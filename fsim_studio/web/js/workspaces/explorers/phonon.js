// Explorer: phonon sidebands, ZPL weight and filter transmission (#/explain/phonon). fsim_core.qd_gf.
//   Plot A  S(omega) (1/meV, log) relative to the ZPL, filter acceptance shaded, Stokes side omega < 0
//   Plot B  Z(T) = exp(-S_total) and S_total on 1..300 K (live, about 2 ms a point)
//   Plot C  transmission vs line offset: ibm_transmission against the Lorentzian-only gap
// Cost rule (brief 3(e)): T, l_xy, l_z, alpha select precomputed FRAMES (one batch job, one worker, so the
// per-process qd_gf caches stay warm); Gamma_zpl, w, kappa, delta, F_cav are live (the first live call at a
// new geometry/T is cold, about 1 s, and says so). No physics here: every number comes from the backend.
import { h, clear } from "../../ui/dom.js";
import { tagChip } from "../../ui/tagChip.js";
import { switchToggle, errorState } from "../../ui/controls.js";
import { fmt, isNum } from "../../ui/format.js";
import { chartCard } from "../../charts/chartCard.js";
import { baseLayout, colors as themeColors, alpha as rgba } from "../../charts/theme.js";
import { ctl, stepper, readoutGroup, caveatBox, qs, jobResult, sourceLine, latest } from "./commonB.js";

const TAG = { v: "A", filter: "A" }; // result / filter-width tags for chart hovers, set from /meta at mount

export async function mount(body, host) {
  const api = host.ctx.api;
  let meta;
  try { meta = await api.getJSON("/api/explore/phonon/meta"); } catch (e) {
    if (host.alive()) clear(body).append(errorState({ title: "The phonon explorer did not answer", message: e.message }));
    return;
  }
  if (!host.alive()) return;
  const D = meta.domain, DEF = meta.defaults, TG = meta.tags, TN = meta.tag_notes;
  const T_frames = meta.T_frames;
  TAG.v = TG.result; TAG.filter = TG.filter;
  const st = { alphaOn: false, kappaOn: true, frames: null, zt: null, cur: null, source: "", framesState: "none" };

  const chartA = chartCard({ title: "Emission spectrum: ZPL plus phonon sidebands", height: 360,
    build: () => buildA(st.cur, st.kappaOn), table: () => tableA(st.cur), filename: () => "fsim-phonon-spectrum" });
  const chartB = chartCard({ title: "ZPL weight and Huang-Rhys factor against T", height: 300,
    build: () => buildB(st.zt, st.cur), table: () => tableB(st.zt), filename: () => "fsim-phonon-zt" });
  const chartC = chartCard({ title: "Filter transmission: IBM against Lorentzian-only", height: 300,
    build: () => buildC(st.cur), table: () => tableC(st.cur), filename: () => "fsim-phonon-transmission" });
  Object.assign(host.charts, { phA: chartA, phB: chartB, phC: chartC });
  for (const c of [chartA, chartB, chartC]) c.setTag(TG.result);

  const caveat = caveatBox([meta.caveat, meta.not_brightness]);
  const status = h("p.xp-status", { "aria-live": "polite" });
  const framesNote = h("p.xb-note", { dataset: { role: "frames" } });

  const P = "xb-ph";
  const S = {
    T: stepper(P, { id: "T", label: "Temperature T (precomputed frame)", values: T_frames, index: T_frames.indexOf(100) >= 0 ? T_frames.indexOf(100) : 5, unit: "K", tag: "A", note: "T selects a precomputed frame", onInput: () => refresh() }),
    l_xy: ctl(P, { id: "l_xy", label: "In-plane extent l_xy", min: D.l_xy[0], max: D.l_xy[1], step: 0.1, default: DEF.l_xy, unit: "nm", tag: TG.geometry, note: TN.geometry }, () => geomChanged()),
    l_z: ctl(P, { id: "l_z", label: "Growth-axis extent l_z", min: D.l_z[0], max: D.l_z[1], step: 0.1, default: DEF.l_z, unit: "nm", tag: TG.geometry, note: TN.geometry }, () => geomChanged()),
    alpha: ctl(P, { id: "alpha", label: "Coupling alpha override", min: D.alpha[0], max: D.alpha[1], step: 0.001, default: meta.alpha_ramsay_ps2, unit: "ps²", tag: TG.alpha, note: TN.alpha }, () => geomChanged()),
    gamma: ctl(P, { id: "gamma", label: "ZPL width Γ_zpl", min: D.gamma_zpl[0], max: D.gamma_zpl[1], default: DEF.gamma_zpl, unit: "meV", scale: "log", tag: "A", note: "phenomenological; the V-a fit Γ(T) is shown beside it [DR]" }, () => filterChanged()),
    w: ctl(P, { id: "w", label: "Filter width w", min: D.w[0], max: D.w[1], step: 0.1, default: DEF.w, unit: "meV", tag: TG.filter, note: TN.filter }, () => filterChanged()),
    kappa: ctl(P, { id: "kappa", label: "Cavity κ", min: D.kappa[0], max: D.kappa[1], step: 0.05, default: 1.0, unit: "meV", tag: TG.cavity, note: TN.filter }, () => filterChanged()),
    F: ctl(P, { id: "F", label: "Purcell factor F_cav", min: D.F_cav[0], max: D.F_cav[1], step: 0.5, default: DEF.F_cav, tag: TG.cavity, note: "total ZPL rate factor supplied by the caller (qd_gf.ibm_purcell_transmission)" }, () => filterChanged()),
    delta: ctl(P, { id: "delta", label: "Line offset δ from filter centre", min: D.delta[0], max: D.delta[1], step: 0.1, default: DEF.delta, unit: "meV", tag: TG.filter }, () => filterChanged()),
  };
  S.alpha.disable(true);
  const alphaSw = switchToggle({ label: "override α (default: computed from the InP-class constants)", checked: false, onChange: (v) => { st.alphaOn = v; S.alpha.disable(!v); geomChanged(); } });
  const kappaSw = switchToggle({ label: "cavity Lorentzian κ in the filter stack", checked: true, onChange: (v) => { st.kappaOn = v; S.kappa.disable(!v); filterChanged(); } });
  const fitBtn = h("button.btn.btn-quiet", { type: "button", dataset: { action: "use-fit" }, on: { click: () => { if (st.cur?.gamma_fit) { S.gamma.set(Math.min(D.gamma_zpl[1], Math.max(D.gamma_zpl[0], Number(st.cur.gamma_fit.toPrecision(3))))); filterChanged(); } } } }, "Use the V-a fit Γ(T)");

  const { el: roEl, ro } = readoutGroup([
    ["Z", "ZPL weight Z(T)"], ["S", "Huang-Rhys S_total"], ["sb", "Sideband fraction 1 - Z"],
    ["t_ibm", "t_IBM at δ"], ["t_lor", "t_Lorentz-only at δ"], ["gap", "Lorentzian optimism gap"],
    ["t_pur", "t with Purcell F_cav"], ["z_eff", "Z_eff = Z F / (Z F + 1 - Z)"], ["rate", "Emission-rate multiplier"],
    ["gfit", "Γ_zpl(T), V-a fit"], ["alpha", "Coupling α"]]);
  const src = h("div");

  body.append(h("div.xp-split",
    h("div.xp-main", caveat.el, chartA.el,
      h("div.xb-grid2", chartB.el, chartC.el),
      h("p.xp-cap", "Plot A: blue is positive; the phonon-emission (Stokes) sideband sits at ω < 0 and dominates at low T. The shaded band is what the filter accepts at the chosen offset. Plot C: a narrow filter loses sideband weight (brightness) while an offset filter admits the other line's sideband (purity); the dashed Lorentzian-only curve is optimistic in both.")),
    h("aside.xp-side.panel",
      h("header.panel-head", h("h2.engrave", "Inputs"), status),
      h("div.xp-ctls", S.T.row, S.l_xy.row, S.l_z.row, h("div.xs-row", alphaSw), S.alpha.row,
        S.gamma.row, h("div.xs-row", fitBtn), S.w.row, h("div.xs-row", kappaSw), S.kappa.row, S.F.row, S.delta.row, framesNote),
      roEl, src)));
  S.kappa.disable(false);

  // ------------------------------------------------------------------ data flow
  const alphaVal = () => (st.alphaOn ? S.alpha.value : null);
  const geom = () => ({ l_xy: S.l_xy.value, l_z: S.l_z.value, alpha: alphaVal() });
  const filt = () => ({ gamma_zpl: S.gamma.value, w: S.w.value, kappa_on: st.kappaOn ? 1 : 0, kappa: st.kappaOn ? S.kappa.value : null, F_cav: S.F.value, delta: S.delta.value });
  const sameInputs = (a, b) => ["l_xy", "l_z", "alpha", "gamma_zpl", "w", "kappa", "F_cav", "delta"].every((k) => (a[k] ?? null) === (b[k] ?? null));
  const currentInputs = () => ({ ...geom(), gamma_zpl: S.gamma.value, w: S.w.value, kappa: st.kappaOn ? S.kappa.value : null, F_cav: S.F.value, delta: S.delta.value });

  const liveJob = latest(async (my, still) => {
    status.textContent = "computing on the job pool";
    try {
      // brief 3(e): the one-frame evaluation runs on the pool (the same phonon_frames job), never in the request thread
      const resp = await api.getJSON(`/api/explore/phonon/live?${qs({ ...geom(), ...filt(), T: S.T.value })}`);
      const res = await jobResult(api, resp);
      if (!still() || !host.alive()) return;
      status.textContent = "";
      show(res.frames[0], "live");
    } catch (e) {
      if (!still() || !host.alive()) return;
      status.textContent = `error: ${e.message}`;
    }
  }, 120);

  function refresh() {
    const f = st.frames;
    const T = S.T.value;
    if (f && sameInputs(f.inputs, currentInputs()) && f.byT[T]) { liveJob.cancel(); status.textContent = ""; show(f.byT[T], "frame"); return; }
    liveJob.run();
  }
  const filterChanged = () => refresh();

  const ztJob = latest(async (my, still) => {
    try {
      const g = geom();
      const r = await api.getJSON(`/api/explore/phonon/zt?${qs(g)}`);
      if (!still() || !host.alive()) return;
      st.zt = r; chartB.render();
    } catch (e) { if (still() && host.alive()) status.textContent = `error: ${e.message}`; }
  }, 200);

  let framesSeq = 0;
  const framesJob = latest(async (my, still) => {
    const mine = ++framesSeq;
    st.framesState = "computing";
    framesNote.textContent = `Computing ${T_frames.length} T frames in one worker (the phonon caches are per process); the live curves answer meanwhile.`;
    try {
      const inputs = currentInputs();
      const resp = await api.postJSON("/api/explore/phonon/frames", { ...geom(), ...filt() });
      const res = await jobResult(api, resp);
      if (mine !== framesSeq || !host.alive()) return;
      st.frames = { inputs, byT: Object.fromEntries(res.frames.map((f) => [f.T, f])) };
      st.framesState = "ready";
      framesNote.textContent = `${res.frames.length} T frames ready (precomputed on the job pool). Moving Γ, w, κ, F or δ switches to live.`;
      refresh();
    } catch (e) {
      if (mine !== framesSeq || !host.alive()) return;
      st.framesState = "error";
      framesNote.textContent = `Frames failed: ${e.message}`;
    }
  }, 500);

  function geomChanged() {
    st.frames = null;
    ztJob.run();
    framesJob.run();
    refresh();
  }

  // ------------------------------------------------------------------ display
  function show(r, source) {
    st.cur = r; st.source = source;
    const t = TG.result;
    const where = `${source === "frame" ? "precomputed frame" : "live"}, T = ${fmt(r.T)} K`;
    ro.Z.update({ value: r.Z, tag: t, caption: `fraction of the line left in the ZPL, ${where}. Not brightness.`, markChange: false });
    ro.S.update({ value: r.S_total, tag: t, caption: "int dE J/E² coth(E/2kT)", markChange: false });
    ro.sb.update({ value: r.sideband_fraction, tag: t, caption: "weight a ZPL-tuned narrow filter can never collect", markChange: false });
    ro.t_ibm.update({ value: r.point.t_ibm, tag: t, caption: `ibm_transmission, δ = ${fmt(r.delta)} meV`, markChange: false });
    ro.t_lor.update({ value: r.point.t_lorentz, tag: t, caption: "spectral.transmission (Lorentzian only)", markChange: false });
    ro.gap.update({ value: r.point.optimism_gap, tag: t, caption: "t_Lorentz-only - t_IBM", markChange: false });
    ro.t_pur.update({ value: r.point.t_purcell, tag: t, caption: `ibm_purcell_transmission, F_cav = ${fmt(r.F_cav)}`, markChange: false });
    ro.z_eff.update({ value: r.point.Z_eff, tag: t, caption: "redistributed ZPL fraction", markChange: false });
    ro.rate.update({ value: r.point.rate_mult, tag: t, caption: "Z F_cav + (1 - Z)", markChange: false });
    ro.gfit.update({ value: r.gamma_fit, unit: "meV", tag: TG.gamma_fit, caption: "effective_gamma_zpl; the a_ac T term partly double-counts acoustic broadening", markChange: false });
    ro.alpha.update({ value: r.alpha_ps2, unit: "ps²", tag: TG.alpha, caption: st.alphaOn ? "override" : "computed from material constants", markChange: false });
    clear(src).append(sourceLine(meta.calls, "Model: independent boson."));
    chartA.setTitle("Emission spectrum: ZPL plus phonon sidebands", `T = ${fmt(r.T)} K · Γ_zpl ${fmt(r.gamma_zpl)} meV · Z = ${fmt(r.Z)} · ${where}`);
    chartC.setTitle("Filter transmission: IBM against Lorentzian-only", `T = ${fmt(r.T)} K · w ${fmt(r.w)} meV${r.kappa != null ? ` · κ ${fmt(r.kappa)} meV` : ""}`);
    chartA.render(); chartC.render(); chartB.render();
  }

  geomChanged();
}

// ---------------------------------------------------------------------- charts
function buildA(r) {
  if (!r) return null;
  const c = themeColors();
  const om = r.spectrum.omega, S = r.spectrum.S;
  const i0 = om.findIndex((x) => x >= -20), i1 = om.findIndex((x) => x > 20);
  const xs = om.slice(i0, i1 < 0 ? undefined : i1), ys = S.slice(i0, i1 < 0 ? undefined : i1).map((v) => Math.max(v, 1e-12));
  const lo = -r.delta - r.w / 2, hi = -r.delta + r.w / 2;
  const layout = baseLayout({
    hovermode: "x unified", margin: { l: 62, r: 18, t: 34, b: 46 },
    legend: { orientation: "h", x: 0, y: -0.2, yanchor: "top", font: { size: 12 } },
    xaxis: { ...baseLayout().xaxis, title: { text: "ω − ω_ZPL (meV), blue positive", standoff: 6 }, range: [-20, 20] },
    yaxis: { ...baseLayout().yaxis, type: "log", title: { text: "S(ω) (1/meV)", standoff: 6 }, exponentformat: "power" },
    shapes: [{ type: "rect", xref: "x", yref: "paper", x0: lo, x1: hi, y0: 0, y1: 1, fillcolor: c.refWash, line: { width: 1, color: c.ref }, layer: "below" }],
    annotations: [
      { x: (lo + hi) / 2, y: 1, yref: "paper", yanchor: "bottom", showarrow: false, text: `filter acceptance, w = ${fmt(r.w)} meV [${TAG.filter}]`, font: { size: 12, color: c.ink2 } },
      { x: -10, y: 0.02, yref: "paper", showarrow: false, text: "Stokes side (ω < 0)", font: { size: 11, color: c.muted } },
    ],
  });
  return { data: [{ x: xs, y: ys, type: "scatter", mode: "lines", name: "S(ω): ZPL + sidebands", line: { color: c.series[0], width: 2 }, hovertemplate: `%{y:.3e} /meV [${TAG.v}]<extra></extra>` }], layout };
}
function tableA(r) {
  if (!r) return null;
  const om = r.spectrum.omega, step = Math.ceil(om.length / 200), rows = [];
  for (let i = 0; i < om.length; i += step) rows.push({ w: om[i], S: r.spectrum.S[i] });
  return { caption: `every ${step}th grid point of ${om.length}`, columns: [{ key: "w", label: "ω − ω_ZPL (meV)" }, { key: "S", label: "S (1/meV)" }], rows };
}

function buildB(zt, r) {
  if (!zt) return null;
  const c = themeColors();
  const data = [
    { x: zt.T, y: zt.Z, type: "scatter", mode: "lines", name: "Z(T) = exp(−S_total)", line: { color: c.series[0], width: 2 }, hovertemplate: `Z %{y:.4f} [${TAG.v}]<extra></extra>` },
    { x: zt.T, y: zt.S_total, type: "scatter", mode: "lines", name: "S_total", yaxis: "y2", line: { color: c.series[1], width: 2, dash: "dot" }, hovertemplate: `S %{y:.4f} [${TAG.v}]<extra></extra>` },
  ];
  if (r) data.push({ x: [r.T], y: [r.Z], type: "scatter", mode: "markers", name: "selected T", marker: { color: c.beam, size: 9, line: { color: c.ink1, width: 1 } }, hoverinfo: "skip" });
  const layout = baseLayout({
    margin: { l: 58, r: 54, t: 34, b: 46 }, legend: { orientation: "h", x: 0, y: -0.22, yanchor: "top", font: { size: 12 } },
    xaxis: { ...baseLayout().xaxis, title: { text: "T (K)", standoff: 6 } },
    yaxis: { ...baseLayout().yaxis, title: { text: "Z (ZPL weight)", standoff: 6 }, range: [0, 1.02] },
    yaxis2: { overlaying: "y", side: "right", title: { text: "S_total", standoff: 6 }, zeroline: false, showgrid: false },
  });
  return { data, layout };
}
function tableB(zt) {
  if (!zt) return null;
  return { columns: [{ key: "T", label: "T (K)" }, { key: "Z", label: "Z" }, { key: "S", label: "S_total" }], rows: zt.T.map((T, i) => ({ T, Z: zt.Z[i], S: zt.S_total[i] })) };
}

function buildC(r) {
  if (!r) return null;
  const c = themeColors();
  const d = r.transmission.delta;
  const data = [
    { x: d, y: r.transmission.ibm, type: "scatter", mode: "lines", name: "IBM (sidebands)", line: { color: c.series[0], width: 2 }, hovertemplate: `IBM %{y:.4f} [${TAG.v}]<extra></extra>` },
    { x: d, y: r.transmission.lorentz, type: "scatter", mode: "lines", name: "Lorentzian only", line: { color: c.series[1], width: 2, dash: "dash" }, hovertemplate: `Lorentz %{y:.4f} [${TAG.v}]<extra></extra>` },
  ];
  const layout = baseLayout({
    margin: { l: 58, r: 18, t: 34, b: 46 }, legend: { orientation: "h", x: 0, y: -0.22, yanchor: "top", font: { size: 12 } },
    xaxis: { ...baseLayout().xaxis, title: { text: "line offset δ from filter centre (meV)", standoff: 6 } },
    yaxis: { ...baseLayout().yaxis, title: { text: "transmitted fraction", standoff: 6 }, range: [0, 1.02] },
    shapes: [{ type: "line", xref: "x", yref: "paper", x0: r.delta, x1: r.delta, y0: 0, y1: 1, line: { color: c.ref, width: 1, dash: "dot" } }],
    annotations: isNum(r.point.optimism_gap) ? [{ x: r.delta, y: r.point.t_lorentz, text: `gap ${fmt(r.point.optimism_gap)} [${TAG.v}]`, showarrow: true, arrowhead: 0, ax: 40, ay: -28, font: { size: 12, color: c.ink1 }, bgcolor: c.surface, borderpad: 2 }] : [],
  });
  return { data, layout };
}
function tableC(r) {
  if (!r) return null;
  return { columns: [{ key: "d", label: "δ (meV)" }, { key: "i", label: "t_IBM" }, { key: "l", label: "t_Lorentz" }],
    rows: r.transmission.delta.map((d, i) => ({ d, i: r.transmission.ibm[i], l: r.transmission.lorentz[i] })) };
}
