// Explorer: injection and Stark tuning (#/explain/transport). fsim_core.transport / nitride_transport / nitride_stark.
//   InP p-i-n   what reaches one dot at (I, T), the background it brings, the field on it
//   Nitride     E_X(V_j), dE_X/dV (stark_derivatives), screening as series, Zhang -10 meV/V as a labelled guide only
// Cost rule (brief 4(e)): single-point readouts and the depletion curve are live; the 50-point I sweep and the
// Stark traces are jobs on the pool. The explorer never labels a slope as agreeing or disagreeing with a
// measurement: that verdict exists only in the committed sweep (CONTRACT rule 4). No physics here.
import { h, clear } from "../../ui/dom.js";
import { segmented, switchToggle, errorState } from "../../ui/controls.js";
import { fmt, isNum } from "../../ui/format.js";
import { chartCard } from "../../charts/chartCard.js";
import { baseLayout, colors as themeColors } from "../../charts/theme.js";
import { ctl, readoutGroup, caveatBox, qs, jobResult, sourceLine, latest } from "./commonB.js";

export async function mount(body, host) {
  const api = host.ctx.api;
  let meta;
  try { meta = await api.getJSON("/api/explore/transport/meta"); } catch (e) {
    if (host.alive()) clear(body).append(errorState({ title: "The transport explorer did not answer", message: e.message }));
    return;
  }
  if (!host.alive()) return;
  const D = meta.domain, DEF = meta.defaults, TG = meta.tags, TN = meta.tag_notes;
  const st = { mode: "inp", preset: DEF.preset, sweep: null, dep: null, pt: null, npt: null, traces: {}, orientation: DEF.orientation, polarity: DEF.polarity, eOn: false, nStatus: "" };

  const caveat = caveatBox([meta.caveat], { id: "transport" });
  const nCaveat = caveatBox([meta.nitride_caveat], { id: "nitride" });
  const status = h("p.xp-status", { "aria-live": "polite" });
  const modeSeg = segmented({ label: "Device", value: "inp", options: [{ value: "inp", label: "InP p-i-n" }, { value: "nit", label: "Nitride Stark" }], onChange: (v) => { st.mode = v; syncMode(); } });

  // ------------------------------------------------------------------ InP panel
  const cV = chartCard({ title: "Junction voltage and injection efficiency against I", height: 300, build: () => buildV(st.sweep, st.pt, TG), table: () => tbl(st.sweep, ["I_uA", "V_j", "V_applied", "eta_inj"]), filename: () => "fsim-transport-voltage" });
  const cM = chartCard({ title: "Dot loading μ and background b_e against I", height: 300, build: () => buildMuB(st.sweep, st.pt, TG), table: () => tbl(st.sweep, ["I_uA", "mu", "b_e"]), filename: () => "fsim-transport-mu-be" });
  const cR = chartCard({ title: "Per-dot capture rate and junction power against I", height: 300, build: () => buildRP(st.sweep, st.pt, TG), table: () => tbl(st.sweep, ["I_uA", "r_dot", "r_captured", "r_matrix", "P_junction_W"]), filename: () => "fsim-transport-rate-power" });
  const cD = chartCard({ title: "Depletion field and capacitance against V_j", height: 300, build: () => buildDep(st.dep, st.pt, TG), table: () => (st.dep ? { columns: [{ key: "V", label: "V_j (V)" }, { key: "F", label: "F (kV/cm)" }, { key: "C", label: "C_dep (pF)" }], rows: st.dep.V_j.map((V, i) => ({ V, F: st.dep.F_kVcm[i], C: st.dep.C_dep_pF[i] })) } : null), filename: () => "fsim-transport-depletion" });
  // ------------------------------------------------------------------ nitride panel
  const cE = chartCard({ title: "Exciton energy E_X against junction voltage V_j", height: 340, build: () => buildE(st, meta), table: () => stTable(st), filename: () => "fsim-nitride-stark-energy" });
  const cS = chartCard({ title: "Stark slope dE_X/dV against V_j", height: 340, build: () => buildS(st, meta), table: () => stTable(st), filename: () => "fsim-nitride-stark-slope" });
  Object.assign(host.charts, { tV: cV, tM: cM, tR: cR, tD: cD, tE: cE, tS: cS });
  for (const c of [cV, cM, cR, cD]) c.setTag(TG.inputs);
  for (const c of [cE, cS]) c.setTag(TG.nitride);

  const P = "xb-tr";
  const A = {
    I: ctl(P, { id: "I", label: "Current I", min: D.I_uA[0], max: D.I_uA[1], default: DEF.I_uA, unit: "µA", scale: "log", tag: "A", note: TN.inputs }, () => pointJob.run()),
    T: ctl(P, { id: "T", label: "Temperature T", min: D.T[0], max: D.T[1], step: 1, default: DEF.T, unit: "K", tag: "A" }, () => { pointJob.run(); sweepJob.run(); depJob.run(); }),
    n_dot: ctl(P, { id: "n_dot", label: "Dot density n_dot", min: D.n_dot[0], max: D.n_dot[1], default: DEF.n_dot, unit: "cm⁻²", scale: "log", tag: "A", note: TN.inputs }, () => { pointJob.run(); sweepJob.run(); }),
    aperture: ctl(P, { id: "aperture", label: "Aperture", min: D.aperture[0], max: D.aperture[1], step: 0.05, default: DEF.aperture, unit: "µm²", tag: "A", note: TN.inputs }, () => { pointJob.run(); sweepJob.run(); }),
    tau_pulse: ctl(P, { id: "tau_pulse", label: "Pulse length τ_pulse", min: D.tau_pulse[0], max: D.tau_pulse[1], default: DEF.tau_pulse, unit: "ns", scale: "log", tag: "A" }, () => { pointJob.run(); sweepJob.run(); }),
    w: ctl(P, { id: "w", label: "Window w", min: D.w[0], max: D.w[1], step: 0.1, default: DEF.w, unit: "meV", tag: "E", note: TN.inputs }, () => { pointJob.run(); sweepJob.run(); }),
    tau_rad: ctl(P, { id: "tau_rad", label: "Radiative lifetime τ_rad", min: D.tau_rad[0], max: D.tau_rad[1], step: 0.05, default: DEF.tau_rad, unit: "ns", tag: "A" }, () => { pointJob.run(); sweepJob.run(); }),
    E_X: ctl(P, { id: "E_X", label: "Dot energy E_X (sub-turn-on suppression)", min: 1.5, max: 2.6, step: 0.01, default: 1.9, unit: "eV", tag: "A", note: "optional: switches on f_qfl = min(1, exp(-(E_X - V_j)/kT))" }, () => pointJob.run()),
  };
  A.E_X.disable(true);
  const eSw = switchToggle({ label: "supply E_X (sub-turn-on loading suppression)", checked: false, onChange: (v) => { st.eOn = v; A.E_X.disable(!v); pointJob.run(); } });
  const presetSeg = segmented({ label: "Diode preset", value: st.preset, options: meta.presets.map((p) => ({ value: p, label: { red: "red AlGaInP", hkust: "HKUST AlGaAs", gaas: "GaAs p-i-n" }[p] || p })), onChange: (v) => { st.preset = v; pointJob.run(); sweepJob.run(); depJob.run(); } });
  const inp = readoutGroup([
    ["V_j", "Junction voltage V_j"], ["V_app", "Terminal voltage"], ["V_bi", "Built-in V_bi"], ["F", "Field F"], ["C", "C_dep"],
    ["eta", "Injection efficiency η_inj"], ["mu", "Loading μ"], ["b_e", "Background b_e"], ["P", "Junction power"],
    ["r_dot", "Per-dot rate r_dot"], ["rcap", "Captured rate"], ["rmat", "Matrix rate"], ["xi", "Window factor ξ"], ["bgfit", "Legacy b_e fit error"]]);

  // nitride controls
  const N = {
    T: ctl(P, { id: "nT", label: "Junction temperature T", min: 230, max: 350, step: 1, default: DEF.T_nitride, unit: "K", tag: "A" }, () => { nPointJob.run(); traceJob.run(); }),
    I: ctl(P, { id: "nI", label: "Current for the point readout", min: 1e-3, max: 10, default: 0.02, unit: "µA", scale: "log", tag: "A", note: "qcap-staged.yaml:54-57 [A]" }, () => nPointJob.run()),
    ext: ctl(P, { id: "next", label: "External field", min: D.ext[0] / 5, max: D.ext[1] / 5, step: 1, default: 0, unit: "kV/cm", tag: "A" }, () => { nPointJob.run(); traceJob.run(); }),
  };
  const oriSeg = segmented({ label: "Orientation", value: st.orientation, options: [{ value: "c_plane", label: "c-plane" }, { value: "a_plane", label: "a-plane" }], onChange: (v) => { st.orientation = v; traceJob.run(); } });
  const polSeg = segmented({ label: "Field polarity", value: String(st.polarity), options: [{ value: "1", label: "+1" }, { value: "-1", label: "−1" }], onChange: (v) => { st.polarity = Number(v); nPointJob.run(); traceJob.run(); } });
  const nit = readoutGroup([
    ["V_j", "Junction voltage V_j"], ["Vt", "Terminal voltage"], ["Fd", "Diode field"], ["Fa", "Applied field"],
    ["mu", "Loading μ"], ["b_e", "Background b_e"], ["slope", "Slope at the trace start"], ["guide", "Zhang 2016 guide"]]);
  const nNote = h("p.xb-note", { dataset: { role: "stark-status" } });

  const inpSide = h("div.xp-ctls", h("div.xs-row", h("span.xs-label", "Diode preset"), presetSeg.el), A.I.row, A.T.row, A.n_dot.row, A.aperture.row, A.tau_pulse.row, A.w.row, A.tau_rad.row, h("div.xs-row", eSw), A.E_X.row);
  const nitSide = h("div.xp-ctls", h("div.xs-row", h("span.xs-label", "Orientation"), oriSeg.el), h("div.xs-row", h("span.xs-label", "Field polarity"), polSeg.el), N.T.row, N.I.row, N.ext.row, nNote);
  const inpMain = h("div.xb-panel", { dataset: { panel: "inp" } }, caveat.el, h("div.xb-grid2", cV.el, cM.el, cR.el, cD.el),
    h("p.xp-cap", "Dot loading μ rises roughly linearly with I until the dot saturates at 1/τ_rad; b_e is the matrix background per X photon and is flat in I below saturation, linear above. Dashed line: μ = 0.5 (cap-2)."));
  const nitMain = h("div.xb-panel", { dataset: { panel: "nit" } }, nCaveat.el, h("div.xb-grid2", cE.el, cS.el),
    h("p.xp-cap", "Screening is fixed along each trace and the three series are the explorer's assumption [A]. The Zhang slope is drawn only as a labelled guide: it comes from another device (10 K, x = 0.15), and this explorer does not compare it with the model."));
  const inpRO = h("div.xb-panel", { dataset: { panel: "inp" } }, inp.el, h("div", { dataset: { role: "inp-src" } }));
  const nitRO = h("div.xb-panel", { dataset: { panel: "nit" } }, nit.el, h("div", { dataset: { role: "nit-src" } }));

  body.append(h("div.xp-split",
    h("div.xp-main", h("div.xb-modebar", modeSeg.el), inpMain, nitMain),
    h("aside.xp-side.panel", h("header.panel-head", h("h2.engrave", "Inputs"), status),
      h("div.xb-panel", { dataset: { panel: "inp" } }, inpSide), h("div.xb-panel", { dataset: { panel: "nit" } }, nitSide), inpRO, nitRO)));

  function syncMode() {
    for (const el of body.querySelectorAll(".xb-panel")) el.hidden = el.dataset.panel !== st.mode;
    for (const c of st.mode === "inp" ? [cV, cM, cR, cD] : [cE, cS]) c.render();
    if (st.mode === "nit") { nPointJob.run(); traceJob.run(); }
  }

  // ------------------------------------------------------------------ InP data flow
  const inpArgs = () => ({ preset: st.preset, I_uA: A.I.value, T: A.T.value, n_dot: A.n_dot.value, aperture: A.aperture.value, tau_pulse: A.tau_pulse.value, w: A.w.value, dE_WL: DEF.dE_WL, tau_rad: A.tau_rad.value, E_X_eV: st.eOn ? A.E_X.value : null });
  const pointJob = latest(async (my, still) => {
    status.textContent = "computing";
    try {
      const r = await api.getJSON(`/api/explore/transport/point?${qs(inpArgs())}`);
      if (!still() || !host.alive()) return;
      status.textContent = ""; st.pt = r; showPoint(r);
    } catch (e) { if (still() && host.alive()) status.textContent = `error: ${e.message}`; }
  }, 100);
  const sweepJob = latest(async (my, still) => {
    status.textContent = "sweep running on the job pool";
    try {
      const a = inpArgs(); delete a.I_uA; delete a.E_X_eV;
      const resp = await api.postJSON("/api/explore/transport/sweep", a);
      const res = await jobResult(api, resp);
      if (!still() || !host.alive()) return;
      status.textContent = ""; st.sweep = res;
      inp.ro.bgfit.update({ value: res.bg_fit.max_rel_err, tag: TG.b_e, caption: `max |fit/computed - 1| of A (I/I_ref)^m exp(-E_act/kT) on ${fmt(res.bg_fit.I_range_uA[0])}-${fmt(res.bg_fit.I_range_uA[1])} µA`, markChange: false });
      for (const c of [cV, cM, cR]) c.render();
    } catch (e) { if (still() && host.alive()) status.textContent = `error: ${e.message}`; }
  }, 400);
  const depJob = latest(async (my, still) => {
    try {
      const r = await api.getJSON(`/api/explore/transport/depletion?${qs({ preset: st.preset, T: A.T.value })}`);
      if (!still() || !host.alive()) return;
      st.dep = r; cD.render();
    } catch (e) { if (still() && host.alive()) status.textContent = `error: ${e.message}`; }
  }, 200);

  function showPoint(r) {
    const R = inp.ro;
    R.V_j.update({ value: r.V_j, unit: "V", tag: TG.V_j, caption: "junction voltage at this I", markChange: false });
    R.V_app.update({ value: r.V_applied, unit: "V", tag: TG.V_j, caption: "V_j + I R_s", markChange: false });
    R.V_bi.update({ value: r.V_bi, unit: "V", tag: TG.V_bi, caption: r.flat_band ? "flat band: V_j has reached V_bi" : "abrupt depletion", markChange: false });
    R.F.update({ value: r.F_kVcm, unit: "kV/cm", tag: TG.F, caption: "uniform in the i-region", markChange: false });
    R.C.update({ value: r.C_dep_pF, unit: "pF", tag: TG.C_dep, markChange: false });
    R.eta.update({ value: r.eta_inj, tag: TG.eta_inj, caption: `thermionic leakage, valleys ${r.leak_valleys} (a LOWER bound on leakage)`, markChange: false });
    R.mu.update({ value: r.mu, tag: TG.mu, caption: r.saturated ? "dot saturated (r_dot > 1/τ_rad)" : "per-pulse loading", markChange: false });
    R.b_e.update({ value: r.b_e, tag: TG.b_e, caption: "background per X photon; excludes neighbour dots and the X filter", markChange: false });
    R.P.update({ value: r.P_junction_W, unit: "W", tag: TG.P_junction, markChange: false });
    R.r_dot.update({ value: r.r_dot, unit: "s⁻¹", tag: TG.r_dot, caption: `saturates at 1/τ_rad = ${fmt(r.r_dot_max)} s⁻¹`, markChange: false });
    R.rcap.update({ value: r.r_captured, unit: "s⁻¹", tag: TG.eta_inj, caption: "r_captured + r_matrix = η_inj I / q", markChange: false });
    R.rmat.update({ value: r.r_matrix, unit: "s⁻¹", tag: TG.eta_inj, caption: `f_qfl = ${fmt(r.f_qfl)}`, markChange: false });
    R.xi.update({ value: r.xi, tag: TG.xi, caption: `Urbach tail E_U = ${fmt(r.E_U_meV)} meV`, markChange: false });
    for (const c of [cV, cM, cR, cD]) c.render();
    const slot = body.querySelector('[data-role="inp-src"]');
    if (slot) clear(slot).append(sourceLine(meta.calls.slice(0, 3), "Compact p-i-n [E/A]."));
  }

  // ------------------------------------------------------------------ nitride data flow
  const nPointJob = latest(async (my, still) => {
    if (st.mode !== "nit") return;
    try {
      const r = await api.getJSON(`/api/explore/transport/nitride_point?${qs({ I_uA: N.I.value, T: N.T.value, n_dot: DEF.n_dot, aperture: DEF.aperture, tau_pulse: DEF.tau_pulse, w: DEF.w, dE_WL: DEF.dE_WL, tau_rad: DEF.tau_rad, polarity: st.polarity, ext: N.ext.value })}`);
      if (!still() || !host.alive()) return;
      st.npt = r;
      const R = nit.ro, b = r.bias;
      R.V_j.update({ value: b.V_j, unit: "V", tag: TG.nitride, caption: "resolve_bias from the current", markChange: false });
      R.Vt.update({ value: b.V_terminal, unit: "V", tag: TG.nitride, markChange: false });
      R.Fd.update({ value: b.diode_field_kVcm, unit: "kV/cm", tag: TG.nitride, caption: `abrupt depletion [DR], regime: ${b.depletion_regime}${b.flat_band ? " (flat band)" : ""}`, markChange: false });
      R.Fa.update({ value: b.applied_field_kVcm, unit: "kV/cm", tag: TG.nitride, caption: `polarity ${st.polarity > 0 ? "+1" : "-1"} · external ${fmt(N.ext.value)} kV/cm`, markChange: false });
      R.mu.update({ value: r.mu, tag: TG.mu, caption: "nitride_transport.evaluate_injection", markChange: false });
      R.b_e.update({ value: r.b_e, tag: TG.b_e, caption: "barrier-limited, doping dependence omitted [A]", markChange: false });
      const slot = body.querySelector('[data-role="nit-src"]');
      if (slot) clear(slot).append(sourceLine(meta.calls.slice(3), "Nitride diode: Zhang 2016 defaults [V source, E transfer]."));
    } catch (e) { if (still() && host.alive()) status.textContent = `error: ${e.message}`; }
  }, 100);

  let traceSeq = 0;
  const traceJob = latest(async (my, still) => {
    if (st.mode !== "nit") return;
    const mine = ++traceSeq;
    st.traces = {};
    nNote.textContent = "Computing the Stark traces on the job pool (E_X from nitride_levels.levels at each bias, 0.2 V steps).";
    try {
      const resp = await api.postJSON("/api/explore/transport/stark", { orientation: st.orientation, polarity: st.polarity, T: N.T.value, ext: N.ext.value, screenings: meta.screenings });
      let done = 0;
      await Promise.all(resp.traces.map(async (t) => {
        const res = await jobResult(api, t);
        if (mine !== traceSeq || !host.alive()) return;
        st.traces[t.screening_fraction] = res; done++;
        nNote.textContent = `${done}/${resp.traces.length} traces ready.`;
        cE.render(); cS.render(); showSlope();
      }));
      if (mine === traceSeq && host.alive()) nNote.textContent = `${resp.traces.length} traces ready (one point per 0.2 V; derivatives never cross an invalid row or a regime kink).`;
    } catch (e) { if (mine === traceSeq && host.alive()) nNote.textContent = `Stark traces failed: ${e.message}`; }
  }, 400);

  function showSlope() {
    const tr = st.traces[0] || Object.values(st.traces)[0];
    if (!tr) return;
    const first = tr.rows.find((r) => r.derivative_valid);
    nit.ro.slope.update({ value: first ? first.dE_X_dV_meV_per_V : null, unit: "meV/V", tag: TG.slope, reason: "no valid derivative on this trace", caption: first ? `V_j = ${fmt(first.V_j)} V, screening ${fmt(tr.screening_fraction)} [A]` : "", markChange: false });
    nit.ro.guide.update({ value: meta.zhang.slope_meV_per_V, unit: "meV/V", tag: meta.zhang.tag, caption: meta.zhang.label, markChange: false });
  }

  syncMode();
  pointJob.run(); sweepJob.run(); depJob.run();
}

// ---------------------------------------------------------------------- charts
function vline(x, c) { return { type: "line", xref: "x", yref: "paper", x0: x, x1: x, y0: 0, y1: 1, line: { color: c.ref, width: 1, dash: "dot" } }; }
const LOGX = (extra = {}) => ({ ...baseLayout().xaxis, type: "log", title: { text: "current I (µA)", standoff: 6 }, ...extra });
const LEG = { orientation: "h", x: 0, y: -0.24, yanchor: "top", font: { size: 12 } };

function buildV(sw, pt, TG) {
  if (!sw) return null;
  const c = themeColors(), k = sw.columns;
  const data = [
    { x: k.I_uA, y: k.V_j, type: "scatter", mode: "lines", name: "V_j", line: { color: c.series[0], width: 2 }, hovertemplate: `V_j %{y:.3f} V [${TG.V_j}]<extra></extra>` },
    { x: k.I_uA, y: k.V_applied, type: "scatter", mode: "lines", name: "V_applied", line: { color: c.series[1], width: 2 }, hovertemplate: `V_app %{y:.3f} V [${TG.V_j}]<extra></extra>` },
    { x: k.I_uA, y: k.eta_inj, type: "scatter", mode: "lines", name: "η_inj", yaxis: "y2", line: { color: c.series[2], width: 2, dash: "dot" }, hovertemplate: `η_inj %{y:.4f} [${TG.eta_inj}]<extra></extra>` },
  ];
  return { data, layout: baseLayout({ margin: { l: 58, r: 54, t: 30, b: 50 }, legend: LEG, xaxis: LOGX(),
    yaxis: { ...baseLayout().yaxis, title: { text: "voltage (V)", standoff: 6 } }, yaxis2: { overlaying: "y", side: "right", title: { text: "η_inj", standoff: 6 }, range: [0, 1.05], showgrid: false, zeroline: false },
    shapes: pt ? [vline(pt.I_uA, c)] : [] }) };
}
function buildMuB(sw, pt, TG) {
  if (!sw) return null;
  const c = themeColors(), k = sw.columns;
  const data = [
    { x: k.I_uA, y: k.mu, type: "scatter", mode: "lines", name: "μ (per-pulse loading)", line: { color: c.series[0], width: 2 }, hovertemplate: `μ %{y:.3e} [${TG.mu}]<extra></extra>` },
    { x: k.I_uA, y: k.b_e, type: "scatter", mode: "lines", name: "b_e (background per X photon)", line: { color: c.series[1], width: 2 }, hovertemplate: `b_e %{y:.3e} [${TG.b_e}]<extra></extra>` },
  ];
  return { data, layout: baseLayout({ margin: { l: 62, r: 18, t: 30, b: 50 }, legend: LEG, xaxis: LOGX(), yaxis: { ...baseLayout().yaxis, type: "log", title: { text: "μ, b_e", standoff: 6 }, exponentformat: "power" },
    shapes: [...(pt ? [vline(pt.I_uA, c)] : []), { type: "line", xref: "paper", yref: "y", x0: 0, x1: 1, y0: 0.5, y1: 0.5, line: { color: c.ink2, width: 1, dash: "dash" } }],
    annotations: [{ xref: "paper", x: 0.01, y: Math.log10(0.5), yref: "y", text: "μ = 0.5 [A]", showarrow: false, yshift: 10, font: { size: 11, color: c.ink2 }, xanchor: "left" }] }) };
}
function buildRP(sw, pt, TG) {
  if (!sw) return null;
  const c = themeColors(), k = sw.columns;
  const data = [
    { x: k.I_uA, y: k.r_dot, type: "scatter", mode: "lines", name: "r_dot", line: { color: c.series[0], width: 2 }, hovertemplate: `r_dot %{y:.3e} /s [${TG.r_dot}]<extra></extra>` },
    { x: k.I_uA, y: k.P_junction_W, type: "scatter", mode: "lines", name: "P_junction", yaxis: "y2", line: { color: c.series[2], width: 2, dash: "dot" }, hovertemplate: `P %{y:.3e} W [${TG.P_junction}]<extra></extra>` },
  ];
  return { data, layout: baseLayout({ margin: { l: 62, r: 62, t: 30, b: 50 }, legend: LEG, xaxis: LOGX(),
    yaxis: { ...baseLayout().yaxis, type: "log", title: { text: "r_dot (1/s)", standoff: 6 }, exponentformat: "power" },
    yaxis2: { overlaying: "y", side: "right", type: "log", title: { text: "P_junction (W)", standoff: 6 }, showgrid: false, zeroline: false, exponentformat: "power" },
    shapes: [...(pt ? [vline(pt.I_uA, c)] : []), { type: "line", xref: "paper", yref: "y", x0: 0, x1: 1, y0: sw.r_dot_max, y1: sw.r_dot_max, line: { color: c.ink2, width: 1, dash: "dash" } }],
    annotations: [{ xref: "paper", x: 0.01, y: Math.log10(sw.r_dot_max), yref: "y", text: "1/τ_rad (dot saturation)", showarrow: false, yshift: 10, font: { size: 11, color: c.ink2 }, xanchor: "left" }] }) };
}
function buildDep(dep, pt, TG) {
  if (!dep) return null;
  const c = themeColors();
  const data = [
    { x: dep.V_j, y: dep.F_kVcm, type: "scatter", mode: "lines", name: "F (kV/cm)", line: { color: c.series[0], width: 2 }, hovertemplate: `F %{y:.2f} kV/cm [${TG.F}]<extra></extra>` },
    { x: dep.V_j, y: dep.C_dep_pF, type: "scatter", mode: "lines", name: "C_dep (pF)", yaxis: "y2", line: { color: c.series[1], width: 2, dash: "dot" }, hovertemplate: `C %{y:.4e} pF [${TG.C_dep}]<extra></extra>` },
  ];
  return { data, layout: baseLayout({ margin: { l: 58, r: 62, t: 30, b: 50 }, legend: LEG,
    xaxis: { ...baseLayout().xaxis, title: { text: "junction voltage V_j (V); flat band beyond V_bi", standoff: 6 } },
    yaxis: { ...baseLayout().yaxis, title: { text: "F (kV/cm)", standoff: 6 }, rangemode: "tozero" },
    yaxis2: { overlaying: "y", side: "right", title: { text: "C_dep (pF)", standoff: 6 }, showgrid: false, zeroline: false, rangemode: "tozero" },
    shapes: [{ type: "line", xref: "x", yref: "paper", x0: dep.V_bi, x1: dep.V_bi, y0: 0, y1: 1, line: { color: c.ink2, width: 1, dash: "dash" } }, ...(pt ? [vline(pt.V_j, c)] : [])],
    annotations: [{ x: dep.V_bi, y: 1, yref: "paper", yanchor: "bottom", showarrow: false, text: `V_bi = ${fmt(dep.V_bi)} V [${TG.V_bi}]`, font: { size: 11, color: c.ink2 } }] }) };
}
function tbl(sw, keys) {
  if (!sw) return null;
  return { columns: keys.map((k) => ({ key: k, label: k })), rows: sw.columns.I_uA.map((_, i) => Object.fromEntries(keys.map((k) => [k, sw.columns[k][i]]))) };
}

const SCR = [0, 0.5, 1];
function series(st) { return SCR.filter((s) => st.traces[s]).map((s) => ({ s, tr: st.traces[s] })); }
function buildE(st, meta) {
  const list = series(st);
  if (!list.length) return null;
  const c = themeColors();
  const data = list.map(({ s, tr }, i) => ({ x: tr.rows.map((r) => r.V_j), y: tr.rows.map((r) => (r.spectroscopy_valid ? r.E_X_eV : null)), type: "scatter", mode: "lines+markers", name: `screening ${fmt(s)} [A]`, line: { color: c.series[i], width: 2 }, marker: { size: 5 }, connectgaps: false, hovertemplate: `E_X %{y:.4f} eV, screening ${fmt(s)} [${meta.tags.E_X}]<extra></extra>` }));
  const base = list[0].tr.rows.find((r) => r.spectroscopy_valid);
  if (base) {
    const m = meta.zhang.slope_meV_per_V / 1000, x1 = list[0].tr.rows[list[0].tr.rows.length - 1].V_j;
    data.push({ x: [base.V_j, x1], y: [base.E_X_eV, base.E_X_eV + m * (x1 - base.V_j)], type: "scatter", mode: "lines", name: "Zhang −10 meV/V guide only", line: { color: c.ink2, width: 1, dash: "dot" }, hoverinfo: "skip" });
  }
  const flat = list[0].tr.rows.filter((r) => r.flat_band).map((r) => r.V_j);
  return { data, layout: baseLayout({ margin: { l: 62, r: 18, t: 30, b: 50 }, legend: LEG,
    xaxis: { ...baseLayout().xaxis, title: { text: "junction voltage V_j (V)", standoff: 6 } },
    yaxis: { ...baseLayout().yaxis, title: { text: "E_X (eV)", standoff: 6 } },
    annotations: flat.length ? [{ x: flat[0], y: 1, yref: "paper", text: "flat band", showarrow: false, font: { size: 11, color: c.ink2 } }] : [] }) };
}
function buildS(st, meta) {
  const list = series(st);
  if (!list.length) return null;
  const c = themeColors();
  const data = list.map(({ s, tr }, i) => ({ x: tr.rows.map((r) => r.V_j), y: tr.rows.map((r) => (r.derivative_valid ? r.dE_X_dV_meV_per_V : null)), type: "scatter", mode: "lines+markers", name: `screening ${fmt(s)} [A]`, line: { color: c.series[i], width: 2 }, marker: { size: 5 }, connectgaps: false, hovertemplate: `dE_X/dV %{y:.2f} meV/V, screening ${fmt(s)} [${meta.tags.slope}]<extra></extra>` }));
  const z = meta.zhang.slope_meV_per_V;
  const layout = baseLayout({ margin: { l: 62, r: 18, t: 30, b: 50 }, legend: LEG,
    xaxis: { ...baseLayout().xaxis, title: { text: "junction voltage V_j (V)", standoff: 6 } },
    yaxis: { ...baseLayout().yaxis, title: { text: "dE_X/dV (meV/V)", standoff: 6 } },
    shapes: [{ type: "line", xref: "paper", yref: "y", x0: 0, x1: 1, y0: z, y1: z, line: { color: c.ink2, width: 1, dash: "dot" } }],
    annotations: [{ xref: "paper", x: 0.99, xanchor: "right", y: z, yref: "y", yshift: 10, showarrow: false, text: "Zhang −10 meV/V: guide only, 10 K, x = 0.15, another device [V]", font: { size: 11, color: c.ink2 } }] });
  return { data, layout };
}
function stTable(st) {
  const list = series(st);
  if (!list.length) return null;
  const rows = [];
  for (const { s, tr } of list) for (const r of tr.rows) rows.push({ s, V: r.V_j, E: r.E_X_eV, d: r.dE_X_dV_meV_per_V, ok: r.derivative_valid ? "valid" : "no slope", reg: r.depletion_regime });
  return { columns: [{ key: "s", label: "screening" }, { key: "V", label: "V_j (V)" }, { key: "E", label: "E_X (eV)" }, { key: "d", label: "dE_X/dV (meV/V)" }, { key: "ok", label: "slope" }, { key: "reg", label: "regime" }], rows };
}
