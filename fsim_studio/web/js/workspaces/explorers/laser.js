// Explorer: Langevin QD laser (Zhao) (#/explain/laser). fsim_core.sde.
//   Plot A  L-I: steady-state S against I/I_th, I_th marked
//   Plot B  g2(0) against I/I_th for F_pump 1 and 0.08 with +-SE, the Zhao card points (V) and the committed fit run
//   Plot C  one simulate() S(t) trace, burn-in shaded
// Cost rule (brief 5(e)): every point is a seeded frame computed on the job pool; sliders SELECT frames and never
// integrate synchronously. The pinned caveat is Zhao's: this is a laser, g2 ~ 1 is correct, the quiet pump does
// NOT reproduce the card (5.8 sigma, FAIL), and no curve is drawn through the quiet point. No physics here.
import { h, clear } from "../../ui/dom.js";
import { segmented, errorState } from "../../ui/controls.js";
import { tagChip } from "../../ui/tagChip.js";
import { fmt, isNum } from "../../ui/format.js";
import { chartCard } from "../../charts/chartCard.js";
import { baseLayout, colors as themeColors, alpha as rgba } from "../../charts/theme.js";
import { ctl, stepper, readoutGroup, caveatBox, jobResult, sourceLine, latest } from "./commonB.js";

export async function mount(body, host) {
  const api = host.ctx.api;
  let meta, card;
  try {
    [meta, card] = await Promise.all([api.getJSON("/api/explore/sde/meta"), api.getJSON("/api/explore/sde/card")]);
  } catch (e) {
    if (host.alive()) clear(body).append(errorState({ title: "The laser explorer did not answer", message: e.message }));
    return;
  }
  if (!host.alive()) return;
  const D = meta.domain, DEF = meta.defaults, TG = meta.tags;
  const st = { F: 1.0, zoom: "near", frames: {}, li: null, trace: null, seq: 0, nDone: 0, nTotal: 0 };

  const chartA = chartCard({ title: "L-I: steady-state photon number against pump", height: 300, build: () => buildLI(st, meta, S.I.value), table: () => (st.li ? { columns: [{ key: "I", label: "I/I_th (sde axis)" }, { key: "S", label: "S (photons)" }, { key: "rg", label: "rho_GS" }, { key: "re", label: "rho_ES" }], rows: st.li.I.map((I, i) => ({ I, S: st.li.S[i], rg: st.li.rho_GS[i], re: st.li.rho_ES[i] })) } : null), filename: () => "fsim-laser-li" });
  const chartB = chartCard({ title: "g²(0) against pump: normal and quiet drive", height: 380, build: () => buildG2(st, card, meta), table: () => tableG2(st), filename: () => "fsim-laser-g2" });
  const chartC = chartCard({ title: "One seeded S(t) trace", height: 280, build: () => buildTrace(st), table: () => (st.trace ? { columns: [{ key: "t", label: "t (ns)" }, { key: "S", label: "S" }], rows: st.trace.t_ns.map((t, i) => ({ t, S: st.trace.S[i] })) } : null), filename: () => "fsim-laser-trace" });
  Object.assign(host.charts, { lA: chartA, lB: chartB, lC: chartC });
  for (const c of [chartA, chartB, chartC]) c.setTag(TG.result);

  const caveat = caveatBox([meta.caveat, meta.axis_note], { id: "zhao" });
  const status = h("p.xp-status", { "aria-live": "polite" });
  const note = h("p.xb-note", { dataset: { role: "frames" } });

  const P = "xb-ls";
  const S = {
    I: stepper(P, { id: "I", label: "Pump I/I_th (selects a frame)", values: meta.I_grid, index: Math.max(0, meta.I_grid.indexOf(DEF.I)), tag: TG.model, note: "Axis is I/I_th, not amperes (sde.py:196-203)", onInput: () => select() }),
    n_runs: ctl(P, { id: "n_runs", label: "Runs per frame", min: D.n_runs[0], max: D.n_runs[1], step: 1, default: DEF.n_runs, tag: "A", note: "more runs: smaller SE, longer job" }, () => dirty()),
    t_end: ctl(P, { id: "t_end", label: "Run length t_end", min: D.t_end[0] * 1e9, max: D.t_end[1] * 1e9, step: 0.5, default: DEF.t_end * 1e9, unit: "ns", tag: "A", note: "first 10% is discarded as SDE transient (BURNIN_FRAC)" }, () => dirty()),
  };
  const seed = h("input.fr-num", { type: "number", min: 0, step: 1, value: DEF.seed, "aria-label": "Random seed", dataset: { ctl: "seed" }, on: { change: () => dirty() } });
  const fSeg = segmented({ label: "Pump statistics shown in the readouts", value: "1", options: [{ value: "1", label: "normal F = 1" }, { value: "0.08", label: "quiet F = 0.08" }], onChange: (v) => { st.F = Number(v); select(); traceJob.run(); } });
  const zoomSeg = segmented({ label: "g2 axis", value: "near", options: [{ value: "near", label: "near 1" }, { value: "full", label: "full range" }], onChange: (v) => { st.zoom = v; chartB.render(); } });
  const runBtn = h("button.btn.btn-outline", { type: "button", dataset: { action: "compute" }, on: { click: () => compute() } }, "Compute frames");
  const { el: roEl, ro } = readoutGroup([
    ["g2", "g²(0), selected frame"], ["S", "Mean S, selected frame"], ["Ith", "Threshold I_th"], ["Sss", "Steady-state S (noise-free)"],
    ["card", "Zhao card point at I/I_th = 4"], ["sigma", "Fit against the card (committed run)"], ["tg2", "Trace g²(0)"], ["fano", "Photon Fano factor"]]);
  const src = h("div");

  body.append(h("div.xp-split",
    h("div.xp-main", caveat.el, chartB.el, h("div.xb-grid2", chartA.el, chartC.el),
      h("p.xp-cap", "Each marker is one seeded run batch with its standard error; the hollow diamonds are the committed fit run (out/zhao). The Zhao card points are drawn at I/I_th = 4 with their error bars and are not connected to the simulation. The quiet pump is expected to sit below the normal pump in sign only; the model does not reach the published magnitude.")),
    h("aside.xp-side.panel", h("header.panel-head", h("h2.engrave", "Inputs"), status),
      h("div.xp-ctls", S.I.row, h("div.xs-row", h("span.xs-label", "Series in the readouts ", tagChip(TG.model, { size: "sm" })), fSeg.el), S.n_runs.row, S.t_end.row,
        h("div.xs-row", h("label.xs-label", { for: "xb-ls-seed" }, "Seed ", tagChip("A", { size: "sm" })), h("div.xs-ctl", Object.assign(seed, { id: "xb-ls-seed" }))),
        h("div.xs-row", h("span.xs-label", "g² axis"), zoomSeg.el), h("div.xs-row", runBtn), note),
      roEl, src)));

  // ------------------------------------------------------------------ data flow
  const dirty = () => { note.textContent = "Inputs changed: press Compute frames to recompute (the plotted frames are from the previous inputs)."; };
  const frameKey = (F, I) => `${F}|${I}`;

  async function compute() {
    const mine = ++st.seq;
    const n_runs = S.n_runs.value, t_end = S.t_end.value * 1e-9, sd = Number(seed.value) || 0;
    st.frames = {}; st.nDone = 0;
    note.textContent = "Frames are running on the job pool; points appear as they finish.";
    status.textContent = "computing frames";
    chartB.render();
    try {
      const resp = await api.postJSON("/api/explore/sde/frames", { I_grid: meta.I_grid, F_pumps: meta.F_pumps, n_runs, t_end, seed: sd });
      st.nTotal = resp.frames.length;
      const liPromise = api.postJSON("/api/explore/sde/li", {}).then((r) => jobResult(api, r)).then((li) => { if (mine === st.seq && host.alive()) { st.li = li; chartA.render(); select(); } });
      traceJob.run();
      await Promise.all([liPromise, ...resp.frames.map(async (f) => {
        const res = await jobResult(api, f);
        if (mine !== st.seq || !host.alive()) return;
        st.frames[frameKey(res.F_pump, res.I_over_Ith)] = res; st.nDone++;
        status.textContent = `frames ${st.nDone}/${st.nTotal}`;
        chartB.render(); select();
      })]);
      if (mine === st.seq && host.alive()) { status.textContent = ""; note.textContent = `${st.nTotal} frames ready: ${S.n_runs.value} runs, ${fmt(S.t_end.value)} ns, seed ${sd}. Frames are cached; the same inputs return instantly.`; }
    } catch (e) { if (mine === st.seq && host.alive()) { status.textContent = ""; note.textContent = `Frames failed: ${e.message}`; } }
  }

  const traceJob = latest(async (my, still) => {
    try {
      const resp = await api.postJSON("/api/explore/sde/trace", { I: S.I.value, F_pump: st.F, t_end: S.t_end.value * 1e-9, seed: Number(seed.value) || 0 });
      const res = await jobResult(api, resp);
      if (!still() || !host.alive()) return;
      st.trace = res; chartC.render(); select();
    } catch (e) { if (still() && host.alive()) note.textContent = `Trace failed: ${e.message}`; }
  }, 300);

  function select() {
    const I = S.I.value, f = st.frames[frameKey(st.F, I)];
    const tag = TG.result;
    if (f) {
      ro.g2.update({ value: f.g2_mean, tag, caption: `± ${fmt(f.g2_se)} SE · ${f.n_runs} runs · seed ${f.seed} · F = ${fmt(st.F)}`, markChange: false });
      ro.S.update({ value: f.mean_S, tag, caption: "mean photon number over the runs", markChange: false });
    } else {
      ro.g2.update({ pending: true }); ro.S.update({ pending: true });
    }
    if (st.li) {
      const k = st.li.I.indexOf(I);
      ro.Ith.update({ value: st.li.I_th, tag, caption: "find_threshold: linear-fit extrapolation of the L-I curve", markChange: false });
      ro.Sss.update({ value: k >= 0 ? st.li.S[k] : null, tag, reason: "not on the L-I grid", caption: k >= 0 ? `steady_state residual ${fmt(st.li.residual[k])}` : "", markChange: false });
    }
    const cp = card.card_points.find((p) => p.pump === (st.F === 1 ? "normal" : "quiet"));
    if (cp) ro.card.update({ value: cp.g2, tag: card.card_tag, caption: `± ${fmt(cp.err)} · ${cp.pump} pump · cards/zhao.yaml`, markChange: false });
    const cmp = card.comparison.find((r) => r.pump === (st.F === 1 ? "normal" : "quiet"));
    if (cmp) ro.sigma.update({ value: cmp.sigma_ratio, unit: "σ", tag: TG.model, caption: `verdict ${cmp.verdict} (out/zhao/zhao_fit_comparison.csv); sim ${fmt(cmp.sim_g2)} vs published ${fmt(cmp.published_g2)}`, markChange: false });
    if (st.trace) {
      ro.tg2.update({ value: st.trace.g2_0, tag, caption: `single trace, burn-in ${fmt(st.trace.burnin_frac * 100)}% discarded; g²(0) is unchanged by collection loss`, markChange: false });
      ro.fano.update({ value: st.trace.fano, tag, caption: `thinned by eta_ext = ${fmt(st.trace.eta_ext)}: ${fmt(st.trace.fano_ext)} (loading.f8b_thin_fano)`, markChange: false });
    }
    clear(src).append(sourceLine(meta.calls, "Tuned Tier-2 parameter set, all [E] except T [V]."));
    chartA.render(); chartC.render();
  }

  select();
  compute();
}

// ---------------------------------------------------------------------- charts
function buildLI(st, meta, I) {
  const li = st.li;
  if (!li) return null;
  const c = themeColors();
  const k = li.I.indexOf(I);
  const data = [{ x: li.I, y: li.S, type: "scatter", mode: "lines+markers", name: "S(I), noise-free", line: { color: c.series[0], width: 2 }, marker: { size: 5 }, hovertemplate: "S %{y:.4g} at I/I_th %{x} [E]<extra></extra>" }];
  if (k >= 0) data.push({ x: [I], y: [li.S[k]], type: "scatter", mode: "markers", name: "selected frame", marker: { color: c.beam, size: 10, line: { color: c.ink1, width: 1 } }, hoverinfo: "skip" });
  return { data, layout: baseLayout({ margin: { l: 62, r: 18, t: 30, b: 50 }, legend: { orientation: "h", x: 0, y: -0.26, yanchor: "top", font: { size: 12 } },
    xaxis: { ...baseLayout().xaxis, title: { text: "I / I_th (sde axis, not amperes)", standoff: 6 } },
    yaxis: { ...baseLayout().yaxis, title: { text: "photon number S", standoff: 6 }, rangemode: "tozero" },
    shapes: [{ type: "line", xref: "x", yref: "paper", x0: li.I_th, x1: li.I_th, y0: 0, y1: 1, line: { color: c.ink2, width: 1, dash: "dash" } }],
    annotations: [{ x: li.I_th, y: 1, yref: "paper", yanchor: "bottom", showarrow: false, text: `I_th = ${fmt(li.I_th)} [E]`, font: { size: 11, color: c.ink2 } }] }) };
}

function frameSeries(st, F) {
  return Object.values(st.frames).filter((f) => f.F_pump === F).sort((a, b) => a.I_over_Ith - b.I_over_Ith);
}
function buildG2(st, card, meta) {
  const c = themeColors();
  const data = [];
  // Markers are dodged horizontally (the model cannot separate normal from quiet, so they would overlap);
  // customdata carries the true I/I_th for the hover. Markers are deliberately NOT joined by a line.
  const DX = { sim1: -0.12, sim008: 0.12, run1: -0.06, run008: 0.06 };
  const style = { 1: { name: "simulation, normal F = 1", color: c.series[0], symbol: "circle", dx: DX.sim1 }, 0.08: { name: "simulation, quiet F = 0.08", color: c.series[1], symbol: "square", dx: DX.sim008 } };
  const offscale = [];
  for (const F of meta.F_pumps) {
    const fr = frameSeries(st, F);
    if (!fr.length) continue;
    const sty = style[F] || { name: `simulation, F = ${F}`, color: c.series[2], symbol: "circle", dx: 0 };
    data.push({ x: fr.map((f) => f.I_over_Ith + sty.dx), customdata: fr.map((f) => f.I_over_Ith), y: fr.map((f) => f.g2_mean),
      error_y: { type: "data", array: fr.map((f) => f.g2_se), visible: true, thickness: 1.5, width: 4, color: sty.color },
      type: "scatter", mode: "markers", name: sty.name, marker: { color: sty.color, symbol: sty.symbol, size: 9 },
      hovertemplate: `I/I_th %{customdata}: g² %{y:.4f} ± SE, F = ${F} [E]<extra></extra>` });
    if (st.zoom === "near") for (const f of fr) if (f.g2_mean > 1.1) offscale.push(`${fmt(f.I_over_Ith)}: ${fmt(f.g2_mean)}`);
  }
  const names = { normal: "Zhao card, normal", quiet: "Zhao card, quiet" };
  for (const p of card.card_points) {
    data.push({ x: [p.I_over_Ith], customdata: [p.I_over_Ith], y: [p.g2], error_y: { type: "data", array: [p.err], visible: true, thickness: 2, width: 6, color: c.ink1 },
      type: "scatter", mode: "markers", name: `${names[p.pump]} [${card.card_tag}]`, marker: { color: c.ink1, symbol: p.pump === "quiet" ? "x" : "cross", size: 12, line: { color: c.ink1, width: 2 } },
      hovertemplate: `I/I_th %{customdata}: ${p.g2} ± ${p.err} [${card.card_tag}], ${p.pump}<extra></extra>` });
  }
  for (const [key, col, label, dx, color] of [["g2_normal", "se_normal", "committed fit run, normal", DX.run1, c.series[0]], ["g2_quiet", "se_quiet", "committed fit run, quiet", DX.run008, c.series[1]]]) {
    data.push({ x: card.curve.map((r) => r.I_over_Ith + dx), customdata: card.curve.map((r) => r.I_over_Ith), y: card.curve.map((r) => r[key]),
      error_y: { type: "data", array: card.curve.map((r) => r[col]), visible: true, thickness: 1, width: 3, color },
      type: "scatter", mode: "markers", name: `${label} (out/zhao) [E]`, marker: { color, line: { color, width: 2 }, symbol: "diamond-open", size: 11 },
      hovertemplate: `I/I_th %{customdata}: ${label} %{y:.4f} [E]<extra></extra>` });
  }
  const layout = baseLayout({ hovermode: "closest", margin: { l: 62, r: 18, t: 34, b: 70 }, legend: { orientation: "h", x: 0, y: -0.22, yanchor: "top", font: { size: 12 } },
    xaxis: { ...baseLayout().xaxis, title: { text: "I / I_th (sde axis, not amperes; markers offset by up to 0.12 for legibility)", standoff: 6 }, range: [1.0, 6.4] },
    yaxis: { ...baseLayout().yaxis, title: { text: "g²(0)", standoff: 6 }, ...(st.zoom === "near" ? { range: [0.95, 1.1] } : {}) },
    shapes: [{ type: "line", xref: "paper", yref: "y", x0: 0, x1: 1, y0: 1, y1: 1, line: { color: c.ref, width: 1, dash: "dot" } }],
    annotations: [{ xref: "paper", x: 0.01, xanchor: "left", y: 1, yref: "y", yshift: 10, showarrow: false, text: "g² = 1 (coherent)", font: { size: 11, color: c.ink2 } },
      ...(offscale.length ? [{ xref: "paper", yref: "paper", x: 0.99, xanchor: "right", y: 0.97, yanchor: "top", showarrow: false, text: `off scale (I/I_th: g²): ${offscale.join(", ")}`, font: { size: 11, color: c.muted } }] : [])] });
  return { data, layout };
}
function tableG2(st) {
  const rows = Object.values(st.frames).sort((a, b) => a.F_pump - b.F_pump || a.I_over_Ith - b.I_over_Ith).map((f) => ({ F: f.F_pump, I: f.I_over_Ith, g: f.g2_mean, se: f.g2_se, n: f.n_runs }));
  return rows.length ? { columns: [{ key: "F", label: "F_pump" }, { key: "I", label: "I/I_th" }, { key: "g", label: "g²(0)" }, { key: "se", label: "± SE" }, { key: "n", label: "runs" }], rows } : null;
}

function buildTrace(st) {
  const tr = st.trace;
  if (!tr) return null;
  const c = themeColors();
  return { data: [{ x: tr.t_ns, y: tr.S, type: "scatter", mode: "lines", name: "S(t)", line: { color: c.series[0], width: 1.5 }, hovertemplate: "S %{y:.4g} [E]<extra></extra>" }],
    layout: baseLayout({ margin: { l: 62, r: 18, t: 30, b: 50 }, showlegend: false,
      xaxis: { ...baseLayout().xaxis, title: { text: "t (ns)", standoff: 6 } },
      yaxis: { ...baseLayout().yaxis, title: { text: "photon number S", standoff: 6 } },
      shapes: [{ type: "rect", xref: "x", yref: "paper", x0: tr.t_ns[0], x1: tr.burn_t_ns ?? tr.t_ns[0], y0: 0, y1: 1, fillcolor: rgba(c.ink2 && c.ink2.startsWith("#") ? c.ink2 : "#808080", 0.15), line: { width: 0 }, layer: "below" }],
      annotations: [{ x: (tr.burn_t_ns ?? 0) / 2, y: 1, yref: "paper", yanchor: "bottom", showarrow: false, text: `burn-in ${fmt(tr.burnin_frac * 100)}% discarded`, font: { size: 11, color: c.ink2 } }] }) };
}
