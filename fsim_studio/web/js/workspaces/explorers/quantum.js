// Explorer: Lindblad quantum cross-check (#/explain/lindblad). Live (brief 6: 0.1-1 ms per call).
//   incoherent  3-level Lindblad g2(tau) on top of cw_g2.g2_cw (must agree to 1e-10)
//   rabi        coherent resonance fluorescence, g2(0) = 0 and Rabi ringing
//   hom         HOM indistinguishability against pure dephasing gamma*
//   filter      cascaded Lorentzian filter, filtered g2(0) against filter width
// No physics here: every curve and number comes from GET /api/explore/lindblad?panel=...
import { h, clear } from "../../ui/dom.js";
import { tagChip, tagDash, tagWidth } from "../../ui/tagChip.js";
import { segmented, switchToggle, errorState } from "../../ui/controls.js";
import { fmt } from "../../ui/format.js";
import { chartCard } from "../../charts/chartCard.js";
import { baseLayout, colors as themeColors } from "../../charts/theme.js";
import { loadControls, ctlSlider, readoutGroup, caveatBox, liveFetch, qs, sourceLine } from "./commonA.js";

const PANELS = [
  { id: "incoherent", label: "Incoherent 3-level", sub: "vs cw_g2", usesT: true,
    ids: ["r", "gamma_X", "eps", "pump_ratio", "tau_max", "T"],
    readouts: [["g2_0_lindblad", "g²(0) Lindblad"], ["g2_0_cw_g2", "g²(0) cw_g2"], ["max_rel_diff", "max relative difference"], ["detected_rate", "detected rate"]],
    title: "g²(τ): Lindblad against the rate equations" },
  { id: "rabi", label: "Resonance fluorescence", sub: "coherent Rabi drive", usesT: false,
    ids: ["omega", "gamma", "deph", "tau_max"],
    readouts: [["rho_ee", "ρ_ee (steady state)"], ["g2_0", "g²(0)"], ["g2_max", "max g²(τ)"], ["tau_at_max", "τ at maximum"], ["detected_rate", "detected rate"]],
    title: "g²(τ) under coherent drive: antibunching and Rabi ringing" },
  { id: "hom", label: "Indistinguishability", sub: "HOM", usesT: false,
    ids: ["gamma", "gstar"],
    readouts: [["indistinguishability", "HOM indistinguishability I"]],
    title: "HOM indistinguishability against pure dephasing γ*" },
  { id: "filter", label: "Cascaded filter", sub: "filtered g²(0)", usesT: true,
    ids: ["r", "gamma_X", "L", "w", "T"],
    readouts: [["g2_filtered_0", "filtered g²(0)"], ["filtered_flux", "filtered output flux"]],
    title: "Filtered g²(0) against filter FWHM w" },
];

export async function mount(body, host) {
  let all;
  try { all = await loadControls(host.ctx.api); } catch (e) {
    if (host.alive()) clear(body).append(errorState({ title: "The Lindblad explorer did not answer", message: e.message }));
    return;
  }
  if (!host.alive()) return;
  const caveat = caveatBox([all.caveats.lindblad[0], all.caveats.lindblad[1]]);
  const seg = segmented({ label: "Lindblad panel", value: "incoherent", options: PANELS.map((p) => ({ value: p.id, label: p.label, sub: p.sub })), onChange: (v) => open(v) });
  const stage = h("div.xa-stage");
  body.append(h("div.xa-quantum", h("div.xa-panelbar", seg.el), caveat.el, stage));
  let current = null;

  function open(id) {
    if (current) current.dispose();
    clear(stage);
    current = panel(PANELS.find((p) => p.id === id), all.controls.lindblad[id], stage, host, caveat);
  }
  open("incoherent");
}

function panel(def, specs, stage, host, caveat) {
  const st = { thermal: false, last: null };
  let dead = false;
  const by = Object.fromEntries(specs.map((s) => [s.id, s]));
  const key = `lb-${def.id}`;
  const chart = chartCard({ title: def.title, height: 440, build: () => build(def, st.last), table: () => table(def, st.last), filename: () => `fsim-lindblad-${def.id}` });
  host.charts[key] = chart;

  const status = h("p.xp-status", { "aria-live": "polite" });
  const S = {};
  const go = liveFetch(host, () => {
    const q = { panel: def.id };
    for (const id of def.ids) if (id !== "T") q[id] = S[id].value;
    if (def.usesT && st.thermal) q.T = S.T.value;
    return `/api/explore/lindblad?${qs(q)}`;
  }, render, { ms: 100, status });
  for (const id of def.ids) S[id] = ctlSlider(`xa-lb-${def.id}`, by[id], () => go());
  if (def.usesT) S.T.disable(true);
  const rows = def.ids.filter((id) => id !== "T").map((id) => S[id].row);
  if (def.usesT) {
    const sw = switchToggle({ label: "thermal escape at T (V-a fit)", checked: false, onChange: (v) => { st.thermal = v; S.T.disable(!v); go(); } });
    rows.push(h("div.xs-row", sw), S.T.row);
  }
  const { el: roEl, ro } = readoutGroup(def.readouts);
  const stateBox = h("p.xa-check");
  const srcBox = h("div");
  const root = h("div.xp-split", { dataset: { panel: def.id } },
    h("div.xp-main", chart.el, h("p.xp-cap", capFor(def.id))),
    h("aside.xp-side.panel", h("header.panel-head", h("h2.engrave", "Inputs"), status), h("div.xp-ctls", rows), roEl, stateBox, srcBox));
  stage.append(root);

  function render(r, err) {
    if (dead) return;
    if (!r) {
      st.last = null;
      chart.setEmpty(errorState({ title: "The Lindblad explorer did not answer", message: err?.message || "" }));
      chart.render();
      return;
    }
    chart.setEmpty(null);
    st.last = r;
    const t = r.tags;
    for (const [k] of def.readouts) {
      const v = r.readouts[k];
      const tag = k === "max_rel_diff" ? t.max_rel_diff : t.readouts;
      ro[k].update({ value: v, tag, caption: capRead(def.id, k, r), markChange: false });
    }
    clear(stateBox);
    if (r.state) stateBox.append(h("span", `state: |Tr ρ − 1| ${fmt(r.state.trace_error)} · min eigenvalue ${fmt(r.state.min_eig)} · Hermiticity ${fmt(r.state.hermiticity)} `), tagChip("DR", { size: "sm" }), h("span", r.state.ok ? " physical" : " NOT physical"));
    stateBox.dataset.ok = String(r.state ? r.state.ok : true);
    chart.setTag(t.curves);
    chart.render();
    clear(srcBox).append(sourceLine(r));
  }
  go();
  return { dispose() { dead = true; try { chart.dispose(); } catch { /* */ } delete host.charts[key]; } };
}

function capFor(id) {
  return {
    incoherent: "With no coherent drive and incoherent pumping, the 3-level Lindblad model must reproduce the rate-equation g²(τ) of cw_g2; the two curves lie on top of each other and their maximum relative difference is the readout.",
    rabi: "A resonantly driven two-level emitter: g²(0) = 0, then damped Rabi oscillation. The rate-equation modules refuse this regime.",
    hom: "Two-photon interference visibility of a single-excitation emitter: pure dephasing γ* lowers it. The marker is the chosen γ*.",
    filter: "A cascaded Lorentzian filter acting on a two-level emitter. If the requested line FWHM L is below the lifetime-limited width, the line is lifetime-limited and L has no effect.",
  }[id];
}
function capRead(id, k, r) {
  if (k === "max_rel_diff") return `tolerance ${fmt(r.tolerance)}: ${r.agrees ? "agrees" : "DISAGREES"}`;
  if (id === "rabi" && k === "rho_ee") return `γ₂ = ${fmt(r.inputs.gamma_2)} 1/ns`;
  if (id === "filter" && k === "g2_filtered_0") return r.inputs.lifetime_limited ? "line is lifetime-limited (L below the population-relaxation width)" : `pure dephasing ${fmt(r.inputs.dephasing_rate)} 1/ns`;
  return "";
}

function build(def, r) {
  if (!r) return null;
  const c = themeColors();
  const tag = r.tags.curves, w = tagWidth(tag), dash = tagDash(tag);
  const L = baseLayout();
  const base = { hovermode: "x unified", margin: { l: 58, r: 18, t: 30, b: 62 }, legend: { orientation: "h", x: 0, y: -0.2, yanchor: "top", font: { size: 12 } } };
  if (def.id === "incoherent") {
    return {
      data: [
        { x: r.curves.tau, y: r.curves.g2_cw_g2, type: "scatter", mode: "lines", name: "cw_g2.g2_cw (rate equations)", line: { color: c.series[0], width: w + 2, dash }, hovertemplate: `rate eq. %{y:.5f} [${tag}]<extra></extra>` },
        { x: r.curves.tau, y: r.curves.g2_lindblad, type: "scatter", mode: "lines", name: "lindblad.g2_tau", line: { color: c.series[1], width: 1.5, dash: "solid" }, hovertemplate: `Lindblad %{y:.5f} [${tag}]<extra></extra>` },
      ],
      layout: baseLayout({ ...base, xaxis: { ...L.xaxis, title: { text: "τ (ns)", standoff: 6 } }, yaxis: { ...L.yaxis, title: { text: "g²(τ)", standoff: 6 }, rangemode: "tozero" } }),
    };
  }
  if (def.id === "rabi") {
    return {
      data: [{ x: r.curves.tau, y: r.curves.g2, type: "scatter", mode: "lines", name: "g²(τ)", line: { color: c.series[0], width: w, dash }, hovertemplate: `g² %{y:.4f} [${tag}]<extra></extra>` },
        { x: [r.readouts.tau_at_max], y: [r.readouts.g2_max], type: "scatter", mode: "markers", name: "maximum", marker: { color: c.beam, size: 9, symbol: "diamond" }, hoverinfo: "skip" }],
      layout: baseLayout({ ...base, shapes: [{ type: "line", xref: "paper", x0: 0, x1: 1, yref: "y", y0: 1, y1: 1, line: { color: c.ref, width: 1, dash: "dot" } }],
        xaxis: { ...L.xaxis, title: { text: "τ (ns)", standoff: 6 } }, yaxis: { ...L.yaxis, title: { text: "g²(τ)", standoff: 6 }, rangemode: "tozero" } }),
    };
  }
  if (def.id === "hom") {
    return {
      data: [{ x: r.curves.gstar, y: r.curves.indistinguishability, type: "scatter", mode: "lines", name: "I(γ*)", line: { color: c.series[0], width: w, dash }, hovertemplate: `I %{y:.4f} [${tag}]<extra></extra>` },
        { x: [r.inputs.gstar > 0 ? r.inputs.gstar : null], y: [r.readouts.indistinguishability], type: "scatter", mode: "markers", name: "chosen γ*", marker: { color: c.beam, size: 9, symbol: "diamond" }, hoverinfo: "skip" }],
      layout: baseLayout({ ...base, xaxis: { ...L.xaxis, type: "log", title: { text: "pure dephasing γ* (1/ns)", standoff: 6 } }, yaxis: { ...L.yaxis, title: { text: "indistinguishability I", standoff: 6 }, range: [0, 1.02] } }),
    };
  }
  return {
    data: [{ x: r.curves.w, y: r.curves.g2_filtered_0, type: "scatter", mode: "lines", name: "filtered g²(0)", line: { color: c.series[0], width: w, dash }, hovertemplate: `g²_f(0) %{y:.4f} [${tag}]<extra></extra>` },
      { x: [r.inputs.w], y: [r.readouts.g2_filtered_0], type: "scatter", mode: "markers", name: "chosen w", marker: { color: c.beam, size: 9, symbol: "diamond" }, hoverinfo: "skip" }],
    layout: baseLayout({ ...base, shapes: [{ type: "line", xref: "paper", x0: 0, x1: 1, yref: "y", y0: 1, y1: 1, line: { color: c.ref, width: 1, dash: "dot" } }],
      xaxis: { ...L.xaxis, type: "log", title: { text: "filter FWHM w (1/ns)", standoff: 6 } }, yaxis: { ...L.yaxis, title: { text: "filtered g²(0)", standoff: 6 }, rangemode: "tozero" } }),
  };
}

function table(def, r) {
  if (!r) return null;
  const cv = r.curves;
  if (def.id === "incoherent") return { columns: [{ key: "t", label: "τ (ns)" }, { key: "a", label: "Lindblad" }, { key: "b", label: "cw_g2" }], rows: cv.tau.map((x, i) => ({ t: x, a: cv.g2_lindblad[i], b: cv.g2_cw_g2[i] })) };
  if (def.id === "rabi") return { columns: [{ key: "t", label: "τ (ns)" }, { key: "g", label: "g²" }], rows: cv.tau.map((x, i) => ({ t: x, g: cv.g2[i] })) };
  if (def.id === "hom") return { columns: [{ key: "t", label: "γ* (1/ns)" }, { key: "g", label: "I" }], rows: cv.gstar.map((x, i) => ({ t: x, g: cv.indistinguishability[i] })) };
  return { columns: [{ key: "t", label: "w (1/ns)" }, { key: "g", label: "filtered g²(0)" }], rows: cv.w.map((x, i) => ({ t: x, g: cv.g2_filtered_0[i] })) };
}
