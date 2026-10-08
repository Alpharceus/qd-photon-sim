// Explorer: pulsed g2 from periodic photon counting (#/explain/pulse-counting). Fully live (brief 2(e)).
//   Plot A  g2, g2_adj vs r with loading.f1b_g2(mu = r tau_on, eps) as the instantaneous-pulse approximation
//   Plot B  mean counts per period vs r
//   Plot C  deterministic-load cycle: g2, mean counts, blocked-load probability vs period
// No physics here: every curve and number comes from GET /api/explore/pulse_counting.
import { h, clear } from "../../ui/dom.js";
import { tagChip, tagDash, tagWidth } from "../../ui/tagChip.js";
import { segmented, switchToggle, errorState } from "../../ui/controls.js";
import { fmt, isNum } from "../../ui/format.js";
import { chartCard } from "../../charts/chartCard.js";
import { baseLayout, colors as themeColors } from "../../charts/theme.js";
import { loadControls, ctlSlider, readoutGroup, caveatBox, liveFetch, qs, sourceLine } from "./commonA.js";

export async function mount(body, host) {
  let specs;
  try { specs = (await loadControls(host.ctx.api)).controls.pulse_counting; } catch (e) {
    if (host.alive()) clear(body).append(errorState({ title: "The pulse-counting explorer did not answer", message: e.message }));
    return;
  }
  if (!host.alive()) return;
  const st = { gate: "none", thermal: false, last: null };
  const by = Object.fromEntries(specs.map((s) => [s.id, s]));

  const mk = (key, title, build, table, file, height) => {
    const c = chartCard({ title, height, build: () => build(st.last), table: () => table(st.last), filename: () => file });
    host.charts[key] = c;
    return c;
  };
  const chartA = mk("pcA", "Pulsed g² against pump rate r", buildA, tableA, "fsim-pulse-counting-g2", 380);
  const chartB = mk("pcB", "Mean counts per period against r", buildB, tableB, "fsim-pulse-counting-counts", 260);
  const chartC = mk("pcC", "Deterministic load cycle against period", buildC, tableC, "fsim-deterministic-cycle", 420);

  const caveat = caveatBox([]);
  const status = h("p.xp-status", { "aria-live": "polite" });
  const S = {};
  const go = liveFetch(host, () => `/api/explore/pulse_counting?${qs({
    r: S.r.value, tau_on: S.tau_on.value, period: S.period.value, gamma_X: S.gamma_X.value, eps: S.eps.value,
    pump_ratio: S.pump_ratio.value, eta_load: S.eta_load.value, gate: st.gate, T: st.thermal ? S.T.value : null })}`,
  render, { ms: 100, status });
  for (const id of ["r", "tau_on", "period", "gamma_X", "eps", "pump_ratio", "eta_load", "T"]) S[id] = ctlSlider("xa-pc", by[id], () => go());
  S.T.disable(true);
  const gateSeg = segmented({ label: "Counting gate", value: "none", options: [{ value: "none", label: "whole period" }, { value: "auto", label: "τ_on + 5/γ_X" }], onChange: (v) => { st.gate = v; go(); } });
  const thermalSw = switchToggle({ label: "thermal escape at T (V-a fit)", checked: false, onChange: (v) => { st.thermal = v; S.T.disable(!v); go(); } });
  const { el: roEl, ro } = readoutGroup([
    ["g2", "g² (long-delay norm.)"], ["g2_adj", "g²_adj (adjacent norm.)"], ["mean_counts", "mean counts / period"],
    ["adj_factor", "adjacent peak factor"], ["f1b", "f1b_g2(μ = r τ_on, ε)"],
    ["d_g2", "deterministic g²"], ["d_mean", "deterministic mean counts"], ["d_block", "blocked-load probability"], ["d_pair", "one-pair valid"]]);
  const srcBox = h("div");

  body.append(h("div.xp-split",
    h("div.xp-main", caveat.el, chartA.el, chartB.el, chartC.el,
      h("p.xp-cap", "Plot A: the dashed instantaneous-pulse curve is the static f1b model; the gap to the solid pulsed g² is re-excitation during the pulse. The 0.5 line is the acceptance gate from /api/gates.")),
    h("aside.xp-side.panel",
      h("header.panel-head", h("h2.engrave", "Inputs"), status),
      h("div.xp-ctls", S.r.row, S.tau_on.row, S.period.row, S.gamma_X.row, S.eps.row, S.pump_ratio.row, S.eta_load.row,
        h("div.xs-row", h("span.xs-label", "Counting gate ", tagChip("A", { size: "sm", title: "pulse_counting.py:118-120 [A]" })), gateSeg.el),
        h("div.xs-row", thermalSw), S.T.row),
      roEl, srcBox)));

  function render(r, err) {
    if (!r) {
      st.last = null;
      chartA.setEmpty(errorState({ title: "The pulse-counting explorer did not answer", message: err?.message || "" }));
      chartA.render(); chartB.render(); chartC.render();
      return;
    }
    chartA.setEmpty(null);
    st.last = r;
    const t = r.tags, p = r.point, d = r.deterministic;
    caveat.set([r.caveat, r.caveat_extra]);
    const at = `at r = ${fmt(r.inputs.r)}/ns, μ = ${fmt(r.mu)}`;
    ro.g2.update({ value: p.g2, tag: t.point, caption: `${at}${p.converged ? "" : " · NOT converged"}`, markChange: false });
    ro.g2_adj.update({ value: p.g2_adj, tag: t.point, caption: "g² divided by the adjacent-peak factor", markChange: false });
    ro.mean_counts.update({ value: p.mean_counts, tag: t.point, caption: `X ${fmt(p.mean_counts_x)} · XX ${fmt(p.mean_counts_xx)}`, markChange: false });
    ro.adj_factor.update({ value: p.adjacent_peak_factor, tag: t.point, caption: "⟨m_n m_n+1⟩ / ⟨m⟩²", markChange: false });
    ro.f1b.update({ value: p.f1b_g2_inst, tag: t.f1b, caption: "instantaneous-pulse approximation, same μ", markChange: false });
    ro.d_g2.update({ value: d.g2, tag: t.deterministic, caption: d.invalid_reason ? `invalid: ${d.invalid_reason}` : `period ${fmt(r.inputs.period)} ns`, reason: d.invalid_reason, markChange: false });
    ro.d_mean.update({ value: d.mean_counts, tag: t.deterministic, caption: `gate ${fmt(d.gate_ns_used)} ns`, markChange: false });
    ro.d_block.update({ value: d.blocked_load_probability, tag: t.deterministic, caption: "residual X blocks the next load", markChange: false });
    ro.d_pair.update({ value: d.one_pair_valid ? 1 : 0, tag: t.deterministic, caption: d.one_pair_valid ? "true: η_load = 1 and no blocking" : "false: η_load < 1 or loads blocked", markChange: false });
    chartA.setTag(t.sweep); chartB.setTag(t.sweep); chartC.setTag(t.deterministic);
    chartA.setTitle("Pulsed g² against pump rate r", `τ_on ${fmt(r.inputs.tau_on)} ns · period ${fmt(r.inputs.period)} ns · ${fmt(r.elapsed_ms)} ms`);
    chartB.setTitle("Mean counts per period against r", "dot-only, no background");
    chartC.setTitle("Deterministic load cycle against period", `η_load ${fmt(r.inputs.eta_load)}`);
    chartA.render(); chartB.render(); chartC.render();
    clear(srcBox).append(sourceLine(r));
  }
  go();
}

function buildA(r) {
  if (!r) return null;
  const c = themeColors();
  const tag = r.tags.sweep, v = r.vs_r;
  const w = tagWidth(tag), dash = tagDash(tag);
  const data = [
    { x: v.r, y: v.g2, type: "scatter", mode: "lines", name: "g² (pulsed, finite pulse)", line: { color: c.series[0], width: w, dash }, connectgaps: false, hovertemplate: `g² %{y:.4f} [${tag}]<extra></extra>` },
    { x: v.r, y: v.g2_adj, type: "scatter", mode: "lines", name: "g²_adj", line: { color: c.series[1], width: 1.5, dash: "solid" }, connectgaps: false, hovertemplate: `g²_adj %{y:.4f} [${tag}]<extra></extra>` },
    { x: v.r, y: v.f1b_g2, type: "scatter", mode: "lines", name: "f1b_g2(μ = r τ_on): instantaneous-pulse approximation", line: { color: c.series[2], width: tagWidth(r.tags.f1b), dash: tagDash(r.tags.f1b) }, hovertemplate: `f1b %{y:.4f} [${r.tags.f1b}]<extra></extra>` },
    { x: [r.inputs.r], y: [r.point.g2], type: "scatter", mode: "markers", name: "this r", marker: { color: c.beam, size: 9, symbol: "diamond" }, hoverinfo: "skip" },
  ];
  const shapes = [], annotations = [];
  const g = r.gate_g2;
  if (g && isNum(g.value)) {
    shapes.push({ type: "line", xref: "paper", x0: 0, x1: 1, yref: "y", y0: g.value, y1: g.value, line: { color: c.ink1, width: 1.2, dash: "dash" } });
    annotations.push({ xref: "paper", x: 0, xanchor: "left", yref: "y", y: Math.log10(g.value), yanchor: "bottom", text: `gate ${fmt(g.value)} [${g.tag}] (/api/gates)`, showarrow: false, font: { size: 11, color: c.ink1 } });
  }
  const layout = baseLayout({
    hovermode: "x unified", margin: { l: 58, r: 18, t: 30, b: 62 },
    legend: { orientation: "h", x: 0, y: -0.22, yanchor: "top", font: { size: 12 } },
    xaxis: { ...baseLayout().xaxis, type: "log", title: { text: "pump rate r (1/ns)", standoff: 6 } },
    yaxis: { ...baseLayout().yaxis, type: "log", title: { text: "g²", standoff: 6 } },
    shapes, annotations,
  });
  return { data, layout };
}

function buildB(r) {
  if (!r) return null;
  const c = themeColors();
  const tag = r.tags.sweep, v = r.vs_r;
  const data = [{ x: v.r, y: v.mean_counts, type: "scatter", mode: "lines", name: "mean counts / period", line: { color: c.series[0], width: tagWidth(tag), dash: tagDash(tag) }, hovertemplate: `⟨m⟩ %{y:.4g} [${tag}]<extra></extra>` },
    { x: [r.inputs.r], y: [r.point.mean_counts], type: "scatter", mode: "markers", name: "this r", marker: { color: c.beam, size: 9, symbol: "diamond" }, hoverinfo: "skip" }];
  const layout = baseLayout({
    hovermode: "x unified", margin: { l: 58, r: 18, t: 30, b: 46 },
    xaxis: { ...baseLayout().xaxis, type: "log", title: { text: "pump rate r (1/ns)", standoff: 6 } },
    yaxis: { ...baseLayout().yaxis, type: "log", title: { text: "⟨m⟩", standoff: 6 } },
  });
  return { data, layout };
}

function buildC(r) {
  if (!r) return null;
  const c = themeColors();
  const tag = r.tags.deterministic, v = r.vs_period;
  const w = tagWidth(tag), dash = tagDash(tag);
  const mark = (yaxis, y, name) => ({ x: [r.inputs.period], y: [y], type: "scatter", mode: "markers", name, yaxis, marker: { color: c.beam, size: 8, symbol: "diamond" }, showlegend: false, hoverinfo: "skip" });
  const d = r.deterministic;
  const data = [
    { x: v.period, y: v.g2, type: "scatter", mode: "lines", name: "g²", line: { color: c.series[0], width: w, dash }, connectgaps: false, hovertemplate: `g² %{y:.4g} [${tag}]<extra></extra>` },
    { x: v.period, y: v.mean_counts, type: "scatter", mode: "lines", name: "mean counts", yaxis: "y2", line: { color: c.series[1], width: w, dash }, hovertemplate: `⟨m⟩ %{y:.4g} [${tag}]<extra></extra>` },
    { x: v.period, y: v.blocked_load_probability, type: "scatter", mode: "lines", name: "blocked-load probability", yaxis: "y3", line: { color: c.series[2], width: w, dash }, hovertemplate: `blocked %{y:.3g} [${tag}]<extra></extra>` },
    mark("y", d.g2, "g² here"), mark("y2", d.mean_counts, "mean here"), mark("y3", d.blocked_load_probability, "blocked here"),
  ];
  const ax = baseLayout().yaxis;
  const layout = baseLayout({
    hovermode: "x unified", margin: { l: 66, r: 18, t: 30, b: 46 },
    xaxis: { ...baseLayout().xaxis, title: { text: "period (ns)", standoff: 6 }, anchor: "y3" },
    yaxis: { ...ax, domain: [0.70, 1], title: { text: "g²", standoff: 4 } },
    yaxis2: { ...ax, domain: [0.36, 0.64], title: { text: "⟨m⟩", standoff: 4 } },
    yaxis3: { ...ax, domain: [0, 0.30], title: { text: "P(blocked)", standoff: 4 }, exponentformat: "e" },
    showlegend: false,
  });
  return { data, layout };
}

function tableA(r) {
  if (!r) return null;
  const v = r.vs_r;
  return { columns: [{ key: "r", label: "r (1/ns)" }, { key: "g2", label: "g²" }, { key: "adj", label: "g²_adj" }, { key: "f", label: "f1b_g2" }, { key: "m", label: "⟨m⟩" }],
    rows: v.r.map((x, i) => ({ r: x, g2: v.g2[i], adj: v.g2_adj[i], f: v.f1b_g2[i], m: v.mean_counts[i] })) };
}
function tableB(r) {
  if (!r) return null;
  const v = r.vs_r;
  return { columns: [{ key: "r", label: "r (1/ns)" }, { key: "m", label: "⟨m⟩" }], rows: v.r.map((x, i) => ({ r: x, m: v.mean_counts[i] })) };
}
function tableC(r) {
  if (!r) return null;
  const v = r.vs_period;
  return { columns: [{ key: "p", label: "period (ns)" }, { key: "g2", label: "g²" }, { key: "m", label: "⟨m⟩" }, { key: "b", label: "P(blocked)" }],
    rows: v.period.map((x, i) => ({ p: x, g2: v.g2[i], m: v.mean_counts[i], b: v.blocked_load_probability[i] })) };
}
