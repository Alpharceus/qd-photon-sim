// Explorer: CW g2(tau), IRF and background (#/explain/cw-g2). Live sliders (brief 1(e)).
//   Plot A  g2_dot / g2_meas / g2_raw vs tau (fsim_core.cw_g2.cw_report), 1 and gate lines
//   Plot B  g2_dot(0) vs r with the r->0 and r->inf limits (cw_g2.g2_cw_zero)
// No physics here: every curve and number comes from GET /api/explore/cw_g2.
import { h, clear } from "../../ui/dom.js";
import { tagChip, tagDash, tagWidth } from "../../ui/tagChip.js";
import { segmented, switchToggle, errorState } from "../../ui/controls.js";
import { fmt, isNum } from "../../ui/format.js";
import { chartCard } from "../../charts/chartCard.js";
import { baseLayout, colors as themeColors } from "../../charts/theme.js";
import { loadControls, ctlSlider, readoutGroup, caveatBox, liveFetch, qs, sourceLine } from "./commonA.js";

export async function mount(body, host) {
  let specs;
  try { specs = (await loadControls(host.ctx.api)).controls.cw_g2; } catch (e) {
    if (host.alive()) clear(body).append(errorState({ title: "The CW g2 explorer did not answer", message: e.message }));
    return;
  }
  if (!host.alive()) return;
  const st = { shape: "exponential", thermal: false, last: null };
  const by = Object.fromEntries(specs.map((s) => [s.id, s]));

  const chartA = chartCard({
    title: "g²(τ): intrinsic, with background, as the HBT measures it", height: 400,
    build: () => buildA(st.last), table: () => tableA(st.last), filename: () => "fsim-cw-g2-tau",
  });
  const chartB = chartCard({
    title: "Intrinsic g²(0) against pump rate r", height: 300,
    build: () => buildB(st.last), table: () => tableB(st.last), filename: () => "fsim-cw-g2-vs-r",
  });
  host.charts.cwA = chartA;
  host.charts.cwB = chartB;

  const caveat = caveatBox([]);
  const status = h("p.xp-status", { "aria-live": "polite" });
  const S = {};
  const go = liveFetch(host, () => `/api/explore/cw_g2?${qs({
    r: S.r.value, gamma_X: S.gamma_X.value, eps: S.eps.value, rho: S.rho.value, irf_fwhm_ps: S.irf_fwhm_ps.value,
    irf_shape: st.shape, pump_ratio: S.pump_ratio.value, tau_max: S.tau_max.value, T: st.thermal ? S.T.value : null })}`,
  render, { ms: 100, status });

  for (const id of ["r", "gamma_X", "eps", "rho", "irf_fwhm_ps", "pump_ratio", "tau_max", "T"]) S[id] = ctlSlider("xa-cw", by[id], () => go());
  S.T.disable(true);
  const shapeSeg = segmented({ label: "IRF shape", value: st.shape, options: [{ value: "exponential", label: "exponential" }, { value: "gaussian", label: "gaussian" }], onChange: (v) => { st.shape = v; go(); } });
  const thermalSw = switchToggle({ label: "thermal escape at T (V-a fit)", checked: false, onChange: (v) => { st.thermal = v; S.T.disable(!v); go(); } });
  const { el: roEl, ro } = readoutGroup([
    ["g2_dot0", "g²_dot(0)"], ["g2_meas0", "g²_meas(0)"], ["g2_raw0", "g²_raw(0)"],
    ["tau_dip", "τ_dip"], ["A_dip", "A_dip"], ["g2_intrinsic_from_raw", "g² intrinsic from raw"], ["k_X", "k_X (escape)"]]);
  const checksBox = h("div.xa-checks", { "aria-label": "Consistency checks" });
  const srcBox = h("div");

  body.append(h("div.xp-split",
    h("div.xp-main", caveat.el, chartA.el, chartB.el,
      h("p.xp-cap", "Plot A: the dashed dot curve is what the dot does; the background law lifts it to g²_meas; the instrument response (IRF) smears the dip into g²_raw, the dc HBT histogram. Plot B: the closed form for g²_dot(0) with its two pump limits.")),
    h("aside.xp-side.panel",
      h("header.panel-head", h("h2.engrave", "Inputs"), status),
      h("div.xp-ctls", S.r.row, S.gamma_X.row, S.eps.row, S.rho.row, S.irf_fwhm_ps.row,
        h("div.xs-row", h("span.xs-label", "IRF shape ", tagChip("V", { size: "sm", title: "Reischle 2008 exponential IRF" })), shapeSeg.el),
        S.pump_ratio.row, S.tau_max.row, h("div.xs-row", thermalSw), S.T.row),
      roEl, checksBox, srcBox)));

  function render(r, err) {
    if (!r) {
      st.last = null;
      chartA.setEmpty(errorState({ title: "The CW g2 explorer did not answer", message: err?.message || "" }));
      chartA.render(); chartB.render();
      return;
    }
    chartA.setEmpty(null);
    st.last = r;
    const t = r.tags;
    const q = r.readouts;
    caveat.set([r.caveat, r.bunching ? r.caveat_bunching : null]);
    ro.g2_dot0.update({ value: q.g2_dot0, tag: t.g2_dot0, qualifier: r.bunching ? "cascade bunching" : "", caption: "closed form, cw_g2.g2_cw_zero via cw_report", markChange: false });
    ro.g2_meas0.update({ value: q.g2_meas0, tag: t.g2_meas0, caption: "background law 1 − ρ²(1 − g²_dot)", markChange: false });
    ro.g2_raw0.update({ value: q.g2_raw0, tag: t.g2_raw0, caption: `what a ${fmt(r.inputs.irf_fwhm_ps)} ps ${r.inputs.irf_shape} IRF measures`, markChange: false });
    ro.tau_dip.update({ value: q.tau_dip, unit: "ns", tag: t.tau_dip, caption: "single-exponential dip fitted to g²_meas", markChange: false });
    ro.A_dip.update({ value: q.A_dip, tag: t.A_dip, caption: "dip depth", markChange: false });
    ro.g2_intrinsic_from_raw.update({ value: q.g2_intrinsic_from_raw, tag: t.g2_intrinsic_from_raw, caption: "deconvolving the raw value with the fitted τ_dip", markChange: false });
    if (r.inputs.T != null) ro.k_X.update({ value: r.k_X, unit: "1/ns", tag: t.k_X, caption: `k_XX ${fmt(r.k_XX)} · T ${fmt(r.inputs.T)} K`, markChange: false });
    else ro.k_X.update({ value: null, reason: "thermal escape off (k_X = 0)", tag: null });
    clear(checksBox).append(h("p.xa-checks-h", "Consistency"), ...r.checks.map((c) => h("p.xa-check", { dataset: { ok: String(c.ok) } },
      h("span", `${c.label}: `), h("span.xa-num", fmt(c.value)), h("span", ` ≤ ${fmt(c.tol)} `), tagChip(c.tag, { size: "sm" }), h("span", c.ok ? " ok" : " FAIL"))));
    chartA.setTag(t.curves);
    chartA.setTitle("g²(τ): intrinsic, with background, as the HBT measures it",
      `${r.n_tau} grid points${r.decimated ? ` (thinned to ${r.curves.tau.length} for display)` : ""} · ${fmt(r.elapsed_ms)} ms${r.stiff ? " · stiff escape: dense grid" : ""}`);
    chartB.setTag(t.limits);
    chartB.setTitle("Intrinsic g²(0) against pump rate r", "closed form, limits r→0 and r→∞");
    chartA.render(); chartB.render();
    clear(srcBox).append(sourceLine(r));
  }
  go();
}

function buildA(r) {
  if (!r) return null;
  const c = themeColors();
  const tag = r.tags.curves;
  const dash = tagDash(tag), w = tagWidth(tag);
  const data = [
    { x: r.curves.tau, y: r.curves.g2_dot, type: "scatter", mode: "lines", name: "g²_dot (dot)", line: { color: c.series[0], width: w, dash }, hovertemplate: `dot %{y:.3f} [${tag}]<extra></extra>` },
    { x: r.curves.tau, y: r.curves.g2_meas, type: "scatter", mode: "lines", name: "g²_meas (with background)", line: { color: c.series[1], width: w, dash }, hovertemplate: `meas %{y:.3f} [${tag}]<extra></extra>` },
    { x: r.curves.tau, y: r.curves.g2_raw, type: "scatter", mode: "lines", name: "g²_raw (HBT histogram)", line: { color: c.series[2], width: w, dash }, hovertemplate: `raw %{y:.3f} [${tag}]<extra></extra>` },
  ];
  const shapes = [{ type: "line", xref: "paper", x0: 0, x1: 1, yref: "y", y0: 1, y1: 1, line: { color: c.ref, width: 1, dash: "dot" } }];
  const annotations = [{ xref: "paper", x: 1, xanchor: "right", yref: "y", y: 1, yanchor: "bottom", text: "g² = 1", showarrow: false, font: { size: 11, color: c.ink2 } }];
  const g = r.gate_g2;
  if (g && isNum(g.value)) {
    shapes.push({ type: "line", xref: "paper", x0: 0, x1: 1, yref: "y", y0: g.value, y1: g.value, line: { color: c.ink1, width: 1.2, dash: "dash" } });
    annotations.push({ xref: "paper", x: 0, xanchor: "left", yref: "y", y: g.value, yanchor: "bottom", text: `gate ${fmt(g.value)} [${g.tag}] (/api/gates)`, showarrow: false, font: { size: 11, color: c.ink1 } });
  }
  const layout = baseLayout({
    hovermode: "x unified", margin: { l: 58, r: 18, t: 30, b: 46 },
    xaxis: { ...baseLayout().xaxis, title: { text: "τ (ns)", standoff: 6 } },
    yaxis: { ...baseLayout().yaxis, title: { text: "g²(τ)", standoff: 6 }, rangemode: "tozero" },
    shapes, annotations,
  });
  return { data, layout };
}

function buildB(r) {
  if (!r) return null;
  const c = themeColors();
  const v = r.vs_r;
  const tag = r.tags.limits;
  const data = [{ x: v.r, y: v.g2_dot0, type: "scatter", mode: "lines", name: "g²_dot(0)", line: { color: c.series[0], width: tagWidth(r.tags.g2_dot0), dash: tagDash(r.tags.g2_dot0) }, hovertemplate: `g²_dot(0) %{y:.3f} [${r.tags.g2_dot0}]<extra></extra>` },
    { x: [r.inputs.r], y: [r.readouts.g2_dot0], type: "scatter", mode: "markers", name: "this r", marker: { color: c.beam, size: 9, symbol: "diamond" }, hoverinfo: "skip" }];
  const shapes = [], annotations = [];
  const lim = (y, text) => {   // annotations on a log axis take log10 coordinates; shapes take data values
    shapes.push({ type: "line", xref: "paper", x0: 0, x1: 1, yref: "y", y0: y, y1: y, line: { color: c.ref, width: 1, dash: "dot" } });
    annotations.push({ xref: "paper", x: 1, xanchor: "right", yref: "y", y: Math.log10(y), yanchor: "bottom", text: `${text} ${fmt(y)} [${tag}]`, showarrow: false, font: { size: 11, color: c.ink2 } });
  };
  if (isNum(v.low_limit)) lim(v.low_limit, "r→0:");
  if (isNum(v.high_limit)) lim(v.high_limit, "r→∞:");
  const layout = baseLayout({
    hovermode: "x unified", margin: { l: 58, r: 18, t: 30, b: 46 },
    xaxis: { ...baseLayout().xaxis, type: "log", title: { text: "pump rate r (1/ns)", standoff: 6 } },
    yaxis: { ...baseLayout().yaxis, type: "log", title: { text: "g²_dot(0)", standoff: 6 } },
    shapes, annotations,
  });
  return { data, layout };
}

function tableA(r) {
  if (!r) return null;
  const n = r.curves.tau.length, step = Math.ceil(n / 240), rows = [];
  for (let i = 0; i < n; i += step) rows.push({ tau: r.curves.tau[i], dot: r.curves.g2_dot[i], meas: r.curves.g2_meas[i], raw: r.curves.g2_raw[i] });
  return { caption: `every ${step}th of ${n} points`, columns: [{ key: "tau", label: "τ (ns)" }, { key: "dot", label: "g²_dot" }, { key: "meas", label: "g²_meas" }, { key: "raw", label: "g²_raw" }], rows };
}
function tableB(r) {
  if (!r) return null;
  const v = r.vs_r;
  return { columns: [{ key: "r", label: "r (1/ns)" }, { key: "g", label: "g²_dot(0)" }], rows: v.r.map((x, i) => ({ r: x, g: v.g2_dot0[i] })) };
}
