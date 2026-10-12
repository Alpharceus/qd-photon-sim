// Paper-vs-model series overlay (Papers workspace). Paper points are ink markers with error bars
// (digitised [DR] points are OUTLINE markers); the model is the engine's curve at the central value
// of every other leaf, drawn with the provenance dash grammar of its tag. Nothing is computed here:
// x/y come from the card (points) and from the reproduce job (model_series).
import { h } from "../ui/dom.js";
import { chartCard } from "./chartCard.js";
import { baseLayout, colors as themeColors } from "./theme.js";
import { normTag, tagDash, tagWidth, widest } from "../ui/tagChip.js";
import { isNum } from "../ui/format.js";
import { displayUnit, displayNotation } from "../ui/units.js";

// Human axis labels (with units) for the ids the ledgers use; anything else falls back to words, never the raw id.
const OBS_AXIS = {
  g2_0: "g²(0)", g2_XX: "g²(0) of the XX line", hom_I: "HOM overlap I", hom_V: "HOM visibility V",
  rate_detected_Hz: "Detected rate (counts/s)", tau_ps: "Lifetime τ (ps)", tau_XX_ps: "XX lifetime τ (ps)",
  linewidth_ueV: "Emitter linewidth (µeV)", Gamma_ns: "Radiative rate Γ (ns⁻¹)", F_P_eff: "Effective Purcell factor F_P·L",
  B_fibre: "Brightness into the fibre", B_first_lens: "First-lens brightness",
};
const X_AXIS = { T_K: "Temperature (K)", "drive.f_rep_MHz": "Repetition rate (MHz)", "drive.fwhm_ps": "Pulse duration (ps)" };
export function yAxisLabel(y) {
  const base = OBS_AXIS[y?.observable] || String(y?.observable || "y").replace(/_/g, " ");
  const c = y?.convention && y.convention !== "intrinsic" ? `${y.convention} ` : "";
  return displayNotation(c + base);
}
export function xAxisLabel(x) {
  const path = x?.path || "x";
  if (X_AXIS[path]) return X_AXIS[path];
  const leaf = path.split(".").pop().replace(/_/g, " ");
  const u = x?.unit ? displayUnit(x.unit) : "";
  return u ? `${leaf} (${u})` : leaf;
}
const reqName = (y) => (y?.convention ? `${y.observable}:${y.convention}` : y?.observable || "y");

/**
 * series: card series entry {id, x:{path,unit}, y:{observable,convention?}, tag, where, digitized, method, points}
 * model:  reproduce model_series entry {x, y, tag, unit, failed, note} | {error} | null (not run yet)
 * inputs: {paths: [leaf paths the series' ledger entry depends on], fromPaper: Set of leaf paths taken from the paper}
 * Returns {el, chart, render()}.
 */
export function paperOverlay({ series, model = null, inputs = null, filename = null } = {}) {
  const ptag = normTag(series.tag);
  const hasModel = !!(model && !model.error && Array.isArray(model.x) && model.x.length);
  const mtag = hasModel ? normTag(model.tag) : null;
  const chart = chartCard({
    title: `${series.id.replace(/_/g, " ")}`,
    subtitle: "",
    height: 320,
    build: () => build(series, hasModel ? model : null),
    table: () => table(series, hasModel ? model : null),
    filename: () => filename || `fsim-paper-${series.id}`,
  });
  chart.setTag(widest([ptag, mtag]) || ptag);
  const sub = [
    `paper ${series.digitized ? "digitised" : "printed"} [${ptag || "?"}]`,
    hasModel ? `model [${mtag || "?"}]` : model?.error ? "model did not resolve" : "model: run Reproduce",
  ].join(" · ");
  chart.setTitle(series.id.replace(/_/g, " "), sub);

  const cap = h("p.pp-ov-cap", { dataset: { role: "params-source" } });
  const paths = inputs?.paths || [];
  const from = inputs?.fromPaper || new Set();
  const nPaper = paths.filter((p) => from.has(p)).length;
  const bits = [];
  if (hasModel) bits.push(model.note || "model at the central value of every other leaf");
  else if (model?.error) bits.push(`model series failed: ${model.error}`);
  else bits.push("model series appears after Reproduce");
  if (paths.length) bits.push(`parameters: ${paths.join(", ")} (${nPaper} of ${paths.length} taken from the paper, the rest are estimates or assumptions)`);
  if (series.where) bits.push(`paper points: ${series.where}`);
  if (series.digitized && series.method) bits.push(`digitised by ${series.method}`);
  if (hasModel && model.failed?.length) bits.push(`${model.failed.length} model point(s) failed and are not drawn`);
  cap.textContent = bits.join(". ") + ".";

  const el = h("div.pp-ov", { dataset: { series: series.id, model: hasModel ? "1" : "0" } }, chart.el, cap);
  return { el, chart, render: () => chart.render() };
}

function pts(series) {
  return (series.points || []).filter((p) => isNum(p[0]) && isNum(p[1])).slice().sort((a, b) => a[0] - b[0]);
}

function build(series, model) {
  const c = themeColors();
  const ink = c.ink1;
  const ptag = normTag(series.tag);
  const T = `[${ptag || "?"}]`;
  const P = pts(series);
  if (!P.length && !model) return null;
  const data = [];
  const open = !!series.digitized;
  if (P.length) {
    const errs = P.map((p) => (isNum(p[2]) ? p[2] : 0));
    data.push({
      x: P.map((p) => p[0]), y: P.map((p) => p[1]), type: "scatter", mode: "markers",
      name: `paper ${series.digitized ? "(digitised) " : ""}${T}`,
      marker: open
        ? { color: "rgba(0,0,0,0)", size: 10, symbol: "circle-open", line: { width: 2, color: ink } }
        : { color: ink, size: 9, symbol: "circle", line: { color: c.surface, width: 1.5 } },
      error_y: errs.some((e) => e > 0) ? { type: "data", array: errs, visible: true, color: ink, thickness: 1.25, width: 5 } : undefined,
      hovertemplate: `paper %{x}: %{y:.4g}${errs.some((e) => e > 0) ? " ± %{error_y.array:.3g}" : ""} ${T}<extra></extra>`,
    });
  }
  if (model) {
    const mtag = normTag(model.tag);
    const order = model.x.map((x, i) => i).sort((a, b) => model.x[a] - model.x[b]);
    data.push({
      x: order.map((i) => model.x[i]), y: order.map((i) => (isNum(model.y[i]) ? model.y[i] : null)),
      type: "scatter", mode: "lines+markers", connectgaps: false,
      name: `model [${mtag || "?"}]`,
      line: { color: c.series[0], width: tagWidth(mtag), dash: tagDash(mtag) },
      marker: { color: c.series[0], size: 6, symbol: "x-thin-open", line: { width: 2, color: c.series[0] } },
      hovertemplate: `model %{x}: %{y:.4g} [${mtag || "?"}]<extra></extra>`,
    });
  }
  const ys = data.flatMap((d) => d.y).filter((v) => isNum(v));
  const logY = ys.length > 1 && ys.every((v) => v > 0) && Math.max(...ys) / Math.min(...ys) > 100;
  const layout = baseLayout({
    hovermode: "closest", margin: { l: 62, r: 18, t: 30, b: 48 },
    xaxis: { ...baseLayout().xaxis, showspikes: false, title: { text: xAxisLabel(series.x), standoff: 6 } },
    yaxis: { ...baseLayout().yaxis, type: logY ? "log" : "linear", title: { text: yAxisLabel(series.y), standoff: 6 } },
  });
  return { data, layout };
}

function table(series, model) {
  const P = pts(series);
  const rows = P.map((p, i) => ({
    x: p[0], paper: p[1], err: isNum(p[2]) ? p[2] : null,
    model: model ? (model.x.indexOf(p[0]) >= 0 ? model.y[model.x.indexOf(p[0])] : null) : null,
  }));
  return {
    caption: `${series.id}: ${rows.length} paper points${model ? " and the model at the same x" : ""}`,
    columns: [
      { key: "x", label: xAxisLabel(series.x) },
      { key: "paper", label: `paper ${yAxisLabel(series.y)}` },
      { key: "err", label: "± err" },
      ...(model ? [{ key: "model", label: "model" }] : []),
    ],
    rows,
  };
}
