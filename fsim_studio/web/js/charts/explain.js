// Explain charts: the spectral window (port of app.py's Spectral explainer, redesigned per
// 02-charts.md section 2). Every array comes from /api/explain/spectral (fsim_core.spectral);
// this file only maps them to marks.
import { baseLayout, colors as themeColors, alpha } from "./theme.js";
import { fmt, isNum } from "../ui/format.js";

/** X slot 1 and XX slot 2 at 2px; filter window as a neutral wash; XX leakage at 25% (the reported quantity). */
export function buildSpectral(r) {
  if (!r?.curves) return null;
  const c = themeColors();
  const tag = r.tags?.chain || "A";
  const data = [
    { x: r.leakage.x, y: r.leakage.y, type: "scatter", mode: "none", fill: "tozeroy", fillcolor: alpha(c.series[1], 0.25),
      name: "t_XX leakage through the window", hoverinfo: "skip" },
    { x: r.curves.x, y: r.curves.X, type: "scatter", mode: "lines", name: "X", line: { color: c.series[0], width: 2 },
      hovertemplate: `X %{y:.3f} [${tag}]<extra></extra>` },
    { x: r.curves.x, y: r.curves.XX, type: "scatter", mode: "lines", name: "XX", line: { color: c.series[1], width: 2 },
      hovertemplate: `XX %{y:.3f} [${tag}]<extra></extra>` },
  ];
  const peakXX = Math.max(...r.curves.XX.filter(isNum));
  const annotations = [
    { x: r.window.center, y: 1, yref: "paper", yanchor: "bottom", showarrow: false, text: `filter w = ${fmt(r.w)} meV${r.published_window ? " (published window)" : ""}`, font: { size: 12, color: c.ink2 } },
    { x: 0, y: 1, text: "X", showarrow: false, yshift: 12, font: { size: 13, color: c.ink1 } },
    { x: -r.delta_xx, y: peakXX, text: "XX", showarrow: false, yshift: 12, font: { size: 13, color: c.ink1 } },
  ];
  if (r.eps_label) {
    annotations.push({ x: r.eps_label.x, y: r.eps_label.y, text: `ε = ${fmt(r.eps)} [${tag}]`, showarrow: true, arrowhead: 0, arrowcolor: c.ref, ax: 48, ay: -34, font: { size: 13, color: c.ink1 }, bgcolor: c.surface, borderpad: 2 });
  }
  const layout = baseLayout({
    hovermode: "x unified", margin: { l: 58, r: 18, t: 34, b: 46 },
    legend: { orientation: "h", x: 0, y: -0.18, yanchor: "top", font: { size: 12 } },
    xaxis: { ...baseLayout().xaxis, title: { text: "energy − E_X (meV)", standoff: 6 } },
    yaxis: { ...baseLayout().yaxis, title: { text: "peak-normalized intensity", standoff: 6 }, range: [0, 1.12] },
    shapes: [{ type: "rect", xref: "x", yref: "paper", x0: r.window.lo, x1: r.window.hi, y0: 0, y1: 1, fillcolor: c.refWash, line: { width: 1, color: c.ref }, layer: "below" }],
    annotations,
  });
  return { data, layout };
}

export function spectralTable(r) {
  if (!r?.curves) return null;
  const step = Math.ceil(r.curves.x.length / 240);
  const rows = [];
  for (let i = 0; i < r.curves.x.length; i += step) rows.push({ x: r.curves.x[i], X: r.curves.X[i], XX: r.curves.XX[i] });
  return { caption: `every ${step}th grid point of ${r.curves.x.length}`, columns: [{ key: "x", label: "E − E_X (meV)" }, { key: "X", label: "X" }, { key: "XX", label: "XX" }], rows };
}
