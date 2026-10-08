// Plotly theming: the fsim template comes from /api/theme/plotly?mode= (generated from
// fsim_theme/tokens.json). Every live chart registers here; a theme switch re-templates all.
import { plotlyTemplate } from "../api.js";
import { cssVar } from "../ui/dom.js";

const live = new Set();
let mode = document.documentElement.dataset.theme || "dark";
const templates = {};

export async function loadTemplate(m = mode) {
  if (!templates[m]) {
    try { templates[m] = await plotlyTemplate(m); } catch { templates[m] = null; }
  }
  return templates[m];
}

export function template() { return templates[mode] || undefined; }
export function themeMode() { return mode; }

export async function setChartTheme(m) {
  mode = m;
  await loadTemplate(m);
  for (const c of live) {
    try { c.render(); } catch (e) { console.warn("chart re-theme failed", e); }
  }
}

/** Re-render every live chart (e.g. Present mode changes the stroke scale). */
export function rerenderCharts() {
  for (const c of live) { try { c.render(); } catch (e) { console.warn("chart re-render failed", e); } }
}

/** Faint minor graticule for instrument screens (the g2(T) screen). */
export function graticule() {
  const c = colors();
  return { showgrid: true, gridcolor: alpha(c.grid && c.grid.startsWith("#") ? c.grid : "#808080", 0.45), gridwidth: 0.5, nticks: 4 };
}

export function register(chart) { live.add(chart); return () => live.delete(chart); }

/** Theme colours read from the generated tokens at render time. */
export function colors() {
  return {
    ink1: cssVar("--ink-1"), ink2: cssVar("--ink-2"), muted: cssVar("--muted"),
    ref: cssVar("--ref"), refWash: cssVar("--ref-wash"), grid: cssVar("--grid"), axis: cssVar("--axis"),
    surface: cssVar("--surface"), beam: cssVar("--beam"),
    series: [1, 2, 3, 4, 5, 6, 7, 8].map((i) => cssVar(`--series-${i}`)),
  };
}

/** "#rrggbb" + alpha -> rgba() */
export function alpha(hex, a) {
  const m = /^#?([0-9a-f]{6})$/i.exec(String(hex).trim());
  if (!m) return hex;
  const n = parseInt(m[1], 16);
  return `rgba(${(n >> 16) & 255},${(n >> 8) & 255},${n & 255},${a})`;
}

export const PLOT_CONFIG = { displayModeBar: false, responsive: true, displaylogo: false };

export function baseLayout(extra = {}) {
  const c = colors();
  const tpl = template();
  return {
    template: tpl,
    autosize: true,
    margin: { l: 54, r: 18, t: 30, b: 42 },
    hovermode: "x unified",
    hoverdistance: 40,
    spikedistance: -1,
    showlegend: true,
    legend: { orientation: "h", x: 0, xanchor: "left", y: 1.0, yanchor: "bottom", font: { size: 12 } },
    xaxis: {
      showspikes: true, spikemode: "across", spikesnap: "cursor", spikethickness: 1,
      spikecolor: c.ref, spikedash: "solid", zeroline: false, ticks: "outside",
    },
    yaxis: { zeroline: false, ticks: "outside" },
    font: { family: "Barlow, Segoe UI, system-ui, sans-serif" },
    ...extra,
  };
}
