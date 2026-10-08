// Hero chart: g2(0) against heatsink temperature. y from 0 to 1 (expanded when the data rises
// above 1, never clipped) with the g2 ceiling from /api/gates labelled, T_c as a rule (point) or
// span (envelope), operating-point marker, bands for envelopes, a faint graticule.
// Also the secondary single-quantity plots (eps, rho^2, dT_J, Gamma): separate plots, never dual axis.
// A run whose result carries no provenance tag is not drawn (CONTRACT rule 2): it is listed in
// spec.withheld and the chart card shows the reason.
import { baseLayout, colors as themeColors, alpha, graticule } from "./theme.js";
import { bandTraces, tcMarks, MID_LABEL } from "./envelope.js";
import { tagDash, tagWidth, normTag } from "../ui/tagChip.js";
import { fmt, isNum } from "../ui/format.js";
import { gatesNow } from "../api.js";

/** Upper y limit for a g2 axis: 1.02 unless the plotted data rises above 1. */
export function g2Top(values) {
  let m = 0;
  for (const v of values) if (isNum(v) && v > m) m = v;
  return m > 1 ? Math.ceil(m * 1.04 * 20) / 20 : 1.02;
}

/** The g2 ceiling reference from /api/gates: {value, tag, source} or null when not loaded. */
export function ceilingRef() {
  const g = gatesNow()?.g2_ceiling;
  return g && isNum(g.value) ? g : null;
}

/**
 * runs: [{result, name, color?, slot?, opT}] ; options: {overlap}
 * Returns {data, layout} or null when nothing is plottable.
 */
export function buildG2T(runs, { compare = false } = {}) {
  const c = themeColors();
  const data = [];
  const shapes = [];
  const annotations = [];
  let xmin = Infinity, xmax = -Infinity;
  const all = runs.filter((r) => r && r.result && r.result.curves && Array.isArray(r.result.curves.T_hs));
  const withheld = all.filter((r) => !normTag(r.result.tag_chain)).map((r) => r.name || r.result.run_id || "run");
  const live = all.filter((r) => normTag(r.result.tag_chain));
  const ys = [];
  if (!live.length) return withheld.length ? { data: [], layout: {}, withheld } : null;
  const overlap = live.filter((r) => r.result.bands).length > 1;

  live.forEach((r, i) => {
    const res = r.result;
    const x = res.curves.T_hs;
    x.forEach((v) => { if (isNum(v)) { xmin = Math.min(xmin, v); xmax = Math.max(xmax, v); } });
    const color = r.color || c.series[0];
    const tag = normTag(res.tag_chain);
    const tagTxt = tag ? `[${tag}]` : "";
    const name = r.name || "g²(0)";
    if (res.bands && res.bands.g2) {
      ys.push(...(res.bands.g2.hi || []), ...(res.bands.g2.lo || []));
      data.push(...bandTraces({ x, band: res.bands.g2, color, name, tag, overlap, valueLabel: "g²(0)" }));
    } else if (x.length > 1) {
      ys.push(...(res.curves.g2 || []));
      data.push({
        x, y: res.curves.g2, type: "scatter", mode: "lines", name,
        line: { color, width: tagWidth(tag), dash: tagDash(tag) },
        hovertemplate: `${compare ? `${name} ` : ""}g²(0) %{y:.4f} ${tagTxt}<extra></extra>`,
      });
    }
    // operating point
    const opRes = r.opResult || res;
    const s = opRes.scalars || {};
    const opT = isNum(r.opT) ? r.opT : (x.length === 1 ? x[0] : null);
    const g2op = isNum(s.g2_op) ? s.g2_op : null;
    const g2b = opRes.scalar_bands?.g2_op;
    ys.push(g2op, ...(Array.isArray(g2b) ? g2b : []));
    if (isNum(opT) && g2op === null && Array.isArray(g2b) && isNum(g2b[0]) && isNum(g2b[1])) {
      // no mid-design value: show the interval itself, never its average
      data.push({
        x: [opT, opT], y: [g2b[0], g2b[1]], type: "scatter", mode: "lines+markers", name: `${name} operating-point interval`,
        line: { color, width: 2 }, marker: { size: 7, color: c.surface, line: { color, width: 2 } }, showlegend: false,
        hovertemplate: `operating point g²(0) ${fmt(g2b[0])} – ${fmt(g2b[1])} ${tagTxt}<extra></extra>`,
      });
    } else if (isNum(opT) && g2op !== null) {
      const y = g2op;
      const err = Array.isArray(g2b) && isNum(g2b[0]) && isNum(g2b[1]) && y !== null
        ? { type: "data", symmetric: false, array: [g2b[1] - y], arrayminus: [y - g2b[0]], color, thickness: 1.25, width: 6 } : undefined;
      data.push({
        x: [opT], y: [y], type: "scatter", mode: "markers", name: `${compare ? `${name} ` : ""}operating point`,
        marker: { size: 11, color: c.surface, line: { color, width: 2.5 }, symbol: "circle" },
        error_y: err, showlegend: !compare,
        hovertemplate: err
          ? `operating point T_hs %{x:.1f} K · g²(0) ${fmt(g2b[0])} – ${fmt(g2b[1])} ${tagTxt}<extra></extra>`
          : `operating point T_hs %{x:.1f} K · g²(0) %{y:.4f} ${tagTxt}<extra></extra>`,
      });
      if (x.length === 1 && !compare) {
        annotations.push({ x: opT, y, xref: "x", yref: "y", text: "single-T point", showarrow: false, yshift: 18, font: { size: 12, color: c.ink2 } });
      }
    }
    if (!compare) {
      const tc = tcMarks({ Tc: res.scalars?.T_c, TcBand: res.scalar_bands?.T_c, colors: c, xRange: [xmin, xmax], singleT: x.length === 1 });
      shapes.push(...tc.shapes);
      annotations.push(...tc.annotations);
    }
  });

  if (!Number.isFinite(xmin)) return null;
  if (xmin === xmax) { xmin -= 40; xmax += 40; }
  // g2 ceiling from /api/gates: a reference, not a verdict
  const ceil = ceilingRef();
  if (ceil) {
    shapes.push({ type: "line", xref: "paper", x0: 0, x1: 1, yref: "y", y0: ceil.value, y1: ceil.value, line: { color: c.ref, width: 1 } });
    annotations.push({ xref: "paper", x: 1, xanchor: "right", yref: "y", y: ceil.value, yanchor: "bottom", text: `g²(0) = ${fmt(ceil.value)} ceiling${ceil.tag ? ` [${String(ceil.tag).replace(/[[\]]/g, "")}]` : ""}`, showarrow: false, font: { size: 12, color: c.ink2 } });
  }

  const top = g2Top(ys);
  const layout = baseLayout({
    xaxis: { ...baseLayout().xaxis, title: { text: "Heatsink temperature T_hs (K)", standoff: 8 }, range: [xmin, xmax], minor: graticule() },
    yaxis: { ...baseLayout().yaxis, title: { text: "g²(0)", standoff: 6 }, range: [0, top], dtick: 0.25, fixedrange: true, minor: graticule() },
    shapes, annotations,
  });
  return { data, layout, withheld, yTop: top };
}

export const SECONDARY = {
  eps: { label: "XX leakage ε", unit: "", curve: "eps", band: "eps" },
  rho2: { label: "Loading purity ρ²", unit: "", curve: "rho2", band: "rho2", yrange: [0, 1.02] },
  dTJ: { label: "Self-heating ΔT_J = T_j − T_hs", unit: "K", curve: "Tj", band: "Tj", minusX: true },
  gamma: { label: "Linewidth Γ(T_j)", unit: "meV", curve: "gamma", band: "gamma" },
};

/** Secondary single-quantity plot for one or more runs. */
export function buildSecondary(key, runs, { compare = false } = {}) {
  const spec = SECONDARY[key];
  const c = themeColors();
  const all = runs.filter((r) => r?.result?.curves?.[spec.curve] && r.result.curves.T_hs?.length > 1);
  const withheld = all.filter((r) => !normTag(r.result.tag_chain)).map((r) => r.name || r.result.run_id || "run");
  const live = all.filter((r) => normTag(r.result.tag_chain));
  if (!live.length) return withheld.length ? { data: [], layout: {}, withheld } : null;
  const data = [];
  live.forEach((r, i) => {
    const res = r.result;
    const x = res.curves.T_hs;
    const sub = (arr) => (spec.minusX ? arr.map((v, j) => (isNum(v) && isNum(x[j]) ? v - x[j] : null)) : arr);
    const color = r.color || c.series[i % 8];
    const tag = normTag(res.tag_chain);
    const name = compare ? r.name : spec.label;
    if (res.bands?.[spec.band]) {
      const b = res.bands[spec.band];
      data.push(...bandTraces({ x, band: { lo: sub(b.lo), hi: sub(b.hi), mid: b.mid ? sub(b.mid) : null }, color, name, tag, overlap: live.length > 1, showMid: !!b.mid, valueLabel: spec.label, showLegend: compare }));
    } else {
      data.push({
        x, y: sub(res.curves[spec.curve]), type: "scatter", mode: "lines", name,
        line: { color, width: tagWidth(tag), dash: tagDash(tag) }, showlegend: compare,
        hovertemplate: `${compare ? `${name} ` : ""}%{y:.4g}${spec.unit ? ` ${spec.unit}` : ""} ${tag ? `[${tag}]` : ""}<extra></extra>`,
      });
    }
  });
  const top = spec.yrange ? g2Top(live.flatMap((r) => r.result.bands?.[spec.band]?.hi || r.result.curves[spec.curve] || [])) : null;
  const layout = baseLayout({
    showlegend: compare,
    margin: { l: 54, r: 18, t: compare ? 30 : 12, b: 40 },
    xaxis: { ...baseLayout().xaxis, title: { text: "T_hs (K)", standoff: 6 } },
    yaxis: { ...baseLayout().yaxis, title: { text: spec.unit ? `${spec.unit}` : "", standoff: 6 }, range: spec.yrange ? [spec.yrange[0], Math.max(spec.yrange[1], top)] : undefined, rangemode: spec.yrange ? undefined : "tozero" },
  });
  return { data, layout, withheld };
}

export function g2TableRows(result) {
  const x = result?.curves?.T_hs || [];
  if (result?.bands?.g2) {
    const b = result.bands.g2;
    return { caption: `g²(0) envelope: ${MID_LABEL}`, columns: [{ key: "T", label: "T_hs (K)" }, { key: "lo", label: "g² lo" }, { key: "mid", label: "g² mid" }, { key: "hi", label: "g² hi" }], rows: x.map((t, i) => ({ T: t, lo: b.lo[i], mid: b.mid?.[i], hi: b.hi[i] })) };
  }
  return { columns: [{ key: "T", label: "T_hs (K)" }, { key: "g2", label: "g²(0)" }, { key: "Tj", label: "T_j (K)" }], rows: x.map((t, i) => ({ T: t, g2: result.curves.g2[i], Tj: result.curves.Tj?.[i] })) };
}

export { alpha };

// ---------------------------------------------------------------- animation playhead (studio-p2a)
/**
 * Pixel position of the playhead marker (T_hs, g²(0)) inside a drawn g²(T) plot div, from
 * Plotly's own axis mapping (no trace is added: the marker is a DOM overlay the workspace moves
 * per frame, so playback never calls Plotly). Returns {x, y, inX, inY} or null before the first
 * draw. Values outside the axes are clamped to the plot edge and flagged.
 */
export function playheadXY(gd, T, g2) {
  const fl = gd?._fullLayout;
  const xa = fl?.xaxis, ya = fl?.yaxis;
  if (!xa || !ya || typeof xa.d2p !== "function" || !isNum(T)) return null;
  const [x0, x1] = xa.range, [y0, y1] = ya.range;
  const inX = T >= Math.min(x0, x1) && T <= Math.max(x0, x1);
  const inY = isNum(g2) && g2 >= Math.min(y0, y1) && g2 <= Math.max(y0, y1);
  const cx = Math.max(Math.min(x0, x1), Math.min(Math.max(x0, x1), T));
  const cy = isNum(g2) ? Math.max(Math.min(y0, y1), Math.min(Math.max(y0, y1), g2)) : Math.min(y0, y1);
  return { x: xa._offset + xa.d2p(cx), y: ya._offset + ya.d2p(cy), inX, inY, hasY: isNum(g2) };
}
