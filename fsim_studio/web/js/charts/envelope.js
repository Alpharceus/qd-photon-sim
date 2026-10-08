// The band rule, defined once (02-charts.md section 3):
// band fill in the series hue at 12% (8% when bands overlap), lo/hi edges at 1.25px and 60%
// opacity, mid curve labelled "all ranges at midpoint (not a prediction)". Never a mean.
import { alpha } from "./theme.js";
import { tagDash, tagWidth } from "../ui/tagChip.js";
import { fmt, isNum } from "../ui/format.js";
import { gatesNow } from "../api.js";

export const MID_LABEL = "all ranges at midpoint (not a prediction)";

const clean = (arr) => (arr || []).map((v) => (isNum(v) ? v : null));

/**
 * x, band {lo, hi, mid}; returns Plotly traces [hi edge, lo edge (fill to hi), mid].
 * name: series name for the legend; tag: provenance tag (line dash grammar).
 */
export function bandTraces({ x, band, color, name, tag, overlap = false, showMid = true, legendgroup, fmtHover = (v) => fmt(v) , valueLabel = "value", showLegend = true }) {
  const lo = clean(band.lo);
  const hi = clean(band.hi);
  const mid = clean(band.mid);
  const fillA = overlap ? 0.08 : 0.12;
  const group = legendgroup || name;
  const hover = `${valueLabel === "value" ? name : valueLabel} band %{customdata[0]:.4g} – %{customdata[1]:.4g} ${tag ? `[${tag}]` : ""}<extra></extra>`;
  // Without a mid curve the hover is anchored on the hi edge itself: no averaged series exists.
  const hiTrace = {
    x, y: hi, type: "scatter", mode: "lines", name: `${name} upper edge`,
    line: { color: alpha(color, 0.6), width: 1.25 }, showlegend: false, legendgroup: group,
    ...(showMid ? { hoverinfo: "skip" } : { customdata: x.map((_, i) => [lo[i], hi[i]]), hovertemplate: hover }),
  };
  const loTrace = {
    x, y: lo, type: "scatter", mode: "lines", name: `${name} band`,
    line: { color: alpha(color, 0.6), width: 1.25 }, fill: "tonexty", fillcolor: alpha(color, fillA),
    hoverinfo: "skip", showlegend: showLegend, legendgroup: group,
  };
  const traces = [hiTrace, loTrace];
  if (showMid) {
    const t = tag ? `[${tag}]` : "";
    traces.push({
      x, y: mid, type: "scatter", mode: "lines", name: `${name}: ${MID_LABEL}`,
      line: { color, width: tagWidth(tag || "A"), dash: tag ? tagDash(tag) : "dash" },
      customdata: x.map((_, i) => [lo[i], hi[i]]),
      hovertemplate: `${valueLabel} band %{customdata[0]:.4g} – %{customdata[1]:.4g} · mid %{y:.4g} ${t}<extra></extra>`,
      legendgroup: group, showlegend: showLegend,
    });
  }
  return traces;
}

/** T_c rendering: a 1px reference rule (point run) or a reference-wash span (envelope). */
export function tcMarks({ Tc, TcBand, colors, xRange, singleT = false }) {
  const shapes = [];
  const annotations = [];
  if (Array.isArray(TcBand)) {
    const [lo, hi] = TcBand;
    const x0 = isNum(lo) ? lo : xRange[0];
    const x1 = isNum(hi) ? hi : xRange[1];
    const text = isNum(lo) && isNum(hi) ? `T_c ${fmt(lo)} – ${fmt(hi)} K`
      : isNum(lo) || isNum(hi) ? `T_c ${isNum(lo) ? fmt(lo) : "not crossed"} – ${isNum(hi) ? fmt(hi) : "not crossed in range"} K`
        : "T_c not crossed in range";
    if (isNum(lo) || isNum(hi)) {
      shapes.push({ type: "rect", xref: "x", yref: "paper", x0, x1, y0: 0, y1: 1, fillcolor: colors.refWash, line: { width: 0 }, layer: "below" });
    }
    annotations.push({ xref: isNum(lo) || isNum(hi) ? "x" : "paper", x: isNum(lo) || isNum(hi) ? (x0 + x1) / 2 : 0.02, xanchor: isNum(lo) || isNum(hi) ? "center" : "left", yref: "paper", y: 1, yanchor: "top", text, showarrow: false, font: { size: 12, color: colors.ink2 }, bgcolor: colors.surface, borderpad: 2 });
  } else if (isNum(Tc)) {
    shapes.push({ type: "line", xref: "x", yref: "paper", x0: Tc, x1: Tc, y0: 0, y1: 1, line: { color: colors.ref, width: 1 } });
    annotations.push({ xref: "x", x: Tc, yref: "paper", y: 1, yanchor: "top", xanchor: "left", text: `T_c = ${fmt(Tc)} K`, showarrow: false, font: { size: 12, color: colors.ink2 }, bgcolor: colors.surface, borderpad: 2, xshift: 4 });
  } else if (Tc === null) {
    annotations.push({ xref: "paper", x: 0.01, xanchor: "left", yref: "paper", y: 1, yanchor: "top", text: singleT ? "T_c needs the full T_hs grid (single-T run)" : `T_c: g²(0) does not cross ${isNum(gatesNow()?.g2_ceiling?.value) ? fmt(gatesNow().g2_ceiling.value) : "the ceiling"} in range`, showarrow: false, font: { size: 12, color: colors.muted } });
  }
  return { shapes, annotations };
}

/** Rows for a Table toggle of a banded series. */
export function bandRows(x, band) {
  return x.map((t, i) => ({ T: t, lo: band.lo?.[i], mid: band.mid?.[i], hi: band.hi?.[i] }));
}
