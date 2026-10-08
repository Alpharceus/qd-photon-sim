// Sweep charts over committed out/ CSVs (Results explorer + Story proofs).
// Rules (05-science-content section 4): invalid rows are hollow grey markers with a count,
// NaN never plots as 0 (counted, with reasons), bounds and scenarios are paired series or
// facets (never averaged), thresholds are references in --ref ink, gates are a wash box.
// No physics: every value plotted is a CSV cell or a baked story number.
import { baseLayout, colors as themeColors, alpha } from "./theme.js";
import { tagDash, tagWidth } from "../ui/tagChip.js";
import { g2Top } from "./g2T.js";
import { fmt, isNum } from "../ui/format.js";

// Display names for sweep columns (labels only).
const LABELS = {
  T_hs_K: "Heatsink T_hs (K)", T_hs: "Heatsink T_hs (K)", T_K: "T (K)",
  g2_pulsed: "Pulsed g²(0)", g2_op: "g²(0) at operating point", g2_cw0: "CW g²(0)", g2_cw0_raw: "CW g²(0), raw",
  collected_flux_pulsed_s: "Collected flux (s⁻¹)", collected_flux_delivered_s: "Collected flux, delivered (s⁻¹)",
  collected_flux_cw_s: "Collected flux, CW (s⁻¹)",
  delta_xx_meV: "Δ_XX (meV)", gamma300_meV: "Γ(300 K) (meV)", irf_ps: "IRF (ps)",
  emission_NA: "Collection NA", emission_R_back: "Back-facet R", emission_L_um: "Ridge length L (µm)",
  height_nm: "Dot height (nm)", radius_nm: "Dot radius (nm)", core_radius_nm: "Core radius (nm)",
  x_in: "Indium fraction x_in", current_uA: "Current (µA)", Q: "Cavity Q", rep_rate_hz: "Repetition rate (Hz)",
  set_EC_over_kT: "E_C / kT", lambda_nm: "λ (nm)", S_X: "Retention S_X", photons_per_cycle: "Photons per cycle",
  screening_fraction: "Screening fraction s", strain_bound: "Strain bound", family: "Family", regime: "Regime",
  orientation: "Orientation", card_id: "Card",
};

/** Axis label for a column; nanowire flux is never a bare "flux". */
export function colLabel(c, { nanowire = false } = {}) {
  if (nanowire && c === "collected_flux_pulsed_s") return "Collected flux, commanded (s⁻¹)";
  return LABELS[c] || c;
}
export const isG2 = (c) => /(^|_)g2(_|$)/.test(String(c)) && !/valid|reason|tag/.test(c);
export const isFlux = (c) => /flux/.test(String(c)) && !/margin|excluded|reason/.test(c);

function distinct(rows, key) {
  const s = new Map();
  for (const r of rows) { const v = r[key]; const k = v == null ? "n/a" : String(v); if (!s.has(k)) s.set(k, v); }
  return [...s.entries()].sort((a, b) => {
    const x = a[1], y = b[1];
    if (isNum(x) && isNum(y)) return x - y;
    return String(a[0]).localeCompare(String(b[0]));
  }).map(([k, v]) => ({ key: k, value: v }));
}

const shortV = (v) => (typeof v === "number" ? fmt(v) : v == null ? "n/a" : String(v));

/**
 * opts: {rows: [obj], x, y | y: [colA, colB] (paired series), color?, facet?, logX?, logY?,
 *        validKey?, reasonKey?, idKey?, greyRows?: [obj], greyLabel?, refs: {ceiling?, floor?,
 *        gate?: {g2: number, flux: number}, identity?}, maxFacets = 4, tag?, nanowire?, gl?}
 * Returns {data, layout, stats} or null.
 */
export function buildSweepScatter(opts) {
  const c = themeColors();
  const rows = opts.rows || [];
  if (!rows.length && !(opts.greyRows || []).length) return null;
  const ys = Array.isArray(opts.y) ? opts.y : [opts.y];
  const paired = ys.length > 1;
  const nanowire = !!opts.nanowire;
  const facets = opts.facet ? distinct(rows, opts.facet) : [{ key: "_all", value: null }];
  const shownFacets = facets.slice(0, opts.maxFacets || 4);
  const nF = shownFacets.length;
  const gl = opts.gl ?? (rows.length + (opts.greyRows?.length || 0)) > 1000;
  const type = gl ? "scattergl" : "scatter";
  const stats = { shown: 0, invalid: 0, invalidDrawn: 0, nan: 0, nonpos: 0, above1: 0, offcanvas: 0, reasons: new Map(), droppedFacets: facets.length - nF, facetValues: facets.map((f) => f.key) };
  const g2Axis = isG2(ys[0]) && !opts.logY;
  const plottedY = [];

  // colour categories (fixed order, never cycled past 8)
  let colorCats = null, colorNumeric = false;
  if (!paired && opts.color) {
    const d = distinct(rows, opts.color);
    if (d.length <= 8) colorCats = d.map((x) => x.key);
    else colorNumeric = d.every((x) => isNum(x.value));
  }
  const tagTxt = opts.tag ? ` [${opts.tag}]` : "";
  const xLab = colLabel(opts.x, { nanowire });
  const yLab = paired ? ys.map((y) => colLabel(y, { nanowire })).join(" vs ") : colLabel(ys[0], { nanowire });

  const data = [];
  const gap = nF > 1 ? 0.05 : 0;
  const w = (1 - gap * (nF - 1)) / nF;
  const layout = baseLayout({ hovermode: "closest", showlegend: true, margin: { l: 62, r: 16, t: nF > 1 ? 46 : 34, b: 50 } });
  const shapes = [];
  const annotations = [];
  const seenLegend = new Set();
  const legendOnce = (name) => { if (seenLegend.has(name)) return false; seenLegend.add(name); return true; };
  const ok = (v, log) => isNum(v) && (!log || v > 0);

  let xmin = Infinity, xmax = -Infinity, ymin = Infinity, ymax = -Infinity;
  const allX = rows.map((r) => r[opts.x]).filter((v) => isNum(v) && (!opts.logX || v > 0));
  const gateX1 = allX.length ? Math.max(...allX) * (opts.logX ? 1.6 : 1.05) : 1;
  const track = (x, y) => { xmin = Math.min(xmin, x); xmax = Math.max(xmax, x); ymin = Math.min(ymin, y); ymax = Math.max(ymax, y); };

  shownFacets.forEach((f, fi) => {
    const ax = fi === 0 ? "" : String(fi + 1);
    const xa = `x${ax}`, ya = `y${ax}`;
    const inF = (r) => (opts.facet ? (r[opts.facet] == null ? "n/a" : String(r[opts.facet])) === f.key : true);
    const fr = rows.filter(inF);
    const grey = (opts.greyRows || []).filter(inF);

    const pushPts = (list, ycol, style, name, showlegend, group) => {
      const X = [], Y = [], T = [], C = [], K = [];
      for (const r of list) {
        const xv = r[opts.x], yv = r[ycol];
        if (!ok(xv, opts.logX) || !ok(yv, opts.logY)) {
          if (style.kind === "valid") {
            if (isNum(xv) && isNum(yv)) stats.nonpos++;
            else {
              stats.nan++;
              const why = (opts.reasonKey && r[opts.reasonKey]) || "value is n/a in the CSV";
              stats.reasons.set(String(why).slice(0, 90), (stats.reasons.get(String(why).slice(0, 90)) || 0) + 1);
            }
          }
          continue;
        }
        X.push(xv); Y.push(yv); K.push(r.__v ?? null); track(xv, yv); plottedY.push(yv);
        if (g2Axis && yv > 1 && style.kind === "valid") stats.above1++;
        if (style.colorCol) C.push(r[style.colorCol]);
        const id = opts.idKey ? `${r[opts.idKey]} · ` : "";
        const col = opts.color && !paired ? ` · ${colLabel(opts.color)} ${shortV(r[opts.color])}` : "";
        const bad = style.kind === "invalid" && opts.reasonKey && r[opts.reasonKey] ? `<br>invalid: ${String(r[opts.reasonKey]).slice(0, 80)}` : "";
        T.push(`${id}${colLabel(opts.x, { nanowire })} ${fmt(xv)}<br>${colLabel(ycol, { nanowire })} ${fmt(yv)}${tagTxt}${col}${bad}`);
        if (style.kind === "valid") stats.shown++;
        else if (style.kind === "invalid") stats.invalidDrawn++;
      }
      if (!X.length) return;
      const marker = style.kind === "invalid"
        ? { symbol: "circle-open", size: 7, color: c.muted, line: { width: 1.25, color: c.muted } }
        : style.kind === "grey"
          ? { symbol: "circle", size: 6, color: alpha(c.muted, 0.28), line: { width: 0 } }
          : { symbol: style.symbol || "circle", size: 8, opacity: 0.85, color: style.colorCol ? C : style.color, line: { width: 1, color: c.surface }, ...(style.colorCol ? { colorscale: seqScale(), showscale: fi === nF - 1, colorbar: { title: { text: colLabel(opts.color), side: "right" }, thickness: 10, len: 0.8 } } : {}) };
      data.push({ x: X, y: Y, type, mode: "markers", name, legendgroup: group || name, showlegend, xaxis: xa, yaxis: ya,
        marker, text: T, customdata: K, hovertemplate: "%{text}<extra></extra>" });
    };

    if (grey.length) {
      const name = opts.greyLabel || "other model combos (never gate)";
      pushPts(grey, ys[0], { kind: "grey" }, name, legendOnce(name), "grey");
    }
    const valid = (r) => !opts.validKey || r[opts.validKey] !== false;
    const inv = fr.filter((r) => !valid(r));
    const good = fr.filter(valid);
    stats.invalid += inv.length;
    ys.forEach((ycol, si) => {
      if (paired) {
        const name = colLabel(ycol, { nanowire });
        pushPts(good, ycol, { kind: "valid", color: c.series[si], symbol: si ? "diamond" : "circle" }, name, legendOnce(name), ycol);
      } else if (colorCats) {
        colorCats.forEach((cat, ci) => {
          const sub = good.filter((r) => (r[opts.color] == null ? "n/a" : String(r[opts.color])) === cat);
          const name = `${colLabel(opts.color)} ${cat}`;
          pushPts(sub, ycol, { kind: "valid", color: c.series[ci] }, name, legendOnce(name), name);
        });
      } else if (colorNumeric) {
        pushPts(good, ycol, { kind: "valid", colorCol: opts.color }, colLabel(ycol, { nanowire }), false, "num");
      } else {
        const name = colLabel(ycol, { nanowire });
        pushPts(good, ycol, { kind: "valid", color: c.series[0] }, name, false, name);
      }
    });
    if (inv.length) pushPts(inv, ys[0], { kind: "invalid" }, "invalid rows (hollow)", legendOnce("invalid rows (hollow)"), "invalid");

    const dom = [fi * (w + gap), fi * (w + gap) + w];
    layout[`xaxis${ax}`] = {
      ...baseLayout().xaxis, showspikes: false, domain: dom, anchor: ya, type: opts.logX ? "log" : "linear",
      title: { text: xLab, standoff: 6, font: { size: 12 } }, automargin: true,
    };
    layout[`yaxis${ax}`] = {
      ...baseLayout().yaxis, anchor: xa, type: opts.logY ? "log" : "linear", automargin: true,
      title: fi === 0 ? { text: yLab, standoff: 6, font: { size: 12 } } : undefined,
      ...(g2Axis ? { dtick: 0.25 } : {}),
      ...(fi > 0 ? { matches: "y", showticklabels: true } : {}),
    };
    if (opts.facet) {
      annotations.push({ xref: "paper", yref: "paper", x: (dom[0] + dom[1]) / 2, y: 1.0, yanchor: "bottom", yshift: 4,
        text: `${colLabel(opts.facet)}: ${f.key}`, showarrow: false, font: { size: 12, color: c.ink2 } });
    }
    // references per panel
    const refs = opts.refs || {};
    if (isNum(refs.ceiling) && isG2(ys[0])) {
      shapes.push({ type: "line", xref: `${xa} domain`, x0: 0, x1: 1, yref: ya, y0: refs.ceiling, y1: refs.ceiling, line: { color: c.ref, width: 1 } });
      if (fi === nF - 1) annotations.push({ xref: `${xa} domain`, x: 1, xanchor: "right", yref: ya, y: refs.ceiling, yanchor: "bottom", text: `g²(0) = ${fmt(refs.ceiling)} ceiling${refs.ceilingTag ? ` [${refs.ceilingTag}]` : ""}`, showarrow: false, font: { size: 11, color: c.ink2 } });
    }
    if (isNum(refs.floor) && isFlux(opts.x)) {
      shapes.push({ type: "line", yref: `${ya} domain`, y0: 0, y1: 1, xref: xa, x0: refs.floor, x1: refs.floor, line: { color: c.ref, width: tagWidth(refs.floorTag || "A", 1.5), dash: tagDash(refs.floorTag || "A") } });
      if (fi === 0) annotations.push({ yref: `${ya} domain`, y: 1, yanchor: "top", xref: xa, x: opts.logX ? Math.log10(refs.floor) : refs.floor, xanchor: "left", xshift: 4, text: `flux floor ${fmt(refs.floor)} s⁻¹ [${refs.floorTag || "A"}]`, showarrow: false, font: { size: 11, color: c.ink2 } });
    }
    if (refs.gate && isNum(refs.gate.g2) && isNum(refs.gate.flux) && isFlux(opts.x) && isG2(ys[0])) {
      shapes.push({ type: "rect", xref: xa, yref: ya, x0: refs.gate.flux, x1: Math.max(gateX1, refs.gate.flux * 1.5), y0: 0, y1: refs.gate.g2, fillcolor: c.refWash, line: { color: c.ref, width: 1 }, layer: "below" });
      if (fi === 0) annotations.push({ xref: xa, yref: ya, x: opts.logX ? Math.log10(refs.gate.flux) : refs.gate.flux, y: refs.gate.g2 / 2, xanchor: "left", xshift: 8, text: "gate box: g²(0) below the ceiling and flux above the floor", showarrow: false, font: { size: 11, color: c.ink2 }, align: "left" });
    }
  });
  if (!stats.shown && !data.length) return { data: [], layout, stats, empty: true };
  // identity line (delivered vs commanded): y = x across the shared range
  if (opts.refs?.identity && Number.isFinite(xmin)) {
    const lo = Math.min(xmin, ymin), hi = Math.max(xmax, ymax);
    shownFacets.forEach((f, fi) => {
      const ax = fi === 0 ? "" : String(fi + 1);
      data.push({ x: [lo, hi], y: [lo, hi], type: "scatter", mode: "lines", xaxis: `x${ax}`, yaxis: `y${ax}`, name: "delivered = commanded (y = x)",
        line: { color: c.ref, width: 1 }, hoverinfo: "skip", showlegend: fi === 0, legendgroup: "identity" });
    });
  }
  // g2 axes run 0..1 unless the data rises above 1: then the axis is extended, never clipped
  if (g2Axis) {
    const top = g2Top([...plottedY, isNum(opts.refs?.ceiling) ? opts.refs.ceiling : 0]);
    shownFacets.forEach((f, fi) => { layout[`yaxis${fi === 0 ? "" : String(fi + 1)}`].range = [0, top]; });
    stats.yTop = top;
    stats.offcanvas = plottedY.filter((v) => v > top || v < 0).length;
  }
  layout.shapes = shapes;
  layout.annotations = annotations;
  layout.legend = { orientation: "h", x: 0, xanchor: "left", y: -0.16, yanchor: "top", font: { size: 12 } };
  layout.margin.b = 92;
  return { data, layout, stats };
}

function seqScale() {
  const c = getComputedStyle(document.documentElement);
  const s = [1, 2, 3, 4, 5, 6, 7].map((i) => c.getPropertyValue(`--seq-${i}`).trim()).filter(Boolean);
  if (!s.length) return "Blues";
  return s.map((col, i) => [i / (s.length - 1), col]);
}

/** Category bars against an optional reference line. bars: [{label, value, tag}] */
export function buildBars({ bars, yTitle, floor = null, floorTag = "A", floorLabel = "", logY = false, color = null }) {
  const c = themeColors();
  const vals = bars.filter((b) => isNum(b.value) && (!logY || b.value > 0));
  if (!vals.length) return null;
  const data = [{
    x: vals.map((b) => b.label), y: vals.map((b) => b.value), type: "bar", name: yTitle,
    marker: { color: color || c.series[0], line: { width: 0 } }, width: 0.46,
    text: vals.map((b) => fmt(b.value)), textposition: "outside", cliponaxis: false, textfont: { color: c.ink1, size: 13 },
    customdata: vals.map((b) => b.tag || ""), hovertemplate: "%{x}: %{y:.4g} [%{customdata}]<extra></extra>",
  }];
  const shapes = [], annotations = [];
  if (isNum(floor)) {
    shapes.push({ type: "line", xref: "paper", x0: 0, x1: 1, yref: "y", y0: floor, y1: floor, line: { color: c.ref, width: tagWidth(floorTag, 1.5), dash: tagDash(floorTag) } });
    annotations.push({ xref: "paper", x: 1, xanchor: "right", yref: "y", y: logY ? Math.log10(floor) : floor, yanchor: "bottom", text: floorLabel || `floor ${fmt(floor)} [${floorTag}]`, showarrow: false, font: { size: 12, color: c.ink2 } });
  }
  const layout = baseLayout({
    hovermode: "closest", showlegend: false, bargap: 0.4, margin: { l: 64, r: 16, t: 26, b: 46 },
    xaxis: { ...baseLayout().xaxis, showspikes: false, type: "category" },
    yaxis: { ...baseLayout().yaxis, title: { text: yTitle, standoff: 6 }, type: logY ? "log" : "linear", rangemode: logY ? undefined : "tozero" },
    shapes, annotations,
  });
  return { data, layout };
}

/** V-a calibration: data with error bars (upper bounds as down-triangles), model curve, residual strip. */
export function buildCalibration(va) {
  const c = themeColors();
  const pts = va.points || [];
  const vals = pts.filter((p) => p.bound !== "upper");
  const ups = pts.filter((p) => p.bound === "upper");
  const tag = String(va.tag_chain || "").replace(/[[\]]/g, "") || "A";
  const data = [
    { x: va.curve.T_K, y: va.curve.g2_model, type: "scatter", mode: "lines", name: `model g²(T) (fit) [${tag}]`,
      line: { color: c.series[0], width: tagWidth(tag), dash: tagDash(tag) }, xaxis: "x", yaxis: "y", hovertemplate: `model %{y:.4f} [${tag}]<extra></extra>` },
    { x: vals.map((p) => p.T_K), y: vals.map((p) => p.g2), type: "scatter", mode: "markers", name: "measured, digitized [V]",
      marker: { color: c.ink1, size: 9, line: { color: c.surface, width: 2 } },
      error_y: { type: "data", array: vals.map((p) => p.err), color: c.ink1, thickness: 1.25, width: 5 }, xaxis: "x", yaxis: "y",
      hovertemplate: "measured %{y:.3f} ± %{error_y.array:.3f} [V]<extra></extra>" },
    { x: ups.map((p) => p.T_K), y: ups.map((p) => p.g2), type: "scatter", mode: "markers", name: "upper bound (≤) [V]",
      marker: { color: c.ink1, size: 12, symbol: "triangle-down-open", line: { width: 2, color: c.ink1 } }, xaxis: "x", yaxis: "y",
      hovertemplate: "upper bound ≤ %{y:.3f} [V]<extra></extra>" },
    { x: pts.map((p) => p.T_K), y: pts.map((p) => p.residual), type: "bar", name: "residual (model − data)",
      marker: { color: alpha(c.series[0], 0.55) }, width: 4, xaxis: "x2", yaxis: "y2", showlegend: false,
      hovertemplate: `residual %{y:+.4f} [${tag}]<extra></extra>` },
  ];
  const layout = baseLayout({
    hovermode: "closest", margin: { l: 62, r: 18, t: 30, b: 46 },
    legend: { orientation: "h", x: 0, y: 1.0, yanchor: "bottom", font: { size: 12 } },
    xaxis: { ...baseLayout().xaxis, showspikes: false, anchor: "y", matches: "x2", showticklabels: false, domain: [0, 1] },
    yaxis: { ...baseLayout().yaxis, domain: [0.34, 1], range: [0, 0.6], title: { text: "g²(0)", standoff: 6 } },
    xaxis2: { ...baseLayout().xaxis, showspikes: false, anchor: "y2", domain: [0, 1], title: { text: "Temperature T (K)", standoff: 6 } },
    yaxis2: { ...baseLayout().yaxis, domain: [0, 0.24], title: { text: "residual", standoff: 6 }, zeroline: true, zerolinecolor: c.ref },
    shapes: isNum(va.Tc_K) ? [{ type: "line", xref: "x", x0: va.Tc_K, x1: va.Tc_K, yref: "y domain", y0: 0, y1: 1, line: { color: c.ref, width: 1 } }] : [],
    annotations: isNum(va.Tc_K) ? [{ xref: "x", x: va.Tc_K, yref: "y domain", y: 1, yanchor: "top", xanchor: "right", xshift: -4, text: `T_c = ${fmt(va.Tc_K)} K (fit) [${tag}]`, showarrow: false, font: { size: 12, color: c.ink2 } }] : [],
  });
  return { data, layout };
}

/** Model envelope intervals at the measured temperatures, with the measured points (no mid curve). */
export function buildIntervals({ T, lo, hi, measured, tagModel = "A", tagData = "V", ceiling = null }) {
  const c = themeColors();
  const data = [
    { x: T, y: hi.map((h, i) => h - lo[i]), base: lo, type: "bar", width: 14, name: `model envelope lo – hi [${tagModel}]`,
      marker: { color: alpha(c.series[0], 0.18), line: { color: alpha(c.series[0], 0.75), width: 1.25 } },
      customdata: T.map((_, i) => [lo[i], hi[i]]), hovertemplate: `envelope %{customdata[0]:.3f} – %{customdata[1]:.3f} [${tagModel}]<extra></extra>` },
    { x: T, y: measured, type: "scatter", mode: "markers", name: `measured [${tagData}]`,
      marker: { color: c.ink1, size: 11, line: { color: c.surface, width: 2 } }, hovertemplate: `measured %{y:.3f} [${tagData}]<extra></extra>` },
  ];
  const layout = baseLayout({
    hovermode: "closest", margin: { l: 62, r: 18, t: 30, b: 46 },
    xaxis: { ...baseLayout().xaxis, showspikes: false, title: { text: "Temperature T (K)", standoff: 6 }, range: [40, 330] },
    yaxis: { ...baseLayout().yaxis, title: { text: "g²(0)", standoff: 6 }, range: [0, g2Top([...hi, ...measured])], dtick: 0.25 },
    shapes: ceiling ? [{ type: "line", xref: "paper", x0: 0, x1: 1, yref: "y", y0: ceiling.value, y1: ceiling.value, line: { color: c.ref, width: 1 } }] : [],
    annotations: ceiling ? [{ xref: "paper", x: 0, xanchor: "left", yref: "y", y: ceiling.value, yanchor: "bottom", text: `g²(0) = ${fmt(ceiling.value)} ceiling [${ceiling.tag || "A"}]`, showarrow: false, font: { size: 11, color: c.ink2 } }] : [],
  });
  return { data, layout };
}

/** Dot plot (markers only; two temperatures are never joined into a curve). series: [{name, points:[{x, y, tag}]}] */
export function buildDots(panels) {
  const c = themeColors();
  const data = [];
  const n = panels.length;
  const gap = 0.12;
  const w = (1 - gap * (n - 1)) / n;
  const layout = baseLayout({ hovermode: "closest", showlegend: false, margin: { l: 62, r: 16, t: 34, b: 46 } });
  const annotations = [];
  panels.forEach((p, i) => {
    const ax = i === 0 ? "" : String(i + 1);
    data.push({ x: p.points.map((q) => q.x), y: p.points.map((q) => q.y), type: "scatter", mode: "markers+text",
      xaxis: `x${ax}`, yaxis: `y${ax}`, marker: { size: 13, color: c.series[i], line: { color: c.surface, width: 2 } },
      text: p.points.map((q) => `${fmt(q.y)} [${q.tag}]`), textposition: "middle right", textfont: { size: 13, color: c.ink1 },
      cliponaxis: false, hovertemplate: `${p.name} at %{x} K: %{y:.4g}<extra></extra>` });
    const dom = [i * (w + gap), i * (w + gap) + w];
    layout[`xaxis${ax}`] = { ...baseLayout().xaxis, showspikes: false, domain: dom, anchor: `y${ax}`, title: { text: "T_hs (K)", standoff: 6 }, range: [200, 330] };
    layout[`yaxis${ax}`] = { ...baseLayout().yaxis, anchor: `x${ax}`, rangemode: "tozero", title: { text: p.yTitle || "", standoff: 6 } };
    annotations.push({ xref: "paper", yref: "paper", x: (dom[0] + dom[1]) / 2, y: 1, yanchor: "bottom", yshift: 6, text: p.name, showarrow: false, font: { size: 12, color: c.ink2 } });
  });
  layout.annotations = annotations;
  return { data, layout };
}
