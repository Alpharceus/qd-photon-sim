// Chart frame shared by every chart: engraved title, tag-chain chip, Table toggle, export
// PNG/SVG (Plotly.toImage at 2x). Owns the Plotly lifecycle, resize and re-theme.
import { h, clear, reducedMotion } from "../ui/dom.js";
import { icon } from "../ui/icons.js";
import { tagChip } from "../ui/tagChip.js";
import { dataTable, toast } from "../ui/controls.js";
import { register, PLOT_CONFIG, loadTemplate, themeMode } from "./theme.js";
import { fmt } from "../ui/format.js";
import { gates } from "../api.js";

/**
 * opts: {title, subtitle?, tag?, height?, build: () => {data, layout} | null,
 *        table: () => {columns, rows}, empty?: Node, filename}
 */
export function chartCard(opts) {
  const titleEl = h("h3.chart-title", opts.title);
  const subEl = h("p.chart-sub", opts.subtitle || "");
  const tagSlot = h("span.chart-tag");
  const plot = h("div.chart-plot", { role: "img", "aria-label": opts.title });
  const tableHost = h("div.chart-table", { hidden: true });
  const emptyHost = h("div.chart-empty", { hidden: true });
  const withheldHost = h("div.chart-withheld", { hidden: true, role: "note" });
  let showTable = false;
  let hasData = false;
  let plotted = false;

  const tableBtn = h("button.icon-btn.chart-tool", {
    type: "button", "aria-pressed": "false", title: "Show the data as a table", "aria-label": `Table view of ${opts.title}`,
    on: { click: () => { showTable = !showTable; tableBtn.setAttribute("aria-pressed", String(showTable)); sync(); } },
  }, icon("table", { size: 18 }));
  const pngBtn = h("button.icon-btn.chart-tool", { type: "button", title: "Export PNG (2x)", "aria-label": `Export ${opts.title} as PNG`, on: { click: () => exportAs("png") } }, icon("download", { size: 18 }), h("span.tool-txt", "PNG"));
  const svgBtn = h("button.icon-btn.chart-tool", { type: "button", title: "Export SVG", "aria-label": `Export ${opts.title} as SVG`, on: { click: () => exportAs("svg") } }, h("span.tool-txt", "SVG"));

  const root = h("section.chart-card", { "aria-label": opts.title },
    h("header.chart-head",
      h("div.chart-titles", h("div.chart-titleline", titleEl, tagSlot), subEl),
      h("div.chart-tools", tableBtn, pngBtn, svgBtn)),
    h("div.chart-screen", plot, tableHost, emptyHost), withheldHost);
  if (opts.height) root.style.setProperty("--chart-h", `${opts.height}px`);

  const chart = {
    el: root,
    plot,
    setTitle(t, sub) { titleEl.textContent = t; plot.setAttribute("aria-label", t); if (sub !== undefined) subEl.textContent = sub || ""; },
    setTag(tag) { clear(tagSlot); const c = tagChip(tag, { title: "tag chain of the plotted run" }); if (c) tagSlot.append(c); },
    async render() {
      if (!window.Plotly) { await waitPlotly(); }
      await loadTemplate(themeMode());
      await gates();
      let spec = opts.build();
      // A series without a provenance tag is refused (CONTRACT rule 2): the builder lists it in
      // spec.withheld and the card says why instead of drawing it.
      const withheld = spec?.withheld || [];
      clear(withheldHost);
      withheldHost.hidden = !withheld.length;
      if (withheld.length) {
        withheldHost.append(icon("warn", { size: 14 }), h("span", `Not drawn: ${withheld.join(", ")} returned no provenance tag (tag_chain is null), so the series is withheld.`));
      }
      if (spec && !spec.data?.length && withheld.length) spec = null;
      root.dataset.withheld = String(withheld.length);
      hasData = !!spec;
      if (!spec) {
        if (plotted) { window.Plotly.purge(plot); plotted = false; }
        sync();
        return;
      }
      sync();
      await window.Plotly.react(plot, spec.data, spec.layout, PLOT_CONFIG);
      plotted = true;
      if (showTable) renderTable();
    },
    /** Phosphor persistence: draw a trace at 15% alpha that decays over 600 ms. */
    async phosphor(x, y, color) {
      if (!plotted || reducedMotion() || !window.Plotly) return;
      const P = window.Plotly;
      await P.addTraces(plot, { x, y, type: "scatter", mode: "lines", line: { color, width: 1.5 }, opacity: 0.15, hoverinfo: "skip", showlegend: false, meta: "phosphor" });
      const idx = plot.data.length - 1;
      const t0 = performance.now();
      const step = () => {
        const k = (performance.now() - t0) / 600;
        const i = plot.data.findIndex((d) => d.meta === "phosphor");
        if (i < 0) return;
        if (k >= 1) { P.deleteTraces(plot, i); return; }
        P.restyle(plot, { opacity: 0.15 * (1 - k) }, [i]);
        requestAnimationFrame(step);
      };
      void idx;
      requestAnimationFrame(step);
    },
    setEmpty(node) { clear(emptyHost); if (node) emptyHost.append(node); sync(); },
    dispose() { unregister(); ro.disconnect(); if (plotted && window.Plotly) window.Plotly.purge(plot); },
  };

  function sync() {
    const emptyShown = !hasData && emptyHost.firstChild;
    emptyHost.hidden = !emptyShown;
    plot.hidden = !hasData || showTable;
    tableHost.hidden = !hasData || !showTable;
    if (showTable && hasData) renderTable();
    safeResize();
  }

  function safeResize() {
    if (!plotted || plot.hidden || !plot.offsetParent || !plot.clientWidth) return;
    try { Promise.resolve(window.Plotly?.Plots.resize(plot)).catch(() => {}); } catch { /* not displayed yet */ }
  }

  function renderTable() {
    const t = opts.table?.();
    clear(tableHost);
    if (!t) { tableHost.append(h("p.chart-sub", "No rows.")); return; }
    tableHost.append(dataTable({ columns: t.columns.map((c) => ({ align: "right", fmt: (v) => (typeof v === "number" ? fmt(v) : v ?? "n/a"), ...c })), rows: t.rows, caption: t.caption }));
  }

  async function exportAs(format) {
    if (!plotted) { toast("Run the design first; there is nothing to export yet."); return; }
    try {
      const w = plot.clientWidth || 800, hgt = plot.clientHeight || 400;
      const url = await window.Plotly.toImage(plot, { format, scale: 2, width: w, height: hgt });
      const a = h("a", { href: url, download: `${opts.filename?.() || "fsim-chart"}.${format}` });
      document.body.append(a); a.click(); a.remove();
      toast(`Exported ${format.toUpperCase()} at 2x.`, { kind: "ok" });
    } catch (e) {
      toast(`Export failed: ${e.message}`, { kind: "error" });
    }
  }

  const unregister = register(chart);
  const ro = new ResizeObserver(() => safeResize());
  ro.observe(plot);
  if (opts.empty) chart.setEmpty(opts.empty);
  return chart;
}

function waitPlotly() {
  return new Promise((resolve) => {
    const t = setInterval(() => { if (window.Plotly) { clearInterval(t); resolve(); } }, 30);
  });
}
