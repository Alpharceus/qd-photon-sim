// Results explorer (#/results/<id>[/<tab>]): committed campaigns in out/.
// Tabs: Verdict (parsed VERDICT/BEST lines as tables) | Explore (faceted filters from the CSV's
// input axes, chart builder x / y / colour / facet with presets; rt_edge pinned to the headline
// model, other model combos only in a greyed "never gates" view; invalid rows hollow grey with a
// count; NaN as n/a + reason; commanded vs delivered paired for nanowires) | 3D (sweep lattice)
// | Figures (committed PNG/SVG, lightbox) | Provenance (manifest commit, hashes, runtime, dirty).
// Row click opens a drawer with the full CSV row and a live re-run beside the CSV value.
// No physics here: rows come from /api/campaigns/<id>/sweep, re-runs from /api/run.

import { h, clear } from "../ui/dom.js";
import { icon } from "../ui/icons.js";
import { tagChip } from "../ui/tagChip.js";
import { segmented, switchToggle, dataTable, skeleton, emptyState, errorState, toast } from "../ui/controls.js";
import { fmt, isNum } from "../ui/format.js";
import { verdictBadge, humanWord } from "../ui/verdictBadge.js";
import { sro, campaignTag, lightbox } from "../ui/figure.js";
import { chartCard } from "../charts/chartCard.js";
import { buildSweepScatter, colLabel, isG2, isFlux } from "../charts/sweep.js";
import { rtEdgeRowToDesign, rtEdgePairs, liveNanReason } from "./rowSeed.js";

const TABS = [
  { id: "verdict", label: "Verdict" },
  { id: "explore", label: "Explore" },
  { id: "3d", label: "3D" },
  { id: "figures", label: "Figures" },
  { id: "provenance", label: "Provenance" },
];

// Explore configuration per campaign: validity / reason / id columns and the chart presets.
// Presets only choose columns and axis scales; nothing is derived from the values.
const EXPLORE = {
  rt_edge: {
    valid: "diagnostic_valid", reason: "invalid_reasons_pulsed", id: "config_id",
    key: ["config_id", "model_finite_pulse", "model_tau_cap_density"], headline: true, rerun: "rt_edge",
    hideFacets: ["model_finite_pulse", "model_tau_cap_density", "config_id", "delta_xx_tag", "gamma300_tag", "irf_tag", "card_class"],
    presets: {
      g2T: { x: "T_hs_K", y: "g2_pulsed", color: "delta_xx_meV", facet: "card_id" },
      geom: { x: "emission_NA", y: "collected_flux_pulsed_s", color: "emission_R_back", facet: "T_hs_K", logY: true },
      pareto: { x: "collected_flux_pulsed_s", y: "g2_pulsed", color: "T_hs_K", facet: "card_id", logX: true },
    },
  },
  nitride_cavity: {
    valid: "valid", reason: "invalid_reasons", id: "row_id", key: ["row_id"],
    presets: {
      g2T: { x: "T_hs", y: "g2_op", color: "height_nm", facet: "regime" },
      geom: { x: "height_nm", y: "collected_flux_pulsed_s", color: "radius_nm", facet: "regime", logY: true },
      pareto: { x: "collected_flux_pulsed_s", y: "g2_op", color: "T_hs", facet: "regime", logX: true },
    },
  },
  nitride_geometry_stark: {
    valid: "valid", reason: "invalid_reasons", id: "row_id", key: ["row_id"],
    hideFacets: ["T_hs_requested"],
    presets: {
      g2T: { x: "T_hs", y: "g2_op", color: "regime", facet: "orientation" },
      geom: { x: "height_nm", y: "collected_flux_pulsed_s", color: "radius_nm", facet: "orientation", logY: true },
      pareto: { x: "collected_flux_pulsed_s", y: "g2_op", color: "screening_fraction", facet: "regime", logX: true },
    },
  },
  "nitride_nanowire/full": {
    valid: "valid", reason: "invalid_reasons", id: "row_id", key: ["row_id"], nanowire: true,
    defaults: { row_kind: ["core"] },
    hideFacets: ["card_file", "sensitivity_value", "bound_role"],
    presets: {
      g2T: { x: "T_hs", y: "g2_op", color: "family", facet: "strain_bound" },
      geom: { x: "core_radius_nm", y: ["collected_flux_pulsed_s", "collected_flux_delivered_s"], facet: "strain_bound", logY: true },
      pareto: { x: "collected_flux_delivered_s", y: "g2_op", color: "family", facet: "strain_bound", logX: true },
    },
  },
};
const PRESET_LABELS = { g2T: "g² vs T_hs", geom: "Flux vs geometry", pareto: "g²–flux Pareto", custom: "Custom" };
const LATTICE = { rt_edge: "rt_edge", nitride_cavity: "nitride_cavity", "nitride_nanowire/full": "nitride_nanowire" };

let R = null;

export function mount(el, ctx) {
  R = { el, ctx, id: null, tab: "verdict", list: [], camp: null, charts: {}, scene: null, disposers: [], ex: null, seq: 0 };
  el.classList.add("ws-results");
  R.disposers.push(ctx.state.on("theme", (m) => R?.scene?.setTheme?.(m)));
  const onKey = (e) => { if (e.key === "Escape" && R?.drawer) { e.stopPropagation(); closeDrawer(); } };
  document.addEventListener("keydown", onKey, true);
  R.disposers.push(() => document.removeEventListener("keydown", onKey, true));
  update(ctx.args);
}

export function unmount() {
  if (!R) return;
  for (const d of R.disposers) d();
  disposeBody();
  closeDrawer();
  R.el.classList.remove("ws-results");
  R = null;
}

function parseArgs(args) {
  const a = [...(args || [])];
  let tab = "verdict";
  if (a.length && TABS.some((t) => t.id === a.at(-1))) tab = a.pop();
  return { id: a.join("/") || null, tab };
}

export async function update(args) {
  if (!R) return;
  const { id, tab } = parseArgs(args);
  if (!R.list.length) {
    renderFrame();
    try { R.list = await R.ctx.api.listCampaigns(); } catch (e) {
      clear(R.el).append(h("div.ws-pad", errorState({ title: "Could not list out/ campaigns", message: e.message, actions: [{ label: "Retry", fn: () => update(args) }] })));
      return;
    }
    if (!R) return;
    renderList();
  }
  const cid = id || "rt_edge";
  if (!id) history.replaceState(null, "", `#/results/${cid}${tab !== "verdict" ? `/${tab}` : ""}`);
  const changed = cid !== R.id;
  R.id = cid;
  R.tab = tab;
  syncList();
  if (changed) {
    R.camp = null; R.ex = null; R.tagInfo = null;
    renderHead(null);
    renderTabs();
    showBodySkeleton();
    const seq = ++R.seq;
    try {
      const [camp, tagInfo] = await Promise.all([
        R.ctx.api.getJSON(`/api/campaigns/${cid}`),
        campaignTag(cid),
      ]);
      if (!R || seq !== R.seq) return;
      R.camp = camp; R.tagInfo = tagInfo;
    } catch (e) {
      if (!R || seq !== R.seq) return;
      clear(R.dom.body).append(errorState({ title: `Could not open ${cid}`, message: e.message, actions: [{ label: "Back to the board", fn: () => R.ctx.navigate("#/overview") }] }));
      return;
    }
    renderHead(R.camp);
  }
  renderTabs();
  renderBody();
}

// ---------------------------------------------------------------- frame
function renderFrame() {
  clear(R.el);
  R.dom = {
    list: h("nav.res-list", { "aria-label": "Campaigns in out/" }, skeleton("block", { lines: 8 })),
    head: h("header.res-head"),
    tabs: h("nav.res-tabs", { "aria-label": "Result views" }),
    body: h("div.res-body"),
  };
  R.el.append(h("div.results", R.dom.list, h("section.res-main", R.dom.head, R.dom.tabs, R.dom.body)));
}

function renderList() {
  const groups = [
    ["Campaigns", R.list.filter((c) => c.kind === "campaign")],
    ["Figures + CSV", R.list.filter((c) => c.kind === "collection")],
    ["Legacy F-series", R.list.filter((c) => c.kind === "legacy")],
  ];
  clear(R.dom.list).append(...groups.filter(([, l]) => l.length).map(([label, l]) => h("div.rl-group",
    h("h2.engrave.rl-h", label),
    h("ul.rl-items", l.map((c) => h("li", h("a.rl-item", { href: `#/results/${c.id}`, dataset: { id: c.id } },
      h("span.rl-id", c.id),
      c.status_word ? h("span.rl-word", { dataset: { status: c.status_word } }, c.status_word === "mixed" ? "mixed" : humanWord(c.status_word).toLowerCase()) : null,
      c.stale ? h("span.stale-dot", { title: c.stale_note || "pre-audit" }, "pre-audit") : null)))))));
}

function syncList() {
  for (const a of R.dom.list.querySelectorAll(".rl-item")) {
    if (a.dataset.id === R.id) a.setAttribute("aria-current", "page"); else a.removeAttribute("aria-current");
  }
}

function renderHead(camp) {
  const meta = R.list.find((c) => c.id === R.id) || {};
  clear(R.dom.head).append(
    h("div.rh-titles",
      h("h1.ws-title", R.id),
      h("p.ws-sub", camp?.md ? `out/${camp.md}` : meta.path || "")),
    h("div.rh-meta",
      meta.status_word && meta.status_word !== "mixed" ? verdictBadge(meta.status_word) : meta.status_word === "mixed" ? h("span.chip", `mixed: ${Object.entries(meta.status_words || {}).map(([w, n]) => `${n} ${w}`).join(", ")}`) : h("span.chip", meta.kind === "legacy" ? "legacy F-series" : "no verdict lines"),
      R.tagInfo?.tag ? h("span.rh-tag", tagChip(R.tagInfo.tag, { title: R.tagInfo.note }), h("span.rh-tag-t", "rows' widest tag")) : null,
      h("span.id-chip.is-static", { title: meta.commit_source || "" }, h("span.id-k", "DATA"), h("span.id-v", meta.commit || "n/a")),
      h("span.chip", `generated ${String(meta.generated || "").slice(0, 10)}`),
      meta.stale ? h("span.stale-badge", { title: meta.stale_note || "" }, icon("warn", { size: 14 }), `pre-audit${meta.stale_files ? `: ${meta.stale_files.join(", ")}` : ""}`) : null));
}

function renderTabs() {
  clear(R.dom.tabs).append(...TABS.map((t) => h("a.res-tab", {
    href: `#/results/${R.id}/${t.id}`, dataset: { tab: t.id },
    "aria-current": t.id === R.tab ? "page" : null,
  }, t.label)));
}

function disposeBody() {
  for (const c of Object.values(R?.charts || {})) { try { c.dispose(); } catch { /* */ } }
  if (R) R.charts = {};
  try { R?.scene?.dispose?.(); } catch { /* */ }
  if (R) R.scene = null;
}

function showBodySkeleton() {
  disposeBody();
  clear(R.dom.body).append(skeleton("chart", { lines: 0, height: 420 }));
}

function renderBody() {
  if (!R.camp) return;
  disposeBody();
  clear(R.dom.body);
  R.dom.body.dataset.tab = R.tab;
  ({ verdict: renderVerdict, explore: renderExplore, "3d": render3D, figures: renderFigures, provenance: renderProvenance })[R.tab]();
}

// ================================================================ VERDICT
function renderVerdict() {
  const c = R.camp;
  const vs = c.verdicts || [];
  const best = c.best || [];
  const body = R.dom.body;
  if (!vs.length && !best.length) {
    body.append(h("div.panel", emptyState({ title: "No VERDICT lines here", body: `${R.id} is a ${c.manifest && Object.keys(c.manifest).length ? "campaign" : "collection"} without parsed VERDICT or BEST lines. Its CSVs open in Explore and its figures in Figures.`, iconName: "info" })));
    return;
  }
  body.append(h("p.res-note", icon("info", { size: 16 }), `Parsed verbatim from out/${c.md || R.id}. Values are the printed text; nothing here is recomputed. Each count keeps the denominator the line prints.`));
  if (vs.length) body.append(lineTable("VERDICT lines", vs, true));
  if (best.length) body.append(lineTable("BEST lines", best, false));
  body.append(h("details.drawer.raw-lines", h("summary", h("span.engrave", "Raw lines"), h("span.drawer-count", `${vs.length + best.length} lines`)),
    h("div.drawer-body", [...vs, ...best].map((v) => h("pre.raw-line", `${v.line}: ${v.raw}`)))));
}

function lineTable(title, lines, isVerdict) {
  const keys = [];
  for (const v of lines) for (const k of Object.keys(v.fields || {})) if (!keys.includes(k)) keys.push(k);
  const tag = R.tagInfo?.tag;
  let table;
  if (lines.length <= 2) {
    // transposed: one row per key, one column per line
    table = h("div.table-wrap", h("table.data-table.vt-table.is-tx",
      h("thead", h("tr", h("th", { scope: "col" }, "field"), lines.map((v) => h("th", { scope: "col" }, `line ${v.line}`)))),
      h("tbody",
        h("tr", h("th", { scope: "row" }, isVerdict ? "word" : "kind"), lines.map((v) => h("td", isVerdict ? verdictBadge(v, { size: "sm" }) : v.kind))),
        keys.map((k) => h("tr", h("th", { scope: "row" }, k), lines.map((v) => h("td", v.fields?.[k] ?? "")))),
        lines.some((v) => v.note) ? h("tr", h("th", { scope: "row" }, "note"), lines.map((v) => h("td", v.note || ""))) : null)));
  } else {
    table = h("div.table-wrap.vt-wide", h("table.data-table.vt-table",
      h("thead", h("tr", h("th", { scope: "col" }, "line"), h("th", { scope: "col" }, isVerdict ? "word" : "kind"), keys.map((k) => h("th", { scope: "col" }, k)), lines.some((v) => v.note) ? h("th", { scope: "col" }, "note") : null)),
      h("tbody", lines.map((v) => h("tr",
        h("td.num", String(v.line)),
        h("td", isVerdict ? verdictBadge(v, { size: "sm" }) : v.kind),
        keys.map((k) => h("td", { class: /^[-+\d.e/:,]+$/i.test(String(v.fields?.[k] ?? "")) ? "num" : "" }, v.fields?.[k] ?? "")),
        lines.some((x) => x.note) ? h("td", v.note || "") : null)))));
  }
  return h("section.panel.vt-panel",
    h("header.panel-head", h("h2.engrave", title),
      h("span.panel-note", `${lines.length} line${lines.length === 1 ? "" : "s"} · numbers carry the rows' widest tag `, tag ? tagChip(tag, { size: "sm", title: R.tagInfo?.note }) : "(no provenance column)")),
    table);
}

// ================================================================ EXPLORE
function exploreCfg() {
  const base = EXPLORE[R.id];
  if (base) return base;
  // generic: first numeric output against first numeric input
  const cols = R.camp.columns || [];
  const num = cols.filter((c) => c.kind === "number");
  const inp = num.filter((c) => c.role === "input");
  const out = num.filter((c) => c.role === "output");
  const x = (inp[0] || num[0])?.name;
  const y = (out.find((c) => isG2(c.name)) || out[0] || num[1] || num[0])?.name;
  return { valid: cols.some((c) => c.name === "valid") ? "valid" : null, reason: cols.some((c) => c.name === "invalid_reasons") ? "invalid_reasons" : null,
    id: cols.some((c) => c.name === "row_id") ? "row_id" : null, key: cols.some((c) => c.name === "row_id") ? ["row_id"] : [],
    presets: x && y ? { custom: { x, y, color: null, facet: null } } : {} };
}

async function renderExplore() {
  const body = R.dom.body;
  const cfg = exploreCfg();
  const camp = R.camp;
  if (!camp.primary_csv) {
    body.append(h("div.panel", emptyState({ title: "No CSV in this folder", body: "This campaign only holds figures; open the Figures tab.", iconName: "table" })));
    return;
  }
  if (!R.ex) {
    const presetId = Object.keys(cfg.presets)[0];
    R.ex = { cfg, preset: presetId, ...(cfg.presets[presetId] || {}), sensitivity: false, filters: {}, rows: null, grey: [], loaded: null };
    for (const [k, v] of Object.entries(cfg.defaults || {})) R.ex.filters[k] = new Set(v);
  }
  const ex = R.ex;
  body.append(h("div.explore", { "aria-busy": "true" }, skeleton("chart", { lines: 0, height: 420 })));
  try { await loadRows(); } catch (e) {
    if (!R || R.tab !== "explore") return;
    clear(body).append(errorState({ title: "Could not read the sweep", message: e.message, actions: [{ label: "Retry", fn: () => { R.ex.loaded = null; renderBody(); } }] }));
    return;
  }
  if (!R || R.tab !== "explore") return;
  clear(body);
  const filters = h("aside.ex-filters", { "aria-label": "Filters" });
  const tools = h("div.ex-tools");
  const count = h("p.ex-count", { "aria-live": "polite" });
  R.dom.ex = { filters, tools, count };
  const chart = chartCard({
    title: "Sweep", height: 440,
    build: () => buildChart(),
    table: () => tableRows(),
    filename: () => `${R.id.replace(/\//g, "_")}-${ex.preset}`,
  });
  R.charts.ex = chart;
  chart.el.classList.add("ex-chart");
  const rowsHost = h("div.ex-rows");
  R.dom.ex.rowsHost = rowsHost;
  body.append(h("div.explore", filters, h("div.ex-main", tools, count, chart.el, rowsHost)));
  renderFilters();
  renderTools();
  await redraw();
  chart.plot.on?.("plotly_click", (ev) => { const i = ev?.points?.[0]?.customdata; if (isNum(i)) openRow(ex.view[i]); });
}

function neededCols() {
  const ex = R.ex, cfg = ex.cfg;
  const cols = new Set();
  for (const c of R.camp.columns || []) if (c.role === "input") cols.add(c.name);
  for (const p of Object.values(cfg.presets)) for (const k of ["x", "color", "facet"]) if (p[k]) cols.add(p[k]);
  for (const p of Object.values(cfg.presets)) for (const y of [].concat(p.y || [])) cols.add(y);
  for (const k of ["x", "color", "facet"]) if (ex[k]) cols.add(ex[k]);
  for (const y of [].concat(ex.y || [])) cols.add(y);
  for (const k of [cfg.valid, cfg.reason, cfg.id, ...(cfg.key || [])]) if (k) cols.add(k);
  const known = new Set((R.camp.columns || []).map((c) => c.name));
  return [...cols].filter((c) => known.has(c));
}

async function loadRows() {
  const ex = R.ex;
  const cols = neededCols();
  const sig = `${cols.join(",")}|${ex.sensitivity}`;
  if (ex.loaded === sig) return;
  const headline = ex.cfg.headline && !ex.sensitivity ? 1 : 0;
  const file = R.camp.primary_csv.split("/").pop();
  const r = await R.ctx.api.getJSON(`/api/campaigns/${R.id}/sweep?file=${encodeURIComponent(file)}&cols=${encodeURIComponent(cols.join(","))}&limit=50000&headline=${headline}`);
  const objs = r.rows.map((row) => Object.fromEntries(r.columns.map((c, i) => [c, row[i]])));
  if (ex.cfg.headline && ex.sensitivity) {
    const isHead = (o) => o.model_finite_pulse === true && o.model_tau_cap_density === false;
    ex.rows = objs.filter(isHead);
    ex.grey = objs.filter((o) => !isHead(o));
  } else {
    ex.rows = objs; ex.grey = [];
  }
  ex.total = r.total;
  ex.filtered = r.filtered;
  ex.columns = r.columns;
  ex.loaded = sig;
  ex.rows.forEach((o, i) => { o.__i = i; });
  // facet candidates: input columns with 2..12 distinct values
  const hide = new Set([...(ex.cfg.hideFacets || []), ex.cfg.id, ...(ex.cfg.key || [])]);
  const facetCols = (R.camp.columns || []).filter((c) => c.role === "input" && !hide.has(c.name) && r.columns.includes(c.name));
  ex.facets = facetCols.map((c) => {
    const vals = [...new Set(ex.rows.map((o) => (o[c.name] == null ? "n/a" : String(o[c.name]))))];
    vals.sort((a, b) => (isNum(Number(a)) && isNum(Number(b)) ? Number(a) - Number(b) : a.localeCompare(b)));
    return { name: c.name, values: vals };
  }).filter((f) => f.values.length >= 2 && f.values.length <= 12);
}

function applyFilters(rows) {
  const f = R.ex.filters;
  const active = Object.entries(f).filter(([, s]) => s && s.size);
  if (!active.length) return rows;
  return rows.filter((o) => active.every(([k, s]) => s.has(o[k] == null ? "n/a" : String(o[k]))));
}

function renderFilters() {
  const ex = R.ex;
  const host = R.dom.ex.filters;
  const items = [];
  if (ex.cfg.headline) {
    items.push(h("section.ex-fgroup.ex-model",
      h("h3.engrave", "Model"),
      h("span.headline-chip", { title: "model_finite_pulse = True, model_tau_cap_density = False" }, icon("pin", { size: 14 }), "headline model"),
      h("p.ex-fnote", "Only the headline model gates PASS. The other three combos render greyed, never coloured."),
      switchToggle({ label: "model sensitivity (never gates)", checked: ex.sensitivity, onChange: async (on) => { ex.sensitivity = on; ex.loaded = null; await loadRows(); renderFilters(); redraw(); } })));
  }
  for (const f of ex.facets) {
    const sel = ex.filters[f.name];
    items.push(h("section.ex-fgroup", { dataset: { facet: f.name } },
      h("div.ex-fhead", h("h3.engrave", colLabel(f.name)),
        sel && sel.size ? h("button.btn-link", { type: "button", on: { click: () => { delete ex.filters[f.name]; renderFilters(); redraw(); } } }, "all") : null),
      h("div.ex-chips", f.values.map((v) => h("button.fchip", {
        type: "button", "aria-pressed": String(!!sel?.has(v)),
        on: {
          click: () => {
            const s = ex.filters[f.name] || new Set();
            if (s.has(v)) s.delete(v); else s.add(v);
            if (s.size) ex.filters[f.name] = s; else delete ex.filters[f.name];
            renderFilters(); redraw();
          },
        },
      }, fmtCat(v))))));
  }
  if (!ex.facets.length) items.push(h("p.ex-fnote", "No input axis with 2 to 12 values to filter on."));
  clear(host).append(h("h2.engrave.ex-ftitle", "Filters"), ...items);
}

const fmtCat = (v) => (isNum(Number(v)) && v !== "" ? fmt(Number(v)) : String(v).replace(/_/g, " "));

function renderTools() {
  const ex = R.ex;
  const cols = (ex.columns || []).filter((c) => !c.startsWith("__"));
  const numCols = cols.filter((c) => ex.rows.some((o) => isNum(o[c])));
  const catCols = ex.facets.map((f) => f.name);
  const presets = Object.keys(ex.cfg.presets);
  const seg = presets.length > 1 ? segmented({
    label: "Chart preset", value: ex.preset,
    options: [...presets.map((p) => ({ value: p, label: PRESET_LABELS[p] || p })), { value: "custom", label: "Custom", disabled: true }],
    onChange: (v) => { Object.assign(ex, { x: null, y: null, color: null, facet: null, logX: false, logY: false }, ex.cfg.presets[v]); ex.preset = v; renderTools(); reloadThenDraw(); },
  }) : null;
  if (seg && ex.preset === "custom") seg.select("custom");
  const sel = (label, key, options, allowNone) => h("label.ex-sel", h("span.ctx-k", label),
    h("select.fr-select", {
      on: {
        change: (e) => {
          const v = e.target.value || null;
          if (key === "y" && Array.isArray(ex.y) && v === "__pair") return;
          ex[key] = v;
          if (key === "x") ex.logX = isFlux(v);
          if (key === "y") ex.logY = isFlux(v);
          ex.preset = "custom";
          renderTools(); reloadThenDraw();
        },
      },
    },
    allowNone ? h("option", { value: "" }, "none") : null,
    key === "y" && Array.isArray(ex.y) ? h("option", { value: "__pair", selected: true }, "commanded + delivered (paired)") : null,
    options.map((c) => h("option", { value: c, selected: !Array.isArray(ex[key]) && ex[key] === c }, colLabel(c, { nanowire: ex.cfg.nanowire })))));
  const logs = h("div.ex-logs",
    switchToggle({ label: "log x", checked: !!ex.logX, onChange: (v) => { ex.logX = v; redraw(); } }),
    switchToggle({ label: "log y", checked: !!ex.logY, onChange: (v) => { ex.logY = v; redraw(); } }));
  clear(R.dom.ex.tools).append(...[
    seg ? seg.el : null,
    h("div.ex-sels",
      sel("x", "x", numCols, false),
      sel("y", "y", numCols, false),
      sel("colour", "color", catCols.concat(numCols.filter((c) => !catCols.includes(c))), true),
      sel("facet", "facet", catCols, true)),
    logs].filter(Boolean));
}

async function reloadThenDraw() {
  try { await loadRows(); } catch (e) { toast(`Could not read the sweep: ${e.message}`, { kind: "error" }); return; }
  renderFilters();
  redraw();
}

/**
 * The campaign's own flux floor (N3): the flux_floor field its VERDICT lines print
 * (e.g. "flux_floor=1000/s"), cited as out/<md>:<line>. Null when no line prints one, or when the
 * lines print different floors (no single floor to draw; the text says so).
 */
function campaignFluxFloor(camp) {
  const lines = (camp?.verdicts || []).filter((v) => v.fields?.flux_floor != null);
  if (!lines.length) return null;
  const parse = (t) => { const m = /^\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)/.exec(String(t)); return m ? Number(m[1]) : NaN; };
  const vals = lines.map((v) => parse(v.fields.flux_floor));
  if (!vals.every(isNum)) return null;
  const md = `out/${camp.md || camp.id}`;
  if (new Set(vals).size > 1) return { conflict: true, source: `${md} (VERDICT flux_floor differs between lines)` };
  return { value: vals[0], raw: lines[0].fields.flux_floor, source: `${md}:${lines[0].line} (VERDICT flux_floor=${lines[0].fields.flux_floor})` };
}

/**
 * Gates: g2 ceiling from GET /api/gates (I3); flux floor from the campaign's own VERDICT
 * flux_floor when it prints one, else from GET /api/gates. Every value carries its source.
 */
async function gateRefs() {
  const g = await R.ctx.api.gates();
  const own = campaignFluxFloor(R.camp);
  if (!g && !(own && !own.conflict)) return null;
  const c = g?.g2_ceiling || {}, f = g?.flux_floor || {};
  const tag = (t) => (t ? String(t).replace(/[[\]]/g, "") : null);
  const useOwn = own && !own.conflict;
  return {
    ceiling: isNum(c.value) ? c.value : null, ceilingTag: tag(c.tag), ceilingSource: c.source || null,
    floor: useOwn ? own.value : isNum(f.value) ? f.value : null,
    // a VERDICT line prints no bracket tag next to flux_floor: the gate convention tag stands
    floorTag: tag(f.tag) || "A",
    floorSource: useOwn ? own.source : f.source ? `${f.source} (GET /api/gates; ${own?.conflict ? own.source : "this campaign prints no flux_floor"})` : null,
    floorOrigin: useOwn ? "campaign" : "gates",
  };
}

async function redraw() {
  if (!R?.ex || !R.charts.ex) return;
  const ex = R.ex;
  ex.view = applyFilters(ex.rows);
  ex.view.forEach((o, i) => { o.__v = i; });
  ex.gates = await gateRefs();
  if (!R?.charts?.ex) return;
  const grey = applyFilters(ex.grey);
  ex.greyView = grey;
  const b = buildChart();
  const st = b?.stats || { shown: 0, invalid: 0, nan: 0, nonpos: 0, reasons: new Map() };
  const parts = [];
  parts.push(h("span.ex-n", { dataset: { rows: String(ex.view.length), total: String(ex.rows.length) } },
    `${ex.rows.length} ${ex.cfg.headline ? "headline " : ""}rows`));
  if (ex.view.length !== ex.rows.length) parts.push(` · ${ex.view.length} match the filters`);
  parts.push(` · ${st.shown} plotted`);
  if (st.invalid) parts.push(h("span.ex-inv", ` · ${st.invalid} invalid rows: ${st.invalidDrawn} drawn hollow grey${st.invalid - st.invalidDrawn ? `, ${st.invalid - st.invalidDrawn} without values` : ""}`));
  if (st.nan) parts.push(h("span.ex-nan", { title: [...st.reasons.entries()].map(([r, n]) => `${n} × ${r}`).join("\n") }, ` · ${st.nan} n/a: ${[...st.reasons.keys()][0] || "value missing"}${st.reasons.size > 1 ? ` (+${st.reasons.size - 1} more reasons)` : ""}`));
  if (st.nonpos) parts.push(` · ${st.nonpos} at ≤ 0 hidden by the log axis`);
  if (st.above1) parts.push(h("span.ex-above", ` · ${st.above1} above g² = 1: axis extended to ${fmt(st.yTop)}${st.offcanvas ? `, ${st.offcanvas} still off-canvas (hidden)` : ", none hidden"}`));
  if (ex.gates && isNum(ex.gates.floor)) parts.push(h("span.ex-floor-src", { dataset: { origin: ex.gates.floorOrigin } }, ` · flux floor ${fmt(ex.gates.floor)} s⁻¹ from ${ex.gates.floorSource || "unknown source"}`));
  if (!ex.gates) parts.push(h("span.ex-nogate", ` · gates unavailable (${R.ctx.api.gatesProblem?.() || "GET /api/gates failed"}): no ceiling, floor or gate box drawn`));
  if (grey.length) parts.push(h("span.ex-grey", ` · ${grey.length} other-model rows greyed (never gate)`));
  if (st.droppedFacets > 0) parts.push(` · facet shows the first 4 of ${st.facetValues.length} values`);
  clear(R.dom.ex.count).append(...parts);
  const ys = [].concat(ex.y || []);
  R.charts.ex.setTitle(`${ys.map((y) => colLabel(y, { nanowire: ex.cfg.nanowire })).join(" and ")} against ${colLabel(ex.x, { nanowire: ex.cfg.nanowire })}`,
    `${R.camp.primary_csv}${ex.cfg.headline ? (ex.sensitivity ? " · headline coloured, other models greyed" : " · headline model only") : ""}${ex.facet ? ` · facet ${ex.facet}` : ""}`);
  R.charts.ex.setTag(R.tagInfo?.tag);
  R.charts.ex.setEmpty(emptyState({ title: "No rows match", body: "Every row is filtered out. Clear a filter on the left.", iconName: "table" }));
  await R.charts.ex.render();
  renderRowTable();
}

function buildChart() {
  const ex = R?.ex;
  if (!ex?.view || !ex.view.length) return null;
  const isPair = Array.isArray(ex.y);
  const res = buildSweepScatter({
    rows: ex.view, x: ex.x, y: ex.y, color: isPair ? null : ex.color, facet: ex.facet,
    logX: ex.logX, logY: ex.logY, validKey: ex.cfg.valid, reasonKey: ex.cfg.reason, idKey: ex.cfg.id,
    greyRows: ex.greyView || [], tag: R.tagInfo?.tag, nanowire: ex.cfg.nanowire,
    refs: {
      ceiling: ex.gates?.ceiling, ceilingTag: ex.gates?.ceilingTag, floor: ex.gates?.floor, floorTag: ex.gates?.floorTag,
      gate: isNum(ex.gates?.ceiling) && isNum(ex.gates?.floor) ? { g2: ex.gates.ceiling, flux: ex.gates.floor } : null,
    },
  });
  return res;
}

function tableRows() {
  const ex = R?.ex;
  if (!ex?.view) return null;
  const ys = [].concat(ex.y || []);
  const cols = [ex.cfg.id, ex.x, ...ys, ex.color, ex.facet, ex.cfg.valid].filter((c, i, a) => c && a.indexOf(c) === i);
  return { caption: `${ex.view.length} rows (filters applied)`, columns: cols.map((c) => ({ key: c, label: c })), rows: ex.view.slice(0, 2000) };
}

function renderRowTable() {
  const ex = R.ex;
  const host = R.dom.ex.rowsHost;
  const ys = [].concat(ex.y || []);
  const cols = [ex.cfg.id, ex.x, ...ys, ex.color, ex.facet, ex.cfg.valid, ex.cfg.reason].filter((c, i, a) => c && a.indexOf(c) === i);
  const shown = ex.view.slice(0, 60);
  const table = h("table.data-table.ex-table",
    h("caption", `First ${shown.length} of ${ex.view.length} rows. Click a row for the full CSV row${ex.cfg.rerun ? " and a live re-run" : ""}.`),
    h("thead", h("tr", cols.map((c) => h("th", { scope: "col", class: typeof shown[0]?.[c] === "number" ? "num" : "" }, colLabel(c, { nanowire: ex.cfg.nanowire }))))),
    h("tbody", shown.map((o) => h("tr.ex-row", {
      tabindex: "0", dataset: { i: String(o.__v) }, class: o[ex.cfg.valid] === false ? "is-invalid" : "",
      on: { click: () => openRow(o), keydown: (e) => { if (e.key === "Enter") openRow(o); } },
    }, cols.map((c) => {
      const v = o[c];
      return h("td", { class: typeof v === "number" ? "num" : "" }, v == null ? h("span.nan", "n/a") : typeof v === "number" ? fmt(v) : String(v).slice(0, 60));
    })))));
  clear(host).append(h("div.panel.ex-tablepanel", h("div.table-wrap", table)));
}

// ---------------------------------------------------------------- row drawer
async function openRow(o) {
  if (!o) return;
  closeDrawer();
  const ex = R.ex;
  const cfg = ex.cfg;
  const keyCols = (cfg.key || []).filter((k) => o[k] !== undefined);
  const idText = keyCols.map((k) => `${k}=${o[k]}`).join(" · ") || `row ${o.__v}`;
  const body = h("div.rd-body", skeleton("block", { lines: 6 }));
  const live = h("div.rd-live");
  const actions = h("div.rd-actions");
  const drawer = h("aside.row-drawer", { role: "dialog", "aria-label": `Sweep row ${idText}`, tabindex: "-1" },
    h("header.rd-head",
      h("div", h("h2.rd-title", o[cfg.id] != null ? String(o[cfg.id]) : "Sweep row"), h("p.rd-sub", idText)),
      h("button.icon-btn", { type: "button", "aria-label": "Close row", on: { click: closeDrawer } }, icon("close", { size: 18 }))),
    actions, live, body);
  document.body.append(drawer);
  R.drawer = drawer;
  drawer.focus();
  // full row from the CSV
  let full = null;
  try {
    const filt = Object.fromEntries(keyCols.map((k) => [k, o[k]]));
    const file = R.camp.primary_csv.split("/").pop();
    const r = await R.ctx.api.getJSON(`/api/campaigns/${R.id}/sweep?file=${encodeURIComponent(file)}&filter=${encodeURIComponent(JSON.stringify(filt))}&limit=2`);
    if (r.rows.length) full = Object.fromEntries(r.columns.map((c, i) => [c, r.rows[0][i]]));
  } catch (e) {
    if (R?.drawer === drawer) clear(body).append(errorState({ title: "Could not read the row", message: e.message }));
    return;
  }
  if (R?.drawer !== drawer) return;
  full = full || o;
  const card = full.card_id || (full.card_file ? String(full.card_file).replace(/\.yaml$/, "") : null);
  clear(actions).append(
    cfg.rerun === "rt_edge" && card ? h("button.btn.btn-outline", { type: "button", on: { click: () => rerunRow(full, card, live) } }, icon("retry", { size: 16 }), "Re-run this row live") : null,
    card ? h("button.btn.btn-quiet", {
      type: "button",
      on: { click: () => { R.ctx.state.set({ rerunSeed: { card, campaign: R.id, row: full, at: Date.now() } }); closeDrawer(); R.ctx.navigate(`#/design/${card}`); } },
    }, "Open card in Designer", icon("chevronRight", { size: 16 })) : null);
  if (cfg.rerun !== "rt_edge") {
    live.append(h("p.rd-note", icon("info", { size: 14 }), "Live re-run is wired for rt_edge rows only: their grid maps one-to-one onto card fields (scripts/run_rt_edge.py). This campaign's runner applies geometry presets the Studio does not replicate, so no live value is shown rather than an unmatched one."));
  }
  const tag = R.tagInfo?.tag;
  clear(body).append(h("dl.rd-dl", Object.entries(full).filter(([k]) => !k.startsWith("__")).map(([k, v]) => h("div.drow",
    h("dt", k),
    h("dd", v == null ? h("span.nan", full[`${k}__nan_reason`] || "n/a in the CSV") : typeof v === "number" ? fmt(v) : String(v))))),
  h("p.rd-foot", "Values from ", h("code", R.camp.primary_csv), " (rounded for display; hover the chart for full precision). Row tag ", tag ? tagChip(tag, { size: "sm" }) : "none recorded", "."));
}

function closeDrawer() {
  if (!R?.drawer) return;
  R.drawer.remove();
  R.drawer = null;
}

/** rt_edge row -> card fields (configuration only, same overrides as scripts/run_rt_edge.py). */
async function rerunRow(row, card, host) {
  clear(host).append(h("div.rd-cmp.is-pending", h("p.rd-note", "Evaluating the row live (single T_hs point, about 0.3 s per edge card)…")));
  try {
    const info = await R.ctx.api.getCard(card);
    const d = rtEdgeRowToDesign(info.design, row);
    const job = R.ctx.api.runJob({ card, design: d, mode: "point", T_grid: [row.T_hs_K] });
    const r = await job.promise;
    if (!R?.drawer) return;
    const s = r.result.scalars || {};
    const tag = r.result.tag_chain;
    const pairs = rtEdgePairs(row, s);
    const labels = r.result.labels || [];
    clear(host).append(h("section.rd-cmp",
      h("h3.engrave", "CSV value vs live value"),
      labels.length ? h("p.rd-labels", labels.map((l) => h("span.label-chip", l))) : null,
      h("table.data-table.cmp-table",
        h("thead", h("tr", h("th", "quantity"), h("th.num", "CSV (committed)"), h("th.num", "live now"), h("th.num", "match"))),
        h("tbody", pairs.map(([l, a, b, u, key]) => h("tr",
          h("th", { scope: "row" }, l),
          h("td.num", sro({ label: "csv", value: a, unit: u, tag: R.tagInfo?.tag, size: "sm" })),
          h("td.num", sro({ label: labels.length ? `live · ${labels[0]}` : "live", value: b, unit: u, tag, size: "sm", caption: b == null ? liveNanReason(s, key) : "" })),
          h("td.num.match", a === b ? h("span.m-ok", icon("check", { size: 14 }), "bit-identical") : isNum(a) && isNum(b) ? h("span.m-diff", `Δ ${fmt(b - a)}`) : "n/a"))))),
      h("p.rd-note", `${r.result.run_id} · fsim_core.device.evaluate at T_hs = ${fmt(row.T_hs_K)} K with the row's grid values and model switches.`)));
  } catch (e) {
    if (R?.drawer) clear(host).append(errorState({ title: "Live re-run failed", message: e.message }));
  }
}

// ================================================================ 3D lattice
async function render3D() {
  const body = R.dom.body;
  const lid = LATTICE[R.id];
  if (!lid) {
    body.append(h("div.panel", emptyState({ title: "No lattice view for this folder", body: "The sweep lattice is built for rt_edge, nitride_cavity and nitride_nanowire/full, whose CSVs are factorial grids.", iconName: "results" })));
    return;
  }
  if (!R.lat || R.lat.id !== R.id) R.lat = { id: R.id, x: null, y: null, z: null, headline: true };
  const L = R.lat;
  const slot = h("div.lat-viewport", { role: "region", "aria-label": "Sweep lattice, 3D" });
  const pickers = h("div.lat-pickers");
  body.append(h("div.lat", h("div.lat-bar", pickers), slot));
  const load = async () => {
    const qs = { campaign: lid, headline: L.headline ? 1 : 0 };
    if (L.x) qs.x = L.x; if (L.y) qs.y = L.y; if (L.z) qs.z = L.z;
    let spec;
    try { spec = await R.ctx.api.sceneSpec("lattice", qs); } catch (e) {
      toast(`Lattice: ${e.message}`, { kind: "error" });
      return;
    }
    if (!R || R.tab !== "3d") return;
    const ax = spec.overlays?.axes || {};
    L.x = ax.x?.name; L.y = ax.y?.name; L.z = ax.z?.name;
    renderPickers(spec);
    if (R.scene) { try { R.scene.update(spec); return; } catch { /* remount */ } }
    try {
      const mod = await import("../viz3d/index.js");
      if (!R || R.tab !== "3d") return;
      clear(slot);
      R.scene = mod.mountScene(slot, "lattice", spec, { theme: R.ctx.state.get("theme") });
    } catch (e) {
      clear(slot).append(h("div.vp-placeholder", h("div.vp-plate", h("span.engrave", "3D view unavailable"), h("span.vp-sub", e.message))));
    }
  };
  function renderPickers(spec) {
    const cands = (R.camp.columns || []).filter((c) => c.role === "input" && c.kind === "number").map((c) => c.name);
    const mk = (axis) => h("label.ex-sel", h("span.ctx-k", `${axis} axis`), h("select.fr-select", {
      on: { change: (e) => { L[axis] = e.target.value; load(); } },
    }, cands.map((c) => h("option", { value: c, selected: L[axis] === c }, colLabel(c)))));
    clear(pickers).append(...[mk("x"), mk("y"), mk("z"),
      R.id === "rt_edge" ? switchToggle({ label: "headline model only", checked: L.headline, onChange: (v) => { L.headline = v; load(); } }) : null,
      latBanner(spec.overlays?.banner)].filter(Boolean));
  }
  clear(slot).append(h("div.vp-placeholder", h("div.vp-plate", h("span.engrave", "Lattice loading"), h("span.vp-sub", R.id))));
  await load();
}

/** Lattice banner from the SceneSpec's banner fields (I4): gating passes only; others "never gates". */
function latBanner(b) {
  if (!b) return null;
  if (typeof b === "string") return h("span.lat-banner", b);
  const main = b.text || (isNum(b.pass) && isNum(b.gating_rows) ? `${b.pass} / ${b.gating_rows} gating rows pass` : "");
  const ng = b.non_gating_passes;
  const ngText = ng == null ? "" : typeof ng === "object" ? (ng.text || `${ng.count ?? ""} other-model passes: never gates`) : `${ng} other-model passes: never gates`;
  return h("span.lat-banner", { dataset: { gating: String(b.pass ?? ""), style: b.style || "" } }, main, ngText ? h("span.lat-ng", ` · ${ngText}`) : null);
}

// ================================================================ FIGURES
function renderFigures() {
  const figs = (R.camp.files?.figures || []).filter((f) => /\.(png|svg|jpe?g|webp)$/i.test(f));
  const pdfs = (R.camp.files?.figures || []).filter((f) => /\.pdf$/i.test(f));
  const body = R.dom.body;
  if (!figs.length) {
    body.append(h("div.panel", emptyState({ title: "No committed figures", body: pdfs.length ? `${pdfs.length} PDF file(s) only; open them from out/.` : "This folder has CSVs only; open Explore.", iconName: "chart" })));
    return;
  }
  // prefer PNG when the same figure also exists as SVG (one tile per figure)
  const stems = new Map();
  for (const f of figs) { const s = f.replace(/\.[^.]+$/, ""); if (!stems.has(s) || /\.png$/i.test(f)) stems.set(s, f); }
  const items = [...stems.values()].map((f) => ({ src: `/out/${f}`, caption: f }));
  body.append(
    h("p.res-note", icon("info", { size: 16 }), `${items.length} committed figure${items.length === 1 ? "" : "s"} from out/${R.id}. These are the files as committed (matplotlib); the Explore tab redraws the data.`),
    h("div.fig-grid", items.map((it, i) => h("button.fig-tile", { type: "button", on: { click: () => lightbox(items, i) }, "aria-label": `Open ${it.caption}` },
      h("span.fig-img", h("img", { src: it.src, alt: it.caption, loading: "lazy" })),
      h("span.fig-cap", it.caption.split("/").pop())))));
}

// ================================================================ PROVENANCE
function renderProvenance() {
  const m = R.camp.manifest || {};
  const meta = R.list.find((c) => c.id === R.id) || {};
  const body = R.dom.body;
  const commit = m.generation_commit || m.report_only_generation_commit || null;
  const dirtyKeys = Object.keys(m).filter((k) => /dirty/.test(k));
  const runtimeKeys = Object.keys(m).filter((k) => /runtime_s|runtime_seconds/.test(k));
  const hashKeys = Object.keys(m).filter((k) => /hash/.test(k) && m[k] && typeof m[k] === "object");
  const kv = (k, v) => h("div.insp-row", h("dt", k), h("dd", v));
  const facts = h("dl.insp-dl.prov-facts",
    kv("Folder", meta.path || `out/${R.id}`),
    kv("Data commit", h("span", h("code", meta.commit || "n/a"), ` (${meta.commit_source || "unknown source"})`)),
    commit ? kv("Generation commit", h("code", commit)) : null,
    kv("Generated", String(meta.generated || m.generated_utc || "n/a")),
    ...runtimeKeys.map((k) => kv(k, isNum(m[k]) ? `${fmt(m[k])} s` : String(m[k]))),
    ...dirtyKeys.map((k) => kv(k, m[k] === true ? h("span.dirty-flag", icon("warn", { size: 14 }), "true: generated from a dirty worktree") : String(m[k]))),
    m.complete != null ? kv("complete", String(m.complete)) : null,
    m.quick != null ? kv("quick", String(m.quick)) : null,
    m.versions || m.runtime_versions ? kv("versions", Object.entries(m.versions || m.runtime_versions).map(([a, b]) => `${a} ${b}`).join(" · ")) : null,
    meta.stale ? kv("Status", h("span.stale-badge", icon("warn", { size: 14 }), meta.stale_note || "pre-audit")) : null);
  body.append(h("div.prov",
    h("section.panel.prov-main", h("header.panel-head", h("h2.engrave", "Manifest")), Object.keys(m).length ? h("div.pad", facts) : h("div.pad", facts, h("p.res-note", "No manifest.json in this folder: the commit is the last git commit touching it."))),
    ...hashKeys.map((k) => h("section.panel.prov-hash", h("header.panel-head", h("h2.engrave", k), h("span.panel-note", "as recorded in manifest.json; the Studio does not re-hash files")),
      h("div.table-wrap", dataTable({ columns: [{ key: "f", label: "file" }, { key: "v", label: "sha256" }], rows: Object.entries(m[k]).map(([f, v]) => ({ f, v: h("code.hash", typeof v === "string" ? v : JSON.stringify(v)) })) })))),
    h("details.drawer", h("summary", h("span.engrave", "Full manifest.json"), h("span.drawer-count", `${Object.keys(m).length} keys`)),
      h("div.drawer-body", h("pre.raw-line.manifest-json", JSON.stringify(m, null, 2).slice(0, 60000))))));
}

