// Library workspace (#/library, #/library/<card>): the literature parameter cards.
//   index  the parameter cards grouped by paper: citation, [V]/[DR]/[E]/[A] mix, datasets + row counts
//   card   header (paper, device, source, class), parameter table (ranges as bracketed spans, never
//          midpoints), every dataset plotted (ink markers, error bars, upper-bound glyphs) with a model
//          overlay only where the physics rules file allows one (GET /api/params/<name>/overlay/<ds>,
//          which calls the named fsim_core function), "Open as design" where a design card derives
//          from the card, and "Edit copy" (PUT /api/params/<name>-edited, the legacy editor schema).
// No physics here: every number comes from GET /api/params*.

import { h, clear } from "../ui/dom.js";
import { icon } from "../ui/icons.js";
import { tagChip, normTag, tagDash, tagWidth, TAG_MEANING } from "../ui/tagChip.js";
import { skeleton, errorState, emptyState, toast } from "../ui/controls.js";
import { fmt, unitText, isNum } from "../ui/format.js";
import { chartCard } from "../charts/chartCard.js";
import { baseLayout, colors as themeColors, alpha } from "../charts/theme.js";

const TAGS = ["V", "DR", "E", "A"];
const SYMBOLS = ["circle", "square", "diamond", "triangle-up", "pentagon", "star"];

let L = null;

export function mount(el, ctx) {
  L = { el, ctx, key: null, charts: [], seq: 0, editing: false };
  el.classList.add("ws-library");
  update(ctx.args);
}

export function unmount() {
  if (!L) return;
  dispose();
  L.el.classList.remove("ws-library");
  L = null;
}

function dispose() {
  for (const c of L.charts) { try { c.dispose(); } catch { /* */ } }
  L.charts = [];
}

export function update(args) {
  if (!L) return;
  const name = args?.[0] || null;
  const key = name || "";
  if (key === L.key && !L.force) return;
  L.key = key;
  L.force = false;
  L.editing = false;
  L.seq++;
  dispose();
  clear(L.el);
  L.el.dataset.ready = "0";
  if (name) cardPage(name); else indexPage();
}

function reload() { if (L) { L.force = true; update(L.key ? [L.key] : []); } }

// ---------------------------------------------------------------- shared bits
function tagMix(counts, { label = "parameters" } = {}) {
  const total = TAGS.reduce((s, t) => s + (counts?.[t] || 0), 0);
  return h("div.lib-mix", { role: "group", "aria-label": `Provenance mix of ${total} ${label}` },
    TAGS.map((t) => h("span.lib-mix-cell", { dataset: { tag: t, n: String(counts?.[t] || 0) }, class: counts?.[t] ? "" : "is-zero" },
      tagChip(t, { size: "sm", title: `${counts?.[t] || 0} ${label} tagged [${t}]` }),
      h("span.lib-mix-n.lib-count", String(counts?.[t] || 0)))));
}

/** Range as a bracketed span: [ lo , hi ] with a bracket bar under it. Never a midpoint. */
function rangeSpan(lo, hi) {
  return h("span.lib-range", { title: `range ${fmt(lo)} to ${fmt(hi)} (swept or fitted, never averaged)` },
    h("span.rb.rb-l", { "aria-hidden": "true" }),
    h("span.lib-num", fmt(lo)), h("span.rb-sep", { "aria-hidden": "true" }, ","), h("span.lib-num", fmt(hi)),
    h("span.rb.rb-r", { "aria-hidden": "true" }),
    h("span.sr-only", ` range ${fmt(lo)} to ${fmt(hi)}`));
}

function valueCell(row) {
  if (Array.isArray(row.range)) return rangeSpan(row.range[0], row.range[1]);
  if (isNum(row.value)) return h("span.lib-num", fmt(row.value));
  return h("span.lib-text", row.text ?? "n/a");
}

function clsChip(cls) {
  if (!cls) return null;
  const code = cls.code ? ` (${cls.code})` : "";
  return h("span.label-chip.lib-cls", { dataset: { cls: cls.code || "none" }, title: cls.reason || "" }, `${cls.label}${code}`);
}

// ---------------------------------------------------------------- index
async function indexPage() {
  const seq = L.seq;
  const body = h("div.lib-body", skeleton("block", { lines: 6, height: 420 }));
  L.el.append(h("div.library",
    h("header.ws-head.lib-head",
      h("div", h("h1.ws-title", "Library"),
        h("p.ws-sub", "Literature parameter cards: what each paper published, with the provenance of every number.")),
      legend()),
    body));
  let list;
  try { list = await L.ctx.api.getJSON("/api/params"); } catch (e) {
    if (!L || seq !== L.seq) return;
    clear(body).append(errorState({ title: "Could not list the parameter cards", message: e.message }));
    return;
  }
  if (!L || seq !== L.seq) return;
  const groups = new Map();
  for (const c of list) {
    const p = c.paper || { id: "other", label: "Other cards", cite: "", order: 99 };
    if (!groups.has(p.id)) groups.set(p.id, { paper: p, cards: [] });
    groups.get(p.id).cards.push(c);
  }
  const ordered = [...groups.values()].sort((x, y) => (x.paper.order ?? 99) - (y.paper.order ?? 99));
  clear(body).append(...ordered.map(({ paper, cards }) => h("section.lib-group", { "aria-label": paper.label },
    h("header.lib-group-head", h("h2.engrave", paper.label), h("span.lib-cite", paper.cite)),
    h("div.lib-plaques", cards.map(plaque)))));
  L.el.dataset.ready = "1";
}

function legend() {
  return h("div.ov-legend.lib-legend", TAGS.map((t) => h("span.lib-leg", tagChip(t, { size: "sm" }), h("span.ov-legend-txt", TAG_MEANING[t]))));
}

function plaque(c) {
  if (c.error) {
    return h("div.cplaque.is-missing.lib-plaque", h("div.cp-top", h("h3.cp-title", c.name)), h("p.cp-sub", c.error));
  }
  return h("a.cplaque.lib-plaque", { href: `#/library/${encodeURIComponent(c.name)}`, dataset: { card: c.name } },
    h("div.cp-top",
      h("div", h("h3.cp-title", c.name), h("p.cp-sub", c.role || "")),
      h("span.cp-go", icon("chevronRight", { size: 18 }))),
    c.device ? h("p.lib-device", c.device) : null,
    h("div.lib-plaque-row",
      ...(c.n_params || !c.n_device_fields
        ? [h("span.cp-rail-k", `${c.n_params} params`), tagMix(c.tag_counts)]
        : [h("span.cp-rail-k", `${c.n_device_fields} device fields`), tagMix(c.device_tag_counts, { label: "device-tier fields" })]),
      clsChip(c.cls),
      c.edited ? h("span.chip.chip-dirty", "edited copy") : null),
    c.datasets.length
      ? h("ul.lib-ds-list", c.datasets.map((d) => h("li.lib-ds", tagChip(d.tag, { size: "sm", title: `dataset ${d.name}` }),
        h("span.lib-ds-name", d.name), h("span.lib-ds-n", h("span.lib-num", String(d.n_rows)), " rows"))))
      : h("p.lib-none", "no datasets: parameters only"));
}

// ---------------------------------------------------------------- card page
async function cardPage(name) {
  const seq = L.seq;
  const host = h("div.library.lib-card", h("div.lib-body", skeleton("block", { lines: 8, height: 520 })));
  L.el.append(host);
  let c;
  try { c = await L.ctx.api.getJSON(`/api/params/${encodeURIComponent(name)}`); } catch (e) {
    if (!L || seq !== L.seq) return;
    clear(host).append(h("header.ws-head.lib-head", h("div", h("a.lib-back", { href: "#/library" }, icon("chevronRight", { size: 16, cls: "flip" }), "Library"), h("h1.ws-title", name))),
      errorState({ title: `No parameter card ${name}`, message: e.message }));
    L.el.dataset.ready = "1";
    return;
  }
  if (!L || seq !== L.seq) return;
  L.card = c;
  clear(host);
  const actions = h("div.lib-actions",
    c.design ? h("a.btn.btn-outline", { href: `#/design/${encodeURIComponent(c.design)}`, title: `Open the design card ${c.design}, which derives from this card` }, icon("design", { size: 16 }), "Open as design") : null,
    c.editable ? h("button.btn#lib-edit", { type: "button", on: { click: () => setEditing(!L.editing) } }, icon("table", { size: 16 }), h("span.btn-txt", "Edit copy")) : null);
  host.append(
    h("header.ws-head.lib-head",
      h("div.lib-titles",
        h("a.lib-back", { href: "#/library" }, icon("chevronRight", { size: 16, cls: "flip" }), "Library"),
        h("div.lib-titleline", h("h1.ws-title", c.name), clsChip(c.cls), c.edited ? h("span.chip.chip-dirty", "edited copy") : null),
        c.paper ? h("p.ws-sub", h("strong", c.paper.label), ` · ${c.paper.cite}`) : null),
      actions),
    h("dl.lib-meta",
      metaRow("Device", c.device), metaRow("Role", c.role), metaRow("Source", c.source),
      c.cls?.reason ? metaRow("Class", c.cls.reason) : null,
      c.loader && !c.loader.startsWith("fsim_core") ? metaRow("Loader", c.loader) : null),
    c.meta?.note ? h("p.lib-note", c.meta.note) : null);

  const hasData = c.datasets.length > 0;
  const left = h("div.lib-col");
  const right = h("div.lib-col");
  host.append(h("div.lib-grid", { class: hasData ? "" : "is-single" }, left, hasData ? right : null));

  left.append(paramPanel(c));
  if (c.device_block?.length) left.append(devicePanel(c));
  if (c.placeholders?.length) left.append(h("p.lib-warn", icon("warn", { size: 14 }), `PLACEHOLDER entries: ${c.placeholders.join(", ")}; no gate can close on this card.`));

  const pending = [];
  for (const ds of c.datasets) pending.push(datasetCard(c, ds, right));
  await Promise.allSettled(pending);
  if (!L || seq !== L.seq) return;
  L.el.dataset.ready = "1";
}

function metaRow(k, v) {
  if (!v) return null;
  return h("div.lib-mrow", h("dt.engrave", k), h("dd", v));
}

function paramPanel(c) {
  const panel = h("section.panel.lib-params", { "aria-label": "Parameters" });
  const head = h("header.panel-head", h("h2.engrave", `Parameters · ${c.params.length}`), tagMix(c.tag_counts));
  const tableHost = h("div.lib-table-host");
  panel.append(head, tableHost);
  L.paramHost = tableHost;
  if (!c.params.length) {
    tableHost.append(h("p.lib-none.pad", "This card carries no params block; its numbers live in the device-tier configuration below."));
    return panel;
  }
  const nConf = c.params.filter((p) => p.conflict).length;
  if (nConf) {
    panel.insertBefore(h("p.lib-conflict-banner", { role: "note", dataset: { role: "conflict-banner" } }, icon("warn", { size: 16 }),
      h("strong", `Source conflict on ${nConf} parameter${nConf > 1 ? "s" : ""}: `),
      c.params.filter((p) => p.conflict).map((p) => p.name).join(", "), ". The card's value and the note below disagree; read the note before using the number."), tableHost);
  }
  renderParamTable(c);
  return panel;
}

function renderParamTable(c) {
  const host = L.paramHost;
  clear(host);
  const editing = L.editing;
  const rows = c.params.map((p) => {
    const cells = [
      h("td.lib-tagcell", tagChip(p.tag, { title: p.source })),
      h("th.lib-pname", { scope: "row" }, p.name),
    ];
    if (editing) {
      const v = h("input.fr-num.lib-in", { type: "number", step: "any", value: isNum(p.value) ? p.value : "", "aria-label": `${p.name} value`, dataset: { k: "value", p: p.name } });
      const lo = h("input.fr-num.lib-in", { type: "number", step: "any", value: p.range ? p.range[0] : "", "aria-label": `${p.name} range lo`, dataset: { k: "lo", p: p.name } });
      const hi = h("input.fr-num.lib-in", { type: "number", step: "any", value: p.range ? p.range[1] : "", "aria-label": `${p.name} range hi`, dataset: { k: "hi", p: p.name } });
      cells.push(h("td.lib-edit", h("span.lib-edit-row", v, h("span.lib-or", "or"), h("span.lib-edit-range", h("span.rb-txt", "["), lo, h("span.rb-txt", ","), hi, h("span.rb-txt", "]")))));
    } else {
      cells.push(h("td.lib-val", valueCell(p)));
    }
    // a `conflict:` note on the card (two sources disagree on this value) is shown verbatim and prominently
    cells.push(h("td.lib-unit", unitText(p.unit)), h("td.lib-src", p.source,
      p.conflict ? h("p.lib-conflict", { role: "note", dataset: { conflict: p.name } }, icon("warn", { size: 14 }), h("span.label-chip", "source conflict"), " ", p.conflict) : null));
    return h("tr", { dataset: { tagged: "1", param: p.name, conflict: p.conflict ? "1" : "0" } }, cells);
  });
  const table = h("div.table-wrap", h("table.data-table.lib-ptable",
    h("thead", h("tr", h("th", { scope: "col" }, "Tag"), h("th", { scope: "col" }, "Parameter"), h("th", { scope: "col" }, editing ? "Value, or [lo, hi]" : "Value / range"), h("th", { scope: "col" }, "Unit"), h("th", { scope: "col" }, "Source"))),
    h("tbody", rows)));
  host.append(table);
  if (editing) {
    const target = `${c.base}-edited`;
    host.append(h("div.lib-editbar",
      h("p.lib-edit-note", icon("info", { size: 14 }), ` Same schema as the legacy card editor: a value replaces the range; clear the value and give lo and hi to make it a range. Saves cards/${target}.yaml; shipped cards are never overwritten.`),
      h("div.lib-edit-btns",
        h("button.btn.btn-quiet", { type: "button", on: { click: () => setEditing(false) } }, "Cancel"),
        h("button.btn.btn-outline#lib-save", { type: "button", on: { click: (e) => save(c, target, e.currentTarget) } }, icon("download", { size: 16 }), `Save ${target}`))));
  }
}

function setEditing(on) {
  if (!L?.card) return;
  L.editing = on;
  const b = document.getElementById("lib-edit");
  if (b) { b.setAttribute("aria-pressed", String(on)); b.querySelector(".btn-txt").textContent = on ? "Editing copy" : "Edit copy"; }
  renderParamTable(L.card);
}

async function save(c, target, btn) {
  const rows = {};
  for (const inp of L.paramHost.querySelectorAll("input.lib-in")) {
    const p = inp.dataset.p;
    rows[p] = rows[p] || {};
    rows[p][inp.dataset.k] = inp.value === "" ? null : Number(inp.value);
  }
  btn.setAttribute("aria-busy", "true");
  try {
    const r = await L.ctx.api.putJSON(`/api/params/${encodeURIComponent(target)}`, { params: rows });
    toast(`Saved ${r.path}`, { kind: "ok" });
    if (L.ctx.args?.[0] === target || L.key === target) reload();
    else L.ctx.navigate(`#/library/${encodeURIComponent(target)}`);
  } catch (e) {
    toast(`Not saved: ${e.message}`, { kind: "error" });
  } finally {
    btn.removeAttribute("aria-busy");
  }
}

function devicePanel(c) {
  const rows = c.device_block.filter((r) => !r.hidden);
  const hidden = c.device_block.find((r) => r.hidden);
  const shown = rows.filter((r) => r.text != null || normTag(r.tag));
  const withheld = rows.filter((r) => r.text == null && !normTag(r.tag));
  return h("section.panel.lib-devblock", { "aria-label": "Device-tier configuration" },
    h("header.panel-head", h("h2.engrave", "Device-tier configuration"), h("span.panel-note", "tags as written in the card's own line comments")),
    h("div.table-wrap", h("table.data-table.lib-ptable",
      h("thead", h("tr", h("th", { scope: "col" }, "Tag"), h("th", { scope: "col" }, "Field"), h("th", { scope: "col" }, "Value / range"), h("th", { scope: "col" }, "Card comment"))),
      h("tbody", shown.map((r) => h("tr", { dataset: { tagged: r.text != null && !normTag(r.tag) ? "text" : "1" } },
        h("td.lib-tagcell", normTag(r.tag) ? tagChip(r.tag, { title: r.comment }) : h("span.lib-notag", "text")),
        h("th.lib-pname", { scope: "row" }, r.path),
        h("td.lib-val", valueCell(r)),
        h("td.lib-src", r.comment || "")))))),
    withheld.length ? h("p.lib-none.pad", `Withheld (no provenance tag in the card): ${withheld.map((r) => r.path).join(", ")}.`) : null,
    hidden ? h("p.lib-none.pad", `Range centres not shown (${hidden.hidden.length}): ${hidden.comment}.`) : null);
}

// ---------------------------------------------------------------- datasets
async function datasetCard(c, ds, host) {
  const ovStatus = c.overlays?.[ds.name] || { available: false };
  const st = { ov: null, ovErr: null };
  const chart = chartCard({
    title: ds.name.replace(/_/g, " "), height: 300,
    build: () => buildDataset(ds, st.ov),
    table: () => ({ caption: ds.excluded?.length ? `${ds.n_rows} rows (+${ds.excluded.length} excluded)` : `${ds.n_rows} rows`, columns: ds.columns.map((k) => ({ key: k, label: k })), rows: [...ds.rows, ...(ds.excluded || []).map((r) => ({ ...r, mechanism: r.mechanism || "excluded" }))] }),
    filename: () => `fsim-${c.name}-${ds.name}`,
  });
  L.charts.push(chart);
  chart.setTag(ds.tag);
  const ovLine = h("p.lib-ov", { dataset: { overlay: ovStatus.available ? "pending" : "none" } });
  const wrap = h("div.lib-ds-card", { dataset: { dataset: ds.name } }, chart.el,
    h("p.lib-ds-src", tagChip(ds.tag, { size: "sm", title: "dataset tag" }), h("span", ds.source), ds.source_note ? h("em", ` (${ds.source_note})`) : null),
    ovLine);
  host.append(wrap);
  chart.setTitle(ds.name.replace(/_/g, " "), `${ds.n_rows} rows · data [${ds.tag}]${ovStatus.available ? "" : " · no model overlay"}`);
  if (ovStatus.available) {
    try {
      st.ov = await L.ctx.api.getJSON(`/api/params/${encodeURIComponent(c.name)}/overlay/${encodeURIComponent(ds.name)}`);
      ovLine.dataset.overlay = "shown";
      const cls = st.ov.cls ? `${st.ov.cls.label}${st.ov.cls.code ? ` (${st.ov.cls.code})` : ""}` : "model";
      clear(ovLine).append(icon("chart", { size: 14 }), tagChip(st.ov.tag, { size: "sm", title: "overlay tag chain" }),
        h("span", ` Overlay: ${cls}. ${st.ov.label || ""} Calls: ${(st.ov.calls || []).join(", ")}.${st.ov.params_source ? ` Parameters: ${st.ov.params_source}.` : ""}`));
      chart.setTitle(ds.name.replace(/_/g, " "), `${ds.n_rows} rows · data [${ds.tag}] · ${cls} overlay [${normTag(st.ov.tag) || "?"}]`);
    } catch (e) {
      st.ovErr = e.message;
      ovLine.dataset.overlay = "error";
      clear(ovLine).append(icon("warn", { size: 14 }), h("span", ` Overlay did not resolve: ${e.message}`));
    }
  } else {
    clear(ovLine).append(icon("info", { size: 14 }), h("span", ` No model overlay: ${ovStatus.reason || "not allowed by the physics rules file"}.`));
  }
  if (!L) return;
  await chart.render();
}

function plotHint(ds) {
  if (ds.plot) return ds.plot;
  const numeric = ds.columns.filter((k) => ds.rows.every((r) => isNum(r[k])));
  return { x: numeric[0], y: numeric.slice(1, 4), x_label: numeric[0], y_label: "" };
}

function buildDataset(ds, ov) {
  const c = themeColors();
  const hint = plotHint(ds);
  const tag = normTag(ds.tag);
  const T = `[${tag}]`;
  const data = [];
  const annotations = [];
  const ink = c.ink1;
  const ring = { color: c.surface, width: 1.5 };
  const sentinel = hint.sentinel || {};
  const isSentinel = (r) => Object.entries(sentinel).some(([k, v]) => r[k] === v);
  const rows = ds.rows.filter((r) => !isSentinel(r));
  const sent = ds.rows.filter(isSentinel);
  const normal = rows.filter((r) => !r.mechanism && r.bound !== "upper");
  const ups = rows.filter((r) => r.bound === "upper");
  const mech = [...rows.filter((r) => r.mechanism), ...(ds.excluded || [])];

  if (hint.interval) {
    const [klo, khi] = hint.interval;
    const x = [], y = [];
    for (const r of rows) { x.push(r[hint.x], r[hint.x], null); y.push(r[klo], r[khi], null); }
    data.push({ x, y, type: "scatter", mode: "lines+markers", name: `${klo} – ${khi} interval ${T}`,
      line: { color: ink, width: 2 }, marker: { symbol: "line-ew-open", size: 12, color: ink, line: { color: ink, width: 2 } },
      hovertemplate: `${hint.x} %{x}: %{y:.3f} ${T}<extra></extra>` });
  } else {
    const groups = new Map();
    for (const r of normal) {
      const g = hint.group ? hint.group.map((k) => `${k} ${r[k]}`).join(", ") : "";
      if (!groups.has(g)) groups.set(g, []);
      groups.get(g).push(r);
    }
    let si = 0;
    for (const ycol of hint.y) {
      for (const [g, rs] of groups) {
        const errCol = hint.err_for?.[ycol] || (hint.err && ycol === hint.y[0] ? hint.err : null);
        const pts = rs.filter((r) => isNum(r[ycol]));
        if (!pts.length) continue;
        const nm = `${hint.y.length > 1 ? ycol : "measured"}${g ? `, ${g}` : ""} ${T}`;
        data.push({
          x: pts.map((r) => r[hint.x]), y: pts.map((r) => r[ycol]), type: "scatter",
          mode: hint.label ? "markers+text" : "markers", name: nm,
          text: hint.label ? pts.map((r) => r[hint.label]) : undefined, textposition: "middle right", textfont: { size: 12, color: c.ink2 }, cliponaxis: false,
          marker: { color: ink, size: 9, symbol: SYMBOLS[si % SYMBOLS.length], line: ring },
          error_y: errCol && pts.some((r) => isNum(r[errCol])) ? { type: "data", array: pts.map((r) => r[errCol] ?? 0), color: ink, thickness: 1.25, width: 5 } : undefined,
          hovertemplate: `${nm}: %{y:.4g}${errCol ? " ± %{error_y.array:.3g}" : ""}<extra></extra>`,
        });
        si++;
      }
    }
    if (ups.length) {
      const ycol = hint.y[0];
      data.push({ x: ups.map((r) => r[hint.x]), y: ups.map((r) => r[ycol]), type: "scatter", mode: "markers", name: `upper bound (≤) ${T}`,
        marker: { color: ink, size: 12, symbol: "triangle-down-open", line: { width: 2, color: ink } },
        hovertemplate: `upper bound ≤ %{y:.3f} ${T}<extra></extra>` });
    }
    if (mech.length) {
      const ycol = hint.y[0];
      data.push({ x: mech.map((r) => r[hint.x]), y: mech.map((r) => r[ycol]), type: "scatter", mode: "markers",
        name: `excluded: ${[...new Set(mech.map((r) => String(r.mechanism || "excluded").split(" (")[0]))].join("; ")} ${T}`,
        marker: { color: ink, size: 10, symbol: "circle-open", line: { width: 1.5, color: c.ink2 } },
        error_y: hint.err && mech.some((r) => isNum(r[hint.err])) ? { type: "data", array: mech.map((r) => r[hint.err] ?? 0), color: c.ink2, thickness: 1, width: 4 } : undefined,
        hovertemplate: `excluded %{y:.3f} ${T}<extra></extra>` });
    }
  }
  if (sent.length) {
    annotations.push({ xref: "paper", yref: "paper", x: 0, y: 0, xanchor: "left", yanchor: "bottom", showarrow: false,
      text: `not plotted: ${sent.map((r) => `${r[hint.label] || "row"} (${r.note || "sentinel value"})`).join("; ")}`, font: { size: 11, color: c.muted } });
  }

  // model overlay (only when the rules file allows it; backend-computed arrays)
  if (ov) {
    const otag = normTag(ov.tag);
    const cls = ov.cls?.code ? ` (${ov.cls.code})` : "";
    for (const b of ov.bands || []) {
      data.push({ x: [...b.x, ...[...b.x].reverse()], y: [...b.hi, ...[...b.lo].reverse()], type: "scatter", mode: "lines", fill: "toself",
        fillcolor: alpha(c.series[0], 0.12), line: { color: alpha(c.series[0], 0.6), width: 1 }, name: `${b.name}${cls} [${otag}]`, hoverinfo: "skip" });
    }
    (ov.series || []).forEach((s, i) => {
      data.push({ x: s.x, y: s.y, type: "scatter", mode: s.mode || "lines", name: `${s.name}${cls} [${otag}]`,
        line: { color: c.series[i % 8], width: tagWidth(otag), dash: tagDash(otag) },
        ...(s.err ? { error_y: { type: "data", array: s.err, visible: true, color: c.series[i % 8], thickness: 1.4, width: 4 } } : {}),
        marker: { color: c.series[i % 8], size: 8, symbol: "x-thin-open", line: { width: 2, color: c.series[i % 8] } },
        hovertemplate: `${s.name} %{y:.4g} [${otag}]<extra></extra>` });
    });
  }

  if (!data.length) return null;
  const layout = baseLayout({
    hovermode: "closest", margin: { l: 58, r: 18, t: 30, b: 46 },
    xaxis: { ...baseLayout().xaxis, showspikes: false, title: { text: hint.x_label || hint.x, standoff: 6 }, type: hint.categorical ? "category" : undefined },
    yaxis: { ...baseLayout().yaxis, title: { text: hint.y_label || "", standoff: 6 }, rangemode: hint.zero === false ? "normal" : "tozero" },
    annotations,
  });
  return { data, layout };
}
