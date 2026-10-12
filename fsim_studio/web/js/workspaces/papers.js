// Papers workspace (#/papers, #/papers/<id>): the paper-device library.
//   left   catalogue of paper plaques (cite + catalog ID, kind chips, read level by line form, class counts,
//          last ledger outcome n/N), filterable by material / structure / drive kind
//   main   device header (Reproduce key = the only beam-filled control, Open in Composer, job progress),
//          model-vs-paper table (one row per ledger item), series overlays (paperOverlay.js)
// No physics here: every number comes from GET /api/papers* and the reproduce job result. The only
// arithmetic is display: delta = model - paper, and the widest tag of the inputs for a model value.

import { h, clear } from "../ui/dom.js";
import { icon } from "../ui/icons.js";
import { tagChip, normTag, widest } from "../ui/tagChip.js";
import { skeleton, errorState, toast } from "../ui/controls.js";
import { fmt, isNum } from "../ui/format.js";
import { paperOverlay } from "../charts/paperOverlay.js";

const AXES = [
  { key: "material", label: "Material" },
  { key: "structure", label: "Structure" },
  { key: "drive", label: "Drive" },
];

const CLASS_INFO = {
  C: { label: "Calibration", tip: "Calibration: the model was fitted to this paper number (listed in the card's fits), so agreement is by construction. It is not a test." },
  M: { label: "Unfitted", tip: "Unfitted comparison: the model was built from other inputs and then compared with this paper number. The only class that counts as a test." },
  T: { label: "Transcription", tip: "Transcription: the observable passes an input taken from this paper straight through, so it can only agree. It is not a test." },
  P: { label: "Prediction", tip: "Prediction: the paper reports no value for this observable. The model number is stated, not tested." },
};

const READ_INFO = {
  "V-full": "Read in full: every number was checked against the paper text",
  "V-abs": "Read from the abstract only",
  "V-snip": "Read from a snippet only",
};

let P = null;

export function mount(el, ctx) {
  P = {
    el, ctx, list: [], id: null, paper: null, result: null, resultCached: false, running: false,
    charts: [], seq: 0, filters: { material: null, structure: null, drive: null }, dom: {},
    job: null, timer: null, t0: 0, eta: 5, listErr: null,
  };
  el.classList.add("ws-papers");
  ensureCss();
  el.dataset.ready = "0";
  buildShell();
  loadList().then(() => { if (P) update(P.ctx.args); });
}

export function unmount() {
  if (!P) return;
  stopJob();
  disposeCharts();
  P.el.classList.remove("ws-papers");
  delete P.el.dataset.ready;
  P = null;
}

export function update(args) {
  if (!P) return;
  const id = args?.[0] || null;
  if (id === P.id && P.paperLoaded) return;
  P.id = id;
  P.paperLoaded = false;
  P.seq++;
  stopJob();
  disposeCharts();
  P.paper = null;
  P.result = null;
  P.resultCached = false;
  syncSelected();
  if (!id) { renderPick(); P.el.dataset.ready = P.listDone ? "1" : "0"; return; }
  openPaper(id);
}

function ensureCss() {
  const href = new URL("../../css/papers.css", import.meta.url).href;
  const have = [...document.querySelectorAll('link[rel="stylesheet"]')].some((l) => l.href === href);
  if (!have) document.head.append(h("link", { rel: "stylesheet", href }));
}

function disposeCharts() {
  for (const c of P.charts) { try { c.dispose(); } catch { /* gone */ } }
  P.charts = [];
}

function stopJob() {
  clearInterval(P.timer);
  P.timer = null;
  try { P.job?.close(); } catch { /* closed */ }
  P.job = null;
  P.running = false;
}

// ---------------------------------------------------------------- shell + catalogue
function buildShell() {
  const d = P.dom;
  d.filters = h("div.pp-filters", { role: "group", "aria-label": "Filter the catalogue by material, structure or drive" });
  d.count = h("p.pp-count", { "aria-live": "polite" });
  d.plaques = h("div.pp-plaques", skeleton("plaque", { lines: 4 }), skeleton("plaque", { lines: 4 }));
  d.cat = h("aside.pp-cat", { "aria-label": "Paper catalogue" }, d.filters, d.count, d.plaques);
  d.main = h("section.pp-main", { "aria-label": "Paper device", dataset: { state: "idle" } });
  P.el.append(h("div.papers",
    h("header.ws-head.pp-head",
      h("div", h("h1.ws-title", "Papers"),
        h("p.ws-sub", "Each card builds a device from only what one paper reports, runs the composed engine, and records per reported number whether the model agrees. A mismatch is a recorded result, never tuned away.")),
      h("div.pp-head-links",
        h("a.btn.btn-quiet", { href: "#/library", title: "The older literature parameter cards (V-a tiers)" }, icon("library", { size: 16 }), "Literature cards (legacy tiers)"))),
    h("div.pp-grid", d.cat, d.main)));
}

async function loadList() {
  try {
    P.list = await P.ctx.api.getJSON("/api/papers");
    P.listErr = null;
  } catch (e) {
    if (!P) return;
    P.list = [];
    P.listErr = e;
  }
  if (!P) return;
  P.listDone = true;
  renderCatalogue();
}

function kindsOf(axis) {
  return [...new Set(P.list.map((p) => p.kinds?.[axis]).filter(Boolean))].sort();
}

function visible() {
  return P.list.filter((p) => AXES.every((a) => !P.filters[a.key] || p.kinds?.[a.key] === P.filters[a.key]));
}

function renderCatalogue() {
  const d = P.dom;
  clear(d.filters);
  if (P.list.length > 1) {
    for (const a of AXES) {
      const kinds = kindsOf(a.key);
      if (kinds.length < 2) continue;
      d.filters.append(h("div.pp-fgroup", { role: "group", "aria-label": a.label },
        h("span.pp-flabel.engrave", a.label),
        h("div.pp-fchips", kinds.map((k) => h("button.chip.pp-filter", {
          type: "button", "aria-pressed": String(P.filters[a.key] === k), dataset: { axis: a.key, kind: k },
          title: `${a.label}: ${k}. Click again to clear.`,
          on: { click: () => { P.filters[a.key] = P.filters[a.key] === k ? null : k; renderCatalogue(); } },
        }, k.replace(/_/g, " "))))));
    }
    if (AXES.some((a) => P.filters[a.key])) {
      d.filters.append(h("button.btn.btn-quiet.pp-clear", { type: "button", on: { click: () => { for (const a of AXES) P.filters[a.key] = null; renderCatalogue(); } } }, icon("close", { size: 14 }), "Clear filters"));
    }
  }
  clear(d.plaques);
  if (P.listErr) {
    d.count.textContent = "";
    d.plaques.append(errorState({ title: "Could not read the paper catalogue", message: P.listErr.message, actions: [{ label: "Retry", fn: () => { loadList(); } }] }));
    return;
  }
  if (!P.list.length) {
    d.count.textContent = "0 papers";
    d.plaques.append(h("div.cplaque.pp-empty", { dataset: { role: "empty" } },
      h("div.cp-top", h("div", h("h2.cp-title", "No paper cards yet"))),
      h("p.pp-empty-body", "A paper card is one YAML file, ", h("code", "cards/papers/<id>.yaml"), ", that builds a device from what a paper reports and lists, per number, what the model is expected to do."),
      h("p.pp-empty-body", "Copy the schema and the five rules for adding a paper from ", h("code", "cards/papers/README.md"), ", then reload this page. Files starting with an underscore are fixtures and are not listed.")));
    return;
  }
  const vis = visible();
  d.count.textContent = vis.length === P.list.length ? `${P.list.length} papers` : `${vis.length} of ${P.list.length} papers`;
  if (!vis.length) {
    d.plaques.append(h("div.cplaque.is-missing", h("p.pp-empty-body", "No paper matches this combination of filters.")));
    return;
  }
  d.plaques.append(...vis.map(plaque));
  syncSelected();
}

function shortCite(cite, id) {
  if (!cite) return id;
  const year = [...String(cite).matchAll(/\((\d{4})\)/g)].map((m) => m[1]).pop() || "";
  const first = String(cite).split(",")[0].replace(/\s+et al\.?/i, "").trim();
  const sur = first.split(/\s+/).pop() || id;
  const many = /et al/i.test(cite) || /^[^,]+,\s*[A-Z]\.\s/.test(cite);
  return `${sur}${many ? " et al." : ""}${year ? ` ${year}` : ""}`;
}

/** Short cite, with the card id appended when two cards share it (e.g. two lines of one paper). */
function titleOf(id, cite) {
  const t = shortCite(cite, id);
  const dup = P.list.filter((p) => shortCite(p.cite, p.id) === t).length > 1;
  return dup ? `${t} · ${id.replace(/^[a-z]+\d{4}_?/i, "") || id}` : t;
}

function readBadge(level) {
  return h("span.pp-read", { dataset: { level: level || "unknown" }, title: READ_INFO[level] || "read level not stated" },
    h("span.pp-read-k", "read"), h("span.pp-read-v", level || "?"));
}

function lastBadge(lo) {
  if (!lo) return h("span.vb.vb-sm.vb-info.pp-last", { dataset: { status: "none" }, title: "No reproduce result is cached for this card and engine yet" }, icon("minus", { size: 12 }), h("span.vb-word", "NOT RUN"));
  if (lo.n_items == null) return h("span.vb.vb-sm.vb-info.pp-last", { dataset: { status: "info" } }, icon("info", { size: 12 }), h("span.vb-word", "CACHED"));
  const ok = !!lo.all_match;
  return h("span.vb.vb-sm.pp-last", { class: ok ? "vb-pass" : "vb-fail", dataset: { status: ok ? "pass" : "fail" },
    title: `${lo.n_match} of ${lo.n_items} ledger items gave the outcome the card expected (this is ledger agreement, not model-equals-paper)` },
  icon(ok ? "check" : "cross", { size: 12 }), h("span.vb-word", ok ? "LEDGER MATCH" : "LEDGER MISMATCH"), h("span.pp-frac", `${lo.n_match}/${lo.n_items}`));
}

function classCounts(c) {
  return h("div.pp-classes", { role: "group", "aria-label": "Items per evidence class" },
    ["M", "C", "T", "P"].map((k) => h("span.pp-cc", { class: c?.[k] ? "" : "is-zero", dataset: { cls: k }, title: `${c?.[k] || 0} ${CLASS_INFO[k].label.toLowerCase()} item(s). ${CLASS_INFO[k].tip}` },
      h("span.pp-cc-k", k), h("span.pp-cc-n", String(c?.[k] || 0)))));
}

function plaque(p) {
  return h("a.cplaque.pp-plaque", { href: `#/papers/${encodeURIComponent(p.id)}`, dataset: { paper: p.id }, "aria-label": `${titleOf(p.id, p.cite)}: open` },
    h("div.cp-top",
      h("div.cp-titles", h("h2.cp-title", titleOf(p.id, p.cite)), h("p.pp-cite", { title: p.cite || "" }, p.cite || "")),
      h("span.cp-go", icon("chevronRight", { size: 18 }))),
    h("div.pp-kinds", AXES.map((a) => p.kinds?.[a.key] ? h("span.chip.pp-kind", { title: `${a.label}: ${p.kinds[a.key]}` }, p.kinds[a.key].replace(/_/g, " ")) : null)),
    h("div.pp-row", readBadge(p.read_level), classCounts(p.classes)),
    h("div.pp-foot",
      lastBadge(p.last_outcome),
      h("span.id-chip.is-static", { title: "Catalog ID" }, h("span.id-k", "ID"), h("span.id-v", `FS-paper-${p.id}`))));
}

function syncSelected() {
  if (!P?.dom.plaques) return;
  for (const a of P.dom.plaques.querySelectorAll(".pp-plaque")) {
    if (a.dataset.paper === P.id) a.setAttribute("aria-current", "true"); else a.removeAttribute("aria-current");
  }
}

// ---------------------------------------------------------------- main: pick / load
function renderPick() {
  const m = P.dom.main;
  m.dataset.state = "idle";
  clear(m);
  if (P.listErr || !P.list.length) { m.append(h("p.pp-hint", P.listErr ? "" : "Add a card to start.")); return; }
  m.append(h("div.cplaque.pp-pick",
    h("div.cp-top", h("div", h("h2.cp-title", "Choose a paper"))),
    h("p.pp-empty-body", "Open a plaque on the left. Reproduce runs the composed engine on the device built from that paper and fills the model column; nothing is run until you press it, except that an already-cached result for the same card bytes and engine is shown at once."),
    h("p.pp-empty-body", "Class badges say what an agreement is worth: ",
      ...["M", "C", "T", "P"].map((k) => [classBadge(k), " "]), "Only M items are tests.")));
}

async function openPaper(id) {
  const seq = P.seq;
  const m = P.dom.main;
  m.dataset.state = "loading";
  P.el.dataset.ready = "0";
  clear(m).append(skeleton("block", { lines: 5, height: 200 }), skeleton("block", { lines: 6, height: 260 }));
  let paper;
  try {
    paper = await P.ctx.api.getJSON(`/api/papers/${encodeURIComponent(id)}`);
  } catch (e) {
    if (!P || seq !== P.seq) return;
    m.dataset.state = "error";
    clear(m).append(errorState({
      title: e.status === 404 ? `No paper card ${id}` : `Paper card ${id} could not be loaded`,
      message: e.message + (e.body?.problems?.length ? `\n${e.body.problems.join("\n")}` : ""),
      actions: [{ label: "Back to the catalogue", fn: () => P.ctx.navigate("#/papers") }],
    }));
    P.el.dataset.ready = "1";
    return;
  }
  if (!P || seq !== P.seq) return;
  P.paper = paper;
  P.paperLoaded = true;
  renderPaper();
  // An exact cached result (same card bytes + engine hash) is shown at once; this POST returns it inline.
  if (paper.last_outcome?.cached) await reproduce({ auto: true });
  if (!P || seq !== P.seq) return;
  P.el.dataset.ready = "1";
}

// ---------------------------------------------------------------- main: paper page
function classBadge(k) {
  const i = CLASS_INFO[k];
  if (!i) return h("span.label-chip.pp-cls", { dataset: { cls: "none" } }, String(k || "?"));
  return h("span.label-chip.pp-cls", { dataset: { cls: k }, title: i.tip, "aria-label": `${k}: ${i.tip}` },
    h("b.pp-cls-k", k), h("span.pp-cls-l", i.label));
}

function outcomeBadge(outcome) {
  const o = String(outcome || "error").toLowerCase();
  const spec = o === "pass" ? ["vb-pass", "check", "PASS"]
    : o === "fail" ? ["vb-fail", "cross", "FAIL"]
      : o === "unavailable" ? ["vb-info", "info", "UNAVAILABLE"]
        : ["vb-fail", "warn", "ERROR"];
  return h("span.vb.vb-sm.pp-outcome", { class: spec[0], dataset: { outcome: o } }, icon(spec[1], { size: 12 }), h("span.vb-word", spec[2]));
}

/** Outcome differs from the card's expectation: a split badge (outcome | LEDGER MISMATCH). */
function mismatchBadge(outcome, expect) {
  const o = String(outcome || "error").toLowerCase();
  const word = o === "pass" ? "PASS" : o === "fail" ? "FAIL" : o === "unavailable" ? "UNAVAILABLE" : "ERROR";
  const ico = o === "pass" ? "check" : o === "unavailable" ? "info" : o === "fail" ? "cross" : "warn";
  const half = o === "pass" ? "is-pass" : o === "unavailable" ? "is-info" : "is-fail";
  return h("span.vb.vb-sm.vb-split.pp-outcome.pp-mismatch", {
    dataset: { outcome: o, mismatch: "1" }, role: "img",
    "aria-label": `${word}; ledger mismatch: the card expected ${String(expect).toUpperCase()}`,
    title: `The card expected ${String(expect).toUpperCase()}; the run gave ${word}. A flip in either direction needs a human review before the card changes.`,
  }, h("span.vb-half", { class: half }, icon(ico, { size: 12 }), word),
  h("span.vb-half.pp-mm", icon("warn", { size: 12 }), `LEDGER MISMATCH · expected ${String(expect).toUpperCase()}`));
}

function renderPaper() {
  const p = P.paper;
  const m = P.dom.main;
  m.dataset.state = "idle";
  clear(m);
  const meta = p.meta || {};
  const d = P.dom;
  d.runBtn = h("button.run-key.pp-reproduce#pp-reproduce", { type: "button", title: "Run the composed engine on this paper's device and fill the model column", on: { click: () => reproduce() } },
    icon("run", { size: 18 }), h("span", "REPRODUCE"));
  d.composerBtn = h("button.btn.btn-outline.pp-composer", { type: "button", title: "Open this device in the Composer at the central value of every ranged leaf", on: { click: (e) => openInComposer(e.currentTarget) } },
    icon("design", { size: 16 }), "Open in Composer");
  d.progress = h("div.progress", { role: "progressbar", "aria-label": "Reproduce progress", "aria-valuemin": "0", "aria-valuemax": "100", hidden: true }, h("span.progress-fill"));
  d.status = h("p.run-status.pp-status", { "aria-live": "polite", role: "status" });
  d.summary = h("div.pp-summary");
  d.err = h("div.pp-err");
  d.tableHost = h("div.pp-table-host");
  d.charts = h("div.pp-charts");

  const doi = meta.doi ? h("a.pp-doi", { href: `https://doi.org/${meta.doi}`, target: "_blank", rel: "noopener noreferrer", title: "Open the paper at doi.org" }, `doi:${meta.doi}`) : null;
  const arxiv = meta.arxiv ? h("a.pp-doi", { href: `https://arxiv.org/abs/${String(meta.arxiv).replace(/v\d+$/, "")}`, target: "_blank", rel: "noopener noreferrer" }, `arXiv:${meta.arxiv}`) : null;

  m.append(
    h("section.panel.pp-device", { "aria-label": "Device" },
      h("header.pp-dev-head",
        h("div.pp-dev-titles",
          h("div.pp-dev-line", h("h2.pp-dev-title", titleOf(p.id, meta.cite)), readBadge(meta.read_level),
            h("span.id-chip.is-static", { title: "Catalog ID" }, h("span.id-k", "ID"), h("span.id-v", `FS-paper-${p.id}`))),
          h("p.pp-dev-cite", meta.cite || ""),
          h("p.pp-dev-links", doi, doi && arxiv ? " · " : null, arxiv)),
        h("div.pp-actions",
          h("div.pp-status-col", d.progress, d.status),
          d.composerBtn, d.runBtn)),
      meta.summary ? h("p.pp-summary-text", meta.summary) : null,
      h("div.pp-kinds.pp-dev-kinds", AXES.map((a) => p.kinds?.[a.key] ? h("span.chip.pp-kind", { title: `${a.label}: ${p.kinds[a.key]}` }, h("span.pp-kind-k", a.label), p.kinds[a.key].replace(/_/g, " ")) : null),
        classCounts(p.classes))),
    d.err, d.summary,
    h("section.panel.pp-table-panel", { "aria-label": "Model versus paper" },
      h("header.panel-head", h("h2.engrave", `Model vs paper · ${p.reproduce.length} item${p.reproduce.length === 1 ? "" : "s"}`),
        h("span.panel-note", "ranges are [min, max] over the ranged inputs, never means")),
      d.tableHost),
    h("section.pp-charts-sec", { "aria-label": "Series overlays" }, d.charts));
  renderTable();
  renderCharts();
}

function reportedOf(id) {
  return P.paper.reported.find((r) => r.id === id) || null;
}
function seriesOf(id) {
  return P.paper.series.find((s) => s.id === id) || null;
}

function allTags() {
  return Object.values(P.paper.resolved || {}).map((r) => r?.tag);
}
/** Widest tag of the leaves an observable depends on (widest wins); falls back to the whole device. */
function modelTag(inputs) {
  const t = widest((inputs || []).map((k) => P.paper.resolved?.[k]?.tag));
  return t || widest(allTags());
}

function rows() {
  const ledger = P.paper.reproduce;
  const res = P.result?.items;
  if (!Array.isArray(res)) return ledger.map((e) => ({ ledger: e, id: e.item, cls: e.class, expect: e.expect || "pass", item: null }));
  const byId = new Map(res.map((i) => [i.id, i]));
  const out = ledger.map((e) => ({ ledger: e, id: e.item, cls: e.class, expect: e.expect || "pass", item: byId.get(e.item) || null }));
  for (const i of res) if (!out.some((r) => r.id === i.id)) out.push({ ledger: null, id: i.id, cls: i.cls, expect: i.expect, item: i });
  return out;
}

function renderTable() {
  const host = P.dom.tableHost;
  clear(host);
  const rs = rows();
  if (!rs.length) { host.append(h("p.pp-hint.pad", "This card has an empty reproduce ledger.")); return; }
  const th = (t, extra = {}) => h("th", { scope: "col", ...extra }, t);
  host.append(h("div.table-wrap", h("table.data-table.pp-table",
    h("thead", h("tr", th("Item"), th("Paper"), th("Convention"), th("Model"), th("Δ and tolerance"), th("Class"), th("Outcome"), th("Inputs from paper"))),
    h("tbody", rs.map(row)))));
}

function paperCell(r) {
  const rep = reportedOf(r.id);
  const ser = seriesOf(r.id);
  if (r.cls === "P" || (!rep && !ser)) {
    return h("td.pp-paper", h("span.pp-muted", "no paper value"));
  }
  if (ser) {
    const tag = normTag(ser.tag);
    return h("td.pp-paper", h("div.pp-val",
      h("span.pp-num", `${ser.points.length} points`), tagChip(tag, { title: ser.where, size: "sm" }),
      ser.digitized ? h("span.label-chip.pp-dig", { title: ser.method || "digitised from a figure" }, "digitised") : null),
    h("p.pp-where", { title: ser.where }, ser.where));
  }
  const tag = normTag(rep.tag);
  const err = isNum(rep.err) ? h("span.pp-err-v", { title: "paper-quoted uncertainty" }, ` ± ${fmt(rep.err)}`) : null;
  return h("td.pp-paper", h("div.pp-val", h("span.pp-num", fmt(rep.value)), err, tagChip(tag, { title: rep.where, size: "sm" })),
    h("p.pp-where", { title: rep.where }, rep.where));
}

function convCell(r) {
  const it = r.item;
  const rep = reportedOf(r.id);
  const ser = seriesOf(r.id);
  const conv = it?.convention || rep?.convention || ser?.y?.convention || r.ledger?.convention || "intrinsic";
  const obs = it?.observable || rep?.observable || ser?.y?.observable || r.ledger?.observable || "";
  return h("td.pp-conv", h("span.pp-obs", obs), h("span.label-chip.pp-conv-chip", { title: "paper and model are compared in the same convention only" }, conv));
}

function bracket(lo, hi) {
  return h("span.pp-range", { title: "model [min, max] over the ranged inputs: an envelope, never a mean" }, `[${fmt(lo)}, ${fmt(hi)}]`);
}

function modelCell(r) {
  const it = r.item;
  if (!it) return h("td.pp-model", h("span.pp-muted", "not run"));
  const ser = it.kind === "series" ? P.result?.model_series?.find((s) => s.id === it.id) : null;
  const tag = normTag(ser?.tag) || modelTag(it.inputs);
  if (it.outcome === "unavailable" || it.outcome === "error" || !isNum(it.model)) {
    return h("td.pp-model", h("span.pp-muted", "n/a"), it.note ? h("p.pp-note", it.note) : null);
  }
  if (!tag) return h("td.pp-model", h("span.pp-muted", "n/a"), h("p.pp-note", "no provenance tag; number withheld"));
  if (it.kind === "series") {
    return h("td.pp-model", h("div.pp-val", h("span.pp-k", "rms rel"), h("span.pp-num", fmt(it.model)), tagChip(tag, { title: "widest tag of the model curve", size: "sm" })),
      isNum(it.lo) && isNum(it.hi) ? h("p.pp-note", "curve ", bracket(it.lo, it.hi)) : null);
  }
  const ranged = isNum(it.lo) && isNum(it.hi) && Math.abs(it.hi - it.lo) > 1e-9 * Math.max(Math.abs(it.hi), Math.abs(it.lo), 1e-300);
  return h("td.pp-model", h("div.pp-val", h("span.pp-num", fmt(it.model)), tagChip(tag, { title: "widest tag of the inputs this observable depends on", size: "sm" })),
    ranged ? h("p.pp-note", bracket(it.lo, it.hi), tagChip(tag, { title: "same tag as the central value", size: "sm" })) : null);
}

function deltaCell(r) {
  const it = r.item;
  const tol = it?.tol_str || tolText(r.ledger?.tol);
  const tolEl = tol ? h("span.pp-tol", { title: it?.tol_reason || r.ledger?.tol_reason || "" }, `tol ${tol}`) : h("span.pp-muted", "no tolerance");
  if (!it || it.kind !== "scalar" || !isNum(it.model) || !isNum(it.paper)) {
    return h("td.pp-delta", it && it.kind === "series" && isNum(it.model) ? h("span.pp-num", `rms ${fmt(it.model)}`) : h("span.pp-muted", "n/a"), h("br"), tolEl);
  }
  const d = it.model - it.paper;
  return h("td.pp-delta", h("span.pp-num", { title: "model minus paper, central point" }, `${d > 0 ? "+" : d < 0 ? "−" : ""}${fmt(Math.abs(d))}`), h("br"), tolEl);
}

function tolText(t) {
  if (!t || typeof t !== "object") return "";
  if (t.metric) return `${t.metric}<=${t.max}`;
  return Object.entries(t).map(([k, v]) => `${k} ${v}`).join(", ");
}

function outcomeCell(r) {
  const it = r.item;
  if (!it) {
    return h("td.pp-outcome-cell", h("span.pp-expect", { title: "the outcome the card records before any run" }, `ledger expects ${String(r.expect).toUpperCase()}`),
      r.ledger?.reason ? h("p.pp-note", r.ledger.reason) : null);
  }
  const mismatch = it.outcome !== it.expect;
  const why = [["Ledger reason", it.reason], ["Tolerance reason", it.tol_reason], ["Engine note", it.note]].filter(([, v]) => v);
  return h("td.pp-outcome-cell", { dataset: { match: mismatch ? "0" : "1" } },
    mismatch ? mismatchBadge(it.outcome, it.expect) : outcomeBadge(it.outcome),
    !mismatch && it.expect !== "pass" ? h("p.pp-note", `as the card expected (${it.expect})`) : null,
    why.length ? h("details.pp-more", h("summary", "why"), why.map(([k, v]) => h("p.pp-why", h("span.pp-why-k", k), v))) : null);
}

function inputsCell(r) {
  const it = r.item;
  const set = new Set(P.paper.inputs_from_paper || []);
  if (!it) return h("td.pp-inputs", h("span.pp-muted", "run to list"));
  const ins = it.inputs || [];
  if (!ins.length) return h("td.pp-inputs", h("span.pp-muted", it.outcome === "unavailable" ? "no model run" : "none"));
  const nPaper = ins.filter((k) => set.has(k)).length;
  const near = ins.length > 0 && nPaper === ins.length;
  return h("td.pp-inputs", h("details.pp-in", { dataset: { fromPaper: String(nPaper), total: String(ins.length) } },
    h("summary", `${nPaper} of ${ins.length} from the paper`),
    near ? h("p.pp-note.pp-near", icon("warn", { size: 12 }), " every input comes from this paper: the item is near-transcription") : null,
    h("ul.pp-in-list", ins.map((k) => {
      const fp = set.has(k);
      return h("li", { dataset: { fromPaper: fp ? "1" : "0" } }, h("span.pp-leaf", k), tagChip(P.paper.resolved?.[k]?.tag, { size: "sm", title: fp ? "taken from the paper" : "not from the paper" }),
        fp ? h("span.label-chip.pp-fp", "paper") : h("span.pp-muted.pp-notfp", "not from paper"));
    }))));
}

function row(r) {
  const it = r.item;
  return h("tr", { dataset: { item: r.id, cls: r.cls, outcome: it?.outcome || "", state: it ? "run" : "ledger" } },
    h("th.pp-item", { scope: "row" }, h("span.pp-item-id", r.id), h("span.pp-item-kind", r.cls === "P" ? "prediction" : seriesOf(r.id) ? "series" : "scalar")),
    paperCell(r), convCell(r), modelCell(r), deltaCell(r),
    h("td.pp-class", classBadge(r.cls)),
    outcomeCell(r), inputsCell(r));
}

// ---------------------------------------------------------------- summary + charts
function renderSummary() {
  const host = P.dom.summary;
  clear(host);
  const res = P.result;
  if (!res || !Array.isArray(res.items)) return;
  const items = res.items;
  const n = items.length;
  const nMatch = items.filter((i) => i.match).length;
  const ok = n > 0 && nMatch === n;
  const Ms = items.filter((i) => i.cls === "M");
  const Mpass = Ms.filter((i) => i.outcome === "pass").length;
  const cnt = (o) => items.filter((i) => i.outcome === o).length;
  const rep = res.report || {};
  const status = ok ? "pass" : "fail";
  host.append(h("section.plaque.pp-plaque-sum", { class: `plaque-${status}`, dataset: { status }, "aria-label": "Reproduce summary" },
    h("div.plaque-top",
      h("div.plaque-word", icon(ok ? "check" : "cross", { size: 18 }), h("span", ok ? "LEDGER MATCHES" : "LEDGER MISMATCH")),
      h("div.plaque-denom", h("span.k", "items as expected"), h("span.v", `${nMatch} / ${n}`))),
    h("div.plaque-reasons",
      h("span.plaque-kv", h("span.k", "M tests"), h("span.v", `${Mpass} / ${Ms.length} pass`)),
      h("span.plaque-kv", h("span.k", "outcomes"), h("span.v", `${cnt("pass")} pass · ${cnt("fail")} fail · ${cnt("unavailable")} unavailable${cnt("error") ? ` · ${cnt("error")} error` : ""}`))),
    rep.valid === false && rep.invalid_reasons?.length
      ? h("p.plaque-reason", icon("warn", { size: 14 }), ` Run flagged invalid: ${rep.invalid_reasons.join("; ")}`) : null,
    h("p.plaque-source", `${P.resultCached ? "cached result for this card and engine" : "fresh run"} · card sha256 ${String(res.card_sha256 || "").slice(0, 12)}${rep.solver ? ` · solver ${rep.solver}` : ""}`),
    rep.notes?.length ? h("details.pp-more.pp-notes", h("summary", `Engine notes (${rep.notes.length})`), h("ul", rep.notes.map((t) => h("li", t)))) : null));
}

function renderCharts() {
  disposeCharts();
  const host = P.dom.charts;
  clear(host);
  const ser = P.paper.series || [];
  if (!ser.length) {
    host.append(h("div.pp-nochart", { dataset: { role: "no-series" } }, icon("info", { size: 14 }),
      h("span", " This card reports no series: its numbers are the scalar items above, so there is no curve to overlay.")));
    return;
  }
  const set = new Set(P.paper.inputs_from_paper || []);
  const grid = h("div.pp-chart-grid");
  host.append(h("h2.engrave.pp-sec-title", `Series overlays · ${ser.length}`), grid);
  for (const s of ser) {
    const model = P.result?.model_series?.find((m) => m.id === s.id) || null;
    const item = P.result?.items?.find((i) => i.id === s.id);
    const ov = paperOverlay({ series: s, model, inputs: item ? { paths: item.inputs || [], fromPaper: set } : null, filename: `fsim-paper-${P.paper.meta.id}-${s.id}` });
    P.charts.push(ov.chart);
    grid.append(ov.el);
    ov.render();
  }
}

// ---------------------------------------------------------------- reproduce job
function setRunState(state, text = "") {
  const d = P.dom;
  if (!d.runBtn) return;
  P.dom.main.dataset.state = state;
  d.runBtn.disabled = state === "running";
  d.runBtn.setAttribute("aria-busy", String(state === "running"));
  d.status.textContent = text;
  d.progress.hidden = state !== "running";
  if (state !== "running") d.progress.removeAttribute("aria-valuenow");
}

function tick() {
  if (!P?.dom.progress) return;
  const el = (performance.now() - P.t0) / 1000;
  const f = Math.min(0.95, el / Math.max(1, P.eta));
  const pct = Math.round(f * 100);
  P.dom.progress.firstChild.style.transform = `scaleX(${f})`;
  P.dom.progress.setAttribute("aria-valuenow", String(pct));
  P.dom.status.textContent = `Running the reproduce ledger · about ${Math.max(1, Math.round(P.eta))} s expected · ${el.toFixed(0)} s elapsed`;
}

async function reproduce({ auto = false } = {}) {
  if (!P?.paper || P.running) return;
  const seq = P.seq;
  const id = P.id;
  P.running = true;
  clear(P.dom.err);
  if (!auto) setRunState("running", "Submitting the reproduce job");
  let resp;
  try {
    resp = await P.ctx.api.postJSON(`/api/papers/${encodeURIComponent(id)}/reproduce`, {});
  } catch (e) {
    if (!P || seq !== P.seq) return;
    return failed(e);
  }
  if (!P || seq !== P.seq) return;
  try {
    let result = resp.result;
    if (!result) {
      P.eta = resp.eta_s || 5;
      P.t0 = performance.now();
      setRunState("running", "Queued");
      P.dom.progress.firstChild.style.transform = "scaleX(0)";
      clearInterval(P.timer);
      P.timer = setInterval(tick, 250);
      P.job = P.ctx.api.awaitJob(resp.job_id, { onState: (s) => { if (P && s && s !== "running") P.dom.status.textContent = `Job ${s}`; } });
      result = await P.job.promise;
    }
    if (!P || seq !== P.seq) return;
    clearInterval(P.timer);
    P.timer = null;
    P.job = null;
    P.running = false;
    P.result = result;
    P.resultCached = !!resp.cached;
    setRunState("done", `${result.items?.length ?? 0} items · ${resp.cached ? "cached result" : "fresh run"}${resp.run_id ? ` · ${resp.run_id}` : ""}`);
    P.dom.main.dataset.state = "done";
    renderSummary();
    renderTable();
    renderCharts();
    refreshListQuiet();
  } catch (e) {
    if (!P || seq !== P.seq) return;
    failed(e);
  }
}

function failed(e) {
  clearInterval(P.timer);
  P.timer = null;
  P.job = null;
  P.running = false;
  setRunState("error", "");
  P.dom.main.dataset.state = "error";
  const msg = e?.body?.error?.message || e?.message || "The reproduce job failed";
  clear(P.dom.err).append(errorState({ title: "Reproduce failed", message: msg, actions: [{ label: "Try again", fn: () => reproduce(), primary: true }] }));
}

async function refreshListQuiet() {
  try {
    const l = await P.ctx.api.getJSON("/api/papers");
    if (!P) return;
    P.list = l;
    renderCatalogue();
  } catch { /* the catalogue keeps its old badges */ }
}

// ---------------------------------------------------------------- Composer hand-off
async function openInComposer(btn) {
  if (!P?.paper) return;
  const id = P.id;
  btn.setAttribute("aria-busy", "true");
  try {
    const r = await P.ctx.api.getJSON(`/api/papers/${encodeURIComponent(id)}/design`);
    const seed = { design: r.design, paper: id, collapsed: r.collapsed || [], note: r.note || "" };
    P.ctx.state.set({ composeSeed: seed });
    try { sessionStorage.setItem("fsim-studio:compose-seed", JSON.stringify(seed)); } catch { /* private window: state.js carries it */ }
    P.ctx.navigate("#/compose");
  } catch (e) {
    toast(`Could not open in the Composer: ${e.message}`, { kind: "error" });
  } finally {
    btn.removeAttribute("aria-busy");
  }
}
