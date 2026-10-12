// Composer workspace: pick ANY material x structure x drive x setup from the registry and run it.
// Board (left 60%): viz3d viewport, the beam chain of mounts (Drive -> Dot -> Structure -> Setup),
// the temperature plate, parameter sheet, run bar. Shelf (right 40%): validity plaque, readouts,
// chart cards. No physics and no per-kind branching here: every kind, parameter, band, tag and
// status comes from GET /api/compose/registry; every number from /api/compose/*.
import { h, clear, svgEl, debounce, reducedMotion } from "../ui/dom.js";
import { icon } from "../ui/icons.js";
import { tagChip, normTag, widest, tagDash, tagWidth } from "../ui/tagChip.js";
import { createReadout } from "../ui/readout.js";
import { segmented, toast, skeleton, emptyState, errorState, popover, closePopover } from "../ui/controls.js";
import { fmt, fmtSeconds, isNum, unitText } from "../ui/format.js";
import { chartCard } from "../charts/chartCard.js";
import { baseLayout, colors as themeColors } from "../charts/theme.js";
import { renderParamForm, readLeaf } from "../ui/paramForm.js";
import { openChipPicker } from "../ui/chipPicker.js";

// Chain order and display names. The registry key is only a lookup into the payload.
const AXES = [
  { key: "drive", label: "Drive", reg: "drives" },
  { key: "material", label: "Dot", reg: "materials" },
  { key: "structure", label: "Structure", reg: "structures" },
  { key: "setup", label: "Setup", reg: "setups" },
];
const SCHEMA = "fsim-composed/1";
// Display limits of the temperature slider (the spec accepts the wider T_K band; the number box is free).
const T_SLIDER = [1, 300];
const TSWEEP_MAX = 31;
const SERIES_MAX = 200;
const T_SOURCE = "operating temperature set in the Composer";
// Readouts on the shelf: observable id, label, unit, which result field carries the number.
const READOUTS = [
  { key: "g2_0", label: "g²(0)", unit: "", curve: "g2" },
  { key: "hom_I", label: "I_HOM", unit: "" },
  { key: "B_first_lens", label: "B, first lens", unit: "" },
  { key: "rate_detected_Hz", label: "Detected rate", unit: "Hz" },
  { key: "F_P_eff", label: "F_P · L", unit: "" },
  { key: "tau_ps", label: "τ", unit: "ps" },
];
const DEFAULT_OBS = ["g2_0", "hom_I", "B_first_lens"];
const DEFAULT_EXAMPLE = "example-pillar-resonant";
const clone = (x) => JSON.parse(JSON.stringify(x));

let S = null;

// ---------------------------------------------------------------- lifecycle
export function mount(el, ctx) {
  S = {
    ctx, el, reg: null, spec: null, entries: {}, val: null, valSeq: 0, run: null, series: null,
    mode: "point", tGrid: { lo: 4, hi: 100, n: 13 }, ser: { path: null, lo: 0, hi: 1, n: 21, obs: [] },
    openAxis: null, forms: {}, charts: {}, readouts: {}, dom: {}, job: null, runState: "idle", timer: null,
    scene: null, sceneSeq: 0, card: null, edits: 0, lastError: null, disposers: [], tt0: 0,
  };
  el.classList.add("ws-compose");
  S.disposers.push(ctx.state.on("theme", (m) => S?.scene?.setTheme?.(m)));
  update(ctx.args);
}

export function unmount() {
  if (!S) return;
  closePopover();
  S.job?.cancel?.();
  clearInterval(S.timer);
  for (const d of S.disposers) d();
  for (const c of Object.values(S.charts)) c.dispose();
  try { S.scene?.dispose?.(); } catch { /* gone */ }
  S.ro?.disconnect();
  S.ctx.state.set({ run: null });
  S.el.classList.remove("ws-compose");
  S = null;
}

export async function update(args) {
  if (!S) return;
  const card = args?.[0] || null;
  if (S.dom.board && card === S.card) return;
  S.card = card;
  S.job?.cancel?.();
  renderSkeleton();
  try {
    S.reg = await S.ctx.api.composeRegistry();
    indexRegistry();
    if (card) {
      S.seed = null;
      const info = await S.ctx.api.getJSON(`/api/cards/${encodeURIComponent(card)}`);
      const composed = info?.design?.composed;
      if (!composed) throw new Error(`card ${card} carries no composed spec`);
      S.spec = normaliseSpec(composed);
      for (const a of AXES) if (!S.entries[a.key].get(S.spec[a.key].kind)) throw new Error(`card ${card}: unknown ${a.label} kind "${S.spec[a.key].kind}"`);
    } else if ((S.seed = takeSeed())) {
      S.spec = normaliseSpec(S.seed.design.composed);
      for (const a of AXES) if (!S.entries[a.key].get(S.spec[a.key].kind)) throw new Error(`paper design: unknown ${a.label} kind "${S.spec[a.key].kind}"`);
    } else {
      S.seed = null;
      S.spec = await startSpec();
    }
    if (!S || S.card !== card) return;
    S.run = null; S.series = null; S.edits = 0; S.lastError = null; S.openAxis = null;
    S.ctx.state.set({ run: null });
    S.ser.path = null;
    renderAll();
    await revalidate();
    loadScene();
  } catch (e) {
    if (!S) return;
    clear(S.el).append(h("div.ws-pad", errorState({
      title: card ? `Could not open ${card}` : "Could not open the Composer", message: e.message,
      actions: [{ label: "Retry", fn: () => { S.card = "\u0000"; update(args); } }],
    })));
  }
}

/** Papers -> "Open in Composer": state.composeSeed first, sessionStorage as fallback; taken once. */
function takeSeed() {
  let seed = S.ctx.state.get("composeSeed");
  if (!seed) { try { seed = JSON.parse(sessionStorage.getItem("fsim-studio:compose-seed") || "null"); } catch { seed = null; } }
  S.ctx.state.set({ composeSeed: null });
  try { sessionStorage.removeItem("fsim-studio:compose-seed"); } catch { /* private window */ }
  return seed?.design?.composed ? seed : null;
}

// ---------------------------------------------------------------- registry + spec
function indexRegistry() {
  for (const a of AXES) {
    const m = new Map();
    for (const e of S.reg[a.reg] || []) m.set(e.kind, e);
    S.entries[a.key] = m;
  }
}
const entry = (axis) => S.entries[axis].get(S.spec[axis].kind);
const axisOf = (key) => AXES.find((a) => a.key === key);

function tLeaf(v, tag = "A", source = T_SOURCE) { return { value: v, unit: "K", tag, source }; }

function normaliseSpec(raw) {
  const spec = { schema: SCHEMA, solver: raw.solver || "auto" };
  const t = raw.T_K;
  spec.T_K = t !== null && typeof t === "object" ? clone(t) : tLeaf(Number(t));
  for (const a of AXES) spec[a.key] = { kind: raw[a.key].kind, params: clone(raw[a.key].params || {}) };
  return spec;
}

function specWith(axis, kind) {
  const s = clone(S.spec);
  s[axis] = { kind, params: {} };
  return s;
}

/** Opening with no card or seed: the default example; the first-valid search when it is missing. */
async function startSpec() {
  try {
    S.examples = await S.ctx.api.getJSON("/api/compose/examples");
  } catch { S.examples = []; }
  const ex = (S.examples || []).find((e) => e.name === DEFAULT_EXAMPLE);
  if (ex) {
    const spec = normaliseSpec(ex.design.composed);
    if (AXES.every((a) => S.entries[a.key].get(spec[a.key].kind))) { S.exampleName = ex.name; return spec; }
  }
  S.exampleName = null;
  return defaultSpec();
}

/** Replace the working spec by a design from the examples row or a tier chip, then redraw. */
async function openDesign(design, name) {
  const spec = normaliseSpec(design.composed);
  for (const a of AXES) if (!S.entries[a.key].get(spec[a.key].kind)) { toast(`${name}: unknown ${a.label} kind "${spec[a.key].kind}"`, { kind: "error" }); return; }
  S.job?.cancel?.();
  S.spec = spec; S.exampleName = name; S.seed = null;
  S.run = null; S.series = null; S.edits = 0; S.lastError = null; S.openAxis = null; S.val = null;
  S.ctx.state.set({ run: null });
  S.ser.path = null;
  renderAll();
  await revalidate();
  loadScene();
}

/** First ready entry per axis, then the first combination the validator accepts (greedy per axis). */
async function defaultSpec() {
  const ready = (a) => (S.reg[a.reg] || []).filter((e) => e.status === "ready");
  const spec = { schema: SCHEMA, solver: "auto", T_K: tLeaf(4.2) };
  for (const a of AXES) spec[a.key] = { kind: ready(a)[0].kind, params: {} };
  let res = await S.ctx.api.composeValidate(spec);
  for (let pass = 0; pass < 2 && !res.ok; pass++) {
    for (const a of AXES) {
      for (const e of ready(a)) {
        const t = clone(spec); t[a.key] = { kind: e.kind, params: {} };
        const r = await S.ctx.api.composeValidate(t);
        if (r.ok) { spec[a.key] = t[a.key]; res = r; break; }
      }
      if (res.ok) break;
    }
  }
  return spec;
}

// ---------------------------------------------------------------- validation
const scheduleValidate = debounce(() => { revalidate(); }, 250);

async function revalidate() {
  if (!S) return;
  const seq = ++S.valSeq;
  try {
    const v = await S.ctx.api.composeValidate(S.spec);
    if (!S || seq !== S.valSeq) return;
    S.val = v;
  } catch (e) {
    if (!S || seq !== S.valSeq) return;
    S.val = { ok: false, issues: [{ level: "error", path: "request", message: e.message }], resolved: {}, solver: null, eta_s: null };
  }
  renderMounts();
  renderPlaque();
  syncRunbar();
  drawBeam();
  sceneSoon();
}

const sceneSoon = debounce(() => loadScene(), 500);

function issuesFor(axis) {
  return (S.val?.issues || []).filter((i) => String(i.path || "").startsWith(`${axis}.`) || i.path === `${axis}.kind` || i.path === axis);
}

// ---------------------------------------------------------------- skeleton + layout
function renderSkeleton() {
  clear(S.el).append(h("div.composer.is-loading", { "aria-busy": "true" },
    h("section.board", h("div.board-head", skeleton("title", { lines: 2 })),
      h("div.board-table", h("div.sk-viewport", skeleton("block", { lines: 0 })),
        h("div.sk-chain", Array.from({ length: 4 }, () => skeleton("mount", { lines: 3 }))))),
    h("section.shelf", skeleton("plaque", { lines: 2 }), skeleton("readouts", { lines: 5 }), skeleton("chart", { lines: 0, height: 260 }))));
}

function renderAll() {
  for (const c of Object.values(S.charts)) c.dispose();
  S.charts = {};
  try { S.scene?.dispose?.(); } catch { /* gone */ }
  S.scene = null;
  S.sceneKey = null;
  S.ro?.disconnect();
  S.dom = {};
  clear(S.el);
  const board = renderBoard();
  const shelf = renderShelf();
  S.el.append(h("div.composer", board, shelf));
  S.dom.board = board;
  S.ro = new ResizeObserver(() => drawBeam());
  S.ro.observe(S.dom.chain);
  renderMounts();
  renderPlaque();
  renderResults();
  syncRunbar();
  requestAnimationFrame(drawBeam);
}

// ================================================================ BOARD
function renderBoard() {
  const tiers = S.reg.tiers || [];
  const head = h("header.board-head",
    h("div.bh-titles",
      h("h1.board-title", S.card ? S.card.replace(/-design$/, "") : "Composer"),
      S.seed ? h("p.seed-note.compose-seed", { dataset: { paper: S.seed.paper || "" } }, icon("info", { size: 14 }), `from paper ${S.seed.paper || ""}: ${S.seed.note || "ranged inputs collapsed to their central values"}${S.seed.collapsed?.length ? ` (${S.seed.collapsed.length} collapsed)` : ""}`) : null,
      h("p.board-desc", "Pick any dot, structure, drive and setup from the registry. Every value keeps its provenance tag; edits are yours and default to [A].")),
    h("div.bh-meta",
      S.dom.dirty = h("span.chip.chip-dirty", { hidden: true, title: "Edits are local until you save them as a card" }, "edited"),
      h("div.bh-keys",
        S.dom.saveBtn = h("button.btn.btn-outline.bh-save", { type: "button", "aria-haspopup": "dialog", on: { click: (e) => saveAsPopover(e.currentTarget) } }, icon("download", { size: 16 }), "Save as composed card"))));

  const tierStrip = tiers.length ? h("div.tier-strip", { role: "group", "aria-label": "Tier presets" },
    h("span.engrave", "Tiers"),
    tiers.map((t) => {
      const design = t.handoff === "design" && (t.cards || []).length;
      const here = t.handoff === "compose";
      const on = design || here;
      const b = h("button.chip.tier-chip", {
        type: "button", dataset: { tier: t.id, handoff: t.handoff || "" }, disabled: !on,
        title: design ? `${t.summary} Opens ${t.cards[0]} in the Designer.`
          : here ? `${t.summary} Opens here as a composed spec.`
          : `${t.summary} This tier has no spec mapping, so it cannot open here.`,
        on: { click: async () => {
          if (design) { S.ctx.navigate(`#/design/${t.cards[0]}`); return; }
          if (!here) return;
          try {
            const d = await S.ctx.api.getJSON(`/api/compose/tier/${encodeURIComponent(t.id)}`);
            await openDesign(d.design, d.design.name);
          } catch (e) { toast(`Could not open ${t.label}: ${e.message}`, { kind: "error" }); }
        } },
      }, t.label, on ? icon("chevronRight", { size: 12 }) : null);
      return b;
    })) : null;

  const exampleStrip = (S.examples || []).length ? h("div.tier-strip.example-strip", { role: "group", "aria-label": "Example devices" },
    h("span.engrave", "Examples"),
    S.examples.map((ex) => h("button.chip.tier-chip.example-chip", {
      type: "button", dataset: { example: ex.name }, "aria-pressed": String(S.exampleName === ex.name),
      title: ex.summary, on: { click: () => openDesign(ex.design, ex.name) },
    }, ex.label))) : null;

  S.dom.viewport = h("div.viewport-slot", { "aria-label": "Device model", role: "region" });
  const viewportMount = h("div.viewport-mount",
    h("div.vm-rail", h("span.engrave", "Device model"), S.dom.vmNote = h("span.vm-note", "")),
    S.dom.viewport, h("span.post", { "aria-hidden": "true" }));

  S.dom.beamSvg = svgEl("svg", { class: "beam-svg", "aria-hidden": "true" });
  S.dom.mountHost = {};
  const cells = AXES.map((a) => (S.dom.mountHost[a.key] = h("div.mount-cell", { role: "listitem", dataset: { axis: a.key } })));
  S.dom.chain = h("div.chain", { role: "list", "aria-label": "Chain: drive to setup" }, S.dom.beamSvg, cells);
  S.dom.beamCap = h("p.beam-cap");

  S.dom.temp = renderTemp();
  S.dom.sheet = h("div.sheet-host");
  S.dom.runbar = renderRunbar();

  const table = h("div.board-table", viewportMount,
    h("div.chain-wrap", h("div.chain-scroll", S.dom.chain), S.dom.beamCap),
    S.dom.temp, S.dom.sheet);
  return h("section.board", { "aria-label": "Breadboard" }, head, tierStrip, exampleStrip, table, S.dom.runbar);
}

// ---------------------------------------------------------------- mounts
function leafText(p, raw) {
  const l = readLeaf(p, raw);
  const u = unitText(p.unit === "enum" ? "" : p.unit);
  if (l.kind === "range") return { text: `[${fmt(l.range[0])}, ${fmt(l.range[1])}]`, unit: u, tag: l.tag, ranged: true };
  const v = l.value;
  return { text: typeof v === "number" ? fmt(v) : typeof v === "boolean" ? (v ? "on" : "off") : String(v), unit: u, tag: l.tag, ranged: false };
}

function renderMounts() {
  for (const a of AXES) {
    const host = S.dom.mountHost[a.key];
    const e = entry(a.key);
    const params = S.spec[a.key].params;
    const main = e.params.filter((p) => p.group !== "advanced").slice(0, 3);
    const tags = e.params.map((p) => readLeaf(p, params[p.name]).tag);
    const wide = widest(tags);
    const nEdit = e.params.filter((p) => params[p.name] !== undefined).length;
    const bad = issuesFor(a.key).filter((i) => i.level === "error");
    const open = S.openAxis === a.key;
    const vals = h("div.mount-vals", main.map((p) => {
      const t = leafText(p, params[p.name]);
      return h("div.mv", { title: `${p.name}: ${p.help || ""}` },
        h("span.mv-k", p.name),
        h("span.mv-line", h("span.mv-v", t.text), h("span.mv-u", t.unit), tagChip(t.tag, { size: "sm", title: p.source })));
    }));
    const mount = h("div.mount", { class: `${open ? "is-open" : ""} ${bad.length ? "has-issue" : ""}`, dataset: { axis: a.key, rim: wide || "" }, role: "group", "aria-label": `${a.label} mount: ${e.label}` },
      h("span.mount-head", h("span.engrave", a.label), wide ? tagChip(wide, { size: "sm", title: "widest tag of this mount's inputs" }) : null),
      h("span.mount-aperture", { "aria-hidden": "true" }),
      h("button.kind-chip", {
        type: "button", dataset: { axis: a.key, kind: e.kind }, "aria-haspopup": "dialog", title: e.summary,
        on: { click: (ev) => openPicker(a, ev.currentTarget) },
      }, h("span.kc-label", e.label), icon("chevronDown", { size: 14 })),
      h("button.mount-body", {
        type: "button", dataset: { axis: a.key }, "aria-expanded": String(open),
        "aria-label": `Open the ${a.label} parameter sheet`,
        on: { click: () => toggleSheet(a.key) },
      }, vals,
      bad.length ? h("span.mount-issue", icon("warn", { size: 13 }), `${bad.length} issue${bad.length === 1 ? "" : "s"}`) : null,
      h("span.mount-foot", `${e.params.length} field${e.params.length === 1 ? "" : "s"}${nEdit ? `, ${nEdit} edited` : ""}`)));
    host.replaceChildren(mount, h("span.post", { "aria-hidden": "true" }));
  }
}

function toggleSheet(axis) {
  S.openAxis = S.openAxis === axis ? null : axis;
  renderMounts();
  renderSheet();
  requestAnimationFrame(drawBeam);
}

function renderSheet() {
  const host = S.dom.sheet;
  clear(host);
  if (!S.openAxis) { host.classList.remove("is-open"); return; }
  const a = axisOf(S.openAxis);
  const e = entry(a.key);
  const form = renderParamForm({
    axis: a.key, entry: e,
    getParams: () => S.spec[a.key].params,
    onEdit: (name, leaf) => editParam(a.key, name, leaf),
    onInspect: (info, o) => S.ctx.inspector.field(info, o),
  });
  S.forms[a.key] = form;
  host.classList.add("is-open");
  host.append(h("section.sheet", { "aria-label": `${a.label} parameters`, dataset: { axis: a.key } },
    h("header.sheet-head",
      h("div", h("h2.engrave.sheet-title", `${a.label}: ${e.label}`),
        h("p.sheet-sub", `${e.params.length} fields from the registry. Edits apply to the next run; bands warn, they never clamp.`)),
      h("button.icon-btn", { type: "button", "aria-label": `Close ${a.label} sheet`, on: { click: () => toggleSheet(a.key) } }, icon("close", { size: 18 }))),
    e.refs?.length ? h("details.sheet-refs", h("summary", `${e.refs.length} reference${e.refs.length === 1 ? "" : "s"}`), h("ul", e.refs.map((r) => h("li", r)))) : null,
    form.el));
  if (!reducedMotion()) host.querySelector(".sheet").classList.add("mount-in");
  requestAnimationFrame(() => host.scrollIntoView({ block: "nearest", behavior: reducedMotion() ? "auto" : "smooth" }));
}

function markEdited() {
  S.edits += 1;
  if (S.dom.dirty) S.dom.dirty.hidden = false;
}

function editParam(axis, name, leaf) {
  const params = S.spec[axis].params;
  if (leaf === null) delete params[name]; else params[name] = leaf;
  markEdited();
  renderMounts();
  scheduleValidate();
  refreshSeriesPaths();
}

// ---------------------------------------------------------------- picker
function openPicker(a, anchor) {
  const entries = S.reg[a.reg] || [];
  const base = new Set((S.val?.issues || []).filter((i) => i.level === "error").map((i) => `${i.path}|${i.message}`));
  openChipPicker(anchor, {
    axisLabel: a.label, entries, current: S.spec[a.key].kind,
    check: async (kind) => {
      const pair = (S.reg.compat || []).find(([ax, k, ax2, k2]) => (ax === a.key && k === kind && S.spec[ax2]?.kind === k2) || (ax2 === a.key && k2 === kind && S.spec[ax]?.kind === k));
      if (pair) return { ok: false, reason: pair[4] };
      try {
        const r = await S.ctx.api.composeValidate(specWith(a.key, kind));
        const fresh = (r.issues || []).filter((i) => i.level === "error" && !base.has(`${i.path}|${i.message}`));
        return fresh.length ? { ok: false, reason: fresh.map((i) => i.message).join("; ") } : { ok: true };
      } catch (e) { return { ok: false, reason: e.message }; }
    },
    onPick: (kind) => {
      if (kind === S.spec[a.key].kind) return;
      S.spec[a.key] = { kind, params: {} };
      markEdited();
      S.run = null; S.series = null;
      S.ser.path = null;
      renderMounts();
      if (S.openAxis === a.key) renderSheet();
      renderResults();
      syncRunbar();
      revalidate();
      requestAnimationFrame(() => S.dom.mountHost[a.key].querySelector(".kind-chip")?.focus());
    },
  });
}

// ---------------------------------------------------------------- beam
function beamAlpha() {
  const rate = S.run?.scalars?.rate_detected_Hz;
  if (!(isNum(rate) && rate > 0)) return { a: 0.45, text: S.run ? "Beam at rest level: no detected rate on this path" : "Beam at rest level until the first run" };
  const k = Math.max(0, Math.min(1, (Math.log10(rate) - 3) / 5));
  return { a: 0.45 + 0.55 * k, text: `Beam intensity follows log₁₀ detected rate: ${fmt(rate)} Hz` };
}

function drawBeam() {
  if (!S?.dom?.chain) return;
  const svg = S.dom.beamSvg;
  const box = S.dom.chain.getBoundingClientRect();
  if (!box.width) return;
  svg.setAttribute("viewBox", `0 0 ${box.width} ${box.height}`);
  svg.setAttribute("width", box.width);
  svg.setAttribute("height", box.height);
  clear(svg);
  const at = (key) => {
    const m = S.dom.mountHost[key]?.querySelector(".mount");
    if (!m) return null;
    const r = m.getBoundingClientRect();
    const ap = m.querySelector(".mount-aperture")?.getBoundingClientRect();
    return { x: r.left - box.left + r.width / 2, y: ap ? ap.top - box.top + ap.height / 2 : r.top - box.top + r.height / 2 };
  };
  const drive = at("drive"), dot = at("material"), last = at("setup");
  if (!drive || !dot || !last) return;
  const y = dot.y;
  svg.append(svgEl("path", { d: `M ${drive.x} ${y} L ${dot.x} ${y}`, class: "lead" }));
  const { a, text } = beamAlpha();
  const running = S.runState === "running";
  const d = `M ${dot.x} ${y} L ${last.x} ${y}`;
  svg.append(svgEl("path", { d, class: "beam-halo", style: `opacity:${(a * 0.4).toFixed(3)}` }));
  svg.append(svgEl("path", { d, class: `beam${running && !reducedMotion() ? " is-running" : ""}`, style: `opacity:${a.toFixed(3)}`, "data-alpha": a.toFixed(3) }));
  svg.append(svgEl("circle", { cx: dot.x, cy: y, r: 7, class: "beam-src", style: `opacity:${a.toFixed(3)}` }));
  svg.append(svgEl("path", { d: `M ${last.x} ${y - 9} L ${last.x} ${y + 9}`, class: "beam-stop" }));
  S.dom.beamCap.textContent = text;
}

// ---------------------------------------------------------------- temperature plate
function renderTemp() {
  const leaf = () => S.spec.T_K;
  const tNow = () => Number(leaf().value ?? (Array.isArray(leaf().range) ? leaf().range[0] : 4.2));
  const valEl = h("input.fr-num.temp-num#temp-num", {
    type: "number", step: "any", value: String(Number(tNow().toPrecision(5))), "aria-label": "Temperature T_K in kelvin",
    on: { change: (e) => { const v = Number(e.target.value); if (Number.isFinite(v) && v > 0) setT(v); } },
  });
  const slider = h("input.hs-slider#temp-slider", {
    type: "range", min: T_SLIDER[0], max: T_SLIDER[1], step: 0.5, value: Math.min(T_SLIDER[1], Math.max(T_SLIDER[0], tNow())),
    "aria-label": "Temperature T_K, kelvin. Dragging revalidates",
    on: { input: (e) => { const v = Number(e.target.value); valEl.value = String(v); setT(v, true); } },
  });
  const sel = h("div.pf-tagsel", { role: "group", "aria-label": "Temperature provenance tag" });
  const syncSel = () => {
    clear(sel);
    for (const t of ["V", "DR", "E", "A"]) {
      const b = tagChip(t, { size: "sm", title: `set the tag of T_K to [${t}]`, onClick: () => { S.spec.T_K = { ...leaf(), tag: t }; markEdited(); syncSel(); scheduleValidate(); } });
      b.classList.toggle("is-on", t === normTag(leaf().tag));
      b.setAttribute("aria-pressed", String(t === normTag(leaf().tag)));
      sel.append(b);
    }
  };
  function setT(v) {
    S.spec.T_K = { ...leaf(), value: v, unit: "K" };
    delete S.spec.T_K.range;
    markEdited();
    slider.value = Math.min(T_SLIDER[1], Math.max(T_SLIDER[0], v));
    scheduleValidate();
  }
  syncSel();
  return h("div.heatsink.temp-plate", { role: "group", "aria-label": "Temperature plate" },
    h("div.hs-label", h("span.engrave", "Temperature"), h("span.temp-chip", tagChip(leaf().tag, { size: "sm", title: leaf().source }))),
    h("div.hs-control", h("span.hs-end", `${T_SLIDER[0]} K`), slider, h("span.hs-end", `${T_SLIDER[1]} K`)),
    h("div.hs-read", valEl, h("span.hs-unit", "K")),
    h("div.hs-keys", sel));
}

// ---------------------------------------------------------------- run bar
function renderRunbar() {
  S.seg = segmented({
    label: "Run mode", value: S.mode,
    options: [
      { value: "point", label: "Point", sub: "one T" },
      { value: "tsweep", label: "T sweep", sub: `n ≤ ${TSWEEP_MAX}` },
      { value: "series", label: "Series", sub: "any parameter" },
    ],
    onChange: (v) => { S.mode = v; if (v === "series") ensureSeriesPath(); syncRunbar(); },
  });
  const num = (key, obj, label, opts = {}) => h("label.rb-f", h("span.rb-fl", label),
    h("input.fr-num.rb-in", { type: "number", step: "any", value: String(obj[key]), "aria-label": label, dataset: { f: `${label}` },
      on: { change: (e) => { const v = Number(e.target.value); if (Number.isFinite(v)) { obj[key] = opts.int ? Math.round(v) : v; clampGrid(); syncRunbar(); } } } }));
  S.dom.tsweepBox = h("div.rb-sub", { dataset: { for: "tsweep" } }, num("lo", S.tGrid, "T from (K)"), num("hi", S.tGrid, "T to (K)"), num("n", S.tGrid, "points", { int: true }));
  S.dom.seriesPath = h("select.fr-select.rb-path", { "aria-label": "Series x parameter", on: { change: (e) => { S.ser.path = e.target.value; seriesDefaults(); renderSeriesBox(); syncRunbar(); } } });
  S.dom.seriesBox = h("div.rb-sub.rb-series", { dataset: { for: "series" } });
  S.dom.eta = h("p.run-eta");
  S.dom.runBtn = h("button.run-key", { type: "button", on: { click: () => run() } }, icon("run", { size: 18 }), h("span", "RUN"));
  S.dom.cancelBtn = h("button.btn.btn-outline.cancel-btn", { type: "button", hidden: true, on: { click: cancelRun } }, icon("stop", { size: 16 }), "Cancel");
  S.dom.progress = h("div.progress", { role: "progressbar", "aria-label": "Run progress", "aria-valuemin": "0", "aria-valuemax": "100", hidden: true }, h("span.progress-fill"));
  S.dom.runStatus = h("p.run-status", { "aria-live": "polite" });
  return h("footer.runbar", { "aria-label": "Run controls" },
    h("div.rb-left", S.seg.el, S.dom.tsweepBox, S.dom.seriesBox, S.dom.eta),
    h("div.rb-right", h("div.rb-status", S.dom.progress, S.dom.runStatus), S.dom.cancelBtn, S.dom.runBtn));
}

function clampGrid() {
  S.tGrid.n = Math.max(1, Math.min(TSWEEP_MAX, S.tGrid.n));
  S.ser.n = Math.max(2, Math.min(SERIES_MAX, S.ser.n));
}

function paramSpecOf(path) {
  const [axis, name] = path.split(".");
  return entry(axis)?.params.find((p) => p.name === name) || null;
}

/** Sweepable numeric leaves of the current spec: T_K plus every numeric sweepable parameter that resolved. */
function sweepPaths() {
  const out = [{ path: "T_K", label: "T_K (K)" }];
  const resolved = S.val?.resolved || {};
  for (const a of AXES) {
    for (const p of entry(a.key).params) {
      const path = `${a.key}.${p.name}`;
      const r = resolved[path];
      if (p.sweepable && typeof p.default === "number" && r && typeof r.value === "number") out.push({ path, label: `${a.label}: ${p.name}${p.unit && p.unit !== "enum" ? ` (${unitText(p.unit)})` : ""}` });
    }
  }
  return out;
}

function ensureSeriesPath() {
  const list = sweepPaths();
  if (!list.some((x) => x.path === S.ser.path)) { S.ser.path = list[0]?.path || null; seriesDefaults(); }
  if (!S.ser.obs.length) S.ser.obs = (S.reg.observables || []).map((o) => o.id).filter((id) => DEFAULT_OBS.includes(id));
}

function seriesDefaults() {
  const path = S.ser.path;
  if (!path) return;
  if (path === "T_K") { S.ser.lo = 2; S.ser.hi = 100; S.ser.log = false; return; }
  const p = paramSpecOf(path);
  const cur = S.val?.resolved?.[path]?.value;
  if (p && isNum(p.lo) && isNum(p.hi)) { S.ser.lo = p.lo; S.ser.hi = p.hi; }
  else if (isNum(cur)) { S.ser.lo = cur * 0.5; S.ser.hi = cur * 1.5; }
  S.ser.log = !!p && p.scale === "log" && S.ser.lo > 0;
}

function refreshSeriesPaths() {
  if (!S?.dom?.seriesPath) return;
  ensureSeriesPath();
  renderSeriesBox();
}

function renderSeriesBox() {
  const list = sweepPaths();
  clear(S.dom.seriesPath);
  S.dom.seriesPath.append(...list.map((x) => h("option", { value: x.path, selected: x.path === S.ser.path }, x.label)));
  const inp = (key, label, int) => h("label.rb-f", h("span.rb-fl", label),
    h("input.fr-num.rb-in", { type: "number", step: "any", value: String(Number(Number(S.ser[key]).toPrecision(6))), "aria-label": label,
      on: { change: (e) => { const v = Number(e.target.value); if (Number.isFinite(v)) { S.ser[key] = int ? Math.round(v) : v; clampGrid(); syncRunbar(); } } } }));
  const obs = (S.reg.observables || []).flatMap((o) => (o.conventions || ["intrinsic"]).map((c, i) => ({ req: i === 0 ? o.id : `${o.id}:${c}`, label: `${o.label}${i ? ` (${c})` : ""}` })));
  const obsBox = h("details.rb-obs", h("summary.btn.btn-quiet", `Observables (${S.ser.obs.length})`),
    h("div.rb-obs-list", obs.map((o) => h("label.rb-check",
      h("input", { type: "checkbox", checked: S.ser.obs.includes(o.req), dataset: { obs: o.req },
        on: { change: (e) => { S.ser.obs = e.target.checked ? [...S.ser.obs, o.req] : S.ser.obs.filter((x) => x !== o.req); obsBox.querySelector("summary").textContent = `Observables (${S.ser.obs.length})`; syncRunbar(); } } }),
      h("span", o.label)))));
  S.dom.seriesBox.replaceChildren(h("label.rb-f", h("span.rb-fl", "x"), S.dom.seriesPath), inp("lo", "from"), inp("hi", "to"), inp("n", "points", true), obsBox);
}

function pointsFor() { return S.mode === "tsweep" ? S.tGrid.n : S.mode === "series" ? S.ser.n : 1; }

function syncRunbar() {
  if (!S?.dom?.runBtn) return;
  S.seg.select(S.mode);
  S.dom.tsweepBox.hidden = S.mode !== "tsweep";
  S.dom.seriesBox.hidden = S.mode !== "series";
  if (S.mode === "series") { ensureSeriesPath(); if (!S.dom.seriesBox.firstChild) renderSeriesBox(); }
  const ok = !!S.val?.ok;
  const running = S.runState === "running";
  S.dom.runBtn.disabled = running || !ok || (S.mode === "series" && (!S.ser.obs.length || !S.ser.path));
  S.dom.runBtn.title = !ok ? "The spec is not valid; the plaque lists why" : S.mode === "series" && !S.ser.obs.length ? "Choose at least one observable" : "";
  const eta = S.val?.eta_s;
  const n = pointsFor();
  S.dom.eta.textContent = !S.val ? "Validating"
    : !ok ? "Not runnable until the issues on the plaque are resolved"
    : isNum(eta) ? `${n} point${n === 1 ? "" : "s"}; estimated ${fmtSeconds(eta * n)} (${fmtSeconds(eta)} per point)` : `${n} point${n === 1 ? "" : "s"}`;
}

// ---------------------------------------------------------------- run
function tGridBody() {
  const { lo, hi, n } = S.tGrid;
  return { lo, hi, n };
}

function seriesValues() {
  const { lo, hi, n, log } = S.ser;
  const out = [];
  for (let i = 0; i < n; i++) {
    const f = n === 1 ? 0 : i / (n - 1);
    out.push(log && lo > 0 && hi > 0 ? lo * Math.pow(hi / lo, f) : lo + (hi - lo) * f);
  }
  return out;
}

function run() {
  if (!S || !S.val?.ok) return;
  S.job?.cancel?.();
  clampGrid();
  S.lastError = null;
  const name = S.card || "composed";
  const spec = clone(S.spec);
  setRunState("running", { eta: (S.val.eta_s || 0) * pointsFor() });
  const hooks = { onProgress: (k, n) => { if (S?.job === job) setProgress(k, n); } };
  let job;
  const mode = S.mode;
  if (mode === "series") {
    const values = seriesValues();
    job = S.ctx.api.composeSeries({ spec, path: S.ser.path, values, observables: S.ser.obs }, hooks);
  } else {
    job = S.ctx.api.composeRun(mode === "tsweep" ? { spec, mode, T_grid: tGridBody(), name } : { spec, mode, name }, hooks);
  }
  S.job = job;
  job.promise.then((r) => {
    if (!S || S.job !== job) return;
    S.job = null;
    if (mode === "series") S.series = { ...r.result, spec_path: S.ser.path, obs: [...S.ser.obs] };
    else { S.run = r.result; S.ctx.state.set({ run: r.result }); }
    setRunState("done", { elapsed: r.elapsed_s, cached: r.cached });
    renderResults({ fresh: true });
    drawBeam();
  }).catch((e) => {
    if (!S || S.job !== job) return;
    S.job = null;
    const b = e.body;
    const err = b?.error && typeof b.error === "object" ? b.error : { kind: b?.kind, message: b?.message || (typeof b?.error === "string" ? b.error : e.message) };
    if (err.kind === "cancelled") { setRunState("cancelled"); return; }
    S.lastError = { ...err, issues: b?.issues };
    setRunState("error");
    renderResults();
  });
}

function cancelRun() {
  S.job?.cancel();
  S.job = null;
  setRunState("cancelled");
}

function setProgress(k, n) {
  const pct = n ? Math.round((100 * k) / n) : 0;
  S.dom.progress.hidden = false;
  S.dom.progress.classList.remove("is-indeterminate");
  S.dom.progress.firstChild.style.transform = `scaleX(${pct / 100})`;
  S.dom.progress.setAttribute("aria-valuenow", String(pct));
  S.dom.runStatus.textContent = `Point ${k} of ${n}`;
}

function setRunState(st, info = {}) {
  S.runState = st;
  const running = st === "running";
  S.dom.cancelBtn.hidden = !running;
  S.dom.progress.hidden = !running;
  clearInterval(S.timer);
  if (running) {
    S.tt0 = performance.now();
    S.dom.progress.classList.add("is-indeterminate");
    S.dom.progress.firstChild.style.transform = "";
    S.dom.progress.removeAttribute("aria-valuenow");
    const tick = () => {
      if (!S) return;
      const el = (performance.now() - S.tt0) / 1000;
      if (!S.dom.runStatus.textContent.startsWith("Point")) {
        S.dom.runStatus.textContent = info.eta > 0 ? `Running, ${fmtSeconds(el)} of about ${fmtSeconds(info.eta)}` : `Running, ${fmtSeconds(el)}`;
        const f = !(info.eta > 0) ? 0.5 : el > info.eta ? 0.9 + 0.09 * (1 - Math.exp(-(el - info.eta) / info.eta)) : 0.9 * el / info.eta;
        S.dom.progress.firstChild.style.transform = `scaleX(${f})`;
      }
    };
    S.dom.runStatus.textContent = "Submitting";
    tick();
    S.timer = setInterval(tick, 250);
  } else if (st === "done") {
    S.dom.runStatus.textContent = info.cached ? `Done from cache in ${fmtSeconds(info.elapsed)}` : `Done in ${fmtSeconds(info.elapsed)}`;
  } else if (st === "cancelled") {
    S.dom.runStatus.textContent = "Cancelled; the last finished result stays on the shelf";
  } else if (st === "error") {
    S.dom.runStatus.textContent = "Run failed; see the shelf";
  }
  renderStatus();
  syncRunbar();
  drawBeam();
}

// ================================================================ SHELF
function renderShelf() {
  S.dom.plaqueHost = h("div.plaque-host");
  S.dom.status = h("div.shelf-status");
  S.dom.alert = h("div.shelf-alert");
  S.dom.labels = h("div.shelf-labels");
  S.readouts = {};
  S.dom.readouts = h("div.readouts", { role: "group", "aria-label": "Key readouts" },
    READOUTS.map((r) => (S.readouts[r.key] = createReadout({ key: r.key, label: r.label, onInspect: inspectKpi })).el));
  S.charts.series = chartCard({
    title: "Series", height: 280,
    build: () => buildSeries(), table: () => seriesTable(), filename: () => `${S?.run?.run_id || "fsim"}-series`,
  });
  S.charts.series.el.dataset.chart = "series";
  S.charts.tcurve = chartCard({
    title: "Against temperature", height: 280,
    build: () => buildTCurve(), table: () => tCurveTable(), filename: () => `${S?.run?.run_id || "fsim"}-T`,
  });
  S.charts.tcurve.el.dataset.chart = "tcurve";
  S.charts.hom = chartCard({
    title: "HOM visibility against photon separation", height: 240,
    build: () => buildHom(), table: () => homTable(), filename: () => `${S?.run?.run_id || "fsim"}-hom-dt`,
  });
  S.charts.hom.el.dataset.chart = "hom";
  S.dom.notes = h("details.drawer.notes-drawer", h("summary", h("span.engrave", "Model notes"), h("span.drawer-count")), h("div.drawer-body"));
  return h("section.shelf", { "aria-label": "Instrument shelf" },
    S.dom.plaqueHost,
    h("div.shelf-bar", S.dom.status),
    S.dom.alert, S.dom.labels, S.dom.readouts,
    S.charts.series.el, S.charts.tcurve.el, S.charts.hom.el, S.dom.notes);
}

function inspectKpi(key, last) {
  const s = S.run?.scalars || {};
  const spec = READOUTS.find((r) => r.key === key);
  const prov = S.run?.provenance?.[key];
  S.ctx.inspector.number({
    label: spec?.label || key, key,
    valueText: last && isNum(last.value) ? `${fmt(last.value)}${spec?.unit ? ` ${spec.unit}` : ""}` : "n/a",
    tag: last?.tag,
    notes: (prov?.inputs || []).slice(0, 40).map((p) => {
      const r = S.val?.resolved?.[p];
      return { key: p, tag: r?.tag, note: r ? `${r.value ?? JSON.stringify(r.range ?? r.choice)}${r.unit ? ` ${r.unit}` : ""}` : "" };
    }).filter((n) => n.tag),
    source: s.solver ? `fsim_core.compose engine, solver ${s.solver}` : "fsim_core.compose engine",
  });
}

// ---------------------------------------------------------------- plaque
function renderPlaque() {
  if (!S?.dom?.plaqueHost) return;
  const v = S.val;
  const run = S.run?.scalars;
  if (!v) { clear(S.dom.plaqueHost).append(h("section.plaque.plaque-none", h("div.plaque-word", icon("minus", { size: 18 }), h("span", "VALIDATING")))); return; }
  const ok = v.ok;
  const errs = (v.issues || []).filter((i) => i.level === "error");
  const warns = (v.issues || []).filter((i) => i.level !== "error");
  const runInvalid = run && run.valid === false;
  const status = ok && !runInvalid ? "pass" : "fail";
  const solver = run?.solver || v.solver;
  const regime = run?.regime || (v.regime && v.regime !== "unknown" ? v.regime : null);
  const reasons = [...errs.map((i) => ({ k: i.path, t: i.message, lvl: "error" })), ...(runInvalid ? (run.invalid_reasons || []).map((r) => ({ k: "run", t: r, lvl: "error" })) : []), ...warns.map((i) => ({ k: i.path, t: i.message, lvl: "warn" }))];
  const notes = (run?.notes || []).filter((n) => /phonon|polaron|Z\^2|Z²/i.test(n));
  clear(S.dom.plaqueHost).append(h("section.plaque", { class: `plaque-${status}`, "aria-label": "Validity", dataset: { status, ok: String(ok) } },
    h("div.plaque-top",
      h("div.plaque-word", icon(status === "pass" ? "check" : "cross", { size: 18 }), h("span", status === "pass" ? "VALID" : "INVALID")),
      h("div.plaque-denom", h("span.k", "solver"), h("span.v", solver || "none"))),
    h("div.plaque-reasons",
      h("span.plaque-kv", h("span.k", "regime"), h("span.v", regime || "known after a run")),
      isNum(v.eta_s) ? h("span.plaque-kv", h("span.k", "estimate"), h("span.v", fmtSeconds(v.eta_s))) : null),
    reasons.length ? h("ul.plaque-issues", reasons.slice(0, 8).map((r) => h("li", { class: r.lvl }, icon(r.lvl === "error" ? "cross" : "warn", { size: 13 }), h("span.pi-k", r.k), h("span.pi-t", r.t)))) : null,
    notes.length ? h("ul.plaque-issues.pi-notes", notes.slice(0, 2).map((n) => h("li", icon("info", { size: 13 }), h("span.pi-t", n)))) : null,
    h("p.plaque-source", ok ? "fsim_core.compose.validate_spec" : "Fix the listed fields or change the kind; RUN stays disabled until this is valid")));
}

function renderStatus() {
  if (!S?.dom?.status) return;
  const r = S.run;
  const word = { idle: "Ready", running: "Running", done: "Result on shelf", error: "Run failed", cancelled: "Cancelled" }[S.runState] || S.runState;
  const ico = { idle: "minus", running: "retry", done: "check", error: "warn", cancelled: "stop" }[S.runState] || "minus";
  // Node.append() stringifies null ("nullnull"): pass only real nodes.
  clear(S.dom.status).append(...[
    h("span.status-word", { dataset: { state: S.runState } }, icon(ico, { size: 16 }), word),
    r ? h("span.status-mode", r.mode === "tsweep" || (r.T_grid || []).length > 1 ? `T sweep, ${r.T_grid.length} points` : "point run") : null,
    r?.run_id ? h("span.id-chip.is-static", h("span.id-k", "RUN"), h("span.id-v", r.run_id)) : null].filter(Boolean));
}

// ---------------------------------------------------------------- results -> shelf
function renderResults() {
  if (!S?.dom?.readouts) return;
  renderStatus();
  const r = S.run;
  clear(S.dom.alert);
  if (S.lastError) {
    const issues = S.lastError.issues;
    S.dom.alert.append(errorState({
      title: "The run failed", message: [S.lastError.message, S.lastError.suggestion].filter(Boolean).join(" ") + (issues ? ` ${issues.map((i) => i.message).join("; ")}` : ""),
      actions: [{ label: "Run again", primary: true, fn: () => run() }],
    }));
  }
  renderPlaque();
  clear(S.dom.labels);
  const s = r?.scalars || {};
  const labels = r?.labels || [];
  if (labels.length) S.dom.labels.append(...labels.map((l) => h("span.label-chip", l)));
  const mid = labels.some((l) => /midpoint/i.test(l));
  for (const ro of READOUTS) {
    const R = S.readouts[ro.key];
    if (!r) { R.update(null); continue; }
    const v = s[ro.key];
    const prov = r.provenance?.[ro.key];
    const tag = normTag(prov?.tag) || normTag(r.tag_chain);
    const conv = captionFor(ro.key, s);
    R.update({
      value: v, unit: ro.unit, tag, tagNote: prov?.inputs ? `${prov.inputs.length} inputs` : null,
      qualifier: mid ? "range midpoint" : "", caption: conv, markChange: true,
      reason: (s.invalid_reasons || [])[0] || `${ro.label} not computed for this spec`,
    });
  }
  const odd = READOUTS.length % 2;
  S.dom.readouts.classList.toggle("is-odd", !!odd);
  // charts
  const tag = r ? r.tag_chain : null;
  S.charts.tcurve.setTag(tag);
  S.charts.hom.setTag(tag);
  S.charts.series.setTag(S.series ? widestOf(S.series.series?.tag) : null);
  const running = S.runState === "running";
  S.charts.series.setEmpty(running && S.mode === "series" ? loadingPlate() : emptyState({ title: "No series yet", body: "Choose Series, pick a parameter and observables, then RUN.", iconName: "chart" }));
  S.charts.tcurve.setEmpty(running && S.mode === "tsweep" ? loadingPlate() : emptyState({
    title: r && (r.curves?.T_hs || []).length <= 1 ? "Single temperature" : "No T sweep yet",
    body: r && (r.curves?.T_hs || []).length <= 1 ? "Run a T sweep to draw the curves against temperature." : "Choose T sweep and RUN.", iconName: "chart" }));
  S.charts.hom.setEmpty(emptyState({ title: "No HOM curve", body: r ? "This drive and solver return no photon-separation curve." : "Runs of pulsed coherent drives fill this plot.", iconName: "chart" }));
  S.charts.series.render();
  S.charts.tcurve.render();
  S.charts.hom.render();
  // notes
  const notes = s.notes || [];
  S.dom.notes.hidden = !notes.length;
  S.dom.notes.querySelector(".drawer-count").textContent = notes.length ? `${notes.length}` : "";
  clear(S.dom.notes.querySelector(".drawer-body")).append(...notes.map((n) => h("p.pf-note", n)));
}

function widestOf(tagMap) { return tagMap ? widest(Object.values(tagMap)) : null; }
function loadingPlate() { return h("div.chart-loading", skeleton("chart", { lines: 0 }), h("p", "Evaluating; the curve draws when the job finishes.")); }

function captionFor(key, s) {
  const obs = (S.reg.observables || []).find((o) => o.id === key);
  const conv = obs?.conventions?.[0] || "intrinsic";
  if (key === "B_first_lens") return `${conv}: ${(s.brightness_convention || "").split(";")[0]}`.slice(0, 120);
  if (key === "hom_I") return `${conv} overlap; γ* repo = 2 × textbook`;
  if (key === "g2_0") return `${conv} (no laser leak, no dark counts)`;
  if (key === "rate_detected_Hz") return "f_rep · B_fibre · η_optics · η_det";
  if (key === "F_P_eff") return "additional rate ratio; total = 1 + F_P·L";
  if (key === "tau_ps") return "1 / (Γ + k_nr)";
  return conv;
}

// ---------------------------------------------------------------- charts
function stacked(rows, { xTitle, xLog = false } = {}) {
  const c = themeColors();
  const n = rows.length;
  const gap = n > 1 ? 0.07 : 0;
  const hh = (1 - gap * (n - 1)) / n;
  const data = [];
  const layout = baseLayout({ showlegend: false, hovermode: "x unified", margin: { l: 66, r: 18, t: 14, b: 44 } });
  layout.xaxis = { ...layout.xaxis, title: { text: xTitle }, type: xLog ? "log" : "linear", anchor: `y${n > 1 ? n : ""}` };
  rows.forEach((r, i) => {
    const tag = normTag(r.tag) || "A";
    const color = c.series[i % c.series.length];
    const id = i === 0 ? "y" : `y${i + 1}`;
    data.push({
      x: r.x, y: r.y.map((v) => (isNum(v) ? v : null)), type: "scatter", mode: "lines+markers", name: r.label, yaxis: id,
      line: { color, width: tagWidth(tag), dash: tagDash(tag) }, marker: { size: 5, color },
      connectgaps: false, hovertemplate: `${r.label} %{y:.4g} [${tag}]<extra></extra>`,
    });
    const top = 1 - i * (hh + gap);
    layout[i === 0 ? "yaxis" : `yaxis${i + 1}`] = { domain: [top - hh, top], title: { text: r.label + (r.unit ? ` (${r.unit})` : ""), font: { size: 11 } }, zeroline: false, ticks: "outside", exponentformat: "e" };
  });
  return { data, layout, nrows: n };
}

function labelOf(req) {
  const [id, conv] = String(req).split(":");
  const o = (S.reg.observables || []).find((x) => x.id === id);
  return `${o?.label || id}${conv ? ` (${conv})` : ""}`;
}
function unitOfReq(req) { return unitText((S.reg.observables || []).find((x) => x.id === String(req).split(":")[0])?.unit); }

function buildSeries() {
  const sr = S?.series?.series;
  if (!sr) return null;
  const reqs = Object.keys(sr.y || {});
  const withheld = reqs.filter((q) => !normTag(sr.tag?.[q]));
  const rows = reqs.filter((q) => normTag(sr.tag?.[q])).map((q) => ({ x: sr.x, y: sr.y[q], label: labelOf(q), unit: unitOfReq(q), tag: sr.tag[q] }));
  if (!rows.length) return withheld.length ? { data: [], layout: {}, withheld: withheld.map(labelOf) } : null;
  const log = S.ser.log && S.series.spec_path === S.ser.path;
  const p = sr.path || S.series.spec_path;
  const out = stacked(rows, { xTitle: `${p}${sr.unit ? ` (${unitText(sr.unit)})` : ""}`, xLog: log });
  out.withheld = withheld.map(labelOf);
  S.charts.series.el.style.setProperty("--chart-h", `${Math.max(240, 118 * rows.length + 60)}px`);
  S.charts.series.setTitle("Series", `${p} · ${sr.x.length} points${sr.failed?.length ? ` · ${sr.failed.length} failed (gaps)` : ""} · other inputs at their central values`);
  return out;
}
function seriesTable() {
  const sr = S?.series?.series;
  if (!sr) return null;
  const reqs = Object.keys(sr.y);
  return { columns: [{ key: "x", label: sr.path }, ...reqs.map((q, i) => ({ key: `y${i}`, label: `${labelOf(q)} [${sr.tag?.[q] || "?"}]` }))], rows: sr.x.map((x, j) => Object.fromEntries([["x", x], ...reqs.map((q, i) => [`y${i}`, sr.y[q][j]])])) };
}

function curveRows() {
  const r = S?.run;
  const T = r?.curves?.T_hs;
  if (!T || T.length < 2) return null;
  const rows = [];
  for (const [key, label, unit] of [["g2", "g²(0)", ""], ["hom_I", "I_HOM", ""], ["B_first_lens", "B, first lens", ""], ["rate_detected_Hz", "Detected rate", "Hz"], ["tau_ps", "τ", "ps"]]) {
    const y = r.curves[key];
    const provKey = key === "g2" ? "g2_0" : key;
    const tag = normTag(r.provenance?.[provKey]?.tag) || normTag(r.tag_chain);
    if (Array.isArray(y) && y.length === T.length && tag) rows.push({ x: T, y, label, unit, tag });
  }
  return rows.slice(0, 4);
}
function buildTCurve() {
  const rows = curveRows();
  if (!rows?.length) return null;
  S.charts.tcurve.el.style.setProperty("--chart-h", `${Math.max(240, 112 * rows.length + 60)}px`);
  S.charts.tcurve.setTitle("Against temperature", `${S.run.run_id} · ${S.run.T_grid.length} points · other inputs at their central values`);
  return stacked(rows, { xTitle: "T (K)" });
}
function tCurveTable() {
  const rows = curveRows();
  if (!rows) return null;
  return { columns: [{ key: "T", label: "T (K)" }, ...rows.map((r, i) => ({ key: `y${i}`, label: `${r.label} [${r.tag}]` }))], rows: rows[0].x.map((t, j) => Object.fromEntries([["T", t], ...rows.map((r, i) => [`y${i}`, r.y[j]])])) };
}

function homSeries() {
  const hs = S?.run?.series?.hom_vs_dt;
  if (!hs || !Array.isArray(hs.x) || !hs.x.length) return null;
  const tag = normTag(S.run.provenance?.hom_I?.tag) || normTag(S.run.tag_chain);
  return tag ? { ...hs, tag } : null;
}
function buildHom() {
  const hs = homSeries();
  if (!hs) return null;
  S.charts.hom.setSub?.(hs.note);
  return stacked([{ x: hs.x, y: hs.y, label: "I_HOM", unit: "", tag: hs.tag }], { xTitle: `photon separation Δt (${hs.x_unit || "ns"})`, xLog: true });
}
function homTable() {
  const hs = homSeries();
  return hs ? { columns: [{ key: "x", label: `dt (${hs.x_unit || "ns"})` }, { key: "y", label: `I_HOM [${hs.tag}]` }], rows: hs.x.map((x, i) => ({ x, y: hs.y[i] })) } : null;
}

// ================================================================ VIEWPORT
async function loadScene() {
  if (!S?.dom?.viewport) return;
  const seq = ++S.sceneSeq;
  const slot = S.dom.viewport;
  let scene = null;
  try {
    scene = await S.ctx.api.composeScene(S.spec);
  } catch (e) {
    if (!S || seq !== S.sceneSeq) return;
    flatPlate(slot, null, e.message);
    return;
  }
  if (!S || seq !== S.sceneSeq) return;
  const note = scene.catalog_id ? scene.catalog_id : "";
  S.dom.vmNote.textContent = note;
  const known = !scene.fallback && scene.builder && scene.builder !== "fallback2d" && scene.variant !== "fallback2d";
  if (!known) {
    try { S.scene?.dispose?.(); } catch { /* gone */ }
    S.scene = null;
    flatPlate(slot, scene);
    return;
  }
  // An identical descriptor (the open-time scene fetch and the post-validate one) must not rebuild the stage.
  const key = JSON.stringify(scene);
  if (S.scene && S.sceneKey === key) return;
  S.sceneKey = key;
  if (S.scene) { try { S.scene.update(scene); return; } catch (e) { console.warn("scene update failed", e); try { S.scene.dispose(); } catch { /* gone */ } S.scene = null; } }
  let mod = null;
  try { mod = await import("../viz3d/index.js"); } catch { mod = null; }
  if (!S || seq !== S.sceneSeq) return;
  if (!mod?.mountScene) { flatPlate(slot, scene); return; }
  try {
    clear(slot);
    S.scene = mod.mountScene(slot, "device", scene, { theme: S.ctx.state.get("theme") });
    // a builder the viewer does not know ends up as a message panel: show the flat plate instead
    if (slot.dataset.mode === "unavailable") { S.scene.dispose(); S.scene = null; flatPlate(slot, scene); }
  } catch (e) {
    console.warn("3D mount failed", e);
    S.scene = null;
    flatPlate(slot, scene, "3D view unavailable on this machine");
  }
}

/** fallback2d: a flat plate naming the chain, with the scene's own labels and catalog ID. */
function flatPlate(slot, scene, why) {
  const chain = AXES.map((a) => ({ a, e: entry(a.key) }));
  clear(slot).append(h("div.vp-placeholder", { dataset: { mode: "fallback2d" } },
    h("div.vp-plate",
      h("span.engrave", "No 3D model for this structure"),
      h("span.vp-sub", why || scene?.reason || "The viewer has no builder for it; the chain below is the device."),
      h("div.vp-chain", chain.map(({ a, e }) => h("div.vp-link", h("span.vp-lk", a.label), h("span.vp-lv", e.label)))),
      scene?.catalog_id ? h("span.vp-sub", scene.catalog_id) : null,
      null)));
}

// ================================================================ SAVE AS COMPOSED CARD
const SLUG_RE = /^[a-z0-9]+(?:-[a-z0-9]+)*$/;
function saveAsPopover(anchor) {
  const suggest = S.card ? `${S.card}-edited` : `composed-${S.spec.structure.kind}-${S.spec.drive.kind}`.replace(/_/g, "-");
  const input = h("input.fr-num.save-name", { type: "text", id: "save-name", value: suggest.slice(0, 70), spellcheck: false, "aria-label": "Card name (lower-case slug)" });
  const msg = h("p.save-msg#save-msg");
  const go = h("button.btn.btn-outline.save-go", { type: "submit" }, icon("download", { size: 16 }), "Save");
  const check = () => {
    const v = input.value.trim();
    const shipped = S.ctx.cards().some((c) => c.name === v && !c.saved);
    const bad = !SLUG_RE.test(v) || v.length > 80;
    msg.textContent = bad ? "Use a lower-case slug: letters, digits and single hyphens." : shipped ? "That is a shipped card; it is never overwritten. Pick a new name." : `Saves cards/studio/${v}.yaml as a composed card.`;
    msg.dataset.state = bad || shipped ? "bad" : "ok";
    go.disabled = bad || shipped || !S.val?.ok;
    if (!S.val?.ok && !bad && !shipped) msg.textContent = "The spec is not valid yet; fix the plaque issues first.";
    return !go.disabled;
  };
  input.addEventListener("input", check);
  const form = h("form.save-form", {
    on: {
      submit: async (e) => {
        e.preventDefault();
        if (!check()) return;
        const name = input.value.trim();
        go.disabled = true;
        msg.textContent = "Saving";
        try {
          const r = await S.ctx.api.putJSON(`/api/cards/${encodeURIComponent(name)}`, { design: { name, platform: "composed", composed: clone(S.spec) } });
          await S.ctx.refreshCards?.();
          closePopover();
          toast(`Saved ${r.path || name}`, { kind: "ok" });
          S.edits = 0;
          if (S.dom.dirty) S.dom.dirty.hidden = true;
          S.ctx.navigate(`#/compose/${name}`);
        } catch (err) {
          msg.textContent = err.status === 409 ? `Refused (409): ${err.message}` : err.message;
          msg.dataset.state = "bad";
          go.disabled = false;
        }
      },
    },
  }, h("h2.why-title", "Save as composed card"), h("label.save-label", { for: "save-name" }, "Card name"), input, msg, h("div.save-actions", go));
  popover(anchor, form, { label: "Save as composed card" });
  check();
  input.focus();
  input.select();
}
