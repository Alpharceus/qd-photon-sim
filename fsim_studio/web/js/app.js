// FSIM Studio shell: context bar, workspace rail, inspector drawer, hash router, theme,
// Present mode, offline banner.
//
// Workspace module contract (js/workspaces/*.js):
//   export function mount(el, ctx)   -- render into el (the <main> workspace host)
//   export function unmount()        -- release charts, timers, subscriptions
//   export function update(args)     -- OPTIONAL: same workspace, new route args (e.g. a new card)
// ctx = { state, api, navigate(hash), inspector: {field, number, close}, toast, platforms,
//         platformOf(card), cards() }
// Routes: #/overview, #/design/<card>, #/compare, #/results/<id>, #/explain/<topic>,
//         #/story/<name>[/<step>], #/library[/<card>]

import * as state from "./state.js";
import * as api from "./api.js";
import { h, clear } from "./ui/dom.js";
import { icon } from "./ui/icons.js";
import { tagChip, normTag } from "./ui/tagChip.js";
import { toast, popover, closePopover } from "./ui/controls.js";
import { initInspector, inspectField, inspectNumber, closeInspector, isInspectorOpen } from "./ui/inspector.js";
import { setChartTheme, loadTemplate, rerenderCharts } from "./charts/theme.js";

const WORKSPACES = [
  { id: "overview", label: "Overview", icon: "overview", load: () => import("./workspaces/overview.js") },
  { id: "design", label: "Designer", icon: "design", load: () => import("./workspaces/designer.js") },
  { id: "compare", label: "Compare", icon: "compare", load: () => import("./workspaces/compare.js") },
  { id: "results", label: "Results", icon: "results", load: () => import("./workspaces/results.js") },
  { id: "explain", label: "Explain", icon: "explain", load: () => import("./workspaces/explain.js") },
  { id: "story", label: "Story", icon: "story", load: () => import("./workspaces/story.js") },
  { id: "library", label: "Library", icon: "library", load: () => import("./workspaces/library.js") },
];

// Fixed, ranked platform roster (DIRECTION: the roster is finite and ranked).
export const PLATFORMS = [
  { id: "edge", label: "InP edge RT", match: (c) => c.name.startsWith("edge-") },
  { id: "planar", label: "Planar nitride", match: (c) => c.platform === "ingan_gan_planar" },
  { id: "nwh", label: "Nanowire H", match: (c) => c.platform === "ingan_gan_nanowire" && c.name.includes("horizontal") },
  { id: "nwv", label: "Nanowire V", match: (c) => c.platform === "ingan_gan_nanowire" && c.name.includes("vertical") },
  { id: "legacy", label: "Legacy / V-a", match: (c) => c.platform === "legacy" && !c.name.startsWith("edge-") },
];
const DEFAULT_CARD = "edge-inp-gainp-design";

let cards = [];
let current = { ws: null, mod: null };
let presentPrevTheme = null;

const els = {
  shell: document.getElementById("shell"),
  ctx: document.getElementById("ctxbar"),
  rail: document.getElementById("rail"),
  main: document.getElementById("workspace"),
  insp: document.getElementById("inspector"),
};

export function platformOf(card) {
  const c = typeof card === "string" ? cards.find((x) => x.name === card) : card;
  if (!c) return null;
  return PLATFORMS.find((p) => p.match(c)) || null;
}

const ctx = {
  state,
  api,
  navigate: (hash) => { if (location.hash !== hash) location.hash = hash; else route(); },
  inspector: { field: inspectField, number: inspectNumber, close: closeInspector },
  toast,
  platforms: PLATFORMS,
  platformOf,
  cards: () => cards,
  // studio-p2a: re-read the card list after "Save as card" so the picker shows it at once
  refreshCards: async () => { cards = await api.getJSON("/api/cards"); syncContextBar(); return cards; },
};

// ------------------------------------------------------------------ theme
function applyTheme(mode) {
  document.documentElement.dataset.theme = mode;
  setChartTheme(mode);
  const btn = document.getElementById("theme-btn");
  if (btn) {
    clear(btn).append(icon(mode === "dark" ? "sun" : "moon", { size: 18 }));
    btn.setAttribute("aria-label", mode === "dark" ? "Switch to light theme (bead-blasted aluminium)" : "Switch to dark theme (anodized)");
    btn.title = btn.getAttribute("aria-label");
  }
}

function setPresent(on) {
  if (on === state.get("present")) return;
  if (on) {
    presentPrevTheme = state.get("theme");
    closeInspector();
    closePopover();
    state.set({ present: true, theme: "light" });
  } else {
    state.set({ present: false, theme: presentPrevTheme || "dark" });
  }
  document.documentElement.toggleAttribute("data-present", on);
  const b = document.getElementById("present-btn");
  b?.setAttribute("aria-pressed", String(on));
  if (b) b.querySelector(".btn-txt").textContent = on ? "Exit present" : "Present";
  window.dispatchEvent(new Event("resize"));
  rerenderCharts();
}

// ------------------------------------------------------------------ context bar
function buildContextBar() {
  const rosterEl = h("div.roster", { role: "tablist", "aria-label": "Platform roster (ranked)" },
    PLATFORMS.map((p, i) => h("button.roster-tab", {
      type: "button", role: "tab", dataset: { platform: p.id }, "aria-selected": "false", "aria-label": `${p.label} (rank ${i + 1})`,
      on: {
        click: () => {
          const first = cards.find((c) => c.kind === "design" && p.match(c));
          if (first) ctx.navigate(`#/design/${first.name}`);
        },
      },
    }, h("span.roster-rank", String(i + 1)), h("span.roster-label", p.label))));

  const picker = h("select.card-picker#card-picker", {
    "aria-label": "Design card",
    // Literature (parameter) cards route to the Library page, design cards to the Designer.
    on: { change: (e) => { const v = e.target.value; ctx.navigate(v.startsWith("lib:") ? `#/library/${v.slice(4)}` : `#/design/${v}`); } },
  });

  const runChip = h("button.id-chip#run-chip", {
    type: "button", title: "Run catalog ID; click to copy", "aria-label": "Run ID: none yet",
    on: {
      click: async () => {
        const id = state.get("run")?.run_id;
        if (!id) return;
        try { await navigator.clipboard.writeText(id); toast(`Copied ${id}`, { kind: "ok", timeout: 1800 }); } catch { /* clipboard blocked */ }
      },
    },
  }, h("span.id-k", "RUN"), h("span.id-v", "—"));

  const tagWrap = h("div.result-tag#result-tag", { hidden: true });

  const themeBtn = h("button.icon-btn.ctx-btn#theme-btn", { type: "button", on: { click: () => state.set({ theme: state.get("theme") === "dark" ? "light" : "dark" }) } });
  const presentBtn = h("button.btn.btn-outline.ctx-present#present-btn", {
    type: "button", "aria-pressed": "false", title: "Present mode: light theme, larger type, no chrome (Esc exits)",
    on: { click: () => setPresent(!state.get("present")) },
  }, icon("present", { size: 18 }), h("span.btn-txt", "Present"));

  els.ctx.append(
    h("a.wordmark", { href: "#/design/" + (state.get("card") || DEFAULT_CARD), "aria-label": "FSIM Studio home" },
      h("span.wm-mark", { "aria-hidden": "true" }), h("span.wm-text", "FSIM"), h("span.wm-sub", "Studio")),
    h("div.ctx-sep", { "aria-hidden": "true" }),
    rosterEl,
    h("div.ctx-card", h("span.ctx-k", "Card"), picker),
    h("div.ctx-spacer"),
    runChip,
    tagWrap,
    themeBtn,
    presentBtn,
  );
}

function syncContextBar() {
  const card = state.get("card");
  const plat = platformOf(card);
  for (const t of els.ctx.querySelectorAll(".roster-tab")) {
    const on = plat && t.dataset.platform === plat.id;
    t.setAttribute("aria-selected", on ? "true" : "false");
    t.tabIndex = 0;
  }
  const picker = document.getElementById("card-picker");
  if (picker) {
    clear(picker);
    const list = cards.filter((c) => c.kind === "design" && (!plat || plat.match(c)));
    const r = state.get("route") || {};
    const libCard = r.ws === "library" ? (r.args?.[0] || "") : null;
    if (libCard === "") picker.append(h("option", { value: "", disabled: true, selected: true, hidden: true }, "Literature card"));
    picker.append(h("optgroup", { label: "Design" }, list.map((c) => h("option", { value: c.name, selected: libCard === null && c.name === card }, c.name.replace(/-design$/, "")))));
    const lit = cards.filter((c) => c.kind === "param");
    if (lit.length) picker.append(h("optgroup", { label: "Literature" }, lit.map((c) => h("option", { value: `lib:${c.name}`, selected: libCard === c.name }, c.name))));
    picker.disabled = !list.length && !lit.length;
    // studio-p2a: an unsaved preset design is not a card; the picker says so
    if (card === UNSAVED_CARD && state.get("draft")) picker.prepend(h("option", { value: "", disabled: true, selected: true }, "unsaved design"));
  }
  syncUnsaved();
  syncRunChips();
}

// ------------------------------------------------------------------ unsaved indicator (studio-p2a)
const UNSAVED_CARD = "unsaved";
function syncUnsaved() {
  const u = state.get("unsaved");
  let chip = document.getElementById("unsaved-chip");
  if (!chip) {
    chip = h("span.chip.chip-dirty.ctx-unsaved#unsaved-chip", { hidden: true, role: "status" });
    document.getElementById("card-picker")?.closest(".ctx-card")?.after(chip);
  }
  chip.hidden = !u;
  chip.textContent = u ? (u.short || u.title) : "";
  chip.title = u ? `${u.title}: not saved as a card. Leaving asks first.` : "";
}

function syncRunChips() {
  const run = state.get("run");
  const runChip = document.getElementById("run-chip");
  if (runChip) {
    runChip.querySelector(".id-v").textContent = run?.run_id || "—";
    runChip.setAttribute("aria-label", run?.run_id ? `Run ID ${run.run_id}; click to copy` : "Run ID: none yet");
    runChip.disabled = !run?.run_id;
  }
  const tagWrap = document.getElementById("result-tag");
  if (!tagWrap) return;
  clear(tagWrap);
  const t = normTag(run?.tag_chain);
  tagWrap.hidden = !t;
  if (!t) return;
  const why = h("button.why-btn", {
    type: "button", "aria-haspopup": "dialog",
    on: { click: (e) => whyPopover(e.currentTarget, run) },
  }, `why [${t}]`);
  tagWrap.append(h("span.ctx-k", "Result"), tagChip(t, { title: "widest tag in the result's tag chain" }), why);
}

function whyPopover(anchor, run) {
  const t = normTag(run.tag_chain);
  const prov = run.provenance || {};
  const entries = Object.entries(prov).filter(([, v]) => normTag(v?.tag) === t);
  const content = h("div.why",
    h("h2.why-title", `Result tag [${t}]`),
    h("p.why-lede", `The widest tag wins: this result is [${t}] because these provenance entries are [${t}].`),
    entries.length
      ? h("ul.why-list", entries.map(([k, v]) => h("li", h("span.why-k", k), h("span.why-v", v.note || ""))))
      : h("p.why-lede", "No individual entry carries this tag in the provenance map; it comes from the design-card inputs."),
    h("p.why-foot", `${Object.keys(prov).length} provenance entries in ${run.run_id}`));
  popover(anchor, content, { label: "Why this result tag" });
}

// ------------------------------------------------------------------ rail
function buildRail() {
  els.rail.append(
    h("ul.rail-list", WORKSPACES.map((w) => h("li",
      h("a.rail-item", { href: `#/${w.id}`, dataset: { ws: w.id }, "aria-label": w.label },
        icon(w.icon, { size: 22 }), h("span.rail-label", w.label))))));
}

function syncRail(ws) {
  for (const a of els.rail.querySelectorAll(".rail-item")) {
    const on = a.dataset.ws === ws;
    if (on) a.setAttribute("aria-current", "page"); else a.removeAttribute("aria-current");
    if (a.dataset.ws === "design") a.href = `#/design/${state.get("card") || DEFAULT_CARD}`;
  }
}

// ------------------------------------------------------------------ router
function parseHash() {
  const raw = (location.hash || "").replace(/^#\/?/, "");
  const parts = raw.split("/").filter(Boolean).map(decodeURIComponent);
  // Landing route: an empty hash opens the verdict board.
  const ws = parts[0] || "overview";
  return { ws, args: parts.slice(1) };
}

let lastHash = null;
async function route() {
  // studio-p2a: a workspace with unsaved edits may hold the navigation (it asks once)
  const target = location.hash;
  if (lastHash && target !== lastHash && current.mod?.beforeLeave) {
    const ok = current.mod.beforeLeave(target);
    if (ok !== true) {
      history.replaceState(null, "", lastHash);
      if (await ok) location.hash = target;
      return;
    }
  }
  const r = parseHash();
  let ws = WORKSPACES.find((w) => w.id === r.ws) ? r.ws : "overview";
  let args = r.args;
  if (!location.hash || location.hash === "#" || location.hash === "#/") history.replaceState(null, "", "#/overview");
  if (ws === "design") {
    const draftOk = args[0] === UNSAVED_CARD && !!state.get("draft");
    if (!args[0] || (!draftOk && cards.length && !cards.some((c) => c.name === args[0] && c.kind === "design"))) {
      if (args[0]) toast(`No design card named ${args[0]}; opened ${state.get("card") || DEFAULT_CARD} instead.`, { kind: "error" });
      const fallback = state.get("card") || DEFAULT_CARD;
      history.replaceState(null, "", `#/design/${fallback}`);
      args = [fallback];
    }
    if (state.get("card") !== args[0]) state.set({ card: args[0] });
  }
  lastHash = location.hash;
  state.set({ route: { ws, args } });
  syncRail(ws);
  syncContextBar();
  document.title = `${WORKSPACES.find((w) => w.id === ws).label}${ws === "design" || (ws === "library" && args[0]) ? `: ${args[0]}` : ""} · FSIM Studio`;

  if (current.ws === ws && current.mod?.update) {
    current.mod.update(args);
    return;
  }
  try { current.mod?.unmount?.(); } catch (e) { console.warn(e); }
  closeInspector();
  clear(els.main);
  els.main.dataset.ws = ws;
  const def = WORKSPACES.find((w) => w.id === ws);
  current = { ws, mod: null };
  try {
    const mod = await def.load();
    if (current.ws !== ws) return;
    current.mod = mod;
    mod.mount(els.main, { ...ctx, args });
  } catch (e) {
    console.error(e);
    els.main.append(h("div.error-state", { role: "alert" }, h("p.err-title", `The ${def.label} workspace failed to load.`), h("p.err-msg", e.message)));
  }
}

// ------------------------------------------------------------------ offline banner
function offlineBanner(onlineNow) {
  let b = document.getElementById("offline-banner");
  if (onlineNow) { b?.remove(); return; }
  if (b) return;
  b = h("div.offline#offline-banner", { role: "alert" },
    icon("warn", { size: 18 }),
    h("span", "Lost contact with the Studio server at ", h("code", location.host), ". Results on screen stay; new runs wait until it is back."),
    h("button.btn.btn-quiet", { type: "button", on: { click: async () => { try { await api.health(); toast("Server is back.", { kind: "ok" }); } catch { /* still down */ } } } }, icon("retry", { size: 16 }), "Retry"));
  document.body.append(b);
}

// ------------------------------------------------------------------ boot
async function boot() {
  applyTheme(state.get("theme"));
  buildContextBar();
  buildRail();
  initInspector(els.insp);
  state.on("theme", applyTheme);
  state.on("run", syncRunChips);
  state.on("card", syncContextBar);
  state.on("unsaved", syncUnsaved);
  api.onReachability(offlineBanner);

  document.addEventListener("keydown", (e) => {
    if (e.key === "Escape") {
      if (state.get("present")) { setPresent(false); e.preventDefault(); return; }
      if (isInspectorOpen()) { closeInspector(); }
    }
  });

  loadTemplate(state.get("theme"));
  api.gates();
  try {
    cards = await api.listCards();
  } catch (e) {
    els.main.append(h("div.boot-fail", { role: "alert" },
      icon("warn", { size: 28 }),
      h("h1.err-title", "The Studio server is not answering"),
      h("p.err-msg", `${e.message} Expected it at ${location.origin}.`),
      h("button.btn.btn-outline", { type: "button", on: { click: () => location.reload() } }, icon("retry", { size: 16 }), "Retry")));
    return;
  }
  window.addEventListener("hashchange", route);
  await route();
}

boot();
