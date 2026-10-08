// Story mode (#/story/<name>/<step>): full screen, the baked scenes of a story JSON
// (fsim_studio/stories/<name>.json via /api/stories/<name>). Each scene: the claim headline,
// the proof visual immediately below it (chart, live 3D scene or board), the numbers as readouts
// with tag chips, the caveat PINNED on screen, catalog ID + data commit in the footer and a
// progress rail. Keys: Right/PageDown/Space next, Left/PageUp prev, Home/End, E opens the scene's
// live workspace (Esc returns), B opens the presenter view, T toggles theme, Esc exits.
// Numbers are the baked values (never typed here); chart data comes from out/ via the API.

import { h, clear, reducedMotion } from "../ui/dom.js";
import { icon } from "../ui/icons.js";
import { tagChip } from "../ui/tagChip.js";
import { toast, errorState, skeleton } from "../ui/controls.js";
import { fmt, isNum } from "../ui/format.js";
import { verdictBadge, lineStatus } from "../ui/verdictBadge.js";
import { sro, valueText, campaignTag } from "../ui/figure.js";
import { chartCard } from "../charts/chartCard.js";
import { buildSweepScatter, buildBars, buildCalibration, buildIntervals, buildDots } from "../charts/sweep.js";
import * as state from "../state.js";

const FADE_MS = 240;

// Live workspace per scene (E). Results/Explain/Designer routes; Esc comes back here.
const LIVE = {
  s01: "#/explain/va-fit", s02: "#/results/validation/figures", s03: "#/overview",
  s04: "#/results/rt_edge/explore", s05: "#/results/rt_edge/explore", s06: "#/explain/cascade",
  s07: "#/results/rt_edge/explore", s08: "#/design/edge-inp-gainp-design", s09: "#/results/nitride_cavity/explore",
  s10: "#/results/nitride_cavity/verdict", s11: "#/results/nitride_nanowire/full/explore", s12: "#/explain/band",
  s13: "#/overview",
};

let S = null;
let exitHash = "#/overview";

// ---------------------------------------------------------------- live-workspace return (module level)
// While a scene is open in its live workspace (E), Esc returns to the story. Registered once, in the
// capture phase so it runs before workspace handlers, and it yields to open drawers, popovers,
// lightboxes and the inspector.
let returnPill = null;
window.addEventListener("keydown", (e) => {
  const ret = state.get("storyReturn");
  if (e.key !== "Escape" || !ret || S) return;
  if (document.querySelector(".row-drawer, .lightbox, .popover") || !document.getElementById("inspector")?.hidden) return;
  e.preventDefault();
  e.stopPropagation();
  goBackToStory();
}, true);

function goBackToStory() {
  const ret = state.get("storyReturn");
  state.set({ storyReturn: null });
  returnPill?.remove(); returnPill = null;
  if (ret) location.hash = ret;
}

function showReturnPill() {
  returnPill?.remove();
  returnPill = h("button.story-return", { type: "button", on: { click: goBackToStory } },
    icon("story", { size: 18 }), h("span", "Back to the story"), h("kbd", "Esc"));
  document.body.append(returnPill);
}

// ---------------------------------------------------------------- lifecycle
export function mount(el, ctx) {
  returnPill?.remove(); returnPill = null;
  state.set({ storyReturn: null });
  // Present default is light; the theme in use before the story is restored on exit
  // (kept while a scene is open in its live workspace).
  if (!ctx.state.get("storyPrevTheme")) ctx.state.set({ storyPrevTheme: ctx.state.get("theme") });
  S = { el, ctx, name: null, story: null, step: 0, camps: [], health: null, layer: null, visualSeq: 0, disposers: [], presenter: null, t0: performance.now(), timer: null };
  el.classList.add("ws-story");
  document.documentElement.dataset.story = "1";
  ctx.state.set({ theme: "light" });
  const onKey = (e) => keys(e);
  document.addEventListener("keydown", onKey);
  S.disposers.push(() => document.removeEventListener("keydown", onKey));
  S.timer = setInterval(tickPresenter, 1000);
  update(ctx.args);
}

export function unmount() {
  if (!S) return;
  for (const d of S.disposers) d();
  clearInterval(S.timer);
  disposeLayer(S.layer);
  try { if (S.presenter && !S.presenter.closed) S.presenter.close(); } catch { /* */ }
  delete document.documentElement.dataset.story;
  if (!S.ctx.state.get("storyReturn")) {
    const prev = S.ctx.state.get("storyPrevTheme");
    S.ctx.state.set({ storyPrevTheme: null, ...(prev ? { theme: prev } : {}) });
  }
  S.el.classList.remove("ws-story");
  S = null;
}

export async function update(args) {
  if (!S) return;
  const name = args?.[0] || "rt-single-photons";
  const step = Math.max(1, parseInt(args?.[1] || "1", 10) || 1);
  if (name !== S.name || !S.story) {
    S.name = name;
    renderShell();
    try {
      const [story, camps, health] = await Promise.all([
        S.ctx.api.getJSON(`/api/stories/${encodeURIComponent(name)}`),
        S.ctx.api.listCampaigns().catch(() => []),
        S.ctx.api.health().catch(() => null),
      ]);
      if (!S) return;
      S.story = story; S.camps = camps; S.health = health;
    } catch (e) {
      if (!S) return;
      clear(S.dom.scene).append(h("div.st-error", errorState({ title: `Story ${name} did not load`, message: e.message, actions: [{ label: "Back to the board", fn: () => S.ctx.navigate("#/overview") }] })));
      return;
    }
    renderChrome();
  }
  const n = S.story.scenes.length;
  const idx = Math.min(step, n) - 1;
  if (step > n) history.replaceState(null, "", `#/story/${name}/${n}`);
  showScene(idx);
}

// ---------------------------------------------------------------- shell
function renderShell() {
  clear(S.el);
  S.dom = {
    top: h("header.st-top"),
    scene: h("main.st-scene", { "aria-live": "polite" }, h("div.st-loading", skeleton("title", { lines: 2 }), skeleton("chart", { lines: 0, height: 420 }))),
    foot: h("footer.st-foot"),
  };
  S.root = h("div.story", { role: "region", "aria-label": "Story mode" }, S.dom.top, S.dom.scene, S.dom.foot);
  S.el.append(S.root);
}

function renderChrome() {
  const st = S.story;
  const bake = st.stale
    ? h("div.bake.is-stale", { role: "status" },
      h("span.bake-tick", { "aria-hidden": "true" }), h("span.bake-txt", "stale bake: sources changed since ", h("time", String(st.baked_at || "").slice(0, 16).replace("T", " "))),
      h("button.btn.btn-outline.bake-btn", { type: "button", on: { click: rebake } }, icon("retry", { size: 16 }), "Re-bake"))
    : h("div.bake", { title: `baked_hash ${String(st.baked_hash || "").slice(0, 12)}` },
      h("span.bake-tick.is-ok", { "aria-hidden": "true" }), h("span.bake-txt", "baked ", h("time", String(st.baked_at || "").slice(0, 16).replace("T", " "))));
  S.dom.counter = h("span.st-counter");
  const kbtn = (label, key, ico, fn, id) => h("button.st-key", { type: "button", id, title: `${label} (${key})`, "aria-label": `${label} (${key})`, on: { click: fn } }, icon(ico, { size: 18 }), h("span.st-key-l", label), h("kbd", key));
  clear(S.dom.top).append(
    h("div.st-brand", h("span.wm-mark", { "aria-hidden": "true" }), h("span.st-wm", "FSIM"), h("span.st-story", st.title)),
    h("div.st-top-r",
      bake,
      S.dom.counter,
      kbtn("Live", "E", "results", openLive, "st-live"),
      kbtn("Presenter", "B", "present", openPresenter, "st-presenter"),
      kbtn("Theme", "T", S.ctx.state.get("theme") === "dark" ? "sun" : "moon", toggleTheme, "st-theme"),
      kbtn("Exit", "Esc", "close", exitStory, "st-exit")));
  if (st.bake_errors?.length) toast(`Bake reported ${st.bake_errors.length} error(s); affected numbers show their error.`, { kind: "error" });
}

function renderRail() {
  const scenes = S.story.scenes;
  return h("nav.st-rail", { "aria-label": "Scenes" },
    scenes.map((sc, i) => h("button.st-tick", {
      type: "button", dataset: { i: String(i) }, "aria-current": i === S.step ? "step" : null,
      "aria-label": `Scene ${i + 1}: ${disp(sc.title)}`, title: `${i + 1}. ${disp(sc.title)}`,
      on: { click: () => go(i) },
    }, h("span.st-tick-n", String(i + 1).padStart(2, "0")))));
}

// ---------------------------------------------------------------- navigation
function go(i) {
  const n = S.story.scenes.length;
  const j = Math.max(0, Math.min(n - 1, i));
  if (j === S.step && S.layer) return;
  S.ctx.navigate(`#/story/${S.name}/${j + 1}`);
}

function keys(e) {
  if (!S?.story || e.altKey || e.ctrlKey || e.metaKey) return;
  const tag = e.target?.tagName || "";
  if (/input|select|textarea/i.test(tag)) return;
  if (document.querySelector(".lightbox, .popover")) return;
  const k = e.key;
  if (k === "ArrowRight" || k === "PageDown" || k === " ") { e.preventDefault(); go(S.step + 1); }
  else if (k === "ArrowLeft" || k === "PageUp") { e.preventDefault(); go(S.step - 1); }
  else if (k === "Home") { e.preventDefault(); go(0); }
  else if (k === "End") { e.preventDefault(); go(S.story.scenes.length - 1); }
  else if (k === "e" || k === "E") { e.preventDefault(); openLive(); }
  else if (k === "b" || k === "B") { e.preventDefault(); openPresenter(); }
  else if (k === "t" || k === "T") { e.preventDefault(); toggleTheme(); }
  else if (k === "Escape") { e.preventDefault(); exitStory(); }
}

function toggleTheme() {
  const m = S.ctx.state.get("theme") === "dark" ? "light" : "dark";
  S.ctx.state.set({ theme: m });
  const b = document.getElementById("st-theme");
  if (b) b.replaceChild(icon(m === "dark" ? "sun" : "moon", { size: 18 }), b.firstChild);
  S.layer?.scene?.setTheme?.(m);
  for (const sc of S.layer?.scenes || []) sc?.setTheme?.(m);
}

function exitStory() { S.ctx.navigate(exitHash); }

function openLive() {
  const sc = S.story.scenes[S.step];
  const target = LIVE[sc.id] || "#/overview";
  S.ctx.state.set({ storyReturn: `#/story/${S.name}/${S.step + 1}` });
  S.ctx.navigate(target);
  showReturnPill();
}

async function rebake() {
  const btn = S.dom.top.querySelector(".bake-btn");
  if (btn) { btn.disabled = true; btn.setAttribute("aria-busy", "true"); btn.lastChild.textContent = "Baking (about 10 s)"; }
  try {
    const st = await S.ctx.api.postJSON(`/api/stories/${encodeURIComponent(S.name)}/bake`, {});
    if (!S) return;
    S.story = { ...st, stale: false };
    const fresh = await S.ctx.api.getJSON(`/api/stories/${encodeURIComponent(S.name)}`).catch(() => null);
    if (fresh) S.story = fresh;
    renderChrome();
    showScene(S.step, { force: true });
    toast(`Re-baked in ${fmt(st.bake_seconds)} s.`, { kind: "ok" });
  } catch (e) {
    toast(`Bake failed: ${e.message}`, { kind: "error" });
    if (btn) { btn.disabled = false; btn.removeAttribute("aria-busy"); btn.lastChild.textContent = "Re-bake"; }
  }
}

// ---------------------------------------------------------------- scene
const disp = (t) => String(t || "").replace(/\bg2\(/g, "g²(").replace(/\bg2\b/g, "g²").replace(/\+-/g, "±");

function dataCommits(sc) {
  const files = new Set();
  let live = false;
  for (const n of sc.numbers || []) {
    if (n.source?.kind === "live") { live = true; continue; }
    const f = n.source?.file || String(n.resolved_from || "").split(":")[0];
    if (f) files.add(f.replace(/\\/g, "/"));
  }
  const spec = sc.visual?.spec || {};
  for (const c of [].concat(spec.campaign || [], spec.campaigns || [])) files.add(`out/${c}/`);
  if (spec.data) files.add(spec.data);
  const out = new Map();
  for (const f of files) {
    const camp = [...S.camps].sort((a, b) => b.path.length - a.path.length).find((c) => f.startsWith(`${c.path}/`) || f === c.path || f.startsWith(c.path));
    if (camp) out.set(camp.path, camp.commit);
  }
  const chips = [...out.entries()].map(([p, c]) => ({ k: "DATA", v: c || "n/a", t: p }));
  if (live) chips.push({ k: "LIVE", v: `fsim_core ${String(S.health?.fsim_core_hash || "").slice(0, 7)}`, t: "live fsim_core.device.evaluate at bake time" });
  return chips;
}

function showScene(i, { force = false } = {}) {
  const sc = S.story.scenes[i];
  const first = !S.layer;
  const same = i === S.step && !first;
  S.step = i;
  if (same && !force) return;
  S.dom.counter.textContent = `${String(i + 1).padStart(2, "0")} / ${String(S.story.scenes.length).padStart(2, "0")}`;
  const isBoard = sc.visual?.kind === "board";
  // claim + caveat swap instantly (text), the proof visual crossfades
  if (!S.dom.claim) {
    S.dom.claim = h("div.st-claim");
    S.dom.proof = h("section.st-proof", { "aria-label": "Proof" });
    S.dom.nums = h("aside.st-numbers", { "aria-label": "Numbers" });
    S.dom.body = h("div.st-body", S.dom.proof, S.dom.nums);
    S.dom.caveat = h("div.st-caveat", { role: "note" });
    clear(S.dom.scene).append(S.dom.claim, S.dom.body, S.dom.caveat);
  }
  clear(S.dom.claim).append(
    h("h1.st-title", disp(sc.title)),
    h("p.st-lede", disp(sc.claim)));
  S.dom.body.classList.toggle("is-board", isBoard);
  clear(S.dom.nums);
  if (!isBoard) S.dom.nums.append(...(sc.numbers || []).map((n) => (isStatusWord(n.value) ? statusRo(n) : numberRo(n))));
  S.dom.nums.classList.toggle("is-dense", (sc.numbers || []).length > 7);
  clear(S.dom.caveat).append(h("span.cv-ico", icon("warn", { size: 20 })), h("span.cv-k", "Caveat"), h("p.cv-text", disp(sc.caveat || "No caveat recorded for this scene.")));
  // footer
  clear(S.dom.foot).append(
    h("div.st-cat",
      h("span.id-chip.is-static", h("span.id-k", "CAT"), h("span.id-v", sc.catalog_id || sc.id)),
      ...dataCommits(sc).map((c) => h("span.id-chip.is-static", { title: c.t }, h("span.id-k", c.k), h("span.id-v", c.v), h("span.id-p", c.t.replace(/^out\//, "")))),
      h("span.st-scene-id", sc.id)),
    renderRail());
  // visual with crossfade
  const old = S.layer;
  const layer = { el: h("div.st-layer"), scenes: [], charts: [] };
  S.layer = layer;
  const seq = ++S.visualSeq;
  const fade = !first && !reducedMotion();
  if (fade) layer.el.classList.add("is-in");
  S.dom.proof.append(layer.el);
  renderVisual(sc, layer, seq).catch((e) => {
    if (S?.visualSeq !== seq) return;
    clear(layer.el).append(h("div.st-error", errorState({ title: "The proof visual did not render", message: e.message })));
  });
  if (old) {
    if (fade) {
      requestAnimationFrame(() => { layer.el.classList.remove("is-in"); old.el.classList.add("is-out"); });
      setTimeout(() => { disposeLayer(old); old.el.remove(); }, FADE_MS + 40);
    } else { disposeLayer(old); old.el.remove(); }
  }
  renderPresenter();
}

function disposeLayer(l) {
  if (!l) return;
  for (const c of l.charts) { try { c.dispose(); } catch { /* */ } }
  for (const s of l.scenes) { try { s.dispose(); } catch { /* */ } }
  l.charts = []; l.scenes = [];
}

function numberRo(n) {
  return sro({ label: disp(n.label), value: n.value, raw: n.raw, unit: n.unit, tag: n.tag, error: n.error || null,
    note: "", title: n.resolved_from ? `${n.resolved_from}${n.raw != null ? ` · raw ${n.raw}` : ""}` : null, key: n.label });
}

// A printed status token (e.g. pass_hardware_infeasible) never renders raw in a readout (R4):
// the readout says it in room language and the token stays in the tooltip. The plate keeps a
// single status chip (F7), so no second badge here.
const STATUS_SAY = { pass_hardware_infeasible: "passes optically, fails in hardware", no_idealized_pass: "no optical pass" };
const sayStatus = (v) => STATUS_SAY[v] || String(v).replace(/_/g, " ");
const isStatusWord = (v) => typeof v === "string" && /^[a-z]+(?:_[a-z]+)+$/.test(v);
function statusRo(n) {
  const el = numberRo(n);
  el.classList.add("sro-status");
  const val = el.querySelector(".sro-value");
  if (val && el.dataset.state === "value") {
    val.textContent = sayStatus(n.value);
    el.title = `${n.resolved_from || ""} · printed word ${n.value}`;
    const st = lineStatus(n.value).status;
    el.dataset.status = st;
  }
  return el;
}

function findNumber(re, sceneId = null) {
  for (const sc of S.story.scenes) {
    if (sceneId && sc.id !== sceneId) continue;
    const n = (sc.numbers || []).find((x) => re.test(x.label));
    if (n) return n;
  }
  return null;
}
// Gates (g2 ceiling, flux floor, gate box) come only from GET /api/gates (I3).
function gate() {
  const gs = S.ctx.api.gatesNow() || {};
  const f = gs.flux_floor || {}, g = gs.g2_ceiling || {};
  const t = (x) => (x ? String(x).replace(/[[\]]/g, "") : "A");
  return { flux: isNum(f.value) ? f.value : null, fluxTag: t(f.tag), g2: isNum(g.value) ? g.value : null, g2Tag: t(g.tag) };
}

// ---------------------------------------------------------------- visuals
function chartIn(layer, opts) {
  const c = chartCard({ ...opts, filename: () => `fsim-story-${S.story.scenes[S.step].id}` });
  c.el.classList.add("st-chart");
  layer.charts.push(c);
  return c;
}

async function sweepRows(campaign, { file = null, filter = null, headline = false, cols }) {
  const qs = new URLSearchParams({ cols: cols.join(","), limit: "50000", headline: headline ? "1" : "0" });
  if (file) qs.set("file", file);
  if (filter) qs.set("filter", JSON.stringify(filter));
  const r = await S.ctx.api.getJSON(`/api/campaigns/${campaign}/sweep?${qs}`);
  return { rows: r.rows.map((row, i) => ({ ...Object.fromEntries(r.columns.map((c, j) => [c, row[j]])), __v: i })), total: r.total, filtered: r.filtered };
}

async function mount3D(host, kind, params, layer, seq) {
  clear(host).append(h("div.vp-placeholder", h("div.vp-plate", h("span.engrave", "3D view loading"), h("span.vp-sub", params.card || ""))));
  const spec = await S.ctx.api.sceneSpec(kind, params);
  if (!S || S.visualSeq !== seq) return null;
  const mod = await import("../viz3d/index.js");
  if (!S || S.visualSeq !== seq) return null;
  clear(host);
  const sc = mod.mountScene(host, kind, spec, { theme: S.ctx.state.get("theme") });
  layer.scenes.push(sc);
  return { sc, spec };
}

async function renderVisual(sc, layer, seq) {
  const v = sc.visual || {};
  const spec = v.spec || {};
  const L = layer.el;
  const nums = sc.numbers || [];
  await S.ctx.api.gates();
  if (S?.visualSeq !== seq) return;
  const g = gate();
  const preset = spec.preset || spec.scene;
  if (v.kind === "chart" && preset === "calibration_g2T") {
    const va = await S.ctx.api.getJSON("/api/explain/va_fit");
    if (S?.visualSeq !== seq) return;
    const c = chartIn(layer, { title: "Measured g²(T) and the calibrated model", subtitle: `${spec.data} · ${spec.model} · residual strip below`, build: () => buildCalibration(va),
      table: () => ({ columns: [{ key: "T_K", label: "T (K)" }, { key: "g2", label: "measured" }, { key: "bound", label: "bound" }, { key: "model_g2", label: "model" }, { key: "residual", label: "residual" }], rows: va.points }) });
    c.setTag(String(va.tag_chain || "").replace(/[[\]]/g, ""));
    L.append(c.el); await c.render();
  } else if (v.kind === "chart" && preset === "band_vs_points") {
    const camp = spec.data.split("/")[1];
    const file = spec.data.split("/").pop();
    const { rows } = await sweepRows(camp, { file, filter: spec.filter, cols: ["T_K", spec.points, spec.lo, spec.hi] });
    if (S?.visualSeq !== seq) return;
    const c = chartIn(layer, { title: "Model envelope at each measured temperature", subtitle: `${spec.data} · intervals lo – hi, no midpoint curve`,
      build: () => buildIntervals({ T: rows.map((r) => r.T_K), lo: rows.map((r) => r[spec.lo]), hi: rows.map((r) => r[spec.hi]), measured: rows.map((r) => r[spec.points]), ceiling: isNum(g.g2) ? { value: g.g2, tag: g.g2Tag } : null }),
      table: () => ({ columns: [{ key: "T_K", label: "T (K)" }, { key: spec.points, label: "measured" }, { key: spec.lo, label: "model lo" }, { key: spec.hi, label: "model hi" }], rows }) });
    c.setTag("A");
    L.append(c.el); await c.render();
  } else if (v.kind === "scene3d" && spec.scene === "device") {
    const cards = spec.cards || [spec.card];
    const host = h("div.st-viewport", { role: "region", "aria-label": "3D device" });
    const cap = h("span.vm-note");
    const tabs = h("div.hero-tabs", { role: "tablist", "aria-label": "Device" }, cards.map((c, i) => h("button.hero-tab", {
      type: "button", role: "tab", dataset: { card: c }, "aria-selected": String(i === 0),
      on: { click: async () => { for (const b of tabs.children) b.setAttribute("aria-selected", String(b.dataset.card === c)); await swap(c); } },
    }, c.replace(/-design$/, "").replace(/^nitride-/, ""))));
    L.append(h("section.viewport-mount.st-vm", h("div.vm-rail", h("span.engrave", "Device, true scale"), tabs, cap), host));
    const swap = async (card) => {
      for (const s of layer.scenes) { try { s.dispose(); } catch { /* */ } }
      layer.scenes = [];
      const r = await mount3D(host, "device", { card }, layer, seq);
      if (r) cap.textContent = r.spec.catalog_id || "";
    };
    await swap(cards[0]);
  } else if (v.kind === "chart" && preset === "bars_vs_floor") {
    const bars = nums.map((n) => ({ m: /at (\d+(?:\.\d+)?) K/.exec(n.label), n })).filter((x) => x.m && /flux max/i.test(x.n.label))
      .map(({ m, n }) => ({ label: `${m[1]} K`, value: n.value, tag: n.tag }));
    const c = chartIn(layer, { title: "Best collected flux per heatsink temperature, headline model", subtitle: "baked from out/rt_edge/sweep.csv (headline rows) · floor from /api/gates",
      build: () => buildBars({ bars, yTitle: "Collected flux, max (s⁻¹)", floor: g.flux, floorTag: g.fluxTag, floorLabel: `flux floor ${fmt(g.flux)} s⁻¹ [${g.fluxTag}]` }),
      table: () => ({ columns: [{ key: "label", label: "T_hs" }, { key: "value", label: "flux max (s⁻¹)" }, { key: "tag", label: "tag" }], rows: bars }) });
    c.setTag("A");
    L.append(c.el); await c.render();
  } else if (v.kind === "chart" && preset === "g2_vs_flux_scatter") {
    const cols = [spec.x, spec.y, spec.color, spec.facet, "valid", "invalid_reasons", "diagnostic_valid"].filter(Boolean);
    const campCols = (await S.ctx.api.getJSON(`/api/campaigns/${spec.campaign}`)).columns.map((c) => c.name);
    const { rows, total, filtered } = await sweepRows(spec.campaign, { headline: !!spec.headline, cols: cols.filter((c) => campCols.includes(c)) });
    if (S?.visualSeq !== seq) return;
    const validKey = campCols.includes("valid") ? "valid" : campCols.includes("diagnostic_valid") ? "diagnostic_valid" : null;
    const c = chartIn(layer, {
      title: `${spec.y === "g2_pulsed" ? "Pulsed g²(0)" : "g²(0)"} against collected flux · ${spec.headline ? `${filtered} headline rows of ${total}` : `${filtered} rows`}`,
      subtitle: `out/${spec.campaign}/sweep.csv${spec.facet ? ` · facet ${spec.facet}` : ""} · gate box from /api/gates`,
      build: () => buildSweepScatter({ rows, x: spec.x, y: spec.y, color: spec.color, facet: spec.facet, logX: true, validKey, reasonKey: "invalid_reasons",
        tag: "A", refs: { ceiling: g.g2, ceilingTag: g.g2Tag, floor: g.flux, floorTag: g.fluxTag, gate: isNum(g.flux) && isNum(g.g2) ? { g2: g.g2, flux: g.flux } : null } }),
      table: () => ({ columns: [spec.x, spec.y, spec.color, spec.facet].filter(Boolean).map((k) => ({ key: k, label: k })), rows: rows.slice(0, 2000) }) });
    c.setTag("A");
    L.append(c.el); await c.render();
  } else if (v.kind === "scene3d" && spec.scene === "cascade") {
    const host = h("div.st-viewport", { role: "region", "aria-label": "3D cascade" });
    const cap = h("span.vm-note");
    L.append(h("section.viewport-mount.st-vm", h("div.vm-rail", h("span.engrave", "XX → X → 0 inside one drive period"), cap), host));
    const r = await mount3D(host, "cascade", { card: spec.card }, layer, seq);
    if (r) cap.textContent = `${spec.card} · ${r.spec.catalog_id || ""}`;
  } else if (v.kind === "board" && preset === "model_2x2") {
    L.append(model2x2(nums));
  } else if (v.kind === "chart" && preset === "hom_vs_T") {
    const pick = (re) => nums.map((n) => ({ m: /at (\d+) K/.exec(n.label), n })).filter((x) => x.m && re.test(x.n.label)).map(({ m, n }) => ({ x: Number(m[1]), y: n.value, tag: n.tag }));
    const panels = [{ name: "HOM indistinguishability M", yTitle: "M", points: pick(/HOM/i) }, { name: "Filtered CW g²(0) shift (corr)", yTitle: "Δg²(0)", points: pick(/filtered CW/i) }];
    const c = chartIn(layer, { title: "Lindblad diagnostics at two heatsink temperatures", subtitle: "live fsim_core at bake time · markers only: two temperatures are not a curve",
      build: () => buildDots(panels), table: () => ({ columns: [{ key: "q", label: "quantity" }, { key: "x", label: "T_hs (K)" }, { key: "y", label: "value" }], rows: panels.flatMap((p) => p.points.map((q) => ({ q: p.name, ...q }))) }) });
    c.setTag("DR");
    L.append(c.el); await c.render();
  } else if (v.kind === "chart" && preset === "paired_bars") {
    const bars = nums.map((n) => ({ m: /s=([\d.]+)/.exec(n.label), n })).filter((x) => x.m).map(({ m, n }) => ({ label: `s = ${m[1]}`, value: n.value, tag: n.tag }));
    const c = chartIn(layer, { title: "SET collected flux at 230 K under three screening scenarios", subtitle: "out/nitride_cavity · log scale · scenarios side by side, never averaged",
      build: () => buildBars({ bars, yTitle: "Collected flux (s⁻¹)", logY: true }), table: () => ({ columns: [{ key: "label", label: "scenario" }, { key: "value", label: "flux (s⁻¹)" }, { key: "tag", label: "tag" }], rows: bars }) });
    c.setTag("A");
    L.append(c.el); await c.render();
  } else if (v.kind === "chart" && preset === "delivered_vs_commanded") {
    const { rows } = await sweepRows("nitride_nanowire/full", { filter: { row_kind: "core" }, cols: [spec.x, spec.y, "family", "strain_bound", "valid", "invalid_reasons", "row_id"] });
    if (S?.visualSeq !== seq) return;
    const c = chartIn(layer, { title: "Delivered against commanded flux, nanowire core rows", subtitle: "out/nitride_nanowire/full/sweep.csv · log–log · y = x is perfect delivery · strain bounds in separate panels",
      build: () => buildSweepScatter({ rows, x: spec.x, y: spec.y, color: "family", facet: "strain_bound", logX: true, logY: true, validKey: "valid", reasonKey: "invalid_reasons", idKey: "row_id", tag: "A", nanowire: true, refs: { identity: true } }),
      table: () => ({ columns: [{ key: "row_id", label: "row" }, { key: spec.x, label: "commanded" }, { key: spec.y, label: "delivered" }, { key: "family", label: "family" }, { key: "strain_bound", label: "bound" }], rows: rows.slice(0, 2000) }) });
    c.setTag("A");
    L.append(c.el); await c.render();
  } else if (v.kind === "scene3d" && spec.scene === "band") {
    const host = h("div.st-viewport", { role: "region", "aria-label": "3D band profile" });
    const bars = [
      ...nums.filter((n) => /planar E_C\/kT/i.test(n.label)).map((n) => ({ label: n.label.replace(/^planar E_C\/kT, island /, "planar ").replace(/, /g, "<br>"), value: n.value, tag: n.tag })),
      ...nums.filter((n) => /nanowire E_C\/kT/i.test(n.label)).map((n) => ({ label: n.label.replace(/^nanowire E_C\/kT \((lo|hi)\)/, "nanowire<br>range $1"), value: n.value, tag: n.tag })),
    ];
    const req = nums.find((n) => /required E_C\/kT/i.test(n.label));
    const c = chartIn(layer, { title: "Charging energy against the thermal margin", subtitle: "E_C / kT · nanowire lo and hi are range edges, not averaged",
      build: () => buildBars({ bars, yTitle: "E_C / kT", floor: req?.value, floorTag: req?.tag || "A", floorLabel: `required ≥ ${fmt(req?.value)} [${req?.tag || "A"}]` }),
      table: () => ({ columns: [{ key: "label", label: "case" }, { key: "value", label: "E_C/kT" }, { key: "tag", label: "tag" }], rows: bars }) });
    c.setTag("A");
    const cap = h("span.vm-note");
    L.append(h("div.st-split", h("section.viewport-mount.st-vm", h("div.vm-rail", h("span.engrave", "Band profile"), cap), host), c.el));
    await c.render();
    const r = await mount3D(host, "band", { card: spec.card }, layer, seq);
    if (r) cap.textContent = `${spec.card} · ${r.spec.catalog_id || ""}`;
  } else if (v.kind === "board" && preset === "scorecard") {
    const ids = ["rt_edge", "nitride_cavity", "nitride_geometry_stark", "nitride_nanowire/full"];
    const tags = Object.fromEntries(await Promise.all(ids.map(async (id) => [id, await campaignTag(id)])));
    if (S?.visualSeq !== seq) return;
    L.append(scorecard(nums, tags));
  } else {
    L.append(h("div.st-error", errorState({ title: "No renderer for this visual", message: `${v.kind} / ${preset || "?"}` })));
  }
}

// ---------------------------------------------------------------- boards
function model2x2(nums) {
  const cells = { "static, no density": [0, 0], "static + density": [0, 1], headline: [1, 0], "finite + density": [1, 1] };
  const by = {};
  for (const n of nums) {
    const [pre, ...rest] = n.label.split(":");
    if (!(pre in cells)) continue;
    (by[pre] ||= []).push({ ...n, label: rest.join(":").trim() });
  }
  const cell = (pre) => {
    const isHead = pre === "headline";
    return h("div.m2-cell", { class: isHead ? "is-headline" : "", dataset: { cell: pre } },
      h("div.m2-head", h("span.m2-name", pre === "headline" ? "headline model" : pre), isHead ? h("span.headline-chip", icon("pin", { size: 14 }), "the gating model") : h("span.m2-never", "never gates")),
      h("div.m2-nums", (by[pre] || []).map((n) => numberRo(n))));
  };
  return h("div.board.m2", { role: "table", "aria-label": "Model switch grid" },
    h("div.m2-corner"),
    h("div.m2-col", "ret.tau_cap_scales_with_density = false"), h("div.m2-col", "= true"),
    h("div.m2-row", "drive.finite_pulse = false"), cell("static, no density"), cell("static + density"),
    h("div.m2-row", "= true"), cell("headline"), cell("finite + density"));
}

function scorecard(nums, tags) {
  const N = (re) => nums.find((n) => re.test(n.label));
  const camp = (id) => S.camps.find((c) => c.id === id);
  // one VERDICT line per plate, selected by its printed scope; fields shown verbatim
  const line = (id, pred) => (camp(id)?.verdicts || []).find((v) => pred(v.fields || {}));
  const SC_KEYS = ["g2_min", "flux_max", "eligible", "paired_optical_pass", "quality_pass", "hardware_qualified", "hardware_pass_count", "coverage"];
  const lineBox = (id, pred) => {
    const v = line(id, pred);
    if (!v) return null;
    const t = tags[id] || {};
    const keys = SC_KEYS.filter((k) => v.fields?.[k] != null).slice(0, 4);
    // one status per plate: the plate's word is the only badge; the line's own word and
    // idealized_status are quoted as text, never as a second (or compound) chip
    const f = v.fields || {};
    return h("section.sc-line",
      h("p.sc-line-h", h("span.sc-line-k", `${id} · VERDICT line ${v.line}`),
        // printed tokens stay verbatim, but only inside neutral token chips (never in prose)
        h("span.sc-line-w", h("span.tok-chip", `word ${v.word}`),
          f.idealized_status && f.idealized_status !== v.word ? h("span.tok-chip", `idealized_status ${f.idealized_status}`) : null)),
      h("div.sc-line-nums", keys.map((k) => sro({ label: k.replace(/_/g, " "), value: isNum(Number(v.fields[k])) && !/\//.test(v.fields[k]) ? Number(v.fields[k]) : v.fields[k], raw: v.fields[k], tag: t.tag, note: t.note, size: "sm", key: `${id}:${k}` }))));
  };
  const plate = (title, sub, wordNum, extra, campIds, lines = []) => h("article.sc-plate",
    h("header.sc-head", h("h2.sc-title", title), h("p.sc-sub", sub)),
    wordNum ? h("div.sc-word", verdictBadge(String(wordNum.value), { size: "lg" }), tagChip(wordNum.tag, { size: "sm", title: wordNum.resolved_from })) : null,
    h("div.sc-nums", extra.filter(Boolean).map((n) => (isStatusWord(n.value) ? statusRo(n) : numberRo(n)))),
    ...lines.filter(Boolean),
    h("footer.sc-foot", campIds.map((id) => camp(id)).filter(Boolean).map((c) => h("span.id-chip.is-static", { title: c.path }, h("span.id-k", "DATA"), h("span.id-v", c.commit), h("span.id-p", c.id)))));
  return h("div.board.sc",
    plate("InP edge emitter, RT", "finite-pulse headline model", N(/RT edge/i), [], ["rt_edge"],
      [lineBox("rt_edge", () => true)]),
    plate("Planar InGaN cavity", "SET loading, c-plane and screened geometry", N(/planar cavity SET$/i), [N(/planar cavity SET hardware/i), N(/geometry\/Stark/i)], ["nitride_cavity", "nitride_geometry_stark"],
      [lineBox("nitride_cavity", (f) => f.regime === "deterministic_pair")]),
    plate("Nitride nanowire", "vertical photonic wire, relaxed bound, SET 80 MHz", N(/nanowire vertical relaxed SET 80 MHz$/i), [N(/hardware qualified/i)], ["nitride_nanowire/full"],
      [lineBox("nitride_nanowire/full", (f) => f.family === "vertical_photonic" && f.regime === "deterministic_pair" && f.strain_bound === "relaxed" && Number(f.rep_rate_hz) === 8e7)]));
}

// ---------------------------------------------------------------- presenter view (B)
function openPresenter() {
  if (S.presenter && !S.presenter.closed) { S.presenter.focus(); renderPresenter(); return; }
  const w = window.open("", "fsim-presenter", "popup,width=1180,height=760");
  if (!w) { toast("The presenter window was blocked. Allow pop-ups for this page, then press B again.", { kind: "error" }); return; }
  S.presenter = w;
  S.t0 = performance.now();
  const base = location.origin;
  w.document.open();
  w.document.write(`<!doctype html><html lang="en" data-theme="dark"><head><meta charset="utf-8"><title>FSIM presenter</title>
<link rel="stylesheet" href="${base}/css/tokens.css"><link rel="stylesheet" href="${base}/css/app.css"></head>
<body class="presenter"><div id="pv"></div></body></html>`);
  w.document.close();
  w.addEventListener("keydown", (e) => keys(e));
  renderPresenter();
}

function tickPresenter() {
  const w = S?.presenter;
  if (!w || w.closed) return;
  const t = w.document.getElementById("pv-timer");
  if (t) t.textContent = mmss((performance.now() - S.t0) / 1000);
  const c = w.document.getElementById("pv-clock");
  if (c) c.textContent = new Date().toLocaleTimeString([], { hour: "2-digit", minute: "2-digit" });
}

const mmss = (s) => `${String(Math.floor(s / 60)).padStart(2, "0")}:${String(Math.floor(s % 60)).padStart(2, "0")}`;

function renderPresenter() {
  const w = S?.presenter;
  if (!w || w.closed || !S.story) return;
  const doc = w.document;
  const root = doc.getElementById("pv");
  if (!root) return;
  const scenes = S.story.scenes;
  const cur = scenes[S.step];
  const nxt = scenes[S.step + 1];
  const el = (tag, cls, text) => { const e = doc.createElement(tag); if (cls) e.className = cls; if (text != null) e.textContent = text; return e; };
  root.replaceChildren();
  const head = el("header", "pv-head");
  head.append(el("span", "pv-timer", mmss((performance.now() - S.t0) / 1000)), el("span", "pv-step", `Scene ${S.step + 1} of ${scenes.length}`), el("span", "pv-clock", ""), el("span", "pv-keys", "→ / Space next · ← previous · E live · T theme"));
  head.querySelector(".pv-timer").id = "pv-timer";
  head.querySelector(".pv-clock").id = "pv-clock";
  const now = el("section", "pv-now");
  now.append(el("h1", "pv-title", disp(cur.title)), el("p", "pv-claim", disp(cur.claim)));
  const cav = el("p", "pv-caveat", `Caveat: ${disp(cur.caveat || "none recorded")}`);
  const notes = el("div", "pv-notes");
  notes.append(el("h2", "pv-h", "Notes"), el("p", "pv-note", cur.notes || "No presenter notes in the story JSON for this scene. The claim, every number and the caveat are on the main screen."));
  const numsBox = el("div", "pv-nums");
  numsBox.append(el("h2", "pv-h", "Numbers on screen"));
  for (const n of cur.numbers || []) {
    const r = el("p", "pv-num");
    r.append(el("span", "pv-num-l", disp(n.label)), el("span", "pv-num-v", `${n.error ? `error: ${n.error}` : isStatusWord(n.value) ? sayStatus(n.value) : valueText(n.value, n.raw)}${n.unit ? ` ${n.unit}` : ""}`), el("span", "pv-num-t", `[${n.tag}]`));
    numsBox.append(r);
  }
  now.append(cav, notes, numsBox);
  const next = el("aside", "pv-next");
  next.append(el("h2", "pv-h", "Next"));
  if (nxt) next.append(el("p", "pv-next-t", `${S.step + 2}. ${disp(nxt.title)}`), el("p", "pv-next-c", disp(nxt.claim)), el("p", "pv-next-k", `proof: ${nxt.visual?.kind} · ${nxt.visual?.spec?.preset || nxt.visual?.spec?.scene || ""}`), el("p", "pv-next-k", `caveat: ${disp(nxt.caveat)}`));
  else next.append(el("p", "pv-next-t", "End of the story."));
  root.append(head, el("div", "pv-grid"));
  root.lastChild.append(now, next);
  doc.title = `FSIM presenter · ${S.step + 1}/${scenes.length}`;
  tickPresenter();
}

// remember where the story was entered from, so Exit returns there
window.addEventListener("hashchange", (e) => {
  const from = String(e.oldURL || "").split("#")[1] || "";
  if (!from.startsWith("/story") && String(e.newURL || "").includes("#/story")) exitHash = from ? `#${from}` : "#/overview";
});
