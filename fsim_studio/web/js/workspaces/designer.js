// Designer workspace: the optical breadboard (left ~60%) and the instrument shelf (right ~40%).
// Breadboard: chain mounts joined by the beam, 3D device viewport at the centre, Heatsink plate
// beneath, block editor sheet opening in place, run controls bottom-right with the RUN key.
// Shelf: verdict plaque, KPI readouts with tag chips, hero g2(T), secondary plots, detail drawer.
// No physics here: every number is from /api/run (fsim_core) or /api/campaigns (out/).

import { h, clear, append, svgEl, debounce, reducedMotion } from "../ui/dom.js";
import { icon } from "../ui/icons.js";
import { tagChip, normTag, widest } from "../ui/tagChip.js";
import { createReadout } from "../ui/readout.js";
import { verdictPlaque } from "../ui/verdictPlaque.js";
import { sro } from "../ui/figure.js";
import { segmented, toast, skeleton, emptyState, errorState, popover, closePopover } from "../ui/controls.js";
import { fieldRow } from "../ui/fieldRow.js";
import { fmt, fmtInterval, fmtSeconds, isNum, unitText, SCALAR_UNITS } from "../ui/format.js";
import { chartCard } from "../charts/chartCard.js";
import { buildG2T, buildSecondary, g2TableRows, SECONDARY } from "../charts/g2T.js";
import { colors as themeColors } from "../charts/theme.js";
import { lineScope, verdictBadge } from "../ui/verdictBadge.js";
import { SEEDABLE, rtEdgeRowToDesign, rtEdgeSeedFields, rtEdgePairs, liveNanReason } from "./rowSeed.js";
import { createTimeline } from "../ui/timeline.js";
import { playheadXY } from "../charts/g2T.js";

// ---------------------------------------------------------------- constants
const STATIC_LABEL = "static (non-headline)";
const NON_HEADLINE_LABEL = "non-headline model";
// I2: the backend adds this label when the scalar b_res background floor sets g2_op.
const FLOOR_RE = /^structural floor set by b_res/i;
// Seconds per evaluate() per platform class: (single T, 120-point grid). Mirrors
// fsim_studio/jobs.py ETA_TABLE (measured, ideas/01 section 3); used for the pre-run estimate only.
const ETA = { legacy: [0.001, 0.006], edge: [0.31, 20.0], planar: [0.001, 0.31], nanowire: [1.0, 295.0] };
const POOL = 8;
// Display T grids the backend runs for nitride cards when no grid is sent (jobs.NITRIDE_DISPLAY_GRID).
const NITRIDE_GRID_N = { ingan_gan_planar: 31, ingan_gan_nanowire: 16 };
const NANOWIRE_ENVELOPE_REASON = "Envelope runs are not available for the nanowire tier: its evaluate() returns g2_op/T_j curves, not the g2/Tj curves evaluate_envelope bands. Use the committed sweep in Results.";

const BLOCKS = [
  { id: "drive", label: "Drive", prefixes: ["drive."], show: ["drive.V", "drive.I_uA", "drive.mu"] },
  { id: "dot", label: "Dot", prefixes: ["dot.", "ret."], show: ["dot.delta_xx", "dot.gamma_scale"] },
  { id: "cavity", label: "Cavity", prefixes: ["cavity."], show: ["cavity.F_P", "cavity.kappa"], enabled: "cavity.enabled" },
  { id: "filter", label: "Filter", prefixes: ["filter."], show: ["filter.w", "filter.dx"], enabled: "filter.enabled" },
  { id: "aperture", label: "Aperture", prefixes: ["aperture."], show: ["aperture.diameter_um", "aperture.density_cm2"] },
  { id: "detection", label: "Detection", prefixes: [], show: [] },
];
const HEATSINK = { id: "heatsink", label: "Heatsink", prefixes: ["thermal."] };

// ---------------------------------------------------------------- helpers
function getPath(obj, path) {
  return path.split(".").reduce((o, k) => (o == null ? undefined : o[k]), obj);
}
function setPath(obj, path, v) {
  const ks = path.split(".");
  let o = obj;
  for (const k of ks.slice(0, -1)) o = o[k];
  o[ks.at(-1)] = v;
}
const clone = (x) => JSON.parse(JSON.stringify(x));

function platformClass(card, design, mode) {
  if (design?.platform === "ingan_gan_nanowire") return "nanowire";
  if (design?.platform === "ingan_gan_planar") return "planar";
  if (mode === "headline" || design?.emission?.type === "edge" || design?.drive?.finite_pulse) return "edge";
  return "legacy";
}
function evalSeconds(cls, nT) {
  const [t1, tn] = ETA[cls];
  return t1 + (tn - t1) * Math.max(nT - 1, 0) / 119;
}

function cavityLabel(design) {
  const t = design?.emission?.type;
  if (t === "edge") return "Waveguide";
  if (t === "nanowire" || design?.platform === "ingan_gan_nanowire") return "Photonic wire";
  return "Cavity";
}

// ---------------------------------------------------------------- module state
let S = null;

export function mount(el, ctx) {
  S = {
    ctx, el, card: null, info: null, design: null, meta: {}, envDefaults: {}, validation: null,
    ranged: {}, mode: "point", grid: "full", run: null, opOverride: null, runState: "idle", job: null,
    openBlock: null, verdict: null, secondary: "eps", disposers: [], scene: null, sceneSeq: 0,
    readouts: {}, charts: {}, dirty: false, progress: null, lastError: null, t0: 0, timer: null,
    head: null, headJob: null, seed: null,
    // studio-p2a: preset draft, sheet edits (unsaved), response animation
    draft: null, edits: 0, leaveOK: false, leaveAsk: null, anim: null, presetSel: null,
  };
  el.classList.add("ws-designer");
  S.disposers.push(ctx.state.on("theme", (m) => S?.scene?.setTheme?.(m)));
  update(ctx.args);
}

export function unmount() {
  if (!S) return;
  closeAnim({ quiet: true });
  closeLeaveDialog();
  S.ctx.state.set({ unsaved: null });
  S.job?.cancel?.();
  S.headJob?.cancel?.();
  clearInterval(S.timer);
  for (const d of S.disposers) d();
  for (const c of Object.values(S.charts)) c.dispose();
  try { S.scene?.dispose?.(); } catch { /* already gone */ }
  S.ro?.disconnect();
  S.el.classList.remove("ws-designer");
  S = null;
}

export async function update(args) {
  const card = args[0];
  if (!S) return;
  const draft = card === DRAFT_CARD ? S.ctx.state.get("draft") : null;
  if (card === S.card && (card !== DRAFT_CARD || draft === S.draft)) return;
  closeAnim({ quiet: true });
  S.draft = draft;
  S.edits = 0;
  S.leaveOK = false;
  S.leaveAsk = null;
  syncUnsaved();
  S.job?.cancel?.();
  S.job = null;
  S.headJob?.cancel?.();
  S.headJob = null;
  S.head = null;
  S.card = card;
  S.run = null;
  S.opOverride = null;
  S.lastError = null;
  S.openBlock = null;
  S.dirty = false;
  S.ctx.state.set({ run: null });
  renderSkeleton();
  try {
    const saved = !draft && !!S.ctx.cards().find((c) => c.name === card)?.saved;
    const [info, meta, campaigns, verdict] = await Promise.all([
      // a preset draft is not a card (studio-p2a); a Studio-saved card is re-read (it may be re-saved)
      draft ? Promise.resolve({ name: DRAFT_CARD, design: draft.design, meta: draft.meta || {}, platform: draft.design.platform, draft: true })
        : saved ? S.ctx.api.getJSON(`/api/cards/${encodeURIComponent(card)}`) : S.ctx.api.getCard(card),
      S.ctx.api.getMeta(),
      S.ctx.api.listCampaigns().catch(() => []),
      draft || saved
        ? Promise.resolve({ ok: true, v: { match: "none", lines: [], reason: draft ? "An unsaved preset design: no committed VERDICT line speaks for it." : "A design saved in the Studio: the committed VERDICT lines were computed for the shipped cards, not for this one." } })
        : S.ctx.api.cardVerdict(card).then((v) => ({ ok: true, v }), (e) => ({ ok: false, error: e.message })),
    ]);
    if (!S || S.card !== card) return;
    S.info = info;
    S.design = clone(info.design);
    S.meta = info.meta || {};
    S.envDefaults = meta.ENV_DEFAULTS || {};
    S.verdict = verdictFromApi(verdict, campaigns);
    S.seed = takeSeed(card);
    S.validation = await S.ctx.api.validate(S.design).catch(() => null);
    if (!S || S.card !== card) return;
    S.ranged = clone(S.validation?.default_ranged || {});
    S.metaAll = meta.META || {};
    S.presets = meta.presets || {};
    const edge = S.design.emission?.type === "edge";
    const slow = ["edge", "nanowire"].includes(platformClass(card, S.design, "point"));
    S.mode = edge ? "headline" : "point";
    S.grid = slow ? "op" : "full";
    if (S.draft) S.dirty = true;
    if (S.seed?.design) {
      // the row's grid values and model switches, evaluated as a point at the row's T_hs
      S.design = S.seed.design;
      S.mode = "point";
      S.grid = "op";
      S.dirty = true;
    }
    renderAll();
    loadScene();
    run({ auto: true });
  } catch (e) {
    if (!S || S.card !== card) return;
    clear(S.el).append(h("div.ws-pad", errorState({
      title: `Could not open ${card}`,
      message: e.message,
      actions: [{ label: "Retry", fn: () => { S.card = null; update([card]); } }],
    })));
  }
}

// ---------------------------------------------------------------- verdict selection (I6)
// The card -> VERDICT line mapping is the backend's (GET /api/cards/<name>/verdict: exact regime,
// family/orientation, strain bound, rep rate). The UI renders what it returns and never picks a
// line itself: an ambiguous match is shown as ambiguous, with every matching line.
function verdictFromApi(res, campaigns) {
  if (!res.ok) return { state: "error", reason: `verdict mapping unavailable: ${res.error}`, lines: [] };
  const v = res.v || {};
  const camp = campaigns.find((c) => c.id === v.campaign) || null;
  const lines = (v.lines || []).map((l) => normLine(l, camp)).filter(Boolean);
  const match = v.match || (lines.length === 1 ? "exact" : lines.length ? "ambiguous" : "none");
  return { state: match, lines, reason: v.reason || "", campaignId: v.campaign || null, source: camp?.path || (v.campaign ? `out/${v.campaign}` : "") };
}

/** A returned line may be a parsed object or the raw text; map raw text onto the parsed campaign line. */
function normLine(l, camp) {
  if (l && typeof l === "object") return l.fields || l.word ? l : null;
  if (typeof l === "string") {
    const hit = (camp?.verdicts || []).find((x) => x.raw === l || l.includes(x.raw) || x.raw.includes(l));
    return hit || { raw: l, word: (/VERDICT\S*\s+(\S+)/.exec(l) || [])[1] || "", fields: {} };
  }
  return null;
}

/** Results -> "Open card in Designer" leaves state.rerunSeed; take it once, for this card only. */
function takeSeed(card) {
  const seed = S.ctx.state.get("rerunSeed");
  if (!seed || seed.card !== card) return null;
  S.ctx.state.set({ rerunSeed: null });
  const out = { ...seed, design: null };
  if (SEEDABLE.has(seed.campaign) && seed.row) {
    try { out.design = rtEdgeRowToDesign(S.info.design, seed.row); } catch { out.design = null; }
  }
  return out;
}

// ---------------------------------------------------------------- skeleton
function renderSkeleton() {
  clear(S.el);
  S.el.append(h("div.designer.is-loading", { "aria-busy": "true" },
    h("section.board", h("div.board-head", skeleton("title", { lines: 2 })),
      h("div.board-table", h("div.sk-viewport", skeleton("block", { lines: 0 })),
        h("div.sk-chain", Array.from({ length: 6 }, () => skeleton("mount", { lines: 3 }))),
        skeleton("plate", { lines: 1 }))),
    h("section.shelf", skeleton("plaque", { lines: 2 }), skeleton("readouts", { lines: 5 }), skeleton("chart", { lines: 0, height: 300 }))));
}

// ---------------------------------------------------------------- full render
function renderAll() {
  for (const c of Object.values(S.charts)) c.dispose();
  S.charts = {};
  try { S.scene?.dispose?.(); } catch { /* noop */ }
  S.scene = null;
  S.ro?.disconnect();
  clear(S.el);

  const board = renderBoard();
  const shelf = renderShelf();
  S.el.append(h("div.designer", board, shelf));
  S.ro = new ResizeObserver(() => { drawBeam(); syncChainScroll(); });
  S.ro.observe(S.dom.chain);
  S.ro.observe(S.dom.chainScroll);
  requestAnimationFrame(syncChainScroll);
  requestAnimationFrame(drawBeam);
  renderResults();
  renderSeed();
}

// ================================================================ BREADBOARD
function renderBoard() {
  const d = S.design;
  const info = S.info;
  const plat = S.ctx.platformOf(S.card);
  S.dom = {};

  const head = h("header.board-head",
    h("div.bh-titles",
      h("h1.board-title", S.draft ? S.draft.title : S.card.replace(/-design$/, "")),
      h("p.board-desc", S.draft ? draftDesc() : info.design?.provenance?.device || cardDesc())),
    h("div.bh-meta",
      plat ? h("span.chip.chip-plat", plat.label) : null,
      h("span.chip.chip-mono", d.platform),
      S.dom.dirty = h("span.chip.chip-dirty", { hidden: !S.draft, title: "Edits are local to this session until you run" }, S.draft ? "unsaved" : "edited"),
      h("div.bh-keys",
        S.dom.presetBtn = h("button.btn.btn-quiet.bh-presets", { type: "button", "aria-expanded": "false", on: { click: () => togglePresetSheet() } }, icon("design", { size: 16 }), "New from presets"),
        S.dom.saveBtn = h("button.btn.btn-outline.bh-save", { type: "button", hidden: !S.draft && !S.edits, "aria-haspopup": "dialog", on: { click: (e) => saveAsPopover(e.currentTarget) } }, icon("download", { size: 16 }), "Save as card"))));

  // 3D viewport slot
  S.dom.viewport = h("div.viewport-slot", { "aria-label": "3D device model", role: "region" });
  const viewportMount = h("div.viewport-mount",
    h("div.vm-rail", h("span.engrave", "Device model"), h("span.vm-note#vm-note", "true scale; exaggerations are badged")),
    S.dom.viewport,
    h("span.post", { "aria-hidden": "true" }));

  // chain mounts + beam
  S.dom.beamSvg = svgEl("svg", { class: "beam-svg", "aria-hidden": "true" });
  S.dom.mounts = {};
  const mounts = BLOCKS.map((b) => {
    const m = renderMount(b);
    S.dom.mounts[b.id] = m;
    return m;
  });
  S.dom.chain = h("div.chain", { role: "list", "aria-label": "Design chain: drive to detection" }, S.dom.beamSvg, mounts);
  S.dom.beamCap = h("p.beam-cap");
  // narrow screens: the chain scrolls sideways; edge fades + a thin track say so (R1)
  S.dom.chainScroll = h("div.chain-scroll", { on: { scroll: syncChainScroll } }, S.dom.chain);
  S.dom.chainTrack = h("div.chain-track", { "aria-hidden": "true" }, h("span.chain-thumb"));

  // heatsink plate (signature)
  S.dom.heatsink = renderHeatsink();

  // sheet host
  S.dom.sheet = h("div.sheet-host");

  // run bar
  S.dom.runbar = renderRunbar();

  // studio-p2a hosts: preset builder sheet (top of the table), animation sheet (under the heatsink)
  S.dom.presetHost = h("div.preset-host");
  S.dom.animHost = h("div.anim-host");

  const table = h("div.board-table",
    S.dom.presetHost,
    viewportMount,
    h("div.chain-wrap", S.dom.chainScroll, S.dom.chainTrack, S.dom.beamCap),
    S.dom.heatsink,
    S.dom.animHost,
    S.dom.sheet);

  return h("section.board", { "aria-label": "Breadboard" }, head, table, S.dom.runbar);
}

function cardDesc() {
  const c = S.ctx.cards().find((x) => x.name === S.card);
  return c?.device || `Design card cards/${S.card}.yaml (${S.design?.platform || "legacy"} platform); no device description recorded.`;
}

function blockMetaPaths(b) {
  return Object.keys(S.meta).filter((p) => b.prefixes.some((pre) => p.startsWith(pre)) && getPath(S.design, p) !== undefined);
}

function renderMount(b) {
  const label = b.id === "cavity" ? cavityLabel(S.design) : b.label;
  const paths = blockMetaPaths(b);
  const wide = widest(paths.map((p) => S.meta[p]?.tag));
  const enabled = b.enabled ? getPath(S.design, b.enabled) !== false : true;
  const vals = h("div.mount-vals");
  if (b.id === "detection") {
    const em = S.design.emission || {};
    vals.append(h("div.mv-text", em.type && em.type !== "none" ? `${em.type.replace(/_/g, " ")} collection` : "free-space collection"));
  } else {
    for (const p of b.show.filter((x) => S.meta[x] && getPath(S.design, x) !== undefined).slice(0, 2)) {
      const v = getPath(S.design, p);
      const m = S.meta[p];
      vals.append(h("div.mv", { title: `${p}: ${m.source || ""}` }, h("span.mv-k", m.label || p.split(".").pop()),
        h("span.mv-line", h("span.mv-v", v === null ? "auto" : typeof v === "number" ? fmt(v) : String(v)),
          h("span.mv-u", unitText(m.unit)), tagChip(m.tag, { size: "sm" }))));
    }
  }
  const btn = h("button.mount", {
    type: "button",
    class: `${!enabled ? "is-bypassed" : ""} ${S.openBlock === b.id ? "is-open" : ""}`,
    dataset: { block: b.id, rim: wide || "" },
    "aria-expanded": S.openBlock === b.id ? "true" : "false",
    "aria-label": `${label} block${!enabled ? ", bypassed" : ""}${wide ? `, widest input tag ${wide}` : ""}. Open parameter sheet`,
    on: { click: () => toggleBlock(b) },
  },
  h("span.mount-head", h("span.engrave", label), wide ? tagChip(wide, { size: "sm", title: "widest tag of this block's inputs" }) : null),
  // clear aperture: the beam passes through every optic mount at this height
  h("span.mount-aperture", { "aria-hidden": "true" }),
  !enabled ? h("span.mount-bypass", "bypassed") : vals,
  h("span.mount-foot", `${paths.length} field${paths.length === 1 ? "" : "s"}`));
  return h("div.mount-cell", { role: "listitem" }, btn, h("span.post", { "aria-hidden": "true" }));
}

function refreshMounts() {
  for (const b of BLOCKS) {
    const fresh = renderMount(b);
    S.dom.mounts[b.id].replaceWith(fresh);
    S.dom.mounts[b.id] = fresh;
  }
  S.dom.chain.prepend(S.dom.beamSvg);
  requestAnimationFrame(drawBeam);
}

function toggleBlock(b) {
  S.openBlock = S.openBlock === b.id ? null : b.id;
  refreshMounts();
  renderSheet();
}

// ---------------------------------------------------------------- beam
// ONE continuous beam: it leaves the device's emission point (the Dot mount's aperture), passes
// through the clear aperture of every optic mount and stops on Detection. Drawn behind the mounts
// (they are windowed at the beam height), >= 4px, visible at rest, dimming with log10 collected flux.
function beamAlpha(s) {
  const flux = isNum(s?.collected_flux_delivered_s) ? s.collected_flux_delivered_s : s?.collected_flux_pulsed_s;
  if (!(isNum(flux) && flux > 0)) return { a: 0.45, flux: null };
  const k = Math.max(0, Math.min(1, (Math.log10(flux) - 2) / 5));
  return { a: 0.45 + 0.55 * k, flux };
}

function beamIntensity() {
  const s = (S.opOverride?.result || S.run)?.scalars || {};
  const { a, flux } = beamAlpha(s);
  if (flux != null) {
    return { a, flux, text: `Beam intensity follows log₁₀ collected flux: ${fmt(flux)} s⁻¹` };
  }
  return { a: 0.45, flux: null, text: S.run ? "Beam at rest level: this path computes no collected flux" : "Beam at rest level until the first run" };
}

/** Scroll affordance for the chain row: data-overflow / data-edge drive the CSS edge fades; the
 *  track thumb mirrors the visible fraction (headless browsers hide native scrollbars). */
function syncChainScroll() {
  const sc = S?.dom?.chainScroll;
  if (!sc) return;
  const wrap = sc.parentElement;
  const over = sc.scrollWidth > sc.clientWidth + 1;
  wrap.dataset.overflow = over ? "1" : "0";
  if (!over) { delete wrap.dataset.edge; return; }
  const max = sc.scrollWidth - sc.clientWidth;
  wrap.dataset.edge = sc.scrollLeft <= 1 ? "start" : sc.scrollLeft >= max - 1 ? "end" : "mid";
  const thumb = S.dom.chainTrack.firstChild;
  thumb.style.width = `${(100 * sc.clientWidth) / sc.scrollWidth}%`;
  thumb.style.left = `${(100 * sc.scrollLeft) / sc.scrollWidth}%`;
}

// Mounts snap to the 25px M6 hole pitch: the mount pitch is a multiple of 25px and each mount's
// width is chosen so its centre (and its post) lands on a hole centre (12.5 + 25n from the page
// origin, where the body's hole grid is anchored).
const HOLE = 25;
function snapChain() {
  const chain = S.dom.chain;
  const reset = () => { chain.style.gridTemplateColumns = ""; chain.style.columnGap = ""; chain.dataset.snap = "0"; };
  if (window.innerWidth <= 760) { reset(); return; }
  const left = chain.getBoundingClientRect().left + parseFloat(getComputedStyle(chain).paddingLeft || "0");
  const avail = chain.clientWidth - parseFloat(getComputedStyle(chain).paddingLeft || "0") - parseFloat(getComputedStyle(chain).paddingRight || "0");
  const r = (((HOLE - 2 * left) % 50) + 50) % 50; // mount width W must satisfy W = r (mod 50)
  let P = (Math.floor(avail / 6 / HOLE) + 1) * HOLE;
  for (; P >= 75; P -= HOLE) {
    const W = r + 50 * Math.floor((P - 10 - r) / 50);
    if (W >= 60 && 5 * P + W <= avail) {
      const cols = `repeat(6, ${W}px)`;
      if (chain.style.gridTemplateColumns !== cols) chain.style.gridTemplateColumns = cols;
      if (chain.style.columnGap !== `${P - W}px`) chain.style.columnGap = `${P - W}px`;
      chain.dataset.snap = String(P);
      return;
    }
  }
  reset();
}

function drawBeam() {
  if (!S?.dom?.chain) return;
  snapChain();
  const svg = S.dom.beamSvg;
  const box = S.dom.chain.getBoundingClientRect();
  if (!box.width) return;
  svg.setAttribute("viewBox", `0 0 ${box.width} ${box.height}`);
  svg.setAttribute("width", box.width);
  svg.setAttribute("height", box.height);
  clear(svg);
  const at = (id) => {
    const m = S.dom.mounts[id].querySelector(".mount");
    const r = m.getBoundingClientRect();
    const ap = m.querySelector(".mount-aperture")?.getBoundingClientRect();
    return { x: r.left - box.left + r.width / 2, y: ap ? ap.top - box.top + ap.height / 2 : r.top - box.top + r.height / 2, l: r.left - box.left, r: r.right - box.left };
  };
  const drive = at("drive"), dot = at("dot"), det = at("detection");
  const y = dot.y;
  // electrical lead drive -> dot (ink, not beam)
  svg.append(svgEl("path", { d: `M ${drive.x} ${y} L ${dot.x} ${y}`, class: "lead" }));
  const { a, text } = beamIntensity();
  const running = S.runState === "running";
  const x1 = det.x;
  const pathD = `M ${dot.x} ${y} L ${x1} ${y}`;
  svg.append(svgEl("path", { d: pathD, class: "beam-halo", style: `opacity:${(a * 0.4).toFixed(3)}` }));
  svg.append(svgEl("path", { d: pathD, class: `beam${running && !reducedMotion() ? " is-running" : ""}`, style: `opacity:${a.toFixed(3)}`, "data-alpha": a.toFixed(3) }));
  // emission point glow (device) and the detector stop
  svg.append(svgEl("circle", { cx: dot.x, cy: y, r: 7, class: "beam-src", style: `opacity:${a.toFixed(3)}` }));
  svg.append(svgEl("path", { d: `M ${x1} ${y - 9} L ${x1} ${y + 9}`, class: "beam-stop" }));
  S.dom.beamCap.textContent = text;
}

// ---------------------------------------------------------------- heatsink (signature interaction)
function renderHeatsink() {
  const p = "thermal.T_hs";
  const m = S.meta[p] || { tag: "E", unit: "K", band: [4, 350], label: "T_hs" };
  const T = getPath(S.design, p);
  const band = m.band || [4, 350];
  const valEl = h("output.hs-val", { for: "hs-slider" }, fmt(T));
  const slider = h("input.hs-slider#hs-slider", {
    type: "range", min: band[0], max: band[1], step: 1, value: T,
    "aria-label": "Heatsink temperature T_hs, kelvin. Dragging re-runs the operating point",
    on: {
      input: (e) => { valEl.textContent = fmt(Number(e.target.value)); signatureRun(Number(e.target.value)); },
    },
  });
  S.dom.hsSlider = slider;
  S.dom.hsVal = valEl;
  const open = h("button.btn.btn-quiet.hs-open", { type: "button", "aria-expanded": "false", on: { click: () => toggleBlock(HEATSINK) } }, "Thermal stack", icon("chevronDown", { size: 16 }));
  const anim = h("button.btn.btn-quiet.hs-anim", { type: "button", "aria-expanded": "false", title: "Animate the response over T_hs: compute the frames once, then play them back", on: { click: () => (S.anim?.param === p ? closeAnim() : openAnim(p)) } }, icon("run", { size: 16 }), "Animate");
  S.dom.hsAnim = anim;
  return h("div.heatsink", { role: "group", "aria-label": "Heatsink plate" },
    h("div.hs-label", h("span.engrave", "Heatsink"), tagChip(m.tag, { size: "sm", title: m.source })),
    h("div.hs-control",
      h("span.hs-end", `${fmt(band[0])} K`), slider, h("span.hs-end", `${fmt(band[1])} K`)),
    h("div.hs-read", valEl, h("span.hs-unit", "K")),
    h("div.hs-keys", anim, open));
}

const signatureRun = debounce((T) => {
  if (!S) return;
  setPath(S.design, "thermal.T_hs", T);
  markDirty();
  const mode = S.mode;
  const body = { card: S.card, design: S.design, mode, T_grid: [T] };
  if (mode === "envelope") body.ranged = S.ranged;
  S.sigJob?.cancel?.();
  S.dom.heatsink.classList.add("is-settling");
  const job = S.ctx.api.runJob(body);
  S.sigJob = job;
  job.promise.then((r) => {
    if (!S || S.sigJob !== job) return;
    S.opOverride = { T, result: r.result };
    S.dom.heatsink.classList.remove("is-settling");
    if ((r.result?.labels || []).includes(STATIC_LABEL) && r.result?.mode === "point" && S.design.emission?.type === "edge") runHeadlineBeside(T);
    renderResults({ signature: true });
    loadScene(T);
  }).catch((e) => {
    if (!S || S.sigJob !== job) return;
    S.dom.heatsink.classList.remove("is-settling");
    if (e.body?.error?.kind !== "cancelled") toast(`Operating point at ${fmt(T)} K: ${e.message}`, { kind: "error" });
  });
}, 280);

function markDirty() {
  S.dirty = true;
  if (S.dom?.dirty) S.dom.dirty.hidden = false;
}

// ---------------------------------------------------------------- block editor sheet
function renderSheet() {
  const host = S.dom.sheet;
  clear(host);
  const hsBtn = S.dom.heatsink.querySelector(".hs-open");
  hsBtn?.setAttribute("aria-expanded", S.openBlock === "heatsink" ? "true" : "false");
  if (!S.openBlock) { host.classList.remove("is-open"); return; }
  const b = S.openBlock === "heatsink" ? HEATSINK : BLOCKS.find((x) => x.id === S.openBlock);
  const label = b.id === "cavity" ? cavityLabel(S.design) : b.label;
  const paths = blockMetaPaths(b);
  const body = h("div.sheet-grid");
  if (b.id === "detection") {
    const em = S.design.emission || {};
    const prov = (S.run?.provenance || {}).emission;
    body.append(h("div.sheet-note",
      h("p", `Collection path: ${em.type && em.type !== "none" ? em.type.replace(/_/g, " ") : "free space (no edge or wire collection model)"}.`),
      h("p", "Detection and collection fields carry no META entry yet, so their numbers stay in the card YAML; the run reports their provenance below."),
      prov ? h("p.sheet-prov", tagChip(prov.tag), h("span", prov.note || "")) : h("p.sheet-prov", "Run the design to see the collection provenance.")));
  } else if (!paths.length) {
    body.append(h("p.sheet-note", "This block has no editable META fields on this card."));
  } else {
    for (const p of paths) {
      const row = fieldRow({
        path: p, meta: S.meta[p], value: getPath(S.design, p),
        ranged: S.ranged[p] || null, envDefault: S.envDefaults[p],
        onChange: (v) => {
          setPath(S.design, p, v);
          markDirty();
          noteEdit();
          if (p === "thermal.T_hs") { S.dom.hsSlider.value = v; S.dom.hsVal.textContent = fmt(v); }
          if (p === "cavity.enabled" || p === "filter.enabled" || b.show?.includes(p)) refreshMountsKeepSheet();
          inspect(p);
          S.ctx.inspector.field({ path: p, meta: S.meta[p], value: v, inBand: inBand(p, v) }, { quiet: true });
        },
        onRange: (r) => {
          if (r) S.ranged[p] = r; else delete S.ranged[p];
          syncModeSubs();
        },
        onInspect: (path, o) => inspect(path, o),
      });
      addAnimGlyph(row, p);
      body.append(row);
    }
  }
  const nRanged = paths.filter((p) => S.ranged[p]).length;
  host.classList.add("is-open");
  host.append(h("section.sheet", { "aria-label": `${label} parameters` },
    h("header.sheet-head",
      h("div", h("h2.engrave.sheet-title", label), h("p.sheet-sub", `${paths.length} fields from design_meta.META${nRanged ? `; ${nRanged} ranged for envelope runs` : ""}. Edits apply to the next run.`)),
      h("button.icon-btn", { type: "button", "aria-label": `Close ${label} sheet`, on: { click: () => { S.openBlock = null; refreshMounts(); renderSheet(); S.dom.mounts[b.id]?.querySelector(".mount")?.focus(); } } }, icon("close", { size: 18 }))),
    body));
  if (!reducedMotion()) host.querySelector(".sheet").classList.add("mount-in");
  requestAnimationFrame(() => host.scrollIntoView({ block: "nearest", behavior: reducedMotion() ? "auto" : "smooth" }));
}

function refreshMountsKeepSheet() { refreshMounts(); }

function inBand(p, v) {
  const band = S.meta[p]?.band;
  if (!isNum(v) || !Array.isArray(band) || !isNum(band[0])) return null;
  return v >= band[0] && v <= band[1];
}

function inspect(path, { quiet = false } = {}) {
  const v = getPath(S.design, path);
  S.ctx.inspector.field({ path, meta: S.meta[path], value: v, inBand: inBand(path, v) }, { quiet });
}

// ---------------------------------------------------------------- run bar
function renderRunbar() {
  const edge = S.design.emission?.type === "edge";
  S.seg = segmented({
    label: "Run mode",
    value: S.mode,
    options: [
      { value: "point", label: "Point" },
      { value: "headline", label: "Headline", disabled: !edge, title: edge ? "Finite-pulse headline model (RT edge)" : "Headline model applies to RT edge cards only" },
      { value: "envelope", label: "Envelope" },
    ],
    onChange: (v) => { S.mode = v; if (v === "headline") S.grid = "op"; syncModeSubs(); },
  });
  const slow = ["edge", "nanowire"].includes(platformClass(S.card, S.design, "point"));
  S.gridSeg = slow ? segmented({
    label: "Temperature grid",
    value: S.grid,
    options: [{ value: "op", label: "T_hs only" }, { value: "full", label: "Full curve" }],
    onChange: (v) => { S.grid = v; syncModeSubs(); },
  }) : null;

  S.dom.runBtn = h("button.run-key", { type: "button", on: { click: () => run() } }, icon("run", { size: 18 }), h("span", "RUN"));
  S.dom.cancelBtn = h("button.btn.btn-outline.cancel-btn", { type: "button", hidden: true, on: { click: cancelRun } }, icon("stop", { size: 16 }), "Cancel");
  S.dom.progress = h("div.progress", { role: "progressbar", "aria-label": "Run progress", "aria-valuemin": "0", "aria-valuemax": "100", hidden: true }, h("span.progress-fill"));
  S.dom.runStatus = h("p.run-status", { "aria-live": "polite" });
  S.dom.eta = h("p.run-eta");
  const bar = h("footer.runbar", { "aria-label": "Run controls" },
    h("div.rb-left", h("div.rb-row", S.seg.el, S.gridSeg ? S.gridSeg.el : null), S.dom.eta),
    h("div.rb-right", h("div.rb-status", S.dom.progress, S.dom.runStatus), S.dom.cancelBtn, S.dom.runBtn));
  queueMicrotask(syncModeSubs);
  return bar;
}

function runPlan() {
  const cls = platformClass(S.card, S.design, S.mode);
  const nT = S.grid === "op" ? 1 : (NITRIDE_GRID_N[S.design?.platform] || 120);
  if (S.mode === "envelope") {
    const ks = Object.values(S.ranged);
    const n = ks.reduce((a, v) => a * v.length, 1);
    const nEval = n + 1 + ks.reduce((a, v) => a + n / v.length, 0) + 1;
    return { cls, nT, n, eta: evalSeconds(cls, nT) * (nEval / POOL + 2) };
  }
  return { cls, nT, n: 1, eta: evalSeconds(cls, nT) + (["edge", "nanowire"].includes(cls) && nT > 1 ? evalSeconds(cls, 1) : 0) };
}

function syncModeSubs() {
  if (!S?.seg) return;
  S.seg.select(S.mode);
  S.gridSeg?.select(S.grid);
  const k = Object.keys(S.ranged).length;
  const plan = runPlan();
  const nEnv = Object.values(S.ranged).reduce((a, v) => a * v.length, 1);
  const envSub = k ? `${nEnv} samples` : "no ranged inputs";
  if (S.design?.platform === "ingan_gan_nanowire") {
    // N1: the backend refuses envelope runs for this tier (400 kind "unsupported"); say why here.
    if (S.mode === "envelope") S.mode = "point";
    S.seg.select(S.mode);
    S.seg.setSub("envelope", "not for nanowire");
    S.seg.setDisabled("envelope", true, NANOWIRE_ENVELOPE_REASON);
  } else {
    S.seg.setSub("envelope", envSub);
    S.seg.setDisabled("envelope", k === 0, k === 0 ? "Range at least one [A] or [E] input in a block sheet" : "Sweep every ranged input at its extremes");
  }
  S.dom.eta.textContent = S.mode === "envelope"
    ? `${k} ranged input${k === 1 ? "" : "s"}, ${plan.n} samples at extremes; estimated ${fmtSeconds(plan.eta)} on this machine`
    : `${S.grid === "op" ? "Operating point only" : `${plan.nT}-point T_hs grid`}; estimated ${fmtSeconds(plan.eta)}`;
}

// ---------------------------------------------------------------- run
function run({ auto = false } = {}) {
  if (!S) return;
  if (S.mode === "envelope" && !Object.keys(S.ranged).length) { toast("Range at least one input before an envelope run."); return; }
  S.job?.cancel?.();
  const T = getPath(S.design, "thermal.T_hs");
  const body = { card: S.card, design: S.design, mode: S.mode };
  const runT = T;
  if (S.grid === "op") body.T_grid = [T];
  if (S.mode === "envelope") body.ranged = S.ranged;
  const plan = runPlan();
  S.lastError = null;
  S.opOverride = null;
  setRunState("running", { eta: plan.eta, n: plan.n });
  const job = S.ctx.api.runJob(body, {
    onProgress: (k, n) => { if (S?.job === job) setProgress(k, n); },
    onPartial: (res) => {
      if (S?.job !== job) return;
      // progressive result: fill the readouts with the single-T values, draw a phosphor trace
      S.partial = res;
      renderReadouts(res, { partial: true });
      const c = themeColors();
      if (res.curves?.T_hs?.length > 1) S.charts.hero?.phosphor(res.curves.T_hs, res.curves.g2, c.series[0]);
    },
  });
  S.job = job;
  job.promise.then((r) => {
    if (!S || S.job !== job) return;
    S.job = null;
    S.run = r.result;
    S.runT = runT;
    S.partial = null;
    S.ctx.state.set({ run: r.result });
    S.head = null;
    S.headJob?.cancel?.();
    S.headJob = null;
    if (r.result?.mode === "point" && (r.result.labels || []).includes(STATIC_LABEL) && S.design.emission?.type === "edge") runHeadlineBeside(runT);
    setRunState("done", { elapsed: r.elapsed_s, cached: r.cached });
    renderResults({ fresh: true });
    renderSeed();
    if (!auto) S.charts.hero?.el.scrollIntoView?.({ block: "nearest" });
  }).catch((e) => {
    if (!S || S.job !== job) return;
    S.job = null;
    // a 400 from POST /api/run carries {error: "<text>", kind, message} at the top level
    const b = e.body;
    const err = b?.error && typeof b.error === "object" ? b.error : b?.kind ? { kind: b.kind, message: b.message || b.error } : { message: e.message };
    if (err.kind === "cancelled") { setRunState("cancelled"); return; }
    S.lastError = err;
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
  S.dom.runStatus.textContent = `Sample ${k} of ${n}`;
}

function setRunState(st, info = {}) {
  S.runState = st;
  S.ctx.state.set({ runState: st });
  const running = st === "running";
  S.dom.runBtn.disabled = running;
  S.dom.runBtn.setAttribute("aria-busy", running ? "true" : "false");
  S.dom.cancelBtn.hidden = !running;
  S.dom.progress.hidden = !running;
  clearInterval(S.timer);
  if (running) {
    S.t0 = performance.now();
    S.dom.progress.classList.add("is-indeterminate");
    S.dom.progress.firstChild.style.transform = "";
    S.dom.progress.removeAttribute("aria-valuenow");
    const tick = () => {
      if (!S) return;
      const el = (performance.now() - S.t0) / 1000;
      if (!S.dom.runStatus.textContent.startsWith("Sample")) {
        const over = info.eta > 0 && el > info.eta;
        S.dom.runStatus.textContent = over
          ? `Running, ${fmtSeconds(el)}; past the ${fmtSeconds(info.eta)} estimate`
          : `Running, ${fmtSeconds(el)} of about ${fmtSeconds(info.eta)}`;
        // time-based fill (estimate), dimmed until the server reports sample progress; past the
        // estimate it keeps creeping toward (never reaching) full instead of parking at one value
        const f = !(info.eta > 0) ? 0.5 : over ? 0.9 + 0.09 * (1 - Math.exp(-(el - info.eta) / info.eta)) : 0.9 * el / info.eta;
        S.dom.progress.firstChild.style.transform = `scaleX(${f})`;
      }
    };
    S.dom.runStatus.textContent = "Submitting";
    tick();
    S.timer = setInterval(tick, 250);
    renderReadouts(null, { pending: true });
  } else if (st === "done") {
    S.dom.runStatus.textContent = info.cached ? `Done from cache in ${fmtSeconds(info.elapsed)}` : `Done in ${fmtSeconds(info.elapsed)}`;
  } else if (st === "cancelled") {
    S.dom.runStatus.textContent = "Cancelled; the last finished result stays on the shelf";
    renderResults();
  } else if (st === "error") {
    S.dom.runStatus.textContent = "Run failed; see the shelf";
  } else {
    S.dom.runStatus.textContent = "";
  }
  renderStatus();
  drawBeam();
}

// ================================================================ SHELF
function renderShelf() {
  S.dom.plaqueHost = h("div.plaque-host");
  S.dom.status = h("div.shelf-status");
  S.dom.pin = renderPin();
  S.dom.alert = h("div.shelf-alert");
  const ro = (key, label) => (S.readouts[key] = createReadout({ key, label, onInspect: inspectKpi })).el;
  S.dom.readouts = h("div.readouts", { role: "group", "aria-label": "Key readouts" },
    ro("g2_op", "g²(0)"),
    ro("g2_head", "g²(0), headline model"),
    ro("flux_delivered", "Collected flux, delivered"),
    ro("flux", "Collected flux"),
    ro("brightness", "Brightness per pulse"),
    ro("T_c", "Ceiling T_c"),
    ro("T_j_op", "Junction T_j"));

  S.charts.hero = chartCard({
    title: "g²(0) against heatsink temperature",
    build: () => (S?.run ? buildG2T([{ result: S.run, opT: heroOpT(), opResult: S.opOverride?.result, name: "g²(0)" }]) : null),
    table: () => (S?.run ? g2TableRows(S.run) : null),
    filename: () => `${S?.run?.run_id || "fsim"}-g2T`,
  });
  S.charts.hero.el.classList.add("chart-hero");

  const tabs = Object.entries(SECONDARY).map(([k, v]) => h("button.sec-tab", {
    type: "button", role: "tab", dataset: { key: k }, "aria-selected": String(S.secondary === k),
    on: { click: () => { S.secondary = k; syncSecTabs(); S.charts.sec.setTitle(v.label); S.charts.sec.render(); } },
  }, { eps: "ε", rho2: "ρ²", dTJ: "ΔT_J", gamma: "Γ" }[k]));
  S.dom.secTabs = h("div.sec-tabs", { role: "tablist", "aria-label": "Secondary quantity" }, tabs);
  S.charts.sec = chartCard({
    title: SECONDARY[S.secondary].label,
    height: 200,
    build: () => (S?.run ? buildSecondary(S.secondary, [{ result: S.run }]) : null),
    table: () => secTable(),
    filename: () => `${S?.run?.run_id || "fsim"}-${S?.secondary}`,
  });
  S.charts.sec.el.classList.add("chart-sec");

  S.dom.seedHost = h("div.seed-host");
  S.dom.details = h("details.drawer", h("summary", h("span.engrave", "All scalars"), h("span.drawer-count")), h("div.drawer-body"));
  return h("section.shelf", { "aria-label": "Instrument shelf" },
    S.dom.plaqueHost,
    h("div.shelf-bar", S.dom.status, S.dom.pin),
    S.dom.alert,
    S.dom.seedHost,
    S.dom.readouts,
    S.charts.hero.el,
    h("div.sec-wrap", S.dom.secTabs, S.charts.sec.el),
    S.dom.details);
}

function syncSecTabs() {
  for (const t of S.dom.secTabs.children) t.setAttribute("aria-selected", String(t.dataset.key === S.secondary));
}

function secTable() {
  const r = S?.run;
  if (!r?.curves?.T_hs) return null;
  const spec = SECONDARY[S.secondary];
  const x = r.curves.T_hs;
  const sub = (arr, i) => (arr?.[i] == null ? null : spec.minusX ? arr[i] - x[i] : arr[i]);
  const b = r.bands?.[spec.band];
  if (b) return { columns: [{ key: "T", label: "T_hs (K)" }, { key: "lo", label: "lo" }, { key: "mid", label: "mid" }, { key: "hi", label: "hi" }], rows: x.map((t, i) => ({ T: t, lo: sub(b.lo, i), mid: sub(b.mid, i), hi: sub(b.hi, i) })) };
  return { columns: [{ key: "T", label: "T_hs (K)" }, { key: "v", label: spec.label }], rows: x.map((t, i) => ({ T: t, v: sub(r.curves[spec.curve], i) })) };
}

function heroOpT() {
  return S.opOverride ? S.opOverride.T : S.runT ?? null;
}

function renderPin() {
  const slots = () => S.ctx.state.get("slots");
  const pinTo = (slot) => {
    if (!S.run) { toast("Run the design before pinning it."); return; }
    const all = { ...slots() };
    all[slot] = {
      slot, card: S.card, run_id: S.run.run_id, mode: S.run.mode, result: S.run,
      opT: heroOpT(), design: S.design, pinned: new Date().toISOString(),
    };
    S.ctx.state.set({ slots: all });
    toast(`Pinned ${S.run.run_id} to slot ${slot}.`, { kind: "ok", action: { label: "Open Compare", fn: () => S.ctx.navigate("#/compare") } });
    closePopover();
    syncPin();
  };
  const main = h("button.btn.btn-outline.pin-main", { type: "button", on: { click: () => { const s = slots(); pinTo(["A", "B", "C"].find((k) => !s[k]) || "A"); } } }, icon("pin", { size: 16 }), h("span.pin-txt", "Pin to A"));
  const more = h("button.btn.btn-outline.pin-more", {
    type: "button", "aria-haspopup": "dialog", "aria-label": "Choose a comparison slot",
    on: {
      click: (e) => {
        const s = slots();
        popover(e.currentTarget, h("div.pin-menu", h("h2.why-title", "Pin this run to"),
          ["A", "B", "C"].map((k) => h("button.pin-opt", { type: "button", dataset: { slot: k }, on: { click: () => pinTo(k) } },
            h("span.slot-badge", { dataset: { slot: k } }, k),
            h("span.pin-opt-txt", s[k] ? `replace ${s[k].run_id}` : "empty slot")))), { label: "Pin to slot" });
      },
    },
  }, icon("chevronDown", { size: 16 }));
  S.syncPin = syncPin;
  function syncPin() {
    const s = slots();
    const next = ["A", "B", "C"].find((k) => !s[k]) || "A";
    main.querySelector(".pin-txt").textContent = `Pin to ${next}`;
    main.dataset.slot = next;
  }
  queueMicrotask(syncPin);
  return h("div.pin-split", { role: "group", "aria-label": "Pin result for comparison" }, main, more);
}

// ---------------------------------------------------------------- results -> shelf
function renderResults({ fresh = false, signature = false } = {}) {
  if (!S?.dom) return;
  // plaque
  clear(S.dom.plaqueHost).append(renderPlaque());
  renderStatus();
  const r = S.run;

  // error / alerts
  clear(S.dom.alert);
  if (S.lastError) S.dom.alert.append(renderError(S.lastError));
  else if (r) {
    const s = r.scalars || {};
    const alerts = [];
    if (s.runaway === true) alerts.push("Thermal runaway at the operating point: the junction temperature did not settle.");
    if (isNum(s.eps_op) && s.eps_op > 1) alerts.push("XX leakage ε exceeds 1 at the operating point.");
    if (s.g2_op_valid === false && s.g2_op_invalid_reason) alerts.push(`g²(0) not valid: ${s.g2_op_invalid_reason}.`);
    if (alerts.length) S.dom.alert.append(h("div.alert-list", { role: "note" }, alerts.map((a) => h("p", icon("warn", { size: 16 }), a))));
  }

  renderReadouts(S.opOverride?.result || r, { markChange: fresh || signature, envelopeRun: r });
  drawBeam();
  // an open animation keeps the shelf on the playhead
  if (S.anim?.tl) queueMicrotask(() => { if (S?.anim?.tl) animApply(S.anim.tl.pos, S.anim.tl.where()); });

  // charts
  if (r) {
    const sub = `${r.run_id} · ${r.mode === "envelope" ? "envelope bands, never averaged" : r.mode}${S.opOverride ? ` · operating point at ${fmt(S.opOverride.T)} K` : ""}`;
    S.charts.hero.setTitle(r.mode === "headline" ? "g²(0) against heatsink temperature, headline model" : "g²(0) against heatsink temperature", sub);
    S.charts.hero.setTag(r.tag_chain);
    S.charts.sec.setTag(r.tag_chain);
  } else {
    S.charts.hero.setTag(null);
  }
  const plan = runPlan();
  const noGrid = r && !(r.curves?.T_hs?.length);
  const emptyNode = noGrid && S.runState !== "running"
    ? emptyState({ title: "No temperature axis returned", body: `${r.run_id} carries no T_hs grid (curves.T_hs is empty), so there is nothing to plot against temperature; the readouts above are the operating point.`, iconName: "chart" })
    : S.runState === "running"
    ? h("div.chart-loading", skeleton("chart", { lines: 0 }), h("p", "Evaluating; the curve draws when the job finishes."))
    : emptyState({ title: "Not run yet", body: `Run this design to draw g²(0) against T_hs (estimated ${fmtSeconds(plan.eta)}).`, action: { label: "Run now", fn: () => run() }, iconName: "chart" });
  S.charts.hero.setEmpty(emptyNode);
  S.charts.sec.setEmpty(r && r.curves?.T_hs?.length === 1
    ? emptyState({ title: "Single temperature point", body: "Secondary curves need the full T_hs grid. Switch the grid to Full curve and run.", iconName: "chart" })
    : h("p.chart-sub.pad", r ? "" : "Runs fill this plot."));
  S.charts.hero.render();
  S.charts.sec.render();
  renderDetails();
}

function renderPlaque() {
  const v = S.verdict || { state: "none", lines: [] };
  if (v.state === "exact" && v.lines.length === 1) {
    return verdictPlaque(v.lines[0], { source: v.source, note: v.reason });
  }
  if (v.state === "ambiguous" || v.lines.length > 1) {
    return h("section.plaque.plaque-info.plaque-ambiguous", { "aria-label": "Campaign verdict", dataset: { status: "ambiguous" } },
      h("div.plaque-top", h("div.plaque-word", icon("info", { size: 18 }), h("span", `AMBIGUOUS: ${v.lines.length} LINES MATCH`))),
      h("p.plaque-reason", v.reason || "The card matches more than one VERDICT line; none is picked for it."),
      h("ul.plaque-lines", v.lines.map((l) => h("li.plaque-line",
        verdictBadge(l, { size: "sm" }),
        h("span.pl-scope", `line ${l.line ?? "?"}${lineScope(l) ? ` · ${lineScope(l)}` : ""}`)))),
      h("p.plaque-source", `${v.source}`));
  }
  return verdictPlaque(null, { note: v.reason || "No VERDICT line in out/ matches this card." });
}

function renderStatus() {
  if (!S?.dom?.status) return;
  const r = S.run;
  const labels = r?.labels || [];
  append(clear(S.dom.status), [
    h("span.status-word", { dataset: { state: S.runState } }, statusIcon(), statusWord()),
    r ? h("span.status-mode", r.mode === "headline" ? "headline model" : r.mode === "envelope" ? `envelope, ${r.n_samples} samples` : "point run") : null,
    S.opOverride ? h("span.status-op", `operating point re-run at ${fmt(S.opOverride.T)} K`) : null,
    labels.filter((l) => l !== "pre-retention" && !l.startsWith("all ranges")).map((l) => h("span.label-chip", l))]);
}

function statusWord() {
  return { idle: "Ready", running: "Running", done: "Result on shelf", error: "Run failed", cancelled: "Cancelled" }[S.runState] || S.runState;
}
function statusIcon() {
  return icon({ idle: "minus", running: "retry", done: "check", error: "warn", cancelled: "stop" }[S.runState] || "minus", { size: 16 });
}

function renderError(err) {
  const actions = [];
  if (err.kind === "f8_domain" && isNum(err.f8_floor)) {
    actions.push({
      label: `Narrow drive.mu to ≥ ${fmt(err.f8_floor)}`, primary: true,
      fn: () => {
        const cur = S.ranged["drive.mu"] || [err.f8_floor, 1];
        S.ranged["drive.mu"] = [Math.max(cur[0], err.f8_floor), Math.max(cur[1], err.f8_floor)];
        syncModeSubs();
        run();
      },
    });
    actions.push({ label: "Run a point instead", fn: () => { S.mode = "point"; syncModeSubs(); run(); } });
    actions.push({ label: "Open the Drive sheet", fn: () => { S.openBlock = "drive"; refreshMounts(); renderSheet(); } });
  } else if (err.kind === "envelope_all_dropped" || err.kind === "unsupported") {
    // honest message only: no drive.mu remedy (cycle_loading cards never read drive.mu)
    actions.push({ label: "Run a point instead", primary: true, fn: () => { S.mode = "point"; syncModeSubs(); run(); } });
  } else if (err.kind === "offline") {
    actions.push({ label: "Retry", primary: true, fn: () => run() });
  } else {
    actions.push({ label: "Run again", primary: true, fn: () => run() });
  }
  return errorState({
    title: err.kind === "f8_domain" ? "Every envelope sample left the F8 loading domain"
      : err.kind === "envelope_all_dropped" ? "Every envelope sample failed to evaluate"
      : err.kind === "unsupported" ? "This run is not available for this card" : "The run failed",
    message: [err.message, err.suggestion].filter(Boolean).join(" "),
    actions,
  });
}

function nanReason(s, key, fallback) {
  return s?.[`${key}__nan_reason`] || (key === "g2_op" && s?.g2_op_invalid_reason) || fallback;
}

function renderReadouts(res, { pending = false, partial = false, markChange = false, envelopeRun = null, anim = null } = {}) {
  if (!S?.readouts?.g2_op) return;
  const R = S.readouts;
  const keys = Object.keys(R);
  if (anim?.gap || anim?.runaway) {
    // the playhead sits on a frame the pool has not returned (no number, no interpolation across it),
    // or on/after a thermal runaway (the displayable range has ended: "runaway", never a blend)
    const why = anim.gap || anim.runaway, q = anim.gap ? "not computed" : "runaway";
    for (const k of keys) { R[k].el.classList.remove("is-pending"); R[k].update({ value: null, reason: why, qualifier: q, markChange: false }); }
    return;
  }
  if (pending && !res) {
    if (!S.run) for (const k of keys) R[k].update({ pending: true });
    for (const k of keys) R[k].el.classList.add("is-pending");
    return;
  }
  for (const k of keys) R[k].el.classList.remove("is-pending");
  if (!res) { for (const k of keys) R[k].update(null); R.flux_delivered.el.hidden = true; R.g2_head.el.hidden = true; syncOdd(); return; }
  const s = res.scalars || {};
  const sb = res.scalar_bands || {};
  const tag = normTag(res.tag_chain);
  const labels = res.labels || envelopeRun?.labels || [];
  const isStatic = labels.includes(STATIC_LABEL);
  const nonHead = labels.includes(NON_HEADLINE_LABEL);
  const qual0 = isStatic ? STATIC_LABEL : nonHead ? NON_HEADLINE_LABEL : res.mode === "headline" ? "headline" : "";
  const floorLabel = labels.find((l) => FLOOR_RE.test(l));
  const floorMarker = floorLabel ? { text: floorLabel.replace(/\s*\[[A-Z/]+\]\s*$/, ""), tag: (/\[([A-Z]+)\]\s*$/.exec(floorLabel) || [])[1] || "A", title: "g²(0) here equals the structural floor set by the scalar background b_res: it is not a dynamics result" } : null;
  // envelope runs: scalars without a band are the mid-design point, never a prediction
  const midQ0 = res.mode === "envelope" ? "mid-design" : qual0;
  // animation: values between two computed frames are display interpolation and say so, and are
  // drawn in the dashed [A] chip form whatever the run's tag (rule 2); only a number that was really
  // blended counts (an n/a, a snapped value and T_c are never "interpolated").
  const iQ = (v, base, k) => (anim?.interpKeys?.has(k) && isNum(v) ? (base ? `${base} · interp.` : "interp.") : base);
  const tg = (v, k, t = tag) => (anim?.interpKeys?.has(k) && isNum(v) ? "A" : t);
  const prov = res.provenance || {};
  const m = markChange !== false;
  const opNote = anim ? anim.note : S.opOverride ? `at T_hs = ${fmt(S.opOverride.T)} K` : partial ? "single-T partial result" : "";

  R.g2_op.update({ value: s.g2_op, lo: sb.g2_op?.[0], hi: sb.g2_op?.[1], tag: tg(s.g2_op, "g2_op"), tagNote: prov.drive?.note, qualifier: iQ(s.g2_op, qual0, "g2_op"), caption: opNote || (sb.g2_op ? "interval over the envelope" : ""), reason: nanReason(s, "g2_op", "g²(0) not computed"), marker: floorMarker, markChange: m });
  renderHeadReadout(isStatic && !anim);

  const hasDelivered = "collected_flux_delivered_s" in s;
  R.flux_delivered.el.hidden = !hasDelivered;
  if (hasDelivered) {
    R.flux_delivered.update({ value: s.collected_flux_delivered_s, unit: "s⁻¹", tag: tg(s.collected_flux_delivered_s, "collected_flux_delivered_s"), qualifier: iQ(s.collected_flux_delivered_s, midQ0, "collected_flux_delivered_s"), caption: "delivered photons after the hardware screens", reason: nanReason(s, "collected_flux_delivered_s", "delivered flux not computed"), markChange: m });
  }
  R.flux.setLabel(hasDelivered ? "Collected flux, commanded" : "Collected flux");
  R.flux.update({ value: s.collected_flux_pulsed_s, unit: "s⁻¹", tag: tg(s.collected_flux_pulsed_s, "collected_flux_pulsed_s"), qualifier: iQ(s.collected_flux_pulsed_s, midQ0, "collected_flux_pulsed_s"), caption: hasDelivered ? "idealized: commanded drive, before hardware screens" : (isNum(s.rep_rate_hz) ? `at ${fmt(s.rep_rate_hz)} Hz repetition` : ""), reason: nanReason(s, "collected_flux_pulsed_s", s.flux_measurable === false ? "no pulsed repetition rate on this path, so no collected flux" : "flux not computed"), markChange: m });

  R.brightness.update({ value: s.brightness_per_pulse, tag: tg(s.brightness_per_pulse, "brightness_per_pulse", normTag(prov.brightness?.tag) ? widest([tag, prov.brightness.tag]) : tag), tagNote: prov.brightness?.note, qualifier: iQ(s.brightness_per_pulse, midQ0, "brightness_per_pulse"), caption: s.brightness_convention || "", reason: nanReason(s, "brightness_per_pulse", "brightness not computed"), markChange: m });

  const tcB = sb.T_c;
  const gv = S.ctx.api.gatesNow()?.g2_ceiling?.value;
  const ceilTxt = isNum(gv) ? fmt(gv) : "the g² ceiling";
  R.T_c.update(tcB
    ? { lo: tcB[0], hi: tcB[1], unit: "K", tag, caption: !isNum(tcB[0]) || !isNum(tcB[1]) ? `one edge does not cross ${ceilTxt} in range` : "interval over the envelope", markChange: m }
    : { value: s.T_c, unit: "K", tag, caption: anim ? anim.tcNote : "", reason: nanReason(s, "T_c", (res.curves?.T_hs?.length || 0) <= 1 ? "needs the full T_hs grid (single-T run)" : `g²(0) never crosses ${ceilTxt} in the T_hs range`), markChange: m });

  R.T_j_op.update({ value: s.T_j_op, unit: "K", tag: tg(s.T_j_op, "T_j_op"), qualifier: iQ(s.T_j_op, res.mode === "envelope" ? "mid-design" : "", "T_j_op"), caption: anim ? anim.note : isNum(s.dT_J) ? `self-heating ΔT_J ${fmt(s.dT_J)} K` : "", reason: nanReason(s, "T_j_op", "junction temperature not computed"), markChange: m });
  syncOdd();
}

/** Last readout spans both columns when an odd number is visible. */
function syncOdd() {
  const vis = [...S.dom.readouts.children].filter((e) => !e.hidden);
  for (const e of S.dom.readouts.children) e.classList.remove("is-span");
  if (vis.length % 2) vis.at(-1).classList.add("is-span");
}

// ---------------------------------------------------------------- seed from Results (rerunSeed)
// "Open card in Designer" from a Results row: the CSV values beside the live values of the seeded
// design (rt_edge rows map one-to-one onto card fields; other campaigns' rows are shown, not mapped).
function renderSeed() {
  const host = S?.dom?.seedHost;
  if (!host) return;
  clear(host);
  const sd = S.seed;
  if (!sd) return;
  const row = sd.row || {};
  const id = row.config_id ?? row.row_id ?? "";
  const close = h("button.icon-btn", { type: "button", "aria-label": "Dismiss the seeded row", on: { click: () => { S.seed = null; renderSeed(); } } }, icon("close", { size: 16 }));
  const head = h("header.seed-head", h("div", h("h2.engrave", "Seeded from a sweep row"), h("p.seed-sub", `${sd.campaign} · row ${id}`)), close);
  if (!sd.design) {
    host.append(h("section.seed", { "aria-label": "Seeded sweep row" }, head,
      h("p.seed-note", icon("info", { size: 14 }), `${sd.campaign} rows are not mapped onto card fields (its runner applies geometry presets the Studio does not replicate), so the card runs as written and no live value is set beside the CSV.`)));
    return;
  }
  const live = S.run;
  const ok = live && S.runT === row.T_hs_K;
  const pairs = rtEdgePairs(row, ok ? live.scalars : {});
  const tagLive = ok ? normTag(live.tag_chain) : null;
  const labels = ok ? (live.labels || []) : [];
  host.append(h("section.seed", { "aria-label": "Seeded sweep row", dataset: { seeded: "1" } }, head,
    h("p.seed-note", `The card carries the row's grid values and model switches (${rtEdgeSeedFields(row).map(([k]) => k.split(".").pop()).join(", ")}); evaluated as a point at T_hs = ${fmt(row.T_hs_K)} K.`),
    labels.length ? h("p.seed-labels", labels.map((l) => h("span.label-chip", l))) : null,
    h("table.data-table.cmp-table.seed-table",
      h("thead", h("tr", h("th", "quantity"), h("th.num", "CSV (committed)"), h("th.num", "live now"), h("th.num", "match"))),
      h("tbody", pairs.map(([l, a, b, u, key]) => h("tr",
        h("th", { scope: "row" }, l),
        h("td.num", sro({ label: "csv", value: a, unit: u, tag: "A", note: `out/${sd.campaign} row ${id}`, size: "sm" })),
        h("td.num", ok ? sro({ label: "live", value: b, unit: u, tag: tagLive, size: "sm", caption: b == null ? liveNanReason(live.scalars, key) : "" }) : h("span.nan", S.runState === "running" ? "running" : "run the seeded design")),
        h("td.num.match", ok && a === b ? h("span.m-ok", icon("check", { size: 14 }), "bit-identical") : ok && isNum(a) && isNum(b) ? h("span.m-diff", `Δ ${fmt(b - a)}`) : "n/a")))))));
}

// ---------------------------------------------------------------- headline beside static (H21)
// A point run on a static RT-edge card shows static numbers; the headline-model g2 at the same
// T_hs is run as a second job (mode headline, T_grid [T_hs]) and shown beside them.
function renderHeadReadout(isStatic) {
  const R = S.readouts.g2_head;
  const show = isStatic && S.run?.mode === "point";
  R.el.hidden = !show;
  if (!show) return;
  const hd = S.head;
  if (!hd) { R.update({ pending: true }); return; }
  if (hd.error) { R.update({ value: null, tag: "A", reason: `headline run failed: ${hd.error}` }); return; }
  const r = hd.result;
  R.update({ value: r.scalars?.g2_op, tag: normTag(r.tag_chain), qualifier: "headline", caption: `headline model at T_hs = ${fmt(hd.T)} K · ${r.run_id}`, reason: nanReason(r.scalars, "g2_op", "g²(0) not computed"), markChange: false });
}

function runHeadlineBeside(T) {
  S.headJob?.cancel?.();
  S.head = null;
  const job = S.ctx.api.runJob({ card: S.card, design: S.design, mode: "headline", T_grid: [T] });
  S.headJob = job;
  job.promise.then((r) => {
    if (!S || S.headJob !== job) return;
    S.headJob = null;
    S.head = { T, result: r.result };
    renderHeadReadout(true);
    syncOdd();
  }).catch((e) => {
    if (!S || S.headJob !== job) return;
    S.headJob = null;
    if (e.body?.error?.kind === "cancelled") return;
    S.head = { T, error: e.message };
    renderHeadReadout(true);
  });
}

function inspectKpi(key, spec) {
  const r = S.opOverride?.result || S.run;
  if (!r) return;
  const prov = r.provenance || {};
  const tagT = normTag(r.tag_chain);
  const notes = Object.entries(prov).filter(([, v]) => v && v.tag).map(([k, v]) => ({ key: k, tag: v.tag, note: v.note }))
    .sort((a, b) => (normTag(b.tag) === tagT) - (normTag(a.tag) === tagT));
  const labelMap = { g2_op: "g²(0) at the operating point", flux: "Collected flux", flux_delivered: "Collected flux, delivered", brightness: "Brightness per pulse", T_c: "Ceiling temperature T_c", T_j_op: "Junction temperature T_j" };
  const valueText = spec ? (isNum(spec.lo) || isNum(spec.hi) ? `${fmtInterval(spec.lo, spec.hi)} ${spec.unit || ""}` : isNum(spec.value) ? `${fmt(spec.value)} ${spec.unit || ""}` : `n/a: ${spec.reason || ""}`) : "n/a";
  S.ctx.inspector.number({ label: labelMap[key] || key, key: `${r.run_id} / ${key}`, valueText, tag: r.tag_chain, notes, source: r.mode === "envelope" ? r.scalars_source || "evaluate_envelope" : "fsim_core.device.evaluate scalars" });
}

// ---------------------------------------------------------------- details drawer
function groupOf(k) {
  if (k.includes(".")) return k.split(".")[0];
  for (const p of ["edge_", "finite_pulse_", "cw_", "aperture_", "photon_budget", "track_", "injection", "background", "retention", "nitride_", "nanowire_", "wire_", "set_"]) if (k.startsWith(p)) return p.replace(/_$/, "");
  return "core";
}

function renderDetails() {
  const r = S.opOverride?.result || S.run;
  const body = S.dom.details.querySelector(".drawer-body");
  const count = S.dom.details.querySelector(".drawer-count");
  clear(body);
  if (!r) { count.textContent = ""; body.append(h("p.chart-sub.pad", "Run the design to list every scalar.")); return; }
  const s = r.scalars || {};
  const skip = (k) => k === "tag_chain" || k === "provenance" || k.endsWith("__nan_reason");
  const keys = Object.keys(s).filter((k) => !skip(k));
  count.textContent = `${keys.length} values, all within tag chain`;
  const groups = {};
  for (const k of keys) (groups[groupOf(k)] ||= []).push(k);
  const tag = normTag(r.tag_chain);
  const order = Object.keys(groups).sort((a, b) => (a === "core" ? -1 : b === "core" ? 1 : a.localeCompare(b)));
  for (const g of order) {
    body.append(h("section.dgroup",
      h("h3.dgroup-h", h("span.engrave", g), tagChip(tag, { size: "sm", title: "every value in this run carries the run's tag chain" })),
      h("dl.dgroup-list", groups[g].map((k) => {
        const v = s[k];
        let text;
        if (v === null) text = h("span.nan", s[`${k}__nan_reason`] || "not computed on this path");
        else if (Array.isArray(v)) text = v.length ? v.join("; ") : h("span.nan", "none");
        else if (typeof v === "object") text = Object.entries(v).map(([a, b]) => `${a} ${typeof b === "number" ? fmt(b) : b}`).join("; ");
        else text = typeof v === "number" ? `${fmt(v)}${SCALAR_UNITS[k] ? ` ${SCALAR_UNITS[k]}` : ""}` : String(v);
        const band = r.scalar_bands?.[k];
        return h("div.drow", h("dt", k), h("dd", band ? `${fmtInterval(band[0], band[1])}${SCALAR_UNITS[k] ? ` ${SCALAR_UNITS[k]}` : ""}` : text));
      }))));
  }
}

// ---------------------------------------------------------------- 3D viewport
async function loadScene(T) {
  if (!S) return;
  const seq = ++S.sceneSeq;
  const slot = S.dom.viewport;
  const theme = S.ctx.state.get("theme");
  if (S.draft) { placeholder2D(null, "The 3D model is built from a card: save this design as a card to view it"); return; }
  if (!S.scene) {
    clear(slot).append(h("div.vp-placeholder", h("div.vp-plate", h("span.engrave", "3D model loading"), h("span.vp-sub", S.card))));
  }
  let spec = null;
  try {
    const params = { card: S.card };
    if (isNum(T)) params.T = T;
    spec = await S.ctx.api.sceneSpec("device", params);
  } catch (e) {
    if (!S || seq !== S.sceneSeq) return;
    if (!S.scene) placeholder2D(null, e.message);
    return;
  }
  if (!S || seq !== S.sceneSeq) return;
  if (S.scene) { try { S.scene.update(spec); } catch (e) { console.warn("scene update failed", e); } return; }
  let mod = null;
  try { mod = await import("../viz3d/index.js"); } catch { mod = null; }
  if (!S || seq !== S.sceneSeq) return;
  if (!mod?.mountScene) { placeholder2D(spec); return; }
  try {
    clear(slot);
    S.scene = mod.mountScene(slot, "device", spec, { theme });
  } catch (e) {
    console.warn("3D mount failed", e);
    placeholder2D(spec, "3D view unavailable on this machine");
  }
}

function placeholder2D(spec, why) {
  const slot = S.dom.viewport;
  const layers = (spec?.layers || []).filter((l) => isNum(l.t_nm) && l.t_nm > 0).slice(0, 10);
  const total = layers.reduce((a, l) => a + Math.log10(1 + l.t_nm), 0) || 1;
  clear(slot).append(h("div.vp-placeholder",
    h("div.vp-plate",
      h("span.engrave", "3D model loading"),
      h("span.vp-sub", why || spec?.title || S.card),
      layers.length ? h("div.vp-stack", { role: "img", "aria-label": "Layer stack, thickness on a log scale" },
        layers.map((l) => h("div.vp-layer", { style: { flexGrow: String(Math.log10(1 + l.t_nm) / total) } },
          h("span.vp-ln", l.name), h("span.vp-lt", `${fmt(l.t_nm)} nm`), tagChip(l.tag, { size: "sm" })))) : null,
      layers.length ? h("span.vp-sub", "Section on a log thickness scale (not true scale)") : null)));
}

// ================================================================ PRESET BUILDER, SAVE AS CARD, UNSAVED (studio-p2a)
// "New from presets" composes a design from the five fsim_core.presets axes through
// POST /api/presets/apply (preset_device + apply_injection_preset, the DPG designer's order) and
// opens it as an unsaved design (#/design/unsaved). "Save as card" PUTs it to cards/studio/ via
// DeviceDesign.save. Picking a preset never upgrades a field's provenance: the summary shows the
// META tag of every field, and an option's chip is the tag printed in its presets.py note or,
// without one, the widest META tag of the fields it sets.
const DRAFT_CARD = "unsaved";
const PRESET_AXES = [
  { id: "dot", label: "Dot" },
  { id: "template", label: "Thermal template" },
  { id: "cavity", label: "Cavity" },
  { id: "drive", label: "Drive" },
  { id: "injection", label: "Injection" },
];
// fsim_gui.designer's combo defaults
const PRESET_DEFAULT = { dot: "chatzarakis-class", template: "GaAs", cavity: "none", drive: "cw-electrical", injection: "standard" };
const SUMMARY_PATHS = ["dot.delta_xx", "dot.gamma_scale", "dot.r_xx", "cavity.F_P", "cavity.kappa", "cavity.G", "drive.V", "drive.I_uA", "drive.duty", "drive.mu", "drive.b_e", "drive.F_p", "drive.eta_capture", "thermal.T_hs"];
const SLUG_RE = /^[a-z0-9]+(?:-[a-z0-9]+)*$/;
const NOTE_TAG = /\[(V|DR|E|A)\b[^\]]*\]/;

function comboText(sel) { return PRESET_AXES.map((a) => sel[a.id]).join(" / "); }

function draftDesc() {
  const p = S.draft?.presets || {};
  return `Built from presets: ${PRESET_AXES.map((a) => `${a.label.toLowerCase()} ${p[a.id]}`).join(", ")}. Every preset is a point choice inside its META band; each field keeps its META tag. Not saved as a card yet.`;
}

function presetFields(axis, p) {
  const skip = new Set(["note", "source"]);
  if (axis === "template") return [["thermal.layers", (p.layers || []).map((L) => L.name).join(" + ")], ["thermal.substrate", p.substrate?.name]];
  if (axis === "injection") {
    const out = Object.entries(p.drive || {}).map(([k, v]) => [`drive.${k}`, v]);
    if (p.halve_b_e) out.push(["drive.b_e", "× 0.5"]);
    if (p.layer) out.push(["thermal.layers", `+ ${p.layer.name}`]);
    return out;
  }
  const block = axis === "drive" ? "drive" : axis === "cavity" ? "cavity" : "dot";
  return Object.entries(p).filter(([k]) => !skip.has(k)).map(([k, v]) => [`${block}.${k}`, v]);
}

function presetInfo(axis, key, p) {
  const META = S.metaAll || {};
  const fields = presetFields(axis, p);
  const noteTag = typeof p.note === "string" ? (NOTE_TAG.exec(p.note) || [])[1] : null;
  const metaTag = widest(fields.map(([path]) => META[path]?.tag).filter(Boolean));
  const desc = typeof p.note === "string"
    ? p.note.replace(/\s*\[[^\]]*\]\s*/g, " ").trim() + (p.source ? ` (${p.source})` : "")
    : fields.length
      ? fields.map(([path, v]) => `${path.split(".")[1]} ${typeof v === "number" ? fmt(v) : v === true ? "on" : v === false ? "off" : v}${typeof v === "number" ? (unitText(META[path]?.unit) ? ` ${unitText(META[path].unit)}` : "") : ""}`).join(" · ")
      : "no change";
  // nothing recorded (e.g. the thermal stacks: no META entry, no note): the widest tag, [A], is
  // assumed, the same rule as a run whose tag_chain is missing (fsim_studio/tags.py)
  const assumed = !noteTag && !metaTag && fields.length > 0;
  return {
    desc, tag: noteTag || metaTag || (assumed ? "A" : null),
    tagTitle: noteTag ? "tag printed in this preset's note (fsim_core/presets.py); the design fields keep their META tag"
      : metaTag ? "widest META tag of the fields this preset sets"
      : assumed ? "no tag recorded for these values in presets.py or META; assumed widest" : null,
  };
}

function togglePresetSheet() {
  const host = S.dom.presetHost;
  const open = !host.classList.contains("is-open");
  S.dom.presetBtn?.setAttribute("aria-expanded", String(open));
  if (!open) { host.classList.remove("is-open"); clear(host); return; }
  S.presetSel = S.presetSel || { ...PRESET_DEFAULT, ...(S.draft?.presets || {}) };
  host.classList.add("is-open");
  const axes = PRESET_AXES.map((ax) => {
    const opts = S.presets?.[ax.id] || {};
    return h("fieldset.preset-axis", { dataset: { axis: ax.id } },
      h("legend.engrave", ax.label),
      h("div.po-list", { role: "radiogroup", "aria-label": `${ax.label} preset` },
        Object.entries(opts).map(([key, p]) => {
          const info = presetInfo(ax.id, key, p);
          return h("button.po", {
            type: "button", role: "radio", dataset: { axis: ax.id, key }, "aria-checked": String(S.presetSel[ax.id] === key),
            on: { click: () => { S.presetSel[ax.id] = key; syncPresetRadios(); presetApply(); } },
          },
          h("span.po-head", h("span.po-name", key), tagChip(info.tag, { size: "sm", title: info.tagTitle })),
          h("span.po-desc", info.desc));
        })));
  });
  S.dom.presetSummary = h("div.preset-summary", { "aria-live": "polite" });
  S.dom.presetOpen = h("button.btn.btn-outline.preset-open", { type: "button", disabled: true, on: { click: openDraft } }, icon("chevronRight", { size: 16 }), "Open in Designer");
  const sheet = h("section.sheet.preset-sheet", { "aria-label": "New device from presets" },
    h("header.sheet-head",
      h("div", h("h2.engrave.sheet-title", "New device from presets"),
        h("p.sheet-sub", "Five preset axes from fsim_core/presets.py, composed by preset_device + apply_injection_preset (the desktop designer's order). Opens as an unsaved design.")),
      h("button.icon-btn", { type: "button", "aria-label": "Close the preset builder", on: { click: () => togglePresetSheet() } }, icon("close", { size: 18 }))),
    h("div.preset-body", h("div.preset-axes", axes),
      h("aside.preset-side", h("h3.engrave", "Resulting design"), S.dom.presetSummary,
        h("div.preset-actions", S.dom.presetOpen))));
  clear(host).append(sheet);
  if (!reducedMotion()) sheet.classList.add("mount-in");
  presetApply();
  requestAnimationFrame(() => host.scrollIntoView({ block: "nearest", behavior: reducedMotion() ? "auto" : "smooth" }));
}

function syncPresetRadios() {
  for (const b of S.dom.presetHost.querySelectorAll(".po")) b.setAttribute("aria-checked", String(S.presetSel[b.dataset.axis] === b.dataset.key));
}

let presetSeq = 0;
async function presetApply() {
  const seq = ++presetSeq;
  const sel = { ...S.presetSel };
  S.dom.presetOpen.disabled = true;
  S.dom.presetSummary.classList.add("is-pending");
  let res;
  try {
    res = await S.ctx.api.postJSON("/api/presets/apply", { ...sel, name: "preset-device" });
  } catch (e) {
    if (!S || seq !== presetSeq) return;
    clear(S.dom.presetSummary).append(h("p.fr-warn", `presets/apply failed: ${e.message}`));
    return;
  }
  if (!S || seq !== presetSeq || !S.dom.presetSummary) return;
  S.presetBuilt = { sel, design: res.design, meta: res.meta || {} };
  const d = res.design, meta = res.meta || {};
  const rows = SUMMARY_PATHS.filter((p) => meta[p] && getPath(d, p) !== undefined && getPath(d, p) !== null).map((p) => {
    const v = getPath(d, p);
    return sro({ label: meta[p].label || p, value: typeof v === "boolean" ? (v ? "on" : "off") : v, unit: unitText(meta[p].unit), tag: meta[p].tag, note: meta[p].source, size: "sm", key: p });
  });
  const stack = (d.thermal?.layers || []).map((L) => L.name).join(" + ");
  clear(S.dom.presetSummary).append(
    h("p.preset-combo", h("span.chip.chip-mono", comboText(sel))),
    h("p.preset-line", `${d.cavity?.enabled ? `cavity on (${d.cavity.type})` : "no cavity"} · stack ${stack || "none"} on ${d.thermal?.substrate?.name || "?"} · platform ${d.platform}`),
    h("div.preset-grid", rows));
  S.dom.presetSummary.classList.remove("is-pending");
  S.dom.presetOpen.disabled = false;
}

function openDraft() {
  const b = S.presetBuilt;
  if (!b) return;
  const draft = { design: b.design, meta: b.meta, presets: { ...b.sel }, title: `unsaved: ${comboText(b.sel)}`, stamp: Date.now() };
  S.ctx.state.set({ draft });
  S.ctx.navigate(`#/design/${DRAFT_CARD}`);
}

function noteEdit() {
  S.edits++;
  syncUnsaved();
}

function syncUnsaved() {
  if (!S) return;
  const u = S.draft ? { title: S.draft.title, short: "unsaved design", kind: "draft" }
    : S.edits ? { title: `${(S.card || "").replace(/-design$/, "")}: edited, not saved`, short: "edited, not saved", kind: "edits" } : null;
  S.ctx.state.set({ unsaved: u });
  if (S.dom?.saveBtn) S.dom.saveBtn.hidden = !u;
}

function suggestName() {
  if (S.draft) {
    const p = S.draft.presets || {};
    return `preset-${p.dot}-${p.cavity}-${p.drive}`.toLowerCase().replace(/[^a-z0-9]+/g, "-").replace(/^-|-$/g, "").slice(0, 70);
  }
  const c = S.ctx.cards().find((x) => x.name === S.card);
  return c?.saved ? S.card : `${S.card.replace(/-design$/, "")}-edited`;
}

function saveAsPopover(anchor) {
  const input = h("input.fr-num.save-name", { type: "text", value: suggestName(), spellcheck: false, "aria-label": "Card name (lower-case slug)", "aria-describedby": "save-msg" });
  const msg = h("p.save-msg#save-msg");
  const go = h("button.btn.btn-outline.save-go", { type: "submit" }, icon("download", { size: 16 }), "Save");
  const check = () => {
    const v = input.value.trim();
    const shipped = S.ctx.cards().some((c) => c.name === v && !c.saved);
    const bad = !SLUG_RE.test(v) || v.length > 80 || v === DRAFT_CARD;
    msg.textContent = bad ? "Use a lower-case slug: letters, digits and single hyphens." : shipped ? "That is a shipped card; it is never overwritten. Pick a new name." : `Saves cards/studio/${v}.yaml through DeviceDesign.save.`;
    msg.dataset.state = bad || shipped ? "bad" : "ok";
    go.disabled = bad || shipped;
    return !(bad || shipped);
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
          const r = await S.ctx.api.putJSON(`/api/cards/${encodeURIComponent(name)}`, { design: S.design, presets: S.draft?.presets || null });
          S.leaveOK = true;
          S.edits = 0;
          S.ctx.state.set({ draft: null, unsaved: null });
          await S.ctx.refreshCards?.();
          closePopover();
          toast(`Saved ${r.path}`, { kind: "ok" });
          if (S.card === name) { S.leaveOK = false; syncUnsaved(); } else S.ctx.navigate(`#/design/${name}`);
        } catch (err) {
          msg.textContent = err.status === 409 ? `Refused (409): ${err.message}` : err.message;
          msg.dataset.state = "bad";
          go.disabled = false;
        }
      },
    },
  }, h("h2.why-title", "Save as card"),
  h("label.save-label", { for: "save-name" }, "Card name"), input, msg, h("div.save-actions", go));
  input.id = "save-name";
  popover(anchor, form, { label: "Save as card" });
  check();
  input.focus();
  input.select();
}

/** app.js asks before leaving the route; an unsaved draft or sheet edits hold it once. */
export function beforeLeave() {
  if (!S || S.leaveOK || !(S.draft || S.edits)) return true;
  if (S.leaveAsk) return S.leaveAsk;
  S.leaveAsk = new Promise((resolve) => openLeaveDialog(resolve));
  return S.leaveAsk;
}

let leaveDlg = null;
function closeLeaveDialog() { leaveDlg?.remove(); leaveDlg = null; }

function openLeaveDialog(resolve) {
  closeLeaveDialog();
  const what = S.draft ? S.draft.title : `${S.card.replace(/-design$/, "")} with ${S.edits} edit${S.edits === 1 ? "" : "s"}`;
  const done = (ok) => {
    closeLeaveDialog();
    if (!S) { resolve(ok); return; }
    S.leaveAsk = null;
    if (ok) { S.leaveOK = true; S.ctx.state.set({ unsaved: null, draft: null }); }
    resolve(ok);
  };
  const stay = h("button.btn.btn-outline.leave-stay", { type: "button", on: { click: () => done(false) } }, "Stay");
  const save = h("button.btn.btn-quiet.leave-save", { type: "button", on: { click: () => { done(false); if (S?.dom?.saveBtn) saveAsPopover(S.dom.saveBtn); } } }, "Save as card first");
  const leave = h("button.btn.leave-go", { type: "button", on: { click: () => done(true) } }, "Leave without saving");
  leaveDlg = h("div.leave-scrim", { on: { keydown: (e) => { if (e.key === "Escape") { e.stopPropagation(); done(false); } } } },
    h("section.leave-dialog.sheet", { role: "alertdialog", "aria-modal": "true", "aria-labelledby": "leave-t", "aria-describedby": "leave-b" },
      h("h2.engrave#leave-t", "Leave without saving?"),
      h("p#leave-b", `${what} is not saved as a card. Leaving discards it; Save as card keeps it in cards/studio/.`),
      h("div.leave-actions", save, stay, leave)));
  document.body.append(leaveDlg);
  stay.focus();
}

// ================================================================ RESPONSE ANIMATION (studio-p2a)
// POST /api/animate computes one evaluate() per parameter value on the CPU pool (cached per frame
// under the /api/run key) and streams `frame` events. Playback happens here, in the browser:
// the shelf readouts, the hero chart's playhead marker and the beam follow the playhead; values
// between two COMPUTED neighbouring frames are linear display interpolation and carry "interp.";
// Honesty rules (physics brief section 7): a blend needs two ADJACENT computed frames that are both
// g2_op_valid, neither runaway and carry identical label sets (interpRefusal); otherwise the
// playhead snaps to the nearer computed frame and says so. T_c, booleans, labels, __nan_reason and
// verdict words are never blended: they come from the NEAREST computed frame (rule 3), T_c from
// that frame's own evaluate() (null, with its reason, on a grid "op" frame). Blended numbers are
// marked "interp." and drawn in the dashed [A] chip form; the chip returns to the run's tag on a
// computed frame. A runaway frame ends the displayable range. The 3D scene follows the playhead
// through viz3d's cheap setPlayhead glow hook and is rebuilt only when the playhead rests.
const ANIM_MAX = 120;
const ANIM_INTERP = ["g2_op", "collected_flux_pulsed_s", "collected_flux_delivered_s", "brightness_per_pulse", "T_j_op", "eps_op", "rho_op"];
const ANIM_REGIME = ["flat_band", "depletion_regime"]; // a regime change blocks a blend (nitride Stark, rule 6)
const A_SLICE_NOTE = "one range at a time, not a prediction"; // rule 7: an animated [A] input

const labelKey = (f) => [...(f.labels || [])].sort().join("\u0001");
const frameOK = (f) => !!f && !f.error;
const frameRunaway = (A, i) => A.runawayIdx != null && i >= A.runawayIdx;

/** Why two adjacent computed frames may NOT be blended (null: they may). Rules 1 and 5. */
function interpRefusal(A, i, j) {
  const a = A.frames[i], b = A.frames[j];
  if (!frameOK(a) || !frameOK(b)) return "a neighbouring frame is missing";
  if (j !== i + 1) return "the frames are not adjacent";
  if (a.scalars?.runaway === true || b.scalars?.runaway === true || frameRunaway(A, j)) return "thermal runaway ends the displayable range";
  if (a.scalars?.g2_op_valid !== true || b.scalars?.g2_op_valid !== true) return "g²(0) is not valid at one of the two frames";
  if (labelKey(a) !== labelKey(b)) return "the two frames carry different labels";
  for (const k of ANIM_REGIME) if ((k in (a.scalars || {}) || k in (b.scalars || {})) && a.scalars?.[k] !== b.scalars?.[k]) return "regime change (flat band / depletion regime) between the frames";
  return null;
}

/** Recompute the runaway cut (the first runaway frame ends the range) and tell the timeline. */
function syncRunaway(A, serverAt = null) {
  let r = isNum(serverAt) ? serverAt : null;
  A.frames.forEach((f, i) => { if (frameOK(f) && f.scalars?.runaway === true && (r == null || i < r)) r = i; });
  A.runawayIdx = r;
  A.tl?.setLimit(r);
}

function animPaths() {
  return Object.keys(S.meta).filter((p) => {
    const m = S.meta[p];
    const v = getPath(S.design, p);
    return Array.isArray(m.band) && isNum(m.band[0]) && isNum(m.band[1]) && m.band[1] > m.band[0]
      && !m.choices && m.unit !== "enum" && m.unit !== "bool" && typeof v !== "boolean" && (isNum(v) || v === null);
  });
}

function addAnimGlyph(row, p) {
  if (!animPaths().includes(p)) return;
  const tools = row.querySelector(".fr-tools");
  tools?.append(h("button.icon-btn.fr-anim", {
    type: "button", title: `Animate the response over ${S.meta[p].label || p}`, "aria-label": `Animate the response over ${p}`,
    on: { click: () => openAnim(p) },
  }, icon("run", { size: 14 })));
}

function paramLabel(A) { return A.meta?.label || A.param.split(".").pop(); }
function paramText(A, v) { const u = unitText(A.meta?.unit); return `${fmt(v)}${u ? ` ${u}` : ""}`; }

function openAnim(param) {
  closeAnim({ quiet: true });
  const m = S.meta[param] || {};
  const mode = S.mode === "headline" ? "headline" : "point";
  const cls = platformClass(S.card, S.design, mode);
  S.anim = {
    param, meta: m, lo: m.band[0], hi: m.band[1], n: cls === "nanowire" ? 16 : 48, mode, cls,
    grid: ["edge", "nanowire"].includes(cls) ? "op" : "full",
    phase: "setup", frames: [], values: [], tl: null, dom: {}, heldShown: undefined,
  };
  S.dom.hsAnim?.setAttribute("aria-expanded", String(param === "thermal.T_hs"));
  renderAnimSheet();
  requestAnimationFrame(() => S?.dom?.animHost?.scrollIntoView({ block: "nearest", behavior: reducedMotion() ? "auto" : "smooth" }));
}

function animEta(A) {
  const nT = A.grid === "op" ? 1 : (NITRIDE_GRID_N[S.design?.platform] || 120);
  const per = evalSeconds(A.cls, nT);
  return { per, total: (per * A.n) / POOL };
}

function renderAnimSheet() {
  const A = S.anim;
  const host = S.dom.animHost;
  const paths = animPaths();
  const unit = unitText(A.meta.unit);
  const sel = h("select.fr-select.anim-param", { "aria-label": "Parameter to animate", on: { change: (e) => openAnim(e.target.value) } },
    paths.map((p) => h("option", { value: p, selected: p === A.param }, p)));
  const num = (cls, v, label, step = "any") => h(`input.fr-num.${cls}`, { type: "number", step, value: String(v), "aria-label": label });
  A.dom.lo = num("anim-lo", A.lo, "From");
  A.dom.hi = num("anim-hi", A.hi, "To");
  A.dom.n = num("anim-n", A.n, "Frames", "1");
  A.dom.n.min = "2"; A.dom.n.max = String(ANIM_MAX);
  // review finding 1: From/To live inside the META validity band (the backend refuses anything outside it)
  A.dom.lo.min = A.dom.hi.min = String(A.meta.band[0]);
  A.dom.lo.max = A.dom.hi.max = String(A.meta.band[1]);
  A.dom.eta = h("p.anim-eta");
  const syncEta = () => {
    A.n = Math.round(Number(A.dom.n.value));
    const { per, total } = animEta(A);
    A.dom.eta.textContent = `${Number.isFinite(A.n) ? A.n : "?"} frames × ${fmtSeconds(per)} per evaluate() on ${POOL} workers: about ${fmtSeconds(total)}; frames already in the cache are free. ${A.grid === "op" ? "Each frame is the operating point only, so T_c is not available." : "Each frame runs the full T_hs grid."}`;
  };
  A.dom.n.addEventListener("input", syncEta);
  syncEta();
  A.dom.go = h("button.btn.btn-outline.anim-go", { type: "button", on: { click: startAnim } }, icon("run", { size: 16 }), "Compute frames");
  A.dom.cancel = h("button.btn.btn-quiet.anim-cancel", { type: "button", hidden: true, on: { click: cancelAnimJob } }, icon("stop", { size: 16 }), "Cancel");
  A.dom.status = h("p.anim-status", { "aria-live": "polite" });
  A.dom.player = h("div.anim-player", { hidden: true });
  A.dom.labels = h("div.anim-labels", { "aria-live": "polite" });
  A.dom.sheet = h("section.sheet.anim-sheet", { "aria-label": "Response animation", dataset: { phase: "setup", param: A.param } },
    h("header.sheet-head",
      h("div", h("h2.engrave.sheet-title", "Animate the response"),
        h("p.sheet-sub", "Frames are computed once on the CPU pool (one evaluate() each), then played back here. Between computed frames the shelf interpolates for display and says interp.")),
      h("button.icon-btn", { type: "button", "aria-label": "Close the animation", on: { click: () => closeAnim() } }, icon("close", { size: 18 }))),
    h("div.anim-setup",
      h("label.anim-field", h("span.engrave", "Parameter"), sel, tagChip(A.meta.tag, { size: "sm", title: A.meta.source })),
      h("label.anim-field", h("span.engrave", "From"), A.dom.lo, h("span.fr-unit", unit)),
      h("label.anim-field", h("span.engrave", "To"), A.dom.hi, h("span.fr-unit", unit)),
      h("label.anim-field", h("span.engrave", "Frames"), A.dom.n),
      A.meta?.tag === "A" ? h("span.label-chip.anim-slice", { title: "an assumed [A] input: each frame is one value of it, run one at a time" }, A_SLICE_NOTE) : null,
      h("span.label-chip.anim-mode", { title: A.mode === "headline" ? "RT-edge headline switches, as the Headline run mode" : "point evaluate(), as the Point run mode" }, `${A.mode} model`),
      h("div.anim-keys", A.dom.cancel, A.dom.go)),
    A.dom.eta,
    A.dom.status,
    A.dom.player);
  A.dom.sheet.__anim = {
    seek: (p) => A.tl?.seek(p), where: () => A.tl?.where(), play: () => A.tl?.play(), pause: () => A.tl?.pause(),
    get pos() { return A.tl?.pos; }, get frames() { return A.frames; }, get values() { return A.values; }, get phase() { return A.phase; },
    states: () => A.tl?.states(),
    // test hooks (not part of the shell contract): the honesty rule and a synthetic-frame patch
    refusal: (i, j) => interpRefusal(A, i, j),
    _patch: (i, scalars, labels, tag) => {
      const f = A.frames[i];
      if (!f) return false;
      f.scalars = { ...f.scalars, ...scalars };
      if (labels) f.labels = labels;
      if (tag) f.tag_chain = tag;
      syncRunaway(A, A.serverRunaway);
      A.patchRev = (A.patchRev || 0) + 1;
      A.heldShown = undefined;
      A.tl.refresh();
      return true;
    },
    get runawayIdx() { return A.runawayIdx; },
  };
  clear(host).append(A.dom.sheet);
  host.classList.add("is-open");
  if (!reducedMotion()) A.dom.sheet.classList.add("mount-in");
}

async function startAnim() {
  const A = S?.anim;
  if (!A) return;
  const lo = Number(A.dom.lo.value), hi = Number(A.dom.hi.value), n = Math.round(Number(A.dom.n.value));
  if (!Number.isFinite(lo) || !Number.isFinite(hi) || lo === hi) { A.dom.status.textContent = "From and To must be numbers that differ."; return; }
  const band = A.meta?.band;
  if (Array.isArray(band) && (Math.min(lo, hi) < band[0] || Math.max(lo, hi) > band[1])) {
    A.dom.status.textContent = `From and To must lie inside the validity band of ${A.param}: ${fmt(band[0])} to ${fmt(band[1])}${unitText(A.meta.unit) ? ` ${unitText(A.meta.unit)}` : ""}.`;
    return;
  }
  if (!(n >= 2 && n <= ANIM_MAX)) { A.dom.status.textContent = `Frames must be between 2 and ${ANIM_MAX}.`; return; }
  stopAnimJob(A);
  A.tl?.dispose();
  Object.assign(A, { lo, hi, n, frames: new Array(n).fill(null), phase: "computing", heldShown: undefined, started: false, jobOver: false, runawayIdx: null, serverRunaway: null, density: null, lastGlow: null, crossShown: "" });
  A.dom.sheet.dataset.phase = "computing";
  A.dom.go.disabled = true;
  A.dom.cancel.hidden = false;
  A.dom.status.textContent = "Submitting the frames";
  let resp;
  try {
    resp = await S.ctx.api.postJSON("/api/animate", { card: S.card, design: S.design, param: A.param, lo, hi, n, mode: A.mode, grid: A.grid });
  } catch (e) {
    if (S?.anim !== A) return;
    A.phase = "setup";
    A.dom.sheet.dataset.phase = "setup";
    A.dom.go.disabled = false;
    A.dom.cancel.hidden = true;
    A.dom.status.textContent = `Animation refused: ${e.message}`;
    return;
  }
  if (S?.anim !== A) return;
  A.values = resp.values;
  A.jobId = resp.job_id;
  A.dom.status.textContent = resp.cached ? "Every frame came from the cache." : `${resp.n_cached || 0} of ${n} frames from the cache; computing the rest (estimated ${fmtSeconds(resp.eta_s)}).`;
  buildPlayer(A);
  if (resp.result) { ingestResult(A, resp.result); finishAnim(A); return; }
  followAnim(A);
}

function followAnim(A) {
  let es;
  try { es = new EventSource(`/api/jobs/${A.jobId}/events`); } catch { pollAnim(A); return; }
  A.es = es;
  const parse = (ev) => { try { return JSON.parse(ev.data); } catch { return {}; } };
  es.addEventListener("frame", (ev) => { if (S?.anim === A) { const d = parse(ev); ingestFrame(A, d.index, d.frame); } });
  es.addEventListener("done", (ev) => { es.close(); if (S?.anim === A) { ingestResult(A, parse(ev).result); finishAnim(A); } });
  es.addEventListener("state", (ev) => { if (parse(ev).state === "cancelled") { es.close(); if (S?.anim === A) finishAnim(A, "cancelled"); } });
  es.addEventListener("error", (ev) => {
    if (ev.data) { es.close(); if (S?.anim === A) finishAnim(A, parse(ev).error?.message || "failed"); return; }
    es.close();
    if (S?.anim === A && A.phase === "computing") pollAnim(A);
  });
}

async function pollAnim(A) {
  while (S?.anim === A && A.phase === "computing") {
    try {
      const snap = await S.ctx.api.getJSON(`/api/jobs/${A.jobId}`);
      if (snap.state === "done") { ingestResult(A, snap.result); finishAnim(A); return; }
      if (snap.state === "error" || snap.state === "cancelled") { finishAnim(A, snap.error?.message || snap.state); return; }
    } catch (e) { if (e.offline) { finishAnim(A, e.message); return; } }
    await new Promise((r) => { A.poll = setTimeout(r, 600); });
  }
}

function ingestFrame(A, i, f) {
  if (!f || i == null || i < 0 || i >= A.n) return;
  A.frames[i] = f;
  syncRunaway(A, A.serverRunaway);
  A.tl.setComputed(i, f.error ? "error" : "computed");
  if (!A.started && !f.error) { A.started = true; A.tl.seek(i, "api"); return; }
  if (!A.tl.playing) animApply(A.tl.pos, A.tl.where());
}

function ingestResult(A, res) {
  for (const f of res?.frames || []) if (f && !A.frames[f.index]) ingestFrame(A, f.index, f);
  if (isNum(res?.runaway_at)) { A.serverRunaway = res.runaway_at; syncRunaway(A, res.runaway_at); }
  A.density = res?.density || null;
  renderDensity(A);
}

/** Rule 6 notice: frames too coarse near T_c (or for the Stark bias) are said so, with a one-click refine. */
function renderDensity(A) {
  const host = A.dom.density;
  if (!host) return;
  const d = A.density;
  clear(host);
  host.hidden = true;
  if (!d || d.rule == null) return;
  if (d.rule === "arrhenius" && d.checked && !d.ok) {
    host.hidden = false;
    host.append(h("p.anim-dens-t", `Frames are ${fmt(d.max_step)} K apart within ±20 K of T_c (${fmt(d.T_c)} K); Arrhenius quantities need no more than ${fmt(d.required_step)} K there, so read the values near T_c as coarse samples, not a curve.`),
      h("button.btn.btn-outline.anim-refine", { type: "button", on: { click: () => refineAround(A, d.suggest) } }, icon("run", { size: 16 }), `Recompute ${fmt(d.suggest.lo)} to ${fmt(d.suggest.hi)} K (${d.suggest.n} frames)`));
  } else if (d.rule === "arrhenius" && !d.checked) {
    host.hidden = false;
    host.append(h("p.anim-dens-t", `Frame density near T_c is not checked: ${d.reason}.`));
  } else if (d.rule === "stark" && !d.ok) {
    host.hidden = false;
    host.append(h("p.anim-dens-t", `Frames are ${fmt(d.max_step)} V apart; the nitride Stark shift needs one frame per ${fmt(d.required_step)} V, and E_X is never interpolated across a flat-band or regime change.`));
  }
}

function refineAround(A, sg) {
  if (!sg || S?.anim !== A) return;
  const [b0, b1] = A.meta?.band || [-Infinity, Infinity];
  A.dom.lo.value = String(Math.max(b0, sg.lo)); // the suggestion never leaves the validity band
  A.dom.hi.value = String(Math.min(b1, sg.hi));
  A.dom.n.value = String(sg.n);
  A.dom.n.dispatchEvent(new Event("input"));
  startAnim();
}

function finishAnim(A, problem = null) {
  A.phase = "ready";
  A.jobOver = true;
  A.dom.sheet.dataset.phase = "ready";
  A.dom.go.disabled = false;
  A.dom.go.lastChild.textContent = "Recompute";
  A.dom.cancel.hidden = true;
  A.tl?.setJobOver(true);
  const c = A.frames.filter((f) => f && !f.error).length;
  const e = A.frames.filter((f) => f?.error).length;
  A.dom.status.textContent = problem
    ? `Stopped: ${problem}. ${c} of ${A.n} frames computed; the rest stay hatched.`
    : `${c} of ${A.n} frames computed${e ? `; ${e} failed (${A.frames.find((f) => f?.error)?.error?.message || "error"})` : ""}. Play, or drag the scrubber.`;
  A.dom.sheet.dataset.computed = String(c);
}

function stopAnimJob(A) {
  A.es?.close();
  A.es = null;
  clearTimeout(A.poll);
  if (A.jobId && A.phase === "computing") S.ctx.api.del(`/api/jobs/${A.jobId}`).catch(() => {});
}

function cancelAnimJob() {
  const A = S?.anim;
  if (!A) return;
  stopAnimJob(A);
  finishAnim(A, "cancelled");
}

function buildPlayer(A) {
  const unit = unitText(A.meta.unit);
  A.tl = createTimeline({
    n: A.n, values: A.values, label: paramLabel(A),
    valueText: (v) => `${fmt(v)}${unit ? ` ${unit}` : ""}`,
    onChange: (pos, w) => animApply(pos, w),
    onSettle: (pos, w) => animSettle(w),
    canInterp: (i, j) => interpRefusal(A, i, j),
  });
  A.dom.cross = h("p.anim-cross", { "aria-live": "polite" });
  A.dom.density = h("div.anim-density", { hidden: true, "aria-live": "polite" });
  clear(A.dom.player).append(A.tl.el, A.dom.cross, A.dom.labels, A.dom.density);
  A.dom.player.hidden = false;
  // hero chart playhead marker: a DOM overlay (no Plotly call per frame)
  const screen = S.charts.hero?.el.querySelector(".chart-screen");
  if (screen) {
    A.dom.marker = h("span.anim-marker", { hidden: true, "aria-hidden": "true" });
    screen.append(A.dom.marker);
  }
  S.el.querySelector(".designer")?.classList.add("is-animating");
  fitHeroRange(A);
}

/** T_hs animations: widen the hero x axis once so every frame's operating point is on it. */
function fitHeroRange(A) {
  const gd = S.charts.hero?.plot;
  const r = gd?._fullLayout?.xaxis?.range;
  if (A.param !== "thermal.T_hs" || !r || !window.Plotly) return;
  const lo = Math.min(r[0], A.lo, A.hi), hi = Math.max(r[1], A.lo, A.hi);
  if (lo < r[0] || hi > r[1]) {
    window.Plotly.relayout(gd, { "xaxis.range": [lo, hi] }).then(() => { if (S?.anim === A && A.tl) animApply(A.tl.pos, A.tl.where()); }).catch(() => {});
  }
}

function valueAtPos(A, pos) {
  const i = Math.floor(pos), f = pos - i;
  const vs = A.values;
  return i >= A.n - 1 ? vs[A.n - 1] : vs[i] + (vs[i + 1] - vs[i]) * f;
}

/** Display blend of two adjacent computed frames. Everything outside ANIM_INTERP (T_c, booleans, text,
 * __nan_reason) is taken whole from the NEAREST frame; `keys` lists what was really blended. */
function interpScalars(a, b, t) {
  const out = { ...(t < 0.5 ? a : b).scalars };
  const keys = new Set();
  for (const k of ANIM_INTERP) {
    if (!(k in a.scalars) && !(k in b.scalars)) continue;
    const va = a.scalars[k], vb = b.scalars[k];
    delete out[`${k}__nan_reason`];
    if (isNum(va) && isNum(vb)) { out[k] = va + (vb - va) * t; keys.add(k); }
    else {
      out[k] = null;
      const why = a.scalars[`${k}__nan_reason`] || b.scalars[`${k}__nan_reason`] || (isNum(va) || isNum(vb) ? "not finite at one end of this interval, so not interpolated" : null);
      if (why) out[`${k}__nan_reason`] = why;
    }
  }
  return { scalars: out, keys };
}

function frameView(A, f, scalars, tag) {
  return { mode: A.mode, scalars, tag_chain: tag ?? f.tag_chain, labels: f.labels || [], provenance: f.provenance || {}, curves: { T_hs: f.curves?.T_hs || [] } };
}

function crossNote(A, i, j) {
  // a gate crossing (g2 ceiling, flux floor) or a T_c change between two computed frames is reported
  // as "between frames i and j", never located inside the interval (rules 3 and 4)
  const g = S.ctx.api.gatesNow?.() || {};
  const gv = g.g2_ceiling?.value, ff = g.flux_floor?.value;
  const a = A.frames[i]?.scalars, b = A.frames[j]?.scalars;
  if (!a || !b) return "";
  const out = [];
  if (isNum(gv) && isNum(a.g2_op) && isNum(b.g2_op) && (a.g2_op - gv) * (b.g2_op - gv) < 0) out.push(`g²(0) ceiling crossed between frames ${i + 1} and ${j + 1}`);
  const fl = (s) => (isNum(s.collected_flux_delivered_s) ? s.collected_flux_delivered_s : s.collected_flux_pulsed_s);
  if (isNum(ff) && isNum(fl(a)) && isNum(fl(b)) && (fl(a) - ff) * (fl(b) - ff) < 0) out.push(`flux floor crossed between frames ${i + 1} and ${j + 1}`);
  if (isNum(a.T_c) !== isNum(b.T_c) || (isNum(a.T_c) && Math.abs(a.T_c - b.T_c) > 1e-9 * Math.max(1, Math.abs(a.T_c)))) out.push(`T_c changes between frames ${i + 1} and ${j + 1}`);
  return out.length ? `${out.join("; ")} (placed in the interval only, not inside it)` : "";
}

function animApply(pos, w) {
  const A = S?.anim;
  if (!A?.tl || !w) return;
  const fr = A.frames;
  const v = valueAtPos(A, pos);
  const lab = paramLabel(A);
  const slice = A.meta?.tag === "A" ? `; ${A_SLICE_NOTE}` : "";
  const near = w.near ?? w.held;
  const interp = w.kind === "interp";
  const snap = w.kind === "snap";
  const pair = interp ? [w.held, w.next] : snap ? [w.a, w.b] : null;
  // the crossing text depends only on the pair of frames: computed once per pair, not per screen frame
  let crossing = "";
  if (pair) {
    const ck = `${pair[0]}:${pair[1]}:${A.patchRev || 0}`;
    if (A.crossKey !== ck) { A.crossKey = ck; A.crossText = crossNote(A, pair[0], pair[1]); }
    crossing = A.crossText;
  }
  let res = null, anim = null;
  if (w.kind !== "gap" && near != null && (frameRunaway(A, near) || fr[near]?.scalars?.runaway === true)) {
    const r = A.runawayIdx ?? near;
    anim = { runaway: `thermal runaway at frame ${r + 1} (${lab} = ${paramText(A, A.values[r])}): the displayable range ends here; nothing is interpolated through it` };
  } else if (w.kind === "frame") {
    res = frameView(A, fr[w.held], fr[w.held].scalars);
    anim = { note: `${lab} = ${paramText(A, v)}; computed frame ${w.held + 1} of ${A.n}${slice}`, tcNote: `T_c from computed frame ${w.held + 1}; never interpolated` };
  } else if (interp) {
    const ip = interpScalars(fr[w.held], fr[w.next], w.f);
    res = frameView(A, fr[near], ip.scalars);
    anim = { interpKeys: ip.keys, note: `${lab} = ${paramText(A, v)}; interp. between computed frames ${w.held + 1} and ${w.next + 1} (display only)${slice}`, tcNote: `snapped to the nearest computed frame ${near + 1}; T_c is never interpolated` };
  } else if (snap) {
    res = frameView(A, fr[near], fr[near].scalars);
    anim = { note: `${lab} = ${paramText(A, A.values[near] ?? v)}; snapped to computed frame ${near + 1}, not interpolated (${w.why})${slice}`, tcNote: `snapped to the nearest computed frame ${near + 1}; T_c is never interpolated` };
  } else {
    anim = { gap: `frame ${(w.at ?? 0) + 1} (${lab} = ${paramText(A, v)}) is not computed yet; nothing is interpolated across it` };
  }
  renderReadouts(res, { markChange: false, anim });
  if (A.dom.cross && A.crossShown !== crossing) { A.crossShown = crossing; A.dom.cross.textContent = crossing; }
  // beam: follows the flux at the playhead (no transition while animating)
  const { a, flux } = beamAlpha(res?.scalars);
  const svg = S.dom.beamSvg;
  const beam = svg.querySelector(".beam");
  if (beam) { beam.style.opacity = a.toFixed(3); beam.dataset.alpha = a.toFixed(3); }
  const halo = svg.querySelector(".beam-halo");
  if (halo) halo.style.opacity = (a * 0.4).toFixed(3);
  const src = svg.querySelector(".beam-src");
  if (src) src.style.opacity = a.toFixed(3);
  S.dom.beamCap.textContent = flux != null ? `Beam follows the playhead: log₁₀ collected flux ${fmt(flux)} s⁻¹${interp ? " (interp.)" : ""}` : "Beam at rest level: no collected flux at the playhead";
  // 3D: the cheap glow follow (no rebuild, no new WebGL context); the geometry rebuild waits for rest
  const glow = Math.max(0, Math.min(1, (a - 0.45) / 0.55));
  if (S.scene?.setPlayhead && (A.lastGlow == null || Math.abs(glow - A.lastGlow) > 0.004)) {
    A.lastGlow = glow;
    S.scene.setPlayhead({ T: A.param === "thermal.T_hs" ? (snap ? (A.values[near] ?? v) : v) : undefined, glow });
  }
  // hero playhead marker: on a snap it sits on the frame it snapped to, never between frames
  const mk = A.dom.marker;
  const gd = S.charts.hero?.plot;
  if (mk && gd) {
    const T = A.param === "thermal.T_hs" ? (snap ? A.values[near] : v) : getPath(S.design, "thermal.T_hs");
    const xy = res ? playheadXY(gd, T, res.scalars?.g2_op) : null;
    if (!xy) mk.hidden = true;
    else {
      mk.hidden = false;
      mk.style.transform = `translate(${(gd.offsetLeft + xy.x).toFixed(1)}px, ${(gd.offsetTop + xy.y).toFixed(1)}px)`;
      mk.dataset.kind = w.kind;
      mk.dataset.out = xy.inX && xy.inY && xy.hasY ? "0" : "1";
    }
  }
  // labels and booleans: only from the NEAREST computed frame (rule 3)
  const shown = w.kind === "gap" ? w.held : near;
  if (A.heldShown !== shown) { A.heldShown = shown; renderAnimLabels(A, shown); }
  A.dom.sheet.dataset.kind = w.kind;
}

function renderAnimLabels(A, idx) {
  const host = A.dom.labels;
  host.dataset.frame = idx == null ? "" : String(idx);
  const f = idx == null ? null : A.frames[idx];
  if (!f) { clear(host).append(h("p.anim-lab-none", "No computed frame at or before the playhead.")); return; }
  const flags = [];
  if (f.scalars?.runaway === true) flags.push("runaway: thermal runaway at this frame");
  if (f.scalars?.g2_op_valid === false) flags.push(`g²(0) not valid: ${f.scalars.g2_op_invalid_reason || "see the run"}`);
  clear(host).append(
    h("span.engrave.anim-lab-k", `Frame ${idx + 1} labels (nearest computed frame)`),
    h("span.anim-lab-list", (f.labels || []).concat(flags).map((l) => h("span.label-chip", l))),
    h("span.chip.chip-mono.anim-lab-id", { title: "catalog ID of this frame's run (same cache entry as /api/run)" }, f.run_id || ""));
}

/** The playhead came to rest: the 3D scene follows, at the computed frame it rests on or after. */
function animSettle(w) {
  const A = S?.anim;
  const idx = w?.near ?? w?.held;
  if (!A || A.param !== "thermal.T_hs" || S.draft || !S.scene || idx == null || frameRunaway(A, idx)) return;
  const T = A.frames[idx]?.value;
  if (!isNum(T) || T === A.sceneT) return;
  A.sceneT = T;
  const note = S.el.querySelector("#vm-note");
  if (note) note.textContent = `3D at computed frame ${idx + 1}: T_hs ${fmt(T)} K (follows the playhead when it rests)`;
  loadScene(T);
}

function closeAnim({ quiet = false } = {}) {
  const A = S?.anim;
  if (!A) return;
  stopAnimJob(A);
  A.tl?.dispose();
  A.dom.marker?.remove();
  S.scene?.setGlow?.(1);
  S.anim = null;
  if (S.dom?.animHost) { clear(S.dom.animHost); S.dom.animHost.classList.remove("is-open"); }
  S.dom?.hsAnim?.setAttribute("aria-expanded", "false");
  S.el?.querySelector(".designer")?.classList.remove("is-animating");
  if (quiet) return;
  if (A.sceneT != null) {
    const note = S.el.querySelector("#vm-note");
    if (note) note.textContent = "true scale; exaggerations are badged";
    loadScene(getPath(S.design, "thermal.T_hs"));
  }
  renderResults();
}
