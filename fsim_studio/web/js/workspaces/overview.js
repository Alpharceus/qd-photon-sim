// Overview workspace (landing): the verdict board. One machined plaque per campaign from
// /api/campaigns (word + reason + counts with their own denominators + generated date + commit),
// the nanowire BEST strip (commanded and delivered side by side, photons-per-cycle caveat
// inline), the figures/CSV collections with pre-audit badges, and a 3D hero of the selected
// platform's device. Verdicts are parsed VERDICT/BEST lines; nothing is recomputed here.

import { h, clear } from "../ui/dom.js";
import { icon } from "../ui/icons.js";
import { skeleton, errorState } from "../ui/controls.js";
import { verdictBadge, lineStatus, lineScope, humanWord, COUNT_KEYS, countWithDenom, conservativeLine } from "../ui/verdictBadge.js";
import { sro, campaignTag } from "../ui/figure.js";

const CAMPAIGNS = [
  { id: "rt_edge", title: "RT edge emitter", sub: "InP dots in GaInP / GaAsP, finite-pulse headline model" },
  { id: "nitride_cavity", title: "Planar InGaN cavity", sub: "c-plane dot in a DBR microcavity, pulse and SET loading" },
  { id: "nitride_geometry_stark", title: "Nitride geometry + Stark", sub: "orientation, screening scenarios, bias response" },
  { id: "nitride_nanowire/full", title: "Nitride nanowires", sub: "horizontal as-built and vertical photonic wires, strain bounds" },
];
const HERO_CARDS = [
  { card: "edge-inp-gainp-design", label: "InP edge RT" },
  { card: "nitride-cavity-set-design", label: "Planar nitride" },
  { card: "nitride-nanowire-horizontal-set-design", label: "Nanowire H" },
  { card: "nitride-nanowire-vertical-set-design", label: "Nanowire V" },
];

let O = null;

export function mount(el, ctx) {
  O = { el, ctx, scene: null, sceneSeq: 0, disposers: [], heroCard: pickHeroCard(ctx) };
  el.classList.add("ws-overview");
  O.disposers.push(ctx.state.on("theme", (m) => O?.scene?.setTheme?.(m)));
  render();
}

export function unmount() {
  if (!O) return;
  for (const d of O.disposers) d();
  try { O.scene?.dispose?.(); } catch { /* gone */ }
  O.el.classList.remove("ws-overview");
  O = null;
}

function pickHeroCard(ctx) {
  const cur = ctx.state.get("card");
  const plat = cur ? ctx.platformOf(cur) : null;
  if (plat) {
    const hit = HERO_CARDS.find((c) => ctx.platformOf(c.card)?.id === plat.id);
    if (hit) return hit.card;
  }
  return HERO_CARDS[0].card;
}

async function render() {
  const el = O.el;
  clear(el);
  const board = h("div.ov-board", { "aria-busy": "true" },
    h("div.ov-plaques", CAMPAIGNS.map(() => skeleton("plaque", { lines: 4 }))),
    skeleton("block", { lines: 2, height: 120 }));
  O.dom = { board, hero: renderHero() };
  el.append(h("div.overview",
    h("header.ws-head.ov-head",
      h("div", h("h1.ws-title", "Verdict board"),
        h("p.ws-sub", "Every word below is a VERDICT or BEST line parsed from out/*.md. The Studio never recomputes a verdict.")),
      h("div.ov-legend", { "aria-label": "Badge legend" },
        verdictBadge("FAIL", { size: "sm" }), verdictBadge("pass_hardware_infeasible", { size: "sm" }),
        h("span.ov-legend-txt", "status is always icon + word"))),
    h("div.ov-grid", board, O.dom.hero)));
  loadScene();
  let camps;
  try {
    camps = await O.ctx.api.listCampaigns();
  } catch (e) {
    if (!O) return;
    clear(board).append(errorState({ title: "Could not read out/", message: e.message, actions: [{ label: "Retry", fn: () => render() }] }));
    return;
  }
  if (!O) return;
  const tags = {};
  await Promise.all(CAMPAIGNS.map(async (c) => { tags[c.id] = await campaignTag(c.id); }));
  if (!O) return;
  board.removeAttribute("aria-busy");
  clear(board).append(...[
    h("div.ov-plaques", CAMPAIGNS.map((c) => {
      const camp = camps.find((x) => x.id === c.id);
      return camp ? campaignPlaque(c, camp, tags[c.id]) : missingPlaque(c);
    })),
    bestStrip(camps.find((x) => x.id === "nitride_nanowire/full"), tags["nitride_nanowire/full"]),
    collections(camps.filter((x) => x.kind !== "campaign"))].filter(Boolean));
}

// ---------------------------------------------------------------- plaques
// Representative line: the conservative (least favourable) line by its own parsed word, never
// the most favourable one; every line stays visible in the rail below it.
function representative(camp) {
  const vs = camp.verdicts || [];
  if (vs.length <= 1) return { v: vs[0] || null, why: "" };
  const v = conservativeLine(vs);
  return { v, why: `line ${v.line}: the least favourable of ${vs.length} lines` };
}

function headerBadge(camp) {
  const rep = representative(camp).v;
  return verdictBadge(rep || camp.status_word, { size: "lg" });
}

function campaignPlaque(c, camp, tagInfo) {
  const { v, why } = representative(camp);
  const f = v?.fields || {};
  const tag = tagInfo?.tag;
  const words = camp.status_words || {};
  const nLines = camp.verdicts?.length || 0;
  const reasonKeys = ["g2_min", "flux_max"].filter((k) => f[k] != null);
  // every count keeps its own denominator (coverage N printed on the same line) and coverage itself
  const counts = COUNT_KEYS.filter((k) => f[k] != null && !["coverage", "headline_coverage", "headline_coverage_pulsed"].includes(k)).slice(0, 4);
  const cover = ["coverage", "headline_coverage_pulsed", "headline_coverage"].find((k) => f[k] != null);
  const plaque = h("a.cplaque", {
    href: `#/results/${camp.id}`, dataset: { campaign: camp.id, status: v ? lineStatus(v).status : "info" },
    "aria-label": `${c.title}: ${camp.status_word}. Open in Results`,
  },
  h("div.cp-top",
    h("div.cp-titles", h("h2.cp-title", c.title), h("p.cp-sub", c.sub)),
    h("span.cp-go", icon("chevronRight", { size: 18 }))),
  h("div.cp-word", headerBadge(camp),
    camp.status_word === "mixed"
      ? h("p.cp-mix", Object.entries(words).map(([w, n]) => `${n} × ${humanWord(w).toLowerCase()}`).join(" · "))
      : null),
  h("div.cp-reason",
    v ? h("p.cp-scope", h("span.cp-line", `VERDICT line ${v.line}`), lineScope(v) ? ` · ${lineScope(v)}` : "", why && nLines > 1 ? ` · ${why}` : "") : h("p.cp-scope", "No VERDICT line in this campaign."),
    reasonKeys.length ? h("div.cp-kv", reasonKeys.map((k) => sro({ label: k.replace(/_/g, " "), value: Number(f[k]), raw: f[k], unit: k.startsWith("flux") ? "/s" : "", tag, note: tagInfo?.note, size: "sm" }))) : null),
  counts.length ? h("div.cp-counts", { role: "group", "aria-label": "Counts with their own denominators" },
    counts.map((k) => {
      const c = countWithDenom(f, k);
      return sro({ label: k.replace(/_/g, " "), value: c.text, raw: f[k], tag, note: tagInfo?.note, size: "sm", key: `count:${k}`, title: c.denomFrom ? `denominator: the N of ${c.denomFrom} on line ${v.line}` : null });
    }),
    cover ? sro({ label: cover.replace(/_/g, " "), value: String(f[cover]).replace(/\s*\/\s*/, " / "), raw: f[cover], tag, note: tagInfo?.note, size: "sm", key: `count:${cover}` }) : null) : null,
  nLines > 1 ? lineRail(camp) : null,
  h("footer.cp-foot",
    h("span.cp-meta", icon("folder", { size: 14 }), camp.path),
    h("span.cp-meta", `generated ${String(camp.generated || "").slice(0, 10)}`),
    h("span.id-chip.is-static", { title: camp.commit_source || "" }, h("span.id-k", "DATA"), h("span.id-v", camp.commit || "n/a")),
    camp.stale ? h("span.stale-badge", { title: camp.stale_note || "" }, icon("warn", { size: 14 }), "pre-audit") : null,
    tag ? null : h("span.cp-meta", "no provenance column")));
  return plaque;
}

function lineRail(camp) {
  const vs = camp.verdicts || [];
  return h("div.cp-rail", { role: "img", "aria-label": `${vs.length} VERDICT lines: ${Object.entries(camp.status_words || {}).map(([w, n]) => `${n} ${w}`).join(", ")}` },
    h("span.cp-rail-k", `${vs.length} lines`),
    h("span.cp-ticks", vs.map((v) => h("span.cp-tick", { dataset: { status: lineStatus(v).status }, title: `line ${v.line}: ${v.word} · ${lineScope(v)}` }))));
}

function missingPlaque(c) {
  return h("div.cplaque.is-missing", h("div.cp-top", h("div.cp-titles", h("h2.cp-title", c.title), h("p.cp-sub", "no manifest found under out/"))));
}

// ---------------------------------------------------------------- BEST strip (nanowire)
function bestStrip(camp, tagInfo) {
  if (!camp) return null;
  const best = (camp.best || []);
  if (!best.length) return null;
  const tag = tagInfo?.tag;
  return h("section.best-strip", { "aria-labelledby": "best-h" },
    h("header.best-head",
      h("h2.engrave#best-h", "Best passing flux · nanowires"),
      h("p.best-note", "BEST_PASSING_FLUX lines from out/nitride_nanowire/full/results.md. Commanded flux is idealized; delivered is what the as-built contact passes. Neither is a hardware pass.")),
    h("div.best-rows", best.map((b) => {
      const f = b.fields || {};
      return h("article.best-row", { dataset: { family: f.family } },
        h("div.best-id",
          h("span.best-fam", (f.family || "").replace(/_/g, " ")),
          h("span.best-meta", `${b.kind === "BEST_PASSING_FLUX_300K" ? "at 300 K" : "any T"} · row ${f.row_id} · T_hs ${f.T_hs} K · ${Number(f.rep_rate_hz) / 1e6} MHz · ${f.strain_bound}`),
          f.hardware_qualified === "False" ? h("span.best-hw", icon("cross", { size: 12 }), "hardware not qualified") : null),
        h("div.best-pair",
          sro({ label: "commanded", value: Number(f.commanded_flux), raw: f.commanded_flux, unit: "/s", tag, note: tagInfo?.note }),
          h("span.best-arrow", { "aria-hidden": "true" }, "→"),
          sro({ label: "delivered", value: Number(f.delivered_flux), raw: f.delivered_flux, unit: "/s", tag, note: tagInfo?.note })),
        h("div.best-ppc",
          sro({ label: "photons / cycle, commanded", value: Number(f.photons_per_cycle_commanded), raw: f.photons_per_cycle_commanded, tag, size: "sm" }),
          sro({ label: "photons / cycle, delivered", value: Number(f.photons_per_cycle_delivered), raw: f.photons_per_cycle_delivered, tag, size: "sm" }),
          b.note ? h("p.best-caveat", icon("info", { size: 14 }), b.note) : null));
    })));
}

// ---------------------------------------------------------------- collections
function collections(list) {
  if (!list.length) return null;
  return h("section.ov-coll", { "aria-labelledby": "coll-h" },
    h("h2.engrave#coll-h", "Figures + CSV collections"),
    h("ul.coll-list", list.map((c) => h("li",
      h("a.coll-item", { href: `#/results/${c.id}` },
        h("span.coll-id", c.id),
        h("span.coll-kind", c.kind === "legacy" ? "legacy F-series" : c.collection_label || "figures + CSV"),
        h("span.coll-date", String(c.generated || "").slice(0, 10)),
        h("span.coll-commit", c.commit || ""),
        c.stale ? h("span.stale-badge", { title: c.stale_note || "" }, icon("warn", { size: 14 }), "pre-audit") : null)))));
}

// ---------------------------------------------------------------- 3D hero
function renderHero() {
  const tabs = h("div.hero-tabs", { role: "tablist", "aria-label": "Device scene" },
    HERO_CARDS.map((c) => h("button.hero-tab", {
      type: "button", role: "tab", dataset: { card: c.card }, "aria-selected": String(c.card === O.heroCard),
      on: { click: () => { O.heroCard = c.card; syncTabs(); loadScene(); } },
    }, c.label)));
  const slot = h("div.hero-viewport", { role: "region", "aria-label": "3D device model" });
  const cap = h("p.hero-cap");
  O.heroDom = { tabs, slot, cap };
  return h("aside.ov-hero", { "aria-label": "Device" },
    h("header.hero-head", h("h2.engrave", "Device"), tabs),
    slot,
    h("footer.hero-foot", cap, h("a.btn.btn-outline", { href: `#/design/${O.heroCard}`, id: "hero-open" }, "Open in Designer", icon("chevronRight", { size: 16 }))));
}

function syncTabs() {
  for (const t of O.heroDom.tabs.children) t.setAttribute("aria-selected", String(t.dataset.card === O.heroCard));
  const a = document.getElementById("hero-open");
  if (a) a.href = `#/design/${O.heroCard}`;
}

async function loadScene() {
  if (!O) return;
  const seq = ++O.sceneSeq;
  const { slot, cap } = O.heroDom;
  cap.textContent = O.heroCard;
  if (!O.scene) clear(slot).append(h("div.vp-placeholder", h("div.vp-plate", h("span.engrave", "3D model loading"), h("span.vp-sub", O.heroCard))));
  let spec;
  try { spec = await O.ctx.api.sceneSpec("device", { card: O.heroCard }); } catch (e) {
    if (!O || seq !== O.sceneSeq) return;
    try { O.scene?.dispose?.(); } catch { /* */ }
    O.scene = null;
    clear(slot).append(h("div.vp-placeholder", h("div.vp-plate", h("span.engrave", "3D view unavailable"), h("span.vp-sub", e.message))));
    return;
  }
  if (!O || seq !== O.sceneSeq) return;
  cap.textContent = `${spec.title || O.heroCard} · ${spec.catalog_id || ""}`;
  if (O.scene) { try { O.scene.update(spec); return; } catch { /* remount */ } }
  try {
    const mod = await import("../viz3d/index.js");
    if (!O || seq !== O.sceneSeq) return;
    clear(slot);
    O.scene = mod.mountScene(slot, "device", spec, { theme: O.ctx.state.get("theme") });
  } catch (e) {
    console.warn("3D mount failed", e);
    clear(slot).append(h("div.vp-placeholder", h("div.vp-plate", h("span.engrave", "3D view unavailable on this machine"), h("span.vp-sub", spec.title || ""))));
  }
}

