// Small readouts for baked / file numbers (story rail, overview strips, drawers), the campaign
// provenance tag lookup, and the figure lightbox.
// A number without a tag does not render (CONTRACT rule 2): sro() withholds it.
import { h, clear } from "./dom.js";
import { icon } from "./icons.js";
import { tagChip, widest, normTag } from "./tagChip.js";
import { fmt, unitText, isNum } from "./format.js";
import { getJSON } from "../api.js";

/** Display text for a value that may be a number, list, bool or string (raw-faithful; rounding only). */
export function valueText(v, raw) {
  if (Array.isArray(v)) return v.map((x) => (isNum(x) ? fmt(x) : String(x))).join(", ");
  if (isNum(v)) return fmt(v);
  if (typeof v === "boolean") return v ? "true" : "false";
  if (v == null) return raw ?? null;
  return String(v);
}

/**
 * Small readout: {label, value, raw?, unit?, tag, note? (tag-chip tooltip), caption? (visible line), error?, size?, title?}
 * Renders "n/a + reason" for null, the error text for a failed resolution, and refuses a value without a tag.
 */
export function sro({ label, value, raw, unit = "", tag, note = "", caption = "", error = null, size = "", title = null, key = null }) {
  const t = normTag(tag);
  const chip = tagChip(t, { title: note || null, size: size === "lg" ? "" : "sm" });
  let text = error ? null : valueText(value, raw);
  let state = "value";
  let cap = caption;
  if (error) { state = "error"; cap = `did not resolve: ${error}`; }
  else if (text == null || text === "") { state = "nan"; text = "n/a"; }
  else if (!chip) { state = "nan"; cap = "no provenance tag; number withheld"; text = "n/a"; }
  const u = unitText(unit === "photons/s" ? "/s" : unit);
  return h("div.sro", { class: size ? `sro-${size}` : "", dataset: { state, key: key || label }, title: title || (raw != null ? `raw: ${raw}` : null) },
    h("div.sro-head", h("span.sro-label", label), chip),
    h("div.sro-reading", h("span.sro-value", state === "error" ? "n/a" : text), u && state === "value" ? h("span.sro-unit", u) : null),
    cap ? h("p.sro-cap", cap) : null);
}

// ------------------------------------------------------------------ campaign provenance tag
const TAG_RE = /\[((?:V|DR|E|A)(?:\/(?:V|DR|E|A))*)\]/g;
const tagMemo = new Map();

/**
 * Widest provenance tag of a campaign's sweep rows, read from the CSV itself: the bracket tags in
 * its provenance column, or [A] when its assumptions column lists assumed inputs (rt_edge).
 * Resolves {tag, note} (tag null when the CSV carries no provenance).
 */
export function campaignTag(id, columns = null) {
  if (tagMemo.has(id)) return tagMemo.get(id);
  const p = (async () => {
    let cols = columns;
    if (!cols) { try { cols = (await getJSON(`/api/campaigns/${id}`)).columns || []; } catch { cols = []; } }
    const names = cols.map((c) => c.name || c);
    const col = names.includes("provenance") ? "provenance" : names.includes("assumptions") ? "assumptions" : null;
    if (!col) return { tag: null, note: "this CSV carries no provenance column" };
    const r = await getJSON(`/api/campaigns/${id}/sweep?cols=${col}&limit=200`);
    const vals = r.rows.map((x) => x[0]).filter((x) => x != null && x !== "");
    if (col === "assumptions") {
      return vals.length ? { tag: "A", note: `assumptions column lists assumed inputs on ${vals.length}/${r.rows.length} sampled rows` } : { tag: null, note: "assumptions column empty" };
    }
    const found = [];
    for (const v of vals) for (const m of String(v).matchAll(TAG_RE)) found.push(...m[1].split("/"));
    const t = widest(found);
    return { tag: t, note: t ? `widest tag in the provenance column (${[...new Set(found)].join(", ")} present)` : "no bracket tags in the provenance column" };
  })().catch(() => ({ tag: null, note: "provenance column not readable" }));
  tagMemo.set(id, p);
  return p;
}

// ------------------------------------------------------------------ lightbox
let lb = null;
/** items: [{src, caption}], start index. Esc closes, arrows step. */
export function lightbox(items, start = 0) {
  closeLightbox();
  let i = start;
  const img = h("img.lb-img", { alt: "" });
  const cap = h("p.lb-cap");
  const count = h("span.lb-count");
  const go = (d) => { i = (i + d + items.length) % items.length; show(); };
  const prev = h("button.icon-btn.lb-nav", { type: "button", "aria-label": "Previous figure", on: { click: () => go(-1) } }, icon("chevronRight", { size: 22, cls: "flip" }));
  const next = h("button.icon-btn.lb-nav", { type: "button", "aria-label": "Next figure", on: { click: () => go(1) } }, icon("chevronRight", { size: 22 }));
  const close = h("button.icon-btn.lb-close", { type: "button", "aria-label": "Close figure", on: { click: closeLightbox } }, icon("close", { size: 20 }));
  const open = h("a.btn.btn-quiet.lb-open", { target: "_blank", rel: "noopener" }, "Open file");
  const dlg = h("div.lightbox", { role: "dialog", "aria-modal": "true", "aria-label": "Figure viewer", tabindex: "-1" },
    h("div.lb-bar", count, cap, open, close),
    h("div.lb-stage", prev, h("figure.lb-fig", img), next));
  function show() {
    const it = items[i];
    img.src = it.src;
    img.alt = it.caption;
    cap.textContent = it.caption;
    count.textContent = `${i + 1} / ${items.length}`;
    open.href = it.src;
    prev.hidden = next.hidden = items.length < 2;
  }
  const onKey = (e) => {
    if (e.key === "Escape") { e.stopPropagation(); e.preventDefault(); closeLightbox(); }
    else if (e.key === "ArrowRight") { go(1); e.preventDefault(); }
    else if (e.key === "ArrowLeft") { go(-1); e.preventDefault(); }
  };
  dlg.addEventListener("click", (e) => { if (e.target === dlg || e.target.classList.contains("lb-stage")) closeLightbox(); });
  document.addEventListener("keydown", onKey, true);
  const prevFocus = document.activeElement;
  lb = { dlg, cleanup: () => { document.removeEventListener("keydown", onKey, true); prevFocus?.focus?.(); } };
  document.body.append(dlg);
  show();
  dlg.focus();
}
export function closeLightbox() { if (!lb) return; lb.cleanup(); lb.dlg.remove(); lb = null; }
export const isLightboxOpen = () => !!lb;

export { clear };
