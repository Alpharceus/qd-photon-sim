// Provenance tag chip. Graded by LINE FORM, never hue:
// [V] solid chip, [DR] solid with an inset rule, [E] outline, [A] dashed outline.
import { h } from "./dom.js";

export const TAG_MEANING = {
  V: "verified against the cited paper",
  DR: "derived from verified inputs",
  E: "estimate or class range",
  A: "assumption",
};
const ORDER = ["V", "DR", "E", "A"];

/** "[A]" | "A" | "[DR]" -> "DR"; anything else -> null */
export function normTag(t) {
  if (t == null) return null;
  const s = String(t).replace(/[[\]\s]/g, "").toUpperCase();
  return ORDER.includes(s) ? s : null;
}

/** Widest tag wins (V < DR < E < A). */
export function widest(tags) {
  let best = -1;
  for (const t of tags) { const i = ORDER.indexOf(normTag(t)); if (i > best) best = i; }
  return best >= 0 ? ORDER[best] : null;
}

export function tagChip(tag, { title = null, onClick = null, size = "" } = {}) {
  const t = normTag(tag);
  if (!t) return null;
  const label = `Provenance [${t}]: ${TAG_MEANING[t]}${title ? `. ${title}` : ""}`;
  const el = h(onClick ? "button.tag" : "span.tag", {
    class: `tag-${t.toLowerCase()} ${size ? `tag-${size}` : ""}`,
    title: label,
    "aria-label": label,
    type: onClick ? "button" : null,
    on: onClick ? { click: onClick } : null,
    dataset: { tag: t },
  }, t);
  return el;
}

/** Plotly line dash for a tag (same grammar as the chips): [A] is dashed, never dotted. */
export function tagDash(tag) {
  return { V: "solid", DR: "solid", E: "longdash", A: "dash" }[normTag(tag)] || "solid";
}

/** Presentation scale for chart strokes: wider in Present and Story mode (projector). */
export function strokeScale() {
  const d = document.documentElement;
  return d.hasAttribute("data-present") || d.dataset.story === "1" ? 1.4 : 1;
}

/** Plotly line width for a model curve with this tag: dashed [A] at >= 2.5px so the dash reads. */
export function tagWidth(tag, base = 2) {
  const t = normTag(tag);
  const w = t === "A" || t === "E" ? Math.max(2.5, base) : base;
  return Math.round(w * strokeScale() * 100) / 100;
}
