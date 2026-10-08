// Right inspector drawer: tag, unit, band and source of the selected field or number,
// straight from /api/cards/<name> meta (design_meta.META) or the run's provenance.
import { h, clear } from "./dom.js";
import { icon } from "./icons.js";
import { tagChip, TAG_MEANING, normTag } from "./tagChip.js";
import { fmt, unitText, isNum } from "./format.js";

let host = null;
let lastFocus = null;

export function initInspector(el) {
  host = el;
  host.addEventListener("keydown", (e) => { if (e.key === "Escape") { e.stopPropagation(); closeInspector(); } });
}

export function isInspectorOpen() { return host && !host.hidden; }

export function closeInspector() {
  if (!host) return;
  host.hidden = true;
  document.documentElement.classList.remove("inspector-open");
  lastFocus?.focus?.();
}

function frame(title, sub, ...body) {
  clear(host);
  host.append(
    h("div.insp-head",
      h("div", h("h2.insp-title", title), sub ? h("p.insp-sub", sub) : null),
      h("button.icon-btn", { type: "button", "aria-label": "Close inspector", on: { click: closeInspector } }, icon("close", { size: 18 }))),
    h("div.insp-body", ...body));
}

function row(k, v) { return h("div.insp-row", h("dt", k), h("dd", v)); }

function tagBlock(tag, note) {
  const t = normTag(tag);
  if (!t) return row("Provenance", "no tag returned");
  return h("div.insp-tag", tagChip(t, { size: "lg" }), h("div", h("strong", `[${t}] ${TAG_MEANING[t]}`), note ? h("p", note) : null));
}

/** Show a design field. quiet: only refresh if already open. */
export function inspectField({ path, meta, value, inBand }, { quiet = false } = {}) {
  if (!host || (quiet && host.hidden)) return;
  const unit = unitText(meta?.unit === "enum" ? "" : meta?.unit);
  const band = meta?.band;
  frame(meta?.label || path, path,
    tagBlock(meta?.tag, meta?.source),
    h("dl.insp-dl",
      row("Value", value === null ? "auto (resolved by the model)" : `${typeof value === "number" ? fmt(value) : String(value)}${unit ? ` ${unit}` : ""}`),
      row("Unit", unit || "dimensionless"),
      row("Known-physical band", Array.isArray(band) && isNum(band[0]) ? `${fmt(band[0])} – ${fmt(band[1])}${unit ? ` ${unit}` : ""}` : "none recorded"),
      row("Band check", inBand == null ? "not applicable" : inBand ? "inside band" : "outside band: result is conditional"),
      row("Source", meta?.source || "no source recorded")),
    h("p.insp-foot", "From fsim_core.design_meta.META. Bands warn only; they never clamp a value."));
  open(quiet);
}

/** Show a result number: {label, valueText, tag, notes:[{key, tag, note}], source} */
export function inspectNumber({ label, key, valueText, tag, notes = [], source }) {
  if (!host) return;
  frame(label, key,
    tagBlock(tag, "Widest tag in the run's tag chain; the number is only as strong as its weakest input."),
    h("dl.insp-dl", row("Value", valueText), row("Source", source || "fsim_core.device.evaluate scalars")),
    notes.length ? h("div.insp-notes", h("h3.insp-h3", "Provenance entries"),
      h("ul", notes.map((n) => h("li", tagChip(n.tag), h("span.insp-note-k", n.key), h("span.insp-note-v", n.note || ""))))) : null);
  open(false);
}

function open(quiet) {
  if (!quiet) lastFocus = document.activeElement;
  host.hidden = false;
  document.documentElement.classList.add("inspector-open");
  if (!quiet) host.querySelector(".icon-btn")?.focus();
}
