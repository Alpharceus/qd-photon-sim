// One parameter row in a block editor, generated from META: label, value input or slider,
// unit, tag chip, and (for [A]/[E] numeric fields) a range toggle that turns the slider into a
// bracketed [lo, hi] span. Out-of-band values get a dashed beam underline and a sentence.
// Validation here is a band comparison of the input only (META lo/hi), never physics.
import { h } from "./dom.js";
import { icon } from "./icons.js";
import { tagChip, normTag } from "./tagChip.js";
import { fmt, unitText, isNum } from "./format.js";

const SLIDER_STEPS = 400;

function sliderable(band) {
  if (!Array.isArray(band) || !isNum(band[0]) || !isNum(band[1]) || band[1] <= band[0]) return false;
  const floor = Math.max(Math.abs(band[0]), Math.abs(band[1]) * 1e-6);
  return band[1] / floor <= 1e4;
}

function niceStep(band) {
  return (band[1] - band[0]) / SLIDER_STEPS;
}

export function fieldRow({ path, meta, value, ranged, envDefault, onChange, onRange, onInspect, signature = false }) {
  const tag = normTag(meta.tag);
  const band = meta.band;
  const unit = unitText(meta.unit === "enum" ? "" : meta.unit);
  const id = `f-${path.replace(/\W/g, "-")}`;
  const label = h("label.fr-label", { for: id }, meta.label || path.split(".").pop());
  const pathEl = h("span.fr-path", path);
  const chip = tagChip(tag, { title: meta.source, onClick: () => onInspect?.(path) });
  const warn = h("p.fr-warn", { hidden: true, id: `${id}-warn` });
  const row = h("div.field-row", { dataset: { path, tag: tag || "" } });
  let control;
  let rangeBtn = null;

  // ----- enum / bool / structured
  if (Array.isArray(meta.choices) || meta.unit === "enum") {
    const choices = meta.choices || (path === "drive.mode" ? ["EL", "PL"] : [String(value)]);
    control = h("select.fr-select", { id, on: { change: (e) => onChange(e.target.value) } },
      choices.map((c) => h("option", { value: c, selected: String(c) === String(value) }, c)));
  } else if (typeof value === "boolean") {
    const input = h("input", { id, type: "checkbox", role: "switch", checked: value, on: { change: (e) => onChange(e.target.checked) } });
    control = h("label.switch.fr-switch", input, h("span.switch-track", { "aria-hidden": "true" }), h("span.switch-label", value ? "on" : "off"));
    input.addEventListener("change", () => { control.querySelector(".switch-label").textContent = input.checked ? "on" : "off"; });
  } else if (value !== null && typeof value === "object") {
    control = h("span.fr-struct", Object.keys(value).length ? `${Object.keys(value).length} entries; edit in the card YAML` : "empty; edit in the card YAML");
  } else if (typeof value === "string") {
    control = h("input.fr-num", { id, type: "text", value, on: { change: (e) => onChange(e.target.value) } });
  } else {
    // ----- numeric: slider (point) or bracket (range)
    const canSlide = sliderable(band);
    const canRange = (tag === "A" || tag === "E") && (value === null || isNum(value));
    const num = h("input.fr-num", {
      id, type: "number", step: "any", value: isNum(value) ? String(value) : "",
      placeholder: value === null ? "auto" : "",
      "aria-describedby": `${id}-warn`,
      on: { change: (e) => { const v = e.target.value === "" ? null : Number(e.target.value); if (v === null || Number.isFinite(v)) { onChange(v); syncSlider(v); checkBand(v); } } },
    });
    let slider = null;
    if (canSlide) {
      slider = h("input.fr-slider", {
        type: "range", min: band[0], max: band[1], step: niceStep(band),
        value: isNum(value) ? value : band[0],
        "aria-label": `${meta.label || path} slider`,
        on: {
          input: (e) => { const v = Number(e.target.value); num.value = String(Number(v.toPrecision(5))); checkBand(v); if (signature) onChange(Number(v.toPrecision(6)), { live: true }); },
          change: (e) => { const v = Number(Number(e.target.value).toPrecision(6)); onChange(v); },
        },
      });
    }
    function syncSlider(v) { if (slider && isNum(v)) slider.value = v; }

    // range (bracket) editor
    const lo0 = ranged?.[0] ?? envDefault?.[0] ?? (canSlide ? band[0] : value);
    const hi0 = ranged?.[1] ?? envDefault?.[1] ?? (canSlide ? band[1] : value);
    const loIn = h("input.fr-num.fr-lo", { type: "number", step: "any", value: lo0, "aria-label": `${meta.label || path} range low` });
    const hiIn = h("input.fr-num.fr-hi", { type: "number", step: "any", value: hi0, "aria-label": `${meta.label || path} range high` });
    const span = h("div.bracket", { "aria-hidden": "true" }, h("span.bracket-fill"));
    const placeSpan = () => {
      if (!canSlide) { span.hidden = true; return; }
      const a = (Number(loIn.value) - band[0]) / (band[1] - band[0]);
      const b = (Number(hiIn.value) - band[0]) / (band[1] - band[0]);
      const l = Math.max(0, Math.min(a, b)), r = Math.min(1, Math.max(a, b));
      span.firstChild.style.left = `${l * 100}%`;
      span.firstChild.style.width = `${Math.max(0.5, (r - l) * 100)}%`;
    };
    const emitRange = () => {
      const lo = Number(loIn.value), hi = Number(hiIn.value);
      if (Number.isFinite(lo) && Number.isFinite(hi)) onRange?.([Math.min(lo, hi), Math.max(lo, hi)]);
      placeSpan();
    };
    loIn.addEventListener("change", emitRange);
    hiIn.addEventListener("change", emitRange);
    const rangeBox = h("div.fr-range", h("span.br", "["), loIn, h("span.br-sep", ","), hiIn, h("span.br", "]"));

    const pointBox = h("div.fr-point", slider, num);
    const rangeWrap = h("div.fr-rangewrap", span, rangeBox);
    control = h("div.fr-control", pointBox, rangeWrap);

    const setRanged = (on) => {
      row.classList.toggle("is-ranged", on);
      pointBox.hidden = on;
      rangeWrap.hidden = !on;
      rangeBtn?.setAttribute("aria-pressed", on ? "true" : "false");
      if (on) placeSpan();
    };
    if (canRange) {
      rangeBtn = h("button.icon-btn.fr-rangebtn", {
        type: "button", "aria-pressed": "false",
        title: "Run this input as a range (envelope)", "aria-label": `Range ${meta.label || path} for envelope runs`,
        on: {
          click: () => {
            const on = rangeBtn.getAttribute("aria-pressed") !== "true";
            setRanged(on);
            onRange?.(on ? [Math.min(Number(loIn.value), Number(hiIn.value)), Math.max(Number(loIn.value), Number(hiIn.value))] : null);
          },
        },
      }, icon("range", { size: 18 }));
    }
    setRanged(!!ranged && canRange);
    checkBand(value);
  }

  function checkBand(v) {
    const out = isNum(v) && Array.isArray(band) && isNum(band[0]) && isNum(band[1]) && (v < band[0] || v > band[1]);
    row.classList.toggle("is-out", !!out);
    warn.hidden = !out;
    if (out) warn.textContent = `Outside the known-physical band ${fmt(band[0])} – ${fmt(band[1])}${unit ? ` ${unit}` : ""}; the result is a conditional number.`;
  }

  row.append(
    h("div.fr-name", label, pathEl),
    h("div.fr-body", control, h("span.fr-unit", unit)),
    h("div.fr-tools", chip, rangeBtn),
    warn,
  );
  row.addEventListener("focusin", () => onInspect?.(path, { quiet: true }));
  return row;
}
