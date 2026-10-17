// Instrument readout: engraved label, tabular value (or lo – hi interval), unit, tag chip,
// optional caption. A changed value holds a dashed beam underline until it is clicked
// (gate-board raise). A number without a tag does not render (CONTRACT rule 2).
import { h, clear } from "./dom.js";
import { tagChip } from "./tagChip.js";
import { fmt, fmtInterval, isNum } from "./format.js";
import { popover } from "./controls.js";
import { displayUnit } from "./units.js";
import { tex, richText, richPlain, symbolFor, plainSymbol } from "./math.js";

// ---------------------------------------------------------------- symbols
// Physics symbols are typeset by KaTeX (ui/math.js): g^{(2)}(0), F_{\mathrm P}L, \eta_{\mathrm{out}} ...
// sym() keeps the old call sites working: it accepts TeX, or the former mini-notation (greek words, F_{P}·L).
const GREEK_RE = /(^|[^\\A-Za-z])(alpha|beta|gamma|Gamma|delta|Delta|epsilon|eta|theta|kappa|lambda|Lambda|mu|nu|rho|sigma|tau|chi|phi|omega|Omega|pi)(?![A-Za-z])/g;
function legacyToTex(markup) {
  return String(markup)
    .replace(GREEK_RE, "$1\\$2")
    .replace(/_\{([A-Za-z]{2,})\}/g, "_{\\mathrm{$1}}")
    .replace(/_P(?![A-Za-z{])/g, "_{\\mathrm P}")
    .replace(/·/g, "\\,\\cdot\\,");
}

/** sym("F_{P}·L") -> typeset F_P L. Pass TeX directly where possible. */
export function sym(markup) {
  return tex(legacyToTex(markup), { cls: "sym" });
}

const unitKey = (u) => String(u || "").toLowerCase().replace(/µ/g, "u").replace(/[^a-z0-9]/g, "");

/**
 * Display name of a registry parameter or scene scalar. `meta` may be the ParamSpec (with .tex) or a bare id.
 * Returns {text, tex}: `tex` is the typeset symbol (the spec's own tex, else derived), `text` its Unicode form.
 */
export function paramName(name, unit = "", texSrc = null) {
  const meta = typeof name === "object" && name ? name : { name, tex: texSrc };
  const t = symbolFor(meta);
  return { text: plainSymbol(t), tex: t, markup: t };
}

/** DOM node for a parameter label: the typeset symbol (never the snake_case id). */
export function paramLabelNode(name, unit = "", texSrc = null) {
  const meta = typeof name === "object" && name ? name : { name, tex: texSrc };
  return tex(symbolFor(meta), { cls: "sym" });
}

/** Replace the content of `el` with a rich string (plain text and $...$ math). */
function setRich(el, str) {
  const t = str == null ? "" : String(str);
  if (el.__rich === t) return;               // unchanged (an animation frame often repeats a caption)
  el.__rich = t;
  if (!/[$_²^*×]/.test(t)) { el.textContent = t; return; }  // plain: no parsing, no nodes
  clear(el);
  el.append(richText(t));
}

export function createReadout({ key, label, labelText, onInspect }) {
  const name = labelText || (typeof label === "string" ? richPlain(label) : label?.textContent || key);
  const valueEl = h("span.ro-value");
  const unitEl = h("span.ro-unit");
  const tagSlot = h("span.ro-tag");
  const capEl = h("span.ro-caption");
  const titleEl = h("span.ro-label", typeof label === "string" ? richText(label) : label);
  const whyBtn = h("button.ro-why", { type: "button", hidden: true, "aria-haspopup": "dialog", on: { click: () => openWhy() } }, "why");
  const capRow = h("div.ro-cap-row");
  const qualEl = h("span.ro-qual");
  const markEl = h("span.ro-marker", { hidden: true });
  const valueBtn = h("button.ro-reading", {
    type: "button",
    "aria-label": `${name}: no value yet`,
    on: {
      click: () => {
        root.classList.remove("is-changed");
        onInspect?.(key, last);
      },
    },
  }, valueEl, unitEl);
  const root = h("div.readout", { dataset: { readout: key } },
    h("div.ro-head", titleEl, qualEl, tagSlot),
    valueBtn,
    markEl,
    capRow);
  capRow.append(capEl, whyBtn);
  let whyText = "";
  function openWhy() {
    if (!whyText) return;
    popover(whyBtn, h("div.why", h("h2.why-title", name), h("p.why-lede", richText(whyText))), { label: `Why ${name}` });
  }
  let last = null;
  let lastText = null;

  function update(spec) {
    // spec: {value, lo, hi, unit, tag, caption, reason, qualifier, pending}
    last = spec;
    root.classList.toggle("is-pending", !!spec?.pending);
    clear(tagSlot);
    qualEl.textContent = spec?.qualifier || "";
    setRich(capEl, spec?.caption || "");
    whyText = spec?.why || "";
    whyBtn.hidden = !whyText;
    // marker: a backend label pinned to this reading (e.g. the b_res structural floor)
    clear(markEl);
    markEl.hidden = !spec?.marker;
    if (spec?.marker) {
      markEl.append(h("span.ro-marker-rule", { "aria-hidden": "true" }), h("span.ro-marker-t", richText(spec.marker.text)));
      const mc = tagChip(spec.marker.tag, { size: "sm", title: spec.marker.title || null });
      if (mc) markEl.append(mc);
      markEl.title = richPlain(spec.marker.title || spec.marker.text);
    }
    if (!spec || spec.pending) {
      valueEl.textContent = spec?.pending ? "" : "—";
      unitEl.textContent = "";
      valueBtn.setAttribute("aria-label", `${name}: ${spec?.pending ? "running" : "not run yet"}`);
      root.classList.remove("is-nan");
      root.dataset.state = spec?.pending ? "pending" : "empty";
      return;
    }
    const chip = tagChip(spec.tag, { title: spec.tagNote });
    const hasInterval = isNum(spec.lo) || isNum(spec.hi);
    let text;
    if (hasInterval) text = fmtInterval(spec.lo, spec.hi);
    else if (isNum(spec.value)) text = fmt(spec.value);
    else text = null;
    if (text != null && !chip) {
      // Provenance missing: refuse to render the number.
      text = null;
      spec = { ...spec, reason: "no provenance tag returned; number withheld" };
    }
    let shownText = "n/a";
    if (text == null) {
      root.classList.add("is-nan");
      valueEl.textContent = "n/a";
      unitEl.textContent = "";
      setRich(capEl, spec.reason || spec.caption || "not computed for this design");
      root.dataset.state = "nan";
    } else {
      root.classList.remove("is-nan");
      setRich(valueEl, text);
      shownText = text;
      unitEl.textContent = displayUnit(spec.unit);
      root.dataset.state = "value";
    }
    if (chip) tagSlot.appendChild(chip);
    valueBtn.setAttribute("aria-label",
      `${name}: ${shownText} ${unitEl.textContent}${spec.tag ? `, tag ${spec.tag}` : ""}`);
    const sig = `${shownText}|${unitEl.textContent}`;
    if (lastText != null && sig !== lastText && spec.markChange !== false) root.classList.add("is-changed");
    lastText = sig;
  }

  function setLabel(text) { clear(titleEl); titleEl.append(typeof text === "string" ? richText(text) : text); }

  return { el: root, update, setLabel, clearChanged: () => root.classList.remove("is-changed") };
}
