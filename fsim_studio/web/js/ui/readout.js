// Instrument readout: engraved label, tabular value (or lo – hi interval), unit, tag chip,
// optional caption. A changed value holds a dashed beam underline until it is clicked
// (gate-board raise). A number without a tag does not render (CONTRACT rule 2).
import { h, clear } from "./dom.js";
import { tagChip } from "./tagChip.js";
import { fmt, fmtInterval, isNum } from "./format.js";
import { popover } from "./controls.js";

// ---------------------------------------------------------------- symbols
// Physics symbols are never engraved: they sit in .sym (text-transform: none) with real notation,
// so g^{(2)}(0) reads g(2)(0) with a superscript and F_{P} keeps its subscript.
const GREEK = {
  alpha: "α", beta: "β", gamma: "γ", Gamma: "Γ", delta: "δ", Delta: "Δ", epsilon: "ε", eta: "η", theta: "θ",
  kappa: "κ", lambda: "λ", Lambda: "Λ", mu: "μ", nu: "ν", rho: "ρ", sigma: "σ", tau: "τ", chi: "χ", phi: "φ",
  omega: "ω", Omega: "Ω", pi: "π",
};
const greekize = (t) => String(t).replace(/[A-Za-z]+/g, (m) => GREEK[m] || m);

/** sym("F_{P}·L") -> <span class="sym">F<sub>P</sub>·L</span>. `_x`, `_{xyz}`, `^x`, `^{(2)}`; Greek names become letters. */
export function sym(markup) {
  const out = h("span.sym");
  const re = /([_^])(\{[^}]*\}|[A-Za-z0-9]+)|([^_^]+)/g;
  let m;
  while ((m = re.exec(String(markup)))) {
    if (m[3] !== undefined) out.append(greekize(m[3]));
    else out.append(h(m[1] === "_" ? "sub" : "sup", greekize(m[2].replace(/^\{|\}$/g, ""))));
  }
  return out;
}

const unitKey = (u) => String(u || "").toLowerCase().replace(/µ/g, "u").replace(/[^a-z0-9]/g, "");

/**
 * Display name of a registry parameter or scene scalar from its snake_case id.
 * Drops a trailing _<unit> that matches `unit` (the unit is shown beside the value anyway).
 * Returns {text, markup}: `markup` is sym() notation for symbols (N_top -> N_{top}, lambda_nm -> λ),
 * null when the name is a plain word (fwhm_ps -> "fwhm"); `text` is the plain-text fallback.
 */
export function paramName(name, unit = "") {
  const parts = String(name).split("_").filter(Boolean);
  const u = unitKey(unit);
  if (parts.length > 1 && u && unitKey(parts[parts.length - 1]) === u) parts.pop();
  const base = parts[0] || String(name);
  const rest = parts.slice(1).join(" ");
  const text = parts.join(" ");
  const symbolic = base.length <= 2 || !!GREEK[base];
  if (!symbolic) return { text, markup: null };
  return { text, markup: rest ? `${base}_{${rest}}` : base };
}

/** DOM node for a parameter label: a symbol in .sym, or the plain words. */
export function paramLabelNode(name, unit = "") {
  const n = paramName(name, unit);
  return n.markup ? sym(n.markup) : document.createTextNode(n.text);
}

export function createReadout({ key, label, labelText, onInspect }) {
  const name = labelText || (typeof label === "string" ? label : label?.textContent || key);
  const valueEl = h("span.ro-value");
  const unitEl = h("span.ro-unit");
  const tagSlot = h("span.ro-tag");
  const capEl = h("span.ro-caption");
  const titleEl = h("span.ro-label", label);
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
    popover(whyBtn, h("div.why", h("h2.why-title", name), h("p.why-lede", whyText)), { label: `Why ${name}` });
  }
  let last = null;
  let lastText = null;

  function update(spec) {
    // spec: {value, lo, hi, unit, tag, caption, reason, qualifier, pending}
    last = spec;
    root.classList.toggle("is-pending", !!spec?.pending);
    clear(tagSlot);
    qualEl.textContent = spec?.qualifier || "";
    capEl.textContent = spec?.caption || "";
    whyText = spec?.why || "";
    whyBtn.hidden = !whyText;
    // marker: a backend label pinned to this reading (e.g. the b_res structural floor)
    clear(markEl);
    markEl.hidden = !spec?.marker;
    if (spec?.marker) {
      markEl.append(h("span.ro-marker-rule", { "aria-hidden": "true" }), h("span.ro-marker-t", spec.marker.text));
      const mc = tagChip(spec.marker.tag, { size: "sm", title: spec.marker.title || null });
      if (mc) markEl.append(mc);
      markEl.title = spec.marker.title || spec.marker.text;
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
    if (text == null) {
      root.classList.add("is-nan");
      valueEl.textContent = "n/a";
      unitEl.textContent = "";
      capEl.textContent = spec.reason || spec.caption || "not computed for this design";
      root.dataset.state = "nan";
    } else {
      root.classList.remove("is-nan");
      valueEl.textContent = text;
      unitEl.textContent = spec.unit || "";
      root.dataset.state = "value";
    }
    if (chip) tagSlot.appendChild(chip);
    valueBtn.setAttribute("aria-label",
      `${name}: ${valueEl.textContent} ${unitEl.textContent}${spec.tag ? `, tag ${spec.tag}` : ""}`);
    const sig = `${valueEl.textContent}|${unitEl.textContent}`;
    if (lastText != null && sig !== lastText && spec.markChange !== false) root.classList.add("is-changed");
    lastText = sig;
  }

  function setLabel(text) { clear(titleEl); titleEl.append(text); }

  return { el: root, update, setLabel, clearChanged: () => root.classList.remove("is-changed") };
}
