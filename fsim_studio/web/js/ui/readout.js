// Instrument readout: engraved label, tabular value (or lo – hi interval), unit, tag chip,
// optional caption. A changed value holds a dashed beam underline until it is clicked
// (gate-board raise). A number without a tag does not render (CONTRACT rule 2).
import { h, clear } from "./dom.js";
import { tagChip } from "./tagChip.js";
import { fmt, fmtInterval, isNum } from "./format.js";

export function createReadout({ key, label, onInspect }) {
  const valueEl = h("span.ro-value");
  const unitEl = h("span.ro-unit");
  const tagSlot = h("span.ro-tag");
  const capEl = h("span.ro-caption");
  const titleEl = h("span.ro-label", label);
  const qualEl = h("span.ro-qual");
  const markEl = h("span.ro-marker", { hidden: true });
  const valueBtn = h("button.ro-reading", {
    type: "button",
    "aria-label": `${label}: no value yet`,
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
    capEl);
  let last = null;
  let lastText = null;

  function update(spec) {
    // spec: {value, lo, hi, unit, tag, caption, reason, qualifier, pending}
    last = spec;
    root.classList.toggle("is-pending", !!spec?.pending);
    clear(tagSlot);
    qualEl.textContent = spec?.qualifier || "";
    capEl.textContent = spec?.caption || "";
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
      valueBtn.setAttribute("aria-label", `${label}: ${spec?.pending ? "running" : "not run yet"}`);
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
      `${label}: ${valueEl.textContent} ${unitEl.textContent}${spec.tag ? `, tag ${spec.tag}` : ""}`);
    const sig = `${valueEl.textContent}|${unitEl.textContent}`;
    if (lastText != null && sig !== lastText && spec.markChange !== false) root.classList.add("is-changed");
    lastText = sig;
  }

  function setLabel(text) { titleEl.textContent = text; }

  return { el: root, update, setLabel, clearChanged: () => root.classList.remove("is-changed") };
}
