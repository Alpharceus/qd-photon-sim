// Parameter sheet body for one registry entry (devlib-p11).
// The form is generated ONLY from the registry's ParamSpec list (to_meta shape plus help, group,
// visible_if, scale): no field, label or branch is typed per kind. Rows are the Designer's
// fieldRow (band warnings, tag chip, inspector, value | range toggle); this adds the per-leaf
// tag selector, the dashed rim on edited leaves, the folded advanced group, and visible_if.
//
// Spec leaves (fsim-composed/1): a bare scalar, or a mapping with exactly one of value | range
// (+ scale) and unit, tag, source. An edit writes a mapping; "reset" deletes the leaf so the
// registry default (with its own tag) applies again.
import { h, clear } from "./dom.js";
import { fieldRow } from "./fieldRow.js";
import { tagChip, normTag } from "./tagChip.js";
import { isNum } from "./format.js";
import { paramLabelNode } from "./readout.js";
import { symbolText, richPlain } from "./math.js";
import { displayUnit } from "./units.js";

const TAGS = ["V", "DR", "E", "A"];
const EDIT_SOURCE = "set in the Composer (user value)";

/** Normalise a stored leaf against its ParamSpec. */
export function readLeaf(p, raw) {
  if (raw === undefined) return { kind: "default", value: p.default, range: null, tag: normTag(p.tag) || "A", explicit: false, source: p.source };
  if (raw !== null && typeof raw === "object" && !Array.isArray(raw)) {
    const tag = normTag(raw.tag) || normTag(p.tag) || "A";
    if (Array.isArray(raw.range)) return { kind: "range", value: null, range: raw.range.map(Number), tag, explicit: true, source: raw.source || EDIT_SOURCE };
    if ("value" in raw) return { kind: "value", value: raw.value, range: null, tag, explicit: true, source: raw.source || EDIT_SOURCE };
  }
  return { kind: "value", value: raw, range: null, tag: normTag(p.tag) || "A", explicit: false, source: p.source };
}

/** True when the leaf value is the registry default at the registry's tag (nothing to store). */
function isPlainDefault(p, leaf) {
  return leaf.kind === "value" && leaf.value === p.default && leaf.tag === (normTag(p.tag) || "A");
}

function visibleNow(p, params, specs) {
  if (!p.visible_if) return true;
  return Object.entries(p.visible_if).every(([name, allowed]) => {
    const dep = specs.find((s) => s.name === name);
    if (!dep) return true;
    const cur = readLeaf(dep, params[name]);
    const v = cur.kind === "range" ? null : cur.value;
    return (allowed || []).some((a) => String(a) === String(v));
  });
}

/**
 * renderParamForm({axis, entry, getParams, onEdit, onInspect})
 *   getParams(): the entry's current `params` object in the spec
 *   onEdit(name, leaf|null): store (or delete) the leaf; the caller revalidates
 *   onInspect({path, meta, value, inBand}, {quiet}): show the field in the inspector drawer
 * Returns {el, render}.
 */
export function renderParamForm({ axis, entry, getParams, onEdit, onInspect }) {
  const root = h("div.pf", { dataset: { axis, kind: entry.kind } });
  const specs = entry.params;

  function itemFor(p) {
    const params = getParams();
    const leaf = readLeaf(p, params[p.name]);
    const path = `${axis}.${p.name}`;
    const band = isNum(p.lo) && isNum(p.hi) ? [p.lo, p.hi] : null;
    const meta = { label: symbolText(p), tex: p.tex, unit: displayUnit(p.unit), tag: leaf.tag, source: leaf.source || p.source, band, choices: p.choices || undefined };
    const shown = leaf.kind === "range" ? (isNum(p.default) ? p.default : leaf.range[0]) : leaf.value;
    const item = h("div.pf-item", { class: leaf.explicit ? "is-edited" : "", dataset: { param: p.name, tag: leaf.tag } });

    const inspect = (_path, o) => onInspect?.({ path, meta, value: curValue(), inBand: inBandNow() }, o || {});
    const curValue = () => { const l = readLeaf(p, getParams()[p.name]); return l.kind === "range" ? `[${l.range[0]}, ${l.range[1]}]` : l.value; };
    const inBandNow = () => {
      const l = readLeaf(p, getParams()[p.name]);
      if (!band || l.kind === "range" || !isNum(l.value)) return null;
      return l.value >= band[0] && l.value <= band[1];
    };
    const rebuild = () => { render(); };
    const dependents = specs.some((s) => s.visible_if && p.name in s.visible_if);

    const store = (patch) => {
      const cur = readLeaf(p, getParams()[p.name]);
      const next = { ...cur, ...patch };
      if (next.kind === "value" && isPlainDefault(p, next)) { onEdit(p.name, null); return; }
      const out = { unit: p.unit, tag: next.tag, source: next.source || EDIT_SOURCE };
      if (next.kind === "range") { out.range = next.range; out.scale = p.scale; } else out.value = next.value;
      onEdit(p.name, out);
    };

    const row = fieldRow({
      path, meta, value: shown, ranged: leaf.kind === "range" ? leaf.range : null, envDefault: band,
      onChange: (v) => {
        if (v === null && typeof p.default !== "string") { onEdit(p.name, null); rebuild(); return; }
        const cur = readLeaf(p, getParams()[p.name]);
        // an edited number is the user's own: [A] unless they pick a tag themselves
        const tag = cur.explicit ? cur.tag : "A";
        store({ kind: "value", value: v, range: null, tag, explicit: true, source: cur.explicit ? cur.source : EDIT_SOURCE });
        syncEdited();
        if (dependents) rebuild();
      },
      onRange: (r) => {
        const cur = readLeaf(p, getParams()[p.name]);
        if (r) {
          const tag = cur.explicit && (cur.tag === "A" || cur.tag === "E") ? cur.tag : "A";
          store({ kind: "range", range: r, value: null, tag, explicit: true, source: cur.explicit ? cur.source : EDIT_SOURCE });
        } else {
          const v = isNum(cur.value) ? cur.value : (isNum(p.default) ? p.default : null);
          store({ kind: "value", value: v, range: null, tag: cur.explicit ? cur.tag : "A", explicit: true, source: cur.source });
        }
        syncEdited();
      },
      onInspect: inspect,
    });

    // per-leaf tag selector (line form, never hue)
    const sel = h("div.pf-tagsel", { role: "group", "aria-label": `${symbolText(p)} provenance tag` },
      TAGS.map((t) => {
        const b = tagChip(t, {
          size: "sm", title: `set the tag of this value to [${t}]`,
          onClick: () => {
            const cur = readLeaf(p, getParams()[p.name]);
            if (cur.tag === t && cur.explicit) return;
            const v = cur.value;
            store({ kind: cur.kind === "default" ? "value" : cur.kind, value: cur.kind === "range" ? null : v, tag: t, explicit: true, source: cur.explicit ? cur.source : EDIT_SOURCE });
            rebuild();
          },
        });
        b.setAttribute("aria-pressed", String(t === leaf.tag));
        b.classList.toggle("is-on", t === leaf.tag);
        return b;
      }));
    row.querySelector(".fr-tools")?.append(sel);
    // label: the registry id is a code name; show the symbol or plain words and keep the unit beside the value only
    const lab = row.querySelector(".fr-label");
    if (lab) { lab.replaceChildren(paramLabelNode(p)); lab.title = `${symbolText(p)} (registry id ${p.name})`; }
    item.dataset.find = `${p.name} ${symbolText(p)} ${richPlain(p.help || "")}`.toLowerCase();

    const reset = h("button.btn-link.pf-reset", {
      type: "button", hidden: !leaf.explicit, title: `back to the registry default (${String(p.default)}) and its tag [${normTag(p.tag)}]`,
      on: { click: () => { onEdit(p.name, null); rebuild(); } },
    }, "reset");

    function syncEdited() {
      const l = readLeaf(p, getParams()[p.name]);
      item.classList.toggle("is-edited", l.explicit);
      reset.hidden = !l.explicit;
      item.dataset.tag = l.tag;
      const chip = row.querySelector(".fr-tools > .tag");
      if (chip) {
        const fresh = tagChip(l.tag, { title: l.source, onClick: () => inspect(path) });
        chip.replaceWith(fresh);
      }
      for (const b of sel.children) { const on = b.dataset.tag === l.tag; b.setAttribute("aria-pressed", String(on)); b.classList.toggle("is-on", on); }
    }

    item.append(row,
      h("div.pf-foot",
        p.help ? h("p.pf-help", p.help) : null,
        h("p.pf-src", { title: leaf.source || "" }, leaf.explicit ? "your value" : `default: ${leaf.source || "no source recorded"}`),
        reset));
    return item;
  }

  // A long sheet (> 12 fields) gets a filter box; groups other than main/advanced get headings.
  let filter = "";
  const FILTER_MIN = 12;

  function applyFilter() {
    const q = filter.trim().toLowerCase();
    for (const it of root.querySelectorAll(".pf-item")) it.hidden = !!q && !it.dataset.find.includes(q);
    for (const g of root.querySelectorAll(".pf-group")) g.hidden = !!q && ![...g.querySelectorAll(".pf-item")].some((i) => !i.hidden);
    const adv = root.querySelector(".pf-adv");
    if (adv) {
      const any = [...adv.querySelectorAll(".pf-item")].some((i) => !i.hidden);
      adv.hidden = !!q && !any;
      if (q && any) adv.open = true;
    }
    const n = root.querySelectorAll(".pf-item:not([hidden])").length;
    const note = root.querySelector(".pf-filter-n");
    if (note) note.textContent = q ? `${n} match${n === 1 ? "" : "es"}` : "";
    const none = root.querySelector(".pf-none");
    if (none) none.hidden = !q || n > 0;
  }

  function render() {
    const params = getParams();
    clear(root);
    const vis = specs.filter((p) => visibleNow(p, params, specs));
    const main = vis.filter((p) => p.group !== "advanced");
    const adv = vis.filter((p) => p.group === "advanced");
    if (vis.length > FILTER_MIN) {
      const box = h("input.fr-num.pf-filter", {
        type: "search", value: filter, placeholder: `Filter ${vis.length} fields`, "aria-label": `Filter the ${vis.length} fields of this sheet`, spellcheck: false,
        on: { input: (e) => { filter = e.target.value; applyFilter(); } },
      });
      root.append(h("div.pf-filter-row", box, h("span.pf-filter-n", { "aria-live": "polite" })));
    }
    const groups = [];
    for (const p of main) {
      const g = p.group || "main";
      let b = groups.find((x) => x.name === g);
      if (!b) groups.push((b = { name: g, items: [] }));
      b.items.push(p);
    }
    if (groups.length > 1) {
      for (const g of groups) root.append(h("section.pf-group", h("h3.engrave.pf-group-h", g.name), h("div.pf-grid", g.items.map(itemFor))));
    } else root.append(h("div.pf-grid", main.map(itemFor)));
    if (adv.length) {
      const edited = adv.filter((p) => params[p.name] !== undefined).length;
      root.append(h("details.pf-adv", { open: edited > 0 },
        h("summary", h("span.engrave", "Advanced"), h("span.pf-adv-n", `${adv.length} fields${edited ? `, ${edited} edited` : ""}`)),
        h("div.pf-grid", adv.map(itemFor))));
    }
    if (!vis.length) root.append(h("p.sheet-note", "This kind has no parameters."));
    root.append(h("p.sheet-note.pf-none", { hidden: true }, "No field matches the filter."));
    if (filter) applyFilter();
  }

  render();
  return { el: root, render };
}
