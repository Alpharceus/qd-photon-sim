// Compare workspace: slots A/B/C pinned from the Designer shelf. Overlay hero chart in the
// comparison set (series 1-3), envelopes at lower opacity, delta table against a chosen
// reference (interval notation when banded; direction-of-good as glyph + word), design diff,
// and the report bundle export (POST /api/compare/report -> fsim_viz.report.designer_report).

import { h, clear } from "../ui/dom.js";
import { icon } from "../ui/icons.js";
import { tagChip, normTag, widest } from "../ui/tagChip.js";
import { segmented, dataTable, emptyState, errorState, toast } from "../ui/controls.js";
import { fmt, fmtInterval, isNum } from "../ui/format.js";
import { chartCard } from "../charts/chartCard.js";
import { buildG2T, buildSecondary, SECONDARY } from "../charts/g2T.js";
import { colors as themeColors } from "../charts/theme.js";

const SLOTS = ["A", "B", "C"];
// Metric rows: direction of good is a display convention for the delta arrows (no physics).
const METRICS = [
  { key: "g2_op", label: "g²(0) at operating point", unit: "", good: "lower" },
  { key: "T_c", label: "Ceiling temperature T_c", unit: "K", good: "higher" },
  { key: "brightness_per_pulse", label: "Brightness per pulse", unit: "", good: "higher" },
  { key: "collected_flux_pulsed_s", label: "Collected flux, commanded", unit: "s⁻¹", good: "higher" },
  { key: "collected_flux_delivered_s", label: "Collected flux, delivered", unit: "s⁻¹", good: "higher" },
  { key: "T_j_op", label: "Junction temperature T_j", unit: "K", good: "lower" },
  { key: "eps_op", label: "XX leakage ε", unit: "", good: "lower" },
  { key: "rho_op", label: "Loading purity ρ", unit: "", good: "higher" },
];

let C = null;

export function mount(el, ctx) {
  C = { el, ctx, charts: {}, secondary: "eps", report: null, disposers: [] };
  el.classList.add("ws-compare");
  C.disposers.push(ctx.state.on("slots", render));
  C.disposers.push(ctx.state.on("compareRef", render));
  render();
}

export function unmount() {
  if (!C) return;
  for (const d of C.disposers) d();
  for (const c of Object.values(C.charts)) c.dispose();
  C.el.classList.remove("ws-compare");
  C = null;
}

function slots() { return C.ctx.state.get("slots") || {}; }
function filled() { return SLOTS.filter((k) => slots()[k]); }

function slotColor(k) { return themeColors().series[SLOTS.indexOf(k)]; }

function render() {
  if (!C) return;
  for (const c of Object.values(C.charts)) c.dispose();
  C.charts = {};
  clear(C.el);
  const s = slots();
  const have = filled();
  let ref = C.ctx.state.get("compareRef");
  if (!s[ref]) ref = have[0] || "A";

  const head = h("header.cmp-head",
    h("div", h("h1.ws-title", "Compare"), h("p.ws-sub", "Pinned runs side by side. Envelopes stay bands; deltas stay intervals.")),
    h("div.cmp-actions", exportButton(have)));

  const plates = h("div.slot-plates", SLOTS.map((k) => slotPlate(k, s[k], k === ref)));

  if (!have.length) {
    C.el.append(h("div.compare", head, plates, h("div.cmp-empty", emptyState({
      title: "No runs pinned yet",
      body: "Run a design, then use Pin to A on the instrument shelf. Pin up to three runs to overlay them here.",
      action: { label: "Open the Designer", fn: () => C.ctx.navigate(`#/design/${C.ctx.state.get("card") || "edge-inp-gainp-design"}`) },
      iconName: "pin",
    }))));
    return;
  }

  const refSeg = segmented({
    label: "Reference slot",
    value: ref,
    options: SLOTS.map((k) => ({ value: k, label: `Ref ${k}`, disabled: !s[k] })),
    onChange: (v) => C.ctx.state.set({ compareRef: v }),
  });

  const runs = have.map((k) => ({ result: s[k].result, opT: s[k].opT, color: slotColor(k), name: `${k}: ${short(s[k].card)}` }));
  C.charts.hero = chartCard({
    title: "g²(0) against heatsink temperature, overlay",
    subtitle: have.map((k) => `${k} ${s[k].run_id}`).join("   "),
    build: () => buildG2T(have.map((k) => ({ result: slots()[k].result, opT: slots()[k].opT, color: slotColor(k), name: `${k}: ${short(slots()[k].card)}` })), { compare: true }),
    table: () => overlayTable(have),
    filename: () => `compare-${have.join("")}-g2T`,
  });
  C.charts.hero.setTag(widestTag(have));
  C.charts.hero.el.classList.add("chart-hero");
  C.charts.hero.setEmpty(emptyState({ title: "No curves to overlay", body: "Pinned runs are single-temperature points; run the full curve to overlay g²(T).", iconName: "chart" }));

  const tabs = h("div.sec-tabs", { role: "tablist", "aria-label": "Secondary quantity" },
    Object.entries(SECONDARY).map(([k]) => h("button.sec-tab", {
      type: "button", role: "tab", "aria-selected": String(C.secondary === k),
      on: { click: () => { C.secondary = k; render(); } },
    }, { eps: "ε", rho2: "ρ²", dTJ: "ΔT_J", gamma: "Γ" }[k])));
  C.charts.sec = chartCard({
    title: SECONDARY[C.secondary].label,
    height: 220,
    build: () => buildSecondary(C.secondary, have.map((k) => ({ result: slots()[k].result, color: slotColor(k), name: `${k}: ${short(slots()[k].card)}` })), { compare: true }),
    table: () => null,
    filename: () => `compare-${C.secondary}`,
  });
  C.charts.sec.setEmpty(h("p.chart-sub.pad", "Secondary curves need full T_hs grids."));
  void runs;

  C.el.append(h("div.compare", head, plates,
    h("div.cmp-grid",
      h("div.cmp-main", C.charts.hero.el, h("div.sec-wrap", tabs, C.charts.sec.el)),
      h("div.cmp-side",
        h("section.panel", { "aria-label": "Delta table" },
          h("header.panel-head", h("h2.engrave", "Deltas against the reference"), refSeg.el),
          deltaTable(have, ref)),
        h("section.panel", { "aria-label": "Design differences" },
          h("header.panel-head", h("h2.engrave", "Design diff"), h("span.panel-note", "fields that differ between slots")),
          designDiff(have)),
        C.report ? reportPlate(C.report) : null))));
  C.charts.hero.render();
  C.charts.sec.render();
}

function short(card) { return String(card || "").replace(/-design$/, ""); }

function widestTag(have) {
  const order = ["V", "DR", "E", "A"];
  let best = -1;
  for (const k of have) best = Math.max(best, order.indexOf(normTag(slots()[k].result.tag_chain)));
  return best >= 0 ? order[best] : null;
}

function slotPlate(k, entry, isRef) {
  if (!entry) {
    return h("article.slot-plate.is-empty", { "aria-label": `Slot ${k}, empty` },
      h("div.slot-top", h("span.slot-badge", { dataset: { slot: k } }, k), h("span.engrave", "Empty")),
      h("p.slot-hint", "Pin a run from the Designer shelf."));
  }
  const r = entry.result;
  return h("article.slot-plate", { class: isRef ? "is-ref" : "", "aria-label": `Slot ${k}: ${entry.card}` },
    h("div.slot-top",
      h("span.slot-badge", { dataset: { slot: k } }, k),
      h("span.slot-card", short(entry.card)),
      isRef ? h("span.chip", "reference") : null,
      h("button.icon-btn.slot-x", {
        type: "button", "aria-label": `Unpin slot ${k}`,
        on: { click: () => { const all = { ...slots() }; all[k] = null; C.ctx.state.set({ slots: all }); } },
      }, icon("close", { size: 16 }))),
    h("div.slot-meta",
      h("span.id-chip.is-static", h("span.id-v", r.run_id)),
      tagChip(r.tag_chain),
      h("span.slot-mode", r.mode === "envelope" ? `envelope, ${r.n_samples} samples` : r.mode === "headline" ? "headline model" : "point"),
      (r.labels || []).includes("static (non-headline)") ? h("span.label-chip", "static (non-headline)") : null),
    h("a.slot-open", { href: `#/design/${entry.card}` }, "Open in Designer"));
}

// An envelope run's scalars without a band are its mid-design point (all ranges at midpoint):
// shown with a "mid-design" qualifier and never compared better/worse.
function scalarSpec(r, key) {
  const band = r.scalar_bands?.[key];
  const tag = cellTag(r, key);
  if (band) return { lo: band[0], hi: band[1], band: true, tag };
  const v = r.scalars?.[key];
  return { value: v, reason: r.scalars?.[`${key}__nan_reason`], mid: r.mode === "envelope", tag };
}

/** Tag of a cell: the run's tag chain, widened by the quantity's own provenance entry when present. */
function cellTag(r, key) {
  const t = normTag(r.tag_chain);
  const pk = key === "brightness_per_pulse" ? "brightness" : null;
  const own = pk ? normTag(r.provenance?.[pk]?.tag) : null;
  return own ? widest([t, own]) : t;
}

function cellText(sp, unit) {
  if (sp.band) return `${fmtInterval(sp.lo, sp.hi)}${unit ? ` ${unit}` : ""}`;
  if (isNum(sp.value)) return `${fmt(sp.value)}${unit ? ` ${unit}` : ""}`;
  return null;
}

function deltaCell(sp, refSp, m) {
  if (sp.mid || refSp.mid) {
    return h("span.delta", { dataset: { verdict: "midpoint" } }, icon("minus", { size: 14 }), h("span.delta-w", "midpoint, not compared"));
  }
  if (!sp.tag || !refSp.tag) return h("span.nan", "no delta: a value has no provenance tag");
  const a = sp.band ? [sp.lo, sp.hi] : [sp.value, sp.value];
  const b = refSp.band ? [refSp.lo, refSp.hi] : [refSp.value, refSp.value];
  if (![...a, ...b].every(isNum)) return h("span.nan", "no delta: a value is missing");
  // interval difference [a.lo - b.hi, a.hi - b.lo]; a point difference when neither is banded
  const lo = a[0] - b[1], hi = a[1] - b[0];
  const banded = sp.band || refSp.band;
  let dir = "same";
  if (lo > 0) dir = "up"; else if (hi < 0) dir = "down"; else if (banded) dir = "overlap"; else dir = "same";
  const better = dir === "up" ? m.good === "higher" : dir === "down" ? m.good === "lower" : null;
  const word = dir === "overlap" ? "intervals overlap" : dir === "same" ? "no change" : better ? "better" : "worse";
  const glyph = dir === "up" ? "arrowUp" : dir === "down" ? "arrowDown" : "minus";
  const txt = banded ? `Δ ${fmtInterval(lo, hi)}` : `Δ ${lo > 0 ? "+" : ""}${fmt(lo)}`;
  return h("span.delta", { dataset: { verdict: better === null ? "neutral" : better ? "better" : "worse" } },
    icon(glyph, { size: 14 }), h("span.delta-v", txt), h("span.delta-w", word));
}

function deltaTable(have, ref) {
  const s = slots();
  const refR = s[ref].result;
  const rows = METRICS.filter((m) => have.some((k) => {
    const sp = scalarSpec(s[k].result, m.key);
    return sp.band || isNum(sp.value);
  }));
  const table = h("table.data-table.delta-table",
    h("caption.sr-only", `Deltas of each slot against slot ${ref}`),
    h("thead", h("tr", h("th", { scope: "col" }, "Metric"), have.map((k) => h("th", { scope: "col" }, h("span.slot-badge.sm", { dataset: { slot: k } }, k))))),
    h("tbody", rows.map((m) => h("tr",
      h("th", { scope: "row" }, m.label, h("span.good", m.good === "lower" ? "lower is better" : "higher is better")),
      have.map((k) => {
        const r = s[k].result;
        const sp = scalarSpec(r, m.key);
        let t = cellText(sp, m.unit);
        const chip = tagChip(sp.tag, { size: "sm" });
        // a number without a tag does not render (CONTRACT rule 2)
        const withheld = t && !chip;
        if (withheld) t = null;
        const refSp = scalarSpec(refR, m.key);
        const conv = m.key === "brightness_per_pulse" && t ? r.scalars?.brightness_convention : null;
        return h("td.num", { dataset: { metric: m.key, slot: k, tagged: chip ? "1" : "0", mid: sp.mid && t ? "1" : "0" } },
          t ? h("span.cell-line", h("span.cell-v", t), chip) : h("span.nan", withheld ? "no provenance tag; number withheld" : sp.reason || "n/a"),
          t && sp.mid ? h("span.cell-qual", { title: "envelope run: this scalar has no band, so it is the design with every range at its midpoint (not a prediction)" }, "mid-design") : null,
          conv ? h("span.cell-conv", conv) : null,
          k === ref ? h("span.delta.is-ref", "reference") : t ? deltaCell(sp, refSp, m) : null);
      })))));
  return h("div.table-wrap", table);
}

function flatten(obj, pre = "", out = {}) {
  for (const [k, v] of Object.entries(obj || {})) {
    const p = pre ? `${pre}.${k}` : k;
    if (v && typeof v === "object" && !Array.isArray(v)) flatten(v, p, out);
    else out[p] = Array.isArray(v) ? JSON.stringify(v) : v;
  }
  return out;
}

function designDiff(have) {
  const s = slots();
  if (have.length < 2) return h("p.panel-note.pad", "Pin a second run to see which fields differ.");
  const flats = have.map((k) => flatten(s[k].design));
  const keys = new Set(flats.flatMap((f) => Object.keys(f)));
  const diff = [...keys].filter((p) => p !== "name" && new Set(flats.map((f) => JSON.stringify(f[p] ?? null))).size > 1).sort();
  if (!diff.length) return h("p.panel-note.pad", "The pinned designs are identical; only the run mode differs.");
  const shown = diff.slice(0, 60);
  return h("div",
    dataTable({
      className: "diff-table",
      columns: [{ key: "path", label: "Field" }, ...have.map((k) => ({ key: k, label: k, align: "right" }))],
      rows: shown.map((p) => Object.fromEntries([["path", p], ...have.map((k, i) => [k, flats[i][p] === undefined ? "absent" : typeof flats[i][p] === "number" ? fmt(flats[i][p]) : String(flats[i][p])])])),
    }),
    diff.length > shown.length ? h("p.panel-note.pad", `${diff.length - shown.length} more differing fields not listed.`) : null);
}

function overlayTable(have) {
  const s = slots();
  const rows = [];
  for (const k of have) {
    const r = s[k].result;
    const x = r.curves?.T_hs || [];
    x.forEach((t, i) => {
      const b = r.bands?.g2;
      rows.push({ slot: k, T: t, g2: b ? `${fmt(b.lo[i])} – ${fmt(b.hi[i])}` : r.curves.g2[i] });
    });
  }
  return { columns: [{ key: "slot", label: "Slot", align: "left" }, { key: "T", label: "T_hs (K)" }, { key: "g2", label: "g²(0)" }], rows };
}

function exportButton(have) {
  return h("button.btn.btn-outline", {
    type: "button", disabled: !have.length,
    title: have.length ? "Write the designer report bundle for the pinned runs" : "Pin a run first",
    on: {
      click: async (e) => {
        const btn = e.currentTarget;
        btn.disabled = true;
        btn.setAttribute("aria-busy", "true");
        try {
          const s = slots();
          const r = await C.ctx.api.compareReport(have.map((k) => ({ label: `${k} ${short(s[k].card)}`, run_id: s[k].run_id })), "FSIM Studio comparison");
          C.report = { ok: true, ...r };
          toast("Report bundle written.", { kind: "ok" });
        } catch (err) {
          C.report = { ok: false, message: err.message };
        }
        render();
      },
    },
  }, icon("folder", { size: 16 }), "Export report bundle");
}

function reportPlate(rep) {
  if (!rep.ok) return errorState({ title: "Report bundle failed", message: rep.message, actions: [] });
  return h("section.panel.report-plate", { "aria-label": "Report bundle" },
    h("header.panel-head", h("h2.engrave", "Report bundle written"), icon("check", { size: 16 })),
    h("p.report-path", h("code", rep.path)),
    rep.files?.length ? h("ul.report-files", rep.files.map((f) => h("li", f))) : null);
}
