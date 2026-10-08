// Shared pieces of the model explorers cwg2 / pulsed / quantum (studio-p2d, Coder E1).
// No physics here: sliders send inputs, GET /api/explore/... returns every number and the
// pinned caveat text. Control ranges, defaults and tags come from GET /api/explore/a/controls.
import { h, clear, debounce } from "../../ui/dom.js";
import { icon } from "../../ui/icons.js";
import { tagChip } from "../../ui/tagChip.js";
import { createReadout } from "../../ui/readout.js";
import { fmt } from "../../ui/format.js";

let controlsPromise = null;
export function loadControls(api) {
  if (!controlsPromise) controlsPromise = api.getJSON("/api/explore/a/controls").catch((e) => { controlsPromise = null; throw e; });
  return controlsPromise;
}

/** Slider with a tag chip and (optionally) a log mapping. spec = one entry of /controls. */
export function ctlSlider(prefix, spec, onInput) {
  const id = `${prefix}-${spec.id}`;
  const log = spec.scale === "log";
  const N = 1000;
  const toPos = (v) => (log ? (Math.log(v / spec.min) / Math.log(spec.max / spec.min)) * N : v);
  const fromPos = (p) => {
    if (log) return spec.min * Math.pow(spec.max / spec.min, p / N);
    return p;
  };
  const roundTo = (v) => {
    if (log) return Number(v.toPrecision(3));
    const st = spec.step || (spec.max - spec.min) / 200;
    return Number((Math.round((v - spec.min) / st) * st + spec.min).toPrecision(6));
  };
  const out = h("output.xs-val", { for: id }, fmt(spec.default));
  const input = h("input.xs-range", {
    id, type: "range", min: log ? 0 : spec.min, max: log ? N : spec.max, step: log ? 1 : (spec.step || "any"), value: toPos(spec.default),
    "aria-label": `${spec.label}${spec.unit ? ` (${spec.unit})` : ""}`, dataset: { ctl: spec.id },
    on: { input: (e) => { const v = roundTo(fromPos(Number(e.target.value))); out.textContent = fmt(v); state.value = v; onInput(v); } },
  });
  const chip = tagChip(spec.tag, { size: "sm", title: spec.note || null });
  const row = h("div.xs-row", { dataset: { ctl: spec.id } },
    h("label.xs-label", { for: id }, spec.label, " ", chip),
    h("div.xs-ctl", input, h("span.xs-read", out, spec.unit ? h("span.xs-unit", spec.unit) : null)));
  if (spec.note) row.title = spec.note;
  const state = { value: spec.default, row, input, out, spec };
  state.set = (v) => { state.value = v; input.value = toPos(v); out.textContent = fmt(v); };
  state.disable = (d) => { input.disabled = d; row.classList.toggle("is-disabled", d); };
  return state;
}

export function readoutGroup(defs) {
  const ro = {};
  const el = h("div.readouts.xp-readouts", { role: "group", "aria-label": "Readouts" },
    defs.map(([k, l]) => (ro[k] = createReadout({ key: k, label: l })).el));
  return { el, ro };
}

/** Pinned caveat: always visible, never collapsible. texts: array of strings from the backend. */
export function caveatBox(texts, { id = "caveat" } = {}) {
  const box = h("aside.xa-caveat", { role: "note", "aria-label": "Model caveat", dataset: { caveat: id } },
    icon("warn", { size: 16 }), h("div.xa-caveat-body"));
  const body = box.querySelector(".xa-caveat-body");
  const set = (list) => { clear(body).append(...list.filter(Boolean).map((t) => h("p", t))); };
  set(texts || []);
  return { el: box, set };
}

/** Debounced fetch with stale-reply protection. fn(json) renders; err(e) shows the failure. */
export function liveFetch(host, build, render, { ms = 100, status } = {}) {
  let inflight = 0;
  const run = debounce(async () => {
    const my = ++inflight;
    let url;
    try { url = build(); } catch (e) { return; }
    if (url == null) return;
    if (status) status.textContent = "computing";
    try {
      const r = await host.ctx.api.getJSON(url);
      if (my !== inflight || !host.alive()) return;
      if (status) status.textContent = "";
      render(r);
    } catch (e) {
      if (my !== inflight || !host.alive()) return;
      if (status) status.textContent = `error: ${e.message}`;
      render(null, e);
    }
  }, ms);
  return run;
}

export function qs(obj) {
  const p = new URLSearchParams();
  for (const [k, v] of Object.entries(obj)) if (v != null) p.set(k, v);
  return p.toString();
}

/** One footer line listing the fsim_core calls and the parameter source. */
export function sourceLine(r) {
  return h("p.xp-src", icon("info", { size: 14 }), ` ${r.params_source ? `${r.params_source}. ` : ""}Calls: ${(r.calls || []).join(", ")}.`);
}
