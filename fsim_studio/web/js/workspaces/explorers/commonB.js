// Shared pieces of the model explorers phonon / transport / laser (studio-p2d, Coder E2).
// No physics here: sliders send inputs, GET/POST /api/explore/{phonon,transport,sde}/... returns every
// number, every tag and the pinned caveat text. Job-backed frames come through the job pool.
import { h, clear } from "../../ui/dom.js";
import { icon } from "../../ui/icons.js";
import { tagChip } from "../../ui/tagChip.js";
import { createReadout } from "../../ui/readout.js";
import { fmt } from "../../ui/format.js";

/** Slider with a tag chip and an optional log mapping. spec: {id,label,min,max,step?,default,unit?,tag,note?,scale?} */
export function ctl(prefix, spec, onInput) {
  const id = `${prefix}-${spec.id}`;
  const log = spec.scale === "log";
  const N = 1000;
  const toPos = (v) => (log ? (Math.log(v / spec.min) / Math.log(spec.max / spec.min)) * N : v);
  const fromPos = (p) => (log ? spec.min * Math.pow(spec.max / spec.min, p / N) : p);
  const roundTo = (v) => {
    if (log) return Number(v.toPrecision(3));
    const st = spec.step || (spec.max - spec.min) / 200;
    return Number((Math.round((v - spec.min) / st) * st + spec.min).toPrecision(6));
  };
  const out = h("output.xs-val", { for: id }, fmt(spec.default));
  const input = h("input.xs-range", {
    id, type: "range", min: log ? 0 : spec.min, max: log ? N : spec.max, step: log ? 1 : (spec.step || "any"), value: toPos(spec.default),
    "aria-label": `${spec.label}${spec.unit ? ` (${spec.unit})` : ""}`, dataset: { ctl: spec.id },
    on: { input: (e) => { const v = roundTo(fromPos(Number(e.target.value))); out.textContent = fmt(v); state.value = v; onInput?.(v); } },
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

/** Discrete slider over a list of values (frame selection). Calls onInput(value, index). */
export function stepper(prefix, { id, label, values, index = 0, unit = "", tag, note, onInput }) {
  const iid = `${prefix}-${id}`;
  const out = h("output.xs-val", { for: iid }, fmt(values[index]));
  const input = h("input.xs-range", {
    id: iid, type: "range", min: 0, max: values.length - 1, step: 1, value: index,
    "aria-label": `${label}${unit ? ` (${unit})` : ""}`, dataset: { ctl: id },
    on: { input: (e) => { const i = Number(e.target.value); state.index = i; out.textContent = fmt(values[i]); onInput?.(values[i], i); } },
  });
  const row = h("div.xs-row", { dataset: { ctl: id } },
    h("label.xs-label", { for: iid }, label, " ", tagChip(tag, { size: "sm", title: note || null })),
    h("div.xs-ctl", input, h("span.xs-read", out, unit ? h("span.xs-unit", unit) : null)));
  if (note) row.title = note;
  const state = { row, input, out, index, get value() { return values[state.index]; }, set(i) { state.index = i; input.value = i; out.textContent = fmt(values[i]); } };
  return state;
}

export function readoutGroup(defs) {
  const ro = {};
  const el = h("div.readouts.xp-readouts", { role: "group", "aria-label": "Readouts" },
    defs.map(([k, l]) => (ro[k] = createReadout({ key: k, label: l })).el));
  return { el, ro };
}

/** Pinned caveat: always visible, never collapsible. texts: strings from the backend (verbatim). */
export function caveatBox(texts, { id = "caveat" } = {}) {
  const body = h("div.xb-caveat-body");
  const box = h("aside.xb-caveat", { role: "note", "aria-label": "Model caveat", dataset: { caveat: id } }, icon("warn", { size: 16 }), body);
  const set = (list) => { clear(body).append(...list.filter(Boolean).map((t) => h("p", t))); };
  set(texts || []);
  return { el: box, set };
}

export function qs(obj) {
  const p = new URLSearchParams();
  for (const [k, v] of Object.entries(obj)) if (v != null) p.set(k, v);
  return p.toString();
}

/** A submit_call response ({result} when cached, else {job_id}) -> promise of the job result. */
export function jobResult(api, resp, hooks) {
  if (resp.result !== undefined && resp.result !== null) return Promise.resolve(resp.result);
  return api.awaitJob(resp.job_id, hooks).promise;
}

export function sourceLine(calls, extra = "") {
  return h("p.xp-src", icon("info", { size: 14 }), ` ${extra}${extra ? " " : ""}Calls: ${(calls || []).join(", ")}.`);
}

/** Run fn after ms of quiet, dropping superseded calls. Returns {run(), cancel()}. */
export function latest(fn, ms = 100) {
  let t = null;
  let seq = 0;
  return {
    run(...a) { clearTimeout(t); const my = ++seq; t = setTimeout(() => fn(my, () => my === seq, ...a), ms); },
    cancel() { clearTimeout(t); seq++; },
  };
}

export function chipText(tag) { return tag ? `[${String(tag).replace(/[[\]]/g, "")}]` : ""; }
