// "Pulse-end occupation" panel of the Composer's Results area: the exciton occupation at the END of the pulse window,
// started from the ground state, with the phonon and radiative-decay switches. No physics here: the number is the
// return value of POST /api/compose/pulse-end (a job, because the phonon-on call takes about two seconds).
import { h } from "./dom.js";
import { icon } from "./icons.js";
import { tex } from "./math.js";
import { switchToggle } from "./controls.js";

const RHO_TEX = String.raw`\rho_{XX}`;
const key = (spec, ph, ra) => JSON.stringify([spec, ph, ra]);

/** api: ctx.api (submitJob); getSpec(): the current composed spec. */
export function pulseEndPanel({ api, getSpec }) {
  const st = { ph: true, ra: true, state: "idle", result: null, error: "", doneKey: null, preset: "", job: null, visible: false };
  const valueEl = h("div.pe-value", { dataset: { role: "pe-value" } });
  const stateEl = h("div.pe-state", { role: "status", "aria-live": "polite", dataset: { role: "pe-state" } });
  const presetEl = h("p.pe-preset", { dataset: { role: "pe-preset" }, hidden: true });
  const btn = h("button.btn.btn-outline.pe-compute", { type: "button", dataset: { role: "pe-compute" }, on: { click: () => compute() } }, icon("run", { size: 14 }), "Compute");
  const sw = (label, role, k) => {
    const t = switchToggle({ label, checked: st[k], onChange: (v) => { st[k] = v; sync(); } });
    t.dataset.role = role;
    return t;
  };
  const swPh = sw("Phonons", "pe-phonons", "ph");
  const swRa = sw("Radiative decay", "pe-radiative", "ra");
  const el = h("section.pe-panel", { dataset: { role: "pulse-end" }, "aria-label": "Pulse-end occupation", hidden: true },
    h("header.pe-head", h("h2.pe-title", "Pulse-end occupation"), stateEl),
    h("div.pe-controls", swPh, swRa, btn),
    presetEl,
    h("div.pe-out", valueEl),
    h("p.pe-cap", "Starts from the ground state and is the occupation at the end of the pulse window, not the steady-state readout above."));

  function currentKey() { return key(getSpec(), st.ph, st.ra); }
  function stale() { return st.state === "done" && st.doneKey !== currentKey(); }

  function sync() {
    el.hidden = !st.visible;
    if (!st.visible) return;
    const loading = st.state === "loading";
    const isStale = stale();
    el.classList.toggle("is-stale", isStale);
    el.classList.toggle("is-loading", loading);
    el.dataset.state = loading ? "loading" : st.state === "error" ? "error" : isStale ? "stale" : st.state;
    btn.disabled = loading;
    swPh.querySelector("input").checked = st.ph;
    swRa.querySelector("input").checked = st.ra;
    presetEl.hidden = !st.preset;
    presetEl.textContent = st.preset;
    valueEl.replaceChildren();
    stateEl.replaceChildren();
    if (loading) {
      stateEl.append(icon("retry", { size: 13 }), "Computing");
    } else if (st.state === "error") {
      valueEl.append(h("p.pe-error", { dataset: { role: "pe-error" } }, icon("warn", { size: 14 }), st.error));
    } else if (st.result) {
      valueEl.append(tex(RHO_TEX), h("span.pe-eq", " at pulse end = "), h("span.pe-num", { dataset: { role: "pe-number" } }, st.result.value.toFixed(3)));
      if (isStale) stateEl.append(h("span.stale-badge", { dataset: { role: "pe-stale" }, title: "The design or a switch changed after this number was computed" }, icon("warn", { size: 13 }), "Stale"));
    } else {
      valueEl.append(h("span.pe-none", "Not computed yet"));
    }
  }

  function compute() {
    if (st.state === "loading") return;
    st.job?.cancel?.();
    const spec = JSON.parse(JSON.stringify(getSpec()));
    const ph = st.ph, ra = st.ra, k = key(spec, ph, ra);
    st.state = "loading"; st.error = "";
    sync();
    const job = api.submitJob("/api/compose/pulse-end", { spec, phonons: ph, radiative: ra });
    st.job = job;
    job.promise.then((r) => {
      if (st.job !== job) return;
      st.job = null; st.result = r.result; st.doneKey = k; st.state = "done"; sync();
    }).catch((e) => {
      if (st.job !== job) return;
      st.job = null;
      const b = e.body;
      st.error = (typeof b?.error === "string" ? b.error : b?.error?.message) || e.message;
      st.state = "error"; st.result = null; sync();
    });
  }

  return {
    el,
    /** visible: the drive is pulsed; called after every spec change and render. */
    sync(visible) { if (visible !== undefined) st.visible = !!visible; sync(); },
    /** Preset the switches from a paper optimum item's solver ({phonons, radiative}: true = on) and say so. */
    preset(solver, text) {
      st.ph = !!solver.phonons; st.ra = !!solver.radiative; st.preset = text || "";
      st.state = "idle"; st.result = null; st.doneKey = null;
      sync();
    },
    dispose() { st.job?.cancel?.(); },
  };
}
