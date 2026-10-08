// Explain workspace (#/explain/<topic>): physics explainers ported from app.py and the 3D module.
//   spectral  X/XX Lorentzians through the filter window, leakage fill, eps readout
//             (GET /api/explain/spectral -> fsim_core.spectral / loading)
//   cascade   3D XX -> X -> 0 cascade (GET /api/scene/cascade, Lindblad populations) beside the
//             static formula g2 = 1 - rho^2 (1 - g2_0) (GET /api/explain/cascade -> integrator.g2_of_T)
//   band      band-profile scenes (GET /api/scene/band): tilted nitride band with psi_e / psi_h, InP flat levels
//   va-fit    the committed V-a calibration (GET /api/explain/va_fit -> out/phase0), residual strip
// No physics in this file: sliders send inputs, the backend returns every number.

import { h, clear, debounce } from "../ui/dom.js";
import { icon } from "../ui/icons.js";
import { tagChip } from "../ui/tagChip.js";
import { segmented, switchToggle, dataTable, skeleton, errorState } from "../ui/controls.js";
import { createReadout } from "../ui/readout.js";
import { fmt } from "../ui/format.js";
import { sro } from "../ui/figure.js";
import { chartCard } from "../charts/chartCard.js";
import { buildSpectral, spectralTable } from "../charts/explain.js";
import { buildCalibration } from "../charts/sweep.js";

const TOPICS = [
  { id: "spectral", label: "Spectral window", sub: "why a filter cannot remove the biexciton completely" },
  { id: "cascade", label: "Cascade", sub: "XX → X → 0, and the closed form for g²(0)" },
  { id: "band", label: "Band profile", sub: "where the electron and hole sit in the dot" },
  { id: "va-fit", label: "V-a calibration", sub: "the model against Chatzarakis 2023 g²(T)" },
  // model explorers (studio-p2d): each topic with a `module` is mounted from ./explorers/<module>.js
  { id: "cw-g2", label: "CW g²(τ)", sub: "dc HBT histogram, IRF and background (explorer)", module: "cwg2" },
  { id: "pulse-counting", label: "Pulse counting", sub: "pulsed g² from periodic photon counting (explorer)", module: "pulsed" },
  { id: "lindblad", label: "Quantum cross-check", sub: "Lindblad: Rabi drive, HOM, cascaded filter (explorer)", module: "quantum" },
  { id: "phonon", label: "Phonon sidebands", sub: "ZPL weight and filter transmission, independent-boson model (explorer)", module: "phonon" },
  { id: "transport", label: "Injection and Stark", sub: "InP p-i-n loading and background, nitride Stark tuning (explorer)", module: "transport" },
  { id: "laser", label: "QD laser (Zhao)", sub: "Langevin laser: g² against pump statistics (explorer)", module: "laser" },
];
const CASCADE_CARDS = ["edge-inp-gainp-design", "nitride-cavity-pulse-design", "nitride-cavity-set-design"];
const BAND_CARDS = [
  { card: "nitride-cavity-set-design", label: "InGaN/GaN (c-plane): tilted band" },
  { card: "edge-inp-gainp-design", label: "InP/GaInP: flat levels" },
];

let X = null;

export function mount(el, ctx) {
  X = { el, ctx, topic: null, charts: {}, scenes: [], disposers: [], seq: 0 };
  el.classList.add("ws-explain");
  X.disposers.push(ctx.state.on("theme", (m) => { for (const s of X?.scenes || []) s?.setTheme?.(m); }));
  update(ctx.args);
}

export function unmount() {
  if (!X) return;
  for (const d of X.disposers) d();
  dispose();
  X.el.classList.remove("ws-explain");
  X = null;
}

function dispose() {
  for (const c of Object.values(X.charts)) { try { c.dispose(); } catch { /* */ } }
  X.charts = {};
  for (const s of X.scenes) { try { s.dispose(); } catch { /* */ } }
  X.scenes = [];
}

export function update(args) {
  if (!X) return;
  let topic = args?.[0];
  if (!TOPICS.some((t) => t.id === topic)) { topic = "spectral"; history.replaceState(null, "", "#/explain/spectral"); }
  if (topic === X.topic) return;
  X.topic = topic;
  X.seq++;
  dispose();
  clear(X.el);
  const t = TOPICS.find((x) => x.id === topic);
  const body = h("div.xp-body", { dataset: { topic } });
  X.el.append(h("div.explain",
    h("header.ws-head.xp-head",
      h("div", h("h1.ws-title", t.label), h("p.ws-sub", t.sub)),
      h("nav.xp-tabs", { "aria-label": "Explainers" }, TOPICS.map((x) => h("a.res-tab", { href: `#/explain/${x.id}`, "aria-current": x.id === topic ? "page" : null }, x.label)))),
    body));
  if (t.module) { mountExplorer(t, body); return; }
  ({ spectral, cascade, band, "va-fit": vaFit })[topic](body);
}

// ---------------------------------------------------------------- explorers (studio-p2d)
// A topic with `module` loads ./explorers/<module>.js and calls mount(body, host); the host hands over
// the app ctx and a liveness test so late replies from a previous topic are dropped.
async function mountExplorer(t, body) {
  const seq = X.seq;
  const host = { ctx: X.ctx, charts: X.charts, topic: t.id, alive: () => !!X && X.seq === seq };
  try {
    const mod = await import(`./explorers/${t.module}.js`);
    if (!host.alive()) return;
    mod.mount(body, host);
  } catch (e) {
    if (!host.alive()) return;
    clear(body).append(errorState({ title: `The ${t.label} explorer failed to load`, message: e.message }));
  }
}

// ---------------------------------------------------------------- helpers
function slider({ label, min, max, step, value, unit = "", onInput, id }) {
  const out = h("output.xs-val", { for: id }, fmt(value));
  const input = h("input.xs-range", { id, type: "range", min, max, step, value, "aria-label": `${label}${unit ? ` (${unit})` : ""}`,
    on: { input: (e) => { out.textContent = fmt(Number(e.target.value)); onInput(Number(e.target.value)); } } });
  const row = h("div.xs-row", h("label.xs-label", { for: id }, label), h("div.xs-ctl", input, h("span.xs-read", out, unit ? h("span.xs-unit", unit) : null)));
  return { row, input, out, set(v) { input.value = v; out.textContent = fmt(v); }, disable(d) { input.disabled = d; row.classList.toggle("is-disabled", d); } };
}

function readoutGroup(defs) {
  const ro = {};
  const el = h("div.readouts.xp-readouts", { role: "group", "aria-label": "Readouts" }, defs.map(([k, l]) => (ro[k] = createReadout({ key: k, label: l })).el));
  return { el, ro };
}

async function mountSceneInto(slot, kind, params) {
  clear(slot).append(h("div.vp-placeholder", h("div.vp-plate", h("span.engrave", "3D view loading"), h("span.vp-sub", params.card || ""))));
  const seq = X.seq;
  const spec = await X.ctx.api.sceneSpec(kind, params);
  if (!X || seq !== X.seq) return null;
  const mod = await import("../viz3d/index.js");
  if (!X || seq !== X.seq) return null;
  clear(slot);
  const s = mod.mountScene(slot, kind, spec, { theme: X.ctx.state.get("theme") });
  X.scenes.push(s);
  return { scene: s, spec };
}

// ================================================================ spectral
function spectral(body) {
  const st = { T: 78, pub: true, w: 2, dx: 0, mu: null, mode: "PL", dg_inj: 2, I_ratio: 1, last: null };
  const chart = chartCard({
    title: "Spectral window: X and XX through the filter", height: 470,
    build: () => buildSpectral(st.last), table: () => spectralTable(st.last), filename: () => "fsim-spectral-window",
  });
  X.charts.spec = chart;
  const status = h("p.xp-status", { "aria-live": "polite" });
  const sT = slider({ id: "xs-T", label: "Temperature T", min: 4, max: 320, step: 1, value: st.T, unit: "K", onInput: (v) => { st.T = v; go(); } });
  const sW = slider({ id: "xs-w", label: "Filter width w", min: 0.1, max: 15, step: 0.1, value: st.w, unit: "meV", onInput: (v) => { st.w = v; go(); } });
  const sDx = slider({ id: "xs-dx", label: "X offset from window centre dx", min: -6, max: 6, step: 0.1, value: st.dx, unit: "meV", onInput: (v) => { st.dx = v; go(); } });
  const sMu = slider({ id: "xs-mu", label: "Mean loading μ (F1b)", min: 0, max: 5, step: 0.05, value: 0.33, unit: "", onInput: (v) => { st.mu = v; go(); } });
  const pubSw = switchToggle({ label: "published window at this T (card g2_vs_T dataset)", checked: true, onChange: (v) => { st.pub = v; syncPub(); go(); } });
  const elBox = h("div.xp-el", { hidden: true });
  const dg = h("input.fr-num", { type: "number", min: 0, max: 20, step: 0.5, value: st.dg_inj, "aria-label": "dGamma_inj (meV)", on: { change: (e) => { st.dg_inj = Number(e.target.value); go(); } } });
  const ir = h("input.fr-num", { type: "number", min: 0, max: 20, step: 0.5, value: st.I_ratio, "aria-label": "I / I_ref", on: { change: (e) => { st.I_ratio = Number(e.target.value); go(); } } });
  elBox.append(h("label.xs-inline", h("span", "ΔΓ_inj"), dg, h("span.xs-unit", "meV")), h("label.xs-inline", h("span", "I / I_ref"), ir));
  const mode = segmented({ label: "Drive mode", value: "PL", options: [{ value: "PL", label: "PL" }, { value: "EL", label: "EL" }], onChange: (v) => { st.mode = v; elBox.hidden = v !== "EL"; go(); } });
  const { el: roEl, ro } = readoutGroup([["gamma", "Linewidth Γ"], ["ratio", "Γ_EL / Γ_PL"], ["eps", "Leakage ε = t_XX / t_X"], ["g2dot", "g²₀(μ)"], ["tx", "t_X (brightness)"], ["window", "Window w, dx"]]);
  const src = h("p.xp-src");
  function syncPub() { sW.disable(st.pub); sDx.disable(st.pub); }
  syncPub();
  body.append(h("div.xp-split",
    h("div.xp-main", chart.el, h("p.xp-cap", "The window passes most of X but also a tail of XX; that tail, relative to X, is ε. Wider lines (higher T, electrical injection) push more XX through any window.")),
    h("aside.xp-side.panel",
      h("header.panel-head", h("h2.engrave", "Inputs"), status),
      h("div.xp-ctls", sT.row, h("div.xs-row", pubSw), sW.row, sDx.row, sMu.row, h("div.xs-row", h("span.xs-label", "Drive"), mode.el), elBox),
      roEl, src)));

  let inflight = 0;
  const go = debounce(async () => {
    const my = ++inflight;
    const qs = new URLSearchParams({ T: st.T, pub: st.pub ? 1 : 0, w: st.w, dx: st.dx, mode: st.mode, dg_inj: st.dg_inj, I_ratio: st.I_ratio });
    if (st.mu != null) qs.set("mu", st.mu);
    status.textContent = "evaluating";
    let r;
    try { r = await X.ctx.api.getJSON(`/api/explain/spectral?${qs}`); } catch (e) {
      if (my !== inflight || !X) return;
      status.textContent = "";
      chart.setEmpty(errorState({ title: "The spectral explainer did not answer", message: e.message }));
      st.last = null; chart.render();
      return;
    }
    if (my !== inflight || !X || X.topic !== "spectral") return;
    st.last = r;
    if (st.mu == null) { st.mu = r.mu; sMu.set(r.mu); }
    if (st.pub) { sW.set(r.w); sDx.set(r.dx); }
    status.textContent = "";
    const t = r.tags || {};
    ro.gamma.update({ value: r.gamma, unit: "meV", tag: t.linewidth, qualifier: r.mode, caption: `PL ${fmt(r.gamma_pl)} · EL ${fmt(r.gamma_el)} meV`, markChange: false });
    ro.ratio.update({ value: r.gamma_ratio, unit: "×", tag: t.linewidth, caption: r.mode === "EL" ? "injection broadening, F5′ (loading.gamma_eff)" : "PL: no injection broadening", markChange: false });
    ro.eps.update({ value: r.eps, tag: t.chain, caption: "F1, spectral.epsilon", markChange: false });
    ro.g2dot.update({ value: r.g2_dot, tag: t.chain, caption: r.g2_dot_formula, markChange: false });
    ro.tx.update({ value: r.t_x, tag: t.chain, caption: "X transmission through the window", markChange: false });
    ro.window.update({ value: r.w, unit: "meV", tag: t.window, caption: `dx ${fmt(r.dx)} meV · ${r.published_window ? "published window, interpolated from the card" : "manual window [A]"}`, markChange: false });
    chart.setTitle("Spectral window: X and XX through the filter", `T = ${fmt(r.T)} K · ${r.mode} · Δ_XX ${fmt(r.delta_xx)} meV · μ ${fmt(r.mu)}`);
    chart.setTag(t.chain);
    chart.render();
    clear(src).append(icon("info", { size: 14 }), ` ${r.params_source}. Calls: ${r.calls.join(", ")}.`);
  }, 70);
  go();
}

// ================================================================ cascade
function cascade(body) {
  const st = { card: CASCADE_CARDS[0], T: 150 };
  const slot = h("div.xp-viewport", { role: "region", "aria-label": "3D cascade" });
  const sceneNote = h("p.xp-cap");
  const tabs = h("div.hero-tabs", { role: "tablist", "aria-label": "Cascade card" }, CASCADE_CARDS.map((c) => h("button.hero-tab", {
    type: "button", role: "tab", dataset: { card: c }, "aria-selected": String(c === st.card),
    on: { click: () => { st.card = c; for (const b of tabs.children) b.setAttribute("aria-selected", String(b.dataset.card === c)); loadScene(); } },
  }, c.replace(/-design$/, ""))));
  const sT = slider({ id: "xc-T", label: "Temperature for the rates", min: 4, max: 320, step: 1, value: st.T, unit: "K", onInput: (v) => { st.T = v; go(); } });
  const formula = h("div.formula", { role: "math", "aria-live": "polite" });
  const { el: roEl, ro } = readoutGroup([["g2", "g²(0)"], ["rho", "ρ (signal fraction)"], ["eps", "ε (XX leakage)"], ["g2dot", "g²₀ = F1b(μ, ε)"], ["p1", "P₁ (one exciton)"], ["p2", "P₂ (two or more)"]]);
  const src = h("p.xp-src");
  body.append(h("div.xp-split",
    h("div.xp-main",
      h("section.viewport-mount.xp-vm",
        h("div.vm-rail", h("span.engrave", "Cascade in time · design card"), tabs),
        slot),
      sceneNote),
    h("aside.xp-side.panel",
      h("header.panel-head", h("h2.engrave", "Closed form · V-a card"), h("span.panel-note", "integrator.g2_of_T")),
      h("div.xp-ctls", sT.row),
      formula, roEl, src)));

  async function loadScene() {
    for (const s of X.scenes) { try { s.dispose(); } catch { /* */ } }
    X.scenes = [];
    try {
      const r = await mountSceneInto(slot, "cascade", { card: st.card });
      if (r) sceneNote.textContent = `${r.spec.title}. ${(r.spec.labels || []).slice(0, 2).join(". ")}. The 3D view is the design card at its operating point; the closed form on the right is the calibrated V-a card. They are different devices.`;
    } catch (e) {
      clear(slot).append(h("div.vp-placeholder", h("div.vp-plate", h("span.engrave", "3D view unavailable"), h("span.vp-sub", e.message))));
    }
  }
  const go = debounce(async () => {
    let r;
    try { r = await X.ctx.api.getJSON(`/api/explain/cascade?T=${st.T}`); } catch (e) {
      clear(formula).append(errorState({ title: "The cascade explainer did not answer", message: e.message }));
      return;
    }
    if (!X || X.topic !== "cascade") return;
    const t = r.tags || {};
    clear(formula).append(
      h("p.f-line", h("span.f-sym", "g²(0)"), " = 1 − ", h("span.f-sym", "ρ²"), " (1 − ", h("span.f-sym", "g²₀"), ")"),
      h("p.f-line.f-num", "= 1 − ", h("span.f-v", fmt(r.rho2)), " × (1 − ", h("span.f-v", fmt(r.g2_dot)), ") = ", h("span.f-res", fmt(r.g2)), tagChip(t.chain, { title: "tag chain of the V-a fit" })),
      h("p.f-note", `At μ → 0 the cascade gives g²₀ = ε, so g²(0) = 1 − ρ²(1 − ε). Here μ = ${fmt(r.mu)} (F1b).`));
    ro.g2.update({ value: r.g2, tag: t.chain, caption: `at T = ${fmt(r.T)} K`, markChange: false });
    ro.rho.update({ value: r.rho, tag: t.chain, caption: `ρ² = ${fmt(r.rho2)}`, markChange: false });
    ro.eps.update({ value: r.eps, tag: t.chain, caption: `window w ${fmt(r.w)}, dx ${fmt(r.dx)} meV (published)`, markChange: false });
    ro.g2dot.update({ value: r.g2_dot, tag: t.chain, caption: `μ = ${fmt(r.mu)}`, markChange: false });
    ro.p1.update({ value: r.P1, tag: t.chain, caption: "cap-2 Poisson loading", markChange: false });
    ro.p2.update({ value: r.P2, tag: t.chain, caption: "loading.loading_probs", markChange: false });
    clear(src).append(icon("info", { size: 14 }), ` ${r.params_source}. Calls: ${r.calls.join(", ")}.`);
  }, 70);
  loadScene();
  go();
}

// ================================================================ band
function band(body) {
  const grid = h("div.xp-bands");
  body.append(grid, h("p.xp-cap", "Band edges and levels come from fsim_core.nitride_levels and dot_levels. The nitride panel shows the built-in polarization field tilting the band, pulling the electron and hole apart (|ψ|² overlap is the readout). dot_levels returns no wavefunctions for InP, so none are drawn."));
  for (const b of BAND_CARDS) {
    const slot = h("div.xp-viewport.is-band", { role: "region", "aria-label": `Band profile, ${b.label}` });
    const notes = h("ul.lat-notes");
    grid.append(h("section.viewport-mount.xp-vm", h("div.vm-rail", h("span.engrave", b.label), h("span.vm-note", b.card)), slot, notes));
    mountSceneInto(slot, "band", { card: b.card }).then((r) => {
      if (r) clear(notes).append(...(r.spec.labels || []).map((t) => h("li", t)));
    }).catch((e) => clear(slot).append(h("div.vp-placeholder", h("div.vp-plate", h("span.engrave", "3D view unavailable"), h("span.vp-sub", e.message)))));
  }
}

// ================================================================ V-a fit
async function vaFit(body) {
  body.append(skeleton("chart", { lines: 0, height: 440 }));
  let va;
  try { va = await X.ctx.api.getJSON("/api/explain/va_fit"); } catch (e) {
    clear(body).append(errorState({ title: "Could not read out/phase0", message: e.message }));
    return;
  }
  if (!X || X.topic !== "va-fit") return;
  const tag = String(va.tag_chain || "").replace(/[[\]]/g, "");
  const chart = chartCard({
    title: "g²(T): Chatzarakis 2023 data and the fitted model", height: 470,
    build: () => buildCalibration(va),
    table: () => ({ columns: [{ key: "T_K", label: "T (K)" }, { key: "g2", label: "measured" }, { key: "err", label: "± err" }, { key: "bound", label: "bound" }, { key: "model_g2", label: "model" }, { key: "residual", label: "residual" }], rows: va.points }),
    filename: () => "fsim-va-calibration",
  });
  X.charts.va = chart;
  const { el: roEl, ro } = readoutGroup([["tc", "T_c (fit)"], ["cost", "Fit cost"]]);
  const passed = sro({ label: "passed ±0.03 (every point)", value: va.passed, tag, note: "fit_params.json passed_pm0.03", caption: "fit_params.json passed_pm0.03", size: "lg" });
  clear(body).append(h("div.xp-split",
    h("div.xp-main",
      h("p.calib-label", icon("warn", { size: 16 }), h("strong", va.label), " · digitized data · no held-out (P) prediction exists yet · the 78 K point is an upper bound (▽)"),
      chart.el),
    h("aside.xp-side.panel",
      h("header.panel-head", h("h2.engrave", "Fit bundle"), h("span.panel-note", "out/phase0, as committed")),
      roEl,
      h("div.pad", passed),
      h("div.pad", dataTable({ caption: "residual = model − data", columns: [{ key: "T", label: "T (K)", align: "right" }, { key: "r", label: "residual", align: "right" }],
        rows: va.points.map((p) => ({ T: fmt(p.T_K), r: h("span", `${p.residual >= 0 ? "+" : ""}${fmt(p.residual)} `, tagChip(tag, { size: "sm" })) })) })),
      h("p.xp-src", icon("info", { size: 14 }), ` Files: ${va.files.join(", ")}. The fit is stored, never re-run here.`))));
  chart.setTitle("g²(T): Chatzarakis 2023 data and the fitted model", "measured points [V] · model [A] · residual strip below");
  chart.setTag(tag);
  chart.render();
  ro.tc.update({ value: va.Tc_K, unit: "K", tag, caption: "where the fitted g²(T) crosses 0.5", markChange: false });
  ro.cost.update({ value: va.cost, tag, caption: "least-squares cost as stored", markChange: false });
}

