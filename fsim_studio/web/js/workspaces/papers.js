// Papers workspace (#/papers, #/papers/<id>): the paper-device library.
//   left   catalogue of paper plaques (cite + catalog ID, kind chips, read level by line form, class counts,
//          last ledger outcome n/N), filterable by material / structure / drive kind
//   main   device header (Reproduce key = the only beam-filled control, Open in Composer, job progress),
//          model-vs-paper table (one row per ledger item), series overlays (paperOverlay.js)
// No physics here: every number comes from GET /api/papers* and the reproduce job result. The only
// arithmetic is display: delta = model - paper, and the widest tag of the inputs for a model value.

import { h, clear } from "../ui/dom.js";
import { icon } from "../ui/icons.js";
import { tagChip, normTag, widest } from "../ui/tagChip.js";
import { skeleton, errorState, toast } from "../ui/controls.js";
import { fmt, isNum } from "../ui/format.js";
import { richText, tex, symbolFor, texIndex, pathSymbol } from "../ui/math.js";
import { paperOverlay } from "../charts/paperOverlay.js";

const AXES = [
  { key: "material", label: "Material" },
  { key: "structure", label: "Structure" },
  { key: "drive", label: "Drive" },
];

const CLASS_INFO = {
  C: { label: "Calibration", tip: "Calibration: the model was fitted to this paper number (listed in the card's fits), so agreement is by construction. It is not a test." },
  M: { label: "Unfitted", tip: "Unfitted comparison: the model was built from other inputs and then compared with this paper number. The only class that counts as a test." },
  T: { label: "Transcription", tip: "Transcription: the observable passes an input taken from this paper straight through, so it can only agree. It is not a test." },
  P: { label: "Prediction", tip: "Prediction: the paper reports no value for this observable. The model number is stated, not tested." },
  R: { label: "Paper reference", tip: "Paper reference: the reported number is the paper's own theory or simulation (or was tuned by the paper), so the comparison is model against another model. It is shown, not counted as a test." },
};

const READ_INFO = {
  "V-full": "Read in full: every number was checked against the paper text",
  "V-abs": "Read from the abstract only",
  "V-snip": "Read from a snippet only",
};

// Display notation of the observables the ledgers use (rich strings: `$...$` is typeset; the registry's own names are sentences).
const OBS_RICH = {
  g2_0: "$g^{(2)}(0)$", g2_XX: "$g^{(2)}(0)$ of the XX line", g2_raw: "$g^{(2)}(0)$", g2_binned: "binned $g^{(2)}(0)$", g2_irf: "IRF-convolved $g^{(2)}(0)$",
  hom_I: "HOM overlap $I_{\\mathrm{HOM}}$", hom_V: "HOM visibility $V_{\\mathrm{HOM}}$", hom_XX: "HOM overlap of the XX line", hom_I_gated: "gated HOM overlap $I_{\\mathrm{HOM}}^{\\mathrm{gate}}$",
  F_P_eff: "Effective Purcell factor $F_{\\mathrm P}L$", F_P_peak: "Purcell factor $F_{\\mathrm P}$ at resonance", tau_ps: "Lifetime $\\tau$", tau_XX_ps: "XX lifetime $\\tau_{XX}$",
  tau_purcell_ps: "Purcell lifetime $\\tau_{\\mathrm P}$", tau_EL_ps: "EL decay time $\\tau_{\\mathrm{EL}}$", Gamma_ns: "Radiative rate $\\Gamma$", beta: "$\\beta$ factor",
  eta_extract: "Extraction $\\beta\\,\\eta_{\\mathrm{out}}$", eta_out: "Extraction $\\eta_{\\mathrm{out}}$", eta_prep: "State preparation $\\eta_{\\mathrm{prep}}$",
  eta_rad: "Radiative efficiency $\\eta_{\\mathrm{rad}}$", eta_rad_eff: "Effective radiative efficiency $\\eta_{\\mathrm{rad}}^{\\mathrm{eff}}$",
  B_fibre: "Brightness into the fibre $B_{\\mathrm{fibre}}$", B_first_lens: "First-lens brightness $B_{\\mathrm{lens}}$", B_lens: "First-lens brightness $B_{\\mathrm{lens}}$",
  rate_detected_Hz: "Detected count rate $R_{\\mathrm{det}}$", Z_phonon: "ZPL weight $Z(T)$", inversion: "Inversion $w$ after the pulse",
  linewidth_ueV: "Emitter linewidth $\\hbar\\Gamma_{\\mathrm{tot}}$", lambda_nm: "Wavelength $\\lambda$", spectral_overlap: "Spectral leakage $\\varepsilon_{XX}$",
};

/** "raw g(2)(0)": the convention in words, then the observable in notation (intrinsic is the default and stays silent). */
function obsNode(obs, conv) {
  const c = conv && conv !== "intrinsic" ? `${String(conv).replace(/_/g, " ")} ` : "";
  const body = OBS_RICH[obs] || (obs ? `$${symbolFor({ name: obs })}$` : "");
  return h("span.pp-obsl", { title: "paper and model are compared in the same convention only" }, c, richText(body));
}

// A ledger item id ("eta_X_QD1", "g2_ps_vs_T") is a file key: shown as a typeset head symbol plus words, never as snake_case.
const ITEM_HEAD = {
  g2: String.raw`$g^{(2)}(0)$`, eta: String.raw`$\eta$`, tau: String.raw`$\tau$`, beta: String.raw`$\beta$`, lambda: String.raw`$\lambda$`,
  gamma: String.raw`$\Gamma$`, F: String.raw`$F_{\mathrm P}$`, FP: String.raw`$F_{\mathrm P}$`, Z: String.raw`$Z$`, T1: String.raw`$T_1$`,
  Fprep: String.raw`$F_{\mathrm{prep}}$`, M: String.raw`$M$`, hom: "HOM", linewidth: "linewidth", rate: "rate", inversion: "inversion",
  extraction: "extraction", brightness: "brightness", sys: "system efficiency", pQD: String.raw`$p_{\mathrm{QD}}$`, CX: String.raw`$C_X$`,
};
const ITEM_WORD = {
  pred: "prediction", prediction: "prediction", corr: "corrected", raw: "raw", meas: "measured", sim: "simulated", sat: "at saturation",
  res: "resonant", env: "envelope", lying: "mislabelled", dc: "dark-count corrected", darkcorr: "dark-count corrected", analyt: "analytic",
  bgc: "background-corrected", irf: "IRF", ps: "pulsed", vs: "vs", frep: String.raw`$f_{\mathrm{rep}}$`, T: String.raw`$T$`, FP: String.raw`$F_{\mathrm P}$`,
  rabi: "Rabi", xcorr: "cross-correlation", EL: "EL", eff: "efficiency", dev: "device", rad: "radiative", cav: "cavity", fit: "fit",
  max: "max", star: String.raw`$\gamma^{*}$`, det: "detected", ebg: "eBG", cbg: "CBG", TPE: "TPE", TPI: "TPI", XX: "XX", X: "X",
};
function itemToken(t) {
  let m;
  if (ITEM_WORD[t] !== undefined) return ITEM_WORD[t];
  if ((m = /^(\d+)p(\d+)([A-Za-z]+)$/.exec(t)) || (m = /^(\d+)()([A-Za-z]+)$/.exec(t))) {
    const num = m[2] ? `${m[1]}.${m[2]}` : m[1];
    const u = m[3].replace(/^u/, "µ");
    return u === "Psat" ? String.raw`${num}$\,P_{\mathrm{sat}}$` : `${num} ${u}`;
  }
  return t;
}
function itemLabel(id) {
  const toks = String(id).split("_").filter(Boolean);
  const head = ITEM_HEAD[toks[0]];
  const rest = toks.slice(head ? 1 : 0).map(itemToken);
  return [head, ...rest].filter(Boolean).join(" · ");
}

function leafNode(path) {
  const s = String(path);
  const axis = s.split(".")[0];
  const word = { material: "material", structure: "structure", drive: "drive", setup: "setup" }[axis] || axis;
  return h("span.pp-leaf", `${word} `, tex(pathSymbol(s, P?.texIdx || null)));
}

/** The ledger's tolerance string in words: "2sigma|abs0.005" -> "±2 sigma or ±0.005 abs (whichever is larger)". */
function tolHuman(str) {
  if (!str) return "";
  const parts = String(str).split("|").map((t) => {
    let m;
    if ((m = /^([\d.]+)sigma$/.exec(t))) return String.raw`±${m[1]}$\,\sigma$`;
    if ((m = /^rel([\d.]+)%$/.exec(t))) return `±${m[1]}%`;
    if ((m = /^abs(.+)$/.exec(t))) return `±${m[1]} abs`;
    if ((m = /^(.+)<=(.+)$/.exec(t))) return `${m[1].replace(/_/g, " ")} at most ${m[2]}`;
    return t;
  });
  return parts.length > 1 ? `${parts.join(" or ")} (whichever is larger)` : parts[0];
}

/** Items that are tests (class M) and how many of them the model agrees with. Calibration, transcription, prediction and paper-reference (R)
 *  items are not tests. Agreement counts the CENTRAL value within the tolerance only: "pass-envelope" (paper value inside the input
 *  envelope but the central value misses) is counted apart as inRange, never as agreement. Weak passes are flagged in the ledger. */
function agreement(items) {
  const tested = items.filter((i) => i.cls === "M");
  const agreed = tested.filter((i) => i.outcome === "pass");
  return { tested: tested.length, agree: agreed.length, inRange: tested.filter((i) => i.outcome === "pass-envelope").length,
    weak: agreed.filter((i) => i.weak).length };
}

const NARROW = window.matchMedia("(max-width: 1100px)");
let P = null;

export function mount(el, ctx) {
  P = {
    el, ctx, list: [], id: null, paper: null, result: null, resultCached: false, running: false,
    charts: [], seq: 0, filters: { material: null, structure: null, drive: null }, dom: {},
    job: null, timer: null, t0: 0, eta: 5, listErr: null, agree: {},
  };
  el.classList.add("ws-papers");
  ensureCss();
  el.dataset.ready = "0";
  buildShell();
  NARROW.addEventListener("change", syncFold);
  loadList().then(() => { if (P) update(P.ctx.args); });
  ctx.api.composeRegistry?.().then((reg) => { if (P) { P.texIdx = texIndex(reg); if (P.paper) renderTable(); } }).catch(() => { /* symbols fall back to the derived ones */ });
}

export function unmount() {
  if (!P) return;
  NARROW.removeEventListener("change", syncFold);
  stopJob();
  disposeCharts();
  P.el.classList.remove("ws-papers");
  delete P.el.dataset.ready;
  P = null;
}

export function update(args) {
  if (!P) return;
  const id = args?.[0] || null;
  if (id === P.id && P.paperLoaded) return;
  P.id = id;
  P.paperLoaded = false;
  P.seq++;
  stopJob();
  disposeCharts();
  P.paper = null;
  P.result = null;
  P.resultCached = false;
  syncSelected();
  P.dom.grid.dataset.sel = id ? "1" : "0";
  syncFold();
  if (id && NARROW.matches) window.scrollTo?.(0, 0);
  if (!id) { renderPick(); P.el.dataset.ready = P.listDone ? "1" : "0"; return; }
  openPaper(id);
}

function ensureCss() {
  const href = new URL("../../css/papers.css", import.meta.url).href;
  const have = [...document.querySelectorAll('link[rel="stylesheet"]')].some((l) => l.href === href);
  if (!have) document.head.append(h("link", { rel: "stylesheet", href }));
}

function disposeCharts() {
  for (const c of P.charts) { try { c.dispose(); } catch { /* gone */ } }
  P.charts = [];
}

function stopJob() {
  clearInterval(P.timer);
  P.timer = null;
  try { P.job?.close(); } catch { /* closed */ }
  P.job = null;
  P.running = false;
}

// ---------------------------------------------------------------- shell + catalogue
function buildShell() {
  const d = P.dom;
  d.filters = h("div.pp-filters", { role: "group", "aria-label": "Filter the catalogue by material, structure or drive" });
  d.count = h("p.pp-count", { "aria-live": "polite" });
  d.key = h("p.pp-key", { hidden: true, title: "What an agreement is worth" },
    h("b", "M"), " unfitted test", h("span.pp-key-sep", "·"), h("b", "C"), " calibration", h("span.pp-key-sep", "·"), h("b", "T"), " transcription", h("span.pp-key-sep", "·"), h("b", "P"), " prediction", h("span.pp-key-sep", "·"), h("b", "R"), " paper reference");
  d.plaques = h("div.pp-plaques", skeleton("plaque", { lines: 4 }), skeleton("plaque", { lines: 4 }));
  // On a phone the catalogue folds behind "Papers (n)" and the chosen paper comes first; on desktop it is always open.
  d.catN = h("span.pp-cat-n");
  d.fold = h("details.pp-cat-fold", { open: true },
    h("summary.pp-cat-sum", h("span.engrave", "Papers"), d.catN),
    h("div.pp-cat-body", d.filters, d.count, d.key, d.plaques));
  d.fold.addEventListener("toggle", () => { if (!NARROW.matches && !d.fold.open) d.fold.open = true; });
  d.cat = h("aside.pp-cat", { "aria-label": "Paper catalogue" }, d.fold);
  d.main = h("section.pp-main", { "aria-label": "Paper device", dataset: { state: "idle" } });
  d.grid = h("div.pp-grid", { dataset: { sel: "0" } }, d.cat, d.main);
  P.el.append(h("div.papers",
    h("header.ws-head.pp-head",
      h("div", h("h1.ws-title", "Papers"),
        h("p.ws-sub", "Each card builds a device from only what one paper reports, runs the composed engine, and records per reported number whether the model agrees. A mismatch is a recorded result, never tuned away.")),
      h("div.pp-head-links",
        h("a.pp-legacy", { href: "#/library", title: "The older literature parameter cards (V-a tiers)" }, "Literature cards (legacy tiers)"))),
    d.grid));
}

/** Narrow screens: the catalogue is a disclosure, closed once a paper is chosen. Wide screens: always open. */
function syncFold() {
  const f = P?.dom?.fold;
  if (!f) return;
  f.open = NARROW.matches ? !P.id : true;
}

async function loadList() {
  try {
    P.list = await P.ctx.api.getJSON("/api/papers");
    P.listErr = null;
  } catch (e) {
    if (!P) return;
    P.list = [];
    P.listErr = e;
  }
  if (!P) return;
  P.listDone = true;
  renderCatalogue();
  prefetchAgreement();
}

/** A cached result comes back inline from the reproduce endpoint (nothing is computed), so the list can lead with agreement. */
async function prefetchAgreement() {
  const seq = P?.seq;
  for (const p of (P?.list || []).filter((x) => x.last_outcome?.cached && !P.agree[x.id]).slice(0, 40)) {
    if (!P || P.listErr) return;
    try {
      const resp = await P.ctx.api.postJSON(`/api/papers/${encodeURIComponent(p.id)}/reproduce`, {});
      if (!P) return;
      if (resp.result?.items) P.agree[p.id] = agreement(resp.result.items);
    } catch { /* the list keeps its ledger-only badge */ }
  }
  if (P && seq === P.seq) renderCatalogue();
}

function kindsOf(axis) {
  return [...new Set(P.list.map((p) => p.kinds?.[axis]).filter(Boolean))].sort();
}

function visible() {
  return P.list.filter((p) => AXES.every((a) => !P.filters[a.key] || p.kinds?.[a.key] === P.filters[a.key]));
}

function renderCatalogue() {
  const d = P.dom;
  clear(d.filters);
  if (P.list.length > 1) {
    for (const a of AXES) {
      const kinds = kindsOf(a.key);
      if (kinds.length < 2) continue;
      d.filters.append(h("div.pp-fgroup", { role: "group", "aria-label": a.label },
        h("span.pp-flabel.engrave", a.label),
        h("div.pp-fchips", kinds.map((k) => h("button.chip.pp-filter", {
          type: "button", "aria-pressed": String(P.filters[a.key] === k), dataset: { axis: a.key, kind: k },
          title: `${a.label}: ${k}. Click again to clear.`,
          on: { click: () => { P.filters[a.key] = P.filters[a.key] === k ? null : k; renderCatalogue(); } },
        }, k.replace(/_/g, " "))))));
    }
    if (AXES.some((a) => P.filters[a.key])) {
      d.filters.append(h("button.btn.btn-quiet.pp-clear", { type: "button", on: { click: () => { for (const a of AXES) P.filters[a.key] = null; renderCatalogue(); } } }, icon("close", { size: 14 }), "Clear filters"));
    }
  }
  clear(d.plaques);
  d.catN.textContent = P.list.length ? ` (${P.list.length})` : "";
  d.key.hidden = !P.list.length;
  if (P.listErr) {
    d.count.textContent = "";
    d.plaques.append(errorState({ title: "Could not read the paper catalogue", message: P.listErr.message, actions: [{ label: "Retry", fn: () => { loadList(); } }] }));
    return;
  }
  if (!P.list.length) {
    d.count.textContent = "0 papers";
    d.plaques.append(h("div.cplaque.pp-empty", { dataset: { role: "empty" } },
      h("div.cp-top", h("div", h("h2.cp-title", "No paper cards yet"))),
      h("p.pp-empty-body", "A paper card is one YAML file, ", h("code", "cards/papers/<id>.yaml"), ", that builds a device from what a paper reports and lists, per number, what the model is expected to do."),
      h("p.pp-empty-body", "Copy the schema and the five rules for adding a paper from ", h("code", "cards/papers/README.md"), ", then reload this page. Files starting with an underscore are fixtures and are not listed.")));
    return;
  }
  const vis = visible();
  d.count.textContent = vis.length === P.list.length ? `${P.list.length} papers` : `${vis.length} of ${P.list.length} papers`;
  if (!vis.length) {
    d.plaques.append(h("div.cplaque.is-missing", h("p.pp-empty-body", "No paper matches this combination of filters.")));
    return;
  }
  d.plaques.append(...vis.map(plaque));
  syncSelected();
}

function shortCite(cite, id) {
  if (!cite) return id;
  const year = [...String(cite).matchAll(/\((\d{4})\)/g)].map((m) => m[1]).pop() || "";
  const first = String(cite).split(",")[0].replace(/\s+et al\.?/i, "").trim();
  const sur = first.split(/\s+/).pop() || id;
  const many = /et al/i.test(cite) || /^[^,]+,\s*[A-Z]\.\s/.test(cite);
  return `${sur}${many ? " et al." : ""}${year ? ` ${year}` : ""}`;
}

/** Short cite, with the card id appended when two cards share it (e.g. two lines of one paper). */
function titleOf(id, cite) {
  const t = shortCite(cite, id);
  const dup = P.list.filter((p) => shortCite(p.cite, p.id) === t).length > 1;
  return dup ? `${t} · ${id.replace(/^[a-z]+\d{4}_?/i, "") || id}` : t;
}

function readBadge(level) {
  return h("span.pp-read", { dataset: { level: level || "unknown" }, title: READ_INFO[level] || "read level not stated" },
    h("span.pp-read-k", "read"), h("span.pp-read-v", level || "?"));
}

/** Last outcome on a list plaque: what was tested first ("5/7 AGREE"), the ledger as a quiet second part. */
function lastBadge(lo, id) {
  if (!lo) return h("span.vb.vb-sm.vb-info.pp-last", { dataset: { status: "none" }, title: "No reproduce result is cached for this card and engine yet" }, icon("minus", { size: 12 }), h("span.vb-word", "NOT RUN"));
  if (lo.n_items == null) return h("span.vb.vb-sm.vb-info.pp-last", { dataset: { status: "info" } }, icon("info", { size: 12 }), h("span.vb-word", "CACHED"));
  const ag = P.agree[id] || (lo.n_tested != null ? { tested: lo.n_tested, agree: lo.n_agree } : null);
  const ledger = h("span.pp-ledger", { title: `${lo.n_match} of ${lo.n_items} ledger items gave the outcome the card expected. This says the run reproduces the recorded results, including recorded failures; it is not model-equals-paper.` }, `ledger ${lo.n_match}/${lo.n_items}`);
  if (!ag || !ag.tested) {
    // agreement not known from the list alone: never show a green word for the ledger alone
    return h("span.pp-last-wrap", h("span.vb.vb-sm.vb-info.pp-last", { dataset: { status: "info" }, title: ag ? "This card has no unfitted (M) test items" : "Open the paper to see how many tested items the model agrees with" },
      icon("info", { size: 12 }), h("span.vb-word", ag ? "NO TESTS" : "OPEN TO SEE")), ledger);
  }
  const ok = ag.agree === ag.tested;
  return h("span.pp-last-wrap", h("span.vb.vb-sm.pp-last", { class: ok ? "vb-pass" : "vb-fail", dataset: { status: ok ? "pass" : "fail" },
    title: `The model agrees with the paper on ${ag.agree} of ${ag.tested} tested (M) items.` },
  icon(ok ? "check" : "cross", { size: 12 }), h("span.vb-word", "AGREE"), h("span.pp-frac", `${ag.agree}/${ag.tested}`)), ledger);
}

function classCounts(c) {
  return h("div.pp-classes", { role: "group", "aria-label": "Items per evidence class" },
    ["M", "C", "T", "P", "R"].filter((k) => k !== "R" || c?.R != null).map((k) => h("span.pp-cc", { class: c?.[k] ? "" : "is-zero", dataset: { cls: k }, title: `${c?.[k] || 0} ${CLASS_INFO[k].label.toLowerCase()} item(s). ${CLASS_INFO[k].tip}` },
      h("span.pp-cc-k", k), h("span.pp-cc-n", String(c?.[k] || 0)))));
}

function plaque(p) {
  return h("a.cplaque.pp-plaque", { href: `#/papers/${encodeURIComponent(p.id)}`, dataset: { paper: p.id }, "aria-label": `${titleOf(p.id, p.cite)}: open` },
    h("div.cp-top",
      h("div.cp-titles", h("h2.cp-title", titleOf(p.id, p.cite)), h("p.pp-cite", { title: p.cite || "" }, p.cite || "")),
      h("span.cp-go", icon("chevronRight", { size: 18 }))),
    h("div.pp-kinds", AXES.map((a) => p.kinds?.[a.key] ? h("span.chip.pp-kind", { title: `${a.label}: ${p.kinds[a.key]}` }, p.kinds[a.key].replace(/_/g, " ")) : null)),
    h("div.pp-row", readBadge(p.read_level), classCounts(p.classes)),
    h("div.pp-foot",
      lastBadge(p.last_outcome, p.id),
      h("span.id-chip.is-static", { title: "Catalog ID" }, h("span.id-k", "ID"), h("span.id-v", `FS-paper-${p.id}`))));
}

function syncSelected() {
  if (!P?.dom.plaques) return;
  for (const a of P.dom.plaques.querySelectorAll(".pp-plaque")) {
    if (a.dataset.paper === P.id) a.setAttribute("aria-current", "true"); else a.removeAttribute("aria-current");
  }
}

// ---------------------------------------------------------------- main: pick / load
function renderPick() {
  const m = P.dom.main;
  m.dataset.state = "idle";
  clear(m);
  if (P.listErr || !P.list.length) { m.append(h("p.pp-hint", P.listErr ? "" : "Add a card to start.")); return; }
  m.append(h("div.cplaque.pp-pick",
    h("div.cp-top", h("div", h("h2.cp-title", "Choose a paper"))),
    h("p.pp-empty-body", "Open a plaque on the left. Reproduce runs the composed engine on the device built from that paper and fills the model column; nothing is run until you press it, except that an already-cached result for the same card bytes and engine is shown at once."),
    h("p.pp-empty-body", "Class badges say what an agreement is worth: ",
      ...["M", "C", "T", "P", "R"].map((k) => [classBadge(k), " "]), "Only M items are tests, and only when the central model value is within the tolerance; a paper value that is merely inside the input envelope is shown as IN RANGE ONLY and is not counted.")));
}

async function openPaper(id) {
  const seq = P.seq;
  const m = P.dom.main;
  m.dataset.state = "loading";
  P.el.dataset.ready = "0";
  clear(m).append(skeleton("block", { lines: 5, height: 200 }), skeleton("block", { lines: 6, height: 260 }));
  let paper;
  try {
    paper = await P.ctx.api.getJSON(`/api/papers/${encodeURIComponent(id)}`);
  } catch (e) {
    if (!P || seq !== P.seq) return;
    m.dataset.state = "error";
    clear(m).append(errorState({
      title: e.status === 404 ? `No paper card ${id}` : `Paper card ${id} could not be loaded`,
      message: e.message + (e.body?.problems?.length ? `\n${e.body.problems.join("\n")}` : ""),
      actions: [{ label: "Back to the catalogue", fn: () => P.ctx.navigate("#/papers") }],
    }));
    P.el.dataset.ready = "1";
    return;
  }
  if (!P || seq !== P.seq) return;
  P.paper = paper;
  P.paperLoaded = true;
  renderPaper();
  // An exact cached result (same card bytes + engine hash) is shown at once; this POST returns it inline.
  if (paper.last_outcome?.cached) await reproduce({ auto: true });
  if (!P || seq !== P.seq) return;
  P.el.dataset.ready = "1";
}

// ---------------------------------------------------------------- main: paper page
function classBadge(k) {
  const i = CLASS_INFO[k];
  if (!i) return h("span.label-chip.pp-cls", { dataset: { cls: "none" } }, String(k || "?"));
  return h("span.label-chip.pp-cls", { dataset: { cls: k }, title: i.tip, "aria-label": `${k}: ${i.tip}` },
    h("b.pp-cls-k", k), h("span.pp-cls-l", i.label));
}

function outcomeBadge(outcome) {
  const o = String(outcome || "error").toLowerCase();
  const spec = o === "pass" ? ["vb-pass", "check", "PASS"]
    : o === "pass-envelope" ? ["vb-info pp-envelope", "info", "IN RANGE ONLY"]
      : o === "fail" ? ["vb-fail", "cross", "FAIL"]
        : o === "unavailable" ? ["vb-info", "info", "UNAVAILABLE"]
          : ["vb-fail", "warn", "ERROR"];
  const b = h("span.vb.vb-sm.pp-outcome", { class: spec[0], dataset: { outcome: o } }, icon(spec[1], { size: 12 }), h("span.vb-word", spec[2]));
  if (o === "pass-envelope") b.title = "The central model value misses the tolerance; the paper value lies inside the model's input envelope only. Not counted as agreement.";
  return b;
}

/** Outcome differs from the card's expectation: a split badge (outcome | LEDGER MISMATCH). */
function mismatchBadge(outcome, expect) {
  const o = String(outcome || "error").toLowerCase();
  const word = o === "pass" ? "PASS" : o === "pass-envelope" ? "IN RANGE ONLY" : o === "fail" ? "FAIL" : o === "unavailable" ? "UNAVAILABLE" : "ERROR";
  const ico = o === "pass" ? "check" : o === "unavailable" || o === "pass-envelope" ? "info" : o === "fail" ? "cross" : "warn";
  const half = o === "pass" ? "is-pass" : o === "unavailable" || o === "pass-envelope" ? "is-info" : "is-fail";
  return h("span.vb.vb-sm.vb-split.pp-outcome.pp-mismatch", {
    dataset: { outcome: o, mismatch: "1" }, role: "img",
    "aria-label": `${word}; ledger mismatch: the card expected ${String(expect).toUpperCase()}`,
    title: `The card expected ${String(expect).toUpperCase()}; the run gave ${word}. A flip in either direction needs a human review before the card changes.`,
  }, h("span.vb-half", { class: half }, icon(ico, { size: 12 }), word),
  h("span.vb-half.pp-mm", icon("warn", { size: 12 }), `LEDGER MISMATCH · expected ${String(expect).toUpperCase()}`));
}

function renderPaper() {
  const p = P.paper;
  const m = P.dom.main;
  m.dataset.state = "idle";
  clear(m);
  const meta = p.meta || {};
  const d = P.dom;
  d.runBtn = h("button.run-key.pp-reproduce#pp-reproduce", { type: "button", title: "Run the composed engine on this paper's device and fill the model column", on: { click: () => reproduce() } },
    icon("run", { size: 18 }), h("span", "REPRODUCE"));
  d.composerBtn = h("button.btn.btn-outline.pp-composer", { type: "button", title: "Open this device in the Composer at the central value of every ranged leaf", on: { click: (e) => openInComposer(e.currentTarget) } },
    icon("design", { size: 16 }), "Open in Composer");
  d.progress = h("div.progress", { role: "progressbar", "aria-label": "Reproduce progress", "aria-valuemin": "0", "aria-valuemax": "100", hidden: true }, h("span.progress-fill"));
  d.status = h("p.run-status.pp-status", { "aria-live": "polite", role: "status" });
  d.summary = h("div.pp-summary");
  d.err = h("div.pp-err");
  d.tableHost = h("div.pp-table-host");
  d.charts = h("div.pp-charts");

  const doi = meta.doi ? h("a.pp-doi", { href: `https://doi.org/${meta.doi}`, target: "_blank", rel: "noopener noreferrer", title: "Open the paper at doi.org" }, `doi:${meta.doi}`) : null;
  const arxiv = meta.arxiv ? h("a.pp-doi", { href: `https://arxiv.org/abs/${String(meta.arxiv).replace(/v\d+$/, "")}`, target: "_blank", rel: "noopener noreferrer" }, `arXiv:${meta.arxiv}`) : null;

  m.append(
    h("section.panel.pp-device", { "aria-label": "Device" },
      h("header.pp-dev-head",
        h("div.pp-dev-titles",
          h("div.pp-dev-line", h("h2.pp-dev-title", titleOf(p.id, meta.cite)), readBadge(meta.read_level),
            h("span.id-chip.is-static", { title: "Catalog ID" }, h("span.id-k", "ID"), h("span.id-v", `FS-paper-${p.id}`))),
          h("p.pp-dev-cite", meta.cite || ""),
          h("p.pp-dev-links", doi, doi && arxiv ? " · " : null, arxiv)),
        h("div.pp-actions",
          h("div.pp-status-col", d.progress, d.status),
          d.composerBtn, d.runBtn)),
      meta.summary ? h("p.pp-summary-text", meta.summary) : null,
      h("div.pp-kinds.pp-dev-kinds", AXES.map((a) => p.kinds?.[a.key] ? h("span.chip.pp-kind", { title: `${a.label}: ${p.kinds[a.key]}` }, h("span.pp-kind-k", a.label), p.kinds[a.key].replace(/_/g, " ")) : null),
        classCounts(p.classes))),
    d.err, d.summary,
    h("section.panel.pp-table-panel", { "aria-label": "Model versus paper" },
      h("header.panel-head", h("h2.engrave", `Model vs paper · ${p.reproduce.length} item${p.reproduce.length === 1 ? "" : "s"}`),
        h("span.panel-note", "ranges are [min, max] over the ranged inputs, never means")),
      d.tableHost),
    h("section.pp-charts-sec", { "aria-label": "Series overlays" }, d.charts));
  renderTable();
  renderCharts();
}

function reportedOf(id) {
  return P.paper.reported.find((r) => r.id === id) || null;
}
function seriesOf(id) {
  return P.paper.series.find((s) => s.id === id) || null;
}

function allTags() {
  return Object.values(P.paper.resolved || {}).map((r) => r?.tag);
}
/** Widest tag of the leaves an observable depends on (widest wins); falls back to the whole device. */
function modelTag(inputs) {
  const t = widest((inputs || []).map((k) => P.paper.resolved?.[k]?.tag));
  return t || widest(allTags());
}

function rows() {
  const ledger = P.paper.reproduce;
  const res = P.result?.items;
  if (!Array.isArray(res)) return ledger.map((e) => ({ ledger: e, id: e.item, cls: e.class, expect: e.expect || "pass", item: null }));
  const byId = new Map(res.map((i) => [i.id, i]));
  const out = ledger.map((e) => ({ ledger: e, id: e.item, cls: e.class, expect: e.expect || "pass", item: byId.get(e.item) || null }));
  for (const i of res) if (!out.some((r) => r.id === i.id)) out.push({ ledger: null, id: i.id, cls: i.cls, expect: i.expect, item: i });
  return out;
}

function renderTable() {
  const host = P.dom.tableHost;
  clear(host);
  const rs = rows();
  if (!rs.length) { host.append(h("p.pp-hint.pad", "This card has an empty reproduce ledger.")); return; }
  const th = (t, extra = {}) => h("th", { scope: "col", ...extra }, t);
  host.append(h("div.table-wrap.pp-table-wrap", h("table.data-table.pp-table",
    h("thead", h("tr", th("Item"), th("Paper"), th("Model"), th("Tolerance"), th("Outcome"), th("Inputs from paper"))),
    h("tbody", rs.map(row)))));
}

function paperCell(r) {
  const rep = reportedOf(r.id);
  const ser = seriesOf(r.id);
  if (r.cls === "P" || (!rep && !ser)) {
    return h("td.pp-paper", { dataset: { label: "Paper" } }, h("span.pp-muted", "no paper value"));
  }
  if (ser) {
    const tag = normTag(ser.tag);
    return h("td.pp-paper", { dataset: { label: "Paper" } }, h("div.pp-val",
      h("span.pp-num", `${ser.points.length} points`), tagChip(tag, { title: ser.where, size: "sm" }),
      ser.digitized ? h("span.label-chip.pp-dig", { title: ser.method || "digitised from a figure" }, "digitised") : null),
    h("p.pp-where", { title: ser.where }, ser.where));
  }
  const tag = normTag(rep.tag);
  const err = isNum(rep.err) ? h("span.pp-err-v", { title: "paper-quoted uncertainty" }, ` ± ${fmt(rep.err)}`) : null;
  return h("td.pp-paper", { dataset: { label: "Paper" } }, h("div.pp-val", h("span.pp-num", fmt(rep.value)), err, tagChip(tag, { title: rep.where, size: "sm" })),
    h("p.pp-where", { title: rep.where }, rep.where));
}

function bracket(lo, hi) {
  return h("span.pp-range", { title: "model [min, max] over the ranged inputs: an envelope, never a mean" }, `[${fmt(lo)}, ${fmt(hi)}]`);
}

/** Model value with ONE tag chip; the difference from the paper sits right under it, in absolute and relative terms. */
function modelCell(r) {
  const it = r.item;
  const td = (...kids) => h("td.pp-model", { dataset: { label: "Model" } }, ...kids);
  if (!it) return td(h("span.pp-muted", "not run"));
  const ser = it.kind === "series" ? P.result?.model_series?.find((s) => s.id === it.id) : null;
  const tag = normTag(ser?.tag) || modelTag(it.inputs);
  if (it.outcome === "unavailable" || it.outcome === "error" || !isNum(it.model)) {
    return td(h("span.pp-muted", "n/a"), it.note ? h("p.pp-note", it.note) : null);
  }
  if (!tag) return td(h("span.pp-muted", "n/a"), h("p.pp-note", "no provenance tag; number withheld"));
  if (it.kind === "series") {
    return td(h("div.pp-val", h("span.pp-k", "rms rel"), h("span.pp-num", fmt(it.model)), tagChip(tag, { title: "widest tag of the model curve", size: "sm" })),
      isNum(it.lo) && isNum(it.hi) ? h("p.pp-note", "curve ", bracket(it.lo, it.hi)) : null);
  }
  const ranged = isNum(it.lo) && isNum(it.hi) && Math.abs(it.hi - it.lo) > 1e-9 * Math.max(Math.abs(it.hi), Math.abs(it.lo), 1e-300);
  let delta = null;
  if (it.kind === "scalar" && isNum(it.paper)) {
    const d = it.model - it.paper;
    const sign = d > 0 ? "+" : d < 0 ? "−" : "";
    const rel = it.paper !== 0 ? ` (${sign}${fmt(Math.abs(d / it.paper) * 100, { sig: 2 })}%)` : "";
    delta = h("p.pp-delta", { title: "model minus paper at the central point; the percentage is relative to the paper value" }, h("span.pp-k", "Δ"), h("span.pp-num.pp-delta-v", `${sign}${fmt(Math.abs(d))}${rel}`));
  }
  return td(h("div.pp-val", h("span.pp-num", fmt(it.model)), tagChip(tag, { title: "widest tag of the inputs this observable depends on", size: "sm" })),
    delta, ranged ? h("p.pp-note", bracket(it.lo, it.hi)) : null);
}

function tolCell(r) {
  const it = r.item;
  const tol = tolHuman(it?.tol_str || tolText(r.ledger?.tol));
  const why = it?.tol_reason || r.ledger?.tol_reason || "";
  return h("td.pp-tolc", { dataset: { label: "Tolerance" } }, tol ? h("span.pp-tol", { title: why }, tol) : h("span.pp-muted", "none stated"));
}

function tolText(t) {
  if (!t || typeof t !== "object") return "";
  if (t.metric) return `${t.metric}<=${t.max}`;
  return Object.entries(t).map(([k, v]) => (k === "k_sigma" ? `${v}sigma` : k === "rel" ? `rel${Number(v) * 100}%` : k === "abs" ? `abs${v}` : `${k} ${v}`)).join("|");
}

function outcomeCell(r) {
  const it = r.item;
  if (!it) {
    return h("td.pp-outcome-cell", { dataset: { label: "Outcome" } }, h("span.pp-expect", { title: "the outcome the card records before any run" }, `ledger expects ${String(r.expect).toUpperCase()}`),
      r.ledger?.reason ? h("p.pp-note", r.ledger.reason) : null);
  }
  const mismatch = it.outcome !== it.expect;
  const why = [["Ledger reason", it.reason], ["Tolerance reason", it.tol_reason], ["Engine note", it.note]].filter(([, v]) => v);
  return h("td.pp-outcome-cell", { dataset: { match: mismatch ? "0" : "1", label: "Outcome" } },
    mismatch ? mismatchBadge(it.outcome, it.expect) : outcomeBadge(it.outcome),
    it.weak && it.outcome === "pass" ? h("span.label-chip.pp-weak", { title: it.weak_reason || "weak pass", dataset: { weak: "1" } }, icon("warn", { size: 12 }), " WEAK PASS") : null,
    it.weak && it.outcome === "pass" && it.weak_reason ? h("p.pp-note.pp-weak-why", it.weak_reason) : null,
    !mismatch && it.expect !== "pass" ? h("p.pp-note", `as the card expected (${it.expect})`) : null,
    why.length ? h("details.pp-more", h("summary", "why"), why.map(([k, v]) => h("p.pp-why", h("span.pp-why-k", k), v))) : null);
}

function inputsCell(r) {
  const it = r.item;
  const set = new Set(P.paper.inputs_from_paper || []);
  const td = (...kids) => h("td.pp-inputs", { dataset: { label: "Inputs from paper" } }, ...kids);
  if (!it) return td(h("span.pp-muted", "run to list"));
  const ins = it.inputs || [];
  if (!ins.length) return td(h("span.pp-muted", it.outcome === "unavailable" ? "no model run" : "none"));
  const nPaper = ins.filter((k) => set.has(k)).length;
  const near = ins.length > 0 && nPaper === ins.length;
  return td(h("details.pp-in", { dataset: { fromPaper: String(nPaper), total: String(ins.length) } },
    h("summary", `${nPaper} of ${ins.length} from the paper`),
    near ? h("p.pp-note.pp-near", icon("warn", { size: 12 }), " every input comes from this paper: the item is near-transcription") : null,
    h("ul.pp-in-list", ins.map((k) => {
      const fp = set.has(k);
      return h("li", { dataset: { fromPaper: fp ? "1" : "0" } }, leafNode(k), tagChip(P.paper.resolved?.[k]?.tag, { size: "sm", title: fp ? "taken from the paper" : "not from the paper" }),
        fp ? h("span.label-chip.pp-fp", "paper") : h("span.pp-muted.pp-notfp", "not from paper"));
    }))));
}

function row(r) {
  const it = r.item;
  const rep = reportedOf(r.id);
  const ser = seriesOf(r.id);
  const obs = it?.observable || rep?.observable || ser?.y?.observable || r.ledger?.observable || "";
  const conv = it?.convention || rep?.convention || ser?.y?.convention || r.ledger?.convention || "intrinsic";
  return h("tr", { dataset: { item: r.id, cls: r.cls, outcome: it?.outcome || "", state: it ? "run" : "ledger" } },
    h("th.pp-item", { scope: "row" },
      h("span.pp-item-obs", obsNode(obs, conv)),
      h("span.pp-item-id", { title: r.id }, itemLabel(r.id)),
      h("span.pp-item-meta", classBadge(r.cls), h("span.pp-item-kind", r.cls === "P" ? "prediction" : r.cls === "R" ? "paper reference" : ser ? "series" : "scalar"))),
    paperCell(r), modelCell(r), tolCell(r),
    outcomeCell(r), inputsCell(r));
}

// ---------------------------------------------------------------- summary + charts
function renderSummary() {
  const host = P.dom.summary;
  clear(host);
  const res = P.result;
  if (!res || !Array.isArray(res.items)) return;
  const items = res.items;
  const n = items.length;
  const nMatch = items.filter((i) => i.match).length;
  const ledgerOk = n > 0 && nMatch === n;
  const { tested, agree, inRange, weak } = agreement(items);
  const cnt = (o) => items.filter((i) => i.outcome === o).length;
  const rep = res.report || {};
  const status = !tested ? "info" : agree === tested ? "pass" : "fail";
  const word = !tested ? "NO TESTED ITEMS" : `MODEL AGREES ON ${agree} OF ${tested} TESTED ITEM${tested === 1 ? "" : "S"}`;
  const notTests = items.length - tested;
  const calib = items.filter((i) => i.cls === "C").length, pred = items.filter((i) => i.cls === "P").length, trans = items.filter((i) => i.cls === "T").length;
  const refs = items.filter((i) => i.cls === "R").length;
  host.append(h("section.plaque.pp-plaque-sum", { class: `plaque-${status}`, dataset: { status, agree: String(agree), tested: String(tested), inRange: String(inRange), weak: String(weak) }, "aria-label": "Reproduce summary" },
    h("div.plaque-top",
      h("div.plaque-word", icon(status === "pass" ? "check" : status === "fail" ? "cross" : "info", { size: 18 }), h("span", word)),
      h("div.plaque-denom", h("span.k", "tested items agree"), h("span.v", tested ? `${agree} / ${tested}` : "0 / 0"))),
    h("div.plaque-reasons",
      h("span.plaque-kv", h("span.k", "tested"), h("span.v", "unfitted (M) items, central value within tolerance")),
      inRange ? h("span.plaque-kv", { dataset: { role: "in-range" } }, h("span.k", "in range only"), h("span.v", `${inRange} tested item${inRange === 1 ? "" : "s"}: paper value inside the input envelope, central value misses (not counted)`)) : null,
      weak ? h("span.plaque-kv", { dataset: { role: "weak" } }, h("span.k", "weak passes"), h("span.v", `${weak} of the ${agree} agreeing item${agree === 1 ? "" : "s"} cannot discriminate (flagged in the card)`)) : null,
      notTests ? h("span.plaque-kv", h("span.k", "not tests"), h("span.v", [calib ? `${calib} calibrated` : "", trans ? `${trans} transcribed` : "", pred ? `${pred} predicted` : "", refs ? `${refs} paper reference${refs === 1 ? "" : "s"}` : ""].filter(Boolean).join(" · "))) : null,
      h("span.plaque-kv", h("span.k", "outcomes"), h("span.v", `${cnt("pass")} pass · ${cnt("pass-envelope")} in range only · ${cnt("fail")} fail · ${cnt("unavailable")} unavailable${cnt("error") ? ` · ${cnt("error")} error` : ""}`))),
    // the ledger is a secondary line: it says the run reproduces what the card recorded, failures included
    h("p.pp-ledger-line", { class: ledgerOk ? "" : "is-off", dataset: { ledger: ledgerOk ? "ok" : "off" } },
      icon(ledgerOk ? "check" : "warn", { size: 14 }),
      h("span", h("b", `Ledger: ${nMatch}/${n} outcomes as the card expected.`), ledgerOk
        ? " The run gives the results the card recorded, recorded failures included. This is a regression check, not agreement with the paper."
        : " At least one outcome differs from the card's record: that flip needs a human review before the card changes.")),
    rep.valid === false && rep.invalid_reasons?.length
      ? h("p.plaque-reason", icon("warn", { size: 14 }), ` Run flagged invalid: ${rep.invalid_reasons.join("; ")}`) : null,
    h("p.plaque-source", `${P.resultCached ? "cached result for this card and engine" : "fresh run"} · card sha256 ${String(res.card_sha256 || "").slice(0, 12)}${rep.solver ? ` · solver ${rep.solver}` : ""}`),
    rep.notes?.length ? h("details.pp-more.pp-notes", h("summary", `Engine notes (${rep.notes.length})`), h("ul", rep.notes.map((t) => h("li", t)))) : null));
}

function renderCharts() {
  disposeCharts();
  const host = P.dom.charts;
  clear(host);
  const ser = P.paper.series || [];
  if (!ser.length) {
    host.append(h("div.pp-nochart", { dataset: { role: "no-series" } }, icon("info", { size: 14 }),
      h("span", " This card reports no series: its numbers are the scalar items above, so there is no curve to overlay.")));
    return;
  }
  const set = new Set(P.paper.inputs_from_paper || []);
  const grid = h("div.pp-chart-grid");
  host.append(h("h2.engrave.pp-sec-title", `Series overlays · ${ser.length}`), grid);
  for (const s of ser) {
    const model = P.result?.model_series?.find((m) => m.id === s.id) || null;
    const item = P.result?.items?.find((i) => i.id === s.id);
    const ov = paperOverlay({ series: s, model, inputs: item ? { paths: item.inputs || [], fromPaper: set } : null, filename: `fsim-paper-${P.paper.meta.id}-${s.id}` });
    P.charts.push(ov.chart);
    grid.append(ov.el);
    ov.render();
  }
}

// ---------------------------------------------------------------- reproduce job
function setRunState(state, text = "") {
  const d = P.dom;
  if (!d.runBtn) return;
  P.dom.main.dataset.state = state;
  d.runBtn.disabled = state === "running";
  d.runBtn.setAttribute("aria-busy", String(state === "running"));
  d.status.textContent = text;
  d.progress.hidden = state !== "running";
  if (state !== "running") d.progress.removeAttribute("aria-valuenow");
}

function tick() {
  if (!P?.dom.progress) return;
  const el = (performance.now() - P.t0) / 1000;
  const f = Math.min(0.95, el / Math.max(1, P.eta));
  const pct = Math.round(f * 100);
  P.dom.progress.firstChild.style.transform = `scaleX(${f})`;
  P.dom.progress.setAttribute("aria-valuenow", String(pct));
  P.dom.status.textContent = `Running the reproduce ledger · about ${Math.max(1, Math.round(P.eta))} s expected · ${el.toFixed(0)} s elapsed`;
}

async function reproduce({ auto = false } = {}) {
  if (!P?.paper || P.running) return;
  const seq = P.seq;
  const id = P.id;
  P.running = true;
  clear(P.dom.err);
  if (!auto) setRunState("running", "Submitting the reproduce job");
  let resp;
  try {
    resp = await P.ctx.api.postJSON(`/api/papers/${encodeURIComponent(id)}/reproduce`, {});
  } catch (e) {
    if (!P || seq !== P.seq) return;
    return failed(e);
  }
  if (!P || seq !== P.seq) return;
  try {
    let result = resp.result;
    if (!result) {
      P.eta = resp.eta_s || 5;
      P.t0 = performance.now();
      setRunState("running", "Queued");
      P.dom.progress.firstChild.style.transform = "scaleX(0)";
      clearInterval(P.timer);
      P.timer = setInterval(tick, 250);
      P.job = P.ctx.api.awaitJob(resp.job_id, { onState: (s) => { if (P && s && s !== "running") P.dom.status.textContent = `Job ${s}`; } });
      result = await P.job.promise;
    }
    if (!P || seq !== P.seq) return;
    clearInterval(P.timer);
    P.timer = null;
    P.job = null;
    P.running = false;
    P.result = result;
    P.resultCached = !!resp.cached;
    if (Array.isArray(result?.items)) P.agree[id] = agreement(result.items);
    setRunState("done", `${result.items?.length ?? 0} items · ${resp.cached ? "cached result" : "fresh run"}${resp.run_id ? ` · ${resp.run_id}` : ""}`);
    P.dom.main.dataset.state = "done";
    renderSummary();
    renderTable();
    renderCharts();
    refreshListQuiet();
  } catch (e) {
    if (!P || seq !== P.seq) return;
    failed(e);
  }
}

function failed(e) {
  clearInterval(P.timer);
  P.timer = null;
  P.job = null;
  P.running = false;
  setRunState("error", "");
  P.dom.main.dataset.state = "error";
  const msg = e?.body?.error?.message || e?.message || "The reproduce job failed";
  clear(P.dom.err).append(errorState({ title: "Reproduce failed", message: msg, actions: [{ label: "Try again", fn: () => reproduce(), primary: true }] }));
}

async function refreshListQuiet() {
  try {
    const l = await P.ctx.api.getJSON("/api/papers");
    if (!P) return;
    P.list = l;
    renderCatalogue();
  } catch { /* the catalogue keeps its old badges */ }
}

// ---------------------------------------------------------------- Composer hand-off
async function openInComposer(btn) {
  if (!P?.paper) return;
  const id = P.id;
  btn.setAttribute("aria-busy", "true");
  try {
    const r = await P.ctx.api.getJSON(`/api/papers/${encodeURIComponent(id)}/design`);
    const seed = { design: r.design, paper: id, collapsed: r.collapsed || [], note: r.note || "" };
    P.ctx.state.set({ composeSeed: seed });
    try { sessionStorage.setItem("fsim-studio:compose-seed", JSON.stringify(seed)); } catch { /* private window: state.js carries it */ }
    P.ctx.navigate("#/compose");
  } catch (e) {
    toast(`Could not open in the Composer: ${e.message}`, { kind: "error" });
  } finally {
    btn.removeAttribute("aria-busy");
  }
}
