// User-facing names for the machine keys that data files, campaigns and scalar dictionaries carry
// (FIN3, review item 5 and 6). Labels only: no value is computed here.
//   fieldLabel(key)   -> rich string ("Minimum $g^{(2)}(0)$"): feed it to h() / richText(), or richPlain() for <option>
//   fieldUnit(key)    -> display unit text ("s⁻¹") or ""
//   humanValue(v, key)-> display text of a printed value: numbers to 4 significant figures with the key's unit,
//                        counts with spaces around the slash, identifiers as words, booleans as yes / no
//   campaignTitle(id) -> "RT edge emitter" for out/rt_edge
//   words(id)         -> "deterministic_pair" -> "deterministic pair"
//   cardTitle(name)   -> "chatzarakis2023" -> "Chatzarakis 2023" (Library card names; hyphens of an id are kept)
//   designTitle(name) -> "edge-inp-gainp-design" -> "Edge InP GaInP design" (scene titles)
// The raw key stays available for a tooltip (title attribute) and for tests (data attributes), never as visible text.
import { fmt, isNum } from "./format.js";
import { displayUnit } from "./units.js";
import { symbolFor } from "./math.js";

/** "deterministic_pair" -> "deterministic pair" (identifier-shaped text only; anything else is returned unchanged). */
export const IDENT_RE = /^[A-Za-z][A-Za-z0-9]*(?:_[A-Za-z0-9]+)+$/;
export const words = (s) => (typeof s === "string" && IDENT_RE.test(s) ? s.replace(/_/g, " ") : s);

// ---------------------------------------------------------------- parameter-card names
const CARD_TOKENS = { inp: "InP", gainp: "GaInP", gaasp: "GaAsP", gaas: "GaAs", ingan: "InGaN", gan: "GaN", qw: "QW", qcap: "QCAP",
  set: "SET", cw: "CW", rt: "RT", hom: "HOM", nm: "nm" };
/** Human name of a parameter card: "chatzarakis2023" -> "Chatzarakis 2023". The id stays in a title attribute. */
export function cardTitle(name) {
  const t = words(String(name)).replace(/([A-Za-z])(\d)/g, "$1 $2").replace(/(\d)([A-Za-z])/g, "$1 $2");
  return t.replace(/^./, (c) => c.toUpperCase());
}
/** Display name of a design card: "edge-inp-gainp-design" -> "Edge InP GaInP design" (material names in their usual case). */
export function designTitle(name) {
  const t = String(name).split(/[-_]+/).filter(Boolean).map((tok) => {
    const m = /^([A-Za-z]+)(\d.*)$/.exec(tok);              // chatzarakis2023 -> chatzarakis 2023; 4K -> 4 K
    if (m) return `${CARD_TOKENS[m[1].toLowerCase()] ?? m[1]} ${m[2]}`;
    const m2 = /^(\d+)([A-Za-z]+)$/.exec(tok);
    if (m2) return `${m2[1]} ${m2[2]}`;
    return CARD_TOKENS[tok.toLowerCase()] ?? tok;
  }).join(" ");
  return t.replace(/^./, (c) => c.toUpperCase());
}

/** Rich title of a SceneSpec: the card id at its head becomes the card's display name, "XX -> X -> 0" is typeset. */
export function sceneTitle(spec) {
  let t = String(spec?.title || spec?.kind || "");
  const card = spec?.card;
  if (card && t.startsWith(`${card} - `)) t = `${designTitle(card)} - ${t.slice(card.length + 3)}`;
  else if (card && t === card) t = designTitle(card);
  return t.replace("XX -> X -> 0", "$XX \\to X \\to 0$");
}

// ---------------------------------------------------------------- campaigns (folders under out/)
const CAMPAIGN_TITLES = {
  rt_edge: "RT edge emitter",
  nitride_cavity: "Planar InGaN cavity",
  nitride_geometry_stark: "Nitride geometry and Stark",
  "nitride_nanowire/full": "Nitride nanowires",
  rt_campaign: "RT campaign (early sweeps)",
  tier_geometry: "Tier sweep: geometry",
  tier_device: "Tier sweep: device",
  phase0: "Phase 0: validation fit",
  phase1: "Phase 1: sweeps",
  phase2: "Phase 2: sweeps",
  phase3: "Phase 3: sweeps",
  spec: "Specification checks",
  validation: "Validation",
  zhao: "QD laser (Zhao)",
};
export function campaignTitle(id, fallback = null) {
  if (CAMPAIGN_TITLES[id]) return CAMPAIGN_TITLES[id];
  const t = fallback && fallback !== id ? fallback : String(id).split("/").pop();
  return words(t).replace(/^./, (c) => c.toUpperCase());
}

// ---------------------------------------------------------------- field dictionary
// [label, unit]. Labels with `$...$` are typeset. Keys of VERDICT / BEST lines, sweep CSV columns and design scalars.
const F = {
  // VERDICT / BEST lines
  regime: ["Regime"], model: ["Model"], complete: ["Complete"], evidence: ["Evidence"], conditional: ["Conditional"],
  g2_min: ["Minimum $g^{(2)}(0)$"], g2_median_eligible: ["Median $g^{(2)}(0)$, eligible rows"],
  diag_g2_min: ["Minimum $g^{(2)}(0)$, diagnostic"], diag_g2_median_diagnostic: ["Median $g^{(2)}(0)$, diagnostic"],
  diag_g2_flux_max: ["$g^{(2)}(0)$ at the maximum flux, diagnostic"],
  flux_max: ["Maximum collected flux", "s⁻¹"], flux_margin: ["Flux margin (maximum over floor)"],
  flux_floor: ["Flux floor", "s⁻¹"], flux_shortfall_deprecated: ["Flux shortfall (deprecated)"],
  flux_floor_excluded: ["Rows below the flux floor"],
  eligible: ["Eligible rows"], eligible_fraction: ["Eligible fraction"], invalid: ["Invalid rows"],
  eligible_dedup: ["Eligible rows, deduplicated"], eligible_dedup_mismatch_groups: ["Deduplication mismatch groups, eligible"],
  headline_dedup_mismatch_groups: ["Deduplication mismatch groups, headline"],
  coverage: ["Coverage"], coverage_over_eligible: ["Coverage of the eligible rows"],
  headline_coverage: ["Headline coverage"], headline_coverage_pulsed: ["Headline coverage, pulsed"],
  headline_by_T: ["Headline coverage by $T_{\\mathrm{hs}}$"], headline_by_T_pulsed: ["Headline coverage by $T_{\\mathrm{hs}}$, pulsed"],
  rows_scheduled: ["Rows scheduled"], cw_raw_coverage: ["CW raw coverage"],
  median_pass: ["Median passes"], T_pass_min: ["Lowest passing $T_{\\mathrm{hs}}$", "K"],
  gamma300_pass_max: ["Largest passing $\\Gamma$ at 300 K", "meV"], gamma300_threshold: ["$\\Gamma$ threshold at 300 K", "meV"],
  idealized_status: ["Idealized status"], idealized_pass_count: ["Idealized passes"],
  hardware_infeasible_count: ["Hardware-infeasible rows"], hardware_pass_count: ["Hardware passes"],
  paired_optical_pass: ["Paired optical passes"], quality_pass: ["Quality passes"],
  hardware_qualified: ["Hardware qualified"], rti_qualified: ["RTI qualified"],
  family: ["Family"], screening: ["Screening"], access: ["Dot access"],
  strain_bound: ["Strain bound"], bound_role: ["Bound role"], bound_reversal: ["Bound reversal"],
  rep_rate_hz: ["Repetition rate", "Hz"], value: ["Value"], row_id: ["Row"], T_hs: ["$T_{\\mathrm{hs}}$", "K"],
  core_radius_nm: ["Core radius", "nm"], height_nm: ["Dot height", "nm"], radius_nm: ["Dot radius", "nm"], x_in: ["Indium fraction $x_{\\mathrm{in}}$"],
  commanded_flux: ["Commanded flux", "s⁻¹"], delivered_flux: ["Delivered flux", "s⁻¹"],
  photons_per_cycle_commanded: ["Photons per cycle, commanded"], photons_per_cycle_delivered: ["Photons per cycle, delivered"],
  headline_eligible: ["Headline eligible"], headline_eligible_candidates: ["Headline-eligible candidates"],
  optical_pass_candidates: ["Optical-pass candidates"], blocked_load_probability: ["Blocked-load probability"],
  unrelaxed_partner_row_id: ["Unrelaxed partner row"], unrelaxed_partner_flux: ["Unrelaxed partner flux", "s⁻¹"],
  // sweep CSV columns not in charts/sweep.js
  invalid_reasons: ["Reasons the row is invalid"], invalid_reasons_pulsed: ["Reasons the row is invalid (pulsed)"],
  config_id: ["Configuration"], diagnostic_valid: ["Diagnostic valid"], flux_below_floor: ["Below the flux floor"],
  card_file: ["Card file"], occupied_dot_access: ["Occupied-dot access"], row_kind: ["Row kind"],
  model_finite_pulse: ["Finite-pulse model"], model_tau_cap_density: ["Density-capped lifetime model"],
  delta_xx_tag: ["$\\Delta_{XX}$ tag"], gamma300_tag: ["$\\Gamma$(300 K) tag"], irf_tag: ["IRF tag"], card_class: ["Card class"],
  g2_q: ["Quantum $g^{(2)}(0)$"], deterministic_pair: ["Deterministic pair"],
  report_only_runtime_s: ["Run time (report only)", "s"], runtime_s: ["Run time", "s"], runtime_seconds: ["Run time", "s"],
  generation_worktree_dirty: ["Worktree had uncommitted changes"], report_only_generation_worktree_dirty: ["Worktree had uncommitted changes (report only)"],
  card_hashes: ["Card hashes"], output_hashes: ["Output hashes"], source_hashes: ["Source hashes"],
  output_hashes_excluded_downstream: ["Output hashes, excluded downstream"],
  // design scalars (Designer "All scalars", result tables)
  mean_counts: ["Mean counts per pulse"], T_j_op: ["Junction temperature $T_{\\mathrm{j}}$", "K"], dT_J: ["Self-heating $\\Delta T_{\\mathrm{J}}$", "K"],
  eps_op: ["XX leakage $\\varepsilon$"], g2_op: ["$g^{(2)}(0)$ at the operating point"], t_x_op: ["X transmission through the filter"],
  brightness_per_pulse: ["Brightness per pulse"], brightness_convention: ["Brightness convention"],
  g2_op_valid: ["$g^{(2)}(0)$ valid"], g2_op_invalid_reason: ["Why $g^{(2)}(0)$ is not valid"],
  F_eff_wire_op: ["Effective Purcell factor, wire"], cavity_detuning_op_meV: ["Cavity detuning", "meV"],
  rho_g2_op: ["Loading purity $\\rho$"], linewidth_source: ["Linewidth source"], drive_source: ["Drive source"],
  eta_capture_resolved: ["Capture efficiency, resolved"], b_e_resolved: ["Electron blocking, resolved"],
  V_j_op: ["Junction voltage $V_{\\mathrm{j}}$", "V"], P_junction_W: ["Junction power", "W"], T_j_transport: ["Junction temperature from transport", "K"],
  tau_pulse_ns: ["Pulse length", "ns"], duty_resolved: ["Duty cycle, resolved"],
  collected_flux_pulsed_s: ["Collected flux, pulsed", "s⁻¹"], collected_flux_delivered_s: ["Collected flux, delivered", "s⁻¹"],
  flux_measurable: ["Flux measurable"], emission_type: ["Emission type"], xx_in_window: ["XX inside the filter window"],
  g2_cw0: ["CW $g^{(2)}(0)$"], g2_cw0_raw: ["CW $g^{(2)}(0)$, raw"],
  aperture_g2_penalty: ["Aperture $g^{(2)}(0)$ penalty"], aperture_compose: ["Aperture composition"], aperture_lambda_op: ["Aperture wavelength", "nm"],
  cw_r_ns: ["CW pump rate", "ns⁻¹"], cw_gamma_X_ns: ["CW X decay rate", "ns⁻¹"], cw_rho_op: ["CW loading purity $\\rho$"],
  edge_lambda_nm: ["Emission wavelength", "nm"], edge_beta: ["Coupling factor $\\beta$"], edge_T_facet: ["Facet transmission"],
  edge_eta_prop: ["Propagation efficiency"], edge_eta_NA: ["NA collection efficiency"], edge_eta_total: ["Total collection efficiency"],
  edge_n_eff: ["Effective index"], edge_n_g: ["Group index"], edge_Gamma_dot: ["Dot confinement factor $\\Gamma$"], edge_A_mode_um2: ["Mode area", "µm²"],
  finite_pulse: ["Finite-pulse model"], finite_pulse_g2_dot: ["Dot $g^{(2)}(0)$, finite pulse"],
  finite_pulse_mean_counts: ["Mean counts, finite pulse"], finite_pulse_mean_counts_x: ["Mean X counts, finite pulse"],
  finite_pulse_r_ns: ["Pump rate, finite pulse", "ns⁻¹"], finite_pulse_gamma_X_ns: ["X decay rate, finite pulse", "ns⁻¹"],
  finite_pulse_gamma_XX_ns: ["XX decay rate, finite pulse", "ns⁻¹"], finite_pulse_k_X: ["X filter fraction, finite pulse"],
  finite_pulse_k_XX: ["XX filter fraction, finite pulse"], finite_pulse_tau_on_ns: ["Pulse on-time", "ns"],
  finite_pulse_tau_dark_ns: ["Dark time", "ns"], finite_pulse_gate_ns_used: ["Gate used", "ns"], finite_pulse_converged: ["Finite-pulse model converged"],
  photon_budget: ["Photon budget"], retention_source: ["Retention source"], retention_params_used: ["Retention parameters used"],
  track_material: ["Track material"], track_material_effective: ["Track material, effective"],
  "aperture.diameter_um": ["Diameter", "µm"], "aperture.density_cm2": ["Dot density", "cm⁻²"],
  diameter_um: ["Aperture diameter", "µm"], density_cm2: ["Dot density", "cm⁻²"],
  // layer names of the 3D and 2D stacks
  barrier_lower: ["Lower barrier"], barrier_upper: ["Upper barrier"], matrix_lower: ["Lower matrix"], matrix_upper: ["Upper matrix"],
  // Explain / explorers
  applied_field_kVcm: ["Applied field", "kV/cm"], dot_levels: ["Dot levels"],
  dE_e_matrix: ["$\\Delta E_e$, matrix", "meV"], dE_e_barrier: ["$\\Delta E_e$, barrier", "meV"],
  dE_h_matrix: ["$\\Delta E_h$, matrix", "meV"], dE_h_barrier: ["$\\Delta E_h$, barrier", "meV"],
  ibm_transmission: ["Independent-boson transmission"], ibm_purcell_transmission: ["Independent-boson transmission with Purcell"],
  effective_gamma_zpl: ["Effective ZPL linewidth $\\Gamma_{\\mathrm{ZPL}}$"],
  // library datasets
  g2_vs_T: ["$g^{(2)}(0)$ against temperature"], gamma_vs_T: ["Linewidth against temperature"], tau_vs_T: ["Lifetime against temperature"],
  g2_vs_ERR: ["$g^{(2)}(0)$ against error"], g2_vs_pump: ["$g^{(2)}$ against pump"], g2_oe2008: ["$g^{(2)}$, Oe 2008 table"], device_points: ["Device points"],
};

const UNIT_SUFFIX = { _nm: "nm", _um: "µm", _ns: "ns", _ps: "ps", _meV: "meV", _ueV: "µeV", _K: "K", _Hz: "Hz", _hz: "Hz", _W: "W", _V: "V" };

/** Rich label of a machine key: the dictionary entry, else the words of the key, first letter capital. */
export function fieldLabel(key) {
  const k = String(key);
  const e = F[k];
  if (e) return e[0];
  const toks = k.split("_").filter(Boolean);
  for (const suf of Object.keys(UNIT_SUFFIX)) if (k.endsWith(suf) && toks.length > 1) { toks.pop(); break; }
  return toks.join(" ").replace(/\bg2\b/g, "$g^{(2)}$").replace(/\bT hs\b/g, "$T_{\\mathrm{hs}}$").replace(/^./, (c) => c.toUpperCase());
}

export const knownField = (key) => Object.prototype.hasOwnProperty.call(F, String(key));

/**
 * A /api/cards meta map whose `label` is the machine name ("diameter_um") gets a readable one: the dictionary entry
 * for the path or its last segment, else no label at all (the symbol is then derived from the registry name).
 */
export function cleanMeta(meta) {
  const out = {};
  for (const [p, m] of Object.entries(meta || {})) {
    const tail = p.split(".").pop();
    if (m && typeof m === "object" && (!m.label || m.label === tail || IDENT_RE.test(m.label))) {
      const key = knownField(p) ? p : knownField(tail) ? tail : null;
      const { label, ...rest } = m;
      out[p] = key ? { ...rest, label: fieldLabel(key) } : rest;
    } else out[p] = m;
  }
  return out;
}

export function fieldUnit(key) {
  const e = F[String(key)];
  if (e && e[1]) return e[1];
  const k = String(key);
  for (const [suf, u] of Object.entries(UNIT_SUFFIX)) if (k.endsWith(suf) && k.length > suf.length) return displayUnit(u);
  return "";
}

/** Symbol-only label (a short TeX symbol) when the registry knows one; used by dense readouts. */
export function fieldSymbol(key) { return symbolFor({ name: String(key) }); }

// ---------------------------------------------------------------- values
const NUM_RE = /^[-+]?(?:\d+\.?\d*|\.\d+)(?:[eE][-+]?\d+)?$/;
const yesNo = { true: "yes", false: "no", True: "yes", False: "no" };

/** Display text of a printed field value. Rounding and wording only: the raw text stays in the caller's tooltip. */
export function humanValue(raw, key = "") {
  if (raw == null) return "";
  if (typeof raw === "boolean") return raw ? "yes" : "no";
  if (typeof raw === "number") return isNum(raw) ? fmt(raw) : "n/a";
  const s = String(raw).trim();
  if (s === "") return "";
  if (s in yesNo) return yesNo[s];
  if (/^(nan|none|n\/a)$/i.test(s)) return /^nan$/i.test(s) ? "n/a" : s;
  if (NUM_RE.test(s)) {
    const n = Number(s);
    const k = String(key);
    if (/_hz$/i.test(k) && Math.abs(n) >= 1e6) return `${fmt(n / 1e6)} MHz`;
    if (/_hz$/i.test(k)) return `${fmt(n)} Hz`;
    const u = fieldUnit(key);
    return u ? `${fmt(n)} ${u}` : fmt(n);
  }
  let m = /^([-+]?[\d.]+(?:e[-+]?\d+)?)\s*\/\s*s$/i.exec(s);                 // "1000/s"
  if (m) return `${fmt(Number(m[1]))} s⁻¹`;
  m = /^(\d+)\s*\/\s*(\d+)$/.exec(s);                                       // "16/768": count and denominator
  if (m) return `${m[1]} / ${m[2]}`;
  if (/_by_T(_pulsed)?$/.test(String(key)) && /^\d+:/.test(s)) {            // "230:0/192,250:0/192"
    return s.split(",").map((p) => { const [t, c] = p.split(":"); return `${t} K: ${(c || "").replace(/\s*\/\s*/, " / ")}`; }).join(" · ");
  }
  if (/^[A-Za-z_]+:[A-Za-z0-9_.-]+(?:,[A-Za-z_]+:[A-Za-z0-9_.-]+)*$/.test(s)) { // "finite_pulse:true,tau_cap_density:false"
    return s.split(",").map((p) => { const [a, b] = p.split(":"); return `${a.replace(/_/g, " ")}: ${yesNo[b] ?? b}`; }).join(", ");
  }
  if (IDENT_RE.test(s)) return s.replace(/_/g, " ");
  return s;
}
