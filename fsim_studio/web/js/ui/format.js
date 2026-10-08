// Display formatting only (rounding for the screen; the server sends full precision).

const UNIT_MAP = { "1": "", uA: "µA", um: "µm", um2: "µm²", cm2: "cm⁻²", "1/cm2": "cm⁻²", s: "s", "/s": "s⁻¹", meV: "meV", K: "K" };

export function unitText(u) {
  if (u == null) return "";
  return UNIT_MAP[u] ?? String(u).replace(/\^-1/g, "⁻¹").replace(/\^-2/g, "⁻²").replace(/\^-3/g, "⁻³").replace(/\^2/g, "²").replace(/\^3/g, "³");
}

export function isNum(v) { return typeof v === "number" && Number.isFinite(v); }

/** 4 significant figures; grouped integers for big values; scientific when tiny/huge. */
export function fmt(v, { sig = 4 } = {}) {
  if (v == null) return "n/a";
  if (typeof v === "boolean") return v ? "true" : "false";
  if (typeof v !== "number") return String(v);
  if (!Number.isFinite(v)) return "n/a";
  if (v === 0) return "0";
  const a = Math.abs(v);
  if (a >= 1e6 || a < 1e-3) {
    const [m, e] = v.toExponential(sig - 1).split("e");
    return `${m}×10${superscript(Number(e))}`;
  }
  if (a >= 1000) return Math.round(v).toLocaleString("en-US");
  return Number(v.toPrecision(sig)).toString();
}

const SUP = { "-": "⁻", 0: "⁰", 1: "¹", 2: "²", 3: "³", 4: "⁴", 5: "⁵", 6: "⁶", 7: "⁷", 8: "⁸", 9: "⁹" };
function superscript(n) { return String(n).split("").map((c) => SUP[c] ?? c).join(""); }

/** Interval "lo – hi" (en dash, thin spaces). */
export function fmtInterval(lo, hi, opts) {
  if (!isNum(lo) && !isNum(hi)) return "n/a";
  return `${isNum(lo) ? fmt(lo, opts) : "n/a"} – ${isNum(hi) ? fmt(hi, opts) : "n/a"}`;
}

export function fmtSeconds(s) {
  if (!isNum(s)) return "";
  if (s < 1) return `${Math.max(1, Math.round(s * 1000))} ms`;
  if (s < 90) return `${s.toFixed(s < 10 ? 1 : 0)} s`;
  return `${Math.round(s / 60)} min`;
}

// Display names for result scalars and curves (labels only, no physics).
export const SCALAR_LABELS = {
  g2_op: "g²(0) at operating point",
  T_c: "Ceiling temperature T_c",
  T_j_op: "Junction temperature T_j",
  dT_J: "Self-heating ΔT_J",
  brightness_per_pulse: "Brightness per pulse",
  collected_flux_pulsed_s: "Collected flux, commanded",
  collected_flux_delivered_s: "Collected flux, delivered",
  eps_op: "XX leakage ε",
  rho_op: "Loading purity ρ",
  gamma_op: "Linewidth Γ",
};

export const SCALAR_UNITS = {
  T_c: "K", T_j_op: "K", dT_J: "K", gamma_op: "meV", collected_flux_pulsed_s: "s⁻¹",
  collected_flux_delivered_s: "s⁻¹", rep_rate_hz: "Hz", tau_pulse_ns: "ns", V_j_op: "V", V_j: "V",
  V_applied: "V", V_bi: "V", P_junction_W: "W", edge_lambda_nm: "nm", edge_A_mode_um2: "µm²",
  w_resolved: "meV", cavity_detuning_op_meV: "meV", aperture_lambda_op: "nm",
};
