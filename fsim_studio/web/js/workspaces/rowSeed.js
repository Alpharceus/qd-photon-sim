// Sweep row -> card configuration (no physics): the rt_edge grid maps one-to-one onto card fields,
// the same overrides scripts/run_rt_edge.py applies. Shared by the Results row drawer (live
// re-run) and the Designer (state.rerunSeed from "Open card in Designer").

/** Campaigns whose rows can be applied to a card exactly. */
export const SEEDABLE = new Set(["rt_edge"]);

/** Copy of `design` with the rt_edge row's grid values and model switches applied. */
export function rtEdgeRowToDesign(design, row) {
  const d = JSON.parse(JSON.stringify(design));
  d.dot.delta_xx = row.delta_xx_meV;
  d.dot.gamma300 = row.gamma300_meV;
  d.thermal.T_hs = row.T_hs_K;
  d.drive.duty = row.duty_pulsed;
  d.drive.cw = false;
  d.drive.diode = { ...(d.drive.diode || {}), tau_pulse_ns: row.pulse_width_ns };
  d.drive.finite_pulse = row.model_finite_pulse;
  d.ret.tau_cap_scales_with_density = row.model_tau_cap_density;
  d.emission.NA = row.emission_NA;
  d.emission.R_back = row.emission_R_back;
  d.emission.L_um = row.emission_L_um;
  return d;
}

/** The fields of a row that the seed changes, for display: [[card path, value]]. */
export function rtEdgeSeedFields(row) {
  return [
    ["dot.delta_xx", row.delta_xx_meV], ["dot.gamma300", row.gamma300_meV], ["thermal.T_hs", row.T_hs_K],
    ["drive.duty", row.duty_pulsed], ["drive.diode.tau_pulse_ns", row.pulse_width_ns],
    ["drive.finite_pulse", row.model_finite_pulse], ["ret.tau_cap_scales_with_density", row.model_tau_cap_density],
    ["emission.NA", row.emission_NA], ["emission.R_back", row.emission_R_back], ["emission.L_um", row.emission_L_um],
  ];
}

/**
 * CSV vs live pairs: [label, csv value, live value, unit, scalar key]. The scalar key is the
 * result key (used for the n/a reason lookup `<key>__nan_reason`), never the display label.
 */
export function rtEdgePairs(row, scalars) {
  const s = scalars || {};
  return [
    ["pulsed g²(0)", row.g2_pulsed, s.g2_op, "", "g2_op"],
    ["collected flux", row.collected_flux_pulsed_s, s.collected_flux_pulsed_s, "/s", "collected_flux_pulsed_s"],
    ["T_j", row.Tj_pulsed_K, s.T_j_op, "K", "T_j_op"],
  ];
}

/** n/a reason for a live value: the backend's `<key>__nan_reason`, or g2_op_invalid_reason. */
export function liveNanReason(scalars, key) {
  const s = scalars || {};
  return s[`${key}__nan_reason`] || (key === "g2_op" ? s.g2_op_invalid_reason : null) || "";
}
