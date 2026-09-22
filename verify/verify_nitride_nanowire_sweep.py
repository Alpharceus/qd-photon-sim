"""Independent structural, grid, VERDICT-recomputation, mutation-sensitivity
and presentation checks for scripts/run_nitride_nanowire.py's output (piece 9).

Without --out-dir this runs entirely in memory: pure combo-builder counts
(never a real evaluate() call except the minimal smoke test at the very
end), independently-transcribed gate predicates exercised at boundary
values, and a fabricated-RT-pass detector exercised on synthetic rows --
no previously generated artifact is required, so this script is always safe
to run as part of the project's `for %F in (verify\\verify_nitride_*.py)`
test-command sweep, before any sweep output exists.

With --out-dir <dir> it additionally validates the real quick/full run
written there, READ-ONLY: row counts against the run module's own pure
build_core/build_reduced_cuts/build_dipole_falsification/
build_deshpande2013/build_planar_reference builders (an independent
re-derivation of the plan, not a re-read of the CSV's own row count),
card/output hash verification, strain-bound pairing and core-coordinate
uniqueness (a dropped bound or a duplicated coordinate fails these), every
VERDICT line's fields recomputed from sweep.csv via the SAME independently
-transcribed predicates the self-test exercises (a flipped verdict boolean
fails this), BEST_PASSING_FLUX(_300K) recomputation, figure/manifest
plot-construction-contract tracing back to sweep.csv by row_id (a merged
trace's fixed coordinates fails this), and the contract's results-text
obligations (bullet 11: RC caveat, access-1.0 lifetime cap, opposite-
endpoint anchor match, E_C/kT wall, dipole falsification, 200/80 MHz
one_pair_valid comparison).

The independently-transcribed predicates below (`optical_pass`,
`hardware_qualified`, `rti_qualified`, `eligible`) are copied from
fsim_core/nitride_nanowire_device.py's own `_one_raw` (the `optical =
bool(valid and g2 < .5 and flux >= 1000 and (regime != "deterministic_pair"
or (one and supply)))` line and the two `hardware_qualified`/
`rti_qualified` dict entries immediately below it) and from
scripts/run_nitride_nanowire.py's own `_row()` (`eligible`), NEVER by
calling the production routine being tested -- a corrupted saved
optical_pass/hardware_qualified/rti_qualified/eligible column is caught by
this recomputation, not trusted.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, re, sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import scripts.run_nitride_nanowire as m  # noqa: E402


# --------------------------------------------------------------- utilities
def rows(p):
    with open(p, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))


def f(x):
    try:
        return float(x)
    except (TypeError, ValueError):
        return float("nan")


def isfin(x):
    try:
        return math.isfinite(float(x))
    except (TypeError, ValueError):
        return False


def as_bool(x):
    """CSV round-trips Python bools as the strings 'True'/'False'; a
    genuine sentinel like 'not_applicable' stays neither True nor False."""
    if isinstance(x, bool):
        return x
    if x == "True":
        return True
    if x == "False":
        return False
    return None


def close_nan_safe(a, b, rtol=1e-6, atol=1e-9):
    fa, fb = f(a), f(b)
    na, nb = math.isnan(fa), math.isnan(fb)
    if na or nb:
        return na and nb
    return math.isclose(fa, fb, rel_tol=rtol, abs_tol=atol)


# ------------------------------------------------- independent predicates
# Transcribed verbatim from fsim_core/nitride_nanowire_device.py's _one_raw
# (optical/hardware_qualified/rti_qualified) and scripts/
# run_nitride_nanowire.py's _row() (eligible) -- see module docstring.
def optical_pass(valid, g2, flux, regime, one_pair_valid, pair_supply_possible):
    if not (valid and isfin(g2) and isfin(flux)):
        return False
    ok = valid and g2 < 0.5 and flux >= 1000
    if regime == "deterministic_pair":
        ok = ok and bool(one_pair_valid) and bool(pair_supply_possible)
    return bool(ok)


def hardware_qualified(optical, regime, set_feasible):
    return bool(optical and regime == "deterministic_pair" and set_feasible is True)


def rti_qualified(optical, regime, rti_feasible):
    return bool(optical and regime == "deterministic_pair" and rti_feasible is True)


def eligible(valid, g2, flux):
    return bool(valid and isfin(g2) and isfin(flux) and flux >= 1000.0)


# H3 fix: quality_pass is a SCRIPT-computed gate orthogonal to the
# contract's own optical_pass (one_pair_valid rewards faster emptying by
# ADDED loss just as readily as by improved device quality -- see the
# script's "Gate anti-monotonicity" results.md section). Transcribed
# independently from scripts/run_nitride_nanowire.py's own
# QUALITY_PHOTONS_PER_CYCLE_FLOOR/_attach_quality_columns, never by calling
# the production routine being tested.
QUALITY_PHOTONS_PER_CYCLE_FLOOR = 0.01


def quality_pass(optical, photons_per_cycle):
    return bool(optical and isfin(photons_per_cycle) and photons_per_cycle >= QUALITY_PHOTONS_PER_CYCLE_FLOOR)


def photons_per_cycle_of(flux, rep_rate_hz):
    if not (isfin(flux) and isfin(rep_rate_hz)) or rep_rate_hz <= 0:
        return float("nan")
    return flux / rep_rate_hz


def headline_pass(optical, family, headline_eligible):
    """H2 fix: a nomination (BEST_PASSING_FLUX(_300K) and the 16
    per-family/regime/bound/rate nominations) additionally requires
    headline_eligible -- vertical_photonic rows above the LP11 single-mode
    cutoff are never nominated even if optical_pass."""
    if not optical:
        return False
    if family == "vertical_photonic":
        return headline_eligible is True
    return True


def rti_fabricated_pass(second_pair_control_known, rti_feasible):
    """Contract Composition-rules bullet 7: 'unsupported second-pair/hole-
    occupation controls produce an honest failed or unknown RT screen
    (rti_feasible=False, rti_status), never a fabricated pass.' Returns True
    (a FINDING, not a pass) iff a row claims rti_feasible=True despite
    second_pair_control_known=False."""
    return (second_pair_control_known is False) and (rti_feasible is True)


CORE_COORD_KEYS = ("family", "regime", "rep_rate_hz", "core_radius_nm", "height_nm", "x_in", "T_hs")


def core_coord(row, key_fn=lambda r, k: r.get(k)):
    return tuple(key_fn(row, k) for k in CORE_COORD_KEYS)


def find_bound_partner(core_rows, row):
    """Independent (non-imported) re-implementation of run_nitride_nanowire.
    _bound_partner: same family/regime/rate/radius/height/x_in/T_hs, the
    OPPOSITE strain_bound. Operates on plain dicts (CSV string rows or
    synthetic fixtures) via close_nan_safe for the float coordinates."""
    other = "relaxed" if row.get("strain_bound") == "unrelaxed" else "unrelaxed"
    for cand in core_rows:
        if cand is row or cand.get("strain_bound") != other:
            continue
        if (cand.get("family") == row.get("family") and cand.get("regime") == row.get("regime")
                and close_nan_safe(cand.get("rep_rate_hz"), row.get("rep_rate_hz"))
                and close_nan_safe(cand.get("core_radius_nm"), row.get("core_radius_nm"))
                and close_nan_safe(cand.get("height_nm"), row.get("height_nm"))
                and close_nan_safe(cand.get("x_in"), row.get("x_in"))
                and close_nan_safe(cand.get("T_hs"), row.get("T_hs"))):
            return cand
    return None


def duplicate_core_coords(core_rows):
    """Returns the set of (family,regime,rate,radius,height,x_in,T_hs,
    strain_bound) tuples that appear more than once in `core_rows` --
    catches a duplicated coordinate (acceptance criterion 2)."""
    seen, dupes = {}, set()
    for r in core_rows:
        key = core_coord(r) + (r.get("strain_bound"),)
        key = tuple(round(v, 6) if isinstance(v, float) else v for v in key)
        seen[key] = seen.get(key, 0) + 1
    for k, n in seen.items():
        if n > 1:
            dupes.add(k)
    return dupes


def trace_shares_fixed_coords(entry, all_rows, x_key):
    """Given one plot-contract trace entry ({"row_ids": [...], "x": [...],
    "y": [...]}), confirm every referenced row agrees on every core
    coordinate OTHER than x_key -- a merged trace (two different fixed
    geometries plotted as one line) fails this."""
    fixed_keys = [k for k in CORE_COORD_KEYS if k != x_key]
    seen_fixed = None
    for rid in entry.get("row_ids", []):
        row = all_rows.get(rid)
        if row is None:
            return False
        vals = tuple(round(f(row.get(k)), 6) if k != "family" and k != "regime" else row.get(k) for k in fixed_keys)
        if seen_fixed is None:
            seen_fixed = vals
        elif vals != seen_fixed:
            return False
    return True


# --------------------------------------------------------- self-test (A)
def check_predicate_boundaries(checks):
    checks.append(("optical_pass: g2 exactly at gate (0.5) is rejected", optical_pass(True, 0.5, 2000., "rectangular", True, True) is False))
    checks.append(("optical_pass: g2 just under gate (0.4999) with flux at floor is accepted (rectangular)", optical_pass(True, 0.4999, 1000., "rectangular", True, True) is True))
    checks.append(("optical_pass: flux just under floor (999.999) is rejected", optical_pass(True, 0.1, 999.999, "rectangular", True, True) is False))
    checks.append(("optical_pass: an invalid row is never optical_pass regardless of g2/flux", optical_pass(False, 0.01, 1e9, "rectangular", True, True) is False))
    checks.append(("optical_pass: deterministic_pair with one_pair_valid=False is rejected even with good g2/flux", optical_pass(True, 0.1, 2000., "deterministic_pair", False, True) is False))
    checks.append(("optical_pass: deterministic_pair with pair_supply_possible=False is rejected even with good g2/flux", optical_pass(True, 0.1, 2000., "deterministic_pair", True, False) is False))
    checks.append(("optical_pass: deterministic_pair with one_pair_valid=True and pair_supply_possible=True is accepted", optical_pass(True, 0.1, 2000., "deterministic_pair", True, True) is True))
    checks.append(("optical_pass: rectangular regime never requires one_pair_valid/pair_supply_possible", optical_pass(True, 0.1, 2000., "rectangular", False, False) is True))
    checks.append(("hardware_qualified: rectangular regime is never hardware_qualified regardless of set_feasible", hardware_qualified(True, "rectangular", True) is False))
    checks.append(("hardware_qualified: deterministic_pair with set_feasible=False is never hardware_qualified", hardware_qualified(True, "deterministic_pair", False) is False))
    checks.append(("hardware_qualified: deterministic_pair with set_feasible=True and optical_pass is hardware_qualified", hardware_qualified(True, "deterministic_pair", True) is True))
    checks.append(("hardware_qualified: a non-optical_pass row is never hardware_qualified", hardware_qualified(False, "deterministic_pair", True) is False))
    checks.append(("rti_qualified: rectangular regime is never rti_qualified regardless of rti_feasible", rti_qualified(True, "rectangular", True) is False))
    checks.append(("rti_qualified: deterministic_pair with rti_feasible=False is never rti_qualified", rti_qualified(True, "deterministic_pair", False) is False))
    checks.append(("rti_qualified: deterministic_pair with rti_feasible=True and optical_pass is rti_qualified", rti_qualified(True, "deterministic_pair", True) is True))
    checks.append(("eligible: flux exactly at floor (1000) with finite g2 is eligible", eligible(True, 0.9, 1000.) is True))
    checks.append(("eligible: flux just under floor (999.999) is not eligible", eligible(True, 0.01, 999.999) is False))
    checks.append(("eligible: invalid row is never eligible", eligible(False, 0.01, 1e9) is False))
    checks.append(("eligible: NaN g2 is never eligible", eligible(True, float("nan"), 2000.) is False))
    checks.append(("quality_pass: photons_per_cycle exactly at the 0.01 floor with optical_pass is accepted",
                    quality_pass(True, 0.01) is True))
    checks.append(("quality_pass: photons_per_cycle just under the 0.01 floor is rejected",
                    quality_pass(True, 0.009999) is False))
    checks.append(("quality_pass: a non-optical_pass row is never quality_pass regardless of photons_per_cycle",
                    quality_pass(False, 100.0) is False))
    checks.append(("quality_pass: NaN photons_per_cycle is never quality_pass", quality_pass(True, float("nan")) is False))
    checks.append(("photons_per_cycle_of: flux/rate division matches a direct value", photons_per_cycle_of(2.0e6, 200.0e6) == 0.01))
    checks.append(("photons_per_cycle_of: non-finite rate is nan-safe", math.isnan(photons_per_cycle_of(2.0e6, 0.0))))
    checks.append(("headline_pass: horizontal_as_built never requires headline_eligible",
                    headline_pass(True, "horizontal_as_built", False) is True))
    checks.append(("headline_pass: vertical_photonic with headline_eligible=False is rejected even if optical_pass "
                    "(H2 fix: never nominate a brighter multimode row)", headline_pass(True, "vertical_photonic", False) is False))
    checks.append(("headline_pass: vertical_photonic with headline_eligible=True and optical_pass is accepted",
                    headline_pass(True, "vertical_photonic", True) is True))
    checks.append(("headline_pass: a non-optical_pass row is never headline_pass", headline_pass(False, "vertical_photonic", True) is False))


def check_rti_fabrication_fixture(checks):
    checks.append(("mutation fixture: second_pair_control_known=False + rti_feasible=True IS flagged as a fabricated RT pass",
                    rti_fabricated_pass(False, True) is True))
    checks.append(("mutation fixture: second_pair_control_known=False + rti_feasible=False is NOT flagged",
                    rti_fabricated_pass(False, False) is False))
    checks.append(("mutation fixture: second_pair_control_known=True + rti_feasible=True is NOT flagged (supported control)",
                    rti_fabricated_pass(True, True) is False))


def check_bound_partner_and_duplicate_fixtures(checks):
    """Synthetic (no artifacts) exercise of find_bound_partner/
    duplicate_core_coords -- directly demonstrates the 'dropping a bound'
    and 'duplicating a coordinate' failure modes acceptance criterion 2
    names, without needing a real corrupted sweep.csv."""
    base = dict(family="horizontal_as_built", regime="deterministic_pair", rep_rate_hz=80.0e6,
                core_radius_nm=12.5, height_nm=2.0, x_in=0.40, T_hs=300.0)
    relaxed = dict(base, strain_bound="relaxed", row_id="A")
    unrelaxed = dict(base, strain_bound="unrelaxed", row_id="B")
    checks.append(("find_bound_partner locates the opposite-strain-bound partner when both rows are present",
                    find_bound_partner([relaxed, unrelaxed], relaxed) is unrelaxed))
    checks.append(("find_bound_partner returns None when the partner bound is DROPPED (acceptance criterion 2)",
                    find_bound_partner([relaxed], relaxed) is None))
    dupes = duplicate_core_coords([relaxed, unrelaxed, dict(relaxed)])
    checks.append(("duplicate_core_coords flags a DUPLICATED coordinate (acceptance criterion 2)",
                    len(dupes) == 1))
    checks.append(("duplicate_core_coords reports no duplicates on a clean (no-repeat) row set",
                    len(duplicate_core_coords([relaxed, unrelaxed])) == 0))
    # Merged-trace-fixed-coordinates fixture.
    rowmap = {"A": dict(relaxed, core_radius_nm=10.0), "B": dict(relaxed, core_radius_nm=20.0, height_nm=3.0)}
    good_entry = {"row_ids": ["A"], "x": [10.0], "y": [1.0]}
    bad_entry = {"row_ids": ["A", "B"], "x": [10.0, 20.0], "y": [1.0, 2.0]}
    checks.append(("trace_shares_fixed_coords accepts a trace whose non-swept coordinates all agree",
                    trace_shares_fixed_coords(good_entry, rowmap, "core_radius_nm") is True))
    checks.append(("trace_shares_fixed_coords rejects a trace MERGING two different fixed heights (acceptance criterion 2)",
                    trace_shares_fixed_coords(bad_entry, rowmap, "core_radius_nm") is False))


def check_grid_self(checks, detail):
    horiz = m.build_core("horizontal_as_built", False)
    vert = m.build_core("vertical_photonic", False)
    checks.append(("full-mode horizontal core grid is exactly 1536 rows (contract literal)", len(horiz) == 1536))
    checks.append(("full-mode vertical core grid is exactly 1024 rows (contract literal)", len(vert) == 1024))
    checks.append(("full-mode core grid totals 2560 rows (contract literal)", len(horiz) + len(vert) == 2560))
    ax = ("core_radius_nm", "height_nm", "x_in", "T_hs", "regime", "strain_bound", "rep_rate_hz")
    combos = {tuple(p[k] for k in ax) for p in horiz}
    checks.append(("full-mode horizontal grid has no duplicate axis combination", len(combos) == len(horiz)))
    combos_v = {tuple(p[k] for k in ax) for p in vert}
    checks.append(("full-mode vertical grid has no duplicate axis combination", len(combos_v) == len(vert)))
    checks.append(("full-mode horizontal core radii match the contract literal set",
                    set(m.CORE_R_NM["horizontal_as_built"]) == {10.0, 12.5, 15.0, 20.0, 25.0, 40.0}))
    checks.append(("full-mode vertical core radii match the contract literal set",
                    set(m.CORE_R_NM["vertical_photonic"]) == {60.0, 80.0, 100.0, 120.0}))
    checks.append(("both regimes present in the module's REGIMES", set(m.REGIMES) == {"rectangular", "deterministic_pair"}))
    checks.append(("both strain bounds present in the module's STRAIN_BOUNDS", set(m.STRAIN_BOUNDS) == {"unrelaxed", "relaxed"}))
    checks.append(("both rep rates present (80/200 MHz core axis)", set(m.REP_RATES) == {80.0e6, 200.0e6}))
    checks.append(("quick mode never claims full coverage (fewer horizontal core rows than full)",
                    len(m.build_core("horizontal_as_built", True)) < 1536))
    checks.append(("quick mode still samples both regimes/bounds/rates",
                    {p["regime"] for p in m.build_core("horizontal_as_built", True)} == {"rectangular", "deterministic_pair"}
                    and {p["strain_bound"] for p in m.build_core("horizontal_as_built", True)} == {"unrelaxed", "relaxed"}
                    and {p["rep_rate_hz"] for p in m.build_core("horizontal_as_built", True)} == {80.0e6, 200.0e6}))

    # H6 radius-coupling invariants (DEVICE-PIECE CONSTRAINT 1), checked on
    # real (in-memory) DeviceDesign objects built by the run module's own
    # _design() -- no evaluate() call.
    ph = dict(horiz[0]); ph.update(core_radius_nm=20.0, regime="deterministic_pair")
    dh, _ = m._design(ph)
    checks.append(("horizontal_as_built: disc radius equals core radius (pre-H6 convention retained)",
                    dh.nitride["dot"]["radius_nm"] == 20.0))
    checks.append(("horizontal_as_built: SET radius equals core radius", dh.drive.set_params["radius_nm"] == 20.0))
    checks.append(("horizontal_as_built: conducting radius equals core radius", dh.drive.diode["conducting_radius_nm"] == 20.0))

    pv = dict(vert[0]); pv.update(core_radius_nm=100.0, regime="deterministic_pair")
    dv, _ = m._design(pv)
    checks.append(("vertical_photonic: disc radius is STRICTLY LESS than core radius (H6)",
                    dv.nitride["dot"]["radius_nm"] < dv.nitride["nanowire"]["core_radius_nm"]))
    checks.append(("vertical_photonic: disc radius stays at its own 12.5 nm default despite the 100 nm core sweep (H6)",
                    dv.nitride["dot"]["radius_nm"] == 12.5))
    checks.append(("vertical_photonic: SET radius follows the DISC radius, never the core radius (Composition rules bullet 1)",
                    dv.drive.set_params["radius_nm"] == dv.nitride["dot"]["radius_nm"] == 12.5))
    checks.append(("vertical_photonic: conducting radius follows the core radius (drive.diode only)",
                    dv.drive.diode["conducting_radius_nm"] == 100.0))

    detail.append("grid_self: horizontal=%d vertical=%d core_total=%d" % (len(horiz), len(vert), len(horiz) + len(vert)))


REQUIRED_CUT_AXES = {
    "screening_fraction": {0.0, 1.0},
    "occupied_dot_access": {1.0},
    "S_cm_s": {1.0e2, 1.0e4},
    "shell": {"AlGaN"},
    "al_fraction": {0.2},
    "growth_tolerance_steps": {2.0},
    "alignment_uncertainty_meV": {30.0},
    "occupation_control_uncertainty": {True},
    "injector_barrier_thickness_nm": {1.0},
    "R_s_ohm": {1.0e6},
}


def check_reduced_cut_axes_declared(checks, detail):
    cuts = m.build_reduced_cuts(False)
    checks.append(("full-mode reduced cuts is non-empty", len(cuts) > 0))
    by_axis = {}
    for p in cuts:
        by_axis.setdefault(p["sensitivity_axis"], set()).add(p["sensitivity_value"])
    for axis, expected_values in REQUIRED_CUT_AXES.items():
        got = by_axis.get(axis, set())
        checks.append((f"reduced-cut axis '{axis}' is present with its declared alternative value(s) {expected_values}",
                        expected_values <= got))
    checks.append(("R_s_ohm reduced cut is evaluated on the horizontal family only (spec: 'R_s designed 1e6 on the horizontal family')",
                    {p["family"] for p in cuts if p.get("sensitivity_axis") == "R_s_ohm"} == {"horizontal_as_built"}))
    checks.append(("every declared-dropped axis in DROPPED_CUT_AXES is actually absent from the full-mode cuts",
                    all(name.split(" ")[0].split("(")[0].strip() not in by_axis for name in
                        ("reservoir_access", "gamma300", "tau_rad0_ns", "tau_cap_ps", "C_parasitic_F", "NA", "bottom_reflectivity"))))
    checks.append(("the RESTORED current/pulse-width axis 'current_pulse_width' is not in DROPPED_CUT_AXES "
                    "(M9-M10 fix: 'restore the current / pulse-width cut')",
                    all("current_pulse_width" not in name and "current_uA" not in name and "tau_pulse_ns" not in name
                        for name in m.DROPPED_CUT_AXES)))
    pulse_cuts = [p for p in cuts if p.get("sensitivity_axis") == "current_pulse_width"]
    checks.append(("current_pulse_width reduced cut is present with exactly 9 rows (I in {0.001,0.002,0.02} uA "
                    "x tau_pulse in {0.01,0.1,1} ns, Cartesian, horizontal reference geometry)",
                    len(pulse_cuts) == 9
                    and {p["I_uA"] for p in pulse_cuts} == {0.001, 0.002, 0.02}
                    and {p["tau_pulse_ns"] for p in pulse_cuts} == {0.01, 0.1, 1.0}
                    and {p["family"] for p in pulse_cuts} == {"horizontal_as_built"}))
    detail.append("reduced_cut_rows(full)=%d axes=%s" % (len(cuts), sorted(by_axis)))


def check_dry_run_cap(checks):
    counts = m.planned_counts(False)
    checks.append(("full-mode planned evaluate-call upper bound stays within the 10000 cap",
                    counts["planned_evaluate_calls_upper_bound"] <= 10000))
    counts_q = m.planned_counts(True)
    checks.append(("quick-mode planned evaluate-call upper bound stays within the 600 s quick command's practical budget "
                    "(<=1000 calls)", counts_q["planned_evaluate_calls_upper_bound"] <= 1000))
    checks.append(("N_WORKERS is a positive integer (process-pool sizing)", isinstance(m.N_WORKERS, int) and m.N_WORKERS >= 1))


def check_evaluate_smoke(checks, detail):
    """The 'minimal evaluate smoke test': exactly two REAL evaluate() calls
    (one per family) at the reference geometry, confirming the pipeline
    returns a well-shaped, valid scalar dict end to end. Not a sweep."""
    ph = m.full_defaults("horizontal_as_built")
    dh, _ = m._design(ph)
    sh = m.device.evaluate(dh)["scalars"]
    checks.append(("horizontal smoke evaluate() returns valid=True at reference geometry", sh.get("valid") is True))
    checks.append(("horizontal smoke evaluate() reports family=horizontal_as_built", sh.get("family") == "horizontal_as_built"))
    checks.append(("horizontal smoke evaluate() carries optical_pass/hardware_qualified/rti_qualified keys",
                    all(k in sh for k in ("optical_pass", "hardware_qualified", "rti_qualified"))))

    pv = m.full_defaults("vertical_photonic")
    dv, _ = m._design(pv)
    sv = m.device.evaluate(dv)["scalars"]
    checks.append(("vertical smoke evaluate() returns valid=True at reference geometry", sv.get("valid") is True))
    checks.append(("vertical smoke evaluate() reports family=vertical_photonic", sv.get("family") == "vertical_photonic"))
    checks.append(("vertical smoke evaluate() carries headline_eligible key", "headline_eligible" in sv))
    detail.append("smoke: horizontal valid=%s g2_op=%s | vertical valid=%s g2_op=%s"
                   % (sh.get("valid"), sh.get("g2_op"), sv.get("valid"), sv.get("g2_op")))


# --------------------------------------------------------- out-dir checks
REQUIRED_ROW_COLUMNS = (
    "platform", "family", "T_hs", "T_j", "cycle_loading", "rep_rate_hz", "strain_bound", "bound_role",
    "screening_fraction", "core_radius_nm", "outer_radius_nm", "conducting_radius_nm",
    "E_X_eV", "lambda_nm", "field_kVcm", "overlap_sq", "gamma_X0_ns", "gamma_XX0_ns", "k_X_ns", "k_XX_ns",
    "V_number", "beta_HE11", "eta_collection_X", "eta_collection_XX", "gamma_X_ns", "gamma_XX_ns",
    "single_mode", "degree_of_linear_polarization", "antenna_rate_factor",
    "k_side_ns", "k_surface_X_ns", "k_surface_XX_ns",
    "area_cm2", "J_A_cm2", "V_j", "eta_inj", "f_capture", "r_supply_s", "r_captured_s", "mu",
    "tau_RC_ns", "delivered_step_fraction", "pulse_delivery_feasible",
    "tau_rad_bare_ns", "tau_rad_photonic_ns", "tau_total_X_ns",
    "collected_flux_pulsed_s", "collected_flux_delivered_s", "g2_op", "one_pair_valid",
    "pair_supply_possible", "valid", "invalid_reasons",
    "set_feasible", "set_E_C_meV", "set_EC_over_kT",
    "rti_feasible", "rti_status", "rti_level_margin_kT",
    "optical_pass", "hardware_qualified", "rti_qualified", "headline_eligible",
    "row_id", "row_kind", "card_file", "card_hash", "cache_hit", "cache_identity", "eligible",
    "photons_per_cycle", "emission_probability_per_cycle", "quality_pass",
)


def check_files_present(out, checks):
    for name in ("sweep.csv", "manifest.json", "results.md"):
        checks.append((f"{name} exists", (out / name).is_file()))
    pngs = list(out.glob("*.png"))
    checks.append(("at least one figure PNG exists", len(pngs) > 0))
    for p in pngs:
        checks.append((f"figure is non-trivial (>1000 bytes): {p.name}", p.stat().st_size > 1000))


def check_row_counts(core, man, quick, checks, detail):
    core_rows = [r for r in core if r["row_kind"] == "core"]
    sens_rows = [r for r in core if r["row_kind"] == "sensitivity"]
    planar_rows = [r for r in core if r["row_kind"] == "planar_reference"]
    replay_rows = [r for r in core if r["row_kind"] == "planar_2014_replay"]

    exp_h = len(m.build_core("horizontal_as_built", quick))
    exp_v = len(m.build_core("vertical_photonic", quick))
    checks.append(("core row count matches independent build_core() for both families",
                    len(core_rows) == exp_h + exp_v))
    if not quick:
        checks.append(("full run's core row count is exactly 2560 (contract literal)", len(core_rows) == 2560))

    exp_cuts = len(m.build_reduced_cuts(quick)) + len(m.build_dipole_falsification())
    checks.append(("sensitivity row count matches independent build_reduced_cuts()+build_dipole_falsification()",
                    len(sens_rows) == exp_cuts + len(m.build_deshpande2013(quick))))
    exp_planar = len(m.build_planar_reference(quick))
    checks.append(("planar_reference row count matches independent build_planar_reference()", len(planar_rows) == exp_planar))
    checks.append(("exactly one planar_2014_replay row", len(replay_rows) == 1))

    all_ids = [r["row_id"] for r in core]
    checks.append(("all row_ids are unique", len(set(all_ids)) == len(all_ids)))
    checks.append(("manifest row_counts.core matches sweep.csv", man.get("row_counts", {}).get("core") == len(core_rows)))
    checks.append(("manifest row_counts.total matches sweep.csv", man.get("row_counts", {}).get("total") == len(core)))
    detail.append("row_counts: core=%d sensitivity=%d planar_reference=%d planar_2014_replay=%d total=%d"
                   % (len(core_rows), len(sens_rows), len(planar_rows), len(replay_rows), len(core)))


def check_column_set(core, checks):
    nanowire_rows = [r for r in core if r["row_kind"] in ("core", "sensitivity")]
    checks.append(("nanowire rows are present to check the column set on", len(nanowire_rows) > 0))
    missing = set()
    for col in REQUIRED_ROW_COLUMNS:
        if not all(col in r for r in nanowire_rows[:1]):
            missing.add(col)
    # csv.DictWriter with restval="" over a sorted key union means every row
    # (once written) carries every column that ANY row produced; check the
    # header itself carries every required column.
    header = set(nanowire_rows[0].keys()) if nanowire_rows else set()
    missing = {c for c in REQUIRED_ROW_COLUMNS if c not in header}
    checks.append((f"every required nanowire row column is present in sweep.csv's header (missing: {sorted(missing)})",
                    len(missing) == 0))


def check_card_hashes(core, man, checks):
    ok = True
    for name in set(m.CARD_NAME.values()) | set(m.PLANAR_CARD_NAME.values()) | {m.DESHPANDE2014_CARD}:
        real = hashlib.sha256((ROOT / "cards" / name).read_bytes()).hexdigest()
        if man.get("card_hashes", {}).get(name) != real:
            ok = False
    checks.append(("manifest card_hashes match sha256 of the actual card files on disk", ok))
    row_ok = True
    for r in core:
        if r["row_kind"] in ("core", "sensitivity"):
            expected = m.CARD_NAME.get((r.get("family"), r.get("regime")))
        elif r["row_kind"] == "planar_2014_replay":
            expected = m.DESHPANDE2014_CARD
        else:
            expected = None
        if expected is not None and r.get("card_file") != expected:
            row_ok = False
        if expected is not None and r.get("card_hash") != hashlib.sha256((ROOT / "cards" / expected).read_bytes()).hexdigest():
            row_ok = False
    checks.append(("every core/sensitivity/2014-replay row's card_file/card_hash matches the real on-disk card it claims", row_ok))


def check_bound_pairing_and_duplicates(core, checks):
    core_rows = [r for r in core if r["row_kind"] == "core"]
    dupes = duplicate_core_coords(core_rows)
    checks.append(("no core row shares its full coordinate (incl. strain_bound) with another core row "
                    "(would indicate a duplicated coordinate)", len(dupes) == 0))
    unpaired = [r for r in core_rows if find_bound_partner(core_rows, r) is None]
    checks.append(("every core row has its opposite-strain-bound partner present in this run "
                    "(a dropped bound would leave a row unpaired)", len(unpaired) == 0))


VERDICT_RE = re.compile(
    r"^VERDICT: idealized_status=(?P<ideal>\S+) family=(?P<family>\S+) regime=(?P<regime>\S+) "
    r"strain_bound=(?P<sb>\S+) bound_role=(?P<role>\S+) rep_rate_hz=(?P<rate>\S+) "
    r"complete=(?P<complete>\S+) eligible=(?P<elig>\d+) paired_optical_pass=(?P<paired>\d+) "
    r"quality_pass=(?P<quality>\d+) "
    r"hardware_qualified=(?P<hw>\d+) rti_qualified=(?P<rti>\d+) coverage=(?P<covn>\d+)/(?P<covd>\d+) "
    r"invalid=(?P<invalid>\d+) flux_floor=1000/s screening=(?P<screening>\S+) access=(?P<access>\S+)$")


def check_verdict_lines(text, core, quick, checks, detail):
    verdicts = [ln for ln in text.splitlines() if ln.startswith("VERDICT:")]
    checks.append(("results.md contains VERDICT lines", len(verdicts) > 0))
    checks.append(("exactly 16 VERDICT lines (family x regime x strain_bound x rep_rate_hz)", len(verdicts) == 16))
    core_rows = [r for r in core if r["row_kind"] == "core"]
    for family in m.FAMILIES:
        for regime in m.REGIMES:
            for sb in m.STRAIN_BOUNDS:
                for rate in m.REP_RATES:
                    checks.append((f"VERDICT line present: family={family} regime={regime} strain_bound={sb} rep_rate_hz={rate:g}",
                                    any(f"family={family} " in v and f"regime={regime} " in v
                                        and f"strain_bound={sb} " in v and f"rep_rate_hz={rate:g} " in v for v in verdicts)))
    n_parsed = 0
    for v in verdicts:
        mo = VERDICT_RE.match(v)
        if mo is None:
            checks.append((f"VERDICT line parses in the documented field order: {v[:90]}", False))
            continue
        n_parsed += 1
        family, regime, sb, rate_s = mo["family"], mo["regime"], mo["sb"], mo["rate"]
        rate = f(rate_s)
        group = [r for r in core_rows if r["family"] == family and r["regime"] == regime
                 and r["strain_bound"] == sb and close_nan_safe(r["rep_rate_hz"], rate)]
        expected = [p for p in m.build_core(family, quick)
                    if p["regime"] == regime and p["strain_bound"] == sb and p["rep_rate_hz"] == rate]

        n_invalid = n_elig = n_paired = n_hw = n_rti = n_quality = 0
        for r in group:
            valid = as_bool(r["valid"]) is True
            g2v, fluxv = f(r.get("g2_op")), f(r.get("collected_flux_pulsed_s"))
            one = as_bool(r.get("one_pair_valid")) is True
            supply = as_bool(r.get("pair_supply_possible")) is True
            opt = optical_pass(valid, g2v, fluxv, regime, one, supply)
            elig = eligible(valid, g2v, fluxv)
            n_elig += int(elig)
            if not valid:
                n_invalid += 1
            if opt and find_bound_partner(core_rows, r) is not None:
                n_paired += 1
            set_feas = as_bool(r.get("set_feasible"))
            rti_feas = as_bool(r.get("rti_feasible"))
            n_hw += int(hardware_qualified(opt, regime, set_feas))
            n_rti += int(rti_qualified(opt, regime, rti_feas))
            # H3 fix: quality_pass recomputed from flux/rep_rate_hz (raw
            # columns), never trusting the CSV's own photons_per_cycle/
            # quality_pass columns.
            ppc = photons_per_cycle_of(fluxv, f(r.get("rep_rate_hz")))
            n_quality += int(quality_pass(opt, ppc))

        ok = (int(mo["elig"]) == n_elig and int(mo["paired"]) == n_paired and int(mo["quality"]) == n_quality
              and int(mo["hw"]) == n_hw and int(mo["rti"]) == n_rti and int(mo["invalid"]) == n_invalid
              and int(mo["covn"]) == len(group) and int(mo["covd"]) == len(expected))
        checks.append((f"VERDICT counts recomputed from sweep.csv match the line "
                        f"({family}/{regime}/{sb}/{rate:g}Hz): eligible={n_elig} paired={n_paired} "
                        f"quality={n_quality} hw={n_hw} rti={n_rti} invalid={n_invalid} "
                        f"coverage={len(group)}/{len(expected)}", ok))
        complete_expected = bool((not quick) and len(group) == len(expected))
        checks.append((f"VERDICT complete flag is consistent with coverage ({family}/{regime}/{sb}/{rate:g}Hz)",
                        (mo["complete"] == "True") == complete_expected or quick))
    detail.append("verdict lines parsed: %d/%d" % (n_parsed, len(verdicts)))


BEST_FLUX_RE = re.compile(r"^BEST_PASSING_FLUX family=(?P<family>\S+) value=(?P<value>\S+) row_id=(?P<row_id>\S+)")
BEST_FLUX_300K_RE = re.compile(r"^BEST_PASSING_FLUX_300K family=(?P<family>\S+) value=(?P<value>\S+) row_id=(?P<row_id>\S+)")


def _best_passing_flux(core_rows, family, restrict_300k, require_headline=True):
    """H2 fix: headline selection requires optical_pass AND (for
    vertical_photonic) headline_eligible -- the brightest optical_pass row
    is never nominated if it is headline-ineligible (single_mode False /
    above the LP11 cutoff)."""
    cand = []
    for r in core_rows:
        if r["family"] != family:
            continue
        if restrict_300k and not close_nan_safe(r.get("T_hs"), 300.0):
            continue
        valid = as_bool(r["valid"]) is True
        g2v, fluxv = f(r.get("g2_op")), f(r.get("collected_flux_pulsed_s"))
        one = as_bool(r.get("one_pair_valid")) is True
        supply = as_bool(r.get("pair_supply_possible")) is True
        opt = optical_pass(valid, g2v, fluxv, r["regime"], one, supply)
        if require_headline and not headline_pass(opt, family, as_bool(r.get("headline_eligible"))):
            continue
        if opt and isfin(fluxv):
            cand.append((fluxv, r["row_id"]))
    return max(cand, key=lambda z: z[0]) if cand else (None, None)


def check_headline_selection_fixture(checks):
    """H2 fix, mutation-sensitivity fixture: a synthetic row set where the
    BRIGHTEST optical_pass row is deliberately headline-ineligible
    (single_mode=False / above the LP11 cutoff); _best_passing_flux must
    never select it, only the dimmer headline_eligible row."""
    bright_ineligible = dict(family="vertical_photonic", regime="deterministic_pair", T_hs="230.0",
                              valid="True", g2_op="0.17", collected_flux_pulsed_s="9000000",
                              one_pair_valid="True", pair_supply_possible="True",
                              headline_eligible="False", row_id="FIXA")
    dim_eligible = dict(bright_ineligible, collected_flux_pulsed_s="5000000",
                         headline_eligible="True", row_id="FIXB")
    rows = [bright_ineligible, dim_eligible]
    best_val, best_rid = _best_passing_flux(rows, "vertical_photonic", restrict_300k=False)
    checks.append(("_best_passing_flux prefers a DIMMER headline_eligible row over a BRIGHTER "
                    "headline_eligible=False row (H2 fix)", best_rid == "FIXB" and best_val == 5000000.0))
    best_val_no_req, best_rid_no_req = _best_passing_flux(rows, "vertical_photonic", restrict_300k=False,
                                                            require_headline=False)
    checks.append(("_best_passing_flux WITHOUT the headline requirement would have picked the brighter "
                    "ineligible row (confirms the fixture actually exercises the filter)",
                    best_rid_no_req == "FIXA"))


def check_best_passing_flux(text, core, checks):
    core_rows = [r for r in core if r["row_kind"] == "core"]
    for family in m.FAMILIES:
        lines_found = [ln for ln in text.splitlines() if ln.startswith("BEST_PASSING_FLUX ") or ln.startswith("BEST_PASSING_FLUX family")]
        mo = next((BEST_FLUX_RE.match(ln) for ln in lines_found if f"family={family} " in ln), None)
        best_val, best_rid = _best_passing_flux(core_rows, family, restrict_300k=False)
        if best_val is None:
            checks.append((f"BEST_PASSING_FLUX family={family} matches independent recomputation (no passing rows)",
                            mo is not None and mo["value"] == "none"))
        else:
            checks.append((f"BEST_PASSING_FLUX family={family} row_id/value match independent recomputation "
                            f"(expected {best_rid}, {best_val:.6g})",
                            mo is not None and mo["row_id"] == best_rid and close_nan_safe(mo["value"], best_val, rtol=1e-4)))

        lines_300 = [ln for ln in text.splitlines() if ln.startswith("BEST_PASSING_FLUX_300K")]
        mo300 = next((BEST_FLUX_300K_RE.match(ln) for ln in lines_300 if f"family={family} " in ln), None)
        best_val3, best_rid3 = _best_passing_flux(core_rows, family, restrict_300k=True)
        if best_val3 is None:
            checks.append((f"BEST_PASSING_FLUX_300K family={family} matches independent recomputation (no passing rows at 300K)",
                            mo300 is not None and mo300["value"] == "none"))
        else:
            checks.append((f"BEST_PASSING_FLUX_300K family={family} row_id/value match independent recomputation "
                            f"(expected {best_rid3}, {best_val3:.6g})",
                            mo300 is not None and mo300["row_id"] == best_rid3 and close_nan_safe(mo300["value"], best_val3, rtol=1e-4)))


def check_invalid_preserved(core, checks):
    invalid_rows = [r for r in core if as_bool(r.get("valid")) is False]
    checks.append(("at least some rows are invalid (screening/E_C/etc. constraints are expected to bite somewhere)",
                    len(invalid_rows) > 0))
    bad = []
    for r in invalid_rows:
        try:
            reasons = json.loads(r.get("invalid_reasons", "[]"))
        except json.JSONDecodeError:
            bad.append(r["row_id"]); continue
        if not isinstance(reasons, list) or len(reasons) == 0:
            bad.append(r["row_id"])
    checks.append(("every invalid row carries a non-empty, decodable invalid_reasons list (invalid rows preserved with reasons)",
                    len(bad) == 0))


def check_rti_fabrication_real(core, checks):
    nanowire_rows = [r for r in core if r["row_kind"] in ("core", "sensitivity") and r.get("regime") == "deterministic_pair"]
    fab = []
    for r in nanowire_rows:
        spck = as_bool(r.get("second_pair_control_known"))
        rti = as_bool(r.get("rti_feasible"))
        if spck is None:
            continue
        if rti_fabricated_pass(spck, rti):
            fab.append(r["row_id"])
    checks.append(("no deterministic_pair row in this run fabricates an RT pass (rti_feasible=True) despite "
                    "second_pair_control_known=False (Composition rules bullet 7)", len(fab) == 0))


def check_composition_counts(text, core, checks):
    """H1 fix fixture: independently recompute per-(family,strain_bound,
    x_in) optical_pass counts from sweep.csv (this file's OWN
    independently-transcribed optical_pass() predicate above, never the
    production routine) and cross-check them against the exact counts
    results.md prints -- a regression guard against the fixed false claim
    "these are the ONLY x_in=0.40 optical passes in either family this
    run" (64 vertical relaxed x_in=0.40 rows pass too)."""
    for fam in ("horizontal_as_built", "vertical_photonic"):
        n_rel25 = n_rel40 = n_unrel40 = 0
        for r in core:
            if r.get("row_kind") != "core" or r.get("family") != fam:
                continue
            valid = as_bool(r.get("valid")) is True
            g2v, fluxv = f(r.get("g2_op")), f(r.get("collected_flux_pulsed_s"))
            one = as_bool(r.get("one_pair_valid")) is True
            supply = as_bool(r.get("pair_supply_possible")) is True
            if not optical_pass(valid, g2v, fluxv, r.get("regime"), one, supply):
                continue
            xin = f(r.get("x_in"))
            if r.get("strain_bound") == "relaxed":
                if close_nan_safe(xin, 0.25):
                    n_rel25 += 1
                elif close_nan_safe(xin, 0.40):
                    n_rel40 += 1
            elif r.get("strain_bound") == "unrelaxed" and close_nan_safe(xin, 0.40):
                n_unrel40 += 1
        n_rel = n_rel25 + n_rel40
        mo = re.search(rf"\({re.escape(fam)}\) RELAXED optical passes: (\d+) total \(x_in=0\.25: (\d+), "
                        rf"x_in=0\.40: (\d+)\)", text)
        line_ok = bool(mo) and int(mo.group(1)) == n_rel and int(mo.group(2)) == n_rel25 and int(mo.group(3)) == n_rel40
        checks.append((f"results.md's {fam} RELAXED composition line matches an independent recount "
                        f"(n_rel={n_rel}, x_in0.25={n_rel25}, x_in0.40={n_rel40})", line_ok))
        mo2 = re.search(rf"\({re.escape(fam)}\) x_in=0\.40 UNRELAXED optical passes: (\d+) total", text)
        line2_ok = bool(mo2) and int(mo2.group(1)) == n_unrel40
        checks.append((f"results.md's {fam} x_in=0.40 UNRELAXED composition line matches an independent "
                        f"recount (n={n_unrel40})", line2_ok))
    checks.append(("results.md no longer claims horizontal-only x_in=0.40 unrelaxed passes are the ONLY ones "
                    "in EITHER family (the unqualified false claim this fix removed)",
                    "these are the ONLY x_in=0.40 optical passes in either family this run" not in text))


def check_qcse_field_label(text, core, checks):
    """H2 fix fixture: field_kVcm on an invalid row is the transport
    DEPLETION field alone (never the built-in piezoelectric/polarization
    field) -- the paragraph must label it as such and must not repeat the
    old false 'the built-in piezoelectric field this strain bound
    carries' framing attached to field_kVcm's own value; it must also
    report a separate levels-only F_pz/total-field replay."""
    parts = text.split("### QCSE excursion physics paragraph", 1)
    ok_present = len(parts) > 1
    body = parts[1] if ok_present else ""
    checks.append(("results.md's QCSE excursion physics paragraph exists", ok_present))
    checks.append(("QCSE paragraph labels field_kVcm as the depletion field, never the piezoelectric field",
                    ok_present and "depletion field" in body
                    and "the built-in piezoelectric field this strain bound carries" not in body))
    checks.append(("QCSE paragraph reports a levels-only F_pz/total-field replay",
                    ok_present and "F_pz_kVcm=" in body and "total_field_kVcm=" in body))
    mo = re.search(r"F_pz_kVcm=(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?), "
                    r"total_field_kVcm=(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)", body)
    checks.append(("the replayed F_pz_kVcm/total_field_kVcm values are finite numbers",
                    bool(mo) and isfin(mo.group(1)) and isfin(mo.group(2))))
    if mo:
        # numeric cross-check: replay fsim_core.nitride_nanowire_levels
        # directly (a DEPENDENCY module, never scripts/run_nitride_
        # nanowire.py itself -- the routine under test) at the paragraph's
        # own cited representative row's geometry.
        rid_mo = re.search(r"representative row \(deterministic rule[^)]*\), (\S+?):", body)
        row = next((r for r in core if r.get("row_id") == (rid_mo.group(1) if rid_mo else None)), None)
        if row is not None:
            try:
                sys.path.insert(0, str(ROOT))
                from fsim_core.nitride_nanowire_levels import NitrideNanowireSystem, levels
                depletion = f(row.get("field_kVcm"))
                sys_ = NitrideNanowireSystem(
                    height_nm=f(row.get("height_nm")), core_radius_nm=f(row.get("core_radius_nm")),
                    outer_radius_nm=f(row.get("outer_radius_nm")) or f(row.get("core_radius_nm")),
                    disc_radius_nm=None, x_in=f(row.get("x_in")), strain_bound="unrelaxed",
                    screening_fraction=f(row.get("screening_fraction")) or 0.0,
                    external_field_kVcm=depletion)
                lv = levels(sys_, T_K=f(row.get("T_j")) if isfin(row.get("T_j")) else f(row.get("T_hs")))
                replay_ok = (lv.valid and close_nan_safe(lv.F_pz_kVcm, f(mo.group(1)), rtol=1e-4)
                             and close_nan_safe(lv.field_kVcm, f(mo.group(2)), rtol=1e-4))
            except Exception:
                replay_ok = False
            checks.append((f"an independent levels-only replay of row {row.get('row_id')} reproduces the "
                            "printed F_pz_kVcm/total_field_kVcm numbers", replay_ok))


def check_figure_captions(out, man, checks):
    """M7 fix fixture: every exported figure carries the SAME one-line
    qualification footer (recorded into plot_row_mapping's own
    "__caption__" entry) and round-trips through matplotlib.image.imread
    without error -- "open the PNGs (matplotlib reads them back) and
    assert the footer text exists"."""
    mapping = man.get("plot_row_mapping", {})
    expected_figs = ("g2_flux_vs_core_radius.png", "g2_flux_vs_disc_thickness.png", "g2_flux_vs_ths.png",
                      "delivered_vs_commanded_flux.png", "hardware_screens.png", "strain_reversal_map.png")
    ok_caption = True
    ok_common = True
    ok_readable = True
    for figname in expected_figs:
        contract = mapping.get(figname) or {}
        caption = contract.get("__caption__")
        if not isinstance(caption, str) or len(caption) < 40:
            ok_caption = False
        if not isinstance(caption, str) or "Qualification (applies to every figure this run)" not in caption:
            ok_common = False
        p = out / figname
        if not p.is_file():
            ok_readable = False
            continue
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.image as mpimg
            img = mpimg.imread(str(p))
            if img is None or img.shape[0] < 10 or img.shape[1] < 10:
                ok_readable = False
        except Exception:
            ok_readable = False
    checks.append(("every exported figure's plot_row_mapping carries a non-empty '__caption__' footer string "
                    "(M7/L18: one-line qualification on every figure)", ok_caption))
    checks.append(("every exported figure's caption carries the SAME common qualification sentence "
                    "(M7: not distributed piecemeal)", ok_common))
    checks.append(("every exported figure PNG round-trips through matplotlib.image.imread without error",
                    ok_readable))


def check_plots(out, man, core, checks):
    mapping = man.get("plot_row_mapping", {})
    checks.append(("plot_row_mapping is non-empty", len(mapping) > 0))
    hashes = man.get("output_hashes", {})
    checks.append(("every recorded output_hash matches the file on disk",
                    len(hashes) > 0 and all((out / n).is_file()
                        and hashlib.sha256((out / n).read_bytes()).hexdigest() == v for n, v in hashes.items())))
    all_rows = {r["row_id"]: {**r, "T_hs": f(r.get("T_hs")), "core_radius_nm": f(r.get("core_radius_nm")),
                               "height_nm": f(r.get("height_nm")), "x_in": f(r.get("x_in")),
                               "rep_rate_hz": f(r.get("rep_rate_hz"))}
                for r in core}
    x_key_for_fig = {"g2_flux_vs_core_radius.png": "core_radius_nm", "g2_flux_vs_disc_thickness.png": "height_nm",
                      "g2_flux_vs_ths.png": "T_hs"}
    n_traced = 0
    n_points_checked = 0
    for figname, x_key in x_key_for_fig.items():
        contract = mapping.get(figname)
        if not contract:
            checks.append((f"plot-construction contract present: {figname}", False))
            continue
        xy_ok = True
        for key, entry in contract.items():
            if key == "__caption__":
                continue
            # L fix: `key` is the contract's OWN "family|y_key|label"
            # construction (see scripts/run_nitride_nanowire.py's
            # plot_vs_axis), so the plotted y column is read directly from
            # the key structure -- never guessed from the label text (the
            # prior `"g2_op" if "g2" in label.lower() or True else None`
            # was a tautology: `or True` makes the condition always True,
            # and the computed `ykey` was then never even used, so no
            # plotted y VALUE was ever checked against sweep.csv here).
            parts = key.split("|", 2)
            ykey = parts[1] if len(parts) >= 2 else None
            for rid, xv, yv in zip(entry["row_ids"], entry["x"], entry["y"]):
                row = all_rows.get(rid)
                if row is None:
                    xy_ok = False
                    continue
                if not close_nan_safe(row.get(x_key), xv):
                    xy_ok = False
                if ykey is not None and not close_nan_safe(f(row.get(ykey)), yv):
                    xy_ok = False
                n_points_checked += 1
            if not trace_shares_fixed_coords({"row_ids": entry["row_ids"]}, all_rows, x_key):
                xy_ok = False
        checks.append((f"plot-construction contract's traces all share fixed (non-swept) coordinates AND every "
                        f"plotted x/y value matches its row_id's sweep.csv columns: {figname} (a merged trace "
                        "or a corrupted plotted value fails this)", xy_ok))
        n_traced += 1
    checks.append(("at least one grouped figure's trace-fixed-coordinate contract was actually checked", n_traced > 0))
    checks.append((f"plotted x/y values checked against sweep.csv for {n_points_checked} points across traced "
                    "figures (L fix: the y value is now actually compared, not just computed and discarded)",
                    n_points_checked > 0))


RESULTS_TEXT_OBLIGATIONS = (
    "CANNOT deliver a 100 ps step",           # RC caveat
    "0.625 ns",                                # access-1.0 lifetime cap (X)
    "0.3125 ns",                                # access-1.0 lifetime cap (XX)
    "OPPOSITE endpoints",                       # opposite-endpoint anchor match
    "never averaged",
    # E_C/kT wall statement (M5 fix: the Coulomb wall is a demonstrated
    # result; the RT injector screen is explicitly NOT, so the text no
    # longer asserts "fails on both ... by EITHER charging mechanism").
    "E_C/kT wall", "fails the Coulomb-blockade screen", "unknown_incomplete",
    "wrong sign against the +70%",             # c-plane dipole prior falsification
    "one_pair_valid",                           # 200 vs 80 MHz comparison
    "quality_pass",                             # H3 gate anti-monotonicity fix
    "photons_per_cycle",                        # H3 gate anti-monotonicity fix
    "photons_per_cycle_delivered",              # M4 fix: delivered per-cycle beside commanded
    "depletion field",                          # H2 fix: field_kVcm on invalid rows is the depletion field
    "F_pz_kVcm",                                # H2 fix: levels-only polarization-field replay
    "DEMONSTRATES the anti-monotonicity",       # M3 fix: the pair that actually flips
    "does NOT repair this loss-induced ordering",  # M5 fix
    "share the SAME insufficient disc E_C/kT", # M6 fix
    "pulsed replay of a CW measurement",        # L fix: renamed heading
    "MEASURED leverage",                        # L fix: sensitivity ranking table
)


def check_results_text_obligations(text, checks):
    for phrase in RESULTS_TEXT_OBLIGATIONS:
        checks.append((f"results.md states the required obligation text: {phrase!r}", phrase in text))
    checks.append(("results.md carries the paired_optical_pass definition paragraph", "Definition: paired_optical_pass" in text))
    checks.append(("results.md carries a Deshpande 2013 comparison section (drive_mismatch, non-gating)",
                    "Deshpande 2013 comparison" in text and "drive_mismatch" in text))
    checks.append(("results.md carries a Deshpande 2014 comparison section", "Deshpande 2014 comparison" in text))
    checks.append(("results.md documents the reduced-cut budget trim and which axes were dropped",
                    "Reduced-cut coverage" in text and "DROPPED" in text))


def check_caps(man, quick, checks):
    checks.append(("evaluate_calls within the 10000 cap", man.get("evaluate_calls", 10**9) <= 10000))
    if quick:
        checks.append(("quick run stays within its own 600 s budget", man.get("runtime_s", 1e9) < 600.0))
        checks.append(("quick run's complete flag is False (never a full-coverage claim)", man.get("complete") is False))
    else:
        checks.append(("full run's complete flag is a bool", isinstance(man.get("complete"), bool)))
        checks.append(("full run reports its own runtime honestly (complete=True only if <1800s AND full "
                        "expected coverage; a budget overrun is reported as complete=False, not hidden)",
                        (man.get("complete") is False) or (man.get("runtime_s", 1e9) < 1800.0)))


# ------------------------------------------------------------------- main
def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--out-dir", default=None)
    a = ap.parse_args(argv)

    checks, detail = [], []
    check_predicate_boundaries(checks)
    check_rti_fabrication_fixture(checks)
    check_bound_partner_and_duplicate_fixtures(checks)
    check_headline_selection_fixture(checks)
    check_grid_self(checks, detail)
    check_reduced_cut_axes_declared(checks, detail)
    check_dry_run_cap(checks)
    check_evaluate_smoke(checks, detail)

    if a.out_dir:
        out = Path(a.out_dir)
        core = rows(out / "sweep.csv")
        man = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
        text = (out / "results.md").read_text(encoding="utf-8")
        quick = bool(man.get("quick"))

        check_files_present(out, checks)
        check_row_counts(core, man, quick, checks, detail)
        check_column_set(core, checks)
        check_card_hashes(core, man, checks)
        check_bound_pairing_and_duplicates(core, checks)
        check_verdict_lines(text, core, quick, checks, detail)
        check_best_passing_flux(text, core, checks)
        check_invalid_preserved(core, checks)
        check_rti_fabrication_real(core, checks)
        check_plots(out, man, core, checks)
        check_figure_captions(out, man, checks)
        check_composition_counts(text, core, checks)
        check_qcse_field_label(text, core, checks)
        check_results_text_obligations(text, checks)
        check_caps(man, quick, checks)

    passed = sum(1 for _, ok in checks if ok)
    total = len(checks)
    for name, ok in checks:
        if not ok:
            print("FAILED:", name)
    for line in detail:
        print("NOTE:", line)
    print("%d/%d nitride nanowire sweep checks passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    raise SystemExit(main())
