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


def check_csv_nan_bool_coercion_fixture(checks):
    """H1 fix, mutation-sensitivity fixture: a boolean column's underlying
    scalar can be a real Python None/float('nan') on rows where the
    predicate does not apply (e.g. rti_transport_feasible on rti_status=
    not_applicable rows); csv.DictWriter's str() rendering turns that into
    the literal CSV token "nan" (or "NaN"/"None"), which the OLD
    scripts.run_nitride_nanowire._coerce_csv_row's `v == ""` guard did NOT
    catch, so it fell through to `(v == "True")` and silently became False.
    Exercised entirely in memory -- a csv.DictReader row is already just a
    dict of strings, so no temp CSV file is needed to reproduce this."""
    checks.append(("_coerce_csv_row maps the CSV token 'nan' in a boolean column to None, never False (H1 fix)",
                    m._coerce_csv_row({"rti_transport_feasible": "nan"})["rti_transport_feasible"] is None))
    for token in ("NaN", "NAN", "None", "none", ""):
        checks.append((f"_coerce_csv_row maps the CSV token {token!r} in a boolean column to None",
                        m._coerce_csv_row({"set_feasible": token})["set_feasible"] is None))
    checks.append(("_coerce_csv_row still correctly parses a genuine 'False' token in a boolean column",
                    m._coerce_csv_row({"optical_pass": "False"})["optical_pass"] is False))
    checks.append(("_coerce_csv_row still correctly parses a genuine 'True' token in a boolean column",
                    m._coerce_csv_row({"valid": "True"})["valid"] is True))
    g2_coerced = m._coerce_csv_row({"g2_op": "nan"})["g2_op"]
    checks.append(("_coerce_csv_row leaves a NON-boolean numeric column's 'nan' token as a float NaN (only "
                    "boolean columns get the None treatment)", isinstance(g2_coerced, float) and math.isnan(g2_coerced)))

    # Regression count: a small synthetic row set, exactly as csv.DictReader
    # would hand it back (every value a string), mixing applicable (False)
    # and not_applicable (nan-token) rti_transport_feasible rows -- the OLD
    # bug counted ALL FOUR as False; the fix must count only the two
    # genuinely-applicable ones.
    synth = [
        {"rti_status": "unknown_incomplete", "rti_transport_feasible": "False"},
        {"rti_status": "unknown_incomplete", "rti_transport_feasible": "False"},
        {"rti_status": "not_applicable", "rti_transport_feasible": "nan"},
        {"rti_status": "not_applicable", "rti_transport_feasible": "nan"},
    ]
    coerced = [m._coerce_csv_row(r) for r in synth]
    n_false = sum(1 for r in coerced if r.get("rti_transport_feasible") is False)
    n_not_na = sum(1 for r in coerced if r.get("rti_status") != "not_applicable")
    checks.append(("synthetic fixture: rti_transport_feasible=False count EXCLUDES the not_applicable/'nan'-"
                    "token rows (2 of 4, not the old bug's 4 of 4)", n_false == 2))
    checks.append(("synthetic fixture: the False count equals the count of rows whose rti_status is NOT "
                    "not_applicable (the H1 denominator invariant)", n_false == n_not_na))


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
    "occupied_dot_access": {0.1, 1.0},  # audit C4 re-pin: 0.1 added (contract sensitivity set {0,0.1,1.0})
    "S_cm_s": {1.0e2, 1.0e4},
    "shell": {"AlGaN"},
    "al_fraction": {0.2},
    "growth_tolerance_steps": {2.0},
    "alignment_uncertainty_meV": {30.0},
    "occupation_control_uncertainty": {True},
    "injector_barrier_thickness_nm": {1.0},
    "R_s_ohm": {1.0e6},
    # audit C4 restored cuts (contract "Sweep grid": linewidth/tau_rad0/capture,
    # parasitic capacitance); C_parasitic_F values are the contract's own
    # Card-schema sensitivity set {1e-18, 1e-16}.
    "tau_rad0_ns": {0.2, 4.0},
    "tau_cap_ps": {1.0, 100.0},
    "C_parasitic_F": {1.0e-18, 1.0e-16},
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
    # audit C4 re-pin: the literal list below now names only the axes still
    # dropped (tau_rad0_ns / tau_cap_ps / C_parasitic_F were restored), and
    # every one of them must still be declared in m.DROPPED_CUT_AXES.
    _still_dropped = ("reservoir_access", "gamma300", "Rth_K_W", "NA", "bottom_reflectivity")
    checks.append(("every declared-dropped axis in DROPPED_CUT_AXES is actually absent from the full-mode cuts",
                    all(name not in by_axis for name in _still_dropped)
                    and all(any(name in d for d in m.DROPPED_CUT_AXES) for name in _still_dropped)))
    checks.append(("audit C4: the restored axes are no longer listed in DROPPED_CUT_AXES",
                    all(not any(name in d for d in m.DROPPED_CUT_AXES)
                        for name in ("tau_rad0_ns", "tau_cap_ps", "C_parasitic_F"))))
    # audit C4: the pre-C4 grid version reproduces the grid every older
    # sweep.csv was generated with (411 full-mode cut rows incl. the 2
    # dipole rows = 2992 - 2560 core - 4 D2013 - 16 planar - 1 replay; the
    # c4 grid adds 7 alternative values x 32 rows = 224).
    checks.append(("audit C4: reduced-cut grid versions differ by exactly the 224 restored-axis rows "
                    "(pre_c4 full grid = 409 + 2 dipole rows, the committed pre-C4 sweep's count)",
                    len(m.build_reduced_cuts(False, grid_version="pre_c4")) == 409
                    and len(cuts) - len(m.build_reduced_cuts(False, grid_version="pre_c4")) == 224))
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


def check_loading_window_axis(checks, detail):
    """Audit C6 (2026-09-23, M3; user decision Q6): the loading window is an
    OPT-IN sweep axis. The default grid must be exactly the pre-C6 grid (no
    row carries loading_window_ns, so every default design and physics-only
    identity is unchanged); the opt-in axis samples 0.1/0.3/1/3/10 ns
    strictly below each row's period (contract bullet 7)."""
    default_jobs = ([p for f in m.FAMILIES for p in m.build_core(f, False)] + m.build_reduced_cuts(False)
                    + m.build_dipole_falsification() + m.build_deshpande2013(False))
    checks.append(("audit C6: no default-grid row carries loading_window_ns (default grid unchanged)",
                    all("loading_window_ns" not in p for p in default_jobs)
                    and "loading_window_ns" not in m.full_defaults("horizontal_as_built")))
    d_def, _ = m._design(m.full_defaults("horizontal_as_built"))
    checks.append(("audit C6: a default design writes no drive.diode.loading_window_ns leaf (device defaults it to tau_pulse_ns)",
                    "loading_window_ns" not in d_def.drive.diode and d_def.drive.diode["tau_pulse_ns"] == 0.1))
    checks.append(("audit C6: planned counts without the flag add zero loading-window rows",
                    m.planned_counts(False)["loading_window_cut_rows"] == 0
                    and m.planned_counts(False) == m.planned_counts(False, loading_window_cut=False)))
    lw = m.build_loading_window_cut(False)
    # independent count: 2 families x 2 regimes x |CUT_T_HS| temperatures x
    # (5 windows below the 12.5 ns 80 MHz period + 4 below the 5 ns 200 MHz one)
    n_exp = 2 * 2 * len(m.CUT_T_HS) * (5 + 4)
    checks.append((f"audit C6: opt-in loading-window axis has {n_exp} full-mode rows, values {{0.1,0.3,1,3,10}} ns, "
                    "each strictly below its row's period",
                    len(lw) == n_exp
                    and {p["loading_window_ns"] for p in lw} == {0.1, 0.3, 1.0, 3.0, 10.0}
                    and all(p["loading_window_ns"] * 1e-9 < 1.0 / p["rep_rate_hz"] for p in lw)
                    and all(p["sensitivity_axis"] == "loading_window_ns" and p["tau_pulse_ns"] == 0.1 for p in lw)
                    and m.planned_counts(False, loading_window_cut=True)["loading_window_cut_rows"] == n_exp))
    p1 = dict(lw[0]); d1, _ = m._design(p1)
    checks.append(("audit C6: an opt-in row writes its window into drive.diode.loading_window_ns",
                    d1.drive.diode.get("loading_window_ns") == p1["loading_window_ns"]))
    checks.append(("audit C6: the quick axis is non-empty and a strict subset of the full axis values",
                    0 < len(m.build_loading_window_cut(True)) < len(lw)
                    and {p["loading_window_ns"] for p in m.build_loading_window_cut(True)} < {0.1, 0.3, 1.0, 3.0, 10.0}))
    detail.append("loading_window_cut rows(full)=%d (opt-in, --loading-window-cut)" % len(lw))


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

    # audit C4: compare against the grid version that PRODUCED this artifact
    # (manifest reduced_cut_grid_version; absent = generated before C4).
    # Review follow-up: default to pre_c4 ONLY for a genuinely pre-C4
    # artifact (no generation_commit either -- C4 full runs record both);
    # a C4-era manifest missing the grid-version key is an error.
    _gv = man.get("reduced_cut_grid_version") if isinstance(man, dict) else None
    _pre_c4_marker = isinstance(man, dict) and "generation_commit" not in man
    checks.append(("reduced_cut_grid_version present, or the manifest is a pre-C4 artifact (no generation_commit)",
                    _gv in m.REDUCED_CUT_GRID_VERSIONS or (_gv is None and _pre_c4_marker)))
    _grid_version = _gv if _gv in m.REDUCED_CUT_GRID_VERSIONS else ("pre_c4" if (_gv is None and _pre_c4_marker) else None)
    exp_cuts = (len(m.build_reduced_cuts(quick, grid_version=_grid_version)) + len(m.build_dipole_falsification())
                if _grid_version is not None else -1)
    # audit C6: the opt-in loading-window rows count only when the
    # manifest says the axis was run (absent = pre-C6 artifact = not run).
    _lw_run = bool(man.get("loading_window_cut", False)) if isinstance(man, dict) else False
    _lw_rows = [r for r in sens_rows if r.get("sensitivity_axis") == "loading_window_ns"]
    exp_lw = len(m.build_loading_window_cut(quick)) if _lw_run else 0
    checks.append(("sensitivity row count matches independent build_reduced_cuts()+build_dipole_falsification()",
                    len(sens_rows) == exp_cuts + len(m.build_deshpande2013(quick)) + exp_lw))
    checks.append(("audit C6: loading-window rows present iff the manifest records loading_window_cut=True",
                    len(_lw_rows) == exp_lw))
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


_SI_DIAG_RE = re.compile(r"si_complex_index: lambda_nm=([0-9.eE+-]+) outside")
_CODATA_EPS0 = 8.8541878128e-12  # F/m [V] CODATA 2018


def _reasons_of(row):
    try:
        reasons = json.loads(row.get("invalid_reasons", "[]"))
    except (TypeError, json.JSONDecodeError):
        return ["<undecodable invalid_reasons>"]
    return reasons if isinstance(reasons, list) else []


def independent_reason_label(reason):
    """Audit C3 item 6, this file's own transcription: a reason's label is
    its text up to the first ':' at parenthesis depth 0 (a ':' inside a
    parenthesised clause does not cut the label)."""
    depth = 0
    for i, ch in enumerate(reason):
        if ch == "(":
            depth += 1
        elif ch == ")":
            depth = max(0, depth - 1)
        elif ch == ":" and depth == 0:
            return reason[:i]
    return reason


def independent_invalid_triggers(core):
    """Audit C7: recompute, from sweep.csv rows alone (this file's own code,
    never the generator's _invalid_paragraph_triggers), the trigger the
    invalid-rows paragraph obeys. A row is invalid unless valid == True
    (a 'nan'/empty token counts as invalid, as the generator's falsy test
    does). Reason key = independent_reason_label (text up to the first ':'
    outside parentheses; audit C3 item 6). Returns dict with the
    per-key counts, the tied-max key set, invalid_total, dominant (> half of
    invalid rows), si (dominant and >= 1 invalid core row carrying the key
    AND a recoverable si_complex_index diagnostic wavelength), fc (dominant
    and key is field_collapse), plus the fc core rows."""
    counts, order = {}, []
    invalid_total = 0
    for r in core:
        if as_bool(r.get("valid")) is True:
            continue
        invalid_total += 1
        for reason in _reasons_of(r):
            key = independent_reason_label(reason)
            if key not in counts:
                order.append(key)
            counts[key] = counts.get(key, 0) + 1
    if not counts:
        return dict(counts={}, top_keys=set(), top_n=0, invalid_total=invalid_total,
                    dominant=False, si=False, fc=False, fc_rows=[])
    top_n = max(counts.values())
    top_keys = {k for k in order if counts[k] == top_n}
    top_key = next(k for k in order if counts[k] == top_n)  # dict max() keeps the first-inserted tie
    dominant = invalid_total > 0 and top_n / invalid_total > 0.5
    inv_core = [r for r in core if r.get("row_kind") == "core" and as_bool(r.get("valid")) is not True]
    carrying = [r for r in inv_core if any(x.startswith(top_key) for x in _reasons_of(r))]
    si = dominant and any(any(_SI_DIAG_RE.search(x) for x in _reasons_of(r)) for r in carrying)
    fc = dominant and top_key.startswith("field_collapse")
    return dict(counts=counts, top_keys=top_keys, top_key=top_key, top_n=top_n, invalid_total=invalid_total,
                dominant=dominant, si=si, fc=fc, fc_rows=carrying if fc else [])


def closed_form_polarization_field_kVcm(x_in, strain_bound, screening_fraction, external_kVcm, T_K):
    """Audit C7 closed form [DR]: the interface sheet-charge field inside an
    InGaN slab between GaN barriers, F = (1-s)*(P_sp,GaN - P_sp,InGaN -
    P_pz,InGaN)/(eps0*eps_r) + F_ext, with P_pz = P_total(strained) -
    P_sp (Bernardini, Fiorentini and Vanderbilt, PRB 56, R10024 (1997) [V]
    polarization constants as tabulated in fsim_core.nitride_materials;
    eps0 CODATA 2018 [V]). 1 V/m = 1e-5 kV/cm. Returns (F_sp, F_pz, F) in
    kV/cm and the strained gap Ec-Ev (eV)."""
    from fsim_core.nitride_materials import band_edges, binary, ingaN
    d, g = ingaN(x_in), binary("GaN")
    sf = 0.0 if strain_bound == "relaxed" else 1.0
    de = band_edges(d, T_K, substrate=g, strain_fraction=sf)
    f_sp = (g.Psp_Cm2 - d.Psp_Cm2) / (_CODATA_EPS0 * d.eps_r) * 1e-5
    f_pz = -(de["P_total_Cm2"] - d.Psp_Cm2) / (_CODATA_EPS0 * d.eps_r) * 1e-5
    F = (1.0 - screening_fraction) * (f_sp + f_pz) + external_kVcm
    return f_sp, f_pz, F, de["Ec_eV"] - de["Ev_eV"]


def check_invalid_paragraph_trigger_fixture(checks):
    """Audit C7 fixture (no artifacts): synthetic rows fed to the
    GENERATOR's trigger/paragraph helpers must yield the same flags as this
    file's independent recomputation, for a field_collapse-dominated, an
    Si-index-dominated and a no-dominant case."""
    zr = "field_collapse (|F|*h_eff >= strained InGaN gap: interband Zener breakdown, unscreened field not self-consistent)"
    si = "si_complex_index: lambda_nm=812.5 outside 380-750 nm"
    pc = "photon counting did not converge: x"
    def mk(i, reason, valid=False):
        return {"row_id": f"CO{i:05d}", "row_kind": "core", "valid": "True" if valid else "False",
                "invalid_reasons": json.dumps([] if valid else [reason])}
    cases = {"field_collapse-dominated": [mk(i, zr) for i in range(3)] + [mk(9, pc), mk(10, None, True)],
             "si-dominated": [mk(i, si) for i in range(3)] + [mk(9, pc)],
             "no-dominant": [mk(0, zr), mk(1, pc)]}
    ok = True
    for tag, rs in cases.items():
        ind = independent_invalid_triggers(rs)
        coerced = [m._coerce_csv_row(r) for r in rs]
        inv_total = sum(1 for r in coerced if not r.get("valid"))
        dom, si_rows, fc = m._invalid_paragraph_triggers(coerced, ind["top_key"], ind["top_n"], inv_total)
        para = "\n".join(m._dominant_invalid_paragraph(ind["top_key"], ind["top_n"], inv_total, dom, bool(si_rows), fc))
        ok = ok and (dom, bool(si_rows), fc) == (ind["dominant"], ind["si"], ind["fc"]) and \
            f"si_index_excursion_trigger={ind['si']}" in para and f"field_collapse_dominant={ind['fc']}" in para
    want = {"field_collapse-dominated": (True, False, True), "si-dominated": (True, True, False),
            "no-dominant": (False, False, False)}
    ok = ok and all((independent_invalid_triggers(cases[k])["dominant"], independent_invalid_triggers(cases[k])["si"],
                     independent_invalid_triggers(cases[k])["fc"]) == v for k, v in want.items())
    checks.append(("C7 fixture: generator's invalid-paragraph trigger flags equal the independent recomputation "
                    "(field_collapse-dominated / Si-dominated / no-dominant synthetic rows)", ok))


def check_invalid_paragraph(text, core, checks):
    """Audit C7: the invalid-rows paragraph is ALWAYS present and states the
    dominant reason and count; its trigger flags, and the presence/absence of
    the Si-index excursion and field-collapse sub-paragraphs, match the
    trigger recomputed here from sweep.csv. When field_collapse dominates,
    the printed field and drop are checked against a closed-form sheet-charge
    field (not the generator's number) and the Zener identity drop >= gap."""
    ind = independent_invalid_triggers(core)
    head = "### QCSE / invalid-rows paragraph (dominant invalid reason)"
    checks.append(("C7: results.md's QCSE / invalid-rows paragraph is present (unconditional)", head in text))
    body = text.split(head, 1)[1] if head in text else ""
    if ind["counts"]:
        mo = re.search(r"Dominant invalid reason: '(.*?)' on (\d+)/(\d+) invalid rows \(\d+ percent\); "
                        r"dominant=(True|False); si_index_excursion_trigger=(True|False); "
                        r"field_collapse_dominant=(True|False)\.", body)
        ok = (bool(mo) and mo.group(1) in ind["top_keys"] and int(mo.group(2)) == ind["top_n"]
              and int(mo.group(3)) == ind["invalid_total"] and mo.group(4) == str(ind["dominant"])
              and mo.group(5) == str(ind["si"]) and mo.group(6) == str(ind["fc"]))
    else:
        ok = "Dominant invalid reason: none (0 invalid rows this run)" in body
    checks.append((f"C7: the dominant invalid reason, its count and the trigger flags match an independent "
                    f"recount from sweep.csv (top_n={ind['top_n']}/{ind['invalid_total']}, dominant="
                    f"{ind['dominant']}, si={ind['si']}, fc={ind['fc']})", ok))
    si_present = "### QCSE excursion physics paragraph" in text
    fc_head = "### Field-collapse (QCSE collapse) paragraph"
    fc_present = fc_head in text
    checks.append(("C7: the Si-index QCSE excursion paragraph is present iff the Si-index trigger holds",
                    si_present == ind["si"]))
    checks.append(("C7: the field-collapse paragraph is present iff field_collapse is the dominant reason",
                    fc_present == ind["fc"]))
    qcse_needed = ind["si"] or ind["fc"]
    for phrase in ("depletion field", "F_pz_kVcm"):
        checks.append((f"results.md states the QCSE obligation text {phrase!r} whenever the Si-index or "
                        f"field_collapse trigger holds (trigger={qcse_needed})",
                        (phrase in text) if qcse_needed else True))
    if not ind["fc"]:
        return
    fbody = text.split(fc_head, 1)[1].split("\n## ", 1)[0] if fc_present else ""
    fc_rows = ind["fc_rows"]
    mo_n = re.search(r"The (\d+) invalid core rows carrying", fbody)
    checks.append((f"C7: field-collapse paragraph's row count matches an independent recount (n={len(fc_rows)})",
                    bool(mo_n) and int(mo_n.group(1)) == len(fc_rows)))
    sample = min(fc_rows, key=lambda r: r["row_id"]) if fc_rows else None
    rid = re.search(r"representative row \(deterministic rule[^)]*\), (\S+?):", fbody)
    checks.append(("C7: field-collapse representative row is the smallest row_id among the rows",
                    sample is not None and bool(rid) and rid.group(1) == sample["row_id"]))
    checks.append(("C7: field-collapse paragraph labels field_kVcm as the depletion field, never the polarization "
                    "field", "depletion field alone, never the built-in polarization field" in fbody))
    num = r"(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)"
    mo = re.search(rf"F_sp_kVcm={num}, F_pz_kVcm={num}, total_field_kVcm={num}; field drop \|F\|\*h={num} eV "
                    rf"across h={num} nm against the strained InGaN gap Ec-Ev={num} eV \(drop >= gap: (True|False)\)",
                    fbody)
    ok_cf = False
    if mo and sample is not None:
        try:
            dep = f(sample.get("field_kVcm"))
            ext = dep if math.isfinite(dep) else 0.0
            T = f(sample.get("T_j")) if isfin(sample.get("T_j")) else f(sample.get("T_hs"))
            h = f(sample.get("height_nm"))
            f_sp, f_pz, F, gap = closed_form_polarization_field_kVcm(
                f(sample.get("x_in")), sample.get("strain_bound"), f(sample.get("screening_fraction")) or 0.0, ext, T)
            drop = abs(F) * 1e-4 * h  # e*(kV/cm)*nm = 1e-4 eV [DR]
            ok_cf = (close_nan_safe(f(mo.group(1)), f_sp, rtol=1e-4) and close_nan_safe(f(mo.group(2)), f_pz, rtol=1e-4)
                     and close_nan_safe(f(mo.group(3)), F, rtol=1e-4) and close_nan_safe(f(mo.group(4)), drop, rtol=1e-4)
                     and close_nan_safe(f(mo.group(5)), h) and close_nan_safe(f(mo.group(6)), gap, rtol=1e-4)
                     and drop >= gap and mo.group(7) == "True")
        except Exception:
            ok_cf = False
    checks.append(("C7 expectation: the printed F_sp/F_pz/total field and drop |F|*h equal the closed-form "
                    "Bernardini sheet-charge field (1-s)(P_sp,GaN-P_sp,InGaN-P_pz)/(eps0 eps_r)+F_ext and "
                    "the Zener identity |F|*h >= strained gap Ec-Ev holds for the representative row", ok_cf))
    ok_lv = False
    if sample is not None:
        try:
            from fsim_core.nitride_nanowire_levels import NitrideNanowireSystem, levels
            dep = f(sample.get("field_kVcm"))
            sys_ = NitrideNanowireSystem(
                height_nm=f(sample.get("height_nm")), core_radius_nm=f(sample.get("core_radius_nm")),
                outer_radius_nm=f(sample.get("outer_radius_nm")) or f(sample.get("core_radius_nm")),
                disc_radius_nm=None, x_in=f(sample.get("x_in")), strain_bound=sample.get("strain_bound"),
                screening_fraction=f(sample.get("screening_fraction")) or 0.0,
                external_field_kVcm=dep if math.isfinite(dep) else 0.0)
            lv = levels(sys_, T_K=f(sample.get("T_j")) if isfin(sample.get("T_j")) else f(sample.get("T_hs")))
            ok_lv = (not lv.valid) and any(str(x).startswith("field_collapse") for x in lv.invalid_reasons)
        except Exception:
            ok_lv = False
    checks.append(("C7: an independent levels-only replay of the representative row is invalid with the "
                    "field_collapse reason", ok_lv))


def check_qcse_field_label(text, core, checks):
    """H2 fix fixture: field_kVcm on an invalid row is the transport
    DEPLETION field alone (never the built-in piezoelectric/polarization
    field) -- the paragraph must label it as such and must not repeat the
    old false 'the built-in piezoelectric field this strain bound
    carries' framing attached to field_kVcm's own value; it must also
    report a separate levels-only F_pz/total-field replay. Audit C7: these
    Si-index-excursion checks apply only when the Si-index trigger
    (independent_invalid_triggers) holds, the same condition under which the
    generator writes the paragraph; otherwise check_invalid_paragraph
    covers the section."""
    check_invalid_paragraph(text, core, checks)
    if not independent_invalid_triggers(core)["si"]:
        return
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


def check_rti_transport_feasible_denominator(text, core, checks):
    """H1 fix, artifact-mode regression check: independently recompute, from
    sweep.csv's own deterministic_pair core rows via THIS file's own
    as_bool() (never scripts/run_nitride_nanowire.py's _coerce_csv_row),
    the rti_transport_feasible=False count and cross-check it against
    results.md's printed 'rti_transport_feasible=False on N/M' line -- a
    regression guard against the fixed bug where the 96 rti_status=
    not_applicable rows (a 'nan' CSV token in that boolean column) were
    silently counted as False, inflating N to M. The printed denominator M
    equals the total deterministic_pair core row count; the correct
    numerator N equals the number of those rows whose rti_status is NOT
    not_applicable (contract: the artifact-mode check the H1 fix requires)."""
    set_rows = [r for r in core if r.get("row_kind") == "core" and r.get("regime") == "deterministic_pair"]
    n_total = len(set_rows)
    n_not_na = sum(1 for r in set_rows if r.get("rti_status") != "not_applicable")
    n_false = sum(1 for r in set_rows if as_bool(r.get("rti_transport_feasible")) is False)
    checks.append(("rti_transport_feasible=False count equals the number of deterministic_pair core rows whose "
                    "rti_status is NOT not_applicable (H1 fix denominator invariant: the not_applicable/'nan' "
                    f"rows must never be counted as False; n_false={n_false} n_not_applicable_excluded="
                    f"{n_total - n_not_na})", n_false == n_not_na))
    mo = re.search(r"rti_transport_feasible=False on (\d+)/(\d+)", text)
    checks.append(("results.md's printed rti_transport_feasible=False count/total matches an independent "
                    f"recount from sweep.csv (expected {n_false}/{n_total})",
                    bool(mo) and int(mo.group(1)) == n_false and int(mo.group(2)) == n_total))


def _footer_bbox_ok(bbox):
    """L4 fix: a recorded footer bbox (scripts/run_nitride_nanowire.py's
    _footer_bbox_fraction, stored per figure at plot_row_mapping[fig]
    ["__footer_bbox__"]) is a fraction of the FINAL saved PNG -- it is only
    valid, i.e. the footer is NOT clipped, if fully inside that image: 0 <=
    x0 < x1 <= 1 and 0 <= y0 < y1 <= 1. This is a real geometric assertion
    on where matplotlib actually rendered the footer, never just an
    imread() round-trip (which succeeds on a badly clipped image too)."""
    if not isinstance(bbox, dict):
        return False
    try:
        x0, x1 = float(bbox["x0_frac"]), float(bbox["x1_frac"])
        y0, y1 = float(bbox["y0_frac"]), float(bbox["y1_frac"])
    except (KeyError, TypeError, ValueError):
        return False
    if not all(math.isfinite(v) for v in (x0, x1, y0, y1)):
        return False
    return 0.0 <= x0 < x1 <= 1.0 and 0.0 <= y0 < y1 <= 1.0


def check_footer_bbox_fixture(checks):
    """L4 fix, mutation-sensitivity fixture (no artifacts needed): a footer
    bbox whose right/left/bottom edge falls outside the saved PNG (a
    clipped footer) must be REJECTED by _footer_bbox_ok, not silently
    accepted the way a bare imread() success would accept it."""
    good = {"x0_frac": 0.08, "x1_frac": 0.92, "y0_frac": 0.95, "y1_frac": 0.995}
    checks.append(("_footer_bbox_ok accepts a footer bbox fully inside the saved PNG", _footer_bbox_ok(good) is True))
    checks.append(("_footer_bbox_ok rejects a footer whose RIGHT edge falls outside the saved PNG (clipped)",
                    _footer_bbox_ok(dict(good, x1_frac=1.15)) is False))
    checks.append(("_footer_bbox_ok rejects a footer whose LEFT edge falls outside the saved PNG (clipped)",
                    _footer_bbox_ok(dict(good, x0_frac=-0.05)) is False))
    checks.append(("_footer_bbox_ok rejects a footer extending past the BOTTOM of the saved PNG (clipped)",
                    _footer_bbox_ok(dict(good, y1_frac=1.02)) is False))
    checks.append(("_footer_bbox_ok rejects a missing/malformed bbox",
                    _footer_bbox_ok(None) is False and _footer_bbox_ok({}) is False))


def check_figure_captions(out, man, checks):
    """M7 fix fixture: every exported figure carries the SAME one-line
    qualification footer (recorded into plot_row_mapping's own
    "__caption__" entry) and round-trips through matplotlib.image.imread
    without error -- "open the PNGs (matplotlib reads them back) and
    assert the footer text exists". L4 fix: additionally assert the
    footer's own recorded bbox (__footer_bbox__) is fully inside the saved
    PNG -- i.e. the footer's rightmost/bottommost pixel column/row is
    inside the PNG's actual pixel dimensions -- never just that the file
    decodes."""
    mapping = man.get("plot_row_mapping", {})
    expected_figs = ("g2_flux_vs_core_radius.png", "g2_flux_vs_disc_thickness.png", "g2_flux_vs_ths.png",
                      "delivered_vs_commanded_flux.png", "hardware_screens.png", "strain_reversal_map.png")
    ok_caption = True
    ok_common = True
    ok_readable = True
    ok_bbox = True
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
            ok_bbox = False
            continue
        try:
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.image as mpimg
            img = mpimg.imread(str(p))
            if img is None or img.shape[0] < 10 or img.shape[1] < 10:
                ok_readable = False
                png_w = png_h = None
            else:
                png_h, png_w = img.shape[0], img.shape[1]
        except Exception:
            ok_readable = False
            png_w = png_h = None
        bbox = contract.get("__footer_bbox__")
        if not _footer_bbox_ok(bbox):
            ok_bbox = False
        elif png_w is not None:
            # Literal L4 ask: "the footer's rightmost pixel column is
            # inside the PNG width from a stored bbox" -- convert the
            # recorded fraction to actual pixel columns/rows of THIS run's
            # own saved PNG and check both edges land inside it.
            x1_px = bbox["x1_frac"] * png_w
            y1_px = bbox["y1_frac"] * png_h
            if not (0.0 <= x1_px <= png_w and 0.0 <= y1_px <= png_h):
                ok_bbox = False
    checks.append(("every exported figure's plot_row_mapping carries a non-empty '__caption__' footer string "
                    "(M7/L18: one-line qualification on every figure)", ok_caption))
    checks.append(("every exported figure's caption carries the SAME common qualification sentence "
                    "(M7: not distributed piecemeal)", ok_common))
    checks.append(("every exported figure PNG round-trips through matplotlib.image.imread without error",
                    ok_readable))
    checks.append(("every exported figure's recorded footer bbox is fully inside its saved PNG's actual pixel "
                    "dimensions (L4 fix: a real clipped-footer detector, not just an imread round-trip)", ok_bbox))


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
            if key in ("__caption__", "__footer_bbox__"):
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

    # L3 fix: the RC figure (plot_delivered_vs_commanded) and the hardware
    # figure (plot_hardware_screens) use their OWN contract key shapes,
    # neither of which is the "family|y_key|label" convention the loop
    # above parses -- they were previously entirely unchecked. Extend the
    # y-value comparison to both so every mapped point across every figure
    # is compared, not just the three plot_vs_axis figures.
    rc_contract = mapping.get("delivered_vs_commanded_flux.png") or {}
    rc_ok = True
    n_rc_checked = 0
    for key, entry in rc_contract.items():
        if key in ("__caption__", "__footer_bbox__"):
            continue
        for rid, xv, yv in zip(entry["row_ids"], entry["x"], entry["y"]):
            row = all_rows.get(rid)
            if row is None:
                rc_ok = False
                continue
            if not close_nan_safe(f(row.get("collected_flux_pulsed_s")), xv):
                rc_ok = False
            if not close_nan_safe(f(row.get("collected_flux_delivered_s")), yv):
                rc_ok = False
            n_rc_checked += 1
    checks.append((f"delivered_vs_commanded_flux.png: every plotted (commanded, delivered) point matches its "
                    f"row_id's own sweep.csv columns ({n_rc_checked} points checked)",
                    rc_ok and n_rc_checked > 0))
    n_points_checked += n_rc_checked

    hw_contract = mapping.get("hardware_screens.png") or {}
    hw_ok = True
    n_hw_checked = 0
    for key, entry in hw_contract.items():
        if key in ("__caption__", "__footer_bbox__"):
            continue
        # hardware_screens' own contract key is "y_key|family" (the
        # y_key/family order is REVERSED from plot_vs_axis's "family|
        # y_key|label" convention -- see plot_hardware_screens).
        ykey = key.split("|", 1)[0]
        for rid, xv, yv in zip(entry["row_ids"], entry["x"], entry["y"]):
            row = all_rows.get(rid)
            if row is None:
                hw_ok = False
                continue
            if not close_nan_safe(row.get("core_radius_nm"), xv):
                hw_ok = False
            if not close_nan_safe(f(row.get(ykey)), yv):
                hw_ok = False
            n_hw_checked += 1
    checks.append((f"hardware_screens.png: every plotted (core_radius_nm, {{set_EC_over_kT,rti_bypass_fraction,"
                    f"rti_alignment_error_e_meV}}) point matches its row_id's own sweep.csv columns "
                    f"({n_hw_checked} points checked)", hw_ok and n_hw_checked > 0))
    n_points_checked += n_hw_checked

    checks.append((f"plotted x/y values checked against sweep.csv for {n_points_checked} points across traced "
                    "figures (L fix: the y value is now actually compared, not just computed and discarded; L3 "
                    "fix: the RC and hardware figures are now included, not just the three plot_vs_axis figures)",
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
    # "depletion field" / "F_pz_kVcm" (H2 fix) moved to check_qcse_field_label
    # (audit C7): required whenever the Si-index or field_collapse trigger
    # holds, i.e. under the same condition the generator writes them.
    "DEMONSTRATES the anti-monotonicity",       # M3 fix: the pair that actually flips
    "does NOT repair this loss-induced ordering",  # M5 fix
    "share the SAME insufficient disc E_C/kT", # M6 fix
    "pulsed replay of a CW measurement",        # L fix: renamed heading
    "MEASURED leverage",                        # L fix: sensitivity ranking table
)


def check_dipole_weights_ranking(text, core, checks):
    """M2 fix: build_dipole_falsification() evaluates DW02970/DW02971 ONLY
    in the rectangular regime (it never overrides full_defaults'
    regime="rectangular" default); the sensitivity-ranking table's
    per-axis filter used to hardcode regime="deterministic_pair" for every
    axis, so dipole_weights matched zero rows there and was printed with a
    leverage of 0. Independently recompute its commanded-flux leverage
    from sweep.csv (against the SAME family's rectangular-regime reference
    row at the reference geometry) and cross-check it is printed, nonzero,
    in the horizontal_as_built commanded-flux ranking line."""
    dw = [r for r in core if r.get("sensitivity_axis") == "dipole_weights"]
    checks.append(("dipole_weights rows are present to check the M2 ranking fix on", len(dw) > 0))
    if not dw:
        return
    ref_row = next((r for r in core if r.get("row_kind") == "core" and r.get("family") == "horizontal_as_built"
                     and r.get("regime") == "rectangular" and r.get("strain_bound") == "relaxed"
                     and close_nan_safe(r.get("T_hs"), 300.0) and close_nan_safe(r.get("rep_rate_hz"), 200.0e6)
                     and close_nan_safe(r.get("core_radius_nm"), f(dw[0].get("core_radius_nm")))
                     and close_nan_safe(r.get("height_nm"), f(dw[0].get("height_nm")))
                     and close_nan_safe(r.get("x_in"), f(dw[0].get("x_in")))), None)
    checks.append(("dipole_weights' own rectangular-regime reference row is present in this run's core coverage",
                    ref_row is not None))
    if ref_row is None or not isfin(ref_row.get("collected_flux_pulsed_s")):
        return
    ref_flux = f(ref_row.get("collected_flux_pulsed_s"))
    fluxes = [f(r.get("collected_flux_pulsed_s")) for r in dw if isfin(r.get("collected_flux_pulsed_s"))]
    lev = max(abs(v / ref_flux - 1.0) for v in fluxes) if fluxes else 0.0
    sec = text.split("### horizontal_as_built", 1)
    body = sec[1].split("### vertical_photonic", 1)[0] if len(sec) > 1 else ""
    mo = re.search(r"Ranked by commanded flux \|ratio-1\| \(measured leverage\): ([^\n]+)", body)
    lev_pat = re.escape(f"{lev:.3g}")
    checks.append((f"results.md's horizontal_as_built commanded-flux ranking line lists dipole_weights with a "
                    f"NONZERO leverage matching an independent recount from sweep.csv (expected {lev:.3g}, "
                    "not the old bug's 0)",
                    bool(mo) and re.search(rf"dipole_weights\({lev_pat}\)", mo.group(1)) is not None
                    and lev > 0.0))
    checks.append(("results.md states the per-axis regime used for the sensitivity ranking (M2 fix)",
                    "Regime used per axis for this ranking" in body and "dipole_weights=rectangular" in body))


def check_results_text_obligations(text, checks):
    for phrase in RESULTS_TEXT_OBLIGATIONS:
        checks.append((f"results.md states the required obligation text: {phrase!r}", phrase in text))
    checks.append(("results.md carries the paired_optical_pass definition paragraph", "Definition: paired_optical_pass" in text))
    checks.append(("results.md carries a Deshpande 2013 comparison section (drive_mismatch, non-gating)",
                    "Deshpande 2013 comparison" in text and "drive_mismatch" in text))
    checks.append(("results.md carries a Deshpande 2014 comparison section", "Deshpande 2014 comparison" in text))
    checks.append(("results.md documents the reduced-cut budget trim and which axes were dropped",
                    "Reduced-cut coverage" in text and "DROPPED" in text))


def check_c4_prose(text, core, checks):
    """Audit C4 pins, recomputed from sweep.csv here (never read back from
    the generator): (1) the Deshpande per-anchor photon-energy residuals are
    printed as computed values hc/lambda_pred - hc/lambda_meas (hc =
    1239.841984 eV nm, CODATA 2018; measured 436.56 nm [Deshpande et al.,
    Nat. Commun. 2013 Fig. 3c] and ~630 nm [Deshpande et al., APL 2014
    abstract]), with no hard-coded 'roughly 0.3 eV'; (2) the 80-vs-200 MHz
    'period effect' wording is stated per family from the rows' own
    one_pair_valid flags."""
    import re as _re
    hc = 1239.841984
    checks.append(("C4: no hard-coded 'roughly 0.3 eV' anchor residual in results.md", "roughly 0.3 eV" not in text))
    d13 = [r for r in core if r.get("sensitivity_axis") == "deshpande2013_comparison" and f(r.get("core_radius_nm")) == 12.5]
    ref = m.REF_GEOM["horizontal_as_built"]
    def _d14(sb):
        c = [r for r in core if r.get("row_kind") == "core" and r.get("family") == "horizontal_as_built"
             and r.get("regime") == "deterministic_pair" and r.get("strain_bound") == sb and f(r.get("T_hs")) == 300.0
             and f(r.get("rep_rate_hz")) == 200.0e6 and f(r.get("x_in")) == 0.40
             and f(r.get("core_radius_nm")) == ref["core_radius_nm"] and f(r.get("height_nm")) == ref["height_nm"]]
        return f(c[0].get("lambda_nm")) if c else float("nan")
    expect = [("2013 relaxed", f(d13[0].get("lambda_nm")) if d13 else float("nan"), 436.56),
              ("2014 relaxed", _d14("relaxed"), 630.0), ("2014 unrelaxed", _d14("unrelaxed"), 630.0)]
    ok = True
    for tag, lam, meas in expect:
        mm = _re.search(_re.escape(tag) + r" ([+-][0-9.]+) eV \(predicted", text)
        ok = ok and lam == lam and mm is not None and abs(float(mm.group(1)) - (hc / lam - hc / meas)) <= 6e-4
    checks.append(("C4: per-anchor residuals in results.md equal hc/lambda_pred - hc/lambda_meas recomputed from sweep.csv", ok))
    ok = True
    for fam in m.FAMILIES:
        g = m.REF_GEOM[fam]
        for t in (230.0, 300.0):
            pair = {}
            for r in core:
                if (r.get("row_kind") == "core" and r.get("family") == fam and r.get("regime") == "deterministic_pair"
                        and r.get("strain_bound") == "relaxed" and f(r.get("T_hs")) == t and f(r.get("x_in")) == 0.40
                        and f(r.get("core_radius_nm")) == g["core_radius_nm"] and f(r.get("height_nm")) == g["height_nm"]):
                    pair.setdefault(f(r.get("rep_rate_hz")), r)
            line = next((l for l in text.splitlines() if l.startswith(f"{fam} T_hs={t:g}K: at 80 MHz")), "")
            if 80.0e6 not in pair or 200.0e6 not in pair:
                ok = ok and "not found" in line
                continue
            v80 = pair[80.0e6].get("one_pair_valid") == "True"; v200 = pair[200.0e6].get("one_pair_valid") == "True"
            want = ("PERIOD EFFECT" if (v80 and not v200) else "NOT a period effect" if (not v80 and not v200)
                    else "no period effect" if (v80 and v200) else "inverted")
            ok = ok and want in line and (want == "PERIOD EFFECT" or "-> PERIOD EFFECT" not in line)
            if want == "NOT a period effect":
                both = f(pair[80.0e6].get("blocked_load_probability")) > 1e-9 and f(pair[200.0e6].get("blocked_load_probability")) > 1e-9
                ok = ok and (("both above the 1e-9 threshold" in line) == both)
    checks.append(("C4: the 80/200 MHz period-effect clause per family/T_hs matches the rows' own one_pair_valid flags", ok))


# ------------------------------------------------ audit C3 prose checks
# Each check below recomputes its expectation from sweep.csv / manifest.json
# with this file's own code (never the generator's helpers).
C3_NON_CUT_AXES = ("", "deshpande2013_comparison")


def _c3_norm(tok):
    tok = str(tok).strip()
    try:
        return ("n", round(float(tok), 12))
    except ValueError:
        return ("s", tok)


def _c3_norm_set(cell, sep=","):
    return {_c3_norm(t) for t in str(cell).split(sep) if t.strip()}


def independent_cut_coverage(core):
    """axis -> coverage sets, from raw sweep.csv strings."""
    cov = {}
    for r in core:
        if r.get("row_kind") != "sensitivity" or (r.get("sensitivity_axis") or "") in C3_NON_CUT_AXES:
            continue
        c = cov.setdefault(r["sensitivity_axis"], dict(values=set(), rows=0, families=set(), regimes=set(),
                                                         strain_bounds=set(), T_hs_K=set(), rep_rate_hz=set()))
        c["values"].add(_c3_norm(r.get("sensitivity_value")))
        c["rows"] += 1
        c["families"].add(_c3_norm(r.get("family")))
        c["regimes"].add(_c3_norm(r.get("regime")))
        c["strain_bounds"].add(_c3_norm(r.get("strain_bound")))
        c["T_hs_K"].add(_c3_norm(r.get("T_hs")))
        c["rep_rate_hz"].add(_c3_norm(r.get("rep_rate_hz")))
    all_fams = {_c3_norm(x) for x in m.FAMILIES}
    for c in cov.values():
        c["full"] = (c["families"] == all_fams and c["regimes"] == {_c3_norm(x) for x in ("rectangular", "deterministic_pair")}
                     and c["strain_bounds"] == {_c3_norm(x) for x in ("relaxed", "unrelaxed")}
                     and c["T_hs_K"] == {_c3_norm(230.0), _c3_norm(300.0)}
                     and c["rep_rate_hz"] == {_c3_norm(80.0e6), _c3_norm(200.0e6)})
    return cov


def _section(text, heading):
    parts = text.split(heading, 1)
    return parts[1].split("\n## ", 1)[0] if len(parts) > 1 else ""


def check_c3_reduced_cut_list(text, core, checks):
    """Item 1: the 'Reduced-cut coverage' section's cut list is the set of
    cuts actually in sweep.csv (axis, values, row count, coverage), not a
    literal that omits the restored C4 cuts."""
    cov = independent_cut_coverage(core)
    body = _section(text, "## Reduced-cut coverage")
    table = {}
    for line in body.splitlines():
        cells = [c.strip() for c in line.strip().strip("|").split("|")]
        if len(cells) == 8 and cells[0] in cov:
            table[cells[0]] = cells
    ok = set(table) == set(cov) and len(cov) > 0
    for ax, c in cov.items():
        cells = table.get(ax)
        ok = ok and cells is not None and _c3_norm_set(cells[1], ", ") == c["values"] and cells[2] == str(c["rows"]) \
            and _c3_norm_set(cells[3], ", ") == c["families"] and _c3_norm_set(cells[4], ", ") == c["regimes"] \
            and _c3_norm_set(cells[5], ", ") == c["strain_bounds"] and _c3_norm_set(cells[6], ", ") == c["T_hs_K"] \
            and _c3_norm_set(cells[7], ", ") == c["rep_rate_hz"]
    ok = ok and "occupied_dot_access's 1.0 conservative partner" not in body
    detail = ", ".join(f"{ax}:{c['rows']}" for ax, c in sorted(cov.items()))
    checks.append(("C3 item 1: results.md's Reduced-cut coverage table lists exactly the cuts in sweep.csv with "
                    f"their values, row counts and coverage (independent recount: {detail})", ok))


def check_c3_ec_wall_denominator(text, core, checks):
    """Item 2: the E_C/kT wall's set_feasible / rti_feasible denominators
    are labelled as VALID rows and the bin's total/invalid counts are
    printed beside them."""
    wall = [r for r in core if r.get("row_kind") == "core" and r.get("regime") == "deterministic_pair"
            and isfin(r.get("core_radius_nm")) and f(r.get("core_radius_nm")) >= 10.0]
    valid_ec = [r for r in wall if as_bool(r.get("valid")) is True and isfin(r.get("set_EC_over_kT"))]
    valid_rti = [r for r in wall if as_bool(r.get("valid")) is True and isfin(r.get("rti_level_margin_kT"))]
    n_set = sum(1 for r in valid_ec if as_bool(r.get("set_feasible")) is True)
    n_rti = sum(1 for r in valid_rti if as_bool(r.get("rti_feasible")) is True)
    n_inv = sum(1 for r in wall if as_bool(r.get("valid")) is not True)
    mo = re.search(r"set_feasible=True count=(\d+)/(\d+) valid rows; rti_feasible=True count=(\d+)/(\d+) valid "
                    r"rows \(of (\d+) deterministic_pair core rows with core_radius_nm>=10 nm in total; the (\d+) "
                    r"invalid rows", text)
    ok = bool(mo) and [int(g) for g in mo.groups()] == [n_set, len(valid_ec), n_rti, len(valid_rti), len(wall), n_inv]
    checks.append(("C3 item 2: E_C/kT wall prints set_feasible/rti_feasible counts over VALID rows plus the bin total "
                    f"and invalid count (expected {n_set}/{len(valid_ec)} valid, {n_rti}/{len(valid_rti)} valid, "
                    f"of {len(wall)}, {n_inv} invalid)", ok))


def check_c3_rc_core_only(text, core, checks):
    """Item 3: the as-built tau_RC/delivered_step_fraction span is labelled
    CORE rows only, and the as-built sensitivity-row span (C_parasitic_F cut
    included) is printed separately with its worst row."""
    as_built = [r for r in core if r.get("row_kind") == "core" and r.get("family") == "horizontal_as_built"
                and r.get("regime") == "deterministic_pair" and isfin(r.get("tau_RC_ns"))]
    rs = {f(r.get("R_s_ohm")) for r in as_built}
    sens = [r for r in core if r.get("row_kind") == "sensitivity" and r.get("family") == "horizontal_as_built"
            and r.get("regime") == "deterministic_pair" and any(close_nan_safe(f(r.get("R_s_ohm")), v) for v in rs)
            and isfin(r.get("tau_RC_ns")) and isfin(r.get("delivered_step_fraction"))]
    body = _section(text, "## RC caveat")
    num = r"(-?\d+(?:\.\d+)?(?:[eE][+-]?\d+)?)"
    ok = False
    if as_built and sens:
        c_tau = [f(r["tau_RC_ns"]) for r in as_built]
        c_dsf = [f(r["delivered_step_fraction"]) for r in as_built if isfin(r.get("delivered_step_fraction"))]
        s_tau = [f(r["tau_RC_ns"]) for r in sens]
        s_dsf = [f(r["delivered_step_fraction"]) for r in sens]
        worst = max(sens, key=lambda r: f(r["tau_RC_ns"]))
        mc = re.search(rf"CORE rows only[^\n]*?\): tau_RC_ns spans {num}-{num}, delivered_step_fraction spans {num}-{num}", body)
        ms = re.search(rf"SENSITIVITY rows[^\n]*?\): tau_RC_ns spans {num}-{num}, delivered_step_fraction spans "
                        rf"{num}-{num}; the longest is row (\S+) ", body)
        near = lambda a, b: close_nan_safe(f(a), b, rtol=1e-3)
        ok = (bool(mc) and bool(ms)
              and near(mc.group(1), min(c_tau)) and near(mc.group(2), max(c_tau))
              and near(mc.group(3), min(c_dsf)) and near(mc.group(4), max(c_dsf))
              and near(ms.group(1), min(s_tau)) and near(ms.group(2), max(s_tau))
              and near(ms.group(3), min(s_dsf)) and near(ms.group(4), max(s_dsf))
              and ms.group(5) == worst["row_id"])
        tag = (f"core {min(c_tau):.4g}-{max(c_tau):.4g}; sensitivity {min(s_tau):.4g}-{max(s_tau):.4g} ns, "
               f"fraction down to {min(s_dsf):.4g} at {worst['row_id']}")
    else:
        tag = "no as-built rows"
    checks.append(("C3 item 3: the as-built RC span is labelled core rows only and the as-built sensitivity-row span "
                    f"is printed separately, both matching sweep.csv ({tag})", ok))


def check_c3_limitations_coverage(text, core, checks):
    """Item 4: the Limitations paragraph states per-cut coverage from the
    data -- no blanket 'both regimes/bounds, T in {230,300} K, both rates'
    claim while narrower cuts exist; every narrower cut is named with its
    actual families/regimes/bounds/temperatures/rates."""
    cov = independent_cut_coverage(core)
    body = _section(text, "## Limitations")
    narrow = {ax: c for ax, c in cov.items() if not c["full"]}
    ok = bool(body)
    if narrow:
        ok = ok and "both regimes/bounds, T in {230,300} K, both rates" not in body
    for ax, c in narrow.items():
        mo = re.search(re.escape(ax) + r" \((families=[^)]*)\)", body)
        if not mo:
            ok = False
            continue
        kv = dict(part.split("=", 1) for part in mo.group(1).split("; ") if "=" in part)
        ok = ok and all(_c3_norm_set(kv.get(k, "")) == c[k]
                        for k in ("families", "regimes", "strain_bounds", "T_hs_K", "rep_rate_hz"))
    for ax, c in cov.items():
        if c["full"]:
            ok = ok and re.search(r"cuts \([^)]*\b" + re.escape(ax) + r"\b[^)]*\) span every family", body) is not None
    checks.append(("C3 item 4: Limitations states per-cut coverage from the rows (narrower cuts: "
                    + (", ".join(sorted(narrow)) or "none") + ")", ok))


def check_c3_family_specific_ranking(text, core, checks):
    """Item 5: the per-family ranking has an entry only for cuts with rows
    in that family, and the page names the family-specific cuts."""
    sens = [r for r in core if r.get("row_kind") == "sensitivity"
            and (r.get("sensitivity_axis") or "") not in C3_NON_CUT_AXES]
    fams = {}
    for r in sens:
        fams.setdefault(r["sensitivity_axis"], set()).add(r.get("family"))
    specific = {ax: fs for ax, fs in fams.items() if fs != set(m.FAMILIES)}
    body = _section(text, "## Numerical sensitivity table")
    mo = re.search(r"Family-specific cuts \([^\n]*?\): ([^\n]*)", body)
    printed = dict((a, set(x.split(", "))) for a, x in re.findall(r"(\w+) \[([^\]]+) only\]", mo.group(1))) if mo else None
    ok = printed == specific if specific else (mo is None and "Family-specific cuts: none" in body)
    for fam in m.FAMILIES:
        sec = body.split(f"### {fam}", 1)
        fb = sec[1].split("\n### ", 1)[0] if len(sec) > 1 else ""
        ml = re.search(r"Ranked by commanded flux \|ratio-1\| \(measured leverage\): ([^\n]+)", fb)
        ranked = set(re.findall(r"(\w+)\([^)]*\)", ml.group(1))) if ml else None
        want = {ax for ax, fs in fams.items() if fam in fs and ax != "current_pulse_width"}
        ok = ok and ranked == want
    dw = _section(text, "## Numerical sensitivity table")
    ok = ok and "ranked, per family/axis above" not in dw
    checks.append(("C3 item 5: per-family ranking entries exist only where the cut has rows for that family and "
                    "the family-specific cuts are named ("
                    + "; ".join(f"{a}: {sorted(v)}" for a, v in sorted(specific.items())) + ")", ok))


def check_c3_reason_label_and_replay(text, core, checks):
    """Item 6: every printed invalid-reason label equals the independent
    label (whole reason when its ':' is parenthesised) with balanced
    parentheses, and the field-collapse replay says its external field is the
    replay input while the row's stored external_field_kVcm is nan."""
    counts = {}
    for r in core:
        if as_bool(r.get("valid")) is True:
            continue
        for reason in _reasons_of(r):
            k = independent_reason_label(reason)
            counts[k] = counts.get(k, 0) + 1
    body = _section(text, "## Invalid rows")
    printed = re.findall(r"^- (.+): (\d+)$", body, flags=re.M)
    labels = [p for p, _ in printed]
    quoted = re.findall(r"(?:Dominant invalid reason: |carrying )'([^']*)'", body)
    bal = lambda sx: sx.count("(") == sx.count(")")
    ok = (dict((p, int(n)) for p, n in printed) == counts and all(bal(x) for x in labels + quoted)
          and all(q in counts for q in quoted) and len(quoted) > 0) if counts else True
    ind = independent_invalid_triggers(core)
    ok_replay = True
    if ind["fc"] and ind["fc_rows"]:
        sample = min(ind["fc_rows"], key=lambda r: r["row_id"])
        fb = _section(text, "### Field-collapse (QCSE collapse) paragraph")
        if not isfin(sample.get("field_kVcm")):
            stored = "nan" if not isfin(sample.get("external_field_kVcm")) else f"{f(sample['external_field_kVcm']):.4g}"
            ok_replay = (f"stored external_field_kVcm is {stored}" in fb and "as its replay input" in fb
                         and "the replay uses external_field_kVcm=0.0" not in fb)
    checks.append(("C3 item 6: invalid-reason labels are printed whole with balanced parentheses and match an "
                    "independent recount, and the field-collapse replay input is distinguished from the row's "
                    "stored external_field_kVcm", ok and ok_replay))


def check_c3_downstream_hashes(out, man, checks):
    """Item 7: graphs.html is written after the manifest by its own builder,
    so the manifest must either exclude it (declared) or carry a hash that
    matches disk -- rebuild order must not matter."""
    hashes = man.get("output_hashes", {})
    excluded = man.get("output_hashes_excluded_downstream", [])
    if "graphs.html" in hashes:
        ok = (out / "graphs.html").is_file() and \
            hashlib.sha256((out / "graphs.html").read_bytes()).hexdigest() == hashes["graphs.html"]
    else:
        ok = "graphs.html" in excluded
    ok = ok and list(getattr(m, "DOWNSTREAM_ARTIFACTS", ())) == ["graphs.html"]
    checks.append(("C3 item 7: manifest.json deliberately excludes the downstream graphs.html from output_hashes "
                    "(or its recorded hash matches disk), so a graphs-page rebuild never stales the manifest", ok))


def check_c3_prose(text, core, checks):
    check_c3_reduced_cut_list(text, core, checks)
    check_c3_ec_wall_denominator(text, core, checks)
    check_c3_rc_core_only(text, core, checks)
    check_c3_limitations_coverage(text, core, checks)
    check_c3_family_specific_ranking(text, core, checks)
    check_c3_reason_label_and_replay(text, core, checks)


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
    check_footer_bbox_fixture(checks)
    check_csv_nan_bool_coercion_fixture(checks)
    check_invalid_paragraph_trigger_fixture(checks)
    check_grid_self(checks, detail)
    check_reduced_cut_axes_declared(checks, detail)
    check_loading_window_axis(checks, detail)
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
        check_rti_transport_feasible_denominator(text, core, checks)
        check_dipole_weights_ranking(text, core, checks)
        check_results_text_obligations(text, checks)
        check_c4_prose(text, core, checks)
        check_c3_prose(text, core, checks)
        check_c3_downstream_hashes(out, man, checks)
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
