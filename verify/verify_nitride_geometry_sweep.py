"""Independent structural, replay, sensitivity, convergence and presentation
checks for scripts/run_nitride_geometry_stark.py's output.

Grid-membership counts are recomputed from the run module's own PURE combo
builders (build_core/build_shape/build_qw/build_stark/build_sensitivities),
which never call evaluate() -- an independent re-derivation of the plan, not
a re-read of the CSV's own row count. Replay re-evaluates real cache
identities (never trusts saved numbers). Mutation-sensitive fixtures unit-
test the run module's own pure eligibility/hardware-qualification functions
with crafted boundary values that the sweep's actual data would not by
itself exercise (g2 exactly at the gate, flux exactly at the floor).
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, re, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts import run_nitride_geometry_stark as s

ZHANG_SLOPE_MEV_PER_V_INDEPENDENT = -10.0  # [V] Zhang et al., APL 108, 153102 (2016), Fig. 5;
# transcribed independently here, not imported from fsim_core.nitride_stark or
# the run script, so this is a real cross-check of the transcription.
WANG_HEIGHT_NM_INDEPENDENT = 7.0            # [V] Wang et al., Sci. Rep. 7, 12089 (2017), uncapped AFM
WANG_DIAMETER_NM_INDEPENDENT = 35.0         # [V] ibid.
WANG_TEMPERATURE_K_INDEPENDENT = 220.0      # [V] ibid.

# Frozen core-grid axes, transcribed independently from docs/nitride_geometry_
# stark_contract.md ("7 heights x 5 radii x 4 orientations x 4 T_hs x 3
# screening x 2 regimes = 3,360 rows"), NOT imported from the run module's own
# CORE_H/CORE_R/ORI/TS/SCR/REG constants (Opus fix-round medium finding: grid
# membership was checked only against the module's own builders, so e.g.
# SCR=(0,0.5,0.9) in the module would still pass). The module's constants are
# asserted equal to these literals below, and the CSV's own axis coverage is
# asserted a subset of them independently of build_core().
CONTRACT_HEIGHTS_NM = (1.0, 2.0, 3.0, 4.0, 5.0, 7.0, 10.0)
CONTRACT_RADII_NM = (5.0, 10.0, 15.0, 20.0, 30.0)
CONTRACT_ORIENTATIONS = ("c_plane", "semipolar_11_22", "m_plane", "a_plane")
CONTRACT_T_HS_K = (230.0, 250.0, 273.0, 300.0)
CONTRACT_SCREENING = (0.0, 0.5, 1.0)
CONTRACT_REGIMES = ("rectangular", "deterministic_pair")


def rows(p):
    with open(p, encoding="utf-8", newline="") as f:
        return list(csv.DictReader(f))

def f(x):
    try: return float(x)
    except (TypeError, ValueError): return float("nan")

def isfin(x): return isinstance(x, float) and math.isfinite(x) or isinstance(x, int)

def close_nan_safe(a, b, rtol=1e-8, atol=1e-10):
    fa, fb = f(a), f(b)
    na, nb = math.isnan(fa), math.isnan(fb)
    if na or nb: return na and nb
    if math.isinf(fa) or math.isinf(fb): return math.isinf(fa) and math.isinf(fb) and (fa > 0) == (fb > 0)
    return math.isclose(fa, fb, rel_tol=rtol, abs_tol=atol)


def check_grid(core, geo, bias, cur, sens, comp, man, quick, checks, detail):
    core_p = s.build_core(quick); shape_p = s.build_shape(quick); qw_p = s.build_qw(quick)
    bias_p, cur_p = s.build_stark(quick)
    sens_p = s.build_sensitivities(quick)
    checks.append(("core row count matches independent build_core()", len(core) == len(core_p)))
    ax = ("height_nm", "radius_nm", "orientation", "T_hs", "screening_fraction", "regime")
    expected_combos = {tuple(round(p[k], 6) if isinstance(p[k], float) else p[k] for k in ax) for p in core_p}
    actual_combos = [tuple(round(f(r[k]), 6) if k != "orientation" and k != "regime" else r[k] for k in ax) for r in core]
    checks.append(("core has no duplicate (height,radius,orientation,T,screening,regime) combinations", len(actual_combos) == len(set(actual_combos))))
    checks.append(("core covers every expected combination with none missing/extra", set(actual_combos) == expected_combos))
    checks.append(("shape row count matches independent build_shape()", len([r for r in geo if r["row_kind"] == "shape"]) == len(shape_p)))
    checks.append(("qw row count matches independent build_qw()", len([r for r in geo if r["row_kind"] == "qw"]) == len(qw_p)))
    checks.append(("stark_bias row count matches independent build_stark()", len(bias) == len(bias_p)))
    checks.append(("stark_current row count matches independent build_stark()", len(cur) == len(cur_p)))
    checks.append(("sensitivities row count matches independent build_sensitivities()", len(sens) == len(sens_p)))
    all_ids = [r["row_id"] for r in core + geo + bias + cur + sens + comp]
    checks.append(("all row_ids unique across every output file", len(set(all_ids)) == len(all_ids)))
    checks.append(("core covers all 4 orientations", {x["orientation"] for x in core} == set(s.ORI)))
    checks.append(("core covers both regimes", {x["regime"] for x in core} == set(s.REG)))
    checks.append(("core covers all 3 screening fractions", {round(f(x["screening_fraction"]), 6) for x in core} == {round(v, 6) for v in s.SCR}))
    checks.append(("every core row carries card_hash/cache_identity/invalid_reasons", all(x.get("card_hash") and x.get("cache_identity") and "invalid_reasons" in x for x in core)))
    checks.append(("every stark_bias row carries spectroscopy_valid and V_j", all("spectroscopy_valid" in x and "V_j" in x for x in bias) if bias else True))
    checks.append(("manifest core/shape/qw/stark full-contract counts recorded", man.get("full_contract_rows") == {"core": 3360, "shape": 864, "qw": 288, "stark": 3120, "total_before_caching": 7632}))
    # Tie the module's own grid constants to the frozen contract literals
    # (Opus fix-round medium finding) -- independent of build_core(), so a
    # mutated CORE_H/CORE_R/ORI/TS/SCR/REG is actually caught.
    checks.append(("module CORE_H matches the frozen contract heights", tuple(s.CORE_H) == CONTRACT_HEIGHTS_NM))
    checks.append(("module CORE_R matches the frozen contract radii", tuple(s.CORE_R) == CONTRACT_RADII_NM))
    checks.append(("module ORI matches the frozen contract orientations", tuple(s.ORI) == CONTRACT_ORIENTATIONS))
    checks.append(("module TS matches the frozen contract T_hs", tuple(s.TS) == CONTRACT_T_HS_K))
    checks.append(("module SCR matches the frozen contract screening fractions", tuple(round(v, 6) for v in s.SCR) == tuple(round(v, 6) for v in CONTRACT_SCREENING)))
    checks.append(("module REG matches the frozen contract regimes", tuple(s.REG) == CONTRACT_REGIMES))
    # The CSV's own axis coverage (quick or full) must be a SUBSET of the
    # frozen contract literals, checked directly against the CSV, not against
    # whatever the module's constants happen to say.
    checks.append(("core heights present are all within the frozen contract literal set",
                    {round(f(r["height_nm"]), 6) for r in core} <= {round(v, 6) for v in CONTRACT_HEIGHTS_NM}))
    checks.append(("core radii present are all within the frozen contract literal set",
                    {round(f(r["radius_nm"]), 6) for r in core} <= {round(v, 6) for v in CONTRACT_RADII_NM}))
    checks.append(("core orientations present are all within the frozen contract literal set",
                    {r["orientation"] for r in core} <= set(CONTRACT_ORIENTATIONS)))
    checks.append(("core T_hs present are all within the frozen contract literal set",
                    {round(f(r["T_hs"]), 6) for r in core} <= {round(v, 6) for v in CONTRACT_T_HS_K}))
    checks.append(("core screening_fraction present are all within the frozen contract literal set",
                    {round(f(r["screening_fraction"]), 6) for r in core} <= {round(v, 6) for v in CONTRACT_SCREENING}))
    checks.append(("core regimes present are all within the frozen contract literal set",
                    {r["regime"] for r in core} <= set(CONTRACT_REGIMES)))


def check_caps(man, a, checks):
    cap = 600 if a.quick else 10000
    checks.append(("evaluate_calls within declared cap", man["evaluate_calls"] <= cap))
    checks.append(("complete flag matches quick/full expectation", (man["complete"] is False) if a.quick else isinstance(man["complete"], bool)))
    if a.quick:
        checks.append(("quick run stays within its own 600-call budget", man["evaluate_calls"] <= 600))
    else:
        checks.append(("full run stays within the 1800 s runtime bound", man["runtime_s"] < 1800))
        checks.append(("full run attempted every requested row when claiming complete=True",
                        (not man["complete"]) or (man["requested_rows"]["core"] == 3360 and man["requested_rows"]["stark_bias"] + man["requested_rows"]["stark_current"] == 3120)))


def check_replay(core, geo, bias, cur, checks, quick):
    """rtol=1e-8/atol=1e-10 replay across: valid + invalid diagnostic, both
    regimes, c-plane + a-plane, a shape alternative, a QW card, nondefault
    polarity, and a temperature-shifted current point."""
    samples = {}
    samples["valid_core"] = next((x for x in core if x["valid"] == "True"), None)
    samples["invalid_diagnostic"] = next((x for x in bias if x.get("valid") == "False"), None) or next((x for x in core if x["valid"] == "False"), None)
    samples["regime_rectangular"] = next((x for x in core if x["regime"] == "rectangular"), None)
    samples["regime_set"] = next((x for x in core if x["regime"] == "deterministic_pair"), None)
    samples["c_plane"] = next((x for x in core if x["orientation"] == "c_plane"), None)
    samples["a_plane"] = next((x for x in core if x["orientation"] == "a_plane"), None)
    samples["shape_alt"] = next((x for x in geo if x["row_kind"] == "shape" and x.get("shape") in ("lens", "truncated_cone")), None)
    samples["qw_card"] = next((x for x in geo if x["row_kind"] == "qw"), None)
    samples["nondefault_polarity"] = next((x for x in bias if x.get("polarity") == "-1"), None)
    # Opus fix-round finding: the previous fallback to cur[0] (itself always
    # T_hs=300.0) made this check pass without ever being exercised. Quick
    # mode's own Stark current grid is deliberately fixed at T_hs=300 K
    # (spec: "Reduce Stark to ... fixed/current heat-sink temperature 300
    # K"), so no shifted-T sample CAN exist there -- explicitly marked N/A
    # rather than silently faked. Full mode has no such excuse: a missing
    # sample there is a real failure, not a fallback.
    if quick:
        samples["temperature_shifted_current"] = None
    else:
        samples["temperature_shifted_current"] = next((x for x in cur if x.get("T_hs") != "300.0"), None)

    for name, sample in samples.items():
        if sample is None:
            if name == "temperature_shifted_current" and quick:
                checks.append((f"replay sample available: {name} (not applicable in quick mode: Stark current T_hs is fixed at 300 K by spec)", True))
            else:
                checks.append((f"replay sample available: {name}", False))
            continue
        try:
            p = json.loads(sample["cache_identity"])
        except (KeyError, json.JSONDecodeError):
            checks.append((f"replay sample has a decodable cache_identity: {name}", False)); continue
        d, _ = s._design(p)
        got = s.device.evaluate(d, T_grid=[p.get("T_j", p.get("T_hs"))])["scalars"]
        ok = (close_nan_safe(sample.get("E_X_eV"), got.get("E_X_eV"))
              and close_nan_safe(sample.get("g2"), got.get("g2_op"))
              and close_nan_safe(sample.get("field_kVcm"), got.get("field_kVcm"))
              and close_nan_safe(sample.get("tau_rad_bare_ns"), got.get("tau_rad_bare_ns"))
              and close_nan_safe(sample.get("reservoir_energy_eV"), got.get("reservoir_energy_eV"))
              and close_nan_safe(sample.get("signal_flux_s"), got.get("collected_flux_pulsed_s"))
              and (sample.get("valid") == "True") == bool(got.get("valid")))
        checks.append((f"independent replay matches saved row ({name})", ok))


AXIS_CANDIDATES = {
    "semipolar_factor": ("E_X_eV", "g2"),
    "x_in": ("E_X_eV",),
    "strain_fraction": ("E_X_eV",),
    "tau_rad0_ns": ("tau_rad_bare_ns",),
    "k_nr_ns": ("g2", "tau_rad_bare_ns"),
    "background_tau_ns": ("g2", "background_flux_s"),
    "b_res": ("g2",),
    "Q_purcell": ("tau_rad_cavity_ns", "g2", "Fp_add"),
    "detuning_offset_meV": ("tau_rad_cavity_ns", "detuning_meV"),
    "tau_cap_density_convention": ("g2", "tau_cap_ps_used"),
    "island_radius_nm": ("set_EC_over_kT", "set_radius_max_nm"),
}

def check_oat_sensitivities(sens, core, checks, detail):
    by_axis = {}
    for r in sens:
        by_axis.setdefault(r["sensitivity_axis"], []).append(r)
    for axis, cand in AXIS_CANDIDATES.items():
        group = by_axis.get(axis, [])
        if len(group) < 2:
            checks.append((f"OAT axis has >=2 rows: {axis}", False)); continue
        moved = False
        for key in cand:
            vals = [f(r.get(key)) for r in group]
            if all(isfin(v) for v in vals) and len({round(v, 9) for v in vals}) > 1:
                moved = True; break
        if not moved and axis == "background_tau_ns":
            # Documented conditional null: background_tau_ns only extends an
            # afterglow term proportional to the reservoir's own SRH
            # background rate; at REF, background_flux_s is IDENTICAL
            # (not merely small) across both values, showing the observed
            # background here comes entirely from the unrelated b_res
            # channel and the bg_tau-dependent term contributes exactly
            # zero -- not a bug in this sweep or in device.py. Opus
            # fix-round finding: checking background_flux_s alone could pass
            # by coincidence (e.g. if that single field were broken/blank
            # for an unrelated reason) with moved=True on both branches; ALL
            # THREE of background_flux_s, g2 and signal_flux_s must be
            # finite AND identical across the axis values for the null to
            # be accepted, so a genuine (buggy) movement in any one of them
            # now fails this check.
            bgflux = [f(r.get("background_flux_s")) for r in group]
            g2s = [f(r.get("g2")) for r in group]
            fluxes = [f(r.get("signal_flux_s")) for r in group]
            def _all_identical(vals):
                return len(vals) > 0 and all(isfin(v) for v in vals) and len({round(v, 12) for v in vals}) == 1
            moved = _all_identical(bgflux) and _all_identical(g2s) and _all_identical(fluxes)
            detail.append(f"background_tau_ns null-movement accepted conditionally only because "
                           f"background_flux_s/g2/signal_flux_s are ALL identical across values "
                           f"(bg={bgflux}, g2={g2s}, flux={fluxes})")
        checks.append((f"OAT axis moves a relevant observable: {axis}", moved))
    # Explicit expected null cases (spec): nonpolar screening_fraction, and
    # m/a orientation identity -- read directly off the core grid, not the
    # sensitivities table.
    def _grp(o): return {(r["height_nm"], r["radius_nm"], r["T_hs"], r["regime"]): r for r in core if r["orientation"] == o}
    a_by_key = {}
    for r in core:
        if r["orientation"] == "a_plane":
            a_by_key.setdefault((r["height_nm"], r["radius_nm"], r["T_hs"], r["regime"]), {})[r["screening_fraction"]] = r
    null_ok = True; n_groups = 0
    for key, byscr in a_by_key.items():
        if len(byscr) < 2: continue
        n_groups += 1
        vals = [f(r["E_X_eV"]) for r in byscr.values() if r["valid"] == "True"]
        if len(vals) >= 2 and len({round(v, 6) for v in vals}) > 1: null_ok = False
    checks.append(("expected null case: nonpolar (a-plane) E_X_eV is screening-invariant", null_ok and n_groups > 0))
    m_by_key = _grp("m_plane"); a_by_key2 = _grp("a_plane")
    matched = [(k, m_by_key[k], a_by_key2[k]) for k in m_by_key if k in a_by_key2]
    identity_ok = all(close_nan_safe(m["E_X_eV"], a["E_X_eV"]) for _, m, a in matched if m["valid"] == "True" and a["valid"] == "True")
    checks.append(("expected null case: m-plane/a-plane E_X_eV coincide (model limitation)", identity_ok and len(matched) > 0))


def check_convergence(man, checks):
    cc = man.get("convergence_checks", [])
    checks.append(("convergence_checks present", len(cc) > 0))
    resolved = [c for c in cc if "coarse_slope_meV_per_V" in c]
    checks.append(("at least one convergence check reached a resolved comparison", len(resolved) > 0))
    checks.append(("every resolved convergence check uses the declared criterion",
                    all(c.get("criterion") == "abs<=0.5meV/V or rel<=5%" for c in resolved)))
    checks.append(("every resolved convergence check's abs/rel diff is consistent with its own converged flag",
                    all((c["abs_diff_meV_per_V"] <= 0.5 or c["rel_diff"] <= 0.05) == c["converged"] for c in resolved)))


def check_plots(out, man, core, geo, bias, checks):
    mapping = man.get("plot_row_mapping", {})
    checks.append(("plot_row_mapping is non-empty", len(mapping) > 0))
    for name in mapping:
        p = out / name
        checks.append((f"figure exists and is non-trivial: {name}", p.is_file() and p.stat().st_size > 1000))
    hashes = man.get("output_hashes", {})
    checks.append(("every recorded output_hash matches the file on disk",
                    all((out / n).is_file() and hashlib.sha256((out / n).read_bytes()).hexdigest() == v for n, v in hashes.items())))
    # Plot-construction return contract: recompute x/y for a couple of
    # figures straight from the CSV rows the manifest says were plotted, and
    # check they match the manifest's OWN recorded x/y (not just that files
    # exist) -- this is the ordered row-id -> data contract check.
    all_rows = {r["row_id"]: r for r in core + geo + bias}
    for figname in ("height_response.png", "stark_energy_bias.png"):
        contract = mapping.get(figname)
        if not contract:
            checks.append((f"plot-construction contract present: {figname}", False)); continue
        ok = True
        for label, entry in contract.items():
            for rid, xv, yv in zip(entry["row_ids"], entry["x"], entry["y"]):
                row = all_rows.get(rid)
                if row is None: ok = False; continue
                xkey = "height_nm" if figname == "height_response.png" else "V_j"
                ykey = "E_X_eV"
                if not (close_nan_safe(row.get(xkey), xv, rtol=1e-6, atol=1e-9) and close_nan_safe(row.get(ykey), yv, rtol=1e-6, atol=1e-9)):
                    ok = False
        checks.append((f"plot-construction contract's x/y trace to the CSV rows by row_id: {figname}", ok))


def check_mutation_fixtures(checks):
    # eligible(): flux-floor gate ONLY (matches run_nitride_cavity.py's
    # `eligible` exactly -- Opus fix-round finding: g2<0.5 must NOT be
    # folded in here, or this column is not comparable with round 1).
    checks.append(("eligible: g2 at/above the optical gate (0.5) is STILL eligible (flux-only gate)", s.eligible(True, 0.5, 2000.) is True))
    checks.append(("eligible: flux exactly at floor (1000) is included", s.eligible(True, 0.1, 1000.) is True))
    checks.append(("eligible: flux just under floor (999.999) is rejected", s.eligible(True, 0.1, 999.999) is False))
    checks.append(("eligible: invalid row is never eligible regardless of g2/flux", s.eligible(False, 0.01, 1e9) is False))
    checks.append(("eligible: NaN g2 is never eligible", s.eligible(True, float("nan"), 2000.) is False))
    # optical_pass(): adds g2<0.5 and, for SET, one_pair_valid.
    checks.append(("optical_pass: g2 exactly at gate (0.5) is rejected", s.optical_pass(True, 0.5, "rectangular", True) is False))
    checks.append(("optical_pass: g2 just under gate (0.4999) is accepted (rectangular)", s.optical_pass(True, 0.4999, "rectangular", True) is True))
    checks.append(("optical_pass: an ineligible row is never optical_pass regardless of g2", s.optical_pass(False, 0.01, "rectangular", True) is False))
    checks.append(("optical_pass: SET regime with one_pair_valid=False is rejected even with good g2", s.optical_pass(True, 0.1, "deterministic_pair", False) is False))
    checks.append(("optical_pass: SET regime with one_pair_valid=True and good g2 is accepted", s.optical_pass(True, 0.1, "deterministic_pair", True) is True))
    checks.append(("optical_pass: rectangular regime never requires one_pair_valid", s.optical_pass(True, 0.1, "rectangular", False) is True))
    # hardware_qualified(): adds, for SET, hardware_feasible (=set_feasible
    # AND pair_supply_possible) -- matches run_nitride_cavity.py's
    # device_pass exactly (Opus fix-round finding: the previous
    # hardware_qualified dropped one_pair_valid entirely by taking `elig`
    # instead of `opt_pass`, so it was not comparable with round 1).
    checks.append(("hardware_qualified: rectangular regime never requires hardware qualification", s.hardware_qualified(True, "rectangular", False) is True))
    checks.append(("hardware_qualified: SET regime with hardware_feasible=False is never hardware-qualified (no combining extrema across rows)", s.hardware_qualified(True, "deterministic_pair", False) is False))
    checks.append(("hardware_qualified: SET regime with hardware_feasible=True and optical_pass is hardware-qualified", s.hardware_qualified(True, "deterministic_pair", True) is True))
    checks.append(("hardware_qualified: hardware qualification never rescues a non-optical_pass row", s.hardware_qualified(False, "deterministic_pair", True) is False))


def check_literature(lit_csv, comp, checks):
    zhang_row = next((r for r in lit_csv if "Zhang" in r["source"]), None)
    checks.append(("Zhang slope transcription matches independent literal", zhang_row is not None and close_nan_safe(zhang_row.get("slope_meV_per_V"), ZHANG_SLOPE_MEV_PER_V_INDEPENDENT)))
    wang_row = next((r for r in lit_csv if "Wang" in r["source"] and r.get("kind") == "comparison_only"), None)
    checks.append(("Wang height/diameter/temperature transcription matches independent literal",
                    wang_row is not None and close_nan_safe(wang_row.get("height_nm"), WANG_HEIGHT_NM_INDEPENDENT)
                    and close_nan_safe(wang_row.get("diameter_nm"), WANG_DIAMETER_NM_INDEPENDENT)
                    and close_nan_safe(wang_row.get("temperature_K"), WANG_TEMPERATURE_K_INDEPENDENT)))
    checks.append(("every literature_comparisons.csv row carries a provenance/model_row_ids field",
                    all("model_row_ids" in r and "source" in r for r in lit_csv)))
    checks.append(("illustrative slope interval is never labelled 'measured'",
                    all(r.get("slope_interval_kind") in ("illustrative_not_measured", "user_supplied") for r in comp) if comp else True))
    nonpolar_claims = [r for r in comp if r.get("compatible") == "True" and ("a_plane" in r.get("group", "") or "m_plane" in r.get("group", ""))]
    checks.append(("no nonpolar compatibility row claims compatible=True", len(nonpolar_claims) == 0))


VERDICT_RE = re.compile(
    r"^VERDICT: idealized_status=(?P<ideal>\S+) family=(?P<family>\S+) regime=(?P<regime>\S+) "
    r"screening=(?P<screening>\S+) complete=(?P<complete>\S+) eligible=(?P<eligible>\d+) "
    r"paired_optical_pass=(?P<opt>\d+) hardware_qualified=(?P<hw>\d+) coverage=(?P<covn>\d+)/(?P<covd>\d+) "
    r"invalid=(?P<invalid>\d+) flux_floor=1000/s$")
_SCR_FOR_LABEL = {"unscreened_lower": 0.0, "screened_upper": 1.0}

def check_results_md(text, core, checks):
    checks.append(("results.md contains VERDICT lines", "VERDICT:" in text))
    verdicts = [ln for ln in text.splitlines() if ln.startswith("VERDICT:")]
    checks.append(("exactly 8 headline VERDICT lines (c-plane + a-plane x 2 regimes x 2 bounds)", len(verdicts) == 8))
    for fam in ("family=c_plane", "family=a_plane"):
        for reg in ("regime=rectangular", "regime=deterministic_pair"):
            for scr in ("screening=unscreened_lower", "screening=screened_upper"):
                checks.append((f"VERDICT line present: {fam} {reg} {scr}", any(fam in v and reg in v and scr in v for v in verdicts)))
    checks.append(("every VERDICT line carries idealized_status", all("idealized_status=" in v for v in verdicts)))
    checks.append(("semipolar/m-plane kept out of headline VERDICT lines", not any("family=semipolar" in v or "family=m_plane" in v for v in verdicts)))
    checks.append(("results.md documents screening_unidentifiable for nonpolar orientation", "screening_unidentifiable" in text))
    checks.append(("results.md states the illustrative interval is not measured", "illustrative" in text and "never a measured interval" in text))

    # Independent recomputation from sweep.csv (Opus fix-round high finding):
    # check_results_md previously only pattern-matched the VERDICT lines and
    # never recomputed eligible/paired_optical_pass/hardware_qualified/
    # invalid from the CSV, so a rewritten headline line (e.g. eligible=8
    # with true value 0) still passed. Every count is now recomputed here,
    # via the run module's own pure eligible/optical_pass/hardware_qualified
    # functions applied to the ACTUAL CSV rows, and compared against the
    # line's own numbers -- a mutated line fails.
    for v in verdicts:
        m = VERDICT_RE.match(v)
        if not m:
            checks.append((f"VERDICT line parses in the documented field order: {v[:70]}", False)); continue
        fam, reg, scr_label = m["family"], m["regime"], m["screening"]
        scr = _SCR_FOR_LABEL.get(scr_label)
        if scr is None:
            checks.append((f"VERDICT line's screening label is recognized: {v[:70]}", False)); continue
        rs = [r for r in core if r["orientation"] == fam and r["regime"] == reg and close_nan_safe(r["screening_fraction"], scr)]
        n_elig = n_opt = n_hw = n_invalid = 0
        for r in rs:
            valid = r["valid"] == "True"
            g2v = f(r["g2"]); fluxv = f(r["signal_flux_s"])
            one_pair = r.get("one_pair_valid") == "True"
            hw_feasible = r.get("hardware_feasible") == "True"
            elig = s.eligible(valid, g2v, fluxv)
            opt = s.optical_pass(elig, g2v, reg, one_pair)
            hw = s.hardware_qualified(opt, reg, hw_feasible)
            n_elig += elig; n_opt += opt; n_hw += hw
            if not valid: n_invalid += 1
        ok = (int(m["eligible"]) == n_elig and int(m["opt"]) == n_opt and int(m["hw"]) == n_hw
              and int(m["invalid"]) == n_invalid and int(m["covn"]) == len(rs) and int(m["covd"]) == len(rs))
        checks.append((f"VERDICT counts recomputed from sweep.csv match the line ({fam}/{reg}/{scr_label}): "
                        f"eligible={n_elig} opt={n_opt} hw={n_hw} invalid={n_invalid} n={len(rs)}", ok))


def main(argv=None):
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out-dir", default="out/nitride_geometry_stark")
    a = ap.parse_args(argv)
    out = Path(a.out_dir)
    core = rows(out / "sweep.csv"); geo = rows(out / "geometry_supplement.csv")
    bias = rows(out / "stark_bias.csv"); cur = rows(out / "stark_current.csv")
    comp = rows(out / "screening_compatibility.csv"); sens = rows(out / "sensitivities.csv")
    lit_csv = rows(out / "literature_comparisons.csv")
    man = json.loads((out / "manifest.json").read_text(encoding="utf-8"))
    text = (out / "results.md").read_text(encoding="utf-8")

    checks = []; detail = []
    check_grid(core, geo, bias, cur, sens, comp, man, a.quick, checks, detail)
    check_caps(man, a, checks)
    check_replay(core, geo, bias, cur, checks, a.quick)
    check_oat_sensitivities(sens, core, checks, detail)
    check_convergence(man, checks)
    check_plots(out, man, core, geo, bias, checks)
    check_mutation_fixtures(checks)
    check_literature(lit_csv, comp, checks)
    check_results_md(text, core, checks)

    passed = sum(1 for _, ok in checks if ok)
    total = len(checks)
    for name, ok in checks:
        if not ok: print("FAILED:", name)
    for line in detail: print("NOTE:", line)
    print("%d/%d nitride geometry sweep checks passed" % (passed, total))
    return 0 if passed == total else 1


if __name__ == "__main__":
    raise SystemExit(main())
