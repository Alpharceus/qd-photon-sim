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
    sens_p = s.build_all_sensitivities(quick); anchor_p = s.build_fixed_anchor(quick)
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
    sens_only = [r for r in sens if r.get("row_kind") == "sensitivity"]
    anchor_only = [r for r in sens if r.get("row_kind") == "fixed_anchor"]
    checks.append(("sensitivities row count matches independent build_all_sensitivities() "
                    "(fix round 3, item 5: all three baselines)", len(sens_only) == len(sens_p)))
    checks.append(("fixed_anchor row count matches independent build_fixed_anchor() (item 1)",
                    len(anchor_only) == len(anchor_p)))
    all_ids = [r["row_id"] for r in core + geo + bias + cur + sens + comp]
    checks.append(("all row_ids unique across every output file", len(set(all_ids)) == len(all_ids)))
    checks.append(("core covers all 4 orientations", {x["orientation"] for x in core} == set(s.ORI)))
    checks.append(("core covers both regimes", {x["regime"] for x in core} == set(s.REG)))
    checks.append(("core covers all 3 screening fractions", {round(f(x["screening_fraction"]), 6) for x in core} == {round(v, 6) for v in s.SCR}))
    checks.append(("every core row carries card_hash/cache_identity/invalid_reasons", all(x.get("card_hash") and x.get("cache_identity") and "invalid_reasons" in x for x in core)))
    checks.append(("every stark_bias row carries spectroscopy_valid and V_j", all("spectroscopy_valid" in x and "V_j" in x for x in bias) if bias else True))
    # Fix round 3, item 1/8: cavity_tracking column values -- every core row
    # (headline, re-tracked) is "per_T_hs"; every Stark bias/current row
    # (fixed cavity reference across the whole V_j/current sweep) is
    # "fixed_300K"; no row carries any other value.
    checks.append(("every row's cavity_tracking is one of per_T_hs/fixed_300K",
                    all(x.get("cavity_tracking") in ("per_T_hs", "fixed_300K") for x in core + geo + bias + cur)))
    checks.append(("every core row has cavity_tracking=per_T_hs (item 1: headline rows re-track the cavity)",
                    all(x.get("cavity_tracking") == "per_T_hs" for x in core)))
    checks.append(("every stark_bias row has cavity_tracking=fixed_300K (fixed cavity reference across the trace)",
                    all(x.get("cavity_tracking") == "fixed_300K" for x in bias) if bias else True))
    checks.append(("every stark_current row has cavity_tracking=fixed_300K",
                    all(x.get("cavity_tracking") == "fixed_300K" for x in cur) if cur else True))
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


def _synthetic_bias_rows(bias_p, corrupt_t_hs=True):
    """Build synthetic (never evaluate()-derived) Stark bias rows that
    reproduce fsim_core/device.py's own reporting convention for
    bias_mode=='junction_voltage' rows: scalars["T_hs"] is the CARD's fixed
    heat-sink setting -- identical across every T_j in the sweep -- while
    T_hs_requested (set by run_nitride_geometry_stark._row()) carries the
    real, uncorrupted sweep value. E_X_eV is a smooth, strictly monotone
    function of V_j and T_hs_requested so stark_derivatives has a
    well-posed trace to differentiate whenever the grouping is correct."""
    out = []
    for i, p in enumerate(bias_p):
        out.append({
            "row_id": "SYN%05d" % i,
            "height_nm": p["height_nm"], "radius_nm": p["radius_nm"],
            "orientation": p["orientation"], "screening_fraction": p["screening_fraction"],
            "regime": p["regime"], "x_in": p.get("x_in"),
            "geometry_type": p.get("geometry_type"), "shape": p.get("shape"),
            "top_radius_fraction": p.get("top_radius_fraction"),
            "field_polarity": p["polarity"], "bias_mode": p["bias_mode"],
            "T_hs": (300.0 if corrupt_t_hs else p["T_hs"]),
            "T_hs_requested": p["T_hs"],
            "V_j": p["V_j"],
            "E_X_eV": 2.2 - 0.01 * p["V_j"] - 1e-4 * p["T_hs"],
            "spectroscopy_valid": True,
        })
    return out


def check_derivative_trace_grouping(checks, detail):
    """Fix round 2: reproduce, without ever calling evaluate(), the exact
    condition that crashed attach_derivatives at scripts/run_nitride_
    geometry_stark.py:383 -> fsim_core/nitride_stark.py:64 after 787 s in
    the full run (ValueError: 'trace valid voltages must be strictly
    increasing and distinct'). Every planned full-mode Stark bias trace
    must group to distinct, strictly increasing V_j on its own (a sanity
    check on the plan itself, independent of evaluate()'s T_hs-corruption
    quirk), and attach_derivatives must run the full-mode plan AS
    evaluate() actually reports it (T_hs corruption reproduced) end to end
    without raising and without any trace failure."""
    bias_p, _ = s.build_stark(False)
    groups = {}
    for p in bias_p:
        groups.setdefault(s._derivative_group_key(p), []).append(p)
    checks.append(("full-mode Stark bias plan forms exactly 120 independent trace groups "
                    "(height x orientation x screening x regime x T_hs x polarity)", len(groups) == 120))
    all_distinct_increasing = True
    for plist in groups.values():
        vjs = sorted(pp["V_j"] for pp in plist)
        if len(vjs) != 17 or len(set(vjs)) != len(vjs) or any(vjs[i + 1] <= vjs[i] for i in range(len(vjs) - 1)):
            all_distinct_increasing = False
    checks.append(("every planned full-mode Stark bias trace has exactly 17 distinct, strictly increasing V_j (no evaluation)",
                    all_distinct_increasing))

    synth = _synthetic_bias_rows(bias_p, corrupt_t_hs=True)
    failures = s.attach_derivatives(synth, "V_j")
    checks.append(("attach_derivatives runs the full-mode synthetic plan (T_hs corruption reproduced) "
                    "without raising and reports zero trace failures", failures == []))
    n_resolved = sum(1 for r in synth if r.get("derivative_valid"))
    checks.append(("synthetic full-mode plan yields resolved (non-NaN) derivatives on interior points", n_resolved > 0))
    detail.append("synthetic full-mode Stark bias plan: %d groups, %d/%d rows derivative_valid, %d attach_derivatives failures"
                   % (len(groups), n_resolved, len(synth), len(failures)))

    # Negative control: a genuinely duplicated-voltage trace (two rows that
    # really do share every trace coordinate including V_j) must degrade
    # gracefully (derivative_valid=False, reason='non_monotonic_trace'),
    # never raise past attach_derivatives.
    dup = [dict(synth[0], row_id="DUP0", derivative_valid=None), dict(synth[0], row_id="DUP1", derivative_valid=None)]
    dup_failures = s.attach_derivatives(dup, "V_j")
    checks.append(("attach_derivatives degrades a genuinely duplicated-voltage trace instead of raising",
                    len(dup_failures) == 1 and dup_failures[0]["reason"] == "non_monotonic_trace"
                    and all(r.get("derivative_valid") is False for r in dup)
                    and all(r.get("derivative_invalid_reason") == "non_monotonic_trace" for r in dup)))

    # Confirms this reproduces the ACTUAL round-2 bug rather than a
    # strawman: a naive key that groups on the row's own (corrupted) "T_hs"
    # field directly (never falling back to T_hs_requested) DOES collide on
    # this same synthetic data.
    naive_groups = {}
    for r in synth:
        naive_groups.setdefault((r["height_nm"], r["orientation"], r["screening_fraction"], r["regime"],
                                  r["T_hs"], r["field_polarity"], r["bias_mode"]), []).append(r)
    naive_collides = any(len({rr["V_j"] for rr in grp}) < len(grp) for grp in naive_groups.values())
    checks.append(("a naive T_hs-keyed grouping on the same synthetic data DOES collide "
                    "(confirms the reproduced bug is real, not a strawman)", naive_collides))


def check_failure_logging(man, checks):
    """Fix round 2: 'any per-trace exception in derivatives, compatibility
    or plotting must be caught, logged in the manifest with the trace id
    and reason, and never abort the sweep'. Each list is always present,
    possibly empty (the quick grid never triggers the underlying bug)."""
    for key in ("derivative_trace_failures", "compatibility_trace_failures", "plot_trace_failures"):
        val = man.get(key)
        checks.append((f"manifest records {key} as a list", isinstance(val, list)))
        checks.append((f"every entry in {key} carries a trace_id and reason",
                        all(isinstance(x, dict) and "trace_id" in x and "reason" in x for x in val) if isinstance(val, list) else False))


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


BEST_FLUX_RE = re.compile(r"^BEST_PASSING_FLUX family=(?P<family>\S+) value=(?P<value>\S+) row_id=(?P<row_id>\S+)")

def check_best_passing_flux(text, core, checks):
    """Fix round 3, item 4/8: independently recompute the best flux among
    optical_pass rows per family straight from sweep.csv (using the run
    module's own pure eligible/optical_pass functions, never trusting the
    row's own saved optical_pass column) and compare against the
    machine-checkable BEST_PASSING_FLUX line in results.md."""
    lines_found = [ln for ln in text.splitlines() if ln.startswith("BEST_PASSING_FLUX")]
    checks.append(("results.md carries a BEST_PASSING_FLUX line per family", len(lines_found) >= 1))
    for fam in ("c_plane", "a_plane"):
        m = next((BEST_FLUX_RE.match(ln) for ln in lines_found if f"family={fam} " in ln), None)
        cand = []
        for r in core:
            if r["orientation"] != fam: continue
            valid = r["valid"] == "True"; g2v = f(r["g2"]); fluxv = f(r["signal_flux_s"])
            one_pair = r.get("one_pair_valid") == "True"
            elig = s.eligible(valid, g2v, fluxv)
            opt = s.optical_pass(elig, g2v, r["regime"], one_pair)
            if opt and isfin(fluxv): cand.append((fluxv, r["row_id"]))
        if not cand:
            checks.append((f"BEST_PASSING_FLUX family={fam} matches independent recomputation (no passing rows)",
                            m is not None and m["value"] == "none"))
            continue
        best_val, best_rid = max(cand, key=lambda z: z[0])
        # results.md prints value with %.6g -- compare by row_id (exact) and
        # value with a tolerance wide enough for that 6-sig-fig rounding.
        ok = (m is not None and m["family"] == fam and m["row_id"] == best_rid
              and close_nan_safe(m["value"], best_val, rtol=1e-4, atol=1e-9))
        checks.append((f"BEST_PASSING_FLUX family={fam} row_id/value match independent max-over-optical_pass-rows "
                        f"recomputation (expected row {best_rid}, value {best_val:.6g})", ok))


def check_distinct_hypotheses(text, comp, checks):
    """Fix round 3, item 3/8: independently collapse screening_compatibility.
    csv's raw fit rows to DISTINCT (orientation,height,polarity,screening)
    hypotheses (regime/T_hs are replicates, not independent hypotheses) and
    compare the count against the 'N raw ... reduce to M DISTINCT
    hypotheses' line in results.md."""
    m = re.search(r"(?P<raw>\d+) raw compatibility-fit rows reduce to (?P<distinct>\d+) DISTINCT hypotheses", text)
    checks.append(("results.md states the raw->distinct compatibility hypothesis count", m is not None))
    if m is None or not comp: return
    seen = set()
    for r in comp:
        try:
            gk = dict(json.loads(r["group"]))
        except (KeyError, ValueError, json.JSONDecodeError):
            continue
        seen.add((gk.get("orientation"), gk.get("height_nm"), gk.get("field_polarity"), r.get("screening")))
    checks.append(("results.md's raw compatibility-row count matches screening_compatibility.csv",
                    int(m["raw"]) == len(comp)))
    checks.append((f"results.md's DISTINCT hypothesis count ({m['distinct']}) matches an independent "
                    f"(orientation,height,polarity,screening) collapse of screening_compatibility.csv "
                    f"(expected {len(seen)})", int(m["distinct"]) == len(seen)))


def check_nonpolar_slope_in_window(comp, man, checks):
    """Fix round 3, item 8: every nonpolar screening_unidentifiable row's
    fitted slope actually lies inside the run's own slope window (the
    override in build_compatibility only applies when the window already
    accepted the shared curve)."""
    si = man.get("slope_interval", {})
    lo, hi = (si.get("range_meV_per_V") or [None, None])[:2]
    nonpolar_unident = [r for r in comp if r.get("identification_status") == "screening_unidentifiable"
                         and ("a_plane" in r.get("group", "") or "m_plane" in r.get("group", ""))]
    if lo is None or not nonpolar_unident:
        checks.append(("nonpolar screening_unidentifiable rows' fitted slope lies within the run's slope window "
                        "(no such rows this run, or window unavailable -- vacuously true)", True))
        return
    ok = all(lo - 1e-6 <= f(r.get("fitted_slope_meV_per_V")) <= hi + 1e-6 for r in nonpolar_unident)
    checks.append((f"every nonpolar screening_unidentifiable row's fitted slope lies within [{lo},{hi}] meV/V "
                    "(the override only applies when the window already accepted the shared curve)", ok))


def check_sensitivity_baselines_pass(sens, core, checks):
    """Fix round 3, item 5/8: the two new passing baselines actually
    produce at least one optical-pass row (independently recomputed via
    s.eligible/s.optical_pass on the sensitivities.csv rows themselves,
    never trusting the saved optical_pass column) -- checked against
    sensitivities.csv directly so it holds in --quick too, where the
    reduced core grid (H in {1,7}, R in {5,30}) does not sample the
    baseline's own H=3/R=10 reference point."""
    baselines = {r.get("sensitivity_baseline") for r in sens}
    checks.append(("sensitivities.csv carries the REF_pass_c_plane_screened_SET baseline",
                    "REF_pass_c_plane_screened_SET" in baselines))
    checks.append(("sensitivities.csv carries the REF_pass_a_plane_SET baseline",
                    "REF_pass_a_plane_SET" in baselines))

    def _baseline_passes(label):
        rs = [r for r in sens if r.get("sensitivity_baseline") == label and r.get("row_kind") == "sensitivity"]
        for r in rs:
            valid = r.get("valid") == "True"; g2v = f(r.get("g2")); fluxv = f(r.get("signal_flux_s"))
            one_pair = r.get("one_pair_valid") == "True"
            elig = s.eligible(valid, g2v, fluxv)
            if s.optical_pass(elig, g2v, r.get("regime"), one_pair): return True
        return False
    checks.append(("REF_pass_c_plane_screened_SET yields at least one independently-recomputed optical_pass row",
                    _baseline_passes("REF_pass_c_plane_screened_SET")))
    checks.append(("REF_pass_a_plane_SET yields at least one independently-recomputed optical_pass row",
                    _baseline_passes("REF_pass_a_plane_SET")))


def check_panel_grouping(man, checks):
    """Fix round 3, item 6/8: no Stark figure's panel_index splits a
    screening triplet or polarity pair -- every contract entry sharing the
    same subtitle (the 'orientation h=..nm regime T=..K' super-group
    prefix before ' | ' in its label) must carry the SAME panel_index."""
    mapping = man.get("plot_row_mapping", {})
    stark_figs = [n for n in mapping if n.startswith("stark_")]
    checks.append(("at least one Stark figure is present to check panel grouping on", len(stark_figs) > 0))
    all_ok = True
    for name in stark_figs:
        by_subtitle = {}
        for label, entry in mapping[name].items():
            subtitle = label.split(" | ", 1)[0]
            by_subtitle.setdefault(subtitle, set()).add(entry.get("panel_index"))
        if any(len(idxs) > 1 for idxs in by_subtitle.values()):
            all_ok = False
    checks.append(("every Stark figure's screening/polarity lines within one trace-identity super-group "
                    "share a single panel_index (no split triplets/pairs)", all_ok))


def check_invalid_completeness(text, man, checks):
    """Fix round 3, item 7/8: results.md's invalid-rows accounting matches
    the manifest's own invalid_counts_by_kind/invalid_total exactly (all
    kinds, not just core)."""
    invalid_total = man.get("invalid_total", 0)
    checks.append((f"results.md states invalid_total={invalid_total} matching the manifest",
                    f"invalid_total={invalid_total}" in text))
    by_kind = man.get("invalid_counts_by_kind", {})
    checks.append(("results.md's invalid-rows section names every kind in the manifest's invalid_counts_by_kind",
                    all(f"{k}={v}" in text for k, v in by_kind.items())))


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
    check_derivative_trace_grouping(checks, detail)
    check_failure_logging(man, checks)
    check_convergence(man, checks)
    check_plots(out, man, core, geo, bias, checks)
    check_panel_grouping(man, checks)
    check_best_passing_flux(text, core, checks)
    check_distinct_hypotheses(text, comp, checks)
    check_nonpolar_slope_in_window(comp, man, checks)
    check_sensitivity_baselines_pass(sens, core, checks)
    check_invalid_completeness(text, man, checks)
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
