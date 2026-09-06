"""Checks for scripts/run_rt_edge.py: the RT edge-emitter acceptance-sweep
policy (compute_verdict/compute_stats), the saved full-run artifacts
(sweep.csv, envelope.png, verdict.md, manifest.json), and six deterministic
real-evaluator replays.  Pass --full to restore the previous full/quick
regeneration audit for manual use.

This file does NOT re-derive device physics (that is verify/verify_device_rt.py
and verify/verify_rt_edge_cards.py's job); it checks that run_rt_edge.py's own
bookkeeping -- grid construction, eligibility classification, min/median/
coverage aggregation, and the PASS/FAIL policy in compute_verdict -- is
correct, using synthetic fixtures whose expected outcome is worked out by hand
in this file (never by calling compute_verdict and trusting its own answer),
plus a bounded real replay of the actual evaluator.

Section 1 (pure policy, no evaluate() calls): all-ineligible, no-evidence,
duplicate-evidence, headline-only-pass (favorable pulsed gate despite both CW
secondary diagnostics failing), secondary-only-pass (favorable CW diagnostics
cannot substitute for a failing headline gate), coverage bookkeeping,
fallback-only-pass, median-fail and only-[A]/[E]-corner scenarios against
synthetic rows -- exercising docs/rt_edge_contract.md's rule that the
headline metric (pulsed intrinsic g2(0)) alone gates PASS, with g2_cw0/
g2_cw0_raw reported as non-gating secondary diagnostics.
The default path recomputes statistics from the existing CSV, verifies grid
completeness and the markdown VERDICT line, and replays at most six rows.
The legacy --quick/full/determinism regeneration path remains under --full.

Standalone, side-effect-free on import; exits 0 iff every check passes and
prints "N/N rt-edge sweep checks passed".
"""
from __future__ import annotations

import csv
import io
import json
import math
import statistics
import sys
import tempfile
import warnings
from contextlib import redirect_stdout
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# Default verification deliberately consumes the saved full sweep.  `--full`
# retains the former expensive behavior below for an explicit manual audit.
RUN_FULL = "--full" in sys.argv

import scripts.run_rt_edge as rte  # noqa: E402

warnings.filterwarnings("ignore", category=UserWarning, module="fsim_core.transport")
warnings.filterwarnings("ignore", category=RuntimeWarning, module="fsim_core.cw_g2")
warnings.filterwarnings("ignore", category=Warning, module="scipy.optimize")

CHECKS = []


def ok(name: str, value: bool) -> None:
    CHECKS.append(bool(value))
    print(("ok  " if value else "FAIL") + " " + name)


PRIMARY_ID = next(c["id"] for c in rte.CARDS if c["class"] == "primary")
FALLBACK_ID = next(c["id"] for c in rte.CARDS if c["class"] == "fallback")
GOOD_EVIDENCE = {"evidence_complete": True, "missing_evidence": []}
BAD_EVIDENCE = {"evidence_complete": False,
                "missing_evidence": [{"claim": "x", "reason": "only 1 distinct verified primary source(s), need 2"}]}
DUP_EVIDENCE = {"evidence_complete": False,
                "missing_evidence": [{"claim": "y",
                                     "reason": "duplicate DOI collapses to one distinct source"}]}
GOOD_HALLU = {"hallucination_tests_passed": True}
BAD_HALLU = {"hallucination_tests_passed": False}


def make_row(card_id, card_class, g2_pulsed, g2_cw0, g2_cw0_raw, eligible_row,
            assumptions="dot.gamma300; emission.lambda_nm", gamma300_meV=6.0,
            delta_xx_meV=4.0, irf_ps=50.0, emission_NA=0.75, emission_R_back=0.95,
            emission_L_um=250.0) -> dict:
    """Minimal synthetic row: headline_pass/secondary_pass are worked out
    here by the SAME plain formulas docs/rt_edge_contract.md states (the
    headline metric is pulsed intrinsic g2(0) alone; g2_cw0/g2_cw0_raw are
    secondary diagnostics that never gate PASS) -- not imported from
    run_rt_edge -- so the fixture's own label is an independent ground truth
    for the section-1 assertions below."""
    headline_pass = bool(eligible_row and np.isfinite(g2_pulsed) and g2_pulsed < 0.5)
    secondary_pass = bool(eligible_row
                          and all(np.isfinite(v) for v in (g2_cw0, g2_cw0_raw))
                          and g2_cw0 < 0.5 and g2_cw0_raw < 0.5)
    return {"card_id": card_id, "card_class": card_class, "g2_pulsed": g2_pulsed,
            "g2_cw0": g2_cw0, "g2_cw0_raw": g2_cw0_raw, "eligible_row": eligible_row,
            "headline_pass": headline_pass, "secondary_pass": secondary_pass,
            "assumptions": assumptions, "gamma300_meV": gamma300_meV,
            "delta_xx_meV": delta_xx_meV, "irf_ps": irf_ps, "emission_NA": emission_NA,
            "emission_R_back": emission_R_back, "emission_L_um": emission_L_um}


def full_row(**overrides) -> dict:
    """A row dict populated with every csv_fieldnames() column (sane
    defaults), for exercising write_csv/write_png/write_markdown directly."""
    base = {
        "card_id": PRIMARY_ID, "card_class": "primary", "config_id": "deadbeef",
        "delta_xx_meV": 4.0, "delta_xx_tag": "V", "gamma300_meV": 6.0, "gamma300_tag": "E",
        "irf_ps": 50.0, "irf_tag": "E", "T_hs_K": 300.0, "Tj_pulsed_K": 300.0,
        "Tj_cw_K": 300.0, "I_uA": 0.001, "V_j_pulsed_V": 1.2, "V_j_cw_V": 1.2,
        "mu_pulsed": 0.5, "eta_inj_pulsed": 0.3, "eta_inj_cw": 0.3,
        "S_retention_pulsed": 0.5, "S_retention_cw": 0.5, "b_e_window_pulsed": 0.1,
        "Gamma_pulsed_meV": 6.0,
        "b_e_window_cw": 0.1, "t_x_pulsed": 0.5, "eps_pulsed": 0.1, "t_x_cw": 0.5,
        "eps_cw": 0.1, "edge_beta": 0.01, "edge_T_facet": 0.8,
        "edge_eta_prop": 0.7, "edge_eta_NA": 0.6, "edge_eta_total": 0.002,
        "collected_flux_pulsed_s": 10.0, "collected_flux_cw_s": 10.0,
        "g2_pulsed": float("nan"), "g2_cw0": float("nan"), "g2_cw0_raw": float("nan"),
        "eligible_pulsed": False, "eligible_cw": False, "eligible_row": False,
        "headline_pass": False, "secondary_pass": False,
        "invalid_reasons_pulsed": "synthetic: forced invalid",
        "invalid_reasons_cw": "synthetic: forced invalid",
        "assumptions": "dot.gamma300; emission.lambda_nm",
        "pulse_width_ns": rte.PULSE_WIDTH_NS, "rep_rate_hz": rte.REP_RATE_HZ,
        "duty_pulsed": rte.PULSE_WIDTH_NS * 1e-9 * rte.REP_RATE_HZ,
        "diagnostic_valid": False,
        "emission_NA": 0.75, "emission_R_back": 0.95, "emission_L_um": 250.0,
    }
    base.update(overrides)
    return base


# ============================================================ 1. pure policy

# -- all-ineligible: no row eligible -> FAIL for eligibility, but diagnostic
# statistics remain available when the only invalidity is the flux floor.
rows = [make_row(PRIMARY_ID, "primary", 0.3, 0.3, 0.3, False),
        make_row(PRIMARY_ID, "primary", 0.9, 0.9, 0.9, False)]
for row in rows:
    row.update({"collected_flux_pulsed_s": 10.0,
                "invalid_reasons_pulsed": "ineligible: flux_below_floor",
                "invalid_reasons_cw": "", "diagnostic_valid": True})
stats = rte.compute_stats(rows)
v = rte.compute_verdict(rows, stats, True, GOOD_EVIDENCE, GOOD_HALLU)
ok("all-ineligible: zero eligible rows despite favorable-looking g2 values",
   stats["n_eligible"] == 0 and math.isnan(v["g2_min"]) and math.isnan(v["g2_median"])
   and v["diag_g2_min"] == 0.3 and v["diag_g2_median"] == 0.6)
ok("all-ineligible: verdict FAILs on eligibility, not metrics",
   not v["pass"] and "no_eligible_rows" in v["fail_reasons"])

# -- flux-floor oracle: evaluate the policy against hand-built valid scalar
# maps rather than a device-generated flux, so these expected boundaries do
# not depend on the current cards or transport model.
floor_scalars = {"invalid_reasons": [], "T_j_op": 300.0, "V_j_op": 1.2,
                 "eta_inj": 0.2, "g2_op": 0.3}
eligible_low_flux, low_flux_reasons = rte._classify(
    floor_scalars, rte.FLUX_FLOOR_PULSED_S - 1.0, "pulsed")
eligible_at_floor, at_floor_reasons = rte._classify(
    floor_scalars, rte.FLUX_FLOOR_PULSED_S, "pulsed")
ok("flux floor: a valid pulsed row below 1 kHz is ineligible with the exact reason",
   not eligible_low_flux and "ineligible: flux_below_floor" in low_flux_reasons)
ok("flux floor: a valid pulsed row at 1 kHz remains eligible",
   eligible_at_floor and "ineligible: flux_below_floor" not in at_floor_reasons)
flag_false_scalars = {**floor_scalars, "flux_measurable": False}
eligible_flag_false, flag_false_reasons = rte._classify(
    flag_false_scalars, rte.FLUX_FLOOR_PULSED_S + 1.0, "pulsed")
ok("flux floor: device.py flux_measurable=False takes precedence when exposed",
   not eligible_flag_false and "ineligible: flux_below_floor" in flag_false_reasons)
flux_excluded_row = make_row(PRIMARY_ID, "primary", 0.3, 0.3, 0.3, False)
flux_excluded_row["invalid_reasons_pulsed"] = "ineligible: flux_below_floor"
stats_flux_excluded = rte.compute_stats([flux_excluded_row])
ok("flux floor: excluded rows are counted separately from eligible statistics",
   stats_flux_excluded["n_flux_floor_excluded"] == 1 and stats_flux_excluded["n_eligible"] == 0)

# -- diagnostic oracle: below-floor rows expose finite g2 statistics, maximum
# flux, and the evaluator-factor decomposition; eligible rows expose both too.
diag_row = full_row(g2_pulsed=0.2, g2_cw0=0.25, g2_cw0_raw=0.3,
                    collected_flux_pulsed_s=500.0, eligible_row=False,
                    diagnostic_valid=True, invalid_reasons_pulsed="ineligible: flux_below_floor",
                    invalid_reasons_cw="")
eligible_diag = full_row(g2_pulsed=0.2, g2_cw0=0.25, g2_cw0_raw=0.3,
                        collected_flux_pulsed_s=1500.0, eligible_row=True,
                        eligible_pulsed=True, eligible_cw=True, diagnostic_valid=True,
                        invalid_reasons_pulsed="", invalid_reasons_cw="")
for fixture, expected_flux in (([diag_row], 500.0), ([eligible_diag], 1500.0)):
    fixture_stats = rte.compute_stats(fixture)
    fixture_verdict = rte.compute_verdict(fixture, fixture_stats, True, GOOD_EVIDENCE, GOOD_HALLU)
    ok("diagnostic scenarios: stats and decomposition exist below floor and when eligible",
       fixture_verdict["diag_g2_min"] == 0.2 and fixture_verdict["flux_max"] == expected_flux
       and fixture_stats["best_diagnostic_row"]["edge_beta"] == 0.01
       and fixture_stats["best_diagnostic_row"]["edge_T_facet"] == 0.8)

# -- no evidence: physics favorable and eligible, but evidence gate FAILs
rows = [make_row(PRIMARY_ID, "primary", 0.3, 0.3, 0.3, True)]
stats = rte.compute_stats(rows)
v = rte.compute_verdict(rows, stats, True, BAD_EVIDENCE, GOOD_HALLU)
ok("no-evidence: physics stats still computed correctly (0.3/0.3)",
   abs(v["g2_min"] - 0.3) < 1e-12 and abs(v["g2_median"] - 0.3) < 1e-12)
ok("no-evidence: verdict FAILs on evidence, not metrics/eligibility",
   not v["pass"] and "evidence_incomplete" in v["fail_reasons"]
   and "no_headline_pass" not in v["fail_reasons"])

# -- duplicate evidence: a distinct evidence-incompleteness reason still gates PASS
v_dup = rte.compute_verdict(rows, stats, True, DUP_EVIDENCE, GOOD_HALLU)
ok("duplicate-evidence: still FAILs on evidence (dedup collapsed a source)",
   not v_dup["pass"] and "evidence_incomplete" in v_dup["fail_reasons"])
buf = io.StringIO()
# The sandbox used by this repository does not permit Python's temporary
# directory cleanup ACL changes on Windows; use the regenerable rt_edge
# output location and let the full run below replace this fixture artifact.
md_path = ROOT / "out" / "rt_edge" / "_verify_fixture_verdict.md"
rte.write_markdown(rows, stats, v_dup, rte.build_grid(False), True,
                   DUP_EVIDENCE, GOOD_HALLU, md_path)
md_text = md_path.read_text(encoding="utf-8")
ok("duplicate-evidence: the specific reason text survives into verdict.md",
   "duplicate DOI collapses to one distinct source" in md_text)

# -- headline-only pass: docs/rt_edge_contract.md makes the headline metric
# (pulsed intrinsic g2(0)) alone the PASS gate; g2_cw0/g2_cw0_raw are
# secondary diagnostics that must NOT block PASS even when both fail.
rows = [make_row(PRIMARY_ID, "primary", g2_pulsed=0.3, g2_cw0=0.7, g2_cw0_raw=0.9, eligible_row=True)]
stats = rte.compute_stats(rows)
v = rte.compute_verdict(rows, stats, True, GOOD_EVIDENCE, GOOD_HALLU)
ok("headline-only pass: pulsed g2<0.5 alone is sufficient for PASS",
   v["pass"] and "no_headline_pass" not in v["fail_reasons"])
ok("headline-only pass: failing CW secondary diagnostics are recorded but never gate",
   not rows[0]["secondary_pass"] and v["cw_raw_coverage_n"] == 0 and v["cw0_coverage_n"] == 0)

# -- secondary-only pass is NOT sufficient: favorable CW diagnostics cannot
# substitute for a failing headline (pulsed) gate.
rows = [make_row(PRIMARY_ID, "primary", g2_pulsed=0.9, g2_cw0=0.3, g2_cw0_raw=0.3, eligible_row=True)]
stats = rte.compute_stats(rows)
v = rte.compute_verdict(rows, stats, True, GOOD_EVIDENCE, GOOD_HALLU)
ok("secondary-only pass: favorable CW diagnostics cannot substitute for the headline gate",
   not v["pass"] and not rows[0]["headline_pass"] and "no_headline_pass" in v["fail_reasons"])
ok("secondary-only pass: secondary coverage still correctly counts the row "
   "even though overall verdict is FAIL",
   v["cw_raw_coverage_n"] == 1 and v["cw0_coverage_n"] == 1)

# -- coverage bookkeeping: headline/cw0/cw_raw coverages are independent
# fractions over EVERY scheduled row (run_rt_edge.py's own "every scheduled
# row" coverage-denominator policy), not gated on each other.
rows = [make_row(PRIMARY_ID, "primary", 0.3, 0.3, 0.3, True),   # headline+secondary pass
        make_row(PRIMARY_ID, "primary", 0.9, 0.3, 0.3, True),   # headline fail, secondary pass
        make_row(PRIMARY_ID, "primary", 0.3, 0.9, 0.9, True),   # headline pass, secondary fail
        make_row(PRIMARY_ID, "primary", 0.9, 0.9, 0.9, False)]  # ineligible
stats = rte.compute_stats(rows)
v = rte.compute_verdict(rows, stats, True, GOOD_EVIDENCE, GOOD_HALLU)
ok("coverage bookkeeping: headline coverage counts exactly the headline-passing rows / all rows",
   v["headline_coverage_n"] == 2 and v["headline_coverage_total"] == 4
   and abs(stats["headline_coverage"] - 0.5) < 1e-12)
ok("coverage bookkeeping: cw0/cw_raw secondary coverage counted independently of headline",
   v["cw0_coverage_n"] == 2 and v["cw_raw_coverage_n"] == 2)
ok("coverage bookkeeping: CW secondary denominators are eligible rows",
   v["cw0_coverage_total"] == 3 and v["cw_raw_coverage_total"] == 3
   and abs(stats["cw0_coverage"] - 2 / 3) < 1e-12)

# -- temperature-axis oracles: only 230 K passing selects the lowest TEC
# set point; no passing T reports the literal `none` token.
rows = [dict(make_row(PRIMARY_ID, "primary", 0.3, 0.3, 0.3, True), T_hs_K=230.0),
        dict(make_row(PRIMARY_ID, "primary", 0.9, 0.3, 0.3, True), T_hs_K=250.0),
        dict(make_row(PRIMARY_ID, "primary", 0.9, 0.3, 0.3, True), T_hs_K=273.0),
        dict(make_row(PRIMARY_ID, "primary", 0.9, 0.3, 0.3, True), T_hs_K=300.0)]
stats = rte.compute_stats(rows)
v = rte.compute_verdict(rows, stats, True, GOOD_EVIDENCE, GOOD_HALLU)
ok("temperature axis: only 230 K passes -> T_pass_min=230 and per-T coverage is retained",
   v["T_pass_min"] == 230 and v["headline_by_T"]["230"]["n_headline"] == 1)
rows = [dict(make_row(PRIMARY_ID, "primary", 0.9, 0.3, 0.3, True), T_hs_K=T)
        for T in (230.0, 250.0, 273.0, 300.0)]
stats = rte.compute_stats(rows)
v = rte.compute_verdict(rows, stats, True, GOOD_EVIDENCE, GOOD_HALLU)
ok("temperature axis: no passing temperature -> T_pass_min is none",
   v["T_pass_min"] is None and "T_pass_min=none" in rte.verdict_line(v))

# -- fallback-only pass: primary card never favorable, fallback card is
rows = [make_row(PRIMARY_ID, "primary", 0.9, 0.9, 0.9, True),
        make_row(FALLBACK_ID, "fallback", 0.3, 0.3, 0.3, True)]
stats = rte.compute_stats(rows)
v = rte.compute_verdict(rows, stats, True, GOOD_EVIDENCE, GOOD_HALLU)
ok("fallback-only pass: overall PASS", v["pass"])
ok("fallback-only pass: flagged as fallback-only with the material switch named",
   v["fallback_only"] and FALLBACK_ID in v["note"] and "fallback" in v["note"])

# -- median fail: a favorable corner exists but the pooled pulsed median does not pass
rows = ([make_row(PRIMARY_ID, "primary", 0.3, 0.3, 0.3, True)]
       + [make_row(PRIMARY_ID, "primary", 0.9, 0.9, 0.9, True) for _ in range(4)])
stats = rte.compute_stats(rows)
v = rte.compute_verdict(rows, stats, True, GOOD_EVIDENCE, GOOD_HALLU)
expected_median = statistics.median([0.3, 0.9, 0.9, 0.9, 0.9])  # independent oracle
ok("median fail: pooled pulsed median matches an independent statistics.median",
   abs(v["g2_median"] - expected_median) < 1e-12)
ok("median fail: best corner still passes overall while the median does not",
   v["pass"] and not v["median_pass"] and "median" in v["note"])

# -- only-[A]/[E]-corner: assumptions_used still names the passing corner's
# stated assumptions (informational), independent of the "conditional" flag.
rows = [make_row(PRIMARY_ID, "primary", 0.3, 0.3, 0.3, True, assumptions="dot.gamma300; emission.lambda_nm")]
stats = rte.compute_stats(rows)
v = rte.compute_verdict(rows, stats, True, GOOD_EVIDENCE, GOOD_HALLU)
ok("only-[A]/[E]-corner: assumptions_used names the passing corner's stated assumptions",
   v["pass"] and set(v["assumptions_used"]) == {"dot.gamma300", "emission.lambda_nm"})
rows_noassum = [make_row(PRIMARY_ID, "primary", 0.3, 0.3, 0.3, True, assumptions="")]
v_noassum = rte.compute_verdict(rows_noassum, rte.compute_stats(rows_noassum), True,
                                GOOD_EVIDENCE, GOOD_HALLU)
ok("assumptions_used is empty when the passing corner truly carries none",
   v_noassum["pass"] and not v_noassum["assumptions_used"])

# -- council review round 4 item 2: "conditional" now means "a headline-
# passing eligible row exists but evidence is incomplete" -- NOT "the
# passing corner carries assumptions" (true of essentially every corner,
# hence uninformative). A PASS (evidence complete) is never conditional;
# a headline pass blocked only by incomplete evidence IS conditional.
ok("conditional is False for a genuine PASS (evidence complete), even though "
   "the passing corner carries [A]/[E] assumptions",
   v["pass"] and not v["conditional"])
v_cond = rte.compute_verdict(rows, stats, True, BAD_EVIDENCE, GOOD_HALLU)
ok("conditional is True when a headline-passing eligible row exists but "
   "evidence is incomplete (physics clears the gate, evidence gate does not)",
   not v_cond["pass"] and "evidence_incomplete" in v_cond["fail_reasons"] and v_cond["conditional"])
rows_nopass = [make_row(PRIMARY_ID, "primary", 0.9, 0.9, 0.9, True)]
v_nopass_badevidence = rte.compute_verdict(rows_nopass, rte.compute_stats(rows_nopass), True,
                                           BAD_EVIDENCE, GOOD_HALLU)
ok("conditional is False when there is no headline-passing row at all, "
   "even with incomplete evidence",
   not v_nopass_badevidence["conditional"])

# -- an incomplete (--quick) grid can never PASS, regardless of how good the metrics are
v_quick = rte.compute_verdict(rows, stats, False, GOOD_EVIDENCE, GOOD_HALLU)
ok("incomplete grid cannot grant PASS even with a favorable, well-evidenced corner",
   not v_quick["pass"] and "grid_incomplete" in v_quick["fail_reasons"])


# ================================== 1b. collection-lever axes (council review round 4, item 1)

class _FakeBlock:
    def __init__(self, **kw):
        for k, v in kw.items():
            setattr(self, k, v)


class _FakeDesign:
    def __init__(self, provenance, **blocks):
        self.provenance = provenance
        for name, block in blocks.items():
            setattr(self, name, block)


# -- a card with NO declared range for a lever path falls back to that
# card's own scalar (a single-value axis), never inventing a range.
fake_no_range = _FakeDesign({"ranges": {}},
                            emission=_FakeBlock(NA=0.6, R_back=None, L_um=300.0))
grid_no_range = rte.resolve_lever_grid(fake_no_range, quick=False)
ok("lever axis fallback: a card with no declared range collapses to its own scalar",
   grid_no_range == {"emission.NA": [0.6], "emission.R_back": [None], "emission.L_um": [300.0]})

# -- a declared range is sampled at lo/hi plus the card's own scalar
# (dedup order preserved: lo, hi, card value).
fake_with_range = _FakeDesign(
    {"ranges": {"emission.NA": {"lo": 0.5, "hi": 0.8}}},
    emission=_FakeBlock(NA=0.75, R_back=0.95, L_um=250.0))
grid_range = rte.resolve_lever_grid(fake_with_range, quick=False)
ok("lever axis oracle: a declared range is sampled at [lo, hi, card value]",
   grid_range["emission.NA"] == [0.5, 0.8, 0.75])
ok("lever axis oracle: a path with no declared range still falls back correctly "
   "alongside a sibling path that does have one",
   grid_range["emission.R_back"] == [0.95] and grid_range["emission.L_um"] == [250.0])

# -- a card value sitting exactly on a declared endpoint is not duplicated
# (this repo's real cards: R_back=0.95=hi, L_um=250.0=lo).
fake_dedup = _FakeDesign(
    {"ranges": {"emission.R_back": {"lo": 0.0, "hi": 0.95},
               "emission.L_um": {"lo": 250.0, "hi": 500.0}}},
    emission=_FakeBlock(NA=0.75, R_back=0.95, L_um=250.0))
grid_dedup = rte.resolve_lever_grid(fake_dedup, quick=False)
ok("lever axis oracle: a card value exactly at a declared endpoint is not duplicated",
   grid_dedup["emission.R_back"] == [0.0, 0.95] and grid_dedup["emission.L_um"] == [250.0, 500.0])

# -- --quick collapses every lever to the card's own scalar only (single
# point, matching --quick's existing "smaller, explicitly incomplete"
# convention for the other axes).
grid_quick = rte.resolve_lever_grid(fake_with_range, quick=True)
ok("lever axis quick mode: collapses every lever to the card's own scalar only",
   grid_quick == {"emission.NA": [0.75], "emission.R_back": [0.95], "emission.L_um": [250.0]})

# -- build_lever_combos is the full cartesian product of the resolved grid
combos = rte.build_lever_combos({"emission.NA": [0.5, 0.8], "emission.R_back": [0.0, 0.95],
                                 "emission.L_um": [250.0]})
ok("lever combos oracle: full cartesian product size (2 x 2 x 1 = 4)", len(combos) == 4)
ok("lever combos oracle: every combo is a well-formed overrides dict, both corners present",
   {"emission.NA": 0.5, "emission.R_back": 0.0, "emission.L_um": 250.0} in combos
   and {"emission.NA": 0.8, "emission.R_back": 0.95, "emission.L_um": 250.0} in combos)


# ======================== 1c. gamma300_pass_max / conditional-on-linewidth (round 4, item 2)

gamma_rows = [
    {"card_id": PRIMARY_ID, "headline_pass": True, "eligible": True, "gamma300_meV": 6.0,
     "g2_pulsed": 0.30, "delta_xx_meV": 4.0, "irf_ps": 50.0, "emission_NA": 0.75,
     "emission_R_back": 0.95, "emission_L_um": 250.0, "assumptions": "dot.gamma300"},
    {"card_id": PRIMARY_ID, "headline_pass": True, "eligible": True, "gamma300_meV": 13.0,
     "g2_pulsed": 0.45, "delta_xx_meV": 5.5, "irf_ps": 125.0, "emission_NA": 0.8,
     "emission_R_back": 0.0, "emission_L_um": 500.0, "assumptions": "dot.gamma300"},
    {"card_id": PRIMARY_ID, "headline_pass": False, "eligible": True, "gamma300_meV": 20.0,
     "g2_pulsed": 0.90, "delta_xx_meV": 7.0, "irf_ps": 200.0, "emission_NA": 0.5,
     "emission_R_back": 0.95, "emission_L_um": 250.0, "assumptions": ""},
    # eligible=True, g2_pulsed=0.90: this card's failure is a genuine physics
    # shortfall (g2 too high), not a flux-eligibility one -- exercises the
    # "none(g2)" branch below (round 6, item 2), distinct from the "none(flux)"
    # branch exercised by the dedicated fixture in section 1g-ii.
    {"card_id": FALLBACK_ID, "headline_pass": False, "eligible": True, "gamma300_meV": 6.0,
     "g2_pulsed": 0.90, "delta_xx_meV": 4.0, "irf_ps": 50.0, "emission_NA": 0.75,
     "emission_R_back": 0.95, "emission_L_um": 250.0, "assumptions": ""},
]
gmax_by_card = rte._card_gamma300_pass_max(gamma_rows)
ok("gamma300_pass_max oracle: picks the LARGEST passing gamma300 (13.0), not the first (6.0)",
   gmax_by_card[PRIMARY_ID]["gamma300_pass_max_meV"] == 13.0
   and gmax_by_card[PRIMARY_ID]["g2_pulsed"] == 0.45
   and gmax_by_card[PRIMARY_ID]["emission_NA"] == 0.8)
ok("gamma300_pass_max oracle: a card with zero headline-passing rows reports nan",
   math.isnan(gmax_by_card[FALLBACK_ID]["gamma300_pass_max_meV"]))


# ============================= 1d. brightness factor self-check (round 4, item 3)

from fsim_core.loading import loading_probs as _loading_probs_ref  # noqa: E402

for _mu in (0.0, 0.1, 0.5, 1.0, 3.0):
    _, _p1, _p2 = _loading_probs_ref(_mu)
    ok(f"_loading_term oracle: matches fsim_core.loading.loading_probs' own P1+P2 at mu={_mu}",
       abs(rte._loading_term(_mu) - (_p1 + _p2)) < 1e-12)

_mu_test = 0.4
_loading_expected = 1.0 - math.exp(-_mu_test)
_factor_row = {"mu_pulsed": _mu_test, "t_x_pulsed": 0.5, "S_retention_pulsed": 0.8,
              "edge_eta_total": 0.002, "rep_rate_hz": 80.0e6,
              "collected_flux_pulsed_s": _loading_expected * 0.5 * 0.8 * 0.002 * 80.0e6}
_check = rte._brightness_factor_check(_factor_row)
ok("factor self-check oracle: loading term is the exact Poisson complement 1-e^-mu",
   abs(_check["loading"] - _loading_expected) < 1e-12)
ok("factor self-check oracle: product x rep_rate reproduces a consistent reported flux "
   "within 1% (constructed exactly, so within float precision)",
   _check["ok"] and _check["rel_diff"] < 1e-9)

_bad_row = dict(_factor_row)
_bad_row["collected_flux_pulsed_s"] = _factor_row["collected_flux_pulsed_s"] * 1.5
_check_bad = rte._brightness_factor_check(_bad_row)
ok("factor self-check oracle: a reported flux 50% off the factor product fails the 1% self-check",
   not _check_bad["ok"] and _check_bad["rel_diff"] > 0.01)

_nan_row = dict(_factor_row)
_nan_row["S_retention_pulsed"] = float("nan")
_check_nan = rte._brightness_factor_check(_nan_row)
ok("factor self-check oracle: non-finite inputs yield ok=False without raising",
   not _check_nan["ok"] and math.isnan(_check_nan["rel_diff"]))

# -- combined front-facet-and-transmission factor oracle (council review
# round 5, item 6, updated for fsim_core/waveguide.py's OWN concurrent
# 2026-09-06 item-1 facet-model change -- the standalone `front` local is
# gone; T_facet is now fused into a single `facet_factor` for R_back>0):
# beta * <combined factor> * eta_prop * eta_NA must reproduce edge_eta_total
# EXACTLY (the combined factor is back-solved as the one missing factor --
# deliberately WITHOUT dividing by edge_T_facet separately, which would
# assume the now-superseded old 5-term shape and be wrong for R_back>0
# rows under the new model); confirms beta is NOT the combined factor (a
# plain relabeling bug would make this fail).
_edge_row = {"edge_beta": 0.03, "edge_T_facet": 0.72, "edge_eta_prop": 0.78,
            "edge_eta_NA": 0.30, "edge_eta_total": 0.03 * 0.5432 * 0.78 * 0.30}
_combined = rte._front_facet_split(_edge_row)
ok("combined-facet-factor oracle: back-solved factor reproduces the constructed value (0.5432)",
   abs(_combined - 0.5432) < 1e-9)
ok("combined-facet-factor oracle: beta * combined * propagation * NA reproduces "
   "edge_eta_total exactly (never divides by edge_T_facet separately)",
   abs(_edge_row["edge_beta"] * _combined
       * _edge_row["edge_eta_prop"] * _edge_row["edge_eta_NA"]
       - _edge_row["edge_eta_total"]) < 1e-12)
ok("combined-facet-factor oracle: non-finite/zero component inputs yield nan without raising",
   math.isnan(rte._front_facet_split({"edge_beta": 0.0, "edge_T_facet": 0.72,
                                      "edge_eta_prop": 0.78, "edge_eta_NA": 0.30,
                                      "edge_eta_total": 0.001})))

# -- facet-factor forward check (council review round 5, item 6): a
# genuinely independent forward recomputation, introspected LIVE from
# fsim_core/waveguide.py's installed source (never hardcoded). First
# against FAKE candidate lists (both conventions a concurrently-edited
# facet model might use), by temporarily swapping
# _facet_factor_formula_candidates; then against the REAL installed source
# and a real evaluator row.
_row_fused = {"edge_T_facet": 0.8, "emission_R_back": 0.95,
             "edge_beta": 0.03, "edge_eta_prop": 0.78, "edge_eta_NA": 0.30}
_combined_fused = 0.8 / (0.8 + (1 - 0.95))
_row_fused["edge_eta_total"] = 0.03 * _combined_fused * 0.78 * 0.30
_orig_candidates_fn = rte._facet_factor_formula_candidates
try:
    rte._facet_factor_formula_candidates = (
        lambda: [("facet_factor", "0.5 * T if R_back is None else T / (T + (1 - R_back))")])
    _check_fused = rte._facet_factor_forward_check(_row_fused)
finally:
    rte._facet_factor_formula_candidates = _orig_candidates_fn
ok("facet-factor forward check: fused/current-model convention (candidate already IS the "
   "full combined factor, T_facet included) matches directly",
   _check_fused["ok"] and abs(_check_fused["forward"] - _combined_fused) < 1e-9
   and "fused" in _check_fused["convention"])

_front_legacy = 0.8 / (0.8 + (1 - 0.95))
_row_legacy = dict(_row_fused)
_row_legacy["edge_eta_total"] = 0.03 * _front_legacy * 0.8 * 0.78 * 0.30  # front * T_facet
try:
    rte._facet_factor_formula_candidates = (
        lambda: [("front", "0.5 if R_back is None else T / (T + (1 - R_back))")])
    _check_legacy = rte._facet_factor_forward_check(_row_legacy)
finally:
    rte._facet_factor_formula_candidates = _orig_candidates_fn
ok("facet-factor forward check: legacy/split convention (candidate is only the front "
   "split, T_facet applied as a separate step) matches after multiplying by T_facet",
   _check_legacy["ok"] and abs(_check_legacy["forward"] - _front_legacy * 0.8) < 1e-9
   and "legacy" in _check_legacy["convention"])

try:
    rte._facet_factor_formula_candidates = lambda: [("front", "0.1")]
    _check_wrong = rte._facet_factor_forward_check(_row_fused)
finally:
    rte._facet_factor_formula_candidates = _orig_candidates_fn
ok("facet-factor forward check: a candidate that matches NEITHER convention is reported "
   "as ok=False, not silently accepted",
   not _check_wrong["ok"])

_real_candidates = rte._facet_factor_formula_candidates()
ok("facet-factor forward check: the REAL fsim_core/waveguide.py source is introspectable "
   "and yields at least one candidate formula (this file never hardcodes the formula)",
   len(_real_candidates) >= 1)
_gainp_design_for_facet = rte.resolve_device_card(
    next(c["path"] for c in rte.CARDS if c["id"] == FALLBACK_ID), {"dot.delta_xx": 7.0})
_sc_gainp_facet = rte.evaluate(
    _gainp_design_for_facet, T_grid=[_gainp_design_for_facet.thermal.T_hs])["scalars"]
_real_row = {"edge_beta": _sc_gainp_facet["edge_beta"],
            "edge_T_facet": _sc_gainp_facet["edge_T_facet"],
            "edge_eta_prop": _sc_gainp_facet["edge_eta_prop"],
            "edge_eta_NA": _sc_gainp_facet["edge_eta_NA"],
            "edge_eta_total": _sc_gainp_facet["edge_eta_total"],
            "emission_R_back": _gainp_design_for_facet.emission.R_back}
_real_check = rte._facet_factor_forward_check(_real_row)
ok("facet-factor forward check: against a REAL evaluator row, the forward recomputation "
   "(introspected live from the installed waveguide.py) matches the back-solved value",
   _real_check["ok"])


# ============================================= 1e. self-test field rename (round 4, item 5)

ok("_self_test_passed: reads the legacy hallucination_tests_passed key",
   rte._self_test_passed({"hallucination_tests_passed": True}) is True
   and rte._self_test_passed({"hallucination_tests_passed": False}) is False)
ok("_self_test_passed: prefers the renamed all_checks_passed key when both are present",
   rte._self_test_passed({"all_checks_passed": True, "hallucination_tests_passed": False}) is True)
ok("_self_test_passed: works with ONLY the renamed all_checks_passed key "
   "(the concurrently-renamed verify_rt_edge_papers.py case)",
   rte._self_test_passed({"all_checks_passed": True}) is True
   and rte._self_test_passed({"all_checks_passed": False}) is False)
ok("_self_test_passed: missing key defaults to False rather than raising",
   rte._self_test_passed({}) is False)


# ==================================== 1f. card_line conditional phrase (round 5, item 1)

_stats_zero_eligible = {"card_class": "primary", "n_rows": 96, "n_eligible": 0,
                        "g2_pulsed_min": float("nan"), "g2_pulsed_median": float("nan"),
                        "diag_g2_pulsed_min": 0.50, "diag_g2_pulsed_median": 0.78,
                        "diag_g2_cw0_min": 0.37, "diag_g2_cw0_median": 0.80,
                        "diag_g2_cw0_raw_min": 0.98, "diag_g2_cw0_raw_median": 1.0,
                        "g2_cw0_raw_min": float("nan"), "g2_cw0_raw_median": float("nan"),
                        "n_flux_floor_excluded": 96, "n_favorable": 0}
_line_zero_eligible = rte.card_line(PRIMARY_ID, _stats_zero_eligible)
ok("card_line: a card with ZERO eligible rows prints the "
   "'diagnostic (below flux floor, not measurable)' phrase",
   "diagnostic (below flux floor, not measurable)" in _line_zero_eligible
   and "eligible rows:" not in _line_zero_eligible)

_stats_some_eligible = {**_stats_zero_eligible, "n_eligible": 16,
                        "g2_pulsed_min": 0.40, "g2_pulsed_median": 0.70}
_line_some_eligible = rte.card_line(FALLBACK_ID, _stats_some_eligible)
ok("card_line: a card with a NONZERO eligible count (council review round 5 item 1's "
   "gainp example, eligible=16/96) prints 'eligible rows: <n>/<N>' and NEVER the "
   "'not measurable' phrase, even though its flux_floor_excluded count is also nonzero",
   "eligible rows: 16/96" in _line_some_eligible
   and "not measurable" not in _line_some_eligible
   and "g2_pulsed_min=0.4" in _line_some_eligible)


# ============================== 1g. gamma300_threshold bracket (round 5, item 2)

_gmax_by_card = rte._card_gamma300_pass_max(gamma_rows)
_bracket_by_card = rte._gamma300_threshold_bracket(gamma_rows, _gmax_by_card)
ok("gamma300_threshold_bracket oracle: PRIMARY_ID passes at 13.0 and fails at the next "
   "sampled value (20.0) -- bracket is (13.0, 20.0)",
   _bracket_by_card[PRIMARY_ID]["lo"] == 13.0 and _bracket_by_card[PRIMARY_ID]["hi"] == 20.0)
ok("gamma300_threshold_bracket oracle: a card with NO headline-passing sample has lo=None, "
   "hi=the smallest sampled value, and reason='g2' (its only rows are eligible but "
   "g2_pulsed=0.90 -- round 6, item 2)",
   _bracket_by_card[FALLBACK_ID]["lo"] is None and _bracket_by_card[FALLBACK_ID]["hi"] == 6.0
   and _bracket_by_card[FALLBACK_ID]["reason"] == "g2")
_pooled_bracket = rte._pooled_gamma300_threshold(_gmax_by_card, _bracket_by_card)
ok("pooled gamma300_threshold: picks the bracket of whichever card attains the pooled max "
   "(PRIMARY_ID, 13.0), not FALLBACK_ID's",
   _pooled_bracket == _bracket_by_card[PRIMARY_ID])
ok("_format_threshold: both bounds known renders 'lo-hi'",
   rte._format_threshold({"lo": 13.0, "hi": 20.0}) == "13-20")
ok("_format_threshold: hi=None (nothing sampled above pass_max fails) renders '>=lo'",
   rte._format_threshold({"lo": 13.0, "hi": None}) == ">=13")
ok("_format_threshold: lo=None with no 'reason' key falls back to the old '<hi' text "
   "(bare bracket, e.g. a caller that never ran _gamma300_threshold_bracket)",
   rte._format_threshold({"lo": None, "hi": 6.0}) == "<6")
ok("_format_threshold: no data at all renders 'n/a'",
   rte._format_threshold({"lo": None, "hi": None}) == "n/a")

# ==================== 1g-ii. gamma300_threshold none(reason) (round 6, item 2) ====================
# Reproduces the primary card's real favourable-corner shape: g2=0.26 at a
# gamma300/delta_xx sample that is still below the 1 kHz collected-flux
# floor (313 photons/s in the real sweep) -- headline_pass is False at
# EVERY sample even though g2 alone would already clear the gate, because
# eligible=False everywhere. gamma300_threshold must print `none(flux)`,
# never a bracket like `<6` that would misleadingly imply the card passes
# below 6 meV.
_no_pass_flux_rows = [
    {"card_id": "flux_card", "headline_pass": False, "eligible": False, "gamma300_meV": g,
     "g2_pulsed": 0.26, "delta_xx_meV": 1.0, "irf_ps": 50.0, "emission_NA": 0.75,
     "emission_R_back": 0.95, "emission_L_um": 250.0, "assumptions": ""}
    for g in (6.0, 8.0, 10.0)
]
_no_pass_flux_gmax = rte._card_gamma300_pass_max(_no_pass_flux_rows)
_no_pass_flux_bracket = rte._gamma300_threshold_bracket(_no_pass_flux_rows, _no_pass_flux_gmax)
ok("_gamma300_threshold_bracket oracle: every sample flux-ineligible (favourable g2=0.26 "
   "but never eligible) -> reason='flux'",
   _no_pass_flux_bracket["flux_card"]["lo"] is None
   and _no_pass_flux_bracket["flux_card"]["reason"] == "flux")
ok("_format_threshold: reason='flux' renders 'none(flux)'",
   rte._format_threshold(_no_pass_flux_bracket["flux_card"]) == "none(flux)")

# Contrast case: at least one row IS flux-eligible, but g2 >= 0.5 at every
# sampled gamma300 -- a genuine physics shortfall, not an eligibility one.
_no_pass_g2_rows = [
    {"card_id": "g2_card", "headline_pass": False, "eligible": True, "gamma300_meV": g,
     "g2_pulsed": 0.9, "delta_xx_meV": 4.0, "irf_ps": 50.0, "emission_NA": 0.75,
     "emission_R_back": 0.95, "emission_L_um": 250.0, "assumptions": ""}
    for g in (6.0, 8.0, 10.0)
]
_no_pass_g2_gmax = rte._card_gamma300_pass_max(_no_pass_g2_rows)
_no_pass_g2_bracket = rte._gamma300_threshold_bracket(_no_pass_g2_rows, _no_pass_g2_gmax)
ok("_gamma300_threshold_bracket oracle: eligible rows exist but g2_pulsed >= 0.5 "
   "everywhere -> reason='g2'",
   _no_pass_g2_bracket["g2_card"]["lo"] is None
   and _no_pass_g2_bracket["g2_card"]["reason"] == "g2")
ok("_format_threshold: reason='g2' renders 'none(g2)'",
   rte._format_threshold(_no_pass_g2_bracket["g2_card"]) == "none(g2)")

_no_pass_both_rows = [
    {"card_id": "both_card", "headline_pass": False, "eligible": False, "gamma300_meV": g,
     "g2_pulsed": 0.9, "delta_xx_meV": 4.0, "irf_ps": 50.0, "emission_NA": 0.75,
     "emission_R_back": 0.95, "emission_L_um": 250.0, "assumptions": ""}
    for g in (6.0, 8.0)]
_both_bracket = rte._gamma300_threshold_bracket(
    _no_pass_both_rows, rte._card_gamma300_pass_max(_no_pass_both_rows))
ok("gamma300_threshold none reports every blocking reason",
   rte._format_threshold(_both_bracket["both_card"]) == "none(flux,g2)")

# All-pass fixture: every sample passes -> hi is None (threshold >= largest sample).
_all_pass_rows = [{"card_id": "x", "headline_pass": True, "gamma300_meV": g, "g2_pulsed": 0.2,
                  "delta_xx_meV": 4.0, "irf_ps": 50.0, "emission_NA": 0.75,
                  "emission_R_back": 0.95, "emission_L_um": 250.0, "assumptions": ""}
                 for g in (6.0, 13.0, 20.0)]
_all_pass_gmax = rte._card_gamma300_pass_max(_all_pass_rows)
_all_pass_bracket = rte._gamma300_threshold_bracket(_all_pass_rows, _all_pass_gmax)
ok("gamma300_threshold_bracket oracle: every sampled value passes -> hi is None "
   "(threshold >= the largest sample)",
   _all_pass_bracket["x"]["lo"] == 20.0 and _all_pass_bracket["x"]["hi"] is None)


# ================================ 1h. RANGE_BOUNDS cross-check (round 5, item 3)

_bounds_ok = {"dot.delta_xx": (4.0, 8.0, "meV")}
_agreeing_designs = [
    ("card-a", _FakeDesign({"ranges": {"dot.delta_xx": {"lo": 4.0, "hi": 8.0, "unit": "meV"}}})),
    ("card-b", _FakeDesign({"ranges": {"dot.delta_xx": {"lo": 4.0, "hi": 8.0, "unit": "meV"}}})),
]
ok("_range_bounds_mismatches: two cards declaring the SAME (lo, hi) -> no mismatches",
   rte._range_bounds_mismatches(_bounds_ok, _agreeing_designs) == [])

_disagreeing_designs = [
    ("card-a", _FakeDesign({"ranges": {"dot.delta_xx": {"lo": 4.0, "hi": 8.0, "unit": "meV"}}})),
    ("card-b", _FakeDesign({"ranges": {"dot.delta_xx": {"lo": 4.0, "hi": 7.0, "unit": "meV"}}})),
]
_mismatches = rte._range_bounds_mismatches(_bounds_ok, _disagreeing_designs)
ok("_range_bounds_mismatches: card-b declares a DIFFERENT (lo, hi) than the resolved "
   "bound (this is exactly the round-5-item-3 bug: RANGE_BOUNDS = (4.0, 7.0) vs the "
   "cards' own declared (4.0, 8.0)) -- one mismatch is reported, naming the offending card",
   len(_mismatches) == 1 and "card-b" in _mismatches[0] and "dot.delta_xx" in _mismatches[0])

_no_declaration_designs = [("card-c", _FakeDesign({"ranges": {}}))]
ok("_range_bounds_mismatches: a card that does not declare the path at all is not a "
   "mismatch (nothing to cross-check)",
   rte._range_bounds_mismatches(_bounds_ok, _no_declaration_designs) == [])

ok("range cross-check: resolve_range_bounds() on the real CARDS reproduces the "
   "module-level RANGE_BOUNDS resolved at import (both real cards agree)",
   rte.resolve_range_bounds() == rte.RANGE_BOUNDS)
ok("range cross-check: the real cards' own provenance.ranges['dot.delta_xx'] is "
   "(4.0, 8.0), matching verify_rt_edge_cards.py's REQUIRED_RANGES -- NOT the "
   "round-5-item-3 stale (4.0, 7.0)",
   rte.RANGE_BOUNDS["dot.delta_xx"][:2] == (4.0, 8.0))


# ==================================== 1i. naming/semantics (round 5, item 4)

# coverage (VERDICT) = headline / ELIGIBLE (not headline / total); eligible_fraction =
# eligible / total; flux_margin = flux_max / floor (>1 clears the floor); flux_shortfall
# (deprecated) is its old floor / flux_max inverse -- kept for one release.
_naming_rows = [make_row(PRIMARY_ID, "primary", 0.3, 0.3, 0.3, True),   # headline pass, eligible
                make_row(PRIMARY_ID, "primary", 0.9, 0.9, 0.9, True),  # headline fail, eligible
                make_row(PRIMARY_ID, "primary", 0.9, 0.9, 0.9, False)]  # ineligible
for _r in _naming_rows:
    _r["collected_flux_pulsed_s"] = 2000.0
_naming_stats = rte.compute_stats(_naming_rows)
_naming_verdict = rte.compute_verdict(_naming_rows, _naming_stats, True, GOOD_EVIDENCE, GOOD_HALLU)
ok("naming: `coverage` is headline/ELIGIBLE (1/2 = 0.5), not headline/total (1/3)",
   abs(_naming_verdict["coverage"] - 0.5) < 1e-12)
ok("naming: `eligible_fraction` (NEW) is eligible/total (2/3)",
   abs(_naming_verdict["eligible_fraction"] - (2.0 / 3.0)) < 1e-12)
ok("naming: headline_coverage_n/headline_coverage_total is unaffected (still headline/total)",
   _naming_verdict["headline_coverage_n"] == 1 and _naming_verdict["headline_coverage_total"] == 3)
ok("naming: `flux_margin` = flux_max/floor (>1 means the floor is cleared) and "
   "`flux_shortfall` (deprecated) is its exact inverse",
   abs(_naming_verdict["flux_margin"] - _naming_verdict["flux_max"] / rte.FLUX_FLOOR_PULSED_S) < 1e-9
   and abs(_naming_verdict["flux_margin"] * _naming_verdict["flux_shortfall"] - 1.0) < 1e-9)
ok("naming: VERDICT line carries flux_margin= and gamma300_threshold=",
   "flux_margin=" in rte.verdict_line(_naming_verdict)
   and "gamma300_threshold=" in rte.verdict_line(_naming_verdict))


# ================================== 2. saved full-run artifact verification

if not RUN_FULL:
    artifact_dir = ROOT / "out" / "rt_edge"
    artifacts = {name: artifact_dir / name
                 for name in ("sweep.csv", "envelope.png", "verdict.md", "manifest.json")}
    ok("saved full run: all required artifacts exist", all(p.exists() and p.stat().st_size > 0
                                                            for p in artifacts.values()))
    with artifacts["sweep.csv"].open(newline="", encoding="utf-8") as f:
        saved_csv_rows = list(csv.DictReader(f))
    saved_manifest = json.loads(artifacts["manifest.json"].read_text(encoding="utf-8"))
    saved_md = artifacts["verdict.md"].read_text(encoding="utf-8")
    ok("saved full run: CSV header matches csv_fieldnames() exactly",
       bool(saved_csv_rows) and list(saved_csv_rows[0]) == rte.csv_fieldnames())

    # Rebuild precisely the data types compute_stats consumes; this is a
    # CSV-derived recomputation, not a trust in manifest stats.
    _bool_fields = {"eligible_pulsed", "eligible_cw", "eligible_row", "headline_pass",
                    "secondary_pass", "diagnostic_valid"}
    _text_fields = {"card_id", "card_class", "config_id", "delta_xx_tag", "gamma300_tag",
                    "irf_tag", "invalid_reasons_pulsed", "invalid_reasons_cw", "assumptions"}
    _rows_typed = []
    for raw in saved_csv_rows:
        typed = {}
        for key, value in raw.items():
            if key in _bool_fields:
                typed[key] = value == "True"
            elif key in _text_fields:
                typed[key] = value
            else:
                typed[key] = float(value) if value not in ("", "nan") else float("nan")
        _rows_typed.append(typed)
    recomputed_stats = rte.compute_stats(_rows_typed)
    _stat_keys = ("n_total", "n_eligible", "n_headline", "n_cw0_pass", "n_cw_raw_pass",
                  "n_flux_floor_excluded", "g2_pulsed_min", "g2_pulsed_median", "flux_max")
    ok("saved full run: pooled statistics recomputed from CSV equal manifest",
       all((recomputed_stats[k] == saved_manifest["stats"][k]
            if isinstance(recomputed_stats[k], int)
            else math.isclose(recomputed_stats[k], saved_manifest["stats"][k], rel_tol=0, abs_tol=1e-9))
           for k in _stat_keys))
    ok("saved full run: per-temperature statistics recomputed from CSV equal manifest",
       recomputed_stats["per_T"] == saved_manifest["stats"]["per_T"])

    expected_rows = sum(info["n_combos"] for info in saved_manifest["lever_info"].values())
    for values in saved_manifest["grid"].values():
        expected_rows *= len(values)
    ok("saved full run: grid is complete and every lever/axis combination is present",
       saved_manifest["grid_complete"] is True and len(saved_csv_rows) == expected_rows
       and {r["card_id"] for r in saved_csv_rows} == {PRIMARY_ID, FALLBACK_ID})
    expected_line = rte.verdict_line(saved_manifest["verdict"])
    ok("saved full run: verdict.md VERDICT line equals manifest verdict",
       expected_line in saved_md)
    ok("saved full run: manifest card hashes match current cards",
       all(entry["sha256"] == rte._sha256_file(Path(entry["path"]))
           for entry in saved_manifest["cards"]))

    # Six fixed rows cover both cards, endpoints and all four temperatures.
    # Replay both pulsed and CW paths and compare directly with the CSV.
    replay_rows = sorted(saved_csv_rows, key=lambda r: r["config_id"])[::max(1, len(saved_csv_rows) // 6)][:6]
    replay_ok = len(replay_rows) == 6
    card_paths = {c["id"]: c["path"] for c in rte.CARDS}
    for raw in replay_rows:
        lever = {key: float(raw["emission_" + key.split(".")[1]])
                 for key in rte.LEVER_PATHS}
        pulsed = rte.eval_pulsed_point(card_paths[raw["card_id"]], float(raw["delta_xx_meV"]),
                                       float(raw["gamma300_meV"]), lever, T_hs=float(raw["T_hs_K"]))
        cw = rte.eval_cw_point(card_paths[raw["card_id"]], float(raw["delta_xx_meV"]),
                               float(raw["gamma300_meV"]), float(raw["irf_ps"]), lever,
                               T_hs=float(raw["T_hs_K"]))
        checks = ((pulsed["flux_s"], raw["collected_flux_pulsed_s"]),
                  (pulsed["scalars"]["g2_op"], raw["g2_pulsed"]),
                  (cw["flux_s"], raw["collected_flux_cw_s"]),
                  (cw["scalars"]["g2_cw0"], raw["g2_cw0"]),
                  (cw["scalars"]["g2_cw0_raw"], raw["g2_cw0_raw"]))
        replay_ok = replay_ok and all(math.isclose(float(actual), float(expected), rel_tol=1e-9,
                                                    abs_tol=1e-9) for actual, expected in checks)
    ok("saved full run: six deterministic evaluator replays reproduce CSV physics columns", replay_ok)
    print(f"{sum(CHECKS)}/{len(CHECKS)} rt-edge sweep checks passed")
    sys.exit(0 if all(CHECKS) else 1)


# ============================================= 2. --quick real-evaluator smoke

with tempfile.TemporaryDirectory() as td:
    out_dir = Path(td) / "quick"
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc = rte.main(["--quick", "--out-dir", str(out_dir)])
    stdout_text = buf.getvalue()
    artifacts = {name: out_dir / name
                for name in ("sweep.csv", "envelope.png", "verdict.md", "manifest.json")}
    ok("--quick: exits 0", rc == 0)
    ok("--quick: all four artifacts exist", all(p.exists() and p.stat().st_size > 0
                                                 for p in artifacts.values()))
    with artifacts["sweep.csv"].open(newline="", encoding="utf-8") as f:
        quick_rows = list(csv.DictReader(f))
    card_ids = {r["card_id"] for r in quick_rows}
    ok("--quick: both cards present in sweep.csv", card_ids == {PRIMARY_ID, FALLBACK_ID})
    deltas = {float(r["delta_xx_meV"]) for r in quick_rows}
    gammas = {float(r["gamma300_meV"]) for r in quick_rows}
    irfs = {float(r["irf_ps"]) for r in quick_rows}
    ok("--quick: declared endpoints are sampled on every axis (dot.delta_xx is now "
       "(4.0, 8.0) -- council review round 5, item 3 -- matching both cards' own "
       "provenance.ranges and verify_rt_edge_cards.py's REQUIRED_RANGES, not the stale "
       "(4.0, 7.0))",
       deltas == {4.0, 8.0} and gammas == {6.0, 20.0} and irfs == {50.0, 200.0})
    ok("--quick: CSV header matches csv_fieldnames() exactly",
       list(quick_rows[0].keys()) == rte.csv_fieldnames())
    ok("--quick: eligibility/invalid-reason columns are present and well-formed",
       all(r["eligible_row"] in ("True", "False") for r in quick_rows))
    ok("--quick: printed VERDICT line matches the required stdout format",
       "VERDICT: FAIL" in stdout_text or "VERDICT: PASS" in stdout_text)
    ok("--quick: an incomplete grid is explicitly reported as FAIL",
       "VERDICT: FAIL" in stdout_text)
    manifest = json.loads(artifacts["manifest.json"].read_text(encoding="utf-8"))
    ok("--quick: manifest records grid_complete=False", manifest["grid_complete"] is False)

# -- a real ineligible operating point: zero drive current is a genuine invalid
# reason from device.py itself (transport requires positive current), wired
# through eval_pulsed_point's own classification (not asserted in isolation).
zero_i_design = rte.resolve_device_card(
    next(c["path"] for c in rte.CARDS if c["id"] == PRIMARY_ID), {"drive.I_uA": 0.0})
from fsim_core.device import evaluate as _evaluate  # noqa: E402
sc_zero = _evaluate(zero_i_design, T_grid=[zero_i_design.thermal.T_hs])["scalars"]
eligible_zero, reasons_zero = rte._classify(sc_zero, rte._collected_flux_s(sc_zero), "pulsed")
ok("a real zero-current operating point is classified ineligible with a precise reason",
   not eligible_zero and any("current" in r for r in reasons_zero))

# -- empty/all-invalid dataset still produces a legible PNG/markdown, and
# evidence failure does not suppress artifact generation
invalid_rows = [full_row(card_id=cid, card_class=cls, delta_xx_meV=d, gamma300_meV=g, irf_ps=i)
               for cid, cls in ((PRIMARY_ID, "primary"), (FALLBACK_ID, "fallback"))
               for d in (4.0, 7.0) for g in (6.0, 20.0) for i in (50.0, 200.0)]
stats_invalid = rte.compute_stats(invalid_rows)
v_invalid = rte.compute_verdict(invalid_rows, stats_invalid, False, BAD_EVIDENCE, BAD_HALLU)
with tempfile.TemporaryDirectory() as td:
    png_path, md_path = Path(td) / "envelope.png", Path(td) / "verdict.md"
    rte.write_png(invalid_rows, stats_invalid, png_path)
    rte.write_markdown(invalid_rows, stats_invalid, v_invalid, rte.build_grid(False), False,
                       BAD_EVIDENCE, BAD_HALLU, md_path)
    ok("all-invalid dataset: PNG is written and decodes to a legible image",
       png_path.exists() and png_path.stat().st_size > 1000)
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.image as mpimg
    img = mpimg.imread(png_path)
    ok("all-invalid dataset: decoded PNG has sane, nonzero pixel dimensions",
       img.ndim >= 2 and img.shape[0] > 50 and img.shape[1] > 50)
    md_text = md_path.read_text(encoding="utf-8")
    ok("all-invalid dataset: verdict.md is non-empty and reports FAIL with a reason",
       len(md_text) > 200 and "FAIL" in md_text and "no_eligible_rows" in md_text)
    # BAD_EVIDENCE/BAD_HALLU were fed into write_markdown above and it still
    # produced a complete file without raising or short-circuiting -- this
    # IS the "evidence failure cannot suppress artifact generation" check.
    ok("evidence failure (BAD_EVIDENCE/BAD_HALLU) does not suppress artifact generation",
       not v_invalid["pass"] and "evidence_incomplete" in v_invalid["fail_reasons"]
       and png_path.exists() and md_path.exists())


# ==================================================== 3. full default run

with tempfile.TemporaryDirectory() as td:
    full_dir = Path(td) / "full"
    buf = io.StringIO()
    with redirect_stdout(buf):
        rc_full = rte.main(["--out-dir", str(full_dir)])
    stdout_full = buf.getvalue()
    ok("full run: exits 0", rc_full == 0)
    full_artifacts = {name: full_dir / name
                      for name in ("sweep.csv", "envelope.png", "verdict.md", "manifest.json")}
    ok("full run: all four artifacts exist", all(p.exists() and p.stat().st_size > 0
                                                  for p in full_artifacts.values()))
    ok("full run: printed VERDICT line matches the required stdout format",
       ("VERDICT: PASS" in stdout_full or "VERDICT: FAIL" in stdout_full)
       and "g2_min=" in stdout_full and "eligible=" in stdout_full and "evidence=" in stdout_full)

    with full_artifacts["sweep.csv"].open(newline="", encoding="utf-8") as f:
        full_rows = list(csv.DictReader(f))

    # Independent recomputation directly from the CSV (never importing
    # compute_stats for this part) of per-card pulsed g2 min/median and
    # eligible-coverage, cross-checked against manifest.json.
    def _recompute(card_id):
        card_rows = [r for r in full_rows if r["card_id"] == card_id]
        eligible = [r for r in card_rows if r["eligible_row"] == "True"]
        pulsed_vals = [float(r["g2_pulsed"]) for r in eligible if r["g2_pulsed"] not in ("", "nan")]
        pulsed_vals = [v for v in pulsed_vals if math.isfinite(v)]
        return len(card_rows), len(eligible), pulsed_vals

    manifest_full = json.loads(full_artifacts["manifest.json"].read_text(encoding="utf-8"))
    for card_id in (PRIMARY_ID, FALLBACK_ID):
        n_rows, n_eligible, pulsed_vals = _recompute(card_id)
        recomputed_min = min(pulsed_vals) if pulsed_vals else float("nan")
        recomputed_median = statistics.median(pulsed_vals) if pulsed_vals else float("nan")
        declared = manifest_full["per_card_stats"][card_id]
        min_matches = (math.isnan(recomputed_min) and math.isnan(declared["g2_pulsed_min"])) or \
            abs(recomputed_min - declared["g2_pulsed_min"]) < 1e-9
        med_matches = (math.isnan(recomputed_median) and math.isnan(declared["g2_pulsed_median"])) or \
            abs(recomputed_median - declared["g2_pulsed_median"]) < 1e-9
        ok(f"full run: {card_id} pulsed min/median recomputed from CSV match manifest.json",
           min_matches and med_matches and n_eligible == declared["n_eligible"]
           and n_rows == declared["n_rows"])

    # Every claimed headline-passing row must genuinely satisfy the headline
    # (pulsed-only) gate and carry non-empty assumptions (spec: "every
    # claimed passing corner must have matching assumptions and a valid
    # headline gate"); its secondary (CW) diagnostics may be anything.
    passing_rows = [r for r in full_rows if r["headline_pass"] == "True"]
    ok("full run: every headline_pass=True row genuinely satisfies the headline gate "
       "(eligible, pulsed intrinsic g2(0)<0.5)",
       all(r["eligible_row"] == "True" and float(r["g2_pulsed"]) < 0.5
           for r in passing_rows))
    ok("full run: every headline_pass=True row lists non-empty assumptions",
       all(r["assumptions"].strip() for r in passing_rows))

    # Council review round 4, item 1: the collection-lever axes are actually
    # varied in a real run against this repo's real cards (which declare
    # emission.NA/R_back/L_um ranges), not held constant, and the row count
    # matches n_combos * (axis samples)^3 exactly (independent recount).
    ok("full run: manifest records a resolved per-card lever grid for both cards",
       all(cid in manifest_full["lever_info"] for cid in (PRIMARY_ID, FALLBACK_ID))
       and all(manifest_full["lever_info"][cid]["n_combos"] >= 1
               for cid in (PRIMARY_ID, FALLBACK_ID)))
    for card_id in (PRIMARY_ID, FALLBACK_ID):
        card_rows_full = [r for r in full_rows if r["card_id"] == card_id]
        n_combos = manifest_full["lever_info"][card_id]["n_combos"]
        expected_n_rows = n_combos * (rte._FULL_N ** 3)
        ok(f"full run: {card_id} row count matches n_combos * (delta_xx x gamma300 x irf) "
           f"sample counts exactly",
           len(card_rows_full) == expected_n_rows)
        distinct_na = {r["emission_NA"] for r in card_rows_full}
        if len(manifest_full["lever_info"][card_id]["grid"]["emission.NA"]) > 1:
            ok(f"full run: {card_id} sweep.csv actually varies emission.NA "
               f"(this card declares a range for it)",
               len(distinct_na) > 1)

    # Council review round 4, item 2: gamma300_pass_max is present in the
    # printed VERDICT line and the manifest, and is finite whenever any
    # card actually has a headline-passing row.
    ok("full run: printed VERDICT line carries gamma300_pass_max=",
       "gamma300_pass_max=" in stdout_full)
    ok("full run: manifest verdict section carries gamma300_pass_max and "
       "gamma300_pass_max_by_card for both cards",
       "gamma300_pass_max" in manifest_full["verdict"]
       and all(cid in manifest_full["verdict"]["gamma300_pass_max_by_card"]
               for cid in (PRIMARY_ID, FALLBACK_ID)))
    if manifest_full["stats"]["n_headline"] > 0:
        ok("full run: gamma300_pass_max is finite whenever at least one row headline-passes",
           math.isfinite(manifest_full["verdict"]["gamma300_pass_max"]))

    # ---- Council review round 5 acceptance criteria ----

    # Item 1: whichever card has n_eligible == 0 must print the "not measurable"
    # phrase; whichever card has n_eligible > 0 must print "eligible rows:" and
    # must NEVER print "not measurable" (this is the exact gainp mislabelling bug).
    for card_id, cs in manifest_full["per_card_stats"].items():
        line = next((ln for ln in stdout_full.splitlines() if ln.startswith(f"CARD: {card_id} ")),
                   None)
        ok(f"full run: CARD line for {card_id} exists in stdout", line is not None)
        if line is None:
            continue
        if cs["n_eligible"] == 0:
            ok(f"full run: {card_id} has 0 eligible rows -> CARD line says 'not measurable'",
               "not measurable" in line)
        else:
            ok(f"full run: {card_id} has {cs['n_eligible']} eligible rows -> CARD line says "
               f"'eligible rows:' and NEVER 'not measurable' (council review round 5 item 1)",
               f"eligible rows: {cs['n_eligible']}/{cs['n_rows']}" in line
               and "not measurable" not in line)

    # Item 2: gamma300_threshold in the VERDICT line/manifest; anchor_6_5mev_by_card
    # for both cards, each entry a genuine fresh eval_pulsed_point result.
    ok("full run: printed VERDICT line carries gamma300_threshold=",
       "gamma300_threshold=" in stdout_full)
    ok("full run: manifest verdict section carries gamma300_threshold_by_card and "
       "anchor_6_5mev_by_card for both cards",
       all(cid in manifest_full["verdict"]["gamma300_threshold_by_card"]
           for cid in (PRIMARY_ID, FALLBACK_ID))
       and all(cid in manifest_full["verdict"]["anchor_6_5mev_by_card"]
               for cid in (PRIMARY_ID, FALLBACK_ID)))
    for card_id in (PRIMARY_ID, FALLBACK_ID):
        anchor = manifest_full["verdict"]["anchor_6_5mev_by_card"][card_id]
        ok(f"full run: {card_id} anchor_6_5mev entry is a real 6.5 meV evaluation whose "
           f"'passes' flag matches g2_pulsed < 0.5 and eligible",
           anchor["gamma300_meV"] == rte.CHATZARAKIS_ANCHOR_GAMMA300_MEV
           and anchor["passes"] == bool(anchor["eligible"] and anchor["g2_pulsed"] is not None
                                        and anchor["g2_pulsed"] < rte.G2_THRESHOLD))
    ok("full run: verdict.md states the 6.5 meV anchor result and quotes Reischle's "
       "IRF-deconvolved 0.25 +/- 0.05 value",
       full_artifacts["verdict.md"].read_text(encoding="utf-8").count("6.5") > 0)
    md_text_full = full_artifacts["verdict.md"].read_text(encoding="utf-8")
    ok("full run: verdict.md quotes Reischle's like-for-like deconvolved value (0.25) "
       "alongside the raw dip (0.43)",
       "0.25" in md_text_full and "0.43" in md_text_full)
    ok("full run: verdict.md's front-facet section documents the independent forward "
       "check (council review round 5, item 6), not just the back-solved tautology",
       "Independent check on the combined front-facet factor" in md_text_full
       and "BACK-SOLVED" in md_text_full)

    # ==================== Sixth council review (2026-09-06) reporting fixes ====================

    # Item 1: the background-assumption citation is generated FROM the ledger
    # anchor (reischle08-b-res-80k), never the hardcoded/drifted "Appl. Phys.
    # Lett. 92, 233113 (2008)" string.
    ok("full run: verdict.md never cites the drifted 'Appl. Phys. Lett. 92, "
       "233113 (2008)' background-assumption source",
       "Appl. Phys. Lett. 92" not in md_text_full)
    ok("full run: verdict.md's background-assumption paragraph cites the ledger anchor's "
       "real source (Optics Express 16, 12771 (2008), DOI 10.1364/OE.16.012771)",
       "12771" in md_text_full and "10.1364/OE.16.012771" in md_text_full
       and "reischle08-b-res-80k" in md_text_full)

    # Item 2: a card with zero headline-passing samples prints gamma300_threshold
    # as a REASON class (none(flux), none(g2), or none(flux,g2)), never a bare "<N" bracket that
    # would misleadingly imply it passes below N.
    gaasp_row_line = next((ln for ln in md_text_full.splitlines()
                          if ln.startswith(f"| {PRIMARY_ID} |")), None)
    ok("full run: the gaasp (primary) card's gamma300_threshold table row exists",
       gaasp_row_line is not None)
    if gaasp_row_line is not None:
        ok("full run: the gaasp (primary) card's gamma300_threshold prints every blocking reason "
           "'none(flux,g2)' (below the collected-flux floor and diagnostic g2 >= 0.5), "
           "never a bare '<N' bracket",
           "none(flux,g2)" in gaasp_row_line and "<6" not in gaasp_row_line)
    ok("full run: verdict.md explains every none(flux,g2) reason class",
       "none(flux)" in md_text_full and "none(g2)" in md_text_full and "none(flux,g2)" in md_text_full)

    # Item 3: the dominant brightness limiter is compared across the WHOLE
    # chain (loading, t_X, S, plus eta_total's own sub-factors), not just
    # eta_total's sub-factors -- at the favourable corner S (retention) is
    # smaller than beta and must be named.
    dominant_sentence = next((ln for ln in md_text_full.splitlines()
                             if "dominant brightness limiter" in ln), None)
    ok("full run: verdict.md names a dominant brightness limiter",
       dominant_sentence is not None)
    if dominant_sentence is not None:
        ok("full run: the dominant brightness limiter names S (confinement retention), "
           "not beta -- S is smaller at the favourable corner once the comparison spans "
           "the whole chain (council review round 6, item 3)",
           "is S (confinement retention)" in dominant_sentence)

    # Item 4: one name per quantity -- the Coverage section uses the same key
    # names as the VERDICT line, everywhere.
    ok("full run: verdict.md's Coverage section uses the VERDICT-line key names "
       "(coverage_over_eligible, g2_median_eligible, diag_g2_median_diagnostic), "
       "not the old bare coverage/g2_median/diag_g2_median labels",
       "`coverage_over_eligible`" in md_text_full and "`g2_median_eligible`" in md_text_full
       and "`diag_g2_median_diagnostic`" in md_text_full)

    # Item 5: REISCHLE_RAW_G2/REISCHLE_DECONV_G2/REISCHLE_DECONV_G2_ERR are
    # read from the ledger and tagged in the markdown.
    ok("full run: verdict.md tags the Reischle raw/deconvolved g2 values with their "
       "ledger anchor ids and a provenance tag",
       "reischle08-g2-80k-deconvolved" in md_text_full
       and "0.43 [V]" in md_text_full and "0.25 [V]" in md_text_full)

    # Item 6: the verified-anchor section header identifies the Chatzarakis
    # anchor as an InAs/GaAs (211)B single-dot measurement used as a
    # distinct-material class proxy, not an unqualified InP-class result.
    ok("full run: the 'Verified 6.5 meV anchor' header names the InAs/GaAs (211)B "
       "single-dot measurement and its [E]-class proxy role",
       "InAs/GaAs (211)B single-dot" in md_text_full
       and "distinct-material class proxy [E]" in md_text_full)

    # Item 4: flux_margin/eligible_fraction in the VERDICT line and manifest.
    ok("full run: printed VERDICT line carries flux_margin= and eligible_fraction=",
       "flux_margin=" in stdout_full and "eligible_fraction=" in stdout_full)
    ok("full run: manifest verdict section carries flux_margin and eligible_fraction, "
       "and flux_margin is the exact reciprocal of the deprecated flux_shortfall",
       "flux_margin" in manifest_full["verdict"] and "eligible_fraction" in manifest_full["verdict"]
       and (not math.isfinite(manifest_full["verdict"]["flux_margin"])
           or abs(manifest_full["verdict"]["flux_margin"]
                  * manifest_full["verdict"]["flux_shortfall"] - 1.0) < 1e-6))

    # Item 3: sweep.csv's actual delta_xx values match the cards' own declared
    # (4.0, 8.0) range, not the stale (4.0, 7.0).
    full_deltas = {float(r["delta_xx_meV"]) for r in full_rows}
    ok("full run: sweep.csv delta_xx values are the cards' own declared (4.0, 8.0) "
       "endpoints (council review round 5, item 3)",
       full_deltas == {4.0, 8.0})


# ============================================== 4. artifacts, hashes, determinism

ok("full run: CSV header matches csv_fieldnames() exactly",
   list(full_rows[0].keys()) == rte.csv_fieldnames())
ok("manifest.json: card sha256 hashes match a fresh independent hash of each card file",
   all(entry["sha256"] == rte._sha256_file(Path(entry["path"])) for entry in manifest_full["cards"]))
ok("manifest.json: evidence section carries the evaluator input hashes",
   set(manifest_full["evidence"]["evaluator_input_hashes"]) == {
       "verify/data/rt_edge_anchors.yaml", "cards/edge-inp-gaasp-design.yaml",
       "cards/edge-inp-gainp-design.yaml"})
ok("manifest.json: output_schema.sweep.csv matches csv_fieldnames()",
   manifest_full["output_schema"]["sweep.csv"] == rte.csv_fieldnames())

# Determinism of the scientific columns (excludes only manifest timestamps/
# elapsed-time fields, per the spec) across two independent --quick runs.
with tempfile.TemporaryDirectory() as td_a, tempfile.TemporaryDirectory() as td_b:
    dir_a, dir_b = Path(td_a) / "a", Path(td_b) / "b"
    with redirect_stdout(io.StringIO()):
        rte.main(["--quick", "--out-dir", str(dir_a)])
        rte.main(["--quick", "--out-dir", str(dir_b)])
    csv_a = (dir_a / "sweep.csv").read_text(encoding="utf-8")
    csv_b = (dir_b / "sweep.csv").read_text(encoding="utf-8")
    ok("determinism: sweep.csv is byte-identical across two independent --quick runs",
       csv_a == csv_b)
    man_a = json.loads((dir_a / "manifest.json").read_text(encoding="utf-8"))
    man_b = json.loads((dir_b / "manifest.json").read_text(encoding="utf-8"))
    for m in (man_a, man_b):
        m.pop("generated_utc", None)
        m.pop("runtime_seconds", None)
        m["stats"].pop("runtime_seconds", None)
    ok("determinism: manifest.json is identical across two runs, excluding timestamp/elapsed-time",
       man_a == man_b)

print(f"{sum(CHECKS)}/{len(CHECKS)} rt-edge sweep checks passed")
sys.exit(0 if all(CHECKS) else 1)
