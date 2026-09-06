"""Checks for scripts/run_rt_edge.py: the RT edge-emitter acceptance-sweep
policy (compute_verdict/compute_stats), the artifacts it writes (sweep.csv,
envelope.png, verdict.md, manifest.json), and the real-evaluator smoke paths
(--quick and the full default grid).

This file does NOT re-derive device physics (that is verify/verify_device_rt.py
and verify/verify_rt_edge_cards.py's job); it checks that run_rt_edge.py's own
bookkeeping -- grid construction, eligibility classification, min/median/
coverage aggregation, and the PASS/FAIL policy in compute_verdict -- is
correct, using synthetic fixtures whose expected outcome is worked out by hand
in this file (never by calling compute_verdict and trusting its own answer),
plus real (but --quick, then once full) runs of the actual script.

Section 1 (pure policy, no evaluate() calls): all-ineligible, no-evidence,
duplicate-evidence, headline-only-pass (favorable pulsed gate despite both CW
secondary diagnostics failing), secondary-only-pass (favorable CW diagnostics
cannot substitute for a failing headline gate), coverage bookkeeping,
fallback-only-pass, median-fail and only-[A]/[E]-corner scenarios against
synthetic rows -- exercising docs/rt_edge_contract.md's rule that the
headline metric (pulsed intrinsic g2(0)) alone gates PASS, with g2_cw0/
g2_cw0_raw reported as non-gating secondary diagnostics.
Section 2: a --quick real-evaluator smoke run, plus a synthetic all-invalid
dataset run through write_png/write_markdown directly (an empty/all-invalid
case must still produce legible artifacts, and evidence failure must not
suppress them).
Section 3: one full (non-quick) real-evaluator run; statistics are recomputed
directly from the written CSV, independently of compute_stats.
Section 4: CSV/PNG/JSON structural checks, hash cross-checks against
manifest.json, and determinism of the scientific columns (excluding elapsed
time/timestamps) across two independent --quick runs.

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
            assumptions="dot.gamma300; emission.lambda_nm") -> dict:
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
            "assumptions": assumptions}


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
        "b_e_window_cw": 0.1, "t_x_pulsed": 0.5, "eps_pulsed": 0.1, "t_x_cw": 0.5,
        "eps_cw": 0.1, "edge_beta": 0.01, "edge_eta_total": 0.002,
        "collected_flux_pulsed_s": 10.0, "collected_flux_cw_s": 10.0,
        "g2_pulsed": float("nan"), "g2_cw0": float("nan"), "g2_cw0_raw": float("nan"),
        "eligible_pulsed": False, "eligible_cw": False, "eligible_row": False,
        "headline_pass": False, "secondary_pass": False,
        "invalid_reasons_pulsed": "synthetic: forced invalid",
        "invalid_reasons_cw": "synthetic: forced invalid",
        "assumptions": "dot.gamma300; emission.lambda_nm",
        "pulse_width_ns": rte.PULSE_WIDTH_NS, "rep_rate_hz": rte.REP_RATE_HZ,
        "duty_pulsed": rte.PULSE_WIDTH_NS * 1e-9 * rte.REP_RATE_HZ,
    }
    base.update(overrides)
    return base


# ============================================================ 1. pure policy

# -- all-ineligible: no row eligible -> FAIL for eligibility, nan stats
rows = [make_row(PRIMARY_ID, "primary", 0.3, 0.3, 0.3, False),
        make_row(PRIMARY_ID, "primary", 0.9, 0.9, 0.9, False)]
stats = rte.compute_stats(rows)
v = rte.compute_verdict(rows, stats, True, GOOD_EVIDENCE, GOOD_HALLU)
ok("all-ineligible: zero eligible rows despite favorable-looking g2 values",
   stats["n_eligible"] == 0 and math.isnan(v["g2_min"]) and math.isnan(v["g2_median"]))
ok("all-ineligible: verdict FAILs on eligibility, not metrics",
   not v["pass"] and "no_eligible_rows" in v["fail_reasons"])

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
with tempfile.TemporaryDirectory() as td:
    md_path = Path(td) / "verdict.md"
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

# -- only-[A]/[E]-corner: a favorable corner's assumptions carry into "conditional"
rows = [make_row(PRIMARY_ID, "primary", 0.3, 0.3, 0.3, True, assumptions="dot.gamma300; emission.lambda_nm")]
stats = rte.compute_stats(rows)
v = rte.compute_verdict(rows, stats, True, GOOD_EVIDENCE, GOOD_HALLU)
ok("only-[A]/[E]-corner: PASS is flagged conditional on its stated assumptions",
   v["pass"] and v["conditional"] and set(v["assumptions_used"]) == {"dot.gamma300", "emission.lambda_nm"})
rows_noassum = [make_row(PRIMARY_ID, "primary", 0.3, 0.3, 0.3, True, assumptions="")]
v_noassum = rte.compute_verdict(rows_noassum, rte.compute_stats(rows_noassum), True,
                                GOOD_EVIDENCE, GOOD_HALLU)
ok("conditional is False only when the passing corner truly carries no assumptions",
   v_noassum["pass"] and not v_noassum["conditional"])

# -- an incomplete (--quick) grid can never PASS, regardless of how good the metrics are
v_quick = rte.compute_verdict(rows, stats, False, GOOD_EVIDENCE, GOOD_HALLU)
ok("incomplete grid cannot grant PASS even with a favorable, well-evidenced corner",
   not v_quick["pass"] and "grid_incomplete" in v_quick["fail_reasons"])


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
    ok("--quick: declared endpoints are sampled on every axis",
       deltas == {4.0, 7.0} and gammas == {6.0, 20.0} and irfs == {50.0, 200.0})
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
