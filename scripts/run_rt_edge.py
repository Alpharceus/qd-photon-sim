"""RT edge-emitter acceptance sweep (integration-c): evaluates both room-
temperature electrical edge-emitter design cards (cards/edge-inp-gaasp-
design.yaml, cards/edge-inp-gainp-design.yaml) over the contract-declared
ranges (docs/rt_edge_contract.md; endpoints match verify/verify_rt_edge_cards.py
REQUIRED_RANGES / each card's own provenance.ranges) and renders a
reproducible 300 K heatsink acceptance envelope: sweep.csv, envelope.png,
verdict.md, manifest.json.

This file is a THIN caller of the validated evaluator (fsim_core.device.
evaluate via resolve_device_card below) -- it contains no confinement,
transport, thermal or optical physics of its own. Every number in a CSV row
traces to a single evaluate() call's scalars; the only arithmetic performed
here is bookkeeping (grid construction, min/median/coverage aggregation,
plotting, an [A] repetition-rate multiplication for reported flux).

Two duty states are evaluated at every (delta_xx, gamma300) grid point,
matching docs/rt_edge_contract.md's b_e/T_j conventions and the spec's
"pair each pulsed operating point with a duty=1 DC evaluation":

  pulsed  -- drive.duty = tau_pulse_ns * rep_rate_Hz, drive.diode.tau_pulse_ns
             = tau_pulse_ns (transport.dot_loading's mu = r_dot * tau_pulse_ns,
             the ONLY place a pulse width enters the model), drive.cw=False.
             Headline metric: g2_op (pulsed intrinsic g2(0)).
  cw      -- drive.duty = 1.0 (continuous), drive.cw=True, drive.cw_irf_fwhm_ps
             = the swept irf_ps value. Uses injection.loading.r_dot directly
             (cw_g2's rate-equation model, wired inside evaluate()) -- never
             the pulsed mu. Secondary metrics: g2_cw0 (intrinsic+background),
             g2_cw0_raw (IRF-convolved).

The pulsed sub-result does not depend on irf_ps at all, so it is cached per
(card, delta_xx, gamma300) and reused across the irf_ps axis -- this is the
"avoid a full expensive CW recalculation for every identical IRF-independent
input" reuse the spec asks for: the pulsed evaluate() call has drive.cw=False
(skips the CW rate-equation/dip-fit block entirely) and is never repeated for
an unchanged (card, delta_xx, gamma300) triple.

pulse_width_ns and rep_rate_hz are themselves [A] assumptions (no published
single-dot InP/GaAsP or InP/GaInP pulsed drive scheme exists at this current
scale; see PULSE_ASSUMPTION below) -- they are recorded in every row's
assumptions column and in manifest.json, and are used ONLY for (a) the pulsed
duty/tau_pulse_ns inputs above and (b) the reported collected_flux_s column
(brightness_per_pulse * rep_rate_hz); they never feed g2 or Tj directly.

CLI: python scripts/run_rt_edge.py [--quick] [--out-dir PATH]
--quick uses a 2-point (endpoints-only) grid, explicitly marked incomplete;
an incomplete grid can never grant VERDICT: PASS (see compute_verdict).

Exit code: 0 whenever the four artifacts were generated (including a
scientific VERDICT: FAIL -- the machine-readable verdict, not the process
exit code, is authoritative for acceptance); nonzero only for a runtime or
artifact-writing error.

Standalone, side-effect-free on import (all work happens under
`if __name__ == "__main__"`); verify/verify_rt_edge_sweep.py imports the
functions below directly (grid, resolve_device_card, eval_pulsed_point,
eval_cw_point, compute_stats, compute_verdict, write_csv, write_png,
write_markdown, write_manifest) rather than parsing this script's stdout.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import platform
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import scipy
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.device import DeviceDesign, evaluate  # noqa: E402
import verify.verify_rt_edge_papers as rt_papers  # noqa: E402

CARDS = [
    {"id": "edge-inp-gaasp-design", "class": "primary",
     "path": ROOT / "cards" / "edge-inp-gaasp-design.yaml"},
    {"id": "edge-inp-gainp-design", "class": "fallback",
     "path": ROOT / "cards" / "edge-inp-gainp-design.yaml"},
]

# Contract ranges (docs/rt_edge_contract.md; exact endpoints checked
# mechanically against every card's own provenance.ranges by
# verify/verify_rt_edge_cards.py REQUIRED_RANGES -- reproduced here, not
# re-derived, so a range typo in either place would disagree and be caught
# by that file, not silently drift).
RANGE_BOUNDS = {
    "dot.delta_xx": (4.0, 7.0, "meV"),
    "dot.gamma300": (6.0, 20.0, "meV"),
    "irf_ps": (50.0, 200.0, "ps"),
}

# [A] pulsed-drive assumption: no published single-dot InP/GaAsP or
# InP/GaInP pulsed excitation scheme exists at this design's pA-nA current
# scale (cards/edge-inp-*-design.yaml provenance.sources["drive.I_uA"]).
# 100 ps is the fast-gain-switched/mode-locked-diode class window used
# elsewhere in this repo's own RT analysis (scripts/run_rt_campaign.py T-5
# pulse-edge section, "100-200 ps windows"); 80 MHz is the standard mode-
# locked-source repetition-rate class (Ti:Sapphire / gain-switched diode).
# Used ONLY for (a) drive.diode.tau_pulse_ns / drive.duty at the pulsed sub-
# point and (b) the reported collected_flux_s column -- never for g2 itself.
PULSE_WIDTH_NS = 0.1
REP_RATE_HZ = 80.0e6
PULSE_ASSUMPTION = {
    "pulse_width_ns": PULSE_WIDTH_NS, "rep_rate_hz": REP_RATE_HZ, "tag": "A",
    "note": ("no published single-dot pulsed-drive scheme exists at this "
             "design's pA-nA current scale; 100 ps / 80 MHz is the fast-"
             "gain-switched/mode-locked class window (cf. scripts/"
             "run_rt_campaign.py T-5), used only for drive.duty/"
             "diode.tau_pulse_ns and the reported collected_flux_s column"),
}
SWEEP_ASSUMPTION_KEYS = ["drive.diode.tau_pulse_ns", "drive.duty (rep-rate-derived)"]

G2_THRESHOLD = 0.5
_ENDPOINT_N = 2
_FULL_N = 3


# --------------------------------------------------------------------- grid

def build_grid(quick: bool) -> dict:
    """{"dot.delta_xx": [...], "dot.gamma300": [...], "irf_ps": [...]}, always
    including both declared endpoints; quick=True is endpoints-only (spec:
    "explicitly marked incomplete and cannot grant scientific PASS"); the
    full grid adds one interior sample per axis (device.py evaluate_envelope's
    own caveat: a bare 2-level factorial only bounds a monotone response)."""
    # n=3 -> endpoints plus one interior (midpoint) sample per axis: a
    # genuine (non-degenerate) grid, catching a non-monotone response the
    # bare 2-level factorial (quick=True) cannot.
    n = _ENDPOINT_N if quick else _FULL_N
    return {path: [float(v) for v in np.linspace(lo, hi, n)]
            for path, (lo, hi, _unit) in RANGE_BOUNDS.items()}


# ------------------------------------------------------------ card plumbing

def resolve_device_card(card_path: Path, overrides: dict) -> DeviceDesign:
    """Load a design card and apply explicit dotted-path overrides (same
    two-segment block.field convention as fsim_core.device._set_path) --
    pure plumbing, no physics: every number that reaches evaluate() below is
    either the card's own value or one of these caller-supplied overrides."""
    design = DeviceDesign.load(card_path)
    for path, value in overrides.items():
        block_name, field_name = path.split(".", 1)
        setattr(getattr(design, block_name), field_name, value)
    return design


def _sha256_file(path: Path) -> str:
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


# --------------------------------------------------------------- evaluation

def eval_pulsed_point(card_path: Path, delta_xx: float, gamma300: float,
                      cache: dict | None = None) -> dict:
    """Pulsed sub-result at (delta_xx, gamma300); irf_ps-independent, so
    cache is keyed without it and reused across the irf_ps sweep axis."""
    key = (str(card_path), "pulsed", delta_xx, gamma300)
    if cache is not None and key in cache:
        return cache[key]
    base = DeviceDesign.load(card_path)
    duty = PULSE_WIDTH_NS * 1e-9 * REP_RATE_HZ
    design = resolve_device_card(card_path, {
        "dot.delta_xx": delta_xx, "dot.gamma300": gamma300,
        "drive.duty": duty, "drive.cw": False,
        "drive.diode": {**base.drive.diode, "tau_pulse_ns": PULSE_WIDTH_NS},
    })
    sc = evaluate(design, T_grid=[design.thermal.T_hs])["scalars"]
    flux = _collected_flux_s(sc)
    eligible, reasons = _classify(sc, flux, "pulsed")
    result = {"scalars": sc, "flux_s": flux, "eligible": eligible,
              "reasons": reasons, "duty": duty}
    if cache is not None:
        cache[key] = result
    return result


def eval_cw_point(card_path: Path, delta_xx: float, gamma300: float,
                  irf_ps: float) -> dict:
    """CW (duty=1 DC) sub-result; genuinely irf_ps-dependent (g2_cw0_raw),
    so every irf_ps grid value gets its own evaluate() call."""
    design = resolve_device_card(card_path, {
        "dot.delta_xx": delta_xx, "dot.gamma300": gamma300,
        "drive.duty": 1.0, "drive.cw": True, "drive.cw_irf_fwhm_ps": irf_ps,
    })
    sc = evaluate(design, T_grid=[design.thermal.T_hs])["scalars"]
    flux = _collected_flux_s(sc)
    eligible, reasons = _classify(sc, flux, "cw")
    return {"scalars": sc, "flux_s": flux, "eligible": eligible, "reasons": reasons}


def _collected_flux_s(sc: dict) -> float:
    """[A] reporting-only quantity: brightness_per_pulse (device.py's own
    fully-composed per-event collected photon yield -- loading, spectral
    transmission, retention and edge collection each exactly once, per
    device.py's Lemma 1 comment) times the [A] repetition rate. Never fed
    back into g2/eps/rho/Tj."""
    b = sc.get("brightness_per_pulse", float("nan"))
    return float(b) * REP_RATE_HZ if np.isfinite(b) else float("nan")


def _classify(sc: dict, flux_s: float, role: str) -> tuple[bool, list]:
    """Eligibility per docs/rt_edge_contract.md / the spec's interface
    constraints: supported confinement/transport/thermal/optical results,
    finite diode V/I with positive operating injection, eta_inj>0 and
    collected_flux_s>0. device.py's own invalid_reasons already covers
    thermal nonconvergence/runaway, zero/negative current, and absorbing/
    unguided/unsupported emission stacks (edge_err); the checks below are
    additional, never a floor or a way to make an otherwise-invalid row
    eligible."""
    reasons = list(sc.get("invalid_reasons", []))
    if not np.isfinite(sc.get("T_j_op", float("nan"))):
        reasons.append(f"{role}: non-finite T_j_op")
    if not np.isfinite(sc.get("V_j_op", float("nan"))):
        reasons.append(f"{role}: non-finite diode V_j_op")
    if not (sc.get("eta_inj", float("nan")) > 0):
        reasons.append(f"{role}: eta_inj not positive")
    if not (np.isfinite(flux_s) and flux_s > 0):
        reasons.append(f"{role}: collected_flux_s not positive")
    if role == "pulsed" and not np.isfinite(sc.get("g2_op", float("nan"))):
        reasons.append("pulsed: g2_op is NaN (zero signal or invalid operating point)")
    if role == "cw":
        if not np.isfinite(sc.get("g2_cw0", float("nan"))):
            reasons.append("cw: g2_cw0 is NaN")
        if not np.isfinite(sc.get("g2_cw0_raw", float("nan"))):
            reasons.append("cw: g2_cw0_raw is NaN")
    # de-duplicate while preserving order (invalid_reasons may already
    # overlap with the CW-requested-but-not-evaluable message above)
    seen, deduped = set(), []
    for r in reasons:
        if r not in seen:
            seen.add(r)
            deduped.append(r)
    return len(deduped) == 0, deduped


# ------------------------------------------------------------------ sweep

def sweep_card(card: dict, grid: dict, pulsed_cache: dict) -> list:
    """One card's full cartesian (delta_xx, gamma300, irf_ps) grid -> rows.
    Every scheduled row is emitted, eligible or not (spec: never drop invalid
    rows)."""
    card_path = card["path"]
    design0 = DeviceDesign.load(card_path)
    provenance = design0.provenance or {}
    ranges = provenance.get("ranges", {})
    card_assumptions = list(provenance.get("assumptions", []))
    i_ua = design0.drive.I_uA
    rows = []
    for delta_xx in grid["dot.delta_xx"]:
        for gamma300 in grid["dot.gamma300"]:
            pulsed = eval_pulsed_point(card_path, delta_xx, gamma300, pulsed_cache)
            for irf_ps in grid["irf_ps"]:
                cw = eval_cw_point(card_path, delta_xx, gamma300, irf_ps)
                rows.append(_build_row(card, ranges, card_assumptions, i_ua,
                                       delta_xx, gamma300, irf_ps, pulsed, cw))
    return rows


def _build_row(card, ranges, card_assumptions, i_ua, delta_xx, gamma300, irf_ps,
              pulsed, cw) -> dict:
    scp, sccw = pulsed["scalars"], cw["scalars"]
    eligible_row = bool(pulsed["eligible"] and cw["eligible"])
    g2_p, g2_cw0, g2_cw0_raw = scp.get("g2_op"), sccw.get("g2_cw0"), sccw.get("g2_cw0_raw")
    # Contract rule (docs/rt_edge_contract.md Acceptance gates): the headline
    # metric is pulsed intrinsic g2(0) alone -- g2_cw0/g2_cw0_raw are
    # secondary diagnostics, reported per-row and as coverage fractions, but
    # they never gate PASS. `secondary_pass` is retained as a diagnostic-only
    # column (both CW gates favorable), never consumed by compute_verdict.
    headline_pass = bool(eligible_row and np.isfinite(g2_p) and g2_p < G2_THRESHOLD)
    secondary_pass = bool(
        eligible_row and np.isfinite(g2_cw0) and np.isfinite(g2_cw0_raw)
        and g2_cw0 < G2_THRESHOLD and g2_cw0_raw < G2_THRESHOLD)
    assumptions = sorted(set(card_assumptions) | set(SWEEP_ASSUMPTION_KEYS))
    config_hash = hashlib.sha256(json.dumps(
        {"card": card["id"], "delta_xx": delta_xx, "gamma300": gamma300,
         "irf_ps": irf_ps, "pulse_width_ns": PULSE_WIDTH_NS, "rep_rate_hz": REP_RATE_HZ},
        sort_keys=True).encode()).hexdigest()[:16]
    return {
        "card_id": card["id"], "card_class": card["class"], "config_id": config_hash,
        "delta_xx_meV": delta_xx, "delta_xx_tag": ranges.get("dot.delta_xx", {}).get("tag", ""),
        "gamma300_meV": gamma300, "gamma300_tag": ranges.get("dot.gamma300", {}).get("tag", ""),
        "irf_ps": irf_ps, "irf_tag": ranges.get("irf_ps", {}).get("tag", ""),
        "T_hs_K": scp.get("T_j_op", float("nan")) - scp.get("dT_J", 0.0),
        "Tj_pulsed_K": scp.get("T_j_op"), "Tj_cw_K": sccw.get("T_j_op"),
        "I_uA": i_ua,
        "V_j_pulsed_V": scp.get("V_j_op"), "V_j_cw_V": sccw.get("V_j_op"),
        "mu_pulsed": scp.get("mu_resolved"),
        "eta_inj_pulsed": scp.get("eta_inj"), "eta_inj_cw": sccw.get("eta_inj"),
        "S_retention_pulsed": scp.get("S_resolved"), "S_retention_cw": sccw.get("S_resolved"),
        "b_e_window_pulsed": scp.get("b_e_resolved"), "b_e_window_cw": sccw.get("b_e_resolved"),
        "t_x_pulsed": scp.get("t_x_op"), "eps_pulsed": scp.get("eps_op"),
        "t_x_cw": sccw.get("t_x_op"), "eps_cw": sccw.get("eps_op"),
        "edge_beta": scp.get("edge_beta"), "edge_eta_total": scp.get("edge_eta_total"),
        "collected_flux_pulsed_s": pulsed["flux_s"], "collected_flux_cw_s": cw["flux_s"],
        "g2_pulsed": g2_p, "g2_cw0": g2_cw0, "g2_cw0_raw": g2_cw0_raw,
        "eligible_pulsed": pulsed["eligible"], "eligible_cw": cw["eligible"],
        "eligible_row": eligible_row, "headline_pass": headline_pass,
        "secondary_pass": secondary_pass,
        "invalid_reasons_pulsed": "; ".join(pulsed["reasons"]),
        "invalid_reasons_cw": "; ".join(cw["reasons"]),
        "assumptions": "; ".join(assumptions),
        "pulse_width_ns": PULSE_WIDTH_NS, "rep_rate_hz": REP_RATE_HZ,
        "duty_pulsed": pulsed["duty"],
    }


def csv_fieldnames() -> list:
    return ["card_id", "card_class", "config_id", "delta_xx_meV", "delta_xx_tag",
            "gamma300_meV", "gamma300_tag", "irf_ps", "irf_tag", "T_hs_K",
            "Tj_pulsed_K", "Tj_cw_K", "I_uA", "V_j_pulsed_V", "V_j_cw_V",
            "mu_pulsed", "eta_inj_pulsed", "eta_inj_cw", "S_retention_pulsed",
            "S_retention_cw", "b_e_window_pulsed", "b_e_window_cw", "t_x_pulsed",
            "eps_pulsed", "t_x_cw", "eps_cw", "edge_beta", "edge_eta_total",
            "collected_flux_pulsed_s", "collected_flux_cw_s", "g2_pulsed",
            "g2_cw0", "g2_cw0_raw", "eligible_pulsed", "eligible_cw",
            "eligible_row", "headline_pass", "secondary_pass", "invalid_reasons_pulsed",
            "invalid_reasons_cw", "assumptions", "pulse_width_ns", "rep_rate_hz",
            "duty_pulsed"]


# -------------------------------------------------------------- statistics

def _finite_stats(values: list) -> tuple:
    arr = np.array([v for v in values if v is not None and np.isfinite(v)], dtype=float)
    if arr.size == 0:
        return float("nan"), float("nan")
    return float(np.min(arr)), float(np.median(arr))


def compute_stats(rows: list) -> dict:
    """Per-card and pooled statistics, computed directly from the row list
    (the same computation the verifier redoes directly from the written
    CSV). Coverage denominator is every scheduled row (spec: "invalid rows
    count as nonpassing"), for the headline metric and for each secondary
    (CW) diagnostic independently -- none of the three gates the others."""
    by_card: dict = {}
    for row in rows:
        by_card.setdefault(row["card_id"], []).append(row)
    per_card = {}
    for card_id, card_rows in by_card.items():
        eligible_rows = [r for r in card_rows if r["eligible_row"]]
        pulsed_vals = [r["g2_pulsed"] for r in eligible_rows]
        cw_raw_vals = [r["g2_cw0_raw"] for r in eligible_rows]
        cw0_vals = [r["g2_cw0"] for r in eligible_rows]
        headline_rows = [r for r in card_rows if r["headline_pass"]]
        p_min, p_med = _finite_stats(pulsed_vals)
        r_min, r_med = _finite_stats(cw_raw_vals)
        c_min, c_med = _finite_stats(cw0_vals)
        per_card[card_id] = {
            "card_class": card_rows[0]["card_class"], "n_rows": len(card_rows),
            "n_eligible": len(eligible_rows),
            "g2_pulsed_min": p_min, "g2_pulsed_median": p_med,
            "g2_cw0_min": c_min, "g2_cw0_median": c_med,
            "g2_cw0_raw_min": r_min, "g2_cw0_raw_median": r_med,
            "n_favorable": len(headline_rows), "favorable_rows": headline_rows,
        }
    all_eligible = [r for r in rows if r["eligible_row"]]
    pooled_min, pooled_median = _finite_stats([r["g2_pulsed"] for r in all_eligible])
    n_total = len(rows)
    n_eligible = len(all_eligible)
    n_headline = sum(1 for r in rows if r["headline_pass"])
    n_cw0_pass = sum(1 for r in rows if r["eligible_row"]
                     and np.isfinite(r["g2_cw0"]) and r["g2_cw0"] < G2_THRESHOLD)
    n_cw_raw_pass = sum(1 for r in rows if r["eligible_row"]
                        and np.isfinite(r["g2_cw0_raw"]) and r["g2_cw0_raw"] < G2_THRESHOLD)
    return {
        "per_card": per_card, "n_total": n_total, "n_eligible": n_eligible,
        "eligible_coverage": (n_eligible / n_total) if n_total else 0.0,
        "n_headline": n_headline,
        "headline_coverage": (n_headline / n_total) if n_total else 0.0,
        "n_cw0_pass": n_cw0_pass,
        "cw0_coverage": (n_cw0_pass / n_total) if n_total else 0.0,
        "n_cw_raw_pass": n_cw_raw_pass,
        "cw_raw_coverage": (n_cw_raw_pass / n_total) if n_total else 0.0,
        "g2_pulsed_min": pooled_min, "g2_pulsed_median": pooled_median,
    }


def compute_verdict(rows: list, stats: dict, grid_complete: bool,
                    evidence_report: dict, hallucination_report: dict) -> dict:
    """Pure policy function over already-computed rows/stats plus the two
    freshly-run paper-check reports -- no card loading, no evaluate() calls,
    so verify_rt_edge_sweep.py can drive it with synthetic fixtures.

    Contract rule (docs/rt_edge_contract.md Acceptance gates, council review
    2026-09-05): the headline metric is pulsed intrinsic g2(0) ALONE.
    g2_cw0/g2_cw0_raw are secondary diagnostics -- reported as coverage
    fractions and per-row, but they never gate PASS. This replaces the prior
    three-gate (pulsed AND g2_cw0 AND g2_cw0_raw) rule, which the contract
    never stated and which made the raw-CW gate (structurally >= ~0.89 for
    every corner) force coverage=0 regardless of the headline metric."""
    evidence_complete = bool(evidence_report.get("evidence_complete", False))
    hallucination_ok = bool(hallucination_report.get("hallucination_tests_passed", False))
    headline_rows = [r for r in rows if r["headline_pass"]]
    headline_by_card: dict = {}
    for r in headline_rows:
        headline_by_card.setdefault(r["card_id"], []).append(r)

    fail_reasons = []
    if not grid_complete:
        fail_reasons.append("grid_incomplete")
    if not evidence_complete:
        fail_reasons.append("evidence_incomplete")
    if not hallucination_ok:
        fail_reasons.append("hallucination_checks_failed")
    if stats["n_eligible"] == 0:
        fail_reasons.append("no_eligible_rows")
    elif not headline_rows:
        fail_reasons.append("no_headline_pass")

    passed = len(fail_reasons) == 0
    g2_min, g2_median = stats["g2_pulsed_min"], stats["g2_pulsed_median"]
    median_pass = bool(np.isfinite(g2_median) and g2_median < G2_THRESHOLD)

    assumptions_used = sorted(set(
        a for r in headline_rows for a in r["assumptions"].split("; ") if a))
    conditional = bool(passed and assumptions_used)

    fallback_only = bool(passed and headline_by_card
                         and all(r["card_class"] != "primary"
                                 for rs in headline_by_card.values() for r in rs))
    primary_id = next((c["id"] for c in CARDS if c["class"] == "primary"), None)
    fallback_id = next((c["id"] for c in CARDS if c["class"] == "fallback"), None)
    note = ""
    if fallback_only and fallback_id:
        note = (f"only the fallback card ({fallback_id}) reaches the acceptance "
                f"gate; the primary card ({primary_id}) does not -- PASS requires "
                f"the material switch to the fallback stack")
    elif passed and not median_pass:
        note = "the best corner passes but the pooled pulsed median does not"

    return {
        "pass": passed, "fail_reasons": fail_reasons,
        "g2_min": g2_min, "g2_median": g2_median, "median_pass": median_pass,
        "coverage": stats["headline_coverage"],
        "headline_coverage_n": stats["n_headline"], "headline_coverage_total": stats["n_total"],
        "cw0_coverage_n": stats["n_cw0_pass"], "cw0_coverage_total": stats["n_total"],
        "cw_raw_coverage_n": stats["n_cw_raw_pass"], "cw_raw_coverage_total": stats["n_total"],
        "eligible_n": stats["n_eligible"], "eligible_total": stats["n_total"],
        "evidence_complete": evidence_complete, "hallucination_ok": hallucination_ok,
        "conditional": conditional, "assumptions_used": assumptions_used,
        "headline_rows": headline_rows, "headline_by_card": headline_by_card,
        "fallback_only": fallback_only, "note": note,
    }


def verdict_line(verdict: dict) -> str:
    return (f"VERDICT: {'PASS' if verdict['pass'] else 'FAIL'} "
            f"g2_min={verdict['g2_min']:.4g} g2_median={verdict['g2_median']:.4g} "
            f"median_pass={'true' if verdict['median_pass'] else 'false'} "
            f"coverage={verdict['coverage']:.4g} "
            f"eligible={verdict['eligible_n']}/{verdict['eligible_total']} "
            f"evidence={'complete' if verdict['evidence_complete'] else 'incomplete'} "
            f"conditional={'true' if verdict['conditional'] else 'false'} "
            f"headline_coverage={verdict['headline_coverage_n']}/{verdict['headline_coverage_total']} "
            f"cw_raw_coverage={verdict['cw_raw_coverage_n']}/{verdict['cw_raw_coverage_total']}")


def card_line(card_id: str, card_stats: dict) -> str:
    return (f"CARD: {card_id} role={card_stats['card_class']} "
            f"g2_pulsed_min={card_stats['g2_pulsed_min']:.4g} "
            f"g2_pulsed_median={card_stats['g2_pulsed_median']:.4g} "
            f"g2_cw_raw_min={card_stats['g2_cw0_raw_min']:.4g} "
            f"g2_cw_raw_median={card_stats['g2_cw0_raw_median']:.4g} "
            f"eligible={card_stats['n_eligible']}/{card_stats['n_rows']} "
            f"favorable_rows={card_stats['n_favorable']}")


# ---------------------------------------------------------------- artifacts

def write_csv(rows: list, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=csv_fieldnames())
        w.writeheader()
        w.writerows(rows)


def write_png(rows: list, stats: dict, path: Path) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # dataviz skill categorical palette (references/palette.md): slot 1 blue
    # for the primary card, slot 2 orange for the fallback card; a neutral
    # gray dashed line for the (non-series) 0.5 acceptance threshold.
    colors = {"edge-inp-gaasp-design": "#2a78d6", "edge-inp-gainp-design": "#eb6834"}
    path.parent.mkdir(parents=True, exist_ok=True)
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(11.0, 4.6))

    for card in CARDS:
        card_id = card["id"]
        card_rows = [r for r in rows if r["card_id"] == card_id]
        if not card_rows:
            continue
        color = colors.get(card_id, "#4a3aa7")
        deltas = sorted({r["delta_xx_meV"] for r in card_rows})
        lo, hi = [], []
        for dxx in deltas:
            vals = [r["g2_pulsed"] for r in card_rows
                    if r["delta_xx_meV"] == dxx and r["eligible_row"] and np.isfinite(r["g2_pulsed"])]
            lo.append(min(vals) if vals else np.nan)
            hi.append(max(vals) if vals else np.nan)
        ax1.plot(deltas, lo, color=color, lw=2, marker="o", ms=5,
                 label=f"{card_id} ({card['class']})")
        ax1.fill_between(deltas, lo, hi, color=color, alpha=0.18, linewidth=0)

        irfs = sorted({r["irf_ps"] for r in card_rows})
        raw_lo, raw_hi, dot_lo, dot_hi = [], [], [], []
        for irf in irfs:
            raw_vals = [r["g2_cw0_raw"] for r in card_rows
                        if r["irf_ps"] == irf and r["eligible_row"] and np.isfinite(r["g2_cw0_raw"])]
            dot_vals = [r["g2_cw0"] for r in card_rows
                        if r["irf_ps"] == irf and r["eligible_row"] and np.isfinite(r["g2_cw0"])]
            raw_lo.append(min(raw_vals) if raw_vals else np.nan)
            raw_hi.append(max(raw_vals) if raw_vals else np.nan)
            dot_lo.append(min(dot_vals) if dot_vals else np.nan)
            dot_hi.append(max(dot_vals) if dot_vals else np.nan)
        ax2.plot(irfs, raw_lo, color=color, lw=2, marker="s", ms=5, linestyle="-",
                 label=f"{card_id} raw (IRF-convolved)")
        ax2.fill_between(irfs, raw_lo, raw_hi, color=color, alpha=0.18, linewidth=0)
        ax2.plot(irfs, dot_lo, color=color, lw=1.4, marker="^", ms=4, linestyle="--",
                 label=f"{card_id} intrinsic (g2_cw0)")

    for ax, xlabel, title in (
        (ax1, "dot.delta_xx (meV)", "Pulsed intrinsic g2(0) envelope"),
        (ax2, "IRF FWHM (ps)", "CW g2(0) envelope: intrinsic (dashed) vs raw (solid)"),
    ):
        ax.axhline(G2_THRESHOLD, color="#52514e", linestyle=":", lw=1.5,
                   label="g2 = 0.5 threshold" if ax is ax1 else None)
        ax.set_xlabel(xlabel)
        ax.set_ylabel("g2(0)")
        ax.set_title(title, fontsize=10)
        ax.set_ylim(bottom=0.0)
        ax.legend(fontsize=6.5, loc="best")

    n_total, n_eligible = stats["n_total"], stats["n_eligible"]
    invalid_frac = 1.0 - (n_eligible / n_total if n_total else 0.0)
    fig.suptitle(
        f"RT edge-emitter acceptance sweep, T_hs=300 K -- "
        f"{n_eligible}/{n_total} rows eligible (invalid coverage {invalid_frac:.1%})",
        fontsize=10)
    fig.tight_layout(rect=(0.0, 0.0, 1.0, 0.94))
    fig.savefig(path, dpi=150)
    plt.close(fig)


def write_markdown(rows: list, stats: dict, verdict: dict, grid: dict,
                   grid_complete: bool, evidence_report: dict,
                   hallucination_report: dict, path: Path) -> None:
    lines = []
    lines.append("# RT edge-emitter acceptance sweep verdict")
    lines.append("")
    lines.append(f"Generated {datetime.now(timezone.utc).isoformat()}; "
                 f"contract: docs/rt_edge_contract.md.")
    lines.append("")
    lines.append(f"```\n{verdict_line(verdict)}\n```")
    for card_id, cs in stats["per_card"].items():
        lines.append(f"```\n{card_line(card_id, cs)}\n```")
    lines.append("")
    lines.append("## Literature ceiling")
    lines.append(
        "No electrically driven III-V single quantum dot g2(0) at 300 K has "
        "been published (docs/rt_edge_contract.md Evidence section; "
        "verify/data/rt_edge_anchors.yaml) -- there is no existing electrical "
        "300 K baseline for this sweep to exceed. The best reported "
        "electrically driven single-dot result at any temperature is "
        "Reischle et al. 2008 at 80 K: g2(0) = 0.43 raw, 0.03 after "
        "background correction. The best reported 300 K single-dot values "
        "are optically pumped: g2(0) ~ 0.5-0.57 (Laferriere et al. 2023, "
        "InAsP/InP nanowire dot, g2(0) = 0.57 at 300 K). This sweep's pooled "
        f"pulsed g2_min={verdict['g2_min']:.4g}, "
        f"g2_median={verdict['g2_median']:.4g} is reported against that "
        "electrical-vs-optical literature picture, not as a claim of an "
        "existing electrical 300 K result to exceed.")
    lines.append("")
    lines.append("## Grid")
    lines.append(f"Grid complete: {grid_complete}"
                 + ("" if grid_complete else " (--quick: endpoints-only, cannot grant PASS)"))
    for path_key, (lo, hi, unit) in RANGE_BOUNDS.items():
        vals = ", ".join(f"{v:g}" for v in grid[path_key])
        lines.append(f"- `{path_key}` in [{lo}, {hi}] {unit}: sampled at {vals}")
    lines.append("")
    lines.append("## Coverage")
    lines.append(f"- eligible coverage: {stats['n_eligible']}/{stats['n_total']} "
                 f"= {stats['eligible_coverage']:.3f}")
    lines.append(f"- **headline coverage** (pulsed intrinsic g2(0) < 0.5, eligible rows -- "
                 f"the contract's PASS metric): "
                 f"{stats['n_headline']}/{stats['n_total']} = {stats['headline_coverage']:.3f}")
    lines.append(f"- secondary coverage, CW intrinsic g2_cw0 < 0.5 (diagnostic only, does "
                 f"not gate PASS): {stats['n_cw0_pass']}/{stats['n_total']} "
                 f"= {stats['cw0_coverage']:.3f}")
    lines.append(f"- secondary coverage, CW IRF-convolved g2_cw0_raw < 0.5 (diagnostic "
                 f"only, does not gate PASS): {stats['n_cw_raw_pass']}/{stats['n_total']} "
                 f"= {stats['cw_raw_coverage']:.3f}")
    lines.append("")
    lines.append("## Per-card statistics")
    lines.append("| card | role | g2_pulsed min/median | g2_cw0 min/median | "
                 "g2_cw0_raw min/median | eligible | headline rows |")
    lines.append("|---|---|---|---|---|---|---|")
    for card_id, cs in stats["per_card"].items():
        lines.append(
            f"| {card_id} | {cs['card_class']} | "
            f"{cs['g2_pulsed_min']:.4g} / {cs['g2_pulsed_median']:.4g} | "
            f"{cs['g2_cw0_min']:.4g} / {cs['g2_cw0_median']:.4g} | "
            f"{cs['g2_cw0_raw_min']:.4g} / {cs['g2_cw0_raw_median']:.4g} | "
            f"{cs['n_eligible']}/{cs['n_rows']} | {cs['n_favorable']} |")
    lines.append("")
    lines.append("## Assumptions required by any headline-passing corner")
    if verdict["assumptions_used"]:
        for a in verdict["assumptions_used"]:
            lines.append(f"- `{a}`")
    else:
        lines.append("- (no row satisfies the headline gate; see fail_reasons below)")
    lines.append("")
    if verdict["note"]:
        lines.append(f"**Note:** {verdict['note']}")
        lines.append("")
    lines.append("## Evidence gate")
    lines.append(f"- evidence_complete: {evidence_report.get('evidence_complete')}")
    lines.append(f"- hallucination_tests_passed (self-test): "
                 f"{hallucination_report.get('hallucination_tests_passed')}")
    missing = evidence_report.get("missing_evidence", [])
    if missing:
        lines.append("- missing/incomplete evidence claims:")
        for m in missing:
            lines.append(f"  - `{m['claim']}`: {m['reason']}")
    lines.append("")
    lines.append("## Fail reasons" if not verdict["pass"] else "## Pass basis")
    if verdict["pass"]:
        lines.append(f"At least one eligible corner has pulsed intrinsic g2(0) < 0.5 -- "
                     f"the contract's headline metric ({stats['n_headline']} such row(s) "
                     f"across both cards), grid complete, evidence complete, hallucination "
                     f"self-test passed. g2_cw0 and g2_cw0_raw are reported above as "
                     f"secondary diagnostics (see Coverage) and do not gate this PASS.")
    else:
        for r in verdict["fail_reasons"]:
            lines.append(f"- `{r}`")
    lines.append("")
    lines.append("## Sources")
    lines.append("- docs/rt_edge_contract.md (acceptance gates, evidence status)")
    lines.append("- verify/data/rt_edge_anchors.yaml (literature evidence ledger)")
    lines.append("- cards/edge-inp-gaasp-design.yaml, cards/edge-inp-gainp-design.yaml "
                 "(design.provenance.sources for every scalar)")
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_manifest(rows: list, stats: dict, verdict: dict, grid: dict,
                   grid_complete: bool, evidence_report: dict,
                   hallucination_report: dict, runtime_s: float, quick: bool,
                   path: Path) -> None:
    manifest = {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "quick": quick, "grid_complete": grid_complete, "grid": grid,
        "range_bounds": RANGE_BOUNDS,
        "pulse_assumption": PULSE_ASSUMPTION,
        "cards": [{"id": c["id"], "class": c["class"], "path": str(c["path"]),
                  "sha256": _sha256_file(c["path"])} for c in CARDS],
        "runtime_versions": {
            "python": platform.python_version(), "numpy": np.__version__,
            "scipy": scipy.__version__, "yaml": getattr(yaml, "__version__", "unknown"),
        },
        "policy": {
            "g2_threshold": G2_THRESHOLD,
            "pass_requires": ["grid_complete", "evidence_complete",
                              "hallucination_tests_passed",
                              "at least one eligible row with pulsed intrinsic g2(0)<0.5 "
                              "(the headline metric; g2_cw0 and g2_cw0_raw are secondary "
                              "diagnostics reported as coverage fractions but do not gate PASS)"],
            "coverage_denominator": "every scheduled row (invalid rows count as nonpassing)",
        },
        "output_schema": {"sweep.csv": csv_fieldnames()},
        "stats": {k: v for k, v in stats.items() if k != "per_card"},
        "per_card_stats": {cid: {k: v for k, v in cs.items() if k not in ("favorable_rows",)}
                          for cid, cs in stats["per_card"].items()},
        "verdict": {k: v for k, v in verdict.items()
                   if k not in ("headline_rows", "headline_by_card")},
        "evidence": {k: evidence_report.get(k) for k in
                    ("evidence_complete", "source_ledger_sha256", "evaluator_input_hashes",
                     "missing_evidence")},
        "hallucination_self_test": {"hallucination_tests_passed":
                                    hallucination_report.get("hallucination_tests_passed")},
        "runtime_seconds": runtime_s,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2, sort_keys=False, default=str), encoding="utf-8")


# ------------------------------------------------------------------- main

def run_sweep(quick: bool) -> tuple:
    grid = build_grid(quick)
    pulsed_cache: dict = {}
    rows = []
    for card in CARDS:
        rows.extend(sweep_card(card, grid, pulsed_cache))
    stats = compute_stats(rows)
    return rows, stats, grid


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--quick", action="store_true",
                        help="small endpoints-only grid; explicitly incomplete, cannot grant PASS")
    parser.add_argument("--out-dir", default=str(ROOT / "out" / "rt_edge"))
    args = parser.parse_args(argv)
    out_dir = Path(args.out_dir)

    t0 = time.time()
    rows, stats, grid = run_sweep(args.quick)

    # Fresh paper-check invocation (spec: "Invoke paper-check run_checks()
    # freshly"); this file never writes evidence.json (out of scope --
    # verify/verify_rt_edge_papers.py owns that file).
    evidence_report = rt_papers.run_checks(self_test=False)
    hallucination_report = rt_papers.run_checks(self_test=True)

    grid_complete = not args.quick
    verdict = compute_verdict(rows, stats, grid_complete, evidence_report, hallucination_report)
    runtime_s = time.time() - t0

    try:
        write_csv(rows, out_dir / "sweep.csv")
        write_png(rows, stats, out_dir / "envelope.png")
        write_markdown(rows, stats, verdict, grid, grid_complete,
                       evidence_report, hallucination_report, out_dir / "verdict.md")
        write_manifest(rows, stats, verdict, grid, grid_complete, evidence_report,
                       hallucination_report, runtime_s, args.quick, out_dir / "manifest.json")
    except OSError as exc:
        print(f"FAIL: could not write output artifacts: {exc}", file=sys.stderr)
        return 1

    for card_id, cs in stats["per_card"].items():
        print(card_line(card_id, cs))
    print(verdict_line(verdict))
    if verdict["note"]:
        print(f"NOTE: {verdict['note']}")
    print(f"artifacts written to {out_dir} ({runtime_s:.1f}s)")
    return 0  # artifacts generated successfully; VERDICT: FAIL is a valid scientific outcome


if __name__ == "__main__":
    sys.exit(main())
