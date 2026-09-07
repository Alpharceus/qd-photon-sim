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

Fourth council review (2026-09-06) additions: three collection-lever axes --
emission.NA, emission.R_back, emission.L_um -- are now genuine sweep axes,
read from EACH card's own provenance.ranges (per-card, unlike the shared
RANGE_BOUNDS below, since the two cards' collection-lever ranges need not be
identical); a card lacking a declared range for one of these falls back to
that card's own scalar (single value). Because these three axes multiply the
grid size, the shared dot.delta_xx/dot.gamma300/irf_ps axes are sampled at
their declared endpoints only for the full (non-quick) grid too (_FULL_N=2,
down from 3) to keep the whole sweep under the ~15 minute budget -- this
interior-sample reduction (and every axis's resolved sample count) is
recorded in manifest.json's "grid"/"lever_info" sections, per the review's
explicit "reduce interior samples on the other axes if needed and say so in
the manifest" instruction.

Fifth council review (2026-09-06) fixes: (1) card_line() no longer prints the
"diagnostic (below flux floor, not measurable)" phrase next to g2_pulsed_min/
median unless a card has zero eligible rows (those two numbers are ALREADY
the eligible-row statistics). (2) gamma300_pass_max/gamma300_threshold are
now derived from a finer, headline-only GAMMA300_REFINE_MEV sample set (8
points, including the verified Chatzarakis 2023 6.5 meV anchor) rather than
just the 2-point sweep.csv grid endpoints; a bracketed gamma300_threshold
states where the true (unsampled) threshold lies, and a per-card
anchor_check() states explicitly whether the verified anchor passes. (3)
RANGE_BOUNDS is now resolve_range_bounds()-derived from the cards' own
provenance.ranges (cross-checked for agreement across cards) instead of a
hardcoded dict that had drifted from them. (4) `coverage` now means
headline-pass fraction over ELIGIBLE rows (a new `eligible_fraction` key is
eligible/total); `flux_margin` (flux_max/floor, >1 clears the floor)
replaces `flux_shortfall` as the preferred key (kept one release,
deprecated). (5) the literature-ceiling paragraph adds Reischle et al.
2008's like-for-like IRF-deconvolved, background-included g2(0) = 0.25 +/-
0.05 anchor alongside the (still IRF-broadened) 0.43 raw dip. (6) the
front-facet-split "self-check" (which only ever reproduced its own
back-solved inputs, a tautology) is replaced with a genuinely independent
forward recomputation read from fsim_core/waveguide.py's OWN source at run
time (never hardcoded, since that file's facet model may change).

CLI: python scripts/run_rt_edge.py [--quick] [--out-dir PATH]
--quick uses a 2-point (endpoints-only) grid, explicitly marked incomplete;
an incomplete grid can never grant VERDICT: PASS (see compute_verdict).

Exit code: 0 whenever the four artifacts were generated (including a
scientific VERDICT: FAIL -- the machine-readable verdict, not the process
exit code, is authoritative for acceptance); nonzero only for a runtime or
artifact-writing error.

Standalone on import (all evaluation happens under
`if __name__ == "__main__"`) except for resolve_range_bounds()'s own two
small provenance.ranges reads off the cards' YAML files at module load
(computing the module-level RANGE_BOUNDS constant -- no evaluate() calls,
no writes); verify/verify_rt_edge_sweep.py imports the functions below
directly (grid, resolve_lever_grid, build_lever_combos, resolve_device_card,
eval_pulsed_point, eval_cw_point, compute_stats, compute_verdict,
_card_gamma300_pass_max, _gamma300_threshold_bracket, resolve_range_bounds,
_range_bounds_mismatches, refine_gamma300, anchor_check, _loading_term,
_brightness_factor_check, _front_facet_split,
_facet_factor_forward_check, _self_test_passed, write_csv, write_png,
write_markdown, write_manifest) rather than parsing this script's stdout.
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import itertools
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
from fsim_core import waveguide  # noqa: E402
import verify.verify_rt_edge_papers as rt_papers  # noqa: E402

CARDS = [
    {"id": "edge-inp-gaasp-design", "class": "primary",
     "path": ROOT / "cards" / "edge-inp-gaasp-design.yaml"},
    {"id": "edge-inp-gainp-design", "class": "fallback",
     "path": ROOT / "cards" / "edge-inp-gainp-design.yaml"},
]

# Contract ranges (docs/rt_edge_contract.md; exact endpoints checked
# mechanically against every card's own provenance.ranges by
# verify/verify_rt_edge_cards.py REQUIRED_RANGES). Fifth council review
# (2026-09-06) item 3: this used to be a bare hardcoded dict --
# ("dot.delta_xx": (4.0, 7.0)) had drifted from both cards' own declared
# provenance.ranges["dot.delta_xx"] = (4.0, 8.0) and from
# verify_rt_edge_cards.py's own REQUIRED_RANGES (4.0, 8.0), so the sweep was
# silently sampling a narrower delta_xx window than the cards/contract
# actually declare. RANGE_BOUNDS is now RESOLVED from the cards' own
# provenance.ranges at run time (see resolve_range_bounds() below); the
# tuples here are only the FALLBACK used for a path a card does not declare
# (currently unused in practice -- both real cards declare all three).
_RANGE_BOUNDS_FALLBACK = {
    "dot.delta_xx": (4.0, 8.0, "meV"),
    "dot.gamma300": (6.0, 20.0, "meV"),
    "irf_ps": (50.0, 200.0, "ps"),
    "thermal.T_hs": (230.0, 300.0, "K"),
}


def _range_bounds_mismatches(bounds: dict, designs: list) -> list:
    """Pure check, no I/O: for each RANGE_BOUNDS path, verify every design's
    own provenance.ranges[path] (lo, hi) agrees with the already-resolved
    `bounds`. `designs` is a list of (card_id, design) pairs -- `design`
    needs only a `.provenance` dict, so a synthetic fixture (no real card
    file) can drive this in verify/verify_rt_edge_sweep.py. Returns a list
    of human-readable mismatch strings (empty when everything agrees)."""
    mismatches = []
    for card_id, design in designs:
        ranges = (getattr(design, "provenance", None) or {}).get("ranges", {})
        for path_key, (lo, hi, _unit) in bounds.items():
            entry = ranges.get(path_key)
            if entry is None or entry.get("lo") is None or entry.get("hi") is None:
                continue  # this design doesn't declare the axis -- nothing to cross-check
            declared = (float(entry["lo"]), float(entry["hi"]))
            if declared != (float(lo), float(hi)):
                mismatches.append(
                    f"{path_key}: card {card_id!r} declares provenance.ranges = "
                    f"{declared} but the resolved RANGE_BOUNDS is ({lo}, {hi})")
    return mismatches


def resolve_range_bounds(cards: list | None = None) -> dict:
    """Council review round 5, item 3: read {lo, hi, unit} for each shared
    RANGE_BOUNDS path from the FIRST card that declares it (in `cards`
    order) instead of a hardcoded constant; fall back to
    _RANGE_BOUNDS_FALLBACK only for a path no card declares. Then
    cross-check with _range_bounds_mismatches that every OTHER card agrees
    with the resolved bound exactly -- RANGE_BOUNDS is a SHARED axis across
    both cards (unlike the per-card LEVER_PATHS ranges below), so a card
    that declares a different (lo, hi) for the same path is a real
    inconsistency, not a value to silently prefer. Raises ValueError on any
    mismatch rather than resolving to one card's number silently."""
    cards = CARDS if cards is None else cards
    designs = [(c["id"], DeviceDesign.load(c["path"])) for c in cards]
    bounds = {}
    for path_key, (lo_fb, hi_fb, unit_fb) in _RANGE_BOUNDS_FALLBACK.items():
        resolved = None
        for _card_id, design in designs:
            ranges = (design.provenance or {}).get("ranges", {})
            entry = ranges.get(path_key)
            if entry is None or entry.get("lo") is None or entry.get("hi") is None:
                continue
            resolved = (float(entry["lo"]), float(entry["hi"]), entry.get("unit", unit_fb))
            break
        bounds[path_key] = resolved if resolved is not None else (lo_fb, hi_fb, unit_fb)
    mismatches = _range_bounds_mismatches(bounds, designs)
    if mismatches:
        raise ValueError("RANGE_BOUNDS disagree with a card's declared "
                         "provenance.ranges: " + "; ".join(mismatches))
    return bounds


RANGE_BOUNDS = resolve_range_bounds()

# Collection-lever axes (council review round 4, item 1): unlike RANGE_BOUNDS
# above, these are declared per-card in each card's own provenance.ranges
# (the two cards' collection-lever ranges need not agree), so there is no
# shared (lo, hi) tuple here -- resolve_lever_grid() reads each card's own
# declared range (or falls back to that card's own scalar) at run time.
LEVER_PATHS = ["emission.NA", "emission.R_back", "emission.L_um"]

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
# [A] Below roughly 1 kHz collected pulsed flux, a g2 measurement is not
# practical within hours at single-photon-detector counting rates.  This is
# an eligibility/reporting floor only; it is never fed into device physics.
FLUX_FLOOR_PULSED_S = 1.0e3
_ENDPOINT_N = 2
# Reduced from 3 to 2 (council review round 4, item 1): adding the three
# emission.NA/R_back/L_um lever axes below multiplies the grid size by up to
# ~12 (this repo's two cards' own declared ranges); keeping the prior
# endpoints+interior (n=3) sampling on dot.delta_xx/dot.gamma300/irf_ps as
# well would push the full sweep to ~20+ minutes. Endpoints-only (n=2) on
# these three axes keeps the full sweep under the ~15 minute budget; this
# is the documented "reduce interior samples on the other axes" trade-off
# the review explicitly allows, recorded in manifest.json's grid section.
_FULL_N = 2

# Council review round 5, item 2: the main grid's dot.gamma300 axis is only
# the two RANGE_BOUNDS endpoints ({6, 20} meV) -- gamma300_pass_max derived
# from it can only ever equal whichever endpoint happens to pass, which is
# NOT the same thing as "the largest linewidth at which the card passes"
# (there could be passing/failing samples in between). This finer,
# headline-only sample set (used ONLY for the gamma300_pass_max /
# gamma300_threshold determination below -- never written to sweep.csv, so
# it does not change the CSV/PNG grid-completeness contract) includes the
# verified Chatzarakis et al., Phys. Rev. Applied 20, 034011 (2023) 6.5 meV
# lower-anchor value. irf_ps is intentionally NOT swept here: g2_op (the
# headline metric) is irf_ps-independent (eval_pulsed_point's own pulsed
# sub-result is cached without it), so no CW/irf evaluation is performed by
# refine_gamma300() at all -- this is the "reduce irf_ps to the two
# endpoints if runtime needs it" trade-off the review allows, taken to its
# natural conclusion (irf_ps contributes nothing to this determination).
GAMMA300_REFINE_MEV = [6.0, 6.5, 7.0, 8.0, 10.0, 12.0, 16.0, 20.0]

# Council review round 5, item 2: the verified Chatzarakis 2023 300 K
# linewidth lower-anchor value (chatzarakis23-gamma300-class; also the
# lo endpoint of RANGE_BOUNDS["dot.gamma300"]) -- report explicitly, per
# card, whether it passes the headline gate at the card's OWN default
# delta_xx (never the sweep's endpoint grid), since both real cards'
# provenance.ranges["dot.gamma300"].note already document a fresh
# evaluation at exactly this value.
CHATZARAKIS_ANCHOR_GAMMA300_MEV = 6.5

# Council review round 5, item 5: Reischle et al., Optics Express 16, 12771
# (2008) reports g2(0) on THREE different conventions (QD C, 80 K); the
# ceiling paragraph below used to compare this sweep's IRF-free intrinsic
# g2_min against Reischle's RAW dip (0.43, still IRF-broadened) -- not a
# like-for-like comparison. The IRF-DECONVOLVED-but-background-included
# value is the correct like-for-like anchor against this sweep's g2_op
# (also background-included via drive.b_res/rho, also never IRF-convolved
# for the pulsed metric): g2_b(0) = 0.25 +/- 0.05 (QD C, 80 K)
# (../_goal/paper_digests.md line ~36).
#
# Council review round 6, item 5: these three numbers (and the background-
# assumption citation below, item 1) are now READ from the ledger
# (verify/data/rt_edge_anchors.yaml) instead of being hardcoded, so they
# cannot drift from their anchors again: RAW/BG_RES from the single-source
# raw-dip/residual-background anchors, DECONV/its error from the deconvolved
# anchor added for this review.
_RAW_G2_ANCHOR_ID = "reischle08-g2-80k"
_DECONV_G2_ANCHOR_ID = "reischle08-g2-80k-deconvolved"
BG_RES_ANCHOR_ID = "reischle08-b-res-80k"
_rt_anchors = rt_papers.load_anchors()
REISCHLE_RAW_G2 = float(_rt_anchors[_RAW_G2_ANCHOR_ID]["value"])
REISCHLE_DECONV_G2 = float(_rt_anchors[_DECONV_G2_ANCHOR_ID]["value"])
REISCHLE_DECONV_G2_ERR = rt_papers.anchor_bound(_rt_anchors[_DECONV_G2_ANCHOR_ID])


# --------------------------------------------------------------------- grid

def build_grid(quick: bool) -> dict:
    """{"dot.delta_xx": [...], "dot.gamma300": [...], "irf_ps": [...]}, always
    including both declared endpoints; quick=True is endpoints-only (spec:
    "explicitly marked incomplete and cannot grant scientific PASS"). The
    full grid was endpoints-plus-one-interior-sample (n=3) before council
    review round 4; it is now endpoints-only (n=2, same as --quick) so the
    three new emission.NA/R_back/L_um lever axes below (up to ~12 combos
    per card) fit the ~15 minute runtime budget -- see _FULL_N above."""
    n = _ENDPOINT_N if quick else _FULL_N
    grid = {path: [float(v) for v in np.linspace(lo, hi, n)]
            for path, (lo, hi, _unit) in RANGE_BOUNDS.items()}
    lo, hi, _unit = RANGE_BOUNDS["thermal.T_hs"]
    grid["thermal.T_hs"] = ([float(lo), 250.0, 273.0, float(hi)]
                             if not quick else [float(lo), float(hi)])
    return grid


def resolve_lever_grid(design0: DeviceDesign, quick: bool) -> dict:
    """Per-card sample set for each of the three collection-lever axes
    (LEVER_PATHS): the card's own declared provenance.ranges endpoints plus
    the card's own scalar value, deduplicated (spec: "sampled at the range
    ends plus the card value"); a card with no declared range for one of
    these paths falls back to that card's single declared scalar (== the
    device.py class default whenever the card does not override it) --
    this file never invents a range the card does not state. --quick
    collapses every lever to the card's own scalar only (a single point),
    matching --quick's existing endpoints-only-but-smaller convention for
    the other axes (still explicitly incomplete; grid_complete governs
    PASS eligibility, not this axis's density)."""
    ranges = (design0.provenance or {}).get("ranges", {})
    grid = {}
    for path in LEVER_PATHS:
        block_name, field_name = path.split(".", 1)
        card_value = getattr(getattr(design0, block_name), field_name)
        if quick:
            grid[path] = [card_value]
            continue
        rng = ranges.get(path)
        if not rng or rng.get("lo") is None or rng.get("hi") is None:
            grid[path] = [card_value]
            continue
        candidates = [float(rng["lo"]), float(rng["hi"])]
        if card_value is not None:
            candidates.append(float(card_value))
        deduped = []
        for v in candidates:
            if not any(abs(v - seen) < 1e-12 for seen in deduped):
                deduped.append(v)
        grid[path] = deduped
    return grid


def build_lever_combos(lever_grid: dict) -> list:
    """Full cartesian product of the per-card lever grid -- LEVER_PATHS are
    genuine sweep axes, exactly like dot.delta_xx/dot.gamma300/irf_ps, just
    resolved per-card rather than from the shared RANGE_BOUNDS. Each entry
    is an overrides dict ({"emission.NA": v, "emission.R_back": v,
    "emission.L_um": v}) directly consumable by resolve_device_card."""
    return [dict(zip(LEVER_PATHS, combo))
            for combo in itertools.product(*(lever_grid[p] for p in LEVER_PATHS))]


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
                      lever: dict, cache: dict | None = None,
                      T_hs: float | None = None) -> dict:
    """Pulsed sub-result at (delta_xx, gamma300, lever); irf_ps-independent,
    so cache is keyed without it and reused across the irf_ps sweep axis.
    g2_op itself does not depend on the emission.NA/R_back/L_um collection
    levers (device.py's edge-out-coupling factor only multiplies into
    brightness_per_pulse, never into op["g2"]) but brightness/collected
    flux do, so the lever values are still part of the cache key -- a
    lever-varying row cannot reuse another lever's cached brightness."""
    base_T_hs = DeviceDesign.load(card_path).thermal.T_hs if T_hs is None else T_hs
    key = (str(card_path), "pulsed", delta_xx, gamma300, base_T_hs,
          lever["emission.NA"], lever["emission.R_back"], lever["emission.L_um"])
    if cache is not None and key in cache:
        return cache[key]
    base = DeviceDesign.load(card_path)
    duty = PULSE_WIDTH_NS * 1e-9 * REP_RATE_HZ
    overrides = {
        "dot.delta_xx": delta_xx, "dot.gamma300": gamma300, "thermal.T_hs": base_T_hs,
        "drive.duty": duty, "drive.cw": False,
        "drive.diode": {**base.drive.diode, "tau_pulse_ns": PULSE_WIDTH_NS},
    }
    overrides.update(lever)
    design = resolve_device_card(card_path, overrides)
    sc = evaluate(design, T_grid=[design.thermal.T_hs])["scalars"]
    flux = _collected_flux_s(sc)
    eligible, reasons = _classify(sc, flux, "pulsed")
    result = {"scalars": sc, "flux_s": flux, "eligible": eligible,
              "reasons": reasons, "duty": duty}
    if cache is not None:
        cache[key] = result
    return result


def eval_cw_point(card_path: Path, delta_xx: float, gamma300: float,
                  irf_ps: float, lever: dict, T_hs: float | None = None) -> dict:
    """CW (duty=1 DC) sub-result; genuinely irf_ps-dependent (g2_cw0_raw),
    so every irf_ps grid value gets its own evaluate() call."""
    overrides = {
        "dot.delta_xx": delta_xx, "dot.gamma300": gamma300,
        "thermal.T_hs": (DeviceDesign.load(card_path).thermal.T_hs if T_hs is None else T_hs),
        "drive.duty": 1.0, "drive.cw": True, "drive.cw_irf_fwhm_ps": irf_ps,
    }
    overrides.update(lever)
    design = resolve_device_card(card_path, overrides)
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


def _flux_measurable(sc: dict, flux_s: float) -> bool:
    """Use device.py's measurability flag when available; otherwise apply
    the sweep's [A] 1 kHz pulsed collected-flux floor."""
    if "flux_measurable" in sc:
        return bool(sc["flux_measurable"])
    return bool(np.isfinite(flux_s) and flux_s >= FLUX_FLOOR_PULSED_S)


def _classify(sc: dict, flux_s: float, role: str) -> tuple[bool, list]:
    """Eligibility per docs/rt_edge_contract.md / the spec's interface
    constraints: supported confinement/transport/thermal/optical results,
    finite diode V/I with positive operating injection, eta_inj>0 and
    collected_flux_s>0. device.py's own invalid_reasons already covers
    thermal nonconvergence/runaway, zero/negative current, and absorbing/
    unguided/unsupported emission stacks (edge_err); the checks below are
    additional. Pulsed rows must also meet the [A] 1 kHz collected-flux floor;
    rows below it remain in the CSV but are ineligible."""
    reasons = list(sc.get("invalid_reasons", []))
    if not np.isfinite(sc.get("T_j_op", float("nan"))):
        reasons.append(f"{role}: non-finite T_j_op")
    if not np.isfinite(sc.get("V_j_op", float("nan"))):
        reasons.append(f"{role}: non-finite diode V_j_op")
    if not (sc.get("eta_inj", float("nan")) > 0):
        reasons.append(f"{role}: eta_inj not positive")
    if not (np.isfinite(flux_s) and flux_s > 0):
        reasons.append(f"{role}: collected_flux_s not positive")
    if role == "pulsed" and not _flux_measurable(sc, flux_s):
        reasons.append("ineligible: flux_below_floor")
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

def sweep_card(card: dict, grid: dict, pulsed_cache: dict, quick: bool) -> tuple:
    """One card's full cartesian (delta_xx, gamma300, irf_ps, lever) grid ->
    rows, plus the resolved per-card lever grid/combo count for manifest
    reporting. Every scheduled row is emitted, eligible or not (spec: never
    drop invalid rows)."""
    card_path = card["path"]
    design0 = DeviceDesign.load(card_path)
    provenance = design0.provenance or {}
    ranges = provenance.get("ranges", {})
    card_assumptions = list(provenance.get("assumptions", []))
    i_ua = design0.drive.I_uA
    alpha_cm = design0.emission.alpha_cm
    lever_grid = resolve_lever_grid(design0, quick)
    combos = build_lever_combos(lever_grid)
    rows = []
    for lever in combos:
        for T_hs in grid["thermal.T_hs"]:
            for delta_xx in grid["dot.delta_xx"]:
                for gamma300 in grid["dot.gamma300"]:
                    pulsed = eval_pulsed_point(card_path, delta_xx, gamma300, lever, pulsed_cache, T_hs)
                    for irf_ps in grid["irf_ps"]:
                        cw = eval_cw_point(card_path, delta_xx, gamma300, irf_ps, lever, T_hs)
                        rows.append(_build_row(card, ranges, card_assumptions, i_ua,
                                               delta_xx, gamma300, irf_ps, T_hs, lever, pulsed, cw,
                                               alpha_cm))
    return rows, lever_grid, len(combos)


def refine_gamma300(card: dict, delta_xx_grid: list, lever_combos: list,
                    pulsed_cache: dict, T_hs_values: list) -> list:
    """Council review round 5, item 2: headline_pass at every (delta_xx,
    lever-combo) x GAMMA300_REFINE_MEV sample -- "keep the lever axes"
    (every lever combo is still searched, since eligibility, via the
    collected-flux floor, DOES depend on the lever even though g2_op does
    not) while never evaluating a CW/irf_ps sub-point at all (headline_pass
    depends only on the pulsed sub-result, which eval_pulsed_point already
    caches independently of irf_ps -- so this is strictly cheaper than
    adding gamma300 samples to the main grid, which would also multiply the
    CW/irf axis). Reuses `pulsed_cache` with the SAME (card, delta_xx,
    gamma300, lever) key scheme as the main sweep, so gamma300 values
    already in RANGE_BOUNDS's endpoints (computed by sweep_card for the
    same delta_xx/lever) are cache hits, not recomputed. Returns a list of
    row-like dicts (never written to sweep.csv) with just the fields
    _card_gamma300_pass_max / _gamma300_threshold_bracket need -- including
    `eligible` (council review round 6, item 2), which _gamma300_threshold_
    bracket uses to classify WHY no sample passes when that happens."""
    card_path = card["path"]
    design0 = DeviceDesign.load(card_path)
    provenance = design0.provenance or {}
    assumptions = "; ".join(sorted(set(provenance.get("assumptions", []))
                                   | set(SWEEP_ASSUMPTION_KEYS)))
    rows = []
    for T_hs in T_hs_values:
      for lever in lever_combos:
        for delta_xx in delta_xx_grid:
            for gamma300 in GAMMA300_REFINE_MEV:
                pulsed = eval_pulsed_point(card_path, delta_xx, gamma300, lever, pulsed_cache, T_hs)
                g2_p = _f_or_none(pulsed["scalars"].get("g2_op"))
                headline_pass = bool(pulsed["eligible"] and g2_p is not None
                                     and g2_p < G2_THRESHOLD)
                rows.append({
                    "card_id": card["id"], "headline_pass": headline_pass,
                    "eligible": bool(pulsed["eligible"]),
                    "T_hs_K": T_hs,
                    "gamma300_meV": gamma300, "g2_pulsed": g2_p,
                    "delta_xx_meV": delta_xx, "irf_ps": None,
                    "emission_NA": lever["emission.NA"],
                    "emission_R_back": lever["emission.R_back"],
                    "emission_L_um": lever["emission.L_um"],
                    "assumptions": assumptions,
                })
    return rows


def anchor_check(card: dict, pulsed_cache: dict) -> dict:
    """Council review round 5, item 2: a fresh, single eval_pulsed_point
    call at the verified Chatzarakis 2023 CHATZARAKIS_ANCHOR_GAMMA300_MEV
    (6.5 meV) anchor, using the CARD'S OWN default delta_xx and
    emission.NA/R_back/L_um (no sweep-grid overrides at all beyond
    dot.gamma300) -- both cards' own provenance.ranges["dot.gamma300"].note
    already document this exact fresh evaluation by hand; this reproduces
    it programmatically instead of quoting the card comment's numbers."""
    card_path = card["path"]
    design0 = DeviceDesign.load(card_path)
    lever = {"emission.NA": design0.emission.NA, "emission.R_back": design0.emission.R_back,
            "emission.L_um": design0.emission.L_um}
    delta_xx = design0.dot.delta_xx
    pulsed = eval_pulsed_point(card_path, delta_xx, CHATZARAKIS_ANCHOR_GAMMA300_MEV,
                              lever, pulsed_cache)
    g2 = _f_or_none(pulsed["scalars"].get("g2_op"))
    passes = bool(pulsed["eligible"] and g2 is not None and g2 < G2_THRESHOLD)
    return {"card_id": card["id"], "delta_xx_meV": float(delta_xx),
           "gamma300_meV": CHATZARAKIS_ANCHOR_GAMMA300_MEV, "g2_pulsed": g2,
           "flux_s": _f_or_none(pulsed["flux_s"]), "eligible": pulsed["eligible"],
           "passes": passes}


def _build_row(card, ranges, card_assumptions, i_ua, delta_xx, gamma300, irf_ps, T_hs,
              lever, pulsed, cw, alpha_cm) -> dict:
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
        {"card": card["id"], "delta_xx": delta_xx, "gamma300": gamma300, "T_hs": T_hs,
         "irf_ps": irf_ps, "pulse_width_ns": PULSE_WIDTH_NS, "rep_rate_hz": REP_RATE_HZ,
         "lever": lever},
        sort_keys=True).encode()).hexdigest()[:16]
    return {
        "card_id": card["id"], "card_class": card["class"], "config_id": config_hash,
        "delta_xx_meV": delta_xx, "delta_xx_tag": ranges.get("dot.delta_xx", {}).get("tag", ""),
        "gamma300_meV": gamma300, "gamma300_tag": ranges.get("dot.gamma300", {}).get("tag", ""),
        "irf_ps": irf_ps, "irf_tag": ranges.get("irf_ps", {}).get("tag", ""),
        "T_hs_K": T_hs,
        "Tj_pulsed_K": scp.get("T_j_op"), "Tj_cw_K": sccw.get("T_j_op"),
        "I_uA": i_ua,
        "V_j_pulsed_V": scp.get("V_j_op"), "V_j_cw_V": sccw.get("V_j_op"),
        "mu_pulsed": scp.get("mu_resolved"),
        "eta_inj_pulsed": scp.get("eta_inj"), "eta_inj_cw": sccw.get("eta_inj"),
        "S_retention_pulsed": scp.get("S_resolved"), "S_retention_cw": sccw.get("S_resolved"),
        "Gamma_pulsed_meV": scp.get("gamma_op"),
        "b_e_window_pulsed": scp.get("b_e_resolved"), "b_e_window_cw": sccw.get("b_e_resolved"),
        "t_x_pulsed": scp.get("t_x_op"), "eps_pulsed": scp.get("eps_op"),
        "t_x_cw": sccw.get("t_x_op"), "eps_cw": sccw.get("eps_op"),
        "edge_beta": scp.get("edge_beta"), "edge_eta_total": scp.get("edge_eta_total"),
        "edge_T_facet": scp.get("edge_T_facet"), "edge_eta_prop": scp.get("edge_eta_prop"),
        "edge_eta_NA": scp.get("edge_eta_NA"),
        "collected_flux_pulsed_s": pulsed["flux_s"], "collected_flux_cw_s": cw["flux_s"],
        "g2_pulsed": g2_p, "g2_cw0": g2_cw0, "g2_cw0_raw": g2_cw0_raw,
        # Per-row diagnostic aliases are deliberately written alongside the
        # legacy metric columns.  They retain finite evaluator values for
        # below-floor rows; diagnostic_valid is the separate validity gate.
        "diag_g2_pulsed": g2_p, "diag_g2_cw0": g2_cw0,
        "diag_g2_cw0_raw": g2_cw0_raw,
        "eligible_pulsed": pulsed["eligible"], "eligible_cw": cw["eligible"],
        "eligible_row": eligible_row, "headline_pass": headline_pass,
        "secondary_pass": secondary_pass,
        "diagnostic_valid": _diagnostic_valid(pulsed, cw),
        "invalid_reasons_pulsed": "; ".join(pulsed["reasons"]),
        "invalid_reasons_cw": "; ".join(cw["reasons"]),
        "assumptions": "; ".join(assumptions),
        "pulse_width_ns": PULSE_WIDTH_NS, "rep_rate_hz": REP_RATE_HZ,
        "duty_pulsed": pulsed["duty"],
        # Collection-lever axes (council review round 4, item 1): which
        # emission.NA/R_back/L_um values produced THIS row.
        "emission_NA": lever["emission.NA"], "emission_R_back": lever["emission.R_back"],
        "emission_L_um": lever["emission.L_um"],
        # emission_alpha_cm (peer-review pkg2 fix3, 2026-09-07, item 3): not
        # a swept lever, but recorded per-row from design.emission.alpha_cm
        # so _facet_factor_forward_check's forward recomputation reads the
        # ACTUAL card value instead of a hardcoded literal.
        "emission_alpha_cm": alpha_cm,
    }


def csv_fieldnames() -> list:
    return ["card_id", "card_class", "config_id", "delta_xx_meV", "delta_xx_tag",
            "gamma300_meV", "gamma300_tag", "irf_ps", "irf_tag", "T_hs_K",
            "Tj_pulsed_K", "Tj_cw_K", "I_uA", "V_j_pulsed_V", "V_j_cw_V",
            "mu_pulsed", "eta_inj_pulsed", "eta_inj_cw", "S_retention_pulsed",
            "S_retention_cw", "Gamma_pulsed_meV", "b_e_window_pulsed", "b_e_window_cw", "t_x_pulsed",
            "eps_pulsed", "t_x_cw", "eps_cw", "edge_beta", "edge_eta_total",
            "collected_flux_pulsed_s", "collected_flux_cw_s", "g2_pulsed",
            "g2_cw0", "g2_cw0_raw", "eligible_pulsed", "eligible_cw",
            "eligible_row", "headline_pass", "secondary_pass", "invalid_reasons_pulsed",
            "invalid_reasons_cw", "assumptions", "pulse_width_ns", "rep_rate_hz",
            "duty_pulsed", "edge_T_facet", "edge_eta_prop", "edge_eta_NA",
            "diagnostic_valid", "diag_g2_pulsed", "diag_g2_cw0",
            "diag_g2_cw0_raw", "emission_NA", "emission_R_back", "emission_L_um",
            "emission_alpha_cm"]


# -------------------------------------------------------------- statistics

def _finite_stats(values: list) -> tuple:
    arr = np.array([v for v in values if v is not None and np.isfinite(v)], dtype=float)
    if arr.size == 0:
        return float("nan"), float("nan")
    return float(np.min(arr)), float(np.median(arr))


def _diagnostic_valid(pulsed: dict, cw: dict) -> bool:
    """A row usable for below-floor diagnostics: only the measurement-floor
    exclusion is tolerated; all evaluator/classification invalidity remains
    excluded.  The g2 finiteness check is applied by compute_stats per metric.
    """
    allowed = {"ineligible: flux_below_floor"}
    reasons = set(pulsed.get("reasons", [])) | set(cw.get("reasons", []))
    return reasons <= allowed


def _row_diagnostic_valid(row: dict) -> bool:
    """CSV/synthetic-row counterpart of _diagnostic_valid."""
    reasons = []
    for key in ("invalid_reasons_pulsed", "invalid_reasons_cw"):
        text = row.get(key, "")
        reasons.extend(r.strip() for r in text.split(";") if r.strip())
    return set(reasons) <= {"ineligible: flux_below_floor"}


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
        diagnostic_rows = [r for r in card_rows if _row_diagnostic_valid(r)]
        # Preserve the existing g2_* statistics as eligible-row statistics.
        # The diagnostic_* statistics below intentionally use all valid rows,
        # including rows excluded only by the collected-flux floor.
        eligible_pulsed_vals = [r["g2_pulsed"] for r in eligible_rows]
        eligible_cw_raw_vals = [r["g2_cw0_raw"] for r in eligible_rows]
        eligible_cw0_vals = [r["g2_cw0"] for r in eligible_rows]
        pulsed_vals = [r.get("diag_g2_pulsed", r.get("g2_pulsed"))
                       for r in diagnostic_rows]
        cw_raw_vals = [r.get("diag_g2_cw0_raw", r.get("g2_cw0_raw"))
                       for r in diagnostic_rows]
        cw0_vals = [r.get("diag_g2_cw0", r.get("g2_cw0"))
                    for r in diagnostic_rows]
        headline_rows = [r for r in card_rows if r["headline_pass"]]
        flux_floor_excluded = sum(
            1 for r in card_rows
            if "ineligible: flux_below_floor" in r.get("invalid_reasons_pulsed", ""))
        p_min, p_med = _finite_stats(eligible_pulsed_vals)
        r_min, r_med = _finite_stats(eligible_cw_raw_vals)
        c_min, c_med = _finite_stats(eligible_cw0_vals)
        dp_min, dp_med = _finite_stats(pulsed_vals)
        dr_min, dr_med = _finite_stats(cw_raw_vals)
        dc_min, dc_med = _finite_stats(cw0_vals)
        per_card[card_id] = {
            "card_class": card_rows[0]["card_class"], "n_rows": len(card_rows),
            "n_eligible": len(eligible_rows),
            "g2_pulsed_min": p_min, "g2_pulsed_median": p_med,
            "g2_cw0_min": c_min, "g2_cw0_median": c_med,
            "g2_cw0_raw_min": r_min, "g2_cw0_raw_median": r_med,
            "n_favorable": len(headline_rows), "favorable_rows": headline_rows,
            "n_flux_floor_excluded": flux_floor_excluded,
            "diag_g2_pulsed_min": dp_min, "diag_g2_pulsed_median": dp_med,
            "diag_g2_cw0_min": dc_min, "diag_g2_cw0_median": dc_med,
            "diag_g2_cw0_raw_min": dr_min, "diag_g2_cw0_raw_median": dr_med,
            "flux_max": _max_finite([r.get("collected_flux_pulsed_s") for r in diagnostic_rows]),
            "best_diagnostic_row": _best_diagnostic_row(diagnostic_rows),
        }
    diagnostic_rows = [r for r in rows if _row_diagnostic_valid(r)]
    all_eligible = [r for r in rows if r["eligible_row"]]
    pooled_min, pooled_median = _finite_stats([r["g2_pulsed"] for r in all_eligible])
    n_total = len(rows)
    n_eligible = len(all_eligible)
    n_headline = sum(1 for r in rows if r["headline_pass"])
    n_flux_floor_excluded = sum(
        1 for r in rows
        if "ineligible: flux_below_floor" in r.get("invalid_reasons_pulsed", ""))
    n_cw0_pass = sum(1 for r in rows if r["eligible_row"]
                     and np.isfinite(r["g2_cw0"]) and r["g2_cw0"] < G2_THRESHOLD)
    n_cw_raw_pass = sum(1 for r in rows if r["eligible_row"]
                        and np.isfinite(r["g2_cw0_raw"]) and r["g2_cw0_raw"] < G2_THRESHOLD)
    # Peer review finding 6 / pkg5 fix, item 1: eval_pulsed_point's own
    # docstring establishes that the pulsed sub-result (g2_pulsed) is
    # irf_ps-independent, yet the main grid repeats every (card, delta_xx,
    # gamma300, lever, T_hs) combination once per irf_ps sample -- so
    # headline_coverage (n_headline/n_total, kept below unchanged for
    # backward compatibility) double-counts an axis that cannot change the
    # pulsed sub-result. headline_coverage_pulsed instead counts each such
    # combination once, via the key that excludes irf_ps.
    #
    # eligible_row = eligible_pulsed AND eligible_cw, and cw eligibility's
    # g2_cw0_raw finiteness check IS IRF-convolved, so eligible_row (and
    # therefore headline_pass, which requires eligible_row) is NOT
    # guaranteed irf-independent -- a group's rows can legitimately
    # disagree. The dedup below is therefore a group reduction, not a
    # first-row pick: a corner counts as passing only if EVERY sampled irf
    # value agrees it passes (`all()`), and every group where the irf axis
    # disagrees is counted and surfaced (headline_dedup_mismatch_groups,
    # eligible_dedup_mismatch_groups) rather than silently resolved by
    # picking whichever row happened to be seen first.
    dedup_groups: dict = {}
    for r in rows:
        key = (r.get("card_id"), r.get("delta_xx_meV"), r.get("gamma300_meV"),
              r.get("emission_NA"), r.get("emission_R_back"), r.get("emission_L_um"),
              r.get("T_hs_K"))
        dedup_groups.setdefault(key, []).append(r)
    n_scheduled_dedup = len(dedup_groups)
    n_headline_dedup = 0
    n_headline_mismatch_groups = 0
    n_eligible_dedup = 0
    n_eligible_mismatch_groups = 0
    for grp in dedup_groups.values():
        headline_flags = [bool(r["headline_pass"]) for r in grp]
        if all(headline_flags):
            n_headline_dedup += 1
        if len(set(headline_flags)) > 1:
            n_headline_mismatch_groups += 1
        eligible_flags = [bool(r["eligible_row"]) for r in grp]
        if all(eligible_flags):
            n_eligible_dedup += 1
        if len(set(eligible_flags)) > 1:
            n_eligible_mismatch_groups += 1
    per_T = {}
    for T_hs in (230.0, 250.0, 273.0, 300.0):
        t_rows = [r for r in rows if float(r.get("T_hs_K", float("nan"))) == T_hs]
        t_eligible = [r for r in t_rows if r.get("eligible_row")]
        t_headline = [r for r in t_rows if r.get("headline_pass")]
        per_T[str(int(T_hs))] = {
            "n_total": len(t_rows), "n_eligible": len(t_eligible),
            "n_headline": len(t_headline),
            "g2_min": _finite_stats([r.get("g2_pulsed") for r in t_eligible])[0],
            "flux_max": _max_finite([r.get("collected_flux_pulsed_s") for r in t_rows]),
        }
    # pkg5-fix2, item 3: per_T_pulsed is built from the T_hs_K values
    # actually present in the dedup groups (sorted; a group whose T_hs_K is
    # None -- a synthetic fixture row that never set it -- is kept under a
    # "T=?" bucket, not silently dropped), rather than a hardcoded tuple of
    # the four T_hs set points that would drop any other T_hs sampled by a
    # future grid change. The numerators/denominators are guarded to sum to
    # the pooled headline_coverage_pulsed count so a partition bug fails
    # loudly instead of silently under/over-counting.
    per_T_pulsed = {}
    t_hs_values_present = sorted(
        {key[-1] for key in dedup_groups}, key=lambda v: (v is None, v))
    for T_hs in t_hs_values_present:
        bucket = "T=?" if T_hs is None else str(int(T_hs))
        t_groups = [grp for key, grp in dedup_groups.items() if key[-1] == T_hs]
        t_n_scheduled_dedup = len(t_groups)
        t_n_headline_dedup = sum(1 for grp in t_groups
                                 if all(bool(r["headline_pass"]) for r in grp))
        per_T_pulsed[bucket] = {
            "n_scheduled_dedup": t_n_scheduled_dedup,
            "n_headline_dedup": t_n_headline_dedup,
        }
    _per_T_pulsed_num = sum(v["n_headline_dedup"] for v in per_T_pulsed.values())
    _per_T_pulsed_den = sum(v["n_scheduled_dedup"] for v in per_T_pulsed.values())
    if _per_T_pulsed_num != n_headline_dedup or _per_T_pulsed_den != n_scheduled_dedup:
        raise ValueError(
            f"per_T_pulsed buckets sum to {_per_T_pulsed_num}/{_per_T_pulsed_den} but the "
            f"pooled headline_coverage_pulsed is {n_headline_dedup}/{n_scheduled_dedup} -- "
            "the per-T partition of dedup_groups does not account for every group")
    return {
        "per_card": per_card, "n_total": n_total, "n_eligible": n_eligible,
        "eligible_coverage": (n_eligible / n_total) if n_total else 0.0,
        "n_headline": n_headline,
        "headline_coverage": (n_headline / n_total) if n_total else 0.0,
        "n_headline_dedup": n_headline_dedup, "n_scheduled_dedup": n_scheduled_dedup,
        "headline_coverage_pulsed": ((n_headline_dedup / n_scheduled_dedup)
                                     if n_scheduled_dedup else 0.0),
        "headline_dedup_mismatch_groups": n_headline_mismatch_groups,
        "n_eligible_dedup": n_eligible_dedup,
        # pkg5-fix2, item 2: named "eligible_dedup" (not "eligible_pulsed",
        # pkg5's original name), which collided with the pre-existing
        # per-row CSV column of the same name (a bool: was the pulsed
        # sub-path eligible for THIS row).
        "eligible_dedup": ((n_eligible_dedup / n_scheduled_dedup)
                           if n_scheduled_dedup else 0.0),
        "eligible_dedup_mismatch_groups": n_eligible_mismatch_groups,
        "per_T_pulsed": per_T_pulsed,
        "n_cw0_pass": n_cw0_pass,
        "cw0_coverage": (n_cw0_pass / n_eligible) if n_eligible else 0.0,
        "n_cw_raw_pass": n_cw_raw_pass,
        "cw_raw_coverage": (n_cw_raw_pass / n_eligible) if n_eligible else 0.0,
        "per_T": per_T,
        "n_flux_floor_excluded": n_flux_floor_excluded,
        "g2_pulsed_min": pooled_min, "g2_pulsed_median": pooled_median,
        "diag_g2_pulsed_min": _finite_stats([r.get("diag_g2_pulsed", r.get("g2_pulsed")) for r in diagnostic_rows])[0],
        "diag_g2_pulsed_median": _finite_stats([r.get("diag_g2_pulsed", r.get("g2_pulsed")) for r in diagnostic_rows])[1],
        "diag_g2_cw0_min": _finite_stats([r.get("diag_g2_cw0", r.get("g2_cw0")) for r in diagnostic_rows])[0],
        "diag_g2_cw0_median": _finite_stats([r.get("diag_g2_cw0", r.get("g2_cw0")) for r in diagnostic_rows])[1],
        "diag_g2_cw0_raw_min": _finite_stats([r.get("diag_g2_cw0_raw", r.get("g2_cw0_raw")) for r in diagnostic_rows])[0],
        "diag_g2_cw0_raw_median": _finite_stats([r.get("diag_g2_cw0_raw", r.get("g2_cw0_raw")) for r in diagnostic_rows])[1],
        "flux_max": _max_finite([r.get("collected_flux_pulsed_s") for r in diagnostic_rows]),
        "best_diagnostic_row": _best_diagnostic_row(diagnostic_rows),
    }


def _max_finite(values):
    vals = []
    for value in values:
        try:
            if value is not None and np.isfinite(float(value)):
                vals.append(float(value))
        except (TypeError, ValueError):
            pass
    return max(vals) if vals else float("nan")


def _best_diagnostic_row(rows):
    """Lowest g2_pulsed wins; ties are broken by the HIGHEST collected flux,
    not the lowest. g2_op does not depend on the emission.NA/R_back/L_um
    collection levers (device.py's edge out-coupling factor only multiplies
    into brightness, never into g2), so council review round 4's lever axes
    (item 1) now routinely produce many rows tied on g2_pulsed that differ
    only in collected flux; the "favourable diagnostic corner" this feeds
    (write_markdown's brightness decomposition) is the best-g2,
    most-measurable operating point, so the highest-flux tie is preferred."""
    valid = []
    for row in rows:
        try:
            if np.isfinite(float(row.get("g2_pulsed", float("nan")))):
                valid.append(row)
        except (TypeError, ValueError):
            pass

    def _key(row):
        g2 = float(row["g2_pulsed"])
        try:
            flux = float(row.get("collected_flux_pulsed_s", float("nan")))
        except (TypeError, ValueError):
            flux = float("nan")
        return (g2, -flux if np.isfinite(flux) else float("inf"))

    return min(valid, key=_key) if valid else None


def _self_test_passed(report: dict) -> bool:
    """verify/verify_rt_edge_papers.py (out of this file's scope) is
    concurrently being renamed by another worker: its self-test flag is
    moving from `hallucination_tests_passed` to `all_checks_passed`. Read
    whichever key the installed module actually returns, preferring the
    new name."""
    if "all_checks_passed" in report:
        return bool(report["all_checks_passed"])
    return bool(report.get("hallucination_tests_passed", False))


def _f_or_none(value):
    try:
        return float(value) if value is not None and np.isfinite(float(value)) else None
    except (TypeError, ValueError):
        return None


def _fmt_or_na(value, spec: str = ".4g") -> str:
    """Markdown-table-cell formatter: "n/a" for None/unparseable, "nan" for
    a genuine non-finite float, else `spec`-formatted."""
    if value is None:
        return "n/a"
    try:
        fval = float(value)
    except (TypeError, ValueError):
        return "n/a"
    return "nan" if not np.isfinite(fval) else format(fval, spec)


def _card_gamma300_pass_max(rows: list) -> dict:
    """Per-card ceiling on the (currently unmeasured, [E]-class) 300 K
    linewidth gamma300 (council review round 4, item 2): the largest
    gamma300_meV at which ANY eligible row (any lever/delta_xx/irf
    combination) still clears the headline pulsed g2(0)<0.5 gate --
    restates the verdict as conditional on this one unmeasured input
    rather than as an absolute platform ceiling. Returns {card_id: {...}},
    "gamma300_pass_max_meV": nan when the card has no headline-passing row."""
    by_card: dict = {}
    for r in rows:
        by_card.setdefault(r["card_id"], []).append(r)
    result = {}
    for card_id, card_rows in by_card.items():
        passing = [r for r in card_rows if r["headline_pass"]]
        if not passing:
            result[card_id] = {"gamma300_pass_max_meV": float("nan")}
            continue
        gmax = max(float(r["gamma300_meV"]) for r in passing)
        at_max = [r for r in passing if float(r["gamma300_meV"]) == gmax]
        best = min(at_max, key=lambda r: float(r["g2_pulsed"]))
        result[card_id] = {
            "gamma300_pass_max_meV": gmax,
            "g2_pulsed": _f_or_none(best.get("g2_pulsed")),
            "delta_xx_meV": _f_or_none(best.get("delta_xx_meV")),
            "irf_ps": _f_or_none(best.get("irf_ps")),
            "emission_NA": _f_or_none(best.get("emission_NA")),
            "emission_R_back": _f_or_none(best.get("emission_R_back")),
            "emission_L_um": _f_or_none(best.get("emission_L_um")),
            "assumptions": str(best.get("assumptions", "")),
        }
    return result


def _no_pass_reason(card_rows: list) -> str:
    """Council review round 6, item 2: WHY a card has zero headline-passing
    gamma300 samples. `"flux"` when every sampled row is below the
    collected-flux eligibility floor (headline_pass can never be True for
    an ineligible row, no matter how good g2 would be -- e.g. the primary
    card's favourable corner reaches g2=0.26 at 1 meV but only 313 photons/s,
    below the 1 kHz floor); `"g2"` when at least one row IS eligible but
    g2(0) >= 0.5 at every sampled gamma300 (a genuine physics shortfall, not
    a flux-eligibility one)."""
    flux_blocks = not any(r.get("eligible") for r in card_rows)
    g2_blocks = not any(_f_or_none(r.get("g2_pulsed")) is not None
                        and _f_or_none(r.get("g2_pulsed")) < G2_THRESHOLD
                        for r in card_rows)
    reasons = [name for name, blocked in (("flux", flux_blocks), ("g2", g2_blocks)) if blocked]
    return ",".join(reasons) or "g2"


def _gamma300_threshold_bracket(rows: list, gamma300_pass_max_by_card: dict) -> dict:
    """Council review round 5, item 2: gamma300_pass_max alone only says the
    true (continuous, unsampled) threshold is AT LEAST this value -- it says
    nothing about how much further it might extend before failing. This
    finds, per card, the SMALLEST sampled gamma300 (from the same `rows`
    _card_gamma300_pass_max was given) strictly above gamma300_pass_max,
    i.e. the first sample already known to fail; the true threshold lies in
    (gamma300_pass_max, that value]. Returns {card_id: {"lo": ..., "hi":
    ...}}: hi is None when every sampled value >= gamma300_pass_max passes
    (the threshold is >= the largest sample); lo is None when no sampled
    value passes at all (the threshold is < the smallest sample) -- in that
    case a `"reason"` key (council review round 6, item 2) is also set,
    `"flux"` or `"g2"` per _no_pass_reason, so _format_threshold can print
    `none(flux)`/`none(g2)` instead of a bracket like `<6` that would
    otherwise misleadingly imply the card passes below that value."""
    by_card_gammas: dict = {}
    by_card_rows: dict = {}
    for r in rows:
        by_card_gammas.setdefault(r["card_id"], set()).add(float(r["gamma300_meV"]))
        by_card_rows.setdefault(r["card_id"], []).append(r)
    result = {}
    for card_id, gammas in by_card_gammas.items():
        gmax = gamma300_pass_max_by_card.get(card_id, {}).get("gamma300_pass_max_meV", float("nan"))
        sorted_g = sorted(gammas)
        if not np.isfinite(gmax):
            result[card_id] = {"lo": None, "hi": (sorted_g[0] if sorted_g else None),
                               "reason": _no_pass_reason(by_card_rows.get(card_id, []))}
            continue
        higher = sorted(g for g in sorted_g if g > gmax)
        result[card_id] = {"lo": gmax, "hi": (higher[0] if higher else None)}
    return result


def _pooled_gamma300_threshold(gamma300_pass_max_by_card: dict, threshold_by_card: dict) -> dict:
    """VERDICT-line pooled bracket (council review round 5, item 2): the
    bracket of whichever card attains the pooled gamma300_pass_max (the max
    over all cards); ties broken by CARDS order. {"lo": None, "hi": None}
    when no card has any headline-passing sample."""
    finite = [(cid, info["gamma300_pass_max_meV"])
             for cid, info in gamma300_pass_max_by_card.items()
             if np.isfinite(info.get("gamma300_pass_max_meV", float("nan")))]
    if not finite:
        return {"lo": None, "hi": None}
    best_cid = max(finite, key=lambda kv: kv[1])[0]
    return threshold_by_card.get(best_cid, {"lo": None, "hi": None})


def _format_threshold(bracket: dict) -> str:
    """Human-readable `gamma300_threshold` text for the VERDICT line/
    verdict.md: "<lo>-<hi>" when both bounds are known, ">=<lo>" when every
    sampled value at/above gamma300_pass_max passes, "n/a" when there is no
    card data at all. When nothing passes, council review round 6 item 2:
    a numeric "<hi" bracket would misleadingly imply the card passes
    somewhere below `hi`, when in fact NO sample passed at all -- print the
    REASON class instead, `none(flux)` (every sample flux-ineligible) or
    `none(g2)` (eligible but g2(0) >= 0.5 everywhere), per the bracket's own
    `reason` (set by _gamma300_threshold_bracket); a bare bracket with no
    `reason` key falls back to the old "<hi" text."""
    lo, hi = (bracket or {}).get("lo"), (bracket or {}).get("hi")
    if lo is None and hi is None:
        return "n/a"
    if lo is None:
        reason = (bracket or {}).get("reason")
        return f"none({reason})" if reason else f"<{hi:g}"
    if hi is None:
        return f">={lo:g}"
    return f"{lo:g}-{hi:g}"


def _loading_term(mu) -> float:
    """The Poisson cap-2 loading multiplier actually used in
    brightness_per_pulse (fsim_core.loading.loading_probs: P0=e^-mu,
    P1+P2 = 1-P0 EXACTLY -- brightness_per_pulse's own P1+P2 factor, loading.py
    lines 49-54/69-72) -- a plain algebraic identity on the already-computed
    mu_resolved scalar, not new physics; council review round 4 item 3's
    complaint was that the old factor table printed bare `mu` as if IT were
    the multiplier, when the actual multiplier is this quantity."""
    try:
        mu_f = float(mu)
    except (TypeError, ValueError):
        return float("nan")
    return 1.0 - np.exp(-mu_f) if np.isfinite(mu_f) else float("nan")


def _brightness_factor_check(row: dict) -> dict:
    """Reconstructs collected_flux_pulsed_s as the plain product of
    device.py's own already-computed scalars -- loading (1-e^-mu), t_X, S
    (confinement retention), edge_eta_total (edge out-coupling, which
    itself already contains beta/eta_facet/eta_NA exactly once each --
    T_facet and single-pass propagation are folded into eta_facet's own
    ray series, not separate factors of eta_total; peer-review pkg2 facet
    fix, 2026-09-07, item 3), times rep_rate_hz -- with NO duty term (duty
    is already folded into rep_rate_hz's own derivation, device.py's
    rep_rate_hz = duty_eff /
    tau_pulse_ns) and no separate beta/facet/propagation/NA terms (already
    inside edge_eta_total): council review round 4 item 3's exact
    complaint was that the old table's `duty` row and standalone `beta` row
    were NOT independent multipliers on top of rep_rate/eta_total, so
    multiplying every printed row together double-counted / mis-stated the
    chain. This IS that same chain, reconstructed field-by-field, so its
    product times rep_rate reproduces the reported flux to float precision
    (checked here at a generous 1% tolerance)."""
    mu, t_x = row.get("mu_pulsed"), row.get("t_x_pulsed")
    S, eta_total = row.get("S_retention_pulsed"), row.get("edge_eta_total")
    rep_rate, reported_flux = row.get("rep_rate_hz"), row.get("collected_flux_pulsed_s")
    loading = _loading_term(mu)
    fields = [loading, t_x, S, eta_total, rep_rate]
    try:
        fields_f = [float(v) for v in fields]
    except (TypeError, ValueError):
        fields_f = [float("nan")] * len(fields)
    if not all(np.isfinite(v) for v in fields_f):
        return {"ok": False, "loading": loading, "product_flux": float("nan"),
                "reported_flux": _f_or_none(reported_flux) or float("nan"),
                "rel_diff": float("nan")}
    loading_f, t_x_f, S_f, eta_f, rep_f = fields_f
    product_flux = loading_f * t_x_f * S_f * eta_f * rep_f
    reported_flux_f = _f_or_none(reported_flux)
    if reported_flux_f is None or reported_flux_f == 0:
        rel_diff = float("nan")
    else:
        rel_diff = abs(product_flux - reported_flux_f) / reported_flux_f
    ok = bool(np.isfinite(rel_diff) and rel_diff < 0.01)
    return {"ok": ok, "loading": loading_f, "t_x": t_x_f, "S": S_f, "eta_total": eta_f,
            "rep_rate": rep_f, "product_flux": product_flux,
            "reported_flux": reported_flux_f if reported_flux_f is not None else float("nan"),
            "rel_diff": rel_diff}


def _front_facet_split(row: dict) -> float:
    """Peer-review pkg2 facet fix (2026-09-07,
    .workers/specs/pr-pkg2-facet-fix.md item 3), updated again from council
    review round 5, item 6: fsim_core/waveguide.py's edge_emission() now
    folds single-pass propagation entirely into the facet ray-series
    (facet_escape_fraction), so `eta_prop` is no longer a factor in
    eta_total at all -- dividing by it here would silently reintroduce a
    propagation term that was never applied on this path. The STRUCTURAL
    equation is now eta_total = beta * eta_facet * eta_NA (no eta_prop);
    this backs out eta_facet as the one factor missing, deliberately never
    dividing by edge_T_facet or edge_eta_prop separately."""
    try:
        beta = float(row.get("edge_beta"))
        eta_na = float(row.get("edge_eta_NA"))
        eta_total = float(row.get("edge_eta_total"))
    except (TypeError, ValueError):
        return float("nan")
    denom = beta * eta_na
    if not (np.isfinite(denom) and denom != 0 and np.isfinite(eta_total)):
        return float("nan")
    return eta_total / denom


def _facet_factor_forward_check(row: dict) -> dict:
    """Peer-review pkg2 facet fix (2026-09-07,
    .workers/specs/pr-pkg2-facet-fix.md item 3): the OLD version of this
    check scraped a facet-factor formula string out of the live
    fsim_core/waveguide.py source and eval()-ed it in a bare {T, R_back}
    namespace -- after the pkg2 checkpoint, edge_emission() assigns
    `facet_factor = facet_escape_fraction(T, R_back_eff, alpha_cm, L_um,
    dot_position)`, a name reference the eval() namespace could never
    resolve, so every row's forward value silently failed to evaluate and
    `ok` was False everywhere.  This calls
    `waveguide.facet_escape_fraction(T, R_back, alpha_cm, L_um)` directly
    instead -- the SAME pure function edge_emission() itself now calls
    (single source of truth), so this is a wiring/regression check, not a
    duplicated physics derivation. `emission_R_back=None` resolves the same
    way edge_emission() resolves it for an uncoated stack, `R_back = 1 - T`
    (item 4's uncoated-Fresnel resolution reduces to this whenever no
    coating override is in play, which this sweep never sets). `alpha_cm`
    is read from the row's own `emission_alpha_cm` column (peer-review pkg2
    fix3, 2026-09-07, item 3 -- it is not a swept lever, but IS recorded
    per-row from `design.emission.alpha_cm`); a stale CSV predating that
    column falls back to the fixed `fsim_core.device.EmissionBlock.alpha_cm`
    class default, noted in `convention` as `"ray-series-midpoint (alpha
    from card)"` so the fallback is visible in the printed verdict rather
    than silently indistinguishable from a genuinely per-row value. Compares
    the forward value against the back-solved `eta_total/(beta*eta_NA)` from
    `_front_facet_split`. Returns {"back_solved", "forward", "convention",
    "ok"}; `convention` is `None` (not the string) whenever `forward`
    could not be evaluated at all -- either a malformed row or a
    facet_escape_fraction exception -- so the verdict text can print
    "eta_facet could not be evaluated" instead of a nonsensical "forward =
    nan via ray-series-midpoint"."""
    back_solved = _front_facet_split(row)
    try:
        T = float(row.get("edge_T_facet"))
        L_um = float(row.get("emission_L_um"))
        r_back_raw = row.get("emission_R_back")
        R_back = None if r_back_raw in (None, "", "None") else float(r_back_raw)
    except (TypeError, ValueError):
        return {"back_solved": back_solved, "forward": float("nan"),
               "convention": None, "ok": False}
    if R_back is None:
        R_back = 1.0 - T  # uncoated Fresnel resolution reduces to this here (item 4)
    alpha_raw = row.get("emission_alpha_cm")
    if alpha_raw in (None, "", "None"):
        # stale CSV predating the emission_alpha_cm column: fall back to
        # the fixed class default and say so in convention.
        alpha_cm = 5.0  # fsim_core.device.EmissionBlock.alpha_cm default
        convention = "ray-series-midpoint (alpha from card)"
    else:
        try:
            alpha_cm = float(alpha_raw)
        except (TypeError, ValueError):
            return {"back_solved": back_solved, "forward": float("nan"),
                   "convention": None, "ok": False}
        convention = "ray-series-midpoint"
    try:
        forward = float(waveguide.facet_escape_fraction(T, R_back, alpha_cm, L_um))
    except Exception:
        return {"back_solved": back_solved, "forward": float("nan"),
               "convention": None, "ok": False}
    ok = bool(np.isfinite(back_solved) and np.isfinite(forward)
             and abs(forward - back_solved) < 1e-6)
    return {"back_solved": back_solved, "forward": forward,
           "convention": convention, "ok": ok}


def compute_verdict(rows: list, stats: dict, grid_complete: bool,
                    evidence_report: dict, hallucination_report: dict,
                    gamma300_refine_rows: list | None = None,
                    anchor_report: dict | None = None) -> dict:
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
    hallucination_ok = _self_test_passed(hallucination_report)
    headline_rows = [r for r in rows if r["headline_pass"]]
    passing_temperatures = [int(T) for T, info in stats.get("per_T", {}).items()
                            if info.get("n_headline", 0) > 0]
    T_pass_min = min(passing_temperatures) if passing_temperatures else None
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
    # Council review round 4, item 2: `conditional` no longer means "PASS
    # relies on some assumption" (true of essentially every corner, since
    # every row is [A]/[E]-tagged -- an uninformative signal). It now marks
    # the specific, actionable state where the physics already clears the
    # headline gate but the run cannot declare PASS only because the
    # evidence gate (verify_rt_edge_papers.py's two-source-per-claim rule)
    # is not yet satisfied -- i.e. a result that is conditional on evidence
    # completion, not on the science itself.
    conditional = bool(headline_rows) and not evidence_complete

    # Council review round 5, item 2: prefer the finer gamma300_refine_rows
    # (GAMMA300_REFINE_MEV, 8 samples) when the caller supplies it; fall
    # back to the main `rows` (RANGE_BOUNDS's 2 endpoints only) so this
    # function stays usable exactly as before when no refinement was run
    # (e.g. --quick, or the existing section-1 synthetic-fixture tests).
    gamma_source_rows = gamma300_refine_rows if gamma300_refine_rows else rows
    gamma300_pass_max_by_card = _card_gamma300_pass_max(gamma_source_rows)
    gamma300_pass_max = _max_finite(
        [v["gamma300_pass_max_meV"] for v in gamma300_pass_max_by_card.values()])
    gamma300_threshold_by_card = _gamma300_threshold_bracket(
        gamma_source_rows, gamma300_pass_max_by_card)
    gamma300_threshold = _pooled_gamma300_threshold(
        gamma300_pass_max_by_card, gamma300_threshold_by_card)
    gamma300_threshold_by_T = {}
    for T in stats.get("per_T", {}):
        t_rows = [r for r in gamma_source_rows if str(int(round(float(r.get("T_hs_K", 300))))) == T]
        if t_rows:
            t_max = _card_gamma300_pass_max(t_rows)
            gamma300_threshold_by_T[T] = _pooled_gamma300_threshold(
                t_max, _gamma300_threshold_bracket(t_rows, t_max))

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

    # Council review round 5, item 4 (naming/semantics -- one word per
    # quantity): `coverage` now means headline-pass fraction OVER ELIGIBLE
    # ROWS ONLY (it used to duplicate headline_coverage_n/total's all-rows
    # fraction under a different-looking float, which is what made it
    # ambiguous); `eligible_fraction` (NEW) is eligible/total. `flux_margin`
    # (NEW, preferred) is flux_max/floor (>1 means the floor is CLEARED);
    # `flux_shortfall` (its old floor/flux_max inverse framing) is kept one
    # release, deprecated.
    eligible_fraction = (stats["n_eligible"] / stats["n_total"]) if stats["n_total"] else 0.0
    coverage_over_eligible = ((stats["n_headline"] / stats["n_eligible"])
                              if stats["n_eligible"] else float("nan"))
    flux_shortfall = (FLUX_FLOOR_PULSED_S / stats["flux_max"]
                      if np.isfinite(stats["flux_max"]) and stats["flux_max"] > 0
                      else float("nan"))
    flux_margin = (stats["flux_max"] / FLUX_FLOOR_PULSED_S
                  if np.isfinite(stats["flux_max"]) else float("nan"))

    return {
        "pass": passed, "fail_reasons": fail_reasons,
        "g2_min": g2_min, "g2_median": g2_median, "median_pass": median_pass,
        "coverage": coverage_over_eligible, "eligible_fraction": eligible_fraction,
        "headline_coverage_n": stats["n_headline"], "headline_coverage_total": stats["n_total"],
        "headline_coverage_pulsed_n": stats["n_headline_dedup"],
        "headline_coverage_pulsed_total": stats["n_scheduled_dedup"],
        "headline_dedup_mismatch_groups": stats.get("headline_dedup_mismatch_groups", 0),
        "eligible_dedup_n": stats.get("n_eligible_dedup", 0),
        "eligible_dedup_total": stats["n_scheduled_dedup"],
        "eligible_dedup": stats.get("eligible_dedup", 0.0),
        "eligible_dedup_mismatch_groups": stats.get("eligible_dedup_mismatch_groups", 0),
        "cw0_coverage_n": stats["n_cw0_pass"], "cw0_coverage_total": stats["n_eligible"],
        "cw_raw_coverage_n": stats["n_cw_raw_pass"], "cw_raw_coverage_total": stats["n_eligible"],
        "eligible_n": stats["n_eligible"], "eligible_total": stats["n_total"],
        "flux_floor_excluded": stats["n_flux_floor_excluded"],
        "diag_g2_min": stats["diag_g2_pulsed_min"],
        "diag_g2_median": stats["diag_g2_pulsed_median"],
        "flux_max": stats["flux_max"],
        "flux_margin": flux_margin,
        "flux_shortfall": flux_shortfall,  # DEPRECATED: use flux_margin (1/flux_shortfall)
        "evidence_complete": evidence_complete, "hallucination_ok": hallucination_ok,
        "conditional": conditional, "assumptions_used": assumptions_used,
        "headline_rows": headline_rows, "headline_by_card": headline_by_card,
        "fallback_only": fallback_only, "note": note,
        "gamma300_pass_max": gamma300_pass_max,
        "gamma300_pass_max_by_card": gamma300_pass_max_by_card,
        "gamma300_threshold_by_card": gamma300_threshold_by_card,
        "gamma300_threshold": gamma300_threshold,
        "anchor_6_5mev_by_card": anchor_report or {},
        "T_pass_min": T_pass_min,
        "headline_by_T": stats.get("per_T", {}),
        "headline_by_T_pulsed": stats.get("per_T_pulsed", {}),
        "gamma300_threshold_by_T": gamma300_threshold_by_T,
    }


def verdict_line(verdict: dict) -> str:
    # Peer review finding 6: appended, not substituted -- headline_coverage
    # (existing key, unchanged) stays for backward compatibility with
    # verdict dicts saved by the pre-fix writer (e.g. the stale saved
    # manifest.json/verdict.md this file's grid has not rerun yet), which
    # do not carry the new keys; omit this segment entirely for those so
    # the reconstructed line still equals the stale saved text exactly.
    pulsed_segment = ""
    by_T_pulsed_segment = ""
    has_pulsed_dedup = "headline_coverage_pulsed_n" in verdict
    if has_pulsed_dedup:
        pulsed_segment = (
            f"headline_coverage_pulsed={verdict['headline_coverage_pulsed_n']}/"
            f"{verdict['headline_coverage_pulsed_total']} "
            f"headline_dedup_mismatch_groups={verdict.get('headline_dedup_mismatch_groups', 0)} "
            f"eligible_dedup={verdict.get('eligible_dedup_n', 0)}/"
            f"{verdict.get('eligible_dedup_total', 0)} "
            f"eligible_dedup_mismatch_groups={verdict.get('eligible_dedup_mismatch_groups', 0)} "
            f"rows_scheduled={verdict['headline_coverage_n']}/{verdict['headline_coverage_total']} ")
        by_T_pulsed_segment = " headline_by_T_pulsed=" + ",".join(
            f"{T}:{info.get('n_headline_dedup', 0)}/{info.get('n_scheduled_dedup', 0)}"
            for T, info in verdict.get('headline_by_T_pulsed', {}).items())
    return (f"VERDICT: {'PASS' if verdict['pass'] else 'FAIL'} "
            f"g2_min={verdict['g2_min']:.4g} g2_median_eligible={verdict['g2_median']:.4g} "
            f"diag_g2_min={verdict['diag_g2_min']:.4g} "
            f"diag_g2_median_diagnostic={verdict['diag_g2_median']:.4g} "
            f"flux_max={verdict['flux_max']:.4g} "
            f"flux_margin={verdict['flux_margin']:.4g} "
            f"flux_shortfall_deprecated={verdict['flux_shortfall']:.4g} "
            f"median_pass={'true' if verdict['median_pass'] else 'false'} "
            f"coverage_over_eligible={verdict['coverage']:.4g} "
            f"eligible_fraction={verdict['eligible_fraction']:.4g} "
            f"eligible={verdict['eligible_n']}/{verdict['eligible_total']} "
            f"flux_floor_excluded={verdict['flux_floor_excluded']} "
            f"evidence={'complete' if verdict['evidence_complete'] else 'incomplete'} "
            f"conditional={'true' if verdict['conditional'] else 'false'} "
            f"headline_coverage={verdict['headline_coverage_n']}/{verdict['headline_coverage_total']} "
            + pulsed_segment +
            f"cw_raw_coverage={verdict['cw_raw_coverage_n']}/{verdict['cw_raw_coverage_total']} "
            f"gamma300_pass_max={verdict['gamma300_pass_max']:.4g} "
            f"gamma300_threshold={_format_threshold(verdict['gamma300_threshold'])} "
            f"T_pass_min={verdict.get('T_pass_min') if verdict.get('T_pass_min') is not None else 'none'} "
            f"headline_by_T=" + ",".join(
                f"{T}:{info.get('n_headline', 0)}/{info.get('n_total', 0)}"
                for T, info in verdict.get('headline_by_T', {}).items())
            + by_T_pulsed_segment)


def card_line(card_id: str, card_stats: dict) -> str:
    """Council review round 5, item 1: the fixed 'diagnostic (below flux
    floor, not measurable)' phrase used to print unconditionally right
    after g2_pulsed_min/median -- which are ALREADY the eligible-row
    statistics (compute_stats' eligible_pulsed_vals) -- wrongly implying
    those numbers themselves were unmeasurable diagnostics whenever a card
    had ANY eligible rows at all (e.g. the gainp card, eligible=16/96,
    flux >= floor). Print that phrase ONLY when the card has zero eligible
    rows (so g2_pulsed_min/median are genuinely nan and only the diag_*
    figures below carry information); otherwise print the eligible count
    with its own min/median."""
    if card_stats["n_eligible"] == 0:
        eligibility_phrase = "diagnostic (below flux floor, not measurable)"
    else:
        eligibility_phrase = (
            f"eligible rows: {card_stats['n_eligible']}/{card_stats['n_rows']} "
            f"g2_pulsed_min={card_stats['g2_pulsed_min']:.4g} "
            f"g2_pulsed_median={card_stats['g2_pulsed_median']:.4g}")
    return (f"CARD: {card_id} role={card_stats['card_class']} "
            f"{eligibility_phrase} "
            f"diag_g2_min={card_stats['diag_g2_pulsed_min']:.4g} "
            f"diag_g2_median={card_stats['diag_g2_pulsed_median']:.4g} "
            f"diag_g2_cw0_min={card_stats['diag_g2_cw0_min']:.4g} "
            f"diag_g2_cw0_median={card_stats['diag_g2_cw0_median']:.4g} "
            f"diag_g2_cw0_raw_min={card_stats['diag_g2_cw0_raw_min']:.4g} "
            f"diag_g2_cw0_raw_median={card_stats['diag_g2_cw0_raw_median']:.4g} "
            f"g2_cw_raw_min={card_stats['g2_cw0_raw_min']:.4g} "
            f"g2_cw_raw_median={card_stats['g2_cw0_raw_median']:.4g} "
            f"eligible={card_stats['n_eligible']}/{card_stats['n_rows']} "
            f"flux_floor_excluded={card_stats['n_flux_floor_excluded']} "
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
                   hallucination_report: dict, path: Path,
                   lever_info: dict | None = None) -> None:
    lever_info = lever_info or {}
    lines = []
    lines.append("# RT edge-emitter acceptance sweep verdict")
    lines.append("")
    lines.append(f"Generated {datetime.now(timezone.utc).isoformat()}; "
                 f"contract: docs/rt_edge_contract.md.")
    lines.append("")
    lines.append(f"```\n{verdict_line(verdict)}\n```")
    for card_id, cs in stats["per_card"].items():
        lines.append(f"```\n{card_line(card_id, cs)}\n```")
    # pkg5 fix, item 1: the irf_ps axis is deduplicated with a group all()
    # reduction (a corner passes only if it passes at every sampled irf_ps),
    # not a first-row pick -- surface any group where the irf_ps axis
    # actually disagrees, rather than resolving it silently.
    n_headline_mismatch = stats.get("headline_dedup_mismatch_groups", 0)
    n_eligible_mismatch = stats.get("eligible_dedup_mismatch_groups", 0)
    if n_headline_mismatch or n_eligible_mismatch:
        warning = (f"**WARNING:** the irf_ps axis disagrees within "
                   f"{n_headline_mismatch} dedup group(s) on headline_pass "
                   f"(headline_dedup_mismatch_groups={n_headline_mismatch}) and "
                   f"{n_eligible_mismatch} dedup group(s) on eligible_row "
                   f"(eligible_dedup_mismatch_groups={n_eligible_mismatch}); a group "
                   f"counts as passing only if every sampled irf_ps value agrees.")
        lines.append("")
        lines.append(warning)
        print(warning.replace("**", ""), file=sys.stderr)
    lines.append("")
    lines.append("## Literature ceiling")
    lines.append(
        "In the reviewed literature set of this repository (six papers, "
        "../_goal/paper_digests.md) and the anchors ledger, no electrically "
        "driven III-V single-dot g2(0) at 300 K is reported; the best "
        "electrical result in that set is Reischle et al., Optics Express 16, "
        "12771 (2008) at 80 K (QD C): g2(0) = "
        f"{REISCHLE_RAW_G2:.2f} raw (still IRF-broadened), 0.03 after "
        "background correction with the ~500 ps detector IRF also deconvolved. "
        "The best reported 300 K single-dot values "
        "are optically pumped: g2(0) ~ 0.5-0.57 (Laferriere et al. 2023, "
        "InAsP/InP nanowire dot, g2(0) = 0.57 at 300 K). This sweep's pooled "
        f"pulsed g2_min={verdict['g2_min']:.4g}, "
        f"g2_median(eligible)={verdict['g2_median']:.4g} is reported against "
        "that electrical-vs-optical literature picture, not as a claim of an "
        "existing electrical 300 K result to exceed.")
    lines.append("")
    lines.append(
        "**Convention-matched comparison (council review round 5, item 5; "
        "peer review finding 9).** The paragraph above juxtaposed this "
        "sweep's IRF-FREE intrinsic g2_min against Reischle's RAW dip "
        f"({REISCHLE_RAW_G2:.2f} [V], still IRF-broadened; ledger anchor "
        f"`{_RAW_G2_ANCHOR_ID}`) -- not a convention-matched comparison. The "
        "correct convention-matched anchor is Reischle's IRF-DECONVOLVED-"
        f"but-background-included value, g2_b(0) = {REISCHLE_DECONV_G2:.2f} "
        f"[V] +/- {REISCHLE_DECONV_G2_ERR:.2f} (QD C, 80 K; ledger anchor "
        f"`{_DECONV_G2_ANCHOR_ID}`, ../_goal/paper_digests.md "
        "line ~36), since this sweep's g2_op is likewise background-included "
        "(via drive.b_res/rho) and never IRF-convolved for the pulsed "
        "metric. Pulsed peak-area g2(0) and CW zero-delay g2(0) remain "
        "distinct observables even after matching background and IRF.")
    # Peer review finding 9, item 1: the temperature that actually produces
    # the pooled g2_min (today 230 K, not 300 K) must be named, and both
    # ratios reported -- taken from stats["per_T"], never hardcoded.
    g2_min = verdict.get("g2_min", float("nan"))
    per_T_stats = stats.get("per_T", {})
    argmin_T = None
    for T, info in per_T_stats.items():
        val = info.get("g2_min")
        if (val is not None and np.isfinite(val) and np.isfinite(g2_min)
                and abs(val - g2_min) < 1e-9):
            argmin_T = T
            break
    if np.isfinite(g2_min) and g2_min > 0:
        ratio = g2_min / REISCHLE_DECONV_G2
        if argmin_T == "300":
            corner_label = "best 300 K corner"
        elif argmin_T is not None:
            # pkg5 fix, item 6: no inner parens around "at N K" -- callers
            # below already wrap the whole label in "(g2_min=... )", so a
            # second, nested pair here rendered as the double parenthesis
            # "best corner (at 230 K) (g2_min=...)".
            corner_label = f"best corner at {argmin_T} K"
        else:
            corner_label = "best corner (heatsink temperature not resolved in this row set)"
        if ratio > 1.0:
            lines.append(
                f"On matching conventions, the 80 K Reischle device is about "
                f"{ratio:.2f}x better (lower g2(0)) than this sweep's "
                f"{corner_label} (g2_min={g2_min:.4g} vs "
                f"{REISCHLE_DECONV_G2:.2f} +/- {REISCHLE_DECONV_G2_ERR:.2f}).")
        else:
            lines.append(
                f"On matching conventions, this sweep's {corner_label} "
                f"(g2_min={g2_min:.4g}) is already at or below the 80 K "
                f"Reischle deconvolved value ({REISCHLE_DECONV_G2:.2f} +/- "
                f"{REISCHLE_DECONV_G2_ERR:.2f}).")
        g2_min_300 = per_T_stats.get("300", {}).get("g2_min")
        if (g2_min_300 is not None and np.isfinite(g2_min_300) and g2_min_300 > 0
                and argmin_T != "300"):
            ratio_300 = g2_min_300 / REISCHLE_DECONV_G2
            # pkg5 fix, item 6: when argmin_T is unresolved, "an unresolved
            # temperature" already says there is no number to give a unit
            # to -- appending " K" after it produced "an unresolved
            # temperature K".
            argmin_T_text = f"{argmin_T} K" if argmin_T is not None else "an unresolved temperature"
            lines.append(
                f"On the 300 K line specifically (not this sweep's pooled "
                f"best, which is at {argmin_T_text}"
                f"): pulsed g2_min={g2_min_300:.4g}, about {ratio_300:.2f}x "
                f"{'worse than' if ratio_300 > 1.0 else 'at or below'} the "
                f"80 K Reischle deconvolved value.")
    else:
        lines.append(
            "No finite pooled g2_min is available in this run to compare "
            "against Reischle's deconvolved value.")
    lines.append("")
    lines.append("## Conditional on the 300 K linewidth")
    gamma300_samples_text = ", ".join(f"{g:g}" for g in GAMMA300_REFINE_MEV)
    lines.append(
        "The 300 K single-dot linewidth gamma300 is not a platform ceiling -- it is an "
        "unmeasured [E]-class quantity (see Assumptions below) whose actual value decides "
        "whether the headline gate passes. `gamma300_pass_max` (council review round 5, "
        f"item 2) is the HIGHEST SAMPLED linewidth, among this determination's own "
        f"{len(GAMMA300_REFINE_MEV)}-point grid ({gamma300_samples_text} meV -- finer than "
        "the sweep.csv grid's own 2-point dot.gamma300 endpoints, and never itself written "
        "to sweep.csv), at which any eligible row (any lever/delta_xx combination) still "
        "clears pulsed intrinsic g2(0) < 0.5; `gamma300_threshold` brackets the true "
        "(continuous, unsampled) threshold between that value and the next sampled value "
        "above it that already fails -- the true threshold lies somewhere inside the "
        "bracket, at the delta_xx shown for that card's gamma300_pass_max row.")
    lines.append("")
    lines.append(
        "**When nothing passes (council review round 6, item 2).** If NO sampled gamma300 "
        "clears the headline gate for a card, `gamma300_threshold` no longer prints a "
        "numeric bracket like `<6` -- that would misleadingly imply the card passes "
        "somewhere below 6 meV, when in fact it never passes anywhere in this grid. It "
        "instead prints the REASON class: `none(flux)` if every sampled row is below the "
        "collected-flux eligibility floor (headline_pass can never be True there, no "
        "matter how favourable g2 would be), or `none(g2)` if at least one sampled row IS "
        "flux-eligible but pulsed intrinsic g2(0) >= 0.5 at every sampled gamma300 (a "
        "genuine physics shortfall, not an eligibility one); `none(flux,g2)` reports both "
        "blocking conditions when neither a flux-eligible row nor a diagnostic g2 < 0.5 row exists.")
    lines.append("")
    lines.append(f"`gamma300_pass_max` (pooled, both cards): {verdict['gamma300_pass_max']:.4g} meV; "
                 f"`gamma300_threshold` (pooled): "
                 f"{_format_threshold(verdict['gamma300_threshold'])} meV")
    lines.append("")
    lines.append("| card | gamma300_pass_max (meV) | gamma300_threshold bracket (meV) | "
                 "delta_xx | lever values | assumptions |")
    lines.append("|---|---:|---:|---:|---|---|")
    for card_id, info in verdict["gamma300_pass_max_by_card"].items():
        gmax = info.get("gamma300_pass_max_meV", float("nan"))
        bracket = verdict["gamma300_threshold_by_card"].get(card_id, {})
        if not np.isfinite(gmax):
            lines.append(f"| {card_id} | n/a (no headline-passing sample in this grid) | "
                         f"{_format_threshold(bracket)} | | | |")
            continue
        def _g(key):
            v = info.get(key)
            return f"{v:g}" if v is not None else "n/a"
        lever_text = f"NA={_g('emission_NA')}, R_back={_g('emission_R_back')}, L_um={_g('emission_L_um')}"
        lines.append(f"| {card_id} | {gmax:.4g} | {_format_threshold(bracket)} | "
                     f"{_g('delta_xx_meV')} meV | {lever_text} | {info.get('assumptions', '')} |")
    lines.append("")
    lines.append(
        "Restated (council review round 5, item 2): for each card, `gamma300_pass_max` is "
        "the highest SAMPLED linewidth at which the card passes; the true threshold lies "
        "between the bracket's two values (at the delta_xx shown above for that row).")
    lines.append("")
    lines.append(f"### Verified {CHATZARAKIS_ANCHOR_GAMMA300_MEV:g} meV anchor "
                 "(Chatzarakis et al., Phys. Rev. Applied 20, 034011, 2023) -- an InAs/GaAs "
                 "(211)B single-dot linewidth measurement, used here as a distinct-material "
                 "class proxy [E] for the InP-dot family, not a direct InP/GaAsP or "
                 "InP/GaInP measurement")
    lines.append(
        "Fresh evaluation, per card, at the card's OWN default delta_xx (never the sweep's "
        "endpoint grid) -- both cards' own provenance.ranges[\"dot.gamma300\"].note already "
        "document this by hand; the numbers below are computed programmatically.")
    lines.append("")
    lines.append("| card | delta_xx (meV, card default) | g2_pulsed | collected flux (photons/s) | "
                 "eligible | passes (<0.5) |")
    lines.append("|---|---:|---:|---:|---|---|")
    anchor_by_card = verdict.get("anchor_6_5mev_by_card", {}) or {}
    for card_id, info in anchor_by_card.items():
        lines.append(
            f"| {card_id} | {_fmt_or_na(info.get('delta_xx_meV'))} | "
            f"{_fmt_or_na(info.get('g2_pulsed'))} | {_fmt_or_na(info.get('flux_s'))} | "
            f"{info.get('eligible')} | {'yes' if info.get('passes') else 'no'} |")
    if anchor_by_card:
        anchor_sentences = []
        for card_id, info in anchor_by_card.items():
            verdict_word = "PASSES" if info.get("passes") else "does not pass"
            anchor_sentences.append(
                f"at gamma300 = {CHATZARAKIS_ANCHOR_GAMMA300_MEV:g} meV (this card's own "
                f"delta_xx = {_fmt_or_na(info.get('delta_xx_meV'))} meV), card `{card_id}` "
                f"{verdict_word} the headline gate (g2_pulsed = "
                f"{_fmt_or_na(info.get('g2_pulsed'))}, collected flux = "
                f"{_fmt_or_na(info.get('flux_s'))} photons/s)")
        lines.append("")
        lines.append("; ".join(anchor_sentences) + ".")
    lines.append("")
    lines.append("## Grid")
    lines.append(f"Grid complete: {grid_complete}"
                 + ("" if grid_complete else " (--quick: endpoints-only, cannot grant PASS)"))
    for path_key, (lo, hi, unit) in RANGE_BOUNDS.items():
        vals = ", ".join(f"{v:g}" for v in grid[path_key])
        lines.append(f"- `{path_key}` in [{lo}, {hi}] {unit}: sampled at {vals}")
    if lever_info:
        lines.append("")
        lines.append("Collection-lever axes (per-card, council review round 4 item 1):")
        for card_id, info in lever_info.items():
            for lever_path, values in info.get("grid", {}).items():
                vals = ", ".join(f"{v:g}" if v is not None else "None" for v in values)
                lines.append(f"- `{card_id}` `{lever_path}`: sampled at {vals}")
            lines.append(f"- `{card_id}` lever combinations evaluated: {info.get('n_combos')}")
    lines.append("")
    lines.append("## Per-temperature acceptance")
    lines.append("| T_hs (K) | eligible | headline passes | deduplicated (IRF axis collapsed) | "
                 "pulsed g2 min | flux max (photons/s) | gamma300 threshold |")
    lines.append("|---:|---:|---:|---:|---:|---:|---|")
    headline_by_T_pulsed = verdict.get("headline_by_T_pulsed", {})
    for T, info in verdict.get("headline_by_T", {}).items():
        threshold = verdict.get("gamma300_threshold_by_T", {}).get(T, {"lo": None, "hi": None})
        dedup_info = headline_by_T_pulsed.get(T, {})
        dedup_cell = (f"{dedup_info['n_headline_dedup']}/{dedup_info['n_scheduled_dedup']}"
                     if dedup_info else "n/a")
        lines.append(f"| {T} | {info['n_eligible']}/{info['n_total']} | {info['n_headline']}/{info['n_total']} | "
                     f"{dedup_cell} | "
                     f"{_fmt_or_na(info['g2_min'])} | {_fmt_or_na(info['flux_max'])} | {_format_threshold(threshold)} |")
    lines.append("")
    lines.append("## What cooling buys")
    for T, info in verdict.get("headline_by_T", {}).items():
        candidates = [r for r in rows if int(round(float(r.get("T_hs_K", -1)))) == int(T)]
        best = _best_diagnostic_row(candidates)
        if best:
            lines.append(f"At T_hs={T} K, the favourable corner has retention S={_fmt_or_na(best.get('S_retention_pulsed'))}, "
                         f"linewidth Gamma(T)={_fmt_or_na(best.get('Gamma_pulsed_meV'))} meV, and window background "
                         f"b_e={_fmt_or_na(best.get('b_e_window_pulsed'))}.")
    lines.append("")
    lines.append("## Coverage")
    lines.append(f"- **eligible fraction** (eligible/total; council review round 5, item 4 -- "
                 f"one word per quantity): {stats['n_eligible']}/{stats['n_total']} "
                 f"= {verdict['eligible_fraction']:.3f}")
    lines.append(f"- pulsed collected-flux eligibility floor [A]: "
                 f"{FLUX_FLOOR_PULSED_S:.0f} photons/s; rows excluded by this floor: "
                 f"{stats['n_flux_floor_excluded']}")
    lines.append(f"- diagnostic pooled g2 (below flux floor, not measurable): "
                 f"pulsed {stats['diag_g2_pulsed_min']:.4g} / {stats['diag_g2_pulsed_median']:.4g}; "
                 f"g2_cw0 {stats['diag_g2_cw0_min']:.4g} / {stats['diag_g2_cw0_median']:.4g}; "
                 f"g2_cw0_raw {stats['diag_g2_cw0_raw_min']:.4g} / {stats['diag_g2_cw0_raw_median']:.4g}")
    lines.append(f"- maximum collected pulsed flux: {stats['flux_max']:.4g} photons/s; "
                 f"**flux_margin** (flux_max/floor; >1 means the floor is CLEARED): "
                 f"{verdict['flux_margin']:.4g} (flux_shortfall, DEPRECATED, its old "
                 f"floor/flux_max inverse framing: {verdict['flux_shortfall']:.4g})")
    lines.append(f"- **rows_scheduled** (pulsed intrinsic g2(0) < 0.5; invalid/ineligible "
                 f"rows count as nonpassing; raw grid count, every irf_ps sample counted "
                 f"separately -- `headline_coverage` in the VERDICT line, kept for backward "
                 f"compatibility): {stats['n_headline']}/{stats['n_total']} = "
                 f"{stats['headline_coverage']:.3f}")
    lines.append(f"- **headline_coverage_pulsed** (peer review finding 6, pkg5 fix item 2: "
                 f"same numerator rule as rows_scheduled, but the irf_ps axis is "
                 f"deduplicated first. The pulsed g2 sub-result is IRF-independent; "
                 f"`headline_pass` also requires CW eligibility, which is IRF-convolved, "
                 f"so the deduplicated count is computed with `all()` over the IRF axis "
                 f"and any disagreement is reported as `headline_dedup_mismatch_groups` "
                 f"(currently {stats.get('headline_dedup_mismatch_groups', 0)}) rather than "
                 f"resolved by picking whichever irf_ps row happened to be seen first): "
                 f"{stats['n_headline_dedup']}/{stats['n_scheduled_dedup']} = "
                 f"{stats['headline_coverage_pulsed']:.3f}. Coverage is the fraction of a "
                 f"chosen endpoint grid that passes, not a fabrication-yield probability or "
                 f"a confidence level.")
    if stats["n_eligible"]:
        lines.append(f"- **headline coverage over ELIGIBLE rows only** (same numerator, "
                     f"denominator restricted to eligible rows -- this is "
                     f"`coverage_over_eligible` in the VERDICT line, council review round 6 "
                     f"item 4, one name per quantity): "
                     f"{stats['n_headline']}/{stats['n_eligible']} = {verdict['coverage']:.3f}")
    else:
        lines.append("- **headline coverage over ELIGIBLE rows only**: n/a (0 eligible rows)")
    lines.append(f"- pooled pulsed g2 median OVER ELIGIBLE ROWS ONLY (`g2_median_eligible` in "
                 f"the VERDICT line; diagnostic only, it never appears in fail_reasons): "
                 f"{stats['g2_pulsed_median']:.4g}")
    lines.append(f"- pooled pulsed g2 median over ALL VALID/DIAGNOSTIC rows (eligible rows "
                 f"plus rows excluded ONLY by the flux floor; `diag_g2_median_diagnostic` in "
                 f"the VERDICT line): {stats['diag_g2_pulsed_median']:.4g}")
    lines.append(f"- secondary coverage, CW intrinsic g2_cw0 < 0.5 (eligible rows, diagnostic "
                 f"only, does not gate PASS): {stats['n_cw0_pass']}/{stats['n_eligible']} "
                 f"= {stats['cw0_coverage']:.3f}")
    lines.append(f"- secondary coverage, CW IRF-convolved g2_cw0_raw < 0.5 (eligible rows, "
                 f"diagnostic only, does not gate PASS): {stats['n_cw_raw_pass']}/{stats['n_eligible']} "
                 f"= {stats['cw_raw_coverage']:.3f}")
    lines.append("")
    lines.append("## Per-card statistics")
    lines.append("| card | role | diagnostic g2_pulsed min/median | diagnostic g2_cw0 min/median | "
                 "diagnostic g2_cw0_raw min/median | eligible | flux-floor excluded | headline rows |")
    lines.append("|---|---|---|---|---|---|---|---|")
    for card_id, cs in stats["per_card"].items():
        lines.append(
            f"| {card_id} | {cs['card_class']} | "
            f"{cs['diag_g2_pulsed_min']:.4g} / {cs['diag_g2_pulsed_median']:.4g} | "
            f"{cs['diag_g2_cw0_min']:.4g} / {cs['diag_g2_cw0_median']:.4g} | "
            f"{cs['diag_g2_cw0_raw_min']:.4g} / {cs['diag_g2_cw0_raw_median']:.4g} | "
            f"{cs['n_eligible']}/{cs['n_rows']} | {cs['n_flux_floor_excluded']} | "
                 f"{cs['n_favorable']} |")
    if stats.get("best_diagnostic_row") is not None:
        best = stats.get("best_diagnostic_row")
        lines.append("")
        if stats["n_eligible"] == 0:
            lines.append("## Why no row is eligible")
            lines.append(f"Every row is below the {FLUX_FLOOR_PULSED_S:.0f} photons/s collected-flux floor; "
                         f"the grid maximum is only {stats['flux_max']:.4g} photons/s "
                         f"(floor/maximum = {verdict['flux_shortfall']:.4g}).")
        else:
            lines.append("## Best diagnostic-g2 row and brightness decomposition")
        # Council review round 4, item 3: print the ACTUAL multiplicative
        # chain (loading term, t_X, S, eta_total, rep rate) with a
        # self-check that their product reproduces the reported flux --
        # not a decorative list that double-counts beta/facet/propagation/
        # NA (already folded into eta_total) or duty (already folded into
        # rep_rate_hz's own derivation).
        check = _brightness_factor_check(best)
        chain = [("loading = 1 - e^-mu (mu={:.4g})".format(
                     float(best["mu_pulsed"]) if best.get("mu_pulsed") is not None
                     and np.isfinite(best.get("mu_pulsed")) else float("nan")),
                 check.get("loading")),
                ("t_X (spectral transmission)", check.get("t_x")),
                ("S (confinement retention)", check.get("S")),
                ("eta_total (edge out-coupling: waveguide coupling x front/back facet "
                 "split x facet transmission x propagation x NA, all in one factor)",
                 check.get("eta_total")),
                ("rep rate (Hz)", check.get("rep_rate"))]
        facet_check = _facet_factor_forward_check(best)
        front_split = facet_check["back_solved"]
        component_factors = [("beta (waveguide coupling / spontaneous-emission factor)",
                              best.get("edge_beta")),
                             ("eta_facet (front/back-facet split, transmission, and "
                              "single-pass propagation all fused into one ray-series factor; "
                              "BACK-SOLVED as the one factor missing from "
                              "eta_total/(beta*eta_NA) -- device.py folds it into eta_total "
                              "but never exposes it on its own; see the independent check "
                              "below)", front_split),
                             ("T_facet (raw Fresnel transmission, diagnostic only -- NOT "
                              "an independent multiplicative step beyond eta_facet above; "
                              "see the independent check below)",
                              best.get("edge_T_facet")),
                             ("NA (numerical aperture)", best.get("edge_eta_NA"))]
        # Council review round 6, item 3: the dominant-limiter comparison
        # used to run over eta_total's sub-factors only (component_factors:
        # beta, the combined facet factor, T_facet, propagation, NA), which
        # can never surface a smaller limiter sitting OUTSIDE eta_total in
        # the chain (loading, t_X, S) -- e.g. S=0.0046 at the favourable
        # corner, well below beta=0.029. Comparing across the WHOLE chain
        # (loading, t_X, S, plus eta_total's own sub-factors -- never
        # eta_total or rep_rate themselves, which are aggregates/a rate,
        # not standalone limiting fractions) makes the printed dominant
        # limiter the smallest factor anywhere in the reported flux product.
        whole_chain_candidates = [chain[0], chain[1], chain[2]] + component_factors
        finite_components = [(n, float(v)) for n, v in whole_chain_candidates
                             if v is not None and np.isfinite(v) and v > 0]
        dominant = (min(finite_components, key=lambda x: x[1])[0]
                   if finite_components else "the collection chain")
        lines.append(f"The dominant brightness limiter at the favourable diagnostic corner -- "
                     f"the smallest factor across the WHOLE chain (loading, t_X, S, and "
                     f"eta_total's own sub-factors; council review round 6, item 3) -- is "
                     f"{dominant}; the multiplicative chain that reproduces the reported "
                     f"collected pulsed flux is:")
        lines.append("")
        lines.append("| factor | value |")
        lines.append("|---|---:|")
        for name, value in chain:
            lines.append(f"| {name} | {value:.6g} |" if value is not None and np.isfinite(value)
                         else f"| {name} | nan |")
        lines.append(f"| **product x rep rate** | {check['product_flux']:.6g} photons/s |"
                     if np.isfinite(check.get("product_flux", float("nan")))
                     else "| **product x rep rate** | nan |")
        lines.append(f"| reported collected_flux_pulsed_s | {check['reported_flux']:.6g} photons/s |"
                     if np.isfinite(check.get("reported_flux", float("nan")))
                     else "| reported collected_flux_pulsed_s | nan |")
        lines.append(
            f"| self-check: relative difference | {check['rel_diff']:.4%} "
            f"({'PASS, within 1%' if check['ok'] else 'FAIL, exceeds 1%'}) |"
            if np.isfinite(check.get("rel_diff", float("nan")))
            else "| self-check: relative difference | nan (non-finite inputs) |")
        lines.append("")
        lines.append(f"`beta`/`eta_facet`/`T_facet`/`NA` are shown below for diagnosis only -- "
                     f"`beta`, `eta_facet` and `NA` are already folded into `eta_total` above "
                     f"exactly once each (their product reproduces eta_total, by construction "
                     f"of eta_facet -- NOT independent evidence, see below) and must NOT also "
                     f"be multiplied into the flux self-check. Single-pass propagation is NOT "
                     f"listed here as a separate multiplicative row (peer-review pkg2 facet "
                     f"fix, 2026-09-07): it is reported below as an informational line only, "
                     f"already folded inside `eta_facet`'s ray series.")
        lines.append("")
        lines.append("| component (already inside eta_total) | value |")
        lines.append("|---|---:|")
        for name, value in component_factors:
            lines.append(f"| {name} | {value:.6g} |" if value is not None and np.isfinite(value)
                         else f"| {name} | nan |")
        lines.append("")
        edge_eta_prop = best.get("edge_eta_prop")
        lines.append(
            f"single-pass propagation (already inside the facet factor): "
            f"{edge_eta_prop:.6g}" if edge_eta_prop is not None and np.isfinite(edge_eta_prop)
            else "single-pass propagation (already inside the facet factor): nan")
        lines.append("")
        lines.append(
            "**Independent check on eta_facet** (council review round 5, item 6, fixing a "
            "tautology; peer-review pkg2 facet fix, 2026-09-07): the value above is "
            "BACK-SOLVED -- the one factor missing once eta_total, beta and NA are all "
            "already known (eta_total = beta * eta_facet * eta_NA) -- so a self-check "
            "comparing it against that same back-solving would only ever reproduce its own "
            "inputs. The genuinely independent check instead recomputes eta_facet FORWARD "
            "from this row's own facet transmission (T_facet), back-facet reflectivity "
            "(R_back), and ridge loss (alpha_cm), by a direct call to "
            "fsim_core.waveguide.facet_escape_fraction -- the SAME pure function "
            "edge_emission() itself calls (single source of truth, never a duplicated or "
            "scraped formula), so this is a wiring/regression check, not a physics "
            "re-derivation: "
            + (f"forward = {facet_check['forward']:.6g} via {facet_check['convention']}, "
               if facet_check.get("convention") else "eta_facet could not be evaluated, ")
            + f"back-solved = {front_split:.6g} -- "
            f"{'MATCH' if facet_check['ok'] else 'MISMATCH'} (agreement confirms the "
            "back-solved value really is the eta_facet factor waveguide.py computes, not "
            "some other quantity folded into eta_total).")
        lines.append("")
        lines.append("## Diagnostic g2 landscape")
        lines.append("| corner | pulsed g2 min | pooled diagnostic median | lever values | assumptions |")
        lines.append("|---|---:|---:|---|---|")
        def _grid_text(key):
            value = best.get(key)
            return f"{value:g}" if value is not None and np.isfinite(value) else "n/a"
        lever_text = (f"NA={_grid_text('emission_NA')}, R_back={_grid_text('emission_R_back')}, "
                     f"L_um={_grid_text('emission_L_um')}")
        lines.append(f"| {best.get('card_id', 'unknown')} ({_grid_text('delta_xx_meV')} meV, "
                     f"gamma300={_grid_text('gamma300_meV')} meV, "
                     f"irf={_grid_text('irf_ps')} ps) | "
                         f"{best['g2_pulsed']:.6g} | {stats['diag_g2_pulsed_median']:.6g} | "
                         f"{lever_text} | {best.get('assumptions', '')} |")
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
    lines.append("## Assumptions in plain words")
    lines.append(
        "**1. The collection window assumption makes flux look independent of the "
        "300 K linewidth.** Every row's spectral filter is set automatically to exactly "
        "match the exciton line's own width at the operating temperature "
        "(device.py's `filter.auto_w` convention; `filter.auto_w_scale` can widen or "
        "narrow it, 1.0 = as-is here). A plain mathematical consequence of matching the "
        "window to the line's own width is that exactly HALF of the X-line's photons fall "
        "inside the window in every single row of this sweep (t_X = 0.5, column "
        "`t_x_pulsed`) -- no matter how wide or narrow the real (currently unmeasured) "
        "300 K linewidth gamma300 turns out to be. That is why the reported collected "
        "photon flux does not move with gamma300 in this sweep: it is not that the real "
        "device would be insensitive to linewidth, it is that this filter-window "
        "convention cancels that sensitivity out by construction.")
    lines.append("")
    # Council review round 6, item 1: the citation below used to be a
    # hardcoded string ("Reischle et al., Appl. Phys. Lett. 92, 233113
    # (2008)") that had drifted from the actual source of the rho ~ 0.88
    # residual-background anchor (Optics Express 16, 12771 (2008), DOI
    # 10.1364/OE.16.012771). It is now generated FROM the ledger anchor
    # (verify/data/rt_edge_anchors.yaml, id `reischle08-b-res-80k`) so it
    # cannot drift again.
    _bres_anchor = _rt_anchors[BG_RES_ANCHOR_ID]
    _rho = float(_bres_anchor["value"])
    lines.append(
        "**2. The background-light assumption is borrowed from a different, colder "
        "device.** Every row carries a constant background term (`drive.b_res`) that is "
        f"not measured on this platform: it is transferred from {_bres_anchor['source']} "
        f"(DOI {_bres_anchor['doi']} [{_bres_anchor['tag']}], ledger anchor "
        f"`{BG_RES_ANCHOR_ID}`), whose 80 K electrically driven single-photon "
        f"source had about {_rho * 100:.0f}% real signal and {(1 - _rho) * 100:.0f}% "
        f"background light (signal fraction rho ~ {_rho:.2f}). This sweep assumes the "
        f"same {(1 - _rho) * 100:.0f}% background fraction still applies "
        "at 300 K, on a different material system (InP/GaAsP or InP/GaInP edge "
        "emitters) than the one actually measured. No 300 K electrical background "
        "measurement exists for either card's platform.")
    lines.append("")
    _plain_best = stats.get("best_diagnostic_row")
    if _plain_best is not None and np.isfinite(_plain_best.get("g2_cw0", float("nan"))) \
            and np.isfinite(_plain_best.get("g2_cw0_raw", float("nan"))):
        # pkg5-fix2, item 5: "ABOVE the 0.5 threshold" was asserted
        # unconditionally, but this best-diagnostic corner's raw CW g2(0)
        # is not guaranteed to actually clear 0.5 (e.g. the live-like
        # fixture's 0.321 does not) -- render the threshold verdict from
        # the value itself instead of hardcoding it.
        _cw_raw = _plain_best['g2_cw0_raw']
        _cw_threshold_phrase = ("ABOVE the 0.5 threshold" if _cw_raw >= G2_THRESHOLD
                                else "still below the 0.5 threshold")
        lines.append(
            "**3. Why pulsed drive, not continuous-wave (CW) drive, is required.** At the "
            f"best diagnostic operating point in this sweep, the intrinsic CW g2(0) is "
            f"{_plain_best['g2_cw0']:.3g} (that alone would already satisfy the g2 < 0.5 "
            f"single-photon criterion), but once a realistic single-photon detector's "
            f"finite timing resolution (instrument response function, IRF) is folded in, "
            f"the measured raw CW g2(0) rises to {_cw_raw:.3g} -- {_cw_threshold_phrase}. "
            f"In plain terms: the antibunching dip this device produces "
            f"under continuous drive is narrower in time than a real detector can resolve, "
            f"so single-photon emission cannot be demonstrated by a CW measurement alone "
            f"at the IRF values sampled here (50-200 ps); a faster detector, a different "
            f"gate or different physical rates could change this. Pulsed (gated) operation "
            f"sidesteps the detector's timing resolution and is therefore required at "
            f"these IRF values.")
    else:
        lines.append(
            "**3. Why pulsed drive, not continuous-wave (CW) drive, is required.** No "
            "diagnostic row with both a finite intrinsic and IRF-convolved CW g2(0) is "
            "available in this run to quote a concrete pair of numbers, but the "
            "underlying reason pulsed drive is used is the same in every corner: the "
            "antibunching dip this device produces under continuous drive is narrower in "
            "time than a real single-photon detector's instrument response, so "
            "single-photon emission cannot be demonstrated by a CW measurement alone on "
            "this platform at the IRF values sampled here (50-200 ps); a faster detector, "
            "a different gate or different physical rates could change this.")
    lines.append("")
    lines.append("## Evidence gate")
    lines.append(f"- evidence_complete: {evidence_report.get('evidence_complete')}")
    lines.append(f"- self_test_passed (self-test): {_self_test_passed(hallucination_report)}")
    missing = evidence_report.get("missing_evidence", [])
    if missing:
        lines.append("- missing/incomplete evidence claims:")
        for m in missing:
            lines.append(f"  - `{m['claim']}`: {m['reason']}")
    lines.append("")
    lines.append("## Fail reasons" if not verdict["pass"] else "## Pass basis")
    if verdict["pass"]:
        # pkg5-fix2, item 1: the PASS gate itself is row-based (compute_verdict
        # requires at least one raw row with headline_pass, never the
        # deduplicated count), so quoting only the deduplicated count here
        # could read as self-contradictory when every dedup group disagrees
        # across the IRF axis (n_headline_dedup == 0 while PASS is still
        # true). Quote the dedup count only when it is actually >= 1;
        # otherwise state plainly that PASS rests on raw rows whose IRF-axis
        # partner disagrees.
        if stats['n_headline_dedup'] >= 1:
            lines.append(f"At least one eligible corner has pulsed intrinsic g2(0) < 0.5 -- "
                         f"the contract's headline metric "
                         f"({stats['n_headline_dedup']} such corner(s), deduplicated over the "
                         f"IRF axis; {stats['n_headline']} raw grid row(s), across both cards), "
                         f"grid complete, evidence complete, hallucination "
                         f"self-test passed. g2_cw0 and g2_cw0_raw are reported above as "
                         f"secondary diagnostics (see Coverage) and do not gate this PASS.")
        else:
            lines.append(
                f"PASS rests on {stats['n_headline']} raw grid row(s) whose IRF-axis partner "
                f"disagrees (headline_dedup_mismatch_groups = "
                f"{stats.get('headline_dedup_mismatch_groups', 0)}); no corner passes at "
                f"every sampled IRF. Grid complete, evidence complete, hallucination "
                f"self-test passed. g2_cw0 and g2_cw0_raw are reported above as secondary "
                f"diagnostics (see Coverage) and do not gate this PASS.")
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
                   path: Path, lever_info: dict | None = None) -> None:
    manifest = {
        "schema_version": 1,
        "generated_utc": datetime.now(timezone.utc).isoformat(),
        "quick": quick, "grid_complete": grid_complete, "grid": grid,
        "range_bounds": RANGE_BOUNDS,
        # Per-card resolved collection-lever grid (council review round 4,
        # item 1) -- {card_id: {"grid": {path: [values]}, "n_combos": int}};
        # unlike range_bounds/grid above these are per-card, not shared.
        "lever_info": lever_info or {},
        "pulse_assumption": PULSE_ASSUMPTION,
        "cards": [{"id": c["id"], "class": c["class"], "path": str(c["path"]),
                  "sha256": _sha256_file(c["path"])} for c in CARDS],
        "runtime_versions": {
            "python": platform.python_version(), "numpy": np.__version__,
            "scipy": scipy.__version__, "yaml": getattr(yaml, "__version__", "unknown"),
        },
        "policy": {
            "g2_threshold": G2_THRESHOLD,
            "pulsed_collected_flux_floor_s": FLUX_FLOOR_PULSED_S,
            "pulsed_collected_flux_floor_basis": "[A] below ~1 kHz a g2 measurement is not practical within hours at single-photon-detector counting rates",
            "pass_requires": ["grid_complete", "evidence_complete",
                              "self_test_passed",
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
        # Field renamed from hallucination_tests_passed to self_test_passed
        # (spec item 5); verify/verify_rt_edge_papers.py (out of this file's
        # scope) is concurrently renaming its OWN key from
        # hallucination_tests_passed to all_checks_passed -- _self_test_passed
        # reads whichever key that module actually returns.
        "hallucination_self_test": {"self_test_passed": _self_test_passed(hallucination_report)},
        "runtime_seconds": runtime_s,
    }
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(manifest, indent=2, sort_keys=False, default=str), encoding="utf-8")


# ------------------------------------------------------------------- main

def run_sweep(quick: bool) -> tuple:
    grid = build_grid(quick)
    pulsed_cache: dict = {}
    rows = []
    lever_info = {}
    for card in CARDS:
        card_rows, lever_grid, n_combos = sweep_card(card, grid, pulsed_cache, quick)
        rows.extend(card_rows)
        lever_info[card["id"]] = {"grid": lever_grid, "n_combos": n_combos}
    stats = compute_stats(rows)
    # Council review round 5, item 2: the finer gamma300_pass_max/threshold
    # determination (GAMMA300_REFINE_MEV) is skipped for --quick (kept fast,
    # matching --quick's existing "smaller, explicitly incomplete"
    # convention for every other axis -- an incomplete grid can never PASS
    # regardless, so the extra resolution buys nothing there). The verified
    # 6.5 meV anchor check is two eval_pulsed_point calls total (negligible
    # cost) and is always computed.
    gamma300_refine_rows = []
    if not quick:
        passing_T = [float(T) for T, info in stats.get("per_T", {}).items()
                     if info.get("n_headline", 0) > 0]
        # The eight-point linewidth refinement is deliberately restricted to
        # the lowest passing TEC set point and the retained 300 K headline.
        refine_T = sorted(set(([min(passing_T)] if passing_T else []) + [300.0]))
        for card in CARDS:
            design0 = DeviceDesign.load(card["path"])
            # Refinement is a threshold diagnostic, not another full
            # collection-lever sweep: use each card's declared operating
            # corner so the eight samples remain affordable at two T values.
            combos = [{"emission.NA": design0.emission.NA,
                       "emission.R_back": design0.emission.R_back,
                       "emission.L_um": design0.emission.L_um}]
            gamma300_refine_rows.extend(
                refine_gamma300(card, [design0.dot.delta_xx], combos, pulsed_cache, refine_T))
    anchor_report = {card["id"]: anchor_check(card, pulsed_cache) for card in CARDS}
    return rows, stats, grid, lever_info, gamma300_refine_rows, anchor_report


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--quick", action="store_true",
                        help="small endpoints-only grid; explicitly incomplete, cannot grant PASS")
    parser.add_argument("--out-dir", default=str(ROOT / "out" / "rt_edge"))
    args = parser.parse_args(argv)
    out_dir = Path(args.out_dir)

    t0 = time.time()
    rows, stats, grid, lever_info, gamma300_refine_rows, anchor_report = run_sweep(args.quick)

    # Fresh paper-check invocation (spec: "Invoke paper-check run_checks()
    # freshly"); this file never writes evidence.json (out of scope --
    # verify/verify_rt_edge_papers.py owns that file).
    evidence_report = rt_papers.run_checks(self_test=False)
    hallucination_report = rt_papers.run_checks(self_test=True)

    grid_complete = not args.quick
    verdict = compute_verdict(rows, stats, grid_complete, evidence_report, hallucination_report,
                              gamma300_refine_rows=(gamma300_refine_rows or None),
                              anchor_report=anchor_report)
    runtime_s = time.time() - t0

    try:
        write_csv(rows, out_dir / "sweep.csv")
        write_png(rows, stats, out_dir / "envelope.png")
        write_markdown(rows, stats, verdict, grid, grid_complete,
                       evidence_report, hallucination_report, out_dir / "verdict.md",
                       lever_info=lever_info)
        write_manifest(rows, stats, verdict, grid, grid_complete, evidence_report,
                       hallucination_report, runtime_s, args.quick, out_dir / "manifest.json",
                       lever_info=lever_info)
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
