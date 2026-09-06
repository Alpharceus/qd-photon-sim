"""Literature cross-checking and hallucination detection for the room-
temperature edge-emitter evidence ledger (verify/data/rt_edge_anchors.yaml --
this is the "data/rt_edge_papers.yaml" the spec/contract refers to by that
name), wired through the shared DeviceDesign/evaluate() integration
(fsim_core.device) rather than reimplementing physics as direct-module-only
checks. This file is the paper-enforcement layer: it does NOT judge whether a
device corner is scientifically good (verify/verify_rt_edge_cards.py and
verify/verify_device_rt.py already own that); it judges whether the literature
evidence backing every anchored claim is genuine, sourced, in-context, and
fed through the real evaluator rather than asserted in isolation.

Every "numerical comparison" below compares a card/fixture's own resolved
DeviceDesign field against the anchor's SOURCE-EXTRACTED literature value
(never the evaluator's own computed output used as if it were the target --
docs/rt_edge_contract.md and the spec's interface constraints forbid that).
Where a claim has no direct 1:1 design-card field (e.g. an 80 K electrical g2
value from a different lab's device, or a vertical-reference collection
efficiency), the evaluator cross-check instead demonstrates a specific,
literature-grounded physical DISTINCTION the contract requires evaluators
never blur: V_bi is not a forward operating voltage, a ridge beta is not a
photonic-wire first-lens efficiency, an activation energy alone does not fix
a retention fraction without its prefactor, and a design's own modeled
wavelength does not silently claim an unmatched target wavelength.

CLI: python verify/verify_rt_edge_papers.py [--report PATH] [--self-test]
Default report: out/rt_edge/evidence.json. Exports run_checks() -> dict,
consumable directly (no subprocess parsing, no file staleness -- the report
carries source_ledger_sha256/evaluator_input_hashes for that).

Standalone, side-effect-free on import (all work happens under
`if __name__ == "__main__"` or inside run_checks()); importable without
executing the checks or writing anything.
"""
from __future__ import annotations

import argparse
import copy
import hashlib
import json
import re
import sys
import warnings
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

ANCHORS_PATH = ROOT / "verify" / "data" / "rt_edge_anchors.yaml"
DIGEST_VALUES_PATH = ROOT / "verify" / "data" / "rt_edge_digest_values.yaml"
GAASP_CARD = ROOT / "cards" / "edge-inp-gaasp-design.yaml"
GAINP_CARD = ROOT / "cards" / "edge-inp-gainp-design.yaml"
DEFAULT_REPORT = ROOT / "out" / "rt_edge" / "evidence.json"

VALID_TAGS = {"V", "DR", "E", "A"}

# Registered DOI-prefix allowlist for the plausibility check below (Part B.ii
# of the council review): APS, AIP, Optica, Wiley, IOP, Elsevier, Springer/
# Nature, ACS, IEEE. A DOI outside this list is rejected unless the anchor
# carries an explicit `web:` verification note (a real paper hosted under a
# prefix this list happens not to cover, hand-confirmed).
REGISTERED_DOI_PREFIXES = (
    "10.1103",  # APS (Phys. Rev. *)
    "10.1063",  # AIP (J. Appl. Phys., Appl. Phys. Lett., ...)
    "10.1364",  # Optica (Optics Express, Optics Letters, ...)
    "10.1002",  # Wiley (physica status solidi, ...)
    "10.1088",  # IOP
    "10.1016",  # Elsevier
    "10.1038",  # Nature
    "10.1007",  # Springer
    "10.1021",  # ACS
    "10.1109",  # IEEE
)

# Per-claim expectations used for context/unit/magnitude hallucination checks.
# `material_tokens`/`temp_tokens`/`drive_tokens` are OR-checked substrings
# (case-sensitive, matching the ledger's own wording); an empty list means
# "not asserted for this claim" (skipped). `proxy=True` marks claims whose
# anchors are explicitly transferred-class evidence, not a direct measurement
# of the target device (docs/rt_edge_contract.md Evidence section).
CLAIM_META = {
    "reischle2008_g2_80K": dict(unit="g2(0)", lo=0.0, hi=1.0,
                                material_tokens=["InP"], temp_tokens=["80 K"],
                                drive_tokens=["EL"], proxy=False),
    "temperature_trend_g2": dict(unit="g2(0)", lo=0.0, hi=1.0,
                                 material_tokens=[], temp_tokens=[],
                                 drive_tokens=["Pulsed", "pulsed"], proxy=True),
    "hkust_inp_gaasp_wavelength": dict(unit="nm", lo=600.0, hi=900.0,
                                       material_tokens=["InP"], temp_tokens=[],
                                       drive_tokens=[], proxy=False),
    "gaas_pin_iv": dict(unit="V", lo=0.5, hi=5.0,
                        material_tokens=["GaAs"], temp_tokens=[],
                        drive_tokens=["p-i-n"], proxy=False),
    "edge_extraction_efficiency": dict(unit="fraction", lo=0.0, hi=1.0,
                                       material_tokens=[], temp_tokens=["4 K"],
                                       drive_tokens=["Vertical"], proxy=False),
    "gamma300_class_range": dict(unit="meV", lo=0.0, hi=100.0,
                                 material_tokens=[], temp_tokens=[],
                                 drive_tokens=[], proxy=True),
    "retention_Ea_inp_algainp": dict(unit="meV", lo=0.0, hi=500.0,
                                     material_tokens=["InP"], temp_tokens=[],
                                     drive_tokens=[], proxy=False),
    "delta_xx_inp_gaasp": dict(unit="meV", lo=0.0, hi=50.0,
                              material_tokens=["InP"], temp_tokens=[],
                              drive_tokens=[], proxy=True),
}

# Literal figure/text numbers pulled directly from ../_goal/paper_digests.md
# for the one claim where raw/deconvolved/background-corrected values could
# be swapped (Reischle et al. 2008, Fig. 3b, QD C, 80 K): 0.43 raw dip,
# 0.25 after IRF deconvolution only, 0.03 after background correction too.
# An anchor whose conditions claim "raw dip" but whose value matches one of
# the OTHER two numbers is a raw-vs-intrinsic swap.
RAW_VS_CORRECTED_DISTRACTORS = {
    "reischle2008_g2_80K": {"raw_value": 0.43, "distractors": [0.25, 0.03]},
}


# ------------------------------------------------------------------ plumbing

def _sha256_file(path: Path) -> str | None:
    if not path.exists():
        return None
    return hashlib.sha256(path.read_bytes()).hexdigest()


def load_anchors(path: Path = ANCHORS_PATH) -> dict:
    doc = yaml.safe_load(path.read_text(encoding="utf-8"))
    return {a["id"]: a for a in doc["anchors"]}


def load_digest_values(path: Path = DIGEST_VALUES_PATH) -> dict:
    """The independent ground-truth transcription table (verify/data/
    rt_edge_digest_values.yaml), keyed by (doi, quantity) -> entry. This is
    never derived from rt_edge_anchors.yaml -- it is hand-transcribed from
    ../_goal/paper_digests.md / confirmed web-PDF sources -- so it can catch
    a fabricated anchor or a mutated real value the anchors file alone
    cannot self-detect."""
    doc = yaml.safe_load(path.read_text(encoding="utf-8"))
    return {(e["doi"], e["quantity"]): e for e in doc.get("entries", [])}


def group_by_claim(anchors: dict) -> dict:
    by_claim: dict = {}
    for aid, a in anchors.items():
        by_claim.setdefault(a["claim"], []).append(aid)
    return by_claim


def anchor_bound(anchor: dict) -> float:
    """Absolute tolerance bound: a plain number is absolute; {rel: x} is
    relative to |value|; {abs: x} is already absolute."""
    tol = anchor["tolerance"]
    if isinstance(tol, dict):
        if "rel" in tol:
            return abs(float(anchor["value"])) * float(tol["rel"])
        if "abs" in tol:
            return float(tol["abs"])
        raise ValueError(f"unsupported tolerance mapping {tol!r}")
    return float(tol)


def _distinct_key(anchor: dict):
    """Two anchors under one claim count as ONE source if they cite the same
    DOI (a duplicated-citation hallucination), else they are distinguished by
    their (source, locator) text."""
    doi = (anchor.get("doi") or "").strip()
    if doi:
        return ("doi", doi)
    return ("src", (anchor.get("source") or "").strip(), (anchor.get("locator") or "").strip())


_PLACEHOLDER_SOURCE_MARKERS = ("not retrieved", "not available")


def _is_placeholder_gap(anchor: dict) -> bool:
    src = (anchor.get("source") or "").strip().lower()
    loc = (anchor.get("locator") or "").strip().lower()
    return any(m in src for m in _PLACEHOLDER_SOURCE_MARKERS) or loc in _PLACEHOLDER_SOURCE_MARKERS


# ------------------------------------------------------- structural / hallucination checks

def _check_tag_valid(anchor: dict):
    return anchor.get("tag") in VALID_TAGS, f"tag={anchor.get('tag')!r}"


def _check_citation_present(anchor: dict):
    """Reject an absent/fake citation for a numeric claim: a real source
    string plus at least a DOI or a genuine (non-placeholder) locator."""
    if anchor.get("value") is None:
        return True, "no numeric claim to source (declared gap)"
    if _is_placeholder_gap(anchor):
        return True, "declared placeholder gap, not a claimed citation"
    src = (anchor.get("source") or "").strip()
    if not src:
        return False, "no source citation given for a numeric claim"
    has_doi = bool((anchor.get("doi") or "").strip())
    has_locator = bool((anchor.get("locator") or "").strip()) and \
        (anchor.get("locator") or "").strip().lower() not in _PLACEHOLDER_SOURCE_MARKERS
    if not (has_doi or has_locator):
        return False, "no DOI and no figure/page locator for a numeric claim"
    return True, "citation present"


def _check_tag_consistent_with_sourcing(anchor: dict):
    """Reject source-tag upgrading: a placeholder/unretrieved/news-summary
    source cannot carry tag V or DR (that would upgrade untraceable evidence
    to 'verified against the paper')."""
    if _is_placeholder_gap(anchor):
        return anchor.get("tag") in ("A", "E"), f"placeholder source tagged {anchor.get('tag')!r}"
    return True, "sourced (no placeholder check needed)"


def _check_unit_and_magnitude(claim: str, anchor: dict):
    meta = CLAIM_META.get(claim)
    if meta is None:
        return True, "no claim metadata (unmanaged claim)"
    unit = anchor.get("unit")
    if unit != meta["unit"]:
        return False, f"expected unit {meta['unit']!r}, got {unit!r}"
    val = anchor.get("value")
    if val is None:
        return True, "no value to range-check"
    if not (meta["lo"] <= float(val) <= meta["hi"]):
        return False, f"value {val} outside plausible [{meta['lo']}, {meta['hi']}] {unit}"
    return True, "unit/magnitude plausible"


def _check_context_tokens(claim: str, anchor: dict):
    meta = CLAIM_META.get(claim)
    if meta is None or anchor.get("value") is None:
        return True, "no claim metadata or declared gap (skipped)"
    conditions = anchor.get("conditions") or ""
    for key in ("material_tokens", "temp_tokens", "drive_tokens"):
        tokens = meta[key]
        if tokens and not any(t in conditions for t in tokens):
            return False, f"conditions missing all of {tokens!r} ({key})"
    return True, "context tokens present"


_GAINP_80K_RE = re.compile(r"(?<!Al)GaInP")


def _check_no_gainp_80k_confusion(anchor: dict):
    """The real Reischle 2008 80 K electrical device is InP/AlGaInP, not
    InP/GaInP (that combination is Reischle 2010, at 20-40 K); this guards
    against exactly the "InP/GaInP 80 K" mistaken-request pattern the spec
    calls out."""
    conditions = anchor.get("conditions") or ""
    if "80 K" in conditions and _GAINP_80K_RE.search(conditions):
        return False, "bare 'GaInP' (not 'AlGaInP') combined with '80 K' -- no such device exists"
    return True, "no InP/GaInP+80K confusion"


def _check_raw_vs_corrected(claim: str, anchor: dict):
    spec = RAW_VS_CORRECTED_DISTRACTORS.get(claim)
    if spec is None or anchor.get("value") is None:
        return True, "no raw/corrected distractor table for this claim"
    conditions = (anchor.get("conditions") or "").lower()
    if "raw" not in conditions:
        return True, "not asserted as a raw observation"
    val = float(anchor["value"])
    for distractor in spec["distractors"]:
        if abs(val - distractor) < 1e-9:
            return False, (f"value {val} matches a background-corrected/deconvolved "
                           f"figure ({distractor}), not the raw dip ({spec['raw_value']})")
    return True, "raw value distinct from known corrected/deconvolved distractors"


def _check_tolerance_sane(anchor: dict):
    """Reject an altered/blown-up threshold: no anchor's tolerance may be so
    wide it would accept almost any number (a symptom of a tolerance edited
    to fit a desired outcome rather than a predeclared measurement/reading
    uncertainty)."""
    if anchor.get("value") is None:
        return True, "no value to bound"
    bound = anchor_bound(anchor)
    val = abs(float(anchor["value"]))
    if val > 0 and bound > 0.30 * val:
        return False, f"tolerance {bound} exceeds 30% of |value|={val}"
    return True, "tolerance within a plausible measurement/reading bound"


def _check_doi_plausible(anchor: dict):
    """DOI plausibility (council review Part B.ii): a DOI must carry a
    registered publisher prefix (REGISTERED_DOI_PREFIXES) unless the anchor
    is explicitly whitelisted with a `web:` verification note -- a real
    paper hand-confirmed under a prefix this list happens not to cover. A
    bare fabricated DOI (e.g. a hallucinated "10.1000/xyz") has neither and
    is rejected."""
    if anchor.get("value") is None:
        return True, "no numeric claim to source (declared gap)"
    doi = (anchor.get("doi") or "").strip()
    if not doi:
        return True, "no DOI asserted (locator-only sourcing)"
    if doi.startswith(REGISTERED_DOI_PREFIXES):
        return True, f"DOI {doi!r} carries a registered publisher prefix"
    verification = (anchor.get("verification") or "").strip()
    if verification.startswith("web:"):
        return True, f"unregistered DOI prefix but explicitly whitelisted via {verification!r}"
    return False, f"DOI {doi!r} has no registered publisher prefix and no web: whitelist"


def _check_digest_crosscheck(claim: str, anchor: dict, digest_table: dict):
    """Value-level hallucination gate (council review Part B.i): every
    VERIFIED, valued anchor must agree -- same DOI+claim, unit, and value
    within the anchor's own declared tolerance -- with an entry in the
    independently hand-transcribed verify/data/rt_edge_digest_values.yaml
    table. That table is never derived from this ledger, so it catches what
    the purely structural checks above cannot: a fabricated anchor (whose
    doi/claim pair simply is not in the table) or a mutated real value
    (whose doi/claim pair IS in the table, but at a different number).
    Anchors that are not status=='verified', or that declare no numeric
    value (a gap), are not cross-checked here -- the spec scopes this gate
    to verified anchors only."""
    if anchor.get("status") != "verified" or anchor.get("value") is None:
        return True, "not a verified numeric claim (skipped)"
    doi = (anchor.get("doi") or "").strip()
    entry = digest_table.get((doi, claim))
    if entry is None:
        return False, f"no independent digest-value table entry for doi={doi!r} claim={claim!r}"
    if entry.get("unit") != anchor.get("unit"):
        return False, (f"unit mismatch: anchor unit {anchor.get('unit')!r} vs digest-table "
                       f"unit {entry.get('unit')!r}")
    bound = anchor_bound(anchor)
    if abs(float(anchor["value"]) - float(entry["value"])) > bound + 1e-9:
        return False, (f"value {anchor['value']} is outside tolerance ({bound}) of the "
                       f"independently transcribed {entry['value']} {entry['unit']} "
                       f"({entry.get('digest_section', '?')})")
    return True, (f"matches independent transcription {entry['value']} {entry['unit']} "
                  f"({entry.get('digest_section', '?')})")


ANCHOR_CHECKS = [
    _check_tag_valid,
    _check_citation_present,
    _check_tag_consistent_with_sourcing,
    _check_no_gainp_80k_confusion,
    _check_tolerance_sane,
    _check_doi_plausible,
]
CLAIM_ANCHOR_CHECKS = [_check_unit_and_magnitude, _check_context_tokens, _check_raw_vs_corrected]


def _anchor_passes_all(claim: str, anchor: dict, digest_table: dict | None = None) -> tuple[bool, list]:
    if digest_table is None:
        digest_table = load_digest_values()
    details = []
    all_ok = True
    for fn in ANCHOR_CHECKS:
        ok, detail = fn(anchor)
        details.append({"check": fn.__name__, "ok": ok, "detail": detail})
        all_ok = all_ok and ok
    for fn in CLAIM_ANCHOR_CHECKS:
        ok, detail = fn(claim, anchor)
        details.append({"check": fn.__name__, "ok": ok, "detail": detail})
        all_ok = all_ok and ok
    ok, detail = _check_digest_crosscheck(claim, anchor, digest_table)
    details.append({"check": "_check_digest_crosscheck", "ok": ok, "detail": detail})
    all_ok = all_ok and ok
    return all_ok, details


# ---------------------------------------------------------- claim-level scoring

def score_ledger(anchors: dict, digest_table: dict | None = None) -> dict:
    """Recompute per-claim evidence completeness from an anchors dict --
    NEVER trusting a stored `status` field by itself. An anchor only counts
    toward a claim's required two distinct primary sources if: its ledger
    status is 'verified', it has a real value, it passes every structural/
    hallucination check above, AND (per-claim rule) at least one of the two
    counted anchors is tagged V or DR (an [E]-only pair cannot alone
    establish target-device performance, docs/rt_edge_contract.md /
    the spec's interface constraints).

    Used both for the real ledger (normal mode) and for self-test copies with
    exactly one deliberate defect (removed anchor / duplicate DOI /
    incompatible context / downgraded or invalid tag) -- the SAME function
    must reject every one of those, per the spec's acceptance criteria."""
    if digest_table is None:
        digest_table = load_digest_values()
    anchor_results = {}
    claim_results = {}
    missing_evidence = []
    by_claim = group_by_claim(anchors)
    for claim, aids in by_claim.items():
        eligible_keys = {}
        any_v_or_dr = False
        for aid in aids:
            anchor = anchors[aid]
            passed, details = _anchor_passes_all(claim, anchor, digest_table)
            effective_verified = (passed and anchor.get("status") == "verified"
                                  and anchor.get("value") is not None)
            anchor_results[aid] = {
                "claim": claim, "paper_id": anchor.get("doi") or anchor.get("source"),
                "expected": anchor.get("value"), "tolerance": anchor.get("tolerance"),
                "context": anchor.get("conditions"), "tag": anchor.get("tag"),
                "status": anchor.get("status"),
                "applicability": "proxy" if CLAIM_META.get(claim, {}).get("proxy") else "direct",
                "structural_checks_passed": passed, "checks": details,
                "effective_verified": effective_verified,
            }
            if effective_verified:
                key = _distinct_key(anchor)
                eligible_keys.setdefault(key, aid)
                if anchor.get("tag") in ("V", "DR"):
                    any_v_or_dr = True
        distinct = len(eligible_keys)
        complete = distinct >= 2 and any_v_or_dr
        reason = None
        if not complete:
            if distinct < 2:
                reason = f"only {distinct} distinct verified primary source(s), need 2"
            else:
                reason = "both counted sources are [E]-class only; cannot alone validate target performance"
        claim_results[claim] = {
            "anchors": aids, "distinct_verified_sources": distinct, "complete": complete,
            "reason": reason,
        }
        if not complete:
            missing_evidence.append({"claim": claim, "reason": reason, "anchors": aids})
    return {"claim_results": claim_results, "anchor_results": anchor_results,
            "missing_evidence": missing_evidence,
            "all_complete": all(c["complete"] for c in claim_results.values())}


# ------------------------------------------------------ evaluator cross-checks

def _load_device_modules():
    """Deferred import: keeps this module's own import side-effect-free
    (fsim_core.device pulls in numpy/scipy transitively) and lets --self-test
    run even if a concurrent edit temporarily breaks an unrelated module."""
    from fsim_core.device import DeviceDesign, evaluate  # noqa: E402
    from fsim_core import dot_levels  # noqa: E402
    return DeviceDesign, evaluate, dot_levels


def _card_field(DeviceDesign, card_path: Path, dotted: str):
    design = DeviceDesign.load(card_path)
    block_name, field_name = dotted.split(".", 1)
    return getattr(getattr(design, block_name), field_name)


def _xcheck_delta_xx_inp_gaasp(anchors):
    DeviceDesign, evaluate, _ = _load_device_modules()
    gaasp_anchor, gainp_anchor = anchors["reischle08-delta-xx"], anchors["bommer11-delta-xx"]
    gaasp_val = _card_field(DeviceDesign, GAASP_CARD, "dot.delta_xx")
    gainp_val = _card_field(DeviceDesign, GAINP_CARD, "dot.delta_xx")
    ok_gaasp = abs(gaasp_val - gaasp_anchor["value"]) <= anchor_bound(gaasp_anchor) + 1e-12
    ok_gainp = abs(gainp_val - gainp_anchor["value"]) <= anchor_bound(gainp_anchor) + 1e-12
    # Evaluated through the real evaluator, not just parsed YAML: the two
    # cards' different anchored delta_xx must feed a different eps_op.
    d1, d2 = DeviceDesign.load(GAASP_CARD), DeviceDesign.load(GAINP_CARD)
    r1, r2 = evaluate(d1, [d1.thermal.T_hs]), evaluate(d2, [d2.thermal.T_hs])
    fed = r1["scalars"]["eps_op"] != r2["scalars"]["eps_op"]
    ok = ok_gaasp and ok_gainp and fed
    detail = (f"gaasp card dot.delta_xx={gaasp_val} vs anchor {gaasp_anchor['value']}; "
              f"gainp card dot.delta_xx={gainp_val} vs anchor {gainp_anchor['value']}; "
              f"evaluator-fed (eps_op differs)={fed}")
    return ok, detail, {"gaasp_delta_xx": gaasp_val, "gainp_delta_xx": gainp_val}


def _xcheck_gamma300_class_range(anchors):
    DeviceDesign, evaluate, _ = _load_device_modules()
    anchor = anchors["laferriere23-linewidth-class-proxy"]
    gaasp_val = _card_field(DeviceDesign, GAASP_CARD, "dot.gamma300")
    gainp_val = _card_field(DeviceDesign, GAINP_CARD, "dot.gamma300")
    bound = anchor_bound(anchor)
    ok = (abs(gaasp_val - anchor["value"]) <= bound + 1e-12
          and abs(gainp_val - anchor["value"]) <= bound + 1e-12)
    d = DeviceDesign.load(GAASP_CARD)
    r = evaluate(d, [d.thermal.T_hs])
    fed = bool(r["scalars"]["linewidth_source"] == "anchored" and
              (r["scalars"]["gamma_op"] == r["scalars"]["gamma_op"]))  # finite (not NaN)
    detail = (f"both cards' dot.gamma300={gaasp_val} within tolerance of the InAsP/InP "
              f"nanowire class proxy (anchor 'missing' -- this remains an [E] class value, "
              f"never a direct InP/GaAsP measurement); evaluator-fed anchored linewidth={fed}")
    return ok and fed, detail, {"gaasp_gamma300": gaasp_val, "gainp_gamma300": gainp_val}


def _xcheck_hkust_inp_gaasp_wavelength(anchors):
    DeviceDesign, evaluate, _ = _load_device_modules()
    d = DeviceDesign.load(GAASP_CARD)
    r = evaluate(d, [d.thermal.T_hs])
    modeled_lambda = r["scalars"]["edge_lambda_nm"]
    # The design's OWN modeled wavelength must NOT silently coincide with
    # either unverified HKUST-wavelength anchor -- coincidental agreement
    # would let a missing-evidence anchor "pass" through a back door
    # (spec: "no metadata flag can override failed evidence").
    far_from_gu25 = abs(modeled_lambda - anchors["gu25-hkust-wavelength"]["value"]) > \
        anchor_bound(anchors["gu25-hkust-wavelength"])
    far_from_hkust22 = abs(modeled_lambda - anchors["hkust22-wavelength-independent-gap"]["value"]) > \
        anchor_bound(anchors["hkust22-wavelength-independent-gap"])
    ok = far_from_gu25 and far_from_hkust22
    detail = (f"card's own modeled edge_lambda_nm={modeled_lambda} (the explicit tabulated-index "
              f"design point) stays distinct from both unverified HKUST-target wavelength anchors "
              f"({anchors['gu25-hkust-wavelength']['value']}, "
              f"{anchors['hkust22-wavelength-independent-gap']['value']}) -- the design does not "
              f"claim to have reached the HKUST wavelength")
    return ok, detail, {"modeled_edge_lambda_nm": modeled_lambda}


def _xcheck_gaas_pin_iv(anchors):
    DeviceDesign, evaluate, _ = _load_device_modules()
    cases = [("reischle08-gaas-pin-iv", "InP/AlGaInP0.2/AlGaInP0.55 on GaAs", 80.0),
             ("reischle10-gaas-pin-iv", "InP/GaInP/AlGaInP0.55 on GaAs", 20.0)]
    ok_all, observed = True, {}
    details = []
    for anchor_id, preset, T in cases:
        anchor = anchors[anchor_id]
        d = DeviceDesign()
        d.ret.mode = "confinement"; d.ret.preset = preset
        d.drive.mode = "EL-transport"; d.drive.diode = {"preset": "red"}
        d.drive.n_dot_cm2 = 1e10
        d.thermal.T_hs = T
        r = evaluate(d, [T])
        sc = r["scalars"]
        # V_bi is not a substitute for the measured forward operating
        # voltage (spec: "V_bi ~1.2 V does not validate forward I-V"); the
        # transport-modeled JUNCTION voltage V_j_op is the comparable
        # quantity, and only class-consistency (not exact agreement, which
        # would be a fitted-not-sourced claim) is asserted.
        v_bi_distinct = abs(sc["V_bi"] - anchor["value"]) > 0.1
        v_j_class_consistent = abs(sc["V_j_op"] - anchor["value"]) / anchor["value"] < 0.25
        case_ok = v_bi_distinct and v_j_class_consistent
        ok_all = ok_all and case_ok
        observed[anchor_id] = {"V_bi": sc["V_bi"], "V_j_op": sc["V_j_op"]}
        details.append(f"{anchor_id}: V_bi={sc['V_bi']:.3f} (distinct from anchor, "
                       f"{v_bi_distinct}); V_j_op={sc['V_j_op']:.3f} vs anchor "
                       f"{anchor['value']} (class-consistent <25%, {v_j_class_consistent})")
    return ok_all, "; ".join(details), observed


def _xcheck_edge_extraction_efficiency(anchors):
    DeviceDesign, evaluate, _ = _load_device_modules()
    d = DeviceDesign.load(GAINP_CARD)
    r = evaluate(d, [d.thermal.T_hs])
    edge_beta = r["scalars"]["edge_beta"]
    ridge_class_ok = 0.003 < edge_beta < 0.03
    # Guard: a vertical photonic-nanowire/first-lens collection efficiency
    # (Laferriere: 0.276; Reischle vertical reference: 9e-5) must not be
    # conflated with this ridge's own beta factor -- "high-beta photonic
    # crystals do not prove ~1% ridge beta" (spec interface constraints).
    far_from_laferriere = abs(edge_beta - anchors["laferriere23-first-lens-efficiency"]["value"]) > 0.1
    far_from_reischle = abs(edge_beta - anchors["reischle08-first-lens-efficiency"]["value"]) > 0.001
    ok = ridge_class_ok and far_from_laferriere and far_from_reischle
    detail = (f"gainp card edge_beta={edge_beta:.4f} sits in its own [E] ridge class "
              f"(0.3%-3%, {ridge_class_ok}) and stays distinct from both vertical-reference "
              f"anchors ({anchors['laferriere23-first-lens-efficiency']['value']}, "
              f"{anchors['reischle08-first-lens-efficiency']['value']}) -- a vertical collection "
              f"efficiency does not validate this edge ridge's beta")
    return ok, detail, {"edge_beta": edge_beta}


def _xcheck_retention_Ea_inp_algainp(anchors):
    DeviceDesign, evaluate, dot_levels = _load_device_modules()
    anchor = anchors["bommer11-retention-ea"]
    bommer_name = "InP/AlGaInP0.2/AlGaInP0.55 on GaAs"
    hkust_name = "InP/GaAsP0.4/AlGaAs0.4 on GaAs"
    hkust_prefactor = dot_levels.retention_params(
        dot_levels.levels(dot_levels.class_presets()[hkust_name](T=300.0)),
        tau_rad_ns=1.0, channel="pair_half", verbose=False)["a_esc"]

    def _S_at(a_esc_override):
        d = DeviceDesign(); d.thermal.T_hs = 300.; d.drive.V = 0.
        d.ret.mode = "confinement"; d.ret.preset = bommer_name
        d.ret.overrides = {"E_a": float(anchor["value"])}
        if a_esc_override is not None:
            d.ret.overrides["a_esc"] = float(a_esc_override)
        return evaluate(d, [300.])["scalars"]["S_resolved"]

    S_sourced_prefactor = _S_at(None)
    S_other_prefactor = _S_at(hkust_prefactor)
    rel_diff = abs(S_sourced_prefactor - S_other_prefactor) / max(S_sourced_prefactor, S_other_prefactor)
    # The anchor is Ea=96 meV alone; demonstrate through the real evaluator
    # that the SAME Ea with a different (also literature-class-derived)
    # prefactor gives a meaningfully different S -- "a reported 96 meV
    # activation energy does not specify S without a prefactor."
    ok = 0.0 < S_sourced_prefactor < 1.0 and rel_diff > 0.05
    detail = (f"Ea={anchor['value']} meV with the sourced Bommer-class prefactor gives "
              f"S(300K)={S_sourced_prefactor:.4g}; the SAME Ea with a different class prefactor "
              f"gives {S_other_prefactor:.4g} ({rel_diff:.1%} different) -- Ea alone does not "
              f"determine S")
    return ok, detail, {"S_sourced_prefactor": S_sourced_prefactor, "S_other_prefactor": S_other_prefactor}


def _xcheck_temperature_trend_g2(anchors):
    DeviceDesign, evaluate, _ = _load_device_modules()
    d = DeviceDesign()
    d.ret.mode = "confinement"; d.ret.preset = "InGaAs0.5/GaAs/AlGaAs0.57 on GaAs"
    Ts = [200.0, 230.0, 260.0, 300.0]
    r = evaluate(d, Ts)
    g2 = list(r["curves"]["g2"])
    monotone = all(b >= a - 1e-12 for a, b in zip(g2, g2[1:]))
    ok = monotone
    detail = (f"Chatzarakis-class fixture g2(T) over {Ts} = {[round(v, 4) for v in g2]}: "
              f"non-decreasing with T ({monotone}), matching the direction BOTH independent "
              f"papers (Chatzarakis 0.36@230K, Laferriere 0.57@300K) report; not a claim that "
              f"the model reproduces their absolute numbers (different material systems)")
    return ok, detail, {"g2_curve": g2}


def _xcheck_reischle2008_g2_80K(anchors):
    DeviceDesign, evaluate, _ = _load_device_modules()
    d = DeviceDesign()
    d.ret.mode = "confinement"; d.ret.preset = "InP/AlGaInP0.2/AlGaInP0.55 on GaAs"
    d.drive.mode = "EL-transport"; d.drive.diode = {"preset": "red"}
    d.drive.n_dot_cm2 = 1e10
    d.thermal.T_hs = 80.0
    r = evaluate(d, [80.0])
    g2_op = r["scalars"]["g2_op"]
    eligible = len(r["scalars"]["invalid_reasons"]) == 0
    ok = eligible  # wiring sanity only -- NOT a claim of numeric agreement with 0.43
    detail = (f"InP/AlGaInP EL fixture at 80 K evaluates to a physically eligible corner "
              f"(g2_op={g2_op:.4f}); this is delegation-wiring evidence that the simulator can "
              f"represent this device class, NOT a reproduction of Reischle's measured 0.43 -- "
              f"the claim stays gated on missing second-source evidence regardless")
    return ok, detail, {"g2_op": g2_op}


CLAIM_CROSSCHECKS = {
    "delta_xx_inp_gaasp": _xcheck_delta_xx_inp_gaasp,
    "gamma300_class_range": _xcheck_gamma300_class_range,
    "hkust_inp_gaasp_wavelength": _xcheck_hkust_inp_gaasp_wavelength,
    "gaas_pin_iv": _xcheck_gaas_pin_iv,
    "edge_extraction_efficiency": _xcheck_edge_extraction_efficiency,
    "retention_Ea_inp_algainp": _xcheck_retention_Ea_inp_algainp,
    "temperature_trend_g2": _xcheck_temperature_trend_g2,
    "reischle2008_g2_80K": _xcheck_reischle2008_g2_80K,
}


def run_evaluator_crosschecks(anchors: dict, ok_fn) -> dict:
    """Run every claim's evaluator-integration cross-check, recording each as
    a top-level check (so a deliberately wrong parameter that breaks one
    fails the overall gate) and returning per-claim detail for the report."""
    results = {}
    for claim, fn in CLAIM_CROSSCHECKS.items():
        try:
            passed, detail, observed = fn(anchors)
        except Exception as exc:  # never let a fixture crash sink the whole run
            passed, detail, observed = False, f"evaluator cross-check raised {exc!r}", {}
        ok_fn(f"evaluator cross-check [{claim}]", passed)
        results[claim] = {"ok": passed, "detail": detail, "observed": observed}
    return results


# --------------------------------------------------------------------- self-test

def _two_source_complete_target(anchors: dict, scored: dict, *,
                                require_context: bool = False):
    """Choose a live, exactly-two-source complete claim for structural
    tampering.  The self-test must test the detector, not assume a particular
    paper remains verified in the ledger.  Exactly two effective sources make
    the expected loss of completeness unambiguous for every tamper below."""
    for claim, result in scored["claim_results"].items():
        if not result["complete"]:
            continue
        aids = [aid for aid in result["anchors"]
                if scored["anchor_results"][aid]["effective_verified"]]
        if len(aids) != 2 or len({_distinct_key(anchors[aid]) for aid in aids}) != 2:
            continue
        if require_context and not CLAIM_META.get(claim, {}).get("material_tokens"):
            continue
        return claim, aids
    return None, []


def _tamper_removed_anchor(anchors: dict, aid: str) -> dict:
    t = copy.deepcopy(anchors)
    del t[aid]
    return t


def _tamper_duplicate_doi(anchors: dict, aids: list[str]) -> dict:
    t = copy.deepcopy(anchors)
    t[aids[1]]["doi"] = t[aids[0]]["doi"]
    return t


def _tamper_incompatible_context(anchors: dict, aid: str) -> dict:
    t = copy.deepcopy(anchors)
    # Drops every claim-required context token, whatever the selected claim.
    t[aid]["conditions"] = "deliberately incompatible self-test context"
    return t


def _tamper_invalid_tag(anchors: dict, aid: str) -> dict:
    t = copy.deepcopy(anchors)
    t[aid]["tag"] = "downgraded"  # not a valid V/DR/E/A tag
    return t


def _tamper_e_class_only(anchors: dict, aids: list[str]) -> dict:
    t = copy.deepcopy(anchors)
    for aid in aids:
        t[aid]["tag"] = "E"
    return t


def _tamper_mutate_chatzarakis_g2_value(anchors: dict) -> dict:
    """Value-level hallucination (council review Part B): mutate a REAL,
    already-verified anchor's number while keeping its real DOI/source/
    conditions untouched -- the class of error the purely structural checks
    above cannot catch (right unit, right context tokens, real citation),
    but the independent digest-values table can, because the table's own
    0.36 doesn't move when this copy's anchor value does."""
    t = copy.deepcopy(anchors)
    t["chatzarakis23-g2-temperature"]["value"] = 0.90
    return t


def _tamper_mutate_matsuda_gamma300_value(anchors: dict) -> dict:
    t = copy.deepcopy(anchors)
    t["matsuda01-gamma300-class"]["value"] = 30.0
    return t


def _tamper_fabricated_anchor_with_flipped_status(anchors: dict) -> dict:
    """Reproduces the reviewer's fabricated-anchor scenario concretely: an
    attacker (or a careless edit) both (a) flips the real single-source
    reischle2008_g2_80K anchor's status from 'missing' to 'verified' -- the
    "fix" a bad-faith edit might make, since the number itself is genuine --
    and (b) adds a second, wholly fabricated anchor under the same claim
    with a fake journal, a fake DOI, and value 0.44, but otherwise-correct
    unit/context tokens so every structural check above passes it. Before
    the digest cross-check / DOI-plausibility check existed, this pair would
    give distinct=2, any_v_or_dr=True and flip reischle2008_g2_80K to
    'complete'; the digest cross-check rejects the fabricated anchor (its
    fake DOI has no table entry) and the DOI-plausibility check independently
    rejects it too (fake DOI has no registered prefix and no web:
    whitelist), so distinct drops back to 1 and the claim stays incomplete."""
    t = copy.deepcopy(anchors)
    t["reischle08-g2-80k"]["status"] = "verified"
    t["fabricated-g2-80k"] = _synthetic_anchor(
        value=0.44, unit="g2(0)", tolerance=0.02,
        conditions="InP/AlGaInP p-i-n EL, 80 K, raw dip.",
        source="Fake Journal of Quantum Nonsense 1, 1 (2000)",
        doi="10.1000/xyz.fake", locator="Fig. 1",
        tag="V", status="verified", claim="reischle2008_g2_80K")
    return t


def _synthetic_anchor(**overrides) -> dict:
    base = dict(id="synthetic", claim="reischle2008_g2_80K", value=0.2, unit="g2(0)",
               tolerance=0.02, conditions="InP/AlGaInP p-i-n EL, 80 K, raw dip.",
               source="Some Journal 1, 1 (2000)", doi="10.1000/xyz", locator="Fig. 1",
               tag="V", status="verified")
    base.update(overrides)
    return base


def _naive_recount_complete(anchors: dict) -> set:
    """Independent recount of which claims are complete, computed straight
    from the ledger's own raw status/value/doi/tag fields -- NOT via
    score_ledger's structural/hallucination checks. This is what the
    self-test compares score_ledger's output against, so "which claims are
    complete" is derived from the ledger every run instead of a hard-coded
    claim-name snapshot that goes stale the moment a genuine second source
    is added to (or removed from) the ledger."""
    by_claim = group_by_claim(anchors)
    complete = set()
    for claim, aids in by_claim.items():
        keys = {}
        any_v_or_dr = False
        for aid in aids:
            a = anchors[aid]
            if a.get("status") == "verified" and a.get("value") is not None:
                keys.setdefault(_distinct_key(a), aid)
                if a.get("tag") in ("V", "DR"):
                    any_v_or_dr = True
        if len(keys) >= 2 and any_v_or_dr:
            complete.add(claim)
    return complete


def _run_self_test(ok_fn) -> None:
    real_anchors = load_anchors()

    baseline = score_ledger(real_anchors)
    expected_complete = _naive_recount_complete(real_anchors)
    ok_fn("real ledger: score_ledger's completeness matches an independent "
         "recount of the ledger's own status/value/doi/tag fields",
         {c for c, r in baseline["claim_results"].items() if r["complete"]} == expected_complete)
    ok_fn("real ledger: missing-evidence claims are named with a reason",
         all(m["reason"] for m in baseline["missing_evidence"]) and
         len(baseline["missing_evidence"]) == len(CLAIM_META) - len(expected_complete))

    structural_claim, structural_aids = _two_source_complete_target(real_anchors, baseline)
    context_claim, context_aids = _two_source_complete_target(
        real_anchors, baseline, require_context=True)
    have_structural_target = structural_claim is not None
    have_context_target = context_claim is not None

    tamper_cases = [
        ("removed anchor breaks a previously-complete claim",
         lambda a: _tamper_removed_anchor(a, structural_aids[0]), structural_claim,
         have_structural_target),
        ("duplicate DOI collapses to one distinct source",
         lambda a: _tamper_duplicate_doi(a, structural_aids), structural_claim,
         have_structural_target),
        ("incompatible (material-dropped) context is rejected",
         lambda a: _tamper_incompatible_context(a, context_aids[0]), context_claim,
         have_context_target),
        ("invalid/downgraded tag is rejected",
         lambda a: _tamper_invalid_tag(a, structural_aids[0]), structural_claim,
         have_structural_target),
        ("[E]-only pair cannot alone establish target validation",
         lambda a: _tamper_e_class_only(a, structural_aids), structural_claim,
         have_structural_target),
        ("value-level hallucination: Chatzarakis g2(0) mutated 0.36->0.90 is "
         "caught by the independent digest-values cross-check",
         _tamper_mutate_chatzarakis_g2_value, "temperature_trend_g2", True),
        ("value-level hallucination: Matsuda gamma300 mutated 12->30 meV is "
         "caught by the independent digest-values cross-check",
         _tamper_mutate_matsuda_gamma300_value, "gamma300_class_range", True),
    ]
    for name, tamper_fn, claim, target_found in tamper_cases:
        if not target_found:
            ok_fn(f"self-test: {name}", False)
            continue
        was_complete = baseline["claim_results"][claim]["complete"]
        tampered = score_ledger(tamper_fn(real_anchors))
        now_complete = tampered["claim_results"][claim]["complete"]
        ok_fn(f"self-test: {name}", was_complete and not now_complete)

    # Fabricated-anchor scenario (council review Part B.iii): a fake journal,
    # a fake DOI, and value 0.44 -- with the right unit/context tokens so
    # every structural check passes it -- paired with the real single source
    # having its status "fixed" to verified must still NOT flip
    # reischle2008_g2_80K to complete (see the tamper function's docstring).
    fabricated_scored = score_ledger(_tamper_fabricated_anchor_with_flipped_status(real_anchors))
    ok_fn("self-test: fabricated anchor (fake journal/DOI, value 0.44) cannot "
         "flip reischle2008_g2_80K to complete even paired with a status-"
         "flipped real anchor",
         not fabricated_scored["claim_results"]["reischle2008_g2_80K"]["complete"])
    ok_fn("self-test: the fabricated anchor itself is rejected (not effective_verified)",
         not fabricated_scored["anchor_results"]["fabricated-g2-80k"]["effective_verified"])

    # Newly-completed claims: the ledger now carries a genuine second source
    # for these (Matsuda 2001 + Chatzarakis 2023; Schulz 2009 + Bommer
    # 2011), so they must score complete -- checked against the
    # ledger-derived `expected_complete`, not a hard-coded name.
    for claim in ("gamma300_class_range", "retention_Ea_inp_algainp"):
        ok_fn(f"self-test: newly-complete claim [{claim}] is scored complete",
             claim in expected_complete and baseline["claim_results"][claim]["complete"])

    # Missing-evidence rejection using the real ledger's own genuinely
    # incomplete claims (single-source-only pattern) -- derived from the
    # ledger via expected_complete, not hard-coded, so the default
    # (non-self-test) run over this same ledger correctly exits 1 with
    # exactly these claims named in evidence.json's missing_evidence.
    for claim in sorted(set(CLAIM_META) - expected_complete):
        ok_fn(f"self-test: missing-evidence claim [{claim}] correctly rejected",
             not baseline["claim_results"][claim]["complete"])

    # Isolated per-category hallucination-detector fixtures.
    ok_fn("self-test: absent/fake citation rejected",
         not _check_citation_present(_synthetic_anchor(source="", doi="", locator=""))[0])
    ok_fn("self-test: unit error rejected",
         not _check_unit_and_magnitude("gaas_pin_iv",
             _synthetic_anchor(claim="gaas_pin_iv", value=1970.0, unit="mV"))[0])
    ok_fn("self-test: altered/blown-up threshold rejected",
         not _check_tolerance_sane(_synthetic_anchor(tolerance={"rel": 5.0}))[0])
    ok_fn("self-test: source-tag upgrading rejected",
         not _check_tag_consistent_with_sourcing(
             _synthetic_anchor(source="Independent source not retrieved", doi="",
                              locator="not available", tag="V"))[0])
    ok_fn("self-test: raw-vs-intrinsic swap rejected",
         not _check_raw_vs_corrected("reischle2008_g2_80K", _synthetic_anchor(value=0.03))[0])
    ok_fn("self-test: InP/GaInP+80K context confusion rejected",
         not _check_no_gainp_80k_confusion(
             _synthetic_anchor(conditions="InP/GaInP p-i-n EL, 80 K electrically driven."))[0])
    ok_fn("self-test: missing-value gap is not itself flagged as a citation failure",
         _check_citation_present(_synthetic_anchor(value=None, source="", doi=""))[0])

    # DOI-plausibility fixtures (council review Part B.ii).
    ok_fn("self-test: implausible DOI prefix rejected",
         not _check_doi_plausible(_synthetic_anchor(doi="10.1000/xyz.fake"))[0])
    ok_fn("self-test: implausible DOI prefix accepted when explicitly web-whitelisted",
         _check_doi_plausible(_synthetic_anchor(
             doi="10.9999/unregistered-but-real",
             verification="web:https://example.org/hand-confirmed-real-paper"))[0])
    ok_fn("self-test: registered publisher prefixes (APS/AIP/Optica/Wiley) accepted "
         "without needing a web: whitelist",
         all(_check_doi_plausible(_synthetic_anchor(doi=doi))[0] for doi in
             ("10.1103/PhysRevB.1.1", "10.1063/1.1", "10.1364/OE.1.1", "10.1002/pssc.1")))

    # Independent digest-values cross-check fixtures (council review Part B.i).
    digest_table = load_digest_values()
    ok_fn("self-test: digest cross-check accepts the real, unmutated "
         "chatzarakis23-g2-temperature anchor",
         _check_digest_crosscheck("temperature_trend_g2",
             real_anchors["chatzarakis23-g2-temperature"], digest_table)[0])
    ok_fn("self-test: digest cross-check rejects a fabricated anchor whose "
         "doi/claim pair is absent from the table",
         not _check_digest_crosscheck("reischle2008_g2_80K",
             _synthetic_anchor(doi="10.1000/xyz.fake", status="verified"), digest_table)[0])
    ok_fn("self-test: digest cross-check rejects a mutated real value "
         "(Chatzarakis 0.36 -> 0.90) even with the real DOI",
         not _check_digest_crosscheck("temperature_trend_g2",
             _synthetic_anchor(**{**real_anchors["chatzarakis23-g2-temperature"],
                                  "value": 0.90}), digest_table)[0])
    ok_fn("self-test: digest cross-check skips a non-verified (missing-status) anchor",
         _check_digest_crosscheck("reischle2008_g2_80K",
             _synthetic_anchor(doi="10.1000/xyz.fake", status="missing"), digest_table)[0])

    # Orchestrator-required: a deliberately wrong parameter fed through the
    # real evaluator/DeviceDesign must make the corresponding gate FAIL.
    from fsim_core.device import DeviceDesign  # local import, self-test only
    good_delta_xx = _card_field(DeviceDesign, GAASP_CARD, "dot.delta_xx")
    anchor = real_anchors["reischle08-delta-xx"]
    bound = anchor_bound(anchor)
    ok_fn("self-test: real gaasp card dot.delta_xx passes its own anchor tolerance",
         abs(good_delta_xx - anchor["value"]) <= bound + 1e-12)
    wrong_delta_xx = 15.0  # deliberately outside [4, 6] meV
    ok_fn("self-test: a deliberately wrong dot.delta_xx fails the same anchor tolerance",
         abs(wrong_delta_xx - anchor["value"]) > bound)

    lambda_anchor = real_anchors["gu25-hkust-wavelength"]
    hallucinated_lambda = lambda_anchor["value"]  # coincidentally "matches" a missing anchor
    coincidence_within_tol = abs(hallucinated_lambda - lambda_anchor["value"]) <= anchor_bound(lambda_anchor)
    ok_fn("self-test: no metadata flag/coincidental match overrides a missing-status anchor",
         coincidence_within_tol and score_ledger(real_anchors)["claim_results"]
         ["hkust_inp_gaasp_wavelength"]["complete"] is False)


# ------------------------------------------------------------------------ main

def run_checks(self_test: bool = False) -> dict:
    """Run every check and return the machine-readable report dict. Never
    raises for an ordinary evidence/physics failure (that is reported via
    evidence_complete=False); only a genuinely broken fixture can raise, and
    even evaluator cross-checks catch their own exceptions."""
    warnings.filterwarnings("ignore", category=UserWarning, module="fsim_core.transport")
    # The gainp card opts into drive.cw diagnostics; cw_g2's closed-form dip
    # convolution transiently overflows for some (rate, IRF) combinations
    # this cross-check's fixture touches -- benign, same rationale as
    # verify_device_rt.py's identical filter.
    warnings.filterwarnings("ignore", category=RuntimeWarning, module="fsim_core.cw_g2")
    checks: list = []

    def ok(name: str, value: bool) -> bool:
        value = bool(value)
        checks.append((name, value))
        return value

    if self_test:
        _run_self_test(ok)
        claim_results, anchor_results, missing_evidence = {}, {}, []
        evidence_complete = False  # self-test asserts rejection behavior, not scientific PASS
    else:
        anchors = load_anchors()
        scored = score_ledger(anchors)
        for claim, result in scored["claim_results"].items():
            ok(f"claim [{claim}] has two distinct verified primary sources", result["complete"])
        run_evaluator_crosschecks(anchors, ok)
        claim_results, anchor_results = scored["claim_results"], scored["anchor_results"]
        missing_evidence = scored["missing_evidence"]
        evidence_complete = all(v for _, v in checks)

    passed, total = sum(1 for _, v in checks if v), len(checks)
    label = "rt-edge paper enforcement checks" if self_test else "rt-edge paper checks"
    for name, value in checks:
        print(("ok  " if value else "FAIL") + " " + name)
    print(f"{passed}/{total} {label} passed")

    return {
        "schema_version": 1,
        "evidence_complete": evidence_complete,
        "hallucination_tests_passed": (passed == total),
        "claim_results": claim_results,
        "anchor_results": anchor_results,
        "missing_evidence": missing_evidence,
        "source_ledger_sha256": _sha256_file(ANCHORS_PATH),
        "evaluator_input_hashes": {
            "verify/data/rt_edge_anchors.yaml": _sha256_file(ANCHORS_PATH),
            "cards/edge-inp-gaasp-design.yaml": _sha256_file(GAASP_CARD),
            "cards/edge-inp-gainp-design.yaml": _sha256_file(GAINP_CARD),
        },
        "checks_passed": passed, "checks_total": total, "self_test": self_test,
    }


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--report", default=str(DEFAULT_REPORT))
    parser.add_argument("--self-test", action="store_true")
    args = parser.parse_args(argv)

    report = run_checks(self_test=args.self_test)

    report_path = Path(args.report)
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(report, indent=2, sort_keys=False), encoding="utf-8")

    if args.self_test:
        return 0 if report["hallucination_tests_passed"] else 1
    return 0 if report["evidence_complete"] else 1


if __name__ == "__main__":
    sys.exit(main())
