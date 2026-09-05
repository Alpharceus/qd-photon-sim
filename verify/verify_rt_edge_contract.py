"""Ledger/contract structural validator for the room-temperature edge-emitter
evidence contract (docs/rt_edge_contract.md + verify/data/rt_edge_anchors.yaml).

Several rt-edge specs and docs/rt_edge_contract.md itself refer to a
"contract validator" that checks the ledger and the contract document are
STRUCTURALLY well-formed, but until now no such script existed. This is NOT
the paper/hallucination-detection gate (that is verify/verify_rt_edge_papers.py,
which judges whether individual anchors are genuinely sourced, in-context,
and evaluator-wired) and it is NOT the design-card checker (that is
verify/verify_rt_edge_cards.py). This script only checks:

  1. every anchor in verify/data/rt_edge_anchors.yaml carries the full
     required key set (id, claim, value, unit, tolerance, conditions,
     source, doi, locator, tag, status);
  2. every anchor id is unique;
  3. every anchor's tag is one of V/DR/E/A and its status is one of
     verified/missing;
  4. the contract-documented claim set (the eight claims docs/rt_edge_
     contract.md and the rt-edge specs anchor against) is present in the
     ledger;
  5. every one of those claims has >= 2 distinct-DOI anchors whose status
     is "verified" (the ledger's own definition of "has two independent
     verified sources" -- docs/rt_edge_contract.md Evidence status);
  6. docs/rt_edge_contract.md carries every one of its seven required
     section headings.

CLI: python verify/verify_rt_edge_contract.py [--allow-missing]
Prints "ok "/"FAIL " per check and a final "N/N contract checks passed"
line; exits 0 iff every check passes. A claim that legitimately still lacks
a second independent verified source fails check 5 in the default (strict)
run; --allow-missing downgrades that specific failure to a printed warning
(and exits 0 if that is the only problem) for use while evidence-gathering
is incomplete -- it never relaxes checks 1-4 or 6.

Standalone, side-effect-free on import (all work happens under
`if __name__ == "__main__"` or inside run_checks()).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
ANCHORS_PATH = ROOT / "verify" / "data" / "rt_edge_anchors.yaml"
CONTRACT_DOC = ROOT / "docs" / "rt_edge_contract.md"

REQUIRED_ANCHOR_KEYS = {
    "id", "claim", "value", "unit", "tolerance", "conditions",
    "source", "doi", "locator", "tag", "status",
}
VALID_TAGS = {"V", "DR", "E", "A"}
VALID_STATUSES = {"verified", "missing"}

# The claim set every rt-edge spec/doc anchors against (docs/rt_edge_
# contract.md "Evidence status"; verify_rt_edge_papers.py CLAIM_META).
REQUIRED_CLAIMS = {
    "reischle2008_g2_80K",
    "temperature_trend_g2",
    "hkust_inp_gaasp_wavelength",
    "gaas_pin_iv",
    "edge_extraction_efficiency",
    "gamma300_class_range",
    "retention_Ea_inp_algainp",
    "delta_xx_inp_gaasp",
}

REQUIRED_HEADINGS = [
    "Block schema",
    "Units and normalization",
    "Temperature and rate mapping",
    "Aperture assumptions",
    "Lemma 1",
    "Acceptance gates",
    "Evidence status",
]


def load_anchor_list(path: Path = ANCHORS_PATH) -> list:
    doc = yaml.safe_load(path.read_text(encoding="utf-8"))
    return doc.get("anchors", []) if doc else []


def _distinct_doi(anchor: dict) -> str:
    return (anchor.get("doi") or "").strip()


def run_checks(allow_missing: bool = False) -> tuple[bool, list]:
    """Run every structural check and return (all_passed, checks), where
    checks is [(name, passed), ...] -- a claim downgraded by --allow-missing
    is printed as a warning and excluded from the pass/fail tally, matching
    verify_rt_edge_papers.py's own --self-test-vs-strict convention."""
    checks: list = []
    warnings: list = []

    def ok(name: str, value: bool) -> bool:
        value = bool(value)
        checks.append((name, value))
        return value

    anchors = load_anchor_list()

    for a in anchors:
        aid = a.get("id", "<missing id>")
        present_keys = set(a.keys())
        missing_keys = REQUIRED_ANCHOR_KEYS - present_keys
        ok(f"anchor [{aid}] has all required keys", not missing_keys)
        ok(f"anchor [{aid}] tag is valid (V/DR/E/A)", a.get("tag") in VALID_TAGS)
        ok(f"anchor [{aid}] status is valid (verified/missing)",
           a.get("status") in VALID_STATUSES)

    ids = [a.get("id") for a in anchors]
    ok("every anchor id is unique", len(ids) == len(set(ids)))

    claims_present = {a.get("claim") for a in anchors}
    missing_claims = REQUIRED_CLAIMS - claims_present
    ok(f"required claim set is present (missing: {sorted(missing_claims)})",
       not missing_claims)

    by_claim: dict = {}
    for a in anchors:
        by_claim.setdefault(a.get("claim"), []).append(a)

    for claim in sorted(REQUIRED_CLAIMS):
        claim_anchors = by_claim.get(claim, [])
        verified = [a for a in claim_anchors if a.get("status") == "verified"]
        distinct_dois = {d for d in (_distinct_doi(a) for a in verified) if d}
        sufficient = len(distinct_dois) >= 2
        name = (f"claim [{claim}] has >= 2 distinct-DOI verified sources "
                f"(has {len(distinct_dois)})")
        if not sufficient and allow_missing:
            warnings.append(name)
            print("WARN " + name + " -- allowed by --allow-missing "
                  "(evidence-gathering incomplete)")
            continue
        ok(name, sufficient)

    doc_text = CONTRACT_DOC.read_text(encoding="utf-8") if CONTRACT_DOC.exists() else ""
    ok("docs/rt_edge_contract.md exists", CONTRACT_DOC.exists())
    for heading in REQUIRED_HEADINGS:
        ok(f"docs/rt_edge_contract.md has heading [{heading}]", heading in doc_text)

    passed, total = sum(1 for _, v in checks if v), len(checks)
    for name, value in checks:
        print(("ok  " if value else "FAIL") + " " + name)
    print(f"{passed}/{total} contract checks passed")
    if warnings:
        print(f"(allowed missing-second-source claims: "
              f"{', '.join(w.split('[')[1].split(']')[0] for w in warnings)})")

    return passed == total, checks


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        description="Structural validator for the rt-edge evidence contract.")
    parser.add_argument("--allow-missing", action="store_true",
                        help="downgrade a still-single-sourced claim's failure "
                             "to a warning; never relaxes the other checks.")
    args = parser.parse_args(argv)

    all_passed, _ = run_checks(allow_missing=args.allow_missing)
    return 0 if all_passed else 1


if __name__ == "__main__":
    sys.exit(main())
