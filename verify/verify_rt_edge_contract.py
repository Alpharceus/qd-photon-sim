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
import re
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
ANCHORS_PATH = ROOT / "verify" / "data" / "rt_edge_anchors.yaml"
CONTRACT_DOC = ROOT / "docs" / "rt_edge_contract.md"
README_DOC = ROOT / "README.md"

# The five validation classes README.md's "Validation record" section must
# sort every check into (finding 10, peer-review-triage.md).
VALIDATION_CLASS_LETTERS = ("N", "T", "C", "M", "P")

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


def _readme_validation_record_section(text: str) -> str:
    """The body of README.md's "## Validation record" section only, up to
    (not including) the next "## " heading -- so class parsing below never
    bleeds into unrelated sections (e.g. the "## Quickstart" verify_fsim.py
    command line)."""
    m = re.search(r"## Validation record\b", text)
    if not m:
        return ""
    rest = text[m.end():]
    nxt = re.search(r"\n## ", rest)
    return rest[: nxt.start()] if nxt else rest


_CLASS_BULLET_RE = re.compile(r"- \*\*\(([NTCMP])\)[^\n]*")
# Matches "<name>.py" wherever it appears in backticks/prose (verify_*.py,
# gate_*.py, verify4.py, audit_physics.py -- not just the "verify_" prefix).
_SCRIPT_NAME_RE = re.compile(r"\b([A-Za-z_][A-Za-z0-9_]*)\.py\b")


def _readme_validation_class_list_block(text: str) -> str:
    """Just the five "- **(X) ..." class-definition bullets, stopping at the
    first blank line after the last one -- so the Arm table and the prose
    that follows it (which also use "(C)"/"(M)"/"(T)" as table-cell labels)
    are never swept into the last class's block and misread as more
    script-name assignments."""
    body = _readme_validation_record_section(text)
    matches = list(_CLASS_BULLET_RE.finditer(body))
    if not matches:
        return ""
    tail = body[matches[-1].start():]
    blank = re.search(r"\n[ \t]*\n", tail)
    end = matches[-1].start() + (blank.start() if blank else len(tail))
    return body[matches[0].start():end]


def _readme_validation_class_sections(text: str) -> dict:
    """Split the class-list block into per-class text blocks keyed by class
    letter (N/T/C/M/P), using the "- **(X) ..." class bullets as boundaries
    -- everything from one bullet up to the next belongs to that class
    (multi-line bullets included)."""
    block = _readme_validation_class_list_block(text)
    matches = list(_CLASS_BULLET_RE.finditer(block))
    sections: dict = {}
    for i, m in enumerate(matches):
        letter = m.group(1)
        start = m.end()
        end = matches[i + 1].start() if i + 1 < len(matches) else len(block)
        sections[letter] = block[start:end]
    return sections


def _readme_validation_record_ok() -> tuple:
    """(N)/(T)/(C)/(M)/(P) all present as headers, no <name>.py script
    claimed by more than one class, and every verify/*.py script (including
    audit_physics.py) named exactly once in the class list -- finding 10."""
    text = README_DOC.read_text(encoding="utf-8") if README_DOC.exists() else ""
    headers_present = all(f"({c})" in text for c in VALIDATION_CLASS_LETTERS)
    sections = _readme_validation_class_sections(text)
    all_classes_found = set(sections) == set(VALIDATION_CLASS_LETTERS)
    seen: dict = {}
    counts: dict = {}
    no_duplicate_script = True
    for letter, block in sections.items():
        for name in _SCRIPT_NAME_RE.findall(block):
            if name in seen and seen[name] != letter:
                no_duplicate_script = False
            seen[name] = letter
            counts[name] = counts.get(name, 0) + 1
    verify_dir = ROOT / "verify"
    all_scripts = {p.stem for p in verify_dir.glob("*.py")} if verify_dir.exists() else set()
    missing = sorted(all_scripts - set(seen))
    repeated = sorted(n for n in all_scripts & set(counts) if counts[n] != 1)
    all_scripts_named_once = not missing and not repeated
    ok = (headers_present and all_classes_found and no_duplicate_script
          and all_scripts_named_once)
    return ok, {"assignments": seen, "missing": missing, "repeated": repeated}


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
        # Optional `verification` key (digest/pdf:<path>/web:<url>/unverified,
        # for anchors citing a paper outside the fact-checked digest set): an
        # anchor that admits it is `unverified` cannot also claim status
        # verified -- that combination would be a self-contradicting anchor.
        verification = a.get("verification")
        if verification is not None:
            ok(f"anchor [{aid}] verification=='unverified' is not paired with status: verified",
               not (str(verification).strip() == "unverified" and a.get("status") == "verified"))

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

    readme_record_ok, readme_detail = _readme_validation_record_ok()
    detail = ""
    if not readme_record_ok:
        detail = (f" (missing: {readme_detail['missing']}, "
                  f"repeated/duplicated: {readme_detail['repeated']})")
    ok("README.md validation record has all five class headers "
       "(N)/(T)/(C)/(M)/(P), no <name>.py script claimed by more than one "
       f"class, and every verify/*.py script named exactly once{detail}",
       readme_record_ok)

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
