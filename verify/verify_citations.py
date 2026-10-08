"""Citation existence check for the project's evidence ledgers
(including verify/data/nitride_nanowire_anchors.yaml): every anchor's
`doi_or_url` is resolved against Crossref / Semantic Scholar / arXiv
(via skills/citation_gate) and cross-checked against the anchor's
`citation` string.

This is a (T) source-transcription check in the README.md five-way split:
it confirms the cited paper *exists* and that title/year agree with the
citation string. It does not re-grade evidence quality (see
skills/source-verification/SKILL.md) and it never edits a ledger -- an
anchor that fails is reported, not "fixed".

Class discipline: no re-verification of fsim_core physics here; this
piece only exercises skills/citation_gate against literal ledger text.

Usage:
    python verify/verify_citations.py            # networked first run
    python verify/verify_citations.py --offline  # cache-only, no network

The cache is verify/data/citation_verification_cache.json (committed
after the first networked run so --offline reproduces the same result).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from skills.citation_gate import VerificationCache, to_ascii, verify_ledger

LEDGERS = [
    ROOT / "verify" / "data" / "nitride_cavity_anchors.yaml",
    ROOT / "verify" / "data" / "nitride_geometry_stark_anchors.yaml",
    ROOT / "verify" / "data" / "nitride_nanowire_anchors.yaml",
]
CACHE_PATH = ROOT / "verify" / "data" / "citation_verification_cache.json"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--offline", action="store_true",
        help="Serve only from the committed cache; a miss fails ('not in cache'). No network call is made.",
    )
    args = parser.parse_args(argv)

    cache = VerificationCache(CACHE_PATH)

    total_passed = 0
    total_count = 0
    for ledger_path in LEDGERS:
        with ledger_path.open("r", encoding="utf-8") as f:
            doc = yaml.safe_load(f)
        anchors = (doc or {}).get("anchors") or {}
        rows, passed, count = verify_ledger(anchors, cache=cache, offline=args.offline)
        for anchor_id, ok, identifier_display, detail in rows:
            # fix-2 finding 11: a recorded gap (null identifier,
            # evidence_status=missing) auto-passes (ok=True, see
            # skills/citation_gate/ledger.py verify_anchor) but is not a
            # DOI-resolved PASS -- print it as GAP so the gate is honest
            # about what was actually checked, while still exiting 0 for
            # every recorded gap (ok stays True; only the printed label
            # changes).
            if ok and "recorded gap" in detail:
                outcome = "GAP"
            else:
                outcome = "PASS" if ok else "FAIL"
            print(f"{anchor_id} | {identifier_display} | {outcome} | {to_ascii(detail)}")
        total_passed += passed
        total_count += count

    # Audit C4 regression pin (no network: every fixture below returns
    # before any resolver call). GAP widening for figure_reading anchors:
    # null identifier + no parseable title -> recorded gap (ok, printed GAP)
    # for evidence_status in {missing, figure_reading}; still a FAIL for
    # full_text / abstract_only / absent status (an anchor claiming the text
    # was read must carry an identifier or a title).
    from skills.citation_gate.ledger import verify_anchor
    _untitled = "Author et al., Some Journal 1, 2 (2023)"
    for status, want_gap in (("missing", True), ("figure_reading", True), ("full_text", False),
                             ("abstract_only", False), (None, False)):
        a = {"citation": _untitled, "doi_or_url": None}
        if status is not None:
            a["evidence_status"] = status
        ok, _disp, detail = verify_anchor("c4_fixture", a, cache=None, offline=True, clients={})
        good = (ok and "recorded gap" in detail) if want_gap else (not ok and "recorded gap" not in detail)
        print(f"c4_fixture_gap_rule_{status} | null | {'PASS' if good else 'FAIL'} | "
              f"{'GAP' if want_gap else 'FAIL'} expected; got ok={ok}")
        total_passed += int(good)
        total_count += 1
    # A V-tagged figure_reading anchor with no identifier must still FAIL.
    ok, _disp, detail = verify_anchor("c4_fixture_v", {"citation": _untitled, "doi_or_url": None,
                                      "evidence_status": "figure_reading", "tag": "V"},
                                      cache=None, offline=True, clients={})
    good = (not ok) and "recorded gap" not in detail
    print(f"c4_fixture_gap_rule_figure_reading_tag_V | null | {'PASS' if good else 'FAIL'} | FAIL expected; got ok={ok}")
    total_passed += int(good)
    total_count += 1

    print(f"{total_passed}/{total_count} citation checks passed")
    return 0 if total_passed == total_count else 1


if __name__ == "__main__":
    sys.exit(main())
