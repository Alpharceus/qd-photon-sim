#!/usr/bin/env python3
# Copyright (c) 2026 Cheng-I Wu. Licensed CC-BY-NC-4.0 (non-commercial
# research use). Original code written for qd-photon-sim; not itself
# derived from an upstream file. See THIRD_PARTY.md.
"""CLI: `python -m skills.citation_gate <ledger.yaml> [--offline]`.

Verifies every anchor in a `{anchors: {...}}`-shaped YAML ledger (the
schema used by this project's verify/data/*.yaml evidence ledgers) and
prints one line per anchor, then `N/N citation checks passed`. Exits 0
iff every anchor passed.

The cache defaults to `verify/data/citation_verification_cache.json`
inside this repository (repo-relative, not a home-directory path).
`--offline` serves only from that cache -- a miss counts as a failure
("not in cache") and no network call is made.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import yaml

from .ledger import to_ascii, verify_ledger
from .verification_cache import VerificationCache

_DEFAULT_CACHE = Path(__file__).resolve().parents[2] / "verify" / "data" / "citation_verification_cache.json"


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(
        prog="python -m skills.citation_gate",
        description="Verify every citation in a ledger's `anchors:` mapping.",
    )
    parser.add_argument("ledger", help="Path to a YAML file with a top-level `anchors:` mapping.")
    parser.add_argument("--offline", action="store_true", help="Serve only from the cache; a miss fails.")
    parser.add_argument("--cache", default=str(_DEFAULT_CACHE), help="Cache file path (default: repo-local).")
    args = parser.parse_args(argv)

    with open(args.ledger, "r", encoding="utf-8") as f:
        doc = yaml.safe_load(f)
    anchors = (doc or {}).get("anchors") or {}

    cache = VerificationCache(args.cache)
    rows, passed, total = verify_ledger(anchors, cache=cache, offline=args.offline)
    for anchor_id, ok, identifier_display, detail in rows:
        outcome = "PASS" if ok else "FAIL"
        print(f"{anchor_id} | {identifier_display} | {outcome} | {to_ascii(detail)}")
    print(f"{passed}/{total} citation checks passed")
    return 0 if passed == total else 1


if __name__ == "__main__":
    sys.exit(main())
