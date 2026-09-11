#!/usr/bin/env python3
# Copyright (c) 2026 Cheng-I Wu. Licensed under the Creative Commons
# Attribution-NonCommercial 4.0 International License (CC-BY-NC-4.0).
#
# This package vendors and adapts pieces of academic-research-skills
# (https://github.com/Imbad0202/academic-research-skills), tag v3.21.2,
# for non-commercial research use inside qd-photon-sim. See
# THIRD_PARTY.md at the repository root for the full file list, upstream
# URL, and licence text.
"""citation-gate: verify that a citation exists via DOI/Crossref,
Semantic Scholar, and arXiv, with a repo-local JSON cache.

Public API:
  - verify_citation(entry, *, cache=None, offline=False) -> outcome dict
      (skills.citation_gate.gate)
  - verify_ledger(anchors, *, cache=None, offline=False) -> (rows, passed, total)
      (skills.citation_gate.ledger) -- for the `{anchors: {...}}` ledger
      schema used by this project's verify/data/*.yaml files.
  - VerificationCache(path) (skills.citation_gate.verification_cache)

See SKILL.md for the `python -m skills.citation_gate <ledger.yaml>` CLI.
"""
from .gate import verify_citation
from .ledger import parse_citation, parse_identifier, to_ascii, verify_anchor, verify_ledger
from .verification_cache import VerificationCache

__all__ = [
    "verify_citation",
    "verify_ledger",
    "verify_anchor",
    "parse_citation",
    "parse_identifier",
    "to_ascii",
    "VerificationCache",
]
