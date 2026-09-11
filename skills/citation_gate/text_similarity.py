#!/usr/bin/env python3
# Copyright (c) 2026 Cheng-I Wu. Licensed CC-BY-NC-4.0 (non-commercial
# research use). Adapted for qd-photon-sim from academic-research-skills
# (v3.21.2) scripts/_text_similarity.py. See THIRD_PARTY.md.
#
# NOTE ON SCOPE: this is a minimal reimplementation of the parts of the
# upstream module that scripts/crossref_client.py depends on (title
# normalization + a Levenshtein-style similarity ratio + the 0.70 match
# threshold). _text_similarity.py itself is not one of the vendored files
# for this task, so the dotted-acronym and CJK-title normalization passes
# (irrelevant to this project's ASCII physics-journal citations) were
# intentionally left out rather than carried over unread.
"""Title-similarity helper shared by the vendored resolver clients.

`_similarity` mirrors the upstream algorithm's core: lowercase, strip
punctuation to whitespace, collapse whitespace, then a difflib
SequenceMatcher ratio. `_TITLE_SIMILARITY_THRESHOLD` (0.70) is the same
match/no-match cutoff the upstream clients use.
"""
from __future__ import annotations

import string
from difflib import SequenceMatcher

_PUNCT_TRANSLATION = str.maketrans({c: " " for c in string.punctuation})

# Per upstream protocol: title-similarity threshold for a "matched" verdict.
_TITLE_SIMILARITY_THRESHOLD = 0.70

# Shared retry budget for the resolver clients (429 backoff).
_BACKOFF_SECONDS = 2.0
_MAX_RETRIES = 3


def normalize_title(s: str) -> str:
    """Case-insensitive, punctuation-stripped, whitespace-collapsed form."""
    cleaned = (s or "").lower().translate(_PUNCT_TRANSLATION)
    return " ".join(cleaned.split())


def similarity(a: str, b: str) -> float:
    """Levenshtein-style ratio (via difflib) of the normalized titles."""
    return SequenceMatcher(None, normalize_title(a), normalize_title(b)).ratio()


def titles_match(a: str, b: str) -> bool:
    """True iff both titles are non-empty and clear the match threshold."""
    if not a or not b:
        return False
    return similarity(a, b) >= _TITLE_SIMILARITY_THRESHOLD
