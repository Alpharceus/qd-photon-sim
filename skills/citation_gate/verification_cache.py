#!/usr/bin/env python3
# Copyright (c) 2026 Cheng-I Wu. Licensed CC-BY-NC-4.0 (non-commercial
# research use). Vendored/adapted for qd-photon-sim from
# academic-research-skills (v3.21.2) scripts/verification_cache.py. See
# THIRD_PARTY.md for the upstream URL, tag, and licence.
#
# Adaptation notes (from the original docstring/behavior):
#   - Storage backend changed from the upstream tool's SQLite database
#     under a caller-home-directory cache folder to a single UTF-8 JSON
#     file inside this repository (default
#     `verify/data/citation_verification_cache.json`), per spec: the
#     cache must default to a repo-local path, not a home-directory
#     path, and must not require an environment variable.
#   - Cache key is (resolver_name, query_form) rather than
#     (citation_key, resolver_name, query_form): this project verifies
#     one citation per query form (a DOI or arXiv id is queried once and
#     shared across every anchor that cites the same paper), so the
#     citation_key dimension was dropped as redundant with query_form.
#   - The 90-day TTL / #541 stale-advisory machinery, SQLite WAL
#     concurrency handling, and the `stale_advisory_days()` env override
#     were not carried over (no cross-model/passport/hook code; see
#     spec). `entry_age_days` / `row_age_days` are omitted for the same
#     reason -- nothing in this project reads them.
#   - `get`/`put`/`invalidate` keep the same "malformed payload = miss,
#     never a hard error" contract as the original.
"""Persistent, repo-local, JSON-backed verification cache.

Cache key: (resolver_name, query_form).
  - resolver_name identifies which resolver/step produced the value
    (e.g. "crossref_doi", "crossref_title", "semantic_scholar", "arxiv_id",
    "arxiv_title").
  - query_form is the exact string the resolver was queried with (a DOI,
    an arXiv id, or a title).
Cache value: a small JSON-serializable dict (the resolver's projected
result, e.g. {"found": true, "title": ..., "year": ...}) plus a
verification_timestamp.

The cache file is written with a deterministic key order (sorted) so
repeated runs against an unchanged cache produce a byte-identical file.
"""
from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


class VerificationCache:
    """JSON-file-backed (resolver_name, query_form) -> response cache.

    Single-process use; the whole cache is read into memory on
    construction and rewritten (sorted, UTF-8) on every `put`.
    """

    def __init__(self, path: str | Path) -> None:
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)
        self._data: dict[str, dict[str, Any]] = self._load()

    def _load(self) -> dict[str, dict[str, Any]]:
        if not self.path.exists():
            return {}
        try:
            with self.path.open("r", encoding="utf-8") as f:
                loaded = json.load(f)
        except (json.JSONDecodeError, OSError, UnicodeDecodeError):
            # A corrupted cache file is a miss, not a hard error (mirrors
            # the upstream "malformed cache payload = miss" contract).
            return {}
        if not isinstance(loaded, dict):
            return {}
        return loaded

    @staticmethod
    def _key(resolver_name: str, query_form: str) -> str:
        return f"{resolver_name}|{query_form}"

    def get(self, resolver_name: str, query_form: str) -> dict[str, Any] | None:
        """Return the cached response dict, or None on a true cache miss
        (no row for this key). A cached row is always a hit -- there is no
        TTL here; the cache is meant to make the committed, offline
        reproduction exact, not to expire."""
        row = self._data.get(self._key(resolver_name, query_form))
        if row is None:
            return None
        response = row.get("response")
        if not isinstance(response, dict):
            return None
        return response

    def put(self, resolver_name: str, query_form: str, response: dict[str, Any]) -> None:
        """Store (or overwrite) the resolver response and persist the whole
        cache file immediately, sorted for a deterministic byte layout."""
        self._data[self._key(resolver_name, query_form)] = {
            "response": response,
            "verification_timestamp": datetime.now(timezone.utc).isoformat(),
        }
        self._write()

    def invalidate(self, resolver_name: str, query_form: str) -> None:
        """Drop one cached row. No-op if it is not present."""
        self._data.pop(self._key(resolver_name, query_form), None)
        self._write()

    def _write(self) -> None:
        with self.path.open("w", encoding="utf-8") as f:
            json.dump(self._data, f, sort_keys=True, indent=2, ensure_ascii=True)
            f.write("\n")
