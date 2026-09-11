#!/usr/bin/env python3
# Copyright (c) 2026 Cheng-I Wu. Licensed CC-BY-NC-4.0 (non-commercial
# research use). Vendored/adapted for qd-photon-sim from
# academic-research-skills (v3.21.2) scripts/crossref_client.py. See
# THIRD_PARTY.md for the upstream URL, tag, and licence.
#
# Adaptation notes (from the original docstring/behavior):
#   - `_text_similarity` import replaced with this package's
#     `text_similarity` (a minimal reimplementation, see that module).
#   - The dual-path (sibling vs. `scripts.`) import fallback was dropped;
#     this package uses a single relative import.
#   - Polite-pool email support (CROSSREF_POLITE_EMAIL) was dropped per
#     spec: no environment variable is required, and the polite-pool
#     email is left empty (anonymous rate only).
#   - Otherwise the DOI-lookup-with-title-check / title-search logic and
#     retry/backoff behavior are unchanged.
"""Minimal Crossref API client wrapper.

DOI-first with title cross-check (DOI_MISMATCH pattern), title-similarity
fallback, 429 -> 2s backoff x 3 retries, 404/5xx -> miss vs. skip.

Crossref-specific: DOI endpoint is /works/{doi} (no doi: prefix); title
search is /works?query.title=...&rows=5; response shape is nested under
`message`; title is a list (multi-language variants).
"""
from __future__ import annotations

import http.client
import json
import time
import urllib.error
import urllib.parse
import urllib.request
from typing import Any, Mapping

from .text_similarity import (
    _BACKOFF_SECONDS,
    _MAX_RETRIES,
    _TITLE_SIMILARITY_THRESHOLD,
    similarity as _similarity,
)

_API_BASE = "https://api.crossref.org"
_API_HOST = "api.crossref.org"

# Crossref polite pool: 10 req/s with mailto, ~5 req/s anonymous. No
# polite-pool email is configured here (spec: left empty), so this client
# always uses the anonymous pacing.
_ANONYMOUS_MIN_INTERVAL = 0.2


def _require_api_url(url: str) -> None:
    parsed = urllib.parse.urlsplit(url)
    if parsed.scheme != "https" or parsed.netloc != _API_HOST:
        raise CrossrefUnavailable(f"Refusing non-Crossref URL: {url}")


def _extract_title(message_or_item: Mapping[str, Any]) -> str:
    """Crossref returns `title` as a list of language variants. Take first or empty."""
    titles = message_or_item.get("title") or []
    return titles[0] if titles else ""


def _extract_year(item: Mapping[str, Any]) -> int | None:
    """Crossref year lives in `issued.date-parts[0][0]` (or `published-print` /
    `published-online`). Prefer `issued` as canonical; fall through to
    alternatives."""
    for key in ("issued", "published-print", "published-online"):
        val = item.get(key)
        if not isinstance(val, dict):
            continue
        date_parts = val.get("date-parts")
        if date_parts and date_parts[0]:
            return date_parts[0][0]
    return None


class CrossrefUnavailable(Exception):
    """Crossref API degraded -- caller must treat this citation's Crossref
    step as inconclusive, not as a hard failure of the whole run."""


class CrossrefClient:
    """Production lookup-by-(doi-with-cross-check-then-title) client for
    Crossref. Concurrency note: rate-limit pacing is per-instance."""

    def __init__(self) -> None:
        self._min_interval = _ANONYMOUS_MIN_INTERVAL
        self._last_request_at: float | None = None
        self._user_agent = "qd-photon-sim-citation-gate/1.0"

    def _throttle(self) -> None:
        if self._last_request_at is None:
            return
        elapsed = time.monotonic() - self._last_request_at
        if elapsed < self._min_interval:
            time.sleep(self._min_interval - elapsed)

    def _get(self, path: str, query: Mapping[str, str]) -> dict[str, Any]:
        url = f"{_API_BASE}{path}"
        if query:
            url += "?" + urllib.parse.urlencode(query)
        _require_api_url(url)
        req = urllib.request.Request(url, headers={"User-Agent": self._user_agent})

        self._throttle()
        self._last_request_at = time.monotonic()

        for attempt in range(_MAX_RETRIES + 1):
            try:
                with urllib.request.urlopen(req, timeout=30) as resp:  # nosec B310
                    try:
                        body = resp.read()
                        return json.loads(body.decode("utf-8"))
                    except (
                        OSError,
                        http.client.HTTPException,
                        UnicodeDecodeError,
                        json.JSONDecodeError,
                    ) as e:
                        raise CrossrefUnavailable(
                            f"Crossref response read/parse failed: {e}"
                        ) from e
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    return {}
                if e.code == 429 and attempt < _MAX_RETRIES:
                    time.sleep(_BACKOFF_SECONDS)
                    self._last_request_at = time.monotonic()
                    continue
                raise CrossrefUnavailable(f"Crossref HTTP {e.code}: {e.reason}") from e
            except (urllib.error.URLError, TimeoutError) as e:
                raise CrossrefUnavailable(f"Crossref network error: {e}") from e

        raise CrossrefUnavailable("Crossref rate limit exhausted after retries")

    def doi_lookup_with_title_check(
        self, doi: str, expected_title: str,
    ) -> dict[str, Any] | None:
        """DOI lookup with title cross-check (skipped when expected_title is
        empty -- callers with no parseable title accept the bare DOI hit).

        Returns the `message` dict if the DOI resolves AND (no expected
        title, or the title cross-check passes); None on 404 (miss) or
        DOI_MISMATCH.
        """
        data = self._get(f"/works/{urllib.parse.quote(doi, safe='')}", {})
        if not data:  # 404 -> empty dict from _get
            return None
        message = data.get("message", {})
        if not expected_title:
            return message
        title = _extract_title(message)
        if _similarity(title, expected_title) >= _TITLE_SIMILARITY_THRESHOLD:
            return message
        return None  # DOI_MISMATCH

    def title_search(
        self, title: str, year: int | None = None,
    ) -> dict[str, Any] | None:
        """Title search: a candidate matches iff it clears the 0.70
        similarity ratio (best-scoring candidate wins; a year match adds a
        small tiebreak bonus)."""
        if not title:
            return None
        data = self._get("/works", {"query.title": title, "rows": "5"})
        candidates = data.get("message", {}).get("items", [])
        scored = []
        for cand in candidates:
            cand_title = _extract_title(cand)
            sim = _similarity(cand_title, title)
            if sim < _TITLE_SIMILARITY_THRESHOLD:
                continue
            year_match = year is not None and _extract_year(cand) == year
            score = sim + (0.05 if year_match else 0.0)
            scored.append((cand, score))
        if not scored:
            return None
        scored.sort(key=lambda cand_score: (-cand_score[1],))
        return scored[0][0]

    @staticmethod
    def extract_title(message: Mapping[str, Any]) -> str:
        return _extract_title(message)

    @staticmethod
    def extract_year(message: Mapping[str, Any]) -> int | None:
        return _extract_year(message)
