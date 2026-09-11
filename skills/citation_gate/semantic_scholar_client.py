#!/usr/bin/env python3
# Copyright (c) 2026 Cheng-I Wu. Licensed CC-BY-NC-4.0 (non-commercial
# research use). Adapted for qd-photon-sim, following the DOI-then-title
# lookup pattern of academic-research-skills (v3.21.2)
# scripts/semantic_scholar_client.py -- rewritten to this project's
# stdlib-only style rather than vendored verbatim, since the source file
# was not one of the pieces named for vendoring (see THIRD_PARTY.md); the
# request shape, throttle, and 429/5xx handling mirror the upstream
# client's documented behavior.
"""Minimal Semantic Scholar API client wrapper.

DOI-first, title-similarity fallback. `S2_API_KEY` is read from the
environment ONLY if already set (raises the request-rate ceiling from
1 req/s to 10 req/s); it is never required and never written.
"""
from __future__ import annotations

import json
import os
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

_API_BASE = "https://api.semanticscholar.org/graph/v1"
_API_HOST = "api.semanticscholar.org"
_API_KEY_ENV = "S2_API_KEY"
_FIELDS = "title,year,externalIds"

_UNAUTHENTICATED_MIN_INTERVAL = 1.0
_AUTHENTICATED_MIN_INTERVAL = 0.1


def _require_api_url(url: str) -> None:
    parsed = urllib.parse.urlsplit(url)
    if parsed.scheme != "https" or parsed.netloc != _API_HOST:
        raise SemanticScholarUnavailable(f"Refusing non-S2 URL: {url}")


class SemanticScholarUnavailable(Exception):
    """S2 API degraded -- caller must treat this citation's S2 step as
    inconclusive, not as a hard failure of the whole run."""


class SemanticScholarClient:
    """Lookup-by-(doi-then-title) client for Semantic Scholar."""

    def __init__(self) -> None:
        api_key = os.environ.get(_API_KEY_ENV)  # optional, read-only
        self._api_key = api_key
        self._min_interval = (
            _AUTHENTICATED_MIN_INTERVAL if api_key else _UNAUTHENTICATED_MIN_INTERVAL
        )
        self._last_request_at: float | None = None

    def _throttle(self) -> None:
        if self._last_request_at is None:
            return
        elapsed = time.monotonic() - self._last_request_at
        if elapsed < self._min_interval:
            time.sleep(self._min_interval - elapsed)

    def _request(self, path: str) -> dict[str, Any]:
        url = f"{_API_BASE}{path}"
        _require_api_url(url)
        headers = {"User-Agent": "qd-photon-sim-citation-gate/1.0"}
        if self._api_key:
            headers["x-api-key"] = self._api_key
        req = urllib.request.Request(url, headers=headers)

        self._throttle()
        self._last_request_at = time.monotonic()

        for attempt in range(_MAX_RETRIES + 1):
            try:
                with urllib.request.urlopen(req, timeout=30) as resp:  # nosec B310
                    return json.loads(resp.read().decode("utf-8"))
            except urllib.error.HTTPError as e:
                if e.code == 404:
                    return {}
                if e.code == 429 and attempt < _MAX_RETRIES:
                    time.sleep(_BACKOFF_SECONDS)
                    self._last_request_at = time.monotonic()
                    continue
                raise SemanticScholarUnavailable(f"S2 API HTTP {e.code}") from e
            except (urllib.error.URLError, TimeoutError) as e:
                raise SemanticScholarUnavailable(f"S2 API network error: {e}") from e
        raise SemanticScholarUnavailable(f"S2 API exhausted {_MAX_RETRIES} retries")

    def lookup_by_doi(self, doi: str, expected_title: str = "") -> dict[str, Any] | None:
        data = self._request(f"/paper/DOI:{urllib.parse.quote(doi, safe='')}?fields={_FIELDS}")
        if not data or not data.get("paperId"):
            return None
        title = data.get("title") or ""
        if expected_title and _similarity(expected_title, title) < _TITLE_SIMILARITY_THRESHOLD:
            return None  # DOI_MISMATCH
        return {"title": title, "year": data.get("year")}

    def lookup_by_title(self, title: str, year: int | None = None) -> dict[str, Any] | None:
        if not title:
            return None
        query = urllib.parse.urlencode({"query": title, "fields": _FIELDS, "limit": "5"})
        data = self._request(f"/paper/search?{query}")
        candidates = data.get("data") or []
        scored = []
        for cand in candidates:
            cand_title = cand.get("title") or ""
            sim = _similarity(cand_title, title)
            if sim < _TITLE_SIMILARITY_THRESHOLD:
                continue
            year_match = year is not None and cand.get("year") == year
            scored.append((cand, sim + (0.05 if year_match else 0.0)))
        if not scored:
            return None
        scored.sort(key=lambda cand_score: -cand_score[1])
        best = scored[0][0]
        return {"title": best.get("title") or "", "year": best.get("year")}

    def lookup(
        self, doi: str | None, title: str, year: int | None = None,
    ) -> dict[str, Any] | None:
        """DOI-first, falling through to title search on miss/mismatch."""
        if doi:
            result = self.lookup_by_doi(doi, title)
            if result is not None:
                return result
        return self.lookup_by_title(title, year)
