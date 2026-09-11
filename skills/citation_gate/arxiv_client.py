#!/usr/bin/env python3
# Copyright (c) 2026 Cheng-I Wu. Licensed CC-BY-NC-4.0 (non-commercial
# research use). Adapted for qd-photon-sim, following the ID-then-title
# lookup pattern of academic-research-skills (v3.21.2)
# scripts/arxiv_client.py -- rewritten to this project's stdlib-only style
# rather than vendored verbatim, since the source file was not one of the
# pieces named for vendoring (see THIRD_PARTY.md); the endpoint, Atom
# parsing, and 429 pacing mirror the upstream client's documented
# behavior (arXiv's own Terms of Use ask for >=3s between requests).
"""Minimal arXiv API client wrapper.

arXiv-ID-first with title cross-check, title-similarity fallback. The
query API returns Atom 1.0 XML (not JSON); the exact-key endpoint is
`?id_list={id}`, the title fallback is `?search_query=ti:"{title}"`.
"""
from __future__ import annotations

import time
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from typing import Any

from .text_similarity import _MAX_RETRIES, _TITLE_SIMILARITY_THRESHOLD, similarity as _similarity

_API_BASE = "http://export.arxiv.org/api/query"
_ATOM_NS = "{http://www.w3.org/2005/Atom}"

# arXiv API Terms of Use ask callers to pace requests ~3s apart.
_ARXIV_MIN_INTERVAL = 3.0


class ArxivUnavailable(Exception):
    """arXiv API degraded -- caller must treat this citation's arXiv step
    as inconclusive, not as a hard failure of the whole run."""


def _extract_title(entry: ET.Element) -> str:
    node = entry.find(f"{_ATOM_NS}title")
    if node is None or node.text is None:
        return ""
    return " ".join(node.text.split())


def _extract_year(entry: ET.Element) -> int | None:
    node = entry.find(f"{_ATOM_NS}published")
    if node is None or not node.text:
        return None
    head = node.text[:4]
    return int(head) if head.isdigit() else None


class ArxivClient:
    """Lookup-by-(arxiv-id-with-cross-check-then-title) client."""

    def __init__(self) -> None:
        self._min_interval = _ARXIV_MIN_INTERVAL
        self._last_request_at: float | None = None
        self._user_agent = "qd-photon-sim-citation-gate/1.0"

    def _throttle(self) -> None:
        if self._last_request_at is None:
            return
        elapsed = time.monotonic() - self._last_request_at
        if elapsed < self._min_interval:
            time.sleep(self._min_interval - elapsed)

    def _get(self, query: dict[str, str]) -> list[ET.Element]:
        url = _API_BASE + "?" + urllib.parse.urlencode(query)
        req = urllib.request.Request(url, headers={"User-Agent": self._user_agent})

        self._throttle()
        self._last_request_at = time.monotonic()

        for attempt in range(_MAX_RETRIES + 1):
            try:
                with urllib.request.urlopen(req, timeout=30) as resp:  # nosec B310
                    try:
                        root = ET.fromstring(resp.read())
                    except (OSError, ET.ParseError) as e:
                        raise ArxivUnavailable(
                            f"arXiv response read/parse failed: {e}"
                        ) from e
                    if root.tag != f"{_ATOM_NS}feed":
                        raise ArxivUnavailable(
                            f"arXiv returned a non-Atom body (root tag {root.tag!r})"
                        )
                    return root.findall(f"{_ATOM_NS}entry")
            except urllib.error.HTTPError as e:
                if e.code == 429 and attempt < _MAX_RETRIES:
                    time.sleep(_ARXIV_MIN_INTERVAL)
                    self._last_request_at = time.monotonic()
                    continue
                raise ArxivUnavailable(f"arXiv HTTP {e.code}: {e.reason}") from e
            except (urllib.error.URLError, TimeoutError) as e:
                raise ArxivUnavailable(f"arXiv network error: {e}") from e
        raise ArxivUnavailable("arXiv rate limit exhausted after retries")

    def id_lookup(self, arxiv_id: str, expected_title: str = "") -> dict[str, Any] | None:
        """arXiv ID lookup with an optional title cross-check (skipped when
        expected_title is empty). Returns {"title", "year"} on a hit that
        clears the check; None on empty feed (miss) or ID_MISMATCH."""
        entries = self._get({"id_list": arxiv_id})
        if not entries:
            return None
        entry = entries[0]
        title = _extract_title(entry)
        if expected_title and _similarity(title, expected_title) < _TITLE_SIMILARITY_THRESHOLD:
            return None  # ID_MISMATCH
        return {"title": title, "year": _extract_year(entry)}

    def title_search(self, title: str, year: int | None = None) -> dict[str, Any] | None:
        if not title:
            return None
        entries = self._get({"search_query": f'ti:"{title}"', "max_results": "5"})
        scored = []
        for cand in entries:
            cand_title = _extract_title(cand)
            sim = _similarity(cand_title, title)
            if sim < _TITLE_SIMILARITY_THRESHOLD:
                continue
            year_match = year is not None and _extract_year(cand) == year
            scored.append((cand, sim + (0.05 if year_match else 0.0)))
        if not scored:
            return None
        scored.sort(key=lambda cand_score: -cand_score[1])
        best = scored[0][0]
        return {"title": _extract_title(best), "year": _extract_year(best)}
