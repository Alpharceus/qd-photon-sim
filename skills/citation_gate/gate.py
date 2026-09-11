#!/usr/bin/env python3
# Copyright (c) 2026 Cheng-I Wu. Licensed CC-BY-NC-4.0 (non-commercial
# research use). Vendored/adapted for qd-photon-sim from
# academic-research-skills (v3.21.2) scripts/verification_gate/__init__.py
# (the `verify_citation` resolver waterfall). See THIRD_PARTY.md for the
# upstream URL, tag, and licence.
#
# Adaptation notes (from the original docstring/behavior):
#   - The upstream `verify_citation` runs all four resolvers (crossref,
#     openalex, semantic_scholar, arxiv) and reduces their outcomes via a
#     3-class lookup_verified vote. This project needs only the simpler
#     waterfall the spec asks for: DOI lookup, then Crossref title
#     search, then Semantic Scholar, then arXiv -- stopping at the first
#     match. OpenAlex is not used (not one of the four public APIs the
#     upstream code exercises for this project; see THIRD_PARTY.md).
#   - No `ref_slug` / `anchor` / passport / prose-join machinery: this
#     project's citations are keyed directly by ledger anchor id, not by
#     a writer-prose marker, so that whole join layer (upstream #332) is
#     out of scope and was not carried over.
#   - `entry` here is a plain dict of {doi, arxiv_id, title, year}; no
#     `citation_key` / `obtained_via` / corpus-schema fields.
"""verify_citation -- citation existence verification for one entry.

Runs the waterfall: DOI lookup (Crossref) -> Crossref title search ->
Semantic Scholar (DOI-then-title) -> arXiv (id-then-title, only when an
arxiv_id is present). Stops at the first resolver step that returns a
result which matches `entry`'s expected title (when the entry has a
parsed title) and expected year (when the entry has a parsed year and the
resolver returned one). Every step's result -- hit or miss -- is
persisted through `cache` so a second run can reproduce the same verdict
with `offline=True` and no resolver call at all.
"""
from __future__ import annotations

from typing import Any

from .arxiv_client import ArxivClient, ArxivUnavailable
from .crossref_client import CrossrefClient, CrossrefUnavailable
from .semantic_scholar_client import SemanticScholarClient, SemanticScholarUnavailable
from .text_similarity import _TITLE_SIMILARITY_THRESHOLD, similarity as _similarity
from .verification_cache import VerificationCache

_NOT_IN_CACHE = {"found": False, "note": "not in cache"}


def _acceptable(
    resolved_title: str | None,
    expected_title: str,
    resolved_year: int | None,
    expected_year: int | None,
) -> bool:
    """A resolver hit is acceptable iff every check the entry actually
    supplies agrees: title similarity >= 0.70 when `expected_title` is
    non-empty, AND year equality when both years are known. A citation
    string with no parseable title (the common case for this project's
    "Author, Journal Vol, Page (Year)" style) only needs the year to
    agree; one with no parseable year only needs the title to agree.
    An entry that supplies neither is accepted on existence alone."""
    if expected_title:
        if not resolved_title:
            return False
        if _similarity(resolved_title, expected_title) < _TITLE_SIMILARITY_THRESHOLD:
            return False
    if expected_year is not None and resolved_year is not None:
        # +/-1 year tolerance: a Crossref "issued" date is often the
        # online-first date, one year ahead of the print-issue year a
        # citation names (see e.g. taylor2010-cavity-q, issued 2009-12,
        # published-print 2010-03). Matches the same tolerance documented
        # in the vendored source_verification_agent.md's Tier 0 S2 rule
        # ("year matches (or within +/-1 year)").
        if abs(resolved_year - expected_year) > 1:
            return False
    return True


def _cached_call(cache, resolver_name, query_form, fetch_fn, offline):
    """Cache-through wrapper shared by every resolver step.

    Returns (response_dict, note) where `response_dict` always has a
    "found" key. `note`, when set, is a short human-readable reason
    string for the trace (e.g. "not in cache", or an *Unavailable
    message) -- None on an ordinary cache/network hit-or-miss.
    """
    cached = cache.get(resolver_name, query_form) if cache is not None else None
    if cached is not None:
        return cached, None
    if offline:
        return dict(_NOT_IN_CACHE), "not in cache"
    try:
        result = fetch_fn()
    except (CrossrefUnavailable, SemanticScholarUnavailable, ArxivUnavailable) as e:
        result = {"found": False, "note": str(e)}
    if cache is not None:
        cache.put(resolver_name, query_form, result)
    return result, (None if result.get("found") else result.get("note"))


def verify_citation(
    entry: dict[str, Any],
    *,
    cache: VerificationCache | None = None,
    offline: bool = False,
    clients: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """Verify one citation's existence across the DOI/Crossref/S2/arXiv
    waterfall.

    `entry`: {"doi": str|None, "arxiv_id": str|None, "title": str,
    "year": int|None}. `clients` optionally injects
    {"crossref", "semantic_scholar", "arxiv"} client instances (for
    tests); defaults to one live instance of each.

    Returns {"matched": bool, "resolver": str|None, "title": str|None,
    "year": int|None, "trace": [str, ...]}. `trace` records one short
    reason per step that did not match, in the order tried; `resolver`
    names the step that matched (None when unmatched).
    """
    clients = clients or {}
    crossref = clients.get("crossref") or CrossrefClient()
    s2 = clients.get("semantic_scholar") or SemanticScholarClient()
    arxiv = clients.get("arxiv") or ArxivClient()

    doi = entry.get("doi")
    arxiv_id = entry.get("arxiv_id")
    title = entry.get("title") or ""
    year = entry.get("year")
    trace: list[str] = []

    def matched(resolver, resolved):
        return {
            "matched": True,
            "resolver": resolver,
            "title": resolved.get("title"),
            "year": resolved.get("year"),
            "trace": trace,
        }

    # Step 1: DOI lookup (Crossref).
    if doi:
        result, note = _cached_call(
            cache, "crossref_doi", doi,
            lambda: _project_crossref(crossref.doi_lookup_with_title_check(doi, title)),
            offline,
        )
        if result.get("found") and _acceptable(result.get("title"), title, result.get("year"), year):
            return matched("crossref_doi", result)
        trace.append(f"crossref DOI lookup: {note or 'no match'}")

    # Step 2: Crossref title search.
    if title:
        result, note = _cached_call(
            cache, "crossref_title", title,
            lambda: _project_crossref(crossref.title_search(title, year)),
            offline,
        )
        if result.get("found") and _acceptable(result.get("title"), title, result.get("year"), year):
            return matched("crossref_title", result)
        trace.append(f"crossref title search: {note or 'no match'}")

    # Step 3: Semantic Scholar (DOI-then-title internally).
    if doi or title:
        query_form = doi or title
        result, note = _cached_call(
            cache, "semantic_scholar", query_form,
            lambda: _project_s2(s2.lookup(doi, title, year)),
            offline,
        )
        if result.get("found") and _acceptable(result.get("title"), title, result.get("year"), year):
            return matched("semantic_scholar", result)
        trace.append(f"semantic scholar: {note or 'no match'}")

    # Step 4: arXiv (only applicable when the entry carries an arXiv id).
    if arxiv_id:
        result, note = _cached_call(
            cache, "arxiv_id", arxiv_id,
            lambda: _project_arxiv(arxiv.id_lookup(arxiv_id, title)),
            offline,
        )
        if result.get("found") and _acceptable(result.get("title"), title, result.get("year"), year):
            return matched("arxiv_id", result)
        trace.append(f"arxiv id lookup: {note or 'no match'}")
        if title:
            result, note = _cached_call(
                cache, "arxiv_title", title,
                lambda: _project_arxiv(arxiv.title_search(title, year)),
                offline,
            )
            if result.get("found") and _acceptable(result.get("title"), title, result.get("year"), year):
                return matched("arxiv_title", result)
            trace.append(f"arxiv title search: {note or 'no match'}")

    if not trace:
        trace.append("no identifier or title to verify")
    return {"matched": False, "resolver": None, "title": None, "year": None, "trace": trace}


def _project_crossref(message):
    if not message:
        return {"found": False}
    return {
        "found": True,
        "title": CrossrefClient.extract_title(message),
        "year": CrossrefClient.extract_year(message),
    }


def _project_s2(result):
    if not result:
        return {"found": False}
    return {"found": True, "title": result.get("title") or "", "year": result.get("year")}


def _project_arxiv(result):
    if not result:
        return {"found": False}
    return {"found": True, "title": result.get("title") or "", "year": result.get("year")}
