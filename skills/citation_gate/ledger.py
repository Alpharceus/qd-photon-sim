#!/usr/bin/env python3
# Copyright (c) 2026 Cheng-I Wu. Licensed CC-BY-NC-4.0 (non-commercial
# research use). Original code written for qd-photon-sim on top of the
# vendored `gate.verify_citation` waterfall (see THIRD_PARTY.md); not
# itself derived from an upstream file.
"""Ledger-shaped citation verification.

A "ledger" here is the evidence-ledger YAML schema already used by this
project's other verify_* scripts: a top-level `anchors:` mapping keyed by
anchor id, each entry carrying at least `citation` (a free-text citation
string) and `doi_or_url` (a DOI, `arXiv:<id>`, `OSTI <id>`, `PMC<id>`, a
`https://doi.org/...` URL, or null). `evidence_status` (when present) is
read only for the null-identifier "recorded gap" rule below.

Citation strings in this project's ledgers follow the common physics-
journal style with NO title ("Author, Author2, Journal Vol, Page
(Year)"); a handful carry a quoted title. Title/year parsing is
therefore deliberately best-effort (see `parse_citation`).
"""
from __future__ import annotations

import re
import unicodedata
from typing import Any

from .arxiv_client import ArxivClient
from .crossref_client import CrossrefClient
from .gate import verify_citation
from .semantic_scholar_client import SemanticScholarClient

_YEAR_RE = re.compile(r"\(~?\s*(\d{4})\)")
_TITLE_RE = re.compile("'([^']{6,})'|\"([^\"]{6,})\"")

_ASCII_REPLACEMENTS = {
    chr(0x2009): ' ',  # thin space
    chr(0x2013): '-',  # en dash
    chr(0x2014): '-',  # em dash
    chr(0x2018): "'",  # left single quote
    chr(0x2019): "'",  # right single quote
    chr(0x201C): '"',  # left double quote
    chr(0x201D): '"',  # right double quote
    chr(0x2026): '...',  # ellipsis
    chr(0x00A0): ' ',  # non-breaking space
}


def to_ascii(s: str | None) -> str:
    """ASCII-safe rendering for stdout (this project prints only ASCII).
    A resolver-returned title can carry Unicode punctuation (thin spaces,
    em dashes, curly quotes); fold the common ones to an ASCII equivalent,
    then drop anything else via NFKD + ascii-ignore. The cache file itself
    stays UTF-8/Unicode -- this is a print-time transform only."""
    if not s:
        return ""
    for u, a in _ASCII_REPLACEMENTS.items():
        s = s.replace(u, a)
    return unicodedata.normalize("NFKD", s).encode("ascii", "ignore").decode("ascii")


def default_clients() -> dict[str, Any]:
    """One shared client per resolver, so rate-limit pacing (per-instance
    state) is respected across an entire ledger run rather than reset for
    every anchor."""
    return {
        "crossref": CrossrefClient(),
        "semantic_scholar": SemanticScholarClient(),
        "arxiv": ArxivClient(),
    }


def parse_citation(citation: str) -> dict[str, Any]:
    """Best-effort {"title", "year", "authors_raw"} from a free-text
    citation string. `title` is "" when no quoted title is present (the
    common case here); `year` is None when no `(YYYY)` / `(~YYYY)` marker
    is found."""
    title = ""
    m = _TITLE_RE.search(citation or "")
    if m:
        title = (m.group(1) or m.group(2) or "").strip()
    year = None
    m = _YEAR_RE.search(citation or "")
    if m:
        year = int(m.group(1))
    authors_raw = (citation or "").split("(")[0].strip().rstrip(",")
    return {"title": title, "year": year, "authors_raw": authors_raw}


def parse_identifier(doi_or_url: str | None) -> tuple[str | None, str | None]:
    """(kind, value) from a `doi_or_url` field.

    kind is one of "doi", "arxiv", "osti", "pmc", "unknown", or None
    (when doi_or_url itself is None).
    """
    if doi_or_url is None:
        return None, None
    s = doi_or_url.strip()
    if s.startswith("arXiv:"):
        return "arxiv", s[len("arXiv:"):]
    if s.startswith("OSTI "):
        return "osti", s[len("OSTI "):]
    if s.startswith("PMC"):
        return "pmc", s
    if s.startswith("http") and "doi.org/" in s:
        return "doi", s.split("doi.org/", 1)[1]
    if s.startswith("10."):
        return "doi", s
    return "unknown", s


def verify_anchor(
    anchor_id: str, anchor: dict[str, Any], *, cache=None, offline: bool = False,
    clients: dict[str, Any] | None = None,
) -> tuple[bool, str, str]:
    """Verify one ledger anchor.

    Returns (passed, identifier_display, detail) where `detail` is either
    the matched title (on a resolver match) or the failure/recorded-gap
    reason, suitable for the 4th column of the report line.
    """
    citation = anchor.get("citation") or ""
    doi_or_url = anchor.get("doi_or_url")
    evidence_status = anchor.get("evidence_status")
    parsed = parse_citation(citation)
    kind, value = parse_identifier(doi_or_url)
    identifier_display = doi_or_url if doi_or_url is not None else "null"

    if kind is None:
        if evidence_status == "missing":
            return True, identifier_display, "recorded gap (null identifier, evidence_status=missing)"
        if parsed["title"]:
            # A null identifier with a digest-printed title (the ledger rule
            # forbids inventing a DOI) is resolved by title search only,
            # exactly like the OSTI/PMC case below.
            entry = {"doi": None, "arxiv_id": None, "title": parsed["title"], "year": parsed["year"]}
            outcome = verify_citation(entry, cache=cache, offline=offline, clients=clients)
            if outcome["matched"] and outcome["resolver"] in ("crossref_title", "semantic_scholar"):
                return True, identifier_display, outcome.get("title") or "title match"
            return False, identifier_display, "; ".join(outcome["trace"]) or "no title match"
        return False, identifier_display, (
            f"null identifier, no parseable title, and evidence_status={evidence_status!r} "
            "(only evidence_status=missing auto-passes with no identifier and no title)"
        )

    if kind == "doi":
        entry = {"doi": value, "arxiv_id": None, "title": parsed["title"], "year": parsed["year"]}
        outcome = verify_citation(entry, cache=cache, offline=offline, clients=clients)
        if outcome["matched"]:
            return True, identifier_display, outcome.get("title") or f"doi resolved (year {outcome.get('year')})"
        return False, identifier_display, "; ".join(outcome["trace"])

    if kind == "arxiv":
        entry = {"doi": None, "arxiv_id": value, "title": parsed["title"], "year": parsed["year"]}
        outcome = verify_citation(entry, cache=cache, offline=offline, clients=clients)
        if outcome["matched"]:
            return True, identifier_display, outcome.get("title") or f"arxiv resolved (year {outcome.get('year')})"
        return False, identifier_display, "; ".join(outcome["trace"])

    if kind in ("osti", "pmc"):
        # Neither OSTI nor PMC is one of the four public APIs this gate
        # calls -- these identifiers are resolved by title search only
        # (Crossref, then Semantic Scholar), per spec.
        if not parsed["title"]:
            return False, identifier_display, (
                f"no title parseable from citation string; cannot title-search "
                f"a bare {kind.upper()} identifier"
            )
        entry = {"doi": None, "arxiv_id": None, "title": parsed["title"], "year": parsed["year"]}
        outcome = verify_citation(entry, cache=cache, offline=offline, clients=clients)
        if outcome["matched"] and outcome["resolver"] in ("crossref_title", "semantic_scholar"):
            return True, identifier_display, outcome.get("title") or "title match"
        return False, identifier_display, "; ".join(outcome["trace"]) or "no title match"

    return False, identifier_display, f"unrecognized identifier form: {value!r}"


def verify_ledger(
    anchors: dict[str, Any], *, cache=None, offline: bool = False,
    clients: dict[str, Any] | None = None,
) -> tuple[list[tuple[str, bool, str, str]], int, int]:
    """Verify every anchor in a `{anchors: {...}}`-shaped mapping (already
    unwrapped -- pass the inner mapping). Returns (rows, passed, total)
    where each row is (anchor_id, passed, identifier_display, detail).

    `clients` defaults to one shared client set for the whole call (see
    `default_clients`), so rate-limit pacing is respected across every
    anchor rather than reset per anchor.
    """
    if clients is None:
        clients = default_clients()
    rows = []
    passed = 0
    for anchor_id, anchor in anchors.items():
        ok, identifier_display, detail = verify_anchor(
            anchor_id, anchor, cache=cache, offline=offline, clients=clients)
        rows.append((anchor_id, ok, identifier_display, detail))
        if ok:
            passed += 1
    return rows, passed, len(rows)
