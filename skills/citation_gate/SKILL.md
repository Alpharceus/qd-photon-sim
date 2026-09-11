---
name: citation-gate
description: "verify that every citation in a ledger or bibliography exists, via DOI/Crossref/Semantic Scholar/arXiv, with a repo-local cache"
---

# citation-gate

Verifies that citations actually resolve to a real, existing paper,
rather than being fabricated or mistyped, by querying Crossref, Semantic
Scholar, and arXiv's public APIs (DOI-first, falling through to a title
search) and cross-checking the resolved title/year against what the
citation string claims.

## How to call it

```
python -m skills.citation_gate <ledger.yaml> [--offline] [--cache PATH]
```

`<ledger.yaml>` must have a top-level `anchors:` mapping keyed by an
arbitrary id, each entry carrying at least:

- `citation`: a free-text citation string (title, author, journal,
  year -- whatever form the ledger already uses; title and year are
  parsed best-effort, so a title-less "Author, Journal Vol, Page (Year)"
  string is fine).
- `doi_or_url`: a DOI (`10.xxxx/...`), `arXiv:<id>`, `OSTI <id>`,
  `PMC<id>`, a `https://doi.org/...` URL, or `null`.
- `evidence_status` (optional): only `missing` is read, to auto-pass a
  citation that is deliberately recorded with no identifier (see
  Outcomes below).

Prints one line per anchor:

```
<anchor_id> | <identifier> | PASS|FAIL | <matched title or reason>
```

followed by `N/N citation checks passed`, and exits 0 iff every anchor
passed.

## Outcomes

- **DOI or arXiv identifier present**: resolved via the DOI-lookup /
  arXiv-id-lookup step (falling through to a Crossref title search, then
  Semantic Scholar, then an arXiv title search on a miss). Passes if a
  resolver hit's title (when the citation string has a parseable title)
  and year (when the citation string has a parseable year) agree with
  the citation string.
- **`OSTI <id>` / `PMC<id>` identifier**: neither OSTI nor PMC is one of
  the public APIs this gate calls, so these are resolved by title search
  only (Crossref, then Semantic Scholar). Requires a parseable title in
  the citation string; passes on a title match.
- **`doi_or_url: null` with `evidence_status: missing`**: passes as a
  "recorded gap" -- the ledger has already documented that no identifier
  is available for this claim, so there is nothing to verify.
- **Anything else** (a null identifier with any other evidence_status,
  or an unrecognized identifier form): fails with a reason. This project
  never edits a ledger to make a citation pass; a failure here is
  reported, not silently fixed.

## Cache

Every resolver call (hit or miss) is cached at
`verify/data/citation_verification_cache.json` (repo-relative; override
with `--cache`), keyed by `(resolver_step, query_form)`, so re-running
against an unchanged ledger costs no new network calls. `--offline`
serves only from that cache -- a cache miss is a failure ("not in
cache"), and no network call is attempted.

## Scope and non-goals

This is the existence-check layer only: does the cited paper exist, and
does its title/year match the citation. It does not grade evidence
quality, detect predatory venues, or verify the factual claim a citation
is attached to -- see `skills/source-verification/SKILL.md` for that.

No environment variable is required. `S2_API_KEY`, if already set in the
environment, raises the Semantic Scholar rate limit; it is read only,
never written. Nothing is read or written outside this repository.

## Provenance

Adapted from `academic-research-skills` (CC-BY-NC-4.0, tag v3.21.2),
non-commercial research use. See `README.md` in this directory and
`THIRD_PARTY.md` at the repository root for the full attribution.
