# citation-gate

## Licence

This package vendors and adapts code from
[`academic-research-skills`](https://github.com/Imbad0202/academic-research-skills)
by Cheng-I Wu, at tag **v3.21.2**, licensed under the **Creative Commons
Attribution-NonCommercial 4.0 International License (CC-BY-NC-4.0)**.

Use of this package within `qd-photon-sim` is **non-commercial research
use**, consistent with that licence. Every adapted or vendored file below
carries a copyright/licence header pointing back here; see
`THIRD_PARTY.md` at the repository root for the complete list, the
upstream URL, and the full licence text (also mirrored at
`skills/citation_gate/LICENSE`).

## What was vendored vs. adapted vs. written fresh

- `crossref_client.py`, `verification_cache.py`, and the resolver
  waterfall in `gate.py` (`verify_citation`) are vendored/adapted from
  the upstream files of the same purpose
  (`scripts/crossref_client.py`, `scripts/verification_cache.py`,
  `scripts/verification_gate/__init__.py`), trimmed to what this project
  needs (see the adaptation notes at the top of each file): a JSON,
  repo-local cache instead of a home-directory SQLite database; no
  environment variable required; no cross-model, passport, inquiry-ledger,
  or hook code; a simple DOI -> Crossref-title -> Semantic-Scholar ->
  arXiv waterfall instead of the upstream's four-way parallel vote
  (OpenAlex is not used here).
- `semantic_scholar_client.py` and `arxiv_client.py` are original code
  written for this project, following the same DOI/id-then-title lookup
  pattern and public API endpoints as the upstream clients of the same
  name, but not vendored verbatim -- those two files were not among the
  pieces named for vendoring into this project.
- `text_similarity.py` is a minimal reimplementation of the title
  normalization + similarity-ratio + 0.70 threshold that
  `crossref_client.py` depends on (the upstream `_text_similarity.py`
  module was likewise not one of the vendored pieces; its dotted-acronym
  and CJK-title normalization passes are irrelevant to this project's
  ASCII physics-journal citations and were intentionally left out).
- `ledger.py` and `__main__.py` are original code: the ledger-shaped
  verification loop, citation-string parsing, and CLI are specific to
  this project's evidence-ledger YAML schema.

## No home-directory or environment footprint

- Cache: `verify/data/citation_verification_cache.json` inside this
  repository (JSON, UTF-8, deterministic key order) -- not the upstream
  tool's SQLite database under the caller's home-directory cache folder.
- No environment variable is required to run. `S2_API_KEY`, if already
  present in the environment, is read (never written) to raise the
  Semantic Scholar rate limit.
- No home-directory path (a user profile cache folder, a Claude
  configuration directory, or otherwise) is read or written.
