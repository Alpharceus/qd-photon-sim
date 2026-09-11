# Third-party code

This repository vendors and adapts a small number of files from
**academic-research-skills** by Cheng-I Wu.

- Upstream repository: https://github.com/Imbad0202/academic-research-skills
- Tag vendored from: **v3.21.2**
- Licence: **Creative Commons Attribution-NonCommercial 4.0 International
  (CC-BY-NC-4.0)** -- full text at `skills/citation_gate/LICENSE`
- Use: non-commercial research use, consistent with the licence.

## Vendored / adapted files

| File in this repo | Adapted from (upstream path) | Notes |
|---|---|---|
| `skills/citation_gate/crossref_client.py` | `scripts/crossref_client.py` | Trimmed: no polite-pool email / env var, no dual-path import fallback; DOI-lookup + title-search logic and retry/backoff unchanged. |
| `skills/citation_gate/verification_cache.py` | `scripts/verification_cache.py` | Storage changed from a SQLite database under `~/.cache/ars/` to a single JSON file inside this repository (`verify/data/citation_verification_cache.json`, no environment variable); TTL/stale-advisory machinery dropped (not needed here). |
| `skills/citation_gate/gate.py` (`verify_citation`) | `scripts/verification_gate/__init__.py` | Simplified from the upstream's four-way parallel resolver vote (crossref/openalex/semantic_scholar/arxiv + ref_slug/anchor prose-join) to a linear waterfall: DOI lookup, then Crossref title search, then Semantic Scholar, then arXiv. OpenAlex and the ref_slug/anchor/passport join layer were not carried over. |
| `skills/source-verification/SKILL.md` | `deep-research/agents/source_verification_agent.md` | Reduced to evidence grading, per-claim verification verdicts, and the predatory-venue flag; pipeline phase-boundary boilerplate removed; the evidence hierarchy is mapped onto this project's own `[V]/[DR]/[E]/[A]` tag convention. |

Every file above carries its own copyright/licence header at the top
pointing back to this document.

## Written fresh (not vendored), following the same public API pattern

`skills/citation_gate/semantic_scholar_client.py`,
`skills/citation_gate/arxiv_client.py`, and
`skills/citation_gate/text_similarity.py` were **not** among the pieces
named for vendoring into this project. They are original code written
for `qd-photon-sim`, following the same DOI/id-then-title lookup pattern,
public API endpoints, and (for `text_similarity.py`) the 0.70
title-similarity threshold as the upstream modules of similar purpose,
but are not copies of upstream source. See the header comment in each
file.

`skills/citation_gate/ledger.py` and `skills/citation_gate/__main__.py`
are original code specific to this project's evidence-ledger YAML
schema and are not derived from any upstream file.

## Attribution files copied verbatim

- `skills/citation_gate/LICENSE` -- the upstream CC-BY-NC-4.0 licence text.
