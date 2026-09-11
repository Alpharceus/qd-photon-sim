---
name: source-verification
description: "grade evidence quality and verify per-claim verdicts against this project's [V]/[DR]/[E]/[A] provenance tags"
---

# source-verification

Reduced from the `source_verification_agent` prompt of
`academic-research-skills` (CC-BY-NC-4.0, tag v3.21.2) to what this
project needs: evidence grading, per-claim verification verdicts, and
predatory-venue flags, with the pipeline phase-boundary boilerplate
(Phase 2 write-scope fencing, cross-model/passport machinery) removed
and the evidence hierarchy mapped onto this project's own `[V]/[DR]/[E]/[A]`
provenance-tag convention (`CLAUDE.md`, "Conventions"). See
`THIRD_PARTY.md` at the repository root for attribution.

Use this when grading or spot-checking a literature number before it
enters `fsim_core/`, a card, or an evidence ledger (`verify/data/*.yaml`)
-- not for checking that a citation merely *exists*; that is
`skills/citation-gate` (`skills/citation_gate/`).

## Evidence hierarchy -> project tag mapping

Upstream's 7-level evidence hierarchy (systematic reviews down to expert
opinion) grades *study design*. This project's tags grade *how a number
got from a paper into code*, which is a different axis; the two compose
rather than substitute for one another. Map a source's hierarchy level
into the project tag it can support:

| Upstream level | What it means | Project tag it can support |
|---|---|---|
| I-IV (systematic review, RCT, controlled/cohort study) | Directly measured, high confidence | `[V]` -- verified against the cited paper, IF the number is read from the paper's own text/table/figure (not re-derived or transferred across platforms) |
| V-VI (review of descriptive studies, single descriptive/qualitative study) | Single measurement or reading, lower confidence | `[V]` for a direct read, `[DR]` if the project's own code derives a further quantity from it |
| VII (expert opinion / committee report) | No new measurement | `[E]` estimate/class range at best, never `[V]` |
| No source located, or an off-platform/off-condition transfer of an otherwise `[V]` number onto a *different* card/platform | -- | `[A]` assumption or `[E]`, even when the underlying measurement is itself tagged `[V]` in its own ledger entry -- never conflate a source's own grade with the grade of transferring it elsewhere (see `verify/data/nitride_cavity_anchors.yaml`'s header comment on this exact point) |

`evidence_status` in this project's ledgers (`full_text`, `abstract_only`,
`figure_reading`, `estimate`, `missing`) is the closest local analogue of
upstream's Tier 0-2 evidence-strength ladder: `full_text` supports `[V]`,
`abstract_only`/`figure_reading` still support `[V]` but with a narrower
excerpt (record what was actually read in `location`), `estimate` caps at
`[E]`, and `missing` records a gap rather than inventing a number (never
silently promote a `missing` entry to a tag by guessing).

## Per-claim verification verdicts

For a specific numeric claim (a card field, a `fsim_core/` constant, a
ledger anchor's `value`), verify against the cited paper and record one
of:

- **Confirmed** -- the paper states this number (or one the project's own
  documented derivation reaches from what the paper states); tag stays
  as recorded.
- **Transcription error** -- the paper states a different number; flag
  for correction, do not silently fix the ledger (this project's own
  policy: report, don't auto-edit an evidence ledger).
- **Unsupported transfer** -- the number is real and correctly read, but
  is being used somewhere the paper did not measure it (different
  platform, orientation, temperature, excitation); the consuming tag
  must be `[A]`/`[E]`, not `[V]`, per the mapping table above.
- **Existence unverified** -- the citation itself could not be confirmed
  to exist; hand off to `skills/citation-gate` rather than grading its
  content.

## Predatory-venue flag

Carried over from upstream essentially unchanged, since it is a pure
existence/reputation check independent of the pipeline machinery: flag a
source if the venue is not indexed in any recognized index (Scopus / Web
of Science / DOAJ), the publisher is not a member of COPE or DOAJ, or the
citation shows other classic markers (no identifiable editorial board,
article-processing charges far below market rate, aggressive
solicitation). A flagged venue does not auto-fail the claim; record it
alongside the verdict and let a human weigh it -- consistent with
upstream's "red flags, not censorship" principle.

## What was intentionally dropped

The upstream prompt's Phase 2 write-scope fence, its `PreToolUse` hook
enforcement note, the multi-phase pipeline vocabulary
(`phase{M}_*/` directories, `bibliography_agent` hand-off), and the
Tier-0-2 WebSearch/DOI-resolution procedure (superseded here by
`skills/citation-gate`, which does the same job through actual API calls
rather than manual WebSearch spot-checks) are all out of scope for this
project and were not carried over.
