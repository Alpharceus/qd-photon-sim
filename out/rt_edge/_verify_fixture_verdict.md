# RT edge-emitter acceptance sweep verdict

Generated 2026-09-06T22:34:20.614924+00:00; contract: docs/rt_edge_contract.md.

```
VERDICT: FAIL g2_min=0.3 g2_median_eligible=0.3 diag_g2_min=0.3 diag_g2_median_diagnostic=0.3 flux_max=nan flux_margin=nan flux_shortfall_deprecated=nan median_pass=true coverage_over_eligible=1 eligible_fraction=1 eligible=1/1 flux_floor_excluded=0 evidence=incomplete conditional=true headline_coverage=1/1 cw_raw_coverage=1/1 gamma300_pass_max=6 gamma300_threshold=>=6 T_pass_min=none headline_by_T=230:0/0,250:0/0,273:0/0,300:0/0
```
```
CARD: edge-inp-gaasp-design role=primary eligible rows: 1/1 g2_pulsed_min=0.3 g2_pulsed_median=0.3 diag_g2_min=0.3 diag_g2_median=0.3 diag_g2_cw0_min=0.3 diag_g2_cw0_median=0.3 diag_g2_cw0_raw_min=0.3 diag_g2_cw0_raw_median=0.3 g2_cw_raw_min=0.3 g2_cw_raw_median=0.3 eligible=1/1 flux_floor_excluded=0 favorable_rows=1
```

## Literature ceiling
In the reviewed literature set of this repository (six papers, ../_goal/paper_digests.md) and the anchors ledger, no electrically driven III-V single-dot g2(0) at 300 K is reported; the best electrical result in that set is Reischle et al., Optics Express 16, 12771 (2008) at 80 K (QD C): g2(0) = 0.43 raw (still IRF-broadened), 0.03 after background correction with the ~500 ps detector IRF also deconvolved. The best reported 300 K single-dot values are optically pumped: g2(0) ~ 0.5-0.57 (Laferriere et al. 2023, InAsP/InP nanowire dot, g2(0) = 0.57 at 300 K). This sweep's pooled pulsed g2_min=0.3, g2_median(eligible)=0.3 is reported against that electrical-vs-optical literature picture, not as a claim of an existing electrical 300 K result to exceed.

**Like-for-like comparison (council review round 5, item 5).** The paragraph above juxtaposed this sweep's IRF-FREE intrinsic g2_min against Reischle's RAW dip (0.43 [V], still IRF-broadened; ledger anchor `reischle08-g2-80k`) -- not a like-for-like convention match. The correct like-for-like anchor is Reischle's IRF-DECONVOLVED-but-background-included value, g2_b(0) = 0.25 [V] +/- 0.05 (QD C, 80 K; ledger anchor `reischle08-g2-80k-deconvolved`, ../_goal/paper_digests.md line ~36), since this sweep's g2_op is likewise background-included (via drive.b_res/rho) and never IRF-convolved for the pulsed metric.
On matching conventions, the 80 K Reischle device is about 1.2x better (lower g2(0)) than this sweep's best 300 K corner (g2_min=0.3 vs 0.25 +/- 0.05).

## Conditional on the 300 K linewidth
The 300 K single-dot linewidth gamma300 is not a platform ceiling -- it is an unmeasured [E]-class quantity (see Assumptions below) whose actual value decides whether the headline gate passes. `gamma300_pass_max` (council review round 5, item 2) is the HIGHEST SAMPLED linewidth, among this determination's own 8-point grid (6, 6.5, 7, 8, 10, 12, 16, 20 meV -- finer than the sweep.csv grid's own 2-point dot.gamma300 endpoints, and never itself written to sweep.csv), at which any eligible row (any lever/delta_xx combination) still clears pulsed intrinsic g2(0) < 0.5; `gamma300_threshold` brackets the true (continuous, unsampled) threshold between that value and the next sampled value above it that already fails -- the true threshold lies somewhere inside the bracket, at the delta_xx shown for that card's gamma300_pass_max row.

**When nothing passes (council review round 6, item 2).** If NO sampled gamma300 clears the headline gate for a card, `gamma300_threshold` no longer prints a numeric bracket like `<6` -- that would misleadingly imply the card passes somewhere below 6 meV, when in fact it never passes anywhere in this grid. It instead prints the REASON class: `none(flux)` if every sampled row is below the collected-flux eligibility floor (headline_pass can never be True there, no matter how favourable g2 would be), or `none(g2)` if at least one sampled row IS flux-eligible but pulsed intrinsic g2(0) >= 0.5 at every sampled gamma300 (a genuine physics shortfall, not an eligibility one); `none(flux,g2)` reports both blocking conditions when neither a flux-eligible row nor a diagnostic g2 < 0.5 row exists.

`gamma300_pass_max` (pooled, both cards): 6 meV; `gamma300_threshold` (pooled): >=6 meV

| card | gamma300_pass_max (meV) | gamma300_threshold bracket (meV) | delta_xx | lever values | assumptions |
|---|---:|---:|---:|---|---|
| edge-inp-gaasp-design | 6 | >=6 | 4 meV | NA=0.75, R_back=0.95, L_um=250 | dot.gamma300; emission.lambda_nm |

Restated (council review round 5, item 2): for each card, `gamma300_pass_max` is the highest SAMPLED linewidth at which the card passes; the true threshold lies between the bracket's two values (at the delta_xx shown above for that row).

### Verified 6.5 meV anchor (Chatzarakis et al., Phys. Rev. Applied 20, 034011, 2023) -- an InAs/GaAs (211)B single-dot linewidth measurement, used here as a distinct-material class proxy [E] for the InP-dot family, not a direct InP/GaAsP or InP/GaInP measurement
Fresh evaluation, per card, at the card's OWN default delta_xx (never the sweep's endpoint grid) -- both cards' own provenance.ranges["dot.gamma300"].note already document this by hand; the numbers below are computed programmatically.

| card | delta_xx (meV, card default) | g2_pulsed | collected flux (photons/s) | eligible | passes (<0.5) |
|---|---:|---:|---:|---|---|

## Grid
Grid complete: True
- `dot.delta_xx` in [4.0, 8.0] meV: sampled at 4, 8
- `dot.gamma300` in [6.0, 20.0] meV: sampled at 6, 20
- `irf_ps` in [50.0, 200.0] ps: sampled at 50, 200
- `thermal.T_hs` in [230.0, 300.0] K: sampled at 230, 250, 273, 300

## Per-temperature acceptance
| T_hs (K) | eligible | headline passes | pulsed g2 min | flux max (photons/s) | gamma300 threshold |
|---:|---:|---:|---:|---:|---|
| 230 | 0/0 | 0/0 | nan | nan | n/a |
| 250 | 0/0 | 0/0 | nan | nan | n/a |
| 273 | 0/0 | 0/0 | nan | nan | n/a |
| 300 | 0/0 | 0/0 | nan | nan | >=6 |

## What cooling buys

## Coverage
- **eligible fraction** (eligible/total; council review round 5, item 4 -- one word per quantity): 1/1 = 1.000
- pulsed collected-flux eligibility floor [A]: 1000 photons/s; rows excluded by this floor: 0
- diagnostic pooled g2 (below flux floor, not measurable): pulsed 0.3 / 0.3; g2_cw0 0.3 / 0.3; g2_cw0_raw 0.3 / 0.3
- maximum collected pulsed flux: nan photons/s; **flux_margin** (flux_max/floor; >1 means the floor is CLEARED): nan (flux_shortfall, DEPRECATED, its old floor/flux_max inverse framing: nan)
- **headline coverage over ALL SCHEDULED rows** (pulsed intrinsic g2(0) < 0.5; invalid/ineligible rows count as nonpassing -- the contract's PASS metric): 1/1 = 1.000
- **headline coverage over ELIGIBLE rows only** (same numerator, denominator restricted to eligible rows -- this is `coverage_over_eligible` in the VERDICT line, council review round 6 item 4, one name per quantity): 1/1 = 1.000
- pooled pulsed g2 median OVER ELIGIBLE ROWS ONLY (`g2_median_eligible` in the VERDICT line, the contract's median gate): 0.3
- pooled pulsed g2 median over ALL VALID/DIAGNOSTIC rows (eligible rows plus rows excluded ONLY by the flux floor; `diag_g2_median_diagnostic` in the VERDICT line): 0.3
- secondary coverage, CW intrinsic g2_cw0 < 0.5 (eligible rows, diagnostic only, does not gate PASS): 1/1 = 1.000
- secondary coverage, CW IRF-convolved g2_cw0_raw < 0.5 (eligible rows, diagnostic only, does not gate PASS): 1/1 = 1.000

## Per-card statistics
| card | role | diagnostic g2_pulsed min/median | diagnostic g2_cw0 min/median | diagnostic g2_cw0_raw min/median | eligible | flux-floor excluded | headline rows |
|---|---|---|---|---|---|---|---|
| edge-inp-gaasp-design | primary | 0.3 / 0.3 | 0.3 / 0.3 | 0.3 / 0.3 | 1/1 | 0 | 1 |

## Best diagnostic-g2 row and brightness decomposition
The dominant brightness limiter at the favourable diagnostic corner -- the smallest factor across the WHOLE chain (loading, t_X, S, and eta_total's own sub-factors; council review round 6, item 3) -- is the collection chain; the multiplicative chain that reproduces the reported collected pulsed flux is:

| factor | value |
|---|---:|
| loading = 1 - e^-mu (mu=nan) | nan |
| t_X (spectral transmission) | nan |
| S (confinement retention) | nan |
| eta_total (edge out-coupling: waveguide coupling x front/back facet split x facet transmission x propagation x NA, all in one factor) | nan |
| rep rate (Hz) | nan |
| **product x rep rate** | nan |
| reported collected_flux_pulsed_s | nan |
| self-check: relative difference | nan (non-finite inputs) |

`beta`/the combined front-facet factor/`T_facet`/`propagation`/`NA` are shown below for diagnosis only -- `beta`, the combined factor, `propagation` and `NA` are already folded into `eta_total` above exactly once each (their product reproduces eta_total, by construction of the combined factor -- NOT independent evidence, see below) and must NOT also be multiplied into the flux self-check:

| component (already inside eta_total) | value |
|---|---:|
| beta (waveguide coupling / spontaneous-emission factor) | nan |
| combined front/back-facet-and-transmission factor; BACK-SOLVED as the one factor missing from eta_total/(beta*eta_prop*eta_NA) -- device.py folds it into eta_total but never exposes it on its own; see the independent check below | nan |
| T_facet (raw Fresnel transmission, diagnostic only -- NOT necessarily an independent multiplicative step beyond the combined factor above; see the independent check below) | nan |
| propagation | nan |
| NA (numerical aperture) | nan |

**Independent check on the combined front-facet factor** (council review round 5, item 6, fixing a tautology): the value above is BACK-SOLVED -- the one factor missing once eta_total, beta, propagation and NA are all already known (eta_total = beta * <combined factor> * eta_prop * eta_NA) -- so a self-check comparing it against that same back-solving would only ever reproduce its own inputs. The genuinely independent check instead recomputes the combined factor FORWARD from this row's own facet transmission (T_facet) and back-facet reflectivity (R_back), trying every candidate formula introspected LIVE from fsim_core/waveguide.py's installed `edge_emission` source right now (never hardcoded, since that file's facet model is being edited concurrently and may change again) -- both as a fused current-model factor and as a legacy front-split-times-T_facet factor, since either grouping may be in effect: no candidate formula could be evaluated, back-solved = nan -- MISMATCH (agreement confirms the back-solved value really is the geometric front/back-facet-and-transmission factor waveguide.py computes, not some other quantity folded into eta_total).

## Diagnostic g2 landscape
| corner | pulsed g2 min | pooled diagnostic median | lever values | assumptions |
|---|---:|---:|---|---|
| edge-inp-gaasp-design (4 meV, gamma300=6 meV, irf=50 ps) | 0.3 | 0.3 | NA=0.75, R_back=0.95, L_um=250 | dot.gamma300; emission.lambda_nm |

## Assumptions required by any headline-passing corner
- `dot.gamma300`
- `emission.lambda_nm`

## Assumptions in plain words
**1. The collection window assumption makes flux look independent of the 300 K linewidth.** Every row's spectral filter is set automatically to exactly match the exciton line's own width at the operating temperature (device.py's `filter.auto_w` convention; `filter.auto_w_scale` can widen or narrow it, 1.0 = as-is here). A plain mathematical consequence of matching the window to the line's own width is that exactly HALF of the X-line's photons fall inside the window in every single row of this sweep (t_X = 0.5, column `t_x_pulsed`) -- no matter how wide or narrow the real (currently unmeasured) 300 K linewidth gamma300 turns out to be. That is why the reported collected photon flux does not move with gamma300 in this sweep: it is not that the real device would be insensitive to linewidth, it is that this filter-window convention cancels that sensitivity out by construction.

**2. The background-light assumption is borrowed from a different, colder device.** Every row carries a constant background term (`drive.b_res`) that is not measured on this platform: it is transferred from Reischle et al., Optics Express 16, 12771 (2008) (DOI 10.1364/OE.16.012771 [V], ledger anchor `reischle08-b-res-80k`), whose 80 K electrically driven single-photon source had about 88% real signal and 12% background light (signal fraction rho ~ 0.88). This sweep assumes the same 12% background fraction still applies at 300 K, on a different material system (InP/GaAsP or InP/GaInP edge emitters) than the one actually measured. No 300 K electrical background measurement exists for either card's platform.

**3. Why pulsed drive, not continuous-wave (CW) drive, is required.** At the best diagnostic operating point in this sweep, the intrinsic CW g2(0) is 0.3 (that alone would already satisfy the g2 < 0.5 single-photon criterion), but once a realistic single-photon detector's finite timing resolution (instrument response function, IRF) is folded in, the measured raw CW g2(0) rises to 0.3 -- ABOVE the 0.5 threshold. In plain terms: the antibunching dip this device produces under continuous drive is narrower in time than a real detector can resolve, so a CW measurement alone could never demonstrate single-photon emission on this platform. Pulsed (gated) operation sidesteps the detector's timing resolution and is therefore required, not optional.

## Evidence gate
- evidence_complete: False
- self_test_passed (self-test): True
- missing/incomplete evidence claims:
  - `y`: duplicate DOI collapses to one distinct source

## Fail reasons
- `evidence_incomplete`

## Sources
- docs/rt_edge_contract.md (acceptance gates, evidence status)
- verify/data/rt_edge_anchors.yaml (literature evidence ledger)
- cards/edge-inp-gaasp-design.yaml, cards/edge-inp-gainp-design.yaml (design.provenance.sources for every scalar)
