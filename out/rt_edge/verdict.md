# RT edge-emitter acceptance sweep verdict

Generated 2026-09-06T22:34:19.910122+00:00; contract: docs/rt_edge_contract.md.

```
VERDICT: FAIL g2_min=0.3214 g2_median_eligible=0.6373 diag_g2_min=0.3214 diag_g2_median_diagnostic=0.6706 flux_max=6890 flux_margin=6.89 flux_shortfall_deprecated=0.1451 median_pass=false coverage_over_eligible=0.3066 eligible_fraction=0.5521 eligible=424/768 flux_floor_excluded=344 evidence=incomplete conditional=true headline_coverage=130/768 cw_raw_coverage=0/424 gamma300_pass_max=10 gamma300_threshold=10-12 T_pass_min=230 headline_by_T=230:60/192,250:34/192,273:24/192,300:12/192
```
```
CARD: edge-inp-gaasp-design role=primary eligible rows: 104/384 g2_pulsed_min=0.4188 g2_pulsed_median=0.6696 diag_g2_min=0.4188 diag_g2_median=0.7201 diag_g2_cw0_min=0.2834 diag_g2_cw0_median=0.6425 diag_g2_cw0_raw_min=0.9269 diag_g2_cw0_raw_median=0.9912 g2_cw_raw_min=0.9269 g2_cw_raw_median=0.987 eligible=104/384 flux_floor_excluded=280 favorable_rows=26
```
```
CARD: edge-inp-gainp-design role=fallback eligible rows: 320/384 g2_pulsed_min=0.3214 g2_pulsed_median=0.6373 diag_g2_min=0.3214 diag_g2_median=0.6373 diag_g2_cw0_min=0.2836 diag_g2_cw0_median=0.6377 diag_g2_cw0_raw_min=0.7345 diag_g2_cw0_raw_median=0.9596 g2_cw_raw_min=0.7345 g2_cw_raw_median=0.9555 eligible=320/384 flux_floor_excluded=64 favorable_rows=104
```

## Literature ceiling
In the reviewed literature set of this repository (six papers, ../_goal/paper_digests.md) and the anchors ledger, no electrically driven III-V single-dot g2(0) at 300 K is reported; the best electrical result in that set is Reischle et al., Optics Express 16, 12771 (2008) at 80 K (QD C): g2(0) = 0.43 raw (still IRF-broadened), 0.03 after background correction with the ~500 ps detector IRF also deconvolved. The best reported 300 K single-dot values are optically pumped: g2(0) ~ 0.5-0.57 (Laferriere et al. 2023, InAsP/InP nanowire dot, g2(0) = 0.57 at 300 K). This sweep's pooled pulsed g2_min=0.3214, g2_median(eligible)=0.6373 is reported against that electrical-vs-optical literature picture, not as a claim of an existing electrical 300 K result to exceed.

**Like-for-like comparison (council review round 5, item 5).** The paragraph above juxtaposed this sweep's IRF-FREE intrinsic g2_min against Reischle's RAW dip (0.43 [V], still IRF-broadened; ledger anchor `reischle08-g2-80k`) -- not a like-for-like convention match. The correct like-for-like anchor is Reischle's IRF-DECONVOLVED-but-background-included value, g2_b(0) = 0.25 [V] +/- 0.05 (QD C, 80 K; ledger anchor `reischle08-g2-80k-deconvolved`, ../_goal/paper_digests.md line ~36), since this sweep's g2_op is likewise background-included (via drive.b_res/rho) and never IRF-convolved for the pulsed metric.
On matching conventions, the 80 K Reischle device is about 1.3x better (lower g2(0)) than this sweep's best 300 K corner (g2_min=0.3214 vs 0.25 +/- 0.05).

## Conditional on the 300 K linewidth
The 300 K single-dot linewidth gamma300 is not a platform ceiling -- it is an unmeasured [E]-class quantity (see Assumptions below) whose actual value decides whether the headline gate passes. `gamma300_pass_max` (council review round 5, item 2) is the HIGHEST SAMPLED linewidth, among this determination's own 8-point grid (6, 6.5, 7, 8, 10, 12, 16, 20 meV -- finer than the sweep.csv grid's own 2-point dot.gamma300 endpoints, and never itself written to sweep.csv), at which any eligible row (any lever/delta_xx combination) still clears pulsed intrinsic g2(0) < 0.5; `gamma300_threshold` brackets the true (continuous, unsampled) threshold between that value and the next sampled value above it that already fails -- the true threshold lies somewhere inside the bracket, at the delta_xx shown for that card's gamma300_pass_max row.

**When nothing passes (council review round 6, item 2).** If NO sampled gamma300 clears the headline gate for a card, `gamma300_threshold` no longer prints a numeric bracket like `<6` -- that would misleadingly imply the card passes somewhere below 6 meV, when in fact it never passes anywhere in this grid. It instead prints the REASON class: `none(flux)` if every sampled row is below the collected-flux eligibility floor (headline_pass can never be True there, no matter how favourable g2 would be), or `none(g2)` if at least one sampled row IS flux-eligible but pulsed intrinsic g2(0) >= 0.5 at every sampled gamma300 (a genuine physics shortfall, not an eligibility one); `none(flux,g2)` reports both blocking conditions when neither a flux-eligible row nor a diagnostic g2 < 0.5 row exists.

`gamma300_pass_max` (pooled, both cards): 10 meV; `gamma300_threshold` (pooled): 10-12 meV

| card | gamma300_pass_max (meV) | gamma300_threshold bracket (meV) | delta_xx | lever values | assumptions |
|---|---:|---:|---:|---|---|
| edge-inp-gaasp-design | n/a (no headline-passing sample in this grid) | none(g2) | | | |
| edge-inp-gainp-design | 10 | 10-12 | 7 meV | NA=0.75, R_back=0.95, L_um=250 | dot.gamma300; drive.b_res; drive.diode.tau_pulse_ns; drive.duty (rep-rate-derived); emission.lambda_nm; ret.system.barrier.x_al |

Restated (council review round 5, item 2): for each card, `gamma300_pass_max` is the highest SAMPLED linewidth at which the card passes; the true threshold lies between the bracket's two values (at the delta_xx shown above for that row).

### Verified 6.5 meV anchor (Chatzarakis et al., Phys. Rev. Applied 20, 034011, 2023) -- an InAs/GaAs (211)B single-dot linewidth measurement, used here as a distinct-material class proxy [E] for the InP-dot family, not a direct InP/GaAsP or InP/GaInP measurement
Fresh evaluation, per card, at the card's OWN default delta_xx (never the sweep's endpoint grid) -- both cards' own provenance.ranges["dot.gamma300"].note already document this by hand; the numbers below are computed programmatically.

| card | delta_xx (meV, card default) | g2_pulsed | collected flux (photons/s) | eligible | passes (<0.5) |
|---|---:|---:|---:|---|---|
| edge-inp-gaasp-design | 5 | 0.6523 | 474.5 | False | no |
| edge-inp-gainp-design | 7 | 0.4563 | 2308 | True | yes |

at gamma300 = 6.5 meV (this card's own delta_xx = 5 meV), card `edge-inp-gaasp-design` does not pass the headline gate (g2_pulsed = 0.6523, collected flux = 474.5 photons/s); at gamma300 = 6.5 meV (this card's own delta_xx = 7 meV), card `edge-inp-gainp-design` PASSES the headline gate (g2_pulsed = 0.4563, collected flux = 2308 photons/s).

## Grid
Grid complete: True
- `dot.delta_xx` in [4.0, 8.0] meV: sampled at 4, 8
- `dot.gamma300` in [6.0, 20.0] meV: sampled at 6, 20
- `irf_ps` in [50.0, 200.0] ps: sampled at 50, 200
- `thermal.T_hs` in [230.0, 300.0] K: sampled at 230, 250, 273, 300

Collection-lever axes (per-card, council review round 4 item 1):
- `edge-inp-gaasp-design` `emission.NA`: sampled at 0.5, 0.8, 0.75
- `edge-inp-gaasp-design` `emission.R_back`: sampled at 0, 0.95
- `edge-inp-gaasp-design` `emission.L_um`: sampled at 250, 500
- `edge-inp-gaasp-design` lever combinations evaluated: 12
- `edge-inp-gainp-design` `emission.NA`: sampled at 0.5, 0.8, 0.75
- `edge-inp-gainp-design` `emission.R_back`: sampled at 0, 0.95
- `edge-inp-gainp-design` `emission.L_um`: sampled at 250, 500
- `edge-inp-gainp-design` lever combinations evaluated: 12

## Per-temperature acceptance
| T_hs (K) | eligible | headline passes | pulsed g2 min | flux max (photons/s) | gamma300 threshold |
|---:|---:|---:|---:|---:|---|
| 230 | 144/192 | 60/192 | 0.3214 | 6890 | 10-12 |
| 250 | 136/192 | 34/192 | 0.376 | 5168 | n/a |
| 273 | 96/192 | 24/192 | 0.3775 | 3639 | n/a |
| 300 | 48/192 | 12/192 | 0.3986 | 2417 | 7-8 |

## What cooling buys
At T_hs=230 K, the favourable corner has retention S=0.01311, linewidth Gamma(T)=3.55 meV, and window background b_e=0.01479.
At T_hs=250 K, the favourable corner has retention S=0.009837, linewidth Gamma(T)=4.213 meV, and window background b_e=0.03107.
At T_hs=273 K, the favourable corner has retention S=0.006935, linewidth Gamma(T)=5.015 meV, and window background b_e=0.01963.
At T_hs=300 K, the favourable corner has retention S=0.004617, linewidth Gamma(T)=6 meV, and window background b_e=0.01247.

## Coverage
- **eligible fraction** (eligible/total; council review round 5, item 4 -- one word per quantity): 424/768 = 0.552
- pulsed collected-flux eligibility floor [A]: 1000 photons/s; rows excluded by this floor: 344
- diagnostic pooled g2 (below flux floor, not measurable): pulsed 0.3214 / 0.6706; g2_cw0 0.2834 / 0.6377; g2_cw0_raw 0.7345 / 0.9845
- maximum collected pulsed flux: 6890 photons/s; **flux_margin** (flux_max/floor; >1 means the floor is CLEARED): 6.89 (flux_shortfall, DEPRECATED, its old floor/flux_max inverse framing: 0.1451)
- **headline coverage over ALL SCHEDULED rows** (pulsed intrinsic g2(0) < 0.5; invalid/ineligible rows count as nonpassing -- the contract's PASS metric): 130/768 = 0.169
- **headline coverage over ELIGIBLE rows only** (same numerator, denominator restricted to eligible rows -- this is `coverage_over_eligible` in the VERDICT line, council review round 6 item 4, one name per quantity): 130/424 = 0.307
- pooled pulsed g2 median OVER ELIGIBLE ROWS ONLY (`g2_median_eligible` in the VERDICT line, the contract's median gate): 0.6373
- pooled pulsed g2 median over ALL VALID/DIAGNOSTIC rows (eligible rows plus rows excluded ONLY by the flux floor; `diag_g2_median_diagnostic` in the VERDICT line): 0.6706
- secondary coverage, CW intrinsic g2_cw0 < 0.5 (eligible rows, diagnostic only, does not gate PASS): 176/424 = 0.415
- secondary coverage, CW IRF-convolved g2_cw0_raw < 0.5 (eligible rows, diagnostic only, does not gate PASS): 0/424 = 0.000

## Per-card statistics
| card | role | diagnostic g2_pulsed min/median | diagnostic g2_cw0 min/median | diagnostic g2_cw0_raw min/median | eligible | flux-floor excluded | headline rows |
|---|---|---|---|---|---|---|---|
| edge-inp-gaasp-design | primary | 0.4188 / 0.7201 | 0.2834 / 0.6425 | 0.9269 / 0.9912 | 104/384 | 280 | 26 |
| edge-inp-gainp-design | fallback | 0.3214 / 0.6373 | 0.2836 / 0.6377 | 0.7345 / 0.9596 | 320/384 | 64 | 104 |

## Best diagnostic-g2 row and brightness decomposition
The dominant brightness limiter at the favourable diagnostic corner -- the smallest factor across the WHOLE chain (loading, t_X, S, and eta_total's own sub-factors; council review round 6, item 3) -- is S (confinement retention); the multiplicative chain that reproduces the reported collected pulsed flux is:

| factor | value |
|---|---:|
| loading = 1 - e^-mu (mu=1.007) | 0.634719 |
| t_X (spectral transmission) | 0.5 |
| S (confinement retention) | 0.0131112 |
| eta_total (edge out-coupling: waveguide coupling x front/back facet split x facet transmission x propagation x NA, all in one factor) | 0.0206985 |
| rep rate (Hz) | 8e+07 |
| **product x rep rate** | 6890.03 photons/s |
| reported collected_flux_pulsed_s | 6890.03 photons/s |
| self-check: relative difference | 0.0000% (PASS, within 1%) |

`beta`/the combined front-facet factor/`T_facet`/`propagation`/`NA` are shown below for diagnosis only -- `beta`, the combined factor, `propagation` and `NA` are already folded into `eta_total` above exactly once each (their product reproduces eta_total, by construction of the combined factor -- NOT independent evidence, see below) and must NOT also be multiplied into the flux self-check:

| component (already inside eta_total) | value |
|---|---:|
| beta (waveguide coupling / spontaneous-emission factor) | 0.0292412 |
| combined front/back-facet-and-transmission factor; BACK-SOLVED as the one factor missing from eta_total/(beta*eta_prop*eta_NA) -- device.py folds it into eta_total but never exposes it on its own; see the independent check below | 0.935012 |
| T_facet (raw Fresnel transmission, diagnostic only -- NOT necessarily an independent multiplicative step beyond the combined factor above; see the independent check below) | 0.719371 |
| propagation | 0.882497 |
| NA (numerical aperture) | 0.857853 |

**Independent check on the combined front-facet factor** (council review round 5, item 6, fixing a tautology): the value above is BACK-SOLVED -- the one factor missing once eta_total, beta, propagation and NA are all already known (eta_total = beta * <combined factor> * eta_prop * eta_NA) -- so a self-check comparing it against that same back-solving would only ever reproduce its own inputs. The genuinely independent check instead recomputes the combined factor FORWARD from this row's own facet transmission (T_facet) and back-facet reflectivity (R_back), trying every candidate formula introspected LIVE from fsim_core/waveguide.py's installed `edge_emission` source right now (never hardcoded, since that file's facet model is being edited concurrently and may change again) -- both as a fused current-model factor and as a legacy front-split-times-T_facet factor, since either grouping may be in effect: forward = 0.935012 via `facet_factor = T / (T + (1 - R_back))` (fused/current-model convention: already includes T_facet), back-solved = 0.935012 -- MATCH (agreement confirms the back-solved value really is the geometric front/back-facet-and-transmission factor waveguide.py computes, not some other quantity folded into eta_total).

## Diagnostic g2 landscape
| corner | pulsed g2 min | pooled diagnostic median | lever values | assumptions |
|---|---:|---:|---|---|
| edge-inp-gainp-design (8 meV, gamma300=6 meV, irf=50 ps) | 0.321359 | 0.670574 | NA=0.8, R_back=0.95, L_um=250 | dot.gamma300; drive.b_res; drive.diode.tau_pulse_ns; drive.duty (rep-rate-derived); emission.lambda_nm; ret.system.barrier.x_al |

## Assumptions required by any headline-passing corner
- `dot.gamma300`
- `drive.b_res`
- `drive.diode.tau_pulse_ns`
- `drive.duty (rep-rate-derived)`
- `emission.lambda_nm`
- `ret.system.barrier.x_al`

## Assumptions in plain words
**1. The collection window assumption makes flux look independent of the 300 K linewidth.** Every row's spectral filter is set automatically to exactly match the exciton line's own width at the operating temperature (device.py's `filter.auto_w` convention; `filter.auto_w_scale` can widen or narrow it, 1.0 = as-is here). A plain mathematical consequence of matching the window to the line's own width is that exactly HALF of the X-line's photons fall inside the window in every single row of this sweep (t_X = 0.5, column `t_x_pulsed`) -- no matter how wide or narrow the real (currently unmeasured) 300 K linewidth gamma300 turns out to be. That is why the reported collected photon flux does not move with gamma300 in this sweep: it is not that the real device would be insensitive to linewidth, it is that this filter-window convention cancels that sensitivity out by construction.

**2. The background-light assumption is borrowed from a different, colder device.** Every row carries a constant background term (`drive.b_res`) that is not measured on this platform: it is transferred from Reischle et al., Optics Express 16, 12771 (2008) (DOI 10.1364/OE.16.012771 [V], ledger anchor `reischle08-b-res-80k`), whose 80 K electrically driven single-photon source had about 88% real signal and 12% background light (signal fraction rho ~ 0.88). This sweep assumes the same 12% background fraction still applies at 300 K, on a different material system (InP/GaAsP or InP/GaInP edge emitters) than the one actually measured. No 300 K electrical background measurement exists for either card's platform.

**3. Why pulsed drive, not continuous-wave (CW) drive, is required.** At the best diagnostic operating point in this sweep, the intrinsic CW g2(0) is 0.284 (that alone would already satisfy the g2 < 0.5 single-photon criterion), but once a realistic single-photon detector's finite timing resolution (instrument response function, IRF) is folded in, the measured raw CW g2(0) rises to 0.734 -- ABOVE the 0.5 threshold. In plain terms: the antibunching dip this device produces under continuous drive is narrower in time than a real detector can resolve, so a CW measurement alone could never demonstrate single-photon emission on this platform. Pulsed (gated) operation sidesteps the detector's timing resolution and is therefore required, not optional.

## Evidence gate
- evidence_complete: False
- self_test_passed (self-test): True
- missing/incomplete evidence claims:
  - `reischle2008_g2_80K`: only 1 distinct verified primary source(s), need 2
  - `hkust_inp_gaasp_wavelength`: only 1 distinct verified primary source(s), need 2
  - `residual_background_80K`: only 1 distinct verified primary source(s), need 2

## Fail reasons
- `evidence_incomplete`

## Sources
- docs/rt_edge_contract.md (acceptance gates, evidence status)
- verify/data/rt_edge_anchors.yaml (literature evidence ledger)
- cards/edge-inp-gaasp-design.yaml, cards/edge-inp-gainp-design.yaml (design.provenance.sources for every scalar)
