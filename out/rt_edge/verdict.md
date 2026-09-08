# RT edge-emitter acceptance sweep verdict

Generated 2026-09-08T02:10:45.950278+00:00; contract: docs/rt_edge_contract.md.

```
VERDICT: FAIL model=finite_pulse:true,tau_cap_density:false g2_min=nan g2_median_eligible=nan diag_g2_min=0.9817 diag_g2_median_diagnostic=0.9993 flux_max=582.8 flux_margin=0.5828 flux_shortfall_deprecated=1.716 median_pass=false coverage_over_eligible=nan eligible_fraction=0 eligible=0/768 flux_floor_excluded=768 evidence=incomplete conditional=false headline_coverage=0/768 headline_coverage_pulsed=0/384 headline_dedup_mismatch_groups=0 eligible_dedup=0/384 eligible_dedup_mismatch_groups=0 rows_scheduled=0/768 cw_raw_coverage=0/0 gamma300_pass_max=nan gamma300_threshold=n/a T_pass_min=none headline_by_T=230:0/192,250:0/192,273:0/192,300:0/192 headline_by_T_pulsed=230:0/96,250:0/96,273:0/96,300:0/96
```
```
CARD: edge-inp-gaasp-design role=primary diagnostic (below flux floor, not measurable) diag_g2_min=0.9981 diag_g2_median=0.9997 diag_g2_cw0_min=0.2928 diag_g2_cw0_median=0.7583 diag_g2_cw0_raw_min=0.9964 diag_g2_cw0_raw_median=0.9998 g2_cw_raw_min=nan g2_cw_raw_median=nan eligible=0/384 flux_floor_excluded=384 favorable_rows=0
```
```
CARD: edge-inp-gainp-design role=fallback diagnostic (below flux floor, not measurable) diag_g2_min=0.9817 diag_g2_median=0.9973 diag_g2_cw0_min=0.3114 diag_g2_cw0_median=0.6689 diag_g2_cw0_raw_min=0.9654 diag_g2_cw0_raw_median=0.998 g2_cw_raw_min=nan g2_cw_raw_median=nan eligible=0/384 flux_floor_excluded=384 favorable_rows=0
```

## Literature ceiling
In the reviewed literature set of this repository (six papers, ../_goal/paper_digests.md) and the anchors ledger, no electrically driven III-V single-dot g2(0) at 300 K is reported; the best electrical result in that set is Reischle et al., Optics Express 16, 12771 (2008) at 80 K (QD C): g2(0) = 0.43 raw (still IRF-broadened), 0.03 after background correction with the ~500 ps detector IRF also deconvolved. The best reported 300 K single-dot values are optically pumped: g2(0) ~ 0.5-0.57 (Laferriere et al. 2023, InAsP/InP nanowire dot, g2(0) = 0.57 at 300 K). This sweep's pooled pulsed g2_min=nan, g2_median(eligible)=nan is reported against that electrical-vs-optical literature picture, not as a claim of an existing electrical 300 K result to exceed.

**Convention-matched comparison (council review round 5, item 5; peer review finding 9).** The paragraph above juxtaposed this sweep's IRF-FREE intrinsic g2_min against Reischle's RAW dip (0.43 [V], still IRF-broadened; ledger anchor `reischle08-g2-80k`) -- not a convention-matched comparison. The correct convention-matched anchor is Reischle's IRF-DECONVOLVED-but-background-included value, g2_b(0) = 0.25 [V] +/- 0.05 (QD C, 80 K; ledger anchor `reischle08-g2-80k-deconvolved`, ../_goal/paper_digests.md line ~36), since this sweep's g2_op is likewise background-included (via drive.b_res/rho) and never IRF-convolved for the pulsed metric. Pulsed peak-area g2(0) and CW zero-delay g2(0) remain distinct observables even after matching background and IRF.
No finite pooled g2_min is available in this run to compare against Reischle's deconvolved value.

## Conditional on the 300 K linewidth
The 300 K single-dot linewidth gamma300 is not a platform ceiling -- it is an unmeasured [E]-class quantity (see Assumptions below) whose actual value decides whether the headline gate passes. `gamma300_pass_max` (council review round 5, item 2) is the HIGHEST SAMPLED linewidth, among this determination's own 8-point grid (6, 6.5, 7, 8, 10, 12, 16, 20 meV -- finer than the sweep.csv grid's own 2-point dot.gamma300 endpoints, and never itself written to sweep.csv), at which any eligible row (any lever/delta_xx combination) still clears pulsed intrinsic g2(0) < 0.5; `gamma300_threshold` brackets the true (continuous, unsampled) threshold between that value and the next sampled value above it that already fails -- the true threshold lies somewhere inside the bracket, at the delta_xx shown for that card's gamma300_pass_max row.

**When nothing passes (council review round 6, item 2).** If NO sampled gamma300 clears the headline gate for a card, `gamma300_threshold` no longer prints a numeric bracket like `<6` -- that would misleadingly imply the card passes somewhere below 6 meV, when in fact it never passes anywhere in this grid. It instead prints the REASON class: `none(flux)` if every sampled row is below the collected-flux eligibility floor (headline_pass can never be True there, no matter how favourable g2 would be), or `none(g2)` if at least one sampled row IS flux-eligible but pulsed intrinsic g2(0) >= 0.5 at every sampled gamma300 (a genuine physics shortfall, not an eligibility one); `none(flux,g2)` reports both blocking conditions when neither a flux-eligible row nor a diagnostic g2 < 0.5 row exists.

`gamma300_pass_max` (pooled, both cards): nan meV; `gamma300_threshold` (pooled): n/a meV

| card | gamma300_pass_max (meV) | gamma300_threshold bracket (meV) | delta_xx | lever values | assumptions |
|---|---:|---:|---:|---|---|
| edge-inp-gaasp-design | n/a (no headline-passing sample in this grid) | none(flux,g2) | | | |
| edge-inp-gainp-design | n/a (no headline-passing sample in this grid) | none(flux,g2) | | | |

Restated (council review round 5, item 2): for each card, `gamma300_pass_max` is the highest SAMPLED linewidth at which the card passes; the true threshold lies between the bracket's two values (at the delta_xx shown above for that row).

### Verified 6.5 meV anchor (Chatzarakis et al., Phys. Rev. Applied 20, 034011, 2023) -- an InAs/GaAs (211)B single-dot linewidth measurement, used here as a distinct-material class proxy [E] for the InP-dot family, not a direct InP/GaAsP or InP/GaInP measurement
Fresh evaluation, per card, at the card's OWN default delta_xx (never the sweep's endpoint grid) -- both cards' own provenance.ranges["dot.gamma300"].note already document this by hand; the numbers below are computed programmatically.

| card | delta_xx (meV, card default) | g2_pulsed | collected flux (photons/s) | eligible | passes (<0.5) |
|---|---:|---:|---:|---|---|
| edge-inp-gaasp-design | 5 | 0.9999 | 2.285 | False | no |
| edge-inp-gainp-design | 7 | 0.9979 | 76.68 | False | no |

at gamma300 = 6.5 meV (this card's own delta_xx = 5 meV), card `edge-inp-gaasp-design` does not pass the headline gate (g2_pulsed = 0.9999, collected flux = 2.285 photons/s); at gamma300 = 6.5 meV (this card's own delta_xx = 7 meV), card `edge-inp-gainp-design` does not pass the headline gate (g2_pulsed = 0.9979, collected flux = 76.68 photons/s).

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

Grid sizing (pkg5b item 1): one row measured at 0.494 s; extrapolated to 25.3 minutes for the full four-model-combination grid (768 rows x 4 combinations), against a 90-minute budget -- the full grid was used for all four model combinations.

## Model sensitivity
Two opt-in physics-correction switches (`drive.finite_pulse`, peer-review-triage.md finding 1; `ret.tau_cap_scales_with_density`, finding 1b) are now genuine sweep axes (pkg5b) instead of both being silently fixed at their card-default (false, false), which every prior package used. The headline above -- headline_coverage_pulsed, T_pass_min, headline_by_T*, g2_min, flux_max -- is computed over ONLY `model=finite_pulse:true,tau_cap_density:false` rows (the corrected, no-cancellation model). The other three combinations never gate PASS; they are reported here for comparison, pooled across both cards at each combination's own best corner.

| model (finite_pulse, tau_cap_density) | eligible | headline passes | pulsed g2 min | flux max (photons/s) |
|---|---:|---:|---:|---:|
| finite_pulse:false,tau_cap_density:false | 0/768 | 0 | nan | 449.4 |
| finite_pulse:false,tau_cap_density:true | 392/768 | 170 | 0.2995 | 4767 |
| finite_pulse:true,tau_cap_density:false (headline) | 0/768 | 0 | nan | 582.8 |
| finite_pulse:true,tau_cap_density:true | 496/768 | 0 | 0.8291 | 6164 |

## Per-temperature acceptance
| T_hs (K) | eligible | headline passes | deduplicated (IRF axis collapsed) | pulsed g2 min | flux max (photons/s) | gamma300 threshold |
|---:|---:|---:|---:|---:|---:|---|
| 230 | 0/192 | 0/192 | 0/96 | nan | 582.8 | n/a |
| 250 | 0/192 | 0/192 | 0/96 | nan | 302.4 | n/a |
| 273 | 0/192 | 0/192 | 0/96 | nan | 156.4 | n/a |
| 300 | 0/192 | 0/192 | 0/96 | nan | 80.43 | n/a |

## What cooling buys
At T_hs=230 K, the favourable corner has retention S=0.001348, linewidth Gamma(T)=3.55 meV, and window background b_e=0.07775.
At T_hs=250 K, the favourable corner has retention S=0.0006996, linewidth Gamma(T)=4.213 meV, and window background b_e=0.06934.
At T_hs=273 K, the favourable corner has retention S=0.0003624, linewidth Gamma(T)=5.015 meV, and window background b_e=0.05962.
At T_hs=300 K, the favourable corner has retention S=0.0001873, linewidth Gamma(T)=6 meV, and window background b_e=0.04887.

## Coverage
- **eligible fraction** (eligible/total; council review round 5, item 4 -- one word per quantity): 0/768 = 0.000
- pulsed collected-flux eligibility floor [A]: 1000 photons/s; rows excluded by this floor: 768
- diagnostic pooled g2 (below flux floor, not measurable): pulsed 0.9817 / 0.9993; g2_cw0 0.2928 / 0.7142; g2_cw0_raw 0.9654 / 0.9995
- maximum collected pulsed flux: 582.8 photons/s; **flux_margin** (flux_max/floor; >1 means the floor is CLEARED): 0.5828 (flux_shortfall, DEPRECATED, its old floor/flux_max inverse framing: 1.716)
- **rows_scheduled** (pulsed intrinsic g2(0) < 0.5; invalid/ineligible rows count as nonpassing; raw grid count, every irf_ps sample counted separately -- `headline_coverage` in the VERDICT line, kept for backward compatibility): 0/768 = 0.000
- **headline_coverage_pulsed** (peer review finding 6, pkg5 fix item 2: same numerator rule as rows_scheduled, but the irf_ps axis is deduplicated first. The pulsed g2 sub-result is IRF-independent; `headline_pass` also requires CW eligibility, which is IRF-convolved, so the deduplicated count is computed with `all()` over the IRF axis and any disagreement is reported as `headline_dedup_mismatch_groups` (currently 0) rather than resolved by picking whichever irf_ps row happened to be seen first): 0/384 = 0.000. Coverage is the fraction of a chosen endpoint grid that passes, not a fabrication-yield probability or a confidence level.
- **headline coverage over ELIGIBLE rows only**: n/a (0 eligible rows)
- pooled pulsed g2 median OVER ELIGIBLE ROWS ONLY (`g2_median_eligible` in the VERDICT line; diagnostic only, it never appears in fail_reasons): nan
- pooled pulsed g2 median over ALL VALID/DIAGNOSTIC rows (eligible rows plus rows excluded ONLY by the flux floor; `diag_g2_median_diagnostic` in the VERDICT line): 0.9993
- secondary coverage, CW intrinsic g2_cw0 < 0.5 (eligible rows, diagnostic only, does not gate PASS): 0/0 = 0.000
- secondary coverage, CW IRF-convolved g2_cw0_raw < 0.5 (eligible rows, diagnostic only, does not gate PASS): 0/0 = 0.000

## Per-card statistics
| card | role | diagnostic g2_pulsed min/median | diagnostic g2_cw0 min/median | diagnostic g2_cw0_raw min/median | eligible | flux-floor excluded | headline rows |
|---|---|---|---|---|---|---|---|
| edge-inp-gaasp-design | primary | 0.9981 / 0.9997 | 0.2928 / 0.7583 | 0.9964 / 0.9998 | 0/384 | 384 | 0 |
| edge-inp-gainp-design | fallback | 0.9817 / 0.9973 | 0.3114 / 0.6689 | 0.9654 / 0.998 | 0/384 | 384 | 0 |

## Why no row is eligible
Every row is below the 1000 photons/s collected-flux floor; the grid maximum is only 582.8 photons/s (floor/maximum = 1.716).
The dominant brightness limiter at the favourable diagnostic corner -- the smallest factor across the WHOLE chain (loading, t_X, S, and eta_total's own sub-factors; council review round 6, item 3) -- is S (confinement retention); the multiplicative chain that reproduces the reported collected pulsed flux is:

| factor | value |
|---|---:|
| loading = 1 - e^-mu (mu=0.5445) | 0.419876 |
| t_X (spectral transmission) | 0.5 |
| S (confinement retention) | 0.00134802 |
| eta_total (edge out-coupling: waveguide coupling x facet escape (mid-ridge ray series, propagation included) x NA) | 0.0198516 |
| rep rate (Hz) | 8e+07 |
| **product x rep rate** | 449.442 photons/s |
| reported collected_flux_pulsed_s | 578.843 photons/s |
| self-check: relative difference | 22.3553% (FAIL, exceeds 1%) |

`beta`/`eta_facet`/`T_facet`/`NA` are shown below for diagnosis only -- `beta`, `eta_facet` and `NA` are already folded into `eta_total` above exactly once each (their product reproduces eta_total, by construction of eta_facet -- NOT independent evidence, see below) and must NOT also be multiplied into the flux self-check. Single-pass propagation is NOT listed here as a separate multiplicative row (peer-review pkg2 facet fix, 2026-09-07): it is reported below as an informational line only, already folded inside `eta_facet`'s ray series.

| component (already inside eta_total) | value |
|---|---:|
| beta (waveguide coupling / spontaneous-emission factor) | 0.0296973 |
| eta_facet (front/back-facet split, transmission, and single-pass propagation all fused into one ray-series factor; BACK-SOLVED as the one factor missing from eta_total/(beta*eta_NA) -- device.py folds it into eta_total but never exposes it on its own; see the independent check below) | 0.784015 |
| T_facet (raw Fresnel transmission, diagnostic only -- NOT an independent multiplicative step beyond eta_facet above; see the independent check below) | 0.719581 |
| NA (numerical aperture) | 0.852617 |

single-pass propagation (already inside the facet factor): 0.882497

**Independent check on eta_facet** (council review round 5, item 6, fixing a tautology; peer-review pkg2 facet fix, 2026-09-07): the value above is BACK-SOLVED -- the one factor missing once eta_total, beta and NA are all already known (eta_total = beta * eta_facet * eta_NA) -- so a self-check comparing it against that same back-solving would only ever reproduce its own inputs. The genuinely independent check instead recomputes eta_facet FORWARD from this row's own facet transmission (T_facet), back-facet reflectivity (R_back), and ridge loss (alpha_cm), by a direct call to fsim_core.waveguide.facet_escape_fraction -- the SAME pure function edge_emission() itself calls (single source of truth, never a duplicated or scraped formula), so this is a wiring/regression check, not a physics re-derivation: forward = 0.784015 via ray-series-midpoint, back-solved = 0.784015 -- MATCH (agreement confirms the back-solved value really is the eta_facet factor waveguide.py computes, not some other quantity folded into eta_total).

## Diagnostic g2 landscape
| corner | pulsed g2 min | pooled diagnostic median | lever values | assumptions |
|---|---:|---:|---|---|
| edge-inp-gainp-design (8 meV, gamma300=6 meV, irf=50 ps) | 0.981747 | 0.999275 | NA=0.8, R_back=0.95, L_um=250 | dot.gamma300; drive.b_res; drive.diode.tau_pulse_ns; drive.duty (rep-rate-derived); emission.lambda_nm; ret.system.barrier.x_al |

## Assumptions required by any headline-passing corner
- (no row satisfies the headline gate; see fail_reasons below)

## Assumptions in plain words
**1. The collection window assumption makes flux look independent of the 300 K linewidth.** Every row's spectral filter is set automatically to exactly match the exciton line's own width at the operating temperature (device.py's `filter.auto_w` convention; `filter.auto_w_scale` can widen or narrow it, 1.0 = as-is here). A plain mathematical consequence of matching the window to the line's own width is that exactly HALF of the X-line's photons fall inside the window in every single row of this sweep (t_X = 0.5, column `t_x_pulsed`) -- no matter how wide or narrow the real (currently unmeasured) 300 K linewidth gamma300 turns out to be. That is why the reported collected photon flux does not move with gamma300 in this sweep: it is not that the real device would be insensitive to linewidth, it is that this filter-window convention cancels that sensitivity out by construction.

**2. The background-light assumption is borrowed from a different, colder device.** Every row carries a constant background term (`drive.b_res`) that is not measured on this platform: it is transferred from Reischle et al., Optics Express 16, 12771 (2008) (DOI 10.1364/OE.16.012771 [V], ledger anchor `reischle08-b-res-80k`), whose 80 K electrically driven single-photon source had about 88% real signal and 12% background light (signal fraction rho ~ 0.88). This sweep assumes the same 12% background fraction still applies at 300 K, on a different material system (InP/GaAsP or InP/GaInP edge emitters) than the one actually measured. No 300 K electrical background measurement exists for either card's platform.

**3. CW versus pulsed measurement at the best diagnostic point.** At the best diagnostic operating point in this sweep, the intrinsic CW g2(0) is 0.311 (that alone would already satisfy the g2 < 0.5 single-photon criterion), but once a realistic single-photon detector's finite timing resolution (instrument response function, IRF) is folded in, the measured raw CW g2(0) rises to 0.965 -- ABOVE the 0.5 threshold, so a CW measurement alone cannot demonstrate single-photon emission at the IRF values sampled here (50-200 ps); pulsed, gated operation is required at these IRF values; a faster detector, a different gate or different physical rates could change this.

## Evidence gate
- evidence_complete: False
- self_test_passed (self-test): True
- missing/incomplete evidence claims:
  - `reischle2008_g2_80K`: only 1 distinct verified primary source(s), need 2
  - `hkust_inp_gaasp_wavelength`: only 1 distinct verified primary source(s), need 2
  - `residual_background_80K`: only 1 distinct verified primary source(s), need 2

## Fail reasons
- `evidence_incomplete`
- `no_eligible_rows`

## Sources
- docs/rt_edge_contract.md (acceptance gates, evidence status)
- verify/data/rt_edge_anchors.yaml (literature evidence ledger)
- cards/edge-inp-gaasp-design.yaml, cards/edge-inp-gainp-design.yaml (design.provenance.sources for every scalar)
