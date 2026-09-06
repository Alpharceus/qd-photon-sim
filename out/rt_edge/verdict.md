# RT edge-emitter acceptance sweep verdict

Generated 2026-09-06T07:06:37.352639+00:00; contract: docs/rt_edge_contract.md.

```
VERDICT: FAIL g2_min=0.4305 g2_median=0.7076 diag_g2_min=0.4305 diag_g2_median=0.7409 flux_max=1111 flux_shortfall=0.9001 median_pass=false coverage=0.02083 eligible=16/192 flux_floor_excluded=176 evidence=incomplete conditional=true headline_coverage=4/192 cw_raw_coverage=0/192 gamma300_pass_max=6
```
```
CARD: edge-inp-gaasp-design role=primary g2_pulsed_min=nan g2_pulsed_median=nan diagnostic (below flux floor, not measurable) diag_g2_min=0.5306 diag_g2_median=0.7856 diag_g2_cw0_min=0.4017 diag_g2_cw0_median=0.8252 diag_g2_cw0_raw_min=0.9784 diag_g2_cw0_raw_median=0.9983 g2_cw_raw_min=nan g2_cw_raw_median=nan eligible=0/96 flux_floor_excluded=96 favorable_rows=0
```
```
CARD: edge-inp-gainp-design role=fallback g2_pulsed_min=0.4305 g2_pulsed_median=0.7076 diagnostic (below flux floor, not measurable) diag_g2_min=0.4305 diag_g2_median=0.7076 diag_g2_cw0_min=0.3915 diag_g2_cw0_median=0.8086 diag_g2_cw0_raw_min=0.9003 diag_g2_cw0_raw_median=0.9898 g2_cw_raw_min=0.9003 g2_cw_raw_median=0.9898 eligible=16/96 flux_floor_excluded=80 favorable_rows=4
```

## Literature ceiling
In the reviewed literature set of this repository (six papers, ../_goal/paper_digests.md) and the anchors ledger, no electrically driven III-V single-dot g2(0) at 300 K is reported; the best electrical result in that set is Reischle et al. 2008 at 80 K: g2(0) = 0.43 raw, 0.03 after background correction. The best reported 300 K single-dot values are optically pumped: g2(0) ~ 0.5-0.57 (Laferriere et al. 2023, InAsP/InP nanowire dot, g2(0) = 0.57 at 300 K). This sweep's pooled pulsed g2_min=0.4305, g2_median=0.7076 is reported against that electrical-vs-optical literature picture, not as a claim of an existing electrical 300 K result to exceed.

## Conditional on the 300 K linewidth
The 300 K single-dot linewidth gamma300 is not a platform ceiling -- it is an unmeasured [E]-class quantity (see Assumptions below) whose actual value decides whether the headline gate passes. The table below reports, per card, the LARGEST gamma300 in this sweep's grid at which any eligible row (any lever/delta_xx/irf combination) still clears pulsed intrinsic g2(0) < 0.5; a smaller real-world gamma300 than this value keeps the corresponding card passing, a larger one does not.

`gamma300_pass_max` (pooled, both cards): 6 meV

| card | gamma300_pass_max (meV) | delta_xx | irf | lever values | assumptions |
|---|---:|---:|---:|---|---|
| edge-inp-gaasp-design | n/a (no headline-passing row in this grid) | | | | |
| edge-inp-gainp-design | 6 | 7 meV | 50 ps | NA=0.8, R_back=0.95, L_um=250 | dot.gamma300; drive.b_res; drive.diode.tau_pulse_ns; drive.duty (rep-rate-derived); emission.lambda_nm; ret.system.barrier.x_al |

## Grid
Grid complete: True
- `dot.delta_xx` in [4.0, 7.0] meV: sampled at 4, 7
- `dot.gamma300` in [6.0, 20.0] meV: sampled at 6, 20
- `irf_ps` in [50.0, 200.0] ps: sampled at 50, 200

Collection-lever axes (per-card, council review round 4 item 1):
- `edge-inp-gaasp-design` `emission.NA`: sampled at 0.5, 0.8, 0.75
- `edge-inp-gaasp-design` `emission.R_back`: sampled at 0, 0.95
- `edge-inp-gaasp-design` `emission.L_um`: sampled at 250, 500
- `edge-inp-gaasp-design` lever combinations evaluated: 12
- `edge-inp-gainp-design` `emission.NA`: sampled at 0.5, 0.8, 0.75
- `edge-inp-gainp-design` `emission.R_back`: sampled at 0, 0.95
- `edge-inp-gainp-design` `emission.L_um`: sampled at 250, 500
- `edge-inp-gainp-design` lever combinations evaluated: 12

## Coverage
- eligible coverage: 16/192 = 0.083
- pulsed collected-flux eligibility floor [A]: 1000 photons/s; rows excluded by this floor: 176
- diagnostic pooled g2 (below flux floor, not measurable): pulsed 0.4305 / 0.7409; g2_cw0 0.3915 / 0.8133; g2_cw0_raw 0.9003 / 0.9958
- maximum collected pulsed flux: 1111 photons/s; shortfall factor (floor/flux_max): 0.9001
- **headline coverage** (pulsed intrinsic g2(0) < 0.5, eligible rows -- the contract's PASS metric): 4/192 = 0.021
- secondary coverage, CW intrinsic g2_cw0 < 0.5 (eligible rows, diagnostic only, does not gate PASS): 4/192 = 0.021
- secondary coverage, CW IRF-convolved g2_cw0_raw < 0.5 (eligible rows, diagnostic only, does not gate PASS): 0/192 = 0.000

## Per-card statistics
| card | role | diagnostic g2_pulsed min/median | diagnostic g2_cw0 min/median | diagnostic g2_cw0_raw min/median | eligible | flux-floor excluded | headline rows |
|---|---|---|---|---|---|---|---|
| edge-inp-gaasp-design | primary | 0.5306 / 0.7856 | 0.4017 / 0.8252 | 0.9784 / 0.9983 | 0/96 | 96 | 0 |
| edge-inp-gainp-design | fallback | 0.4305 / 0.7076 | 0.3915 / 0.8086 | 0.9003 / 0.9898 | 16/96 | 80 | 4 |

## Best diagnostic-g2 row and brightness decomposition
The dominant brightness limiter at the favourable diagnostic corner is beta (waveguide coupling / spontaneous-emission factor); the multiplicative chain that reproduces the reported collected pulsed flux is:

| factor | value |
|---|---:|
| loading = 1 - e^-mu (mu=1) | 0.632121 |
| t_X (spectral transmission) | 0.5 |
| S (confinement retention) | 0.00461745 |
| eta_total (edge out-coupling: waveguide coupling x front/back facet split x facet transmission x propagation x NA, all in one factor) | 0.00951594 |
| rep rate (Hz) | 8e+07 |
| **product x rep rate** | 1111 photons/s |
| reported collected_flux_pulsed_s | 1111 photons/s |
| self-check: relative difference | 0.0000% (PASS, within 1%) |

`beta`/`front`/`facet`/`propagation`/`NA` are shown below for diagnosis only -- they are already folded into `eta_total` above exactly once each (their product reproduces eta_total) and must NOT also be multiplied into the flux self-check:

| component (already inside eta_total) | value |
|---|---:|
| beta (waveguide coupling / spontaneous-emission factor) | 0.0292412 |
| front (front/back facet split: 0.5 with no R_back, else T_facet/(T_facet+(1-R_back)); back-solved -- device.py folds it into eta_total but never exposes it on its own) | 0.935012 |
| facet (Fresnel transmission) | 0.719371 |
| propagation | 0.882497 |
| NA (numerical aperture) | 0.548244 |

## Diagnostic g2 landscape
| corner | pulsed g2 min | pooled diagnostic median | lever values | assumptions |
|---|---:|---:|---|---|
| edge-inp-gainp-design (7 meV, gamma300=6 meV, irf=50 ps) | 0.430515 | 0.740852 | NA=0.8, R_back=0.95, L_um=250 | dot.gamma300; drive.b_res; drive.diode.tau_pulse_ns; drive.duty (rep-rate-derived); emission.lambda_nm; ret.system.barrier.x_al |

## Assumptions required by any headline-passing corner
- `dot.gamma300`
- `drive.b_res`
- `drive.diode.tau_pulse_ns`
- `drive.duty (rep-rate-derived)`
- `emission.lambda_nm`
- `ret.system.barrier.x_al`

## Assumptions in plain words
**1. The collection window assumption makes flux look independent of the 300 K linewidth.** Every row's spectral filter is set automatically to exactly match the exciton line's own width at the operating temperature (device.py's `filter.auto_w` convention; `filter.auto_w_scale` can widen or narrow it, 1.0 = as-is here). A plain mathematical consequence of matching the window to the line's own width is that exactly HALF of the X-line's photons fall inside the window in every single row of this sweep (t_X = 0.5, column `t_x_pulsed`) -- no matter how wide or narrow the real (currently unmeasured) 300 K linewidth gamma300 turns out to be. That is why the reported collected photon flux does not move with gamma300 in this sweep: it is not that the real device would be insensitive to linewidth, it is that this filter-window convention cancels that sensitivity out by construction.

**2. The background-light assumption is borrowed from a different, colder device.** Every row carries a constant background term (`drive.b_res`) that is not measured on this platform: it is transferred from Reischle et al., Appl. Phys. Lett. 92, 233113 (2008), whose 80 K electrically driven single-photon source had about 88% real signal and 12% background light (signal fraction rho ~ 0.88). This sweep assumes the same 12% background fraction still applies at 300 K, on a different material system (InP/GaAsP or InP/GaInP edge emitters) than the one actually measured. No 300 K electrical background measurement exists for either card's platform.

**3. Why pulsed drive, not continuous-wave (CW) drive, is required.** At the best diagnostic operating point in this sweep, the intrinsic CW g2(0) is 0.391 (that alone would already satisfy the g2 < 0.5 single-photon criterion), but once a realistic single-photon detector's finite timing resolution (instrument response function, IRF) is folded in, the measured raw CW g2(0) rises to 0.9 -- ABOVE the 0.5 threshold. In plain terms: the antibunching dip this device produces under continuous drive is narrower in time than a real detector can resolve, so a CW measurement alone could never demonstrate single-photon emission on this platform. Pulsed (gated) operation sidesteps the detector's timing resolution and is therefore required, not optional.

## Evidence gate
- evidence_complete: False
- self_test_passed (self-test): True
- missing/incomplete evidence claims:
  - `reischle2008_g2_80K`: only 0 distinct verified primary source(s), need 2
  - `hkust_inp_gaasp_wavelength`: only 1 distinct verified primary source(s), need 2
  - `residual_background_80K`: only 1 distinct verified primary source(s), need 2

## Fail reasons
- `evidence_incomplete`

## Sources
- docs/rt_edge_contract.md (acceptance gates, evidence status)
- verify/data/rt_edge_anchors.yaml (literature evidence ledger)
- cards/edge-inp-gaasp-design.yaml, cards/edge-inp-gainp-design.yaml (design.provenance.sources for every scalar)
