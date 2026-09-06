# RT edge-emitter acceptance sweep verdict

Generated 2026-09-06T04:16:57.223992+00:00; contract: docs/rt_edge_contract.md.

```
VERDICT: FAIL g2_min=nan g2_median=nan diag_g2_min=0.2154 diag_g2_median=0.8683 flux_max=43.12 flux_shortfall=23.19 median_pass=false coverage=0 eligible=0/54 flux_floor_excluded=54 evidence=incomplete conditional=false headline_coverage=0/54 cw_raw_coverage=0/54
```
```
CARD: edge-inp-gaasp-design role=primary g2_pulsed_min=nan g2_pulsed_median=nan diagnostic (below flux floor, not measurable) diag_g2_min=0.2437 diag_g2_median=0.8735 diag_g2_cw0_min=0.2413 diag_g2_cw0_median=0.9126 diag_g2_cw0_raw_min=0.9726 diag_g2_cw0_raw_median=0.9987 g2_cw_raw_min=nan g2_cw_raw_median=nan eligible=0/27 flux_floor_excluded=27 favorable_rows=0
```
```
CARD: edge-inp-gainp-design role=fallback g2_pulsed_min=nan g2_pulsed_median=nan diagnostic (below flux floor, not measurable) diag_g2_min=0.2154 diag_g2_median=0.8631 diag_g2_cw0_min=0.2133 diag_g2_cw0_median=0.903 diag_g2_cw0_raw_min=0.8694 diag_g2_cw0_raw_median=0.9933 g2_cw_raw_min=nan g2_cw_raw_median=nan eligible=0/27 flux_floor_excluded=27 favorable_rows=0
```

## Literature ceiling
In the reviewed literature set of this repository (six papers, ../_goal/paper_digests.md) and the anchors ledger, no electrically driven III-V single-dot g2(0) at 300 K is reported; the best electrical result in that set is Reischle et al. 2008 at 80 K: g2(0) = 0.43 raw, 0.03 after background correction. The best reported 300 K single-dot values are optically pumped: g2(0) ~ 0.5-0.57 (Laferriere et al. 2023, InAsP/InP nanowire dot, g2(0) = 0.57 at 300 K). This sweep's pooled pulsed g2_min=nan, g2_median=nan is reported against that electrical-vs-optical literature picture, not as a claim of an existing electrical 300 K result to exceed.

## Grid
Grid complete: True
- `dot.delta_xx` in [4.0, 7.0] meV: sampled at 4, 5.5, 7
- `dot.gamma300` in [6.0, 20.0] meV: sampled at 6, 13, 20
- `irf_ps` in [50.0, 200.0] ps: sampled at 50, 125, 200

## Coverage
- eligible coverage: 0/54 = 0.000
- pulsed collected-flux eligibility floor [A]: 1000 photons/s; rows excluded by this floor: 54
- diagnostic pooled g2 (below flux floor, not measurable): pulsed 0.2154 / 0.8683; g2_cw0 0.2133 / 0.9078; g2_cw0_raw 0.8694 / 0.9973
- maximum collected pulsed flux: 43.12 photons/s; shortfall factor (floor/flux_max): 23.19
- **headline coverage** (pulsed intrinsic g2(0) < 0.5, eligible rows -- the contract's PASS metric): 0/54 = 0.000
- secondary coverage, CW intrinsic g2_cw0 < 0.5 (diagnostic only, does not gate PASS): 0/54 = 0.000
- secondary coverage, CW IRF-convolved g2_cw0_raw < 0.5 (diagnostic only, does not gate PASS): 0/54 = 0.000

## Per-card statistics
| card | role | diagnostic g2_pulsed min/median | diagnostic g2_cw0 min/median | diagnostic g2_cw0_raw min/median | eligible | flux-floor excluded | headline rows |
|---|---|---|---|---|---|---|---|
| edge-inp-gaasp-design | primary | 0.2437 / 0.8735 | 0.2413 / 0.9126 | 0.9726 / 0.9987 | 0/27 | 27 | 0 |
| edge-inp-gainp-design | fallback | 0.2154 / 0.8631 | 0.2133 / 0.903 | 0.8694 / 0.9933 | 0/27 | 27 | 0 |

## Why no row is eligible
Every row is below the 1000 photons/s collected-flux floor; the grid maximum is only 43.12 photons/s (floor/maximum = 23.19).
The dominant brightness limiter at the favourable diagnostic corner is S; the factor decomposition is:

| factor | value |
|---|---:|
| mu | 0.100002 |
| S | 0.00461745 |
| t_X | 0.5 |
| beta | 0.0292412 |
| facet | 0.719371 |
| propagation | 0.778801 |
| NA | 0.299495 |
| eta_total | 0.0024532 |
| rep rate | 8e+07 |
| duty | 0.008 |

## Diagnostic g2 landscape
| corner | pulsed g2 min | pooled diagnostic median | assumptions |
|---|---:|---:|---|
| edge-inp-gainp-design (7 meV, gamma300=6 meV, irf=50 ps) | 0.2154 | 0.868329 | dot.gamma300; drive.diode.tau_pulse_ns; drive.duty (rep-rate-derived); emission.lambda_nm; ret.system.barrier.x_al |

## Assumptions required by any headline-passing corner
- (no row satisfies the headline gate; see fail_reasons below)

## Evidence gate
- evidence_complete: False
- hallucination_tests_passed (self-test): True
- missing/incomplete evidence claims:
  - `reischle2008_g2_80K`: only 0 distinct verified primary source(s), need 2
  - `hkust_inp_gaasp_wavelength`: only 1 distinct verified primary source(s), need 2

## Fail reasons
- `evidence_incomplete`
- `no_eligible_rows`

## Sources
- docs/rt_edge_contract.md (acceptance gates, evidence status)
- verify/data/rt_edge_anchors.yaml (literature evidence ledger)
- cards/edge-inp-gaasp-design.yaml, cards/edge-inp-gainp-design.yaml (design.provenance.sources for every scalar)
