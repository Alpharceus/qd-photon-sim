# RT edge-emitter acceptance sweep verdict

Generated 2026-09-05T21:20:46.482341+00:00; contract: docs/rt_edge_contract.md.

```
VERDICT: FAIL g2_min=0.7723 g2_median=0.9509 median_pass=false coverage=0 eligible=54/54 evidence=incomplete conditional=false headline_coverage=0/54 cw_raw_coverage=0/54
```
```
CARD: edge-inp-gaasp-design role=primary g2_pulsed_min=0.9763 g2_pulsed_median=0.9941 g2_cw_raw_min=0.9978 g2_cw_raw_median=0.9997 eligible=27/27 favorable_rows=0
```
```
CARD: edge-inp-gainp-design role=fallback g2_pulsed_min=0.7723 g2_pulsed_median=0.8747 g2_cw_raw_min=0.9628 g2_cw_raw_median=0.9909 eligible=27/27 favorable_rows=0
```

## Literature ceiling
No electrically driven III-V single quantum dot g2(0) at 300 K has been published (docs/rt_edge_contract.md Evidence section; verify/data/rt_edge_anchors.yaml) -- there is no existing electrical 300 K baseline for this sweep to exceed. The best reported electrically driven single-dot result at any temperature is Reischle et al. 2008 at 80 K: g2(0) = 0.43 raw, 0.03 after background correction. The best reported 300 K single-dot values are optically pumped: g2(0) ~ 0.5-0.57 (Laferriere et al. 2023, InAsP/InP nanowire dot, g2(0) = 0.57 at 300 K). This sweep's pooled pulsed g2_min=0.7723, g2_median=0.9509 is reported against that electrical-vs-optical literature picture, not as a claim of an existing electrical 300 K result to exceed.

## Grid
Grid complete: True
- `dot.delta_xx` in [4.0, 7.0] meV: sampled at 4, 5.5, 7
- `dot.gamma300` in [6.0, 20.0] meV: sampled at 6, 13, 20
- `irf_ps` in [50.0, 200.0] ps: sampled at 50, 125, 200

## Coverage
- eligible coverage: 54/54 = 1.000
- **headline coverage** (pulsed intrinsic g2(0) < 0.5, eligible rows -- the contract's PASS metric): 0/54 = 0.000
- secondary coverage, CW intrinsic g2_cw0 < 0.5 (diagnostic only, does not gate PASS): 0/54 = 0.000
- secondary coverage, CW IRF-convolved g2_cw0_raw < 0.5 (diagnostic only, does not gate PASS): 0/54 = 0.000

## Per-card statistics
| card | role | g2_pulsed min/median | g2_cw0 min/median | g2_cw0_raw min/median | eligible | headline rows |
|---|---|---|---|---|---|---|
| edge-inp-gaasp-design | primary | 0.9763 / 0.9941 | 0.94 / 0.9823 | 0.9978 / 0.9997 | 27/27 | 0 |
| edge-inp-gainp-design | fallback | 0.7723 / 0.8747 | 0.7762 / 0.8773 | 0.9628 / 0.9909 | 27/27 | 0 |

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
- `no_headline_pass`

## Sources
- docs/rt_edge_contract.md (acceptance gates, evidence status)
- verify/data/rt_edge_anchors.yaml (literature evidence ledger)
- cards/edge-inp-gaasp-design.yaml, cards/edge-inp-gainp-design.yaml (design.provenance.sources for every scalar)
