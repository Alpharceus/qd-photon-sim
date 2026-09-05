# RT edge-emitter acceptance sweep verdict

Generated 2026-09-05T20:24:08.189904+00:00; contract: docs/rt_edge_contract.md.

```
VERDICT: FAIL g2_min=0.3648 g2_median=0.773 median_pass=false coverage=0 eligible=54/54 evidence=incomplete conditional=false
```
```
CARD: edge-inp-gaasp-design role=primary g2_pulsed_min=0.4117 g2_pulsed_median=0.8555 g2_cw_raw_min=0.9784 g2_cw_raw_median=0.9992 eligible=27/27 favorable_rows=0
```
```
CARD: edge-inp-gainp-design role=fallback g2_pulsed_min=0.3648 g2_pulsed_median=0.7325 g2_cw_raw_min=0.8907 g2_cw_raw_median=0.9926 eligible=27/27 favorable_rows=0
```

## Literature ceiling
No published III-V quantum dot has demonstrated g2(0) < 0.5 at 300 K under electrical driving; the best reported room-temperature values sit around 0.5-0.57 (docs/rt_edge_contract.md Evidence section; verify/data/rt_edge_anchors.yaml). This sweep's pooled pulsed g2_min=0.3648, g2_median=0.773 is reported against that ceiling, not as a claim of having exceeded it.

## Grid
Grid complete: True
- `dot.delta_xx` in [4.0, 7.0] meV: sampled at 4, 5.5, 7
- `dot.gamma300` in [6.0, 20.0] meV: sampled at 6, 13, 20
- `irf_ps` in [50.0, 200.0] ps: sampled at 50, 125, 200

## Coverage
- eligible coverage: 54/54 = 1.000
- secondary-gate coverage (CW intrinsic AND raw < 0.5, eligible): 0.000 combined / 0.000 secondary-only
- combined passing coverage (all three gates, eligible): 0/54 = 0.000

## Per-card statistics
| card | role | g2_pulsed min/median | g2_cw0 min/median | g2_cw0_raw min/median | eligible | favorable rows |
|---|---|---|---|---|---|---|
| edge-inp-gaasp-design | primary | 0.4117 / 0.8555 | 0.4026 / 0.9443 | 0.9784 / 0.9992 | 27/27 | 0 |
| edge-inp-gainp-design | fallback | 0.3648 / 0.7325 | 0.3413 / 0.905 | 0.8907 / 0.9926 | 27/27 | 0 |

## Assumptions required by any favorable corner
- (no row satisfies all three same-row gates; see fail_reasons below)

## Evidence gate
- evidence_complete: False
- hallucination_tests_passed (self-test): True
- missing/incomplete evidence claims:
  - `reischle2008_g2_80K`: only 0 distinct verified primary source(s), need 2
  - `hkust_inp_gaasp_wavelength`: only 0 distinct verified primary source(s), need 2

## Fail reasons
- `evidence_incomplete`
- `no_favorable_corner_metrics`

## Sources
- docs/rt_edge_contract.md (acceptance gates, evidence status)
- verify/data/rt_edge_anchors.yaml (literature evidence ledger)
- cards/edge-inp-gaasp-design.yaml, cards/edge-inp-gainp-design.yaml (design.provenance.sources for every scalar)
