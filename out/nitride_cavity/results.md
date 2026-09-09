# Nitride cavity sweep results

Assumed planar model, evaluator outputs only; no held-out prediction.

VERDICT: FAIL regime=rectangular model=finite electrical pulse complete=True g2_min=none g2_median_eligible=none diag_g2_min=0.9769926994149587 diag_g2_flux_max=4.226009050399168e-18 flux_max=92.48444325385057 flux_margin=0.09248444325385057 eligible_fraction=0.0 eligible=0 flux_floor_excluded=1620 coverage_over_eligible=none headline_coverage=1620/1620 idealized_pass_count=0 hardware_infeasible_count=0 hardware_pass_count=0 evidence=incomplete conditional=False T_pass_min=none

VERDICT: FAIL regime=deterministic_pair model=idealized deterministic SET complete=True g2_min=0.0 g2_median_eligible=0.0 diag_g2_min=0.0 diag_g2_flux_max=5600.493440683151 flux_max=5600.493440683151 flux_margin=5.600493440683151 eligible_fraction=0.016666666666666666 eligible=27 flux_floor_excluded=1593 coverage_over_eligible=0.0 headline_coverage=1620/1620 idealized_pass_count=27 hardware_infeasible_count=1620 hardware_pass_count=0 evidence=incomplete conditional=False T_pass_min=none

## Deshpande comparison
Measured g2=0.29 [V abstract-only; Deshpande et al., APL 105, 141109 (2014), DOI 10.1063/1.4897640]. Evaluated card g2=0.9976537308395578, flux=0.018999124213812685/s; count rate unavailable. CONDITIONS INCOMPLETE.
| transfer | value |
|---|---|
| geometry/x | 2 nm / 12.5 nm [A], x=0.40 [V abstract-only] |
| rate/waveform | 200 MHz reported maximum [V], waveform/current [A] |
| field/cavity | planar QCSE, Q/V/outcoupling [A] |

## Sensitivity
Full mode contains all one-at-a-time contract axes at 230/300 K; quick contains named diagnostics. screening_fraction 0 through 1 SET rows are in sensitivities.csv.

## Limitations
Q/V realization, oscillator/QCSE uncertainty, omitted field-assisted tunnelling, optimistic nonradiative/background assumptions, SET pair delivery and hardware screen limit this model. Idealized optical pass is not hardware demonstration.
