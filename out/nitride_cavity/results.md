# Nitride cavity sweep results

Assumed planar model, evaluator outputs only; no held-out prediction.

VERDICT: FAIL regime=rectangular model=finite electrical pulse complete=True g2_min=none g2_median_eligible=none diag_g2_min=0.9809857019958337 diag_g2_flux_max=4.226009050399168e-18 flux_max=92.48444325385057 flux_margin=0.09248444325385057 eligible_fraction=0.0 eligible=0 flux_floor_excluded=1620 coverage_over_eligible=none headline_coverage=1620/1620 idealized_pass_count=0 idealized_status=no_idealized_pass hardware_infeasible_count=0 hardware_pass_count=0 evidence=incomplete conditional=False T_pass_min=none

VERDICT: FAIL regime=deterministic_pair model=idealized deterministic SET complete=True g2_min=0.17355371900826433 g2_median_eligible=0.17355371900826433 diag_g2_min=0.17355371900826433 diag_g2_flux_max=5543.33605560159 flux_max=5600.493440683151 flux_margin=5.600493440683151 eligible_fraction=0.016666666666666666 eligible=27 flux_floor_excluded=1593 coverage_over_eligible=0.0 headline_coverage=1620/1620 idealized_pass_count=27 idealized_status=pass_hardware_infeasible hardware_infeasible_count=1620 hardware_pass_count=0 evidence=incomplete conditional=False T_pass_min=none

## Bounds (unscreened lower bound / screening_fraction=1 upper bound)
Headline rows use the card default screening_fraction=0.0 -- an explicit CONSERVATIVE LOWER BOUND (no polarization-field screening). screening_fraction=1.0 rows below are the UPPER BOUND, evaluated the same way (row_kind='bound' in sensitivities.csv), never hand-entered.

| regime | T_hs K | screening_fraction | g2_op | flux /s | S_X | E_C/kT | hardware_feasible |
|---|---|---|---|---|---|---|---|
| deterministic_pair | 230 | 0 | 0.173554 | 0.829668 | 1.646e-06 | 1.53 | False |
| deterministic_pair | 230 | 0.5 | 0.173554 | 328.736 | 0.0005512 | 1.53 | False |
| deterministic_pair | 230 | 1 | 0.173554 | 400172 | 0.6234 | 1.53 | False |
| deterministic_pair | 300 | 0 | 0.173554 | 0.137422 | 5.703e-07 | 1.173 | False |
| deterministic_pair | 300 | 0.5 | 0.173554 | 30.903 | 0.0001073 | 1.173 | False |
| deterministic_pair | 300 | 1 | 0.173554 | 25642.8 | 0.08233 | 1.173 | False |
| rectangular | 230 | 0 | 0.998901 | 0.0256475 | 1.646e-06 |  | n/a |
| rectangular | 230 | 0.5 | 0.987036 | 0.00043571 | 0.0005512 |  | n/a |
| rectangular | 230 | 1 | 0.588879 | 6.5708e-06 | 0.6234 |  | n/a |
| rectangular | 300 | 0 | 0.999679 | 0.00424714 | 5.703e-07 |  | n/a |
| rectangular | 300 | 0.5 | 0.997844 | 8.44093e-05 | 0.0001073 |  | n/a |
| rectangular | 300 | 1 | 0.895883 | 1.18216e-05 | 0.08233 |  | n/a |

Note: with drive.b_res wired into this branch (fix round 2), SET g2_op is NOT exactly 0: under ideal one-pair loading the exact-one-pair counting result is 0 by construction (cnt_g2=0), so g2_op=1-rho^2 depends only on rho=signal/(signal+background); since bg_counts includes b_res*(collected X counts) and one-pair loading has negligible XX, signal~=collected X counts, so rho asymptotes to 1/(1+b_res) and g2_op asymptotes to 1-(1/(1+b_res))^2 nearly independent of absolute flux (b_res=0.1 on the shipped cards -> g2_op~0.174, matching the bounds table below across screening/T/Q). This is a STRUCTURAL floor set by the assumed residual background channel, not a demonstrated device number.

Cavity re-tuning disclosure: every headline row re-tracks the cavity resonance to T_track=T_hs (the SAME dot's own E_X at that operating point) -- headline results assume a cavity re-tuned per dot/temperature, not a single fixed cavity swept across all conditions. The fixed-anchor (T_track=300 K) sensitivity set below holds the cavity fixed while T_hs varies, to show the resulting detuning as a labelled, separate effect.

### Fixed-anchor (T_track=300 K) sensitivity

| regime | T_hs K | T_track K | detuning_meV | g2_op | flux /s |
|---|---|---|---|---|---|
| deterministic_pair | 230 | 300 | 22.21 | 0.173554 | 0.00125519 |
| deterministic_pair | 250 | 300 | 16.31 | 0.173554 | 0.00357734 |
| deterministic_pair | 273 | 300 | 9.07 | 0.173554 | 0.136923 |
| deterministic_pair | 300 | 300 | -2.159e-07 | 0.173554 | 0.137422 |
| rectangular | 230 | 300 | 22.21 | 0.999138 | 3.88017e-05 |
| rectangular | 250 | 300 | 16.31 | 0.999236 | 0.000110584 |
| rectangular | 273 | 300 | 9.07 | 0.999296 | 0.00423232 |
| rectangular | 300 | 300 | -2.159e-07 | 0.999679 | 0.00424714 |

### SET island-radius sensitivity (classical charging-energy screen only; not demonstrated feasible manufacture)

| island radius nm | T_hs K | E_C meV | E_C/kT | allowed max radius nm | R_T/R_Q | hardware_feasible | f_max Hz |
|---|---|---|---|---|---|---|---|
| 0.5 | 230 | 303.2 | 15.3 | 0.7648 | 38.74 | True | 1.892e+11 |
| 0.5 | 250 | 303.2 | 14.07 | 0.7036 | 38.74 | True | 1.892e+11 |
| 0.5 | 273 | 303.2 | 12.89 | 0.6443 | 38.74 | True | 1.892e+11 |
| 0.5 | 300 | 303.2 | 11.73 | 0.5863 | 38.74 | True | 1.892e+11 |
| 1 | 230 | 151.6 | 7.648 | 0.7648 | 38.74 | False | 9.461e+10 |
| 1 | 250 | 151.6 | 7.036 | 0.7036 | 38.74 | False | 9.461e+10 |
| 1 | 273 | 151.6 | 6.443 | 0.6443 | 38.74 | False | 9.461e+10 |
| 1 | 300 | 151.6 | 5.863 | 0.5863 | 38.74 | False | 9.461e+10 |
| 5 | 230 | 30.32 | 1.53 | 0.7648 | 38.74 | False | 1.892e+10 |
| 5 | 250 | 30.32 | 1.407 | 0.7036 | 38.74 | False | 1.892e+10 |
| 5 | 273 | 30.32 | 1.289 | 0.6443 | 38.74 | False | 1.892e+10 |
| 5 | 300 | 30.32 | 1.173 | 0.5863 | 38.74 | False | 1.892e+10 |

## Deshpande comparison
Measured g2=0.29 [V abstract-only; Deshpande et al., APL 105, 141109 (2014), DOI 10.1063/1.4897640]. Evaluated card g2=0.9980609340138185, flux=0.018999124213812685/s; count rate unavailable. CONDITIONS INCOMPLETE.
| transfer | value |
|---|---|
| geometry/x | 2 nm / 12.5 nm [A], x=0.40 [V abstract-only] |
| rate/waveform | 200 MHz reported maximum [V], waveform/current [A] |
| field/cavity | planar QCSE, Q/V/outcoupling [A] |

## Sensitivity
Full mode contains all one-at-a-time contract axes at 230/300 K, both regimes (quick contains the named diagnostic subset plus Q=167 and purcell_enabled=False as two independent one-at-a-time rows). Values below are the actual evaluated g2/flux for each axis value.

| axis | regime | T_hs K | value | g2_op | flux /s |
|---|---|---|---|---|---|
| Q | deterministic_pair | 230 | 167.0 | 0.173554 | 3.30937 |
| Q | deterministic_pair | 300 | 167.0 | 0.173554 | 0.958533 |
| Q | rectangular | 230 | 167.0 | 0.999173 | 0.102303 |
| Q | rectangular | 300 | 167.0 | 0.999782 | 0.0296243 |
| b_res | deterministic_pair | 230 | 0.0 | 0 | 0.829668 |
| b_res | deterministic_pair | 230 | 0.1 | 0.173554 | 0.829668 |
| b_res | deterministic_pair | 230 | 1.0 | 0.75 | 0.829668 |
| b_res | deterministic_pair | 300 | 0.0 | 0 | 0.137422 |
| b_res | deterministic_pair | 300 | 0.1 | 0.173554 | 0.137422 |
| b_res | deterministic_pair | 300 | 1.0 | 0.75 | 0.137422 |
| b_res | rectangular | 230 | 0.0 | 0.99867 | 0.0256475 |
| b_res | rectangular | 230 | 0.1 | 0.998901 | 0.0256475 |
| b_res | rectangular | 230 | 1.0 | 0.999668 | 0.0256475 |
| b_res | rectangular | 300 | 0.0 | 0.999612 | 0.00424714 |
| b_res | rectangular | 300 | 0.1 | 0.999679 | 0.00424714 |
| b_res | rectangular | 300 | 1.0 | 0.999903 | 0.00424714 |
| background_tau_ns | deterministic_pair | 230 | 0.0 | 0.173554 | 0.829668 |
| background_tau_ns | deterministic_pair | 230 | 1.0 | 0.173554 | 0.829668 |
| background_tau_ns | deterministic_pair | 300 | 0.0 | 0.173554 | 0.137422 |
| background_tau_ns | deterministic_pair | 300 | 1.0 | 0.173554 | 0.137422 |
| background_tau_ns | rectangular | 230 | 0.0 | 0.998901 | 0.0256475 |
| background_tau_ns | rectangular | 230 | 1.0 | 0.998901 | 0.0256475 |
| background_tau_ns | rectangular | 300 | 0.0 | 0.999679 | 0.00424714 |
| background_tau_ns | rectangular | 300 | 1.0 | 0.999679 | 0.00424714 |
| delta_xx | deterministic_pair | 230 | -10.0 | 0.173554 | 0.829668 |
| delta_xx | deterministic_pair | 230 | 10.0 | 0.173554 | 0.829668 |
| delta_xx | deterministic_pair | 230 | 16.0 | 0.173554 | 0.829668 |
| delta_xx | deterministic_pair | 300 | -10.0 | 0.173554 | 0.137422 |
| delta_xx | deterministic_pair | 300 | 10.0 | 0.173554 | 0.137422 |
| delta_xx | deterministic_pair | 300 | 16.0 | 0.173554 | 0.137422 |
| delta_xx | rectangular | 230 | -10.0 | 0.998901 | 0.0256475 |
| delta_xx | rectangular | 230 | 10.0 | 0.998901 | 0.0256475 |
| delta_xx | rectangular | 230 | 16.0 | 0.998699 | 0.0256474 |
| delta_xx | rectangular | 300 | -10.0 | 0.999679 | 0.00424714 |
| delta_xx | rectangular | 300 | 10.0 | 0.999679 | 0.00424714 |
| delta_xx | rectangular | 300 | 16.0 | 0.999455 | 0.00424712 |
| detuning_offset_meV | deterministic_pair | 230 | 0.0 | 0.173554 | 0.829668 |
| detuning_offset_meV | deterministic_pair | 230 | 10.0 | 0.173554 | 0.0251003 |
| detuning_offset_meV | deterministic_pair | 300 | 0.0 | 0.173554 | 0.137422 |
| detuning_offset_meV | deterministic_pair | 300 | 10.0 | 0.173554 | 0.086737 |
| detuning_offset_meV | rectangular | 230 | 0.0 | 0.998901 | 0.0256475 |
| detuning_offset_meV | rectangular | 230 | 10.0 | 0.998893 | 0.000775923 |
| detuning_offset_meV | rectangular | 300 | 0.0 | 0.999679 | 0.00424714 |
| detuning_offset_meV | rectangular | 300 | 10.0 | 0.999495 | 0.00268066 |
| eta_out | deterministic_pair | 230 | 0.01 | 0.173554 | 0.0829668 |
| eta_out | deterministic_pair | 230 | 0.1 | 0.173554 | 0.829668 |
| eta_out | deterministic_pair | 230 | 0.5 | 0.173554 | 4.14834 |
| eta_out | deterministic_pair | 300 | 0.01 | 0.173554 | 0.0137422 |
| eta_out | deterministic_pair | 300 | 0.1 | 0.173554 | 0.137422 |
| eta_out | deterministic_pair | 300 | 0.5 | 0.173554 | 0.687111 |
| eta_out | rectangular | 230 | 0.01 | 0.998901 | 0.00256475 |
| eta_out | rectangular | 230 | 0.1 | 0.998901 | 0.0256475 |
| eta_out | rectangular | 230 | 0.5 | 0.998901 | 0.128237 |
| eta_out | rectangular | 300 | 0.01 | 0.999679 | 0.000424714 |
| eta_out | rectangular | 300 | 0.1 | 0.999679 | 0.00424714 |
| eta_out | rectangular | 300 | 0.5 | 0.999679 | 0.0212357 |
| gamma300 | deterministic_pair | 230 | 25.0 | 0.173554 | 1.22686 |
| gamma300 | deterministic_pair | 230 | 35.0 | 0.173554 | 0.829668 |
| gamma300 | deterministic_pair | 230 | 55.0 | 0.173554 | 0.462512 |
| gamma300 | deterministic_pair | 300 | 25.0 | 0.173554 | 0.240419 |
| gamma300 | deterministic_pair | 300 | 35.0 | 0.173554 | 0.137422 |
| gamma300 | deterministic_pair | 300 | 55.0 | 0.173554 | 0.0665605 |
| gamma300 | rectangular | 230 | 25.0 | 0.998778 | 0.0379257 |
| gamma300 | rectangular | 230 | 35.0 | 0.998901 | 0.0256475 |
| gamma300 | rectangular | 230 | 55.0 | 0.999153 | 0.0142977 |
| gamma300 | rectangular | 300 | 25.0 | 0.999512 | 0.00743032 |
| gamma300 | rectangular | 300 | 35.0 | 0.999679 | 0.00424714 |
| gamma300 | rectangular | 300 | 55.0 | 0.999847 | 0.00205711 |
| k_nr_ns | deterministic_pair | 230 | 0.0 | 0.173554 | 0.829668 |
| k_nr_ns | deterministic_pair | 230 | 1.0 | 0.173554 | 0.829597 |
| k_nr_ns | deterministic_pair | 300 | 0.0 | 0.173554 | 0.137422 |
| k_nr_ns | deterministic_pair | 300 | 1.0 | 0.173554 | 0.137415 |
| k_nr_ns | rectangular | 230 | 0.0 | 0.998901 | 0.0256475 |
| k_nr_ns | rectangular | 230 | 1.0 | 0.998901 | 0.0256453 |
| k_nr_ns | rectangular | 300 | 0.0 | 0.999679 | 0.00424714 |
| k_nr_ns | rectangular | 300 | 1.0 | 0.999679 | 0.00424691 |
| mode_volume_norm | deterministic_pair | 230 | 1.0 | 0.173554 | 1.51718 |
| mode_volume_norm | deterministic_pair | 230 | 10.0 | 0.173554 | 0.279657 |
| mode_volume_norm | deterministic_pair | 230 | 2.0 | 0.173554 | 0.829668 |
| mode_volume_norm | deterministic_pair | 300 | 1.0 | 0.173554 | 0.23319 |
| mode_volume_norm | deterministic_pair | 300 | 10.0 | 0.173554 | 0.0608081 |
| mode_volume_norm | deterministic_pair | 300 | 2.0 | 0.173554 | 0.137422 |
| mode_volume_norm | rectangular | 230 | 1.0 | 0.998875 | 0.0469005 |
| mode_volume_norm | rectangular | 230 | 10.0 | 0.999015 | 0.00864503 |
| mode_volume_norm | rectangular | 230 | 2.0 | 0.998901 | 0.0256475 |
| mode_volume_norm | rectangular | 300 | 1.0 | 0.999659 | 0.00720691 |
| mode_volume_norm | rectangular | 300 | 10.0 | 0.99974 | 0.00187932 |
| mode_volume_norm | rectangular | 300 | 2.0 | 0.999679 | 0.00424714 |
| purcell_enabled | deterministic_pair | 230 | False | 0.173554 | 0.142154 |
| purcell_enabled | deterministic_pair | 230 | False | 0.173554 | 0.142154 |
| purcell_enabled | deterministic_pair | 230 | True | 0.173554 | 0.829668 |
| purcell_enabled | deterministic_pair | 300 | False | 0.173554 | 0.0416546 |
| purcell_enabled | deterministic_pair | 300 | False | 0.173554 | 0.0416546 |
| purcell_enabled | deterministic_pair | 300 | True | 0.173554 | 0.137422 |
| purcell_enabled | rectangular | 230 | False | 0.999181 | 0.00439441 |
| purcell_enabled | rectangular | 230 | False | 0.999181 | 0.00439441 |
| purcell_enabled | rectangular | 230 | True | 0.998901 | 0.0256475 |
| purcell_enabled | rectangular | 300 | False | 0.99979 | 0.00128737 |
| purcell_enabled | rectangular | 300 | False | 0.99979 | 0.00128737 |
| purcell_enabled | rectangular | 300 | True | 0.999679 | 0.00424714 |
| screening_fraction | deterministic_pair | 230 | 0.0 | 0.173554 | 0.829668 |
| screening_fraction | deterministic_pair | 230 | 0.5 | 0.173554 | 328.736 |
| screening_fraction | deterministic_pair | 230 | 1.0 | 0.173554 | 400172 |
| screening_fraction | deterministic_pair | 300 | 0.0 | 0.173554 | 0.137422 |
| screening_fraction | deterministic_pair | 300 | 0.5 | 0.173554 | 30.903 |
| screening_fraction | deterministic_pair | 300 | 1.0 | 0.173554 | 25642.8 |
| screening_fraction | rectangular | 230 | 0.0 | 0.998901 | 0.0256475 |
| screening_fraction | rectangular | 230 | 0.5 | 0.987036 | 0.00043571 |
| screening_fraction | rectangular | 230 | 1.0 | 0.588879 | 6.5708e-06 |
| screening_fraction | rectangular | 300 | 0.0 | 0.999679 | 0.00424714 |
| screening_fraction | rectangular | 300 | 0.5 | 0.997844 | 8.44093e-05 |
| screening_fraction | rectangular | 300 | 1.0 | 0.895883 | 1.18216e-05 |
| strain_fraction | deterministic_pair | 230 | 0.5 | 0.173554 | 2085.64 |
| strain_fraction | deterministic_pair | 230 | 1.0 | 0.173554 | 0.829668 |
| strain_fraction | deterministic_pair | 300 | 0.5 | 0.173554 | 134.485 |
| strain_fraction | deterministic_pair | 300 | 1.0 | 0.173554 | 0.137422 |
| strain_fraction | rectangular | 230 | 0.5 | 0.936075 | 0.387644 |
| strain_fraction | rectangular | 230 | 1.0 | 0.998901 | 0.0256475 |
| strain_fraction | rectangular | 300 | 0.5 | 0.992492 | 0.0161284 |
| strain_fraction | rectangular | 300 | 1.0 | 0.999679 | 0.00424714 |
| tau_cap_ps | deterministic_pair | 230 | 1.0 | 0.173554 | 0.0829669 |
| tau_cap_ps | deterministic_pair | 230 | 10.0 | 0.173554 | 0.829668 |
| tau_cap_ps | deterministic_pair | 230 | 100.0 | 0.173554 | 8.29656 |
| tau_cap_ps | deterministic_pair | 300 | 1.0 | 0.173554 | 0.0137422 |
| tau_cap_ps | deterministic_pair | 300 | 10.0 | 0.173554 | 0.137422 |
| tau_cap_ps | deterministic_pair | 300 | 100.0 | 0.173554 | 1.37422 |
| tau_cap_ps | rectangular | 230 | 1.0 | 0.99989 | 0.0025648 |
| tau_cap_ps | rectangular | 230 | 10.0 | 0.998901 | 0.0256475 |
| tau_cap_ps | rectangular | 230 | 100.0 | 0.989098 | 0.256424 |
| tau_cap_ps | rectangular | 300 | 1.0 | 0.999968 | 0.000424716 |
| tau_cap_ps | rectangular | 300 | 10.0 | 0.999679 | 0.00424714 |
| tau_cap_ps | rectangular | 300 | 100.0 | 0.996807 | 0.0424689 |
| tau_cap_scales_with_density | deterministic_pair | 230 | False | 0.173554 | 0.0829669 |
| tau_cap_scales_with_density | deterministic_pair | 230 | True | 0.173554 | 0.829668 |
| tau_cap_scales_with_density | deterministic_pair | 300 | False | 0.173554 | 0.0137422 |
| tau_cap_scales_with_density | deterministic_pair | 300 | True | 0.173554 | 0.137422 |
| tau_cap_scales_with_density | rectangular | 230 | False | 0.99989 | 0.00235495 |
| tau_cap_scales_with_density | rectangular | 230 | True | 0.998901 | 0.0235491 |
| tau_cap_scales_with_density | rectangular | 300 | False | 0.999968 | 0.000389967 |
| tau_cap_scales_with_density | rectangular | 300 | True | 0.999679 | 0.00389965 |
| tau_rad0_ns | deterministic_pair | 230 | 0.2 | 0.173554 | 4.14831 |
| tau_rad0_ns | deterministic_pair | 230 | 1.0 | 0.173554 | 0.829668 |
| tau_rad0_ns | deterministic_pair | 230 | 4.0 | 0.173554 | 0.207417 |
| tau_rad0_ns | deterministic_pair | 300 | 0.2 | 0.173554 | 0.68711 |
| tau_rad0_ns | deterministic_pair | 300 | 1.0 | 0.173554 | 0.137422 |
| tau_rad0_ns | deterministic_pair | 300 | 4.0 | 0.173554 | 0.0343556 |
| tau_rad0_ns | rectangular | 230 | 0.2 | 0.998901 | 0.128236 |
| tau_rad0_ns | rectangular | 230 | 1.0 | 0.998901 | 0.0256475 |
| tau_rad0_ns | rectangular | 230 | 4.0 | 0.998901 | 0.00641187 |
| tau_rad0_ns | rectangular | 300 | 0.2 | 0.999679 | 0.0212356 |
| tau_rad0_ns | rectangular | 300 | 1.0 | 0.999679 | 0.00424714 |
| tau_rad0_ns | rectangular | 300 | 4.0 | 0.999679 | 0.00106179 |

## Limitations
Q/V realization, oscillator/QCSE uncertainty, omitted field-assisted tunnelling, optimistic nonradiative/background assumptions, SET pair delivery and hardware screen limit this model. Idealized optical pass is not hardware demonstration.
Reservoir-vs-diode-SRH-layer disagreement: the background reservoir energy (fsim_core.device._nitride_reservoir_energy_eV) is read from nitride.dot's own wetting-layer thickness (0.0 nm on shipped cards -> GaN barrier edge), while the diode's SRH background (fsim_core.nitride_transport.evaluate_injection) lives in drive.diode's separate 0.5 nm InGaN layer -- the two blocks disagree about which material hosts the background carriers. No numerical consequence today (the reservoir energy only sets the cavity/slit spectral ACCEPTANCE of the SRH rate, not the rate itself), but the two should eventually be unified.
