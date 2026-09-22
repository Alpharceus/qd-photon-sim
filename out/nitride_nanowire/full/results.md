# Nitride nanowire sweep results

Independent predictions (design brief User answer 2): every number below is an evaluate() output or a declared, traceable reduction of one (row_id in sweep.csv); no plotted point or headline number is fitted to the Deshpande 2013/2014 held-out lifetime/g2/wavelength anchors. `family` is horizontal_as_built (Deshpande 2013/2014 as-built dispersed wire) or vertical_photonic (designed HE11 wire, headline-eligible only when single_mode AND approximation_error==0); `regime` is rectangular (100 ps electrical pulse) or deterministic_pair (idealized one-pair-per-cycle SET loading); `strain_bound` unrelaxed is the conservative_lower scenario, relaxed the headline_upper scenario -- neither is asserted a rigorous flux bound, and a reversal is reported (bound_reversal_pair) rather than hidden. Two independent, non-gating hardware screens (Coulomb-blockade `set_feasible`/`set_EC_over_kT` and resonant-tunnelling-injector `rti_feasible`) are reported alongside the idealized optical statistics; neither overwrites `g2_op`/`collected_flux_pulsed_s`.

## VERDICT lines (16: family x regime x strain_bound x rep_rate_hz)

Definition: paired_optical_pass counts an optical_pass row whose strain-bound PARTNER row (same family/regime/rep_rate_hz/core_radius_nm/height_nm/x_in/T_hs, opposite strain_bound) is present in this run -- NOT an AND requiring both strain scenarios to individually pass.

VERDICT: idealized_status=no_idealized_pass family=horizontal_as_built regime=rectangular strain_bound=unrelaxed bound_role=conservative_lower rep_rate_hz=8e+07 complete=True eligible=0 paired_optical_pass=0 quality_pass=0 hardware_qualified=0 rti_qualified=0 coverage=192/192 invalid=48 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=no_idealized_pass family=horizontal_as_built regime=rectangular strain_bound=unrelaxed bound_role=conservative_lower rep_rate_hz=2e+08 complete=True eligible=0 paired_optical_pass=0 quality_pass=0 hardware_qualified=0 rti_qualified=0 coverage=192/192 invalid=48 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=no_idealized_pass family=horizontal_as_built regime=rectangular strain_bound=relaxed bound_role=headline_upper rep_rate_hz=8e+07 complete=True eligible=192 paired_optical_pass=0 quality_pass=0 hardware_qualified=0 rti_qualified=0 coverage=192/192 invalid=0 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=no_idealized_pass family=horizontal_as_built regime=rectangular strain_bound=relaxed bound_role=headline_upper rep_rate_hz=2e+08 complete=True eligible=192 paired_optical_pass=0 quality_pass=0 hardware_qualified=0 rti_qualified=0 coverage=192/192 invalid=0 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=pass_hardware_infeasible family=horizontal_as_built regime=deterministic_pair strain_bound=unrelaxed bound_role=conservative_lower rep_rate_hz=8e+07 complete=True eligible=3 paired_optical_pass=3 quality_pass=0 hardware_qualified=0 rti_qualified=0 coverage=192/192 invalid=48 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=pass_hardware_infeasible family=horizontal_as_built regime=deterministic_pair strain_bound=unrelaxed bound_role=conservative_lower rep_rate_hz=2e+08 complete=True eligible=14 paired_optical_pass=14 quality_pass=0 hardware_qualified=0 rti_qualified=0 coverage=192/192 invalid=48 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=pass_hardware_infeasible family=horizontal_as_built regime=deterministic_pair strain_bound=relaxed bound_role=headline_upper rep_rate_hz=8e+07 complete=True eligible=192 paired_optical_pass=58 quality_pass=34 hardware_qualified=0 rti_qualified=0 coverage=192/192 invalid=0 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=pass_hardware_infeasible family=horizontal_as_built regime=deterministic_pair strain_bound=relaxed bound_role=headline_upper rep_rate_hz=2e+08 complete=True eligible=192 paired_optical_pass=19 quality_pass=0 hardware_qualified=0 rti_qualified=0 coverage=192/192 invalid=0 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=no_idealized_pass family=vertical_photonic regime=rectangular strain_bound=unrelaxed bound_role=conservative_lower rep_rate_hz=8e+07 complete=True eligible=0 paired_optical_pass=0 quality_pass=0 hardware_qualified=0 rti_qualified=0 coverage=128/128 invalid=10 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=no_idealized_pass family=vertical_photonic regime=rectangular strain_bound=unrelaxed bound_role=conservative_lower rep_rate_hz=2e+08 complete=True eligible=0 paired_optical_pass=0 quality_pass=0 hardware_qualified=0 rti_qualified=0 coverage=128/128 invalid=3 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=no_idealized_pass family=vertical_photonic regime=rectangular strain_bound=relaxed bound_role=headline_upper rep_rate_hz=8e+07 complete=True eligible=128 paired_optical_pass=0 quality_pass=0 hardware_qualified=0 rti_qualified=0 coverage=128/128 invalid=0 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=no_idealized_pass family=vertical_photonic regime=rectangular strain_bound=relaxed bound_role=headline_upper rep_rate_hz=2e+08 complete=True eligible=128 paired_optical_pass=0 quality_pass=0 hardware_qualified=0 rti_qualified=0 coverage=128/128 invalid=0 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=no_idealized_pass family=vertical_photonic regime=deterministic_pair strain_bound=unrelaxed bound_role=conservative_lower rep_rate_hz=8e+07 complete=True eligible=0 paired_optical_pass=0 quality_pass=0 hardware_qualified=0 rti_qualified=0 coverage=128/128 invalid=0 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=no_idealized_pass family=vertical_photonic regime=deterministic_pair strain_bound=unrelaxed bound_role=conservative_lower rep_rate_hz=2e+08 complete=True eligible=0 paired_optical_pass=0 quality_pass=0 hardware_qualified=0 rti_qualified=0 coverage=128/128 invalid=0 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=pass_hardware_infeasible family=vertical_photonic regime=deterministic_pair strain_bound=relaxed bound_role=headline_upper rep_rate_hz=8e+07 complete=True eligible=128 paired_optical_pass=128 quality_pass=87 hardware_qualified=0 rti_qualified=0 coverage=128/128 invalid=0 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=pass_hardware_infeasible family=vertical_photonic regime=deterministic_pair strain_bound=relaxed bound_role=headline_upper rep_rate_hz=2e+08 complete=True eligible=128 paired_optical_pass=62 quality_pass=21 hardware_qualified=0 rti_qualified=0 coverage=128/128 invalid=0 flux_floor=1000/s screening=0 access=0.05

## Best passing flux per family

BEST_PASSING_FLUX family=horizontal_as_built value=2.27654e+06 row_id=CO00407 core_radius_nm=12.5 height_nm=3 x_in=0.25 T_hs=273 rep_rate_hz=8e+07 strain_bound=relaxed screening=0 regime=deterministic_pair commanded_flux=2.27654e+06 delivered_flux=84941.3 photons_per_cycle_commanded=0.0284568 photons_per_cycle_delivered=0.00106177 (one delivered photon per 941.8 cycles) headline_eligible=True hardware_qualified=False rti_qualified=False blocked_load_probability=7.40345e-10 quality_pass=True unrelaxed_partner_row_id=CO00405 unrelaxed_partner_flux=1.7782 bound_reversal=not_comparable optical_pass_candidates=94 headline_eligible_candidates=94
BEST_PASSING_FLUX_300K family=horizontal_as_built value=1.97467e+06 row_id=CO00223 core_radius_nm=10 height_nm=4 x_in=0.25 T_hs=300 rep_rate_hz=8e+07 strain_bound=relaxed screening=0 regime=deterministic_pair commanded_flux=1.97467e+06 delivered_flux=116045 photons_per_cycle_commanded=0.0246834 photons_per_cycle_delivered=0.00145056 (one delivered photon per 689.4 cycles) headline_eligible=True hardware_qualified=False rti_qualified=False blocked_load_probability=4.29242e-11 quality_pass=True unrelaxed_partner_row_id=CO00221 unrelaxed_partner_flux=0.0450371 bound_reversal=not_comparable optical_pass_candidates=39 headline_eligible_candidates=39
BEST_PASSING_FLUX family=vertical_photonic value=5.97635e+06 row_id=CO01992 core_radius_nm=80 height_nm=4 x_in=0.25 T_hs=230 rep_rate_hz=2e+08 strain_bound=relaxed screening=0 regime=deterministic_pair commanded_flux=5.97635e+06 delivered_flux=5.46023e+06 photons_per_cycle_commanded=0.0298817 photons_per_cycle_delivered=0.0273012 (one delivered photon per 36.63 cycles) headline_eligible=True hardware_qualified=False rti_qualified=False blocked_load_probability=5.77172e-11 quality_pass=True unrelaxed_partner_row_id=CO01990 unrelaxed_partner_flux=0.00530137 bound_reversal=not_comparable optical_pass_candidates=190 headline_eligible_candidates=79
BEST_PASSING_FLUX_300K family=vertical_photonic value=4.35616e+06 row_id=CO02239 core_radius_nm=100 height_nm=3 x_in=0.4 T_hs=300 rep_rate_hz=8e+07 strain_bound=relaxed screening=0 regime=deterministic_pair commanded_flux=4.35616e+06 delivered_flux=3.5296e+06 photons_per_cycle_commanded=0.054452 photons_per_cycle_delivered=0.04412 (one delivered photon per 22.67 cycles) headline_eligible=True hardware_qualified=False rti_qualified=False blocked_load_probability=1.67158e-16 quality_pass=True unrelaxed_partner_row_id=CO02237 unrelaxed_partner_flux=0.0674613 bound_reversal=not_comparable optical_pass_candidates=48 headline_eligible_candidates=22

## Bounds table (conservative unrelaxed vs headline relaxed, reference geometry)

| family | regime | rep_rate_hz | T_hs K | unrelaxed g2 | unrelaxed flux/s | relaxed g2 | relaxed flux/s | bound_reversal | row_ids |
|---|---|---|---|---|---|---|---|---|---|
| horizontal_as_built | rectangular | 8e+07 | 230 | 0.9997970764002508 | 17.434502542020805 | 0.690346256114617 | 3454359.421252767 | not_comparable | CO00353,CO00355 |
| horizontal_as_built | rectangular | 8e+07 | 300 | 0.9999799925881296 | 5.03217836049003 | 0.6904497439787254 | 3463381.2718472197 | not_comparable | CO00377,CO00379 |
| horizontal_as_built | rectangular | 2e+08 | 230 | 0.9997986248418571 | 43.40659404724062 | 0.675598053999424 | 8318636.925812878 | not_comparable | CO00354,CO00356 |
| horizontal_as_built | rectangular | 2e+08 | 300 | 0.9999800961852406 | 12.545291309052995 | 0.6757369210026002 | 8341278.017765789 | not_comparable | CO00378,CO00380 |
| horizontal_as_built | deterministic_pair | 8e+07 | 230 | 0.17355371900826433 | 375.87290058394416 | 0.179324136421946 | 3457024.318770842 | not_comparable | CO00357,CO00359 |
| horizontal_as_built | deterministic_pair | 8e+07 | 300 | 0.17355371900826455 | 117.97525883530035 | 0.1792890886443651 | 3467714.52640254 | not_comparable | CO00381,CO00383 |
| horizontal_as_built | deterministic_pair | 2e+08 | 230 | 0.17355371900826433 | 935.9896504626765 | 0.3168346309252512 | 8630544.613331493 | not_comparable | CO00358,CO00360 |
| horizontal_as_built | deterministic_pair | 2e+08 | 300 | 0.17355371900826455 | 294.1896059761477 | 0.3165284147670757 | 8657131.143435951 | not_comparable | CO00382,CO00384 |
| vertical_photonic | rectangular | 8e+07 | 230 | 0.9999999190333605 | 0.7323628286480333 | 0.6990997824990145 | 4224580.724281467 | not_comparable | CO01889,CO01891 |
| vertical_photonic | rectangular | 8e+07 | 300 | 0.999999992048697 | 0.2152402640389279 | 0.7002722427091184 | 3947403.003283787 | not_comparable | CO01913,CO01915 |
| vertical_photonic | rectangular | 2e+08 | 230 | 0.9999999190356622 | 1.8308783417608352 | 0.698095799405805 | 10532068.02444331 | not_comparable | CO01890,CO01892 |
| vertical_photonic | rectangular | 2e+08 | 300 | 0.9999999920488613 | 0.5380951210996248 | 0.6995223665864898 | 9847688.958371958 | not_comparable | CO01914,CO01916 |
| vertical_photonic | deterministic_pair | 8e+07 | 230 | 0.17355371900826455 | 15.745684531805457 | 0.1735647853023734 | 4159250.1925281966 | not_comparable | CO01893,CO01895 |
| vertical_photonic | deterministic_pair | 8e+07 | 300 | 0.17355371900826455 | 5.070544439430814 | 0.17355930780084994 | 3886469.141773661 | not_comparable | CO01917,CO01919 |
| vertical_photonic | deterministic_pair | 2e+08 | 230 | 0.17355371900826433 | 39.36362623729463 | 0.1867292085013761 | 10398117.084807843 | not_comparable | CO01894,CO01896 |
| vertical_photonic | deterministic_pair | 2e+08 | 300 | 0.17355371900826455 | 12.676243623410885 | 0.18360081169386988 | 9716146.139613492 | not_comparable | CO01918,CO01920 |

Comparable pairs (relaxed count): 17/1280 (17 reversed, 0 not reversed); not_comparable=1263/1280 (no pair where BOTH strain-bound rows clear the 1000/s optical floor -- the legend value for every 'not_comparable' cell/entry elsewhere in this report and in strain_reversal_map.png). Reversal-trigger tally among the 17 reversed pairs: mu=0, collected_flux_pulsed_s=0, g2_op=17 -- every reversed pair triggers on g2_op ONLY, never on mu or collected_flux_pulsed_s.

## Repetition-rate sensitivity: 80 MHz vs 200 MHz SET one_pair_valid (family-asymmetric)

one_pair_valid requires blocked_load_probability<=1e-9 (the disc must empty between cycles); a shorter 200 MHz period leaves less time per cycle for that reset than 80 MHz, so a one_pair_valid failure at 200 MHz where it holds at 80 MHz is a loading-window/period effect, not a fit. collected_flux_pulsed_s is the idealized/commanded flux, collected_flux_delivered_s is the RC-limited delivered flux (bullet 6). No row at any OFF-GRID repetition rate was evaluated this run -- only this run's own main-grid values, 80 MHz and 200 MHz, were ever sampled.

horizontal_as_built T_hs=230K: at 80 MHz one_pair_valid=False optical_pass=False blocked_load_probability=5.752516105036039e-08 (row CO00359, g2=0.179324136421946, flux=3457024.318770842, delivered=122866.09454133948); at 200 MHz one_pair_valid=False optical_pass=False blocked_load_probability=0.0014040041293181974 (row CO00360, g2=0.3168346309252512, flux=8630544.613331493, delivered=306776.87976046297).
horizontal_as_built T_hs=300K: at 80 MHz one_pair_valid=False optical_pass=False blocked_load_probability=5.647793631516201e-08 (row CO00383, g2=0.1792890886443651, flux=3467714.52640254, delivered=128792.56555700889); at 200 MHz one_pair_valid=False optical_pass=False blocked_load_probability=0.0013932024440156097 (row CO00384, g2=0.3165284147670757, flux=8657131.143435951, delivered=321566.51162706804).
vertical_photonic T_hs=230K: at 80 MHz one_pair_valid=True optical_pass=True blocked_load_probability=3.994511988600227e-16 (row CO01895, g2=0.1735647853023734, flux=4159250.1925281966, delivered=3767861.7438239465); at 200 MHz one_pair_valid=False optical_pass=False blocked_load_probability=6.988915925889203e-07 (row CO01896, g2=0.1867292085013761, flux=10398117.084807843, delivered=9419648.437798565).
vertical_photonic T_hs=300K: at 80 MHz one_pair_valid=True optical_pass=True blocked_load_probability=5.145228695925724e-17 (row CO01919, g2=0.17355930780084994, flux=3886469.141773661, delivered=3574053.9491063515); at 200 MHz one_pair_valid=False optical_pass=False blocked_load_probability=3.072275330426634e-07 (row CO01920, g2=0.18360081169386988, flux=9716146.139613492, delivered=8935111.502103966).

## Deshpande 2014 comparison (x_in=0.40, T_hs=300 K, 200 MHz, deterministic_pair)

Measured (abstract-only [V], CONDITIONS INCOMPLETE): lambda~630 nm, lifetime 1.3+/-0.3 ns, g2=0.29.

| strain_bound | lambda_nm (predicted) | tau_rad_bare_ns | tau_rad_photonic_ns | tau_total_X_ns | g2_op | collected_flux_pulsed_s | row_id |
|---|---|---|---|---|---|---|---|
| unrelaxed | 622.8766000859592 | 24.450986850391498 | 62.259743487195735 | 0.0017032738525559954 | 0.17355371900826455 | 294.1896059761477 | CO00382 |
| relaxed | 542.6613131119037 | 1.0712970826414099 | 2.7454237780315776 | 2.2463730632050725 | 0.3165284147670757 | 8657131.143435951 | CO00384 |

Untouched OLD planar 2014 comparison replay (non-gating, unchanged from the prior piece): row LI02992, g2=0.9980609340138185, flux=0.018999124213812685 vs measured g2=0.29.

## Deshpande 2013 comparison (pulsed replay of a CW measurement, 10 K, drive_mismatch)

Measured (CW electrical, 10 K, [V]): X g2 raw/corrected 0.30/0.16, XX 0.38/0.25 at 1 nA; g2-fit lifetimes X 1.1 ns, XX 0.7 ns; TRPL XX 711 ps; emission X=2.84 eV (436.56 nm). This evaluator has no CW/HBT drive path (pulsed rectangular / deterministic_pair only); the rows below are the pulsed model's own prediction at the paper's geometry/current, published as drive_mismatch -- never substituted for a predicted CW g2.

| R nm | I uA | lambda_nm (predicted) | tau_rad_bare_ns | tau_rad_photonic_ns | tau_total_X_ns | g2_op (drive_mismatch) | commanded_flux/s | time_avg_current_pA | valid | row_id |
|---|---|---|---|---|---|---|---|---|---|---|
| 12.5 | 0.001 | 446.41772793972837 | 1.066861783858344 | 2.7732036580958352 | 2.269664341693179 | 0.7917008239538932 (drive_mismatch) | 1390068.9541082678 | 8 | True | D1302972 |
| 12.5 | 0.002 | 446.4177252000232 | 1.0668421845337484 | 2.773152713169029 | 2.269630217520449 | 0.6903016283384055 (drive_mismatch) | 2536737.5901067983 | 16 | True | D1302973 |
| 15.0 | 0.001 | 446.2330005160856 | 1.0666577502540837 | 2.772782270671441 | 2.34019262863002 | 0.7910213051220907 (drive_mismatch) | 1369801.4318615806 | 8 | True | D1302974 |
| 15.0 | 0.002 | 446.2329979828811 | 1.0666381870417456 | 2.772731417496704 | 2.3401564050816375 | 0.6895194267777083 (drive_mismatch) | 2498202.7016585306 | 16 | True | D1302975 |

tau_rad_bare_ns / tau_rad0_ns (this run's own [A] input) ratio: 1.067-1.067 across these rows -- the bare radiative lifetime's closeness to the 1.1 ns anchor reflects this fixed [A] input times a near-unity geometry factor, NOT an independently verified radiative-rate prediction; tau_rad_photonic_ns (after antenna suppression) and tau_total_X_ns (after nonradiative/surface competition) are the physically relevant, larger, device-level lifetimes and are reported alongside it, never substituted for it.
Drive AND heating mismatch: these rows run 100 ps rectangular pulses at 80 MHz (this run's own I_uA/tau_pulse_ns/rep_rate_hz columns above give the time-average current column), roughly 2 orders of magnitude below the paper's CW 1 nA -- both the counting statistics (pulsed vs CW g2) and the junction/thermal operating point (heating scales with time-average power) differ from the measured device; g2_op above is published as drive_mismatch, never a predicted CW g2.

## Opposite-endpoint anchor match (bullet 11 obligation)

2013 relaxed predicted lambda_nm=446.41772793972837 vs measured ~437 nm (X=2.84 eV, [V] Fig. 3c). 2014 relaxed predicted lambda_nm=542.6613131119037 / unrelaxed predicted lambda_nm=622.8766000859592 vs measured ~630 nm. The two strain-bound anchors are matched by OPPOSITE endpoints (2013 by relaxed, 2014 by whichever endpoint lands closer) and are never averaged: one of x_in transfer, disc thickness, VBO/bowing, or lateral localization carries roughly 0.3 eV of the remaining discrepancy -- stated here, not resolved.

## RC caveat (bullet 6 obligation)

As-built R_s_ohm=2.38e9 (horizontal_as_built, deterministic_pair rows): tau_RC_ns spans 1.651-26.47, delivered_step_fraction spans 0.00377-0.05877 -- this device CANNOT deliver a 100 ps step.
Designed-contact R_s_ohm=1e6 (sensitivity rows): tau_RC_ns spans 0.00111-0.001161, delivered_step_fraction spans 1-1 -- only this designed contact can deliver the 100 ps step.
These numbers are THIS RUN's own tau_RC_ns/delivered_step_fraction, recomputed from sweep.csv rows, not restated from the contract: docs/nitride_nanowire_contract.md bullet 6 quotes tau_RC_ns 0.82-1.19 / delivered_step_fraction 0.08-0.12 for the as-built device, which this run's own as-built range above SUPERSEDES (the contract text is out of scope for this piece and is not edited here).

## Access-1.0 lifetime cap (bullet 8 obligation)

occupied_dot_access=1.0 at core_radius_nm=12.5 (row SN02753): 1/k_surface_X_ns=0.625 ns, 1/k_surface_XX_ns=0.3125 ns (contract's stated cap: 0.625 ns X, 0.3125 ns XX) -- this caps tau_X/tau_XX regardless of any other physics; the headline default occupied_dot_access=0.05 is reported next to it, never chosen to reproduce the held-out anchors.

## E_C/kT wall (bullet 11 obligation)

Across every deterministic_pair core row with core_radius_nm>=10 nm at 230-300 K: set_EC_over_kT spans 0.1464-0.7643 (required ec_margin=10), set_feasible=True count=0/1184; rti_feasible=True count=0/1136. Deterministic loading at 230-300 K fails the Coulomb-blockade screen for every core_radius_nm>=10 nm priced in this tier.
The RT injector screen is SEPARATELY reported, not folded into the Coulomb wall above: rti_status=unknown_incomplete on 1184/1280 deterministic_pair rows, rti_transport_feasible=False on 1280/1280 -- these are conditional engineering screens on an unsupported occupation/second-pair control, NOT a demonstrated hardware failure by either charging mechanism; deterministic loading at 230-300 K fails on the Coulomb-blockade wall (set_feasible) and separately carries an incomplete/unresolved RT-injector screen, and neither screen overwrites optical_pass/g2_op/collected_flux_pulsed_s computed upstream of it.
Both priced charging-based loading mechanisms share the SAME insufficient disc E_C/kT: the RT-injector screen's own second_pair_addition_meV input is set_E_C_meV -- the identical Coulomb charging energy the set_EC_over_kT wall above already reports as an E_C/kT<10 failure at every core_radius_nm>=10 nm priced here. This is a CONDITIONAL MODEL LIMITATION shared by both screens' inputs (contract bullet 11's E_C/kT wall statement covers 'any charging mechanism priced in this tier'), never a demonstrated hardware failure of the RT-injector mechanism specifically (contract bullet 7: rti_feasible=False must not be reported as a demonstrated physics result).

## c-plane dipole prior falsification (Composition rules bullet 9)

isotropic default (row DW02970): degree_of_linear_polarization=0.8399435321115863 (anchor: +70% axial, deshpande2013_polarization [V]). CPLANE_ONLY_DIPOLE_WEIGHTS=(0,0.5,0.5) (row DW02971): degree_of_linear_polarization=-0.9588380316429788 -- the wrong sign against the +70% anchor, confirming this sensitivity stays a named, explicitly falsified alternative and never the headline default.

## Gate anti-monotonicity (H3 obligation)

one_pair_valid can be satisfied by ADDING sidewall loss (faster occupied-dot emptying) as readily as by improving device quality -- optical_pass/paired_optical_pass alone therefore do NOT certify throughput. quality_pass=optical_pass AND photons_per_cycle>=0.01 [A, orchestrator threshold: one collected photon per hundred cycles] is reported beside optical_pass/paired_optical_pass in every VERDICT line and the per-temperature tables below. M4 fix: photons_per_cycle (sweep.csv column) and quality_pass BOTH use COMMANDED flux (collected_flux_pulsed_s/rep_rate_hz), never the RC-limited delivered flux -- the BEST lines below additionally print a report-derived photons_per_cycle_delivered (collected_flux_delivered_s/rep_rate_hz) beside it so the commanded per-cycle number is never mistaken for what an RC-limited detector would actually see.

S_cm_s pair that actually DEMONSTRATES the anti-monotonicity (selected programmatically: the (T_hs, rep_rate_hz) combination at the horizontal reference geometry where one_pair_valid flips False->True while collected_flux_pulsed_s FALLS; the 300K/200MHz reference-condition triple below does NOT demonstrate it -- one_pair_valid is False at all three S_cm_s values there):
- S_cm_s=100.0 (lower surface loss, T_hs=230K, 80MHz) (row SN02681): one_pair_valid=False g2_op=0.18766432896768592 flux=4124930.72172916 quality_pass=False
- S_cm_s=10000.0 (higher surface loss, T_hs=230K, 80MHz) (row SN02682): one_pair_valid=True g2_op=0.17355443393845815 flux=1319877.5375752521 quality_pass=True
quality_pass (photons_per_cycle>=0.01) removes only ABSOLUTELY DIM optical passes; it does NOT repair this loss-induced ordering: row SN02682 (S_cm_s=10000.0) keeps quality_pass=True even though the 100x surface-recombination increase (S_cm_s 100.0->10000.0) cut its commanded flux 3.1x (4.125e+06->1.32e+06 /s).

S_cm_s and occupied_dot_access triple/pair at the 300K/200MHz reference condition, for context only (S_cm_s here does NOT flip; occupied_dot_access DOES):
- S_cm_s=100 (lower surface loss) (row SN02687): one_pair_valid=False g2_op=0.3656546763984454 flux=10297783.614345277 quality_pass=False
- S_cm_s=1000 (main-grid default) (row CO00384): one_pair_valid=False g2_op=0.3165284147670757 flux=8657131.143435951 quality_pass=False
- S_cm_s=10000 (higher surface loss) (row SN02688): one_pair_valid=False g2_op=0.17797379959021353 flux=3312159.633154987 quality_pass=False
- occupied_dot_access=0.05 (main-grid default, lower loss) (row CO00384): one_pair_valid=False g2_op=0.3165284147670757 flux=8657131.143435951 quality_pass=False
- occupied_dot_access=1.0 (higher loss) (row SN02768): one_pair_valid=True g2_op=0.1736349147052857 flux=1963808.9310572965 quality_pass=False

Per-T_hs optical_pass counts, horizontal_as_built (both regimes, split by strain_bound), with x_in composition (read from the rows, never hardcoded):
- horizontal_as_built/rectangular/unrelaxed: T=230K:0(x_in0.25=0,x_in0.40=0) T=250K:0(x_in0.25=0,x_in0.40=0) T=273K:0(x_in0.25=0,x_in0.40=0) T=300K:0(x_in0.25=0,x_in0.40=0)
- horizontal_as_built/rectangular/relaxed: T=230K:0(x_in0.25=0,x_in0.40=0) T=250K:0(x_in0.25=0,x_in0.40=0) T=273K:0(x_in0.25=0,x_in0.40=0) T=300K:0(x_in0.25=0,x_in0.40=0)
- horizontal_as_built/deterministic_pair/unrelaxed: T=230K:7(x_in0.25=0,x_in0.40=7) T=250K:5(x_in0.25=0,x_in0.40=5) T=273K:3(x_in0.25=0,x_in0.40=3) T=300K:2(x_in0.25=0,x_in0.40=2)
- horizontal_as_built/deterministic_pair/relaxed: T=230K:3(x_in0.25=3,x_in0.40=0) T=250K:11(x_in0.25=11,x_in0.40=0) T=273K:26(x_in0.25=26,x_in0.40=0) T=300K:37(x_in0.25=37,x_in0.40=0)

Composition of optical passes, per family (counted from the rows; see the definition of 'entirely'/'both compositions' inline):
HORIZONTAL (horizontal_as_built) RELAXED optical passes: 77 total (x_in=0.25: 77, x_in=0.40: 0) -- ENTIRELY x_in=0.25, no relaxed x_in=0.40 optical pass this run.
HORIZONTAL (horizontal_as_built) x_in=0.40 UNRELAXED optical passes: 17 total (the only x_in=0.40 passes in THIS family), collected_flux_pulsed_s spans 1014-5290 /s, row_ids=CO00037,CO00038,CO00045,CO00046,CO00054,CO00062,CO00102,CO00110....
VERTICAL (vertical_photonic) RELAXED optical passes: 190 total (x_in=0.25: 126, x_in=0.40: 64) -- this family passes at BOTH compositions this run.
VERTICAL (vertical_photonic) x_in=0.40 UNRELAXED optical passes: 0 total.
Across BOTH families this run: 81/284 optical passes are x_in=0.40 -- the horizontal family's x_in=0.40 passes are UNRELAXED-only, while the vertical family also passes at x_in=0.40 under the RELAXED bound (see the per-family lines above); the earlier 'these are the only x_in=0.40 passes in either family' claim was HORIZONTAL-only and is not repeated here.

## Per-temperature counts (M9-M10 obligation)

| family | regime | T_hs K | n_core | n_eligible | n_optical_pass | n_quality_pass |
|---|---|---|---|---|---|---|
| horizontal_as_built | rectangular | 230 | 192 | 96 | 0 | 0 |
| horizontal_as_built | rectangular | 250 | 192 | 96 | 0 | 0 |
| horizontal_as_built | rectangular | 273 | 192 | 96 | 0 | 0 |
| horizontal_as_built | rectangular | 300 | 192 | 96 | 0 | 0 |
| horizontal_as_built | deterministic_pair | 230 | 192 | 103 | 10 | 2 |
| horizontal_as_built | deterministic_pair | 250 | 192 | 101 | 16 | 8 |
| horizontal_as_built | deterministic_pair | 273 | 192 | 99 | 29 | 14 |
| horizontal_as_built | deterministic_pair | 300 | 192 | 98 | 39 | 10 |
| vertical_photonic | rectangular | 230 | 128 | 64 | 0 | 0 |
| vertical_photonic | rectangular | 250 | 128 | 64 | 0 | 0 |
| vertical_photonic | rectangular | 273 | 128 | 64 | 0 | 0 |
| vertical_photonic | rectangular | 300 | 128 | 64 | 0 | 0 |
| vertical_photonic | deterministic_pair | 230 | 128 | 64 | 46 | 44 |
| vertical_photonic | deterministic_pair | 250 | 128 | 64 | 48 | 32 |
| vertical_photonic | deterministic_pair | 273 | 128 | 64 | 48 | 16 |
| vertical_photonic | deterministic_pair | 300 | 128 | 64 | 48 | 16 |

## 16 per-family/regime/bound/rate reference-headline nominations (M9-M10 obligation)

One nominated row per (family,regime,strain_bound,rep_rate_hz) group -- the brightest optical_pass row in that group (headline_eligible additionally required for vertical_photonic, H2/bullet 10); 'none' if the group has no such row.

| family | regime | strain_bound | rep_rate_hz | row_id | commanded_flux/s | photons_per_cycle_delivered | blocked_load_probability | g2_op | headline_eligible | quality_pass |
|---|---|---|---|---|---|---|---|---|---|---|
| horizontal_as_built | rectangular | unrelaxed | 8e+07 | none | | | | | False | |
| horizontal_as_built | rectangular | unrelaxed | 2e+08 | none | | | | | False | |
| horizontal_as_built | rectangular | relaxed | 8e+07 | none | | | | | False | |
| horizontal_as_built | rectangular | relaxed | 2e+08 | none | | | | | False | |
| horizontal_as_built | deterministic_pair | unrelaxed | 8e+07 | CO00037 | 2124.32 | 1.421e-06 | 0 | 0.17355371900826455 | True | False |
| horizontal_as_built | deterministic_pair | unrelaxed | 2e+08 | CO00038 | 5290.03 | 1.415e-06 | 0 | 0.17355371900826455 | True | False |
| horizontal_as_built | deterministic_pair | relaxed | 8e+07 | CO00407 | 2.27654e+06 | 0.001062 | 7.403e-10 | 0.17491140858775645 | True | True |
| horizontal_as_built | deterministic_pair | relaxed | 2e+08 | CO00992 | 1.91507e+06 | 0.0001517 | 6.815e-10 | 0.17487446045397337 | True | False |
| vertical_photonic | rectangular | unrelaxed | 8e+07 | none | | | | | False | |
| vertical_photonic | rectangular | unrelaxed | 2e+08 | none | | | | | False | |
| vertical_photonic | rectangular | relaxed | 8e+07 | none | | | | | False | |
| vertical_photonic | rectangular | relaxed | 2e+08 | none | | | | | False | |
| vertical_photonic | deterministic_pair | unrelaxed | 8e+07 | none | | | | | False | |
| vertical_photonic | deterministic_pair | unrelaxed | 2e+08 | none | | | | | False | |
| vertical_photonic | deterministic_pair | relaxed | 8e+07 | CO02279 | 4.59824e+06 | 0.04571 | 1.747e-13 | 0.17363770828071723 | True | True |
| vertical_photonic | deterministic_pair | relaxed | 2e+08 | CO01992 | 5.97635e+06 | 0.0273 | 5.772e-11 | 0.17413411911590182 | True | True |

## Ensemble-yield proxy at 10K/300K vs measured 0.52 (M9-M10 obligation, non-gating)

Inputs (this run's own gamma_X_ns radiative / k_X_ns nonradiative additive rate pair, rows CO00347 300K / D1302973 10K): gamma_300_ns=0.35922, loss_300_ns=0.630736, gamma_10_ns=0.3606, loss_10_ns=0.08.
yield_300K=0.362864, yield_10K=0.81843, MODEL PROXY ratio=0.443367 -- vs measured ensemble PL ratio 0.52 (Deshpande 2013, ledger deshpande2013_thermal_and_pl). This proxy assumes equal absorption/capture/collection between 10K and 300K [A] and is a SINGLE-DOT proxy, never verification of the ensemble (many-wire) PL measurement; no channel in this model was fit to reproduce 0.52, and none does.

## Planar-reference rows and cross-round reconciliation (M9-M10 obligation)

| family_tag | regime | T_hs K | screening | g2_op | collected_flux_pulsed_s | valid | row_id |
|---|---|---|---|---|---|---|---|
| c_plane | deterministic_pair | 230.0 | 0.0 | 0.17355371900826433 | 0.001255189394466358 | True | PR02980 |
| c_plane | deterministic_pair | 230.0 | 1.0 | 0.17355371900826433 | 1731.5228538626084 | True | PR02981 |
| c_plane | deterministic_pair | 300.0 | 0.0 | 0.17355371900826455 | 0.13742222076335409 | True | PR02982 |
| c_plane | deterministic_pair | 300.0 | 1.0 | 0.17355371900826433 | 25642.788985855634 | True | PR02983 |
| c_plane | rectangular | 230.0 | 0.0 | 0.9991379030238006 | 3.880173957537539e-05 | True | PR02976 |
| c_plane | rectangular | 230.0 | 1.0 | 0.5763774545747093 | 2.8431507378029667e-08 | True | PR02977 |
| c_plane | rectangular | 300.0 | 0.0 | 0.9996791936846767 | 0.004247138693754914 | True | PR02978 |
| c_plane | rectangular | 300.0 | 1.0 | 0.8958834914876267 | 1.1821590535478755e-05 | True | PR02979 |
| nonpolar | deterministic_pair | 230.0 | 0.0 | 0.17355371900826455 | 144.8264154955226 | True | PR02988 |
| nonpolar | deterministic_pair | 230.0 | 1.0 | 0.17355371900826455 | 144.8264154955226 | True | PR02989 |
| nonpolar | deterministic_pair | 300.0 | 0.0 | 0.17355371900826455 | 2151.3242842275577 | True | PR02990 |
| nonpolar | deterministic_pair | 300.0 | 1.0 | 0.17355371900826455 | 2151.3242842275577 | True | PR02991 |
| nonpolar | rectangular | 230.0 | 0.0 | 0.8762316964092836 | 5.298015525163914e-10 | True | PR02984 |
| nonpolar | rectangular | 230.0 | 1.0 | 0.8762316964092836 | 5.298015525163914e-10 | True | PR02985 |
| nonpolar | rectangular | 300.0 | 0.0 | 0.9895138559015328 | 3.142264791100479e-07 | True | PR02986 |
| nonpolar | rectangular | 300.0 | 1.0 | 0.9895138559015328 | 3.142264791100479e-07 | True | PR02987 |

Reconciliation: planar c-plane SET, screening=1, 300 K (row PR02983): 25642.788985855634 /s (round 2's own optimized headline was a DIFFERENT, separately optimized design point ~28 kHz/s, not this screened-default card row -- the two must not be conflated) vs this run's nanowire optical_pass collected_flux_pulsed_s range 1014-6.154e+06 /s. The nanowire tier's higher flux is attributed to: a single-wire supply without the planar aperture partition, relaxed (piezoelectric-field-free) strain raising overlap and suppressing surface escape, and antenna/waveguide collection geometry differing from the planar cavity -- qualitative mechanism attributions, not a controlled one-parameter comparison between the two platforms.

## Numerical sensitivity table, ranked by MEASURED leverage (M9-M10/L obligation)

Ranking is measured from this run's own rows, separately per family and per observable -- never a single mixed narrative ordering.

### horizontal_as_built

Reference row: CO00384 g2_op=0.3165284147670757 commanded_flux=8657131.143435951 delivered_flux=321566.51162706804

Ranked by commanded flux |ratio-1| (measured leverage): occupied_dot_access(0.773), S_cm_s(0.617), shell(0.175), screening_fraction(0.00171), R_s_ohm(8.83e-06), al_fraction(0), alignment_uncertainty_meV(0), b_res(0), dipole_weights(0), growth_tolerance_steps(0), injector_barrier_thickness_nm(0), occupation_control_uncertainty(0)
Ranked by delivered flux |ratio-1| (measured leverage): R_s_ohm(25.9), occupied_dot_access(0.773), S_cm_s(0.617), shell(0.175), screening_fraction(0.00171), al_fraction(0), alignment_uncertainty_meV(0), b_res(0), dipole_weights(0), growth_tolerance_steps(0), injector_barrier_thickness_nm(0), occupation_control_uncertainty(0)
Ranked by g2_op |delta| (measured leverage): occupied_dot_access(0.143), S_cm_s(0.139), b_res(0.098), shell(0.0491), screening_fraction(0.00118), R_s_ohm(5.53e-06), al_fraction(0), alignment_uncertainty_meV(0), dipole_weights(0), growth_tolerance_steps(0), injector_barrier_thickness_nm(0), occupation_control_uncertainty(0)

(the contact R_s_ohm row has ratio~1.0 on commanded flux -- it only changes delivered flux via the RC time constant, never the idealized/commanded flux.)

| axis | value | g2_op | commanded_flux/s | flux_ratio_to_reference | delivered_flux/s | delivered_flux_ratio_to_reference | one_pair_valid | optical_pass | quality_pass | row_id |
|---|---|---|---|---|---|---|---|---|---|---|
| (reference) | main-grid default | 0.3165284147670757 | 8657131.143435951 | 1.0 | 321566.51162706804 | 1.0 | False | False | False | CO00384 |
| occupied_dot_access | 1.0 | 0.1736349147052857 | 1963808.9310572965 | 0.2268 | 72945.08735044302 | 0.2268 | True | True | False | SN02768 |
| S_cm_s | 100.0 | 0.3656546763984454 | 10297783.614345277 | 1.19 | 382508.05024089216 | 1.19 | False | False | False | SN02687 |
| S_cm_s | 10000.0 | 0.17797379959021353 | 3312159.633154987 | 0.3826 | 123029.16538271548 | 0.3826 | False | False | False | SN02688 |
| shell | AlGaN | 0.3656546763984454 | 10171773.983798828 | 1.175 | 377827.46071823477 | 1.175 | False | False | False | SN02736 |
| screening_fraction | 0.0 | 0.3165284147670757 | 8657131.143435951 | 1 | 321566.51162706804 | 1 | False | False | False | SN02591 |
| screening_fraction | 1.0 | 0.3153438449500987 | 8671952.281091604 | 1.002 | 322117.0382682204 | 1.002 | False | False | False | SN02592 |
| R_s_ohm | 1000000.0 | 0.31653394196232165 | 8657207.554713154 | 1 | 8657207.554713154 | 26.92 | False | False | False | SN02960 |
| al_fraction | 0.2 | 0.3165284147670757 | 8657131.143435951 | 1 | 321566.51162706804 | 1 | False | False | False | SN02800 |
| alignment_uncertainty_meV | 30.0 | 0.3165284147670757 | 8657131.143435951 | 1 | 321566.51162706804 | 1 | False | False | False | SN02864 |
| b_res | 0.02 | 0.2184923908909121 | 8657131.143435951 | 1 | 321566.51162706804 | 1 | False | False | False | SN02640 |
| growth_tolerance_steps | 2.0 | 0.3165284147670757 | 8657131.143435951 | 1 | 321566.51162706804 | 1 | False | False | False | SN02832 |
| injector_barrier_thickness_nm | 1.0 | 0.3165284147670757 | 8657131.143435951 | 1 | 321566.51162706804 | 1 | False | False | False | SN02928 |
| occupation_control_uncertainty | True | 0.3165284147670757 | 8657131.143435951 | 1 | 321566.51162706804 | 1 | False | False | False | SN02896 |

### vertical_photonic

Reference row: CO01920 g2_op=0.18360081169386988 commanded_flux=9716146.139613492 delivered_flux=8935111.502103966

Ranked by commanded flux |ratio-1| (measured leverage): occupied_dot_access(0.192), S_cm_s(0.101), shell(0.0425), screening_fraction(0.0205), al_fraction(0), alignment_uncertainty_meV(0), b_res(0), growth_tolerance_steps(0), injector_barrier_thickness_nm(0), occupation_control_uncertainty(0)
Ranked by delivered flux |ratio-1| (measured leverage): occupied_dot_access(0.192), S_cm_s(0.101), shell(0.0425), screening_fraction(0.0205), al_fraction(0), alignment_uncertainty_meV(0), b_res(0), growth_tolerance_steps(0), injector_barrier_thickness_nm(0), occupation_control_uncertainty(0)
Ranked by g2_op |delta| (measured leverage): b_res(0.132), occupied_dot_access(0.00697), S_cm_s(0.00431), screening_fraction(0.000585), shell(0.000577), al_fraction(0), alignment_uncertainty_meV(0), growth_tolerance_steps(0), injector_barrier_thickness_nm(0), occupation_control_uncertainty(0)



| axis | value | g2_op | commanded_flux/s | flux_ratio_to_reference | delivered_flux/s | delivered_flux_ratio_to_reference | one_pair_valid | optical_pass | quality_pass | row_id |
|---|---|---|---|---|---|---|---|---|---|---|
| (reference) | main-grid default | 0.18360081169386988 | 9716146.139613492 | 1.0 | 8935111.502103966 | 1.0 | False | False | False | CO01920 |
| occupied_dot_access | 1.0 | 0.17663266881313455 | 7851673.319151593 | 0.8081 | 7220514.757253847 | 0.8081 | False | False | False | SN02784 |
| S_cm_s | 100.0 | 0.1841778750087586 | 9826678.12954954 | 1.011 | 9036758.352659397 | 1.011 | False | False | False | SN02719 |
| S_cm_s | 10000.0 | 0.17929547780047272 | 8733757.21437212 | 0.8989 | 8031692.136098741 | 0.8989 | False | False | False | SN02720 |
| shell | AlGaN | 0.18417787500875837 | 10128721.810067073 | 1.042 | 9314522.182592627 | 1.042 | False | False | False | SN02752 |
| screening_fraction | 0.0 | 0.18360081169386988 | 9716146.139613492 | 1 | 8935111.502103966 | 1 | False | False | False | SN02623 |
| screening_fraction | 1.0 | 0.18418554198214543 | 9915424.030337261 | 1.021 | 9118370.39384097 | 1.021 | False | False | False | SN02624 |
| al_fraction | 0.2 | 0.18360081169386988 | 9716146.139613492 | 1 | 8935111.502103966 | 1 | False | False | False | SN02816 |
| alignment_uncertainty_meV | 30.0 | 0.18360081169386988 | 9716146.139613492 | 1 | 8935111.502103966 | 1 | False | False | False | SN02880 |
| b_res | 0.02 | 0.051433710524236065 | 9716146.139613492 | 1 | 8935111.502103966 | 1 | False | False | False | SN02656 |
| growth_tolerance_steps | 2.0 | 0.18360081169386988 | 9716146.139613492 | 1 | 8935111.502103966 | 1 | False | False | False | SN02848 |
| injector_barrier_thickness_nm | 1.0 | 0.18360081169386988 | 9716146.139613492 | 1 | 8935111.502103966 | 1 | False | False | False | SN02944 |
| occupation_control_uncertainty | True | 0.18360081169386988 | 9716146.139613492 | 1 | 8935111.502103966 | 1 | False | False | False | SN02912 |

tau_rad0_ns/dipole-prior leverage proxy (from the c-plane dipole-prior falsification rows above): switching dipole_weights from isotropic to CPLANE_ONLY materially changes antenna_rate_factor/tau_rad_photonic_ns/commanded flux (see that section's row values) -- a qualitative leverage indicator, not a numeric ranking-table entry, since it sweeps an orientation prior rather than the tau_rad0_ns magnitude (not sampled this run).

Restored current/pulse-width cut (I x tau_pulse, 9 rows, horizontal reference geometry, rectangular regime, relaxed, 300K, 200MHz -- preserves the 100 ps headline point):

| I_uA | tau_pulse_ns | g2_op | commanded_flux/s | valid | row_id |
|---|---|---|---|---|---|
| 0.001 | 0.01 | 0.9662555392101329 | 493872.55163792474 | True | SN02961 |
| 0.001 | 0.1 | 0.7760225567041406 | 4678461.242980294 | True | SN02962 |
| 0.001 | 1.0 | 0.7221832835909678 | 18863436.764750484 | True | SN02963 |
| 0.002 | 0.01 | 0.935680972607993 | 985834.8098184542 | True | SN02964 |
| 0.002 | 0.1 | 0.6757369210026002 | 8341278.017765789 | True | SN02965 |
| 0.002 | 1.0 | 0.7574275474958108 | 20513980.71067336 | True | SN02966 |
| 0.02 | 0.01 | 0.665518902987269 | 8313735.120467398 | True | SN02967 |
| 0.02 | 0.1 | 0.5749764791782161 | 16001423.511084666 | True | SN02968 |
| 0.02 | 1.0 | 0.9226680814012316 | 12444539.18653733 | True | SN02969 |

## Column/gate definitions (L17 obligation)

`eligible` = valid AND finite g2_op AND finite collected_flux_pulsed_s AND collected_flux_pulsed_s>=1000/s (flux_floor) -- NOT the same as `headline_eligible` (the vertical_photonic single_mode AND approximation_error==0 check, bullet 10) or `optical_pass` (which additionally requires g2_op<0.5 and, for deterministic_pair rows, one_pair_valid AND pair_supply_possible). `device_pass` (a device-module output column, sweep.csv) mirrors the planar convention's Coulomb-screen alias; it is reported as a row column alongside `set_feasible`/`hardware_qualified` but this results text does not gate any VERDICT/BEST/quality_pass computation on it beyond the `set_feasible`/`hardware_qualified` values already discussed above.

## Fitted electrical/thermal transport inputs (L14 obligation)

Distinct from the unfitted optical prediction (levels/photonics/surface: no parameter there is chosen to reproduce a held-out optical anchor): fsim_core/nitride_nanowire_transport.py's NitrideWireDiode.tau_SRH_ns default (0.01 ns, 10 ps) is [A, back-solved] so the GaN-kernel dark current reproduces the empirically inferred 2.6-3.1 V junction-voltage window at 1 nA/300 K (Deshpande 2013 terminal range minus the IR drop) -- used in EVERY row's V_j solve this run, not a measured SRH lifetime. WIRE_RTH_PRESETS_K_W['deshpande_fig4_replay']=2.8e9 K/W is [DR, re-fit] to both Deshpande 2013 Fig.4 heating rises under the H2 GaN-kernel V_j, at the paper's 10 K bath -- a SEPARATE replay anchor from this sweep's own Rth_K_W card defaults (1.0e9 horizontal, 1.0e7 vertical, see the Card schema), which are NOT re-fit to Fig.4 and are the values this sweep's core/sensitivity rows actually use.

## Invalid rows (all kinds)

invalid_total=205 across kinds: core=205

Breakdown by reason (leading token before ':', counted per row -- a row may carry more than one reason):

- si_complex_index: 192
- photon counting did not converge: 13

NOTE: 'si_complex_index' alone accounts for 192/205 (94 percent) of all invalid rows this run -- a single dependency-module data-table coverage gap, generated from the rows themselves: fsim_core/nitride_nanowire_photonics.py's Si-substrate complex-index table covers 380-750 nm (read from the module's own _SI_INDEX_ANCHORS_NM table). The affected core rows (recovered from their own invalid_reasons diagnostic, never a repaired/assumed value) span lambda_nm=898.156-1744.47 nm (red-shifted ABOVE the table) and are composed of family=['horizontal_as_built'], strain_bound=['unrelaxed'], x_in=[0.4], height_nm=[3.0, 4.0] -- not a physics failure this sweep introduces or can fix within its own three scoped files. Reported as the interface coordination item for the orchestrator/reviewer, not silently repaired (no n_substrate override was invented to force these rows valid).

### QCSE excursion physics paragraph (H4 obligation)

The 192 rows above are the unrelaxed (conservative_lower), x_in=0.40, height_nm in [3.0, 4.0] disc geometries. field_kVcm on an INVALID row like these is the transport DEPLETION field alone (fsim_core/nitride_nanowire_device.py's invalid-row fallback captures inj['depletion_field_kVcm'] before the row goes invalid), NOT the built-in piezoelectric field -- this run's own field_kVcm on the representative row (deterministic rule: smallest row_id among the 192 affected rows), CO00161: 43.43 kV/cm depletion field. The actual polarization field driving QCSE is reported separately: a levels-only replay (fsim_core.nitride_nanowire_levels, unrelaxed bound, this row's own x_in/height_nm/core_radius_nm, external_field_kVcm set to this row's OWN recorded field_kVcm as the resolved depletion field, per device.py's own post-feedback convention) gives F_pz_kVcm=-6582.28, total_field_kVcm=-6426.02. This field drives the quantum-confined Stark effect (QCSE) far enough to red-shift emission to 0.898-1.74 um, past the photonics module's own 750 nm table ceiling -- an excursion this run reports as MISSING (invalid), never silently repaired. Their relaxed strain-bound partner (strain_fraction=0, no piezoelectric field) is valid for 192/192 of these rows and invalid/absent for 0/192 -- i.e. the conservative (unrelaxed) bound is ABSENT at exactly these x_in=0.40, thick-disc geometries where only the relaxed bound is computable this run, an asymmetry the bounds table and headline selection must not paper over by silently reporting only the relaxed side.

## Decisions (conservative choices under ambiguity)

- Measured evaluate() cost on this machine: 6.1-9.8 s/call serially (fsim_core.nitride_nanowire_injector.injector_feasibility's resonance search dominates, ~85 percent of wall time, run on EVERY row regardless of regime -- a dependency-module property, not something this piece may patch). At that cost the mandated 2560-row core grid cannot finish serially inside 1800 s (a 436-call serial quick run took 1587 s). Rejected alternative: leave evaluation serial and let the full run overrun the budget. Chosen: dispatch every unique evaluate() call across a process pool (N_WORKERS=os.cpu_count(), see _evaluate_all) -- identical inputs/outputs per row, only wall-clock parallelized; benchmarked throughput ~1.0-2.0 calls/s depending on core radius.
- Reduced-cut axis LIST trimmed to the orchestrator task's explicit reduced-cut deliverable list (screening {0,1}; occupied_dot_access's 1.0 partner; S {1e2,1e4}; shell AlGaN; the spec-named injector cuts -- barrier thickness/Al fraction/alignment/growth tolerance/occupation-control uncertainty; R_s designed 1e6 on the horizontal family) plus one b_res=0.02 point, each sampling only its alternative (non-default) value(s), PLUS this fix round's restored current_pulse_width cut (I in {0.001,0.002,0.02} uA x tau_pulse in {0.01,0.1,1} ns, 9 rows, horizontal reference geometry, rectangular regime, relaxed bound, 300 K, 200 MHz) -- see DROPPED_CUT_AXES in scripts/run_nitride_nanowire.py and results.md's 'Reduced-cut coverage' section for the axes still dropped relative to the contract's fuller reduced-cut list (reservoir_access, gamma300, tau_rad0_ns, tau_cap_ps, C_parasitic_F, Rth_K_W, vertical R_s_ohm, NA/bottom_reflectivity). Rejected alternative: the contract's full one-at-a-time axis list (~1350 extra calls), which single-benchmark projections put close to or beyond the 40-minute stop threshold once added to the core grid's own runtime; the core/main grid itself is NEVER reduced. The job/cache identity (_physics_only) was also changed this round to exclude the bookkeeping sensitivity_axis/sensitivity_value labels, so a reduced-cut row whose sampled value equals the main-grid default (e.g. the screening_fraction=0.0 cut) now correctly deduplicates onto its core-row twin's real evaluate() call instead of hashing separately -- rejected alternative: leave the identity as-is and simply accept the ~32 redundant real evaluate() calls it caused.
- bound_reversal is computed as a declared post-hoc transform over paired core rows (bound_reversal_pair column) rather than via a second evaluate_strain_pair() call, to avoid doubling the core-grid evaluate() count; the device module's own per-row bound_reversal field stays its 'not_computed' sentinel on every row here.
- The 2013 CW comparison is published as a pulsed-model drive_mismatch (this evaluator has no CW/HBT path), per the contract's own instruction, never substituted for a predicted CW g2.

## Reduced-cut coverage (budget trim)

evaluate() measured 6.1-9.8 s/call serially on this machine (fsim_core.nitride_nanowire_injector.injector_feasibility's resonance search, run on every row regardless of regime, dominates); the mandated 2560-row core grid alone consumes most of the 1800 s/40-minute operational budget even parallelized across all cores (see manifest n_workers). The reduced-cut axis list is therefore trimmed to screening_fraction {0,1}, occupied_dot_access's 1.0 conservative partner, S_cm_s {1e2,1e4}, shell AlGaN, b_res 0.02, the spec-named injector cuts (barrier thickness, al_fraction, alignment_uncertainty_meV, growth_tolerance_steps, occupation-control uncertainty), and R_s_ohm=1e6 on the horizontal family only -- each sampling only its non-default alternative value(s). DROPPED entirely this run (a real, reported coverage gap against the contract's fuller reduced-cut list, never silently omitted): reservoir_access, gamma300 (linewidth), tau_rad0_ns, tau_cap_ps (capture), C_parasitic_F, Rth_K_W (both families), R_s_ohm (vertical family), NA (horizontal)/bottom_reflectivity (vertical) collection envelope.

## Limitations

Reduced cuts are one-at-a-time and moderately sampled (not exhaustive), always at the family's reference geometry (12.5/80 nm core, 2 nm height, x_in=0.40), both regimes/bounds, T in {230,300} K, both rates -- never the Cartesian product. bound_reversal_pair is a declared post-hoc transform over already-evaluated row pairs (see attach_bound_reversal's docstring), not a separate evaluate_strain_pair() call. No fitting to any Deshpande anchor occurred anywhere in this sweep. Evaluate calls are dispatched across a process pool (N_WORKERS, manifest.json) for wall-clock only -- every row's evaluate() input/output is identical to a serial run.

runtime_s=1285.9 evaluate_calls=2959 complete=True
