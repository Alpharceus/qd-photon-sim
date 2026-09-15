# Nitride nanowire sweep results

Independent predictions (design brief User answer 2): every number below is an evaluate() output or a declared, traceable reduction of one (row_id in sweep.csv); no plotted point or headline number is fitted to the Deshpande 2013/2014 held-out lifetime/g2/wavelength anchors. `family` is horizontal_as_built (Deshpande 2013/2014 as-built dispersed wire) or vertical_photonic (designed HE11 wire, headline-eligible only when single_mode AND approximation_error==0); `regime` is rectangular (100 ps electrical pulse) or deterministic_pair (idealized one-pair-per-cycle SET loading); `strain_bound` unrelaxed is the conservative_lower scenario, relaxed the headline_upper scenario -- neither is asserted a rigorous flux bound, and a reversal is reported (bound_reversal_pair) rather than hidden. Two independent, non-gating hardware screens (Coulomb-blockade `set_feasible`/`set_EC_over_kT` and resonant-tunnelling-injector `rti_feasible`) are reported alongside the idealized optical statistics; neither overwrites `g2_op`/`collected_flux_pulsed_s`.

## VERDICT lines (16: family x regime x strain_bound x rep_rate_hz)

Definition: paired_optical_pass counts an optical_pass row whose strain-bound PARTNER row (same family/regime/rep_rate_hz/core_radius_nm/height_nm/x_in/T_hs, opposite strain_bound) is present in this run -- NOT an AND requiring both strain scenarios to individually pass.

VERDICT: idealized_status=no_idealized_pass family=horizontal_as_built regime=rectangular strain_bound=unrelaxed bound_role=conservative_lower rep_rate_hz=8e+07 complete=True eligible=0 paired_optical_pass=0 hardware_qualified=0 rti_qualified=0 coverage=192/192 invalid=96 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=no_idealized_pass family=horizontal_as_built regime=rectangular strain_bound=unrelaxed bound_role=conservative_lower rep_rate_hz=2e+08 complete=True eligible=0 paired_optical_pass=0 hardware_qualified=0 rti_qualified=0 coverage=192/192 invalid=96 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=no_idealized_pass family=horizontal_as_built regime=rectangular strain_bound=relaxed bound_role=headline_upper rep_rate_hz=8e+07 complete=True eligible=168 paired_optical_pass=0 hardware_qualified=0 rti_qualified=0 coverage=192/192 invalid=24 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=no_idealized_pass family=horizontal_as_built regime=rectangular strain_bound=relaxed bound_role=headline_upper rep_rate_hz=2e+08 complete=True eligible=168 paired_optical_pass=0 hardware_qualified=0 rti_qualified=0 coverage=192/192 invalid=24 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=pass_hardware_infeasible family=horizontal_as_built regime=deterministic_pair strain_bound=unrelaxed bound_role=conservative_lower rep_rate_hz=8e+07 complete=True eligible=6 paired_optical_pass=6 hardware_qualified=0 rti_qualified=0 coverage=192/192 invalid=96 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=pass_hardware_infeasible family=horizontal_as_built regime=deterministic_pair strain_bound=unrelaxed bound_role=conservative_lower rep_rate_hz=2e+08 complete=True eligible=15 paired_optical_pass=15 hardware_qualified=0 rti_qualified=0 coverage=192/192 invalid=96 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=pass_hardware_infeasible family=horizontal_as_built regime=deterministic_pair strain_bound=relaxed bound_role=headline_upper rep_rate_hz=8e+07 complete=True eligible=168 paired_optical_pass=40 hardware_qualified=0 rti_qualified=0 coverage=192/192 invalid=24 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=pass_hardware_infeasible family=horizontal_as_built regime=deterministic_pair strain_bound=relaxed bound_role=headline_upper rep_rate_hz=2e+08 complete=True eligible=168 paired_optical_pass=12 hardware_qualified=0 rti_qualified=0 coverage=192/192 invalid=24 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=no_idealized_pass family=vertical_photonic regime=rectangular strain_bound=unrelaxed bound_role=conservative_lower rep_rate_hz=8e+07 complete=True eligible=0 paired_optical_pass=0 hardware_qualified=0 rti_qualified=0 coverage=128/128 invalid=10 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=no_idealized_pass family=vertical_photonic regime=rectangular strain_bound=unrelaxed bound_role=conservative_lower rep_rate_hz=2e+08 complete=True eligible=0 paired_optical_pass=0 hardware_qualified=0 rti_qualified=0 coverage=128/128 invalid=3 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=no_idealized_pass family=vertical_photonic regime=rectangular strain_bound=relaxed bound_role=headline_upper rep_rate_hz=8e+07 complete=True eligible=128 paired_optical_pass=0 hardware_qualified=0 rti_qualified=0 coverage=128/128 invalid=0 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=no_idealized_pass family=vertical_photonic regime=rectangular strain_bound=relaxed bound_role=headline_upper rep_rate_hz=2e+08 complete=True eligible=128 paired_optical_pass=0 hardware_qualified=0 rti_qualified=0 coverage=128/128 invalid=0 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=no_idealized_pass family=vertical_photonic regime=deterministic_pair strain_bound=unrelaxed bound_role=conservative_lower rep_rate_hz=8e+07 complete=True eligible=0 paired_optical_pass=0 hardware_qualified=0 rti_qualified=0 coverage=128/128 invalid=0 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=no_idealized_pass family=vertical_photonic regime=deterministic_pair strain_bound=unrelaxed bound_role=conservative_lower rep_rate_hz=2e+08 complete=True eligible=0 paired_optical_pass=0 hardware_qualified=0 rti_qualified=0 coverage=128/128 invalid=0 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=pass_hardware_infeasible family=vertical_photonic regime=deterministic_pair strain_bound=relaxed bound_role=headline_upper rep_rate_hz=8e+07 complete=True eligible=128 paired_optical_pass=128 hardware_qualified=0 rti_qualified=0 coverage=128/128 invalid=0 flux_floor=1000/s screening=0 access=0.05
VERDICT: idealized_status=pass_hardware_infeasible family=vertical_photonic regime=deterministic_pair strain_bound=relaxed bound_role=headline_upper rep_rate_hz=2e+08 complete=True eligible=128 paired_optical_pass=62 hardware_qualified=0 rti_qualified=0 coverage=128/128 invalid=0 flux_floor=1000/s screening=0 access=0.05

## Best passing flux per family

BEST_PASSING_FLUX family=horizontal_as_built value=2.80164e+06 row_id=CO00407 T_hs=273 screening=0 regime=deterministic_pair strain_bound=relaxed rep_rate_hz=8e+07
BEST_PASSING_FLUX_300K family=horizontal_as_built value=2.42252e+06 row_id=CO00223 T_hs=300 screening=0 regime=deterministic_pair strain_bound=relaxed rep_rate_hz=8e+07
BEST_PASSING_FLUX family=vertical_photonic value=6.15418e+06 row_id=CO01928 T_hs=230 screening=0 regime=deterministic_pair strain_bound=relaxed rep_rate_hz=2e+08
BEST_PASSING_FLUX_300K family=vertical_photonic value=4.35616e+06 row_id=CO02239 T_hs=300 screening=0 regime=deterministic_pair strain_bound=relaxed rep_rate_hz=8e+07

## Bounds table (conservative unrelaxed vs headline relaxed, reference geometry)

| family | regime | rep_rate_hz | T_hs K | unrelaxed g2 | unrelaxed flux/s | relaxed g2 | relaxed flux/s | bound_reversal | row_ids |
|---|---|---|---|---|---|---|---|---|---|
| horizontal_as_built | rectangular | 8e+07 | 230 | 0.9997970764002506 | 17.560624769775856 | 0.6903462561146173 | 3529045.4836993823 | not_comparable | CO00353,CO00355 |
| horizontal_as_built | rectangular | 8e+07 | 300 | 0.9999799925881296 | 5.051288531581465 | 0.6904497439787256 | 3539837.19708995 | not_comparable | CO00377,CO00379 |
| horizontal_as_built | rectangular | 2e+08 | 230 | 0.9997986248418573 | 43.72023750846927 | 0.6755980539994237 | 8498502.140028652 | not_comparable | CO00354,CO00356 |
| horizontal_as_built | rectangular | 2e+08 | 300 | 0.9999800961852413 | 12.592803801596693 | 0.6757369210026001 | 8525427.047682272 | not_comparable | CO00378,CO00380 |
| horizontal_as_built | deterministic_pair | 8e+07 | 230 | 0.17355371900826455 | 443.6425259744916 | 0.1793241364219459 | 4148676.236885604 | not_comparable | CO00357,CO00359 |
| horizontal_as_built | deterministic_pair | 8e+07 | 300 | 0.17355371900826433 | 138.9118555618235 | 0.1792890886443651 | 4158022.514176035 | not_comparable | CO00381,CO00383 |
| horizontal_as_built | deterministic_pair | 2e+08 | 230 | 0.17355371900826455 | 1104.741358244114 | 0.3168346309252511 | 10357246.372222234 | True | CO00358,CO00360 |
| horizontal_as_built | deterministic_pair | 2e+08 | 300 | 0.17355371900826455 | 346.39580604219583 | 0.31652841476707594 | 10380458.669112533 | not_comparable | CO00382,CO00384 |
| vertical_photonic | rectangular | 8e+07 | 230 | 0.9999999190333605 | 0.7323628286480333 | 0.6990997824990145 | 4224580.724281467 | not_comparable | CO01889,CO01891 |
| vertical_photonic | rectangular | 8e+07 | 300 | 0.999999992048697 | 0.2152402640389279 | 0.7002722427091184 | 3947403.003283787 | not_comparable | CO01913,CO01915 |
| vertical_photonic | rectangular | 2e+08 | 230 | 0.9999999190356622 | 1.8308783417608352 | 0.698095799405805 | 10532068.02444331 | not_comparable | CO01890,CO01892 |
| vertical_photonic | rectangular | 2e+08 | 300 | 0.9999999920488613 | 0.5380951210996248 | 0.6995223665864898 | 9847688.958371958 | not_comparable | CO01914,CO01916 |
| vertical_photonic | deterministic_pair | 8e+07 | 230 | 0.17355371900826455 | 15.745684531805457 | 0.1735647853023734 | 4159250.1925281966 | not_comparable | CO01893,CO01895 |
| vertical_photonic | deterministic_pair | 8e+07 | 300 | 0.17355371900826455 | 5.070544439430814 | 0.17355930780084994 | 3886469.141773661 | not_comparable | CO01917,CO01919 |
| vertical_photonic | deterministic_pair | 2e+08 | 230 | 0.17355371900826433 | 39.36362623729463 | 0.1867292085013761 | 10398117.084807843 | not_comparable | CO01894,CO01896 |
| vertical_photonic | deterministic_pair | 2e+08 | 300 | 0.17355371900826455 | 12.676243623410885 | 0.18360081169386988 | 9716146.139613492 | not_comparable | CO01918,CO01920 |

## Repetition-rate sensitivity: 80 MHz vs 200 MHz SET one_pair_valid (family-asymmetric)

vertical_photonic's 80 MHz SET row clears one_pair_valid (and hence optical_pass) at BOTH 230 K and 300 K; horizontal_as_built's 80 MHz SET row fails one_pair_valid at both temperatures too (it would need an off-grid rate near ~60 MHz to clear the loading window), passing only its own g2<0.5 and flux>=1000/s thresholds in isolation -- collected_flux_pulsed_s is the idealized/commanded flux, collected_flux_delivered_s is the RC-limited delivered flux (bullet 6).

horizontal_as_built T_hs=230K: at 80 MHz one_pair_valid=False optical_pass=False (row CO00359, g2=0.1793241364219459, flux=4148676.236885604, delivered=147448.09400815328); at 200 MHz one_pair_valid=False optical_pass=False (row CO00360, g2=0.3168346309252511, flux=10357246.372222234, delivered=368153.32836269424). The 200 MHz period leaves less time per cycle for the deterministic-pair loading window, so a one_pair_valid failure there is a loading-window/period effect, not a fit.
horizontal_as_built T_hs=300K: at 80 MHz one_pair_valid=False optical_pass=False (row CO00383, g2=0.1792890886443651, flux=4158022.514176035, delivered=154430.93229479156); at 200 MHz one_pair_valid=False optical_pass=False (row CO00384, g2=0.31652841476707594, flux=10380458.669112533, delivered=385578.9900845425). The 200 MHz period leaves less time per cycle for the deterministic-pair loading window, so a one_pair_valid failure there is a loading-window/period effect, not a fit.
vertical_photonic T_hs=230K: at 80 MHz one_pair_valid=True optical_pass=True (row CO01895, g2=0.1735647853023734, flux=4159250.1925281966, delivered=3767861.7438239465); at 200 MHz one_pair_valid=False optical_pass=False (row CO01896, g2=0.1867292085013761, flux=10398117.084807843, delivered=9419648.437798565). The 200 MHz period leaves less time per cycle for the deterministic-pair loading window, so a one_pair_valid failure there is a loading-window/period effect, not a fit.
vertical_photonic T_hs=300K: at 80 MHz one_pair_valid=True optical_pass=True (row CO01919, g2=0.17355930780084994, flux=3886469.141773661, delivered=3574053.9491063515); at 200 MHz one_pair_valid=False optical_pass=False (row CO01920, g2=0.18360081169386988, flux=9716146.139613492, delivered=8935111.502103966). The 200 MHz period leaves less time per cycle for the deterministic-pair loading window, so a one_pair_valid failure there is a loading-window/period effect, not a fit.

## Deshpande 2014 comparison (x_in=0.40, T_hs=300 K, 200 MHz, deterministic_pair)

Measured (abstract-only [V], CONDITIONS INCOMPLETE): lambda~630 nm, lifetime 1.3+/-0.3 ns, g2=0.29.

| strain_bound | lambda_nm (predicted) | tau_rad_bare_ns | tau_rad_photonic_ns | tau_total_X_ns | g2_op | collected_flux_pulsed_s | row_id |
|---|---|---|---|---|---|---|---|
| unrelaxed | 622.8766000859592 | 24.450986850391498 | 62.259743487195735 | 0.0017032738525559954 | 0.17355371900826455 | 346.39580604219583 | CO00382 |
| relaxed | 542.6613131119037 | 1.0712970826414099 | 2.7454237780315776 | 2.2463730632050725 | 0.31652841476707594 | 10380458.669112533 | CO00384 |

Untouched OLD planar 2014 comparison replay (non-gating, unchanged from the prior piece): row LI02983, g2=0.9980609340138185, flux=0.018999124213812685 vs measured g2=0.29.

## Deshpande 2013 comparison (10 K CW-equivalent, drive_mismatch)

Measured (CW electrical, 10 K, [V]): X g2 raw/corrected 0.30/0.16, XX 0.38/0.25 at 1 nA; g2-fit lifetimes X 1.1 ns, XX 0.7 ns; TRPL XX 711 ps; emission X=2.84 eV (436.56 nm). This evaluator has no CW/HBT drive path (pulsed rectangular / deterministic_pair only); the rows below are the pulsed model's own prediction at the paper's geometry/current, published as drive_mismatch -- never substituted for a predicted CW g2.

| R nm | I uA | lambda_nm (predicted) | tau_rad_bare_ns | g2_op (drive_mismatch) | flux/s | valid | row_id |
|---|---|---|---|---|---|---|---|
| 12.5 | 0.001 | nan | nan | nan (drive_mismatch) | nan | False | D1302963 |
| 12.5 | 0.002 | nan | nan | nan (drive_mismatch) | nan | False | D1302964 |
| 15.0 | 0.001 | nan | nan | nan (drive_mismatch) | nan | False | D1302965 |
| 15.0 | 0.002 | nan | nan | nan (drive_mismatch) | nan | False | D1302966 |

MISSING (not silently repaired): 4/4 2013 comparison rows are invalid this run. Reason(s): si_complex_index: lambda_nm=446.2329979828811 outside the tabulated 450.0-630.0 nm anchor range; supply an explicit n_substrate override instead; si_complex_index: lambda_nm=446.2330005160856 outside the tabulated 450.0-630.0 nm anchor range; supply an explicit n_substrate override instead; si_complex_index: lambda_nm=446.4177252000232 outside the tabulated 450.0-630.0 nm anchor range; supply an explicit n_substrate override instead; si_complex_index: lambda_nm=446.41772793972837 outside the tabulated 450.0-630.0 nm anchor range; supply an explicit n_substrate override instead.

## Opposite-endpoint anchor match (bullet 11 obligation)

2013 relaxed predicted lambda_nm=446.418 (DIAGNOSTIC ONLY: the row itself is invalid -- its bare-dot emission wavelength, already computed upstream by the levels module, is recovered from the row's own invalid_reasons message, not a repaired/assumed value; the downstream photonics step fails because fsim_core/nitride_nanowire_photonics.py's Si-substrate index table covers only 450-630 nm and 446 nm falls just below it -- see the 'Deshpande 2013 comparison' section's MISSING line and this run's STATUS notes for the interface coordination item) vs measured ~437 nm (X=2.84 eV, [V] Fig. 3c). 2014 relaxed predicted lambda_nm=542.6613131119037 / unrelaxed predicted lambda_nm=622.8766000859592 vs measured ~630 nm. The two strain-bound anchors are matched by OPPOSITE endpoints (2013 by relaxed, 2014 by whichever endpoint lands closer) and are never averaged: one of x_in transfer, disc thickness, VBO/bowing, or lateral localization carries roughly 0.3 eV of the remaining discrepancy -- stated here, not resolved.

## RC caveat (bullet 6 obligation)

As-built R_s_ohm=2.38e9 (horizontal_as_built, deterministic_pair rows): tau_RC_ns spans 1.651-26.47, delivered_step_fraction spans 0.00377-0.05877 -- this device CANNOT deliver a 100 ps step.
Designed-contact R_s_ohm=1e6 (sensitivity rows): tau_RC_ns spans 0.00111-0.001161, delivered_step_fraction spans 1-1 -- only this designed contact can deliver the 100 ps step.

## Access-1.0 lifetime cap (bullet 8 obligation)

occupied_dot_access=1.0 at core_radius_nm=12.5 (row SN02753): 1/k_surface_X_ns=0.625 ns, 1/k_surface_XX_ns=0.3125 ns (contract's stated cap: 0.625 ns X, 0.3125 ns XX) -- this caps tau_X/tau_XX regardless of any other physics; the headline default occupied_dot_access=0.05 is reported next to it, never chosen to reproduce the held-out anchors.

## E_C/kT wall (bullet 11 obligation)

Across every deterministic_pair core row with core_radius_nm>=10 nm at 230-300 K: set_EC_over_kT spans 0.1464-0.7643 (required ec_margin=10), set_feasible=True count=0/1040; rti_feasible=True count=0/992. Deterministic loading at 230-300 K fails on both the Coulomb-blockade and resonant-tunnelling-injector screens for every core_radius_nm>=10 nm priced in this tier -- by EITHER charging mechanism.

## c-plane dipole prior falsification (Composition rules bullet 9)

isotropic default (row DW02961): degree_of_linear_polarization=0.8401003547559116 (anchor: +70% axial, deshpande2013_polarization [V]). CPLANE_ONLY_DIPOLE_WEIGHTS=(0,0.5,0.5) (row DW02962): degree_of_linear_polarization=-0.9610335721927784 -- the wrong sign against the +70% anchor, confirming this sensitivity stays a named, explicitly falsified alternative and never the headline default.

## Invalid rows (all kinds)

invalid_total=497 across kinds: core=493, sensitivity=4

Breakdown by reason (leading token before ':', counted per row -- a row may carry more than one reason):

- si_complex_index: 484
- photon counting did not converge: 13

NOTE: 'si_complex_index' alone accounts for 484/497 (97 percent) of all invalid rows this run -- a single dependency-module data-table coverage gap (fsim_core/nitride_nanowire_photonics.py's Si-substrate complex-index table covers only 450-630 nm; the horizontal_as_built family's x_in=0.25 rows, whose predicted emission blue-shifts below 450 nm, hit this on every T_hs/height/regime/bound/rate combination), not a physics failure this sweep introduces or can fix within its own three scoped files. This is reported as the interface coordination item for the orchestrator/reviewer, not silently repaired by this run (no n_substrate override was invented to force these rows valid).

## Decisions (conservative choices under ambiguity)

- Measured evaluate() cost on this machine: 6.1-9.8 s/call serially (fsim_core.nitride_nanowire_injector.injector_feasibility's resonance search dominates, ~85 percent of wall time, run on EVERY row regardless of regime -- a dependency-module property, not something this piece may patch). At that cost the mandated 2560-row core grid cannot finish serially inside 1800 s (a 436-call serial quick run took 1587 s). Rejected alternative: leave evaluation serial and let the full run overrun the budget. Chosen: dispatch every unique evaluate() call across a process pool (N_WORKERS=os.cpu_count(), see _evaluate_all) -- identical inputs/outputs per row, only wall-clock parallelized; benchmarked throughput ~1.0-2.0 calls/s depending on core radius.
- Reduced-cut axis LIST trimmed to the orchestrator task's explicit reduced-cut deliverable list (screening {0,1}; occupied_dot_access's 1.0 partner; S {1e2,1e4}; shell AlGaN; the spec-named injector cuts -- barrier thickness/Al fraction/alignment/growth tolerance/occupation-control uncertainty; R_s designed 1e6 on the horizontal family) plus one b_res=0.02 point, each sampling only its alternative (non-default) value(s) -- see DROPPED_CUT_AXES in scripts/run_nitride_nanowire.py and results.md's 'Reduced-cut coverage' section for the axes this drops relative to the contract's fuller reduced-cut list (reservoir_access, gamma300, tau_rad0_ns, tau_cap_ps, C_parasitic_F, current/tau_pulse sensitivity, Rth_K_W, vertical R_s_ohm, NA/bottom_reflectivity). Rejected alternative: the contract's full one-at-a-time axis list (~1350 extra calls), which single-benchmark projections put close to or beyond the 40-minute stop threshold once added to the core grid's own runtime; the core/main grid itself is NEVER reduced.
- bound_reversal is computed as a declared post-hoc transform over paired core rows (bound_reversal_pair column) rather than via a second evaluate_strain_pair() call, to avoid doubling the core-grid evaluate() count; the device module's own per-row bound_reversal field stays its 'not_computed' sentinel on every row here.
- The 2013 CW comparison is published as a pulsed-model drive_mismatch (this evaluator has no CW/HBT path), per the contract's own instruction, never substituted for a predicted CW g2.

## Reduced-cut coverage (budget trim)

evaluate() measured 6.1-9.8 s/call serially on this machine (fsim_core.nitride_nanowire_injector.injector_feasibility's resonance search, run on every row regardless of regime, dominates); the mandated 2560-row core grid alone consumes most of the 1800 s/40-minute operational budget even parallelized across all cores (see manifest n_workers). The reduced-cut axis list is therefore trimmed to screening_fraction {0,1}, occupied_dot_access's 1.0 conservative partner, S_cm_s {1e2,1e4}, shell AlGaN, b_res 0.02, the spec-named injector cuts (barrier thickness, al_fraction, alignment_uncertainty_meV, growth_tolerance_steps, occupation-control uncertainty), and R_s_ohm=1e6 on the horizontal family only -- each sampling only its non-default alternative value(s). DROPPED entirely this run (a real, reported coverage gap against the contract's fuller reduced-cut list, never silently omitted): reservoir_access, gamma300 (linewidth), tau_rad0_ns, tau_cap_ps (capture), C_parasitic_F, current_uA/tau_pulse_ns pulse sensitivity, Rth_K_W (both families), R_s_ohm (vertical family), NA (horizontal)/bottom_reflectivity (vertical) collection envelope.

## Limitations

Reduced cuts are one-at-a-time and moderately sampled (not exhaustive), always at the family's reference geometry (12.5/80 nm core, 2 nm height, x_in=0.40), both regimes/bounds, T in {230,300} K, both rates -- never the Cartesian product. bound_reversal_pair is a declared post-hoc transform over already-evaluated row pairs (see attach_bound_reversal's docstring), not a separate evaluate_strain_pair() call. No fitting to any Deshpande anchor occurred anywhere in this sweep. Evaluate calls are dispatched across a process pool (N_WORKERS, manifest.json) for wall-clock only -- every row's evaluate() input/output is identical to a serial run.

runtime_s=1370.9 evaluate_calls=2983 complete=True
