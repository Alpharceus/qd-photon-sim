# Nitride geometry and Stark sweep results

Model-only planar predictions; every number below is an evaluate() output or a traceable reduction of one (row_id in manifest.json's plot_row_mapping / the CSVs). Unscreened (screening_fraction=0) is labelled the conventional lower bound and screened (screening_fraction=1) the upper bound -- conventional polarization-screening SCENARIOS, not a rigorous ordered bound on every metric under every polarity; crossings are not suppressed. `eligible` is the flux-floor gate only (matches scripts/run_nitride_cavity.py's `eligible` exactly); `paired_optical_pass` additionally requires g2<0.5 and, for the SET regime, one_pair_valid; `hardware_qualified` further requires set_feasible AND pair_supply_possible for the SET regime. Every headline core row carries cavity_tracking=per_T_hs: the cavity is re-tracked to that row's own operating temperature (nitride.cavity.T_track=T_hs), exactly as scripts/run_nitride_cavity.py's round-1 headline does. A separate cavity_tracking=fixed_300K sensitivity set (below) holds the cavity fixed at the card's 300 K default while T_hs varies, to show the resulting detuning as a labelled, quantified effect rather than an undisclosed artifact.

VERDICT: idealized_status=no_idealized_pass family=c_plane regime=rectangular screening=unscreened_lower complete=True eligible=0 paired_optical_pass=0 hardware_qualified=0 coverage=140/140 invalid=40 flux_floor=1000/s
VERDICT: idealized_status=no_idealized_pass family=c_plane regime=rectangular screening=screened_upper complete=True eligible=0 paired_optical_pass=0 hardware_qualified=0 coverage=140/140 invalid=0 flux_floor=1000/s
VERDICT: idealized_status=no_idealized_pass family=c_plane regime=deterministic_pair screening=unscreened_lower complete=True eligible=0 paired_optical_pass=0 hardware_qualified=0 coverage=140/140 invalid=40 flux_floor=1000/s
VERDICT: idealized_status=pass_hardware_infeasible family=c_plane regime=deterministic_pair screening=screened_upper complete=True eligible=139 paired_optical_pass=139 hardware_qualified=0 coverage=140/140 invalid=0 flux_floor=1000/s
VERDICT: idealized_status=no_idealized_pass family=a_plane regime=rectangular screening=unscreened_lower complete=True eligible=0 paired_optical_pass=0 hardware_qualified=0 coverage=140/140 invalid=0 flux_floor=1000/s
VERDICT: idealized_status=no_idealized_pass family=a_plane regime=rectangular screening=screened_upper complete=True eligible=0 paired_optical_pass=0 hardware_qualified=0 coverage=140/140 invalid=0 flux_floor=1000/s
VERDICT: idealized_status=pass_hardware_infeasible family=a_plane regime=deterministic_pair screening=unscreened_lower complete=True eligible=110 paired_optical_pass=110 hardware_qualified=0 coverage=140/140 invalid=0 flux_floor=1000/s
VERDICT: idealized_status=pass_hardware_infeasible family=a_plane regime=deterministic_pair screening=screened_upper complete=True eligible=110 paired_optical_pass=110 hardware_qualified=0 coverage=140/140 invalid=0 flux_floor=1000/s

## Supplementary orientation and geometry families (not headline coverage)

| family | regime | screening | paired_optical_pass | hardware_qualified | total | invalid |
|---|---|---|---|---|---|---|
| semipolar_11_22 | rectangular | 0 | 0 | 0 | 140 | 0 |
| semipolar_11_22 | rectangular | 0.5 | 0 | 0 | 140 | 0 |
| semipolar_11_22 | rectangular | 1 | 0 | 0 | 140 | 0 |
| semipolar_11_22 | deterministic_pair | 0 | 85 | 0 | 140 | 0 |
| semipolar_11_22 | deterministic_pair | 0.5 | 122 | 0 | 140 | 0 |
| semipolar_11_22 | deterministic_pair | 1 | 139 | 0 | 140 | 0 |
| m_plane | rectangular | 0 | 0 | 0 | 140 | 0 |
| m_plane | rectangular | 0.5 | 0 | 0 | 140 | 0 |
| m_plane | rectangular | 1 | 0 | 0 | 140 | 0 |
| m_plane | deterministic_pair | 0 | 110 | 0 | 140 | 0 |
| m_plane | deterministic_pair | 0.5 | 110 | 0 | 140 | 0 |
| m_plane | deterministic_pair | 1 | 110 | 0 | 140 | 0 |

## Maximum flux among optically-passing rows, and eligible-g2 statistics

MACHINE-CHECKABLE (verifier recomputes these from sweep.csv):
BEST_PASSING_FLUX family=c_plane value=422153 row_id=CO01350 T_hs=230 screening=1 regime=deterministic_pair
BEST_PASSING_FLUX family=a_plane value=88746.2 row_id=CO02378 T_hs=230 screening=0 regime=deterministic_pair

Core-grid paired-optical-pass rows (g2<0.5, flux>=1000/s, plus one_pair_valid for SET): 1159 of 3360; g2_min=0.173554, g2_median=0.173554.
SET g2 floor: g2 = 1-(1/(1+b_res))^2 = 0.173554 (b_res from the evaluated card's own drive.b_res). Under idealized deterministic one-pair loading, exact-one-pair counting gives zero coincidence probability by construction, so rho asymptotes to 1/(1+b_res) and g2_op to this floor nearly independent of geometry -- distinct SET g2 values observed among optical-pass rows: [0.17355]. The g2 gate therefore carries no geometry information in the SET regime; only the flux/eligibility gates discriminate geometry there.

## Fixed-anchor (cavity_tracking=fixed_300K) sensitivity

Headline core rows re-track the cavity to T_hs every row (cavity_tracking=per_T_hs). These rows instead hold the cavity fixed at the card's 300 K default while T_hs varies, at a reference geometry sampled by the core grid (H=3nm, R=10nm, screening=0), matching each row against the tracked core row at the SAME coordinates.

| orientation | regime | T_hs K | detuning_meV (fixed) | flux/s (fixed) | flux/s (tracked, matching core row) | ratio tracked/fixed | fixed row_id | tracked row_id |
|---|---|---|---|---|---|---|---|---|
| a_plane | deterministic_pair | 230 | 20.56424536789825 | 144.8264154955226 | 57181.920600857135 | 394.8 | FI07789 | CO01130 |
| a_plane | deterministic_pair | 250 | 15.098602513977433 | 272.43978739384926 | 19969.753824237505 | 73.3 | FI07790 | CO01136 |
| a_plane | deterministic_pair | 273 | 8.397437473822045 | 3814.352713434182 | 6660.280018813191 | 1.746 | FI07791 | CO01142 |
| a_plane | deterministic_pair | 300 | -2.0118662291679357e-07 | 2151.3242842275577 | 2151.3242842275577 | 1 | FI07792 | CO01148 |
| a_plane | rectangular | 230 | 20.56424536789825 | 5.298015525163914e-10 | 2.0918193829867042e-07 | 394.8 | FI07785 | CO01129 |
| a_plane | rectangular | 250 | 15.098602513977433 | 3.609719802484194e-09 | 2.6459136721420436e-07 | 73.3 | FI07786 | CO01135 |
| a_plane | rectangular | 273 | 8.397437473822045 | 1.7213278191626666e-07 | 3.0056279901492125e-07 | 1.746 | FI07787 | CO01141 |
| a_plane | rectangular | 300 | -2.0118662291679357e-07 | 3.142264791100479e-07 | 3.142264791100479e-07 | 1 | FI07788 | CO01147 |
| c_plane | deterministic_pair | 230 | 22.209883596023605 | 0.001255189394466358 | 0.8296682312254734 | 661 | FI07781 | CO01058 |
| c_plane | deterministic_pair | 250 | 16.307205251468115 | 0.003577343699706345 | 0.46962353789324307 | 131.3 | FI07782 | CO01064 |
| c_plane | deterministic_pair | 273 | 9.07031298689187 | 0.1369233171876589 | 0.25775337506045704 | 1.882 | FI07783 | CO01070 |
| c_plane | deterministic_pair | 300 | -2.1593571375433385e-07 | 0.13742222076335409 | 0.13742222076335409 | 1 | FI07784 | CO01076 |
| c_plane | rectangular | 230 | 22.209883596023605 | 3.880173957537539e-05 | 0.02564746663466914 | 661 | FI07777 | CO01057 |
| c_plane | rectangular | 250 | 16.307205251468115 | 0.00011058353261280492 | 0.014517082852426663 | 131.3 | FI07778 | CO01063 |
| c_plane | rectangular | 273 | 9.07031298689187 | 0.004232320588079984 | 0.007967221511740175 | 1.882 | FI07779 | CO01069 |
| c_plane | rectangular | 300 | -2.1593571375433385e-07 | 0.004247138693754914 | 0.004247138693754914 | 1 | FI07780 | CO01075 |

## SET hardware feasibility (item 2)

hardware_qualified=0 in every headline VERDICT line because set_EC_over_kT is 1.173-1.53 across all 1620 SET core rows at the card's assumed island radius ([5.0] nm; allowed max 0.5863-0.7648 nm), versus the required margin ec_margin=10 -- the sole failing criterion. R_T/R_Q (38.74-38.74) and f_max (1.892e+10-1.892e+10 Hz) both pass. The island-radius sensitivity rows below carry set_EC_over_kT across an explicit radius sweep (0.5/1/5 nm) to show how strongly this single criterion depends on the assumed island size.

| baseline | island radius nm | T_hs K | E_C/kT | allowed max radius nm | R_T/R_Q | f_max Hz | row_id |
|---|---|---|---|---|---|---|---|
| REF_pass_a_plane_SET | 0.5 | 300 | 11.726382211781894 | 0.5863191105890947 | 38.74045864977526 | 189211616679.18253 | SE07774 |
| REF_pass_a_plane_SET | 1 | 300 | 5.863191105890947 | 0.5863191105890947 | 38.74045864977526 | 94605808339.59126 | SE07775 |
| REF_pass_a_plane_SET | 5 | 300 | 1.1726382211781896 | 0.5863191105890947 | 38.74045864977526 | 18921161667.91825 | SE07776 |
| REF_pass_c_plane_screened_SET | 0.5 | 300 | 11.726382211781894 | 0.5863191105890947 | 38.74045864977526 | 189211616679.18253 | SE07727 |
| REF_pass_c_plane_screened_SET | 1 | 300 | 5.863191105890947 | 0.5863191105890947 | 38.74045864977526 | 94605808339.59126 | SE07728 |
| REF_pass_c_plane_screened_SET | 5 | 300 | 1.1726382211781896 | 0.5863191105890947 | 38.74045864977526 | 18921161667.91825 | SE07729 |
| REF_unscreened_rectangular | 0.5 | 300 | 11.726382211781894 | 0.5863191105890947 | 38.74045864977526 | 189211616679.18253 | SE07680 |
| REF_unscreened_rectangular | 1 | 300 | 5.863191105890947 | 0.5863191105890947 | 38.74045864977526 | 94605808339.59126 | SE07681 |
| REF_unscreened_rectangular | 5 | 300 | 1.1726382211781896 | 0.5863191105890947 | 38.74045864977526 | 18921161667.91825 | SE07682 |

## Per-card reservoir choice

| card_file | orientation | geometry_type | reservoir_kind (sample) |
|---|---|---|---|
| nitride-cavity-pulse-design.yaml | c_plane | isolated_dot | gan_barrier |
| nitride-cavity-set-design.yaml | c_plane | isolated_dot | gan_barrier |
| nitride-cavity-pulse-design.yaml | semipolar_11_22 | isolated_dot | gan_barrier |
| nitride-cavity-set-design.yaml | semipolar_11_22 | isolated_dot | gan_barrier |
| nitride-cavity-pulse-design.yaml | m_plane | isolated_dot | gan_barrier |
| nitride-cavity-set-design.yaml | m_plane | isolated_dot | gan_barrier |
| nitride-nonpolar-pulse-design.yaml | a_plane | isolated_dot | gan_barrier |
| nitride-nonpolar-set-design.yaml | a_plane | isolated_dot | gan_barrier |
| nitride-qw-fluctuation-pulse-design.yaml | c_plane | qw_fluctuation | gan_barrier |
| nitride-qw-fluctuation-set-design.yaml | c_plane | qw_fluctuation | gan_barrier |
| nitride-qw-fluctuation-pulse-design.yaml | a_plane | qw_fluctuation | ingan_qw |
| nitride-qw-fluctuation-set-design.yaml | a_plane | qw_fluctuation | ingan_qw |
| nitride-deshpande2014-comparison-design.yaml | c_plane | isolated_dot | gan_barrier |

## Coverage and geometry change versus round 1 (scripts/run_nitride_cavity.py)

Round 1's headline grid was c-plane only, radius fixed at 10 nm (height in {1,2,3,4,5} nm, x_in in {0.15,0.25,0.4}, Q in {500,2000,10000}, current in {0.002,0.02,0.2} uA), no orientation axis, and no Stark bias/current traces. This piece adds a 4-orientation x 5-radius x 7-height x 4-T_hs x 3-screening x 2-regime core grid (3360 rows this run), lens/truncated_cone shape (864 rows) and QW-fluctuation (288 rows) geometry supplements, and bias/current Stark diagnostic traces (2040+1080 rows) with derivative-refinement convergence checks and screening-compatibility fits. Round 1's own results.md and artifacts under out/nitride_cavity/ are unchanged and are not reinterpreted here (out of scope for this piece). Comparability (fix round 3, item 7): round 1's headline used cavity_tracking=per_T_hs (T_track=T_hs) and screening_fraction=0 exclusively. This round's headline now ALSO uses cavity_tracking=per_T_hs for every core row (item 1), so the family=c_plane screening=unscreened_lower VERDICT line above IS on comparable tracking-convention and screening terms with round 1's headline; the screening=screened_upper line is a round-2-only addition with no round-1 counterpart and is NOT comparable to round 1 on screening.

## Radius-degeneracy disclosure

Checked directly against sweep.csv: 648 of 648 (orientation, height, T_hs, screening, regime) groups have IDENTICAL overlap_sq across every sampled radius (5,10,15,20,30 nm) -- this planar model's overlap_sq and tau_rad_bare_ns do not depend on radius_nm; the radius axis enters only through E_a_meV (escape barrier, ~11 meV over 5-30 nm) and the escape/counting prefactors. Most core rows at fixed (orientation,height,T_hs,screening,regime) are therefore near-duplicates in E_X/overlap/tau_rad_bare, differing materially only through the escape-limited flux and validity at large radius.


## Shape supplement (mapping sensitivity, not a demonstrated shape accuracy)
864 rows (lens / truncated_cone, native and full-height variants). The numerical spread across shape variants at matched height/radius is a mapping-sensitivity diagnostic; systematic shape accuracy remains unquantified (no 3-D solver).

## QW-fluctuation supplement
288 rows; diode.wl_thickness_nm/x_in and dot.wl_thickness_nm/x_in are matched per card contract; d_i_nm=24+w, d_active_nm=height [A]. No lens-in-QW cases were run (unsupported combination).

## Literature comparisons (non-gating)

Wang et al., Sci. Rep. 7, 12089 (2017): uncapped a-plane AFM ~7 nm/~35 nm dots, E_X=2.54 eV [V] at 220 K, 19.0+/-0.4 meV linewidth [V], raw/corrected g2=0.47/0.21 [V], optical excitation (76 MHz, 1 ps, 800 nm two-photon) -- NOT an electrical SPS or a 300 K validation. Composition transfer [A]: this model's cards use x_in=0.25 (In0.25Ga0.75N) at the same nominal geometry; Wang's paper does not report a composition, so the model's own default is retained rather than fitted. Model rows (comparison-only, excluded from headline coverage):

- row LI04513: T_hs=220 K, x_in=0.25, E_X_eV=2.790254882472666, g2=0.6792796752386774, valid=True (comparison_only_reported_temperature)
- row LI04514: T_hs=300 K, x_in=0.25, E_X_eV=2.7644024701218513, g2=0.9822939349276038, valid=True (comparison_only_model_headline_temperature)
- Discrepancy at 220 K (never fitted to force agreement): model E_X=2.79 eV vs measured 2.54 eV [V]; model g2=0.6792796752386774 vs measured raw/corrected 0.47/0.21 [V]. The gap reflects the assumed x_in=0.25 transfer, uncapped-vs-capped surface treatment and this model's idealized single-band overlap, none of which were tuned to close it.
- Deshpande et al., APL 105, 141109 (2014) [V abstract-only; CONDITIONS INCOMPLETE], replay row LI04515: measured g2=0.29, model g2=0.9980609340138185.

Zhang et al., APL 108, 153102 (2016), Fig. 5: -10 meV/V below 2 V [V], measured at T=10 K on a 3 nm In0.15Ga0.85N QW [V Fig. 5 / device section] -- both the temperature and the composition (x=0.15 vs this model's x_in=0.25 curves) are disclosed differences, non-gating (different device, temperature and excitation from this model's 230-300 K planar-dot cards); shown only as a labelled guide anchored to a model row on the Stark figures, never a fitted target.

## Screening-compatibility table

The [-12,-8] meV/V window used when --slope-range is not supplied is EXPLICITLY illustrative around Zhang's approximate value, never a measured interval; the same --bias-window applies to any user-supplied slope. screening_compatibility is now called ONCE per (orientation,height,polarity,T_hs) geometry group with ALL THREE screening hypotheses together (fix round 3, item 3), so a genuine degeneracy across screening_fraction is detected instead of being hidden by testing each hypothesis in isolation. Nonpolar (m-plane/a-plane) rows are marked screening_unidentifiable ONLY when the shared curve already matches the window (polarization_factor=0 makes screening_fraction physically inert there, so an ACCEPTED nonpolar slope cannot distinguish a screening hypothesis); a nonpolar slope the window rejects stays incompatible.

120 raw compatibility-fit rows reduce to 30 DISTINCT hypotheses (orientation, height, polarity, screening) after collapsing regime/T_hs replicates: 2 compatible, 28 incompatible, 0 screening_unidentifiable.

| orientation | height_nm | polarity | screening | fitted_slope meV/V | lifetime_ratio(s0/s1) | row_id |
|---|---|---|---|---|---|---|
| c_plane | 2.0 | 1 | 0.0 | -8.541033800799307 | 12.88422779657124 | CT07841 |
| c_plane | 3.0 | 1 | 0.5 | -10.57221037117655 | 183.97227797648063 | CT07866 |

Height-screening degeneracy at orientation=c_plane, polarity=1: [(2.0, 0.0), (3.0, 0.5)] all fit the same illustrative window -- the window does NOT identify screening_fraction from height alone. The discriminating observable is the bias-resolved lifetime (see lifetime_ratio_s0_over_s1 column above, computed from stark_bias.csv's tau_rad_bare_ns at the unscreened vs screened hypotheses).

## One-at-a-time sensitivities at multiple baselines (item 5)

The original REF baseline (rectangular, unscreened) never passes the optical gate; REF_pass_c_plane_screened_SET and REF_pass_a_plane_SET (H=3nm, R=10nm, T_hs=300K, deterministic_pair regime) DO pass, so every [A] assumption below is also bounded at a configuration that produces a verdict.

### REF_unscreened_rectangular (optical_pass at baseline: False)

| axis | value | T_hs K | g2 | flux/s | optical_pass | row_id |
|---|---|---|---|---|---|---|
| Q_purcell | 10000.0 | 230 | 0.998874238514352 | 0.005690607601621742 | False | SE07668 |
| Q_purcell | 10000.0 | 300 | 0.999668584358809 | 0.0008886884921624619 | False | SE07669 |
| Q_purcell | 167.0 | 230 | 0.9991729883153779 | 0.10230272846147936 | False | SE07664 |
| Q_purcell | 167.0 | 300 | 0.9997815927440364 | 0.02962426339338993 | False | SE07665 |
| Q_purcell | 500.0 | 230 | 0.9989953799307202 | 0.07201897999749715 | False | SE07666 |
| Q_purcell | 500.0 | 230 | 0.9992526492989291 | 0.014230183139320951 | False | SE07670 |
| Q_purcell | 500.0 | 300 | 0.9997145796704188 | 0.014456284713042571 | False | SE07667 |
| Q_purcell | 500.0 | 300 | 0.9998118256421988 | 0.004651655297391359 | False | SE07671 |
| b_res | 0.0 | 230 | 0.9986704369648789 | 0.02564746663466914 | False | SE07660 |
| b_res | 0.0 | 300 | 0.9996118251057948 | 0.004247138693754914 | False | SE07661 |
| b_res | 0.5 | 230 | 0.9994090808002184 | 0.02564746663466914 | False | SE07662 |
| b_res | 0.5 | 300 | 0.9998274766069121 | 0.004247138693754914 | False | SE07663 |
| background_tau_ns | 0.0 | 230 | 0.9989011864102583 | 0.02564746663466914 | False | SE07656 |
| background_tau_ns | 0.0 | 300 | 0.9996791936846767 | 0.004247138693754914 | False | SE07657 |
| background_tau_ns | 1.0 | 230 | 0.9989011864102583 | 0.02564746663466914 | False | SE07658 |
| background_tau_ns | 1.0 | 300 | 0.9996791936846767 | 0.004247138693754914 | False | SE07659 |
| detuning_offset_meV | 0.0 | 230 | 0.9989011864102583 | 0.02564746663466914 | False | SE07672 |
| detuning_offset_meV | 0.0 | 300 | 0.9996791936846767 | 0.004247138693754914 | False | SE07673 |
| detuning_offset_meV | 10.0 | 230 | 0.9988929950189479 | 0.0007759227944418132 | False | SE07674 |
| detuning_offset_meV | 10.0 | 300 | 0.9994951045143491 | 0.002680663150882631 | False | SE07675 |
| island_radius_nm | 0.5 | 300 | 0.17355371900826455 | 0.13742222076335409 | False | SE07680 |
| island_radius_nm | 1.0 | 300 | 0.17355371900826455 | 0.13742222076335409 | False | SE07681 |
| island_radius_nm | 5.0 | 300 | 0.17355371900826455 | 0.13742222076335409 | False | SE07682 |
| k_nr_ns | 0.0 | 230 | 0.9989011864102583 | 0.02564746663466914 | False | SE07652 |
| k_nr_ns | 0.0 | 300 | 0.9996791936846767 | 0.004247138693754914 | False | SE07653 |
| k_nr_ns | 1.0 | 230 | 0.998901280159109 | 0.025645276589994107 | False | SE07654 |
| k_nr_ns | 1.0 | 300 | 0.9996792109019386 | 0.0042469106385417165 | False | SE07655 |
| semipolar_factor | 0.1 | 230 | 0.6061277704722288 | 5.733868746510371e-06 | False | SE07636 |
| semipolar_factor | 0.1 | 300 | 0.9314341260827127 | 7.598092260246392e-06 | False | SE07637 |
| semipolar_factor | 0.3 | 230 | 0.9272838867016583 | 1.116450709650972e-05 | False | SE07638 |
| semipolar_factor | 0.3 | 300 | 0.991598778800592 | 6.838818575934152e-06 | False | SE07639 |
| strain_fraction | 0.0 | 230 | 0.5998581126239653 | 1.6852257894762432 | False | SE07644 |
| strain_fraction | 0.0 | 300 | 0.7944633953659482 | 0.6948968499953352 | False | SE07645 |
| strain_fraction | 0.5 | 230 | 0.9360753047771412 | 0.38764442609137595 | False | SE07646 |
| strain_fraction | 0.5 | 300 | 0.9924924046426596 | 0.016128383823171638 | False | SE07647 |
| tau_cap_density_convention | False | 230 | 0.9998900328830732 | 0.0023549511918210496 | False | SE07676 |
| tau_cap_density_convention | False | 300 | 0.9999679041937176 | 0.00038996687724695585 | False | SE07677 |
| tau_cap_density_convention | True | 230 | 0.9989011854677752 | 0.02354907714289014 | False | SE07678 |
| tau_cap_density_convention | True | 300 | 0.9996791940083959 | 0.0038996474417952887 | False | SE07679 |
| tau_rad0_ns | 0.5 | 230 | 0.998901188459612 | 0.0512948488517463 | False | SE07648 |
| tau_rad0_ns | 0.5 | 300 | 0.9996791939204358 | 0.008494272543473626 | False | SE07649 |
| tau_rad0_ns | 2.0 | 230 | 0.9989011853858228 | 0.012823743869555839 | False | SE07650 |
| tau_rad0_ns | 2.0 | 300 | 0.9996791935669387 | 0.0021235699523821355 | False | SE07651 |
| x_in | 0.15 | 230 | 0.999509320280985 | 0.023133881803611527 | False | SE07640 |
| x_in | 0.15 | 300 | 0.999825218680176 | 0.004527042996645048 | False | SE07641 |
| x_in | 0.4 | 230 | 0.9889723241767102 | 0.04769243032115655 | False | SE07642 |
| x_in | 0.4 | 300 | 0.9981685901672467 | 0.004750533928600966 | False | SE07643 |

### REF_pass_c_plane_screened_SET (optical_pass at baseline: True)

| axis | value | T_hs K | g2 | flux/s | optical_pass | row_id |
|---|---|---|---|---|---|---|
| Q_purcell | 10000.0 | 230 | 0.17355371900826455 | 88406.31439103976 | True | SE07715 |
| Q_purcell | 10000.0 | 300 | 0.17355371900826433 | 5437.181532723568 | True | SE07716 |
| Q_purcell | 167.0 | 230 | 0.17355371900826433 | 1651115.4461871076 | True | SE07711 |
| Q_purcell | 167.0 | 300 | 0.17355371900826455 | 155585.25131821044 | True | SE07712 |
| Q_purcell | 500.0 | 230 | 0.17355371900826455 | 1137032.4352100824 | True | SE07713 |
| Q_purcell | 500.0 | 230 | 0.17355371900826433 | 367917.5667745471 | True | SE07717 |
| Q_purcell | 500.0 | 300 | 0.17355371900826455 | 83427.84112599432 | True | SE07714 |
| Q_purcell | 500.0 | 300 | 0.17355371900826455 | 24115.392004527985 | True | SE07718 |
| b_res | 0.0 | 230 | 0.0 | 400171.62811885105 | True | SE07707 |
| b_res | 0.0 | 300 | 0.0 | 25642.788985855634 | True | SE07708 |
| b_res | 0.5 | 230 | 0.5555555555555555 | 400171.62811885105 | False | SE07709 |
| b_res | 0.5 | 300 | 0.5555555555555556 | 25642.788985855634 | False | SE07710 |
| background_tau_ns | 0.0 | 230 | 0.17355371900826455 | 400171.62811885105 | True | SE07703 |
| background_tau_ns | 0.0 | 300 | 0.17355371900826433 | 25642.788985855634 | True | SE07704 |
| background_tau_ns | 1.0 | 230 | 0.17355371900826455 | 400171.62811885105 | True | SE07705 |
| background_tau_ns | 1.0 | 300 | 0.17355371900826433 | 25642.788985855634 | True | SE07706 |
| detuning_offset_meV | 0.0 | 230 | 0.17355371900826455 | 400171.62811885105 | True | SE07719 |
| detuning_offset_meV | 0.0 | 300 | 0.17355371900826433 | 25642.788985855634 | True | SE07720 |
| detuning_offset_meV | 10.0 | 230 | 0.17355371900826455 | 22310.87181166967 | True | SE07721 |
| detuning_offset_meV | 10.0 | 300 | 0.17355371900826455 | 16296.641000198546 | True | SE07722 |
| island_radius_nm | 0.5 | 300 | 0.17355371900826433 | 25642.788985855634 | True | SE07727 |
| island_radius_nm | 1.0 | 300 | 0.17355371900826433 | 25642.788985855634 | True | SE07728 |
| island_radius_nm | 5.0 | 300 | 0.17355371900826433 | 25642.788985855634 | True | SE07729 |
| k_nr_ns | 0.0 | 230 | 0.17355371900826455 | 400171.62811885105 | True | SE07699 |
| k_nr_ns | 0.0 | 300 | 0.17355371900826433 | 25642.788985855634 | True | SE07700 |
| k_nr_ns | 1.0 | 230 | 0.17355371900826433 | 365371.0019426988 | True | SE07701 |
| k_nr_ns | 1.0 | 300 | 0.17355371900826455 | 25073.25855486469 | True | SE07702 |
| semipolar_factor | 0.1 | 230 | 0.17355371900826455 | 400171.62811885105 | True | SE07683 |
| semipolar_factor | 0.1 | 300 | 0.17355371900826433 | 25642.788985855634 | True | SE07684 |
| semipolar_factor | 0.3 | 230 | 0.17355371900826455 | 400171.62811885105 | True | SE07685 |
| semipolar_factor | 0.3 | 300 | 0.17355371900826433 | 25642.788985855634 | True | SE07686 |
| strain_fraction | 0.0 | 230 | 0.17355371900826455 | 584724.0790347588 | True | SE07691 |
| strain_fraction | 0.0 | 300 | 0.17355371900826455 | 171710.28413773407 | True | SE07692 |
| strain_fraction | 0.5 | 230 | 0.17355371900826433 | 563754.0177845515 | True | SE07693 |
| strain_fraction | 0.5 | 300 | 0.17355371900826455 | 80157.36992368038 | True | SE07694 |
| tau_cap_density_convention | False | 230 | 0.17355371900826455 | 91168.62376676891 | True | SE07723 |
| tau_cap_density_convention | False | 300 | 0.17355371900826455 | 2769.4993849051384 | True | SE07724 |
| tau_cap_density_convention | True | 230 | 0.17355371900826455 | 400171.62811885093 | True | SE07725 |
| tau_cap_density_convention | True | 300 | 0.17355371900826433 | 25642.788985855634 | True | SE07726 |
| tau_rad0_ns | 0.5 | 230 | 0.17355371900826455 | 493002.8923843233 | True | SE07695 |
| tau_rad0_ns | 0.5 | 300 | 0.17355371900826455 | 47384.26262052851 | True | SE07696 |
| tau_rad0_ns | 2.0 | 230 | 0.17355371900826455 | 290696.65851941006 | True | SE07697 |
| tau_rad0_ns | 2.0 | 300 | 0.17355371900826433 | 13371.871427560782 | True | SE07698 |
| x_in | 0.15 | 230 | 0.17355371900826433 | 17732.554550375757 | True | SE07687 |
| x_in | 0.15 | 300 | 0.17355371900826455 | 1271.7733239861022 | True | SE07688 |
| x_in | 0.4 | 230 | 0.17355371900826433 | 567547.1247777955 | True | SE07689 |
| x_in | 0.4 | 300 | 0.17355371900826455 | 258752.24049906997 | True | SE07690 |

### REF_pass_a_plane_SET (optical_pass at baseline: True)

| axis | value | T_hs K | g2 | flux/s | optical_pass | row_id |
|---|---|---|---|---|---|---|
| Q_purcell | 10000.0 | 230 | 0.17355371900826455 | 13053.236789716952 | True | SE07762 |
| Q_purcell | 10000.0 | 300 | 0.17355371900826455 | 457.2997096285444 | False | SE07763 |
| Q_purcell | 167.0 | 230 | 0.17355371900826455 | 180878.68866609386 | True | SE07758 |
| Q_purcell | 167.0 | 300 | 0.17355371900826433 | 12744.002351576071 | True | SE07759 |
| Q_purcell | 500.0 | 230 | 0.17355371900826433 | 146965.27256956187 | True | SE07760 |
| Q_purcell | 500.0 | 230 | 0.17355371900826455 | 26091.505176731993 | True | SE07764 |
| Q_purcell | 500.0 | 300 | 0.17355371900826455 | 6941.9125037904205 | True | SE07761 |
| Q_purcell | 500.0 | 300 | 0.17355371900826455 | 1891.6709013662366 | True | SE07765 |
| b_res | 0.0 | 230 | 0.0 | 57181.920600857135 | True | SE07754 |
| b_res | 0.0 | 300 | 0.0 | 2151.3242842275577 | True | SE07755 |
| b_res | 0.5 | 230 | 0.5555555555555555 | 57181.920600857135 | False | SE07756 |
| b_res | 0.5 | 300 | 0.5555555555555556 | 2151.3242842275577 | False | SE07757 |
| background_tau_ns | 0.0 | 230 | 0.17355371900826455 | 57181.920600857135 | True | SE07750 |
| background_tau_ns | 0.0 | 300 | 0.17355371900826455 | 2151.3242842275577 | True | SE07751 |
| background_tau_ns | 1.0 | 230 | 0.17355371900826455 | 57181.920600857135 | True | SE07752 |
| background_tau_ns | 1.0 | 300 | 0.17355371900826455 | 2151.3242842275577 | True | SE07753 |
| detuning_offset_meV | 0.0 | 230 | 0.17355371900826455 | 57181.920600857135 | True | SE07766 |
| detuning_offset_meV | 0.0 | 300 | 0.17355371900826455 | 2151.3242842275577 | True | SE07767 |
| detuning_offset_meV | 10.0 | 230 | 0.17355371900826455 | 2348.298741773773 | True | SE07768 |
| detuning_offset_meV | 10.0 | 300 | 0.17355371900826455 | 1348.786846431506 | True | SE07769 |
| island_radius_nm | 0.5 | 300 | 0.17355371900826455 | 2151.3242842275577 | True | SE07774 |
| island_radius_nm | 1.0 | 300 | 0.17355371900826455 | 2151.3242842275577 | True | SE07775 |
| island_radius_nm | 5.0 | 300 | 0.17355371900826455 | 2151.3242842275577 | True | SE07776 |
| k_nr_ns | 0.0 | 230 | 0.17355371900826455 | 57181.920600857135 | True | SE07746 |
| k_nr_ns | 0.0 | 300 | 0.17355371900826455 | 2151.3242842275577 | True | SE07747 |
| k_nr_ns | 1.0 | 230 | 0.17355371900826455 | 56470.59121734534 | True | SE07748 |
| k_nr_ns | 1.0 | 300 | 0.17355371900826455 | 2147.5359524866785 | True | SE07749 |
| semipolar_factor | 0.1 | 230 | 0.17355371900826455 | 262119.7530452213 | True | SE07730 |
| semipolar_factor | 0.1 | 300 | 0.17355371900826455 | 13183.392890562816 | True | SE07731 |
| semipolar_factor | 0.3 | 230 | 0.17355371900826433 | 6993.311170239615 | True | SE07732 |
| semipolar_factor | 0.3 | 300 | 0.17355371900826433 | 438.5540861159112 | False | SE07733 |
| strain_fraction | 0.0 | 230 | 0.17355371900826433 | 465511.7827716468 | True | SE07738 |
| strain_fraction | 0.0 | 300 | 0.17355371900826433 | 27750.412451603923 | True | SE07739 |
| strain_fraction | 0.5 | 230 | 0.17355371900826433 | 227819.75707326696 | True | SE07740 |
| strain_fraction | 0.5 | 300 | 0.17355371900826455 | 7913.174978671014 | True | SE07741 |
| tau_cap_density_convention | False | 230 | 0.17355371900826455 | 6211.576496416844 | True | SE07770 |
| tau_cap_density_convention | False | 300 | 0.17355371900826455 | 216.46488105836696 | False | SE07771 |
| tau_cap_density_convention | True | 230 | 0.17355371900826455 | 57181.920600857135 | True | SE07772 |
| tau_cap_density_convention | True | 300 | 0.17355371900826455 | 2151.3242842275577 | True | SE07773 |
| tau_rad0_ns | 0.5 | 230 | 0.17355371900826433 | 105089.15944389762 | True | SE07742 |
| tau_rad0_ns | 0.5 | 300 | 0.17355371900826455 | 4273.420678012138 | True | SE07743 |
| tau_rad0_ns | 2.0 | 230 | 0.17355371900826455 | 29910.857051552517 | True | SE07744 |
| tau_rad0_ns | 2.0 | 300 | 0.17355371900826433 | 1079.3532388476697 | True | SE07745 |
| x_in | 0.15 | 230 | 0.17355371900826455 | 1150.3280722100662 | True | SE07734 |
| x_in | 0.15 | 300 | 0.17355371900826455 | 104.59868607157303 | False | SE07735 |
| x_in | 0.4 | 230 | 0.17355371900826455 | 568554.1794343806 | True | SE07736 |
| x_in | 0.4 | 300 | 0.17355371900826455 | 165679.72860721534 | True | SE07737 |

## Invalid rows (all kinds)

invalid_total=311 across kinds: core=120, shape=72, stark_bias=119.

Reasons: photon counting did not converge=119, unbound dot: nonphysical E_X=64, unbound dot: overlap_unresolved (QCSE-separated pair)=128

c-plane screening=0.5 invalid core rows: 40 (included in the core count above).


## Derivative-refinement convergence checks

| check | orientation | height_nm | screening | polarity | T_hs | anchor V_j | coarse meV/V | fine (half-step) meV/V | abs diff | rel diff | converged | anchor row_id | probe row_ids |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 0 | c_plane | 3.0 | 0.0 | 1 | 300.0 | 1 | -14.32 | -14.31 | 0.01005 | 0.0007013 | True | ST05354 | ['RF07793', 'RF07794'] |
| 1 | c_plane | 3.0 | 1.0 | 1 | 300.0 | 1 | 2.925 | 2.925 | 0.0003173 | 0.0001085 | True | ST05490 | ['RF07795', 'RF07796'] |
| 2 | a_plane | 3.0 | 0.0 | 1 | 300.0 | 1 | 1.913 | 1.914 | 0.0002305 | 0.0001205 | True | ST05558 | ['RF07797', 'RF07798'] |
| 3 | c_plane | 3.0 | 0.5 | 1 | 300.0 | 1 | -10.6 | -10.6 | 0.007743 | 0.0007302 | True | ST05456 | ['RF07799', 'RF07800'] |
| 4 | c_plane | 1.0 | 0.0 | 1 | 300.0 | 1 | -3.568 | -3.565 | 0.002579 | 0.000723 | True | ST04538 | ['RF07801', 'RF07802'] |
| 5 | c_plane | 5.0 | 1.0 | 1 | 300.0 | 1.4 | 10.15 | 10.15 | 0.0007964 | 7.847e-05 | True | ST05900 | ['RF07803', 'RF07804'] |
| 6 | a_plane | 2.0 | 0.5 | 1 | 300.0 | 1 | 0.8297 | 0.8298 | 0.000102 | 0.0001229 | True | ST05218 | ['RF07805', 'RF07806'] |
| 7 | c_plane | 3.0 | 0.0 | 1 | 230.0 | 1 | -14.06 | -14.05 | 0.009253 | 0.0006581 | True | ST05337 | ['RF07807', 'RF07808'] |
| 8 | a_plane | 3.0 | 1.0 | 1 | 300.0 | 1 | 1.913 | 1.914 | 0.0002305 | 0.0001205 | True | ST05728 | ['RF07809', 'RF07810'] |
| 9 | c_plane | 3.0 | 1.0 | -1 | 300.0 | 1 | 2.925 | 2.925 | 0.0003173 | 0.0001085 | True | ST06306 | ['RF07811', 'RF07812'] |
| 10 | c_plane | 5.0 | 0.0 | 1 | 300.0 | 1.6 | -30.28 | -30.25 | 0.038 | 0.001255 | True | ST05799 | ['RF07813', 'RF07814'] |
| 11 | a_plane | 1.0 | 1.0 | 1 | 300.0 | 1 | 0.3747 | 0.3747 | 5.384e-05 | 0.0001437 | True | ST04878 | ['RF07815', 'RF07816'] |

## Limitations
Nonpolar strain, valence-band ordering and FSS are not modelled; shape mapping error is unquantified; finite-dot lateral fields, field-assisted escape and injection-dependent screening are not solved (screening_fraction is fixed along every current trace). m-plane and a-plane rows coincide under this scalar model (both polarization_factor=0) -- a model limitation, not independent evidence for either. SET hardware feasibility (island charging energy / RC bandwidth) is reported independently of the idealized optical pass and never substitutes for it. No field-free lifetime or screening fraction was ever fitted to force agreement with Wang or Zhang.

runtime_s=1078.1 evaluate_calls=7816 complete=True
