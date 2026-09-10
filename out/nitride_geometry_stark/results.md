# Nitride geometry and Stark sweep results

Model-only planar predictions; every number below is an evaluate() output or a traceable reduction of one (row_id in manifest.json's plot_row_mapping / the CSVs). Unscreened (screening_fraction=0) is labelled the conventional lower bound and screened (screening_fraction=1) the upper bound -- conventional polarization-screening SCENARIOS, not a rigorous ordered bound on every metric under every polarity; crossings are not suppressed. `eligible` is the flux-floor gate only (matches scripts/run_nitride_cavity.py's `eligible` exactly); `paired_optical_pass` additionally requires g2<0.5 and, for the SET regime, one_pair_valid; `hardware_qualified` further requires set_feasible AND pair_supply_possible for the SET regime.

VERDICT: idealized_status=no_idealized_pass family=c_plane regime=rectangular screening=unscreened_lower complete=True eligible=0 paired_optical_pass=0 hardware_qualified=0 coverage=140/140 invalid=20 flux_floor=1000/s
VERDICT: idealized_status=no_idealized_pass family=c_plane regime=rectangular screening=screened_upper complete=True eligible=0 paired_optical_pass=0 hardware_qualified=0 coverage=140/140 invalid=0 flux_floor=1000/s
VERDICT: idealized_status=no_idealized_pass family=c_plane regime=deterministic_pair screening=unscreened_lower complete=True eligible=0 paired_optical_pass=0 hardware_qualified=0 coverage=140/140 invalid=20 flux_floor=1000/s
VERDICT: idealized_status=pass_hardware_infeasible family=c_plane regime=deterministic_pair screening=screened_upper complete=True eligible=112 paired_optical_pass=112 hardware_qualified=0 coverage=140/140 invalid=0 flux_floor=1000/s
VERDICT: idealized_status=no_idealized_pass family=a_plane regime=rectangular screening=unscreened_lower complete=True eligible=0 paired_optical_pass=0 hardware_qualified=0 coverage=140/140 invalid=0 flux_floor=1000/s
VERDICT: idealized_status=no_idealized_pass family=a_plane regime=rectangular screening=screened_upper complete=True eligible=0 paired_optical_pass=0 hardware_qualified=0 coverage=140/140 invalid=0 flux_floor=1000/s
VERDICT: idealized_status=pass_hardware_infeasible family=a_plane regime=deterministic_pair screening=unscreened_lower complete=True eligible=45 paired_optical_pass=45 hardware_qualified=0 coverage=140/140 invalid=0 flux_floor=1000/s
VERDICT: idealized_status=pass_hardware_infeasible family=a_plane regime=deterministic_pair screening=screened_upper complete=True eligible=45 paired_optical_pass=45 hardware_qualified=0 coverage=140/140 invalid=0 flux_floor=1000/s

## Supplementary orientation and geometry families (not headline coverage)

| family | regime | screening | paired_optical_pass | hardware_qualified | total | invalid |
|---|---|---|---|---|---|---|
| semipolar_11_22 | rectangular | 0 | 0 | 0 | 140 | 0 |
| semipolar_11_22 | rectangular | 0.5 | 0 | 0 | 140 | 0 |
| semipolar_11_22 | rectangular | 1 | 0 | 0 | 140 | 0 |
| semipolar_11_22 | deterministic_pair | 0 | 34 | 0 | 140 | 0 |
| semipolar_11_22 | deterministic_pair | 0.5 | 73 | 0 | 140 | 0 |
| semipolar_11_22 | deterministic_pair | 1 | 112 | 0 | 140 | 0 |
| m_plane | rectangular | 0 | 0 | 0 | 140 | 0 |
| m_plane | rectangular | 0.5 | 0 | 0 | 140 | 0 |
| m_plane | rectangular | 1 | 0 | 0 | 140 | 0 |
| m_plane | deterministic_pair | 0 | 45 | 0 | 140 | 0 |
| m_plane | deterministic_pair | 0.5 | 45 | 0 | 140 | 0 |
| m_plane | deterministic_pair | 1 | 45 | 0 | 140 | 0 |

## Maximum flux and eligible-g2 statistics

Maximum valid signal_flux_s across all core/shape/QW/Stark rows: 171666/s (row ST04753).
Core-grid paired-optical-pass rows (g2<0.5, flux>=1000/s, plus one_pair_valid for SET): 601 of 3360; g2_min=0.173554, g2_median=0.173554.

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

Round 1's headline grid was c-plane only, radius fixed at 10 nm (height in {1,2,3,4,5} nm, x_in in {0.15,0.25,0.4}, Q in {500,2000,10000}, current in {0.002,0.02,0.2} uA), no orientation axis, and no Stark bias/current traces. This piece adds a 4-orientation x 5-radius x 7-height x 4-T_hs x 3-screening x 2-regime core grid (3360 rows this run), lens/truncated_cone shape (864 rows) and QW-fluctuation (288 rows) geometry supplements, and bias/current Stark diagnostic traces (2040+1080 rows) with derivative-refinement convergence checks and screening-compatibility fits. Round 1's own results.md and artifacts under out/nitride_cavity/ are unchanged and are not reinterpreted here (out of scope for this piece).


## Shape supplement (mapping sensitivity, not a demonstrated shape accuracy)
864 rows (lens / truncated_cone, native and full-height variants). The numerical spread across shape variants at matched height/radius is a mapping-sensitivity diagnostic; systematic shape accuracy remains unquantified (no 3-D solver).

## QW-fluctuation supplement
288 rows; diode.wl_thickness_nm/x_in and dot.wl_thickness_nm/x_in are matched per card contract; d_i_nm=24+w, d_active_nm=height [A]. No lens-in-QW cases were run (unsupported combination).

## Literature comparisons (non-gating)

Wang et al., Sci. Rep. 7, 12089 (2017): uncapped a-plane AFM ~7 nm/~35 nm dots, 2.54 eV at 220 K, 19.0+/-0.4 meV linewidth, raw/corrected g2 0.47/0.21, optical excitation (76 MHz, 1 ps, 800 nm two-photon) -- NOT an electrical SPS or a 300 K validation. Model rows at the same nominal geometry (comparison-only, excluded from headline coverage):

- row LI04513: T_hs=220 K, E_X_eV=2.790254882472666, g2=0.7219685638585884, valid=True (comparison_only_reported_temperature)
- row LI04514: T_hs=300 K, E_X_eV=2.7644024701218513, g2=0.9822939349276038, valid=True (comparison_only_model_headline_temperature)
- Deshpande et al., APL 105, 141109 (2014) [V abstract-only; CONDITIONS INCOMPLETE], replay row LI04515: measured g2=0.29, model g2=0.9980609340138185.

Zhang et al., APL 108, 153102 (2016), Fig. 5: -10 meV/V below 2 V, non-gating (different device, excitation and unspecified temperature); shown only as a labelled guide anchored to a model row on the Stark figures, never a fitted target.

## Screening-compatibility table

The [-12,-8] meV/V window used when --slope-range is not supplied is EXPLICITLY illustrative around Zhang's approximate value, never a measured interval; the same --bias-window applies to any user-supplied slope. Nonpolar (m-plane/a-plane) rows are marked screening_unidentifiable: polarization_factor=0 makes screening_fraction physically inert for those orientations, so no compatibility test there can resolve a screening fraction -- this is a structural degeneracy, not a measurement result.

120 compatibility rows: 8 compatible, 112 incompatible, 0 screening_unidentifiable (nonpolar or degenerate).

## Derivative-refinement convergence checks

| check | anchor V_j | coarse meV/V | fine (half-step) meV/V | abs diff | rel diff | converged |
|---|---|---|---|---|---|---|
| 0 | 1 | -14.32 | -14.31 | 0.01005 | 0.0007013 | True |
| 1 | 1 | 2.925 | 2.925 | 0.0003173 | 0.0001085 | True |
| 2 | 1 | 1.913 | 1.914 | 0.0002305 | 0.0001205 | True |
| 3 | 1 | -10.6 | -10.6 | 0.007743 | 0.0007302 | True |
| 4 | 1 | -3.568 | -3.565 | 0.002579 | 0.000723 | True |
| 5 | 1.4 | 10.15 | 10.15 | 0.0007964 | 7.847e-05 | True |
| 6 | 1 | 0.8297 | 0.8298 | 0.000102 | 0.0001229 | True |
| 7 | 1 | -14.06 | -14.05 | 0.009253 | 0.0006581 | True |
| 8 | 1 | 1.913 | 1.914 | 0.0002305 | 0.0001205 | True |
| 9 | 1 | 2.925 | 2.925 | 0.0003173 | 0.0001085 | True |
| 10 | 1.6 | -30.28 | -30.25 | 0.038 | 0.001255 | True |
| 11 | 1 | 0.3747 | 0.3747 | 5.384e-05 | 0.0001437 | True |

## Limitations
Nonpolar strain, valence-band ordering and FSS are not modelled; shape mapping error is unquantified; finite-dot lateral fields, field-assisted escape and injection-dependent screening are not solved (screening_fraction is fixed along every current trace). m-plane and a-plane rows coincide under this scalar model (both polarization_factor=0) -- a model limitation, not independent evidence for either. SET hardware feasibility (island charging energy / RC bandwidth) is reported independently of the idealized optical pass and never substitutes for it. No field-free lifetime or screening fraction was ever fitted to force agreement with Wang or Zhang.

runtime_s=802.2 evaluate_calls=7706 complete=True
