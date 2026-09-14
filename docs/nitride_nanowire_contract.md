# Nitride nanowire contract

This is the piece-1 contract for the opt-in `ingan_gan_nanowire` tier
(`.workers/briefs/nitride-nanowire.md`, `.workers/specs/nitride-nanowire-*.md`,
nine pieces). It does not change `ingan_gan_planar`, legacy evaluation, or
the planar contracts (`docs/nitride_cavity_contract.md`,
`docs/nitride_geometry_stark_contract.md`), which remain unchanged. All
nanowire predictions are independent predictions: the Deshpande 2013/2014
results below are comparison anchors, never fit targets. This revision
(round 2 of piece 1) replaces the previous 97-line version after an Opus
FAIL review; it TRANSCRIBES every module, function signature, card leaf,
row key, VERDICT field, and output path from the nine specs' own
"Interface or signature constraints" and "Acceptance criteria" sections --
it never invents or renames what a spec already names. Where this document
and a spec appear to disagree, the spec wins; report the conflict rather
than silently resolving it in a future revision.

## Families and common card contract

`platform` is `ingan_gan_nanowire`; `emission.type` is `nanowire`;
`drive.diode.preset` is `nitride-nanowire`; and `cavity.enabled` is `false`
for BOTH the `horizontal_as_built` and `vertical_photonic` families. There
is no Q, DBR, cavity-Q axis, or planar aperture supply factor. Required
blocks are `nitride.nanowire`, `nitride.surface`, `nitride.photonics`,
`nitride.wire_thermal`, and `nitride.injector`; `nitride.dot` retains axial
height and composition (see "Card schema" below for every leaf). Headline
`nitride.tau_rad0_ns=1.0` [A, deliberately not the held-out 2014 lifetime]
and `drive.b_res=0.1` [A] are explicit defaults. `b_res=0.02` is only a
sensitivity, motivated by the 2013 corrected HBT value (0.16), never a
headline substitution.

`nitride.nanowire` leaves are `family` (`horizontal_as_built` or
`vertical_photonic`), `core_radius_nm`, `outer_radius_nm`, `strain_bound`
(`unrelaxed` or `relaxed`), `shell` (`none` or `AlGaN`), `barrier_left_nm`,
and `barrier_right_nm`. `nitride.dot.radius_nm` equals `core_radius_nm`;
`outer_radius_nm >= core_radius_nm`; shell `none` requires equality
(`outer_radius_nm == core_radius_nm`); shell `AlGaN` uses a declared 3 nm
shell thickness [A] (`outer_radius_nm = core_radius_nm + 3`). The disc
fills the semiconductor core. The independently explicit
`drive.diode.conducting_radius_nm` is in `(0, core_radius_nm]` and defaults
to the core; it is never an optical radius alias (piece 5's
`0 < conducting_radius_nm <= core_radius_nm`). A radius sweep updates
dot/core/outer/conducting radii, the SET island radius
(`drive.set_params.radius_nm`), and
`d_i_nm = barrier_left_nm + height_nm + barrier_right_nm` together.

The two strain endpoints are scenarios, not a universal relaxation curve:
`strain_bound=unrelaxed` sets `nitride.dot.strain_fraction=1` (piece 2's
`NitrideNanowireSystem.strain_bound` field, not a separate duplicate
input), `strain_bound=relaxed` sets `strain_fraction=0`. Both change
strain-deformation shifts and `P_pz`; `P_sp` remains (piece 2's own
materials-transfer rule: "Materials strain_fraction=0 for relaxed removes
both piezoelectric strain and strain-induced band shifts, while
spontaneous polarization discontinuity survives"). Screening is a separate
`screening_fraction=0` main-grid input, with `{0,1}` only on a reduced cut
(never an alias of the strain bound). Unrelaxed is labelled
conservative/lower and relaxed headline/upper (the design brief's user
answer 4), without asserting monotonic flux ordering between them -- a
reversal must be reported, not hidden. The radial model is a hard-wall
cylinder [A boundary envelope] with GaN axial reservoirs (piece 2's
`E_perp(m,n)=hbar^2*j_(m,n)^2/(2*m_xy*R^2)` applied on both the InGaN disc
and GaN axial sections); it must not model a full-core disc as laterally
surrounded by infinite GaN. Dielectric self energy, sidewall band bending,
elastic spatial variation, and lateral alloy localization remain
unresolved systematic error (piece 2's own docstring carries the
quantitative caveat).

## Family-specific interfaces

Horizontal-as-built means the dispersed, end-contacted wire lying flat on
100 nm thermal SiO2 on (001) Si (`deshpande2013_device_geometry` below),
normal-incidence objective collection, no mirror/taper/top contact, and a
subwavelength dielectric-antenna envelope (piece 3's horizontal family).
`NA=0.5` [A] and its collection scale are card leaves. The reported 70
percent axial degree of linear polarization (`deshpande2013_polarization`
below) constrains orientation/antenna weighting only; it is never an
automatic 70 percent collection efficiency.

Vertical-photonic means a designed HE11 wire (piece 3's vertical family):
`n_wire`, `n_ambient`, taper transmission, bottom reflectivity, top-contact
transmission, propagation transmission, and collection envelope are
explicit [A]/[E] card leaves. HE11 confinement, beta, first-lens
extraction, and radiative-rate factor remain separately normalized;
baseline radiative-rate factor is one [A] (a high beta never multiplies
gamma by that same fraction -- piece 3's own conjunction rule). Beta sums
both propagation directions; mirror return is incoherent (no uncontrolled
repeated-count mirror gain).

`nitride.surface` declares `S_cm_s=1e3` [E secondary attribution: Deshpande
2013 citing its own ref. 35, itself never read -- see "Provenance
corrections" below; this is NOT [V]], `shell`/`shell_multiplier` [A], and
`reservoir_access`/`occupied_dot_access` [A] (piece 4's
`NitrideNanowireSurfaceParams`). Surface loss is `k_side_ns =
2*S_cm_s/(core_radius_nm*1e-7)*1e-9` [DR] only in its declared channel
(reservoir vs. occupied-dot) and is never counted twice (piece 5 uses
reservoir loss in supply competition; the device piece adds occupied-dot
loss once to intrinsic `k_nr`).

`nitride.wire_thermal` declares series resistance (`R_s_ohm`; 2.38e9 ohm is
[V] for the 2013 device, 1e6 ohm is the [A] designed-contact alternative),
thermal resistance `Rth_K_W` [A], `f_Rs_local`/`eta_total` [A], and duty
averaging (piece 5's `wire_operating_point`/`pulse_delivery`).

`nitride.injector` declares barrier topology/material/thickness/alignment,
tunnelling and thermionic controls (piece 6's
`NitrideNanowireInjectorParams`); unsupported second-pair/hole-occupation
controls produce an honest failed or `unknown` RT screen
(`rti_feasible=False`, `rti_status`), never a fabricated pass.

## Module table

Every symbol below is transcribed VERBATIM (the "Signature" column is a
literal substring) from the named spec's own "Interface or signature
constraints" section. A rename, a dropped parameter, or a changed default
in either this table or the corresponding production module breaks that
substring match, and `verify/verify_nitride_nanowire_contract.py` fails.

| Module | Symbol | Signature | Verifier |
| --- | --- | --- | --- |
| `fsim_core/nitride_nanowire_levels.py` | `NitrideNanowireSystem` | `NitrideNanowireSystem carries height_nm, core_radius_nm, outer_radius_nm, x_in, strain_bound, screening_fraction, external_field_kVcm, vbo_InN_GaN_eV and strain_c_fraction` | `verify/verify_nitride_nanowire_levels.py` |
| `fsim_core/nitride_nanowire_levels.py` | `NanowireLevels` | `frozen NitrideNanowireSystem and NanowireLevels result types` | `verify/verify_nitride_nanowire_levels.py` |
| `fsim_core/nitride_nanowire_levels.py` | `levels` | `levels(system, T_K=300.0, *, z_points=1201, exterior_nm=45.0)` | `verify/verify_nitride_nanowire_levels.py` |
| `fsim_core/nitride_nanowire_levels.py` | `rates` | `rates(lv, T_K, *, tau_rad0_ns, tau_cap_ps, reservoir_length_nm, channel='min', k_nr_ns=0.0)` | `verify/verify_nitride_nanowire_levels.py` |
| `fsim_core/nitride_nanowire_photonics.py` | `NitrideNanowirePhotonicsParams` | `Create frozen NitrideNanowirePhotonicsParams and response(params, *, lambda_nm, outer_radius_nm, gamma_X0_ns, gamma_XX0_ns) -> dict` | `verify/verify_nitride_nanowire_photonics.py` |
| `fsim_core/nitride_nanowire_photonics.py` | `response` | `response(params, *, lambda_nm, outer_radius_nm, gamma_X0_ns, gamma_XX0_ns) -> dict` | `verify/verify_nitride_nanowire_photonics.py` |
| `fsim_core/nitride_nanowire_surface.py` | `NitrideNanowireSurfaceParams` | `Create frozen NitrideNanowireSurfaceParams and surface_rates(params, *, core_radius_nm, T_K) -> dict` | `verify/verify_nitride_nanowire_surface.py` |
| `fsim_core/nitride_nanowire_surface.py` | `surface_rates` | `surface_rates(params, *, core_radius_nm, T_K) -> dict` | `verify/verify_nitride_nanowire_surface.py` |
| `fsim_core/nitride_nanowire_surface.py` | `yield_ratio` | `yield_ratio(*, gamma_300_ns, loss_300_ns, gamma_10_ns, loss_10_ns)` | `verify/verify_nitride_nanowire_surface.py` |
| `fsim_core/nitride_nanowire_transport.py` | `NitrideWireDiode` | `Create NitrideWireDiode, wire_pin(**overrides)` | `verify/verify_nitride_nanowire_transport.py` |
| `fsim_core/nitride_nanowire_transport.py` | `wire_pin` | `wire_pin(**overrides)` | `verify/verify_nitride_nanowire_transport.py` |
| `fsim_core/nitride_nanowire_transport.py` | `evaluate_injection` | `evaluate_injection(diode, *, I_uA, T_K, tau_pulse_ns, E_X_eV, reservoir_energy_eV, barrier_e_eV, barrier_h_eV, surface_reservoir_ns, tau_cap_ps, S_dot, w_meV, eta_rad_matrix, eta_total) -> dict` | `verify/verify_nitride_nanowire_transport.py` |
| `fsim_core/nitride_nanowire_transport.py` | `wire_operating_point` | `wire_operating_point(diode, *, I_uA, T_hs_K, duty, Rth_K_W, eta_total, h_nu_eV) -> dict` | `verify/verify_nitride_nanowire_transport.py` |
| `fsim_core/nitride_nanowire_transport.py` | `pulse_delivery` | `pulse_delivery(diode, *, V_j, T_K, tau_pulse_ns, rep_rate_hz, C_parasitic_F=0.0) -> dict` | `verify/verify_nitride_nanowire_transport.py` |
| `fsim_core/nitride_nanowire_injector.py` | `NitrideNanowireInjectorParams` | `Create frozen NitrideNanowireInjectorParams, transmission(params, energy_eV, *, bias_V=0.0, field_kVcm=0.0) and injector_feasibility` | `verify/verify_nitride_nanowire_injector.py` |
| `fsim_core/nitride_nanowire_injector.py` | `transmission` | `transmission(params, energy_eV, *, bias_V=0.0, field_kVcm=0.0)` | `verify/verify_nitride_nanowire_injector.py` |
| `fsim_core/nitride_nanowire_injector.py` | `injector_feasibility` | `injector_feasibility(params, *, T_K, rep_rate_hz, loading_window_ns, electron_level_eV, hole_level_eV, electron_spacing_meV, hole_spacing_meV, second_pair_addition_meV, available_pair_rate_Hz, field_kVcm=0.0) -> dict` | `verify/verify_nitride_nanowire_injector.py` |
| `fsim_core/device.py` + `fsim_core/nitride_nanowire_device.py` | `platform` | `recognize platform='ingan_gan_nanowire'` | `verify/verify_nitride_nanowire_device.py` |
| `fsim_core/device.py` + `fsim_core/nitride_nanowire_device.py` | `evaluate_nanowire` | `dispatch to evaluate_nanowire(design, T_grid=None) in the new module` | `verify/verify_nitride_nanowire_device.py` |

### Module evidence map

Every module below has at least one ledger anchor (`verify/data/
nitride_nanowire_anchors.yaml`) carrying a non-null `value`, satisfying the
original piece-1 acceptance criterion "every required new module has at
least one published benchmark identified in the ledger." A module may also
list a `missing`/null anchor as a documented gap; that never substitutes
for the required non-null companion.

| Module | Ledger anchors |
| --- | --- |
| `nitride_nanowire_levels` | `deshpande2013_geometry` (non-null), `deshpande2013_emission_energy` (non-null) |
| `nitride_nanowire_photonics` | `maslov2004_he11` (null, documented gap), `claudon2010_extraction` (null, documented gap), `deshpande2013_polarization` (non-null), `deshpande2013_device_geometry` (non-null) |
| `nitride_nanowire_surface` | `deshpande2013_surface_velocity` (non-null), `deshpande2013_thermal_and_pl` (non-null) |
| `nitride_nanowire_transport` | `deshpande2013_geometry` (non-null), `deshpande2013_thermal_and_pl` (non-null) |
| `nitride_nanowire_injector` | `kitamura2026_architecture` (null, documented gap), `deshpande2013_coulomb_reference` (non-null) |
| `nitride_nanowire_device` | `deshpande2014_abstract` (non-null), `deshpande2013_hbt` (non-null), `deshpande2013_lifetimes` (non-null) |

## Card schema

Every leaf below carries a unit, a default, an allowed range, and a
provenance tag from `{V, DR, E, A}` (the card-field TRANSFER tag, distinct
from a ledger anchor's own measurement tag). "Both families" means the
leaf and its default/range apply identically to `horizontal_as_built` and
`vertical_photonic` unless a family-specific default is given.

### `nitride.nanowire` (both families)

| Leaf | Unit | Default | Range | Tag |
| --- | --- | --- | --- | --- |
| `family` | enum | required | `horizontal_as_built`, `vertical_photonic` | A |
| `core_radius_nm` (horizontal) | nm | 12.5 | `{10,12.5,15,20,25,40}` | V (2013 optical diameter/2) |
| `core_radius_nm` (vertical) | nm | 100.0 | `{60,80,100,120}` | A (design) |
| `outer_radius_nm` | nm | equal to `core_radius_nm` (shell none); `core_radius_nm+3` (shell AlGaN) | `>= core_radius_nm` | A |
| `strain_bound` | enum | required | `unrelaxed`, `relaxed` | DR (scenario label) |
| `shell` | enum | `none` | `none`, `AlGaN` | A |
| `barrier_left_nm` | nm | 15.0 | `>0`, sweep fixed at 15.0 | V (2013 geometry) |
| `barrier_right_nm` | nm | 15.0 | `>0`, sweep fixed at 15.0 | V (2013 geometry) |

### `nitride.dot` (nanowire-relevant leaves; full enum inherited from the planar accepted-field table, `docs/nitride_geometry_stark_contract.md`)

| Leaf | Unit | Default | Range | Tag |
| --- | --- | --- | --- | --- |
| `height_nm` | nm | 2.0 | sweep `{1.5,2,3,4}` | V (2013 disc thickness) |
| `radius_nm` | nm | equal to `nitride.nanowire.core_radius_nm` | derived, not independently set | A (equality constraint) |
| `x_in` | dimensionless | 0.40 | sweep `{0.25,0.40}` | V (0.25 2013; 0.40 2014) |
| `strain_fraction` | dimensionless | derived from `strain_bound` (1.0 unrelaxed, 0.0 relaxed) | `{0,1}` | DR |
| `screening_fraction` | dimensionless | 0.0 | main grid fixed 0; reduced cut `{0,1}` | A |
| `external_field_kVcm` | kV/cm | 0.0 | finite | A |
| `vbo_InN_GaN_eV` | eV | 1.15 | fixed | V (Rinke 2008) |
| `strain_c_fraction` | dimensionless | 0.7 | fixed | DR |

### `nitride.surface` (both families)

| Leaf | Unit | Default | Range | Tag |
| --- | --- | --- | --- | --- |
| `S_cm_s` | cm/s | 1.0e3 | sensitivity `{1e2,1e3,1e4}` | E (secondary attribution, ref. 35 unread) |
| `shell` | enum | `none` | `none`, `AlGaN` (mirrors `nitride.nanowire.shell`) | A |
| `shell_multiplier` | dimensionless | 1.0 (`none`), 0.1 (`AlGaN`) | `[0,1]` | A |
| `reservoir_access` | dimensionless | 1.0 | sensitivity `{0,0.1}` (occupied-dot only, per piece 8) | A |
| `occupied_dot_access` | dimensionless | 1.0 | sensitivity `{0,0.1}` | A (conservative full access) |

### `nitride.photonics` (both families; some leaves apply to one family only, see Notes)

| Leaf | Unit | Default | Range | Tag | Notes |
| --- | --- | --- | --- | --- | --- |
| `family` | enum | mirrors `nitride.nanowire.family` | `horizontal_as_built`, `vertical_photonic` | A | both |
| `NA` | dimensionless | 0.5 | `(0,1]` | A | horizontal (objective NA) |
| `n_wire` | dimensionless | GaN index at emission wavelength | provenance-tagged dispersion source | E | both |
| `n_ambient` | dimensionless | 1.0 (air) | fixed | A | both |
| `oxide_thickness_nm` | nm | 100.0 | fixed | V (2013 substrate) | horizontal |
| `collection_scale` | dimensionless | 1.0 | `[0,1]`, quantified assumption range | A | horizontal |
| `radiative_rate_factor` | dimensionless | 1.0 | independent envelope on sensitivity cuts | A | both |
| `taper_transmission` | dimensionless | 1.0 | `[0,1]` | A | vertical |
| `bottom_reflectivity` | dimensionless | 0.0 (no mirror baseline) | `[0,1]` | A | vertical |
| `top_contact_transmission` | dimensionless | 1.0 | `[0,1]` | A | vertical |
| `propagation_transmission` | dimensionless | 1.0 | `[0,1]` | A | vertical |

### `nitride.wire_thermal` (both families)

| Leaf | Unit | Default | Range | Tag |
| --- | --- | --- | --- | --- |
| `R_s_ohm` (horizontal) | ohm | 2.38e9 | sensitivity 1e6 [A designed contact] | V (2013 device) |
| `R_s_ohm` (vertical) | ohm | 1.0e6 | sensitivity 2.38e9 | A (designed contact) |
| `Rth_K_W` | K/W | 1.0e9 (horizontal), 1.0e7 (vertical) | one decade each way | A |
| `f_Rs_local` | dimensionless | 1.0 | fixed | A |
| `eta_total` | dimensionless | 0.01 | fixed baseline | A |
| `C_parasitic_F` | F | 0.0 | sensitivity `{1e-18,1e-16}` | A (optimistic default) |

### `nitride.injector` (both families; SET/`deterministic_pair` regime only)

| Leaf | Unit | Default | Range | Tag |
| --- | --- | --- | --- | --- |
| `topology` | enum | `double_barrier` | `single_barrier`, `double_barrier` | A |
| `al_fraction` | dimensionless | 0.30 | sensitivity `{0.2,0.3,0.4}` | A |
| `barrier_thickness_nm` | nm | 2.0 (each barrier) | sensitivity `{1,2,3}` | A |
| `well_width_nm` | nm | 4.0 (`double_barrier` only) | fixed | A |

### `drive.set_params` (SET/`deterministic_pair` regime only, both families)

| Leaf | Unit | Default | Range | Tag |
| --- | --- | --- | --- | --- |
| `radius_nm` | nm | equal to `nitride.nanowire.core_radius_nm` | derived, never independently resized by a geometry sweep | A |
| `eps_r` | dimensionless | 9.5 | fixed | E |
| `ec_margin` | dimensionless | 10.0 | fixed | A |
| `R_T_ohm` | ohm | 1.0e6 | fixed (unchanged SET convention) | A |

## Row columns

Every scalar key below is produced by `evaluate()` through
`fsim_core/nitride_nanowire_device.py`'s `evaluate_nanowire`. Names already
established by the planar evaluator (`fsim_core/device.py` `_evaluate_nitride`,
line 1333 and 1357 in the current tree) are reused unchanged where the
same quantity applies; nanowire-only names are transcribed verbatim from
pieces 2-6's own "Required ... outputs" bullets.

Identity / regime: `platform`, `family`, `T_hs`, `T_j`, `cycle_loading`,
`rep_rate_hz`, `strain_bound`, `bound_role`, `screening_fraction`,
`core_radius_nm`, `outer_radius_nm`, `conducting_radius_nm`.

Levels/optical (piece 2): `E_X_eV`, `lambda_nm`, `field_kVcm`,
`overlap_sq`, `electron_bound`, `hole_bound`, `E_a_meV`, `gamma_X0_ns`,
`gamma_XX0_ns`, `k_X_ns`, `k_XX_ns`, `escape_prefactor_ns`,
`tau_cap_ps_used`, `reservoir_state_count_e`, `reservoir_state_count_h`.

Photonics (piece 3): `radius_over_lambda`, `V_number`, `beta_HE11`,
`eta_collection_X`, `eta_collection_XX`, `radiative_rate_factor`,
`gamma_X_ns`, `gamma_XX_ns`, `approximation_error`.

Surface (piece 4): `k_side_ns`, `k_surface_reservoir_ns`,
`k_surface_X_ns`, `k_surface_XX_ns`, `shell_multiplier_used`,
`reservoir_access_used`, `occupied_dot_access_used`.

Transport (piece 5): `area_cm2`, `J_A_cm2`, `V_j`, `V_terminal`,
`depletion_field_kVcm`, `C_dep_F`, `eta_inj`, `f_capture`, `f_qfl_dot`,
`f_qfl_background`, `r_supply_s`, `r_captured_s`, `r_matrix_radiative_s`,
`r_matrix_nonradiative_s`, `r_surface_reservoir_s`, `r_leakage_s`, `mu`,
`power_on_W`, `accounting_residual_s`, `tau_RC_ns`,
`delivered_step_fraction`, `pulse_delivery_feasible`.

Lifetimes (device, piece 7): `tau_rad_bare_ns`, `tau_rad_photonic_ns`,
`tau_total_X_ns`.

Counting: `g2_op`, `collected_flux_pulsed_s`, `counting_converged`,
`one_pair_valid`, `pair_supply_possible`, `valid`, `invalid_reasons`,
`provenance`.

Coulomb screen (unchanged `drive_mech.set_feasibility` convention):
`set_feasible`, `set_E_C_meV`, `set_EC_over_kT`, `set_radius_nm`,
`set_radius_max_nm`, `set_C_sigma_F`, `set_R_T_over_RQ`, `set_f_max_Hz`.

RT injector screen (piece 6, transcribed verbatim): `rti_feasible`,
`rti_status`, `rti_transport_feasible`, `rti_level_margin_kT`,
`rti_orbital_margin_kT`, `rti_alignment_error_meV`, `rti_linewidth_meV`,
`rti_rate_Hz`, `rti_e_rate_Hz`, `rti_h_rate_Hz`, `rti_bypass_fraction`,
`rti_missed_load_probability`, `rti_second_pair_probability`,
`rti_growth_feasible`, `rti_failed_checks`, `rti_evidence_status`.

Gates (piece 7): `optical_pass`, `hardware_qualified`, `rti_qualified`,
`device_pass`, `rti_device_pass`.

`optical_pass = valid AND counting_converged AND g2_op < 0.5 AND
collected_flux_pulsed_s >= 1000 AND (one_pair_valid AND
pair_supply_possible for deterministic-pair rows)`. `hardware_qualified =
optical_pass AND cycle_loading==deterministic_pair AND set_feasible`.
`rti_qualified = optical_pass AND cycle_loading==deterministic_pair AND
rti_feasible`. Pulse-regime rows report the hardware fields
`not_applicable`/`False`. `device_pass` uses the Coulomb screen (mirrors
the planar convention); `rti_device_pass` is a separate alias of
`rti_qualified`. Neither hardware failure overwrites the idealized optical
statistics (`g2_op`, `collected_flux_pulsed_s`) computed upstream of it.

## Sweep grid, VERDICT format, and output paths

### Main grid (piece 9)

| Family | Axis | Values | Count |
| --- | --- | --- | --- |
| horizontal | core_radius_nm | 10, 12.5, 15, 20, 25, 40 | 6 |
| horizontal | height_nm | 1.5, 2, 3, 4 | 4 |
| horizontal | x_in | 0.25, 0.40 | 2 |
| horizontal | T_hs_K | 230, 250, 273, 300 | 4 |
| horizontal | regime | rectangular, deterministic_pair | 2 |
| horizontal | strain_bound | unrelaxed, relaxed | 2 |
| horizontal | rep_rate_hz | 80e6, 200e6 | 2 |
| vertical | core_radius_nm | 60, 80, 100, 120 | 4 |
| vertical | height_nm | 1.5, 2, 3, 4 | 4 |
| vertical | x_in | 0.25, 0.40 | 2 |
| vertical | T_hs_K | 230, 250, 273, 300 | 4 |
| vertical | regime | rectangular, deterministic_pair | 2 |
| vertical | strain_bound | unrelaxed, relaxed | 2 |
| vertical | rep_rate_hz | 80e6, 200e6 | 2 |

Horizontal main grid total: 6 * 4 * 2 * 4 * 2 * 2 * 2 = 1536 rows. Vertical
main grid total: 4 * 4 * 2 * 4 * 2 * 2 * 2 = 1024 rows. Core grid
(horizontal + vertical): 2560 rows. Reduced cuts (screening, `b_res`, `S`,
shell, reservoir/occupied access, NA/collection, linewidth/`tau_rad0`/
capture, `R_s`/`Rth` decades, parasitic capacitance, injector barrier
thickness/Al fraction/alignment/growth tolerance, current/pulse
sensitivity) are one-at-a-time, not Cartesian, and are counted separately
from the 2560-row core grid; the full run's total `evaluate()` call count
(core + reduced cuts + experiment-comparison + planar-reference rows) MUST
stay `<= 10000` and finish under 1800 s.

### VERDICT line format (piece 9, exact field list)

```
VERDICT: idealized_status=<idealized_status> family=<family> regime=<regime> strain_bound=<strain_bound> bound_role=<bound_role> rep_rate_hz=<rep_rate_hz> complete=<complete> eligible=<eligible> paired_optical_pass=<paired_optical_pass> hardware_qualified=<hardware_qualified> rti_qualified=<rti_qualified> coverage=<n>/<total> invalid=<invalid> flux_floor=1000/s
```

Every field name in this template (`idealized_status`, `family`, `regime`,
`strain_bound`, `bound_role`, `rep_rate_hz`, `complete`, `eligible`,
`paired_optical_pass`, `hardware_qualified`, `rti_qualified`, `coverage`,
`invalid`, `flux_floor`) is the exact list piece 9 prescribes ("VERDICT
lines retain established fields and add family, regime, strain_bound,
bound_role, rep_rate_hz, complete, eligible, paired_optical_pass,
hardware_qualified, rti_qualified, coverage, invalid, idealized_status and
flux_floor=1000/s"). `paired_optical_pass` counts an optical pass whose
matched strain-bound partner row is also present, not an AND requiring
both strain scenarios to individually pass (piece 9's own definition,
printed alongside the VERDICT lines). One VERDICT line is emitted per
`family` x `regime` x `strain_bound` combination in the main grid, plus
supplementary lines for planar-reference and sensitivity rows, mirroring
`scripts/run_nitride_geometry_stark.py`'s existing per-combination pattern.

### Output paths (piece 9)

`scripts/run_nitride_nanowire.py --out-dir <dir>` writes, under `<dir>`
(which must resolve inside `out/nitride_nanowire/`):

- `<dir>/sweep.csv` -- every evaluated row's coordinates and scalars.
- `<dir>/manifest.json` -- input hashes, source/card hashes, row kinds,
  grid specification, actual call count, runtime, completeness, and the
  plot-to-row/column/transform mapping.
- `<dir>/results.md` -- VERDICT lines, per-temperature tables, and the
  nominated headline row per family/regime/bound/rate.
- `<dir>/*.png` -- figures, flat under `<dir>` (matching
  `out/nitride_geometry_stark/`'s existing layout).

`--out-dir out/nitride_nanowire/quick` and `--out-dir
out/nitride_nanowire/full` are the two invocations the project test
command uses; `--dry-run` writes nothing and only prints the enumerated
counts.

## Evidence and verdict semantics

The separate ledger is `verify/data/nitride_nanowire_anchors.yaml`. Its
required fields are `citation`, `doi_or_url`, `location`, `evidence_status`,
`evidence_kind`, `tag`, `platform`, `excitation`, `observable`, `value`,
`unit`, `tolerance`, `transfer_notes`, `temperature`, `geometry`,
`raw_corrected`, and `observable_definition`. Status is `full_text`,
`abstract_only`, `figure_reading`, or `missing`; unknown values are null,
never zero. Every anchor with `value: null` carries `evidence_status:
missing`. No anchor tagged `V` cites a source marked unread or secondary
(that tag is reserved for a number read directly from the cited paper's
own full text, abstract, or figure).

2013 is CW electrical: X g2 raw/corrected `0.30/0.16`, XX `0.38/0.25`, at
1 nA. Its TRPL XX decay is 711 ps while HBT-fit X/XX correlation times are
1.1/0.7 ns; these are not interchangeable. Optical diameter is 25 +/- 5 nm;
the paper's thermal/current-density simulation cross section is 30 nm.
Thus 1 nA is 203.7 A/cm2 at 25 nm and 141.5 A/cm2 at 30 nm [DR], versus its
rounded 142 A/cm2. The reported 15/49 K rises are thermal-simulation
estimates, not thermometry or a universal thermal resistance. The 52
percent 10-to-300 K ensemble PL ratio is a loss comparison, not an S or
lifetime fit. The 2014 record is abstract-only: x=0.40, about 630 nm,
300 K, g2=0.29, 1.3 +/- 0.3 ns, up to 200 MHz. Geometry, shell, contact
layout and waveform are inherited assumptions, not verified 2014 details.
Supplement retrieval was not established and is marked missing.

### Provenance corrections (this revision)

- `S_cm_s=1e3` is `[E]`, not `[V]`: Deshpande 2013 (p.1) attributes the
  value to its own ref. 35 (GaN nanowire surface-recombination-velocity
  literature), which this project has never read. The number is a
  secondary attribution inside a primary source, not a number Deshpande
  2013 itself measured. See `deshpande2013_surface_velocity` below.
- The 2013 device's substrate stack (100 nm thermal SiO2 on (001) Si),
  Ti/Au 5/45 nm end contacts, and ~600 nm contacted length are `[V]`
  (Deshpande 2013, p.3 device description and p.6 methods). See
  `deshpande2013_device_geometry` below.
- The 70 percent axial degree of linear polarization and its
  dielectric-contrast attribution (thin wire in air acting as an antenna,
  not dot anisotropy) are `[V]` (Deshpande 2013, p.5). See
  `deshpande2013_polarization` below.
- The single-dot emission energy (X = 2.84 eV = 436.56 nm, XX 10 meV ABOVE
  X, antibinding, 10 K ambient / ~25 K junction from Joule heating at 1 nA)
  is `[V]` (Deshpande 2013, Fig. 3c). See `deshpande2013_emission_energy`
  below, which cross-references the existing
  `deshpande2013-xx-splitting` anchor in
  `verify/data/nitride_cavity_anchors.yaml` (same measurement, already
  reused unmodified by that ledger as `dot.delta_xx=-10.0 meV`).
- The Coulomb-screen reference value (design brief user answer 3) is `[DR]`:
  at r=12.5 nm, eps_r=9.5, `fsim_core.drive_mech.set_feasibility`'s own
  isolated-SPHERE convention (`C_sigma = 4*pi*eps0*eps_r*R`, matching
  `fsim_core/loading.py island_radius_nm`) gives E_C = 12.126 meV,
  E_C/kT = 0.612 at 230 K and 0.469 at 300 K -- the numbers the design
  brief and the device spec (piece 7 acceptance criterion 4) both state.
  The alternative isolated-DISC convention (`C = 8*eps0*eps_r*R`) gives a
  DIFFERENT E_C = 19.05 meV, E_C/kT = 0.961 at 230 K and 0.737 at 300 K.
  Both are recorded in `deshpande2013_coulomb_reference` below with their
  formulas; the sphere convention is the one `set_feasibility` (and hence
  piece 7's device evaluator) actually uses.
- `deshpande2013_supplement` no longer carries the main article's DOI: its
  `doi_or_url` is `null` with `evidence_status: missing`, so the citation
  gate reports a recorded gap (an explicit unretrieved-evidence note), not
  a DOI-resolves-fine PASS that misleadingly implies the supplement itself
  was checked.
- `kitamura2026_architecture`'s `value` is `null`; consistent with the
  null-value rule its `evidence_status` is `missing` (the abstract was
  read, but it names no numeric benchmark for this project's use, so no
  quantity is transcribed) and its tag is `A` (architecture-only
  reference, not a measured number), matching the same pattern already
  used by `maslov2004_he11` and `claudon2010_extraction`.

## Evidence ledger schema

`verify/data/nitride_nanowire_anchors.yaml` is a SEPARATE ledger from
`verify/data/nitride_cavity_anchors.yaml` and `verify/data/
nitride_geometry_stark_anchors.yaml` (both unmodified by this piece).
`doi_or_url` uses only identifiers actually printed in `../_goal/
nitride_digests.md` or independently resolvable via Crossref (never a
fabricated DOI). New anchors added by this revision:

- `deshpande2013_surface_velocity` -- `S_cm_s=1000.0` [E], secondary
  attribution to ref. 35 (unread), Deshpande 2013 p.1.
- `deshpande2013_device_geometry` -- substrate stack, contact metal/
  thickness, contacted length, collection objective [V], Deshpande 2013
  p.3, p.6.
- `deshpande2013_polarization` -- 70 percent axial DOLP and its
  dielectric-contrast attribution [V], Deshpande 2013 p.5.
- `deshpande2013_emission_energy` -- X/XX energies, antibinding ordering,
  ambient/junction temperatures [V], Deshpande 2013 Fig. 3c; cross-
  references `deshpande2013-xx-splitting` in the cavity ledger.
- `deshpande2013_coulomb_reference` -- E_C/E_C/kT at r=12.5 nm under both
  the sphere and disc conventions [DR], derived from the 2013 geometry and
  a stated formula, not a number Deshpande 2013 itself reports.

`deshpande2013_thermal_and_pl`'s `value` no longer carries
`surface_S_cm_s_secondary` (moved to `deshpande2013_surface_velocity` with
its own correct `[E]` tag, so the `[V]`-tagged thermal/PL anchor no longer
mixes a secondary-attributed number into a directly-measured one).

## Round trip and verifier scope

`verify/verify_nitride_nanowire_contract.py` parses this document (markdown
tables and the fenced VERDICT block; a small stdlib+yaml parser, no new
dependency) and the ledger, and asserts: every module/function signature
cell in the "Module table" appears verbatim in the corresponding spec
file's "Interface or signature constraints" section; every card leaf has a
unit, a default, and a tag in `{V, DR, E, A}`; the sweep grid arithmetic in
this document equals the product of its own axis counts and equals the
independently expected 1536/1024; the VERDICT template contains every
field piece 9's own field list requires; every ledger anchor tagged `V`
has a non-null `doi_or_url` or a `title`; no `V`-tagged anchor's
`transfer_notes` reads as unread/secondary; every `value: null` anchor has
`evidence_status: missing`; no anchor has `value: 0`/`0.0` with
`evidence_status: missing`; every module in the "Module evidence map" has
at least one anchor with a non-null value; `deshpande2013_supplement` has
a null `doi_or_url`; and `cavity.enabled` is documented `false` for both
families. It prints `N/N nitride nanowire contract checks passed` and
exits 0 iff all pass. A corrupted copy of this document (one module-table
function name changed) is exercised once, in memory, against the same
parser to confirm the substring check fails; the on-disk document is never
modified by that self-test.
