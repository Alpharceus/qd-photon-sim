# Nitride nanowire contract

This is the piece-1 contract for the opt-in `ingan_gan_nanowire` tier
(`.workers/briefs/nitride-nanowire.md`, `.workers/specs/nitride-nanowire-*.md`,
nine pieces). It does not change `ingan_gan_planar`, legacy evaluation, or
the planar contracts (`docs/nitride_cavity_contract.md`,
`docs/nitride_geometry_stark_contract.md`), which remain unchanged. All
nanowire predictions are independent predictions: the Deshpande 2013/2014
results below are comparison anchors, never fit targets. This revision
(round 3 of piece 1, "contract-fix2") replaces round 2 after a second Opus
FAIL review (1 high, 5 medium, 7 low -- see `.workers/review/
nitride-nanowire-contract-opus-findings.md`) and folds in the Opus physics
coherence audit's H6 vertical-family correction and DEVICE-PIECE
CONSTRAINTS 1-10 (`.workers/review/
nitride-nanowire-physics-opus-coherence-findings.md`). It TRANSCRIBES every
module, function signature, card leaf, row key, VERDICT field, and output
path from the nine specs' own "Interface or signature constraints" and
"Acceptance criteria" sections -- it never invents or renames what a spec
already names. Where this document and a spec appear to disagree, the spec
wins; report the conflict rather than silently resolving it in a future
revision. Where a committed module's fix/directive round has drifted past
its own frozen spec text (the injector's `carrier` keyword, its new
`reflection` function, and its `gate_ns` keyword; transport's
GaN-reservoir `evaluate_injection` keywords), the Module table below
records BOTH the frozen spec-verbatim row (still checked against the
spec) and a separate, clearly labelled "module-only" row bound solely to
the live module (see "Module table" below); neither replaces the other.

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
repeated-count mirror gain). Above the LP11 cutoff (`V_number > 2.404826`)
`beta_HE11` carries an explicit multimode penalty (`beta_multimode_penalty`,
Bleuse 2011 / Claudon 2010 [E]-calibrated relative decline) and
`single_mode=False`; a vertical headline row is never nominated above
cutoff or with a nonzero `approximation_error`.

Coherence-review correction H6 (design correction, binding for this and
every later piece): `vertical_photonic` is a designed dot INSIDE a large
photonic wire, not a full-core disc. `nitride.dot.radius_nm` (the disc
radius, `NitrideNanowireSystem.disc_radius_nm` in
`fsim_core/nitride_nanowire_levels.py`) is STRICTLY LESS than
`nitride.nanowire.core_radius_nm` (the wire/reservoir radius) for this
family -- default disc radius 12.5 nm inside a 60-120 nm core, laterally
confined by the InGaN/GaN finite radial barrier (`nitride_levels.
finite_disk_2d`-style matching), never by a full-core hard wall. The
`horizontal_as_built` family keeps the pre-H6 convention, disc radius
equal to core radius (the as-measured 2013 device has no separate
dot-in-wire structure). Surface loss uses the core radius (the dot sits
far from the sidewall; access weight captures the evanescent overlap via
`sidewall_overlap`, a `NanowireLevels` field); photonics uses the outer
radius. This replaces the pre-H6 vertical grid (which forced a 120-240 nm
diameter InGaN quantum WELL, `E_perp` 0.11-0.44 meV vs kT, to stand in for
a two-level emitter) and is the reason `nitride.dot.radius_nm`'s vertical
default/range differs from `core_radius_nm` in the Card schema below.

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

## Composition rules for the device piece

Binding on piece 7 (device) and piece 9 (sweep); transcribed from the
physics coherence review's DEVICE-PIECE CONSTRAINTS 1-10.

1. **Radius plumbing.** `core_radius_nm` feeds `NitrideNanowireSystem`
   (the disc AND the axial reservoir for `horizontal_as_built`; for
   `vertical_photonic` the disc radius is the separate, smaller
   `nitride.dot.radius_nm` -- see H6 above), `surface_rates(core_radius_nm=
   )`, and `set_feasibility(radius_nm=)` (the disc radius). `drive.diode.
   conducting_radius_nm` (`<= core_radius_nm`) feeds `NitrideWireDiode`
   only. `nitride.nanowire.outer_radius_nm` feeds `photonics.response`
   only. The device piece must assert `nitride.dot.radius_nm ==
   drive.set_params.radius_nm` and reject an explicit `set_params.
   radius_nm` / `C_sigma_F` override that disagrees. Changing only
   `outer_radius_nm` must move `V_number`/`eta_collection` and nothing
   else.
2. **Units.** Rates crossing module boundaries are `ns^-1` EXCEPT
   transport's `r_*_s` outputs (`s^-1`) and the injector's `rti_*_Hz`
   outputs (`Hz`, i.e. `s^-1`). `available_pair_rate_Hz` (injector input)
   equals transport's `r_captured_s`. `surface_reservoir_ns` (transport's
   `evaluate_injection` keyword, aliased `k_surface_reservoir_per_ns`) is a
   RATE in `ns^-1` despite the `_ns` suffix -- never invert it as a
   lifetime. `tau_cap_ps` / `tau_matrix_ns` are TIMES, not rates. Always
   pass `evaluate_injection`'s surface-reservoir argument by keyword
   (`surface_reservoir_ns=` or `k_surface_reservoir_per_ns=`); never rely
   on `NitrideWireDiode.reservoir_surface_ns` (dead/absent on the
   dataclass).
3. **Temperature flow.** Solve `wire_operating_point(T_hs_K=...)` first,
   take its `T_j_K`, and evaluate `levels`, `rates`, `surface_rates`,
   `photonics.response` (`lambda_nm` from `levels` at `T_j`),
   `evaluate_injection`, `set_feasibility`, and `injector_feasibility` all
   at that `T_j`. Report `T_hs`, `T_j`, `thermal_iterations`, and
   `thermal_converged` (the convergence reason on failure) as row columns.
   A row whose `T_j` sits below the diode kernel's numerical floor is
   TRANSPORT-INVALID (an explicit invalid row, never NaN-propagated
   silently); the 2013 10 K replay is the one row where the kernel is
   valid (log-space `vbi`/`vj_of_j`, `_GaNJunctionKernel`).
4. **Surface, no 4x biexciton double count.** Either (a) call
   `rates(k_nr_ns=k_intrinsic + k_surface_X_ns)` and use the returned
   `k_X_ns`/`k_XX_ns` UNCHANGED (`k_XX_ns` already carries `2 *
   k_surface_X_ns` as `k_surface_XX_ns`), or (b) call `rates(k_nr_ns=
   k_intrinsic)` and then add `k_surface_X_ns` to `k_X_ns` and
   `k_surface_XX_ns` to `k_XX_ns` by hand. Never both in the same row.
   `k_surface_reservoir_ns` goes ONLY into `evaluate_injection(
   surface_reservoir_ns=)`, never into post-capture survival.
5. **Single strain-bound switch.** `nitride.nanowire.strain_bound` is the
   single source of truth, feeding `NitrideNanowireSystem.strain_bound`
   directly. `nitride.dot.strain_fraction` is DISPLAY-ONLY: validated
   (`1.0` for `unrelaxed`, `0.0` for `relaxed`), rejected on contradiction,
   and never passed to `levels` (it is not a `NitrideNanowireSystem`
   constructor field). Emit `bound_role` as a label (`conservative_lower`
   for unrelaxed, `headline_upper` for relaxed) and never assert relaxed is
   the flux upper bound; emit `bound_reversal=True` when `mu`,
   `collected_flux_pulsed_s`, or `g2_op` orders the strain-bound pair
   against that label. A row whose strain-bound partner is invalid is
   carried as unpaired, never dropped.
6. **Two non-gating hardware screens plus the RC diagnostic.** Report
   `set_feasible`/`set_EC_over_kT` (Coulomb screen), the injector's
   `rti_*` columns, and `tau_RC_ns` / `delivered_step_fraction` /
   `pulse_delivery_feasible` (from `pulse_delivery`) as independent
   columns. `optical_pass` reads none of them. `hardware_qualified` and
   `rti_qualified` are separate conjunctions with `optical_pass`;
   `rti_device_pass` is a distinct alias of `rti_qualified`. Never set
   `eta_load` in `deterministic_cycle_g2` from `rti_missed_load_
   probability`; never fold `delivered_step_fraction` into
   `collected_flux_pulsed_s`. Instead emit `collected_flux_delivered_s =
   collected_flux_pulsed_s * delivered_step_fraction` as its own column.
   Results text must state that the as-built 2.38 Gohm device cannot
   deliver a 100 ps step (`tau_RC_ns` 0.82-1.19, `delivered_step_fraction`
   0.08-0.12) and only the `R_s_ohm=1e6` designed-contact card can (the RC
   caveat).
7. **Injector energy zero and inputs.** `electron_level_eV =
   -dE_e_meV/1000`, `hole_level_eV = -dE_h_meV/1000` (both carriers'
   levels referenced to THEIR OWN bulk GaN band edge, per
   `fsim_core/nitride_nanowire_injector.py`'s own module docstring --
   `NanowireLevels.dE_e_meV`/`dE_h_meV` are escape depths below that edge,
   so the sign flip is required). `electron_spacing_meV`/
   `hole_spacing_meV = sp_split_e_meV`/`sp_split_h_meV` (NaN-guarded).
   `available_pair_rate_Hz = r_captured_s`. `field_kVcm =
   depletion_field_kVcm`. `carrier` must be passed EXPLICITLY for holes
   (`transmission`/`reflection` default to `carrier='electron'`; an
   omitted keyword silently prices holes with electron parameters).
   `loading_window_ns` is strictly `< 1/rep_rate_hz`.
   `second_pair_addition_meV = set_E_C_meV`, labelled the Coulomb charging
   energy, NEVER the optical XX splitting. `rti_feasible=False` must not
   be reported as a demonstrated physics result until the physics
   coherence review's H1/H4/H5 findings are resolved in the injector
   module.
8. **Surface access defaults.** `reservoir_access=1.0` [A]. `occupied_dot_
   access` headline default is `0.05` [DR] (see the Card schema below) with
   `1.0` [A] recorded as the declared conservative partner side by side,
   and `0.1` kept only as a sensitivity; never chosen to reproduce the
   held-out 0.71/1.1/1.3 ns anchors. Results text states the access-1.0
   lifetime cap (0.625 ns X, 0.3125 ns XX at `core_radius_nm=12.5`) next to
   whichever default is used.
9. **Photonics weights.** `dipole_weights` default is isotropic
   `(1/3,1/3,1/3)` [A, orientation prior unknown]. An earlier directive
   round set the default to `(0.0, 0.5, 0.5)` [DR] on the reasoning that a
   c-plane disc's transition dipole lies IN the c-plane, perpendicular to
   `along_wire`; the following photonics fix round found that this
   c-plane-only prior predicts a NEGATIVE `degree_of_linear_polarization`
   at this geometry -- the WRONG SIGN against the `deshpande2013_
   polarization` 70 percent anchor -- and reverted the default to
   isotropic, keeping `(0.0, 0.5, 0.5)` (`CPLANE_ONLY_DIPOLE_WEIGHTS`) as a
   named, explicitly falsified sensitivity rather than the headline.
   `eta_collection` is applied once per X/XX/background channel; `gamma`
   is never multiplied by `beta` or by `eta_collection` (they are
   independent envelopes); the antenna screening factor is applied via its
   own `antenna_rate_factor` output, never folded into `gamma` directly.
10. **Results text obligations.** Every results table/figure states: the
    RC caveat (bullet 6), the access-1.0 lifetime cap (bullet 8), that the
    two strain-bound anchors (2013 relaxed-matched, 2014 unrelaxed-matched)
    are matched by OPPOSITE endpoints and are never averaged (one of
    `x_in` transfer, disc thickness, VBO/bowing, or lateral localization
    carries roughly 0.3 eV of the discrepancy -- state it, do not resolve
    it), and an explicit E_C/kT wall statement: deterministic loading at
    230-300 K fails on `set_EC_over_kT`/E_C/kT by any charging mechanism
    priced in this tier at `core_radius_nm >= 10` nm.

## Module table

Every symbol below is transcribed VERBATIM (the "Signature" column is a
literal substring) from the named spec's own "Interface or signature
constraints" section. A rename, a dropped parameter, or a changed default
in either this table or the corresponding production module breaks that
substring match, and `verify/verify_nitride_nanowire_contract.py` fails.
This is now TRUE, not aspirational: for each of the five committed modules
(levels, photonics, surface, transport, injector) the verifier additionally
imports the live module and asserts, independently of the spec text, (a)
the Symbol cell resolves via `getattr` (a Symbol rename is caught even if
the Signature text is untouched -- fix-1 finding 12), and (b) for a
function-symbol row the module's own `inspect.signature` (annotations
stripped, defaults formatted, whitespace normalised) equals the row's
Signature cell; for a dataclass-symbol row every constructor field name
(`dataclasses.fields`) appears as a leaf in the corresponding Card schema
table below (see "Card schema"). Four rows below intentionally have NO
spec text at all (marked "module-only" in the Verifier cell): a
fix/directive round added a keyword or a whole function to the live
module after its spec was frozen, and the frozen row above each one is
kept verbatim so the original spec-substring check keeps passing
unchanged.

| Module | Symbol | Signature | Verifier |
| --- | --- | --- | --- |
| `fsim_core/nitride_nanowire_levels.py` | `NitrideNanowireSystem` | `NitrideNanowireSystem carries height_nm, core_radius_nm, outer_radius_nm, disc_radius_nm, x_in, strain_bound, screening_fraction, external_field_kVcm, vbo_InN_GaN_eV and strain_c_fraction` | `verify/verify_nitride_nanowire_levels.py` |
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
| `fsim_core/nitride_nanowire_transport.py` | `evaluate_injection` (GaN-reservoir directive round) | `evaluate_injection(diode, *, I_uA, T_K, tau_pulse_ns, E_X_eV, reservoir_energy_eV, barrier_e_eV, barrier_h_eV, surface_reservoir_ns=None, tau_cap_ps, S_dot, w_meV, eta_rad_matrix, eta_total, k_surface_reservoir_per_ns=None)` | `verify/verify_nitride_nanowire_transport.py` (module-only: introspected live; the transport spec predates the `surface_reservoir_ns=None` default and the `k_surface_reservoir_per_ns` alias keyword added by the piece-5 directive round, and also predates the `f_qfl_dot_thermodynamic_limit` output key -- see Row columns) |
| `fsim_core/nitride_nanowire_transport.py` | `wire_operating_point` | `wire_operating_point(diode, *, I_uA, T_hs_K, duty, Rth_K_W, eta_total, h_nu_eV) -> dict` | `verify/verify_nitride_nanowire_transport.py` |
| `fsim_core/nitride_nanowire_transport.py` | `pulse_delivery` | `pulse_delivery(diode, *, V_j, T_K, tau_pulse_ns, rep_rate_hz, C_parasitic_F=0.0) -> dict` | `verify/verify_nitride_nanowire_transport.py` |
| `fsim_core/nitride_nanowire_injector.py` | `NitrideNanowireInjectorParams` | `Create frozen NitrideNanowireInjectorParams, transmission(params, energy_eV, *, bias_V=0.0, field_kVcm=0.0) and injector_feasibility` | `verify/verify_nitride_nanowire_injector.py` |
| `fsim_core/nitride_nanowire_injector.py` | `transmission` | `transmission(params, energy_eV, *, bias_V=0.0, field_kVcm=0.0)` | `verify/verify_nitride_nanowire_injector.py` |
| `fsim_core/nitride_nanowire_injector.py` | `transmission` (carrier keyword, fix-1 directive round) | `transmission(params, energy_eV, *, bias_V=0.0, field_kVcm=0.0, carrier='electron')` | `verify/verify_nitride_nanowire_injector.py` (module-only: introspected live; the injector spec predates the `carrier` keyword -- a device-piece caller that omits it silently prices holes with electron parameters, fix-1 finding 6) |
| `fsim_core/nitride_nanowire_injector.py` | `reflection` (fix-1 directive round) | `reflection(params, energy_eV, *, bias_V=0.0, field_kVcm=0.0, carrier='electron')` | `verify/verify_nitride_nanowire_injector.py` (module-only: introspected live; public `reflection()` was added by the fix-1 round and has no spec text, fix-1 finding 6) |
| `fsim_core/nitride_nanowire_injector.py` | `injector_feasibility` | `injector_feasibility(params, *, T_K, rep_rate_hz, loading_window_ns, electron_level_eV, hole_level_eV, electron_spacing_meV, hole_spacing_meV, second_pair_addition_meV, available_pair_rate_Hz, field_kVcm=0.0) -> dict` | `verify/verify_nitride_nanowire_injector.py` |
| `fsim_core/nitride_nanowire_injector.py` | `injector_feasibility` (`gate_ns` keyword, directive round) | `injector_feasibility(params, *, T_K, rep_rate_hz, loading_window_ns, electron_level_eV, hole_level_eV, electron_spacing_meV, hole_spacing_meV, second_pair_addition_meV, available_pair_rate_Hz, field_kVcm=0.0, gate_ns=None) -> dict` | `verify/verify_nitride_nanowire_injector.py` (module-only: introspected live; commit `e696bdd` added an explicit `gate_ns` keyword separating the counting gate from `loading_window_ns`, postdating the frozen injector spec) |
| `fsim_core/device.py` + `fsim_core/nitride_nanowire_device.py` | `platform` | `recognize platform='ingan_gan_nanowire'` | `verify/verify_nitride_nanowire_device.py` |
| `fsim_core/device.py` + `fsim_core/nitride_nanowire_device.py` | `evaluate_nanowire` | `dispatch to evaluate_nanowire(design, T_grid=None) in the new module` | `verify/verify_nitride_nanowire_device.py` |

The device/sweep rows (piece 7/9) are not yet committed
(`fsim_core/nitride_nanowire_device.py` does not exist on disk at this
revision); the verifier binds them to their spec text only, exactly as
before, and skips the live-module checks for those two rows until the
module exists.

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
| `nitride_nanowire_photonics` | `maslov2004_he11` (null, documented gap), `claudon2010_extraction` (null, documented gap), `bleuse2011_claudon2010_beta_envelope` (non-null), `deshpande2013_polarization` (non-null), `deshpande2013_device_geometry` (non-null) |
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
| `radius_nm` (horizontal) | nm | equal to `nitride.nanowire.core_radius_nm` | derived, not independently set | A (equality constraint) |
| `radius_nm` (vertical) | nm | 12.5 (H6 disc-in-wire default, maps to `NitrideNanowireSystem.disc_radius_nm`) | `(0, core_radius_nm)`, sweep fixed 12.5 unless a cut says otherwise | A (H6 design correction: coherence review, see "Family-specific interfaces") |
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
| `reservoir_access` | dimensionless | 1.0 | sensitivity `{0,0.1}` | A (conservative full access) |
| `occupied_dot_access` | dimensionless | 0.05 (declared conservative partner 1.0 [A], see "Composition rules for the device piece" bullet 8) | sensitivity `{0,0.1,1.0}`; 0.1 kept only as a sensitivity; never tuned to the held-out 0.71/1.1/1.3 ns anchors | DR (M1: J0 hard-wall ground-state probability in a 2 nm sidewall capture layer over the uniform-density value, at `core_radius_nm=12.5`) |

### `nitride.photonics` (both families; some leaves apply to one family only, see Notes)

| Leaf | Unit | Default | Range | Tag | Notes |
| --- | --- | --- | --- | --- | --- |
| `family` | enum | mirrors `nitride.nanowire.family` | `horizontal_as_built`, `vertical_photonic` | A | both |
| `NA` | dimensionless | 0.5 | `(0,1]` | A | horizontal (objective NA) |
| `n_wire` | dimensionless | None (falls back to `gan_ordinary_index(lambda_nm)`) | provenance-tagged dispersion source | E | both |
| `n_ambient` | dimensionless | 1.0 (air) | fixed | A | both |
| `n_oxide` | dimensionless | None (falls back to `sio2_index(lambda_nm)`, Malitson 1965) | provenance-tagged override | V (default)/A (override) | horizontal |
| `oxide_thickness_nm` | nm | 100.0 | fixed | V (2013 substrate) | horizontal |
| `n_substrate` | dimensionless, complex | None (falls back to a 2-point Si n,k anchor table) | provenance-tagged override | E (default)/A (override) | horizontal |
| `dipole_weights` | tuple (along_wire, transverse_inplane, vertical), dimensionless | `(1/3, 1/3, 1/3)` (isotropic) | `CPLANE_ONLY_DIPOLE_WEIGHTS=(0.0,0.5,0.5)` [DR] kept as a named, explicitly falsified sensitivity (predicts the wrong-sign DOLP against the `deshpande2013_polarization` anchor) | A (orientation prior unknown; see "Composition rules for the device piece" bullet 9) | both |
| `emitter_height_nm` | nm | None (falls back to `outer_radius_nm`) | provenance-tagged override | A | horizontal |
| `collection_scale` | dimensionless | 1.0 | `[0,1]`, quantified assumption range | A | horizontal |
| `radiative_rate_factor` | dimensionless | 1.0 | independent envelope on sensitivity cuts | A | both |
| `beta_scale` | dimensionless | 1.0 | independent envelope on sensitivity cuts | A | vertical |
| `taper_transmission` | dimensionless | 1.0 | `[0,1]` | A | vertical |
| `bottom_reflectivity` | dimensionless | 0.0 (no mirror baseline) | `[0,1]` | A | vertical |
| `top_contact_transmission` | dimensionless | 1.0 | `[0,1]` | A | vertical |
| `propagation_transmission` | dimensionless | 1.0 | `[0,1]` | A | vertical |
| `unguided_collection_scale` | dimensionless | 0.0 | `[0,1]` | A | vertical |
| `taper_output_mfr_nm` | nm | None (falls back to the bare-wire Marcuse mode-field radius) | positive if set | A | vertical |

### `nitride.wire_thermal` (both families)

| Leaf | Unit | Default | Range | Tag |
| --- | --- | --- | --- | --- |
| `R_s_ohm` (horizontal) | ohm | 2.38e9 | sensitivity 1e6 [A designed contact] | V (2013 device) |
| `R_s_ohm` (vertical) | ohm | 1.0e6 | sensitivity 2.38e9 | A (designed contact) |
| `Rth_K_W` | K/W | 1.0e9 (horizontal), 1.0e7 (vertical) | one decade each way | A |
| `f_Rs_local` | dimensionless | 1.0 | fixed | A |
| `eta_total` | dimensionless | 0.01 | fixed baseline | A |
| `C_parasitic_F` | F | 0.0 | sensitivity `{1e-18,1e-16}` | A (optimistic default) |

### `drive.diode` (both families; `fsim_core.nitride_nanowire_transport.NitrideWireDiode` constructor fields not already carried by `nitride.nanowire`/`nitride.dot`/`nitride.wire_thermal`)

`NitrideWireDiode` has 15 constructor fields; six are already exact-name
leaves elsewhere in this schema (`core_radius_nm`, `barrier_left_nm`,
`barrier_right_nm` mirror `nitride.nanowire`; `x_in` mirrors `nitride.dot`;
`R_s_ohm`, `f_Rs_local` mirror `nitride.wire_thermal`) and are not
repeated here; this table carries the remaining nine so the verifier's
live field-set check has a named leaf for every constructor field.

| Leaf | Unit | Default | Range | Tag |
| --- | --- | --- | --- | --- |
| `N_A` | cm^-3 | 5.0e17 (p-GaN reservoir) | fixed | V (Deshpande et al. 2013; DIRECTIVE ROUND H2) |
| `N_D` | cm^-3 | 3.0e18 (n-GaN reservoir) | fixed | V (Deshpande et al. 2013; DIRECTIVE ROUND H2) |
| `n_ideality` | dimensionless | 2.0 | fixed | A |
| `tau_SRH_ns` | ns | 0.01 | fixed | DR (back-solved so the GaN-kernel dark current reproduces the inferred 2.6-3.1 V junction-voltage window at 1 nA/300 K; not an independent SRH-lifetime measurement) |
| `eps_r` | dimensionless | 10.28 | fixed | V/E (GaN dielectric constant) |
| `T` | K | 300.0 | overridden every call by the device piece's `T_j` from `wire_operating_point` (see "Composition rules for the device piece" bullet 3); not an independently swept card default | A |
| `conducting_radius_nm` | nm | equal to `core_radius_nm` | `(0, core_radius_nm]`, never an optical radius alias | A (DEVICE-PIECE CONSTRAINT 1) |
| `d_active_nm` | nm | equal to `nitride.dot.height_nm` | derived, not independently set | DR |
| `tau_matrix_ns` | ns | 1.0 | fixed | A |

### `nitride.injector` (both families; SET/`deterministic_pair` regime only)

Every leaf below is a constructor field of
`fsim_core.nitride_nanowire_injector.NitrideNanowireInjectorParams` as
committed after the piece-6 fix-1 round (electron/hole split, growth
tolerance, alignment gate) AND the subsequent directive round (commit
`e696bdd`: polarization sheet charges, emitter-quasi-Fermi-level
alignment, Mg-acceptor ionization, disc-E_C second-pair exclusion, a
separate counting gate). The verifier asserts this leaf set equals
`dataclasses.fields(NitrideNanowireInjectorParams)` exactly, live.

| Leaf | Unit | Default | Range | Tag |
| --- | --- | --- | --- | --- |
| `al_fraction` | dimensionless | 0.30 | sensitivity `{0.2,0.3,0.4}` | A |
| `electron_topology` | enum | `double_barrier` | `single_barrier`, `double_barrier` | A |
| `hole_topology` | enum | `double_barrier` | `single_barrier`, `double_barrier` | A |
| `electron_barrier_thickness_nm` | nm | 2.0 (each barrier) | sensitivity `{1,2,3}` | A |
| `hole_barrier_thickness_nm` | nm | 2.0 (each barrier) | sensitivity `{1,2,3}` | A |
| `electron_well_width_nm` | nm | 4.0 (`double_barrier` only) | fixed | A |
| `hole_well_width_nm` | nm | 4.0 (`double_barrier` only) | fixed | A |
| `growth_step_nm` | nm | 0.2593 (GaN c-axis bilayer spacing c/2) | fixed | V (Bernardini PRB 56, R10024 (1997), lattice constant) |
| `growth_tolerance_steps` | growth steps | 1.0 | sensitivity, injector barrier growth-tolerance cut | A |
| `me_barrier_override` | m0 | None (VCA-interpolated barrier electron mass used) | optional override, positive if set | A |
| `mh_barrier_override` | m0 | None (VCA-interpolated barrier hole mass used) | optional override, positive if set | A |
| `dEc_eV_override` | eV | None (VCA-interpolated conduction offset used) | optional override | A |
| `dEv_eV_override` | eV | None (VCA-interpolated valence offset used) | optional override | A |
| `me_well` | m0 | 0.209 (GaN electron effective mass) | fixed | V (nitride_materials) |
| `mh_well` | m0 | 1.88 (GaN hole effective mass) | fixed | V (nitride_materials) |
| `delta_Ev_GaN_AlN_eV` | eV | 0.70 | sensitivity `{0.30 Tsai and Bayram partition}` | V (Martin, Yu, Waldrop, APL 68, 2541 (1996), 0.70 +/- 0.24 eV) |
| `n_cm3` | cm^-3 | 3.0e18 (n-GaN reservoir) | fixed | V (Deshpande et al. 2013) |
| `p_cm3` | cm^-3 | 5.0e17 (p-GaN reservoir) | fixed | V (Deshpande et al. 2013) |
| `alignment_uncertainty_meV` | meV | 15.0 | fixed | A |
| `degeneracy` | dimensionless | 2.0 (spin/valley) | fixed | A |
| `reservoir_state_count_e` | dimensionless | 2.0 | fixed | A (directive round, `e696bdd`) |
| `reservoir_state_count_h` | dimensionless | 2.0 | fixed | A (directive round, `e696bdd`) |
| `mg_acceptor_energy_meV` | meV | 170.0 | fixed | V (directive round, `e696bdd`; Mg acceptor ionization energy, p-GaN free-hole fraction) |
| `include_polarization` | bool | `True` | `{True, False}` | DR (directive round, `e696bdd`; pseudomorphic fixed-D polarization sheet-charge envelope at every AlGaN/GaN interface) |
| `bypass_prefactor` | dimensionless | 1.0 | fixed | A |
| `field_leverarm` | dimensionless | 1.0 | fixed | A |
| `alignment_tunable` | bool | `False` | `{False, True}` | A |
| `bias_tuning_range_meV` | meV | 0.0 | sensitivity when `alignment_tunable=True` | A |
| `occupancy_control_known` | bool | `False` | `{False, True}` | A (honest-failure flag, not a physics input) |
| `second_pair_control_known` | bool | `False` | `{False, True}` | A (honest-failure flag, not a physics input) |

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

Identity / regime: `platform`, `family`, `T_hs`, `T_j`, `thermal_iterations`,
`thermal_converged`, `cycle_loading`, `rep_rate_hz`, `strain_bound`,
`bound_role`, `bound_reversal`, `screening_fraction`, `core_radius_nm`,
`outer_radius_nm`, `conducting_radius_nm`.

Levels/optical (piece 2): `E_X_eV`, `lambda_nm`, `field_kVcm`,
`overlap_sq`, `electron_bound`, `hole_bound`, `E_a_meV`, `gamma_X0_ns`,
`gamma_XX0_ns`, `k_X_ns`, `k_XX_ns`, `escape_prefactor_ns`,
`tau_cap_ps_used`, `reservoir_state_count_e`, `reservoir_state_count_h`,
`sidewall_overlap` (coherence column, piece 2 levels directive round;
`NanowireLevels.sidewall_overlap`).

Photonics (piece 3): `radius_over_lambda`, `V_number`, `beta_HE11`,
`eta_collection_X`, `eta_collection_XX`, `radiative_rate_factor`,
`gamma_X_ns`, `gamma_XX_ns`, `approximation_error`, `single_mode`,
`beta_multimode_penalty` (both piece 3 directive round, H7 fix),
`degree_of_linear_polarization` (piece 3 fix-1 round, non-gating against
`deshpande2013_polarization`'s 70 percent anchor), `antenna_rate_factor`
(piece 3 fix round, commit `5dccc4b`: the wire-antenna screening factor
applied to the radiative rate, orientation-weighted by the card's own
`dipole_weights`; never folded into `gamma` directly).

Surface (piece 4): `k_side_ns`, `k_surface_reservoir_ns`,
`k_surface_X_ns`, `k_surface_XX_ns`, `shell_multiplier_used`,
`reservoir_access_used`, `occupied_dot_access_used`.

Transport (piece 5): `area_cm2`, `J_A_cm2`, `V_j`, `V_terminal`,
`depletion_field_kVcm`, `C_dep_F`, `eta_inj`, `f_capture`, `f_qfl_dot`,
`f_qfl_dot_thermodynamic_limit` (coherence column, piece 5 directive round
H2 -- the delivered-vs-thermodynamic-ceiling QFL pair, see "Composition
rules for the device piece" bullet 3; never averaged or substituted for
`f_qfl_dot`), `f_qfl_background`, `r_supply_s`, `r_captured_s`,
`r_matrix_radiative_s`, `r_matrix_nonradiative_s`, `r_surface_reservoir_s`,
`r_leakage_s`, `mu`, `power_on_W`, `accounting_residual_s`, `tau_RC_ns`,
`delivered_step_fraction`, `pulse_delivery_feasible`.

Lifetimes (device, piece 7): `tau_rad_bare_ns`, `tau_rad_photonic_ns`,
`tau_total_X_ns`.

Counting: established count/flux/background/filter/thermal columns,
transcribed by name from `fsim_core/device.py`'s planar `_evaluate_nitride`
(around line 1333 in the current tree) and reused unchanged for the
nanowire tier where the same quantity applies -- `collected_flux_pulsed_s`,
`collected_flux_x_s`, `collected_flux_xx_s`, `background_flux_s`,
`total_detected_flux_s`, `mean_counts`, `mean_counts_x`, `mean_counts_xx`,
`gate_ns_used`, `S_X`, `S_XX`, `rho_pulsed`, `blocked_load_probability`,
`counting_converged`, and `eta_out` (the filter-transmission key the
planar path uses); plus `g2_op`, `one_pair_valid`, `pair_supply_possible`,
`valid`, `invalid_reasons`, `provenance`; plus the coherence column
`collected_flux_delivered_s = collected_flux_pulsed_s *
delivered_step_fraction` ("Composition rules for the device piece" bullet
6; never folded back into `collected_flux_pulsed_s` itself).

Coulomb screen (unchanged `drive_mech.set_feasibility` convention):
`set_feasible`, `set_E_C_meV`, `set_EC_over_kT`, `set_radius_nm`,
`set_radius_max_nm`, `set_C_sigma_F`, `set_R_T_over_RQ`, `set_f_max_Hz`.

RT injector screen (piece 6, transcribed verbatim): `rti_feasible`,
`rti_status`, `rti_transport_feasible`, `rti_level_margin_kT`,
`rti_orbital_margin_kT`, `rti_alignment_error_meV`, `rti_linewidth_meV`,
`rti_rate_Hz`, `rti_e_rate_Hz`, `rti_h_rate_Hz`, `rti_bypass_fraction`,
`rti_missed_load_probability`, `rti_second_pair_probability`,
`rti_growth_feasible`, `rti_failed_checks`, `rti_evidence_status`.

RT injector screen, fix-1 and directive round additions (module-only:
bound live below by calling `injector_feasibility` at its card defaults
and listing every `rti_`-prefixed key it returns, not by spec text --
`fsim_core/nitride_nanowire_injector.py`'s spec predates all of these):
`rti_level_margin_e_kT`, `rti_level_margin_h_kT`,
`rti_alignment_error_e_meV`, `rti_alignment_error_h_meV`,
`rti_required_bias_shift_meV`, `rti_bypass_fraction_tsai_partition`,
`rti_second_carrier_probability`, `rti_growth_nearest_commensurate_nm`,
`rti_growth_perturbed_margins_kT` (piece-6 fix-1 round), and
`rti_barrier_polarization_tilt_eV`, `rti_well_to_dot_drop_meV`,
`rti_p_free_cm3`, `rti_numerics_ok`, `rti_gate_ns`,
`rti_reservoir_state_count_e`, `rti_reservoir_state_count_h` (the
subsequent directive round, commit `e696bdd`: polarization sheet charges,
alignment against the emitter quasi-Fermi level with the required bias
shift and well-to-dot drop, Mg-acceptor free-hole fraction, a numerics
sanity flag, and a `gate_ns` input separate from `loading_window_ns`).
This directive round landed as `fsim_core/nitride_nanowire_injector.py`
mid-way through this contract fix-2 round (it was reported in flight at
dispatch time); this document was updated to match once it was observed
committed, and `verify/verify_nitride_nanowire_injector.py` itself was
still being updated by that round's own worker when this document was
finished -- see this worker's STATUS `notes:` for the exact `git status`
snapshot to reconcile against.

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
VERDICT: idealized_status=<idealized_status> family=<family> regime=<regime> strain_bound=<strain_bound> bound_role=<bound_role> rep_rate_hz=<rep_rate_hz> complete=<complete> eligible=<eligible> paired_optical_pass=<paired_optical_pass> hardware_qualified=<hardware_qualified> rti_qualified=<rti_qualified> coverage=<n>/<total> invalid=<invalid> flux_floor=1000/s screening=<screening_fraction> access=<occupied_dot_access>
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

Two additional fields, `screening=` and `access=`, are ADDED by this
contract revision per the physics coherence review (fix-1 finding 8 noted
the established `screening=` field was missing from the template
entirely; "Composition rules for the device piece" bullets 8-9 motivate
surfacing the access default alongside it). They report the row's
`screening_fraction` and `occupied_dot_access` card inputs so a reader can
tell, from the VERDICT line alone, which reduced-cut branch and which
surface-access default a given line's numbers used. Unlike the field list
above, `screening=`/`access=` are NOT yet named in piece 9's own field-list
bullet (they postdate that spec's freeze); they are checked directly
against this document's own VERDICT template, not cross-checked against
`nitride-nanowire-sweep.md`.

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
  Ti/Au 5/45 nm end contacts, and ~600 nm WIRE LENGTH are `[V]` (Deshpande
  2013, p.3 device description and p.6 methods). Relabelled this revision
  (fix-1 finding 10): the paper's ~600 nm figure is the total dispersed
  wire length, not the inter-contact (electrically active) length, which
  is unstated and necessarily smaller than 600 nm -- the ledger key is now
  `wire_length_nm` (was `contacted_length_nm`), and no code in this tier
  reads the old key name for an active-length computation. See
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
- `occupied_dot_access`'s headline card default changes `1.0 [A]` ->
  `0.05 [DR]` this revision (fix-2, coherence finding M1): the previous
  `1.0` default caps `tau_X` at 0.625 ns / `tau_XX` at 0.3125 ns at
  `core_radius_nm=12.5` regardless of any other physics, which the held-out
  2013/2014 lifetime anchors (1.1 ns X, 0.71/0.70 ns XX) falsify. `0.05` is
  the J0 hard-wall ground-state probability in a 2 nm sidewall capture
  layer over the uniform-density value at that radius; `1.0` is retained as
  the declared conservative `[A]` partner (never removed), and `0.1` as a
  third sensitivity point. Neither default was tuned to reproduce the
  held-out lifetimes.
- `dipole_weights`' default reverts to isotropic `(1/3,1/3,1/3)` `[A]` in
  the Card schema this revision: a photonics directive round had set the
  default to `(0.0, 0.5, 0.5)` `[DR]` on the geometric reasoning that a
  c-plane exciton's transition dipole lies IN the c-plane, excluding the
  `along_wire` component; the SUBSEQUENT photonics fix round (commit
  `5dccc4b`, landed mid-way through this contract fix-2 round) found that
  prior predicts a NEGATIVE `degree_of_linear_polarization` -- the WRONG
  SIGN against the `deshpande2013_polarization` 70 percent anchor -- and
  reverted the default to isotropic, keeping `(0.0,0.5,0.5)`
  (`CPLANE_ONLY_DIPOLE_WEIGHTS`) as a named, explicitly falsified
  sensitivity. This document was updated to track that reversal once it
  was observed committed; see this worker's STATUS `notes:`.
- `deshpande2013_emission_energy`'s `tolerance.value` changes `0.0` ->
  `null` this revision (fix-1 finding 7): a literal `0.0` tolerance reads
  as "zero error, exact value," conflating "no reported error bar" with a
  measured zero uncertainty. The ledger's own null-is-unknown convention
  (see the opening paragraph of this section) now applies to `tolerance`
  values as well as to `value`, checked by `verify/
  verify_nitride_nanowire_contract.py`'s `check_ledger_rules`.
- New anchor `bleuse2011_claudon2010_beta_envelope` (`[E]`, non-null)
  records the two photonic-wire beta points the photonics directive round
  transcribed to calibrate `beta_multimode_penalty`'s RELATIVE decline
  (beta ~0.95 near d/lambda 0.22-0.24, falling to beta ~0.7 by d/lambda
  ~0.4; Bleuse et al., PRL 106, 103601 (2011), Fig. 2; Claudon et al., Nat.
  Photon. 4, 174 (2010)) -- see `fsim_core/nitride_nanowire_photonics.py`'s
  own `multimode_penalty` provenance string. `maslov2004_he11` and the
  original `claudon2010_extraction` (first-lens extraction efficiency, a
  DIFFERENT quantity from this beta-envelope point) remain separate,
  still-`missing` anchors; this addition does not resolve either gap.

## Evidence ledger schema

`verify/data/nitride_nanowire_anchors.yaml` is a SEPARATE ledger from
`verify/data/nitride_cavity_anchors.yaml` and `verify/data/
nitride_geometry_stark_anchors.yaml` (both unmodified by this piece).
`doi_or_url` uses only identifiers actually printed in `../_goal/
nitride_digests.md` or independently resolvable via Crossref (never a
fabricated DOI). Anchors added since piece 1's first round (round 2, fix-1):

- `deshpande2013_surface_velocity` -- `S_cm_s=1000.0` [E], secondary
  attribution to ref. 35 (unread), Deshpande 2013 p.1.
- `deshpande2013_device_geometry` -- substrate stack, contact metal/
  thickness, wire length, collection objective [V], Deshpande 2013
  p.3, p.6.
- `deshpande2013_polarization` -- 70 percent axial DOLP and its
  dielectric-contrast attribution [V], Deshpande 2013 p.5.
- `deshpande2013_emission_energy` -- X/XX energies, antibinding ordering,
  ambient/junction temperatures [V], Deshpande 2013 Fig. 3c; cross-
  references `deshpande2013-xx-splitting` in the cavity ledger.
- `deshpande2013_coulomb_reference` -- E_C/E_C/kT at r=12.5 nm under both
  the sphere and disc conventions [DR], derived from the 2013 geometry and
  a stated formula, not a number Deshpande 2013 itself reports.
- `encomendero2023_resonant_tunneling` -- injector fix-1 round's AlN/GaN/AlN
  double-barrier resonant-tunnelling benchmark (device structure and I-V
  peak, non-gating); this revision's networked `verify/verify_citations.py`
  run attempted to resolve `arXiv:2303.08352` and hit a transient arXiv
  network timeout (not a not-found result); still a recorded gap
  (`doi_or_url: null`, `evidence_status: missing`) -- see the anchor's own
  `transfer_notes` for the retry procedure.

New this revision (round 3, fix-2):

- `bleuse2011_claudon2010_beta_envelope` -- the two photonic-wire beta
  points (`[E]`, Fig. 2 references) the photonics directive round used to
  calibrate `beta_multimode_penalty`'s relative decline; see "Provenance
  corrections" above.

`deshpande2013_thermal_and_pl`'s `value` no longer carries
`surface_S_cm_s_secondary` (moved to `deshpande2013_surface_velocity` with
its own correct `[E]` tag, so the `[V]`-tagged thermal/PL anchor no longer
mixes a secondary-attributed number into a directly-measured one).
`deshpande2013_device_geometry`'s `value` carries `wire_length_nm` (was
`contacted_length_nm`; see "Provenance corrections" above) --
`verify/verify_nitride_nanowire_contract.py` reads the new key name.

## Round trip and verifier scope

`verify/verify_nitride_nanowire_contract.py` parses this document (markdown
tables and the fenced VERDICT block; a small stdlib+yaml parser, no new
dependency), the ledger, and (new this revision) LIVE-IMPORTS the five
committed production modules (`fsim_core/nitride_nanowire_{levels,
photonics,surface,transport,injector}.py`), and asserts:

- Spec binding (unchanged since round 1): every module/function signature
  cell in the "Module table" appears verbatim in the corresponding spec
  file's "Interface or signature constraints" section.
- Module binding (new this revision, fix-2 required change 5): for every
  row, the Symbol cell resolves via `getattr` on the live, imported module
  (a Symbol rename is now caught, fix-1 finding 12); for a row whose symbol
  is a plain function, the module's own `inspect.signature` (annotations
  stripped -- this project's `from __future__ import annotations` makes
  every module's own annotations opaque source-text strings, not the
  simplified names this table uses -- defaults formatted, whitespace
  normalised) equals the row's Signature cell exactly; for a row whose
  symbol is a `@dataclass`, its `dataclasses.fields()` name set is checked
  against the Card schema instead (see below), since a full constructor
  signature would duplicate the Card schema table without adding
  information. Four rows per fix/directive-round drift (the injector's
  `transmission`/`reflection` carrier keyword, the injector's
  `injector_feasibility` `gate_ns` keyword, transport's
  `evaluate_injection` GaN-reservoir keywords) are marked "module-only" in
  the Verifier cell and are bound ONLY to the live module, not to any spec
  text (see "Module table" above); the device/sweep rows (piece 7/9) are
  bound to spec text only, since `fsim_core/nitride_nanowire_device.py`
  does not exist yet.
- Card-to-dataclass binding (new this revision): for
  `NitrideNanowirePhotonicsParams`, `NitrideNanowireSurfaceParams`, and
  `NitrideNanowireInjectorParams`, the leaf-name set of the corresponding
  single Card schema subsection (`nitride.photonics`, `nitride.surface`,
  `nitride.injector`) equals `dataclasses.fields()` exactly. For
  `NitrideNanowireSystem` (levels) and `NitrideWireDiode` (transport),
  whose card leaves are split across multiple subsections with some
  renaming (`nitride.dot.radius_nm` for `disc_radius_nm`, per H6; the rest
  exact-name), every field name is asserted present SOMEWHERE in the Card
  schema section instead of in one named subsection.
- Row-column binding (new this revision): `injector_feasibility` is called
  once at its card defaults and every `rti_`-prefixed key it returns is
  asserted present in the "Row columns" section text (catches a key this
  document has not yet documented, live, independent of spec text); the
  established count/flux/background/filter/thermal columns named in
  required change 3 are asserted present both in the "Row columns" section
  and, by literal dict-key substring, in `fsim_core/device.py`'s own
  source text.
- Structural checks (unchanged): every card leaf has a unit, a default,
  and a tag in `{V, DR, E, A}`; the sweep grid arithmetic in this document
  equals the product of its own axis counts and equals the independently
  expected 1536/1024; the VERDICT template contains every field piece 9's
  own field list requires, plus (new) the `screening=`/`access=` fields
  checked against this document only; every ledger anchor tagged `V` has a
  non-null `doi_or_url` or a `title`; no `V`-tagged anchor's
  `transfer_notes` reads as unread/secondary; every `value: null` anchor
  has `evidence_status: missing`; no anchor has `value: 0`/`0.0` with
  `evidence_status: missing`; (new) no anchor has `tolerance.value: 0`/
  `0.0` (unknown tolerance is null, never zero, matching the `value`
  convention); every module in the "Module evidence map" has at least one
  anchor with a non-null value; `deshpande2013_supplement` has a null
  `doi_or_url`; and `cavity.enabled` is documented `false` for both
  families.
- Current-density cross-check (rewritten this revision, fix-1 finding 4):
  the 25 nm and 30 nm current densities are recomputed from the RADII READ
  LIVE from `deshpande2013_geometry`'s `optical_diameter_nm` and
  `deshpande2013_thermal_and_pl`'s `thermal_diameter_nm` ledger fields (not
  from a hardcoded area literal), cross-checked against the 203.7/141.5
  A/cm2 figures parsed out of this document's own "Evidence and verdict
  semantics" prose, and the 30 nm figure is additionally cross-checked
  against Deshpande 2013's own independently reported rounded
  `J_1nA_reported_A_cm2=142.0` value.

It prints `N/N nitride nanowire contract checks passed` and exits 0 iff all
pass. A corrupted copy of this document (one module-table function name
changed) is exercised once, in memory, against the same parser to confirm
the substring check fails; the on-disk document is never modified by that
self-test.
