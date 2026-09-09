# Planar InGaN/GaN cavity contract

This is the opt-in contract for the planar InGaN/GaN vertical-cavity,
electrically driven single-photon design tier
(`fsim_core.device.DeviceDesign.platform="ingan_gan_planar"`, implemented by
`fsim_core/device.py`'s `_evaluate_nitride`). It is a model-and-evidence
assessment, not a device demonstration: no InGaN/GaN electrically driven
planar vertical-cavity single-photon source has been published (the closest
published devices are a nanowire dot-in-wire diode, Deshpande et al. 2013/
2014, and separate optical-only planar cavity/QD-on-Si papers -- see
"Evidence status" below). Literature anchors live in
`verify/data/nitride_cavity_anchors.yaml`; `[V]` means verified published
value, `[DR]` derived from published values, `[E]` class estimate, `[A]`
assumption. Copying a verified off-platform measurement into this planar
design makes that TRANSFER `[E]` or `[A]`, never `[V]`, even when the
underlying number is `[V]` in the ledger.

This piece (nitride-cards-contract) ships three design cards
(`cards/nitride-cavity-pulse-design.yaml`, `cards/
nitride-cavity-set-design.yaml`, `cards/
nitride-deshpande2014-comparison-design.yaml`), this contract, the evidence
ledger, and `verify/verify_nitride_cards.py`. No production code
(`fsim_core/`) is edited by this piece; the schema below is the one
`fsim_core/device.py`'s `_evaluate_nitride` already implements (pieces 1-5).
Piece 7 (`nitride-sweep-results`, not part of this piece) consumes the grid
and sensitivity ranges this contract freezes.

## Platform and geometry

`platform: ingan_gan_planar` opts a `DeviceDesign` into the nitride
evaluator. It requires, at `DeviceDesign.load()` time, that `dot.gamma0`,
`dot.a_ac`, `dot.E_LO`, `dot.gamma300`, `dot.r_xx`, `dot.delta_xx` are all
EXPLICITLY present in the card's own YAML (`DotBlock`'s class defaults are
InP class numbers -- a nitride card that leaves any of these unset would
silently inherit a wrong-material default with no error). It also requires,
at `evaluate()` time: `dot.linewidth="anchored"`, `dot.lineshape=
"lorentzian"`, `ret.mode="nitride_confinement"`, `ret.channel` explicitly
resolvable (`"min"` on all three cards here -- `RetentionBlock`'s own class
default `"pair_half"` is REJECTED by this evaluator, which requires a bound
wetting-layer continuum this finite-well model does not carry), `drive.mode=
"EL-transport"` with `drive.diode.preset="nitride-planar"`,
`drive.finite_pulse=true` with an empty `drive.mechanism`, `emission.type=
"vertical_cavity"`, `cavity.enabled=true` with `cavity.type=
"nitride_planar"`, and legacy `cavity.F_P`/`cavity.G` left at their class
defaults (10.0/8.0) -- `nitride.cavity` owns Q/mode-volume/collection
instead. `ret.system`/`ret.preset` (the legacy inline-stack mechanism) must
stay empty.

The `nitride` block carries the planar-specific geometry and background
knobs not covered by the shared `dot`/`ret`/`drive`/`cavity` blocks:
`nitride.dot` (a `fsim_core.nitride_levels.NitrideDotSystem`: `height_nm`,
`radius_nm`, `x_in` required; `wl_thickness_nm`, `strain_fraction`,
`screening_fraction`, `external_field_kVcm`, `vbo_InN_GaN_eV`,
`strain_c_fraction` optional, class defaults apply if omitted),
`nitride.cavity` (a `fsim_core.nitride_cavity.NitrideCavityParams`: `Q`,
`mode_volume_norm`, `spatial_overlap`, `eta_out`, `T_track`,
`dEdT_cav_meV_K`, `detuning_offset_meV`, `purcell_enabled`), and three flat
scalars: `nitride.tau_rad0_ns`, `nitride.k_nr_ns`, `nitride.
background_tau_ns`, `nitride.eta_background`. Any other key under `nitride`
is rejected (`_evaluate_nitride`'s own `allowed` set).

Both headline cards (`cards/nitride-cavity-pulse-design.yaml`, `cards/
nitride-cavity-set-design.yaml`) share IDENTICAL explicit geometry/optics/
material parameters -- `dot`, `ret`, `drive.diode`/`finite_pulse`/
`rep_rate_hz`/`duty`/`gate_ns`/`I_uA`/`cw`/`b_res`, `thermal`, `cavity`,
`emission`, `aperture`, `nitride.dot`, `nitride.cavity`, `nitride.
tau_rad0_ns`/`k_nr_ns`/`background_tau_ns`/`eta_background` -- only
`meta.name`/`design.name`, `drive.cycle_loading`, and the SET pricing
metadata (`drive.eta_load`, `drive.set_params`) differ. This is checked
mechanically (`verify/verify_nitride_cards.py`, AC1). Geometric headline
values (`nitride.dot.height_nm=3.0`, `radius_nm=10.0`,
`thermal.mesa_diameter_um`, aperture geometry, thermal-stack thicknesses)
are `[A]` design choices except where a leaf's own
`design.provenance.sources` entry documents a literature transfer
(`nitride.dot.x_in=0.25` and `dot.delta_xx=-10.0 meV`, both Deshpande et
al. 2013; `dot.E_LO=91.5 meV`, Seguin et al. 2006; `nitride.dot.
vbo_InN_GaN_eV=1.15 eV`, Tsai & Bayram 2020) -- a geometric value is never
promoted to `[V]` merely because it happens to coincide with a published
number.

## Electrical regimes

`drive.cycle_loading` selects one of two full-cycle loading regimes
(`fsim_core.pulse_counting`), deliberately independent of the legacy
`drive.loading_model`/`F8` Poisson machinery:

- `"rectangular"` (`cards/nitride-cavity-pulse-design.yaml`): a real
  finite-electrical-pulse waveform loading calculation (`pulse_counting.
  pulse_g2`), propagating the exact state-resolved factorial-moment
  hierarchy through the pump/dark windows of one 100 ps/80 MHz pulse
  period (12.5 ns), counted over the full period (`drive.gate_ns=null`).
- `"deterministic_pair"` (`cards/nitride-cavity-set-design.yaml`): an
  idealized single-electron-turnstile (SET) one-pair-per-cycle load
  (`pulse_counting.deterministic_cycle_g2`, `drive.eta_load`), priced
  against a hardware feasibility screen (`drive_mech.set_feasibility`/
  `mech_set`, `drive.set_params`: `radius_nm`, `eps_r`, `R_T_ohm`,
  `ec_margin`). The SET island radius is INDEPENDENT of the confined QD
  radius (`nitride.dot.radius_nm`) -- it is a separate classical
  charging-energy island, never a nanowire dimension, and never a
  substitute Fano factor for deterministic photon delivery. The evaluator
  reports BOTH the idealized photon-counting result (`mean_counts`,
  `g2_op` etc., always computed) and the hardware diagnostic
  (`set_feasible`, `set_priced_F_p`, `set_EC_over_kT`, `set_radius_max_nm`,
  `one_pair_valid`, `blocked_load_probability`) side by side; an
  infeasible SET screen does not overwrite the idealized counting result,
  but `scalars["device_pass"]` requires BOTH the idealized gates AND
  `set_feasible`/`pair_supply_possible` for the deterministic regime --
  it cannot report a hardware pass when the island is infeasible.

Both regimes require `drive.diode.tau_pulse_ns` and `drive.rep_rate_hz`
explicit and positive, with `tau_pulse_ns <= 1e9/rep_rate_hz` (the period).
The transport diode (`nitride_transport.planar_pin` ->
`NitrideDiode`) reads `drive.diode` (minus `preset`/`tau_pulse_ns`): both
headline cards set `N_A=1e17`, `N_D=1e18 cm^-3` `[V Zhang 2016 source, E
planar transfer]` (Zhang et al., APL 108, 153102 (2016), their own p-GaN/
n-GaN doping), `d_i_nm=24.5` `[A adapted Zhang geometry]` (their 12 nm
undoped-GaN spacers either side of THIS card's own `nitride.dot.
wl_thickness_nm` reservoir, not their 3 nm QW), `x_in` set EQUAL to
`nitride.dot.x_in`, and `drive.diode.wl_thickness_nm` = 0.5 nm kept as the
diode's SRH reservoir thickness only (since commit 59014bc the confinement
solver has no wetting-layer model and `nitride.dot.wl_thickness_nm` is 0;
the two fields are therefore NOT equal, by design). All remaining `NitrideDiode` fields (`d_active_nm`, `area_um2`,
`R_s_ohm`, `n_ideality`, `tau_SRH_ns`, `tau_cap0_ps`, `tau_matrix_ns`,
`eps_r`, `f_Rs_local`) are left at their module class defaults; the
resolved values are recorded in each card's `design.provenance.
diode_unset_defaults` rather than silently left implicit. `drive.diode.
area_um2` (0.785 um^2 default) is independent of `aperture.diameter_um`
for this evaluator (unlike the legacy `EL-transport` diode path's own
`_diode_from_drive`, which DOES default the injection area to the aperture
area) -- a pre-existing module behaviour, not a defect introduced by these
cards.

## Cavity and tracking

`nitride.cavity` models a planar dielectric-DBR vertical cavity:
`Q`, `mode_volume_norm` (`V/(lambda/n)^3`), `spatial_overlap`, `eta_out`,
`T_track` (fixed at 300 K across the full `T_hs` sweep, per
`nitride-sweep-results.md`), `dEdT_cav_meV_K`, `detuning_offset_meV`,
`purcell_enabled`. Both headline cards assume `Q=2000`, well above the
measured InGaN-QD cavity ceiling (Taylor et al., Nanoscale Res. Lett. 5,
608 (2010): Q=167 best spot, up to 260 across the digest's nitride-cavity
literature) -- flagged in "Known limitations" below, not presented as
demonstrated. `mode_volume_norm` is not reported for any InGaN-QD cavity in
the digest (`[A]`). `purcell_enabled=true` is likewise an assumption Taylor
et al. 2010's own measurement does not support at their much lower Q (they
report NO observed Purcell enhancement). Legacy `cavity.kappa`/`T_track`/
`E_X0`/`dEdT_cav`/`F_P`/`G`/`beta_sin`/`purcell_wire` (plain `CavityBlock`
fields) are inert for this platform except the `F_P`/`G` guard above; they
are left at class defaults and never set on any of these three cards.

## Counting and background

`nitride.tau_rad0_ns` is the field-free radiative lifetime; both headline
cards set `1.0 ns [A]`, deliberately distinct from Deshpande et al., Appl.
Phys. Lett. 105, 141109 (2014)'s own abstract-reported 300 K measured
recombination lifetime (`1.3 +/- 0.3 ns`, `[V]` abstract-only) -- NOT a
transferred measured field-suppressed lifetime (that citation already
backs `fsim_core/nitride_levels.py`'s own `TAU_RAD0_DEFAULT_NS=1.3` module
constant, out of scope for this piece). `nitride.k_nr_ns=0 [A optimistic]`
and `nitride.background_tau_ns=0 [A optimistic]` are both explicit
best-case assumptions with mandatory sensitivity coverage
(`nitride-sweep-results.md`: `k_nr_ns=[0,1]`, `background_tau_ns=[0,1]`).
`nitride.eta_background=0.1 [A]` and `drive.b_res=0.1 [A]` are both round,
explicitly-stated design assumptions -- no InGaN/GaN-specific residual-
background or background-acceptance measurement exists in the digest.

`fsim_core.nitride_transport.evaluate_injection`'s `eta_rad_matrix`
(default `0.1 [A]`), `E_urbach_meV` (default `None`, module resolves its
own `[E]` class value), and `eta_total` (default `0.01 [A]`) are NOT
forwarded from any card field by `_evaluate_nitride`'s own call site (it
does not pass them through) -- they are resolved module defaults, not
accepted nitride card overrides. Each card's `design.provenance.
background_model_defaults` records this explicitly rather than inventing a
new YAML field to expose them.

The retention model uses `ret.channel="min"` (the faster of the electron/
hole 2D-reservoir escape rates) and `ret.tau_cap_ps=10.0 [A]`, with `ret.
tau_cap_scales_with_density=false` -- the shipped "no-cancellation"
convention (matching the RT-edge InP tier's own default,
`docs/rt_edge_contract.md`), not the opt-in "full-cancellation"
alternative. `ret.b0`/`ret.beta` are set to `0.0` explicitly (a no-op
restatement: the nitride evaluator never reads these legacy proxy-
Arrhenius fields) to document "no proxy Arrhenius override" rather than
leaving them silently unset.

A single selected dot: `aperture.density_cm2=1e10 cm^-2` and
`aperture.diameter_um=0.05 um` give an expected dot count under this
aperture of `pi*(0.05/2 um)^2 * 1e-8 cm^2/um^2 * 1e10 cm^-2 = 0.19635`
`[DR from these two A inputs]` (`design.provenance.single_dot_selection` on
each headline card). `drive.n_dot_cm2` is left unset (0), so the transport
dot density falls back to `aperture.density_cm2`.

## Evidence status

The evidence ledger (`verify/data/nitride_cavity_anchors.yaml`) is a
top-level `anchors` MAPPING keyed by anchor id; each entry carries
`citation, doi_or_url, location, evidence_status, tag, platform,
excitation, observable, value, unit, tolerance, transfer_notes`.
`evidence_status` is one of `full_text`, `abstract_only`, `figure_reading`,
`estimate`, `missing`; a missing numerical value is `null`. Every
`tolerance` names what its number means (a reported paper uncertainty, or
an explicit "no reported uncertainty" / "not applicable" for a bare
transcription). `doi_or_url` uses only the identifier `../_goal/
nitride_digests.md` actually prints for that paper (a DOI where stated,
else the arXiv/OSTI/PMC id given) -- never a fabricated DOI.

Retained anchors: Bernardini et al. 1997 (GaN polarization, `[V]`, full
text); Rinke et al. 2008 (GaN electron mass, `[V]`, full text); Wu et al.
~2003 (InGaN bandgap bowing b=1.43 eV, `[V]`, full text); Tsai & Bayram
2020 (InN/GaN VBO=1.15 eV, `[V]`, full text) -- these four back internal
`fsim_core/nitride_materials.py`/`nitride_levels.py` physics (out of scope
for this piece) rather than a single card leaf, except `vbo_InN_GaN_eV`
which both headline cards set directly. Deshpande et al. 2013 (Nat.
Commun. 4, 1675): dot thickness 2.0 nm, diameter 25+/-5 nm, `x_in=0.25`,
Varshni alpha/beta, and the signed XX splitting `-10.0 meV` (antibinding),
all `[V]`, full text, 10 K nanowire dot-in-wire. Deshpande et al. 2014
(Appl. Phys. Lett. 105, 141109, DOI 10.1063/1.4897640): `x_in=0.40`,
`T=300 K`, `g2=0.29`, `lifetime=1.3+/-0.3 ns`, `max repetition=200 MHz`,
all `[V]` **abstract only** -- no primary full text was retrieved as of
this digest's 2026-09-09 compilation (AIP paywalled); dot size, drive
waveform/amplitude, measured count rate, and the repetition rate the
g2=0.29 measurement itself was taken at are all recorded `missing`
(`value: null`) rather than invented. Taylor et al. 2010 (Nanoscale Res.
Lett. 5, 608): cavity Q=167 and "no Purcell effect observed", both `[V]`,
full text. Zhang et al. 2016 (APL 108, 153102): 3 nm active layer,
`x_in=0.15`, `N_A=1e17`/`N_D=1e18 cm^-3`, and BOTH raw (0.42) and
dark-count-corrected (0.38) g2, all `[V]`, full text, planar electrical --
raw and corrected values are never conflated. Seguin et al. 2006
(arXiv:cond-mat/0610425): GaN LO phonon 91.5 meV, `[V]`, full text.

External-context anchors (a different material system and/or pumping
scheme from the headline InGaN/GaN electrical platform; never transferred
onto any card field): Kako et al., Nature Materials 5, 887 (2006), 200 K
triggered GaN/AlN single-photon emission, `[V]` abstract only. Holmes et
al., Nano Lett. 14, 982 (2014), GaN/AlGaN nanowire optical g2=0.13 at
300 K, `[V]` abstract only. Tamariz et al., ACS Photonics 7, 1515 (2020),
optical GaN/AlN g2=0.17+/-0.08 at 300 K and linewidth 55+/-7 meV, `[V]`
full text.

This ledger does NOT adopt the RT-edge InP tier's two-independent-source
policy (`docs/rt_edge_contract.md` "Evidence status") as a gating rule for
this piece: no InGaN/GaN electrical single-dot measurement exists at all
outside the Deshpande/Zhang family, so a two-source requirement would be
unenforceable here. Every anchor's own `evidence_status` and `tag` still
distinguish full-text from abstract-only from missing, and every card
leaf's `tag` still distinguishes a measurement from its transfer.

## Sweep grid and cost

The headline Cartesian grid is FROZEN for piece 7 (`nitride-sweep-
results.md`), both loading regimes: `height_nm=[1,2,3,4,5]`,
`radius_nm=[5,10,15]`, `x_in=[0.15,0.25,0.40]`, `T_hs=[230,250,273,300]`,
`Q=[500,2000,10000]`, `current_uA=[0.002,0.02,0.2]`. All six axes are
`[A]` design axes; 5*3*3*4*3*3*2 = 3240 headline rows. The fixed reference
point when an axis is not being varied is `radius_nm=10`, `x_in=0.25`,
`current_uA=0.02`, `Q=2000`, `height_nm=3`, `T_hs` as swept -- exactly the
two headline cards' own explicit values. `cards/
nitride-deshpande2014-comparison-design.yaml` is explicitly NOT part of
this 3240-row grid (`headline_grid_row_count: 0` in its own provenance).

One-at-a-time sensitivity ranges (fixed reference geometry, NOT another
full Cartesian product, `nitride-sweep-results.md` "Sensitivity rows"):
`screening_fraction=[0,0.5,1]`; `strain_fraction=[0.5,1]`;
`gamma300=[25,35,55] meV`; `delta_xx=[-10,+10,+16] meV`;
`tau_cap_ps=[1,10,100]`; `tau_cap_scales_with_density=[false,true]` at
density 1e9 cm^-2 (including `false` at that same density);
`tau_rad0_ns=[0.2,1,4] ns`; `k_nr_ns=[0,1] ns^-1`;
`background_tau_ns=[0,1] ns`; `b_res=[0,0.1,1]`; `eta_out=[0.01,0.1,0.5]`;
`mode_volume_norm=[1,2,10]`; `detuning_offset_meV=[0,+10] meV`;
`purcell_enabled=[false,true]`; plus an observed-cavity-scale reference
`Q=167, purcell_enabled=false` `[V Q Taylor2010, A device transfer]`. All
are `[A]` envelopes except the separately cited source inputs
(`gamma300=55` is NOT a measured InGaN value -- it is inside the union of
the Wang 2017 220 K extrapolation and the off-platform Tamariz 2020 GaN/AlN
number, never a fit to either). SET feasibility sensitivity: island
`radius_nm=[0.5,1,5] nm` at four `T_hs`, independent of the QD radius;
sub-nm entries are classical isolated-sphere self-capacitance estimates,
explicitly beyond a validated hardware design, not demonstrated
manufacture.

## Acceptance gates

The headline metric is pulsed intrinsic `g2(0)` (`scalars["g2_op"]`), and
the headline brightness metric is collected pulsed flux
(`scalars["collected_flux_pulsed_s"]`) -- following the RT-edge InP tier's
own convention (`docs/rt_edge_contract.md`), reused here: `g2(0) < 0.5`
AND collected flux `>= 1000 photons/s` (the same 1 kHz floor). For the
`deterministic_pair` regime, `scalars["device_pass"]` additionally
requires `one_pair_valid` AND `set_feasible` AND `pair_supply_possible` --
an idealized-optical pass with an infeasible SET island is reported as
exactly that (idealized pass, hardware infeasible), never as an
unconditional device pass. Bound/bright/pure results are NOT card-validity
gates: a physical failure (`scalars["valid"]=false` with populated
`invalid_reasons`) is valid, honestly-reported model output, checked
mechanically by `verify/verify_nitride_cards.py`'s round-trip/evaluate
checks (AC4) -- never silently converted into a passing or a missing row.
Absent optical pumping is checked mechanically (`drive.cw=false`, no
`emission.type` value overlapping with an optical-excitation convention).
This piece ships design cards and their evidence/sweep contract only; it
does NOT itself run or gate the full acceptance sweep (`nitride-sweep-
results.md`, piece 7, not part of this piece's scope) -- no PASS/FAIL
device verdict is asserted here.

## Deshpande comparison

`cards/nitride-deshpande2014-comparison-design.yaml` evaluates the SAME
`ingan_gan_planar` model at Deshpande et al. 2014's own published
conditions where the abstract states them (`x_in=0.40`, `T_hs=300 K`,
`rep_rate_hz=2.0e8`, the paper's own reported MAXIMUM excitation rate --
the abstract does NOT state that the g2=0.29 measurement was itself taken
at that rate, recorded as its own `missing` anchor,
`deshpande2014-g2-measurement-rate-missing`), with finite electrical
loading (`drive.cycle_loading=rectangular`, `drive.finite_pulse=true`),
distinct from the headline 80 MHz grid. Dot size and drive waveform/
amplitude are not published in the retrieved abstract, so the card uses
explicit `[A]` placeholders motivated ONLY by the same group's 2013
nanowire family: `height_nm=2.0` (Deshpande 2013's own dot thickness),
`radius_nm=12.5` (half Deshpande 2013's own 25 nm diameter),
`diode.tau_pulse_ns=0.1`/`I_uA=0.02` (the headline cards' own design
values) -- never presented as verified 2014 geometry or drive. `design.
provenance.comparison_status: conditions_incomplete` records this (a
provenance note, not a new/unsupported design schema field); the ledger
entries for the missing published size, waveform, count rate, and g2
measurement rate are all `null` (`verify/data/nitride_cavity_
anchors.yaml`), never invented -- an unreported measured count rate is
recorded as unavailable, never as zero.

This card produces a REQUIRED, non-gating model evaluation
(`verify/verify_nitride_cards.py` evaluates it, AC4), labelled an
assumed-geometry transfer at the reported maximum rate, NOT a reproduction
of the published experiment: the model's own evaluated g2/flux are never
adjusted, fitted, or forced to match the published g2=0.29, and this
comparison never gates the headline cards' own acceptance verdict.
Deshpande et al. 2013's 10 K nanowire measurement is never presented as
the 2014 device (they are cited as two distinct, separately dated,
separately conditioned anchors). Transfer assumptions this card states
explicitly (`design.provenance.transfer_assumptions`): nanowire strain
relaxation removed; planar polarization/QCSE restored (the nanowire's own
near-zero blueshift-with-current finding does not apply to this
unscreened planar disk); cavity out-coupling in place of the nanowire's
own dielectric-contrast end-facet extraction; assumed optical mode volume/
Q/collection (the same `[A]` headline cavity assumptions, not calibrated
to either Deshpande paper, which reports no cavity at all); unknown drive
waveform/contact heating/background; and the unspecified raw/corrected
status of the published g2=0.29 (the abstract does not say which).

## Artifacts and reproduction

Load and evaluate any of the three cards with `fsim_core.device.
DeviceDesign.load`/`evaluate` (never the unrelated params-style
`fsim_core.card.Card` schema): e.g. `DeviceDesign.load("cards/
nitride-cavity-pulse-design.yaml")` then `evaluate(design,
T_grid=[design.thermal.T_hs])`. `verify/verify_nitride_cards.py` is the
mechanical checker for this piece (round-trip, provenance completeness,
anchor cross-references, the frozen grid, common-geometry equality between
the two headline cards, and one evaluate() call per card at its own
`T_hs`); run it with `python verify/verify_nitride_cards.py` from the
repository root. It prints `N/N ... passed` and exits 0 iff every scored
check passes. This piece ships no sweep script or output directory --
`scripts/run_nitride_cavity.py` and `out/nitride_cavity/` belong to piece 7
(`nitride-sweep-results`, not part of this piece).

## Known limitations

- The measured InGaN-QD cavity Q ceiling is 167-260 (Taylor et al. 2010
  and the digest's own Sec. 5/10d nitride-cavity literature summary); both
  headline cards assume `Q=2000`, well above any published InGaN-QD
  cavity, and are not presented as a demonstrated cavity.
- Taylor et al. 2010 explicitly report NO observed Purcell enhancement at
  their much lower measured Q; both headline cards set `purcell_enabled=
  true`, an assumption that paper does not support at this higher,
  assumed Q.
- No InGaN/GaN electrical single-dot 300 K linewidth measurement exists in
  the digest; `dot.gamma300=35 meV` sits inside the union of the Wang et
  al. 2017 220 K InGaN class extrapolation and the off-platform Tamariz et
  al. 2020 GaN/AlN optical number, but is fit to neither -- do not read a
  hidden fit into this value.
- No InGaN/GaN cavity thermal-conductivity/package-resistance measurement
  exists in the digest; the thermal stack (`thermal.layers`,
  `thermal.substrate`) uses generic `[A]` crystalline-semiconductor-class
  GaN/Si values and a deliberately low-conductivity lumped `[A]` DBR/
  contact effective layer, explicitly not the legacy GaAs class fallback.
- Deshpande et al. 2014's own abstract is the ONLY evidence for the
  comparison card's published-condition fields; no primary full text was
  retrieved as of this digest's 2026-09-09 compilation, so dot size, drive
  waveform, count rate, and the g2 measurement's own repetition rate
  remain unverifiable and are recorded `missing`, never estimated as a
  number.
- The single-band, separable-disk confinement solver
  (`fsim_core/nitride_levels.py`) is a deliberately limited model (see
  that module's own docstring); it is not fit to or validated against any
  InGaN/GaN dot's measured transition energy in this digest.
- No held-out prediction is claimed anywhere in this piece: every
  literature-sourced card leaf is either a direct-platform verified value
  transferred with an explicit `[E]`/`[A]` note, or an external-context
  anchor never transferred onto any card field at all.
