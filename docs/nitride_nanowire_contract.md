# Nitride nanowire contract

This contract freezes the opt-in `ingan_gan_nanowire` tier.  It does not
change `ingan_gan_planar`, legacy evaluation, or the planar cavity contract.
All nanowire predictions are independent predictions: the Deshpande results
below are comparison anchors, never fit targets.

## Families and common card contract

`platform` is `ingan_gan_nanowire`; `emission.type` is `nanowire`;
`drive.diode.preset` is `nitride-nanowire`; and `cavity.enabled` is false.
There is no Q, DBR, cavity-Q axis, or planar aperture supply factor.  Required
blocks are `nitride.nanowire`, `nitride.surface`, `nitride.photonics`,
`nitride.wire_thermal`, and `nitride.injector`; `nitride.dot` retains axial
height and composition. `tau_rad0_ns=1.0` [A, deliberately not the held-out
2014 lifetime] and `drive.b_res=0.1` [A] are explicit.  `b_res=0.02` is only
a sensitivity motivated by the 2013 corrected HBT value.

`nitride.nanowire` leaves are `family` (`horizontal_as_built` or
`vertical_photonic`), `core_radius_nm`, `outer_radius_nm`, `strain_bound`
(`unrelaxed` or `relaxed`), `shell` (`none` or `AlGaN`),
`barrier_left_nm`, and `barrier_right_nm`.  `nitride.dot.radius_nm` equals
`core_radius_nm`; `outer_radius_nm >= core_radius_nm`; and shell `none`
requires equality. The disc fills the semiconductor core.  The independently
explicit `drive.diode.conducting_radius_nm` is in `(0, core_radius_nm]` and
defaults to the core; it is never an optical radius alias. A radius sweep
updates dot/core/outer/conducting radii, SET radius, and
`d_i_nm=barrier_left_nm+height_nm+barrier_right_nm` together.

The two endpoints are scenarios, not a universal relaxation curve:
unrelaxed has `strain_fraction=1`, relaxed has `strain_fraction=0`. Both
change strain deformation shifts and `P_pz`; `P_sp` remains. Screening is a
separate `screening_fraction=0` main-grid input, with `{0,1}` only on a
reduced cut. Unrelaxed is labelled conservative/lower and relaxed
headline/upper, without asserting monotonic flux ordering. The radial model
is a hard-wall cylinder [A] with GaN axial reservoirs; it must not model a
full-core disc as laterally surrounded by infinite GaN. Dielectric self
energy, sidewall band bending, elastic spatial variation, and lateral alloy
localization remain unresolved.

## Family-specific interfaces

Horizontal-as-built means the dispersed, end-contacted wire on SiO2/Si,
normal-incidence objective collection, no mirror/taper/top contact, and a
subwavelength antenna envelope. `NA=0.5` and its collection scale are [A].
The reported 70 percent axial polarization constrains orientation only; it
is not collection efficiency. Vertical-photonic means a designed HE11 wire:
`n_wire`, `n_ambient`, taper transmission, bottom reflectivity, top-contact
transmission, propagation transmission, and collection envelope are explicit
[A]/[E] inputs. HE11 confinement, beta, first-lens extraction, and radiative
rate factor remain separately normalized; baseline radiative-rate factor is
one [A]. Beta is both propagation directions, and mirror return is
incoherent. A high beta never changes the radiative rate by itself.

`nitride.surface` declares `S_cm_s=1e3` [V secondary attribution in
Deshpande 2013, ref. 35], shell multiplier [A], and whether reservoir or
occupied-dot access is being priced. Surface loss is `2S/R` only in its
declared channel and is not counted twice. `nitride.wire_thermal` declares
series resistance (2.38 GOhm is [V] for the 2013 device), thermal resistance
[A], and duty averaging. `nitride.injector` declares barrier material,
thickness, alignment, tunnelling and thermionic controls; unsupported
second-pair/hole controls produce an honest failed RT screen.

## Evidence and verdict semantics

The separate ledger is `verify/data/nitride_nanowire_anchors.yaml`. Its
required fields are citation, doi_or_url, location, evidence_status,
evidence_kind, tag, platform, excitation, observable, value, unit,
tolerance, transfer_notes, temperature, geometry, raw_corrected, and
observable_definition. Status is `full_text`, `abstract_only`,
`figure_reading`, or `missing`; unknown values are null, never zero.

2013 is CW electrical: X g2 raw/corrected `0.30/0.16`, XX `0.38/0.25`, at
1 nA. Its TRPL XX decay is 711 ps while HBT-fit X/XX correlation times are
1.1/0.7 ns; these are not interchangeable. Optical diameter is 25 +/- 5 nm;
the paper's thermal/current-density simulation cross section is 30 nm.
Thus 1 nA is 203.7 A/cm2 at 25 nm and 141.5 A/cm2 at 30 nm [DR], versus its
rounded 142 A/cm2. The reported 15/49 K rises are thermal-simulation
estimates, not thermometry or universal thermal resistance. The 52 percent
10-to-300 K ensemble PL ratio is a loss comparison, not an S or lifetime fit.
The 2014 record is abstract-only: x=0.40, about 630 nm, 300 K, g2=0.29,
1.3 +/- 0.3 ns, up to 200 MHz. Geometry, shell, contact layout and waveform
are inherited assumptions, not verified 2014 details. Supplement retrieval
was not established and is marked missing.

`optical_pass` is valid AND converged AND `g2_op < 0.5` AND collected flux
at least 1000/s, plus one-pair validity and supply for deterministic rows.
`hardware_qualified` is an optical SET pass with the unchanged Coulomb
screen; `rti_qualified` is an optical SET pass with the independent injector
screen. Pulse hardware fields are not applicable and false. `device_pass`
uses the Coulomb screen; `rti_device_pass` is separate. Neither failure
overwrites idealized optical statistics.

The main sweep is 1536 horizontal plus 1024 vertical rows: radii
`{10,12.5,15,20,25,40}` and `{60,80,100,120}` nm respectively, heights
`{1.5,2,3,4}`, x `{0.25,0.40}`, T `{230,250,273,300}` K, two regimes, two
bounds and `{80,200}` MHz. It permits no more than 10000 evaluations.
