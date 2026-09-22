# COMSOL work order for qd-photon-sim

Written 2026-09-22 by the orchestrating session for a COMSOL agent that has never seen this project. Read sections 0 to 2 in full before opening any model. Sections 3 to 6 are the simulation items, one per numbered heading; each is self-contained. Section 7 is what COMSOL cannot settle. Section 8 is the delivery checklist.

## 0. What this project is and why COMSOL is being asked

`qd-photon-sim` (folder `C:\Users\rama4\Downloads\quantum dot\qd-photon-sim`, git branch `rt-edge-emitter`, never pushed) is a Python simulator (numpy/scipy) that predicts whether an electrically driven, self-assembled quantum-dot diode can work as a single-photon source at elevated temperature. It computes the second-order correlation at zero delay g2 and the collected photon flux, and gates a design on g2 below 0.5 and at least 1000 collected photons per second, over heat-sink temperatures 230 to 300 K. Every literature number in the code carries a provenance tag: `[V]` verified against a cited paper, `[DR]` derived from published data, `[E]` class estimate, `[A]` assumption (see `CLAUDE.md`).

The simulator is built from analytic and one-dimensional models. Each module states its approximations in its docstring. Several of those approximations decide the verdicts, and several docstrings say outright that a finite-element field solve should replace them ("COMSOL replaces it at Phase-1 exit" in `fsim_core/thermal.py`; "not a replacement for a full-vector FEM calculation" in `fsim_core/waveguide.py`; "[A -> computed from MEEP/COMSOL]" in `fsim_core/cavity.py`; "that is the later MEEP/COMSOL job" in `fsim_core/spec.py`). This work order lists every such approximation that COMSOL's headless batch mode can check, tier by tier, with the geometry, the quantity to extract, the number the Python model currently uses, and the pass criterion.

The four device tiers and their current verdicts:

| Tier | Device | Where | Verdict today |
|---|---|---|---|
| 1 legacy | InAs/GaAs class proxy, cryogenic, no cavity | `cards/staged-device-design.yaml`, `cards/chatzarakis*.yaml` | calibration tier (fits), no gate |
| 2 rt-edge | InP dot in a GaAs0.6P0.4 well, AlGaInP claddings on GaAs, ridge waveguide, edge emission, 230 to 300 K | `docs/rt_edge_tier.md`, `out/rt_edge/verdict.md` | FAIL: no corner clears 1000 photons/s (flux_max 583/s) |
| 3 planar nitride | In0.25Ga0.75N/GaN disc in a planar p-i-n mesa with a dielectric DBR microcavity; round 2 adds nonpolar orientations and a bias (Stark) prediction | `docs/nitride_cavity_contract.md`, `docs/nitride_geometry_stark_contract.md`, `out/nitride_cavity/results.md`, `out/nitride_geometry_stark/results.md` | pulse regime FAIL; one-pair-per-cycle (SET) regime passes optically (best 28 kHz at 300 K) but is hardware-infeasible (charging energy far below 10 kT) |
| 4 nanowire | InGaN disc in a GaN nanowire p-i-n, as-built horizontal wire on SiO2/Si (Deshpande 2013/2014) and a designed vertical photonic wire | `docs/nitride_nanowire_contract.md`, `out/nitride_nanowire/full/results.md` | pulse regime FAIL; SET regime passes optically (millions of photons/s at the relaxed-strain bound) but is hardware-infeasible on both the Coulomb and the resonant-tunnelling screens |

How to read a module's current number: run Python from the repository root, e.g. `python -c "from fsim_core import nitride_materials as m; print(m.binary('GaN'))"`. Every module has a verifier `verify/verify_<module>.py` that prints the published anchors it checks against; run it to see the numbers the module is already held to. Never edit `fsim_core/`, `verify/`, `cards/` or `docs/` from the COMSOL side; results come back as files under a new folder (section 2).

## 1. What "COMSOL CLI" means for this work order

All work is headless. No GUI session, no interactive model building. Two routes are allowed:

1. **Java Model API files.** Write `Model.java` using the COMSOL Model API (`ModelUtil.create`, `model.component().create(...)`, `model.geom()`, `model.material()`, `model.physics()`, `model.mesh()`, `model.study()`, `model.result().export()`), compile with `comsol compile Model.java` (produces `Model.class`), run with `comsol batch -inputfile Model.class -outputfile Model_out.mph -batchlog run.log`. Every model in this work order is authored this way so it is reproducible from text.
2. **Existing .mph with parameters.** `comsol batch -inputfile model.mph -outputfile out.mph -study std1 -pname R -plist "10,12.5,15" -batchlog run.log` for parametric sweeps of a model already built by route 1.

Conventions that apply to every item:

- Export numbers with `model.result().export()` nodes of type `Data` (point/line/global evaluation) to CSV; never read numbers off plots. One CSV per item, columns named exactly as in the item's "Extract" list, one row per parameter combination, plus a `mesh_level` column.
- Mesh convergence: run every item at two mesh levels (the production mesh and one uniform refinement). Report both; the item passes convergence when the extracted quantities change by less than 2 percent (energies: less than 1 meV). Record the element count and the run time.
- Units: SI in the model; report in the units named in each item (nm, meV, kV/cm, ns, K, photons/s).
- Licence modules used: Wave Optics (`emw`, frequency domain and mode analysis), Semiconductor (`semi`, drift-diffusion with Fermi statistics and incomplete ionization; and the Schrodinger Equation interface `schr` with the Schrodinger-Poisson multiphysics coupling), Heat Transfer in Solids (`ht`), Solid Mechanics (`solid`, with initial strain), AC/DC Electrostatics (`es`), Electric Currents (`ec`). If a module is not licensed, mark the item `BLOCKED: licence` in the note and move on; do not substitute a different physics.
- Materials: take the constants from the project's own modules so the comparison is like-for-like. Nitrides: `fsim_core/nitride_materials.py` (wurtzite GaN / InN / AlN; Bernardini 1997 polarization; Rinke 2008 masses; Ioffe elastic constants; Wu 2003 bowing). Arsenides and phosphides: `fsim_core/materials.py` (Vurgaftman 2001 model-solid theory). Refractive indices: `fsim_core/nitride_nanowire_photonics.py` (GaN Barker-Ilegems Sellmeier, SiO2 Malitson, Si tabulated n,k 380 to 750 nm) and `fsim_core/materials.MATERIAL_EXTRA` (Schubert 1995 for the phosphides). Quote the constants you use in the note with their tags.
- Never tune a model to the published measurement it is compared with. The project's rule is that measurements are held-out comparisons.

Wurtzite constants for quick reference (from `nitride_materials.py`, tags in the module): GaN a 3.189 A, c 5.186 A, C13 106 GPa, C33 398 GPa, e31 -0.49 C/m2, e33 0.73 C/m2, P_sp -0.029 C/m2, eps_r 10.28; AlN a 3.112, c 4.982, C13 99, C33 389, e31 -0.60, e33 1.46, P_sp -0.081, eps_r 10.31; InN a 3.54, C13 92 to 121, C33 182 to 224, e31 -0.57, e33 0.97, P_sp -0.032, eps_r 14.61. InGaN by linear (virtual-crystal) interpolation with the Wu bowing on the gap; InN/GaN valence-band offset 1.15 eV (Tsai 2020, the module's [V] value); GaN/AlN valence offset 0.30 eV in the materials module (Tsai) versus 0.70 eV used by the injector module (Martin 1996): run both where the item says so.

## 2. Where results go and how they enter the project

Create `qd-photon-sim/comsol/<tier>/<item-id>/` with: `model.java` (the source you compiled), `params.csv` (the sweep table), `results.csv` (extracted quantities), `convergence.csv` (both mesh levels), `run.log`, and `note.md` (one page: what was solved, boundary conditions, mesh, licence modules, run time, the comparison table against the module value named in the item, and a one-line PASS / DEVIATION / BLOCKED statement with the reason). Do not put .mph files under git if they exceed 50 MB; keep them beside the folder and record their hash in the note.

Each results row becomes a ledger anchor in `verify/data/<tier>_anchors.yaml` with `evidence_kind: numerical_reference`, tag `[DR]` (it is derived from a field solve, not measured), or `[V]` only when the item reproduces a published measurement. The orchestrator writes the anchors and a verifier `verify/verify_comsol_anchors.py` that compares each module output with its anchor at the item's pass criterion; the COMSOL agent only delivers the folders. Item ids are used as anchor ids.

## 3. Tier 2, room-temperature edge-emitting InP dot (branch rt-edge-emitter)

Geometry (from `cards/edge-inp-gaasp-design.yaml`; print any leaf with `python -c "import yaml;print(yaml.safe_load(open('cards/edge-inp-gaasp-design.yaml'))['design'])"`): InP quantum dot, height 4.0 nm, radius 12.0 nm, in a GaAs0.60P0.40 well; inner cladding (Al0.50Ga0.50)0.51In0.49P; outer cladding (Al0.70Ga0.30)0.51In0.49P; GaAs substrate; ridge waveguide with the width and etch depth of the card's `waveguide` block; edge emission at the model's 815.7 nm (the HKUST reference device emits at 750 to 755 nm, a known 9 percent model deviation); p-i-n diode with the card's dopings and a mesa; heat sink at the substrate back; 100 ps current pulses at 80 MHz; T_hs 230 to 300 K. Current verdict FAIL on flux (583 photons/s maximum); the two levers are the waveguide collection and the thermal/electrical operating point.

### E1. Ridge waveguide full-vector mode (Wave Optics, mode analysis, 2-D cross-section)

- Replaces: `fsim_core/waveguide.py`, a scalar 1-D slab effective-index solve with a Gaussian far-field; its confinement factor Gamma_dot, the guided-mode Purcell factor F_wg (Lecamp 2007 form) and the NA collection eta_NA set the tier's flux. A collection bug in this module already moved flux by 4.46x once (`docs/rt_edge_tier.md`).
- Model: 2-D cross-section of the ridge (layers and widths from the card; indices from `materials.refractive_index` at 815.7 nm and at 750 nm), scattering boundaries or PML around, mode analysis for the fundamental TE and TM modes.
- Sweep: ridge width and etch depth over the card's sensitivity range; both wavelengths.
- Extract: effective index, group index (from two wavelengths 5 nm apart), effective mode area at the dot position, the fraction of the mode energy inside the well layer (Gamma), the mode field at the dot centre normalised to the mode power (for beta), and the far-field half-angles in both planes.
- Compare to: the module's `slab_modes` / `effective_index_ridge` / `beta_factor` / `na_collection` outputs for the same stack (run `python -c "from fsim_core import waveguide as w; print(w.edge_emission(w.hkust_ridge_stack(815.7), <width>, <etch>, 815.7, <L>, <NA>))"` with the card's values).
- Pass: n_eff within 0.01, Gamma within 10 percent, far-field half-angles within 15 percent, beta within a factor 1.3.

### E2. Facet transmission and NA collection (Wave Optics, frequency domain, 2-D or 3-D)

- Replaces: `facet_transmission` and `na_collection` in `waveguide.py` (Fresnel facet with an optional coating; Gaussian far-field integrated over the objective cone) and `facet_escape_fraction` (ridge loss 5/cm [E], back-facet reflection, dot position).
- Model: the E1 mode launched at a cleaved facet into air (and with the card's coating if any), PML outside; far-field transform.
- Extract: facet power transmission, far-field pattern, fraction within NA 0.5 and NA 0.7.
- Compare to: the module's facet transmission and eta_NA for the same NA; pass within 15 percent.

### E3. Junction self-heating of the mesa (Heat Transfer in Solids, 3-D or axisymmetric)

- Replaces: `fsim_core/thermal.py`, an analytic spreading-resistance network (pillar, cone, substrate half-space) with a fixed-point junction temperature; the docstring names COMSOL as its replacement. Feeds T_j into every temperature-dependent rate.
- Model: mesa of the card's diameter and height on the epitaxial stack and a 350 um GaAs substrate, heat-sink boundary at the back at T_hs, the on-state power of the card's operating point applied as a volumetric source in the junction layer, duty-averaged (0.1 ns times 80 MHz) and also as the instantaneous on-state power.
- Sweep: mesa diameter over the card's range, T_hs 230 and 300 K, thermal conductivities from the module.
- Extract: junction temperature rise, thermal resistance, the 1/e thermal time constant (transient study).
- Compare to: `thermal.py`'s T_j - T_hs for the same inputs; pass within 20 percent. If the transient time constant is shorter than the 12.5 ns period the duty-averaged treatment is wrong and the note must say so.

### E4. p-i-n injection and junction field (Semiconductor, 1-D and 2-D axisymmetric)

- Replaces: `fsim_core/transport.py` (abrupt-junction depletion electrostatics, 3-D effective density of states, SRH / radiative / diffusion current channels with mobilities 2000 / 100 cm2/Vs and 1 ns lifetimes, all [E]), which sets the injected current per dot, the junction field and the background photon rate.
- Model: the card's p-i-n stack with the InP dot layer treated as a thin well (no dot); drift-diffusion with Fermi statistics, SRH with the module's lifetimes, radiative recombination with the module's B, at 230 and 300 K; then 2-D axisymmetric with the mesa and the current-confinement aperture the card declares.
- Extract: I-V from 1 nA to 1 mA, junction voltage at the card's current, depletion field at the well, injection efficiency (fraction of current recombining in the well), lateral current spreading (fraction of current inside the aperture radius).
- Compare to: `transport.py`'s V_j, field and eta_inj at the same point; pass within 0.1 V, 20 percent, 20 percent. The lateral spreading number is new information: the module assumes the aperture partition.

### E5. Dot confinement with a real strain field (Solid Mechanics then Schrodinger Equation, 3-D or axisymmetric)

- Replaces: `fsim_core/dot_levels.py`, a separable disc model (finite square well in z, finite circular well in-plane, adiabatic decoupling) that its docstring says "honestly misses by 100 to 300 meV" for large-mismatch dots; it sets the escape activation energies that shape g2 versus temperature.
- Model: axisymmetric InP lens or disc (4 nm by 12 nm radius) in GaAs0.6P0.4; solid mechanics with the lattice-mismatch initial strain, free top surface; then the single-band Schrodinger equation for electron and hole on the strain-shifted band edges (deformation potentials and offsets from `materials.py`; no k.p mixing available, state this).
- Extract: electron and hole ground-state energies, transition energy, escape depths to the well continuum, in-plane excited-state spacing, hydrostatic and biaxial strain at the dot centre.
- Compare to: `dot_levels.py` levels for the card's dot; pass when the transition energy is within 50 meV and escape depths within 30 percent. A larger miss is a result, not a failure: it quantifies the docstring's own caveat.

## 4. Tier 3, planar InGaN/GaN cavity dot (rounds 1 and 2)

Geometry (from `cards/nitride-cavity-pulse-design.yaml` and the round-2 cards): In0.25Ga0.75N/GaN disc, height 3 nm (1 to 10 nm swept), radius 10 nm (5 to 30 nm), on a c-plane GaN p-i-n with N_A 1e17 and N_D 1e18 cm-3 and an intrinsic region 24.5 nm; mesa 20 um; vertical dielectric DBR microcavity with an assumed Q of 2000 (the measured InGaN-dot cavity ceiling is 167 to 260); GaN on Si; T_hs 230 to 300 K; round 2 adds a-plane / m-plane (nonpolar) orientations and quantum-well-fluctuation dots, and a bias-dependent Stark-shift prediction (slopes -12.95, -9.55, +2.88 meV/V for screening 0, 0.5, 1 at 3 nm; Zhang 2016 measured about -10 meV/V).

### P1. Dielectric DBR microcavity: Q, mode volume, Purcell, out-coupling (Wave Optics, eigenfrequency, 2-D axisymmetric)

- Replaces: `fsim_core/dbr.py` (1-D normal-incidence transfer matrix, "SCOPE WALL: no lateral confinement, sidewall loss, absorption") and `fsim_core/nitride_cavity.py` (Purcell factor from an assumed normalised mode volume, spatial overlap, out-coupling efficiency and cavity temperature coefficient, all [A]; Q 2000 [A] against Taylor 2010's measured 167). The cavity response sets every planar flux number.
- Model: the card's DBR pairs (indices and thicknesses from the card's cavity block, GaN spacer with the dot layer, mesa diameter 20 um and also 2 um), axisymmetric eigenfrequency solve with PML; then a dipole excitation at the dot position for the emitted-power spectrum.
- Extract: resonant wavelength, Q, mode volume in units of (lambda/n)^3, the Purcell factor for a dipole at the dot position, the fraction of the emitted power leaving through the top DBR into NA 0.5, the resonance shift per kelvin with the module's dn/dT.
- Compare to: the card's Q, mode_volume_norm, spatial_overlap, eta_out and dEdT_cav; pass when Q is within a factor 2 of the assumed value (or the note states the achievable Q), mode volume within 30 percent, eta_out within 20 percent.

### P2. Polarization field, screening and Stark shift in the planar disc (Semiconductor, Schrodinger-Poisson, 1-D and axisymmetric)

- Replaces: `fsim_core/nitride_levels.py` (1-D finite-volume BenDaniel-Duke solver with the polarization field confined to the disc and a screening_fraction scenario knob) plus `fsim_core/nitride_stark.py` (Sze-Ng depletion electrostatics converting bias to field). The screening fraction is the single largest lever in the planar tier: it moves the SET flux from 0.83 to 400,000 photons/s at 230 K.
- Model: c-plane GaN p-i-n with the InGaN disc as a strained layer, spontaneous and piezoelectric sheet charges at the interfaces, Poisson with the card's dopings under bias from -2 to +3 V, Schrodinger for electron and hole in the disc; then with a sheet of free carriers in the disc (0, 0.5, 1 times the polarization charge) to reproduce the screening scenarios; repeat for the nonpolar case (polarization sheet charge zero).
- Sweep: disc height 1, 3, 5, 10 nm; x_in 0.25 and 0.40; bias; screening.
- Extract: field in the disc versus bias, electron-hole transition energy versus bias (the Stark slope in meV/V), overlap versus bias, escape depths, the bias at which the built-in field is compensated.
- Compare to: `nitride_stark.py`'s pinned slopes and `nitride_levels.py`'s E_X, overlap and escape depths at screening 0, 0.5 and 1; pass within 2 meV/V on the slope and 20 meV on E_X. Any COMSOL-derived relation between injected carrier density and screening becomes the replacement for the scenario knob.

### P3. Mesa self-heating (Heat Transfer in Solids, axisymmetric)

- As E3, for the 20 um mesa on GaN on Si with the card's DBR stack as a thermal barrier; compare to `thermal.py` at the planar operating point; pass within 20 percent.

### P4. Planar p-i-n injection and aperture partition (Semiconductor, axisymmetric)

- As E4, for the nitride stack with Mg incomplete ionization (170 meV) and the card's current-confinement aperture; compare to `nitride_transport.py`'s eta_inj, depletion field and the aperture supply fraction; pass within 20 percent.

### P5. Charging energy of the injection island (Electrostatics, 3-D)

- Replaces: `drive_mech.set_feasibility`, an isolated-sphere capacitance E_C = e^2/C with C = 4 pi eps0 eps_r R, which is the sole reason every SET regime in tiers 3 and 4 is hardware-infeasible (E_C/kT 1.2 to 1.5 at a 5 nm island; 0.47 at the 12.5 nm disc).
- Model: a conducting island of radius 0.5, 1, 5, 12.5 nm in the actual epitaxial environment (GaN eps_r 10.28, the doped contacts at the card's distances, an InGaN disc of the tier's dimensions), self-capacitance from the electrostatics solve with the island held at 1 V.
- Extract: total capacitance and E_C in meV, its ratio to kT at 230 and 300 K.
- Compare to: the module's E_C (12.126 meV at 12.5 nm); pass within 30 percent. A capacitance larger than the isolated-sphere value strengthens the infeasibility verdict; a smaller one would be a real finding.

## 5. Tier 4, nitride nanowire (round 3)

Reference geometry (the Deshpande 2013 as-built wire, `cards/nitride-nanowire-horizontal-*-design.yaml`): GaN wire, core radius 12.5 nm (also 15 nm), In0.25Ga0.75N disc 2 nm thick filling the core (x 0.40 for the 2014 replay), 15 nm undoped GaN each side, n-GaN 300 nm (Si 3e18 cm-3), p-GaN 250 nm (Mg 5e17 cm-3), wire lying on 100 nm thermal SiO2 on (001) Si, Ti/Au (5/45 nm) end contacts, wire length about 600 nm, emission 446 nm (x 0.25) and 543 nm (x 0.40 relaxed strain). Vertical family (`cards/nitride-nanowire-vertical-*-design.yaml`): the same 12.5 nm disc inside a 60 to 120 nm GaN core, bottom mirror, tapered top, 1 MOhm contact. Verdicts: the SET regime passes optically at the relaxed-strain bound (vertical 128 of 128 rows at 80 MHz, best 5.98e6 photons/s at 230 K; horizontal 58 of 192, best 2.28e6) and fails both hardware screens; the pulse regime fails on re-excitation.

### N1. Strain relaxation and polarization field in the disc (Solid Mechanics then Electrostatics, 3-D or axisymmetric)

- Replaces: the relaxed-versus-unrelaxed strain endpoints in `nitride_nanowire_levels.py` (74.5 kV/cm spontaneous-only versus 4038 kV/cm pseudomorphic at x 0.25; 112.8 versus 6469 at x 0.40). The model has no radius-dependent relaxation law; the whole round is presented as bounds. Deshpande measured under 2 meV blueshift with current, which implies a field near the relaxed endpoint.
- Model: axisymmetric wire with the InGaN disc, initial strain from the lattice mismatch, wurtzite stiffness, free sidewall; then electrostatics with the polarization as body and interface charge (P_sp discontinuity plus e31, e33 times the computed strain), potential pinned in the doped reservoirs so the undoped stack carries zero net drop.
- Sweep: core radius 10, 12.5, 15, 20, 25, 40 nm; disc thickness 1.5, 2, 3, 4 nm; x_in 0.25 and 0.40; no shell and a 3 nm AlGaN shell.
- Extract: volume-averaged in-plane strain in the disc relative to the pseudomorphic value (relaxation fraction), axial field at the disc centre, axial field profile.
- Compare to: the two endpoints; pass criterion is a relaxation fraction and field versus radius from which a `strain_fraction(R, h)` law is derived and tagged [DR]. If the 12.5 nm disc retains under 30 percent of the pseudomorphic strain the relaxed headline stands; if over 70 percent the conservative bound becomes the headline and the sweep is re-run.

### N2. Disc-in-wire single-particle states (Schrodinger Equation, axisymmetric, on the N1 field)

- Replaces: `nitride_nanowire_levels.py` (hard-wall radial envelope for the full-core disc; finite-barrier BenDaniel-Duke disc for the vertical family; adiabatic decoupling of axial and radial motion; analytic Coulomb envelope).
- Model: electron and hole on the N1 field profile with the InGaN/GaN offsets (VBO 1.15 eV), effective masses from `nitride_materials.py` (in-plane and axial), hard sidewall for the full-core disc, finite lateral barrier for the 12.5 nm disc in a 60 to 120 nm core.
- Extract: ground-state energies, transition energy with the same analytic Coulomb envelope the module applies (state it), electron-hole overlap, escape depths to the GaN reservoirs, first transverse excited state.
- Compare to: relaxed E_X 2.777 eV (446.4 nm) at x 0.25 and 25 K, overlap 0.94, transverse spacings 13.93 meV (electron) and 3.31 meV (hole); vertical 12.5 nm disc in a 100 nm core E_perp 7.98 meV (the hard wall would give 9.05). Pass: transition energy within 15 meV, overlap within 5 percent, spacings within 10 percent.

### N3. Horizontal wire on SiO2/Si: radiative rate, collection, polarization (Wave Optics, frequency domain, 3-D)

- Replaces: `nitride_nanowire_photonics.py`'s horizontal family, a point dipole over an air/SiO2/Si layered substrate (Fresnel image, objective cone NA 0.5) with a quasi-static cylinder screening factor (2/(n^2+1))^2 on the transverse dipole components; the wire body is not in the electromagnetic model. These set the horizontal collection 0.079 (450 nm) and the antenna rate factor 0.385 (isotropic dipole) or 0.077 (c-plane-only dipoles), and the predicted degree of linear polarization +84 percent (isotropic) or -100 percent (c-plane-only) against the measured +70 percent along the wire axis.
- Model: electric point dipole at the wire centre, GaN cylinder 25 nm diameter and 600 nm long lying on 100 nm SiO2 on Si (complex index from the module's table), Ti/Au pads at the ends, PML; 446 and 543 nm; the three dipole orientations (along the wire, transverse in the substrate plane, vertical).
- Extract: total radiated power for each orientation relative to the same dipole in bulk GaN and in vacuum (the rate factor), far-field power within NA 0.5 above the substrate (collection), and the degree of linear polarization of the collected light for an isotropic population and for the c-plane-only population (equal transverse and vertical, zero axial).
- Compare to: the numbers above; pass when rate factors are within a factor 1.5 and collection within 1.3. The sign of the c-plane-only polarization is decisive: the module predicts the wrong sign against the measurement, which is why the headline uses an isotropic prior; a COMSOL result with the measured sign for c-plane dipoles would overturn that decision.

### N4. Vertical photonic wire: modes, beta, extraction (Wave Optics, mode analysis then frequency domain, axisymmetric)

- Replaces: the vertical family of `nitride_nanowire_photonics.py`: a Marcuse Gaussian surrogate for the HE11 confinement, a beta estimator from the confinement and a group-index ratio (1.02 [A]), a multimode penalty s(V) [E from Bleuse 2011 and Claudon 2010], an LP11 cutoff at V 2.405 that decides which rows are headline-eligible, mirror, taper and contact factors.
- Model: GaN cylinder in air, radius 60 to 120 nm, 446 and 543 nm: mode analysis for HE11 and the next mode (effective index, group index from two wavelengths, confinement); then a dipole in the wire with a bottom Au mirror (and a DBR variant), a linear taper to a 1.5 um top, PML above.
- Extract: beta (fraction of emission into HE11 both directions), effective and group index, the radius at which the second guided mode appears, extraction into NA 0.5 and 0.9 with and without the taper.
- Compare to: beta 0.81 to 0.86 near d/lambda 0.24, the module's single-mode boundary (100 nm core is above cutoff at 543 nm, 80 nm is not); pass when beta is within 0.1 across the radius sweep and the cutoff radius within 5 nm.

### N5. Electro-thermal operating point of the as-built diode (Semiconductor, Heat Transfer, Electric Currents; 2-D axisymmetric wire plus a 3-D thermal model of the lying wire)

- Replaces: `nitride_nanowire_transport.py`: a GaN homojunction kernel giving V_j 3.09 V at 300 K and 3.51 V at 10 K, a wire thermal resistance of 2.8e9 K/W fitted [DR] to Deshpande's Fig. 4 junction rises, an SRH lifetime refit to 0.01 ns [A], an RC pulse-delivery diagnostic (tau_RC 1.65 ns, delivered step fraction 0.04 for the 2.38 GOhm as-built contact).
- Model: axial p-i-n drift-diffusion with Fermi statistics and Mg incomplete ionization, SRH and radiative recombination, the disc as a heterojunction well; a 3-D thermal model of the wire on SiO2/Si with the Ti/Au pads as sinks and the substrate at the bath temperature; Joule heating from the 2.38 GOhm series path.
- Extract: I-V 0.5 to 5 nA at 10 and 300 K; V_j at 1 and 2 nA; junction rise at 1 and 2 nA from a 10 K bath; depletion field across the disc; junction capacitance.
- Compare to: the paper's +15 K at 142 A/cm2 and +49 K at 283 A/cm2; the model's +16.5 and +46.1 K; V_j 3.09 V. Pass: rises within 20 percent without any fitted thermal resistance, V_j within 0.2 V. If the drift-diffusion V_j at 10 K is not within 0.1 V of E_g/q the low-temperature treatment in the transport module is wrong.

### N6. Resonant-tunnelling injector stack (Semiconductor, Schrodinger-Poisson, 1-D)

- Replaces: `nitride_nanowire_injector.py`: a 1-D transfer-matrix tunnelling model with the polarization sheet charges solved by a zero-net-drop closure (1.4755 MV/cm per 2 nm Al0.30 barrier, 0.295 eV, well tilted -0.59 eV), resonances at 268.9 meV (electron) and 48.7 / 101.2 / 166.0 meV (hole), misalignment 233 meV (electron) and 192 to 217 meV (hole) from the emitter Fermi level, thermionic bypass 0.0021 / 0.0151 at 230 / 300 K, a T = 0 degenerate reservoir formula [A].
- Model: GaN / Al0.30Ga0.70N 2 nm / GaN 4 nm / AlGaN 2 nm / GaN, Ga-polar, n-GaN emitter (3e18) and p-GaN collector (5e17, incomplete ionization), self-consistent Schrodinger-Poisson at 230 and 300 K, applied field 0, 50, 200 kV/cm; run both valence partitions (0.70 and 0.30 eV).
- Extract: barrier and well fields, resonance energies for electrons and holes relative to their band edges, their position relative to the emitter quasi-Fermi level, transmission versus energy, the thermionic fraction above the true barrier top.
- Compare to: the numbers above; pass when barrier fields are within 10 percent and resonances within 20 meV. If a resonance sits within 3 kT of the emitter Fermi level under a realistic bias, the injector verdict changes and the round is re-run.

### N7. Surface recombination access weight (Semiconductor, transient carrier continuity, axisymmetric)

- Replaces: `nitride_nanowire_surface.py`: a uniform-cylinder rate 2S/R with S 1e3 cm/s [E] and an occupied-dot access weight 0.05 [DR, from the hard-wall overlap] as the headline, 1.0 [A] as the conservative partner; at 1.0 the exciton lifetime is capped at 0.625 ns at 12.5 nm regardless of other physics, against Deshpande's measured 0.71 to 1.1 ns.
- Model: carrier continuity in the disc and reservoirs with sidewall surface recombination velocity S 1e2, 1e3, 1e4 cm/s and the confined-state density from N2; decay of a single injected pair; also a delocalised reservoir carrier.
- Extract: effective nonradiative rate for the localized exciton and for the reservoir carrier versus core radius; the effective access weight (localized rate divided by the uniform 2S/R rate).
- Compare to: 0.05; pass within a factor 2. Above 0.3 the headline flux drops by up to 4x and the sweep is re-run with the new default.

## 6. Cross-cutting checks

### X1. Pseudomorphic polarization field in a planar InGaN well (Electrostatics, 1-D)

- Validates `nitride_materials.polarization_field` (linear interpolation, Bernardini constants) at x 0.15, 0.25, 0.40 for a 3 nm well in GaN; compare the field to the module's value at strain fraction 1 and 0; pass within 5 percent (this is arithmetic; a miss means a constant is mistranscribed).

### X2. Material dispersion tables (no solve)

- Confirm the Sellmeier and tabulated indices the electromagnetic items use match the module's functions at 446, 450, 543, 630, 750 nm to 1e-3; record the values in the note.

## 7. What COMSOL cannot settle (do not attempt)

- Photon-counting statistics, g2, the one-pair-per-cycle loading and the finite-pulse re-excitation (`pulse_counting.py`, `cw_g2.py`, `loading.py`, `integrator.py`): these are rate and master-equation problems on top of the field-solved inputs.
- Exciton-phonon dephasing and the polaron lineshape (`linewidth.py`, `qd_gf.py`): phonon coupling constants, not classical fields.
- Multiband k.p mixing and valence-band ordering: not available as a native interface; single-band Schrodinger results must be labelled as such.
- The CW Hanbury Brown-Twiss measurement of Deshpande 2013: needs a rate-equation model of continuous injection with the detector response.
- The 52 percent ensemble yield ratio 300 K over 10 K: a defect-kinetics question the project has no model for.

## 8. Delivery checklist and order

Run in this order; each item is independent and can be delivered alone. Tier 4 first (it is the active round), then tier 3, then tier 2.

1. N1, N2 (the strain endpoints and the disc states: decide the headline bound).
2. N3, N4 (the two extraction models: decide the collection numbers and the dipole prior).
3. N5, N6, N7 (operating point, injector, surface).
4. P1, P2, P5 (cavity, Stark and screening, charging energy).
5. E1, E2, E3, E4, E5.
6. P3, P4, X1, X2.

For every item deliver the folder of section 2 with the note's PASS / DEVIATION / BLOCKED line. A DEVIATION is a result, not a failure: report the number, the module's number and the ratio, and stop. Never adjust a model input to close a gap with the module or with a measurement. If a licence blocks an item, say which interface and stop that item.
