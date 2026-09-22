# COMSOL work order: proofs for the nitride nanowire tier

Written 2026-09-22 by the orchestrator for itself. Purpose: every module in the nanowire tier (`fsim_core/nitride_nanowire_*.py`) is an analytic or 1-D model with declared approximations. The sweep (`out/nitride_nanowire/full`, commit after 5f35479) rests on a handful of those approximations. This work order lists the COMSOL simulations that would replace each [A]/[E] with a numerical reference, in the order of their leverage on the headline, with the exact geometry, the quantity to extract, the number in the current model it is compared to, and the pass criterion. Nothing here changes the simulator; each result becomes a `[V]`-class anchor in `verify/data/nitride_nanowire_anchors.yaml` with `evidence_kind: numerical_reference` and a verifier check against the analytic module.

Reference geometry throughout (the Deshpande 2013 as-built wire unless stated): GaN wire, core radius 12.5 nm (also 15 nm), In0.25Ga0.75N disc 2 nm thick filling the core, 15 nm undoped GaN each side, n-GaN 300 nm (Si 3e18 cm-3), p-GaN 250 nm (Mg 5e17 cm-3), wire lying on 100 nm thermal SiO2 on (001) Si, Ti/Au (5/45 nm) end contacts about 600 nm apart, emission 446 nm (x 0.25) and 543 nm (x 0.40 relaxed). Vertical family: same disc (12.5 nm) inside a 60 to 120 nm GaN core, bottom mirror, tapered top, 1 MOhm contact. Materials: GaN/InN/AlN constants as in `fsim_core/nitride_materials.py` (Bernardini 1997 polarization, Rinke 2008 masses, Ioffe elastic constants).

## Priority 1. Strain relaxation and polarization field in the disc (Solid Mechanics + Electrostatics)

Why first: the whole round hinges on the relaxed-versus-unrelaxed strain endpoints. The model has no radius-dependent relaxation law; it reports both endpoints. A single 3-D strain solve gives the actual number.

- Physics: 3-D linear elasticity with initial strain (lattice mismatch InGaN/GaN, wurtzite anisotropic stiffness), free sidewall; then electrostatics with spontaneous plus piezoelectric polarization as a body/interface charge (P_sp discontinuity plus e31, e33 times strain), zero net drop across the undoped stack between the doped reservoirs.
- Geometry sweep: core radius {10, 12.5, 15, 20, 25, 40} nm; disc thickness {1.5, 2, 3, 4} nm; x_in {0.25, 0.40}; no shell and 3 nm AlGaN shell.
- Extract: the volume-averaged in-plane strain in the disc relative to the pseudomorphic value (relaxation fraction), the axial polarization field at the disc centre, and the field profile along the axis.
- Compare to: relaxed endpoint 74.5 kV/cm at x 0.25 (112.8 at x 0.40, spontaneous only) and unrelaxed 4038 kV/cm (6469 at x 0.40) in `nitride_nanowire_levels.py`. Deshpande's less than 2 meV blueshift with current implies a field near the relaxed endpoint.
- Pass criterion: the relaxation fraction and field versus radius, with which a single `strain_fraction(R, h)` law can be derived and tagged [DR]; if the 12.5 nm disc sits below 30 percent of the pseudomorphic strain, the relaxed headline stands.

## Priority 2. Disc-in-wire single-particle states and radiative rate (Semiconductor / Schrodinger, axisymmetric)

- Physics: 2-D axisymmetric effective-mass Schrodinger equation for electron and hole on the field profile from Priority 1, BenDaniel-Duke mass discontinuities, finite InGaN/GaN offsets (VBO 1.15 eV [V Tsai 2020]), hard sidewall for the full-core disc and the finite lateral barrier for the 12.5 nm disc inside a 60 to 120 nm core (vertical family).
- Extract: ground-state energies, exciton transition energy (single-band plus the analytic Coulomb envelope the module uses, applied identically), electron-hole overlap, escape depths to the GaN reservoirs, the first transverse excited state.
- Compare to: relaxed E_X 2.777 eV (446.4 nm) at x 0.25 / 25 K; overlap 0.94; transverse spacing 13.93 meV (e) and 3.31 meV (h); vertical 12.5 nm disc in a 100 nm core E_perp 7.98 meV (hard wall would give 9.05).
- Pass criterion: transition energy within 15 meV, overlap within 5 percent, spacings within 10 percent. A miss on the vertical disc points at the adiabatic-decoupling [E] in the disc branch.

## Priority 3. Optical extraction, horizontal wire on SiO2/Si (Wave Optics, 3-D)

Why: the horizontal collection efficiency 0.079 (isotropic) and the antenna rate factor 0.385 come from a point dipole over a layered substrate plus a quasi-static cylinder screening; the wire body itself (n 2.49, 25 nm diameter, 600 nm long) is not in the electromagnetic model.

- Physics: 3-D frequency-domain electromagnetic wave, electric point dipole at the wire centre, GaN cylinder 25 nm diameter and 600 nm long lying on 100 nm SiO2 on Si (complex index), Ti/Au contact pads at the ends, PML boundaries, 446 and 543 nm.
- Extract: total radiated power relative to the same dipole in bulk GaN and in vacuum (Purcell / LDOS factor) for the three dipole orientations; far-field power within NA 0.5 above the substrate; degree of linear polarization of the collected light for an isotropic dipole population and for the c-plane-only population.
- Compare to: antenna rate factor 0.385 (isotropic) and 0.077 (c-plane only); eta_collection 0.079 (450 nm), 0.107 (630 nm); DOLP +84 percent (isotropic) and -100 percent (c-plane only) against the measured +70 percent axial.
- Pass criterion: rate factors within a factor 1.5 and collection within a factor 1.3; the sign of the DOLP for the c-plane population is the decisive check of the orientation prior.

## Priority 4. Vertical photonic wire (Wave Optics, 2-D axisymmetric, mode analysis plus dipole)

- Physics: mode analysis of the GaN cylinder in air at 446 and 543 nm for radius 60 to 120 nm (effective index, confinement, group index of HE11, LP11 cutoff); then a dipole-in-wire frequency-domain solve with a bottom Au or DBR mirror, a linear taper to a 1.5 um top, PML above.
- Extract: beta (fraction of emission into HE11), effective index and group index, extraction into NA 0.5 and 0.9 with and without the taper, the radius where the second mode appears.
- Compare to: beta 0.81 to 0.86 at d/lambda 0.24, the Marcuse confinement surrogate, the multimode penalty s(V) [E from Bleuse 2011], the LP11 cutoff at V 2.405 (100 nm core above cutoff at 543 nm), the group-index ratio 1.02 [A].
- Pass criterion: beta within 0.1 of the module across the radius sweep; the module's headline_eligible boundary (single-mode) within 5 nm of the COMSOL cutoff radius.

## Priority 5. Electro-thermal operating point of the as-built diode (Semiconductor + Heat Transfer + Electric Currents)

Why: the junction voltage (3.09 V at 300 K, 3.51 V at 10 K), the 2.8e9 K/W thermal resistance [DR fit to Fig. 4] and the tau_SRH 0.01 ns refit [A] are the fitted electrical and thermal inputs the results review asked to be disclosed.

- Physics: 1-D or 2-D axisymmetric drift-diffusion with Fermi statistics and incomplete ionization (Mg 170 meV), SRH and radiative recombination, the InGaN disc as a heterojunction; coupled heat transfer in the lying wire with the SiO2/Si substrate and the Ti/Au pads as heat sinks, Joule heating from the series resistance path.
- Extract: I-V from 0.5 to 5 nA at 10 K and 300 K; junction voltage at 1 and 2 nA; junction temperature rise at 1 and 2 nA from a 10 K bath; the depletion field across the disc at the operating point.
- Compare to: Deshpande 2013 rises +15 K (142 A/cm2) and +49 K (283 A/cm2); the model's +16.5 / +46.1 K with Rth 2.8e9; V_j 3.09 V (300 K); the RC time 1.65 ns with C_dep 0.34 aF.
- Pass criterion: junction rise within 20 percent without a fitted Rth; V_j within 0.2 V; if the drift-diffusion V_j at 10 K is not within 0.1 V of E_g/q the transport module's low-temperature V_bi treatment is wrong.

## Priority 6. Resonant-tunnelling injector stack (Semiconductor, 1-D, plus Schrodinger)

- Physics: self-consistent Poisson-Schrodinger of the GaN / Al0.30Ga0.70N (2 nm) / GaN (4 nm) / AlGaN (2 nm) / GaN stack with spontaneous and piezoelectric sheet charges, Ga-polar, n-GaN emitter (3e18) and p-GaN collector; transmission versus energy at 0, 50, 200 kV/cm applied field.
- Extract: the built-in field in each barrier and the well (zero-net-drop closure gives 1.4755 MV/cm per barrier, 0.295 eV, well minus 0.59 eV), the resonance energies (electron 268.9 meV, hole 48.7 / 101.2 / 166.0 meV above the respective band edges), the resonance versus the emitter Fermi level (misalignment 233 meV electron, 192 to 217 meV hole), thermionic bypass fraction at 230 and 300 K (0.0021 / 0.0151).
- Pass criterion: barrier fields within 10 percent, resonances within 20 meV; if the self-consistent solution shows a resonance within 3 kT of the emitter Fermi level under a realistic bias, the injector verdict changes.

## Priority 7. Surface recombination in the disc (Semiconductor, 2-D axisymmetric, transient)

- Physics: carrier continuity in the disc and reservoirs with a sidewall surface recombination velocity S {1e2, 1e3, 1e4} cm/s and the confined-state density from Priority 2; decay of a single injected pair.
- Extract: effective nonradiative rate for a localized exciton versus a delocalized reservoir carrier as a function of core radius; the effective access weight (ratio of the confined-state loss to the uniform-density 2S/R value).
- Compare to: occupied_dot_access 0.05 [DR from the hard-wall overlap] against the conservative 1.0; the lifetime cap 0.625 ns at access 1.0; Deshpande's TRPL 0.71 ns and g2-fit 1.1 ns read as no nonradiative loss at 25 K.
- Pass criterion: the extracted access weight within a factor 2 of 0.05; if it is above 0.3 the headline flux figures drop by up to 4x and the sweep must be re-run with the new default.

## Deliverables per simulation

For each item: the COMSOL model file, a one-page note with the mesh convergence statement (result change under one refinement below 2 percent), the extracted numbers in a CSV named by the item, and the comparison table against the module values above. Each CSV row becomes a ledger anchor with `evidence_kind: numerical_reference`, tag `[V]` only for quantities that reproduce a published measurement (Priority 5 heating; Priority 3 DOLP), `[DR]` otherwise. The analytic module keeps its verifier; a new `verify/verify_nitride_nanowire_comsol_anchors.py` compares module outputs to the anchors with the pass criteria above and prints N/N.

## What this does not cover

The one-pair-per-cycle loading itself (the Coulomb and resonant-tunnelling screens) is a transport-statistics question, not a field solve; COMSOL confirms the charging energy (Electrostatics: capacitance of the disc in the wire, model value 12.13 meV at 12.5 nm) but not the loading statistics. The CW g2 the 2013 paper measured needs a rate-equation or master-equation model, not COMSOL.
