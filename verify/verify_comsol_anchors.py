"""COMSOL anchor intake harness (audit Phase C item 3).

Walks ``comsol/<tier>/<item-id>/`` for the 19 items of the COMSOL
specification and validates every delivered folder:

* file set (model.java, params.csv, results.csv, convergence.csv, run.log,
  note.md; X2 is a no-solve item and needs only results.csv and note.md);
* results.csv schema: every column of the item's schema present, with the
  unit written in the header as ``name [unit]``, plus ``mesh_level`` and
  ``evidence_tag`` ([DR] field solve, or [V] when the item reproduces a
  published measurement);
* convergence.csv: every results row present at mesh_level ``production`` and
  ``refined``, the refined mesh has more elements, and every numeric quantity
  changes by less than 2 percent (energies: less than 1 meV), "Mesh
  convergence";
* note.md provenance: solver version, geometry, boundary conditions, mesh,
  licence modules, run time, operator, date, constants with evidence tags,
  a comparison table and the one-line PASS / DEVIATION / BLOCKED verdict;
* cheap cross-checks against the surrogate numbers (fsim_core calls or closed
  forms).  A mismatch is a DEVIATION, which is noted but does not fail the
  harness; it fails only when note.md claims PASS for a row the harness finds
  outside the item's pass criterion (an inconsistent claim).

A folder holding only README.md is PENDING and passes.  A folder whose
note.md says BLOCKED (licence) and that has no results is BLOCKED and passes.

Usage:
    python verify/verify_comsol_anchors.py              # walk comsol/
    python verify/verify_comsol_anchors.py --self-test  # fixture self-test
    python verify/verify_comsol_anchors.py --write-readmes  # regenerate READMEs

Prints ``N/N ... passed`` and exits 0 iff every check passes.
"""
from __future__ import annotations

import csv
import math
import re
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Callable

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
COMSOL_DIR = ROOT / "comsol"
WORK_ORDER = "docs/comsol_work_order_nanowire.md"

# CODATA 2018 exact / recommended SI constants [V] (Tiesinga et al., Rev. Mod.
# Phys. 93, 025010 (2021)).
E_CHARGE = 1.602176634e-19      # C, exact
EPS0 = 8.8541878128e-12         # F/m
KB = 1.380649e-23               # J/K, exact

# Convergence rule [A: project acceptance rule].
CONV_REL = 0.02                 # 2 percent for non-energy quantities
CONV_ENERGY_MEV = 1.0           # 1 meV for energies
ENERGY_UNITS_MEV = {"meV": 1.0, "eV": 1000.0}
MESH_LEVELS = ("production", "refined")
EVIDENCE_TAGS = {"[DR]", "[V]"}
NA_TOKEN = "NA"                 # not-applicable cell (e.g. 1-D row of a 2-D-only quantity)

NOTE_FIELDS = ("Solver version", "Geometry", "Boundary conditions", "Mesh",
               "Licence modules", "Run time", "Operator", "Date", "Constants",
               "Verdict")
NOTE_FIELDS_BLOCKED = ("Operator", "Date", "Verdict")


# ---------------------------------------------------------------- schema


@dataclass(frozen=True)
class Col:
    name: str
    unit: str           # "label" for categorical parameters
    desc: str

    @property
    def header(self) -> str:
        return f"{self.name} [{self.unit}]"

    @property
    def is_label(self) -> bool:
        return self.unit == "label"


@dataclass
class Ref:
    """Cheap cross-check of one quantity on the rows matching ``where``.

    mode: "abs" |v-r| <= tol; "rel" |v-r| <= tol*|r|; "factor" max(v/r, r/v)
    <= tol; "mag_rel" as rel on magnitudes (sign conventions differ).
    ``where`` values may be a tuple (any of).  A ``diagnostic`` ref tests the
    module, not COMSOL: a miss is reported but never counts against a PASS
    verdict."""
    quantity: str
    where: dict
    ref: float | Callable[[dict], float]
    mode: str
    tol: float
    source: str
    diagnostic: bool = False


@dataclass
class Item:
    id: str
    tier: str
    section: str
    title: str
    physics: str
    replaces: str
    params: list
    quantities: list
    compare_to: str
    pass_rule: str
    refs: list = field(default_factory=list)
    solve: bool = True
    extra: Callable | None = None

    @property
    def columns(self):
        return self.params + self.quantities

    @property
    def anchor(self) -> str:
        return f"verify/data/{self.tier}_anchors.yaml, anchor id `{self.id}` (written by the orchestrator, work order section 2)"


SEC3 = "3. Tier 2, room-temperature edge-emitting InP dot"
SEC4 = "4. Tier 3, planar InGaN/GaN cavity dot"
SEC5 = "5. Tier 4, nitride nanowire"
SEC6 = "6. Cross-cutting checks"

C = Col


# -- surrogate reference callables (lazy fsim_core imports) -----------------

def sphere_EC_meV(radius_nm: float, eps_r: float) -> float:
    """Closed form isolated-sphere charging energy E_C = e^2 / (4 pi eps0
    eps_r R), the convention of drive_mech.set_feasibility ([V] Pekola et al.,
    Rev. Mod. Phys. 85, 1421 (2013), Sec. II.B)."""
    C_F = 4.0 * math.pi * EPS0 * eps_r * radius_nm * 1e-9
    return E_CHARGE / C_F * 1e3


def _module_EC(row):
    from fsim_core import drive_mech
    # eps_r 9.5: the value the nanowire device path passes for the 12.5 nm
    # disc [A, per verify_nitride_nanowire_device Coulomb anchor 12.126 meV].
    return drive_mech.set_feasibility(300.0, radius_nm=row["island_radius"], eps_r=9.5)["E_C_meV"]


def _module_pol_field(row):
    from fsim_core import nitride_materials as nm
    return nm.polarization_field(nm.ingaN(row["x_in"]), nm.binary("GaN"), 300.0,
                                 strain_fraction=row["strain_fraction"])


def _module_n_gan(row):
    from fsim_core import nitride_nanowire_photonics as ph
    return ph.gan_ordinary_index(row["lambda"])


def _module_n_sio2(row):
    from fsim_core import nitride_nanowire_photonics as ph
    return ph.sio2_index(row["lambda"])


def _module_n_si(row):
    from fsim_core import nitride_nanowire_photonics as ph
    return ph.si_complex_index(row["lambda"]).real


def _module_k_si(row):
    from fsim_core import nitride_nanowire_photonics as ph
    return ph.si_complex_index(row["lambda"]).imag


def _lp11_cutoff_radius_nm(row):
    """Module single-mode boundary: V = 2 pi R / lambda sqrt(n^2 - 1) = 2.405
    (Snyder and Love 1983 [V]) with the module's Barker-Ilegems GaN index."""
    from fsim_core import nitride_nanowire_photonics as ph
    n = ph.gan_ordinary_index(row["lambda"])
    return 2.405 * row["lambda"] / (2.0 * math.pi * math.sqrt(n * n - 1.0))


def _gan_gap_eV(row):
    """E_g(T_bath)/q of GaN from nitride_materials.bandgap (the module's
    Varshni form; its tags apply). A drift-diffusion diode's V_j approaches it
    at 10 K."""
    from fsim_core import nitride_materials as nm
    return nm.bandgap(nm.binary("GaN"), row["T_bath"])


def _uniform_surface_rate_per_ns(row):
    """Uniform-cylinder surface rate 2S/R of nitride_nanowire_surface.py
    [DR, closed form]; S in cm/s, R in nm, result in 1/ns."""
    return 2.0 * row["S"] * 1e-2 / (row["core_radius"] * 1e-9) * 1e-9


def _tau_note_check(item, rows, note_text):
    """E3/P3: if the 1/e thermal time constant is below the 12.5 ns period
    (80 MHz), the note must say the duty-averaged treatment is wrong."""
    probs = []
    taus = [r["tau_th"] for r in rows if isinstance(r.get("tau_th"), float)]
    if taus and min(taus) < 12.5 and "duty" not in note_text.lower():
        probs.append("tau_th < 12.5 ns but note.md does not discuss the duty-averaged treatment (work order E3)")
    return probs


ITEMS: list[Item] = [
    # ---------------------------------------------------------------- tier 2
    Item("E1", "rt_edge", SEC3, "Ridge waveguide full-vector mode",
         "Wave Optics, mode analysis, 2-D cross-section",
         "fsim_core/waveguide.py scalar 1-D slab effective-index solve with a Gaussian far-field",
         [C("ridge_width", "um", "ridge width"), C("etch_depth", "um", "etch depth"),
          C("lambda", "nm", "wavelength (815.7 and 750)"), C("mode", "label", "TE or TM fundamental")],
         [C("n_eff", "1", "effective index"), C("n_group", "1", "group index from two wavelengths 5 nm apart"),
          C("A_eff", "um^2", "effective mode area at the dot position"),
          C("Gamma_well", "1", "fraction of mode energy in the well layer"),
          C("beta", "1", "beta factor from the normalised field at the dot centre"),
          C("theta_par", "deg", "far-field half-angle, plane parallel to the junction"),
          C("theta_perp", "deg", "far-field half-angle, plane perpendicular to the junction")],
         "waveguide.slab_modes / effective_index_ridge / beta_factor / na_collection via waveguide.edge_emission(hkust_ridge_stack(815.7), ...) with the card's values",
         "n_eff within 0.01, Gamma within 10 percent, far-field half-angles within 15 percent, beta within a factor 1.3"),
    Item("E2", "rt_edge", SEC3, "Facet transmission and NA collection",
         "Wave Optics, frequency domain, 2-D or 3-D",
         "waveguide.facet_transmission, na_collection and facet_escape_fraction",
         [C("lambda", "nm", "wavelength"), C("coating", "label", "none or the card's coating")],
         [C("T_facet", "1", "facet power transmission"),
          C("eta_NA_0p5", "1", "far-field fraction within NA 0.5"),
          C("eta_NA_0p7", "1", "far-field fraction within NA 0.7")],
         "waveguide facet transmission and eta_NA for the same NA",
         "within 15 percent"),
    Item("E3", "rt_edge", SEC3, "Junction self-heating of the mesa",
         "Heat Transfer in Solids, 3-D or axisymmetric, stationary and transient",
         "fsim_core/thermal.py analytic spreading-resistance network",
         [C("mesa_diameter", "um", "mesa diameter"), C("T_hs", "K", "heat-sink temperature (230, 300)"),
          C("power_mode", "label", "duty_avg or on_state")],
         [C("dT_junction", "K", "junction temperature rise"), C("R_th", "K/W", "thermal resistance"),
          C("tau_th", "ns", "1/e thermal time constant (transient study)")],
         "thermal.py T_j - T_hs for the same inputs",
         "within 20 percent; if tau_th < 12.5 ns the note must say the duty-averaged treatment is wrong",
         extra=_tau_note_check),
    Item("E4", "rt_edge", SEC3, "p-i-n injection and junction field",
         "Semiconductor (semi), drift-diffusion with Fermi statistics, 1-D and 2-D axisymmetric",
         "fsim_core/transport.py abrupt-junction electrostatics and SRH / radiative / diffusion channels",
         [C("T", "K", "temperature (230, 300)"), C("geometry", "label", "1D or 2D_axi"),
          C("I", "mA", "terminal current of the operating point")],
         [C("V_j", "V", "junction voltage"), C("F_well", "kV/cm", "depletion field at the well"),
          C("eta_inj", "1", "fraction of current recombining in the well"),
          C("f_aperture", "1", "fraction of current inside the aperture radius (NA for 1D rows)")],
         "transport.py V_j, field and eta_inj at the same point (full I-V 1 nA to 1 mA may be added as iv.csv)",
         "V_j within 0.1 V, field within 20 percent, eta_inj within 20 percent; f_aperture is new information"),
    Item("E5", "rt_edge", SEC3, "Dot confinement with a real strain field",
         "Solid Mechanics then Schrodinger Equation (single band, no k.p), 3-D or axisymmetric",
         "fsim_core/dot_levels.py separable disc model",
         [C("shape", "label", "lens or disc (4 nm by 12 nm radius)")],
         [C("E_e", "meV", "electron ground state above the strained CB edge"),
          C("E_h", "meV", "hole ground state below the strained VB edge"),
          C("E_trans", "meV", "transition energy"),
          C("escape_e", "meV", "electron escape depth to the well continuum"),
          C("escape_h", "meV", "hole escape depth to the well continuum"),
          C("dE_inplane_e", "meV", "electron in-plane excited-state spacing"),
          C("dE_inplane_h", "meV", "hole in-plane excited-state spacing"),
          C("eps_hydro", "1", "hydrostatic strain at the dot centre"),
          C("eps_biax", "1", "biaxial strain at the dot centre")],
         "dot_levels.py levels for the card's dot",
         "transition energy within 50 meV, escape depths within 30 percent (a larger miss is a result)"),
    # ---------------------------------------------------------------- tier 3
    Item("P1", "nitride_cavity", SEC4, "Dielectric DBR microcavity: Q, mode volume, Purcell, out-coupling",
         "Wave Optics, eigenfrequency then dipole excitation, 2-D axisymmetric",
         "fsim_core/dbr.py 1-D transfer matrix and fsim_core/nitride_cavity.py assumed Purcell inputs",
         [C("mesa_diameter", "um", "mesa diameter (20 and 2)")],
         [C("lambda_res", "nm", "resonant wavelength"), C("Q", "1", "quality factor"),
          C("V_mode", "(lambda/n)^3", "mode volume"), C("F_P", "1", "Purcell factor at the dot position"),
          C("eta_out_NA_0p5", "1", "fraction of emitted power through the top DBR into NA 0.5"),
          C("dlambda_dT", "nm/K", "resonance shift per kelvin with the module's dn/dT")],
         "the card's Q (2000 [A]), mode_volume_norm, spatial_overlap, eta_out and dEdT_cav",
         "Q within a factor 2 (or the note states the achievable Q), mode volume within 30 percent, eta_out within 20 percent"),
    Item("P2", "nitride_cavity", SEC4, "Polarization field, screening and Stark shift in the planar disc",
         "Semiconductor, Schrodinger-Poisson, 1-D and axisymmetric",
         "fsim_core/nitride_levels.py (screening_fraction knob) and fsim_core/nitride_stark.py",
         [C("h_disc", "nm", "disc height (1, 3, 5, 10)"), C("x_in", "1", "indium fraction (0.25, 0.40)"),
          C("bias", "V", "applied bias (-2 to +3)"), C("screening", "1", "free-carrier sheet / polarization charge (0, 0.5, 1)"),
          C("orientation", "label", "c_plane or nonpolar")],
         [C("F_disc", "kV/cm", "field in the disc"), C("E_X", "meV", "electron-hole transition energy"),
          C("overlap", "1", "electron-hole overlap"), C("escape_e", "meV", "electron escape depth"),
          C("escape_h", "meV", "hole escape depth"),
          C("stark_slope", "meV/V", "dE_X/dV of the (h, x, screening, orientation) series, repeated on each row"),
          C("V_flatband", "V", "bias at which the built-in field is compensated, repeated on each row")],
         "nitride_stark.py pinned slopes (-12.95, -9.55, +2.88 meV/V at 3 nm, screening 0 / 0.5 / 1) and nitride_levels.py E_X, overlap, escape depths",
         "Stark slope within 2 meV/V, E_X within 20 meV",
         refs=[Ref("stark_slope", {"h_disc": 3.0, "x_in": 0.25, "screening": s, "orientation": "c_plane"},
                   v, "abs", 2.0, "nitride_stark.py pinned slope, work order section 4 header")
               for s, v in ((0.0, -12.95), (0.5, -9.55), (1.0, 2.88))]),
    Item("P3", "nitride_cavity", SEC4, "Mesa self-heating",
         "Heat Transfer in Solids, axisymmetric (as E3)",
         "fsim_core/thermal.py at the planar operating point",
         [C("mesa_diameter", "um", "mesa diameter (20)"), C("T_hs", "K", "heat-sink temperature (230, 300)"),
          C("power_mode", "label", "duty_avg or on_state")],
         [C("dT_junction", "K", "junction temperature rise"), C("R_th", "K/W", "thermal resistance"),
          C("tau_th", "ns", "1/e thermal time constant")],
         "thermal.py T_j - T_hs at the planar operating point (GaN on Si, DBR as thermal barrier)",
         "within 20 percent", extra=_tau_note_check),
    Item("P4", "nitride_cavity", SEC4, "Planar p-i-n injection and aperture partition",
         "Semiconductor, axisymmetric, Mg incomplete ionization 170 meV",
         "fsim_core/nitride_transport.py eta_inj, depletion field and aperture supply fraction",
         [C("T", "K", "temperature"), C("geometry", "label", "1D or 2D_axi"), C("I", "mA", "terminal current")],
         [C("V_j", "V", "junction voltage"), C("F_depl", "kV/cm", "depletion field at the disc"),
          C("eta_inj", "1", "injection efficiency"),
          C("f_aperture", "1", "aperture supply fraction (NA for 1D rows)")],
         "nitride_transport.py eta_inj, depletion field and aperture supply fraction",
         "within 20 percent"),
    Item("P5", "nitride_cavity", SEC4, "Charging energy of the injection island",
         "AC/DC Electrostatics, 3-D, island held at 1 V",
         "drive_mech.set_feasibility isolated-sphere capacitance E_C = e^2/C, C = 4 pi eps0 eps_r R",
         [C("island_radius", "nm", "island radius (0.5, 1, 5, 12.5)"),
          C("environment", "label", "planar or nanowire epitaxial environment")],
         [C("C_total", "aF", "total self-capacitance"), C("E_C", "meV", "charging energy e^2/C"),
          C("EC_over_kT_230K", "1", "E_C / kT at 230 K"), C("EC_over_kT_300K", "1", "E_C / kT at 300 K")],
         "drive_mech.set_feasibility E_C_meV (12.126 meV at 12.5 nm, eps_r 9.5)",
         "E_C within 30 percent (a larger capacitance strengthens the infeasibility verdict)",
         refs=[Ref("E_C", {}, _module_EC, "rel", 0.30, "drive_mech.set_feasibility(radius_nm=R, eps_r=9.5)")]),
    # ---------------------------------------------------------------- tier 4
    Item("N1", "nitride_nanowire", SEC5, "Strain relaxation and polarization field in the disc",
         "Solid Mechanics then Electrostatics, 3-D or axisymmetric",
         "relaxed / pseudomorphic endpoints of nitride_nanowire_levels.py (74.5 vs 4038 kV/cm at x 0.25; 112.8 vs 6469 at x 0.40)",
         [C("core_radius", "nm", "core radius (10, 12.5, 15, 20, 25, 40)"),
          C("h_disc", "nm", "disc thickness (1.5, 2, 3, 4)"), C("x_in", "1", "indium fraction (0.25, 0.40)"),
          C("shell", "label", "none or AlGaN_3nm")],
         [C("strain_fraction", "1", "volume-averaged in-plane strain / pseudomorphic value"),
          C("F_axial_center", "kV/cm", "axial field at the disc centre (module sign convention, +c positive)")],
         "the two endpoints; the result derives a strain_fraction(R, h) law tagged [DR] (axial profiles may be added as profile.csv)",
         "under 30 percent retained at 12.5 nm: relaxed headline stands; over 70 percent: conservative bound becomes the headline"),
    Item("N2", "nitride_nanowire", SEC5, "Disc-in-wire single-particle states",
         "Schrodinger Equation, axisymmetric, on the N1 field",
         "nitride_nanowire_levels.py hard-wall / finite-barrier disc with adiabatic decoupling",
         [C("x_in", "1", "indium fraction"), C("family", "label", "horizontal or vertical"),
          C("core_radius", "nm", "GaN core radius"), C("disc_radius", "nm", "disc radius"),
          C("T", "K", "temperature")],
         [C("E_e", "meV", "electron ground state"), C("E_h", "meV", "hole ground state"),
          C("E_trans", "meV", "transition energy with the module's analytic Coulomb envelope"),
          C("overlap", "1", "electron-hole overlap"),
          C("escape_e", "meV", "electron escape depth to the GaN reservoir"),
          C("escape_h", "meV", "hole escape depth to the GaN reservoir"),
          C("dE_perp_e", "meV", "first transverse excited-state spacing, electron; on family=vertical rows this is "
            "the finite-lateral-barrier E_perp of the disc inside the wider core (work order: 7.98 meV for the "
            "12.5 nm disc in a 100 nm core; a hard wall would give 9.05)"),
          C("dE_perp_h", "meV", "first transverse excited-state spacing, hole")],
         "relaxed E_X 2.777 eV at x 0.25, 25 K; overlap 0.94; transverse spacings 13.93 / 3.31 meV; vertical 12.5 nm disc in a 100 nm core E_perp 7.98 meV (dE_perp_e on family=vertical rows)",
         "transition energy within 15 meV, overlap within 5 percent, spacings within 10 percent",
         refs=[Ref(q, {"x_in": 0.25, "family": "horizontal", "core_radius": 12.5, "disc_radius": 12.5, "T": 25.0},
                   v, m, t, "nitride_nanowire_levels.py relaxed bound, work order N2")
               for q, v, m, t in (("E_trans", 2777.0, "abs", 15.0), ("overlap", 0.94, "rel", 0.05),
                                  ("dE_perp_e", 13.93, "rel", 0.10), ("dE_perp_h", 3.31, "rel", 0.10))]
         # Vertical family: 12.5 nm disc in a core of radius 100 nm, the module's
         # own case in verify_nitride_nanowire_levels (H6/L1: core_radius_nm=100,
         # disc_radius_nm=12.5, x 0.25, relaxed) [DR].
         + [Ref("dE_perp_e", {"x_in": 0.25, "family": "vertical", "core_radius": 100.0, "disc_radius": 12.5},
                7.98, "rel", 0.10, "nitride_nanowire_levels finite-barrier disc-in-wire E_perp, work order N2")]),
    Item("N3", "nitride_nanowire", SEC5, "Horizontal wire on SiO2/Si: radiative rate, collection, polarization",
         "Wave Optics, frequency domain, 3-D",
         "nitride_nanowire_photonics.py horizontal family (point dipole over air/SiO2/Si, quasi-static cylinder screening)",
         [C("lambda", "nm", "wavelength (446, 543)"),
          C("dipole", "label", "axial, transverse, vertical, isotropic or cplane_only")],
         [C("rate_factor_bulk", "1", "radiated power relative to the same dipole in bulk GaN"),
          C("rate_factor_vacuum", "1", "radiated power relative to the same dipole in vacuum"),
          C("eta_NA_0p5", "1", "far-field power within NA 0.5 above the substrate"),
          C("DOLP", "1", "degree of linear polarization of collected light, + along the wire axis")],
         "collection 0.079 (450 nm), antenna rate factor 0.385 (isotropic) / 0.077 (c-plane-only), DOLP +84 / -100 percent vs measured +70",
         "rate factors within a factor 1.5, collection within a factor 1.3; the c-plane-only DOLP sign is decisive",
         refs=[Ref("eta_NA_0p5", {"lambda": 446.0, "dipole": "isotropic"}, 0.079, "factor", 1.3,
                   "nitride_nanowire_photonics horizontal collection 0.079 at 450 nm")]),
    Item("N4", "nitride_nanowire", SEC5, "Vertical photonic wire: modes, beta, extraction",
         "Wave Optics, mode analysis then frequency domain, axisymmetric",
         "vertical family of nitride_nanowire_photonics.py (Marcuse Gaussian HE11, beta estimator, s(V), LP11 cutoff V 2.405)",
         [C("core_radius", "nm", "GaN cylinder radius (60 to 120)"), C("lambda", "nm", "wavelength (446, 543)"),
          C("mirror", "label", "Au or DBR"), C("taper", "label", "none or linear_1p5um")],
         [C("n_eff_HE11", "1", "HE11 effective index"), C("n_group_HE11", "1", "HE11 group index"),
          C("beta", "1", "fraction of emission into HE11, both directions"),
          C("R_cutoff_2nd", "nm", "radius at which the second guided mode appears"),
          C("eta_NA_0p5", "1", "extraction into NA 0.5"), C("eta_NA_0p9", "1", "extraction into NA 0.9")],
         "beta 0.81 to 0.86 near d/lambda 0.24; module single-mode boundary V = 2.405",
         "beta within 0.1 across the radius sweep, cutoff radius within 5 nm",
         refs=[Ref("R_cutoff_2nd", {}, _lp11_cutoff_radius_nm, "abs", 5.0,
                   "module LP11 cutoff V = 2.405 with gan_ordinary_index")]),
    Item("N5", "nitride_nanowire", SEC5, "Electro-thermal operating point of the as-built diode",
         "Semiconductor, Heat Transfer, Electric Currents; 2-D axisymmetric wire plus 3-D thermal model",
         "nitride_nanowire_transport.py (V_j 3.09 V at 300 K, fitted R_th 2.8e9 K/W, RC diagnostic)",
         [C("T_bath", "K", "bath temperature (10, 300)"), C("I", "nA", "current (0.5 to 5)")],
         [C("V_j", "V", "junction voltage"), C("dT_junction", "K", "junction rise from the bath"),
          C("F_disc", "kV/cm", "depletion field across the disc"), C("C_j", "aF", "junction capacitance")],
         "Deshpande 2013 +15 K at 142 A/cm2 and +49 K at 283 A/cm2 (1 and 2 nA); model +16.5 / +46.1 K; V_j 3.09 V",
         "rises within 20 percent without a fitted thermal resistance, V_j within 0.2 V",
         refs=[Ref("V_j", {"T_bath": 300.0, "I": (1.0, 2.0)}, 3.09, "abs", 0.2,
                   "nitride_nanowire_transport V_j(300 K) at 1 and 2 nA"),
               Ref("V_j", {"T_bath": 10.0, "I": (1.0, 2.0)}, _gan_gap_eV, "abs", 0.1,
                   "E_g(10 K)/q of GaN, nitride_materials.bandgap; a miss means the transport "
                   "module's low-temperature treatment is wrong (work order N5)", diagnostic=True),
               Ref("dT_junction", {"T_bath": 10.0, "I": 1.0}, 15.0, "rel", 0.20,
                   "Deshpande et al., Nat. Commun. 4, 1675 (2013), Fig. 4 [V]"),
               Ref("dT_junction", {"T_bath": 10.0, "I": 2.0}, 49.0, "rel", 0.20,
                   "Deshpande et al., Nat. Commun. 4, 1675 (2013), Fig. 4 [V]")]),
    Item("N6", "nitride_nanowire", SEC5, "Resonant-tunnelling injector stack",
         "Semiconductor, Schrodinger-Poisson, 1-D",
         "nitride_nanowire_injector.py transfer-matrix tunnelling with the zero-net-drop closure",
         [C("T", "K", "temperature (230, 300)"), C("F_applied", "kV/cm", "applied field (0, 50, 200)"),
          C("vbo_GaN_AlN", "eV", "GaN/AlN valence-band offset partition (0.70, 0.30)")],
         [C("F_barrier", "kV/cm", "field in the AlGaN barrier"), C("F_well", "kV/cm", "field in the GaN well"),
          C("E_res_e1", "meV", "first electron resonance above the CB edge"),
          C("E_res_h1", "meV", "first hole resonance"), C("E_res_h2", "meV", "second hole resonance"),
          C("E_res_h3", "meV", "third hole resonance"),
          C("dE_res_e_mu", "meV", "electron resonance minus emitter quasi-Fermi level"),
          C("dE_res_h_mu", "meV", "hole resonance minus emitter quasi-Fermi level"),
          C("f_thermionic", "1", "thermionic fraction above the true barrier top")],
         "barrier field 1.4755 MV/cm; resonances 268.9 (e), 48.7 / 101.2 / 166.0 meV (h); misalignment 233 / 192 to 217 meV; bypass 0.0021 / 0.0151",
         "barrier fields within 10 percent, resonances within 20 meV",
         refs=[Ref("F_barrier", {"F_applied": 0.0, "vbo_GaN_AlN": 0.70}, 1475.5, "mag_rel", 0.10,
                   "nitride_nanowire_injector zero-net-drop closure")]
         + [Ref(q, {"F_applied": 0.0, "vbo_GaN_AlN": 0.70}, v, "abs", 20.0, "nitride_nanowire_injector resonances")
            for q, v in (("E_res_e1", 268.9), ("E_res_h1", 48.7), ("E_res_h2", 101.2), ("E_res_h3", 166.0))]),
    Item("N7", "nitride_nanowire", SEC5, "Surface recombination access weight",
         "Semiconductor, transient carrier continuity, axisymmetric",
         "nitride_nanowire_surface.py uniform-cylinder rate 2S/R with access weight 0.05 [DR]",
         [C("core_radius", "nm", "core radius"), C("S", "cm/s", "surface recombination velocity (1e2, 1e3, 1e4)"),
          C("carrier", "label", "localized or reservoir")],
         [C("rate_nonrad", "1/ns", "effective nonradiative rate"),
          C("access_weight", "1", "localized rate / uniform 2S/R rate")],
         "access weight 0.05 (headline) / 1.0 (conservative); uniform rate 2S/R",
         "access weight within a factor 2 of 0.05; above 0.3 the sweep is re-run",
         refs=[Ref("access_weight", {"carrier": "localized", "core_radius": 12.5, "S": 1000.0}, 0.05, "factor", 2.0,
                   "nitride_nanowire_surface occupied-dot access weight"),
               Ref("rate_nonrad", {"carrier": "reservoir"}, _uniform_surface_rate_per_ns, "factor", 2.0,
                   "closed form 2S/R")]),
    # ---------------------------------------------------------------- cross
    Item("X1", "cross_cutting", SEC6, "Pseudomorphic polarization field in a planar InGaN well",
         "AC/DC Electrostatics, 1-D, 3 nm well in GaN",
         "validates nitride_materials.polarization_field (Bernardini constants, linear interpolation)",
         [C("x_in", "1", "indium fraction (0.15, 0.25, 0.40)"), C("strain_fraction", "1", "1 (pseudomorphic) or 0")],
         [C("F_well", "kV/cm", "field in the well (magnitude compared; state the sign convention)")],
         "nitride_materials.polarization_field(ingaN(x), binary('GaN'), 300, strain_fraction=...)",
         "within 5 percent",
         refs=[Ref("F_well", {}, _module_pol_field, "mag_rel", 0.05, "nitride_materials.polarization_field")]),
    Item("X2", "cross_cutting", SEC6, "Material dispersion tables",
         "no solve", "confirms the indices the electromagnetic items use",
         [C("lambda", "nm", "wavelength (446, 450, 543, 630, 750)")],
         [C("n_GaN", "1", "GaN ordinary index"), C("n_SiO2", "1", "SiO2 index"),
          C("n_Si", "1", "Si real index"), C("k_Si", "1", "Si extinction coefficient")],
         "nitride_nanowire_photonics.gan_ordinary_index / sio2_index / si_complex_index",
         "match to 1e-3",
         refs=[Ref("n_GaN", {}, _module_n_gan, "abs", 1e-3, "gan_ordinary_index"),
               Ref("n_SiO2", {}, _module_n_sio2, "abs", 1e-3, "sio2_index"),
               Ref("n_Si", {}, _module_n_si, "abs", 1e-3, "si_complex_index real"),
               Ref("k_Si", {}, _module_k_si, "abs", 1e-3, "si_complex_index imag")],
         solve=False),
]
ITEM_BY_ID = {it.id: it for it in ITEMS}


# ---------------------------------------------------------------- parsing

_HDR = re.compile(r"^\s*([A-Za-z0-9_]+)\s*(?:\[([^\]]*)\])?\s*$")


def _parse_header(cells):
    out = {}
    for i, c in enumerate(cells):
        m = _HDR.match(c)
        if m:
            out[m.group(1)] = (i, (m.group(2) or "").strip())
        else:
            out[c.strip()] = (i, None)
    return out


def _read_csv(path: Path):
    with path.open(newline="", encoding="utf-8") as fh:
        rows = [r for r in csv.reader(fh) if any(c.strip() for c in r)]
    if not rows:
        return None, []
    return rows[0], rows[1:]


def _check_table(item: Item, path: Path, extra_cols, probs):
    """Parse a results/convergence table into dicts; append problems."""
    label = path.name
    header, body = _read_csv(path)
    if header is None:
        probs.append(f"{label}: empty")
        return []
    hmap = _parse_header(header)
    ok = True
    for col in item.columns:
        if col.name not in hmap:
            probs.append(f"{label}: missing column '{col.header}'")
            ok = False
        elif hmap[col.name][1] != col.unit:
            probs.append(f"{label}: wrong units for '{col.name}': [{hmap[col.name][1]}] != [{col.unit}]")
            ok = False
    for name, unit in extra_cols:
        if name not in hmap:
            probs.append(f"{label}: missing column '{name}{' [' + unit + ']' if unit else ''}'")
            ok = False
        elif unit and hmap[name][1] != unit:
            probs.append(f"{label}: wrong units for '{name}': [{hmap[name][1]}] != [{unit}]")
            ok = False
    if not ok:
        return []
    if not body:
        probs.append(f"{label}: no data rows")
        return []
    out = []
    for ln, raw in enumerate(body, start=2):
        row = {}
        for col in item.columns:
            idx = hmap[col.name][0]
            cell = raw[idx].strip() if idx < len(raw) else ""
            if col.is_label:
                if not cell:
                    probs.append(f"{label}:{ln}: empty label '{col.name}'")
                row[col.name] = cell
            elif cell == NA_TOKEN and col in item.quantities:
                row[col.name] = None
            else:
                try:
                    v = float(cell)
                    if not math.isfinite(v):
                        raise ValueError
                    row[col.name] = v
                except ValueError:
                    probs.append(f"{label}:{ln}: '{col.name}' is not a finite number: {cell!r}")
                    row[col.name] = None
        for name, _unit in extra_cols:
            idx = hmap[name][0]
            row[name] = raw[idx].strip() if idx < len(raw) else ""
        out.append(row)
    return out


def _key(item, row):
    k = []
    for p in item.params:
        v = row.get(p.name)
        k.append(v if p.is_label or v is None else round(v, 9))
    return tuple(k)


def _matches(row, where):
    for k, v in where.items():
        rv = row.get(k)
        if isinstance(v, tuple):
            if not any(_matches(row, {k: vv}) for vv in v):
                return False
            continue
        if isinstance(v, str):
            if rv != v:
                return False
        elif rv is None or isinstance(rv, str) or not math.isclose(rv, v, rel_tol=1e-6, abs_tol=1e-9):
            return False
    return True


def _compare(value, ref, mode, tol):
    if mode == "abs":
        return abs(value - ref) <= tol, f"|diff| {abs(value - ref):.4g} (tol {tol:g})"
    if mode == "rel":
        err = abs(value - ref) / abs(ref)
        return err <= tol, f"rel {err:.3%} (tol {tol:.0%})"
    if mode == "mag_rel":
        err = abs(abs(value) - abs(ref)) / abs(ref)
        return err <= tol, f"|rel| {err:.3%} (tol {tol:.0%})"
    if mode == "factor":
        if value <= 0 or ref <= 0:
            return False, "non-positive value in a factor comparison"
        f = max(value / ref, ref / value)
        return f <= tol, f"factor {f:.3g} (tol {tol:g})"
    raise ValueError(mode)


def _note_fields(text):
    fields = {}
    for line in text.splitlines():
        m = re.match(r"^\s*[-*]?\s*\**([A-Za-z ]+?)\**\s*:\s*\**\s*(.*)$", line)
        if m:
            fields.setdefault(m.group(1).strip().lower(), m.group(2).strip())
    return fields


def _verdict(fields):
    v = fields.get("verdict", "")
    m = re.match(r"^\**\s*(PASS|DEVIATION|BLOCKED)\b", v)
    return m.group(1) if m else None


# ---------------------------------------------------------------- validator

@dataclass
class Outcome:
    item: str
    state: str                 # PENDING | BLOCKED | VALID | INVALID
    problems: list
    info: list


DELIVERY_FILES = ("results.csv", "convergence.csv", "note.md", "model.java", "params.csv", "run.log")


def validate_item(item: Item, folder: Path) -> Outcome:
    probs, info = [], []
    if not (folder / "README.md").is_file():
        probs.append("README.md missing")
    present = [f for f in DELIVERY_FILES if (folder / f).is_file()]
    if not present:
        return Outcome(item.id, "INVALID" if probs else "PENDING", probs, info)

    note_path = folder / "note.md"
    note_text = note_path.read_text(encoding="utf-8") if note_path.is_file() else ""
    fields = _note_fields(note_text)
    verdict = _verdict(fields)

    # BLOCKED: note only, no results.
    if verdict == "BLOCKED" and not (folder / "results.csv").is_file():
        for f in NOTE_FIELDS_BLOCKED:
            if not fields.get(f.lower()):
                probs.append(f"note.md: missing field '{f}'")
        if "date" in fields and not re.search(r"\d{4}-\d{2}-\d{2}", fields["date"]):
            probs.append("note.md: Date is not YYYY-MM-DD")
        if len(re.sub(r"^\**\s*BLOCKED\W*", "", fields.get("verdict", ""))) < 3:
            probs.append("note.md: BLOCKED verdict gives no reason")
        return Outcome(item.id, "INVALID" if probs else "BLOCKED", probs, info)

    required = ("results.csv", "note.md") if not item.solve else DELIVERY_FILES
    for f in required:
        if not (folder / f).is_file():
            probs.append(f"{f} missing")

    # note.md provenance
    if note_path.is_file():
        for f in NOTE_FIELDS:
            if not fields.get(f.lower()):
                probs.append(f"note.md: missing field '{f}'")
        if fields.get("date") and not re.search(r"\d{4}-\d{2}-\d{2}", fields["date"]):
            probs.append("note.md: Date is not YYYY-MM-DD")
        if fields.get("constants") and not re.search(r"\[(V|DR|E|A)\]", fields["constants"]):
            probs.append("note.md: Constants carries no evidence tag [V]/[DR]/[E]/[A]")
        if fields.get("verdict") and verdict is None:
            probs.append("note.md: Verdict is not PASS / DEVIATION / BLOCKED")
        if not re.search(r"^\s*\|.*\|\s*$", note_text, re.M):
            probs.append("note.md: no comparison table against the module value")

    # results.csv
    rows = []
    if (folder / "results.csv").is_file():
        rows = _check_table(item, folder / "results.csv", [("mesh_level", ""), ("evidence_tag", "")], probs)
        for i, r in enumerate(rows, start=2):
            if r["evidence_tag"] not in EVIDENCE_TAGS:
                probs.append(f"results.csv:{i}: evidence_tag {r['evidence_tag']!r} not in {sorted(EVIDENCE_TAGS)}")
            if r["mesh_level"] not in MESH_LEVELS and item.solve:
                probs.append(f"results.csv:{i}: mesh_level {r['mesh_level']!r} not in {MESH_LEVELS}")
        for q in item.quantities:
            if rows and all(r.get(q.name) is None for r in rows):
                probs.append(f"results.csv: quantity '{q.name}' has no numeric value")

    # convergence.csv
    if item.solve and (folder / "convergence.csv").is_file() and rows:
        crow = _check_table(item, folder / "convergence.csv",
                            [("mesh_level", ""), ("n_elements", "1"), ("run_time", "s")], probs)
        by = {}
        for r in crow:
            by.setdefault(_key(item, r), {})[r["mesh_level"]] = r
        for r in rows:
            k = _key(item, r)
            pair = by.get(k, {})
            if not all(lv in pair for lv in MESH_LEVELS):
                probs.append(f"convergence.csv: parameter point {k} lacks both mesh levels {MESH_LEVELS}")
                continue
            a, b = pair["production"], pair["refined"]
            try:
                if float(b["n_elements"]) <= float(a["n_elements"]):
                    probs.append(f"convergence.csv: {k} refined mesh has no more elements than production")
            except ValueError:
                probs.append(f"convergence.csv: {k} n_elements not numeric")
            for q in item.quantities:
                va, vb = a.get(q.name), b.get(q.name)
                if (va is None) != (vb is None):
                    probs.append(f"convergence.csv: {k} '{q.name}' NA at one mesh level only")
                    continue
                if va is None:
                    continue
                if q.unit in ENERGY_UNITS_MEV:
                    d = abs(va - vb) * ENERGY_UNITS_MEV[q.unit]
                    if d >= CONV_ENERGY_MEV:
                        probs.append(f"convergence.csv: {k} '{q.name}' not converged: {d:.3g} meV >= {CONV_ENERGY_MEV} meV")
                else:
                    den = max(abs(va), abs(vb))
                    rel = 0.0 if den == 0 else abs(va - vb) / den
                    if rel >= CONV_REL:
                        probs.append(f"convergence.csv: {k} '{q.name}' not converged: {rel:.2%} >= {CONV_REL:.0%}")

    # cross-checks against the surrogate numbers
    deviations = 0
    if rows and not probs:
        for ref in item.refs:
            for r in rows:
                if not _matches(r, ref.where) or r.get(ref.quantity) is None:
                    continue
                try:
                    rv = ref.ref(r) if callable(ref.ref) else ref.ref
                except Exception as exc:  # module import / range error is reported, not hidden
                    probs.append(f"cross-check {ref.quantity}: surrogate call failed: {exc!r}")
                    continue
                ok, msg = _compare(r[ref.quantity], rv, ref.mode, ref.tol)
                if ref.diagnostic:
                    tag = "diagnostic agree" if ok else "DIAGNOSTIC-MISS (module finding, not a COMSOL failure)"
                else:
                    tag = "agree" if ok else "DEVIATION"
                    deviations += 0 if ok else 1
                info.append(f"{ref.quantity} @ {_key(item, r)}: {r[ref.quantity]:.6g} vs {rv:.6g} "
                            f"({ref.source}) {msg} -> {tag}")
        if deviations and verdict == "PASS":
            probs.append(f"note.md claims PASS but {deviations} cross-check(s) fall outside the pass criterion")
        if item.extra:
            probs.extend(item.extra(item, rows, note_text))
    return Outcome(item.id, "INVALID" if probs else "VALID", probs, info)


def walk(root: Path):
    """Validate every work-order item folder under root; report strays."""
    outcomes, strays = [], []
    for it in ITEMS:
        folder = root / it.tier / it.id
        if not folder.is_dir():
            outcomes.append(Outcome(it.id, "INVALID", [f"folder {folder} missing"], []))
        else:
            outcomes.append(validate_item(it, folder))
    known = {(it.tier, it.id) for it in ITEMS}
    if root.is_dir():
        for tier_dir in sorted(p for p in root.iterdir() if p.is_dir()):
            for d in sorted(p for p in tier_dir.iterdir() if p.is_dir()):
                if (tier_dir.name, d.name) not in known:
                    strays.append(str(d.relative_to(root)))
    return outcomes, strays


# ---------------------------------------------------------------- README

def readme_text(it: Item) -> str:
    L = [f"# {it.id}: {it.title}", "",
         f"Work order: `{WORK_ORDER}`, section {it.section}, item {it.id}.",
         f"Status: PENDING until the COMSOL agent delivers the files below "
         f"(`verify/verify_comsol_anchors.py` treats a README-only folder as PENDING).", "",
         f"- Physics / interfaces: {it.physics}",
         f"- Replaces or validates: {it.replaces}",
         f"- Compared against: {it.compare_to}",
         f"- Anchor: {it.anchor}",
         f"- Pass criterion (work order): {it.pass_rule}", "",
         "## Files to deliver (work order section 2)", ""]
    if it.solve:
        L += ["`model.java`, `params.csv`, `results.csv`, `convergence.csv`, `run.log`, `note.md`. "
              "Keep .mph files over 50 MB outside git and record their hash in the note.", ""]
    else:
        L += ["No solve: `results.csv` and `note.md` only (no mesh, so no `convergence.csv`).", ""]
    L += ["## results.csv schema", "",
          "One row per parameter combination. Header cells are written `name [unit]` exactly as below; "
          "the harness fails a missing column or a different unit. Quantity cells may be `NA` only where "
          "the quantity does not apply to that row.", "",
          "| column header | role | meaning |", "|---|---|---|"]
    for c in it.params:
        L.append(f"| `{c.header}` | parameter | {c.desc} |")
    for c in it.quantities:
        L.append(f"| `{c.header}` | quantity | {c.desc} |")
    L += ["| `mesh_level` | mesh | `production` or `refined` |",
          "| `evidence_tag` | provenance | `[DR]` (field solve) or `[V]` (reproduces a published measurement) |", ""]
    if it.solve:
        L += ["## convergence.csv schema", "",
              "The same parameter and quantity columns as results.csv, plus `mesh_level` "
              "(`production` and `refined`, one uniform refinement), `n_elements [1]` and `run_time [s]`. "
              "Every results row must appear at both mesh levels; the refined mesh must have more elements.", "",
              "Convergence criterion (work order section 1): every quantity changes by less than 2 percent "
              "between the two levels; energies (meV, eV) by less than 1 meV.", ""]
    L += ["## note.md contents", "",
          "One page with these `Field: value` lines (the harness checks each):", ""]
    L += [f"- {f}:" for f in NOTE_FIELDS]
    L += ["",
          "`Solver version` is the COMSOL version and build; `Operator` is who ran it; `Date` is YYYY-MM-DD; "
          "`Constants` quotes the material constants used with their tags ([V]/[DR]/[E]/[A]) and source module; "
          "`Verdict` is one line `PASS`, `DEVIATION` or `BLOCKED` followed by the reason "
          "(`BLOCKED: licence <interface>` needs only Operator, Date and Verdict). "
          "Include a markdown comparison table: COMSOL value, module value, ratio.", ""]
    if it.refs:
        L += ["## Automated cross-checks in the harness", ""]
        for r in it.refs:
            where = ", ".join(f"{k}={v}" for k, v in r.where.items()) or "every row"
            val = "computed per row" if callable(r.ref) else f"{r.ref:g}"
            kind = " [diagnostic of the module; never fails a PASS]" if r.diagnostic else ""
            L.append(f"- `{r.quantity}` at {where}: {val} ({r.source}), {r.mode} tolerance {r.tol:g}{kind}")
        L += ["", "A mismatch is reported as DEVIATION (a result); it fails the harness only if note.md claims PASS.", ""]
    return "\n".join(L)


TOP_README = """# COMSOL anchor intake

Delivery folders for `docs/comsol_work_order_nanowire.md`, one per work-order item:
`comsol/<tier>/<item-id>/`, with tiers `rt_edge` (section 3, E1-E5), `nitride_cavity`
(section 4, P1-P5), `nitride_nanowire` (section 5, N1-N7) and `cross_cutting`
(section 6, X1-X2). Each item's README states the results.csv, convergence.csv and
note.md contract. `python verify/verify_comsol_anchors.py` validates every delivered
folder; README-only folders are PENDING. The READMEs are generated from the harness
table with `python verify/verify_comsol_anchors.py --write-readmes`.
"""


def write_readmes(root: Path = COMSOL_DIR):
    root.mkdir(exist_ok=True)
    (root / "README.md").write_text(TOP_README, encoding="utf-8")
    for it in ITEMS:
        d = root / it.tier / it.id
        d.mkdir(parents=True, exist_ok=True)
        (d / "README.md").write_text(readme_text(it), encoding="utf-8")
    print(f"wrote {len(ITEMS)} item READMEs under {root}")


# ---------------------------------------------------------------- checks

RESULTS: list[bool] = []


def check(name, ok, detail=""):
    RESULTS.append(bool(ok))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f" -- {detail}" if detail else ""))
    return bool(ok)


def surrogate_sanity():
    """New expectation checks: closed forms the code did not produce."""
    print("== surrogate numbers the work order names (closed form vs fsim_core) ==")
    ec = sphere_EC_meV(12.5, 9.5)
    check("closed-form E_C(12.5 nm, eps_r 9.5) = 12.126 meV (work order P5)", abs(ec - 12.126) < 5e-3,
          f"e^2/(4 pi eps0 eps_r R) = {ec:.4f} meV")
    for T, ref in ((230.0, 0.612), (300.0, 0.469)):
        r = ec * 1e-3 * E_CHARGE / (KB * T)
        check(f"closed-form E_C/kT at {T:.0f} K = {ref} (work order P5)", abs(r - ref) < 1e-3, f"{r:.4f}")
    try:
        mod = _module_EC({"island_radius": 12.5})
        check("drive_mech.set_feasibility E_C agrees with the closed form", abs(mod - ec) < 1e-6,
              f"{mod:.6f} vs {ec:.6f} meV")
    except Exception as exc:
        check("drive_mech.set_feasibility E_C agrees with the closed form", False, repr(exc))
    # Spontaneous-only field dP_sp / (eps0 eps_r) with the Bernardini constants
    # quoted in work order section 1 [V] (Bernardini et al., PRB 56, R10024
    # (1997)): P_sp GaN -0.029, InN -0.032 C/m2; eps_r GaN 10.28, InN 14.61;
    # linear interpolation [A].
    for x, wo in ((0.25, 74.5), (0.40, 112.8)):
        psp = (1 - x) * -0.029 + x * -0.032
        eps = (1 - x) * 10.28 + x * 14.61
        f = (-0.029 - psp) / (EPS0 * eps) * 1e-5
        check(f"closed-form spontaneous-only field x={x} = {wo} kV/cm (work order N1)", abs(f - wo) < 0.1,
              f"{f:.3f} kV/cm")
        try:
            mod = _module_pol_field({"x_in": x, "strain_fraction": 0.0})
            check(f"nitride_materials.polarization_field(x={x}, strain 0) matches the closed form",
                  abs(mod - f) / f < 1e-6, f"{mod:.4f} vs {f:.4f} kV/cm")
        except Exception as exc:
            check(f"nitride_materials.polarization_field(x={x}, strain 0) matches the closed form", False, repr(exc))
    # 2S/R at R 12.5 nm, S 1e3 cm/s caps the lifetime at 0.625 ns (work order N7).
    tau = 1.0 / _uniform_surface_rate_per_ns({"S": 1000.0, "core_radius": 12.5})
    check("closed-form 2S/R lifetime cap 0.625 ns at 12.5 nm, S 1e3 cm/s (work order N7)",
          abs(tau - 0.625) < 1e-9, f"{tau:.4f} ns")


def run_default():
    print(f"COMSOL anchor intake: {COMSOL_DIR}")
    print("== item folders ==")
    outcomes, strays = walk(COMSOL_DIR)
    for it, oc in zip(ITEMS, outcomes):
        rd = COMSOL_DIR / it.tier / it.id / "README.md"
        text = rd.read_text(encoding="utf-8") if rd.is_file() else ""
        names = (f"section {it.section}" in text and f"item {it.id}" in text
                 and all(c.header in text for c in it.columns))
        check(f"{it.tier}/{it.id} README names section, item id and schema", names)
        check(f"{it.tier}/{it.id} state {oc.state}", oc.state != "INVALID", "; ".join(oc.problems))
        for line in oc.info:
            print(f"        {line}")
    check("no stray item folders outside the work order", not strays, ", ".join(strays))
    surrogate_sanity()
    states = {}
    for oc in outcomes:
        states.setdefault(oc.state, []).append(oc.item)
    print("== summary ==")
    for s in ("PENDING", "BLOCKED", "VALID", "INVALID"):
        if s in states:
            print(f"  {s}: {', '.join(states[s])}")


# ---------------------------------------------------------------- self-test

GOOD_NOTE = """# P5 fixture note

- Solver version: COMSOL 6.2.0.290 (fixture)
- Geometry: conducting sphere in GaN, contacts at card distances
- Boundary conditions: island at 1 V, contacts grounded, infinite element domain
- Mesh: physics-controlled finer, one uniform refinement
- Licence modules: AC/DC (es)
- Run time: 42 s
- Operator: self-test fixture
- Date: 2026-09-23
- Constants: GaN eps_r 9.5 [A] (drive_mech convention); e, eps0 CODATA [V]
- Verdict: {verdict} fixture

| quantity | COMSOL | module | ratio |
|---|---|---|---|
| E_C (12.5 nm) | x | 12.126 | x |
"""


def _write_fixture(root: Path, *, scale=1.05, refined_drift=0.001, drop_col=None,
                   unit_override=None, note=True, verdict="PASS"):
    it = ITEM_BY_ID["P5"]
    d = root / it.tier / it.id
    d.mkdir(parents=True, exist_ok=True)
    (d / "README.md").write_text(readme_text(it), encoding="utf-8")
    for f in ("model.java", "params.csv", "run.log"):
        (d / f).write_text("fixture\n", encoding="utf-8")

    def hdr(cols):
        out = []
        for c in cols:
            if c.name == drop_col:
                continue
            out.append(f"{c.name} [{unit_override}]" if unit_override and c.name == "C_total" else c.header)
        return out

    def vals(R, env, s):
        ec = sphere_EC_meV(R, 9.5) * s
        c_aF = E_CHARGE / (ec * 1e-3) * 1e18
        v = {"island_radius": R, "environment": env, "C_total": c_aF, "E_C": ec,
             "EC_over_kT_230K": ec * 1e-3 * E_CHARGE / (KB * 230.0),
             "EC_over_kT_300K": ec * 1e-3 * E_CHARGE / (KB * 300.0)}
        return [f"{v[c.name]:.9g}" if not c.is_label else v[c.name] for c in it.columns if c.name != drop_col]

    points = [(r, "nanowire") for r in (0.5, 1.0, 5.0, 12.5)]
    with (d / "results.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(hdr(it.columns) + ["mesh_level", "evidence_tag"])
        for R, env in points:
            w.writerow(vals(R, env, scale) + ["production", "[DR]"])
    with (d / "convergence.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(hdr(it.columns) + ["mesh_level", "n_elements [1]", "run_time [s]"])
        for R, env in points:
            w.writerow(vals(R, env, scale) + ["production", "120000", "20"])
            w.writerow(vals(R, env, scale * (1 + refined_drift)) + ["refined", "960000", "150"])
    if note:
        (d / "note.md").write_text(GOOD_NOTE.format(verdict=verdict), encoding="utf-8")
    return it, d


def _write_rows_fixture(root: Path, item: Item, rows, verdict="PASS"):
    """Generic well-formed folder for ``item`` with the given production rows
    (dicts over item.columns; unspecified quantities default to 1.0)."""
    d = root / item.tier / item.id
    d.mkdir(parents=True, exist_ok=True)
    (d / "README.md").write_text(readme_text(item), encoding="utf-8")
    for f in ("model.java", "params.csv", "run.log"):
        (d / f).write_text("fixture\n", encoding="utf-8")

    def cells(r):
        return [r[c.name] if c.is_label else f"{r.get(c.name, 1.0):.9g}" for c in item.columns]

    with (d / "results.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([c.header for c in item.columns] + ["mesh_level", "evidence_tag"])
        for r in rows:
            w.writerow(cells(r) + ["production", "[DR]"])
    with (d / "convergence.csv").open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow([c.header for c in item.columns] + ["mesh_level", "n_elements [1]", "run_time [s]"])
        for r in rows:
            w.writerow(cells(r) + ["production", "120000", "20"])
            w.writerow(cells(r) + ["refined", "960000", "150"])
    (d / "note.md").write_text(GOOD_NOTE.replace("P5", item.id).format(verdict=verdict), encoding="utf-8")
    return d


def _n5_rows(vj300=3.09, vj300_off=2.40, vj10=None):
    """N5 fixture: 300 K rows at 0.5, 1, 2, 5 nA (only 1 and 2 nA carry the
    3.09 V criterion) and 10 K rows at 1 and 2 nA with the Deshpande rises."""
    eg10 = _gan_gap_eV({"T_bath": 10.0}) if vj10 is None else vj10
    rows = [{"T_bath": 300.0, "I": i, "V_j": vj300 if i in (1.0, 2.0) else vj300_off, "dT_junction": 1.0}
            for i in (0.5, 1.0, 2.0, 5.0)]
    rows += [{"T_bath": 10.0, "I": 1.0, "V_j": eg10, "dT_junction": 15.5},
             {"T_bath": 10.0, "I": 2.0, "V_j": eg10, "dT_junction": 47.0}]
    return rows


def _n2_rows(dperp_vertical=7.98 * 1.03):
    return [{"x_in": 0.25, "family": "horizontal", "core_radius": 12.5, "disc_radius": 12.5, "T": 25.0,
             "E_trans": 2777.0, "overlap": 0.94, "dE_perp_e": 13.93, "dE_perp_h": 3.31},
            {"x_in": 0.25, "family": "vertical", "core_radius": 100.0, "disc_radius": 12.5, "T": 25.0,
             "E_trans": 2777.0, "overlap": 0.94, "dE_perp_e": dperp_vertical, "dE_perp_h": 3.0}]


def run_self_test():
    print("== self-test: fixture item folders (P5) ==")
    cases = [
        ("well-formed item", {}, "VALID", None),
        ("missing column", {"drop_col": "E_C"}, "INVALID", "missing column"),
        ("non-converged mesh", {"refined_drift": 0.05}, "INVALID", "not converged"),
        ("missing note.md", {"note": False}, "INVALID", "note.md missing"),
        ("wrong units", {"unit_override": "fF"}, "INVALID", "wrong units"),
        ("PASS claimed on a deviating result", {"scale": 2.0}, "INVALID", "claims PASS"),
        ("DEVIATION honestly reported", {"scale": 2.0, "verdict": "DEVIATION"}, "VALID", None),
    ]
    with tempfile.TemporaryDirectory() as tmp:
        for i, (name, kw, want, reason) in enumerate(cases):
            root = Path(tmp) / f"case{i}"
            it, d = _write_fixture(root, **kw)
            oc = validate_item(it, d)
            ok = oc.state == want and (reason is None or any(reason in p for p in oc.problems))
            if want == "VALID" and "scale" not in kw:
                ok = ok and oc.info and all(line.endswith("agree") for line in oc.info)
            check(f"{name}: expect {want}" + (f" for '{reason}'" if reason else ""), ok,
                  f"got {oc.state}; " + "; ".join(oc.problems[:2]))
        # README-only folder is PENDING; BLOCKED note is BLOCKED; stray folder is caught.
        root = Path(tmp) / "pending"
        write_readmes_quiet(root)
        outcomes, strays = walk(root)
        check("README-only tree: every item PENDING", all(o.state == "PENDING" for o in outcomes) and not strays,
              ", ".join(f"{o.item}={o.state}" for o in outcomes if o.state != "PENDING"))
        bd = root / "nitride_nanowire" / "N3"
        (bd / "note.md").write_text("- Operator: fixture\n- Date: 2026-09-23\n- Verdict: BLOCKED: licence Wave Optics (emw)\n",
                                    encoding="utf-8")
        check("BLOCKED note without results: state BLOCKED", validate_item(ITEM_BY_ID["N3"], bd).state == "BLOCKED")
        (bd / "note.md").write_text("- Verdict: BLOCKED\n", encoding="utf-8")
        check("BLOCKED note missing operator/date/reason: INVALID",
              validate_item(ITEM_BY_ID["N3"], bd).state == "INVALID")
        (root / "nitride_nanowire" / "N9").mkdir()
        _, strays = walk(root)
        check("stray folder N9 is reported", strays == [str(Path("nitride_nanowire") / "N9")], str(strays))
        # N5 V_j filter: off-criterion V_j at 0.5 and 5 nA must not be compared.
        n5 = ITEM_BY_ID["N5"]
        oc = validate_item(n5, _write_rows_fixture(Path(tmp) / "n5a", n5, _n5_rows()))
        vj300 = [ln for ln in oc.info if ln.startswith("V_j @ (300.0")]
        diag_ok = [ln for ln in oc.info if ln.endswith("diagnostic agree")]
        check("N5 V_j(300 K) cross-check only at I = 1 and 2 nA (0.5 / 5 nA rows 0.69 V off are ignored)",
              oc.state == "VALID" and len(vj300) == 2 and all(ln.endswith("agree") for ln in vj300),
              f"got {oc.state}; {len(vj300)} V_j(300 K) lines; " + "; ".join(oc.problems[:2]))
        oc = validate_item(n5, _write_rows_fixture(Path(tmp) / "n5b", n5, _n5_rows(vj300=2.70)))
        check("N5 V_j(300 K) 0.39 V off at 1 and 2 nA with a PASS note: INVALID",
              oc.state == "INVALID" and any("claims PASS" in p for p in oc.problems), f"got {oc.state}")
        # N5 diagnostic: drift-diffusion V_j(10 K) within 0.1 V of E_g(10 K)/q.
        oc = validate_item(n5, _write_rows_fixture(Path(tmp) / "n5c", n5,
                                                   _n5_rows(vj10=_gan_gap_eV({"T_bath": 10.0}) - 0.25)))
        miss = [ln for ln in oc.info if "DIAGNOSTIC-MISS" in ln]
        check("N5 diagnostic V_j(10 K) vs E_g/q: agrees at E_g, flags a 0.25 V miss without failing a PASS note",
              len(diag_ok) == 2 and len(miss) == 2 and oc.state == "VALID",
              f"agree {len(diag_ok)}, miss {len(miss)}, state {oc.state}")
        # N2 vertical-family E_perp vs 7.98 meV.
        n2 = ITEM_BY_ID["N2"]
        oc = validate_item(n2, _write_rows_fixture(Path(tmp) / "n2a", n2, _n2_rows()))
        vert = [ln for ln in oc.info if ln.startswith("dE_perp_e") and "'vertical'" in ln]
        check("N2 vertical dE_perp_e 8.22 meV vs 7.98 meV (3 percent): agree, VALID",
              oc.state == "VALID" and len(vert) == 1 and vert[0].endswith("agree"),
              f"got {oc.state}; " + "; ".join(oc.problems[:2]))
        oc = validate_item(n2, _write_rows_fixture(Path(tmp) / "n2b", n2, _n2_rows(dperp_vertical=9.05)))
        check("N2 vertical dE_perp_e at the 9.05 meV hard-wall value (13 percent off) with a PASS note: INVALID",
              oc.state == "INVALID" and any("claims PASS" in p for p in oc.problems), f"got {oc.state}")
        # E3 duty-average rule: tau < 12.5 ns without a note sentence fails.
        check("E3 tau_th < 12.5 ns without a duty-average statement is flagged",
              bool(_tau_note_check(ITEM_BY_ID["E3"], [{"tau_th": 3.0}], "no statement")))


def write_readmes_quiet(root: Path):
    root.mkdir(parents=True, exist_ok=True)
    for it in ITEMS:
        d = root / it.tier / it.id
        d.mkdir(parents=True, exist_ok=True)
        (d / "README.md").write_text(readme_text(it), encoding="utf-8")


def main(argv):
    if "--write-readmes" in argv:
        write_readmes()
        return 0
    if "--self-test" in argv:
        run_self_test()
        label = "comsol anchor harness self-test checks"
    else:
        run_default()
        label = "comsol anchor checks"
    n, p = len(RESULTS), sum(RESULTS)
    print(f"\n{p}/{n} {label} passed")
    return 0 if p == n else 1


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
