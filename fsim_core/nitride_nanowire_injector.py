"""Piece 6 -- GaN/AlGaN tunnel injector: an independent, second deterministic-
pair feasibility screen (M-5a, Kitamura-class resonant-tunnelling hardware,
see fsim_core.drive_mech.mech_rti) for the nitride nanowire disc-in-wire
single-photon diode.

This module is standalone by design: it imports material constants from
fsim_core.nitride_materials (piece 1's frozen wurtzite parameters) but does
NOT import fsim_core.device, fsim_core.nitride_nanowire_levels, or
fsim_core.drive_mech, and does not touch DeviceDesign or nanowire result
types.  It prices a thin GaN/AlGaN tunnel barrier (single_barrier) or double
barrier with a GaN well (double_barrier) sitting in the axial i-segment that
feeds electrons and holes into the disc, as a SECOND, INDEPENDENT screen next
to the Coulomb-blockade island test (fsim_core.drive_mech.set_feasibility).
It answers a narrower question than that screen: "can this specific tunnel
structure deliver one electron and one hole into the dot's target levels,
selectively and on time, at the stated junction temperature" -- not whether
the disc itself blocks a second occupant electrostatically.

Physics and conventions
------------------------
* Energy zero: every energy this module consumes or returns for a given
  carrier (electron or hole) is referenced to THAT carrier's own bulk GaN
  band edge (Ec for electrons, Ev for holes, with the hole-energy axis
  defined so that hole energy increases moving away from the gap, i.e. hole
  potential = -Ev(x)).  This matches the contract's shared-zero requirement
  (docs/nitride_nanowire_contract.md) and lets injector_feasibility receive
  `electron_level_eV` / `hole_level_eV` from any dot-levels solver without
  this module importing that solver.
* Transmission: a numerically stable 1-D effective-mass transfer-matrix
  method (flux/mass-normalized, BenDaniel-Duke boundary condition), NOT a
  bare exp(-2 kappa d) WKB factor -- the WKB exponent only appears as the
  correct ASYMPTOTIC decay rate of this engine's own thick-barrier limit
  (checked in the verifier), and the engine natively resolves double-barrier
  resonances, which a bare WKB prefactor cannot represent at all.
* Rate normalization: capture/escape rate through the barrier stack is the
  single-channel Landauer/Buttiker rate, Gamma = (g_s/h) INTEGRAL T(E) f(E)
  dE [V; standard mesoscopic-transport result, e.g. Datta, "Electronic
  Transport in Mesoscopic Systems" (Cambridge, 1995), Ch. 2 -- the 1/h
  normalization is the same single-channel quantum that gives the
  conductance quantum e^2/h and is independent of the reservoir's effective
  mass or dispersion].  g_s (spin/valley degeneracy) is params.degeneracy
  [A, default 2].  f(E) is the injecting reservoir's Fermi-Dirac occupation;
  the receiving (dot) state is treated as empty (forward/loading rate) or
  full (reverse/reload rate) as documented per call site.
* Reservoir electrochemical energy: a bulk 3-D degenerate free-electron-gas
  Fermi level, E_F - E_band_edge = (hbar^2/2 m*)(3 pi^2 n)^(2/3) [DR;
  Ashcroft & Mermin, "Solid State Physics" (1976) Ch. 2 Eq. 2.33 -- a
  standard free-electron-gas result, not GaN-specific].  Using the T=0
  degenerate formula at finite T is itself an explicit [A] approximation
  (good to within roughly kT at these dopings, since E_F - E_edge is tens of
  meV and kT is 20-26 meV over 230-300 K); no non-degenerate/finite-T
  Fermi-Dirac inversion is attempted.
* Barrier material: an Al_xGa_1-xN barrier's conduction/valence BARRIER
  HEIGHTS come from a linear (no-bowing [A]) virtual-crystal interpolation
  of the GaN/AlN bandgap difference (nitride_materials' [V] Wu et al.
  bandgap()), PARTITIONED into electron/hole shares by the explicit
  `delta_Ev_GaN_AlN_eV` knob on NitrideNanowireInjectorParams -- default
  0.70 eV [V, Martin, Yu, Waldrop, APL 68, 2541 (1996), 0.70 +/- 0.24 eV;
  Rinke et al., PRB 77, 075202 (2008), ~0.8 eV].  This module does NOT use
  nitride_materials' own pinned Tsai & Bayram 0.30 eV valence-band offset
  (that module is read-only and its 0.30 eV value is a cited [V] number
  pinned by verify_nitride_materials.py) -- 0.30 eV is retained here only
  as the ALTERNATIVE partition, reported as the non-gating sensitivity
  output `rti_bypass_fraction_tsai_partition`.  Effective masses use a
  separate linear mass VCA between the same GaN/AlN endpoints.  Barrier
  bowing is neglected [A].  The barrier is evaluated UNSTRAINED
  (strain_fraction=0): this piece does not model the barrier as
  pseudomorphically strained to the GaN reservoir lattice; that is a
  conservative simplifying choice (coherency strain would shift the offset
  by an amount this piece does not attempt to sign) and is recorded under
  the module's `provenance` string and the worker's STATUS `decisions:`
  field.
* Growth realism: growth_tolerance_steps is a TOLERANCE (in units of the
  GaN c-axis bilayer spacing c/2, taken from nitride_materials' [V] lattice
  constant), not a commensurability requirement -- a nominal thickness is
  growth-feasible only if it lies within growth_tolerance_steps growth
  steps of its nearest integer-bilayer thickness AND the feasibility
  conjunction still holds when the barrier stack is perturbed by +/-
  growth_tolerance_steps growth steps in BOTH directions.  A design that
  only clears its screens at an exact nominal thickness, or whose nominal
  thickness itself sits far from any growable bilayer count, is reported
  fragile (rti_growth_feasible=False), never claimed manufacturable.
* Carrier default: transmission() and reflection() default to
  carrier="electron" -- a call that omits `carrier` prices ELECTRONS ONLY.
  Any device-integration caller pricing the hole path must pass
  carrier="hole" explicitly; injector_feasibility itself always prices both
  paths explicitly and does not rely on this default.

Honesty tags: no branch in this module can set rti_feasible=True without
computed, finite, out-of-branch-independent energetic/rate numbers passing
their stated thresholds; params.occupancy_control_known and
params.second_pair_control_known are necessary but not sufficient flags --
even with both True, rti_feasible still requires the numeric second-pair and
missed-load probabilities computed from real inputs to clear their
thresholds.  This module never reuses fsim_core.drive_mech.mech_rti's F_p
boost as evidence of determinism, and does not touch photon counting.
"""
from __future__ import annotations

import math
import warnings
from collections import namedtuple
from dataclasses import dataclass, replace

import numpy as np
from scipy.integrate import quad, IntegrationWarning
from scipy.optimize import brentq, minimize_scalar

from . import nitride_materials as NM

# ------------------------------------------------------------------ constants

HBAR_JS = 1.054571817e-34      # [V] CODATA 2018
H_EVS = 4.135667696e-15        # [V] CODATA 2018, h in eV s
KB_EV = 8.617333262e-5         # [V] CODATA 2018, k_B in eV/K
M0_KG = 9.1093837015e-31       # [V] CODATA 2018
EV_J = 1.602176634e-19         # [V] CODATA 2018, exact by SI definition
H_JS = H_EVS * EV_J             # [V] CODATA 2018, h in J s (derived from H_EVS)
_MG_ACCEPTOR_DEGENERACY = 4.0   # [E] standard shallow-acceptor degeneracy (spin x
                                 # valence-band-top degeneracy), e.g. Gotz et al.,
                                 # APL 68, 667 (1996); Kozodoy et al., JAP 87, 1832 (2000)

_GAN = NM.binary("GaN")
_ALN = NM.binary("AlN")
_REF_T_K = 300.0
# [A] Barrier band-offset reference temperature for the flat-barrier
# envelope. transmission() takes no T_K argument; the Varshni shift of the
# ~2.4 eV AlN-GaN offset changes by well under 0.1% over 230-300 K, which is
# negligible next to the [A] linear-interpolation and [A] partition
# uncertainties already carried by the offset itself.
_EG_GAN_EV = NM.bandgap(_GAN, _REF_T_K)   # [V] Wu et al., APL 2002, via nitride_materials
_EG_ALN_EV = NM.bandgap(_ALN, _REF_T_K)   # [V] Wu et al., APL 2002, via nitride_materials
_GAN_C_NM = _GAN.c_A / 10.0
_GROWTH_STEP_NM_DEFAULT = _GAN_C_NM / 2.0
# [V] GaN c-axis bilayer spacing c/2; c_A from nitride_materials' Bernardini
# PRB 56, R10024 (1997) Table II lattice constant.

_SLICE_LENGTH_M = 2.0e-10        # 0.2 nm spatial staircase step for field tilt
_MIN_SLICES_PER_SEGMENT = 8


_TSAI_BAYRAM_DELTA_EV_EV = 0.30
# [V] Tsai & Bayram, ACS Omega 5, 3917 (2020), Table 2: AlN valence 0.30 eV
# below GaN. This is the SAME number nitride_materials._VBO_ALN_EV pins for
# the shared strained-band-edge solver (that module is read-only and its
# 0.30 eV value is not touched here); it is retained below only as the
# ALTERNATIVE partition for the rti_bypass_fraction_tsai_partition
# sensitivity output, never as this module's default.


def _layer_polarization_Cm2_and_eps(al_fraction: float):
    """Natural-sign (Bernardini '+c positive', never abs()) pseudomorphic
    Al_xGa_1-xN-on-GaN polarization/permittivity, and the bulk unstrained
    GaN polarization/permittivity used by the undoped GaN well/reservoir.

    [DR] The AlGaN P_sp and P_pz are VCA interpolations of Bernardini,
    Fiorentini & Vanderbilt, PRB 56, R10024 (1997), Table II.  The alloy is
    constrained to the GaN in-plane lattice constant (pseudomorphic barrier
    on relaxed GaN reservoirs/well)."""
    aln = _ALN
    a = _GAN.a_A + al_fraction * (aln.a_A - _GAN.a_A)
    c13 = _GAN.C13_GPa + al_fraction * (aln.C13_GPa - _GAN.C13_GPa)
    c33 = _GAN.C33_GPa + al_fraction * (aln.C33_GPa - _GAN.C33_GPa)
    e31 = _GAN.e31_Cm2 + al_fraction * (aln.e31_Cm2 - _GAN.e31_Cm2)
    e33 = _GAN.e33_Cm2 + al_fraction * (aln.e33_Cm2 - _GAN.e33_Cm2)
    psp = _GAN.Psp_Cm2 + al_fraction * (aln.Psp_Cm2 - _GAN.Psp_Cm2)
    eps_algan = _GAN.eps_r + al_fraction * (aln.eps_r - _GAN.eps_r)
    ep = (_GAN.a_A - a) / a
    ez = -2.0 * c13 / c33 * ep
    p_algan = psp + 2.0 * e31 * ep + e33 * ez
    return p_algan, eps_algan, _GAN.Psp_Cm2, _GAN.eps_r


def _stack_polarization_fields_eV_per_m(al_fraction: float, segments, polarity: str):
    """Signed polarization field (eV/m, i.e. numerically V/m) for each
    (length_m, V_eV, m_ratio) segment of ONE carrier path's stack (HIGH
    3/4/5 fix, Opus re-review 2026-09-14): the undoped stack sits between
    two reservoirs that PIN THE POTENTIAL AT BOTH ENDS, so D (displacement)
    is a single constant across every interface (no free charge anywhere
    in the undoped stack -- Gauss's law) and the net potential drop across
    the WHOLE stack is zero: sum_i d_i E_i = 0 with E_i = (D - P_i)/(eps0
    eps_i) per layer [DR].  AlGaN segments (V_eV > 0, i.e. this path's
    barrier layers) use the pseudomorphic-on-GaN P/eps; GaN segments
    (V_eV == 0, the well) use bulk unstrained GaN P/eps.  For the symmetric
    2 nm / 4 nm / 2 nm Al0.30 design this gives E_b = +/-1.4755 MV/cm
    (0.2951 eV per 2 nm barrier) and E_w = -/+1.4755 MV/cm (-/+0.5902 eV
    across the 4 nm well) -- verified against the hand-computed closure in
    verify_nitride_nanowire_injector.py.  A single_barrier stack (one
    segment only) has nothing to balance against, so the closure trivially
    gives E = 0 (no net drop possible across a single, already-pinned-at-
    both-ends layer); this is a real consequence of the "pinned at both
    ends" assumption, not a special case coded separately.

    Sign: the natural (+c-positive, Bernardini) sign of this closure
    already gives the physically correct LOCAL (entry-to-exit) slope for
    BOTH carrier paths under Ga-polar material, with no further
    carrier-dependent flip.  Proof sketch (module docstring "Growth
    polarity and transport direction" has the full derivation): for
    Ga-polar growth the electron potential Ec(z) and the valence edge
    Ev(z) shift together with the SAME sign of slope in absolute +c
    coordinates (a uniform electrostatic tilt moves both bands together);
    the hole's own energy axis is defined as -Ev(z) (this module's
    convention, see module docstring "Energy zero"), and the hole's own
    LOCAL transport coordinate runs opposite to +c (holes travel toward
    -c) -- two sign flips that cancel, leaving the hole's local slope
    equal in sign to the electron's.  polarity="N" reverses the bound
    sheet-charge sign and therefore this closure's sign, for both paths
    identically.  never abs()."""
    if polarity not in ("Ga", "N"):
        raise ValueError("polarity must be 'Ga' or 'N'")
    p_algan, eps_algan, p_gan, eps_gan = _layer_polarization_Cm2_and_eps(al_fraction)
    layer_info = []
    num = 0.0
    den = 0.0
    for (length_m, V_eV, _m_ratio) in segments:
        if V_eV > 0.0:
            p_layer, eps_layer = p_algan, eps_algan
        else:
            p_layer, eps_layer = p_gan, eps_gan
        layer_info.append((length_m, p_layer, eps_layer))
        if length_m > 0.0:
            num += length_m * p_layer / eps_layer
            den += length_m / eps_layer
    if den <= 0.0:
        return [0.0 for _ in segments]
    D = num / den
    sign = 1.0 if polarity == "Ga" else -1.0
    return [sign * (D - p_layer) / (NM.EPS0_SI * eps_layer) for (_l, p_layer, eps_layer) in layer_info]


def _n_v_valence_dos_m3(mh_ratio: float, T_K: float) -> float:
    """3-D valence-band effective density of states, N_V =
    2 (2 pi m_h* k_B T / h^2)^(3/2) [V; standard semiconductor-physics
    result, e.g. Sze & Ng, "Physics of Semiconductor Devices," 3rd ed.
    (2007), Ch. 1]. mh_ratio in units of the free electron mass."""
    m = mh_ratio * M0_KG
    kT_J = KB_EV * T_K * EV_J
    return 2.0 * (2.0 * math.pi * m * kT_J / (H_JS ** 2)) ** 1.5


def _hole_quasi_fermi_eV(params: "NitrideNanowireInjectorParams", T_K: float):
    """Ionised free-hole density and the p-GaN reservoir's quasi-Fermi
    level (this module's own positive-away-from-the-gap hole-energy axis,
    see module docstring "Energy zero") from Mg-acceptor charge neutrality
    -- MEDIUM 6 fix, replaces the earlier Na*exp(-E_A/kT), which is not an
    ionization fraction.  Mass-action law for a single acceptor level with
    no compensating donors (standard incomplete-ionization result, e.g.
    Sze & Ng, "Physics of Semiconductor Devices," 3rd ed. (2007), Ch. 1):
    p^2 / (Na - p) = (N_V/g) exp(-E_A/kT), solved in closed form as a
    quadratic in p.  N_V uses mh_well (GaN heavy-hole-axis mass); E_A =
    params.mg_acceptor_energy_meV [E, Gotz et al., APL 68, 667 (1996);
    Kozodoy et al., JAP 87, 1832 (2000)]; g = 4 [E] standard shallow
    acceptor degeneracy.  Returns (p_free_cm3, mu_h_eV) where mu_h_eV =
    kT ln(p_free/N_V) (negative for a non-degenerate, partially-ionized
    reservoir -- this is expected and is NOT clamped positive: a negative
    value here correctly represents a dilute hole gas, unlike the old
    always-positive degenerate-formula estimate)."""
    kT_eV = KB_EV * T_K
    Na_m3 = params.p_cm3 * 1.0e6
    NV_m3 = _n_v_valence_dos_m3(params.mh_well, T_K)
    g = _MG_ACCEPTOR_DEGENERACY
    K = (NV_m3 / g) * math.exp(-params.mg_acceptor_energy_meV / 1000.0 / kT_eV)
    p_m3 = (-K + math.sqrt(K * K + 4.0 * K * Na_m3)) / 2.0
    p_m3 = max(p_m3, 1e-6)   # numerical floor only; physically p > 0 always
    mu_h_eV = kT_eV * math.log(p_m3 / NV_m3)
    return p_m3 / 1.0e6, float(mu_h_eV)


def _bandgap_diff_eV(al_fraction: float) -> float:
    """Linear (no-bowing [A]) virtual-crystal Al_xGa_1-xN - GaN bandgap
    difference at the barrier reference temperature, from nitride_materials'
    [V] Wu et al. bandgap() -- independent of how that gap is PARTITIONED
    into conduction/valence offsets (see _default_barrier_material)."""
    if not math.isfinite(al_fraction) or not 0.0 <= al_fraction <= 1.0:
        raise ValueError("al_fraction must be finite and in [0, 1]")
    return (_EG_ALN_EV - _EG_GAN_EV) * al_fraction


def _default_barrier_material(al_fraction: float, delta_Ev_eV: float) -> dict:
    """Al_xGa_1-xN barrier conduction/valence BARRIER HEIGHTS (relative to
    GaN, unstrained) and effective masses.

    The valence (hole) barrier height is x * delta_Ev_eV, an EXPLICIT
    partition of the VCA bandgap difference (see _bandgap_diff_eV) that
    this module owns -- NOT nitride_materials' own pinned Tsai & Bayram
    0.30 eV VBO (nitride_materials.py is read-only; see module docstring
    and NitrideNanowireInjectorParams.delta_Ev_GaN_AlN_eV). The conduction
    (electron) barrier height takes the remainder of the gap difference,
    dEc = dEg(x) - x*delta_Ev_eV, so Ec-Ev reproduces the VCA gap exactly
    for any choice of delta_Ev_eV. Effective masses are an unrelated linear
    mass VCA between the same GaN/AlN endpoints. [A] interpolation of [V]
    endpoints (see module docstring)."""
    if not math.isfinite(delta_Ev_eV) or delta_Ev_eV < 0.0:
        raise ValueError("delta_Ev_eV must be non-negative and finite")
    dEg = _bandgap_diff_eV(al_fraction)
    hole_barrier_eV = delta_Ev_eV * al_fraction
    dEc = dEg - hole_barrier_eV
    me = _GAN.me_z + (_ALN.me_z - _GAN.me_z) * al_fraction
    mh = _GAN.mh_z + (_ALN.mh_z - _GAN.mh_z) * al_fraction
    return dict(dEc_eV=float(dEc), dEv_eV=float(-hole_barrier_eV), me=float(me), mh=float(mh))


def _degenerate_mu_eV(n_m3: float, m_ratio: float) -> float:
    """Bulk 3-D degenerate free-electron-gas Fermi level above the band
    edge [DR; Ashcroft & Mermin, "Solid State Physics" (1976) Ch. 2 Eq.
    2.33]. n_m3 in carriers/m^3, m_ratio in units of the free electron mass."""
    if not math.isfinite(n_m3) or n_m3 <= 0:
        raise ValueError("n_m3 must be positive and finite")
    if not math.isfinite(m_ratio) or m_ratio <= 0:
        raise ValueError("m_ratio must be positive and finite")
    k_f = (3.0 * math.pi ** 2 * n_m3) ** (1.0 / 3.0)
    return (HBAR_JS ** 2 * k_f ** 2) / (2.0 * m_ratio * M0_KG) / EV_J


# --------------------------------------------------------------------- params

@dataclass(frozen=True)
class NitrideNanowireInjectorParams:
    """Frozen injector design. Electron and hole injection paths are priced
    independently (own topology, thickness, well width); by default they
    share the same physical Al_xGa_1-xN barrier composition (al_fraction),
    which is the ordinary case of one physical barrier stack probed by two
    carrier species with different mass/offset, but every per-path field can
    be overridden to model a genuinely asymmetric design.

    Starting designed stack [A]: Al fraction 0.30, each barrier 2 nm,
    optional GaN well 4 nm; sensitivities barrier thickness {1,2,3} nm and
    Al fraction {0.2,0.3,0.4} are exercised by the caller (device/sweep
    piece), not enumerated here.
    """

    al_fraction: float = 0.30                        # [A] designed stack Al content
    electron_topology: str = "double_barrier"         # "single_barrier" | "double_barrier"
    hole_topology: str = "double_barrier"
    electron_barrier_thickness_nm: float = 2.0        # [A]
    hole_barrier_thickness_nm: float = 2.0            # [A]
    electron_well_width_nm: float = 4.0               # [A] GaN well (double_barrier only)
    hole_well_width_nm: float = 4.0                   # [A]

    growth_step_nm: float = _GROWTH_STEP_NM_DEFAULT   # [V] GaN c/2 bilayer spacing
    growth_tolerance_steps: float = 1.0               # [A] +/- N growth increments swept for robustness

    me_barrier_override: float | None = None          # None => derived from al_fraction
    mh_barrier_override: float | None = None
    dEc_eV_override: float | None = None
    dEv_eV_override: float | None = None
    me_well: float = _GAN.me_z                        # [V] Rinke PRB 2008 GaN conduction mass
    mh_well: float = _GAN.mh_z                        # [V] Rinke PRB 2008 GaN valence mass (z)

    delta_Ev_GaN_AlN_eV: float = 0.70                 # [V] Martin, Yu, Waldrop, APL 68, 2541 (1996),
                                                       # 0.70 +/- 0.24 eV; Rinke et al. PRB 77, 075202
                                                       # (2008) ~0.8 eV. Explicit valence-band-offset
                                                       # partition knob for the AlGaN barrier (see
                                                       # _default_barrier_material); the Tsai & Bayram
                                                       # 0.30 eV alternative pinned by nitride_materials
                                                       # (read-only) is reported separately as
                                                       # rti_bypass_fraction_tsai_partition.

    n_cm3: float = 3.0e18                             # [V] Deshpande et al. 2013 n-GaN doping anchor
    p_cm3: float = 5.0e17                             # [V] Deshpande et al. 2013 p-GaN doping anchor

    alignment_uncertainty_meV: float = 15.0           # [A] growth/doping-limited level-alignment uncertainty
    degeneracy: float = 2.0                           # [A] spin/valley degeneracy in the Landauer rate
    reservoir_state_count_e: float = 1.0              # [A] MEDIUM 6 fix: spin lives in
                                                       # degeneracy alone (see `degeneracy`
                                                       # above) -- default 1 per carrier;
                                                       # caller should replace with the
                                                       # levels module's own transverse
                                                       # (non-spin) channel count.
    reservoir_state_count_h: float = 1.0              # [A] see reservoir_state_count_e
    mg_acceptor_energy_meV: float = 170.0              # [E] Gotz et al., APL 68, 667 (1996); Kozodoy et al., JAP 87, 1832 (2000)
    include_polarization: bool = True                  # [A] pseudomorphic-barrier envelope enabled by default
    polarity: str = "Ga"                              # [A] Deshpande et al. 2013 catalyst-free PA-MBE c-axis
                                                       # growth is ASSUMED Ga-polar [A]; "N" is the opposite-sign
                                                       # sensitivity. Combined with the fixed axial p-i-n growth
                                                       # order (n-GaN first, electrons toward +c, holes toward
                                                       # -c) this sets the SIGN of the polarization tilt (HIGH
                                                       # 3/4/5 fix; see _stack_polarization_fields_eV_per_m and
                                                       # module docstring "Growth polarity and transport
                                                       # direction"); never abs().
    bypass_prefactor: float = 1.0                     # [A] declared prefactor on the thermionic (over-barrier) integral
    field_leverarm: float = 1.0                       # [A] fraction of the applied bias dropping across the injector
    slice_length_nm: float = 0.2                      # [A numerical budget] field-tilt staircase step (see
                                                       # _transmission_reflection_scalar); default matches the
                                                       # module's original fixed 0.2 nm slice.
    min_slices_per_segment: int = 8                   # [A numerical budget] floor on slices per flat-band
                                                       # segment; both knobs exist so the verifier can force a
                                                       # deliberately coarse staircase for the rti_numerics_ok
                                                       # falsification case (MEDIUM 7) without perturbing the
                                                       # default (bit-identical) resolution.

    alignment_tunable: bool = False                   # [A] can a post-growth bias retune the resonance onto the target level
    bias_tuning_range_meV: float = 0.0                # [A] maximum level shift deliverable by that bias, if alignment_tunable

    occupancy_control_known: bool = False             # necessary, not sufficient (see module docstring)
    second_pair_control_known: bool = False

    def __post_init__(self):
        for name in ("electron_topology", "hole_topology"):
            v = getattr(self, name)
            if v not in ("single_barrier", "double_barrier"):
                raise ValueError(f"{name} must be 'single_barrier' or 'double_barrier'")
        if self.polarity not in ("Ga", "N"):
            raise ValueError("polarity must be 'Ga' or 'N'")
        if not math.isfinite(self.al_fraction) or not 0.0 <= self.al_fraction <= 1.0:
            raise ValueError("al_fraction must be finite and in [0, 1]")
        positive = (
            "electron_barrier_thickness_nm", "hole_barrier_thickness_nm",
            "electron_well_width_nm", "hole_well_width_nm", "growth_step_nm",
            "me_well", "mh_well", "n_cm3", "p_cm3", "alignment_uncertainty_meV",
            "degeneracy", "reservoir_state_count_e", "reservoir_state_count_h", "mg_acceptor_energy_meV", "bypass_prefactor", "field_leverarm", "delta_Ev_GaN_AlN_eV",
            "slice_length_nm",
        )
        for name in positive:
            v = getattr(self, name)
            if not math.isfinite(v) or v <= 0:
                raise ValueError(f"{name} must be positive and finite")
        if not isinstance(self.min_slices_per_segment, int) or self.min_slices_per_segment <= 0:
            raise ValueError("min_slices_per_segment must be a positive integer")
        if not math.isfinite(self.growth_tolerance_steps) or self.growth_tolerance_steps < 0:
            raise ValueError("growth_tolerance_steps must be non-negative and finite")
        if not math.isfinite(self.bias_tuning_range_meV) or self.bias_tuning_range_meV < 0:
            raise ValueError("bias_tuning_range_meV must be non-negative and finite")
        for name in ("me_barrier_override", "mh_barrier_override"):
            v = getattr(self, name)
            if v is not None and (not math.isfinite(v) or v <= 0):
                raise ValueError(f"{name} must be positive and finite if given")
        for name in ("dEc_eV_override", "dEv_eV_override"):
            v = getattr(self, name)
            if v is not None and not math.isfinite(v):
                raise ValueError(f"{name} must be finite if given")


_ResolvedPath = namedtuple(
    "_ResolvedPath",
    "carrier topology thickness_nm well_nm m_barrier m_well barrier_height_eV mu_eV",
)


def _resolve_path(params: NitrideNanowireInjectorParams, carrier: str) -> _ResolvedPath:
    if carrier not in ("electron", "hole"):
        raise ValueError("carrier must be 'electron' or 'hole'")
    default = _default_barrier_material(params.al_fraction, params.delta_Ev_GaN_AlN_eV)
    me_b = params.me_barrier_override if params.me_barrier_override is not None else default["me"]
    mh_b = params.mh_barrier_override if params.mh_barrier_override is not None else default["mh"]
    dEc = params.dEc_eV_override if params.dEc_eV_override is not None else default["dEc_eV"]
    dEv = params.dEv_eV_override if params.dEv_eV_override is not None else default["dEv_eV"]
    if carrier == "electron":
        mu = _degenerate_mu_eV(params.n_cm3 * 1e6, params.me_well)
        return _ResolvedPath("electron", params.electron_topology,
                              params.electron_barrier_thickness_nm,
                              params.electron_well_width_nm,
                              me_b, params.me_well, dEc, mu)
    # LOW 11 fix (Opus re-review of 12b39cd, 2026-09-14): the hole path's own
    # degenerate free-electron-gas mu is DEAD -- every caller of a hole
    # emitter quasi-Fermi level uses _hole_quasi_fermi_eV's Mg-acceptor
    # mass-action solve instead (MEDIUM 6, injector_feasibility's
    # _screen_at), never this field.  Kept as NaN (not computed at all) so
    # nothing can mistake an unused degenerate-gas estimate for a
    # meaningful number; the field name is unchanged (mu_eV) since the
    # electron branch above still genuinely uses it.
    return _ResolvedPath("hole", params.hole_topology,
                          params.hole_barrier_thickness_nm,
                          params.hole_well_width_nm,
                          mh_b, params.mh_well, -dEv, float("nan"))


def _stack_segments(path: _ResolvedPath):
    """List of (length_m, V_eV, m_ratio) flat-band segments for one carrier
    path's barrier stack. A single tunnel barrier has no claimed resonant
    well; a double barrier brackets a GaN well of the path's own mass."""
    if path.topology == "single_barrier":
        return [(path.thickness_nm * 1e-9, path.barrier_height_eV, path.m_barrier)]
    if path.topology == "double_barrier":
        return [(path.thickness_nm * 1e-9, path.barrier_height_eV, path.m_barrier),
                (path.well_nm * 1e-9, 0.0, path.m_well),
                (path.thickness_nm * 1e-9, path.barrier_height_eV, path.m_barrier)]
    raise ValueError("topology must be 'single_barrier' or 'double_barrier'")


def _tilt_eV_per_m(field_kVcm: float, bias_V: float, field_leverarm: float, total_len_m: float) -> float:
    """Uniform additive field (eV/m, numerically equal to V/m) tilting the
    barrier stack. Positive values lower the potential on the OUTPUT (dot)
    side [A sign convention: positive bias/field is forward/injecting]."""
    f_field = field_kVcm * 1.0e5   # kV/cm -> V/m
    f_bias = (bias_V * field_leverarm / total_len_m) if total_len_m > 0 else 0.0
    return f_field + f_bias


def _segment_k(E_eV: float, V_eV: float, m_ratio: float) -> complex:
    """Complex wavevector (1/m) in a constant-potential segment. For E > V
    this is real and positive (propagating); for E < V numpy's principal
    complex sqrt branch places the result on the positive imaginary axis
    (evanescent), which the backward recursion below relies on."""
    val = 2.0 * m_ratio * M0_KG * (E_eV - V_eV) * EV_J
    return complex(np.sqrt(complex(val))) / HBAR_JS


def _transmission_reflection_scalar(segments, energy_eV: float, m_left: float,
                                     m_right: float, tilt_eV_per_m: float,
                                     v_right_eV: float = 0.0,
                                     pol_fields_eV_per_m=None,
                                     slice_length_m: float = _SLICE_LENGTH_M,
                                     min_slices_per_segment: int = _MIN_SLICES_PER_SEGMENT):
    """Flux/mass-normalized transmission and reflection through a stack of
    (length_m, V_eV, m_ratio) flat-band segments. The LEFT (emitter)
    reservoir stays at the module's energy zero (V=0); the RIGHT
    (collector) reservoir is referenced at v_right_eV, the field/bias-
    tilted band edge at the far end of the stack (see _transmission_scalar)
    -- pinning both reservoirs at V=0 regardless of tilt would reproduce
    the interior ramp but then snap back by tilt*L at the exit, an
    unphysical discontinuity that also means bias_V could never lower the
    collector (MEDIUM 5 fix). A linear field tilt is applied by subdividing
    each segment into thin slices (staircase approximation of a linear
    ramp; a standard treatment, e.g. Ridley, "Quantum Processes in
    Semiconductors").  `pol_fields_eV_per_m`, one signed field per segment
    (see _stack_polarization_fields_eV_per_m), is accumulated into a
    RUNNING potential offset that is NEVER reset between segments (HIGH
    3/4/5 fix: the polarization potential is continuous at every interface,
    matching D continuity -- only its SLOPE changes segment to segment,
    exactly as a piecewise-constant-field integral must).  Numerically
    stable backward coefficient recursion (BenDaniel-Duke effective-mass
    boundary condition: continuity of psi and (1/m) dpsi/dx), not a
    chained-matrix product, to avoid overflow for thick/deep barriers.
    Returns RAW (unclamped) T, R -- the public transmission()/reflection()
    wrappers clip to the physical [0, 1] range; keeping this internal
    engine unclamped lets the rti_numerics_ok diagnostic (MEDIUM 7) see
    genuine numerical overshoot instead of a value already forced sane."""
    if pol_fields_eV_per_m is None:
        pol_fields_eV_per_m = [0.0] * len(segments)
    slices = []
    x0 = 0.0
    pol_offset = 0.0   # running polarization potential; continuous across segments
    for (length_m, V_eV, m_ratio), E_pol in zip(segments, pol_fields_eV_per_m):
        if length_m <= 0.0:
            continue
        n = max(min_slices_per_segment, int(math.ceil(length_m / slice_length_m)))
        dl = length_m / n
        for i in range(n):
            xc = x0 + (i + 0.5) * dl
            local_x = (i + 0.5) * dl
            pol_here = pol_offset + E_pol * local_x
            # V_electron(x) = -e*phi(x) + flat-band offset, and E = -dphi/dx
            # (standard physics sign), so dV_electron/dx = +e*E(x): the
            # polarization contribution is ADDED (not subtracted, unlike
            # the independent bias/field tilt term above, which uses its
            # own established "positive tilt lowers the collector side"
            # convention) -- this is what gives a Ga-polar barrier's exit
            # face ABOVE its entry face (HIGH 3/4/5 fix; see module
            # docstring "Growth polarity and transport direction").
            slices.append((dl, V_eV - tilt_eV_per_m * xc + pol_here, m_ratio))
        pol_offset += E_pol * length_m
        x0 += length_m

    k_right = _segment_k(energy_eV, v_right_eV, m_right)
    A, B = 1.0 + 0j, 0.0 + 0j
    k_prev, m_prev = k_right, m_right
    for (dl, V_eV, m_ratio) in reversed(slices):
        k = _segment_k(energy_eV, V_eV, m_ratio)
        u = A + B
        v = (k_prev / m_prev) * (A - B)
        theta = k * dl
        if abs(k) < 1e-3:
            k = 1e-3 + 0j  # degenerate-energy guard; physically negligible width
        A_new = (u + v * m_ratio / k) / 2.0 * np.exp(-1j * theta)
        B_new = (u - v * m_ratio / k) / 2.0 * np.exp(1j * theta)
        A, B = A_new, B_new
        k_prev, m_prev = k, m_ratio

    k_left = _segment_k(energy_eV, 0.0, m_left)
    u = A + B
    v = (k_prev / m_prev) * (A - B)
    if abs(k_left) < 1e-12:
        return 0.0, 1.0
    A0 = (u + v * m_left / k_left) / 2.0
    B0 = (u - v * m_left / k_left) / 2.0
    if abs(A0) < 1e-300:
        return 0.0, 1.0
    v_left = k_left.real / m_left
    v_right = k_right.real / m_right
    if v_left <= 0.0 or v_right <= 0.0:
        return 0.0, 1.0
    T = (v_right / v_left) * (1.0 / abs(A0)) ** 2
    R = abs(B0 / A0) ** 2
    return float(T), float(R)


def _transmission_scalar(params, energy_eV, bias_V, field_kVcm, carrier):
    path = _resolve_path(params, carrier)
    segs = _stack_segments(path)
    total_len_m = sum(s[0] for s in segs)
    tilt = _tilt_eV_per_m(field_kVcm, bias_V, params.field_leverarm, total_len_m)
    # MEDIUM 5 fix: the collector (output) reservoir band edge is lowered
    # (for a positive/forward tilt) by the same tilt that ramps the barrier
    # interior, so the ramp is continuous into the collector instead of
    # snapping back to the emitter's V=0 at the exit face. The polarization
    # closure's net drop across the WHOLE stack is zero by construction
    # (see _stack_polarization_fields_eV_per_m), so v_right_eV needs no
    # separate polarization term.
    v_right_eV = -tilt * total_len_m
    if params.include_polarization:
        pol_fields = _stack_polarization_fields_eV_per_m(params.al_fraction, segs, params.polarity)
    else:
        pol_fields = [0.0] * len(segs)
    T, R = _transmission_reflection_scalar(
        segs, float(energy_eV), path.m_well, path.m_well, tilt, v_right_eV=v_right_eV,
        pol_fields_eV_per_m=pol_fields,
        slice_length_m=params.slice_length_nm * 1e-9,
        min_slices_per_segment=params.min_slices_per_segment)
    return T, R


def transmission(params: NitrideNanowireInjectorParams, energy_eV, *,
                  bias_V: float = 0.0, field_kVcm: float = 0.0,
                  carrier: str = "electron"):
    """1-D effective-mass transmission probability through the named
    carrier's barrier stack, at energy_eV referenced to that carrier's own
    GaN reservoir band edge (see module docstring). Both ends of the stack
    are taken as GaN reservoirs of the path's own well/reservoir mass.
    energy_eV may be a scalar or array-like; returns the matching shape."""
    if carrier not in ("electron", "hole"):
        raise ValueError("carrier must be 'electron' or 'hole'")
    if not math.isfinite(bias_V) or not math.isfinite(field_kVcm):
        raise ValueError("bias_V and field_kVcm must be finite")
    scalar_input = np.ndim(energy_eV) == 0
    energies = np.atleast_1d(np.asarray(energy_eV, dtype=float))
    if not np.all(np.isfinite(energies)):
        raise ValueError("energy_eV must be finite")
    out = np.empty_like(energies)
    for i, E in enumerate(energies):
        out[i], _ = _transmission_scalar(params, float(E), bias_V, field_kVcm, carrier)
    out = np.clip(out, 0.0, 1.0)   # MEDIUM 7 fix: physical clamp <= 1, not 1+1e-9;
                                    # the internal engine stays unclamped (see
                                    # rti_numerics_ok / _numerics_scan)
    return float(out[0]) if scalar_input else out


def reflection(params: NitrideNanowireInjectorParams, energy_eV, *,
               bias_V: float = 0.0, field_kVcm: float = 0.0,
               carrier: str = "electron"):
    """Companion to transmission(): R = 1 - T by flux conservation for this
    engine (see verify_nitride_nanowire_injector.py's T+R=1 unitarity
    check). Exposed for verification and diagnostics, not part of the
    minimal required interface."""
    if carrier not in ("electron", "hole"):
        raise ValueError("carrier must be 'electron' or 'hole'")
    if not math.isfinite(bias_V) or not math.isfinite(field_kVcm):
        raise ValueError("bias_V and field_kVcm must be finite")
    scalar_input = np.ndim(energy_eV) == 0
    energies = np.atleast_1d(np.asarray(energy_eV, dtype=float))
    if not np.all(np.isfinite(energies)):
        raise ValueError("energy_eV must be finite")
    out = np.empty_like(energies)
    for i, E in enumerate(energies):
        _, out[i] = _transmission_scalar(params, float(E), bias_V, field_kVcm, carrier)
    out = np.clip(out, 0.0, 1.0)   # MEDIUM 7 fix: physical clamp <= 1, not 1+1e-9
    return float(out[0]) if scalar_input else out


# ------------------------------------------------------------- rate/feasibility

def _fermi_dirac(E_eV: float, mu_eV: float, kT_eV: float) -> float:
    x = (E_eV - mu_eV) / kT_eV
    if x > 40.0:
        return math.exp(-x)
    if x < -40.0:
        return 1.0
    return 1.0 / (1.0 + math.exp(x))


def _with_thickness(params: NitrideNanowireInjectorParams, delta_nm: float) -> NitrideNanowireInjectorParams:
    """A params variant with both carrier paths' barrier thickness shifted
    by delta_nm (a single epitaxial run grows both barriers together),
    floored at one growth step so a perturbation cannot request a negative
    or sub-monolayer thickness."""
    floor_nm = params.growth_step_nm
    e_t = max(params.electron_barrier_thickness_nm + delta_nm, floor_nm)
    h_t = max(params.hole_barrier_thickness_nm + delta_nm, floor_nm)
    return replace(params, electron_barrier_thickness_nm=e_t, hole_barrier_thickness_nm=h_t)


def _analytic_well_seed_eV(m_w: float, m_b: float, V0: float, L_m: float):
    """Seed energy for the lowest quasi-bound resonance of a SYMMETRIC,
    mass-mismatched finite square well, from the standard finite-well
    transcendental equation (e.g. Griffiths, "Introduction to Quantum
    Mechanics", symmetric well, with the BenDaniel-Duke mass-mismatched
    matching condition already used throughout this module): even parity
    k*tan(k L/2) = (m_w/m_b)*kappa, k=sqrt(2 m_w E)/hbar,
    kappa=sqrt(2 m_b (V0-E))/hbar.  This is a LOCATOR only -- some designed
    double barriers give an extremely narrow true resonance (linewidth many
    orders of magnitude below the barrier height), which a blind uniform
    energy scan can straddle without ever sampling; the lowest root of this
    well-known equation is guaranteed to fall strictly before the tangent's
    first pole and is used only to seed the local transmission scan/refine
    below.  The reported resonance energy and width always come from the
    transfer-matrix transmission itself, never from this estimate.  Returns
    None if no root is found (falls back to a plain coarse scan)."""
    if not (math.isfinite(V0) and V0 > 0 and math.isfinite(L_m) and L_m > 0):
        return None

    def k_of(E):
        return math.sqrt(2.0 * m_w * M0_KG * E * EV_J) / HBAR_JS

    def kappa_of(E):
        return math.sqrt(2.0 * m_b * M0_KG * max(V0 - E, 1e-30) * EV_J) / HBAR_JS

    def g(E):
        k = k_of(E)
        return k * math.tan(k * L_m / 2.0) - (m_w / m_b) * kappa_of(E)

    n = 20000
    Es = np.linspace(V0 * 1e-6, V0 * (1.0 - 1e-6), n)
    prev = g(Es[0])
    for j in range(1, n):
        cur = g(Es[j])
        if math.isfinite(prev) and math.isfinite(cur) and prev < 0.0 <= cur:
            try:
                return float(brentq(g, Es[j - 1], Es[j]))
            except ValueError:
                return float(Es[j])
        prev = cur
    return None


def _tilted_faces_eV(params: NitrideNanowireInjectorParams, path: "_ResolvedPath",
                      bias_V: float, field_kVcm: float):
    """Per-segment (V_eV, entry_face_eV, exit_face_eV) for this path's
    stack in the ACTUAL (bias/field tilt + polarization) profile,
    referenced to the SAME (emitter) reservoir zero as transmission()'s
    energy_eV, plus the tilt (eV/m). Shared by _tilted_window_eV (resonance
    -search bounds) and _tilted_profile_max_eV (thermionic/over-barrier
    partition bound) so the two can never independently disagree about
    where the tilted profile's segment endpoints sit.  Each segment's
    potential is linear in x (tilt and the running polarization offset are
    both linear within one flat-band segment), so the segment's extremum
    is exactly one of its two endpoints -- the polarization contribution
    is accumulated continuously (never reset) exactly as in
    _transmission_reflection_scalar."""
    segs = _stack_segments(path)
    total_len_m = sum(s[0] for s in segs)
    tilt = _tilt_eV_per_m(field_kVcm, bias_V, params.field_leverarm, total_len_m)
    if params.include_polarization:
        pol_fields = _stack_polarization_fields_eV_per_m(params.al_fraction, segs, params.polarity)
    else:
        pol_fields = [0.0] * len(segs)
    x = 0.0
    pol_offset = 0.0
    faces = []
    for (length_m, V_eV, m_ratio), E_pol in zip(segs, pol_fields):
        x0, x1 = x, x + length_m
        pol0, pol1 = pol_offset, pol_offset + E_pol * length_m
        v0 = V_eV - tilt * x0 + pol0
        v1 = V_eV - tilt * x1 + pol1
        faces.append((V_eV, v0, v1))
        pol_offset = pol1
        x = x1
    return faces, tilt


def _tilted_window_eV(params: NitrideNanowireInjectorParams, path: "_ResolvedPath",
                       bias_V: float, field_kVcm: float):
    """(well_floor_eV, lower_barrier_top_eV, tilt_eV_per_m) bounding the
    energy window a resonance search must cover for this path's stack at
    the given bias/field, referenced to the SAME (emitter) reservoir zero
    as transmission()'s energy_eV.  A nonzero field/bias tilts the barrier
    stack (see _tilt_eV_per_m): the true transmission resonance can sit far
    from its flat-band location -- HIGH 2 fix; Opus review 2026-09-14 found
    the previous (flat-band-seeded, +/-10%-window) finder stuck at 62 meV /
    T 7e-7 for the electron path at 50 kV/cm, against a dense scan of this
    engine's own transmission() peaking at 252 meV / T 0.98.  The coarse
    scan below must therefore span the FULL field-tilted window: from the
    lowest point either barrier's tilted top reaches (past which the
    structure no longer presents a barrier at all) down to the lowest
    point the tilted well floor (or the reservoir zero) reaches.

    HIGH 1 / MEDIUM 11 fix: this window must be built from the ACTUAL
    (polarization + bias + field) profile the transmission engine uses --
    previously this function ignored the polarization tilt entirely, so a
    barrier "top" computed here could sit well above where the true
    (polarization-tilted) barrier top actually is, and _find_resonance
    could search a window with no relation to the real transmission
    landscape.

    HIGH 1 fix (Opus re-review of 12b39cd, 2026-09-14): a tilted barrier's
    own top is its HIGHEST face, max(v0, v1), not min(v0, v1) -- the
    previous min() picked the barrier's LOWEST face, which for the hole
    path collapsed the window to (well_floor, barrier_height - tilt) and
    made _find_resonance bail entirely ("no hole resonance"), when a
    200k-point scan of this engine's own transmission() finds three real
    hole resonances once the window is built correctly.  The window's
    upper bound stays the LOWER of the two (now correctly computed)
    barrier TOPS -- past that point at least one barrier no longer
    presents any barrier at all, so no bound state can exist above it."""
    faces, tilt = _tilted_faces_eV(params, path, bias_V, field_kVcm)
    barrier_tops = [max(v0, v1) for (V_eV, v0, v1) in faces if V_eV > 0.0]
    well_floors = [0.0] + [min(v0, v1) for (V_eV, v0, v1) in faces if V_eV <= 0.0]
    lower_barrier_top = min(barrier_tops) if barrier_tops else float("nan")
    well_floor = min(well_floors)
    return well_floor, lower_barrier_top, tilt


def _tilted_profile_max_eV(params: NitrideNanowireInjectorParams, path: "_ResolvedPath",
                            bias_V: float, field_kVcm: float) -> float:
    """HIGH 2 fix (Opus re-review of 12b39cd, 2026-09-14): the ACTUAL
    tilted-profile maximum a carrier must clear to be genuinely over-
    barrier/thermionic -- the MAXIMUM of max(v0, v1) across every barrier
    segment (not the lower of the two barrier tops _tilted_window_eV uses
    to bound the resonance-search window: that bound exists so the search
    does not run past the point where the structure stops being a barrier
    AT ALL, whereas the thermionic partition needs the highest point a
    carrier must actually clear).  Must replace the flat-band
    path.barrier_height_eV wherever the thermionic/over-barrier partition
    of the Landauer integral, the orbital/continuum separation, or that
    integral's upper limit is computed -- the polarization-tilted profile
    peaks well above the flat-band barrier height (e.g. the designed hole
    stack's flat-band barrier 0.210 eV vs. tilted top 0.505 eV), so using
    the flat-band value there flips the thermionic-bypass screen from a
    genuine ~40-95x fail to a spurious pass.  Falls back to the flat-band
    path.barrier_height_eV only if the path has no barrier segment at all
    (should not occur for single_barrier/double_barrier topologies)."""
    faces, _tilt = _tilted_faces_eV(params, path, bias_V, field_kVcm)
    barrier_tops = [max(v0, v1) for (V_eV, v0, v1) in faces if V_eV > 0.0]
    if not barrier_tops:
        return path.barrier_height_eV
    top = max(barrier_tops)
    return top if math.isfinite(top) else path.barrier_height_eV


def _find_resonance(params: NitrideNanowireInjectorParams, carrier: str,
                     bias_V: float, field_kVcm: float):
    """Locate the lowest double-barrier resonance under the given bias/
    field. HIGH 2 fix: scans the module's own transmission() on a coarse
    grid over the FULL field-tilted window (_tilted_window_eV above, never
    a +/-10% band around the flat-band analytic seed, which a field can
    move the true resonance well outside of), then refines with a
    multi-level local zoom on the transfer-matrix transmission itself, so
    the REPORTED energy and height are always numerically computed, never
    the analytic estimate.  The flat-band estimate is folded in only as one
    extra sample point injected into the coarse grid (some designed double
    barriers give an extremely narrow true resonance that a blind uniform
    grid can straddle without ever sampling); a coarse pass that finds
    nothing resolvable is retried once with a much denser full-window scan
    before giving up.  Returns (E_res_eV, T_res) or None (single_barrier
    topology, or no resolvable peak anywhere in the window)."""
    path = _resolve_path(params, carrier)
    if path.topology != "double_barrier":
        return None
    well_floor, lower_top, _tilt = _tilted_window_eV(params, path, bias_V, field_kVcm)
    if not (math.isfinite(lower_top) and lower_top > well_floor):
        return None

    def T_at(E):
        return _transmission_scalar(params, E, bias_V, field_kVcm, carrier)[0]

    e_lo = max(well_floor, 1e-9)
    e_hi_scan = lower_top - abs(lower_top) * 1e-6 - 1e-12
    if e_hi_scan <= e_lo:
        return None
    seed = _analytic_well_seed_eV(path.m_well, path.m_barrier,
                                   path.barrier_height_eV, path.well_nm * 1e-9)

    def _first_local_max(Ts, threshold):
        """Index of the FIRST (lowest-energy) sampled local maximum above
        threshold, scanning low-to-high -- a real double barrier generally
        supports several quasi-bound states (this design's flat-band
        window has three), and the lowest one is the physically relevant
        target; a plain global argmax would instead return whichever
        excited state happens to transmit best, silently mispricing
        alignment against the wrong state."""
        for j in range(1, len(Ts) - 1):
            if Ts[j] > threshold and Ts[j] >= Ts[j - 1] and Ts[j] >= Ts[j + 1]:
                return j
        return None

    def _coarse_scan(n, threshold):
        grid = np.linspace(e_lo, e_hi_scan, n)
        if seed is not None and e_lo < seed < e_hi_scan:
            grid = np.sort(np.append(grid, seed))
        Ts = np.array([T_at(E) for E in grid])
        i = _first_local_max(Ts, threshold)
        return grid, Ts, i

    grid, Ts, i = _coarse_scan(2000, 1e-9)
    if i is None:
        grid, Ts, i = _coarse_scan(20000, 1e-12)   # narrow-resonance safety net
    if i is None:
        return None

    best_E, best_T = float(grid[i]), float(Ts[i])
    lo = grid[max(i - 2, 0)]
    hi = grid[min(i + 2, len(grid) - 1)]
    for _ in range(4):
        if hi <= lo:
            break
        local = np.linspace(lo, hi, 60)
        Ts_local = np.array([T_at(E) for E in local])
        j = int(np.argmax(Ts_local))
        if Ts_local[j] > best_T:
            best_E, best_T = float(local[j]), float(Ts_local[j])
        span = (hi - lo) / 30.0
        lo = max(local[j] - span, e_lo)
        hi = min(local[j] + span, e_hi_scan)

    if hi <= lo:
        return best_E, best_T
    res = minimize_scalar(lambda E: -T_at(E), bounds=(lo, hi), method="bounded",
                           options={"xatol": 1e-12})
    E_res, T_res = float(res.x), float(T_at(res.x))
    if T_res < best_T:
        E_res, T_res = best_E, best_T
    return E_res, T_res


def _resonance_width_eV(params: NitrideNanowireInjectorParams, carrier: str,
                         bias_V: float, field_kVcm: float,
                         E_res: float, T_res: float, e_hi: float) -> float:
    """FWHM of the resonance found by _find_resonance, via half-max
    bracketing on each side. Returns NaN if it cannot be bracketed (missing
    information, never assumed zero or infinite)."""
    half = T_res / 2.0

    def f(E):
        return _transmission_scalar(params, E, bias_V, field_kVcm, carrier)[0] - half

    step = max(e_hi * 1e-4, 1e-6)
    lo = hi = None
    x = E_res
    for _ in range(200):
        x -= step
        if x <= 0:
            break
        if f(x) < 0:
            lo = brentq(f, x, x + step)
            break
    x = E_res
    for _ in range(200):
        x += step
        if x >= e_hi:
            break
        if f(x) < 0:
            hi = brentq(f, x - step, x)
            break
    if lo is None or hi is None:
        return float("nan")
    return float(hi - lo)


def _forward_rate_hz(params: NitrideNanowireInjectorParams, carrier: str,
                      bias_V: float, field_kVcm: float, E_center: float,
                      combined_width_eV: float, kT_eV: float, mu_eV: float):
    """Single-channel Landauer capture rate into an EMPTY target state,
    Gamma = g_s * channels * (1/h) INTEGRAL T(E) f_reservoir(E) dE (see
    module docstring) -- MEDIUM 8 fix: params.degeneracy (g_s, spin/valley)
    and the caller-supplied transverse-channel count are BOTH applied (the
    degeneracy factor was previously dead/unused).  `mu_eV` is the
    injecting reservoir's quasi-Fermi level in the carrier's own energy
    axis, passed in explicitly by the caller (electron: the degenerate
    free-electron-gas estimate; hole: the Mg-acceptor-neutrality quasi-
    Fermi level, MEDIUM 6) rather than re-derived from path.mu_eV here.
    Also returns the same integral restricted to energies above the
    ACTUAL tilted-profile maximum (HIGH 2 fix, Opus re-review of 12b39cd,
    2026-09-14: the flat-band path.barrier_height_eV understates the true
    peak by the full polarization tilt, e.g. 0.2951 eV per barrier for the
    designed stack, which alone flipped the hole thermionic-bypass screen
    from a genuine ~40-95x fail to a spurious pass -- see
    _tilted_profile_max_eV) -- the thermionic/over-barrier share, used for
    bypass pricing without double-counting: it is the SAME T(E) integral,
    partitioned by energy domain, not a second exp(-barrier/kT) estimate),
    and whether scipy's quad raised an integration-quality warning
    (MEDIUM 7).

    MEDIUM 4 fix (Opus re-review of 12b39cd, 2026-09-14): a non-finite
    E_center means the caller has NO resolved center to price a rate
    around (a double-barrier path whose resonance could not be found --
    resonance_unresolved) -- this must return NaN/NaN, never silently
    substitute mu_eV and integrate as if the resonance sat right at the
    reservoir's own quasi-Fermi level (a fabricated, not missing, rate)."""
    if not math.isfinite(E_center):
        return float("nan"), float("nan"), False
    path = _resolve_path(params, carrier)
    true_top = _tilted_profile_max_eV(params, path, bias_V, field_kVcm)
    mu = mu_eV
    width = combined_width_eV if math.isfinite(combined_width_eV) else 0.0
    center = E_center
    e_lo = max(1e-6, min(center, mu) - 15.0 * kT_eV - 5.0 * width)
    # HIGH 2 fix: e_hi is bounded by the ACTUAL tilted profile maximum, not
    # the flat-band barrier height, and so extends at least 15 kT above the
    # true top (not merely above the understated flat-band value).
    e_hi = max(true_top, center, mu) + 15.0 * kT_eV

    def integrand(E):
        T = _transmission_scalar(params, E, bias_V, field_kVcm, carrier)[0]
        return T * _fermi_dirac(E, mu, kT_eV)

    hint_points = sorted({p for p in (center, true_top) if e_lo < p < e_hi})
    quad_warned = False
    with warnings.catch_warnings(record=True) as caught:
        warnings.simplefilter("always", category=IntegrationWarning)
        total, _ = quad(integrand, e_lo, e_hi, points=hint_points or None,
                         limit=200, epsabs=1e-18, epsrel=1e-7)
        if true_top < e_hi:
            above, _ = quad(integrand, true_top, e_hi,
                             limit=200, epsabs=1e-18, epsrel=1e-7)
        else:
            above = 0.0
        quad_warned = any(issubclass(w.category, IntegrationWarning) for w in caught)
    # [A] Caller-provided transverse-channel count from the levels module;
    # MEDIUM 6 fix: the default is now one (spin lives in `degeneracy`
    # alone), not two -- old callers relying on the previous default get
    # an explicit sensitivity, never a silent double-count of spin.
    channels = (params.reservoir_state_count_e if carrier == "electron"
                else params.reservoir_state_count_h)
    rate_total = params.degeneracy * channels * total / H_EVS
    rate_above = params.degeneracy * channels * above / H_EVS
    return float(rate_total), float(rate_above), quad_warned


def _scan_required_bias_shift_meV(params: NitrideNanowireInjectorParams, carrier: str,
                                   field_kVcm: float, mu_emitter_eV: float):
    """MEDIUM 10 fix: the ENGINE's own answer for how much bias is needed to
    bring the resonance onto the emitter quasi-Fermi level, replacing the
    previous 1:1 unit-lever-arm assumption the engine itself falsifies (a
    10 meV drop can move the resonance by much less, or non-monotonically,
    per the applied field). Scans bias_V over the declared
    bias_tuning_range_meV (converted through field_leverarm), finds the
    resonance at each bias via the production _find_resonance, and reports
    the bias point (converted back to an equivalent meV shift) and the
    residual |E_res - mu_emitter| there closest to zero. Single-barrier
    paths have no resonance to move in this model; returns (NaN, NaN).
    Returns (required_bias_shift_meV, best_residual_meV)."""
    if params.bias_tuning_range_meV <= 0.0 or params.field_leverarm <= 0.0:
        return float("nan"), float("nan")
    path = _resolve_path(params, carrier)
    if path.topology != "double_barrier":
        return float("nan"), float("nan")
    v_max = params.bias_tuning_range_meV / 1000.0 / params.field_leverarm
    best_err_eV = float("inf")
    best_bias_V = float("nan")
    for bV in np.linspace(-v_max, v_max, 41):
        res = _find_resonance(params, carrier, float(bV), field_kVcm)
        if res is None:
            continue
        err = abs(res[0] - mu_emitter_eV)
        if err < best_err_eV:
            best_err_eV, best_bias_V = err, float(bV)
    if not math.isfinite(best_bias_V):
        return float("nan"), float("nan")
    return best_bias_V * params.field_leverarm * 1000.0, best_err_eV * 1000.0


def _carrier_screen(params: NitrideNanowireInjectorParams, carrier: str, *,
                     kT_eV: float, level_eV: float, mu_emitter_eV: float, spacing_meV: float,
                     second_pair_addition_meV: float, available_pair_rate_Hz: float,
                     loading_window_ns: float, gate_ns: float, rep_rate_hz: float, field_kVcm: float):
    """All diagnostics for one carrier path at params' current (possibly
    growth-perturbed) thickness. Held at zero series bias: the interface
    constraint exposes only field_kVcm to injector_feasibility, and
    params.field_leverarm lets a card fold an equivalent series-bias
    contribution into field_kVcm [A, documented simplification -- bias_V
    itself remains an independent transmission()/reflection() argument for
    direct callers]. `mu_emitter_eV` is this carrier's own emitter
    reservoir quasi-Fermi level (HIGH 2); `level_eV` (the dot level) is
    used ONLY for the informational well-to-dot drop, never for gating."""
    bias_V = 0.0
    path = _resolve_path(params, carrier)

    res = _find_resonance(params, carrier, bias_V, field_kVcm)
    resonance_unresolved = False
    if path.topology == "double_barrier":
        if res is not None:
            E_res, T_res = res
            _well_floor, lower_top, _tilt = _tilted_window_eV(params, path, bias_V, field_kVcm)
            e_hi_width = lower_top if math.isfinite(lower_top) else path.barrier_height_eV
            linewidth_eV = _resonance_width_eV(params, carrier, bias_V, field_kVcm,
                                                E_res, T_res, e_hi_width)
            # HIGH 2 fix: ALIGNMENT is versus the carrier's own EMITTER
            # quasi-Fermi level, not the dot level -- the dot level only
            # enters as the informational well-to-dot drop below.
            alignment_error_meV = abs(E_res - mu_emitter_eV) * 1000.0
            E_center = E_res
            linewidth_na = False
            well_to_dot_drop_meV = (E_res - level_eV) * 1000.0 if math.isfinite(level_eV) else float("nan")
        else:
            # HIGH 1 / MEDIUM 11 fix: a double-barrier path with NO
            # resolvable resonance under the ACTUAL (polarization + bias +
            # field) profile is MISSING information, not an aligned/
            # zero-error pass -- both the linewidth and the alignment
            # error are unknown (NaN), the combined width is unknown (NaN,
            # same as the "found but unresolved width" branch below), and
            # injector_feasibility must record
            # "<carrier>_resonance_unresolved" and force valid=False (the
            # critical-NaN rule), never silently fall into the
            # single-barrier "aligned, zero error" branch.
            linewidth_eV = float("nan")
            alignment_error_meV = float("nan")
            E_center = float("nan")
            linewidth_na = False
            resonance_unresolved = True
            well_to_dot_drop_meV = float("nan")
    else:
        # Single tunnel barrier: no claimed resonant well, so no resonance
        # to misalign against [A]; a NaN width here is NOT APPLICABLE (no
        # resonance exists to have a width), distinct from a double-barrier
        # resonance whose width genuinely could not be bracketed -- MEDIUM
        # 6/7 fix: this must not be double-counted into combined_width_eV
        # alongside the alignment uncertainty below, and must not trip the
        # "unresolved linewidth" invalidity check the other case does.
        linewidth_eV = float("nan")
        alignment_error_meV = 0.0
        E_center = level_eV
        linewidth_na = True
        well_to_dot_drop_meV = float("nan")   # no resonance state to report a drop from

    linewidth_meV = linewidth_eV * 1000.0 if math.isfinite(linewidth_eV) else float("nan")
    if linewidth_na:
        # No resonance -> no separate linewidth to add in quadrature; the
        # combined width is the alignment uncertainty alone (MEDIUM 7 fix,
        # previously this branch quadrature-added alignment_uncertainty
        # with itself, over-broadening by sqrt(2)).
        combined_width_eV = params.alignment_uncertainty_meV / 1000.0
    elif math.isfinite(linewidth_eV):
        combined_width_eV = math.sqrt(linewidth_eV ** 2 + (params.alignment_uncertainty_meV / 1000.0) ** 2)
    else:
        # A double-barrier resonance was found but its width could not be
        # bracketed (or no resonance was found at all, resonance_unresolved
        # above): missing information, not zero broadening (MEDIUM 6 fix)
        # -- propagate NaN instead of assuming an infinitely sharp
        # resonance; see the linewidth_unresolved/resonance_unresolved
        # invalidity checks in injector_feasibility.
        combined_width_eV = float("nan")

    # HIGH 3 fix: allowed-state ALIGNMENT error (delivery into the
    # REQUESTED level) is a distinct, separately-gated screen from the
    # unwanted-state SEPARATION priced into margin_kT/orbital_margin_kT
    # below -- a high-transmission resonance that simply is not where the
    # dot's level sits must not be reported as a passing delivery rate.
    _combined_width_for_gate_eV = combined_width_eV if math.isfinite(combined_width_eV) else 0.0
    align_threshold_meV = max(2.0 * kT_eV, _combined_width_for_gate_eV) * 1000.0
    if params.alignment_tunable:
        # MEDIUM 10 fix: the required bias shift is the ENGINE's own
        # answer (a genuine bias sweep of the production resonance
        # finder), not an assumed 1:1 lever arm.
        required_bias_shift_meV, best_residual_meV = _scan_required_bias_shift_meV(
            params, carrier, field_kVcm, mu_emitter_eV)
        alignment_ok = math.isfinite(best_residual_meV) and best_residual_meV <= align_threshold_meV
    else:
        required_bias_shift_meV = float("nan")
        alignment_ok = (not resonance_unresolved) and math.isfinite(alignment_error_meV) \
            and alignment_error_meV <= align_threshold_meV

    rate_hz, rate_above_hz, quad_warned = _forward_rate_hz(
        params, carrier, bias_V, field_kVcm, E_center, combined_width_eV, kT_eV, mu_emitter_eV)
    rate_below_hz = max(rate_hz - rate_above_hz, 0.0)
    bypassed = params.bypass_prefactor * rate_above_hz
    denom = rate_below_hz + bypassed
    # MEDIUM 4 fix companion: an unresolved rate (NaN, see _forward_rate_hz)
    # must not silently read back as bypass_fraction 0.0 (a fabricated
    # "no bypass" pass) via the denom>0 comparison's False-on-NaN behavior.
    if not math.isfinite(rate_hz):
        bypass_fraction = float("nan")
    else:
        bypass_fraction = float(bypassed / denom) if denom > 0 else 0.0

    # HIGH 2 fix: continuum onset uses the ACTUAL tilted profile maximum
    # (see _tilted_profile_max_eV), not the flat-band path.barrier_height_eV
    # -- the true tilted top sits above the flat-band value by the full
    # polarization tilt, so using the flat-band value here understated the
    # separation to the continuum.
    true_top = _tilted_profile_max_eV(params, path, bias_V, field_kVcm)
    orbital_candidates = [c for c in (
        spacing_meV / 1000.0 if math.isfinite(spacing_meV) else float("nan"),
        true_top - E_center if math.isfinite(E_center) else float("nan"),
    ) if math.isfinite(c) and c > -1e9]
    orbital_sep_eV = min(orbital_candidates) if orbital_candidates else float("nan")
    full_candidates = list(orbital_candidates)
    if math.isfinite(second_pair_addition_meV):
        full_candidates.append(second_pair_addition_meV / 1000.0)
    full_sep_eV = min(full_candidates) if full_candidates else float("nan")

    def _margin_kT(sep_eV):
        if not math.isfinite(sep_eV) or not math.isfinite(combined_width_eV):
            return float("nan")
        return max(sep_eV - combined_width_eV, 0.0) / kT_eV

    margin_kT = _margin_kT(full_sep_eV)
    orbital_margin_kT = _margin_kT(orbital_sep_eV)

    t_load_s = loading_window_ns * 1e-9
    missed_load_p = math.exp(-rate_hz * t_load_s)

    period_s = 1.0 / rep_rate_hz
    # [A] Reload is priced over the explicit counting gate, not implicitly
    # over the unused remainder of an electrical period.
    t_remaining_s = min(gate_ns * 1e-9, period_s)
    if math.isfinite(second_pair_addition_meV) and math.isfinite(available_pair_rate_Hz):
        # LOW 16 fix: suppress reload by the EXCESS separation beyond the
        # combined resonance/alignment width, floored at 0 -- a charging
        # energy smaller than the state's own width cannot suppress reload
        # at all (exp(-x/kT) with x<0 would be > 1, unphysical).
        width_for_suppression_eV = combined_width_eV if math.isfinite(combined_width_eV) else 0.0
        net_eV = max(second_pair_addition_meV / 1000.0 - width_for_suppression_eV, 0.0)
        suppression = math.exp(-net_eV / kT_eV)
        rate_second_hz = min(rate_hz, available_pair_rate_Hz) * suppression
        second_pair_p = 1.0 - math.exp(-rate_second_hz * t_remaining_s)
    else:
        rate_second_hz = float("nan")
        second_pair_p = float("nan")

    return dict(
        rate_hz=rate_hz, alignment_error_meV=alignment_error_meV,
        linewidth_meV=linewidth_meV, linewidth_na=linewidth_na,
        resonance_unresolved=resonance_unresolved,
        alignment_ok=alignment_ok, required_bias_shift_meV=required_bias_shift_meV,
        bypass_fraction=bypass_fraction,
        margin_kT=margin_kT, orbital_margin_kT=orbital_margin_kT,
        missed_load_p=missed_load_p, second_pair_p=second_pair_p,
        rate_second_hz=rate_second_hz, barrier_height_eV=path.barrier_height_eV,
        thickness_nm=path.thickness_nm, well_to_dot_drop_meV=well_to_dot_drop_meV,
        quad_warned=quad_warned, second_pair_addition_meV=second_pair_addition_meV,
        E_center=E_center,   # LOW 10: internal-only, lets rti_numerics_ok refine
                              # its scan around the actual found resonance
    )


def _transport_ok(screen: dict, margin_floor_kT: float, bypass_ceiling: float,
                   missed_ceiling: float) -> bool:
    # HIGH 3 fix: a path whose delivery is not aligned onto the requested
    # level (screen["alignment_ok"] False) is not transport-feasible, even
    # if its resonant transmission/rate numbers look otherwise excellent --
    # resonant transmission is not the same claim as one-pair loading into
    # the intended state.
    return (not screen["resonance_unresolved"]
            and math.isfinite(screen["margin_kT"]) and screen["margin_kT"] >= margin_floor_kT
            and math.isfinite(screen["bypass_fraction"]) and screen["bypass_fraction"] <= bypass_ceiling
            and math.isfinite(screen["missed_load_p"]) and screen["missed_load_p"] <= missed_ceiling
            and screen["alignment_ok"])


def _nan_safe_min(a: float, b: float) -> float:
    """min(a, b), but NaN if EITHER is NaN -- plain min()/max() silently
    drop a NaN that appears as the SECOND argument (Python: min(nan, 5) is
    nan, but min(5, nan) is 5), which would hide a genuinely missing/
    unresolved carrier path behind the other, finite, path's number."""
    if not math.isfinite(a) or not math.isfinite(b):
        return float("nan")
    return min(a, b)


def _nan_safe_max(a: float, b: float) -> float:
    if not math.isfinite(a) or not math.isfinite(b):
        return float("nan")
    return max(a, b)


def _barrier_polarization_tilt_eV(params: NitrideNanowireInjectorParams, carrier: str) -> float:
    """Signed potential drop across ONE (the first) barrier of this
    carrier's stack, from the zero-net-drop closure (HIGH 3/4/5 fix) --
    replaces the old abs()-based fixed-D single value.  0.0 if
    include_polarization is False or the path has no barrier segment."""
    if not params.include_polarization:
        return 0.0
    path = _resolve_path(params, carrier)
    segs = _stack_segments(path)
    if not segs:
        return 0.0
    fields = _stack_polarization_fields_eV_per_m(params.al_fraction, segs, params.polarity)
    return float(fields[0] * segs[0][0])


def _numerics_scan_ok(params: NitrideNanowireInjectorParams, carrier: str,
                       field_kVcm: float, mu_eV: float, width_eV: float, kT_eV: float,
                       E_center: float = float("nan")) -> bool:
    """MEDIUM 7 fix: rti_numerics_ok is a REAL, falsifiable computation --
    max |T+R-1| <= 1e-8 and no RAW (unclamped) T > 1+1e-9 over a
    representative energy grid spanning this carrier's rate-integration
    window (same bounds logic as _forward_rate_hz).  Uses _transmission_
    scalar directly (bypassing the public API's [0,1] clamp) so genuine
    numerical overshoot from e.g. a deliberately coarse staircase
    (params.slice_length_nm / min_slices_per_segment) is visible instead
    of hidden by the clamp.

    LOW 10 fix (Opus re-review of 12b39cd, 2026-09-14): the broad 300-point
    scan above alone has ~3 meV spacing over a window that can be a volt or
    more wide, and so cannot land on a resonance whose own FWHM is a few
    meV or less -- a genuine numerical overshoot confined to the resonance
    itself would be invisible to it.  When the caller supplies a resolved
    resonance center (E_center finite, double-barrier path only), fold in
    a SECOND, densely-sampled scan of 50 points spanning +/- 5 FWHM
    (width_eV) around that center, so the diagnostic actually samples the
    one place a transfer-matrix resonance is most likely to misbehave."""
    path = _resolve_path(params, carrier)
    w = width_eV if math.isfinite(width_eV) else 0.0
    e_lo = max(1e-6, mu_eV - 15.0 * kT_eV - 5.0 * w)
    e_hi = max(path.barrier_height_eV, mu_eV) + 15.0 * kT_eV
    grids = []
    if e_hi > e_lo and math.isfinite(e_lo) and math.isfinite(e_hi):
        grids.append(np.linspace(e_lo, e_hi, 300))
    if math.isfinite(E_center) and w > 0.0:
        r_lo = max(E_center - 5.0 * w, 1e-9)
        r_hi = E_center + 5.0 * w
        if r_hi > r_lo:
            grids.append(np.linspace(r_lo, r_hi, 50))
    if not grids:
        return True
    max_resid = 0.0
    max_T = 0.0
    for grid in grids:
        for E in grid:
            T, R = _transmission_scalar(params, float(E), 0.0, field_kVcm, carrier)
            if not (math.isfinite(T) and math.isfinite(R)):
                return False
            max_resid = max(max_resid, abs(T + R - 1.0))
            max_T = max(max_T, T)
    return bool(max_resid <= 1e-8 and max_T <= 1.0 + 1e-9)


# --------------------------------------------------------- engineering thresholds
# [A] Proposed conditional engineering thresholds (design-brief round 3).

RTI_MARGIN_FLOOR_KT = 10.0
RTI_BYPASS_CEILING = 0.01
RTI_MISSED_LOAD_CEILING = 0.01
RTI_SECOND_PAIR_CEILING = 0.01

_PROVENANCE = (
    "[V] material endpoints (GaN/AlN masses, bandgaps, lattice constant) from "
    "fsim_core.nitride_materials (Bernardini PRB 1997; Bernardini & "
    "Fiorentini pss(b) 1999; Rinke PRB 2008; Wu et al. bandgap, APL 2002); "
    "[A] linear Al_xGa_1-xN VCA (no bowing) of that gap and of masses, "
    "unstrained barrier; [V] default valence-band-offset partition Martin, "
    "Yu, Waldrop, APL 68, 2541 (1996), 0.70 +/- 0.24 eV (Rinke PRB 2008 "
    "~0.8 eV), with nitride_materials' own pinned Tsai & Bayram ACS Omega "
    "2020 0.30 eV VBO reported as the alternative partition "
    "(rti_bypass_fraction_tsai_partition) -- [A] which partition is correct "
    "for this device; [V] n/p-GaN doping anchors Deshpande et al., Nat. "
    "Commun. 4, 1675 (2013); [DR] degenerate free-electron-gas reservoir "
    "electrochemical energy (Ashcroft & Mermin, 'Solid State Physics', "
    "1976, Ch. 2 Eq. 2.33); [V] single-channel Landauer/Buttiker rate "
    "normalization 1/h (e.g. Datta, 'Electronic Transport in Mesoscopic "
    "Systems', 1995, Ch. 2); [A] spin/valley degeneracy, bypass prefactor, "
    "field lever-arm, alignment uncertainty, alignment-gate threshold "
    "max(2 kT, combined width) unless alignment_tunable, growth tolerance "
    "in bilayer steps, engineering thresholds (margin 10 kT, "
    "bypass/missed-load/second-pair 0.01). Conditional engineering "
    "screening only; not a demonstrated hardware claim. [DR] barrier "
    "polarization uses B97 P_sp/e31/e33 VCA pseudomorphic on GaN and "
    "fixed-D interface sheets; [E] Mg activation 170 meV (Gotz APL 1996; "
    "Kozodoy JAP 2000), so rti_p_free_cm3 is not nominal Mg doping; [A] "
    "reservoir_state_count_e/h defaults to 1 (spin lives in degeneracy "
    "alone) until supplied by levels and "
    "gate_ns defaults to loading_window_ns for legacy pulse callers."
)


def _unknown_result(failed_checks, status: str) -> dict:
    """MEDIUM 12 fix: returns EVERY key the valid path (injector_feasibility's
    final return dict) returns, NaN/None/False-filled -- a caller destructuring
    any of those keys off an invalid-input result must not KeyError."""
    nan = float("nan")
    return dict(
        rti_feasible=False, rti_status=status, rti_transport_feasible=False,
        rti_level_margin_kT=nan, rti_orbital_margin_kT=nan,
        rti_alignment_error_meV=nan, rti_linewidth_meV=nan,
        rti_rate_Hz=nan, rti_e_rate_Hz=nan, rti_h_rate_Hz=nan,
        rti_bypass_fraction=nan, rti_missed_load_probability=nan,
        rti_second_pair_probability=nan, rti_second_carrier_probability=nan,
        rti_growth_feasible=False,
        rti_failed_checks=list(failed_checks),
        rti_evidence_status="conditional_engineering_screen_no_measured_device",
        valid=False, provenance=_PROVENANCE,
        rti_level_margin_e_kT=nan, rti_level_margin_h_kT=nan,
        rti_alignment_error_e_meV=nan, rti_alignment_error_h_meV=nan,
        rti_required_bias_shift_meV=nan, rti_bypass_fraction_tsai_partition=nan,
        rti_growth_nearest_commensurate_nm=None, rti_growth_perturbed_margins_kT=None,
        rti_barrier_polarization_tilt_eV={"electron": nan, "hole": nan},
        rti_well_to_dot_drop_meV={"electron": nan, "hole": nan},
        rti_p_free_cm3=nan, rti_reservoir_state_count_e=nan, rti_reservoir_state_count_h=nan,
        rti_gate_ns=nan, rti_numerics_ok=False, rti_polarity=None,
    )


def injector_feasibility(params: NitrideNanowireInjectorParams, *, T_K: float,
                          rep_rate_hz: float, loading_window_ns: float,
                          electron_level_eV: float, hole_level_eV: float,
                          electron_spacing_meV: float, hole_spacing_meV: float,
                          second_pair_addition_meV: float,
                          available_pair_rate_Hz: float,
                          field_kVcm: float = 0.0,
                          gate_ns: float | None = None) -> dict:
    """Price this injector design as a SECOND, independent deterministic-
    pair feasibility screen (see module docstring). Every energy argument
    is referenced to its own carrier's bulk GaN band edge (see
    'Energy zero' above). Returns a dict; NaN/missing critical diagnostics
    cannot pass (rti_feasible is False whenever a required numeric input is
    missing, non-finite, or out of range) -- this function never raises on
    bad SCREENING inputs, it fails the screen instead ("fail safely").

    occupancy_control_known / second_pair_control_known on params are
    necessary, not sufficient: rti_feasible still requires the computed
    rti_second_pair_probability (from actual second_pair_addition_meV and
    available_pair_rate_Hz numbers) to clear its threshold. There is no
    branch that can set rti_feasible=True without a passing numeric
    computation behind it.

    Each carrier path is also gated on ALIGNMENT: a resonance (or, for a
    single tunnel barrier, the requested level itself) more than
    max(2 kT, combined linewidth) away from the requested
    electron_level_eV/hole_level_eV fails that path's transport screen
    (params.alignment_tunable=True substitutes a bias-tuning-range gate on
    the reported rti_required_bias_shift_meV instead) -- resonant
    transmission is not itself evidence of delivery into the REQUESTED
    state. rti_second_pair_probability is P(second electron OR second
    hole), a conservative upper bound on the true second-PAIR probability;
    rti_second_carrier_probability is the same number under its honest
    name. rti_bypass_fraction_tsai_partition reports the same design's
    bypass fraction under the alternative (Tsai & Bayram 0.30 eV) valence
    partition, non-gating."""
    if gate_ns is None:
        gate_ns = loading_window_ns  # [A] compatibility: pulse width is its own gate
    raw = dict(T_K=T_K, rep_rate_hz=rep_rate_hz, loading_window_ns=loading_window_ns, gate_ns=gate_ns,
               electron_level_eV=electron_level_eV, hole_level_eV=hole_level_eV,
               electron_spacing_meV=electron_spacing_meV, hole_spacing_meV=hole_spacing_meV,
               second_pair_addition_meV=second_pair_addition_meV,
               available_pair_rate_Hz=available_pair_rate_Hz, field_kVcm=field_kVcm)

    def _nonfinite(x):
        try:
            return x is None or not math.isfinite(x)
        except TypeError:
            return True

    # Orbital/continuum selectivity tolerates a missing excited-state
    # spacing, and second-pair pricing tolerates a missing addition energy
    # or supply rate (all treated as unknown, not infinite, downstream via
    # the "second_pair_exclusion_unknown" / orbital-only diagnostics) --
    # those are NOT hard requirements to run the screen at all.  A missing
    # hole (or electron) target level, temperature, repetition rate, or
    # loading window IS: there is no energy to center the rate integral on.
    unknown_ok = ("electron_spacing_meV", "hole_spacing_meV",
                  "second_pair_addition_meV", "available_pair_rate_Hz")
    required = {k: v for k, v in raw.items() if k not in unknown_ok}
    missing = [k for k, v in required.items() if _nonfinite(v)]
    range_bad = []
    if "T_K" not in missing and T_K <= 0:
        range_bad.append("out_of_range:T_K")
    if "rep_rate_hz" not in missing and rep_rate_hz <= 0:
        range_bad.append("out_of_range:rep_rate_hz")
    if "loading_window_ns" not in missing and loading_window_ns <= 0:
        range_bad.append("out_of_range:loading_window_ns")
    if "gate_ns" not in missing and gate_ns <= 0:
        range_bad.append("out_of_range:gate_ns")
    if not _nonfinite(available_pair_rate_Hz) and available_pair_rate_Hz < 0:
        range_bad.append("out_of_range:available_pair_rate_Hz")
    if ("loading_window_ns" not in missing and "rep_rate_hz" not in missing
            and loading_window_ns * 1e-9 > 1.0 / rep_rate_hz):
        range_bad.append("out_of_range:loading_window_exceeds_period")

    failed = [f"missing:{k}" for k in missing] + range_bad
    if failed:
        return _unknown_result(failed, "missing_or_invalid_input")

    kT_eV = KB_EV * T_K
    # MEDIUM 6 fix: the ionised free-hole density (and the reservoir's own
    # quasi-Fermi level) come from Mg-acceptor charge neutrality, not the
    # earlier Na*exp(-E_A/kT) mislabelled "ionisation fraction" -- the
    # number is USED for mu_h below, not just reported.
    p_free_cm3, mu_h_eV = _hole_quasi_fermi_eV(params, T_K)

    def _screen_at(p):
        mu_e = _resolve_path(p, "electron").mu_eV
        _p_free, mu_h = _hole_quasi_fermi_eV(p, T_K)
        e = _carrier_screen(p, "electron", kT_eV=kT_eV, level_eV=electron_level_eV,
                             mu_emitter_eV=mu_e, spacing_meV=electron_spacing_meV,
                             second_pair_addition_meV=second_pair_addition_meV,
                             available_pair_rate_Hz=available_pair_rate_Hz,
                             loading_window_ns=loading_window_ns, gate_ns=gate_ns, rep_rate_hz=rep_rate_hz,
                             field_kVcm=field_kVcm)
        h = _carrier_screen(p, "hole", kT_eV=kT_eV, level_eV=hole_level_eV,
                             mu_emitter_eV=mu_h, spacing_meV=hole_spacing_meV,
                             second_pair_addition_meV=second_pair_addition_meV,
                             available_pair_rate_Hz=available_pair_rate_Hz,
                             loading_window_ns=loading_window_ns, gate_ns=gate_ns, rep_rate_hz=rep_rate_hz,
                             field_kVcm=field_kVcm)
        return e, h

    e, h = _screen_at(params)

    failed_checks = []
    e_ok = _transport_ok(e, RTI_MARGIN_FLOOR_KT, RTI_BYPASS_CEILING, RTI_MISSED_LOAD_CEILING)
    h_ok = _transport_ok(h, RTI_MARGIN_FLOOR_KT, RTI_BYPASS_CEILING, RTI_MISSED_LOAD_CEILING)
    if not e_ok:
        failed_checks.append("electron_transport")
        if e["resonance_unresolved"]:
            failed_checks.append("electron_resonance_unresolved")
        elif not e["alignment_ok"]:
            failed_checks.append("electron_alignment")
    if not h_ok:
        failed_checks.append("hole_transport")
        if h["resonance_unresolved"]:
            failed_checks.append("hole_resonance_unresolved")
        elif not h["alignment_ok"]:
            failed_checks.append("hole_alignment")

    # MEDIUM 6 fix: a double-barrier resonance that was FOUND but whose
    # width could not be bracketed is missing information, not zero
    # broadening -- it must invalidate the screen. A single-barrier path's
    # NaN linewidth is "not applicable" (no resonance exists at all) and
    # must NOT trip this (see _carrier_screen's linewidth_na). A path whose
    # resonance was never found at all is tagged "<carrier>_
    # resonance_unresolved" above instead (a more specific diagnosis of the
    # same missing-information invalidity), so it is excluded here to avoid
    # a redundant/confusing double tag.
    e_linewidth_unresolved = ((not e["linewidth_na"]) and not e["resonance_unresolved"]
                               and not math.isfinite(e["linewidth_meV"]))
    h_linewidth_unresolved = ((not h["linewidth_na"]) and not h["resonance_unresolved"]
                               and not math.isfinite(h["linewidth_meV"]))
    if e_linewidth_unresolved:
        failed_checks.append("electron_linewidth_unresolved")
    if h_linewidth_unresolved:
        failed_checks.append("hole_linewidth_unresolved")
    linewidth_unresolved = (e_linewidth_unresolved or h_linewidth_unresolved
                             or e["resonance_unresolved"] or h["resonance_unresolved"])

    missed_load_probability = 1.0 - (1.0 - e["missed_load_p"]) * (1.0 - h["missed_load_p"])
    if missed_load_probability > RTI_MISSED_LOAD_CEILING:
        failed_checks.append("pair_missed_load_window")

    rti_transport_feasible = bool(
        e_ok and h_ok and missed_load_probability <= RTI_MISSED_LOAD_CEILING
        and not linewidth_unresolved
    )

    # Growth realism (HIGH 1 fix): growth_tolerance_steps is a TOLERANCE
    # around the nearest commensurate (integer-bilayer) thickness, not a
    # commensurability test in itself -- a nominal thickness is
    # growth-feasible only if (a) it lies within
    # growth_tolerance_steps*growth_step_nm of its nearest integer-bilayer
    # thickness AND (b) the feasibility conjunction still holds when the
    # barrier stack is perturbed by +/- growth_tolerance_steps bilayers in
    # BOTH directions (a design that only clears exactly at nominal is
    # fragile, never a manufacturable claim).
    def _nearest_commensurate_nm(thickness_nm):
        n = round(thickness_nm / params.growth_step_nm)
        return n * params.growth_step_nm

    nearest_e_nm = _nearest_commensurate_nm(params.electron_barrier_thickness_nm)
    nearest_h_nm = _nearest_commensurate_nm(params.hole_barrier_thickness_nm)
    tolerance_nm = params.growth_tolerance_steps * params.growth_step_nm
    commensurate_ok = (abs(params.electron_barrier_thickness_nm - nearest_e_nm) <= tolerance_nm
                        and abs(params.hole_barrier_thickness_nm - nearest_h_nm) <= tolerance_nm)

    # MEDIUM 9 fix: the growth_tolerance TAG is a diagnosis of the
    # growth-specific criteria alone (commensurability + perturbed-
    # thickness robustness), decoupled from whether NOMINAL transport
    # already failed for an unrelated reason (bypass, margin, alignment,
    # ...) -- previously growth_ok = transport_feasible and commensurate_ok
    # meant "growth_tolerance" was named whenever transport failed at all,
    # misattributing the cause. rti_growth_feasible (the boolean folded
    # into the overall rti_feasible AND-gate) still requires nominal
    # transport to pass too.
    #
    # MEDIUM 5 fix (Opus re-review of 12b39cd, 2026-09-14): a perturbed
    # screen re-running the FULL _transport_ok conjunction reproduces any
    # nominal failure regardless of thickness (bypass/margin/alignment do
    # not improve just because the barrier moved by one growth step), so
    # ANDing the perturbed pass/fail into growth_examined_ok tagged
    # "growth_tolerance" on every one of those unrelated nominal failures
    # too.  The tag must instead mean a growth-specific REGRESSION: a
    # per-check pass at nominal that becomes a fail under perturbation.
    # A check that already failed at nominal is never re-counted here (it
    # is already reported under its own nominal failure tag above).
    missed_load_ok_nominal = missed_load_probability <= RTI_MISSED_LOAD_CEILING
    growth_regression = False
    perturbed_margins_kT = {}
    if params.growth_tolerance_steps > 0:
        for sign, tag in ((-1.0, "minus"), (1.0, "plus")):
            delta = sign * params.growth_tolerance_steps * params.growth_step_nm
            p_pert = _with_thickness(params, delta)
            e_p, h_p = _screen_at(p_pert)
            e_p_ok = _transport_ok(e_p, RTI_MARGIN_FLOOR_KT, RTI_BYPASS_CEILING, RTI_MISSED_LOAD_CEILING)
            h_p_ok = _transport_ok(h_p, RTI_MARGIN_FLOOR_KT, RTI_BYPASS_CEILING, RTI_MISSED_LOAD_CEILING)
            missed_pert = 1.0 - (1.0 - e_p["missed_load_p"]) * (1.0 - h_p["missed_load_p"])
            missed_pert_ok = missed_pert <= RTI_MISSED_LOAD_CEILING
            perturbed_margins_kT[tag] = dict(
                electron_margin_kT=e_p["margin_kT"], hole_margin_kT=h_p["margin_kT"],
                missed_load_probability=missed_pert, ok=bool(e_p_ok and h_p_ok and missed_pert_ok),
            )
            if (e_ok and not e_p_ok) or (h_ok and not h_p_ok) or (missed_load_ok_nominal and not missed_pert_ok):
                growth_regression = True
    growth_examined_ok = commensurate_ok and not growth_regression
    rti_growth_feasible = bool(rti_transport_feasible and growth_examined_ok)
    if not growth_examined_ok:
        failed_checks.append("growth_tolerance")

    second_pair_known = bool(math.isfinite(e["second_pair_p"]) and math.isfinite(h["second_pair_p"]))
    if not second_pair_known:
        failed_checks.append("second_pair_exclusion_unknown")
        second_pair_probability = float("nan")
    else:
        # P(second electron OR second hole): a conservative UPPER BOUND on
        # true second-PAIR probability, not the pair probability itself
        # (LOW 10 fix -- rti_second_carrier_probability is the honestly
        # named quantity; rti_second_pair_probability is retained,
        # documented as that same conservative bound, for contract
        # compatibility).
        second_pair_probability = 1.0 - (1.0 - e["second_pair_p"]) * (1.0 - h["second_pair_p"])
        if second_pair_probability > RTI_SECOND_PAIR_CEILING:
            # LOW 16 fix: an honest, distinguishing tag when the charging
            # energy itself is smaller than kT (no floored suppression can
            # help; this is a fundamental design limit at this
            # temperature), vs. the plain tag when E_C > kT but the
            # combined-width-subtracted margin still isn't enough.
            if math.isfinite(second_pair_addition_meV) and second_pair_addition_meV / 1000.0 < kT_eV:
                failed_checks.append("second_pair_reload (E_C below kT)")
            else:
                failed_checks.append("second_pair_reload")

    if not params.occupancy_control_known:
        failed_checks.append("occupancy_control_unspecified")
    if not params.second_pair_control_known:
        failed_checks.append("second_pair_control_unspecified")

    rti_feasible = bool(
        rti_transport_feasible
        and rti_growth_feasible
        and second_pair_known
        and second_pair_probability <= RTI_SECOND_PAIR_CEILING
        and params.occupancy_control_known
        and params.second_pair_control_known
    )

    if rti_feasible:
        rti_status = "feasible"
    elif not (params.occupancy_control_known and params.second_pair_control_known
              and second_pair_known):
        rti_status = "unknown_incomplete"
    else:
        rti_status = "infeasible"

    rti_level_margin_kT = _nan_safe_min(e["margin_kT"], h["margin_kT"])
    rti_orbital_margin_kT = _nan_safe_min(e["orbital_margin_kT"], h["orbital_margin_kT"])
    rti_alignment_error_meV = _nan_safe_max(e["alignment_error_meV"], h["alignment_error_meV"])
    rti_linewidth_meV = _nan_safe_max(e["linewidth_meV"], h["linewidth_meV"])
    rti_bypass_fraction = _nan_safe_max(e["bypass_fraction"], h["bypass_fraction"])
    rti_rate_Hz = _nan_safe_min(e["rate_hz"], h["rate_hz"])   # slower completed-pair rate
    if params.alignment_tunable:
        _bias_shift_candidates = [v for v in (e["required_bias_shift_meV"], h["required_bias_shift_meV"])
                                   if math.isfinite(v)]
        rti_required_bias_shift_meV = max(_bias_shift_candidates) if _bias_shift_candidates else float("nan")
    else:
        rti_required_bias_shift_meV = float("nan")

    # HIGH 4 sensitivity: the same designed stack's bypass fraction under
    # the Tsai & Bayram 0.30 eV valence-offset partition (nitride_materials'
    # own pinned VBO) instead of this module's default 0.70 eV partition --
    # reported, never gated on, so the choice of partition is visible
    # without silently changing the pass/fail verdict above.
    p_tsai = replace(params, delta_Ev_GaN_AlN_eV=_TSAI_BAYRAM_DELTA_EV_EV)
    e_tsai, h_tsai = _screen_at(p_tsai)
    rti_bypass_fraction_tsai_partition = float(max(e_tsai["bypass_fraction"], h_tsai["bypass_fraction"]))

    critical = [rti_level_margin_kT, rti_alignment_error_meV, e["rate_hz"], h["rate_hz"],
                rti_bypass_fraction, missed_load_probability]
    valid = bool(all(math.isfinite(v) for v in critical) and not linewidth_unresolved)
    if not valid:
        rti_feasible = False
        if rti_status == "feasible":
            rti_status = "infeasible"

    # MEDIUM 7 fix: a real, falsifiable numerical-quality diagnostic (see
    # _numerics_scan_ok) instead of a hardcoded True -- also folds in
    # whether scipy's quad raised an integration-quality warning while
    # computing either carrier's rate.
    mu_e = _resolve_path(params, "electron").mu_eV
    e_width_for_scan = e["linewidth_meV"] / 1000.0 if math.isfinite(e["linewidth_meV"]) else params.alignment_uncertainty_meV / 1000.0
    h_width_for_scan = h["linewidth_meV"] / 1000.0 if math.isfinite(h["linewidth_meV"]) else params.alignment_uncertainty_meV / 1000.0
    rti_numerics_ok = bool(
        _numerics_scan_ok(params, "electron", field_kVcm, mu_e, e_width_for_scan, kT_eV,
                           E_center=e["E_center"])
        and _numerics_scan_ok(params, "hole", field_kVcm, mu_h_eV, h_width_for_scan, kT_eV,
                               E_center=h["E_center"])
        and not e["quad_warned"] and not h["quad_warned"]
    )

    return dict(
        rti_feasible=rti_feasible,
        rti_status=rti_status,
        rti_transport_feasible=rti_transport_feasible,
        rti_level_margin_kT=float(rti_level_margin_kT),
        rti_orbital_margin_kT=float(rti_orbital_margin_kT),
        rti_alignment_error_meV=float(rti_alignment_error_meV),
        rti_linewidth_meV=float(rti_linewidth_meV),
        rti_rate_Hz=float(rti_rate_Hz),
        rti_e_rate_Hz=float(e["rate_hz"]),
        rti_h_rate_Hz=float(h["rate_hz"]),
        rti_bypass_fraction=float(rti_bypass_fraction),
        rti_missed_load_probability=float(missed_load_probability),
        rti_second_pair_probability=float(second_pair_probability),
        rti_second_carrier_probability=float(second_pair_probability),
        rti_growth_feasible=rti_growth_feasible,
        rti_failed_checks=failed_checks,
        rti_evidence_status="conditional_engineering_screen_no_measured_device",
        valid=valid,
        provenance=_PROVENANCE,
        # extra, non-required diagnostics for transparency (device piece may ignore)
        rti_level_margin_e_kT=float(e["margin_kT"]),
        rti_level_margin_h_kT=float(h["margin_kT"]),
        rti_alignment_error_e_meV=float(e["alignment_error_meV"]),
        rti_alignment_error_h_meV=float(h["alignment_error_meV"]),
        rti_required_bias_shift_meV=float(rti_required_bias_shift_meV),
        rti_bypass_fraction_tsai_partition=rti_bypass_fraction_tsai_partition,
        rti_growth_nearest_commensurate_nm={"electron": float(nearest_e_nm), "hole": float(nearest_h_nm)},
        rti_growth_perturbed_margins_kT=perturbed_margins_kT,
        rti_barrier_polarization_tilt_eV={
            "electron": _barrier_polarization_tilt_eV(params, "electron"),
            "hole": _barrier_polarization_tilt_eV(params, "hole"),
        },
        rti_well_to_dot_drop_meV={
            "electron": float(e["well_to_dot_drop_meV"]),
            "hole": float(h["well_to_dot_drop_meV"]),
        },
        rti_p_free_cm3=float(p_free_cm3),
        rti_reservoir_state_count_e=float(params.reservoir_state_count_e),
        rti_reservoir_state_count_h=float(params.reservoir_state_count_h),
        rti_gate_ns=float(gate_ns),
        rti_numerics_ok=rti_numerics_ok,
        rti_polarity=params.polarity,
    )
