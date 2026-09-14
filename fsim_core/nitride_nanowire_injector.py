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
* Barrier material: an Al_xGa_1-xN barrier's conduction/valence offsets and
  effective masses come from a linear (Vegard-law) virtual-crystal
  interpolation between the GaN and AlN binaries already parameterized and
  cited in fsim_core.nitride_materials (Bernardini PRB 1997, Bernardini &
  Fiorentini pss(b) 1999, Rinke PRB 2008, Tsai & Bayram ACS Omega 2020) --
  [A] interpolation of [V] endpoints, exactly the same epistemic split that
  module already uses for its In_xGa_1-xN alloy.  Barrier bowing is
  neglected [A].  The barrier is evaluated UNSTRAINED (strain_fraction=0):
  this piece does not model the barrier as pseudomorphically strained to
  the GaN reservoir lattice; that is a conservative simplifying choice
  (coherency strain would shift the offset by an amount this piece does not
  attempt to sign) and is recorded under the module's `provenance` string
  and the worker's STATUS `decisions:` field.
* Growth realism: barrier thickness is checked against the GaN c-axis
  bilayer spacing (c/2, taken from nitride_materials' [V] lattice constant)
  as the growth increment; a design that only clears its screens at an
  exact nominal thickness and fails at +/- one growth step is reported
  fragile (rti_growth_feasible=False), never claimed manufacturable.

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
from collections import namedtuple
from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, minimize_scalar

from . import nitride_materials as NM

# ------------------------------------------------------------------ constants

HBAR_JS = 1.054571817e-34      # [V] CODATA 2018
H_EVS = 4.135667696e-15        # [V] CODATA 2018, h in eV s
KB_EV = 8.617333262e-5         # [V] CODATA 2018, k_B in eV/K
M0_KG = 9.1093837015e-31       # [V] CODATA 2018
EV_J = 1.602176634e-19         # [V] CODATA 2018, exact by SI definition

_GAN = NM.binary("GaN")
_ALN = NM.binary("AlN")
_REF_T_K = 300.0
# [A] Barrier band-offset reference temperature for the flat-barrier
# envelope. transmission() takes no T_K argument; the Varshni shift of the
# ~2.4 eV AlN-GaN offset changes by well under 0.1% over 230-300 K, which is
# negligible next to the [A] linear-interpolation and [A] partition
# uncertainties already carried by the offset itself.
_EDGES_GAN = NM.band_edges(_GAN, _REF_T_K, substrate=_GAN, strain_fraction=0.0)
_EDGES_ALN = NM.band_edges(_ALN, _REF_T_K, substrate=_ALN, strain_fraction=0.0)
_GAN_C_NM = _GAN.c_A / 10.0
_GROWTH_STEP_NM_DEFAULT = _GAN_C_NM / 2.0
# [V] GaN c-axis bilayer spacing c/2; c_A from nitride_materials' Bernardini
# PRB 56, R10024 (1997) Table II lattice constant.

_SLICE_LENGTH_M = 2.0e-10        # 0.2 nm spatial staircase step for field tilt
_MIN_SLICES_PER_SEGMENT = 8


def _default_barrier_material(al_fraction: float) -> dict:
    """Al_xGa_1-xN barrier conduction/valence offsets (relative to GaN,
    unstrained) and effective masses via linear virtual-crystal
    interpolation between the GaN and AlN endpoints of nitride_materials.
    [A] interpolation of [V] endpoints (see module docstring)."""
    if not math.isfinite(al_fraction) or not 0.0 <= al_fraction <= 1.0:
        raise ValueError("al_fraction must be finite and in [0, 1]")
    dEc = (_EDGES_ALN["Ec_eV"] - _EDGES_GAN["Ec_eV"]) * al_fraction
    dEv = (_EDGES_ALN["Ev_eV"] - _EDGES_GAN["Ev_eV"]) * al_fraction
    me = _GAN.me_z + (_ALN.me_z - _GAN.me_z) * al_fraction
    mh = _GAN.mh_z + (_ALN.mh_z - _GAN.mh_z) * al_fraction
    return dict(dEc_eV=float(dEc), dEv_eV=float(dEv), me=float(me), mh=float(mh))


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

    n_cm3: float = 3.0e18                             # [V] Deshpande et al. 2013 n-GaN doping anchor
    p_cm3: float = 5.0e17                             # [V] Deshpande et al. 2013 p-GaN doping anchor

    alignment_uncertainty_meV: float = 15.0           # [A] growth/doping-limited level-alignment uncertainty
    degeneracy: float = 2.0                           # [A] spin/valley degeneracy in the Landauer rate
    bypass_prefactor: float = 1.0                     # [A] declared prefactor on the thermionic (over-barrier) integral
    field_leverarm: float = 1.0                       # [A] fraction of the applied bias dropping across the injector

    occupancy_control_known: bool = False             # necessary, not sufficient (see module docstring)
    second_pair_control_known: bool = False

    def __post_init__(self):
        for name in ("electron_topology", "hole_topology"):
            v = getattr(self, name)
            if v not in ("single_barrier", "double_barrier"):
                raise ValueError(f"{name} must be 'single_barrier' or 'double_barrier'")
        if not math.isfinite(self.al_fraction) or not 0.0 <= self.al_fraction <= 1.0:
            raise ValueError("al_fraction must be finite and in [0, 1]")
        positive = (
            "electron_barrier_thickness_nm", "hole_barrier_thickness_nm",
            "electron_well_width_nm", "hole_well_width_nm", "growth_step_nm",
            "me_well", "mh_well", "n_cm3", "p_cm3", "alignment_uncertainty_meV",
            "degeneracy", "bypass_prefactor", "field_leverarm",
        )
        for name in positive:
            v = getattr(self, name)
            if not math.isfinite(v) or v <= 0:
                raise ValueError(f"{name} must be positive and finite")
        if not math.isfinite(self.growth_tolerance_steps) or self.growth_tolerance_steps < 0:
            raise ValueError("growth_tolerance_steps must be non-negative and finite")
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
    default = _default_barrier_material(params.al_fraction)
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
    mu = _degenerate_mu_eV(params.p_cm3 * 1e6, params.mh_well)
    return _ResolvedPath("hole", params.hole_topology,
                          params.hole_barrier_thickness_nm,
                          params.hole_well_width_nm,
                          mh_b, params.mh_well, -dEv, mu)


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
                                     m_right: float, tilt_eV_per_m: float):
    """Flux/mass-normalized transmission and reflection through a stack of
    (length_m, V_eV, m_ratio) flat-band segments, ends referenced to GaN
    reservoirs of mass m_left / m_right at V=0 (before the field tilt).  A
    linear field tilt is applied by subdividing each segment into thin
    slices (staircase approximation of a linear ramp; a standard treatment,
    e.g. Ridley, "Quantum Processes in Semiconductors").  Numerically stable
    backward coefficient recursion (BenDaniel-Duke effective-mass boundary
    condition: continuity of psi and (1/m) dpsi/dx), not a chained-matrix
    product, to avoid overflow for thick/deep barriers."""
    slices = []
    x0 = 0.0
    for (length_m, V_eV, m_ratio) in segments:
        if length_m <= 0.0:
            continue
        n = max(_MIN_SLICES_PER_SEGMENT, int(math.ceil(length_m / _SLICE_LENGTH_M)))
        dl = length_m / n
        for i in range(n):
            xc = x0 + (i + 0.5) * dl
            slices.append((dl, V_eV - tilt_eV_per_m * xc, m_ratio))
        x0 += length_m

    k_right = _segment_k(energy_eV, 0.0, m_right)
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
    return float(min(max(T, 0.0), 1.0 + 1e-9)), float(min(max(R, 0.0), 1.0 + 1e-9))


def _transmission_scalar(params, energy_eV, bias_V, field_kVcm, carrier):
    path = _resolve_path(params, carrier)
    segs = _stack_segments(path)
    total_len_m = sum(s[0] for s in segs)
    tilt = _tilt_eV_per_m(field_kVcm, bias_V, params.field_leverarm, total_len_m)
    T, R = _transmission_reflection_scalar(segs, float(energy_eV), path.m_well, path.m_well, tilt)
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
    from dataclasses import replace
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


def _find_resonance(params: NitrideNanowireInjectorParams, carrier: str,
                     bias_V: float, field_kVcm: float):
    """Locate the lowest double-barrier resonance below the barrier top.
    Seeds the search from the analytic finite-well estimate above (a real
    double barrier's true transmission resonance can be far narrower than a
    blind uniform grid can resolve for a deep/thick design), then refines
    with a multi-level local zoom on the transfer-matrix transmission
    itself, so the REPORTED energy and height are always numerically
    computed, never the analytic estimate. Falls back to a coarse full-range
    scan if the analytic seed is unavailable. Returns (E_res_eV, T_res) or
    None (single_barrier topology, or no resolvable peak -- "a single
    tunnel barrier has no claimed resonant well")."""
    path = _resolve_path(params, carrier)
    if path.topology != "double_barrier":
        return None
    e_hi = path.barrier_height_eV
    if not math.isfinite(e_hi) or e_hi <= 0:
        return None

    def T_at(E):
        return _transmission_scalar(params, E, bias_V, field_kVcm, carrier)[0]

    seed = _analytic_well_seed_eV(path.m_well, path.m_barrier, e_hi, path.well_nm * 1e-9)
    if seed is None or not (0.0 < seed < e_hi):
        grid = np.linspace(e_hi * 1e-4, e_hi * 0.999, 600)
        Ts = np.array([T_at(E) for E in grid])
        i = int(np.argmax(Ts))
        if Ts[i] < 1e-12:
            return None
        lo, hi = grid[max(i - 2, 0)], grid[min(i + 2, len(grid) - 1)]
        best_E, best_T = float(grid[i]), float(Ts[i])
    else:
        lo, hi = seed * 0.9, min(seed * 1.1, e_hi * 0.999999)
        best_E, best_T = seed, T_at(seed)
        for _ in range(4):
            local = np.linspace(lo, hi, 60)
            Ts_local = np.array([T_at(E) for E in local])
            j = int(np.argmax(Ts_local))
            if Ts_local[j] > best_T:
                best_E, best_T = float(local[j]), float(Ts_local[j])
            span = (hi - lo) / 30.0
            lo = max(local[j] - span, 1e-12)
            hi = min(local[j] + span, e_hi * 0.999999)
            if hi <= lo:
                break

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
                      combined_width_eV: float, kT_eV: float):
    """Single-channel Landauer capture rate into an EMPTY target state,
    Gamma = (g_s/h) INTEGRAL T(E) f_reservoir(E) dE (see module docstring).
    Also returns the same integral restricted to energies above the
    flat-band barrier top (the thermionic/over-barrier share, used for
    bypass pricing without double-counting: it is the SAME T(E) integral,
    partitioned by energy domain, not a second exp(-barrier/kT) estimate)."""
    path = _resolve_path(params, carrier)
    mu = path.mu_eV
    width = combined_width_eV if math.isfinite(combined_width_eV) else 0.0
    e_lo = max(1e-6, min(E_center, mu) - 15.0 * kT_eV - 5.0 * width)
    e_hi = max(path.barrier_height_eV, E_center, mu) + 15.0 * kT_eV

    def integrand(E):
        T = _transmission_scalar(params, E, bias_V, field_kVcm, carrier)[0]
        return T * _fermi_dirac(E, mu, kT_eV)

    hint_points = sorted({p for p in (E_center, path.barrier_height_eV) if e_lo < p < e_hi})
    total, _ = quad(integrand, e_lo, e_hi, points=hint_points or None,
                     limit=200, epsabs=1e-18, epsrel=1e-7)
    if path.barrier_height_eV < e_hi:
        above, _ = quad(integrand, path.barrier_height_eV, e_hi,
                         limit=200, epsabs=1e-18, epsrel=1e-7)
    else:
        above = 0.0
    rate_total = params.degeneracy * total / H_EVS
    rate_above = params.degeneracy * above / H_EVS
    return float(rate_total), float(rate_above)


def _carrier_screen(params: NitrideNanowireInjectorParams, carrier: str, *,
                     kT_eV: float, level_eV: float, spacing_meV: float,
                     second_pair_addition_meV: float, available_pair_rate_Hz: float,
                     loading_window_ns: float, rep_rate_hz: float, field_kVcm: float):
    """All diagnostics for one carrier path at params' current (possibly
    growth-perturbed) thickness. Held at zero series bias: the interface
    constraint exposes only field_kVcm to injector_feasibility, and
    params.field_leverarm lets a card fold an equivalent series-bias
    contribution into field_kVcm [A, documented simplification -- bias_V
    itself remains an independent transmission()/reflection() argument for
    direct callers]."""
    bias_V = 0.0
    path = _resolve_path(params, carrier)

    res = _find_resonance(params, carrier, bias_V, field_kVcm)
    if path.topology == "double_barrier" and res is not None:
        E_res, T_res = res
        linewidth_eV = _resonance_width_eV(params, carrier, bias_V, field_kVcm,
                                            E_res, T_res, path.barrier_height_eV)
        alignment_error_meV = abs(E_res - level_eV) * 1000.0
        E_center = E_res
    else:
        # Single tunnel barrier: no claimed resonant well, so no resonance
        # to misalign against [A]; the only broadening is the declared
        # alignment uncertainty itself.
        linewidth_eV = params.alignment_uncertainty_meV / 1000.0
        alignment_error_meV = 0.0
        E_center = level_eV

    linewidth_meV = linewidth_eV * 1000.0 if math.isfinite(linewidth_eV) else float("nan")
    combined_width_eV = math.sqrt(
        (linewidth_eV if math.isfinite(linewidth_eV) else 0.0) ** 2
        + (params.alignment_uncertainty_meV / 1000.0) ** 2
    )

    rate_hz, rate_above_hz = _forward_rate_hz(params, carrier, bias_V, field_kVcm,
                                               E_center, combined_width_eV, kT_eV)
    rate_below_hz = max(rate_hz - rate_above_hz, 0.0)
    bypassed = params.bypass_prefactor * rate_above_hz
    denom = rate_below_hz + bypassed
    bypass_fraction = float(bypassed / denom) if denom > 0 else 0.0

    orbital_candidates = [c for c in (
        spacing_meV / 1000.0 if math.isfinite(spacing_meV) else float("nan"),
        path.barrier_height_eV - E_center,
    ) if math.isfinite(c) and c > -1e9]
    orbital_sep_eV = min(orbital_candidates) if orbital_candidates else float("nan")
    full_candidates = list(orbital_candidates)
    if math.isfinite(second_pair_addition_meV):
        full_candidates.append(second_pair_addition_meV / 1000.0)
    full_sep_eV = min(full_candidates) if full_candidates else float("nan")

    def _margin_kT(sep_eV):
        if not math.isfinite(sep_eV):
            return float("nan")
        return max(sep_eV - combined_width_eV, 0.0) / kT_eV

    margin_kT = _margin_kT(full_sep_eV)
    orbital_margin_kT = _margin_kT(orbital_sep_eV)

    t_load_s = loading_window_ns * 1e-9
    missed_load_p = math.exp(-rate_hz * t_load_s)

    period_s = 1.0 / rep_rate_hz
    t_remaining_s = max(period_s - t_load_s, 0.0)
    if math.isfinite(second_pair_addition_meV) and math.isfinite(available_pair_rate_Hz):
        suppression = math.exp(-second_pair_addition_meV / 1000.0 / kT_eV)
        rate_second_hz = min(rate_hz, available_pair_rate_Hz) * suppression
        second_pair_p = 1.0 - math.exp(-rate_second_hz * t_remaining_s)
    else:
        rate_second_hz = float("nan")
        second_pair_p = float("nan")

    return dict(
        rate_hz=rate_hz, alignment_error_meV=alignment_error_meV,
        linewidth_meV=linewidth_meV, bypass_fraction=bypass_fraction,
        margin_kT=margin_kT, orbital_margin_kT=orbital_margin_kT,
        missed_load_p=missed_load_p, second_pair_p=second_pair_p,
        rate_second_hz=rate_second_hz, barrier_height_eV=path.barrier_height_eV,
        thickness_nm=path.thickness_nm,
    )


def _transport_ok(screen: dict, margin_floor_kT: float, bypass_ceiling: float,
                   missed_ceiling: float) -> bool:
    return (math.isfinite(screen["margin_kT"]) and screen["margin_kT"] >= margin_floor_kT
            and math.isfinite(screen["bypass_fraction"]) and screen["bypass_fraction"] <= bypass_ceiling
            and math.isfinite(screen["missed_load_p"]) and screen["missed_load_p"] <= missed_ceiling)


# --------------------------------------------------------- engineering thresholds
# [A] Proposed conditional engineering thresholds (design-brief round 3).

RTI_MARGIN_FLOOR_KT = 10.0
RTI_BYPASS_CEILING = 0.01
RTI_MISSED_LOAD_CEILING = 0.01
RTI_SECOND_PAIR_CEILING = 0.01

_PROVENANCE = (
    "[V] material endpoints (GaN/AlN masses, offsets, lattice constant) from "
    "fsim_core.nitride_materials (Bernardini PRB 1997; Bernardini & "
    "Fiorentini pss(b) 1999; Rinke PRB 2008; Tsai & Bayram ACS Omega 2020); "
    "[A] linear Al_xGa_1-xN interpolation of those endpoints, unstrained "
    "barrier; [V] n/p-GaN doping anchors Deshpande et al., Nat. Commun. 4, "
    "1675 (2013); [DR] degenerate free-electron-gas reservoir "
    "electrochemical energy (Ashcroft & Mermin, 'Solid State Physics', "
    "1976, Ch. 2 Eq. 2.33); [V] single-channel Landauer/Buttiker rate "
    "normalization 1/h (e.g. Datta, 'Electronic Transport in Mesoscopic "
    "Systems', 1995, Ch. 2); [A] spin/valley degeneracy, bypass prefactor, "
    "field lever-arm, alignment uncertainty, engineering thresholds "
    "(margin 10 kT, bypass/missed-load/second-pair 0.01). Conditional "
    "engineering screening only; not a demonstrated hardware claim."
)


def _unknown_result(failed_checks, status: str) -> dict:
    nan = float("nan")
    return dict(
        rti_feasible=False, rti_status=status, rti_transport_feasible=False,
        rti_level_margin_kT=nan, rti_orbital_margin_kT=nan,
        rti_alignment_error_meV=nan, rti_linewidth_meV=nan,
        rti_rate_Hz=nan, rti_e_rate_Hz=nan, rti_h_rate_Hz=nan,
        rti_bypass_fraction=nan, rti_missed_load_probability=nan,
        rti_second_pair_probability=nan, rti_growth_feasible=False,
        rti_failed_checks=list(failed_checks),
        rti_evidence_status="conditional_engineering_screen_no_measured_device",
        valid=False, provenance=_PROVENANCE,
    )


def injector_feasibility(params: NitrideNanowireInjectorParams, *, T_K: float,
                          rep_rate_hz: float, loading_window_ns: float,
                          electron_level_eV: float, hole_level_eV: float,
                          electron_spacing_meV: float, hole_spacing_meV: float,
                          second_pair_addition_meV: float,
                          available_pair_rate_Hz: float,
                          field_kVcm: float = 0.0) -> dict:
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
    computation behind it."""
    raw = dict(T_K=T_K, rep_rate_hz=rep_rate_hz, loading_window_ns=loading_window_ns,
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
    if not _nonfinite(available_pair_rate_Hz) and available_pair_rate_Hz < 0:
        range_bad.append("out_of_range:available_pair_rate_Hz")
    if ("loading_window_ns" not in missing and "rep_rate_hz" not in missing
            and loading_window_ns * 1e-9 > 1.0 / rep_rate_hz):
        range_bad.append("out_of_range:loading_window_exceeds_period")

    failed = [f"missing:{k}" for k in missing] + range_bad
    if failed:
        return _unknown_result(failed, "missing_or_invalid_input")

    kT_eV = KB_EV * T_K

    def _screen_at(p):
        e = _carrier_screen(p, "electron", kT_eV=kT_eV, level_eV=electron_level_eV,
                             spacing_meV=electron_spacing_meV,
                             second_pair_addition_meV=second_pair_addition_meV,
                             available_pair_rate_Hz=available_pair_rate_Hz,
                             loading_window_ns=loading_window_ns, rep_rate_hz=rep_rate_hz,
                             field_kVcm=field_kVcm)
        h = _carrier_screen(p, "hole", kT_eV=kT_eV, level_eV=hole_level_eV,
                             spacing_meV=hole_spacing_meV,
                             second_pair_addition_meV=second_pair_addition_meV,
                             available_pair_rate_Hz=available_pair_rate_Hz,
                             loading_window_ns=loading_window_ns, rep_rate_hz=rep_rate_hz,
                             field_kVcm=field_kVcm)
        return e, h

    e, h = _screen_at(params)

    failed_checks = []
    e_ok = _transport_ok(e, RTI_MARGIN_FLOOR_KT, RTI_BYPASS_CEILING, RTI_MISSED_LOAD_CEILING)
    h_ok = _transport_ok(h, RTI_MARGIN_FLOOR_KT, RTI_BYPASS_CEILING, RTI_MISSED_LOAD_CEILING)
    if not e_ok:
        failed_checks.append("electron_transport")
    if not h_ok:
        failed_checks.append("hole_transport")

    missed_load_probability = 1.0 - (1.0 - e["missed_load_p"]) * (1.0 - h["missed_load_p"])
    if missed_load_probability > RTI_MISSED_LOAD_CEILING:
        failed_checks.append("pair_missed_load_window")

    rti_transport_feasible = bool(
        e_ok and h_ok and missed_load_probability <= RTI_MISSED_LOAD_CEILING
    )

    # Growth realism: perturb both barriers together by +/- N growth steps
    # and require the necessary transport conditions to survive at every
    # perturbation; a design that only clears at the exact nominal
    # thickness is fragile, never a manufacturable claim.
    growth_ok = rti_transport_feasible
    if params.growth_tolerance_steps > 0:
        for sign in (-1.0, 1.0):
            delta = sign * params.growth_tolerance_steps * params.growth_step_nm
            p_pert = _with_thickness(params, delta)
            e_p, h_p = _screen_at(p_pert)
            e_p_ok = _transport_ok(e_p, RTI_MARGIN_FLOOR_KT, RTI_BYPASS_CEILING, RTI_MISSED_LOAD_CEILING)
            h_p_ok = _transport_ok(h_p, RTI_MARGIN_FLOOR_KT, RTI_BYPASS_CEILING, RTI_MISSED_LOAD_CEILING)
            missed_pert = 1.0 - (1.0 - e_p["missed_load_p"]) * (1.0 - h_p["missed_load_p"])
            growth_ok = growth_ok and e_p_ok and h_p_ok and missed_pert <= RTI_MISSED_LOAD_CEILING
    nominal_steps_e = params.electron_barrier_thickness_nm / params.growth_step_nm
    nominal_steps_h = params.hole_barrier_thickness_nm / params.growth_step_nm
    submonolayer = (abs(nominal_steps_e - round(nominal_steps_e)) > 0.05
                    or abs(nominal_steps_h - round(nominal_steps_h)) > 0.05)
    rti_growth_feasible = bool(growth_ok and not submonolayer)
    if not rti_growth_feasible:
        failed_checks.append("growth_tolerance")

    second_pair_known = bool(math.isfinite(e["second_pair_p"]) and math.isfinite(h["second_pair_p"]))
    if not second_pair_known:
        failed_checks.append("second_pair_exclusion_unknown")
        second_pair_probability = float("nan")
    else:
        second_pair_probability = 1.0 - (1.0 - e["second_pair_p"]) * (1.0 - h["second_pair_p"])
        if second_pair_probability > RTI_SECOND_PAIR_CEILING:
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

    rti_level_margin_kT = min(e["margin_kT"], h["margin_kT"])
    rti_orbital_margin_kT = min(e["orbital_margin_kT"], h["orbital_margin_kT"])
    rti_alignment_error_meV = max(e["alignment_error_meV"], h["alignment_error_meV"])
    rti_linewidth_meV = (max(e["linewidth_meV"], h["linewidth_meV"])
                          if math.isfinite(e["linewidth_meV"]) and math.isfinite(h["linewidth_meV"])
                          else float("nan"))
    rti_bypass_fraction = max(e["bypass_fraction"], h["bypass_fraction"])
    rti_rate_Hz = min(e["rate_hz"], h["rate_hz"])   # slower completed-pair rate

    critical = [rti_level_margin_kT, rti_alignment_error_meV, e["rate_hz"], h["rate_hz"],
                rti_bypass_fraction, missed_load_probability]
    valid = bool(all(math.isfinite(v) for v in critical))
    if not valid:
        rti_feasible = False
        if rti_status == "feasible":
            rti_status = "infeasible"

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
    )
