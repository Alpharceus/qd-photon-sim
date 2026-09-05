"""Module T -- Transport: p-i-n diode model of the actual layer stack
(RT edge-emitter tier, 2026-09-02).

WHY THIS MODULE EXISTS. Before this tier the drive block carried three free
[A] inputs -- mu (per-pulse cap-2 loading), eta_capture (F8b thinning) and
b_e (injection background amplitude, plus its exponents b_e_m / b_e_Eact) --
and the thermal model took the diode voltage V as another free input. The
code audit found b_e used with three inconsistent meanings across the chain.
This module computes all of them from a p-i-n diode built from
`materials.Material` layers, and it defines b_e ONCE:

    b_e(I, T)  =  [matrix/WL background photons per second that fall inside
                   the X spectral window]  /  [X photons per second of the
                   selected dot]           at the operating point (I, T).

    Dimensionless.  Same normalization as loading.b_injection (background per
    signal photon, signal at full transmission).  Neighbour dots inside the
    aperture are NOT included here: their contribution is the F5 aperture
    lemma (loading.aperture_g2 / loading.n_window_competitors).  Adding them
    here too would double count.  The window transmission of the X line
    itself (t_X) is also not applied here: b_e is "per X photon emitted into
    the line", the filter acts on both downstream.

PROVENANCE TAGS.  [DR] derived from the materials database (Vurgaftman 2001
compilation) or textbook closed forms; [E] estimates -- every [E] default is
named in the `Diode` field comments and in ASSUMPTIONS_E below.

PHYSICS (all closed forms; one implementation each):

  1. Carrier statistics.  3-D effective densities of states
         N_c = 2 (2 pi m_e k T / h^2)^{3/2},   N_v likewise with the DOS hole
         mass  m_dos,h = (m_hh^{3/2} + m_lh^{3/2})^{2/3},
     with m_hh, m_lh the spherically averaged Luttinger masses
     1/(g1 -+ 2 gbar), gbar = (2 g2 + 3 g3)/5 (Chuang ch. 4; GaAs -> 0.55 /
     0.082, the standard DOS values).  Indirect materials (AlInP, AlGaInP
     x > 0.53): the conduction DOS is the valley sum
         N_c,eff = sum_v N_c,v exp(-(E_c,v - E_c,min)/kT)
     with [E] X/L DOS masses (XL_DOS_MASS) -- the materials module carries
     no X/L masses.  n_i = sqrt(N_c,eff N_v) exp(-E_g,min/2kT), evaluated in
     log space so that 20 K does not underflow.

  2. Built-in voltage of the heterojunction p-i-n:
         V_bi = (E_Fn,n-side - E_Fp,p-side)/q,
         E_Fn = E_c,n,min - kT ln(N_c,eff,n/N_D),  E_Fp = E_v,p + kT ln(N_v,p/N_A)
     on the absolute (V01 VBO) scale from materials.band_edges.  Boltzmann
     statistics; a doping above the effective DOS triggers a warning
     (degenerate: the Boltzmann Fermi level is then off by ~kT).  For a
     homojunction this collapses exactly to (kT/q) ln(N_A N_D / n_i^2).

  3. Electrostatics.  Abrupt p-i-n, i-region fully depleted, doped-side
     depletion widths from the standard quadratic
         V_bi - V_j = (q N_A x_p / eps) [ d_i + x_p (1 + N_A/N_D)/2 ],
         x_n = x_p N_A/N_D,   F_i = q N_A x_p / eps  (uniform in the i-region)
     If V_j >= V_bi: flat band, F = 0.  Field returned in kV/cm, plus the
     potential drop across the active layer, total depletion width and the
     depletion capacitance (feeds the F9 granularity block).  The QCSE
     placeholder qcse_shift_meV = -p F^2 / 2 has p = 0 by default: the dot
     polarizability is NOT invented here.

  4. I-V.  J(V_j) = sum over three channels, V_applied = V_j + I R_s:
       * SRH generation-recombination in the depleted i-region, ideality
         n_ideality (default 2):  J0,SRH = q sum_layers n_i,layer d_layer /
         (2 tau_SRH).  The spec's one-liner "q n_i,active d_i / (2 tau)" is
         the special case barrier == active; resolving the layers matters
         because n_i,barrier/n_i,active ~ exp(-dE_g/2kT) ~ 1e-4.
       * Radiative recombination in the active layer, n = 1:
         J0,rad = q B d_active n_i,active^2  (B = B_rad_cm3s [E]).  This is
         the physical n = 1 channel of an LED.
       * Shockley diffusion of minority carriers into the neutral claddings,
         n = 1:  J0,diff = q [ n_i,p^2 D_n/(L_n N_A) + n_i,n^2 D_p/(L_p N_D) ]
         with each cladding's OWN n_i (the spec text wrote n_i,active, which
         is the homojunction special case; for an Al-rich cladding this
         channel is negligible, exactly as the spec remarks, because
         n_i,clad^2/n_i,active^2 ~ exp(-2 dE_g/kT)).  D = mu kT/q,
         L = sqrt(D tau), mu_n = 2000, mu_p = 100 cm^2/Vs, tau = 1 ns [E].
     V_j(I) is solved by brentq on ln J.  V_on = V_j at J = 10 A/cm^2.

  5. Injection efficiency (electron leakage over the p-cladding, the classic
     AlGaInP problem; Bour et al., IEEE JQE 30, 593 (1994); Coldren, Corzine
     & Mashanovitch, "Diode Lasers and Photonic Integrated Circuits", Sec.
     4.4-4.5 "carrier leakage"):
         eta_inj = 1 / (1 + r_e + r_h)
         r_e = sum_v (N_c,v,pclad / N_c,active) (D_n tau_active / (L_n d_active))
               exp(-(E_c,v,pclad - E_c,active - E_Fn,offset) / kT)
     where the valley sum v runs over the p-cladding valleys selected by
     `leak_valleys`.  Default "G": the Gamma-only expression with dE_c from
     materials.offsets ([DR] strained model-solid edges) -- a LOWER BOUND on
     the leakage.  "GX" adds the cladding X valley with an [E] DOS mass; for
     AlInP the X edge lies ~0.33 eV below its Gamma edge and only 0.08-0.12
     eV above the active CB in the V01 model-solid picture, and the
     Gamma-X leakage is the known AlGaInP limitation (Bour 1994).  "GX" is
     not the default because the result moves 50x per 0.1 eV of X-edge
     uncertainty (literature GaInP/AlInP Gamma-X offsets span 0.13-0.28 eV)
     and the X-electron mobility/DOS are [E]; verify_transport prints both
     modes so the choice is visible, never silent.
     E_Fn,offset = qfl_e_meV = 0 [E] (electrons at the active band edge).
     L_n = sqrt(D_n tau_n), capped at the p-cladding thickness d_pclad_nm if
     given (the GaInP/GaAs contact layers are a perfect sink).  Hole leakage
     r_h symmetric with dE_v (top valence band, strain-resolved).
     The X/L edges of the cladding are UNSTRAINED (no X/L deformation
     potentials in materials); the claddings are lattice matched so this is
     a < 10 meV effect.

  6. Carrier budget at the dot layer.  Carriers/s reaching the active layer
     = eta_inj I/q.  Fraction captured by dots vs recombining in the matrix
     / wetting layer:  f_QD = 1/(1 + tau_cap,QD/tau_matrix),
     tau_cap,QD = tau_cap0 (1e10 cm^-2 / n_dot) [E: capture rate proportional
     to dot density, tau_cap0 = 5 ps at 1e10], tau_matrix = 1 ns [E].
     N_dots = n_dot A_ap;  r_dot = f_QD eta_inj (I/q)/N_dots;
     pulsed: mu = r_dot tau_pulse (cap-2 truncation happens downstream in
     loading.f1b_g2 / f8_g2);  CW: occupancy estimate r_dot tau_rad.
     THE EYE-OPENER: 1 uA, 1 ns pulses, one dot -> mu ~ 6000.  Single-dot
     cap-2 operation (mu ~ 0.5) needs ~80 pA into ONE dot, or ~0.6 nA for the
     7 dots under a 1 um aperture at 7e8 cm^-2.  Every published electrically
     driven single-dot result at mA currents therefore runs either with
     millions of dots sharing the current (Reischle 2010: 20 um aperture,
     1e10 cm^-2 -> 3e6 dots) or on a filamentary / aperture-limited current
     path.  A designed device must make the current path that small on
     purpose (oxide aperture, nanowire p-i-n) -- it is a hard number, not a
     tuning knob.

  7. Injection background in the X window (definition at the top).
     numerator = (1 - f_QD) eta_inj (I/q) eta_rad,matrix xi(w, dE_WL, T)
     with xi the fraction of the matrix/WL luminescence inside a window of
     full width w centred dE_WL below the WL/matrix peak.  Low-energy tail
     modelled as a normalized one-sided Urbach/thermal exponential
         f(x) = (1/E_U) exp(-x/E_U),  x = E_WL - E >= 0,  int_0^inf f = 1,
         E_U = max(kT, E_urbach)   [E]
     so that (window = [dE_WL - w/2, dE_WL + w/2] in x, clipped at x = 0)
         xi = exp(-max(0, dE_WL - w/2)/E_U) - exp(-(dE_WL + w/2)/E_U)
            = exp(-dE_WL/E_U) * 2 sinh(w/(2 E_U))     when dE_WL >= w/2,
     which is the Arrhenius-in-T form loading.BackgroundChannel assumed
     (E_act ~ dE_WL, with a slowly varying prefactor 2 sinh(w/2kT) ~ w/kT).
     The clip guarantees xi <= 1 (w -> inf gives xi -> 1); for dE_WL = 0 and
     w << E_U the clipped form gives w/(2 E_U) -- half of the unclipped
     w/E_U, because half the window lies above the peak where the one-sided
     tail carries no weight.
     denominator = min(r_dot, 1/tau_rad) S_dot(T)  (X photons/s of the dot;
     the X-line share of the dot's photons is taken as 1).
     Current dependence: below dot saturation (r_dot < 1/tau_rad) b_e is
     independent of I (both scale with I); above it b_e grows linearly
     (m = 1).  `to_background_channel` least-squares fits ln b = ln A +
     m ln(I/I_ref) - E_act/kT over a grid and returns the legacy
     BackgroundChannel plus the fit error, so the fit quality is visible.

  8. Heating.  P_junction = I V_j + f_Rs,local I^2 R_s - (I/q) eta_total h nu.
     The thermal module (thermal.t_junction) takes P = duty * P_junction.

References: S. M. Sze & K. K. Ng, "Physics of Semiconductor Devices" 3rd ed.
(n_i, V_bi, ideal diode, SRH); Bour et al. 1994 (AlGaInP leakage); Coldren-
Corzine-Mashanovitch Sec. 4.4 (leakage current expression); Vurgaftman et
al. JAP 89, 5815 (2001) via materials.py.
"""
from __future__ import annotations

import warnings
from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import brentq
from scipy.special import logsumexp

from .loading import BackgroundChannel
from .materials import (
    KB_EV,
    AlGaAs,
    AlGaInP,
    AlInP,
    GaAsP,
    GaInP,
    Material,
    band_edges,
    bandgap,
    binary,
    offsets,
    strain_shifts,
)

# ------------------------------------------------------------------ SI constants
Q_SI = 1.602176634e-19      # C
KB_SI = 1.380649e-23        # J/K
H_SI = 6.62607015e-34       # J s
M0_SI = 9.1093837015e-31    # kg
EPS0_SI = 8.8541878128e-12  # F/m

# 2 (2 pi m0 k_B 300 K / h^2)^{3/2} in cm^-3  (= 2.509e19)
NC300_CM3 = 2.0 * (2.0 * np.pi * M0_SI * KB_SI * 300.0 / H_SI**2) ** 1.5 * 1e-6

# [E] density-of-states masses (in m0) of the satellite valleys, INCLUDING the
# valley multiplicity: X ~ 0.85 (GaP-class 0.79, AlAs 0.71-0.85; Ioffe NSM),
# L ~ 0.56 (GaAs L, Ioffe).  materials.py carries only the Gamma mass.
XL_DOS_MASS = {"X": 0.85, "L": 0.56}

ASSUMPTIONS_E = {
    "n_ideality": "2.0 -- SRH-dominated i-region recombination at low bias",
    "tau_SRH_ns": "1.0 -- SRH lifetime in the i-region (quality dependent, 1-100 ns)",
    "R_s_ohm": "50 -- series resistance of a small mesa incl. contacts",
    "f_Rs_local": "0.3 -- fraction of the I^2 R_s heat dissipated at the mesa",
    "mu_n_cm2 / mu_p_cm2": "2000 / 100 cm^2/Vs -- cladding minority mobilities (an "
                           "Al-rich, X-valley cladding is likely 5-10x lower for "
                           "electrons: r_leak scales as sqrt(mu_n))",
    "tau_n_ns / tau_p_ns": "1.0 -- cladding minority lifetimes (L = sqrt(D tau))",
    "tau_active_ns": "1.0 -- active-layer recombination time in the leakage ratio",
    "qfl_e_meV": "0 -- electron quasi-Fermi level at the active band edge",
    "B_rad_cm3s": "1e-10 -- bimolecular coefficient of the active (GaInP/GaAs class)",
    "eps_r": "12.0 -- static dielectric constant of the i-region (GaInP 11.8, AlInP ~11.3)",
    "XL_DOS_MASS": "X: 0.85, L: 0.56 -- satellite-valley DOS masses incl. multiplicity",
    "X/L edges unstrained": "no X/L deformation potentials in materials.py",
    "tau_cap0_ps": "5 ps -- dot capture time at 1e10 cm^-2, rate proportional to density",
    "tau_matrix_ns": "1.0 -- matrix / wetting-layer recombination lifetime",
    "eta_rad_matrix": "0.1 -- radiative efficiency of matrix/WL recombination",
    "E_urbach": "E_U = max(kT, E_urbach) -- one-sided exponential tail of the WL/matrix line",
    "X-line share": "1 -- all dot photons counted in the X line for the b_e denominator",
    "eta_total": "0.01 -- wall-plug photon efficiency in the heat balance (negligible)",
    "3-D DOS for the active": "the 8 nm well/dot layer is treated with bulk N_c, N_v",
}


# ========================================================== carrier statistics

def hole_dos_mass(m: Material) -> float:
    """m_dos,h = (m_hh^{3/2} + m_lh^{3/2})^{2/3} with spherically averaged
    Luttinger masses m_hh = 1/(g1 - 2 gbar), m_lh = 1/(g1 + 2 gbar),
    gbar = (2 g2 + 3 g3)/5.  GaAs -> 0.55 / 0.082 -> m_dos,h = 0.60 [DR]."""
    g1, g2, g3 = m.p["g1"], m.p["g2"], m.p["g3"]
    # Spherical Luttinger average used for the standard DOS masses [DR].
    # The former 2:3 weighting was inconsistent with the GaAs reference
    # masses (m_hh ~= 0.50, m_lh ~= 0.084) and inflated N_v and n_i.
    gbar = (g2 + g3) / 2.0
    m_hh = 1.0 / (g1 - 2.0 * gbar)
    m_lh = 1.0 / (g1 + 2.0 * gbar)
    return float((m_hh**1.5 + m_lh**1.5) ** (2.0 / 3.0))


def effective_dos(mass: float, T: float) -> float:
    """3-D effective density of states [cm^-3] for a DOS mass (units of m0)."""
    return NC300_CM3 * mass**1.5 * (T / 300.0) ** 1.5


def conduction_valleys(m: Material, T: float = 300.0, dE_c_gamma: float = 0.0) -> list:
    """[(valley, E_c,v absolute [eV], N_c,v [cm^-3])] for G, X, L.  The Gamma
    edge may carry a strain shift dE_c_gamma (from materials.strain_shifts);
    X and L are unstrained (see module docstring)."""
    E_v0, E_c0 = band_edges(m, T)
    out = [("G", E_c0 + dE_c_gamma, effective_dos(m.p["m_e"], T))]
    for v in ("X", "L"):
        out.append((v, E_v0 + bandgap(m, T, v), effective_dos(XL_DOS_MASS[v], T)))
    return out


def n_c_eff(m: Material, T: float = 300.0, dE_c_gamma: float = 0.0) -> tuple[float, float]:
    """(N_c,eff, E_c,min): valley-summed conduction DOS referred to the lowest
    conduction edge, N_c,eff = sum_v N_c,v exp(-(E_c,v - E_c,min)/kT)."""
    vals = conduction_valleys(m, T, dE_c_gamma)
    E_min = min(E for _, E, _ in vals)
    kT = KB_EV * T
    return float(sum(N * np.exp(-(E - E_min) / kT) for _, E, N in vals)), float(E_min)


def n_v(m: Material, T: float = 300.0) -> float:
    return effective_dos(hole_dos_mass(m), T)


def ln_n_i(m: Material, T: float = 300.0, Eg: float | None = None) -> float:
    """ln n_i [ln cm^-3].  Eg (eV) overrides the material's minimum gap (used
    for a strained active layer); the conduction DOS is the valley sum."""
    Nc, E_cmin = n_c_eff(m, T)
    if Eg is None:
        Eg = E_cmin - band_edges(m, T)[0]
    return float(0.5 * (np.log(Nc) + np.log(n_v(m, T))) - Eg / (2.0 * KB_EV * T))


def n_i(m: Material, T: float = 300.0, Eg: float | None = None) -> float:
    """Intrinsic carrier density [cm^-3] (exp of ln_n_i; may underflow to 0
    below ~30 K for wide gaps -- use ln_n_i in that regime)."""
    return float(np.exp(ln_n_i(m, T, Eg)))


def fermi_levels(n_side: Material, N_D: float, p_side: Material, N_A: float,
                 T: float = 300.0) -> tuple[float, float]:
    """(E_Fn on the n side, E_Fp on the p side) [eV, absolute V01 scale],
    Boltzmann statistics; warns when the doping exceeds the effective DOS."""
    kT = KB_EV * T
    Nc, E_cmin = n_c_eff(n_side, T)
    Nv = n_v(p_side, T)
    if N_D > Nc:
        warnings.warn(f"n-side {n_side.label}: N_D = {N_D:.2e} > N_c,eff = {Nc:.2e} "
                      f"cm^-3 at {T:.0f} K -- degenerate, Boltzmann E_Fn off by ~kT",
                      stacklevel=2)
    if N_A > Nv:
        warnings.warn(f"p-side {p_side.label}: N_A = {N_A:.2e} > N_v = {Nv:.2e} "
                      f"cm^-3 at {T:.0f} K -- degenerate, Boltzmann E_Fp off by ~kT",
                      stacklevel=2)
    E_Fn = E_cmin - kT * np.log(Nc / N_D)
    E_Fp = band_edges(p_side, T)[0] + kT * np.log(Nv / N_A)
    return float(E_Fn), float(E_Fp)


def homojunction_vbi(m: Material, N_A: float, N_D: float, T: float = 300.0) -> float:
    """(kT/q) ln(N_A N_D / n_i^2) [V] -- the textbook homojunction form, used
    as the verify cross-check of the Fermi-level route."""
    return float(KB_EV * T * (np.log(N_A) + np.log(N_D) - 2.0 * ln_n_i(m, T)))


def qcse_shift_meV(F_kVcm, height_nm: float = 0.0, polarizability: float = 0.0):
    """[E] placeholder second-order Stark shift  dE = -p F^2 / 2  [meV] with
    p the exciton polarizability in meV/(kV/cm)^2, a USER parameter (default
    0: no shift).  No value is invented here: for a self-assembled dot p
    scales roughly as height^4 and must be measured or computed (8-band k.p);
    `height_nm` is accepted for that future scaling and is otherwise unused."""
    return -0.5 * polarizability * np.asarray(F_kVcm, dtype=float) ** 2


def xi_window(w_meV, dE_WL_meV, E_U_meV):
    """Fraction of a normalized one-sided exponential tail (decay energy E_U)
    inside a window of full width w centred dE_WL below the peak, clipped at
    the peak:  xi = exp(-max(0, dE - w/2)/E_U) - exp(-(dE + w/2)/E_U) <= 1.
    Reduces to exp(-dE/E_U) 2 sinh(w/2E_U) when dE >= w/2 (module docstring)."""
    w = np.asarray(w_meV, dtype=float)
    dE = np.asarray(dE_WL_meV, dtype=float)
    lo = np.maximum(0.0, dE - 0.5 * w)
    return np.exp(-lo / E_U_meV) - np.exp(-(dE + 0.5 * w) / E_U_meV)


def xi_window_unclipped(w_meV, dE_WL_meV, E_U_meV):
    """The spec's closed form exp(-dE/E_U) 2 sinh(w/2E_U) without the clip
    (exceeds 1 for a window wider than the tail; kept for the derivation
    check only -- xi_window is what b_e uses)."""
    return np.exp(-np.asarray(dE_WL_meV, float) / E_U_meV) * 2.0 * np.sinh(
        0.5 * np.asarray(w_meV, float) / E_U_meV)


# ================================================================ result types

@dataclass
class Depletion:
    V_j: float
    V_bi: float
    F_kVcm: float          # field in the i-region
    dV_active_mV: float    # potential drop across d_active
    x_p_nm: float
    x_n_nm: float
    W_nm: float            # x_p + d_i + x_n
    C_dep_pF: float        # eps A / W
    flat_band: bool


@dataclass
class Leakage:
    T: float
    dE_c_eff_eV: float     # lowest electron barrier used (active -> p-cladding)
    dE_v_eff_eV: float     # hole barrier (active -> n-cladding)
    r_leak_e: float
    r_leak_h: float
    eta_e: float           # 1/(1 + r_e)
    eta_h: float
    eta_inj: float         # 1/(1 + r_e + r_h)
    valleys: str
    per_valley: dict       # valley -> (barrier_eV, r_leak contribution)


@dataclass
class DotLoading:
    I_uA: float
    T: float
    eta_inj: float
    f_QD: float
    N_dots: float
    r_dot: float           # captured carriers per dot per second
    mu: float | None       # per-pulse loading (pulsed) -- cap-2 applied downstream
    n_occ_cw: float        # r_dot * tau_rad (CW mean-occupancy estimate, uncapped)
    saturated: bool        # r_dot > 1/tau_rad
    eta_capture_dot: float  # per-dot thinning probability of the injected stream


@dataclass
class Background:
    b_e: float
    T: float
    I_uA: float
    xi: float
    E_U_meV: float
    eta_inj: float
    f_QD: float
    rate_bg_window: float   # background photons/s inside the window
    rate_x: float           # dot X photons/s
    r_dot: float
    S_dot: float
    saturated: bool
    note: str = ("neighbour dots enter only via loading.aperture_g2 / "
                 "n_window_competitors -- not counted here")


@dataclass
class InjectionResult:
    I_uA: float
    T: float
    V_j: float
    V_applied: float
    V_bi: float
    V_on: float
    J_Acm2: float
    depletion: Depletion
    leakage: Leakage
    loading: DotLoading
    background: Background
    P_junction_W: float
    mu: float | None
    eta_capture: float
    b_e: float
    tags: dict = field(default_factory=lambda: {
        "V_j": "DR", "V_bi": "DR", "F": "DR(eps_r E)", "eta_inj": "E",
        "mu": "E", "b_e": "E", "P_junction": "E"})


# ===================================================================== Diode

@dataclass
class Diode:
    """A p-i-n diode built from Material layers.  Field comments carry the
    provenance tag; every [E] default is also listed in ASSUMPTIONS_E."""
    p_cladding: Material
    n_cladding: Material
    active: Material            # lowest-gap i-layer (GaAsP well / GaInP matrix)
    barrier: Material           # i-region material next to the claddings
    substrate: Material
    N_A: float = 5e17           # cm^-3 (Reischle 2008: Zn ~5e17)
    N_D: float = 5e17           # cm^-3 (Si ~5e17)
    d_i_nm: float = 120.0       # intrinsic region thickness
    d_active_nm: float = 8.0
    n_ideality: float = 2.0     # [E] SRH channel ideality
    tau_SRH_ns: float = 1.0     # [E]
    R_s_ohm: float = 50.0       # [E] small mesa incl. contacts
    area_um2: float = 0.785     # mesa / current-aperture area (1 um diameter)
    f_Rs_local: float = 0.3     # [E] share of I^2 R_s dissipated at the mesa
    T: float = 300.0            # default temperature when a method gets T=None
    # ---- transport parameters
    mu_n_cm2: float = 2000.0    # [E] cladding electron mobility
    mu_p_cm2: float = 100.0     # [E] cladding hole mobility
    tau_n_ns: float = 1.0       # [E] cladding minority-electron lifetime
    tau_p_ns: float = 1.0       # [E]
    tau_active_ns: float = 1.0  # [E] active recombination time (leakage ratio)
    d_pclad_nm: float | None = None   # p-cladding thickness: caps L_n if given
    d_nclad_nm: float | None = None
    qfl_e_meV: float = 0.0      # [E] electron QFL above E_c,active
    leak_valleys: str = "G"     # cladding valleys in the leakage sum: "G" (spec
                                # formula, materials.offsets [DR]) | "GX" | "GXL"
                                # -- see module docstring item 5 on AlInP
    B_rad_cm3s: float = 1e-10   # [E]
    eps_r: float = 12.0         # [E]
    # ---- dot-layer parameters
    tau_cap0_ps: float = 5.0    # [E] dot capture time at 1e10 cm^-2
    tau_matrix_ns: float = 1.0  # [E]

    # ------------------------------------------------------------ helpers
    def _T(self, T):
        return self.T if T is None else float(T)

    @property
    def area_cm2(self) -> float:
        return self.area_um2 * 1e-8

    def active_edges(self, T=None):
        """Strained band edges of the active layer on the substrate."""
        return strain_shifts(self.active, self.substrate, self._T(T))

    def active_gap(self, T=None) -> float:
        """Strained minimum gap of the active layer [eV] (top valence band =
        hh under compression, lh under tension)."""
        e = self.active_edges(T)
        return float(min(e.Eg_hh, e.Eg_lh))

    def ln_ni_active(self, T=None) -> float:
        return ln_n_i(self.active, self._T(T), Eg=self.active_gap(T))

    # -------------------------------------------------------- built-in V
    def vbi(self, T=None) -> float:
        """Heterojunction built-in voltage [V] from the cladding Fermi levels."""
        T = self._T(T)
        E_Fn, E_Fp = fermi_levels(self.n_cladding, self.N_D, self.p_cladding, self.N_A, T)
        return E_Fn - E_Fp

    # ------------------------------------------------------ electrostatics
    def depletion(self, V_j: float, T=None) -> Depletion:
        T = self._T(T)
        V_bi = self.vbi(T)
        eps = self.eps_r * EPS0_SI
        d_i = self.d_i_nm * 1e-9
        NA, ND = self.N_A * 1e6, self.N_D * 1e6          # m^-3
        c = V_bi - V_j
        if c <= 0.0:
            W = d_i
            return Depletion(V_j, V_bi, 0.0, 0.0, 0.0, 0.0, W * 1e9,
                             eps * self.area_um2 * 1e-12 / W * 1e12, True)
        a = (Q_SI * NA / (2.0 * eps)) * (1.0 + NA / ND)
        b = Q_SI * NA * d_i / eps
        x_p = (-b + np.sqrt(b * b + 4.0 * a * c)) / (2.0 * a)
        x_n = x_p * NA / ND
        F = Q_SI * NA * x_p / eps                          # V/m
        W = x_p + d_i + x_n
        C = eps * self.area_um2 * 1e-12 / W                # F
        return Depletion(float(V_j), float(V_bi), float(F * 1e-5),
                         float(F * self.d_active_nm * 1e-9 * 1e3), float(x_p * 1e9),
                         float(x_n * 1e9), float(W * 1e9), float(C * 1e12), False)

    # ------------------------------------------------------------- I-V
    def _diffusion(self, T):
        kT = KB_EV * T
        D_n, D_p = self.mu_n_cm2 * kT, self.mu_p_cm2 * kT        # cm^2/s
        L_n = np.sqrt(D_n * self.tau_n_ns * 1e-9)
        L_p = np.sqrt(D_p * self.tau_p_ns * 1e-9)
        if self.d_pclad_nm is not None:
            L_n = min(L_n, self.d_pclad_nm * 1e-7)
        if self.d_nclad_nm is not None:
            L_p = min(L_p, self.d_nclad_nm * 1e-7)
        return D_n, D_p, L_n, L_p

    def ln_j0(self, T=None) -> dict:
        """{channel: (ln J0 [ln A/cm^2], ideality n)} for the three channels;
        a channel switched off (tau_SRH = inf, B = 0) is omitted."""
        T = self._T(T)
        out = {}
        ln_q = np.log(Q_SI)
        ln_ni_a = self.ln_ni_active(T)
        d_a = self.d_active_nm * 1e-7
        d_b = max(self.d_i_nm - self.d_active_nm, 0.0) * 1e-7
        if np.isfinite(self.tau_SRH_ns) and self.tau_SRH_ns > 0:
            terms = [ln_ni_a + np.log(d_a)]
            if d_b > 0:
                terms.append(ln_n_i(self.barrier, T) + np.log(d_b))
            out["srh"] = (ln_q + logsumexp(terms) - np.log(2.0 * self.tau_SRH_ns * 1e-9),
                          self.n_ideality)
        if self.B_rad_cm3s > 0:
            out["rad"] = (ln_q + np.log(self.B_rad_cm3s) + np.log(d_a) + 2.0 * ln_ni_a, 1.0)
        D_n, D_p, L_n, L_p = self._diffusion(T)
        out["diff"] = (ln_q + logsumexp([
            2.0 * ln_n_i(self.p_cladding, T) + np.log(D_n / (L_n * self.N_A)),
            2.0 * ln_n_i(self.n_cladding, T) + np.log(D_p / (L_p * self.N_D))]), 1.0)
        return out

    def ln_j(self, V_j: float, T=None) -> float:
        """ln J(V_j) [A/cm^2] for V_j > 0 (log-space sum of the channels)."""
        T = self._T(T)
        kT = KB_EV * T
        terms = []
        for lnJ0, n in self.ln_j0(T).values():
            x = V_j / (n * kT)
            terms.append(lnJ0 + x + np.log(-np.expm1(-x)))
        return float(logsumexp(terms))

    def j_of_vj(self, V_j: float, T=None) -> float:
        """J(V_j) [A/cm^2], both signs of V_j."""
        T = self._T(T)
        if V_j > 0:
            return float(np.exp(self.ln_j(V_j, T)))
        kT = KB_EV * T
        return float(sum(np.exp(lnJ0) * np.expm1(V_j / (n * kT))
                         for lnJ0, n in self.ln_j0(T).values()))

    def vj_of_j(self, J: float, T=None) -> float:
        """Junction voltage [V] at current density J [A/cm^2] (brentq on ln J)."""
        T = self._T(T)
        if J <= 0:
            return 0.0
        lnJt = np.log(J)
        kT = KB_EV * T
        chans = self.ln_j0(T).values()
        V_hi = max(n * kT * (lnJt - lnJ0) for lnJ0, n in chans) + 0.5
        V_lo = 1e-9
        f = lambda V: self.ln_j(V, T) - lnJt
        if f(V_lo) >= 0:
            return V_lo
        return float(brentq(f, V_lo, V_hi, xtol=1e-12, rtol=1e-13, maxiter=200))

    def v_of_i(self, I: float, T=None) -> tuple[float, float]:
        """(V_applied, V_j) [V] at current I [A]: V_applied = V_j + I R_s."""
        V_j = self.vj_of_j(I / self.area_cm2, T)
        return V_j + I * self.R_s_ohm, V_j

    def iv(self, V_applied: float, T=None) -> tuple[float, float]:
        """(I [A], V_j [V]) at an applied voltage (solves V_j + I(V_j) R_s = V)."""
        T = self._T(T)
        if V_applied <= 0 or self.R_s_ohm == 0:
            V_j = V_applied
            return self.area_cm2 * self.j_of_vj(V_j, T), V_j
        g = lambda Vj: Vj + self.area_cm2 * self.j_of_vj(Vj, T) * self.R_s_ohm - V_applied
        V_j = float(brentq(g, 0.0, V_applied, xtol=1e-12, maxiter=200))
        return self.area_cm2 * self.j_of_vj(V_j, T), V_j

    def v_on(self, T=None, J_on: float = 10.0) -> float:
        """Turn-on voltage: V_j at J = J_on A/cm^2 (default 10)."""
        return self.vj_of_j(J_on, T)

    # ------------------------------------------------- injection efficiency
    def eta_inj(self, T=None, dE_c_override_eV: float | None = None,
                dE_v_override_eV: float | None = None) -> Leakage:
        """Injection efficiency from thermionic electron leakage over the
        p-cladding and hole leakage over the n-cladding (module docstring,
        item 5).  Overrides replace the barrier energies (sanity checks)."""
        T = self._T(T)
        kT = KB_EV * T
        D_n, D_p, L_n, L_p = self._diffusion(T)
        d_a = self.d_active_nm * 1e-7
        tau_a = self.tau_active_ns * 1e-9
        pref_e = D_n * tau_a / (L_n * d_a)
        pref_h = D_p * tau_a / (L_p * d_a)
        a = self.active_edges(T)
        Nc_a = effective_dos(self.active.p["m_e"], T)
        Nv_a = n_v(self.active, T)
        # electrons: valleys of the p-cladding (Gamma strained via offsets)
        o_p = offsets(self.active, self.p_cladding, self.substrate, T)
        vals = conduction_valleys(self.p_cladding, T, dE_c_gamma=o_p.barrier.dE_c)
        per, r_e = {}, 0.0
        for v, E_c, N_c in vals:
            if v not in self.leak_valleys:
                continue
            dE = E_c - a.E_c - self.qfl_e_meV * 1e-3
            if dE_c_override_eV is not None:
                dE = dE_c_override_eV
            contrib = (N_c / Nc_a) * pref_e * np.exp(-dE / kT)
            per[v] = (float(dE), float(contrib))
            r_e += contrib
        dE_c_eff = min(d for d, _ in per.values())
        # holes: top valence band of the active vs top band of the n-cladding
        o_n = offsets(self.active, self.n_cladding, self.substrate, T)
        E_v_top_a = max(o_n.well.E_hh, o_n.well.E_lh)
        E_v_top_b = max(o_n.barrier.E_hh, o_n.barrier.E_lh)
        dE_v = E_v_top_a - E_v_top_b
        if dE_v_override_eV is not None:
            dE_v = dE_v_override_eV
        r_h = (n_v(self.n_cladding, T) / Nv_a) * pref_h * np.exp(-dE_v / kT)
        return Leakage(T, float(dE_c_eff), float(dE_v), float(r_e), float(r_h),
                       float(1.0 / (1.0 + r_e)), float(1.0 / (1.0 + r_h)),
                       float(1.0 / (1.0 + r_e + r_h)), self.leak_valleys, per)

    # ------------------------------------------------------ dot loading
    def f_qd(self, n_dot_cm2: float) -> float:
        """Fraction of active-layer carriers captured by dots vs recombining
        in the matrix/WL: 1/(1 + tau_cap/tau_matrix), tau_cap = tau_cap0
        (1e10/n_dot).  1e10 -> 0.995, 7e8 -> 0.93."""
        tau_cap = self.tau_cap0_ps * 1e-12 * (1e10 / n_dot_cm2)
        return 1.0 / (1.0 + tau_cap / (self.tau_matrix_ns * 1e-9))

    def dot_loading(self, I_uA: float, T=None, n_dot_cm2: float = 1e10,
                    aperture_um2: float = 0.785, tau_pulse_ns: float | None = None,
                    f_capture: float | None = None, tau_rad_ns: float = 1.0,
                    leakage: Leakage | None = None) -> DotLoading:
        """Carrier budget at the dot layer (module docstring, item 6)."""
        T = self._T(T)
        lk = leakage or self.eta_inj(T)
        f_QD = self.f_qd(n_dot_cm2) if f_capture is None else float(f_capture)
        N_dots = n_dot_cm2 * aperture_um2 * 1e-8
        r_dot = f_QD * lk.eta_inj * (I_uA * 1e-6 / Q_SI) / N_dots
        mu = r_dot * tau_pulse_ns * 1e-9 if tau_pulse_ns is not None else None
        n_occ = float(r_dot * tau_rad_ns * 1e-9)
        return DotLoading(float(I_uA), T, float(lk.eta_inj), float(f_QD), float(N_dots),
                          float(r_dot), None if mu is None else float(mu), n_occ,
                          bool(n_occ > 1.0), float(f_QD * lk.eta_inj / N_dots))

    def current_for_mu(self, mu_target: float, T=None, n_dot_cm2: float = 1e10,
                       aperture_um2: float = 0.785, tau_pulse_ns: float = 1.0,
                       f_capture: float | None = None) -> float:
        """Current [uA] giving per-pulse loading mu_target (exact inverse of
        dot_loading: mu is linear in I)."""
        per_uA = self.dot_loading(1.0, T, n_dot_cm2, aperture_um2, tau_pulse_ns,
                                  f_capture).mu
        return mu_target / per_uA

    # ------------------------------------------------- injection background
    def b_e(self, I_uA: float, T=None, w_meV: float = 1.0, dE_WL_meV: float = 100.0,
            n_dot_cm2: float = 1e10, aperture_um2: float = 0.785, S_dot=1.0,
            eta_rad_matrix: float = 0.1, E_urbach_meV: float | None = None,
            tau_rad_ns: float = 1.0, f_capture: float | None = None,
            leakage: Leakage | None = None) -> Background:
        """b_e(I, T): background photons inside the X window per dot X photon
        (module docstring, definition and item 7).  S_dot: the dot's thermal
        retention factor, a float or a callable S_dot(T)."""
        T = self._T(T)
        kT_meV = 1e3 * KB_EV * T
        E_U = max(kT_meV, E_urbach_meV if E_urbach_meV is not None else kT_meV)
        ld = self.dot_loading(I_uA, T, n_dot_cm2, aperture_um2, None, f_capture,
                              tau_rad_ns, leakage)
        S = float(S_dot(T)) if callable(S_dot) else float(S_dot)
        xi = float(xi_window(w_meV, dE_WL_meV, E_U))
        rate_bg = (1.0 - ld.f_QD) * ld.eta_inj * (I_uA * 1e-6 / Q_SI) * eta_rad_matrix * xi
        rate_x = min(ld.r_dot, 1.0 / (tau_rad_ns * 1e-9)) * S
        b = rate_bg / rate_x if rate_x > 0 else np.inf
        return Background(float(b), T, I_uA, xi, E_U, ld.eta_inj, ld.f_QD, rate_bg,
                          rate_x, ld.r_dot, S, ld.saturated)

    # ------------------------------------------------------------ heating
    def junction_power(self, I_uA: float, T=None, eta_total: float = 0.01,
                       h_nu_eV: float | None = None) -> float:
        """P = I V_j + f_Rs,local I^2 R_s - (I/q) eta_total h nu  [W]."""
        T = self._T(T)
        I = I_uA * 1e-6
        _, V_j = self.v_of_i(I, T)
        hnu = self.active_gap(T) if h_nu_eV is None else h_nu_eV
        return float(I * V_j + self.f_Rs_local * I * I * self.R_s_ohm - I * eta_total * hnu)


# ============================================================ convenience API

def evaluate_injection(diode: Diode, I_uA: float, T: float, n_dot_cm2: float,
                       aperture_um2: float, tau_pulse_ns: float | None, w_meV: float,
                       dE_WL_meV: float, S_dot=1.0, tau_rad_ns: float = 1.0,
                       eta_rad_matrix: float = 0.1, E_urbach_meV: float | None = None,
                       eta_total: float = 0.01) -> InjectionResult:
    """Everything for one operating point (I, T): junction voltage, field,
    leakage, dot loading, background and heat.  The drive-block mapping is
        DriveBlock.mu          <- result.mu
        DriveBlock.eta_capture <- result.eta_capture  (per-dot thinning)
        DriveBlock.b_e/_m/_Eact <- to_background_channel(...)
        DriveBlock.V           <- result.V_j (thermal.t_junction wants P = duty * P_junction)
    """
    I = I_uA * 1e-6
    V_app, V_j = diode.v_of_i(I, T)
    dep = diode.depletion(V_j, T)
    lk = diode.eta_inj(T)
    ld = diode.dot_loading(I_uA, T, n_dot_cm2, aperture_um2, tau_pulse_ns, None,
                           tau_rad_ns, lk)
    bg = diode.b_e(I_uA, T, w_meV, dE_WL_meV, n_dot_cm2, aperture_um2, S_dot,
                   eta_rad_matrix, E_urbach_meV, tau_rad_ns, None, lk)
    P = diode.junction_power(I_uA, T, eta_total)
    return InjectionResult(I_uA, T, V_j, V_app, dep.V_bi, diode.v_on(T),
                           I / diode.area_cm2, dep, lk, ld, bg, P, ld.mu,
                           ld.eta_capture_dot, bg.b_e)


def to_background_channel(diode: Diode, I_grid_uA, T_grid, I_ref_uA: float = 1.0,
                          name: str = "injection", **b_e_kwargs):
    """Least-squares fit of the computed b_e(I, T) to the legacy form
    b = A (I/I_ref)^m exp(-E_act/kT) over the grid (log-linear in ln A, m,
    E_act).  Returns (BackgroundChannel, info) with info['max_rel_err'] the
    worst-case |fit/computed - 1| on the grid, info['drive_fields'] the
    DriveBlock entries, and info['grid'] the computed values -- so a legacy
    card is populated WITH its fit error, never silently."""
    from .spectral import KB as KB_meV
    Is = np.asarray(I_grid_uA, float)
    Ts = np.asarray(T_grid, float)
    rows, y, vals = [], [], np.zeros((len(Is), len(Ts)))
    for i, I in enumerate(Is):
        for j, T in enumerate(Ts):
            b = diode.b_e(I, T, **b_e_kwargs).b_e
            vals[i, j] = b
            rows.append([1.0, np.log(I / I_ref_uA), -1.0 / (KB_meV * T)])
            y.append(np.log(b))
    coef, *_ = np.linalg.lstsq(np.asarray(rows), np.asarray(y), rcond=None)
    lnA, m, E_act = coef
    ch = BackgroundChannel(name, float(np.exp(lnA)), float(m), float(E_act), float(I_ref_uA))
    fit = np.array([[ch.rate(I, T) for T in Ts] for I in Is])
    err = np.abs(fit / vals - 1.0)
    info = {"max_rel_err": float(err.max()), "rel_err": err, "grid": vals,
            "I_grid_uA": Is, "T_grid": Ts,
            "drive_fields": {"b_e": ch.A, "b_e_m": ch.m, "b_e_Eact": ch.E_act,
                             "I_ref_uA": ch.I_ref}}
    return ch, info


# ==================================================================== presets

def red_diode_preset(active: str = "GaAsP", **overrides) -> Diode:
    """Reischle-2008-like red p-i-n on GaAs: p-Al0.52In0.48P 5e17 / i:
    (Al0.55Ga0.45)0.51In0.49P 50 nm + active 8 nm + (Al0.55Ga0.45)InP 50 nm /
    n-Al0.52In0.48P 5e17.  active = "GaAsP" (GaAs0.6P0.4 well, 1.92 eV
    unstrained, tensile on GaAs) or "GaInP" (Ga0.51In0.49P matrix of the InP
    dots).  Reischle's 10 nm (Al0.2Ga0.8)InP dot-adjacent layers are folded
    into the active/barrier split.  p-cladding thickness 200 nm (their
    AlInP:Zn), which caps the electron diffusion length."""
    act = GaAsP(0.4) if active == "GaAsP" else GaInP()
    kw = dict(p_cladding=AlInP(0.52), n_cladding=AlInP(0.52), active=act,
              barrier=AlGaInP(0.55), substrate=binary("GaAs"), N_A=5e17, N_D=5e17,
              d_i_nm=108.0, d_active_nm=8.0, d_pclad_nm=200.0, d_nclad_nm=200.0)
    kw.update(overrides)
    return Diode(**kw)


def hkust_preset(x_p: float = 0.4, **overrides) -> Diode:
    """Gu et al. 2025 (HKUST)-like stack: Al0.4Ga0.6As claddings and barriers
    (70 nm each side of the well) around a GaAs_1-x P_x well (they use
    x = 0.35; default here 0.40 to match the materials preset).  NOTE: the
    materials model puts GaAsP(0.4)/AlGaAs(0.4) at dE_v,hh ~ 0 (type II for
    holes) -- hole leakage is then large; this is a model statement about
    the VBOs, flagged, not a design claim."""
    kw = dict(p_cladding=AlGaAs(0.4), n_cladding=AlGaAs(0.4), active=GaAsP(x_p),
              barrier=AlGaAs(0.4), substrate=binary("GaAs"), N_A=5e17, N_D=5e17,
              d_i_nm=148.0, d_active_nm=8.0)
    kw.update(overrides)
    return Diode(**kw)


def gaas_homojunction(N_A: float = 1e17, N_D: float = 1e17, **overrides) -> Diode:
    """All-GaAs p-i-n (verify reference: Sze's n_i / V_bi / Shockley diode)."""
    g = binary("GaAs")
    kw = dict(p_cladding=g, n_cladding=g, active=g, barrier=g, substrate=g,
              N_A=N_A, N_D=N_D, d_i_nm=120.0, d_active_nm=8.0)
    kw.update(overrides)
    return Diode(**kw)
