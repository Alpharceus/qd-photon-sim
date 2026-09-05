"""Module L -- Confinement levels of a self-assembled dot computed from the
ACTUAL dot / matrix / barrier materials, and their mapping onto the
two-channel Arrhenius retention of Module E (integrator.retention):

    S(T) = 1 / (1 + a_esc e^{-E_a/kT} + b_p e^{-E_b/kT}).

WHY THIS MODULE EXISTS. Until this tier every thermal-escape energy in the
simulator was an arsenide-class fit proxy (Chatzarakis InAs/GaAs numbers
re-used for InP/GaAsP and InP/GaInP dots). This module derives E_a and E_b
from the band offsets (materials.offsets, Vurgaftman 2001 model-solid VBOs +
strain), the dot geometry, and the effective masses, so that changing the
matrix from GaAs to GaAsP or AlGaInP changes the retention curve for a
stated physical reason.

MODEL (tag [E]). HONEST ACCURACY: tens of meV for a dot close to the
  lattice-matched, non-piezoelectric, single-parameter-source regime this
  separable-disk + linear-strain treatment was built for. Three classes of
  input honestly miss by 100-300 meV, documented with the numbers and the
  mechanism in verify/verify_dot_levels.py's "known deviations" table rather
  than forced to pass: (i) large-mismatch 3D islands compared to an ensemble
  PL/lasing line (Gu et al. HKUST InP/GaAsP, ~160 meV, see below); (ii)
  piezoelectric (211)B-grown dots, whose built-in field this model omits
  entirely (Chatzarakis et al. InAs/GaAs, ~90-300 meV); (iii) combining one
  paper's isolated, measured band-offset number (e.g. Pryor's InP/GaInP
  unstrained VBO) with this module's OWN bulk band gaps and linear
  deformation-potential strain, which are not the ingredients that
  paper's own (typically k.p) calculation used internally -- the isolated
  number is not separable from the rest of that calculation, so re-injecting
  it here does not reproduce that paper's own strained offset (Pryor/
  Reischle InP/GaInP, ~130-260 meV).
  Lens / truncated-pyramid dots are replaced by a SEPARABLE DISK:
    * along z: a symmetric finite square well of width = dot height, with
      BenDaniel-Duke matching (mass discontinuity) -> E_z;
    * in-plane: a finite circular well of radius R whose depth is the
      REMAINING barrier V - E_z (adiabatic decoupling, valid for flat dots,
      height << 2R; it keeps E_z + E_r < V so that an escape energy
      V - E_z - E_r is never negative) -> E_r (l = 0) and the p-shell
      (l = 1);
    * the same biaxial (pseudomorphic-film) strain as materials.offsets is
      used for the dot. A true 3D island partly relaxes, whereas the amount
      and spatial distribution of that relaxation require continuum
      elasticity plus multiband k.p (Grundmann and Stier, Phys. Rev. B 60,
      1999) [E]. This is the leading systematic of the module; no
      unmeasured relaxation fraction is fitted here;
    * no k.p band mixing, no piezoelectric field, no Coulomb correlation
      beyond a Gaussian-orbital exciton binding (Section 7 below);
    * electrons escape to the Gamma edge of the matrix / barrier even when
      that layer is X-indirect (AlGaAs x > 0.45): the X continuum lies
      lower and would make the electron escape energy SMALLER -- flagged
      in DotLevels.notes.
  Literature classes reproduced (or honestly missed) by these choices are
  listed in verify/verify_dot_levels.py with the numbers.

HKUST InP/GaAsP FINDING. Gu et al., Opt. Express 33, 23732 (2025) [V]
report 4--7 nm high InP islands and approximately 750--755 nm ensemble QD
PL/lasing in their OWN GaAs0.65P0.35 well (class_presets()
"InP/GaAs0.65P0.35/AlGaAs0.4 on GaAs"). For that composition at the
midpoint 5.5 nm height, this model gives about 1.484 eV (V_e = 175 meV,
dE_e_matrix = 105 meV); the design-card variant at the simulator's own
GaAs0.60P0.40 well [A] gives about 1.520 eV (V_e = 200 meV). Neither
difference from the ensemble line is a missing Varshni correction:
materials.bandgap already applies the Vurgaftman et al., J. Appl. Phys. 89,
5815 (2001) [V] temperature dependence at the requested T. Nor can it be
removed honestly by changing the tabulated InP electron mass. It is the
expected uncertainty of applying a fully coherent, single-band disk to a
large-mismatch 3D island, compounded by unknown alloy ordering/composition
and the fact that the measured laser line is an ensemble observable. Gu et
al. also observe both type-I and type-II alignment depending on growth, so
their ensemble wavelength cannot uniquely determine either test structure;
because a bound electron has a positive confinement energy, its escape
energy is in both cases strictly below its own V_e, not a solver defect.

UNITS. Energies in meV inside the solvers and the level table, eV for band
edges and transition energies (suffix _eV), nm for lengths, masses in m0.
hbar^2 / (2 m0) = 38.09982 meV nm^2;  e^2/(4 pi eps0) = 1439.964 meV nm.
"""
from __future__ import annotations

from dataclasses import dataclass, field, replace

import numpy as np
from scipy.optimize import brentq
from scipy.special import jv, jvp, kve

from . import materials as M

HB2_2M0 = 38.09982        # hbar^2/(2 m0)  [meV nm^2]
E2_4PIEPS0 = 1439.964     # e^2/(4 pi eps0) [meV nm]
KB_MEV = 0.08617333262    # Boltzmann constant [meV/K] (same as spectral.KB)
HC_EV_NM = 1239.84198     # h c [eV nm]

# Static dielectric constants of the binaries (Ioffe NSM / Adachi tables).
# InP 12.5, GaAs 12.9, GaP 11.1, InAs 15.15 [DR]; AlAs 10.06, AlP 9.8 [E].
# Alloys: linear in the binary fractions (Material.composition).
EPS_R_BINARY = {"InP": 12.5, "GaAs": 12.9, "GaP": 11.1, "InAs": 15.15,
                "AlAs": 10.06, "AlP": 9.8}


def eps_r_static(m: M.Material) -> float:
    """Static dielectric constant, linear alloy interpolation [DR/E]."""
    return float(sum(f * EPS_R_BINARY[b] for b, f in m.composition.items()))


# ------------------------------------------------------------------ geometry / system

@dataclass
class DotGeometry:
    """Dot dimensions (nm).

    height_nm        : dot height (z well width)
    radius_nm        : in-plane radius (disk well radius; for a lens dot use
                       ~ base radius, the ground state sits well inside it)
    wl_thickness_nm  : wetting-layer thickness; 0.5 nm ~ 2 ML of InP / InAs
    shape            : note on the real shape -> separable-disk mapping [E]:
                       lens or truncated pyramid -> disk of the same height
                       and base radius. Expected accuracy: tens of meV on
                       each confinement energy (the height sets E_z to
                       ~10 %, the in-plane shape/size sets E_r to ~30 %).
    """
    height_nm: float = 3.0
    radius_nm: float = 10.0
    wl_thickness_nm: float = 0.5
    shape: str = "lens/truncated-pyramid -> separable disk [E]"
    # Thickness of the matrix layer between the dot/WL plane and the outer
    # barrier ON EACH SIDE (nm). None = matrix semi-infinite (the plain
    # separable-disk model). A finite value makes the z-problems of the dot
    # and of the wetting layer three-region wells (dot | matrix | barrier)
    # solved by the layered BenDaniel-Duke finite-difference solver -- this
    # is what a thin GaAs spacer between an InGaAs WL and an AlGaAs / SSL
    # barrier does to E_WL (Chatzarakis class), [E].
    matrix_thickness_nm: float | None = None


@dataclass
class DotSystem:
    """A dot of `dot` material buried in `matrix` (the layer it sits in: a
    GaAsP well, GaInP, GaAs, ...), which is itself clad by `barrier`
    (AlGaAs, AlGaInP, ...), all pseudomorphic on `substrate`.

    vbo_override_eV : optional unstrained valence-band-offset overrides
        {"dot-matrix": E_v0(dot) - E_v0(matrix) [eV],
         "matrix-barrier": E_v0(matrix) - E_v0(barrier) [eV]}.
        Applied by rewriting the `vbo` parameter of the dot (resp. barrier)
        so that the UNSTRAINED offset to the matrix equals the override;
        the strain terms of materials.strain_shifts are then added on top
        unchanged (both E_c0 and E_v0 of the layer move together, so its
        gap is preserved). Example: {"dot-matrix": -0.045} = Pryor, PRB 56,
        10404 (1997) InP/GaInP VBO (InP VB 45 meV BELOW GaInP, from Fe
        impurity-level spectra) in place of the model-solid +0.17 eV.
    eps_r : static dielectric constant for the exciton binding; default =
        eps_r_static(dot).
    """
    dot: M.Material
    matrix: M.Material
    barrier: M.Material
    substrate: M.Material
    T: float = 300.0
    geometry: DotGeometry = field(default_factory=DotGeometry)
    vbo_override_eV: dict = field(default_factory=dict)
    eps_r: float | None = None
    name: str = ""


def _with_vbo(layer: M.Material, reference: M.Material, offset_eV: float | None,
              sign: float) -> M.Material:
    """Copy of `layer` whose unstrained VBO is reference.vbo + sign*offset."""
    if offset_eV is None:
        return layer
    p = dict(layer.p)
    p["vbo"] = reference.p["vbo"] + sign * float(offset_eV)
    return replace(layer, p=p, label=layer.label + "[vbo-override]")


# ------------------------------------------------------------------- 1D finite well

@dataclass
class Well1D:
    V_meV: float
    width_nm: float
    m_in: float
    m_out: float
    energies_meV: list          # all bound states found, ascending (<= 4)
    parity: list                # "even" / "odd" for each entry
    E_meV: float | None         # the requested state n (None if absent)
    rms_z_nm: float | None      # <z^2>^{1/2} of the ground state
    bound: bool


def _scan_roots(f, xs, n_max=4):
    """Sign-change scan + brentq on a grid; skips non-finite samples."""
    vals = np.asarray(f(xs), float)
    roots = []
    for i in range(len(xs) - 1):
        a, b = vals[i], vals[i + 1]
        if not (np.isfinite(a) and np.isfinite(b)):
            continue
        if a == 0.0:
            roots.append(float(xs[i]))
        elif a * b < 0.0:
            roots.append(float(brentq(f, xs[i], xs[i + 1], xtol=1e-12, rtol=1e-12)))
        if len(roots) >= n_max:
            break
    return roots


def finite_well_1d(V_meV: float, width_nm: float, m_in: float, m_out: float,
                   n: int = 0) -> Well1D:
    """Bound states of a symmetric finite square well of depth V and width w
    with BenDaniel-Duke matching (psi and psi'/m continuous):
        even:  k tan(k a) = (m_in/m_out) kappa
        odd:  -k cot(k a) = (m_in/m_out) kappa
    a = w/2, k = sqrt(2 m_in E)/hbar, kappa = sqrt(2 m_out (V - E))/hbar.
    Solved in the pole-free forms  k sin(ka) - r kappa cos(ka) = 0  and
    k cos(ka) + r kappa sin(ka) = 0 (r = m_in/m_out) on a grid uniform in k
    (spacing k a < 0.02) followed by brentq. A symmetric 1D well always
    binds one even state for V > 0; V <= 0 -> bound = False.
    The ground-state rms extent <z^2>^{1/2} is integrated numerically from
    the normalized wavefunction cos(kz) / cos(ka) e^{-kappa(|z|-a)}."""
    V = float(V_meV)
    a = 0.5 * float(width_nm)
    if V <= 0.0 or width_nm <= 0.0:
        return Well1D(V, width_nm, m_in, m_out, [], [], None, None, False)
    r = m_in / m_out

    def k_of(E):
        return np.sqrt(np.maximum(m_in * E, 0.0) / HB2_2M0)

    def q_of(E):
        return np.sqrt(np.maximum(m_out * (V - E), 0.0) / HB2_2M0)

    def f_even(E):
        k, q = k_of(E), q_of(E)
        return k * np.sin(k * a) - r * q * np.cos(k * a)

    def f_odd(E):
        k, q = k_of(E), q_of(E)
        return k * np.cos(k * a) + r * q * np.sin(k * a)

    kmax = float(k_of(V))
    N = max(4000, int(kmax * a / 0.02) + 2)
    ks = np.linspace(0.0, kmax, N)[1:-1]
    Es = HB2_2M0 * ks**2 / m_in
    Es = Es[(Es > 0) & (Es < V * (1 - 1e-9))]
    ev = _scan_roots(f_even, Es)
    od = _scan_roots(f_odd, Es)
    states = sorted([(E, "even") for E in ev] + [(E, "odd") for E in od])[:4]
    if not states:
        return Well1D(V, width_nm, m_in, m_out, [], [], None, None, False)
    energies = [s[0] for s in states]
    parity = [s[1] for s in states]
    E0 = energies[0]
    k0, q0 = float(k_of(E0)), float(q_of(E0))
    rms = None
    if q0 > 1e-6:
        zout = min(40.0 / q0, 500.0)
        z_in = np.linspace(0.0, a, 20001)
        z_out = a + np.linspace(0.0, zout, 20001)[1:]
        z = np.concatenate([z_in, z_out])
        psi = np.concatenate([np.cos(k0 * z_in),
                              np.cos(k0 * a) * np.exp(-q0 * (z_out - a))])
        w = psi**2
        rms = float(np.sqrt(np.trapezoid(z**2 * w, z) / np.trapezoid(w, z)))
    E_n = energies[n] if n < len(energies) else None
    return Well1D(V, width_nm, m_in, m_out, energies, parity, E_n, rms, True)


def well_1d_layered(widths_nm, V_meV, masses, h_nm: float = 0.01,
                    pad_nm: float = 20.0) -> Well1D:
    """Symmetric multilayer 1D well with the BenDaniel-Duke kinetic operator
    -(hbar^2/2) d/dz (1/m(z)) d/dz + V(z), solved by finite differences
    (harmonic-mean interface masses, Dirichlet walls at +-L).
    widths_nm[0] = FULL width of the central region; widths_nm[i > 0] =
    thickness of the i-th shell ON EACH SIDE; the last region (V_meV[-1],
    masses[-1]) extends pad_nm beyond the last interface. Bound states are
    those below the outer potential. Used for the dot | matrix | barrier
    z-problem when DotGeometry.matrix_thickness_nm is finite; reduces to
    finite_well_1d for a single shell of large thickness (checked in
    verify_dot_levels)."""
    from scipy.linalg import eigh_tridiagonal
    widths = [float(w) for w in widths_nm]
    edges = [0.5 * widths[0]]
    for w in widths[1:]:
        edges.append(edges[-1] + w)
    L = edges[-1] + pad_nm
    N = int(2 * L / h_nm) + 1
    if N > 40001:
        N = 40001
    z = np.linspace(-L, L, N)
    h = z[1] - z[0]
    idx = np.searchsorted(edges, np.abs(z), side="right")   # region index
    Vz = np.asarray(V_meV, float)[idx]
    mz = np.asarray(masses, float)[idx]
    inv_m_half = 2.0 / (mz[:-1] + mz[1:])
    c = HB2_2M0 * inv_m_half / h**2
    diag = np.concatenate([[0.0], c]) + np.concatenate([c, [0.0]]) + Vz
    off = -c
    E, U = eigh_tridiagonal(diag, off, select="i", select_range=(0, 3))
    V_out = float(V_meV[-1])
    bound = [(float(e), U[:, i]) for i, e in enumerate(E) if e < V_out]
    if not bound:
        return Well1D(V_out, widths[0], masses[0], masses[-1], [], [], None, None, False)
    energies = [b[0] for b in bound]
    parity = ["even" if abs(b[1][N // 2]) > 1e-3 * np.abs(b[1]).max() else "odd" for b in bound]
    u0 = bound[0][1]
    rms = float(np.sqrt(np.sum(z**2 * u0**2) / np.sum(u0**2)))
    return Well1D(V_out, widths[0], masses[0], masses[-1], energies, parity, energies[0], rms, True)


# --------------------------------------------------------------- 2D circular well

@dataclass
class Disk2D:
    V_meV: float
    radius_nm: float
    m_in: float
    m_out: float
    E0_meV: float | None        # l = 0 ground state
    E1_meV: float | None        # l = 1 ground state (p shell), None if unbound
    E_meV: float | None         # ground state of the requested l
    rms_r_nm: float | None      # <r^2>^{1/2} of the l = 0 state
    bound: bool
    p_bound: bool


def _disk_ground(V, R, m_in, m_out, l):
    """Lowest root of the BenDaniel-Duke matching for angular momentum l:
    (k/m_in) J_l'(kR) K_l(kappa R) = (kappa/m_out) K_l'(kappa R) J_l(kR),
    written with exponentially scaled K (kve) so that huge kappa R does not
    underflow; K_l' = -(K_{l-1} + K_{l+1})/2."""
    def k_of(E):
        return np.sqrt(np.maximum(m_in * E, 0.0) / HB2_2M0)

    def q_of(E):
        return np.sqrt(np.maximum(m_out * (V - E), 0.0) / HB2_2M0)

    def f(E):
        k, q = k_of(E), q_of(E)
        x, y = k * R, q * R
        Jl, Jlp = jv(l, x), jvp(l, x)
        Kl = kve(l, y)
        Klp = -0.5 * (kve(l - 1, y) + kve(l + 1, y))
        return (k / m_in) * Jlp * Kl - (q / m_out) * Klp * Jl

    kmax = float(k_of(V))
    N = max(4000, int(kmax * R / 0.02) + 2)
    ks = np.linspace(0.0, kmax, N)[1:-1]
    Es = HB2_2M0 * ks**2 / m_in
    Es = Es[(Es > 0) & (Es < V * (1 - 1e-9))]
    roots = _scan_roots(f, Es, n_max=1)
    if not roots and l == 0:
        # The l = 0 state of a 2D finite circular well is bound for EVERY
        # V > 0 (no threshold depth, unlike l >= 1): a shallow/narrow well
        # (small m_in V R^2) still binds, but with an exponentially small
        # binding energy V - E0 (kappa R << 1), which can fall below the
        # resolution of the E-grid above (itself cut off at V*(1 - 1e-9) to
        # avoid the kappa -> 0 singular point of Kl). Recover it by bracketing
        # geometrically in the residual (V - E) down toward E -> V^-, using
        # the last (still V > E) grid sample as the starting point.
        e_lo = float(Es[-1]) if len(Es) else 0.0
        f_lo = f(e_lo)
        if np.isfinite(f_lo) and f_lo != 0.0:
            for resid in V * np.logspace(-1.0, -15.0, 57):
                e_hi = V - resid
                if e_hi <= e_lo:
                    continue
                f_hi = f(e_hi)
                if not np.isfinite(f_hi):
                    continue
                if f_lo * f_hi < 0.0:
                    roots = [brentq(f, e_lo, e_hi, xtol=1e-14 * V, rtol=1e-13)]
                    break
                e_lo, f_lo = e_hi, f_hi
    return roots[0] if roots else None


def finite_disk_2d(V_meV: float, radius_nm: float, m_in: float, m_out: float,
                   l: int = 0) -> Disk2D:
    """2D circular finite well of depth V, radius R, BenDaniel-Duke matching.
    Inside J_l(k r), outside K_l(kappa r). Returns the l = 0 ground state
    (always bound for V > 0 in 2D), the l = 1 ground state (p shell, may be
    unbound -> None, p_bound = False), the ground state of the requested l,
    and the rms radius <r^2>^{1/2} of the l = 0 state integrated numerically
    from J_0(kr) / J_0(kR) K_0(kappa r)/K_0(kappa R)."""
    V, R = float(V_meV), float(radius_nm)
    if V <= 0.0 or R <= 0.0:
        return Disk2D(V, R, m_in, m_out, None, None, None, None, False, False)
    E0 = _disk_ground(V, R, m_in, m_out, 0)
    if E0 is None:
        return Disk2D(V, R, m_in, m_out, None, None, None, None, False, False)
    E1 = _disk_ground(V, R, m_in, m_out, 1)
    El = E0 if l == 0 else (E1 if l == 1 else _disk_ground(V, R, m_in, m_out, l))
    k0 = np.sqrt(m_in * E0 / HB2_2M0)
    q0 = np.sqrt(m_out * (V - E0) / HB2_2M0)
    rms = None
    if q0 > 1e-6:
        rout = min(40.0 / q0, 500.0)
        r_in = np.linspace(0.0, R, 20001)
        r_out = R + np.linspace(0.0, rout, 20001)[1:]
        r = np.concatenate([r_in, r_out])
        psi_in = jv(0, k0 * r_in)
        psi_out = jv(0, k0 * R) * kve(0, q0 * r_out) / kve(0, q0 * R) * np.exp(-q0 * (r_out - R))
        w = np.concatenate([psi_in, psi_out])**2
        rms = float(np.sqrt(np.trapezoid(r**3 * w, r) / np.trapezoid(r * w, r)))
    return Disk2D(V, R, m_in, m_out, E0, E1, El, rms, True, E1 is not None)


# ------------------------------------------------------------------------ levels

@dataclass
class DotLevels:
    """All energies in meV unless the name ends in _eV. Confinement energies
    are measured from the STRAINED dot band edges into the gap (positive =
    deeper transition). For an unbound hole (type II) E_h = V_h <= 0: the
    hole sits at the matrix VB top, which lies |V_h| above the dot HH edge."""
    system: DotSystem
    # band offsets (well depths seen by the dot carriers), meV
    V_e: float                  # dot -> matrix, conduction band
    V_h: float                  # dot -> matrix, top valence band (HH of dot vs max(HH,LH) of matrix)
    V_e_mb: float               # matrix -> barrier step, conduction band
    V_h_mb: float               # matrix -> barrier step, valence band
    Eg_dot_strained_eV: float
    Eg_matrix_eV: float
    # confinement
    E_e: float
    E_e_z: float
    E_e_r: float
    E_h: float
    E_h_z: float
    E_h_r: float
    sp_split_e: float
    sp_split_h: float
    electron_bound: bool
    hole_bound: bool
    p_bound_e: bool
    p_bound_h: bool
    type: str                   # "I" | "II"
    # exciton
    E_bind: float
    E_bind_uncapped: float
    E_X_eV: float
    lambda_nm: float
    E_WL_eV: float
    E_WL_e_z: float
    E_WL_h_z: float
    # escape energies, meV
    dE_e_matrix: float
    dE_h_matrix: float
    dE_e_barrier: float
    dE_h_barrier: float
    dE_pair_WL: float
    dE_pair_half: float
    # Gaussian sizes for qd_gf.PhononParams
    l_xy_nm: float
    l_z_nm: float
    rms: dict
    masses: dict
    eps_r: float
    notes: list


def _vtop(e: M.StrainedEdges) -> float:
    return max(e.E_hh, e.E_lh)


def _gauss_binding(l_e_xy, l_h_xy, eps_r):
    """Exciton binding of two in-plane Gaussian orbitals [E].
    Derivation: with |psi_e|^2 ~ exp(-r^2/l_e^2), |psi_h|^2 ~ exp(-r^2/l_h^2)
    (the qd_gf convention, per-axis variance l^2/2), the relative coordinate
    rho = r_e - r_h is Gaussian with per-axis variance (l_e^2 + l_h^2)/2,
    i.e. |psi_rel|^2 ~ exp(-rho^2/L^2), L^2 = l_e^2 + l_h^2. The 2D radial
    density is p(rho) = (2 rho/L^2) exp(-rho^2/L^2) and
        <1/rho> = int_0^inf (2/L^2) exp(-rho^2/L^2) drho = sqrt(pi)/L.
    Written with per-axis rms widths s = l/sqrt(2) this is the familiar
    sqrt(pi/2)/sqrt(s_e^2 + s_h^2). So
        E_bind = e^2/(4 pi eps0 eps_r) * sqrt(pi) / sqrt(l_e^2 + l_h^2).
    Rydberg check: for l -> large the true 2D exciton binds with 4 Ry* on the
    Bohr radius, which this first-order (frozen-orbital) expectation value
    underestimates; for strongly confined dots (l < a_B ~ 10 nm) the
    frozen-orbital estimate is the standard one and adequate at the
    10 meV level. The z extent is neglected (overestimates E_bind slightly)."""
    L = np.sqrt(l_e_xy**2 + l_h_xy**2)
    return float(np.sqrt(np.pi) * E2_4PIEPS0 / (eps_r * L))


def levels(system: DotSystem) -> DotLevels:
    """Solve the separable-disk levels for a DotSystem. See module docstring.

    A zero-argument DotSystem constructor is also accepted for convenience,
    so ``levels(class_presets()[name])`` is equivalent to
    ``levels(class_presets()[name]())`` [A].

    Gaussian-size mapping for qd_gf.PhononParams: the exciton form factor of
    qd_gf, F2(E) = <exp(-(E/hbar)^2 (l_xy^2 sin^2 + l_z^2 cos^2)/(2 c_s^2))>,
    is the squared Fourier transform of |psi|^2 ~ exp(-r^2/l_xy^2 - z^2/l_z^2)
    (FT of exp(-x^2/l^2) is ~ exp(-q^2 l^2/4); squared -> exp(-q^2 l^2/2)).
    For that density  <x^2> = l_xy^2/2 per axis, hence
        <r^2>_in-plane = l_xy^2      ->  l_xy = <r^2>^{1/2},
        <z^2>          = l_z^2 / 2   ->  l_z  = sqrt(2) <z^2>^{1/2}.
    The single exciton form factor is fed the e/h AVERAGE of <r^2> and
    <z^2> (the standard same-form-factor IBM approximation, [E])."""
    S = system() if callable(system) else system
    g = S.geometry
    T = float(S.T)
    notes = []
    sub = S.substrate
    dot = _with_vbo(S.dot, S.matrix, S.vbo_override_eV.get("dot-matrix"), +1.0)
    barrier = _with_vbo(S.barrier, S.matrix, S.vbo_override_eV.get("matrix-barrier"), -1.0)
    if "dot-matrix" in S.vbo_override_eV:
        notes.append(f"unstrained dot-matrix VBO overridden to {S.vbo_override_eV['dot-matrix']*1e3:+.0f} meV "
                     f"(model-solid would be {(S.dot.p['vbo'] - S.matrix.p['vbo'])*1e3:+.0f} meV)")

    om = M.offsets(dot, S.matrix, sub, T)          # dot in matrix
    ob = M.offsets(S.matrix, barrier, sub, T)      # matrix in barrier
    V_e = 1e3 * om.dE_c
    V_h = 1e3 * (om.well.E_hh - _vtop(om.barrier))
    V_e_mb = 1e3 * ob.dE_c
    V_h_mb = 1e3 * (_vtop(ob.well) - _vtop(ob.barrier))
    if om.barrier.E_lh > om.barrier.E_hh + 1e-6:
        notes.append("matrix is tensile: its LH edge is the top VB, hole depth measured to LH continuum")
    if not M.is_direct(S.matrix, T):
        notes.append(f"matrix {S.matrix.label} is indirect: electron escape to its X/L valley would be "
                     f"{1e3*(M.bandgap(S.matrix,T,'G') - M.bandgap(S.matrix,T,'min')):.0f} meV smaller than quoted")
    if not M.is_direct(barrier, T):
        notes.append(f"barrier {S.barrier.label} is indirect: Gamma-edge barrier quoted; X continuum is "
                     f"{1e3*(M.bandgap(barrier,T,'G') - M.bandgap(barrier,T,'min')):.0f} meV lower")
    Eg_dot = om.well.Eg_hh
    Eg_mat = M.bandgap(S.matrix, T, "G")

    m_e_d, m_e_m, m_e_b = dot.p["m_e"], S.matrix.p["m_e"], barrier.p["m_e"]
    m_hz_d, m_hz_m, m_hz_b = M.hh_mass_z(dot), M.hh_mass_z(S.matrix), M.hh_mass_z(barrier)
    m_hr_d, m_hr_m = M.hh_mass_inplane(dot), M.hh_mass_inplane(S.matrix)
    masses = dict(m_e_dot=m_e_d, m_e_matrix=m_e_m, m_hz_dot=m_hz_d, m_hz_matrix=m_hz_m,
                  m_hr_dot=m_hr_d, m_hr_matrix=m_hr_m, m_e_barrier=m_e_b, m_hz_barrier=m_hz_b)
    d_mat = g.matrix_thickness_nm

    def zwell(V, V_mb, width, m_d, m_m, m_b):
        """z-problem: plain finite well (semi-infinite matrix) or the layered
        dot | matrix(d) | barrier well when matrix_thickness_nm is set."""
        if V <= 0.0:
            return finite_well_1d(V, width, m_d, m_m)
        if d_mat is None:
            return finite_well_1d(V, width, m_d, m_m)
        return well_1d_layered([width, d_mat], [0.0, V, V + V_mb], [m_d, m_m, m_b])

    if d_mat is not None:
        notes.append(f"matrix layer {d_mat} nm thick on each side of the dot/WL: z-levels from the "
                     f"layered dot|matrix|barrier well [E]")

    # ---- electron
    ez = zwell(V_e, V_e_mb, g.height_nm, m_e_d, m_e_m, m_e_b)
    if ez.bound:
        E_e_z = ez.energies_meV[0]
        er = finite_disk_2d(V_e - E_e_z, g.radius_nm, m_e_d, m_e_m)
        E_e_r = er.E0_meV if er.bound else 0.0
        if not er.bound:
            # V_e - E_e_z <= 0: the layered z-solution already uses up (or
            # exceeds) the dot-matrix step, e.g. when it is really a
            # matrix-layer state riding on the extra matrix-barrier step
            # V_e_mb (matrix_thickness_nm finite). No residual in-plane well
            # is left for the disk problem, so there is no p-shell either.
            sp_e = 0.0
        else:
            sp_e = (er.E1_meV - er.E0_meV) if er.p_bound else (er.V_meV - er.E0_meV)
        rms_e_r = er.rms_r_nm if (er.bound and er.rms_r_nm) else 1.5 * g.radius_nm
        rms_e_z = ez.rms_z_nm
        e_bound, p_e = True, er.p_bound
    else:
        E_e_z, E_e_r, sp_e = V_e, 0.0, 0.0
        rms_e_r, rms_e_z = 1.5 * g.radius_nm, g.height_nm
        e_bound, p_e = False, False
        notes.append("electron NOT bound to the dot (V_e <= 0)")
    E_e = E_e_z + E_e_r

    # ---- hole
    hz = zwell(V_h, V_h_mb, g.height_nm, m_hz_d, m_hz_m, m_hz_b)
    if hz.bound:
        E_h_z = hz.energies_meV[0]
        hr = finite_disk_2d(V_h - E_h_z, g.radius_nm, m_hr_d, m_hr_m)
        E_h_r = hr.E0_meV if hr.bound else 0.0
        if not hr.bound:
            sp_h = 0.0    # see the matching electron comment above
        else:
            sp_h = (hr.E1_meV - hr.E0_meV) if hr.p_bound else (hr.V_meV - hr.E0_meV)
        rms_h_r = hr.rms_r_nm if (hr.bound and hr.rms_r_nm) else 1.5 * g.radius_nm
        rms_h_z = hz.rms_z_nm
        h_bound, p_h = True, hr.p_bound
    else:
        E_h_z, E_h_r, sp_h = V_h, 0.0, 0.0
        rms_h_r, rms_h_z = 2.0 * g.radius_nm, g.height_nm
        h_bound, p_h = False, False
        notes.append("hole NOT bound to the dot (V_h <= 0): type II, hole at the matrix VB top; "
                     "l_h in-plane set to 2R [E]")
    E_h = E_h_z + E_h_r
    typ = "I" if (e_bound and h_bound) else "II"

    # ---- Gaussian sizes (mapping in the docstring)
    l_e_xy, l_h_xy = rms_e_r, rms_h_r
    l_xy = float(np.sqrt(0.5 * (rms_e_r**2 + rms_h_r**2)))
    l_z = float(np.sqrt(rms_e_z**2 + rms_h_z**2))   # sqrt(2) * rms of the average <z^2>

    # ---- exciton
    eps_r = S.eps_r if S.eps_r is not None else eps_r_static(S.dot)
    Eb0 = _gauss_binding(l_e_xy, l_h_xy, eps_r)
    # No cap is applied: an earlier version capped Eb0 at half the summed
    # in-plane confinement (0.5*(E_e_r+E_h_r)), with no citation, and it
    # fired on 3 of 8 presets, e.g. shifting the HKUST binding 22.0 ->
    # 10.9 meV. The frozen-orbital Gaussian estimate _gauss_binding already
    # has a stated accuracy (see its docstring) and is used uncapped.
    Eb = Eb0
    E_X = Eg_dot + (E_e + E_h - Eb) * 1e-3
    lam = HC_EV_NM / E_X

    # ---- wetting layer: 1D well of the dot material in the matrix, in-plane free
    wz_e = zwell(V_e, V_e_mb, g.wl_thickness_nm, m_e_d, m_e_m, m_e_b)
    wz_h = zwell(V_h, V_h_mb, g.wl_thickness_nm, m_hz_d, m_hz_m, m_hz_b)
    E_WL_e = wz_e.energies_meV[0] if wz_e.bound else V_e
    E_WL_h = wz_h.energies_meV[0] if wz_h.bound else V_h
    E_WL = Eg_dot + (E_WL_e + E_WL_h - 10.0) * 1e-3   # 10 meV 2D-QW exciton binding [E]

    dE_e_matrix = V_e - E_e
    dE_h_matrix = V_h - E_h if h_bound else 0.0
    dE_e_barrier = V_e + V_e_mb - E_e
    dE_h_barrier = (V_h if h_bound else 0.0) + V_h_mb - (E_h if h_bound else 0.0)
    dE_pair_WL = (E_WL - E_X) * 1e3
    dE_pair_half = 0.5 * dE_pair_WL

    rms = dict(e_r=rms_e_r, e_z=rms_e_z, h_r=rms_h_r, h_z=rms_h_z)
    return DotLevels(S, V_e, V_h, V_e_mb, V_h_mb, Eg_dot, Eg_mat,
                     E_e, E_e_z, E_e_r, E_h, E_h_z, E_h_r, sp_e, sp_h,
                     e_bound, h_bound, p_e, p_h, typ,
                     Eb, Eb0, E_X, lam, E_WL, E_WL_e, E_WL_h,
                     dE_e_matrix, dE_h_matrix, dE_e_barrier, dE_h_barrier,
                     dE_pair_WL, dE_pair_half, l_xy, l_z, rms, masses, eps_r, notes)


# ------------------------------------------------------------- retention parameters

def retention_params(lv: DotLevels, tau_rad_ns: float, channel: str = "pair_half",
                     T_ref: float = 300.0, n_dot_cm2: float = 1e10,
                     tau_cap_ps: float = 10.0, verbose: bool = True) -> dict:
    """Map the level table onto integrator.retention(T, a_esc, E_a, b_p, E_b).

    E_a by channel:
      "pair_half" : 0.5 (E_WL - E_QD), correlated e-h pair escape, Gelinas et
                    al. arXiv:0910.0480 (their Ea ~ 1/2 Delta E_WL)   [DR]
      "pair"      : E_WL - E_QD (exciton escape to the WL, the Chatzarakis
                    Delta E identification)                            [DR]
      "electron"  : V_e - E_e (single electron to the matrix continuum)
      "hole"      : V_h - E_h (single hole to the matrix continuum)
      "min"       : the smallest of electron / hole / pair_half.
    a_esc from detailed balance for escape into the 2D (WL/matrix) continuum:
      a_esc = tau_rad * nu_esc0,   nu_esc0 = (1/tau_cap) * N2D / N_dot,
      N2D / N_dot = [m* kT / (pi hbar^2)] / n_dot  = number of 2D states
      within kT per dot (spin included; m* = exciton translational in-plane
      mass m_e + m_hh,xy of the matrix for the pair channels, the carrier
      mass for single-carrier channels) and tau_cap = 10 ps the capture time
      [E]. Typical result 1e3-1e6; the 1e7-1e9 fitted by Chatzarakis et al.
      is larger because a fitted Arrhenius prefactor absorbs the
      multiplicity of final states above the nominal E_a (higher WL
      subbands, barrier continuum) and the WL-transport / recapture
      dynamics (Gelinas): the fit E_a is then effectively degenerate with
      a_esc, which is why the two are not separately meaningful.
    E_b = the smaller of the two carriers' s-p splitting (p-shell channel),
      counting only carriers that are actually bound with a meaningful
      p-shell energy (sp_split_e/h is 0 exactly when that carrier is not
      bound to the dot at all, or -- for a finite matrix_thickness_nm --
      when its z-solution already exceeds the dot-matrix step and no
      residual in-plane well is left; such a carrier contributes no
      p-channel). If NEITHER carrier has a meaningful p-shell, the p-channel
      is absent (b_p = 0, E_b = 0) rather than a spurious E_b = 0 with
      b_p != 0, which would make the p-term T-independent and floor S(T)
      away from 1 even at T = 0.
      b_p = 100, carried directly as the Chatzarakis et al., Phys. Rev.
      Applied 20, 034011 (2023) Fig. 2a Arrhenius-fit class value (b ~ 100,
      Eb ~ 35 +/- 5 meV) [E]; an earlier version derived b_p = 4 from "two p
      states x an s/p degeneracy ratio of 2", a self-referential argument
      with no independent citation, 25x below the fitted class value. As
      with a_esc (see above), the fitted b likely also absorbs additional
      physics (multiple p-sublevels, phonon-assisted capture) beyond the
      bare level degeneracy, so it is carried as-is rather than re-derived.
    All [E] tags are repeated in the returned note."""
    kT = KB_MEV * T_ref
    ms = lv.masses
    if channel in ("pair_half", "pair", "min"):
        m_star = ms["m_e_matrix"] + ms["m_hr_matrix"]
    elif channel == "electron":
        m_star = ms["m_e_matrix"]
    elif channel == "hole":
        m_star = ms["m_hr_matrix"]
    else:
        raise ValueError(f"unknown channel {channel!r}")
    Ea_table = {"pair_half": lv.dE_pair_half, "pair": lv.dE_pair_WL,
                "electron": lv.dE_e_matrix, "hole": lv.dE_h_matrix if lv.hole_bound else 0.0}
    if channel == "min":
        cands = {k: v for k, v in Ea_table.items() if k in ("electron", "hole", "pair_half")}
        ch_min = min(cands, key=cands.get)
        E_a = cands[ch_min]
        ch_label = f"min({ch_min})"
        m_star = {"electron": ms["m_e_matrix"], "hole": ms["m_hr_matrix"]}.get(ch_min, m_star)
    else:
        E_a = Ea_table[channel]
        ch_label = channel
    # 2D DOS m/(pi hbar^2) x kT  [1/nm^2]  with hbar^2/m0 = 2*HB2_2M0
    n2d_nm2 = m_star * kT / (np.pi * 2.0 * HB2_2M0)
    n2d_cm2 = n2d_nm2 * 1e14
    states_per_dot = n2d_cm2 / n_dot_cm2
    nu_esc0_ps = states_per_dot / tau_cap_ps          # 1/ps
    a_esc = tau_rad_ns * 1e3 * nu_esc0_ps
    p_splits = [s for s in (lv.sp_split_e, lv.sp_split_h) if s > 0.0]
    if p_splits:
        E_b = min(p_splits)
        b_p = 100.0
    else:
        E_b = 0.0
        b_p = 0.0    # no carrier has a meaningful p-shell: the p-channel is absent, not T-independent
    note = (f"E_a = {E_a:.1f} meV ({ch_label}); a_esc = tau_rad ({tau_rad_ns} ns) x nu_esc0 "
            f"({nu_esc0_ps:.3g}/ps) with N2D/N_dot = {states_per_dot:.1f} states per dot within kT "
            f"at {T_ref:.0f} K (m* = {m_star:.3f}, n_dot = {n_dot_cm2:.1e} cm^-2), tau_cap = {tau_cap_ps} ps [E]; "
            f"E_b = min bound-carrier s-p splitting = {E_b:.1f} meV, b_p = {b_p:g} "
            f"[E, Chatzarakis et al. PRA 20, 034011 (2023) Fig. 2a fit class]; separable-disk levels [E]; "
            f"biaxial-film strain in the dot [E]; 2D-QW WL exciton binding 10 meV [E]; "
            f"Gaussian-orbital exciton binding [E]")
    if verbose:
        print(f"  retention_params[{lv.system.name or 'dot'}]: channel={ch_label} E_a={E_a:.1f} meV "
              f"a_esc={a_esc:.3g} (nu_esc0={nu_esc0_ps:.3g}/ps, {states_per_dot:.1f} states/dot) "
              f"E_b={E_b:.1f} meV b_p={b_p:g}  [Chatzarakis fit class: a_esc 1e7-1e9, "
              f"difference absorbed by E_a degeneracy]")
    return dict(E_a=float(E_a), a_esc=float(a_esc), E_b=float(E_b), b_p=b_p, channel=ch_label,
                nu_esc0_per_ps=float(nu_esc0_ps), states_per_dot=float(states_per_dot),
                m_star=float(m_star), note=note)


# ------------------------------------------------------------------------ presets

def _quaternary_Q(lambda_um: float, T: float = 300.0, label: str | None = None) -> M.Material:
    """InGaAsP lattice-matched to InP with gap hc/lambda at T, approximated
    as the linear mix of In0.532Ga0.468As and InP (both lattice-matched to
    InP, so the mix is too); the weight is solved so that bandgap(., T) hits
    the target. Quaternary bowing is neglected [E] -- adequate for the
    Q1.15 barrier used as a dot matrix (its role is a ~0.3 eV step)."""
    Eg_target = HC_EV_NM / (lambda_um * 1e3)
    t_ingaas, t_inp = M.InGaAs(0.532), M.binary("InP")
    E1, E2 = M.bandgap(t_ingaas, T), M.bandgap(t_inp, T)
    w = (E2 - Eg_target) / (E2 - E1)
    w = float(np.clip(w, 0.0, 1.0))
    q = M.quaternary_from_ternaries([t_ingaas, t_inp], [w, 1.0 - w],
                                    label or f"InGaAsP-Q{lambda_um:.2f}")
    q.tag = "E"
    return q


def class_presets() -> dict:
    """Constructors for the literature dot classes. Each entry is a callable
    f(T=..., height_nm=..., radius_nm=...) -> DotSystem. Geometries are the
    stated class values; heights from the papers where given.
      "InP/GaAs0.65P0.35/AlGaAs0.4 on GaAs" HKUST, Gu et al. Opt. Express 33,
            23732 (2025) [V]: the paper's OWN composition, GaAs0.65P0.35, and
            InP dot heights 4-7 nm (default 5.5 nm, the midpoint) [V]; radius
            not given by the paper, carried over at 12 nm [A]. Use this
            fixture, not the 0.40 variant below, for literature comparisons.
      "InP/GaAsP0.4/AlGaAs0.4 on GaAs"      the same HKUST class at the
            simulator's design-card composition GaAs0.60P0.40 [A] (chosen to
            keep the well direct-gap with margin, materials_research.md
            Section 0(c)) and radius 12 nm [A]; default height 4 nm.
            Gu et al. give no single-photon/single-dot data for either
            composition -- both are ensemble laser material [V].
      "InP/GaInP/AlGaInP0.55 on GaAs"       Reischle class (Stuttgart): dot
            InP 3 nm x R 10 nm in Ga0.51In0.49P, (Al0.55Ga0.45)0.51In0.49P
            cladding; dot-matrix VBO -45 meV (Pryor 1997).
      "InP/AlGaInP0.2/AlGaInP0.55 on GaAs"  Bommer class: same dot in
            (Al0.2Ga0.8)InP; Pryor VBO transferred through the model-solid
            GaInP -> AlGaInP0.2 step.
      "InGaAs0.5/GaAs/AlGaAs0.57 on GaAs"   Chatzarakis class: In0.5Ga0.5As
            3 nm x R 12 nm in GaAs with Al0.57Ga0.43As SSL barriers.
      "InGaAs0.5/GaAs/GaAs on GaAs"         same dot, plain GaAs barriers.
      "InGaAs0.5/AlGaAs0.57/AlGaAs0.57 on GaAs"  the SSL barrier as the
            matrix (the WL sits in AlGaAs).
      "InAs/InP/InP"                        telecom: InAs 2.5 nm x R 15 nm.
      "InAs/InGaAsP-Q1.15/InP"              InAs in a Q1.15 quaternary on InP
            (quaternary approximated, see _quaternary_Q) [E].
    """
    InP, GaAs = M.binary("InP"), M.binary("GaAs")

    def hkust_gu(T=300.0, height_nm=5.5, radius_nm=12.0):
        return DotSystem(InP, M.GaAsP(0.35), M.AlGaAs(0.4), GaAs, T,
                         DotGeometry(height_nm, radius_nm), name="InP/GaAs0.65P0.35/AlGaAs0.4 on GaAs")

    def hkust(T=300.0, height_nm=4.0, radius_nm=12.0):
        return DotSystem(InP, M.GaAsP(0.4), M.AlGaAs(0.4), GaAs, T,
                         DotGeometry(height_nm, radius_nm), name="InP/GaAsP0.4/AlGaAs0.4 on GaAs")

    def reischle(T=5.0, height_nm=3.0, radius_nm=10.0):
        return DotSystem(InP, M.GaInP(0.51), M.AlGaInP(0.55), GaAs, T,
                         DotGeometry(height_nm, radius_nm), {"dot-matrix": -0.045},
                         name="InP/GaInP/AlGaInP0.55 on GaAs")

    def bommer(T=5.0, height_nm=3.0, radius_nm=10.0):
        gainp, mat = M.GaInP(0.51), M.AlGaInP(0.2)
        vbo = -0.045 + (gainp.p["vbo"] - mat.p["vbo"])   # Pryor + model-solid Al step
        return DotSystem(InP, mat, M.AlGaInP(0.55), GaAs, T,
                         DotGeometry(height_nm, radius_nm), {"dot-matrix": vbo},
                         name="InP/AlGaInP0.2/AlGaInP0.55 on GaAs")

    # Chatzarakis et al. give the wetting layer of this class as 0.8 nm
    # (Section 2: "Intermixed: modelled as In0.5Ga0.5As, 3 nm QD, 0.8 nm
    # WL"), not the module's generic 0.5 nm (~2 ML InP/InAs) default [V].
    CHATZ_WL_NM = 0.8

    def chatz(T=78.0, height_nm=3.0, radius_nm=12.0):
        return DotSystem(M.InGaAs(0.5), GaAs, M.AlGaAs(0.57), GaAs, T,
                         DotGeometry(height_nm, radius_nm, CHATZ_WL_NM),
                         name="InGaAs0.5/GaAs/AlGaAs0.57 on GaAs")

    def chatz_gaas(T=78.0, height_nm=3.0, radius_nm=12.0):
        return DotSystem(M.InGaAs(0.5), GaAs, GaAs, GaAs, T,
                         DotGeometry(height_nm, radius_nm, CHATZ_WL_NM),
                         name="InGaAs0.5/GaAs/GaAs on GaAs")

    def chatz_ssl(T=78.0, height_nm=3.0, radius_nm=12.0):
        return DotSystem(M.InGaAs(0.5), M.AlGaAs(0.57), M.AlGaAs(0.57), GaAs, T,
                         DotGeometry(height_nm, radius_nm, CHATZ_WL_NM),
                         name="InGaAs0.5/AlGaAs0.57/AlGaAs0.57 on GaAs")

    def telecom(T=300.0, height_nm=2.5, radius_nm=15.0):
        return DotSystem(M.binary("InAs"), InP, InP, InP, T,
                         DotGeometry(height_nm, radius_nm), name="InAs/InP/InP")

    def telecom_q(T=300.0, height_nm=2.5, radius_nm=15.0):
        return DotSystem(M.binary("InAs"), _quaternary_Q(1.15, T), InP, InP, T,
                         DotGeometry(height_nm, radius_nm), name="InAs/InGaAsP-Q1.15/InP")

    return {
        "InP/GaAs0.65P0.35/AlGaAs0.4 on GaAs": hkust_gu,
        "InP/GaAsP0.4/AlGaAs0.4 on GaAs": hkust,
        "InP/GaInP/AlGaInP0.55 on GaAs": reischle,
        "InP/AlGaInP0.2/AlGaInP0.55 on GaAs": bommer,
        "InGaAs0.5/GaAs/AlGaAs0.57 on GaAs": chatz,
        "InGaAs0.5/GaAs/GaAs on GaAs": chatz_gaas,
        "InGaAs0.5/AlGaAs0.57/AlGaAs0.57 on GaAs": chatz_ssl,
        "InAs/InP/InP": telecom,
        "InAs/InGaAsP-Q1.15/InP": telecom_q,
    }


def report(lv: DotLevels) -> str:
    """One-block human-readable summary."""
    S = lv.system
    g = S.geometry
    lines = [
        f"{S.name or 'dot'}: {S.dot.label} in {S.matrix.label} / {S.barrier.label} on {S.substrate.label}, "
        f"T = {S.T:.0f} K, h = {g.height_nm} nm, R = {g.radius_nm} nm, WL = {g.wl_thickness_nm} nm",
        f"  offsets [meV]: V_e = {lv.V_e:.0f}, V_h = {lv.V_h:.0f} (dot->matrix); "
        f"V_e_mb = {lv.V_e_mb:.0f}, V_h_mb = {lv.V_h_mb:.0f} (matrix->barrier); type {lv.type}",
        f"  Eg_dot(strained) = {lv.Eg_dot_strained_eV:.3f} eV, Eg_matrix = {lv.Eg_matrix_eV:.3f} eV",
        f"  E_e = {lv.E_e:.0f} (z {lv.E_e_z:.0f} + r {lv.E_e_r:.0f}), s-p {lv.sp_split_e:.0f} meV"
        f"{'' if lv.p_bound_e else ' (p unbound: continuum)'};  "
        f"E_h = {lv.E_h:.0f} (z {lv.E_h_z:.0f} + r {lv.E_h_r:.0f}), s-p {lv.sp_split_h:.0f} meV"
        f"{'' if lv.p_bound_h else ' (p unbound)'}",
        f"  E_bind = {lv.E_bind:.1f} meV (Gaussian {lv.E_bind_uncapped:.1f}, eps_r {lv.eps_r:.1f}); "
        f"E_X = {lv.E_X_eV:.3f} eV ({lv.lambda_nm:.0f} nm); E_WL = {lv.E_WL_eV:.3f} eV",
        f"  escape [meV]: dE_e_matrix = {lv.dE_e_matrix:.0f}, dE_h_matrix = {lv.dE_h_matrix:.0f}, "
        f"dE_e_barrier = {lv.dE_e_barrier:.0f}, dE_h_barrier = {lv.dE_h_barrier:.0f}, "
        f"dE_pair_WL = {lv.dE_pair_WL:.0f}, dE_pair_half = {lv.dE_pair_half:.0f}",
        f"  Gaussian sizes: l_xy = {lv.l_xy_nm:.2f} nm, l_z = {lv.l_z_nm:.2f} nm "
        f"(rms e: r {lv.rms['e_r']:.2f} z {lv.rms['e_z']:.2f}; h: r {lv.rms['h_r']:.2f} z {lv.rms['h_z']:.2f})",
    ]
    for n in lv.notes:
        lines.append(f"  note: {n}")
    return "\n".join(lines)
