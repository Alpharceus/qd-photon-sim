"""Module B-analytic -- 1-D waveguide/DBR semi-analytic tier (Phase T3).

Millisecond-fast transfer-matrix physics under the cavity card: DBR mirror
spectra, stopband, penetration depth, cavity mode position/linewidth/Q,
dE_cav/dT from layer dn/dT, a planar (1-D) Purcell estimate, and the
inverse direction (kappa target -> pair counts).

SCOPE WALL (stated up front, per the T3 plan): this tier is 1-D transfer
matrices at normal incidence for lossless, isotropic, nonmagnetic layers
with real refractive indices. It knows NOTHING of lateral confinement,
sidewall loss, absorption, or the real collection efficiency G into an
objective -- micropillar/CBG numbers stay COMSOL/MEEP property (SIM-B2).
Its job is pre-screening and analytic inversion. Every number produced
here carries the tag [E->analytic-1D, confirm: SIM-B]; the first COMSOL
SIM-B1 result becomes a permanent cross-check the day it lands. The
planar_purcell() output in particular is a 1-D LDOS estimate, NOT a 3-D
V_eff / Purcell factor (energy-normalized conventions only; no
participation-ratio proxies are used anywhere in this module).

Conventions (uniform through the module):
  * Units: energies meV, lengths nm, temperatures K.
    E[meV] = 1.239841984e6 / lambda[nm]  (h*c = 1239.841984 eV nm).
  * Phasor convention exp(-i*omega*t); forward wave exp(+i*k*z).
  * Field amplitude vectors (E_forward, E_backward). The transfer matrix M
    maps OUTPUT-side amplitudes to INPUT-side amplitudes,
        (a_in, b_in)^T = M (a_out, b_out)^T,
    built as M = D(n_in,n_1) P_1 D(n_1,n_2) P_2 ... P_N D(n_N,n_out) with
    interface matrices D(n_i,n_j) = (1/t_ij) [[1, r_ij], [r_ij, 1]],
    r_ij = (n_i-n_j)/(n_i+n_j), t_ij = 2 n_i/(n_i+n_j) (Fresnel, normal
    incidence), and propagation P = diag(exp(-i delta), exp(+i delta)),
    delta = 2 pi n d / lambda. Then r = M21/M11, t = 1/M11, and the power
    coefficients are R = |r|^2, T = (n_out/n_in) |t|^2 with R + T = 1 for
    lossless stacks (Poynting-flux normalization).
  * Textbook provenance: transfer-matrix thin-film optics per Born & Wolf,
    *Principles of Optics* ch. 1.6 / Hecht *Optics*; Bragg-mirror closed
    forms per Yariv & Yeh, *Optical Waves in Crystals*; penetration-depth
    phase-slope convention per Babic & Corzine, IEEE JQE 28, 514 (1992)
    class of definitions (exact factor derived below in penetration_depth).

Every closed form here is exercised against an independent second method
in verify/verify_dbr.py (Airy multiple-beam summation, direct 1-D
Helmholtz ODE integration, quarter-wave admittance closed forms, numeric
phase slopes). This module is fresh code and must never import the
quarantined S-1 test oracle package (oracle-only rule).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# h*c in meV*nm: E[meV] = HC_MEV_NM / lambda[nm]. (1239.841984 eV nm.)
HC_MEV_NM = 1.239841984e6


def energy_meV(lambda_nm):
    """Photon energy in meV from vacuum wavelength in nm."""
    return HC_MEV_NM / np.asarray(lambda_nm, dtype=float)


def wavelength_nm(E_meV):
    """Vacuum wavelength in nm from photon energy in meV."""
    return HC_MEV_NM / np.asarray(E_meV, dtype=float)


# ------------------------------------------------------------------- layers

@dataclass(frozen=True)
class Layer:
    """One homogeneous dielectric layer.

    n     : real refractive index at the reference temperature T_ref
            (lossless tier: complex/absorbing indices are out of scope,
            see the module scope wall).
    d_nm  : physical thickness in nm.
    dn_dT : thermo-optic coefficient dn/dT in 1/K (linear model
            n(T) = n + dn_dT*(T - T_ref); thermal expansion of d is
            neglected -- it is an order of magnitude weaker than dn/dT
            for the III-V pairs this tier targets).
    """
    n: float
    d_nm: float
    dn_dT: float = 0.0

    def n_at(self, T_K: float, T_ref: float = 300.0) -> float:
        """Index at temperature T_K under the linear thermo-optic model."""
        return self.n + self.dn_dT * (T_K - T_ref)

    def at_T(self, T_K: float, T_ref: float = 300.0) -> "Layer":
        """Copy of this layer with n evaluated at T_K."""
        return Layer(self.n_at(T_K, T_ref), self.d_nm, self.dn_dT)


# ---------------------------------------------------------- transfer matrix

def _m_elems(layers, lambda_nm, n_in, n_out):
    """Transfer-matrix elements (m11, m12, m21, m22) for scalar or array
    lambda_nm, vectorized over wavelength. Internal kernel; see the module
    docstring for the exact convention."""
    lam = np.asarray(lambda_nm, dtype=float)
    if np.any(lam <= 0):
        raise ValueError("wavelength must be positive")
    one = np.ones_like(lam, dtype=complex)
    m11, m12, m21, m22 = one.copy(), 0 * one, 0 * one, one.copy()
    n_prev = float(n_in)
    for L in layers:
        n_j = float(L.n)
        if n_j <= 0:
            raise ValueError("refractive indices must be positive")
        r = (n_prev - n_j) / (n_prev + n_j)
        t = 2.0 * n_prev / (n_prev + n_j)
        # M <- M @ D, D = (1/t) [[1, r], [r, 1]]
        m11, m12 = (m11 + m12 * r) / t, (m11 * r + m12) / t
        m21, m22 = (m21 + m22 * r) / t, (m21 * r + m22) / t
        # M <- M @ P, P = diag(exp(-i delta), exp(+i delta))
        delta = 2.0 * np.pi * n_j * L.d_nm / lam
        pf, pb = np.exp(-1j * delta), np.exp(+1j * delta)
        m11, m12 = m11 * pf, m12 * pb
        m21, m22 = m21 * pf, m22 * pb
        n_prev = n_j
    n_o = float(n_out)
    r = (n_prev - n_o) / (n_prev + n_o)
    t = 2.0 * n_prev / (n_prev + n_o)
    m11, m12 = (m11 + m12 * r) / t, (m11 * r + m12) / t
    m21, m22 = (m21 + m22 * r) / t, (m21 * r + m22) / t
    return m11, m12, m21, m22


def transfer_matrix(layers, lambda_nm, n_in=1.0, n_out=1.0):
    """2x2 complex transfer matrix M of the stack at scalar lambda_nm,
    (a_in, b_in)^T = M (a_out, b_out)^T (convention in module docstring).
    An empty layer list gives the bare n_in|n_out interface matrix."""
    m11, m12, m21, m22 = _m_elems(layers, float(lambda_nm), n_in, n_out)
    return np.array([[m11, m12], [m21, m22]], dtype=complex)


def rt_amplitudes(layers, lambda_nm, n_in=1.0, n_out=1.0):
    """Complex reflection and transmission amplitudes (r, t) for a wave
    incident from n_in; scalar or array lambda_nm. r is referenced to the
    front interface, t to the field just inside n_out at the back
    interface (no extra exit-medium propagation phase)."""
    m11, _, m21, _ = _m_elems(layers, lambda_nm, n_in, n_out)
    return m21 / m11, 1.0 / m11


def power_RT(layers, lambda_nm, n_in=1.0, n_out=1.0):
    """Power reflectance and transmittance (R, T): R = |r|^2,
    T = (n_out/n_in) |t|^2 (Poynting normalization). R + T = 1 for the
    lossless stacks this tier admits."""
    r, t = rt_amplitudes(layers, lambda_nm, n_in, n_out)
    return np.abs(r) ** 2, (float(n_out) / float(n_in)) * np.abs(t) ** 2


# ------------------------------------------------------------- DBR mirrors

def quarter_wave_stack(n_h, n_l, n_pairs, lambda0_nm, dn_dT_h=0.0, dn_dT_l=0.0):
    """(H L)^n_pairs quarter-wave stack, HIGH-index layer first (i.e.
    adjacent to the incidence medium), each layer of optical thickness
    lambda0/4: d = lambda0/(4 n). Returns a list of Layer."""
    d_h = lambda0_nm / (4.0 * n_h)
    d_l = lambda0_nm / (4.0 * n_l)
    out = []
    for _ in range(int(n_pairs)):
        out.append(Layer(n_h, d_h, dn_dT_h))
        out.append(Layer(n_l, d_l, dn_dT_l))
    return out


def qw_peak_reflectivity(n_in, n_out, n_h, n_l, n_pairs):
    """Textbook closed-form peak (Bragg-wavelength) power reflectivity of
    an (H L)^N quarter-wave stack, H adjacent to the incidence medium n_in
    and L adjacent to the exit medium n_out.

    Derivation (admittance transformation, Yariv & Yeh / Born & Wolf): a
    quarter-wave layer of index n maps the load admittance Y -> n^2/Y.
    Starting from Y = n_out and applying the 2N layers back-to-front
    (L then H per pair, N times) gives Y_in = (n_h/n_l)^(2N) * n_out, so

        r = (n_in - Y_in)/(n_in + Y_in),
        R = [ (1 - (n_out/n_in)(n_h/n_l)^(2N)) /
              (1 + (n_out/n_in)(n_h/n_l)^(2N)) ]^2.

    Exact at lambda0 for lossless quarter-wave layers; the general-N
    transfer matrix must reproduce it to machine precision."""
    b = (float(n_out) / float(n_in)) * (float(n_h) / float(n_l)) ** (2 * int(n_pairs))
    return ((1.0 - b) / (1.0 + b)) ** 2


def stopband_edges_analytic(n_h, n_l, lambda0_nm):
    """Analytic stopband edges (lambda_lo_nm, lambda_hi_nm) of an
    infinite quarter-wave Bragg mirror. In frequency the band is
    symmetric: f_pm/f0 = 1 +/- (2/pi) asin((n_h-n_l)/(n_h+n_l))
    (Yariv & Yeh, |trace| = 2 condition of the pair matrix), so
    lambda_pm = lambda0 / (f_pm/f0)."""
    x = (2.0 / np.pi) * np.arcsin((n_h - n_l) / (n_h + n_l))
    return lambda0_nm / (1.0 + x), lambda0_nm / (1.0 - x)


def dbr_reflectivity(n_h, n_l, n_pairs, lambda0_nm, n_in=1.0, n_out=1.0,
                     n_points=2001, span_frac=0.35):
    """R(lambda) spectrum of an (H L)^n_pairs quarter-wave mirror, plus
    stopband edges and peak R.

    Scans a frequency-uniform grid E = E0*(1 -/+ span_frac). The numeric
    stopband edges are taken as the contiguous half-peak region around
    lambda0 (first crossings of R_peak/2 walking outward) -- for
    many-pair mirrors R falls off a cliff at the band edge, so this
    converges onto the analytic (infinite-stack) edges, also returned.

    Returns dict: lambda_nm, R (arrays); R_peak, R_at_lambda0;
    stopband_lo_nm/stopband_hi_nm (numeric half-peak edges);
    stopband_analytic_nm (infinite-stack closed form, (lo, hi)).
    Tag: [E->analytic-1D, confirm: SIM-B]."""
    layers = quarter_wave_stack(n_h, n_l, n_pairs, lambda0_nm)
    E0 = HC_MEV_NM / lambda0_nm
    E = np.linspace(E0 * (1.0 - span_frac), E0 * (1.0 + span_frac), int(n_points))
    lam = HC_MEV_NM / E
    R, _ = power_RT(layers, lam, n_in, n_out)
    i0 = int(np.argmin(np.abs(lam - lambda0_nm)))
    R_peak = float(np.max(R))
    half = 0.5 * R_peak
    i_lo = i0
    while i_lo > 0 and R[i_lo - 1] >= half:
        i_lo -= 1
    i_hi = i0
    while i_hi < len(R) - 1 and R[i_hi + 1] >= half:
        i_hi += 1
    lam_edges = sorted((float(lam[i_lo]), float(lam[i_hi])))
    return {
        "lambda_nm": lam,
        "R": R,
        "R_peak": R_peak,
        "R_at_lambda0": float(R[i0]),
        "stopband_lo_nm": lam_edges[0],
        "stopband_hi_nm": lam_edges[1],
        "stopband_analytic_nm": stopband_edges_analytic(n_h, n_l, lambda0_nm),
    }


def penetration_depth(n_h, n_l, lambda0_nm, n_c=None, first_layer="H"):
    """Phase penetration depth L_pen (nm) of a semi-infinite quarter-wave
    Bragg mirror, seen from a medium of index n_c.

    CONVENTION (stated per project law): L_pen is the equivalent
    hard-mirror displacement defined by the reflection phase slope at the
    Bragg frequency,

        L_pen = (c / (2 n_c)) * |d phi_r / d omega|_{omega0}
              = (lambda0^2 / (4 pi n_c)) * |d phi_r / d lambda|_{lambda0},

    i.e. the mirror responds, to first order in detuning, like a
    fixed-phase mirror a distance L_pen into an index-n_c medium. This is
    the Babic & Corzine phase (tau) penetration class of definitions --
    NOT the energy-density penetration, which differs by O(1) factors.

    Closed forms (derived from the first-order admittance fixed point of
    the pair matrix; verified against the numeric phase slope in
    verify/verify_dbr.py):

      first_layer='H' (mirror starts with the high-index layer, the
      n_c <= n_l spacer geometry):
          L_pen = lambda0 / (4 (n_h - n_l))            [n_c-independent]
      first_layer='L' (mirror starts with the low-index layer, the
      high-index spacer geometry, e.g. n_c = n_h):
          L_pen = lambda0 / (4 (n_h - n_l)) * n_h n_l / n_c^2

    Both are exact to first order in the detuning for the semi-infinite
    stack. Tag: [E->analytic-1D, confirm: SIM-B]."""
    dn = float(n_h) - float(n_l)
    if dn <= 0:
        raise ValueError("need n_h > n_l")
    base = lambda0_nm / (4.0 * dn)
    if first_layer.upper() == "H":
        return base
    if first_layer.upper() == "L":
        n_c = float(n_h) if n_c is None else float(n_c)
        return base * (float(n_h) * float(n_l)) / n_c**2
    raise ValueError("first_layer must be 'H' or 'L'")


# ------------------------------------------------------------- lambda cavity

def cavity_stack(n_h, n_l, n_c, n_pairs_top, n_pairs_bottom, lambda0_nm,
                 dn_dT_h=0.0, dn_dT_l=0.0, dn_dT_c=0.0):
    """Planar lambda-cavity: top DBR + spacer + bottom DBR.

    Structure, listed from the incidence side:
        n_in | (H L)^n_pairs_top | spacer (n_c, d = lambda0/n_c) |
             (L H)^n_pairs_bottom | n_out
    The spacer physical thickness lambda0/n_c is one optical wavelength
    (a "lambda cavity"); both mirrors present their LOW-index layer to
    the spacer, so for n_c > n_l their on-resonance reflection phase is 0
    and the mode sits at lambda0 by construction (round-trip phase
    2 k n_c d + 0 + 0 = 4 pi).

    Returns (top_layers, spacer_layer, bottom_layers); the flat stack for
    transfer-matrix calls is top + [spacer] + bottom."""
    top = quarter_wave_stack(n_h, n_l, n_pairs_top, lambda0_nm, dn_dT_h, dn_dT_l)
    spacer = Layer(n_c, lambda0_nm / n_c, dn_dT_c)
    bottom = []
    for _ in range(int(n_pairs_bottom)):
        bottom.append(Layer(n_l, lambda0_nm / (4.0 * n_l), dn_dT_l))
        bottom.append(Layer(n_h, lambda0_nm / (4.0 * n_h), dn_dT_h))
    return top, spacer, bottom


def cavity_mode(layers, lambda0_nm, n_in=1.0, n_out=1.0, span_frac=0.05):
    """Locate the cavity resonance in T(E) near the design energy and
    measure its linewidth.

    Method (adaptive, second-method-friendly): coarse frequency-uniform
    scan of the transmission over E0*(1 +/- span_frac); iteratively zoom
    the grid onto the argmax until the spacing resolves the peak to
    ~1e-11*E0; then find the two half-maximum crossings by bracketed
    bisection (march outward with doubling steps to bracket, 80
    bisections each side). Works from broad few-pair modes to Q ~ 1e6
    without tuning.

    Returns dict:
      E_cav_meV : resonance energy (transmission peak)
      kappa_meV : full width at half maximum of T(E)
      Q         : E_cav/kappa
      T_peak    : peak power transmission (=1 for symmetric lossless)
    kappa_meV is NaN if a half-max crossing is not found within +/-45% of
    E_cav (no resolvable resonance). Tag: [E->analytic-1D, confirm: SIM-B]."""
    E0 = HC_MEV_NM / lambda0_nm

    def T_of(E):
        return power_RT(layers, HC_MEV_NM / np.asarray(E, dtype=float),
                        n_in, n_out)[1]

    lo, hi = E0 * (1.0 - span_frac), E0 * (1.0 + span_frac)
    npts = 801
    Eg = np.linspace(lo, hi, npts)
    for _ in range(64):
        T = T_of(Eg)
        i = int(np.argmax(T))
        dE = Eg[1] - Eg[0]
        if dE < 1e-11 * E0:
            break
        lo, hi = Eg[i] - 2.0 * dE, Eg[i] + 2.0 * dE
        Eg = np.linspace(lo, hi, npts)
    # sub-grid parabolic refinement of the peak position
    E_p, T_p = float(Eg[i]), float(T[i])
    if 0 < i < npts - 1:
        y0, y1, y2 = float(T[i - 1]), float(T[i]), float(T[i + 1])
        denom = y0 - 2.0 * y1 + y2
        if denom < 0.0:
            E_p = float(Eg[i]) + 0.5 * dE * (y0 - y2) / denom
            T_p = float(T_of(E_p))

    def half_cross(sign):
        """Distance from E_p to the half-max crossing on one side."""
        step = max(dE, 1e-9 * E0)
        delta = step
        while float(T_of(E_p + sign * delta)) > 0.5 * T_p:
            delta *= 2.0
            if delta > 0.45 * E_p:
                return np.nan
        a, b = (delta / 2.0 if delta > step else 0.0), delta
        for _ in range(80):
            mid = 0.5 * (a + b)
            if float(T_of(E_p + sign * mid)) > 0.5 * T_p:
                a = mid
            else:
                b = mid
        return 0.5 * (a + b)

    kappa = half_cross(+1.0) + half_cross(-1.0)
    return {
        "E_cav_meV": E_p,
        "kappa_meV": float(kappa),
        "Q": float(E_p / kappa),
        "T_peak": T_p,
    }


def dEdT_cav(layers, lambda0_nm, n_in=1.0, n_out=1.0, T0=300.0, dT=10.0,
             T_ref=300.0, span_frac=0.05):
    """dE_cav/dT in meV/K by central finite difference: re-solve
    cavity_mode with every layer's index shifted by its own thermo-optic
    slope, n(T) = n + dn_dT*(T - T_ref), at T0 +/- dT.

    Physics/limits: linear thermo-optic model only; thermal expansion and
    dn_dT dispersion neglected. Sanity anchor: if all indices scaled by a
    common relative amount, dE/E = -dn/n exactly (optical paths scale);
    with per-layer additive dn/dT the result is the mode-weighted
    average of that estimate. Typical III-V dn/dT > 0 gives dE/dT < 0
    (red shift), the slow mode drift the F6(iv) tracking rule relies on.
    Tag: [E->analytic-1D, confirm: SIM-B]."""
    def E_at(T):
        mod = [L.at_T(T, T_ref) for L in layers]
        return cavity_mode(mod, lambda0_nm, n_in, n_out, span_frac)["E_cav_meV"]

    return (E_at(T0 + dT) - E_at(T0 - dT)) / (2.0 * dT)


# --------------------------------------------------------- planar Purcell

def planar_purcell(top_layers, bottom_layers, n_c, d_c_nm, z_dot_nm,
                   lambda_nm, n_in=1.0, n_out=1.0):
    """1-D (planar) LDOS/spontaneous-emission enhancement for a dipole in
    the spacer between two mirror stacks.

    F(z) = Re[ (1 + r1 e^{2i phi1}) (1 + r2 e^{2i phi2})
               / (1 - r1 r2 e^{2i(phi1+phi2)}) ]

    where r1 (r2) is the complex reflection amplitude of the top (bottom)
    stack as seen FROM the dot plane inside the spacer -- computed by
    transfer matrix of each half-stack with the spacer index as incidence
    medium -- and phi1 = k n_c z_dot, phi2 = k n_c (d_c - z_dot) are the
    one-way propagation phases from the dot to each mirror
    (k = 2 pi/lambda). This is the standard two-mirror multiple-round-trip
    (image-sum) result for a 1-D emitter, e.g. Bjork et al. / Benisty et
    al. planar-microcavity treatments, normalized so F = 1 in free space.

    Geometry: z_dot_nm measured from the top-mirror/spacer interface into
    the spacer (0 <= z_dot_nm <= d_c_nm). top_layers are ordered from the
    outside (n_in side) toward the spacer, exactly as returned by
    cavity_stack; they are reversed internally for the looking-up
    reflection. bottom_layers are ordered spacer -> n_out.

    HONESTY WALL: this is the 1-D/planar estimate of the emission-rate
    modification -- it integrates no transverse mode profile and is NOT a
    3-D energy-normalized V_eff Purcell factor. Use it for trend
    pre-screening only. Tag: [E->analytic-1D, confirm: SIM-B].
    Limits: r1 = r2 = 0 gives F = 1 exactly; a symmetric high-R cavity on
    resonance gives F ~ (1+r)/(1-r) ~ 4/(1-R) at an antinode and F < 1 at
    a node (both exercised in verify/verify_dbr.py)."""
    if not 0.0 <= z_dot_nm <= d_c_nm:
        raise ValueError("dot must sit inside the spacer")
    r1, _ = rt_amplitudes(list(reversed(top_layers)), lambda_nm,
                          n_in=n_c, n_out=n_in)
    r2, _ = rt_amplitudes(bottom_layers, lambda_nm, n_in=n_c, n_out=n_out)
    k = 2.0 * np.pi / lambda_nm
    phi1 = k * n_c * z_dot_nm
    phi2 = k * n_c * (d_c_nm - z_dot_nm)
    e1, e2 = np.exp(2j * phi1), np.exp(2j * phi2)
    F = (1.0 + r1 * e1) * (1.0 + r2 * e2) / (1.0 - r1 * r2 * e1 * e2)
    return float(np.real(F))


# ---------------------------------------------------------------- inversion

def invert_kappa(kappa_target_meV, n_h, n_l, n_c, lambda0_nm,
                 n_top_max=30, n_bottom_extra=6, n_in=1.0, n_out=None):
    """INVERSE direction (the point of this tier): smallest top-mirror
    pair count whose lambda-cavity linewidth meets kappa <= target.

    Brute-force sweep n_top = 1..n_top_max with n_bottom = n_top +
    n_bottom_extra (bottom mirror kept stronger so the top port
    dominates the loss, the extraction geometry), each evaluated by a
    direct cavity_mode solve -- milliseconds per configuration, which is
    the point of the analytic tier. Substrate defaults to n_out = n_h.

    Returns dict {n_top, n_bottom, kappa_meV, Q, E_cav_meV, L_pen_nm}
    (L_pen for the mirror seen from the spacer, low-index layer first).
    Raises ValueError if the target is not met by n_top_max.
    Tag: [E->analytic-1D, confirm: SIM-B]."""
    n_out = float(n_h) if n_out is None else float(n_out)
    for n_top in range(1, int(n_top_max) + 1):
        n_bottom = n_top + int(n_bottom_extra)
        top, spacer, bottom = cavity_stack(n_h, n_l, n_c, n_top, n_bottom,
                                           lambda0_nm)
        mode = cavity_mode(top + [spacer] + bottom, lambda0_nm,
                           n_in=n_in, n_out=n_out)
        if np.isfinite(mode["kappa_meV"]) and mode["kappa_meV"] <= kappa_target_meV:
            return {
                "n_top": n_top,
                "n_bottom": n_bottom,
                "kappa_meV": mode["kappa_meV"],
                "Q": mode["Q"],
                "E_cav_meV": mode["E_cav_meV"],
                "L_pen_nm": penetration_depth(n_h, n_l, lambda0_nm,
                                              n_c=n_c, first_layer="L"),
            }
    raise ValueError(
        f"kappa <= {kappa_target_meV} meV not reachable with n_top <= {n_top_max} "
        f"at contrast {n_h}/{n_l} (raise n_top_max or the contrast)"
    )
