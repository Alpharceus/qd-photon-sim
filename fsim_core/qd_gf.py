"""Module C Tier-3 (Phase T2 core) -- independent-boson-model (IBM) polaron
lineshape with a dot-geometry form factor.

Physics:
  * IBM: a two-level exciton coupled diagonally to bulk LA phonons through the
    deformation potential is exactly solvable (polaron / cumulant Green's
    function). The emission dipole correlator, with the polaron shift absorbed
    into the ZPL position, is

        C(t) = exp(phi(t) - phi(0)),
        phi(t) = int_0^inf dE  J(E)/E^2
                 * [coth(E / (2 kB T)) cos(E t / hbar) - i sin(E t / hbar)],

    and the unit-area emission spectrum (photon energy measured RELATIVE TO
    THE ZPL, blue positive) is

        S(omega) = (1 / (pi hbar)) Re int_0^inf dt
                   e^{-i omega t / hbar} C(t) e^{-Gamma_zpl t / (2 hbar)}.

  * SIGN / SIDE CONVENTION (derived, and checked in verify/verify_gf.py):
    with the -i sin cumulant above and the e^{-i omega t/hbar} transform,
    expanding C(t) to one phonon gives (n+1) e^{-iEt/hbar} + n e^{+iEt/hbar},
    i.e. weight (n+1) lands at omega = -E and weight n at omega = +E. So on
    this axis omega < 0 (photon RED of the ZPL) is the phonon-EMISSION
    (Stokes) side -- the photon leaves with less energy because an LA phonon
    carries the difference -- and it dominates as T -> 0; omega > 0 (blue) is
    phonon absorption (anti-Stokes) with detailed balance
    S(+E)/S(-E) = exp(-E / kB T).

  * Spectral density (super-ohmic, bulk LA + Gaussian confinement):
        J(E) = alpha_meV * E^3 * F2(E)      [J in meV when E is in meV]
    with the anisotropic-Gaussian exciton form factor, angular average done
    numerically (Gauss-Legendre over u = cos(theta)):
        F2(E) = <exp(-(E/hbar)^2 (l_xy^2 sin^2 th + l_z^2 cos^2 th)
                     / (2 c_s^2))>_solid angle.
    CONVENTION: the /2 in the exponent, so the spherical limit l_xy = l_z = l
    collapses to F2 = exp(-(E/E_c)^2) with E_c = hbar * sqrt(2) * c_s / l --
    the standard omega_c = sqrt(2) c_s / l cutoff. (Some papers omit the /2
    and quote omega_c = c_s/l; we do not.)

  * Coupling constant: in SI, alpha_SI = (D_e - D_h)^2 / (4 pi^2 rho_m hbar
    c_s^5) with units of s^2 (J(omega) = alpha omega^3 for omega in rad/s).
    Literature class for InAs/GaAs dots: alpha ~ 0.01-0.1 ps^2 (0.027 ps^2 in
    Ramsay et al.). Internally we use the energy-axis constant
    alpha_meV = alpha_ps2 / hbar^2  [meV^-2], so that J(E) = alpha_meV E^3 F2
    is in meV and int dE J/E^2 is dimensionless.
    Both entry paths are exposed: material constants {D_e_eV, D_h_eV,
    rho_m_kg_m3, c_s_m_s} or a direct alpha_ps2 override.

  * ZPL weight (Franck-Condon / Debye-Waller): phi(t) -> 0 as t -> inf (the
    oscillatory integrand dephases; there is no residual linear phase because
    J/E -> 0 at E -> 0), so C(t) -> exp(-phi(0)) = Z with the total
    Huang-Rhys factor
        S_total(T) = phi(0) = int_0^inf dE J(E)/E^2 coth(E / (2 kB T)),
        Z(T) = exp(-S_total(T)),
    and at T = 0, S_total -> S_0 = int dE J/E^2 (coth -> 1).

  * Gamma_zpl is a PHENOMENOLOGICAL Lorentzian ZPL width (radiative + pure
    dephasing) that the strict IBM does not produce (the IBM ZPL is a delta
    function); it multiplies C(t) as exp(-Gamma_zpl t / 2 hbar), i.e. the ZPL
    becomes a Lorentzian of FWHM Gamma_zpl sitting on the phonon-sideband
    pedestal. See effective_gamma_zpl for the bridge to the validated
    Gamma(T) model -- and its honest double-counting caveat.

Numerics: the spectrum is computed exactly as ZPL + sideband,
    S = Z * lorentzian(omega; Gamma_zpl) + FFT[(C(t) - Z) * envelope],
which is algebraically identical to transforming C(t) whole but needs only a
~80 ps time grid (the sideband correlator C - Z dies on the phonon timescale,
a few ps; thermal decay rate ~ 2 pi kB T / hbar limits validity to T >~ 0.5 K
for the fixed grid used here). Time grid: dt = 5 fs, 16384 points,
zero-padded to 131072 (frequency span ~ +-413 meV, resolution ~ 6 ueV).

Units: energies/omega in meV, T in K, lengths in nm, time in ps.
hbar = 0.6582119569 meV ps; kB = 0.08617333262 meV/K (from fsim_core.spectral).
Headless: numpy only (three-layer rule).
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np

from fsim_core.spectral import (
    KB,
    cavity_transmission,
    combined_transmission,
    gamma_of_T,
    lorentzian,
    tophat_transmission,
)

HBAR = 0.6582119569       # meV*ps
_EV_J = 1.602176634e-19   # J per eV (exact, SI 2019)
_HBAR_SI = 1.054571817e-34  # J*s

# time/frequency grid for the polaron correlator (see module docstring)
_DT_PS = 0.005
_NT = 16384
_NPAD = 131072

# Gauss-Legendre nodes/weights for the angular average over u = cos(theta),
# mapped to [0, 1] (the integrand is even in u); weights sum to 1.
_GL_X, _GL_W = np.polynomial.legendre.leggauss(48)
_U_NODES = 0.5 * (_GL_X + 1.0)
_U_WEIGHTS = 0.5 * _GL_W


# ------------------------------------------------------------------ parameters

@dataclass(frozen=True)
class PhononParams:
    """Exciton--LA-phonon coupling card (IBM inputs).

    Material constants ([DR] tags: literature-class values, honest status):
      D_e_eV, D_h_eV : electron/hole deformation potentials [eV]. Defaults are
          the InP conduction/valence hydrostatic potentials a_c = -6.0 eV,
          a_v = -0.6 eV (Vurgaftman 2001 compilation class) used as a
          CLASS-PROXY for the InP/GaAsP visible-platform dots [DR]; only the
          difference D_e - D_h enters. The same-form-factor approximation for
          electron and hole (single exciton form factor) is standard IBM
          practice [E].
      rho_m_kg_m3 : mass density [kg/m^3]; 4810 = bulk InP [DR].
      c_s_m_s : LA sound velocity [m/s]; 4600 = InP-class longitudinal [DR].

    Geometry (Gaussian confinement lengths of the exciton wavefunction,
    psi ~ exp(-r_xy^2/(2 l_xy^2) - z^2/(2 l_z^2)) -- NOT the physical dot
    dimensions; the ground state is smaller than the island):
      l_xy_nm, l_z_nm : defaults 4.5 / 1.5 nm, a [DR] class-proxy mapping of
          the tier-plan default lens dot (base ~25 nm, height ~3 nm:
          l_z ~ height/2, l_xy ~ base/5.5). Replace when AFM/TEM (test #8)
          and a confinement calculation exist; [A, null] until then.

    alpha_ps2 : direct coupling override [ps^2]. None (default) -> compute
          from the material constants. Set it to pin alpha to a measured or
          literature value (e.g. 0.027 ps^2, InGaAs/GaAs, Ramsay-class).
    """
    D_e_eV: float = -6.0        # [DR] InP a_c, class-proxy
    D_h_eV: float = -0.6        # [DR] InP a_v, class-proxy
    rho_m_kg_m3: float = 4810.0  # [DR] InP density
    c_s_m_s: float = 4600.0     # [DR] InP LA sound velocity
    l_xy_nm: float = 4.5        # [DR] class-proxy, 25 nm-base lens dot
    l_z_nm: float = 1.5         # [DR] class-proxy, 3 nm-height lens dot
    alpha_ps2: float | None = None  # direct override; None -> from materials


def coupling_alpha_ps2(params: PhononParams) -> float:
    """IBM coupling alpha in ps^2 (J(omega) = alpha omega^3 exp(...) with
    omega in rad/ps).

    alpha_SI = (D_e - D_h)^2 / (4 pi^2 rho_m hbar c_s^5)  [s^2], converted to
    ps^2 (1 ps^2 = 1e-24 s^2). If params.alpha_ps2 is set, that override is
    returned unchanged. Plausibility anchor: InAs/GaAs literature class is
    ~0.01-0.1 ps^2; the default InP-class card lands at ~0.018 ps^2.
    """
    if params.alpha_ps2 is not None:
        return float(params.alpha_ps2)
    d_j = (params.D_e_eV - params.D_h_eV) * _EV_J
    alpha_si = d_j**2 / (
        4.0 * np.pi**2 * params.rho_m_kg_m3 * _HBAR_SI * params.c_s_m_s**5
    )
    return alpha_si * 1e24


def coupling_alpha_meV(params: PhononParams) -> float:
    """alpha on the energy axis [meV^-2]: alpha_meV = alpha_ps2 / hbar^2, so
    that J(E) = alpha_meV E^3 F2(E) is in meV for E in meV and int dE J/E^2
    (the Huang-Rhys integral) is dimensionless."""
    return coupling_alpha_ps2(params) / HBAR**2


# -------------------------------------------------------------- spectral density

def _form_factor_sq(E_meV: np.ndarray, params: PhononParams) -> np.ndarray:
    """Angular-averaged squared exciton form factor F2(E), dimensionless.

    F2(E) = int_0^1 du exp(-(E/hbar)^2 [l_xy^2 (1-u^2) + l_z^2 u^2]/(2 c_s^2)),
    u = cos(theta); Gauss-Legendre (48 nodes, effectively exact for these
    smooth integrands). Spherical limit l_xy = l_z = l: F2 = exp(-(E/E_c)^2),
    E_c = hbar sqrt(2) c_s / l  (the /2-exponent convention; see module
    docstring). F2(0) = 1; F2 -> 0 super-exponentially at large E.
    """
    c_nmps = params.c_s_m_s * 1e-3  # nm/ps
    l_eff2 = (params.l_xy_nm**2 * (1.0 - _U_NODES**2)
              + params.l_z_nm**2 * _U_NODES**2)          # (48,)
    x = np.multiply.outer((E_meV / HBAR) ** 2, l_eff2) / (2.0 * c_nmps**2)
    return np.exp(-x) @ _U_WEIGHTS


def spectral_density(omega_meV, params: PhononParams):
    """IBM spectral density J(E) [meV] on the energy axis E = hbar*omega [meV].

    J(E) = alpha_meV * E^3 * F2(E): super-ohmic (E^3, bulk LA deformation
    potential -- no piezoelectric or LO term, an [E] scope limit) with the
    anisotropic-Gaussian geometry cutoff. J(E <= 0) = 0; J -> 0 at both ends,
    peaking near the geometric cutoff E_c. Scalar in -> scalar out.
    """
    E = np.atleast_1d(np.asarray(omega_meV, dtype=float))
    J = np.zeros_like(E)
    pos = E > 0.0
    if np.any(pos):
        J[pos] = coupling_alpha_meV(params) * E[pos] ** 3 * _form_factor_sq(E[pos], params)
    if np.ndim(omega_meV) == 0:
        return float(J[0])
    return J


# ------------------------------------------------------------ cumulant / weights

def _thermal_weights(params: PhononParams, T: float):
    """Energy grid E plus the two cumulant weights:
    g(E) = J/E^2 (the -i sin weight) and p(E) = g(E) coth(E/2kBT) (the cos
    weight), with the E -> 0 limits taken analytically (g -> 0;
    p -> 2 alpha_meV kB T since coth -> 2kBT/E and F2(0) = 1)."""
    if T < 0.0:
        raise ValueError("T must be >= 0")
    c_nmps = params.c_s_m_s * 1e-3
    e_c_max = HBAR * np.sqrt(2.0) * c_nmps / min(params.l_xy_nm, params.l_z_nm)
    e_hi = 12.0 * e_c_max
    n_e = max(2001, int(round(e_hi / 0.01)) + 1)
    E = np.linspace(0.0, e_hi, n_e)
    J = spectral_density(E, params)
    g = np.zeros_like(E)
    g[1:] = J[1:] / E[1:] ** 2
    if T > 0.0:
        p = np.empty_like(E)
        p[1:] = g[1:] / np.tanh(E[1:] / (2.0 * KB * T))
        p[0] = 2.0 * coupling_alpha_meV(params) * KB * T
    else:
        p = g.copy()
    return E, p, g


def _phi(t_ps: np.ndarray, params: PhononParams, T: float) -> np.ndarray:
    """Cumulant phi(t) = int dE [p(E) cos(Et/hbar) - i g(E) sin(Et/hbar)],
    trapezoid over the _thermal_weights grid, chunked over t for memory.
    The E-grid spacing (~0.01 meV) aliases the time axis with period
    2 pi hbar / dE ~ 400 ps >> the 82 ps window: alias error is negligible
    because phi has decayed by then (validity T >~ 0.5 K, see module doc)."""
    E, p, g = _thermal_weights(params, T)
    out = np.empty(t_ps.size, dtype=complex)
    chunk = 256
    for i in range(0, t_ps.size, chunk):
        ph = np.multiply.outer(t_ps[i:i + chunk], E) / HBAR
        re = np.trapezoid(p * np.cos(ph), E, axis=1)
        im = np.trapezoid(g * np.sin(ph), E, axis=1)
        out[i:i + chunk] = re - 1j * im
    return out


def huang_rhys(params: PhononParams, T: float) -> float:
    """Total Huang-Rhys factor S_total(T) = int_0^inf dE J(E)/E^2
    coth(E/(2 kB T)) = phi(0). Monotone increasing in T (coth is); T = 0
    limit is S_0 = int dE J/E^2 (spherical closed form: alpha_meV E_c^2 / 2).
    High-T asymptote: S ~ 2 alpha_meV kB T int F2 dE (linear in T)."""
    E, p, _ = _thermal_weights(params, float(T))
    return float(np.trapezoid(p, E))


def zpl_weight(params: PhononParams, T: float) -> float:
    """Franck-Condon / Debye-Waller ZPL weight Z(T) = exp(-S_total(T)),
    in (0, 1]: the fraction of the emission that stays in the zero-phonon
    line. Derivation: phi(t->inf) -> 0 (Riemann-Lebesgue on the smooth
    integrand; no residual phase since J/E -> 0), so C(inf) = exp(-phi(0))."""
    return float(np.exp(-huang_rhys(params, T)))


def sideband_fraction(params: PhononParams, T: float) -> float:
    """1 - Z(T): fraction of the line emitted into the phonon sidebands.
    This is the weight a ZPL-tuned narrow filter can never collect."""
    return 1.0 - zpl_weight(params, T)


# ------------------------------------------------------------------- spectrum

@lru_cache(maxsize=32)
def _phi_window(params: PhononParams, T: float):
    """phi(t) on the module's fixed FFT time window, cached by (params, T).
    Split out of _spectrum_parts so the X and XX lines at the same junction
    temperature (same params, different ZPL widths) share the expensive
    cumulant integral -- the FFT below is cheap by comparison."""
    t = np.arange(_NT) * _DT_PS
    return t, _phi(t, params, T)


@lru_cache(maxsize=16)
def _spectrum_parts(params: PhononParams, T: float, gamma_zpl_meV: float):
    """(omega_axis, S_sideband(omega), Z): the FFT sideband spectrum on its
    native axis (fftshifted, ascending, ~+-413 meV span, ~6 ueV spacing) and
    the ZPL weight. Exact decomposition C = Z + (C - Z): the ZPL term
    transforms analytically to Z * lorentzian(omega; Gamma_zpl); only the
    fast-decaying sideband correlator is transformed numerically (trapezoid
    endpoint handled via the -f(0)/2 correction). Cached; treat the returned
    arrays as read-only."""
    t, phi = _phi_window(params, T)
    s_total = phi[0].real
    z = float(np.exp(-s_total))
    c_sb = np.exp(phi - phi[0]) - z
    f = c_sb * np.exp(-gamma_zpl_meV * t / (2.0 * HBAR))
    spec = (_DT_PS / (np.pi * HBAR)) * (np.fft.fft(f, n=_NPAD) - 0.5 * f[0]).real
    omega = np.fft.fftshift(2.0 * np.pi * HBAR * np.fft.fftfreq(_NPAD, d=_DT_PS))
    return omega, np.fft.fftshift(spec), z


def ibm_spectrum(omega_grid_meV, params: PhononParams, T: float,
                 gamma_zpl_meV: float) -> np.ndarray:
    """Unit-area IBM emission spectrum on the given grid [1/meV].

    omega_grid_meV: photon energy relative to the (polaron-shifted) ZPL,
    blue positive; the phonon-emission (Stokes) sideband therefore sits at
    omega < 0 and dominates at low T (see module docstring for the
    derivation). S = Z * Lorentzian(FWHM Gamma_zpl) + sidebands, normalized
    numerically on the grid (which must be dense vs Gamma_zpl and wide enough
    to cover the line: asserted at the 10% level; tiny FFT ringing < 1e-10 of
    the peak is clipped to zero and asserted).
    """
    grid = np.asarray(omega_grid_meV, dtype=float)
    omega, s_sb, z = _spectrum_parts(params, float(T), float(gamma_zpl_meV))
    s = z * lorentzian(grid, 0.0, gamma_zpl_meV) + np.interp(grid, omega, s_sb,
                                                             left=0.0, right=0.0)
    peak = float(s.max())
    assert float(s.min()) > -1e-10 * peak, "FFT ringing above tolerance"
    s = np.clip(s, 0.0, None)
    area = float(np.trapezoid(s, grid))
    assert abs(area - 1.0) < 0.1, f"grid does not cover the line (raw area {area:.4f})"
    return s / area


def ibm_transmission(delta_meV, params: PhononParams, T: float,
                     gamma_zpl_meV: float, w_meV=None, kappa_meV=None):
    """Fraction of the IBM line transmitted through the filter stack --
    the sideband-aware replacement for the Lorentzian t_X/t_XX in the F1
    eps path.

    The line center sits at delta_meV from the filter center (same delta
    convention as spectral.tophat_transmission / cavity_transmission).
    Acceptance mirrors fsim_core.spectral EXACTLY: top-hat of FULL width
    w_meV (|x| <= w/2), peak-normalized Lorentzian cavity of FWHM kappa_meV,
    or their product; both None -> 1.0 (unit-area line fully collected).

    Computed as [Z * t_closed + int S_sb * acceptance] / [Z + int S_sb]:
    the ZPL Lorentzian goes through spectral.py's own closed forms (so the
    coupling -> 0 limit reproduces them identically), the sideband through a
    numeric trapezoid on the FFT axis. Physics: sidebands put weight far from
    the ZPL, so narrow filters lose it (brightness) while wide/offset filters
    admit the other line's sideband (purity) -- the 'Lorentzian optimism'
    correction. Scalar or array delta_meV.
    """
    scalar = np.ndim(delta_meV) == 0
    deltas = np.atleast_1d(np.asarray(delta_meV, dtype=float))
    if w_meV is None and kappa_meV is None:
        out = np.ones_like(deltas)
        return float(out[0]) if scalar else out
    omega, s_sb, z = _spectrum_parts(params, float(T), float(gamma_zpl_meV))
    sb_area = float(np.trapezoid(s_sb, omega))
    out = np.empty_like(deltas)
    for i, d in enumerate(deltas):
        if w_meV is not None and kappa_meV is not None:
            t_zpl = combined_transmission(d, gamma_zpl_meV, w_meV, kappa_meV)
        elif w_meV is not None:
            t_zpl = tophat_transmission(d, gamma_zpl_meV, w_meV)
        else:
            t_zpl = cavity_transmission(d, gamma_zpl_meV, kappa_meV)
        x = omega + d  # photon energy relative to the filter center
        acc = np.ones_like(x)
        if w_meV is not None:
            acc *= np.abs(x) <= w_meV / 2.0
        if kappa_meV is not None:
            acc *= (kappa_meV / 2.0) ** 2 / (x**2 + (kappa_meV / 2.0) ** 2)
        t_sb = float(np.trapezoid(s_sb * acc, omega))
        out[i] = (z * t_zpl + t_sb) / (z + sb_area)
    return float(out[0]) if scalar else out


# ------------------------------------------------------------ broadening bridge

def effective_gamma_zpl(T, gamma0, a_ac, b_lo, E_lo):
    """ZPL width to feed ibm_spectrum/ibm_transmission: the standing,
    V-a-validated Gamma(T) = gamma0 + a_ac*T + b_lo*n_BE(E_lo,T) form,
    imported unchanged from fsim_core.spectral.gamma_of_T (core->core).

    HONESTY NOTE [E]: the a_ac*T linear term was fitted on data whose acoustic
    dephasing a full phonon treatment (IBM + beyond-IBM real dephasing) would
    partly generate microscopically, so pairing it with IBM sidebands
    double-counts some acoustic broadening. Kept deliberately: the IBM alone
    produces NO finite ZPL width (its ZPL is a delta), so the fitted form is
    still the best carrier of the observed width. The V-a refit gate (T2
    integrator task, not this module) decides whether the double-count
    matters at the decision points.
    """
    return gamma_of_T(T, gamma0, a_ac, b_lo, E_lo)
