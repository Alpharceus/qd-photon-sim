"""Module C -- Spectral closed forms (F1/F1a) and their MC second methods.

Physics (F-series):
  * X and XX emission lines are unit-area Lorentzians; the filter is centered on X.
  * F1 (filtered-cascade identity):  g2_0 = eps = t_XX / t_X, where t_i is the
    transmission of line i through the acceptance function.
  * Acceptance functions: top-hat window of full width w, Lorentzian cavity of
    FWHM kappa (peak-normalized), or their product (numeric).

Units: energies in meV, temperatures in K, throughout.
Every closed form here has an MC second method (mc_*) exercised by verify/.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.integrate import quad

KB = 0.08617333262  # Boltzmann constant, meV/K


# ----------------------------------------------------------------------------- lineshape

def lorentzian(x, center, fwhm):
    """Unit-area Lorentzian."""
    hg = fwhm / 2.0
    return (hg / np.pi) / ((x - center) ** 2 + hg**2)


def gamma_of_T(T, gamma0, a_ac, b_lo, E_lo):
    """Linewidth model: Gamma(T) = Gamma0 + a_ac*T + b_lo * n_BE(E_lo, T).

    Gamma0: zero-T width; a_ac: acoustic-phonon linear term; b_lo: LO-phonon
    coupling with Bose-Einstein occupation of the E_lo phonon.
    """
    T = np.asarray(T, dtype=float)
    x = np.minimum(E_lo / (KB * T), 700.0)
    return gamma0 + a_ac * T + b_lo / np.expm1(x)


# ------------------------------------------------------------------------- transmissions

def tophat_transmission(delta, gamma, w):
    """Fraction of a Lorentzian (FWHM gamma, centered at delta from the filter
    center) transmitted by a top-hat window of full width w."""
    hg = gamma / 2.0
    return (np.arctan((w / 2.0 - delta) / hg) + np.arctan((w / 2.0 + delta) / hg)) / np.pi


def cavity_transmission(delta, gamma, kappa):
    """Lorentzian line through a peak-normalized Lorentzian cavity filter (FWHM
    kappa). Closed form via Lorentzian-Lorentzian convolution: the result is a
    Lorentzian of FWHM gamma+kappa. At delta=0 this reduces to kappa/(gamma+kappa)."""
    s = (gamma + kappa) / 2.0
    return (kappa * (gamma + kappa) / 4.0) / (delta**2 + s**2)


def combined_transmission(delta, gamma, w, kappa):
    """Top-hat window times cavity Lorentzian (numeric; no closed form kept)."""
    val, _ = quad(
        lambda x: lorentzian(x, delta, gamma) * (kappa / 2.0) ** 2 / (x**2 + (kappa / 2.0) ** 2),
        -w / 2.0, w / 2.0, limit=200,
    )
    return val


def transmission(delta, gamma, w=None, kappa=None):
    if w is not None and kappa is not None:
        return combined_transmission(delta, gamma, w, kappa)
    if w is not None:
        return tophat_transmission(delta, gamma, w)
    if kappa is not None:
        return cavity_transmission(delta, gamma, kappa)
    return 1.0  # no filter: unit-area line fully collected


# ------------------------------------------------------------------------------- epsilon

@dataclass
class SpectralResult:
    eps: float   # leakage ratio t_XX/t_X = g2_0 (F1)
    t_x: float   # X brightness factor through the filter
    t_xx: float


def epsilon(delta_xx, gamma_x, gamma_xx, w=None, kappa=None) -> SpectralResult:
    """F1: eps = t_XX/t_X for an XX line detuned by delta_xx from the X-centered filter."""
    t_x = transmission(0.0, gamma_x, w, kappa)
    t_xx = transmission(delta_xx, gamma_xx, w, kappa)
    return SpectralResult(eps=t_xx / t_x, t_x=t_x, t_xx=t_xx)


def epsilon_narrow_filter(delta_xx, gamma_x, gamma_xx):
    """F1a narrow-filter bound (w -> 0): ratio of on-center spectral densities."""
    return (gamma_x * gamma_xx / 4.0) / (delta_xx**2 + (gamma_xx / 2.0) ** 2)


def optimal_width(delta_xx, gamma_x, gamma_xx, eps_max, w_hi=None):
    """Widest top-hat w with eps(w) <= eps_max -- maximizes brightness t_X at the
    allowed leakage. Returns np.nan if even w->0 violates the bound."""
    from scipy.optimize import brentq

    if epsilon_narrow_filter(delta_xx, gamma_x, gamma_xx) > eps_max:
        return np.nan
    w_hi = w_hi if w_hi is not None else 20.0 * (abs(delta_xx) + gamma_x + gamma_xx)
    f = lambda w: epsilon(delta_xx, gamma_x, gamma_xx, w=w).eps - eps_max
    if f(w_hi) < 0:
        return w_hi  # bound never binds in range
    return brentq(f, 1e-9, w_hi)


# ---------------------------------------------------------------- MC second methods (Yuki)

def mc_transmission(delta, gamma, w=None, kappa=None, n=1_000_000, seed=0):
    """Monte-Carlo second method for every transmission closed form: draw photon
    frequencies from the Lorentzian (Cauchy) line, apply the acceptance."""
    rng = np.random.default_rng(seed)
    x = delta + (gamma / 2.0) * rng.standard_cauchy(n)
    acc = np.ones(n, dtype=bool)
    if w is not None:
        acc &= np.abs(x) <= w / 2.0
    if kappa is not None:
        p = (kappa / 2.0) ** 2 / (x**2 + (kappa / 2.0) ** 2)
        acc &= rng.random(n) < p
    t = acc.mean()
    return t, np.sqrt(t * (1.0 - t) / n)
