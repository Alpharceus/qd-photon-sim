"""Module B -- Cavity/EM, analytic tier (F6).

A cavity of linewidth kappa acts as (F6):
  (i)   a filter: Lorentzian acceptance multiplying the line transmissions
        (implemented in spectral.transmission via kappa) -- XX suppressed by
        ~((Gamma+kappa)/2)^2/Delta^2 at large Delta;
  (ii)  a collection/rate boost -> rho (integrator collection_gain);
  (iii) selective Purcell limited by spectral overlap: F_eff = F_P kappa/(kappa+Gamma);
  (iv)  the mode-tracking design rule: the emitter Varshni-redshifts with T while
        the cavity mode moves only via dn/dT -- place the mode red of the
        cryogenic line so that detuning vanishes at T_target.

The SiN evanescent-coupling variant (beta) is brightness-only by Lemma 1: the
schema forbids beta from entering the g2 path (enforced by a regression check).

The numeric kappa/F_P/collection-gain values for a real design come from
MEEP/COMSOL [A -> computed]; this module supplies the design rules and the
closed forms they plug into.
"""
from __future__ import annotations

import numpy as np

from .spectral import epsilon

# Varshni parameters, E(T) = E0 - alpha T^2/(T+beta). The (211)B InAs/GaAs dots
# follow bulk GaAs closely (Chatzarakis S1(a) [V]).
VARSHNI = {
    "GaAs": (5.405e-4, 204.0),   # eV/K, K
    "InP": (3.63e-4, 162.0),
}


def varshni_shift(T, material="GaAs"):
    """Emitter redshift E(T)-E(0) in eV (negative)."""
    a, b = VARSHNI[material]
    T = np.asarray(T, dtype=float)
    return -a * T**2 / (T + b)


def emitter_energy(T, E0, material="GaAs"):
    return E0 + varshni_shift(T, material)


def cavity_mode_energy(T, E_cav0, dEdT=-4.0e-5):
    """Cavity mode vs T: refractive-index drift only (dE/dT ~ -0.04 meV/K for
    GaAs-based cavities [E]) -- much slower than the emitter Varshni shift."""
    return E_cav0 + dEdT * np.asarray(T, dtype=float)


def tracking_detuning(T, E_X0, T_target, material="GaAs", dEdT=-4.0e-5,
                      E_cav0=None):
    """Detuning delta(T) = E_X(T) - E_cav(T) in eV under the mode-tracking rule:
    by default the mode is placed so delta(T_target) = 0 (red of the cryogenic
    line by the Varshni shift, F6 iv)."""
    if E_cav0 is None:
        E_cav0 = (emitter_energy(T_target, E_X0, material)
                  - dEdT * T_target)
    return emitter_energy(T, E_X0, material) - cavity_mode_energy(T, E_cav0, dEdT)


def purcell_eff(F_P, kappa, gamma):
    """F6 iii: spectral-overlap-limited effective Purcell factor."""
    return F_P * kappa / (kappa + gamma)


def cavity_eps(T, p, kappa, detuning_meV):
    """epsilon with the cavity acceptance centered at detuning_meV from X
    (positive = X above the mode). Uses the spectral Module C kernel: the
    cavity Lorentzian replaces/multiplies the slit."""
    from .spectral import gamma_of_T

    gam = float(gamma_of_T(T, p["gamma0"], p["a_ac"], p["b_lo"], p["E_lo"]))
    gam_xx = p.get("r_xx", 1.0) * gam
    w = p.get("w")  # optional additional slit
    return epsilon(p["delta_xx"], gam, gam_xx, w=w, kappa=kappa,
                   dx=detuning_meV)


def sin_beta_brightness(t_x, beta):
    """Lemma 1: SiN evanescent coupling multiplies brightness only. Returns the
    waveguide-coupled brightness; deliberately returns nothing g2-related."""
    return beta * t_x
