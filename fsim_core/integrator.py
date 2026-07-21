"""Module E -- Integrator: assembles g2(T; design), solves the master ceiling
for T_c, and runs one-at-a-time sensitivities.

F-series governing math:
  * background law      g2 = 1 - rho^2 (1 - eps)
  * master ceiling      rho(Tc)^2 [1 - eps(Tc)] = 1/2   (i.e. g2(Tc) = 0.5)
  * rho(T) from Arrhenius retention against a constant background channel:
        S(T)   = 1 / (1 + C_esc * exp(-E_a / kB T))     (signal retention)
        rho(T) = S / (S + b0)                           (signal fraction)
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .spectral import KB, epsilon, gamma_of_T


def rho_of_T(T, C_esc, E_a, b0):
    T = np.asarray(T, dtype=float)
    S = 1.0 / (1.0 + C_esc * np.exp(-E_a / (KB * T)))
    return S / (S + b0)


def g2_from(eps, rho):
    """Background law (F-series)."""
    return 1.0 - rho**2 * (1.0 - eps)


@dataclass
class ModelPoint:
    T: float
    gamma: float
    eps: float
    rho: float
    g2: float
    t_x: float


# parameter dict keys expected by the assembler
PARAM_KEYS = ("delta_xx", "gamma0", "a_ac", "b_lo", "E_lo", "w", "C_esc", "E_a", "b0")


def g2_of_T(T, p: dict) -> ModelPoint:
    """Assemble one model point from a flat parameter dict (values in meV/K)."""
    gam = float(gamma_of_T(T, p["gamma0"], p["a_ac"], p["b_lo"], p["E_lo"]))
    spec = epsilon(p["delta_xx"], gam, p.get("gamma_xx", gam), w=p["w"], kappa=p.get("kappa"))
    rho = float(rho_of_T(T, p["C_esc"], p["E_a"], p["b0"]))
    return ModelPoint(T=float(T), gamma=gam, eps=spec.eps, rho=rho,
                      g2=g2_from(spec.eps, rho), t_x=spec.t_x)


def g2_curve(Ts, p: dict) -> list[ModelPoint]:
    return [g2_of_T(T, p) for T in np.asarray(Ts, dtype=float)]


def solve_Tc(p: dict, T_lo=4.0, T_hi=400.0):
    """Master ceiling: T where rho^2(1-eps) drops to 1/2 (g2 crosses 0.5).
    Returns np.nan if the ceiling is not crossed inside [T_lo, T_hi]."""
    f = lambda T: g2_of_T(T, p).g2 - 0.5
    if f(T_lo) >= 0:
        return T_lo  # already above the ceiling at the coldest point
    if f(T_hi) < 0:
        return np.nan
    return brentq(f, T_lo, T_hi, xtol=1e-3)


def oat_sensitivity(p: dict, rel=0.05, keys=None) -> dict[str, float]:
    """One-at-a-time sensitivity of T_c: dTc for a +rel fractional bump of each
    parameter (first pass before Sobol, per the planning doc)."""
    base = solve_Tc(p)
    out = {}
    for k in keys or [k for k in PARAM_KEYS if k in p]:
        q = dict(p)
        q[k] = p[k] * (1.0 + rel) if p[k] != 0 else rel
        out[k] = solve_Tc(q) - base
    return out
