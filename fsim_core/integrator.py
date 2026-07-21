"""Module E -- Integrator: assembles g2(T; design), solves the master ceiling
for T_c, and runs one-at-a-time sensitivities.

F-series governing math:
  * background law      g2 = 1 - rho^2 (1 - eps)
  * master ceiling      rho(Tc)^2 [1 - eps(Tc)] = 1/2   (i.e. g2(Tc) = 0.5)
  * rho(T) from the paper's two-channel Arrhenius retention with a coupled
    background (escaped carriers re-emit in the WL/SSL):
        S(T)   = 1 / (1 + a_esc e^{-E_a/kT} + b_p e^{-E_b/kT})   (retention)
        B(T)   = b0 + beta (1 - S(T))                            (background)
        rho(T) = S / (S + B)                                     (signal fraction)
    E_a: WL escape (dominant, high T); E_b: p-shell channel (weak, low T).
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import brentq

from .loading import b_injection, f1b_g2
from .spectral import KB, epsilon, gamma_of_T


def retention(T, a_esc, E_a, b_p, E_b):
    T = np.asarray(T, dtype=float)
    kT = KB * T
    return 1.0 / (1.0 + a_esc * np.exp(-E_a / kT) + b_p * np.exp(-E_b / kT))


def rho_of_T(T, a_esc, E_a, b_p, E_b, b0, beta):
    S = retention(T, a_esc, E_a, b_p, E_b)
    B = b0 + beta * (1.0 - S)
    return S / (S + B)


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


RHO_KEYS = ("a_esc", "E_a", "b_p", "E_b", "b0", "beta")
PARAM_KEYS = ("delta_xx", "gamma0", "a_ac", "b_lo", "E_lo", "w", "dx") + RHO_KEYS


def g2_of_T(T, p: dict) -> ModelPoint:
    """Assemble one model point from a flat parameter dict (meV/K units).

    Filter geometry: p["w"] window full width, p["dx"] X offset from window
    center. For temperature sweeps with T-dependent windows, callables are
    accepted for w and dx.

    Drive realism (Module D, optional keys):
      p["mu"]        mean cap-2 loading -> F1b replaces the F1 identity
      p["I"], p["channels"]  injection current + BackgroundChannel set ->
                     b_e(I,T) added to the background before rho
    """
    w = p["w"](T) if callable(p["w"]) else p["w"]
    dx = p.get("dx", 0.0)
    dx = dx(T) if callable(dx) else dx
    gam = float(gamma_of_T(T, p["gamma0"], p["a_ac"], p["b_lo"], p["E_lo"]))
    spec = epsilon(p["delta_xx"], gam, p.get("gamma_xx", gam),
                   w=w, kappa=p.get("kappa"), dx=dx)

    S = float(retention(T, p["a_esc"], p["E_a"], p["b_p"], p["E_b"]))
    B = p["b0"] + p["beta"] * (1.0 - S)
    if p.get("channels") and p.get("I") is not None:
        B += float(b_injection(p["channels"], p["I"], T))
    rho = S / (S + B)

    mu = p.get("mu")
    g2_dot = float(f1b_g2(mu, spec.eps)) if mu else spec.eps
    return ModelPoint(T=float(T), gamma=gam, eps=spec.eps, rho=rho,
                      g2=g2_from(g2_dot, rho), t_x=spec.t_x)


def g2_electrical(T_hs, I, V, p: dict, a, stack, duty=1.0) -> ModelPoint:
    """Electrical-separation theorem (F-series): electrical drive differs from
    optical only through Delta T_J (junction heating, Module A) and Delta rho
    (injection background channels, Module D). The eps path is untouched.

    Evaluates the device at T_j = T_hs + duty * I * V * R_th and returns the
    ModelPoint at the junction temperature (point.T == T_j)."""
    from .thermal import t_junction

    Tj = t_junction(duty * I * V, a, stack, T_hs)
    return g2_of_T(Tj, {**p, "I": I})


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
    scalar parameter (first pass before Sobol, per the planning doc)."""
    base = solve_Tc(p)
    out = {}
    for k in keys or [k for k in PARAM_KEYS if k in p and not callable(p[k])]:
        q = dict(p)
        q[k] = p[k] * (1.0 + rel) if p[k] != 0 else rel
        out[k] = solve_Tc(q) - base
    return out
