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

from .loading import b_injection, f1b_g2, f8_g2, f8b_thin_fano
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
      p["mu"]        mean cap-2 loading -> F1b (or F8, see below) replaces
                     the F1 identity
      p["I"], p["channels"]  injection current + BackgroundChannel set ->
                     b_e(I,T) added to the background before rho
      p["dg_inj"], p["p_inj"], p["I_ref"]  F5' injection broadening (v1.1):
                     when p["dg_inj"] and p["I"] are both given, the X/XX
                     linewidth becomes Gamma_eff = Gamma(T) + dg_inj *
                     (I/I_ref)^p_inj BEFORE epsilon (and hence before any
                     auto-w filter width derived from Gamma downstream, e.g.
                     device.evaluate's filter.auto_w). p_inj defaults to 1.0,
                     I_ref defaults to 1.0.
      p["F_p"], p["eta_capture"]  F8/F8b (v1.1): a pump Fano factor F_p != 1.0
                     switches g2_dot from f1b_g2(mu, eps) (Poisson cap-2) to
                     f8_g2(mu, F_eff, eps), the moment-matched (mu, Fano)
                     cap-2 family, with F_eff = f8b_thin_fano(eta_capture,
                     F_p) (F8b dot-capture thinning; eta_capture defaults to
                     1.0, i.e. F_eff = F_p). F_p == 1.0 (the default) keeps
                     the f1b_g2 path EXACTLY -- the two cap-2 conventions are
                     not identical at finite mu (see loading.py module
                     docstring), so every pre-v1.1 result is preserved
                     bit-for-bit unless F_p is explicitly set away from 1.0.
    """
    w = p["w"](T) if callable(p["w"]) else p["w"]
    dx = p.get("dx", 0.0)
    dx = dx(T) if callable(dx) else dx
    gam = float(gamma_of_T(T, p["gamma0"], p["a_ac"], p["b_lo"], p["E_lo"]))
    if p.get("dg_inj") and p.get("I") is not None:
        gam = gam + p["dg_inj"] * (p["I"] / p.get("I_ref", 1.0)) ** p.get("p_inj", 1.0)
    gam_xx = p.get("gamma_xx", p.get("r_xx", 1.0) * gam)  # S2: Gamma_XX < Gamma_X
    spec = epsilon(p["delta_xx"], gam, gam_xx,
                   w=w, kappa=p.get("kappa"), dx=dx)

    S = float(retention(T, p["a_esc"], p["E_a"], p["b_p"], p["E_b"]))
    B = p["b0"] + p["beta"] * (1.0 - S)
    if p.get("channels") and p.get("I") is not None:
        B += float(b_injection(p["channels"], p["I"], T))
    G = p.get("collection_gain", 1.0)  # cavity/waveguide gain acts on signal only (F6 ii)
    rho = G * S / (G * S + B)

    mu = p.get("mu")
    F_p = p.get("F_p", 1.0)
    if mu:
        if F_p != 1.0:
            F_eff = f8b_thin_fano(p.get("eta_capture", 1.0), F_p)
            g2_dot = float(f8_g2(mu, F_eff, spec.eps))
        else:
            g2_dot = float(f1b_g2(mu, spec.eps))
    else:
        g2_dot = spec.eps
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
