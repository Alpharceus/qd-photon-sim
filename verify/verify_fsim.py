"""FSIM regression suite (Yuki standard: units, limits, second methods).

Every closed form in Module C/E is checked against an independent method:
analytic limits, numeric quadrature, and Monte-Carlo sampling.
Run: python verify/verify_fsim.py   (exit code 0 iff all pass)
"""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core.integrator import g2_from, g2_of_T, rho_of_T, solve_Tc
from fsim_core.spectral import (
    KB,
    cavity_transmission,
    combined_transmission,
    epsilon,
    epsilon_narrow_filter,
    gamma_of_T,
    lorentzian,
    mc_transmission,
    optimal_width,
    tophat_transmission,
)

CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


# ---------------------------------------------------------------- top-hat closed form

@check("tophat: wide-window limit t->1")
def _():
    assert abs(tophat_transmission(0.0, 1.0, 1e6) - 1.0) < 1e-5


@check("tophat: narrow-window limit t ~ w * L(0)")
def _():
    w = 1e-6
    t = tophat_transmission(0.3, 1.0, w)
    assert abs(t - w * lorentzian(0.0, 0.3, 1.0)) < 1e-9


@check("tophat: detuning symmetry t(+d)=t(-d)")
def _():
    assert abs(tophat_transmission(1.7, 0.9, 2.3) - tophat_transmission(-1.7, 0.9, 2.3)) < 1e-14


@check("tophat: closed form vs numeric quadrature (20 random cases)")
def _():
    rng = np.random.default_rng(1)
    for _ in range(20):
        d, g, w = rng.uniform(0, 10), rng.uniform(0.05, 8), rng.uniform(0.1, 20)
        num, _err = quad(lorentzian, -w / 2, w / 2, args=(d, g), limit=200)
        assert abs(tophat_transmission(d, g, w) - num) < 1e-7, (d, g, w)


@check("tophat: MC second method (1e6 samples, 4 sigma)")
def _():
    for d, g, w in [(0.0, 1.0, 2.0), (9.0, 6.0, 2.0), (3.0, 0.5, 1.0)]:
        t, s = mc_transmission(d, g, w=w)
        assert abs(t - tophat_transmission(d, g, w)) < 4 * s + 1e-4, (d, g, w)


# ----------------------------------------------------------------- cavity closed form

@check("cavity: on-resonance value kappa/(gamma+kappa)")
def _():
    g, k = 1.3, 0.7
    assert abs(cavity_transmission(0.0, g, k) - k / (g + k)) < 1e-14


@check("cavity: closed form vs numeric quadrature (20 random cases)")
def _():
    rng = np.random.default_rng(2)
    for _ in range(20):
        d, g, k = rng.uniform(0, 10), rng.uniform(0.05, 8), rng.uniform(0.05, 8)
        num, _err = quad(
            lambda x: lorentzian(x, d, g) * (k / 2) ** 2 / (x**2 + (k / 2) ** 2),
            -np.inf, np.inf, limit=400,
        )
        assert abs(cavity_transmission(d, g, k) - num) < 1e-6, (d, g, k)


@check("cavity: MC second method")
def _():
    for d, g, k in [(0.0, 1.0, 1.0), (5.0, 2.0, 0.5)]:
        t, s = mc_transmission(d, g, kappa=k)
        assert abs(t - cavity_transmission(d, g, k)) < 4 * s + 1e-4


@check("combined tophat*cavity: numeric vs MC")
def _():
    d, g, w, k = 2.0, 1.0, 3.0, 1.5
    t, s = mc_transmission(d, g, w=w, kappa=k)
    assert abs(t - combined_transmission(d, g, w, k)) < 4 * s + 1e-4


# ---------------------------------------------------------------------------- epsilon

@check("epsilon: F1 identity eps=1 at zero splitting, equal widths")
def _():
    assert abs(epsilon(0.0, 1.0, 1.0, w=2.0).eps - 1.0) < 1e-14


@check("epsilon: narrow-filter bound is the w->0 limit")
def _():
    d, gx, gxx = 4.0, 1.0, 1.5
    r = epsilon(d, gx, gxx, w=1e-7)
    assert abs(r.eps - epsilon_narrow_filter(d, gx, gxx)) < 1e-6


@check("epsilon: monotone non-decreasing in w (leakage grows with window)")
def _():
    d, g = 6.0, 2.0
    ws = np.linspace(0.05, 30, 200)
    e = [epsilon(d, g, g, w=w).eps for w in ws]
    assert np.all(np.diff(e) > -1e-12)


@check("optimal_width: eps at w_opt equals the target")
def _():
    d, g = 9.0, 2.0
    tgt = 0.05
    w = optimal_width(d, g, g, tgt)
    assert abs(epsilon(d, g, g, w=w).eps - tgt) < 1e-9


# ------------------------------------------------------------------------ Module E

@check("Gamma(T): T->0 recovers Gamma0; monotone increasing")
def _():
    g = gamma_of_T(np.array([1e-3, 50, 150, 300]), 0.5, 1e-3, 20.0, 36.6)
    assert abs(g[0] - (0.5 + 1e-3 * 1e-3)) < 1e-12  # LO term frozen out at T->0
    assert np.all(np.diff(g) > 0)


@check("rho(T): frozen-escape limit 1/(1+b0); high-T limit below it")
def _():
    b0 = 0.02
    assert abs(rho_of_T(1.0, 1e7, 265.0, b0) - 1 / (1 + b0)) < 1e-12
    assert rho_of_T(400.0, 1e7, 265.0, b0) < 1 / (1 + b0)


@check("background law identities: rho=1 -> g2=eps; eps=1 -> g2=1; eps=0 -> 1-rho^2")
def _():
    assert abs(g2_from(0.3, 1.0) - 0.3) < 1e-14
    assert abs(g2_from(1.0, 0.7) - 1.0) < 1e-14
    assert abs(g2_from(0.0, 0.9) - (1 - 0.81)) < 1e-14


@check("master ceiling: g2(Tc) = 0.5 at the solved Tc")
def _():
    p = dict(delta_xx=9.0, gamma0=0.5, a_ac=2e-3, b_lo=25.0, E_lo=36.6,
             w=2.0, C_esc=1e7, E_a=265.0, b0=0.01)
    Tc = solve_Tc(p)
    assert np.isfinite(Tc) and abs(g2_of_T(Tc, p).g2 - 0.5) < 1e-6


@check("fit round-trip: synthetic g2(T) from known params recovered within +-0.03")
def _():
    from fsim_core.card import Card, DataSet, Param, Tag
    from fsim_core.fitting import fit_phase0

    truth = dict(delta_xx=9.0, gamma0=0.5, a_ac=2e-3, b_lo=25.0, E_lo=36.6,
                 w=2.0, C_esc=1e7, E_a=265.0, b0=0.01)
    Ts = [78.0, 110.0, 140.0, 170.0, 200.0, 230.0]
    rows = [{"T": T, "g2": g2_of_T(T, truth).g2, "err": 0.03} for T in Ts]

    def P(name, **kw):
        kw.setdefault("unit", "meV"); kw.setdefault("tag", Tag.A); kw.setdefault("source", "synthetic")
        return Param(name=name, **kw)

    card = Card(
        name="synthetic", meta={"name": "synthetic"},
        params={
            "delta_xx": P("delta_xx", value=9.0), "w": P("w", value=2.0),
            "E_lo": P("E_lo", value=36.6),
            "gamma0": P("gamma0", range=(0.05, 2.0)),
            "a_ac": P("a_ac", range=(0.0, 0.02), unit="meV/K"),
            "b_lo": P("b_lo", range=(1.0, 80.0)),
            "C_esc": P("C_esc", range=(1e4, 1e10), unit="1"),
            "E_a": P("E_a", range=(200.0, 330.0)),
            "b0": P("b0", range=(1e-4, 0.2), unit="1"),
            "gamma_anchor": P("gamma_anchor", value=6.78),
            "gamma_anchor_tol": P("gamma_anchor_tol", value=0.5),
            "gamma_anchor_T": P("gamma_anchor_T", value=250.0, unit="K"),
            "E_a_prior": P("E_a_prior", value=265.0),
            "E_a_prior_tol": P("E_a_prior_tol", value=30.0),
        },
        datasets={"g2_vs_T": DataSet("g2_vs_T", Tag.A, "synthetic", ["T", "g2", "err"], rows)},
    )
    fit = fit_phase0(card)
    assert fit.passed, f"max |resid| = {np.max(np.abs(fit.residuals)):.4f}"


@check("units: KB*300K ~ 25.85 meV (meV/K convention holds)")
def _():
    assert abs(KB * 300.0 - 25.852) < 0.01


def main():
    failed = 0
    for name, fn in CHECKS:
        try:
            fn()
            print(f"  PASS  {name}")
        except AssertionError as e:
            failed += 1
            print(f"* FAIL  {name}  {e}")
    n = len(CHECKS)
    print(f"\n{n - failed}/{n} checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
