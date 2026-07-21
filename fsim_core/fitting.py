"""Phase-0 fit driver: fit {eps(T), rho(T)} to the g2(T) validation series.

V-a criterion (planning doc, section 1): one physically constrained
parameterization must reproduce the six-point series within +-0.03 per point,
with Gamma(T) hitting the published >=250 K anchor and E_a consistent with the
published Arrhenius activation energy. Anchors enter as pseudo-residuals so the
fit is constrained, not merely regularized after the fact.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

from .card import Card, Tag, widest
from .integrator import g2_of_T
from .spectral import gamma_of_T

FREE = ("gamma0", "a_ac", "b_lo", "C_esc", "E_a", "b0")
LOG_FREE = {"C_esc", "b0"}  # fitted in log10 space (span decades)

V_A_TOL = 0.03  # per-point absolute tolerance, planning doc section 1


def _pack(values: dict) -> np.ndarray:
    return np.array([np.log10(values[k]) if k in LOG_FREE else values[k] for k in FREE])


def _unpack(x: np.ndarray) -> dict:
    return {k: (10.0**v if k in LOG_FREE else v) for k, v in zip(FREE, x)}


@dataclass
class FitResult:
    params: dict                 # full flat parameter dict (fitted + fixed)
    data: list[dict]             # T, g2, err rows
    model_g2: np.ndarray
    residuals: np.ndarray        # model - data, absolute
    passed: bool                 # all |resid| <= V_A_TOL
    gamma_anchor: tuple          # (model Gamma(T_anchor), anchor value, anchor tol)
    tag: Tag                     # widest tag in the input chain
    cost: float = 0.0
    notes: list[str] = field(default_factory=list)


def fit_phase0(card: Card, dataset="g2_vs_T", seed=0) -> FitResult:
    ds = card.datasets[dataset]
    rows = sorted(ds.rows, key=lambda r: r["T"])
    Ts = np.array([r["T"] for r in rows], float)
    g2s = np.array([r["g2"] for r in rows], float)
    errs = np.array([r.get("err", V_A_TOL) for r in rows], float)

    fixed = {k: card[k].fixed for k in ("delta_xx", "E_lo", "w") if k in card.params}
    T_anchor = card["gamma_anchor_T"].fixed
    g_anchor = card["gamma_anchor"].fixed
    g_anchor_tol = card["gamma_anchor_tol"].fixed
    Ea_prior = card["E_a_prior"].fixed
    Ea_tol = card["E_a_prior_tol"].fixed

    def model_all(theta: dict) -> np.ndarray:
        p = {**fixed, **theta}
        return np.array([g2_of_T(T, p).g2 for T in Ts])

    def residuals(x: np.ndarray) -> np.ndarray:
        th = _unpack(x)
        r_data = (model_all(th) - g2s) / errs
        g250 = gamma_of_T(T_anchor, th["gamma0"], th["a_ac"], th["b_lo"], fixed["E_lo"])
        r_anchor = (g250 - g_anchor) / g_anchor_tol
        r_ea = (th["E_a"] - Ea_prior) / Ea_tol
        return np.concatenate([r_data, [r_anchor, r_ea]])

    lo = _pack({k: card[k].bounds[0] for k in FREE})
    hi = _pack({k: card[k].bounds[1] for k in FREE})

    # multistart over the fit box: ranges are swept, never averaged
    rng = np.random.default_rng(seed)
    best = None
    for i in range(24):
        x0 = lo + (hi - lo) * (0.5 if i == 0 else rng.random(len(FREE)))
        try:
            sol = least_squares(residuals, x0, bounds=(lo, hi), method="trf")
        except Exception:
            continue
        if best is None or sol.cost < best.cost:
            best = sol

    th = _unpack(best.x)
    p_full = {**fixed, **th}
    m = model_all(th)
    resid = m - g2s
    g_at_anchor = float(gamma_of_T(T_anchor, th["gamma0"], th["a_ac"], th["b_lo"], fixed["E_lo"]))

    tag = widest(
        ds.tag,
        *(card[k].tag for k in FREE),
        *(card[k].tag for k in ("delta_xx", "E_lo", "w") if k in card.params),
    )

    notes = []
    ph = card.placeholders()
    if ph:
        notes.append(
            "PLACEHOLDER inputs present (%s): fit exercises the machinery only -- "
            "the V-a gate CANNOT close on this run." % ", ".join(ph)
        )

    return FitResult(
        params=p_full,
        data=[{"T": float(T), "g2": float(g), "err": float(e)} for T, g, e in zip(Ts, g2s, errs)],
        model_g2=m, residuals=resid,
        passed=bool(np.all(np.abs(resid) <= V_A_TOL)),
        gamma_anchor=(g_at_anchor, g_anchor, g_anchor_tol),
        tag=tag, cost=float(best.cost), notes=notes,
    )
