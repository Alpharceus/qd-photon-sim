"""Phase-0 fit driver: fit {eps(T), rho(T)} to the g2(T) validation series.

V-a criterion (planning doc, section 1): one physically constrained
parameterization must reproduce the six-point series within +-0.03 per point,
with Gamma(T) hitting the published >=250 K anchor and the extracted E_a
consistent with the published activation energy. Anchors and priors enter as
pseudo-residuals so the fit is constrained, not regularized after the fact.

Data rows carry their own filter geometry {w, dx} (the published windows move
and widen with T) and an optional bound: "upper" flag (the 78 K point is an
upper bound; it only penalizes the model from above).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

from .card import Card, Tag, widest
from .integrator import g2_of_T
from .spectral import gamma_of_T

FREE = ("gamma0", "a_ac", "b_lo", "a_esc", "E_a", "b_p", "b0", "beta")
LOG_FREE = {"a_esc", "b_p", "b0", "beta"}  # span decades: fitted in log10 space

V_A_TOL = 0.03  # per-point absolute tolerance, planning doc section 1


def _pack(values: dict) -> np.ndarray:
    return np.array([np.log10(values[k]) if k in LOG_FREE else values[k] for k in FREE])


def _unpack(x: np.ndarray) -> dict:
    return {k: (10.0**v if k in LOG_FREE else v) for k, v in zip(FREE, x)}


@dataclass
class FitResult:
    params: dict                 # full flat parameter dict (fitted + fixed)
    data: list[dict]             # T, g2, err, w, dx, bound rows
    model_g2: np.ndarray
    residuals: np.ndarray        # model - data (0 for satisfied upper bounds)
    passed: bool                 # all |resid| <= V_A_TOL
    gamma_anchor: tuple          # (model Gamma(T_anchor), anchor value, anchor tol, T_anchor)
    tag: Tag                     # widest tag in the input chain
    cost: float = 0.0
    notes: list[str] = field(default_factory=list)


def fit_phase0(card: Card, dataset="g2_vs_T", seed=0, n_starts=40) -> FitResult:
    ds = card.datasets[dataset]
    rows = sorted((dict(r) for r in ds.rows), key=lambda r: r["T"])
    for r in rows:
        r.setdefault("err", V_A_TOL)
        r.setdefault("bound", "value")

    fixed = {k: card[k].fixed for k in ("delta_xx", "E_lo", "E_b") if k in card.params}
    if "mu_op" in card.params:
        fixed["mu"] = card["mu_op"].fixed  # F1b operating point (Module D)
    T_anchor = card["gamma_anchor_T"].fixed
    g_anchor = card["gamma_anchor"].fixed
    g_anchor_tol = card["gamma_anchor_tol"].fixed
    Ea_prior, Ea_tol = card["E_a_prior"].fixed, card["E_a_prior_tol"].fixed
    bp_prior, bp_tol = card["b_p_prior_log10"].fixed, card["b_p_prior_log10_tol"].fixed

    def model_all(theta: dict) -> np.ndarray:
        return np.array([
            g2_of_T(r["T"], {**fixed, **theta, "w": r["w"], "dx": r["dx"]}).g2
            for r in rows
        ])

    def residuals(x: np.ndarray) -> np.ndarray:
        th = _unpack(x)
        m = model_all(th)
        r_data = np.array([
            max(0.0, mi - r["g2"]) / r["err"] if r["bound"] == "upper"
            else (mi - r["g2"]) / r["err"]
            for mi, r in zip(m, rows)
        ])
        g_at = gamma_of_T(T_anchor, th["gamma0"], th["a_ac"], th["b_lo"], fixed["E_lo"])
        return np.concatenate([r_data, [
            (g_at - g_anchor) / g_anchor_tol,
            (th["E_a"] - Ea_prior) / Ea_tol,
            (np.log10(th["b_p"]) - bp_prior) / bp_tol,
        ]])

    lo = _pack({k: card[k].bounds[0] for k in FREE})
    hi = _pack({k: card[k].bounds[1] for k in FREE})

    # multistart over the fit box: ranges are swept, never averaged
    rng = np.random.default_rng(seed)
    best = None
    for i in range(n_starts):
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
    resid = np.array([
        max(0.0, mi - r["g2"]) if r["bound"] == "upper" else mi - r["g2"]
        for mi, r in zip(m, rows)
    ])
    g_at_anchor = float(gamma_of_T(T_anchor, th["gamma0"], th["a_ac"], th["b_lo"], fixed["E_lo"]))

    tag = widest(
        ds.tag,
        *(card[k].tag for k in FREE),
        *(card[k].tag for k in ("delta_xx", "E_lo", "E_b") if k in card.params),
    )

    notes = []
    ph = card.placeholders()
    if ph:
        notes.append(
            "PLACEHOLDER inputs present (%s): fit exercises the machinery only -- "
            "the V-a gate CANNOT close on this run." % ", ".join(ph)
        )

    return FitResult(
        params=p_full, data=rows,
        model_g2=m, residuals=resid,
        passed=bool(np.all(np.abs(resid) <= V_A_TOL)),
        gamma_anchor=(g_at_anchor, g_anchor, g_anchor_tol, T_anchor),
        tag=tag, cost=float(best.cost), notes=notes,
    )
