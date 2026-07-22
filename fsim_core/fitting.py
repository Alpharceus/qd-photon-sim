"""Phase-0/1 fit driver: joint fit of {eps(T), rho(T)} to the validation card.

Residual blocks (each optional block activates when the card carries it):
  1. g2_vs_T      -- the published g2(0) series; rows carry {w, dx} filter
                     geometry and an optional bound:"upper" flag (one-sided).
  2. gamma_vs_T   -- digitized Gamma_X(T) [and Gamma_XX(T)] linewidths (S2);
                     constrains the Gamma model and the XX/X width ratio r_xx.
  3. tau_vs_T     -- published X decay times (S4). tau(T)/tau(T0) is an
                     independent measure of the retention ratio S(T)/S(T0)
                     (assumes Gamma_r const -- assumption ledger).
Anchors/priors: Gamma(250 K) anchor (only when no gamma_vs_T data), E_a prior,
log10 b_p prior. V-a criterion: +-0.03 per g2 point (planning doc section 1).
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.optimize import least_squares

from .card import Card, Tag, widest
from .integrator import g2_of_T, retention
from .spectral import gamma_of_T

BASE_FREE = ("gamma0", "a_ac", "b_lo", "a_esc", "E_a", "b_p", "b0", "beta")
LOG_FREE = {"a_esc", "b_p", "b0", "beta"}  # span decades: fitted in log10 space

V_A_TOL = 0.03          # per-point absolute tolerance on g2
GAMMA_RELERR = 0.15     # digitization-class relative error on linewidths
TAU_RELERR = 0.10       # relative error on lifetime ratios


@dataclass
class FitResult:
    params: dict                 # full flat parameter dict (fitted + fixed)
    data: list[dict]             # g2 rows: T, g2, err, w, dx, bound
    model_g2: np.ndarray
    residuals: np.ndarray        # g2 block: model - data (0 for satisfied bounds)
    passed: bool                 # all g2 |resid| <= V_A_TOL
    gamma_anchor: tuple          # (model Gamma(T_anchor), anchor, tol, T_anchor)
    tag: Tag
    cost: float = 0.0
    notes: list[str] = field(default_factory=list)
    aux: dict = field(default_factory=dict)   # per-block worst residuals etc.


def fit_phase0(card: Card, dataset="g2_vs_T", seed=0, n_starts=40) -> FitResult:
    # free set = every candidate the card declares as a range (swept, never averaged)
    CAND = BASE_FREE + ("r_xx", "E_lo")
    FREE = [k for k in CAND if k in card.params and card[k].is_free]

    def pack(values):
        return np.array([np.log10(values[k]) if k in LOG_FREE else values[k] for k in FREE])

    def unpack(x):
        return {k: (10.0**v if k in LOG_FREE else v) for k, v in zip(FREE, x)}

    ds = card.datasets[dataset]
    rows = sorted((dict(r) for r in ds.rows), key=lambda r: r["T"])
    for r in rows:
        r.setdefault("err", V_A_TOL)
        r.setdefault("bound", "value")

    gam_rows = sorted(card.datasets["gamma_vs_T"].rows, key=lambda r: r["T"]) \
        if "gamma_vs_T" in card.datasets else []
    tau_rows = sorted(card.datasets["tau_vs_T"].rows, key=lambda r: r["T"]) \
        if "tau_vs_T" in card.datasets else []

    fixed = {k: card[k].fixed for k in ("delta_xx", "E_lo", "E_b")
             if k in card.params and not card[k].is_free}
    if "mu_op" in card.params:
        fixed["mu"] = card["mu_op"].fixed  # F1b operating point (Module D)
    T_anchor = card["gamma_anchor_T"].fixed
    g_anchor = card["gamma_anchor"].fixed
    g_anchor_tol = card["gamma_anchor_tol"].fixed
    Ea_prior, Ea_tol = card["E_a_prior"].fixed, card["E_a_prior_tol"].fixed
    bp_prior, bp_tol = card["b_p_prior_log10"].fixed, card["b_p_prior_log10_tol"].fixed

    def model_g2_all(th):
        return np.array([
            g2_of_T(r["T"], {**fixed, **th, "w": r["w"], "dx": r["dx"]}).g2
            for r in rows
        ])

    def residuals(x):
        th = unpack(x)
        out = []
        m = model_g2_all(th)
        out.extend(
            max(0.0, mi - r["g2"]) / r["err"] if r["bound"] == "upper"
            else (mi - r["g2"]) / r["err"]
            for mi, r in zip(m, rows))
        gm = lambda T: gamma_of_T(T, th["gamma0"], th["a_ac"], th["b_lo"],
                                  th.get("E_lo", fixed.get("E_lo")))
        for r in gam_rows:
            out.append((gm(r["T"]) - r["gamma_x"]) / (GAMMA_RELERR * r["gamma_x"]))
            if "gamma_xx" in r and "r_xx" in th:
                out.append((th["r_xx"] * gm(r["T"]) - r["gamma_xx"])
                           / (GAMMA_RELERR * r["gamma_xx"]))
        if tau_rows:
            S = lambda T: retention(T, th["a_esc"], th["E_a"], th["b_p"], fixed["E_b"])
            T0, tau0 = tau_rows[0]["T"], tau_rows[0]["tau_ns"]
            S0 = S(T0)
            for r in tau_rows[1:]:
                ratio_d = r["tau_ns"] / tau0
                out.append((S(r["T"]) / S0 - ratio_d) / (TAU_RELERR * ratio_d))
        if not gam_rows:  # anchor only when no real linewidth data
            out.append((gm(T_anchor) - g_anchor) / g_anchor_tol)
        out.append((th["E_a"] - Ea_prior) / Ea_tol)
        out.append((np.log10(th["b_p"]) - bp_prior) / bp_tol)
        return np.array(out)

    lo = pack({k: card[k].bounds[0] for k in FREE})
    hi = pack({k: card[k].bounds[1] for k in FREE})

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

    th = unpack(best.x)
    p_full = {**fixed, **th}
    m = model_g2_all(th)
    resid = np.array([
        max(0.0, mi - r["g2"]) if r["bound"] == "upper" else mi - r["g2"]
        for mi, r in zip(m, rows)
    ])
    E_lo_fit = th.get("E_lo", fixed.get("E_lo"))
    g_at_anchor = float(gamma_of_T(T_anchor, th["gamma0"], th["a_ac"], th["b_lo"], E_lo_fit))

    aux = {}
    if gam_rows:
        gr = [abs(gamma_of_T(r["T"], th["gamma0"], th["a_ac"], th["b_lo"], E_lo_fit)
                  - r["gamma_x"]) / r["gamma_x"] for r in gam_rows]
        aux["gamma_worst_relerr"] = float(max(gr))
    if tau_rows:
        S = lambda T: retention(T, th["a_esc"], th["E_a"], th["b_p"], fixed["E_b"])
        S0, tau0 = S(tau_rows[0]["T"]), tau_rows[0]["tau_ns"]
        tr = [abs(S(r["T"]) / S0 - r["tau_ns"] / tau0) / (r["tau_ns"] / tau0)
              for r in tau_rows[1:]]
        aux["tau_worst_relerr"] = float(max(tr))

    tag = widest(
        ds.tag,
        *(card.datasets[d].tag for d in ("gamma_vs_T", "tau_vs_T") if d in card.datasets),
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
        tag=tag, cost=float(best.cost), notes=notes, aux=aux,
    )
