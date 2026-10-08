"""Explain endpoints for FSIM Studio (ports of the Streamlit app.py explainers).

Three-layer rule: every number here comes from an fsim_core call
(spectral.epsilon / lorentzian / gamma_of_T, loading.gamma_eff / f1b_g2 /
loading_probs, integrator.g2_of_T / g2_from) or from a committed out/ file
(out/phase0/fit_*.csv, fit_params.json). The browser only plots what this
module returns; it never evaluates a model.

Routes (registered by register(app) from server.create_app):
  GET /api/explain/spectral?T=&w=&dx=&mu=&mode=PL|EL&dg_inj=&I_ratio=&pub=0|1
  GET /api/explain/cascade?T=
  GET /api/explain/va_fit

The parameterization is the one app.py uses for the V-a card: the saved fit
bundle out/phase0/fit_params.json (params_source says so), and for the
published window the card's g2_vs_T dataset (w, dx interpolated in T).
No physics fallbacks: parameters resolve through
fsim_core.scene_support.explain_fit_params; a parameter missing from the
bundle is a 400 {error, missing: [...]}. p_inj (not fitted) resolves to
fsim_core.loading.gamma_eff's own default and `param_sources` says so.
"""
from __future__ import annotations

import csv
import json

import numpy as np
from flask import Response, request

from . import CARDS, OUT
from .serialize import dumps

CARD_PATH = CARDS / "chatzarakis.yaml"
FIT_JSON = OUT / "phase0" / "fit_params.json"
FIT_POINTS = OUT / "phase0" / "fit_data_points.csv"
FIT_CURVE = OUT / "phase0" / "fit_model_curve.csv"

# slider domains copied from app.py's spectral explainer (input validation only)
DOMAIN = {"T": (4.0, 320.0), "w": (0.1, 15.0), "dx": (-6.0, 6.0), "mu": (0.0, 5.0),
          "dg_inj": (0.0, 20.0), "I_ratio": (0.0, 20.0)}
N_GRID = 1200


class ExplainError(Exception):
    def __init__(self, message, status=400, **extra):
        super().__init__(message)
        self.message = message
        self.status = status
        self.extra = extra


def _fit_params(required):
    """Fit-bundle parameters through fsim_core.scene_support.explain_fit_params:
    no silent physics fallbacks; a missing key is a 400 naming the list."""
    from fsim_core.scene_support import MissingFitParams, explain_fit_params
    try:
        return explain_fit_params(_fit()["params"], required)
    except MissingFitParams as exc:
        raise ExplainError(f"fit bundle {FIT_JSON.relative_to(OUT.parent).as_posix()} lacks {exc.missing}",
                           missing=exc.missing)


def _card():
    from fsim_core.card import load_card
    return load_card(CARD_PATH)


def _fit():
    return json.loads(FIT_JSON.read_text(encoding="utf-8"))


def _window_fns(card):
    rows = sorted(card.datasets["g2_vs_T"].rows, key=lambda r: r["T"])
    Ts = np.array([r["T"] for r in rows])
    ws = np.array([r["w"] for r in rows])
    dxs = np.array([r["dx"] for r in rows])
    return (lambda T: float(np.interp(T, Ts, ws))), (lambda T: float(np.interp(T, Ts, dxs)))


def _tag_name(tag) -> str:
    return getattr(tag, "name", str(tag))


def _arg(name, default):
    v = request.args.get(name)
    if v in (None, ""):
        return float(default)
    try:
        f = float(v)
    except ValueError:
        raise ExplainError(f"{name} must be a number")
    if not np.isfinite(f):
        raise ExplainError(f"{name} must be finite")
    lo, hi = DOMAIN.get(name, (-np.inf, np.inf))
    if f < lo or f > hi:
        raise ExplainError(f"{name}={f} is outside the explainer range [{lo}, {hi}]")
    return f


def spectral(T=78.0, w=2.0, dx=0.0, mu=None, mode="PL", dg_inj=2.0, I_ratio=1.0, pub=True):
    """app.py spectral tab: X/XX Lorentzians through a top-hat filter window."""
    from fsim_core.card import Tag, widest
    from fsim_core.loading import f1b_g2, gamma_eff
    from fsim_core.spectral import epsilon, gamma_of_T, lorentzian

    from fsim_core.scene_support import EXPLAIN_SPECTRAL_KEYS

    card = _card()
    resolved = _fit_params(EXPLAIN_SPECTRAL_KEYS)
    p, psrc = resolved["params"], resolved["sources"]
    w_of_T, dx_of_T = _window_fns(card)
    if pub:
        w, dx = w_of_T(T), dx_of_T(T)
    if mu is None:
        mu = float(p["mu"])
    if mode not in ("PL", "EL"):
        raise ExplainError("mode must be PL or EL")
    if mode == "PL":
        dg_inj, I_ratio = 0.0, 1.0
    delta = p["delta_xx"]
    p_inj = p["p_inj"]
    gam_pl = float(gamma_of_T(T, p["gamma0"], p["a_ac"], p["b_lo"], p["E_lo"]))
    gam_el = float(gamma_eff(gam_pl, dg_inj, I_ratio, p_inj))
    gam = gam_el if mode == "EL" else gam_pl
    spec = epsilon(delta, gam, gam, w=w, dx=dx)
    g2_dot = float(f1b_g2(mu, spec.eps)) if mu > 0 else float(spec.eps)

    x = np.linspace(-delta - 6 * gam - 4, 6 * gam + 4, N_GRID)
    LX = lorentzian(x, 0.0, gam)
    LXX = lorentzian(x, -delta, gam)
    peak = LX.max()
    cen = -dx
    inw = (x >= cen - w / 2) & (x <= cen + w / 2)
    leak_label = None
    if inw.any():
        k = int(np.argmax(np.where(inw, LXX, -np.inf)))
        leak_label = {"x": float(x[k]), "y": float(LXX[k] / peak)}
    delta_tag = card.params["delta_xx"].tag if "delta_xx" in card.params else Tag.A
    chain = widest(delta_tag, Tag.A)  # linewidth params are fitted [A] (app.py)
    return {
        "T": T, "mode": mode, "w": w, "dx": dx, "mu": mu, "published_window": bool(pub),
        "dg_inj": dg_inj, "I_ratio": I_ratio, "p_inj": p_inj,
        "delta_xx": delta, "gamma_pl": gam_pl, "gamma_el": gam_el, "gamma": gam,
        "gamma_ratio": (gam_el / gam_pl) if gam_pl > 0 else None,
        "eps": float(spec.eps), "t_x": float(spec.t_x), "t_xx": float(spec.t_xx),
        "g2_dot": g2_dot, "g2_dot_formula": "F1: g2_0 = eps (mu = 0)" if mu <= 0 else "F1b: g2_0 = f1b_g2(mu, eps)",
        "curves": {"x": x.tolist(), "X": (LX / peak).tolist(), "XX": (LXX / peak).tolist()},
        "window": {"center": cen, "lo": cen - w / 2, "hi": cen + w / 2},
        "leakage": {"x": x[inw].tolist(), "y": (LXX / peak)[inw].tolist()},
        "eps_label": leak_label,
        "tags": {"chain": _tag_name(chain), "delta_xx": _tag_name(delta_tag), "linewidth": "A",
                 "window": "V" if pub else "A", "mu": "A"},
        "params_source": f"saved bundle ({FIT_JSON.relative_to(OUT.parent).as_posix()}); V-a card cards/chatzarakis.yaml",
        "param_sources": psrc,
        "calls": ["fsim_core.spectral.gamma_of_T", "fsim_core.loading.gamma_eff",
                  "fsim_core.spectral.epsilon", "fsim_core.spectral.lorentzian", "fsim_core.loading.f1b_g2"],
    }


def cascade(T=150.0):
    """app.py cascade tab: one integrator.g2_of_T point plus cap-2 loading probabilities."""
    from fsim_core.integrator import g2_from, g2_of_T
    from fsim_core.loading import f1b_g2, loading_probs

    from fsim_core.scene_support import EXPLAIN_CASCADE_KEYS

    card = _card()
    _fit_params(EXPLAIN_CASCADE_KEYS)  # 400 with the missing list before any call
    p = dict(_fit()["params"])
    w_of_T, dx_of_T = _window_fns(card)
    p_full = {**p, "w": w_of_T(T), "dx": dx_of_T(T)}
    pt = g2_of_T(T, p_full)
    mu = p["mu"]  # a fitted mu of 0 stays 0 (no `or` fallback)
    P0, P1, P2 = loading_probs(mu)
    g2_dot = float(f1b_g2(p_full["mu"], pt.eps)) if p_full.get("mu") else float(pt.eps)
    return {
        "T": T, "w": p_full["w"], "dx": p_full["dx"], "mu": mu,
        "g2": float(pt.g2), "rho": float(pt.rho), "rho2": float(pt.rho) ** 2, "eps": float(pt.eps),
        "g2_dot": g2_dot, "g2_from_check": float(g2_from(g2_dot, pt.rho)),
        "t_x": float(pt.t_x), "t_xx": float(pt.eps * pt.t_x), "gamma": float(pt.gamma),
        "P0": float(P0), "P1": float(P1), "P2": float(P2),
        "formula": "g2(0) = 1 - rho^2 (1 - g2_0), g2_0 = F1b(mu, eps); at mu -> 0, g2_0 = eps",
        "tags": {"chain": (_fit().get("tag_chain") or "[A]").strip("[]"), "window": "V"},
        "params_source": f"saved bundle ({FIT_JSON.relative_to(OUT.parent).as_posix()}); window from cards/chatzarakis.yaml g2_vs_T",
        "calls": ["fsim_core.integrator.g2_of_T", "fsim_core.loading.loading_probs",
                  "fsim_core.loading.f1b_g2", "fsim_core.integrator.g2_from"],
    }


def _read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return v


def va_fit():
    """The committed V-a calibration bundle, verbatim (numbers parsed, never refit)."""
    pts = _read_csv(FIT_POINTS)
    cur = _read_csv(FIT_CURVE)
    fit = _fit()
    return {
        "points": [{k: _num(v) for k, v in r.items()} for r in pts],
        "curve": {k: [_num(r[k]) for r in cur] for k in (cur[0].keys() if cur else [])},
        "params": fit.get("params", {}), "Tc_K": fit.get("Tc_K"), "tag_chain": fit.get("tag_chain"),
        "passed": fit.get("passed_pm0.03"), "cost": fit.get("cost"), "notes": fit.get("notes", []),
        "files": [p.relative_to(OUT.parent).as_posix() for p in (FIT_POINTS, FIT_CURVE, FIT_JSON)],
        "label": "calibration (C), not held-out",
    }


def register(app):
    def _out(obj, status=200):
        return Response(dumps(obj), status=status, mimetype="application/json")

    @app.errorhandler(ExplainError)
    def _explain_error(exc):
        return _out({"error": exc.message, **exc.extra}, exc.status)

    @app.get("/api/explain/spectral")
    def explain_spectral():
        pub = request.args.get("pub", "1") not in ("0", "false", "False")
        mu = request.args.get("mu")
        return _out(spectral(
            T=_arg("T", 78.0), w=_arg("w", 2.0), dx=_arg("dx", 0.0),
            mu=None if mu in (None, "") else _arg("mu", 0.0),
            mode=(request.args.get("mode") or "PL").upper(),
            dg_inj=_arg("dg_inj", 2.0), I_ratio=_arg("I_ratio", 1.0), pub=pub))

    @app.get("/api/explain/cascade")
    def explain_cascade():
        return _out(cascade(T=_arg("T", 150.0)))

    @app.get("/api/explain/va_fit")
    def explain_va_fit():
        return _out(va_fit())

    return app


__all__ = ["register", "spectral", "cascade", "va_fit", "ExplainError"]
