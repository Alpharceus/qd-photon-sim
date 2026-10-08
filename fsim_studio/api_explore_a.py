"""Model explorers A for FSIM Studio (studio-p2d, Coder E1): cw_g2, pulse_counting, lindblad.

Three-layer rule: every number here is the return value of a public fsim_core
call (cw_g2.cw_report / g2_cw_zero / g2_cw_zero_low_pump, pulse_counting.pulse_g2 /
deterministic_cycle_g2, loading.f1b_g2, lindblad.build_system / g2_tau /
steady_state / indistinguishability / filtered_g2_zero / filtered_flux).  This
module only parses and range-checks inputs, picks grids, thins arrays for the
wire and attaches provenance tags and the pinned caveat text.  The browser plots
what it gets; it never evaluates a model.

All three explorers are LIVE (brief cost lines: 1.8 ms / 0.3 ms / 0.1-1 ms per
call), so no job-pool entry is needed.

Routes (registered by register(app) from server.create_app):
  GET /api/explore/a/controls                 control specs + pinned caveats (ranges, defaults, tags)
  GET /api/explore/cw_g2?r=&gamma_X=&eps=&rho=&irf_fwhm_ps=&irf_shape=&pump_ratio=&T=&tau_max=
  GET /api/explore/pulse_counting?r=&tau_on=&period=&gamma_X=&eps=&pump_ratio=&T=&gate=none|auto&eta_load=
  GET /api/explore/lindblad?panel=incoherent|rabi|hom|filter&...

T is optional in every route: omitted means no thermal escape (k_X = k_XX = 0).
"""
from __future__ import annotations

import json
import math
import time
import warnings

import numpy as np
from flask import Response, request

from . import OUT
from .serialize import dumps
from .tags import widest

FIT_JSON = OUT / "phase0" / "fit_params.json"
FIT_KEYS = ("a_esc", "E_a", "b_p", "E_b")
WIRE_MAX = 2001          # curve points sent to the browser (display decimation only)

# ------------------------------------------------------------------ pinned caveats (brief 1(d), 2(d), 6)
CAVEAT_CW = (
    "CW dc model, rate equations, incoherent pump only; not valid for resonant (Rabi) drive "
    "(cw_g2.py:18-23). g2_raw is what a 500 ps HBT measures, g2_dot is what the dot does; the "
    "designer's pulsed peak-area g2 is a different observable (cw_g2.py:439-502).")
CAVEAT_CW_BUNCH = (
    "g2_dot(0) > 1 is cascade bunching of the leaked XX partner, not a model failure.")
CAVEAT_PULSED = (
    "Rectangular pump [A]; incoherent capture, not a pi-pulse (:101-106). g2 is normalised to the "
    "long-delay peak, not the adjacent one (audit D2, :79-99); g2_adj is the adjacent-normalised "
    "value. The static f1b curve is exact only for an instantaneous pulse; the gap is the "
    "re-excitation the headline model sees.")
CAVEAT_PULSED_EXTRA = (
    "Background is not in this module (:122-146): the explorer shows dot-only g2, rho is a "
    "device.py addition. Deshpande 1.3 ns is [V abstract-only].")
CAVEAT_LINDBLAD = "Markovian baths, Lorentzian ZPL; sidebands are the separate qd_gf layer (lindblad.py:91-92)."
CAVEAT_LINDBLAD_EXTRA = (
    "This is the quantum cross-check: it covers what the rate modules refuse (coherent Rabi drive, "
    "cw_g2.py:18-23, pulse_counting.py:101-106), HOM indistinguishability and the cascaded filter.")


class ExploreAError(Exception):
    def __init__(self, message, status=400, **extra):
        super().__init__(message)
        self.message = message
        self.status = status
        self.extra = extra


# ------------------------------------------------------------------ control specs
# tag = the tag of the cited line of the brief / module for the DEFAULT; "scale" is a slider
# mapping hint for the browser (lin | log), not physics.
def _c(id_, label, lo, hi, default, unit="", tag="A", scale="lin", step=None, note=""):
    return {"id": id_, "label": label, "min": lo, "max": hi, "default": default, "unit": unit,
            "tag": tag, "scale": scale, "step": step, "note": note}


_T_CTL = _c("T", "Escape temperature T", 4.0, 350.0, 230.0, "K", "DR", "lin", 1.0,
            "feeds escape_rates_from_retention with the V-a fit [DR]; off = no escape")

CONTROLS = {
    "cw_g2": [
        _c("r", "Pump rate r", 1e-3, 100.0, 0.5, "1/ns", "A", "log", note="brief 1(b) default 0.5 [A]"),
        _c("gamma_X", "Radiative rate gamma_X", 0.2, 5.0, 1.0, "1/ns", "E", "lin", 0.05,
           "Deshpande tau 1.3 +- 0.3 ns [E] (pulse_counting.py:368-370); gamma_XX = 2 gamma_X [DR]"),
        _c("eps", "Leakage eps = t_XX / t_X", 0.0, 1.0, 0.1, "", "A", "lin", 0.005, "t_X = 1"),
        _c("rho", "Signal fraction rho", 0.5, 1.0, 0.9, "", "A", "lin", 0.005),
        _c("irf_fwhm_ps", "IRF FWHM", 0.0, 1000.0, 500.0, "ps", "V", "lin", 5.0,
           "Reischle 2008 exponential IRF (cw_g2.py:57-61)"),
        _c("pump_ratio", "Pump ratio p (X to XX)", 0.0, 1.0, 1.0, "", "A", "lin", 0.01),
        _c("tau_max", "Display range tau_max", 0.5, 20.0, 10.0, "ns", "A", "lin", 0.5),
        _T_CTL,
    ],
    "pulse_counting": [
        _c("r", "Pump rate r (readouts)", 0.01, 1e4, 10.0, "1/ns", "A", "log",
           note="Plot A sweeps r; the readouts are at this r"),
        _c("tau_on", "Pulse length tau_on", 0.01, 2.0, 0.1, "ns", "A", "log",
           note="RT cards drive.diode.tau_pulse_ns 0.1 [A]"),
        _c("period", "Period", 1.0, 50.0, 12.5, "ns", "A", "lin", 0.5, "12.5 ns [A]"),
        _c("gamma_X", "Radiative rate gamma_X", 0.2, 5.0, 1.0, "1/ns", "E", "lin", 0.05),
        _c("eps", "Leakage eps = t_XX / t_X", 0.0, 1.0, 0.1, "", "A", "lin", 0.005),
        _c("pump_ratio", "Pump ratio p (X to XX)", 0.0, 1.0, 1.0, "", "A", "lin", 0.01),
        _c("eta_load", "Deterministic load eta_load", 0.0, 1.0, 1.0, "", "A", "lin", 0.01,
           "pulse_counting.py:318-321 [A]"),
        _T_CTL,
    ],
    "lindblad": {
        "incoherent": [
            _c("r", "Pump rate r", 1e-3, 100.0, 0.5, "1/ns", "A", "log"),
            _c("gamma_X", "Radiative rate gamma_X", 0.2, 5.0, 1.0, "1/ns", "E", "lin", 0.05),
            _c("eps", "Leakage eps", 0.0, 1.0, 0.1, "", "A", "lin", 0.005),
            _c("pump_ratio", "Pump ratio p", 0.0, 1.0, 1.0, "", "A", "lin", 0.01),
            _c("tau_max", "Display range tau_max", 1.0, 20.0, 10.0, "ns", "A", "lin", 0.5),
            _T_CTL,
        ],
        "rabi": [
            _c("omega", "Rabi frequency Omega", 0.1, 20.0, 3.0, "1/ns", "A", "lin", 0.1, "verify (b3): Omega 3, gamma 1"),
            _c("gamma", "Decay rate gamma_1", 0.2, 5.0, 1.0, "1/ns", "A", "lin", 0.05),
            _c("deph", "Pure dephasing c", 0.0, 5.0, 0.0, "1/ns", "A", "lin", 0.05, "gamma_2 = gamma_1/2 + c/2"),
            _c("tau_max", "Display range tau_max", 1.0, 20.0, 8.0, "ns", "A", "lin", 0.5),
        ],
        "hom": [
            _c("gamma", "Decay rate gamma", 0.2, 5.0, 1.0, "1/ns", "A", "lin", 0.05),
            _c("gstar", "Pure dephasing gamma*", 0.0, 50.0, 10.0, "1/ns", "A", "lin", 0.1, "verify (d1): gamma* = 10"),
        ],
        "filter": [
            _c("r", "Pump rate r", 1e-3, 100.0, 0.5, "1/ns", "A", "log"),
            _c("gamma_X", "Radiative rate gamma_X", 0.2, 5.0, 1.0, "1/ns", "E", "lin", 0.05),
            _c("L", "Line FWHM L", 0.1, 100.0, 2.0, "1/ns", "A", "log", note="angular, 1/ns (lindblad.mev_to_rate)"),
            _c("w", "Filter FWHM w", 0.05, 100.0, 2.0, "1/ns", "A", "log"),
            _T_CTL,
        ],
    },
}



# ------------------------------------------------------------------ small helpers
def _spec_of(name, panel=None):
    c = CONTROLS[name]
    if panel is not None:
        c = c[panel]
    return {s["id"]: s for s in c}


def _arg(specs, name):
    s = specs[name]
    raw = request.args.get(name)
    if raw in (None, ""):
        return float(s["default"])
    try:
        v = float(raw)
    except ValueError:
        raise ExploreAError(f"{name} must be a number")
    if not math.isfinite(v):
        raise ExploreAError(f"{name} must be finite")
    if v < s["min"] or v > s["max"]:
        raise ExploreAError(f"{name}={v} is outside the explorer range [{s['min']}, {s['max']}]")
    return v


def _arg_T(specs):
    """T is optional: None means no thermal escape."""
    raw = request.args.get("T")
    if raw in (None, "", "off", "none"):
        return None
    return _arg(specs, "T")


def _fit():
    try:
        d = json.loads(FIT_JSON.read_text(encoding="utf-8"))
    except OSError as exc:
        raise ExploreAError(f"cannot read {FIT_JSON.name}: {exc}", status=500)
    p = d.get("params", {})
    missing = [k for k in FIT_KEYS if p.get(k) is None]
    if missing:
        raise ExploreAError(f"fit bundle out/phase0/fit_params.json lacks {missing}", missing=missing)
    return {k: float(p[k]) for k in FIT_KEYS}


def escape_rates(gamma_X, T):
    """(k_X, k_XX) from the V-a retention fit; (0, 0) when T is None.  The call is
    cw_g2.escape_rates_from_retention (gamma_XX defaults to 2 gamma_X)."""
    if T is None:
        return 0.0, 0.0
    from fsim_core.cw_g2 import escape_rates_from_retention
    f = _fit()
    return escape_rates_from_retention(gamma_X, f["a_esc"], f["E_a"], f["b_p"], f["E_b"], T)


def _thin(arrs, n_max=WIRE_MAX, keep=None):
    """Stride-thin parallel arrays for the wire, always keeping index `keep`."""
    n = len(arrs[0])
    if n <= n_max:
        return arrs, False
    stride = int(math.ceil(n / n_max))
    idx = np.arange(0, n, stride)
    if keep is not None and keep not in idx:
        idx = np.unique(np.concatenate([idx, [keep]]))
    return [a[idx] for a in arrs], True


def _chain(*tags):
    return widest([t for t in tags if t])


def _gates_g2():
    """The 0.5 (g2 ceiling) line comes from campaigns.gates() only."""
    try:
        from . import campaigns
        g = campaigns.gates().get("g2_ceiling") or {}
        return g.get("value"), g.get("tag"), g.get("source")
    except Exception:
        return None, None, None


def _tags(specs, *ids):
    return [specs[i]["tag"] for i in ids]


def _tag_k(T):
    return "DR" if T is not None else None


# ================================================================== 1. cw_g2
def cw_g2_explore(r, gamma_X, eps, rho, irf_fwhm_ps, irf_shape, pump_ratio, T, tau_max):
    from fsim_core import cw_g2

    sp = _spec_of("cw_g2")
    gamma_XX = 2.0 * gamma_X
    k_X, k_XX = escape_rates(gamma_X, T)
    t0 = time.perf_counter()
    rep = cw_g2.cw_report(r, gamma_X, gamma_XX, k_X, k_XX, 1.0, eps, rho, irf_fwhm_ps,
                          tau_max, pump_ratio, irf_shape)
    ms = (time.perf_counter() - t0) * 1e3
    tau = rep["curves"]["tau"]
    n_full = len(tau)
    (tau_w, g2_dot, g2_meas, g2_raw), decimated = _thin(
        [tau, rep["curves"]["g2_dot"], rep["curves"]["g2_meas"], rep["curves"]["g2_raw"]], keep=n_full // 2)

    # Plot B: g2_dot(0) vs r on a log grid, closed form, with its two limits
    r_grid = np.logspace(-3, 2, 101)
    g2_vs_r = cw_g2.g2_cw_zero(r_grid, gamma_X, gamma_XX, eps, k_X, k_XX, pump_ratio)
    low = float(cw_g2.g2_cw_zero_low_pump(gamma_X, gamma_XX, eps, k_X, k_XX, pump_ratio))
    # r -> inf limit gamma_X/(eps gamma_XX) (cw_g2.py:65-77); None for eps = 0 or p = 0
    from fsim_core.studio_support import g2_cw_high_pump_limit
    high = g2_cw_high_pump_limit(gamma_X, gamma_XX, eps, pump_ratio)

    # live consistency rows (comparisons between fsim_core calls, not new physics)
    g2_matrix0 = float(cw_g2.g2_cw(np.array([0.0]), r, gamma_X, gamma_XX, 1.0, eps, k_X, k_XX, pump_ratio)[0])
    M = cw_g2.generator(r, gamma_X, gamma_XX, k_X, k_XX, pump_ratio)
    P = cw_g2.steady_state(r, gamma_X, gamma_XX, k_X, k_XX, pump_ratio)
    balance = float(np.max(np.abs(M @ P)))
    checks = [
        {"id": "closed_vs_matrix", "label": "g2_cw_zero vs g2_cw(0)", "value": abs(rep["g2_dot0"] - g2_matrix0),
         "tol": 1e-10, "tag": "DR"},
        {"id": "detailed_balance", "label": "|generator @ steady_state|", "value": balance, "tol": 1e-12, "tag": "DR"},
    ]
    for c in checks:
        c["ok"] = bool(c["value"] <= c["tol"])

    t_r = _chain(*_tags(sp, "r", "gamma_X", "eps", "pump_ratio"), _tag_k(T))
    # r -> 0 / r -> inf limits: no r, rho, tau_max, IRF; eps [A], p [A], gamma_X [E], k [DR] (review finding 2)
    t_lim = _chain("DR", *_tags(sp, "gamma_X", "eps", "pump_ratio"), _tag_k(T))
    t_rho = _chain(t_r, sp["rho"]["tag"])
    t_irf = _chain(t_rho, "V")
    gate, gate_tag, gate_src = _gates_g2()
    return {
        "inputs": {"r": r, "gamma_X": gamma_X, "gamma_XX": gamma_XX, "eps": eps, "rho": rho,
                   "irf_fwhm_ps": irf_fwhm_ps, "irf_shape": irf_shape, "pump_ratio": pump_ratio,
                   "T": T, "tau_max": tau_max},
        "k_X": k_X, "k_XX": k_XX,
        "readouts": {
            "g2_dot0": rep["g2_dot0"], "g2_meas0": rep["g2_meas0"], "g2_raw0": rep["g2_raw0"],
            "tau_dip": rep["tau_dip"], "A_dip": rep["A_dip"],
            "g2_intrinsic_from_raw": rep["g2_intrinsic_from_raw"],
            "I_X": rep["I_X"], "I_XX": rep["I_XX"], "I_det": rep["I_det"],
        },
        "tags": {"g2_dot0": t_r, "g2_meas0": t_rho, "g2_raw0": t_irf, "tau_dip": t_rho, "A_dip": t_rho,
                 "g2_intrinsic_from_raw": t_irf, "k_X": _chain(sp["gamma_X"]["tag"], "DR") if T is not None else None,
                 "intensity": t_r, "curves": t_irf, "limits": t_lim},
        "bunching": bool(rep["g2_dot0"] > 1.0),
        "curves": {"tau": tau_w, "g2_dot": g2_dot, "g2_meas": g2_meas, "g2_raw": g2_raw},
        "n_tau": n_full, "decimated": decimated,
        "vs_r": {"r": r_grid, "g2_dot0": g2_vs_r, "low_limit": low, "high_limit": high},
        "gate_g2": {"value": gate, "tag": gate_tag, "source": gate_src},
        "checks": checks,
        "caveat": CAVEAT_CW, "caveat_bunching": CAVEAT_CW_BUNCH,
        "notes": rep["notes"],
        "stiff": bool(k_X > 100.0),
        "elapsed_ms": ms,
        "calls": ["fsim_core.cw_g2.escape_rates_from_retention", "fsim_core.cw_g2.cw_report",
                  "fsim_core.cw_g2.g2_cw_zero", "fsim_core.cw_g2.g2_cw_zero_low_pump",
                  "fsim_core.cw_g2.g2_cw", "fsim_core.cw_g2.steady_state",
                  "fsim_core.studio_support.g2_cw_high_pump_limit"],
        "params_source": ("escape rates: out/phase0/fit_params.json V-a fit [DR]" if T is not None
                          else "no thermal escape (k_X = k_XX = 0)"),
    }


# ================================================================== 2. pulse_counting
def pulsed_explore(r, tau_on, period, gamma_X, eps, pump_ratio, T, gate_mode, eta_load):
    from fsim_core import pulse_counting as pc
    from fsim_core.loading import f1b_g2

    sp = _spec_of("pulse_counting")
    if period <= tau_on:
        raise ExploreAError(f"period={period} must exceed tau_on={tau_on}")
    gamma_XX = 2.0 * gamma_X
    k_X, k_XX = escape_rates(gamma_X, T)
    tau_dark = period - tau_on
    gate = None if gate_mode == "none" else tau_on + 5.0 / gamma_X       # [A], pulse_counting.py:118-120
    t0 = time.perf_counter()
    pt = pc.pulse_g2(r, gamma_X, gamma_XX, k_X, k_XX, 1.0, eps, tau_on, tau_dark, pump_ratio,
                     split=True, gate_ns=gate, adjacent=True)
    det = pc.deterministic_cycle_g2(gamma_X, gamma_XX, k_X, k_XX, 1.0, eps, period,
                                    eta_load=eta_load, gate_ns=gate, adjacent=True)
    ms_point = (time.perf_counter() - t0) * 1e3

    # Plot A / B: sweep r (the instantaneous-pulse approximation is f1b_g2(mu = r tau_on, eps))
    r_grid = np.logspace(-2, 4, 49)
    g2s, g2a, mc, conv = [], [], [], []
    for rv in r_grid:
        o = pc.pulse_g2(float(rv), gamma_X, gamma_XX, k_X, k_XX, 1.0, eps, tau_on, tau_dark, pump_ratio,
                        split=False, gate_ns=gate, adjacent=True)
        g2s.append(o["g2"])
        g2a.append(o.get("g2_adj", float("nan")))
        mc.append(o["mean_counts"])
        conv.append(o["converged"])
    f1b = np.asarray(f1b_g2(r_grid * tau_on, eps), dtype=float)
    mu_point = r * tau_on
    f1b_point = float(f1b_g2(mu_point, eps))

    # Plot C: deterministic cycle vs period 1..50 ns
    per_grid = np.linspace(1.0, 50.0, 50)
    d_g2, d_mc, d_blk = [], [], []
    for pv in per_grid:
        o = pc.deterministic_cycle_g2(gamma_X, gamma_XX, k_X, k_XX, 1.0, eps, float(pv),
                                      eta_load=eta_load, gate_ns=gate)
        d_g2.append(o["g2"])
        d_mc.append(o["mean_counts"])
        d_blk.append(o["blocked_load_probability"])
    ms = (time.perf_counter() - t0) * 1e3

    t_dyn = _chain(*_tags(sp, "r", "tau_on", "period", "gamma_X", "eps", "pump_ratio"), _tag_k(T))
    t_det = _chain(*_tags(sp, "period", "gamma_X", "eps", "eta_load"), _tag_k(T))
    # f1b_g2(mu = r tau_on, eps) is a [DR] closed form of r, tau_on and eps; eps and r tau_on are [A]
    t_f1b = _chain("DR", *_tags(sp, "r", "tau_on", "eps"))
    gate_v, gate_tag, gate_src = _gates_g2()
    return {
        "inputs": {"r": r, "tau_on": tau_on, "tau_dark": tau_dark, "period": period, "gamma_X": gamma_X,
                   "gamma_XX": gamma_XX, "eps": eps, "pump_ratio": pump_ratio, "T": T,
                   "gate_mode": gate_mode, "gate_ns": gate, "eta_load": eta_load},
        "k_X": k_X, "k_XX": k_XX, "mu": mu_point,
        "point": {
            "g2": pt["g2"], "g2_adj": pt.get("g2_adj"), "mean_counts": pt["mean_counts"],
            "mean_counts_x": pt.get("mean_counts_x"), "mean_counts_xx": pt.get("mean_counts_xx"),
            "adjacent_peak_factor": pt.get("adjacent_peak_factor"), "converged": pt["converged"],
            "p_period": pt["p_period"], "f1b_g2_inst": f1b_point,
        },
        "deterministic": {
            "g2": det["g2"], "g2_adj": det.get("g2_adj"), "mean_counts": det["mean_counts"],
            "blocked_load_probability": det["blocked_load_probability"],
            "mean_loaded_pairs": det["mean_loaded_pairs"], "one_pair_valid": det["one_pair_valid"],
            "converged": det["converged"], "invalid_reason": det["invalid_reason"],
            "gate_ns_used": det["gate_ns_used"],
        },
        "tags": {"point": t_dyn, "f1b": t_f1b, "deterministic": t_det, "sweep": t_dyn},
        "vs_r": {"r": r_grid, "g2": g2s, "g2_adj": g2a, "mean_counts": mc, "f1b_g2": f1b, "converged": conv},
        "vs_period": {"period": per_grid, "g2": d_g2, "mean_counts": d_mc, "blocked_load_probability": d_blk},
        "gate_g2": {"value": gate_v, "tag": gate_tag, "source": gate_src},
        "caveat": CAVEAT_PULSED, "caveat_extra": CAVEAT_PULSED_EXTRA,
        "elapsed_ms": ms, "point_ms": ms_point,
        "calls": ["fsim_core.cw_g2.escape_rates_from_retention", "fsim_core.pulse_counting.pulse_g2",
                  "fsim_core.pulse_counting.deterministic_cycle_g2", "fsim_core.loading.f1b_g2"],
        "params_source": ("escape rates: out/phase0/fit_params.json V-a fit [DR]" if T is not None
                          else "no thermal escape (k_X = k_XX = 0)"),
    }


# ================================================================== 6. lindblad cross-check
def _phys(rho):
    from fsim_core import lindblad as lb
    tr, mineig, herm = lb.check_physical(rho)
    return {"trace_error": tr, "min_eig": mineig, "hermiticity": herm,
            "ok": bool(tr <= 1e-12 and mineig >= -1e-12 and herm <= 1e-12)}


def lindblad_incoherent(r, gamma_X, eps, pump_ratio, T, tau_max):
    from fsim_core import cw_g2
    from fsim_core import lindblad as lb

    sp = _spec_of("lindblad", "incoherent")
    gamma_XX = 2.0 * gamma_X
    k_X, k_XX = escape_rates(gamma_X, T)
    taus = np.linspace(0.0, tau_max, 201)
    s = lb.build_system(levels=3, r_ns=r, gamma_X_ns=gamma_X, gamma_XX_ns=gamma_XX, k_X=k_X, k_XX=k_XX,
                        pump_ratio=pump_ratio)
    g2L, rho_ss, I = lb.g2_tau(s, taus, weights=[1.0, eps])
    g2C = cw_g2.g2_cw(taus, r, gamma_X, gamma_XX, t_X=1.0, t_XX=eps, k_X=k_X, k_XX=k_XX, pump_ratio=pump_ratio)
    rel = float(np.max(np.abs(g2L - g2C) / np.maximum(np.abs(g2C), 1e-300)))
    t_in = _chain(*_tags(sp, "r", "gamma_X", "eps", "pump_ratio"), _tag_k(T))
    return {
        "panel": "incoherent",
        "inputs": {"r": r, "gamma_X": gamma_X, "gamma_XX": gamma_XX, "eps": eps, "pump_ratio": pump_ratio,
                   "T": T, "tau_max": tau_max},
        "k_X": k_X, "k_XX": k_XX,
        "curves": {"tau": taus, "g2_lindblad": g2L, "g2_cw_g2": g2C},
        "readouts": {"g2_0_lindblad": float(g2L[0]), "g2_0_cw_g2": float(g2C[0]), "max_rel_diff": rel,
                     "detected_rate": float(I)},
        "tags": {"curves": t_in, "readouts": t_in, "max_rel_diff": "DR"},
        "state": _phys(rho_ss),
        "agrees": bool(rel <= 1e-10),
        "tolerance": 1e-10,
    }


def lindblad_rabi(omega, gamma, deph, tau_max):
    from fsim_core import lindblad as lb

    sp = _spec_of("lindblad", "rabi")
    dl = [(deph, (0.0, 1.0))] if deph > 0 else []
    s = lb.build_system(levels=2, gamma_X_ns=gamma, omega_ns=omega, deph=dl)
    taus = np.linspace(0.0, tau_max, 401)
    g2, rho_ss, I = lb.g2_tau(s, taus)
    rho = lb.steady_state(s.L, s.dim)
    i_max = int(np.argmax(g2))
    t_in = _chain(*_tags(sp, "omega", "gamma", "deph"))
    return {
        "panel": "rabi",
        "inputs": {"omega": omega, "gamma": gamma, "deph": deph, "tau_max": tau_max,
                   "gamma_2": gamma / 2.0 + deph / 2.0},
        "curves": {"tau": taus, "g2": g2},
        "readouts": {"rho_ee": float(rho[1, 1].real), "g2_0": float(g2[0]), "g2_max": float(g2[i_max]),
                     "tau_at_max": float(taus[i_max]), "detected_rate": float(I)},
        "tags": {"curves": t_in, "readouts": t_in},
        "state": _phys(rho_ss),
    }


def lindblad_hom(gamma, gstar):
    from fsim_core import lindblad as lb

    sp = _spec_of("lindblad", "hom")
    rho0 = np.diag([0.0, 1.0]).astype(complex)

    def one(gs):
        s = lb.build_system(levels=2, gamma_X_ns=gamma, deph=[(gs, (0.0, 1.0))] if gs > 0 else [])
        return lb.indistinguishability(s.L, s.dim, rho0, s.ops["c_X"])

    gs_grid = np.logspace(-2, 2, 41)
    sweep = [float(one(float(g))) for g in gs_grid]
    t_in = _chain(*_tags(sp, "gamma", "gstar"))
    return {
        "panel": "hom",
        "inputs": {"gamma": gamma, "gstar": gstar},
        "curves": {"gstar": gs_grid, "indistinguishability": sweep},
        "readouts": {"indistinguishability": float(one(gstar))},
        "tags": {"curves": t_in, "readouts": t_in},
    }


def lindblad_filter(r, gamma_X, L, w, T):
    from fsim_core import lindblad as lb

    sp = _spec_of("lindblad", "filter")
    k_X, _ = escape_rates(gamma_X, T)

    deph = lb.dephasing_for_fwhm(L, r_ns=r, gamma_X_ns=gamma_X, k_X=k_X, levels=2)

    def build(wv):
        return lb.build_system(levels=2, r_ns=r, gamma_X_ns=gamma_X, k_X=k_X, deph=deph,
                               filt=dict(fwhm_ns=wv, det_ns=0.0, n_max=4))

    # The similarity-scaled steady-state solve (lindblad.py:430) is deliberately ill-conditioned
    # in the raw basis; scipy's LinAlgWarning is noise, the result is checked by verify_lindblad (e1, e2).
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        s = build(w)
        g2f, rho_f = lb.filtered_g2_zero(s, return_state=True)
        flux = lb.filtered_flux(s)
        w_grid = np.logspace(-1.3, 2, 21)
        sweep = [float(lb.filtered_g2_zero(build(float(wv)))) for wv in w_grid]
    t_in = _chain(*_tags(sp, "r", "gamma_X", "L", "w"), _tag_k(T))
    return {
        "panel": "filter",
        "inputs": {"r": r, "gamma_X": gamma_X, "L": L, "w": w, "T": T, "k_X": k_X, "n_max": 4,
                   "dephasing_rate": float(deph[0][0]) if deph else 0.0,
                   "lifetime_limited": not deph},
        "curves": {"w": w_grid, "g2_filtered_0": sweep},
        "readouts": {"g2_filtered_0": float(g2f), "filtered_flux": float(flux)},
        "tags": {"curves": t_in, "readouts": t_in},
        "state": _phys(rho_f),
    }


# ================================================================== routes
def _out(obj, status=200):
    return Response(dumps(obj), status=status, mimetype="application/json")


def register(app):
    @app.errorhandler(ExploreAError)
    def _explore_a_error(exc):
        return _out({"error": exc.message, **exc.extra}, exc.status)

    @app.get("/api/explore/a/controls")
    def explore_a_controls():
        return _out({
            "controls": CONTROLS,
            "caveats": {"cw_g2": [CAVEAT_CW], "cw_g2_bunching": CAVEAT_CW_BUNCH,
                        "pulse_counting": [CAVEAT_PULSED, CAVEAT_PULSED_EXTRA],
                        "lindblad": [CAVEAT_LINDBLAD, CAVEAT_LINDBLAD_EXTRA]},
        })

    @app.get("/api/explore/cw_g2")
    def explore_cw_g2():
        sp = _spec_of("cw_g2")
        shape = request.args.get("irf_shape", "exponential")
        if shape not in ("gaussian", "exponential"):
            raise ExploreAError("irf_shape must be gaussian or exponential")
        return _out(cw_g2_explore(
            _arg(sp, "r"), _arg(sp, "gamma_X"), _arg(sp, "eps"), _arg(sp, "rho"),
            _arg(sp, "irf_fwhm_ps"), shape, _arg(sp, "pump_ratio"), _arg_T(sp), _arg(sp, "tau_max")))

    @app.get("/api/explore/pulse_counting")
    def explore_pulse_counting():
        sp = _spec_of("pulse_counting")
        gate = request.args.get("gate", "none")
        if gate not in ("none", "auto"):
            raise ExploreAError("gate must be none or auto")
        return _out(pulsed_explore(
            _arg(sp, "r"), _arg(sp, "tau_on"), _arg(sp, "period"), _arg(sp, "gamma_X"), _arg(sp, "eps"),
            _arg(sp, "pump_ratio"), _arg_T(sp), gate, _arg(sp, "eta_load")))

    @app.get("/api/explore/lindblad")
    def explore_lindblad():
        panel = request.args.get("panel", "incoherent")
        if panel not in CONTROLS["lindblad"]:
            raise ExploreAError(f"panel must be one of {sorted(CONTROLS['lindblad'])}")
        sp = _spec_of("lindblad", panel)
        if panel == "incoherent":
            body = lindblad_incoherent(_arg(sp, "r"), _arg(sp, "gamma_X"), _arg(sp, "eps"),
                                       _arg(sp, "pump_ratio"), _arg_T(sp), _arg(sp, "tau_max"))
        elif panel == "rabi":
            body = lindblad_rabi(_arg(sp, "omega"), _arg(sp, "gamma"), _arg(sp, "deph"), _arg(sp, "tau_max"))
        elif panel == "hom":
            body = lindblad_hom(_arg(sp, "gamma"), _arg(sp, "gstar"))
        else:
            body = lindblad_filter(_arg(sp, "r"), _arg(sp, "gamma_X"), _arg(sp, "L"), _arg(sp, "w"), _arg_T(sp))
        body["caveat"] = CAVEAT_LINDBLAD
        body["caveat_extra"] = CAVEAT_LINDBLAD_EXTRA
        body["calls"] = ["fsim_core.lindblad.build_system", "fsim_core.lindblad.g2_tau",
                         "fsim_core.lindblad.steady_state", "fsim_core.lindblad.indistinguishability",
                         "fsim_core.lindblad.filtered_g2_zero", "fsim_core.lindblad.filtered_flux",
                         "fsim_core.lindblad.check_physical"]
        return _out(body)

    return app


__all__ = ["register", "cw_g2_explore", "pulsed_explore", "lindblad_incoherent", "lindblad_rabi",
           "lindblad_hom", "lindblad_filter", "ExploreAError", "CONTROLS"]
