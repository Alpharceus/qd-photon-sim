"""Physics helpers behind the FSIM Studio 3D scenes and explainers.

Two kinds of helper live here.

Moved from the Flask layer (fsim_studio/scene.py, fsim_studio/api_explain.py)
because the three-layer rule puts every physical construction in fsim_core:
per_s_to_per_ns, pulse_timing_ns, cascade_trajectory, dbr_defaults,
planar_cavity_optics, band_tilt_eV.  verify/verify_scene_support.py proves
the Studio SceneSpecs they feed are identical to the pre-move output.

New in the move (not a relocation of older code): background_g2_floor,
is_background_floor and explain_fit_params (with MissingFitParams).  They
replace logic the Flask layer used to inline -- the b_res floor test and the
explainer's silent fit-parameter fallbacks (a fitted mu of 0 used to become
0.33; that behaviour change is documented in the verifier) -- and
verify/verify_scene_support.py checks them against independent values.

Nothing here is new physics: each helper composes existing fsim_core calls
(lindblad, nitride_cavity.dbr_diagnostic,
nitride_nanowire_photonics.gan_ordinary_index, integrator.g2_from,
loading.gamma_eff) or performs a unit conversion.
"""
from __future__ import annotations

import inspect
import math

import numpy as np

from . import lindblad, nitride_cavity
from .integrator import g2_from
from .loading import gamma_eff
from .nitride_nanowire_photonics import gan_ordinary_index

# ------------------------------------------------------------------ units

NS_PER_S = 1e-9      # 1/s -> 1/ns (unit conversion only)
KVCM_TO_V_PER_NM = 1e-4  # 1 kV/cm = 1e3 V / 1e7 nm (unit conversion only)


def per_s_to_per_ns(rate_per_s):
    """Rate in 1/s -> 1/ns."""
    return rate_per_s * NS_PER_S


def pulse_timing_ns(rep_rate_hz, *, duty=None, tau_pulse_ns=None):
    """(period_ns, tau_on_ns, tau_dark_ns) for a pulsed drive.  The on-time
    is duty * period (planar cards) or the stated pulse width (nanowire sweep
    rows); exactly one of duty / tau_pulse_ns must be given."""
    if (duty is None) == (tau_pulse_ns is None):
        raise ValueError("give exactly one of duty / tau_pulse_ns")
    period = 1e9 / rep_rate_hz
    tau_on = float(duty) * period if duty is not None else tau_pulse_ns
    return period, tau_on, period - tau_on


# --------------------------------------------------------------- cascade

def cascade_trajectory(*, r_ns, gamma_X_ns, gamma_XX_ns, k_X, k_XX, tau_on_ns, tau_dark_ns,
                       pump_ratio, deterministic_pair, n_on=121, n_off=280, settle_periods=3):
    """XX -> X -> 0 populations over one drive period from the 3-level
    Lindblad ladder (lindblad.build_system / evolve).

    Drive construction (as the Studio cascade view always used it):
    pump r_ns during the on-window, r = 0 in the dark window, the card's
    drive.cw_pump_ratio reused as the XX/X pump ratio for the pulsed drive
    [A] (no separate pulsed ratio exists on the card).  The start state is
    |XX> for deterministic_pair loading (the deterministic load map resets
    each period to |XX>), otherwise the ground state settled for
    `settle_periods` periods toward the periodic state [A: 3 periods].

    Returns {t_ns, P_G, P_X, P_XX, rate_X, rate_XX, rate_esc (1/ns,
    gamma*P and k*P), per_pulse {X, XX, esc} (time integrals of the rates
    over the period), tau_on_ns, period_ns, start_state}.
    """
    on = lindblad.build_system(levels=3, r_ns=float(r_ns), gamma_X_ns=float(gamma_X_ns),
                               gamma_XX_ns=float(gamma_XX_ns), k_X=float(k_X), k_XX=float(k_XX),
                               pump_ratio=pump_ratio)
    off = lindblad.build_system(levels=3, r_ns=0.0, gamma_X_ns=float(gamma_X_ns),
                                gamma_XX_ns=float(gamma_XX_ns), k_X=float(k_X), k_XX=float(k_XX),
                                pump_ratio=pump_ratio)
    det_pair = bool(deterministic_pair)
    rho = np.zeros((3, 3), dtype=complex)
    rho[2 if det_pair else 0, 2 if det_pair else 0] = 1.0
    t_on = np.linspace(0.0, float(tau_on_ns), n_on)
    t_off = np.concatenate([[0.0], np.geomspace(1e-5, float(tau_dark_ns), n_off - 1)])
    ops = on.ops
    for _ in range(settle_periods):  # settle to the periodic state
        r_on = lindblad.evolve(on.L, rho, t_on)
        r_off = lindblad.evolve(off.L, r_on[-1], t_off)
        rho = r_off[-1]
        if det_pair:  # deterministic load map: the next period starts from |XX>
            rho = np.zeros((3, 3), dtype=complex)
            rho[2, 2] = 1.0
    traj = r_on + r_off[1:]
    t = np.concatenate([t_on, float(tau_on_ns) + t_off[1:]])
    P = {k: np.array([lindblad.expect(ops[k], x).real for x in traj]) for k in ("P_G", "P_X", "P_XX")}
    rate_X = float(gamma_X_ns) * P["P_X"]
    rate_XX = float(gamma_XX_ns) * P["P_XX"]
    rate_esc = float(k_X) * P["P_X"] + float(k_XX) * P["P_XX"]
    per_pulse = {k: float(np.trapezoid(rr, t))
                 for k, rr in (("X", rate_X), ("XX", rate_XX), ("esc", rate_esc))}
    return {"t_ns": t, "P_G": P["P_G"], "P_X": P["P_X"], "P_XX": P["P_XX"],
            "rate_X": rate_X, "rate_XX": rate_XX, "rate_esc": rate_esc, "per_pulse": per_pulse,
            "tau_on_ns": float(tau_on_ns), "period_ns": float(tau_on_ns) + float(tau_dark_ns),
            "start_state": "|XX>" if det_pair else f"periodic ({settle_periods} periods settled)"}


# ---------------------------------------------------------- planar cavity

def dbr_defaults() -> dict:
    """n_high, n_low, pairs as nitride_cavity.dbr_diagnostic defines them
    (its keyword defaults are the single source; exploratory [E/A])."""
    sig = inspect.signature(nitride_cavity.dbr_diagnostic).parameters
    return {"n_high": sig["n_high"].default, "n_low": sig["n_low"].default, "pairs": sig["pairs"].default}


def planar_cavity_optics(lambda_nm, mode_volume_norm) -> dict:
    """Optical glyph inputs for the planar nitride cavity view.

    n_GaN = gan_ordinary_index(lambda); mode volume V = V_norm (lambda/n)^3
    (the [A] normalized-volume convention of nitride_cavity, no 3-D
    geometry), the equivalent-sphere radius, the lambda/2n standing-wave
    period, and the dbr_diagnostic stack at its own default indices / pairs.
    """
    lam = float(lambda_nm)
    n_gan = gan_ordinary_index(lam)
    V_um3 = mode_volume_norm * (lam / 1000.0 / n_gan) ** 3
    dd = dbr_defaults()
    dbr = nitride_cavity.dbr_diagnostic(lam, n_high=dd["n_high"], n_low=dd["n_low"], pairs=dd["pairs"])
    return {
        "n_gan": n_gan, "V_um3": V_um3,
        "radius_um": (3.0 * V_um3 / (4.0 * math.pi)) ** (1.0 / 3.0),
        "standing_wave_period_nm": lam / (2.0 * n_gan),
        "dbr": {**dbr, **dd, "total_nm": dd["pairs"] * (dbr["d_high_nm"] + dbr["d_low_nm"])},
    }


# ----------------------------------------------------------------- bands

def band_tilt_eV(field_kVcm, height_nm):
    """Electrostatic drop |F| * h across the dot (eV): kV/cm -> V/nm."""
    return abs(field_kVcm) * KVCM_TO_V_PER_NM * height_nm


# ------------------------------------------------------- background floor

def background_g2_floor(b_res):
    """g2(0) of an ideal emitter (g2_dot = 0) seen through a residual
    background b_res per collected X photon: rho = 1/(1 + b_res),
    g2 = integrator.g2_from(0, rho) = 1 - rho^2.  This is the structural
    floor a deterministic (pair) source cannot go below at that b_res."""
    rho = 1.0 / (1.0 + float(b_res))
    return float(g2_from(0.0, rho))


def is_background_floor(g2_op, b_res, tol=1e-9) -> bool:
    """True when g2_op equals the b_res background floor within tol."""
    try:
        g, b = float(g2_op), float(b_res)
    except (TypeError, ValueError):
        return False
    if not (math.isfinite(g) and math.isfinite(b)) or b <= 0.0:
        return False
    return abs(g - background_g2_floor(b)) <= tol


# ------------------------------------------------------ explainer params

EXPLAIN_SPECTRAL_KEYS = ("delta_xx", "gamma0", "a_ac", "b_lo", "E_lo", "mu")
EXPLAIN_CASCADE_KEYS = ("delta_xx", "gamma0", "a_ac", "b_lo", "E_lo", "mu", "a_esc", "E_a", "b_p",
                        "E_b", "b0", "beta")


class MissingFitParams(KeyError):
    def __init__(self, missing):
        super().__init__(f"fit bundle lacks {missing}")
        self.missing = list(missing)


def explain_fit_params(params: dict, required) -> dict:
    """Resolve the explainer's linewidth / loading parameters from a fit
    bundle without any silent physics fallback: every `required` key must be
    present (else MissingFitParams lists them).  p_inj, which the V-a bundle
    does not fit, resolves to fsim_core.loading.gamma_eff's own default
    exponent [A], and the returned sources say so.

    Returns {"params": {...}, "sources": {key: text}}.
    """
    missing = [k for k in required if params.get(k) is None]
    if missing:
        raise MissingFitParams(missing)
    out = {k: params[k] for k in required}
    src = {k: "fit bundle" for k in required}
    if params.get("p_inj") is not None:
        out["p_inj"] = float(params["p_inj"])
        src["p_inj"] = "fit bundle"
    else:
        out["p_inj"] = float(inspect.signature(gamma_eff).parameters["p_inj"].default)
        src["p_inj"] = "fsim_core.loading.gamma_eff default p_inj [A] (not in the fit bundle)"
    return {"params": out, "sources": src}


__all__ = ["per_s_to_per_ns", "pulse_timing_ns", "cascade_trajectory", "dbr_defaults",
           "planar_cavity_optics", "band_tilt_eV", "background_g2_floor", "is_background_floor",
           "explain_fit_params", "MissingFitParams", "EXPLAIN_SPECTRAL_KEYS", "EXPLAIN_CASCADE_KEYS"]
