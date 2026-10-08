"""Model explorers, set B (studio-p2d, Coder E2): phonon spectra (qd_gf), transport + nitride
Stark, and the Langevin QD-laser (sde, Zhao).

Three-layer rule: nothing here is physics. Every curve and readout is a public fsim_core call
(qd_gf.huang_rhys / zpl_weight / ibm_spectrum / ibm_transmission / ibm_purcell_transmission /
effective_gamma_zpl, spectral.transmission, transport.evaluate_injection / Diode.depletion,
nitride_transport.evaluate_injection, nitride_stark.resolve_bias / stark_derivatives,
nitride_levels.levels, sde.steady_state / find_threshold / g2_vs_pump / simulate).

Cost rule (brief 3(e), 4(e), 5(e)): live endpoints only where the brief says live; everything else
is a job on the pool through JobManager.submit_call with exact allowlist entries (jobs.py):
  fsim_studio.api_explore_b:phonon_frames   one batch of T frames per task (the qd_gf lru caches
                                            are per process, so a batch keeps one worker warm)
  fsim_studio.api_explore_b:transport_sweep the 50-point I sweep (one brentq per point)
  fsim_studio.api_explore_b:stark_trace     one nitride Stark trace (levels 144 ms cold per point)
  fsim_studio.api_explore_b:sde_frame       one seeded g2_vs_pump frame per (I/I_th, F_pump)
  fsim_studio.api_explore_b:sde_li          steady-state L-I curve plus find_threshold
  fsim_studio.api_explore_b:sde_trace       one simulate() S(t) trace

Routes (registered by register(app) from server.create_app):
  GET  /api/explore/phonon/meta | /zt | /live        POST /api/explore/phonon/frames
  GET  /api/explore/transport/meta | /point | /depletion | /nitride_point
  POST /api/explore/transport/sweep | /stark
  GET  /api/explore/sde/meta | /card                 POST /api/explore/sde/frames | /li | /trace

The Zhao/sde endpoints refuse disable_stim, beta_sp and every tuned Tier-2 parameter (brief 5(b)
and 5(c): the default-beta estimator control gave g2 = -1.85 +- 0.27 and must never be shown).
"""
from __future__ import annotations

import csv
import math

import numpy as np

from . import CARDS, OUT
from .tags import widest

# ----------------------------------------------------------------- pinned caveats (verbatim from the brief)
PHONON_CAVEAT = (
    "Independent-boson model, bulk LA deformation potential only, no LO or piezo term (:191-197) [E]. "
    "Gamma_zpl is phenomenological; the fitted a_ac T term partly double-counts acoustic broadening "
    "(:436-443). Lengths are wavefunction extents, not dot sizes (:127-134).")
PHONON_NOT_BRIGHTNESS = ("The ZPL-weight panel is not brightness (retention and collection missing, "
                         "README \"Model scope\").")
TRANSPORT_CAVEAT = (
    "Compact p-i-n [E/A]: Boltzmann statistics, abrupt depletion, SRH n=2; leak_valleys='G' is a LOWER "
    "bound on leakage (docstring item 5, 'GX' moves 50x per 0.1 eV); b_e excludes neighbour dots (F5) "
    "and the X filter (transport.py:1-30). 1 uA, 1 ns, one dot gives mu ~ 6000; cap-2 needs ~80 pA into "
    "one dot (item 6).")
NITRIDE_CAVEAT = (
    "barrier-limited, doping dependence omitted [A]; screening fixed along a trace; m- and a-plane "
    "coincide by construction (out/nitride_geometry_stark/results.md:328); Zhang is 10 K, x=0.15, "
    "another device (results.md:129)")
SDE_CAVEAT = (
    "Langevin model of Zhao's QD LASER (thousands of dots, lasing mode), not a single-photon source; "
    "g2 ~ 1 is correct physics. Normal pump reproduces 1.0224+-0.003 to 1.7 sigma; quiet pump does NOT: "
    "1.017 vs 0.982, 5.8 sigma (out/zhao/zhao_fit_comparison.csv, verdict FAIL both rows). Mechanism "
    "and 12 tuning experiments: sde.py:59-91; Supplement parameters unavailable (sde.py:15-23).")
SDE_AXIS_NOTE = "Axis is I/I_th, not amperes (sde.py:196-203). No curve is drawn through the quiet point."


class ExploreError(Exception):
    def __init__(self, message, status=400, **extra):
        super().__init__(message)
        self.message = message
        self.status = status
        self.extra = extra


# ----------------------------------------------------------------- argument handling
def _num(name, v, lo, hi, default=None, allow_none=False):
    if v in (None, ""):
        if allow_none:
            return None
        if default is None:
            raise ExploreError(f"{name} is required")
        v = default
    try:
        f = float(v)
    except (TypeError, ValueError):
        raise ExploreError(f"{name} must be a number")
    if not math.isfinite(f):
        raise ExploreError(f"{name} must be finite")
    if f < lo or f > hi:
        raise ExploreError(f"{name}={f:g} is outside the explorer range [{lo:g}, {hi:g}]")
    return f


def _num_list(name, v, default, max_len):
    """A JSON array of numbers (or `default` when absent/empty); anything else is a 400, never a 500."""
    if v in (None, [], ""):
        return list(default)
    if not isinstance(v, (list, tuple)):
        raise ExploreError(f"{name} must be an array of numbers")
    if len(v) > max_len:
        raise ExploreError(f"{name}: at most {max_len} values")
    return list(v)


def _reject(body: dict, forbidden: tuple):
    bad = [k for k in forbidden if k in body]
    if bad:
        raise ExploreError(f"{bad} are not exposed by this explorer (tuned Tier-2 parameters and the "
                           "stimulated-emission-off estimator control are verify-only)")


# =================================================================== phonon (qd_gf)
PHONON_DOMAIN = {"T": (1.0, 300.0), "l_xy": (2.0, 8.0), "l_z": (1.0, 4.0), "alpha": (0.005, 0.1),
                 "gamma_zpl": (0.05, 10.0), "w": (0.5, 20.0), "kappa": (0.3, 3.0),
                 "F_cav": (1.0, 20.0), "delta": (-12.0, 12.0)}
PHONON_DEFAULTS = {"T": 100.0, "l_xy": 4.5, "l_z": 1.5, "alpha": None, "gamma_zpl": 0.5, "w": 2.0,
                   "kappa": None, "kappa_on": 1.0, "F_cav": 1.0, "delta": 0.0}
# V-a fit (out/phase0/fit_params.json [DR]) feeds effective_gamma_zpl
FIT_JSON = OUT / "phase0" / "fit_params.json"
PHONON_T_FRAMES = (4.0, 20.0, 40.0, 60.0, 80.0, 100.0, 130.0, 160.0, 200.0, 230.0, 260.0, 300.0)
OMEGA_GRID = np.linspace(-25.0, 25.0, 5001)   # photon energy relative to the ZPL, blue positive
DELTA_GRID = np.linspace(-12.0, 12.0, 49)     # line offset from the filter centre


def _params(l_xy, l_z, alpha):
    from fsim_core.qd_gf import PhononParams
    return PhononParams(l_xy_nm=float(l_xy), l_z_nm=float(l_z),
                        alpha_ps2=None if alpha is None else float(alpha))


def phonon_state(l_xy, l_z, alpha, T, gamma_zpl, w, kappa, F_cav, delta):
    """One (geometry, T) state: the IBM spectrum, the transmission curves against the
    Lorentzian-only gap, and the Purcell-reweighted point. Used by the job and the live route."""
    from fsim_core import qd_gf, spectral
    P = _params(l_xy, l_z, alpha)
    S_total = qd_gf.huang_rhys(P, T)
    Z = qd_gf.zpl_weight(P, T)
    spec = qd_gf.ibm_spectrum(OMEGA_GRID, P, T, gamma_zpl)
    t_ibm = qd_gf.ibm_transmission(DELTA_GRID, P, T, gamma_zpl, w_meV=w, kappa_meV=kappa)
    t_lor = np.array([float(spectral.transmission(float(d), gamma_zpl, w=w, kappa=kappa)) for d in DELTA_GRID])
    t_pt = float(qd_gf.ibm_transmission(float(delta), P, T, gamma_zpl, w_meV=w, kappa_meV=kappa))
    t_lor_pt = float(spectral.transmission(float(delta), gamma_zpl, w=w, kappa=kappa))
    tp, z_eff, rate_mult = qd_gf.ibm_purcell_transmission(float(delta), P, T, gamma_zpl, float(F_cav),
                                                          w_meV=w, kappa_meV=kappa)
    return {
        "T": float(T), "S_total": float(S_total), "Z": float(Z), "sideband_fraction": float(1.0 - Z),
        "gamma_zpl": float(gamma_zpl), "w": float(w), "kappa": None if kappa is None else float(kappa),
        "F_cav": float(F_cav), "delta": float(delta),
        "spectrum": {"omega": OMEGA_GRID.tolist(), "S": spec.tolist()},
        "transmission": {"delta": DELTA_GRID.tolist(), "ibm": np.asarray(t_ibm).tolist(), "lorentz": t_lor.tolist()},
        "point": {"t_ibm": t_pt, "t_lorentz": t_lor_pt, "optimism_gap": t_lor_pt - t_pt, "t_purcell": tp,
                  "Z_eff": z_eff, "rate_mult": rate_mult},
        "gamma_fit": phonon_gamma_fit(float(T)), "alpha_ps2": float(qd_gf.coupling_alpha_ps2(P)),
    }


def phonon_frames(l_xy, l_z, alpha, T_list, gamma_zpl, w, kappa, F_cav, delta):
    """Job body: one batch of T frames in ONE worker so the per-process qd_gf caches stay warm."""
    frames = [phonon_state(l_xy, l_z, alpha, float(T), gamma_zpl, w, kappa, F_cav, delta) for T in T_list]
    return {"frames": frames, "T": [float(T) for T in T_list],
            "inputs": {"l_xy": l_xy, "l_z": l_z, "alpha": alpha, "gamma_zpl": gamma_zpl, "w": w,
                       "kappa": kappa, "F_cav": F_cav, "delta": delta}}


def _phonon_args(src):
    """Validate geometry / filter arguments from a query string or JSON body."""
    g = src.get
    kappa_on = str(g("kappa_on", "1")) not in ("0", "false", "False")
    alpha = g("alpha")
    return {
        "l_xy": _num("l_xy", g("l_xy"), *PHONON_DOMAIN["l_xy"], default=PHONON_DEFAULTS["l_xy"]),
        "l_z": _num("l_z", g("l_z"), *PHONON_DOMAIN["l_z"], default=PHONON_DEFAULTS["l_z"]),
        "alpha": _num("alpha", alpha, *PHONON_DOMAIN["alpha"], allow_none=True),
        "gamma_zpl": _num("gamma_zpl", g("gamma_zpl"), *PHONON_DOMAIN["gamma_zpl"], default=PHONON_DEFAULTS["gamma_zpl"]),
        "w": _num("w", g("w"), *PHONON_DOMAIN["w"], default=PHONON_DEFAULTS["w"]),
        "kappa": _num("kappa", g("kappa"), *PHONON_DOMAIN["kappa"], default=1.0) if kappa_on else None,
        "F_cav": _num("F_cav", g("F_cav"), *PHONON_DOMAIN["F_cav"], default=1.0),
        "delta": _num("delta", g("delta"), *PHONON_DOMAIN["delta"], default=0.0),
    }


def phonon_meta():
    from fsim_core import qd_gf
    P = qd_gf.PhononParams()
    out = _phonon_meta(qd_gf, P)
    tg = out["tags"]
    tg["result"] = widest([tg[k] for k in ("geometry", "alpha", "material", "model", "filter", "cavity")])
    return out


def _phonon_meta(qd_gf, P):
    return {
        "domain": {k: list(v) for k, v in PHONON_DOMAIN.items()}, "defaults": PHONON_DEFAULTS,
        "T_frames": list(PHONON_T_FRAMES), "alpha_default_ps2": float(qd_gf.coupling_alpha_ps2(P)),
        "alpha_ramsay_ps2": 0.027,
        "caveat": PHONON_CAVEAT, "not_brightness": PHONON_NOT_BRIGHTNESS,
        "tags": {"geometry": "A", "alpha": "DR", "material": "DR", "model": "E", "gamma_fit": "DR",
                 "filter": "A", "cavity": "A", "result": None},
        "tag_notes": {"geometry": "l_xy 4.5, l_z 1.5 nm: [DR] class-proxy, [A] until AFM/TEM",
                      "alpha": "0.0181 ps^2 from InP-class constants [DR]; Ramsay 0.027 [DR] (cards/chatzarakis2023.yaml:21)",
                      "filter": "w, kappa: qcap-cavity.yaml:35-39 [A]",
                      "result": "widest tag of the inputs (geometry [A])"},
        "calls": ["fsim_core.qd_gf.huang_rhys", "fsim_core.qd_gf.zpl_weight", "fsim_core.qd_gf.ibm_spectrum",
                  "fsim_core.qd_gf.ibm_transmission", "fsim_core.qd_gf.ibm_purcell_transmission",
                  "fsim_core.spectral.transmission", "fsim_core.qd_gf.effective_gamma_zpl"],
    }


def phonon_zt(l_xy, l_z, alpha):
    """Z(T) and S_total(T) on 1..300 K (about 2 ms per point: live)."""
    from fsim_core import qd_gf
    P = _params(l_xy, l_z, alpha)
    Ts = np.linspace(1.0, 300.0, 60)
    S = [qd_gf.huang_rhys(P, float(T)) for T in Ts]
    Z = [qd_gf.zpl_weight(P, float(T)) for T in Ts]
    return {"T": Ts.tolist(), "S_total": S, "Z": Z, "alpha_ps2": float(qd_gf.coupling_alpha_ps2(P)),
            "alpha_override": alpha is not None}


def phonon_gamma_fit(T):
    import json
    from fsim_core import qd_gf
    p = json.loads(FIT_JSON.read_text(encoding="utf-8"))["params"]
    need = ("gamma0", "a_ac", "b_lo", "E_lo")
    miss = [k for k in need if k not in p]
    if miss:
        raise ExploreError(f"fit bundle {FIT_JSON.name} lacks {miss}")
    return float(qd_gf.effective_gamma_zpl(T, p["gamma0"], p["a_ac"], p["b_lo"], p["E_lo"]))


# =================================================================== transport (InP + nitride)
TRANSPORT_DOMAIN = {"I_uA": (1e-3, 200.0), "T": (230.0, 300.0), "n_dot": (1e8, 2e10), "aperture": (0.1, 3.0),
                    "tau_pulse": (0.01, 2.0), "w": (1.0, 3.0), "dE_WL": (100.0, 100.0), "tau_rad": (0.2, 5.0),
                    "V_j": (0.0, 3.2), "E_X_eV": (1.0, 4.0), "ext": (-500.0, 500.0)}
TRANSPORT_DEFAULTS = {"preset": "red", "I_uA": 1.0, "T": 230.0, "n_dot": 2e9, "aperture": 1.0,
                      "tau_pulse": 0.1, "w": 2.0, "dE_WL": 100.0, "tau_rad": 1.0, "T_nitride": 300.0,
                      "orientation": "c_plane", "polarity": 1, "ext": 0.0}
PRESETS = ("red", "hkust", "gaas")
STARK_V = tuple(round(0.2 * i, 4) for i in range(17))   # 0..3.2 V, one point per 0.2 V (brief section 7 rule 6)
SCREENINGS = (0.0, 0.5, 1.0)


def _diode(preset):
    from fsim_core import transport
    if preset == "red":
        return transport.red_diode_preset()
    if preset == "hkust":
        return transport.hkust_preset()
    if preset == "gaas":
        return transport.gaas_homojunction()
    raise ExploreError(f"preset must be one of {PRESETS}")


def _inj_point(diode, I, T, n_dot, aperture, tau_pulse, w, dE_WL, tau_rad, E_X_eV=None):
    from fsim_core import transport
    r = transport.evaluate_injection(diode, float(I), float(T), float(n_dot), float(aperture), float(tau_pulse),
                                     float(w), float(dE_WL), tau_rad_ns=float(tau_rad), E_X_eV=E_X_eV)
    return r


def _inj_row(r, tau_rad):
    d, ld, lk, bg = r.depletion, r.loading, r.leakage, r.background
    return {
        "I_uA": r.I_uA, "T": r.T, "V_j": r.V_j, "V_applied": r.V_applied, "V_bi": r.V_bi, "V_on": r.V_on,
        "eta_inj": lk.eta_inj, "leak_valleys": lk.valleys, "mu": r.mu, "b_e": r.b_e,
        "P_junction_W": r.P_junction_W, "r_dot": ld.r_dot, "r_dot_max": 1.0 / (tau_rad * 1e-9),
        "r_captured": ld.r_captured, "r_matrix": ld.r_matrix, "f_qfl": r.f_qfl, "saturated": bool(ld.saturated),
        "F_kVcm": d.F_kVcm, "C_dep_pF": d.C_dep_pF, "W_nm": d.W_nm, "flat_band": bool(d.flat_band),
        "f_QD": ld.f_QD, "xi": bg.xi, "E_U_meV": bg.E_U_meV, "rate_bg_window": bg.rate_bg_window, "rate_x": bg.rate_x,
    }


def transport_point(preset, I, T, n_dot, aperture, tau_pulse, w, dE_WL, tau_rad, E_X_eV=None):
    diode = _diode(preset)
    r = _inj_point(diode, I, T, n_dot, aperture, tau_pulse, w, dE_WL, tau_rad, E_X_eV)
    return _inj_row(r, tau_rad)


def transport_sweep(preset, T, n_dot, aperture, tau_pulse, w, dE_WL, tau_rad, I_lo, I_hi, n):
    """Job body: the I sweep (one brentq per point, about 18 ms each)."""
    from fsim_core import transport
    diode = _diode(preset)
    Is = np.geomspace(I_lo, I_hi, int(n))
    rows = [_inj_row(_inj_point(diode, float(I), T, n_dot, aperture, tau_pulse, w, dE_WL, tau_rad), tau_rad) for I in Is]
    keys = ("I_uA", "V_j", "V_applied", "eta_inj", "mu", "b_e", "P_junction_W", "r_dot", "r_captured", "r_matrix", "saturated")
    # the legacy-form fit b = A (I/I_ref)^m exp(-E_act/kT), shown WITH its error (never silently)
    _, info = transport.to_background_channel(diode, np.geomspace(0.01, 10.0, 8), [230.0, 250.0, 273.0, 300.0],
                                              w_meV=w, dE_WL_meV=dE_WL, n_dot_cm2=n_dot, aperture_um2=aperture,
                                              tau_rad_ns=tau_rad)
    return {"preset": preset, "T": T, "n": int(n), "tau_rad": tau_rad, "r_dot_max": 1.0 / (tau_rad * 1e-9),
            "columns": {k: [row[k] for row in rows] for k in keys},
            "bg_fit": {"max_rel_err": float(info["max_rel_err"]), **{k: float(v) for k, v in info["drive_fields"].items()},
                       "I_range_uA": [0.01, 10.0], "T_grid": [230.0, 250.0, 273.0, 300.0]}}


def transport_depletion(preset, T, n=60):
    """F(V_j), C_dep(V_j), W(V_j) from V_j = 0 to V_bi (flat band beyond), via Diode.depletion."""
    diode = _diode(preset)
    V_bi = float(diode.vbi(T))
    Vs = np.concatenate([np.linspace(0.0, V_bi, int(n)), [V_bi + 0.1]])
    deps = [diode.depletion(float(v), T) for v in Vs]
    return {"preset": preset, "T": T, "V_bi": V_bi, "V_j": Vs.tolist(),
            "F_kVcm": [d.F_kVcm for d in deps], "C_dep_pF": [d.C_dep_pF for d in deps],
            "W_nm": [d.W_nm for d in deps], "flat_band": [bool(d.flat_band) for d in deps]}


def nitride_point(I_uA, T, n_dot, aperture, tau_pulse, w, dE_WL, tau_rad, polarity, ext):
    """Nitride injection point: nitride_transport.evaluate_injection plus nitride_stark.resolve_bias."""
    from fsim_core import nitride_stark, nitride_transport
    diode = nitride_transport.NitrideDiode()
    r = nitride_transport.evaluate_injection(diode, float(I_uA), float(T), float(n_dot), float(aperture),
                                             float(tau_pulse), float(w), float(dE_WL), tau_rad_ns=float(tau_rad))
    b = nitride_stark.resolve_bias(diode, T_j_K=float(T), current_uA=float(I_uA),
                                   field_polarity=int(polarity), external_field_kVcm=float(ext))
    row = _inj_row(r, tau_rad)
    row["bias"] = {k: b[k] for k in ("V_j", "V_terminal", "current_uA", "diode_field_kVcm", "applied_field_kVcm",
                                     "flat_band", "depletion_regime", "bias_valid", "invalid_reasons")}
    return row


def stark_trace(orientation, screening_fraction, polarity, T, ext, V_list):
    """Job body: one nitride Stark trace. Each row is resolve_bias -> levels; the slopes come
    from stark_derivatives (never across invalid rows or regime kinks)."""
    from fsim_core import nitride_levels, nitride_stark, nitride_transport
    diode = nitride_transport.NitrideDiode()
    rows = []
    for i, V in enumerate(V_list):
        b = nitride_stark.resolve_bias(diode, T_j_K=float(T), junction_voltage_V=float(V),
                                       field_polarity=int(polarity), external_field_kVcm=float(ext))
        row = {"row_id": i, "V_j": float(V), "E_X_eV": float("nan"), "field_kVcm": float("nan"),
               "spectroscopy_valid": False, "depletion_regime": b["depletion_regime"], "flat_band": b["flat_band"],
               "applied_field_kVcm": b["applied_field_kVcm"], "diode_field_kVcm": b["diode_field_kVcm"],
               "current_uA": b["current_uA"], "bias_valid": b["bias_valid"]}
        if b["bias_valid"]:
            sysm = nitride_levels.NitrideDotSystem(orientation=orientation, screening_fraction=float(screening_fraction),
                                                   external_field_kVcm=float(b["applied_field_kVcm"]))
            lv = nitride_levels.levels(sysm, float(T))
            row.update(E_X_eV=float(lv.E_X_eV), field_kVcm=float(lv.field_kVcm), overlap_sq=float(lv.overlap_sq),
                       spectroscopy_valid=bool(lv.valid and b["bias_valid"]))
        rows.append(row)
    rows = nitride_stark.stark_derivatives(rows)
    cols = ("row_id", "V_j", "E_X_eV", "field_kVcm", "overlap_sq", "spectroscopy_valid", "depletion_regime", "flat_band",
            "applied_field_kVcm", "current_uA", "dE_X_dV_meV_per_V", "derivative_valid")
    return {"orientation": orientation, "screening_fraction": float(screening_fraction), "polarity": int(polarity),
            "T": float(T), "ext": float(ext),
            "rows": [{k: r.get(k) for k in cols} for r in rows],
            "zhang_slope_meV_per_V": nitride_stark.ZHANG2016_SLOPE_MEV_PER_V}


# Tags of the transport readouts (review finding 2: the widest-tag rule). An output carries the widest
# tag of its own model line and of every control it depends on; the controls are the tags of
# transport.js (I, T, n_dot, aperture, tau_pulse, tau_rad [A]; w, dE_WL [E]).
TRANSPORT_INPUT_TAGS = {"I": "A", "T": "A", "n_dot": "A", "aperture": "A", "tau_pulse": "A", "tau_rad": "A",
                        "w": "E", "dE_WL": "E"}
# output: (tag of the model line, controls it depends on)
TRANSPORT_OUTPUT_DEPS = {
    "V_j": ("DR", ("I", "T")), "V_bi": ("DR", ("T",)), "F": ("DR", ("I", "T")), "C_dep": ("DR", ("I", "T")),
    "eta_inj": ("E", ("I", "T")), "mu": ("E", ("I", "T", "n_dot", "aperture", "tau_pulse", "tau_rad", "w", "dE_WL")),
    "b_e": ("E", ("I", "T", "n_dot", "aperture", "tau_pulse", "tau_rad", "w", "dE_WL")),
    "P_junction": ("E", ("I", "T")), "r_dot": ("E", ("I", "T", "n_dot", "aperture", "tau_rad")),
    "xi": ("E", ("T", "w", "dE_WL")),
}


def transport_output_tags():
    return {k: widest([model] + [TRANSPORT_INPUT_TAGS[d] for d in deps])
            for k, (model, deps) in TRANSPORT_OUTPUT_DEPS.items()}


def transport_meta():
    return {
        "domain": {k: list(v) for k, v in TRANSPORT_DOMAIN.items()}, "defaults": TRANSPORT_DEFAULTS,
        "presets": list(PRESETS), "stark_V": list(STARK_V), "screenings": list(SCREENINGS),
        "caveat": TRANSPORT_CAVEAT, "nitride_caveat": NITRIDE_CAVEAT,
        "zhang": {"slope_meV_per_V": -10.0, "label": "Zhang 2016, -10 meV/V: labelled guide only, non-gating",
                  "tag": "V"},
        "tags": {**transport_output_tags(), "inputs": "A", "nitride": "A", "E_X": "A", "slope": "A"},
        "tag_notes": {"inputs": "I, n_dot, aperture: qcap-staged.yaml:54-57, 73-77 [A]; w, dE_WL [E] (transport docstring item 7)",
                      "nitride": "NitrideDiode: Zhang 2016 defaults [V source, E transfer]",
                      "E_X": "levels(): geometry and screening are [A]"},
        "calls": ["fsim_core.transport.evaluate_injection", "fsim_core.transport.Diode.depletion",
                  "fsim_core.nitride_transport.evaluate_injection", "fsim_core.nitride_stark.resolve_bias",
                  "fsim_core.nitride_levels.levels", "fsim_core.nitride_stark.stark_derivatives"],
    }


# =================================================================== sde (Zhao)
SDE_DOMAIN = {"I": (1.2, 6.0), "F_pump": (0.0, 1.0), "n_runs": (8, 60), "t_end": (1e-9, 3e-9)}
SDE_DT = 1.0e-13
SDE_I_GRID = (1.2, 1.6, 2.0, 3.0, 4.0, 5.0, 6.0)
SDE_F_PUMPS = (1.0, 0.08)
SDE_FORBIDDEN = ("disable_stim", "beta_sp", "N_dots", "E_a_ES_meV", "E_a_GS_meV", "eps_gain", "params")
ZHAO_CARD_POINTS = {"normal": (1.0224, 0.003), "quiet": (0.9823, 0.006)}


def sde_frame(I_over_Ith, F_pump, n_runs, t_end, dt, seed):
    """Job body: one seeded g2_vs_pump frame with the tuned default parameters."""
    from fsim_core import sde
    r = sde.g2_vs_pump(sde.QDLaserParams(), float(I_over_Ith), n_runs=int(n_runs), t_end=float(t_end),
                       dt=float(dt), seed=int(seed), F_pump=float(F_pump))
    return {"I_over_Ith": float(I_over_Ith), "F_pump": float(F_pump), "n_runs": int(n_runs), "t_end": float(t_end),
            "dt": float(dt), "seed": int(seed), "g2_mean": float(r["g2_mean"]), "g2_se": float(r["g2_se"]),
            "mean_S": float(r["mean_S"]), "R_pump": float(r["R_pump"])}


def sde_li(I_grid):
    """Job body: noise-free L-I curve (steady_state at each I) and the operational threshold."""
    from fsim_core import sde
    P = sde.QDLaserParams()
    rows = [sde.steady_state(P, float(I)) for I in I_grid]
    return {"I": [float(I) for I in I_grid], "S": [r["S"] for r in rows], "rho_ES": [r["rho_ES"] for r in rows],
            "rho_GS": [r["rho_GS"] for r in rows], "residual": [r["residual"] for r in rows],
            "I_th": float(sde.find_threshold(P))}


def sde_trace(I_over_Ith, F_pump, t_end, dt, seed, max_points=900):
    """Job body: one simulate() S(t) trace, decimated for display; burn-in reported, not hidden."""
    from fsim_core import sde
    P = sde.QDLaserParams()
    r = sde.simulate(P, float(I_over_Ith), float(t_end), float(dt), seed=int(seed), F_pump=float(F_pump))
    t, S = np.asarray(r["t"]), np.asarray(r["S"])
    stride = max(1, int(math.ceil(len(t) / max_points)))
    mean_S, var_S = float(r["mean_S"]), float(r["var_S"])
    return {"t_ns": (t[::stride] * 1e9).tolist(), "S": S[::stride].tolist(), "n_steps": int(len(t)),
            "n_burn": int(r["n_burn"]), "burnin_frac": float(sde.BURNIN_FRAC), "burn_t_ns": float(t[int(r["n_burn"])] * 1e9) if int(r["n_burn"]) < len(t) else None,
            "mean_S": mean_S, "var_S": var_S, "g2_0": float(r["g2_0"]), "g2_0_ext": float(r["g2_0_ext"]),
            "fano": var_S / mean_S if mean_S else None,
            "fano_ext": float(r["var_S_ext"]) / float(r["mean_S_ext"]) if r["mean_S_ext"] else None,
            "eta_ext": float(P.eta_ext), "I_over_Ith": float(I_over_Ith), "F_pump": float(F_pump), "seed": int(seed)}


def _read_csv(path):
    with open(path, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def sde_card():
    """The Zhao card points and the committed fit files, verbatim (never refit here)."""
    import yaml
    card = yaml.safe_load((CARDS / "zhao.yaml").read_text(encoding="utf-8"))
    pts = card["data"]["g2_vs_pump"]["points"]
    cmp_rows = _read_csv(OUT / "zhao" / "zhao_fit_comparison.csv")
    curve = _read_csv(OUT / "zhao" / "zhao_g2_vs_pump_curve.csv")
    num = lambda v: float(v) if v not in (None, "") else None
    return {
        "card_points": [{"pump": p["pump"], "I_over_Ith": float(p["I_over_Ith"]), "g2": float(p["g2"]), "err": float(p["err"])} for p in pts],
        "card_tag": card["data"]["g2_vs_pump"].get("tag", "V"),
        "comparison": [{k: (v if k in ("pump", "verdict") else num(v)) for k, v in r.items()} for r in cmp_rows],
        "curve": [{k: num(v) for k, v in r.items()} for r in curve],
        "files": ["out/zhao/zhao_fit_comparison.csv", "out/zhao/zhao_g2_vs_pump_curve.csv", "cards/zhao.yaml"],
    }


def sde_meta():
    return {
        "domain": {k: list(v) for k, v in SDE_DOMAIN.items()}, "dt": SDE_DT, "I_grid": list(SDE_I_GRID),
        "F_pumps": list(SDE_F_PUMPS), "defaults": {"I": 4.0, "n_runs": 8, "t_end": 1e-9, "seed": 1, "F_pump": 1.0},
        "caveat": SDE_CAVEAT, "axis_note": SDE_AXIS_NOTE, "burnin_frac": 0.1,
        "tags": {"model": "E", "card": "V", "temperature": "V", "result": "E"},
        "tag_notes": {"model": "N_dots, E_a_ES, tau_p, I_ref, F_pump: all [E] (sde.py:137, 159, 172, 196, 205); T 300 K [V]",
                      "card": "cards/zhao.yaml:34-35 [V]"},
        "calls": ["fsim_core.sde.steady_state", "fsim_core.sde.find_threshold", "fsim_core.sde.g2_vs_pump",
                  "fsim_core.sde.simulate"],
    }


# =================================================================== routes
def register(app):
    from flask import Response, request

    from .serialize import dumps

    def _out(obj, status=200):
        return Response(dumps(obj), status=status, mimetype="application/json")

    @app.errorhandler(ExploreError)
    def _explore_error(exc):
        return _out({"error": exc.message, **exc.extra}, exc.status)

    def jm():
        return app.config["JOB_MANAGER"]

    def _submit(fn, kwargs, eta, label):
        try:
            return jm().submit_call(f"fsim_studio.api_explore_b:{fn}", kwargs, eta, label=label)
        except ValueError as exc:
            raise ExploreError(str(exc))

    def _body():
        data = request.get_json(silent=True)
        if not isinstance(data, dict):
            raise ExploreError("request body must be a JSON object")
        return data

    # ---------------------------------------------------------------- phonon
    @app.get("/api/explore/phonon/meta")
    def phonon_meta_route():
        return _out(phonon_meta())

    @app.get("/api/explore/phonon/zt")
    def phonon_zt_route():
        a = _phonon_args(request.args)
        return _out(phonon_zt(a["l_xy"], a["l_z"], a["alpha"]))

    @app.get("/api/explore/phonon/live")
    def phonon_live_route():
        """Brief 3(e): T, l_xy, l_z and alpha select PRECOMPUTED frames (the phi cache is per process, ~1 s
        cold per new (params, T)); only Gamma_zpl, w, kappa, F_cav and delta move per request. So this route
        never integrates in the Flask thread: T must be one of the frame temperatures and the one-frame
        evaluation runs on the pool through the same phonon_frames job (cached by its arguments)."""
        a = _phonon_args(request.args)
        T = _num("T", request.args.get("T"), *PHONON_DOMAIN["T"], default=PHONON_DEFAULTS["T"])
        if T not in PHONON_T_FRAMES:
            raise ExploreError(f"T={T:g} K is not a precomputed frame temperature; use one of "
                               f"{[int(t) for t in PHONON_T_FRAMES]} (POST /api/explore/phonon/frames)")
        resp = _submit("phonon_frames", {**a, "T_list": [float(T)]}, 1.2, "phonon")
        return _out({**resp, "T": T, "inputs": a}, 200 if resp.get("cached") else 202)

    @app.post("/api/explore/phonon/frames")
    def phonon_frames_route():
        a = _phonon_args(_body())
        Ts = [float(t) for t in PHONON_T_FRAMES]
        resp = _submit("phonon_frames", {**a, "T_list": Ts}, 1.2 * len(Ts), "phonon")
        return _out({**resp, "T": Ts, "inputs": a}, 200 if resp.get("cached") else 202)

    # ---------------------------------------------------------------- transport
    def _transport_args(src):
        g = src.get
        d = TRANSPORT_DEFAULTS
        return {
            "I": _num("I_uA", g("I_uA"), *TRANSPORT_DOMAIN["I_uA"], default=d["I_uA"]),
            "T": _num("T", g("T"), *TRANSPORT_DOMAIN["T"], default=d["T"]),
            "n_dot": _num("n_dot", g("n_dot"), *TRANSPORT_DOMAIN["n_dot"], default=d["n_dot"]),
            "aperture": _num("aperture", g("aperture"), *TRANSPORT_DOMAIN["aperture"], default=d["aperture"]),
            "tau_pulse": _num("tau_pulse", g("tau_pulse"), *TRANSPORT_DOMAIN["tau_pulse"], default=d["tau_pulse"]),
            "w": _num("w", g("w"), *TRANSPORT_DOMAIN["w"], default=d["w"]),
            "dE_WL": _num("dE_WL", g("dE_WL"), *TRANSPORT_DOMAIN["dE_WL"], default=d["dE_WL"]),
            "tau_rad": _num("tau_rad", g("tau_rad"), *TRANSPORT_DOMAIN["tau_rad"], default=d["tau_rad"]),
        }

    def _preset(src):
        p = src.get("preset") or "red"
        if p not in PRESETS:
            raise ExploreError(f"preset must be one of {list(PRESETS)}")
        return p

    @app.get("/api/explore/transport/meta")
    def transport_meta_route():
        return _out(transport_meta())

    @app.get("/api/explore/transport/point")
    def transport_point_route():
        a = _transport_args(request.args)
        ex = _num("E_X_eV", request.args.get("E_X_eV"), *TRANSPORT_DOMAIN["E_X_eV"], allow_none=True)
        return _out(transport_point(_preset(request.args), a["I"], a["T"], a["n_dot"], a["aperture"], a["tau_pulse"],
                                    a["w"], a["dE_WL"], a["tau_rad"], ex))

    @app.get("/api/explore/transport/depletion")
    def transport_depletion_route():
        T = _num("T", request.args.get("T"), *TRANSPORT_DOMAIN["T"], default=TRANSPORT_DEFAULTS["T"])
        return _out(transport_depletion(_preset(request.args), T))

    @app.get("/api/explore/transport/nitride_point")
    def nitride_point_route():
        a = _transport_args(request.args)
        T = _num("T", request.args.get("T"), 230.0, 350.0, default=TRANSPORT_DEFAULTS["T_nitride"])
        pol = int(_num("polarity", request.args.get("polarity"), -1, 1, default=1))
        if pol not in (-1, 1):
            raise ExploreError("polarity must be +1 or -1")
        ext = _num("ext", request.args.get("ext"), *TRANSPORT_DOMAIN["ext"], default=0.0)
        return _out(nitride_point(a["I"], T, a["n_dot"], a["aperture"], a["tau_pulse"], a["w"], a["dE_WL"],
                                  a["tau_rad"], pol, ext))

    @app.post("/api/explore/transport/sweep")
    def transport_sweep_route():
        body = _body()
        a = _transport_args(body)
        kw = {"preset": _preset(body), "T": a["T"], "n_dot": a["n_dot"], "aperture": a["aperture"],
              "tau_pulse": a["tau_pulse"], "w": a["w"], "dE_WL": a["dE_WL"], "tau_rad": a["tau_rad"],
              "I_lo": 1e-3, "I_hi": 100.0, "n": 50}
        resp = _submit("transport_sweep", kw, 1.0, "tsweep")
        return _out({**resp, "inputs": kw}, 200 if resp.get("cached") else 202)

    @app.post("/api/explore/transport/stark")
    def transport_stark_route():
        body = _body()
        ori = body.get("orientation") or "c_plane"
        if ori not in ("c_plane", "a_plane"):
            raise ExploreError("orientation must be c_plane or a_plane")
        pol = int(_num("polarity", body.get("polarity"), -1, 1, default=1))
        if pol not in (-1, 1):
            raise ExploreError("polarity must be +1 or -1")
        T = _num("T", body.get("T"), 230.0, 350.0, default=TRANSPORT_DEFAULTS["T_nitride"])
        ext = _num("ext", body.get("ext"), *TRANSPORT_DOMAIN["ext"], default=0.0)
        scr = _num_list("screenings", body.get("screenings"), SCREENINGS, len(SCREENINGS))
        jobs = []
        for s in scr:
            s = _num("screening_fraction", s, 0.0, 1.0)
            if s not in SCREENINGS:
                raise ExploreError(f"screening_fraction must be one of {list(SCREENINGS)}")
            kw = {"orientation": ori, "screening_fraction": s, "polarity": pol, "T": T, "ext": ext,
                  "V_list": list(STARK_V)}
            resp = _submit("stark_trace", kw, 2.5, "stark")
            jobs.append({"screening_fraction": s, **resp})
        return _out({"traces": jobs, "orientation": ori, "polarity": pol, "T": T, "ext": ext})

    # ---------------------------------------------------------------- sde
    @app.get("/api/explore/sde/meta")
    def sde_meta_route():
        return _out(sde_meta())

    @app.get("/api/explore/sde/card")
    def sde_card_route():
        return _out(sde_card())

    def _sde_common(body):
        _reject(body, SDE_FORBIDDEN)
        n_runs = int(_num("n_runs", body.get("n_runs"), *SDE_DOMAIN["n_runs"], default=8))
        t_end = _num("t_end", body.get("t_end"), *SDE_DOMAIN["t_end"], default=1e-9)
        seed = int(_num("seed", body.get("seed"), 0, 2**31 - 1, default=1))
        return n_runs, t_end, seed

    @app.post("/api/explore/sde/frames")
    def sde_frames_route():
        body = _body()
        n_runs, t_end, seed = _sde_common(body)
        grid = _num_list("I_grid", body.get("I_grid"), SDE_I_GRID, 12)
        fps = _num_list("F_pumps", body.get("F_pumps"), SDE_F_PUMPS, 3)
        eta = 0.7 * (n_runs / 8.0) * (t_end / 1e-9)
        frames = []
        for F in fps:
            F = _num("F_pump", F, *SDE_DOMAIN["F_pump"])
            for I in grid:
                I = _num("I", I, *SDE_DOMAIN["I"])
                kw = {"I_over_Ith": I, "F_pump": F, "n_runs": n_runs, "t_end": t_end, "dt": SDE_DT, "seed": seed}
                resp = _submit("sde_frame", kw, eta, "sde")
                frames.append({"I": I, "F_pump": F, **resp})
        return _out({"frames": frames, "n_runs": n_runs, "t_end": t_end, "seed": seed})

    @app.post("/api/explore/sde/li")
    def sde_li_route():
        body = _body()
        _reject(body, SDE_FORBIDDEN)
        grid = [_num("I", x, 0.5, SDE_DOMAIN["I"][1])
                for x in _num_list("I_grid", body.get("I_grid"), [1.0, 1.2, 1.6, 2.0, 2.5, 3.0, 4.0, 5.0, 6.0], 16)]
        resp = _submit("sde_li", {"I_grid": grid}, 0.2 * len(grid) + 0.4, "sdeli")
        return _out({**resp, "I_grid": grid})

    @app.post("/api/explore/sde/trace")
    def sde_trace_route():
        body = _body()
        n_runs, t_end, seed = _sde_common(body)
        I = _num("I", body.get("I"), *SDE_DOMAIN["I"], default=4.0)
        F = _num("F_pump", body.get("F_pump"), *SDE_DOMAIN["F_pump"], default=1.0)
        resp = _submit("sde_trace", {"I_over_Ith": I, "F_pump": F, "t_end": t_end, "dt": SDE_DT, "seed": seed},
                       0.7 * (t_end / 1e-9), "sdetr")
        return _out({**resp, "I": I, "F_pump": F, "t_end": t_end, "seed": seed})

    return app


__all__ = ["register", "ExploreError", "phonon_state", "phonon_frames", "transport_sweep", "stark_trace",
           "sde_frame", "sde_li", "sde_trace"]
