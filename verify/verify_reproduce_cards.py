"""Gating reproduction of published measurements by the PRODUCTION model, against an
EXPECTATION LEDGER (verify/data/reproduce_ledger.yaml).  devlib-vf2.

Closes p0-verifiability defects H3/H4/H5 and M4/M5/M6: device-level model outputs for the
existing literature cards are compared with the published numbers under tolerances that
were written down in the ledger BEFORE the first run (class C calibration, U unfitted,
H held-out, KD known deviation).  Each item carries `expect: pass | fail`; a documented
model miss is a gated `expect: fail`, so an unexpected pass (the model changed) and an
unexpected fail (a regression) both turn the exit code to 1.

Evaluators call production code only (fsim_core.fitting.fit_phase0, fsim_core.device,
fsim_core.sde, fsim_core.dot_levels, fsim_core.nitride_levels, scripts.run_nitride_cavity,
fsim_core.nitride_nanowire_*); they do not re-implement the model.  scripts/run_v1_validation.py
runs its computation at import time (writes out/validation), so its CALL SEQUENCE is
replicated here instead of imported.

Usage:
    python verify/verify_reproduce_cards.py            # gating run (non-slow items)
    python verify/verify_reproduce_cards.py --all      # also the items marked `slow: true`
    python verify/verify_reproduce_cards.py --self-test

Exit 0 iff every item's outcome equals its ledger `expect`.
"""
from __future__ import annotations

import copy
import json
import math
import sys
import time
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

LEDGER_PATH = ROOT / "verify" / "data" / "reproduce_ledger.yaml"
HC_EV_NM = 1239.841984  # eV nm

# ----------------------------------------------------------------------------- evaluators
# Each returns {key: float}; a "<key>__se" entry carries a Monte-Carlo standard error.
# Groups listed in SLOW_GROUPS run only with --all.
SLOW_GROUPS: set[str] = set()


def _chatz_card():
    from fsim_core.card import load_card
    return load_card(ROOT / "cards" / "chatzarakis.yaml")


def ev_chatz_va():
    """V-a joint fit in memory (production fit_phase0), per-point model g2."""
    from fsim_core.fitting import fit_phase0
    fit = fit_phase0(_chatz_card())
    return {f"g2_T{int(r['T'])}": float(m) for r, m in zip(fit.data, fit.model_g2)}


def _ensure_class_proxy():
    """device.class_proxy_params() resolves the V-a bundle through fsim_core.fit_bundle (out/ if
    present, else the tracked snapshot fsim_core/data/phase0_fit_params.json), so a clean clone
    runs every item; this only reports which file was used."""
    from fsim_core.fit_bundle import fit_params_source
    return fit_params_source()


def _base_design(DeviceDesign, lineshape, geom, alpha=None, rho_one=True):
    # call sequence of scripts/run_v1_validation.py::base_design
    d = DeviceDesign()
    d.drive.V = 0.0
    d.drive.I_uA = 0.0
    d.drive.duty = 0.0
    d.drive.b_e = 0.0
    d.drive.mode = "PL"
    d.cavity.enabled = False
    d.filter.enabled = True
    d.filter.auto_w = False
    d.dot.lineshape = lineshape
    ph = dict(geom)
    if alpha is not None:
        ph["alpha_ps2"] = alpha
    d.dot.phonon = ph
    if rho_one:
        d.ret.b0 = 1e-12
        d.ret.beta = 1e-12
    return d


def _run_point(DeviceDesign, evaluate, d, T, w, dx, mu, delta, gamma_target=None):
    # call sequence of scripts/run_v1_validation.py::run_point
    d.thermal.T_hs = T
    d.filter.w = w
    d.filter.dx = dx
    d.drive.mu = mu
    d.dot.delta_xx = delta
    if gamma_target is not None:
        from fsim_core.device import class_proxy_params
        from fsim_core.spectral import gamma_of_T
        p = class_proxy_params()
        g_proxy = float(gamma_of_T(T, p["gamma0"], p["a_ac"], p["b_lo"], p["E_lo"]))
        d.dot.gamma_scale = gamma_target / g_proxy
    return float(evaluate(d, T_grid=[T])["scalars"]["g2_op"])


def ev_chatz_va_ibm():
    """devlib-r2-pha2 D1: the opt-in IBM-calibrated V-a refit in memory (production fit_phase0, lineshape="ibm")."""
    from fsim_core.fitting import fit_phase0
    fit = fit_phase0(_chatz_card(), lineshape="ibm")
    return {f"g2_T{int(r['T'])}": float(m) for r, m in zip(fit.data, fit.model_g2)}


def _chatz_dev(ibm_bundle):
    from fsim_core.device import DeviceDesign, evaluate
    _ensure_class_proxy()
    card = yaml.safe_load((ROOT / "cards" / "chatzarakis2023.yaml").read_text())
    dev = card["device"]
    delta_xx, r_xx = dev["delta_xx"], None
    if ibm_bundle:
        from fsim_core.fit_bundle import load_fit_params
        bp = load_fit_params("ibm")["params"]
        delta_xx, r_xx = bp["delta_xx"], bp["r_xx"]   # the bundle's own XX offset / width ratio
    out = {}
    for pt in card["data"]["g2_vs_T"]["points"]:
        d = _base_design(DeviceDesign, "ibm", dev["geometry"], alpha=dev["phonon_alpha_ps2"], rho_one=False)
        if ibm_bundle:
            d.dot.fit_bundle = "ibm"   # non-field opt-in attribute read by device.evaluate (D1)
            d.dot.r_xx = r_xx
        out[f"g2_T{int(pt['T'])}"] = _run_point(DeviceDesign, evaluate, d, float(pt["T"]), float(pt["w"]),
                                                float(pt["dx"]), dev["mu"], delta_xx)
    return out


def ev_chatz_dev_ibm():
    """Device tier on the IBM-calibrated bundle (devlib-r2-pha2: re-pointed from the Lorentzian bundle)."""
    return _chatz_dev(True)


def ev_chatz_dev_ibm_lor():
    """LEGACY wiring kept visible: Lorentzian bundle, full measured FWHM used as the IBM ZPL width."""
    return _chatz_dev(False)


def _centre(rng):
    return 0.5 * (float(rng[0]) + float(rng[1]))


def ev_laf_mid():
    """Pre-registered mid-point: centre of every [E] range (delta_xx, Gamma, w, mu)."""
    from fsim_core.device import DeviceDesign, evaluate
    _ensure_class_proxy()
    card = yaml.safe_load((ROOT / "cards" / "laferriere2023.yaml").read_text())
    dev = card["device"]
    out = {}
    for pt in card["data"]["g2_vs_T"]["points"]:
        T = float(pt["T"])
        key = str(int(T))
        mu_rng = dev["mu"]["range_77"] if T < 100 else dev["mu"]["range_high"]
        d = _base_design(DeviceDesign, "ibm", dev["geometry"])
        out[f"g2_T{int(T)}"] = _run_point(
            DeviceDesign, evaluate, d, T, _centre(dev["w_meV"][key]["range"]), 0.0, _centre(mu_rng),
            _centre(dev["delta_xx"]["range"]), _centre(dev["gamma_meV"][key]["range"]))
    return out


def ev_reischle():
    """Call-for-call copy of verify_rt_edge_papers._xcheck_reischle2008_g2_80K's fixture."""
    from fsim_core.device import DeviceDesign, evaluate
    d = DeviceDesign()
    d.ret.mode = "confinement"
    d.ret.preset = "InP/AlGaInP0.2/AlGaInP0.55 on GaAs"
    d.drive.mode = "EL-transport"
    d.drive.diode = {"preset": "red", "tau_pulse_ns": 0.1}
    d.drive.duty = 0.008
    d.drive.n_dot_cm2 = 1e10
    d.thermal.T_hs = 80.0
    r = evaluate(d, [80.0])
    return {"g2_op": float(r["scalars"]["g2_op"])}


# devlib-rc2: faithful DC comparison for Reischle 2008 (QD C, 80 K).  Constants are paper inputs, fixed before the first run.
REISCHLE_TAU_C_NS = 1.0       # [V] Reischle 2008: system recovery time ~1 ns from the rise of the antibunching dip (prose, rounded)
REISCHLE_TAU_C_DR_NS = 1.6    # [DR] r2-physics-a Q4: inverting the paper's own pair (raw 0.43, deconvolved 0.25, tau_i 0.5 ns)
REISCHLE_R_OVER_D = 0.2       # [A] r/D_eff, sub-saturation (moves g2_dot0 by 0.004 only); D_eff + r = 1/tau_c
REISCHLE_RHO_MEASURED = 0.88  # [V] Reischle 2008: measured rho = S/(S+B), QD C (80 K); TRANSCRIBED, not predicted
REISCHLE_IRF_FWHM_PS = 693.0  # [V] 2 ln2 x 0.5 ns: the paper's fit-model IRF C exp(-|tau/0.5 ns|), exponential shape
REISCHLE_APERTURE_UM = 0.46   # [V] nanosphere shadow-mask aperture
REISCHLE_R_FLOOR_NS = 0.01    # [A] smallest realisable pump rate (r = 1/tau_c - gamma_X - k_X < 0 is not physical)
REISCHLE_I_PROBE_UA = 0.001   # [A] low-current probe point (the paper's 11.9 mA runs the 1 um default mesa to thermal runaway)


def ev_reischle_dc():
    """Card-faithful DC comparison.  device.evaluate (EL-transport, drive.cw) supplies gamma_X, the 80 K retention S
    (-> sink k_X, k_XX = 2 k_X as cw_g2.escape_rates_from_retention) and eps; the sink extraction is cross-checked
    against device.evaluate's own g2_cw0 / g2_cw0_raw at the probe point (raises if it disagrees).  The comparison
    itself (devlib-r2-pha2, r2-physics-a Q4) maps S -> (Gamma_eff, k_eff) by cw_g2.rates_fast_exchange at the dip-width
    D_eff and gives the intrinsic, with-background and IRF-convolved g2(0) via cw_g2.cw_report."""
    from fsim_core import cw_g2
    from fsim_core.device import DeviceDesign, evaluate
    d = DeviceDesign()
    d.ret.mode = "confinement"
    d.ret.preset = "InP/AlGaInP0.2/AlGaInP0.55 on GaAs"
    d.drive.mode = "EL-transport"
    d.drive.diode = {"preset": "red"}
    d.drive.cw = True
    d.drive.cw_irf_shape = "exponential"
    d.drive.cw_irf_fwhm_ps = REISCHLE_IRF_FWHM_PS
    d.drive.duty = 1.0
    d.drive.n_dot_cm2 = 1e10
    d.drive.I_uA = REISCHLE_I_PROBE_UA
    d.aperture.diameter_um = REISCHLE_APERTURE_UM
    d.thermal.T_hs = 80.0
    sc = evaluate(d, [80.0])["scalars"]
    gX, S, eps = float(sc["cw_gamma_X_ns"]), float(sc["S_resolved"]), float(sc["eps_op"])
    r_dev, rho_dev = float(sc["cw_r_ns"]), float(sc["cw_rho_op"])
    kX = gX * (1.0 / S - 1.0)   # perfect-sink mapping: used only for the device.evaluate cross-check below
    kXX = 2.0 * kX

    def rep(r, rho):
        return cw_g2.cw_report(r, gX, 2.0 * gX, kX, kXX, 1.0, eps, rho, REISCHLE_IRF_FWHM_PS,
                               tau_max_ns=10.0, irf_shape="exponential")

    chk = rep(r_dev, rho_dev)
    for mine, theirs, name in ((chk["g2_meas0"], float(sc["g2_cw0"]), "g2_cw0"),
                               (chk["g2_raw0"], float(sc["g2_cw0_raw"]), "g2_cw0_raw")):
        if abs(mine - theirs) > 1e-3 * max(1.0, abs(theirs)):
            raise RuntimeError(f"reischle_dc extraction disagrees with device.evaluate on {name}: {mine} vs {theirs}")
    # r2-physics-a Q4: the sink mapping (kX above) is kept ONLY for the cross-check against device.evaluate.  The
    # comparison itself uses the ch09 fast-exchange mapping: same S, Gamma_eff = S D_eff, k_eff = (1 - S) D_eff, with the
    # dip width fixing D_eff + r = 1/tau_c at r/D_eff = 0.2 [A]; tau_c is the bracket [1.0 (prose, V), 1.6 (DR)].
    def fx(tau_c, rho):
        D = 1.0 / (tau_c * (1.0 + REISCHLE_R_OVER_D))
        G_eff, k_eff = cw_g2.rates_fast_exchange(gX, S, D)
        return cw_g2.cw_report(REISCHLE_R_OVER_D * D, G_eff, 2.0 * G_eff, k_eff, 2.0 * k_eff, 1.0, eps, rho,
                               REISCHLE_IRF_FWHM_PS, tau_max_ns=10.0, irf_shape="exponential"), D, G_eff, k_eff

    f10_hi, D10, G10, k10 = fx(REISCHLE_TAU_C_NS, 1.0)
    f10_m, _, _, _ = fx(REISCHLE_TAU_C_NS, REISCHLE_RHO_MEASURED)
    f16_m, D16, _, _ = fx(REISCHLE_TAU_C_DR_NS, REISCHLE_RHO_MEASURED)
    f10_d, _, _, _ = fx(REISCHLE_TAU_C_NS, rho_dev)
    return {"g2_dot0": float(f10_hi["g2_dot0"]), "g2_meas0_rho088": float(f10_m["g2_meas0"]),
            "g2_raw0_rho088_tauc16": float(f16_m["g2_raw0"]), "g2_raw0_rho088_tauc10": float(f10_m["g2_raw0"]),
            "g2_raw0_rho_model": float(f10_d["g2_raw0"]),
            "D_eff_tauc10_ns": float(D10), "D_eff_tauc16_ns": float(D16), "Gamma_eff_tauc10_ns": float(G10),
            "k_eff_tauc10_ns": float(k10), "gamma_X_ns": gX, "k_X_sink_ns": float(kX), "S": S,
            "eps": eps, "rho_model": rho_dev, "tau_dip_ns": float(f10_m["tau_dip"])}


def ev_zhao():
    """scripts/run_zhao_fit.py main computation through production sde.py."""
    from fsim_core.sde import QDLaserParams, find_threshold, g2_vs_pump
    params = QDLaserParams()
    I_op = 4.0 * find_threshold(params)
    out = {}
    for label, F_pump in (("normal", 1.0), ("quiet", 0.08)):
        r = g2_vs_pump(params, I_op, n_runs=60, t_end=2.0e-9, dt=1.0e-13, seed=42, F_pump=F_pump)
        out[label] = float(r["g2_mean"])
        out[label + "__se"] = float(r["g2_se"])
    return out


SSL_GAAS_NM = 0.86  # [DR] <Al> = y/(x+y) = 16/(x+16) = 0.65 -> x = 8.6 A (Chatzarakis 2023 SSL sample)


def ev_dot_levels():
    from dataclasses import replace

    from fsim_core import dot_levels as D
    from fsim_core import materials as M
    P = D.class_presets()

    def lv(name, **kw):
        return D.levels(P[name](**kw))

    L = lv("InGaAs0.5/GaAs/GaAs on GaAs")
    Lc = lv("InGaAs0.5/AlGaAs0.57/AlGaAs0.57 on GaAs")
    S = P["InGaAs0.5/GaAs/AlGaAs0.57 on GaAs"]()
    S.geometry.matrix_thickness_nm = 2.0
    shift = D.levels(S).dE_pair_WL - L.dE_pair_WL
    Lr = lv("InP/GaInP/AlGaInP0.55 on GaAs")
    Lt = lv("InAs/InP/InP")
    Lg = lv("InP/GaAs0.65P0.35/AlGaAs0.4 on GaAs", height_nm=5.5)
    # r2-physics-b Q2: [film, sphere] strain-state brackets (an island lies between; V_e identical)
    Lc_f, Lc_s = D.levels_bracket(P["InGaAs0.5/GaAs/GaAs on GaAs"]())
    Lg_f, Lg_s = D.levels_bracket(P["InP/GaAs0.65P0.35/AlGaAs0.4 on GaAs"](height_nm=5.5))
    Lr_f, Lr_s = D.levels_bracket(P["InP/GaInP/AlGaInP0.55 on GaAs"]())
    Lb = lv("InP/AlGaInP0.2/AlGaInP0.55 on GaAs")   # Bommer class: Pryor VBO transferred through the Al step
    # devlib-rc2 faithful SSL: the paper's <Al> = 65 % sample = GaAs 8.6 A / AlAs 16 A -> 0.86 nm GaAs layer, Al0.65 barrier [DR]
    Sf = P["InGaAs0.5/GaAs/AlGaAs0.57 on GaAs"]()
    Sf = replace(Sf, barrier=M.AlGaAs(0.65))
    Sf.geometry.matrix_thickness_nm = SSL_GAAS_NM
    shift_f = D.levels(Sf).dE_pair_WL - L.dE_pair_WL
    # devlib-rc2 Pryor model-solid VBO (E_v0(dot) - E_v0(matrix) of the module's own materials, +0.168 eV)
    Sr = P["InP/GaInP/AlGaInP0.55 on GaAs"]()
    vbo_ms = Sr.dot.p["vbo"] - Sr.matrix.p["vbo"]
    Lms = D.levels(replace(Sr, vbo_override_eV={"dot-matrix": vbo_ms}))
    # devlib-rc2 Laferriere faithful preset (presets.LEVEL_PRESETS), 300 K, axial-rod vs radial-InP matrix bracket
    from fsim_core import presets as PR
    e_rod = D.levels(PR.level_system("laferriere-inasp-0.68", T=300.0, matrix="rod")).E_X_eV
    e_inp = D.levels(PR.level_system("laferriere-inasp-0.68", T=300.0, matrix="inp")).E_X_eV
    return {"chatz_EX_film_eV": float(Lc_f.E_X_eV), "chatz_EX_sphere_eV": float(Lc_s.E_X_eV),
            "chatz_dEWL_film_meV": float(Lc_f.dE_pair_WL), "chatz_dEWL_sphere_meV": float(Lc_s.dE_pair_WL),
            "gu_EX_film_eV": float(Lg_f.E_X_eV), "gu_EX_sphere_eV": float(Lg_s.E_X_eV),
            "reischle_EX_film_eV": float(Lr_f.E_X_eV), "reischle_EX_sphere_eV": float(Lr_s.E_X_eV),
            "bommer_Vh_meV": float(Lb.V_h),
            "ssl_shift_faithful_meV": float(shift_f), "pryor_Ve_modelsolid_meV": float(Lms.V_e),
            "laf_faithful_lo_eV": float(min(e_rod, e_inp)), "laf_faithful_hi_eV": float(max(e_rod, e_inp)),"chatz_EX_eV": float(L.E_X_eV), "chatz_dEWL_meV": float(L.dE_pair_WL),
            "ssl0_dEWL_meV": float(Lc.dE_pair_WL), "ssl_shift_meV": float(shift),
            "reischle_EX_eV": float(Lr.E_X_eV), "pryor_Ve_meV": float(Lr.V_e),
            "laf_EX_eV": float(Lt.E_X_eV), "gu_EX_eV": float(Lg.E_X_eV)}


ZHANG_T_K = 10.0   # [V] Zhang et al., APL 103, 192114 (2013), Fig. 2: 10 K lifetime and ZPL


def _zhang_extra(NL, zh):
    """r2-physics-b Q3b: strain x screening brackets and the one-parameter partial-relaxation fit.
    Brackets: lifetime between the field-free end (relaxed, or fully strained with screening 1; the smaller tau)
    and the fully strained unscreened planar field; E_X between the strained-unscreened and the strained-screened
    energies.  Fit: strain_fraction sf solved on the 3.40 ns lifetime ONLY (class C, exact by construction),
    then E_X evaluated at that sf and compared with the 2.95 eV ZPL as a held-out number.  sf is a pillar-
    relaxation [A], never fitted on E_X."""
    from dataclasses import replace

    from scipy.optimize import brentq

    def lt(sf, scr=0.0):
        lv = NL.levels(replace(zh, strain_fraction=sf, screening_fraction=scr), ZHANG_T_K)
        if not lv.valid:
            raise RuntimeError(f"zhang levels invalid at sf={sf}, scr={scr}: {lv.invalid_reasons}")
        r = NL.rates(lv, ZHANG_T_K, tau_rad0_ns=NL.TAU_RAD0_DEFAULT_NS)
        return 1.0 / r["gamma_X0_ns"], lv.E_X_eV

    tau_relaxed, _ = lt(0.0)
    tau_screened, ex_screened = lt(1.0, 1.0)
    tau_unscreened, ex_unscreened = lt(1.0, 0.0)
    sf = brentq(lambda f: lt(f)[0] - NL.ZHANG2013_LIFETIME_NS, 0.0, 1.0, xtol=1e-9)
    tau_fit, ex_fit = lt(sf)
    return {"zhang2013_tau_lo_ns": float(min(tau_relaxed, tau_screened)), "zhang2013_tau_hi_ns": float(tau_unscreened),
            "zhang2013_EX_lo_eV": float(min(ex_unscreened, ex_screened)),
            "zhang2013_EX_hi_eV": float(max(ex_unscreened, ex_screened)),
            "zhang_fit_tau_ns": float(tau_fit), "zhang_fit_lambda_nm": float(HC_EV_NM / ex_fit),
            "zhang_fit_sf": float(sf)}


def ev_nitride():
    from fsim_core import nitride_levels as NL
    from fsim_core.nitride_levels import NitrideDotSystem, levels, rates
    desh = NitrideDotSystem(height_nm=NL.DESHPANDE2013_HEIGHT_NM,
                            radius_nm=NL.DESHPANDE2013_DIAMETER_NM / 2., x_in=NL.DESHPANDE2013_X_IN)
    zh = NitrideDotSystem(height_nm=NL.ZHANG2013_HEIGHT_NM,
                          radius_nm=NL.ZHANG2013_DIAMETER_NM / 2., x_in=NL.ZHANG2013_X_IN)
    lv_d, lv_z300 = levels(desh, 300.), levels(zh, 300.)
    lv_z = levels(zh, ZHANG_T_K)   # r2-physics-b: Zhang's lifetime and ZPL are 10 K measurements
    lv_d25 = levels(desh, 25.)  # devlib-rc2: temperature-matched to the 25 K measurement
    if not (lv_d.valid and lv_z.valid):
        raise RuntimeError(f"planar nitride levels not bound: {lv_d.invalid_reasons} {lv_z.invalid_reasons}")
    r_z = rates(lv_z, ZHANG_T_K, tau_rad0_ns=NL.TAU_RAD0_DEFAULT_NS)
    r_z300 = rates(lv_z300, 300., tau_rad0_ns=NL.TAU_RAD0_DEFAULT_NS)
    zh_ext = _zhang_extra(NL, zh)
    from scripts import run_nitride_cavity as sweep
    row = sweep._row("C000", "comparison", "rectangular", {"T_hs": 300.}, {}, True)
    if not lv_d25.valid:
        raise RuntimeError(f"planar nitride 25 K levels not bound: {lv_d25.invalid_reasons}")
    return {"desh2013_lambda_nm": float(lv_d.lambda_nm), "desh2013_lambda_nm_25K": float(lv_d25.lambda_nm),
            "zhang2013_tau_ns": float(1. / r_z["gamma_X0_ns"]),
            "zhang2013_tau_ns_300K": float(1. / r_z300["gamma_X0_ns"]),
            **zh_ext,
            "desh2014_g2": float(row["g2"])}


DESH2014_K_NR_XX_NS = 3.8   # [A] = 2.4 gamma_XX with gamma_XX = 2/tau_X, tau_X ~ 1.3 ns (r2-physics-b Q3c, [DR]); never tuned


def _desh2014_nanowire():
    """r2-physics-b Q3c: Deshpande-2014-like conditions on the NANOWIRE device path (horizontal_as_built, x_in 0.40,
    relaxed, 300 K, 200 MHz, the sweep script's reference geometry), production device.evaluate.  Baseline, and the
    same row with the opt-in XX-only loss k_nr_XX_ns injected into nitride_nanowire_levels.rates by a wrapper (the
    device module itself exposes no such leaf).  Both counting regimes are reported (the paper text, abstract only,
    does not say which).  Diagnostic, NOT a validation: k_nr_XX was solved from g2 ~ S_XX = 0.29."""
    from fsim_core import device, nitride_nanowire_device as ND
    from scripts import run_nitride_nanowire as R
    out = {}
    orig = ND.rates
    for regime, tag in (("rectangular", "rect"), ("deterministic_pair", "set")):
        p = R.full_defaults("horizontal_as_built")
        p["regime"] = regime
        p["rep_rate_hz"] = 200.0e6
        d, _ = R._design(p)
        out[f"desh2014_nw_g2_{tag}"] = float(device.evaluate(d)["scalars"]["g2_op"])
        ND.rates = lambda *a, **k: orig(*a, **{**k, "k_nr_XX_ns": DESH2014_K_NR_XX_NS})
        try:
            d, _ = R._design(p)
            out[f"desh2014_nw_g2_{tag}_knrXX"] = float(device.evaluate(d)["scalars"]["g2_op"])
        finally:
            ND.rates = orig
    return out


def ev_nanowire():
    from fsim_core.nitride_nanowire_levels import NitrideNanowireSystem, levels
    from fsim_core.nitride_nanowire_photonics import NitrideNanowirePhotonicsParams, response
    desh_sys = NitrideNanowireSystem(height_nm=2., core_radius_nm=12.5, x_in=.25, strain_bound="relaxed")
    lt = levels(desh_sys, 25., z_points=4807, exterior_nm=45.)  # as verify_nitride_nanowire_levels.py (M3)
    if not lt.valid:
        raise RuntimeError(f"nanowire levels not valid: {lt.invalid_reasons}")
    # Claudon-matched GaAs reference, parameters as verify_nitride_nanowire_photonics.py:821-826
    p = NitrideNanowirePhotonicsParams(family="vertical_photonic", n_wire=3.45, n_group_override=4.3,
                                       NA=0.75, bottom_reflectivity=1.0, taper_output_mfr_nm=750.0)
    r = response(p, lambda_nm=950.0, outer_radius_nm=0.22 * 950.0 / 2.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
    # r2-physics-b Q3a: Deshpande 2013 field-free 2D bracket on the nanowire path (25 K)
    def exr(sb, scr):
        sysx = NitrideNanowireSystem(height_nm=2., core_radius_nm=12.5, x_in=.25, strain_bound=sb,
                                     screening_fraction=scr)
        lx = levels(sysx, 25., z_points=4807, exterior_nm=45.)
        if not lx.valid:
            raise RuntimeError(f"nanowire bracket levels not valid: {lx.invalid_reasons}")
        return float(lx.E_X_eV)

    out = {"desh2013_lambda_nm": float(lt.lambda_nm), "claudon_eta": float(r["eta_collection_X"]),
           "desh2013_EX_relaxed_F0_eV": exr("relaxed", 1.0), "desh2013_EX_unrelaxed_scr1_eV": exr("unrelaxed", 1.0)}
    out.update(_desh2014_nanowire())
    return out


EVALUATORS = {"chatz_va": ev_chatz_va, "chatz_va_ibm": ev_chatz_va_ibm, "chatz_dev_ibm": ev_chatz_dev_ibm,
              "chatz_dev_ibm_lor": ev_chatz_dev_ibm_lor, "laf_mid": ev_laf_mid,
              "reischle": ev_reischle, "reischle_dc": ev_reischle_dc, "zhao": ev_zhao, "dot_levels": ev_dot_levels,
              "nitride": ev_nitride, "nanowire": ev_nanowire}


# ----------------------------------------------------------------------------- judging
def judge(item, model, se=None, hi=None):
    """Return (outcome 'pass'|'fail', paper-text, tol-text).  Tolerance logic only; no model-derived numbers."""
    paper, tol = item["paper"], item["tol"]
    kind = tol["kind"]
    if kind == "abs":
        ok = abs(model - paper["value"]) <= tol["value"] + 1e-12
        return ok, f"{paper['value']}", f"+-{tol['value']}"
    if kind == "upper_bound":
        ok = model <= paper["value"] + tol["value"] + 1e-12
        return ok, f"<{paper['value']}", f"<= {paper['value'] + tol['value']:.4g}"
    if kind == "range":
        ok = paper["lo"] - 1e-12 <= model <= paper["hi"] + 1e-12
        return ok, f"[{paper['lo']},{paper['hi']}]", "in range"
    if kind == "energy_meV":
        dE = abs(HC_EV_NM / model - HC_EV_NM / paper["value"]) * 1e3
        return dE <= tol["value"] + 1e-9, f"{paper['value']}", f"+-{tol['value']} meV (dE={dE:.0f})"
    if kind == "rel":
        ok = abs(model - paper["value"]) / abs(paper["value"]) <= tol["value"] + 1e-12
        return ok, f"{paper['value']}", f"+-{tol['value']:.0%}"
    if kind == "combined_sigma":
        if se is None:
            raise ValueError(f"{item['id']}: combined_sigma needs a model standard error")
        lim = tol.get("k", 1.0) * math.hypot(se, paper["err"])
        return abs(model - paper["value"]) <= lim, f"{paper['value']}+-{paper['err']}", f"+-{lim:.4f} ({tol.get('k', 1.0):g} sigma)"
    if kind == "contained":  # measured value inside the model bracket [model, hi]
        if hi is None:
            raise ValueError(f"{item['id']}: contained needs source_hi")
        ok = model - 1e-12 <= paper["value"] <= hi + 1e-12
        return ok, f"{paper['value']} in", f"[{model:.3f}, {hi:.3f}]"
    if kind == "regression":
        ok = abs(model - paper["pin"]) <= tol["value"] + 1e-12
        return ok, f"pin {paper['pin']}", f"+-{tol['value']}"
    raise ValueError(f"{item['id']}: unknown tol kind {kind!r}")


def load_ledger(path=LEDGER_PATH):
    led = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    ids = [it["id"] for it in led["items"]]
    assert len(ids) == len(set(ids)), "duplicate ledger ids"
    for it in led["items"]:
        for k in ("class", "tol", "tol_reason", "expect", "paper", "source"):
            assert k in it, f"{it['id']}: ledger item missing {k}"
        assert it["expect"] in ("pass", "fail"), it["id"]
        assert it["class"] in ("C", "U", "H", "KD"), it["id"]
        if it["expect"] == "fail" and not str(it.get("reason", "")).strip():
            raise AssertionError(f"{it['id']}: expect: fail needs a reason")
    # devlib-r2-pha2: `unavailable:` entries (no judgement possible); a reason is mandatory, never counted as a pass
    for it in led.get("unavailable", []):
        for k in ("id", "class", "desc", "source", "reason"):
            assert str(it.get(k, "")).strip(), f"unavailable item missing {k}: {it}"
        assert it["id"] not in ids, f"{it['id']}: unavailable id clashes with an item"
    return led


def evaluate_ledger(led, run_all=False, verbose=True):
    """Run evaluators lazily per group, judge every item.  Returns list of result dicts."""
    cache, errors = {}, {}
    results = []
    for it in led["items"]:
        group, key = it["source"].split(".", 1)
        res = {"id": it["id"], "class": it["class"], "expect": it["expect"], "item": it}
        if (it.get("slow") or group in SLOW_GROUPS) and not run_all:
            res.update(outcome="skipped", model=None, detail="slow (use --all)")
            results.append(res)
            continue
        if group not in cache and group not in errors:
            t0 = time.time()
            try:
                cache[group] = EVALUATORS[group]()
                if verbose:
                    print(f"  [group {group}: {time.time() - t0:.1f} s]", flush=True)
            except Exception as exc:  # an evaluator crash is a mismatch, never silently skipped
                errors[group] = f"{type(exc).__name__}: {exc}"
        if group in errors:
            res.update(outcome="error", model=None, detail=errors[group])
        else:
            try:
                model = cache[group][key]
                hi = cache[group][it["source_hi"].split(".", 1)[1]] if "source_hi" in it else None
                ok, ptxt, ttxt = judge(it, model, cache[group].get(key + "__se"), hi)
                res.update(outcome="pass" if ok else "fail", model=model, paper_txt=ptxt, tol_txt=ttxt)
            except Exception as exc:
                res.update(outcome="error", model=None, detail=f"{type(exc).__name__}: {exc}")
        results.append(res)
    return results


def exit_code(results):
    return 0 if all(r["outcome"] == r["expect"] for r in results if r["outcome"] != "skipped") else 1


def report(results, led=None):
    print(f"{'id':26s} {'cl':>2s} {'paper':>16s} {'model':>10s} {'tol':>26s} {'outcome':>7s} {'expect':>6s}  match")
    for r in results:
        if r["outcome"] in ("skipped", "error"):
            m = "SKIP" if r["outcome"] == "skipped" else "MISMATCH"
            print(f"{r['id']:26s} {r['class']:>2s} {'-':>16s} {'-':>10s} {'-':>26s} {r['outcome']:>7s} {r['expect']:>6s}  {m}  ({r['detail']})")
            continue
        m = "MATCH" if r["outcome"] == r["expect"] else "MISMATCH"
        print(f"{r['id']:26s} {r['class']:>2s} {r['paper_txt']:>16s} {r['model']:10.4f} {r['tol_txt']:>26s} "
              f"{r['outcome']:>7s} {r['expect']:>6s}  {m}")
    ran = [r for r in results if r["outcome"] != "skipped"]
    skipped = len(results) - len(ran)
    matched = [r for r in ran if r["outcome"] == r["expect"]]
    P = sum(r["outcome"] == "pass" for r in matched)
    F = sum(r["outcome"] == "fail" for r in matched)
    print("\nper-class table (matched/total; expected-pass / expected-fail among matched):")
    for c, name in (("C", "calibration"), ("U", "unfitted"), ("H", "held-out"), ("KD", "known-deviation")):
        sub = [r for r in ran if r["class"] == c]
        if not sub:
            continue
        mt = [r for r in sub if r["outcome"] == r["expect"]]
        print(f"  {c:>2s} {name:16s} {len(mt)}/{len(sub)}   pass {sum(r['outcome'] == 'pass' for r in mt)}   "
              f"known-fail {sum(r['outcome'] == 'fail' for r in mt)}")
    if skipped:
        print(f"  ({skipped} slow items skipped, not counted as passed)")
    for r in ran:
        if r["outcome"] != r["expect"]:
            print(f"MISMATCH {r['id']}: outcome {r['outcome']} but ledger expects {r['expect']}")
    for it in (led or {}).get("unavailable", []):
        print(f"UNAVAILABLE {it['id']} (not counted): {it['reason'][:150]}...")
    print(f"{len(matched)}/{len(ran)} reproduce checks matched ledger ({P} pass, {F} known-fail)")


def self_test(led):
    results = evaluate_ledger(led, run_all=False, verbose=False)
    base = exit_code(results)
    flipped = copy.deepcopy(results)
    victim = next(r for r in flipped if r["outcome"] != "skipped")
    victim["expect"] = "fail" if victim["expect"] == "pass" else "pass"
    code = exit_code(flipped)
    ok = base == 0 and code == 1
    print(f"self-test: baseline exit {base}; after flipping expect of {victim['id']} in memory exit {code}")
    print("self-test " + ("OK (flip detected)" if ok else "FAILED"))
    return 0 if ok else 1


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    led = load_ledger()
    t0 = time.time()
    if "--self-test" in argv:
        rc = self_test(led)
        print(f"runtime {time.time() - t0:.1f} s")
        return rc
    results = evaluate_ledger(led, run_all="--all" in argv)
    report(results, led)
    print(f"runtime {time.time() - t0:.1f} s")
    return exit_code(results)


if __name__ == "__main__":
    sys.exit(main())
