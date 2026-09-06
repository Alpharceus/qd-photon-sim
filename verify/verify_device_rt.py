"""Regression and anchor checks for device.py's RT opt-in paths (integration-a:
anchored linewidth, confinement retention, EL-transport drive). Every check
has an independent physical/analytical anchor or a delegation/wiring anchor
against the module it calls into -- never against device.py's own arithmetic.

Legacy comparison uses the literal fixture verify/data/legacy_device_rt.npz,
confirmed bit-identical to a fresh capture from the committed pre-RT evaluator
(git show HEAD:fsim_core/device.py) before this file's edits; see the spec's
"Previous attempt failed because" note. No new fixture file is added.
"""
from __future__ import annotations
import copy
import sys
import tempfile
import os
import warnings
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import fsim_core.device as devmod
from fsim_core.device import (DeviceDesign, EmissionBlock, evaluate, evaluate_envelope,
                              _stack, _resolve_edge, _compose_aperture_g2, _aperture_lambda,
                              _tracked_material, _diode_from_drive)
from fsim_core.linewidth import LinewidthParams, gamma_anchor
from fsim_core.integrator import retention, g2_from
from fsim_core.loading import loading_probs
from fsim_core.thermal import t_junction
from fsim_core.cavity import tracking_detuning
from fsim_core import dot_levels, transport, waveguide, cw_g2, materials

# transport's fermi_levels() warns at some (preset, T) combinations that this
# file exercises deliberately (e.g. hkust at 77 K); the warning documents an
# approximation already priced into the model, not a defect under test here.
warnings.filterwarnings("ignore", category=UserWarning, module="fsim_core.transport")
# cw_g2.fit_dip's scipy curve_fit transiently probes negative/near-zero decay
# times while converging on some (rate, IRF) combinations this file exercises
# deliberately (e.g. the synthetic zero-width-IRF anchor below); the eventual
# fitted A_dip/tau_dip are unaffected and are not what these checks assert on.
warnings.filterwarnings("ignore", category=RuntimeWarning, module="fsim_core.cw_g2")

checks = []


def ok(name, value):
    checks.append(bool(value))
    print(("ok  " if value else "FAIL") + " " + name)


def raises(name, fn, exc=ValueError):
    try:
        fn()
        ok(name, False)
    except exc:
        ok(name, True)


# =============================================================== 1. legacy regression
# The fixture's design list covers EL/PL, mechanism-free vs staged, cavity
# on/off, hold/mode tracking, and IBM/Lorentzian lineshapes -- none of them
# touch an opt-in RT field, so every value must replay exactly.

def legacy_designs():
    staged = DeviceDesign.load(ROOT / "cards" / "staged-device-design.yaml")
    cavity = DeviceDesign(name="legacy-cavity"); cavity.cavity.enabled = True
    held = DeviceDesign(name="legacy-held"); held.cavity.enabled = True; held.filter.track = "hold"
    ibm = DeviceDesign(name="legacy-ibm"); ibm.dot.lineshape = "ibm"
    pl = DeviceDesign(name="legacy-pl"); pl.drive.mode = "PL"
    el = DeviceDesign(name="legacy-el"); el.drive.mode = "EL"
    return [staged, cavity, held, ibm, pl, el]


fixture = np.load(ROOT / "verify" / "data" / "legacy_device_rt.npz")
T = fixture["T_grid"]
for i, design in enumerate(legacy_designs()):
    result = evaluate(design, T)
    for name, actual in result["curves"].items():
        key = f"{i}/curve/{name}"
        if key in fixture:
            ok(f"legacy {i} curve {name}", np.array_equal(actual, fixture[key], equal_nan=True))
    for name, actual in result["scalars"].items():
        key = f"{i}/scalar/{name}"
        if key in fixture:
            expected = fixture[key].item()
            ok(f"legacy {i} scalar {name}", actual == expected or
               (isinstance(actual, float) and np.isnan(actual) and np.isnan(expected)))

# ------------------------------------------------------------- legacy envelope
# evaluate_envelope with an empty ranged dict must degenerate exactly to a
# single evaluate() call (module docstring); this is a self-consistency
# invariant, true regardless of the RT additions, and exercises the envelope
# path the point-only checks above never touch.
staged = legacy_designs()[0]
env0 = evaluate_envelope(staged, {}, T_grid=T)
pt = evaluate(staged, T)
for name in ("g2", "eps", "rho2", "Tj", "gamma"):
    lo, hi = env0["bands"][name]
    ok(f"legacy envelope degenerate band {name} == point curve",
       np.array_equal(lo, pt["curves"][name], equal_nan=True) and
       np.array_equal(hi, pt["curves"][name], equal_nan=True))
    ok(f"legacy envelope degenerate mid {name} == point curve",
       np.array_equal(env0["mid"][name], pt["curves"][name], equal_nan=True))

# A ranged sweep of a pre-existing (non-opt-in) field must also stay
# legacy-bit-identical; values captured once from the current evaluator and
# pinned here as a literal fixture (spec: "no extra fixture file needed").
env = evaluate_envelope(staged, {"dot.delta_xx": (3.0, 4.0)}, T_grid=T)
EXP_SCALAR_BANDS = {
    "T_c": (171.4122260596644, 188.4510619638275),
    "g2_op": (0.04173098804308051, 0.0521265047237498),
    "eps_op": (0.011859660016308384, 0.02118115106804508),
    "rho_op": (0.9857099383542735, 0.9857099383542735),
    "dT_J": (1.2624912838191733, 1.2624912838191733),
}
for name, (lo_exp, hi_exp) in EXP_SCALAR_BANDS.items():
    lo, hi = env["scalar_bands"][name]
    ok(f"legacy envelope scalar_bands {name}",
       abs(lo - lo_exp) < 1e-9 and abs(hi - hi_exp) < 1e-9)
EXP_G2_LO = np.array([0.02527392, 0.04173099, 0.14223014, 0.66489531, 0.94389968, 0.98969636])
EXP_G2_HI = np.array([0.02545349, 0.0521265, 0.21528127, 0.76907704, 0.95516474, 0.99121054])
lo, hi = env["bands"]["g2"]
ok("legacy envelope g2 band lo", np.allclose(lo, EXP_G2_LO, atol=1e-8))
ok("legacy envelope g2 band hi", np.allclose(hi, EXP_G2_HI, atol=1e-8))
ok("legacy envelope tornado dot.delta_xx",
   abs(env["tornado"]["dot.delta_xx"] - 17.038835904163108) < 1e-6)
ok("legacy envelope n_samples", env["n_samples"] == 2)

# ================================================================ 2. anchored linewidth
d = DeviceDesign(); d.thermal.T_hs = 300.; d.drive.V = 0.; d.dot.linewidth = "anchored"
r = evaluate(d, [300.])
expected = gamma_anchor(300., LinewidthParams(d.dot.gamma0, d.dot.a_ac, d.dot.E_LO, d.dot.gamma300))
ok("anchored linewidth at 300 K", abs(r["scalars"]["gamma_op"] - expected) < 1e-9)
ok("anchored linewidth does not inherit the old class E_LO=18 meV", d.dot.E_LO != 18.0)

# Independent analytical anchor: gamma_anchor's own b_lo_from_anchor rejects
# an anchor below Gamma_0 + a_ac*T (module contract) -- an invalid anchoring
# parameter set must not silently clip or NaN.
raises("anchored linewidth rejects an anchor below Gamma_0 + a_ac*T",
       lambda: gamma_anchor(300., LinewidthParams(gamma0=50.0, a_ac=0.0, E_LO=43.0, gamma300=1.0)))


def _boom(*_a, **_k):
    raise AssertionError("class_proxy_params (legacy fit JSON) must not be read")


d = DeviceDesign(); d.dot.linewidth = "anchored"; d.ret.mode = "confinement"
d.ret.preset = "InP/GaAsP0.4/AlGaAs0.4 on GaAs"
_orig_proxy = devmod.class_proxy_params
devmod.class_proxy_params = _boom
try:
    evaluate(d, [300.])
    ok("fully anchored+confinement RT device reads no implicit proxy", True)
except AssertionError:
    ok("fully anchored+confinement RT device reads no implicit proxy", False)
finally:
    devmod.class_proxy_params = _orig_proxy

# ============================================================ 3. confinement retention
d = DeviceDesign(); d.thermal.T_hs = 300.; d.drive.V = 0.; d.ret.mode = "confinement"
d.ret.preset = "InP/GaAsP0.4/AlGaAs0.4 on GaAs"
r = evaluate(d, [300.])
p = dot_levels.retention_params(dot_levels.levels(dot_levels.class_presets()[d.ret.preset]()),
                                d.ret.tau_rad_ns, d.ret.channel, verbose=False)
S = retention(300., p["a_esc"], p["E_a"], p["b_p"], p["E_b"])
ok("InP confinement E_a published 96+/-25 meV band", 71. <= p["E_a"] <= 121.)
ok("confinement retention Arrhenius value uses all four derived params",
   abs(r["scalars"]["S_resolved"] / S - 1.) < 1e-6)

# S(300 K) anchor for the Bommer et al., JAP 110, 063108 (2011) InP/AlGaInP
# activation energy (96 +/- 7 meV, verify/data/rt_edge_anchors.yaml id
# bommer11-retention-ea): energy alone is not a physical retention -- pair it
# with an independently derived prefactor (dot_levels' detailed-balance
# escape-attempt rate for the matching Bommer-class dot geometry, [E]) and a
# stated radiative lifetime (1 ns, the RetentionBlock.tau_rad_ns class
# convention, [E]), not the unrelated V-a class-proxy fit prefactor.
bommer_name = "InP/AlGaInP0.2/AlGaInP0.55 on GaAs"
p_bommer = dot_levels.retention_params(
    dot_levels.levels(dot_levels.class_presets()[bommer_name](T=300.0)),
    tau_rad_ns=1.0, channel="pair_half", verbose=False)
S_bommer_96 = retention(300., p_bommer["a_esc"], 96.0, p_bommer["b_p"], p_bommer["E_b"])
ok("S(300K) at Bommer Ea=96meV with sourced prefactor+lifetime is a physical fraction",
   0.0 < S_bommer_96 < 1.0)

# Explicit-zero override: ret.overrides is a presence-based map, so an
# explicitly supplied 0.0 must win over the confinement-derived value (unlike
# the proxy path's 0="use class fit" truthy-or shorthand).
d = DeviceDesign(); d.thermal.T_hs = 300.; d.drive.V = 0.; d.ret.mode = "confinement"
d.ret.preset = "InP/GaAsP0.4/AlGaAs0.4 on GaAs"
d.ret.overrides = {"a_esc": 0.0, "b_p": 0.0}
r = evaluate(d, [300.])
ok("explicit-zero retention overrides are preserved (S=1 when a_esc=b_p=0)",
   r["scalars"]["S_resolved"] == 1.0)

# preset + overrides is allowed (overrides lives outside ret.system, so it
# never collides with the preset-vs-inline-stack ambiguity check)
d = DeviceDesign(); d.ret.mode = "confinement"
d.ret.preset = "InP/GaAsP0.4/AlGaAs0.4 on GaAs"; d.ret.overrides = {"b0": 0.001}
try:
    evaluate(d, [300.])
    ok("preset plus overrides is not an ambiguity", True)
except ValueError:
    ok("preset plus overrides is not an ambiguity", False)

# preset + inline stack IS an ambiguity
d = DeviceDesign(); d.ret.mode = "confinement"; d.ret.preset = "x"; d.ret.system = {"dot": "InP"}
raises("reject preset plus inline system", lambda: evaluate(d))

# unsupported confinement: unknown preset name, and an unsupported override key
d = DeviceDesign(); d.ret.mode = "confinement"; d.ret.preset = "does-not-exist"
raises("reject unsupported confinement preset", lambda: evaluate(d, [300.]))
d = DeviceDesign(); d.ret.mode = "confinement"
d.ret.preset = "InP/GaAsP0.4/AlGaAs0.4 on GaAs"; d.ret.overrides = {"not_a_field": 1.0}
raises("reject unsupported ret.overrides key", lambda: evaluate(d, [300.]))

# Inline DotSystem must resolve to exactly the same stack as its matching
# preset (shared serialization; "no silent preset may change the inline
# stack").
d_preset = DeviceDesign(); d_preset.thermal.T_hs = 300.; d_preset.drive.V = 0.
d_preset.ret.mode = "confinement"; d_preset.ret.preset = "InP/GaAsP0.4/AlGaAs0.4 on GaAs"
d_inline = DeviceDesign(); d_inline.thermal.T_hs = 300.; d_inline.drive.V = 0.
d_inline.ret.mode = "confinement"
d_inline.ret.system = {
    "dot": "InP",
    "matrix": {"kind": "GaAsP", "x_p": 0.4},
    "barrier": {"kind": "AlGaAs", "x_al": 0.4},
    "substrate": "GaAs",
    "geometry": {"height_nm": 4.0, "radius_nm": 12.0},
}
r_preset = evaluate(d_preset, [300.])
r_inline = evaluate(d_inline, [300.])
ok("inline DotSystem matches its equivalent preset exactly",
   r_preset["scalars"]["S_resolved"] == r_inline["scalars"]["S_resolved"])

# =============================================================== 4. EL-transport drive
d = DeviceDesign(); d.drive.mode = "EL-transport"; d.drive.diode = {"preset": "hkust"}
d.drive.n_dot_cm2 = 1e10
r = evaluate(d, [d.thermal.T_hs])
sc = r["scalars"]
ok("EL transport injection efficiency positive", sc["eta_inj"] > 0)
ok("EL transport loading.r_dot positive", sc["loading.r_dot"] > 0)
ok("EL transport normalized background finite", np.isfinite(sc["b_e_resolved"]) and sc["b_e_resolved"] >= 0)
ok("GaAs-class built-in voltage", 1.10 <= transport.gaas_homojunction().vbi() <= 1.35)

# GaAs homojunction V_bi with declared doping/temperature, distinguished from
# the forward turn-on voltage (Sze's ideal-diode V_bi is not the operating
# I-V voltage).
gaas = transport.gaas_homojunction(N_A=1e17, N_D=1e17)
V_bi_300 = gaas.vbi(300.0)
V_on_300 = gaas.v_on(300.0)
ok("GaAs homojunction V_bi(N_A=N_D=1e17, 300K) ~1.2 V", abs(V_bi_300 - 1.2) < 0.15)
ok("GaAs homojunction V_bi is distinct from its forward turn-on voltage",
   abs(V_bi_300 - V_on_300) > 0.05)

# Self-heating residual: the operating Tj must be self-consistent with
# t_junction(duty * P_junction_W, ...) at that same Tj (the fixed point the
# evaluate() loop solves for).
a = 0.5 * d.thermal.mesa_diameter_um * 1e-6
st = _stack(d.thermal)
Tj_check = t_junction(d.drive.duty * sc["P_junction_W"], a, st, d.thermal.T_hs)
ok("EL transport self-heating residual is converged", abs(Tj_check - sc["T_j_op"]) < 1e-6)

# zero current: no injected flux is an invalid operating point, not a
# silently-coerced zero.
d0 = DeviceDesign(); d0.drive.mode = "EL-transport"; d0.drive.diode = {"preset": "hkust"}
d0.drive.I_uA = 0.0
r0 = evaluate(d0, [d0.thermal.T_hs])
ok("zero-current EL-transport is an invalid operating point",
   len(r0["scalars"]["invalid_reasons"]) > 0 and np.isnan(r0["scalars"]["g2_op"]))

# conflicting mechanism/transport modes
d = DeviceDesign(); d.drive.mode = "EL-transport"; d.drive.diode = {"preset": "hkust"}
d.drive.mechanism = "poisson-rail"
raises("reject drive.mechanism combined with EL-transport", lambda: evaluate(d, [300.]))

d = DeviceDesign(); d.drive.mode = "PL"; d.drive.diode = {"preset": "red"}
raises("reject PL diode", lambda: evaluate(d))

# ==================================================== 5. photon accounting (once each)
d = DeviceDesign(); d.drive.mode = "EL-transport"; d.drive.diode = {"preset": "hkust"}
d.drive.n_dot_cm2 = 1e10
r = evaluate(d, [d.thermal.T_hs])
sc = r["scalars"]
rate_bg, rate_x = sc["background.rate_bg_window"], sc["background.rate_x"]
ok("background rates are finite and positive", np.isfinite(rate_bg) and np.isfinite(rate_x)
   and rate_bg > 0 and rate_x > 0)
ok("b_e_resolved is the background/X ratio applied exactly once (no unit swap)",
   abs(sc["b_e_resolved"] - rate_bg / rate_x) < 1e-9 * max(1.0, sc["b_e_resolved"]))
ok("background.rate_bg_window and background.rate_x are not swapped with each other",
   rate_bg != rate_x)

# Brightness must apply the resolved loading probabilities, spectral
# transmission and retention S each exactly once -- an independent
# recomputation from the exposed scalars must reproduce it bit-for-bit.
_, P1, P2 = loading_probs(sc["mu_resolved"])
recomputed = (P1 + P2) * sc["t_x_op"] * 1.0 * 1.0 * sc["S_resolved"]  # G=1, beta_sin=1 (cavity off)
ok("pulsed brightness applies loading/spectral/retention factors exactly once",
   abs(recomputed - sc["brightness_per_pulse"]) < 1e-9 * max(1.0, sc["brightness_per_pulse"]))

# ========================================================= 6. round-trip and immutability
d = DeviceDesign(); d.dot.linewidth = "anchored"; d.ret.mode = "confinement"; d.ret.preset = "x"
d.ret.overrides = {"b0": 0.001}
d.drive.mode = "EL-transport"; d.drive.diode = {"preset": "red", "R_s_ohm": 7.}
d.drive.n_dot_cm2 = 3e9
d.provenance = {"manual": {"tag": "V", "note": "round-trip sentinel"}}
with tempfile.NamedTemporaryFile(dir=ROOT, suffix=".yaml", delete=False) as tmp:
    path = Path(tmp.name)
try:
    d.save(path); loaded = DeviceDesign.load(path)
finally:
    os.unlink(path)
ok("YAML preserves RT fields, overrides and provenance", loaded == d)

d2 = DeviceDesign(); d2.dot.linewidth = "anchored"; d2.ret.mode = "confinement"
d2.ret.preset = "InP/GaAsP0.4/AlGaAs0.4 on GaAs"; d2.ret.overrides = {"b0": 0.001}
d2.drive.mode = "EL-transport"; d2.drive.diode = {"preset": "red", "R_s_ohm": 7.}
d2.drive.n_dot_cm2 = 3e9
before = copy.deepcopy(d2)
evaluate(d2, [300.])
ok("evaluate() does not mutate its input design", d2 == before)

# =============================================================== 7. edge emission (emission.type)
EDGE_SYSTEM = {
    "dot": "InP",
    "matrix": {"kind": "GaInP", "x_ga": 0.51},
    "barrier": {"kind": "AlGaInP", "x_al": 0.50, "y_III": 0.51},
    "substrate": "GaAs",
    "geometry": {"height_nm": 3.0, "radius_nm": 10.0},
}
# The confinement-derived (dot_levels.levels) wavelength for a realistic
# InP-dot geometry falls well outside MATERIAL_EXTRA's sparse tabulated
# points (650/668 nm for the phosphide alloys used here); emission.lambda_nm
# pins the explicit 668 nm design point where the stack's own indices are
# actually measured (Schubert et al., JAP 77, 3416 (1995)) [E], per the
# EmissionBlock field's documented fallback.
d_edge = DeviceDesign(); d_edge.thermal.T_hs = 300.; d_edge.drive.V = 0.
d_edge.ret.mode = "confinement"; d_edge.ret.system = EDGE_SYSTEM
d_edge.emission.type = "edge"; d_edge.emission.lambda_nm = 668.0
d_edge.emission.ridge_width_nm = 2000.0; d_edge.emission.etch_depth_nm = 1200.0
r_edge = evaluate(d_edge, [300.])
sc_edge = r_edge["scalars"]
ok("edge emission is eligible at the pinned 668 nm design point",
   not sc_edge["invalid_reasons"] and np.isfinite(sc_edge["edge_eta_total"]))

# Independent delegation anchor: device.py's OWN _resolve_edge, exercised
# directly (mirrors verify_waveguide.py's ".003 < b < .03" ridge beta class
# and its own facet_transmission()/beta_factor() checks -- never against
# device.py's own arithmetic in-line).
edge_direct, lam_direct = _resolve_edge(d_edge.ret, d_edge.emission, sc_edge["T_j_op"])
ok("2 um ridge fixture beta sits in the supported ~1% class envelope [E]",
   0.003 < edge_direct.beta < 0.03)
ok("device.py edge scalars match a fresh _resolve_edge call bit-for-bit",
   edge_direct.eta_total == sc_edge["edge_eta_total"] and lam_direct == sc_edge["edge_lambda_nm"])
T_recomputed = waveguide.facet_transmission(sc_edge["edge_n_eff"])
ok("T_facet matches an independent facet_transmission(n_eff) recomputation",
   abs(T_recomputed - sc_edge["edge_T_facet"]) < 1e-9)
# Facet convention: R_back=None -> front=0.5 exactly (waveguide.edge_emission
# docstring); eta_total = beta * front * eta_prop * eta_NA -- T_facet is NOT
# a second factor (spec: "do not multiply beta, T_facet ... again").
eta_recomputed = sc_edge["edge_beta"] * 0.5 * sc_edge["edge_eta_prop"] * sc_edge["edge_eta_NA"]
ok("eta_total = beta * front(=0.5) * eta_prop * eta_NA, T_facet not reapplied",
   abs(eta_recomputed - sc_edge["edge_eta_total"]) < 1e-9 * sc_edge["edge_eta_total"])

# Lemma 1: emission.type "none" -> "edge" changes brightness only, never the
# intrinsic multiphoton probability (g2_op/eps_op/rho_op untouched).
d_none = DeviceDesign(); d_none.thermal.T_hs = 300.; d_none.drive.V = 0.
d_none.ret.mode = "confinement"; d_none.ret.system = EDGE_SYSTEM
r_none = evaluate(d_none, [300.])
ok("Lemma 1: edge emission leaves g2_op/eps_op/rho_op exactly unchanged",
   r_none["scalars"]["g2_op"] == sc_edge["g2_op"] and
   r_none["scalars"]["eps_op"] == sc_edge["eps_op"] and
   r_none["scalars"]["rho_op"] == sc_edge["rho_op"])
ok("edge emission changes collected brightness",
   r_none["scalars"]["brightness_per_pulse"] != sc_edge["brightness_per_pulse"])

# Scaling common collection alone (propagation length) changes brightness but
# leaves pulsed intrinsic g2 identical.
d_scaled = copy.deepcopy(d_edge); d_scaled.emission.L_um = 900.0
r_scaled = evaluate(d_scaled, [300.])
ok("scaling propagation length changes flux but not intrinsic g2",
   r_scaled["scalars"]["brightness_per_pulse"] != sc_edge["brightness_per_pulse"] and
   r_scaled["scalars"]["g2_op"] == sc_edge["g2_op"])

# Reject conflicting SiN/edge and resonant-cavity coexistence.
d_conflict = DeviceDesign(); d_conflict.emission.type = "edge"; d_conflict.cavity.enabled = True
raises("reject emission.type='edge' combined with cavity.enabled", lambda: evaluate(d_conflict))

# Missing optical-index evidence at the natural (un-pinned) wavelength is an
# ineligible point with a stated reason, not a crash.
d_missing = DeviceDesign(); d_missing.thermal.T_hs = 300.; d_missing.drive.V = 0.
d_missing.ret.mode = "confinement"; d_missing.ret.system = EDGE_SYSTEM
d_missing.emission.type = "edge"  # no lambda_nm override -> natural, out-of-table wavelength
r_missing = evaluate(d_missing, [300.])
ok("missing optical index data at the natural wavelength is ineligible with a reason, not a crash",
   len(r_missing["scalars"]["invalid_reasons"]) > 0
   and np.isnan(r_missing["scalars"]["brightness_per_pulse"]))

# emission.type='edge' requires the shared inline stack (no bare preset-less design).
d_nostack = DeviceDesign(); d_nostack.emission.type = "edge"
r_nostack = evaluate(d_nostack, [300.])
ok("emission.type='edge' with no shared inline stack is ineligible with a reason, not a crash",
   any("ret.preset or ret.system" in reason for reason in r_nostack["scalars"]["invalid_reasons"]))

# ==================================================== 8. continuous aperture composition
BASE_T = 200.0


def _set_thermal(design, T_hs):
    design.thermal.T_hs = T_hs
    return design


d_ap0 = _set_thermal(DeviceDesign(), BASE_T)
r_ap0 = evaluate(d_ap0, [BASE_T])
g2_legacy = r_ap0["scalars"]["g2_op"]

d_apN0 = _set_thermal(DeviceDesign(), BASE_T)
d_apN0.aperture.compose = True; d_apN0.aperture.density_cm2 = 0.0
ok("continuous aperture N=0 reproduces the unmixed g2 exactly",
   evaluate(d_apN0, [BASE_T])["scalars"]["g2_op"] == g2_legacy)

d_apC0 = _set_thermal(DeviceDesign(), BASE_T)
d_apC0.aperture.compose = True; d_apC0.aperture.comp_brightness = 0.0
ok("continuous aperture c=0 reproduces the unmixed g2 exactly",
   evaluate(d_apC0, [BASE_T])["scalars"]["g2_op"] == g2_legacy)

d_apTiny = _set_thermal(DeviceDesign(), BASE_T)
d_apTiny.aperture.compose = True; d_apTiny.aperture.density_cm2 = 1e-6
r_tiny = evaluate(d_apTiny, [BASE_T])
ok("arbitrarily small N gives a g2 shift that vanishes continuously (no floor)",
   0.0 < abs(r_tiny["scalars"]["g2_op"] - g2_legacy) < 1e-6)

# Independent factorial-moment recomputation of the composed g2_op from the
# device's own reported N_w/eps/rho (redoes the <n(n-1)>/<n>^2 bookkeeping,
# never against device.py's own formula in-line).
d_apOn = _set_thermal(DeviceDesign(), BASE_T); d_apOn.aperture.compose = True
sc_on = evaluate(d_apOn, [BASE_T])["scalars"]
lam = sc_on["aperture_lambda_op"]
# DriveBlock.mu defaults to 0.5 (finite loading), so g2_dot upstream of the
# aperture mix is F1b's cap-2 loading result, not eps directly.
_, P1_on, P2_on = loading_probs(sc_on["mu_resolved"])
denom_on = P1_on + P2_on * (1.0 + sc_on["eps_op"])
g_target = 2.0 * P2_on * sc_on["eps_op"] / denom_on ** 2 if denom_on > 0 else 0.0
g_mix_independent = (g_target + 2.0 * lam + lam * lam) / (1.0 + lam) ** 2
g2_op_independent = g2_from(g_mix_independent, sc_on["rho_op"])
ok("independent factorial-moment recomputation matches the composed g2_op",
   abs(g2_op_independent - sc_on["g2_op"]) < 1e-9)

# Monte-Carlo second method: an exactly brightness-1 target (mean = 1,
# <n(n-1)> = g_target, the two moments _compose_aperture_g2 assumes -- the
# unique 3-outcome distribution {0,1,2} with P0=P2=g_target/2, P1=1-g_target
# achieves both exactly) mixed with an independent Poisson competitor bath
# of mean lam_t.
rng = np.random.default_rng(0)
for g_target, lam_t in ((0.02, 0.25), (0.3, 1.0), (0.05, 3.0)):
    n_pulses = 4_000_000
    P2 = P0 = g_target / 2.0
    P1 = 1.0 - g_target
    target = rng.choice(np.array([0, 1, 2]), size=n_pulses, p=[P0, P1, P2])
    m = target + rng.poisson(lam_t, n_pulses)
    mean = m.mean()
    g2_mc = (m * (m - 1)).mean() / mean ** 2
    g2_formula = float(_compose_aperture_g2(g_target, lam_t))
    ok(f"MC aperture composition matches the closed form (g_target={g_target}, lam={lam_t})",
       abs(g2_mc - g2_formula) < 0.01)


def _nw_for_density(dens):
    d_ = _set_thermal(DeviceDesign(), BASE_T)
    d_.aperture.density_cm2 = dens
    return evaluate(d_, [BASE_T])["scalars"]["N_w"]


def _ap_pen_for(dens):
    d_ = _set_thermal(DeviceDesign(), BASE_T)
    d_.aperture.density_cm2 = dens
    return evaluate(d_, [BASE_T])["scalars"]["aperture_g2_penalty"]


def _g2_op_for(dens, compose):
    d_ = _set_thermal(DeviceDesign(), BASE_T)
    d_.aperture.density_cm2 = dens
    d_.aperture.compose = compose
    return evaluate(d_, [BASE_T])["scalars"]["g2_op"]


# Density continuity across the OLD integer-rounding threshold (Nw = 0.5):
# the legacy informational aperture_g2_penalty scalar (never wired into g2,
# unchanged by this spec) still jumps discretely there; the NEW composed
# g2_op (compose=True) must not. N_w is exactly linear in density_cm2 (all
# other evaluate() factors fixed at a given T_hs), so the Nw=0.5 threshold
# density is found analytically (bisection to machine precision would instead
# straddle the SAME float64 value of Nw on both sides -- not a meaningful
# separation across the round()=0/1 boundary) and the two probe densities are
# placed a comfortable 0.1% to either side of it.
nw_probe, dens_probe = _nw_for_density(1e8), 1e8
dens_at_half = dens_probe * 0.5 / nw_probe
lo, hi = dens_at_half * 0.999, dens_at_half * 1.001
assert _nw_for_density(lo) < 0.5 < _nw_for_density(hi), "Nw=0.5 threshold bracketing failed"
legacy_pen_jump = abs(_ap_pen_for(hi) - _ap_pen_for(lo))
continuous_jump = abs(_g2_op_for(hi, True) - _g2_op_for(lo, True))
ok("legacy informational aperture_g2_penalty still jumps at the Nw=0.5 rounding threshold",
   legacy_pen_jump > 1e-4)
# The 0.2%-in-density probe separation still carries the ordinary smooth
# slope of the continuous composition (a few e-4 in g2_op here) -- the
# absence of a DISCONTINUITY is what matters, so compare against the legacy
# integer-rounding jump's scale (aperture_g2 of a whole extra competitor,
# ~0.3+), not against zero.
ok("continuous composition removes the Nw=0.5 discontinuity from the composed g2_op "
   "(only the ordinary smooth slope remains)",
   continuous_jump < 0.01 and continuous_jump < 0.1 * legacy_pen_jump)

# Losses composed exactly once: an independent recomputation of the full
# op-point chain (spectral eps -> loading g2_dot -> aperture mix ->
# background law) from the exposed scalars alone reproduces g2_op bit-for-bit
# for a loaded (mu>0) design.
staged = DeviceDesign.load(ROOT / "cards" / "staged-device-design.yaml")
staged.thermal.T_hs = 250.0; staged.aperture.compose = True
sc_s = evaluate(staged, [250.])["scalars"]
_, P1, P2 = loading_probs(sc_s["mu_resolved"])
denom = P1 + P2 * (1.0 + sc_s["eps_op"])
g2_dot_recomputed = 2.0 * P2 * sc_s["eps_op"] / denom ** 2 if denom > 0 else 0.0
g2_dot_mixed = _compose_aperture_g2(g2_dot_recomputed, sc_s["aperture_lambda_op"])
g2_op_recomputed = g2_from(g2_dot_mixed, sc_s["rho_op"])
ok("filter/loading/aperture/background losses are each applied exactly once",
   abs(g2_op_recomputed - sc_s["g2_op"]) < 1e-9)

# ============================================================ 9. stack-based material tracking
d_track = DeviceDesign(); d_track.cavity.enabled = True; d_track.ret.mode = "confinement"
d_track.ret.system = EDGE_SYSTEM; d_track.filter.track_material = "dot"
d_track.thermal.T_hs = 300.
r_track = evaluate(d_track, [300.])
d_legacy_track = DeviceDesign(); d_legacy_track.cavity.enabled = True
d_legacy_track.ret.mode = "confinement"; d_legacy_track.ret.system = EDGE_SYSTEM
d_legacy_track.thermal.T_hs = 300.
r_legacy_track = evaluate(d_legacy_track, [300.])
ok("stack-based tracking (InP dot) numerically differs from legacy GaAs tracking",
   r_track["scalars"]["eps_op"] != r_legacy_track["scalars"]["eps_op"])

Tj_track = r_track["scalars"]["T_j_op"]
T_track_target = d_track.cavity.T_track
InP_mat = materials.binary("InP")
mat_resolved = _tracked_material(d_track.ret, "dot")
ok("filter.track_material='dot' resolves to the stack's own dot material (InP)",
   mat_resolved.label == "InP")
dx_stack = (1e3 * (materials.bandgap(InP_mat, Tj_track) - materials.bandgap(InP_mat, T_track_target))
           - d_track.cavity.dEdT_cav * (Tj_track - T_track_target))
dx_legacy_gaas = 1e3 * float(tracking_detuning(Tj_track, d_track.cavity.E_X0, T_track_target,
                                               "GaAs", d_track.cavity.dEdT_cav * 1e-3))
ok("stack-based bandgap-shift anchor: materials.bandgap(InP,T) reproduces the tracking offset "
   "and demonstrably avoids the legacy GaAs tracking function",
   abs(dx_stack) > 0 and abs(dx_stack - dx_legacy_gaas) > 1.0)

d_bad_track = DeviceDesign(); d_bad_track.filter.track_material = "bogus"
raises("reject unknown filter.track_material", lambda: evaluate(d_bad_track, [300.]))

d_hold_track = DeviceDesign(); d_hold_track.cavity.enabled = True; d_hold_track.ret.mode = "confinement"
d_hold_track.ret.system = EDGE_SYSTEM; d_hold_track.filter.track_material = "matrix"
d_hold_track.filter.track = "hold"; d_hold_track.thermal.T_hs = 300.
r_hold_track = evaluate(d_hold_track, [300.])
ok("stack-based tracking composes with filter.track='hold' semantics without crashing",
   np.isfinite(r_hold_track["scalars"]["eps_op"]))

# ===================================================================== 10. CW diagnostics
# (a) Zero-width IRF / no-background limits (module-level anchors, cw_g2's
# own public API, independent of device.py's wiring).
r_ns0, gX0, gXX0, kX0, kXX0, tX0, tXX0 = 1.0, 2.0, 4.0, 0.3, 0.6, 1.0, 0.15
report0 = cw_g2.cw_report(r_ns0, gX0, gXX0, kX0, kXX0, tX0, tXX0, rho=1.0, irf_fwhm_ps=0.0)
ok("zero-width IRF leaves g2_raw0 identical to g2_meas0",
   abs(report0["g2_raw0"] - report0["g2_meas0"]) < 1e-9)
ok("rho=1 (no background) recovers the exact closed-form intrinsic g2_dot0",
   abs(report0["g2_meas0"] - report0["g2_dot0"]) < 1e-12)

A_dip, tau_d_ns, fwhm_ps = 0.7, 1.2, 500.0
sigma_ns = cw_g2.irf_width(fwhm_ps, "gaussian")
tau_grid = np.linspace(-15 * max(tau_d_ns, sigma_ns), 15 * max(tau_d_ns, sigma_ns), 8001)
g_dip = 1.0 - A_dip * np.exp(-np.abs(tau_grid) / tau_d_ns)
g_conv_numeric = cw_g2.convolve_irf(tau_grid, g_dip, fwhm_ps, "gaussian")
g_raw0_numeric = float(np.interp(0.0, tau_grid, g_conv_numeric))
g_raw0_closed = cw_g2.dip_convolved(A_dip, tau_d_ns, sigma_ns)
ok("finite-width single-exponential dip: numeric IRF convolution matches the hand-derived closed form",
   abs(g_raw0_numeric - g_raw0_closed) < 1e-3)

# (b) device.py wiring: independent reconstruction from exposed scalars, and
# no double-quenching (rho_cw is not the pulsed retention S applied again).
d_cw = DeviceDesign(); d_cw.thermal.T_hs = 200.; d_cw.drive.mode = "EL-transport"
d_cw.drive.diode = {"preset": "hkust"}; d_cw.drive.n_dot_cm2 = 1e10; d_cw.drive.cw = True
d_cw.ret.mode = "confinement"; d_cw.ret.preset = "InP/GaAsP0.4/AlGaAs0.4 on GaAs"
sc_cw = evaluate(d_cw, [200.])["scalars"]
Tj_cw = sc_cw["T_j_op"]
p_cw = dot_levels.retention_params(
    dot_levels.levels(dot_levels.class_presets()[d_cw.ret.preset](T=Tj_cw)),
    d_cw.ret.tau_rad_ns, d_cw.ret.channel, verbose=False)
gamma_X_ns = 1.0 / d_cw.ret.tau_rad_ns
k_X, k_XX = cw_g2.escape_rates_from_retention(gamma_X_ns, p_cw["a_esc"], p_cw["E_a"],
                                              p_cw["b_p"], p_cw["E_b"], Tj_cw)
r_ns = sc_cw["loading.r_dot"] / 1e9
t_X, t_XX = sc_cw["t_x_op"], sc_cw["eps_op"] * sc_cw["t_x_op"]
I_X, I_XX = cw_g2.photon_rates(r_ns, gamma_X_ns, 2.0 * gamma_X_ns, k_X, k_XX, 1.0)
sig = t_X * I_X + t_XX * I_XX
# item 3 fix (council review 2026-09-05): the absolute in-window background
# rate is transport's OWN accounting (Background.rate_bg_window, photons/s
# -- the same physical operating point regardless of pulsed/CW framing);
# forming b_e_cw = abs_bg_in_window / I_X (the SAME CW rate-equation X rate
# that feeds sig) and then multiplying back by t_X * I_X cancels I_X
# exactly, leaving bg = abs_bg_in_window * t_X -- one consistent ratio,
# never the pulsed model's own rate_x (b_e_resolved's denominator) chained
# onto this I_X (that was the bug: two different X-rate models mixed
# across the multiplication, 13x too large at a saturated operating point).
abs_bg_in_window_ns = sc_cw["background.rate_bg_window"] / 1e9
b_e_cw = abs_bg_in_window_ns / I_X
bg = b_e_cw * t_X * I_X
rho_cw_expected = sig / (sig + bg)
ok("reconstruction: abs_bg_in_window == b_e * X_rate_collected, pulsed path",
   abs(sc_cw["background.rate_bg_window"] - sc_cw["b_e_resolved"] * sc_cw["background.rate_x"])
   < 1e-6 * max(1.0, sc_cw["background.rate_bg_window"]))
ok("reconstruction: abs_bg_in_window == b_e * X_rate_collected, CW path",
   abs(abs_bg_in_window_ns - b_e_cw * I_X) < 1e-9 * max(1.0, abs_bg_in_window_ns))
report_cw = cw_g2.cw_report(r_ns, gamma_X_ns, 2.0 * gamma_X_ns, k_X, k_XX, t_X, t_XX,
                            rho_cw_expected, d_cw.drive.cw_irf_fwhm_ps,
                            tau_max_ns=d_cw.drive.cw_tau_max_ns,
                            pump_ratio=d_cw.drive.cw_pump_ratio, irf_shape=d_cw.drive.cw_irf_shape)
ok("g2_cw0/g2_cw0_raw match an independent cw_report reconstruction from exposed scalars",
   abs(report_cw["g2_meas0"] - sc_cw["g2_cw0"]) < 1e-9
   and abs(report_cw["g2_raw0"] - sc_cw["g2_cw0_raw"]) < 1e-9)
ok("CW rho is derived from total X+XX signal, not the pulsed retention S applied a second time",
   abs(rho_cw_expected - sc_cw["cw_rho_op"]) < 1e-9 and sc_cw["cw_rho_op"] != sc_cw["S_resolved"])

# aperture-composed CW: same rates, mixed pointwise before background/IRF.
d_cw_ap = copy.deepcopy(d_cw); d_cw_ap.aperture.compose = True
sc_cw_ap = evaluate(d_cw_ap, [200.])["scalars"]
lam_cw = sc_cw_ap["aperture_lambda_op"]
tau_cw = report_cw["curves"]["tau"]
g_mix_cw = _compose_aperture_g2(report_cw["curves"]["g2_dot"], lam_cw)
g_meas_cw = cw_g2.g2_with_background(g_mix_cw, rho_cw_expected)
g_raw_cw = cw_g2.convolve_irf(tau_cw, g_meas_cw, d_cw.drive.cw_irf_fwhm_ps, d_cw.drive.cw_irf_shape)
g2_cw0_expected = float(np.interp(0.0, tau_cw, g_meas_cw))
g2_cw0_raw_expected = float(np.interp(0.0, tau_cw, g_raw_cw))
ok("CW generalizes the same continuous aperture composition, pointwise before background/IRF",
   abs(g2_cw0_expected - sc_cw_ap["g2_cw0"]) < 1e-9
   and abs(g2_cw0_raw_expected - sc_cw_ap["g2_cw0_raw"]) < 1e-9)

# (c) T-aligned curves, mode/domain guards, invalid-rate/zero-signal cases.
r_cw_curve = evaluate(d_cw, np.linspace(150., 350., 5))
ok("drive.cw exposes T-aligned g2_cw0/g2_cw0_raw curves",
   len(r_cw_curve["curves"]["g2_cw0"]) == 5 and len(r_cw_curve["curves"]["g2_cw0_raw"]) == 5)

d_cw_bad_mode = DeviceDesign(); d_cw_bad_mode.drive.cw = True
raises("drive.cw=True requires drive.mode='EL-transport'", lambda: evaluate(d_cw_bad_mode, [300.]))

d_cw_zero = DeviceDesign(); d_cw_zero.drive.mode = "EL-transport"
d_cw_zero.drive.diode = {"preset": "hkust"}; d_cw_zero.drive.I_uA = 0.0; d_cw_zero.drive.cw = True
r_cw_zero = evaluate(d_cw_zero, [d_cw_zero.thermal.T_hs])
ok("zero-current CW is ineligible (nan g2_cw0), not a crash",
   np.isnan(r_cw_zero["scalars"]["g2_cw0"]))

# ============================================================ 11. new-block round-trip
d_rt2 = DeviceDesign()
d_rt2.emission.type = "edge"; d_rt2.emission.lambda_nm = 668.0
d_rt2.emission.ridge_width_nm = 3000.0; d_rt2.emission.R_back = 0.9
d_rt2.aperture.compose = True
d_rt2.filter.track_material = "matrix"
d_rt2.drive.cw = True; d_rt2.drive.cw_irf_fwhm_ps = 300.0
with tempfile.NamedTemporaryFile(dir=ROOT, suffix=".yaml", delete=False) as tmp2:
    path2 = Path(tmp2.name)
try:
    d_rt2.save(path2); loaded2 = DeviceDesign.load(path2)
finally:
    os.unlink(path2)
ok("YAML preserves emission/aperture.compose/track_material/drive.cw fields", loaded2 == d_rt2)
ok("legacy card YAML (no 'emission' key) loads with the default EmissionBlock",
   DeviceDesign.load(ROOT / "cards" / "staged-device-design.yaml").emission == EmissionBlock())

env_new = evaluate_envelope(d_edge, {"emission.ridge_width_nm": (1500.0, 2500.0)}, T_grid=[300.])
ok("evaluate_envelope sweeps a new-block field without error and preserves Lemma 1",
   env_new["n_samples"] == 2
   and np.isfinite(env_new["scalar_bands"]["g2_op"]).all()
   and env_new["scalar_bands"]["g2_op"][0] == env_new["scalar_bands"]["g2_op"][1] == sc_edge["g2_op"])

# ============================================== 12. council review 2026-09-05 fixes
# Item 1: carrier conservation at the gainp card's own EL-transport +
# confinement operating point -- a single dot cannot receive more than the
# supply-limited flux f_QD*eta_inj*(I/q) (the reviewer's own bound, 0.53).
d_gainp = DeviceDesign.load(ROOT / "cards" / "edge-inp-gainp-design.yaml")
sc_gainp = evaluate(d_gainp)["scalars"]
Tj_gainp = sc_gainp["T_j_op"]
diode_gainp = _diode_from_drive(d_gainp.drive)
n_dot_cm2_gainp = d_gainp.drive.n_dot_cm2 or d_gainp.aperture.density_cm2
lk_gainp = diode_gainp.eta_inj(Tj_gainp)
f_QD_gainp = diode_gainp.f_qd(n_dot_cm2_gainp)
supply_bound_gainp = f_QD_gainp * lk_gainp.eta_inj * (d_gainp.drive.I_uA * 1e-6 / transport.Q_SI)
ok("gainp card: mu stays at or below the reviewer's supply-limited bound (0.53)",
   sc_gainp["mu_resolved"] <= 0.53 + 1e-9)
ok("gainp card: loading.r_dot never exceeds f_QD*eta_inj*(I/q) (carrier conservation)",
   sc_gainp["loading.r_dot"] <= supply_bound_gainp * (1.0 + 1e-6))

# Item 3: CW rho must match pulsed rho within 2% in the linear (unsaturated)
# pump regime -- both derive from the SAME transport background/X-rate
# accounting once the ratio is formed with a single, consistent denominator
# (device.py bg = injection.background.rate_bg_window * t_X, item 3 fix).
d_linear = copy.deepcopy(d_gainp)
d_linear.drive.I_uA = 1e-6  # deep linear regime: mu << 1, well off dot saturation
sc_linear = evaluate(d_linear)["scalars"]
ok("CW rho matches pulsed rho within 2% in the linear pump regime",
   sc_linear["mu_resolved"] < 0.01
   and abs(sc_linear["cw_rho_op"] - sc_linear["rho_op"]) < 0.02 * sc_linear["rho_op"])

# Item 4: filter.track_material must no longer be inert with cavity.enabled
# False (both edge-emitter cards run cavity-less) -- the stack-derived
# Varshni walk of the filter window must move g2_op at a T away from
# cavity.T_track, where the InP and GaAs walks differ by more than the
# window-width tolerance.
d_notrack = DeviceDesign(); d_notrack.ret.mode = "confinement"; d_notrack.ret.system = EDGE_SYSTEM
d_notrack.thermal.T_hs = 400.0
r_notrack = evaluate(d_notrack, [400.0])
d_track4 = copy.deepcopy(d_notrack); d_track4.filter.track_material = "dot"
r_track4 = evaluate(d_track4, [400.0])
ok("filter.track_material is no longer inert with cavity.enabled=False (g2_op moves)",
   r_notrack["scalars"]["g2_op"] != r_track4["scalars"]["g2_op"])

print(f"{sum(checks)}/{len(checks)} device RT checks passed")
sys.exit(0 if all(checks) else 1)
