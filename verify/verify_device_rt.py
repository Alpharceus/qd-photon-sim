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
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import fsim_core.device as devmod
from fsim_core import lindblad as lb
from fsim_core.device import (DeviceDesign, EmissionBlock, evaluate, evaluate_envelope,
                              _stack, _resolve_edge, _compose_aperture_g2, _aperture_lambda,
                              _tracked_material, _diode_from_drive, _confinement_params)
from fsim_core.linewidth import LinewidthParams, gamma_anchor
from fsim_core.integrator import retention, g2_from
from fsim_core.loading import loading_probs
from fsim_core.thermal import t_junction
from fsim_core.cavity import tracking_detuning
from fsim_core.spectral import epsilon
from fsim_core import dot_levels, transport, waveguide, cw_g2, materials, pulse_counting

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


def _cw_reduction_guard(design, label):
    """pr-pkg4-fix3 item 3, corrected by pr-pkg4-fix4 item 1 (Opus review):
    device.py's rho_pulsed reduces to cw_rho_op in the tau_dark_ns -> 0
    limit ONLY when the RESOLVED retention params evaluate() actually used
    -- scalars["retention_params_used"]["b0"]/["beta"] -- are both 0.0 and
    the cavity is disabled -- cw_rho_op's own background never carries a
    b0/beta term and its own signal never carries a cavity gain factor, so
    either active makes the two rho definitions differ BY CONSTRUCTION
    (see device.py's finite_pulse block, item 3).

    The raw design.ret.b0/design.ret.beta fields the previous version of
    this guard tested are NOT the right condition: in ret.mode="proxy"
    (device.py's default) an unset (0.0) ret.b0/ret.beta resolves through
    the class-proxy Arrhenius fit to a nonzero b0/beta, so a proxy-mode
    card can fail the reduction while its own raw fields still read 0.0 --
    the old guard would have wrongly let such a card through checks E/F.
    In ret.mode="confinement" an explicit ret.overrides entry can likewise
    set a nonzero b0/beta on top of raw fields that stay 0.0, breaking the
    reduction the same way. Skip with an explanatory message rather than
    asserting a false equivalence if a future card under test resolves
    either to nonzero."""
    # pr-final-nits item 3: restrict the sweep to the single T_hs point the
    # guard actually reads (identical result to the full default T_grid,
    # since retention_params_used at T_hs is unaffected by which other
    # points are also evaluated) -- 59 s -> ~1 s per call.
    rpu = evaluate(design, T_grid=[design.thermal.T_hs])["scalars"]["retention_params_used"]
    b0_used = rpu["b0"] if rpu is not None else float("nan")
    beta_used = rpu["beta"] if rpu is not None else float("nan")
    if not (b0_used == 0.0 and beta_used == 0.0) or design.cavity.enabled:
        print(f"skip {label}: resolved retention params b0={b0_used!r}/"
              f"beta={beta_used!r} are not both 0.0, or the cavity is "
              "enabled -- rho_pulsed and cw_rho_op differ by construction "
              "(pr-pkg4-fix3 item 3), not comparable here")
        return True
    return False


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

# pr-pkg1-fix2 item 1 (peer-review-triage.md finding 1b): legacy invariance.
# A confinement design with no explicit aperture density set (the design
# above: ApertureBlock.density_cm2 defaults to None) must still resolve
# dot_levels.retention_params' own 1e10 cm^-2 default -- ApertureBlock.
# density_cm2's OLD literal 7.0e8 default silently fed _confinement_params
# too, so no legacy design could ever reach that 1e10 default; this is the
# exact regression that made the check above fail before the fix.
ok("ApertureBlock.density_cm2 defaults to None (no aperture density set)",
   d.aperture.density_cm2 is None)
ok("legacy confinement design (no explicit aperture density) reproduces "
   "S(300K) = 9.63774e-4 (dot_levels.retention_params' 1e10 cm^-2 default, "
   "not a silently substituted aperture density)",
   abs(r["scalars"]["S_resolved"] / 9.63774e-4 - 1.) < 1e-6)

# pr-pkg1-fix3 item 2: drive.n_dot_cm2 is TRANSPORT-only and must never
# reach confinement. The same legacy design above (no explicit aperture
# density), but now with drive.n_dot_cm2=1e9 set (as a legacy card that
# configures only the transport-side density would), must resolve the
# EXACT SAME S(300K) = 9.637738e-4 as the drive.n_dot_cm2-unset case above
# -- the pre-fix `d.drive.n_dot_cm2 or d.aperture.density_cm2` call site
# forwarded drive.n_dot_cm2 into confinement instead, moving this design's
# S to ~1.021141e-4 (dot_levels.retention_params' states_per_dot then used
# 1e9 instead of falling to its own 1e10 default).
d_drive_only = DeviceDesign(); d_drive_only.thermal.T_hs = 300.; d_drive_only.drive.V = 0.
d_drive_only.ret.mode = "confinement"
d_drive_only.ret.preset = "InP/GaAsP0.4/AlGaAs0.4 on GaAs"
d_drive_only.drive.n_dot_cm2 = 1e9
r_drive_only = evaluate(d_drive_only, [300.])
ok("legacy design with only drive.n_dot_cm2=1e9 set (no aperture density) "
   "still resolves confinement's S(300K) = 9.637738e-4 rel 1e-6 -- "
   "drive.n_dot_cm2 is never forwarded to confinement",
   d_drive_only.aperture.density_cm2 is None
   and abs(r_drive_only["scalars"]["S_resolved"] / 9.637738e-4 - 1.) < 1e-6)

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
# council review 2026-09-05 item 1: EL-transport pulsed drive (drive.cw is
# False, the default) now REQUIRES an explicit drive.diode['tau_pulse_ns'] and
# either drive.duty or drive.rep_rate_hz -- every fixture below states both
# explicitly instead of relying on the (now-rejected) untouched defaults.
d = DeviceDesign(); d.drive.mode = "EL-transport"
d.drive.diode = {"preset": "hkust", "tau_pulse_ns": 0.1}; d.drive.duty = 0.008
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
d0 = DeviceDesign(); d0.drive.mode = "EL-transport"
d0.drive.diode = {"preset": "hkust", "tau_pulse_ns": 0.1}; d0.drive.duty = 0.008
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
d = DeviceDesign(); d.drive.mode = "EL-transport"
d.drive.diode = {"preset": "hkust", "tau_pulse_ns": 0.1}; d.drive.duty = 0.008
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
d2.drive.mode = "EL-transport"
d2.drive.diode = {"preset": "red", "R_s_ohm": 7., "tau_pulse_ns": 0.1}; d2.drive.duty = 0.008
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
# directly (mirrors verify_waveguide.py's ridge beta class check
# and its own facet_transmission()/beta_factor() checks -- never against
# device.py's own arithmetic in-line).
edge_direct, lam_direct = _resolve_edge(d_edge.ret, d_edge.emission, sc_edge["T_j_op"])
# Audit H4 (2026-09-22):
# F_wg now uses the energy-normalised Purcell area (Lecamp, Lalanne &
# Hugonin, PRL 99, 023902 (2007)), which for a Gaussian is exactly half the
# old nonlinear (int I)^2/int I^2 area [DR] (this fixture: beta 0.0261 ->
# 0.0454). Spec audit-edge-cards-label task 2 (2026-09-23): the bound is the
# ONE shared [E] envelope verify/data/rt_edge_anchors.yaml
# envelopes.ridge_guided_beta (0.0060, 0.067; derivation there, re-derived
# in verify_waveguide.py), replacing this file's separate (0.006, 0.06).
_BETA_ENV = yaml.safe_load(
    (ROOT / "verify" / "data" / "rt_edge_anchors.yaml").read_text(encoding="utf-8")
)["envelopes"]["ridge_guided_beta"]
ok(f"2 um ridge fixture beta sits in the shared ({_BETA_ENV['lo']}, {_BETA_ENV['hi']}) "
   "ridge class envelope [E] (envelopes.ridge_guided_beta, Purcell-area rescaled)",
   _BETA_ENV["lo"] < edge_direct.beta < _BETA_ENV["hi"])
ok("device.py edge scalars match a fresh _resolve_edge call bit-for-bit",
   edge_direct.eta_total == sc_edge["edge_eta_total"] and lam_direct == sc_edge["edge_lambda_nm"])
T_recomputed = waveguide.facet_transmission(sc_edge["edge_n_eff"])
ok("T_facet matches an independent facet_transmission(n_eff) recomputation",
   abs(T_recomputed - sc_edge["edge_T_facet"]) < 1e-9)
# Facet convention (peer-review fix, 2026-09-07): the old structural identity
# eta_total == beta * 0.5 * T_facet * eta_prop * eta_NA is obsolete --
# propagation is now folded entirely into the facet ray-series (eta_facet),
# never applied again as a separate eta_prop factor. The new decomposition
# is eta_total == beta * eta_facet(...) * eta_NA, with eta_facet recomputed
# HERE by hand from the row's own T_facet, R_back, alpha_cm, L_um (never by
# calling facet_escape_fraction or edge_emission a second time -- see
# fsim_core/waveguide.py's facet_escape_fraction docstring for the same
# closed form). d_edge.emission carries no coating override, so
# R_back=None resolves (item 4) to the uncoated Fresnel value, which here
# equals 1-T_facet exactly (same as the pre-item-4 R_front resolution).
ok("edge_T_facet is not unity (so the item-2 check below is not vacuous)",
   sc_edge["edge_T_facet"] < 1.0)
_T_dr = sc_edge["edge_T_facet"]
_alpha_dr = d_edge.emission.alpha_cm
_L_dr = d_edge.emission.L_um
_Rback_dr = d_edge.emission.R_back
ok(f"d_edge.emission carries no coating override (coating={d_edge.emission.coating!r}) "
   "before applying the R_back=None -> 1-T_facet substitution below",
   not d_edge.emission.coating)
if _Rback_dr is None:
    _Rback_dr = 1.0 - _T_dr  # no coating override on this stack (item 4)
_a_dr = _alpha_dr * 1e-4
_x_dr = 0.5  # facet_escape_fraction's dot_position default
_Rfront_dr = 1.0 - _T_dr
_prop_rt_dr = np.exp(-2.0 * _a_dr * _L_dr)
_eta_facet_dr = (0.5 * _T_dr * np.exp(-_a_dr * _x_dr * _L_dr)
                * (1.0 + _Rback_dr * np.exp(-2.0 * _a_dr * (1.0 - _x_dr) * _L_dr))
                / (1.0 - _Rback_dr * _Rfront_dr * _prop_rt_dr))
eta_recomputed = sc_edge["edge_beta"] * _eta_facet_dr * sc_edge["edge_eta_NA"]
ok("eta_total = beta * eta_facet(T_facet, R_back, alpha_cm, L_um) * eta_NA, each exactly once",
   abs(eta_recomputed - sc_edge["edge_eta_total"]) < 1e-9)

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
#
# Item 8 (council review 2026-09-05): EDGE_SYSTEM's own natural (3/10 nm dot)
# confinement wavelength (~788 nm) now falls INSIDE MATERIAL_EXTRA's 700-850
# nm Ga0.51In0.49P/(AlGa)InP table added since this check was written, so it
# is no longer a genuine "missing evidence" case. A larger dot (lower
# confinement energy -> longer wavelength) is retargeted just past the
# table's 850 nm edge (~860 nm) instead, so the check still exercises the
# ineligible-not-crash path against data that is genuinely absent.
MISSING_SYSTEM = dict(EDGE_SYSTEM, geometry={"height_nm": 10.0, "radius_nm": 24.0})
d_missing = DeviceDesign(); d_missing.thermal.T_hs = 300.; d_missing.drive.V = 0.
d_missing.ret.mode = "confinement"; d_missing.ret.system = MISSING_SYSTEM
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
# Re-pinned (audit 2026-09-23, device.py:2128): target counts s, Poisson
# bath lam*s and uncorrelated background b are three independent sources, so
# the exact factorial-moment result is 1 - s^2 (1-g_t)/((1+lam)s + b)^2 --
# written with s = rho_op, b = 1 - rho_op (rho_op is the target-only signal
# fraction). The pre-fix g2_from(g_mix, rho_op) double-counted the bath x
# background cross term.
def _three_source_g2(g_t, lam_, s_, b_):
    return 1.0 - s_ * s_ * (1.0 - g_t) / ((1.0 + lam_) * s_ + b_) ** 2
g2_op_independent = _three_source_g2(g_target, lam, sc_on["rho_op"], 1.0 - sc_on["rho_op"])
ok("independent factorial-moment recomputation matches the composed g2_op "
   "(exact three-source form, audit 2026-09-23)",
   abs(g2_op_independent - sc_on["g2_op"]) < 1e-9)

# Audit 2026-09-23 acceptance 2: device.py's aperture composition
# (_compose_aperture_g2 mix, then g2_from at rho_eff = (1+lam)rho/(1+lam rho))
# equals the closed form 1 - s^2(1-g_t)/((1+lam)s+b)^2 on synthetic cases to
# 1e-12, and the exposed rho_g2_op is that rho_eff at the op point.
_ap_synth_ok = True
for _gt, _lam, _s, _b in ((0.02, 0.0339, 0.75, 0.25), (0.3, 1.0, 0.5, 0.5),
                          (0.0, 0.25, 2.0e-4, 7.0e-5), (0.6, 3.0, 1.0, 0.0)):
    _rho = _s / (_s + _b)
    _dev = g2_from(_compose_aperture_g2(_gt, _lam), devmod._aperture_rho_eff(_rho, _lam))
    _ap_synth_ok &= abs(_dev - _three_source_g2(_gt, _lam, _s, _b)) < 1e-12
_rho_eff_op = (1.0 + lam) * sc_on["rho_op"] / ((1.0 + lam) * sc_on["rho_op"] + 1.0 - sc_on["rho_op"])
ok("aperture composition equals the closed form 1 - s^2(1-g_t)/((1+lam)s+b)^2 "
   "on synthetic cases to 1e-12 (audit 2026-09-23); rho_g2_op == (1+lam)s/((1+lam)s+b)",
   _ap_synth_ok and abs(sc_on["rho_g2_op"] - _rho_eff_op) < 1e-12)

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
# re-pinned (audit 2026-09-23): exact three-source background law.
g2_op_recomputed = _three_source_g2(g2_dot_recomputed, sc_s["aperture_lambda_op"],
                                    sc_s["rho_op"], 1.0 - sc_s["rho_op"])
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
d_cw.drive.diode = {"preset": "hkust"}; d_cw.drive.n_dot_cm2 = 1e9; d_cw.drive.cw = True
d_cw.ret.mode = "confinement"; d_cw.ret.preset = "InP/GaAsP0.4/AlGaAs0.4 on GaAs"
sc_cw = evaluate(d_cw, [200.])["scalars"]
Tj_cw = sc_cw["T_j_op"]
# pr-pkg1-fix2 item 4: device.py._confinement_params now forwards
# T_ref=Tj (not retention_params' own 300 K default) and the card's
# resolved n_dot_cm2 (not retention_params' own 1e10 default) -- this
# independent reconstruction must match, or it is no longer reconstructing
# what evaluate() actually computed.
p_cw = dot_levels.retention_params(
    dot_levels.levels(dot_levels.class_presets()[d_cw.ret.preset](T=Tj_cw)),
    d_cw.ret.tau_rad_ns, d_cw.ret.channel, T_ref=Tj_cw,
    # pr-pkg1-fix4 item 2: the resolution device.py actually uses is
    # ret.n_dot_cm2 (explicit) else aperture.density_cm2 (explicit) else
    # retention_params' own 1e10 default -- drive.n_dot_cm2 is
    # TRANSPORT-only and must never reach confinement (see
    # _confinement_params' own docstring); the old `d_cw.drive.n_dot_cm2
    # or ...` chain leaked it and only happened to match here because the
    # d_cw fixture set drive.n_dot_cm2 to the same 1e10 default.
    n_dot_cm2=(d_cw.ret.n_dot_cm2 if d_cw.ret.n_dot_cm2 is not None else (
        d_cw.aperture.density_cm2 if d_cw.aperture.density_cm2 is not None else 1e10)),
    verbose=False)
gamma_X_ns = 1.0 / d_cw.ret.tau_rad_ns
k_X, k_XX = cw_g2.escape_rates_from_retention(gamma_X_ns, p_cw["a_esc"], p_cw["E_a"],
                                              p_cw["b_p"], p_cw["E_b"], Tj_cw)
r_ns = sc_cw["loading.r_dot"] / 1e9
t_X, t_XX = sc_cw["t_x_op"], sc_cw["eps_op"] * sc_cw["t_x_op"]
I_X, I_XX = cw_g2.photon_rates(r_ns, gamma_X_ns, 2.0 * gamma_X_ns, k_X, k_XX, 1.0)
sig = t_X * I_X + t_XX * I_XX
# item 5 fix (round 2, council review 2026-09-05): the absolute in-window
# background rate is transport's OWN accounting (Background.rate_bg_window,
# photons/s -- the same physical operating point regardless of pulsed/CW
# framing); a flat continuum's collected rate scales with the collection
# (filter) window width, never with the X LINE's own transmission t_X
# (round 1's `abs_bg_in_window_ns * t_X`, a ~940x undercount at the sweep's
# worst corner) -- d_cw's design carries no drive.diode.w_meV override, so
# the filter window equals the transport window and the scale is 1 exactly.
abs_bg_in_window_ns = sc_cw["background.rate_bg_window"] / 1e9
bg = abs_bg_in_window_ns
rho_cw_expected = sig / (sig + bg)

# Item 6 (council review 2026-09-05): the round-1 "reconstruction" checks
# here were tautologies -- they re-derived abs_bg_in_window from b_e_resolved
# and background.rate_x, device.py's OWN pre-computed ratio, so they could
# not have caught a wrong rate_bg formula. Reconstruct rate_bg INDEPENDENTLY
# instead, straight from transport's public primitives -- (1-f_QD) eta_inj
# I/q eta_rad_matrix xi PLUS the routed (1-f_qfl) share (item 2) -- never
# from device.py's own Background object.
diode_cw = _diode_from_drive(d_cw.drive, d_cw.aperture.diameter_um)
n_dot_cw = d_cw.drive.n_dot_cm2 or d_cw.aperture.density_cm2
lk_cw = diode_cw.eta_inj(Tj_cw)
f_QD_cw = diode_cw.f_qd(n_dot_cw)
E_X_cw = _confinement_params(d_cw.ret, Tj_cw)["E_X_eV"]
V_j_cw = diode_cw.vj_of_j(d_cw.drive.I_uA * 1e-6 / diode_cw.area_cm2, Tj_cw)
kT_eV_cw = materials.KB_EV * Tj_cw
supply_cw = lk_cw.eta_inj * (d_cw.drive.I_uA * 1e-6 / transport.Q_SI)
f_qfl_cw = transport.qfl_suppression(E_X_cw, V_j_cw, kT_eV_cw)
r_matrix_cw = (1.0 - f_QD_cw) * supply_cw + f_QD_cw * supply_cw * (1.0 - f_qfl_cw)
w_f_cw = sc_cw["gamma_op"]  # filter.auto_w default: w = Gamma(Tj); no drive.diode.w_meV here
kT_meV_cw = 1e3 * kT_eV_cw
xi_cw = float(transport.xi_window(w_f_cw, 100.0, kT_meV_cw))  # dE_WL_meV default 100.0
f_qfl_bg_cw = transport.qfl_suppression(E_X_cw + 100.0e-3, V_j_cw, kT_eV_cw)
rate_bg_independent_ns = r_matrix_cw * 0.1 * xi_cw * f_qfl_bg_cw / 1e9  # eta_rad_matrix default 0.1
ok("independent reconstruction of the in-window background rate matches device.py's "
   "(from transport primitives directly, not device.py's own Background object)",
   abs(rate_bg_independent_ns - abs_bg_in_window_ns) < 1e-6 * max(1.0, abs_bg_in_window_ns))
ok("CW rho reconstructed from the independently-derived background matches device.py's cw_rho_op",
   abs(rho_cw_expected - sc_cw["cw_rho_op"]) < 1e-9)

# Item 5 (council review 2026-09-05): a flat background continuum must be
# independent of the X line's own spectral transmission t_X (dx) -- moving
# the filter off the line must degrade rho_cw (background stays fixed while
# signal shrinks), which round 1's `bg = abs_bg_in_window * t_X` masked
# (bg shrank in lockstep with the signal, leaving rho roughly dx-insensitive).
d_cw_dx = copy.deepcopy(d_cw); d_cw_dx.filter.dx = 10.0
sc_cw_dx = evaluate(d_cw_dx, [200.])["scalars"]
ok("item 5: background.rate_bg_window (flat continuum) is independent of the filter's dx",
   abs(sc_cw_dx["background.rate_bg_window"] - sc_cw["background.rate_bg_window"])
   < 1e-9 * max(1.0, sc_cw["background.rate_bg_window"]))
ok("item 5: moving the filter off the X line degrades CW rho (background no longer "
   "shrinks alongside t_X, unlike round 1) -- dx=10 meV clears the XX-line crossover "
   "region (dx near dot.delta_xx=3.5 briefly INCREASES total X+XX transmission)",
   sc_cw_dx["t_x_op"] < sc_cw["t_x_op"] and sc_cw_dx["cw_rho_op"] < sc_cw["cw_rho_op"])

# Item 5: a flat continuum's collected rate scales as min(1, w_f/w_transport)
# when drive.diode.w_meV gives the transport window a different width than
# the filter's own -- never by t_X. Checked on the simpler pulsed path
# (b_e_resolved); pulsed and CW share the identical win_scale computation.
d_cw_wide = copy.deepcopy(d_cw)
d_cw_wide.drive.diode = dict(d_cw.drive.diode)
d_cw_wide.drive.diode["w_meV"] = w_f_cw * 4.0  # transport window 4x WIDER than the filter
sc_cw_wide = evaluate(d_cw_wide, [200.])["scalars"]
raw_ratio_wide = sc_cw_wide["background.rate_bg_window"] / sc_cw_wide["background.rate_x"]
ok("item 5: b_e_resolved applies min(1, w_f/w_transport) window scaling, not t_X, when the "
   "transport window differs from the filter's own",
   abs(sc_cw_wide["b_e_resolved"] - raw_ratio_wide * 0.25) < 1e-9 * max(1.0, sc_cw_wide["b_e_resolved"]))

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
# re-pinned (audit 2026-09-23): exact three-source law pointwise in tau.
g_meas_cw = _three_source_g2(np.asarray(report_cw["curves"]["g2_dot"]), lam_cw,
                             rho_cw_expected, 1.0 - rho_cw_expected)
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

# (d) pr-pkg1-fix3 item 5: frozen CW-path regression fixture. cw_g2.
# convolve_irf switched from np.convolve to scipy.signal.fftconvolve
# (pr-pkg1-fix2, runtime item 2) for O(N log N) instead of O(N * kernel)
# on a stiff escape rate's tau grid; that is the ONLY path in this file
# that exercises convolve_irf against a design card's own (as-shipped)
# operating point rather than a synthetic fixture, so it is the one place
# that would have caught a future convolve_irf regression on a real card.
# Literals captured 2026-09-07 on this tree (edge-inp-gainp-design.yaml,
# drive.cw=true as shipped, T_hs=300 K as shipped); measured fftconvolve
# vs np.convolve agreement on this exact operating point is 1.111e-16
# relative on g2_cw0_raw (see cw_g2.convolve_irf's docstring) -- float64
# round-off, so 1e-9 relative is a tight but safe regression tolerance,
# nowhere near loose enough to hide a real algorithm change.
d_cw_gainp = DeviceDesign.load(ROOT / "cards" / "edge-inp-gainp-design.yaml")
assert d_cw_gainp.drive.cw, "edge-inp-gainp-design.yaml must ship drive.cw=true for this fixture"
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    sc_cw_gainp = evaluate(d_cw_gainp)["scalars"]
# Re-pinned 2026-09-23 (audit device composition, device.py:2117): the CW
# aperture composition now uses rho_eff = (1+lam)rho/(1+lam rho) (exact
# three-source law); was 0.7624234413164181 / 0.9991637853703624. cw_rho_op
# (target-only) is unchanged.
CW_GAINP_G2_CW0 = 0.7601106076524684
CW_GAINP_G2_CW0_RAW = 0.9991556447299034
CW_GAINP_CW_RHO_OP = 0.8527446400441377
ok("CW-path regression fixture: edge-inp-gainp-design.yaml's own operating "
   "point reproduces g2_cw0 (1e-9 relative)",
   abs(sc_cw_gainp["g2_cw0"] / CW_GAINP_G2_CW0 - 1.0) < 1e-9)
ok("CW-path regression fixture: edge-inp-gainp-design.yaml's own operating "
   "point reproduces g2_cw0_raw (1e-9 relative -- the fftconvolve-carrying "
   "quantity)",
   abs(sc_cw_gainp["g2_cw0_raw"] / CW_GAINP_G2_CW0_RAW - 1.0) < 1e-9)
ok("CW-path regression fixture: edge-inp-gainp-design.yaml's own operating "
   "point reproduces cw_rho_op (1e-9 relative)",
   abs(sc_cw_gainp["cw_rho_op"] / CW_GAINP_CW_RHO_OP - 1.0) < 1e-9)

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

# pr-pkg4-fix item 9: drive.loading_model and drive.gate_ns are now
# validated at DeviceDesign.load() itself, not only once evaluate() runs.
d_load_base = DeviceDesign(); d_load_base.drive.mode = "EL-transport"
d_load_base.drive.diode = {"preset": "hkust"}
with tempfile.NamedTemporaryFile(dir=ROOT, suffix=".yaml", delete=False) as tmp_lm:
    path_lm = Path(tmp_lm.name)
try:
    d_load_base.save(path_lm)
    raw_lm = yaml.safe_load(path_lm.read_text(encoding="utf-8"))
    raw_lm["design"]["drive"]["loading_model"] = "bogus"
    path_lm.write_text(yaml.safe_dump(raw_lm, sort_keys=False), encoding="utf-8")
    raises("item 9: DeviceDesign.load() rejects an unknown drive.loading_model "
          "(not only evaluate())", lambda: DeviceDesign.load(path_lm))
finally:
    os.unlink(path_lm)

with tempfile.NamedTemporaryFile(dir=ROOT, suffix=".yaml", delete=False) as tmp_gn:
    path_gn = Path(tmp_gn.name)
try:
    d_load_base.save(path_gn)
    raw_gn = yaml.safe_load(path_gn.read_text(encoding="utf-8"))
    raw_gn["design"]["drive"]["gate_ns"] = -1.0
    path_gn.write_text(yaml.safe_dump(raw_gn, sort_keys=False), encoding="utf-8")
    raises("item 9: DeviceDesign.load() rejects drive.gate_ns <= 0 (not only evaluate())",
          lambda: DeviceDesign.load(path_gn))
finally:
    os.unlink(path_gn)

env_new = evaluate_envelope(d_edge, {"emission.ridge_width_nm": (1500.0, 2500.0)}, T_grid=[300.])
ok("evaluate_envelope sweeps a new-block field without error and preserves Lemma 1",
   env_new["n_samples"] == 2
   and np.isfinite(env_new["scalar_bands"]["g2_op"]).all()
   and env_new["scalar_bands"]["g2_op"][0] == env_new["scalar_bands"]["g2_op"][1] == sc_edge["g2_op"])

# ============================================== 12. council review 2026-09-05 fixes
# Item 1: carrier conservation at the gainp card's own EL-transport +
# confinement operating point -- a single dot cannot receive more than the
# supply-limited flux f_QD*eta_inj*(I/q).
d_gainp = DeviceDesign.load(ROOT / "cards" / "edge-inp-gainp-design.yaml")
sc_gainp = evaluate(d_gainp)["scalars"]
Tj_gainp = sc_gainp["T_j_op"]
diode_gainp = _diode_from_drive(d_gainp.drive, d_gainp.aperture.diameter_um)
n_dot_cm2_gainp = d_gainp.drive.n_dot_cm2 or d_gainp.aperture.density_cm2
lk_gainp = diode_gainp.eta_inj(Tj_gainp)
f_QD_gainp = diode_gainp.f_qd(n_dot_cm2_gainp)
supply_bound_gainp = f_QD_gainp * lk_gainp.eta_inj * (d_gainp.drive.I_uA * 1e-6 / transport.Q_SI)
ok("gainp card: loading.r_dot never exceeds f_QD*eta_inj*(I/q) (carrier conservation)",
   sc_gainp["loading.r_dot"] <= supply_bound_gainp * (1.0 + 1e-6))

# pr-pkg4-fix4 item 2 (Opus review): regression checks for the hoisted
# density validation (pr-pkg6-fix item 1) -- previously untested. All
# three raise ValueError (not the underlying division's ZeroDivisionError
# -- an uncaught ZeroDivisionError here would crash this script rather
# than being silently missed) on the gainp card, whose ret.mode=
# "confinement" and drive.mode="EL-transport" both satisfy the mode gate.
d_dens0 = copy.deepcopy(d_gainp); d_dens0.aperture.density_cm2 = 0.0
raises("item 2: aperture.density_cm2=0.0 raises ValueError (not "
       "ZeroDivisionError) on the gainp card", lambda: evaluate(d_dens0))
d_densneg = copy.deepcopy(d_gainp); d_densneg.aperture.density_cm2 = -1.0
raises("item 2: aperture.density_cm2<0 raises ValueError (not "
       "ZeroDivisionError) on the gainp card", lambda: evaluate(d_densneg))
d_ndot0 = copy.deepcopy(d_gainp); d_ndot0.ret.n_dot_cm2 = 0.0
raises("item 2: ret.n_dot_cm2=0.0 raises ValueError (not "
       "ZeroDivisionError) on the gainp card", lambda: evaluate(d_ndot0))

# pr-pkg4-fix4 item 2: a NEGATIVE aperture.density_cm2 must raise in EVERY
# mode, unlike the `== 0.0` case above (which stays mode-gated to protect
# aperture.compose's own N=0 use -- see device.py's comment at the
# validation site) -- construct a design where NEITHER mode-gate condition
# holds (ret.mode="proxy", drive.mode="EL", both untouched defaults) and
# confirm the unconditional check still catches the negative density.
d_density_proxy_neg = DeviceDesign()
d_density_proxy_neg.aperture.density_cm2 = -1.0
raises("item 2: NEGATIVE aperture.density_cm2 raises ValueError even in "
       "proxy mode / drive.mode='EL' (neither mode-gate condition holds)",
       lambda: evaluate(d_density_proxy_neg))

# Item 6 (council review 2026-09-05): the round-1 mu bound (0.53) was a
# hard-coded magic number, not computed from the card -- mu = r_dot *
# tau_pulse, and r_dot is bounded by the (per-dot, N_eff-capped) supply, so
# mu's own supply bound is f_QD * eta_inj * (I/q) * tau_pulse at the card's
# current (never tighter than dividing by N_eff first, so this is always a
# valid, if sometimes loose, bound).
tau_pulse_gainp = float(d_gainp.drive.diode.get("tau_pulse_ns", 1.0))
mu_bound_gainp = supply_bound_gainp * tau_pulse_gainp * 1e-9
ok("gainp card: mu stays at or below the supply bound f_QD*eta_inj*(I/q)*tau_pulse",
   sc_gainp["mu_resolved"] <= mu_bound_gainp * (1.0 + 1e-9))

# Item 6 (council review 2026-09-05): the round-1 CW/pulsed rho comparison
# ran at I=1e-6 uA, deep enough sub-turn-on that BOTH rho values were ~5e-7
# -- a vacuous "0 == 0"-class agreement that could not have caught the item
# 3/13x CW bug. Run it instead at the highest current under this card's own
# transport parameters where mu is still inside the cap-2 domain (0.01 < mu
# < 1) -- I=5e-4 uA, mu=0.64, the dot's own f_qfl (loading suppression at
# E_X) = 0.24, i.e. meaningfully turned on rather than the previous probe's
# ~1e-8 -- both rho values must still derive from the SAME transport
# background/X-rate accounting. (f_qfl -> 1 only above mu ~ 5 for this
# card's diode, outside the cap-2 domain -- 0.24 is the closest to full
# turn-on reachable while 0.01 < mu < 1 holds.)
d_linear = copy.deepcopy(d_gainp)
d_linear.drive.I_uA = 5e-4
sc_linear = evaluate(d_linear)["scalars"]
Tj_linear = sc_linear["T_j_op"]
diode_linear = _diode_from_drive(d_linear.drive, d_linear.aperture.diameter_um)
E_X_linear = _confinement_params(d_linear.ret, Tj_linear)["E_X_eV"]
V_j_linear = diode_linear.vj_of_j(d_linear.drive.I_uA * 1e-6 / diode_linear.area_cm2, Tj_linear)
f_qfl_linear = transport.qfl_suppression(E_X_linear, V_j_linear, materials.KB_EV * Tj_linear)
ok("linear-regime probe is not vacuous: 0.01 < mu < 1 and the dot's own f_qfl is well off "
   "the ~1e-8 previous probe (meaningfully turned on, not deep sub-turn-on)",
   0.01 < sc_linear["mu_resolved"] < 1.0 and f_qfl_linear > 0.1)
ok("CW rho matches pulsed rho within 2% in the linear pump regime",
   abs(sc_linear["cw_rho_op"] - sc_linear["rho_op"]) < 0.02 * sc_linear["rho_op"])

# Item 4, round 2 (council review 2026-09-05): round 1 made cavity-less
# track_material do SOMETHING (it had been gated behind cavity.enabled and
# was a no-op) -- but the something was wrong: it held the filter window
# FIXED at cavity.T_track and let the dot's own Varshni-walked line drift
# out of it (t_X 0.5 -> 0.00106 at the worst sweep corner, the regression's
# largest single term). A TRACKING filter is centred ON the dot at every T
# (hold_window default False): dx = filter.dx, so t_X must equal the
# independently-computed centred-window value, and must be UNCHANGED from
# no track_material at all. filter.hold_window=True is the explicit,
# differently-named opt-in for the fixed-window behaviour.
d_notrack = DeviceDesign(); d_notrack.ret.mode = "confinement"; d_notrack.ret.system = EDGE_SYSTEM
d_notrack.thermal.T_hs = 300.0; d_notrack.drive.V = 0.0
r_notrack = evaluate(d_notrack, [300.0])
d_track4 = copy.deepcopy(d_notrack); d_track4.filter.track_material = "dot"
r_track4 = evaluate(d_track4, [300.0])
w_track4 = r_track4["scalars"]["gamma_op"]  # auto_w default: w = Gamma(Tj)
spec_centred = epsilon(d_track4.dot.delta_xx, w_track4, d_track4.dot.r_xx * w_track4,
                      w=w_track4, kappa=None, dx=d_track4.filter.dx)
ok("filter.track_material (hold_window=False, default) centres the window on the dot: "
   "t_X at 300 K equals the independently-computed centred-window value",
   abs(r_track4["scalars"]["t_x_op"] - spec_centred.t_x) < 1e-6)
ok("centred tracking leaves t_X unchanged from no track_material at all (no longer an artefact)",
   r_track4["scalars"]["t_x_op"] == r_notrack["scalars"]["t_x_op"])

d_hold4 = copy.deepcopy(d_track4); d_hold4.filter.hold_window = True
r_hold4 = evaluate(d_hold4, [300.0])
Tj_hold4 = r_hold4["scalars"]["T_j_op"]
mat_hold4 = _tracked_material(d_hold4.ret, "dot")
dx_hold4 = 1e3 * (materials.bandgap(mat_hold4, Tj_hold4)
                  - materials.bandgap(mat_hold4, d_hold4.cavity.T_track))
w_hold4 = r_hold4["scalars"]["gamma_op"]
spec_hold4 = epsilon(d_hold4.dot.delta_xx, w_hold4, d_hold4.dot.r_xx * w_hold4,
                     w=w_hold4, kappa=None, dx=dx_hold4)
ok("filter.hold_window=True is the explicit, differently-named opt-in for a fixed window: "
   "t_X at 300 K matches the held Varshni-walk offset, not the centred value",
   abs(r_hold4["scalars"]["t_x_op"] - spec_hold4.t_x) < 1e-6
   and r_hold4["scalars"]["t_x_op"] != r_track4["scalars"]["t_x_op"])

# ============================================ 13. council review 2026-09-05 round 3
# Item 1: repetition rate must be explicit under EL-transport pulsed drive
# (drive.cw=False); the untouched defaults (duty=1.0, tau_pulse_ns=1.0) used
# to silently report a 1 GHz/100%-duty DC drive.
d_rep_default = DeviceDesign(); d_rep_default.drive.mode = "EL-transport"
d_rep_default.drive.diode = {"preset": "hkust"}
raises("item 1: EL-transport pulsed drive with default duty/tau_pulse/rep_rate_hz raises",
       lambda: evaluate(d_rep_default, [d_rep_default.thermal.T_hs]))

d_rep_tau_only = DeviceDesign(); d_rep_tau_only.drive.mode = "EL-transport"
d_rep_tau_only.drive.diode = {"preset": "hkust", "tau_pulse_ns": 0.1}
raises("item 1: tau_pulse_ns alone (duty and rep_rate_hz both still default) still raises",
       lambda: evaluate(d_rep_tau_only, [d_rep_tau_only.thermal.T_hs]))

d_rep_duty_only = DeviceDesign(); d_rep_duty_only.drive.mode = "EL-transport"
d_rep_duty_only.drive.diode = {"preset": "hkust"}; d_rep_duty_only.drive.duty = 0.008
raises("item 1: duty alone (tau_pulse_ns still default) still raises",
       lambda: evaluate(d_rep_duty_only, [d_rep_duty_only.thermal.T_hs]))

d_rep_cw_exempt = DeviceDesign(); d_rep_cw_exempt.drive.mode = "EL-transport"
d_rep_cw_exempt.drive.diode = {"preset": "hkust"}; d_rep_cw_exempt.drive.cw = True
try:
    evaluate(d_rep_cw_exempt, [d_rep_cw_exempt.thermal.T_hs])
    ok("item 1: drive.cw=True is exempt from the explicit tau_pulse/rep-rate requirement", True)
except ValueError:
    ok("item 1: drive.cw=True is exempt from the explicit tau_pulse/rep-rate requirement", False)

# pr-pkg4-fix item 5: drive.cw=True is exempt from the explicit tau_pulse/
# rep-rate requirement above ONLY while finite_pulse stays False (the CW
# path itself has no pulse period to fabricate) -- a CW card that ALSO
# opts into drive.finite_pulse needs a real pulse period exactly like the
# non-CW pulsed path, and must raise rather than silently reading the
# untouched duty=1.0/tau_pulse_ns=1.0 defaults as a fabricated 1 GHz drive.
d_cw_fp_default = DeviceDesign(); d_cw_fp_default.drive.mode = "EL-transport"
d_cw_fp_default.drive.diode = {"preset": "hkust"}
d_cw_fp_default.drive.cw = True; d_cw_fp_default.drive.finite_pulse = True
raises("item 5: drive.cw=True with drive.finite_pulse=True and default "
      "duty/tau_pulse/rep_rate_hz raises (no fabricated 1 GHz drive)",
      lambda: evaluate(d_cw_fp_default, [d_cw_fp_default.thermal.T_hs]))

d_cw_fp_ok = DeviceDesign(); d_cw_fp_ok.drive.mode = "EL-transport"
d_cw_fp_ok.drive.diode = {"preset": "hkust", "tau_pulse_ns": 0.1}
d_cw_fp_ok.drive.rep_rate_hz = 8.0e7; d_cw_fp_ok.drive.n_dot_cm2 = 1e10
d_cw_fp_ok.drive.cw = True; d_cw_fp_ok.drive.finite_pulse = True
try:
    evaluate(d_cw_fp_ok, [d_cw_fp_ok.thermal.T_hs])
    ok("item 5: drive.cw=True with drive.finite_pulse=True and an explicit "
      "tau_pulse_ns/rep_rate_hz does not raise", True)
except ValueError:
    ok("item 5: drive.cw=True with drive.finite_pulse=True and an explicit "
      "tau_pulse_ns/rep_rate_hz does not raise", False)

# pr-pkg4-fix item 6: drive.mechanism resolves its own (mu, F_p, eta), which
# the finite-pulse rates below never consult -- combining the two would
# silently ignore the mechanism, so it must raise instead. Left at
# drive.mode's "EL" default (not "EL-transport") so this exercises ONLY the
# new finite_pulse+mechanism guard, not the pre-existing EL-transport+
# mechanism guard above.
d_fp_mech = DeviceDesign()
d_fp_mech.drive.finite_pulse = True
d_fp_mech.drive.mechanism = "poisson-rail"
raises("item 6: drive.finite_pulse=True combined with drive.mechanism raises",
      lambda: evaluate(d_fp_mech, [d_fp_mech.thermal.T_hs]))

d_rep_ok = DeviceDesign(); d_rep_ok.drive.mode = "EL-transport"
d_rep_ok.drive.diode = {"preset": "hkust", "tau_pulse_ns": 0.1}
d_rep_ok.drive.rep_rate_hz = 8.0e7
d_rep_ok.drive.n_dot_cm2 = 1e10
sc_rep = evaluate(d_rep_ok, [d_rep_ok.thermal.T_hs])["scalars"]
ok("item 1: explicit tau_pulse_ns=0.1 ns / rep_rate_hz=80 MHz resolves rep_rate_hz = 8e7",
   abs(sc_rep["rep_rate_hz"] - 8.0e7) < 1.0)
ok("item 1: resolved tau_pulse_ns matches the card's explicit value",
   abs(sc_rep["tau_pulse_ns"] - 0.1) < 1e-9)
diode_rep = _diode_from_drive(d_rep_ok.drive, d_rep_ok.aperture.diameter_um)
ld_rep = diode_rep.dot_loading(d_rep_ok.drive.I_uA, sc_rep["T_j_op"], d_rep_ok.drive.n_dot_cm2,
                               float(np.pi * (d_rep_ok.aperture.diameter_um / 2) ** 2),
                               tau_pulse_ns=0.1)
ok("item 1: mu is consistent with the explicit tau_pulse_ns (mu is independent of the "
   "duty/rep-rate bookkeeping -- mu = r_dot * tau_pulse_ns)",
   abs(sc_rep["mu_resolved"] - ld_rep.mu) < 1e-9 * max(1.0, ld_rep.mu))

# duty=0.008 (explicit) and rep_rate_hz=8e7 (== 0.1ns*0.008/1e-9... i.e. the
# SAME physical repetition rate stated the other way) must resolve to the
# same rep_rate_hz and duty_resolved either way.
d_rep_duty = DeviceDesign(); d_rep_duty.drive.mode = "EL-transport"
d_rep_duty.drive.diode = {"preset": "hkust", "tau_pulse_ns": 0.1}
d_rep_duty.drive.duty = 8.0e-3
d_rep_duty.drive.n_dot_cm2 = 1e10
sc_rep_duty = evaluate(d_rep_duty, [d_rep_duty.thermal.T_hs])["scalars"]
ok("item 1: an equivalent explicit duty=0.008 (tau_pulse_ns=0.1 ns) resolves the SAME rep_rate_hz",
   abs(sc_rep_duty["rep_rate_hz"] - sc_rep["rep_rate_hz"]) < 1.0)
ok("item 1: duty_resolved reproduces the card's explicit duty",
   abs(sc_rep_duty["duty_resolved"] - 8.0e-3) < 1e-12)

# Item 2: the two rho values share the per-collected-X convention, but their
# signal rates are only equal in the linear-loading limit.  The gaasp card is
# deliberately at mu=1 (the cap-2 flux maximum), so its pulsed
# P(n>=1)-per-pulse rate and CW steady-state rate must not be compared as if
# they were the same observable.
d_gaasp = DeviceDesign.load(ROOT / "cards" / "edge-inp-gaasp-design.yaml")
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    sc_gaasp = evaluate(d_gaasp)["scalars"]

# In the linear regime the cap-2 pulsed X rate and the CW X rate have the
# same first-order loading.  Keep the card's b_res in this fixture: it is a
# per-collected-X channel and must survive the current reduction unchanged.
d_gaasp_linear = copy.deepcopy(d_gaasp)
d_gaasp_linear.drive.I_uA = 1.0e-5
with warnings.catch_warnings():
    warnings.simplefilter("ignore")
    sc_gaasp_linear = evaluate(d_gaasp_linear)["scalars"]
ok("item 2: pulsed rho_op equals cw_rho_op within 1e-5 relative in a dedicated "
   "gaasp linear-loading fixture (mu <= 0.01; b_res retained)",
   sc_gaasp_linear["mu_resolved"] <= 0.01 and
   abs(sc_gaasp_linear["rho_op"] - sc_gaasp_linear["cw_rho_op"])
   < 1e-5 * sc_gaasp_linear["rho_op"])

# The pulsed closed form includes both the transport background and the
# residual channel.  The latter was added after the original check [DR,
# device.py round-3 residual-background convention].
rho_closed_gaasp = 1.0 / (1.0 + sc_gaasp["b_e_resolved"] / sc_gaasp["t_x_op"]
                          + d_gaasp.drive.b_res)
ok("item 2: the fixed pulsed rho_op matches the closed form "
   "1/(1+b_e_resolved/t_x_op+b_res)",
   abs(sc_gaasp["rho_op"] - rho_closed_gaasp) < 1e-9)

# Quantify the card-point pulsed/CW relationship analytically, at whatever
# mu this card's own (re-solved) drive.I_uA currently resolves to.  Let
# f_cw be the CW collected X+XX signal divided by the low-pump collected-X
# rate, and f_p=(1-exp(-mu))/mu the cap-2 pulsed loading factor.  Thus the
# actual pulsed/CW signal-rate factor is f_p/f_cw.  For rho itself, the
# residual channel follows the X-only CW fraction; after removing that shared
# channel, the injection-background odds must transform by 1/f_cw.
#
# pr-pkg1-fix3 item 6: the previous version of this check asserted
# abs(ratio - 1) < 1.0 with no precondition -- vacuous (a ratio has to be
# wildly wrong, or non-finite compared against 1.0's own scale, to fail
# that). Restoring the ORIGINAL abs(rho_op - cw_rho_op) > 1e-2 precondition
# +1e-2 tolerance pairing does not work on this tree: swept over
# drive.I_uA at this card's own geometry (T_hs/dot/emission unchanged),
# the precondition holds (rho_op and cw_rho_op differ by >1e-2) almost
# everywhere EXCEPT within a couple of % of this card's own resolved
# I_uA -- but the ratio itself is only within 1e-2 of 1 in that SAME
# narrow near-linear neighborhood (e.g. mu=0.155961 at I_uA=0.006, +16%
# off this card's own I_uA=0.005161206, already misses 1e-2: ratio=
# 0.641187). The two conditions are not jointly satisfiable by any
# meaningful margin -- pairing them the original way would either never
# fire (mask a real regression) or fire and immediately fail (a false
# positive on an already-passing tree). This is finding 1b's density fix
# genuinely narrowing the identity's regime of validity, not a check that
# lost a precondition by accident.
#
# Corrected analytic expectation (spec: "replace it with the corrected
# analytic expectation"): AT this card's own resolved operating point
# (mu=0.0999999987, where rho_op and cw_rho_op agree to ~1e-7 already),
# the loading-factor identity holds to float-precision exactly (measured
# abs(ratio - 1) = 2.220446e-16 on this tree) -- assert that tight bound
# directly, scoped to the operating point it is actually true at, instead
# of a "must differ meaningfully" precondition this tree's own I_uA no
# longer satisfies.
Tj_gaasp = sc_gaasp["T_j_op"]
# pr-pkg1-fix2 item 4: forward T_ref=Tj (via Tj_gaasp, already the second
# positional arg) and this card's own resolved n_dot_cm2 (aperture.
# density_cm2=3e8, not _confinement_params' internal 1e10 default) so this
# reconstruction matches what evaluate() actually used for this card.
p_gaasp = _confinement_params(
    d_gaasp.ret, Tj_gaasp,
    # pr-pkg1-fix4 item 2: _confinement_params itself already checks
    # ret.n_dot_cm2 first, so the caller-supplied fallback here must be
    # ONLY the explicit aperture.density_cm2 -- drive.n_dot_cm2 is
    # TRANSPORT-only and must never reach confinement (see
    # _confinement_params' own docstring); the old `d_gaasp.drive.
    # n_dot_cm2 or ...` chain leaked it (this card's drive.n_dot_cm2 is
    # unset/0.0, so it happened not to matter here).
    n_dot_cm2=d_gaasp.aperture.density_cm2)
gamma_gaasp = 1.0 / d_gaasp.ret.tau_rad_ns
kx_gaasp, kxx_gaasp = cw_g2.escape_rates_from_retention(
    gamma_gaasp, p_gaasp["a_esc"], p_gaasp["E_a"], p_gaasp["b_p"], p_gaasp["E_b"], Tj_gaasp)
r_gaasp_ns = sc_gaasp["loading.r_dot"] / 1e9
tx_gaasp = sc_gaasp["t_x_op"]
txx_gaasp = sc_gaasp["eps_op"] * tx_gaasp
ix_gaasp, ixx_gaasp = cw_g2.photon_rates(
    r_gaasp_ns, gamma_gaasp, 2.0 * gamma_gaasp, kx_gaasp, kxx_gaasp,
    d_gaasp.drive.cw_pump_ratio)
sig_gaasp = tx_gaasp * ix_gaasp + txx_gaasp * ixx_gaasp
f_cw_gaasp = sig_gaasp / (tx_gaasp * r_gaasp_ns * sc_gaasp["S_resolved"])
mu_gaasp = sc_gaasp["mu_resolved"]
P_ge1_gaasp = 1.0 - loading_probs(mu_gaasp)[0]
f_p_gaasp = P_ge1_gaasp / mu_gaasp
loading_rate_factor_gaasp = f_p_gaasp / f_cw_gaasp
bp_inj_gaasp = (1.0 / sc_gaasp["rho_op"] - 1.0 - d_gaasp.drive.b_res)
cw_x_fraction_gaasp = tx_gaasp * ix_gaasp / sig_gaasp
bc_inj_gaasp = (1.0 / sc_gaasp["cw_rho_op"] - 1.0
                - d_gaasp.drive.b_res * cw_x_fraction_gaasp)
ratio_gaasp = (bc_inj_gaasp / bp_inj_gaasp) / (1.0 / f_cw_gaasp)
# pr-pkg1-fix4 item 6: tightened 1e-6 -> 1e-8 -- measured abs(ratio - 1) =
# 2.2e-16 at this operating point, while dropping the 1/f_cw factor gives
# 2.9e-6, so 1e-8 discriminates the two by ~300x (1e-6 did not).
ok("item 6 (pr-pkg1-fix3): gaasp card-point rho difference follows the analytic cap-2/CW "
   "loading factor (P(n>=1)/mu divided by CW saturation factor) to float precision AT "
   "this card's own resolved operating point (tolerance 1e-8, not the previous vacuous "
   "<1.0)",
   np.isfinite(loading_rate_factor_gaasp) and loading_rate_factor_gaasp < 1.0
   and np.isfinite(ratio_gaasp) and abs(ratio_gaasp - 1.0) < 1e-8)

# Item 6: photon_budget replaces the tautological carrier_budget_closure --
# it must equal 1 within 1e-6 at a real operating point, AND the self-test
# below must show the check has teeth: dropping one of the four exposed
# channels from the numerator moves the ratio measurably off 1.
sc_pb = evaluate(d_gainp)["scalars"]
ok("item 6: photon_budget closes to 1 within 1e-6 at the gainp card's operating point",
   abs(sc_pb["photon_budget"] - 1.0) < 1e-6)
channels_pb = ("photon_budget.dot_radiative", "photon_budget.dot_nonradiative",
              "photon_budget.matrix_radiative", "photon_budget.matrix_nonradiative")
supply_pb = sc_pb["photon_budget.supply_active"]
full_sum_pb = sum(sc_pb[k] for k in channels_pb)
ok("item 6: the four exposed channels independently sum to supply_active (consistent with "
   "the reported photon_budget)", abs(full_sum_pb / supply_pb - sc_pb["photon_budget"]) < 1e-9)
# pr-pkg1-fix2 items 1+4: the confinement escape-prefactor/T_ref fixes
# collapse S (hence dot_radiative's share of supply_active) at the gainp
# card's own operating point -- dot_radiative is now ~1.6e-4 of the budget
# (was ~1e-2 order before those fixes), so the self-test's threshold is
# lowered from 1e-3 to 1e-5 to stay well below every channel's real share
# (dot_radiative ~1.6e-4, matrix_radiative ~1.4e-2, matrix_nonradiative
# ~1.3e-1, dot_nonradiative ~0.86) while staying far above numerical noise
# -- it still catches a channel silently contributing ~0 to the sum.
for dropped in channels_pb:
    dropped_sum = sum(sc_pb[k] for k in channels_pb if k != dropped) / supply_pb
    ok(f"item 6 self-test: dropping {dropped} from the sum moves the ratio measurably off 1 "
       "(the diagnostic has teeth -- the OLD carrier_budget_closure could never fail like this)",
       abs(dropped_sum - 1.0) > 1e-5)

# Item 7: track_material provenance must say "inert" when cavity is disabled
# and hold_window is False (the setting changes nothing), and must claim
# "materials.bandgap ..." only when it is actually effective.
d_prov_inert = DeviceDesign(); d_prov_inert.ret.mode = "confinement"
d_prov_inert.ret.system = EDGE_SYSTEM; d_prov_inert.filter.track_material = "dot"
r_prov_inert = evaluate(d_prov_inert, [300.0])
ok("item 7: inert track_material (cavity disabled, hold_window false) reports the honest "
   "'inert' provenance note, not a false 'materials.bandgap ...' claim",
   r_prov_inert["scalars"]["provenance"]["tracking"]["note"]
   == "track_material set but inert (cavity disabled, hold_window false)")
ok("item 7: inert track_material is NOT reported as track_material_effective",
   r_prov_inert["scalars"]["track_material_effective"] is False)

d_prov_active = DeviceDesign(); d_prov_active.cavity.enabled = True
d_prov_active.ret.mode = "confinement"; d_prov_active.ret.system = EDGE_SYSTEM
d_prov_active.filter.track_material = "dot"; d_prov_active.thermal.T_hs = 300.0
r_prov_active = evaluate(d_prov_active, [300.0])
ok("item 7: effective track_material (cavity enabled) reports the 'materials.bandgap ...' note",
   r_prov_active["scalars"]["provenance"]["tracking"]["note"]
   == "materials.bandgap on stack layer 'dot'")
ok("item 7: effective track_material IS reported as track_material_effective",
   r_prov_active["scalars"]["track_material_effective"] is True)

d_prov_hold = copy.deepcopy(d_prov_inert); d_prov_hold.filter.hold_window = True
r_prov_hold = evaluate(d_prov_hold, [300.0])
ok("item 7: cavity-less hold_window=True is ALSO effective (materials.bandgap note)",
   r_prov_hold["scalars"]["provenance"]["tracking"]["note"]
   == "materials.bandgap on stack layer 'dot'"
   and r_prov_hold["scalars"]["track_material_effective"] is True)

# Item 8: the residual background channel drive.b_res reproduces the
# Reischle 2008 80 K electrical anchor rho ~ 0.88 when the Urbach injection
# channel is negligible (rho = 1/(1+b_res) exactly, independent of G, S,
# whenever the other background terms vanish).
d_bres = DeviceDesign(); d_bres.thermal.T_hs = 80.0; d_bres.drive.mode = "EL-transport"
d_bres.drive.diode = {"preset": "hkust", "tau_pulse_ns": 0.1}; d_bres.drive.duty = 0.008
d_bres.drive.n_dot_cm2 = 1e10; d_bres.drive.I_uA = 1e-5
# confinement mode with RetentionBlock's own (0.0) b0/beta defaults, so the
# ONLY background channels present are the Urbach injection tail and b_res --
# the proxy-fit b0/beta terms (unrelated to either) would otherwise
# contaminate the "Urbach negligible" isolation this check needs.
d_bres.ret.mode = "confinement"; d_bres.ret.preset = "InP/GaAsP0.4/AlGaAs0.4 on GaAs"
d_bres.drive.b_res = 1.0 / 0.88 - 1.0
sc_bres = evaluate(d_bres, [80.0])["scalars"]
ok("item 8: at 80 K the Urbach injection channel is negligible (b_e_resolved << b_res)",
   sc_bres["b_e_resolved"] < 1e-3 * d_bres.drive.b_res)
ok("item 8: b_res = 1/0.88-1 reproduces rho_op = 0.88 within 1e-3",
   abs(sc_bres["rho_op"] - 0.88) < 1e-3)
d_bres0 = copy.deepcopy(d_bres); d_bres0.drive.b_res = 0.0
sc_bres0 = evaluate(d_bres0, [80.0])["scalars"]
ok("item 8: drive.b_res=0.0 (default) leaves rho_op at its legacy (no residual channel) value",
   sc_bres0["rho_op"] > 0.999)
ok("item 8: b_res provenance is tagged [E] with an anchor note when non-zero",
   sc_bres["provenance"]["b_res"]["tag"] == "E"
   and "anchor" in sc_bres["provenance"]["b_res"]["note"])
ok("item 8: b_res provenance says 'not requested' when 0.0 (legacy)",
   sc_bres0["provenance"]["b_res"]["note"] == "not requested (0.0, legacy)")

# Item 9: filter.auto_w_scale halves the resolved collection window and
# lowers eps at Gamma300 = 20 meV (the XX line sits partly inside the full-
# linewidth window; a narrower window excludes more of it).
d_w1 = DeviceDesign(); d_w1.dot.gamma300 = 20.0; d_w1.dot.linewidth = "anchored"
d_w1.dot.delta_xx = 4.0; d_w1.thermal.T_hs = 300.0; d_w1.drive.V = 0.0
r_w1 = evaluate(d_w1, [300.0])
d_w05 = copy.deepcopy(d_w1); d_w05.filter.auto_w_scale = 0.5
r_w05 = evaluate(d_w05, [300.0])
ok("item 9: auto_w_scale=0.5 exactly halves the resolved w",
   abs(r_w05["scalars"]["w_resolved"] - 0.5 * r_w1["scalars"]["w_resolved"]) < 1e-9)
ok("item 9: halving the window lowers eps at Gamma300=20 meV",
   r_w05["scalars"]["eps_op"] < r_w1["scalars"]["eps_op"])
ok("item 9: xx_in_window is a finite fraction and equals eps_op * t_x_op",
   abs(r_w1["scalars"]["xx_in_window"] - r_w1["scalars"]["eps_op"] * r_w1["scalars"]["t_x_op"]) < 1e-12)
ok("item 9: auto_w_scale=1.0 (default) reproduces the legacy w exactly",
   DeviceDesign().filter.auto_w_scale == 1.0)
ok("item 9: the resolved provenance documents the auto_w operating convention",
   "auto_w" in r_w1["scalars"]["provenance"]["filter_window"]["note"])

# slow: peer-review-triage.md finding 4 / pr-pkg4-fix item 1 -- gate-
# consistent counting. The pr-pkg4-fix (17f2c2f) round's Check A used a
# degenerate duty=1.0/tau_pulse_ns=12.5 construction (tau_dark_ns forced to
# exactly 0.0) that only ever exercised the CW-limit identity; it is
# replaced below by checks at the gainp card's own real 100 ps pulse.

# Check A (item 4a): pulse_counting.pulse_g2's gate-restricted X-line count
# (mean_counts_x, split=True) reduces to the legacy static per-pulse
# X-photon count loading.brightness_per_pulse == (P1+P2)*t_x -- an ordinary
# two-step Poisson-arrival argument (cw_g2.generator's 0->X->XX chain at
# rate r is memoryless, so with escape off and decay negligible over a
# vanishing on-window, P(0 events)=P0, P(1 event)=P1, P(>=2 events)=P2
# exactly, matching loading.loading_probs' own cap-2 Poisson convention) --
# ONLY in the short-pulse limit (mu = r_ns*tau_on held fixed at the card's
# own value while tau_on -> 0). The comparison is against mean_counts_x,
# not the TOTAL mean_counts finite_pulse_mean_counts/rho actually use:
# (P1+P2)*t_x is itself X-photon-only by construction (loading.
# brightness_per_pulse's docstring), so it omits the cascade's XX-line
# photon (mean_counts_xx) entirely -- that XX share is real detected signal
# (rho's own "signal = G*(n_X+n_XX)" counts it), not a finite-pulse
# artefact of the counting window, so the TOTAL mean_counts never converges
# to (P1+P2)*t_x even as tau_on -> 0 whenever eps > 0 (verified: at this
# card's eps_op ~0.65 the gap stays ~16% in that limit too). Escape is
# forced to zero (monkeypatching cw_g2.escape_rates_from_retention) so the
# comparison isolates the finite-pulse re-excitation effect from escape.
d_a = copy.deepcopy(d_gainp)
d_a.drive.finite_pulse = True
period_a = d_a.drive.diode["tau_pulse_ns"] / d_a.drive.duty
d_a.drive.gate_ns = period_a
_orig_escape_a = cw_g2.escape_rates_from_retention
cw_g2.escape_rates_from_retention = lambda *a, **k: (0.0, 0.0)
try:
    sc_a = evaluate(d_a)["scalars"]
finally:
    cw_g2.escape_rates_from_retention = _orig_escape_a
mu_a, t_X_a, eps_a = sc_a["mu_resolved"], sc_a["t_x_op"], sc_a["eps_op"]
gamma_X_a = sc_a["finite_pulse_gamma_X_ns"]; gamma_XX_a = sc_a["finite_pulse_gamma_XX_ns"]
tau_dark_a = sc_a["finite_pulse_tau_dark_ns"]
_, P1_a, P2_a = loading_probs(mu_a)
static_signal_a = (P1_a + P2_a) * t_X_a

tau_on_tiny = 1e-4
pc_tiny = pulse_counting.pulse_g2(mu_a / tau_on_tiny, gamma_X_a, gamma_XX_a, k_X=0.0, k_XX=0.0,
                                  t_X=t_X_a, t_XX=eps_a * t_X_a, tau_on_ns=tau_on_tiny,
                                  tau_dark_ns=tau_dark_a,
                                  # item 4 (pr-pkg4-fix2): explicit, not pulse_g2's own
                                  # default -- this check must not depend on the card's
                                  # cw_pump_ratio happening to equal pulse_g2's default.
                                  pump_ratio=d_a.drive.cw_pump_ratio, split=True,
                                  gate_ns=tau_on_tiny + tau_dark_a)
ok("item 1/finding 4 check A: pulse_g2's gate-restricted X-line count "
   "(mean_counts_x) reproduces the legacy static (P1+P2)*t_x to rtol 1e-3 "
   f"in the short-pulse limit (tau_on={tau_on_tiny:g} ns, mu held at the "
   f"card's own {mu_a:.6f}, escape forced to 0)",
   abs(pc_tiny["mean_counts_x"] - static_signal_a) < 1e-3 * static_signal_a)

# At the card's REAL 100 ps pulse this no longer holds exactly (re-excitation
# is a real, computable correction) -- frozen 2026-09-07 against a number
# derived with pulse_g2 at the card's own parameters (escape forced to 0,
# gate_ns=period), not asserted as an identity.
FP_MEAN_COUNTS_X_100PS_20260907 = 0.21043095348245136
dev_a = abs(sc_a["finite_pulse_mean_counts_x"] - static_signal_a) / static_signal_a
ok("item 1/finding 4 check A: at the card's real 100 ps pulse, "
   "finite_pulse_mean_counts_x reproduces the frozen 2026-09-07 value to "
   f"rtol 1e-9 (deviation from the short-pulse static value is {dev_a:.4%})",
   abs(sc_a["finite_pulse_mean_counts_x"] - FP_MEAN_COUNTS_X_100PS_20260907)
   < 1e-9 * FP_MEAN_COUNTS_X_100PS_20260907)

# Check B (item 4b): gate_ns=period vs. gate_ns=1 ns give DIFFERENT
# mean_counts (the pre-fix bug held mean_counts constant across gate_ns and
# only varied the background) -- the ratio of the device path's own exposed
# finite_pulse_mean_counts between the two gates must equal the ratio of an
# INDEPENDENT pulse_g2 call's mean_counts at the same two gates (same
# extracted rates), to 1e-9, proving finite_pulse_mean_counts really is
# pulse_g2's gate-restricted count and not some other, gate-blind quantity.
# Escape forced to zero: at this card's real (fast) escape rates essentially
# all counting is already resolved well inside 1 ns, so gate=period and
# gate=1 ns agree almost trivially (ratio 1.0) -- zero escape (slower,
# purely radiative decay) makes the two gates meaningfully different.
d_b = copy.deepcopy(d_gainp)
d_b.drive.finite_pulse = True
period_b = d_b.drive.diode["tau_pulse_ns"] / d_b.drive.duty
_orig_escape_b = cw_g2.escape_rates_from_retention
cw_g2.escape_rates_from_retention = lambda *a, **k: (0.0, 0.0)
try:
    d_b1 = copy.deepcopy(d_b); d_b1.drive.gate_ns = period_b
    d_b2 = copy.deepcopy(d_b); d_b2.drive.gate_ns = 1.0
    sc_b1 = evaluate(d_b1)["scalars"]
    sc_b2 = evaluate(d_b2)["scalars"]
finally:
    cw_g2.escape_rates_from_retention = _orig_escape_b
ratio_device = sc_b1["finite_pulse_mean_counts"] / sc_b2["finite_pulse_mean_counts"]
pc_b1 = pulse_counting.pulse_g2(
    sc_b1["finite_pulse_r_ns"], sc_b1["finite_pulse_gamma_X_ns"], sc_b1["finite_pulse_gamma_XX_ns"],
    sc_b1["finite_pulse_k_X"], sc_b1["finite_pulse_k_XX"], sc_b1["t_x_op"], sc_b1["eps_op"] * sc_b1["t_x_op"],
    sc_b1["finite_pulse_tau_on_ns"], sc_b1["finite_pulse_tau_dark_ns"],
    # item 4 (pr-pkg4-fix2): explicit, not pulse_g2's own default -- same
    # reasoning as check A above.
    pump_ratio=d_b.drive.cw_pump_ratio, gate_ns=period_b)
pc_b2 = pulse_counting.pulse_g2(
    sc_b1["finite_pulse_r_ns"], sc_b1["finite_pulse_gamma_X_ns"], sc_b1["finite_pulse_gamma_XX_ns"],
    sc_b1["finite_pulse_k_X"], sc_b1["finite_pulse_k_XX"], sc_b1["t_x_op"], sc_b1["eps_op"] * sc_b1["t_x_op"],
    sc_b1["finite_pulse_tau_on_ns"], sc_b1["finite_pulse_tau_dark_ns"],
    pump_ratio=d_b.drive.cw_pump_ratio, gate_ns=1.0)
ratio_direct = pc_b1["mean_counts"] / pc_b2["mean_counts"]
ok(f"item 1/finding 4 check B: mean_counts(gate=period)/mean_counts(gate=1 ns) "
   f"from the device path ({ratio_device!r}) matches an independent pulse_g2 "
   f"call at the same rates ({ratio_direct!r}) to 1e-9, proving signal is "
   "gate-restricted",
   abs(ratio_device - ratio_direct) < 1e-9 * ratio_direct)

# Check C (item 1/4c, rewritten -- pr-pkg4-fix2, orchestrator decision):
# item 1's background-window fix caps n_bg at min(gate_ns, tau_pulse_ns),
# so B_fp can only grow while the gate is still inside the pump window --
# once gate_ns >= tau_pulse_ns, background stops growing entirely.
# mean_counts (signal) keeps growing a little past tau_pulse_ns too, since
# the cascade's own radiative/escape decay tail has not yet fully resolved
# at gate_ns == tau_pulse_ns exactly -- so rho_pulsed briefly RISES while
# background is already capped but signal is not. Once the cascade has
# also fully resolved (a few tau_rad past the pulse -- well inside 1 ns at
# this card's real, fast escape rates) both numerator and denominator are
# constant and rho_pulsed is flat. rho_pulsed can therefore only rise then
# plateau in gate_ns -- never fall -- which replaces the pre-item-1
# "background keeps growing forever, so rho keeps falling forever"
# property the original (pkg4, pre-pkg4-fix2) check C asserted.
d_gate = copy.deepcopy(d_gainp)
d_gate.drive.finite_pulse = True
gate_grid = (0.02, 0.05, 0.1, 0.15, 1.0, 3.0, 6.0, 12.5)
rho_gate = []
for g in gate_grid:
    dg = copy.deepcopy(d_gate)
    dg.drive.gate_ns = g
    rho_gate.append(evaluate(dg)["scalars"]["rho_pulsed"])
ok(f"item 1/finding 4 check C(a): rho_pulsed is non-decreasing in gate_ns "
   f"over {gate_grid} ns (got {rho_gate})",
   all(np.isfinite(r) for r in rho_gate)
   and all(rho_gate[i] <= rho_gate[i + 1] for i in range(len(rho_gate) - 1)))

# gate_grid's own last point IS the card's real period (tau_pulse_ns/duty
# = 0.1/0.008 = 12.5 ns) -- reuse rho_gate[-1] rather than re-evaluating.
tau_pulse_c = d_gate.drive.diode["tau_pulse_ns"]
period_c = tau_pulse_c / d_gate.drive.duty
ok("item 1/finding 4 check C(b) precondition: gate_grid's last point must be "
   "the card's own period", abs(period_c - gate_grid[-1]) < 1e-9)
rho_period_c = rho_gate[-1]
# item 5 (pr-pkg4-fix3), corrected by item 4 (pr-pkg4-fix4, Opus review):
# the previous absolute-1e-12 comparison implicitly assumed tau_rad_ns
# alone sets the cascade's decay tail, and a since-removed "escape-aware
# effective decay time" max(tau_rad_ns, 1/(k_X+gamma_X)) was computed here
# to guard against a slower escape channel lengthening it -- but escape
# (k_X > 0) only ADDS to the total decay rate, gamma_X + k_X >= gamma_X,
# so 1/(k_X+gamma_X) <= 1/gamma_X = tau_rad_ns ALWAYS (rm >= 1 keeps this
# true with a cavity: gamma_X here is the Purcell-enhanced rate rm*gamma_X_bare,
# and rm >= 1 -- enhancement, never suppression, the only regime these
# cards ship -- means gamma_X >= gamma_X_bare too, so 1/gamma_X stays
# <= tau_rad_ns = 1/gamma_X_bare): the max() could
# never actually pick the escape-derived timescale, so that computation
# silently reduced to plain tau_rad_ns in every case anyway. Escape only
# ever SHORTENS the true decay time, so gate = tau_pulse_ns + 5*tau_rad_ns
# (DriveBlock.gate_ns's own recommendation) is a conservative (wider than
# the true tail needs) cover of the cascade's decay tail regardless of
# escape, never too narrow -- a relative (not absolute) tolerance still
# lets the claim survive whichever escape rate this or a future card
# ships with.
gate_tail_c = tau_pulse_c + 5.0 * d_gate.ret.tau_rad_ns
d_tail = copy.deepcopy(d_gate); d_tail.drive.gate_ns = gate_tail_c
rho_tail_c = evaluate(d_tail)["scalars"]["rho_pulsed"]
ok(f"item 1/finding 4 check C(b): rho_pulsed(gate=period={period_c:g} ns) "
   f"({rho_period_c!r}) equals rho_pulsed(gate=tau_pulse+5*tau_rad="
   f"{gate_tail_c:g} ns) ({rho_tail_c!r}) to rtol 1e-9 -- both gates fully "
   "cover the cascade's decay tail, so both sit on the same plateau",
   np.isfinite(rho_period_c) and np.isfinite(rho_tail_c)
   and abs(rho_period_c - rho_tail_c) < 1e-9 * abs(rho_period_c))

ok(f"item 1/finding 4 check C(c): rho_pulsed(gate=0.02 ns, still inside the "
   f"pump window) ({rho_gate[0]!r}) is strictly less than rho_pulsed(gate="
   f"period) ({rho_period_c!r}) -- the rise while the pump is still open",
   np.isfinite(rho_gate[0]) and np.isfinite(rho_period_c)
   and rho_gate[0] < rho_period_c)

# Check D (item 4d): with every background channel forced to exactly zero
# (injection background via diode.eta_rad_matrix=0.0, and b_res=0.0 --
# ret.b0/ret.beta are already 0.0, this card's own defaults), B_fp above is
# identically 0 and rho = signal/(signal+0) = 1.0 to float round-off.
d_d = copy.deepcopy(d_gainp)
d_d.drive.finite_pulse = True
d_d.drive.gate_ns = d_d.drive.diode["tau_pulse_ns"] / d_d.drive.duty
d_d.drive.diode = dict(d_d.drive.diode); d_d.drive.diode["eta_rad_matrix"] = 0.0
d_d.drive.b_res = 0.0
sc_d = evaluate(d_d)["scalars"]
ok("item 1/finding 4 check D: rho at gate=period with every background "
   f"channel forced to zero equals 1.0 to 1e-12 (got {sc_d['rho_pulsed']!r})",
   sc_d["finite_pulse_converged"] and np.isfinite(sc_d["finite_pulse_mean_counts"])
   and sc_d["finite_pulse_mean_counts"] > 0
   and abs(sc_d["rho_pulsed"] - 1.0) < 1e-12)

# Check E (item 2, pr-pkg4-fix2): the rho_pulsed == cw_rho_op anchor
# (pr-pkg4-fix's original Check A, deleted when checks A-D above replaced
# it with the gainp card's real 100 ps pulse) is restored here as its own
# check: a test-only, degenerate 100%-duty construction (drive.duty=1.0,
# drive.diode['tau_pulse_ns'] = the card's own 12.5 ns period, so
# tau_dark_ns resolves to exactly 0.0) with drive.gate_ns covering that
# whole period. With tau_dark_ns=0 the periodic map Phi reduces to
# expm(M_on*T) alone, so the periodic steady state IS the ordinary CW
# steady state and pulse_counting's per-period m1 (mean_counts) is exactly
# sig*T -- an exact identity, independent of escape, that proves the
# gated-count rho reduces to the existing CW rho in the CW limit rather
# than merely resembling it. Escape is ALSO forced to zero (monkeypatching
# cw_g2.escape_rates_from_retention) to keep the check independent of
# whichever retention prefactor a future package lands on.
#
# pr-pkg4-fix's own item 5 guard (added after this check was first written
# and before it was deleted) requires an EXPLICIT duty or rep_rate_hz away
# from DriveBlock's own untouched defaults whenever finite_pulse needs a
# real pulse period -- setting drive.duty = 1.0 is indistinguishable from
# never having touched it (1.0 IS DriveBlock.duty's default), so the
# 100%-duty construction is expressed via an explicit rep_rate_hz instead
# (1e9/12.5 = 8e7 Hz at tau_pulse_ns=12.5 gives duty_eff = 1.0 and
# tau_dark_ns = 0.0 exactly, the same construction as before).
d_cwlimit = copy.deepcopy(d_cw_gainp)
d_cwlimit.drive.finite_pulse = True
d_cwlimit.drive.gate_ns = 12.5
d_cwlimit.drive.diode = dict(d_cwlimit.drive.diode)
d_cwlimit.drive.diode["tau_pulse_ns"] = 12.5
d_cwlimit.drive.rep_rate_hz = 1e9 / 12.5
if not _cw_reduction_guard(d_cwlimit, "item 2/finding 4 check E"):
    _orig_escape_e = cw_g2.escape_rates_from_retention
    cw_g2.escape_rates_from_retention = lambda *a, **k: (0.0, 0.0)
    try:
        sc_cwlimit = evaluate(d_cwlimit)["scalars"]
    finally:
        cw_g2.escape_rates_from_retention = _orig_escape_e
    ok("item 2/finding 4 check E: gate_ns covering the whole (100%-duty, "
       "zero-escape) period makes rho_pulsed exactly reproduce the CW branch's "
       "own rho (rtol 1e-3)",
       sc_cwlimit["finite_pulse_tau_dark_ns"] == 0.0
       and np.isfinite(sc_cwlimit["rho_pulsed"]) and np.isfinite(sc_cwlimit["cw_rho_op"])
       and abs(sc_cwlimit["rho_pulsed"] - sc_cwlimit["cw_rho_op"]) < 1e-3 * sc_cwlimit["cw_rho_op"])

# Check F (item 2, pr-pkg4-fix2): item 1's actual consequence at the
# card's REAL operating point -- 0.8% duty (tau_dark_ns >> tau_pulse_ns,
# not check E's degenerate duty=1.0 limit) and the card's own (fast,
# non-zero) escape rates, unforced. gate_ns=period, so signal
# (finite_pulse_mean_counts, integrated over the whole period) and the
# item-1-fixed background (integrated over only the min(gate, tau_pulse)
# pump window) sit on the same footing the CW branch's own rate ratio
# does -- the only check that anchors rho_pulsed against cw_rho_op OUTSIDE
# check E's artificial construction, at the card's actual pulsed drive.
d_f = copy.deepcopy(d_gainp)
d_f.drive.finite_pulse = True
d_f.drive.gate_ns = d_f.drive.diode["tau_pulse_ns"] / d_f.drive.duty
if not _cw_reduction_guard(d_f, "item 2/finding 4 check F"):
    sc_f = evaluate(d_f)["scalars"]
    ok("item 2/finding 4 check F: at the card's real 0.8% duty with escape "
       "rates as they are, rho_pulsed equals cw_rho_op to rtol 1e-3 for "
       f"gate=period (rho_pulsed={sc_f['rho_pulsed']!r}, "
       f"cw_rho_op={sc_f['cw_rho_op']!r})",
       sc_f["finite_pulse_converged"] and np.isfinite(sc_f["rho_pulsed"])
       and np.isfinite(sc_f["cw_rho_op"])
       and abs(sc_f["rho_pulsed"] - sc_f["cw_rho_op"]) < 1e-3 * sc_f["cw_rho_op"])

# pr-pkg4-fix4 item 1 (Opus review): _cw_reduction_guard itself is under
# test here -- a proxy-mode variant of the gainp card resolves ret.b0/beta
# to nonzero through the class-proxy Arrhenius fit (rp in device.py) even
# though the card's own raw RetentionBlock fields stay untouched at their
# proxy-mode default of 0.0. The OLD guard (design.ret.b0/design.ret.beta)
# would have wrongly let this card through checks E/F; the fixed guard
# reads the RESOLVED params and must skip.
d_proxy_variant = copy.deepcopy(d_gainp)
d_proxy_variant.ret.mode = "proxy"
_proxy_skipped = _cw_reduction_guard(d_proxy_variant, "item 1 proxy-mode self-test")
# pr-final-nits item 3: assert the RESOLVED b0 is actually nonzero
# explicitly, rather than only checking the guard skipped and the raw
# fields stayed 0.0 -- either of those two alone is also true if the
# guard skipped for the unrelated reason of the cavity being enabled, so
# without this the self-test could pass without ever exercising the
# proxy-fit resolution path it claims to.
_proxy_b0_used = evaluate(
    d_proxy_variant, T_grid=[d_proxy_variant.thermal.T_hs]
)["scalars"]["retention_params_used"]["b0"]
ok("item 1: a proxy-mode variant of the gainp card resolves b0/beta "
   "nonzero via the class-proxy fit even though design.ret.b0=="
   f"{d_proxy_variant.ret.b0!r}/design.ret.beta=={d_proxy_variant.ret.beta!r} "
   f"stay 0.0 -- the guard skips E/F on the RESOLVED value (b0_used={_proxy_b0_used!r} "
   "!= 0.0), not the raw fields",
   _proxy_skipped and d_proxy_variant.ret.b0 == 0.0 and d_proxy_variant.ret.beta == 0.0
   and _proxy_b0_used != 0.0)

# Item 7: an unconverged pulse_counting.pulse_g2 periodic steady state must
# nan g2_dot/rho and mark the row invalid the way the other early-return
# operating points do, rather than reporting a number from a bad fixed
# point -- monkeypatch pulse_g2 to flip a real, otherwise-valid result's
# own `converged` flag to False and confirm device.py actually reads it.
d_unconv = copy.deepcopy(d_gainp)
d_unconv.drive.finite_pulse = True
_orig_pulse_g2 = devmod.pulse_counting.pulse_g2
def _fake_unconverged(*a, **k):
    r = _orig_pulse_g2(*a, **k)
    r["converged"] = False
    return r
devmod.pulse_counting.pulse_g2 = _fake_unconverged
try:
    sc_unconv = evaluate(d_unconv)["scalars"]
finally:
    devmod.pulse_counting.pulse_g2 = _orig_pulse_g2
ok("item 7: an unconverged pulse_counting.pulse_g2 result nans g2_dot/rho, "
   "marks the row invalid, and exposes finite_pulse_converged=False",
   sc_unconv["finite_pulse_converged"] is False
   and np.isnan(sc_unconv["finite_pulse_g2_dot"])
   and np.isnan(sc_unconv["rho_pulsed"])
   and np.isnan(sc_unconv["g2_op"])
   and any("did not converge" in r for r in sc_unconv["invalid_reasons"]))
# pr-pkg4-fix3 item 2: finite_pulse_gate_ns_used must ALSO stay NaN in the
# unconverged case (assigned only after the convergence test succeeds), so
# the "finite_pulse_waveform" provenance note actually falls back to the
# "requested; operating point invalid (<reason>)" text instead of silently
# reporting a bogus gate-restricted-counting note.
ok("item 2 (pr-pkg4-fix3): unconverged finite_pulse leaves "
   "finite_pulse_gate_ns_used NaN and the finite_pulse_waveform "
   "provenance note reports 'requested; operating point invalid ...' "
   f"(got note={sc_unconv['provenance']['finite_pulse_waveform']['note']!r})",
   np.isnan(sc_unconv["finite_pulse_gate_ns_used"])
   and "requested; operating point invalid" in
       sc_unconv["provenance"]["finite_pulse_waveform"]["note"])

# ============================================ audit 2026-09-23 (device composition)
# Q1 (user decision): the legacy non-diode brightness omits S BY DESIGN and is
# labelled "pre-retention" (module docstring note + result label); no number
# moves (the legacy fixture comparisons above stay bit-identical).
_q1 = DeviceDesign(); _q1.thermal.T_hs = 300.0
_q1s = evaluate(_q1, [300.0])["scalars"]
_q1_P = loading_probs(_q1s["mu_resolved"])
ok("Q1: device.py module docstring notes that retention S is intentionally omitted "
   "on the legacy non-diode brightness path (diode path carries it)",
   "retention factor S is intentionally omitted on the legacy non-diode brightness path"
   in (devmod.__doc__ or "") and "diode path carries it" in (devmod.__doc__ or ""))
ok("Q1: legacy non-diode brightness is labelled pre-retention (brightness_convention "
   "and provenance), and still equals (P1+P2)*t_x with no S factor",
   _q1s["brightness_convention"].startswith("pre-retention")
   and _q1s["provenance"]["brightness"]["note"].startswith("pre-retention")
   and _q1s["S_resolved"] < 1.0
   and _q1s["brightness_per_pulse"] == float(_q1_P[1] + _q1_P[2]) * _q1s["t_x_op"])

# Audit item 4 (device.py:1860): a static per-pulse g2_op > 1 carries an
# explicit invalid flag. (a) the real static f8 path at a super-Poissonian
# pump (mu=0.05, F_p=1.9): g2_load = 1 + (F_p-1)/mu = 19 pushes g2_op past 1;
# (b) the pre-fix small-mu artefact form of f1b (P2 = 1 - P0 - P1, catastrophic
# cancellation) reproduced by monkeypatch at mu = 5.4e-10 (audit: g2_op ~ 111).
_i4 = DeviceDesign(); _i4.thermal.T_hs = 80.0; _i4.drive.b_e = 0.0
_i4.drive.mu = 0.05; _i4.drive.F_p = 1.9
_i4s = evaluate(_i4, [80.0])["scalars"]
_orig_f1b = devmod.f1b_g2
def _f1b_prefix(mu, eps):
    P0 = np.exp(-mu); P1 = mu * np.exp(-mu); P2 = 1.0 - P0 - P1
    den = P1 + P2 * (1.0 + eps)
    return np.where(den > 0, 2.0 * P2 * eps / den ** 2, 0.0)
_i4b = DeviceDesign(); _i4b.thermal.T_hs = 80.0; _i4b.drive.b_e = 0.0; _i4b.drive.mu = 5.4e-10
devmod.f1b_g2 = _f1b_prefix
try:
    _i4bs = evaluate(_i4b, [80.0])["scalars"]
finally:
    devmod.f1b_g2 = _orig_f1b
_i4c = evaluate(_set_thermal(DeviceDesign(), 80.0), [80.0])["scalars"]
ok("item 4: a static-path g2_op > 1 carries g2_op_valid=False and a reason "
   f"(f8 F_p=1.9: g2_op={_i4s['g2_op']:.4g}; pre-fix f1b artefact: g2_op={_i4bs['g2_op']:.4g}); "
   "an ordinary point stays valid with an empty reason",
   _i4s["g2_op"] > 1.0 and _i4s["g2_op_valid"] is False and "g2_op > 1" in _i4s["g2_op_invalid_reason"]
   and _i4bs["g2_op"] > 1.0 and _i4bs["g2_op_valid"] is False
   and _i4c["g2_op"] <= 1.0 and _i4c["g2_op_valid"] is True and _i4c["g2_op_invalid_reason"] == "")

# ============================================ audit Phase C item 2 (optional Lindblad path)
# DeviceDesign.quantum (default {}, OFF) wires fsim_core.lindblad into
# evaluate(): the cascaded-filter CW g2(0) and the HOM indistinguishability
# at the operating point, as NEW quantum_* scalars only. Expectations are
# the Phase A audit's independent prototype numbers (edge300 exact
# 0.8879 (indep) / 0.9821 (corr) against a memoryless 0.7959; filter-memory
# shift +0.05..+0.08 at 230 K, +0.09..+0.19 at 300 K) and the
# closed form I = Gamma/(Gamma + 2 gamma*) = hbar Gamma / FWHM (Grange et al.,
# PRL 114, 193601 (2015), Eq. 1), never this code's own output.
_HBAR_MEV_NS_V = 6.582119569e-4   # [V] CODATA 2018 hbar (Tiesinga et al., RMP 93, 025010 (2021))


def _canon(x):
    """Canonical, NaN-safe form of a result dict for exact equality."""
    if isinstance(x, dict):
        return {k: _canon(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_canon(v) for v in x]
    if isinstance(x, np.ndarray):
        return ["ndarray", str(x.dtype), [repr(v) for v in x.ravel().tolist()]]
    if isinstance(x, (float, np.floating)):
        return ["float", repr(float(x))]
    return [type(x).__name__, repr(x)]


def _strip_quantum(result):
    out = copy.deepcopy(result)
    sc = out["scalars"]
    for k in [k for k in sc if k.startswith("quantum_")]:
        del sc[k]
    sc["provenance"].pop("quantum", None)
    return out


_qc = DeviceDesign.load(ROOT / "cards" / "edge-inp-gainp-design.yaml")
ok("C2: DeviceDesign.quantum defaults to {} (optional quantum path off)",
   DeviceDesign().quantum == {} and _qc.quantum == {})
_q_off = evaluate(_qc, [300.0])
_qc_on = copy.deepcopy(_qc); _qc_on.quantum = {"enabled": True}
_q_on = evaluate(_qc_on, [300.0])
_qs = _q_on["scalars"]
ok("C2: flag off adds no quantum_* key and no provenance['quantum'] entry (gainp card)",
   not any(k.startswith("quantum_") for k in _q_off["scalars"])
   and "quantum" not in _q_off["scalars"]["provenance"])
ok("C2: flag on changes NO existing key or value: the flag-on result minus its quantum_* "
   "keys equals the flag-off result exactly (curves and scalars, NaN-safe)",
   _canon(_strip_quantum(_q_on)) == _canon(_q_off))
ok("C2: flag on reports the quantum_* keys with a [DR] provenance entry",
   _qs["provenance"]["quantum"]["tag"] == "DR"
   and all(k in _qs for k in ("quantum_g2_cw0_filtered", "quantum_g2_cw0_filtered_dot",
                              "quantum_hom_indistinguishability", "quantum_hom_gamma_ns",
                              "quantum_hom_gamma_star_ns"))
   and _qs["quantum_invalid_reasons"] == [])
ok("C2: quantum filter FWHM is the SAME window the rate path uses (w_resolved)",
   _qs["quantum_filter_fwhm_meV"] == _qs["w_resolved"] and _qs["quantum_linewidth_X_meV"] == _qs["gamma_op"])
_f, _fd = _qs["quantum_g2_cw0_filtered"], _qs["quantum_g2_cw0_filtered_dot"]
ok(f"C2 AC1: gainp 300 K filtered CW g2(0) in [0,1] and above the rate-equation g2_cw0 "
   f"(filter memory fills the dip): filtered {_f:.4f} (dot {_fd:.4f}) vs g2_cw0 {_qs['g2_cw0']:.4f}",
   0.0 <= _f <= 1.0 and 0.0 <= _fd <= 1.0 and _f > _qs["g2_cw0"]
   and _fd > _qs["quantum_g2_cw0_rate_dot"])
ok("C2 AC1: gainp 300 K exact filtered cascade reproduces the audit prototype "
   f"(corr {_qs['quantum_g2_cw0_filtered_dot_corr']:.4f} vs 0.9821, indep "
   f"{_qs['quantum_g2_cw0_filtered_dot_indep']:.4f} vs 0.8879, memoryless Lorentzian "
   f"{_qs['quantum_g2_cw0_rate_dot_lorentzian']:.4f} vs 0.7959; tol 1e-3)",
   abs(_qs["quantum_g2_cw0_filtered_dot_corr"] - 0.9821) <= 1e-3
   and abs(_qs["quantum_g2_cw0_filtered_dot_indep"] - 0.8879) <= 1e-3
   and abs(_qs["quantum_g2_cw0_rate_dot_lorentzian"] - 0.7959) <= 1e-3)
ok("C2 AC1: gainp 300 K filter shift (filtered - memoryless Lorentzian) lies in the audit's "
   f"+0.05..+0.19 range (corr {_qs['quantum_g2_cw0_filter_shift_corr']:+.4f}, indep "
   f"{_qs['quantum_g2_cw0_filter_shift_indep']:+.4f})",
   all(0.05 <= _qs[f"quantum_g2_cw0_filter_shift_{m}"] <= 0.19 for m in ("corr", "indep")))
ok("C2: gainp 300 K correction against the REPORTED top-hat g2_cw0: "
   f"quantum_g2_cw0_shift_vs_reported = {_qs['quantum_g2_cw0_shift_vs_reported']:+.6f} "
   "== filtered - g2_cw0 exactly, pinned at +0.227575 (tol 1e-5), larger than the audit's "
   "memory-only +0.05..+0.19 (it includes the top-hat -> Lorentzian shape change); dot level "
   f"{_qs['quantum_g2_cw0_filtered_dot'] - _qs['quantum_g2_cw0_rate_dot']:+.4f} (about +0.33)",
   _qs["quantum_g2_cw0_shift_vs_reported"] == _qs["quantum_g2_cw0_filtered"] - _qs["g2_cw0"]
   and abs(_qs["quantum_g2_cw0_shift_vs_reported"] - 0.227575) <= 1e-5
   and _qs["quantum_g2_cw0_shift_vs_reported"] > 0.19
   and abs((_qs["quantum_g2_cw0_filtered_dot"] - _qs["quantum_g2_cw0_rate_dot"]) - 0.33) <= 0.01)
# wiring: the measured filtered value uses the SAME background law and aperture
# bath as g2_cw0 (recomputed from exposed scalars through cw_g2/device helpers).
_lam = _qs["aperture_lambda_op"]
_gmix = _compose_aperture_g2(_fd, _lam) if _qs["aperture_compose"] else _fd
_rho_e = devmod._aperture_rho_eff(_qs["cw_rho_op"], _lam) if _qs["aperture_compose"] else _qs["cw_rho_op"]
ok("C2: quantum_g2_cw0_filtered == g2_with_background(aperture-composed filtered dot g2, "
   "the SAME cw_rho_op) (rel 1e-12)",
   abs(float(cw_g2.g2_with_background(_gmix, _rho_e)) / _f - 1.0) <= 1e-12)
_G, _gs, _I = _qs["quantum_hom_gamma_ns"], _qs["quantum_hom_gamma_star_ns"], _qs["quantum_hom_indistinguishability"]
_I_closed = _HBAR_MEV_NS_V * _G / _qs["gamma_op"]
ok(f"C2 AC1: HOM I = {_I:.6f} in [0,1], equals Gamma/(Gamma+2 gamma*) for the reported numbers "
   f"(rel 1e-9) and hbar Gamma/FWHM = {_I_closed:.6f} (rel 1e-9; audit ~0.29 at 300 K)",
   0.0 <= _I <= 1.0 and abs(_I / (_G / (_G + 2.0 * _gs)) - 1.0) <= 1e-9
   and abs(_I / _I_closed - 1.0) <= 1e-9 and abs(_I - 0.29) <= 0.01
   and _G > _qs["cw_gamma_X_ns"])
_qc230 = copy.deepcopy(_qc_on); _qc230.thermal.T_hs = 230.0
_qs230 = evaluate(_qc230, [230.0])["scalars"]
ok("C2: gainp 230 K filter shift lies in the audit's +0.05..+0.08 range "
   f"(corr {_qs230['quantum_g2_cw0_filter_shift_corr']:+.4f}, indep "
   f"{_qs230['quantum_g2_cw0_filter_shift_indep']:+.4f}) and HOM matches hbar Gamma/FWHM (audit ~0.071)",
   all(0.05 <= _qs230[f"quantum_g2_cw0_filter_shift_{m}"] <= 0.08 for m in ("corr", "indep"))
   and abs(_qs230["quantum_hom_indistinguishability"]
           / (_HBAR_MEV_NS_V * _qs230["quantum_hom_gamma_ns"] / _qs230["gamma_op"]) - 1.0) <= 1e-9
   and abs(_qs230["quantum_hom_indistinguishability"] - 0.071) <= 0.005)
_qnocw = copy.deepcopy(_qc_on); _qnocw.drive.cw = False
_qsn = evaluate(_qnocw, [300.0])["scalars"]
ok("C2: flag on without drive.cw reports NaN quantum values with an explicit reason (never a guess)",
   np.isnan(_qsn["quantum_g2_cw0_filtered"]) and np.isnan(_qsn["quantum_hom_indistinguishability"])
   and any("drive.cw" in r for r in _qsn["quantum_invalid_reasons"]))


def _bad_quantum(block):
    def _f():
        dq = copy.deepcopy(_qc); dq.quantum = block
        evaluate(dq, [300.0])
    return _f


raises("C2: unknown quantum key rejected", _bad_quantum({"enabled": True, "filter": 1}))
raises("C2: bad quantum.dephasing rejected", _bad_quantum({"enabled": True, "dephasing": "x"}))
raises("C2: non-bool quantum.enabled rejected", _bad_quantum({"enabled": "yes"}))
_dq = DeviceDesign(); _dq.quantum = {"enabled": True, "dephasing": "indep"}
with tempfile.NamedTemporaryFile(dir=ROOT, suffix=".yaml", delete=False) as tmp:
    _qpath = Path(tmp.name)
try:
    _dq.save(_qpath); _dq2 = DeviceDesign.load(_qpath)
finally:
    os.unlink(_qpath)
ok("C2: YAML round trip preserves the quantum block", _dq2.quantum == _dq.quantum and _dq2 == _dq)
_dq0 = DeviceDesign()
with tempfile.NamedTemporaryFile(dir=ROOT, suffix=".yaml", delete=False) as tmp:
    _qpath0 = Path(tmp.name)
try:
    _dq0.save(_qpath0); _qtext0 = _qpath0.read_text(encoding="utf-8"); _dq0b = DeviceDesign.load(_qpath0)
finally:
    os.unlink(_qpath0)
ok("C2: a flag-off design saves no 'quantum:' key (saved text unchanged) and reloads with quantum == {}",
   "quantum" not in yaml.safe_load(_qtext0)["design"] and _dq0b == _dq0 and _dq0b.quantum == {})

# Planar nitride path (same optional wiring, NEW keys only).
_pn = DeviceDesign.load(ROOT / "cards" / "nitride-cavity-pulse-design.yaml")
_pn_off = evaluate(_pn, [_pn.thermal.T_hs])
_pn_on_d = copy.deepcopy(_pn); _pn_on_d.quantum = {"enabled": True}
_pn_on = evaluate(_pn_on_d, [_pn.thermal.T_hs])
_pq = _pn_on["scalars"]
_pG = _pq["gamma_X_ns"] + _pq["k_X_ns"]
ok("C2: planar nitride flag off adds no quantum_* key; flag on changes no existing key or value",
   not any(k.startswith("quantum_") for k in _pn_off["scalars"])
   and _canon(_strip_quantum(_pn_on)) == _canon(_pn_off))
ok(f"C2: planar nitride pulse card: filtered CW g2(0) {_pq['quantum_g2_cw0_filtered_dot']:.4f} in [0,1] "
   f"above the memoryless rate value {_pq['quantum_g2_cw0_rate_dot']:.4f}; HOM I = "
   f"{_pq['quantum_hom_indistinguishability']:.4f} == Gamma/(Gamma+2 gamma*) == hbar Gamma/FWHM "
   "(Gamma = the row's gamma_X_ns + k_X_ns; rel 1e-9)",
   _pq["valid"] is True and 0.0 <= _pq["quantum_g2_cw0_filtered_dot"] <= 1.0
   and _pq["quantum_g2_cw0_filtered_dot"] > _pq["quantum_g2_cw0_rate_dot"]
   and abs(_pq["quantum_hom_gamma_ns"] / _pG - 1.0) <= 1e-12
   and abs(_pq["quantum_hom_indistinguishability"]
           / (_pG / (_pG + 2.0 * _pq["quantum_hom_gamma_star_ns"])) - 1.0) <= 1e-9
   and abs(_pq["quantum_hom_indistinguishability"]
           / (_HBAR_MEV_NS_V * _pG / _pq["quantum_linewidth_X_meV"]) - 1.0) <= 1e-9)

# =========================================================== quantum tier Q2/Q3/Q4
# Cascade-weighted, filtered and phonon-partitioned HOM.
# Separate block; every check below tests a number the diagnostic did not produce.

# ---- Q2: cascade-weighted HOM (exact identity I_XX->X = I_X Gamma_XX/(Gamma_XX+Gamma_X))
_I_x, _I_xx = _qs["quantum_hom_indistinguishability"], _qs["quantum_hom_indistinguishability_from_XX"]
_cf = _qs["quantum_hom_cascade_factor"]
_Gx = _qs["quantum_hom_gamma_ns"]                       # gamma_X + k_X
_gX_q = _qs["cw_gamma_X_ns"]                            # gamma_X; gamma_XX = 2 gamma_X, k_XX = 2 k_X on this path
_GXX_q = 2.0 * _Gx
ok(f"Q2: gainp 300 K: I(XX->X)/I(X) = {_I_xx / _I_x:.9f} == quantum_hom_cascade_factor = {_cf:.9f} (rel 1e-9)",
   abs(_I_xx / _I_x / _cf - 1.0) <= 1e-9)
ok(f"Q2: cascade factor == Gamma_XX/(Gamma_XX+Gamma_X) from the reported rates "
   f"(Gamma_X = {_Gx:.3f} ns^-1, Gamma_XX = 2 Gamma_X; rel 1e-12)",
   abs(_cf / (_GXX_q / (_GXX_q + _Gx)) - 1.0) <= 1e-12)
_, _P1q, _P2q = loading_probs(_qs["mu_resolved"])
_Iload = (_P1q * _I_x + _P2q * _I_xx) / (_P1q + _P2q)
ok(f"Q2: loaded HOM = {_qs['quantum_hom_indistinguishability_loaded']:.6f} == [P1 I_X + P2 I_XX->X]/(P1+P2) "
   f"with the cap-2 loading_probs(mu_resolved = {_qs['mu_resolved']:.4f}) (rel 1e-12) and lies between "
   "I_XX->X and I_X",
   abs(_qs["quantum_hom_indistinguishability_loaded"] / _Iload - 1.0) <= 1e-12
   and _I_xx < _qs["quantum_hom_indistinguishability_loaded"] < _I_x)
_Ie = _qs["quantum_hom_indistinguishability_loaded_exact"]
_kXq = _Gx - _gX_q
_s3q = lb.build_system(levels=3, gamma_X_ns=_gX_q, gamma_XX_ns=2.0 * _gX_q, k_X=_kXq, k_XX=2.0 * _kXq,
                       deph=[(2.0 * _qs["quantum_hom_gamma_star_ns"], (0.0, 1.0, 0.0))])
_rmix = np.zeros((_s3q.dim, _s3q.dim), dtype=complex)
_rmix[1, 1] = _P1q / (_P1q + _P2q); _rmix[2, 2] = _P2q / (_P1q + _P2q)
_Idir = float(lb.indistinguishability(_s3q.L, _s3q.dim, _rmix, np.sqrt(_gX_q) * _s3q.ops["c_X"]))
ok(f"F-review 5: exact loaded HOM = {_Ie:.6f} <= linear mixture {_qs['quantum_hom_indistinguishability_loaded']:.6f} "
   f"(upper bound) and == a direct lindblad.indistinguishability call on (P1 rho_X + P2 rho_XX)/(P1+P2) (rel 1e-12)",
   np.isfinite(_Ie) and _Ie <= _qs["quantum_hom_indistinguishability_loaded"] * (1.0 + 1e-12)
   and abs(_Ie / _Idir - 1.0) <= 1e-12)
ok(f"Q2: gainp 230 K audit values: I_X = {_qs230['quantum_hom_indistinguishability']:.4f} (0.0712), "
   f"I_XX->X = {_qs230['quantum_hom_indistinguishability_from_XX']:.4f} (0.0475); 300 K: "
   f"{_I_x:.4f} (0.2928), {_I_xx:.4f} (0.1952) (the work order figures are 4-decimal; tol 1e-4)",
   abs(_qs230["quantum_hom_indistinguishability"] - 0.0712) <= 1e-4
   and abs(_qs230["quantum_hom_indistinguishability_from_XX"] - 0.0475) <= 1e-4
   and abs(_I_x - 0.2928) <= 1e-4 and abs(_I_xx - 0.1952) <= 1e-4)
ok("Q2: flag on changes no existing key (re-check with the new keys present: stripped result == flag-off)",
   _canon(_strip_quantum(_q_on)) == _canon(_q_off))

# ---- Q3: filtered HOM
_qw = copy.deepcopy(_qc_on); _qw.filter.auto_w_scale = 1.0e3     # w = 1e3 x line FWHM
_qsw = evaluate(_qw, [300.0])["scalars"]
ok(f"Q3: w = 1e3 x line FWHM ({_qsw['quantum_filter_fwhm_meV']:.0f} meV): filtered HOM "
   f"{_qsw['quantum_hom_indistinguishability_filtered']:.6f} == unfiltered "
   f"{_qsw['quantum_hom_indistinguishability']:.6f} (rel 1e-2)",
   abs(_qsw["quantum_filter_fwhm_meV"] / _qsw["quantum_linewidth_X_meV"] - 1e3) <= 1e-6
   and abs(_qsw["quantum_hom_indistinguishability_filtered"]
           / _qsw["quantum_hom_indistinguishability"] - 1.0) <= 1e-2)
ok(f"Q3: own w: filtered HOM {_qs['quantum_hom_indistinguishability_filtered']:.6f} >= unfiltered "
   f"{_I_x:.6f}; photons per excitation {_qs['quantum_hom_filtered_photons_per_excitation']:.3e} in (0, 1] "
   f"(230 K: {_qs230['quantum_hom_indistinguishability_filtered']:.6f} >= "
   f"{_qs230['quantum_hom_indistinguishability']:.6f})",
   _qs["quantum_hom_indistinguishability_filtered"] >= _I_x
   and 0.0 < _qs["quantum_hom_filtered_photons_per_excitation"] <= 1.0
   and _qs230["quantum_hom_indistinguishability_filtered"] >= _qs230["quantum_hom_indistinguishability"]
   and 0.0 < _qs230["quantum_hom_filtered_photons_per_excitation"] <= 1.0)
_Nmax = _gX_q / _Gx                                      # gamma_X/Gamma_X: wide-filter photon number
ok(f"Q3: wide filter photons per excitation {_qsw['quantum_hom_filtered_photons_per_excitation']:.6e} == "
   f"gamma_X/(gamma_X+k_X) = {_Nmax:.6e} (rel 1e-2) and the own-w filter collects fewer",
   abs(_qsw["quantum_hom_filtered_photons_per_excitation"] / _Nmax - 1.0) <= 1e-2
   and _qs["quantum_hom_filtered_photons_per_excitation"] < _qsw["quantum_hom_filtered_photons_per_excitation"])
_qn = copy.deepcopy(_qc_on); _qn.filter.auto_w_scale = 0.25
_qsn = evaluate(_qn, [300.0])["scalars"]
ok("Q3: trade-off at the operating point: a narrower window (auto_w_scale 0.25 vs 1) raises the "
   "filtered HOM and lowers photons per excitation",
   _qsn["quantum_hom_indistinguishability_filtered"] > _qs["quantum_hom_indistinguishability_filtered"]
   and _qsn["quantum_hom_filtered_photons_per_excitation"]
   < _qs["quantum_hom_filtered_photons_per_excitation"])
ok("Q3: planar nitride pulse card: filtered HOM >= unfiltered and photons per excitation in (0, 1]",
   _pq["quantum_hom_indistinguishability_filtered"] >= _pq["quantum_hom_indistinguishability"]
   and 0.0 < _pq["quantum_hom_filtered_photons_per_excitation"] <= 1.0)

# ---- Q4: consistent ZPL/sideband partition
from fsim_core import qd_gf as _qdgf
_pp = _qdgf.PhononParams(**_qc_on.dot.phonon)
ok(f"Q4: quantum_hom_zpl_weight = {_qs['quantum_hom_zpl_weight']:.4f} == exp(-huang_rhys(PhononParams, T_j)) "
   "at 300 K and 230 K (rel 1e-12; ~0.092 and ~0.161)",
   abs(_qs["quantum_hom_zpl_weight"] / np.exp(-_qdgf.huang_rhys(_pp, _qs["T_j_op"])) - 1.0) <= 1e-12
   and abs(_qs230["quantum_hom_zpl_weight"] / np.exp(-_qdgf.huang_rhys(_pp, _qs230["T_j_op"])) - 1.0) <= 1e-12
   and abs(_qs["quantum_hom_zpl_weight"] - 0.092) <= 5e-4 and abs(_qs230["quantum_hom_zpl_weight"] - 0.161) <= 5e-4)
ok("Q4: default dot.width_kind == total -> phonon HOM is NaN and the reason string is reported "
   "(quantum_hom_phonon_reason; quantum_invalid_reasons stays empty on a valid point)",
   _qc.dot.width_kind == "total" and np.isnan(_qs["quantum_hom_indistinguishability_phonon"])
   and "cannot be partitioned" in _qs["quantum_hom_phonon_reason"]
   and _qs["quantum_invalid_reasons"] == [])
_qz = copy.deepcopy(_qc_on); _qz.dot.width_kind = "zpl"
_qsz = evaluate(_qz, [300.0])["scalars"]
ok(f"Q4: width_kind zpl: phonon HOM {_qsz['quantum_hom_indistinguishability_phonon']:.6f} == Z^2 I_X "
   "(rel 1e-12); the other quantum keys equal the total-width result",
   abs(_qsz["quantum_hom_indistinguishability_phonon"]
       / (_qsz["quantum_hom_zpl_weight"] ** 2 * _qsz["quantum_hom_indistinguishability"]) - 1.0) <= 1e-12
   and _qsz["quantum_hom_phonon_reason"] == ""
   and _qsz["quantum_hom_indistinguishability"] == _I_x
   and _qsz["quantum_g2_cw0_filtered"] == _qs["quantum_g2_cw0_filtered"])
_qz_off = copy.deepcopy(_qc); _qz_off.dot.width_kind = "zpl"
ok("Q4: width_kind does not touch the legacy outputs: zpl flag-off result == default flag-off result",
   _canon(evaluate(_qz_off, [300.0])) == _canon(_q_off))
raises("Q4: unknown dot.width_kind rejected at construction (card parsing)",
       lambda: devmod.DotBlock(width_kind="both"))


def _bad_width_eval():
    dq = copy.deepcopy(_qc_on); dq.dot.width_kind = "both"; evaluate(dq, [300.0])


raises("Q4: unknown dot.width_kind rejected at evaluate (design mutated after construction)", _bad_width_eval)
_wfd, _wname = tempfile.mkstemp(suffix=".yaml"); os.close(_wfd); _wpath = Path(_wname)
try:
    _qz.save(_wpath); _wtext = _wpath.read_text(encoding="utf-8"); _qz2 = DeviceDesign.load(_wpath)
    _qc.save(_wpath); _wtext0 = _wpath.read_text(encoding="utf-8")
finally:
    os.unlink(_wpath)
ok("Q4: YAML round trip keeps width_kind zpl; the default total is not written (saved text unchanged)",
   _qz2.dot.width_kind == "zpl" and "width_kind: zpl" in _wtext and "width_kind" not in _wtext0)
ok("Q4: planar nitride: phonon HOM NaN by default, zpl weight finite in (0, 1]",
   np.isnan(_pq["quantum_hom_indistinguishability_phonon"]) and 0.0 < _pq["quantum_hom_zpl_weight"] <= 1.0)

print(f"{sum(checks)}/{len(checks)} device RT checks passed")
sys.exit(0 if all(checks) else 1)
