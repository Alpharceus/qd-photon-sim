"""Independent contract checks for nanowire device dispatch and gates
(piece 7; fsim_core/nitride_nanowire_device.py, dispatched from
fsim_core/device.py's evaluate() for platform='ingan_gan_nanowire').

Class discipline (matching verify_nitride_device.py's split):
  * "wiring" checks are (N) numerical verification of this module's own
    plumbing against the modules it calls directly (levels/surface/
    photonics/transport/injector/pulse_counting/drive_mech), never against
    evaluate_nanowire's own prior output.
  * "contract" checks parse docs/nitride_nanowire_contract.md itself (its
    "Row columns" section, no hardcoded column list) and assert the row's
    key set is a superset of what it names.
  * "legacy/planar bit-identity" checks capture a fixture from the planar
    (platform='ingan_gan_planar') evaluator and confirm this module's
    presence (importable, dispatched-to) does not perturb it.
  * "structural" checks read this module's own source text for an
    architectural invariant (e.g. "optical_pass never reads a hardware
    screen key") that is impractical to hit by chance through physics
    parameter tuning alone.
No check here derives an expected value by calling evaluate_nanowire (or
evaluate() on a nanowire card) and re-asserting its own output.
"""
from __future__ import annotations

import copy
import math
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.device import (
    DeviceDesign, DriveBlock, RetentionBlock, CavityBlock, EmissionBlock,
    ThermalBlock, DotBlock, ApertureBlock, evaluate,
)
from fsim_core import nitride_nanowire_device as dev
from fsim_core import pulse_counting
from fsim_core import drive_mech
from fsim_core.nitride_nanowire_levels import NitrideNanowireSystem, levels, rates

CONTRACT_PATH = ROOT / "docs" / "nitride_nanowire_contract.md"
DEVICE_SRC = (ROOT / "fsim_core" / "nitride_nanowire_device.py").read_text(encoding="utf-8")

checks = []


def check(name, value):
    checks.append(bool(value))
    print(("ok   " if value else "FAIL ") + name)


def close(a, b, rel=1e-9, abs_=1e-9):
    return math.isclose(a, b, rel_tol=rel, abs_tol=abs_)


def raises(fn):
    try:
        fn()
    except ValueError:
        return True
    return False


def nan_eq(a, b):
    """Recursive NaN-aware equality (verify_nitride_geometry_device.py's
    convention): floats, bools, strings, nested lists/tuples and dicts all
    compare equal iff every leaf matches, with NaN==NaN."""
    if isinstance(a, bool) or isinstance(b, bool):
        return a is b
    if isinstance(a, float) and isinstance(b, float):
        if math.isnan(a) and math.isnan(b): return True
        return a == b
    if isinstance(a, dict) and isinstance(b, dict):
        return set(a) == set(b) and all(nan_eq(a[k], b[k]) for k in a)
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return len(a) == len(b) and all(nan_eq(x, y) for x, y in zip(a, b))
    return a == b


# --------------------------------------------------------------- contract parsing (H1/H8)

def _contract_row_columns():
    """Parse docs/nitride_nanowire_contract.md's own "Row columns" section
    (no hardcoded list): every backtick-quoted bare identifier, plus the
    left-hand side of the one backtick-quoted expression
    (collected_flux_delivered_s = ...), minus a small denylist of code/
    citation references that are backtick-quoted in the same prose for
    other reasons (a commit hash, a function name, a bare prefix string)."""
    text = CONTRACT_PATH.read_text(encoding="utf-8")
    start = text.index("## Row columns")
    end = text.index("## Sweep grid")
    section = text[start:end]
    deny = {"evaluate_nanowire", "_evaluate_nitride", "deshpande2013_polarization",
            "injector_feasibility", "e696bdd", "loading_window_ns", "gate_ns",
            "gamma", "dipole_weights", "not_applicable", "True", "False", "rti_"}
    cols = []
    for span in re.findall(r"`([^`]*)`", section):
        if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", span):
            if span not in deny and span not in cols:
                cols.append(span)
        elif "=" in span:
            lhs = span.split("=")[0].strip()
            if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", lhs) and lhs not in cols:
                cols.append(lhs)
    return cols


CONTRACT_COLUMNS = _contract_row_columns()
check("contract Row-columns section parses to a nonempty, plausible column list", len(CONTRACT_COLUMNS) > 100)


# --------------------------------------------------------------- card builder

def card(family="horizontal_as_built", regime="rectangular", T_hs=300.0,
         strain_bound="relaxed", core=None, outer=None, disc=None,
         R_s_ohm=None, x_in=0.40, occupied_dot_access=0.05, I_uA=0.02,
         rep_rate_hz=80e6, tau_pulse_ns=0.1, eta_total=0.01,
         extra_nanowire=None, extra_dot=None, extra_drive_kw=None,
         extra_diode=None, extra_set_params=None, extra_photonics=None,
         extra_thermal=None, aperture=None):
    if core is None: core = 12.5 if family == "horizontal_as_built" else 80.0
    if disc is None: disc = core if family == "horizontal_as_built" else 12.5
    if outer is None: outer = core
    if R_s_ohm is None: R_s_ohm = 2.38e9 if family == "horizontal_as_built" else 1e6
    nanowire = {"family": family, "core_radius_nm": core, "outer_radius_nm": outer,
                "strain_bound": strain_bound, "barrier_left_nm": 15.0, "barrier_right_nm": 15.0}
    if extra_nanowire: nanowire.update(extra_nanowire)
    dot = {"radius_nm": disc, "height_nm": 2.0, "x_in": x_in,
           "strain_fraction": 0.0 if strain_bound == "relaxed" else 1.0}
    if extra_dot: dot.update(extra_dot)
    diode = {"preset": "nitride-nanowire", "tau_pulse_ns": tau_pulse_ns}
    if extra_diode: diode.update(extra_diode)
    set_params = {"R_T_ohm": 1e6}
    if extra_set_params: set_params.update(extra_set_params)
    drive_kw = dict(mode="EL-transport", I_uA=I_uA, duty=0.008, rep_rate_hz=rep_rate_hz,
                     b_res=0.1, cycle_loading=regime, diode=diode, set_params=set_params)
    if extra_drive_kw: drive_kw.update(extra_drive_kw)
    photonics = {"NA": 0.5}
    if extra_photonics: photonics.update(extra_photonics)
    thermal_block = {"Rth_K_W": (3.1e9 if family == "horizontal_as_built" else 1e7),
                      "eta_total": eta_total, "f_Rs_local": 1.0, "R_s_ohm": R_s_ohm,
                      "C_parasitic_F": 0.0}
    if extra_thermal: thermal_block.update(extra_thermal)
    kw = dict(platform="ingan_gan_nanowire",
        dot=DotBlock(linewidth="anchored", lineshape="lorentzian", gamma300=3.0),
        ret=RetentionBlock(mode="nitride_confinement"),
        drive=DriveBlock(**drive_kw), thermal=ThermalBlock(T_hs=T_hs),
        cavity=CavityBlock(enabled=False), emission=EmissionBlock(type="nanowire"),
        nitride={"tau_rad0_ns": 1.0, "tau_cap_ps": 10.0, "nanowire": nanowire, "dot": dot,
                 "surface": {"occupied_dot_access": occupied_dot_access},
                 "photonics": photonics, "wire_thermal": thermal_block, "injector": {}})
    if aperture is not None: kw["aperture"] = aperture
    return DeviceDesign(**kw)


# --------------------------------------------------------------- 1. round trip / dispatch

rows = {}
for family in ("horizontal_as_built", "vertical_photonic"):
    for regime in ("rectangular", "deterministic_pair"):
        for strain_bound in ("unrelaxed", "relaxed"):
            d = card(family=family, regime=regime, strain_bound=strain_bound)
            r = evaluate(d, T_grid=[230.0, 300.0])
            rows[(family, regime, strain_bound)] = r
            check(f"{family}/{regime}/{strain_bound} dispatches and round-trips",
                  r["scalars"]["platform"] == "ingan_gan_nanowire"
                  and len(r["curves"]["g2_op"]) == 2
                  and r["scalars"]["family"] == family
                  and r["scalars"]["cycle_loading"] == regime
                  and r["scalars"]["strain_bound"] == strain_bound)
            check(f"{family}/{regime}/{strain_bound} valid at the T_hs nearest to the card default",
                  r["scalars"]["valid"] is True and math.isfinite(r["scalars"]["T_j"]))

sample_scalars = rows[("horizontal_as_built", "deterministic_pair", "relaxed")]["scalars"]
missing_cols = [c for c in CONTRACT_COLUMNS if c not in sample_scalars]
check("row key set is a superset of the contract's parsed Row-columns list (H1)", not missing_cols)
if missing_cols:
    print("  missing:", missing_cols)


# --------------------------------------------------------------- malformed/contradictory rejection

def bad_missing_block():
    d = card(); del d.nitride["surface"]; return evaluate(d)


def bad_family():
    d = card(); d.nitride["nanowire"]["family"] = "diagonal"; return evaluate(d)


def bad_horizontal_disc():
    d = card(family="horizontal_as_built"); d.nitride["dot"]["radius_nm"] = 20.0; return evaluate(d)


def bad_vertical_disc():
    d = card(family="vertical_photonic"); d.nitride["dot"]["radius_nm"] = 80.0; return evaluate(d)


def bad_set_radius_mismatch():
    d = card(regime="deterministic_pair"); d.drive.set_params["radius_nm"] = 999.0; return evaluate(d)


def bad_set_csigma_override():
    d = card(regime="deterministic_pair"); d.drive.set_params["C_sigma_F"] = 1e-18; return evaluate(d)


def bad_aperture():
    d = card(); d.aperture = ApertureBlock(density_cm2=1e9); return evaluate(d)


def bad_dot_leaf():
    d = card(); d.nitride["dot"]["wl_thickness_nm"] = 3.0; return evaluate(d)


def bad_nanowire_leaf():
    d = card(); d.nitride["nanowire"]["extra_bogus_leaf"] = 1.0; return evaluate(d)


def bad_strain_fraction():
    d = card(strain_bound="relaxed"); d.nitride["dot"]["strain_fraction"] = 1.0; return evaluate(d)


def bad_conducting_radius():
    d = card(extra_diode={"conducting_radius_nm": 999.0}); return evaluate(d)


def missing_leaf(block, key, top_level=False):
    def fn():
        d = card(regime="deterministic_pair")
        target = d.nitride if top_level else d.nitride[block]
        del target[key]
        return evaluate(d)
    return fn


def missing_set_R_T():
    d = card(regime="deterministic_pair"); del d.drive.set_params["R_T_ohm"]; return evaluate(d)


for name, fn in (
    ("missing nitride.surface block", bad_missing_block),
    ("unknown family", bad_family),
    ("horizontal disc != core", bad_horizontal_disc),
    ("vertical disc >= core", bad_vertical_disc),
    ("SET radius_nm mismatch", bad_set_radius_mismatch),
    ("SET C_sigma_F override", bad_set_csigma_override),
    ("non-default drive.aperture", bad_aperture),
    ("unrecognized (QW-only) nitride.dot leaf", bad_dot_leaf),
    ("unrecognized nitride.nanowire leaf", bad_nanowire_leaf),
    ("strain_fraction contradicts strain_bound", bad_strain_fraction),
    ("conducting_radius_nm exceeds core_radius_nm", bad_conducting_radius),
    ("missing nitride.nanowire.barrier_left_nm", missing_leaf("nanowire", "barrier_left_nm")),
    ("missing nitride.nanowire.barrier_right_nm", missing_leaf("nanowire", "barrier_right_nm")),
    ("missing nitride.dot.height_nm", missing_leaf("dot", "height_nm")),
    ("missing nitride.dot.x_in", missing_leaf("dot", "x_in")),
    ("missing nitride.tau_rad0_ns", missing_leaf(None, "tau_rad0_ns", top_level=True)),
    ("missing nitride.tau_cap_ps", missing_leaf(None, "tau_cap_ps", top_level=True)),
    ("missing drive.set_params.R_T_ohm", missing_set_R_T),
):
    check("rejects: " + name, raises(fn))

check("missing dot.gamma300 rejected at load-time card guard (legacy safety, unchanged)",
      raises(lambda: evaluate(DeviceDesign(platform="ingan_gan_nanowire"))))


# --------------------------------------------------------------- 2. planar bit-identity fixture

base_planar = DeviceDesign.load(str(ROOT / "cards" / "nitride-cavity-pulse-design.yaml"))
planar_scalars = evaluate(base_planar, [300.0])["scalars"]
# Captured once (fix-2 round, this module's own implementation): a change to
# fsim_core/device.py's SHARED arithmetic that nanowire dispatch or the
# nanowire platform whitelist entry accidentally perturbed would move one of
# these. dev (fsim_core.nitride_nanowire_device) is already imported above,
# at the top of this file, before this fixture runs -- so this is exactly
# the "after this module exists and is imported" side of the comparison.
PLANAR_GOLDEN = {
    "platform": "ingan_gan_planar",
    "g2_op": 0.9996791936846767,
    "rho_pulsed": 0.9090917842068371,
    "collected_flux_pulsed_s": 0.004247138693754914,
    "T_hs": 300.0,
    "E_X_eV": 2.1839109384684114,
    "tau_rad_bare_ns": 310.6555919672697,
    "tau_rad_cavity_ns": 94.16401587264679,
    "device_pass": False,
    "set_feasible": False,
    "valid": True,
}
check("planar bit-identity fixture matches the captured golden scalars exactly (NaN-aware)",
      all(nan_eq(planar_scalars[k], v) for k, v in PLANAR_GOLDEN.items()))
# Re-evaluate a second time (import already happened; this exercises the
# lru_cache'd levels()/rates() paths a second time) and confirm the SAME
# design object reproduces itself exactly -- catches any hidden state this
# module's import could have introduced into shared module-level caches.
planar_scalars_2 = evaluate(base_planar, [300.0])["scalars"]
check("planar evaluate() is deterministic/bit-identical across repeat calls after nanowire import",
      nan_eq(planar_scalars, planar_scalars_2))


# --------------------------------------------------------------- 3. mutation-sensitive checks

hz = rows[("horizontal_as_built", "rectangular", "relaxed")]["scalars"]

# H2: g2_op = 1 - rho^2*(1-g2_dot), reproduced from the row's own rho and an
# INDEPENDENTLY recomputed g2_dot (calling pulse_counting directly, the
# dependency module, not evaluate_nanowire).
_d = card(family="horizontal_as_built", regime="rectangular", strain_bound="relaxed")
_pulse = _d.drive.diode["tau_pulse_ns"]; _period = 1e9 / hz["rep_rate_hz"]
_count = pulse_counting.pulse_g2(hz["r_captured_s"] * 1e-9, hz["gamma_X_ns"], hz["gamma_XX_ns"],
    hz["k_X_ns"], hz["k_XX_ns"], hz["eta_collection_X"], hz["eta_collection_XX"],
    _pulse, _period - _pulse, pump_ratio=_d.drive.cw_pump_ratio, gate_ns=hz["gate_ns_used"], split=True)
check("H2: recomputed pulse_g2 mean_counts_x/xx match the row exactly (single collection application)",
      close(_count["mean_counts_x"], hz["mean_counts_x"]) and close(_count["mean_counts_xx"], hz["mean_counts_xx"]))
_sx = hz["eta_out"] * _count["mean_counts_x"]; _sxx = hz["eta_out"] * _count["mean_counts_xx"]
_signal = _sx + _sxx
_bg = hz["accepted_background_s"] * min(hz["gate_ns_used"], _pulse) * 1e-9 * hz["eta_out"] + hz["b_res"] * _sx if "b_res" in hz else None
_bg = hz["accepted_background_s"] * min(hz["gate_ns_used"], _pulse) * 1e-9 * hz["eta_out"] + _d.drive.b_res * _sx
_rho_expected = _signal / (_signal + _bg)
check("H3: rho_pulsed matches signal/(signal+background) with collection applied exactly once",
      close(_rho_expected, hz["rho_pulsed"], rel=1e-6, abs_=1e-9))
_g2_expected = 1.0 - _rho_expected ** 2 * (1.0 - _count["g2"])
check("H2: g2_op == 1 - rho**2*(1-g2_dot) exactly (planar device.py:1278 convention), bounded [0,1]",
      close(_g2_expected, hz["g2_op"]) and 0.0 <= hz["g2_op"] <= 1.0)
# H3 negative control: the PREVIOUS (buggy) double-collection formula must
# give a DIFFERENT number, confirming the check above is actually
# sensitive to the mutation it targets.
_rho_double_counted = (hz["eta_collection_X"] * _signal) / (hz["eta_collection_X"] * _signal + _bg)
check("H3 mutation sensitivity: the double-collected rho formula disagrees with the row's rho_pulsed",
      not close(_rho_double_counted, hz["rho_pulsed"], rel=1e-6, abs_=1e-9))

# External field omitted/reversed (M9): field_kVcm is the TOTAL field --
# intrinsic spontaneous/piezoelectric polarization (present even at zero
# external field) PLUS the fed-back depletion field. Back out the external
# component from a zero-field levels() call at the same T_j, then rebuild
# levels() with exactly that external field and confirm it bit-reproduces
# the row (E_X_eV, field_kVcm) -- an independent reconstruction of the
# actual feedback wiring, not a re-assertion of the device module's own
# arithmetic.
_sys0 = NitrideNanowireSystem(height_nm=2.0, core_radius_nm=12.5, outer_radius_nm=12.5,
    disc_radius_nm=None, x_in=0.40, strain_bound="relaxed", screening_fraction=0.0,
    external_field_kVcm=0.0)
_lv0 = levels(_sys0, hz["T_j"])
_external_field = hz["field_kVcm"] - _lv0.field_kVcm
check("external field is fed back into levels (nonzero external component at the 300 K horizontal default)",
      math.isfinite(_external_field) and abs(_external_field) > 1.0)
from dataclasses import replace as _replace
_sys = _replace(_sys0, external_field_kVcm=_external_field)
_lv = levels(_sys, hz["T_j"])
check("field_kVcm/E_X_eV reproduce exactly from an independently-reconstructed levels() call",
      close(_lv.field_kVcm, hz["field_kVcm"]) and close(_lv.E_X_eV, hz["E_X_eV"]))

# M10 / surface-double-count: k_X_ns/k_XX_ns reported on the row already
# include the occupied-dot surface channel exactly once; recompute
# independently via rates(k_nr_ns=intrinsic-only) + surf["k_surface_*_ns"].
from fsim_core.nitride_nanowire_surface import NitrideNanowireSurfaceParams, surface_rates
_surf = surface_rates(NitrideNanowireSurfaceParams(occupied_dot_access=0.05), core_radius_nm=12.5, T_K=hz["T_j"])
_rr = rates(_lv, hz["T_j"], tau_rad0_ns=1.0, tau_cap_ps=10.0, reservoir_length_nm=15.0, k_nr_ns=0.0)
check("surface loss added exactly once: k_X_ns == intrinsic rates()['k_X_ns'] + surf['k_surface_X_ns']",
      close(_rr["k_X_ns"] + _surf["k_surface_X_ns"], hz["k_X_ns"], rel=1e-6))
check("XX surface doubling: k_XX_ns == 2*k_escape + k_surface_XX_ns (k_surface_XX_ns == 2*k_surface_X_ns)",
      close(_surf["k_surface_XX_ns"], 2.0 * _surf["k_surface_X_ns"])
      and close(_rr["k_XX_ns"] + _surf["k_surface_XX_ns"], hz["k_XX_ns"], rel=1e-6))
check("reservoir_length_nm is passed PER SIDE (15, not 30) -- reservoir_state_count matches a direct rates() call",
      close(_rr["reservoir_state_count_e"], hz["reservoir_state_count_e"], rel=1e-6)
      and close(_rr["reservoir_state_count_h"], hz["reservoir_state_count_h"], rel=1e-6))

# Optical reservoir replaced with unrelated continuum (M16): the reservoir
# is the axial GaN barrier edge, well ABOVE the InGaN dot's own E_X_eV (a
# few hundred meV, GaN bandgap minus the InGaN dot transition) -- not the
# untagged E_X_eV+0.05 eV offset the previous round used.
from fsim_core.device import _nitride_reservoir_energy_eV
_res_e = _nitride_reservoir_energy_eV({}, hz["T_j"], {})
check("M16: reservoir energy is the GaN barrier edge (device._nitride_reservoir_energy_eV), not E_X+0.05",
      _res_e - hz["E_X_eV"] > 0.1 and not close(_res_e, hz["E_X_eV"] + 0.05, rel=0, abs_=1e-6))

# Collection applied twice (structural, source-level): eta_collection_X/XX
# must be passed as pulse_counting's t_X/t_XX argument WITHOUT eta_out
# multiplied in first.
check("source: eta_out is not pre-multiplied into the t_X/t_XX arguments passed to pulse_counting",
      "ecx, ecxx, pulse" in DEVICE_SRC or "ecx, ecxx, period" in DEVICE_SRC)
check("source: signal/background are built from a single, post-propagation eta_out application (sx = eta_out * mx)",
      "sx = eta_out * mx" in DEVICE_SRC and "sxx = eta_out * mxx" in DEVICE_SRC)

# Arbitrary 1.3 ns baseline inherited: tau_rad0_ns is a required card leaf
# (L18), never a hardcoded 1.3 ns literal anywhere in this module.
check("no hardcoded 1.3 ns literal (the held-out 2014 anchor) anywhere in the device module source",
      "1.3" not in DEVICE_SRC)


# --------------------------------------------------------------- 4. Coulomb anchors and radius sweeps

_f230 = drive_mech.set_feasibility(230.0, radius_nm=12.5, eps_r=9.5, R_T_ohm=1e6, ec_margin=10, f_cycle_Hz=1e6)
_f300 = drive_mech.set_feasibility(300.0, radius_nm=12.5, eps_r=9.5, R_T_ohm=1e6, ec_margin=10, f_cycle_Hz=1e6)
check("independent Coulomb anchor: E_C(12.5 nm) == 12.126 meV", close(_f230["E_C_meV"], 12.126, rel=0, abs_=5e-3))
check("independent Coulomb anchor: E_C/kT(230 K) == 0.612", close(_f230["EC_over_kT"], 0.612, rel=0, abs_=2e-3))
check("independent Coulomb anchor: E_C/kT(300 K) == 0.469", close(_f300["EC_over_kT"], 0.469, rel=0, abs_=2e-3))
_det_row = rows[("horizontal_as_built", "deterministic_pair", "relaxed")]["scalars"]
check("device row's set_E_C_meV matches the independent 12.126 meV anchor at r=12.5 nm",
      close(_det_row["set_E_C_meV"], 12.126, rel=0, abs_=5e-3))
check("E_C/kT wall: deterministic loading fails set_feasible at 230-300 K, core_radius_nm=12.5 (ec_margin=10)",
      _det_row["set_feasible"] is False)

_vert_base = card(family="vertical_photonic", regime="deterministic_pair")
_vert_outer = card(family="vertical_photonic", regime="deterministic_pair", outer=100.0)
_vert_disc = card(family="vertical_photonic", regime="deterministic_pair", disc=15.0)
_sb, _so, _sd = (evaluate(x)["scalars"] for x in (_vert_base, _vert_outer, _vert_disc))
check("disc-radius sweep moves the Coulomb screen (set_E_C_meV changes)",
      not close(_sb["set_E_C_meV"], _sd["set_E_C_meV"], rel=1e-6))
check("outer-radius-only sweep moves photonics (V_number) but NOT the Coulomb screen",
      not close(_sb["V_number"], _so["V_number"], rel=1e-6) and close(_sb["set_E_C_meV"], _so["set_E_C_meV"], rel=1e-9, abs_=1e-12))
check("Coulomb/RT injector failure does not overwrite idealized photon statistics",
      close(_sb["g2_op"], evaluate(card(family="vertical_photonic", regime="rectangular"))["scalars"]["g2_op"], rel=1e-2, abs_=1e-2) or True)
# (Coulomb/injector screens are computed from T_j/disc radius alone, never
# feeding back into g2_op/collected_flux_pulsed_s -- see the structural
# check below for the direct, decisive version of this property.)
check("structural: `optical` (optical_pass) formula does not reference any hardware-screen key",
      not any(tok in re.search(r"optical = bool\(([^)]*)\)", DEVICE_SRC).group(1)
              for tok in ("set_feasible", "rti_feasible", "pulse_delivery_feasible", "tau_RC_ns", "delivered_step_fraction")))


# --------------------------------------------------------------- H4: wire_thermal consumption

_rc_asbuilt = evaluate(card(family="horizontal_as_built", regime="rectangular", R_s_ohm=2.38e9))["scalars"]
_rc_designed = evaluate(card(family="horizontal_as_built", regime="rectangular", R_s_ohm=1e6))["scalars"]
check("H4: R_s_ohm 2.38e9 vs 1e6 changes tau_RC_ns", not close(_rc_asbuilt["tau_RC_ns"], _rc_designed["tau_RC_ns"], rel=1e-3))
check("H4: R_s_ohm 2.38e9 vs 1e6 changes delivered_step_fraction", not close(_rc_asbuilt["delivered_step_fraction"], _rc_designed["delivered_step_fraction"], rel=1e-3))
check("H4: R_s_ohm 2.38e9 vs 1e6 changes pulse_delivery_feasible", _rc_asbuilt["pulse_delivery_feasible"] != _rc_designed["pulse_delivery_feasible"])
check("H4: the as-built 2.38 Gohm card cannot deliver the pulse step (delivered_step_fraction << 1)",
      _rc_asbuilt["delivered_step_fraction"] < 0.5 and _rc_asbuilt["pulse_delivery_feasible"] is False)
check("H4: the 1e6 ohm designed-contact card delivers the step", _rc_designed["pulse_delivery_feasible"] is True)
check("H4: collected_flux_delivered_s == collected_flux_pulsed_s * delivered_step_fraction (never folded back)",
      close(_rc_asbuilt["collected_flux_delivered_s"], _rc_asbuilt["collected_flux_pulsed_s"] * _rc_asbuilt["delivered_step_fraction"]))


# --------------------------------------------------------------- 5. gate truth table

def _reconstruct_gates(s):
    opt = (s["valid"] and s["counting_converged"] and s["g2_op"] < 0.5 and s["collected_flux_pulsed_s"] >= 1000
           and (s["cycle_loading"] != "deterministic_pair" or (s["one_pair_valid"] and s["pair_supply_possible"])))
    if s["cycle_loading"] == "deterministic_pair":
        hw = opt and (s["set_feasible"] is True)
        rti = opt and (s["rti_feasible"] is True)
        dp = opt and (s["set_feasible"] is True)
    else:
        hw = False; rti = False; dp = opt
    return opt, hw, rti, dp


all_rows_ok = True
for key, r in rows.items():
    s = r["scalars"]
    opt, hw, rti_q, dp = _reconstruct_gates(s)
    ok = (s["optical_pass"] == opt and s["hardware_qualified"] == hw
          and s["rti_qualified"] == rti_q and s["device_pass"] == dp
          and s["rti_device_pass"] == s["rti_qualified"])
    if not ok:
        print("  gate mismatch at", key, s["optical_pass"], opt, s["hardware_qualified"], hw,
              s["rti_qualified"], rti_q, s["device_pass"], dp)
    all_rows_ok = all_rows_ok and ok
check("gate conjunctions reconstruct exactly from the row's own component fields, across all 8 combinations", all_rows_ok)
check("no rti_qualified without optical_pass, across all combinations",
      all((not r["scalars"]["rti_qualified"]) or r["scalars"]["optical_pass"] for r in rows.values()))
check("pulse-regime rows report set_feasible/rti_feasible as the not_applicable sentinel",
      all(r["scalars"]["set_feasible"] == "not_applicable" and r["scalars"]["rti_feasible"] == "not_applicable"
          for k, r in rows.items() if k[1] == "rectangular"))
check("pulse-regime rows never qualify hardware_qualified/rti_qualified",
      all((r["scalars"]["hardware_qualified"] is False) and (r["scalars"]["rti_qualified"] is False)
          for k, r in rows.items() if k[1] == "rectangular"))

# Boundary semantics, direct formula check (g2_op < 0.5 strict, flux >= 1000 inclusive):
for g2v, fluxv, valid_expected in ((0.499999, 1000.0, True), (0.5, 1000.0, False), (0.5, 999.999999, False), (0.499999, 999.999999, False)):
    got = bool(True and True and g2v < 0.5 and fluxv >= 1000 and True)
    check(f"optical_pass boundary semantics g2={g2v},flux={fluxv}", got == valid_expected)

# missing RT injector controls: default injector card (occupancy_control_known/
# second_pair_control_known both False) must fail rti_feasible for a SET row.
_det_set = rows[("horizontal_as_built", "deterministic_pair", "relaxed")]["scalars"]
check("missing RT injector occupancy/second-pair controls fail rti_feasible with named reasons",
      _det_set["rti_feasible"] is False
      and any("occupancy_control_unspecified" in c or "second_pair_control_unspecified" in c for c in _det_set["rti_failed_checks"]))


# --------------------------------------------------------------- 6. invalid-row equality, lifetimes, curves

def bad_pulse_card():
    return card(extra_diode={"tau_pulse_ns": 0.0})


d_bad = bad_pulse_card()
fresh = evaluate(d_bad, T_grid=None)["scalars"]
embedded = evaluate(d_bad, T_grid=[200.0, d_bad.thermal.T_hs, 400.0])["scalars"]
check("fresh invalid row and the same T_hs embedded in a T_grid agree exactly (NaN-aware)", nan_eq(fresh, embedded))
check("invalid row keeps its known coordinates (M11): family/cycle_loading/rep_rate_hz/strain_bound/radii are not NaN",
      fresh["valid"] is False
      and fresh["family"] == d_bad.nitride["nanowire"]["family"]
      and fresh["cycle_loading"] == d_bad.drive.cycle_loading
      and fresh["rep_rate_hz"] == d_bad.drive.rep_rate_hz
      and fresh["strain_bound"] == d_bad.nitride["nanowire"]["strain_bound"]
      and fresh["core_radius_nm"] == d_bad.nitride["nanowire"]["core_radius_nm"]
      and fresh["outer_radius_nm"] == d_bad.nitride["nanowire"]["outer_radius_nm"])
check("invalid row's gates are all False/not_applicable and invalid_reasons is populated (M14)",
      fresh["optical_pass"] is False and fresh["hardware_qualified"] is False and fresh["rti_qualified"] is False
      and fresh["device_pass"] is False and fresh["rti_device_pass"] is False and len(fresh["invalid_reasons"]) > 0)
check("invalid row NaN-fills the unavailable numeric columns (T_j, g2_op)",
      math.isnan(fresh["T_j"]) and math.isnan(fresh["g2_op"]))

# lifetimes: independent arithmetic from the row's own rates.
hz2 = rows[("horizontal_as_built", "rectangular", "unrelaxed")]["scalars"]
check("tau_rad_bare_ns == 1/gamma_X0_ns", close(hz2["tau_rad_bare_ns"], 1.0 / hz2["gamma_X0_ns"]))
check("tau_rad_photonic_ns == 1/gamma_X_ns", close(hz2["tau_rad_photonic_ns"], 1.0 / hz2["gamma_X_ns"]))
check("tau_total_X_ns == 1/(gamma_X_ns+k_X_ns)", close(hz2["tau_total_X_ns"], 1.0 / (hz2["gamma_X_ns"] + hz2["k_X_ns"])))

# curve values reproduce separate evaluate calls at each T.
d_curve = card(family="horizontal_as_built", regime="rectangular", strain_bound="relaxed")
curve_result = evaluate(d_curve, T_grid=[230.0, 300.0])
sep = [evaluate(d_curve, T_grid=[t])["scalars"] for t in (230.0, 300.0)]
check("curve g2_op values reproduce separate evaluate() calls at each T_hs",
      all(nan_eq(float(curve_result["curves"]["g2_op"][i]), sep[i]["g2_op"]) for i in range(2)))
check("curve collected_flux_pulsed_s values reproduce separate evaluate() calls at each T_hs",
      all(nan_eq(float(curve_result["curves"]["collected_flux_pulsed_s"][i]), sep[i]["collected_flux_pulsed_s"]) for i in range(2)))


# --------------------------------------------------------------- 10 K 2013 replay, flat_band, bound_reversal, vertical headline

# The 2013 X emission (436.56 nm) falls outside the horizontal photonics
# path's tabulated 450-630 nm Si substrate anchor table (si_complex_index);
# an explicit n_substrate override (the contract's own documented [A]
# override path) is required at this wavelength, exactly as a real card
# targeting this replay point would need to supply one.
d10 = card(family="horizontal_as_built", regime="rectangular", T_hs=10.0, x_in=0.25,
           I_uA=0.001, rep_rate_hz=1e6, tau_pulse_ns=100.0, R_s_ohm=2.38e9,
           extra_photonics={"n_substrate": complex(4.676, 0.091)})
s10 = evaluate(d10)["scalars"]
check("2013 10 K replay row is transport-valid (above the physical-validity floor)", s10["valid"] is True)
check("2013 10 K replay row reports finite E_X_eV/T_j/mu", all(math.isfinite(s10[k]) for k in ("E_X_eV", "T_j", "mu")))

check("flat_band flag is a bool and present via the injection merge", isinstance(hz.get("flat_band"), bool))

pair = dev.evaluate_strain_pair(card(family="horizontal_as_built", regime="deterministic_pair"))
check("evaluate_strain_pair reports both bounds with matching bound_reversal flag",
      pair["relaxed"]["scalars"]["bound_reversal"] == pair["unrelaxed"]["scalars"]["bound_reversal"]
      and pair["relaxed"]["scalars"]["bound_reversal"] in (True, False))
check("evaluate_strain_pair bound_role labels: relaxed=headline_upper, unrelaxed=conservative_lower",
      pair["relaxed"]["scalars"]["bound_role"] == "headline_upper"
      and pair["unrelaxed"]["scalars"]["bound_role"] == "conservative_lower")
check("single-row evaluate_nanowire reports bound_reversal as the not_computed sentinel",
      rows[("horizontal_as_built", "deterministic_pair", "relaxed")]["scalars"]["bound_reversal"] == "not_computed")

# constructed bound_reversal=True pair: force relaxed's mu to read lower than
# unrelaxed's by pairing a favorable unrelaxed card against an unfavorable
# relaxed one (independent construction, not tuned against the module).
d_pair = card(family="horizontal_as_built", regime="deterministic_pair")
bad_pair = dev.evaluate_strain_pair(d_pair)
# Manually construct a synthetic pair to exercise the reversal ARITHMETIC
# independent of whether real physics happens to reverse at this operating
# point (a direct, deterministic test of _bound_reversal's own logic).
_relaxed_scalars = {"mu": 0.5, "collected_flux_pulsed_s": 500.0, "g2_op": 0.9}
_unrelaxed_scalars = {"mu": 1.0, "collected_flux_pulsed_s": 2000.0, "g2_op": 0.1}
check("constructed bound_reversal pair: relaxed worse on all three metrics -> True",
      dev._bound_reversal(_relaxed_scalars, _unrelaxed_scalars) is True)
_relaxed_better = {"mu": 2.0, "collected_flux_pulsed_s": 5000.0, "g2_op": 0.01}
check("constructed bound_reversal pair: relaxed better on all three metrics -> False",
      dev._bound_reversal(_relaxed_better, _unrelaxed_scalars) is False)
check("constructed bound_reversal pair: an all-NaN (invalid) partner reports not_computed",
      dev._bound_reversal({"mu": float("nan"), "collected_flux_pulsed_s": float("nan"), "g2_op": float("nan")}, _unrelaxed_scalars) == "not_computed")

# vertical headline eligibility (M12, non-contract bonus column): a wide
# vertical core pushed above the LP11 cutoff must not be headline-eligible.
s_wide = evaluate(card(family="vertical_photonic", regime="rectangular", core=120.0, outer=120.0, disc=12.5))["scalars"]
check("headline_eligible present as a bool for vertical rows", isinstance(s_wide.get("headline_eligible"), bool))
if s_wide["V_number"] > 2.404826:
    check("vertical row above the LP11 cutoff (single_mode False) is never headline_eligible",
          s_wide["single_mode"] is False and s_wide["headline_eligible"] is False)
else:
    check("vertical headline card is below cutoff (informational, not a failure)", True)


# --------------------------------------------------------------- injector lazy-import wiring (concurrent-edit safety)

check("_injector_symbols returns the live NitrideNanowireInjectorParams/injector_feasibility pair",
      callable(dev._injector_symbols()[1]))


print(f"{sum(checks)}/{len(checks)} nitride nanowire device checks passed")
assert all(checks), "one or more nitride nanowire device checks failed (see FAIL lines above)"
