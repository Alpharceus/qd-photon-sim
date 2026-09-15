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


def bad_shell_none_outer_mismatch():
    d = card(outer=20.0); d.nitride["nanowire"]["shell"] = "none"; return evaluate(d)


def bad_shell_algan_thin():
    d = card(outer=13.0); d.nitride["nanowire"]["shell"] = "AlGaN"; return evaluate(d)


def bad_shell_contradiction():
    d = card(outer=15.5); d.nitride["nanowire"]["shell"] = "AlGaN"; d.nitride["surface"]["shell"] = "none"; return evaluate(d)


def bad_eta_load():
    d = card(); d.drive.eta_load = 0.5; return evaluate(d)


def bad_R_s_contradiction():
    d = card(extra_diode={"R_s_ohm": 1.0}); return evaluate(d)


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
    ("HIGH 3: shell='none' with outer_radius_nm != core_radius_nm", bad_shell_none_outer_mismatch),
    ("HIGH 3: shell='AlGaN' with outer_radius_nm < core_radius_nm + 3", bad_shell_algan_thin),
    ("HIGH 3: nitride.surface.shell contradicts nitride.nanowire.shell", bad_shell_contradiction),
    ("LOW 14: drive.eta_load != 1.0", bad_eta_load),
    ("LOW 13: drive.diode.R_s_ohm contradicts nitride.wire_thermal.R_s_ohm", bad_R_s_contradiction),
):
    check("rejects: " + name, raises(fn))

check("missing dot.gamma300 rejected at load-time card guard (legacy safety, unchanged)",
      raises(lambda: evaluate(DeviceDesign(platform="ingan_gan_nanowire"))))


def _card_missing_gamma300_rejected():
    """MEDIUM 11 (fix-2 round 2): device.py:566-575's own presence check
    fires inside DeviceDesign.load() itself (against the raw YAML mapping,
    before DotBlock's class defaults fill gamma300 in) -- not evaluate()'s
    finiteness check above, which never actually reaches load().  Build an
    in-memory YAML doc for a nanowire card with dot.gamma300 omitted and
    feed it through DeviceDesign.load() by monkeypatching Path.read_text
    for one call, mirroring verify_nitride_device.py's own
    _card_missing_gamma0_rejected -- no file is written."""
    from dataclasses import asdict
    import yaml
    d = card()
    doc = asdict(d)
    doc["dot"].pop("gamma300")
    text = yaml.safe_dump({"meta": {"name": "scratch"}, "design": doc}, sort_keys=False)
    orig_read_text = Path.read_text
    Path.read_text = lambda self, *a, **kw: text
    try:
        DeviceDesign.load("<in-memory nanowire card, gamma300 omitted>")
        return False
    except ValueError:
        return True
    finally:
        Path.read_text = orig_read_text


check("MEDIUM 11: missing dot.gamma300 rejected by DeviceDesign.load() itself on a dict card (not only evaluate())",
      _card_missing_gamma300_rejected())


# --------------------------------------------------------------- 2. planar bit-identity fixture

# MEDIUM 12 (fix-2 round 2): the previous fixture covered 11 scalars of one
# card at one T_hs.  This covers EVERY scalar of EVERY planar
# (platform=ingan_gan_planar) card under cards/, at T_grid=[230, 300] each
# (14 rows total, ~88 scalars per row) -- captured once, after this module's
# own implementation was finished (dev is already imported above, at the top
# of this file, before this fixture runs), as the golden reference going
# forward. The reviewer independently confirmed a 0-line diff in
# fsim_core/device.py's planar-relevant code between 605c00c and this
# commit (".workers/review/nitride-nanowire-device-opus-findings.md",
# "Opus re-review of f877b46" section, "legacy PASS byte-identical"); this
# fixture is the CLAUDE.md-documented alternative to re-deriving that
# historical comparison at runtime ("record the values the reviewer
# verified"), and additionally protects every planar card/T_hs going
# forward against any future perturbation this module's presence, import,
# or dispatch might introduce.
nan = float("nan")
PLANAR_CARDS = [
    "nitride-cavity-pulse-design.yaml",
    "nitride-cavity-set-design.yaml",
    "nitride-deshpande2014-comparison-design.yaml",
    "nitride-nonpolar-pulse-design.yaml",
    "nitride-nonpolar-set-design.yaml",
    "nitride-qw-fluctuation-pulse-design.yaml",
    "nitride-qw-fluctuation-set-design.yaml",
]
PLANAR_GOLDEN = {
('nitride-cavity-pulse-design.yaml', 230.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'rectangular', 'T_hs': 230.0, 'T_j': 230.00000068065262, 'g2_op': 0.9991379030238006, 'rho_pulsed': 0.90909175650759, 'collected_flux_pulsed_s': 3.880173957537539e-05, 'collected_flux_x_s': 3.880134171272575e-05, 'collected_flux_xx_s': 3.9786264964846e-10, 'background_flux_s': 3.880134171272575e-06, 'total_detected_flux_s': 4.268187374664797e-05, 'mean_counts': 4.8502174469219245e-12, 'mean_counts_x': 4.8501677140907184e-12, 'mean_counts_xx': 4.9732831206057494e-17, 'E_X_eV': 2.208920822280462, 'lambda_nm': 561.2885584192206, 'field_kVcm': -3901.8335667870415, 'overlap_sq': 0.003301931629964929, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 16.79288759736658, 'k_X_ns': 11709.673382194089, 'k_XX_ns': 23419.346764388178, 'gamma_X0_ns': 0.003301931629964929, 'gamma_XX0_ns': 0.006603863259929858, 'gamma_X_ns': 0.0053942794486723095, 'gamma_XX_ns': 0.008741852349353886, 'S_X': 4.6066843947164527e-07, 'S_XX': 3.73274676454032e-07, 'Q': 2000.0, 'kappa_meV': 1.0933554693422192, 'detuning_meV': 22.209883596023605, 'Fp_add': 75.99088773175333, 'F_eff_X': 1.6336738773508779, 'F_eff_XX': 1.3237482372472287, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 309135540.64414084, 'mu_resolved': 0.030913554064414085, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.452810661041165, 'V_terminal': 2.4528126610411647, 'power_W': 3.884081592483214e-10, 'counting_converged': True, 'blocked_load_probability': nan, 'one_pair_valid': False, 'set_feasible': False, 'set_priced_F_p': nan, 'set_E_C_meV': nan, 'set_EC_over_kT': nan, 'set_radius_nm': nan, 'set_radius_max_nm': nan, 'set_C_sigma_F': nan, 'set_R_T_over_RQ': nan, 'set_f_max_Hz': nan, 'pair_supply_possible': True, 'ideal_load_F_p': nan, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 302.8530303065728, 'tau_rad_cavity_ns': 185.38157125807217, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 136.70048584617953, 'applied_field_kVcm': 136.70048584617953, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 3.4396357544776124, 'reservoir_offset_meV': 1230.7149321971503, 'optical_reservoir_energy_eV': 3.4396357544776124, 'optical_reservoir_kind': 'gan_barrier', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'isolated_dot', 'effective_height_nm': 3.0, 'effective_radius_nm': 10.0, 'cavity_reference_V_j_V': 2.386891298935645, 'cavity_reference_transition_eV': 2.1839109387116644, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-cavity-pulse-design.yaml', 300.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'rectangular', 'T_hs': 300.0, 'T_j': 300.0000006829321, 'g2_op': 0.9996791936846767, 'rho_pulsed': 0.9090917842068371, 'collected_flux_pulsed_s': 0.004247138693754914, 'collected_flux_x_s': 0.0042470937212596965, 'collected_flux_xx_s': 4.497249521672975e-08, 'background_flux_s': 0.00042470937212596974, 'total_detected_flux_s': 0.004671848065880883, 'mean_counts': 5.308923367193642e-10, 'mean_counts_x': 5.308867151574621e-10, 'mean_counts_xx': 5.6215619020912184e-15, 'E_X_eV': 2.1839109384684114, 'lambda_nm': 567.7163670737909, 'field_kVcm': -3902.748396421988, 'overlap_sq': 0.0032189988716036337, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 16.77901297723189, 'k_X_ns': 18622.169199976695, 'k_XX_ns': 37244.33839995339, 'gamma_X0_ns': 0.0032189988716036337, 'gamma_XX0_ns': 0.006437997743207267, 'gamma_X_ns': 0.010619767973282507, 'gamma_XX_ns': 0.01776219933327007, 'S_X': 5.702752350196195e-07, 'S_XX': 4.769098237579364e-07, 'Q': 2000.0, 'kappa_meV': 1.0919554693421736, 'detuning_meV': -2.1593571375433385e-07, 'Fp_add': 75.99088773175333, 'F_eff_X': 3.2990903062951293, 'F_eff_XX': 2.7589632742588286, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 309059491.7248231, 'mu_resolved': 0.030905949172482315, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.3868912982517547, 'V_terminal': 2.3868932982517546, 'power_W': 3.778979081712954e-10, 'counting_converged': True, 'blocked_load_probability': nan, 'one_pair_valid': False, 'set_feasible': False, 'set_priced_F_p': nan, 'set_E_C_meV': nan, 'set_EC_over_kT': nan, 'set_radius_nm': nan, 'set_radius_max_nm': nan, 'set_C_sigma_F': nan, 'set_R_T_over_RQ': nan, 'set_f_max_Hz': nan, 'pair_supply_possible': True, 'ideal_load_F_p': nan, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 310.6555919672697, 'tau_rad_cavity_ns': 94.16401587264679, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 135.7856562112328, 'applied_field_kVcm': 135.7856562112328, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 3.4126017696256388, 'reservoir_offset_meV': 1228.6908311572274, 'optical_reservoir_energy_eV': 3.4126017696256388, 'optical_reservoir_kind': 'gan_barrier', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'isolated_dot', 'effective_height_nm': 3.0, 'effective_radius_nm': 10.0, 'cavity_reference_V_j_V': 2.386891298935645, 'cavity_reference_transition_eV': 2.1839109387116644, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-cavity-set-design.yaml', 230.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'deterministic_pair', 'T_hs': 230.0, 'T_j': 230.00000068065262, 'g2_op': 0.17355371900826433, 'rho_pulsed': 0.9090909090909092, 'collected_flux_pulsed_s': 0.001255189394466358, 'collected_flux_x_s': 0.001255189394466358, 'collected_flux_xx_s': 0.0, 'background_flux_s': 0.00012551893944663583, 'total_detected_flux_s': 0.0013807083339129938, 'mean_counts': 1.5689867430829476e-10, 'mean_counts_x': 1.5689867430829476e-10, 'mean_counts_xx': 0.0, 'E_X_eV': 2.208920822280462, 'lambda_nm': 561.2885584192206, 'field_kVcm': -3901.8335667870415, 'overlap_sq': 0.003301931629964929, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 16.79288759736658, 'k_X_ns': 11709.673382194089, 'k_XX_ns': 23419.346764388178, 'gamma_X0_ns': 0.003301931629964929, 'gamma_XX0_ns': 0.006603863259929858, 'gamma_X_ns': 0.0053942794486723095, 'gamma_XX_ns': 0.008741852349353886, 'S_X': 4.6066843947164527e-07, 'S_XX': 3.73274676454032e-07, 'Q': 2000.0, 'kappa_meV': 1.0933554693422192, 'detuning_meV': 22.209883596023605, 'Fp_add': 75.99088773175333, 'F_eff_X': 1.6336738773508779, 'F_eff_XX': 1.3237482372472287, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 309135540.64414084, 'mu_resolved': 0.030913554064414085, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.452810661041165, 'V_terminal': 2.4528126610411647, 'power_W': 3.884081592483214e-10, 'counting_converged': True, 'blocked_load_probability': 0.0, 'one_pair_valid': True, 'set_feasible': False, 'set_priced_F_p': 1.0, 'set_E_C_meV': 30.31504311247509, 'set_EC_over_kT': 1.5295281135357035, 'set_radius_nm': 5.0, 'set_radius_max_nm': 0.7647640567678519, 'set_C_sigma_F': 5.285087763377385e-18, 'set_R_T_over_RQ': 38.74045864977526, 'set_f_max_Hz': 18921161667.91825, 'pair_supply_possible': True, 'ideal_load_F_p': 0.0, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 302.8530303065728, 'tau_rad_cavity_ns': 185.38157125807217, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 136.70048584617953, 'applied_field_kVcm': 136.70048584617953, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 3.4396357544776124, 'reservoir_offset_meV': 1230.7149321971503, 'optical_reservoir_energy_eV': 3.4396357544776124, 'optical_reservoir_kind': 'gan_barrier', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'isolated_dot', 'effective_height_nm': 3.0, 'effective_radius_nm': 10.0, 'cavity_reference_V_j_V': 2.386891298935645, 'cavity_reference_transition_eV': 2.1839109387116644, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-cavity-set-design.yaml', 300.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'deterministic_pair', 'T_hs': 300.0, 'T_j': 300.0000006829321, 'g2_op': 0.17355371900826455, 'rho_pulsed': 0.9090909090909091, 'collected_flux_pulsed_s': 0.13742222076335409, 'collected_flux_x_s': 0.13742222076335409, 'collected_flux_xx_s': 0.0, 'background_flux_s': 0.01374222207633541, 'total_detected_flux_s': 0.1511644428396895, 'mean_counts': 1.717777759541926e-08, 'mean_counts_x': 1.717777759541926e-08, 'mean_counts_xx': 0.0, 'E_X_eV': 2.1839109384684114, 'lambda_nm': 567.7163670737909, 'field_kVcm': -3902.748396421988, 'overlap_sq': 0.0032189988716036337, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 16.77901297723189, 'k_X_ns': 18622.169199976695, 'k_XX_ns': 37244.33839995339, 'gamma_X0_ns': 0.0032189988716036337, 'gamma_XX0_ns': 0.006437997743207267, 'gamma_X_ns': 0.010619767973282507, 'gamma_XX_ns': 0.01776219933327007, 'S_X': 5.702752350196195e-07, 'S_XX': 4.769098237579364e-07, 'Q': 2000.0, 'kappa_meV': 1.0919554693421736, 'detuning_meV': -2.1593571375433385e-07, 'Fp_add': 75.99088773175333, 'F_eff_X': 3.2990903062951293, 'F_eff_XX': 2.7589632742588286, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 309059491.7248231, 'mu_resolved': 0.030905949172482315, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.3868912982517547, 'V_terminal': 2.3868932982517546, 'power_W': 3.778979081712954e-10, 'counting_converged': True, 'blocked_load_probability': 0.0, 'one_pair_valid': True, 'set_feasible': False, 'set_priced_F_p': 1.0, 'set_E_C_meV': 30.31504311247509, 'set_EC_over_kT': 1.1726382211781896, 'set_radius_nm': 5.0, 'set_radius_max_nm': 0.5863191105890947, 'set_C_sigma_F': 5.285087763377385e-18, 'set_R_T_over_RQ': 38.74045864977526, 'set_f_max_Hz': 18921161667.91825, 'pair_supply_possible': True, 'ideal_load_F_p': 0.0, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 310.6555919672697, 'tau_rad_cavity_ns': 94.16401587264679, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 135.7856562112328, 'applied_field_kVcm': 135.7856562112328, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 3.4126017696256388, 'reservoir_offset_meV': 1228.6908311572274, 'optical_reservoir_energy_eV': 3.4126017696256388, 'optical_reservoir_kind': 'gan_barrier', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'isolated_dot', 'effective_height_nm': 3.0, 'effective_radius_nm': 10.0, 'cavity_reference_V_j_V': 2.386891298935645, 'cavity_reference_transition_eV': 2.1839109387116644, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-deshpande2014-comparison-design.yaml', 230.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'rectangular', 'T_hs': 230.0, 'T_j': 230.000001368483, 'g2_op': 0.9907533560206871, 'rho_pulsed': 0.9090914420928513, 'collected_flux_pulsed_s': 0.0009098505843099325, 'collected_flux_x_s': 0.0009098447163926163, 'collected_flux_xx_s': 5.86791731634541e-09, 'background_flux_s': 9.098447163926162e-05, 'total_detected_flux_s': 0.0010008350559491941, 'mean_counts': 4.5492529215496624e-11, 'mean_counts_x': 4.5492235819630806e-11, 'mean_counts_xx': 2.933958658172705e-16, 'E_X_eV': 2.0277826554563214, 'lambda_nm': 611.4274528701857, 'field_kVcm': -6293.470230589916, 'overlap_sq': 0.04466462043248278, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 63.53024748309088, 'k_X_ns': 1107.7252740910933, 'k_XX_ns': 2215.4505481821866, 'gamma_X0_ns': 0.04466462043248278, 'gamma_XX0_ns': 0.08932924086496556, 'gamma_X_ns': 0.07473926506248701, 'gamma_XX_ns': 0.11884680015425103, 'S_X': 6.746638758024502e-05, 'S_XX': 5.364165095591524e-05, 'Q': 2000.0, 'kappa_meV': 1.0037158188254973, 'detuning_meV': 20.351017805326777, 'Fp_add': 75.99088773175333, 'F_eff_X': 1.6733437861733658, 'F_eff_XX': 1.33043557746008, 'eta_out': 0.1, 'gate_ns_used': 5.0, 'rep_rate_hz': 200000000.0, 'r_dot_s': 19245073.702180263, 'mu_resolved': 0.0019245073702180266, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 32.04353268, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 1.9727518717532737, 'V_terminal': 1.9727538717532738, 'power_W': 7.809122516323953e-10, 'counting_converged': True, 'blocked_load_probability': nan, 'one_pair_valid': False, 'set_feasible': False, 'set_priced_F_p': nan, 'set_E_C_meV': nan, 'set_EC_over_kT': nan, 'set_radius_nm': nan, 'set_radius_max_nm': nan, 'set_C_sigma_F': nan, 'set_R_T_over_RQ': nan, 'set_f_max_Hz': nan, 'pair_supply_possible': True, 'ideal_load_F_p': nan, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 22.38908537265303, 'tau_rad_cavity_ns': 13.379847917475951, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 175.98467092825504, 'applied_field_kVcm': 175.98467092825504, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 3.439635754235719, 'reservoir_offset_meV': 1411.853098779398, 'optical_reservoir_energy_eV': 3.439635754235719, 'optical_reservoir_kind': 'gan_barrier', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'isolated_dot', 'effective_height_nm': 2.0, 'effective_radius_nm': 12.5, 'cavity_reference_V_j_V': 1.909032479167014, 'cavity_reference_transition_eV': 2.004631637705734, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-deshpande2014-comparison-design.yaml', 300.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'rectangular', 'T_hs': 300.0, 'T_j': 300.0000013653411, 'g2_op': 0.9980609340138185, 'rho_pulsed': 0.9090910413777987, 'collected_flux_pulsed_s': 0.018999124213812685, 'collected_flux_x_s': 0.018999093802463018, 'collected_flux_xx_s': 3.041134966553999e-08, 'background_flux_s': 0.001899909380246302, 'total_detected_flux_s': 0.020899033594058987, 'mean_counts': 9.499562106906342e-10, 'mean_counts_x': 9.499546901231509e-10, 'mean_counts_xx': 1.520567483276999e-15, 'E_X_eV': 2.0046316372313706, 'lambda_nm': 618.4886843910963, 'field_kVcm': -6294.381633328483, 'overlap_sq': 0.04385504869483824, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 63.51348173291438, 'k_X_ns': 3054.398888765918, 'k_XX_ns': 6108.797777531836, 'gamma_X0_ns': 0.04385504869483824, 'gamma_XX0_ns': 0.08771009738967649, 'gamma_X_ns': 0.13663523865049135, 'gamma_XX_ns': 0.22951054272167265, 'S_X': 4.473191998479668e-05, 'S_XX': 3.756908127906067e-05, 'Q': 2000.0, 'kappa_meV': 1.00231581882556, 'detuning_meV': -4.197495684366004e-07, 'Fp_add': 75.99088773175333, 'F_eff_X': 3.115610236833995, 'F_eff_XX': 2.616694651495007, 'eta_out': 0.1, 'gate_ns_used': 5.0, 'rep_rate_hz': 200000000.0, 'r_dot_s': 7658837.350678194, 'mu_resolved': 0.0007658837350678195, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 32.04353268, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 1.9090324778473238, 'V_terminal': 1.9090344778473238, 'power_W': 7.555063307970113e-10, 'counting_converged': True, 'blocked_load_probability': nan, 'one_pair_valid': False, 'set_feasible': False, 'set_priced_F_p': nan, 'set_E_C_meV': nan, 'set_EC_over_kT': nan, 'set_radius_nm': nan, 'set_radius_max_nm': nan, 'set_C_sigma_F': nan, 'set_R_T_over_RQ': nan, 'set_f_max_Hz': nan, 'pair_supply_possible': True, 'ideal_load_F_p': nan, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 22.80239173734404, 'tau_rad_cavity_ns': 7.318756199913908, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 175.07326818968758, 'applied_field_kVcm': 175.07326818968758, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 3.412601769339992, 'reservoir_offset_meV': 1407.9701321086216, 'optical_reservoir_energy_eV': 3.412601769339992, 'optical_reservoir_kind': 'gan_barrier', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'isolated_dot', 'effective_height_nm': 2.0, 'effective_radius_nm': 12.5, 'cavity_reference_V_j_V': 1.909032479167014, 'cavity_reference_transition_eV': 2.004631637705734, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-nonpolar-pulse-design.yaml', 230.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'rectangular', 'T_hs': 230.0, 'T_j': 230.00000068065262, 'g2_op': 0.8762316964092836, 'rho_pulsed': 0.9090909090909217, 'collected_flux_pulsed_s': 5.298015525163914e-10, 'collected_flux_x_s': 5.298015525163105e-10, 'collected_flux_xx_s': 8.095685039830635e-23, 'background_flux_s': 5.298015525163105e-11, 'total_detected_flux_s': 5.827817077680225e-10, 'mean_counts': 6.622519406454893e-17, 'mean_counts_x': 6.622519406453881e-17, 'mean_counts_xx': 1.0119606299788292e-29, 'E_X_eV': 2.905843141084041, 'lambda_nm': 426.6720272923851, 'field_kVcm': 136.70048584617953, 'overlap_sq': 0.9681191226443326, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 152.08648008312636, 'k_X_ns': 72.38103748845266, 'k_XX_ns': 144.76207497690532, 'gamma_X0_ns': 0.9681191226443326, 'gamma_XX0_ns': 1.9362382452886653, 'gamma_X_ns': 1.9048200049803783, 'gamma_XX_ns': 2.863550396205122, 'S_X': 0.02564175832726662, 'S_XX': 0.019397380291990344, 'Q': 2000.0, 'kappa_meV': 1.4426394478580713, 'detuning_meV': 20.56424536789825, 'Fp_add': 75.99088773175333, 'F_eff_X': 1.967547133846018, 'F_eff_XX': 1.4789246122851003, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 0.036581831477622624, 'mu_resolved': 3.658183147762263e-12, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.452810661041165, 'V_terminal': 2.4528126610411647, 'power_W': 3.884081592483214e-10, 'counting_converged': True, 'blocked_load_probability': nan, 'one_pair_valid': False, 'set_feasible': False, 'set_priced_F_p': nan, 'set_E_C_meV': nan, 'set_EC_over_kT': nan, 'set_radius_nm': nan, 'set_radius_max_nm': nan, 'set_C_sigma_F': nan, 'set_R_T_over_RQ': nan, 'set_f_max_Hz': nan, 'pair_supply_possible': True, 'ideal_load_F_p': nan, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 1.0329307381808424, 'tau_rad_cavity_ns': 0.5249839866157333, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 136.70048584617953, 'applied_field_kVcm': 136.70048584617953, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 3.4396357544776124, 'reservoir_offset_meV': 533.7926133935715, 'optical_reservoir_energy_eV': 3.4396357544776124, 'optical_reservoir_kind': 'gan_barrier', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'isolated_dot', 'effective_height_nm': 3.0, 'effective_radius_nm': 10.0, 'cavity_reference_V_j_V': 2.386891298935645, 'cavity_reference_transition_eV': 2.8824788957433687, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-nonpolar-pulse-design.yaml', 300.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'rectangular', 'T_hs': 300.0, 'T_j': 300.0000006829321, 'g2_op': 0.9895138559015328, 'rho_pulsed': 0.9090909090910417, 'collected_flux_pulsed_s': 3.142264791100479e-07, 'collected_flux_x_s': 3.1422647910954366e-07, 'collected_flux_xx_s': 5.042146470683037e-19, 'background_flux_s': 3.142264791095437e-08, 'total_detected_flux_s': 3.4564912702100227e-07, 'mean_counts': 3.927830988875599e-14, 'mean_counts_x': 3.927830988869296e-14, 'mean_counts_xx': 6.302683088353796e-26, 'E_X_eV': 2.882478895514865, 'lambda_nm': 430.13046372314926, 'field_kVcm': 135.7856562112328, 'overlap_sq': 0.9679806914019865, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 152.21130004799122, 'k_X_ns': 563.0045391479651, 'k_XX_ns': 1126.0090782959303, 'gamma_X0_ns': 0.9679806914019865, 'gamma_XX0_ns': 1.935961382803973, 'gamma_X_ns': 3.8771648256202624, 'gamma_XX_ns': 6.407457086794771, 'S_X': 0.006839460152696908, 'S_XX': 0.005658215759476905, 'Q': 2000.0, 'kappa_meV': 1.4412394478580257, 'detuning_meV': -2.0118662291679357e-07, 'Fp_add': 75.99088773175333, 'F_eff_X': 4.005415459274011, 'F_eff_XX': 3.309702943306882, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 1.4606188449323623, 'mu_resolved': 1.4606188449323623e-10, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.3868912982517547, 'V_terminal': 2.3868932982517546, 'power_W': 3.778979081712954e-10, 'counting_converged': True, 'blocked_load_probability': nan, 'one_pair_valid': False, 'set_feasible': False, 'set_priced_F_p': nan, 'set_E_C_meV': nan, 'set_EC_over_kT': nan, 'set_radius_nm': nan, 'set_radius_max_nm': nan, 'set_C_sigma_F': nan, 'set_R_T_over_RQ': nan, 'set_f_max_Hz': nan, 'pair_supply_possible': True, 'ideal_load_F_p': nan, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 1.033078457951096, 'tau_rad_cavity_ns': 0.2579204250982602, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 135.7856562112328, 'applied_field_kVcm': 135.7856562112328, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 3.4126017696256388, 'reservoir_offset_meV': 530.122874110774, 'optical_reservoir_energy_eV': 3.4126017696256388, 'optical_reservoir_kind': 'gan_barrier', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'isolated_dot', 'effective_height_nm': 3.0, 'effective_radius_nm': 10.0, 'cavity_reference_V_j_V': 2.386891298935645, 'cavity_reference_transition_eV': 2.8824788957433687, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-nonpolar-set-design.yaml', 230.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'deterministic_pair', 'T_hs': 230.0, 'T_j': 230.00000068065262, 'g2_op': 0.17355371900826455, 'rho_pulsed': 0.9090909090909091, 'collected_flux_pulsed_s': 144.8264154955226, 'collected_flux_x_s': 144.8264154955226, 'collected_flux_xx_s': 0.0, 'background_flux_s': 14.48264154955226, 'total_detected_flux_s': 159.30905704507487, 'mean_counts': 1.8103301936940324e-05, 'mean_counts_x': 1.8103301936940324e-05, 'mean_counts_xx': 0.0, 'E_X_eV': 2.905843141084041, 'lambda_nm': 426.6720272923851, 'field_kVcm': 136.70048584617953, 'overlap_sq': 0.9681191226443326, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 152.08648008312636, 'k_X_ns': 72.38103748845266, 'k_XX_ns': 144.76207497690532, 'gamma_X0_ns': 0.9681191226443326, 'gamma_XX0_ns': 1.9362382452886653, 'gamma_X_ns': 1.9048200049803783, 'gamma_XX_ns': 2.863550396205122, 'S_X': 0.02564175832726662, 'S_XX': 0.019397380291990344, 'Q': 2000.0, 'kappa_meV': 1.4426394478580713, 'detuning_meV': 20.56424536789825, 'Fp_add': 75.99088773175333, 'F_eff_X': 1.967547133846018, 'F_eff_XX': 1.4789246122851003, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 0.036581831477622624, 'mu_resolved': 3.658183147762263e-12, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.452810661041165, 'V_terminal': 2.4528126610411647, 'power_W': 3.884081592483214e-10, 'counting_converged': True, 'blocked_load_probability': 0.0, 'one_pair_valid': True, 'set_feasible': False, 'set_priced_F_p': 1.0, 'set_E_C_meV': 30.31504311247509, 'set_EC_over_kT': 1.5295281135357035, 'set_radius_nm': 5.0, 'set_radius_max_nm': 0.7647640567678519, 'set_C_sigma_F': 5.285087763377385e-18, 'set_R_T_over_RQ': 38.74045864977526, 'set_f_max_Hz': 18921161667.91825, 'pair_supply_possible': True, 'ideal_load_F_p': 0.0, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 1.0329307381808424, 'tau_rad_cavity_ns': 0.5249839866157333, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 136.70048584617953, 'applied_field_kVcm': 136.70048584617953, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 3.4396357544776124, 'reservoir_offset_meV': 533.7926133935715, 'optical_reservoir_energy_eV': 3.4396357544776124, 'optical_reservoir_kind': 'gan_barrier', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'isolated_dot', 'effective_height_nm': 3.0, 'effective_radius_nm': 10.0, 'cavity_reference_V_j_V': 2.386891298935645, 'cavity_reference_transition_eV': 2.8824788957433687, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-nonpolar-set-design.yaml', 300.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'deterministic_pair', 'T_hs': 300.0, 'T_j': 300.0000006829321, 'g2_op': 0.17355371900826455, 'rho_pulsed': 0.9090909090909091, 'collected_flux_pulsed_s': 2151.3242842275577, 'collected_flux_x_s': 2151.3242842275577, 'collected_flux_xx_s': 0.0, 'background_flux_s': 215.13242842275577, 'total_detected_flux_s': 2366.4567126503134, 'mean_counts': 0.00026891553552844467, 'mean_counts_x': 0.00026891553552844467, 'mean_counts_xx': 0.0, 'E_X_eV': 2.882478895514865, 'lambda_nm': 430.13046372314926, 'field_kVcm': 135.7856562112328, 'overlap_sq': 0.9679806914019865, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 152.21130004799122, 'k_X_ns': 563.0045391479651, 'k_XX_ns': 1126.0090782959303, 'gamma_X0_ns': 0.9679806914019865, 'gamma_XX0_ns': 1.935961382803973, 'gamma_X_ns': 3.8771648256202624, 'gamma_XX_ns': 6.407457086794771, 'S_X': 0.006839460152696908, 'S_XX': 0.005658215759476905, 'Q': 2000.0, 'kappa_meV': 1.4412394478580257, 'detuning_meV': -2.0118662291679357e-07, 'Fp_add': 75.99088773175333, 'F_eff_X': 4.005415459274011, 'F_eff_XX': 3.309702943306882, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 1.4606188449323623, 'mu_resolved': 1.4606188449323623e-10, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.3868912982517547, 'V_terminal': 2.3868932982517546, 'power_W': 3.778979081712954e-10, 'counting_converged': True, 'blocked_load_probability': 0.0, 'one_pair_valid': True, 'set_feasible': False, 'set_priced_F_p': 1.0, 'set_E_C_meV': 30.31504311247509, 'set_EC_over_kT': 1.1726382211781896, 'set_radius_nm': 5.0, 'set_radius_max_nm': 0.5863191105890947, 'set_C_sigma_F': 5.285087763377385e-18, 'set_R_T_over_RQ': 38.74045864977526, 'set_f_max_Hz': 18921161667.91825, 'pair_supply_possible': True, 'ideal_load_F_p': 0.0, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 1.033078457951096, 'tau_rad_cavity_ns': 0.2579204250982602, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 135.7856562112328, 'applied_field_kVcm': 135.7856562112328, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 3.4126017696256388, 'reservoir_offset_meV': 530.122874110774, 'optical_reservoir_energy_eV': 3.4126017696256388, 'optical_reservoir_kind': 'gan_barrier', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'isolated_dot', 'effective_height_nm': 3.0, 'effective_radius_nm': 10.0, 'cavity_reference_V_j_V': 2.386891298935645, 'cavity_reference_transition_eV': 2.8824788957433687, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-qw-fluctuation-pulse-design.yaml', 230.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'rectangular', 'T_hs': 230.0, 'T_j': 230.00000066073818, 'g2_op': 0.9990557073844453, 'rho_pulsed': 0.9050750289732515, 'collected_flux_pulsed_s': 6.181405871099647e-06, 'collected_flux_x_s': 6.181335798776849e-06, 'collected_flux_xx_s': 7.007232279850572e-11, 'background_flux_s': 6.483106421401976e-07, 'total_detected_flux_s': 6.829716513239845e-06, 'mean_counts': 7.726757338874558e-13, 'mean_counts_x': 7.72666974847106e-13, 'mean_counts_xx': 8.759040349813213e-18, 'E_X_eV': 2.017501132340375, 'lambda_nm': 614.5433894065467, 'field_kVcm': -3898.585227106328, 'overlap_sq': 0.0005526434495131081, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 18.778379381758192, 'k_X_ns': 10593.477013431197, 'k_XX_ns': 21186.954026862393, 'gamma_X0_ns': 0.0005526434495131081, 'gamma_XX0_ns': 0.0011052868990262161, 'gamma_X_ns': 0.0008843519129597551, 'gamma_XX_ns': 0.0014403721718059723, 'S_X': 8.348079086895748e-08, 'S_XX': 6.798391463244896e-08, 'Q': 2000.0, 'kappa_meV': 0.9979064484480625, 'detuning_meV': 21.688235444249806, 'Fp_add': 75.99088773175333, 'F_eff_X': 1.6002214696272796, 'F_eff_XX': 1.3031658776331958, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 309135540.64414716, 'mu_resolved': 0.03091355406441472, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.3817857939512925, 'V_terminal': 2.3817877939512924, 'power_W': 3.770441805139323e-10, 'counting_converged': True, 'blocked_load_probability': nan, 'one_pair_valid': False, 'set_feasible': False, 'set_priced_F_p': nan, 'set_E_C_meV': nan, 'set_EC_over_kT': nan, 'set_radius_nm': nan, 'set_radius_max_nm': nan, 'set_C_sigma_F': nan, 'set_R_T_over_RQ': nan, 'set_f_max_Hz': nan, 'pair_supply_possible': True, 'ideal_load_F_p': nan, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 1809.484941658177, 'tau_rad_cavity_ns': 1130.7715688126834, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 139.94882552689288, 'applied_field_kVcm': 139.94882552689288, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 2.223792086291998, 'reservoir_offset_meV': 206.29095395162312, 'optical_reservoir_energy_eV': 2.223792086291998, 'optical_reservoir_kind': 'ingan_qw', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'qw_fluctuation', 'effective_height_nm': 3.5, 'effective_radius_nm': 20.0, 'cavity_reference_V_j_V': 2.294250192018214, 'cavity_reference_transition_eV': 1.9930128969225547, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-qw-fluctuation-pulse-design.yaml', 300.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'rectangular', 'T_hs': 300.0, 'T_j': 300.00000065614495, 'g2_op': 0.9996560535427652, 'rho_pulsed': 0.908184731646119, 'collected_flux_pulsed_s': 0.0006627641648061131, 'collected_flux_x_s': 0.0006627565515291735, 'collected_flux_xx_s': 7.613276939675203e-09, 'background_flux_s': 6.700384572278901e-05, 'total_detected_flux_s': 0.0007297680105289021, 'mean_counts': 8.284552060076414e-11, 'mean_counts_x': 8.284456894114668e-11, 'mean_counts_xx': 9.516596174594003e-16, 'E_X_eV': 1.9930128966771705, 'lambda_nm': 622.0943106123965, 'field_kVcm': -3897.5757525515864, 'overlap_sq': 0.0005395759492898227, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 18.793971050948045, 'k_X_ns': 17225.842121654183, 'k_XX_ns': 34451.684243308366, 'gamma_X0_ns': 0.0005395759492898227, 'gamma_XX0_ns': 0.0010791518985796454, 'gamma_X_ns': 0.0016746752569436426, 'gamma_XX_ns': 0.002813846917772047, 'S_X': 9.721876482477394e-08, 'S_XX': 8.167515608463179e-08, 'Q': 2000.0, 'kappa_meV': 0.9965064484481544, 'detuning_meV': -2.1913826309116757e-07, 'Fp_add': 75.99088773175333, 'F_eff_X': 3.103687736912313, 'F_eff_XX': 2.6074613976730863, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 309059491.724884, 'mu_resolved': 0.0309059491724884, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.29425019115853, 'V_terminal': 2.2942521911585296, 'power_W': 3.630753310363642e-10, 'counting_converged': True, 'blocked_load_probability': nan, 'one_pair_valid': False, 'set_feasible': False, 'set_priced_F_p': nan, 'set_E_C_meV': nan, 'set_EC_over_kT': nan, 'set_radius_nm': nan, 'set_radius_max_nm': nan, 'set_C_sigma_F': nan, 'set_R_T_over_RQ': nan, 'set_f_max_Hz': nan, 'pair_supply_possible': True, 'ideal_load_F_p': nan, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 1853.3072152607558, 'tau_rad_cavity_ns': 597.1306949533874, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 140.95830008163486, 'applied_field_kVcm': 140.95830008163486, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 2.199240884373482, 'reservoir_offset_meV': 206.22798769631157, 'optical_reservoir_energy_eV': 2.199240884373482, 'optical_reservoir_kind': 'ingan_qw', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'qw_fluctuation', 'effective_height_nm': 3.5, 'effective_radius_nm': 20.0, 'cavity_reference_V_j_V': 2.294250192018214, 'cavity_reference_transition_eV': 1.9930128969225547, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-qw-fluctuation-set-design.yaml', 230.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'deterministic_pair', 'T_hs': 230.0, 'T_j': 230.00000066073818, 'g2_op': 0.17378044090749156, 'rho_pulsed': 0.9089662034930168, 'collected_flux_pulsed_s': 0.00019996135028807528, 'collected_flux_x_s': 0.00019996135028807528, 'collected_flux_xx_s': 0.0, 'background_flux_s': 2.002631209107004e-05, 'total_detected_flux_s': 0.00021998766237914531, 'mean_counts': 2.499516878600941e-11, 'mean_counts_x': 2.499516878600941e-11, 'mean_counts_xx': 0.0, 'E_X_eV': 2.017501132340375, 'lambda_nm': 614.5433894065467, 'field_kVcm': -3898.585227106328, 'overlap_sq': 0.0005526434495131081, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 18.778379381758192, 'k_X_ns': 10593.477013431197, 'k_XX_ns': 21186.954026862393, 'gamma_X0_ns': 0.0005526434495131081, 'gamma_XX0_ns': 0.0011052868990262161, 'gamma_X_ns': 0.0008843519129597551, 'gamma_XX_ns': 0.0014403721718059723, 'S_X': 8.348079086895748e-08, 'S_XX': 6.798391463244896e-08, 'Q': 2000.0, 'kappa_meV': 0.9979064484480625, 'detuning_meV': 21.688235444249806, 'Fp_add': 75.99088773175333, 'F_eff_X': 1.6002214696272796, 'F_eff_XX': 1.3031658776331958, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 309135540.64414716, 'mu_resolved': 0.03091355406441472, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.3817857939512925, 'V_terminal': 2.3817877939512924, 'power_W': 3.770441805139323e-10, 'counting_converged': True, 'blocked_load_probability': 0.0, 'one_pair_valid': True, 'set_feasible': False, 'set_priced_F_p': 1.0, 'set_E_C_meV': 30.31504311247509, 'set_EC_over_kT': 1.5295281136681371, 'set_radius_nm': 5.0, 'set_radius_max_nm': 0.7647640568340686, 'set_C_sigma_F': 5.285087763377385e-18, 'set_R_T_over_RQ': 38.74045864977526, 'set_f_max_Hz': 18921161667.91825, 'pair_supply_possible': True, 'ideal_load_F_p': 0.0, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 1809.484941658177, 'tau_rad_cavity_ns': 1130.7715688126834, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 139.94882552689288, 'applied_field_kVcm': 139.94882552689288, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 2.223792086291998, 'reservoir_offset_meV': 206.29095395162312, 'optical_reservoir_energy_eV': 2.223792086291998, 'optical_reservoir_kind': 'ingan_qw', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'qw_fluctuation', 'effective_height_nm': 3.5, 'effective_radius_nm': 20.0, 'cavity_reference_V_j_V': 2.294250192018214, 'cavity_reference_transition_eV': 1.9930128969225547, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-qw-fluctuation-set-design.yaml', 300.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'deterministic_pair', 'T_hs': 300.0, 'T_j': 300.00000065614495, 'g2_op': 0.17360474096971346, 'rho_pulsed': 0.9090628465789846, 'collected_flux_pulsed_s': 0.021444687945288258, 'collected_flux_x_s': 0.021444687945288258, 'collected_flux_xx_s': 0.0, 'background_flux_s': 0.0021451969850986977, 'total_detected_flux_s': 0.023589884930386957, 'mean_counts': 2.680585993161032e-09, 'mean_counts_x': 2.680585993161032e-09, 'mean_counts_xx': 0.0, 'E_X_eV': 1.9930128966771705, 'lambda_nm': 622.0943106123965, 'field_kVcm': -3897.5757525515864, 'overlap_sq': 0.0005395759492898227, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 18.793971050948045, 'k_X_ns': 17225.842121654183, 'k_XX_ns': 34451.684243308366, 'gamma_X0_ns': 0.0005395759492898227, 'gamma_XX0_ns': 0.0010791518985796454, 'gamma_X_ns': 0.0016746752569436426, 'gamma_XX_ns': 0.002813846917772047, 'S_X': 9.721876482477394e-08, 'S_XX': 8.167515608463179e-08, 'Q': 2000.0, 'kappa_meV': 0.9965064484481544, 'detuning_meV': -2.1913826309116757e-07, 'Fp_add': 75.99088773175333, 'F_eff_X': 3.103687736912313, 'F_eff_XX': 2.6074613976730863, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 309059491.724884, 'mu_resolved': 0.0309059491724884, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.29425019115853, 'V_terminal': 2.2942521911585296, 'power_W': 3.630753310363642e-10, 'counting_converged': True, 'blocked_load_probability': 0.0, 'one_pair_valid': True, 'set_feasible': False, 'set_priced_F_p': 1.0, 'set_E_C_meV': 30.31504311247509, 'set_EC_over_kT': 1.172638221282895, 'set_radius_nm': 5.0, 'set_radius_max_nm': 0.5863191106414475, 'set_C_sigma_F': 5.285087763377385e-18, 'set_R_T_over_RQ': 38.74045864977526, 'set_f_max_Hz': 18921161667.91825, 'pair_supply_possible': True, 'ideal_load_F_p': 0.0, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 1853.3072152607558, 'tau_rad_cavity_ns': 597.1306949533874, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 140.95830008163486, 'applied_field_kVcm': 140.95830008163486, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 2.199240884373482, 'reservoir_offset_meV': 206.22798769631157, 'optical_reservoir_energy_eV': 2.199240884373482, 'optical_reservoir_kind': 'ingan_qw', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'qw_fluctuation', 'effective_height_nm': 3.5, 'effective_radius_nm': 20.0, 'cavity_reference_V_j_V': 2.294250192018214, 'cavity_reference_transition_eV': 1.9930128969225547, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
}

for (_card_name, _T), _golden in PLANAR_GOLDEN.items():
    _d = DeviceDesign.load(str(ROOT / "cards" / _card_name))
    _actual = evaluate(_d, [_T])["scalars"]
    check(f"planar bit-identity fixture matches golden scalars exactly (NaN-aware): {_card_name} @ {_T:.0f} K",
          nan_eq(_actual, _golden))

# Re-evaluate a second time (import already happened; this exercises the
# lru_cache'd levels()/rates() paths a second time) and confirm the SAME
# design object reproduces itself exactly -- catches any hidden state this
# module's import could have introduced into shared module-level caches.
base_planar = DeviceDesign.load(str(ROOT / "cards" / "nitride-cavity-pulse-design.yaml"))
planar_scalars = evaluate(base_planar, [300.0])["scalars"]
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
# LOW 16 (fix-2 round 2): the reservoir background channel is attenuated by
# eta_collection_X once, in addition to eta_out (rule 9) -- independently
# reconstructed here, not re-derived from the module's own formula.
_bg = hz["accepted_background_s"] * min(hz["gate_ns_used"], _pulse) * 1e-9 * hz["eta_out"] * hz["eta_collection_X"] + _d.drive.b_res * _sx
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
# MEDIUM 4 (fix-2 round 2): the backed-out external field must equal
# +(depletion_field_kVcm + external_field_kVcm) -- not merely be nonzero --
# so a reversed-sign feedback (e.g. subtracting the depletion field instead
# of adding it) fails this check even though the row would still
# "reproduce itself" internally.
check("MEDIUM 4: backed-out external field equals depletion_field_kVcm + card external_field_kVcm (sign not reversed)",
      close(_external_field, hz["depletion_field_kVcm"] + hz["external_field_kVcm"], rel=1e-6))

# HIGH 2: nitride.dot.external_field_kVcm at 500 kV/cm moves E_X_eV and is
# reported separately (not folded silently into field_kVcm alone).
_d_ext = card(family="horizontal_as_built", regime="rectangular", strain_bound="relaxed")
_d_ext.nitride["dot"]["external_field_kVcm"] = 500.0
_s_ext = evaluate(_d_ext)["scalars"]
check("HIGH 2: nitride.dot.external_field_kVcm=500 kV/cm changes E_X_eV",
      not close(_s_ext["E_X_eV"], hz["E_X_eV"], rel=1e-6, abs_=1e-6))
check("HIGH 2: nitride.dot.external_field_kVcm=500 kV/cm shifts field_kVcm by exactly +500 relative to the same T_j baseline",
      close(_s_ext["field_kVcm"] - hz["field_kVcm"], 500.0, rel=1e-6, abs_=1e-3))
check("HIGH 2: external_field_kVcm is reported separately on the row",
      close(_s_ext["external_field_kVcm"], 500.0) and close(hz["external_field_kVcm"], 0.0))

# HIGH 3: nitride.nanowire.shell propagates to surface (shell_multiplier_used,
# k_side_ns, k_X_ns all move for an AlGaN shell vs the 'none' default).
_s_none = evaluate(card(family="horizontal_as_built", regime="rectangular", strain_bound="relaxed"))["scalars"]
_s_algan = evaluate(card(family="horizontal_as_built", regime="rectangular", strain_bound="relaxed",
                          outer=15.5, extra_nanowire={"shell": "AlGaN"}))["scalars"]
check("HIGH 3: AlGaN shell changes shell_multiplier_used vs shell='none'",
      not close(_s_algan["shell_multiplier_used"], _s_none["shell_multiplier_used"], rel=1e-9))
check("HIGH 3: AlGaN shell changes k_X_ns vs shell='none'",
      not close(_s_algan["k_X_ns"], _s_none["k_X_ns"], rel=1e-6))
check("HIGH 3: k_side_ns itself is geometry-only (unaffected by the shell factor, only the occupied-dot/reservoir channels are)",
      close(_s_algan["k_side_ns"], _s_none["k_side_ns"], rel=1e-9))

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

# MEDIUM 8 (fix-2 round 2): eta_out is the genuine FilterBlock spectral
# transmission (spectral.epsilon, kappa=None), gated on d.filter.enabled --
# never the previous ad-hoc 1/(1+(2dx/w)**2) formula, and never the planar
# nitride path's cav.eta_out (a fixed cavity constant, no nanowire
# equivalent).
from fsim_core.spectral import epsilon as _epsilon
check("no ad-hoc (2*dx/w)**2 Lorentzian formula anywhere in the device module source",
      "2.0 * d.filter.dx" not in DEVICE_SRC and "(2.0 * dx" not in DEVICE_SRC)
_d_filt_on = card(family="horizontal_as_built", regime="rectangular", strain_bound="relaxed")
_d_filt_on.filter.enabled = True; _d_filt_on.filter.auto_w = False; _d_filt_on.filter.w = 1.0; _d_filt_on.filter.dx = 0.5
_s_filt_on = evaluate(_d_filt_on)["scalars"]
_expected_eta_out = float(_epsilon(0.0, _d_filt_on.dot.gamma300, _d_filt_on.dot.gamma300, w=1.0, kappa=None, dx=0.5).t_x)
check("MEDIUM 8: eta_out with an explicit fixed filter window matches an independent spectral.epsilon call",
      close(_s_filt_on["eta_out"], _expected_eta_out, rel=1e-9))
_d_filt_off = card(family="horizontal_as_built", regime="rectangular", strain_bound="relaxed")
_d_filt_off.filter.enabled = False
_s_filt_off = evaluate(_d_filt_off)["scalars"]
check("MEDIUM 8: eta_out is 1.0 (full transmission) when drive.filter.enabled is False",
      close(_s_filt_off["eta_out"], 1.0))
check("MEDIUM 8: eta_out with the filter enabled and narrower than the line differs from the disabled case",
      not close(_s_filt_on["eta_out"], _s_filt_off["eta_out"], rel=1e-6))


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

# HIGH 3 (fix-2 round 2): outer_radius_nm is now tied to the declared shell
# geometry (== core for 'none', >= core+3 for 'AlGaN'), so an "outer only"
# sweep independent of disc now goes through shell='AlGaN' with outer at
# its contract-declared core+3 -- still leaves disc (the Coulomb-charged
# island) untouched, which is the property this sweep tests.
_vert_base = card(family="vertical_photonic", regime="deterministic_pair")
_vert_outer = card(family="vertical_photonic", regime="deterministic_pair", outer=83.0, extra_nanowire={"shell": "AlGaN"})
_vert_disc = card(family="vertical_photonic", regime="deterministic_pair", disc=15.0)
_sb, _so, _sd = (evaluate(x)["scalars"] for x in (_vert_base, _vert_outer, _vert_disc))
check("disc-radius sweep moves the Coulomb screen (set_E_C_meV changes)",
      not close(_sb["set_E_C_meV"], _sd["set_E_C_meV"], rel=1e-6))
check("outer-radius-only sweep moves photonics (V_number) but NOT the Coulomb screen",
      not close(_sb["V_number"], _so["V_number"], rel=1e-6) and close(_sb["set_E_C_meV"], _so["set_E_C_meV"], rel=1e-9, abs_=1e-12))
# MEDIUM 5 (fix-2 round 2): the previous check had an unconditional `or
# True` (unfalsifiable) tacked onto a comparison between two DIFFERENT
# regimes (deterministic_pair vs rectangular), which is not even expected to
# agree. Replaced with two decisive checks: (a) a numeric one -- varying
# ONLY drive.set_params.ec_margin (a Coulomb-screen-only knob, never read by
# the optical pipeline) flips set_feasible while g2_op/collected_flux_pulsed_s
# stay bit-identical; (b) a structural one -- the source computes g2/flux
# BEFORE calling injector_feasibility() at all, so an RT injector failure
# cannot possibly feed back into them.
_vert_set_loose_margin = card(family="vertical_photonic", regime="deterministic_pair", extra_set_params={"ec_margin": 0.1})
_s_loose_margin = evaluate(_vert_set_loose_margin)["scalars"]
check("MEDIUM 5: varying drive.set_params.ec_margin flips set_feasible without changing g2_op/collected_flux_pulsed_s",
      _s_loose_margin["set_feasible"] != _sb["set_feasible"]
      and close(_s_loose_margin["g2_op"], _sb["g2_op"])
      and close(_s_loose_margin["collected_flux_pulsed_s"], _sb["collected_flux_pulsed_s"]))
check("MEDIUM 5: structural -- g2_op is computed before injector_feasibility() runs (RT injector failure cannot overwrite it)",
      DEVICE_SRC.index('g2 = 1.0 - rho * rho') < DEVICE_SRC.index('rti = injector_feasibility('))
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

# MEDIUM 6 (fix-2 round 2): the previous boundary loop computed `got` from
# the same four literals it compared against, calling nothing -- replaced
# with cards tuned (via a monkeypatched pulse_counting.pulse_g2) to sit at
# g2=0.5 and flux=1000 +/- eps, evaluated through the real module, asserting
# its own optical_pass. b_res=0 isolates the background to its fixed,
# unmocked component (inj["accepted_background_s"]-derived), so the SAME
# probe row's own reported eta_out/rep_rate_hz/background_flux_s can be
# inverted to solve for the mean_counts_x/g2_dot that hit an exact target.
_orig_pulse_g2 = pulse_counting.pulse_g2
_mock_counts = {}


def _mock_pulse_g2(*args, **kwargs):
    return dict(_mock_counts)


def _boundary_row(mx, g2dot):
    _mock_counts.clear()
    _mock_counts.update(mean_counts_x=mx, mean_counts_xx=0.0, mean_counts=mx,
                         g2=g2dot, converged=True, one_pair_valid=True,
                         blocked_load_probability=0.0)
    dev.pulse_counting.pulse_g2 = _mock_pulse_g2
    try:
        d = card(family="horizontal_as_built", regime="rectangular", strain_bound="relaxed",
                  extra_drive_kw={"b_res": 0.0})
        return evaluate(d)["scalars"]
    finally:
        dev.pulse_counting.pulse_g2 = _orig_pulse_g2


_probe = _boundary_row(1.0, 0.0)
_eta_out_probe = _probe["eta_out"]; _rep_probe = _probe["rep_rate_hz"]
_bg_fixed = _probe["background_flux_s"] / _rep_probe  # constant regardless of mx/g2dot (b_res=0)


def _solve_boundary_row(flux_target, g2_target):
    signal = flux_target / _rep_probe
    mx = signal / _eta_out_probe
    rho = signal / (signal + _bg_fixed)
    g2dot = 1.0 - (1.0 - g2_target) / (rho * rho)
    return _boundary_row(mx, g2dot)


_eps = 1e-6
for _flux_t, _g2_t, _expect in (
    (1000.0 + _eps, 0.5 - _eps, True),   # both strictly inside -> pass
    (1000.0 - _eps, 0.5 - _eps, False),  # flux just below 1000 -> fail
    (1000.0 + _eps, 0.5 + _eps, False),  # g2 just above 0.5 -> fail
    (1000.0, 0.5, False),                # flux==1000 (>=, ok) but g2==0.5 (< strict, fails) -> fail
):
    _r = _solve_boundary_row(_flux_t, _g2_t)
    check(f"MEDIUM 6: module optical_pass boundary (monkeypatched counts) flux~={_flux_t:.6f}, g2~={_g2_t:.6f}",
          bool(_r["optical_pass"]) == _expect
          and close(_r["g2_op"], _g2_t, abs_=1e-6)
          and close(_r["collected_flux_pulsed_s"], _flux_t, abs_=1e-4))

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

# HIGH 1: invalid and valid rows carry EXACTLY the same key set (this
# module's own dynamically-derived _row_keys(), not a hardcoded list) --
# for both families and both regimes, not merely "close enough".
check("HIGH 1: invalid row key set equals the module's own full row-key set",
      set(fresh) == set(dev._row_keys()))
for _key, _r in rows.items():
    check(f"HIGH 1: valid row key set equals the module's own full row-key set: {_key}",
          set(_r["scalars"]) == set(dev._row_keys()))

# LOW 20 (fix-2 round 2): not_applicable is reserved for pulse (rectangular)
# rows; an invalid deterministic_pair (SET) row must report set_feasible
# False, not the pulse-regime sentinel.
d_bad_det = card(regime="deterministic_pair", extra_diode={"tau_pulse_ns": 0.0})
fresh_det = evaluate(d_bad_det, T_grid=None)["scalars"]
check("LOW 20: invalid deterministic_pair row reports set_feasible False (not not_applicable)",
      fresh_det["valid"] is False and fresh_det["cycle_loading"] == "deterministic_pair"
      and fresh_det["set_feasible"] is False and fresh_det["rti_feasible"] is False)
check("LOW 20: invalid rectangular (pulse) row still reports the not_applicable sentinel",
      fresh["cycle_loading"] == "rectangular" and fresh["set_feasible"] == "not_applicable"
      and fresh["rti_feasible"] == "not_applicable")

# MEDIUM 7 (fix-2 round 2): thermal_converged must be a real False on a
# thermal-solve failure (never left at the row-level NaN default).
d_thermal_fail = card(family="horizontal_as_built", regime="rectangular", I_uA=1e6, extra_thermal={"Rth_K_W": 1e15})
s_thermal_fail = evaluate(d_thermal_fail)["scalars"]
check("MEDIUM 7: thermal-solve failure reports thermal_converged False (not NaN)",
      s_thermal_fail["valid"] is False and s_thermal_fail["thermal_converged"] is False)
check("MEDIUM 7: a converged operating point reports thermal_converged True",
      hz["thermal_converged"] is True)

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
# MEDIUM 9 (fix-2 round 2): the real default horizontal SET card's unrelaxed
# partner emits ~82 photons/s (well below the 1000/s optical floor) --
# exactly the false-alarm case the finding names -- so the guarded
# _bound_reversal now reports "not_comparable" for it instead of a
# (spurious) True/False; "not_comparable" is a legitimate outcome here, not
# a test failure.
check("evaluate_strain_pair reports both bounds with a matching bound_reversal flag (True/False/not_comparable)",
      pair["relaxed"]["scalars"]["bound_reversal"] == pair["unrelaxed"]["scalars"]["bound_reversal"]
      and pair["relaxed"]["scalars"]["bound_reversal"] in (True, False, "not_comparable"))
check("MEDIUM 9: the real default horizontal SET pair (unrelaxed flux ~82/s, below the 1000/s floor) reports not_comparable",
      pair["relaxed"]["scalars"]["bound_reversal"] == "not_comparable")
check("evaluate_strain_pair bound_role labels: relaxed=headline_upper, unrelaxed=conservative_lower",
      pair["relaxed"]["scalars"]["bound_role"] == "headline_upper"
      and pair["unrelaxed"]["scalars"]["bound_role"] == "conservative_lower")
check("single-row evaluate_nanowire reports bound_reversal as the not_computed sentinel",
      rows[("horizontal_as_built", "deterministic_pair", "relaxed")]["scalars"]["bound_reversal"] == "not_computed")

# Manually construct synthetic pairs to exercise the reversal ARITHMETIC
# independent of whether real physics happens to reverse at this operating
# point (a direct, deterministic test of _bound_reversal's own logic).
# MEDIUM 9: both partners must be valid AND clear the 1000/s optical floor
# to be compared at all -- flux bumped to 5000/20000 (both well above 1000)
# so these constructed cases exercise the ORDERING logic, not the floor.
_relaxed_scalars = {"valid": True, "mu": 0.5, "collected_flux_pulsed_s": 5000.0, "g2_op": 0.9}
_unrelaxed_scalars = {"valid": True, "mu": 1.0, "collected_flux_pulsed_s": 20000.0, "g2_op": 0.1}
check("constructed bound_reversal pair: relaxed worse on all three metrics -> True",
      dev._bound_reversal(_relaxed_scalars, _unrelaxed_scalars) is True)
_relaxed_better = {"valid": True, "mu": 2.0, "collected_flux_pulsed_s": 50000.0, "g2_op": 0.01}
check("constructed bound_reversal pair: relaxed better on all three metrics -> False",
      dev._bound_reversal(_relaxed_better, _unrelaxed_scalars) is False)
check("constructed bound_reversal pair: an all-NaN (invalid) partner reports not_comparable",
      dev._bound_reversal({"valid": False, "mu": float("nan"), "collected_flux_pulsed_s": float("nan"), "g2_op": float("nan")}, _unrelaxed_scalars) == "not_comparable")
check("MEDIUM 9: a valid partner below the 1000/s optical floor reports not_comparable (the false-alarm case)",
      dev._bound_reversal({"valid": True, "mu": 0.1, "collected_flux_pulsed_s": 200.0, "g2_op": 0.5}, _unrelaxed_scalars) == "not_comparable")

# vertical headline eligibility (M12, non-contract bonus column): a wide
# vertical core pushed above the LP11 cutoff must not be headline-eligible.
s_wide = evaluate(card(family="vertical_photonic", regime="rectangular", core=120.0, outer=120.0, disc=12.5))["scalars"]
check("headline_eligible present as a bool for vertical rows", isinstance(s_wide.get("headline_eligible"), bool))
if s_wide["V_number"] > 2.404826:
    check("vertical row above the LP11 cutoff (single_mode False) is never headline_eligible",
          s_wide["single_mode"] is False and s_wide["headline_eligible"] is False)
else:
    check("vertical headline card is below cutoff (informational, not a failure)", True)


# --------------------------------------------------------------- injector import wiring

# MEDIUM 10 (fix-2 round 2): the import-retry scaffolding (_injector_symbols,
# six attempts with sleep(0.5)) was removed -- this module now imports
# NitrideNanowireInjectorParams/injector_feasibility as a plain top-level
# import, same as its other five dependency modules.
check("MEDIUM 10: injector_feasibility is imported as a plain top-level symbol (no retry scaffolding)",
      callable(dev.injector_feasibility) and "_injector_symbols" not in DEVICE_SRC and "importlib" not in DEVICE_SRC)


print(f"{sum(checks)}/{len(checks)} nitride nanowire device checks passed")
assert all(checks), "one or more nitride nanowire device checks failed (see FAIL lines above)"
