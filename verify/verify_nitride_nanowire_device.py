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


def _catch(fn):
    """Call fn() with no arguments and return the exception it raised (any
    type), or None if it did not raise -- for asserting on an error
    MESSAGE (MEDIUM 3's named-error checks), not just "raised/did not"."""
    try:
        fn()
    except Exception as exc:
        return exc
    return None


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

# --------------------------------------------------------------- MEDIUM 3: contract-parsing robustness


class _FakeContractPath:
    """A duck-typed stand-in for _CONTRACT_PATH's Path interface (no real
    temp file -- in-memory fixture per the spec's sandbox rule)."""
    def __init__(self, text=None, missing=False):
        self._text = text; self._missing = missing

    def read_text(self, encoding="utf-8"):
        if self._missing:
            raise FileNotFoundError(str(self))
        return self._text

    def __str__(self):
        return "<fake contract path>"


import warnings as _warnings

_orig_contract_path = dev._CONTRACT_PATH
try:
    dev._CONTRACT_PATH = _FakeContractPath(missing=True)
    with _warnings.catch_warnings(record=True) as _w:
        _warnings.simplefilter("always")
        _fallback_cols = dev._contract_row_columns()
    check("MEDIUM 3: a missing contract doc falls back to the frozen column list (not a bare FileNotFoundError)",
          _fallback_cols == list(dev._FROZEN_CONTRACT_COLUMNS))
    check("MEDIUM 3: a missing contract doc warns once (not silent)",
          any(issubclass(x.category, RuntimeWarning) for x in _w))

    dev._CONTRACT_PATH = _FakeContractPath(text="# Contract\n\nNo Row columns heading here.\n\n## Sweep grid\n")
    check("MEDIUM 3: a doc missing the '## Row columns' heading fails loudly, naming that heading",
          "Row columns" in str(_err) if (_err := _catch(dev._contract_row_columns)) else False)

    dev._CONTRACT_PATH = _FakeContractPath(text="# Contract\n\n## Row columns\n\n`foo_col`, `bar_col`.\n\n(No closing heading here.)\n")
    check("MEDIUM 3: a doc missing the closing '## Sweep grid' heading fails loudly, naming that heading",
          "Sweep grid" in str(_err) if (_err := _catch(dev._contract_row_columns)) else False)
finally:
    dev._CONTRACT_PATH = _orig_contract_path

# MEDIUM 3: the frozen fallback must not have silently drifted from the
# live contract doc (a real content change to the Row-columns section
# should be reflected here too, not just in CONTRACT_COLUMNS above).
check("MEDIUM 3: the frozen fallback column list matches a live parse of the real contract doc",
      list(dev._FROZEN_CONTRACT_COLUMNS) == dev._contract_row_columns())

# MEDIUM 3 (generalization, not just today's fixed list): a synthetic
# section demonstrates the STRUCTURAL depth-based rule -- a real column at
# paren depth 0 is kept, a citation/reference nested inside parens
# (including one spanning multiple lines, and one whose own closing colon
# sits inside the parenthetical) is excluded, with no per-word denylist
# entry required for the new citation.
_synthetic_section = (
    "Some intro sentence naming `evaluate_nanowire` before any column list.\n\n"
    "Identity: `real_column_one`, `real_column_two` (see `some_future_citation`,\n"
    "commit `abcdef1`: added in a later round) and `real_column_three`.\n\n"
    "`formula_column = real_column_one + real_column_two`.\n"
)
_synthetic_cols = dev._parse_row_columns_section(_synthetic_section)
check("MEDIUM 3: depth-based parser keeps depth-0 columns and the formula LHS",
      {"real_column_one", "real_column_two", "real_column_three", "formula_column"} <= set(_synthetic_cols))
check("MEDIUM 3: depth-based parser excludes a multi-line citation nested in parens (no denylist entry needed)",
      "some_future_citation" not in _synthetic_cols and "abcdef1" not in _synthetic_cols)
check("MEDIUM 3: depth-based parser excludes the opening sentence's function-name citation",
      "evaluate_nanowire" not in _synthetic_cols)
# (The "broken row-key fixture names itself" check lives further down,
# after card() is defined -- see MEDIUM 3 section 2.)


# --------------------------------------------------------------- card builder

def card(family="horizontal_as_built", regime="rectangular", T_hs=300.0,
         strain_bound="relaxed", core=None, outer=None, disc=None,
         R_s_ohm=None, x_in=0.40, occupied_dot_access=0.05, I_uA=0.02,
         rep_rate_hz=80e6, tau_pulse_ns=0.1, eta_total=0.01,
         shell_multiplier=1.0,
         extra_nanowire=None, extra_dot=None, extra_drive_kw=None,
         extra_diode=None, extra_set_params=None, extra_photonics=None,
         extra_thermal=None, extra_surface=None, extra_nitride=None, aperture=None):
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
    # HIGH 1 (fix-3): shell_multiplier is set EXPLICITLY here (not merely
    # left to NitrideNanowireSurfaceParams's own class default) so that a
    # test which asks for shell='AlGaN' (via extra_nanowire) exercises the
    # real card shape the fix actually guards -- a card that also states
    # nitride.surface.shell_multiplier, not one relying on an unstated
    # class default to accidentally look right.
    surface = {"occupied_dot_access": occupied_dot_access, "shell_multiplier": shell_multiplier}
    if extra_surface: surface.update(extra_surface)
    nitride = {"tau_rad0_ns": 1.0, "tau_cap_ps": 10.0, "nanowire": nanowire, "dot": dot,
               "surface": surface, "photonics": photonics, "wire_thermal": thermal_block,
               "injector": {}}
    if extra_nitride: nitride.update(extra_nitride)
    kw = dict(platform="ingan_gan_nanowire",
        dot=DotBlock(linewidth="anchored", lineshape="lorentzian", gamma300=3.0),
        ret=RetentionBlock(mode="nitride_confinement"),
        drive=DriveBlock(**drive_kw), thermal=ThermalBlock(T_hs=T_hs),
        cavity=CavityBlock(enabled=False), emission=EmissionBlock(type="nanowire"),
        nitride=nitride)
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


# MEDIUM 4 (fix-3): top-level nitride allow-list -- QW/wetting-layer/cavity
# leaves have no meaning for a nanowire card and must be rejected by name,
# not silently ignored.
def bad_nitride_qw_leaf():
    return evaluate(card(extra_nitride={"qw": {"wl_thickness_nm": 3.0}}))


def bad_nitride_wl_leaf():
    return evaluate(card(extra_nitride={"wl": {"thickness_nm": 3.0}}))


def bad_nitride_cavity_leaf():
    return evaluate(card(extra_nitride={"cavity": {"Q": 5000.0}}))


def bad_nitride_unknown_leaf():
    return evaluate(card(extra_nitride={"bogus_leaf": 1.0}))


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
    ("MEDIUM 4: nitride.qw is unrecognized on a nanowire card", bad_nitride_qw_leaf),
    ("MEDIUM 4: nitride.wl is unrecognized on a nanowire card", bad_nitride_wl_leaf),
    ("MEDIUM 4: nitride.cavity is unrecognized on a nanowire card", bad_nitride_cavity_leaf),
    ("MEDIUM 4: unknown top-level nitride leaf is rejected", bad_nitride_unknown_leaf),
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
# commit (byte-identical with legacy code); this
# fixture is the CLAUDE.md-documented alternative to re-deriving that
# historical comparison at runtime ("record the values the reviewer
# verified"), and additionally protects every planar card/T_hs going
# forward against any future perturbation this module's presence, import,
# or dispatch might introduce.
# RE-PINNED 2026-09-23 (audit Phase B, spec audit-device-composition, audit
# H3): the planar background acceptance is now evaluated on the X-centred
# window (single spectral selection), which moves background_flux_s,
# total_detected_flux_s, rho_pulsed and g2_op (QW-fluctuation cards: SET g2_op
# 0.1736 -> 0.4542, pulse rho 0.908 -> 0.107; other cards background-only
# changes); new keys background_acceptance / background_rate_window_s /
# background_window_ns / eta_background were added. Captured from the
# post-fix evaluator; everything else in this block is unchanged.
# RE-PINNED 2026-09-23 (audit Phase B, spec audit-nitride-coulomb,
# nitride_levels:339): the exciton Coulomb term is now the exact frozen-
# orbital sqrt(pi)/L*erfcx(z_sep/L) instead of the quadrature form, so
# E_X_eV rises on every planar card (+5.9 meV cavity, +3.8/+3.8 meV
# Deshpande, +1.1 meV nonpolar, +2.8 meV QW-fluctuation) and every E_X-
# dependent scalar (lambda, detuning, kappa, Purcell/F_eff, gamma, S, flux,
# counts, background, reservoir offset, rho, g2_op; QW-fluctuation SET
# g2_op @300 K 0.4542 -> 0.4769) moves with it. All 14 rows re-captured in
# full from the post-fix evaluator; the comparison is still exact
# full-dict NaN-aware equality (no tolerance, no key dropped).
# RE-PINNED AGAIN 2026-09-23 (spec audit-nitride-strain-mass): Rinke 2008
# electron-mass axes corrected (me_z 0.186 along c) and the anisotropic Yan
# 2014 A-gap strain shift replaces the volume-only a_V form (coherent x=0.25
# gap opening +114 meV instead of +250 meV). E_X_eV falls on every planar
# card (-75 meV cavity, -134 meV Deshpande 2014 comparison, -134 meV
# nonpolar, -74 meV QW-fluctuation) and every E_X-dependent scalar moves;
# the QW-fluctuation cards' reservoir_kind becomes 'mixed' (was
# 'gan_barrier'; SET g2_op @300 K 0.4769 -> 0.2289). All 14 rows re-captured
# in full from the post-fix evaluator; comparison still exact NaN-aware.
# RE-PINNED 2026-09-23 (spec audit-c0-d1-validity, audit C0/D2): the
# planar flat-background acceptance is now the closed-form arctan average
# (fsim_core/device.py _nitride_flat_background_acceptance) instead of a
# scipy quad; background_acceptance moves by <= 4.4e-15 rel on 12 rows and
# the downstream background_flux_s / total_detected_flux_s / rho_pulsed by
# <= 6.9e-16 rel on the QW-fluctuation pulse card (floating-point
# re-association only). Those values were re-captured; comparison still exact.
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
('nitride-cavity-pulse-design.yaml', 230.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'rectangular', 'T_hs': 230.0, 'T_j': 230.00000068034873, 'g2_op': 0.9971190400746142, 'rho_pulsed': 0.9090937125128973, 'collected_flux_pulsed_s': 0.0003067272546533103, 'collected_flux_x_s': 0.0003067168500656418, 'collected_flux_xx_s': 1.0404587668432894e-08, 'background_flux_s': 3.0671685006564186e-05, 'total_detected_flux_s': 0.00033739893965987447, 'mean_counts': 3.834090683166378e-11, 'mean_counts_x': 3.833960625820523e-11, 'mean_counts_xx': 1.3005734585541118e-15, 'E_X_eV': 2.1397504749566387, 'lambda_nm': 579.4329752515301, 'field_kVcm': -3901.840077787776, 'overlap_sq': 0.007979076053709917, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 40.66796905483216, 'k_X_ns': 3510.6869944332298, 'k_XX_ns': 7021.3739888664595, 'gamma_X0_ns': 0.007979076053709917, 'gamma_XX0_ns': 0.015958152107419833, 'gamma_X_ns': 0.013037833686522969, 'gamma_XX_ns': 0.021080516525137875, 'S_X': 3.713741865352201e-06, 'S_XX': 3.0023259362032105e-06, 'Q': 2000.0, 'kappa_meV': 1.0589856741063381, 'detuning_meV': 21.77912674396243, 'Fp_add': 75.99088773175333, 'F_eff_X': 1.6340029345203388, 'F_eff_XX': 1.3209873162780776, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 309141449.63791746, 'mu_resolved': 0.03091414496379175, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.451726857949691, 'V_terminal': 2.451728857949691, 'power_W': 3.882347507536855e-10, 'counting_converged': True, 'blocked_load_probability': nan, 'one_pair_valid': False, 'set_feasible': False, 'set_priced_F_p': nan, 'set_E_C_meV': nan, 'set_EC_over_kT': nan, 'set_radius_nm': nan, 'set_radius_max_nm': nan, 'set_C_sigma_F': nan, 'set_R_T_over_RQ': nan, 'set_f_max_Hz': nan, 'pair_supply_possible': True, 'ideal_load_F_p': nan, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 125.32779400379877, 'tau_rad_cavity_ns': 76.69985858415163, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 136.69397484544533, 'applied_field_kVcm': 136.69397484544533, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 3.4396357544777194, 'reservoir_offset_meV': 1299.8852795210807, 'background_acceptance': 0.0006860483088728302, 'background_rate_window_s': 1.9304259684586937e-45, 'background_window_ns': 0.1, 'eta_background': 0.1, 'optical_reservoir_energy_eV': 3.4396357544777194, 'optical_reservoir_kind': 'gan_barrier', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'isolated_dot', 'effective_height_nm': 3.0, 'effective_radius_nm': 10.0, 'cavity_reference_V_j_V': 2.385477642681538, 'cavity_reference_transition_eV': 2.1151713482398904, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-cavity-pulse-design.yaml', 300.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'rectangular', 'T_hs': 300.0, 'T_j': 300.00000068252336, 'g2_op': 0.9991946624029615, 'rho_pulsed': 0.9090931135917045, 'collected_flux_pulsed_s': 0.024776568686833426, 'collected_flux_x_s': 0.02477590778685501, 'collected_flux_xx_s': 6.608999784111582e-07, 'background_flux_s': 0.002477590778685501, 'total_detected_flux_s': 0.027254159465518926, 'mean_counts': 3.097071085854178e-09, 'mean_counts_x': 3.0969884733568764e-09, 'mean_counts_xx': 8.261249730139476e-14, 'E_X_eV': 2.1151713479951915, 'lambda_nm': 586.166215410941, 'field_kVcm': -3902.756933475466, 'overlap_sq': 0.00785853599655373, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 40.64893353929639, 'k_X_ns': 7396.6190682708875, 'k_XX_ns': 14793.238136541775, 'gamma_X0_ns': 0.00785853599655373, 'gamma_XX0_ns': 0.01571707199310746, 'gamma_X_ns': 0.02537401736480599, 'gamma_XX_ns': 0.04250618033984151, 'S_X': 3.43047682807857e-06, 'S_XX': 2.873343740744442e-06, 'Q': 2000.0, 'kappa_meV': 1.0575856741062948, 'detuning_meV': -2.173980995223701e-07, 'Fp_add': 75.99088773175333, 'F_eff_X': 3.228847889216703, 'F_eff_XX': 2.7044592248786605, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 309125153.2350401, 'mu_resolved': 0.030912515323504017, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.3854776419948416, 'V_terminal': 2.3854796419948414, 'power_W': 3.776717231701891e-10, 'counting_converged': True, 'blocked_load_probability': nan, 'one_pair_valid': False, 'set_feasible': False, 'set_priced_F_p': nan, 'set_E_C_meV': nan, 'set_EC_over_kT': nan, 'set_radius_nm': nan, 'set_radius_max_nm': nan, 'set_C_sigma_F': nan, 'set_R_T_over_RQ': nan, 'set_f_max_Hz': nan, 'pair_supply_possible': True, 'ideal_load_F_p': nan, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 127.25016471751717, 'tau_rad_cavity_ns': 39.41039314440644, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 135.7771191577549, 'applied_field_kVcm': 135.7771191577549, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 3.4126017696258097, 'reservoir_offset_meV': 1297.4304216306182, 'background_acceptance': 0.0465515605136271, 'background_rate_window_s': 4.014083382097431e-34, 'background_window_ns': 0.1, 'eta_background': 0.1, 'optical_reservoir_energy_eV': 3.4126017696258097, 'optical_reservoir_kind': 'gan_barrier', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'isolated_dot', 'effective_height_nm': 3.0, 'effective_radius_nm': 10.0, 'cavity_reference_V_j_V': 2.385477642681538, 'cavity_reference_transition_eV': 2.1151713482398904, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-cavity-set-design.yaml', 230.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'deterministic_pair', 'T_hs': 230.0, 'T_j': 230.00000068034873, 'g2_op': 0.17355371900826455, 'rho_pulsed': 0.9090909090909091, 'collected_flux_pulsed_s': 0.009922441109602753, 'collected_flux_x_s': 0.009922441109602753, 'collected_flux_xx_s': 0.0, 'background_flux_s': 0.0009922441109602752, 'total_detected_flux_s': 0.010914685220563027, 'mean_counts': 1.240305138700344e-09, 'mean_counts_x': 1.240305138700344e-09, 'mean_counts_xx': 0.0, 'E_X_eV': 2.1397504749566387, 'lambda_nm': 579.4329752515301, 'field_kVcm': -3901.840077787776, 'overlap_sq': 0.007979076053709917, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 40.66796905483216, 'k_X_ns': 3510.6869944332298, 'k_XX_ns': 7021.3739888664595, 'gamma_X0_ns': 0.007979076053709917, 'gamma_XX0_ns': 0.015958152107419833, 'gamma_X_ns': 0.013037833686522969, 'gamma_XX_ns': 0.021080516525137875, 'S_X': 3.713741865352201e-06, 'S_XX': 3.0023259362032105e-06, 'Q': 2000.0, 'kappa_meV': 1.0589856741063381, 'detuning_meV': 21.77912674396243, 'Fp_add': 75.99088773175333, 'F_eff_X': 1.6340029345203388, 'F_eff_XX': 1.3209873162780776, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 309141449.63791746, 'mu_resolved': 0.03091414496379175, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.451726857949691, 'V_terminal': 2.451728857949691, 'power_W': 3.882347507536855e-10, 'counting_converged': True, 'blocked_load_probability': 0.0, 'one_pair_valid': True, 'set_feasible': False, 'set_priced_F_p': 1.0, 'set_E_C_meV': 30.31504311247509, 'set_EC_over_kT': 1.5295281135377246, 'set_radius_nm': 5.0, 'set_radius_max_nm': 0.7647640567688623, 'set_C_sigma_F': 5.285087763377385e-18, 'set_R_T_over_RQ': 38.74045864977526, 'set_f_max_Hz': 18921161667.91825, 'pair_supply_possible': True, 'ideal_load_F_p': 0.0, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 125.32779400379877, 'tau_rad_cavity_ns': 76.69985858415163, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 136.69397484544533, 'applied_field_kVcm': 136.69397484544533, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 3.4396357544777194, 'reservoir_offset_meV': 1299.8852795210807, 'background_acceptance': 0.0006860483088728302, 'background_rate_window_s': 1.9304259684586937e-45, 'background_window_ns': 0.1, 'eta_background': 0.1, 'optical_reservoir_energy_eV': 3.4396357544777194, 'optical_reservoir_kind': 'gan_barrier', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'isolated_dot', 'effective_height_nm': 3.0, 'effective_radius_nm': 10.0, 'cavity_reference_V_j_V': 2.385477642681538, 'cavity_reference_transition_eV': 2.1151713482398904, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-cavity-set-design.yaml', 300.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'deterministic_pair', 'T_hs': 300.0, 'T_j': 300.00000068252336, 'g2_op': 0.17355371900826455, 'rho_pulsed': 0.9090909090909091, 'collected_flux_pulsed_s': 0.8015181423879548, 'collected_flux_x_s': 0.8015181423879548, 'collected_flux_xx_s': 0.0, 'background_flux_s': 0.08015181423879549, 'total_detected_flux_s': 0.8816699566267503, 'mean_counts': 1.0018976779849435e-07, 'mean_counts_x': 1.0018976779849435e-07, 'mean_counts_xx': 0.0, 'E_X_eV': 2.1151713479951915, 'lambda_nm': 586.166215410941, 'field_kVcm': -3902.756933475466, 'overlap_sq': 0.00785853599655373, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 40.64893353929639, 'k_X_ns': 7396.6190682708875, 'k_XX_ns': 14793.238136541775, 'gamma_X0_ns': 0.00785853599655373, 'gamma_XX0_ns': 0.01571707199310746, 'gamma_X_ns': 0.02537401736480599, 'gamma_XX_ns': 0.04250618033984151, 'S_X': 3.43047682807857e-06, 'S_XX': 2.873343740744442e-06, 'Q': 2000.0, 'kappa_meV': 1.0575856741062948, 'detuning_meV': -2.173980995223701e-07, 'Fp_add': 75.99088773175333, 'F_eff_X': 3.228847889216703, 'F_eff_XX': 2.7044592248786605, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 309125153.2350401, 'mu_resolved': 0.030912515323504017, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.3854776419948416, 'V_terminal': 2.3854796419948414, 'power_W': 3.776717231701891e-10, 'counting_converged': True, 'blocked_load_probability': 0.0, 'one_pair_valid': True, 'set_feasible': False, 'set_priced_F_p': 1.0, 'set_E_C_meV': 30.31504311247509, 'set_EC_over_kT': 1.1726382211797872, 'set_radius_nm': 5.0, 'set_radius_max_nm': 0.5863191105898936, 'set_C_sigma_F': 5.285087763377385e-18, 'set_R_T_over_RQ': 38.74045864977526, 'set_f_max_Hz': 18921161667.91825, 'pair_supply_possible': True, 'ideal_load_F_p': 0.0, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 127.25016471751717, 'tau_rad_cavity_ns': 39.41039314440644, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 135.7771191577549, 'applied_field_kVcm': 135.7771191577549, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 3.4126017696258097, 'reservoir_offset_meV': 1297.4304216306182, 'background_acceptance': 0.0465515605136271, 'background_rate_window_s': 4.014083382097431e-34, 'background_window_ns': 0.1, 'eta_background': 0.1, 'optical_reservoir_energy_eV': 3.4126017696258097, 'optical_reservoir_kind': 'gan_barrier', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'isolated_dot', 'effective_height_nm': 3.0, 'effective_radius_nm': 10.0, 'cavity_reference_V_j_V': 2.385477642681538, 'cavity_reference_transition_eV': 2.1151713482398904, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-deshpande2014-comparison-design.yaml', 230.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'rectangular', 'T_hs': 230.0, 'T_j': 230.00000136776424, 'g2_op': 0.9107689718420242, 'rho_pulsed': 0.909172476443062, 'collected_flux_pulsed_s': 0.23784322739983843, 'collected_flux_x_s': 0.23760850552837748, 'collected_flux_xx_s': 0.00023472187146094287, 'background_flux_s': 0.02376085055283775, 'total_detected_flux_s': 0.2616040779526762, 'mean_counts': 1.189216136999192e-08, 'mean_counts_x': 1.1880425276418873e-08, 'mean_counts_xx': 1.1736093573047144e-11, 'E_X_eV': 1.8971966560645306, 'lambda_nm': 653.5126340417861, 'field_kVcm': -6293.479898112817, 'overlap_sq': 0.0740023443460749, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 110.20113633540434, 'k_X_ns': 105.14191050730909, 'k_XX_ns': 210.28382101461818, 'gamma_X0_ns': 0.0740023443460749, 'gamma_XX0_ns': 0.1480046886921498, 'gamma_X_ns': 0.12301726450292626, 'gamma_XX_ns': 0.19539950609353837, 'S_X': 0.0011686443633875554, 'S_XX': 0.0009283553293770895, 'Q': 2000.0, 'kappa_meV': 0.938737623563374, 'detuning_meV': 19.721408937782627, 'Fp_add': 75.99088773175333, 'F_eff_X': 1.6623428026500233, 'F_eff_XX': 1.320225107867832, 'eta_out': 0.1, 'gate_ns_used': 5.0, 'rep_rate_hz': 200000000.0, 'r_dot_s': 309142319.55755633, 'mu_resolved': 0.030914231955755638, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 32.04353268, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 1.9717264989480943, 'V_terminal': 1.9717284989480943, 'power_W': 7.805021025103226e-10, 'counting_converged': True, 'blocked_load_probability': nan, 'one_pair_valid': False, 'set_feasible': False, 'set_priced_F_p': nan, 'set_E_C_meV': nan, 'set_EC_over_kT': nan, 'set_radius_nm': nan, 'set_radius_max_nm': nan, 'set_C_sigma_F': nan, 'set_R_T_over_RQ': nan, 'set_f_max_Hz': nan, 'pair_supply_possible': True, 'ideal_load_F_p': nan, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 13.513085414206072, 'tau_rad_cavity_ns': 8.128940308018413, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 175.9750034053537, 'applied_field_kVcm': 175.9750034053537, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 3.4396357542359723, 'reservoir_offset_meV': 1542.4390981714419, 'background_acceptance': 0.0006816298774212105, 'background_rate_window_s': 2.837855999186949e-61, 'background_window_ns': 0.1, 'eta_background': 0.1, 'optical_reservoir_energy_eV': 3.4396357542359723, 'optical_reservoir_kind': 'gan_barrier', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'isolated_dot', 'effective_height_nm': 2.0, 'effective_radius_nm': 12.5, 'cavity_reference_V_j_V': 1.9076950363847969, 'cavity_reference_transition_eV': 1.8746752471814585, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-deshpande2014-comparison-design.yaml', 300.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'rectangular', 'T_hs': 300.0, 'T_j': 300.0000013643743, 'g2_op': 0.9884633086673388, 'rho_pulsed': 0.9091229055954397, 'collected_flux_pulsed_s': 6.981809444275543, 'collected_flux_x_s': 6.979106478089611, 'collected_flux_xx_s': 0.002702966185933, 'background_flux_s': 0.6979106478089612, 'total_detected_flux_s': 7.679720092084505, 'mean_counts': 3.4909047221377717e-07, 'mean_counts_x': 3.489553239044805e-07, 'mean_counts_xx': 1.3514830929665e-10, 'E_X_eV': 1.8746752467200014, 'lambda_nm': 661.3636074670915, 'field_kVcm': -6294.394296864753, 'overlap_sq': 0.07324514799402232, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 110.18100308909551, 'k_X_ns': 502.28002009939814, 'k_XX_ns': 1004.5600401987963, 'gamma_X0_ns': 0.07324514799402232, 'gamma_XX0_ns': 0.14649029598804464, 'gamma_X_ns': 0.2184197095223622, 'gamma_XX_ns': 0.36817836502160717, 'S_X': 0.0004346674381823319, 'S_XX': 0.00036637279978841194, 'Q': 2000.0, 'kappa_meV': 0.9373376235634417, 'detuning_meV': -4.068820835811948e-07, 'Fp_add': 75.99088773175333, 'F_eff_X': 2.9820365649365317, 'F_eff_XX': 2.513329381569786, 'eta_out': 0.1, 'gate_ns_used': 5.0, 'rep_rate_hz': 200000000.0, 'r_dot_s': 309142295.18966424, 'mu_resolved': 0.030914229518966427, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 32.04353268, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 1.9076950350599584, 'V_terminal': 1.9076970350599585, 'power_W': 7.54971353682064e-10, 'counting_converged': True, 'blocked_load_probability': nan, 'one_pair_valid': False, 'set_feasible': False, 'set_priced_F_p': nan, 'set_E_C_meV': nan, 'set_EC_over_kT': nan, 'set_radius_nm': nan, 'set_radius_max_nm': nan, 'set_C_sigma_F': nan, 'set_R_T_over_RQ': nan, 'set_f_max_Hz': nan, 'pair_supply_possible': True, 'ideal_load_F_p': nan, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 13.652781479554276, 'tau_rad_cavity_ns': 4.578341406033315, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 175.06060465341756, 'applied_field_kVcm': 175.06060465341756, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 3.4126017693403967, 'reservoir_offset_meV': 1537.9265226203954, 'background_acceptance': 0.04135055906578536, 'background_rate_window_s': 3.4439704130636314e-46, 'background_window_ns': 0.1, 'eta_background': 0.1, 'optical_reservoir_energy_eV': 3.4126017693403967, 'optical_reservoir_kind': 'gan_barrier', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'isolated_dot', 'effective_height_nm': 2.0, 'effective_radius_nm': 12.5, 'cavity_reference_V_j_V': 1.9076950363847969, 'cavity_reference_transition_eV': 1.8746752471814585, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-nonpolar-pulse-design.yaml', 230.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'rectangular', 'T_hs': 230.0, 'T_j': 230.00000068034873, 'g2_op': 0.6381740859262495, 'rho_pulsed': 0.9090909091212144, 'collected_flux_pulsed_s': 2.3058096107964176e-06, 'collected_flux_x_s': 2.3058096099508915e-06, 'collected_flux_xx_s': 8.45526008355162e-16, 'background_flux_s': 2.3058096099508915e-07, 'total_detected_flux_s': 2.5363905717915063e-06, 'mean_counts': 2.8822620134955214e-13, 'mean_counts_x': 2.882262012438614e-13, 'mean_counts_xx': 1.0569075104439525e-22, 'E_X_eV': 2.7728756966594266, 'lambda_nm': 447.13219041649717, 'field_kVcm': 136.69397484544533, 'overlap_sq': 0.9745248424410299, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 191.0827864773979, 'k_X_ns': 10.11893153266569, 'k_XX_ns': 20.23786306533138, 'gamma_X0_ns': 0.9745248424410299, 'gamma_XX0_ns': 1.9490496848820598, 'gamma_X_ns': 1.8778874666472583, 'gamma_XX_ns': 2.8411949002781305, 'S_X': 0.15653211628472546, 'S_XX': 0.12310705681799676, 'Q': 2000.0, 'kappa_meV': 1.3761976938475826, 'detuning_meV': 20.48030896426134, 'Fp_add': 75.99088773175333, 'F_eff_X': 1.9269775226493442, 'F_eff_XX': 1.4577334391811854, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 28.38693569399813, 'mu_resolved': 2.8386935693998133e-09, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.451726857949691, 'V_terminal': 2.451728857949691, 'power_W': 3.882347507536855e-10, 'counting_converged': True, 'blocked_load_probability': nan, 'one_pair_valid': False, 'set_feasible': False, 'set_priced_F_p': nan, 'set_E_C_meV': nan, 'set_EC_over_kT': nan, 'set_radius_nm': nan, 'set_radius_max_nm': nan, 'set_C_sigma_F': nan, 'set_R_T_over_RQ': nan, 'set_f_max_Hz': nan, 'pair_supply_possible': True, 'ideal_load_F_p': nan, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 1.0261411063623158, 'tau_rad_cavity_ns': 0.532513272366304, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 136.69397484544533, 'applied_field_kVcm': 136.69397484544533, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 3.4396357544777194, 'reservoir_offset_meV': 666.7600578182928, 'background_acceptance': 0.0013374148418590784, 'background_rate_window_s': 1.4556747923502226e-29, 'background_window_ns': 0.1, 'eta_background': 0.1, 'optical_reservoir_energy_eV': 3.4396357544777194, 'optical_reservoir_kind': 'gan_barrier', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'isolated_dot', 'effective_height_nm': 3.0, 'effective_radius_nm': 10.0, 'cavity_reference_V_j_V': 2.385477642681538, 'cavity_reference_transition_eV': 2.7495953877223793, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-nonpolar-pulse-design.yaml', 300.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'rectangular', 'T_hs': 300.0, 'T_j': 300.00000068252336, 'g2_op': 0.9568992350651753, 'rho_pulsed': 0.9090909091802536, 'collected_flux_pulsed_s': 0.00020879234084069921, 'collected_flux_x_s': 0.00020879234061498036, 'collected_flux_xx_s': 2.257188483549681e-13, 'background_flux_s': 2.0879234061498036e-05, 'total_detected_flux_s': 0.00022967157490219722, 'mean_counts': 2.60990426050874e-11, 'mean_counts_x': 2.6099042576872543e-11, 'mean_counts_xx': 2.821485604437101e-20, 'E_X_eV': 2.7495953874638475, 'lambda_nm': 450.9179749328852, 'field_kVcm': 135.7771191577549, 'overlap_sq': 0.9744872825376502, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 191.20931036672715, 'k_X_ns': 124.55872637302448, 'k_XX_ns': 249.11745274604897, 'gamma_X0_ns': 0.9744872825376502, 'gamma_XX0_ns': 1.9489745650753003, 'gamma_X_ns': 3.7733132809001364, 'gamma_XX_ns': 6.247208632907322, 'S_X': 0.02940273754765918, 'S_XX': 0.02446387295396595, 'Q': 2000.0, 'kappa_meV': 1.3747976938475392, 'detuning_meV': -2.312310343199897e-07, 'Fp_add': 75.99088773175333, 'F_eff_X': 3.8721011023089993, 'F_eff_XX': 3.205382330202937, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 236.17039549591823, 'mu_resolved': 2.3617039549591824e-08, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.3854776419948416, 'V_terminal': 2.3854796419948414, 'power_W': 3.776717231701891e-10, 'counting_converged': True, 'blocked_load_probability': nan, 'one_pair_valid': False, 'set_feasible': False, 'set_priced_F_p': nan, 'set_E_C_meV': nan, 'set_EC_over_kT': nan, 'set_radius_nm': nan, 'set_radius_max_nm': nan, 'set_C_sigma_F': nan, 'set_R_T_over_RQ': nan, 'set_f_max_Hz': nan, 'pair_supply_possible': True, 'ideal_load_F_p': nan, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 1.0261806571717513, 'tau_rad_cavity_ns': 0.26501907622190507, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 135.7771191577549, 'applied_field_kVcm': 135.7771191577549, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 3.4126017696258097, 'reservoir_offset_meV': 663.0063821619623, 'background_acceptance': 0.06015865540469076, 'background_rate_window_s': 1.8440054140575044e-21, 'background_window_ns': 0.1, 'eta_background': 0.1, 'optical_reservoir_energy_eV': 3.4126017696258097, 'optical_reservoir_kind': 'gan_barrier', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'isolated_dot', 'effective_height_nm': 3.0, 'effective_radius_nm': 10.0, 'cavity_reference_V_j_V': 2.385477642681538, 'cavity_reference_transition_eV': 2.7495953877223793, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-nonpolar-set-design.yaml', 230.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'deterministic_pair', 'T_hs': 230.0, 'T_j': 230.00000068034873, 'g2_op': 0.17355371900826433, 'rho_pulsed': 0.9090909090909092, 'collected_flux_pulsed_s': 812.2784498774813, 'collected_flux_x_s': 812.2784498774813, 'collected_flux_xx_s': 0.0, 'background_flux_s': 81.22784498774813, 'total_detected_flux_s': 893.5062948652294, 'mean_counts': 0.00010153480623468517, 'mean_counts_x': 0.00010153480623468517, 'mean_counts_xx': 0.0, 'E_X_eV': 2.7728756966594266, 'lambda_nm': 447.13219041649717, 'field_kVcm': 136.69397484544533, 'overlap_sq': 0.9745248424410299, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 191.0827864773979, 'k_X_ns': 10.11893153266569, 'k_XX_ns': 20.23786306533138, 'gamma_X0_ns': 0.9745248424410299, 'gamma_XX0_ns': 1.9490496848820598, 'gamma_X_ns': 1.8778874666472583, 'gamma_XX_ns': 2.8411949002781305, 'S_X': 0.15653211628472546, 'S_XX': 0.12310705681799676, 'Q': 2000.0, 'kappa_meV': 1.3761976938475826, 'detuning_meV': 20.48030896426134, 'Fp_add': 75.99088773175333, 'F_eff_X': 1.9269775226493442, 'F_eff_XX': 1.4577334391811854, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 28.38693569399813, 'mu_resolved': 2.8386935693998133e-09, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.451726857949691, 'V_terminal': 2.451728857949691, 'power_W': 3.882347507536855e-10, 'counting_converged': True, 'blocked_load_probability': 0.0, 'one_pair_valid': True, 'set_feasible': False, 'set_priced_F_p': 1.0, 'set_E_C_meV': 30.31504311247509, 'set_EC_over_kT': 1.5295281135377246, 'set_radius_nm': 5.0, 'set_radius_max_nm': 0.7647640567688623, 'set_C_sigma_F': 5.285087763377385e-18, 'set_R_T_over_RQ': 38.74045864977526, 'set_f_max_Hz': 18921161667.91825, 'pair_supply_possible': True, 'ideal_load_F_p': 0.0, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 1.0261411063623158, 'tau_rad_cavity_ns': 0.532513272366304, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 136.69397484544533, 'applied_field_kVcm': 136.69397484544533, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 3.4396357544777194, 'reservoir_offset_meV': 666.7600578182928, 'background_acceptance': 0.0013374148418590784, 'background_rate_window_s': 1.4556747923502226e-29, 'background_window_ns': 0.1, 'eta_background': 0.1, 'optical_reservoir_energy_eV': 3.4396357544777194, 'optical_reservoir_kind': 'gan_barrier', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'isolated_dot', 'effective_height_nm': 3.0, 'effective_radius_nm': 10.0, 'cavity_reference_V_j_V': 2.385477642681538, 'cavity_reference_transition_eV': 2.7495953877223793, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-nonpolar-set-design.yaml', 300.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'deterministic_pair', 'T_hs': 300.0, 'T_j': 300.00000068252336, 'g2_op': 0.17355371900826433, 'rho_pulsed': 0.9090909090909092, 'collected_flux_pulsed_s': 8840.749939502472, 'collected_flux_x_s': 8840.749939502472, 'collected_flux_xx_s': 0.0, 'background_flux_s': 884.0749939502472, 'total_detected_flux_s': 9724.824933452717, 'mean_counts': 0.0011050937424378088, 'mean_counts_x': 0.0011050937424378088, 'mean_counts_xx': 0.0, 'E_X_eV': 2.7495953874638475, 'lambda_nm': 450.9179749328852, 'field_kVcm': 135.7771191577549, 'overlap_sq': 0.9744872825376502, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 191.20931036672715, 'k_X_ns': 124.55872637302448, 'k_XX_ns': 249.11745274604897, 'gamma_X0_ns': 0.9744872825376502, 'gamma_XX0_ns': 1.9489745650753003, 'gamma_X_ns': 3.7733132809001364, 'gamma_XX_ns': 6.247208632907322, 'S_X': 0.02940273754765918, 'S_XX': 0.02446387295396595, 'Q': 2000.0, 'kappa_meV': 1.3747976938475392, 'detuning_meV': -2.312310343199897e-07, 'Fp_add': 75.99088773175333, 'F_eff_X': 3.8721011023089993, 'F_eff_XX': 3.205382330202937, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 236.17039549591823, 'mu_resolved': 2.3617039549591824e-08, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.3854776419948416, 'V_terminal': 2.3854796419948414, 'power_W': 3.776717231701891e-10, 'counting_converged': True, 'blocked_load_probability': 0.0, 'one_pair_valid': True, 'set_feasible': False, 'set_priced_F_p': 1.0, 'set_E_C_meV': 30.31504311247509, 'set_EC_over_kT': 1.1726382211797872, 'set_radius_nm': 5.0, 'set_radius_max_nm': 0.5863191105898936, 'set_C_sigma_F': 5.285087763377385e-18, 'set_R_T_over_RQ': 38.74045864977526, 'set_f_max_Hz': 18921161667.91825, 'pair_supply_possible': True, 'ideal_load_F_p': 0.0, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 1.0261806571717513, 'tau_rad_cavity_ns': 0.26501907622190507, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 135.7771191577549, 'applied_field_kVcm': 135.7771191577549, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 3.4126017696258097, 'reservoir_offset_meV': 663.0063821619623, 'background_acceptance': 0.06015865540469076, 'background_rate_window_s': 1.8440054140575044e-21, 'background_window_ns': 0.1, 'eta_background': 0.1, 'optical_reservoir_energy_eV': 3.4126017696258097, 'optical_reservoir_kind': 'gan_barrier', 'reservoir_kind': 'gan_barrier', 'geometry_type': 'isolated_dot', 'effective_height_nm': 3.0, 'effective_radius_nm': 10.0, 'cavity_reference_V_j_V': 2.385477642681538, 'cavity_reference_transition_eV': 2.7495953877223793, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-qw-fluctuation-pulse-design.yaml', 230.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'rectangular', 'T_hs': 230.0, 'T_j': 230.0000006604343, 'g2_op': 0.9971698034700771, 'rho_pulsed': 0.8524074935609247, 'collected_flux_pulsed_s': 5.695222929422368e-05, 'collected_flux_x_s': 5.69500664262736e-05, 'collected_flux_xx_s': 2.162867950083938e-09, 'background_flux_s': 9.861154825977145e-06, 'total_detected_flux_s': 6.681338412020083e-05, 'mean_counts': 7.11902866177796e-12, 'mean_counts_x': 7.1187583032841995e-12, 'mean_counts_xx': 2.7035849376049224e-16, 'E_X_eV': 1.9460116952479543, 'lambda_nm': 637.1194926667814, 'field_kVcm': -3898.5914804067743, 'overlap_sq': 0.0015560238184151272, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 42.889389365046156, 'k_X_ns': 3138.456929072946, 'k_XX_ns': 6276.913858145892, 'gamma_X0_ns': 0.0015560238184151272, 'gamma_XX0_ns': 0.0031120476368302545, 'gamma_X_ns': 0.0024854047671060723, 'gamma_XX_ns': 0.004042548557538472, 'S_X': 7.919187215360824e-07, 'S_XX': 6.440340022754197e-07, 'Q': 2000.0, 'kappa_meV': 0.9623664776139798, 'detuning_meV': 21.278740019994746, 'Fp_add': 75.99088773175333, 'F_eff_X': 1.5972793846032234, 'F_eff_XX': 1.2989995749730778, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 309141449.6379184, 'mu_resolved': 0.030914144963791846, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.380701990860128, 'V_terminal': 2.3807039908601277, 'power_W': 3.7687077201934576e-10, 'counting_converged': True, 'blocked_load_probability': nan, 'one_pair_valid': False, 'set_feasible': False, 'set_priced_F_p': nan, 'set_E_C_meV': nan, 'set_EC_over_kT': nan, 'set_radius_nm': nan, 'set_radius_max_nm': nan, 'set_C_sigma_F': nan, 'set_R_T_over_RQ': nan, 'set_f_max_Hz': nan, 'pair_supply_possible': True, 'ideal_load_F_p': nan, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 642.663684299216, 'tau_rad_cavity_ns': 402.3489506557794, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 139.94257222644688, 'applied_field_kVcm': 139.94257222644688, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 2.1502906551273027, 'reservoir_offset_meV': 204.2789598793484, 'background_acceptance': 0.000598193317399888, 'background_rate_window_s': 8.705689411280952, 'background_window_ns': 0.1, 'eta_background': 0.1, 'optical_reservoir_energy_eV': 2.1502906551273027, 'optical_reservoir_kind': 'ingan_qw', 'reservoir_kind': 'mixed', 'geometry_type': 'qw_fluctuation', 'effective_height_nm': 3.5, 'effective_radius_nm': 20.0, 'cavity_reference_V_j_V': 2.292836535807325, 'cavity_reference_transition_eV': 1.921932955254377, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-qw-fluctuation-pulse-design.yaml', 300.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'rectangular', 'T_hs': 300.0, 'T_j': 300.0000006557362, 'g2_op': 0.9998097399428482, 'rho_pulsed': 0.42442975389363125, 'collected_flux_pulsed_s': 0.004517145671265698, 'collected_flux_x_s': 0.004517013674792185, 'collected_flux_xx_s': 1.3199647351304614e-07, 'background_flux_s': 0.006125712492720999, 'total_detected_flux_s': 0.010642858163986697, 'mean_counts': 5.646432089082122e-10, 'mean_counts_x': 5.646267093490231e-10, 'mean_counts_xx': 1.6499559189130768e-14, 'E_X_eV': 1.9219329550210902, 'lambda_nm': 645.1015789915494, 'field_kVcm': -3897.583864360978, 'overlap_sq': 0.0015344407420688783, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 42.910444857839515, 'k_X_ns': 6777.063184672321, 'k_XX_ns': 13554.126369344642, 'gamma_X0_ns': 0.0015344407420688783, 'gamma_XX0_ns': 0.0030688814841377567, 'gamma_X_ns': 0.00465037667413206, 'gamma_XX_ns': 0.007828534040646446, 'S_X': 6.861930243757757e-07, 'S_XX': 5.775753675119778e-07, 'Q': 2000.0, 'kappa_meV': 0.9609664776140737, 'detuning_meV': -2.0705726022640647e-07, 'Fp_add': 75.99088773175333, 'F_eff_X': 3.030665536071456, 'F_eff_XX': 2.5509404912213407, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 309125153.2350552, 'mu_resolved': 0.030912515323505523, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.2928365349450868, 'V_terminal': 2.2928385349450866, 'power_W': 3.6284914604221315e-10, 'counting_converged': True, 'blocked_load_probability': nan, 'one_pair_valid': False, 'set_feasible': False, 'set_priced_F_p': nan, 'set_E_C_meV': nan, 'set_EC_over_kT': nan, 'set_radius_nm': nan, 'set_radius_max_nm': nan, 'set_C_sigma_F': nan, 'set_R_T_over_RQ': nan, 'set_f_max_Hz': nan, 'pair_supply_possible': True, 'ideal_load_F_p': nan, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 651.7032379182693, 'tau_rad_cavity_ns': 215.0363443809933, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 140.95018827224288, 'applied_field_kVcm': 140.95018827224288, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 2.1261495373577244, 'reservoir_offset_meV': 204.2165823366342, 'background_acceptance': 0.04237442164433611, 'background_rate_window_s': 167.3772438968552, 'background_window_ns': 0.1, 'eta_background': 0.1, 'optical_reservoir_energy_eV': 2.1261495373577244, 'optical_reservoir_kind': 'ingan_qw', 'reservoir_kind': 'mixed', 'geometry_type': 'qw_fluctuation', 'effective_height_nm': 3.5, 'effective_radius_nm': 20.0, 'cavity_reference_V_j_V': 2.292836535807325, 'cavity_reference_transition_eV': 1.921932955254377, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-qw-fluctuation-set-design.yaml', 230.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'deterministic_pair', 'T_hs': 230.0, 'T_j': 230.0000006604343, 'g2_op': 0.1769411420851732, 'rho_pulsed': 0.9072259133836659, 'collected_flux_pulsed_s': 0.0018423818136572687, 'collected_flux_x_s': 0.0018423818136572687, 'collected_flux_xx_s': 0.0, 'background_flux_s': 0.00018840432954907667, 'total_detected_flux_s': 0.0020307861432063454, 'mean_counts': 2.3029772670715858e-10, 'mean_counts_x': 2.3029772670715858e-10, 'mean_counts_xx': 0.0, 'E_X_eV': 1.9460116952479543, 'lambda_nm': 637.1194926667814, 'field_kVcm': -3898.5914804067743, 'overlap_sq': 0.0015560238184151272, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 42.889389365046156, 'k_X_ns': 3138.456929072946, 'k_XX_ns': 6276.913858145892, 'gamma_X0_ns': 0.0015560238184151272, 'gamma_XX0_ns': 0.0031120476368302545, 'gamma_X_ns': 0.0024854047671060723, 'gamma_XX_ns': 0.004042548557538472, 'S_X': 7.919187215360824e-07, 'S_XX': 6.440340022754197e-07, 'Q': 2000.0, 'kappa_meV': 0.9623664776139798, 'detuning_meV': 21.278740019994746, 'Fp_add': 75.99088773175333, 'F_eff_X': 1.5972793846032234, 'F_eff_XX': 1.2989995749730778, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 309141449.6379184, 'mu_resolved': 0.030914144963791846, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.380701990860128, 'V_terminal': 2.3807039908601277, 'power_W': 3.7687077201934576e-10, 'counting_converged': True, 'blocked_load_probability': 0.0, 'one_pair_valid': True, 'set_feasible': False, 'set_priced_F_p': 1.0, 'set_E_C_meV': 30.31504311247509, 'set_EC_over_kT': 1.529528113670158, 'set_radius_nm': 5.0, 'set_radius_max_nm': 0.7647640568350789, 'set_C_sigma_F': 5.285087763377385e-18, 'set_R_T_over_RQ': 38.74045864977526, 'set_f_max_Hz': 18921161667.91825, 'pair_supply_possible': True, 'ideal_load_F_p': 0.0, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 642.663684299216, 'tau_rad_cavity_ns': 402.3489506557794, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 139.94257222644688, 'applied_field_kVcm': 139.94257222644688, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 2.1502906551273027, 'reservoir_offset_meV': 204.2789598793484, 'background_acceptance': 0.000598193317399888, 'background_rate_window_s': 8.705689411280952, 'background_window_ns': 0.1, 'eta_background': 0.1, 'optical_reservoir_energy_eV': 2.1502906551273027, 'optical_reservoir_kind': 'ingan_qw', 'reservoir_kind': 'mixed', 'geometry_type': 'qw_fluctuation', 'effective_height_nm': 3.5, 'effective_radius_nm': 20.0, 'cavity_reference_V_j_V': 2.292836535807325, 'cavity_reference_transition_eV': 1.921932955254377, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
('nitride-qw-fluctuation-set-design.yaml', 300.0): {'platform': 'ingan_gan_planar', 'cycle_loading': 'deterministic_pair', 'T_hs': 300.0, 'T_j': 300.0000006557362, 'g2_op': 0.22894890250797317, 'rho_pulsed': 0.8780951528689968, 'collected_flux_pulsed_s': 0.1461291442822906, 'collected_flux_x_s': 0.1461291442822906, 'collected_flux_xx_s': 0.0, 'background_flux_s': 0.02028692555347084, 'total_detected_flux_s': 0.16641606983576146, 'mean_counts': 1.8266143035286323e-08, 'mean_counts_x': 1.8266143035286323e-08, 'mean_counts_xx': 0.0, 'E_X_eV': 1.9219329550210902, 'lambda_nm': 645.1015789915494, 'field_kVcm': -3897.583864360978, 'overlap_sq': 0.0015344407420688783, 'electron_bound': True, 'hole_bound': True, 'E_a_meV': 42.910444857839515, 'k_X_ns': 6777.063184672321, 'k_XX_ns': 13554.126369344642, 'gamma_X0_ns': 0.0015344407420688783, 'gamma_XX0_ns': 0.0030688814841377567, 'gamma_X_ns': 0.00465037667413206, 'gamma_XX_ns': 0.007828534040646446, 'S_X': 6.861930243757757e-07, 'S_XX': 5.775753675119778e-07, 'Q': 2000.0, 'kappa_meV': 0.9609664776140737, 'detuning_meV': -2.0705726022640647e-07, 'Fp_add': 75.99088773175333, 'F_eff_X': 3.030665536071456, 'F_eff_XX': 2.5509404912213407, 'eta_out': 0.1, 'gate_ns_used': 12.5, 'rep_rate_hz': 80000000.0, 'r_dot_s': 309125153.2350552, 'mu_resolved': 0.030912515323505523, 'n_dot_cm2_used': 10000000000.0, 'tau_cap_ps_used': 10.0, 'I_pair_pA': 12.817413071999999, 'transport_current_uA': 0.02, 'resolved_current_uA': 0.02, 'V_j': 2.2928365349450868, 'V_terminal': 2.2928385349450866, 'power_W': 3.6284914604221315e-10, 'counting_converged': True, 'blocked_load_probability': 0.0, 'one_pair_valid': True, 'set_feasible': False, 'set_priced_F_p': 1.0, 'set_E_C_meV': 30.31504311247509, 'set_EC_over_kT': 1.1726382212844928, 'set_radius_nm': 5.0, 'set_radius_max_nm': 0.5863191106422464, 'set_C_sigma_F': 5.285087763377385e-18, 'set_R_T_over_RQ': 38.74045864977526, 'set_f_max_Hz': 18921161667.91825, 'pair_supply_possible': True, 'ideal_load_F_p': 0.0, 'valid': True, 'invalid_reasons': [], 'provenance': {'nitride': '[A/E/DR] planar integration; cavity/transport inputs retain module provenance'}, 'tau_rad_bare_ns': 651.7032379182693, 'tau_rad_cavity_ns': 215.0363443809933, 'spectroscopy_valid': True, 'spectroscopy_invalid_reasons': [], 'temperature_mode': 'self_consistent', 'evaluation_kind': 'source', 'field_polarity': 1, 'diode_field_kVcm': 140.95018827224288, 'applied_field_kVcm': 140.95018827224288, 'flat_band': False, 'depletion_regime': 'depleted', 'reservoir_energy_eV': 2.1261495373577244, 'reservoir_offset_meV': 204.2165823366342, 'background_acceptance': 0.04237442164433611, 'background_rate_window_s': 167.3772438968552, 'background_window_ns': 0.1, 'eta_background': 0.1, 'optical_reservoir_energy_eV': 2.1261495373577244, 'optical_reservoir_kind': 'ingan_qw', 'reservoir_kind': 'mixed', 'geometry_type': 'qw_fluctuation', 'effective_height_nm': 3.5, 'effective_radius_nm': 20.0, 'cavity_reference_V_j_V': 2.292836535807325, 'cavity_reference_transition_eV': 1.921932955254377, 'cavity_reference_convention': 'current_controlled', 'device_pass': False},
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
# HIGH 1 (fix-3): shell_multiplier is now stated EXPLICITLY on both cards
# (via card()'s own shell_multiplier= knob, default 1.0) so this exercises
# the real card shape the fix-3 guard checks, not an unstated class default.
_s_none = evaluate(card(family="horizontal_as_built", regime="rectangular", strain_bound="relaxed"))["scalars"]
_s_algan = evaluate(card(family="horizontal_as_built", regime="rectangular", strain_bound="relaxed",
                          outer=15.5, extra_nanowire={"shell": "AlGaN"}, shell_multiplier=0.1))["scalars"]
check("HIGH 3: AlGaN shell changes shell_multiplier_used vs shell='none'",
      not close(_s_algan["shell_multiplier_used"], _s_none["shell_multiplier_used"], rel=1e-9))
check("HIGH 3: AlGaN shell changes k_X_ns vs shell='none'",
      not close(_s_algan["k_X_ns"], _s_none["k_X_ns"], rel=1e-6))
check("HIGH 3: k_side_ns itself is geometry-only (unaffected by the shell factor, only the occupied-dot/reservoir channels are)",
      close(_s_algan["k_side_ns"], _s_none["k_side_ns"], rel=1e-9))

# HIGH 1 (fix-3): nitride.surface.shell_multiplier is cross-checked against
# the declared shell instead of silently accepted regardless of value.
check("HIGH 1: shell='AlGaN' with shell_multiplier=1.0 (no passivation effect) is rejected",
      raises(lambda: evaluate(card(family="horizontal_as_built", regime="rectangular", strain_bound="relaxed",
                                    outer=15.5, extra_nanowire={"shell": "AlGaN"}, shell_multiplier=1.0))))
check("HIGH 1: shell='none' with shell_multiplier!=1.0 is rejected",
      raises(lambda: evaluate(card(family="horizontal_as_built", regime="rectangular", strain_bound="relaxed",
                                    shell_multiplier=0.5))))
check("HIGH 1: AlGaN shell with the contract's 0.1 default multiplier gives k_X_ns ~= 0.0112",
      close(_s_algan["k_X_ns"], 0.0112, rel=0, abs_=2e-4))
check("HIGH 1: shell_multiplier_used reports the effective 0.1 for the AlGaN card",
      close(_s_algan["shell_multiplier_used"], 0.1, rel=1e-9))

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

# MEDIUM 4 (fix-3): nitride.reservoir_energy_eV is an allow-listed top-level
# override, honored through the SAME planar helper/"background" dict
# convention device.py's own nitride path uses -- previously dropped
# unconditionally (background={} regardless of the card).
_d_res_override = card(family="horizontal_as_built", regime="rectangular", strain_bound="relaxed",
                        extra_nitride={"reservoir_energy_eV": 2.9})
_s_res_override = evaluate(_d_res_override)["scalars"]
check("MEDIUM 4: nitride.reservoir_energy_eV overrides the reservoir energy used by evaluate_injection (accepted_background_s moves)",
      not close(_s_res_override["accepted_background_s"], hz["accepted_background_s"], rel=1e-6))
check("MEDIUM 4: an invalid (non-positive) nitride.reservoir_energy_eV yields an explicit invalid row, not a silent value",
      evaluate(card(extra_nitride={"reservoir_energy_eV": -1.0}))["scalars"]["valid"] is False)

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

# Coordinator addition (cards review, fix-3): NitrideNanowirePhotonicsParams
# requires dipole_weights to be an actual Python tuple, but a YAML/JSON
# card can only express a list -- coerce a 3-element list/tuple to a tuple
# of floats at the card boundary so the contract's cplane_only sensitivity
# (dipole_weights=[0.0, 0.5, 0.5]) is reachable from a card at all, and
# reject a wrong-length list/tuple by name.
_s_dipole_default = evaluate(card(family="horizontal_as_built", regime="rectangular", strain_bound="relaxed"))["scalars"]
_s_dipole_cplane = evaluate(card(family="horizontal_as_built", regime="rectangular", strain_bound="relaxed",
                                  extra_photonics={"dipole_weights": [0.0, 0.5, 0.5]}))["scalars"]
check("dipole_weights: a 3-element list from a card is accepted (coerced to a tuple), not rejected as 'must be a 3-tuple'",
      _s_dipole_cplane["valid"] is True)
check("dipole_weights=[0.0, 0.5, 0.5] moves gamma_X_ns relative to the isotropic default",
      not close(_s_dipole_cplane["gamma_X_ns"], _s_dipole_default["gamma_X_ns"], rel=1e-6))
check("dipole_weights=[0.0, 0.5, 0.5] moves degree_of_linear_polarization relative to the isotropic default",
      not close(_s_dipole_cplane["degree_of_linear_polarization"], _s_dipole_default["degree_of_linear_polarization"], rel=1e-6))
check("dipole_weights: a wrong-length list (2 elements) is rejected with a named error",
      raises(lambda: evaluate(card(extra_photonics={"dipole_weights": [0.3, 0.7]}))))

# LOW 8 (fix-3, structural): a non-finite signal must report a NaN flux, not
# max(0.0, nan)==0.0 (Python's nan-comparisons-always-False behavior would
# otherwise silently report "zero photons collected"). Unreachable through
# a full row today (a non-finite signal already makes rho/g2 non-finite,
# which independently marks the row invalid and NaN-fills every column
# through the invalid-row path -- confirmed by the review), so this is a
# structural source check of the specific guard, not a row-level probe.
check("LOW 8: collected_flux_pulsed_s uses an explicit NaN guard instead of a bare max(0.0, signal)",
      "_nan() if not math.isfinite(signal) else rep_v * max(0.0, signal)" in DEVICE_SRC)

# Audit single-selection fix (2026-09-23; re-pins the LOW 6 fix-3 block,
# which applied a SECOND spectral window at the reservoir offset on top of
# transport.xi_window's X-centred window). The acceptance applied to the
# already-windowed reservoir background is the flat-spectrum average of the
# Lorentzian kernel (kappa = w [A]) over the X-centred window of width w,
# closed form [DR] 0.5*(atan(1 - 2 dx/w) + atan(1 + 2 dx/w)) -- a number
# computed here from arctan alone, not by the module's quad.
_accept_closed = 0.5 * (math.atan(1.0 - 2.0 * 0.5 / 1.0) + math.atan(1.0 + 2.0 * 0.5 / 1.0))
check("single selection: _background_acceptance(w=1, dx=0.5) equals the X-centred closed form 0.5*(atan(0)+atan(2))",
      close(dev._background_acceptance(1.0, 0.5), _accept_closed, rel=1e-9))
check("single selection: _background_acceptance(w, dx=0) equals pi/4 for any w (flat reservoir, X-centred window)",
      all(close(dev._background_acceptance(w_, 0.0), math.pi / 4.0, rel=1e-9) for w_ in (0.3, 1.0, 7.0)))
# Row-level: with b_res = 0 the background is the reservoir channel alone,
# so background_flux_s must equal rep * accepted_background_s * window *
# (closed-form acceptance) * eta_collection_X EXACTLY once -- for two
# different reservoir energies (the acceptance must not depend on the
# reservoir offset, since xi_window already selected the window).
for _res_override in (None, 2.9):
    _d_sel = card(family="horizontal_as_built", regime="rectangular", strain_bound="relaxed",
                  extra_drive_kw={"b_res": 0.0},
                  extra_nitride=({"reservoir_energy_eV": _res_override} if _res_override else None))
    _d_sel.filter.enabled = True; _d_sel.filter.auto_w = False; _d_sel.filter.w = 1.0; _d_sel.filter.dx = 0.5
    _s_sel = evaluate(_d_sel)["scalars"]
    _gate_sel = min(_s_sel["gate_ns_used"], _d_sel.drive.diode["tau_pulse_ns"])
    _bg_sel = _s_sel["rep_rate_hz"] * _s_sel["accepted_background_s"] * _gate_sel * 1e-9 * _accept_closed * _s_sel["eta_collection_X"]
    check("single selection: background_flux_s (b_res=0, reservoir %s) == rep*accepted_background_s*window*A_X-centred*eta_collection_X"
          % ("default" if _res_override is None else "2.9 eV"),
          _s_sel["valid"] is True and _s_sel["background_flux_s"] > 0.0 and close(_s_sel["background_flux_s"], _bg_sel, rel=1e-9, abs_=0.0))
check("single selection: no reservoir-offset acceptance call remains in the device module",
      "reservoir_offset_meV" not in DEVICE_SRC)
check("with the filter disabled, eta_out is full transmission (1.0)",
      close(_s_filt_off["eta_out"], 1.0))

# Audit H1 (2026-09-23) NEW EXPECTATION: capture into the dot does not
# depend on the electron-hole overlap; the overlap sets only the radiative
# rate. Two device evaluations differing ONLY in overlap_sq (levels()
# wrapped to scale overlap_sq by 0.5, everything else bit-identical) must
# give different gamma_X_ns and identical r_captured_s (rel 1e-12).
_orig_levels = dev.levels
def _levels_half_overlap(system, T_K=300.0, **kw):
    lv_ = _orig_levels(system, T_K, **kw)
    return _replace(lv_, overlap_sq=0.5 * lv_.overlap_sq) if lv_.valid else lv_
_d_h1 = card(family="horizontal_as_built", regime="rectangular", strain_bound="relaxed")
_s_h1_ref = evaluate(_d_h1)["scalars"]
try:
    dev.levels = _levels_half_overlap
    _s_h1_half = evaluate(card(family="horizontal_as_built", regime="rectangular", strain_bound="relaxed"))["scalars"]
finally:
    dev.levels = _orig_levels
check("H1: halving ONLY the e-h overlap halves gamma_X0_ns (overlap enters the radiative rate)",
      close(_s_h1_half["gamma_X0_ns"], 0.5 * _s_h1_ref["gamma_X0_ns"], rel=1e-12) and not close(_s_h1_half["gamma_X_ns"], _s_h1_ref["gamma_X_ns"], rel=1e-6))
check("H1: halving ONLY the e-h overlap leaves r_captured_s unchanged (rel 1e-12)",
      _s_h1_ref["valid"] is True and close(_s_h1_half["r_captured_s"], _s_h1_ref["r_captured_s"], rel=1e-12, abs_=0.0))
check("H1: the device no longer passes the e-h overlap as S_dot",
      "S_dot=lv.overlap_sq" not in DEVICE_SRC and "S_dot=lv2.overlap_sq" not in DEVICE_SRC)

# Audit M4 (2026-09-23) NEW EXPECTATION: a row whose V_j exceeds the
# ceiling min(V_bi, E_g(GaN,T)/q) is reported invalid with an out_of_model
# reason. At 10 uA / 300 K on the as-built horizontal wire the ideal SRH
# diode gives V_j above E_g/q (audit: 3.56 V vs 3.44 V). Rth_K_W = 0 and
# the 1e6 ohm designed contact isolate the junction ceiling from thermal
# runaway (which otherwise invalidates this current first).
_s_m4 = evaluate(card(family="horizontal_as_built", regime="rectangular", strain_bound="relaxed", I_uA=10.0,
                      R_s_ohm=1e6, extra_thermal={"Rth_K_W": 0.0}))["scalars"]
check("M4: the same card at the default 0.02 uA stays valid (ceiling not hit)",
      evaluate(card(family="horizontal_as_built", regime="rectangular", strain_bound="relaxed",
                    R_s_ohm=1e6, extra_thermal={"Rth_K_W": 0.0}))["scalars"]["valid"] is True)
check("M4: a 10 uA / 300 K row (V_j above min(V_bi, E_g/q)) is invalid with an out_of_model reason",
      _s_m4["valid"] is False and any("out_of_model" in str(x) for x in _s_m4["invalid_reasons"]))


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
_vert_outer = card(family="vertical_photonic", regime="deterministic_pair", outer=83.0, extra_nanowire={"shell": "AlGaN"}, shell_multiplier=0.1)
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

# LOW 12 (fix-3, documentation/structural): the injector is deliberately
# priced at the transport depletion field ALONE (DEVICE-PIECE addendum),
# unlike sysB's level feedback which adds ext_field_kVcm on top -- confirm
# the injector call site does not also add ext_field_kVcm.
check("LOW 12: injector_feasibility's field_kVcm argument is the depletion field alone (not ext_field_kVcm + depletion_field_kVcm)",
      "field_kVcm=inj[\"depletion_field_kVcm\"],\n            gate_ns=" in DEVICE_SRC)


# --------------------------------------------------------------- 6. invalid-row equality, lifetimes, curves

def bad_pulse_card():
    return card(extra_diode={"tau_pulse_ns": 0.0})


# LOW 9 (fix-3): the pulse==period boundary itself (not just pulse<=0) was
# untested -- LOW 15 (fix-2 round 2) made this platform's loading window
# strictly < the period (stricter than the planar `>`), but nothing
# exercised the equality boundary specifically.
def bad_pulse_equals_period():
    return card(extra_diode={"tau_pulse_ns": 1e9 / 80e6})  # == the default period exactly


check("LOW 9: tau_pulse_ns == period exactly is rejected as invalid (not accepted at the boundary)",
      evaluate(bad_pulse_equals_period())["scalars"]["valid"] is False)

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

# MEDIUM 3 (fix-3), section 2: a broken row-key fixture raises a NAMED
# error identifying the fixture (not a bare traceback pointing at
# _validate()/_one_raw() with no hint that the FIXTURE, not a real card,
# needs fixing). Placed here (after card() is defined) rather than beside
# the doc-parsing MEDIUM 3 checks above.
_orig_fixture_design = dev._row_key_fixture_design
_orig_row_keys_cache = dev._ROW_KEYS_CACHE
try:
    dev._row_key_fixture_design = lambda: card(extra_dot={"height_nm": -1.0})  # malformed: fails card-shape validation
    dev._ROW_KEYS_CACHE = None
    _fixture_err = _catch(dev._row_keys)
    check("MEDIUM 3: a broken row-key fixture names _row_key_fixture_design in the raised error",
          _fixture_err is not None and "_row_key_fixture_design" in str(_fixture_err))
finally:
    dev._row_key_fixture_design = _orig_fixture_design
    dev._ROW_KEYS_CACHE = _orig_row_keys_cache

# MEDIUM 2 (fix-3): the two HIGH-1 key-SET checks above are tautologies by
# construction for VALUE regressions -- dev._row_keys() is contract-columns
# UNION fixture-keys, and _one() always pre-fills that exact key set before
# merging in the raw success-path dict, so dropping a genuinely-computed
# column from the row literal just leaves it NaN under the SAME key (the
# key set never shrinks) -- invisible to a key-set-only check. Replaced
# with (kept alongside the still-legitimate key-set checks, which do catch
# a row producing an UNEXPECTED extra key) a VALUE-level check: every
# contract column the contract marks as always-present must be non-NaN on
# a genuinely valid row, with three named, independently-confirmed exceptions
# from dependency modules this piece does not own:
#   - blocked_load_probability: only meaningful for deterministic_pair (SET)
#     rows; pulse_counting.pulse_g2 does not return it at all.
#   - rti_required_bias_shift_meV: piece 6's own diagnostic (nitride_
#     nanowire_injector.py), NaN whenever no bias-shift candidate applies.
#   - rti_bypass_fraction_tsai_partition: piece 6's NON-GATING sensitivity
#     (the same stack under the alternative Tsai & Bayram 0.30 eV valence
#     partition). Audit L6 fix (2026-09-23): it is now a NaN-safe max over
#     the two carriers, so it is NaN whenever either carrier's resonance is
#     unresolvable under that alternative partition (on the vertical rows
#     the Tsai-partition hole resonance does not exist at the ~80 kV/cm
#     depletion field); the old plain max() hid that by returning the
#     electron value. Pinned behaviourally in verify_nitride_nanowire_
#     injector.py (L6 check). It never gates `valid`. Checker condition
#     (2026-09-23): the allowance is restricted to the vertical family, the
#     only rows where the Tsai-partition hole resonance is absent, so a NaN
#     on a horizontal row fails this guard.
# Strain-mass audit (2026-09-23): the injector's resonance finder now
# rejects shallow non-resonant transmission bumps (peak/valley < 2); under
# the Tsai partition the HORIZONTAL rows' hole path then has no resolvable
# resonance either (the ~54 meV "resonance" accepted before was a 1.06x
# ripple). A horizontal NaN is allowed only where the production finder,
# run here on that row's own card injector params at that row's own
# depletion field, confirms the Tsai-partition hole path is unresolved.
import fsim_core.nitride_nanowire_injector as _inj_mod  # noqa: E402
from fsim_core import nitride_nanowire_device as _nwdev_mod  # noqa: E402
from dataclasses import replace  # noqa: E402


def _tsai_hole_unresolved(s):
    _d_t = card(family=s["family"], regime=s["cycle_loading"] if s["cycle_loading"] == "deterministic_pair" else "rectangular",
                strain_bound=s["strain_bound"])
    _ip = _nwdev_mod._params(_d_t.nitride["injector"], _inj_mod.NitrideNanowireInjectorParams)
    _ip = replace(_ip, delta_Ev_GaN_AlN_eV=_inj_mod._TSAI_BAYRAM_DELTA_EV_EV)
    return _inj_mod._find_resonance(_ip, "hole", 0.0, s["field_kVcm"]) is None


_LEGITIMATELY_NAN_ON_VALID_ROWS = {
    "blocked_load_probability": lambda s: s["cycle_loading"] != "deterministic_pair",
    "rti_required_bias_shift_meV": lambda s: True,
    "rti_bypass_fraction_tsai_partition": lambda s: s["family"] == "vertical_photonic" or _tsai_hole_unresolved(s),
}
check("MEDIUM 2: the Tsai-partition NaN allowance covers vertical rows, and horizontal rows only where the finder confirms the hole path unresolved",
      _LEGITIMATELY_NAN_ON_VALID_ROWS["rti_bypass_fraction_tsai_partition"]({"family": "vertical_photonic"})
      and _inj_mod._find_resonance(replace(_inj_mod.NitrideNanowireInjectorParams(),
                                           delta_Ev_GaN_AlN_eV=0.70), "hole", 0.0, 0.0) is not None)
_valid_scalars_all = [r["scalars"] for r in rows.values() if r["scalars"]["valid"]]
_unexpected_nan = set()
for _s in _valid_scalars_all:
    for _col in CONTRACT_COLUMNS:
        _v = _s.get(_col)
        if isinstance(_v, float) and math.isnan(_v):
            _allow = _LEGITIMATELY_NAN_ON_VALID_ROWS.get(_col)
            if _allow is None or not _allow(_s):
                _unexpected_nan.add((_col, _s["family"], _s["cycle_loading"], _s["strain_bound"]))
check("MEDIUM 2: every always-present contract column is genuinely non-NaN on a valid row (three named, dependency-owned exceptions aside)",
      not _unexpected_nan)
if _unexpected_nan:
    print("  unexpectedly NaN:", sorted(_unexpected_nan))

# MEDIUM 2: dropping a genuinely-computed column from _one_raw's success
# path (monkeypatched) IS caught by the value-level guard above, even
# though it is invisible to the key-set checks (the key stays, only the
# value goes NaN).
_orig_one_raw = dev._one_raw
_DROPPED_KEYS = ("mu", "V_terminal", "eta_inj", "degree_of_linear_polarization")


def _dropping_one_raw(d, T_hs):
    raw = _orig_one_raw(d, T_hs)
    if raw.get("valid"):
        for _k in _DROPPED_KEYS:
            raw.pop(_k, None)
    return raw


dev._one_raw = _dropping_one_raw
try:
    _dropped_row = dev._one(card(family="horizontal_as_built", regime="deterministic_pair", strain_bound="relaxed"), 300.0)
finally:
    dev._one_raw = _orig_one_raw
_dropped_still_present = set(_DROPPED_KEYS) <= set(_dropped_row)
_dropped_all_nan = all(isinstance(_dropped_row.get(k), float) and math.isnan(_dropped_row[k]) for k in _DROPPED_KEYS)
check("MEDIUM 2: dropping mu/V_terminal/eta_inj/degree_of_linear_polarization from _one_raw leaves the key set unchanged (the tautological check's blind spot)",
      _dropped_row["valid"] is True and _dropped_still_present and _dropped_all_nan)
check("MEDIUM 2: the SAME drop is caught by the value-level guard as unexpectedly-NaN",
      set(_DROPPED_KEYS) <= {c for c in CONTRACT_COLUMNS
                              if isinstance(_dropped_row.get(c), float) and math.isnan(_dropped_row[c])
                              and c not in _LEGITIMATELY_NAN_ON_VALID_ROWS})

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

# MEDIUM 5 (fix-3): the catch-all exception handler keeps every coordinate
# already genuinely resolved before the exception hit (T_j,
# thermal_iterations, thermal_converged, field_kVcm), instead of dropping
# them to the row-level NaN/False defaults regardless of progress. x_in=1.0
# (pure InN) converges the thermal solve at T_j~325.2 K, then fails later
# (photonics response()'s si_complex_index lambda-range check) -- the
# review's own confirmed repro of the previous round's bug.
d_x_in_fail = card(family="horizontal_as_built", regime="rectangular", strain_bound="relaxed", x_in=1.0)
s_x_in_fail = evaluate(d_x_in_fail)["scalars"]
check("MEDIUM 5: a chain failure after the thermal solve still reports the real T_j (not NaN)",
      s_x_in_fail["valid"] is False and math.isfinite(s_x_in_fail["T_j"])
      and close(s_x_in_fail["T_j"], 325.2, rel=0, abs_=0.1))
check("MEDIUM 5: the same row keeps thermal_iterations/thermal_converged from the successful thermal solve",
      s_x_in_fail["thermal_converged"] is True and s_x_in_fail["thermal_iterations"] > 0)
check("MEDIUM 5: the same row keeps field_kVcm from the successful injection solve (failure is later, in photonics)",
      math.isfinite(s_x_in_fail["field_kVcm"]))
check("MEDIUM 5: invalid_reasons names the actual downstream failure (photonics lambda range), not a generic message",
      any("si_complex_index" in r or "lambda_nm" in r for r in s_x_in_fail["invalid_reasons"]))

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
# LOW 7 (fix-3): removed the unconditional check(..., True) else-branch
# (a silent pass if V_number ever dropped at/below the cutoff for this
# configuration) -- the precondition is asserted explicitly instead, so a
# future photonics change that moved V_number below the cutoff here would
# FAIL this check loudly rather than silently taking the informational
# branch.
s_wide = evaluate(card(family="vertical_photonic", regime="rectangular", core=120.0, outer=120.0, disc=12.5))["scalars"]
check("headline_eligible present as a bool for vertical rows", isinstance(s_wide.get("headline_eligible"), bool))
check("LOW 7: this configuration's V_number is above the LP11 cutoff (test precondition)",
      s_wide["V_number"] > 2.404826)
check("vertical row above the LP11 cutoff (single_mode False) is never headline_eligible",
      s_wide["single_mode"] is False and s_wide["headline_eligible"] is False)


# --------------------------------------------------------------- injector import wiring

# MEDIUM 10 (fix-2 round 2): the import-retry scaffolding (_injector_symbols,
# six attempts with sleep(0.5)) was removed -- this module now imports
# NitrideNanowireInjectorParams/injector_feasibility as a plain top-level
# import, same as its other five dependency modules.
check("MEDIUM 10: injector_feasibility is imported as a plain top-level symbol (no retry scaffolding)",
      callable(dev.injector_feasibility) and "_injector_symbols" not in DEVICE_SRC and "importlib" not in DEVICE_SRC)



# --------------------------------------------------- audit Phase C item 2 (optional Lindblad path)
# DeviceDesign.quantum (default {}, OFF) adds the cascaded-filter CW g2(0)
# and the HOM indistinguishability as NEW quantum_* keys only. Expectations
# are closed forms the Lindblad code did not produce: I = Gamma/(Gamma +
# 2 gamma*) = hbar Gamma / FWHM (Grange et al., PRL 114, 193601 (2015),
# Eq. 1; Gamma = the row's own gamma_X_ns + k_X_ns) and the two-level
# filtered g2_f(0) = 2D(w+L)/((D+w)(3w+L)) [DR] (verify_lindblad e2), which
# is D/(D+w) at w = L.
import fsim_core.device as _devmod
from fsim_core.device import DeviceDesign as _DD
_HBAR_V = 6.582119569e-4   # [V] CODATA 2018 hbar in meV ns (Tiesinga et al., RMP 93, 025010 (2021))


def _canon_q(x):
    if isinstance(x, dict):
        return {k: _canon_q(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_canon_q(v) for v in x]
    if hasattr(x, "dtype") and hasattr(x, "ravel"):
        return ["ndarray", [repr(v) for v in x.ravel().tolist()]]
    if isinstance(x, float):
        return ["float", repr(x)]
    return [type(x).__name__, repr(x)]


def _strip_q(result):
    out = copy.deepcopy(result)
    sc = out["scalars"]
    for k in [k for k in sc if k.startswith("quantum_")]:
        del sc[k]
    if isinstance(sc.get("provenance"), dict):
        sc["provenance"].pop("quantum", None)
    return out


_nwq = _DD.load(str(ROOT / "cards" / "nitride-nanowire-horizontal-pulse-design.yaml"))
_nw_off = evaluate(_nwq, [_nwq.thermal.T_hs])
_nwq_on = copy.deepcopy(_nwq); _nwq_on.quantum = {"enabled": True}
_nw_on = evaluate(_nwq_on, [_nwq_on.thermal.T_hs])
_q = _nw_on["scalars"]
check("C2: nanowire flag off (default) adds no quantum_* key",
      not any(k.startswith("quantum_") for k in _nw_off["scalars"]))
check("C2: nanowire flag on changes NO existing key or value (flag-on result minus quantum_* "
      "keys == flag-off result, NaN-safe)", _canon_q(_strip_q(_nw_on)) == _canon_q(_nw_off))
check("C2 AC3: nanowire horizontal pulse row (valid, rectangular): filtered CW g2(0) in [0,1] "
      f"and above the memoryless rate-equation value ({_q['quantum_g2_cw0_filtered_dot']:.4f} "
      f"vs {_q['quantum_g2_cw0_rate_dot']:.4f}); [DR] provenance entry present",
      _q["valid"] is True and _q["quantum_invalid_reasons"] == []
      and 0.0 <= _q["quantum_g2_cw0_filtered_dot"] <= 1.0
      and _q["quantum_g2_cw0_filtered_dot"] > _q["quantum_g2_cw0_rate_dot"]
      and "[DR]" in _q["provenance"]["quantum"]
      and _q["quantum_filter_fwhm_meV"] == _nwq.dot.gamma300)
_Gq = _q["gamma_X_ns"] + _q["k_X_ns"]
_Iq = _q["quantum_hom_indistinguishability"]
check(f"C2 AC3: nanowire HOM I = {_Iq:.4e} in [0,1], Gamma = gamma_X_ns + k_X_ns of the row, "
      "I == Gamma/(Gamma + 2 gamma*) for the reported numbers and == hbar Gamma/FWHM (rel 1e-9)",
      0.0 <= _Iq <= 1.0 and close(_q["quantum_hom_gamma_ns"], _Gq, rel=1e-12)
      and close(_Iq, _Gq / (_Gq + 2.0 * _q["quantum_hom_gamma_star_ns"]), rel=1e-9, abs_=0.0)
      and close(_Iq, _HBAR_V * _Gq / _nwq.dot.gamma300, rel=1e-9, abs_=0.0))
# Two-level reduction on THIS row's own rates: w = L = 35 meV, D = r + gamma_X + k_X.
_r_ns = _q["r_captured_s"] * 1e-9
_w_ns = _nwq.dot.gamma300 / _HBAR_V
for _wfac in (1.0, 0.25):
    _q2 = _devmod._quantum_diagnostics(
        r_ns=_r_ns, gamma_X_ns=_q["gamma_X_ns"], gamma_XX_ns=_q["gamma_XX_ns"],
        k_X=_q["k_X_ns"], k_XX=_q["k_XX_ns"], pump_ratio=1.0,
        fwhm_X_meV=_nwq.dot.gamma300, fwhm_XX_meV=_nwq.dot.gamma300, delta_xx_meV=0.0,
        dx_meV=0.0, w_meV=_wfac * _nwq.dot.gamma300, eps_rate=float("nan"), levels=2)
    _D = _r_ns + _q["gamma_X_ns"] + _q["k_X_ns"]
    _W, _L = _wfac * _w_ns, _w_ns
    _ref = 2.0 * _D * (_W + _L) / ((_D + _W) * (3.0 * _W + _L))
    check(f"C2: two-level reduction of the wired filter path on the nanowire row's rates "
          f"(w = {_wfac:g} L): g2_f(0) = {_q2['quantum_g2_cw0_filtered_dot']:.6e} == "
          f"2D(w+L)/((D+w)(3w+L)) = {_ref:.6e} (rel 1e-6)",
          close(_q2["quantum_g2_cw0_filtered_dot"], _ref, rel=1e-6, abs_=0.0))
_nws = _DD.load(str(ROOT / "cards" / "nitride-nanowire-vertical-set-design.yaml"))
_nws.quantum = {"enabled": True}
_qs_set = evaluate(_nws, [_nws.thermal.T_hs])["scalars"]
check("C2: deterministic_pair (SET) row: no continuous pump -> CW filtered g2 NaN with an explicit "
      "reason; HOM still reported and equal to hbar Gamma/FWHM",
      math.isnan(_qs_set["quantum_g2_cw0_filtered_dot"])
      and any("deterministic_pair" in r for r in _qs_set["quantum_invalid_reasons"])
      and close(_qs_set["quantum_hom_indistinguishability"],
                _HBAR_V * (_qs_set["gamma_X_ns"] + _qs_set["k_X_ns"]) / _nws.dot.gamma300, rel=1e-9, abs_=0.0))
_nwq_bad = copy.deepcopy(_nwq); _nwq_bad.quantum = {"enabled": True, "x": 1}
check("C2: malformed quantum block rejected on the nanowire platform",
      raises(lambda: evaluate(_nwq_bad, [_nwq.thermal.T_hs])))

# --------------------------------------------------- C0 review lows (device.py, 2026-09-23)
# (ii) _nitride_flat_background_acceptance: a NaN kappa must propagate (it used
# to be silently mapped to 0.0 by `if not g > 0.0`); zero/negative kappa keeps
# the documented zero acceptance.
_acc = _devmod._nitride_flat_background_acceptance
check("C0 low (ii): NaN kappa -> NaN flat-background acceptance (not a silent 0.0); "
      "kappa 0 and kappa < 0 -> 0.0",
      math.isnan(_acc(float("nan"), 35.0, 0.0, 0.0, 0.0))
      and _acc(0.0, 35.0, 0.0, 0.0, 0.0) == 0.0 and _acc(-1.0, 35.0, 0.0, 0.0, 0.0) == 0.0
      and close(_acc(35.0, 35.0, 0.0, 0.0, 0.0), math.pi / 4.0, rel=1e-12))
# (i) a field_collapse (Zener) dot is bound but invalid: its planar-path reason
# is prefixed "invalid dot:", never "unbound dot:". The levels solver is
# wrapped for the operating-point call only (the first call, the cavity-
# tracking reference, passes through) so the check does not depend on which
# geometry currently collapses.
import dataclasses as _dc
_orig_levels = _devmod.nitride_levels.levels
_calls = {"n": 0}


def _collapsing_levels(system, T):
    lv = _orig_levels(system, T)
    _calls["n"] += 1
    if _calls["n"] == 1:
        return lv
    return _dc.replace(lv, valid=False, invalid_reasons=(
        "field_collapse (|F|*h_eff >= strained InGaN gap: interband Zener breakdown, "
        "unscreened field not self-consistent)",))


_pc = _DD.load(str(ROOT / "cards" / "nitride-cavity-pulse-design.yaml"))
_devmod.nitride_levels.levels = _collapsing_levels
try:
    _pcs = evaluate(_pc, [_pc.thermal.T_hs])["scalars"]
finally:
    _devmod.nitride_levels.levels = _orig_levels
check("C0 low (i): a field_collapse row reads 'invalid dot: field_collapse ...', not "
      f"'unbound dot: ...' (got {_pcs['invalid_reasons']})",
      _pcs["valid"] is False and len(_pcs["invalid_reasons"]) == 1
      and _pcs["invalid_reasons"][0].startswith("invalid dot: field_collapse")
      and not any(r.startswith("unbound dot") for r in _pcs["invalid_reasons"]))

# --------------------------------------------------------------- audit C6: explicit loading window
# Audit C6 (2026-09-23, M3; user decision Q6: the loading window is
# arbitrary). drive.diode.loading_window_ns is a card leaf defaulting to
# tau_pulse_ns; evaluate_nanowire's loading_window_ns keyword overrides it.
import copy as _copy_c6  # noqa: E402
_C6_CARDS = ("nitride-nanowire-horizontal-pulse-design.yaml", "nitride-nanowire-horizontal-set-design.yaml",
             "nitride-nanowire-vertical-pulse-design.yaml", "nitride-nanowire-vertical-set-design.yaml")
for _cn in _C6_CARDS:
    _d_c6 = DeviceDesign.load(str(ROOT / "cards" / _cn))
    _tau_c6 = float(_d_c6.drive.diode.get("tau_pulse_ns", 0.1))
    _a_c6 = dev.evaluate_nanowire(_copy_c6_copy := _copy_c6.deepcopy(_d_c6))
    _b_c6 = dev.evaluate_nanowire(_copy_c6.deepcopy(_d_c6), loading_window_ns=_tau_c6)
    check(f"C6 {_cn}: no card sets loading_window_ns, and the default equals an explicit "
          f"loading_window_ns=tau_pulse_ns ({_tau_c6} ns) bit for bit (scalars and curves)",
          "loading_window_ns" not in _d_c6.drive.diode
          and "loading_window_ns" not in _copy_c6_copy.drive.diode
          and nan_eq(_a_c6["scalars"], _b_c6["scalars"])
          and all(nan_eq(_a_c6["curves"][k].tolist(), _b_c6["curves"][k].tolist()) for k in _a_c6["curves"]))

_d_c6p = DeviceDesign.load(str(ROOT / "cards" / "nitride-nanowire-horizontal-pulse-design.yaml"))
_period_c6 = 1e9 / float(_d_c6p.drive.rep_rate_hz)
_c6_rect = {t: dev.evaluate_nanowire(_copy_c6.deepcopy(_d_c6p), loading_window_ns=t)["scalars"]
            for t in (0.1, 1.0, 3.0) if t < _period_c6}  # 200 MHz card: period 5 ns
_d_c6leaf = _copy_c6.deepcopy(_d_c6p); _d_c6leaf.drive.diode["loading_window_ns"] = 1.0
_c6_leaf = dev.evaluate_nanowire(_d_c6leaf)["scalars"]
check("C6 card leaf drive.diode.loading_window_ns=1 and the keyword loading_window_ns=1 give identical rows",
      1.0 in _c6_rect and nan_eq(_c6_leaf, _c6_rect[1.0]))


def _c6_missed_closed_form(s, t_ns):
    # Poisson loading of each carrier over the window: P(miss) =
    # 1 - (1 - e^{-r_e t})(1 - e^{-r_h t}) (independent carriers).
    return 1.0 - (1.0 - math.exp(-s["rti_e_rate_Hz"] * t_ns * 1e-9)) * (1.0 - math.exp(-s["rti_h_rate_Hz"] * t_ns * 1e-9))


check("C6 rectangular: the second-carrier gate IS the loading window (rti_gate_ns == loading_window_ns) "
      "and the second-carrier probability grows strictly with it (1 - exp(-r_2 t))",
      len(_c6_rect) == 3 and all(s["rti_gate_ns"] == t for t, s in _c6_rect.items())
      and _c6_rect[0.1]["rti_second_pair_probability"] < _c6_rect[1.0]["rti_second_pair_probability"]
      < _c6_rect[3.0]["rti_second_pair_probability"])
check("C6 rectangular: missed load equals the Poisson closed form over the loading window and the "
      "pair_missed_load_window failure is window-independent (hole supply ~kHz-MHz, 0.1-3 ns)",
      all(abs(s["rti_missed_load_probability"] - _c6_missed_closed_form(s, t)) < 1e-12
          and "pair_missed_load_window" in s["rti_failed_checks"] for t, s in _c6_rect.items()))
check("C6 rectangular: the optical/pulse-counting path does not depend on the loading window "
      "(g2_op, flux, gate_ns_used identical at 0.1/1/3 ns); only rti_* moves",
      len({(s["g2_op"], s["collected_flux_pulsed_s"], s["gate_ns_used"]) for s in _c6_rect.values()}) == 1)

_d_c6s = DeviceDesign.load(str(ROOT / "cards" / "nitride-nanowire-horizontal-set-design.yaml"))
_period_c6s = 1e9 / float(_d_c6s.drive.rep_rate_hz)
_c6_set = {t: dev.evaluate_nanowire(_copy_c6.deepcopy(_d_c6s), loading_window_ns=t)["scalars"] for t in (0.1, 1.0)}
check("C6 deterministic_pair: the counting gate stays the card gate (period default), so the "
      "second-carrier probability is window-independent while the missed load follows the window",
      all(s["rti_gate_ns"] == _c6_set[0.1]["gate_ns_used"] for s in _c6_set.values())
      and _c6_set[0.1]["rti_second_pair_probability"] == _c6_set[1.0]["rti_second_pair_probability"]
      and all(abs(s["rti_missed_load_probability"] - _c6_missed_closed_form(s, t)) < 1e-12 for t, s in _c6_set.items()))
_c6_bad = dev.evaluate_nanowire(_copy_c6.deepcopy(_d_c6p), loading_window_ns=_period_c6)["scalars"]
_c6_bad0 = dev.evaluate_nanowire(_copy_c6.deepcopy(_d_c6p), loading_window_ns=0.0)["scalars"]
check("C6 a loading window >= the period (or <= 0) is an explicit invalid row (contract bullet 7: strictly < 1/rep_rate)",
      _c6_bad["valid"] is False and _c6_bad["invalid_reasons"] == ["invalid loading window"]
      and _c6_bad0["valid"] is False and _c6_bad0["invalid_reasons"] == ["invalid loading window"])

print(f"{sum(checks)}/{len(checks)} nitride nanowire device checks passed")
assert all(checks), "one or more nitride nanowire device checks failed (see FAIL lines above)"
