"""Nanowire device integration (piece 7).  This module deliberately imports
no device module at module scope: :mod:`fsim_core.device` dispatches here
locally, avoiding a cycle; the reverse direction (this module reusing
device.py's ``_nitride_reservoir_energy_eV`` helper, and the row-key
fixture's own DeviceDesign construction below) is also a local, function-
scope import for the same reason.

The 1.0 ns bare lifetime and the 0.1 ns electrical pulse are card-level [A]
design inputs, not fits to Deshpande et al. (APL 2014).  Every other
material/geometry/rate input is supplied by the card or one of pieces 2-6
and retains ITS OWN provenance tag; nothing here invents a class default --
see ``_REQUIRED_*`` below for the leaves this module refuses to silently
default (fix-2 finding L18).

Conventions reused verbatim from the planar evaluator (``fsim_core/device.py``
``_evaluate_nitride``, contract "Composition rules for the device piece"):
  - g2_op = 1 - rho**2 * (1 - g2_dot) (device.py:1278 convention), bounded
    in [0, 1] by construction whenever rho is in [0, 1] and g2_dot is.
  - rho = collected signal / (collected signal + collected background),
    where "collected" means the collection factor (eta_collection_X/XX) is
    applied EXACTLY ONCE -- inside the pulse_counting propagation (the t_X/
    t_XX argument), never again afterward -- and the filter transmission
    eta_out is then applied ONCE, after propagation, exactly like the
    planar path's separate cr["eta_out"] * cnt["mean_counts"] step (fix-2
    finding H3: an earlier round applied ph["eta_collection_X"] both
    inside AND outside pulse_g2/deterministic_cycle_g2).
  - Scalars are the row nearest thermal.T_hs in a supplied T_grid (never
    rows[-1]).
  - eta_out (fix-2 round 2, MEDIUM 8): the planar nitride path's own
    "eta_out" (device.py's cr["eta_out"]) is NitrideCavityParams.eta_out, a
    FIXED [A] cavity out-coupling constant -- unrelated to any dx/w
    Lorentzian window -- and nanowire cards forbid a cavity entirely
    (cavity.enabled must be False).  This platform's photonics collection
    (ph["eta_collection_X/XX"]) is explicitly geometric/waveguide-only
    ("no cavity Q, DBR, or spectral linewidth filter" --
    nitride_nanowire_photonics.py's own module docstring), so the genuine
    SPECTRAL FilterBlock transmission is applied here instead, exactly as
    the LEGACY (non-nitride) evaluate() path applies it: spectral.epsilon
    with kappa=None (no cavity mode to net against), gated on
    d.filter.enabled (device.py:1643) -- never an invented closed form.
"""
from __future__ import annotations

import math
import re
import warnings
from dataclasses import replace
from pathlib import Path

import numpy as np

from . import pulse_counting
from .drive_mech import set_feasibility
from .loading import E_SI, KB_SI, island_radius_nm
from .nitride_nanowire_injector import NitrideNanowireInjectorParams, injector_feasibility
from .nitride_nanowire_levels import NitrideNanowireSystem, levels, rates
from .nitride_nanowire_surface import NitrideNanowireSurfaceParams, surface_rates
from .nitride_nanowire_transport import (
    NitrideWireDiode, wire_operating_point, evaluate_injection, pulse_delivery,
)
from .nitride_nanowire_photonics import NitrideNanowirePhotonicsParams, response
from .spectral import epsilon

_CONTRACT_PATH = Path(__file__).resolve().parents[1] / "docs" / "nitride_nanowire_contract.md"

# nitride.dot / nitride.nanowire leaf allow-lists (fix-2 finding M17 "no
# QW"/"no silently ignored options"): the surface/photonics/injector/diode
# blocks are already fed straight into their own frozen dataclasses below
# (``Cls(**mapping)``), so an unrecognized key there already raises a loud
# TypeError from the dataclass constructor itself; these two blocks are
# hand-parsed by this module instead, so nothing would otherwise catch a
# stray planar-only leaf (e.g. a QW fluctuation's wl_thickness_nm/
# geometry_type) silently doing nothing.
_ALLOWED_DOT_KEYS = frozenset((
    "radius_nm", "height_nm", "x_in", "strain_fraction", "screening_fraction",
    "external_field_kVcm", "vbo_InN_GaN_eV", "strain_c_fraction",
    "k_intrinsic_ns",  # [A] intrinsic occupied-dot loss beyond surface/escape; not a contract card leaf, defaults 0.0 (neutral)
))
_ALLOWED_NANOWIRE_KEYS = frozenset((
    "family", "core_radius_nm", "outer_radius_nm", "strain_bound", "shell",
    "barrier_left_nm", "barrier_right_nm",
))
# MEDIUM 4 (fix-3): the top-level nitride dict's own leaves -- the contract's
# sub-blocks (nanowire/dot/surface/photonics/wire_thermal/injector), the two
# bare scalar leaves it requires directly (tau_rad0_ns, tau_cap_ps), and the
# single optional scalar override this platform DOES honor
# (reservoir_energy_eV, via the planar _nitride_reservoir_energy_eV helper's
# own "background" dict convention -- see _one_raw). No QW/cavity/wetting-
# layer leaf (nitride.qw, .wl, .cavity, ...) has any meaning for a nanowire
# card -- previously silently accepted and ignored (fix-3 finding MEDIUM 4).
_ALLOWED_TOP_LEVEL_NITRIDE_KEYS = frozenset((
    "tau_rad0_ns", "tau_cap_ps", "nanowire", "dot", "surface", "photonics",
    "wire_thermal", "injector", "reservoir_energy_eV",
))


def _nan(): return float("nan")


def _lifetime_ns(rate_ns):
    """1/rate, distinguishing a genuinely zero rate (infinite lifetime) from
    a non-finite/not-computed one (NaN), matching device.py's
    _nitride_lifetime_ns convention (interface constraints: "infinity only
    for a genuinely zero decay rate")."""
    if not isinstance(rate_ns, (int, float)) or not math.isfinite(rate_ns):
        return _nan()
    return float("inf") if rate_ns <= 0 else 1.0 / rate_ns


# --------------------------------------------------------------- row-key set (HIGH 1, fix-2 round 2)
#
# The scalar schema is deliberately NOT a hand-maintained literal list: it is
# the UNION of (a) every column parsed independently from the contract's own
# "Row columns" section (docs/nitride_nanowire_contract.md) and (b) every
# key a REAL, successful evaluation of this module's own physics chain
# actually returns (harvested once, from an in-memory fixture card, and
# cached).  This guarantees invalid rows carry EXACTLY the valid-row key set
# by construction -- including any key a dependency module (most notably the
# injector's rti_* diagnostics) adds later, with no second list to keep in
# sync (fix-2 round 1 finding: a hardcoded 145-entry list omitted
# rti_polarity; the generalized failure mode is "any future rti_* key is
# missed on invalid rows").


# MEDIUM 3 (fix-3): a hand-maintained per-citation denylist (the previous
# round's "evaluate_nanowire", "e696bdd", "rti_", ...) is fragile -- any
# NEW backtick-quoted citation the contract prose ever grows (a function
# name, a commit hash, a file path) that is not yet in the list silently
# becomes a bogus extra "column" (a permanently-NaN key on every row,
# undetected).  Replaced with a STRUCTURAL rule instead: a genuine column
# name is only ever backtick-quoted at PAREN DEPTH 0 in the section's own
# text (every citation/aside actually in this doc today sits inside a
# parenthetical, including ones that hide their own closing colon inside
# the parenthetical -- see _parse_row_columns_section); the residual
# denylist below is not citations but a small, closed set of Python
# value-literals ("True"/"False"/"None"/"not_applicable") that legitimately
# appear backtick-quoted in top-level prose as EXAMPLE VALUES, plus the two
# function names the section's own un-parenthesized opening sentence names
# before any column list starts. This reproduces today's 149-column list
# exactly (verified against the old denylist-based parser) without needing
# to grow with every new citation.
_STABLE_VALUE_LITERAL_DENY = frozenset((
    "True", "False", "None", "not_applicable",
    "evaluate_nanowire", "_evaluate_nitride",
))

# Frozen fallback (MEDIUM 3): captured from a real, successful parse of
# docs/nitride_nanowire_contract.md's "Row columns" section (2026-09-14,
# via the depth-based parser below) -- used ONLY when the contract file
# itself cannot be read at import time (moved/deleted), so a missing doc
# degrades to "possibly stale row-key set" instead of bricking every
# evaluation. Regenerate by running the parser against the live doc and
# pasting its output back in here if the contract's Row-columns section is
# intentionally changed.
_FROZEN_CONTRACT_COLUMNS = (
    "platform", "family", "T_hs", "T_j", "thermal_iterations", "thermal_converged",
    "cycle_loading", "rep_rate_hz", "strain_bound", "bound_role", "bound_reversal",
    "screening_fraction", "core_radius_nm", "outer_radius_nm", "conducting_radius_nm",
    "E_X_eV", "lambda_nm", "field_kVcm", "overlap_sq", "electron_bound", "hole_bound",
    "E_a_meV", "gamma_X0_ns", "gamma_XX0_ns", "k_X_ns", "k_XX_ns", "escape_prefactor_ns",
    "tau_cap_ps_used", "reservoir_state_count_e", "reservoir_state_count_h",
    "sidewall_overlap", "provenance", "radius_over_lambda", "V_number", "beta_HE11",
    "eta_collection_X", "eta_collection_XX", "radiative_rate_factor", "gamma_X_ns",
    "gamma_XX_ns", "approximation_error", "single_mode", "beta_multimode_penalty",
    "degree_of_linear_polarization", "antenna_rate_factor", "k_side_ns",
    "k_surface_reservoir_ns", "k_surface_X_ns", "k_surface_XX_ns",
    "shell_multiplier_used", "reservoir_access_used", "occupied_dot_access_used",
    "area_cm2", "J_A_cm2", "V_j", "V_terminal", "depletion_field_kVcm", "C_dep_F",
    "V_bi_minus_Eg_mV", "flat_band", "depletion_regime", "eta_inj", "f_capture",
    "f_qfl_dot", "f_qfl_dot_thermodynamic_limit", "f_qfl_background", "r_supply_s",
    "r_captured_s", "r_matrix_radiative_s", "r_matrix_nonradiative_s",
    "r_surface_reservoir_s", "r_leakage_s", "mu", "power_on_W", "accounting_residual_s",
    "tau_RC_ns", "delivered_step_fraction", "pulse_delivery_feasible", "tau_rad_bare_ns",
    "tau_rad_photonic_ns", "tau_total_X_ns", "collected_flux_pulsed_s",
    "collected_flux_x_s", "collected_flux_xx_s", "background_flux_s",
    "total_detected_flux_s", "mean_counts", "mean_counts_x", "mean_counts_xx",
    "gate_ns_used", "S_X", "S_XX", "rho_pulsed", "blocked_load_probability",
    "counting_converged", "eta_out", "g2_op", "one_pair_valid", "pair_supply_possible",
    "valid", "invalid_reasons", "collected_flux_delivered_s", "set_feasible",
    "set_E_C_meV", "set_EC_over_kT", "set_radius_nm", "set_radius_max_nm",
    "set_C_sigma_F", "set_R_T_over_RQ", "set_f_max_Hz", "rti_feasible", "rti_status",
    "rti_transport_feasible", "rti_level_margin_kT", "rti_orbital_margin_kT",
    "rti_alignment_error_meV", "rti_linewidth_meV", "rti_rate_Hz", "rti_e_rate_Hz",
    "rti_h_rate_Hz", "rti_bypass_fraction", "rti_missed_load_probability",
    "rti_second_pair_probability", "rti_growth_feasible", "rti_failed_checks",
    "rti_evidence_status", "rti_level_margin_e_kT", "rti_level_margin_h_kT",
    "rti_alignment_error_e_meV", "rti_alignment_error_h_meV",
    "rti_required_bias_shift_meV", "rti_bypass_fraction_tsai_partition",
    "rti_second_carrier_probability", "rti_growth_nearest_commensurate_nm",
    "rti_growth_perturbed_margins_kT", "rti_barrier_polarization_tilt_eV",
    "rti_well_to_dot_drop_meV", "rti_p_free_cm3", "rti_numerics_ok", "rti_gate_ns",
    "rti_polarity", "rti_reservoir_state_count_e", "rti_reservoir_state_count_h",
    "optical_pass", "hardware_qualified", "rti_qualified", "device_pass",
    "rti_device_pass", "headline_eligible",
)


def _parse_row_columns_section(section):
    """Depth-based extraction (MEDIUM 3): a backtick span is a genuine
    column name only when its OPENING backtick sits at parenthesis depth 0
    in ``section`` (parens are tracked char-by-char across the whole
    multi-line section, so a citation/aside spanning several lines --
    including one that hides its own closing colon inside the
    parenthetical, as the "RT injector screen, fix-1..." paragraph does --
    is excluded regardless of nesting). A depth-0 span containing "="
    still uses the established formula-LHS convention (e.g.
    "collected_flux_delivered_s = ...")."""
    cols = []
    depth = 0
    i, n = 0, len(section)
    while i < n:
        c = section[i]
        if c == "(":
            depth += 1; i += 1; continue
        if c == ")":
            depth = max(0, depth - 1); i += 1; continue
        if c == "`":
            j = section.find("`", i + 1)
            if j == -1:
                break
            span = section[i + 1:j]
            if depth == 0:
                if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", span):
                    if span not in _STABLE_VALUE_LITERAL_DENY and span not in cols:
                        cols.append(span)
                elif "=" in span:
                    lhs = span.split("=")[0].strip()
                    if re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", lhs) and lhs not in cols:
                        cols.append(lhs)
            i = j + 1
            continue
        i += 1
    return cols


def _contract_row_columns():
    """Parse the contract's "Row columns" section for the identifiers it
    names.  Independent of verify/verify_nitride_nanowire_device.py's own
    parser (which exists to check this module's output from the outside,
    per that file's docstring) -- this copy exists to seed this module's OWN
    row-key set, not to be re-checked by it.

    MEDIUM 3 (fix-3): fails LOUDLY, naming the missing heading, when the
    doc is present but malformed (renamed/removed heading) -- silently
    parsing an empty or truncated section would silently shrink every
    row's key set. A genuinely MISSING file instead degrades to the frozen
    fallback list with a one-time warning (this function runs exactly once,
    at import, to populate _CONTRACT_ROW_COLUMNS below, so "once" is
    automatic here -- not a manual dedup flag)."""
    try:
        text = _CONTRACT_PATH.read_text(encoding="utf-8")
    except FileNotFoundError:
        warnings.warn(
            "nitride_nanowire_device: contract doc not found at %s -- falling "
            "back to the frozen Row-columns list captured in this module "
            "(stale if the contract has changed since)." % _CONTRACT_PATH,
            RuntimeWarning, stacklevel=2)
        return list(_FROZEN_CONTRACT_COLUMNS)
    try:
        start = text.index("## Row columns")
    except ValueError:
        raise RuntimeError(
            "nitride_nanowire_device: %s is missing the '## Row columns' "
            "heading required to parse the row-column contract" % _CONTRACT_PATH)
    try:
        end = text.index("## Sweep grid", start)
    except ValueError:
        raise RuntimeError(
            "nitride_nanowire_device: %s is missing the '## Sweep grid' "
            "heading that closes the '## Row columns' section" % _CONTRACT_PATH)
    return _parse_row_columns_section(text[start:end])


# Frozen at import (MEDIUM 3): parsed once here, not re-parsed on every
# _row_keys() call.
_CONTRACT_ROW_COLUMNS = _contract_row_columns()


def _row_key_fixture_design():
    """A minimal, self-consistent, always-valid nanowire card (horizontal_
    as_built / deterministic_pair / relaxed, at the same operating point the
    review's independent numbers were confirmed against) used ONLY to
    harvest the full set of keys the physics chain actually returns on its
    success path -- it reaches every downstream module (levels/rates/
    surface/photonics/transport/injector) including the SET and RT-injector
    screens.  Local import of fsim_core.device: a module-scope import here
    would create the exact cycle the module docstring's first paragraph
    documents avoiding (device.py dispatches to this module locally)."""
    from .device import (
        DeviceDesign, DriveBlock, RetentionBlock, CavityBlock, EmissionBlock,
        ThermalBlock, DotBlock,
    )
    nanowire = {"family": "horizontal_as_built", "core_radius_nm": 12.5,
                "outer_radius_nm": 12.5, "strain_bound": "relaxed",
                "barrier_left_nm": 15.0, "barrier_right_nm": 15.0}
    dot = {"radius_nm": 12.5, "height_nm": 2.0, "x_in": 0.40, "strain_fraction": 0.0}
    diode = {"preset": "nitride-nanowire", "tau_pulse_ns": 0.1}
    set_params = {"R_T_ohm": 1e6}
    drive_kw = dict(mode="EL-transport", I_uA=0.02, duty=0.008, rep_rate_hz=80e6,
                     b_res=0.1, cycle_loading="deterministic_pair", diode=diode,
                     set_params=set_params)
    thermal_block = {"Rth_K_W": 3.1e9, "eta_total": 0.01, "f_Rs_local": 1.0,
                      "R_s_ohm": 2.38e9, "C_parasitic_F": 0.0}
    return DeviceDesign(
        platform="ingan_gan_nanowire",
        dot=DotBlock(linewidth="anchored", lineshape="lorentzian", gamma300=3.0),
        ret=RetentionBlock(mode="nitride_confinement"),
        drive=DriveBlock(**drive_kw), thermal=ThermalBlock(T_hs=300.0),
        cavity=CavityBlock(enabled=False), emission=EmissionBlock(type="nanowire"),
        nitride={"tau_rad0_ns": 1.0, "tau_cap_ps": 10.0, "nanowire": nanowire, "dot": dot,
                 "surface": {"occupied_dot_access": 0.05},
                 "photonics": {"NA": 0.5}, "wire_thermal": thermal_block, "injector": {}})


_ROW_KEYS_CACHE = None


def _row_keys():
    global _ROW_KEYS_CACHE
    if _ROW_KEYS_CACHE is None:
        try:
            fixture_row = _one_raw(_row_key_fixture_design(), 300.0)
        except Exception as exc:
            # MEDIUM 3: name the fixture, not just "something failed" -- a
            # broken _row_key_fixture_design() card raises OUTSIDE
            # _one_raw's own try/except (card-shape errors are not caught
            # there), so without this wrapper the traceback would point at
            # _validate()/_one_raw() with no hint that the FIXTURE itself
            # (not a real card) is what needs fixing.
            raise RuntimeError(
                "nitride_nanowire_device._row_key_fixture_design() produced a "
                "card _one_raw() could not evaluate at all (cannot derive the "
                "module's own row-key set from it): %s" % exc) from exc
        if not fixture_row.get("valid"):
            raise RuntimeError(
                "nitride_nanowire_device._row_key_fixture_design()'s fixture "
                "card evaluated but was physically invalid (cannot derive the "
                "module's own row-key set from it): " + "; ".join(fixture_row.get("invalid_reasons", [])))
        _ROW_KEYS_CACHE = frozenset(_CONTRACT_ROW_COLUMNS) | frozenset(fixture_row)
    return _ROW_KEYS_CACHE


def _invalid_raw(T, reasons, **extra):
    """The invalid-row fields this module itself imposes (never a dependency
    module's output) -- callers add coordinate/regime context via **extra.
    Prefilling the FULL row-key set with NaN is the caller's job (_one), not
    this function's, so this same helper can also serve the raw success
    path's error branches without depending on _row_keys() (which itself
    depends on a SUCCESSFUL call here -- see _row_keys()'s fixture)."""
    row = {"platform": "ingan_gan_nanowire", "T_hs": float(T), "T_j": _nan(),
           "valid": False, "invalid_reasons": list(reasons), "g2_op": _nan(),
           "collected_flux_pulsed_s": _nan(), "collected_flux_delivered_s": _nan(),
           "optical_pass": False, "hardware_qualified": False,
           "rti_qualified": False, "device_pass": False, "rti_device_pass": False,
           "headline_eligible": False,
           "provenance": {}, "set_feasible": "not_applicable", "rti_feasible": "not_applicable",
           "rti_status": "not_applicable", "rti_failed_checks": [], "bound_reversal": "not_computed"}
    row.update(extra)
    return row


def _params(mapping, cls):
    if not isinstance(mapping, dict): raise ValueError("nanowire card blocks must be mappings")
    return cls(**mapping)


def _req(mapping, key, where):
    """L18: a card-schema leaf this module refuses to silently default --
    reject a missing value with a NAMED ValueError instead of an untagged
    hardcoded fallback (fix-1 finding L18: x_in 0.4, height 2, barriers 15,
    tau_rad0 1, tau_cap 10, R_T 1e6 were all previously silent defaults)."""
    if not isinstance(mapping, dict) or key not in mapping:
        raise ValueError(f"{where}.{key} is required for nanowire cards")
    return float(mapping[key])


def _validate(d):
    n = d.nitride
    # MEDIUM 4 (fix-3): top-level allow-list -- reject nitride.qw / .wl /
    # .cavity / any other planar-only or unknown leaf by name instead of
    # silently ignoring it (a nanowire card has no QW, no wetting layer, no
    # cavity block to configure at all).
    unknown_top = [k for k in n if k not in _ALLOWED_TOP_LEVEL_NITRIDE_KEYS]
    if unknown_top: raise ValueError("nitride has unrecognized (e.g. QW/wetting-layer/cavity-only) top-level leaves for a nanowire card: %s" % unknown_top)
    for k in ("nanowire", "surface", "photonics", "wire_thermal", "injector", "dot"):
        if k not in n: raise ValueError("nitride.%s is required for nanowire cards" % k)
    nw, dot = n["nanowire"], n["dot"]
    unknown_nw = [k for k in nw if k not in _ALLOWED_NANOWIRE_KEYS]
    if unknown_nw: raise ValueError("nitride.nanowire has unrecognized leaves: %s" % unknown_nw)
    unknown_dot = [k for k in dot if k not in _ALLOWED_DOT_KEYS]
    if unknown_dot: raise ValueError("nitride.dot has unrecognized (e.g. QW-only) leaves for a nanowire card: %s" % unknown_dot)
    family = nw.get("family")
    if family not in ("horizontal_as_built", "vertical_photonic"): raise ValueError("unknown nanowire family")
    if d.drive.mode != "EL-transport" or d.drive.diode.get("preset") != "nitride-nanowire": raise ValueError("nanowire requires EL-transport nitride-nanowire diode")
    if d.dot.linewidth != "anchored" or d.dot.lineshape != "lorentzian" or d.ret.mode != "nitride_confinement" or d.cavity.enabled or d.emission.type != "nanowire": raise ValueError("nanowire requires anchored Lorentzian, nitride_confinement, cavity disabled, emission nanowire")
    if d.ret.tau_cap_scales_with_density: raise ValueError("nanowire forbids density-scaled capture")
    if d.drive.loading_model != "auto": raise ValueError("nanowire uses pulse-counting, not planar loading model")
    if d.drive.cycle_loading not in ("rectangular", "deterministic_pair"): raise ValueError("drive.cycle_loading must be rectangular or deterministic_pair")
    if d.aperture != type(d.aperture)(): raise ValueError("nanowire forbids non-default planar drive.aperture options (no Q/DBR/aperture-supply axis)")
    # LOW 14 (fix-2 round 2): eta_load was silently hardcoded to 1.0 downstream
    # regardless of what the card said, so a card that set drive.eta_load to
    # anything else was accepted and silently ignored; reject it by name
    # instead -- this platform models SET/pulse loading explicitly (drive_mech,
    # pulse_counting) and has no separate legacy "loading efficiency" axis.
    if d.drive.eta_load != 1.0: raise ValueError("nanowire cards do not support drive.eta_load != 1.0 (SET/pulse loading is modeled explicitly, not via a scalar loading efficiency)")
    if nw.get("strain_bound") not in ("relaxed", "unrelaxed"): raise ValueError("nanowire strain_bound required")
    sf = dot.get("strain_fraction")
    expected = 0.0 if nw["strain_bound"] == "relaxed" else 1.0
    if sf is not None and float(sf) != expected: raise ValueError("nitride.dot.strain_fraction contradicts strain_bound")
    core = float(nw["core_radius_nm"]); disc = float(dot.get("radius_nm", core))
    if family == "horizontal_as_built" and disc != core: raise ValueError("horizontal disc radius must equal core radius")
    if family == "vertical_photonic" and disc >= core: raise ValueError("vertical photonic disc radius must be less than core radius")
    outer = float(nw.get("outer_radius_nm", core))
    if outer < core: raise ValueError("outer radius below core")
    # HIGH 3 (fix-2 round 2): shell was allow-listed but never read or
    # cross-checked against the contract's own outer-radius rule ("shell
    # 'none' requires outer_radius_nm == core_radius_nm; shell 'AlGaN' uses a
    # declared 3 nm shell thickness").
    shell = nw.get("shell", "none")
    if shell not in ("none", "AlGaN"): raise ValueError("nitride.nanowire.shell must be 'none' or 'AlGaN'")
    if shell == "none" and outer != core:
        raise ValueError("nitride.nanowire.shell='none' requires outer_radius_nm == core_radius_nm")
    if shell == "AlGaN" and outer < core + 3.0:
        raise ValueError("nitride.nanowire.shell='AlGaN' requires outer_radius_nm >= core_radius_nm + 3 (contract's declared 3 nm shell thickness)")
    setp = d.drive.set_params
    if setp.get("radius_nm", disc) != disc or "C_sigma_F" in setp: raise ValueError("SET radius must be the disc radius; C_sigma override forbidden")
    conducting = float(d.drive.diode.get("conducting_radius_nm", core))
    if not (0.0 < conducting <= core): raise ValueError("drive.diode.conducting_radius_nm must be in (0, core_radius_nm]")
    bound_role = "headline_upper" if nw["strain_bound"] == "relaxed" else "conservative_lower"
    ext_field_kVcm = float(dot.get("external_field_kVcm", 0.0))
    return (family, core, disc, outer, conducting, d.drive.cycle_loading, float(d.drive.rep_rate_hz),
            nw["strain_bound"], bound_role, float(dot.get("screening_fraction", 0.0)), shell, ext_field_kVcm)


def _one_raw(d, T_hs):
    """The physics chain, returning EITHER the raw success-path merge
    (exactly and only the keys the six dependency modules plus this
    module's own explicit fields produce -- no pre-fill) or an
    _invalid_raw(...) dict on any card-validated-but-physically-invalid
    branch.  _one() (below) is the public per-row entry point: it adds the
    NaN pre-fill for the full row-key set on top of this function's result."""
    family, core, disc, outer, conducting, regime, rep, strain_bound, bound_role, screening, shell, ext_field_kVcm = _validate(d)
    n = d.nitride
    thermal = n["wire_thermal"]

    def inv(reasons, **extra):
        # LOW 20 (fix-2 round 2): not_applicable is reserved for pulse
        # (rectangular) rows; an invalid deterministic_pair (SET) row must
        # report the hardware screens as False (unresolved, not "N/A"),
        # exactly like a valid deterministic_pair row would if the screen
        # itself failed -- never the pulse-regime sentinel.
        hw_default = False if regime == "deterministic_pair" else "not_applicable"
        return _invalid_raw(T_hs, reasons, family=family, cycle_loading=regime, rep_rate_hz=rep,
                             strain_bound=strain_bound, bound_role=bound_role, screening_fraction=screening,
                             core_radius_nm=core, outer_radius_nm=outer, conducting_radius_nm=conducting,
                             set_feasible=hw_default, rti_feasible=hw_default, **extra)

    # Card-shape / required-leaf validation and dataclass construction from
    # raw card leaves happen OUTSIDE the numerical try/except below: a
    # missing leaf or a malformed sub-block value must be REJECTED (raise),
    # never silently downgraded to an invalid row (acceptance criterion 1).
    # Only the PHYSICS chain (thermal solve, levels, rates, injection,
    # counting, ...) -- operating on already-validated inputs -- is allowed
    # to fail into an explicit invalid row.
    rep_v = rep; pulse = float(d.drive.diode.get("tau_pulse_ns", .1))  # [A] 0.1 ns headline electrical pulse (module docstring); a card leaf, not a fitted physics constant
    duty = float(d.drive.duty)
    if rep_v <= 0 or pulse <= 0:
        return inv(["invalid pulse/rep-rate"])
    period = 1e9 / rep_v
    # LOW 15 (fix-2 round 2): the loading window must be STRICTLY shorter
    # than the period (pulse == period leaves no off-time at all), stricter
    # than the planar path's `>` (which admits equality) by design choice.
    if pulse >= period:
        return inv(["invalid pulse/rep-rate"])
    gate = float(d.drive.gate_ns if d.drive.gate_ns is not None else period)

    height_nm = _req(n["dot"], "height_nm", "nitride.dot")
    x_in = _req(n["dot"], "x_in", "nitride.dot")
    barrier_left = _req(n["nanowire"], "barrier_left_nm", "nitride.nanowire")
    barrier_right = _req(n["nanowire"], "barrier_right_nm", "nitride.nanowire")
    tau_rad0_ns = _req(n, "tau_rad0_ns", "nitride")
    tau_cap_ps = _req(n, "tau_cap_ps", "nitride")
    setp = d.drive.set_params
    R_T_ohm = _req(setp, "R_T_ohm", "drive.set_params")

    # Contact/thermal controls are card inputs, never silently inherited
    # from the Deshpande comparison preset [V/E/A; see wire_thermal
    # schema]; H4: R_s_ohm, f_Rs_local, eta_total and C_parasitic_F are
    # all read from nitride.wire_thermal and threaded through below.
    diode_raw = dict(d.drive.diode)
    for k in ("preset", "tau_pulse_ns", "conducting_radius_nm"): diode_raw.pop(k, None)
    # LOW 13 (fix-2 round 2): drive.diode.R_s_ohm silently overrode
    # nitride.wire_thermal.R_s_ohm with no contradiction check; reject a
    # card that states both, disagreeing, instead of silently picking one.
    if "R_s_ohm" in diode_raw and float(diode_raw["R_s_ohm"]) != float(thermal["R_s_ohm"]):
        raise ValueError("drive.diode.R_s_ohm contradicts nitride.wire_thermal.R_s_ohm")
    diode_raw.setdefault("R_s_ohm", float(thermal["R_s_ohm"]))
    diode_raw.setdefault("f_Rs_local", float(thermal["f_Rs_local"]))
    diode = NitrideWireDiode(core_radius_nm=core, conducting_radius_nm=conducting,
        barrier_left_nm=barrier_left, barrier_right_nm=barrier_right,
        d_active_nm=height_nm, x_in=x_in, **diode_raw)

    # HIGH 2 (fix-2 round 2): nitride.dot.external_field_kVcm is an
    # independent junction-field card leaf (e.g. an externally applied bias
    # beyond the resolved depletion field); it was allow-listed but
    # discarded (sys_kw hardcoded external_field_kVcm=0.0).  It now seeds
    # BOTH the pre-feedback system (sysA) and the post-feedback one (sysB,
    # below), which adds the resolved depletion field ON TOP of it.
    sys_kw = dict(height_nm=height_nm, core_radius_nm=core, outer_radius_nm=outer,
                  disc_radius_nm=(None if family == "horizontal_as_built" else disc),
                  x_in=x_in, strain_bound=strain_bound, screening_fraction=screening,
                  external_field_kVcm=ext_field_kVcm)
    for opt in ("vbo_InN_GaN_eV", "strain_c_fraction"):
        if opt in n["dot"]: sys_kw[opt] = float(n["dot"][opt])
    sysA = NitrideNanowireSystem(**sys_kw)

    eta_total = float(thermal["eta_total"])
    # HIGH 3 (fix-2 round 2): nitride.nanowire.shell is the card's single
    # source of truth for the shell material; propagate it into the surface
    # params instead of leaving NitrideNanowireSurfaceParams at its own
    # 'none' class default regardless of the card, and reject a card that
    # states nitride.surface.shell differently (a genuine contradiction,
    # never silently resolved either way).
    surface_kw = dict(n["surface"])
    if "shell" in surface_kw:
        if surface_kw["shell"] != shell:
            raise ValueError("nitride.surface.shell contradicts nitride.nanowire.shell")
    else:
        surface_kw["shell"] = shell
    # HIGH 1 (fix-3): nitride.surface.shell_multiplier must agree with the
    # declared shell instead of silently falling through to
    # NitrideNanowireSurfaceParams's own class default (0.1) regardless of
    # what the card says -- the previous round allow-listed shell_multiplier
    # but never read or cross-checked it, so an AlGaN card that (like every
    # committed card, which all declare shell='none'/multiplier=1.0) states
    # shell_multiplier=1.0 was silently accepted as an AlGaN shell with NO
    # passivation effect, k_X_ns bit-identical to shell='none' (latent: no
    # committed card is AlGaN today). shell='none' is physically unpassivated
    # (multiplier must be 1.0); shell='AlGaN' must state a multiplier
    # strictly below 1.0 (1.0 would mean "AlGaN shell, no passivation",
    # indistinguishable from shell='none' but mislabeled -- a card
    # contradiction, not a silent no-op), defaulting to the contract's own
    # 0.1 [A] (docs/nitride_nanowire_contract.md "Row columns", Surface
    # piece 4 / shell_multiplier default table) when the card omits it.
    if "shell_multiplier" in surface_kw:
        _mult = float(surface_kw["shell_multiplier"])
        if shell == "none" and _mult != 1.0:
            raise ValueError("nitride.surface.shell_multiplier must be 1.0 when shell='none'")
        if shell == "AlGaN" and _mult >= 1.0:
            raise ValueError("nitride.surface.shell_multiplier must be < 1.0 when shell='AlGaN' (1.0 has no passivation effect, indistinguishable from shell='none')")
    else:
        surface_kw["shell_multiplier"] = 1.0 if shell == "none" else 0.1  # [A] contract default
    surface_params = _params(surface_kw, NitrideNanowireSurfaceParams)
    # Coordinator addition (cards review, fix-3): NitrideNanowirePhotonicsParams
    # requires dipole_weights to be an actual Python tuple (its own
    # __post_init__ does `isinstance(w, tuple)`), but YAML/JSON cards can
    # only express a list -- so nitride.photonics.dipole_weights was
    # unreachable from a card (always raised "must be a 3-tuple", even for
    # a well-formed 3-element list). Coerce a 3-element list/tuple to a
    # tuple of floats here, at the card boundary; a wrong LENGTH is rejected
    # with a named error naming the leaf (an off-length tuple would
    # otherwise still hit the class's own generic message, since a length-
    # check only fires there after the isinstance(tuple) gate this coercion
    # satisfies). A non-sequence value is left untouched, so the class's own
    # "must be a 3-tuple" error still fires for it.
    photonics_kw = dict(n["photonics"])
    if "dipole_weights" in photonics_kw:
        _dw = photonics_kw["dipole_weights"]
        if isinstance(_dw, (list, tuple)):
            if len(_dw) != 3:
                raise ValueError("nitride.photonics.dipole_weights must have exactly 3 elements (along_wire, transverse_inplane, vertical)")
            photonics_kw["dipole_weights"] = tuple(float(x) for x in _dw)
    photonics_params = _params({**photonics_kw, "family": family}, NitrideNanowirePhotonicsParams)
    injector_params_base = _params(n["injector"], NitrideNanowireInjectorParams)

    # MEDIUM 5 (fix-3): tracked outside the try block so the catch-all
    # exception handler below can report whatever of these was ALREADY
    # genuinely resolved before the exception hit, instead of dropping them
    # to the row-level NaN/False defaults regardless of progress (previous
    # round's finding: a card failing partway through (e.g. x_in=1.0) could
    # report T_j=NaN even though wire_operating_point had already converged
    # at a real T_j).
    tj = _nan(); _op_iterations = _nan(); _thermal_converged = False; _field_kVcm = _nan()
    try:
        # H4: h_nu_eV for the electrical heat balance must be the row's own
        # E_X, not a fixed 2.0 eV literal -- seed it with a levels() solve at
        # T_hs/field=ext_field_kVcm (the wire's own photon energy, unknown
        # only up to the (weak) T_hs->T_j shift resolved next), then solve
        # the thermal operating point with that seed.
        lv_seed = levels(sysA, T_hs)
        if not lv_seed.valid: return inv(list(lv_seed.invalid_reasons))
        op = wire_operating_point(diode, I_uA=d.drive.I_uA, T_hs_K=T_hs, duty=duty,
                                   Rth_K_W=float(thermal["Rth_K_W"]), eta_total=eta_total,
                                   h_nu_eV=lv_seed.E_X_eV)
        _op_iterations = op["iterations"]
        # MEDIUM 7 (fix-2 round 2): thermal_converged must be a real False on
        # a thermal-solve failure, never left at the row-level NaN default.
        if not op["valid"]: return inv(op["reasons"], thermal_iterations=op["iterations"], thermal_converged=False)
        tj = op["T_j_K"]
        _thermal_converged = True

        lv = levels(sysA, tj)
        if not lv.valid: return inv(list(lv.invalid_reasons), T_j=tj, thermal_iterations=op["iterations"], thermal_converged=bool(op["valid"]))
        surf = surface_rates(surface_params, core_radius_nm=core, T_K=tj)
        # LOW 10 (fix-3): no pre-feedback rates() gate here any more -- this
        # used to call rates(lv, ...) purely to bail out early on an invalid
        # rate solve before evaluate_injection() ran, but evaluate_injection
        # never consumes that result (only lv.dE_e/h_meV and lv.overlap_sq,
        # both already available), and the SAME validity check runs again,
        # unconditionally, on rr2 below (from the post-field-feedback lv2)
        # moments later. Removed the redundant duplicate solve; rr2's own
        # check is now the single gate ("use the gate result").
        #
        # M16: the optical reservoir is the axial GaN barrier's own bulk
        # band edge (Varshni, minus the same 25 meV localization/Urbach
        # downshift the planar path uses), reusing device.py's already-
        # tagged helper -- not an untagged E_X+0.05 offset.  Local import:
        # avoids a module-load-time cycle with fsim_core.device, which
        # dispatches here locally in the other direction.
        from .device import _nitride_reservoir_energy_eV, _nitride_flat_background_acceptance
        # MEDIUM 4 (fix-3): an explicit nitride.reservoir_energy_eV card
        # override is honored through the SAME planar helper/"background"
        # dict convention device.py's own nitride path uses (including its
        # own finite/positive validation) -- previously allow-listed at the
        # top level but always discarded (background={} unconditionally).
        _bg_kw = {"reservoir_energy_eV": float(n["reservoir_energy_eV"])} if "reservoir_energy_eV" in n else {}
        reservoir_energy_eV = _nitride_reservoir_energy_eV({}, tj, _bg_kw)

        inj = evaluate_injection(diode, I_uA=d.drive.I_uA, T_K=tj, tau_pulse_ns=pulse,
            E_X_eV=lv.E_X_eV, reservoir_energy_eV=reservoir_energy_eV,
            barrier_e_eV=lv.dE_e_meV / 1000, barrier_h_eV=lv.dE_h_meV / 1000,
            surface_reservoir_ns=surf["k_surface_reservoir_ns"], tau_cap_ps=tau_cap_ps,
            S_dot=lv.overlap_sq, w_meV=float(d.dot.gamma300),
            eta_rad_matrix=1.0,  # [A] no separate card leaf; conservative fully-radiative reservoir matrix assumption
            eta_total=eta_total)
        if not inj["valid"]: return inv(list(inj["reasons"]), T_j=tj, thermal_iterations=op["iterations"], thermal_converged=bool(op["valid"]))
        _field_kVcm = inj["depletion_field_kVcm"]

        # M9/HIGH 2: feed the resolved depletion field, PLUS the card's own
        # independent external_field_kVcm baseline, back into the level
        # Hamiltonian and iterate once (one feedback solve is sufficient --
        # transport is already evaluated at T_j; see the module docstring).
        sysB = replace(sysA, external_field_kVcm=ext_field_kVcm + inj["depletion_field_kVcm"])
        lv2 = levels(sysB, tj)
        if not lv2.valid: return inv(list(lv2.invalid_reasons), T_j=tj, field_kVcm=inj["depletion_field_kVcm"], thermal_iterations=op["iterations"], thermal_converged=bool(op["valid"]))
        # M10 (constraint-4 option b): rates() k_nr_ns is INTRINSIC occupied-
        # dot loss only; surface's occupied-dot channel is added by hand
        # below (k_X_ns/k_XX_ns), never passed into rates() itself, since
        # the committed levels.rates() does not double k_nr for XX. This is
        # now the ONLY rates() solve in the row (LOW 10) -- its own validity
        # check just below is the single gate for the whole physics chain.
        rr2 = rates(lv2, tj, tau_rad0_ns=tau_rad0_ns, tau_cap_ps=tau_cap_ps,
                    reservoir_length_nm=15.0,  # [A DEVICE-PIECE CONSTRAINT 2] PER SIDE; rates() doubles it internally
                    k_nr_ns=float(n["dot"].get("k_intrinsic_ns", 0.0)))
        if not rr2["valid"]: return inv(list(rr2["invalid_reasons"]), T_j=tj, thermal_iterations=op["iterations"], thermal_converged=bool(op["valid"]))

        ph = response(photonics_params,
                      lambda_nm=lv2.lambda_nm, outer_radius_nm=outer,
                      gamma_X0_ns=rr2["gamma_X0_ns"], gamma_XX0_ns=rr2["gamma_XX0_ns"])
        if not ph["valid"]: return inv(list(ph["invalid_reasons"]), T_j=tj, thermal_iterations=op["iterations"], thermal_converged=bool(op["valid"]))

        gamma_x, gamma_xx = ph["gamma_X_ns"], ph["gamma_XX_ns"]
        kx = rr2["k_X_ns"] + surf["k_surface_X_ns"]
        kxx = rr2["k_XX_ns"] + surf["k_surface_XX_ns"]

        # MEDIUM 8 (fix-2 round 2): the genuine FilterBlock spectral
        # transmission, gated on d.filter.enabled exactly as the legacy
        # (non-nitride) evaluate() path gates it (device.py:1643) -- never
        # the planar NITRIDE path's cr["eta_out"] (a fixed cavity
        # out-coupling constant with no nanowire equivalent) and never an
        # invented closed form.  delta_xx=0.0: this platform has no
        # nitride.dot.delta_xx leaf (no modeled X/XX spectral offset), so
        # t_x == t_xx and a single eta_out serves both channels below.
        if not d.filter.enabled:
            eta_out = 1.0  # [DR] filter disabled: full transmission (legacy path's w=None convention)
            bg_accept = 1.0  # LOW 6: no filter -> no detuning-dependent suppression either
        else:
            w_val = (d.dot.gamma300 * d.filter.auto_w_scale) if d.filter.auto_w else d.filter.w
            eta_out = float(epsilon(0.0, d.dot.gamma300, d.dot.gamma300, w=w_val, kappa=None, dx=d.filter.dx).t_x)
            # LOW 6 (fix-3): the reservoir background sits reservoir_offset_meV
            # away from the X line (typically hundreds of meV -- the GaN
            # barrier edge vs. the confined transition), so it must NOT be
            # attenuated by the on-resonance X-line eta_out (that made rho/g2
            # filter-invariant to the reservoir's own detuning, i.e. wrong
            # whenever the filter is on). Reused from the planar path's own
            # detuning-aware flat-spectrum acceptance helper
            # (device._nitride_flat_background_acceptance) instead of an
            # invented closed form. This platform's own filter has no cavity
            # kappa (kappa=None, a plain top-hat slit above); [A] substitute
            # kappa_meV=w_val (the slit's own width, the only spectral scale
            # this platform's filter has) and detuning_meV=0.0 (no cavity
            # mode to detune from) so the SAME Lorentzian-kernel average-
            # acceptance math the planar cavity path uses still correctly
            # rejects the far-detuned reservoir once the filter narrows below
            # the reservoir offset.
            reservoir_offset_meV = (reservoir_energy_eV - lv2.E_X_eV) * 1e3
            bg_accept = float(_nitride_flat_background_acceptance(
                w_val, w_val, d.filter.dx, 0.0, reservoir_offset_meV))
        ecx, ecxx = ph["eta_collection_X"], ph["eta_collection_XX"]

        # H3 fix: eta_collection_X/XX is passed as the pulse-counting
        # t_X/t_XX argument -- applied ONCE, inside the ODE propagation
        # (mean_counts_x/xx below are therefore already the COLLECTED
        # counts) -- eta_out is then applied ONCE more, afterward, exactly
        # like the planar path's separate cr["eta_out"] * cnt["mean_counts"]
        # step. Neither factor is applied twice.
        if regime == "deterministic_pair":
            count = pulse_counting.deterministic_cycle_g2(gamma_x, gamma_xx, kx, kxx, ecx, ecxx, period, eta_load=d.drive.eta_load, gate_ns=gate, split=True)
        else:
            count = pulse_counting.pulse_g2(inj["r_captured_s"] * 1e-9, gamma_x, gamma_xx, kx, kxx, ecx, ecxx, pulse, period - pulse, pump_ratio=d.drive.cw_pump_ratio, gate_ns=gate, split=True)

        mx = count.get("mean_counts_x", _nan()); mxx = count.get("mean_counts_xx", _nan())
        sx = eta_out * mx; sxx = eta_out * mxx
        signal = sx + sxx
        bg_window_ns = min(gate, pulse)
        # LOW 16 (fix-2 round 2): the reservoir background channel gets the
        # SAME (geometric) collection factor as the X channel it is compared
        # against, once, in addition to the spectral filter transmission --
        # "Composition rules for the device piece" bullet 9.  Numerically
        # near-inert at contract-default cards (the reservoir rate is
        # already many orders below the X line), but no longer silently
        # skipped.
        bg = inj["accepted_background_s"] * bg_window_ns * 1e-9 * bg_accept * ecx
        if math.isfinite(sx): bg += float(d.drive.b_res) * sx
        rho = signal / (signal + bg) if signal + bg > 0 else _nan()
        g2 = 1.0 - rho * rho * (1.0 - count["g2"]) if math.isfinite(rho) and math.isfinite(count.get("g2", _nan())) else _nan()
        # LOW 8 (fix-3): a non-finite signal must yield a NaN (not-computed)
        # flux, never a silent 0 -- max(0.0, nan) == 0.0 in Python (nan
        # comparisons are always False), which previously reported "zero
        # photons collected" instead of "not computed" for an unreachable-
        # in-practice non-finite-signal state.
        flux = _nan() if not math.isfinite(signal) else rep_v * max(0.0, signal)

        st = set_feasibility(tj, radius_nm=disc, eps_r=float(setp.get("eps_r", 9.5)),  # [A] GaN static relative permittivity default; card may override via drive.set_params.eps_r
                              R_T_ohm=R_T_ohm, ec_margin=float(setp.get("ec_margin", 10.0)),  # [A] charge-quantization margin default (drive_mech.set_feasibility's own convention)
                              f_cycle_Hz=rep_v)
        # Max island radius still satisfying the charging-energy screen,
        # inverted through the SAME e^2/C isolated-sphere convention
        # set_feasibility itself uses (mirrors the planar path).
        ec_margin = float(setp.get("ec_margin", 10.0))  # [A] see above
        C_sigma_max = E_SI ** 2 / (ec_margin * KB_SI * tj)
        st["radius_max_nm"] = float(island_radius_nm(C_sigma_max, float(setp.get("eps_r", 9.5))))  # [A] see above

        ip = replace(injector_params_base,
                     reservoir_state_count_e=rr2["reservoir_state_count_e"],
                     reservoir_state_count_h=rr2["reservoir_state_count_h"])
        rti = injector_feasibility(ip, T_K=tj, rep_rate_hz=rep_v, loading_window_ns=pulse,
            electron_level_eV=-lv2.dE_e_meV / 1000, hole_level_eV=-lv2.dE_h_meV / 1000,
            electron_spacing_meV=lv2.sp_split_e_meV, hole_spacing_meV=lv2.sp_split_h_meV,
            second_pair_addition_meV=st["E_C_meV"], available_pair_rate_Hz=inj["r_captured_s"],
            # LOW 12 (fix-3, documentation only): deliberately the TRANSPORT
            # DEPLETION FIELD ALONE, not ext_field_kVcm + depletion_field_kVcm
            # (unlike sysB's level feedback above, which DOES add the two --
            # see HIGH 2/MEDIUM 4). DEVICE-PIECE addendum: "field_kVcm = the
            # transport depletion field" for the injector call specifically;
            # not an omission to "fix" into consistency with the levels
            # feedback convention.
            field_kVcm=inj["depletion_field_kVcm"],
            gate_ns=(gate if regime == "deterministic_pair" else pulse))

        delivery = pulse_delivery(diode, V_j=inj["V_j"], T_K=tj, tau_pulse_ns=pulse, rep_rate_hz=rep_v, C_parasitic_F=float(thermal["C_parasitic_F"]))

        one = bool(count.get("one_pair_valid", False)); supply = bool(inj["r_captured_s"] >= rep_v)
        counting_converged = bool(count.get("converged", False))
        reasons = []
        if not counting_converged: reasons.append("photon counting did not converge")
        if not math.isfinite(g2): reasons.append("non-finite g2 at operating point")
        valid = bool(counting_converged and math.isfinite(g2))
        optical = bool(valid and g2 < .5 and flux >= 1000 and (regime != "deterministic_pair" or (one and supply)))

        # M12 (non-contract bonus diagnostic): a vertical row above the
        # LP11 single-mode cutoff, or with nonzero approximation_error, is
        # never nominated as a headline row (contract "Family-specific
        # interfaces" / DEVICE-PIECE CONSTRAINT re: photonics).
        if family == "vertical_photonic":
            headline_eligible = bool(valid and ph.get("single_mode") is True and ph.get("approximation_error", 1.0) == 0.0)
        else:
            headline_eligible = bool(valid)

        row = {}
        row.update(lv2.__dict__); row.update(rr2); row.update(surf); row.update(ph)
        row.update(inj); row.update(rti)
        row.update({
            "platform": "ingan_gan_nanowire", "family": family, "T_hs": float(T_hs), "T_j": tj,
            "thermal_iterations": op["iterations"], "thermal_converged": bool(op["valid"]),
            "cycle_loading": regime, "rep_rate_hz": rep_v,
            "strain_bound": strain_bound, "bound_role": bound_role, "bound_reversal": "not_computed",
            "screening_fraction": screening, "core_radius_nm": core, "outer_radius_nm": outer,
            "conducting_radius_nm": conducting,
            "external_field_kVcm": ext_field_kVcm,
            "reservoir_state_count_e": rr2["reservoir_state_count_e"], "reservoir_state_count_h": rr2["reservoir_state_count_h"],
            # M10 (constraint-4 option b): the reported k_X_ns/k_XX_ns are the
            # FULL post-capture survival rates (intrinsic + occupied-dot
            # surface loss, the latter added once by hand above) -- rr2's own
            # k_X_ns/k_XX_ns (intrinsic-only, from the row.update(rr2) merge
            # above) are deliberately overridden here, not left stale.
            "k_X_ns": kx, "k_XX_ns": kxx,
            "valid": valid, "invalid_reasons": reasons, "g2_op": g2,
            "counting_converged": counting_converged,
            "collected_flux_pulsed_s": flux, "collected_flux_delivered_s": flux * delivery["delivered_step_fraction"],
            "collected_flux_x_s": rep_v * sx, "collected_flux_xx_s": rep_v * sxx,
            "background_flux_s": rep_v * bg, "total_detected_flux_s": rep_v * (signal + bg),
            "mean_counts": count.get("mean_counts", _nan()), "mean_counts_x": mx, "mean_counts_xx": mxx,
            "gate_ns_used": gate,
            "S_X": gamma_x / (gamma_x + kx) if (gamma_x + kx) > 0 else _nan(),
            "S_XX": gamma_xx / (gamma_xx + kxx) if (gamma_xx + kxx) > 0 else _nan(),
            "rho_pulsed": rho, "blocked_load_probability": count.get("blocked_load_probability", _nan()),
            "eta_out": eta_out,
            "one_pair_valid": one, "pair_supply_possible": supply,
            "optical_pass": optical,
            "set_feasible": st["feasible"] if regime == "deterministic_pair" else "not_applicable",
            "set_E_C_meV": st["E_C_meV"], "set_EC_over_kT": st["EC_over_kT"],
            "set_radius_nm": st["radius_nm"], "set_radius_max_nm": st["radius_max_nm"],
            "set_C_sigma_F": st["C_sigma_F"], "set_R_T_over_RQ": st["R_T_over_RQ"],
            "set_f_max_Hz": st["f_max_Hz"],
            "hardware_qualified": bool(optical and regime == "deterministic_pair" and st["feasible"]),
            "rti_feasible": rti["rti_feasible"] if regime == "deterministic_pair" else "not_applicable",
            "rti_status": rti.get("rti_status") if regime == "deterministic_pair" else "not_applicable",
            "rti_qualified": bool(optical and regime == "deterministic_pair" and rti["rti_feasible"]),
            "device_pass": bool(optical and (st["feasible"] if regime == "deterministic_pair" else True)),
            "rti_device_pass": bool(optical and regime == "deterministic_pair" and rti["rti_feasible"]),
            "tau_rad_bare_ns": _lifetime_ns(rr2["gamma_X0_ns"]),
            "tau_rad_photonic_ns": _lifetime_ns(gamma_x),
            "tau_total_X_ns": _lifetime_ns(gamma_x + kx),
            "tau_RC_ns": delivery["tau_RC_s"] * 1e9,
            "delivered_step_fraction": delivery["delivered_step_fraction"],
            "pulse_delivery_feasible": delivery["pulse_delivery_feasible"],
            "headline_eligible": headline_eligible,
            "provenance": {"levels": lv2.provenance, "surface": surf["provenance"],
                           "transport": inj["provenance"], "photonics": ph["provenance"],
                           "injector": rti.get("provenance"), "set": st.get("notes")},
        })
        return row
    except (ValueError, ZeroDivisionError) as exc:
        # Numerical/physics invalidity ONLY (unbound state, non-convergence,
        # a downstream module's own runtime domain check) -- card-shape and
        # required-leaf errors were already raised above, outside this
        # try/except, and propagate to the caller unchanged.
        #
        # MEDIUM 5 (fix-3): keep whatever of T_j/thermal_iterations/
        # thermal_converged/field_kVcm was ALREADY genuinely resolved before
        # the exception hit (tracked above), instead of unconditionally
        # dropping to the row-level NaN/False defaults regardless of how far
        # the chain got.
        _extra = {"thermal_iterations": _op_iterations, "thermal_converged": _thermal_converged}
        if math.isfinite(tj): _extra["T_j"] = tj
        if math.isfinite(_field_kVcm): _extra["field_kVcm"] = _field_kVcm
        return inv([str(exc)], **_extra)


def _one(d, T_hs):
    """Public per-row entry point: the raw physics result (_one_raw),
    pre-filled with NaN for the full row-key set (_row_keys()) so that
    invalid rows carry EXACTLY the same keys a valid row would (HIGH 1)."""
    raw = _one_raw(d, T_hs)
    row = {k: _nan() for k in _row_keys()}
    row.update(raw)
    return row


def evaluate_nanowire(design, T_grid=None):
    Ts = np.asarray(T_grid if T_grid is not None else [design.thermal.T_hs], dtype=float)
    rows = [_one(design, float(t)) for t in Ts]
    scalars = min(rows, key=lambda r: abs(r["T_hs"] - design.thermal.T_hs)).copy()
    curves = {k: np.asarray([r.get(k, _nan()) for r in rows], dtype=float)
              for k in ("g2_op", "collected_flux_pulsed_s", "collected_flux_delivered_s", "T_j", "gamma_X_ns")}
    return {"curves": curves, "scalars": scalars}


def _finite(x):
    return isinstance(x, (int, float)) and math.isfinite(x)


def _bound_reversal(relaxed_scalars, unrelaxed_scalars):
    """M13/DEVICE-PIECE CONSTRAINT 5: bound_reversal is True when mu,
    collected_flux_pulsed_s, or g2_op orders the strain-bound pair AGAINST
    the headline_upper (relaxed) / conservative_lower (unrelaxed) label --
    never derived from E_X_eV.

    MEDIUM 9 (fix-2 round 2): a reversal is only DECLARED between two rows
    that both clear the optical validity floor this module itself gates
    on (valid AND collected_flux_pulsed_s >= 1000/s) -- below that floor a
    photon-count difference is noise, not a physical reversal (the
    default horizontal SET card's False-alarm case: the unrelaxed partner
    emits 204 photons/s, far below the 1000/s optical floor, yet an
    unguarded comparison still calls its g2_op ordering a "reversal").
    Returns "not_comparable" (not the single-row "not_computed" sentinel,
    which means "never attempted") whenever either partner fails that
    floor or is missing metrics altogether."""
    def _above_floor(s):
        return (s.get("valid") is True and _finite(s.get("collected_flux_pulsed_s"))
                and s["collected_flux_pulsed_s"] >= 1000.0)
    if not (_above_floor(relaxed_scalars) and _above_floor(unrelaxed_scalars)):
        return "not_comparable"
    checks = []
    for key, higher_is_headline in (("mu", True), ("collected_flux_pulsed_s", True), ("g2_op", False)):
        rv, uv = relaxed_scalars.get(key), unrelaxed_scalars.get(key)
        if not (_finite(rv) and _finite(uv)): continue
        checks.append(rv < uv if higher_is_headline else rv > uv)
    return any(checks) if checks else "not_comparable"


def evaluate_strain_pair(design, T_grid=None):
    """Evaluate both physical strain bounds without mutating the card [DR].
    Sets bound_reversal True/False/"not_comparable" on both scalars dicts
    (see _bound_reversal); the per-row single-evaluator bound_reversal stays
    the "not_computed" string sentinel (M13)."""
    import copy
    out = {}
    for bound in ("unrelaxed", "relaxed"):
        d = copy.deepcopy(design)
        d.nitride["nanowire"]["strain_bound"] = bound
        d.nitride["dot"]["strain_fraction"] = (0.0 if bound == "relaxed" else 1.0)
        out[bound] = evaluate_nanowire(d, T_grid)
    reversal = _bound_reversal(out["relaxed"]["scalars"], out["unrelaxed"]["scalars"])
    for value in out.values():
        value["scalars"]["bound_reversal"] = reversal
    return out
