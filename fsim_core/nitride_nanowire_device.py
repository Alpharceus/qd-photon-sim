"""Nanowire device integration (piece 7).  This module deliberately imports
no device module at module scope: :mod:`fsim_core.device` dispatches here
locally, avoiding a cycle; the reverse direction (this module reusing
device.py's ``_nitride_reservoir_energy_eV`` helper) is also a local,
function-scope import for the same reason.

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
    finding H3: the previous round applied ph["eta_collection_X"] both
    inside AND outside pulse_g2/deterministic_cycle_g2).
  - Scalars are the row nearest thermal.T_hs in a supplied T_grid (never
    rows[-1]).
"""
from __future__ import annotations

import importlib
import math
import time
from dataclasses import replace

import numpy as np

from . import pulse_counting
from .drive_mech import set_feasibility
from .loading import E_SI, KB_SI, island_radius_nm
from .nitride_nanowire_levels import NitrideNanowireSystem, levels, rates
from .nitride_nanowire_surface import NitrideNanowireSurfaceParams, surface_rates
from .nitride_nanowire_transport import (
    NitrideWireDiode, wire_operating_point, evaluate_injection, pulse_delivery,
)
from .nitride_nanowire_photonics import NitrideNanowirePhotonicsParams, response

# The scalar schema is deliberately local rather than inferred from whatever
# happened to run.  Invalid points therefore retain their coordinates and all
# contract columns, preventing a previous-temperature result leaking into a
# failed point.  NaN is the unavailable-computation sentinel [DR].  This is
# the contract's exact 144-column "Row columns" list (docs/
# nitride_nanowire_contract.md; the verifier parses that section
# independently and asserts this set is a superset of it -- H1/H8) plus one
# bonus, non-contract diagnostic (headline_eligible, M12) appended last.
_ROW_KEYS = ("platform family T_hs T_j thermal_iterations thermal_converged cycle_loading rep_rate_hz strain_bound bound_role bound_reversal screening_fraction core_radius_nm outer_radius_nm conducting_radius_nm E_X_eV lambda_nm field_kVcm overlap_sq electron_bound hole_bound E_a_meV gamma_X0_ns gamma_XX0_ns k_X_ns k_XX_ns escape_prefactor_ns tau_cap_ps_used reservoir_state_count_e reservoir_state_count_h sidewall_overlap radius_over_lambda V_number beta_HE11 eta_collection_X eta_collection_XX radiative_rate_factor gamma_X_ns gamma_XX_ns approximation_error single_mode beta_multimode_penalty degree_of_linear_polarization antenna_rate_factor k_side_ns k_surface_reservoir_ns k_surface_X_ns k_surface_XX_ns shell_multiplier_used reservoir_access_used occupied_dot_access_used area_cm2 J_A_cm2 V_j V_terminal depletion_field_kVcm C_dep_F eta_inj f_capture f_qfl_dot f_qfl_dot_thermodynamic_limit f_qfl_background r_supply_s r_captured_s r_matrix_radiative_s r_matrix_nonradiative_s r_surface_reservoir_s r_leakage_s mu power_on_W accounting_residual_s tau_RC_ns delivered_step_fraction pulse_delivery_feasible tau_rad_bare_ns tau_rad_photonic_ns tau_total_X_ns collected_flux_pulsed_s collected_flux_delivered_s collected_flux_x_s collected_flux_xx_s background_flux_s total_detected_flux_s mean_counts mean_counts_x mean_counts_xx gate_ns_used S_X S_XX rho_pulsed blocked_load_probability counting_converged eta_out g2_op one_pair_valid pair_supply_possible valid invalid_reasons provenance set_feasible set_E_C_meV set_EC_over_kT set_radius_nm set_radius_max_nm set_C_sigma_F set_R_T_over_RQ set_f_max_Hz rti_feasible rti_status rti_transport_feasible rti_level_margin_kT rti_orbital_margin_kT rti_alignment_error_meV rti_linewidth_meV rti_rate_Hz rti_e_rate_Hz rti_h_rate_Hz rti_bypass_fraction rti_missed_load_probability rti_second_pair_probability rti_growth_feasible rti_failed_checks rti_evidence_status rti_level_margin_e_kT rti_level_margin_h_kT rti_alignment_error_e_meV rti_alignment_error_h_meV rti_required_bias_shift_meV rti_bypass_fraction_tsai_partition rti_second_carrier_probability rti_growth_nearest_commensurate_nm rti_growth_perturbed_margins_kT rti_barrier_polarization_tilt_eV rti_well_to_dot_drop_meV rti_p_free_cm3 rti_numerics_ok rti_gate_ns rti_reservoir_state_count_e rti_reservoir_state_count_h optical_pass hardware_qualified rti_qualified device_pass rti_device_pass headline_eligible").split()

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


def _nan(): return float("nan")


def _lifetime_ns(rate_ns):
    """1/rate, distinguishing a genuinely zero rate (infinite lifetime) from
    a non-finite/not-computed one (NaN), matching device.py's
    _nitride_lifetime_ns convention (interface constraints: "infinity only
    for a genuinely zero decay rate")."""
    if not isinstance(rate_ns, (int, float)) or not math.isfinite(rate_ns):
        return _nan()
    return float("inf") if rate_ns <= 0 else 1.0 / rate_ns


def _invalid(T, reasons, **extra):
    n = _nan()
    row = {k: n for k in _ROW_KEYS}
    row.update({"platform": "ingan_gan_nanowire", "T_hs": float(T), "T_j": n,
           "valid": False, "invalid_reasons": list(reasons), "g2_op": n,
           "collected_flux_pulsed_s": n, "collected_flux_delivered_s": n,
           "optical_pass": False, "hardware_qualified": False,
           "rti_qualified": False, "device_pass": False, "rti_device_pass": False,
           "headline_eligible": False,
           "provenance": {}, "set_feasible": "not_applicable", "rti_feasible": "not_applicable",
           "rti_status": "not_applicable", "rti_failed_checks": [], "bound_reversal": "not_computed"})
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


# --------------------------------------------------------------- injector
# The injector module (fsim_core/nitride_nanowire_injector.py) is being
# edited concurrently by another worker at the time this piece is
# implemented.  Its public injector_feasibility(params, ...) signature is
# stable (contract Module table) and handles both carriers internally, so
# no other change is needed here -- but the import itself can transiently
# fail mid-edit, so it is done lazily (first call, not module load) and
# retried briefly instead of failing the whole device module at import
# time.  Card leaves for the injector are passed through generically
# (``NitrideNanowireInjectorParams(**mapping)``) so a new field the
# concurrent round adds does not require a change here.
_injector_cache: dict = {}


def _injector_symbols():
    if "syms" in _injector_cache:
        return _injector_cache["syms"]
    last_exc = None
    for attempt in range(6):
        try:
            mod = importlib.import_module(".nitride_nanowire_injector", __package__)
            syms = (mod.NitrideNanowireInjectorParams, mod.injector_feasibility)
            _injector_cache["syms"] = syms
            return syms
        except Exception as exc:  # noqa: BLE001 -- transient mid-edit import failure, retried below
            last_exc = exc
            if attempt < 5:
                time.sleep(0.5)
    raise ImportError(
        "fsim_core.nitride_nanowire_injector failed to import after retries "
        f"(concurrent edit?): {last_exc!r}")


def _validate(d):
    n = d.nitride
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
    if nw.get("strain_bound") not in ("relaxed", "unrelaxed"): raise ValueError("nanowire strain_bound required")
    sf = dot.get("strain_fraction")
    expected = 0.0 if nw["strain_bound"] == "relaxed" else 1.0
    if sf is not None and float(sf) != expected: raise ValueError("nitride.dot.strain_fraction contradicts strain_bound")
    core = float(nw["core_radius_nm"]); disc = float(dot.get("radius_nm", core))
    if family == "horizontal_as_built" and disc != core: raise ValueError("horizontal disc radius must equal core radius")
    if family == "vertical_photonic" and disc >= core: raise ValueError("vertical photonic disc radius must be less than core radius")
    outer = float(nw.get("outer_radius_nm", core))
    if outer < core: raise ValueError("outer radius below core")
    setp = d.drive.set_params
    if setp.get("radius_nm", disc) != disc or "C_sigma_F" in setp: raise ValueError("SET radius must be the disc radius; C_sigma override forbidden")
    conducting = float(d.drive.diode.get("conducting_radius_nm", core))
    if not (0.0 < conducting <= core): raise ValueError("drive.diode.conducting_radius_nm must be in (0, core_radius_nm]")
    bound_role = "headline_upper" if nw["strain_bound"] == "relaxed" else "conservative_lower"
    return family, core, disc, outer, conducting, d.drive.cycle_loading, float(d.drive.rep_rate_hz), nw["strain_bound"], bound_role, float(dot.get("screening_fraction", 0.0))


def _one(d, T_hs):
    family, core, disc, outer, conducting, regime, rep, strain_bound, bound_role, screening = _validate(d)
    n = d.nitride
    thermal = n["wire_thermal"]

    def inv(reasons, **extra):
        return _invalid(T_hs, reasons, family=family, cycle_loading=regime, rep_rate_hz=rep,
                         strain_bound=strain_bound, bound_role=bound_role, screening_fraction=screening,
                         core_radius_nm=core, outer_radius_nm=outer, conducting_radius_nm=conducting, **extra)

    # Card-shape / required-leaf validation and dataclass construction from
    # raw card leaves happen OUTSIDE the numerical try/except below: a
    # missing leaf or a malformed sub-block value must be REJECTED (raise),
    # never silently downgraded to an invalid row (acceptance criterion 1).
    # Only the PHYSICS chain (thermal solve, levels, rates, injection,
    # counting, ...) -- operating on already-validated inputs -- is allowed
    # to fail into an explicit invalid row.
    rep_v = rep; pulse = float(d.drive.diode.get("tau_pulse_ns", .1)); duty = float(d.drive.duty)
    if rep_v <= 0 or pulse <= 0 or pulse * rep_v * 1e-9 > 1:
        return inv(["invalid pulse/rep-rate"])
    period = 1e9 / rep_v; gate = float(d.drive.gate_ns if d.drive.gate_ns is not None else period)

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
    diode_raw.setdefault("R_s_ohm", float(thermal["R_s_ohm"]))
    diode_raw.setdefault("f_Rs_local", float(thermal["f_Rs_local"]))
    diode = NitrideWireDiode(core_radius_nm=core, conducting_radius_nm=conducting,
        barrier_left_nm=barrier_left, barrier_right_nm=barrier_right,
        d_active_nm=height_nm, x_in=x_in, **diode_raw)

    sys_kw = dict(height_nm=height_nm, core_radius_nm=core, outer_radius_nm=outer,
                  disc_radius_nm=(None if family == "horizontal_as_built" else disc),
                  x_in=x_in, strain_bound=strain_bound, screening_fraction=screening,
                  external_field_kVcm=0.0)
    for opt in ("vbo_InN_GaN_eV", "strain_c_fraction"):
        if opt in n["dot"]: sys_kw[opt] = float(n["dot"][opt])
    sysA = NitrideNanowireSystem(**sys_kw)

    eta_total = float(thermal["eta_total"])
    surface_params = _params(n["surface"], NitrideNanowireSurfaceParams)
    photonics_params = _params({**n["photonics"], "family": family}, NitrideNanowirePhotonicsParams)
    InjectorParams, injector_feasibility = _injector_symbols()
    injector_params_base = _params(n["injector"], InjectorParams)

    try:
        # H4: h_nu_eV for the electrical heat balance must be the row's own
        # E_X, not a fixed 2.0 eV literal -- seed it with a levels() solve at
        # T_hs/field=0 (the wire's own photon energy, unknown only up to the
        # (weak) T_hs->T_j shift resolved next), then solve the thermal
        # operating point with that seed.
        lv_seed = levels(sysA, T_hs)
        if not lv_seed.valid: return inv(list(lv_seed.invalid_reasons))
        op = wire_operating_point(diode, I_uA=d.drive.I_uA, T_hs_K=T_hs, duty=duty,
                                   Rth_K_W=float(thermal["Rth_K_W"]), eta_total=eta_total,
                                   h_nu_eV=lv_seed.E_X_eV)
        if not op["valid"]: return inv(op["reasons"], thermal_iterations=op["iterations"])
        tj = op["T_j_K"]

        lv = levels(sysA, tj)
        if not lv.valid: return inv(list(lv.invalid_reasons), T_j=tj, thermal_iterations=op["iterations"], thermal_converged=bool(op["valid"]))
        surf = surface_rates(surface_params, core_radius_nm=core, T_K=tj)
        # M10 (constraint-4 option b): rates() k_nr_ns is INTRINSIC occupied-
        # dot loss only; surface's occupied-dot channel is added by hand
        # below (k_X_ns/k_XX_ns), never passed into rates() itself, since
        # the committed levels.rates() does not double k_nr for XX.
        rr = rates(lv, tj, tau_rad0_ns=tau_rad0_ns, tau_cap_ps=tau_cap_ps,
                   reservoir_length_nm=15.0,  # [A DEVICE-PIECE CONSTRAINT 2] PER SIDE; rates() doubles it internally
                   k_nr_ns=float(n["dot"].get("k_intrinsic_ns", 0.0)))
        if not rr["valid"]: return inv(list(rr["invalid_reasons"]), T_j=tj, thermal_iterations=op["iterations"], thermal_converged=bool(op["valid"]))

        # M16: the optical reservoir is the axial GaN barrier's own bulk
        # band edge (Varshni, minus the same 25 meV localization/Urbach
        # downshift the planar path uses), reusing device.py's already-
        # tagged helper -- not an untagged E_X+0.05 offset.  Local import:
        # avoids a module-load-time cycle with fsim_core.device, which
        # dispatches here locally in the other direction.
        from .device import _nitride_reservoir_energy_eV
        reservoir_energy_eV = _nitride_reservoir_energy_eV({}, tj, {})

        inj = evaluate_injection(diode, I_uA=d.drive.I_uA, T_K=tj, tau_pulse_ns=pulse,
            E_X_eV=lv.E_X_eV, reservoir_energy_eV=reservoir_energy_eV,
            barrier_e_eV=lv.dE_e_meV / 1000, barrier_h_eV=lv.dE_h_meV / 1000,
            surface_reservoir_ns=surf["k_surface_reservoir_ns"], tau_cap_ps=tau_cap_ps,
            S_dot=lv.overlap_sq, w_meV=float(d.dot.gamma300),
            eta_rad_matrix=1.0,  # [A] no separate card leaf; conservative fully-radiative reservoir matrix assumption
            eta_total=eta_total)
        if not inj["valid"]: return inv(list(inj["reasons"]), T_j=tj, thermal_iterations=op["iterations"], thermal_converged=bool(op["valid"]))

        # M9: feed the resolved depletion field back into the level
        # Hamiltonian and iterate once (one feedback solve is sufficient --
        # transport is already evaluated at T_j; see the module docstring).
        sysB = replace(sysA, external_field_kVcm=inj["depletion_field_kVcm"])
        lv2 = levels(sysB, tj)
        if not lv2.valid: return inv(list(lv2.invalid_reasons), T_j=tj, field_kVcm=inj["depletion_field_kVcm"], thermal_iterations=op["iterations"], thermal_converged=bool(op["valid"]))
        rr2 = rates(lv2, tj, tau_rad0_ns=tau_rad0_ns, tau_cap_ps=tau_cap_ps,
                    reservoir_length_nm=15.0, k_nr_ns=float(n["dot"].get("k_intrinsic_ns", 0.0)))
        if not rr2["valid"]: return inv(list(rr2["invalid_reasons"]), T_j=tj, thermal_iterations=op["iterations"], thermal_converged=bool(op["valid"]))

        ph = response(photonics_params,
                      lambda_nm=lv2.lambda_nm, outer_radius_nm=outer,
                      gamma_X0_ns=rr2["gamma_X0_ns"], gamma_XX0_ns=rr2["gamma_XX0_ns"])
        if not ph["valid"]: return inv(list(ph["invalid_reasons"]), T_j=tj, thermal_iterations=op["iterations"], thermal_converged=bool(op["valid"]))

        gamma_x, gamma_xx = ph["gamma_X_ns"], ph["gamma_XX_ns"]
        kx = rr2["k_X_ns"] + surf["k_surface_X_ns"]
        kxx = rr2["k_XX_ns"] + surf["k_surface_XX_ns"]

        # M17: anchored-Lorentzian collection window transmission, applied
        # ONCE per channel below (H3/M17) -- the same closed-form Lorentzian
        # (2*dx/w) window the planar path's cavity/filter response also
        # reduces to when Q/kappa play no role (there is no cavity here).
        eta_out = 1.0 if not d.filter.enabled else 1.0 / (1.0 + (2.0 * d.filter.dx / (d.filter.w if not d.filter.auto_w else d.dot.gamma300)) ** 2)
        ecx, ecxx = ph["eta_collection_X"], ph["eta_collection_XX"]

        # H3 fix: eta_collection_X/XX is passed as the pulse-counting
        # t_X/t_XX argument -- applied ONCE, inside the ODE propagation
        # (mean_counts_x/xx below are therefore already the COLLECTED
        # counts) -- eta_out is then applied ONCE more, afterward, exactly
        # like the planar path's separate cr["eta_out"] * cnt["mean_counts"]
        # step. Neither factor is applied twice.
        if regime == "deterministic_pair":
            count = pulse_counting.deterministic_cycle_g2(gamma_x, gamma_xx, kx, kxx, ecx, ecxx, period, eta_load=1.0, gate_ns=gate, split=True)
        else:
            count = pulse_counting.pulse_g2(inj["r_captured_s"] * 1e-9, gamma_x, gamma_xx, kx, kxx, ecx, ecxx, pulse, period - pulse, pump_ratio=d.drive.cw_pump_ratio, gate_ns=gate, split=True)

        mx = count.get("mean_counts_x", _nan()); mxx = count.get("mean_counts_xx", _nan())
        sx = eta_out * mx; sxx = eta_out * mxx
        signal = sx + sxx
        bg_window_ns = min(gate, pulse)
        bg = inj["accepted_background_s"] * bg_window_ns * 1e-9 * eta_out
        if math.isfinite(sx): bg += float(d.drive.b_res) * sx
        rho = signal / (signal + bg) if signal + bg > 0 else _nan()
        g2 = 1.0 - rho * rho * (1.0 - count["g2"]) if math.isfinite(rho) and math.isfinite(count.get("g2", _nan())) else _nan()
        flux = rep_v * max(0.0, signal)

        st = set_feasibility(tj, radius_nm=disc, eps_r=float(setp.get("eps_r", 9.5)),
                              R_T_ohm=R_T_ohm, ec_margin=float(setp.get("ec_margin", 10.0)), f_cycle_Hz=rep_v)
        # Max island radius still satisfying the charging-energy screen,
        # inverted through the SAME e^2/C isolated-sphere convention
        # set_feasibility itself uses (mirrors the planar path).
        ec_margin = float(setp.get("ec_margin", 10.0))
        C_sigma_max = E_SI ** 2 / (ec_margin * KB_SI * tj)
        st["radius_max_nm"] = float(island_radius_nm(C_sigma_max, float(setp.get("eps_r", 9.5))))

        ip = replace(injector_params_base,
                     reservoir_state_count_e=rr2["reservoir_state_count_e"],
                     reservoir_state_count_h=rr2["reservoir_state_count_h"])
        rti = injector_feasibility(ip, T_K=tj, rep_rate_hz=rep_v, loading_window_ns=pulse,
            electron_level_eV=-lv2.dE_e_meV / 1000, hole_level_eV=-lv2.dE_h_meV / 1000,
            electron_spacing_meV=lv2.sp_split_e_meV, hole_spacing_meV=lv2.sp_split_h_meV,
            second_pair_addition_meV=st["E_C_meV"], available_pair_rate_Hz=inj["r_captured_s"],
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

        row = {k: _nan() for k in _ROW_KEYS}
        row.update(lv2.__dict__); row.update(rr2); row.update(surf); row.update(ph)
        row.update(inj); row.update(rti)
        row.update({
            "platform": "ingan_gan_nanowire", "family": family, "T_hs": float(T_hs), "T_j": tj,
            "thermal_iterations": op["iterations"], "thermal_converged": bool(op["valid"]),
            "cycle_loading": regime, "rep_rate_hz": rep_v,
            "strain_bound": strain_bound, "bound_role": bound_role, "bound_reversal": "not_computed",
            "screening_fraction": screening, "core_radius_nm": core, "outer_radius_nm": outer,
            "conducting_radius_nm": conducting,
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
        return inv([str(exc)])


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
    never derived from E_X_eV. "not_computed" (the single-row sentinel) is
    preserved when either partner is invalid (all-NaN), so an unpaired row
    is never scored as a false reversal."""
    checks = []
    for key, higher_is_headline in (("mu", True), ("collected_flux_pulsed_s", True), ("g2_op", False)):
        rv, uv = relaxed_scalars.get(key), unrelaxed_scalars.get(key)
        if not (_finite(rv) and _finite(uv)): continue
        checks.append(rv < uv if higher_is_headline else rv > uv)
    return any(checks) if checks else "not_computed"


def evaluate_strain_pair(design, T_grid=None):
    """Evaluate both physical strain bounds without mutating the card [DR].
    Sets bound_reversal True/False on both scalars dicts (see
    _bound_reversal); the per-row single-evaluator bound_reversal stays the
    "not_computed" string sentinel (M13)."""
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
