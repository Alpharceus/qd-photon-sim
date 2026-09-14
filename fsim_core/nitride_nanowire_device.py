"""Nanowire device integration.  This module deliberately imports no device
module: :mod:`fsim_core.device` dispatches here locally, avoiding a cycle.

The 1.0 ns bare lifetime and the 0.1 ns electrical pulse are [A] design
inputs, not fits to Deshpande et al. (APL 2014).  Material/geometry inputs
are supplied by the card and retain the provenance of pieces 2--6.
"""
from dataclasses import replace
import math
import numpy as np

from . import pulse_counting
from .drive_mech import set_feasibility
from .nitride_nanowire_levels import NitrideNanowireSystem, levels, rates
from .nitride_nanowire_surface import NitrideNanowireSurfaceParams, surface_rates
from .nitride_nanowire_transport import NitrideWireDiode, wire_operating_point, evaluate_injection, pulse_delivery
from .nitride_nanowire_photonics import NitrideNanowirePhotonicsParams, response
from .nitride_nanowire_injector import NitrideNanowireInjectorParams, injector_feasibility


def _nan(): return float("nan")


def _invalid(T, reasons, **extra):
    n = _nan()
    row = {"platform":"ingan_gan_nanowire", "T_hs":float(T), "T_j":n,
           "valid":False, "invalid_reasons":list(reasons), "g2_op":n,
           "collected_flux_pulsed_s":n, "collected_flux_delivered_s":n,
           "optical_pass":False, "hardware_qualified":False,
           "rti_qualified":False, "device_pass":False, "rti_device_pass":False}
    row.update(extra); return row


def _params(mapping, cls):
    if not isinstance(mapping, dict): raise ValueError("nanowire card blocks must be mappings")
    return cls(**mapping)


def _validate(d):
    n = d.nitride
    for k in ("nanowire", "surface", "photonics", "wire_thermal", "injector", "dot"):
        if k not in n: raise ValueError("nitride.%s is required for nanowire cards" % k)
    nw, dot = n["nanowire"], n["dot"]
    family = nw.get("family")
    if family not in ("horizontal_as_built", "vertical_photonic"): raise ValueError("unknown nanowire family")
    if d.drive.mode != "EL-transport" or d.drive.diode.get("preset") != "nitride-nanowire": raise ValueError("nanowire requires EL-transport nitride-nanowire diode")
    if d.ret.mode != "nitride_confinement" or d.cavity.enabled or d.emission.type != "nanowire": raise ValueError("nanowire requires nitride_confinement, cavity disabled, emission nanowire")
    if d.ret.tau_cap_scales_with_density: raise ValueError("nanowire forbids density-scaled capture")
    if d.drive.loading_model != "auto": raise ValueError("nanowire uses pulse-counting, not planar loading model")
    if nw.get("strain_bound") not in ("relaxed", "unrelaxed"): raise ValueError("nanowire strain_bound required")
    sf = dot.get("strain_fraction")
    expected = 0.0 if nw["strain_bound"] == "relaxed" else 1.0
    if sf is not None and float(sf) != expected: raise ValueError("nitride.dot.strain_fraction contradicts strain_bound")
    core = float(nw["core_radius_nm"]); disc = float(dot.get("radius_nm", core))
    if family == "horizontal_as_built" and disc != core: raise ValueError("horizontal disc radius must equal core radius")
    if family == "vertical_photonic" and disc >= core: raise ValueError("vertical photonic disc radius must be less than core radius")
    if n.get("set_params", {}).get("radius_nm", disc) != disc or "C_sigma_F" in n.get("set_params", {}): raise ValueError("SET radius must be the disc radius; C_sigma override forbidden")
    if nw.get("outer_radius_nm", core) < core: raise ValueError("outer radius below core")
    return family, core, disc


def _one(d, T_hs):
    family, core, disc = _validate(d); n=d.nitride; diode_raw=dict(d.drive.diode); diode_raw.pop("preset",None); diode_raw.pop("tau_pulse_ns",None)
    diode = NitrideWireDiode(core_radius_nm=core, conducting_radius_nm=float(diode_raw.pop("conducting_radius_nm",core)),
        barrier_left_nm=float(n["nanowire"].get("barrier_left_nm",15.0)), barrier_right_nm=float(n["nanowire"].get("barrier_right_nm",15.0)),
        d_active_nm=float(n["dot"].get("height_nm",2.0)), x_in=float(n["dot"].get("x_in",.4)), **diode_raw)
    thermal=n["wire_thermal"]; rep=float(d.drive.rep_rate_hz); pulse=float(d.drive.diode.get("tau_pulse_ns",.1)); duty=float(d.drive.duty)
    if rep <= 0 or pulse <= 0 or pulse*rep*1e-9 > 1: return _invalid(T_hs,["invalid pulse/rep-rate"])
    op=wire_operating_point(diode,I_uA=d.drive.I_uA,T_hs_K=T_hs,duty=duty,Rth_K_W=float(thermal["Rth_K_W"]),eta_total=1.0,h_nu_eV=2.0)
    if not op["valid"]: return _invalid(T_hs,op["reasons"], thermal_iterations=op["iterations"])
    tj=op["T_j_K"]
    sys=NitrideNanowireSystem(height_nm=float(n["dot"].get("height_nm",2.0)),core_radius_nm=core,outer_radius_nm=float(n["nanowire"].get("outer_radius_nm",core)),disc_radius_nm=(None if family=="horizontal_as_built" else disc),x_in=float(n["dot"].get("x_in",.4)),strain_bound=n["nanowire"]["strain_bound"],screening_fraction=float(n["dot"].get("screening_fraction",0.0)),external_field_kVcm=0.0)
    lv=levels(sys,tj)
    if not lv.valid: return _invalid(T_hs,lv.invalid_reasons,T_j=tj)
    surf=surface_rates(_params(n["surface"],NitrideNanowireSurfaceParams),core_radius_nm=core,T_K=tj)
    rr=rates(lv,tj,tau_rad0_ns=float(n["dot"].get("tau_rad0_ns",1.0)),tau_cap_ps=float(n["dot"].get("tau_cap_ps",10.0)),reservoir_length_nm=15.0,k_nr_ns=surf["k_surface_X_ns"])
    ph=response(_params({**n["photonics"],"family":family},NitrideNanowirePhotonicsParams),lambda_nm=lv.lambda_nm,outer_radius_nm=sys.outer_radius_nm,gamma_X0_ns=rr["gamma_X0_ns"],gamma_XX0_ns=rr["gamma_XX0_ns"])
    inj=evaluate_injection(diode,I_uA=d.drive.I_uA,T_K=tj,tau_pulse_ns=pulse,E_X_eV=lv.E_X_eV,reservoir_energy_eV=lv.E_X_eV+0.05,barrier_e_eV=lv.dE_e_meV/1000,barrier_h_eV=lv.dE_h_meV/1000,surface_reservoir_ns=surf["k_surface_reservoir_ns"],tau_cap_ps=float(n["dot"].get("tau_cap_ps",10.0)),S_dot=lv.overlap_sq,w_meV=float(d.dot.gamma300),eta_rad_matrix=1.0,eta_total=ph["eta_collection_X"])
    if not inj["valid"]: return _invalid(T_hs,inj["reasons"],T_j=tj)
    gamma_x,gamma_xx=ph["gamma_X_ns"],ph["gamma_XX_ns"]; kx,kxx=rr["k_X_ns"],rr["k_XX_ns"]
    tx=ph["eta_collection_X"]; txx=ph["eta_collection_XX"]
    period=1e9/rep; gate=float(d.drive.gate_ns if d.drive.gate_ns is not None else period)
    regime=n.get("cycle_loading", "rectangular")
    if regime == "deterministic_pair": count=pulse_counting.deterministic_cycle_g2(gamma_x,gamma_xx,kx,kxx,tx,txx,period,eta_load=1.0,gate_ns=gate,split=True)
    else: count=pulse_counting.pulse_g2(inj["r_captured_s"]*1e-9,gamma_x,gamma_xx,kx,kxx,tx,txx,pulse,period-pulse,pump_ratio=d.drive.cw_pump_ratio,gate_ns=gate,split=True)
    b_res=float(d.drive.b_res); mx=count.get("mean_counts_x",_nan()); rho=(tx*mx)/(tx*mx+b_res*mx) if mx>0 else _nan(); g2=count["g2"]/(rho*rho) if rho>0 else _nan()
    flux=rep*max(0.0,count.get("mean_counts",0.0))
    setp=n.get("set_params",{}); st=set_feasibility(tj,radius_nm=disc,eps_r=9.5,R_T_ohm=float(setp.get("R_T_ohm",1e6)),ec_margin=10,f_cycle_Hz=rep)
    ip=replace(_params(n["injector"],NitrideNanowireInjectorParams),reservoir_state_count_e=rr["reservoir_state_count_e"],reservoir_state_count_h=rr["reservoir_state_count_h"])
    rti=injector_feasibility(ip,T_K=tj,rep_rate_hz=rep,loading_window_ns=pulse,electron_level_eV=-lv.dE_e_meV/1000,hole_level_eV=-lv.dE_h_meV/1000,electron_spacing_meV=lv.sp_split_e_meV,hole_spacing_meV=lv.sp_split_h_meV,second_pair_addition_meV=st["E_C_meV"],available_pair_rate_Hz=inj["r_captured_s"],field_kVcm=inj["depletion_field_kVcm"],gate_ns=(gate if regime=="deterministic_pair" else pulse))
    delivery=pulse_delivery(diode,V_j=inj["V_j"],T_K=tj,tau_pulse_ns=pulse,rep_rate_hz=rep)
    one=count.get("one_pair_valid",False); supply=inj["r_captured_s"]>=rep
    valid=bool(count.get("converged",False) and math.isfinite(g2)); optical=valid and g2<.5 and flux>=1000 and (regime!="deterministic_pair" or (one and supply))
    row={"platform":"ingan_gan_nanowire","family":family,"T_hs":float(T_hs),"T_j":tj,"cycle_loading":regime,"valid":valid,"invalid_reasons":[],"g2_op":g2,"g2_dot":count["g2"],"counting_converged":count.get("converged",False),"collected_flux_pulsed_s":flux,"collected_flux_delivered_s":flux*delivery["delivered_step_fraction"],"one_pair_valid":one,"pair_supply_possible":supply,"optical_pass":optical,"set_feasible":st["feasible"],"hardware_qualified":optical and regime=="deterministic_pair" and st["feasible"],"rti_feasible":rti["rti_feasible"],"rti_qualified":optical and regime=="deterministic_pair" and rti["rti_feasible"],"device_pass":optical and (st["feasible"] if regime=="deterministic_pair" else True),"rti_device_pass":optical and regime=="deterministic_pair" and rti["rti_feasible"],"E_X_eV":lv.E_X_eV,"lambda_nm":lv.lambda_nm,"overlap_sq":lv.overlap_sq,"k_X_ns":kx,"k_XX_ns":kxx,"gamma_X0_ns":rr["gamma_X0_ns"],"gamma_X_ns":gamma_x,"tau_rad_bare_ns":1/rr["gamma_X0_ns"],"tau_rad_photonic_ns":1/gamma_x,"tau_total_X_ns":1/(gamma_x+kx),"set_E_C_meV":st["E_C_meV"],"set_EC_over_kT":st["EC_over_kT"],"tau_RC_ns":delivery["tau_RC_s"]*1e9,"delivered_step_fraction":delivery["delivered_step_fraction"],"pulse_delivery_feasible":delivery["pulse_delivery_feasible"],"bound_role":"headline" if sys.strain_bound=="relaxed" else "conservative","bound_reversal":False,"screening_fraction":sys.screening_fraction,"core_radius_nm":core,"outer_radius_nm":sys.outer_radius_nm,"conducting_radius_nm":diode.conducting_radius_nm,"sidewall_overlap":lv.sidewall_overlap,"f_qfl_dot":inj["f_qfl_dot"],"f_qfl_dot_thermodynamic_limit":inj["f_qfl_dot_thermodynamic_limit"],"flat_band":inj["flat_band"],"depletion_regime":inj["depletion_regime"],"thermal_iterations":op["iterations"],"thermal_converged":op["valid"],"eta_collection_X":tx,"background_s":inj["accepted_background_s"]}
    row.update({k:v for k,v in rti.items() if k.startswith("rti_")}); return row


def evaluate_nanowire(design, T_grid=None):
    Ts=np.asarray(T_grid if T_grid is not None else [design.nitride.get("T_hs_K",300.0)],dtype=float)
    rows=[_one(design,float(t)) for t in Ts]; scalars=rows[-1].copy()
    curves={k:np.asarray([r.get(k,_nan()) for r in rows],dtype=float) for k in ("g2_op","collected_flux_pulsed_s","collected_flux_delivered_s","T_j","gamma_X_ns")}
    return {"curves":curves,"scalars":scalars}
