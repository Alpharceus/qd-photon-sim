"""Independent arithmetic checks for single-wire p-i-n transport."""
import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fsim_core.nitride_nanowire_transport import wire_pin, evaluate_injection, wire_operating_point, pulse_delivery
from fsim_core.transport import Q_SI

checks=[]
def ok(name, value):
    checks.append(bool(value))
    if not value: print("FAIL", name)

d=wire_pin()
area15=math.pi*(15e-7)**2; area125=math.pi*(12.5e-7)**2
ok("30 nm current density", math.isclose(1e-9/area15,141.4710605261292,rel_tol=1e-12))
ok("25 nm current density", math.isclose(1e-9/area125,203.71832715762605,rel_tol=1e-12))
ok("fixed measured resistance", math.isclose(1e-9*d.R_s_ohm,2.38,rel_tol=1e-12))
args=dict(T_K=300.,tau_pulse_ns=.1,E_X_eV=2.0,reservoir_energy_eV=2.1,barrier_e_eV=.2,barrier_h_eV=.1,surface_reservoir_ns=1.,tau_cap_ps=10.,S_dot=1.,w_meV=2.,eta_rad_matrix=.2,eta_total=.01)
for current in (0.,.001,1.):
    r=evaluate_injection(d,I_uA=current,**args)
    parts=sum(r[k] for k in ("r_captured_s","r_matrix_radiative_s","r_matrix_nonradiative_s","r_surface_reservoir_s","r_leakage_s","r_other_declared_loss_s"))
    ok("accounting", abs(r["r_supply_s"]-parts)<=max(1.,r["r_supply_s"])*1e-10)
    ok("single pair supply", math.isclose(r["r_supply_s"],current*1e-6/Q_SI,rel_tol=1e-12,abs_tol=1e-6))
r=evaluate_injection(d,I_uA=1.,**args)
ok("surface competes", 0.<r["f_capture"]<1.)
ok("raw accepted distinct", r["raw_background_radiative_s"]>=r["accepted_background_s"]>=0.)
low_res=evaluate_injection(d,I_uA=1.,**dict(args,reservoir_energy_eV=3.0))
ok("reservoir QFL is separate", low_res["f_qfl_background"] < low_res["f_qfl_dot"] and low_res["r_matrix_radiative_s"] < r["r_matrix_radiative_s"])
r_surface=evaluate_injection(d,I_uA=1.,**dict(args,eta_rad_matrix=0.))
ok("surface-only no photons", r_surface["r_matrix_radiative_s"]==0. and r_surface["accepted_background_s"]==0.)
op=wire_operating_point(d,I_uA=.001,T_hs_K=300.,duty=.1,Rth_K_W=1e7,eta_total=.01,h_nu_eV=2.)
I=1e-9; independent=I*op["V_j"]+I*I*d.R_s_ohm-I*.01*2.
ok("heat formula", math.isclose(op["P_on_W"],independent,rel_tol=1e-9))
ok("duty heat", math.isclose(op["P_average_W"],.1*op["P_on_W"],rel_tol=1e-12))
zero=wire_operating_point(d,I_uA=0.,T_hs_K=230.,duty=1.,Rth_K_W=1e9,eta_total=.01,h_nu_eV=2.)
ok("zero heat", zero["T_j_K"]==230.)
p0=pulse_delivery(d,V_j=5.,T_K=300.,tau_pulse_ns=.1,rep_rate_hz=80e6,C_parasitic_F=0.)
p1=pulse_delivery(d,V_j=5.,T_K=300.,tau_pulse_ns=.1,rep_rate_hz=200e6,C_parasitic_F=1e-16)
ok("rc positive", p0["tau_RC_s"]>0. and p1["tau_RC_s"]>p0["tau_RC_s"])
expected=-math.expm1(-.1e-9/p1["tau_RC_s"])
ok("rc exponential", math.isclose(p1["delivered_step_fraction"],expected,rel_tol=1e-12))
ok("diagnostic does not mutate diode", d.R_s_ohm==2.38e9)
print(f"{sum(checks)}/{len(checks)} nitride nanowire transport checks passed")
raise SystemExit(0 if all(checks) else 1)
