"""Independent checks for wire transport."""
import math,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from fsim_core.nitride_nanowire_transport import wire_pin,evaluate_injection,wire_operating_point,pulse_delivery,WIRE_RTH_PRESETS_K_W
from fsim_core.transport import Q_SI
C=[]
def ok(n,x):
 C.append(bool(x))
 if not x: print("FAIL",n)
a=dict(T_K=300.,tau_pulse_ns=.1,E_X_eV=2.,reservoir_energy_eV=2.1,barrier_e_eV=.21,barrier_h_eV=.08,surface_reservoir_ns=1.,tau_cap_ps=10.,S_dot=.8,w_meV=2.,eta_rad_matrix=.2,eta_total=.01)
d15=wire_pin(preset="deshpande_2013_30nm");d125=wire_pin(core_radius_nm=20.,conducting_radius_nm=12.5)
for d,target in ((d15,141.4710605261292),(d125,203.71832715762605)):
 r=evaluate_injection(d,I_uA=.001,**a); area=math.pi*(d.conducting_radius_nm*1e-7)**2
 ok("module conducting area",math.isclose(r["area_cm2"],area,rel_tol=1e-13));ok("module J anchor",math.isclose(r["J_A_cm2"],target,rel_tol=.001))
ok("core differs from conducting",d125.area_cm2 != math.pi*(d125.core_radius_nm*1e-7)**2);ok("2.38V IR",math.isclose(1e-9*d15.R_s_ohm,2.38,rel_tol=1e-13))
for x in (.25,.4):
 for T in (230.,300.):
  for u in (0.,.001,1.):
   r=evaluate_injection(wire_pin(x_in=x),I_uA=u,**dict(a,T_K=T));s=sum(r[k] for k in ("r_captured_s","r_matrix_radiative_s","r_matrix_nonradiative_s","r_surface_reservoir_s","r_leakage_s","r_other_declared_loss_s"))
   ok("independent accounting",abs(r["r_supply_s"]-s)<=max(1.,r["r_supply_s"])*1e-10);ok("one serial pair",math.isclose(r["r_supply_s"],u*1e-6/Q_SI,rel_tol=1e-12,abs_tol=1e-6))
r=evaluate_injection(d15,I_uA=1.,**a);ok("surface capture competition",math.isclose(r["f_capture"],100/(100+.2+.8+1),rel_tol=1e-12));ok("dot qfl capture",math.isclose(r["r_captured_s"],(r["r_supply_s"]-r["r_leakage_s"])*r["f_capture"]*.8*r["f_qfl_dot"],rel_tol=1e-12))
swap=evaluate_injection(d15,I_uA=1.,**dict(a,barrier_e_eV=.08,barrier_h_eV=.21));ok("barrier paths observable",r["r_e_leak_ratio"]!=r["r_h_leak_ratio"] and r["r_e_leak_ratio"]==swap["r_h_leak_ratio"])
z=evaluate_injection(d15,I_uA=1.,**dict(a,eta_rad_matrix=0.));ok("surface loss no photon",z["r_matrix_radiative_s"]==z["accepted_background_s"]==0.)
below=evaluate_injection(d15,I_uA=1.,**dict(a,reservoir_energy_eV=1.9));ok("below exciton rejected",below["background_window_fraction"]==0.)
ok("minimum pair current",math.isclose(r["ideal_min_pair_current_A"],Q_SI/1e-10,rel_tol=1e-12));low=evaluate_injection(d15,I_uA=.001,**dict(a,T_K=10.));ok("10K invalid no crash",not low["valid"] and "kernel_invalid_below_T" in low["reasons"])
for u,target in ((.001,15.),(.002,49.)):
 o=wire_operating_point(d15,I_uA=u,T_hs_K=60.,duty=1.,Rth_K_W=WIRE_RTH_PRESETS_K_W["deshpande_fig4_replay"],eta_total=0.,h_nu_eV=2.);I=u*1e-6
 ok("independent heat",math.isclose(o["P_on_W"],I*o["V_j"]+I*I*d15.R_s_ohm,rel_tol=1e-10));ok("Fig4 thermal anchor",math.isclose(o["T_j_K"]-60.,target,rel_tol=.15))
o=wire_operating_point(d15,I_uA=.001,T_hs_K=300.,duty=1.,Rth_K_W=0.,eta_total=.01,h_nu_eV=2.);ok("zero Rth electrical retained",o["T_j_K"]==300. and o["P_on_W"]>0 and o["V_j"]>0)
o=wire_operating_point(d15,I_uA=.001,T_hs_K=300.,duty=.1,Rth_K_W=1e7,eta_total=.01,h_nu_eV=2.);ok("duty not squared average current",math.isclose(o["P_average_W"],.1*o["P_on_W"],rel_tol=1e-12))
for cp,rate in ((0.,80e6),(1e-18,80e6),(1e-16,200e6)):
 p=pulse_delivery(d15,V_j=5.,T_K=300.,tau_pulse_ns=.1,rep_rate_hz=rate,C_parasitic_F=cp);tau=d15.R_s_ohm*(p["C_dep_F"]+cp)
 ok("RC independently summed",math.isclose(p["tau_RC_s"],tau,rel_tol=1e-12));ok("RC exponential independent",math.isclose(p["delivered_step_fraction"],1-math.exp(-1e-10/tau),rel_tol=1e-12))
ok("RC diagnostic",not pulse_delivery(d15,V_j=5.,T_K=300.,tau_pulse_ns=.1,rep_rate_hz=80e6)["pulse_delivery_feasible"])
print(f"{sum(C)}/{len(C)} nitride nanowire transport checks passed")
raise SystemExit(0 if all(C) else 1)
