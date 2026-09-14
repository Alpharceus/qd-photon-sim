"""Single-wire axial p-i-n transport; no aperture or ensemble partition."""
from __future__ import annotations
from dataclasses import dataclass
import math
from .nitride_materials import KB_EV
from .nitride_transport import NitrideDiode
from .transport import Q_SI, qfl_suppression, xi_window
from .thermal import Layer, Stack, t_junction

_PROVENANCE={
 "R_s_ohm":"[V] Deshpande et al., Nat. Commun. 4, 1675 (2013), p.4: 2.38 Gohm measured wire resistance.",
 "geometry":"[V] Deshpande et al., Nat. Commun. 4, 1675 (2013), pp.2-4: disc/spacer/wire geometry; [A] default transfer.",
 "capture":"[A] Explicit capture/reservoir competition, not a fitted lifetime.",
 "thermal":"[A] constant-Rth adapter via thermal.t_junction; [DR] 3e9 K/W replays Deshpande Fig.4 heating.",
 "eta_total":"[A] emitted optical power efficiency, not collection efficiency.",
 "pulse":"[A] 90% charge and 10% residual are diagnostic thresholds; loading remains commanded, not delivered.",
 "reservoir_qfl":"[A] Reservoir QFL occupancy multiplies every reservoir recombination rate and cancels in the competing-rate branching ratio."}
WIRE_PRESETS={"deshpande_2013_30nm":{"core_radius_nm":15.,"conducting_radius_nm":15.},"horizontal_designed_contact":{"R_s_ohm":1e6},"vertical_designed_contact":{"R_s_ohm":1e6}}
WIRE_RTH_PRESETS_K_W={"horizontal":1e9,"vertical":1e7,"deshpande_fig4_replay":3.1e9} # [DR] approximately 3e9 Fig.4 replay

def _pos(n,v,z=False):
 v=float(v)
 if not math.isfinite(v) or (v<0 if z else v<=0): raise ValueError(n+" must be finite and "+("nonnegative" if z else "positive"))
 return v

@dataclass(frozen=True)
class NitrideWireDiode:
 N_A:float=1e17; N_D:float=1e18; n_ideality:float=2.; tau_SRH_ns:float=1.; eps_r:float=10.28; T:float=300.
 core_radius_nm:float=12.5; conducting_radius_nm:float=12.5; barrier_left_nm:float=15.; barrier_right_nm:float=15.; d_active_nm:float=2.; x_in:float=.40
 R_s_ohm:float=2.38e9; f_Rs_local:float=1.; tau_matrix_ns:float=1.; tau_cap_ps:float=10.; reservoir_surface_ns:float=0.
 def __post_init__(self):
  for n in ("N_A","N_D","n_ideality","tau_SRH_ns","eps_r","T","core_radius_nm","conducting_radius_nm","barrier_left_nm","barrier_right_nm","d_active_nm","R_s_ohm","tau_matrix_ns","tau_cap_ps"): _pos(n,getattr(self,n))
  _pos("reservoir_surface_ns",self.reservoir_surface_ns,True)
  if not 0<=self.x_in<=1 or not 0<=self.f_Rs_local<=1: raise ValueError("x_in or f_Rs_local outside physical range")
  if self.conducting_radius_nm>self.core_radius_nm: raise ValueError("conducting_radius_nm must not exceed core_radius_nm")
 @property
 def d_i_nm(self): return self.barrier_left_nm+self.d_active_nm+self.barrier_right_nm # [DR] axial sum
 @property
 def area_cm2(self): return math.pi*(self.conducting_radius_nm*1e-7)**2 # [DR]
 def _kernel(self): return NitrideDiode(N_A=self.N_A,N_D=self.N_D,d_i_nm=self.d_i_nm,d_active_nm=self.d_active_nm,area_um2=self.area_cm2*1e8,R_s_ohm=self.R_s_ohm,n_ideality=self.n_ideality,tau_SRH_ns=self.tau_SRH_ns,eps_r=self.eps_r,T=self.T,x_in=self.x_in,wl_thickness_nm=self.d_active_nm,f_Rs_local=self.f_Rs_local)

def wire_pin(**overrides):
 if "area_cm2" in overrides or "area_um2" in overrides: raise TypeError("wire area is derived from conducting_radius_nm")
 p=overrides.pop("preset",None)
 if p is not None:
  if p not in WIRE_PRESETS: raise ValueError("unknown wire preset")
  x=dict(WIRE_PRESETS[p]); x.update(overrides); overrides=x
 return NitrideWireDiode(**overrides)

def _kernel(d,I,T):
 try:
  k=d._kernel()
  if k._j0(T)==0.: return None
  vt,v=k.v_of_i(I,T); dep=k.depletion(v,T)
  return (vt,v,dep) if all(math.isfinite(x) for x in (vt,v,dep.F_kVcm,dep.C_dep_pF)) else None
 except (ArithmeticError,ValueError,OverflowError): return None

def _invalid(d,I):
 return {"valid":False,"reasons":["kernel_invalid_below_T"],"provenance":dict(_PROVENANCE),"area_cm2":d.area_cm2,"J_A_cm2":I/d.area_cm2,"V_j":float("nan"),"V_terminal":float("nan"),"depletion_field_kVcm":float("nan"),"C_dep_F":float("nan"),"r_supply_s":I/Q_SI,"accounting_residual_s":float("nan")}

def evaluate_injection(diode,*,I_uA,T_K,tau_pulse_ns,E_X_eV,reservoir_energy_eV,barrier_e_eV,barrier_h_eV,surface_reservoir_ns,tau_cap_ps,S_dot,w_meV,eta_rad_matrix,eta_total):
 if not isinstance(diode,NitrideWireDiode): raise TypeError("diode must be NitrideWireDiode")
 Iu=_pos("I_uA",I_uA,True); T=_pos("T_K",T_K); tp=_pos("tau_pulse_ns",tau_pulse_ns,True)
 for n,x in (("E_X_eV",E_X_eV),("reservoir_energy_eV",reservoir_energy_eV),("barrier_e_eV",barrier_e_eV),("barrier_h_eV",barrier_h_eV),("surface_reservoir_ns",surface_reservoir_ns),("tau_cap_ps",tau_cap_ps),("S_dot",S_dot),("w_meV",w_meV),("eta_rad_matrix",eta_rad_matrix),("eta_total",eta_total)): _pos(n,x,True)
 if any(float(x)>1 for x in (S_dot,eta_rad_matrix,eta_total)): raise ValueError("efficiencies must be in [0, 1]")
 I=Iu*1e-6; kv=_kernel(diode,I,T)
 if kv is None:return _invalid(diode,I)
 vt,v,dep=kv; re=math.exp(-float(barrier_e_eV)/(KB_EV*T)); rh=math.exp(-float(barrier_h_eV)/(KB_EV*T)); eta=1/(1+re+rh)
 supply=I/Q_SI; leak=supply*(1-eta); usable=supply-leak; fd=qfl_suppression(float(E_X_eV),v,KB_EV*T); fbg=qfl_suppression(float(reservoir_energy_eV),v,KB_EV*T)
 kc=1000/float(tau_cap_ps); kr=float(eta_rad_matrix)/diode.tau_matrix_ns; knr=(1-float(eta_rad_matrix))/diode.tau_matrix_ns; ks=float(surface_reservoir_ns); fc=kc/(kc+kr+knr+ks)
 cap=usable*fc*float(S_dot)*fd; res=usable-cap; kres=kr+knr+ks
 rad=res*kr/kres if kres else 0.; nr=res*knr/kres if kres else 0.; surf=res*ks/kres if kres else 0.; other=res if not kres else 0.
 det=(float(reservoir_energy_eV)-float(E_X_eV))*1e3; xi=float(xi_window(float(w_meV),det,KB_EV*T*1e3)) if det>=0 else 0.; residual=supply-(leak+cap+rad+nr+surf+other)
 ideal=float("inf") if tp==0 else Q_SI/(tp*1e-9); power=I*v+diode.f_Rs_local*I*I*diode.R_s_ohm-I*float(eta_total)*float(E_X_eV)
 return {"valid":math.isfinite(residual),"reasons":[],"provenance":dict(_PROVENANCE),"area_cm2":diode.area_cm2,"J_A_cm2":I/diode.area_cm2,"V_j":v,"V_terminal":vt,"depletion_field_kVcm":dep.F_kVcm,"C_dep_F":dep.C_dep_pF*1e-12,"eta_inj":eta,"r_e_leak_ratio":re,"r_h_leak_ratio":rh,"leakage_e_ratio":re,"leakage_h_ratio":rh,"f_capture":fc,"f_qfl_dot":fd,"f_qfl_background":fbg,"r_supply_s":supply,"r_captured_s":cap,"r_matrix_radiative_s":rad,"r_matrix_nonradiative_s":nr,"r_surface_reservoir_s":surf,"r_leakage_s":leak,"r_other_declared_loss_s":other,"raw_background_radiative_s":rad,"accepted_background_s":rad*xi,"background_window_fraction":xi,"mu":cap*tp*1e-9,"mu_commanded_not_delivered":cap*tp*1e-9,"ideal_min_pair_current_A":ideal,"ideal_min_pair_current_uA":ideal*1e6,"power_on_W":power,"accounting_residual_s":residual}

def wire_operating_point(diode,*,I_uA,T_hs_K,duty,Rth_K_W,eta_total,h_nu_eV):
 if not isinstance(diode,NitrideWireDiode): raise TypeError("diode must be NitrideWireDiode")
 Iu=_pos("I_uA",I_uA,True); th=_pos("T_hs_K",T_hs_K); rth=_pos("Rth_K_W",Rth_K_W,True); _pos("h_nu_eV",h_nu_eV,True)
 if not 0<=float(duty)<=1 or not 0<=float(eta_total)<=1: raise ValueError("duty and eta_total must be in [0, 1]")
 I=Iu*1e-6; kv=_kernel(diode,I,th)
 if kv is None:return {"valid":False,"reasons":["kernel_invalid_below_T"],"T_j_K":th,"P_on_W":float("nan"),"P_average_W":float("nan"),"V_j":float("nan"),"V_terminal":float("nan"),"iterations":0,"provenance":dict(_PROVENANCE)}
 def out(v,vt,t,po,it,valid=True,reasons=[]):return {"valid":valid,"reasons":reasons,"T_j_K":t,"P_on_W":po,"P_average_W":float(duty)*po,"V_j":v,"V_terminal":vt,"iterations":it,"provenance":dict(_PROVENANCE)}
 vt,v,_=kv; p=lambda x:I*x+diode.f_Rs_local*I*I*diode.R_s_ohm-I*float(eta_total)*float(h_nu_eV); po=p(v)
 if I==0 or rth==0 or duty==0:return out(v,vt,th,po,0)
 a=diode.conducting_radius_nm*1e-9; stack=Stack([Layer(rth*math.pi*a*a,1.,0.,False)],1e300,0.); tj=th
 for n in range(1,81):
  q=_kernel(diode,I,tj)
  if q is None:return out(v,vt,tj,po,n,False,["kernel_invalid_below_T"])
  vt,v,_=q; po=p(v); new=t_junction(float(duty)*po,a,stack,th)
  if not math.isfinite(new):return out(v,vt,float("inf"),po,n,False,["thermal runaway/nonconvergence"])
  if abs(new-tj)<1e-5:return out(v,vt,new,po,n)
  tj=.5*(tj+new)
 return out(v,vt,tj,po,80,False,["thermal nonconvergence"])

def pulse_delivery(diode,*,V_j,T_K,tau_pulse_ns,rep_rate_hz,C_parasitic_F=0.):
 if not isinstance(diode,NitrideWireDiode):raise TypeError("diode must be NitrideWireDiode")
 v=_pos("V_j",V_j,True); t=_pos("T_K",T_K); width=_pos("tau_pulse_ns",tau_pulse_ns,True)*1e-9; rate=_pos("rep_rate_hz",rep_rate_hz); cp=_pos("C_parasitic_F",C_parasitic_F,True)
 try: dep=diode._kernel().depletion(v,t)
 except (ArithmeticError,ValueError,OverflowError):return {"valid":False,"reasons":["kernel_invalid_below_T"],"provenance":dict(_PROVENANCE)}
 cd=dep.C_dep_pF*1e-12; ct=cd+cp; tau=diode.R_s_ohm*ct
 step=1. if tau==0 else -math.expm1(-width/tau); f3=float("inf") if tau==0 else 1/(2*math.pi*tau); residual=0. if tau==0 else (math.exp(-(1/rate-width)/tau) if width<1/rate else 1.)
 return {"valid":True,"reasons":[],"provenance":dict(_PROVENANCE),"C_dep_F":cd,"C_total_F":ct,"tau_RC_s":tau,"f_3dB_Hz":f3,"delivered_step_fraction":step,"inter_pulse_residual_fraction":residual,"pulse_delivery_feasible":step>=.9 and residual<=.1,"loading_convention":"[A] mu is commanded, not delivered; RC diagnostic only"}
