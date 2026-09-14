"""Single-wire axial p-i-n transport; no aperture or ensemble partition."""
from __future__ import annotations
from dataclasses import dataclass
import math
from .nitride_materials import KB_EV, binary, bandgap
from .nitride_transport import NitrideDiode
from .transport import Q_SI, qfl_suppression, xi_window, effective_dos
from .thermal import Layer, Stack, t_junction

_PROVENANCE={
 "R_s_ohm":"[V] Deshpande et al., Nat. Commun. 4, 1675 (2013), p.4: 2.38 Gohm measured wire resistance.",
 "designed_R_s_ohm":"[A] 1e6 ohm designed-contact sensitivity (horizontal_designed_contact/vertical_designed_contact presets); not a resistivity-derived geometry law, and R_s_ohm stays fixed at its preset value in geometry sweeps unless this contact model is selected.",
 "geometry":"[V] Deshpande et al., Nat. Commun. 4, 1675 (2013), pp.2-4: disc/spacer/wire geometry; [A] default transfer.",
 "doping":"[V] Deshpande et al., Nat. Commun. 4, 1675 (2013): n-GaN reservoir Si-doped N_D=3e18 cm^-3, p-GaN reservoir Mg-doped N_A=5e17 cm^-3 (NitrideWireDiode defaults; DIRECTIVE ROUND H2 -- previously 1e18/1e17 generic placeholders inherited from the planar diode's own defaults, not this wire's doping).",
 "gan_kernel":"[V/DR H2] The wire junction is a GaN p-n homojunction with a 2 nm InGaN disc as a capture target BELOW the GaN quasi-Fermi separation, not an InGaN-alloy junction: the diode kernel (_GaNJunctionKernel below) is built with x_in=0.0 for j0/V_bi/ideality regardless of NitrideWireDiode.x_in, which remains only a descriptive label of the disc/capture-target composition consumed elsewhere (E_X_eV is supplied explicitly to evaluate_injection, not read off diode.x_in). Root-cause fix for the coherence review's H2/H3: with the InGaN-alloy kernel, V_j at 1 nA was 2.06 V, giving f_qfl_dot ~1e-13 for the measured 2.84 eV dot at the exact bias where Deshpande observed bright single-dot EL (g2 0.30) -- physically incoherent. The GaN kernel gives V_j ~2.6-3.1 V (consistent with the paper's 2-20 V terminal range minus the 2.38 V IR drop [V]) and f_qfl_dot saturates near 1, as expected for a device known to emit at that bias.",
 "tau_SRH_wire_ns":"[DR] NitrideWireDiode.tau_SRH_ns default lowered to 0.01 ns (10 ps) from the planar diode's generic 1.0 ns [E] estimate: no wire-specific SRH-lifetime measurement exists in the digest, so this value is back-solved so the GaN-kernel dark current (H2) reproduces the empirically inferred 2.6-3.1 V junction-voltage window (Deshpande 2013 terminal 2-20 V range minus the 2.38 V IR drop [V]) at 1 nA/300 K; a nanowire's much larger surface-to-volume ratio and dislocation density than a planar layer plausibly shortens SRH lifetime, but this specific number is a fit to the voltage window, not an independent measurement.",
 "capture":"[A] Explicit capture/reservoir competition, not a fitted lifetime.",
 "thermal":"[A] constant-Rth adapter via thermal.t_junction: WIRE_RTH_PRESETS_K_W['horizontal']=1e9 K/W envelope baseline, ['vertical']=1e7 K/W design baseline (one decade sensitivity each way); ['deshpande_fig4_replay']=2.8e9 K/W [DR, re-fit to both Deshpande 2013 Fig.4 heating rises under the H2 GaN-kernel V_j] is a separate replay anchor, not a third baseline. [E, M8] Sign of the 10 K -> 300 K transfer error: this Rth is anchored at the Fig.4 replay's 10 K bath and reused unchanged up to 230-300 K, but GaN's own thermal conductivity FALLS with T above its phonon-scattering peak while amorphous SiO2's conductivity RISES with T; for a horizontal wire whose dominant heat path is lateral into the SiO2 substrate rather than axial conduction along the wire itself, the SiO2-limited path therefore likely has a LOWER Rth at 300 K than at 10 K, so reusing the 10 K-anchored value at 230-300 K plausibly OVERESTIMATES the true junction heating there; no independent 300 K Rth measurement exists in the digest to correct this.",
 "eta_total":"[A] emitted optical power efficiency, not collection efficiency.",
 "pulse":"[A] 90% charge and 10% residual are diagnostic thresholds; loading remains commanded, not delivered.",
 "reservoir_qfl":"[A] f_qfl_background (reservoir-energy occupancy) is reported for diagnostic completeness only in this round; it does NOT suppress the reservoir radiative/nonradiative/surface/other split, which is partitioned by capture-competition rate ratios alone -- f_qfl_dot (source QFL) is the channel that actually suppresses r_captured_s.",
 "qfl_pair":"[DR, H2] f_qfl_dot is evaluated at the diode's actual, current-limited junction voltage V_j -- the DELIVERED quasi-Fermi separation. f_qfl_dot_thermodynamic_limit is the identical Boltzmann-tail expression evaluated at V_bi instead -- the maximum quasi-Fermi separation the junction can thermodynamically support. This is a delivered-vs-thermodynamic-ceiling pair, not two independent physical channels; only f_qfl_dot suppresses r_captured_s (see 'reservoir_qfl').",
 "naming":"[L1] surface_reservoir_ns (and its alias keyword k_surface_reservoir_per_ns) is a RATE in ns^-1, not a lifetime, despite the _ns suffix reading like a time constant next to tau_pulse_ns/tau_matrix_ns/tau_cap_ps (which ARE times) elsewhere in this module's keyword list -- a naming hazard carried over from piece 4's k_surface_reservoir_ns convention.",
 "background":"[A] raw_background_radiative_s equals r_matrix_radiative_s by definition -- the QFL-suppressed matrix radiative rate before the spectral window; accepted_background_s = raw * background_window_fraction is the distinct, window-filtered quantity."}
WIRE_PRESETS={"deshpande_2013_30nm":{"core_radius_nm":15.,"conducting_radius_nm":15.},"horizontal_designed_contact":{"R_s_ohm":1e6},"vertical_designed_contact":{"R_s_ohm":1e6}}
WIRE_RTH_PRESETS_K_W={"horizontal":1e9,"vertical":1e7,"deshpande_fig4_replay":2.8e9} # [DR, re-fit to both Fig.4 rises under the H2 GaN-kernel V_j] -- see _PROVENANCE['thermal'] for the value's derivation and the M8 sign-of-transfer-error note (previously 3.1e9 under the pre-H2 InGaN-alloy kernel).

def _pos(n,v,z=False):
 v=float(v)
 if not math.isfinite(v) or (v<0 if z else v<=0): raise ValueError(n+" must be finite and "+("nonnegative" if z else "positive"))
 return v

class _GaNJunctionKernel(NitrideDiode):
 """GaN p-n homojunction kernel (H2): the wire junction, not the InGaN
 disc alloy, sets j0/V_bi/ideality -- NitrideWireDiode._kernel() always
 constructs this with x_in=0.0 (see _PROVENANCE['gan_kernel']).

 vbi() and vj_of_j() are overridden here in LOG SPACE (the same
 technique transport.ln_n_i uses, and for the same reason) so that
 ln(ni) stays finite for any T > 0: NitrideDiode's own literal
 ni/ni**2/j0 (inherited unmodified, and still reachable through this
 subclass's _j0()/_ni()) underflow to exactly 0.0 in double precision
 below roughly 25-60 K -- see the module-level _kernel() docstring
 below for the M9 consequence. depletion() and its (V_bi - V_j)
 quadratic are inherited from NitrideDiode UNCHANGED; only the vbi()
 it calls internally is replaced.
 """
 def _ln_ni_gan(self,T):
  # [DR] ln(ni) [ln cm^-3] for GaN, computed directly in log space (same
  # DOS-mass geometric mean as nitride_transport._dos/_nv) so that ni
  # itself never has to be formed as a literal, possibly-underflowing
  # number.
  T=self._T(T); gan=binary("GaN")
  Nc=effective_dos((gan.me_xy*gan.me_xy*gan.me_z)**(1/3),T); Nv=effective_dos((gan.mh_xy*gan.mh_xy*gan.mh_z)**(1/3),T)
  return .5*(math.log(Nc)+math.log(Nv))-bandgap(gan,T)/(2.*KB_EV*T)
 def vbi(self,T=None):
  # [DR, M9] Degenerate/frozen-carrier V_bi = (kT/q)*[ln(N_A*N_D) - 2 ln(ni)]
  # -- the textbook ln(N_A N_D / ni^2) form, expanded so ni is never
  # squared (or even formed) as a literal -- valid at any T > 0, unlike
  # the ni**2 formula it replaces (NitrideDiode.vbi, unmodified, still
  # available to any caller wanting the original planar behaviour).
  T=self._T(T); return KB_EV*T*(math.log(self.N_A)+math.log(self.N_D)-2.*self._ln_ni_gan(T))
 def _ln_j0(self,T):
  # [DR] SRH j0 in log space: with x_in=0.0 the "active" and "spacer"
  # regions are the same GaN, so the two-region ni_a*d_a+ni_b*d_b sum
  # NitrideDiode._j0 performs collapses to a single ni*d_i term.
  T=self._T(T); return math.log(Q_SI)+self._ln_ni_gan(T)+math.log(self.d_i_nm*1e-7)-math.log(2.*self.tau_SRH_ns*1e-9)
 def vj_of_j(self,J,T=None):
  T=self._T(T)
  if J<=0.: return 0.
  x=math.log(J)-self._ln_j0(T)
  sp=x+math.log1p(math.exp(-x)) if x>0 else math.log1p(math.exp(x)) # [DR] numerically stable log1p(exp(x)) ("softplus"): replaces n*kT*log1p(J/j0), which would need the literal (underflowing) j0.
  return self.n_ideality*KB_EV*T*sp

@dataclass(frozen=True)
class NitrideWireDiode:
 N_A:float=5e17; N_D:float=3e18; n_ideality:float=2.; tau_SRH_ns:float=.01; eps_r:float=10.28; T:float=300.
 core_radius_nm:float=12.5; conducting_radius_nm:float=12.5; barrier_left_nm:float=15.; barrier_right_nm:float=15.; d_active_nm:float=2.; x_in:float=.40
 R_s_ohm:float=2.38e9; f_Rs_local:float=1.; tau_matrix_ns:float=1.
 def __post_init__(self):
  for n in ("N_A","N_D","n_ideality","tau_SRH_ns","eps_r","T","core_radius_nm","conducting_radius_nm","barrier_left_nm","barrier_right_nm","d_active_nm","R_s_ohm","tau_matrix_ns"): _pos(n,getattr(self,n))
  if not 0<=self.x_in<=1 or not 0<=self.f_Rs_local<=1: raise ValueError("x_in or f_Rs_local outside physical range")
  if self.conducting_radius_nm>self.core_radius_nm: raise ValueError("conducting_radius_nm must not exceed core_radius_nm")
 @property
 def d_i_nm(self): return self.barrier_left_nm+self.d_active_nm+self.barrier_right_nm # [DR] axial sum
 @property
 def area_cm2(self): return math.pi*(self.conducting_radius_nm*1e-7)**2 # [DR]
 def _kernel(self):
  # [V/DR H2] x_in is FORCED to 0.0 here regardless of self.x_in -- the
  # wire junction is GaN p-n, not the InGaN disc alloy. self.x_in is
  # retained only as the disc/capture-target composition label for
  # other callers (see _PROVENANCE['gan_kernel']).
  return _GaNJunctionKernel(N_A=self.N_A,N_D=self.N_D,d_i_nm=self.d_i_nm,d_active_nm=self.d_active_nm,area_um2=self.area_cm2*1e8,R_s_ohm=self.R_s_ohm,n_ideality=self.n_ideality,tau_SRH_ns=self.tau_SRH_ns,eps_r=self.eps_r,T=self.T,x_in=0.,wl_thickness_nm=self.d_active_nm,f_Rs_local=self.f_Rs_local)

def wire_pin(**overrides):
 if "area_cm2" in overrides or "area_um2" in overrides: raise TypeError("wire area is derived from conducting_radius_nm")
 p=overrides.pop("preset",None)
 if p is not None:
  if p not in WIRE_PRESETS: raise ValueError("unknown wire preset")
  x=dict(WIRE_PRESETS[p]); x.update(overrides); overrides=x
 return NitrideWireDiode(**overrides)

def _kernel(d,I,T):
 """Returns (vt,v,dep,reason); reason is None on success.

 H2/M9 fix: d._kernel() builds a _GaNJunctionKernel (GaN p-n
 homojunction, x_in=0.0) whose vbi()/vj_of_j() are LOG-SPACE
 reformulations of NitrideDiode's ni/ni**2/j0 formulas (see
 _GaNJunctionKernel's docstring above) -- ln(ni) stays finite for any
 T > 0, so V_j and V_bi no longer underflow to exactly 0.0 (and the
 division/log(0) errors that used to follow) the way the base class's
 literal ni/ni**2/j0 path does below roughly 25-60 K. The previous
 pre-check `k._j0(T)==0.` (which forced every T below that floor
 invalid, before M9) is REMOVED: the 10 K Deshpande replay now returns
 a valid row.
 'kernel_invalid_below_T' (an exception reaching here) and
 'kernel_error' (a non-finite result reaching here without raising, e.g.
 an absurd R_s_ohm large enough to overflow V_terminal) now cover only
 genuinely pathological inputs, not routine low-temperature operation.
 """
 k=d._kernel()
 try:
  vt,v=k.v_of_i(I,T); dep=k.depletion(v,T)
 except (ArithmeticError,ValueError,OverflowError): return None,None,None,"kernel_invalid_below_T"
 if all(math.isfinite(x) for x in (vt,v,dep.F_kVcm,dep.C_dep_pF)): return vt,v,dep,None
 return None,None,None,"kernel_error"

def _invalid(d,I,reason):
 n=float("nan")
 return {"valid":False,"reasons":[reason],"provenance":dict(_PROVENANCE),
  "area_cm2":d.area_cm2,"J_A_cm2":I/d.area_cm2,"V_j":n,"V_bi":n,"V_terminal":n,"depletion_field_kVcm":n,"C_dep_F":n,
  "eta_inj":n,"r_e_leak_ratio":n,"r_h_leak_ratio":n,
  "f_capture":n,"f_qfl_dot":n,"f_qfl_dot_thermodynamic_limit":n,"f_qfl_background":n,
  "r_supply_s":I/Q_SI,"r_captured_s":n,"r_matrix_radiative_s":n,"r_matrix_nonradiative_s":n,"r_surface_reservoir_s":n,"r_leakage_s":n,"r_other_declared_loss_s":n,
  "raw_background_radiative_s":n,"accepted_background_s":n,"background_window_fraction":n,
  "mu":n,"ideal_min_pair_current_A":n,"ideal_min_pair_current_uA":n,
  "power_on_W":n,"accounting_residual_s":n}

def evaluate_injection(diode,*,I_uA,T_K,tau_pulse_ns,E_X_eV,reservoir_energy_eV,barrier_e_eV,barrier_h_eV,surface_reservoir_ns=None,tau_cap_ps,S_dot,w_meV,eta_rad_matrix,eta_total,k_surface_reservoir_per_ns=None):
 """Evaluate one wire operating point. surface_reservoir_ns is a RATE in
 ns^-1 despite the _ns suffix (piece 4's k_surface_reservoir_ns
 convention -- see _PROVENANCE['naming'], L1); k_surface_reservoir_per_ns
 is an unambiguous alias keyword for the same quantity. Pass exactly one
 of the two, or both if they agree -- passing both with different values
 raises. f_qfl_dot (delivered, at V_j) and f_qfl_dot_thermodynamic_limit
 (ceiling, at V_bi) are a paired pair of columns, not two channels; see
 _PROVENANCE['qfl_pair'].
 """
 if not isinstance(diode,NitrideWireDiode): raise TypeError("diode must be NitrideWireDiode")
 if surface_reservoir_ns is None and k_surface_reservoir_per_ns is None: raise TypeError("surface_reservoir_ns (or its alias k_surface_reservoir_per_ns) is required")
 if surface_reservoir_ns is not None and k_surface_reservoir_per_ns is not None and float(surface_reservoir_ns)!=float(k_surface_reservoir_per_ns): raise ValueError("surface_reservoir_ns and k_surface_reservoir_per_ns disagree")
 surface_reservoir_ns=surface_reservoir_ns if surface_reservoir_ns is not None else k_surface_reservoir_per_ns
 Iu=_pos("I_uA",I_uA,True); T=_pos("T_K",T_K); tp=_pos("tau_pulse_ns",tau_pulse_ns,True)
 for n,x in (("E_X_eV",E_X_eV),("reservoir_energy_eV",reservoir_energy_eV),("barrier_e_eV",barrier_e_eV),("barrier_h_eV",barrier_h_eV),("surface_reservoir_ns",surface_reservoir_ns),("tau_cap_ps",tau_cap_ps),("S_dot",S_dot),("w_meV",w_meV),("eta_rad_matrix",eta_rad_matrix),("eta_total",eta_total)): _pos(n,x,True)
 if any(float(x)>1 for x in (S_dot,eta_rad_matrix,eta_total)): raise ValueError("efficiencies must be in [0, 1]")
 I=Iu*1e-6; vt,v,dep,reason=_kernel(diode,I,T)
 if reason: return _invalid(diode,I,reason)
 re=math.exp(-float(barrier_e_eV)/(KB_EV*T)); rh=math.exp(-float(barrier_h_eV)/(KB_EV*T)); eta=1/(1+re+rh)
 supply=I/Q_SI; leak=supply*(1-eta); usable=supply-leak
 fd=qfl_suppression(float(E_X_eV),v,KB_EV*T); fd_lim=qfl_suppression(float(E_X_eV),dep.V_bi,KB_EV*T) # [DR, H2] delivered (V_j) vs thermodynamic-ceiling (V_bi) QFL pair; see _PROVENANCE['qfl_pair'].
 fbg=qfl_suppression(float(reservoir_energy_eV),v,KB_EV*T)
 kc=0. if float(tau_cap_ps)<=0 else 1000/float(tau_cap_ps) # [DR] tau_cap_ps<=0 guarded to k_cap=0 (no dot-capture channel) instead of a ZeroDivisionError; conservative -- does not assume instantaneous, certain capture at the limit.
 kr=float(eta_rad_matrix)/diode.tau_matrix_ns; knr=(1-float(eta_rad_matrix))/diode.tau_matrix_ns
 ks=float(surface_reservoir_ns) # [DR] surface_reservoir_ns is already a RATE in ns^-1 (piece 4's k_surface_reservoir_ns convention), not a lifetime to invert
 fc=kc/(kc+kr+knr+ks)
 cap=usable*fc*float(S_dot)*fd; res=usable-cap; kres=kr+knr+ks
 # kres=kr+knr+ks is always >0: kr+knr=1/tau_matrix_ns alone is >0
 # (tau_matrix_ns is validated positive), so the reservoir split never
 # reaches kres==0 and r_other_declared_loss_s stays 0 by construction,
 # not by an unreachable branch.
 rad=res*kr/kres; nr=res*knr/kres; surf=res*ks/kres; other=0.
 det=(float(reservoir_energy_eV)-float(E_X_eV))*1e3; xi=float(xi_window(float(w_meV),det,KB_EV*T*1e3)) if det>=0 else 0.; residual=supply-(leak+cap+rad+nr+surf+other)
 ideal=float("inf") if tp==0 else Q_SI/(tp*1e-9); power=I*v+diode.f_Rs_local*I*I*diode.R_s_ohm-I*float(eta_total)*float(E_X_eV)
 return {"valid":math.isfinite(residual),"reasons":[],"provenance":dict(_PROVENANCE),"area_cm2":diode.area_cm2,"J_A_cm2":I/diode.area_cm2,"V_j":v,"V_bi":dep.V_bi,"V_terminal":vt,"depletion_field_kVcm":dep.F_kVcm,"C_dep_F":dep.C_dep_pF*1e-12,"eta_inj":eta,"r_e_leak_ratio":re,"r_h_leak_ratio":rh,"f_capture":fc,"f_qfl_dot":fd,"f_qfl_dot_thermodynamic_limit":fd_lim,"f_qfl_background":fbg,"r_supply_s":supply,"r_captured_s":cap,"r_matrix_radiative_s":rad,"r_matrix_nonradiative_s":nr,"r_surface_reservoir_s":surf,"r_leakage_s":leak,"r_other_declared_loss_s":other,"raw_background_radiative_s":rad,"accepted_background_s":rad*xi,"background_window_fraction":xi,"mu":cap*tp*1e-9,"ideal_min_pair_current_A":ideal,"ideal_min_pair_current_uA":ideal*1e6,"power_on_W":power,"accounting_residual_s":residual}

def wire_operating_point(diode,*,I_uA,T_hs_K,duty,Rth_K_W,eta_total,h_nu_eV):
 if not isinstance(diode,NitrideWireDiode): raise TypeError("diode must be NitrideWireDiode")
 Iu=_pos("I_uA",I_uA,True); th=_pos("T_hs_K",T_hs_K); rth=_pos("Rth_K_W",Rth_K_W,True); _pos("h_nu_eV",h_nu_eV,True)
 if not 0<=float(duty)<=1 or not 0<=float(eta_total)<=1: raise ValueError("duty and eta_total must be in [0, 1]")
 I=Iu*1e-6; vt,v,_,reason=_kernel(diode,I,th)
 if reason:return {"valid":False,"reasons":[reason],"T_j_K":th,"P_on_W":float("nan"),"P_average_W":float("nan"),"V_j":float("nan"),"V_terminal":float("nan"),"iterations":0,"provenance":dict(_PROVENANCE)}
 def out(v,vt,t,po,it,valid=True,reasons=[]):return {"valid":valid,"reasons":reasons,"T_j_K":t,"P_on_W":po,"P_average_W":float(duty)*po,"V_j":v,"V_terminal":vt,"iterations":it,"provenance":dict(_PROVENANCE)}
 p=lambda x:I*x+diode.f_Rs_local*I*I*diode.R_s_ohm-I*float(eta_total)*float(h_nu_eV); po=p(v)
 if I==0 or rth==0 or duty==0:return out(v,vt,th,po,0)
 a=diode.conducting_radius_nm*1e-9; stack=Stack([Layer(rth*math.pi*a*a,1.,0.,False)],1e300,0.); tj=th
 for n in range(1,81):
  qvt,qv,_,qreason=_kernel(diode,I,tj)
  if qreason:return out(v,vt,tj,po,n,False,[qreason])
  vt,v=qvt,qv; po=p(v); new=t_junction(float(duty)*po,a,stack,th)
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
