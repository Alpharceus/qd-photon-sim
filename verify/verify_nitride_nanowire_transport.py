"""Independent checks for wire transport."""
import math,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from fsim_core.nitride_nanowire_transport import wire_pin,evaluate_injection,wire_operating_point,pulse_delivery,WIRE_RTH_PRESETS_K_W,NitrideWireDiode,_GaNJunctionKernel,T_FLOOR_PHYSICAL_K
from fsim_core.transport import Q_SI
from fsim_core.nitride_materials import KB_EV,bandgap,binary
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
   # F: the module's OWN accounting_residual_s field (not the independent
   # sum above) must itself close within the 1e-10 numerical budget -- a
   # mutant dropping "leak" (or any channel) from the module's internal
   # residual formula leaves accounting_residual_s pinned near that
   # channel's rate, which fails this check even though it never touches s.
   ok("residual field small",abs(r["accounting_residual_s"])<=max(1.,r["r_supply_s"])*1e-10)
r=evaluate_injection(d15,I_uA=1.,**a);ok("surface capture competition",math.isclose(r["f_capture"],100/(100+.2+.8+1),rel_tol=1e-12));ok("dot qfl capture",math.isclose(r["r_captured_s"],(r["r_supply_s"]-r["r_leakage_s"])*r["f_capture"]*.8*r["f_qfl_dot"],rel_tol=1e-12))

# B: the module computes r_e_leak_ratio/r_h_leak_ratio as
# exp(-barrier_{e,h}_eV / (KB_EV*T)) (nitride_nanowire_transport.py,
# evaluate_injection: re/rh lines).  Compare against that literal formula
# evaluated independently with math.exp from the SAME barrier_e_eV=0.21 /
# barrier_h_eV=0.08 literals passed in `a`, so a barrier swap inside the
# module changes which value math.exp(-0.21/kT) is compared against and
# the check fails (unlike comparing r_e to a swapped-run's r_h, which is
# invariant under a global e/h swap).
kT300=KB_EV*a["T_K"]
ok("r_e formula",math.isclose(r["r_e_leak_ratio"],math.exp(-.21/kT300),rel_tol=1e-12))
ok("r_h formula",math.isclose(r["r_h_leak_ratio"],math.exp(-.08/kT300),rel_tol=1e-12))
swap=evaluate_injection(d15,I_uA=1.,**dict(a,barrier_e_eV=.08,barrier_h_eV=.21));ok("barrier paths observable",r["r_e_leak_ratio"]!=r["r_h_leak_ratio"] and r["r_e_leak_ratio"]==swap["r_h_leak_ratio"])

# MEDIUM 5: propagate NitrideDiode.depletion's own punch-through/flat-band
# boolean. r is I_uA 1.0 / T_K 300 (the `a` default) on the 30 nm preset --
# the same point the Opus re-review computed independently (V_j 3.4435 V
# > V_bi 3.3657 V): above flat band the row must stay valid with a
# reason-free flag, depletion_field_kVcm pinned to 0, and the SAME words
# the planar Stark module (nitride_stark.py) uses for the identical
# condition.
ok("M5 flat_band True above V_bi at I_uA 1.0/300K",r["flat_band"] is True and r["depletion_regime"]=="flat_band" and r["V_j"]>r["V_bi"])
ok("M5 flat_band field F is 0",r["depletion_field_kVcm"]==0.)
below_flat=evaluate_injection(d15,I_uA=.001,**a)
ok("M5 depleted regime below flat band",below_flat["flat_band"] is False and below_flat["depletion_regime"]=="depleted" and below_flat["V_j"]<below_flat["V_bi"])

# A: f_qfl_dot uses transport.qfl_suppression's formula,
# min(1, exp(-(E_X_eV - V_j)/(KB_EV*T))), evaluated at the dot line E_X_eV
# against the diode's own V_j.  DIRECTIVE ROUND H2 rebuilt the kernel on
# the GaN reservoir (x_in=0.0), which raised V_j at 1 nA/300 K on the
# 15 nm preset from 2.06 V to ~3.09 V -- so the real single-dot line,
# 2.84 eV (Deshpande 2013 Fig.3c), is no longer sub-turn-on (that is
# H2's whole point: see the H2 block below).  Pick a synthetic E_X_eV
# ABOVE the new V_j instead (3.3 eV, a few kT over 3.09 V at 300 K) so
# f_qfl_dot still lands well below 0.5 here, independent of the H2 fix.
r_sub=evaluate_injection(d15,I_uA=.001,**dict(a,E_X_eV=3.3))
fd_expected=min(1.,math.exp(-(3.3-r_sub["V_j"])/kT300))
ok("sub-turn-on present",r_sub["f_qfl_dot"]<0.5)
ok("sub-turn-on qfl formula",math.isclose(r_sub["f_qfl_dot"],fd_expected,rel_tol=1e-9))
ok("r_captured scales with f_qfl_dot",math.isclose(r_sub["r_captured_s"],(r_sub["r_supply_s"]-r_sub["r_leakage_s"])*r_sub["f_capture"]*a["S_dot"]*r_sub["f_qfl_dot"],rel_tol=1e-9))

# H2 (root cause): the wire junction is GaN p-n with the InGaN disc as a
# capture target BELOW the GaN quasi-Fermi separation, not an InGaN-alloy
# junction -- the kernel must be built with x_in=0.0 for j0/V_bi/ideality.
# Consequence: at 1 nA/300 K on the 30 nm preset, V_j lands in the
# empirically inferred window (Deshpande 2013 terminal 2-20 V range minus
# the 2.38 Gohm*1nA=2.38 V IR drop [V]), and f_qfl_dot for the MEASURED
# 2.84 eV single-dot line saturates near 1 (both at 300 K and at the
# Fig.4 replay's 60 K), instead of the pre-fix kernel's 2.06 V / ~1e-13.
r284_300=evaluate_injection(d15,I_uA=.001,**dict(a,E_X_eV=2.84))
ok("H2 V_j at 1nA in Deshpande window",2.5<=r284_300["V_j"]<=3.2)
ok("H2 f_qfl_dot near 1 at 300K for measured 2.84eV dot",r284_300["f_qfl_dot"]>=0.5)
r284_60=evaluate_injection(d15,I_uA=.001,**dict(a,T_K=60.,E_X_eV=2.84))
ok("H2 f_qfl_dot near 1 at 60K for measured 2.84eV dot",r284_60["f_qfl_dot"]>=0.5)
# Mutant guard: a kernel rebuilt on the disc's own x_in (e.g. 0.25, or the
# 0.40 default) instead of forcing GaN (0.0) would make V_j depend on
# that composition; the fixed kernel does not, at the same I/T.
d_x25=wire_pin(preset="deshpande_2013_30nm",x_in=.25); d_x40=wire_pin(preset="deshpande_2013_30nm",x_in=.40)
ok("H2 mutant guard: V_j independent of disc x_in",evaluate_injection(d_x25,I_uA=.001,**a)["V_j"]==evaluate_injection(d_x40,I_uA=.001,**a)["V_j"])

# MEDIUM 2: make the x_in guard REAL. Before this round x_in=0.0 was the
# only guard, and it was inert -- _GaNJunctionKernel hardcoded
# binary("GaN") internally, so the "V_j independent of disc x_in" check
# above was true even for a kernel that never actually looked at any
# argument meant to force GaN. reservoir_material is now threaded through
# _GaNJunctionKernel (binary(reservoir_material) inside _ln_ni_gan): (a)
# the diode's OWN kernel-building method forces the literal "GaN"
# regardless of self.x_in (introspected directly, not inferred from V_j
# alone); (b) reservoir_material is demonstrably load-bearing -- swapping
# it for a narrow-gap binary standing in for what a mutant that let the
# disc composition leak into the reservoir kernel would do changes V_bi
# and V_j by volts, not millivolts, and pushes V_j for the 1 nA/300 K
# anchor completely outside the [2.5, 3.2] V H2 window (the same window
# the "H2 V_j at 1nA in Deshpande window" check above enforces on the
# real, correctly-forced kernel) -- so such a mutant would fail that
# window check.
ok("M2 production kernel forces reservoir_material='GaN' regardless of disc x_in",d_x25._kernel().reservoir_material=="GaN" and d_x40._kernel().reservoir_material=="GaN")
_kw_m2=dict(N_A=d15.N_A,N_D=d15.N_D,d_i_nm=d15.d_i_nm,d_active_nm=d15.d_active_nm,area_um2=d15.area_cm2*1e8,R_s_ohm=d15.R_s_ohm,n_ideality=d15.n_ideality,tau_SRH_ns=d15.tau_SRH_ns,eps_r=d15.eps_r,T=300.,x_in=0.,wl_thickness_nm=d15.d_active_nm,f_Rs_local=d15.f_Rs_local)
k_gan=_GaNJunctionKernel(reservoir_material="GaN",**_kw_m2); k_mutant=_GaNJunctionKernel(reservoir_material="InN",**_kw_m2)
J_1nA=1e-9/d15.area_cm2
ok("M2 reservoir_material is real (not dead like the old x_in argument)",not math.isclose(k_gan.vbi(300.),k_mutant.vbi(300.),rel_tol=1e-3))
ok("M2 mutant reservoir_material would fail the H2 window check",not (2.5<=k_mutant.vj_of_j(J_1nA,300.)<=3.2))
# f_qfl_dot_thermodynamic_limit is the SAME qfl_suppression formula
# evaluated at V_bi (the module's own reported built-in potential)
# instead of V_j -- an independent recomputation from the module's V_bi
# output, distinct from f_qfl_dot at a point (3.3 eV) chosen so V_j <
# E_X_eV < V_bi and the two columns visibly differ.
r_pair=evaluate_injection(d15,I_uA=.001,**dict(a,E_X_eV=3.3))
fd_lim_expected=min(1.,math.exp(-(3.3-r_pair["V_bi"])/kT300))
ok("f_qfl_dot_thermodynamic_limit formula",math.isclose(r_pair["f_qfl_dot_thermodynamic_limit"],fd_lim_expected,rel_tol=1e-9))
ok("f_qfl_dot_thermodynamic_limit is the ceiling (>= delivered f_qfl_dot)",r_pair["f_qfl_dot_thermodynamic_limit"]>=r_pair["f_qfl_dot"] and r_pair["f_qfl_dot_thermodynamic_limit"]>r_pair["f_qfl_dot"])

z=evaluate_injection(d15,I_uA=1.,**dict(a,eta_rad_matrix=0.));ok("surface loss no photon",z["r_matrix_radiative_s"]==z["accepted_background_s"]==0.)
below=evaluate_injection(d15,I_uA=1.,**dict(a,reservoir_energy_eV=1.9));ok("below exciton rejected",below["background_window_fraction"]==0.)

# C: tau_cap_ps<=0 is guarded to k_cap=0 (no capture channel) instead of
# the ZeroDivisionError a naive 1000/tau_cap_ps would raise; f_capture and
# r_captured_s must both read exactly 0, and the row stays valid (this is
# a degenerate but well-defined capture-competition limit, not a
# temperature-floor failure).
zero_cap=evaluate_injection(d15,I_uA=1.,**dict(a,tau_cap_ps=0.))
ok("tau_cap_ps 0 guarded",zero_cap["valid"] and zero_cap["f_capture"]==0. and zero_cap["r_captured_s"]==0.)

ok("minimum pair current",math.isclose(r["ideal_min_pair_current_A"],Q_SI/1e-10,rel_tol=1e-12))

# M9: V_bi is now a degenerate/frozen-carrier band-edge/doping formula
# (no ni**2), so the kernel is valid at every Deshpande anchor
# temperature from the as-measured 10 K bath up through 300 K --
# previously invalid below ~25-60 K (the pre-M9 ni/ni**2 underflow
# floors).  Check both evaluate_injection (V_j, mu) and
# wire_operating_point (T_j) directly through the public API.
for T9 in (10.,25.,60.,230.,300.):
 r9=evaluate_injection(d15,I_uA=.001,**dict(a,T_K=T9,E_X_eV=2.))
 ok("M9 evaluate_injection valid at %gK"%T9,r9["valid"] and math.isfinite(r9["V_j"]) and math.isfinite(r9["mu"]))
 o9=wire_operating_point(d15,I_uA=.001,T_hs_K=T9,duty=1.,Rth_K_W=WIRE_RTH_PRESETS_K_W["horizontal"],eta_total=.01,h_nu_eV=2.)
 ok("M9 wire_operating_point valid at %gK"%T9,o9["valid"] and math.isfinite(o9["T_j_K"]) and math.isfinite(o9["V_j"]))

# MEDIUM 3: T_FLOOR_PHYSICAL_K is a disclosed floor (default 10 K), not a
# silent underflow -- rows at T < 10 K must be invalid with the
# floor-specific reason, both through evaluate_injection and
# wire_operating_point, and the 10 K row itself (the boundary, and the
# Deshpande bath) must stay valid (M9, re-asserted above).
ok("M3 T_FLOOR_PHYSICAL_K is 10 K",T_FLOOR_PHYSICAL_K==10.)
r_below_floor=evaluate_injection(d15,I_uA=.001,**dict(a,T_K=5.))
ok("M3 evaluate_injection invalid below floor",not r_below_floor["valid"] and "below_physical_validity_floor" in r_below_floor["reasons"])
o_below_floor=wire_operating_point(d15,I_uA=.001,T_hs_K=5.,duty=1.,Rth_K_W=WIRE_RTH_PRESETS_K_W["horizontal"],eta_total=.01,h_nu_eV=2.)
ok("M3 wire_operating_point invalid below floor",not o_below_floor["valid"] and "below_physical_validity_floor" in o_below_floor["reasons"])
p_below_floor=pulse_delivery(d15,V_j=5.,T_K=5.,tau_pulse_ns=.1,rep_rate_hz=80e6)
ok("M3 pulse_delivery invalid below floor",not p_below_floor["valid"] and "below_physical_validity_floor" in p_below_floor["reasons"])
# LOW 7: pulse_delivery's invalid row NaN-fills every contract column the
# valid row carries, like evaluate_injection's _invalid() does, instead of
# returning a partial dict.
for k in ("C_dep_F","C_total_F","tau_RC_s","f_3dB_Hz","delivered_step_fraction","inter_pulse_residual_fraction"):
 ok("L7 invalid pulse row column "+k,math.isnan(p_below_floor[k]))
ok("L7 invalid pulse row feasible flag is False, not NaN",p_below_floor["pulse_delivery_feasible"] is False)

# MEDIUM 3: V_bi_minus_Eg_mV is (V_bi - E_g(reservoir_material, T)) * 1000,
# computed independently here from the module's own reported V_bi and
# nitride_materials.bandgap on the same GaN binary the kernel is forced to
# use -- not a re-derivation of V_bi itself. At 10 K it must be positive
# (V_bi exceeds E_g/q, the non-degenerate-formula bias the class docstring
# now discloses), in the 6-11 mV range the Opus re-review measured.
eg10=bandgap(binary("GaN"),10.)
r10_m3=evaluate_injection(d15,I_uA=.001,**dict(a,T_K=10.))
ok("M3 V_bi_minus_Eg_mV formula",math.isclose(r10_m3["V_bi_minus_Eg_mV"],(r10_m3["V_bi"]-eg10)*1e3,rel_tol=1e-9))
ok("M3 V_bi exceeds Eg/q by 6-11 mV at 10K",6.<=r10_m3["V_bi_minus_Eg_mV"]<=11.)

# MEDIUM 4: the H2 window check above already runs at 300 K (the `a`
# default). This is the second half: at 10 K (the as-measured Deshpande
# bath) V_j must sit within 20 mV of E_g(10K)/q REGARDLESS of tau_SRH_ns
# (the SRH-lifetime fit only sets the 300 K window position, not the
# low-T pinning) -- computed here independently from
# nitride_materials.bandgap, not from the module's own V_bi.
ok("M4 V_j at 10K within 20mV of Eg(10K)/q",abs(r10_m3["V_j"]-eg10)*1e3<=20.)

# The 10 K row is now genuinely valid (M9), so it can no longer exercise
# the invalid-row / NaN-column path; demonstrate that path still works
# on a genuinely pathological input instead -- an R_s_ohm so large that
# V_terminal itself overflows to a non-finite value, caught by
# _kernel()'s isfinite check ("kernel_error"), independent of any T floor.
low=evaluate_injection(wire_pin(preset="deshpande_2013_30nm",R_s_ohm=1e308),I_uA=1e7,**a)
ok("genuinely failing input (V_terminal overflow) still invalid",not low["valid"] and "kernel_error" in low["reasons"])
for k in ("eta_inj","f_capture","f_qfl_dot","f_qfl_dot_thermodynamic_limit","r_captured_s","r_matrix_radiative_s","r_matrix_nonradiative_s","r_surface_reservoir_s","r_leakage_s","mu","power_on_W","accounting_residual_s"):
 ok("invalid row column "+k,math.isnan(low[k]))

# M7: NitrideWireDiode carries no dead reservoir_surface_ns / tau_cap_ps
# dataclass fields -- both are keyword-only arguments to
# evaluate_injection instead, wired explicitly (surface_reservoir_ns,
# tau_cap_ps), never read off the diode object.
ok("M7 no dead reservoir_surface_ns/tau_cap_ps dataclass fields","reservoir_surface_ns" not in NitrideWireDiode.__dataclass_fields__ and "tau_cap_ps" not in NitrideWireDiode.__dataclass_fields__)

# M8: Rth provenance must state which heat path dominates and the SIGN
# of the 10K->300K transfer error (GaN falls, SiO2 rises with T).
thermal_doc=r9["provenance"]["thermal"]
ok("M8 Rth sign-of-transfer-error documented","GaN" in thermal_doc and "SiO2" in thermal_doc and "LOWER" in thermal_doc)

# L1: k_surface_reservoir_per_ns is an explicit, unambiguous alias for
# surface_reservoir_ns (a RATE despite the _ns suffix); either keyword
# alone suffices, both together must agree, and the hazard is documented.
r_direct=evaluate_injection(d15,I_uA=1.,**a)
r_alias=evaluate_injection(d15,I_uA=1.,**dict(a,surface_reservoir_ns=None,k_surface_reservoir_per_ns=a["surface_reservoir_ns"]))
ok("L1 alias matches direct keyword",r_alias["f_capture"]==r_direct["f_capture"] and r_alias["r_surface_reservoir_s"]==r_direct["r_surface_reservoir_s"])
try:
 evaluate_injection(d15,I_uA=1.,**dict(a,k_surface_reservoir_per_ns=a["surface_reservoir_ns"]*2.));ok("L1 contradiction raises",False)
except ValueError: ok("L1 contradiction raises",True)
try:
 evaluate_injection(d15,I_uA=1.,**dict(a,surface_reservoir_ns=None));ok("L1 missing both raises",False)
except TypeError: ok("L1 missing both raises",True)
ok("L1 naming hazard documented","RATE" in r_direct["provenance"]["naming"])

# HIGH 1: the Fig.4 anchor's bath is 10 K (Deshpande 2013 anchors.yaml:71),
# not 60 K -- the module's own _PROVENANCE["thermal"] already claimed the
# 10 K bath before this fix, but this replay ran at 60 K; both now agree.
# Print and check both the rise and the signed percent deviation from the
# paper's own +15 K / +49 K at that bath.
for u,target in ((.001,15.),(.002,49.)):
 o=wire_operating_point(d15,I_uA=u,T_hs_K=10.,duty=1.,Rth_K_W=WIRE_RTH_PRESETS_K_W["deshpande_fig4_replay"],eta_total=0.,h_nu_eV=2.);I=u*1e-6
 rise=o["T_j_K"]-10.;dev_pct=100.*(rise-target)/target
 print(f"Fig4 thermal anchor (10K bath): I={u} uA predicted +{rise:.2f} K vs Deshpande +{target:.0f} K ({dev_pct:+.1f} pct)")
 ok("independent heat",math.isclose(o["P_on_W"],I*o["V_j"]+I*I*d15.R_s_ohm,rel_tol=1e-10));ok("Fig4 thermal anchor",math.isclose(rise,target,rel_tol=.15))
# HIGH 1: the pre-H2 3.1e9 K/W value must FAIL the 15% band at 1 nA on the
# SAME 10 K bath -- it is a relic of the pre-H2 InGaN-alloy kernel, not an
# alternative anchor at this bath.
o_31=wire_operating_point(d15,I_uA=.001,T_hs_K=10.,duty=1.,Rth_K_W=3.1e9,eta_total=0.,h_nu_eV=2.)
rise_31=o_31["T_j_K"]-10.
print(f"Fig4 thermal anchor (10K bath, 3.1e9 relic): predicted +{rise_31:.2f} K vs Deshpande +15 K ({100.*(rise_31-15.)/15.:+.1f} pct)")
ok("HIGH1 3.1e9 K/W fails the 15pct band at 1nA/10K bath",not math.isclose(rise_31,15.,rel_tol=.15))
o=wire_operating_point(d15,I_uA=.001,T_hs_K=300.,duty=1.,Rth_K_W=0.,eta_total=.01,h_nu_eV=2.);ok("zero Rth electrical retained",o["T_j_K"]==300. and o["P_on_W"]>0 and o["V_j"]>0)
o=wire_operating_point(d15,I_uA=.001,T_hs_K=300.,duty=.1,Rth_K_W=1e7,eta_total=.01,h_nu_eV=2.)
ok("duty not squared average current",math.isclose(o["P_average_W"],.1*o["P_on_W"],rel_tol=1e-12))
# T1: independently pin electrical, optical, and local-heating terms.
_p = evaluate_injection(d15, I_uA=1., **dict(a, eta_total=.5, E_X_eV=2.84))
_ii = 1e-6
ok("T1 P_on includes optical subtraction", math.isclose(_p["power_on_W"], _ii*_p["V_j"] + _ii*_ii*d15.R_s_ohm - _ii*.5*2.84, rel_tol=1e-10))
_dhalf = wire_pin(preset="deshpande_2013_30nm", f_Rs_local=.5)
_ph = evaluate_injection(_dhalf, I_uA=1., **dict(a, eta_total=0.))
ok("T1 f_Rs_local halves series heat", math.isclose(_ph["power_on_W"], _ii*_ph["V_j"] + .5*_ii*_ii*d15.R_s_ohm, rel_tol=1e-10))
# E: T_j must equal T_hs + Rth*(duty*P_on), computed independently here
# with the LINEAR duty factor -- a mutant heating with duty**2 (a factor
# of 10 off at duty=0.1) fails this even though P_average alone (checked
# above from the module's own P_on) would not catch it.
ok("Tj matches independent duty-scaled rise",math.isclose(o["T_j_K"],300.+1e7*(.1*o["P_on_W"]),rel_tol=1e-6))

# LOW 8: wire_operating_point's internal out() helper no longer defaults
# reasons to a shared mutable list -- two independent calls must not
# return the same list object (which a `reasons=[]` default would, once
# something appended to it across calls).
o_l8a=wire_operating_point(d15,I_uA=.001,T_hs_K=300.,duty=1.,Rth_K_W=WIRE_RTH_PRESETS_K_W["horizontal"],eta_total=.01,h_nu_eV=2.)
o_l8b=wire_operating_point(d15,I_uA=.001,T_hs_K=300.,duty=1.,Rth_K_W=WIRE_RTH_PRESETS_K_W["horizontal"],eta_total=.01,h_nu_eV=2.)
ok("L8 no mutable-default reasons list shared across calls",o_l8a["reasons"] is not o_l8b["reasons"])

for cp,rate in ((0.,80e6),(1e-18,80e6),(1e-16,200e6)):
 p=pulse_delivery(d15,V_j=5.,T_K=300.,tau_pulse_ns=.1,rep_rate_hz=rate,C_parasitic_F=cp);tau=d15.R_s_ohm*(p["C_dep_F"]+cp)
 ok("RC independently summed",math.isclose(p["tau_RC_s"],tau,rel_tol=1e-12));ok("RC exponential independent",math.isclose(p["delivered_step_fraction"],1-math.exp(-1e-10/tau),rel_tol=1e-12))
ok("RC diagnostic",not pulse_delivery(d15,V_j=5.,T_K=300.,tau_pulse_ns=.1,rep_rate_hz=80e6)["pulse_delivery_feasible"])
print(f"{sum(C)}/{len(C)} nitride nanowire transport checks passed")
raise SystemExit(0 if all(C) else 1)
