"""Source transcription, numerical and non-gating comparison checks."""
import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fsim_core.nitride_transport import planar_pin, evaluate_injection
from fsim_core.nitride_materials import KB_EV, EPS0_SI
from fsim_core.transport import Q_SI, qfl_suppression, xi_window

c=[]
def check(label, value):
    c.append(bool(value))
    if not value: print("FAIL",label)
d=planar_pin()
# Source transcription: Zhang et al., Appl. Phys. Lett. 108, 153102 (2016),
# p.2 and Fig.2 [V].  Planar transfer is [E], so this is not device validation.
check("source transcription Zhang doping",d.N_A==1e17 and d.N_D==1e18)
check("source transcription Zhang active and spacers",d.d_active_nm==3 and d.d_i_nm==27)
# Numerical verification: analytic log-domain Shockley inverse [E/A].
for T in (230.,250.,273.,300.):
    I=1e-6; va,vj=d.v_of_i(I,T); back=d.j_of_vj(vj,T)*d.area_cm2
    check("numerical IV inverse",abs(back/I-1)<1e-7 and math.isfinite(va))
    r=evaluate_injection(d,1.,T,1e10,.785,1.,1.,100.)
    check("numerical finite operating point",math.isfinite(r.mu) and r.P_junction_W>=0 and r.b_e>=0)
    r_e=evaluate_injection(d,1.,T,1e10,.785,1.,1.,100.,E_X_eV=3.06)
    check("numerical finite operating point with E_X_eV",math.isfinite(r_e.mu) and math.isfinite(r_e.b_e) and r_e.P_junction_W>=0 and r_e.b_e>=0)
r=evaluate_injection(d,1.,300.,1e10,.785,2.,1.,100.)
check("numerical mu units",abs(r.mu-r.loading.r_dot*2e-9)<1e-15)
check("numerical current conservation",r.loading.r_captured+r.loading.r_matrix<=1e-6/Q_SI*(1+1e-12))
small=d.dot_loading(1.,300.,1e10,.1,1.).r_dot; large=d.dot_loading(1.,300.,1e10,10.,1.).r_dot
check("numerical aperture scaling",large<small)

# --- Regression checks for the Opus review findings on commit b91a3a5 ---

# High finding 1 (nitride_transport.py:63): b_e hardcoded the background's
# OWN sub-turn-on suppression to 1.0 while its denominator (r_dot, inside
# r_matrix's ancestry) already carried f_qfl -- a one-sided asymmetry that
# blew b_e up to O(1e3) instead of the correctly-suppressed O(1e-4).  This
# check passes a non-trivial (neither 0 nor 1) E_X_eV, exercises
# dot_loading's f_qfl path, and independently recomputes b_e using
# transport.py's OWN qfl_suppression / xi_window (not a private
# reimplementation) so the SAME convention is being compared, not a
# tautology of the module's internals.
E_X_eV=3.06; T_be=300.; I_be=1.; w_meV_be=1.; dE_WL_meV_be=100.; eta_rad_matrix_be=.1
ld_be=d.dot_loading(I_be,T_be,1e10,.785,None,None,1.,None,E_X_eV)
check("numerical b_e sub-turn-on point has non-trivial f_qfl",0.<ld_be.f_qfl<1.)
kT_meV_be=1e3*KB_EV*T_be
eu_be=max(kT_meV_be,0.)
xi_be=float(xi_window(w_meV_be,dE_WL_meV_be,eu_be))
f_qfl_bg_expected=qfl_suppression(E_X_eV+dE_WL_meV_be*1e-3,ld_be.V_j,KB_EV*T_be)
rate_bg_expected=ld_be.r_matrix*eta_rad_matrix_be*xi_be*f_qfl_bg_expected
rate_x_expected=min(ld_be.r_dot,1./(1.*1e-9))*1.
b_e_expected=rate_bg_expected/rate_x_expected if rate_x_expected>0 else float('inf')
bg_actual=d.b_e(I_be,T_be,w_meV_be,dE_WL_meV_be,1e10,.785,1.,eta_rad_matrix_be,None,1.,None,None,E_X_eV)
check("numerical b_e applies transport.py's own f_qfl_bg convention",math.isclose(bg_actual.b_e,b_e_expected,rel_tol=1e-9) and bg_actual.f_qfl_bg<1.)
# The specific regression the finding reported: b_e must stay O(1), not blow
# up to O(1e3), once the background's own suppression is applied.
check("numerical b_e sub-turn-on magnitude is suppressed, not amplified",bg_actual.b_e<1.)

# High finding 2 (nitride_transport.py:45): dV_active_mV was passed as
# F * d_active_nm with NO unit conversion (V/m times nm, missing 1e-9 m/nm
# and 1e3 mV/V), reporting ~1e6x too large.  Hand-compute the abrupt p-i-n
# depletion field from base SI constants and this diode's own N_A/N_D/eps_r/
# d_i_nm (the same physics as transport.Diode.depletion), independently of
# the depletion() method's own arithmetic, and compare dV_active_mV to it.
T_dep=300.; V_j_dep=1.0
vbi_dep=d.vbi(T_dep)
eps_dep=d.eps_r*EPS0_SI; d_i_dep=d.d_i_nm*1e-9
NA_dep,ND_dep=d.N_A*1e6,d.N_D*1e6
c_dep=vbi_dep-V_j_dep
a_dep=(Q_SI*NA_dep/(2.*eps_dep))*(1.+NA_dep/ND_dep)
b_dep=Q_SI*NA_dep*d_i_dep/eps_dep
x_p_dep=(-b_dep+math.sqrt(b_dep*b_dep+4.*a_dep*c_dep))/(2.*a_dep)
F_dep=Q_SI*NA_dep*x_p_dep/eps_dep  # V/m
dV_active_mV_hand=F_dep*d.d_active_nm*1e-9*1e3  # V/m -> m -> mV
dep_actual=d.depletion(V_j_dep,T_dep)
check("numerical dV_active_mV unit check vs hand-computed field",math.isclose(dep_actual.dV_active_mV,dV_active_mV_hand,rel_tol=1e-9))
# Sanity bound on the unit bug itself: an mV drop across a 3 nm layer at a
# kV/cm-scale field cannot be anywhere near 1e6 mV (~1e3 V).
check("numerical dV_active_mV is a physically sane millivolt-scale number",0.<=dep_actual.dV_active_mV<1e4)

print("non-gating model comparison: computed 10 A/cm2 turn-on %.3f V; Zhang 2016 reports approximately 4 V" % d.vj_of_j(10.,300.))
print(f"{sum(c)}/{len(c)} nitride transport checks passed")
raise SystemExit(0 if all(c) else 1)
