"""Source transcription, numerical and non-gating comparison checks."""
import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fsim_core.nitride_transport import planar_pin, evaluate_injection
from fsim_core.transport import Q_SI

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
r=evaluate_injection(d,1.,300.,1e10,.785,2.,1.,100.)
check("numerical mu units",abs(r.mu-r.loading.r_dot*2e-9)<1e-15)
check("numerical current conservation",r.loading.r_captured+r.loading.r_matrix<=1e-6/Q_SI*(1+1e-12))
small=d.dot_loading(1.,300.,1e10,.1,1.).r_dot; large=d.dot_loading(1.,300.,1e10,10.,1.).r_dot
check("numerical aperture scaling",large<small)
print("non-gating model comparison: computed 10 A/cm2 turn-on %.3f V; Zhang 2016 reports approximately 4 V" % d.vj_of_j(10.,300.))
print(f"{sum(c)}/{len(c)} nitride transport checks passed")
raise SystemExit(0 if all(c) else 1)
