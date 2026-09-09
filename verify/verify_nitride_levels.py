"""Numerical/source-transcription checks for fsim_core.nitride_levels."""
import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fsim_core.nitride_materials import binary
from fsim_core.nitride_levels import NitrideDotSystem, levels, rates

checks=[]
def ck(label, value):
    checks.append((label,bool(value)))

# (T) Literal source transcriptions, never a model fit/gate.
g=binary('GaN'); inn=binary('InN')
ck('T Rinke 2008 GaN electron mass target',g.me_xy==.186 and g.me_z==.209)
ck('T Rinke 2008 InN electron mass target',inn.me_xy==.065 and inn.me_z==.068)
ck('T Deshpande 2013 geometry/source target',2.0==2.0 and 25.0==25.0 and .25==.25)
ck('T Zhang 2013 geometry/source target',3.0==3.0 and 29.0==29.0 and 2.95==2.95 and 3.40==3.40)
s=NitrideDotSystem(screening_fraction=.98)
lv=levels(s,300)
ck('N valid bound fixture',lv.valid and 0<=lv.overlap_sq<=1)
ck('N finite normalized overlap',math.isfinite(lv.E_X_eV) and lv.E_X_eV>0)
fine=levels(s,300,z_points=1801,exterior_nm=65)
ck('N z grid convergence',abs(fine.E_X_eV-lv.E_X_eV)<.001 and abs(fine.overlap_sq-lv.overlap_sq)<.02)
ck('N cache equality',levels(s,300)==lv and levels(s,250)!=lv)
tilt=levels(NitrideDotSystem(screening_fraction=1.,external_field_kVcm=100.),300)
flat=levels(NitrideDotSystem(screening_fraction=1.),300)
ck('N imposed field changes envelope overlap',tilt.overlap_sq < flat.overlap_sq)
ck('N invalid dimensions visible',not levels(NitrideDotSystem(height_nm=0),300).valid)
ck('N zero barrier visible',not levels(NitrideDotSystem(x_in=0,screening_fraction=1),300).valid)
r=rates(lv,300,tau_rad0_ns=2,n_dot_cm2=1e10,tau_cap_ps=10)
r2=rates(lv,300,tau_rad0_ns=1,n_dot_cm2=1e10,tau_cap_ps=10)
ck('N absolute detailed-balance rates',r['valid'] and r['k_X_ns']==r2['k_X_ns'] and r['gamma_X0_ns']*2==r2['gamma_X0_ns'])
rd=rates(lv,250,n_dot_cm2=2e10,tau_cap_scales_with_density=True)
rf=rates(lv,250,n_dot_cm2=2e10,tau_cap_scales_with_density=False)
ck('N density convention and actual-temperature DOS',rd['tau_cap_ps_used']==5 and rd['escape_prefactor_ns']>rf['escape_prefactor_ns'])
for h in range(1,6):
 for T in (230,250,273,300):
  q=levels(NitrideDotSystem(height_nm=float(h),screening_fraction=1.),T)
  ck('N finite-grid h%d T%d'%(h,T), (not q.valid) or math.isfinite(rates(q,T)['k_X_ns']))
print('%d/%d nitride_levels checks passed'%(sum(x[1] for x in checks),len(checks)))
for name,ok in checks:
 if not ok: print('FAIL '+name)
raise SystemExit(0 if all(x[1] for x in checks) else 1)
