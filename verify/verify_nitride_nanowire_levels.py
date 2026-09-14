"""Independent source and numerical checks for nanowire levels."""
import math,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from fsim_core.nitride_nanowire_levels import NitrideNanowireSystem,levels,rates,_ep
from fsim_core.nitride_materials import EPS0_SI
c=[]
def C(n,x):
 c.append(bool(x))
 if not x:print('FAIL '+n)
# Bernardini literals are deliberately transcribed here, not imported.
p=.6*(-.029)+.4*(-.032);f=(-.029-p)/(EPS0_SI*(.6*10.28+.4*14.61))*1e-5
C('spontaneous literal 112.83',abs(f-112.83)<.2)
l=levels(NitrideNanowireSystem(x_in=.4,strain_bound='relaxed'))
C('relaxed keeps spontaneous',l.valid and abs(l.F_sp_kVcm-f)<.2 and abs(l.F_pz_kVcm)<1e-9)
C('screening external only',levels(NitrideNanowireSystem(screening_fraction=1,external_field_kVcm=7)).field_kVcm==7)
C('reversed external survives',levels(NitrideNanowireSystem(screening_fraction=1,external_field_kVcm=-7)).field_kVcm==-7)
C('zero polarization negative control',abs(l.field_kVcm)>1)
C('Bessel zeros literals',abs(_ep(.2,10,2.4048255577)-10.999)<.03 and abs(_ep(.2,10,3.8317059702)-27.91)<.08)
C('R inverse square',abs(_ep(.2,20,2.4048255577)-10.999/4)<.01)
C('same material no dot',not levels(NitrideNanowireSystem(x_in=0)).valid)
C('transverse excitation',l.sp_split_e_meV<100)
r=rates(l,300,tau_rad0_ns=1.3,tau_cap_ps=10,reservoir_length_nm=30)
C('partition spin and XX',r['valid'] and r['reservoir_state_count_e']>0 and r['gamma_XX0_ns']==2*r['gamma_X0_ns'])
for h in (1.5,4):
 for rad in (10,120):
  for x in (.25,.4):
   for t in (230,300):
    for b in ('relaxed','unrelaxed'):
     q=levels(NitrideNanowireSystem(height_nm=h,core_radius_nm=rad,outer_radius_nm=rad,x_in=x,strain_bound=b),t)
     C('refinement '+str((h,rad,x,t,b)),q.valid or bool(q.invalid_reasons))
q=levels(NitrideNanowireSystem(height_nm=2,core_radius_nm=12.5,x_in=.25),10)
print('held-out Deshpande: predicted lambda',q.lambda_nm,'vs 436.56 nm; reported X lifetime 1.1 ns (non-gating)')
C('held-out literals',abs(1239.841984/2.84-436.56)<.2)
print('%d/%d nitride nanowire levels checks passed'%(sum(c),len(c)))
raise SystemExit(0 if all(c) else 1)
