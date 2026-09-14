import math,sys
from pathlib import Path
sys.path.insert(0,str(Path(__file__).resolve().parents[1]))
from fsim_core.nitride_nanowire_levels import NitrideNanowireSystem,levels,rates,_ep,J01,J11
from fsim_core.nitride_materials import EPS0_SI
c=[]
def C(n,x):
 c.append(bool(x))
 if not x:print('FAIL',n)
s=NitrideNanowireSystem(x_in=.4,strain_bound='relaxed');l=levels(s)
psp=.6*(-.029)+.4*(-.032);f=(-.029-psp)/(EPS0_SI*(.6*10.28+.4*14.61))*1e-5
C('Bernardini spontaneous relaxed field',math.isfinite(l.F_sp_kVcm) and abs(l.F_pz_kVcm)<1e-12)
C('screening retains external',levels(NitrideNanowireSystem(screening_fraction=1,external_field_kVcm=7)).field_kVcm==7)
C('Bessel tabulation',math.isclose(_ep(.2,10),10.999,rel_tol=.02) and math.isclose(_ep(.2,10,J11)/_ep(.2,10),(J11/J01)**2))
C('R inverse square',math.isclose(_ep(.2,20),_ep(.2,10)/4))
C('invalid geometry',not levels(NitrideNanowireSystem(core_radius_nm=-1)).valid)
if l.valid:
 r=rates(l,300,tau_rad0_ns=1,tau_cap_ps=10,reservoir_length_nm=30);C('reservoir and XX rate',r['valid'] and r['reservoir_state_count_e']>0 and r['gamma_XX0_ns']==2*r['gamma_X0_ns'])
print(f'{sum(c)}/{len(c)} nitride nanowire levels checks passed')
raise SystemExit(0 if all(c) else 1)
