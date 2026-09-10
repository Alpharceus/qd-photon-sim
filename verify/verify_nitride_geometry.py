"""Independent effective-shape and QW-fluctuation checks.

Volumes use elementary profile integration [DR], not production helpers.
BenDaniel-Duke roots are independently tested in verify_nitride_levels.py.
"""
import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fsim_core.nitride_levels import NitrideDotSystem, levels, rates

checks=[]
def ck(name, value): checks.append((name,bool(value)))

H,R=6.,12.
base=dict(height_nm=H,radius_nm=R,screening_fraction=1.)
disc=levels(NitrideDotSystem(**base),300.)
lens=levels(NitrideDotSystem(**{**base,'shape':'lens'}),300.)
cone=levels(NitrideDotSystem(**{**base,'shape':'truncated_cone','top_radius_fraction':0.}),300.)
frustum=levels(NitrideDotSystem(**{**base,'shape':'truncated_cone','top_radius_fraction':.5}),300.)
unit=levels(NitrideDotSystem(**{**base,'shape':'truncated_cone','top_radius_fraction':1.}),300.)
ck('G independently integrated disc volume',abs(disc.shape_volume_nm3-math.pi*R*R*H)<1e-10)
ck('G lens volume and effective height H/2',abs(lens.shape_volume_nm3-math.pi*R*R*H/2)<1e-10 and lens.effective_height_nm==H/2)
ck('G cone volume and effective height H/3',abs(cone.shape_volume_nm3-math.pi*R*R*H/3)<1e-10 and cone.effective_height_nm==H/3)
ck('G frustum integral',abs(frustum.shape_volume_nm3-math.pi*R*R*H*(1+.5+.25)/3)<1e-10)
ck('G t=1 frustum exactly reduces to disc',disc.E_X_eV==unit.E_X_eV and disc.overlap_sq==unit.overlap_sq and disc.dE_e_meV==unit.dE_e_meV)
ck('G shape changes confinement',lens.valid and disc.valid and lens.E_X_eV!=disc.E_X_eV)
mapped=levels(NitrideDotSystem(height_nm=H/2,radius_nm=R,screening_fraction=1.),300.)
ck('G equal-volume mapped disc agrees exactly',lens.E_X_eV==mapped.E_X_eV and lens.overlap_sq==mapped.overlap_sq)

old=levels(NitrideDotSystem(height_nm=3.,radius_nm=10.),300.)
explicit=levels(NitrideDotSystem(height_nm=3.,radius_nm=10.,orientation='c_plane',shape='disc',geometry_type='isolated_dot',top_radius_fraction=1.),300.)
ck('G omitted defaults preserve legacy scalars',old.E_X_eV==explicit.E_X_eV and old.overlap_sq==explicit.overlap_sq and old.dE_e_meV==explicit.dE_e_meV and old.field_kVcm==explicit.field_kVcm)

qargs=dict(geometry_type='qw_fluctuation',height_nm=3.5,radius_nm=20.,wl_thickness_nm=3.,screening_fraction=1.)
q=levels(NitrideDotSystem(**qargs),300.)
qwide=levels(NitrideDotSystem(**{**qargs,'radius_nm':30.}),300.)
qsmall=levels(NitrideDotSystem(**{**qargs,'radius_nm':5.}),300.)
ck('Q QW reservoir and positive escape',q.valid and q.reservoir_kind=='ingan_qw' and q.dE_e_meV>0 and q.dE_h_meV>0)
ck('Q reservoir is radius independent',q.reservoir_energy_eV==qwide.reservoir_energy_eV==qsmall.reservoir_energy_eV)
ck('Q local energy varies with radius',q.E_X_eV!=qsmall.E_X_eV)
ck('Q electron hole and min rates valid',all(rates(q,300.,channel=c)['valid'] for c in ('electron','hole','min')))
ck('Q pair rate rejected',not rates(q,300.,channel='pair')['valid'])
equal=levels(NitrideDotSystem(geometry_type='qw_fluctuation',height_nm=3.,radius_nm=5.,wl_thickness_nm=3.,screening_fraction=1.),300.)
ck('Q H equals well is invalid',not equal.valid)

invalid=[]
for kw in (dict(shape='bad'),dict(shape='lens',top_radius_fraction=.5),dict(shape='truncated_cone',top_radius_fraction=1.1),dict(shape='lens',shape_height_fraction=.4),dict(geometry_type='qw_fluctuation',shape='lens',wl_thickness_nm=2.)):
    invalid.append(not levels(NitrideDotSystem(**kw),300.).valid)
ck('G malformed geometry rejected',all(invalid))

coarse=levels(NitrideDotSystem(height_nm=3.,radius_nm=10.,screening_fraction=1.),300.,z_points=1201,exterior_nm=45.)
fine=levels(NitrideDotSystem(height_nm=3.,radius_nm=10.,screening_fraction=1.),300.,z_points=2401,exterior_nm=90.)
ck('N refined bound fixture within 0.5 meV',coarse.valid and fine.valid and abs(coarse.E_X_eV-fine.E_X_eV)*1000<=.5 and abs(coarse.dE_e_meV-fine.dE_e_meV)<=.5 and abs(coarse.dE_h_meV-fine.dE_h_meV)<=.5)

print('%d/%d nitride_geometry checks passed'%(sum(ok for _,ok in checks),len(checks)))
for name,ok in checks:
    if not ok: print('FAIL '+name)
raise SystemExit(0 if all(ok for _,ok in checks) else 1)
