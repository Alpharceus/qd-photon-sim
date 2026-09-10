"""Independent checks for fsim_core.nitride_stark."""
import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fsim_core.nitride_transport import planar_pin
from fsim_core.nitride_stark import (resolve_bias, stark_derivatives,
    screening_compatibility, ZHANG2016_SLOPE_MEV_PER_V, ZHANG2016_COMPARISON)

checks = []
def ok(name, condition):
    checks.append((name, bool(condition)))

d = planar_pin()
# Independent SI depletion transcription [DR Sze & Ng, 2007].
T = 300.0; I = 0.02
r = resolve_bias(d, T_j_K=T, current_uA=I, external_field_kVcm=3.0)
eps = d.eps_r * 8.8541878128e-12; q = 1.602176634e-19
NA, ND, di = d.N_A*1e6, d.N_D*1e6, d.d_i_nm*1e-9
F = r['diode_field_kVcm'] * 1e5
drop = F*di + eps*F*F*(1/NA + 1/ND)/(2*q)
ok('abrupt depletion equation', abs((d.vbi(T)-r['V_j'])-drop) < 1e-10)
ok('flat band field', resolve_bias(d, T_j_K=T, junction_voltage_V=d.vbi(T))['diode_field_kVcm'] == 0.0)
minus = resolve_bias(d, T_j_K=T, current_uA=I, field_polarity=-1, external_field_kVcm=3.0)
ok('polarity only changes applied field', r['diode_field_kVcm'] == minus['diode_field_kVcm'] and r['applied_field_kVcm'] == -minus['applied_field_kVcm'] + 6.0)
for temp in (230., 300.):
    for current in (0., .002, .02, .2):
        a = resolve_bias(d, T_j_K=temp, current_uA=current)
        b = resolve_bias(d, T_j_K=temp, junction_voltage_V=a['V_j'])
        ok('roundtrip %.0f %.3g' % (temp,current), a['bias_valid'] and abs(b['current_uA']-current) < 1e-11 and abs(a['V_terminal']-(a['V_j']+current*1e-6*d.R_s_ohm)) < 1e-14)
for kwargs in ({}, {'current_uA':1,'junction_voltage_V':1}, {'current_uA':-1}, {'T_j_K':float('nan'),'current_uA':1}):
    args = dict(T_j_K=300, current_uA=None, junction_voltage_V=None); args.update(kwargs)
    ok('invalid bias', not resolve_bias(d, **args)['bias_valid'])

def poly_rows(grid, bad=None):
    return [{'row_id':i, 'V_j':v, 'E_X_eV':1.2-.003*v+.0004*v*v,
             'spectroscopy_valid': i != bad} for i,v in enumerate(grid)]
for grid in ([0,1,2,3,4], [0,.2,1.3,2.1,4.7]):
    out=stark_derivatives(poly_rows(grid))
    ok('quadratic derivative', all(abs(x['dE_X_dV_meV_per_V']-1000*(-.003+.0008*x['V_j'])) < 1e-9 for x in out))
bad=stark_derivatives(poly_rows([0,1,2,3,4], 2))
ok('invalid gaps not bridged', not any(x['derivative_valid'] for x in bad))
for grid in ([0,1,1], [0,2,1]):
    try: stark_derivatives(poly_rows(grid)); passed=False
    except ValueError: passed=True
    ok('strict trace ordering', passed)

def fixture(slopes, valid=True):
    rows=[]
    for s in slopes:
        for i,v in enumerate((0.,1.,2.,3.)):
            rows.append({'row_id':'%s-%s'%(s,i),'geometry':'g','temperature_mode':'fixed_junction',
                'screening_fraction':s,'V_j':v,'E_X_eV':1.5+(s*.001)*v,
                'spectroscopy_valid':valid, 'tau_rad_bare_ns':1.,'tau_rad_cavity_ns':.5})
    return rows
one=screening_compatibility(fixture((-10.,-5.,2.)), slope_range_meV_per_V=(-10.,-10.), voltage_window_V=(0.,3.))
ok('single compatible hypothesis', sum(x['compatible'] for x in one)==1)
many=screening_compatibility(fixture((-10.,-9.,2.)), slope_range_meV_per_V=(-10.,-9.), voltage_window_V=(0.,3.))
ok('inclusive slope bounds', sum(x['compatible'] for x in many)==2)
none=screening_compatibility(fixture((-10.,-5.,2.)), slope_range_meV_per_V=(4.,5.), voltage_window_V=(0.,3.))
ok('no match', not any(x['compatible'] for x in none))
same_rows=fixture((-10.,-10.,-10.))
for number, row in enumerate(same_rows): row['screening_fraction'] = (0., .5, 1.)[number // 4]
same=screening_compatibility(same_rows, slope_range_meV_per_V=(-10.,-10.), voltage_window_V=(0.,3.))
ok('nonpolar degeneracy', all(x['compatible'] and x['identification_status']=='screening_unidentifiable' for x in same))
life=screening_compatibility(fixture((-10.,)), slope_range_meV_per_V=(-10.,-10.), voltage_window_V=(0.,3.), lifetime_range_ns={'voltage_V':2.,'min_ns':.5,'max_ns':1.,'kind':'bare'})
ok('lifetime exact sample', life[0]['compatible'])
gap=fixture((-10.,), False)
inc=screening_compatibility(gap, slope_range_meV_per_V=(-10.,-10.), voltage_window_V=(0.,3.))
ok('invalid coverage reported', inc[0]['identification_status']=='incomplete_model_coverage')
ok('Zhang transcription', ZHANG2016_SLOPE_MEV_PER_V == -10.0 and 'non-gating' in ZHANG2016_COMPARISON['use'])
print('Zhang -10 meV/V is a non-gating PL comparison; screening remains conditional.')
failed=[n for n,p in checks if not p]
for name in failed: print('FAILED:',name)
print('%d/%d nitride Stark diagnostics checks passed' % (len(checks)-len(failed),len(checks)))
raise SystemExit(1 if failed else 0)
