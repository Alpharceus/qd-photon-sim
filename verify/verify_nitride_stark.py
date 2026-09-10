"""Independent checks for fsim_core.nitride_stark."""
import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fsim_core.nitride_transport import planar_pin
from fsim_core.nitride_levels import NitrideDotSystem, levels
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

# AC3: negative junction-voltage rejection, transcribed explicitly (not just
# folded into the generic 'invalid bias' loop above).
neg_vj = resolve_bias(d, T_j_K=300., junction_voltage_V=-0.1)
ok('negative junction_voltage_V rejected', not neg_vj['bias_valid'] and
   any('nonnegative' in reason or 'reverse bias' in reason for reason in neg_vj['invalid_reasons']))

# AC3: R_s term reaches V_terminal in the junction-voltage-controlled branch
# specifically (mutation: deleting + current_A*R_s there survives if this
# is not checked separately from the current-controlled branch above).
# V_j=2.5 V is chosen (not 0.5 V) so the Shockley current, and hence the
# R_s*I term, is not so small it rounds away in float64.
b_jv = resolve_bias(d, T_j_K=300., junction_voltage_V=2.5)
current_A = b_jv['current_uA'] * 1e-6
ok('R_s term present in V_terminal (junction-voltage branch)',
   abs(b_jv['V_terminal'] - (b_jv['V_j'] + current_A*d.R_s_ohm)) < 1e-12)
d_hi_rs = planar_pin(R_s_ohm=d.R_s_ohm*5.)
b_hi_rs = resolve_bias(d_hi_rs, T_j_K=300., junction_voltage_V=2.5)
ok('R_s change alters V_terminal but not V_j (junction-voltage branch)',
   b_hi_rs['V_j'] == b_jv['V_j'] and b_hi_rs['V_terminal'] != b_jv['V_terminal'])

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

# AC4: 0/1/2-row traces must return derivative_valid=False rows without
# raising (regression guard for the len(rows)>=3 indexing bug).
ok('empty trace, no exception', stark_derivatives([]) == [])
one_row = stark_derivatives([{'row_id':0,'V_j':0.,'E_X_eV':1.,'spectroscopy_valid':True}])
ok('single-row trace is safe', len(one_row)==1 and not one_row[0]['derivative_valid'])
two_row = stark_derivatives([{'row_id':0,'V_j':0.,'E_X_eV':1.,'spectroscopy_valid':True},
                              {'row_id':1,'V_j':1.,'E_X_eV':.9,'spectroscopy_valid':True}])
ok('two-row trace is safe', len(two_row)==2 and not any(x['derivative_valid'] for x in two_row))

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
ok('nonpolar degeneracy', all(x['compatible'] is None and x['identification_status']=='screening_unidentifiable' for x in same))
life=screening_compatibility(fixture((-10.,)), slope_range_meV_per_V=(-10.,-10.), voltage_window_V=(0.,3.), lifetime_range_ns={'voltage_V':2.,'min_ns':.5,'max_ns':1.,'kind':'bare'})
ok('lifetime exact sample', life[0]['compatible'])
missing_life=screening_compatibility(fixture((-10.,)), slope_range_meV_per_V=(-10.,-10.), voltage_window_V=(0.,3.), lifetime_range_ns={'voltage_V':2.5,'min_ns':.5,'max_ns':1.,'kind':'bare'})
ok('missing lifetime is incomplete coverage', missing_life[0]['identification_status']=='incomplete_model_coverage')
gap=fixture((-10.,), False)
inc=screening_compatibility(gap, slope_range_meV_per_V=(-10.,-10.), voltage_window_V=(0.,3.))
ok('invalid coverage reported', inc[0]['identification_status']=='incomplete_model_coverage')

# AC5: narrow window (fewer than 3 samples land inside voltage_window_V) is
# reported incomplete_model_coverage, not fit from <3 points.
narrow=screening_compatibility(fixture((-10.,)), slope_range_meV_per_V=(-10.,-10.), voltage_window_V=(0.,1.5))
ok('narrow window is incomplete, not fit', narrow[0]['identification_status']=='incomplete_model_coverage' and narrow[0]['compatible'] is False)

short_trace=[{'row_id':i, 'V_j':v, 'E_X_eV':1.5-.010*v,
              'spectroscopy_valid':True} for i,v in enumerate((.45,.50,.55))]
short_result=screening_compatibility(short_trace, slope_range_meV_per_V=(-11.,-9.), voltage_window_V=(.25,.75))
ok('requested voltage window outside valid trace coverage is incomplete',
   short_result[0]['compatible'] is False and short_result[0]['identification_status']=='incomplete_model_coverage')

# AC5: bare-vs-cavity dispatch, window containing only ONE of the two
# lifetime values (mutation: swapping the bare/cavity key selection).
bc_rows=fixture((-10.,))
for row in bc_rows: row['tau_rad_bare_ns'], row['tau_rad_cavity_ns'] = .3, .8
bare_res=screening_compatibility(bc_rows, slope_range_meV_per_V=(-10.,-10.), voltage_window_V=(0.,3.),
    lifetime_range_ns={'voltage_V':2.,'min_ns':.5,'max_ns':1.,'kind':'bare'})
cavity_res=screening_compatibility(bc_rows, slope_range_meV_per_V=(-10.,-10.), voltage_window_V=(0.,3.),
    lifetime_range_ns={'voltage_V':2.,'min_ns':.5,'max_ns':1.,'kind':'cavity'})
ok('bare dispatch excluded by its own (out-of-window) lifetime', not bare_res[0]['compatible'])
ok('cavity dispatch matches its own (in-window) lifetime', cavity_res[0]['compatible'])

# AC2/finding fix: a tie the window definitively EXCLUDES must stay
# compatible=False/"incompatible", never screening_unidentifiable.
excluded_rows=fixture((-10.,-10.))
for number, row in enumerate(excluded_rows):
    row['screening_fraction'] = (0., 1.)[number // 4]; row['row_id'] = 'excl-%d' % number
excluded=screening_compatibility(excluded_rows, slope_range_meV_per_V=(4.,5.), voltage_window_V=(0.,3.))
ok('excluded tie stays compatible=False (not unidentifiable)',
   all(x['compatible'] is False and x['identification_status']=='incompatible' for x in excluded))

# AC2/finding fix: a tie where lifetime coverage is missing for BOTH members
# must stay incomplete_model_coverage, never get rewritten to unidentifiable.
incomplete_rows=fixture((-10.,-10.))
for number, row in enumerate(incomplete_rows):
    row['screening_fraction'] = (0., 1.)[number // 4]; row['row_id'] = 'inc-%d' % number
incomplete_rows[2]['tau_rad_bare_ns'] = float('nan')
incomplete_rows[6]['tau_rad_bare_ns'] = float('nan')
incomplete=screening_compatibility(incomplete_rows, slope_range_meV_per_V=(-10.,-10.), voltage_window_V=(0.,3.),
    lifetime_range_ns={'voltage_V':2.,'min_ns':.5,'max_ns':1.,'kind':'bare'})
ok('incomplete tie stays incomplete_model_coverage (not unidentifiable)',
   all(x['identification_status']=='incomplete_model_coverage' and x['compatible'] is False for x in incomplete))

# AC5: row-id reproducibility -- recompute each fitted slope from its own
# row_ids and voltage_window_V, independently of the module's internal
# ordinary-least-squares helper, and compare to 1e-12.
def _index_by_row_id(rows):
    return {row['row_id']: row for row in rows}
def _independent_slope(rows_subset):
    x=[float(row['V_j']) for row in rows_subset]; y=[1000.*float(row['E_X_eV']) for row in rows_subset]
    n=len(x); xm=sum(x)/n; ym=sum(y)/n
    num=sum((a-xm)*(b-ym) for a,b in zip(x,y)); den=sum((a-xm)**2 for a in x)
    return num/den
def check_reproducibility(label, records, index):
    for rec in records:
        if rec['branch_id'] is None:
            continue
        ok('%s row_ids has >=3 entries' % label, len(rec['row_ids']) >= 3)
        subset=[index[rid] for rid in rec['row_ids']]
        vlo,vhi=rec['voltage_window_V']
        ok('%s row_ids all inside stated window' % label, all(vlo<=float(row['V_j'])<=vhi for row in subset))
        indep=_independent_slope(subset)
        ok('%s slope reproducible from row_ids' % label, abs(indep - rec['fitted_slope_meV_per_V']) < 1e-12)
check_reproducibility('single-hypothesis', one, _index_by_row_id(fixture((-10.,-5.,2.))))
check_reproducibility('multi-hypothesis', many, _index_by_row_id(fixture((-10.,-9.,2.))))
check_reproducibility('bare-cavity', bare_res + cavity_res, _index_by_row_id(bc_rows))

# AC6: independently transcribed Zhang 2016 slope, compared to the module's
# own metadata; provenance/transcription gate only -- the production
# prediction is never required to reproduce it (see the pinned SLOPES
# checks below, which are all model predictions and do not target -10).
INDEPENDENT_ZHANG2016_SLOPE_MEV_PER_V = -10.0  # [V] Zhang et al., Appl. Phys. Lett. 108, 153102 (2016), Fig. 5
ok('Zhang transcription gate', abs(ZHANG2016_SLOPE_MEV_PER_V - INDEPENDENT_ZHANG2016_SLOPE_MEV_PER_V) < 1e-12
   and 'Zhang' in ZHANG2016_COMPARISON['citation'] and '2016' in ZHANG2016_COMPARISON['citation']
   and 'non-gating' in ZHANG2016_COMPARISON['use'])

flat=[{'row_id':i,'V_j':v,'E_X_eV':1.0-.01*v,'spectroscopy_valid':True,
       'flat_band': i >= 2} for i,v in enumerate((3.1,3.2,3.3,3.4))]
flat_out=stark_derivatives(flat)
ok('flat-band kink masks derivatives', not any(x['derivative_valid'] for x in flat_out))
invalid_row=stark_derivatives([{'V_j':0.,'E_X_eV':1.,'spectroscopy_valid':True},
    {'V_j':None,'E_X_eV':None,'spectroscopy_valid':False},
    {'V_j':2.,'E_X_eV':.98,'spectroscopy_valid':True}])
ok('invalid voltage row is skipped', not any(x['derivative_valid'] for x in invalid_row))

# Physics-coupled checks: 3 nm / 10 nm / x_in 0.25 / 300 K planar_pin dot,
# resolved through this module's OWN resolve_bias + stark_derivatives at
# V_j = 0.25/0.50/0.75 V, and levels() from fsim_core.nitride_levels.
# Reference slopes at V_j=0.50 V are the Opus re-review's independently
# computed SLOPES (nitride-stark-opus-rereview-findings.md), pinned to
# 1e-3 meV/V; recomputing them here from scratch (see the coder's own
# session notes) reproduced the same values to within about 2e-4 meV/V.
def stark_trace_slope(polarity, screening):
    diode = planar_pin()
    rows=[]
    for i, V in enumerate((0.25, 0.50, 0.75)):
        b = resolve_bias(diode, T_j_K=300.0, junction_voltage_V=V, field_polarity=polarity)
        sysd = NitrideDotSystem(height_nm=3.0, radius_nm=10.0, x_in=0.25,
                                 screening_fraction=screening, external_field_kVcm=b['applied_field_kVcm'])
        lv = levels(sysd, T_K=300.0)
        ok('physics-coupled trace point is valid', lv.valid and b['bias_valid'])
        rows.append({'row_id':i, 'V_j':V, 'E_X_eV':lv.E_X_eV, 'spectroscopy_valid':lv.valid})
    out=stark_derivatives(rows)
    ok('physics-coupled derivative at V_j=0.50 is valid', out[1]['derivative_valid'])
    return out[1]['dE_X_dV_meV_per_V']

PINNED_SLOPES_MEV_PER_V = {
    (1, 0.0): -12.9527, (1, 0.5): -9.5502, (1, 1.0): 2.8844,
    (-1, 0.0): 13.5229, (-1, 0.5): 10.6844, (-1, 1.0): 2.8844,
}
computed_slopes = {}
for (polarity, screening), expect in PINNED_SLOPES_MEV_PER_V.items():
    got = stark_trace_slope(polarity, screening)
    computed_slopes[(polarity, screening)] = got
    ok('pinned slope polarity=%+d screening=%.1f' % (polarity, screening), abs(got - expect) < 1e-3)

# Zero-polarization limit: at full screening the intrinsic polarization
# contribution vanishes, so the slope is set by the depletion field alone
# and must be the SAME regardless of field_polarity's sign convention.
ok('zero-polarization limit is polarity-independent',
   abs(computed_slopes[(1, 1.0)] - computed_slopes[(-1, 1.0)]) < 1e-6)
ok('zero-polarization limit matches pinned +2.8844 meV/V',
   abs(computed_slopes[(1, 1.0)] - 2.8844) < 1e-3)
# Polarity sign rule (screening=0, unscreened intrinsic field dominates):
# field_polarity=+1 gives a negative slope, -1 gives a positive slope.
ok('polarity sign rule at screening=0', computed_slopes[(1, 0.0)] < 0 < computed_slopes[(-1, 0.0)])
# Screening ordering at fixed polarity +1: increasing screening moves the
# slope from strongly negative toward positive, monotonically.
ok('screening ordering at polarity +1',
   computed_slopes[(1, 0.0)] < computed_slopes[(1, 0.5)] < computed_slopes[(1, 1.0)])

print('Zhang %.1f meV/V (%s) is a non-gating PL comparison, not a fit target; '
      'production Stark slopes above are model predictions and screening '
      'identification remains conditional on the supplied hypotheses.'
      % (ZHANG2016_SLOPE_MEV_PER_V, ZHANG2016_COMPARISON['citation']))
failed=[n for n,p in checks if not p]
for name in failed: print('FAILED:',name)
print('%d/%d nitride Stark diagnostics checks passed' % (len(checks)-len(failed),len(checks)))
raise SystemExit(1 if failed else 0)
