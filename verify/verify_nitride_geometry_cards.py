import os, sys, yaml
ROOT=os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0,ROOT)
from fsim_core.nitride_levels import NitrideDotSystem, levels

CARDS=['nitride-nonpolar-pulse-design.yaml','nitride-nonpolar-set-design.yaml','nitride-qw-fluctuation-pulse-design.yaml','nitride-qw-fluctuation-set-design.yaml']
def main():
    checks=0; failures=[]
    def ok(cond,msg):
        nonlocal checks
        checks+=1
        if not cond: failures.append(msg)
    for name in CARDS:
        with open(os.path.join(ROOT,'cards',name),encoding='utf-8') as f: d=yaml.safe_load(f)['design']
        n=d['nitride']['dot']; ok(n['orientation'] in ('a_plane','c_plane'),name+' orientation'); ok(n['geometry_type'] in ('isolated_dot','qw_fluctuation'),name+' geometry')
        ok(float(d['drive']['rep_rate_hz'])==8.0e7 and float(d['drive']['diode']['tau_pulse_ns'])==0.1,name+' timing')
        ok(all(k in n for k in ('height_nm','radius_nm','x_in','screening_fraction')),name+' geometry fields')
        s=NitrideDotSystem(**n); lv=levels(s,300.0); ok(lv.reservoir_kind in ('gan_barrier','ingan_qw'),name+' reservoir')
    anchors=yaml.safe_load(open(os.path.join(ROOT,'verify','data','nitride_geometry_stark_anchors.yaml'),encoding='utf-8'))['anchors']
    ok(anchors['wang2017_linewidth_g2']['value'].startswith('19.0'), 'Wang transcription')
    ok(anchors['zhang2016_stark_slope']['value']==-10.0, 'Zhang slope')
    ok(anchors['schade2011_orientation']['tag']=='A', 'Schade caveat')
    if failures:
        for f in failures: print('FAIL '+f)
        print(f'{checks-len(failures)}/{checks} nitride geometry card checks passed'); return 1
    print(f'{checks}/{checks} nitride geometry card checks passed'); return 0
if __name__=='__main__': sys.exit(main())
