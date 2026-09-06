import os,csv,numpy as np
import matplotlib.pyplot as plt
o=os.path.join(os.path.dirname(__file__),'out','02_22_landscape.png');os.makedirs(os.path.dirname(o),exist_ok=True);p=os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..','out','rt_edge','sweep.csv'))
with open(p,encoding='utf-8') as h:r=list(csv.DictReader(h))
g=sorted({float(x['gamma300_meV']) for x in r});d=sorted({float(x['delta_xx_meV']) for x in r});z=np.full((len(d),len(g)),np.nan)
for x in r:z[d.index(float(x['delta_xx_meV'])),g.index(float(x['gamma300_meV']))]=float(x['g2_pulsed'])
f,a=plt.subplots(figsize=(16,9),dpi=100);im=a.imshow(z,origin='lower',aspect='auto',extent=[min(g),max(g),min(d),max(d)],vmin=0,vmax=1,cmap='magma_r');f.colorbar(im,ax=a,label='pulsed g²(0)');a.set(xlabel='Γ (meV)',ylabel='ΔXX (meV)',title='Repository sweep: purity is jointly limited by linewidth and ladder separation');f.tight_layout();f.savefig(o,facecolor='white');plt.close(f)

