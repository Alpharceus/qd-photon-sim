import os,sys,numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0,os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..')))
from fsim_core import cw_g2
o=os.path.join(os.path.dirname(__file__),'out','02_14_cw.png');os.makedirs(os.path.dirname(o),exist_ok=True)
t=np.linspace(-8,8,1001);f,a=plt.subplots(figsize=(16,9),dpi=100)
for T,c,r in [(80,'#1769aa',.995),(300,'#c43d3d',.9)]:
 q=cw_g2.cw_report(.05,1,2,0,0,.5,.2,r,irf_fwhm_ps=0,tau_max_ns=8,n_tau=1001);a.plot(t,q['curves']['g2_meas'],lw=3,c=c,label=f'{T} K')
a.axhline(1,c='.5',ls='--');a.set(xlabel='delay τ (ns)',ylabel='g²(τ)',title='CW HBT correlation from the cap-2 rate-equation propagator');a.legend(frameon=False,fontsize=16);a.grid(alpha=.2);f.tight_layout();f.savefig(o,facecolor='white');plt.close(f)

