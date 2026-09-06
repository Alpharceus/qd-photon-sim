import os,sys,numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0,os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..')))
from fsim_core import cw_g2
o=os.path.join(os.path.dirname(__file__),'out','02_15_irf.png');os.makedirs(os.path.dirname(o),exist_ok=True)
t=np.linspace(-.08,.08,1601);f,a=plt.subplots(figsize=(16,9),dpi=100)
for w,c in [(0,'#1769aa'),(20,'#e08e0b'),(100,'#c43d3d')]:
 g=1-.95*np.exp(-np.abs(t)/.005);g=cw_g2.convolve_irf(t,g,w) if w else g;a.plot(t*1000,g,lw=3,c=c,label=f'{w} ps IRF')
a.set(xlabel='delay (ps)',ylabel='g²',title='A 5 ps antibunching dip is erased by a 100 ps detector response');a.legend(frameon=False,fontsize=16);a.grid(alpha=.2);f.tight_layout();f.savefig(o,facecolor='white');plt.close(f)

