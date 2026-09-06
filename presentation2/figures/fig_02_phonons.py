import os,sys,numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0,os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..')))
from fsim_core.qd_gf import PhononParams,ibm_spectrum,zpl_weight
o=os.path.join(os.path.dirname(__file__),'out','02_16_phonons.png');os.makedirs(os.path.dirname(o),exist_ok=True)
x=np.linspace(-8,8,1201);p=PhononParams();f,a=plt.subplots(figsize=(16,9),dpi=100)
for T,c in [(4,'#1769aa'),(80,'#e08e0b'),(300,'#c43d3d')]:
    y=ibm_spectrum(x,p,T,.15);a.plot(x,y/np.trapezoid(y,x),lw=3,c=c,label=f'{T} K, ZPL={zpl_weight(p,T):.2f}')
a.set(xlabel='detuning (meV)',ylabel='normalized spectrum',title='Independent-boson lineshape: ZPL plus acoustic-phonon sideband');a.legend(frameon=False,fontsize=16);a.grid(alpha=.2);f.tight_layout();f.savefig(o,facecolor='white');plt.close(f)
