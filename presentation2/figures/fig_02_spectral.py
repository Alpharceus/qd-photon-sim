import os,sys,numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0,os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..')))
from fsim_core.spectral import epsilon
o=os.path.join(os.path.dirname(__file__),'out','02_19_spectral.png');os.makedirs(os.path.dirname(o),exist_ok=True)
G=np.linspace(.2,20,250);f,a=plt.subplots(figsize=(16,9),dpi=100)
for d,c in zip([4,6,8],['#1769aa','#e08e0b','#c43d3d']):a.plot(G,[epsilon(d,g,g,w=g).eps for g in G],lw=3,c=c,label=f'ΔXX={d} meV')
a.set(xlabel='FWHM Γ (meV)',ylabel='ε=tXX/tX',title='Spectral overlap increases as thermal linewidth grows');a.legend(frameon=False,fontsize=16);a.grid(alpha=.2);f.tight_layout();f.savefig(o,facecolor='white');plt.close(f)

