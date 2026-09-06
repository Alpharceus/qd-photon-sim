import os,sys,numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0,os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..')))
from fsim_core.linewidth import gamma_envelope
o=os.path.join(os.path.dirname(__file__),'out','02_17_linewidth.png');os.makedirs(os.path.dirname(o),exist_ok=True)
T=np.linspace(4,300,300);lo,hi=gamma_envelope(T);f,a=plt.subplots(figsize=(16,9),dpi=100);a.fill_between(T,lo,hi,color='#c43d3d',alpha=.25);a.plot(T,lo,lw=3,label='low edge');a.plot(T,hi,lw=3,label='high edge');a.set(xlabel='temperature (K)',ylabel='FWHM Γ (meV)',title='Thermal broadening moves the linewidth into the meV range');a.legend(frameon=False,fontsize=16);a.grid(alpha=.2);f.tight_layout();f.savefig(o,facecolor='white');plt.close(f)

