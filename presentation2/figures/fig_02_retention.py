import os,sys,numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0,os.path.abspath(os.path.join(os.path.dirname(__file__),'..','..')))
from fsim_core import dot_levels
o=os.path.join(os.path.dirname(__file__),'out','02_20_retention.png');os.makedirs(os.path.dirname(o),exist_ok=True)
T=np.linspace(4,300,300);f,a=plt.subplots(figsize=(16,9),dpi=100)
for k,c in [('InP/GaAsP0.4/AlGaAs0.4 on GaAs','#1769aa'),('InP/GaInP/AlGaInP0.55 on GaAs','#c43d3d')]:
 lv=dot_levels.levels(dot_levels.class_presets()[k]);p=dot_levels.retention_params(lv,1,verbose=False);S=1/(1+p['a_esc']*np.exp(-p['E_a']/(dot_levels.KB_MEV*T))+p['b_p']*np.exp(-p['E_b']/(dot_levels.KB_MEV*T)));a.semilogy(T,S,lw=3,c=c,label=f'{k}: Ea={p["E_a"]:.0f} meV')
a.set(xlabel='temperature (K)',ylabel='retention S(T)',title='Thermal escape converts confinement into usable brightness');a.legend(frameon=False,fontsize=14);a.grid(alpha=.2,which='both');f.tight_layout();f.savefig(o,facecolor='white');plt.close(f)
