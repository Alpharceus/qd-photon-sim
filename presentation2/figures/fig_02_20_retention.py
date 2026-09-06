import os, sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from fsim_core import dot_levels

out = os.path.join(os.path.dirname(__file__), 'out', '02_20_retention.png'); os.makedirs(os.path.dirname(out), exist_ok=True)
T = np.linspace(4, 300, 300); fig, ax = plt.subplots(figsize=(16, 9), dpi=150)
for k, c in [('InP/GaAsP0.4/AlGaAs0.4 on GaAs', '#1769aa'), ('InP/GaInP/AlGaInP0.55 on GaAs', '#c43d3d')]:
    lv = dot_levels.levels(dot_levels.class_presets()[k]); p = dot_levels.retention_params(lv, 1, verbose=False)
    S = 1/(1+p['a_esc']*np.exp(-p['E_a']/(dot_levels.KB_MEV*T))+p['b_p']*np.exp(-p['E_b']/(dot_levels.KB_MEV*T)))
    ax.semilogy(T, S, lw=3, color=c, label=f'{k}: $E_a$={p["E_a"]:.0f} meV')
ax.set(xlabel='temperature (K)', ylabel='retention $S(T)$', title='Thermal escape converts confinement into brightness loss')
ax.legend(frameon=False, fontsize=14); ax.grid(alpha=.2, which='both'); fig.tight_layout(); fig.savefig(out, facecolor='white'); plt.close(fig)
