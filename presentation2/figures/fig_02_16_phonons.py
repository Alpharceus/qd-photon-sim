import os, sys
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 14
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from fsim_core.qd_gf import PhononParams, ibm_spectrum, zpl_weight

out = os.path.join(os.path.dirname(__file__), 'out', '02_16_phonons.png'); os.makedirs(os.path.dirname(out), exist_ok=True)
x = np.linspace(-8, 8, 1201); p = PhononParams(); fig, ax = plt.subplots(figsize=(1600/150, 900/150), dpi=150)
for T, c in [(4, '#1769aa'), (80, '#e08e0b'), (300, '#c43d3d')]:
    y = ibm_spectrum(x, p, T, .15); ax.plot(x, y/np.trapezoid(y, x), lw=3, color=c, label=f'{T} K; ZPL={zpl_weight(p,T):.2f}')
ax.set(xlabel='detuning (meV)', ylabel='normalized spectrum', title='Independent-boson spectrum: ZPL plus acoustic sideband')
ax.legend(frameon=False, fontsize=16); ax.grid(alpha=.2); fig.tight_layout(); fig.savefig(out, facecolor='white'); plt.close(fig)
