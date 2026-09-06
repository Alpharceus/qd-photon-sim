import os, sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from fsim_core.spectral import epsilon

out = os.path.join(os.path.dirname(__file__), 'out', '02_19_spectral.png'); os.makedirs(os.path.dirname(out), exist_ok=True)
G = np.linspace(.2, 20, 250); fig, ax = plt.subplots(figsize=(16, 9), dpi=150)
for d, c in zip([4, 6, 8], ['#1769aa', '#e08e0b', '#c43d3d']):
    ax.plot(G, [epsilon(d, g, g, w=g).eps for g in G], lw=3, color=c, label=f'$\\Delta_{{XX}}={d}$ meV')
ax.set(xlabel='FWHM $\\Gamma$ (meV)', ylabel='$\\varepsilon=t_{XX}/t_X$', title='Spectral overlap increases as thermal linewidth grows')
ax.legend(frameon=False, fontsize=16); ax.grid(alpha=.2); fig.tight_layout(); fig.savefig(out, facecolor='white'); plt.close(fig)
