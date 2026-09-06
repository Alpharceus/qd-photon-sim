import os, sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from fsim_core.cw_g2 import dip_convolved

out = os.path.join(os.path.dirname(__file__), 'out', '02_15_irf.png')
os.makedirs(os.path.dirname(out), exist_ok=True)
t = np.linspace(-.08, .08, 1601); raw = 1 - .95*np.exp(-np.abs(t)/.005)
fig, ax = plt.subplots(figsize=(16, 9), dpi=150)
for w, c in [(0, '#1769aa'), (20, '#e08e0b'), (100, '#c43d3d')]:
    if w == 0: y = raw
    else:
        sigma = w*1e-3/(2*np.sqrt(2*np.log(2)))
        y = np.full_like(t, dip_convolved(.95, .005, sigma))
    ax.plot(t*1000, y, lw=3, color=c, label=f'{w} ps IRF')
ax.set(xlabel='delay (ps)', ylabel='$g^{(2)}$', title='A 5 ps antibunching dip is hidden by a 100 ps IRF')
ax.legend(frameon=False, fontsize=16); ax.grid(alpha=.2); fig.tight_layout(); fig.savefig(out, facecolor='white'); plt.close(fig)
