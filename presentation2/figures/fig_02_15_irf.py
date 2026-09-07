import os, sys
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams['font.size'] = 14
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from fsim_core.cw_g2 import convolve_irf

out = os.path.join(os.path.dirname(__file__), 'out', '02_15_irf.png')
os.makedirs(os.path.dirname(out), exist_ok=True)
t = np.linspace(-.3, .3, 6001); raw = 1 - .95*np.exp(-np.abs(t)/.005)
fig, ax = plt.subplots(figsize=(1600/150, 900/150), dpi=150)
for w, c in [(0, '#1769aa'), (20, '#e08e0b'), (100, '#c43d3d')]:
    if w == 0: y = raw
    else:
        # w is detector FWHM in ps; convolve the complete trace on the
        # uniform delay grid so recovery to one remains visible.
        y = convolve_irf(t, raw, w, shape="gaussian")
    ax.plot(t*1000, y, lw=3, color=c, label=f'{w} ps IRF')
ax.set(xlabel='delay (ps)', ylabel='$g^{(2)}$', title='A 5 ps antibunching dip is hidden by a 100 ps IRF')
ax.legend(frameon=False, fontsize=16); ax.grid(alpha=.2); fig.tight_layout(); fig.savefig(out, facecolor='white'); plt.close(fig)
