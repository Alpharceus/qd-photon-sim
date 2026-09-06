import os, sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from fsim_core import cw_g2

out = os.path.join(os.path.dirname(__file__), 'out', '02_14_cw.png')
os.makedirs(os.path.dirname(out), exist_ok=True)
t = np.linspace(-8, 8, 1001)
fig, ax = plt.subplots(figsize=(16, 9), dpi=150)
for T, c, rho in [(80, '#1769aa', .995), (300, '#c43d3d', .90)]:
    q = cw_g2.cw_report(.05, 1, 2, 0, 0, .5, .2, rho, irf_fwhm_ps=0, tau_max_ns=8, n_tau=1001)
    ax.plot(t, q['curves']['g2_meas'], lw=3, color=c, label=f'{T} K proxy')
ax.axhline(1, ls='--', color='.45'); ax.set(xlabel='delay $\\tau$ (ns)', ylabel='$g^{(2)}(\\tau)$', title='CW regression on the cap-2 rate-equation ladder')
ax.legend(frameon=False, fontsize=16); ax.grid(alpha=.2); fig.tight_layout(); fig.savefig(out, facecolor='white'); plt.close(fig)
