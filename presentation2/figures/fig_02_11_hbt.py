import os, sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from fsim_core.cw_g2 import g2_two_level

out = os.path.join(os.path.dirname(__file__), 'out', '02_11_hbt.png')
os.makedirs(os.path.dirname(out), exist_ok=True)
t = np.linspace(-8, 8, 1201)
g = g2_two_level(t, 0.05, 1.0)
fig, ax = plt.subplots(figsize=(16, 9), dpi=150)
ax.plot(t, g, lw=3, color='#1769aa', label='single-emitter antibunching')
ax.axhline(1, ls='--', color='.45', label='Poisson reference')
ax.set(xlabel='delay $\\tau$ (ns)', ylabel='$g^{(2)}(\\tau)$', title='HBT: split, delay, and count coincidences')
ax.legend(frameon=False, fontsize=16); ax.grid(alpha=.2); fig.tight_layout()
fig.savefig(out, facecolor='white'); plt.close(fig)
