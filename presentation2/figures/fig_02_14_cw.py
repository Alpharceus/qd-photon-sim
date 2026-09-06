import os, sys
from pathlib import Path
import numpy as np
import matplotlib.pyplot as plt
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..', '..')))
from fsim_core import cw_g2
from fsim_core.device import DeviceDesign, evaluate, _confinement_params

out = os.path.join(os.path.dirname(__file__), 'out', '02_14_cw.png')
os.makedirs(os.path.dirname(out), exist_ok=True)
t = np.linspace(-8, 8, 1001)
fig, ax = plt.subplots(figsize=(1600/150, 900/150), dpi=150)
card = Path(__file__).resolve().parents[2] / 'cards' / 'edge-inp-gainp-design.yaml'
for T, c in [(80, '#1769aa'), (300, '#c43d3d')]:
    design = DeviceDesign.load(card)
    design.thermal.T_hs = float(T)
    sc = evaluate(design, T_grid=[float(T)])['scalars']
    params = _confinement_params(design.ret, float(T))
    gamma = float(sc['cw_gamma_X_ns'])
    kx, kxx = cw_g2.escape_rates_from_retention(
        gamma, params['a_esc'], params['E_a'], params['b_p'], params['E_b'], float(T),
        gamma_XX_ns=2.0 * gamma)
    q = cw_g2.cw_report(
        float(sc['cw_r_ns']), gamma, 2.0 * gamma, kx, kxx, 1.0,
        float(sc['eps_op']), float(sc['cw_rho_op']), irf_fwhm_ps=0,
        tau_max_ns=8, n_tau=1001)
    ax.plot(q['curves']['tau'], q['curves']['g2_meas'], lw=3, color=c, label=f'{T} K gainp card')
ax.axhline(1, ls='--', color='.45'); ax.set(xlabel='delay $\\tau$ (ns)', ylabel='$g^{(2)}(\\tau)$', title='CW regression on the cap-2 rate-equation ladder')
ax.legend(frameon=False, fontsize=16); ax.grid(alpha=.2); fig.tight_layout(); fig.savefig(out, facecolor='white'); plt.close(fig)
