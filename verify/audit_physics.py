"""Physics audit: known parameters in -> known results out.

Unlike verify_fsim.py (internal second methods), every item here checks the
implementation against an EXTERNAL known number: the foundation doc's published
check values, the F7 table, published Varshni shifts, the Reischle
background-corrected values, and closed-loop extractions using measured (not
fitted) linewidths. Discrepancies are FLAGGED, not hidden.

Run: python verify/audit_physics.py
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core.cavity import VARSHNI, varshni_shift
from fsim_core.integrator import retention
from fsim_core.loading import aperture_g2, f1b_g2, loading_probs
from fsim_core.spectral import epsilon, epsilon_narrow_filter, tophat_transmission

ROOT = Path(__file__).resolve().parents[1]
FITP = json.loads((ROOT / "out" / "phase0" / "fit_params.json").read_text())["params"]

results = []


def item(name, computed, known, tol, unit="", note=""):
    ok = abs(computed - known) <= tol
    results.append(ok)
    mark = "PASS " if ok else "FLAG*"
    print(f"  {mark} {name}")
    print(f"        computed {computed:.4g}{unit}  vs known {known:.4g}{unit}"
          f"  (tol {tol:g}){('  -- ' + note) if note else ''}")
    return ok


print("== 1. Foundation F1b check triplet (Delta=5.4, Gamma=6.5, w_half=Gamma/2) ==")
# foundation convention: half-width w; our code: full width -> w_full = Gamma
eps = epsilon(5.4, 6.5, 6.5, w=6.5).eps
item("epsilon at the foundation check point", eps, 0.3993, 0.002)
for mu, known in ((0.05, 0.3979), (0.3, 0.3929), (1.0, 0.3877)):
    item(f"g2(mu={mu}) exact closed form", float(f1b_g2(mu, eps)), known, 0.0015)

print("\n== 2. F1a convention equivalence (half-width doc form vs full-width code) ==")
D, G, wh = 5.4, 6.5, 3.25
t_x_doc = (2 / np.pi) * np.arctan(2 * wh / G)
t_xx_doc = (1 / np.pi) * (np.arctan(2 * (D + wh) / G) - np.arctan(2 * (D - wh) / G))
item("t_X doc form vs code", tophat_transmission(0.0, G, 2 * wh), t_x_doc, 1e-12)
item("t_XX doc form vs code", tophat_transmission(D, G, 2 * wh), t_xx_doc, 1e-12)

print("\n== 3. F7 requirement table (Gamma(300K) = 6.5-7 meV [V]) ==")
for D, elo_known, ehi_known in ((1.8, 0.77, 0.79), (5.4, 0.27, 0.30), (8.0, 0.14, 0.16),
                                (30.0, 0.0, 0.02)):
    e65 = epsilon_narrow_filter(D, 6.5, 6.5)
    e70 = epsilon_narrow_filter(D, 7.0, 7.0)
    item(f"eps_min(Delta={D}) low edge", min(e65, e70), elo_known, 0.015)
    item(f"eps_min(Delta={D}) high edge", max(e65, e70), ehi_known, 0.015)
rho_req = lambda e: np.sqrt(0.5 / (1 - e))
r_lo = rho_req(epsilon_narrow_filter(5.4, 6.5, 6.5))
r_hi = rho_req(epsilon_narrow_filter(5.4, 7.0, 7.0))
print(f"  NOTE  rho_req(Delta=5.4) recomputed: {r_lo:.3f}-{r_hi:.3f}; the foundation "
      "table states '>=0.84-0.86'.")
print("        The table's upper edge appears ~0.015-0.02 high (0.86 needs eps=0.32, "
      "outside its own eps row 0.27-0.30). Minor erratum in the doc table; "
      "the eps row and every downstream map use the recomputed values.")

print("\n== 4. Theorem-0 and the two-emitter threshold ==")
item("cascade saturates at the pair-source value", float(f1b_g2(500.0, 1.0)), 0.5, 1e-4)
item("two equal emitters g2 (OE2008's quoted 0.5 threshold)",
     aperture_g2([1.0, 1.0]), 0.5, 1e-12)

print("\n== 5. Reischle OE2008: F2 inversion vs the paper's published corrections ==")
for dot, g2b, rho, known in (("B/4K", 0.15, 0.94, 0.04), ("C/80K", 0.25, 0.88, 0.03)):
    g2s = (g2b - (1 - rho**2)) / rho**2
    item(f"dot {dot} background-corrected g2", g2s, known, 0.011,
         note="paper's Eq.(1) == F2 law")

print("\n== 6. Varshni (GaAs) vs digitized Chatzarakis X positions (NOT fitted) ==")
# measured X energies digitized from Fig. 4: 78/120/170/230 K
meas = {78: 1.3310, 120: 1.3200, 170: 1.3040, 230: 1.2800}
E78 = meas[78]
worst = 0.0
for T, E in meas.items():
    pred = E78 + float(varshni_shift(T, "GaAs") - varshni_shift(78, "GaAs"))
    worst = max(worst, abs(pred - E) * 1e3)
    print(f"        X({T:3d} K): predicted {pred:.4f} eV  vs digitized {E:.4f} eV "
          f"({(pred - E) * 1e3:+.1f} meV)")
item("worst Varshni-GaAs residual over 78-230 K", worst, 0.0, 4.0, unit=" meV",
     note="dot redshifts ~3 meV/150 K slower than bulk GaAs ('practically "
          "identical' per paper; ~6%)")

print("\n== 7. Varshni (InP) vs Laferriere published wavelengths ==")
# published: 1301.28 nm at 4 K -> 1397.8 nm at 300 K
dE_meas = 1239.842 / 1301.28 - 1239.842 / 1397.8  # eV redshift
dE_pred = -float(varshni_shift(300.0, "InP") - varshni_shift(4.0, "InP"))
item("4->300 K redshift", dE_pred * 1e3, dE_meas * 1e3, 8.0, unit=" meV",
     note="InAsP dot in InP nanowire vs bulk-InP Varshni; 7% class agreement")

print("\n== 8. Chatzarakis 78 K bound from MEASURED inputs (no fit parameters) ==")
# Gamma_X(78)=1.06 [S2 digitized], r_xx from S2 ratio ~0.72 at 78: 0.67/1.06
e78 = epsilon(5.9, 1.06, 0.67, w=1.40, dx=-0.37).eps
g2_78 = float(f1b_g2(0.33, e78))
print(f"        eps(78 K) = {e78:.4f} -> g2_dot(mu=0.33) = {g2_78:.4f}")
item("g2_dot(78 K) vs published bound (<= 0.02)", g2_78, 0.0, 0.02,
     note="pure cascade leakage sits well under the bound; remaining budget is rho")

print("\n== 9. Three independent routes to rho(230 K) agree ==")
# route A: the V-a joint fit
S230 = float(retention(230.0, FITP["a_esc"], FITP["E_a"], FITP["b_p"], FITP["E_b"]))
B230 = FITP["b0"] + FITP["beta"] * (1 - S230)
rho_fit = S230 / (S230 + B230)
# route B: invert F2 from the measured g2(230)=0.36 using MEASURED Gamma(230)
gam230 = np.interp(230.0, [200.0, 240.0], [5.58, 6.22])  # S2 digitized
e230 = epsilon(5.9, gam230, 0.72 * gam230, w=7.66, dx=-3.41).eps
g2dot = float(f1b_g2(0.33, e230))
rho_inv = np.sqrt((1 - 0.36) / (1 - g2dot))
# route C: lifetime retention (S4 printed) + fitted background only
S78 = float(retention(78.0, FITP["a_esc"], FITP["E_a"], FITP["b_p"], FITP["E_b"]))
S_tau = S78 * 0.21 / 1.77
B_tau = FITP["b0"] + FITP["beta"] * (1 - S_tau)
rho_tau = S_tau / (S_tau + B_tau)
print(f"        A joint fit: {rho_fit:.3f}   B measured-Gamma F2 inversion: {rho_inv:.3f}"
      f"   C lifetime route: {rho_tau:.3f}")
item("spread of the three routes", max(rho_fit, rho_inv, rho_tau)
     - min(rho_fit, rho_inv, rho_tau), 0.0, 0.02,
     note="eps(230) from measured Gamma = %.3f -> the 230 K point is ~%d%% eps, "
          "~%d%% background" % (e230, 100 * g2dot * (1 - 0) / 0.36,
                                100 - 100 * g2dot / 0.36))

print("\n== 10. Thermal magnitude vs VCSEL literature ==")
# oxide-aperture VCSEL class: R_th ~ 1-3 K/mW for ~8 um aperture on GaAs
from fsim_core.thermal import Stack, stack_rth
R = stack_rth(4e-6, Stack(layers=[], k_sub300=55.0), 300.0)
item("bare-disk R_th, 8 um aperture on GaAs", R / 1000.0, 2.0, 1.5, unit=" K/mW",
     note="published VCSEL thermal resistances ~1-3 K/mW; spreading term dominates")

print("\n== 11. Mode-tracking magnitude vs the papers' remarks ==")
dl_78_300 = 931.0 * (float(varshni_shift(300, "GaAs") - varshni_shift(78, "GaAs"))
                     / -1.3313) * -1.0
print(f"        computed emitter walk 78->300 K at 931 nm: {abs(dl_78_300):.0f} nm; "
      "paper remark: 'shifted ~50 nm -> RT achievable'; foundation doc: '35-50 nm "
      "from 20->300 K'")
print("  NOTE  computed ~59 nm (78->300 K) / ~67 nm (20->300 K) is larger than the "
      "foundation's 35-50 nm figure; the paper's own ~50 nm remark is closer. "
      "Design consequence unchanged (mode red of cryo line by the shift to "
      "T_target), but the doc's nm figure should be revised upward.")

n_pass = sum(results)
print(f"\n{n_pass}/{len(results)} quantitative audit items PASS "
      f"({len(results) - n_pass} flagged); plus 2 documented notes (F7 rho_req "
      "erratum; Varshni-walk nm figure).")
sys.exit(0 if n_pass == len(results) else 1)
