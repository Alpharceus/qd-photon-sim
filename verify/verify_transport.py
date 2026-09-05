"""Transport checks against textbook and published semiconductor references."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core import materials as M  # noqa: E402
from fsim_core.transport import (  # noqa: E402
    binary, evaluate_injection, fermi_levels, gaas_homojunction,
    hkust_preset, homojunction_vbi, n_c_eff, n_i, n_v,
    qcse_shift_meV, red_diode_preset, to_background_channel,
    xi_window, xi_window_unclipped,
)

RESULTS = []


def check(name, value, ref, tol, note="", rel=False):
    err = abs(value - ref) / (abs(ref) if rel and ref != 0 else 1.0)
    ok = err <= tol
    unit = "%" if rel else ""
    RESULTS.append(ok)
    print(f"  {'PASS' if ok else 'FAIL'}  {name}\n        computed {value:.6g} vs ref {ref:.6g} "
          f"(tol {tol * (100 if rel else 1):.3g}{unit}) -- {note}")
    return ok


def bool_check(name, ok, note=""):
    RESULTS.append(bool(ok))
    print(f"  {'PASS' if ok else 'FAIL'}  {name} -- {note}")
    return bool(ok)


gaas, inp, gap = (M.binary(x) for x in ("GaAs", "InP", "GaP"))

print("== 1. Carrier statistics (Sze; Ioffe NSM) ==")
check("GaAs n_i(300 K), cm^-3", n_i(gaas), 2.1e6, 0.3e6,
      "Sze & Ng 3rd ed.; Ioffe NSM 2.1e6 cm^-3")
check("GaAs N_c(300 K), cm^-3", n_c_eff(gaas)[0], 4.7e17, 0.7e17,
      "Ioffe NSM 4.7e17 cm^-3")
check("GaAs N_v(300 K), cm^-3", n_v(gaas), 9.0e18, 0.6e18,
      "Ioffe NSM 9e18 cm^-3")
check("InP n_i(300 K), cm^-3", n_i(inp), 1.3e7, 0.7e7,
      "Ioffe NSM 1.3e7 cm^-3")
bool_check("GaP n_i is a few cm^-3", n_i(gap) < 100.0,
           "Ioffe NSM GaP intrinsic density, order 1-10 cm^-3")

print("== 2. Built-in voltage and depletion (Sze) ==")
vbi_ref = homojunction_vbi(gaas, 1e17, 1e17)
check("GaAs homojunction V_bi", vbi_ref, 1.25, 0.05,
      "Sze & Ng 3rd ed., kT ln(N_A N_D/n_i^2), 1.20-1.30 V")
d0 = gaas_homojunction(d_i_nm=0.0)
check("GaAs homojunction W(0 V), nm", d0.depletion(0.0).W_nm, 160.0, 40.0,
      "Sze & Ng abrupt-junction depletion formula, 0.12-0.20 um")
check("gaas_homojunction V_bi", gaas_homojunction().vbi(), vbi_ref, 0.02,
      "Sze & Ng homojunction Fermi-level equality")

print("== 3. Diode I-V and injection ==")
red = red_diode_preset()
V = red.vbi() - 0.2
parts = {}
for channel, (ln_j0, ideality) in red.ln_j0().items():
    parts[channel] = np.exp(ln_j0 + V / (ideality * M.KB_EV * 300.0))
bool_check("n=1 diode current dominates n=2 at V_bi-0.2 V",
           parts.get("rad", 0.0) + parts.get("diff", 0.0) > parts.get("srh", 0.0),
           "Sze & Ng ideal diode terms; radiative/diffusion n=1 vs SRH n=2")
hk = hkust_preset()
check("HKUST-like diode V_j at 1 uA", hk.v_of_i(1e-6, 300.0)[1], 1.9, 0.4,
      "AlGaInP red-LED class turn-on, OSRAM OSLON datasheets; Chen/Kish LED papers")
eta_200, eta_300, eta_350 = (hk.eta_inj(T).eta_inj for T in (200.0, 300.0, 350.0))
bool_check("HKUST eta_inj(300 K) in the 0.3-class range", 0.25 <= eta_300 <= 1.0,
           "Bour et al., IEEE JQE 30 (1994), thermionic leakage model; class-range edge")
bool_check("eta_inj decreases from 200 K to 350 K", eta_200 > eta_300 > eta_350,
           "Bour et al. 1994; Coldren et al., carrier leakage Arrhenius dependence")

print("== 4. Dot loading and background ==")
lo = evaluate_injection(hk, 0.05, 300.0, 1e10, 0.785, 1.0, 1.0, 100.0)
hi = evaluate_injection(hk, 0.10, 300.0, 1e10, 0.785, 1.0, 1.0, 100.0)
check("mu(2I)/mu(I) at low current", hi.mu / lo.mu, 2.0, 0.10,
      "Poisson carrier budget; linear low-current injection model")
bool_check("dot loading is explicitly linear rather than falsely saturated",
           hi.loading.r_dot == 2.0 * lo.loading.r_dot,
           "transport.py dot_loading docstring; capture saturation is downstream cap-2")
kw = dict(w_meV=1.0, dE_WL_meV=100.0, n_dot_cm2=1e10, aperture_um2=0.785)
direct = hk.b_e(1.0, 300.0, **kw).b_e
channel, info = to_background_channel(hk, [1.0], [300.0], I_ref_uA=1.0, **kw)
check("background channel at (I_ref,300 K)", channel.rate(1.0, 300.0), direct, 1e-9,
      "self-consistency of transport.to_background_channel with direct b_e")
bool_check("background b_e grows with temperature", hk.b_e(0.001, 350.0, **kw).b_e >
           hk.b_e(0.001, 300.0, **kw).b_e,
           "Urbach/thermal tail; Sze & Ng recombination spectrum treatment")

print("== 5. Window, QCSE, and Fermi-level invariants ==")
widths = np.array([0.5, 1.0, 2.0, 4.0])
xis = xi_window(widths, 100.0, 25.0)
bool_check("xi_window is clipped to [0,1] and monotone in width",
           np.all((xis >= 0.0) & (xis <= 1.0)) and np.all(np.diff(xis) >= 0),
           "normalized one-sided Urbach tail; transport.py analytic integral")
bool_check("xi_window_unclipped is at least xi_window",
           np.all(xi_window_unclipped(widths, 100.0, 25.0) + 1e-15 >= xis),
           "unclipped exponential-tail derivation in transport.py")
check("QCSE shift at zero field", float(qcse_shift_meV(0.0, polarizability=1.0)), 0.0, 1e-12,
      "second-order Stark shift; standard perturbation theory")
bool_check("positive-field QCSE shift is red", qcse_shift_meV(10.0, polarizability=1.0) <= 0.0,
           "second-order Stark shift sign from perturbation theory")
efn, efp = fermi_levels(gaas, 1e17, gaas, 1e17)
midgap = (M.band_edges(gaas, 300.0)[0] + M.band_edges(gaas, 300.0)[1]) / 2.0
bool_check("n-side E_F is above midgap", efn > midgap,
           "Sze & Ng Boltzmann Fermi level for n-type doping")
bool_check("p-side E_F is below midgap", efp < midgap,
           "Sze & Ng Boltzmann Fermi level for p-type doping")

n_pass = sum(RESULTS)
print(f"\n{n_pass}/{len(RESULTS)} transport checks passed")
sys.exit(0 if n_pass == len(RESULTS) else 1)
