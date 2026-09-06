"""Transport checks against textbook and published semiconductor references."""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core import materials as M  # noqa: E402
from fsim_core.transport import (  # noqa: E402
    Q_SI, binary, evaluate_injection, fermi_levels, gaas_homojunction,
    hkust_preset, homojunction_vbi, n_c_eff, n_i, n_v,
    qcse_shift_meV, qfl_suppression, red_diode_preset, to_background_channel,
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

print("== 6. Carrier conservation and sub-turn-on loading (council review 2026-09-05) ==")
# Item 1: a single physical dot cannot receive more than the total captured
# supply f_QD*eta_inj*(I/q), regardless of how fractional the Poisson-
# expected N_dots under the aperture is (N_eff = max(N_dots, 1.0)).
I_test_uA = 1.0
lk_test = hk.eta_inj(300.0)
f_QD_test = hk.f_qd(1e10)
supply = f_QD_test * lk_test.eta_inj * (I_test_uA * 1e-6 / Q_SI)
cons_ok = True
for n_dots_target in (0.1, 0.377, 1.0, 10.0):
    aperture_um2 = n_dots_target / (1e10 * 1e-8)  # N_dots = n_dot_cm2 * aperture_um2 * 1e-8
    ld = hk.dot_loading(I_test_uA, 300.0, 1e10, aperture_um2, leakage=lk_test)
    cons_ok = cons_ok and (ld.N_dots > 0) and (ld.r_dot <= supply * (1.0 + 1e-9))
bool_check("r_dot never exceeds f_QD*eta_inj*(I/q) for N_dots in {0.1, 0.377, 1, 10}",
           cons_ok, "carrier conservation; N_eff = max(N_dots, 1.0) [DR]")
ld_frac = hk.dot_loading(I_test_uA, 300.0, 1e10, 0.377 / 100.0, leakage=lk_test)
ld_one = hk.dot_loading(I_test_uA, 300.0, 1e10, 1.0 / 100.0, leakage=lk_test)
bool_check("a fractional N_dots (0.377) clamps to the SAME per-dot share as N_dots=1 (N_eff=max(N_dots,1))",
           ld_frac.r_dot == ld_one.r_dot,
           "council review item 1: N_dots is a Poisson aperture-occupancy statistic, not a divisor < 1")

# Item 2: sub-turn-on Boltzmann-tail loading suppression f_qfl =
# min(1, exp(-(E_X - V_j)/kT)) -- E_X_eV=None (every check above) keeps
# f_qfl = 1 exactly (legacy numerics unchanged).
V_j_ref = hk.vj_of_j(I_test_uA * 1e-6 / hk.area_cm2, 300.0)
kT_eV = M.KB_EV * 300.0
ld_below = hk.dot_loading(I_test_uA, 300.0, 1e10, 0.785, leakage=lk_test,
                          E_X_eV=V_j_ref + 0.3)
check("f_qfl at E_X - qV_j = 0.3 eV, 300 K", ld_below.f_qfl, np.exp(-11.6), 0.02, rel=True,
      note="Boltzmann tail of the carrier reservoirs, Sze & Ng 3rd ed. ch. 12")
ld_above = hk.dot_loading(I_test_uA, 300.0, 1e10, 0.785, leakage=lk_test,
                          E_X_eV=V_j_ref - 0.1)
check("f_qfl = 1 exactly when qV_j >= E_X", ld_above.f_qfl, 1.0, 1e-12,
      note="min(1, exp(...)) clip; loading unsuppressed above turn-on")

print("== 7. Council review 2026-09-05, round 2 (background suppression, budget, area, eta_capture) ==")
# Item 1: the WL/matrix reservoir sits dE_WL ABOVE E_X and needs its OWN,
# generally much stronger, suppression than the dot's own f_qfl -- a one-
# sided bug (rate_bg unsuppressed while rate_x carried f_qfl) drove the
# 300 K sweep's pulsed g2 minimum artefact.
dE_WL_meV = 100.0
E_X_sub = V_j_ref + 0.05
bg_sub = hk.b_e(I_test_uA, 300.0, w_meV=1.0, dE_WL_meV=dE_WL_meV, n_dot_cm2=1e10,
                aperture_um2=0.1257, leakage=lk_test, E_X_eV=E_X_sub)
ld_sub = hk.dot_loading(I_test_uA, 300.0, 1e10, 0.1257, leakage=lk_test, E_X_eV=E_X_sub)
check("background suppression f_qfl_bg matches the closed form at E_X + dE_WL",
      bg_sub.f_qfl_bg,
      qfl_suppression(E_X_sub + dE_WL_meV * 1e-3, ld_sub.V_j, M.KB_EV * 300.0),
      1e-9, rel=True, note="council item 1: background suppressed at its OWN, higher, energy")
bool_check("background suppression is strictly stronger than the dot's own f_qfl when sub-turn-on",
           bg_sub.f_qfl_bg < ld_sub.f_qfl,
           "council item 1: E_X + dE_WL sits further above qV_j than E_X alone")
bg_legacy = hk.b_e(I_test_uA, 300.0, w_meV=1.0, dE_WL_meV=dE_WL_meV, n_dot_cm2=1e10,
                   aperture_um2=0.1257, leakage=lk_test)
bool_check("E_X_eV=None leaves the background suppression at 1.0 (legacy numerics unchanged)",
           bg_legacy.f_qfl_bg == 1.0, "council item 1: legacy default path untouched")

# Item 2: the carrier budget closes to eta_inj*I/q EXACTLY (r_captured +
# r_matrix), at three currents spanning sub-turn-on to near-saturated -- the
# (1-f_qfl) share that dot_loading used to just multiply away no longer
# vanishes; it is routed into r_matrix instead.
for I_budget_uA in (0.001, 0.05, 1.0):
    ld_b = hk.dot_loading(I_budget_uA, 300.0, 1e10, 0.1257, leakage=lk_test, E_X_eV=E_X_sub)
    supply_budget = lk_test.eta_inj * (I_budget_uA * 1e-6 / Q_SI)
    check(f"carrier budget r_captured+r_matrix == eta_inj*I/q at I={I_budget_uA} uA",
          ld_b.r_captured + ld_b.r_matrix, supply_budget, 1e-9, rel=True,
          note="council item 2: sub-turn-on-suppressed carriers routed, not discarded")

# Item 3 (area consistency): at fixed current, a SMALLER injection area
# gives a HIGHER current density and therefore a HIGHER V_j (Sze & Ng ideal-
# diode J(V) is monotone increasing) -- the mechanism device.py's aperture-
# area default leans on to fix V_j/N_dots/per-dot-share area consistency.
small_area = hkust_preset(area_um2=0.1257)
big_area = hkust_preset(area_um2=0.785)
V_j_small = small_area.vj_of_j(I_test_uA * 1e-6 / small_area.area_cm2, 300.0)
V_j_big = big_area.vj_of_j(I_test_uA * 1e-6 / big_area.area_cm2, 300.0)
bool_check("a smaller injection area gives a higher V_j at fixed current (area consistency)",
           V_j_small > V_j_big,
           "council item 3: J = I/area, ideal-diode J(V) monotone increasing")

# Item 7: eta_capture_dot and mu (via r_dot) must carry the SAME f_qfl
# suppression -- drive.F_p != 1 must not feed f8b_thin_fano an inconsistent
# thinning probability.
expected_eta_capture = ld_sub.f_QD * lk_test.eta_inj * ld_sub.f_qfl / max(
    1e10 * 0.1257 * 1e-8, 1.0)
check("eta_capture_dot includes the same f_qfl suppression as r_dot/mu",
      ld_sub.eta_capture_dot, expected_eta_capture, 1e-9, rel=True,
      note="council item 7: eta_capture_dot must not omit f_qfl while mu/r_dot include it")

n_pass = sum(RESULTS)
print(f"\n{n_pass}/{len(RESULTS)} transport checks passed")
sys.exit(0 if n_pass == len(RESULTS) else 1)
