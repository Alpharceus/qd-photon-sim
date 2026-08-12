"""Drive-mechanism library regression suite (fsim_core.drive_mech, tier plan
Phase T1). Same @check style/registry as verify_sde.py / verify_fsim.py.

Covers: the DriveInterface contract, mechanism blocks M-1..M-5a, SET
feasibility pricing against the F9 wall, the fano_regulated [E] limits, and
the WP-M2' re-excitation tier -- the exact counting-moment ODE solution
cross-checked against an independent Gillespie Monte-Carlo second method
(the Yuki standard), plus the thinning-invariance claim the module docstring
makes for filter transmission.

Run: python verify/verify_drive.py   (exit code 0 iff all pass)
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core.drive_mech import (
    R_Q_OHM,
    DriveInterface,
    fano_regulated,
    mc_reexc_g2,
    mech_poisson_rail,
    mech_pulsed,
    mech_quiet_rail,
    mech_rti,
    mech_set,
    reexc_g2,
    roster,
    set_feasibility,
)
from fsim_core.loading import c_dep_for_N, f8_g2_load, fano_pump, granularity_N

CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


# ============================================================ mechanism blocks

@check("M-1 Poisson rail: F_p = 1 exactly, g2_load = 1 identically, any mu")
def _():
    for mu in (0.05, 0.4, 1.0):
        m = mech_poisson_rail(mu)
        assert m.F_p == 1.0
        assert abs(m.g2_load() - 1.0) < 1e-15, m.g2_load()


@check("M-2 quiet rail: F_p equals loading.fano_pump exactly (one "
       "implementation rule -- imported, not re-derived) and g2_load < 1")
def _():
    Zd, Zj = 3e3, 7e3          # F_p = 0.7; mu = 0.5 sits inside the cap-2 domain
    m = mech_quiet_rail(mu=0.5, Z_drive_ohm=Zd, Z_junction_ohm=Zj)
    assert abs(m.F_p - fano_pump(Zd, Zj)) < 1e-15
    assert abs(m.g2_load() - 0.4) < 1e-12   # 1 + (0.7-1)/0.5
    # F8b thinning washes the advantage: eta -> 0 restores g2_load -> 1
    m_thin = mech_quiet_rail(mu=0.5, Z_drive_ohm=Zd, Z_junction_ohm=Zj,
                             eta_capture=1e-6)
    assert abs(m_thin.g2_load() - 1.0) < 1e-4
    # the domain wall itself: an over-quiet rail at low mu must RAISE, not clamp
    try:
        mech_quiet_rail(mu=0.3, Z_drive_ohm=1e6, Z_junction_ohm=1e3).g2_load()
        assert False, "expected ValueError at the cap-2 domain wall"
    except ValueError:
        pass


@check("fano_regulated [E]: recovers F_p0 when f_reg >> f_sample, recovers 1 "
       "when f_reg << f_sample, monotone in between")
def _():
    F0 = 0.1
    assert abs(fano_regulated(F0, 1e12, 1e6) - F0) < 1e-5
    assert abs(fano_regulated(F0, 1e0, 1e9) - 1.0) < 1e-6
    fs = [float(fano_regulated(F0, fr, 1e8)) for fr in (1e6, 1e7, 1e8, 1e9, 1e10)]
    assert all(fs[i] > fs[i + 1] for i in range(len(fs) - 1)), fs


@check("M-3 pulsed: mu = eta * I*tau/e; RC rise time shrinks the window; "
       "window shorter than rise time -> mu = 0, never negative")
def _():
    m = mech_pulsed(I_uA=16.0, tau_on_ns=1.0, eta_capture=0.1)
    # 16 uA * 1 ns / e = 1e5 carriers approx: 16e-6*1e-9/1.602e-19
    n_c = 16e-6 * 1e-9 / 1.602176634e-19
    assert abs(m.notes["carriers_per_window"] - n_c) / n_c < 1e-12
    assert abs(m.mu - 0.1 * n_c) / (0.1 * n_c) < 1e-12
    m2 = mech_pulsed(I_uA=16.0, tau_on_ns=1.0, eta_capture=0.1,
                     R_drive_ohm=100.0, C_dep_F=3.5e-12)  # t_rise = 0.35 ns
    assert m2.mu < m.mu and abs(m2.notes["t_rise_ns"] - 0.35) < 1e-9
    m3 = mech_pulsed(I_uA=16.0, tau_on_ns=0.1, eta_capture=0.1,
                     R_drive_ohm=1e5, C_dep_F=3.5e-12)    # rise 350 ns > window
    assert m3.mu == 0.0


@check("SET pricing: reproduces the F9 room-temperature wall (N=1 island at "
       "300 K is a few nm / few aF; 20 nm island infeasible at 300 K, "
       "feasible at 4 K)")
def _():
    # the F9 wall: C for N=1 at 300 K ~ 6.2 aF (memory-anchored v1.1 number)
    C1 = c_dep_for_N(1.0, 300.0)
    assert 5e-18 < C1 < 8e-18, C1
    f300 = set_feasibility(300.0, radius_nm=20.0)
    f4 = set_feasibility(4.0, radius_nm=20.0, R_T_ohm=2e5)
    assert not f300["feasible"], f300
    assert f4["feasible"], f4
    assert f4["EC_over_kT"] > 10.0
    # granularity consistency with loading.granularity_N (imported, same number)
    assert abs(f4["N_granularity"] - granularity_N(f4["C_sigma_F"], 4.0)) < 1e-12


@check("SET pricing: R_T < R_Q flags infeasible regardless of E_C")
def _():
    f = set_feasibility(4.0, radius_nm=20.0, R_T_ohm=0.5 * R_Q_OHM)
    assert not f["feasible"]
    assert f["R_T_over_RQ"] < 1.0


@check("M-4 turnstile: mu = 1 by construction; g2_load = F_p = eps_cycle "
       "exactly (f8_g2_load(1, F) = F); eps_cycle -> 0 gives g2_load -> 0")
def _():
    m = mech_set("turnstile", T_K=4.0, radius_nm=20.0, R_T_ohm=2e5,
                 eps_cycle=0.02)
    assert m.mu == 1.0 and m.feasible
    assert abs(m.g2_load() - 0.02) < 1e-12
    assert abs(float(f8_g2_load(1.0, 0.02)) - 0.02) < 1e-15
    m0 = mech_set("turnstile", T_K=4.0, radius_nm=20.0, R_T_ohm=2e5,
                  eps_cycle=1e-6)
    assert m0.g2_load() < 1e-5


@check("M-4 infeasible SET degrades honestly to F_p = 1 (no free quiet "
       "statistics from hardware that cannot exist)")
def _():
    m = mech_set("metallic", T_K=300.0, radius_nm=50.0, eps_cycle=0.01)
    assert not m.feasible and m.F_p == 1.0
    assert abs(m.g2_load() - 1.0) < 1e-12


@check("M-5a RTI: boosts eta (capped at 1), carries dg_inj reduction factor")
def _():
    m = mech_rti(mu=0.4, eta_capture_base=0.05, eta_boost=3.0)
    assert abs(m.eta_capture - 0.15) < 1e-15
    assert mech_rti(mu=0.4, eta_capture_base=0.6, eta_boost=3.0).eta_capture == 1.0
    assert 0.0 < m.notes["dg_inj_factor"] < 1.0


@check("roster: every mechanism yields a valid interface (F_p in [0, 2], "
       "mu >= 0, eta in (0, 1], g2_load computable on feasible entries)")
def _():
    rs = roster()
    assert len(rs) >= 7
    for m in rs:
        assert isinstance(m, DriveInterface)
        assert 0.0 <= m.F_p <= 2.0 and m.mu >= 0.0
        assert 0.0 < m.eta_capture <= 1.0
        if (m.feasible and m.mu > 0
                and m.mu >= 1.0 - m.F_p and m.F_p <= 2.0 - m.mu):
            m.g2_load()


# ========================================================= WP-M2' re-excitation

@check("WP-M2' limits: weak-pump g2 plateaus (pump-power-independent: p and "
       "p/10 agree to <1%); g2 -> 0 as the window shrinks (the pulsed-source "
       "design rule); loaded dot + zero window -> exactly one photon, g2 = 0")
def _():
    a = reexc_g2(p_per_ps=1e-6, tau_on_ps=1000.0, tau_x_ps=1000.0)["g2"]
    b = reexc_g2(p_per_ps=1e-7, tau_on_ps=1000.0, tau_x_ps=1000.0)["g2"]
    assert a > 0.01 and abs(a - b) / a < 0.01, (a, b)
    short = reexc_g2(p_per_ps=1e-6, tau_on_ps=5.0, tau_x_ps=1000.0)["g2"]
    assert short < 0.02 * a, (short, a)
    r2 = reexc_g2(p_per_ps=0.0, tau_on_ps=0.0, tau_x_ps=1000.0, mu0=1.0,
                  t_end_factor=40.0)   # exp(-40) undecayed tail ~ 4e-18
    assert abs(r2["mean_photons"] - 1.0) < 1e-8
    assert r2["g2"] < 1e-10


@check("WP-M2' limits: hard long drive (p*tau, gX*tau >> 1) approaches "
       "Poisson-like g2 ~ 1; g2 rises monotonically with window length")
def _():
    hard = reexc_g2(p_per_ps=0.01, tau_on_ps=2e5, tau_x_ps=1000.0)
    assert 0.7 < hard["g2"] < 1.1, hard
    g2s = [reexc_g2(p_per_ps=0.002, tau_on_ps=t, tau_x_ps=1000.0)["g2"]
           for t in (100.0, 1000.0, 5000.0, 2e4)]
    assert all(g2s[i] < g2s[i + 1] for i in range(len(g2s) - 1)), g2s


@check("WP-M2' second method: moment-ODE g2 agrees with Gillespie MC within "
       "4 sigma at three parameter points (the Yuki standard)")
def _():
    pts = [(0.002, 1500.0, 1000.0, 0.0),
           (0.005, 3000.0, 800.0, 0.0),
           (0.001, 2000.0, 1000.0, 1.0)]
    for p, ton, tx, mu0 in pts:
        exact = reexc_g2(p_per_ps=p, tau_on_ps=ton, tau_x_ps=tx, mu0=mu0)["g2"]
        mc, se = mc_reexc_g2(p_per_ps=p, tau_on_ps=ton, tau_x_ps=tx, mu0=mu0,
                             n_pulses=150_000, seed=7)
        assert se > 0 and abs(mc - exact) < 4.0 * se, (p, ton, exact, mc, se)


@check("WP-M2' thinning invariance: filter transmission t_filter does not "
       "change g2 (binomial thinning preserves <m(m-1)>/<m>^2) -- MC at "
       "t=1.0 vs t=0.3 agree within combined 4 sigma")
def _():
    a, sa = mc_reexc_g2(0.003, 2000.0, 1000.0, t_filter=1.0,
                        n_pulses=150_000, seed=11)
    b, sb = mc_reexc_g2(0.003, 2000.0, 1000.0, t_filter=0.3,
                        n_pulses=150_000, seed=12)
    assert abs(a - b) < 4.0 * np.hypot(sa, sb), (a, b, sa, sb)


@check("WP-M2' cascade bookkeeping: XX photons are not counted -- with the "
       "pump confined to a short strong window (near-deterministic double "
       "load) mean X photons stays <= 1 + re-excitation, never ~2")
def _():
    # strong short window loads the dot up to XX; the cascade emits ONE X
    r = reexc_g2(p_per_ps=1.0, tau_on_ps=10.0, tau_x_ps=1000.0)
    assert r["mean_photons"] < 1.15, r
    assert r["mean_photons"] > 0.9


@check("card factory: mech_from_card resolves every roster name; unknown "
       "name raises; params override context defaults")
def _():
    from fsim_core.drive_mech import mech_from_card
    for name in ("poisson-rail", "quiet-rail", "pulsed", "set-metallic",
                 "set-gated", "set-turnstile", "rti", "m1", "m4c"):
        m = mech_from_card(name, T_K=77.0, mu=0.9)
        assert isinstance(m, DriveInterface), name
    assert mech_from_card("m1", {"mu": 0.33}).mu == 0.33
    try:
        mech_from_card("warp-drive")
        assert False, "expected ValueError on unknown mechanism"
    except ValueError:
        pass


@check("device wiring: mechanism='' keeps evaluate() BIT-IDENTICAL to legacy "
       "(staged preset, full curve compare); mechanism='poisson-rail' with "
       "matching mu also reproduces it exactly (F_p=1 branch)")
def _():
    from fsim_core.device import DeviceDesign
    from fsim_core.device import evaluate as dev_eval
    base = DeviceDesign()
    ref = dev_eval(base)
    d2 = DeviceDesign()
    assert d2.drive.mechanism == ""          # default stays legacy
    same = dev_eval(d2)
    assert np.array_equal(ref["curves"]["g2"], same["curves"]["g2"],
                          equal_nan=True)
    d3 = DeviceDesign()
    d3.drive.mechanism = "poisson-rail"      # F_p=1 -> same f1b path
    m3 = dev_eval(d3)
    assert np.array_equal(ref["curves"]["g2"], m3["curves"]["g2"],
                          equal_nan=True)


@check("device wiring: a feasible cryogenic turnstile card LOWERS g2 vs the "
       "Poisson rail at the same design (g2_load = eps_cycle multiplier); an "
       "infeasible RT SET card degrades to the legacy curve, never below")
def _():
    from fsim_core.device import DeviceDesign
    from fsim_core.device import evaluate as dev_eval
    base = DeviceDesign()
    ref = dev_eval(base)["scalars"] if "scalars" in dev_eval(base) else None
    lo = DeviceDesign()
    lo.thermal.T_hs = 10.0
    ref_op = dev_eval(lo)["curves"]
    turn = DeviceDesign()
    turn.thermal.T_hs = 10.0
    turn.drive.mechanism = "set-turnstile"
    turn.drive.mech_params = {"radius_nm": 10.0, "R_T_ohm": 2e5,
                              "eps_cycle": 0.02}
    t_op = dev_eval(turn)["curves"]
    i = 0  # the coldest grid point (4 K): E_C/kT clears the margin-10 bar
    # note: at ~10 K a 15 nm island already FAILS the margin and honestly
    # degrades to F_p=1 -- that near-miss is the pricing doing its job
    assert t_op["g2"][i] < ref_op["g2"][i], (t_op["g2"][i], ref_op["g2"][i])
    # infeasible at RT: the same island priced at a 300 K junction gives F_p=1
    hot = DeviceDesign()
    hot.drive.mechanism = "set-metallic"
    hot.drive.mech_params = {"radius_nm": 50.0, "eps_cycle": 0.02}
    h = dev_eval(hot)["curves"]
    legacy = dev_eval(DeviceDesign())["curves"]
    hi = -5  # a hot grid point
    assert abs(h["g2"][hi] - legacy["g2"][hi]) < 1e-12


@check("hygiene: fsim_core.drive_mech imports no GUI framework and no "
       "duplicate F8/F9 closed forms (must import them from loading)")
def _():
    src = (Path(__file__).resolve().parents[1] / "fsim_core" / "drive_mech.py"
           ).read_text()
    for banned in ("tkinter", "PyQt", "PySide", "wx", "dearpygui", "streamlit"):
        assert banned not in src, banned
    # one-implementation rule: the F8/F9 formulas must be imported, not typed
    assert "from .loading import" in src
    for forbidden in ("1.0 + (F_p - 1.0) / mu", "Z_junction / (Z_junction"):
        assert forbidden not in src, f"duplicate closed form: {forbidden}"


def main():
    failed = 0
    for name, fn in CHECKS:
        try:
            fn()
            print(f"  PASS  {name}")
        except AssertionError as e:
            failed += 1
            print(f"* FAIL  {name}  {e}")
    n = len(CHECKS)
    print(f"\n{n - failed}/{n} checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
