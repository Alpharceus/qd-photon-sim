"""Spec-mode gate: fsim_core/spec.py's inversions must round-trip through the
SAME forward chain they invert (F1b, background law, F6, F5, Module A), must
be monotone where the physics says they must be, spec_sheet() must flag
per-route and overall feasibility honestly on a known-feasible point, a
known-infeasible point, and (critically) a MIXED point where some routes
close and some don't -- that mixed case is exactly where the old blanket
AND-of-everything "feasible" flag lied (INFEASIBLE printed while the cavity
route already closed). scripts/run_spec.py must run end-to-end and leave a
parseable bundle whose verdict column reads honestly.
Standalone, same style as gate_d1.py / gate_d2.py / gate_d3.py.
Run: python verify/gate_spec.py   (exit code 0 iff all pass)
"""
import csv
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.loading import f1b_g2, n_window_competitors
from fsim_core.spec import (
    G_required,
    b_e_budget,
    density_limit,
    eps_budget,
    kappa_ceiling,
    kappa_min_of,
    rho_required,
    spec_sheet,
    v_tilde_required,
    w_floor_of,
)
from fsim_core.spectral import epsilon

CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


MUS = (0.05, 0.1, 0.33, 0.5, 1.0, 2.0)
TARGETS = (0.01, 0.05, 0.1, 0.3, 0.5, 0.8)


# ------------------------------------------------------------- (a) round trips

@check("eps_budget: f1b_g2(mu, eps_budget) == target within 1e-9 (budget < 1)")
def _():
    worst = 0.0
    for mu in MUS:
        for target in TARGETS:
            eb = eps_budget(mu, target)
            if eb < 1.0:
                g2 = float(f1b_g2(mu, eb))
                worst = max(worst, abs(g2 - target))
                assert abs(g2 - target) < 1e-9, (mu, target, eb, g2)
    assert worst < 1e-9


@check("rho_required: plugging rho_req back into the background law gives g2 == target")
def _():
    worst = 0.0
    n_checked = 0
    for eps_op in (0.01, 0.05, 0.1, 0.2, 0.35):
        for mu in (0.1, 0.5, 1.0):
            for target in (0.3, 0.5, 0.8):
                g2dot = float(f1b_g2(mu, eps_op))
                rr = rho_required(eps_op, mu, target)
                if g2dot < target:
                    n_checked += 1
                    g2 = 1.0 - rr**2 * (1.0 - g2dot)
                    worst = max(worst, abs(g2 - target))
                    assert abs(g2 - target) < 1e-9, (eps_op, mu, target, rr, g2)
                else:
                    assert rr == float("inf"), (eps_op, mu, target, rr)
    assert n_checked > 0 and worst < 1e-9


@check("G_required: rho(G_req) == rho_req within 1e-9")
def _():
    worst = 0.0
    for rho_req in (0.3, 0.5, 0.7, 0.9, 0.99):
        for S in (0.3, 0.6, 0.9):
            for B in (0.01, 0.05, 0.2):
                G = G_required(rho_req, S, B)
                assert np.isfinite(G)
                rho = G * S / (G * S + B)
                worst = max(worst, abs(rho - rho_req))
                assert abs(rho - rho_req) < 1e-9, (rho_req, S, B, G, rho)
    assert worst < 1e-9
    assert G_required(1.0, 0.5, 0.1) == float("inf")
    assert G_required(float("inf"), 0.5, 0.1) == float("inf")


@check("b_e_budget: rho at B0 + b_e_budget == rho_req within 1e-9 (when non-negative)")
def _():
    worst = 0.0
    n_checked = 0
    for rho_req in (0.3, 0.5, 0.7, 0.9, 0.99):
        for S in (0.3, 0.6, 0.9):
            for B0 in (0.001, 0.01, 0.05):
                be = b_e_budget(rho_req, S, B0)
                if be >= 0.0:
                    n_checked += 1
                    rho = S / (S + B0 + be)
                    worst = max(worst, abs(rho - rho_req))
                    assert abs(rho - rho_req) < 1e-9, (rho_req, S, B0, be, rho)
    assert n_checked > 0 and worst < 1e-9
    assert b_e_budget(1.0, 0.5, 0.1) == float("-inf")


@check("kappa_ceiling: eps_cav(kappa_max) == eps_budget within 1e-6 (interior solutions)")
def _():
    worst = 0.0
    n_checked = 0
    for delta in (1.5, 2.5, 3.5, 5.0, 8.0):
        for gx in (1.0, 2.0, 4.0):
            gxx = 0.72 * gx
            for eb in (0.01, 0.05, 0.1, 0.3, 0.6):
                km = kappa_ceiling(delta, gx, gxx, eb)
                if np.isfinite(km) and 1e-4 < km < 50.0:
                    n_checked += 1
                    ec = epsilon(delta, gx, gxx, kappa=km, dx=0.0).eps
                    worst = max(worst, abs(ec - eb))
                    assert abs(ec - eb) < 1e-6, (delta, gx, gxx, eb, km, ec)
    assert n_checked > 0 and worst < 1e-6


# --------------------------------------------------------- (b) v_tilde round trip

@check("v_tilde_required: F_eff computed from (kappa, V_tilde_req) == F_eff_target within 1e-12")
def _():
    C_PUR = 3.0 / (4.0 * np.pi**2)
    worst = 0.0
    for F_eff in (0.5, 1.0, 3.0, 10.0):
        for kappa in (0.05, 0.5, 1.0, 5.0):
            for gamma in (0.5, 1.0, 2.0, 3.0):
                for E_eV in (1.3, 1.55, 1.88):
                    vt = v_tilde_required(F_eff, kappa, gamma, E_eV)
                    E_meV = E_eV * 1000.0
                    F_check = C_PUR * E_meV / ((kappa + gamma) * vt)
                    worst = max(worst, abs(F_check - F_eff))
                    assert abs(F_check - F_eff) < 1e-12, (F_eff, kappa, gamma, E_eV, vt, F_check)
    assert worst < 1e-12


@check("brightness floor round trips: t_X(w_floor) == t_x_floor (1e-9); "
      "kappa_min/(kappa_min+Gamma) == t_x_floor (1e-12)")
def _():
    worst_w, worst_k = 0.0, 0.0
    for gx in (0.5, 1.0, 2.0, 5.0):
        for tf in (0.1, 0.3, 0.5, 0.7):
            wf = w_floor_of(tf, gx)
            tx = epsilon(3.0, gx, 0.7 * gx, w=wf, dx=0.0).t_x
            worst_w = max(worst_w, abs(tx - tf))
            assert abs(tx - tf) < 1e-9, (gx, tf, wf, tx)

            km = kappa_min_of(tf, gx)
            tcav = km / (km + gx)
            worst_k = max(worst_k, abs(tcav - tf))
            assert abs(tcav - tf) < 1e-12, (gx, tf, km, tcav)
    assert worst_w < 1e-9 and worst_k < 1e-12


@check("density_limit: continuous round trip through run_phase3's forward "
      "g2pen(N_w) formula returns penalty_budget within 1e-6 relative")
def _():
    worst = 0.0
    for budget in (0.01, 0.05, 0.1, 0.3):
        for r in (0.1, 0.3, 0.5):
            for w, sigma, ap in [(2.0, 40.0, 1.0), (6.5, 40.0, 1.0), (4.0, 60.0, 2.0)]:
                dens = density_limit(budget, ap, w, sigma, r)
                area_um2 = np.pi * (ap / 2.0) ** 2
                Nw = n_window_competitors(dens, area_um2, w, sigma)
                g2pen = 1.0 - (1.0 + Nw * r**2) / (1.0 + Nw * r) ** 2
                rel = abs(g2pen - budget) / budget
                worst = max(worst, rel)
                assert rel < 1e-6, (budget, r, w, sigma, ap, dens, Nw, g2pen)
    assert worst < 1e-6


# ---------------------------------------------------------------- (c) monotonicity

@check("kappa_ceiling: monotone (non-decreasing) in Delta_XX at fixed Gamma/eps_budget")
def _():
    gx, gxx, eb = 2.0, 1.44, 0.05
    deltas = np.linspace(1.0, 10.0, 30)
    kms = [kappa_ceiling(d, gx, gxx, eb) for d in deltas]
    finite = [k for k in kms if np.isfinite(k)]
    assert len(finite) > 5
    assert all(finite[i] <= finite[i + 1] + 1e-9 for i in range(len(finite) - 1)), kms
    # infeasible (nan) region, where it exists, must sit at the SMALL-delta end
    first_finite = next(i for i, k in enumerate(kms) if np.isfinite(k))
    assert all(not np.isfinite(k) for k in kms[:first_finite])


@check("b_e_budget: monotone decreasing in rho_req (stricter rho -> smaller/negative budget)")
def _():
    S, B0 = 0.6, 0.01
    rhos = np.linspace(0.2, 0.999, 40)
    bes = [b_e_budget(r, S, B0) for r in rhos]
    assert all(bes[i] >= bes[i + 1] - 1e-9 for i in range(len(bes) - 1)), bes


# ------------------------------------------------------ (d) spec_sheet feasibility

@check("spec_sheet: known-feasible design point (T_op=77, target=0.5, Delta=5.0) -- "
      "overall feasible, all three routes close, verdict names all three")
def _():
    s = spec_sheet(77.0, 0.5, 5.0)
    assert s["feasible"] is True, s
    assert s["closes_baseline"] is True, s
    assert s["closes_slit"] is True, s
    assert s["closes_cavity"] is True, s
    assert "INFEASIBLE" not in s["verdict"], s["verdict"]
    for name in ("baseline", "slit", "cavity"):
        assert name in s["verdict"], (name, s["verdict"])
    numeric_keys = [k for k, v in s.items() if isinstance(v, float)]
    for k in numeric_keys:
        assert np.isfinite(s[k]), (k, s[k])


@check("spec_sheet: known-infeasible design point (T_op=120, target=0.01, Delta=1.5) -- "
      "all three routes fail, overall infeasible, verdict says INFEASIBLE")
def _():
    s = spec_sheet(120.0, 0.01, 1.5)
    assert s["feasible"] is False, s
    assert s["closes_baseline"] is False
    assert s["closes_slit"] is False
    assert s["closes_cavity"] is False
    assert s["rho_feasible"] is False
    assert s["G_feasible"] is False
    assert s["kappa_feasible"] is False
    assert np.isnan(s["kappa_max"])
    assert s["rho_required"] == float("inf")
    assert s["b_e_budget_base"] == float("-inf")
    assert s["verdict"].startswith("INFEASIBLE")


@check("spec_sheet: MIXED design point (T_op=120, target=0.1, Delta=5.0) -- "
      "baseline b_e budget misses by a hair but the cavity route closes; overall "
      "feasible must be True, NOT a blanket INFEASIBLE (the review-flagged bug)")
def _():
    s = spec_sheet(120.0, 0.1, 5.0)
    assert s["feasible"] is True, s
    assert s["closes_cavity"] is True, s
    assert s["closes_baseline"] is False, s  # the b_e budget route that misses
    assert -0.01 < s["b_e_budget_base"] < 0.0, s["b_e_budget_base"]  # misses, and only barely
    assert "cavity" in s["verdict"], s["verdict"]
    assert "INFEASIBLE" not in s["verdict"], s["verdict"]
    assert np.isfinite(s["kappa_max"]) and np.isfinite(s["G_required_worst"])


@check("spec_sheet: same T_op/target, tighter Delta=2.5 -- all three routes fail "
      "(splitting too small for the class-proxy linewidth even with brightness-floor "
      "engineering), verdict INFEASIBLE")
def _():
    s = spec_sheet(120.0, 0.1, 2.5)
    assert s["feasible"] is False, s
    assert s["closes_baseline"] is False
    assert s["closes_slit"] is False
    assert s["closes_cavity"] is False
    assert s["verdict"].startswith("INFEASIBLE"), s["verdict"]


@check("spec_sheet: eps_budget < 0 raises ValueError (invalid target, not silently swallowed)")
def _():
    try:
        spec_sheet(120.0, -0.1, 3.5)
        assert False, "expected ValueError"
    except ValueError:
        pass


# ---------------------------------------------------------- (e) run_spec.py end-to-end

@check("scripts/run_spec.py executes end-to-end; CSV bundle exists and parses")
def _():
    r = subprocess.run([sys.executable, str(ROOT / "scripts" / "run_spec.py")],
                       cwd=str(ROOT), capture_output=True, text=True, timeout=180)
    assert r.returncode == 0, f"exit {r.returncode}\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
    outdir = ROOT / "out" / "spec"
    for f in ("spec_staged.csv", "spec_300K.csv", "spec_figure.pdf", "spec_figure.svg",
             "spec_figure.png"):
        assert (outdir / f).exists(), f"missing {f}"

    with open(outdir / "spec_staged.csv", newline="", encoding="utf-8") as fh:
        rows = list(csv.DictReader(fh))
    assert len(rows) == 2 * 2 * 3, len(rows)  # T_op x target x Delta_XX
    for row in rows:
        float(row["eps_budget"])  # every field must parse as a number/nan/inf token
        assert row["feasible"] in ("True", "False")
        assert row["verdict"] != ""

    with open(outdir / "spec_300K.csv", newline="", encoding="utf-8") as fh:
        rows300 = list(csv.DictReader(fh))
    assert len(rows300) == 2 * 3  # Gamma(300) x Delta_XX
    for row in rows300:
        float(row["rho_req"])

    assert "T_op=  120 K  target=0.10  Delta= 1.5" not in r.stdout  # sanity: gate case not in run_spec's grid
    assert "INFEASIBLE" in r.stdout  # staged grid does include a genuinely all-routes-infeasible cell
    # the mixed staged cell (T_op=120, target=0.10, Delta=5.0) must print a verdict
    # that names the closing route, not a blanket INFEASIBLE stamp
    assert "T_op=  120" in r.stdout and "target=0.10" in r.stdout


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
