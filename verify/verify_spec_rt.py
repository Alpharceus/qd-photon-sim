"""Regression + behavior suite for fsim_core.spec.spec_sheet_from_device (the
device-resolved inverse-design entry point added for the room-temperature
edge tier, docs/rt_edge_contract.md).

Every check here anchors on a value computed independently inside this file
(a closed-form inversion re-derived from the contract/module docstrings, or a
literal pre-edit reference captured before fsim_core/spec.py was touched) --
never on trusting spec_sheet_from_device's own output as ground truth for
itself.

Standalone, same style as verify_cw_g2.py / gate_spec.py.
Run: python verify/verify_spec_rt.py   (exit code 0 iff all pass)
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import fsim_core.device as device_mod
import fsim_core.spec as spec_mod
from fsim_core.device import DeviceDesign
from fsim_core.loading import f1b_g2
from fsim_core.spec import spec_sheet, spec_sheet_from_device

CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


def nan_eq(a, b):
    """Exact float equality, NaN/inf-aware (NaN == NaN, +-inf == +-inf)."""
    if isinstance(a, float) and isinstance(b, float):
        if np.isnan(a) and np.isnan(b):
            return True
        return a == b
    return a == b


def assert_dict_exact(got, ref, ctx):
    assert set(got) == set(ref), f"{ctx}: key set differs -- {set(got) ^ set(ref)}"
    for k in ref:
        assert nan_eq(got[k], ref[k]), f"{ctx}: key {k!r} differs: {got[k]!r} != {ref[k]!r}"


def assert_raises_naming(fn, field_substr, *args, **kwargs):
    try:
        fn(*args, **kwargs)
    except ValueError as e:
        assert field_substr in str(e), f"message {e!r} does not name {field_substr!r}"
        return
    raise AssertionError(f"{fn.__name__} did not raise ValueError")


# ----------------------------------------------------------------- fixtures

# Literal reference values for spec_sheet() at two design points -- captured
# from fsim_core/spec.py BEFORE this task's edit (integration-a/-b, the
# rt-edge-inverse-accounting spec). (77.0, 0.5, 5.0) is the point
# verify/gate_spec.py itself exercises (its check "spec_sheet: known-feasible
# point closes all three routes").
REF_1 = dict(
    T_op=77.0, target_g2=0.5, delta_xx=5.0, mu=0.5,
    gamma_op=0.8705814822595512, gamma_xx_op=0.6274899214857147,
    S_op=0.8897368818807402, B0_op=0.012772755597174409,
    eps_budget=0.5424270310822167, t_x_floor=0.3,
    eps_op_gamma_window=0.00698060948542608, g2dot_base=0.008108461837463062,
    rho_required=0.709991107851761, rho_feasible=True,
    b_e_budget_base=0.35065659522722553, b_e_feasible=True, closes_baseline=True,
    w_floor=0.4435834210698239, eps_slit=0.0058949732683787354,
    g2dot_slit=0.006850824378822933, b_e_budget_slit=0.35145079977633337,
    slit_feasible=True, closes_slit=True,
    kappa_min=0.3731063495398077, eps_cav_min=0.012320937714864819,
    g2dot_cav=0.01427670027954407, kappa_max=9.997927150105113,
    kappa_feasible=True, G_required_best=0.038308011872943956,
    G_required_worst=0.8699568035178158, G_feasible=True,
    v_tilde_required=38.29011999176762, F_eff_target=3.0, closes_cavity=True,
    density_limit_cm2=598475959.3154613, penalty_budget=0.05, aperture_um=1.0,
    T_j_op=79.34760509403647, dT_J_op=2.347605094036467,
    mesa_min_um=0.88484962819109, mesa_feasible=True, dT_max=3.0,
    feasible=True,
    verdict="closes via: baseline; slit (b_e<=0.351); cavity (G>=0.0383..0.87, "
            "kappa 0.373..10 meV, V~<=38.3)",
    tag="[A]",
)

REF_2 = dict(
    T_op=120.0, target_g2=0.1, delta_xx=2.5, mu=0.5,
    gamma_op=2.3956715212521034, gamma_xx_op=1.7267306569333505,
    S_op=0.549265057085544, B0_op=0.012806802779688722,
    eps_budget=0.08936791728545013, t_x_floor=0.3,
    eps_op_gamma_window=0.22670439744947535, g2dot_base=0.23871795828036402,
    rho_required=float("inf"), rho_feasible=False,
    b_e_budget_base=float("-inf"), b_e_feasible=False, closes_baseline=False,
    w_floor=1.22065560870698, eps_slit=0.16743354910530342,
    g2dot_slit=0.1809505982343796, b_e_budget_slit=float("-inf"),
    slit_feasible=False, closes_slit=False,
    kappa_min=1.0267163662509016, eps_cav_min=0.28922462963372075,
    g2dot_cav=0.2964189330679933, kappa_max=float("nan"), kappa_feasible=False,
    G_required_best=float("inf"), G_required_worst=float("inf"), G_feasible=False,
    v_tilde_required=13.914540922082125, F_eff_target=3.0, closes_cavity=False,
    density_limit_cm2=217484777.49789664, penalty_budget=0.05, aperture_um=1.0,
    T_j_op=122.94208398341067, dT_J_op=2.9420839834106687,
    mesa_min_um=0.9902693910690893, mesa_feasible=True, dT_max=3.0,
    feasible=False,
    verdict="INFEASIBLE: eps floor exceeds budget at any brightness >= 0.3",
    tag="[A]",
)


@check("legacy regression: spec_sheet(77.0, 0.5, 5.0) bit-identical to pre-edit values "
       "(the point verify/gate_spec.py exercises)")
def _():
    assert_dict_exact(spec_sheet(77.0, 0.5, 5.0), REF_1, "spec_sheet(77.0, 0.5, 5.0)")


@check("legacy regression: spec_sheet(120.0, 0.1, 2.5) bit-identical to pre-edit values")
def _():
    assert_dict_exact(spec_sheet(120.0, 0.1, 2.5), REF_2, "spec_sheet(120.0, 0.1, 2.5)")


# ----------------------------------------------------- device-resolved fixture

def _device(T_hs=300.0):
    """Plain legacy-drive (non-transport) RT design: class-proxy linewidth/
    retention INSIDE device.evaluate() (its default), evaluated at a hot
    heatsink so the operating point is not the cryogenic legacy default."""
    d = DeviceDesign()
    d.thermal.T_hs = T_hs
    return d


def _confinement_device(T_hs=300.0):
    """Anchored linewidth + confinement retention: the one combination under
    which device.evaluate() itself never calls class_proxy_params()."""
    d = DeviceDesign()
    d.dot.linewidth = "anchored"
    d.ret.mode = "confinement"
    d.ret.preset = "InP/GaAs0.65P0.35/AlGaAs0.4 on GaAs"
    d.thermal.T_hs = T_hs
    return d


@check("spec_sheet_from_device: round trip -- rho_required recovered from f1b_g2 "
       "at the resolved mu and eps (closed-form check, independent of rho_required())")
def _():
    d = _device()
    target = 0.99
    s = spec_sheet_from_device(d, target)
    g2dot = float(f1b_g2(s["mu"], s["eps_op"]))
    assert abs(g2dot - s["g2dot"]) < 1e-12
    assert g2dot < target, "fixture must be in the feasible-by-cascade regime"
    rho_req_expected = np.sqrt((1.0 - target) / (1.0 - g2dot))
    assert abs(rho_req_expected - s["rho_required"]) < 1e-9
    g2_closed = 1.0 - rho_req_expected**2 * (1.0 - g2dot)
    assert abs(g2_closed - target) < 1e-9


@check("spec_sheet_from_device: round trip -- the raw background budget recovered "
       "through b_e_budget's own closing condition")
def _():
    d = _device()
    target = 0.99
    s = spec_sheet_from_device(d, target)
    S, rho_req = s["S_op"], s["rho_required"]
    B0_raw = S * (1.0 - rho_req) / rho_req - s["b_e_budget"]
    B_max = B0_raw + s["b_e_budget"]
    rho_at_Bmax = S / (S + B_max)
    assert abs(rho_at_Bmax - rho_req) < 1e-9
    assert s["feasible"] is True and s["b_e_budget"] > 0.0


@check("normalization: raw-unit and contract-unit (b_e_budget_norm) budgets differ "
       "by exactly the evaluator's collected-X scalar S_resolved")
def _():
    d = _device()
    s = spec_sheet_from_device(d, 0.99)
    assert s["b_e_units"]  # non-empty convention string
    assert abs(s["b_e_budget"] / s["b_e_budget_norm"] - s["S_op"]) < 1e-9 * abs(s["S_op"])
    assert abs(s["b_e_budget"] - s["b_e_budget_norm"] * s["S_op"]) < 1e-12


@check("Lemma 1: doubling the evaluator's collected-X scalar leaves eps, g2, "
       "rho_required, G_required and the contract-unit background budget unchanged, "
       "and scales only the raw-unit budget by two")
def _():
    d = _device()

    base_scalars = dict(
        T_j_op=300.0, gamma_op=10.0, eps_op=0.4, mu_resolved=0.5,
        rho_op=0.6, b_e_resolved=0.02, invalid_reasons=[],
    )

    def make_fake(S):
        scalars = dict(base_scalars, S_resolved=S)

        def fake_evaluate(design, T_grid=None):
            return {"curves": {}, "scalars": scalars}
        return fake_evaluate

    orig_evaluate = device_mod.evaluate
    try:
        device_mod.evaluate = make_fake(0.05)
        s1 = spec_sheet_from_device(d, 0.9)
        device_mod.evaluate = make_fake(0.10)
        s2 = spec_sheet_from_device(d, 0.9)
    finally:
        device_mod.evaluate = orig_evaluate

    for key in ("eps_op", "g2dot", "rho_required", "G_required", "b_e_budget_norm"):
        assert abs(s1[key] - s2[key]) < 1e-12, (key, s1[key], s2[key])
    assert s1["b_e_budget"] > 0.0
    assert abs(s2["b_e_budget"] - 2.0 * s1["b_e_budget"]) < 1e-9 * abs(s1["b_e_budget"])


@check("rejection: metric other than 'g2_pulsed' raises ValueError naming 'metric', "
       "for g2_cw0, g2_cw0_raw and an unknown metric string")
def _():
    d = _device()
    for bad in ("g2_cw0", "g2_cw0_raw", "not-a-metric"):
        assert_raises_naming(spec_sheet_from_device, "metric", d, 0.5, metric=bad)


@check("rejection: a design requesting CW operation raises ValueError naming the "
       "field ('drive.cw'/'cw'), never approximated with the pulsed formulas")
def _():
    d = _device()
    d.drive.cw = True
    assert_raises_naming(spec_sheet_from_device, "cw", d, 0.5)


@check("rejection: an evaluator-invalid design (thermal runaway) raises ValueError "
       "surfacing the evaluator's own invalid reason")
def _():
    d = _device()
    d.drive.I_uA = 1.0e5
    d.drive.V = 50.0
    assert_raises_naming(spec_sheet_from_device, "thermal runaway", d, 0.5)


@check("no-approximation guard: every rejection path above raised (none returned "
       "a finite dict) -- rerun each and assert the call never completes")
def _():
    d_cw = _device()
    d_cw.drive.cw = True
    d_hot = _device()
    d_hot.drive.I_uA, d_hot.drive.V = 1.0e5, 50.0
    for fn, args, kwargs in (
        (spec_sheet_from_device, (_device(), 0.5), dict(metric="g2_cw0")),
        (spec_sheet_from_device, (_device(), 0.5), dict(metric="g2_cw0_raw")),
        (spec_sheet_from_device, (d_cw, 0.5), {}),
        (spec_sheet_from_device, (d_hot, 0.5), {}),
    ):
        try:
            result = fn(*args, **kwargs)
        except ValueError:
            continue
        raise AssertionError(f"expected ValueError, got a return value: {result!r}")


@check("the adapter never consults class_proxy_params() on the resolved path: "
       "spec_sheet_from_device still succeeds and spec_sheet still fails when it "
       "is monkeypatched to raise")
def _():
    d = _confinement_device()  # linewidth='anchored' + ret.mode='confinement':
                               # the one combination under which device.evaluate()
                               # itself never touches class_proxy_params() either,
                               # isolating whether spec_sheet_from_device's OWN
                               # body calls it.
    orig = device_mod.class_proxy_params

    def boom():
        raise RuntimeError("class_proxy_params called on the resolved path")

    device_mod.class_proxy_params = boom
    try:
        s = spec_sheet_from_device(d, 0.99)
        assert s["source"] == "device-resolved"
        try:
            spec_sheet(77.0, 0.5, 5.0)
        except RuntimeError:
            pass
        else:
            raise AssertionError("spec_sheet should have failed with class_proxy_params patched")
    finally:
        device_mod.class_proxy_params = orig


def main():
    failed = 0
    for name, fn in CHECKS:
        try:
            fn()
            print(f"  PASS  {name}")
        except AssertionError as e:
            failed += 1
            print(f"* FAIL  {name}  {e}")
        except Exception as e:  # noqa: BLE001 -- surface unexpected errors as failures too
            failed += 1
            print(f"* FAIL  {name}  unexpected {type(e).__name__}: {e}")
    n = len(CHECKS)
    print(f"\n{n - failed}/{n} spec RT checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
