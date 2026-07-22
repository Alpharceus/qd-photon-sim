"""Phase D2 gate: envelope mode (fsim_core.device.evaluate_envelope) must
reproduce the existing staged-device envelope (scripts/run_phase3.run_staged),
degenerate correctly to evaluate(), respect the sampling cap, and produce a
sane/deterministic tornado -- with the designer GUI still booting headless in
its (now default-on) envelope path. Standalone, same style as gate_d1.py.
Run: python verify/gate_d2.py   (exit code 0 iff all pass)
"""
import copy
import json
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "scripts"))

from fsim_core.design_meta import default_ranged
from fsim_core.device import DeviceDesign, evaluate, evaluate_envelope

CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


def _staged_design() -> DeviceDesign:
    """The DeviceDesign equivalent of run_phase3.run_staged()'s hardcoded
    setup (mesa 1 um, 1.5 um/k=8.0/alpha=0.5 pillar epi on GaAs 55/1.25,
    1.9 V / 10 uA CW, auto-w filter, cavity off, class-proxy retention).

    Two overrides are required for a bit-for-bit match, beyond what a literal
    reading of the staged card gives:
      * dot.r_xx -- run_staged uses the *exact* V-a joint-fit r_xx
        (fit_params.json), not DotBlock's rounded default (0.72). r_xx isn't
        swept in this comparison, so it must match exactly.
      * drive.b_e_Eact = 0.0 -- run_staged's b_e_total is a flat additive
        background (rho = S/(S + b0 + beta(1-S) + b_e)). DeviceDesign's
        default (100 meV) would run that background through the injection
        channel's Arrhenius activation instead, which run_staged's simpler
        placeholder does not model. I_uA == I_ref_uA (both 10) keeps the
        channel's current-scaling factor at 1 either way.
    """
    fitp = json.loads((ROOT / "out" / "phase0" / "fit_params.json").read_text())["params"]
    d = DeviceDesign(name="gate-d2-staged")
    d.dot.r_xx = fitp["r_xx"]
    d.drive.V = 1.9
    d.drive.I_uA = 10.0
    d.drive.duty = 1.0
    d.drive.b_e_Eact = 0.0
    d.thermal.mesa_diameter_um = 1.0
    d.thermal.layers = [{"name": "epi/cladding", "t_um": 1.5, "k300": 8.0,
                         "alpha": 0.5, "spread": False}]
    d.thermal.substrate = {"name": "GaAs", "k300": 55.0, "alpha": 1.25}
    d.filter.auto_w = True
    d.cavity.enabled = False
    # retention block left at its all-zero default -> evaluate() falls back
    # to the class proxy, exactly what run_staged uses.
    return d


@check("band-for-band vs run_phase3.run_staged (77 K & 120 K, 15 deltas, <1e-9)")
def _():
    from run_phase3 import run_staged

    ref = run_staged()
    deltas = ref["deltas"]
    base = _staged_design()
    ranged = {"dot.gamma_scale": (0.7, 1.0, 1.3),
             "drive.b_e": (1e-3, 0.03, 0.3),
             "drive.mu": (0.1, 0.5, 1.0)}

    worst = 0.0
    for T_hs in (77.0, 120.0):
        band_ref = ref["envelopes"][T_hs]["band"]
        for i, delta in enumerate(deltas):
            d = copy.deepcopy(base)
            d.dot.delta_xx = float(delta)
            env = evaluate_envelope(d, ranged, T_grid=[T_hs])
            lo, hi = env["bands"]["g2"]
            lo_ref, hi_ref = band_ref[i]
            worst = max(worst, abs(float(lo[0]) - lo_ref), abs(float(hi[0]) - hi_ref))
            assert abs(float(lo[0]) - lo_ref) < 1e-9, (T_hs, delta, lo[0], lo_ref)
            assert abs(float(hi[0]) - hi_ref) < 1e-9, (T_hs, delta, hi[0], hi_ref)


@check("envelope degeneracy: ranged={} matches evaluate() exactly (curves + scalars)")
def _():
    for d in (DeviceDesign(), _staged_design()):
        ref = evaluate(d, T_grid=[77.0, 120.0, 200.0])
        env = evaluate_envelope(d, {}, T_grid=[77.0, 120.0, 200.0])
        assert env["n_samples"] == 1
        for k, v in ref["curves"].items():
            assert np.array_equal(v, env["mid"][k], equal_nan=True), k
        for k, v in ref["scalars"].items():
            lo, hi = env["scalar_bands"].get(k, (None, None))
            if lo is None:
                continue  # scalar_bands only covers T_c/g2_op/eps_op/rho_op/dT_J
            if v != v:
                assert lo != lo and hi != hi, (k, v, lo, hi)
            else:
                assert lo == v == hi, (k, v, lo, hi)


@check("extremes cap: >4096 cartesian-product combinations raises ValueError")
def _():
    d = DeviceDesign()
    big = {
        "dot.delta_xx": list(np.linspace(1.5, 5.0, 13)),
        "dot.gamma_scale": list(np.linspace(0.7, 1.3, 13)),
        "drive.b_e": list(np.linspace(1e-3, 0.3, 13)),
        "drive.mu": list(np.linspace(0.1, 1.0, 13)),
        "cavity.kappa": list(np.linspace(0.3, 3.0, 3)),
    }
    n = 1
    for v in big.values():
        n *= len(v)
    assert n > 4096
    try:
        evaluate_envelope(d, big, T_grid=[77.0])
        assert False, "expected ValueError for an over-cap sweep"
    except ValueError:
        pass


@check("tornado sanity on the staged design: widths >= 0, ordering deterministic")
def _():
    d = _staged_design()
    ranged = default_ranged(d)
    assert set(ranged) == {"dot.delta_xx", "dot.gamma_scale", "drive.b_e", "drive.mu"}, ranged

    env1 = evaluate_envelope(d, ranged)
    env2 = evaluate_envelope(d, ranged)
    assert env1["tornado"] == env2["tornado"]
    for path, red in env1["tornado"].items():
        assert red == red, f"{path}: width_reduction is NaN"
        assert red >= -1e-9, (path, red)
    order1 = sorted(env1["tornado"], key=lambda p: -env1["tornado"][p])
    order2 = sorted(env2["tornado"], key=lambda p: -env2["tornado"][p])
    assert order1 == order2


@check("designer smoke test: python fsim_gui/designer.py --frames 5 exits 0 (envelope path)")
def _():
    r = subprocess.run([sys.executable, str(ROOT / "fsim_gui" / "designer.py"),
                        "--frames", "5"], cwd=str(ROOT),
                       capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, f"exit {r.returncode}\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
    assert " - " in r.stdout, f"expected an interval in smoke output:\n{r.stdout}"


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
