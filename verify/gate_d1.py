"""Phase D1 gate: presets + design_meta + designer must be consistent and the
GUI must still boot headless. Standalone, same style as verify_fsim.py.
Run: python verify/gate_d1.py   (exit code 0 iff all pass)
"""
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core.design_meta import missing_meta
from fsim_core.device import DeviceDesign, evaluate
from fsim_core.presets import preset_device

CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


COMBOS = [
    ("chatzarakis-class", "GaAs", "planar-lambda", "cw-electrical"),
    ("staged-inp-gaasp", "GaAs-Si", "none", "cw-electrical"),
    ("piezo-variant", "GaAs", "sin-waveguide", "pulsed-electrical"),
]


def _scalars_match(a, b):
    assert set(a) == set(b)
    for k in a:
        va, vb = a[k], b[k]
        if isinstance(va, float) and va != va:  # NaN
            assert isinstance(vb, float) and vb != vb, (k, va, vb)
            continue
        if isinstance(va, float):
            assert abs(va - vb) < 1e-12, (k, va, vb)
        else:
            assert va == vb, (k, va, vb)


@check("sin-waveguide preset: physically sane operating point (no phantom cavity)")
def _():
    # Opus D1 finding: sin mode must not inherit mode-tracking detuning -- a
    # waveguide has no mode to track. dx must follow filter.dx, eps stays < 1.
    d = preset_device("piezo-variant", "GaAs", "sin-waveguide", "pulsed-electrical")
    s = evaluate(d)["scalars"]
    assert s["eps_op"] < 1.0 and s["g2_op"] <= 1.0, (s["eps_op"], s["g2_op"])
    d2 = preset_device("piezo-variant", "GaAs", "sin-waveguide", "pulsed-electrical")
    d2.filter.dx = 2.0  # user's manual window offset must be honored in sin mode
    assert abs(evaluate(d2)["scalars"]["eps_op"] - s["eps_op"]) > 1e-9


@check("presets: round-trip through YAML preserves every scalar (3 combos)")
def _():
    with tempfile.TemporaryDirectory() as td:
        for dot, template, cavity, drive in COMBOS:
            d = preset_device(dot, template, cavity, drive, name="gate-d1")
            s0 = evaluate(d)["scalars"]
            p = Path(td) / "d.yaml"
            d.save(p)
            e = DeviceDesign.load(p)
            s1 = evaluate(e)["scalars"]
            _scalars_match(s0, s1)


@check("Lemma 1: beta_sin (0.3 vs 3.0) leaves g2/eps/rho2/T_c untouched, brightness x10")
def _():
    d1, d2 = DeviceDesign(), DeviceDesign()
    for d in (d1, d2):
        d.cavity.enabled = True
        d.cavity.type = "sin_waveguide"
    d1.cavity.beta_sin, d2.cavity.beta_sin = 0.3, 3.0
    r1, r2 = evaluate(d1), evaluate(d2)
    c1, c2 = r1["curves"], r2["curves"]
    for k in ("g2", "eps", "rho2"):
        assert np.allclose(c1[k], c2[k], equal_nan=True), k
    s1, s2 = r1["scalars"], r2["scalars"]
    assert (s1["T_c"] == s2["T_c"]) or (s1["T_c"] != s1["T_c"] and s2["T_c"] != s2["T_c"])
    assert abs(s2["brightness_per_pulse"] / s1["brightness_per_pulse"] - 10.0) < 1e-9


@check("META completeness: missing_meta() is empty")
def _():
    m = missing_meta()
    assert m == [], m


@check("designer smoke test: python fsim_gui/designer.py --frames 5 exits 0")
def _():
    root = Path(__file__).resolve().parents[1]
    r = subprocess.run([sys.executable, str(root / "fsim_gui" / "designer.py"),
                        "--frames", "5"], cwd=str(root),
                       capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, f"exit {r.returncode}\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"


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
