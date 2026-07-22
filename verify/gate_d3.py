"""Phase D3 gate: design comparison + one-click reporting. Checks the
fsim_viz.report core path (bundle completeness, re-readable CSVs, exact
numeric round-trip), the three-layer rule on the new module, the real GUI
selftest end-to-end, and bundle determinism. Standalone, same style as
gate_d1.py / gate_d2.py.
Run: python verify/gate_d3.py   (exit code 0 iff all pass)
"""
import csv
import filecmp
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.design_meta import default_ranged
from fsim_core.device import evaluate, evaluate_envelope
from fsim_core.presets import preset_device
from fsim_viz.report import designer_report

CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


def _entries():
    """Two entries built directly against fsim_core (no GUI): A = the staged
    preset device under its default envelope, B = a delta_xx=5.0 point
    variant -- covers both the banded and point-only shapes designer_report
    must accept."""
    dA = preset_device("staged-inp-gaasp", "GaAs", "none", "cw-electrical", name="gate-d3-A")
    ranged = default_ranged(dA)
    env = evaluate_envelope(dA, ranged)
    entry_a = {"label": "A", "name": dA.name, "design": dA, "curves": env["mid"],
              "bands": env["bands"], "scalars": evaluate(dA)["scalars"],
              "scalar_bands": env["scalar_bands"], "ranged": ranged}

    dB = preset_device("staged-inp-gaasp", "GaAs", "none", "cw-electrical", name="gate-d3-B")
    dB.dot.delta_xx = 5.0
    res_b = evaluate(dB)
    entry_b = {"label": "B", "name": dB.name, "design": dB, "curves": res_b["curves"],
              "scalars": res_b["scalars"]}
    return [entry_a, entry_b], env


@check("designer_report: full bundle written, curves_A.csv exact to 1e-12, scalars parse")
def _():
    entries, env = _entries()
    with tempfile.TemporaryDirectory() as td:
        outdir = designer_report(entries, Path(td), title="gate-d3")

        expected = ["report.pdf", "report.svg", "report.png",
                   "curves_A.csv", "curves_B.csv", "scalars_comparison.csv",
                   "design_A.yaml", "design_B.yaml", "ranged_A.csv"]
        for f in expected:
            assert (outdir / f).exists(), f"missing {f}"
        assert not (outdir / "ranged_B.csv").exists(), "B has no ranged dict -- no file expected"

        with open(outdir / "curves_A.csv", newline="", encoding="utf-8") as fh:
            rows = list(csv.DictReader(fh))
        mid = env["mid"]
        g2_lo, g2_hi = env["bands"]["g2"]
        assert len(rows) == len(mid["T_hs"])
        cols = {"T_hs": mid["T_hs"], "Tj": mid["Tj"], "g2": mid["g2"], "eps": mid["eps"],
               "rho2": mid["rho2"], "gamma": mid["gamma"], "g2_lo": g2_lo, "g2_hi": g2_hi}
        for name, arr in cols.items():
            for row, ref in zip(rows, arr):
                assert abs(float(row[name]) - float(ref)) < 1e-12, (name, row[name], ref)

        with open(outdir / "scalars_comparison.csv", newline="", encoding="utf-8") as fh:
            srows = list(csv.DictReader(fh))
        assert set(srows[0].keys()) == {"scalar", "A", "B"}, srows[0].keys()
        assert len(srows) == 10, len(srows)  # SCALAR_ROWS
        for row in srows:
            for label in ("A", "B"):
                v = row[label]
                if ".." in v:
                    lo, hi = v.split("..")
                    float(lo); float(hi)
                else:
                    float(v) if v != "nan" else None


@check("three-layer rule: fsim_viz/report.py never references streamlit/plotly/dearpygui")
def _():
    src = (ROOT / "fsim_viz" / "report.py").read_text(encoding="utf-8")
    for banned in ("streamlit", "plotly", "dearpygui"):
        assert banned not in src, banned


@check("GUI selftest: python fsim_gui/designer.py --selftest <dir> exits 0, bundle complete")
def _():
    with tempfile.TemporaryDirectory() as td:
        outdir = Path(td) / "bundle"
        r = subprocess.run([sys.executable, str(ROOT / "fsim_gui" / "designer.py"),
                            "--selftest", str(outdir)], cwd=str(ROOT),
                           capture_output=True, text=True, timeout=120)
        assert r.returncode == 0, f"exit {r.returncode}\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
        for f in ("report.png", "report.pdf", "curves_A.csv", "curves_B.csv",
                 "scalars_comparison.csv", "design_A.yaml", "design_B.yaml"):
            assert (outdir / f).exists(), f"selftest bundle missing {f}"


@check("determinism: two designer_report() calls on the same entries -> byte-identical CSVs")
def _():
    entries, _env = _entries()
    with tempfile.TemporaryDirectory() as td:
        out1 = designer_report(entries, Path(td) / "run1", title="gate-d3")
        out2 = designer_report(entries, Path(td) / "run2", title="gate-d3")
        for name in ("curves_A.csv", "curves_B.csv", "scalars_comparison.csv",
                    "ranged_A.csv", "design_A.yaml", "design_B.yaml"):
            assert filecmp.cmp(out1 / name, out2 / name, shallow=False), name


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
