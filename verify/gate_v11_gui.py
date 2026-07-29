"""v1.1 GUI gate: designer drive-node controls, injection-engineering presets,
the Streamlit EL/PL toggle, and the single-implementation rule for the F5'
injection-broadening formula. Standalone, same style as gate_d1.py /
gate_d2.py / gate_d3.py.
Run: python verify/gate_v11_gui.py   (exit code 0 iff all pass)
"""
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.device import DeviceDesign, evaluate  # noqa: E402
from fsim_core.presets import INJECTION_PRESETS, apply_injection_preset  # noqa: E402

CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


# ------------------------------------------------------------- (a)/(b)/(c) designer subprocess

@check("designer smoke test: python fsim_gui/designer.py --frames 5 exits 0")
def _():
    r = subprocess.run([sys.executable, str(ROOT / "fsim_gui" / "designer.py"),
                        "--frames", "5"], cwd=str(ROOT),
                       capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, f"exit {r.returncode}\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"


@check("designer selftest: python fsim_gui/designer.py --selftest <dir> exits 0")
def _():
    with tempfile.TemporaryDirectory() as td:
        outdir = Path(td) / "bundle"
        r = subprocess.run([sys.executable, str(ROOT / "fsim_gui" / "designer.py"),
                            "--selftest", str(outdir)], cwd=str(ROOT),
                           capture_output=True, text=True, timeout=120)
        assert r.returncode == 0, f"exit {r.returncode}\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"


@check("designer v1.1 round-trip THROUGH WIDGETS: --roundtrip-check exits 0")
def _():
    r = subprocess.run([sys.executable, str(ROOT / "fsim_gui" / "designer.py"),
                        "--roundtrip-check"], cwd=str(ROOT),
                       capture_output=True, text=True, timeout=120)
    assert r.returncode == 0, f"exit {r.returncode}\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
    assert "roundtrip-check: OK" in r.stdout, r.stdout


# --------------------------------------------------------------- (d) injection presets

@check("INJECTION_PRESETS: apply each to a fresh design, evaluate() runs clean")
def _():
    for key in INJECTION_PRESETS:
        d = DeviceDesign()
        apply_injection_preset(d, key)
        res = evaluate(d)
        assert res["scalars"]["g2_op"] == res["scalars"]["g2_op"], key  # not NaN


@check("rti-quiet: F_p=0.5, and lower g2_op than standard at equal eps (mu=0.6, PL mode)")
def _():
    # PL mode disables dg_inj (F5' broadening never fires), so eps is held
    # identical between the two designs -- an apples-to-apples g2_op
    # comparison that isolates the F8/F8b quiet-pump effect. mu=0.6 keeps
    # both inside the F8 domain (floor = 1 - f8b_thin_fano(0.8, 0.5) = 0.4).
    d_std = DeviceDesign()
    d_std.drive.mode = "PL"
    d_std.drive.mu = 0.6
    apply_injection_preset(d_std, "standard")
    s_std = evaluate(d_std)["scalars"]

    d_rti = DeviceDesign()
    d_rti.drive.mode = "PL"
    d_rti.drive.mu = 0.6
    apply_injection_preset(d_rti, "rti-quiet")
    s_rti = evaluate(d_rti)["scalars"]

    assert d_rti.drive.F_p == 0.5, d_rti.drive.F_p
    assert abs(s_std["eps_op"] - s_rti["eps_op"]) < 1e-12, (s_std["eps_op"], s_rti["eps_op"])
    assert s_rti["g2_op"] < s_std["g2_op"], (s_rti["g2_op"], s_std["g2_op"])


# ------------------------------------------------------------- (e) Streamlit AppTest

@check("Streamlit app.py boots with zero exceptions; PL/EL radio exists")
def _():
    from streamlit.testing.v1 import AppTest
    at = AppTest.from_file(str(ROOT / "fsim_gui" / "app.py"), default_timeout=180)
    at.run()
    assert not at.exception, [str(e) for e in at.exception]
    radios = [r for r in at.radio if set(r.options) == {"PL", "EL"}]
    assert radios, "no PL/EL radio found"


# ------------------------------------------------------- (f) one-implementation check

@check("gamma_eff is the ONE implementation: not re-derived inline in fsim_gui/*.py")
def _():
    app_src = (ROOT / "fsim_gui" / "app.py").read_text(encoding="utf-8")
    designer_src = (ROOT / "fsim_gui" / "designer.py").read_text(encoding="utf-8")

    assert "gamma_eff" in app_src, "app.py must import/call fsim_core.loading.gamma_eff"
    assert re.search(r"from fsim_core\.loading import[^\n]*gamma_eff", app_src), \
        "app.py must import gamma_eff from fsim_core.loading"

    # the formula itself must not appear inline in any GUI file -- only
    # fsim_core.loading.gamma_eff may compute it. Opus v1.1 review found the
    # original regex under-strength (missed attribute/dict-accessed forms like
    # `d.drive.dg_inj * (...) ** d.drive.p_inj` and `p['dg_inj'] * ...`), so
    # match ANY dg_inj access multiplied into ANY exponentiation, and ANY
    # exponentiation by a p_inj access.
    inline_patterns = [
        re.compile(r"dg_inj[\"'\]\)\s]*\s*\*[^\n]*\*\*"),   # dg_inj * (...)**...
        re.compile(r"\*\*\s*[^\n]*p_inj"),                   # ... ** ...p_inj
    ]
    for name, src in (("app.py", app_src), ("designer.py", designer_src)):
        for pat in inline_patterns:
            assert not pat.search(src), \
                f"inline gamma_eff formula found in {name} (pattern {pat.pattern!r})"


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
