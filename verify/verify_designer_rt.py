"""RT-edge GUI preservation checker (spec rt-edge-gui-preservation).

Standalone: drives fsim_gui/designer.py ONLY through subprocess, and never
imports dearpygui itself (this file only needs fsim_core.device for building
the probe design and comparing YAML field-by-field). Exits 0 iff every check
passes and prints "N/N designer RT checks passed".

Run: python verify/verify_designer_rt.py
"""
import dataclasses
import re
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.device import DeviceDesign  # noqa: E402

PY = sys.executable
DESIGNER = str(ROOT / "fsim_gui" / "designer.py")

CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


def _run(args, timeout=180):
    return subprocess.run([PY, DESIGNER] + args, cwd=str(ROOT),
                          capture_output=True, text=True, timeout=timeout)


def _run_script(path, args=(), timeout=300):
    return subprocess.run([PY, str(path), *args], cwd=str(ROOT),
                          capture_output=True, text=True, timeout=timeout)


# --------------------------------------------------------------- (1) round-trip

# Float paths that pass through a Dear PyGui add_input_float widget (float32
# storage, ~7 significant figures -- see fsim_gui/designer.py WIDGET_TAG):
# these get a relative tolerance, everything else must be exact (nothing else
# is ever touched by a widget, so it round-trips through the preservation
# baseline byte/value-exact).
WIDGET_FLOAT_PATHS = {
    "dot.delta_xx", "dot.gamma_scale", "dot.r_xx",
    "drive.V", "drive.I_uA", "drive.duty", "drive.mu", "drive.b_e",
    "drive.b_e_m", "drive.b_e_Eact", "drive.dg_inj", "drive.p_inj",
    "drive.F_p", "drive.eta_capture", "drive.C_dep_pF",
    "thermal.mesa_diameter_um", "thermal.T_hs",
    "cavity.kappa", "cavity.T_track", "cavity.E_X0", "cavity.F_P",
    "cavity.G", "cavity.beta_sin",
    "filter.w", "filter.dx",
    "aperture.density_cm2", "aperture.diameter_um", "aperture.sigma_inh",
    "aperture.comp_brightness",
}
TOL = 1e-5

# Legal non-default values for the enumerated-string fields (per each
# field's own docstring in fsim_core/device.py -- drive.mode's docstring
# lists only "EL"/"PL", so "EL-transport" -- a third value evaluate() also
# accepts -- is intentionally not used here).
ENUM_CHOICES = {
    "dot.lineshape": "ibm",               # default "lorentzian"
    "dot.linewidth": "anchored",          # default "class"
    "ret.mode": "confinement",            # default "proxy"
    "ret.channel": "electron",            # default "pair_half"
    "drive.mode": "PL",                   # default "EL"
    "drive.mechanism": "poisson-rail",    # default ""
    "drive.cw_irf_shape": "exponential",  # default "gaussian"
    "cavity.type": "sin_waveguide",       # default "planar"
    "filter.track": "hold",               # default "mode"
    "filter.track_material": "dot",       # default "" (choices "" | "dot" | "matrix")
    "emission.type": "edge",              # default "none"
    # pkg4 fix's DeviceDesign.load()/evaluate() now reject any
    # drive.loading_model outside ("auto", "capped_poisson",
    # "moment_matched") -- without this entry _non_default's generic
    # string bump turned the "auto" default into "auto-probe", which
    # --collect-dump's own DeviceDesign.load() then rejected (not a
    # pr-pkg6-stale-text C5 item; a pre-existing gap this fixture's
    # ENUM_CHOICES needs regardless, or this check cannot run at all).
    "drive.loading_model": "capped_poisson",  # default "auto"
}


def _bump_float(v: float) -> float:
    return v + max(abs(v), 1.0) * 0.2537 + 0.0091


def _non_default(path: str, value):
    if isinstance(value, bool):
        return not value
    if isinstance(value, float):
        return _bump_float(value)
    if isinstance(value, str):
        if path in ENUM_CHOICES:
            return ENUM_CHOICES[path]
        return (value + "-probe") if value else "probe-value"
    if isinstance(value, dict):
        return {"probe_flag": True, "probe_num": 1.5}
    if value is None:  # EmissionBlock.R_back: float | None
        return 0.42
    raise TypeError(f"{path}: don't know how to bump a default of type {type(value)}")


def build_non_default_design() -> DeviceDesign:
    """Every scalar/dict field of every DeviceDesign block, set to something
    other than its dataclass default -- walks dataclasses.fields recursively
    (acceptance criterion 2), no hand-written per-field list. thermal.layers/
    substrate are handled separately (a list of dicts / one dict edited
    through their own stack-table widget, not a generic scalar)."""
    d = DeviceDesign()
    d.name = "rt-edge-preservation-check"
    for block_field in dataclasses.fields(d):
        block_name = block_field.name
        if block_name in ("name", "provenance"):
            continue
        block = getattr(d, block_name)
        if not dataclasses.is_dataclass(block):
            continue
        for f in dataclasses.fields(block):
            if block_name == "thermal" and f.name in ("layers", "substrate"):
                continue
            path = f"{block_name}.{f.name}"
            value = getattr(block, f.name)
            setattr(block, f.name, _non_default(path, value))
    d.provenance = {"probe_note": "rt-edge preservation check"}
    d.thermal.layers = [
        {"name": "probe-cap", "t_um": 0.73, "k300": 21.5, "alpha": 0.9, "spread": True},
        {"name": "probe-spacer", "t_um": 1.21, "k300": 44.0, "alpha": 1.4, "spread": False},
    ]
    d.thermal.substrate = {"name": "InP", "k300": 68.0, "alpha": 0.9}
    return d


def _cmp_value(path: str, expected, actual, tol_ok: bool, errors: list, tmp_out: Path):
    if isinstance(expected, float) or isinstance(actual, float):
        if expected is None or actual is None:
            if expected != actual:
                errors.append(f"{tmp_out}: field {path} changed: expected {expected!r} got {actual!r}")
            return
        if tol_ok:
            denom = abs(expected) if expected else 1.0
            ok = abs(actual - expected) / denom < TOL
        else:
            ok = expected == actual
        if not ok:
            errors.append(f"{tmp_out}: field {path} reverted to default or mismatched: "
                          f"expected {expected!r} got {actual!r}")
    else:
        if expected != actual:
            errors.append(f"{tmp_out}: field {path} reverted to default or mismatched: "
                          f"expected {expected!r} got {actual!r}")


def _cmp_layers(expected: list, actual: list, errors: list, tmp_out: Path):
    if len(expected) != len(actual):
        errors.append(f"{tmp_out}: thermal.layers length changed: "
                      f"expected {len(expected)} got {len(actual)}")
        return
    for i, (exp, act) in enumerate(zip(expected, actual)):
        if exp.get("name") != act.get("name"):
            errors.append(f"{tmp_out}: thermal.layers[{i}].name mismatch: "
                          f"{exp.get('name')!r} vs {act.get('name')!r}")
        if bool(exp.get("spread", False)) != bool(act.get("spread", False)):
            errors.append(f"{tmp_out}: thermal.layers[{i}].spread mismatch")
        for key in ("t_um", "k300"):
            e, a = float(exp[key]), float(act[key])
            if abs(a - e) / (abs(e) if e else 1.0) > TOL:
                errors.append(f"{tmp_out}: thermal.layers[{i}].{key} mismatch: {e} vs {a}")
        # alpha carries no widget at all -- must be exact
        if exp.get("alpha") != act.get("alpha"):
            errors.append(f"{tmp_out}: thermal.layers[{i}].alpha changed though it has "
                          f"no widget: expected {exp.get('alpha')!r} got {act.get('alpha')!r}")


def _compare_designs(expected: DeviceDesign, actual: DeviceDesign, tmp_out: Path) -> list:
    errors = []
    if expected.name != actual.name:
        errors.append(f"{tmp_out}: field name changed: expected {expected.name!r} "
                      f"got {actual.name!r}")
    for block_field in dataclasses.fields(expected):
        block_name = block_field.name
        if block_name == "name":
            continue
        if block_name == "provenance":
            _cmp_value("provenance", expected.provenance, actual.provenance, False, errors, tmp_out)
            continue
        exp_block, act_block = getattr(expected, block_name), getattr(actual, block_name)
        for f in dataclasses.fields(exp_block):
            if block_name == "thermal" and f.name in ("layers", "substrate"):
                continue
            path = f"{block_name}.{f.name}"
            _cmp_value(path, getattr(exp_block, f.name), getattr(act_block, f.name),
                      path in WIDGET_FLOAT_PATHS, errors, tmp_out)
    _cmp_layers(expected.thermal.layers, actual.thermal.layers, errors, tmp_out)
    if expected.thermal.substrate != actual.thermal.substrate:
        errors.append(f"{tmp_out}: thermal.substrate changed though it has no widget: "
                      f"expected {expected.thermal.substrate!r} got {actual.thermal.substrate!r}")
    return errors


@check("collect_design() round-trip preserves every DeviceDesign field "
       "(widget-backed or not) through --design/--collect-dump")
def _():
    with tempfile.TemporaryDirectory() as td:
        tmp_in = Path(td) / "in-design.yaml"
        tmp_out = Path(td) / "out-design.yaml"
        expected = build_non_default_design()
        expected.save(tmp_in)
        r = _run(["--design", str(tmp_in), "--collect-dump", str(tmp_out)])
        assert r.returncode == 0, (f"--collect-dump exited {r.returncode}\n"
                                   f"STDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}")
        assert tmp_out.exists(), f"{tmp_out} was not written"
        actual = DeviceDesign.load(tmp_out)
        errors = _compare_designs(expected, actual, tmp_out)
        assert not errors, "\n".join(errors)


# ------------------------------------------------------------------ (2) screenshot

DERIVED_LABELS = [
    "junction temperature (transport):",
    "junction voltage V_j:",
    "applied voltage V_applied:",
    "built-in voltage V_bi:",
    "junction power P_junction:",
    "injection efficiency eta_inj:",
    "capture efficiency eta_capture:",
    "loading mu (resolved):",
    "background per collected X photon:",
    "g2_cw0 (intrinsic CW):",
    "g2_cw0_raw (IRF-convolved CW):",
]
_VALUE_RE = r"\s*(n/a|[+-]?\d)"


@check("--frames 10 --screenshot out/rt_edge/gui-smoke.png writes a real PNG "
       "+ a .txt with every derived-transport/CW label")
def _():
    out_png = ROOT / "out" / "rt_edge" / "gui-smoke.png"
    r = _run(["--frames", "10", "--screenshot", str(out_png)])
    assert r.returncode == 0, (f"--screenshot exited {r.returncode}\n"
                               f"STDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}")
    assert out_png.exists(), f"{out_png} was not written"
    data = out_png.read_bytes()
    assert data[:8] == b"\x89PNG\r\n\x1a\n", f"{out_png} does not start with the PNG magic bytes"
    assert len(data) > 2048, f"{out_png} is only {len(data)} bytes (expected > 2 kB)"
    txt_path = Path(str(out_png) + ".txt")
    assert txt_path.exists(), f"{txt_path} was not written"
    content = txt_path.read_text(encoding="utf-8")
    missing = [lbl for lbl in DERIVED_LABELS if not re.search(re.escape(lbl) + _VALUE_RE, content)]
    assert not missing, f"{txt_path} missing derived label(s) (or non-numeric/n/a value): {missing}"


# ------------------------------------------------------------- (3)-(4) pre-existing CLI

@check("pre-existing CLI unaffected: --roundtrip-check still exits 0")
def _():
    r = _run(["--roundtrip-check"])
    assert r.returncode == 0, f"exit {r.returncode}\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
    assert "roundtrip-check: OK" in r.stdout, r.stdout


@check("pr-pkg1-fix4 item 3: a real edit of aperture.density_cm2 to 7.0e8 on a "
       "None-density design is written back as 7.0e8, not silently discarded to "
       "None -- within float32 log-slider quantisation (rel 1.05e-6: the widget "
       "stores log10(7.0e8) as float32, so 10.0**value on read-back does not "
       "reproduce 7.0e8 bit-exactly), and the untouched case still yields None")
def _():
    r = _run(["--roundtrip-check"])
    assert r.returncode == 0, f"exit {r.returncode}\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
    # untouched case (populated from a None-density design, never edited)
    assert "aperture.density_cm2=None" in r.stdout, r.stdout
    # edited case: the SAME widget, on the same run, edited to log10(7.0e8)
    # via its own callback (a real edit, not a value-only comparison)
    m = re.search(r"density_edit_to_7e8=([\d.eE+-]+)", r.stdout)
    assert m, f"density_edit_to_7e8= not found in stdout:\n{r.stdout}"
    edited = float(m.group(1))
    assert abs(edited - 7.0e8) / 7.0e8 < 1e-5, f"expected ~7.0e8, got {edited}\n{r.stdout}"


@check("verify/gate_v11_gui.py still exits 0 with its full N/N")
def _():
    r = _run_script(ROOT / "verify" / "gate_v11_gui.py")
    assert r.returncode == 0, f"exit {r.returncode}\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
    m = re.search(r"(\d+)/(\d+) checks passed", r.stdout)
    assert m and m.group(1) == m.group(2), r.stdout


# --------------------------------------------------------- (5)-(6) fsim_core untouched

@check("verify/verify_fsim.py still prints 51/51 (no fsim_core file modified)")
def _():
    r = _run_script(ROOT / "verify" / "verify_fsim.py")
    assert r.returncode == 0, f"exit {r.returncode}\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
    assert "51/51" in r.stdout, r.stdout


@check("verify/audit_physics.py still prints 23/23 (no fsim_core file modified)")
def _():
    r = _run_script(ROOT / "verify" / "audit_physics.py")
    assert r.returncode == 0, f"exit {r.returncode}\nSTDOUT:\n{r.stdout}\nSTDERR:\n{r.stderr}"
    assert "23/23" in r.stdout, r.stdout


def main():
    failed = 0
    for name, fn in CHECKS:
        try:
            fn()
            print(f"  PASS  {name}")
        except AssertionError as e:
            failed += 1
            print(f"* FAIL  {name}\n        {e}")
        except Exception as e:  # noqa: BLE001 -- surface subprocess/setup errors too
            failed += 1
            print(f"* ERROR {name}\n        {type(e).__name__}: {e}")
    n = len(CHECKS)
    print(f"\n{n - failed}/{n} designer RT checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
