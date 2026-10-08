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
    # work-order Q4: DotBlock.width_kind is validated ("total" | "zpl"), so the generic
    # string bump ("total-probe") would be rejected the same way.
    "dot.width_kind": "zpl",  # default "total"
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
        if not dataclasses.is_dataclass(exp_block):
            # Non-dataclass top-level fields (e.g. `platform` str, `nitride`
            # dict added with the nitride tier) are compared exactly;
            # build_non_default_design() leaves them at their defaults, so
            # they must round-trip unchanged.
            _cmp_value(block_name, exp_block, act_block, False, errors, tmp_out)
            continue
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


# ------------------------------------- (7)-(9) static (non-headline) labelling
#
# Spec audit-edge-cards-label (user decision Q2, 2026-09-23): the RT edge
# cards keep drive.finite_pulse=false; every card-level g2/brightness the
# designer emits for them must carry "static (non-headline)" and name the
# headline finite-pulse g2 in the same panel. The real _run_point /
# _run_envelope code runs headlessly in a subprocess (this file still never
# imports dearpygui); only the T grid is narrowed to the card's own T_hs
# (g2_op/brightness are evaluated AT T_hs either way, so the emitted
# operating-point numbers are identical to a full-grid run).

EDGE_CARDS = {
    "edge-inp-gainp-design": ROOT / "cards" / "edge-inp-gainp-design.yaml",
    "edge-inp-gaasp-design": ROOT / "cards" / "edge-inp-gaasp-design.yaml",
}
# Card-level (static, drive.finite_pulse=false) g2_op and brightness_per_pulse
# at each card's own thermal.T_hs, captured with DeviceDesign.load + evaluate
# on 2026-09-23 immediately BEFORE the labelling edit (after the reviewed
# waveguide Purcell-area and aperture cross-term fixes). The labelling must
# not move them.
PRE_LABEL_CARD_VALUES = {
    "edge-inp-gainp-design": {"g2_op": 0.7698721127665349,
                              "brightness_per_pulse": 1.4199538816091535e-06},
    "edge-inp-gaasp-design": {"g2_op": 0.9787218623195735,
                              "brightness_per_pulse": 5.41092896804776e-08},
}
# Headline finite-pulse g2_op at the gainp card point (300 K): finite_pulse
# gives g2_op 0.99911 against 0.77363 for the static model; the static half
# has since moved -0.004 with the aperture cross-term fix, the finite-pulse
# half is saturated near 1.
AUDIT_HEADLINE_G2_GAINP_300K = 0.99911

_PROBE = r"""
import sys, warnings
warnings.filterwarnings("ignore")
sys.path.insert(0, sys.argv[1]); sys.path.insert(0, sys.argv[1] + "/fsim_gui")
import designer as g
import dearpygui.dearpygui as dpg
from fsim_core.device import DeviceDesign
_ev, _env = g.evaluate, g.evaluate_envelope
g.evaluate = lambda d, T_grid=None: _ev(d, T_grid=[d.thermal.T_hs])
g.evaluate_envelope = lambda d, r, T_grid=None, **kw: _env(d, r, T_grid=[d.thermal.T_hs], **kw)
dpg.create_context(); g.build_ui()
for spec in sys.argv[2:]:
    mode, path = spec.split("=", 1)
    g.apply_design(DeviceDesign.load(path))
    if mode == "point":
        g._run_point(g.collect_design())
    else:
        g._run_envelope(g.collect_design(), g._collect_ranged())
    print("<<<" + mode + " " + path + ">>>")
    print(dpg.get_value("results_text"))
dpg.destroy_context()
"""

_PROBE_CACHE = {}


def _probe(specs):
    key = tuple(specs)
    if key not in _PROBE_CACHE:
        r = subprocess.run([PY, "-c", _PROBE, str(ROOT), *specs], cwd=str(ROOT),
                           capture_output=True, text=True, timeout=900)
        assert r.returncode == 0, f"probe exited {r.returncode}\nSTDERR:\n{r.stderr[-3000:]}"
        panels = {}
        for chunk in r.stdout.split("<<<")[1:]:
            head, _, body = chunk.partition(">>>")
            panels[head.strip()] = body
        _PROBE_CACHE[key] = panels
    return _PROBE_CACHE[key]


def _edge_panels():
    specs = [f"point={p}" for p in EDGE_CARDS.values()]
    specs.append(f"envelope={EDGE_CARDS['edge-inp-gainp-design']}")
    return _probe(specs)


_HEADLINE_RE = re.compile(
    r"headline model = drive\.finite_pulse=true .*headline finite-pulse g2_op at T_hs "
    r"\d+ K = ([0-9.]+)")


@check("audit-edge-cards-label: designer point output for BOTH RT edge cards labels "
       "g2(0) and brightness 'static (non-headline)' and names the headline "
       "finite-pulse g2 (drive.finite_pulse=true) in the same panel; the gainp "
       "envelope output carries the same label; the gainp headline g2 equals the "
       "device audit's independent 0.99911 (1e-3)")
def _():
    panels = _edge_panels()
    for name, path in EDGE_CARDS.items():
        text = panels[f"point {path}"]
        for key in ("g2(0) at operating point", "brightness/pulse (t_X)"):
            line = next((ln for ln in text.splitlines() if ln.startswith(key)), "")
            assert "[static (non-headline)]" in line, f"{name}: {key!r} line unlabelled: {line!r}"
        m = _HEADLINE_RE.search(text)
        assert m, f"{name}: no headline-model pointer with a numeric g2 in:\n{text}"
        if name == "edge-inp-gainp-design":
            assert abs(float(m.group(1)) - AUDIT_HEADLINE_G2_GAINP_300K) < 1e-3, m.group(0)
    env = panels[f"envelope {EDGE_CARDS['edge-inp-gainp-design']}"]
    line = next((ln for ln in env.splitlines() if ln.startswith("g2(0) at op:")), "")
    assert "[static (non-headline)]" in line, f"envelope g2 line unlabelled: {line!r}"
    assert _HEADLINE_RE.search(env), f"envelope panel has no headline pointer:\n{env}"


@check("audit-edge-cards-label: labelling changes no number -- card-level g2_op and "
       "brightness_per_pulse of both edge cards equal the values captured before the "
       "edit (rel 1e-12), and the designer's printed g2(0) equals them to its 3 decimals")
def _():
    import warnings
    from fsim_core.device import evaluate
    panels = _edge_panels()
    for name, path in EDGE_CARDS.items():
        d = DeviceDesign.load(path)
        assert d.drive.finite_pulse is False, f"{name}: drive.finite_pulse changed"
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            sc = evaluate(d, T_grid=[d.thermal.T_hs])["scalars"]
        for key, pinned in PRE_LABEL_CARD_VALUES[name].items():
            got = float(sc[key])
            assert abs(got - pinned) <= 1e-12 * abs(pinned), f"{name}: {key} {got!r} != {pinned!r}"
        m = re.search(r"g2\(0\) at operating point ([0-9.]+)", panels[f"point {path}"])
        assert m and m.group(1) == f"{PRE_LABEL_CARD_VALUES[name]['g2_op']:.3f}", (
            f"{name}: designer printed {m.group(1) if m else None}")


@check("audit H5 / audit-edge-cards-label task 4: the F_eff result line is labelled as "
       "the single-mode F_P input, never a planar-DBR total rate -- planar-lambda "
       "preset design shows 'F_eff (mode-only F_P [A]; planar DBR total rate ~1.0)', "
       "edge cards show 'F_eff (single-mode F_P input, not a planar total rate)'; a "
       "non-edge (legacy) design gets no static label")
def _():
    from fsim_core.presets import preset_device
    with tempfile.TemporaryDirectory() as td:
        planar = Path(td) / "planar-lambda.yaml"
        preset_device("staged-inp-gaasp", "GaAs", "planar-lambda", "cw-electrical").save(planar)
        panels = _probe([f"point={planar}"])
    text = panels[f"point {planar}"]
    assert "F_eff (mode-only F_P [A]; planar DBR total rate ~1.0)" in text, text
    assert "static (non-headline)" not in text, "legacy planar design wrongly labelled static"
    for name, path in EDGE_CARDS.items():
        assert ("F_eff (single-mode F_P input, not a planar total rate)"
                in _edge_panels()[f"point {path}"]), name


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
