"""verify_theme -- checks for the fsim_theme token package and its emitters.

Run: python verify/verify_theme.py   (exit 0 iff all PASS)

Checks: tokens.json parses with identical light/dark key sets; every colour
hex is valid; generated files are byte-identical to a fresh build into a temp
dir; the matplotlib style loads and saves a PNG with Barlow resolved from
fsim_theme/fonts; the Plotly templates load into plotly.io.templates; the
Dear PyGui theme binds in a headless context.
"""
from __future__ import annotations

import re
import sys
import tempfile
import traceback
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import fsim_theme  # noqa: E402
from fsim_theme import build as B  # noqa: E402
from fsim_theme import dpg as fsim_theme_dpg  # noqa: E402

RESULTS = []
HEX = re.compile(r"^#[0-9a-fA-F]{6}$")


def check(name, fn):
    try:
        detail = fn()
        ok = True
    except Exception as exc:  # report, never crash
        ok = False
        detail = f"{type(exc).__name__}: {exc}"
        tb = traceback.format_exc().strip().splitlines()[-3:]
        detail += " | " + " / ".join(s.strip() for s in tb)
    RESULTS.append(ok)
    print(f"  {'PASS' if ok else 'FAIL'}  {name}" + (f"\n        {detail}" if detail else ""))


def _keys(d, prefix=""):
    out = set()
    for k, v in d.items():
        out.add(prefix + k)
        if isinstance(v, dict):
            out |= _keys(v, prefix + k + ".")
        elif isinstance(v, list):
            out.add(f"{prefix}{k}[{len(v)}]")
    return out


def c_tokens_modes():
    tok = fsim_theme.load_tokens()
    light, dark = tok["color"]["light"], tok["color"]["dark"]
    kl, kd = _keys(light), _keys(dark)
    assert kl == kd, f"key mismatch: light-only {sorted(kl - kd)}, dark-only {sorted(kd - kl)}"
    need = {"ground", "plate", "rim", "ink-1", "ink-2", "muted", "grid", "axis", "ref",
            "ref-wash", "beam", "focus", "surface", "series[8]", "seq[7]", "div[3]",
            "ordinal_tag[4]", "status.pass", "status.warn", "status.fail", "status.info"}
    missing = need - kl
    assert not missing, f"missing keys {sorted(missing)}"
    for sec in ("font", "size", "mark", "tag", "radius", "space"):
        assert sec in tok, f"missing section {sec}"
    assert set(tok["tag"]) == {"V", "DR", "E", "A"}
    return f"{len(kl)} keys per mode"


def c_hex():
    tok = fsim_theme.load_tokens()
    n = 0
    for mode in fsim_theme.MODES:
        c = tok["color"][mode]
        for k, v in c.items():
            vals = v if isinstance(v, list) else (list(v.values()) if isinstance(v, dict) else [v])
            for h in vals:
                if k == "ref-wash":
                    assert h.startswith("rgba("), f"{mode}.{k} {h}"
                    continue
                assert HEX.match(h), f"{mode}.{k} invalid hex {h!r}"
                n += 1
    return f"{n} hex values valid"


def c_series_verbatim():
    # 02-charts.md section 4, validated palette (light / dark), verbatim.
    light = "#2a78d6 #eb6834 #1baf7a #eda100 #e87ba4 #008300 #4a3aa7 #e34948".split()
    dark = "#3987e5 #d95926 #199e70 #c98500 #d55181 #008300 #9085e9 #e66767".split()
    tok = fsim_theme.load_tokens()
    assert tok["color"]["light"]["series"] == light
    assert tok["color"]["dark"]["series"] == dark
    assert tok["color"]["dark"]["beam"].upper() == "#FF3B2F"
    assert tok["color"]["light"]["beam"].upper() == "#D2281E"
    assert tok["color"]["dark"]["beam"].lower() not in dark
    assert tok["color"]["light"]["beam"].lower() not in light


def c_generated_in_sync():
    with tempfile.TemporaryDirectory() as td:
        B.build(td)
        bad = []
        for rel in B.generate():
            a = (ROOT / rel)
            b = Path(td) / rel
            if not a.is_file() or a.read_bytes() != b.read_bytes():
                bad.append(rel)
        assert not bad, f"out of sync (run python -m fsim_theme.build): {bad}"
    return f"{len(B.generate())} generated files in sync"


def c_tag_helpers():
    assert fsim_theme.tag_dash("V") == "solid"
    assert fsim_theme.tag_dash("[DR]") == "solid"
    assert fsim_theme.tag_dash("E") == "dash"
    assert fsim_theme.tag_dash("A") == "dot"
    assert fsim_theme.tag_dash("A", backend="mpl") == ":"
    assert fsim_theme.TAG_FORM["E"]["chip"] == "outline"


def c_matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib import font_manager
    out = []
    for mode in fsim_theme.MODES:
        plt.style.use("default")
        fsim_theme.apply_matplotlib(mode)
        fpath = Path(font_manager.findfont(
            font_manager.FontProperties(family="Barlow"), fallback_to_default=False))
        assert fpath.resolve().parent == fsim_theme.FONT_DIR.resolve(), f"Barlow -> {fpath}"
        assert "fsim_seq" in matplotlib.colormaps and "fsim_div" in matplotlib.colormaps
        surf = fsim_theme.load_tokens()["color"][mode]["surface"].lower()
        assert matplotlib.colors.to_hex(plt.rcParams["figure.facecolor"]) == surf
        fig, ax = plt.subplots(figsize=(3, 2))
        for i in range(3):
            ax.plot([0, 1, 2], [i, i + 1, i], label=f"s{i}")
        ax.set_title("theme test")
        ax.set_ylabel("g$^{(2)}$(0)")
        ax.legend()
        with tempfile.TemporaryDirectory() as td:
            p = Path(td) / f"t_{mode}.png"
            fig.savefig(p)
            assert p.stat().st_size > 1000
        plt.close(fig)
        out.append(f"{mode}: {fpath.name}")
    plt.style.use("default")
    return ", ".join(out)


def c_plotly():
    import plotly.io as pio
    import plotly.graph_objects as go
    for mode in fsim_theme.MODES:
        name = fsim_theme.register_plotly(mode)
        assert name in pio.templates and f"fsim_{mode}" in pio.templates
        tpl = pio.templates[name]
        c = fsim_theme.load_tokens()["color"][mode]
        assert list(tpl.layout.colorway) == c["series"]
        assert tpl.layout.hovermode == "x unified"
        fig = go.Figure(go.Scatter(x=[0, 1], y=[1, 0]), layout={"template": name})
        fig.to_dict()
    return "light + dark registered"


def c_dpg():
    import dearpygui.dearpygui as dpg
    out = []
    for mode in ("dark", "light"):
        dpg.create_context()
        try:
            ids = fsim_theme_dpg.bind_theme(dpg, mode=mode)
            assert ids["theme"] is not None
            assert ids["ui"] is not None and ids["mono"] is not None, ids
        finally:
            dpg.destroy_context()
        out.append(mode)
    return "bound headless: " + ", ".join(out)


def main() -> int:
    print("verify_theme: fsim_theme tokens and emitters")
    check("tokens.json parses; light/dark key sets identical", c_tokens_modes)
    check("every colour token is a valid hex", c_hex)
    check("series palette verbatim from 02-charts.md; beam reserved", c_series_verbatim)
    check("generated files in sync with tokens (temp rebuild byte-compare)", c_generated_in_sync)
    check("TAG_FORM / tag_dash helpers", c_tag_helpers)
    check("matplotlib style loads, Barlow resolved from fsim_theme/fonts, PNG saves", c_matplotlib)
    check("Plotly templates load into plotly.io.templates", c_plotly)
    check("Dear PyGui theme binds in a headless context", c_dpg)
    n, k = len(RESULTS), sum(RESULTS)
    print(f"{k}/{n} theme checks passed")
    return 0 if k == n else 1


if __name__ == "__main__":
    sys.exit(main())
