"""fsim_theme -- the single design-token source for every FSIM renderer.

`tokens.json` is the only place colours, fonts, sizes and mark specs live.
`python -m fsim_theme.build` emits from it:
  - fsim_studio/web/css/tokens.css          (web Studio CSS custom properties)
  - fsim_theme/fsim.mplstyle, fsim-dark.mplstyle (matplotlib)
  - fsim_theme/plotly_fsim_{light,dark}.json (Plotly templates)
and `fsim_theme.dpg.bind_theme()` binds a Dear PyGui theme at runtime.

The beam colour is reserved for the active signal path / active control: it is
never a data series and never a verdict. Provenance tags ([V]/[DR]/[E]/[A]) are
encoded by line form (TAG_FORM / tag_dash), never by hue.
"""
from __future__ import annotations

import json
from functools import lru_cache
from pathlib import Path

THEME_DIR = Path(__file__).resolve().parent
TOKENS_PATH = THEME_DIR / "tokens.json"
FONT_DIR = THEME_DIR / "fonts"
MODES = ("light", "dark")

_MPL_DASH = {"solid": "-", "dash": "--", "dot": ":"}


@lru_cache(maxsize=1)
def _tokens_cached() -> str:
    return TOKENS_PATH.read_text(encoding="utf-8")


def load_tokens() -> dict:
    """Return a fresh dict parsed from tokens.json (safe to mutate)."""
    return json.loads(_tokens_cached())


def _check_mode(mode: str) -> str:
    if mode not in MODES:
        raise ValueError(f"mode must be one of {MODES}, got {mode!r}")
    return mode


TAG_FORM = load_tokens()["tag"]


def _norm_tag(tag: str) -> str:
    t = str(tag).strip().strip("[]").upper()
    if t not in TAG_FORM:
        raise KeyError(f"unknown provenance tag {tag!r}; expected one of {list(TAG_FORM)}")
    return t


def tag_dash(tag: str, backend: str = "plotly") -> str:
    """Line dash for a provenance tag ('V', '[DR]', ...).

    backend="plotly" (default) returns the Plotly dash name (solid/dash/dot);
    backend="mpl" returns the matplotlib linestyle (-, --, :).
    """
    dash = TAG_FORM[_norm_tag(tag)]["dash"]
    if backend == "plotly":
        return dash
    if backend == "mpl":
        return _MPL_DASH[dash]
    raise ValueError(f"backend must be 'plotly' or 'mpl', got {backend!r}")


def tag_width(tag: str) -> float:
    """Line width (px) for a provenance tag."""
    return float(TAG_FORM[_norm_tag(tag)]["width"])


# ---------------------------------------------------------------- matplotlib
def mplstyle_path(mode: str = "light") -> Path:
    _check_mode(mode)
    return THEME_DIR / ("fsim.mplstyle" if mode == "light" else "fsim-dark.mplstyle")


def register_fonts_matplotlib() -> list[str]:
    """Add the vendored TTFs to matplotlib's font manager; returns paths added."""
    from matplotlib import font_manager

    added = []
    for face in load_tokens()["font"]["faces"]:
        ttf = face.get("ttf")
        if not ttf:
            continue
        p = FONT_DIR / ttf
        if p.is_file():
            font_manager.fontManager.addfont(str(p))
            added.append(str(p))
    return added


def register_colormaps(mode: str = "light") -> None:
    """Register `fsim_seq` and `fsim_div` (overwriting any earlier mode)."""
    import matplotlib
    from matplotlib.colors import LinearSegmentedColormap

    c = load_tokens()["color"][_check_mode(mode)]
    for name, colors in (("fsim_seq", c["seq"]), ("fsim_div", c["div"])):
        cmap = LinearSegmentedColormap.from_list(name, colors)
        if name in matplotlib.colormaps:
            matplotlib.colormaps.unregister(name)
        matplotlib.colormaps.register(cmap, name=name)


def apply_matplotlib(mode: str = "light") -> Path:
    """Register fonts, apply the fsim style for `mode`, register colormaps."""
    import matplotlib.pyplot as plt

    path = mplstyle_path(mode)
    register_fonts_matplotlib()
    plt.style.use(str(path))
    register_colormaps(mode)
    return path


# -------------------------------------------------------------------- plotly
def plotly_template_path(mode: str = "light") -> Path:
    return THEME_DIR / f"plotly_fsim_{_check_mode(mode)}.json"


def plotly_template(mode: str = "light"):
    """Return a plotly.graph_objects.layout.Template for `mode`."""
    import plotly.graph_objects as go

    data = json.loads(plotly_template_path(mode).read_text(encoding="utf-8"))
    return go.layout.Template(data)


def register_plotly(mode: str = "light", set_default: bool = False) -> str:
    """Register the template as plotly.io.templates['fsim'] (and 'fsim_<mode>').

    Returns the registered name. With set_default=True it also becomes
    plotly.io.templates.default.
    """
    import plotly.io as pio

    tpl = plotly_template(mode)
    pio.templates[f"fsim_{mode}"] = tpl
    pio.templates["fsim"] = tpl
    if set_default:
        pio.templates.default = "fsim"
    return "fsim"


__all__ = [
    "THEME_DIR", "TOKENS_PATH", "FONT_DIR", "MODES", "TAG_FORM",
    "load_tokens", "tag_dash", "tag_width",
    "mplstyle_path", "register_fonts_matplotlib", "register_colormaps", "apply_matplotlib",
    "plotly_template_path", "plotly_template", "register_plotly",
]
