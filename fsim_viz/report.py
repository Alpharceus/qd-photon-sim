"""Design-comparison report (phase D3): one-click bundle for 1-4 device
designs side by side. Consumes fsim-core result dicts only (three-layer
rule: matplotlib here, no GUI imports) -- the designer GUI is a thin client
that hands this module the same curves/scalars it already renders itself.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.patheffects
import matplotlib.pyplot as plt

import fsim_theme

from .figures import _band, _ref_h, _themed, _title, _write_csv

# Slot colours come from the token palette (light, report PNGs): A/B/C are the
# categorical comparison set (series slots 1-3, the same slots the designer
# GUI's SLOT_COLOR uses), "current" (the live/unstored run) is slot 7; any
# other label falls back to the remaining slots.
_SERIES = fsim_theme.load_tokens()["color"]["light"]["series"]
_FIXED_COLORS = {"current": _SERIES[6], "A": _SERIES[0], "B": _SERIES[1], "C": _SERIES[2]}
_FALLBACK_COLORS = [_SERIES[3], _SERIES[4], _SERIES[5], _SERIES[7]]

SCALAR_ROWS = ["T_j_op", "dT_J", "eps_op", "rho_op", "g2_op", "brightness_per_pulse",
              "T_c", "F_eff", "N_w", "aperture_g2_penalty"]

CURVE_HEADER = ["T_hs", "Tj", "g2", "eps", "rho2", "gamma"]


def _color_for(label: str, i: int) -> str:
    return _FIXED_COLORS.get(label, _FALLBACK_COLORS[i % len(_FALLBACK_COLORS)])


def _fnum(x) -> str:
    """Deterministic point-value formatting for a possibly-NaN float/bool/str."""
    if isinstance(x, bool):
        return str(x)
    if isinstance(x, float):
        return "nan" if x != x else f"{x:.6g}"
    return str(x)


def _fmt_scalar(entry: dict, name: str) -> str:
    sb = entry.get("scalar_bands")
    if sb and name in sb:
        lo, hi = sb[name]
        if lo != lo or hi != hi:  # NaN
            return "nan"
        return f"{lo:.6g}..{hi:.6g}"
    return _fnum(entry["scalars"].get(name, float("nan")))


def _write_curves_csv(path: Path, curves: dict, bands: dict | None) -> None:
    header = list(CURVE_HEADER)
    cols = [list(map(float, curves[k])) for k in CURVE_HEADER]
    if bands and "g2" in bands:
        lo, hi = bands["g2"]
        header += ["g2_lo", "g2_hi"]
        cols += [list(map(float, lo)), list(map(float, hi))]
    _write_csv(path, header, zip(*cols))


def _write_ranged_csv(path: Path, ranged: dict) -> None:
    _write_csv(path, ["path", "lo", "hi"],
              ((p, lo, hi) for p, (lo, hi) in ranged.items()))


def _write_scalars_csv(path: Path, entries: list) -> None:
    labels = [e["label"] for e in entries]
    _write_csv(path, ["scalar"] + labels,
              ([name] + [_fmt_scalar(e, name) for e in entries] for name in SCALAR_ROWS))


def _plot(entries: list, outdir: Path, title: str) -> None:
    """Three rows on one shared x (no twin axes): g2 with envelope bands and
    T_c references, dT_J, Gamma. Envelope rule from fsim_viz.figures._band."""
    with _themed("light") as P:
        fig, (ax1, ax2, ax3) = plt.subplots(
            3, 1, figsize=(8.5, 9.6), sharex=True,
            gridspec_kw={"height_ratios": [1.6, 1.0, 1.0]})
        n_banded = sum(1 for e in entries if e.get("bands") and "g2" in e["bands"])
        n_tc = sum(1 for e in entries if e.get("scalar_bands") and "T_c" in e["scalar_bands"])
        # stacked T_c intervals share one wash budget, so overlaps never go grey-black
        tc_wash = (*P.ref_wash[:3], P.ref_wash[3] / max(n_tc, 1))
        halo = [matplotlib.patheffects.withStroke(linewidth=2.5, foreground=P.surface)]
        g2_top = 1.0

        for i, e in enumerate(entries):
            label, c = e["label"], e["curves"]
            color = _color_for(label, i)
            Ts = c["T_hs"]

            bands = e.get("bands")
            if bands and "g2" in bands:
                lo, hi = bands["g2"]
                _band(ax1, Ts, lo, hi, color, P, overlap=n_banded > 1,
                      label=f"{label} band (swept inputs, lo/hi edges)")
                ax1.plot(Ts, c["g2"], color=color, lw=2,
                         label=f"{label}: {e['name']} (design at range midpoints)")
                g2_top = max(g2_top, max((float(v) for v in hi if v == v), default=0.0))
            else:
                ax1.plot(Ts, c["g2"], color=color, lw=2, label=f"{label}: {e['name']}")
            g2_top = max(g2_top, max((float(v) for v in c["g2"] if v == v), default=0.0))

            # T_c as a reference (not a verdict): interval for banded entries,
            # a rule for point entries; labels stacked so entries never collide
            sb = e.get("scalar_bands")
            y_lab = 0.97 - 0.075 * i
            if sb and "T_c" in sb:
                tlo, thi = sb["T_c"]
                if tlo == tlo and thi == thi:
                    ax1.axvspan(tlo, thi, color=tc_wash, lw=0, zorder=0.5)
                    txt, x_lab = f"{label}: T$_c$ {tlo:.0f}–{thi:.0f} K", thi
                else:
                    txt, x_lab = f"{label}: T$_c$ not crossed in range", None
            else:
                tc = e["scalars"].get("T_c", float("nan"))
                if tc == tc:
                    ax1.axvline(tc, color=P.ref, lw=P.ref_w, zorder=1.5)
                    txt, x_lab = f"{label}: T$_c$ = {tc:.0f} K", tc
                else:
                    txt, x_lab = f"{label}: T$_c$ not crossed in range", None
            if x_lab is None:
                ax1.annotate(txt, xy=(0.99, y_lab), xycoords="axes fraction", ha="right",
                             va="top", fontsize=8, color=P.ink2, path_effects=halo, zorder=8)
            else:
                ax1.annotate(txt, xy=(x_lab, y_lab), xycoords=("data", "axes fraction"),
                             xytext=(3, 0), textcoords="offset points", ha="left", va="top",
                             fontsize=8, color=P.ink1, path_effects=halo, zorder=8)

            dTj = [tj - t for tj, t in zip(c["Tj"], Ts)]
            ax2.plot(Ts, dTj, color=color, lw=1.8, label=f"{label}")
            ax3.plot(Ts, c["gamma"], color=color, lw=1.8, label=f"{label}")

        _ref_h(ax1, 0.5, "$g^{(2)}(0)$ = 0.5 ceiling", P, va="top")
        ax1.set_ylim(0, g2_top * 1.02 if g2_top > 1.0 else 1.0)
        ax1.set_ylabel("$g^{(2)}(0)$")
        ax1.legend(fontsize=7.5, loc="upper left")
        ax1.set_title("g$^{(2)}(0)$ vs heatsink T (shaded: swept-input band)", fontsize=10)

        ax2.set_ylabel("dT_J = T_j - T_hs (K)")
        ax2.set_title("junction overheat vs heatsink T", fontsize=10)
        ax3.set_ylabel(r"$\Gamma(T_j)$ (meV)")
        ax3.set_xlabel("heatsink T (K)")
        ax3.set_title("linewidth at the junction vs heatsink T", fontsize=10)
        if len(entries) > 1:
            ax2.legend(fontsize=7.5, ncol=min(len(entries), 4))

        fig.tight_layout()
        _title(fig, f"{title}   |   tag chain [A]: unmeasured inputs shape every curve above",
               "A", P, y=1.0)
        for ext in ("pdf", "svg", "png"):
            fig.savefig(outdir / f"report.{ext}", bbox_inches="tight", dpi=200)
        plt.close(fig)


def designer_report(entries: list, outdir, title: str = "design review") -> Path:
    """One-click comparison bundle for 1-N device designs.

    entries: [{"label", "name", "design" (DeviceDesign), "curves" (mid-or-point
    curves dict: T_hs/Tj/g2/eps/rho2/gamma), "bands" (optional, from
    evaluate_envelope), "scalars", "scalar_bands" (optional), "ranged"
    (optional dict of swept paths)}, ...].

    Writes report.pdf/svg/png (g2 vs T + dT_J/Gamma vs T, one line per entry),
    curves_<label>.csv, scalars_comparison.csv, design_<label>.yaml, and
    ranged_<label>.csv (only for entries carrying a "ranged" dict). Returns
    outdir.
    """
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    for e in entries:
        label = e["label"]
        _write_curves_csv(outdir / f"curves_{label}.csv", e["curves"], e.get("bands"))
        e["design"].save(outdir / f"design_{label}.yaml")
        if "ranged" in e:
            _write_ranged_csv(outdir / f"ranged_{label}.csv", e["ranged"])

    _write_scalars_csv(outdir / "scalars_comparison.csv", entries)
    _plot(entries, outdir, title)
    return outdir
