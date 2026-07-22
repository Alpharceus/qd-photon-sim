"""Design-comparison report (phase D3): one-click bundle for 1-4 device
designs side by side. Consumes fsim-core result dicts only (three-layer
rule: matplotlib here, no GUI imports) -- the designer GUI is a thin client
that hands this module the same curves/scalars it already renders itself.
"""
from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from .figures import _write_csv

# fixed slot colors match the designer GUI's SLOT_COLOR/GREEN palette exactly
# (A amber, B blue, C purple, "current" the live/unstored run in green); any
# other label falls back to a small qualitative cycle.
_FIXED_COLORS = {"current": "#56a632", "A": "#e3a21a", "B": "#3c78d8", "C": "#785ac8"}
_FALLBACK_COLORS = ["#2166ac", "#e66101", "#5e3c99", "#d7191c", "#1a9641", "#e3a21a"]

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
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(8.5, 8.0))
    ax2b = ax2.twinx()

    for i, e in enumerate(entries):
        label, c = e["label"], e["curves"]
        color = _color_for(label, i)
        name_tag = f"{label}: {e['name']}"
        Ts = c["T_hs"]

        bands = e.get("bands")
        if bands and "g2" in bands:
            lo, hi = bands["g2"]
            ax1.fill_between(Ts, lo, hi, color=color, alpha=0.22, lw=0)
        ax1.plot(Ts, c["g2"], color=color, lw=2, label=name_tag)

        dTj = [tj - t for tj, t in zip(c["Tj"], Ts)]
        ax2.plot(Ts, dTj, color=color, lw=1.8, ls="-", label=f"{label} dT_J")
        ax2b.plot(Ts, c["gamma"], color=color, lw=1.8, ls="--", label=f"{label} Gamma")

    ax1.axhline(0.5, color="#888888", ls=":", lw=1)
    ax1.set_xlabel("heatsink T (K)")
    ax1.set_ylabel("$g^{(2)}(0)$")
    ax1.legend(frameon=False, fontsize=8)
    ax1.set_title("g$^{(2)}(0)$ vs heatsink T (shaded: swept-input band)", fontsize=10)

    ax2.set_xlabel("heatsink T (K)")
    ax2.set_ylabel("dT_J = T_j - T_hs (K)  [solid]")
    ax2b.set_ylabel(r"$\Gamma(T_j)$ (meV)  [dashed]")
    h1, l1 = ax2.get_legend_handles_labels()
    h2, l2 = ax2b.get_legend_handles_labels()
    ax2.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8, ncol=2)
    ax2.set_title("junction overheat and linewidth vs heatsink T", fontsize=10)

    fig.suptitle(f"{title}   |   tag chain [A]: unmeasured inputs shape every curve above",
                fontsize=11, color="#d7191c", y=1.0)
    fig.tight_layout()
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
