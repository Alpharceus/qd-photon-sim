"""POST /api/compare/report: hand cached runs to fsim_viz.report.designer_report
unchanged (bundle written under .cache/studio/reports/, never under out/)."""
from __future__ import annotations

import hashlib
import time
from pathlib import Path

import numpy as np

from . import CACHE_DIR, cache
from .serialize import nan_from_null

_CURVES = ("T_hs", "Tj", "g2", "eps", "rho2", "gamma")


def _arr(v):
    return np.asarray(nan_from_null(v if v is not None else []), dtype=float)


def _num(v):
    return float("nan") if v is None else v


def build_report(entries: list, title: str = "design review") -> Path:
    from fsim_viz.report import designer_report

    from .api_cards import design_from_dict
    from .jobs import headline_design

    rep_entries = []
    for e in entries:
        label, run_id = str(e.get("label", "")), str(e.get("run_id", ""))
        rec = cache.find_run(run_id)
        if rec is None:
            raise LookupError(f"no cached run {run_id!r} (run it first)")
        r = rec["result"]
        design = design_from_dict(rec["design"])
        if rec.get("mode") == "headline":
            design = headline_design(design)
        entry = {"label": label or run_id, "name": rec.get("card") or design.name,
                 "design": design,
                 "curves": {k: _arr(r["curves"].get(k)) for k in _CURVES},
                 "scalars": {k: _num(v) for k, v in r.get("scalars", {}).items()}}
        if r.get("bands"):
            entry["bands"] = {k: (_arr(b["lo"]), _arr(b["hi"])) for k, b in r["bands"].items()}
        if r.get("scalar_bands"):
            entry["scalar_bands"] = {k: (_num(v[0]), _num(v[1]))
                                     for k, v in r["scalar_bands"].items()}
        if r.get("ranged"):
            entry["ranged"] = {k: (min(v), max(v)) for k, v in r["ranged"].items()}
        rep_entries.append(entry)
    tag = hashlib.sha256("|".join(str(e.get("run_id")) for e in entries).encode()).hexdigest()[:6]
    outdir = CACHE_DIR / "reports" / f"{time.strftime('%Y%m%d-%H%M%S')}-{tag}"
    return Path(designer_report(rep_entries, outdir, title=title))
