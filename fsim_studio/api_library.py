"""Library workspace backend: the literature parameter cards (cards/*.yaml without a
`design:` block). Configuration and bookkeeping only (three-layer rule): every
parameter and data row comes from the card file through fsim_core.card.load_card,
and every model overlay is the return value of the named fsim_core function
(see OVERLAYS).

Routes (registered by register(app) from server.create_app):
  GET /api/params                          -> [{name, title, paper, device, source, role, cls,
                                                tag_counts, datasets: [{name, rows, tag}], ...}]
  GET /api/params/<name>                   -> {name, meta, params, datasets, tag_counts, ...}
  GET /api/params/<name>/overlay/<dataset> -> model curve per the rules (404 + reason otherwise)
  PUT /api/params/<name>                   -> writes cards/<name>.yaml; <name> must end in
                                              "-edited" (409 for a shipped card name)
"""
from __future__ import annotations

import csv
import math
import re
from pathlib import Path

import numpy as np
import yaml
from flask import Response, request

from . import CARDS, OUT
from .serialize import dumps

# Patched by verify/verify_studio_library.py to a temporary copy (PUT round-trip
# without writing into the repository's cards/).
CARDS_DIR: Path = CARDS

NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,120}$")
TAGS = ("V", "DR", "E", "A")
EDITED = "-edited"

# Display grouping (bibliographic configuration, not physics): the paper each card reads.
PAPERS = [
    {"id": "chatzarakis2023", "label": "Chatzarakis et al. 2023",
     "cite": "Phys. Rev. Applied 20, 034011 (2023)", "cards": ["chatzarakis", "chatzarakis2023"]},
    {"id": "laferriere2023", "label": "Laferrière et al. 2023",
     "cite": "Nano Lett. 23, 962 (2023)", "cards": ["laferriere", "laferriere2023"]},
    {"id": "reischle2010", "label": "Reischle et al. 2010 / 2008",
     "cite": "Appl. Phys. Lett. 97, 143513 (2010); Opt. Express 16, 12771 (2008)", "cards": ["reischle"]},
    {"id": "zhao2024", "label": "Zhao et al. 2024",
     "cite": "Phys. Rev. Research 6, L032021 (2024)", "cards": ["zhao"]},
    {"id": "kitamura2025", "label": "Kitamura et al. 2025",
     "cite": "CLEO 2025", "cards": ["kitamura"]},
    {"id": "qcap", "label": "QCAP program cards",
     "cite": "no paper: requirement and prediction cards (anchors cited per entry)",
     "cards": ["qcap-cavity", "qcap-staged", "qcap-piezo-variant"]},
]

# Design cards that derive from a literature card ("Open as design").
DERIVED_DESIGN = {"qcap-staged": "staged-device-design"}

# Card class per the README "Validation record" and the physics brief section 8: C calibration
# (fitted to its own data), M unfitted comparison, T transcription (anchor / formula card),
# P prediction card. Labels and reasons are the rules file's; the overlay registry follows the
# overlay builders below (OVERLAYS is filled there, NO_OVERLAY with it).
CLASS_LABEL = {"C": "calibration", "M": "unfitted comparison", "T": "transcription", "P": "prediction"}
CLASSES = {
    "chatzarakis": {"code": "C", "label": "calibration", "reason":
                    "V-a: the model is fitted to this card's own g2/Gamma/tau series (tolerance V_A_TOL 0.03); "
                    "agreement with the points is by construction, not a held-out test"},
    "chatzarakis2023": {"code": "C", "label": "calibration (derived)", "reason":
                        "V1 device tier recomposed at the V-a fitted [DR] parameters with the published windows; "
                        "inherits the calibration, so it is not a held-out test"},
    "laferriere": {"code": "M", "label": "unfitted comparison", "reason":
                   "V-b: nothing is fitted to this card; the model envelope over unpublished ranges is compared "
                   "with the published points"},
    "laferriere2023": {"code": "M", "label": "unfitted comparison", "reason":
                       "V1 device tier: an envelope over the 16 corners of the unpublished [E] ranges, compared "
                       "with the published points; nothing is fitted"},
    "reischle": {"code": "M", "label": "unfitted comparison", "reason":
                 "V-c: trion line, eps = 0 structurally, g2 = 1 - rho^2; the spectral rho is compared with the rho "
                 "the measured g2 requires; nothing is fitted"},
    "kitamura": {"code": "T", "label": "transcription", "reason":
                 "an anchor (EL/PL linewidth ratio 4 [V]) and a formula card; there is no dataset to compare"},
    "zhao": {"code": "C", "label": "calibration, verdict FAIL", "reason":
             "Tier-2 Langevin fit to the published g2 values: the normal pump is 1.7 sigma off and the quiet pump "
             "5.8 sigma off (verdict FAIL on both rows); a QD LASER model, not a single-photon source"},
    "qcap-cavity": {"code": "P", "label": "prediction", "reason":
                    "requirement envelope, every input [A], awaiting MEEP/COMSOL; no data"},
    "qcap-staged": {"code": "P", "label": "prediction", "reason":
                    "prediction card for the staged InP/GaAsP device; no data"},
    "qcap-piezo-variant": {"code": "P", "label": "prediction", "reason":
                           "(Delta, rho) requirement map for T_c >= 300 K over [A]-swept inputs, with measured "
                           "devices placed on it"},
}
OVERLAYS = {}       # (card, dataset) -> callable(card_obj, raw_doc) -> overlay dict
NO_OVERLAY = {}     # (card, dataset) or card -> reason (from the rules file)

# Axis hints for the chart system (which columns are x / y / error / bound). Labels only.
PLOTS = {
    "g2_vs_T": {"x": "T", "y": ["g2"], "err": "err", "bound": "bound", "x_label": "T (K)", "y_label": "g²(0)"},
    "gamma_vs_T": {"x": "T", "y": ["gamma_x", "gamma_xx"], "x_label": "T (K)", "y_label": "FWHM Γ (meV)"},
    "tau_vs_T": {"x": "T", "y": ["tau_ns"], "x_label": "T (K)", "y_label": "X decay time τ (ns)"},
    "g2_vs_ERR": {"x": "ERR_MHz", "y": ["g2"], "err": "err", "group": ["device", "position"],
                  "x_label": "excitation repetition rate (MHz)", "y_label": "g²(0)"},
    "g2_oe2008": {"x": "T", "y": ["g2_raw", "g2_b", "g2_corr"], "err_for": {"g2_b": "err"},
                  "x_label": "T (K)", "y_label": "g²(0)"},
    "rho_spectral": {"x": "w_meV", "interval": ["rho_lo", "rho_hi"],
                     "x_label": "window full width w (meV)", "y_label": "signal fraction ρ(w)"},
    "g2_vs_pump": {"x": "pump", "y": ["g2"], "err": "err", "categorical": True, "zero": False,
                   "x_label": "pump", "y_label": "g²(0)"},
    "device_points": {"x": "delta", "y": ["rho"], "label": "name", "sentinel": {"rho": -1.0},
                      "x_label": "Δ_XX (meV)", "y_label": "ρ"},
}


class LibraryError(Exception):
    def __init__(self, status: int, message: str, **extra):
        super().__init__(message)
        self.status = status
        self.message = message
        self.extra = extra


# ------------------------------------------------------------------ discovery + loading
def _path(name: str) -> Path:
    if not NAME_RE.match(name or "") or name.endswith(".yaml"):
        raise LibraryError(400, f"invalid card name {name!r}")
    return CARDS_DIR / f"{name}.yaml"


def _raw(path: Path) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def _is_param_doc(doc: dict) -> bool:
    # A card whose meta.role names a design (fock-source-4K-design is loaded by resonant_source,
    # not fsim_core.card) is a design-tier card even without a `design:` block.
    return (isinstance(doc, dict) and not isinstance(doc.get("design"), dict) and isinstance(doc.get("meta"), dict)
            and not str(doc["meta"].get("role", "")).endswith("-design"))


def param_names() -> list:
    out = []
    for p in sorted(CARDS_DIR.glob("*.yaml")):
        try:
            doc = _raw(p)
        except (OSError, yaml.YAMLError):
            continue
        if _is_param_doc(doc):
            out.append(p.stem)
    return out


def base_of(name: str) -> str:
    return name[: -len(EDITED)] if name.endswith(EDITED) else name


def _load(name: str):
    """(card or None, raw doc, loader note). fsim_core.card.load_card is the loader; a card
    it cannot load at all is read raw and the note says so."""
    from fsim_core.card import load_card
    path = _path(name)
    if not path.is_file():
        raise LibraryError(404, f"no card {name!r}")
    doc = _raw(path)
    if not _is_param_doc(doc):
        raise LibraryError(422, f"card {name!r} is a design card, not a parameter card")
    try:
        return load_card(path), doc, "fsim_core.card.load_card"
    except (KeyError, ValueError, TypeError) as exc:
        return None, doc, (f"read raw: fsim_core.card.load_card raises {type(exc).__name__}: {exc} "
                           "(datasets without a per-dataset source)")


def _num(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool)


def _conflict_of(doc, name):
    """The card's own `conflict:` note on a parameter (a known disagreement between sources), verbatim."""
    raw = (doc.get("params") or {}).get(name)
    c = raw.get("conflict") if isinstance(raw, dict) else None
    return " ".join(str(c).split()) if c else None


def _param_entries(card, doc) -> list:
    if card is not None:
        return [{"name": n, "tag": p.tag.name, "unit": p.unit, "source": p.source,
                 "value": p.value, "range": list(p.range) if p.range is not None else None,
                 "conflict": _conflict_of(doc, n)}
                for n, p in card.params.items()]
    from fsim_core.card import _parse_param
    out = []
    for n, raw in (doc.get("params") or {}).items():
        p = _parse_param(n, raw)
        out.append({"name": n, "tag": p.tag.name, "unit": p.unit, "source": p.source,
                    "value": p.value, "range": list(p.range) if p.range is not None else None,
                    "conflict": _conflict_of(doc, n)})
    return out


def _dataset_entries(card, doc) -> list:
    meta_src = str((doc.get("meta") or {}).get("source", ""))
    out = []
    for n, raw in (doc.get("data") or {}).items():
        ds = card.datasets.get(n) if card is not None else None
        if ds is not None:
            tag, source, rows, cols = ds.tag.name, ds.source, ds.rows, ds.columns
            src_note = None
            if raw.get("source") is None:  # load_card fell back to meta.source (device-tier cards)
                src_note = "the dataset carries no source of its own; this is the card's meta.source"
        else:
            rows = [dict(r) for r in raw.get("points") or []]
            cols = sorted({k for r in rows for k in r})
            tag = str(raw["tag"])
            source = raw.get("source")
            src_note = None
            if source is None:
                source = meta_src
                src_note = "the dataset carries no source of its own; this is the card's meta.source"
        excluded = [dict(r) for r in raw.get("excluded") or []]
        out.append({"name": n, "tag": tag, "source": " ".join(str(source).split()),
                    "source_note": src_note, "columns": cols, "rows": rows, "n_rows": len(rows),
                    "excluded": excluded, "plot": PLOTS.get(n)})
    return out


_DEVICE_LINE = re.compile(r"^(\s+)([A-Za-z_][\w]*|\"[^\"]+\")\s*:\s*(.*?)\s*(?:#\s*(.*))?$")
_TAG_IN = re.compile(r"\[(V|DR|E|A)(?:[ /][^\]]*)?\]")


def _device_block(name: str, doc: dict) -> list:
    """Device-tier configuration rows (cards with a `device:` block). Values come from the
    parsed YAML; the provenance tag is the bracket tag the card itself writes in that line's
    comment (or the nearest tagged parent line). A number with no tag in the card is listed
    with tag None so the UI withholds it (CONTRACT rule 2)."""
    dev = doc.get("device")
    if not isinstance(dev, dict):
        return []
    text = _path(name).read_text(encoding="utf-8").splitlines()
    comments = {}
    stack = []  # (indent, key, tag)
    inside = False
    for line in text:
        if line.startswith("device:"):
            inside = True
            continue
        if inside and line and not line[0].isspace():
            break
        if not inside:
            continue
        m = _DEVICE_LINE.match(line)
        if not m:
            continue
        ind, key, _, com = len(m.group(1)), m.group(2).strip('"'), m.group(3), m.group(4) or ""
        while stack and stack[-1][0] >= ind:
            stack.pop()
        t = _TAG_IN.search(com)
        tag = t.group(1) if t else (stack[-1][2] if stack else None)
        path = ".".join([s[1] for s in stack] + [key])
        comments[path] = (tag, com.strip())
        stack.append((ind, key, tag))

    rows = []
    hidden = []

    def walk(prefix, obj):
        if isinstance(obj, dict):
            for k, v in obj.items():
                walk(f"{prefix}.{k}" if prefix else str(k), v)
            return
        leaf = prefix.split(".")[-1]
        if leaf == "mid" or leaf.startswith("mid_"):
            hidden.append(prefix)  # a range's centre: ranges are never collapsed to a value
            return
        tag, com = comments.get(prefix, (None, ""))
        parts = prefix.split(".")
        while tag is None and len(parts) > 1:  # inherit the nearest tagged parent line
            parts = parts[:-1]
            tag, pcom = comments.get(".".join(parts), (None, ""))
            com = com or pcom
        parts = prefix.split(".")
        while not com and len(parts) > 1:  # the comment that carries the inherited tag
            parts = parts[:-1]
            com = comments.get(".".join(parts), (None, ""))[1]
        if isinstance(obj, list) and len(obj) == 2 and all(_num(x) for x in obj) and prefix.split(".")[-1].startswith("range"):
            rows.append({"path": prefix, "range": [float(obj[0]), float(obj[1])], "value": None,
                         "tag": tag, "comment": com})
        elif _num(obj):
            rows.append({"path": prefix, "value": float(obj), "range": None, "tag": tag, "comment": com})
        else:
            rows.append({"path": prefix, "text": str(obj), "value": None, "range": None, "tag": tag, "comment": com})
    walk("", dev)
    if hidden:
        rows.append({"path": "(range centres)", "hidden": hidden, "value": None, "range": None, "tag": None,
                     "comment": "the card's mid values are range centres; ranges are swept, never collapsed, so they are not shown"})
    return rows


def tag_counts(params: list) -> dict:
    counts = {t: 0 for t in TAGS}
    for p in params:
        counts[p["tag"]] += 1
    return counts


def _paper_of(name: str) -> dict | None:
    b = base_of(name)
    for i, p in enumerate(PAPERS):
        if b in p["cards"]:
            return {**{k: p[k] for k in ("id", "label", "cite")}, "order": i}
    return None


def _class_of(name: str) -> dict:
    c = CLASSES.get(base_of(name))
    if not c:
        return {"code": None, "label": "unclassified", "reason": "the physics rules file does not classify this card yet"}
    return dict(c)


def _summary(name: str, card, doc, loader: str) -> dict:
    meta = doc.get("meta") or {}
    params = _param_entries(card, doc)
    dsets = _dataset_entries(card, doc)
    return {
        "name": name, "base": base_of(name), "edited": name.endswith(EDITED),
        "title": meta.get("name", name), "device": " ".join(str(meta.get("device", "")).split()) or None,
        "role": meta.get("role"), "source": " ".join(str(meta.get("source", "")).split()) or None,
        "paper": _paper_of(name), "cls": _class_of(name),
        "tag_counts": tag_counts(params), "n_params": len(params),
        "n_conflicts": sum(1 for p in params if p.get("conflict")),
        "datasets": [{"name": d["name"], "tag": d["tag"], "n_rows": d["n_rows"]} for d in dsets],
        "dataset_tag_counts": tag_counts(dsets),
        "loader": loader, "design": DERIVED_DESIGN.get(base_of(name)),
        "has_device_block": isinstance(doc.get("device"), dict),
        "device_tag_counts": tag_counts([r for r in _device_block(name, doc) if r.get("tag") in TAGS]),
        "n_device_fields": sum(1 for r in _device_block(name, doc) if r.get("tag") in TAGS),
    }


def list_params() -> list:
    out = []
    for n in param_names():
        try:
            card, doc, loader = _load(n)
            out.append(_summary(n, card, doc, loader))
        except (LibraryError, ValueError, KeyError) as exc:
            out.append({"name": n, "error": str(exc)})
    return out


def get_params(name: str) -> dict:
    card, doc, loader = _load(name)
    out = _summary(name, card, doc, loader)
    meta = doc.get("meta") or {}
    out["meta"] = {k: (" ".join(str(v).split()) if isinstance(v, str) else v) for k, v in meta.items()}
    out["params"] = _param_entries(card, doc)
    out["datasets"] = _dataset_entries(card, doc)
    out["device_block"] = _device_block(name, doc)
    out["placeholders"] = card.placeholders() if card is not None else []
    out["overlays"] = {d["name"]: _overlay_status(name, d["name"]) for d in out["datasets"]}
    out["editable"] = bool(out["params"])
    out["shipped"] = not name.endswith(EDITED)
    return out


# ------------------------------------------------------------------ overlay builders
# Every curve below is the return value of the named fsim_core function; this section only
# arranges the arguments the rules file (physics-brief.md section 8) and the producing
# scripts (run_phase0/1/2/3.py, run_v1_validation.py, run_zhao_fit.py) use. Overlays are
# computed lazily on first request and cached for the process.
_CACHE: dict = {}
IRF_OE2008_PS = 500.0          # Reischle OE 16, 12771 (2008): 0.5 ns exponential IRF [V]
# Reischle's dip recovery time tau_d is not printed. It is the value the two published numbers per point
# imply (studio_support.implied_dip_recovery_time: raw 0.41 -> 0.15 and 0.43 -> 0.25 through the 0.5 ns
# exponential IRF give 1.13 and 1.58 ns [DR]); the overlay band spans that implied range, no invented one.
VC_RHO_WINDOW_CSV = OUT / "phase1" / "vc_rho_window.csv"
MAP_RHO_CSV = OUT / "phase3" / "map_rho_required_300K.csv"
ZHAO_COMPARISON_CSV = OUT / "zhao" / "zhao_fit_comparison.csv"
ZHAO_N_RUNS, ZHAO_T_END, ZHAO_DT, ZHAO_SEED = 60, 2.0e-9, 1.0e-13, 42   # run_zhao_fit.py
ZHAO_F_QUIET = 0.08            # [E] residual-suppression band 0.05-0.1


def _read_csv(path: Path) -> list:
    """A committed out/ CSV as a list of dicts (strings); the producing script owns the numbers."""
    try:
        with open(path, newline="", encoding="utf-8") as f:
            return list(csv.DictReader(f))
    except OSError as exc:
        raise LibraryError(500, f"cannot read {path.name}: {exc}")


def _f(v):
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return x if math.isfinite(x) else None


def _cached(key, fn):
    if key not in _CACHE:
        _CACHE[key] = fn()
    return _CACHE[key]


def _cls(card_name: str) -> dict:
    c = CLASSES[card_name]
    return {k: c[k] for k in ("code", "label", "reason")}


def _ov(card_name: str, tag: str, calls: list, label: str, *, series=(), bands=(), params_source=None,
        readouts=None) -> dict:
    out = {"card": card_name, "cls": _cls(card_name), "tag": tag, "calls": list(calls), "label": label,
           "series": list(series), "bands": list(bands), "params_source": params_source}
    if readouts:
        out["readouts"] = readouts
    return out


def _rows(doc: dict, dataset: str) -> list:
    return [dict(r) for r in ((doc.get("data") or {}).get(dataset) or {}).get("points") or []]


def _grid(rng, n):
    lo, hi = rng
    return np.linspace(lo, hi, n)


# ---- chatzarakis (V-a, class C): fitting.fit_phase0 and the integrator on the fit
def _vA_fit(card):
    from fsim_core.fitting import fit_phase0
    return _cached(("vA_fit", str(card.path)), lambda: fit_phase0(card))


def _vA_fit_shipped():
    from fsim_core.card import load_card
    return _vA_fit(load_card(CARDS / "chatzarakis.yaml"))


def _chatz_g2(card, doc):
    from fsim_core.fitting import V_A_TOL
    from fsim_core.integrator import g2_of_T, solve_Tc
    fit = _vA_fit(card)
    Td = np.array([r["T"] for r in fit.data])
    ws = np.array([r["w"] for r in fit.data])
    dxs = np.array([r["dx"] for r in fit.data])
    p = dict(fit.params)
    p["w"] = lambda T: np.interp(T, Td, ws)
    p["dx"] = lambda T: np.interp(T, Td, dxs)
    Ts = np.linspace(Td.min(), Td.max(), 121)
    g2 = [g2_of_T(T, p).g2 for T in Ts]
    tc_fit = solve_Tc(p)
    tc_fix = solve_Tc(dict(fit.params, w=3.0, dx=-0.7))
    label = (f"calibration, not held out (V_A_TOL {V_A_TOL}, fitting.py:44): the model is fitted to these points, "
             f"so it passes near them by construction. Published windows (w, dx) interpolated per point [E]; "
             f"the 78 K point is an upper bound. T_c = {tc_fit:.1f} K with the fit windows ({tc_fix:.1f} K at fixed "
             f"w = 3 meV, dx = -0.7 meV) lies beyond the data range and is not drawn.")
    return _ov("chatzarakis", fit.tag.name, ["fitting.fit_phase0", "integrator.g2_of_T", "integrator.solve_Tc"], label,
               series=[{"name": "F-series model", "x": Ts, "y": g2, "mode": "lines"}],
               params_source="joint V-a fit, seed 0 (reproduces out/phase0/fit_params.json)",
               readouts={"T_c_fit_windows_K": tc_fit, "T_c_fixed_w3_dx-0.7_K": tc_fix})


def _chatz_gamma(card, doc):
    from fsim_core.spectral import gamma_of_T
    fit = _vA_fit(card)
    P = fit.params
    Td = [r["T"] for r in _rows(doc, "gamma_vs_T")]
    Ts = np.linspace(min(Td), max(Td), 121)
    g = [float(gamma_of_T(T, P["gamma0"], P["a_ac"], P["b_lo"], P["E_lo"])) for T in Ts]
    gxx = [P["r_xx"] * v for v in g]
    label = ("calibration, not held out: Gamma(T) is fitted to these linewidths (15 % digitisation error). "
             "Gamma_XX = r_xx Gamma_X with r_xx fitted.")
    return _ov("chatzarakis", fit.tag.name, ["fitting.fit_phase0", "spectral.gamma_of_T"], label,
               series=[{"name": "Gamma_X (model)", "x": Ts, "y": g, "mode": "lines"},
                       {"name": "Gamma_XX = r_xx Gamma_X (model)", "x": Ts, "y": gxx, "mode": "lines"}],
               params_source="joint V-a fit, seed 0")


def _chatz_tau(card, doc):
    from fsim_core.studio_support import tau_anchored
    fit = _vA_fit(card)
    P = fit.params
    rows = sorted(_rows(doc, "tau_vs_T"), key=lambda r: r["T"])
    Ts = np.linspace(rows[0]["T"], rows[-1]["T"], 121)

    tau = [float(tau_anchored(T, rows[0]["T"], rows[0]["tau_ns"], P["a_esc"], P["E_a"], P["b_p"], P["E_b"]))
           for T in Ts]
    label = (f"calibration: the fit constrains the lifetime RATIO S(T)/S(T0), so the curve is anchored to the "
             f"first point ({rows[0]['T']:.0f} K) by construction; absolute tau is not predicted.")
    return _ov("chatzarakis", fit.tag.name, ["fitting.fit_phase0", "studio_support.tau_anchored"], label,
               series=[{"name": "tau0 S(T)/S(T0) (model)", "x": Ts, "y": tau, "mode": "lines"}],
               params_source="joint V-a fit, seed 0")


# ---- V1 device tier (device.evaluate), scripts/run_v1_validation.py
def _v1_design(lineshape, geom, alpha=None, rho_one=True):
    from fsim_core.device import DeviceDesign
    d = DeviceDesign()
    d.drive.V = 0.0
    d.drive.I_uA = 0.0
    d.drive.duty = 0.0
    d.drive.b_e = 0.0
    d.drive.mode = "PL"
    d.cavity.enabled = False
    d.filter.enabled = True
    d.filter.auto_w = False
    d.dot.lineshape = lineshape
    ph = dict(geom)
    if alpha is not None:
        ph["alpha_ps2"] = alpha
    d.dot.phonon = ph
    if rho_one:
        d.ret.b0 = 1e-12      # rho -> 1: the authors background-correct their correlations
        d.ret.beta = 1e-12
    return d


def _v1_point(d, T, w, dx, mu, delta, gamma_target=None):
    from fsim_core.device import class_proxy_params, evaluate
    from fsim_core.spectral import gamma_of_T
    d.thermal.T_hs = T
    d.filter.w = w
    d.filter.dx = dx
    d.drive.mu = mu
    d.dot.delta_xx = delta
    if gamma_target is not None:
        p = class_proxy_params()
        g_proxy = float(gamma_of_T(T, p["gamma0"], p["a_ac"], p["b_lo"], p["E_lo"]))
        d.dot.gamma_scale = gamma_target / g_proxy
    return float(evaluate(d, T_grid=[T])["scalars"]["g2_op"])


def _chatz2023_g2(card, doc):
    dev = doc["device"]
    pts = _rows(doc, "g2_vs_T")
    ibm, lor = [], []
    for pt in pts:
        res = {}
        for ls in ("lorentzian", "ibm"):
            d = _v1_design(ls, dev["geometry"], alpha=dev["phonon_alpha_ps2"], rho_one=False)
            res[ls] = _v1_point(d, float(pt["T"]), float(pt["w"]), float(pt["dx"]), dev["mu"], dev["delta_xx"])
        ibm.append(res["ibm"])
        lor.append(res["lorentzian"])
    Ts = [float(p["T"]) for p in pts]
    label = ("fitted [DR] parameters + published windows; IBM sideband correction. The Lorentzian line is a wiring "
             "check (it must recompose the validated V-a fit); the IBM line is the new content (sideband "
             "correction bounded by the standing V-a gate). Evaluated at the published points only.")
    return _ov("chatzarakis2023", "DR", ["device.evaluate"], label,
               series=[{"name": "IBM sideband correction", "x": Ts, "y": ibm, "mode": "lines+markers"},
                       {"name": "Lorentzian (wiring check)", "x": Ts, "y": lor, "mode": "lines+markers"}],
               params_source="cards/chatzarakis2023.yaml device: block")


def _laf2023_g2(card, doc):
    import itertools
    dev = doc["device"]
    pts = _rows(doc, "g2_vs_T")
    Ts, mid, lo, hi = [], [], [], []
    for pt in pts:
        T = float(pt["T"])
        key = str(int(T))
        gam, win = dev["gamma_meV"][key], dev["w_meV"][key]
        mu = dev["mu"]["mid_77"] if T < 100 else dev["mu"]["mid_high"]
        mu_rng = dev["mu"]["range_77"] if T < 100 else dev["mu"]["range_high"]
        d = _v1_design("ibm", dev["geometry"])
        mid.append(_v1_point(d, T, win["mid"], 0.0, mu, dev["delta_xx"]["mid"], gam["mid"]))
        corners = [_v1_point(d, T, ww, 0.0, mm, dd, gg)
                   for dd, gg, ww, mm in itertools.product(dev["delta_xx"]["range"], gam["range"],
                                                           win["range"], mu_rng)]
        lo.append(float(np.min(corners)))
        hi.append(float(np.max(corners)))
        Ts.append(T)
    label = ("envelope over the 16 corners of the unpublished Delta/Gamma/w/mu [E] ranges, with a mid-range line; "
             "the 4 K point is excluded (re-excitation, outside the spectral model). rho ~ 1 by the authors' "
             "background removal; IBM sideband lineshape. A band, never a fit.")
    return _ov("laferriere2023", "E", ["device.evaluate"], label,
               series=[{"name": "mid-range model", "x": Ts, "y": mid, "mode": "lines+markers"}],
               bands=[{"name": "16-corner envelope", "x": Ts, "lo": lo, "hi": hi}],
               params_source="cards/laferriere2023.yaml device: block")


# ---- laferriere (V-b, class M): the eps->1 envelope, scripts/run_phase2.py run_vb
def _laf_envelope(card, doc):
    from fsim_core.loading import f1b_g2
    from fsim_core.spectral import epsilon
    rows = sorted((r for r in _rows(doc, "g2_vs_T") if "mechanism" not in r), key=lambda r: r["T"])
    Ts = [r["T"] for r in rows]                    # the 4 K re-excitation point (it names a mechanism) is excluded
    deltas = _grid(card["delta_xx"].bounds, 9)
    lo, hi = [], []
    for T in Ts:
        g = []
        for d in deltas:
            if T < 100:
                g += [float(f1b_g2(mu, epsilon(d, gam, gam, w=w).eps))
                      for gam in _grid(card["gamma_77"].bounds, 5) for w in _grid(card["w_77"].bounds, 5)
                      for mu in _grid(card["mu_77"].bounds, 3)]
            elif T < 260:
                g += [float(f1b_g2(mu, epsilon(d, gam, gam, w=w).eps))
                      for gam in _grid(card["gamma_220"].bounds, 6) for w in _grid(card["w_220"].bounds, 5)
                      for mu in _grid(card["mu_high"].bounds, 6)]
            else:
                w300 = card["w_300"].fixed
                g += [float(f1b_g2(mu, epsilon(d, gam, gam, w=w300).eps))
                      for gam in _grid(card["gamma_300"].bounds, 6) for mu in _grid(card["mu_high"].bounds, 6)]
        lo.append(min(g))
        hi.append(max(g))
    label = ("envelope over the unpublished Delta/Gamma/w/mu [E] (eps -> 1 limit); the 4 K point is excluded "
             "(re-excitation). A band over swept ranges, not a fit and not a line.")
    return _ov("laferriere", "E", ["loading.f1b_g2", "spectral.epsilon"], label,
               bands=[{"name": "reachable envelope", "x": Ts, "lo": lo, "hi": hi}],
               params_source="card ranges (delta_xx, gamma_*, w_*, mu_*)")


# ---- reischle (V-c, class M)
def _reischle_rho(card, doc):
    # rho_required = sqrt(1 - g2) of the 100 MHz point is read from the committed
    # out/phase1/vc_rho_window.csv (scripts/run_phase1.py), never recomputed here.
    rows = card.datasets["g2_vs_ERR"].rows
    base = next(r for r in rows if r["device"] == 1 and r["position"] == 2 and r["ERR_MHz"] == 100.0)
    win = sorted(_read_csv(VC_RHO_WINDOW_CSV), key=lambda r: float(r["w_meV"]))
    ws = [float(r["w_meV"]) for r in win]
    req = [float(r["rho_required"]) for r in win]
    label = (f"trion: eps = 0 structurally, so g2 = 1 - rho^2 (integrator.g2_from at eps = 0). The line is the rho the "
             f"100 MHz point (g2 = {base['g2']}) requires, read from out/phase1/vc_rho_window.csv; the 200/500 MHz "
             f"points (0.53, 0.49) are temporal refilling, outside the spectral model.")
    return _ov("reischle", "DR", ["integrator.g2_from", "out/phase1/vc_rho_window.csv"], label,
               series=[{"name": "rho required by the 100 MHz point", "x": ws, "y": req, "mode": "lines"}],
               params_source="out/phase1/vc_rho_window.csv (run_phase1.py); cards/reischle.yaml g2_vs_ERR, device 1 position 2, 100 MHz")


def _reischle_oe(card, doc):
    from fsim_core.cw_g2 import intrinsic_from_raw
    from fsim_core.studio_support import implied_dip_recovery_time, inverted_background_g2
    rows = card.datasets["g2_oe2008"].rows
    Ts = [r["T"] for r in rows]
    g2s = [float(inverted_background_g2(r["g2_b"], r["rho"])) for r in rows]
    # tau_d implied by each published (raw, intrinsic) pair through the 0.5 ns exponential IRF [DR]
    tau_irf = IRF_OE2008_PS * 1e-3
    implied = [float(implied_dip_recovery_time(r["g2_raw"], r["g2_b"], tau_irf)) for r in rows]
    tds = np.linspace(min(implied), max(implied), 9)
    lo, hi = [], []
    for r in rows:
        v = [float(intrinsic_from_raw(r["g2_raw"], td, IRF_OE2008_PS, "exponential")) for td in tds]
        lo.append(min(v))
        hi.append(max(v))
    label = ("trion: eps = 0 structurally, g2 = 1 - rho^2. Markers: F2 inversion g2_s = (g2_b - (1 - rho^2))/rho^2 "
             "of the published g2_b with the published rho (matches g2_corr to 0.011). Band: cw_g2.intrinsic_from_raw "
             f"of the raw g2 through the {IRF_OE2008_PS:.0f} ps exponential IRF for dip recovery tau_d "
             f"{min(implied):.2f}-{max(implied):.2f} ns, the values the published raw/intrinsic pairs imply "
             f"(anchor {implied[0]:.2f} ns for raw 0.41 -> 0.15 [DR]; tau_d itself is not printed). Each point "
             "recovers its published g2_b at its own implied tau_d by construction.")
    return _ov("reischle", "DR", ["cw_g2.intrinsic_from_raw (cw_g2.deconvolve_dip)", "studio_support.inverted_background_g2",
                                  "studio_support.implied_dip_recovery_time"], label,
               series=[{"name": "g2_s (F2 inversion of g2_b)", "x": Ts, "y": g2s, "mode": "markers"}],
               bands=[{"name": "IRF deconvolution of g2_raw (implied tau_d range)", "x": Ts, "lo": lo, "hi": hi}],
               params_source="cards/reischle.yaml g2_oe2008; IRF 500 ps exponential [V]; tau_d implied [DR]",
               readouts={"tau_d_implied_ns": implied})


# ---- zhao (class C, verdict FAIL): sde.g2_vs_pump, scripts/run_zhao_fit.py settings
def _zhao_runs():
    from fsim_core.sde import QDLaserParams, find_threshold, g2_vs_pump

    def go():
        params = QDLaserParams()
        I_th = find_threshold(params)
        return {lab: g2_vs_pump(params, 4.0 * I_th, n_runs=ZHAO_N_RUNS, t_end=ZHAO_T_END, dt=ZHAO_DT,
                                seed=ZHAO_SEED, F_pump=F)
                for lab, F in (("normal", 1.0), ("quiet", ZHAO_F_QUIET))}
    return _cached("zhao_runs", go)


def _zhao_g2(card, doc):
    res = _zhao_runs()
    pub = {r["pump"]: r for r in card.datasets["g2_vs_pump"].rows}
    pumps = [r["pump"] for r in card.datasets["g2_vs_pump"].rows]
    y = [res[p]["g2_mean"] for p in pumps]
    se = [res[p]["g2_se"] for p in pumps]
    sig = {p: abs(res[p]["g2_mean"] - pub[p]["g2"]) / float(np.hypot(res[p]["g2_se"], pub[p]["err"])) for p in pumps}
    # CONTRACT rule 4: the verdict word is the committed CSV's `verdict`, never a re-threshold of this run's sigma
    cmp_rows = {r["pump"]: r for r in _read_csv(ZHAO_COMPARISON_CSV)}
    verdict = {p: cmp_rows[p]["verdict"] for p in pumps}
    sig_csv = {p: _f(cmp_rows[p]["sigma_ratio"]) for p in pumps}
    label = ("Langevin model of Zhao's QD LASER (thousands of dots, lasing mode), not a single-photon source; g2 ~ 1 is "
             f"correct physics. Verdict {verdict['normal']}/{verdict['quiet']} (normal/quiet), read from "
             f"out/zhao/zhao_fit_comparison.csv at {sig_csv['normal']:.1f} / {sig_csv['quiet']:.1f} sigma; this live run "
             f"gives {sig['normal']:.1f} / {sig['quiet']:.1f} sigma. Points only, with SE; no curve is drawn "
             f"through the quiet point. Axis is I/I_th (I/I_th = 4, {ZHAO_N_RUNS} runs, {ZHAO_T_END * 1e9:g} ns); "
             f"quiet = F_pump {ZHAO_F_QUIET} [E].")
    return _ov("zhao", "E", ["sde.QDLaserParams", "sde.find_threshold", "sde.g2_vs_pump"], label,
               series=[{"name": "sde.g2_vs_pump (SE)", "x": pumps, "y": y, "err": se, "mode": "markers"}],
               params_source="QDLaserParams() defaults, seed 42",
               readouts={"sigma": sig, "sigma_committed": sig_csv, "verdict": verdict,
                         "verdict_source": "out/zhao/zhao_fit_comparison.csv"})


# ---- qcap-piezo-variant (class P): rho requirement for T_c >= 300 K, scripts/run_phase3.py run_map
def _piezo_map(card, doc):
    # The requirement band rho >= sqrt(0.5/(1 - eps_min(300 K))) is read from the committed
    # out/phase3/map_rho_required_300K.csv (scripts/run_phase3.py run_map): NaN where rho > 1 is needed
    # (unreachable). The band is drawn where both Gamma(300) edges are reachable.
    from fsim_core.integrator import rho_of_T
    rows = _read_csv(MAP_RHO_CSV)
    g_lo, g_hi = card["gamma_300_lo"].fixed, card["gamma_300_hi"].fixed
    cols = [k for k in rows[0] if k.startswith("rho_req_gamma")]
    xs_b, lo, hi = [], [], []
    for r in rows:
        v = [_f(r[k]) for k in cols]
        if all(x is not None for x in v):
            xs_b.append(float(r["delta_meV"]))
            lo.append(min(v))
            hi.append(max(v))
    P = _vA_fit_shipped().params
    xs, ys = [], []
    for p in card.datasets["device_points"].rows:
        rho = p["rho"] if p["rho"] > 0 else float(rho_of_T(p["T"], P["a_esc"], P["E_a"], P["b_p"], P["E_b"],
                                                           P["b0"], P["beta"]))
        xs.append(p["delta"])
        ys.append(rho)
    label = ("requirement envelope: rho >= sqrt(0.5/(1 - eps_min(300 K))) over Gamma(300 K) 6-7 meV [V] "
             "(narrow-filter bound), read from out/phase3/map_rho_required_300K.csv. Chatzarakis 230 K rho computed "
             "from the V-a fit at runtime; Laferriere rho ~ 1 for contrast; the Gamma(300) rescale is [A]. The band is "
             "drawn where both Gamma(300) edges are reachable (rho <= 1).")
    return _ov("qcap-piezo-variant", "A", ["out/phase3/map_rho_required_300K.csv", "integrator.rho_of_T", "fitting.fit_phase0"],
               label,
               series=[{"name": "device points (rho from the card; Chatzarakis from the fit)", "x": xs, "y": ys,
                        "mode": "markers"}],
               bands=[{"name": "rho required for T_c >= 300 K", "x": xs_b, "lo": lo, "hi": hi}],
               params_source="out/phase3/map_rho_required_300K.csv (Gamma(300) %g / %g meV [V]); V-a fit rho(230 K)" % (g_lo, g_hi))


OVERLAYS.update({
    ("chatzarakis", "g2_vs_T"): _chatz_g2,
    ("chatzarakis", "gamma_vs_T"): _chatz_gamma,
    ("chatzarakis", "tau_vs_T"): _chatz_tau,
    ("chatzarakis2023", "g2_vs_T"): _chatz2023_g2,
    ("laferriere", "g2_vs_T"): _laf_envelope,
    ("laferriere2023", "g2_vs_T"): _laf2023_g2,
    ("reischle", "rho_spectral"): _reischle_rho,
    ("reischle", "g2_oe2008"): _reischle_oe,
    ("zhao", "g2_vs_pump"): _zhao_g2,
    ("qcap-piezo-variant", "device_points"): _piezo_map,
})

NO_OVERLAY.update({
    ("reischle", "g2_vs_ERR"): "no Gamma(T) fit exists for this device, so no g2_of_T overlay is allowed; the 200/500 MHz "
                               "points (0.53, 0.49) are temporal refilling, outside the spectral model",
    "kitamura": "anchor only; formula card Gamma_eff = Gamma + dg_inj (I/I_ref)^p_inj (EL/PL ratio 4 [V] feeds "
                "loading.gamma_eff as dg_inj); the card carries no dataset",
    "qcap-cavity": "prediction card, all inputs [A]: no data to overlay (requirement envelope awaiting MEEP/COMSOL)",
    "qcap-staged": "prediction card: no data to overlay (thermal maps over epi_k 8..40 are a Designer output, "
                   "not a Library overlay)",
})


# ------------------------------------------------------------------ overlays
def _overlay_status(name: str, dataset: str) -> dict:
    key = (base_of(name), dataset)
    if name.endswith(EDITED):
        return {"available": False, "reason": "overlays are drawn for shipped cards only (an edited copy is not the card the rules file reviewed)"}
    if key in OVERLAYS:
        return {"available": True}
    reason = NO_OVERLAY.get(key) or NO_OVERLAY.get(base_of(name))
    return {"available": False, "reason": reason or "no model overlay is allowed for this dataset by the physics rules file"}


def overlay(name: str, dataset: str) -> dict:
    card, doc, _ = _load(name)
    names = [d for d in (doc.get("data") or {})]
    if dataset not in names:
        raise LibraryError(404, f"card {name!r} has no dataset {dataset!r}")
    st = _overlay_status(name, dataset)
    if not st["available"]:
        raise LibraryError(404, st["reason"], reason=st["reason"])
    return OVERLAYS[(base_of(name), dataset)](card, doc)


# ------------------------------------------------------------------ save (Edit copy)
def _finite(x, what):
    try:
        f = float(x)
    except (TypeError, ValueError):
        raise LibraryError(400, f"{what}: not a number ({x!r})")
    if not math.isfinite(f):
        raise LibraryError(400, f"{what}: not finite")
    return f


def save_params(name: str, body: dict) -> dict:
    """Edit copy: the legacy fsim_gui/app.py card-editor schema. Each row is
    {value, lo, hi}: a non-null value replaces the range; otherwise lo and hi replace the
    value with a range. Only -edited names are written; a shipped name is a 409."""
    path = _path(name)
    if not name.endswith(EDITED):
        raise LibraryError(409, f"refusing to overwrite shipped card {name!r}; "
                                f"save under {name + EDITED!r}")
    from . import api_cards
    if (Path(api_cards.USER_CARDS) / f"{name}.yaml").is_file():
        # _find_card lets a cards/ file win over cards/studio/, so this save would hide a Studio card
        raise LibraryError(409, f"refusing to write {name!r}: a design saved by FSIM Studio already has that name "
                                f"(cards/studio/{name}.yaml); choose another name")
    base = base_of(name)
    base_path = _path(base)
    if not base_path.is_file():
        raise LibraryError(404, f"no shipped card {base!r} to copy")
    rows = body.get("params")
    if not isinstance(rows, dict) or not rows:
        raise LibraryError(400, "body must carry params: {name: {value} | {lo, hi}}")
    src_doc = _raw(path) if path.is_file() else _raw(base_path)
    if not _is_param_doc(src_doc) or not isinstance(src_doc.get("params"), dict) or not src_doc["params"]:
        raise LibraryError(422, f"card {base!r} has no editable params block")
    out = dict(src_doc)
    out["params"] = {k: dict(v) for k, v in src_doc["params"].items()}
    from fsim_core.card import _parse_param
    for n, row in rows.items():
        if n not in out["params"]:
            raise LibraryError(400, f"unknown parameter {n!r} (cards keep the shipped parameter set)")
        if not isinstance(row, dict):
            raise LibraryError(400, f"{n}: row must be an object")
        entry = dict(out["params"][n])
        if row.get("value") is not None:
            entry["value"] = _finite(row["value"], f"{n}.value")
            entry.pop("range", None)
        elif row.get("lo") is not None and row.get("hi") is not None:
            lo, hi = _finite(row["lo"], f"{n}.lo"), _finite(row["hi"], f"{n}.hi")
            if lo > hi:
                raise LibraryError(400, f"{n}: range must be (lo, hi), got [{lo}, {hi}]")
            entry["range"] = [lo, hi]
            entry.pop("value", None)
        else:
            raise LibraryError(400, f"{n}: give a value, or both lo and hi")
        try:
            _parse_param(n, entry)
        except (ValueError, KeyError) as exc:
            raise LibraryError(400, f"{n}: {exc}")
        out["params"][n] = entry
    path.write_text(yaml.safe_dump(out, sort_keys=False, allow_unicode=True), encoding="utf-8")
    try:
        rel = str(path.relative_to(CARDS.parent)).replace("\\", "/")
    except ValueError:
        rel = str(path)
    return {"name": name, "base": base, "path": rel}


# ------------------------------------------------------------------ routes
def register(app):
    def _out(obj, status=200):
        return Response(dumps(obj), status=status, mimetype="application/json")

    @app.errorhandler(LibraryError)
    def _library_error(exc):
        return _out({"error": exc.message, **exc.extra}, exc.status)

    @app.get("/api/params")
    def params_list():
        return _out(list_params())

    @app.get("/api/params/<name>")
    def params_get(name):
        return _out(get_params(name))

    @app.get("/api/params/<name>/overlay/<dataset>")
    def params_overlay(name, dataset):
        return _out(overlay(name, dataset))

    @app.put("/api/params/<name>")
    def params_put(name):
        body = request.get_json(silent=True)
        if not isinstance(body, dict):
            return _out({"error": "request body must be a JSON object"}, 400)
        return _out(save_params(name, body))

    return app


__all__ = ["register", "list_params", "get_params", "overlay", "save_params", "LibraryError",
           "OVERLAYS", "CLASSES", "PAPERS", "DERIVED_DESIGN"]
