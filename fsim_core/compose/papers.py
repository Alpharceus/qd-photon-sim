"""Paper-device library: loader, schema/class validation and the reproduce runner (design 3).

A paper card (cards/papers/<id>.yaml) is a composed device spec built ONLY from what one paper reports
(`device`), the numbers the paper reports (`reported`, `series`) and a ledger (`reproduce`) that says,
per reported item, what class of evidence the comparison is and what outcome is expected. The model is
never tuned to the paper; a mismatch is a recorded result, not a bug to hide.

Classes (design 3.1): C calibration (fitted, listed in `fits[].fitted_to`), M unfitted comparison
(the only class that counts as a test), T transcription (observable is a passthrough of an input that
was itself taken from this paper), P prediction (no reported value; the model number is just stated).

Outcome of one item: pass | fail | unavailable | error (error never matches a ledger).
  scalar pass iff |model - paper| <= max(k_sigma*err, abs, rel*|paper|) OR (ranged inputs) paper lies in
  the model [min, max] envelope; series pass iff rms_rel over the points <= tol.max.
  unavailable iff a chosen drive is `needs_physics`, the solver is unavailable, the central run is flagged
  invalid (e.g. strong coupling, `valid` False) or the model value is not finite.
Raw vs corrected: the reported convention is compared with the SAME model convention only.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from pathlib import Path
from typing import Optional

import numpy as np
import yaml

from . import engine
from . import registry as _reg
from .params import TAGS
from .solvers import SolverUnavailable
from .spec import AXIS_KEYS, SCHEMA, ComposedSpec, validate_spec

READ_LEVELS = ("V-full", "V-abs", "V-snip")
CLASSES = ("C", "M", "T", "P")
EXPECTS = ("pass", "fail", "unavailable")
REPORTED_TAGS = ("V", "DR")
SCALAR_TOL_KEYS = ("k_sigma", "abs", "rel")
DEFAULT_DIR = "cards/papers"


class PaperSchemaError(ValueError):
    """Raised with `.problems`, the full list of schema / class-rule violations."""

    def __init__(self, problems, paper_id=""):
        self.problems = list(problems)
        self.paper_id = paper_id
        head = f"paper card {paper_id!r}: " if paper_id else "paper card: "
        super().__init__(head + f"{len(self.problems)} problem(s):\n  - " + "\n  - ".join(self.problems))


# ------------------------------------------------------------------ data classes
@dataclass
class Paper:
    meta: dict
    device: dict                       # composed spec dict (schema key added)
    reported: list
    series: list
    reproduce: list
    fits: list
    inputs_from_paper: set             # leaf paths tagged V/DR with `where` in this card
    raw: dict
    needs_physics: list = field(default_factory=list)   # ["drive:arp", ...] chosen kinds awaiting physics
    path: Optional[str] = None

    @property
    def id(self) -> str:
        return self.meta["id"]


@dataclass
class ItemResult:
    id: str
    cls: str
    kind: str                          # scalar | series | predict
    observable: str
    convention: str
    paper: Optional[float] = None
    err: Optional[float] = None
    where: str = ""
    model: float = float("nan")        # series: rms_rel
    lo: float = float("nan")
    hi: float = float("nan")
    tol: dict = field(default_factory=dict)
    tol_str: str = ""
    outcome: str = "error"             # pass | fail | unavailable | error
    expect: str = "pass"
    match: bool = False
    inputs: list = field(default_factory=list)   # device leaves (from this paper) the observable depends on
    note: str = ""
    tol_reason: str = ""
    reason: str = ""

    @property
    def request(self) -> str:
        return self.observable if self.convention == "intrinsic" else f"{self.observable}:{self.convention}"


@dataclass
class ReproduceReport:
    paper_id: str
    items: list
    valid: bool = True
    invalid_reasons: list = field(default_factory=list)
    solver: str = ""
    notes: list = field(default_factory=list)

    @property
    def matched(self) -> bool:
        return all(i.match for i in self.items)

    def class_summary(self) -> dict:
        out = {c: {"n": 0, "pass": 0, "fail": 0, "unavailable": 0, "error": 0} for c in CLASSES}
        for i in self.items:
            out[i.cls]["n"] += 1
            out[i.cls][i.outcome] += 1
        return out


# ------------------------------------------------------------------ schema
def _is_num(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def _nonempty(x) -> bool:
    return isinstance(x, str) and bool(x.strip())


def _device_leaf_paths(device: dict) -> dict:
    """{leaf path: raw leaf} for T_K and every axis param actually written in the device."""
    out = {}
    if "T_K" in device:
        out["T_K"] = device["T_K"]
    for a in AXIS_KEYS:
        blk = device.get(a)
        if isinstance(blk, dict) and isinstance(blk.get("params"), dict):
            for k, v in blk["params"].items():
                out[f"{a}.{k}"] = v
    return out


def _observable(oid):
    try:
        return _reg.get("observable", oid)
    except Exception:
        return None


def _validate(raw: dict, stem: Optional[str]) -> tuple:
    """-> (problems, inputs_from_paper, needs_physics, device_spec_dict)."""
    P: list = []
    if not isinstance(raw, dict):
        return ["top level must be a mapping"], set(), [], {}
    _reg.ensure_loaded()
    meta = raw.get("meta")
    if not isinstance(meta, dict):
        P.append("meta: missing or not a mapping")
        meta = {}
    if not _nonempty(meta.get("id")):
        P.append("meta.id: required")
    elif stem is not None and meta["id"] != stem:
        P.append(f"meta.id {meta['id']!r} != file stem {stem!r}")
    if meta.get("role") != "paper-device":
        P.append(f"meta.role must be 'paper-device', got {meta.get('role')!r}")
    if not _nonempty(meta.get("cite")):
        P.append("meta.cite: required (full citation)")
    if meta.get("read_level") not in READ_LEVELS:
        P.append(f"meta.read_level must be one of {list(READ_LEVELS)}, got {meta.get('read_level')!r}")

    # ---- device
    device = raw.get("device")
    spec_dict: dict = {}
    inputs: set = set()
    needs: list = []
    if not isinstance(device, dict):
        P.append("device: missing or not a mapping")
        device = {}
    else:
        spec_dict = {"schema": SCHEMA, **device}
        leaves = _device_leaf_paths(device)
        for path, leaf in leaves.items():
            if not isinstance(leaf, dict):
                P.append(f"device.{path}: leaf must be a mapping with tag and where|source (bare scalars "
                         "carry no provenance)")
                continue
            if leaf.get("tag") not in TAGS:
                P.append(f"device.{path}: tag must be one of {list(TAGS)}, got {leaf.get('tag')!r}")
            has_where, has_src = _nonempty(leaf.get("where")), _nonempty(leaf.get("source"))
            if not (has_where or has_src):
                P.append(f"device.{path}: needs where (paper location) or source")
            if leaf.get("tag") in REPORTED_TAGS and not has_where:
                P.append(f"device.{path}: tag {leaf.get('tag')} (from the paper) needs a non-empty where")
            if leaf.get("tag") in REPORTED_TAGS and has_where:
                inputs.add(path)
        try:
            spec = ComposedSpec.from_dict(spec_dict)
        except ValueError as e:
            P.append(f"device: {e}")
            spec = None
        if spec is not None:
            for axis in AXIS_KEYS:
                try:
                    ent = _reg.get(axis, getattr(spec, axis)["kind"])
                    if ent.status == "needs_physics":
                        needs.append(f"{axis}:{ent.kind}")
                except Exception:
                    pass
            res = validate_spec(spec)
            for i in res["issues"]:
                if i["level"] != "error":
                    continue
                if "awaiting physics" in i["message"]:
                    continue           # recorded in `needs_physics`; items become unavailable
                P.append(f"device {i['path']}: {i['message']}")

    # ---- reported / series
    ids: dict = {}

    def claim(kind, i, ent):
        if not isinstance(ent, dict):
            P.append(f"{kind}[{i}]: must be a mapping")
            return None
        eid = ent.get("id")
        if not _nonempty(eid):
            P.append(f"{kind}[{i}]: id required")
            return None
        if eid in ids:
            P.append(f"{kind}[{i}]: duplicate id {eid!r} (ids are unique across reported and series)")
        ids[eid] = (kind, ent)
        return eid

    for i, r in enumerate(_aslist(raw.get("reported"), "reported", P)):
        eid = claim("reported", i, r)
        if eid is None:
            continue
        w = f"reported[{eid}]"
        o = _observable(r.get("observable"))
        if o is None:
            P.append(f"{w}: observable {r.get('observable')!r} is not a registered observable id")
        else:
            conv = r.get("convention", "intrinsic")
            if conv not in o.conventions:
                P.append(f"{w}: convention {conv!r} invalid for {o.id}; valid: {list(o.conventions)}")
        if not _is_num(r.get("value")):
            P.append(f"{w}: value must be a number")
        if r.get("err") is not None and not (_is_num(r["err"]) and r["err"] >= 0):
            P.append(f"{w}: err must be a non-negative number")
        if r.get("tag") not in REPORTED_TAGS:
            P.append(f"{w}: tag must be V or DR, got {r.get('tag')!r}")
        if not _nonempty(r.get("where")):
            P.append(f"{w}: where required")

    series_ids = set()
    resolved = {}
    if spec_dict:
        try:
            resolved = validate_spec(spec_dict)["resolved"]
        except Exception:
            resolved = {}
    for i, s in enumerate(_aslist(raw.get("series"), "series", P)):
        eid = claim("series", i, s)
        if eid is None:
            continue
        series_ids.add(eid)
        w = f"series[{eid}]"
        x = s.get("x")
        if not (isinstance(x, dict) and _nonempty(x.get("path"))):
            P.append(f"{w}: x must be a mapping with a leaf path")
        else:
            if resolved and x["path"] not in resolved:
                P.append(f"{w}: x.path {x['path']!r} is not a device leaf path")
            elif resolved and x.get("unit") != resolved[x["path"]]["unit"]:
                P.append(f"{w}: x.unit {x.get('unit')!r} != leaf unit {resolved[x['path']]['unit']!r}")
        y = s.get("y")
        o = _observable(y.get("observable")) if isinstance(y, dict) else None
        if o is None:
            P.append(f"{w}: y.observable must be a registered observable id")
        elif (y.get("convention", "intrinsic")) not in o.conventions:
            P.append(f"{w}: y.convention {y.get('convention')!r} invalid for {o.id}")
        if s.get("tag") not in REPORTED_TAGS:
            P.append(f"{w}: tag must be V or DR, got {s.get('tag')!r}")
        if not _nonempty(s.get("where")):
            P.append(f"{w}: where required")
        if s.get("digitized") is True:
            if s.get("tag") != "DR":
                P.append(f"{w}: digitized series must be tag DR")
            if not _nonempty(s.get("method")):
                P.append(f"{w}: digitized series needs a non-empty `method` note")
        pts = s.get("points")
        if not (isinstance(pts, list) and pts and all(
                isinstance(p, (list, tuple)) and len(p) in (2, 3) and all(_is_num(v) for v in p) for p in pts)):
            P.append(f"{w}: points must be a non-empty list of [x, y] or [x, y, yerr] numbers")

    # ---- fits
    fits = _aslist(raw.get("fits"), "fits", P)
    fitted: set = set()
    for i, f in enumerate(fits):
        if not isinstance(f, dict) or not _nonempty(f.get("param")) or not isinstance(f.get("fitted_to"), list):
            P.append(f"fits[{i}]: needs param and a fitted_to list")
            continue
        if resolved and f["param"] not in resolved:
            P.append(f"fits[{i}]: param {f['param']!r} is not a device leaf path")
        for t in f["fitted_to"]:
            if t not in ids:
                P.append(f"fits[{i}]: fitted_to {t!r} is not a reported/series id")
            fitted.add(t)

    # ---- reproduce + class rules
    seen = set()
    for i, e in enumerate(_aslist(raw.get("reproduce"), "reproduce", P)):
        if not isinstance(e, dict):
            P.append(f"reproduce[{i}]: must be a mapping")
            continue
        item = e.get("item")
        w = f"reproduce[{item}]"
        if not _nonempty(item):
            P.append(f"reproduce[{i}]: item required")
            continue
        if item in seen:
            P.append(f"{w}: duplicate reproduce entry")
        seen.add(item)
        cls = e.get("class")
        if cls not in CLASSES:
            P.append(f"{w}: class must be one of {list(CLASSES)}, got {cls!r}")
        if e.get("expect") not in EXPECTS:
            P.append(f"{w}: expect must be one of {list(EXPECTS)}, got {e.get('expect')!r}")
        elif e["expect"] != "pass" and not _nonempty(e.get("reason")):
            P.append(f"{w}: reason required when expect is {e['expect']!r}")
        if item in fitted and cls != "C":
            P.append(f"{w}: appears in fits[].fitted_to so its class must be C, got {cls}")
        if cls == "P":
            if item in ids:
                P.append(f"{w}: class P means no reported value, but {item!r} is a reported/series id")
            o = _observable(e.get("observable"))
            if o is None:
                P.append(f"{w}: class P entries name a registered `observable`")
            elif e.get("convention", "intrinsic") not in o.conventions:
                P.append(f"{w}: convention {e.get('convention')!r} invalid for {o.id}")
            continue
        if item not in ids:
            P.append(f"{w}: item is not a reported or series id (only class P may be free-standing)")
            continue
        kind, ent = ids[item]
        oid = ent.get("observable") if kind == "reported" else (ent.get("y") or {}).get("observable")
        o = _observable(oid)
        pt = getattr(o, "passthrough_of", None) if o is not None else None
        pt_from_paper = bool(pt) and pt in inputs
        if pt_from_paper and cls != "T":
            P.append(f"{w}: observable {oid} is a passthrough of {pt}, an input taken from this paper, "
                     f"so the class must be T, got {cls}")
        if cls == "T" and not pt_from_paper:
            P.append(f"{w}: class T requires an observable that passes through a paper-sourced input "
                     f"leaf (passthrough_of={pt!r})")
        if cls == "M" and item in fitted:
            P.append(f"{w}: class M item is fitted")
        # tolerance
        tol = e.get("tol")
        if not isinstance(tol, dict) or not tol:
            P.append(f"{w}: tol required")
        elif kind == "series":
            if tol.get("metric") != "rms_rel" or not (_is_num(tol.get("max")) and tol["max"] > 0):
                P.append(f"{w}: series tol must be {{metric: rms_rel, max: <positive number>}}")
        else:
            bad = set(tol) - set(SCALAR_TOL_KEYS)
            if bad:
                P.append(f"{w}: unknown tol keys {sorted(bad)}; allowed {list(SCALAR_TOL_KEYS)}")
            if not (set(tol) & set(SCALAR_TOL_KEYS)):
                P.append(f"{w}: tol needs at least one of {list(SCALAR_TOL_KEYS)}")
            for k in SCALAR_TOL_KEYS:
                if k in tol and not (_is_num(tol[k]) and tol[k] >= 0):
                    P.append(f"{w}: tol.{k} must be a non-negative number")
            if "k_sigma" in tol and ent.get("err") is None:
                P.append(f"{w}: tol.k_sigma needs the reported item to carry err")
        if not _nonempty(e.get("tol_reason")):
            P.append(f"{w}: tol_reason required (why this tolerance, written before the run)")

    hist = raw.get("history")
    if hist is not None:
        if not isinstance(hist, list) or not all(
                isinstance(h, dict) and h.get("date") and _nonempty(h.get("reason")) for h in hist):
            P.append("history: must be a list of {date, reason, ...}")
    return P, inputs, needs, spec_dict


def _aslist(v, name, P) -> list:
    if v is None:
        return []
    if not isinstance(v, list):
        P.append(f"{name}: must be a list")
        return []
    return v


def paper_from_dict(raw: dict, stem: Optional[str] = None, path: Optional[str] = None) -> Paper:
    problems, inputs, needs, spec_dict = _validate(raw, stem)
    pid = (raw.get("meta") or {}).get("id", stem or "") if isinstance(raw, dict) else (stem or "")
    if problems:
        raise PaperSchemaError(problems, pid)
    return Paper(meta=raw["meta"], device=spec_dict, reported=list(raw.get("reported") or []),
                 series=list(raw.get("series") or []), reproduce=list(raw["reproduce"]) if raw.get("reproduce") else [],
                 fits=list(raw.get("fits") or []), inputs_from_paper=inputs, raw=raw,
                 needs_physics=needs, path=path)


def load_paper(path) -> Paper:
    p = Path(path)
    with open(p, encoding="utf-8") as fh:
        raw = yaml.safe_load(fh)
    return paper_from_dict(raw, stem=p.stem, path=str(p))


def list_papers(dir=DEFAULT_DIR) -> list:
    d = Path(dir)
    if not d.is_dir():
        return []
    return [load_paper(f) for f in sorted(d.glob("*.yaml")) if not f.name.startswith("_")]


# ------------------------------------------------------------------ reproduce
def _tol_str(tol: dict, kind: str) -> str:
    if kind == "series":
        return f"rms_rel<={tol['max']:g}"
    parts = []
    if "k_sigma" in tol:
        parts.append(f"{tol['k_sigma']:g}sigma")
    if "abs" in tol:
        parts.append(f"abs{tol['abs']:g}")
    if "rel" in tol:
        parts.append(f"rel{tol['rel'] * 100:g}%")
    return "|".join(parts)


def _scalar_outcome(paper, err, model, lo, hi, tol) -> str:
    allow = max((tol.get("k_sigma", 0.0) * err) if err is not None else 0.0,
                tol.get("abs", 0.0), tol.get("rel", 0.0) * abs(paper))
    if abs(model - paper) <= allow:
        return "pass"
    if hi > lo and lo <= paper <= hi:
        return "pass"
    return "fail"


def _rms_rel(model, paper) -> float:
    r = [((m - p) / p) ** 2 for m, p in zip(model, paper) if p != 0.0]
    return math.sqrt(sum(r) / len(r)) if r else float("nan")


def reproduce(paper: Paper, *, envelope: bool = True) -> ReproduceReport:
    """Run the engine at the paper's device and judge every ledger entry (module doc)."""
    ids = {r["id"]: ("scalar", r) for r in paper.reported}
    ids.update({s["id"]: ("series", s) for s in paper.series})
    items: list = []

    def base(e) -> ItemResult:
        item = e["item"]
        if e["class"] == "P":
            kind, ent = "predict", {}
            oid, conv = e["observable"], e.get("convention", "intrinsic")
        else:
            kind, ent = ids[item]
            if kind == "scalar":
                oid, conv = ent["observable"], ent.get("convention", "intrinsic")
            else:
                oid, conv = ent["y"]["observable"], ent["y"].get("convention", "intrinsic")
        r = ItemResult(id=item, cls=e["class"], kind=kind, observable=oid, convention=conv,
                       expect=e["expect"], where=ent.get("where", ""), tol=dict(e.get("tol") or {}),
                       tol_reason=e.get("tol_reason", ""), reason=e.get("reason", ""))
        if kind == "scalar":
            r.paper, r.err = float(ent["value"]), (None if ent.get("err") is None else float(ent["err"]))
        if r.tol:
            r.tol_str = _tol_str(r.tol, kind)
        return r

    for e in paper.reproduce:
        items.append(base(e))

    def finish(note_all: str, outcome: str):
        for r in items:
            r.outcome, r.note = outcome, note_all
            r.match = r.outcome == r.expect
        return ReproduceReport(paper.id, items, valid=False, invalid_reasons=[note_all],
                               notes=list(paper.needs_physics))

    if paper.needs_physics:
        return finish("awaiting physics: " + ", ".join(paper.needs_physics), "unavailable")

    spec = paper.device
    try:
        central = engine.evaluate_spec(spec)["scalars"]
    except SolverUnavailable as e:
        return finish(f"solver unavailable: {e}", "unavailable")
    except Exception as e:                                       # noqa: BLE001
        return finish(f"engine error: {type(e).__name__}: {e}", "error")
    rep_valid, rep_reasons = bool(central["valid"]), list(central["invalid_reasons"])
    report = ReproduceReport(paper.id, items, valid=rep_valid, invalid_reasons=rep_reasons,
                             solver=str(central.get("solver", "")), notes=list(central.get("notes", [])))

    scal = [r for r in items if r.kind in ("scalar", "predict")]
    env: dict = {}
    if scal:
        reqs = sorted({r.request for r in scal})
        try:
            env = engine.envelope(spec, reqs) if envelope else engine.envelope(spec, reqs, max_corners=0, n_lhs=0)
        except SolverUnavailable as e:
            for r in scal:
                r.outcome, r.note = "unavailable", f"solver unavailable: {e}"
        except Exception as e:                                   # noqa: BLE001
            for r in scal:
                r.outcome, r.note = "error", f"engine error: {type(e).__name__}: {e}"
    for r in scal:
        if r.outcome in ("unavailable", "error") and r.note:
            continue
        ev = env.get(r.request)
        if ev is None:
            r.outcome, r.note = "error", "no envelope result"
            continue
        r.model, (r.lo, r.hi) = float(ev["central"]), ev["range"]
        r.inputs = [p for p in ev["inputs"] if p in paper.inputs_from_paper]
        if not rep_valid:
            r.outcome, r.note = "unavailable", "model run flagged invalid: " + "; ".join(rep_reasons)
        elif not math.isfinite(r.model):
            r.outcome = "unavailable"
            r.note = scalars_reason(central, r.observable) or "model value not finite"
        elif r.kind == "predict":
            r.outcome, r.note = "pass", "prediction (no reported value); model run valid"
        else:
            r.outcome = _scalar_outcome(r.paper, r.err, r.model, r.lo, r.hi, r.tol)

    for r in items:
        if r.kind != "series":
            continue
        s = ids[r.id][1]
        pts = s["points"]
        xs = [float(p[0]) for p in pts]
        py = [float(p[1]) for p in pts]
        try:
            sr = engine.series(spec, s["x"]["path"], xs, [r.request])
        except SolverUnavailable as e:
            r.outcome, r.note = "unavailable", f"solver unavailable: {e}"
            continue
        except Exception as e:                                   # noqa: BLE001
            r.outcome, r.note = "error", f"engine error: {type(e).__name__}: {e}"
            continue
        my = [float(v) for v in sr["y"][r.request]]
        if sr["failed"] or not all(math.isfinite(v) for v in my):
            r.outcome, r.note = "unavailable", "model undefined at some series points"
        elif not rep_valid:
            r.outcome, r.note = "unavailable", "model run flagged invalid: " + "; ".join(rep_reasons)
        else:
            r.model = _rms_rel(my, py)
            r.lo, r.hi = min(my), max(my)
            r.outcome = "pass" if r.model <= r.tol["max"] else "fail"
        try:
            ev = engine.envelope(spec, [r.request], max_corners=0, n_lhs=0)[r.request]
            r.inputs = [p for p in ev["inputs"] if p in paper.inputs_from_paper]
        except Exception:                                        # noqa: BLE001
            pass
    for r in items:
        r.match = r.outcome == r.expect
    return report


def scalars_reason(scalars: dict, oid: str) -> str:
    return str(scalars.get(f"{oid}_invalid_reason", "") or "")
