"""Composer and Papers workspaces: HTTP side of fsim_core.compose (devlib-p10, design 4.1).

Three-layer rule: nothing here is physics. Every number is the return value of a public
fsim_core.compose call (registry_payload, validate_spec, engine.estimate / series, papers.reproduce,
tiers.tiers_payload) or of the existing device.evaluate path through jobs.JobManager.submit_run.

Routes (registered by register(app) from server.create_app):
  GET  /api/compose/registry          registry_payload() + tiers_payload(), cached per fsim_core hash
  POST /api/compose/validate          {spec} -> {ok, solver, regime, issues, resolved, eta_s}
  POST /api/compose/run               {spec, mode: point|tsweep, T_grid?, name?} -> same body as POST /api/run
  POST /api/compose/series            {spec, path, values[], observables[]} -> job (series_job), cached
  GET  /api/compose/examples          class-example composed designs from cards/composed/*.yaml
  GET  /api/compose/tier/fock         fsim_core.compose.tiers.fock_card_to_spec() as a composed design
  GET|POST /api/compose/scene         fsim_studio.scene_compose.describe(spec) (501 while the module is absent)
  GET  /api/papers                    catalogue (fsim_core.compose.papers.list_papers)
  GET  /api/papers/<id>               card JSON + resolved device
  POST /api/papers/<id>/reproduce     job (reproduce_job), cache key = card sha256 + fsim_core hash
  GET  /api/papers/<id>/design        composed design at the central values ("Open in Composer")

Job functions (exact entries in jobs.ALLOWED_CALLS):
  fsim_studio.api_compose:series_job      fsim_studio.api_compose:reproduce_job
"""
from __future__ import annotations

import copy
import hashlib
import importlib
import json
import math
import re

from . import CARDS, fsim_core_hash
from . import cache

MAX_SERIES_POINTS = 200
MAX_TSWEEP_POINTS = 31
MAX_OBSERVABLES = 16
PAPER_ID_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9_-]{0,80}$")
PAPERS_DIR = CARDS / "papers"
EXAMPLES_DIR = CARDS / "composed"
SERIES_FN = "fsim_studio.api_compose:series_job"
REPRODUCE_FN = "fsim_studio.api_compose:reproduce_job"
CENTRAL_NOTE = "ranged leaves replaced by their declared central value (lin midpoint, log geometric mean); not a prediction"

_REGISTRY_CACHE: dict = {}


class ComposeError(Exception):
    def __init__(self, message, status=400, **extra):
        super().__init__(message)
        self.message = message
        self.status = status
        self.extra = extra


# ------------------------------------------------------------------ optional modules (P6, P7)
def _optional(name: str):
    """fsim_core.compose.<name> or None while that module is not written yet."""
    try:
        return importlib.import_module(f"fsim_core.compose.{name}")
    except ModuleNotFoundError as exc:
        if exc.name == f"fsim_core.compose.{name}":
            return None
        raise


def _need_papers():
    papers = _optional("papers")
    if papers is None:
        raise ComposeError("fsim_core.compose.papers is not available yet", 501)
    return papers


# ------------------------------------------------------------------ pure helpers
def registry_payload_cached() -> dict:
    """registry_payload() with the tier list merged in (tiers_payload() when tiers.py exists)."""
    h = fsim_core_hash()
    hit = _REGISTRY_CACHE.get(h)
    if hit is None:
        from fsim_core.compose import registry
        registry.ensure_loaded()
        payload = registry.registry_payload()
        tiers = _optional("tiers")
        if tiers is not None and hasattr(tiers, "tiers_payload"):
            payload = dict(payload, tiers=tiers.tiers_payload())
        _REGISTRY_CACHE.clear()
        _REGISTRY_CACHE[h] = hit = payload
    return hit


def validate_payload(spec) -> dict:
    """validate_spec(spec) plus eta_s from engine.estimate (null when the spec cannot run)."""
    from fsim_core.compose import engine, validate_spec
    res = validate_spec(spec)
    out = dict(res)
    out["issues"] = list(res["issues"])
    out["eta_s"] = None
    if res["ok"]:
        try:
            out["eta_s"] = float(engine.estimate(spec, 1))
        except Exception as exc:  # noqa: BLE001 -- SolverUnavailable etc.: reported as an issue
            out["ok"] = False
            out["issues"].append({"level": "error", "path": "solver", "message": f"{type(exc).__name__}: {exc}"})
    return out


def _check_observables(observables) -> list:
    from fsim_core.compose import observables as obs_mod
    if not isinstance(observables, list) or not observables:
        raise ComposeError("observables must be a non-empty list of observable requests ('id' or 'id:convention')")
    if len(observables) > MAX_OBSERVABLES:
        raise ComposeError(f"at most {MAX_OBSERVABLES} observables per request")
    out = []
    for r in observables:
        if not isinstance(r, str):
            raise ComposeError(f"observable request {r!r} must be a string")
        try:
            oid, conv = obs_mod.parse_request(r)
            ob = obs_mod.get_observable(oid)
            convs = list(getattr(ob, "conventions", ()) or ())
            if ":" not in r and conv not in convs and len(convs) == 1:
                conv = convs[0]
            if convs and conv not in convs:
                raise ValueError(f"observable {oid!r} has no convention {conv!r}; known: {convs}; "
                                 f"ask for '{oid}:{convs[0]}'")
        except Exception as exc:  # noqa: BLE001
            raise ComposeError(f"bad observable request {r!r}: {exc}") from exc
        out.append(r)
    return out


def _valid_spec(body: dict):
    spec = body.get("spec")
    if not isinstance(spec, dict):
        raise ComposeError("body must carry spec: a composed spec object")
    res = validate_payload(spec)
    if not res["ok"]:
        raise ComposeError("invalid composed spec", 400, issues=res["issues"])
    return spec, res


def series_eta(spec, n: int) -> float:
    from fsim_core.compose import engine
    try:
        return float(engine.estimate(spec, n))
    except Exception:  # noqa: BLE001 -- advisory only
        return 0.0


def _card_path(paper_id: str):
    if not isinstance(paper_id, str) or not PAPER_ID_RE.match(paper_id):
        raise ComposeError(f"invalid paper id {paper_id!r}", 400)
    return PAPERS_DIR / f"{paper_id}.yaml"


def _card_sha(path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _load_paper(paper_id: str):
    papers = _need_papers()
    path = _card_path(paper_id)
    if paper_id.startswith("_") or not path.is_file():
        raise ComposeError(f"no paper {paper_id!r}", 404)
    try:
        return papers.load_paper(path), path
    except Exception as exc:  # noqa: BLE001 -- PaperSchemaError carries a problem list
        problems = getattr(exc, "problems", None)
        raise ComposeError(f"paper {paper_id!r} failed to load: {exc}", 422,
                           **({"problems": list(problems)} if problems else {})) from exc


def _pget(paper, name, default=None):
    if hasattr(paper, name):
        return getattr(paper, name)
    raw = getattr(paper, "raw", None) or {}
    return raw.get(name, default)


def _reproduce_kwargs(paper_id: str, path) -> dict:
    return {"paper_id": paper_id, "card_sha256": _card_sha(path), "papers_dir": str(path.parent)}


def _reproduce_key(paper_id: str, path) -> str:
    return cache.make_key({"call": REPRODUCE_FN, "kwargs": _reproduce_kwargs(paper_id, path)})


def _last_outcome(paper_id: str, path):
    """Summary of the cached reproduce report for this exact card + core hash, else None."""
    rec = cache.get(_reproduce_key(paper_id, path))
    if rec is None:
        return None
    result = rec.get("result") or {}
    items = result.get("items")
    if not isinstance(items, list):
        report = result.get("report")
        items = report.get("items") if isinstance(report, dict) else None
    if not isinstance(items, list):
        return {"cached": True, "n_items": None, "n_match": None, "all_match": None}
    flags = [it.get("match") for it in items if isinstance(it, dict) and "match" in it]
    return {"cached": True, "n_items": len(items), "n_match": sum(1 for f in flags if f is True),
            "all_match": bool(flags) and all(f is True for f in flags)}


def paper_summary(paper, path) -> dict:
    meta = _pget(paper, "meta", {}) or {}
    device = _pget(paper, "device", {}) or {}
    repro = _pget(paper, "reproduce", []) or []
    classes = {"C": 0, "M": 0, "T": 0, "P": 0, "R": 0}
    for it in repro:
        c = it.get("class") if isinstance(it, dict) else None
        if c in classes:
            classes[c] += 1
    kinds = {a: ((device.get(a) or {}).get("kind") if isinstance(device.get(a), dict) else None)
             for a in ("material", "structure", "drive")}
    return {"id": meta.get("id"), "cite": meta.get("cite"), "read_level": meta.get("read_level"),
            "summary": meta.get("summary"), "kinds": kinds, "classes": classes, "n_items": len(repro),
            "last_outcome": _last_outcome(meta.get("id") or path.stem, path)}


def _param_spec(axis, kind, name):
    """The ParamSpec of `name` on a registered axis entry (T_K for axis None), or None when unknown."""
    try:
        if axis is None:
            from fsim_core.compose.spec import T_K_SPEC
            return T_K_SPEC
        from fsim_core.compose import registry
        registry.ensure_loaded()
        return registry.REGISTRY.get(axis, kind).param(name)
    except Exception:  # noqa: BLE001 -- unknown kind/param: validation reports it, the collapse stays best-effort
        return None


def _collapse_leaf(leaf, axis=None, kind=None, name=None):
    """A ranged / choice leaf mapping -> a value leaf at its declared central value; others unchanged.
    Ranged leaves go through params.resolve_leaf + params.central, so a log-scale ParamSpec gives the geometric mean
    exactly as the engine does (R2 finding 8)."""
    if not isinstance(leaf, dict):
        return leaf, False
    if "range" in leaf:
        pspec = _param_spec(axis, kind, name)
        c = None
        if pspec is not None:
            try:
                from fsim_core.compose import params as _params
                c = _params.central(_params.resolve_leaf(leaf, pspec, name or "T_K"))
            except ValueError:
                c = None
        if c is None:
            lo, hi = (float(x) for x in leaf["range"])
            c = math.sqrt(lo * hi) if leaf.get("scale") == "log" and lo > 0 and hi > 0 else 0.5 * (lo + hi)
        out = {k: v for k, v in leaf.items() if k not in ("range", "scale")}
        out["value"] = c
        return out, True
    if "choice" in leaf:
        out = {k: v for k, v in leaf.items() if k != "choice"}
        out["value"] = list(leaf["choice"])[0]
        return out, True
    return leaf, False


_HINT_CACHE: dict = {}


def _leaf_value(leaf, axis=None, kind=None, name=None):
    """The central value of a card/spec leaf (value, or the central value of a range, or the first choice)."""
    leaf, _hit = _collapse_leaf(leaf, axis, kind, name)
    return leaf.get("value") if isinstance(leaf, dict) else leaf


def _same_value(a, b) -> bool:
    if isinstance(a, bool) or isinstance(b, bool) or not isinstance(a, (int, float)) or not isinstance(b, (int, float)):
        return a == b
    return math.isclose(float(a), float(b), rel_tol=1e-9, abs_tol=1e-12)


def _display_hint_cards() -> list:
    """[(material kind, structure kind, {param: central value}, display_as, card id)] for the paper cards that carry
    meta.display_as (an illustrative 3D geometry for a lumped structure, R2-3D); re-read when a card file changes."""
    import yaml
    files = sorted(PAPERS_DIR.glob("*.yaml"))
    sig = tuple((f.name, f.stat().st_mtime_ns) for f in files)
    if _HINT_CACHE.get("sig") == sig:
        return _HINT_CACHE["cards"]
    out = []
    for f in files:
        if f.name.startswith("_"):
            continue
        try:
            doc = yaml.safe_load(f.read_text(encoding="utf-8")) or {}
            da = (doc.get("meta") or {}).get("display_as")
            dev = doc.get("device") or {}
            if not da:
                continue
            st = dev.get("structure") or {}
            out.append(((dev.get("material") or {}).get("kind"), st.get("kind"),
                        {k: _leaf_value(v, "structure", st.get("kind"), k) for k, v in (st.get("params") or {}).items()},
                        da, f.stem))
        except Exception:  # noqa: BLE001 -- a broken card must not break the scene route
            continue
    _HINT_CACHE.update(sig=sig, cards=out)
    return out


def display_hint(spec):
    """{"display_as": X} when `spec` is (the central design of) a paper card whose meta.display_as names an illustrative
    geometry for its lumped structure; else None.  The composer posts only the spec, so the hint is recovered by matching
    the material kind and the structure kind + parameters against the cards.  Editing a structure parameter drops the
    hint (the lumped schematic is shown); when several cards match, they must agree."""
    try:
        st = spec.get("structure") or {}
        mat = (spec.get("material") or {}).get("kind")
        post = {k: _leaf_value(v, "structure", st.get("kind"), k) for k, v in (st.get("params") or {}).items()}
        hits = set()
        for m, k, params, da, _cid in _display_hint_cards():
            if m == mat and k == st.get("kind") and params and all(n in post and _same_value(post[n], v) for n, v in params.items()):
                hits.add(da)
        return {"display_as": hits.pop()} if len(hits) == 1 else None
    except Exception:  # noqa: BLE001
        return None


def central_design(paper_id: str, device: dict) -> dict:
    spec = copy.deepcopy(device)
    collapsed = []
    spec["T_K"], hit = _collapse_leaf(spec["T_K"], None, None, "T_K")
    if hit:
        collapsed.append("T_K")
    for axis in ("material", "structure", "drive", "setup"):
        params = (spec.get(axis) or {}).get("params") or {}
        for name in list(params):
            params[name], hit = _collapse_leaf(params[name], axis, spec[axis].get("kind"), name)
            if hit:
                collapsed.append(f"{axis}.{name}")
    return {"design": {"name": f"paper-{paper_id}", "platform": "composed", "composed": spec},
            "paper": paper_id, "collapsed": collapsed, "note": CENTRAL_NOTE}


def list_examples() -> list:
    """[{name, label, summary, kinds, design}] for every cards/composed/*.yaml (sorted by file name)."""
    import yaml
    out = []
    for path in sorted(EXAMPLES_DIR.glob("*.yaml")):
        doc = yaml.safe_load(path.read_text(encoding="utf-8"))
        design = doc["design"]
        meta = doc.get("meta") or {}
        spec = design["composed"]
        out.append({"name": design.get("name", path.stem), "label": meta.get("label", path.stem),
                    "summary": meta.get("summary", ""),
                    "kinds": {a: spec[a]["kind"] for a in ("material", "structure", "drive", "setup")},
                    "design": {"name": design.get("name", path.stem), "platform": "composed", "composed": spec}})
    return out


def fock_design() -> dict:
    """The 4 K Fock source card mapped onto a composed spec (tiers.fock_card_to_spec)."""
    tiers = _optional("tiers")
    if tiers is None or not hasattr(tiers, "fock_card_to_spec"):
        raise ComposeError("fsim_core.compose.tiers.fock_card_to_spec is not available", 501)
    return {"design": {"name": "fock-source-4K", "platform": "composed", "composed": tiers.fock_card_to_spec()},
            "tier": "fock", "note": "Fock 4 K resonant source card mapped onto the generic engine"}


# ------------------------------------------------------------------ job bodies (run in the pool)
def series_job(spec, path, values, observables):
    from fsim_core.compose import engine
    res = engine.series(spec, path, [float(v) for v in values], list(observables))
    return {"series": res, "spec_path": path, "n": len(values)}


def _model_series(paper, spec) -> list:
    """The engine's model curve for every paper series entry (central point of every other leaf)."""
    from fsim_core.compose import engine
    out = []
    for s in _pget(paper, "series", []) or []:
        sid = s.get("id") if isinstance(s, dict) else None
        try:
            xpath = s["x"]["path"]
            yobs = s["y"]["observable"]
            conv = s["y"].get("convention")
            req = f"{yobs}:{conv}" if conv else yobs
            xs = [float(p[0]) for p in s.get("points", [])]
            if not xs:
                out.append({"id": sid, "error": "series has no x points"})
                continue
            res = engine.series(spec, xpath, xs, [req])
            out.append({"id": sid, "path": xpath, "request": req, "x": res["x"], "y": res["y"][req],
                        "tag": res["tag"][req], "unit": res.get("unit"), "failed": res["failed"],
                        "note": "model at the central value of every other leaf"})
        except Exception as exc:  # noqa: BLE001 -- one bad series must not lose the report
            out.append({"id": sid, "error": f"{type(exc).__name__}: {exc}"})
    return out


def reproduce_job(paper_id, card_sha256=None, papers_dir=None):
    import dataclasses
    from pathlib import Path
    papers = _need_papers()
    path = (Path(papers_dir) if papers_dir else PAPERS_DIR) / _card_path(paper_id).name
    if card_sha256 is not None and _card_sha(path) != card_sha256:
        raise ValueError(f"card {paper_id!r} changed since the job was submitted")
    paper = papers.load_paper(path)
    report = papers.reproduce(paper)
    if dataclasses.is_dataclass(report) and not isinstance(report, type):
        rep = dataclasses.asdict(report)
    elif hasattr(report, "to_dict"):
        rep = report.to_dict()
    else:
        rep = report
    return {"paper": paper_id, "card_sha256": _card_sha(path), "report": rep,
            "items": rep.get("items") if isinstance(rep, dict) else None,
            "model_series": _model_series(paper, _pget(paper, "device"))}


# ------------------------------------------------------------------ routes
def register(app):
    from flask import Response, request

    from .jobs import UnsupportedRun
    from .serialize import dumps

    def _out(obj, status=200):
        return Response(dumps(obj), status=status, mimetype="application/json")

    @app.errorhandler(ComposeError)
    def _compose_error(exc):
        return _out({"error": exc.message, **exc.extra}, exc.status)

    def jm():
        return app.config["JOB_MANAGER"]

    def _body():
        data = request.get_json(silent=True)
        if not isinstance(data, dict):
            raise ComposeError("request body must be a JSON object")
        return data

    def _submit(fn, kwargs, eta, label):
        try:
            return jm().submit_call(fn, kwargs, eta, label=label)
        except ValueError as exc:
            raise ComposeError(str(exc)) from exc

    # ---------------------------------------------------------------- composer
    @app.get("/api/compose/registry")
    def compose_registry():
        return _out(registry_payload_cached())

    @app.post("/api/compose/validate")
    def compose_validate():
        body = _body()
        if "spec" not in body:
            raise ComposeError("body must carry spec")
        return _out(validate_payload(body["spec"]))

    @app.post("/api/compose/run")
    def compose_run():
        from .api_cards import design_from_dict
        body = _body()
        spec, _res = _valid_spec(body)
        mode = body.get("mode", "point")
        if mode not in ("point", "tsweep"):
            raise ComposeError("mode must be point or tsweep")
        run_body = {"mode": "point"}
        if mode == "tsweep":
            grid = body.get("T_grid")
            if grid is None:
                raise ComposeError("tsweep needs T_grid: [T, ...] or {lo, hi, n}")
            n = grid.get("n") if isinstance(grid, dict) else (len(grid) if isinstance(grid, list) else None)
            if not isinstance(n, int) or isinstance(n, bool) or not 1 <= n <= MAX_TSWEEP_POINTS:
                raise ComposeError(f"tsweep T_grid must carry 1..{MAX_TSWEEP_POINTS} points")
            run_body["T_grid"] = grid
        name = body.get("name") or "composed"
        design = design_from_dict({"name": name, "platform": "composed", "composed": spec})
        try:
            resp = jm().submit_run(run_body, design, name)
        except UnsupportedRun as exc:
            return _out({"error": str(exc), **exc.payload}, 400)
        except (ValueError, KeyError, TypeError) as exc:
            raise ComposeError(str(exc)) from exc
        return _out(resp)

    @app.post("/api/compose/series")
    def compose_series():
        body = _body()
        spec, res = _valid_spec(body)
        path = body.get("path")
        if not isinstance(path, str) or path not in res["resolved"]:
            raise ComposeError(f"path must be one of the spec's leaf paths: {sorted(res['resolved'])}")
        values = body.get("values")
        if (not isinstance(values, list) or not 1 <= len(values) <= MAX_SERIES_POINTS
                or any(isinstance(v, bool) or not isinstance(v, (int, float)) or not math.isfinite(v)
                       for v in values)):
            raise ComposeError(f"values must be a list of 1..{MAX_SERIES_POINTS} finite numbers")
        observables = _check_observables(body.get("observables"))
        resp = _submit(SERIES_FN, {"spec": spec, "path": path, "values": [float(v) for v in values],
                                   "observables": observables}, series_eta(spec, len(values)), "series")
        return _out(resp, 200 if resp.get("cached") else 202)

    def _scene(spec, display):
        try:
            mod = importlib.import_module("fsim_studio.scene_compose")
        except ModuleNotFoundError as exc:
            if exc.name == "fsim_studio.scene_compose":
                raise ComposeError("fsim_studio.scene_compose is not available yet", 501) from exc
            raise
        if not isinstance(spec, dict):
            raise ComposeError("spec must be a composed spec object")
        if display is not None and not isinstance(display, dict):
            raise ComposeError("display must be an object")
        if not (display or {}).get("display_as"):
            hint = display_hint(spec)                 # R2-3D: a paper card's illustrative geometry for a lumped structure
            if hint:
                display = {**(display or {}), **hint}
        return _out(mod.describe(spec, display))

    @app.get("/api/compose/scene")
    def compose_scene_get():
        raw = request.args.get("spec")
        if not raw:
            raise ComposeError("spec query parameter (JSON) is required")
        try:
            spec = json.loads(raw)
            raw_display = request.args.get("display")
            display = json.loads(raw_display) if raw_display else None
        except ValueError as exc:
            raise ComposeError("spec and display must be JSON") from exc
        return _scene(spec, display)

    @app.post("/api/compose/scene")
    def compose_scene_post():
        body = _body()
        return _scene(body.get("spec"), body.get("display"))

    @app.get("/api/compose/examples")
    def compose_examples():
        return _out(list_examples())

    @app.get("/api/compose/tier/fock")
    def compose_tier_fock():
        return _out(fock_design())

    # ---------------------------------------------------------------- papers
    @app.get("/api/papers")
    def papers_list():
        papers = _need_papers()
        out = []
        for p in papers.list_papers(str(PAPERS_DIR)):
            pid = (_pget(p, "meta", {}) or {}).get("id")
            if not pid:
                continue
            try:
                out.append(paper_summary(p, _card_path(pid)))
            except ComposeError:
                continue
        return _out(out)

    @app.get("/api/papers/<paper_id>")
    def paper_get(paper_id):
        from fsim_core.compose import validate_spec
        paper, path = _load_paper(paper_id)
        device = _pget(paper, "device", {})
        return _out({**paper_summary(paper, path), "meta": _pget(paper, "meta"), "device": device,
                     "resolved": validate_spec(device)["resolved"], "reported": _pget(paper, "reported", []),
                     "series": _pget(paper, "series", []), "reproduce": _pget(paper, "reproduce", []),
                     "fits": _pget(paper, "fits", []),
                     "inputs_from_paper": sorted(_pget(paper, "inputs_from_paper", []) or []),
                     "card_sha256": _card_sha(path)})

    @app.post("/api/papers/<paper_id>/reproduce")
    def paper_reproduce(paper_id):
        _paper, path = _load_paper(paper_id)
        resp = _submit(REPRODUCE_FN, _reproduce_kwargs(paper_id, path), 5.0, f"paper-{paper_id}")
        return _out(resp, 200 if resp.get("cached") else 202)

    @app.get("/api/papers/<paper_id>/design")
    def paper_design(paper_id):
        paper, _path = _load_paper(paper_id)
        return _out(central_design(paper_id, _pget(paper, "device")))


__all__ = ["register", "ComposeError", "display_hint", "series_job", "reproduce_job", "validate_payload",
           "registry_payload_cached", "central_design", "list_examples", "fock_design"]
