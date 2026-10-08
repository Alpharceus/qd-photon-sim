"""Flask app factory for FSIM Studio. Routes only; logic lives in api_cards,
jobs, campaigns, stories (and scene, owned by the 3D coder)."""
from __future__ import annotations

import importlib
import json
import sys
import time

from flask import Flask, Response, abort, request, send_from_directory, stream_with_context

from . import OUT, ROOT, WEB, __version__, fsim_core_hash
from .serialize import dumps

if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))


def _json(obj, status: int = 200) -> Response:
    return Response(dumps(obj), status=status, mimetype="application/json")


def _err(status: int, message: str, **extra) -> Response:
    return _json({"error": message, **extra}, status)


def _body() -> dict:
    data = request.get_json(silent=True)
    if not isinstance(data, dict):
        abort(_err(400, "request body must be a JSON object"))
    return data


MUTATING = ("POST", "PUT", "DELETE", "PATCH")
STUDIO_HEADER = "X-FSIM-Studio"


def allowed_hosts(port: int | None) -> frozenset:
    """Host header allowlist (H13, DNS-rebinding guard): 127.0.0.1:<port>
    and localhost:<port>; without a port (test client) the bare names."""
    if port:
        return frozenset({f"127.0.0.1:{port}", f"localhost:{port}"})
    return frozenset({"127.0.0.1", "localhost"})


def create_app(job_manager=None, port: int | None = None) -> Flask:
    from . import api_cards, campaigns, stories
    from .jobs import JobManager, UnsupportedRun

    app = Flask(__name__, static_folder=None)
    jm = job_manager or JobManager()
    app.config["JOB_MANAGER"] = jm
    hosts = allowed_hosts(port)
    app.config["ALLOWED_HOSTS"] = hosts

    @app.before_request
    def _guard():
        host = (request.headers.get("Host") or "").strip().lower()
        if host not in hosts:
            return _err(403, f"host {host!r} not allowed", allowed=sorted(hosts))
        if request.method in MUTATING and request.headers.get(STUDIO_HEADER) != "1":
            return _err(403, f"{request.method} requires the header {STUDIO_HEADER}: 1")
        return None

    @app.errorhandler(api_cards.CardError)
    def _card_error(exc):
        return _err(exc.status, exc.message)

    @app.errorhandler(campaigns.CampaignError)
    def _campaign_error(exc):
        return _err(exc.status, exc.message)

    @app.errorhandler(stories.StoryError)
    def _story_error(exc):
        return _err(exc.status, exc.message)

    # ------------------------------------------------------------ basics
    @app.get("/api/health")
    def health():
        return _json({"ok": True, "version": __version__, "fsim_core_hash": fsim_core_hash()})

    @app.get("/api/theme/plotly")
    def theme_plotly():
        import fsim_theme
        mode = request.args.get("mode", "light")
        if mode not in ("light", "dark"):
            return _err(400, "mode must be light or dark")
        try:
            tpl = fsim_theme.plotly_template(mode).to_plotly_json()
        except ImportError:
            tpl = json.loads(fsim_theme.plotly_template_path(mode).read_text(encoding="utf-8"))
        return _json(tpl)

    # ------------------------------------------------------------ cards
    @app.get("/api/cards")
    def cards():
        return _json(api_cards.list_cards())

    @app.get("/api/cards/<name>")
    def card_get(name):
        return _json(api_cards.get_card(name))

    @app.get("/api/cards/<name>/verdict")
    def card_verdict(name):
        return _json(campaigns.card_verdict(name))

    @app.get("/api/gates")
    def gates():
        return _json(campaigns.gates())

    @app.put("/api/cards/<name>")
    def card_put(name):
        body = _body()
        if "design" not in body:
            return _err(400, "body must carry {design}")
        return _json(api_cards.save_card(name, body["design"], body.get("presets")))

    @app.get("/api/meta")
    def meta():
        return _json(api_cards.meta_payload())

    @app.post("/api/presets/apply")
    def presets_apply():
        return _json(api_cards.apply_presets(_body()))

    @app.post("/api/validate")
    def validate():
        body = _body()
        if "design" not in body:
            return _err(400, "body must carry {design}")
        return _json(api_cards.validate(body["design"]))

    # ------------------------------------------------------------ runs / jobs
    @app.post("/api/run")
    def run():
        body = _body()
        card = body.get("card")
        if body.get("design") is not None:
            design = api_cards.design_from_dict(body["design"])
            card = card or design.name
        elif card:
            design = api_cards.load_design(card)
        else:
            return _err(400, "body must carry card or design")
        try:
            resp = jm.submit_run(body, design, card)
        except UnsupportedRun as exc:
            return _json({"error": str(exc), **exc.payload}, 400)
        except (ValueError, KeyError, TypeError) as exc:
            return _err(400, str(exc))
        return _json(resp)

    # ------------------------------------------------------------ response animation (studio-p2a)
    @app.post("/api/animate")
    def animate():
        from . import animate as animate_mod
        body = _body()
        card = body.get("card")
        if body.get("design") is not None:
            design = api_cards.design_from_dict(body["design"])
            card = card or design.name
        elif card:
            design = api_cards.load_design(card)
        else:
            return _err(400, "body must carry card or design")
        try:
            resp = animate_mod.submit_animate(jm, body, design, card)
        except UnsupportedRun as exc:
            return _json({"error": str(exc), **exc.payload}, 400)
        except (ValueError, KeyError, TypeError) as exc:
            return _err(400, str(exc))
        return _json(resp)

    @app.get("/api/jobs/<job_id>")
    def job_get(job_id):
        job = jm.get(job_id)
        if job is None:
            return _err(404, f"no job {job_id!r}")
        return _json(job.snapshot())

    @app.delete("/api/jobs/<job_id>")
    def job_cancel(job_id):
        snap = jm.cancel(job_id)
        if snap is None:
            return _err(404, f"no job {job_id!r}")
        return _json(snap)

    @app.get("/api/jobs/<job_id>/events")
    def job_events(job_id):
        job = jm.get(job_id)
        if job is None:
            return _err(404, f"no job {job_id!r}")
        resp = Response(stream_with_context(jm.stream(job)), mimetype="text/event-stream")
        resp.headers["Cache-Control"] = "no-cache"
        resp.headers["X-Accel-Buffering"] = "no"
        return resp

    @app.post("/api/compare/report")
    def compare_report():
        from .report import build_report
        body = _body()
        entries = body.get("entries")
        if not isinstance(entries, list) or not entries:
            return _err(400, "body must carry entries: [{label, run_id}]")
        try:
            path = build_report(entries, title=body.get("title", "design review"))
        except LookupError as exc:
            return _err(404, str(exc))
        return _json({"path": str(path), "files": sorted(p.name for p in path.iterdir())})

    # ------------------------------------------------------------ campaigns
    @app.get("/api/campaigns")
    def campaign_list():
        return _json(campaigns.list_campaigns())

    @app.get("/api/campaigns/<path:cid>/sweep")
    def campaign_sweep(cid):
        filt = request.args.get("filter")
        try:
            filt = json.loads(filt) if filt else None
        except ValueError:
            return _err(400, "filter must be JSON")
        cols = request.args.get("cols")
        cols = [c for c in cols.split(",") if c] if cols else None
        try:
            limit = int(request.args.get("limit", 5000))
            offset = int(request.args.get("offset", 0))
        except ValueError:
            return _err(400, "limit/offset must be integers")
        headline = request.args.get("headline", "0") in ("1", "true", "True")
        return _json(campaigns.query_sweep(cid, file=request.args.get("file"), filt=filt,
                                           cols=cols, limit=limit, headline=headline,
                                           offset=offset))

    @app.get("/api/campaigns/<path:cid>")
    def campaign_get(cid):
        return _json(campaigns.get_campaign(cid))

    @app.get("/out/<path:path>")
    def out_file(path):
        return send_from_directory(OUT, path)

    # ------------------------------------------------------------ scenes (3D)
    def _scene(fn_name: str, *args, **kwargs):
        try:
            scene = importlib.import_module("fsim_studio.scene")
        except ModuleNotFoundError as exc:
            if exc.name == "fsim_studio.scene":
                return _err(501, "fsim_studio.scene is not available yet")
            raise
        fn = getattr(scene, fn_name, None)
        if fn is None:
            return _err(501, f"fsim_studio.scene.{fn_name} is not implemented")
        try:
            spec = fn(*args, **kwargs)
        except (api_cards.CardError, campaigns.CampaignError):
            raise
        except (ValueError, KeyError, TypeError, FileNotFoundError, LookupError) as exc:
            return _err(400, f"{type(exc).__name__}: {exc}")
        if isinstance(spec, dict) and spec.get("job") and spec.get("fn"):
            # slow scene: the builder returned a job descriptor; run it in the pool
            try:
                resp = jm.submit_call(spec["fn"], spec.get("kwargs") or {}, spec.get("eta_s", 0.0),
                                      label=f"scene-{spec.get('kind', 'x')}")
            except ValueError as exc:
                return _err(400, str(exc))
            return _json({**resp, "scene_job": True, "descriptor": spec}, 202 if not resp["cached"] else 200)
        return _json(spec)

    def _float_arg(name):
        v = request.args.get(name)
        if v in (None, ""):
            return None
        try:
            return float(v)
        except ValueError:
            abort(_err(400, f"{name} must be a number"))

    @app.get("/api/scene/device")
    def scene_device():
        return _scene("device_scene", request.args.get("card"), _float_arg("T"))

    @app.get("/api/scene/band")
    def scene_band():
        return _scene("band_scene", request.args.get("card"), _float_arg("T"))

    @app.get("/api/scene/cascade")
    def scene_cascade():
        return _scene("cascade_scene", request.args.get("card"), _float_arg("T"))

    @app.get("/api/scene/surface")
    def scene_surface():
        return _scene("surface_scene", request.args.get("card"), request.args.get("param") or "w",
                      request.args.get("mode") or "envelope")

    @app.get("/api/scene/lattice")
    def scene_lattice():
        return _scene("lattice_scene", request.args.get("campaign"), request.args.get("x"),
                      request.args.get("y"), request.args.get("z"),
                      headline=request.args.get("headline", "1") not in ("0", "false", "False"))

    # ------------------------------------------------------------ stories
    @app.get("/api/stories")
    def story_list():
        return _json(stories.list_stories())

    @app.get("/api/stories/<name>")
    def story_get(name):
        return _json(stories.get_story(name))

    @app.put("/api/stories/<name>")
    def story_put(name):
        return _json(stories.save_story(name, _body()))

    @app.post("/api/stories/<name>/bake")
    def story_bake(name):
        t0 = time.time()
        story = stories.bake(name, run_live=request.args.get("live", "1") != "0")
        return _json({**story, "bake_seconds": time.time() - t0})

    # ------------------------------------------------------------ static web
    @app.get("/")
    def index():
        return send_from_directory(WEB, "index.html")

    @app.get("/<path:path>")
    def static_file(path):
        if path.startswith("api/"):
            return _err(404, f"no endpoint /{path}")
        return send_from_directory(WEB, path)

    # ------------------------------------------------------------ explain (api_explain.py)
    from . import api_explain
    api_explain.register(app)
    # ------------------------------------------------------------ library (api_library.py)
    from . import api_library
    api_library.register(app)
    # ------------------------------------------------------------ explorers A (api_explore_a.py)
    from . import api_explore_a
    api_explore_a.register(app)
    # ------------------------------------------------------------ explorers B (api_explore_b.py)
    from . import api_explore_b
    api_explore_b.register(app)
    return app
