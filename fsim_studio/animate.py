"""Response animation: one evaluate() per parameter value, computed once on the
process pool and played back in the browser (POST /api/animate).

No physics here. A frame is exactly the /api/run point (or headline) result of
the design with one numeric META path set to the frame's value: the same
switch selection, the same T grid rule and the same cache key, so a frame
computed here is a cache hit for /api/run and vice versa, and a repeated
animation costs nothing. The browser interpolates between COMPUTED frames for
display only; every number a frame carries comes from fsim_core.

This module must not import flask (it shares the job objects with jobs.py).
"""
from __future__ import annotations

import concurrent.futures as cf
import copy
import math
import threading

import numpy as np

from . import cache
from .jobs import (NITRIDE_DISPLAY_GRID, Job, UnsupportedRun, _w_evaluate, classify_error,
                   eval_seconds, headline_design, make_run_id, platform_class, point_result)

MAX_FRAMES = 120
DEFAULT_FRAMES = 48
DEFAULT_FRAMES_NANOWIRE = 16
# The frame scalars the player reads (spec studio-p2a section 2). Each is
# copied verbatim from the run's scalars, never interpolated server-side.
FRAME_SCALARS = ("g2_op", "collected_flux_pulsed_s", "collected_flux_delivered_s",
                 "brightness_per_pulse", "T_j_op", "eps_op", "rho_op", "T_c")
# Text that travels with the numbers (shown verbatim, never parsed for physics).
FRAME_TEXT = ("brightness_convention", "g2_op_valid", "g2_op_invalid_reason", "runaway",
              "flux_measurable", "rep_rate_hz", "one_pair_valid", "set_feasible",
              "flat_band", "depletion_regime")
# Animation honesty rules (physics brief section 7), enforced here and in the browser player:
#  - T_c is never taken from a grid "op" frame (rule 4); the frame says why it is null.
#  - runaway ends the displayable range (rule 5): result.runaway_at, frames after it flagged.
#  - frame density (rule 6): Arrhenius step near T_c, nitride Stark one frame per STARK_STEP_V.
OP_GRID_TC_REASON = "operating-point frame (grid op): T_c needs the full T_hs grid and is never read from an op frame"
ARRHENIUS_STEP_K = 5.0      # rule 6: no coarser than 5 K within +-20 K of T_c
ARRHENIUS_WINDOW_K = 20.0
STARK_STEP_V = 0.2          # rule 6: nitride Stark, one frame per 0.2 V
STARK_PARAM = "drive.V"
A_SLICE_NOTE = "one range at a time, not a prediction"   # rule 7


def numeric_meta(path: str) -> dict:
    """META entry for an animatable path: numeric with a [lo, hi] band."""
    from fsim_core import design_meta
    m = design_meta.META.get(path)
    if m is None:
        raise ValueError(f"param {path!r} is not a design_meta.META path")
    if m.get("unit") in ("bool", "enum") or "choices" in m or "lo" not in m or "hi" not in m:
        raise ValueError(f"param {path!r} is not a numeric META field with a [lo, hi] band")
    return m


def default_frames(cls: str, platform: str = "", param: str = "", lo: float = 0.0, hi: float = 0.0) -> int:
    """Frame count when the request names none. Nitride Stark bias (drive.V on an InGaN/GaN
    platform): one frame per STARK_STEP_V (brief rule 6); otherwise 16 for nanowire, 48 else."""
    if param == STARK_PARAM and platform.startswith("ingan_gan") and hi != lo:
        return max(2, min(MAX_FRAMES, int(round(abs(hi - lo) / STARK_STEP_V)) + 1))
    return DEFAULT_FRAMES_NANOWIRE if cls == "nanowire" else DEFAULT_FRAMES


def frame_values(lo: float, hi: float, n: int) -> list:
    return [float(x) for x in np.linspace(lo, hi, n)]


def frame_T_grid(d, param: str, value: float, grid: str):
    """The T grid one frame is evaluated on. "op": the operating point only
    ([thermal.T_hs] of the frame design); "full": what /api/run uses when no
    grid is sent (fsim_core's default, or the nitride display grid)."""
    if grid == "op":
        return [float(value) if param == "thermal.T_hs" else float(d.thermal.T_hs)]
    if d.platform in NITRIDE_DISPLAY_GRID:
        lo, hi, n = NITRIDE_DISPLAY_GRID[d.platform]
        return [float(x) for x in np.linspace(lo, hi, n)]
    return None


def frame_of(index: int, value: float, result: dict, grid: str = "full") -> dict:
    """The player's view of one /api/run-shaped result. `grid` is the T grid the frame was
    evaluated on; an "op" frame never carries a T_c (rule 4), whatever the result holds."""
    s = result.get("scalars") or {}
    scalars = {}
    for k in FRAME_SCALARS:
        if k in s:
            scalars[k] = s[k]
            if f"{k}__nan_reason" in s:
                scalars[f"{k}__nan_reason"] = s[f"{k}__nan_reason"]
    if grid == "op":
        scalars["T_c"] = None
        scalars["T_c__nan_reason"] = OP_GRID_TC_REASON
    for k in FRAME_TEXT:
        if k in s:
            scalars[k] = s[k]
    out = {"index": index, "value": value, "grid": grid, "scalars": scalars,
           "tag_chain": result.get("tag_chain"), "tag_chain_source": result.get("tag_chain_source"),
           "labels": result.get("labels", []), "run_id": result.get("run_id"),
           "provenance": result.get("provenance", {})}
    curves = result.get("curves") or {}
    if len(curves.get("T_hs") or []) > 1 and "g2" in curves:
        out["curves"] = {"T_hs": curves["T_hs"], "g2": curves["g2"]}
    return out


def _num(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def density_report(param: str, platform: str, values: list, frames: list) -> dict:
    """Frame-density check of rule 6 (arithmetic on the frame values, no physics). Advisory:
    the player shows it and offers the narrowed lo/hi; nothing is refused."""
    if param == "thermal.T_hs":
        tc = next((f["scalars"]["T_c"] for f in frames
                   if f and "error" not in f and f.get("grid") == "full" and _num(f["scalars"].get("T_c"))), None)
        if tc is None:
            return {"rule": "arrhenius", "checked": False, "required_step": ARRHENIUS_STEP_K,
                    "reason": "no frame carries a T_c from the full T_hs grid, so the +-20 K window is unknown"}
        w0, w1 = tc - ARRHENIUS_WINDOW_K, tc + ARRHENIUS_WINDOW_K
        steps = [(values[i], values[i + 1]) for i in range(len(values) - 1)
                 if max(values[i], values[i + 1]) >= w0 and min(values[i], values[i + 1]) <= w1]
        worst = max((abs(b - a) for a, b in steps), default=None)
        ok = worst is not None and worst <= ARRHENIUS_STEP_K + 1e-9
        n_sug = int(round(2 * ARRHENIUS_WINDOW_K / ARRHENIUS_STEP_K)) + 1
        return {"rule": "arrhenius", "checked": True, "T_c": tc, "window": [w0, w1], "required_step": ARRHENIUS_STEP_K,
                "max_step": worst, "ok": ok, "suggest": None if ok else {"lo": w0, "hi": w1, "n": n_sug}}
    if param == STARK_PARAM and platform.startswith("ingan_gan"):
        worst = max((abs(values[i + 1] - values[i]) for i in range(len(values) - 1)), default=None)
        return {"rule": "stark", "checked": True, "required_step": STARK_STEP_V, "max_step": worst,
                "ok": worst is not None and worst <= STARK_STEP_V + 1e-9}
    return {"rule": None, "checked": False, "reason": "no density rule for this parameter"}


def annotate_frames(frames: list) -> dict:
    """Rule 5. Returns {runaway_at, frames}: the first frame whose evaluate() reports runaway
    ends the displayable range; that frame and every later one is flagged after_runaway
    (copies, the cached frames are untouched)."""
    r = next((i for i, f in enumerate(frames) if f and "error" not in f and (f.get("scalars") or {}).get("runaway") is True), None)
    out = []
    for i, f in enumerate(frames):
        out.append({**f, "after_runaway": True} if (r is not None and f is not None and i >= r) else f)
    return {"runaway_at": r, "frames": out}


def submit_animate(jm, body: dict, design, card: str) -> dict:
    """POST /api/animate body: {card?|design, param, lo?, hi?, n?, mode: point|headline,
    grid?: op|full}. Returns {job_id, cached, eta_s, n, values, result?}."""
    from fsim_core.device import _set_path
    from .api_cards import design_to_dict

    mode = body.get("mode", "point")
    if mode not in ("point", "headline"):
        raise ValueError(f"animate mode must be point or headline, got {mode!r}")
    param = body.get("param") or "thermal.T_hs"
    m = numeric_meta(param)
    block, field = param.split(".", 1)
    cur = getattr(getattr(design, block, None), field, "missing")
    if cur == "missing" or isinstance(cur, bool) or not (cur is None or isinstance(cur, (int, float))):
        raise ValueError(f"param {param!r} is not a numeric field of this design")
    lo = float(body["lo"]) if body.get("lo") is not None else float(m["lo"])
    hi = float(body["hi"]) if body.get("hi") is not None else float(m["hi"])
    if not (math.isfinite(lo) and math.isfinite(hi)) or lo == hi:
        raise ValueError("lo and hi must be finite and differ")
    # Review finding 1: the animated range must sit inside the META validity band (a frame
    # outside it would still come back flagged valid). Refused, not clamped.
    if min(lo, hi) < float(m["lo"]) or max(lo, hi) > float(m["hi"]):
        raise ValueError(f"lo/hi must lie inside the {param} validity band [{m['lo']}, {m['hi']}]")
    d_cls = headline_design(design) if mode == "headline" else design
    cls = platform_class(d_cls)
    n_req = body.get("n")
    if n_req is None:
        n = default_frames(cls, design.platform, param, lo, hi)
    else:
        # n must be an integer: 0, 1.5, "7", True are refused (n:0 used to fall back to the default)
        if isinstance(n_req, bool) or not isinstance(n_req, (int, float)) or not math.isfinite(n_req) or n_req != int(n_req):
            raise ValueError(f"n must be an integer in 2..{MAX_FRAMES}, got {n_req!r}")
        n = int(n_req)
    if n < 2 or n > MAX_FRAMES:
        raise ValueError(f"n must be in 2..{MAX_FRAMES}")
    grid = body.get("grid") or ("op" if cls in ("edge", "nanowire") else "full")
    if grid not in ("op", "full"):
        raise ValueError("grid must be op or full")

    values = frame_values(lo, hi, n)
    plan = []  # (index, value, frame design, T_grid, key, run_id)
    for i, v in enumerate(values):
        d = copy.deepcopy(design)
        _set_path(d, param, v)
        T_grid = frame_T_grid(d, param, v, grid)
        key = cache.make_key({"design": design_to_dict(d), "mode": mode, "ranged": None, "T_grid": T_grid})
        plan.append((i, v, d, T_grid, key, make_run_id(card, key)))

    meta = {"param": param, "lo": lo, "hi": hi, "n": n, "values": values, "mode": mode, "grid": grid,
            "card": card, "platform": design.platform, "unit": m.get("unit"), "param_tag": m.get("tag"),
            "slice_note": A_SLICE_NOTE if m.get("tag") == "A" else None}
    frames = [None] * n
    todo = []
    for i, v, d, T_grid, key, run_id in plan:
        rec = cache.get(key)
        if rec is not None and isinstance(rec.get("result"), dict):
            frames[i] = frame_of(i, v, rec["result"], grid)
        else:
            todo.append((i, v, d, T_grid, key, run_id))

    n_T = len(todo[0][3]) if todo and todo[0][3] is not None else 120
    eta = (eval_seconds(cls, n_T) * len(todo) / jm.max_workers) if todo else 0.0
    anim_key = cache.make_key({"animate": meta, "design": design_to_dict(design)})
    job = jm._register(Job("animate", eta, anim_key, f"FS-anim-{anim_key[:6]}"))
    job.progress = (n - len(todo), n)

    def result_of():
        ann = annotate_frames(frames)
        return {**meta, "frames": ann["frames"], "runaway_at": ann["runaway_at"],
                "density": density_report(param, design.platform, values, frames),
                "computed": [i for i, f in enumerate(frames) if f is not None and "error" not in f],
                "errors": [i for i, f in enumerate(frames) if f is not None and "error" in f]}

    if not todo:
        res = result_of()
        job.state = "done"
        job.result = res
        for f in frames:
            job.events.append(("frame", {"job_id": job.id, "index": f["index"], "frame": f, "cached": True}))
        job.events.append(("state", {"state": "done", "job_id": job.id, "cached": True}))
        job.events.append(("done", {"job_id": job.id, "result": res, "cached": True}))
        return {"job_id": job.id, "cached": True, "eta_s": 0.0, "n": n, "values": values,
                "n_cached": n, "result": res}

    for f in frames:
        if f is not None:
            job.events.append(("frame", {"job_id": job.id, "index": f["index"], "frame": f, "cached": True}))

    def run():
        try:
            if not job.set_state("running"):
                return
            pool = jm.pool()
            futs = {}
            for i, v, d, T_grid, key, run_id in todo:
                d_eval = headline_design(d) if mode == "headline" else d
                fut = pool.submit(_w_evaluate, d_eval, T_grid)
                futs[fut] = (i, v, d, d_eval, T_grid, key, run_id)
                job.futures.append(fut)
            done_n = n - len(todo)
            for fut in cf.as_completed(futs):
                if job.cancelled:
                    return
                i, v, d, d_eval, T_grid, key, run_id = futs[fut]
                try:
                    res = fut.result()
                except cf.CancelledError:
                    return
                except cf.process.BrokenProcessPool:
                    raise
                except Exception as exc:  # noqa: BLE001 -- a failed frame is reported, not fatal
                    frames[i] = {"index": i, "value": v, "error": classify_error(d_eval, exc)}
                else:
                    result = point_result(d, res, mode=mode, run_id=run_id, card=card, T_grid=T_grid)
                    cache.put(key, {"key": key, "run_id": run_id, "card": card, "mode": mode,
                                    "design": design_to_dict(d), "ranged": None, "T_grid": T_grid,
                                    "result": result})
                    frames[i] = frame_of(i, v, result, grid)
                done_n += 1
                job.progress = (done_n, n)
                job.emit("frame", {"job_id": job.id, "index": i, "frame": frames[i], "cached": False})
                job.emit("progress", {"job_id": job.id, "k": done_n, "n": n})
            if job.cancelled:
                return
            job.finish(result_of())
        except cf.process.BrokenProcessPool as exc:
            jm._reset_pool()
            job.fail({"kind": "pool", "message": f"worker pool broke: {exc}"})
        except Exception as exc:  # noqa: BLE001 -- reported to the client
            job.fail({"kind": "error", "message": f"{type(exc).__name__}: {exc}"})

    threading.Thread(target=run, daemon=True).start()
    return {"job_id": job.id, "cached": False, "eta_s": eta, "n": n, "values": values,
            "n_cached": n - len(todo)}


__all__ = ["submit_animate", "UnsupportedRun", "FRAME_SCALARS", "MAX_FRAMES", "annotate_frames",
           "density_report", "default_frames", "frame_of"]
