"""Job runner: process pool, SSE events, cancel, progressive results, cache.

No physics here. A run is fsim_core.device.evaluate (point / headline) or
evaluate_envelope (envelope) on a DeviceDesign; "headline" is only the RT
edge switch selection (copy the design, drive.finite_pulse=True,
ret.tau_cap_scales_with_density=False, drive.cw=False), the same
configuration fsim_gui.designer._headline_g2 applies.

Worker functions are top-level so a Windows (spawn) process pool can pickle
them; this module must not import flask (workers import it).
"""
from __future__ import annotations

import concurrent.futures as cf
import copy
import os
import threading
import time
import uuid

import numpy as np

from . import cache
from .serialize import dumps, to_jsonable
from .tags import chain_of

STATIC_LABEL = "static (non-headline)"
HEADLINE_LABEL = ("headline model (drive.finite_pulse=true, "
                  "ret.tau_cap_scales_with_density=false)")
PRE_RETENTION_LABEL = "pre-retention"
COMMANDED_DELIVERED_LABEL = "flux: commanded (collected_flux_pulsed_s) vs delivered (collected_flux_delivered_s)"
MID_LABEL = "all ranges at midpoint (not a prediction)"
NON_HEADLINE_MODEL_LABEL = "non-headline model"
FLOOR_LABEL = "structural floor set by b_res [A]"

# The only job functions a scene descriptor may name (H16: exact allowlist).
ALLOWED_CALLS = frozenset({"fsim_studio.scene:surface_compute"})
# studio-p2d explorers B (api_explore_b.py): exact job functions, no wildcard.
ALLOWED_CALLS = ALLOWED_CALLS | frozenset({
    "fsim_studio.api_explore_b:phonon_frames", "fsim_studio.api_explore_b:transport_sweep",
    "fsim_studio.api_explore_b:stark_trace", "fsim_studio.api_explore_b:sde_frame",
    "fsim_studio.api_explore_b:sde_li", "fsim_studio.api_explore_b:sde_trace"})
MAX_T_POINTS = 400

# Measured seconds per evaluate() on this machine (ideas/01-ux-architecture.md
# section 3): (single T point, default 120-point grid); eval_seconds is linear
# in the point count between the two.
# nanowire (re-measured 2026-10-07, nitride-nanowire-vertical-pulse-design,
# warm process): 1 point 3.6 s, 2 points 7.9 s, the 16-point display grid
# 104-109 s (~6.8 s per point). The old (3.6, 8.9) entry assumed a cheap
# vectorised grid and estimated ~4 s for that 16-point run. The 120-point
# value is the same per-point slope extrapolated: 3.6 + 6.76 * 119.
ETA_TABLE = {
    "legacy": (0.001, 0.006),
    "edge": (0.31, 20.0),
    "planar": (0.001, 0.31),
    "nanowire": (1.0, 295.0),  # re-measured 2026-10-08 after the vectorized injector (1 pt 1.0 s, 120 pt 295 s)
}
INLINE_MAX_T = 200
DEFAULT_N_T = 120
SLOW_CLASSES = ("edge", "nanowire")

TERMINAL = ("done", "error", "cancelled")


# ----------------------------------------------------------- pure helpers

def platform_class(d) -> str:
    if d.platform == "ingan_gan_nanowire":
        return "nanowire"
    if d.platform == "ingan_gan_planar":
        return "planar"
    if getattr(d.emission, "type", "none") == "edge" or bool(getattr(d.drive, "finite_pulse", False)):
        return "edge"
    return "legacy"


def eval_seconds(cls: str, n_T: int) -> float:
    t1, tn = ETA_TABLE[cls]
    return t1 + (tn - t1) * max(n_T - 1, 0) / (DEFAULT_N_T - 1)


def resolve_T_grid(spec):
    """None | [T, ...] | {lo, hi, n} -> None or a list of floats."""
    if spec is None:
        return None
    if isinstance(spec, dict):
        lo, hi, n = float(spec["lo"]), float(spec["hi"]), int(spec["n"])
        if n < 1 or n > MAX_T_POINTS:
            raise ValueError(f"T_grid.n must be in 1..{MAX_T_POINTS}")
        return [float(x) for x in np.linspace(lo, hi, n)]
    vals = [float(x) for x in spec]
    if not vals:
        raise ValueError("T_grid must not be empty")
    if len(vals) > MAX_T_POINTS:
        raise ValueError(f"T_grid may carry at most {MAX_T_POINTS} points")
    return vals


def is_edge(d) -> bool:
    return getattr(d.emission, "type", "none") == "edge"


def is_headline_model(d) -> bool:
    """The RT-edge headline model switches (docs/rt_edge_contract.md model
    axes): drive.finite_pulse True and ret.tau_cap_scales_with_density
    False. drive.cw only adds CW diagnostics (g2_op is unchanged by it), so
    it is not part of the model identity."""
    return (bool(getattr(d.drive, "finite_pulse", False))
            and not bool(getattr(d.ret, "tau_cap_scales_with_density", False)))


def is_static_edge(d) -> bool:
    """Same predicate as fsim_gui.designer._is_static_edge (re-implemented:
    importing the designer would import dearpygui)."""
    return (getattr(d.emission, "type", "none") == "edge"
            and not bool(getattr(d.drive, "finite_pulse", False)))


def headline_design(d):
    """Switch selection only (designer._headline_g2): no physics."""
    dh = copy.deepcopy(d)
    dh.drive.finite_pulse = True
    dh.ret.tau_cap_scales_with_density = False
    dh.drive.cw = False
    return dh


def labels_for(d_orig, mode: str, scalars: dict) -> list:
    """Mandatory result labels. Configuration bookkeeping only; the b_res
    floor test is fsim_core.scene_support.is_background_floor."""
    from fsim_core.scene_support import is_background_floor
    labels = []
    if mode == "headline":
        labels.append(HEADLINE_LABEL)
    else:
        if is_static_edge(d_orig):
            labels.append(STATIC_LABEL)
        if is_edge(d_orig) and not is_headline_model(d_orig):
            labels.append(NON_HEADLINE_MODEL_LABEL)
    if is_background_floor(scalars.get("g2_op"), getattr(d_orig.drive, "b_res", 0.0)):
        labels.append(FLOOR_LABEL)
    conv = scalars.get("brightness_convention")
    if isinstance(conv, str) and conv.startswith("pre-retention"):
        labels.append(PRE_RETENTION_LABEL)
    if "collected_flux_delivered_s" in scalars:
        labels.append(COMMANDED_DELIVERED_LABEL)
    if mode == "envelope":
        labels.append(MID_LABEL)
    return labels


def card_short(name: str) -> str:
    s = name[:-len("-design")] if name.endswith("-design") else name
    return s or "design"


def make_run_id(card: str, key: str) -> str:
    return f"FS-{card_short(card)}-{key[:6]}"


def f8_error(d, exc: Exception) -> dict:
    from fsim_core.loading import f8b_thin_fano
    floor = 1.0 - f8b_thin_fano(d.drive.eta_capture, d.drive.F_p)
    return {"kind": "f8_domain", "message": str(exc), "f8_floor": floor,
            "suggestion": (f"Every sample sits outside the F8 loading domain (mu >= 1 - F_eff). "
                           f"Narrow drive.mu to >= {floor:.6g}, or raise drive.F_p / "
                           "drive.eta_capture.")}


NITRIDE_ALL_DROPPED_MESSAGE = (
    "Every envelope sample failed to evaluate. evaluate_envelope drops any sample that raises "
    "ValueError and reports the whole box with F8-domain wording, but this card's dot loading "
    "is resolved by the nitride transport model (drive.cycle_loading), not by drive.mu, so "
    "narrowing drive.mu is not a remedy. Run a point at the box edges to see the actual error.")


def classify_error(d, exc: Exception) -> dict:
    msg = str(exc)
    if isinstance(exc, ValueError) and "F8" in msg:
        if getattr(d, "platform", "legacy") in NITRIDE_DISPLAY_GRID:
            # cycle_loading cards: drive.mu is never read, so the F8 floor and its
            # "narrow drive.mu" remedy would be fake; report the honest message only.
            return {"kind": "envelope_all_dropped", "message": NITRIDE_ALL_DROPPED_MESSAGE,
                    "detail": msg}
        return f8_error(d, exc)
    return {"kind": "value_error" if isinstance(exc, ValueError) else "error",
            "message": f"{type(exc).__name__}: {msg}"}


class UnsupportedRun(ValueError):
    """A run the backend refuses up front (HTTP 400 {kind: "unsupported", message})."""

    def __init__(self, message: str):
        super().__init__(message)
        self.payload = {"kind": "unsupported", "message": message}


NANOWIRE_ENVELOPE_UNSUPPORTED = (
    "envelope runs are not available for the nanowire tier (its evaluate() returns g2_op/T_j "
    "curves, not the g2/Tj curves evaluate_envelope bands); use the committed sweep in Results")


# ------------------------------------------------- picklable worker functions

def _w_evaluate(design, T_grid):
    from fsim_core.device import evaluate
    return evaluate(design, T_grid=T_grid)


def _w_call(fn_path: str, kwargs: dict):
    """Generic job body for a slow scene builder: "module:function" restricted
    to the exact ALLOWED_CALLS set (the scene job descriptors)."""
    import importlib
    if fn_path not in ALLOWED_CALLS:
        raise ValueError(f"job function {fn_path!r} not allowed")
    mod_name, _, fn_name = fn_path.partition(":")
    return to_jsonable(getattr(importlib.import_module(mod_name), fn_name)(**kwargs))


# ------------------------------------------------------------- result shaping

NITRIDE_DISPLAY_GRID = {"ingan_gan_planar": (200.0, 350.0, 31), "ingan_gan_nanowire": (200.0, 350.0, 16)}


def point_result(d_orig, res: dict, *, mode, run_id, card, T_grid, partial=False) -> dict:
    scalars = res["scalars"]
    curves = dict(res["curves"])
    for std, alt in (("g2", "g2_op"), ("Tj", "T_j")):  # nanowire tier names its curves differently
        if std not in curves and alt in curves:
            curves[std] = curves[alt]
    if "T_hs" not in curves:  # nitride tiers return no T axis: label it from the grid that was run
        n = len(next(iter(curves.values()), []))
        if T_grid is not None and len(T_grid) == n:
            curves["T_hs"] = list(T_grid)
        elif n == 1:
            curves["T_hs"] = [float(d_orig.thermal.T_hs)]
    res = dict(res, curves=curves)
    tag, tag_src = chain_of(scalars)
    out = {
        "curves": res["curves"],
        "scalars": scalars,
        "tag_chain": tag,
        "tag_chain_source": tag_src,
        "provenance": scalars.get("provenance", {}),
        "labels": labels_for(d_orig, mode, scalars),
        "run_id": run_id,
        "mode": mode,
        "card": card,
        "platform": d_orig.platform,
        "T_grid": T_grid,
    }
    if partial:
        out["partial"] = True
    return to_jsonable(out)


def envelope_result(d_orig, env: dict, mid_res, *, ranged, run_id, card, T_grid) -> dict:
    Ts = np.asarray(T_grid if T_grid is not None else np.linspace(4.0, 350.0, DEFAULT_N_T))
    mid = dict(env["mid"])
    mid.setdefault("T_hs", Ts)
    bands = {name: {"lo": lo, "hi": hi, "mid": mid.get(name)}
             for name, (lo, hi) in env["bands"].items()}
    scalars = mid_res["scalars"] if mid_res is not None else {}
    tag, tag_src = chain_of(scalars)
    out = {
        "curves": mid,
        "scalars": scalars,
        "scalars_source": "mid-design point evaluate(): " + MID_LABEL,
        "tag_chain": tag,
        "tag_chain_source": tag_src,
        "provenance": scalars.get("provenance", {}),
        "bands": bands,
        "scalar_bands": {k: [lo, hi] for k, (lo, hi) in env["scalar_bands"].items()},
        "tornado": env["tornado"],
        "n_samples": env["n_samples"],
        "ranged": ranged,
        "labels": labels_for(d_orig, "envelope", scalars),
        "run_id": run_id,
        "mode": "envelope",
        "card": card,
        "platform": d_orig.platform,
        "T_grid": T_grid,
    }
    return to_jsonable(out)


def mid_design(d, ranged: dict):
    """Every ranged path pinned to the midpoint of its value set -- the same
    point evaluate_envelope's 'mid' curve uses (bookkeeping, not physics)."""
    from fsim_core.device import _set_path
    dm = copy.deepcopy(d)
    for p, spec in ranged.items():
        vals = list(spec)
        _set_path(dm, p, 0.5 * (min(vals) + max(vals)))
    return dm


# --------------------------------------------------------------------- jobs

class _Cancelled(Exception):
    pass


class Job:
    def __init__(self, kind: str, eta_s: float, key: str, run_id: str):
        self.id = uuid.uuid4().hex[:12]
        self.kind = kind
        self.eta_s = eta_s
        self.key = key
        self.run_id = run_id
        self.state = "queued"
        self.progress = None
        self.result = None
        self.partial = None
        self.error = None
        self.cancelled = False
        self.created = time.time()
        self.events = []
        self.futures = []
        self.cond = threading.Condition()

    @property
    def terminal(self) -> bool:
        return self.state in TERMINAL

    def emit(self, etype: str, data: dict) -> None:
        with self.cond:
            self.events.append((etype, data))
            self.cond.notify_all()

    def set_state(self, state: str) -> bool:
        """No-op once the job was cancelled or reached a terminal state (a
        worker setting "running" after a DELETE must not revive it)."""
        with self.cond:
            if self.cancelled or self.terminal:
                return False
            self.state = state
            self.events.append(("state", {"state": state, "job_id": self.id}))
            self.cond.notify_all()
        return True

    def finish(self, result: dict) -> None:
        with self.cond:
            if self.cancelled:
                return
            self.result = result
            self.state = "done"
            self.events.append(("state", {"state": "done", "job_id": self.id}))
            self.events.append(("done", {"job_id": self.id, "result": result}))
            self.cond.notify_all()

    def fail(self, error: dict) -> None:
        with self.cond:
            if self.cancelled:
                return
            self.error = error
            self.state = "error"
            self.events.append(("error", {"job_id": self.id, "error": error}))
            self.cond.notify_all()

    def snapshot(self) -> dict:
        out = {"job_id": self.id, "kind": self.kind, "state": self.state,
               "progress": ({"k": self.progress[0], "n": self.progress[1]}
                            if self.progress else None),
               "eta_s": self.eta_s, "run_id": self.run_id}
        if self.partial is not None and self.result is None:
            out["partial"] = self.partial
        if self.result is not None:
            out["result"] = self.result
        if self.error is not None:
            out["error"] = self.error
        return out


class JobManager:
    def __init__(self, max_workers: int | None = None):
        self.max_workers = max_workers or max(1, (os.cpu_count() or 1) - 2)
        self._pool = None
        self._pool_lock = threading.Lock()
        self.jobs: dict[str, Job] = {}
        self._jobs_lock = threading.Lock()

    # -- pool
    def pool(self) -> cf.ProcessPoolExecutor:
        with self._pool_lock:
            if self._pool is None:
                self._pool = cf.ProcessPoolExecutor(max_workers=self.max_workers)
            return self._pool

    def _reset_pool(self) -> None:
        with self._pool_lock:
            old, self._pool = self._pool, None
        if old is not None:
            old.shutdown(wait=False, cancel_futures=True)

    def shutdown(self) -> None:
        with self._pool_lock:
            old, self._pool = self._pool, None
        if old is not None:
            old.shutdown(wait=True, cancel_futures=True)

    def terminate(self) -> None:
        """Server exit: stop the pool without waiting for running evaluations
        (a cancelled job's future may still be running in a worker)."""
        with self._pool_lock:
            old, self._pool = self._pool, None
        if old is None:
            return
        procs = list(getattr(old, "_processes", {}).values())
        old.shutdown(wait=False, cancel_futures=True)
        for proc in procs:
            try:
                proc.terminate()
            except Exception:  # noqa: BLE001 -- best effort at exit
                pass

    def _register(self, job: Job) -> Job:
        with self._jobs_lock:
            self.jobs[job.id] = job
        return job

    def get(self, job_id: str):
        return self.jobs.get(job_id)

    # -- public entry point
    def submit_run(self, body: dict, design, card: str) -> dict:
        """design: a DeviceDesign already loaded (from a card or a dict)."""
        from fsim_core import design_meta
        from .api_cards import design_to_dict

        mode = body.get("mode", "point")
        if mode not in ("point", "envelope", "headline"):
            raise ValueError(f"unknown mode {mode!r}")
        if mode == "envelope" and getattr(design, "platform", "legacy") == "ingan_gan_nanowire":
            # N1: evaluate_envelope needs g2/Tj curves; refuse immediately instead of
            # failing with KeyError 'g2' after minutes of sampling.
            raise UnsupportedRun(NANOWIRE_ENVELOPE_UNSUPPORTED)
        T_grid = resolve_T_grid(body.get("T_grid"))
        if T_grid is None and getattr(design, "platform", "legacy") in NITRIDE_DISPLAY_GRID:
            # nitride evaluate() defaults to the single T_hs point; the Designer's "full curve"
            # asks for no grid, so give it an explicit display grid (a sampling choice, not physics).
            lo, hi, n = NITRIDE_DISPLAY_GRID[design.platform]
            T_grid = [float(x) for x in np.linspace(lo, hi, n)]
        ranged = None
        if mode == "envelope":
            raw = body.get("ranged")
            if raw is None:
                ranged = {k: list(v) for k, v in design_meta.default_ranged(design).items()}
            else:
                ranged = {str(k): [float(x) for x in v] for k, v in raw.items()}
        key = cache.make_key({"design": design_to_dict(design), "mode": mode,
                              "ranged": ranged, "T_grid": T_grid})
        run_id = make_run_id(card, key)
        cls = platform_class(headline_design(design) if mode == "headline" else design)
        n_T = len(T_grid) if T_grid is not None else DEFAULT_N_T

        if mode == "envelope":
            n_samples = 1
            for v in ranged.values():
                n_samples *= len(v)
            n_eval = n_samples + 1 + sum(n_samples // len(v) for v in ranged.values()) + 1
            eta = eval_seconds(cls, n_T) * (n_eval / self.max_workers + 2)
        else:
            eta = eval_seconds(cls, n_T)
            if cls in SLOW_CLASSES and n_T > 1:
                eta += eval_seconds(cls, 1)

        cached = cache.get(key)
        if cached is not None:
            job = self._register(Job(mode, 0.0, key, run_id))
            job.state = "done"
            job.result = cached["result"]
            job.events.append(("state", {"state": "done", "job_id": job.id, "cached": True}))
            job.events.append(("done", {"job_id": job.id, "result": job.result, "cached": True}))
            return {"job_id": job.id, "cached": True, "eta_s": 0.0, "run_id": run_id,
                    "result": job.result}

        record_base = {"key": key, "run_id": run_id, "card": card, "mode": mode,
                       "design": design_to_dict(design), "ranged": ranged, "T_grid": T_grid}

        if mode == "point" and cls == "legacy" and n_T <= INLINE_MAX_T:
            job = self._register(Job(mode, eta, key, run_id))
            from fsim_core.device import evaluate
            try:
                res = evaluate(design, T_grid=T_grid)
            except Exception as exc:  # noqa: BLE001 -- reported to the client
                err = classify_error(design, exc)
                job.fail(err)
                return {"job_id": job.id, "cached": False, "eta_s": eta, "run_id": run_id,
                        "error": err}
            result = point_result(design, res, mode=mode, run_id=run_id, card=card, T_grid=T_grid)
            cache.put(key, {**record_base, "result": result})
            job.finish(result)
            return {"job_id": job.id, "cached": False, "eta_s": eta, "run_id": run_id,
                    "result": result}

        job = self._register(Job(mode, eta, key, run_id))
        target = self._run_envelope if mode == "envelope" else self._run_point
        args = (job, design, card, T_grid, record_base) + ((ranged,) if mode == "envelope" else (mode, cls))
        threading.Thread(target=target, args=args, daemon=True).start()
        return {"job_id": job.id, "cached": False, "eta_s": eta, "run_id": run_id}

    def submit_call(self, fn_path: str, kwargs: dict, eta_s: float = 0.0, label: str = "scene") -> dict:
        """Run a scene job descriptor {fn: "fsim_studio.scene:...", kwargs} in
        the pool; cached by (fn, kwargs, fsim_core hash). The job result is
        whatever the function returns (a SceneSpec)."""
        if fn_path not in ALLOWED_CALLS:
            raise ValueError(f"job function {fn_path!r} not allowed")
        key = cache.make_key({"call": fn_path, "kwargs": kwargs})
        run_id = f"FS-{label}-{key[:6]}"
        cached = cache.get(key)
        if cached is not None:
            job = self._register(Job(label, 0.0, key, run_id))
            job.state = "done"
            job.result = cached["result"]
            job.events.append(("state", {"state": "done", "job_id": job.id, "cached": True}))
            job.events.append(("done", {"job_id": job.id, "result": job.result, "cached": True}))
            return {"job_id": job.id, "cached": True, "eta_s": 0.0, "run_id": run_id,
                    "result": job.result}
        if fn_path not in ALLOWED_CALLS:
            raise ValueError(f"job function {fn_path!r} not allowed")
        job = self._register(Job(label, float(eta_s or 0.0), key, run_id))

        def body():
            try:
                if not job.set_state("running"):
                    return
                fut = self.pool().submit(_w_call, fn_path, dict(kwargs))
                job.futures.append(fut)
                try:
                    result = fut.result()
                except cf.CancelledError:
                    return
                if job.cancelled:
                    return
                cache.put(key, {"key": key, "run_id": run_id, "mode": label, "call": fn_path,
                                "kwargs": kwargs, "result": result})
                job.finish(result)
            except cf.process.BrokenProcessPool as exc:
                self._reset_pool()
                job.fail({"kind": "pool", "message": f"worker pool broke: {exc}"})
            except Exception as exc:  # noqa: BLE001 -- reported to the client
                job.fail({"kind": "error", "message": f"{type(exc).__name__}: {exc}"})

        threading.Thread(target=body, daemon=True).start()
        return {"job_id": job.id, "cached": False, "eta_s": job.eta_s, "run_id": run_id}

    def cancel(self, job_id: str) -> dict | None:
        job = self.get(job_id)
        if job is None:
            return None
        with job.cond:
            if job.terminal:
                return job.snapshot()
            job.cancelled = True
            job.state = "cancelled"
            job.events.append(("state", {"state": "cancelled", "job_id": job.id}))
            job.cond.notify_all()
        for f in job.futures:
            f.cancel()
        return job.snapshot()

    # -- workers (threads in the server process; heavy work in the pool)
    def _run_point(self, job, design, card, T_grid, record_base, mode, cls):
        d_eval = headline_design(design) if mode == "headline" else design
        try:
            if not job.set_state("running"):
                return
            pool = self.pool()
            f_partial = None
            progressive = cls in SLOW_CLASSES and (T_grid is None or len(T_grid) > 1)
            if progressive:
                f_partial = pool.submit(_w_evaluate, d_eval, [float(d_eval.thermal.T_hs)])
                job.futures.append(f_partial)
            f_full = pool.submit(_w_evaluate, d_eval, T_grid)
            job.futures.append(f_full)
            if f_partial is not None:
                try:
                    pres = f_partial.result()
                except cf.CancelledError:
                    return
                except Exception:  # noqa: BLE001 -- the full run reports errors
                    pres = None
                if pres is not None and not job.cancelled and not f_full.done():
                    partial = point_result(design, pres, mode=mode, run_id=job.run_id, card=card,
                                           T_grid=[float(d_eval.thermal.T_hs)], partial=True)
                    job.partial = partial
                    job.emit("partial", {"job_id": job.id, "result": partial})
            try:
                res = f_full.result()
            except cf.CancelledError:
                return
            if job.cancelled:
                return
            result = point_result(design, res, mode=mode, run_id=job.run_id, card=card,
                                  T_grid=T_grid)
            cache.put(job.key, {**record_base, "result": result})
            job.finish(result)
        except cf.process.BrokenProcessPool as exc:
            self._reset_pool()
            job.fail({"kind": "pool", "message": f"worker pool broke: {exc}"})
        except Exception as exc:  # noqa: BLE001 -- reported to the client
            job.fail(classify_error(d_eval, exc))

    def _run_envelope(self, job, design, card, T_grid, record_base, ranged):
        from fsim_core.device import evaluate_envelope

        def progress(k, n):
            if job.cancelled:
                raise _Cancelled()
            job.progress = (k, n)
            job.emit("progress", {"job_id": job.id, "k": k, "n": n})

        try:
            if not job.set_state("running"):
                return
            pool = self.pool()
            f_mid = pool.submit(_w_evaluate, mid_design(design, ranged), T_grid)
            job.futures.append(f_mid)
            env = evaluate_envelope(design, {k: tuple(v) for k, v in ranged.items()},
                                    T_grid=T_grid, map_fn=pool.map, progress=progress)
            try:
                mid_res = f_mid.result()
            except ValueError:
                mid_res = None
            if job.cancelled:
                return
            result = envelope_result(design, env, mid_res, ranged=ranged, run_id=job.run_id,
                                     card=card, T_grid=T_grid)
            cache.put(job.key, {**record_base, "result": result})
            job.finish(result)
        except _Cancelled:
            return
        except cf.process.BrokenProcessPool as exc:
            self._reset_pool()
            job.fail({"kind": "pool", "message": f"worker pool broke: {exc}"})
        except Exception as exc:  # noqa: BLE001 -- reported to the client
            job.fail(classify_error(design, exc))

    # -- SSE
    def stream(self, job: Job, keepalive_s: float = 15.0):
        i = 0
        while True:
            with job.cond:
                if i >= len(job.events) and not job.terminal:
                    job.cond.wait(timeout=keepalive_s)
                new = job.events[i:]
                i += len(new)
                term = job.terminal and i >= len(job.events)
            if not new and not term:
                yield ": keepalive\n\n"
                continue
            for etype, data in new:
                yield f"event: {etype}\ndata: {dumps(data)}\n\n"
            if term:
                return

