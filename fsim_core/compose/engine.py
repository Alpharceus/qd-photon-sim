"""Composed-engine entry points (design 1.4): evaluate_spec / evaluate_composed / envelope / series /
estimate.

Chain per temperature: material.build -> setup.chain -> structure.response (gets the setup NA) ->
rates.combine -> drive.plan -> solvers.select(...).run -> observables. Only registry interfaces are used.

Ranged leaves: a point run uses the declared central value (lin midpoint / log geometric mean / first
choice) and sets scalars["range_note"] = "midpoint, not a prediction"; `envelope` reports
[min, max] over range corners + Latin-hypercube samples + the central point, never a mean.
"""
from __future__ import annotations

import math
from dataclasses import replace
from typing import Optional

import numpy as np

from .. import lindblad as lb
from .. import resonant_source as rs
from . import observables as obs_mod
from . import rates as rates_mod
from . import solvers as solvers_mod
from .params import ParamSpec, central, corners, lhs, resolve_leaf, widest_tag
from .registry import get as _get
from .solvers import SolverUnavailable
from .spec import AXIS_KEYS, T_K_SPEC, ComposedSpec, validate_spec

RANGE_NOTE = "midpoint, not a prediction"
HEADLINE = ("g2_0", "hom_I", "B_lens", "rate_detected_Hz")
CURVE_OBS = (("g2", "g2_0"), ("g2_raw", "g2_raw"), ("hom_I", "hom_I"), ("B_first_lens", "B_lens"),
             ("B_fibre", "B_fibre"), ("rate_detected_Hz", "rate_detected_Hz"), ("Gamma_ns", "Gamma_ns"),
             ("tau_ps", "tau_ps"), ("F_P_eff", "F_P_eff"), ("Z_phonon", "Z_phonon"), ("inversion", "inversion"))
SCALAR_OBS = (  # every registered chain value reported at the operating point (intrinsic convention)
    "lambda_nm", "F_P_eff", "beta", "eta_out", "eta_extract", "T_NA", "eta_pol", "eta_prep", "eta_rad",
    "eta_prep_x_eta_rad", "inversion", "Z_phonon", "Z_filt", "gamma_star_ns", "Gamma_ns", "n_ph", "tau_ps",
    "tau_purcell_ps", "linewidth_ueV", "detuning_ueV", "g2_0", "g2_raw", "hom_I", "hom_V_raw", "hom_V_corr",
    "B_first_lens", "B_lens", "B_fibre", "rate_detected_Hz")


# ------------------------------------------------------------------ spec handling
def _as_spec(spec) -> ComposedSpec:
    if isinstance(spec, dict):
        spec = ComposedSpec.from_dict(spec)
    res = validate_spec(spec)
    if not res["ok"]:
        lines = [f"{i['level']} {i['path']}: {i['message']}" for i in res["issues"]]
        raise ValueError("invalid composed spec:\n" + "\n".join(lines))
    return spec


def _entries(spec: ComposedSpec) -> dict:
    return {a: _get(a, getattr(spec, a)["kind"]) for a in AXIS_KEYS}


def _leaves(spec: ComposedSpec, entries: dict) -> dict:
    out = {"T_K": resolve_leaf(spec.T_K, T_K_SPEC, "T_K")}
    for a in AXIS_KEYS:
        blk = getattr(spec, a)
        for p in entries[a].params:
            path = f"{a}.{p.name}"
            if p.name in blk["params"]:
                out[path] = resolve_leaf(blk["params"][p.name], p, path)
            else:
                out[path] = replace(resolve_leaf(p.default, p, path), defaulted=True)
    return out


def _central_values(leaves: dict) -> dict:
    return {p: central(l) for p, l in leaves.items()}


# ------------------------------------------------------------------ one point
def _prepare(entries, values, T_K):
    reads = {a: set() for a in AXIS_KEYS}
    vals = {a: {p.name: values[f"{a}.{p.name}"] for p in entries[a].params} for a in AXIS_KEYS}

    def tr(a):
        return obs_mod.TrackedVals(vals[a], reads[a])

    em = entries["material"].build(tr("material"), T_K)
    setup = entries["setup"].chain(tr("setup"))
    sr = entries["structure"].response(tr("structure"), em, T_K, setup.NA)
    rt = rates_mod.combine(em, sr, T_K)
    plan = entries["drive"].plan(tr("drive"), em, rt)
    return em, setup, sr, rt, plan, vals, reads


def _run_point(entries, values, T_K, solver_req):
    em, setup, sr, rt, plan, vals, reads = _prepare(entries, values, T_K)
    solver, issues = solvers_mod.select(plan, em, rt, solver_req, want_hom=True)
    if solver is None:
        raise SolverUnavailable("; ".join(issues))
    res = solver.run(plan, em, rt)
    ctx = obs_mod.Context(em, sr, rt, plan, res, setup, T_K, vals, reads)
    ctx.solver = solver
    ctx.select_issues = issues
    return ctx


def _validity(ctx) -> tuple:
    rt, res = ctx.rt, ctx.res
    reasons = list(res.reasons) if not res.valid else []
    valid = bool(res.valid)
    if rt.regime == "strong" and ctx.solver.name != "jc":
        valid = False
        reasons.append("strong coupling: weak-coupling Purcell model is not valid; use solver 'jc'")
    return valid, reasons


def _notes(ctx) -> list:
    n = list(ctx.rt.notes) + list(ctx.sr.notes) + list(ctx.select_issues)
    n += [str(x) for x in ctx.res.extra.get("notes", [])]
    n += [r for r in ctx.res.reasons]
    n.append("phonons (tier P0): HOM factor Z^2 unfiltered, efficiency x Z ZPL-filtered, never both")
    return n


def _series(ctx) -> dict:
    out = {}
    for k, v in ctx.res.extra.get("series", {}).items():
        out[k] = v
    plan, em, rt = ctx.plan, ctx.em, ctx.rt
    if plan.timing == "pulsed" and math.isfinite(ctx.res.I_single) and ctx.solver.name != "jc":
        sig = float(lb.mev_to_rate(em.sigma_sd_ueV * 1e-3))
        dts = np.logspace(-1, 3, 25)
        y = [rs.hom_vs_separation(float(dt), float(rt.Gamma_ns),
                                  float(rt.gamma_star_ns), sig, float(em.tau_sd_ns)) for dt in dts]
        out["hom_vs_dt"] = {"x": [float(v) for v in dts], "y": [float(v) for v in y], "x_unit": "ns",
                            "note": "spectral diffusion (sigma_sd, tau_sd) between photons dt apart; "
                                    "single-emitter overlap Gamma/(Gamma+gamma*) at dt -> 0"}
    return out


# ------------------------------------------------------------------ evaluate
def evaluate_spec(spec, T_grid=None) -> dict:
    """spec: ComposedSpec or its dict. Returns {'curves', 'scalars', 'series'}."""
    spec = _as_spec(spec)
    entries = _entries(spec)
    leaves = _leaves(spec, entries)
    values = _central_values(leaves)
    T_op = float(values["T_K"])
    grid = sorted({float(t) for t in (T_grid if T_grid is not None else [T_op])} | {T_op})
    ctxs = {T: _run_point(entries, {**values, "T_K": T}, T, spec.solver) for T in grid}
    ctx = ctxs[T_op]
    Ts = [float(t) for t in (T_grid if T_grid is not None else [T_op])]

    def obs_val(c, oid):
        return obs_mod.get_observable(oid).compute(c, "intrinsic")

    curves = {"T_hs": Ts, "Tj": list(Ts)}
    for name, oid in CURVE_OBS:
        curves[name] = [obs_val(ctxs[float(t)], oid) for t in Ts]

    reads = {a: set() for a in AXIS_KEYS}
    for c in ctxs.values():
        for a in AXIS_KEYS:
            reads[a] |= c.reads[a]
    scalars: dict = {}
    prov = {}
    for oid in SCALAR_OBS:
        o = obs_mod.get_observable(oid)
        v = o.compute(ctx, "intrinsic")
        scalars[oid] = v
        tag, inputs = obs_mod.tag_of(o, leaves, reads)
        prov[oid] = {"tag": tag, "inputs": inputs}
        if not math.isfinite(v):
            r = o.invalid_reason(ctx)
            if r:
                scalars[f"{oid}_invalid_reason"] = r
    for oid, key in (("hom_V", "hom_V_raw"),):
        scalars[oid] = scalars[key]
    valid, reasons = _validity(ctx)
    scalars.update({
        "solver": ctx.solver.name, "regime": ctx.rt.regime, "valid": valid, "invalid_reasons": reasons,
        "T_K": T_op, "brightness_convention": obs_mod.BRIGHTNESS_CONVENTION,
        "gamma_star_convention": rates_mod.GAMMA_STAR_CONVENTION, "notes": _notes(ctx),
        "tag_chain": widest_tag("DR", *[prov[h]["tag"] for h in HEADLINE]),
        "provenance": prov, "laser_leak_used": float(ctx.d["laser_leak_used"]),
        "jc_recommended": bool(ctx.rt.jc_recommended), "L": float(ctx.rt.L),
    })
    if "tau_jc_ps" in ctx.res.extra:
        scalars["tau_jc_ps"] = float(ctx.res.extra["tau_jc_ps"])
    for k in ("coherent_fraction", "flux_ns", "tpe_pi_area_rad"):
        if k in ctx.res.extra:
            scalars[k] = float(ctx.res.extra[k])
    if any(l.is_free for l in leaves.values()):
        scalars["range_note"] = RANGE_NOTE
    return {"curves": curves, "scalars": scalars, "series": _series(ctx)}


def evaluate_composed(design, T_grid=None) -> dict:
    """device.evaluate() arm for platform == 'composed'."""
    return evaluate_spec(ComposedSpec.from_dict(design.composed), T_grid)


# ------------------------------------------------------------------ envelope / series / estimate
def _value_of(ctx, request):
    oid, conv = obs_mod.parse_request(request)
    return obs_mod.get_observable(oid).compute(ctx, conv)


def envelope(spec, observables, max_corners: int = 16, n_lhs: int = 32) -> dict:
    """{request: {"range": [min, max], "central": c, "tag": t, "inputs": [...], "n_points": n}}, plus
    "_meta". min/max over range corners, LHS samples and the central point (so min <= central <= max);
    a mean is never computed. A request is 'id' or 'id:convention'."""
    spec = _as_spec(spec)
    entries = _entries(spec)
    leaves = _leaves(spec, entries)
    base = _central_values(leaves)
    pts = [dict(base)]
    if any(l.kind == "range" for l in leaves.values()):
        pts += [{**base, **p} for p in corners(leaves, max_corners)]
        pts += [{**base, **p} for p in lhs(leaves, n_lhs)]
    elif any(l.kind == "choice" for l in leaves.values()):
        pts += [{**base, **p} for p in lhs(leaves, n_lhs)]
    vals = {r: [] for r in observables}
    central_val, nfail, reads = {}, 0, None
    for i, p in enumerate(pts):
        try:
            c = _run_point(entries, p, float(p["T_K"]), spec.solver)
        except SolverUnavailable:
            if i == 0:
                raise
            nfail += 1
            continue
        if reads is None:
            reads = c.reads
        for r in observables:
            v = _value_of(c, r)
            if i == 0:
                central_val[r] = v
            if math.isfinite(v):
                vals[r].append(v)
    out = {}
    for r in observables:
        oid, _ = obs_mod.parse_request(r)
        o = obs_mod.get_observable(oid)
        tag, inputs = obs_mod.tag_of(o, leaves, reads)
        v = vals[r]
        out[r] = {"range": [min(v), max(v)] if v else [float("nan")] * 2, "central": central_val[r],
                  "tag": tag, "inputs": inputs, "n_points": len(v)}
    out["_meta"] = {"n_points": len(pts), "n_failed": nfail, "range_note": RANGE_NOTE}
    return out


def series(spec, path: str, values, observables) -> dict:
    """Sweep one leaf path (e.g. 'structure.F_P', or 'T_K') over `values` at the central point of every
    other leaf. -> {"path", "x", "y": {request: [..]}, "tag": {request: tag}, "failed": [indices]}."""
    spec = _as_spec(spec)
    entries = _entries(spec)
    leaves = _leaves(spec, entries)
    if path not in leaves:
        raise KeyError(f"unknown leaf path {path!r}; known: {sorted(leaves)}")
    base = _central_values(leaves)
    y = {r: [] for r in observables}
    failed, reads = [], None
    for i, x in enumerate(values):
        p = {**base, path: x}
        try:
            c = _run_point(entries, p, float(p["T_K"]), spec.solver)
        except SolverUnavailable:
            failed.append(i)
            for r in observables:
                y[r].append(float("nan"))
            continue
        reads = reads or c.reads
        for r in observables:
            y[r].append(_value_of(c, r))
    tags = {}
    for r in observables:
        o = obs_mod.get_observable(obs_mod.parse_request(r)[0])
        tags[r] = obs_mod.tag_of(o, leaves, reads or {a: set() for a in AXIS_KEYS})[0]
    return {"path": path, "unit": leaves[path].unit, "x": [float(v) for v in values],
            "y": {r: np.asarray(v, dtype=float) for r, v in y.items()}, "tag": tags, "failed": failed}


def estimate(spec, n_points: int = 1) -> float:
    """Seconds for n_points point runs: the selected solver's eta_seconds at the central point."""
    spec = _as_spec(spec)
    entries = _entries(spec)
    values = _central_values(_leaves(spec, entries))
    em, setup, sr, rt, plan, _v, _r = _prepare(entries, values, float(values["T_K"]))
    solver, issues = solvers_mod.select(plan, em, rt, spec.solver, want_hom=True)
    if solver is None:
        raise SolverUnavailable("; ".join(issues))
    return float(solver.eta_seconds(plan, n_points))
