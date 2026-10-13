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
    "lambda_nm", "F_P_eff", "F_P_peak", "beta", "eta_out", "eta_extract", "T_NA", "eta_pol", "eta_prep", "eta_rad",
    "eta_prep_x_eta_rad", "inversion", "Z_phonon", "Z_filt", "gamma_star_ns", "Gamma_ns", "n_ph", "tau_ps",
    "tau_1e_ps", "tau_purcell_ps", "linewidth_ueV", "detuning_ueV", "g2_0", "g2_raw", "hom_I", "hom_V_raw", "hom_V_corr",
    "B_first_lens", "B_lens", "B_fibre", "rate_detected_Hz", "tau_ratio", "s_0", "s_sat", "coherent_fraction")


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
def combine(em, sr, T_K, pump_ns: float = 0.0):
    """rates.combine plus the cw incoherent-pump hook (R2 FIX-1a, Fable-A PHA-1 deferred note): with pump_ns = P = r (1 + p) > 0 the
    pump broadens the dot line the cavity filter sees, so L, Gamma and beta are the PUMPED values (rates.apply_pump, ch08 eq.
    dp-lor); pump_ns = 0 is rates.combine bit for bit. The planner needs the unpumped Gamma to set r, so the engine calls
    combine(pump_ns=0) for the plan and again with the plan's P for the solver and the observables."""
    rt = rates_mod.combine(em, sr, T_K)
    return rates_mod.apply_pump(rt, float(pump_ns)) if float(pump_ns) > 0.0 else rt


def _cw_pump_ns(plan, rt) -> float:
    """P = r (1 + p_xx) of a cw incoherent / electrical plan (0 otherwise); same definition the rate solver uses."""
    if plan.timing != "cw" or plan.family not in ("incoherent", "electrical"):
        return 0.0
    r, _ = solvers_mod._interp_pump(plan, rt)
    p_xx = float(plan.inputs.get("pump_ratio_xx", 1.0)) if "cascade" in plan.needs else 0.0
    return float(r) * (1.0 + p_xx)


def _prepare(entries, values, T_K):
    reads = {a: set() for a in AXIS_KEYS}
    vals = {a: {p.name: values[f"{a}.{p.name}"] for p in entries[a].params} for a in AXIS_KEYS}

    def tr(a):
        return obs_mod.TrackedVals(vals[a], reads[a])

    em = entries["material"].build(tr("material"), T_K)
    setup = entries["setup"].chain(tr("setup"))
    sr = entries["structure"].response(tr("structure"), em, T_K, setup.NA)
    rt = combine(em, sr, T_K)
    plan = entries["drive"].plan(tr("drive"), em, rt)
    if "T_K" not in plan.inputs:                  # P9: phonon-aware solvers need the bath temperature
        plan = replace(plan, inputs={**plan.inputs, "T_K": float(T_K)})
    return em, setup, sr, rt, plan, vals, reads


def _run_point(entries, values, T_K, solver_req):
    em, setup, sr, rt, plan, vals, reads = _prepare(entries, values, T_K)
    solver, issues = solvers_mod.select(plan, em, rt, solver_req, want_hom=True)
    if solver is None:
        raise SolverUnavailable("; ".join(issues))
    P = 0.0 if solver.name == "jc" else _cw_pump_ns(plan, rt)      # the jc solver puts the pump in its Lindbladian
    if P > 0.0:
        rt = combine(em, sr, T_K, pump_ns=P)                       # observables and solver see the pumped Gamma / beta
    res = solver.run(plan, em, rt)
    ctx = obs_mod.Context(em, sr, rt, plan, res, setup, T_K, vals, reads)
    ctx.solver = solver
    ctx.structure_kind = entries["structure"].kind
    ctx.select_issues = issues
    return ctx


JC_REMEDY_UNSUPPORTED = ("not supported in the strong-coupling regime: weaken the coupling (lower Q or larger mode volume V "
                         "until g < (kappa + Gamma_bg + k_nr + gamma*)/4) or choose a drive the Jaynes-Cummings solver supports "
                         "(resonant_pulse, cw_resonant, above_band, p_shell, electrical_pulsed, electrical_cw)")


def _strong_coupling_reason(ctx) -> str:
    """Why a strong-coupling run on a non-jc solver is invalid. Never recommends the jc solver when jc refuses the plan
    (R2 FIX-1a, audit 4): it asks the jc solver itself and states its refusal plus a UI-facing suggestion."""
    jc = _get("solver", "jc")
    ok, why = jc.supports(ctx.plan, ctx.em, ctx.rt)
    if ok:
        return ("strong coupling: the weak-coupling Purcell model is not valid for this solver; use solver 'jc' "
                "(auto selects it)")
    return ("strong coupling: the weak-coupling Purcell model is not valid and the Jaynes-Cummings solver cannot run this "
            f"drive ({why}); " + JC_REMEDY_UNSUPPORTED)


def _validity(ctx) -> tuple:
    rt, res = ctx.rt, ctx.res
    reasons = list(res.reasons) if not res.valid else []
    valid = bool(res.valid)
    try:
        from . import structures as _st
        why = _st.incompatibility(getattr(ctx, "structure_kind", ""), ctx.em, ctx.vals.get("structure", {}))
    except Exception:                                   # noqa: BLE001
        why = None
    if why:
        valid = False
        reasons.append(f"incompatible material / structure pair: {why}")
    if rt.regime == "strong" and ctx.solver.name != "jc":
        valid = False
        reasons.append(_strong_coupling_reason(ctx))
    return valid, reasons


def _notes(ctx) -> list:
    n = list(ctx.rt.notes) + list(ctx.sr.notes) + list(ctx.select_issues)
    n += [str(x) for x in ctx.res.extra.get("notes", [])]
    n += [r for r in ctx.res.reasons]
    n += list(ctx.d.get("_warnings", []))
    n.append("phonons (tier P0): HOM factor Z^2 unfiltered, efficiency x Z ZPL-filtered, never both")
    ep = ctx.d.get("eta_prep", float("nan"))
    if ctx.plan.timing == "pulsed" and math.isfinite(ep) and ep < 0.05:
        n.append(f"preparation probability {ep:.3g} < 5 %: the drive barely excites the emitter at these inputs (off resonance "
                 "or wrong pulse area); g2 and HOM of such a run are a symptom of the residual re-excitation, not a result")
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


def _series_unavailable(ctx) -> dict:
    """{series id: reason} of an advertised series the run does not produce on purpose (R2 N1, a named reason instead of a
    silently missing curve). hom_vs_dt is the Gaussian average over the relative spectral-diffusion detuning of two photons dt
    apart, which needs the overlap of two photons with a DETUNED second emitter. The jc solver has only the two-time overlap of
    two identical photons (lb.indistinguishability), so it does not compute the series; the scalar hom_I is still the jc value."""
    if ctx.plan.timing != "pulsed" and "hom_vs_dt" not in ctx.res.extra.get("series", {}):
        # R2 F6: the cw variants of above_band / p_shell no longer advertise it; a request is answered with this reason
        return {"hom_vs_dt": "a cw stream has no per-pulse HOM overlap, so there is no overlap-against-separation series"}
    if ctx.plan.timing == "pulsed" and ctx.solver.name == "jc" and "hom_vs_dt" not in ctx.res.extra.get("series", {}):
        return {"hom_vs_dt": "not computed with the Jaynes-Cummings solver: the separation series needs the overlap of photons with a "
                             "detuned second emitter, the jc two-time machinery has only identical photons (hom_I is the jc value "
                             "at dt -> 0)"}
    return {}


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
        tag, inputs = obs_mod.tag_of(o, leaves, reads, ctx.sr.tag_floor)
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
    for k in ("P_X_end", "tau_XX_ps", "g2_XX", "hom_XX", "n_XX"):      # FX2b: only when the solver reports them
        if k in ctx.res.extra:
            scalars[k] = float(ctx.d[k])
    scalars["captions"] = captions(ctx)
    # R2 F56: an observable the drive advertises (beyond SCALAR_OBS) that is NaN at the operating point comes with its reason
    # (hom_I_gated with no HOM gate set): never a bare NaN
    adv = entries["drive"]
    adv = adv.observables_for({k.split(".", 1)[1]: v for k, v in values.items() if k.startswith("drive.")})         if hasattr(adv, "observables_for") else adv.observables
    gone_o = {}
    for oid in sorted(set(adv) - set(SCALAR_OBS)):
        try:
            o = obs_mod.get_observable(oid)
            if not math.isfinite(o.compute(ctx, o.conventions[0])):
                gone_o[oid] = str(o.invalid_reason(ctx, o.conventions[0]) or f"{oid} is not computed for this spec")
        except Exception:  # noqa: BLE001 -- a series id (g2_tau, ...) is not a registered scalar
            continue
    if gone_o:
        scalars["observables_unavailable"] = gone_o
    gone = _series_unavailable(ctx)
    if gone:
        scalars["series_unavailable"] = gone
    if any(l.is_free for l in leaves.values()):
        scalars["range_note"] = RANGE_NOTE
    return {"curves": curves, "scalars": scalars, "series": _series(ctx)}


def captions(ctx) -> dict:
    """Formula captions that match what was computed (audit 17): the lifetime formula depends on the solver and the
    detected-rate formula on pulsed vs cw timing. The front end shows these instead of hard-coded text."""
    tau = ("tau = slowest decay time of the one-excitation Liouvillian sector of the dot + cavity system (JC solver); "
           "it is NOT 1/(Gamma + k_nr) in the strong-coupling regime" if ctx.solver.name == "jc"
           else "tau = 1/(Gamma + k_nr), Gamma = Gamma0 (F_bg + F_P L) (Purcell lifetime)")
    if ctx.plan.timing == "pulsed":
        rate = "Detected rate = f_rep . B_fibre . eta_optics . eta_det"
    else:
        rate = ("Detected rate (cw, no f_rep) = photon flux (solver steady state) . beta . eta_out . eta_pol . Z_filt . eta_fibre "
                ". eta_optics . eta_det")
    return {"tau_ps": tau, "rate_detected_Hz": rate}


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
    central_val, nfail, ninvalid, floor = {}, 0, 0, ""
    reads = {a: set() for a in AXIS_KEYS}
    central_valid = True
    for i, p in enumerate(pts):
        try:
            c = _run_point(entries, p, float(p["T_K"]), spec.solver)
        except SolverUnavailable:
            if i == 0:
                raise
            nfail += 1
            continue
        for a in AXIS_KEYS:                       # union of the leaves read at ANY point (a point may take another branch)
            reads[a] |= c.reads[a]
        if c.sr.tag_floor:
            floor = obs_mod.widest_tag(floor, c.sr.tag_floor) if floor else c.sr.tag_floor
        valid, _why = _validity(c)
        if not valid:
            if i == 0:
                central_valid = False             # the central point is always reported (min <= central <= max)
            else:
                ninvalid += 1                     # other invalid points are left out of the range
                continue
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
        tag, inputs = obs_mod.tag_of(o, leaves, reads, floor)
        v = vals[r]
        out[r] = {"range": [min(v), max(v)] if v else [float("nan")] * 2, "central": central_val[r],
                  "tag": tag, "inputs": inputs, "n_points": len(v)}
    out["_meta"] = {"n_points": len(pts), "n_failed": nfail, "n_invalid": ninvalid,
                   "central_valid": central_valid, "range_note": RANGE_NOTE}
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
    failed, floor = [], ""
    last_ctx = None
    reads = {a: set() for a in AXIS_KEYS}
    for i, x in enumerate(values):
        p = {**base, path: x}
        try:
            c = _run_point(entries, p, float(p["T_K"]), spec.solver)
        except SolverUnavailable:
            failed.append(i)
            for r in observables:
                y[r].append(float("nan"))
            continue
        for a in AXIS_KEYS:
            reads[a] |= c.reads[a]
        if c.sr.tag_floor:
            floor = obs_mod.widest_tag(floor, c.sr.tag_floor) if floor else c.sr.tag_floor
        last_ctx = c
        for r in observables:
            y[r].append(_value_of(c, r))
    tags = {}
    for r in observables:
        o = obs_mod.get_observable(obs_mod.parse_request(r)[0])
        tags[r] = obs_mod.tag_of(o, leaves, reads, floor)[0]
    out = {"path": path, "unit": leaves[path].unit, "x": [float(v) for v in values],
           "y": {r: [float(t) for t in v] for r, v in y.items()}, "tag": tags, "failed": failed}
    # R2 F56: a request that is NaN at every point comes back WITH the reason the engine states for it (never a bare NaN)
    gone = {}
    if last_ctx is not None:
        for r in observables:
            if out["y"][r] and not any(math.isfinite(t) for t in out["y"][r]):
                oid, conv = obs_mod.parse_request(r)
                reason = obs_mod.get_observable(oid).invalid_reason(last_ctx, conv)
                gone[r] = str(reason) if reason else f"{r} is not computed for this spec (NaN at every point)"
    if gone:
        out["unavailable"] = gone
    return out


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
