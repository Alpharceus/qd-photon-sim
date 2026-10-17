"""Registry-driven coverage of the composable engine: material x structure x drive x solver (R2 test audit 3, 4, 5, 6, 10, 11, 12).

Kinds are ENUMERATED from fsim_core.compose.registry, never from a hand list. verify/data/compose_compat.yaml only DECLARES
what is supported (and why a combination is refused). Slices (a full cube would take hours):
  S x D   every structure x every drive (+ drive variants), reference material (per structure)
  M x S   every material x every structure, resonant_pulse
  M x D   every material x every drive (+ variants), bulk_planar
  Sv x D  every solver x every drive (declared accepted or cleanly refused), auto-solver choice
For every run: outcome as declared; physical bounds; every advertised observable of the drive finite (R2 6);
a "use solver 'X'" remedy in an invalid reason must really run (R2 4); a drive does what its summary claims (R2 11);
the Composer's observable requests are valid (R2 5); every ParamSpec moves an output (R2 10, dead controls).

KNOWN_OPEN lists R2-audit findings that are being fixed by other agents. A failure that matches an entry prints
KNOWN-OPEN and does not fail the exit code; an entry that matches nothing prints "FIXED? delete the entry". Everything
else must pass.

Run: python verify/verify_compose_matrix.py [--workers N] [--quick]   (exit 0 iff every non-KNOWN-OPEN check passes)
"""
from __future__ import annotations

import fnmatch
import math
import multiprocessing as mp
import os
import sys
import time
import warnings
from pathlib import Path

for _v in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS"):
    os.environ.setdefault(_v, "2")            # user thermal rule 2026-10-12: few threads per process

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")

import yaml                                                                    # noqa: E402

COMPAT = ROOT / "verify" / "data" / "compose_compat.yaml"

# ---------------------------------------------------------------------------------------------------------------------
# KNOWN_OPEN: {id, finding, desc, match}. `match` maps a key of the failing check (slice, check, material, structure,
# drive, solver, param, request) to a value or a list of values; all keys must match. Delete an entry when its fix lands.
KNOWN_OPEN = []


# ---------------------------------------------------------------------------------------------------------------------
# worker side (top-level functions: picklable under the spawn start method)
def _spec(m, s, d, solver="auto", over=None, T=4.2):
    spec = {"schema": "fsim-composed/1", "T_K": T,
            "material": {"kind": m, "params": {}}, "structure": {"kind": s, "params": {}},
            "drive": {"kind": d, "params": {}}, "setup": {"kind": "lens_fibre_det", "params": {}}, "solver": solver}
    for path, v in (over or {}).items():
        if path == "T_K":
            spec["T_K"] = v
        else:
            axis, name = path.split(".", 1)
            spec[axis]["params"][name] = v
    return spec


def _fl(v):
    try:
        return float(v)
    except Exception:                                                          # noqa: BLE001
        return float("nan")


def run_point(job):
    """job: {m, s, d, solver, over, T, series} -> a plain dict (picklable) describing the run."""
    from fsim_core.compose import engine, observables as obs_mod, registry as R
    from fsim_core.compose.solvers import SolverUnavailable
    from fsim_core.compose.spec import validate_spec
    out = {"job": job, "outcome": None}
    spec = _spec(job["m"], job["s"], job["d"], job.get("solver", "auto"), job.get("over"), job.get("T", 4.2))
    t0 = time.time()
    try:
        res = validate_spec(spec)
        if not res["ok"]:
            out.update(outcome="refused", detail="; ".join(f"{i['path']}: {i['message']}" for i in res["issues"] if i["level"] == "error"))
            return out
        sp = engine._as_spec(spec)
        entries = engine._entries(sp)
        leaves = engine._leaves(sp, entries)
        values = engine._central_values(leaves)
        ctx = engine._run_point(entries, values, float(values["T_K"]), sp.solver)
    except SolverUnavailable as e:
        out.update(outcome="unavailable", detail=str(e))
        return out
    except Exception as e:                                                     # noqa: BLE001
        out.update(outcome="crash", detail=f"{type(e).__name__}: {str(e)[:300]}")
        return out
    valid, reasons = engine._validity(ctx)
    obs, obs_reason = {}, {}
    for oid in R.kinds("observable"):
        o = obs_mod.get_observable(oid)
        for conv in o.conventions:
            try:
                obs[f"{oid}:{conv}"] = _fl(o.compute(ctx, conv))
                if obs[f"{oid}:{conv}"] != obs[f"{oid}:{conv}"]:               # NaN: the reason the engine states for it (R2 F6)
                    obs_reason[f"{oid}:{conv}"] = str(o.invalid_reason(ctx, conv) or "")
            except Exception as e:                                             # noqa: BLE001
                obs[f"{oid}:{conv}"] = f"ERR {type(e).__name__}: {str(e)[:120]}"
    series = {}
    try:
        for k, v in engine._series(ctx).items():
            series[k] = [_fl(y) for y in v.get("y", [])]
    except Exception as e:                                                     # noqa: BLE001
        series["__error__"] = [str(e)[:120]]
    try:
        notes = [str(n) for n in engine._notes(ctx)]
    except Exception:                                                          # noqa: BLE001
        notes = []
    try:
        gone = {k: str(v) for k, v in engine._series_unavailable(ctx).items()}   # R2 N1: an advertised series left out WITH a reason
    except Exception:                                                          # noqa: BLE001
        gone = {}
    plan = ctx.plan
    out.update(outcome="ok" if valid else "invalid", valid=bool(valid), reasons=[str(r) for r in reasons],
               solver=ctx.solver.name, regime=str(ctx.rt.regime), family=str(plan.family), timing=str(plan.timing),
               obs=obs, obs_reason=obs_reason, series=series, series_unavailable=gone, notes=notes, seconds=time.time() - t0)
    return out


def _vec(r):
    if r["outcome"] not in ("ok", "invalid"):
        return None
    v = dict(r["obs"])
    v["__valid"], v["__solver"], v["__regime"] = r["valid"], r["solver"], r["regime"]
    for k, ys in r["series"].items():
        v[f"series:{k}"] = tuple(ys)
    return v


def _same(a, b):
    if a.keys() != b.keys():
        return False
    for k in a:
        x, y = a[k], b[k]
        if isinstance(x, float) and isinstance(y, float):
            if math.isnan(x) and math.isnan(y):
                continue
            if not math.isclose(x, y, rel_tol=1e-9, abs_tol=1e-15):
                return False
        elif isinstance(x, tuple) and isinstance(y, tuple):
            if len(x) != len(y) or not all((math.isnan(p) and math.isnan(q)) or math.isclose(p, q, rel_tol=1e-9, abs_tol=1e-15)
                                           for p, q in zip(x, y)):
                return False
        elif x != y:
            return False
    return True


def probe_entry(task):
    """Dead-control probe of every param of ONE registry entry. Each param is moved on the first chain; only a param that
    moves nothing there escalates to the next chain (other drive / structure / temperature), so a parameter that acts only
    for a phonon solver, a TPE cascade or a hot dot is still found."""
    from fsim_core.compose import registry as R
    axis, kind, chains, context = task["axis"], task["kind"], task["chains"], task["context"]
    extra_cands = task.get("candidates", {})
    entry = R.get(axis, kind)
    res = {}
    base_cache = {}

    def base_of(chain, over):
        k = (chain, tuple(sorted(over.items(), key=lambda kv: kv[0])))
        if k not in base_cache:
            m, s, d = chain
            r = run_point({"m": m, "s": s, "d": d, "over": over})
            base_cache[k] = (_vec(r), f"{r['outcome']}: {r.get('detail', '')[:100]}")
        return base_cache[k]

    for p in entry.params:
        path = f"{axis}.{p.name}"
        notes, verdict = [], None
        for chain, chain_over in chains:
            over = dict(chain_over)
            over.update(context.get(f"{axis}.{kind}.{p.name}", context.get(path, {})))      # a per-kind key wins over the axis-wide one
            if p.visible_if:
                for gname, allowed in p.visible_if.items():
                    allowed = [allowed] if isinstance(allowed, (str, bool, int, float)) else list(allowed or [])
                    if any(q.name == gname for q in entry.params) and allowed:
                        over[f"{axis}.{gname}"] = allowed[0]
            base, why = base_of(chain, over)
            if base is None:
                notes.append(f"chain {'|'.join(chain)}: base {why}")
                continue
            dflt = over.get(path, p.default)
            if isinstance(p.default, bool):
                cands = [not dflt]
            elif p.choices:
                cands = [c for c in p.choices if c != dflt]
            elif isinstance(p.default, str):
                verdict = ("skip", "free-text parameter")
                break
            else:
                lo, hi = p.lo, p.hi
                if lo is not None and hi is not None and hi > lo:
                    if p.scale == "log" and lo > 0:
                        g = lambda f: math.exp(math.log(lo) + f * (math.log(hi) - math.log(lo)))   # noqa: E731
                    else:
                        g = lambda f: lo + f * (hi - lo)                                            # noqa: E731
                    cands = [g(0.75), g(0.25), hi, lo]
                else:
                    x = float(dflt)
                    cands = [x * 2.0 if x else 1.0, x * 0.5 if x else -1.0, x + 1.0]
                if isinstance(p.default, int):
                    cands = [int(round(c)) for c in cands]
            cands = [c for c in extra_cands.get(path, []) + cands if c != dflt]      # declared regimes first (R2 F10b)
            if not cands:
                verdict = ("skip", "a single choice: nothing to move")
                break
            n_ok, refused = 0, []
            for c in cands:
                m, s, d = chain
                r = run_point({"m": m, "s": s, "d": d, "over": {**over, path: c}})
                v = _vec(r)
                if v is None:
                    refused.append((c, r["outcome"], r.get("detail", "")[:90]))
                    continue
                n_ok += 1
                if not _same(base, v):
                    verdict = ("ok", f"changes an output at {c!r} on {'|'.join(chain)}")
                    break
            if verdict:
                break
            if not n_ok and refused and all(o in ("unavailable", "refused") for _, o, _ in refused):
                verdict = ("ok", f"every other value is refused with a named reason ({refused[0][2]})")
                break
            notes.append(f"{'|'.join(chain)}: " + (f"no change for {cands!r}" if n_ok else f"{refused[:1]}"))
        if verdict is None:
            crashes = [n for n in notes if "crash" in n]
            verdict = ("error" if crashes and len(crashes) == len(notes) else "dead", "; ".join(notes)[:420])
        res[p.name] = verdict
    return {"axis": axis, "kind": kind, "chains": [c for c, _ in chains], "result": res}


# ---------------------------------------------------------------------------------------------------------------------
class Report:
    def __init__(self):
        self.lines = []            # (ok, known_id, text)
        self.hits = {}
        self.notes = []

    def add(self, name, failures, key_base):
        """failures: list of (check, extra_key_dict, detail). One reported line per name."""
        if not failures:
            self.lines.append((True, None, f"[PASS] {name}"))
            return
        kn, bad = [], []
        for chk, extra, detail in failures:
            key = {**key_base, **extra, "check": chk}
            ent = known(key)
            if ent:
                self.hits[ent["id"]] = self.hits.get(ent["id"], 0) + 1
                kn.append((ent["id"], chk, detail))
            else:
                bad.append((chk, detail))
        if bad:
            self.lines.append((False, None, f"[FAIL] {name}  ({'; '.join(f'{c}: {d}' for c, d in bad)[:600]})"))
        else:
            ids = sorted({k[0] for k in kn})
            self.lines.append((True, ids, f"[KNOWN-OPEN {','.join(ids)}] {name}  ({'; '.join(f'{c}: {d}' for _, c, d in kn)[:300]})"))


def known(key):
    for ent in KNOWN_OPEN:
        ok = True
        for k, want in ent["match"].items():
            have = key.get(k)
            if isinstance(want, str) and want.startswith("re:"):
                import re as _re
                if not (isinstance(have, str) and _re.fullmatch(want[3:], have)):
                    ok = False
                    break
            elif isinstance(want, (list, tuple, set)):
                if have not in want:
                    ok = False
                    break
            elif have != want:
                ok = False
                break
        if ok:
            return ent
    return None


def match_axes(rule, key):
    for ax in ("material", "structure", "drive"):
        if ax in rule:
            w = rule[ax]
            w = w if isinstance(w, list) else [w]
            if key.get(ax) not in w:
                return False
    return True


def excluded_reason(table, m, s, d):
    key = {"material": m, "structure": s, "drive": d}
    for rule in table:
        if match_axes(rule, key):
            return rule["reason"]
    return None


def is_num(x):
    return isinstance(x, float) and math.isfinite(x)


def val(r, oid, conv=None):
    """Observable value from a run result: first declared convention unless `conv`."""
    for k, v in r["obs"].items():
        i, c = k.split(":", 1)
        if i == oid and (conv is None or c == conv):
            return v
    return float("nan")


def first_conv_key(r, oid):
    for k in r["obs"]:
        if k.split(":", 1)[0] == oid:
            return k
    return None


def combo_checks(r, ctx):
    """Physics checks of one run: list of (check, extra, detail) failures."""
    fails = []
    drive_obs = ctx["drive_obs"]
    g2 = val(r, "g2_0")
    if not is_num(g2) or g2 < 0.0:
        fails.append(("bounds", {}, f"g2_0 = {g2!r} (need finite >= 0)"))
    elif g2 > 1.0 + 1e-9:
        fails.append(("bounds", {}, f"g2_0 = {g2:.4g} > 1 for a single emitter"))
    tau = val(r, "tau_ps")
    if not (is_num(tau) and tau > 0.0):
        fails.append(("bounds", {}, f"tau_ps = {tau!r} (need > 0)"))
    for oid in ("B_first_lens", "B_fibre"):
        b = val(r, oid)                                    # NaN is the 'advertised' check's business
        if is_num(b) and not (-1e-12 <= b <= 1.0 + 1e-9):
            fails.append(("bounds", {}, f"{oid} = {b!r} outside [0, 1]"))
    h = val(r, "hom_I")
    if is_num(h) and not (-1e-12 <= h <= 1.0 + 1e-9):
        fails.append(("bounds", {}, f"hom_I = {h!r} outside [0, 1]"))
    # The raw visibility of the engine is V = [(1 - g2) I - g2](1 - eps): it is >= -g2 by construction (an indistinguishable-
    # photon overlap I ~ 0 gives V = -g2 ~ 0-), so the physical floor of THIS convention is -g2, not 0. A V below that floor, or
    # above 1, is a defect (decisions: floor -g2 instead of the spec's 0, which fails on ~60 ordinary runs by < 1e-2).
    g2f = g2 if is_num(g2) else 0.0
    for conv in ("raw", "corrected"):
        v = val(r, "hom_V", conv)
        if is_num(v) and not (-max(g2f, 0.0) - 1e-9 <= v <= 1.0 + 1e-9):
            fails.append(("bounds", {}, f"hom_V ({conv}) = {v:.4g} outside [-g2, 1] (g2 = {g2f:.3g})"))
    ep = val(r, "eta_prep")
    if is_num(ep) and not (-1e-12 <= ep <= 1.0 + 1e-9):
        fails.append(("bounds", {}, f"eta_prep = {ep!r} outside [0, 1]"))
    bl, ee, ez, nph = val(r, "B_first_lens"), val(r, "eta_extract"), val(r, "eta_pol"), val(r, "n_ph")
    if all(is_num(x) for x in (bl, ee, ez, nph)) and nph <= 1.0 + 1e-12 and bl > ee * ez + 1e-12:
        fails.append(("bounds", {}, f"B_first_lens {bl:.4g} > eta_extract x eta_pol = {ee * ez:.4g}"))
    if r["family"] == "coherent" and r["timing"] == "pulsed" and is_num(g2) and g2 >= 1.0:
        fails.append(("bounds", {}, f"coherent pulsed two-level drive with g2 = {g2:.4g} >= 1"))
    # every advertised observable is finite (first declared convention), series present with finite values
    for oid in sorted(drive_obs):
        k = first_conv_key(r, oid)
        if k is not None:
            v = r["obs"][k]
            if not is_num(v) and r.get("obs_reason", {}).get(k):
                pass                               # unavailable WITH a stated reason (R2 F6: e.g. hom_I_gated with no gate)
            elif not is_num(v):
                fails.append(("advertised", {"request": oid}, f"{k}={v if isinstance(v, str) else 'nan'}"[:90]))
        elif oid in r["series"]:
            ys = r["series"][oid]
            if not ys or not all(is_num(y) for y in ys):
                fails.append(("advertised", {"request": oid}, f"series {oid} not finite"))
        elif r.get("series_unavailable", {}).get(oid):
            pass                                   # left out on purpose, with the reason the run reports (R2 N1)
        else:
            fails.append(("advertised", {"request": oid}, f"{oid} not produced"))
    return fails


def purpose_check(r, rule):
    metric = rule["metric"]
    v = val(r, metric)
    if not (is_num(v) and v > rule["min"]):
        return [("purpose", {}, f"{metric} = {v:.4g} <= {rule['min']} at the drive defaults ({rule['source'][:70]})")]
    return []


def prep_note_check(r):
    ep = val(r, "eta_prep")
    if r["family"] == "coherent" and r["timing"] == "pulsed" and r["valid"] and is_num(ep) and ep < 0.05:
        extra = [n for n in r["notes"] if "phonons (tier P0)" not in n]
        if not any(("prepar" in n.lower() or "inversion" in n.lower() or "off-resonan" in n.lower() or "excitation" in n.lower())
                   for n in extra):
            return [("prep_note", {}, f"valid run with eta_prep = {ep:.3g} < 0.05 carries no note about the poor preparation")]
    return []


def main():
    t_start = time.time()
    args = sys.argv[1:]
    workers = int(args[args.index("--workers") + 1]) if "--workers" in args else 2          # thermal rule: at most 2 workers
    quick = "--quick" in args
    compat = yaml.safe_load(COMPAT.read_text(encoding="utf-8"))

    from fsim_core.compose import registry as R, observables as obs_mod
    mats, structs, drives, solvers = (R.kinds(a) for a in ("material", "structure", "drive", "solver"))
    rep = Report()

    # -------- the table must only name registered kinds (a typo or a removed kind fails, a new kind needs no edit)
    bad_names = []
    reg_kinds = {"material": mats, "structure": structs, "drive": drives, "solver": solvers}
    for rule in compat["excluded"]:
        for ax in ("material", "structure", "drive"):
            for k in (rule.get(ax) if isinstance(rule.get(ax), list) else [rule.get(ax)] if rule.get(ax) else []):
                if k not in reg_kinds[ax]:
                    bad_names.append(f"excluded: {ax} {k}")
    for sv, ds in compat["solver_drives"].items():
        if sv not in solvers:
            bad_names.append(f"solver_drives: solver {sv}")
        bad_names += [f"solver_drives[{sv}]: drive {d}" for d in ds if d not in drives]
    bad_names += [f"auto_solver_bulk: drive {d}" for d in compat["auto_solver_bulk"] if d not in drives]
    bad_names += [f"purpose: drive {d}" for d in compat["purpose"] if d not in drives]
    for key in ("reference_material",):
        if compat[key] not in mats:
            bad_names.append(key)
    for key, reg in (("reference_structure", structs), ("reference_drive", drives), ("jc_cavity_structure", structs)):
        if compat[key] not in reg:
            bad_names.append(key)
    for s, m in compat["structure_reference_material"].items():
        if s not in structs or m not in mats:
            bad_names.append(f"structure_reference_material {s}: {m}")
    rep.add(f"compat table names only registered kinds (registry: {len(mats)} materials, {len(structs)} structures, "
            f"{len(drives)} drives, {len(solvers)} solvers)", [("table", {}, ", ".join(bad_names))] if bad_names else [], {})
    undeclared = [d for d in drives if d not in compat["auto_solver_bulk"]]
    rep.add("every registered drive has a declared auto solver at bulk_planar",
            [("table", {}, "undeclared: " + ", ".join(undeclared))] if undeclared else [], {})

    # -------- jobs
    ref_m, ref_s, ref_d = compat["reference_material"], compat["reference_structure"], compat["reference_drive"]
    variants = {v["name"]: v for v in compat.get("variants", [])}
    jobs = {}

    def add_job(m, s, d, solver="auto", over=None, tag=""):
        key = (m, s, d, solver, tag)
        jobs.setdefault(key, {"m": m, "s": s, "d": d, "solver": solver, "over": over or {}, "tag": tag})
        return key

    slice_of = {}
    for s in structs:
        for d in drives:
            slice_of[add_job(compat["structure_reference_material"].get(s, ref_m), s, d)] = "SxD"
        for vn, v in variants.items():
            slice_of[add_job(compat["structure_reference_material"].get(s, ref_m), s, v["drive"], over={f"drive.{k}": x for k, x in v["params"].items()}, tag=vn)] = "SxD"
    for m in mats:
        for s in structs:
            slice_of.setdefault(add_job(m, s, ref_d), "MxS")
        for d in drives:
            slice_of.setdefault(add_job(m, ref_s, d), "MxD")
        for vn, v in variants.items():
            slice_of.setdefault(add_job(m, ref_s, v["drive"], over={f"drive.{k}": x for k, x in v["params"].items()}, tag=vn), "MxD")
    chg_off = add_job(ref_m, ref_s, ref_d, over={"drive.target": "X+"}, tag="chg-off")
    chg_on = add_job(ref_m, ref_s, ref_d, over={"drive.target": "X+", "material.charge_state": "X+"}, tag="chg-on")
    solver_keys = {}
    for sv in solvers:
        for d in drives:
            solver_keys[(sv, d)] = add_job(ref_m, ref_s, d, solver=sv)
        for d in compat["jc_cavity_drives"]:
            solver_keys[(sv, d, "cavity")] = add_job(ref_m, compat["jc_cavity_structure"], d, solver=sv)
    if quick:
        keep = set(list(jobs)[: max(40, len(jobs) // 8)]) | set(solver_keys.values())
        jobs = {k: v for k, v in jobs.items() if k in keep}

    # -------- dead-control tasks (one per registry entry with parameters; escalating chains, see probe_entry)
    tasks = []
    ctxs = compat["dead_control"]["context"]
    for axis in ("material", "structure", "drive", "setup"):
        for kind in R.kinds(axis):
            ent = R.get(axis, kind)
            if not ent.params or ent.status == "needs_physics":
                continue
            m, s, d = ref_m, "micropillar", ref_d
            if axis == "material":
                m = kind
            elif axis == "structure":
                s, m = kind, compat["structure_reference_material"].get(kind, ref_m)
            elif axis == "drive":
                d = kind
            chains = [((m, s, d), {})]
            if axis in ("material", "setup"):
                chains += [((m, "bulk_planar", ref_d), {}), ((m, "phc_waveguide", ref_d), {}), ((m, s, "phonon_assisted"), {}), ((m, "bulk_planar", "tpe"), {}),
                           ((m, "bulk_planar", "above_band"), {}), ((m, "bulk_planar", "cw_resonant"), {}), ((m, s, ref_d), {"T_K": 50.0})]
            elif axis == "structure":
                chains += [((m, kind, "phonon_assisted"), {}), ((m, kind, "above_band"), {}), ((m, kind, "tpe"), {}), ((m, kind, ref_d), {"T_K": 50.0})]
            elif axis == "drive":
                chains += [((m, "bulk_planar", d), {}), ((m, "bulk_planar", d), {"T_K": 50.0})]
                mech = next((q for q in ent.params if q.name == "mechanism" and q.choices), None)
                if mech:
                    chains += [((m, s, d), {"drive.mechanism": c}) for c in mech.choices]
            tasks.append({"axis": axis, "kind": kind, "chains": chains, "context": ctxs,
                          "candidates": compat["dead_control"].get("candidates", {})})
    if quick:
        tasks = tasks[:6]

    print(f"matrix: {len(jobs)} runs + {len(tasks)} dead-control entries on {workers} workers", flush=True)
    load_path = args[args.index("--load") + 1] if "--load" in args else None
    save_path = args[args.index("--save") + 1] if "--save" in args else None
    if load_path:                                       # re-judge a saved raw run without running the engine again
        import pickle
        with open(load_path, "rb") as fh:
            results, remedy_res, probe_results = pickle.load(fh)
        print(f"loaded {len(results)} runs, {len(probe_results)} probes from {load_path}", flush=True)
        miss = {k: j for k, j in jobs.items() if k not in results}
        if miss:                                        # jobs added since the save: run just those
            with mp.get_context("spawn").Pool(workers) as pool:
                for r in pool.imap_unordered(run_point, list(miss.values())):
                    j = r["job"]
                    results[(j["m"], j["s"], j["d"], j["solver"], j["tag"])] = r
            print(f"ran {len(miss)} jobs missing from the saved run", flush=True)
        if "--reprobe-dead" in args:                    # rerun only the entries that had a non-ok parameter
            redo = {(p["axis"], p["kind"]) for p in probe_results if any(v[0] not in ("ok", "skip") for v in p["result"].values())}
            keep = [p for p in probe_results if (p["axis"], p["kind"]) not in redo]
            with mp.get_context("spawn").Pool(workers) as pool:
                probe_results = keep + pool.map(probe_entry, [t for t in tasks if (t["axis"], t["kind"]) in redo], chunksize=1)
            print(f"re-probed {len(redo)} entries", flush=True)
            if save_path:
                import pickle
                with open(save_path, "wb") as fh:
                    pickle.dump((results, remedy_res, probe_results), fh)
    else:
        ctxm = mp.get_context("spawn")
        with ctxm.Pool(workers) as pool:
            results = {}
            probe_results = []
            asyncs = [pool.apply_async(probe_entry, (t,)) for t in tasks]
            for r in pool.imap_unordered(run_point, list(jobs.values()), chunksize=2):
                j = r["job"]
                results[(j["m"], j["s"], j["d"], j["solver"], j["tag"])] = r
                if len(results) % 40 == 0:
                    print(f"  ... {len(results)}/{len(jobs)} runs, {time.time() - t_start:.0f} s", flush=True)
            # phase 2: remedies named in invalid reasons ("use solver 'X'")
            import re
            remedy_jobs = {}
            for key, r in results.items():
                if r["outcome"] == "invalid":
                    for reason in r["reasons"]:
                        mm = re.search(r"use solver '([a-z_]+)'", reason)
                        if mm and mm.group(1) in solvers:
                            j = jobs[key]
                            remedy_jobs[key] = {**j, "solver": mm.group(1), "tag": "remedy:" + j["tag"]}
            remedy_res = {}
            for r in pool.imap_unordered(run_point, list(remedy_jobs.values()), chunksize=1):
                j = r["job"]
                remedy_res[(j["m"], j["s"], j["d"], j["tag"].replace("remedy:", "", 1))] = r
            for a in asyncs:
                probe_results.append(a.get(timeout=3600))

        if save_path:
            import pickle
            with open(save_path, "wb") as fh:
                pickle.dump((results, remedy_res, probe_results), fh)

    # -------- evaluate the matrix
    excl = compat["excluded"]
    n_note = 0
    for key in sorted(results):
        m, s, d, solver, tag = key
        job = jobs[key]
        dname = tag if tag and not tag.startswith("remedy") else d
        if solver != "auto" or tag.startswith("chg"):
            continue                               # solver slice and the charge-state pair are judged below
        r = results[key]
        why = excluded_reason(excl, m, s, d)
        kb = {"slice": slice_of.get(key, "?"), "material": m, "structure": s, "drive": dname, "solver": r.get("solver", "")}
        name = f"{slice_of.get(key, '?')} {m} | {s} | {dname}"
        fails = []
        out = r["outcome"]
        if why is not None:
            if out in ("refused", "unavailable", "invalid"):
                pass                                           # a clean, stated refusal
            elif out == "crash":
                fails.append(("runs", {}, f"declared excluded ({why[:60]}...) but crashed: {r['detail']}"))
            else:
                n_note += 1
                rep.notes.append(f"excluded combination now runs: {name} (update compose_compat.yaml)")
        else:
            if out != "ok":
                fails.append(("runs", {}, f"{out}: {r.get('detail') or '; '.join(r.get('reasons', []))}"[:300]))
        if out in ("ok", "invalid") and why is None and out == "ok":
            ent = R.get("drive", job["d"])
            dvals = {k.split(".", 1)[1]: x for k, x in (job.get("over") or {}).items() if k.startswith("drive.")}
            ctx = {"drive_obs": set(ent.observables_for(dvals) if hasattr(ent, "observables_for") else ent.observables)}
            fails += combo_checks(r, ctx)
            rule = compat["purpose"].get(job["d"])
            if rule and s == ref_s and not tag and slice_of.get(key) in ("SxD", "MxD") and m == ref_m:
                fails += purpose_check(r, rule)
            if not tag:
                fails += prep_note_check(r)
        if out == "invalid":
            rk = (m, s, job["d"], tag)
            rr = remedy_res.get(rk)
            if rr is not None and (rr["outcome"] != "ok"):
                fails.append(("remedy", {}, f"invalid run recommends solver '{rr['job']['solver']}', which then gives {rr['outcome']}: {rr.get('detail', '')[:140]}"))
            elif rr is not None and not rr["valid"]:
                fails.append(("remedy", {}, f"recommended solver '{rr['job']['solver']}' still invalid: {'; '.join(rr['reasons'])[:140]}"))
        rep.add(name, fails, kb)

    # -------- material.charge_state only ENABLES a transition (it moves no output by itself): the drive that targets it must
    # be refused without it and run with it
    ro, rn = results.get(chg_off), results.get(chg_on)
    if ro is not None and rn is not None:
        rep.add("charge_state enables drive.target X+: refused (named) without it, runs with it",
                ([] if ro["outcome"] in ("refused", "unavailable", "crash") and "X+" in ro.get("detail", "") else
                 [("charge_state", {}, f"without charge_state X+: {ro['outcome']} {ro.get('detail', '')[:100]}")])
                + ([] if rn["outcome"] == "ok" else [("charge_state", {}, f"with charge_state X+: {rn['outcome']} {rn.get('detail', '')[:100]}")]),
                {"slice": "charge"})

    # -------- solver x drive
    for sv in solvers:
        declared = set(compat["solver_drives"].get(sv, []))
        for d in drives:
            r = results.get(solver_keys[(sv, d)])
            if r is None:
                continue
            name = f"SvxD solver {sv} | {ref_s} | {d}"
            fails = []
            if d in declared:
                if r["outcome"] != "ok":
                    fails.append(("solver_runs", {}, f"declared accepted, got {r['outcome']}: {r.get('detail', '')[:140]}"))
                else:
                    fails += [(c, e, x) for c, e, x in combo_checks(r, {"drive_obs": set()})]   # no 'advertised' here: a solver may lack coherence
            elif r["outcome"] not in ("unavailable", "refused"):
                fails.append(("solver_runs", {}, f"not declared for {sv}, yet outcome {r['outcome']}"))
            rep.add(name, fails, {"slice": "SvxD", "solver": sv, "drive": d})
        for d in compat["jc_cavity_drives"]:
            if sv == "jc":
                r = results[solver_keys[(sv, d, "cavity")]]
                rep.add(f"SvxD solver jc | {compat['jc_cavity_structure']} | {d}",
                        [] if r["outcome"] == "ok" else [("solver_runs", {}, f"{r['outcome']}: {r.get('detail', '')[:140]}")],
                        {"slice": "SvxD", "solver": "jc", "drive": d})
    for d in drives:
        r = results.get(add_job(ref_m, ref_s, d)) if not quick else None
        if r is None:
            continue
        want = compat["auto_solver_bulk"].get(d)
        rep.add(f"auto solver for {d} at {ref_s} is {want}",
                [] if r.get("solver") == want else [("auto_solver", {}, f"got {r.get('solver')!r}")], {"slice": "auto", "drive": d})

    # -------- Composer observable requests (R2 5): re-derive composer.js line 681 (since 8073001, FIX2):
    #   req: c === "intrinsic" ? o.id : `${o.id}:${c}`      (the bare id ONLY for the convention 'intrinsic')
    bad = []
    for oid in R.kinds("observable"):
        o = obs_mod.get_observable(oid)
        for conv in o.conventions:
            req = oid if conv == "intrinsic" else f"{oid}:{conv}"
            try:
                pid, pconv = obs_mod.parse_request(req)
            except Exception as e:                                             # noqa: BLE001
                bad.append(f"{req}: {e}")
                continue
            if pid != oid or pconv != conv:
                bad.append(f"{req} parses to {pid}:{pconv}, the Composer meant {oid}:{conv}")
    rep.add(f"Composer requests parse to the convention they mean, all {len(R.kinds('observable'))} observables",
            [("composer_request", {}, "; ".join(bad)[:400])] if bad else [], {"slice": "contract"})
    try:
        from fsim_studio import api_compose
        errs = []
        for req in ("g2_0:nonsense", "hom_V", "hom_V:nonsense"):
            try:
                api_compose._check_observables([req])
                errs.append(req)
            except api_compose.ComposeError:
                pass
        rep.add("API _check_observables rejects an observable with a convention it does not have (g2_0:nonsense, hom_V, hom_V:nonsense)",
                [("api_convention", {}, "accepted: " + ", ".join(errs))] if errs else [], {"slice": "contract"})
        refused = []                                       # ... and accepts every request composer.js can send (line 681)
        for oid in R.kinds("observable"):
            for conv in obs_mod.get_observable(oid).conventions:
                req = oid if conv == "intrinsic" else f"{oid}:{conv}"
                try:
                    api_compose._check_observables([req])
                except api_compose.ComposeError as e:
                    refused.append(f"{req}: {e}")
        rep.add("API _check_observables accepts every request the Composer can send",
                [("api_convention", {}, "refused: " + "; ".join(refused)[:300])] if refused else [], {"slice": "contract"})
    except ImportError as e:
        rep.notes.append(f"api_compose not importable, API convention check skipped: {e}")

    # -------- R2 F6: a cw variant does not advertise the per-pulse observables (the timing variant decides the advertised set)
    per_pulse = {"hom_I", "hom_V", "B_first_lens", "B_fibre", "hom_vs_dt"}
    fl6 = []
    for v in variants.values():
        ent = R.get("drive", v["drive"])
        cw = set(ent.observables_for(v["params"]))
        pu = set(ent.observables_for({"cw": False}))
        if cw & per_pulse:
            fl6.append(f"{v['name']} still advertises {sorted(cw & per_pulse)}")
        if not per_pulse <= pu:
            fl6.append(f"{v['drive']} pulsed lost {sorted(per_pulse - pu)}")
        if not {"B_first_lens_Hz", "B_fibre_Hz", "rate_detected_Hz"} <= cw:
            fl6.append(f"{v['name']} lost its cw rates")
    rep.add("cw variants advertise no per-pulse observable (hom_I, hom_V, B_first_lens, B_fibre, hom_vs_dt), the pulsed drives keep them",
            [("variant_obs", {}, "; ".join(fl6))] if fl6 else [], {"slice": "contract"})

    # -------- dead controls
    for pr in sorted(probe_results, key=lambda p: (p["axis"], p["kind"])):
        for pname, (status, detail) in sorted(pr["result"].items()):
            pid = f"{pr['axis']}.{pr['kind']}.{pname}"
            ex = next((why for pat, why in compat["dead_control"]["exempt"].items() if fnmatch.fnmatch(pid, pat)), None)
            name = f"DEAD {pid} moves an output (first chain {' | '.join(pr['chains'][0])}, {len(pr['chains'])} chains)"
            if status == "ok" or status == "skip":
                rep.add(name, [], {})
            elif ex:
                rep.lines.append((True, None, f"[PASS] {name}  (exempt: {ex})"))
            else:
                rep.add(name, [("dead_control", {"param": pid}, f"{status}: {detail}")], {"slice": "dead", "param": pid})

    # -------- report
    for ok, kid, text in rep.lines:
        print(text)
    for n in rep.notes:
        print(f"NOTE {n}")
    unused = [e["id"] for e in KNOWN_OPEN if e["id"] not in rep.hits]
    for e in KNOWN_OPEN:
        if e["id"] not in rep.hits:
            print(f"[FIXED?] KNOWN_OPEN {e['id']} (finding {e['finding']}) matched nothing: delete the entry. {e['desc'][:80]}")
    n_known = sum(1 for ok, kid, _ in rep.lines if ok and kid)
    for e in KNOWN_OPEN:
        if e["id"] in rep.hits:
            print(f"KNOWN-OPEN {e['id']} finding {e['finding']}: {rep.hits[e['id']]} failing check(s) - {e['desc']}")
    total = len(rep.lines)
    passed = sum(1 for ok, _, _ in rep.lines if ok)
    print(f"{passed}/{total} compose matrix checks passed ({n_known} KNOWN-OPEN, {len(unused)} stale entries, "
          f"{time.time() - t_start:.0f} s)")
    return 0 if passed == total else 1


if __name__ == "__main__":
    mp.freeze_support()
    sys.exit(main())
