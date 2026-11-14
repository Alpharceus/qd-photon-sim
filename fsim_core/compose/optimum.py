"""Optimum items of the paper-device ledger (devlib R2 SUPCMP).

A reproduce entry with `mode: optimum` reports the MAXIMUM of one observable over a scanned parameter window
instead of one operating point. It exists because some paper numbers are themselves optima ("maximum exciton
occupation about 97 % in the considered parameter window", Karli et al. 2022); comparing them with one point is
not like-for-like. Entry layout (validated by `validate_entry`):

  - item: <id>               # class M (compared with a reported value, `compare_to` or the item id) or P
    class: M | P
    mode: optimum
    compare_to: <reported id>   # optional; default: the item id itself
    tol: {abs: ..}               # class M only, same keys as any scalar item
    optimum:
      observable: inversion      # exciton occupation after the pulses, started from |G><G| (see below)
      reduction: max
      solver: {phonons: on|off, radiative: on|off}
      overrides:                  # leaf path -> {value, unit, tag, where|source}: replaces the card leaf for this item
      scan:                       # 1-2 scanned variables, each with its window and provenance
        - {name, window: [lo, hi], n: grid points, tag, where|source, path: <leaf>   # sets the leaf directly
                                                                   | (no path: used by a `derive` rule)}
      derive:                     # leaves computed from the scan after it is applied, in order
        - {path: <leaf>, rule: swing_up_area_gauss | sqrt_power_ratio, tag, where|source, ...rule keys}

  - item: <id>
    class: P
    mode: at_optimum              # a second observable evaluated by the engine at the argmax of an earlier optimum item
    of: <optimum item id>
    observable: g2_0
    convention: raw

Ladder variant (R3-LADDER). An entry with `solver: ladder` is evaluated by the phonon-free four-level exciton-biexciton ladder
(ladder.py, Piccinini et al. 2025 Eqs 1-4) instead of the two-level polaron engine:

  - item: <id>
    class: P
    mode: optimum
    solver: ladder
    optimum:
      observable: N_X            # X-line photons per pulse pair, maximised over a FIXED grid (no refinement, no tuning)
      reduction: max
      ladder: {fss_ueV: {value, unit, tag, where}}      # the fine-structure splitting; E_b, Gamma, the pulses come from the device
      scan: [...]                # as above (area1, area2)
  - item: <id>
    class: M | P
    mode: at_optimum
    solver: ladder
    of: <ladder optimum item>
    observable: g2_0 | N_X | n_ph | pXX_end
    compare_to: <reported id>    # class M only (with tol), as for optimum items

Observable. `inversion` here is the exciton occupation rho_XX at the END of the pulse window starting from |G><G|, i.e.
the "maximum exciton occupation" of a pulse-pair map. It is NOT the engine's `inversion` scalar: that one is the
periodic fixed point and, with the engine's long Gaussian window (16 sigma), includes the radiative decay of the
whole post-pulse half window (about 18 % at 165 ps). radiative: off removes the decay (the paper's model is a von
Neumann equation); phonons: off removes the polaron rates and the sqrt(Z) renormalisation (alpha -> 0).

Search. (1) a deterministic grid, evaluated by a batched exciton-frame SU(2) propagator on the engine's own pulse
train (PulseTrain.omega) and the engine's step count (phonon-free, dissipation-free two-level model; an accelerator,
checked against the engine at the argmax); (2) a bounded zoom refinement of that grid; (3) the engine's own propagation
(polaron Generator, phonons and decay as stated) at the argmax, and, when the item keeps phonons or decay, a bounded
Nelder-Mead from the two best grid seeds. The REPORTED value is always an engine value; the maximum of (engine at
the zoomed argmax, Nelder-Mead results) is taken.
"""
from __future__ import annotations

import copy
import json
import math
from dataclasses import replace
from typing import Optional

import numpy as np

from .. import lindblad as lb
from .. import resonant_source as rs
from . import engine as E
from .params import TAGS
from . import ladder as _lad
from .polaron import PolaronSolver, sinc_gauss_area

MODES = ("optimum", "at_optimum")
RULES = ("swing_up_area_gauss", "sqrt_power_ratio")
_HBAR_MEV_PS = lb.HBAR_MEV_NS * 1e3          # meV ps
_SHAPE_KEYS = ("drive.delta1_meV", "drive.delta2_meV", "drive.fwhm1_ps", "drive.fwhm2_ps", "drive.shape", "T_K",
               "drive.slm_alpha_pi", "drive.slm_bb_fwhm_meV", "drive.slm_bb_centre_meV", "drive.slm_mask_sigma_meV", "drive.slm_floor",
               "drive.spectral_fwhm_meV", "drive.r")
LADDER_OBSERVABLES = ("g2_0", "N_X", "n_ph", "pXX_end")
_AMP_PATHS = ("drive.area1_pi", "drive.area2_pi", "drive.slm_T1", "drive.slm_T2")   # the paths that scale the field linearly
ZOOM_LEVELS = 3
ZOOM_N = 7
NZ_SEEDS = 3
NM_SEEDS = 1
NM_MAXFEV = 30
_CACHE: dict = {}


# ------------------------------------------------------------------ rules
def swing_up_area_gauss(v: dict, spec: dict) -> float:
    """Area (units of pi) of Gaussian pulse 1 from the SUPER swing-up (Rabi) condition (Karli 2022 text p.2; Bracht 2021 Eq. 6):
    |Delta2 - Delta1| = sqrt((hbar Omega1)^2 + Delta1^2) -> hbar Omega1_peak, and area = Omega_pk sigma sqrt(2 pi) with
    sigma = fwhm / (2 sqrt(ln 2)) (intensity FWHM) [DR]."""
    d1, d2 = float(v["drive.delta1_meV"]), float(v["drive.delta2_meV"])
    fw = float(v["drive.fwhm1_ps"])
    om_meV = math.sqrt(max((d2 - d1) ** 2 - d1 * d1, 0.0))
    return om_meV / _HBAR_MEV_PS * fw * math.sqrt(2.0 * math.pi) / (2.0 * math.sqrt(math.log(2.0))) / math.pi


def sqrt_power_ratio(v: dict, spec: dict) -> float:
    """area2 = area(ref) sqrt(P(x) / P1), P(x) the piecewise-linear power table [[x, P_uW], ...] [DR] (equal pulse widths:
    area scales as sqrt(power))."""
    x = float(v[spec["from"]])
    tab = sorted((float(a), float(b)) for a, b in spec["table"])
    p = float(np.interp(x, [a for a, _ in tab], [b for _, b in tab]))
    return float(v[spec["ref"]]) * math.sqrt(p / float(spec["P1_uW"]))


_RULE_FN = {"swing_up_area_gauss": swing_up_area_gauss, "sqrt_power_ratio": sqrt_power_ratio}


# ------------------------------------------------------------------ schema
def _num(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def _onoff(v):
    """YAML 1.1 reads an unquoted on / off as a bool: accept both spellings."""
    return "on" if v is True else "off" if v is False else v


def _ne(x) -> bool:
    return isinstance(x, str) and bool(x.strip())


def validate_entry(w: str, e: dict, ids: dict, resolved: dict, P: list, order: list) -> None:
    """Append the problems of an optimum / at_optimum reproduce entry to P. `order` = optimum item ids seen so far."""
    mode, cls = e.get("mode"), e.get("class")
    if mode not in MODES:
        P.append(f"{w}: mode must be one of {list(MODES)}, got {mode!r}")
        return
    if cls not in ("M", "P"):
        P.append(f"{w}: an {mode} entry must be class M or P, got {cls!r}")
    if not _ne(e.get("tol_reason")):
        P.append(f"{w}: tol_reason required (why this tolerance / why no comparison, written before the run)")
    ladder = e.get("solver") == "ladder"
    if e.get("solver") is not None and not ladder:
        P.append(f"{w}: solver must be 'ladder' when given, got {e.get('solver')!r}")

    def check_M(ref):
        if ref not in ids or ids[ref][0] != "reported":
            P.append(f"{w}: class M {mode} item compares with a reported scalar; {ref!r} is not a reported id")
            return
        tol = e.get("tol")
        if not isinstance(tol, dict) or not tol:
            P.append(f"{w}: tol required")
            return
        bad = set(tol) - {"k_sigma", "abs", "rel"}
        if bad:
            P.append(f"{w}: unknown tol keys {sorted(bad)}")
        if "k_sigma" in tol and ids[ref][1].get("err") is None:
            P.append(f"{w}: tol.k_sigma needs the reported item to carry err")
        for k, v in tol.items():
            if not (_num(v) and v >= 0):
                P.append(f"{w}: tol.{k} must be a non-negative number")

    if mode == "at_optimum":
        if cls == "M":
            if not ladder:
                P.append(f"{w}: a class M at_optimum item needs solver: ladder (the polaron at_optimum items are class P)")
            check_M(e.get("compare_to", e.get("item")))
        elif e.get("compare_to") is not None:
            P.append(f"{w}: class P carries no comparison (compare_to); state the reference in tol_reason")
        if e.get("of") not in order:
            P.append(f"{w}: `of` must name an EARLIER optimum item, got {e.get('of')!r}")
        if not _ne(e.get("observable")):
            P.append(f"{w}: at_optimum needs an `observable`")
        elif ladder and e["observable"] not in LADDER_OBSERVABLES:
            P.append(f"{w}: a ladder at_optimum observable must be one of {list(LADDER_OBSERVABLES)}")
        return
    order.append(e.get("item"))
    ref = e.get("compare_to", e.get("item"))
    if cls == "M":
        check_M(ref)
    elif cls == "P" and e.get("compare_to") is not None:
        P.append(f"{w}: class P carries no comparison (compare_to); state the reference in tol_reason")
    o = e.get("optimum")
    if not isinstance(o, dict):
        P.append(f"{w}: optimum block required")
        return
    if ladder:
        if o.get("observable") != "N_X":
            P.append(f"{w}: a ladder optimum.observable must be 'N_X' (X-line photons per pulse pair)")
    elif o.get("observable") != "inversion":
        P.append(f"{w}: optimum.observable must be 'inversion' (exciton occupation after the pulses)")
    if o.get("reduction") != "max":
        P.append(f"{w}: optimum.reduction must be 'max'")
    s = o.get("solver")
    if not ladder and not (isinstance(s, dict) and _onoff(s.get("phonons")) in ("on", "off") and _onoff(s.get("radiative")) in ("on", "off")):
        P.append(f"{w}: optimum.solver must be {{phonons: on|off, radiative: on|off}}")

    def leaf_ok(where, d):
        if not (isinstance(d, dict) and d.get("tag") in TAGS and (_ne(d.get("where")) or _ne(d.get("source")))):
            P.append(f"{where}: needs tag and where|source")
        elif d["tag"] in ("V", "DR") and not _ne(d.get("where")):
            P.append(f"{where}: tag {d['tag']} needs where")

    if ladder:
        lad = o.get("ladder")
        fss = lad.get("fss_ueV") if isinstance(lad, dict) else None
        if not (isinstance(fss, dict) and _num(fss.get("value")) and fss["value"] >= 0):
            P.append(f"{w}: a ladder optimum needs ladder.fss_ueV with a non-negative value")
        else:
            leaf_ok(f"{w}: ladder.fss_ueV", fss)
    if ladder:
        for sv in (o.get("scan") or []):
            if isinstance(sv, dict) and sv.get("path") not in ("drive.area1_pi", "drive.area2_pi"):
                P.append(f"{w}: a ladder scan may vary only drive.area1_pi / drive.area2_pi, got {sv.get('path')!r}")
        if o.get("derive"):
            P.append(f"{w}: a ladder optimum takes no derive rules")
    for path, leaf in (o.get("overrides") or {}).items():
        if path not in resolved:
            P.append(f"{w}: overrides {path!r} is not a device leaf path")
        elif not (isinstance(leaf, dict) and (_num(leaf.get("value")) or _ne(leaf.get("value")))):
            P.append(f"{w}: overrides {path} needs a numeric (or choice) value")
        else:
            leaf_ok(f"{w}: overrides {path}", leaf)
    scan = o.get("scan")
    names = []
    if not (isinstance(scan, list) and 1 <= len(scan) <= 2):
        P.append(f"{w}: optimum.scan must list 1 or 2 scanned variables")
        scan = []
    for i, sv in enumerate(scan):
        ws = f"{w}: scan[{i}]"
        if not isinstance(sv, dict) or not _ne(sv.get("name")):
            P.append(f"{ws}: needs a name")
            continue
        names.append(sv["name"])
        win = sv.get("window")
        if not (isinstance(win, list) and len(win) == 2 and all(_num(x) for x in win) and win[0] < win[1]):
            P.append(f"{ws}: window must be [lo, hi] numbers with lo < hi")
        if not (isinstance(sv.get("n"), int) and sv["n"] >= 3):
            P.append(f"{ws}: n must be an integer >= 3")
        leaf_ok(ws, sv)
        if sv.get("path") is not None and sv["path"] not in resolved:
            P.append(f"{ws}: path {sv['path']!r} is not a device leaf path")
    for i, d in enumerate(o.get("derive") or []):
        wd = f"{w}: derive[{i}]"
        if not isinstance(d, dict) or d.get("rule") not in RULES:
            P.append(f"{wd}: rule must be one of {list(RULES)}")
            continue
        if d.get("path") not in resolved:
            P.append(f"{wd}: path {d.get('path')!r} is not a device leaf path")
        leaf_ok(wd, d)
        if d["rule"] == "sqrt_power_ratio":
            tab = d.get("table")
            if not (d.get("from") in names and _ne(d.get("ref")) and _num(d.get("P1_uW")) and d["P1_uW"] > 0
                    and isinstance(tab, list) and len(tab) >= 2
                    and all(isinstance(r, list) and len(r) == 2 and all(_num(x) for x in r) for r in tab)):
                P.append(f"{wd}: sqrt_power_ratio needs from (a scan name), ref (a leaf path), P1_uW > 0, table [[x, P_uW], ..]")
    free = {sv.get("path") for sv in scan if sv.get("path")} | {d.get("path") for d in (o.get("derive") or [])}
    if not (free & set(_AMP_PATHS)):
        P.append(f"{w}: the scan must move at least one pulse amplitude (a pulse area, or the slm_T1 / slm_T2 transmission)")


# ------------------------------------------------------------------ point construction
def _base_values(entries, leaves, overrides: dict) -> dict:
    vals = E._central_values(leaves)
    for p, leaf in overrides.items():
        vals[p] = float(leaf["value"]) if _num(leaf["value"]) else leaf["value"]
    return vals


def _point(base: dict, o: dict, x: dict) -> dict:
    v = dict(base)
    for sv in o["scan"]:
        if sv.get("path"):
            v[sv["path"]] = float(x[sv["name"]])
        else:
            v[sv["name"]] = float(x[sv["name"]])
    for d in o.get("derive") or []:
        v[d["path"]] = float(_RULE_FN[d["rule"]](v, d))
    for sv in o["scan"]:
        v.setdefault(sv["name"], float(x[sv["name"]]))
    return v


def _axis(sv: dict, lo: float, hi: float, n: int) -> np.ndarray:
    return np.linspace(lo, hi, n)


# ------------------------------------------------------------------ engine evaluator
def _engine_px_end(entries, vals: dict, phonons: bool, radiative: bool) -> float:
    """rho_XX at the end of the pulse window from |G><G|, by the engine's own polaron Generator (same step count as
    PolaronSolver.run)."""
    em, setup, sr, rt, plan, v, reads = E._prepare(entries, vals, float(vals["T_K"]))
    em2 = em if phonons else replace(em, phonon=None)
    rt2 = rt if radiative else replace(rt, Gamma_ns=0.0, k_nr_ns=0.0, gamma_star_ns=0.0)
    sol = PolaronSolver()
    ok, why = sol.supports(plan, em2, rt2)
    if not ok:
        raise RuntimeError(why)
    pm, train, gen = sol._setup(plan, em2, rt2)
    _, nq0, dt = sol._grids(train, pm.B)
    m = max(1, int(math.ceil(train.W / nq0 / dt)))
    n_f = nq0 * m
    _, _, U = gen.grid(n_f)
    r = np.zeros(4, dtype=complex)
    r[0] = 1.0                                  # vec(|G><G|), column stacking, basis (G, X)
    for u in U:
        r = u @ r
    return float(r[3].real)


def pulse_end_check(spec, phonons: bool = True):
    """Raise ValueError (user-facing physics wording) when `pulse_end_occupation` cannot read out this design; else return
    (entries, vals). Cheap (builds the drive plan only), so a server can refuse before it queues the real computation."""
    spec = E._as_spec(spec)
    entries = E._entries(spec)
    vals = E._central_values(E._leaves(spec, entries))
    em, _setup, _sr, rt, plan, _v, _reads = E._prepare(entries, vals, float(vals["T_K"]))
    if plan.family != "coherent" or plan.timing != "pulsed":
        raise ValueError("The pulse-end occupation is defined for pulsed coherent drives only: this drive is not a "
                         "pulse sequence that starts from the ground state, so there is no end-of-pulse state to read out.")
    ok, why = PolaronSolver().supports(plan, em if phonons else replace(em, phonon=None), rt)
    if not ok:
        raise ValueError("The pulse-end occupation cannot be computed for this design: "
                         + str(why).replace("polaron: ", "") + ".")
    return entries, vals


def pulse_end_occupation(spec, phonons: bool = True, radiative: bool = True) -> dict:
    """rho_XX at the END of the pulse window, started from |G><G|, for a composed design (the paper-style readout of the
    optimum items) with the LA-phonon polaron rates and the radiative / non-radiative decay switched on or off [DR].
    Ranged leaves collapse to their central values, as everywhere in the Composer. The occupation is the engine's own
    `_engine_px_end`; this function only turns a spec into its inputs. Raises ValueError (user-facing wording) when the
    drive is not a pulsed coherent drive or the polaron solver does not cover the design."""
    entries, vals = pulse_end_check(spec, phonons)
    value = _engine_px_end(entries, vals, bool(phonons), bool(radiative))
    return {"value": float(value), "phonons": bool(phonons), "radiative": bool(radiative), "T_K": float(vals["T_K"]),
            "note": "Exciton occupation at the end of the pulse window, started from the ground state; "
                    "not the steady-state (periodic) readout."}


# ------------------------------------------------------------------ fast batched evaluator (phonon-free, lossless)
def _unit_series(train, n_f: int):
    """(t_mid, [unit complex Rabi series per component for peak amplitude 1 rad/ns], [peak amp per rad of area])."""
    edges = np.linspace(0.0, train.W, n_f + 1)
    tm = 0.5 * (edges[:-1] + edges[1:])
    series, norms = [], []
    for j, c in enumerate(train.comps):
        t2 = copy.copy(train)
        t2.comps = tuple(replace(cc, amp=1.0 if k == j else 0.0) for k, cc in enumerate(train.comps))
        series.append(t2.omega(tm))
        norms.append(c.dur_ns if c.shape == "rect" else sinc_gauss_area(c.tp_ns, c.r) if c.shape == "sinc_gauss"
                     else c.sigma_ns * math.sqrt(2.0 * math.pi))
    return float(train.W / n_f), series, norms


def _fast_batch(entries, pts: list) -> np.ndarray:
    """rho_XX(end) of a list of full value dicts (phonons off, no decay), batched by pulse-shape group."""
    out = np.full(len(pts), np.nan)
    groups: dict = {}
    for i, v in enumerate(pts):
        groups.setdefault(tuple(v[k] for k in _SHAPE_KEYS if k in v), []).append(i)
    sol = PolaronSolver()
    for idx in groups.values():
        vmax = dict(pts[idx[0]])
        for k in _AMP_PATHS:
            vmax[k] = max(pts[i][k] for i in idx)
        em, setup, sr, rt, plan, v, reads = E._prepare(entries, vmax, float(vmax["T_K"]))
        pm, train, gen = sol._setup(plan, replace(em, phonon=None), rt)
        n_f, _, _ = sol._grids(train, 1.0)
        if train.field:                                  # pulse-shaper pair: Omega = sqrt(T1) F1 + sqrt(T2) F2 + F0
            tm = 0.5 * (np.linspace(0.0, train.W, n_f + 1)[:-1] + np.linspace(0.0, train.W, n_f + 1)[1:])
            series = list(train.field.series(tm - 0.5 * train.W))
            dt = float(train.W / n_f)
            a = [np.sqrt(np.array([pts[i]["drive.slm_T1"] for i in idx])), np.sqrt(np.array([pts[i]["drive.slm_T2"] for i in idx])),
                 np.ones(len(idx))]
        else:
            dt, series, norms = _unit_series(train, n_f)
            a = [np.array([pts[i][f"drive.area{k + 1}_pi"] for i in idx]) * math.pi / norms[k] for k in range(len(series))]
        delta = float(train.det0_ns)
        g = np.ones(len(idx), dtype=complex)
        x = np.zeros(len(idx), dtype=complex)
        for s in range(n_f):
            om = a[0] * series[0][s]
            for k in range(1, len(series)):
                om = om + a[k] * series[k][s]
            h2 = delta * delta + (om.real ** 2 + om.imag ** 2)
            wv = 0.5 * np.sqrt(h2)
            c = np.cos(wv * dt)
            sn = np.where(wv > 0, np.sin(wv * dt) / np.where(wv > 0, wv, 1.0), dt)
            g, x = (c * g - 0.5j * sn * (-delta * g + np.conj(om) * x),
                    c * x - 0.5j * sn * (om * g + delta * x))
        out[idx] = np.abs(x) ** 2
    return out


# ------------------------------------------------------------------ search
def _grid_points(o, windows, ns):
    axes = [_axis(sv, windows[sv["name"]][0], windows[sv["name"]][1], ns[sv["name"]]) for sv in o["scan"]]
    mesh = np.meshgrid(*axes, indexing="ij")
    xs = [{sv["name"]: float(m.flat[k]) for sv, m in zip(o["scan"], mesh)} for k in range(mesh[0].size)]
    return axes, xs


def optimise(paper_device: dict, o: dict) -> dict:
    """Run the optimum search. -> {value, argmax(scan), leaves(full), info}."""
    if o.get("observable") == "N_X":                    # R3-LADDER: the four-level ladder has its own fixed-grid search
        return optimise_ladder(paper_device, o)
    spec = E._as_spec(paper_device)
    entries = E._entries(spec)
    leaves = E._leaves(spec, entries)
    base = _base_values(entries, leaves, o.get("overrides") or {})
    phon, rad = _onoff(o["solver"]["phonons"]) == "on", _onoff(o["solver"]["radiative"]) == "on"
    info: dict = {"n_fast": 0, "n_engine": 0}
    scan = o["scan"]
    windows = {sv["name"]: (float(sv["window"][0]), float(sv["window"][1])) for sv in scan}
    full = {sv["name"]: (float(sv["window"][0]), float(sv["window"][1])) for sv in scan}
    ns = {sv["name"]: int(sv["n"]) for sv in scan}

    def fast_search():
        """Level-0 grid, then a zoom refinement of its NZ_SEEDS best well-separated cells; the best refined cell wins."""
        axes, xs = _grid_points(o, windows, ns)
        pts = [_point(base, o, x) for x in xs]
        vals = _fast_batch(entries, pts)
        info["n_fast"] += len(pts)
        order = np.argsort(-np.nan_to_num(vals, nan=-1.0))
        span = {nm: full[nm][1] - full[nm][0] for nm in full}
        seeds = []
        for j in order:                                 # best cells separated by > 10 % of a window in some variable
            if all(any(abs(xs[j][nm] - sx[nm]) > 0.1 * span[nm] for nm in span) for _, sx in seeds):
                seeds.append((float(vals[j]), xs[j]))
            if len(seeds) == NZ_SEEDS:
                break
        info["seeds"] = []
        best = None
        for v0, x0 in seeds:
            cur = (v0, x0)
            cur_w = {nm: (max(full[nm][0], x0[nm] - (axes[i][1] - axes[i][0])), min(full[nm][1], x0[nm] + (axes[i][1] - axes[i][0])))
                     for i, nm in enumerate(span)}
            for level in range(ZOOM_LEVELS):
                zn = {nm: ZOOM_N for nm in span}
                ax2, xs2 = _grid_points(o, cur_w, zn)
                p2 = [_point(base, o, x) for x in xs2]
                v2 = _fast_batch(entries, p2)
                info["n_fast"] += len(p2)
                k = int(np.argmax(np.nan_to_num(v2, nan=-1.0)))
                if v2[k] >= cur[0]:
                    cur = (float(v2[k]), xs2[k])
                cur_w = {nm: (max(full[nm][0], cur[1][nm] - (ax2[i][1] - ax2[i][0])), min(full[nm][1], cur[1][nm] + (ax2[i][1] - ax2[i][0])))
                         for i, nm in enumerate(span)}
            info["seeds"].append(cur)
            if best is None or cur[0] > best[0]:
                best = cur
        info["seeds"].sort(key=lambda c: -c[0])
        return best

    best_fast = fast_search()
    info["fast_value"], info["fast_argmax"] = best_fast
    phon_free_lossless = (not phon) and (not rad)
    x0 = best_fast[1]
    v0 = _point(base, o, x0)
    val0 = _engine_px_end(entries, v0, phon, rad)
    info["n_engine"] += 1
    results = [(val0, x0)]
    if not phon_free_lossless:
        # the lossless accelerator ranks nearly degenerate maxima by noise: score every refined seed with the engine
        # (phonons / decay as stated), then run a bounded Nelder-Mead from the best of them
        for _, xs in info["seeds"]:
            if xs is not x0:
                results.append((_engine_px_end(entries, _point(base, o, xs), phon, rad), xs))
                info["n_engine"] += 1
        from scipy.optimize import minimize
        lo = np.array([full[s["name"]][0] for s in scan])
        hi = np.array([full[s["name"]][1] for s in scan])

        def f(u):
            u = np.clip(u, 0.0, 1.0)
            x = {s["name"]: float(lo[i] + u[i] * (hi[i] - lo[i])) for i, s in enumerate(scan)}
            info["n_engine"] += 1
            return -_engine_px_end(entries, _point(base, o, x), phon, rad)

        for _, xs in sorted(results, key=lambda r: -r[0])[:NM_SEEDS]:
            u0 = np.array([(xs[s["name"]] - lo[i]) / (hi[i] - lo[i]) for i, s in enumerate(scan)])
            r = minimize(f, u0, method="Nelder-Mead", bounds=[(0.0, 1.0)] * len(scan),
                         options={"maxfev": NM_MAXFEV, "xatol": 2e-4, "fatol": 1e-5, "initial_simplex":
                                  np.vstack([u0] + [np.clip(u0 + 0.03 * np.eye(len(scan))[i], 0, 1) for i in range(len(scan))])})
            u = np.clip(r.x, 0.0, 1.0)
            results.append((-float(r.fun), {s["name"]: float(lo[i] + u[i] * (hi[i] - lo[i])) for i, s in enumerate(scan)}))
    val, xb = max(results, key=lambda r: r[0])
    vfull = _point(base, o, xb)
    edge = [s["name"] for s in scan if min(abs(xb[s["name"]] - full[s["name"]][0]), abs(xb[s["name"]] - full[s["name"]][1]))
            < 1e-6 * (full[s["name"]][1] - full[s["name"]][0])]
    # agreement of the accelerator with the engine at the argmax (phonon-free, lossless)
    chk = _engine_px_end(entries, vfull, False, False) if not phon_free_lossless else val
    fast_at = float(_fast_batch(entries, [vfull])[0])
    info["fast_vs_engine_phonon_free"] = abs(fast_at - chk)
    info["on_boundary"] = edge
    return {"value": float(val), "argmax": xb, "leaves": {k: vfull[k] for k in vfull if k.startswith("drive.area") or
                                                           k.startswith("drive.delta") or k.startswith("drive.fwhm")},
            "values": vfull, "info": info, "entries": entries}


def optimise_ladder(paper_device: dict, o: dict) -> dict:
    """Ladder optimum (R3-LADDER): N_X maximised over a FIXED grid of the two pulse areas (no zoom, no refinement: the argmax is
    a grid point chosen before any g2 is read). -> {value, argmax(scan), leaves, values, info, entries, setup}."""
    spec = E._as_spec(paper_device)
    entries = E._entries(spec)
    leaves = E._leaves(spec, entries)
    base = _base_values(entries, leaves, o.get("overrides") or {})
    scan = o["scan"]
    bad = [sv.get("path") for sv in scan if sv.get("path") not in ("drive.area1_pi", "drive.area2_pi")] + \
          [d.get("path") for d in o.get("derive") or []]
    if bad:
        raise ValueError(f"ladder optimum: the scan may vary only drive.area1_pi / drive.area2_pi (got {bad})")
    windows = {sv["name"]: (float(sv["window"][0]), float(sv["window"][1])) for sv in scan}
    ns = {sv["name"]: int(sv["n"]) for sv in scan}
    _, xs = _grid_points(o, windows, ns)
    pts = [_point(base, o, x) for x in xs]
    em, setup, sr, rt, plan, v, reads = E._prepare(entries, pts[0], float(pts[0]["T_K"]))
    fss = float(o["ladder"]["fss_ueV"]["value"])
    S = _lad.prepare(plan, em, rt, fss_ueV=fss)
    a = [np.array([pt[f"drive.area{k + 1}_pi"] for pt in pts]) for k in range(len(S.series))]
    res = _lad.scan(S, a)
    k = int(np.argmax(res["N_X"]))
    xb = xs[k]
    vfull = pts[k]
    edge = [s["name"] for s in scan if min(abs(xb[s["name"]] - windows[s["name"]][0]), abs(xb[s["name"]] - windows[s["name"]][1]))
            < 1e-6 * (windows[s["name"]][1] - windows[s["name"]][0])]
    info = {"n_fast": 0, "n_engine": 0, "n_ladder": len(pts), "on_boundary": edge, "solver": "ladder",
            "fast_vs_engine_phonon_free": float("nan"), "ladder": dict(S.info), "pXX_end": float(res["pXX_end"][k]),
            "pX_end": float(res["pX_end"][k]), "grid_max_N_X": float(res["N_X"][k])}
    return {"value": float(res["N_X"][k]), "argmax": xb, "leaves": {kk: vfull[kk] for kk in vfull if kk.startswith("drive.area")
                                                                     or kk.startswith("drive.delta") or kk.startswith("drive.fwhm")},
            "values": vfull, "info": info, "entries": entries, "setup": S, "ladder": True}


def _fmt_argmax(res: dict) -> dict:
    out = {k: round(v, 5) for k, v in res["argmax"].items()}
    slm = res["values"].get("drive.shape") == "slm_mask"           # the pulse-shaper pulses ignore the area leaves
    for k in ("drive.area1_pi", "drive.area2_pi", "drive.delta2_meV"):
        if k in res["leaves"] and not (slm and k.startswith("drive.area")):
            out[k.split(".")[1]] = round(res["leaves"][k], 5)
    return out


def run_item(paper, e: dict, item, store: dict) -> None:
    """Fill `item` (an ItemResult) for an optimum / at_optimum entry. `store` carries argmax results between items."""
    key = (paper.id, e["item"], json.dumps(paper.device, sort_keys=True, default=str),
           json.dumps(e.get("optimum") or {"of": e.get("of")}, sort_keys=True, default=str))
    if e["mode"] == "optimum":
        ladder = e.get("solver") == "ladder"
        if key not in _CACHE:
            _CACHE[key] = (optimise_ladder if ladder else optimise)(paper.device, e["optimum"])
        res = _CACHE[key]
        store[e["item"]] = res
        item.model = res["value"]
        item.argmax = _fmt_argmax(res)
        info = res["info"]
        if ladder:
            item.note = (f"max over {[sv['name'] for sv in e['optimum']['scan']]} of the X-line photon number N_X, four-level ladder "
                         f"(no phonons, FSS {info['ladder']['fss_ueV']:g} ueV), fixed grid of {info['n_ladder']} points; argmax "
                         f"{item.argmax}; P_XX(end) {info['pXX_end']:.3f}"
                         + (f"; ON the window boundary {info['on_boundary']}" if info["on_boundary"] else ""))
        else:
            s = e["optimum"]["solver"]
            item.note = (f"max over {[sv['name'] for sv in e['optimum']['scan']]} (phonons {_onoff(s['phonons'])}, decay {_onoff(s['radiative'])}); "
                         f"argmax {item.argmax}" + (f"; ON the window boundary {info['on_boundary']}" if info["on_boundary"] else "")
                         + f"; {info['n_fast']} grid + {info['n_engine']} engine evaluations; accelerator vs engine "
                           f"|diff| {info['fast_vs_engine_phonon_free']:.1e}")
        if item.kind == "scalar":
            from . import papers as _pp
            item.outcome = _pp._scalar_outcome(item.paper, item.err, item.model, 0.0, 0.0, item.tol)
        else:
            item.outcome = "pass"
        return
    src = store.get(e["of"])
    if src is None:
        item.outcome, item.note = "error", f"optimum item {e['of']!r} was not computed"
        return
    key2 = (key, e["of"])
    if e.get("solver") == "ladder":
        _run_ladder_at(e, item, src)
        return
    req = e["observable"] if e.get("convention", "intrinsic") == "intrinsic" else f"{e['observable']}:{e['convention']}"
    try:
        ctx = E._run_point(src["entries"], src["values"], float(src["values"]["T_K"]), "auto")
        valid, why = E._validity(ctx)
        item.model = float(E._value_of(ctx, req))
    except Exception as ex:                                       # noqa: BLE001
        item.outcome, item.note = "error", f"engine error at the argmax: {type(ex).__name__}: {ex}"
        return
    item.argmax = _fmt_argmax(src)
    if not valid:
        item.outcome, item.note = "unavailable", "model run at the argmax flagged invalid: " + "; ".join(why)
    elif not math.isfinite(item.model):
        item.outcome, item.note = "unavailable", "model value not finite at the argmax"
    else:
        item.outcome = "pass"
        item.note = f"{req} by the engine at the argmax of {e['of']}: {item.argmax}"
    _ = key2


def _run_ladder_at(e: dict, item, src: dict) -> None:
    """at_optimum item evaluated by the ladder at the argmax of a ladder optimum item (class M: compared with the reported value)."""
    if not src.get("ladder"):
        item.outcome, item.note = "error", f"optimum item {e['of']!r} was not computed by the ladder solver"
        return
    v = src["values"]
    try:
        r = _lad.g2_at(src["setup"], [float(v[f"drive.area{k + 1}_pi"]) for k in range(len(src["setup"].series))])
        item.model = float(r[e["observable"]])
    except Exception as ex:                                       # noqa: BLE001
        item.outcome, item.note = "error", f"ladder error at the argmax: {type(ex).__name__}: {ex}"
        return
    item.argmax = _fmt_argmax(src)
    if not math.isfinite(item.model):
        item.outcome, item.note = "unavailable", "ladder value not finite at the argmax"
        return
    if item.kind == "scalar":
        from . import papers as _pp
        item.outcome = _pp._scalar_outcome(item.paper, item.err, item.model, 0.0, 0.0, item.tol)
    else:
        item.outcome = "pass"
    item.note = (f"{e['observable']} by the four-level ladder (X line, no phonons) at the argmax of {e['of']}: {item.argmax}; "
                 f"N_X {r['N_X']:.4f}, P_XX(end) {r['pXX_end']:.3f}, g2(0) {r['g2_0']:.5f}")
