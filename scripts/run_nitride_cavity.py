"""Planar InGaN/GaN cavity sweep; this is deliberately only evaluator plumbing.

All model numbers in the artifacts are returned by ``fsim_core.device.evaluate``.
The grid endpoints are design assumptions [A], while the Q=167 diagnostic is
Taylor et al., Nanoscale Research Letters 5, 608 (2010), DOI
10.1007/s11671-009-9514-4 [V] (its transfer to this design is [A]).
"""
from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import itertools
import json
import os
import platform
import sys
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import fsim_core.device as device  # public module attribute: verifier instruments this

G2, FLUX = .5, 1000.0  # [A] brief/design collected first-lens floor
HEAD = {"height_nm": [1.,2.,3.,4.,5.], "radius_nm": [5.,10.,15.],
        "x_in": [.15,.25,.40], "T_hs": [230.,250.,273.,300.],
        "Q": [500.,2000.,10000.], "current_uA": [.002,.02,.2],
        "regime": ["rectangular", "deterministic_pair"]}
REF = dict(height_nm=3., radius_nm=10., x_in=.25, T_hs=300., Q=2000., current_uA=.02)

def build_grid(quick=False):
    """Frozen headline axes [A brief/design]; quick remains a 48-row subset."""
    if not quick: return copy.deepcopy(HEAD)
    return {"height_nm":[1.,3.,5.], "radius_nm":[10.], "x_in":[.25],
            "T_hs":HEAD["T_hs"], "Q":[500.,10000.], "current_uA":[.02],
            "regime":HEAD["regime"]}

def _set(d, path, value):
    if path.startswith("nitride."):
        obj = d.nitride
        bits = path.split(".")[1:]
        for bit in bits[:-1]: obj = obj[bit]
        obj[bits[-1]] = value
    else:
        obj, leaf = path.rsplit(".", 1)
        setattr(getattr(d, obj), leaf, value)

def _design(regime, p):
    card = ROOT / "cards" / ("nitride-cavity-set-design.yaml" if regime == "deterministic_pair" else "nitride-cavity-pulse-design.yaml")
    d = device.DeviceDesign.load(card)
    for k, v in p.items():
        paths = {"height_nm":"nitride.dot.height_nm", "radius_nm":"nitride.dot.radius_nm",
                 "x_in":"nitride.dot.x_in", "T_hs":"thermal.T_hs", "Q":"nitride.cavity.Q",
                 "current_uA":"drive.I_uA"}
        if k in paths: _set(d, paths[k], v)
        elif k == "detuning_offset_meV": _set(d, "nitride.cavity.detuning_offset_meV", v)
        elif k == "background_tau_ns": _set(d, "nitride.background_tau_ns", v)
        elif k == "purcell_enabled": _set(d, "nitride.cavity.purcell_enabled", v)
    d.drive.cycle_loading = regime
    return d

def _finite(v):
    return isinstance(v, (int,float)) and not isinstance(v,bool) and v == v and abs(v) != float("inf")

def _row(row_id, kind, regime, p, cache):
    d = _design(regime, p)
    identity = json.dumps({"regime":regime, "p":p}, sort_keys=True, separators=(",",":"))
    if identity not in cache:
        cache[identity] = device.evaluate(d, T_grid=[p["T_hs"]])["scalars"]
    s = cache[identity]
    r = {"row_id":row_id, "row_kind":kind, "regime":regime, **p,
         "cache_identity":identity, "g2":s.get("g2_op"), "signal_flux_s":s.get("collected_flux_pulsed_s"),
         "valid":bool(s.get("valid")), "invalid_reasons":json.dumps(s.get("invalid_reasons", [])),
         "provenance":"[A/E/DR] evaluator output; Taylor et al. NRL 2010 Q=167 [V] where used"}
    # Keep every scalar diagnostic flat enough for independent CSV replay.
    for k,v in s.items():
        if isinstance(v, (str,int,float,bool)) or v is None: r[k] = v
    r["eligible"] = bool(r["valid"] and _finite(r["g2"]) and _finite(r["signal_flux_s"]) and r["signal_flux_s"] >= FLUX)
    r["one_pair_valid"] = bool(s.get("one_pair_valid", regime != "deterministic_pair"))
    r["optical_pass"] = bool(r["eligible"] and r["g2"] < G2 and (regime != "deterministic_pair" or r["one_pair_valid"]))
    r["hardware_feasible"] = bool(s.get("set_feasible", True) and s.get("pair_supply_possible", True)) if regime == "deterministic_pair" else True
    r["device_pass"] = bool(r["optical_pass"] and r["hardware_feasible"])
    return r

def compute_stats(rows):
    """Pure same-row aggregation; invalid values are never promoted to passes."""
    valid = [r for r in rows if r.get("valid")]
    eligible = [r for r in rows if r.get("eligible")]
    gs = [r["g2"] for r in eligible if _finite(r.get("g2"))]
    fs = [r["signal_flux_s"] for r in valid if _finite(r.get("signal_flux_s"))]
    gs.sort()
    return {"total":len(rows), "valid":len(valid), "eligible":len(eligible),
            "g2_min":min(gs) if gs else None,
            "g2_median_eligible":gs[len(gs)//2] if gs else None,
            "diag_g2_min":min((r["g2"] for r in valid if _finite(r.get("g2"))), default=None),
            "flux_max":max(fs, default=None), "idealized_pass_count":sum(r.get("optical_pass",False) for r in rows),
            "hardware_pass_count":sum(r.get("device_pass",False) for r in rows),
            "hardware_infeasible_count":sum(r.get("optical_pass",False) and not r.get("hardware_feasible",False) for r in rows)}

def compute_verdict(rows, *, complete):
    by = {}
    for regime in HEAD["regime"]:
        rs = [r for r in rows if r["regime"] == regime and r["row_kind"] == "headline"]
        st = compute_stats(rs); n = len(rs)
        passed = st["hardware_pass_count"]
        status = "FAIL" if not passed else ("CONDITIONAL" if complete else "FAIL")
        if regime == "deterministic_pair" and st["idealized_pass_count"] and not passed: status = "FAIL"
        by[regime] = {**st, "regime":regime, "complete":complete, "status":status,
                      "headline_coverage":f"{n}/{n}", "eligible_fraction":(st["eligible"]/n if n else 0),
                      "coverage_over_eligible":(passed/st["eligible"] if st["eligible"] else None),
                      "flux_floor_excluded":n-st["eligible"], "evidence":"incomplete", "conditional":status=="CONDITIONAL",
                      "T_pass_min":min((r["T_hs"] for r in rs if r.get("device_pass")), default=None)}
    return by

def _write_csv(path, rows):
    keys = sorted({k for r in rows for k in r})
    with path.open("w", newline="", encoding="utf-8") as f:
        w=csv.DictWriter(f, fieldnames=keys, extrasaction="ignore"); w.writeheader(); w.writerows(rows)

def _plot(path, rows, x, title):
    os.environ.setdefault("MPLBACKEND", "Agg")
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots(figsize=(7,4)); ax2=ax.twinx()
    for regime, marker in (("rectangular","o"),("deterministic_pair","s")):
        candidates=[r for r in rows if r["regime"]==regime and _finite(r.get(x))]
        groups={}
        for r in candidates:
            key=() if x == "T_hs" else (r["T_hs"],)
            groups.setdefault(key,[]).append(r)
        for key, rs in groups.items():
            rs.sort(key=lambda r:r[x]); suffix="" if not key else f" T_hs={key[0]:g}K"
            ax.plot([r[x] for r in rs],[r["g2"] for r in rs], marker=marker,label=regime+suffix+" g2")
            flux=[max(r["signal_flux_s"], 1e-12) if _finite(r.get("signal_flux_s")) else 1e-12 for r in rs]
            ax2.plot([r[x] for r in rs],flux,marker=marker,ls="--",label=regime+suffix+" flux")
    ax.axhline(G2,color="k",lw=.8); ax2.axhline(FLUX,color="k",lw=.8); ax2.set_yscale("log")
    ax.set_xlabel(x); ax.set_ylabel("g2"); ax2.set_ylabel("first-lens useful flux / s"); ax.set_title(title+" assumed-planar model")
    if ax.lines: ax.legend(loc="upper left",fontsize=7)
    if ax2.lines: ax2.legend(loc="lower right",fontsize=7)
    fig.tight_layout(); fig.savefig(path,dpi=140); plt.close(fig)

def _safe_out(value):
    base=(ROOT/"out"/"nitride_cavity").resolve(); target=Path(value).resolve()
    if target != base and base not in target.parents: raise ValueError("out-dir must remain under out/nitride_cavity")
    return target

def main(argv=None):
    ap=argparse.ArgumentParser(); ap.add_argument("--quick",action="store_true"); ap.add_argument("--out-dir",default=str(ROOT/"out"/"nitride_cavity")); ap.add_argument("--max-evaluations",type=int,default=5000); ap.add_argument("--dry-run",action="store_true"); a=ap.parse_args(argv)
    grid=build_grid(a.quick); headline=list(itertools.product(*(grid[k] for k in ("height_nm","radius_nm","x_in","T_hs","Q","current_uA","regime"))))
    # Auxiliary inputs are reference slices, never headline denominator rows.
    aux=[]
    for name, vals in (("detuning_offset_meV",[10.]),("background_tau_ns",[1.]),("q167",[167.])):
        for T, reg in itertools.product([230.,300.], HEAD["regime"]): aux.append((name,T,reg,vals[0]))
    islands=list(itertools.product([.5,1.,5.],[230.,250.,273.,300.]))
    planned=len(headline)+len(aux)+len(islands)+1
    if planned>a.max_evaluations: raise ValueError("planned calls exceed --max-evaluations")
    if a.dry_run:
        print(f"headline_rows={len(headline)} planned_rows={planned} unique_evaluate_calls<={planned}"); return 0
    out=_safe_out(a.out_dir); out.mkdir(parents=True,exist_ok=True); os.environ["MPLCONFIGDIR"]=str(out/"mplconfig")
    start=time.time(); cache={}; rows=[]
    for i, vals in enumerate(headline):
        p=dict(zip(("height_nm","radius_nm","x_in","T_hs","Q","current_uA","regime"),vals)); reg=p.pop("regime")
        rows.append(_row(f"H{i:04d}","headline",reg,p,cache))
    sens=[]
    for i,(name,T,reg,val) in enumerate(aux):
        p=dict(REF,T_hs=T)
        if name=="detuning_offset_meV": p["detuning_offset_meV"]=val
        elif name=="background_tau_ns": p["background_tau_ns"]=val
        else: p["Q"]=val; p["purcell_enabled"]=False
        sens.append(_row(f"S{i:03d}","sensitivity",reg,p,cache))
    for radius,T in islands:
        p=dict(REF,T_hs=T); r=_row(f"I{len(sens):03d}","island", "deterministic_pair",p,cache); r["island_radius_nm_input"]=radius; sens.append(r)
    comp=_row("C000","comparison","rectangular",dict(REF,T_hs=300.),cache); comp["measured_g2"]=.29; comp["comparison_note"]="Deshpande et al., Applied Physics Letters 105, 141109 (2014), DOI 10.1063/1.4897640, abstract-only [V]; unfitted, non-gating transfer comparison; CONDITIONS INCOMPLETE"
    _write_csv(out/"sweep.csv",rows); _write_csv(out/"sensitivities.csv",sens); _write_csv(out/"deshpande_comparison.csv",[comp])
    fixed=lambda r: r["radius_nm"]==10 and r["x_in"]==.25 and r["current_uA"]==.02
    height_rows=[r for r in rows if fixed(r) and r["Q"]==2000 and r["T_hs"] in (230.,300.)]
    temp_rows=[r for r in rows if fixed(r) and r["height_nm"]==3 and r["Q"]==2000]
    q_rows=[r for r in rows if fixed(r) and r["height_nm"]==3 and r["T_hs"] in (230.,300.)]
    _plot(out/"height_response.png",height_rows,"height_nm","height response")
    _plot(out/"temperature_response.png",temp_rows,"T_hs","temperature response")
    _plot(out/"cavity_q_response.png",q_rows,"Q","cavity Q response")
    _plot(out/"pulse_vs_set.png",temp_rows,"T_hs","pulse versus idealized SET")
    _plot(out/"envelope_pulse.png",[r for r in rows if r["regime"]=="rectangular"],"height_nm","pulse envelope")
    _plot(out/"envelope_set.png",[r for r in rows if r["regime"]=="deterministic_pair"],"height_nm","idealized SET envelope")
    _plot(out/"set_feasibility.png",[r for r in sens if r["row_kind"]=="island"],"T_hs","SET feasibility")
    verdict=compute_verdict(rows,complete=not a.quick)
    lines=["# Nitride cavity sweep results", "", "Assumed planar InGaN/GaN model. [A] endpoints; no device demonstration.", ""]
    for reg,v in verdict.items():
        fields=" ".join(f"{k}={('none' if v[k] is None else v[k])}" for k in ("regime","complete","g2_min","g2_median_eligible","diag_g2_min","flux_max","eligible_fraction","eligible","flux_floor_excluded","coverage_over_eligible","headline_coverage","idealized_pass_count","hardware_infeasible_count","hardware_pass_count","evidence","conditional","T_pass_min"))
        lines += [f"VERDICT: {v['status']} {fields}", ""]
    lines += ["## Limitations", "Q/V realization limits, QCSE/oscillator uncertainty, omitted field-assisted tunnelling, optimistic nonradiative/background assumptions, and SET pair-delivery/hardware constraints remain. Idealized optical feasibility is not hardware demonstration."]
    (out/"results.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
    hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir()
            if p.is_file() and p.name != "manifest.json" and p.suffix != ".log"}
    manifest={"quick":a.quick,"complete":not a.quick,"axes":grid,"requested_headline_rows":len(headline),"completed_headline_rows":len(rows),"evaluate_calls":len(cache),"runtime_s":time.time()-start,"resolved_out_dir":str(out),"versions":{"python":platform.python_version()},"output_hashes":hashes,"card_hashes":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/"cards").glob("nitride-*.yaml")},"plot_row_mapping":{"height_response.png":"sweep.csv fixed-reference rows","temperature_response.png":"sweep.csv fixed-reference rows","cavity_q_response.png":"sweep.csv fixed-reference rows","envelope_pulse.png":"sweep.csv rectangular","envelope_set.png":"sweep.csv deterministic_pair","pulse_vs_set.png":"sweep.csv fixed-reference rows","set_feasibility.png":"sensitivities.csv island rows"}}
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True),encoding="utf-8")
    print(f"generated {len(rows)} headline rows, {len(cache)} evaluate calls in {manifest['runtime_s']:.1f}s")
    return 0
if __name__ == "__main__": raise SystemExit(main())
