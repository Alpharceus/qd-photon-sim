"""Numerical verification and source-transcription checks for the nitride sweep.

This checker independently replays evaluator rows; it is not a held-out
prediction.  The 0.29 literal is a source transcription (Deshpande et al.,
APL 105, 141109 (2014), DOI 10.1063/1.4897640, abstract-only [V]).
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from scripts import run_nitride_cavity as sweep

def _rows(p):
    with Path(p).open(encoding="utf-8", newline="") as f: return list(csv.DictReader(f))
def _num(x):
    try: return float(x)
    except (TypeError,ValueError): return float("nan")
def _bool(x): return str(x).lower()=="true"
def _close(a,b):
    return (math.isnan(a) and math.isnan(b)) or math.isclose(a,b,rel_tol=1e-8,abs_tol=1e-10)

def main(argv=None):
    ap=argparse.ArgumentParser(); ap.add_argument("--out-dir",default="out/nitride_cavity"); ap.add_argument("--quick",action="store_true"); a=ap.parse_args(argv)
    out=Path(a.out_dir); rows=_rows(out/"sweep.csv"); sens=_rows(out/"sensitivities.csv"); checks=[]
    grid=sweep.build_grid(a.quick); expected=1
    for v in grid.values(): expected*=len(v)
    checks += [len(rows)==expected, len({r["row_id"] for r in rows})==len(rows),
               {r["regime"] for r in rows}==set(sweep.HEAD["regime"])]
    for key, vals in grid.items():
        if key=="regime": continue
        checks.append({float(r[key]) for r in rows}==set(vals))
    # strict g2 and inclusive brightness gates, including invalid and independent-extrema fixtures.
    fixture=[{"row_kind":"headline","regime":"rectangular","valid":True,"g2":.49,"signal_flux_s":900.,"eligible":False,"optical_pass":False,"device_pass":False}, {"row_kind":"headline","regime":"rectangular","valid":True,"g2":.6,"signal_flux_s":1001.,"eligible":True,"optical_pass":False,"device_pass":False}]
    st=sweep.compute_stats(fixture); checks += [st["eligible"]==1, st["hardware_pass_count"]==0]
    setfail=[dict(fixture[0], g2=.4, signal_flux_s=1000., eligible=True, optical_pass=True, hardware_feasible=False, device_pass=False, regime="deterministic_pair")]
    checks.append(sweep.compute_verdict(setfail,complete=True)["deterministic_pair"]["status"]=="FAIL")
    checks.append(sweep.compute_verdict([],complete=False)["rectangular"]["complete"] is False)
    # Policy/source-transcription checks: rows are retained even when a model
    # happens to make every headline point valid; this is not a tautology.
    checks += [all("invalid_reasons" in r for r in rows), len(sens)>=12,
               abs(_num(_rows(out/"deshpande_comparison.csv")[0]["measured_g2"])-.29)<1e-15]
    islands=[r for r in sens if r.get("row_kind")=="island"]
    checks += [len(islands)==12, all(abs(_num(r["island_radius_nm"])-_num(r["set_radius_nm"]))<1e-12 for r in islands),
               all(_num(r["set_E_C_meV"])>0 for r in islands)]
    # Rectangular rows must not masquerade as SET-screened hardware data.
    checks.append(all(r.get("set_feasible","") in ("",None) for r in rows if r["regime"]=="rectangular"))
    # Numerical replay: six distinct real rows, same evaluator and exact T grid.
    for r in rows[:6]:
        p={k:_num(r[k]) for k in ("height_nm","radius_nm","x_in","T_hs","Q","current_uA")}
        d=sweep._design(r["regime"],p); got=sweep.device.evaluate(d,T_grid=[p["T_hs"]])["scalars"]
        checks.append(_close(_num(r["g2"]),float(got["g2_op"])))
        checks.append(_close(_num(r["signal_flux_s"]),float(got["collected_flux_pulsed_s"])))
    man=json.loads((out/"manifest.json").read_text(encoding="utf-8"))
    checks += [man["completed_headline_rows"]==expected, man["evaluate_calls"]<= (200 if a.quick else 5000)]
    for name,digest in man["output_hashes"].items(): checks.append(hashlib.sha256((out/name).read_bytes()).hexdigest()==digest)
    comp=_rows(out/"deshpande_comparison.csv")[0]
    checks.append(comp.get("comparison_card_hash")==man["card_hashes"].get("nitride-deshpande2014-comparison-design.yaml"))
    if not a.quick:
        checks += [all(any(r.get(key,"")!="" for r in sens) for key in sweep.SENS),
                   any(r.get("screening_fraction")=="1.0" and r["regime"]=="deterministic_pair" for r in sens)]
    for name in ("height_response.png","temperature_response.png","cavity_q_response.png","envelope_pulse.png","envelope_set.png","pulse_vs_set.png","set_feasibility.png"):
        checks.append((out/name).is_file() and (out/name).stat().st_size>1000)
    print("numerical verification: evaluator replay [A rtol=1e-8, atol=1e-10]")
    print("source transcription: Deshpande 2014 g2=0.29 [V abstract-only]")
    print("non-gating model comparison: no held-out prediction claim")
    print(f"{sum(checks)}/{len(checks)} nitride sweep checks passed")
    return 0 if all(checks) else 1
if __name__=="__main__": raise SystemExit(main())
