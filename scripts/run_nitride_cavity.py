"""Evaluator-only planar InGaN cavity sweep; grid endpoints are [A].

Q=167 is Taylor et al., Nanoscale Research Letters 5, 608 (2010) [V];
its planar-device transfer is [A].  All output physics comes from evaluate.
"""
from __future__ import annotations
import argparse,csv,hashlib,itertools,json,math,os,platform,sys,time
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1];sys.path.insert(0,str(ROOT))
import fsim_core.device as device
G2,FLUX=.5,1000. # [A] brief gates
HEAD={"height_nm":[1.,2.,3.,4.,5.],"radius_nm":[5.,10.,15.],"x_in":[.15,.25,.4],"T_hs":[230.,250.,273.,300.],"Q":[500.,2000.,10000.],"current_uA":[.002,.02,.2],"regime":["rectangular","deterministic_pair"]}
REF={"height_nm":3.,"radius_nm":10.,"x_in":.25,"T_hs":300.,"Q":2000.,"current_uA":.02}
SENS={"screening_fraction":[0.,.5,1.],"strain_fraction":[.5,1.],"gamma300":[25.,35.,55.],"delta_xx":[-10.,10.,16.],"tau_cap_ps":[1.,10.,100.],"tau_cap_scales_with_density":[False,True],"tau_rad0_ns":[.2,1.,4.],"k_nr_ns":[0.,1.],"background_tau_ns":[0.,1.],"b_res":[0.,.1,1.],"eta_out":[.01,.1,.5],"mode_volume_norm":[1.,2.,10.],"detuning_offset_meV":[0.,10.],"purcell_enabled":[False,True]}
def build_grid(quick=False):
 return {"height_nm":[1.,3.,5.],"radius_nm":[10.],"x_in":[.25],"T_hs":HEAD["T_hs"],"Q":[500.,10000.],"current_uA":[.02],"regime":HEAD["regime"]} if quick else {k:list(v) for k,v in HEAD.items()}
def _finite(x):return isinstance(x,(int,float)) and not isinstance(x,bool) and math.isfinite(x)
def _put(o,k,v):
 if isinstance(o,dict):o[k]=v
 else:setattr(o,k,v)
def _design(reg,p,comparison=False):
 name="nitride-deshpande2014-comparison-design.yaml" if comparison else ("nitride-cavity-set-design.yaml" if reg=="deterministic_pair" else "nitride-cavity-pulse-design.yaml")
 d=device.DeviceDesign.load(ROOT/"cards"/name)
 for k,v in p.items():
  if k in ("height_nm","radius_nm","x_in","screening_fraction","strain_fraction"):_put(d.nitride["dot"],k,v)
  elif k=="T_hs":d.thermal.T_hs=v;_put(d.nitride["cavity"],"T_track",v) # [A] re-track operating temperature
  elif k=="Q":_put(d.nitride["cavity"],k,v)
  elif k in ("mode_volume_norm","eta_out","detuning_offset_meV","purcell_enabled","T_track"):_put(d.nitride["cavity"],k,v)
  elif k in ("tau_rad0_ns","k_nr_ns","background_tau_ns"):_put(d.nitride,k,v)
  elif k=="current_uA":d.drive.I_uA=v
  elif k in ("gamma300","delta_xx"):setattr(d.dot,k,v)
  elif k in ("tau_cap_ps","tau_cap_scales_with_density"):setattr(d.ret,k,v)
  elif k=="b_res":d.drive.b_res=v
  elif k=="island_radius_nm":_put(d.drive.set_params,"radius_nm",v)
 # contract: swept composition must drive transport too
 dot=d.nitride["dot"];_put(d.drive.diode,"x_in",dot["x_in"] if isinstance(dot,dict) else dot.x_in)
 d.drive.cycle_loading=reg;return d
def _row(i,kind,reg,p,cache,comparison=False):
 ident=json.dumps({"regime":reg,"comparison":comparison,"inputs":p},sort_keys=True,separators=(",",":"))
 if ident not in cache:cache[ident]=device.evaluate(_design(reg,p,comparison),T_grid=[p["T_hs"]])["scalars"]
 s=cache[ident];r={"row_id":i,"row_kind":kind,"regime":reg,"cache_identity":ident,**p,"g2":s.get("g2_op"),"signal_flux_s":s.get("collected_flux_pulsed_s"),"valid":bool(s.get("valid")),"invalid_reasons":json.dumps(s.get("invalid_reasons",[])),"provenance":"[A/E/DR] evaluator output; Taylor et al. 2010 Q=167 [V] where used"}
 for k,v in s.items():
  if isinstance(v,(str,int,float,bool)) or v is None:r[k]=v
 r["eligible"]=bool(r["valid"] and _finite(r["g2"]) and _finite(r["signal_flux_s"]) and r["signal_flux_s"]>=FLUX)
 if reg=="deterministic_pair":
  r["one_pair_valid"]=bool(s.get("one_pair_valid"));r["set_feasible"]=bool(s.get("set_feasible"));r["pair_supply_possible"]=bool(s.get("pair_supply_possible"));r["hardware_feasible"]=bool(r["set_feasible"] and r["pair_supply_possible"])
 else:r["one_pair_valid"]=r["set_feasible"]=r["pair_supply_possible"]=r["hardware_feasible"]=None
 r["optical_pass"]=bool(r["eligible"] and r["g2"]<G2 and (reg!="deterministic_pair" or r["one_pair_valid"]))
 r["device_pass"]=bool(r["optical_pass"] and (reg!="deterministic_pair" or r["hardware_feasible"]))
 return r
def compute_stats(rows):
 va=[r for r in rows if r.get("valid")];el=[r for r in rows if r.get("eligible")];gs=sorted(r["g2"] for r in el if _finite(r.get("g2")));di=[r for r in va if _finite(r.get("g2"))];dm=min((r["g2"] for r in di),default=None)
 return {"total":len(rows),"valid":len(va),"eligible":len(el),"g2_min":min(gs) if gs else None,"g2_median_eligible":gs[len(gs)//2] if gs else None,"diag_g2_min":dm,"diag_g2_flux_max":max((r["signal_flux_s"] for r in di if r["g2"]==dm and _finite(r.get("signal_flux_s"))),default=None),"flux_max":max((r["signal_flux_s"] for r in va if _finite(r.get("signal_flux_s"))),default=None),"idealized_pass_count":sum(bool(r.get("optical_pass")) for r in rows),"hardware_infeasible_count":sum(r.get("regime")=="deterministic_pair" and r.get("set_feasible") is False for r in rows),"hardware_pass_count":sum(bool(r.get("device_pass")) for r in rows)}
def compute_verdict(rows,*,complete):
 out={}
 for reg in HEAD["regime"]:
  rs=[r for r in rows if r.get("row_kind")=="headline" and r.get("regime")==reg];s=compute_stats(rs);n=len(rs);status="CONDITIONAL" if complete and s["hardware_pass_count"] else "FAIL"
  out[reg]={**s,"regime":reg,"model":"finite electrical pulse" if reg=="rectangular" else "idealized deterministic SET","complete":complete,"status":status,"flux_margin":s["flux_max"]/FLUX if s["flux_max"] is not None else None,"eligible_fraction":s["eligible"]/n if n else 0,"flux_floor_excluded":n-s["eligible"],"coverage_over_eligible":s["hardware_pass_count"]/s["eligible"] if s["eligible"] else None,"headline_coverage":f"{n}/{n}","evidence":"incomplete","conditional":status=="CONDITIONAL","T_pass_min":min((r["T_hs"] for r in rs if r.get("device_pass")),default=None)}
 return out
def _csv(p,rows):
 with p.open("w",newline="",encoding="utf-8") as f:
  w=csv.DictWriter(f,fieldnames=sorted({k for r in rows for k in r}));w.writeheader();w.writerows(rows)
def _plots(out,rows,sens):
 import matplotlib.pyplot as plt
 def xy(name,q,x,title):
  fig,a=plt.subplots(figsize=(7,4));b=a.twinx()
  for reg,m in (("rectangular","o"),("deterministic_pair","s")):
   for t in sorted({r["T_hs"] for r in q if r["regime"]==reg}):
    z=sorted([r for r in q if r["regime"]==reg and r["T_hs"]==t],key=lambda r:r[x])
    if z:a.plot([r[x] for r in z],[r["g2"] for r in z],m+"-",label=f"{reg} g2 {t:g}K");b.plot([r[x] for r in z],[max(r["signal_flux_s"],1e-12) for r in z],m+"--")
  a.axhline(G2,color="k");b.axhline(FLUX,color="k");b.set_yscale("log");a.set(xlabel=x,ylabel="g2",title=title+"; assumed-planar model");b.set_ylabel("first-lens useful flux / s");fig.tight_layout();fig.savefig(out/name,dpi=130);plt.close(fig)
 fixed=lambda r:r["radius_nm"]==10 and r["x_in"]==.25 and r["current_uA"]==.02
 xy("height_response.png",[r for r in rows if fixed(r) and r["Q"]==2000 and r["T_hs"] in (230.,300.)],"height_nm","height response")
 xy("temperature_response.png",[r for r in rows if fixed(r) and r["height_nm"]==3 and r["Q"]==2000],"T_hs","temperature response")
 xy("cavity_q_response.png",[r for r in rows if fixed(r) and r["height_nm"]==3 and r["T_hs"] in (230.,300.)],"Q","cavity Q response")
 xy("pulse_vs_set.png",[r for r in rows if fixed(r) and r["height_nm"]==3 and r["Q"]==2000],"T_hs","pulse versus idealized SET; paired physical inputs")
 for reg,name in (("rectangular","envelope_pulse.png"),("deterministic_pair","envelope_set.png")):
  fig,a=plt.subplots(figsize=(7,4));mat=[];labs=[]
  for t in HEAD["T_hs"]:
   z=[]
   for h in HEAD["height_nm"]:
    q=[r for r in rows if r["regime"]==reg and r["height_nm"]==h and r["T_hs"]==t];z.append(sum(r["optical_pass"] for r in q)/len(q) if q else float("nan"));labs.append((h,t,f"{sum(r['device_pass'] if reg=='deterministic_pair' else r['optical_pass'] for r in q)}/{len(q)}" if q else "n/a"))
   mat.append(z)
  im=a.imshow(mat,vmin=0,vmax=1,aspect="auto");fig.colorbar(im,ax=a,label="same-row optical pass fraction")
  for h,t,label in labs:a.text(HEAD["height_nm"].index(h),HEAD["T_hs"].index(t),label,ha="center",va="center",fontsize=8)
  a.set(xticks=range(5),xticklabels=HEAD["height_nm"],yticks=range(4),yticklabels=HEAD["T_hs"],xlabel="height nm",ylabel="T_hs K",title=("idealized SET optical; numbers hardware-pass/total" if reg!="rectangular" else "finite pulse same-row pass counts"));fig.tight_layout();fig.savefig(out/name,dpi=130);plt.close(fig)
 z=[r for r in sens if r["row_kind"]=="island"];fig,aa=plt.subplots(1,3,figsize=(11,3.5))
 for rad in (.5,1.,5.):
  q=sorted([r for r in z if r["island_radius_nm"]==rad],key=lambda r:r["T_hs"]);x=[r["T_hs"] for r in q];aa[0].plot(x,[r["set_E_C_meV"] for r in q],"o-",label=f"{rad:g} nm");aa[1].plot(x,[rad]*len(x),"o-");aa[1].plot(x,[r["set_radius_max_nm"] for r in q],"x--");aa[2].plot(x,[r["set_f_max_Hz"] for r in q],"o-")
 aa[0].set(title="E_C versus 10 kT",xlabel="T_hs K",ylabel="E_C meV");aa[1].set(title="assumed / allowed island radius",xlabel="T_hs K",ylabel="radius nm");aa[2].set(title="RC rate versus 80 MHz",xlabel="T_hs K",ylabel="f_max Hz",yscale="log");aa[2].axhline(8e7,color="k",ls="--");fig.tight_layout();fig.savefig(out/"set_feasibility.png",dpi=130);plt.close(fig)
def _safe(v):
 b=(ROOT/"out"/"nitride_cavity").resolve();p=Path(v).resolve()
 if p!=b and b not in p.parents:raise ValueError("out-dir must remain under out/nitride_cavity")
 return p
def main(argv=None):
 ap=argparse.ArgumentParser();ap.add_argument("--quick",action="store_true");ap.add_argument("--out-dir",default=str(ROOT/"out"/"nitride_cavity"));ap.add_argument("--max-evaluations",type=int,default=5000);ap.add_argument("--dry-run",action="store_true");a=ap.parse_args(argv);g=build_grid(a.quick);head=[dict(zip(g,z)) for z in itertools.product(*g.values())]
 axes={"detuning_offset_meV":[10.],"background_tau_ns":[1.]} if a.quick else SENS;aux=[(k,v,t,r) for k,vs in axes.items() for v in vs for t in (230.,300.) for r in HEAD["regime"]]+[("q167",167.,t,r) for t in (230.,300.) for r in HEAD["regime"]];islands=[(r,t) for r in (.5,1.,5.) for t in HEAD["T_hs"]];plan=len(head)+len(aux)+len(islands)+1
 if plan>a.max_evaluations:raise ValueError("planned calls exceed --max-evaluations")
 if a.dry_run:print(f"headline_rows={len(head)} planned_rows={plan} unique_evaluate_calls<={plan}");return 0
 out=_safe(a.out_dir);out.mkdir(parents=True,exist_ok=True);os.environ["MPLBACKEND"]="Agg";os.environ["MPLCONFIGDIR"]=str(out/"mplconfig");cache={};start=time.time();rows=[]
 for n,p in enumerate(head):reg=p.pop("regime");rows.append(_row(f"H{n:04d}","headline",reg,p,cache))
 sens=[]
 for n,(k,v,t,reg) in enumerate(aux):
  p={**REF,"T_hs":t,k:v}
  if k=="q167":p.pop(k);p.update(Q=167.,purcell_enabled=False)
  sens.append(_row(f"S{n:03d}","sensitivity",reg,p,cache))
 for rad,t in islands:sens.append(_row(f"I{len(sens):03d}","island","deterministic_pair",{**REF,"T_hs":t,"island_radius_nm":rad},cache))
 comp=_row("C000","comparison","rectangular",{"T_hs":300.},cache,True);ch=hashlib.sha256((ROOT/"cards"/"nitride-deshpande2014-comparison-design.yaml").read_bytes()).hexdigest();comp.update(measured_g2=.29,measured_count_rate="unavailable",comparison_card_hash=ch,comparison_note="Deshpande et al., APL 105, 141109 (2014), DOI 10.1063/1.4897640, abstract-only [V]; CONDITIONS INCOMPLETE")
 _csv(out/"sweep.csv",rows);_csv(out/"sensitivities.csv",sens);_csv(out/"deshpande_comparison.csv",[comp]);_plots(out,rows,sens);vs=compute_verdict(rows,complete=not a.quick)
 lines=["# Nitride cavity sweep results","","Assumed planar model, evaluator outputs only; no held-out prediction.",""]
 keys=("regime","model","complete","g2_min","g2_median_eligible","diag_g2_min","diag_g2_flux_max","flux_max","flux_margin","eligible_fraction","eligible","flux_floor_excluded","coverage_over_eligible","headline_coverage","idealized_pass_count","hardware_infeasible_count","hardware_pass_count","evidence","conditional","T_pass_min")
 for v in vs.values():lines += ["VERDICT: "+v["status"]+" "+" ".join(f"{k}={v[k] if v[k] is not None else 'none'}" for k in keys),""]
 lines += ["## Deshpande comparison",f"Measured g2=0.29 [V abstract-only; Deshpande et al., APL 105, 141109 (2014), DOI 10.1063/1.4897640]. Evaluated card g2={comp['g2']}, flux={comp['signal_flux_s']}/s; count rate unavailable. CONDITIONS INCOMPLETE.","| transfer | value |","|---|---|","| geometry/x | 2 nm / 12.5 nm [A], x=0.40 [V abstract-only] |","| rate/waveform | 200 MHz reported maximum [V], waveform/current [A] |","| field/cavity | planar QCSE, Q/V/outcoupling [A] |","","## Sensitivity","Full mode contains all one-at-a-time contract axes at 230/300 K; quick contains named diagnostics. screening_fraction 0 through 1 SET rows are in sensitivities.csv.","","## Limitations","Q/V realization, oscillator/QCSE uncertainty, omitted field-assisted tunnelling, optimistic nonradiative/background assumptions, SET pair delivery and hardware screen limit this model. Idealized optical pass is not hardware demonstration."]
 (out/"results.md").write_text("\n".join(lines)+"\n",encoding="utf-8");hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file() and p.name!="manifest.json"};man={"quick":a.quick,"complete":not a.quick,"axes":g,"requested_headline_rows":len(head),"completed_headline_rows":len(rows),"invalid_headline_rows":sum(not r["valid"] for r in rows),"evaluate_calls":len(cache),"runtime_s":time.time()-start,"resolved_out_dir":str(out),"versions":{"python":platform.python_version()},"output_hashes":hashes,"card_hashes":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/"cards").glob("nitride-*.yaml")},"source_hashes":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/"fsim_core").glob("*.py")},"plot_row_mapping":{n:("sensitivities.csv island row_id" if n=="set_feasibility.png" else "sweep.csv evaluator row_id") for n in ("height_response.png","temperature_response.png","cavity_q_response.png","envelope_pulse.png","envelope_set.png","pulse_vs_set.png","set_feasibility.png")}}
 (out/"manifest.json").write_text(json.dumps(man,indent=2,sort_keys=True),encoding="utf-8");print(f"generated {len(rows)} headline rows, {len(cache)} evaluate calls in {man['runtime_s']:.1f}s");return 0
if __name__=="__main__":raise SystemExit(main())
