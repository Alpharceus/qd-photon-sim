"""Evaluator-only planar InGaN cavity sweep; grid endpoints are [A].

Q=167 is Taylor et al., Nanoscale Research Letters 5, 608 (2010) [V];
its planar-device transfer is [A].  All output physics comes from evaluate.

Fix round 2 (2026-09-09, Opus re-review "physics closed, presentation
FAIL"): headline rows use the card default screening_fraction=0.0, an
explicit CONSERVATIVE LOWER BOUND (no polarization-field screening).
screening_fraction=1.0 counterpart rows (row_kind='bound') are evaluated
alongside them, at the same inputs otherwise, as the UPPER BOUND -- never
hand-entered, always additional evaluate() calls through the same cache.
Every headline figure plots both (solid=unscreened lower bound,
dashed=screened upper bound) with a legend; results.md carries a bounds
table under the VERDICT lines.
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
# [A] density used for the tau_cap_scales_with_density one-at-a-time axis,
# per this round's binding rule: both False and True run AT the same
# (sparse) 1e9 cm^-2 density so the density-scaling channel is actually
# exercised (it is a no-op at any single density when the flag is False).
DENSITY_AXIS_CM2=1e9
BOUND_SCREEN=1. # [A] upper-bound screening_fraction paired with every unscreened headline figure
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
  elif k=="n_dot_cm2":d.drive.n_dot_cm2=v
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
 else:
  # rectangular rows carry no SET hardware screen; leave these unset
  # (None -> blank CSV cell), never False, so they can never masquerade
  # as hardware-screened SET data.
  r["one_pair_valid"]=r["set_feasible"]=r["pair_supply_possible"]=r["hardware_feasible"]=None
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
  # idealized_status (fix round 2, Opus re-review high finding): the
  # spec-mandated field distinguishing "idealized optical pass, hardware
  # infeasible" from an unconditional device pass or no pass at all.
  if reg=="deterministic_pair":
   idealized_status=("pass_hardware_infeasible" if s["idealized_pass_count"]>0 and s["hardware_pass_count"]==0
                      else "pass_hardware_qualified" if s["hardware_pass_count"]>0
                      else "no_idealized_pass")
  else:
   idealized_status="pass" if s["idealized_pass_count"]>0 else "no_idealized_pass"
  out[reg]={**s,"regime":reg,"model":"finite electrical pulse" if reg=="rectangular" else "idealized deterministic SET","complete":complete,"status":status,"idealized_status":idealized_status,"flux_margin":s["flux_max"]/FLUX if s["flux_max"] is not None else None,"eligible_fraction":s["eligible"]/n if n else 0,"flux_floor_excluded":n-s["eligible"],"coverage_over_eligible":s["hardware_pass_count"]/s["eligible"] if s["eligible"] else None,"headline_coverage":f"{n}/{n}","evidence":"incomplete","conditional":status=="CONDITIONAL","T_pass_min":min((r["T_hs"] for r in rs if r.get("device_pass")),default=None)}
 return out
def _csv(p,rows):
 with p.open("w",newline="",encoding="utf-8") as f:
  w=csv.DictWriter(f,fieldnames=sorted({k for r in rows for k in r}));w.writeheader();w.writerows(rows)
def _bound_rows(grid,cache):
 """Fix round 2: additional evaluate() rows (never hand-entered) that back
 every headline figure's screening bounds and the results.md bounds
 table. row_kind='reference' is the SAME screening_fraction=0.0 lower
 bound as headline (needed as genuine auxiliary calls in --quick, whose
 headline grid omits Q=2000); row_kind='bound' is the screening_fraction
 upper bound (1.0 for the figures, plus 0.5 for the bounds table only).
 All combos are at the fixed reference radius/x/current, matching REF."""
 combos=set()
 for h in grid["height_nm"]:
  for t in (230.,300.):combos.add((h,t,2000.))
 for t in HEAD["T_hs"]:combos.add((3.,t,2000.))
 for t in (230.,300.):
  for q in HEAD["Q"]:combos.add((3.,t,q))
 out=[]
 for h,t,q in sorted(combos):
  for scr in (0.,BOUND_SCREEN):
   for reg in HEAD["regime"]:
    p={"height_nm":h,"radius_nm":10.,"x_in":.25,"current_uA":.02,"T_hs":t,"Q":q,"screening_fraction":scr}
    kind="reference" if scr==0. else "bound"
    row=_row(f"{'R' if kind=='reference' else 'B'}{len(out):03d}",kind,reg,p,cache)
    row["bound_kind"]="unscreened_lower" if scr==0. else "screened_upper"
    out.append(row)
 # bounds-table-only midpoint (screening_fraction=0.5), reference geometry, 230/300 K
 for t in (230.,300.):
  for reg in HEAD["regime"]:
   p={"height_nm":3.,"radius_nm":10.,"x_in":.25,"current_uA":.02,"T_hs":t,"Q":2000.,"screening_fraction":.5}
   row=_row(f"BT{len(out):03d}","bound",reg,p,cache);row["bound_kind"]="screened_mid_0.5";out.append(row)
 return out
def _anchor_rows(cache):
 """Fix round 2 (Opus re-review): a labelled fixed-anchor (T_track=300 K)
 sensitivity row set across the full T_hs axis -- headline rows re-track
 the cavity to T_track=T_hs every row (see _design); these rows instead
 hold T_track fixed at 300 K while T_hs varies, so the resulting
 detuning is a visible, named diagnostic rather than an unlabelled
 artifact of the temperature-response figure."""
 out=[]
 for t in HEAD["T_hs"]:
  for reg in HEAD["regime"]:
   p={**REF,"T_hs":t,"T_track":300.}
   row=_row(f"A{len(out):03d}","anchor",reg,p,cache);row["sensitivity_axis"]="T_track_fixed_300K";out.append(row)
 return out
def _plots(out,rows,sens):
 import matplotlib.pyplot as plt
 from matplotlib.patches import Patch
 rowmap={};legends={}
 def _leg(fig,name,ax,handles=None,outside=False):
  kw=dict(fontsize=7,loc="upper center",bbox_to_anchor=(.5,-.18),ncol=2) if outside else dict(fontsize=7,loc="best")
  lg=ax.legend(handles=handles,**kw) if handles else ax.legend(**kw)
  legends[name]=lg is not None and len(lg.get_texts())>0
 def bound_pair(pool,reg,x,fixed):
  ref=sorted([r for r in pool if r["row_kind"]=="reference" and r["regime"]==reg and fixed(r)],key=lambda r:r[x])
  bnd=sorted([r for r in pool if r["row_kind"]=="bound" and r["bound_kind"]=="screened_upper" and r["regime"]==reg and fixed(r)],key=lambda r:r[x])
  return ref,bnd
 def xy(name,x,fixed,title,pool):
  fig,a=plt.subplots(figsize=(7.5,5.2));b=a.twinx();trace_ids={}
  for reg,m in (("rectangular","o"),("deterministic_pair","s")):
   ref,bnd=bound_pair(pool,reg,x,fixed)
   if ref:
    a.plot([r[x] for r in ref],[r["g2"] for r in ref],m+"-",label=f"{reg} g2 (unscreened lower bound)")
    b.plot([r[x] for r in ref],[max(r["signal_flux_s"],1e-12) for r in ref],m+"--",alpha=.55,label=f"{reg} flux (unscreened lower bound)")
    trace_ids[f"{reg}_unscreened_g2"]=[r["row_id"] for r in ref]
   if bnd:
    a.plot([r[x] for r in bnd],[r["g2"] for r in bnd],m+":",label=f"{reg} g2 (screened upper bound)")
    b.plot([r[x] for r in bnd],[max(r["signal_flux_s"],1e-12) for r in bnd],m+"-.",alpha=.55,label=f"{reg} flux (screened upper bound)")
    trace_ids[f"{reg}_screened_g2"]=[r["row_id"] for r in bnd]
  a.axhline(G2,color="k",lw=.8,label="g2<0.5 optical guide");b.axhline(FLUX,color="gray",lw=.8,ls="--",label="flux>=1000/s guide")
  b.set_yscale("log");a.set(xlabel=x,ylabel="g2");a.set_title(title+"; assumed-planar model\nsolid/dashed=unscreened lower bound, dotted/dash-dot=screening=1 upper bound",fontsize=8);b.set_ylabel("first-lens useful flux / s (log)")
  h1,l1=a.get_legend_handles_labels();h2,l2=b.get_legend_handles_labels();_leg(fig,name,a,h1+h2,outside=True)
  fig.tight_layout();fig.savefig(out/name,dpi=130);plt.close(fig);rowmap[name]=trace_ids
 fixed_h=lambda r:r["radius_nm"]==10 and r["x_in"]==.25 and r["current_uA"]==.02 and r["Q"]==2000 and r["T_hs"] in (230.,300.)
 fixed_t=lambda r:r["radius_nm"]==10 and r["x_in"]==.25 and r["current_uA"]==.02 and r["height_nm"]==3 and r["Q"]==2000
 fixed_q=lambda r:r["radius_nm"]==10 and r["x_in"]==.25 and r["current_uA"]==.02 and r["height_nm"]==3 and r["T_hs"] in (230.,300.)
 xy("height_response.png","height_nm",fixed_h,"height response",sens)
 xy("temperature_response.png","T_hs",fixed_t,"temperature response",sens)
 xy("cavity_q_response.png","Q",fixed_q,"cavity Q response",sens)
 # pulse_vs_set.png: a genuine contrast, not a relabelled temperature plot --
 # grouped bars, same paired rows (unscreened+screened) as temperature_response,
 # g2 and flux side by side per T_hs, with the loading-mechanism difference annotated.
 ref_p,bnd_p=bound_pair(sens,"rectangular","T_hs",fixed_t);ref_s,bnd_s=bound_pair(sens,"deterministic_pair","T_hs",fixed_t)
 fig,(ax1,ax2)=plt.subplots(1,2,figsize=(10,4.2));ts=HEAD["T_hs"];wdt=.18;xs=range(len(ts))
 def _val(rr,t,key):
  m=[r for r in rr if r["T_hs"]==t];return m[0][key] if m else float("nan")
 for i,(rr,lab,off) in enumerate([(ref_p,"pulse unscreened",-1.5),(bnd_p,"pulse screened=1",-.5),(ref_s,"SET unscreened",.5),(bnd_s,"SET screened=1",1.5)]):
  ax1.bar([x+off*wdt for x in xs],[_val(rr,t,"g2") for t in ts],width=wdt,label=lab)
  ax2.bar([x+off*wdt for x in xs],[max(_val(rr,t,"signal_flux_s"),1e-12) for t in ts],width=wdt,label=lab)
 ax1.axhline(G2,color="k",lw=.8);ax1.set(xticks=list(xs),xticklabels=[f"{t:g}" for t in ts],xlabel="T_hs K",ylabel="g2");ax1.set_title("pulse vs idealized SET: g2",fontsize=8)
 ax2.axhline(FLUX,color="gray",lw=.8,ls="--");ax2.set_yscale("log");ax2.set(xticks=list(xs),xticklabels=[f"{t:g}" for t in ts],xlabel="T_hs K",ylabel="flux /s (log)");ax2.set_title("same contrast: collected flux",fontsize=8)
 fig.suptitle("pulse_vs_set: same physical inputs, two loading mechanisms (finite\nelectrical pulse vs idealized deterministic pair)",fontsize=8)
 _leg(fig,"pulse_vs_set.png",ax1);ax2.legend(fontsize=7,loc="best")
 fig.tight_layout();fig.savefig(out/"pulse_vs_set.png",dpi=130);plt.close(fig)
 rowmap["pulse_vs_set.png"]={"pulse_unscreened":[r["row_id"] for r in ref_p],"pulse_screened":[r["row_id"] for r in bnd_p],"set_unscreened":[r["row_id"] for r in ref_s],"set_screened":[r["row_id"] for r in bnd_s]}
 # envelope plots: left panel = full headline aggregate (all other axes,
 # unscreened lower bound, as before); right panel = reference-geometry
 # screened (screening_fraction=1) upper bound, one row per cell.
 refpool={(r["height_nm"],r["T_hs"]):r for r in sens if r["row_kind"]=="bound" and r["bound_kind"]=="screened_upper" and r["radius_nm"]==10 and r["x_in"]==.25 and r["current_uA"]==.02 and r["Q"]==2000}
 for reg,name in (("rectangular","envelope_pulse.png"),("deterministic_pair","envelope_set.png")):
  fig,(a0,a1)=plt.subplots(1,2,figsize=(11.5,4.4));mat=[];labs=[];cellmap={}
  for t in HEAD["T_hs"]:
   zrow=[]
   for h in HEAD["height_nm"]:
    q=[r for r in rows if r["regime"]==reg and r["height_nm"]==h and r["T_hs"]==t]
    zrow.append(sum(r["optical_pass"] for r in q)/len(q) if q else float("nan"))
    idl=sum(r["optical_pass"] for r in q);hw=sum(r["device_pass"] for r in q)
    labs.append((h,t,f"{idl}/{hw}/{len(q)}" if q else "n/a"))
    cellmap[f"h{h:g}_T{t:g}"]=[r["row_id"] for r in q]
   mat.append(zrow)
  im=a0.imshow(mat,vmin=0,vmax=1,aspect="auto",origin="lower");fig.colorbar(im,ax=a0,label="idealized optical-pass fraction")
  for h,t,label in labs:
   xi,yi=HEAD["height_nm"].index(h),HEAD["T_hs"].index(t)
   a0.text(xi,yi,label,ha="center",va="center",fontsize=7)
   if label=="n/a":continue
   idl_n,hw_n,tot_n=(int(v) for v in label.split("/"))
   if reg=="deterministic_pair" and idl_n>0 and hw_n==0:
    a0.add_patch(plt.Rectangle((xi-.5,yi-.5),1,1,fill=False,hatch="///",edgecolor="k",lw=.5))
  a0.set(xticks=range(5),xticklabels=[f"{v:g}" for v in HEAD["height_nm"]],yticks=range(4),yticklabels=[f"{v:g}" for v in HEAD["T_hs"]],xlabel="height nm",ylabel="T_hs K");a0.set_title("unscreened lower bound: idealized/hardware/total\nper cell (all other headline axes)",fontsize=8)
  handles=[Patch(facecolor="none",edgecolor="k",hatch="///",label="idealized pass, hardware infeasible")] if reg=="deterministic_pair" else [Patch(facecolor="none",edgecolor="none",label="text = idealized/hardware/total same-row passes")]
  _leg(fig,name,a0,handles)
  mat2=[];labs2=[]
  for t in HEAD["T_hs"]:
   zrow=[]
   for h in HEAD["height_nm"]:
    r=refpool.get((h,t));reg_r=r if (r and r["regime"]==reg) else None
    if reg_r is None:
     rb=[rr for rr in sens if rr["row_kind"]=="bound" and rr["bound_kind"]=="screened_upper" and rr["height_nm"]==h and rr["T_hs"]==t and rr["regime"]==reg and rr["radius_nm"]==10 and rr["x_in"]==.25 and rr["current_uA"]==.02 and rr["Q"]==2000]
     reg_r=rb[0] if rb else None
    zrow.append(1. if reg_r and reg_r["optical_pass"] else (0. if reg_r else float("nan")))
    labs2.append((h,t,("pass" if reg_r and reg_r["device_pass"] else ("opt" if reg_r and reg_r["optical_pass"] else "fail")) if reg_r else "n/a"))
   mat2.append(zrow)
  im2=a1.imshow(mat2,vmin=0,vmax=1,aspect="auto",origin="lower",cmap="RdYlGn");fig.colorbar(im2,ax=a1,label="screened (screening_fraction=1) optical pass")
  for h,t,label in labs2:
   xi,yi=HEAD["height_nm"].index(h),HEAD["T_hs"].index(t);a1.text(xi,yi,label,ha="center",va="center",fontsize=7)
   if reg=="deterministic_pair" and label=="opt":a1.add_patch(plt.Rectangle((xi-.5,yi-.5),1,1,fill=False,hatch="///",edgecolor="k",lw=.5))
  a1.set(xticks=range(5),xticklabels=[f"{v:g}" for v in HEAD["height_nm"]],yticks=range(4),yticklabels=[f"{v:g}" for v in HEAD["T_hs"]],xlabel="height nm",ylabel="T_hs K");a1.set_title("screened upper bound (reference geometry, 1 row/cell)\n'opt'=idealized pass, hardware infeasible",fontsize=8)
  fig.tight_layout();fig.savefig(out/name,dpi=130);plt.close(fig)
  rowmap[name]={"unscreened_cells":cellmap,"screened_cells":{f"h{h:g}_T{t:g}":([refpool[(h,t)]["row_id"]] if (h,t) in refpool else []) for t in HEAD["T_hs"] for h in HEAD["height_nm"]}}
 # set_feasibility.png: E_C vs 10 kT, assumed vs allowed island radius, RC rate vs 80 MHz.
 z=[r for r in sens if r["row_kind"]=="island"];fig,aa=plt.subplots(1,3,figsize=(12,3.8));ids={}
 ts_i=sorted({r["T_hs"] for r in z})
 kt10=[10.*device.KB_SI*t/device.E_SI*1e3 for t in ts_i] # 10 kT margin threshold, meV, vs T_hs
 for rad in (.5,1.,5.):
  q=sorted([r for r in z if r["island_radius_nm"]==rad],key=lambda r:r["T_hs"]);x=[r["T_hs"] for r in q]
  aa[0].plot(x,[r["set_E_C_meV"] for r in q],"o-",label=f"E_C, {rad:g} nm island")
  aa[1].plot(x,[rad]*len(x),"o-",label=f"assumed radius {rad:g} nm");aa[1].plot(x,[r["set_radius_max_nm"] for r in q],"x--",label=f"allowed max ({rad:g} nm island)")
  aa[2].plot(x,[r["set_f_max_Hz"] for r in q],"o-",label=f"f_max, {rad:g} nm island")
  ids[f"{rad:g}nm"]=[r["row_id"] for r in q]
 aa[0].plot(ts_i,kt10,"k--",label="10 kT margin threshold")
 aa[0].set(title="E_C versus 10 kT",xlabel="T_hs K",ylabel="E_C meV");aa[0].legend(fontsize=6,loc="best")
 aa[1].set(title="assumed / allowed island radius",xlabel="T_hs K",ylabel="radius nm");aa[1].legend(fontsize=6,loc="best")
 aa[2].axhline(8e7,color="k",ls="--",label="80 MHz rep rate");aa[2].set(title="RC rate versus 80 MHz",xlabel="T_hs K",ylabel="f_max Hz",yscale="log");aa[2].legend(fontsize=6,loc="best")
 legends["set_feasibility.png"]=True
 fig.tight_layout();fig.savefig(out/"set_feasibility.png",dpi=130);plt.close(fig)
 rowmap["set_feasibility.png"]=ids
 return rowmap,legends
def _safe(v):
 b=(ROOT/"out"/"nitride_cavity").resolve();p=Path(v).resolve()
 if p!=b and b not in p.parents:raise ValueError("out-dir must remain under out/nitride_cavity")
 return p
def _bounds_table_md(sens):
 tbl=[r for r in sens if r["row_kind"] in ("reference","bound") and r["height_nm"]==3 and r["Q"]==2000 and r["T_hs"] in (230.,300.)]
 tbl=sorted(tbl,key=lambda r:(r["regime"],r["T_hs"],r["screening_fraction"]))
 lines=["| regime | T_hs K | screening_fraction | g2_op | flux /s | S_X | E_C/kT | hardware_feasible |","|---|---|---|---|---|---|---|---|"]
 for r in tbl:
  ec=r.get("set_EC_over_kT");hf=r.get("hardware_feasible")
  ec_s='' if ec in (None,'') or (isinstance(ec,float) and math.isnan(ec)) else f"{float(ec):.4g}"
  hf_s='n/a' if hf in (None,'') else str(hf)
  lines.append(f"| {r['regime']} | {r['T_hs']:g} | {r['screening_fraction']:g} | {r['g2']:.6g} | {r['signal_flux_s']:.6g} | {r.get('S_X',float('nan')):.4g} | {ec_s} | {hf_s} |")
 return "\n".join(lines)
def _sensitivity_md(sens):
 axes=sorted({r.get("sensitivity_axis") for r in sens if r.get("row_kind")=="sensitivity" and r.get("sensitivity_axis")})
 lines=["| axis | regime | T_hs K | value | g2_op | flux /s |","|---|---|---|---|---|---|"]
 for ax in axes:
  rs=sorted([r for r in sens if r.get("sensitivity_axis")==ax],key=lambda r:(r["regime"],r["T_hs"],str(r.get(ax))))
  for r in rs:
   lines.append(f"| {ax} | {r['regime']} | {r['T_hs']:g} | {r.get(ax)} | {r['g2']:.6g} | {r['signal_flux_s']:.6g} |")
 return "\n".join(lines)
def _islands_md(sens):
 rs=sorted([r for r in sens if r["row_kind"]=="island"],key=lambda r:(r["island_radius_nm"],r["T_hs"]))
 lines=["| island radius nm | T_hs K | E_C meV | E_C/kT | allowed max radius nm | R_T/R_Q | hardware_feasible | f_max Hz |","|---|---|---|---|---|---|---|---|"]
 for r in rs:
  lines.append(f"| {r['island_radius_nm']:g} | {r['T_hs']:g} | {r['set_E_C_meV']:.4g} | {r['set_EC_over_kT']:.4g} | {r['set_radius_max_nm']:.4g} | {r['set_R_T_over_RQ']:.4g} | {r['hardware_feasible']} | {r['set_f_max_Hz']:.4g} |")
 return "\n".join(lines)
def _anchor_md(sens):
 rs=sorted([r for r in sens if r["row_kind"]=="anchor"],key=lambda r:(r["regime"],r["T_hs"]))
 lines=["| regime | T_hs K | T_track K | detuning_meV | g2_op | flux /s |","|---|---|---|---|---|---|"]
 for r in rs:
  lines.append(f"| {r['regime']} | {r['T_hs']:g} | {r['T_track']:g} | {r.get('detuning_meV',float('nan')):.4g} | {r['g2']:.6g} | {r['signal_flux_s']:.6g} |")
 return "\n".join(lines)
def main(argv=None):
 ap=argparse.ArgumentParser();ap.add_argument("--quick",action="store_true");ap.add_argument("--out-dir",default=str(ROOT/"out"/"nitride_cavity"));ap.add_argument("--max-evaluations",type=int,default=5000);ap.add_argument("--dry-run",action="store_true");a=ap.parse_args(argv);g=build_grid(a.quick);head=[dict(zip(g,z)) for z in itertools.product(*g.values())]
 axes={"detuning_offset_meV":[10.],"background_tau_ns":[1.]} if a.quick else SENS;aux_combos=[(k,v,t,r) for k,vs in axes.items() for v in vs for t in (230.,300.) for r in HEAD["regime"]]
 q167_combos=[("Q",167.,t,r) for t in (230.,300.) for r in HEAD["regime"]]+[("purcell_enabled",False,t,r) for t in (230.,300.) for r in HEAD["regime"]]
 islands=[(r,t) for r in (.5,1.,5.) for t in HEAD["T_hs"]]
 bound_combo_ct=(len(g["height_nm"])*2+len(HEAD["T_hs"])+2*len(HEAD["Q"]))*2*len(HEAD["regime"])+4
 anchor_ct=len(HEAD["T_hs"])*len(HEAD["regime"])
 plan=len(head)+len(aux_combos)+len(q167_combos)+len(islands)+1+bound_combo_ct+anchor_ct
 if plan>a.max_evaluations:raise ValueError("planned calls exceed --max-evaluations")
 if a.dry_run:print(f"headline_rows={len(head)} planned_rows={plan} unique_evaluate_calls<={plan}");return 0
 out=_safe(a.out_dir);out.mkdir(parents=True,exist_ok=True);os.environ["MPLBACKEND"]="Agg";os.environ["MPLCONFIGDIR"]=str(out/"mplconfig");cache={};start=time.time();rows=[]
 for n,p in enumerate(head):reg=p.pop("regime");rows.append(_row(f"H{n:04d}","headline",reg,p,cache))
 sens=[]
 for k,v,t,reg in aux_combos:
  p={**REF,"T_hs":t,k:v}
  if k=="tau_cap_scales_with_density":p["n_dot_cm2"]=DENSITY_AXIS_CM2
  row=_row(f"S{len(sens):03d}","sensitivity",reg,p,cache);row["sensitivity_axis"]=k;sens.append(row)
 for k,v,t,reg in q167_combos:
  p={**REF,"T_hs":t,k:v};row=_row(f"S{len(sens):03d}","sensitivity",reg,p,cache);row["sensitivity_axis"]=k;sens.append(row)
 for rad,t in islands:
  row=_row(f"I{len(sens):03d}","island","deterministic_pair",{**REF,"T_hs":t,"island_radius_nm":rad},cache);row["sensitivity_axis"]="island_radius_nm";sens.append(row)
 sens += _anchor_rows(cache)
 sens += _bound_rows(g,cache)
 comp=_row("C000","comparison","rectangular",{"T_hs":300.},cache,True);ch=hashlib.sha256((ROOT/"cards"/"nitride-deshpande2014-comparison-design.yaml").read_bytes()).hexdigest();comp.update(measured_g2=.29,measured_count_rate="unavailable",comparison_card_hash=ch,comparison_note="Deshpande et al., APL 105, 141109 (2014), DOI 10.1063/1.4897640, abstract-only [V]; CONDITIONS INCOMPLETE")
 _csv(out/"sweep.csv",rows);_csv(out/"sensitivities.csv",sens);_csv(out/"deshpande_comparison.csv",[comp]);rowmap,legends=_plots(out,rows,sens);vs=compute_verdict(rows,complete=not a.quick)
 lines=["# Nitride cavity sweep results","","Assumed planar model, evaluator outputs only; no held-out prediction.",""]
 keys=("regime","model","complete","g2_min","g2_median_eligible","diag_g2_min","diag_g2_flux_max","flux_max","flux_margin","eligible_fraction","eligible","flux_floor_excluded","coverage_over_eligible","headline_coverage","idealized_pass_count","idealized_status","hardware_infeasible_count","hardware_pass_count","evidence","conditional","T_pass_min")
 for v in vs.values():lines += ["VERDICT: "+v["status"]+" "+" ".join(f"{k}={v[k] if v[k] is not None else 'none'}" for k in keys),""]
 set_diag=vs["deterministic_pair"]["diag_g2_min"]
 if set_diag is not None and set_diag<1e-6:
  set_note="Note: SET g2_op=0.0 wherever it appears is STRUCTURAL under this model's idealized deterministic one-pair loading (exact-one-pair counting has zero coincidence probability by construction), not evidence of any device suppressing multi-photon emission; it carries no additional device information beyond one_pair_valid/eligible."
 else:
  set_note=("Note: with drive.b_res wired into this branch (fix round 2), SET g2_op is NOT exactly 0: under ideal one-pair loading the exact-one-pair counting result is 0 by construction (cnt_g2=0), so g2_op=1-rho^2 depends only on rho=signal/(signal+background); since bg_counts includes b_res*(collected X counts) and one-pair loading has negligible XX, signal~=collected X counts, so rho asymptotes to 1/(1+b_res) and g2_op asymptotes to 1-(1/(1+b_res))^2 nearly independent of absolute flux (b_res=0.1 on the shipped cards -> g2_op~0.174, matching the bounds table below across screening/T/Q). This is a STRUCTURAL floor set by the assumed residual background channel, not a demonstrated device number.")
 lines += ["## Bounds (unscreened lower bound / screening_fraction=1 upper bound)",
           "Headline rows use the card default screening_fraction=0.0 -- an explicit CONSERVATIVE LOWER BOUND (no polarization-field screening). screening_fraction=1.0 rows below are the UPPER BOUND, evaluated the same way (row_kind='bound' in sensitivities.csv), never hand-entered.",
           "",_bounds_table_md(sens),"",
           set_note,
           "",
           "Cavity re-tuning disclosure: every headline row re-tracks the cavity resonance to T_track=T_hs (the SAME dot's own E_X at that operating point) -- headline results assume a cavity re-tuned per dot/temperature, not a single fixed cavity swept across all conditions. The fixed-anchor (T_track=300 K) sensitivity set below holds the cavity fixed while T_hs varies, to show the resulting detuning as a labelled, separate effect.",
           "",
           "### Fixed-anchor (T_track=300 K) sensitivity","",_anchor_md(sens),"",
           "### SET island-radius sensitivity (classical charging-energy screen only; not demonstrated feasible manufacture)","",_islands_md(sens),""]
 lines += ["## Deshpande comparison",f"Measured g2=0.29 [V abstract-only; Deshpande et al., APL 105, 141109 (2014), DOI 10.1063/1.4897640]. Evaluated card g2={comp['g2']}, flux={comp['signal_flux_s']}/s; count rate unavailable. CONDITIONS INCOMPLETE.","| transfer | value |","|---|---|","| geometry/x | 2 nm / 12.5 nm [A], x=0.40 [V abstract-only] |","| rate/waveform | 200 MHz reported maximum [V], waveform/current [A] |","| field/cavity | planar QCSE, Q/V/outcoupling [A] |",""]
 lines += ["## Sensitivity","Full mode contains all one-at-a-time contract axes at 230/300 K, both regimes (quick contains the named diagnostic subset plus Q=167 and purcell_enabled=False as two independent one-at-a-time rows). Values below are the actual evaluated g2/flux for each axis value.","",_sensitivity_md(sens),""]
 lines += ["## Limitations","Q/V realization, oscillator/QCSE uncertainty, omitted field-assisted tunnelling, optimistic nonradiative/background assumptions, SET pair delivery and hardware screen limit this model. Idealized optical pass is not hardware demonstration.","Reservoir-vs-diode-SRH-layer disagreement: the background reservoir energy (fsim_core.device._nitride_reservoir_energy_eV) is read from nitride.dot's own wetting-layer thickness (0.0 nm on shipped cards -> GaN barrier edge), while the diode's SRH background (fsim_core.nitride_transport.evaluate_injection) lives in drive.diode's separate 0.5 nm InGaN layer -- the two blocks disagree about which material hosts the background carriers. No numerical consequence today (the reservoir energy only sets the cavity/slit spectral ACCEPTANCE of the SRH rate, not the rate itself), but the two should eventually be unified."]
 (out/"results.md").write_text("\n".join(lines)+"\n",encoding="utf-8")
 hashes={p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in out.iterdir() if p.is_file() and p.name!="manifest.json"}
 man={"quick":a.quick,"complete":not a.quick,"axes":g,"requested_headline_rows":len(head),"completed_headline_rows":len(rows),"invalid_headline_rows":sum(not r["valid"] for r in rows),"evaluate_calls":len(cache),"runtime_s":time.time()-start,"resolved_out_dir":str(out),"versions":{"python":platform.python_version()},"output_hashes":hashes,"card_hashes":{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/"cards").glob("nitride-*.yaml")},"source_hashes":{str(p.relative_to(ROOT)):hashlib.sha256(p.read_bytes()).hexdigest() for p in (ROOT/"fsim_core").glob("*.py")},"plot_row_mapping":rowmap,"figure_legends":legends}
 (out/"manifest.json").write_text(json.dumps(man,indent=2,sort_keys=True),encoding="utf-8");print(f"generated {len(rows)} headline rows, {len(cache)} evaluate calls in {man['runtime_s']:.1f}s");return 0
if __name__=="__main__":raise SystemExit(main())
