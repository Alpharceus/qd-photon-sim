"""Bias conversion and conditional Stark diagnostics for planar nitride traces.

Depletion equations use [DR Sze & Ng, Physics of Semiconductor Devices, 3rd
ed. (2007)].  ``field_polarity=+1`` adds the depletion field to the evaluator
static field; in the c-plane reference convention it is antiparallel to the
intrinsic polarization field and normally gives a negative exciton slope [A].
Zhang's value is a non-gating transcription [V Zhang et al., Appl. Phys.
Lett. 108, 153102 (2016), Fig. 5].
"""
from __future__ import annotations
import math
import numbers

ZHANG2016_SLOPE_MEV_PER_V = -10.0
ZHANG2016_COMPARISON = {"slope_meV_per_V": -10.0,
 "citation": "Zhang et al., Appl. Phys. Lett. 108, 153102 (2016), Fig. 5",
 "url": "https://arxiv.org/abs/1602.02325", "measurement": "bias-dependent PL below 2 V",
 "temperature": "unknown in supplied digest", "excitation": "unknown in supplied digest",
 "use": "non-gating comparison; not a planar-card fitted target"}

def _finite(x):
    return isinstance(x, numbers.Real) and not isinstance(x, bool) and math.isfinite(float(x))

def resolve_bias(diode, *, T_j_K, current_uA=None, junction_voltage_V=None, field_polarity=1, external_field_kVcm=0.0):
    """Resolve exactly one forward control; flat-band clipping is inherited.

    Negative junction bias/reverse leakage and breakdown are out of scope.
    """
    why=[]
    if (current_uA is None) == (junction_voltage_V is None): why.append("supply exactly one of current_uA or junction_voltage_V")
    if not _finite(T_j_K) or T_j_K <= 0: why.append("T_j_K must be finite and positive")
    if isinstance(field_polarity,bool) or not isinstance(field_polarity,numbers.Integral) or int(field_polarity) not in (-1,1): why.append("field_polarity must be integer +1 or -1")
    if not _finite(external_field_kVcm): why.append("external_field_kVcm must be finite")
    if current_uA is not None and (not _finite(current_uA) or current_uA < 0): why.append("current_uA must be finite and nonnegative")
    if junction_voltage_V is not None and (not _finite(junction_voltage_V) or junction_voltage_V < 0): why.append("junction_voltage_V must be finite and nonnegative; reverse bias unsupported")
    out={"T_j_K":T_j_K,"current_uA":current_uA,"V_j":junction_voltage_V,"V_terminal":float("nan"),"diode_field_kVcm":float("nan"),"applied_field_kVcm":float("nan"),"flat_band":False,"depletion_regime":None,"bias_valid":False,"invalid_reasons":why}
    if why: return out
    try:
        if current_uA is not None:
            terminal,vj=diode.v_of_i(float(current_uA)*1e-6,float(T_j_K)); out["V_j"],out["V_terminal"]=vj,terminal
        else:
            vj=float(junction_voltage_V); current_A=diode.j_of_vj(vj,float(T_j_K))*diode.area_cm2
            out["current_uA"]=current_A*1e6; out["V_terminal"]=vj+current_A*diode.R_s_ohm
        if not all(_finite(out[k]) for k in ("V_j","V_terminal","current_uA")): raise OverflowError("unsupported operating point")
        dep=diode.depletion(float(out["V_j"]),float(T_j_K)); field=dep.F_kVcm
        if not _finite(field): raise OverflowError("depletion field is nonfinite")
        flat=bool(getattr(dep,"flat_band",False))
        out.update(diode_field_kVcm=float(field),applied_field_kVcm=int(field_polarity)*float(field)+float(external_field_kVcm),flat_band=flat,depletion_regime="flat_band" if flat else "depleted",bias_valid=True)
    except (ArithmeticError,OverflowError,ValueError) as exc: out["invalid_reasons"].append("unsupported operating point: %s" % exc)
    return out

def _rid(row,i): return row.get("row_id",row.get("id",i))
def _marker(row): return tuple((k,row.get(k)) for k in ("depletion_regime","depletion_clipped","flat_band","trace_branch") if k in row)
def _d3(x,y,n):
    return sum(y[j]*(2*x[n]-x[a]-x[b])/((x[j]-x[a])*(x[j]-x[b])) for j in range(3) for a,b in [[k for k in range(3) if k!=j]])

def stark_derivatives(rows, *, voltage_key="V_j", energy_key="E_X_eV", valid_key="spectroscopy_valid"):
    """Same-order adjacent, three-point, nonuniform derivatives.

    Invalid rows split the trace. Valid finite points must be strictly
    increasing; no derivative crosses invalid rows or a regime/flat-band kink.
    """
    rows=list(rows); valid=[float(r[voltage_key]) for r in rows if bool(r.get(valid_key,False)) and _finite(r.get(energy_key)) and _finite(r.get(voltage_key))]
    if any(valid[n+1]<=valid[n] for n in range(len(valid)-1)): raise ValueError("trace valid voltages must be strictly increasing and distinct")
    ans=[]
    for i,row in enumerate(rows):
        r=dict(row); r.update(dE_X_dV_meV_per_V=float("nan"),derivative_valid=False,derivative_row_ids=[],derivative_policy="adjacent_three_point_nonuniform")
        own=bool(row.get(valid_key,False)) and _finite(row.get(energy_key)) and _finite(row.get(voltage_key))
        ix=(0,1,2) if i==0 else (len(rows)-3,len(rows)-2,len(rows)-1) if i==len(rows)-1 else (i-1,i,i+1)
        chosen=[rows[j] for j in ix]
        good=own and len(rows)>=3 and all(bool(q.get(valid_key,False)) and _finite(q.get(energy_key)) and _finite(q.get(voltage_key)) for q in chosen)
        if good and len({_marker(q) for q in chosen})==1:
            x=[float(q[voltage_key]) for q in chosen]; y=[float(q[energy_key]) for q in chosen]
            r.update(dE_X_dV_meV_per_V=1000*_d3(x,y,ix.index(i)),derivative_valid=True,derivative_row_ids=[_rid(rows[j],j) for j in ix])
        ans.append(r)
    return ans

# Explicit fixed-input whitelist: outputs never create separate one-row traces.
_TRACE_KEYS=("geometry","orientation","shape","composition","surrounding_well","regime","temperature_mode","controlled_temperature","field_polarity","applied_static_field_kVcm","external_field_kVcm","T_hs","cavity_reference_V_j_V","cavity_reference_settings","height_nm","radius_nm","x_in","strain_fraction","strain_c_fraction","wl_thickness_nm")
def _freeze(x):
    if isinstance(x,dict): return tuple(sorted((k,_freeze(v)) for k,v in x.items()))
    if isinstance(x,(list,tuple)): return tuple(_freeze(v) for v in x)
    return x
def _gkey(r): return tuple((k,_freeze(r[k])) for k in _TRACE_KEYS if k in r)
def _screen(r): return r.get("screening_fraction",r.get("screening"))
def _vkey(r):
    for k in ("V_j","V_j_V","junction_voltage_V"):
        if k in r:return r[k],k
    return None,None
def _fit(items,key):
    x=[float(r[key]) for _,r in items]; y=[1000*float(r["E_X_eV"]) for _,r in items]; xm=sum(x)/len(x); ym=sum(y)/len(y)
    slope=sum((a-xm)*(b-ym) for a,b in zip(x,y))/sum((a-xm)**2 for a in x); intercept=ym-slope*xm
    return slope,math.sqrt(sum((b-(intercept+slope*a))**2 for a,b in zip(x,y))/len(x))

def screening_compatibility(rows, *, slope_range_meV_per_V, voltage_window_V, lifetime_range_ns=None):
    """Fit contiguous valid same-window slopes per supplied screening hypothesis."""
    try: lo,hi=slope_range_meV_per_V; vlo,vhi=voltage_window_V
    except (TypeError,ValueError): raise ValueError("ranges must each contain two values")
    if not all(_finite(x) for x in (lo,hi,vlo,vhi)) or lo>hi or vlo>=vhi: raise ValueError("invalid slope range or voltage window")
    if lifetime_range_ns is not None:
        if set(lifetime_range_ns)!={"voltage_V","min_ns","max_ns","kind"} or lifetime_range_ns["kind"] not in ("bare","cavity"): raise ValueError("lifetime_range_ns must exactly specify voltage_V, min_ns, max_ns, kind")
        if not all(_finite(lifetime_range_ns[k]) and lifetime_range_ns[k]>0 for k in ("voltage_V","min_ns","max_ns")) or lifetime_range_ns["min_ns"]>lifetime_range_ns["max_ns"]: raise ValueError("invalid lifetime range")
    grouped={}
    for i,r in enumerate(rows): grouped.setdefault((_gkey(r),_screen(r)),[]).append((i,r))
    result=[]
    for (group,hyp),numbered in grouped.items():
        branches=[]; branch=[]; direct=0
        for item in numbered:
            val,_=_vkey(item[1]); good=bool(item[1].get("spectroscopy_valid",item[1].get("valid",False))) and _finite(val) and _finite(item[1].get("E_X_eV"))
            if not good:
                if branch: branches.append(branch)
                branch=[]; direct=0; continue
            if branch:
                prev,_=_vkey(branch[-1][1]); delta=float(val)-float(prev); turning=delta==0 or (direct and (delta>0)!=(direct>0)) or _marker(item[1])!=_marker(branch[-1][1])
                if turning: branches.append(branch); branch=[]; direct=0
                elif branch: direct=1 if delta>0 else -1
            branch.append(item)
        if branch: branches.append(branch)
        made=False
        for bid,branch in enumerate(branches):
            win=[q for q in branch if vlo<=float(_vkey(q[1])[0])<=vhi]
            if len(win)<3: continue
            _,key=_vkey(win[0][1]); slope,rms=_fit(win,key); coverage=True; lifeval=None; lifeok=True
            if lifetime_range_ns is not None:
                lk="tau_rad_bare_ns" if lifetime_range_ns["kind"]=="bare" else "tau_rad_cavity_ns"; hits=[r for _,r in branch if float(_vkey(r)[0])==float(lifetime_range_ns["voltage_V"])]
                if not hits or not _finite(hits[0].get(lk)): coverage=False; lifeok=False
                else: lifeval=float(hits[0][lk]); lifeok=lifetime_range_ns["min_ns"]<=lifeval<=lifetime_range_ns["max_ns"]
            eps=1e-12*max(1,abs(lo),abs(hi),abs(slope)); slopeok=lo-eps<=slope<=hi+eps; status="incomplete_model_coverage" if not coverage else "compatible" if slopeok and lifeok else "incompatible"
            result.append({"group":group,"screening":hyp,"branch_id":bid,"voltage_window_V":(vlo,vhi),"fitted_slope_meV_per_V":slope,"slope_interval_meV_per_V":(slope,slope),"residual_rms_meV":rms,"nonlinearity_diagnostic":"rms_residual_meV","row_ids":[_rid(r,i) for i,r in win],"compatible":slopeok and lifeok if coverage else False,"lifetime_value_ns":lifeval,"identification_status":status,"numerical_sensitivity":"not_quantified","shape_sensitivity":"rms_residual_meV","nuisance_sensitivity":"conditional_fixed_settings"}); made=True
        if not made: result.append({"group":group,"screening":hyp,"branch_id":None,"voltage_window_V":(vlo,vhi),"fitted_slope_meV_per_V":float("nan"),"slope_interval_meV_per_V":None,"residual_rms_meV":float("nan"),"nonlinearity_diagnostic":"incomplete_model_coverage","row_ids":[],"compatible":False,"lifetime_value_ns":None,"identification_status":"incomplete_model_coverage","numerical_sensitivity":"not_quantified","shape_sensitivity":"not_available","nuisance_sensitivity":"conditional_fixed_settings"})
    bybranch={}
    for r in result: bybranch.setdefault((r["group"],r["branch_id"]),[]).append(r)
    for members in bybranch.values():
        usable=[r for r in members if _finite(r["fitted_slope_meV_per_V"])]
        if len(usable)>=2 and len({round(r["fitted_slope_meV_per_V"],12) for r in usable})==1 and len({(r["compatible"],r["identification_status"]) for r in usable})==1:
            for r in usable: r["compatible"]=None; r["identification_status"]="screening_unidentifiable"
    return result
