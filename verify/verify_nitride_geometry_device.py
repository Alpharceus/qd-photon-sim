"""Independent integration checks for geometry-aware nitride device rows.

Electrical checks use the abrupt p-i-n equations [DR Sze & Ng, Physics of
Semiconductor Devices, 3rd ed. (2007)].
"""
from __future__ import annotations
import copy, math, sys
from dataclasses import asdict
from pathlib import Path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
import fsim_core.device as device_mod
from fsim_core.device import DeviceDesign, evaluate
from fsim_core import nitride_transport, nitride_levels
import yaml
checks=[]
def check(name, value):
    checks.append(bool(value)); print(("ok   " if value else "FAIL ")+name)
def close(a,b): return math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-10)
def raises(fn):
    try: fn()
    except ValueError: return True
    return False
def nan_eq(a, b):
    """Recursive NaN-aware equality for evaluate() scalar dicts: floats,
    bools, strings, nested lists/tuples (invalid_reasons) and dicts
    (provenance) all compare equal iff every leaf matches, with NaN==NaN."""
    if isinstance(a, bool) or isinstance(b, bool):
        return a is b
    if isinstance(a, float) and isinstance(b, float):
        if math.isnan(a) and math.isnan(b): return True
        return a == b
    if isinstance(a, dict) and isinstance(b, dict):
        return set(a) == set(b) and all(nan_eq(a[k], b[k]) for k in a)
    if isinstance(a, (list, tuple)) and isinstance(b, (list, tuple)):
        return len(a) == len(b) and all(nan_eq(x, y) for x, y in zip(a, b))
    return a == b

base=DeviceDesign.load(str(ROOT/"cards"/"nitride-cavity-pulse-design.yaml"))
s0=evaluate(base,[300.])["scalars"]
check("default source temperature mode",s0["evaluation_kind"]=="source" and s0["temperature_mode"]=="self_consistent")
check("default applied diode field is finite",s0["field_polarity"]==1 and math.isfinite(s0["applied_field_kVcm"]))
check("bare reciprocal radiative lifetime",close(s0["tau_rad_bare_ns"],1/s0["gamma_X0_ns"]))
check("cavity reciprocal radiative lifetime",close(s0["tau_rad_cavity_ns"],1/s0["gamma_X_ns"]))
dv=copy.deepcopy(base)
dv.nitride["bias"]={"mode":"junction_voltage","V_j_V":.5,"T_j_K":300.,"field_polarity":1,"cavity_reference_V_j_V":.5}
sv=evaluate(dv,[300.])["scalars"]
diode=nitride_transport.planar_pin(**{k:v for k,v in dv.drive.diode.items() if k not in ("preset","tau_pulse_ns")})
expected_I=diode.j_of_vj(.5,300.)*diode.area_cm2*1e6
check("fixed requested temperature",sv["T_j"]==300. and sv["temperature_mode"]=="fixed_junction")
check("resolved voltage current",sv["V_j"]==.5 and close(sv["resolved_current_uA"],expected_I))
check("series terminal voltage",close(sv["V_terminal"],.5+expected_I*1e-6*diode.R_s_ohm))
check("independent depletion field",close(sv["diode_field_kVcm"],diode.depletion(.5,300.).F_kVcm))
check("diagnostic excluded from device pass",sv["evaluation_kind"]=="stark_diagnostic" and sv["device_pass"] is False)
check("reference offset at anchor",abs(sv["detuning_meV"]-dv.nitride["cavity"]["detuning_offset_meV"])<1e-8)
dp=copy.deepcopy(dv); dp.nitride["bias"]["V_j_V"]=.8
dm=copy.deepcopy(dv); dm.nitride["bias"]["field_polarity"]=-1
sp=evaluate(dp,[300.])["scalars"]; sm=evaluate(dm,[300.])["scalars"]
check("frozen reference detects voltage Stark detuning",sp["detuning_meV"]!=sv["detuning_meV"])
check("polarity changes field",sm["applied_field_kVcm"]!=sv["applied_field_kVcm"])
npur=copy.deepcopy(dv); npur.nitride["cavity"]["purcell_enabled"]=False
snp=evaluate(npur,[300.])["scalars"]
check("purcell disabled lifetime equality",close(snp["tau_rad_bare_ns"],snp["tau_rad_cavity_ns"]))
check("escape is cavity independent",close(snp["k_X_ns"],sv["k_X_ns"]))
check("contradictory diagnostic temperature rejected",raises(lambda:evaluate(dv,[299.])))
def bad_current():
    q=copy.deepcopy(base); q.nitride["bias"]={"mode":"current","V_j_V":.5}; evaluate(q,[300.])
def bad_voltage():
    q=copy.deepcopy(base); q.nitride["bias"]={"mode":"junction_voltage","V_j_V":-1.,"T_j_K":300.}; evaluate(q,[300.])
def qw_mismatch():
    q=copy.deepcopy(base); q.nitride["dot"].update({"geometry_type":"qw_fluctuation","wl_thickness_nm":1.}); evaluate(q,[300.])
check("current mode contradictory voltage rejected",raises(bad_current))
check("negative junction voltage rejected",raises(bad_voltage))
check("QW diode and well mismatch rejected",raises(qw_mismatch))
text=yaml.safe_dump({"design":asdict(dv)},sort_keys=False)
old_read=Path.read_text
Path.read_text=lambda self,*a,**kw:text
try:
    rt=DeviceDesign.load("<in-memory diagnostic round trip>")
finally:
    Path.read_text=old_read
srt=evaluate(rt,[300.])["scalars"]
check("bias mapping survives in-memory load round trip",rt.nitride["bias"]==dv.nitride["bias"] and close(srt["V_j"],sv["V_j"]))

# ---------------------------------------------------------------------------
# Fix-round-1 additions below. Each block is an INDEPENDENT integration check
# (acceptance criteria 1, 2, 5, 6, 7 of .workers/specs/nitride-geometry-device.md's
# fix round) plus the three named mutation-sensitivity checks the Opus
# re-review's Required list called for by name.
# ---------------------------------------------------------------------------

# --- Acceptance 1: NaN-aware old-card equality, absent vs explicit defaults ---
d_dot_explicit=copy.deepcopy(base)
d_dot_explicit.nitride["dot"].update(orientation="c_plane",polarization_factor=None,shape="disc",
    top_radius_fraction=1.0,geometry_type="isolated_dot",shape_height_fraction=None)
s_dot_explicit=evaluate(d_dot_explicit,[300.])["scalars"]
check("old card + explicit default geometry settings reproduce absent-settings scalars exactly (NaN-aware)",
      nan_eq(s0,s_dot_explicit))
d_bias_explicit=copy.deepcopy(base)
d_bias_explicit.nitride["bias"]={"mode":"current","field_polarity":1}
s_bias_explicit=evaluate(d_bias_explicit,[300.])["scalars"]
check("old card + explicit default bias settings reproduce absent-bias scalars exactly (NaN-aware)",
      nan_eq(s0,s_bias_explicit))

# --- Acceptance 2: polarity direction, external field, nonpolar row, screening ---
d_field=copy.deepcopy(base); d_field.nitride["dot"]["external_field_kVcm"]=50.0
d_field_plus=copy.deepcopy(d_field); d_field_plus.nitride["bias"]={"mode":"current","field_polarity":1}
d_field_minus=copy.deepcopy(d_field); d_field_minus.nitride["bias"]={"mode":"current","field_polarity":-1}
s_field_plus=evaluate(d_field_plus,[300.])["scalars"]
s_field_minus=evaluate(d_field_minus,[300.])["scalars"]
check("static external field reaches the operating-point solve",s_field_plus["E_X_eV"]!=s0["E_X_eV"])
check("static external field reaches the cavity reference solve",
      s_field_plus["cavity_reference_transition_eV"]!=s0["cavity_reference_transition_eV"])
check("polarity +1 vs -1 differ in the documented sign: applied field gap is exactly 2*diode_field_kVcm "
      "(same diode field magnitude at the same current-mode operating point, only its added sign flips)",
      close(s_field_plus["diode_field_kVcm"],s_field_minus["diode_field_kVcm"]) and
      close(s_field_plus["applied_field_kVcm"]-s_field_minus["applied_field_kVcm"],2*s_field_plus["diode_field_kVcm"]))
d_np0=copy.deepcopy(base); d_np0.nitride["dot"].update(orientation="m_plane",screening_fraction=0.0)
d_np1=copy.deepcopy(base); d_np1.nitride["dot"].update(orientation="m_plane",screening_fraction=1.0)
s_np0=evaluate(d_np0,[300.])["scalars"]; s_np1=evaluate(d_np1,[300.])["scalars"]
check("nonpolar (m_plane) row: field stays nonzero below flat band (bias-driven depletion field only)",
      s_np0["field_kVcm"]!=0.0 and math.isfinite(s_np0["field_kVcm"]))
check("nonpolar (m_plane) row: E_X is screening-invariant to 1e-9 eV (zero intrinsic polarization field)",
      abs(s_np0["E_X_eV"]-s_np1["E_X_eV"])<1e-9)
d_screen=copy.deepcopy(base); d_screen.nitride["dot"]["screening_fraction"]=0.5
s_screen=evaluate(d_screen,[300.])["scalars"]
check("c-plane row: screening changes E_X at a valid reference",s_screen["E_X_eV"]!=s0["E_X_eV"])

# --- Required: polarity-direction / external-field / reservoir mutation sensitivity ---
def _with_patched_resolve(patched_fn, fn):
    real=device_mod.resolve_bias; device_mod.resolve_bias=patched_fn
    try: return fn()
    finally: device_mod.resolve_bias=real
def _negate_polarity(diode,*,T_j_K,current_uA=None,junction_voltage_V=None,field_polarity=1,external_field_kVcm=0.0):
    from fsim_core.nitride_stark import resolve_bias as _real
    return _real(diode,T_j_K=T_j_K,current_uA=current_uA,junction_voltage_V=junction_voltage_V,
                 field_polarity=-int(field_polarity),external_field_kVcm=external_field_kVcm)
s_patched_polarity=_with_patched_resolve(_negate_polarity, lambda: evaluate(d_field_plus,[300.])["scalars"])
check("mutation caught: negating field_polarity inside resolve_bias changes the operating-point field",
      s_patched_polarity["E_X_eV"]!=s_field_plus["E_X_eV"] and close(s_patched_polarity["E_X_eV"],s_field_minus["E_X_eV"]))
def _zero_external(diode,*,T_j_K,current_uA=None,junction_voltage_V=None,field_polarity=1,external_field_kVcm=0.0):
    from fsim_core.nitride_stark import resolve_bias as _real
    return _real(diode,T_j_K=T_j_K,current_uA=current_uA,junction_voltage_V=junction_voltage_V,
                 field_polarity=field_polarity,external_field_kVcm=0.0)
s_patched_extfield=_with_patched_resolve(_zero_external, lambda: evaluate(d_field_plus,[300.])["scalars"])
check("mutation caught: forcing external_field_kVcm=0 inside resolve_bias changes the operating-point field",
      s_patched_extfield["E_X_eV"]!=s_field_plus["E_X_eV"])
check("mutation caught: forcing external_field_kVcm=0 inside resolve_bias also changes the cavity reference solve",
      s_patched_extfield["cavity_reference_transition_eV"]!=s_field_plus["cavity_reference_transition_eV"])

# --- Acceptance 5 / Required: QW reservoir sourced from the levels object, not a bulk edge ---
def _qw_card(w,height=3.5,radius=20.,screening=1.0):
    q=copy.deepcopy(base)
    q.nitride["dot"].update(geometry_type="qw_fluctuation",height_nm=height,radius_nm=radius,
                             wl_thickness_nm=w,screening_fraction=screening)
    q.drive.diode["wl_thickness_nm"]=w
    return q
qw_rows={w:evaluate(_qw_card(w),[300.])["scalars"] for w in (1.0,1.5,2.0,2.5)}
check("QW row is valid and reports the ingan_qw reservoir kind",
      all(r["valid"] and r["reservoir_kind"]=="ingan_qw" for r in qw_rows.values()))
# Independent recomputation: call nitride_levels.levels() directly on the SAME
# bias-resolved system the device row used, and compare its reservoir energy --
# not merely internal self-consistency, a genuine second computation path.
_qw_diode=nitride_transport.planar_pin(**{k:v for k,v in _qw_card(1.0).drive.diode.items() if k not in ("preset","tau_pulse_ns")})
from fsim_core.nitride_stark import resolve_bias as _resolve_bias_direct
_qw_dotkw=dict(_qw_card(1.0).nitride["dot"])
_qw_bias=_resolve_bias_direct(_qw_diode,T_j_K=qw_rows[1.0]["T_j"],current_uA=base.drive.I_uA,field_polarity=1,external_field_kVcm=0.0)
_qw_lv=nitride_levels.levels(nitride_levels.NitrideDotSystem(**{**_qw_dotkw,"external_field_kVcm":_qw_bias["applied_field_kVcm"]}),qw_rows[1.0]["T_j"])
check("QW reservoir_energy_eV agrees with an independent nitride_levels.levels() call",
      close(qw_rows[1.0]["reservoir_energy_eV"],_qw_lv.reservoir_energy_eV))
check("QW reservoir moves with w (surrounding-well thickness), not w-independent",
      len({round(r["reservoir_energy_eV"],6) for r in qw_rows.values()})==4 and
      qw_rows[1.0]["reservoir_energy_eV"]>qw_rows[1.5]["reservoir_energy_eV"]>qw_rows[2.0]["reservoir_energy_eV"]>qw_rows[2.5]["reservoir_energy_eV"])
check("QW transport/background observable (reservoir_offset_meV) responds to w",
      len({round(r["reservoir_offset_meV"],3) for r in qw_rows.values()})==4)
qw_h1=evaluate(_qw_card(1.0,height=3.5),[300.])["scalars"]; qw_h2=evaluate(_qw_card(1.0,height=5.0),[300.])["scalars"]
check("fixed w: changing fluctuation height moves the local dot energy E_X",qw_h1["E_X_eV"]!=qw_h2["E_X_eV"])
check("fixed w: changing fluctuation height leaves the surrounding-QW reservoir unmoved",
      close(qw_h1["reservoir_energy_eV"],qw_h2["reservoir_energy_eV"]))
_expected_iso=device_mod._nitride_reservoir_energy_eV(base.nitride["dot"],s0["T_j"],base.nitride.get("background",{}))
check("isolated-dot card reservoir behavior is unchanged (still the bulk-edge helper)",
      close(s0["reservoir_energy_eV"],_expected_iso) and s0["reservoir_kind"]=="gan_barrier")

# --- Required: constant-reservoir mutation caught for isolated dots, QW stays decoupled ---
_real_reservoir_helper=device_mod._nitride_reservoir_energy_eV
device_mod._nitride_reservoir_energy_eV=lambda *a,**kw: 3.0
try:
    s_iso_patched=evaluate(base,[300.])["scalars"]
    s_qw_patched=evaluate(_qw_card(1.0),[300.])["scalars"]
finally:
    device_mod._nitride_reservoir_energy_eV=_real_reservoir_helper
check("mutation caught: replacing the bulk reservoir helper with a constant changes an isolated-dot row",
      s_iso_patched["reservoir_energy_eV"]==3.0 and s0["reservoir_energy_eV"]!=3.0)
check("QW row stays decoupled from the (mutated) bulk helper: reservoir_energy_eV is unaffected",
      close(s_qw_patched["reservoir_energy_eV"],qw_rows[1.0]["reservoir_energy_eV"]))

# --- Acceptance 3 / Required: V_j=0 point and a near-flat-band point, marker present ---
_diode0=nitride_transport.planar_pin(**{k:v for k,v in base.drive.diode.items() if k not in ("preset","tau_pulse_ns")})
def _find_flat_band_v(diode,T=300.,hi=6.0,steps=6000):
    for i in range(steps+1):
        vj=hi*i/steps
        if diode.depletion(vj,T).flat_band: return vj
    raise RuntimeError("flat_band not reached in scan range")
_v_flat=_find_flat_band_v(_diode0)
def _diag(vj,Tj=300.):
    q=copy.deepcopy(base); q.nitride["bias"]={"mode":"junction_voltage","V_j_V":vj,"T_j_K":Tj}
    return evaluate(q,[Tj])["scalars"]
r_vj0=_diag(0.0)
check("V_j=0 point: flat_band/depletion_regime marker present and not flat-band",
      "flat_band" in r_vj0 and r_vj0["flat_band"] is False and r_vj0["depletion_regime"]=="depleted")
check("V_j=0 point: bound-state spectrum is spectroscopically valid at zero injection",
      r_vj0["spectroscopy_valid"] is True and math.isfinite(r_vj0["E_X_eV"]))
check("V_j=0 point: zero current makes the source/counting row invalid, diagnostic device_pass False",
      r_vj0["valid"] is False and r_vj0["device_pass"] is False)
r_near_flat=_diag(_v_flat-0.01)
r_flat=_diag(_v_flat)
check("near-flat-band point: depletion_regime marker reads depleted just below V_bi",
      r_near_flat["depletion_regime"]=="depleted" and r_near_flat["flat_band"] is False)
check("at/above V_bi: depletion_regime marker flips to flat_band",
      r_flat["depletion_regime"]=="flat_band" and r_flat["flat_band"] is True)

# --- Acceptance 6: current-mode T_j responds to thermal resistance/current; fixed_junction clamped ---
d_i_lo=copy.deepcopy(base); d_i_lo.drive.I_uA=0.02
d_i_hi=copy.deepcopy(base); d_i_hi.drive.I_uA=0.2
r_i_lo=evaluate(d_i_lo,[300.])["scalars"]; r_i_hi=evaluate(d_i_hi,[300.])["scalars"]
check("current mode: T_j increases with drive current (self-consistent self-heating)",
      r_i_hi["T_j"]>r_i_lo["T_j"])
d_therm_lo=copy.deepcopy(base); d_therm_lo.thermal.mesa_diameter_um=20.0
d_therm_hi=copy.deepcopy(base); d_therm_hi.thermal.mesa_diameter_um=40.0
r_therm_lo=evaluate(d_therm_lo,[300.])["scalars"]; r_therm_hi=evaluate(d_therm_hi,[300.])["scalars"]
check("current mode: T_j responds to thermal resistance (mesa diameter)",
      r_therm_lo["T_j"]!=r_therm_hi["T_j"])
dj_lo=copy.deepcopy(base); dj_lo.nitride["bias"]={"mode":"junction_voltage","V_j_V":2.0,"T_j_K":350.}
dj_hi=copy.deepcopy(dj_lo); dj_hi.thermal.mesa_diameter_um=40.0
rj_lo=evaluate(dj_lo,[350.])["scalars"]; rj_hi=evaluate(dj_hi,[350.])["scalars"]
check("fixed_junction mode: T_j stays clamped to the requested value regardless of thermal resistance",
      rj_lo["T_j"]==350.0 and rj_hi["T_j"]==350.0)
check("fixed_junction mode: T_hs reports the card's own heat-sink setting, not overwritten by T_j",
      rj_lo["T_hs"]==base.thermal.T_hs and rj_hi["T_hs"]==base.thermal.T_hs)

# --- Required: NaN (not inf) lifetimes and no locals() leak on a non-converged row ---
d_runaway=copy.deepcopy(base); d_runaway.thermal.T_hs=3000.
r_runaway=evaluate(d_runaway,[300.,3000.])["scalars"]
check("non-converged row: valid is False and device_pass is False",
      r_runaway["valid"] is False and r_runaway["device_pass"] is False)
check("non-converged row: tau_rad_bare_ns/tau_rad_cavity_ns are NaN, not inf (rate not computed, not zero)",
      math.isnan(r_runaway["tau_rad_bare_ns"]) and math.isnan(r_runaway["tau_rad_cavity_ns"]))
check("non-converged row: reservoir_energy_eV/reservoir_offset_meV are NaN, not a leaked prior row's value",
      math.isnan(r_runaway["reservoir_energy_eV"]) and math.isnan(r_runaway["reservoir_offset_meV"]))
check("mutation-sensitivity regression: a leaked-locals reservoir value would equal the 300K row's ~3.41 eV; it does not",
      not close(r_runaway["reservoir_energy_eV"] if math.isfinite(r_runaway["reservoir_energy_eV"]) else -1.0, s0["reservoir_energy_eV"]))
check("restored pre-9d7ebd3 semantics: pair_supply_possible reflects the REQUESTED current, "
      "independent of the thermal failure that invalidated this row",
      r_runaway["pair_supply_possible"] is True)

# --- Required: polarization_factor / shape_height_fraction under the ValueError card-schema contract ---
def _bad_polarization_str():
    q=copy.deepcopy(base); q.nitride["dot"]["polarization_factor"]="0.2e0"; evaluate(q,[300.])
def _bad_polarization_bool():
    q=copy.deepcopy(base); q.nitride["dot"]["polarization_factor"]=True; evaluate(q,[300.])
def _bad_shape_height_nan():
    q=copy.deepcopy(base); q.nitride["dot"]["shape_height_fraction"]=float("nan"); evaluate(q,[300.])
check("malformed polarization_factor (non-numeric string) raises ValueError",raises(_bad_polarization_str))
check("malformed polarization_factor (bool) raises ValueError",raises(_bad_polarization_bool))
check("malformed shape_height_fraction (NaN) raises ValueError",raises(_bad_shape_height_nan))

print(f"{sum(checks)}/{len(checks)} nitride geometry device checks passed")
raise SystemExit(0 if all(checks) else 1)
