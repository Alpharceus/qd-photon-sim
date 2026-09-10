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
from fsim_core.device import DeviceDesign, evaluate
from fsim_core import nitride_transport
import yaml
checks=[]
def check(name, value):
    checks.append(bool(value)); print(("ok   " if value else "FAIL ")+name)
def close(a,b): return math.isclose(a,b,rel_tol=1e-10,abs_tol=1e-10)
def raises(fn):
    try: fn()
    except ValueError: return True
    return False
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
print(f"{sum(checks)}/{len(checks)} nitride geometry device checks passed")
raise SystemExit(0 if all(checks) else 1)
