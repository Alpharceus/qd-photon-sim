"""Source transcription and numerical verification for nitride_materials."""
import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fsim_core.nitride_materials import binary, ingaN, bandgap, band_edges, polarization_field

checks=[]
def check(label, ok):
    checks.append(bool(ok))
    if not ok: print("FAIL", label)

# Source transcription: Bernardini, Fiorentini & Vanderbilt, PRB 56, R10024
# (1997), Table II [V].  These literals deliberately do not come from code.
g,i,a=binary("GaN"),binary("InN"),binary("AlN")
check("source transcription B97 GaN", abs(g.Psp_Cm2-(-.029))<=1e-12 and abs(g.e31_Cm2+.49)<=1e-12 and abs(g.e33_Cm2-.73)<=1e-12)
check("source transcription B97 InN", abs(i.Psp_Cm2+.032)<=1e-12 and abs(i.e31_Cm2+.57)<=1e-12 and abs(i.e33_Cm2-.97)<=1e-12)
check("source transcription B97 AlN", abs(a.Psp_Cm2+.081)<=1e-12 and abs(a.e31_Cm2+.60)<=1e-12 and abs(a.e33_Cm2-1.46)<=1e-12)
# Numerical verification: Rinke et al., PRB 77, 075202 (2008), Table V [V].
check("source transcription Rinke masses", abs(g.me_xy-.186)<=1e-12 and abs(g.me_z-.209)<=1e-12 and abs(i.me_xy-.065)<=1e-12 and abs(a.me_z-.329)<=1e-12)
check("numerical alloy endpoints", math.isclose(ingaN(0).a_A,g.a_A,rel_tol=1e-12) and math.isclose(ingaN(1).c_A,i.c_A,rel_tol=1e-12))
check("numerical bowing sign", bandgap(ingaN(.5),300)<.5*(bandgap(g,300)+bandgap(i,300))) # Wu et al. APL 2002, 1.43 eV [V]
check("numerical zero mismatch", abs(band_edges(g,300,substrate=g)["eps_parallel"])<1e-15)
check("numerical polarization sign", polarization_field(ingaN(.15),g,300)<0 and polarization_field(ingaN(.15),g,300,screening_fraction=1)==0)
try: ingaN(1.01); good=False
except ValueError: good=True
check("invalid composition",good)
check("model comparison GaN gap",3.39<=bandgap(g,300)<=3.44) # Guo & Yoshida class cited by Ioffe NSM [V]
print(f"{sum(checks)}/{len(checks)} nitride materials checks passed")
raise SystemExit(0 if all(checks) else 1)
