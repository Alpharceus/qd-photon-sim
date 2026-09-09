"""Source transcription and numerical verification for nitride_materials.

Every literal target below is transcribed independently of the module
under test (see the digest at ../_goal/nitride_digests.md); none is a
value the module produced itself.
"""
import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fsim_core.nitride_materials import binary, ingaN, bandgap, band_edges, polarization_field, EPS0_SI

checks=[]
def check(label, ok):
    checks.append(bool(ok))
    if not ok: print("FAIL", label)

g,i,a=binary("GaN"),binary("InN"),binary("AlN")

# Source transcription: Bernardini, Fiorentini & Vanderbilt, PRB 56, R10024
# (1997), Table II [V].  These literals deliberately do not come from code.
check("source transcription B97 GaN", abs(g.Psp_Cm2-(-.029))<=1e-12 and abs(g.e31_Cm2+.49)<=1e-12 and abs(g.e33_Cm2-.73)<=1e-12)
check("source transcription B97 InN", abs(i.Psp_Cm2+.032)<=1e-12 and abs(i.e31_Cm2+.57)<=1e-12 and abs(i.e33_Cm2-.97)<=1e-12)
check("source transcription B97 AlN", abs(a.Psp_Cm2+.081)<=1e-12 and abs(a.e31_Cm2+.60)<=1e-12 and abs(a.e33_Cm2-1.46)<=1e-12)
# Source transcription: Bernardini PRB 1997, p.38 [V] -- one consistent
# static dielectric-constant family, used with the SAME along-c polarization
# formula (not a mix of Ioffe's perpendicular/parallel/static values).
check("source transcription B97 eps_r", abs(g.eps_r-10.28)<=1e-12 and abs(i.eps_r-14.61)<=1e-12 and abs(a.eps_r-10.31)<=1e-12)
# Numerical verification: Rinke et al., PRB 77, 075202 (2008), Table V and
# Sec. IV.B [V] -- electron masses (all three binaries) and GaN's A-band
# hole masses.
check("source transcription Rinke electron masses", abs(g.me_xy-.186)<=1e-12 and abs(g.me_z-.209)<=1e-12 and abs(i.me_xy-.065)<=1e-12 and abs(i.me_z-.068)<=1e-12 and abs(a.me_xy-.322)<=1e-12 and abs(a.me_z-.329)<=1e-12)
check("source transcription Rinke GaN A-band hole masses", abs(g.mh_z-1.88)<=1e-12 and abs(g.mh_xy-.33)<=1e-12)
# Source transcription: Ioffe NSM archive [V] -- lattice constants and
# elastic constants (InN elastic is the Sheleg & Savastenko 1979 family).
check("source transcription Ioffe lattice constants", abs(g.a_A-3.189)<=1e-12 and abs(g.c_A-5.186)<=1e-12 and abs(a.a_A-3.112)<=1e-12 and abs(a.c_A-4.982)<=1e-12 and abs(i.a_A-3.545)<=1e-12 and abs(i.c_A-5.703)<=1e-12)
check("source transcription Ioffe/Sheleg elastic C13/C33", abs(g.C13_GPa-106.)<=1e-12 and abs(g.C33_GPa-398.)<=1e-12 and abs(a.C13_GPa-99.)<=1e-12 and abs(a.C33_GPa-389.)<=1e-12 and abs(i.C13_GPa-121.)<=1e-12 and abs(i.C33_GPa-182.)<=1e-12)
check("source transcription Ioffe AlN A-band hole masses", abs(a.mh_z-3.53)<=1e-12 and abs(a.mh_xy-10.4)<=1e-12)
# Source transcription: Rinke et al., PRB 77, 075202 (2008), Table IV [V] --
# volume band-gap deformation potentials a_V.
check("source transcription Rinke a_V", abs(g.aV_eV+7.6)<=1e-12 and abs(i.aV_eV+4.2)<=1e-12 and abs(a.aV_eV+9.8)<=1e-12)
# Source transcription: Tsai & Bayram, ACS Omega 5, 3917 (2020), Table 2 [V]
# -- valence-band offsets vs GaN.  substrate=material itself zeroes the
# coherent strain (a_sub == a_layer), isolating the raw VBO term used by
# band_edges from any deformation-potential correction.
check("source transcription Tsai VBO InN above GaN", math.isclose(band_edges(i,300,substrate=i)["Ev_eV"],1.15,abs_tol=1e-12))
check("source transcription Tsai VBO AlN below GaN", math.isclose(band_edges(a,300,substrate=a)["Ev_eV"],-0.30,abs_tol=1e-12))
check("source transcription Tsai VBO GaN reference", abs(band_edges(g,300,substrate=g)["Ev_eV"])<=1e-12)

# Numerical verification: alloy endpoints must match the underlying binary
# fields (not merely return an identical object; every field is compared).
_endpoint_fields=("a_A","c_A","C13_GPa","C33_GPa","e31_Cm2","e33_Cm2","Psp_Cm2","eps_r","me_z","me_xy","mh_z","mh_xy","Eg0_eV","alpha_eVK","beta_K","aV_eV")
check("numerical alloy endpoint x=0 matches GaN fields", all(math.isclose(getattr(ingaN(0.),k),getattr(g,k),rel_tol=1e-12) for k in _endpoint_fields))
check("numerical alloy endpoint x=1 matches InN fields", all(math.isclose(getattr(ingaN(1.),k),getattr(i,k),rel_tol=1e-12) for k in _endpoint_fields))
# Numerical verification: mid-alloy virtual-crystal interpolation reproduces
# the literal linear combination of the transcribed GaN/InN Table II
# constants above, independently recomputed here (not calling ingaN twice).
x=0.15
mid=ingaN(x)
psp_hand=(1-x)*(-.029)+x*(-.032); e31_hand=(1-x)*(-.49)+x*(-.57); e33_hand=(1-x)*.73+x*.97
check("numerical alloy VCA interpolation vs hand-computed B97 blend", math.isclose(mid.Psp_Cm2,psp_hand,rel_tol=1e-12) and math.isclose(mid.e31_Cm2,e31_hand,rel_tol=1e-12) and math.isclose(mid.e33_Cm2,e33_hand,rel_tol=1e-12))
check("numerical bowing sign", bandgap(ingaN(.5),300)<.5*(bandgap(g,300)+bandgap(i,300))) # Wu et al. APL 2002, 1.43 eV [V]
check("numerical zero mismatch", abs(band_edges(g,300,substrate=g)["eps_parallel"])<1e-15)
check("numerical polarization sign", polarization_field(ingaN(.15),g,300)<0 and polarization_field(ingaN(.15),g,300,screening_fraction=1)==0)

# Numerical verification: polarization-field MAGNITUDE at x=0.15 on GaN,
# independently re-derived here from the transcribed B97 constants and the
# documented coherent-strain formula (eps_parallel, eps_zz, P_total), rather
# than compared to itself via the module's own band_edges internals.
a_A_hand=(1-x)*g.a_A+x*i.a_A; C13_hand=(1-x)*g.C13_GPa+x*i.C13_GPa; C33_hand=(1-x)*g.C33_GPa+x*i.C33_GPa
ep_hand=(g.a_A-a_A_hand)/a_A_hand
ez_hand=-2*C13_hand/C33_hand*ep_hand
P_dot_hand=psp_hand+2*e31_hand*ep_hand+e33_hand*ez_hand
eps_r_hand=(1-x)*g.eps_r+x*i.eps_r
F_hand=(g.Psp_Cm2-P_dot_hand)/(EPS0_SI*eps_r_hand)*1e-5
check("numerical polarization field magnitude vs hand-computed B97 formula", math.isclose(polarization_field(mid,g,300), F_hand, rel_tol=1e-9))

try: ingaN(1.01); good=False
except ValueError: good=True
check("invalid composition",good)
check("model comparison GaN gap",3.39<=bandgap(g,300)<=3.44) # Guo & Yoshida class cited by Ioffe NSM [V]
print(f"{sum(checks)}/{len(checks)} nitride materials checks passed")
raise SystemExit(0 if all(checks) else 1)
