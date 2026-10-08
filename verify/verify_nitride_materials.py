"""Source transcription and numerical verification for nitride_materials.

Every literal target below is transcribed independently of the module
under test; none is a value the module produced itself.
"""
import math
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fsim_core.nitride_materials import binary, ingaN, bandgap, band_edges, polarization_field, strained_gap_shift_eV, EPS0_SI

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
# Source transcription: Bernardini & Fiorentini, phys. stat. sol. (b) 216,
# 391 (1999), Sec. III [V] -- one consistent static dielectric-constant
# family, used with the SAME along-c polarization formula.
check("source transcription Bernardini-Fiorentini 1999 eps_r", abs(g.eps_r-10.28)<=1e-12 and abs(i.eps_r-14.61)<=1e-12 and abs(a.eps_r-10.31)<=1e-12)
# Numerical verification: Rinke et al., PRB 77, 075202 (2008), Table V and
# Sec. IV.B [V] -- electron masses (all three binaries) and GaN's A-band
# hole masses.  Axis (strain-mass audit, 2026-09-23): Table V row m_e-par
# (0.186 / 0.065 / 0.322) is ALONG c (App. A Eq. (A5): 1/m_e-par multiplies
# kz^2), so it is me_z; row m_e-perp (0.209 / 0.068 / 0.329) is me_xy.
check("source transcription Rinke electron masses (m_e-par along c = me_z)", abs(g.me_z-.186)<=1e-12 and abs(g.me_xy-.209)<=1e-12 and abs(i.me_z-.065)<=1e-12 and abs(i.me_xy-.068)<=1e-12 and abs(a.me_z-.322)<=1e-12 and abs(a.me_xy-.329)<=1e-12)
# Axis cross-check independent of the electron rows: the hole masses the
# module already stores on the c axis follow m0/m_par = -(A1+A3) with the
# Rinke App. B Luttinger-like GaN parameters A1 = -5.947, A3 = 5.414 [V]
# (1/0.533 = 1.876 -> mh_z 1.88), i.e. "par" means along c for the holes too.
check("axis cross-check: Rinke -(A1+A3) GaN gives the stored c-axis hole mass", abs(1./(5.947-5.414)-g.mh_z)<0.01)
# Yan et al., PRB 90, 125118 (2014), Table III recommended rows [V]:
# acz-D1, act-D2, D3, D4 (eV) for GaN (HSE 0.30), InN (HSE 0.25), AlN (HSE 0.34).
_Y14={"GaN":(-6.07,-8.88,5.38,-2.69),"InN":(-3.64,-4.58,2.68,-1.78),"AlN":(-4.36,-12.35,9.17,-3.72)}
check("source transcription Yan 2014 Table III deformation potentials",
      all(abs(getattr(m,"acz_D1_eV")-v[0])<=1e-12 and abs(m.act_D2_eV-v[1])<=1e-12 and abs(m.D3_eV-v[2])<=1e-12 and abs(m.D4_eV-v[3])<=1e-12
          for m,v in ((g,_Y14["GaN"]),(i,_Y14["InN"]),(a,_Y14["AlN"]))))
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
_endpoint_fields=("a_A","c_A","C13_GPa","C33_GPa","e31_Cm2","e33_Cm2","Psp_Cm2","eps_r","me_z","me_xy","mh_z","mh_xy","Eg0_eV","alpha_eVK","beta_K","aV_eV","acz_D1_eV","act_D2_eV","D3_eV","D4_eV")
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

# Regression for the strained-gap construction.  Independently reconstruct
# the anisotropic A-gap strain shift from the x=0.25 VCA inputs and the Yan
# 2014 Table III literals above, linearly interpolated GaN->InN [V]/[A]; the
# conduction/valence partition [A] may move each edge but cannot change their
# separation.  The superseded volume-only Rinke a_V form (H7, 249.6 meV) is
# recomputed by hand only to document the size of the fix.
x_gap=.25; dot_gap=ingaN(x_gap); T_gap=300.
ep_gap=(g.a_A-dot_gap.a_A)/dot_gap.a_A
ez_gap=-2.*dot_gap.C13_GPa/dot_gap.C33_GPa*ep_gap
aV_hand=(1-x_gap)*(-7.6)+x_gap*(-4.2)
check("numerical x=0.25 superseded volume-only shift is 249.6 meV class (hand)",
      abs(aV_hand*(2.*ep_gap+ez_gap)*1000.-249.5950584942605)<1e-6)
_yx=[(1-x_gap)*gv+x_gap*iv for gv,iv in zip(_Y14["GaN"],_Y14["InN"])]
dEg_gap=(_yx[0]-_yx[2])*ez_gap+(_yx[1]-_yx[3])*2.*ep_gap
expected_gap=bandgap(dot_gap,T_gap)+dEg_gap
check("H7 regression: x=0.25 coherent A-gap shift is +114.0 meV (Yan 2014), not +249.6",
      abs(dEg_gap*1000.-114.0)<0.1)
partition_gaps=[]
for f_c in (0.,.2,.7,1.):
    edges=band_edges(dot_gap,T_gap,substrate=g,strain_c_fraction=f_c)
    partition_gaps.append(edges["Ec_eV"]-edges["Ev_eV"])
check("numerical strained gap equals unstrained gap plus anisotropic A-gap shift",
      all(abs(q-expected_gap)<1e-6 for q in partition_gaps))
check("numerical strained gap is partition-independent",
      max(partition_gaps)-min(partition_gaps)<1e-12)

# H7 closed form (strain-mass audit, 2026-09-23): pure c-plane biaxial strain
# on GaN, dEg_A = (acz-D1-D3) ezz + (act-D2-D4) 2 exx with ezz = -2 C13/C33
# exx [Yan et al., PRB 90, 125118 (2014), Eqs. (3), (9); Chuang & Chang,
# PRB 54, 2491 (1996)], from the literals -6.07, -8.88, 5.38, -2.69 eV and
# C13/C33 = 106/398 GPa, never from module fields.
for exx in (-0.01,0.005):
    ezz=-2.*106./398.*exx
    closed=(-6.07-5.38)*ezz+(-8.88+2.69)*2.*exx
    check(f"H7 biaxial GaN A-gap shift equals Yan 2014 closed form to 1e-9 (exx={exx})",
          abs(strained_gap_shift_eV(g,exx,ezz)-closed)<1e-9)
# Same closed form through band_edges: GaN coherent on AlN (exx = (a_AlN-a_GaN)/a_GaN).
exx_al=(3.112-3.189)/3.189; ezz_al=-2.*106./398.*exx_al
e_al=band_edges(g,300.,substrate=a)
check("H7 band_edges GaN-on-AlN gap equals Eg + Yan 2014 closed form to 1e-9",
      abs((e_al["Ec_eV"]-e_al["Ev_eV"])-(bandgap(g,300.)+(-6.07-5.38)*ezz_al+(-8.88+2.69)*2.*exx_al))<1e-9)
# Hydrostatic strain exx=eyy=ezz=e reduces to a volume potential:
# dEg = 3 a_V,eff e with a_V,eff = [(acz-D1-D3) + 2(act-D2-D4)]/3 = -7.943 eV
# for GaN [DR from the literals], which must agree with the independent
# Rinke 2008 Table IV a_V = -7.6 eV [V] to within 5 % (different methods:
# HSE vs G0W0; Yan 2014 itself compares them).
e_h=1e-3
aV_eff=strained_gap_shift_eV(g,e_h,e_h)/(3.*e_h)
check("hydrostatic strain reduces to the volume form, a_V,eff = [(acz-D1-D3)+2(act-D2-D4)]/3 = -7.943 eV (GaN) to 1e-9",
      abs(aV_eff-((-6.07-5.38)+2.*(-8.88+2.69))/3.)<1e-9 and abs(aV_eff+7.9433333333)<1e-9)
check("hydrostatic a_V,eff agrees with Rinke 2008 a_V = -7.6 eV within 5 %", abs(aV_eff/(-7.6)-1.)<0.05)
# Electron DOS mass of bulk GaN with the corrected axes: (0.209^2 * 0.186)^(1/3).
check("GaN electron DOS mass (me_xy^2 me_z)^(1/3) = 0.2010 with the Rinke axes",
      abs((g.me_xy**2*g.me_z)**(1./3.)-(0.209**2*0.186)**(1./3.))<1e-12 and abs((0.209**2*0.186)**(1./3.)-0.2010)<1e-4)

try: ingaN(1.01); good=False
except ValueError: good=True
check("invalid composition",good)
check("model comparison GaN gap",3.39<=bandgap(g,300)<=3.44) # Guo & Yoshida class cited by Ioffe NSM [V]
print(f"{sum(checks)}/{len(checks)} nitride materials checks passed")
raise SystemExit(0 if all(checks) else 1)
