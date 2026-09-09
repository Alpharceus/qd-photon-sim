"""Wurtzite InGaN/GaN material inputs, kept separate from cubic materials.

The polarization convention is positive along +c.  Alloy interpolation is a
virtual-crystal approximation [A]; it is not a claim that an alloy value was
measured.  Thermal expansion is omitted [A].
"""
from dataclasses import dataclass
import math

EPS0_SI = 8.8541878128e-12
KB_EV = 8.617333262145e-5

@dataclass(frozen=True)
class NitrideMaterial:
    name: str; x_in: float; a_A: float; c_A: float; C13_GPa: float; C33_GPa: float
    e31_Cm2: float; e33_Cm2: float; Psp_Cm2: float; eps_r: float
    me_z: float; me_xy: float; mh_z: float; mh_xy: float; provenance: str
    Eg0_eV: float; alpha_eVK: float; beta_K: float; aV_eV: float

# Bernardini, Fiorentini & Vanderbilt, PRB 56, R10024 (1997), Table II [V]
# (Psp, e31, e33) and p.38 [V] (eps_r: a single internally-consistent static
# dielectric-constant family used with the SAME polarization formula below --
# NOT Ioffe's mixed perpendicular/parallel/static set, which mismatches the
# along-c polarization screening by about 10%).
# Rinke et al., PRB 77, 075202 (2008), Table V and Sec. IV.B [V] (masses).
# Lattice and elastic constants: Ioffe NSM summary [V]; InN elastic choice is
# the Sheleg family [E].  Varshni pairs: Vurgaftman & Meyer, JAP 94, 3675
# (2003) framework [E] pending primary-text confirmation.
_B = {
 "GaN": NitrideMaterial("GaN",0.,3.189,5.186,106.,398.,-.49,.73,-.029,10.28,.209,.186,1.88,.33,"[V] Bernardini PRB 1997 Table II + p.38 eps_r; Rinke PRB 2008 Table V/Sec IV.B masses",3.51,9.09e-4,830.,-7.6),
 "InN": NitrideMaterial("InN",1.,3.545,5.703,121.,182.,-.57,.97,-.032,14.61,.068,.065,1.63,1.63,"[V/E] Bernardini PRB 1997 Table II + p.38 eps_r; Rinke PRB 2008 electron masses [V]; Ioffe NSM hole mass 1.63 applied isotropically [E]",.78,2.45e-4,624.,-4.2),
 "AlN": NitrideMaterial("AlN",0.,3.112,4.982,99.,389.,-.60,1.46,-.081,10.31,.329,.322,3.53,10.4,"[V] Bernardini PRB 1997 Table II + p.38 eps_r; Rinke PRB 2008 electron masses; Ioffe NSM hole masses 3.53/10.4",6.25,1.799e-3,1462.,-9.8),
}

# Tsai & Bayram, ACS Omega 5, 3917 (2020), Table 2 [V]: AlN valence 0.30 eV
# below GaN.  AlN is not a point on the In_xGa_1-xN composition axis (its
# x_in is a placeholder 0.0 that only means "not an InGaN alloy"), so its
# valence-band offset cannot be read off x_in*vbo_InN_GaN_eV the way the
# InGaN alloy family's can -- doing so silently aliased AlN with GaN
# (x_in=0.0 for both) and put the AlN valence edge above GaN instead of
# below it.
_VBO_ALN_EV = -0.30

def binary(name):
    """Return an immutable binary; accepted names are GaN, InN, and AlN."""
    if name not in _B: raise ValueError("nitride binary must be GaN, InN, or AlN")
    return _B[name]

def ingaN(x_in):
    """In_xGa_(1-x)N virtual crystal; 1.4 eV bowing is Wu et al. [V]."""
    if not math.isfinite(x_in) or not 0. <= x_in <= 1.: raise ValueError("x_in must be finite and in [0, 1]")
    if x_in == 0.: return binary("GaN")
    if x_in == 1.: return binary("InN")
    g,i=_B["GaN"],_B["InN"]; q=lambda k:(1-x_in)*getattr(g,k)+x_in*getattr(i,k)
    return NitrideMaterial("InGaN",x_in,*[q(k) for k in ("a_A","c_A","C13_GPa","C33_GPa","e31_Cm2","e33_Cm2","Psp_Cm2","eps_r","me_z","me_xy","mh_z","mh_xy")],"[A] linear virtual-crystal alloy; Wu et al. bowing [V]",q("Eg0_eV"),q("alpha_eVK"),q("beta_K"),q("aV_eV"))

def bandgap(material, T_K):
    if not math.isfinite(T_K) or T_K <= 0: raise ValueError("T_K must be positive and finite")
    return material.Eg0_eV-material.alpha_eVK*T_K*T_K/(T_K+material.beta_K)-1.4*material.x_in*(1-material.x_in) # Wu et al., APL 2002, reported by OSTI 835986 [V]

def band_edges(material, T_K, *, substrate=None, strain_fraction=1., vbo_InN_GaN_eV=1.15, strain_c_fraction=.7):
    if not 0 <= strain_fraction <= 1 or not 0 <= strain_c_fraction <= 1: raise ValueError("strain fractions must be in [0, 1]")
    sub = binary("GaN") if substrate is None else substrate
    ep=(sub.a_A-material.a_A)/material.a_A*strain_fraction
    ez=-2*material.C13_GPa/material.C33_GPa*ep
    p=material.Psp_Cm2+2*material.e31_Cm2*ep+material.e33_Cm2*ez
    dEg=material.aV_eV*(2*ep+ez) # Rinke et al., PRB 77, 075202 (2008), Table IV [V]
    # Valence reference: the InN/GaN family uses the linear x_in transfer of
    # the 1.15 eV Tsai & Bayram VBO [V]/[A]; AlN is off that axis and uses
    # its own Tsai & Bayram Table 2 literal [V] instead (see _VBO_ALN_EV).
    ev0 = _VBO_ALN_EV if material.name == "AlN" else material.x_in*vbo_InN_GaN_eV
    ev=ev0-(1-strain_c_fraction)*dEg # VBO Tsai & Bayram ACS Omega 2020 Table 2 [V]; partition [A]
    return dict(Ec_eV=ev+bandgap(material,T_K)+strain_c_fraction*dEg,Ev_eV=ev,eps_parallel=ep,eps_zz=ez,P_total_Cm2=p,provenance="[V] B97 polarization; [V] Rinke 2008 volume deformation; [V] Tsai 2020 VBO endpoints; [A] 0.7 conduction partition, linear alloy VBO transfer, no thermal expansion")

def polarization_field(dot,matrix,T_K,*,strain_fraction=1.,screening_fraction=0.,external_field_kVcm=0.):
    if not 0 <= screening_fraction <= 1: raise ValueError("screening_fraction must be in [0, 1]")
    pd=band_edges(dot,T_K,substrate=matrix,strain_fraction=strain_fraction)["P_total_Cm2"]
    pm=band_edges(matrix,T_K,substrate=matrix,strain_fraction=0.)["P_total_Cm2"]
    return (1-screening_fraction)*(pm-pd)/(EPS0_SI*dot.eps_r)*1e-5+external_field_kVcm # fixed-D thick GaN reservoirs [A]
