"""Wurtzite InGaN/GaN material inputs, kept separate from cubic materials.

The polarization convention is positive along +c; returned polarization
fields are the c-axis projection in this sign convention.  Orientation
factors are non-negative and do not change that polarity, so consumers add
the depletion field directly.  Alloy interpolation is a
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
# (Psp, e31, e33).  Static eps_r = 10.28/14.61/10.31 for GaN/InN/AlN is
# Bernardini & Fiorentini, phys. stat. sol. (b) 216, 391 (1999), Sec. III [V]:
# one internally-consistent dielectric family used with the SAME polarization
# formula below, not Ioffe's mixed perpendicular/parallel/static set.
# Rinke et al., PRB 77, 075202 (2008), Table V and Sec. IV.B [V] (masses).
# Lattice and elastic constants: Ioffe NSM summary [V]; InN elastic choice is
# the Sheleg family [E].  Varshni pairs: Vurgaftman & Meyer, JAP 94, 3675
# (2003) framework [E] pending primary-text confirmation.
_B = {
 "GaN": NitrideMaterial("GaN",0.,3.189,5.186,106.,398.,-.49,.73,-.029,10.28,.209,.186,1.88,.33,"[V] Bernardini PRB 1997 Table II polarization; Bernardini & Fiorentini pss(b) 1999 Sec. III eps_r; Rinke PRB 2008 Table V/Sec IV.B masses",3.51,9.09e-4,830.,-7.6),
 "InN": NitrideMaterial("InN",1.,3.545,5.703,121.,182.,-.57,.97,-.032,14.61,.068,.065,1.63,1.63,"[V] Bernardini PRB 1997 Table II polarization; Bernardini & Fiorentini pss(b) 1999 Sec. III eps_r; Rinke PRB 2008 electron masses [V]; Ioffe NSM hole mass 1.63 applied isotropically [E]",.78,2.45e-4,624.,-4.2),
 "AlN": NitrideMaterial("AlN",0.,3.112,4.982,99.,389.,-.60,1.46,-.081,10.31,.329,.322,3.53,10.4,"[V] Bernardini PRB 1997 Table II polarization; Bernardini & Fiorentini pss(b) 1999 Sec. III eps_r; Rinke PRB 2008 electron masses; Ioffe NSM hole masses 3.53/10.4",6.25,1.799e-3,1462.,-9.8),
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
    # Both strained edges are constructed from the same unshifted reference.
    # Therefore Ec-Ev = Eg+dEg for every partition, as required by the volume
    # deformation potential; using `ev` here would apply the valence share to
    # Ec a second time. [DR] Rinke et al., PRB 77, 075202 (2008), Table IV.
    ec=ev0+bandgap(material,T_K)+strain_c_fraction*dEg
    return dict(Ec_eV=ec,Ev_eV=ev,eps_parallel=ep,eps_zz=ez,P_total_Cm2=p,provenance="[V] B97 polarization; [V] Rinke 2008 volume deformation; [V] Tsai 2020 VBO endpoints; [A] 0.7 conduction partition, linear alloy VBO transfer, no thermal expansion")

_ORIENTATION_FACTORS = {
    "c_plane": 1.0,
    "semipolar_11_22": 0.2,
    "m_plane": 0.0,
    "a_plane": 0.0,
}

def orientation_factor(orientation, polarization_factor=None):
    """Return the reduced-model normal-polarization factor.

    The c-plane scalar strain and band-edge partition (`band_edges` above)
    are retained for every orientation [A/E transfer; Bernardini et al., PRB
    1997; Rinke et al., PRB 2008] -- this module does not itself rotate the
    strain tensor or the valence partition by orientation.  Growth-axis
    MASSES are handled separately, in `nitride_levels._growth_masses`: for
    c_plane and semipolar_11_22 growth the confinement axis is the crystal
    c-axis (GaN me_z=0.209, mh_z=1.88), while for nonpolar m/a-plane growth
    the confinement axis is perpendicular to c and nitride_levels swaps in
    the PERPENDICULAR (xy) masses (me_xy, mh_xy) for z-confinement and the
    c-axis masses for the lateral direction instead -- since commit 17a333f
    (hardening-round fix; previously this swap was missing and every
    orientation used c-plane masses throughout, differing by up to about
    5.7x in the hole mass) [A; Rinke et al., PRB 77, 075202 (2008)].
    Nonpolar valence-band ordering itself is still not modelled
    [A; Schade et al., phys. status solidi (b) (2011)].  Factors 1, 0.2, and
    0 are assumptions,
    not measurements.  Schade et al., phys. status solidi (b) (2011),
    discusses orientation-dependent band structure and matrix elements that
    this scalar model does not rotate or reproduce.  A zero normal component
    here does not establish that a finite real dot has no lateral fields.
    """
    if not isinstance(orientation, str) or orientation not in _ORIENTATION_FACTORS:
        raise ValueError("orientation must be c_plane, semipolar_11_22, m_plane, or a_plane")
    if polarization_factor is None:
        return _ORIENTATION_FACTORS[orientation]
    if isinstance(polarization_factor, bool) or not isinstance(polarization_factor, (int, float)):
        raise ValueError("polarization_factor must be a finite real number")
    factor = float(polarization_factor)
    if not math.isfinite(factor):
        raise ValueError("polarization_factor must be finite")
    if orientation == "semipolar_11_22":
        if not 0.0 <= factor <= 1.0:
            raise ValueError("semipolar polarization_factor must be in [0, 1]")
        return factor
    if factor != _ORIENTATION_FACTORS[orientation]:
        raise ValueError("polarization_factor contradicts the selected orientation")
    return factor

def polarization_field(dot,matrix,T_K,*,strain_fraction=1.,screening_fraction=0.,external_field_kVcm=0.,orientation="c_plane",polarization_factor=None):
    """Normal polarization field plus an independently applied junction field.

    The returned value is the c-axis projection with the module's "+c
    positive" sign convention.  Orientation factors are non-negative and do
    not change this polarity; nitride_levels._z_potential and
    device._evaluate_nitride therefore add the depletion field directly.
    The single scalar factor multiplies the summed spontaneous and
    piezoelectric discontinuity [A].  At (11-22) those contributions partially
    cancel, so sign-reversed semipolar components are outside the [0,1]
    clamp by construction [A; Romanov et al., J. Appl. Phys. 100, 023522
    (2006)].

    Electrostatic normalization follows Bernardini et al., PRB 56, R10024
    (1997), Table II and Bernardini & Fiorentini, phys. status solidi (b)
    216, 391 (1999), Sec. III/Eqs. 7-8 [V].  Orientation factors are [A],
    as documented by :func:`orientation_factor`; screening is an independent
    scenario parameter and is not a current-dependent law.
    """
    if not 0 <= screening_fraction <= 1: raise ValueError("screening_fraction must be in [0, 1]")
    pd=band_edges(dot,T_K,substrate=matrix,strain_fraction=strain_fraction)["P_total_Cm2"]
    pm=band_edges(matrix,T_K,substrate=matrix,strain_fraction=0.)["P_total_Cm2"]
    # Preserve the pre-orientation c-plane arithmetic path byte-for-byte in
    # numerical operation order for legacy callers and explicit factor=1.
    factor=orientation_factor(orientation,polarization_factor)
    if orientation == "c_plane" and factor == 1.0:
        return (1-screening_fraction)*(pm-pd)/(EPS0_SI*dot.eps_r)*1e-5+external_field_kVcm # fixed-D thick GaN reservoirs [A]
    intrinsic=factor*(1-screening_fraction)*(pm-pd)/(EPS0_SI*dot.eps_r)*1e-5
    return intrinsic+external_field_kVcm # external junction field is never orientation/screening scaled [A]
