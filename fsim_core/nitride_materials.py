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
    # Anisotropic A-gap deformation potentials (eV), Yan et al. PRB 90,
    # 125118 (2014) Eqs. (2)-(3) convention; see strained_gap_shift_eV.
    # acz_D1_eV is acz - D1 (= VM03 a1), act_D2_eV is act - D2 (= VM03 a2).
    acz_D1_eV: float; act_D2_eV: float; D3_eV: float; D4_eV: float

# Bernardini, Fiorentini & Vanderbilt, PRB 56, R10024 (1997), Table II [V]
# (Psp, e31, e33).  Static eps_r = 10.28/14.61/10.31 for GaN/InN/AlN is
# Bernardini & Fiorentini, phys. stat. sol. (b) 216, 391 (1999), Sec. III [V]:
# one internally-consistent dielectric family used with the SAME polarization
# formula below, not Ioffe's mixed perpendicular/parallel/static set.
# Rinke et al., PRB 77, 075202 (2008), Table V and Sec. IV.B [V] (masses).
# Electron-mass AXES (audit fix 2026-09-23; arXiv:0801.0421 Table V rows
# m_e-par / m_e-perp, axis fixed by App. A Eq. (A5), where 1/m_e-par
# multiplies kz^2): m_e along c (me_z) = 0.186 / 0.065 / 0.322 and m_e in the
# c plane (me_xy) = 0.209 / 0.068 / 0.329 for GaN / InN / AlN [V].  Before
# 2026-09-23 the two electron axes were swapped here; the hole masses
# (m0/mh_z = -(A1+A3), c axis) were and are on the correct axis [V; Rinke
# 2008 App. B].
# Gap deformation potentials acz-D1, act-D2, D3, D4: Yan, Rinke, Janotti,
# Scheffler & Van de Walle, PRB 90, 125118 (2014), Table III "(recommended)"
# rows [V] (GaN HSE alpha=0.30, InN HSE alpha=0.25, AlN HSE alpha=0.34).
# The Vurgaftman & Meyer, JAP 94, 3675 (2003) set (GaN a1 -4.90, a2 -11.30,
# D3 8.20, D4 -4.10, reused for InN) is known here only via Yan 2014's
# reprint [E] and is NOT used.
# Lattice and elastic constants: Ioffe NSM summary [V]; InN elastic choice is
# the Sheleg family [E].  Varshni pairs: Vurgaftman & Meyer, JAP 94, 3675
# (2003) framework [E] pending primary-text confirmation.
_B = {
 "GaN": NitrideMaterial("GaN",0.,3.189,5.186,106.,398.,-.49,.73,-.029,10.28,.186,.209,1.88,.33,"[V] Bernardini PRB 1997 Table II polarization; Bernardini & Fiorentini pss(b) 1999 Sec. III eps_r; Rinke PRB 2008 Table V/Sec IV.B masses; Yan PRB 2014 Table III deformation potentials",3.51,9.09e-4,830.,-7.6,-6.07,-8.88,5.38,-2.69),
 "InN": NitrideMaterial("InN",1.,3.545,5.703,121.,182.,-.57,.97,-.032,14.61,.065,.068,1.63,1.63,"[V] Bernardini PRB 1997 Table II polarization; Bernardini & Fiorentini pss(b) 1999 Sec. III eps_r; Rinke PRB 2008 electron masses [V]; Ioffe NSM hole mass 1.63 applied isotropically [E]; Yan PRB 2014 Table III deformation potentials [V]",.78,2.45e-4,624.,-4.2,-3.64,-4.58,2.68,-1.78),
 "AlN": NitrideMaterial("AlN",0.,3.112,4.982,99.,389.,-.60,1.46,-.081,10.31,.322,.329,3.53,10.4,"[V] Bernardini PRB 1997 Table II polarization; Bernardini & Fiorentini pss(b) 1999 Sec. III eps_r; Rinke PRB 2008 electron masses; Ioffe NSM hole masses 3.53/10.4; Yan PRB 2014 Table III deformation potentials [V]",6.25,1.799e-3,1462.,-9.8,-4.36,-12.35,9.17,-3.72),
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
    """In_xGa_(1-x)N virtual crystal; 1.4 eV bowing is Wu et al. [V].

    Every field, including the Yan 2014 deformation potentials [V] and the
    Rinke 2008 a_V [V], is linear in x between the GaN and InN endpoints
    (Vegard-type virtual crystal, no bowing of the potentials) [A]; only
    the gap itself carries the Wu et al. bowing.
    """
    if not math.isfinite(x_in) or not 0. <= x_in <= 1.: raise ValueError("x_in must be finite and in [0, 1]")
    if x_in == 0.: return binary("GaN")
    if x_in == 1.: return binary("InN")
    g,i=_B["GaN"],_B["InN"]; q=lambda k:(1-x_in)*getattr(g,k)+x_in*getattr(i,k)
    return NitrideMaterial("InGaN",x_in,*[q(k) for k in ("a_A","c_A","C13_GPa","C33_GPa","e31_Cm2","e33_Cm2","Psp_Cm2","eps_r","me_z","me_xy","mh_z","mh_xy")],"[A] linear virtual-crystal alloy; Wu et al. bowing [V]",q("Eg0_eV"),q("alpha_eVK"),q("beta_K"),q("aV_eV"),*[q(k) for k in ("acz_D1_eV","act_D2_eV","D3_eV","D4_eV")])

def bandgap(material, T_K):
    if not math.isfinite(T_K) or T_K <= 0: raise ValueError("T_K must be positive and finite")
    return material.Eg0_eV-material.alpha_eVK*T_K*T_K/(T_K+material.beta_K)-1.4*material.x_in*(1-material.x_in) # Wu et al., APL 2002, reported by OSTI 835986 [V]

def strained_gap_shift_eV(material, eps_parallel, eps_zz):
    """Strain shift of the conduction-to-A-valence (HH) gap, eV.

    Wurtzite Bir-Pikus / k.p strain Hamiltonian (Chuang & Chang, PRB 54,
    2491 (1996)) in the form and sign convention of Yan et al., PRB 90,
    125118 (2014), Eqs. (2)-(3) [V]:

        dE_CB = acz*ezz + act*eperp
        dE_A  = (D1 + D3)*ezz + (D2 + D4)*eperp        (A = HH valence edge)
        dEg_A = dE_CB - dE_A
              = (acz - D1 - D3)*ezz + (act - D2 - D4)*eperp,

    with eperp = exx + eyy = 2*eps_parallel for c-plane biaxial strain and
    ezz along c; dEg_A > 0 opens the gap.  acz-D1, act-D2, D3, D4 are the
    Yan 2014 Table III recommended values [V] (only the gap combinations
    acz-D1 and act-D2 are published there, not acz, act, D1, D2
    separately).  Hydrostatic strain (exx=eyy=ezz=e) gives
    dEg = [(acz-D1-D3) + 2(act-D2-D4)] * e, i.e. an effective volume
    potential a_V,eff = [(acz-D1-D3) + 2(act-D2-D4)]/3 = -7.94 eV for GaN
    [DR], consistent with the independent Rinke et al., PRB 77, 075202
    (2008) Table IV a_V = -7.6 eV [V] (kept as `aV_eV`, cross-check only).
    Shear strains (D5, D6) vanish for c-plane biaxial strain.

    Validity caveat [A]: the -(D3*ezz + D4*eperp) valence term is the shift
    of the A (heavy-hole, Gamma9) edge, so this form is the strained GAP only
    while HH is the top valence band (GaN, InN and InGaN).  It does not hold
    for AlN, whose top valence band is the crystal-field split-off (CH)
    band (negative crystal-field splitting); there the gap shift would be
    (acz-D1)*ezz + (act-D2)*eperp with no D3/D4 term.  The AlN literals are
    carried for completeness and the module applies this form unchanged
    for strained AlN.
    """
    return ((material.acz_D1_eV-material.D3_eV)*eps_zz
            +(material.act_D2_eV-material.D4_eV)*2.*eps_parallel)

def band_edges(material, T_K, *, substrate=None, strain_fraction=1., vbo_InN_GaN_eV=1.15, strain_c_fraction=.7):
    """Strained c-plane band edges (eV) of `material` coherent on `substrate`.

    Biaxial strain ezz = -2 C13/C33 eps_parallel [V; Yan et al. PRB 2014
    Eq. (9)].  The A-gap strain shift is the anisotropic Yan 2014 form
    (strained_gap_shift_eV) since the 2026-09-23 audit fix H7; before that
    it was the volume-only a_V*(2 eps_parallel + ezz), which overstated the
    x=0.25 gap opening (+249.6 meV against +114.0 meV).  Because the
    published potentials give only the gap combinations, the split of dEg
    between Ec and Ev stays the explicit strain_c_fraction partition [A]
    (0.7 to the conduction band by default); Ec-Ev = Eg + dEg for every
    partition.
    """
    if not 0 <= strain_fraction <= 1 or not 0 <= strain_c_fraction <= 1: raise ValueError("strain fractions must be in [0, 1]")
    sub = binary("GaN") if substrate is None else substrate
    ep=(sub.a_A-material.a_A)/material.a_A*strain_fraction
    ez=-2*material.C13_GPa/material.C33_GPa*ep
    p=material.Psp_Cm2+2*material.e31_Cm2*ep+material.e33_Cm2*ez
    dEg=strained_gap_shift_eV(material,ep,ez) # Yan et al., PRB 90, 125118 (2014), Eq. (3), Table III [V]; Ec/Ev partition below [A]
    # Valence reference: the InN/GaN family uses the linear x_in transfer of
    # the 1.15 eV Tsai & Bayram VBO [V]/[A]; AlN is off that axis and uses
    # its own Tsai & Bayram Table 2 literal [V] instead (see _VBO_ALN_EV).
    ev0 = _VBO_ALN_EV if material.name == "AlN" else material.x_in*vbo_InN_GaN_eV
    ev=ev0-(1-strain_c_fraction)*dEg # VBO Tsai & Bayram ACS Omega 2020 Table 2 [V]; partition [A]
    # Both strained edges are constructed from the same unshifted reference.
    # Therefore Ec-Ev = Eg+dEg for every partition, as required by the volume
    # gap deformation; using `ev` here would apply the valence share to
    # Ec a second time. [DR] Yan et al., PRB 90, 125118 (2014), Eq. (3).
    ec=ev0+bandgap(material,T_K)+strain_c_fraction*dEg
    return dict(Ec_eV=ec,Ev_eV=ev,eps_parallel=ep,eps_zz=ez,P_total_Cm2=p,provenance="[V] B97 polarization; [V] Yan 2014 anisotropic A-gap deformation potentials; [V] Tsai 2020 VBO endpoints; [A] 0.7 conduction partition, linear alloy VBO and deformation-potential transfer, no thermal expansion")

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
    c-axis (GaN me_z=0.186, mh_z=1.88; me_z read 0.209 before the
    2026-09-23 electron-axis fix), while for nonpolar m/a-plane growth
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
