"""Finite-well c-plane InGaN/GaN dot levels.

This is a single-band, anisotropic-effective-mass, separable disk model [E],
with thick GaN reservoirs and a real-energy (not tunnelling) retention
treatment [A].  Consequently a bound state in a tilted well may still
tunnel; Arrhenius retention is optimistic in that situation.  Material
response is supplied by :mod:`nitride_materials`.

FIELD GEOMETRY [fix round 2, 2026-09-09]. A single InGaN dot embedded in
thick GaN carries polarization SHEET CHARGES only at its own two z
interfaces.  Those sheet charges are equal and opposite (a dipole layer):
by Gauss's law the field they produce is confined to the region BETWEEN
them (inside the dot) and is exactly zero outside, so the exterior GaN
band edges are FLAT and EQUAL on both sides (the model's documented
"thick GaN reservoir" boundary condition below pins the far bands rather
than letting an isolated dipole's own capacitor-like plateau offset stand,
which is the semiconductor-reservoir-screening picture, not vacuum
electrostatics).  `_z_state` therefore applies field (the caller's
polarization field, PLUS `external_field_kVcm` if supplied -- see its
docstring) only to |z| <= height/2 and holds the exterior flat at the
barrier level, the SAME value on both sides.  The previous version instead
held the tilt's boundary value into the exterior, so the two exterior
plateaus differed by the full interior voltage drop, fabricating a large,
artificial reduction of the barrier on one side -- not a physical
depletion effect, a modelling bug -- which is why almost every unscreened
dot came back 'hole unbound'.  `screening_fraction` remains an explicit
[A] parameter on top of this fix.
"""
from dataclasses import dataclass
from functools import lru_cache
import math
import numpy as np
from scipy.linalg import eigh_tridiagonal
from .nitride_materials import binary, ingaN, band_edges, polarization_field, KB_EV
from .dot_levels import finite_disk_2d

_HBAR2_2M0 = 0.0380998212 # eV nm2 [DR] CODATA 2018 constants
_HBAR2_2M0_MEV = 38.0998212 # meV nm2, same constant [DR] CODATA 2018
_HBAR_SI = 1.054571817e-34 # J s [V] CODATA 2018
_M0 = 9.1093837015e-31 # kg [V] CODATA 2018
_KB_SI = 1.380649e-23 # J K-1 [V] SI 2019
_HC_EV_NM = 1239.841984 # h c [V] CODATA 2018, eV nm
_E2_4PIEPS0_EV_NM = 1.439964 # e^2/(4 pi eps0) [DR] CODATA 2018, eV nm

# Source-transcription targets: named module constants (compared, in the
# verifier, to INDEPENDENTLY typed literals so the check is a real
# transcription cross-check rather than a literal-equals-itself tautology).
DESHPANDE2013_HEIGHT_NM = 2.0         # [V] Deshpande et al., Nat. Commun. 4, 1675 (2013), p.2, Fig. 1d
DESHPANDE2013_DIAMETER_NM = 25.0      # [V] ibid., p.2, Fig. 1c (25 +/- 5 nm)
DESHPANDE2013_X_IN = 0.25             # [V] ibid., p.2 (single-wire EDX)
DESHPANDE2013_X_NM = 436.56           # [V] ibid., p.2-3, Fig. 3c (single-dot EL, X line)
DESHPANDE2013_XX_ANTIBIND_MEV = -10.0 # [V] ibid., p.2-3, Fig. 3c (XX above X: antibinding)
ZHANG2013_HEIGHT_NM = 3.0             # [V] Zhang et al., APL 103, 192114 (2013), Fig. 2
ZHANG2013_DIAMETER_NM = 29.0          # [V] ibid., Fig. 2 (D = 29 nm pillar)
ZHANG2013_X_IN = 0.15                 # [V] ibid., Fig. 2 (In0.15Ga0.85N)
ZHANG2013_ZPL_EV = 2.95               # [V] ibid., Fig. 2 (10 K)
ZHANG2013_LIFETIME_NS = 3.40          # [V] ibid., Fig. 2 (mono-exp lifetime, 10 K)
ZHANG2016_SLOPE_MEV_PER_V = -10.0     # [V] Zhang et al., APL 108, 153102 (2016), Fig. 5 (below 2 V)
INGAN_DOT_LIFETIME_RANGE_NS = (1.0, 10.0) # [E] digest class range, nitride_digests.md Sec.11 item 5
TAU_RAD0_DEFAULT_NS = 1.3              # [E] Deshpande et al., APL 105, 141109 (2014), abstract:
                                        # RT electrically-driven recombination lifetime 1.3+/-0.3 ns;
                                        # transferred here as a field-free radiative-lifetime DEFAULT,
                                        # not a verified measurement of this planar design's own
                                        # baseline (their number mixes radiative and any non-radiative
                                        # loss at RT) -- [E], overridable by callers.

@dataclass(frozen=True)
class NitrideDotSystem:
    height_nm: float = 3.0; radius_nm: float = 10.0; x_in: float = .25
    wl_thickness_nm: float = .5; strain_fraction: float = 1.0
    screening_fraction: float = 0.0; external_field_kVcm: float = 0.0
    vbo_InN_GaN_eV: float = 1.15; strain_c_fraction: float = .7

@dataclass(frozen=True)
class NitrideLevels:
    E_X_eV: float; lambda_nm: float; electron_bound: bool; hole_bound: bool
    overlap_sq: float; field_kVcm: float; E_e_meV: float; E_h_meV: float
    dE_e_meV: float; dE_h_meV: float; dE_pair_meV: object
    sp_split_e_meV: float; sp_split_h_meV: float; m_e_matrix_xy: float
    m_h_matrix_xy: float; valid: bool; invalid_reasons: tuple; provenance: str

def _validate(s):
    bad=[]
    for n in ('height_nm','radius_nm','x_in','wl_thickness_nm','strain_fraction',
              'screening_fraction','strain_c_fraction','external_field_kVcm','vbo_InN_GaN_eV'):
        v=getattr(s,n)
        if not math.isfinite(v): bad.append(n+' is not finite')
    if s.height_nm<=0: bad.append('height_nm must be positive')
    if s.radius_nm<=0: bad.append('radius_nm must be positive')
    if not 0<=s.x_in<=1: bad.append('x_in must be in [0,1]')
    if not 0<=s.wl_thickness_nm<s.height_nm: bad.append('wl_thickness_nm must be >=0 and < height_nm')
    if not 0<=s.strain_fraction<=1: bad.append('strain_fraction must be in [0,1]')
    if not 0<=s.screening_fraction<=1: bad.append('screening_fraction must be in [0,1]')
    return bad

def _z_grid(half, pad, target_n):
    """Node grid with +/-half landing EXACTLY on a node (finite-volume cell
    face), so the mass/potential step at the dot interface is not smeared
    to first order by an arbitrary offset between a linspace node and the
    physical edge. `target_n` sets the INTERIOR half-well resolution
    (n_half = target_n/2, floored at 200): the well is nm-scale while
    `pad` (exterior_nm) is tens of nm, so splitting a single uniform dz
    proportionally across the whole padded domain (the previous scheme)
    starves the interior of points for a thin dot -- e.g. z_points=1201,
    exterior_nm=45 gave only ~7 half-well nodes at height_nm=1, a first-
    order interface error far above the 0.5 meV budget. Fixing dz from the
    interior requirement and extending the SAME dz into the pad instead
    keeps the interface well resolved regardless of pad size."""
    n_half = max(200, int(round(target_n / 2.0)))
    dz = half / n_half
    n_pad = max(1, int(round(pad / dz)))
    idx = np.arange(-(n_half + n_pad), n_half + n_pad + 1)
    return idx * dz, dz, n_half

def _z_state(height, barrier, md, mb, field, sign, n=1201, pad=45.):
    """BDD finite-volume state; field is kV/cm, sign is carrier charge sign.

    `field` tilts the potential ONLY inside the dot (|z| <= height/2); the
    exterior plateau is flat AND EQUAL on both sides (the model's own
    "thick GaN reservoir" assumption pins the far bands, exactly as a
    dipole layer of polarization sheet charge at the two dot interfaces
    produces a field confined between them and none outside). The
    previous version instead HELD the tilt's boundary value into the
    exterior, so the two exterior plateaus differed by the full interior
    voltage drop -- that fabricated a large, artificial reduction of the
    barrier seen by the carrier on one side, which is why every unscreened
    (default screening_fraction=0) dot came back 'hole unbound'. `field`
    is the caller's COMBINED internal-polarization + external-depletion
    field (see NitrideDotSystem.external_field_kVcm): both are folded into
    this same interior-confined tilt rather than the external contribution
    being extended into the exterior reservoir, because the exterior
    padding (`pad`/`exterior_nm`) is a NUMERICAL convergence control, not a
    physical depletion width -- letting any field ramp through it would
    make results depend on the padding choice. This is a documented [A]
    simplification: the true external depletion field would, in principle,
    also weakly tilt the exterior bands over the diode's actual (much
    longer) intrinsic-region length, which this nm-scale dot solver does
    not model. BenDaniel & Duke, PR 152, 683 (1966) [V] supplies the
    flux-matching boundary condition (mass-weighted hopping).
    """
    half = height / 2.
    z, dz, n_half = _z_grid(half, pad, n)
    n_pts = len(z)
    node_idx = np.round(z / dz).astype(int)
    inside = np.abs(node_idx) <= n_half
    mass = np.where(inside, md, mb)
    # e*(kV/cm)*(nm) = 1e-4 eV; sign gives opposite electron/hole tilt;
    # zero outside |z|<=half -> flat, equal exterior plateaus both sides.
    tilt = sign * field * 1e-4 * np.where(inside, z, 0.0)
    V = np.where(inside, 0., barrier) + tilt
    invface=2/(mass[:-1]+mass[1:]); a=_HBAR2_2M0*invface/dz**2
    diag=np.empty(n_pts-2); off=-a[1:-1]
    diag[:] = a[:-1]+a[1:]+V[1:-1]
    vals, vecs=eigh_tridiagonal(diag,off,select='i',select_range=(0,1))
    psi=np.zeros(n_pts); psi[1:-1]=vecs[:,0]; psi/=math.sqrt(np.trapezoid(psi*psi,z))
    cont=min(V[0],V[-1]); loc=float(np.trapezoid(psi[inside]**2,z[inside]))
    e1 = vals[1] if len(vals)>1 else float('nan')
    return vals[0], e1, cont, loc, z, psi

def _radial(remaining_offset, radius, md, mb):
    """Finite circular well (BenDaniel-Duke matching, barrier mass mb),
    reusing dot_levels' material-independent disk kernel [E]; `remaining_offset`
    is the barrier headroom left after the z-confinement energy is removed
    (adiabatic separable decoupling: the disk sees V - E_z, not the full
    offset). Returns (E0_eV, E1_eV or nan, bound, p_bound, rms_r_nm or None)."""
    if remaining_offset <= 0.:
        return float('nan'), float('nan'), False, False, None
    d = finite_disk_2d(remaining_offset * 1000., radius, md, mb)
    if not d.bound:
        return float('nan'), float('nan'), False, False, None
    e0 = d.E0_meV / 1000.
    e1 = d.E1_meV / 1000. if d.p_bound else float('nan')
    return e0, e1, True, d.p_bound, d.rms_r_nm

def _coulomb_binding_eV(l_e_xy, l_h_xy, z_sep_nm, eps_r):
    """Screened e-h Coulomb binding, frozen-orbital Gaussian envelopes [E]
    (same in-plane method as dot_levels._gauss_binding, Bernardini/BenDaniel
    geometry aside) with the field-driven z separation added in quadrature
    to the relative-coordinate width, so QCSE separation SUPPRESSES (never
    hard-zeroes) the binding as the carriers are pulled apart:
        E_b = e^2/(4 pi eps0 eps_r) * sqrt(pi) / sqrt(l_e^2 + l_h^2 + 2 z_sep^2)
    Evaluated from the normalized envelopes actually solved above, not
    fitted to any target wavelength."""
    L2 = l_e_xy*l_e_xy + l_h_xy*l_h_xy + 2.*z_sep_nm*z_sep_nm
    if L2 <= 0.: return float('nan')
    return math.sqrt(math.pi) * _E2_4PIEPS0_EV_NM / (eps_r * math.sqrt(L2))

@lru_cache(maxsize=256)
def _levels_cached(s,T_K,n,pad):
    bad=_validate(s)
    if not math.isfinite(T_K) or T_K<=0: bad.append('T_K must be positive and finite')
    if bad: return _invalid(s,bad)
    d,m=ingaN(s.x_in),binary('GaN')
    de=band_edges(d,T_K,substrate=m,strain_fraction=s.strain_fraction,vbo_InN_GaN_eV=s.vbo_InN_GaN_eV,strain_c_fraction=s.strain_c_fraction)
    be=band_edges(m,T_K,substrate=m)
    Ve=be['Ec_eV']-de['Ec_eV']; Vh=de['Ev_eV']-be['Ev_eV']
    F=polarization_field(d,m,T_K,strain_fraction=s.strain_fraction,screening_fraction=s.screening_fraction,external_field_kVcm=s.external_field_kVcm)
    if Ve<=0: bad.append('electron offset is nonpositive')
    if Vh<=0: bad.append('hole offset is nonpositive')
    if bad: return _invalid(s,bad,F)
    ee,ee1,ce,le,ze,pe=_z_state(s.height_nm,Ve,d.me_z,m.me_z,F,-1,n,pad)
    eh,eh1,ch,lh,zh,ph=_z_state(s.height_nm,Vh,d.mh_z,m.mh_z,F,+1,n,pad)
    re,rpe,ok_e,pb_e,rms_e=_radial(ce-ee,s.radius_nm,d.me_xy,m.me_xy)
    rh,rph,ok_h,pb_h,rms_h=_radial(ch-eh,s.radius_nm,d.mh_xy,m.mh_xy)
    eb=ee+(re if ok_e else float('nan')); hb=eh+(rh if ok_h else float('nan'))
    if not (le>.50 and ok_e and eb<ce): bad.append('electron unbound or laterally exhausted offset')
    if not (lh>.50 and ok_h and hb<ch): bad.append('hole unbound or laterally exhausted offset')
    if bad: return _invalid(s,bad,F)
    # Envelopes are separable; normalized radial ground envelopes cancel in
    # overlap ratio under the common-disk approximation [E].
    ov=float(np.trapezoid(pe*ph,ze)**2); ov=max(0.,min(1.,ov))
    eps=(d.eps_r+m.eps_r)/2
    z_sep=abs(float(np.trapezoid(ze*pe*pe,ze))-float(np.trapezoid(zh*ph*ph,zh)))
    l_e_xy = rms_e if rms_e else s.radius_nm; l_h_xy = rms_h if rms_h else s.radius_nm
    coul=_coulomb_binding_eV(l_e_xy,l_h_xy,z_sep,eps) # eV, screened Gaussian-envelope proxy [E]
    ex=(de['Ec_eV']-de['Ev_eV'])+eb+hb-coul
    # First excited state = min(z-excitation, radial p-shell excitation),
    # each admitted only if it is itself bound below the local continuum.
    zg_e = (ee1-ee) if (math.isfinite(ee1) and ee1<ce) else float('inf')
    rg_e = (rpe-re) if (ok_e and math.isfinite(rpe)) else float('inf')
    sp_e = min(zg_e,rg_e)
    zg_h = (eh1-eh) if (math.isfinite(eh1) and eh1<ch) else float('inf')
    rg_h = (rph-rh) if (ok_h and math.isfinite(rph)) else float('inf')
    sp_h = min(zg_h,rg_h)
    valid=True
    return NitrideLevels(ex,_HC_EV_NM/ex if ex>0 else float('nan'),True,True,ov,F,
                          eb*1000,hb*1000,(ce-eb)*1000,(ch-hb)*1000,None,
                          sp_e*1000 if math.isfinite(sp_e) else float('nan'),
                          sp_h*1000 if math.isfinite(sp_h) else float('nan'),
                          m.me_xy,m.mh_xy,valid,tuple(bad),
                          '[V] Rinke et al., PRB 77, 075202 (2008), Table V; '
                          '[V] Bernardini et al., PRB 56, R10024 (1997); '
                          '[V] BenDaniel & Duke, PR 152, 683 (1966); '
                          '[E] separable disk, finite-barrier radial confinement, and screened '
                          'Gaussian-envelope Coulomb approximation; '
                          '[A] thick GaN reservoirs (flat exterior band edges outside the dot), '
                          'real-energy retention')

def _invalid(s,reasons,F=0.):
    return NitrideLevels(float('nan'),float('nan'),False,False,0.,F,float('nan'),float('nan'),float('nan'),float('nan'),None,float('nan'),float('nan'),float('nan'),float('nan'),False,tuple(reasons),'[A] invalid geometry/offset rejected before model evaluation')

def levels(system, T_K=300.0, *, z_points=1201, exterior_nm=45.0):
    """Return immutable cached levels. Numerical controls are cache keys."""
    if not isinstance(system,NitrideDotSystem): raise TypeError('system must be NitrideDotSystem')
    return _levels_cached(system,float(T_K),int(z_points),float(exterior_nm))

def rates(lv,T_K,*,tau_rad0_ns=TAU_RAD0_DEFAULT_NS,n_dot_cm2=1e10,tau_cap_ps=10.,tau_cap_scales_with_density=False,channel='min',k_nr_ns=0.):
    """Absolute detailed-balance escape rates in 1/ns; no cavity/Purcell input.

    tau_rad0_ns default is TAU_RAD0_DEFAULT_NS [E], see module docstring
    constants; callers (fsim_core.device) always pass an explicit value."""
    bad=list(lv.invalid_reasons)
    if not lv.valid: bad.append('levels are invalid')
    if T_K<=0 or tau_rad0_ns<=0 or n_dot_cm2<=0 or tau_cap_ps<=0 or k_nr_ns<0: bad.append('invalid rate input')
    if channel not in ('min','electron','hole','pair','pair_half'): bad.append('unknown channel')
    if channel in ('pair','pair_half'): bad.append('pair channel requires a bound wetting-layer continuum')
    if bad: return dict(gamma_X0_ns=float('nan'),gamma_XX0_ns=float('nan'),k_X_ns=float('nan'),k_XX_ns=float('nan'),escape_prefactor_ns=float('nan'),E_a_meV=float('nan'),S0=float('nan'),tau_cap_ps_used=float('nan'),valid=False,invalid_reasons=tuple(bad),provenance='[A] invalid rate request')
    tc=tau_cap_ps*(1e10/n_dot_cm2) if tau_cap_scales_with_density else tau_cap_ps
    choices={'electron':(lv.dE_e_meV,lv.m_e_matrix_xy),'hole':(lv.dE_h_meV,lv.m_h_matrix_xy)}
    name=min(choices,key=lambda q:choices[q][0]) if channel=='min' else channel
    ea,mass=choices[name]
    if ea<0: return dict(gamma_X0_ns=float('nan'),gamma_XX0_ns=float('nan'),k_X_ns=float('nan'),k_XX_ns=float('nan'),escape_prefactor_ns=float('nan'),E_a_meV=ea,S0=float('nan'),tau_cap_ps_used=tc,valid=False,invalid_reasons=('negative escape barrier',),provenance='[A] negative barrier is not floored')
    n2d=mass*_M0*_KB_SI*T_K/(math.pi*_HBAR_SI**2)/1e4
    attempt=(1000/tc)*(n2d/n_dot_cm2)
    k=attempt*math.exp(-(ea/1000)/(KB_EV*T_K))+k_nr_ns
    gx=lv.overlap_sq/tau_rad0_ns; gxx=2*gx
    return dict(gamma_X0_ns=gx,gamma_XX0_ns=gxx,k_X_ns=k,k_XX_ns=2*k,escape_prefactor_ns=attempt,E_a_meV=ea,S0=gx/(gx+k),tau_cap_ps_used=tc,valid=True,invalid_reasons=(),provenance='[A] 2D-reservoir detailed balance N2D=m kT/(pi hbar2); [A] XX=2X cascade transfer, Reischle et al., OE 16, 12771 (2008); selected '+name+' channel')
