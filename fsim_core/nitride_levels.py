"""Finite-well c-plane InGaN/GaN dot levels.

This is a single-band, separable disk model [E], with thick GaN reservoirs
and a real-energy (not tunnelling) retention treatment [A].  Consequently a
bound state in a tilted well may still tunnel; Arrhenius retention is optimistic
in that situation.  Material response is supplied by :mod:`nitride_materials`.
"""
from dataclasses import dataclass
from functools import lru_cache
import math
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.special import jn_zeros
from .nitride_materials import binary, ingaN, band_edges, polarization_field, KB_EV

_HBAR2_2M0 = 0.0380998212 # eV nm2 [DR] CODATA 2018 constants
_HBAR_SI = 1.054571817e-34 # J s [V] CODATA 2018
_M0 = 9.1093837015e-31 # kg [V] CODATA 2018
_KB_SI = 1.380649e-23 # J K-1 [V] SI 2019

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
    for n in ('height_nm','radius_nm','x_in','wl_thickness_nm','strain_fraction','screening_fraction','strain_c_fraction'):
        v=getattr(s,n)
        if not math.isfinite(v): bad.append(n+' is not finite')
    if s.height_nm<=0: bad.append('height_nm must be positive')
    if s.radius_nm<=0: bad.append('radius_nm must be positive')
    if not 0<=s.x_in<=1: bad.append('x_in must be in [0,1]')
    if not 0<=s.wl_thickness_nm<s.height_nm: bad.append('wl_thickness_nm must be >=0 and < height_nm')
    if not 0<=s.strain_fraction<=1: bad.append('strain_fraction must be in [0,1]')
    if not 0<=s.screening_fraction<=1: bad.append('screening_fraction must be in [0,1]')
    return bad

def _z_state(height, barrier, md, mb, field, sign, n=801, pad=25.):
    """BDD finite-volume state; field is kV/cm, sign is carrier charge sign.

    The exterior potential is held at its interface value, so no artificial
    infinite ramp/box continuum is introduced. BenDaniel & Duke, PR 152, 683
    (1966) [V] supplies the flux-matching boundary condition.
    """
    half=height/2; z=np.linspace(-half-pad,half+pad,n); dz=z[1]-z[0]
    inside=np.abs(z)<=half
    mass=np.where(inside,md,mb)
    # e*(kV/cm)*(nm) = 1e-4 eV; sign gives opposite electron/hole tilt.
    tilt=sign*field*1e-4*np.clip(z,-half,half)
    V=np.where(inside,0.,barrier)+tilt
    invface=2/(mass[:-1]+mass[1:]); a=_HBAR2_2M0*invface/dz**2
    diag=np.empty(n-2); off=-a[1:-1]
    diag[:] = a[:-1]+a[1:]+V[1:-1]
    vals, vecs=eigh_tridiagonal(diag,off,select='i',select_range=(0,1))
    psi=np.zeros(n); psi[1:-1]=vecs[:,0]; psi/=math.sqrt(np.trapezoid(psi*psi,z))
    cont=min(V[0],V[-1]); loc=float(np.trapezoid(psi[inside]**2,z[inside]))
    return vals[0], vals[1] if len(vals)>1 else float('nan'), cont, loc, z, psi

def _radial(radius, barrier, md, mb):
    # Infinite-cylinder Bessel energies are an upper-bound separable estimate
    # [E]; reject when this already exhausts the finite lateral offset.
    a0,a1=jn_zeros(0,2); e0=_HBAR2_2M0*a0*a0/(md*radius*radius)
    e1=_HBAR2_2M0*a1*a1/(md*radius*radius)
    return e0,e1, e0 < barrier

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
    re,rpe,ok_e=_radial(s.radius_nm,Ve,d.me_xy,m.me_xy); rh,rph,ok_h=_radial(s.radius_nm,Vh,d.mh_xy,m.mh_xy)
    eb=ee+re; hb=eh+rh
    if not (eb<ce and le>.50 and ok_e): bad.append('electron unbound or laterally exhausted offset')
    if not (hb<ch and lh>.50 and ok_h): bad.append('hole unbound or laterally exhausted offset')
    # Envelopes are separable; normalized radial ground envelopes cancel in
    # overlap ratio under the common-disk approximation [E].
    ov=float(np.trapezoid(pe*ph,ze)**2); ov=max(0.,min(1.,ov))
    eps=(d.eps_r+m.eps_r)/2
    # Gaussian-envelope Coulomb estimate [E/A], not fitted to any wavelength.
    r_eff=math.sqrt(s.radius_nm*s.radius_nm+(s.height_nm/2)**2)
    coul=1.439964/(eps*r_eff)*ov # eV, screened point/Gaussian proxy [E]
    ex=band_edges(d,T_K,substrate=m,strain_fraction=s.strain_fraction,vbo_InN_GaN_eV=s.vbo_InN_GaN_eV,strain_c_fraction=s.strain_c_fraction)['Ec_eV']-de['Ev_eV']+eb+hb-coul
    valid=not bad
    return NitrideLevels(ex,1239.841984/ex if ex>0 else float('nan'),not any('electron' in x for x in bad),not any('hole' in x for x in bad),ov,F,eb*1000,hb*1000,(ce-eb)*1000,(ch-hb)*1000,None,(ee1-ee+rpe-re)*1000,(eh1-eh+rph-rh)*1000,m.me_xy,m.mh_xy,valid,tuple(bad),'[V] Rinke et al., PRB 77, 075202 (2008), Table V; [V] Bernardini et al., PRB 56, R10024 (1997); [V] BenDaniel & Duke, PR 152, 683 (1966); [E] separable disk and screened Coulomb approximation; [A] thick GaN reservoirs, real-energy retention')

def _invalid(s,reasons,F=0.):
    return NitrideLevels(float('nan'),float('nan'),False,False,0.,F,float('nan'),float('nan'),float('nan'),float('nan'),None,float('nan'),float('nan'),float('nan'),float('nan'),False,tuple(reasons),'[A] invalid geometry/offset rejected before model evaluation')

def levels(system, T_K=300.0, *, z_points=1201, exterior_nm=45.0):
    """Return immutable cached levels. Numerical controls are cache keys."""
    if not isinstance(system,NitrideDotSystem): raise TypeError('system must be NitrideDotSystem')
    return _levels_cached(system,float(T_K),int(z_points),float(exterior_nm))

def rates(lv,T_K,*,tau_rad0_ns=1.,n_dot_cm2=1e10,tau_cap_ps=10.,tau_cap_scales_with_density=False,channel='min',k_nr_ns=0.):
    """Absolute detailed-balance escape rates in 1/ns; no cavity/Purcell input."""
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
