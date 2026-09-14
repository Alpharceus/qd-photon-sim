"""C-axis InGaN disc-in-GaN-wire levels; hard-wall radial envelope [A]."""
from dataclasses import dataclass
from functools import lru_cache
import math
import numpy as np
from scipy.special import jn_zeros,jv
from .nitride_materials import EPS0_SI,band_edges,binary,ingaN,KB_EV
from .nitride_levels import _z_state,_coulomb_binding_eV
HB=1.054571817e-34 # [V] CODATA 2018
M0=9.1093837015e-31 # [V] CODATA 2018
KB=1.380649e-23 # [V] SI 2019
HC=1239.841984 # [V] CODATA 2018 eV nm
J01=2.4048255577;J11=3.8317059702 # [V] Abramowitz and Stegun 1964 Table 9.5
@dataclass(frozen=True)
class NitrideNanowireSystem:
 height_nm:float=2.0 # [V] Deshpande et al., Nat Commun 2013 Fig. 1
 core_radius_nm:float=12.5 # [V] Deshpande et al., Nat Commun 2013, 25 nm diameter
 outer_radius_nm:float=12.5 # [A] no shell in levels model
 x_in:float=.4 # [V] Deshpande et al., APL 2014 abstract
 strain_bound:str='relaxed' # [A] scenario label
 screening_fraction:float=0. # [A] independent polarization screening
 external_field_kVcm:float=0. # [A] independent junction field
 vbo_InN_GaN_eV:float=1.15 # [V] Tsai and Bayram, ACS Omega 2020 Table 2
 strain_c_fraction:float=.7 # [A] deformation-potential partition
@dataclass(frozen=True)
class NanowireLevels:
 E_X_eV:float;lambda_nm:float;electron_bound:bool;hole_bound:bool;overlap_sq:float;field_kVcm:float;F_sp_kVcm:float;F_pz_kVcm:float;P_sp_dot_Cm2:float;P_sp_GaN_Cm2:float;P_pz_dot_Cm2:float;P_pz_GaN_Cm2:float;E_e_meV:float;E_h_meV:float;dE_e_meV:float;dE_h_meV:float;dE_pair_meV:object;sp_split_e_meV:float;sp_split_h_meV:float;m_e_matrix_xy:float;m_h_matrix_xy:float;axial_threshold_e_meV:float;axial_threshold_h_meV:float;transverse_e_meV:float;transverse_h_meV:float;rms_radius_nm:float;core_radius_nm:float;outer_radius_nm:float;geometry:str;approximation_metadata:str;valid:bool;invalid_reasons:tuple;provenance:str
def _ep(m,r,j=J01): return HB*HB*j*j/(2*m*M0*(r*1e-9)**2)/1.602176634e-22 # [DR] cylinder Bessel mode
@lru_cache(maxsize=1)
def _rf():
 u=np.linspace(0,1,4001);w=u*jv(0,J01*u)**2;return math.sqrt(np.trapezoid(u*u*w,u)/np.trapezoid(w,u))
def _bad(b,F=0.):
 n=float('nan');return NanowireLevels(n,n,False,False,0,F,n,n,n,n,n,n,n,n,n,n,None,n,n,n,n,n,n,n,n,n,n,n,'full-core hard-wall cylinder','[A] rejected state',False,tuple(b),'[A] invalid state')
def _once(s,T,n,pad):
 d,g=ingaN(s.x_in),binary('GaN');sf=0 if s.strain_bound=='relaxed' else 1 # [DR] relaxation removes strain shifts and P_pz
 de=band_edges(d,T,substrate=g,strain_fraction=sf,vbo_InN_GaN_eV=s.vbo_InN_GaN_eV,strain_c_fraction=s.strain_c_fraction);ge=band_edges(g,T,substrate=g)
 pp=de['P_total_Cm2']-d.Psp_Cm2;fs=(g.Psp_Cm2-d.Psp_Cm2)/(EPS0_SI*d.eps_r)*1e-5;fp=-pp/(EPS0_SI*d.eps_r)*1e-5;F=(1-s.screening_fraction)*(fs+fp)+s.external_field_kVcm # [V/DR] Bernardini 1997; [A] external unscreened
 te,Te=_ep(d.me_xy,s.core_radius_nm)/1000,_ep(g.me_xy,s.core_radius_nm)/1000;th,Th=_ep(d.mh_xy,s.core_radius_nm)/1000,_ep(g.mh_xy,s.core_radius_nm)/1000
 ve=ge['Ec_eV']-de['Ec_eV']+Te-te;vh=de['Ev_eV']-ge['Ev_eV']+Th-th
 if ve<=0 or vh<=0:return None,['nonpositive ground-channel axial barrier'],F
 ee,e1,ce,_,z,pe=_z_state(s.height_nm,ve,d.me_z,g.me_z,F,-1,n,pad);eh,h1,ch,_,zh,ph=_z_state(s.height_nm,vh,d.mh_z,g.mh_z,F,1,n,pad)
 if not ee<ce or not eh<ch:return None,['axial state unbound'],F
 ov=max(0,min(1,float(np.trapezoid(pe*ph,z)**2)))
 if ov<1e-12:return None,['overlap_unresolved'],F
 sep=abs(np.trapezoid(z*pe*pe,z)-np.trapezoid(zh*ph*ph,zh));r=s.core_radius_nm*_rf();ex=de['Ec_eV']-de['Ev_eV']+te+th+ee+eh-_coulomb_binding_eV(r,r,sep,d.eps_r) # [E] dot eps envelope
 bare_e=ge['Ec_eV']-de['Ec_eV'];bare_h=de['Ev_eV']-ge['Ev_eV'] # [DR] bare band offset, no transverse zero-point term
 def spl(mxy_d,mxy_b,t,bare,e,eone,c):
  # J11 transverse level [fix A]: this is a separable product-state model
  # (V(z) does not depend on the transverse quantum numbers), so the SAME
  # axial eigenstate `e` applies whichever transverse mode is occupied and
  # the J01->J11 spacing is the pure in-plane Bessel difference t1-t, using
  # the IN-PLANE masses (mxy_d/mxy_b) exactly like the ground subband `t`
  # above -- never the axial masses. Boundedness is checked against the
  # J11-SPECIFIC continuum `bare+T1` (the bare band offset plus the J11
  # transverse zero-point energy), because a transverse-excited carrier
  # autoionizes through its own, higher-energy escape channel even while
  # the ground-transverse continuum `c` (=ve/vh) remains closed. A
  # non-positive or unbound spacing means no bound J11 excited state at
  # this geometry: return nan, never a negative meV value with valid=True.
  t1=_ep(mxy_d,s.core_radius_nm,J11)/1000;T1=_ep(mxy_b,s.core_radius_nm,J11)/1000
  rad=t1-t if bare+T1-t1>e else float('inf')
  ax=eone-e if math.isfinite(eone) and eone<c else float('inf');sp=min(ax,rad)
  return sp if sp>0 else float('nan')
 se=spl(d.me_xy,g.me_xy,te,bare_e,ee,e1,ce);sh=spl(d.mh_xy,g.mh_xy,th,bare_h,eh,h1,ch)
 return (ex,ov,F,fs,fp,d,g,de,ee,eh,ce,ch,te,th,r,se,sh),[],F
@lru_cache(maxsize=256)
def _solve(s,T,n,pad):
 bad=[]
 for k in ('height_nm','core_radius_nm','outer_radius_nm','x_in','screening_fraction','external_field_kVcm','vbo_InN_GaN_eV','strain_c_fraction'):
  if not isinstance(getattr(s,k),(int,float)) or isinstance(getattr(s,k),bool) or not math.isfinite(getattr(s,k)):bad.append(k+' must be finite')
 if s.height_nm<=0 or s.core_radius_nm<=0:bad.append('positive geometry required')
 if s.outer_radius_nm<s.core_radius_nm:bad.append('outer radius below core radius')
 if s.strain_bound not in ('relaxed','unrelaxed') or not 0<=s.x_in<=1 or not 0<=s.screening_fraction<=1 or not 0<=s.strain_c_fraction<=1:bad.append('invalid bound or fraction')
 # z_points floor of 401 [A fix C/D]: nitride_levels._z_grid clamps
 # n_half=max(200,round(n/2)), so a base n<=401 could double to a refined
 # grid whose n_half is STILL 200 -- a bit-identical grid that would let a
 # false convergence pass silently. n>=401 guarantees 2*n+1 always lands
 # above the floor (n_half>=402), so base and refined solves are always on
 # different grids.
 if not math.isfinite(T) or T<=0 or n<401 or n%2==0 or not math.isfinite(pad) or pad<=0:bad.append('invalid numerical control')
 if bad:return _bad(bad)
 a,b,F=_once(s,T,n,pad)
 if b:return _bad(b,F)
 qz,bz,_=_once(s,T,2*n+1,pad) # z_points refined alone [fix C/E]
 if bz:return _bad(bz+['z-refinement unresolved'],F)
 qp,bp,_=_once(s,T,n,2*pad) # exterior padding refined alone [fix C/E]
 if bp:return _bad(bp+['pad-refinement unresolved'],F)
 # Rate-equivalent gate [fix C]: k=pref*exp(-Ea/kT) with pref set only by
 # the T/R/reservoir-length transverse partition (unaffected by z_points or
 # pad), so d(k)/k ~= -d(Ea)/kT to first order; gating Ea (=dE_e/dE_h, the
 # axial escape depth already computed here) to 0.02*kT bounds k_X/k_XX
 # drift under refinement to <=2% without inventing a tau_rad0_ns/
 # tau_cap_ps default inside levels() -- tau_rad0_ns has no implicit
 # default per the interface constraints, so rates() itself cannot be
 # called from this gate.  gamma_X0/XX0's refinement stability is already
 # covered by the overlap gate above (gamma_X0=overlap_sq/tau_rad0_ns).
 tol_mev=.02*KB_EV*1000.*T
 for tag,q in (('z-points',qz),('exterior padding',qp)):
  if abs(q[0]-a[0])*1000>.5:return _bad([tag+' refinement: E_X exceeds 0.5 meV'],F)
  if abs(q[1]-a[1])/a[1]>.02:return _bad([tag+' refinement: overlap exceeds 2 percent'],F)
  if abs((q[10]-q[8])-(a[10]-a[8]))*1000>tol_mev:return _bad([tag+' refinement: electron escape depth exceeds rate tolerance'],F)
  if abs((q[11]-q[9])-(a[11]-a[9]))*1000>tol_mev:return _bad([tag+' refinement: hole escape depth exceeds rate tolerance'],F)
 ex,ov,F,fs,fp,d,g,de,ee,eh,ce,ch,te,th,r,se,sh=a
 meta='[E] separable BDD axial/cylinder; [A] vacuum wall, no dielectric images or alloy localization; refinement passed (z_points, exterior_nm and dE_e/dE_h rate-equivalent tolerance gated separately)'
 meta+='; [A] dE_pair_meV is None: no pair-correlation correction to the single-particle escape depths is modeled'
 if not (math.isfinite(se) and math.isfinite(sh)):meta+='; [A] sp_split_e/h nan marks an absent J11 transverse excited state at this geometry, not zero RT-injector selectivity'
 return NanowireLevels(ex,HC/ex,True,True,ov,F,fs,fp,d.Psp_Cm2,g.Psp_Cm2,de['P_total_Cm2']-d.Psp_Cm2,0.,(te+ee)*1000,(th+eh)*1000,(ce-ee)*1000,(ch-eh)*1000,None,se*1000 if math.isfinite(se) else float('nan'),sh*1000 if math.isfinite(sh) else float('nan'),d.me_xy,d.mh_xy,ce*1000,ch*1000,te*1000,th*1000,r,s.core_radius_nm,s.outer_radius_nm,'full-core hard-wall cylinder',meta,True,(),'[V] Bernardini PRB 1997; Rinke PRB 2008; BenDaniel and Duke PR 1966')
def levels(system,T_K=300.,*,z_points=1201,exterior_nm=45.):
 if not isinstance(system,NitrideNanowireSystem):raise TypeError('system must be NitrideNanowireSystem')
 if isinstance(z_points,bool) or int(z_points)!=z_points:raise ValueError('z_points must be integer')
 return _solve(system,float(T_K),int(z_points),float(exterior_nm))
@lru_cache(maxsize=32)
def _part(mz,mxy,r,L,T):
 base=L*1e-9*math.sqrt(mz*M0*KB*T/(2*math.pi*HB**2));tot=0.
 for m in range(400):
  add=0.;zs=jn_zeros(m,300)
  for j in zs:
   x=math.exp(-(_ep(mxy,r,j)-_ep(mxy,r))/1000/(KB_EV*T));add+=(1 if m==0 else 2)*x
   if x<max(tot+add,1)*1e-10:break
  tot+=add
  if m>8 and add<max(tot,1)*1e-10:return base*tot,True
 return base*tot,False
def rates(lv,T_K,*,tau_rad0_ns,tau_cap_ps,reservoir_length_nm,channel='min',k_nr_ns=0.):
 bad=list(lv.invalid_reasons);v=(T_K,tau_rad0_ns,tau_cap_ps,reservoir_length_nm,k_nr_ns)
 if not lv.valid:bad.append('levels are invalid')
 if any(not isinstance(x,(int,float)) or not math.isfinite(x) for x in v) or min(v[:4])<=0 or k_nr_ns<0:bad.append('invalid rate request')
 if channel not in ('min','electron','hole'):bad.append('unknown channel')
 if bad:return dict(gamma_X0_ns=float('nan'),gamma_XX0_ns=float('nan'),k_X_ns=float('nan'),k_XX_ns=float('nan'),E_a_meV=float('nan'),escape_prefactor_ns=float('nan'),tau_cap_ps_used=float('nan'),reservoir_state_count_e=float('nan'),reservoir_state_count_h=float('nan'),valid=False,invalid_reasons=tuple(bad),provenance='[A] invalid rate request')
 g=binary('GaN');ne,ok1=_part(g.me_z,g.me_xy,lv.core_radius_nm,reservoir_length_nm,T_K);nh,ok2=_part(g.mh_z,g.mh_xy,lv.core_radius_nm,reservoir_length_nm,T_K)
 if not ok1 or not ok2:bad.append('transverse partition unconverged')
 if bad:return dict(gamma_X0_ns=float('nan'),gamma_XX0_ns=float('nan'),k_X_ns=float('nan'),k_XX_ns=float('nan'),E_a_meV=float('nan'),escape_prefactor_ns=float('nan'),tau_cap_ps_used=tau_cap_ps,reservoir_state_count_e=2*ne,reservoir_state_count_h=2*nh,valid=False,invalid_reasons=tuple(bad),provenance='[A] partition rejected')
 ne*=2;nh*=2;name='electron' if channel=='min' and lv.dE_e_meV<=lv.dE_h_meV else ('hole' if channel=='min' else channel);ea=lv.dE_e_meV if name=='electron' else lv.dE_h_meV;pref=1000/tau_cap_ps*((ne if name=='electron' else nh)/2);gx=lv.overlap_sq/tau_rad0_ns
 # [fix L] k_nr_ns is intrinsic occupied-dot loss, not a thermal escape
 # attempt: only the escape term doubles for XX (two independent carriers
 # each attempting thermal escape), so k_XX_ns=2*k_escape+k_nr_XX with
 # k_nr_XX=k_nr_ns (there is no separate k_nr_XX_ns parameter in this
 # frozen signature for the caller to override).
 k_escape=pref*math.exp(-ea/1000/(KB_EV*T_K));kx=k_escape+k_nr_ns;kxx=2*k_escape+k_nr_ns
 return dict(gamma_X0_ns=gx,gamma_XX0_ns=2*gx,k_X_ns=kx,k_XX_ns=kxx,E_a_meV=ea,escape_prefactor_ns=pref,tau_cap_ps_used=tau_cap_ps,reservoir_state_count_e=ne,reservoir_state_count_h=nh,valid=True,invalid_reasons=(),provenance='[DR] nondegenerate cylindrical reservoir, spin/angular states; [A] XX=2X thermal-escape doubling, occupied-dot k_nr_ns NOT doubled (k_XX_ns=2*k_escape+k_nr_ns)')
