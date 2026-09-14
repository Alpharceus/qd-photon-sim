"""C-axis InGaN disc-in-GaN-wire levels; hard-wall radial envelope [A] for
the full-core (horizontal, as-built) family, finite-barrier InGaN/GaN disk
[A] for the vertical disc-in-wire family (H6: disc_radius_nm < core_radius_nm,
a real dot laterally confined inside a much larger wire core, reusing
dot_levels.finite_disk_2d's BenDaniel-Duke radial matching)."""
from dataclasses import dataclass
from functools import lru_cache
import math
import numpy as np
from scipy.special import jn_zeros,jv,kve
from .nitride_materials import EPS0_SI,band_edges,binary,ingaN,KB_EV
from .nitride_levels import _z_state,_coulomb_binding_eV
from .dot_levels import finite_disk_2d,HB2_2M0
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
 disc_radius_nm:object=None # [A] H6: None (default) = full-core horizontal family (disc fills the core, hard wall at core_radius_nm); a finite value < core_radius_nm selects the vertical disc-in-wire family (finite InGaN/GaN radial barrier, core radius unchanged for the GaN reservoir)
 x_in:float=.4 # [V] Deshpande et al., APL 2014 abstract
 strain_bound:str='relaxed' # [A] scenario label
 screening_fraction:float=0. # [A] independent polarization screening
 external_field_kVcm:float=0. # [A] independent junction field
 vbo_InN_GaN_eV:float=1.15 # [V] Tsai and Bayram, ACS Omega 2020 Table 2
 strain_c_fraction:float=.7 # [A] deformation-potential partition
@dataclass(frozen=True)
class NanowireLevels:
 E_X_eV:float;lambda_nm:float;electron_bound:bool;hole_bound:bool;overlap_sq:float;field_kVcm:float;F_sp_kVcm:float;F_pz_kVcm:float;P_sp_dot_Cm2:float;P_sp_GaN_Cm2:float;P_pz_dot_Cm2:float;P_pz_GaN_Cm2:float;E_e_meV:float;E_h_meV:float;dE_e_meV:float;dE_h_meV:float;dE_pair_meV:object;sp_split_e_meV:float;sp_split_h_meV:float;m_e_matrix_xy:float;m_h_matrix_xy:float;axial_threshold_e_meV:float;axial_threshold_h_meV:float;transverse_e_meV:float;transverse_h_meV:float;rms_radius_nm:float;core_radius_nm:float;outer_radius_nm:float;geometry:str;approximation_metadata:str;valid:bool;invalid_reasons:tuple;provenance:str;sidewall_overlap:float
def _ep(m,r,j=J01): return HB*HB*j*j/(2*m*M0*(r*1e-9)**2)/1.602176634e-22 # [DR] cylinder Bessel mode
@lru_cache(maxsize=1)
def _rf():
 u=np.linspace(0,1,4001);w=u*jv(0,J01*u)**2;return math.sqrt(np.trapezoid(u*u*w,u)/np.trapezoid(w,u))
def _hw_shell_overlap(R_nm,shell_nm=2.0):
 # H6: hard-wall (full-core) ground-mode radial density fraction found in
 # the outer shell_nm of the hard wall at R_nm [A geometric diagnostic]:
 # mass-independent since the J0 hard-wall mode SHAPE does not depend on
 # mass (only its energy scale does), so a single call covers both carriers.
 if R_nm<=shell_nm:return 1.0
 u=np.linspace(0.,1.,4001);w=u*jv(0,J01*u)**2;lo=1.-shell_nm/R_nm;mask=u>=lo
 return float(np.trapezoid(w[mask],u[mask])/np.trapezoid(w,u))
def _disc_shell_overlap(V_meV,E0_meV,m_in,m_out,R_disc_nm,core_radius_nm,shell_nm=2.0):
 # H6: evanescent overlap of the disc's finite-barrier ground state with a
 # shell_nm-wide sidewall layer AT core_radius_nm [A], i.e. the "sidewall
 # access weight" for a real dot sitting well inside a much larger wire.
 # Reuses the SAME analytic exterior tail dot_levels.finite_disk_2d builds
 # for its own rms radius (infinite-barrier exterior, valid while
 # core_radius_nm >> R_disc_nm); the integral is truncated (not
 # renormalized past) core_radius_nm because the true wire has a hard wall
 # there, not more GaN.
 if core_radius_nm<=R_disc_nm or not math.isfinite(V_meV) or not math.isfinite(E0_meV):return float('nan')
 q0=math.sqrt(max(m_out*(V_meV-E0_meV),0.)/HB2_2M0)
 if q0<=0.:return float('nan')
 k0=math.sqrt(max(m_in*E0_meV,0.)/HB2_2M0)
 r_in=np.linspace(0.,R_disc_nm,4001);r_out=np.linspace(R_disc_nm,core_radius_nm,4001)
 psi_in=jv(0,k0*r_in);denom=kve(0,q0*R_disc_nm)
 if denom==0.:return float('nan')
 psi_out=jv(0,k0*R_disc_nm)/denom*kve(0,q0*r_out)*np.exp(-q0*(r_out-R_disc_nm))
 norm=np.trapezoid(r_in*psi_in**2,r_in)+np.trapezoid(r_out*psi_out**2,r_out)
 if norm<=0.:return float('nan')
 lo=max(R_disc_nm,core_radius_nm-shell_nm);mask=r_out>=lo
 if mask.sum()<2:return 0.0
 return float(np.trapezoid(r_out[mask]*psi_out[mask]**2,r_out[mask])/norm)
def _bad(b,F=0.):
 n=float('nan');return NanowireLevels(n,n,False,False,0,F,n,n,n,n,n,n,n,n,n,n,None,n,n,n,n,n,n,n,n,n,n,n,'full-core hard-wall cylinder','[A] rejected state',False,tuple(b),'[A] invalid state',n)
def _once(s,T,n,pad):
 d,g=ingaN(s.x_in),binary('GaN');sf=0 if s.strain_bound=='relaxed' else 1 # [DR] relaxation removes strain shifts and P_pz
 de=band_edges(d,T,substrate=g,strain_fraction=sf,vbo_InN_GaN_eV=s.vbo_InN_GaN_eV,strain_c_fraction=s.strain_c_fraction);ge=band_edges(g,T,substrate=g)
 pp=de['P_total_Cm2']-d.Psp_Cm2;fs=(g.Psp_Cm2-d.Psp_Cm2)/(EPS0_SI*d.eps_r)*1e-5;fp=-pp/(EPS0_SI*d.eps_r)*1e-5;F=(1-s.screening_fraction)*(fs+fp)+s.external_field_kVcm # [V/DR] Bernardini 1997; [A] external unscreened
 bare_e=ge['Ec_eV']-de['Ec_eV'];bare_h=de['Ev_eV']-ge['Ev_eV'] # [DR] bare band offset, no transverse zero-point term
 disc=s.disc_radius_nm is not None and s.disc_radius_nm<s.core_radius_nm # H6: vertical disc-in-wire family
 if not disc:
  te,Te=_ep(d.me_xy,s.core_radius_nm)/1000,_ep(g.me_xy,s.core_radius_nm)/1000
  th,Th=_ep(d.mh_xy,s.core_radius_nm)/1000,_ep(g.mh_xy,s.core_radius_nm)/1000
  ve=bare_e+Te-te;vh=bare_h+Th-th
  if ve<=0 or vh<=0:return None,['nonpositive ground-channel axial barrier'],F
  ee,e1,ce,_,z,pe=_z_state(s.height_nm,ve,d.me_z,g.me_z,F,-1,n,pad);eh,h1,ch,_,zh,ph=_z_state(s.height_nm,vh,d.mh_z,g.mh_z,F,1,n,pad)
  if not ee<ce or not eh<ch:return None,['axial state unbound'],F
  ov=max(0,min(1,float(np.trapezoid(pe*ph,z)**2)))
  if ov<1e-12:return None,['overlap_unresolved'],F
  sep=abs(np.trapezoid(z*pe*pe,z)-np.trapezoid(zh*ph*ph,zh));r=s.core_radius_nm*_rf();r_e,r_h=r,r
  ex=de['Ec_eV']-de['Ev_eV']+te+th+ee+eh-_coulomb_binding_eV(r,r,sep,d.eps_r) # [E] dot eps envelope
  de_e,de_h=ce-ee,ch-eh # [A] axial escape depth already carries the E_perp_GaN-E_perp_InGaN correction folded into ve/vh above
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
  sw=_hw_shell_overlap(s.core_radius_nm);geom='full-core hard-wall cylinder'
 else:
  R=s.disc_radius_nm
  if bare_e<=0 or bare_h<=0:return None,['nonpositive ground-channel axial barrier'],F
  ee,e1,ce,_,z,pe=_z_state(s.height_nm,bare_e,d.me_z,g.me_z,F,-1,n,pad);eh,h1,ch,_,zh,ph=_z_state(s.height_nm,bare_h,d.mh_z,g.mh_z,F,1,n,pad)
  if not ee<ce or not eh<ch:return None,['axial state unbound'],F
  ov=max(0,min(1,float(np.trapezoid(pe*ph,z)**2)))
  if ov<1e-12:return None,['overlap_unresolved'],F
  sep=abs(np.trapezoid(z*pe*pe,z)-np.trapezoid(zh*ph*ph,zh))
  # [E] adiabatic decoupling, same order as nitride_levels._levels_cached /
  # dot_levels' own dot solve: axial solved first against the BARE offset,
  # then the disc's finite-barrier RADIAL well sees the REMAINING escape
  # depth (V - E_z), not the bare offset again and not the horizontal
  # branch's Te-te axial correction (which assumed InGaN and GaN share the
  # SAME hard-wall radius everywhere -- no longer true once the disc is
  # laterally smaller than the core).
  roff_e,roff_h=ce-ee,ch-eh
  if roff_e<=0 or roff_h<=0:return None,['nonpositive disc-radial headroom'],F
  de_disk=finite_disk_2d(roff_e*1000.,R,d.me_xy,g.me_xy);dh_disk=finite_disk_2d(roff_h*1000.,R,d.mh_xy,g.mh_xy)
  if not de_disk.bound or not dh_disk.bound:return None,['disc radial state unbound'],F
  te,th=de_disk.E0_meV/1000.,dh_disk.E0_meV/1000.
  e1r=de_disk.E1_meV/1000. if de_disk.p_bound else float('nan');h1r=dh_disk.E1_meV/1000. if dh_disk.p_bound else float('nan')
  r_e=de_disk.rms_r_nm if de_disk.rms_r_nm else R;r_h=dh_disk.rms_r_nm if dh_disk.rms_r_nm else R
  ex=de['Ec_eV']-de['Ev_eV']+te+th+ee+eh-_coulomb_binding_eV(r_e,r_h,sep,d.eps_r)
  de_e,de_h=roff_e-te,roff_h-th # true net escape depth after the disc's own radial confinement is removed
  zg_e=(e1-ee) if (math.isfinite(e1) and e1<ce) else float('inf');rg_e=(e1r-te) if math.isfinite(e1r) else float('inf')
  se=min(zg_e,rg_e);se=se if (math.isfinite(se) and se>0) else float('nan')
  zg_h=(h1-eh) if (math.isfinite(h1) and h1<ch) else float('inf');rg_h=(h1r-th) if math.isfinite(h1r) else float('inf')
  sh=min(zg_h,rg_h);sh=sh if (math.isfinite(sh) and sh>0) else float('nan')
  sw_e=_disc_shell_overlap(roff_e*1000.,de_disk.E0_meV,d.me_xy,g.me_xy,R,s.core_radius_nm)
  sw_h=_disc_shell_overlap(roff_h*1000.,dh_disk.E0_meV,d.mh_xy,g.mh_xy,R,s.core_radius_nm)
  finite_sw=[x for x in (sw_e,sw_h) if math.isfinite(x)]
  sw=max(finite_sw) if finite_sw else float('nan') # [A] conservative: report the more-leaked (worse-isolated) carrier
  geom='disc-in-wire finite-barrier radial'
 return (ex,ov,F,fs,fp,d,g,de,ee,eh,ce,ch,te,th,r_e,r_h,se,sh,de_e,de_h,sw,geom),[],F
@lru_cache(maxsize=256)
def _solve(s,T,n,pad):
 bad=[]
 for k in ('height_nm','core_radius_nm','outer_radius_nm','x_in','screening_fraction','external_field_kVcm','vbo_InN_GaN_eV','strain_c_fraction'):
  if not isinstance(getattr(s,k),(int,float)) or isinstance(getattr(s,k),bool) or not math.isfinite(getattr(s,k)):bad.append(k+' must be finite')
 if s.height_nm<=0 or s.core_radius_nm<=0:bad.append('positive geometry required')
 if s.outer_radius_nm<s.core_radius_nm:bad.append('outer radius below core radius')
 if s.disc_radius_nm is not None:
  if not isinstance(s.disc_radius_nm,(int,float)) or isinstance(s.disc_radius_nm,bool) or not math.isfinite(s.disc_radius_nm) or s.disc_radius_nm<=0:
   bad.append('disc_radius_nm must be None or a positive finite number')
  elif s.disc_radius_nm>s.core_radius_nm:
   bad.append('disc_radius_nm must not exceed core_radius_nm')
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
 # true net escape depth already computed by _once for either family) to
 # 0.02*kT bounds k_X/k_XX drift under refinement to <=2% without inventing
 # a tau_rad0_ns/tau_cap_ps default inside levels() -- tau_rad0_ns has no
 # implicit default per the interface constraints, so rates() itself cannot
 # be called from this gate. gamma_X0/XX0's refinement stability is already
 # covered by the overlap gate above (gamma_X0=overlap_sq/tau_rad0_ns).
 tol_mev=.02*KB_EV*1000.*T
 for tag,q in (('z-points',qz),('exterior padding',qp)):
  if abs(q[0]-a[0])*1000>.5:return _bad([tag+' refinement: E_X exceeds 0.5 meV'],F)
  if abs(q[1]-a[1])/a[1]>.02:return _bad([tag+' refinement: overlap exceeds 2 percent'],F)
  if abs(q[18]-a[18])*1000>tol_mev:return _bad([tag+' refinement: electron escape depth exceeds rate tolerance'],F)
  if abs(q[19]-a[19])*1000>tol_mev:return _bad([tag+' refinement: hole escape depth exceeds rate tolerance'],F)
 ex,ov,F,fs,fp,d,g,de,ee,eh,ce,ch,te,th,r_e,r_h,se,sh,de_e,de_h,sw,geom=a
 if geom.startswith('disc'):
  meta='[E] separable BDD axial well + finite-barrier radial disk (H6: adiabatic decoupling, radial well sees the axial escape depth as its barrier, dot_levels.finite_disk_2d BenDaniel-Duke matching); [A] evanescent sidewall tail truncated (not renormalized) at core_radius_nm, no dielectric images or alloy localization; refinement passed (z_points, exterior_nm and dE_e/dE_h rate-equivalent tolerance gated separately)'
 else:
  meta='[E] separable BDD axial/cylinder; [A] vacuum wall, no dielectric images or alloy localization; refinement passed (z_points, exterior_nm and dE_e/dE_h rate-equivalent tolerance gated separately)'
 meta+='; [A] dE_pair_meV is None: no pair-correlation correction to the single-particle escape depths is modeled'
 if not (math.isfinite(se) and math.isfinite(sh)):meta+='; [A] sp_split_e/h nan marks an absent transverse excited state at this geometry, not zero RT-injector selectivity'
 meta+='; z_points_used=%d'%n
 return NanowireLevels(ex,HC/ex,True,True,ov,F,fs,fp,d.Psp_Cm2,g.Psp_Cm2,de['P_total_Cm2']-d.Psp_Cm2,0.,(te+ee)*1000,(th+eh)*1000,de_e*1000,de_h*1000,None,se*1000 if math.isfinite(se) else float('nan'),sh*1000 if math.isfinite(sh) else float('nan'),d.me_xy,d.mh_xy,ce*1000,ch*1000,te*1000,th*1000,r_e,s.core_radius_nm,s.outer_radius_nm,geom,meta,True,(),'[V] Bernardini PRB 1997; Rinke PRB 2008; BenDaniel and Duke PR 1966; [E] dot_levels.finite_disk_2d radial matching (disc-in-wire family)',sw)
def levels(system,T_K=300.,*,z_points=1201,exterior_nm=45.):
 if not isinstance(system,NitrideNanowireSystem):raise TypeError('system must be NitrideNanowireSystem')
 if isinstance(z_points,bool) or int(z_points)!=z_points:raise ValueError('z_points must be integer')
 n=int(z_points);T=float(T_K);pad=float(exterior_nm)
 cap=9601 # [A] resolution-policy cap: double z_points on a gate failure and re-test; give up only here
 while True:
  res=_solve(system,T,n,pad)
  if res.valid or n>=cap:return res
  n=min(cap,2*n+1) # stay odd (n%2==0 is rejected by _solve)
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
 # [A fix L3] reservoir_length_nm is a PER-SIDE GaN reservoir extent (e.g.
 # Deshpande's 15 nm on each side of the disc); the axial partition sum
 # (Interface constraints: "L=barrier_left_nm+barrier_right_nm") counts
 # BOTH flanking reservoirs, i.e. the prefactor uses L_total=2*reservoir_
 # length_nm, not reservoir_length_nm alone.
 g=binary('GaN');L_total=2.*reservoir_length_nm
 ne,ok1=_part(g.me_z,g.me_xy,lv.core_radius_nm,L_total,T_K);nh,ok2=_part(g.mh_z,g.mh_xy,lv.core_radius_nm,L_total,T_K)
 if not ok1 or not ok2:bad.append('transverse partition unconverged')
 if bad:return dict(gamma_X0_ns=float('nan'),gamma_XX0_ns=float('nan'),k_X_ns=float('nan'),k_XX_ns=float('nan'),E_a_meV=float('nan'),escape_prefactor_ns=float('nan'),tau_cap_ps_used=tau_cap_ps,reservoir_state_count_e=2*ne,reservoir_state_count_h=2*nh,valid=False,invalid_reasons=tuple(bad),provenance='[A] partition rejected')
 ne*=2;nh*=2;name='electron' if channel=='min' and lv.dE_e_meV<=lv.dE_h_meV else ('hole' if channel=='min' else channel);ea=lv.dE_e_meV if name=='electron' else lv.dE_h_meV;pref=1000/tau_cap_ps*((ne if name=='electron' else nh)/2);gx=lv.overlap_sq/tau_rad0_ns
 # [fix L] k_nr_ns is intrinsic occupied-dot loss, not a thermal escape
 # attempt: only the escape term doubles for XX (two independent carriers
 # each attempting thermal escape), so k_XX_ns=2*k_escape+k_nr_XX with
 # k_nr_XX=k_nr_ns (there is no separate k_nr_XX_ns parameter in this
 # frozen signature for the caller to override).
 k_escape=pref*math.exp(-ea/1000/(KB_EV*T_K));kx=k_escape+k_nr_ns;kxx=2*k_escape+k_nr_ns
 return dict(gamma_X0_ns=gx,gamma_XX0_ns=2*gx,k_X_ns=kx,k_XX_ns=kxx,E_a_meV=ea,escape_prefactor_ns=pref,tau_cap_ps_used=tau_cap_ps,reservoir_state_count_e=ne,reservoir_state_count_h=nh,valid=True,invalid_reasons=(),provenance='[DR] nondegenerate cylindrical reservoir, spin/angular states, L_total=2*reservoir_length_nm (both flanking GaN reservoirs, fix L3); [A] XX=2X thermal-escape doubling, occupied-dot k_nr_ns NOT doubled (k_XX_ns=2*k_escape+k_nr_ns)')
