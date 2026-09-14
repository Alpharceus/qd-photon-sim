"""Full-core c-axis InGaN disc-in-GaN-wire envelope model.

Inputs use Bernardini et al., PRB 56, R10024 (1997) polarization and Rinke
et al., PRB 77, 075202 (2008) masses [V].  A hard vacuum cylinder and
separable single-band envelope are [A/E], not a fit to Deshpande data.
"""
from dataclasses import dataclass
from functools import lru_cache
import math
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.special import jn_zeros,jv
from .nitride_materials import EPS0_SI,band_edges,binary,ingaN
HB=1.054571817e-34; M0=9.1093837015e-31; KB=1.380649e-23; KBE=8.617333262e-5; HC=1239.841984; E2=1.43996448
J01=2.4048255577; J11=3.8317059702 # [V] Abramowitz & Stegun 1964, Table 9.5
@dataclass(frozen=True)
class NitrideNanowireSystem:
 height_nm:float=2.; core_radius_nm:float=12.5; outer_radius_nm:float=12.5; x_in:float=.4; strain_bound:str='relaxed'; screening_fraction:float=0.; external_field_kVcm:float=0.; vbo_InN_GaN_eV:float=1.15; strain_c_fraction:float=.7
@dataclass(frozen=True)
class NanowireLevels:
 E_X_eV:float; lambda_nm:float; electron_bound:bool; hole_bound:bool; overlap_sq:float; field_kVcm:float; F_sp_kVcm:float; F_pz_kVcm:float; P_sp_dot_Cm2:float; P_sp_GaN_Cm2:float; P_pz_dot_Cm2:float; P_pz_GaN_Cm2:float; E_e_meV:float; E_h_meV:float; dE_e_meV:float; dE_h_meV:float; dE_pair_meV:float|None; sp_split_e_meV:float; sp_split_h_meV:float; m_e_matrix_xy:float; m_h_matrix_xy:float; axial_threshold_e_meV:float; axial_threshold_h_meV:float; transverse_e_meV:float; transverse_h_meV:float; rms_radius_nm:float; geometry:str; approximation_metadata:str; valid:bool; invalid_reasons:tuple; provenance:str
def _ep(m,r,j=J01): return HB*HB*j*j/(2*m*M0*(r*1e-9)**2)/1.602176634e-22
def _rms(r):
 u=np.linspace(0,1,2001);f=u*jv(0,J01*u)**2;return r*math.sqrt(np.trapezoid(u*u*f,u)/np.trapezoid(f,u))
def _bad(b,F=0.):
 n=float('nan');return NanowireLevels(E_X_eV=n,lambda_nm=n,electron_bound=False,hole_bound=False,overlap_sq=0,field_kVcm=F,F_sp_kVcm=n,F_pz_kVcm=n,P_sp_dot_Cm2=n,P_sp_GaN_Cm2=n,P_pz_dot_Cm2=n,P_pz_GaN_Cm2=n,E_e_meV=n,E_h_meV=n,dE_e_meV=n,dE_h_meV=n,dE_pair_meV=None,sp_split_e_meV=n,sp_split_h_meV=n,m_e_matrix_xy=n,m_h_matrix_xy=n,axial_threshold_e_meV=n,axial_threshold_h_meV=n,transverse_e_meV=n,transverse_h_meV=n,rms_radius_nm=n,geometry='full-core cylinder',approximation_metadata='[A] rejected input',valid=False,invalid_reasons=tuple(b),provenance='[A] invalid state')
def _z(h,e,V,md,mb,F,q,n):
 z=np.linspace(-h/2-e,h/2+e,n);dz=z[1]-z[0];inside=np.abs(z)<=h/2;m=np.where(inside,md,mb);zc=np.clip(z,-h/2,h/2);p=np.where(inside,0,V)+q*F*zc*1e-4
 inv=2/(m[:-1]+m[1:]);c=.0380998212/dz**2;o=-c*inv;d=np.empty(n);d[0]=c/m[0]+p[0];d[-1]=c/m[-1]+p[-1];d[1:-1]=c*(inv[:-1]+inv[1:])+p[1:-1]
 w,v=eigh_tridiagonal(d,o,select='i',select_range=(0,3));return w,v,z,min(V-q*F*h/2e4,V+q*F*h/2e4)
@lru_cache(maxsize=256)
def _solve(s,T,n,ext):
 b=[]
 for k in ('height_nm','core_radius_nm','outer_radius_nm','x_in','screening_fraction','external_field_kVcm','vbo_InN_GaN_eV','strain_c_fraction'):
  if not math.isfinite(getattr(s,k)):b.append(k+' must be finite')
 if s.height_nm<=0 or s.core_radius_nm<=0:b.append('positive geometry required')
 if s.outer_radius_nm<s.core_radius_nm:b.append('outer radius below core radius')
 if not 0<=s.x_in<=1 or not 0<=s.screening_fraction<=1 or not 0<=s.strain_c_fraction<=1:b.append('fraction out of range')
 if s.strain_bound not in ('relaxed','unrelaxed'):b.append('strain_bound must be relaxed or unrelaxed')
 if not math.isfinite(T) or T<=0 or n<101 or n%2==0 or ext<=0:b.append('invalid numerical control')
 if b:return _bad(b)
 d=ingaN(s.x_in);g=binary('GaN');sf=0 if s.strain_bound=='relaxed' else 1
 de=band_edges(d,T,substrate=g,strain_fraction=sf,vbo_InN_GaN_eV=s.vbo_InN_GaN_eV,strain_c_fraction=s.strain_c_fraction);be=band_edges(g,T,substrate=g)
 ppz=de['P_total_Cm2']-d.Psp_Cm2;fs=(g.Psp_Cm2-d.Psp_Cm2)/(EPS0_SI*d.eps_r)*1e-5;fp=-ppz/(EPS0_SI*d.eps_r)*1e-5;F=(1-s.screening_fraction)*(fs+fp)+s.external_field_kVcm
 ted=_ep(d.me_xy,s.core_radius_nm)/1000;teg=_ep(g.me_xy,s.core_radius_nm)/1000;thd=_ep(d.mh_xy,s.core_radius_nm)/1000;thg=_ep(g.mh_xy,s.core_radius_nm)/1000;ve=be['Ec_eV']-de['Ec_eV']+teg-ted;vh=de['Ev_eV']-be['Ev_eV']+thg-thd
 if ve<=0 or vh<=0:return _bad(['nonpositive ground-channel axial barrier'],F)
 we,pe,z,ce=_z(s.height_nm,ext,ve,d.me_z,g.me_z,F,-1,n);wh,ph,_,ch=_z(s.height_nm,ext,vh,d.mh_z,g.mh_z,F,1,n);we2,_,_,_=_z(s.height_nm,2*ext,ve,d.me_z,g.me_z,F,-1,n);wh2,_,_,_=_z(s.height_nm,2*ext,vh,d.mh_z,g.mh_z,F,1,n)
 if not we[0]<ce or not wh[0]<ch:return _bad(['axial state unbound'],F)
 pe=pe[:,0]/math.sqrt(np.trapezoid(pe[:,0]**2,z));ph=ph[:,0]/math.sqrt(np.trapezoid(ph[:,0]**2,z));ov=max(0,min(1,float(np.trapezoid(pe*ph,z)**2)));sep=abs(np.trapezoid(z*pe*pe,z)-np.trapezoid(z*ph*ph,z));r=_rms(s.core_radius_nm);bind=math.sqrt(math.pi)*E2/((d.eps_r+g.eps_r)/2*math.sqrt(2*r*r+2*sep*sep));ex=de['Ec_eV']-de['Ev_eV']+ted+thd+we[0]+wh[0]-bind
 return NanowireLevels(ex,HC/ex,True,True,ov,F,fs,fp,d.Psp_Cm2,g.Psp_Cm2,ppz,0.,(ted+we[0])*1000,(thd+wh[0])*1000,(ce-we[0])*1000,(ch-wh[0])*1000,None,(we[1]-we[0])*1000 if we[1]<ce else float('nan'),(wh[1]-wh[0])*1000 if wh[1]<ch else float('nan'),d.me_xy,d.mh_xy,ce*1000,ch*1000,ted*1000,thd*1000,r,'full-core hard-wall cylinder','[E] separable finite-volume axial well; [A] vacuum sidewall; dielectric image, band bending and alloy localization unquantified',True,(),'[V] Bernardini 1997; Rinke 2008; BenDaniel & Duke PR 152, 683 (1966)')
def levels(system,T_K=300.,*,z_points=1201,exterior_nm=45.):
 if not isinstance(system,NitrideNanowireSystem):raise TypeError('system must be NitrideNanowireSystem')
 if int(z_points)!=z_points:raise ValueError('z_points must be integer')
 return _solve(system,float(T_K),int(z_points),float(exterior_nm))
def _N(mz,mxy,r,L,T):
 base=L*1e-9*math.sqrt(mz*M0*KB*T/(2*math.pi*HB**2));tot=0
 for m in range(30):
  q=jn_zeros(m,20);add=sum((1 if m==0 else 2)*math.exp(-(_ep(mxy,r,x)-_ep(mxy,r))/1000/(KBE*T)) for x in q);tot+=add
  if m>4 and add<tot*1e-10:break
 return base*tot
def rates(lv,T_K,*,tau_rad0_ns,tau_cap_ps,reservoir_length_nm,channel='min',k_nr_ns=0.):
 if not lv.valid or min(T_K,tau_rad0_ns,tau_cap_ps,reservoir_length_nm)<=0 or k_nr_ns<0 or channel not in ('min','electron','hole'):return dict(gamma_X0_ns=float('nan'),gamma_XX0_ns=float('nan'),k_X_ns=float('nan'),k_XX_ns=float('nan'),E_a_meV=float('nan'),escape_prefactor_ns=float('nan'),tau_cap_ps_used=float('nan'),reservoir_state_count_e=float('nan'),reservoir_state_count_h=float('nan'),valid=False,invalid_reasons=('invalid rate request',))
 g=binary('GaN');r=lv.rms_radius_nm/.468;ne=2*_N(g.me_z,g.me_xy,r,reservoir_length_nm,T_K);nh=2*_N(g.mh_z,g.mh_xy,r,reservoir_length_nm,T_K);name='electron' if channel=='min' and lv.dE_e_meV<=lv.dE_h_meV else ('hole' if channel=='min' else channel);ea=lv.dE_e_meV if name=='electron' else lv.dE_h_meV;pref=1000/tau_cap_ps*((ne if name=='electron' else nh)/2);gx=lv.overlap_sq/tau_rad0_ns;k=pref*math.exp(-ea/1000/(KBE*T_K))+k_nr_ns
 return dict(gamma_X0_ns=gx,gamma_XX0_ns=2*gx,k_X_ns=k,k_XX_ns=2*k,E_a_meV=ea,escape_prefactor_ns=pref,tau_cap_ps_used=tau_cap_ps,reservoir_state_count_e=ne,reservoir_state_count_h=nh,valid=True,invalid_reasons=(),provenance='[DR] nondegenerate cylindrical reservoir; [A] XX=2X')
