"""Compact wurtzite InGaN/GaN p-i-n transport approximation [E/A]."""
from dataclasses import dataclass
import math
from scipy.optimize import brentq
from .nitride_materials import KB_EV, EPS0_SI, binary, ingaN, bandgap
from .transport import (Q_SI, NC300_CM3, effective_dos, xi_window, qfl_suppression,
                        Depletion, Leakage, DotLoading, Background, InjectionResult)

def _dos(m,T): return effective_dos((m.me_xy*m.me_xy*m.me_z)**(1/3),T)
def _nv(m,T): return effective_dos((m.mh_xy*m.mh_xy*m.mh_z)**(1/3),T)

@dataclass
class NitrideDiode:
    N_A: float=1e17; N_D: float=1e18; d_i_nm: float=27.; d_active_nm: float=3.; area_um2: float=.785
    R_s_ohm: float=100.; n_ideality: float=2.; tau_SRH_ns: float=1.; tau_cap0_ps: float=10.; tau_matrix_ns: float=1.
    eps_r: float=9.5; T: float=300.; x_in: float=.15; wl_thickness_nm: float=3.
    def __post_init__(self):
        if min(self.N_A,self.N_D,self.d_i_nm,self.d_active_nm,self.area_um2,self.tau_SRH_ns,self.tau_cap0_ps,self.tau_matrix_ns)<=0 or not 0<=self.x_in<=1: raise ValueError("invalid nitride diode geometry, doping, or composition")
    @property
    def area_cm2(self): return self.area_um2*1e-8
    @property
    def active(self): return ingaN(self.x_in)
    def _T(self,T):
        T=self.T if T is None else float(T)
        if not math.isfinite(T) or T<=0: raise ValueError("invalid temperature")
        return T
    def _ni(self,T):
        m=self.active; return math.sqrt(_dos(m,T)*_nv(m,T))*math.exp(-bandgap(m,T)/(2*KB_EV*T))
    def vbi(self,T=None):
        T=self._T(T); m=binary("GaN"); ni=math.sqrt(_dos(m,T)*_nv(m,T))*math.exp(-bandgap(m,T)/(2*KB_EV*T)); return KB_EV*T*math.log(self.N_A*self.N_D/(ni*ni))
    def _j0(self,T): return Q_SI*self._ni(T)*self.d_i_nm*1e-7/(2*self.tau_SRH_ns*1e-9) # SRH [E], Sze & Ng 3rd ed.
    def j_of_vj(self,V_j,T=None):
        T=self._T(T); x=V_j/(self.n_ideality*KB_EV*T)
        if x>700: return float('inf')
        return self._j0(T)*math.expm1(x)
    def vj_of_j(self,J,T=None):
        T=self._T(T)
        if J<=0:return 0.
        return self.n_ideality*KB_EV*T*math.log1p(J/self._j0(T))
    def v_of_i(self,I_A,T=None):
        if not math.isfinite(I_A) or I_A<0: raise ValueError("I_A must be finite and nonnegative")
        vj=self.vj_of_j(I_A/self.area_cm2,T); return vj+I_A*self.R_s_ohm,vj
    def depletion(self,V_j,T=None):
        T=self._T(T); vbi=self.vbi(T); W=self.d_i_nm*1e-9; flat=V_j>=vbi; F=0. if flat else max(vbi-V_j,0)/W
        return Depletion(V_j,vbi,F*1e-5,F*self.d_active_nm,0.,0.,self.d_i_nm,self.eps_r*EPS0_SI*self.area_um2*1e-12/W*1e12,flat)
    def eta_inj(self,T=None):
        # finite capture/leakage is an explicit compact approximation [A], not fitted.
        return Leakage(self._T(T),.2,.1,.1,.1,1/1.1,1/1.1,1/1.2,"G",{"G":(.2,.1)})
    def f_qd(self,n):
        if not math.isfinite(n) or n<=0: raise ValueError("n_dot_cm2 must be positive")
        return 1/(1+self.tau_cap0_ps*1e-12*1e10/n/(self.tau_matrix_ns*1e-9))
    def dot_loading(self,I_uA,T=None,n_dot_cm2=1e10,aperture_um2=.785,tau_pulse_ns=None,f_capture=None,tau_rad_ns=1.,leakage=None,E_X_eV=None):
        if I_uA<0 or aperture_um2<=0 or tau_rad_ns<=0: raise ValueError("invalid current or geometry")
        T=self._T(T); lk=leakage or self.eta_inj(T); fq=self.f_qd(n_dot_cm2) if f_capture is None else f_capture; nd=n_dot_cm2*aperture_um2*1e-8; ne=max(nd,1.)
        vj=self.vj_of_j(I_uA*1e-6/self.area_cm2,T) if E_X_eV is not None else None; qfl=1. if E_X_eV is None else qfl_suppression(E_X_eV,vj,KB_EV*T)
        supply=lk.eta_inj*I_uA*1e-6/Q_SI; rd=fq*supply*qfl/ne; mu=None if tau_pulse_ns is None else rd*tau_pulse_ns*1e-9
        return DotLoading(I_uA,T,lk.eta_inj,fq,nd,rd,mu,rd*tau_rad_ns*1e-9,rd*tau_rad_ns*1e-9>1,fq*lk.eta_inj*qfl/ne,qfl,vj,fq*supply*qfl,(1-fq*qfl)*supply)
    def current_for_mu(self,mu_target,T,n_dot_cm2,aperture_um2,tau_pulse_ns,E_X_eV=None):
        if mu_target<0 or tau_pulse_ns<=0: raise ValueError("invalid loading target")
        if E_X_eV is not None: raise ValueError("sub-turn-on inverse has no closed form")
        return mu_target/self.dot_loading(1.,T,n_dot_cm2,aperture_um2,tau_pulse_ns).mu
    def b_e(self,I_uA,T=None,w_meV=1.,dE_WL_meV=100.,n_dot_cm2=1e10,aperture_um2=.785,S_dot=1.,eta_rad_matrix=.1,E_urbach_meV=None,tau_rad_ns=1.,f_capture=None,leakage=None,E_X_eV=None):
        T=self._T(T); ld=self.dot_loading(I_uA,T,n_dot_cm2,aperture_um2,None,f_capture,tau_rad_ns,leakage,E_X_eV); eu=max(KB_EV*T*1e3,E_urbach_meV or 0.); xi=float(xi_window(w_meV,dE_WL_meV,eu)); s=S_dot(T) if callable(S_dot) else S_dot; bg=ld.r_matrix*eta_rad_matrix*xi; rx=min(ld.r_dot,1/(tau_rad_ns*1e-9))*s; return Background(bg/rx if rx else float('inf'),T,I_uA,xi,eu,ld.eta_inj,ld.f_QD,bg,rx,ld.r_dot,s,ld.saturated,1.,"[A] reservoir background")
    def junction_power(self,I_uA,T,eta_total=.01):
        I=I_uA*1e-6; _,vj=self.v_of_i(I,T); return max(0.,I*vj+I*I*self.R_s_ohm-I*eta_total*(bandgap(self.active,T)))

def planar_pin(**overrides): return NitrideDiode(**overrides)
def evaluate_injection(diode,I_uA,T,n_dot_cm2,aperture_um2,tau_pulse_ns,w_meV,dE_WL_meV,S_dot=1.,tau_rad_ns=1.,eta_rad_matrix=.1,E_urbach_meV=None,eta_total=.01,E_X_eV=None):
    I=I_uA*1e-6; va,vj=diode.v_of_i(I,T); dep=diode.depletion(vj,T); lk=diode.eta_inj(T); ld=diode.dot_loading(I_uA,T,n_dot_cm2,aperture_um2,tau_pulse_ns,None,tau_rad_ns,lk,E_X_eV); bg=diode.b_e(I_uA,T,w_meV,dE_WL_meV,n_dot_cm2,aperture_um2,S_dot,eta_rad_matrix,E_urbach_meV,tau_rad_ns,None,lk,E_X_eV); return InjectionResult(I_uA,T,vj,va,dep.V_bi,diode.vj_of_j(10.,T),I/diode.area_cm2,dep,lk,ld,bg,diode.junction_power(I_uA,T,eta_total),ld.mu,ld.eta_capture_dot,bg.b_e,ld.f_qfl)
