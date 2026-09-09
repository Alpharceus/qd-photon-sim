"""Compact wurtzite InGaN/GaN p-i-n transport approximation [E/A].

Reuses transport.py's material-independent kernels (QFL suppression, window
integrals, result containers) and follows its own conventions where they
overlap: the Depletion unit contract (kV/cm field, mV active-layer drop),
the b_e definition (background photons per X photon inside the window,
module docstring of transport.py), and the r_captured + r_matrix == supply
carrier-conservation identity.  See NITRIDE_ASSUMPTIONS_E below for every
[A]/[E] default, mirroring transport.py's ASSUMPTIONS_E table.
"""
from dataclasses import dataclass
import math
from .nitride_materials import KB_EV, EPS0_SI, binary, ingaN, bandgap, band_edges
from .transport import (Q_SI, effective_dos, xi_window, qfl_suppression,
                        Depletion, Leakage, DotLoading, Background, InjectionResult)

def _dos(m,T): return effective_dos((m.me_xy*m.me_xy*m.me_z)**(1/3),T)
def _nv(m,T): return effective_dos((m.mh_xy*m.mh_xy*m.mh_z)**(1/3),T)

# [E] order-unity thermionic-leakage prefactors: no cladding minority
# mobility/diffusion-length data for nitrides is in the digest (unlike the
# InP tier's mu_n_cm2/mu_p_cm2 from transport.ASSUMPTIONS_E), so the barrier
# exponential is used directly rather than inventing a diffusion-length
# ratio.  Doping (N_A, N_D) dependence is likewise omitted [A]: this is a
# barrier-limited, not diffusion-limited, compact approximation.
_LEAK_PREF_E = 1.0
_LEAK_PREF_H = 1.0

NITRIDE_ASSUMPTIONS_E = {
    "N_A / N_D": "1e17 / 1e18 cm^-3 -- Zhang et al., APL 108, 153102 (2016), p.2 [V source, E planar transfer]",
    "d_i_nm / d_active_nm": "27 (two 12 nm GaN spacers + 3 nm active) / 3 -- Zhang 2016 [V source]",
    "wl_thickness_nm": "3 -- SRH recombination-volume thickness of the InGaN reservoir, distinct from d_active_nm (field drop layer); same Zhang 2016 QW thickness [E transfer]",
    "R_s_ohm": "100 -- small planar nitride mesa incl. contacts; no measured value in the digest [E], 2x transport.py's InP default (nitride ohmic contacts are typically more resistive)",
    "f_Rs_local": "0.3 -- same rationale as transport.Diode.f_Rs_local (ASSUMPTIONS_E, transport.py): fraction of I^2 R_s dissipated at the mesa itself [E]",
    "n_ideality": "2.0 -- SRH-dominated i-region recombination at low bias [E], same as transport.py",
    "tau_SRH_ns": "1.0 -- SRH lifetime, quality dependent [E], same class as transport.py",
    "tau_cap0_ps / tau_matrix_ns": "10 / 1.0 -- dot capture time at 1e10 cm^-2 and matrix/WL lifetime [E]; no nitride-specific capture-time measurement in the digest",
    "eta_rad_matrix": "0.1 -- radiative efficiency of matrix/WL recombination [E], same default as transport.py",
    "eps_r (diode i-region)": "10.28 -- Bernardini & Fiorentini, pss(b) 216, 391 (1999), Sec. III GaN static dielectric constant [V], used for the (GaN-dominated) 27 nm i-region rather than a doping- or QW-weighted average [A]",
    "_LEAK_PREF_E / _LEAK_PREF_H": "1.0 / 1.0 -- order-unity thermionic-leakage prefactors [E]; barrier energies come from nitride_materials.band_edges (T, x_in dependent), not literal constants",
    "doping dependence of eta_inj": "omitted [A] -- barrier-limited compact model, no cladding diffusion data",
}

@dataclass
class NitrideDiode:
    N_A: float=1e17; N_D: float=1e18; d_i_nm: float=27.; d_active_nm: float=3.; area_um2: float=.785
    R_s_ohm: float=100.; n_ideality: float=2.; tau_SRH_ns: float=1.; tau_cap0_ps: float=10.; tau_matrix_ns: float=1.
    eps_r: float=10.28; T: float=300.; x_in: float=.15; wl_thickness_nm: float=3.
    f_Rs_local: float=.3  # [E] share of I^2 R_s dissipated at the mesa (see NITRIDE_ASSUMPTIONS_E)
    def __post_init__(self):
        vals=(self.N_A,self.N_D,self.d_i_nm,self.d_active_nm,self.area_um2,self.tau_SRH_ns,
              self.tau_cap0_ps,self.tau_matrix_ns,self.R_s_ohm,self.eps_r,self.n_ideality,
              self.f_Rs_local,self.wl_thickness_nm)
        if not all(math.isfinite(v) and v>0 for v in vals) or not 0<=self.x_in<=1:
            raise ValueError("invalid nitride diode geometry, doping, composition, or transport parameter")
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
    def _j0(self,T):
        # SRH [E], Sze & Ng 3rd ed.  Two-region recombination volume, as
        # transport.py's ln_j0 splits d_active (active) + (d_i - d_active)
        # (barrier): the InGaN reservoir (wl_thickness_nm) has the small
        # nitride ni; the surrounding GaN spacers (d_i_nm - wl_thickness_nm)
        # have GaN's much larger gap and far smaller ni.  Using the whole
        # 27 nm i-region at the InGaN ni (the previous behaviour) overstated
        # the SRH dark current by treating the GaN spacers as if they had
        # the QW's narrow gap.
        ni_a=self._ni(T)
        gan=binary("GaN"); ni_b=math.sqrt(_dos(gan,T)*_nv(gan,T))*math.exp(-bandgap(gan,T)/(2*KB_EV*T))
        d_a=self.wl_thickness_nm*1e-7; d_b=max(self.d_i_nm-self.wl_thickness_nm,0.)*1e-7  # cm
        return Q_SI*(ni_a*d_a+ni_b*d_b)/(2*self.tau_SRH_ns*1e-9)
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
        # Abrupt p-i-n quadratic, same formula and unit contract as
        # transport.Diode.depletion (transport.py ~line 527): the whole
        # (V_bi - V_j) drop is NOT assumed to fall across the i-region alone
        # -- x_p, x_n solve the standard depletion quadratic in the doped
        # sides, using this diode's OWN N_A/N_D (previously the doped-side
        # widths were hardcoded to zero, over-stating F_kVcm and
        # C_dep_pF).  dV_active_mV = F [V/m] * d_active_nm [nm->m] * [V->mV]
        # (previously missing the 1e-9 * 1e3 = 1e-6 unit factor entirely).
        T=self._T(T); vbi=self.vbi(T)
        eps=self.eps_r*EPS0_SI; d_i=self.d_i_nm*1e-9
        NA,ND=self.N_A*1e6,self.N_D*1e6  # cm^-3 -> m^-3
        c=vbi-V_j
        if c<=0.:
            W=d_i
            return Depletion(V_j,vbi,0.,0.,0.,0.,W*1e9,eps*self.area_um2*1e-12/W*1e12,True)
        a=(Q_SI*NA/(2.*eps))*(1.+NA/ND)
        b=Q_SI*NA*d_i/eps
        x_p=(-b+math.sqrt(b*b+4.*a*c))/(2.*a)
        x_n=x_p*NA/ND
        F=Q_SI*NA*x_p/eps  # V/m
        W=x_p+d_i+x_n
        C=eps*self.area_um2*1e-12/W
        return Depletion(V_j,vbi,F*1e-5,F*self.d_active_nm*1e-9*1e3,x_p*1e9,x_n*1e9,W*1e9,C*1e12,False)
    def eta_inj(self,T=None):
        # Band-aware, compact thermionic leakage [A/E]: electron barrier is
        # the InGaN-active-to-GaN-cladding conduction offset, hole barrier
        # is the corresponding valence offset, both from
        # nitride_materials.band_edges (T- and x_in-dependent, replacing the
        # previous literal 0.2/0.1 eV constants that never varied with
        # composition or temperature).  Doping dependence is not modeled
        # (NITRIDE_ASSUMPTIONS_E).
        T=self._T(T); kT=KB_EV*T
        gan=binary("GaN")
        e_a=band_edges(self.active,T,substrate=gan)   # InGaN active, strained on GaN
        e_c=band_edges(gan,T,substrate=gan)            # GaN cladding, unstrained (ep=0)
        dE_c=e_c["Ec_eV"]-e_a["Ec_eV"]  # electron confinement barrier, active -> p-GaN cladding
        dE_v=e_a["Ev_eV"]-e_c["Ev_eV"]  # hole confinement barrier, active -> n-GaN cladding
        def r_leak(dE,pref):
            x=-dE/kT
            return pref*math.exp(x) if x<700 else float('inf')
        r_e=r_leak(dE_c,_LEAK_PREF_E); r_h=r_leak(dE_v,_LEAK_PREF_H)
        eta_e=1./(1.+r_e); eta_h=1./(1.+r_h); eta=1./(1.+r_e+r_h) if math.isfinite(r_e+r_h) else 0.
        return Leakage(T,float(dE_c),float(dE_v),float(r_e),float(r_h),float(eta_e),float(eta_h),float(eta),"G",{"G":(float(dE_c),float(r_e))})
    def f_qd(self,n):
        if not math.isfinite(n) or n<=0: raise ValueError("n_dot_cm2 must be positive")
        return 1/(1+self.tau_cap0_ps*1e-12*1e10/n/(self.tau_matrix_ns*1e-9))
    def dot_loading(self,I_uA,T=None,n_dot_cm2=1e10,aperture_um2=.785,tau_pulse_ns=None,f_capture=None,tau_rad_ns=1.,leakage=None,E_X_eV=None):
        """Per-dot loading in the optically selected aperture.

        The electrical supply is first apportioned by aperture_um2/area_um2,
        then shared across max(N_dots, 1).  This mirrors
        transport.Diode.dot_loading's max(N_dots, 1) fractional-occupancy
        convention, while additionally respecting this planar diode's mesa
        current geometry [DR].  Thus a fractional expected dot cannot receive
        more than the aperture's share of the mesa current.
        """
        if not math.isfinite(I_uA) or I_uA<0 or aperture_um2<=0 or tau_rad_ns<=0: raise ValueError("invalid current or geometry")
        T=self._T(T); lk=leakage or self.eta_inj(T); fq=self.f_qd(n_dot_cm2) if f_capture is None else f_capture; nd=n_dot_cm2*aperture_um2*1e-8; ne=max(nd,1.)
        vj=self.vj_of_j(I_uA*1e-6/self.area_cm2,T) if E_X_eV is not None else None; qfl=1. if E_X_eV is None else qfl_suppression(E_X_eV,vj,KB_EV*T)
        aperture_fraction=min(aperture_um2/self.area_um2,1.)
        supply=lk.eta_inj*I_uA*1e-6/Q_SI*aperture_fraction; rd=fq*supply*qfl/ne; mu=None if tau_pulse_ns is None else rd*tau_pulse_ns*1e-9
        return DotLoading(I_uA,T,lk.eta_inj,fq,nd,rd,mu,rd*tau_rad_ns*1e-9,rd*tau_rad_ns*1e-9>1,fq*lk.eta_inj*qfl*aperture_fraction/ne,qfl,vj,fq*supply*qfl,(1-fq*qfl)*supply)
    def current_for_mu(self,mu_target,T,n_dot_cm2,aperture_um2,tau_pulse_ns,E_X_eV=None):
        if mu_target<0 or tau_pulse_ns<=0: raise ValueError("invalid loading target")
        if E_X_eV is not None: raise ValueError("sub-turn-on inverse has no closed form")
        return mu_target/self.dot_loading(1.,T,n_dot_cm2,aperture_um2,tau_pulse_ns).mu
    def b_e(self,I_uA,T=None,w_meV=1.,dE_WL_meV=100.,n_dot_cm2=1e10,aperture_um2=.785,S_dot=1.,eta_rad_matrix=.1,E_urbach_meV=None,tau_rad_ns=1.,f_capture=None,leakage=None,E_X_eV=None):
        # Same f_qfl_bg convention as transport.Diode.b_e (transport.py
        # ~line 795, council review 2026-09-05 item 1): the background
        # reservoir (WL/matrix, at E_X + dE_WL) needs its OWN sub-turn-on
        # suppression, not a hardcoded 1.0, while the denominator (r_dot)
        # already carries f_qfl -- the one-sided asymmetry previously
        # reintroduced the bug transport.py fixes.
        T=self._T(T); ld=self.dot_loading(I_uA,T,n_dot_cm2,aperture_um2,None,f_capture,tau_rad_ns,leakage,E_X_eV); eu=max(KB_EV*T*1e3,E_urbach_meV or 0.); xi=float(xi_window(w_meV,dE_WL_meV,eu)); s=S_dot(T) if callable(S_dot) else S_dot
        f_qfl_bg=1.
        if E_X_eV is not None: f_qfl_bg=qfl_suppression(E_X_eV+dE_WL_meV*1e-3,ld.V_j,KB_EV*T)
        bg=ld.r_matrix*eta_rad_matrix*xi*f_qfl_bg; rx=min(ld.r_dot,1/(tau_rad_ns*1e-9))*s
        return Background(bg/rx if rx else float('inf'),T,I_uA,xi,eu,ld.eta_inj,ld.f_QD,bg,rx,ld.r_dot,s,ld.saturated,f_qfl_bg,"[A] reservoir background")
    def junction_power(self,I_uA,T,eta_total=.01,h_nu_eV=None):
        # P = I V_j + f_Rs_local I^2 R_s - (I/q) eta_total h_nu, matching
        # transport.Diode.junction_power exactly (previously charged 100% of
        # I^2 R_s to the junction and ignored an h_nu_eV override, so the
        # same R_s produced a different heat load than the InP tier's model
        # for no stated reason; no clamp here either, matching transport.py).
        I=I_uA*1e-6; _,vj=self.v_of_i(I,T); hnu=bandgap(self.active,T) if h_nu_eV is None else h_nu_eV
        return I*vj+self.f_Rs_local*I*I*self.R_s_ohm-I*eta_total*hnu

def planar_pin(**overrides): return NitrideDiode(**overrides)
def evaluate_injection(diode,I_uA,T,n_dot_cm2,aperture_um2,tau_pulse_ns,w_meV,dE_WL_meV,S_dot=1.,tau_rad_ns=1.,eta_rad_matrix=.1,E_urbach_meV=None,eta_total=.01,E_X_eV=None):
    I=I_uA*1e-6; va,vj=diode.v_of_i(I,T); dep=diode.depletion(vj,T); lk=diode.eta_inj(T); ld=diode.dot_loading(I_uA,T,n_dot_cm2,aperture_um2,tau_pulse_ns,None,tau_rad_ns,lk,E_X_eV); bg=diode.b_e(I_uA,T,w_meV,dE_WL_meV,n_dot_cm2,aperture_um2,S_dot,eta_rad_matrix,E_urbach_meV,tau_rad_ns,None,lk,E_X_eV); return InjectionResult(I_uA,T,vj,va,dep.V_bi,diode.vj_of_j(10.,T),I/diode.area_cm2,dep,lk,ld,bg,diode.junction_power(I_uA,T,eta_total),ld.mu,ld.eta_capture_dot,bg.b_e,ld.f_qfl)
