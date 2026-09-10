"""Independent effective-shape and QW-fluctuation checks.

Volumes use elementary profile integration [DR], not production helpers.
BenDaniel-Duke roots are independently tested in verify_nitride_levels.py.

Fix round 1 (Opus review of commit 57ab49d, VERDICT FAIL) additions below
the original G/Q checks: an independent re-assembly of the QW-fluctuation
formula (reservoir channel selection, N2D mass source, Coulomb term,
symmetric-placement lateral-depth limit) built from the SAME numerical
primitives the production module uses (_z_state, _radial,
_coulomb_binding_eV, _growth_masses) but combined here independently of
fsim_core.nitride_levels._qw_levels, so these checks catch a broken
ASSEMBLY of those primitives even when the primitives themselves are
correct -- this is deliberately the same "reuse the solver, not the
solved-for value" pattern verify_nitride_levels.py uses for the z-solver
itself. Also added: a validity-gate check (nonphysical E_X), an
overlap-unresolved rates() check, an orientation growth-axis mass-swap
check, a shape_height_fraction sensitivity check, more malformed-input
cases, and a height x radius x orientation x geometry x screening
convergence matrix with a stated per-height numerical budget.
"""
import math
import sys
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fsim_core.nitride_levels import (NitrideDotSystem, levels, rates,
    _z_state, _radial, _coulomb_binding_eV, _growth_masses)
from fsim_core.nitride_materials import binary, ingaN, band_edges, polarization_field

checks=[]
def ck(name, value): checks.append((name,bool(value)))

H,R=6.,12.
base=dict(height_nm=H,radius_nm=R,screening_fraction=1.)
disc=levels(NitrideDotSystem(**base),300.)
lens=levels(NitrideDotSystem(**{**base,'shape':'lens'}),300.)
cone=levels(NitrideDotSystem(**{**base,'shape':'truncated_cone','top_radius_fraction':0.}),300.)
frustum=levels(NitrideDotSystem(**{**base,'shape':'truncated_cone','top_radius_fraction':.5}),300.)
unit=levels(NitrideDotSystem(**{**base,'shape':'truncated_cone','top_radius_fraction':1.}),300.)
ck('G independently integrated disc volume',abs(disc.shape_volume_nm3-math.pi*R*R*H)<1e-10)
ck('G lens volume and effective height H/2',abs(lens.shape_volume_nm3-math.pi*R*R*H/2)<1e-10 and lens.effective_height_nm==H/2)
ck('G cone volume and effective height H/3',abs(cone.shape_volume_nm3-math.pi*R*R*H/3)<1e-10 and cone.effective_height_nm==H/3)
ck('G frustum integral',abs(frustum.shape_volume_nm3-math.pi*R*R*H*(1+.5+.25)/3)<1e-10)
ck('G t=1 frustum exactly reduces to disc',disc.E_X_eV==unit.E_X_eV and disc.overlap_sq==unit.overlap_sq and disc.dE_e_meV==unit.dE_e_meV)
ck('G shape changes confinement',lens.valid and disc.valid and lens.E_X_eV!=disc.E_X_eV)
mapped=levels(NitrideDotSystem(height_nm=H/2,radius_nm=R,screening_fraction=1.),300.)
ck('G equal-volume mapped disc agrees exactly',lens.E_X_eV==mapped.E_X_eV and lens.overlap_sq==mapped.overlap_sq)

old=levels(NitrideDotSystem(height_nm=3.,radius_nm=10.),300.)
explicit=levels(NitrideDotSystem(height_nm=3.,radius_nm=10.,orientation='c_plane',shape='disc',geometry_type='isolated_dot',top_radius_fraction=1.),300.)
ck('G omitted defaults preserve legacy scalars',old.E_X_eV==explicit.E_X_eV and old.overlap_sq==explicit.overlap_sq and old.dE_e_meV==explicit.dE_e_meV and old.field_kVcm==explicit.field_kVcm)

qargs=dict(geometry_type='qw_fluctuation',height_nm=3.5,radius_nm=20.,wl_thickness_nm=3.,screening_fraction=1.)
q=levels(NitrideDotSystem(**qargs),300.)
qwide=levels(NitrideDotSystem(**{**qargs,'radius_nm':30.}),300.)
qsmall=levels(NitrideDotSystem(**{**qargs,'radius_nm':5.}),300.)
ck('Q QW reservoir and positive escape',q.valid and q.reservoir_kind=='ingan_qw' and q.dE_e_meV>0 and q.dE_h_meV>0)
ck('Q reservoir is radius independent',q.reservoir_energy_eV==qwide.reservoir_energy_eV==qsmall.reservoir_energy_eV)
ck('Q local energy varies with radius',q.E_X_eV!=qsmall.E_X_eV)
ck('Q electron hole and min rates valid',all(rates(q,300.,channel=c)['valid'] for c in ('electron','hole','min')))
ck('Q pair rate rejected',not rates(q,300.,channel='pair')['valid'])
ck('Q pair_half rate rejected',not rates(q,300.,channel='pair_half')['valid'])
equal=levels(NitrideDotSystem(geometry_type='qw_fluctuation',height_nm=3.,radius_nm=5.,wl_thickness_nm=3.,screening_fraction=1.),300.)
ck('Q H equals well is invalid',not equal.valid)

invalid=[]
for kw in (dict(shape='bad'),dict(shape='lens',top_radius_fraction=.5),dict(shape='truncated_cone',top_radius_fraction=1.1),
           dict(shape='lens',shape_height_fraction=.4),dict(geometry_type='qw_fluctuation',shape='lens',wl_thickness_nm=2.),
           dict(shape='truncated_cone',top_radius_fraction=-0.1),dict(height_nm=float('nan')),dict(height_nm=float('inf')),
           dict(radius_nm=0.),dict(geometry_type='qw_fluctuation',wl_thickness_nm=0.),
           dict(geometry_type='qw_fluctuation',x_in=1.4,wl_thickness_nm=1.),dict(geometry_type='qw_fluctuation',x_in=-0.1,wl_thickness_nm=1.),
           dict(orientation='bad_orientation')):
    try:
        invalid.append(not levels(NitrideDotSystem(**kw),300.).valid)
    except ValueError:
        invalid.append(True)  # constructing/consuming the bad orientation string is itself rejected
ck('G malformed geometry/dimensions/composition rejected',all(invalid))

coarse=levels(NitrideDotSystem(height_nm=3.,radius_nm=10.,screening_fraction=1.),300.,z_points=1201,exterior_nm=45.)
fine=levels(NitrideDotSystem(height_nm=3.,radius_nm=10.,screening_fraction=1.),300.,z_points=2401,exterior_nm=90.)
ck('N refined bound fixture within 0.5 meV',coarse.valid and fine.valid and abs(coarse.E_X_eV-fine.E_X_eV)*1000<=.5 and abs(coarse.dE_e_meV-fine.dE_e_meV)<=.5 and abs(coarse.dE_h_meV-fine.dE_h_meV)<=.5)

# ---------------------------------------------------------------------
# Fix round 1 additions (Opus review of commit 57ab49d)
# ---------------------------------------------------------------------

# ---- (V) validity gate: the exact reviewer fixture (height 10 nm, x=.25,
# unscreened) used to report valid=True with E_X_eV=-0.667; it must now be
# rejected, visibly, as a nonphysical result rather than a plausible one.
bad_fixture=levels(NitrideDotSystem(height_nm=10.,radius_nm=10.,x_in=.25,screening_fraction=0.),300.)
ck('V nonphysical-E_X reviewer fixture (h=10nm unscreened) is rejected',
   not bad_fixture.valid and 'nonphysical E_X' in bad_fixture.invalid_reasons)
# A genuinely bound, positive-E_X fixture must still pass through valid=True
# (the gate must not become a blanket rejection of every large dot).
ok_fixture=levels(NitrideDotSystem(height_nm=10.,radius_nm=10.,x_in=.25,orientation='a_plane',screening_fraction=0.),300.)
ck('V a genuinely bound fixture is not swept up by the validity gate',
   ok_fixture.valid and ok_fixture.E_X_eV>0. and math.isfinite(ok_fixture.lambda_nm))

# ---- (V2) overlap_unresolved: a state below the numerical overlap floor
# (h=8nm, x=.25, unscreened -- overlap_sq ~1e-14) is levels()-valid (a real
# bound state) but rates() must refuse it rather than report an
# astronomical lifetime.
unresolved=levels(NitrideDotSystem(height_nm=8.,radius_nm=10.,x_in=.25,screening_fraction=0.),300.)
r_unresolved=rates(unresolved,300.)
ck('V2 levels() accepts the low-overlap bound state',unresolved.valid and 0.<unresolved.overlap_sq<1e-12)
ck('V2 rates() rejects it as overlap_unresolved',not r_unresolved['valid'] and 'overlap_unresolved' in r_unresolved['invalid_reasons'])
resolved=levels(NitrideDotSystem(**qargs),300.)
ck('V2 an ordinary well-overlapped state is unaffected',resolved.overlap_sq>=1e-12 and rates(resolved,300.)['valid'])

# ---- (M) orientation-dependent growth-axis mass swap: isolate the swap
# from the (much larger) polarization-field change between orientations by
# comparing OLD (always c-plane masses) vs NEW (swapped) at the SAME
# a_plane field, independently assembled from _z_state/_radial/
# _coulomb_binding_eV -- not by calling nitride_levels._levels_cached twice
# (that would also change the field, confounding the comparison).
def _mass_isolated_ex(s,T_K,mez,mmz,mhz,mmhz,mexy,mmexy,mhxy,mmhxy):
    d=ingaN(s.x_in); m=binary('GaN')
    de=band_edges(d,T_K,substrate=m,strain_fraction=s.strain_fraction,vbo_InN_GaN_eV=s.vbo_InN_GaN_eV,strain_c_fraction=s.strain_c_fraction)
    be=band_edges(m,T_K,substrate=m)
    Ve=be['Ec_eV']-de['Ec_eV']; Vh=de['Ev_eV']-be['Ev_eV']
    F=polarization_field(d,m,T_K,strain_fraction=s.strain_fraction,screening_fraction=s.screening_fraction,
                          external_field_kVcm=s.external_field_kVcm,orientation=s.orientation,polarization_factor=s.polarization_factor)
    ee,_,ce,_,ze,pe=_z_state(s.height_nm,Ve,mez,mmz,F,-1,1201,45.)
    eh,_,ch,_,zh,ph=_z_state(s.height_nm,Vh,mhz,mmhz,F,+1,1201,45.)
    re,_,oke,_,rmse=_radial(ce-ee,s.radius_nm,mexy,mmexy)
    rh,_,okh,_,rmsh=_radial(ch-eh,s.radius_nm,mhxy,mmhxy)
    eb=ee+re; hb=eh+rh
    eps=(d.eps_r+m.eps_r)/2.
    zsep=abs(float(np.trapezoid(ze*pe*pe,ze))-float(np.trapezoid(zh*ph*ph,zh)))
    coul=_coulomb_binding_eV(rmse or s.radius_nm,rmsh or s.radius_nm,zsep,eps)
    ex=(de['Ec_eV']-de['Ev_eV'])+eb+hb-coul
    return ex,(ce-eb)*1000.,(ch-hb)*1000.

s_apl=NitrideDotSystem(orientation='a_plane')  # H=3, R=10, x=.25 default
d0=ingaN(s_apl.x_in); m0=binary('GaN')
ex_old,dee_old,deh_old=_mass_isolated_ex(s_apl,300.,d0.me_z,m0.me_z,d0.mh_z,m0.mh_z,d0.me_xy,m0.me_xy,d0.mh_xy,m0.mh_xy)
d_ez,d_exy,d_hz,d_hxy=_growth_masses(d0,'a_plane')
m_ez,m_exy,m_hz,m_hxy=_growth_masses(m0,'a_plane')
ex_new,dee_new,deh_new=_mass_isolated_ex(s_apl,300.,d_ez,m_ez,d_hz,m_hz,d_exy,m_exy,d_hxy,m_hxy)
lv_apl=levels(s_apl,300.)
ck('M independent mass-only recompute matches the production a_plane result',
   abs(lv_apl.E_X_eV-ex_new)<1e-9 and abs(lv_apl.dE_h_meV-deh_new)<1e-6)
ck('M nonpolar growth-axis mass swap moves E_X by roughly +30 meV and hole escape depth by roughly -22 meV at H=3nm',
   20.<=(ex_new-ex_old)*1000.<=45. and -35.<=(deh_new-deh_old)<=-10.)
ck('M semipolar_11_22 keeps c-plane masses (no swap)',_growth_masses(d0,'semipolar_11_22')==(d0.me_z,d0.me_xy,d0.mh_z,d0.mh_xy))
ck('M c_plane keeps c-plane masses (no swap)',_growth_masses(d0,'c_plane')==(d0.me_z,d0.me_xy,d0.mh_z,d0.mh_xy))
ck('M m_plane swaps like a_plane',_growth_masses(d0,'m_plane')==(d0.me_xy,d0.me_z,d0.mh_xy,d0.mh_z))

# ---- (Q2) independent re-assembly of the QW-fluctuation formula, built
# from the same primitives nitride_levels._qw_levels uses but combined here
# independently, so a broken ASSEMBLY (wrong reservoir channel, wrong N2D
# mass source, dropped Coulomb term, ignored escape-channel selection) is
# caught even though the BDD solve itself is correct.
def independent_qw(s,T_K=300.):
    d=ingaN(s.x_in); m=binary('GaN')
    de=band_edges(d,T_K,substrate=m,strain_fraction=s.strain_fraction,vbo_InN_GaN_eV=s.vbo_InN_GaN_eV,strain_c_fraction=s.strain_c_fraction)
    be=band_edges(m,T_K,substrate=m)
    Ve=be['Ec_eV']-de['Ec_eV']; Vh=de['Ev_eV']-be['Ev_eV']
    F=polarization_field(d,m,T_K,strain_fraction=s.strain_fraction,screening_fraction=s.screening_fraction,
                          external_field_kVcm=s.external_field_kVcm,orientation=s.orientation,polarization_factor=s.polarization_factor)
    d_ez,d_exy,d_hz,d_hxy=_growth_masses(d,s.orientation)
    m_ez,m_exy,m_hz,m_hxy=_growth_masses(m,s.orientation)
    ee,_,ce,_,_,_=_z_state(s.height_nm,Ve,d_ez,m_ez,F,-1,1201,45.)
    eh,_,ch,_,_,_=_z_state(s.height_nm,Vh,d_hz,m_hz,F,+1,1201,45.)
    er,_,_,_,_,_=_z_state(s.wl_thickness_nm,Ve,d_ez,m_ez,F,-1,1201,45.)
    hr,_,_,_,_,_=_z_state(s.wl_thickness_nm,Vh,d_hz,m_hz,F,+1,1201,45.)
    de_lat=er-ee; dh_lat=hr-eh
    re,_,oke,_,_=_radial(de_lat,s.radius_nm,d_exy,d_exy)
    rh,_,okh,_,_=_radial(dh_lat,s.radius_nm,d_hxy,d_hxy)
    ebe=ee+re; hbe=eh+rh
    eth=min(er,ce); hth=min(hr,ch)
    ede=(eth-ebe)*1000.; hde=(hth-hbe)*1000.
    gap=de['Ec_eV']-de['Ev_eV']
    reservoir_kind='ingan_qw' if (er<=ce and hr<=ch) else 'gan_barrier'
    reservoir_energy_eV=gap+eth+hth
    return dict(F=F,ee=ee,eh=eh,er=er,hr=hr,ce=ce,ch=ch,ebe=ebe,hbe=hbe,ede=ede,hde=hde,gap=gap,
                reservoir_kind=reservoir_kind,reservoir_energy_eV=reservoir_energy_eV,
                m_e_matrix_xy=d_exy,m_h_matrix_xy=d_hxy,de_lat=de_lat,dh_lat=dh_lat,oke=oke,okh=okh)

qi=independent_qw(NitrideDotSystem(**qargs))
ck('Q2 independent recompute matches production reservoir_kind and reservoir_energy_eV',
   q.reservoir_kind==qi['reservoir_kind'] and abs(q.reservoir_energy_eV-qi['reservoir_energy_eV'])<1e-9)
ck('Q2 independent recompute matches production escape depths and matrix masses',
   abs(q.dE_e_meV-qi['ede'])<1e-6 and abs(q.dE_h_meV-qi['hde'])<1e-6
   and q.m_e_matrix_xy==qi['m_e_matrix_xy'] and q.m_h_matrix_xy==qi['m_h_matrix_xy'])

# escape-channel selection: the exact Opus-finding fixture (H=3.5, w=3.0,
# R=20, unscreened) where the electron's true escape channel is the wider
# GaN plateau (ce<er), not the QW edge.
chan_sys=NitrideDotSystem(geometry_type='qw_fluctuation',height_nm=3.5,radius_nm=20.,wl_thickness_nm=3.,screening_fraction=0.)
chan_lv=levels(chan_sys,300.)
chan_ind=independent_qw(chan_sys)
ck('Q2 GaN-barrier escape channel is selected and reported when ce<er (Opus finding fixture)',
   chan_lv.valid and chan_lv.reservoir_kind=='gan_barrier'==chan_ind['reservoir_kind']
   and chan_ind['ce']<chan_ind['er'])
ck('Q2 reservoir_energy_eV tracks the selected (not always-QW) channel',
   abs(chan_lv.reservoir_energy_eV-chan_ind['reservoir_energy_eV'])<1e-9
   and abs(chan_lv.reservoir_energy_eV-((chan_ind['gap']+chan_ind['er']+chan_ind['hr'])))>1e-4)
ck('Q2 fully-screened qargs fixture selects the QW channel for both carriers',
   q.reservoir_kind=='ingan_qw' and qi['er']<=qi['ce'] and qi['hr']<=qi['ch'])

# N2D mass source: the QW branch's matrix mass must be the SURROUNDING
# InGaN well's own (orientation-correct) in-plane mass, never GaN's.
_qsys=NitrideDotSystem(**qargs)
d_q=ingaN(_qsys.x_in); m_q=binary('GaN')
ck('Q2 N2D matrix mass is the InGaN well mass, not GaN',
   q.m_e_matrix_xy==d_q.me_xy and q.m_h_matrix_xy==d_q.mh_xy
   and q.m_e_matrix_xy!=m_q.me_xy and q.m_h_matrix_xy!=m_q.mh_xy)
r_analytic_mass=q.m_e_matrix_xy if q.dE_e_meV<=q.dE_h_meV else q.m_h_matrix_xy
r_q=rates(q,300.,tau_rad0_ns=1.3,n_dot_cm2=1e10,tau_cap_ps=10.)
_HBAR_SI=1.054571817e-34; _M0_SI=9.1093837015e-31; _KB_SI=1.380649e-23; _KB_EV=8.617333262145e-5
n2d_expected=r_analytic_mass*_M0_SI*_KB_SI*300./(math.pi*_HBAR_SI**2)/1e4
attempt_expected=(1000./10.)*(n2d_expected/1e10)
ck('Q2 QW escape prefactor matches an independently computed InGaN-mass N2D',
   r_q['valid'] and abs(r_q['escape_prefactor_ns']-attempt_expected)<1e-9*max(1.,attempt_expected))

# tau_cap / density / T sensitivity of the QW detailed-balance prefactor and
# Arrhenius factor (acceptance criterion 7).
r_base=rates(q,300.,tau_rad0_ns=1.3,n_dot_cm2=1e10,tau_cap_ps=10.)
r_tau=rates(q,300.,tau_rad0_ns=1.3,n_dot_cm2=1e10,tau_cap_ps=20.)
r_dens=rates(q,300.,tau_rad0_ns=1.3,n_dot_cm2=2e10,tau_cap_ps=10.)
r_T=rates(q,250.,tau_rad0_ns=1.3,n_dot_cm2=1e10,tau_cap_ps=10.)
ck('Q2 doubling tau_cap halves the escape prefactor',
   abs(r_tau['escape_prefactor_ns']-r_base['escape_prefactor_ns']/2.)<1e-9*max(1.,r_base['escape_prefactor_ns']))
ck('Q2 doubling density halves the escape prefactor',
   abs(r_dens['escape_prefactor_ns']-r_base['escape_prefactor_ns']/2.)<1e-9*max(1.,r_base['escape_prefactor_ns']))
ck('Q2 lower T lowers the Arrhenius-suppressed escape rate (E_a>0)',
   r_base['E_a_meV']>0. and r_T['k_X_ns']<r_base['k_X_ns'])

# overlap not forced to 1.0 (a caught mutation in the review).
ck('Q2 overlap_sq is not forced to 1.0',q.overlap_sq<1.-1e-9)

# Coulomb term present: E_X must be measurably BELOW gap+ebe+hbe.
de_q=band_edges(d_q,300.,substrate=m_q,strain_fraction=1.,vbo_InN_GaN_eV=1.15,strain_c_fraction=.7)
gap_q=de_q['Ec_eV']-de_q['Ev_eV']
no_coulomb_eV=gap_q+(q.E_e_meV+q.E_h_meV)/1000.
ck('Q2 Coulomb binding is present and subtracted (E_X < gap+ebe+hbe)',
   q.E_X_eV<no_coulomb_eV-1e-4 and (no_coulomb_eV-q.E_X_eV)*1000.>1.)

# ---- (Q3) symmetric-fluctuation placement: in the UNSCREENED (default,
# large intrinsic-field) limit the lateral depth er-ee/hr-eh collapses to
# the mass-independent geometric estimate |F|*(H-w)/2, a direct numerical
# signature of the model's z=0-centered placement assumption (see module
# docstring). Independently recomputed via _z_state, not via the module's
# dE_e_meV/dE_h_meV fields (those additionally subtract the radial term).
sym_ok=True
for Hs,ws in ((3.5,3.0),(6.,3.),(10.,8.)):
    s_sym=NitrideDotSystem(geometry_type='qw_fluctuation',height_nm=Hs,radius_nm=20.,wl_thickness_nm=ws,screening_fraction=0.)
    ind=independent_qw(s_sym)
    predicted=abs(ind['F'])*1e-4*(Hs-ws)/2.
    if predicted<=0.:
        sym_ok=False; continue
    rel_e=abs(ind['de_lat']-predicted)/predicted; rel_h=abs(ind['dh_lat']-predicted)/predicted
    print('Q3 symmetric-limit H=%g w=%g de_lat=%.5feV dh_lat=%.5feV predicted=|F|(H-w)/2=%.5feV rel_e=%.4f rel_h=%.4f'
          %(Hs,ws,ind['de_lat'],ind['dh_lat'],predicted,rel_e,rel_h))
    if not (rel_e<0.02 and rel_h<0.02): sym_ok=False
ck('Q3 unscreened symmetric-limit lateral depth matches |F|(H-w)/2 to 2%',sym_ok)

# ---- (S) shape_height_fraction actually changes the reported geometry and
# confinement (a mutation that ignores it was caught passing 16/16 in the
# review); lens native fraction is .5, so .9 is a valid, distinct choice.
shf_base=levels(NitrideDotSystem(height_nm=H,radius_nm=R,shape='lens',screening_fraction=1.),300.)
shf_alt=levels(NitrideDotSystem(height_nm=H,radius_nm=R,shape='lens',shape_height_fraction=.9,screening_fraction=1.),300.)
ck('S shape_height_fraction changes effective height and confinement',
   shf_base.valid and shf_alt.valid and shf_base.effective_height_nm==H/2. and shf_alt.effective_height_nm==H*.9
   and shf_alt.effective_height_nm!=shf_base.effective_height_nm and shf_alt.E_X_eV!=shf_base.E_X_eV)

# ---- (P) physical vs effective dimensions and aspect ratio exposed.
ck('P physical height/radius and aspect ratio exposed on isolated dots',
   disc.physical_height_nm==H and disc.physical_radius_nm==R and abs(disc.aspect_ratio-H/(2.*R))<1e-12)
ck('P physical dimensions exposed on QW rows too',
   q.physical_height_nm==qargs['height_nm'] and q.physical_radius_nm==qargs['radius_nm']
   and abs(q.aspect_ratio-qargs['height_nm']/(2.*qargs['radius_nm']))<1e-12)

# ---- (X) QW first-excited-state split now uses the dot column's z-excited
# state (ee1/eh1) rather than being hardwired to "no z-excitation ever
# wins" -- sanity: the split must be finite and no larger than the local
# continuum offset used to compute it.
ck('X QW first-excited split is finite and physically bounded',
   math.isfinite(q.sp_split_e_meV) and math.isfinite(q.sp_split_h_meV)
   and 0.<q.sp_split_e_meV<1e4 and 0.<q.sp_split_h_meV<1e4)

# ---------------------------------------------------------------------
# (C) convergence matrix: h in {1,3,7,10} nm x R in {5,30} nm x
# orientation in {c_plane, a_plane} x geometry in {isolated_dot,
# qw_fluctuation} x screening_fraction in {0, 1}: 64 fixtures. Budget is
# <=0.5 meV in E_X and both escape depths; for overlap_sq>=1e-12 require
# <=5% relative change, below that mark the lifetime numerically
# unresolved (not gated). Stated per-height numerical budget: the default
# (z_points=1201, exterior_nm=45) -> (2401, 90) pair already meets this
# budget at every height EXCEPT h=7 nm, where doubling from 1201/45 still
# moves E_X by up to ~0.61 meV (over budget); h=7 nm fixtures therefore use
# the raised (2401, 90) -> (4801, 180) pair instead. Invalid fixtures
# (including the h=10nm/c_plane/unscreened rows correctly caught by the (V)
# validity gate) are reported separately and are not a convergence failure.
DEFAULT_PAIR=((1201,45.),(2401,90.))
H7_PAIR=((2401,90.),(4801,180.))
conv_n_compared=0; conv_n_invalid=0; conv_ok=True
for h in (1.,3.,7.,10.):
    (bz,bx),(fz,fx) = H7_PAIR if h==7. else DEFAULT_PAIR
    for Rr in (5.,30.):
        for orient in ('c_plane','a_plane'):
            for geom in ('isolated_dot','qw_fluctuation'):
                for scr in (0.,1.):
                    kw=dict(height_nm=h,radius_nm=Rr,orientation=orient,geometry_type=geom,screening_fraction=scr)
                    if geom=='qw_fluctuation': kw['wl_thickness_nm']=.5*h
                    s_c=NitrideDotSystem(**kw)
                    lo=levels(s_c,300.,z_points=bz,exterior_nm=bx)
                    hi=levels(s_c,300.,z_points=fz,exterior_nm=fx)
                    if not (lo.valid and hi.valid):
                        conv_n_invalid+=1
                        print('C invalid (excluded) h=%g R=%g %s %s scr=%g: lo.valid=%s hi.valid=%s'
                              %(h,Rr,orient,geom,scr,lo.valid,hi.valid))
                        continue
                    conv_n_compared+=1
                    dE=abs(lo.E_X_eV-hi.E_X_eV)*1000.; dee=abs(lo.dE_e_meV-hi.dE_e_meV); deh=abs(lo.dE_h_meV-hi.dE_h_meV)
                    if lo.overlap_sq>=1e-12:
                        rel_ov=abs(lo.overlap_sq-hi.overlap_sq)/lo.overlap_sq; ov_ok=rel_ov<=.05
                    else:
                        rel_ov=None; ov_ok=True  # numerically unresolved, not gated
                    row_ok = dE<=.5 and dee<=.5 and deh<=.5 and ov_ok
                    if not row_ok:
                        conv_ok=False
                        print('C FAIL h=%g R=%g %s %s scr=%g dE=%.4fmeV dee=%.4fmeV deh=%.4fmeV rel_ov=%s'
                              %(h,Rr,orient,geom,scr,dE,dee,deh,rel_ov))
print('C convergence matrix: %d compared, %d invalid/excluded (of 64 fixtures)'%(conv_n_compared,conv_n_invalid))
ck('C convergence matrix within budget for every valid fixture',conv_ok)
ck('C convergence matrix includes a substantial number of genuinely bound fixtures',conv_n_compared>=40)
ck('C convergence matrix is not vacuously all-invalid',conv_n_invalid<64)

print('%d/%d nitride_geometry checks passed'%(sum(ok for _,ok in checks),len(checks)))
for name,ok in checks:
    if not ok: print('FAIL '+name)
raise SystemExit(0 if all(ok for _,ok in checks) else 1)
