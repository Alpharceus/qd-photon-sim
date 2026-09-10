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
import copy
import math
import sys
from pathlib import Path
import numpy as np
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from fsim_core.nitride_levels import (NitrideDotSystem, levels, rates,
    _z_state, _radial, _coulomb_binding_eV, _growth_masses)
from fsim_core.nitride_materials import binary, ingaN, band_edges, polarization_field
from fsim_core import nitride_materials as nitride_materials_mod
import fsim_core.device as device_mod
from fsim_core.device import DeviceDesign, evaluate

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
    # Trinary channel-selection rule (hardening round after commit 17a333f):
    # 'ingan_qw' only when BOTH carriers select the well, 'gan_barrier' only
    # when BOTH select the GaN plateau, 'mixed' when they select DIFFERENT
    # channels (the prior binary rule mislabelled that split case
    # 'gan_barrier').
    e_is_qw=er<=ce; h_is_qw=hr<=ch
    if e_is_qw and h_is_qw: reservoir_kind='ingan_qw'
    elif e_is_qw or h_is_qw: reservoir_kind='mixed'
    else: reservoir_kind='gan_barrier'
    reservoir_energy_eV=gap+eth+hth
    return dict(F=F,ee=ee,eh=eh,er=er,hr=hr,ce=ce,ch=ch,ebe=ebe,hbe=hbe,ede=ede,hde=hde,gap=gap,
                reservoir_kind=reservoir_kind,reservoir_energy_eV=reservoir_energy_eV,
                m_e_matrix_xy=d_exy,m_h_matrix_xy=d_hxy,de_lat=de_lat,dh_lat=dh_lat,oke=oke,okh=okh,
                eth=eth,hth=hth,e_is_qw=e_is_qw,h_is_qw=h_is_qw)

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

# ---------------------------------------------------------------------
# Hardening round after PASS re-review (commit 17a333f): close findings
# (a)-(g) without moving any default-setting isolated-dot/nonpolar physics
# output. The grid below is captured from commit 17a333f itself (loaded as
# a standalone scratch package, independent of the current working tree)
# BEFORE this round's edits, so it is a real regression pin against a
# frozen historical output, not a tautological re-derivation of the code
# under test.
# ---------------------------------------------------------------------

def _nan_eq(a, b):
    fa, fb = float(a), float(b)
    if math.isnan(fa) and math.isnan(fb): return True
    return fa == fb

# h_nm, R_nm, orientation, screening_fraction, E_X_eV, overlap_sq, dE_e_meV,
# dE_h_meV, valid, reservoir_kind, gamma_X0_ns, k_X_ns, S0, E_a_meV, rates_valid
PIN_GRID_17A333F = [
  (1.0,5.0,'c_plane',0.0, 2.9854444405583824,0.26886426817644143,7.140173001705497,8.745657404253475, True,'gan_barrier', 0.20681866782803188,15238.936375242862,1.3571541732784111e-05,7.140173001705497, True),
  (1.0,5.0,'c_plane',1.0, 3.118904444932365,0.8322086313626568,131.45395644445478,128.05261923264536, True,'gan_barrier', 0.6401604856635821,251.60308120866017,0.002537869721954131,128.05261923264536, True),
  (1.0,5.0,'a_plane',0.0, 3.1878944558087285,0.9962343619112091,124.0025693554127,66.6644455713944, True,'gan_barrier', 0.7663341245470838,15404.08097349792,4.9746297982966326e-05,66.6644455713944, True),
  (1.0,5.0,'a_plane',1.0, 3.1878944558087285,0.9962343619112091,124.0025693554127,66.6644455713944, True,'gan_barrier', 0.7663341245470838,15404.08097349792,4.9746297982966326e-05,66.6644455713944, True),
  (1.0,30.0,'c_plane',0.0, 2.983389895408833,0.26886426817644143,22.922165358933626,16.514582533493833, True,'gan_barrier', 0.20681866782803188,18813.625838443248,1.0992904614247245e-05,16.514582533493833, True),
  (1.0,30.0,'c_plane',1.0, 3.1223393641567427,0.8322086313626568,164.69761865325611,139.00647301737573, True,'gan_barrier', 0.6401604856635821,164.70216607276575,0.003871727820627706,139.00647301737573, True),
  (1.0,30.0,'a_plane',0.0, 3.2015503231620412,0.9962343619112091,154.05110304349938,70.45568769461968, True,'gan_barrier', 0.7663341245470838,13302.881681212599,5.760330727809496e-05,70.45568769461968, True),
  (1.0,30.0,'a_plane',1.0, 3.2015503231620412,0.9962343619112091,154.05110304349938,70.45568769461968, True,'gan_barrier', 0.7663341245470838,13302.881681212599,5.760330727809496e-05,70.45568769461968, True),
  (3.0,5.0,'c_plane',0.0, 2.1624787805772376,0.0025456662371087487,21.683708249353202,9.238168459395624, True,'gan_barrier', 0.0019582047977759606,24929.34563826353,7.855018227807029e-08,9.238168459395624, True),
  (3.0,5.0,'c_plane',1.0, 2.8626310529575254,0.9559569723028649,326.7452012484394,185.45723929806672, True,'gan_barrier', 0.73535151715605,27.31245772901792,0.02621778801695038,185.45723929806672, True),
  (3.0,5.0,'a_plane',0.0, 2.883174908429698,0.9838185368314954,322.5912125515136,168.56163006498454, True,'gan_barrier', 0.7567834898703811,299.1151665400365,0.0025236888271640795,168.56163006498454, True),
  (3.0,5.0,'a_plane',1.0, 2.883174908429698,0.9838185368314954,322.5912125515136,168.56163006498454, True,'gan_barrier', 0.7567834898703811,299.1151665400365,0.0025236888271640795,168.56163006498454, True),
  (3.0,30.0,'c_plane',0.0, 2.1547764245511436,0.0025456662371087487,43.26860224254969,17.082168499493456, True,'gan_barrier', 0.0019582047977759606,18405.07019116673,1.0639484495815067e-07,17.082168499493456, True),
  (3.0,30.0,'c_plane',1.0, 2.8635400547366605,0.9559569723028649,365.68950590012827,196.71634117156933, True,'gan_barrier', 0.73535151715605,17.669171881269083,0.03995493397123088,196.71634117156933, True),
  (3.0,30.0,'a_plane',0.0, 2.895319239125514,0.9838185368314954,357.95824490310775,172.65597742741903, True,'gan_barrier', 0.7567834898703811,255.3032985348733,0.002955491866933212,172.65597742741903, True),
  (3.0,30.0,'a_plane',1.0, 2.895319239125514,0.9838185368314954,357.95824490310775,172.65597742741903, True,'gan_barrier', 0.7567834898703811,255.3032985348733,0.002955491866933212,172.65597742741903, True),
  (7.0,5.0,'c_plane',0.0, 0.5588768913897915,6.993479097544514e-12,21.985081807827655,9.52963234635007, True,'gan_barrier', 5.3795993058034715e-12,24649.86250791069,2.1824054004670564e-16,9.52963234635007, True),
  (7.0,5.0,'c_plane',1.0, 2.7758991598998928,0.9884671514849829,400.9035436547549,197.41108163815528, True,'gan_barrier', 0.7603593472961406,17.200658333950148,0.042333867756839734,197.41108163815528, True),
  (7.0,5.0,'a_plane',0.0, 2.7729225743404613,0.9924346438089483,401.8634076184131,198.8168510504523, True,'gan_barrier', 0.7634112644684218,92.80536688641533,0.008158824765643382,198.8168510504523, True),
  (7.0,5.0,'a_plane',1.0, 2.7729225743404613,0.9924346438089483,401.8634076184131,198.8168510504523, True,'gan_barrier', 0.7634112644684218,92.80536688641533,0.008158824765643382,198.8168510504523, True),
  (7.0,30.0,'c_plane',0.0, 0.5395091096981972,6.993479097544514e-12,43.651270757072645,17.416188502642484, True,'gan_barrier', 5.3795993058034715e-12,18168.797685787624,2.96089999945985e-16,17.416188502642484, True),
  (7.0,30.0,'c_plane',1.0, 2.7761956214375996,0.9884671514849829,441.0147335220955,208.7180290745468, True,'gan_barrier', 0.7603593472961406,11.106998482578803,0.06407149410983533,208.7180290745468, True),
  (7.0,30.0,'a_plane',0.0, 2.784625658800511,0.9924346438089483,438.33135614126127,202.9552510774334, True,'gan_barrier', 0.7634112644684218,79.07715414290227,0.0095616966209253,202.9552510774334, True),
  (7.0,30.0,'a_plane',1.0, 2.784625658800511,0.9924346438089483,438.33135614126127,202.9552510774334, True,'gan_barrier', 0.7634112644684218,79.07715414290227,0.0095616966209253,202.9552510774334, True),
  (10.0,5.0,'c_plane',0.0, float('nan'),0.0,float('nan'),float('nan'), False,'gan_barrier', float('nan'),float('nan'),float('nan'),float('nan'), False),
  (10.0,5.0,'c_plane',1.0, 2.7607401823899,0.9938026469487302,414.26204873668246,199.11851242492003, True,'gan_barrier', 0.7644635745759463,16.101320489828733,0.045326299189929184,199.11851242492003, True),
  (10.0,5.0,'a_plane',0.0, 2.7536855039474144,0.995447267555754,416.42559645016,203.3837631919083, True,'gan_barrier', 0.76572866735058,77.77720721948545,0.009749172967685281,203.3837631919083, True),
  (10.0,5.0,'a_plane',1.0, 2.7536855039474144,0.995447267555754,416.42559645016,203.3837631919083, True,'gan_barrier', 0.76572866735058,77.77720721948545,0.009749172967685281,203.3837631919083, True),
  (10.0,30.0,'c_plane',0.0, float('nan'),0.0,float('nan'),float('nan'), False,'gan_barrier', float('nan'),float('nan'),float('nan'),float('nan'), False),
  (10.0,30.0,'c_plane',1.0, 2.7609383050315426,0.9938026469487302,454.55545900520525,210.43197548582563, True,'gan_barrier', 0.7644635745759463,10.394501831143588,0.0685066712532436,210.43197548582563, True),
  (10.0,30.0,'a_plane',0.0, 2.765314396547345,0.995447267555754,453.06671410160766,207.52800590761424, True,'gan_barrier', 0.76572866735058,66.2570554319813,0.011424900914526663,207.52800590761424, True),
  (10.0,30.0,'a_plane',1.0, 2.765314396547345,0.995447267555754,453.06671410160766,207.52800590761424, True,'gan_barrier', 0.76572866735058,66.2570554319813,0.011424900914526663,207.52800590761424, True),
]

pin_ok = True
_be300 = band_edges(binary('GaN'), 300., substrate=binary('GaN'))
_expected_iso_reservoir_energy_eV = _be300['Ec_eV'] - _be300['Ev_eV'] - .025
for (h, Rr, orient, scr, E_X_eV, overlap_sq, dE_e_meV, dE_h_meV, valid, reservoir_kind,
     gamma_X0_ns, k_X_ns, S0, E_a_meV, rates_valid) in PIN_GRID_17A333F:
    s_pin = NitrideDotSystem(height_nm=h, radius_nm=Rr, orientation=orient, screening_fraction=scr)
    lv_pin = levels(s_pin, 300.)
    rr_pin = rates(lv_pin, 300., tau_rad0_ns=1.3, n_dot_cm2=1e10, tau_cap_ps=10.)
    row_ok = (_nan_eq(lv_pin.E_X_eV, E_X_eV) and _nan_eq(lv_pin.overlap_sq, overlap_sq)
              and _nan_eq(lv_pin.dE_e_meV, dE_e_meV) and _nan_eq(lv_pin.dE_h_meV, dE_h_meV)
              and lv_pin.valid == valid and lv_pin.reservoir_kind == reservoir_kind
              and _nan_eq(rr_pin['gamma_X0_ns'], gamma_X0_ns) and _nan_eq(rr_pin['k_X_ns'], k_X_ns)
              and _nan_eq(rr_pin['S0'], S0) and _nan_eq(rr_pin['E_a_meV'], E_a_meV)
              and rr_pin['valid'] == rates_valid)
    if not row_ok:
        pin_ok = False
        print('PIN FAIL h=%g R=%g %s scr=%g: got E_X=%r overlap=%r dEe=%r dEh=%r valid=%r kind=%r gX0=%r kX=%r S0=%r Ea=%r rvalid=%r'
              % (h, Rr, orient, scr, lv_pin.E_X_eV, lv_pin.overlap_sq, lv_pin.dE_e_meV, lv_pin.dE_h_meV,
                 lv_pin.valid, lv_pin.reservoir_kind, rr_pin['gamma_X0_ns'], rr_pin['k_X_ns'], rr_pin['S0'],
                 rr_pin['E_a_meV'], rr_pin['valid']))
    # reservoir_energy_eV is the ONE isolated-dot output this round is
    # explicitly allowed (and required, finding d) to change: it must now
    # equal the bulk GaN Varshni edge minus the SAME 25 meV offset
    # fsim_core.device._nitride_reservoir_energy_eV applies, independent of
    # this row's own height/radius/orientation/screening (a bulk material
    # edge cannot depend on the dot's own confinement), and it must NOT
    # still equal the pre-hardening value (the dot's own strained gap).
    if lv_pin.valid:
        if abs(lv_pin.reservoir_energy_eV - _expected_iso_reservoir_energy_eV) > 1e-9:
            pin_ok = False
            print('PIN FAIL reservoir_energy_eV h=%g R=%g %s scr=%g: got %r expected %r'
                  % (h, Rr, orient, scr, lv_pin.reservoir_energy_eV, _expected_iso_reservoir_energy_eV))
ck('PIN isolated-dot c-plane/nonpolar E_X/overlap/escape-depths/valid/reservoir_kind/rates '
   'bit-identical to commit 17a333f across the (h,R,orientation,screening) grid', pin_ok)
ck('PIN isolated-dot reservoir_energy_eV is now the T-only bulk GaN edge (device.py agreement, finding d)',
   pin_ok)

# ---- (B) QW ee1/eh1 admission now uses this branch's own threshold
# eth=min(er,ce)/hth=min(hr,ch), not the isolated-dot continuum ce/ch (the
# exact Opus re-review fixture: h=7, R=5, w=3.5, screening=1, where the
# z-excited state sits 33 meV ABOVE eth but 88.69 meV below ce, so the old
# ee1<ce criterion wrongly admitted it).
b_sys = NitrideDotSystem(geometry_type='qw_fluctuation', height_nm=7., radius_nm=5.,
                          wl_thickness_nm=3.5, screening_fraction=1.)
b_lv = levels(b_sys, 300.)
d_b = ingaN(b_sys.x_in); m_b = binary('GaN')
de_b = band_edges(d_b, 300., substrate=m_b, strain_fraction=b_sys.strain_fraction,
                   vbo_InN_GaN_eV=b_sys.vbo_InN_GaN_eV, strain_c_fraction=b_sys.strain_c_fraction)
be_b = band_edges(m_b, 300., substrate=m_b)
Ve_b = be_b['Ec_eV'] - de_b['Ec_eV']
F_b = polarization_field(d_b, m_b, 300., strain_fraction=b_sys.strain_fraction,
                          screening_fraction=b_sys.screening_fraction, external_field_kVcm=b_sys.external_field_kVcm,
                          orientation=b_sys.orientation, polarization_factor=b_sys.polarization_factor)
d_ez_b, d_exy_b, d_hz_b, d_hxy_b = _growth_masses(d_b, b_sys.orientation)
m_ez_b, m_exy_b, m_hz_b, m_hxy_b = _growth_masses(m_b, b_sys.orientation)
ee_b, ee1_b, ce_b, _, _, _ = _z_state(b_sys.height_nm, Ve_b, d_ez_b, m_ez_b, F_b, -1, 1201, 45.)
er_b, _, _, _, _, _ = _z_state(b_sys.wl_thickness_nm, Ve_b, d_ez_b, m_ez_b, F_b, -1, 1201, 45.)
eth_b = min(er_b, ce_b)
ck('B fixture reproduces the Opus re-review geometry (ee1 above eth, below the old ce threshold)',
   b_lv.valid and ee1_b > eth_b and (ee1_b - ee_b) * 1000. > 50. and (ee1_b - eth_b) * 1000. < 50.)
ck('B fixed: ee1 above eth is NOT admitted as the z-excited split (sp_split_e_meV is NaN, not the old 88.69 meV)',
   math.isnan(b_lv.sp_split_e_meV))

# ---- (C) reservoir_electron_edge_meV/reservoir_hole_edge_meV report the
# SELECTED per-carrier channel (eth,hth), not the raw surrounding-well edge
# (er,hr) regardless of which channel was actually chosen.
ck('C QW row (both carriers select the well): edge fields equal the (QW) selected channel',
   abs(q.reservoir_electron_edge_meV - qi['eth'] * 1000.) < 1e-6
   and abs(q.reservoir_hole_edge_meV - qi['hth'] * 1000.) < 1e-6
   and abs(q.reservoir_electron_edge_meV - qi['er'] * 1000.) < 1e-6)
ck('C GaN-barrier row (both carriers select the plateau): edge fields equal the (barrier) selected '
   'channel, not the raw surrounding-well edge',
   abs(chan_lv.reservoir_electron_edge_meV - chan_ind['eth'] * 1000.) < 1e-6
   and abs(chan_lv.reservoir_hole_edge_meV - chan_ind['hth'] * 1000.) < 1e-6
   and abs(chan_lv.reservoir_electron_edge_meV - chan_ind['ce'] * 1000.) < 1e-6
   and abs(chan_lv.reservoir_electron_edge_meV - chan_ind['er'] * 1000.) > 1.)

# ---- (E) mixed-carrier reservoir_kind: electron and hole select DIFFERENT
# channels (found by scanning H/R/w/screening for e_is_qw != h_is_qw).
mixed_sys = NitrideDotSystem(geometry_type='qw_fluctuation', height_nm=3.5, radius_nm=20.,
                              wl_thickness_nm=1.75, screening_fraction=0.7)
mixed_lv = levels(mixed_sys, 300.)
mixed_ind = independent_qw(mixed_sys)
ck('E mixed fixture: production and independent recompute agree the carriers select different channels',
   mixed_lv.valid and mixed_ind['e_is_qw'] != mixed_ind['h_is_qw'])
ck("E mixed fixture: reservoir_kind is the new 'mixed' label (old rule mislabelled this split 'gan_barrier')",
   mixed_lv.reservoir_kind == 'mixed' == mixed_ind['reservoir_kind'])
ck('E mixed fixture: the two per-carrier edge fields report the two DIFFERENT selected channels',
   abs(mixed_lv.reservoir_electron_edge_meV - mixed_ind['eth'] * 1000.) < 1e-6
   and abs(mixed_lv.reservoir_hole_edge_meV - mixed_ind['hth'] * 1000.) < 1e-6
   and abs(mixed_lv.reservoir_electron_edge_meV - mixed_lv.reservoir_hole_edge_meV) > 1.)

# device.py's consumer: reservoir_kind is only ever forwarded as a label,
# reservoir_energy_eV only as a number (see nitride_levels module
# docstring); a 'mixed' row must evaluate() cleanly through the SAME code
# path a real QW card uses, not raise or silently fall back. Load an actual
# shipped QW card as the base design (never mutate the on-disk file) so the
# rest of the schema -- diode, thermal, cavity, drive -- is a realistic,
# already-validated device rather than a hand-built fixture; only the dot's
# own geometry (and the matched drive.diode.wl_thickness_nm the QW
# consistency guard requires) is overridden to reach the mixed regime.
# resolve_bias is patched to zero the diode's own physical field (the same
# pattern verify_nitride_geometry_device.py uses for other bias mutations)
# so the operating-point field matches the intrinsic-only field the mixed
# fixture above was found at, rather than depending on the diode's bias
# solve landing at exactly 300 K.
_mixed_design = DeviceDesign.load(str(ROOT / "cards" / "nitride-qw-fluctuation-pulse-design.yaml"))
_mixed_design = copy.deepcopy(_mixed_design)
_mixed_design.nitride["dot"].update(height_nm=3.5, radius_nm=20., wl_thickness_nm=1.75, screening_fraction=0.7)
_mixed_design.drive.diode["wl_thickness_nm"] = 1.75

def _zero_applied_field(diode, *, T_j_K, current_uA=None, junction_voltage_V=None, field_polarity=1, external_field_kVcm=0.0):
    from fsim_core.nitride_stark import resolve_bias as _real_resolve_bias
    out = dict(_real_resolve_bias(diode, T_j_K=T_j_K, current_uA=current_uA, junction_voltage_V=junction_voltage_V,
                                   field_polarity=field_polarity, external_field_kVcm=external_field_kVcm))
    if out["bias_valid"]:
        out["applied_field_kVcm"] = 0.0
    return out

_real_resolve_bias_fn = device_mod.resolve_bias
device_mod.resolve_bias = _zero_applied_field
try:
    _mixed_scalars = evaluate(_mixed_design, [300.])["scalars"]
finally:
    device_mod.resolve_bias = _real_resolve_bias_fn
ck("E device.evaluate() on a scratch QW design in the mixed regime returns a valid row, not a crash",
   _mixed_scalars["valid"] and _mixed_scalars["geometry_type"] == "qw_fluctuation")
ck("E device.evaluate() forwards the 'mixed' label unmodified (a string device.py never branches on)",
   _mixed_scalars["reservoir_kind"] == "mixed")

# ---- (F) nitride_materials.py's orientation_factor docstring no longer
# claims nitride_levels uses c-plane masses for nonpolar growth (false
# since commit 17a333f's _growth_masses swap).
_of_doc = nitride_materials_mod.orientation_factor.__doc__ or ''
ck('F orientation_factor docstring no longer makes the stale c-plane-masses-for-nonpolar claim',
   'nitride_levels uses the c-plane GaN masses' not in _of_doc)
ck('F orientation_factor docstring now documents the growth-axis mass swap nitride_levels performs',
   '_growth_masses' in _of_doc and '17a333f' in _of_doc)

# ---- (G) the two QW cards' provenance text no longer claims
# reservoir_kind='ingan_qw' at their own default (screening_fraction=0.0)
# geometry, where the resolved value is 'gan_barrier'.
for _card_path in (ROOT / "cards" / "nitride-qw-fluctuation-pulse-design.yaml",
                    ROOT / "cards" / "nitride-qw-fluctuation-set-design.yaml"):
    _card_text = _card_path.read_text(encoding='utf-8')
    _card_flat = _card_text.replace('\n', ' ')
    while '  ' in _card_flat:
        _card_flat = _card_flat.replace('  ', ' ')
    ck('G %s provenance text states the correct gan_barrier default (not the stale ingan_qw claim)'
       % _card_path.name,
       "reservoir_kind='gan_barrier'" in _card_flat and 'screening_fraction=0.0) BOTH' in _card_flat)

# Astra review fixtures: selected plateau masses enter the bulk-channel DOS
# prefactor, and radial p states must lie below that same selected continuum.
_chan_rates = rates(chan_lv, 300., tau_rad0_ns=1.3, n_dot_cm2=1e10, tau_cap_ps=10.)
_chan_mass = chan_lv.escape_e_matrix_xy if chan_lv.dE_e_meV <= chan_lv.dE_h_meV else chan_lv.escape_h_matrix_xy
_chan_prefactor = (1000./10.) * (_chan_mass * _M0_SI * _KB_SI * 300. /
    (math.pi * _HBAR_SI**2) / 1e4 / 1e10)
ck('Astra GaN escape channel exposes GaN reservoir masses for detailed balance',
   chan_lv.escape_h_matrix_xy == m_q.mh_xy and chan_lv.escape_h_matrix_xy != d_q.mh_xy)
ck('Astra detailed-balance prefactor uses the selected escape-channel mass',
   _chan_rates['valid'] and abs(_chan_rates['escape_prefactor_ns']-_chan_prefactor) < 1e-9*max(1.,_chan_prefactor))
_radial_fixture = levels(NitrideDotSystem(geometry_type='qw_fluctuation', height_nm=3.5,
    radius_nm=5., wl_thickness_nm=3., screening_fraction=0.), 300.)
ck('Astra radial excited states above selected escape thresholds are rejected',
   _radial_fixture.valid and math.isnan(_radial_fixture.sp_split_e_meV)
   and math.isnan(_radial_fixture.sp_split_h_meV))

print('%d/%d nitride_geometry checks passed'%(sum(ok for _,ok in checks),len(checks)))
for name,ok in checks:
    if not ok: print('FAIL '+name)
raise SystemExit(0 if all(ok for _,ok in checks) else 1)
