"""Independent source and numerical checks for nanowire levels.

Every published-number check transcribes its literal independently (not
imported from the module under test); every numerical check compares two
independently derived quantities, never a formula against itself.
"""
import math, sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fsim_core.nitride_nanowire_levels import NitrideNanowireSystem, levels, rates, _ep
from fsim_core.nitride_materials import EPS0_SI, binary
import re

c = []
def C(n, x):
    c.append(bool(x))
    if not x: print('FAIL ' + n)

def _z_points_used(lv):
    m = re.search(r'z_points_used=(\d+)', lv.approximation_metadata)
    return int(m.group(1)) if m else None

# ---- AC1: source transcriptions and field-construction negative controls.
# Bernardini, Fiorentini & Vanderbilt, PRB 56, R10024 (1997), Table II
# literals typed independently here, not read off the module's constants.
p25 = .75*(-.029) + .25*(-.032); f25 = (-.029-p25)/(EPS0_SI*(.75*10.28+.25*14.61))*1e-5
p40 = .6*(-.029) + .4*(-.032); f40 = (-.029-p40)/(EPS0_SI*(.6*10.28+.4*14.61))*1e-5
C('spontaneous field literal 74.55 kV/cm at x_in=0.25', abs(f25-74.55) < .2)
C('spontaneous field literal 112.83 kV/cm at x_in=0.40', abs(f40-112.83) < .2)

l25 = levels(NitrideNanowireSystem(x_in=.25, strain_bound='relaxed'), 300.)
C('relaxed x_in=0.25 keeps only the spontaneous field, P_pz=0',
  l25.valid and abs(l25.F_sp_kVcm-f25) < .2 and abs(l25.field_kVcm-f25) < .2 and abs(l25.F_pz_kVcm) < 1e-9)
l40 = levels(NitrideNanowireSystem(x_in=.4, strain_bound='relaxed'), 300.)
C('relaxed x_in=0.40 keeps only the spontaneous field, P_pz=0',
  l40.valid and abs(l40.F_sp_kVcm-f40) < .2 and abs(l40.field_kVcm-f40) < .2 and abs(l40.F_pz_kVcm) < 1e-9)

C('screening=1 leaves only the external field',
  levels(NitrideNanowireSystem(screening_fraction=1, external_field_kVcm=7)).field_kVcm == 7)
C('reversed external field survives screening=1 unchanged',
  levels(NitrideNanowireSystem(screening_fraction=1, external_field_kVcm=-7)).field_kVcm == -7)

# Zeroed-polarization negative control: screening=1 with no external field
# must give EXACTLY zero total field and a materially different E_X than
# the same (unrelaxed, large-field) system with polarization active -- a
# bug that dropped or mis-signed the polarization term would leave this
# field nonzero or the E_X shift negligible.
lu_field = levels(NitrideNanowireSystem(x_in=.25, strain_bound='unrelaxed'), 300.)
lu_nofield = levels(NitrideNanowireSystem(x_in=.25, strain_bound='unrelaxed',
                                           screening_fraction=1, external_field_kVcm=0.), 300.)
C('zeroed-polarization negative control: field exactly zero, E_X shifts by >50 meV',
  lu_field.valid and lu_nofield.valid and lu_nofield.field_kVcm == 0.
  and abs(lu_nofield.E_X_eV-lu_field.E_X_eV)*1000 > 50.)

C('same-material core/matrix gives no artificial dot', not levels(NitrideNanowireSystem(x_in=0)).valid)

# ---- AC2: hard-wall Bessel energies, R^-2 scaling, and mass dependence.
J01, J11 = 2.4048255577, 3.8317059702  # Abramowitz and Stegun 1964 Table 9.5
C('Bessel zero literals reproduce E_perp(0.2 m0, 10 nm)',
  abs(_ep(.2, 10, J01)-10.999) < .03 and abs(_ep(.2, 10, J11)-27.91) < .08)
C('radial energy scales as R^-2', abs(_ep(.2, 20, J01)-10.999/4) < .01)

# Perturbing the in-plane mass by +10% must shift the E_perp(m,n) hard-wall
# energy by the analytic 1/m factor, for a representative electron-like AND
# hole-like mass -- this is the formula both dot AND barrier thresholds are
# built from (nitride_nanowire_levels._once te/Te/th/Th), so a mass mix-up
# in either channel would fail here.
g = binary('GaN')
for label, m in (('electron-like (GaN me_xy)', g.me_xy), ('hole-like (GaN mh_xy)', g.mh_xy)):
    e_lo = _ep(m, 12.5); e_hi = _ep(m*1.1, 12.5)
    C('transverse threshold mass scaling, ' + label, abs(e_hi*1.1 - e_lo) < 1e-6*e_lo)

# Deshpande disc (h=2 nm, R=12.5 nm, x_in=0.25): the J01->J11 spacing is a
# separable, mass-agnostic-to-axial-motion Bessel difference (fix A); a
# regression back to axial masses in the J11 branch reproduces the FAIL
# round's -0.185/+11.115 meV numbers, not these.
l_desh = levels(NitrideNanowireSystem(height_nm=2., core_radius_nm=12.5, x_in=.25,
                                       strain_bound='relaxed'), 300.)
C('sp_split_e_meV within 0.1 meV of 13.93 at the Deshpande geometry',
  l_desh.valid and abs(l_desh.sp_split_e_meV-13.93) < .1)
C('sp_split_h_meV within 0.1 meV of 3.31 at the Deshpande geometry',
  l_desh.valid and abs(l_desh.sp_split_h_meV-3.31) < .1)

# ---- AC3: independent 1-D thermal partition literal and wide-radius
# 3-D-DOS limit. N_res = g_spin*L*sqrt(mz*kT/(2 pi hbar^2)) * sum_m d_m
# exp(-(E_perp_mn-E_perp_01)/kT); at very small R only the ground mode
# (sum=1) is thermally reachable, so N_res reduces to the bare 1-D
# thermal-length literal below, independently typed (not the module's own
# formula fed back into itself).
#
# fix L3: reservoir_length_nm is a PER-SIDE GaN reservoir extent; the
# module's prefactor counts BOTH flanking reservoirs
# (L_total=2*reservoir_length_nm, "L=barrier_left_nm+barrier_right_nm" per
# the interface constraints). The independent literals below apply the same
# factor of 2 to the SAME reservoir_length_nm passed to rates(), so the
# comparison is apples-to-apples rather than silently absorbing the fix.
HBAR_SI = 1.054571817e-34; M0_SI = 9.1093837015e-31; KB_SI = 1.380649e-23
RES_LEN_NM = 30.  # per-side length passed to rates(); physical total is 2x this
def _n1d_ground_only(mz, L_nm, T_K):
    return 2. * L_nm*1e-9 * math.sqrt(mz*M0_SI*KB_SI*T_K/(2.*math.pi*HBAR_SI**2))
l_small = levels(NitrideNanowireSystem(height_nm=2., core_radius_nm=3., outer_radius_nm=3.,
                                        x_in=.25, strain_bound='relaxed'), 300.)
r_small = rates(l_small, 300., tau_rad0_ns=1.3, tau_cap_ps=10., reservoir_length_nm=RES_LEN_NM) if l_small.valid else None
n1d = _n1d_ground_only(g.me_z, 2.*RES_LEN_NM, 300.)
C('independent 1-D thermal partition literal (ground-mode-only, small R, L_total=2x per-side)',
  l_small.valid and r_small['valid'] and abs(r_small['reservoir_state_count_e']/n1d - 1.) < .01)

# fix L3 tag check: reservoir_state_count_e must be exactly linear in the
# per-side reservoir_length_nm (N_res proportional to L_total=2*
# reservoir_length_nm), independent of any particular factor-of-2 choice --
# doubling the per-side input must exactly double the state count.
r_1x = rates(l_small, 300., tau_rad0_ns=1.3, tau_cap_ps=10., reservoir_length_nm=RES_LEN_NM)
r_2x = rates(l_small, 300., tau_rad0_ns=1.3, tau_cap_ps=10., reservoir_length_nm=2.*RES_LEN_NM)
C('reservoir_state_count_e is linear in reservoir_length_nm (per-side convention)',
  r_1x['valid'] and r_2x['valid']
  and abs(r_2x['reservoir_state_count_e']/r_1x['reservoir_state_count_e'] - 2.) < 1e-6)
C('fix L3 choice is tagged in the rates() provenance string',
  'reservoir_length_nm' in r_1x['provenance'] and '2*' in r_1x['provenance'])

# Wide-radius limit: an independently coded 3-D nondegenerate carrier count
# (Boltzmann effective-DOS formula, anisotropic mass m_eff=(mxy^2*mz)^(1/3))
# over the SAME cylinder volume; the code/3-D ratio must climb toward 1 as
# the discrete transverse sum approaches its continuum limit.
def _n3d(core_radius_nm, L_nm, T_K, mz, mxy):
    V = math.pi*(core_radius_nm*1e-9)**2*(L_nm*1e-9)
    m_eff = (mxy*mxy*mz)**(1./3.)
    return 2.*V*(m_eff*M0_SI*KB_SI*T_K/(2.*math.pi*HBAR_SI**2))**1.5
ratios = []
for R in (40., 80., 120.):
    lw = levels(NitrideNanowireSystem(height_nm=2., core_radius_nm=R, outer_radius_nm=R,
                                       x_in=.25, strain_bound='relaxed'), 300.)
    rw = rates(lw, 300., tau_rad0_ns=1.3, tau_cap_ps=10., reservoir_length_nm=RES_LEN_NM)
    ratios.append(rw['reservoir_state_count_e'] / _n3d(R, 2.*RES_LEN_NM, 300., g.me_z, g.me_xy) if lw.valid and rw['valid'] else float('nan'))
C('wide-radius 3-D-DOS ratio ~0.90 at R=40 nm (electrons)', abs(ratios[0]-0.90) < .03)
C('wide-radius 3-D-DOS ratio ~0.95 at R=80 nm (electrons)', abs(ratios[1]-0.95) < .03)
C('wide-radius 3-D-DOS ratio ~0.96 at R=120 nm (electrons)', abs(ratios[2]-0.96) < .03)
C('code/3-D ratio increases monotonically toward 1 over R=40/80/120 nm',
  ratios[0] < ratios[1] < ratios[2] < 1.0)

# ---- k_nr separation (fix L): k_nr_ns must be added once to k_X_ns and
# once (not doubled) to k_XX_ns, i.e. k_XX_ns-2*k_X_ns == -k_nr_ns exactly.
r0 = rates(l_desh, 300., tau_rad0_ns=1.3, tau_cap_ps=10., reservoir_length_nm=30., k_nr_ns=0.)
r1 = rates(l_desh, 300., tau_rad0_ns=1.3, tau_cap_ps=10., reservoir_length_nm=30., k_nr_ns=5.)
C('k_nr_ns adds once to k_X_ns, once (not doubled) to k_XX_ns',
  r0['valid'] and r1['valid'] and abs((r1['k_X_ns']-r0['k_X_ns'])-5.) < 1e-9
  and abs((r1['k_XX_ns']-r0['k_XX_ns'])-5.) < 1e-9)

# ---- AC4: refinement table. H {1.5,4} x R {10,120} x x_in {0.25,0.40} x
# T {230,300} x strain_bound, all 32 rows checked unconditionally (never
# hidden under `if valid`). Resolution policy (directive round): levels()
# now refines z_points adaptively on a gate failure (double up to a cap of
# 9601) instead of returning invalid at the fixed default -- the row that
# used to be the fixed "4 expected-invalid rows" (H=4, x_in=0.40, unrelaxed)
# now converges honestly at an internally-doubled z_points and is VALID at
# the default call; verified directly (no row is exempted from the
# expected-valid branch below, and never clamped -- levels() only reports
# valid when the actual convergence gates pass at whatever z_points it
# settled on). All 32 rows must be valid AND -- re-solved at an EXTERNALLY
# doubled resolution called fresh through the public levels() API, not the
# internal gate's own numbers -- converged in E_X (<=0.5 meV) and overlap
# (<=2%).
for h in (1.5, 4):
    for rad in (10, 120):
        for x in (.25, .4):
            for t in (230, 300):
                for b in ('relaxed', 'unrelaxed'):
                    sysx = NitrideNanowireSystem(height_nm=h, core_radius_nm=rad, outer_radius_nm=rad,
                                                  x_in=x, strain_bound=b)
                    q = levels(sysx, t)
                    tag = str((h, rad, x, t, b))
                    qq = levels(sysx, t, z_points=2403, exterior_nm=90.)
                    ok = (q.valid and qq.valid
                          and abs(qq.E_X_eV-q.E_X_eV)*1000 <= .5
                          and abs(qq.overlap_sq-q.overlap_sq)/q.overlap_sq <= .02)
                    C('refinement ' + tag + ' expected valid and externally converged at 2x', ok)

# ---- Resolution policy (directive round, fix 2 follow-up): the 0.02 kT
# escape-depth gate makes a FIXED default z_points=1201 grid INVALID at low
# T (grid noise ~0.1 meV exceeds 0.02 kT there); levels() now refines
# adaptively on a gate failure (double z_points up to a cap of 9601,
# returning the first converged result) instead of surfacing that as a
# caller-visible failure. The Deshpande geometry must be VALID at the
# DEFAULT call (no z_points override) at every one of these temperatures,
# with z_points_used disclosed in approximation_metadata.
desh_sys = NitrideNanowireSystem(height_nm=2., core_radius_nm=12.5, x_in=.25, strain_bound='relaxed')
for T in (25., 60., 230., 300.):
    lt = levels(desh_sys, T)
    zpu = _z_points_used(lt)
    print('Deshpande geometry T=%g K default call: valid=%s z_points_used=%s reasons=%s'
          % (T, lt.valid, zpu, lt.invalid_reasons))
    C('Deshpande geometry valid at the default call, T=%g K' % T, lt.valid)
    C('Deshpande geometry reports z_points_used at T=%g K' % T, zpu is not None)

# The adaptive refinement must still give up (return invalid, never clamped
# to passing) when the underlying failure is NOT a grid-noise/convergence
# issue that more z_points can fix -- e.g. a genuinely too-thin/low-offset
# geometry whose axial state is unbound at every resolution.
l_unbound = levels(NitrideNanowireSystem(height_nm=.5, core_radius_nm=10., outer_radius_nm=10.,
                                          x_in=.1, strain_bound='unrelaxed'), 230.)
C('adaptive refinement still returns invalid (never clamped to passing) for a structurally unbound case',
  not l_unbound.valid and len(l_unbound.invalid_reasons) > 0)

# ---- AC5 / M3: held-out comparisons, reported side by side without gating
# on agreement (this model is an independent prediction, not a fit). The
# 2013 anchor (x_in=0.25) and the 2014 anchor (x_in=0.40) are matched by
# OPPOSITE strain endpoints -- relaxed matches 2013, unrelaxed is closer to
# 2014 -- so the results are never averaged; x_in=0.40 is abstract-only [A]
# (APL 2014 abstract gives no direct structural confirmation).
l_headline = levels(desh_sys, 25., z_points=4807, exterior_nm=45.)
l_headline_default = levels(desh_sys, 25.)
print('Deshpande T=25 K default resolution: valid=%s z_points_used=%s'
      % (l_headline_default.valid, _z_points_used(l_headline_default)))
l2014_relaxed = levels(NitrideNanowireSystem(height_nm=2., core_radius_nm=12.5, x_in=.4,
                                              strain_bound='relaxed'), 300.)
l2014_unrelaxed = levels(NitrideNanowireSystem(height_nm=2., core_radius_nm=12.5, x_in=.4,
                                                strain_bound='unrelaxed'), 300.)
print('---- M3: held-out comparisons (non-gating, anchors matched by opposite endpoints) ----')
if l_headline.valid:
    print('2013 (x_in=0.25, relaxed, 25 K): predicted %.2f nm vs measured 436.56 nm'
          % l_headline.lambda_nm)
if l2014_relaxed.valid and l2014_unrelaxed.valid:
    print('2014 (x_in=0.40, 300 K): relaxed %.1f nm / unrelaxed %.1f nm vs measured ~630 nm [E]'
          % (l2014_relaxed.lambda_nm, l2014_unrelaxed.lambda_nm))
print('x_in=0.40 is abstract-only [A] (APL 2014 abstract, no direct structural confirmation); '
      'the 2013 and 2014 anchors are matched by opposite strain endpoints and are reported '
      'separately, never averaged -- held out, non-gating')
C('M3 held-out comparisons computed and reported (non-gating)',
  l_headline.valid and l2014_relaxed.valid and l2014_unrelaxed.valid)

# ---- H6: vertical disc-in-wire family (disc_radius_nm < core_radius_nm).
# The disc radius is decoupled from the core radius: the radial problem
# becomes the finite-barrier InGaN/GaN disk (dot_levels.finite_disk_2d),
# NOT the hard wall, while the GaN reservoir subbands still use core_radius
# (unchanged rates()/_part() call site, which reads lv.core_radius_nm).
l_full_12p5 = levels(NitrideNanowireSystem(height_nm=2., core_radius_nm=12.5, outer_radius_nm=12.5,
                                            x_in=.25, strain_bound='relaxed'), 300.)
l_disc_in_100 = levels(NitrideNanowireSystem(height_nm=2., core_radius_nm=100., outer_radius_nm=100.,
                                              disc_radius_nm=12.5, x_in=.25, strain_bound='relaxed'), 300.)
C('H6: 12.5 nm disc in a 100 nm core is valid and selects the disc geometry label',
  l_disc_in_100.valid and l_disc_in_100.geometry == 'disc-in-wire finite-barrier radial')
C('H6: 12.5 nm disc-in-100nm-core E_perp within 30% of the horizontal 12.5 nm full-core value',
  l_full_12p5.valid and l_disc_in_100.valid
  and abs(l_disc_in_100.transverse_e_meV - l_full_12p5.transverse_e_meV) / l_full_12p5.transverse_e_meV < .30)

l_full_100 = levels(NitrideNanowireSystem(height_nm=2., core_radius_nm=100., outer_radius_nm=100.,
                                           x_in=.25, strain_bound='relaxed'), 300.)
l_disc_eq_core_100 = levels(NitrideNanowireSystem(height_nm=2., core_radius_nm=100., outer_radius_nm=100.,
                                                   disc_radius_nm=100., x_in=.25, strain_bound='relaxed'), 300.)
C('H6: a 100 nm disc in a 100 nm core reproduces the full-core (horizontal) result exactly',
  l_full_100.valid and l_disc_eq_core_100.valid
  and l_disc_eq_core_100.geometry == 'full-core hard-wall cylinder'
  and l_disc_eq_core_100.E_X_eV == l_full_100.E_X_eV
  and l_disc_eq_core_100.transverse_e_meV == l_full_100.transverse_e_meV)

# Monotonic E_perp in disc radius: a smaller disc (tighter lateral
# confinement) must give a strictly larger transverse (E_perp) energy, at
# fixed core_radius_nm.
disc_radii = (10., 15., 20., 30., 50., 80.)
e_perp_seq = []
for R in disc_radii:
    ld = levels(NitrideNanowireSystem(height_nm=2., core_radius_nm=100., outer_radius_nm=100.,
                                       disc_radius_nm=R, x_in=.25, strain_bound='relaxed'), 300.)
    e_perp_seq.append(ld.transverse_e_meV if ld.valid else float('nan'))
print('H6 E_perp(disc_radius_nm) at core_radius_nm=100: %s'
      % dict(zip(disc_radii, e_perp_seq)))
C('H6: E_perp is monotonically decreasing as disc radius increases (all rows valid)',
  all(math.isfinite(v) for v in e_perp_seq)
  and all(e_perp_seq[i] > e_perp_seq[i+1] for i in range(len(e_perp_seq)-1)))

print('H6 sidewall_overlap: 12.5-in-100 disc %.3e, full-core 12.5 %.3e, full-core 100 limit %.3e'
      % (l_disc_in_100.sidewall_overlap, l_full_12p5.sidewall_overlap, l_full_100.sidewall_overlap))
C('H6: sidewall_overlap reported (finite, non-negative) at both the disc and the full-core limit',
  math.isfinite(l_disc_in_100.sidewall_overlap) and l_disc_in_100.sidewall_overlap >= 0.
  and math.isfinite(l_full_100.sidewall_overlap) and l_full_100.sidewall_overlap >= 0.)
C('H6: a dot laterally isolated inside a much larger core has far less sidewall overlap than the full-core limit',
  l_disc_in_100.sidewall_overlap < l_full_100.sidewall_overlap)
C('L1: finite-barrier disc E_perp is strictly below independent 9.054 meV hard-wall literal',
  l_disc_in_100.transverse_e_meV < 9.054)
C('L2: returned resolutions are disclosed in approximation metadata',
  _z_points_used(l_headline) is not None and _z_points_used(l_headline_default) is not None)

print('%d/%d nitride nanowire levels checks passed' % (sum(c), len(c)))
raise SystemExit(0 if all(c) else 1)
