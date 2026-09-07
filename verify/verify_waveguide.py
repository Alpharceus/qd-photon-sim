"""Independent numerical checks for fsim_core.waveguide (published-class values)."""
import inspect
import sys
from pathlib import Path
from math import pi, tan, exp
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import fsim_core.waveguide as wgmod
from fsim_core.waveguide import (Layer, slab_modes, effective_index_ridge, beta_factor,
                                 facet_transmission, na_collection, hkust_ridge_stack,
                                 edge_emission, na_collection_numeric, facet_escape_fraction)

checks = []
def ck(ok, name):
    checks.append(bool(ok))
    if not ok: print("FAIL", name)

# Independent even-TE dispersion: u tan(u)=w, u^2+w^2=V^2.
lam, nc, nf, d = 1000., 1.45, 1.50, 1000.
k0 = 2*pi/lam
V = k0*d/2*np.sqrt(nf*nf-nc*nc)
u = brentq(lambda q: q*tan(q)-np.sqrt(V*V-q*q), 1e-10, min(V-1e-10, pi/2-1e-10))
analytic = np.sqrt(nf*nf-(2*u/(k0*d))**2)
mode = slab_modes([Layer('lo',nc,2000),Layer('core',nf,d),Layer('hi',nc,2000)],lam)[0]
ck(abs(mode.n_eff-analytic)<1e-4, 'analytic TE')
thin1 = slab_modes([Layer('lo',1.45,2000),Layer('dot',1.50,4,True),Layer('core',1.50,500),Layer('hi',1.45,2000)],1000)[0]
thin2 = slab_modes([Layer('lo',1.45,2000),Layer('dot',1.50,8,True),Layer('core',1.50,500),Layer('hi',1.45,2000)],1000)[0]
ck(1.7 < thin2.gamma_layer('dot')/thin1.gamma_layer('dot') < 2.3, 'thin gamma')
ridge = effective_index_ridge([Layer('lo',3.4,1000),Layer('core',3.5,300),Layer('hi',3.4,1000)],670,2000,1000)
F,b = beta_factor(ridge.A_mode_um2,670,3.5,ridge.n_eff)
ck(.003 < b < .03, 'ridge beta class')
A=(3/(4*pi))*(.67/3.5)**2*(3.6/3.5)
F1,b1=beta_factor(A,670,3.5,3.6)
ck(abs(F1-1)<1e-12, 'F unity')
ck(abs(b1-.5)<1e-12, 'beta half')
T=facet_transmission(3.4)
ck(.68<T<.72, 'facet')
ck(facet_transmission(3.4,.9)==.9, 'coating')
n0=na_collection(.5,.4,670,0)
n1=na_collection(.5,.4,670,1)
nmid=na_collection(.5,.4,670,.5)
ck(n0==0, 'NA zero'); ck(n1>.99, 'NA one'); ck(0<nmid<n1, 'NA monotonic')

# Round-NA aperture regression (waveguide-na-collection-fix, 2026-09-06).
# For an isotropic 1/e^2 intensity radius w0, eta=1-exp(-2 asin(NA)^2 /
# (lambda/(pi*w0))^2), independently of the implementation's quadrature.
w0, lam0, na0 = 1.5, 770.0, 0.1
eta_gaussian = na_collection(w0, w0, lam0, na0)
ck(abs(eta_gaussian - 0.527) < 0.01 * 0.527,
   'round-cone Gaussian NA=0.1 gives published review anchor 0.527')
xg = (np.arange(4096) - 2048) * 0.025
fg = np.exp(-xg**2 / w0**2)  # amplitude: I=exp(-2 x^2/w0^2)
eta_gaussian_numeric = na_collection_numeric(fg, fg, 0.025, 0.025, lam0, na0)
ck(abs(eta_gaussian_numeric - eta_gaussian) < 0.01 * eta_gaussian,
   'FFT Gaussian far field agrees with circular closed form within 1%')
ck(abs(na_collection(w0, w0, lam0, .8)
       - (1 - np.exp(-2 * np.arcsin(.8)**2 / (lam0 * 1e-3 / (pi*w0))**2))) < 1e-12,
   'equal-width Gaussian uses the analytic round-cone closed form at NA=0.8')

# The effective-width definition is (int I)^2/int I^2.  For a Gaussian it
# yields sqrt(pi)*w per axis, hence A_mode=pi*wx*wy for 1/e^2 radii.
xw = np.linspace(-8.0, 8.0, 20001)
iwx, iwy = np.exp(-2*xw**2 / 1.2**2), np.exp(-2*xw**2 / .8**2)
wx = np.trapezoid(iwx, xw)**2 / np.trapezoid(iwx**2, xw) / np.sqrt(pi)
wy = np.trapezoid(iwy, xw)**2 / np.trapezoid(iwy**2, xw) / np.sqrt(pi)
ck(abs((np.trapezoid(iwx, xw)**2 / np.trapezoid(iwx**2, xw))
       * (np.trapezoid(iwy, xw)**2 / np.trapezoid(iwy**2, xw)) - pi*wx*wy) < 1e-10,
   'Gaussian effective area equals pi wx wy for 1/e^2 intensity radii')
s=hkust_ridge_stack()
rm=effective_index_ridge(s,668,2000,1200)
g=rm.vertical.gamma_layer('dot')
ck(3.05<rm.n_eff<3.22, 'HKUST n'); ck(.005<g<.05, 'HKUST gamma')
r=edge_emission(s,2000,1200,668,500,.5)
ck(0<r.eta_total<1 and r.n_g>0, 'edge result')
ck(any('NA collection method: gaussian' in note for note in r.notes)
   and any('1/e^2 intensity radii' in note for note in r.notes),
   'edge notes state NA method and width convention')
r_na_numeric = edge_emission(s, 2000, 1200, 668, 500, .8, na_method='numeric')
r_na_gaussian = edge_emission(s, 2000, 1200, 668, 500, .8)
ck(abs(r_na_gaussian.eta_NA - r_na_numeric.eta_NA) < 0.03 * r_na_numeric.eta_NA,
   'solved GaInP ridge Gaussian NA estimate agrees with FFT oracle within 3%')
ck(any('NA collection method: numeric' in note for note in r_na_numeric.notes),
   'numeric NA method is recorded in edge notes')

# council review 2026-09-05 item 3, updated for peer-review pkg2 facet fix
# (2026-09-07, .workers/specs/pr-pkg2-facet-fix.md item 1): T_facet and the
# propagation loss are now folded into eta_total exactly once each, via the
# single continuous ray-probability series at the dot_position=0.5 default
# (0.5 * T * exp(-a*L/2) * (1 + R_back * exp(-a*L)) / (1 - R_back*R_front*
# exp(-2*a*L)), a = alpha_cm*1e-4); reproduce that by hand -- NOT via
# exp(-a*L) == sqrt(prop_rt), the pkg2-checkpoint bug used the FULL
# round-trip prop_rt for the returning term instead -- for
# r=edge_emission(..., alpha_cm=5.0 default, L_um=500 default, R_back=None)
# -- R_back=None resolves to R_back_eff=R_front=1-T_facet (no coating
# override on this stack, so the uncoated-Fresnel resolution of item 4 is
# numerically identical to R_front here).
_alpha_r, _L_r = 5.0, 500.0
_a_r = _alpha_r * 1e-4
_prop_half_r = exp(-_a_r * (_L_r / 2.0))
_return_r = exp(-_a_r * _L_r)  # single pass back facet -> front facet at x=0.5
_prop_rt_r = exp(-2.0 * _a_r * _L_r)
_R_front_r = 1.0 - r.T_facet
_eta_facet_r = (0.5 * r.T_facet * _prop_half_r * (1.0 + _R_front_r * _return_r)
               / (1.0 - _R_front_r * _R_front_r * _prop_rt_r))
expected_total = r.beta * _eta_facet_r * r.eta_NA
ck(abs(expected_total - r.eta_total) < 1e-9,
   'facet transmission (once) and propagation (folded via the ray-series at dot_position=0.5) reproduce eta_total by hand')
ck(abs(_eta_facet_r - facet_escape_fraction(r.T_facet, _R_front_r, _alpha_r, _L_r)) < 1e-12,
   'hand-derived eta_facet matches facet_escape_fraction (single source of truth) at the default dot_position')
ck(r.T_facet < 1.0, 'facet transmission is not unity (so omitting it was not a no-op)')

# Peer-review finding 3 (2026-09-07): the old two-branch facet model
# stepped discontinuously by +16.33% at T_facet=0.719371 crossing R_back=0
# (0.5*T=0.3596855 vs the escape-rate limit T/(T+1)=0.4183916 -- ABOVE
# 0.5*T, not "strictly below" as the old comment claimed:
# T/(T+1) - 0.5*T = T(1-T)/(2(T+1)) > 0). Replaced by one continuous
# ray-probability series. Every check below compares against numbers
# computed by hand right here, never by calling the facet formula inside
# edge_emission a second time.
def _facet_only(res):
    # beta and eta_NA do not depend on R_back/alpha_cm, so eta_total is
    # exactly proportional to eta_facet for a fixed stack/ridge/NA/lambda.
    return res.eta_total / (res.beta * res.eta_NA)

alpha_c, L_c = 5.0, 250.0

# (a) continuity through R_back -> 0 (this is the assertion that used to
# encode the defect: it asserted R_back=0.0 == R_back=None EXACTLY, which
# was only true because both were aliased onto the same special-cased
# branch, not because the underlying escape-rate formula is continuous).
e_tiny = _facet_only(edge_emission(s, 2000, 1200, 668, L_c, .5, alpha_cm=alpha_c, R_back=1e-9))
e_zero = _facet_only(edge_emission(s, 2000, 1200, 668, L_c, .5, alpha_cm=alpha_c, R_back=0.0))
ck(abs(e_tiny - e_zero) < 1e-8,
   'eta_facet is continuous through R_back=0 (no discontinuous branch switch)')

# (b) alpha=0, R_back=0: the series reduces exactly to the forward term 0.5*T
r_alpha0 = edge_emission(s, 2000, 1200, 668, L_c, .5, alpha_cm=0.0, R_back=0.0)
ck(abs(_facet_only(r_alpha0) - 0.5 * r_alpha0.T_facet) < 1e-12,
   'eta_facet(R_back=0, alpha=0) equals 0.5*T_facet by hand')

# (c) monotone non-decreasing in R_back on a 20-point grid, alpha=5/cm, L=250um
rb_grid = np.linspace(0.0, 0.99, 20)
facet_grid = [_facet_only(edge_emission(s, 2000, 1200, 668, L_c, .5, alpha_cm=alpha_c, R_back=rb))
             for rb in rb_grid]
ck(all(a <= b + 1e-12 for a, b in zip(facet_grid, facet_grid[1:])),
   'eta_facet is monotone non-decreasing in R_back on a 20-point grid in [0, 0.99]')

# (d) alpha=0 closed form: eta = 0.5*T*(1+R_back) / (1 - R_back*(1-T))
T0 = r_alpha0.T_facet
for rb in (0.3, 0.7, 0.999):
    r_d = edge_emission(s, 2000, 1200, 668, L_c, .5, alpha_cm=0.0, R_back=rb)
    closed_form = 0.5 * T0 * (1.0 + rb) / (1.0 - rb * (1.0 - T0))
    ck(abs(_facet_only(r_d) - closed_form) < 1e-12,
       f'eta_facet matches the hand-derived alpha=0 closed form at R_back={rb}')

# (e) R_back=None resolves to R_back=R_front=1-T_facet exactly
r_none = edge_emission(s, 2000, 1200, 668, L_c, .5, alpha_cm=alpha_c)
r_explicit = edge_emission(s, 2000, 1200, 668, L_c, .5, alpha_cm=alpha_c, R_back=1.0 - r_none.T_facet)
ck(abs(r_none.eta_total - r_explicit.eta_total) < 1e-12,
   'R_back=None equals R_back=1-T_facet (cleaved back facet resolves to the front facet Fresnel reflectivity)')

# item 4: position factor -- centred symmetric stack keeps beta near the
# antinode value (pos > 0.95); a dot moved to the cladding edge (far from the
# vertical antinode) must collapse pos well below 1.  Reuses the SAME hkust
# stack, only moving which layer is marked is_dot=True, so every layer name
# stays unique (no confinement-dict name collision).
r_centre = edge_emission(s, 2000, 1200, 668, 500, .5)
edge_stack = [Layer(x.name, x.n, x.thickness_nm, x.name == "lower_cladding") for x in s]
r_edge_dot = edge_emission(edge_stack, 2000, 1200, 668, 500, .5)
ck(r_centre.beta > 0, 'centred-dot beta positive (position factor sanity)')
ck(r_edge_dot.beta < r_centre.beta, 'dot at the cladding edge has lower beta than at the antinode')

from fsim_core.waveguide import effective_index_ridge as _eir
mode_centre = _eir(s, 668, 2000, 1200)
bounds_c = np.cumsum([0.0] + [x.thickness_nm for x in s])
dot_i = [i for i, x in enumerate(s) if x.is_dot][0]
z_dot_c = 0.5 * (bounds_c[dot_i] + bounds_c[dot_i + 1])
pos_centre = float(np.interp(z_dot_c, mode_centre.vertical.z_nm, mode_centre.vertical.field) ** 2
                   / np.max(mode_centre.vertical.field ** 2))
ck(pos_centre > 0.95, 'position factor at a symmetric core centre is near the antinode (>0.95)')

mode_edge = _eir(edge_stack, 668, 2000, 1200)
bounds_e = np.cumsum([0.0] + [x.thickness_nm for x in edge_stack])
dot_i_e = [i for i, x in enumerate(edge_stack) if x.is_dot][0]
z_dot_e = 0.5 * (bounds_e[dot_i_e] + bounds_e[dot_i_e + 1])
pos_edge = float(np.interp(z_dot_e, mode_edge.vertical.z_nm, mode_edge.vertical.field) ** 2
                 / np.max(mode_edge.vertical.field ** 2))
ck(pos_edge < 0.2, 'position factor at the cladding edge is well off the antinode (<0.2)')

# item 5: group index falls in the published-class range for the hkust stack
# at 668 nm, and a numeric-only (non-dispersive) stack keeps the n_eff
# fallback (n_g == n_eff exactly).
ck(3.6 <= r.n_g <= 4.4, 'HKUST group index n_g in the published-class range at 668 nm')
r_numeric = edge_emission([Layer('lo', 3.4, 1000), Layer('core', 3.5, 300, True),
                          Layer('hi', 3.4, 1000)], 2000, 1000, 670, 500, .5)
ck(r_numeric.n_g == r_numeric.n_eff, 'numeric-only (non-dispersive) layers keep the n_g=n_eff fallback')

# Peer-review pkg2 facet fix (2026-09-07, .workers/specs/pr-pkg2-facet-fix.md
# item 4): a `coating` override on the FRONT transmission must not leak into
# what R_back=None resolves to -- the back facet is a separate, uncoated
# semiconductor/air interface with its own Fresnel reflectivity, computed by
# hand from n_eff the same way facet_transmission computes the uncoated T
# (NOT `1 - T_coated`, the checkpoint's bug).
mode_coat = effective_index_ridge(s, 668, 2000, 1200)
R_uncoated_hand = ((mode_coat.n_eff - 1.0) / (mode_coat.n_eff + 1.0)) ** 2
alpha_coat, L_coat = 5.0, 250.0
r_coated_none = edge_emission(s, 2000, 1200, 668, L_coat, .5, alpha_cm=alpha_coat,
                              coating={'T_facet': 0.95})
eta_facet_coated_none = r_coated_none.eta_total / (r_coated_none.beta * r_coated_none.eta_NA)
eta_facet_hand_uncoated = facet_escape_fraction(0.95, R_uncoated_hand, alpha_coat, L_coat)
ck(abs(eta_facet_coated_none - eta_facet_hand_uncoated) < 1e-9,
   'coating override on T_facet: R_back=None resolves to the uncoated Fresnel '
   'reflectivity of the bare back facet (hand-computed from n_eff), not 1-T_coated')
eta_facet_wrong_1mT = facet_escape_fraction(0.95, 1.0 - 0.95, alpha_coat, L_coat)
ck(abs(eta_facet_coated_none - eta_facet_wrong_1mT) > 1e-6,
   'sanity: the uncoated-Fresnel resolution actually differs from the old (wrong) '
   '1-T_coated resolution for this coating (not a vacuously-passing check)')

# item 5: loss-sensitive independent check -- an explicit finite sum of the
# first 200 round trips, with path lengths accumulated term by term (a
# separate derivation from the closed-form geometric series under test,
# not a restatement of it), at R_back=1 (worst-case, no truncation slack),
# alpha=20/cm, L=300um, T=0.6, dot_position=0.5 default. Forward-emitted
# escape after n round trips travels x*L + 2*n*L; backward-first escape
# (extra R_back factor) travels (2-x)*L + 2*n*L (dot -> back facet
# (1-x)*L, then a FULL pass L back to the front facet, per the
# facet_escape_fraction docstring's L/2+L worked example at x=0.5).
T_fs, Rb_fs, alpha_fs, L_fs, x_fs = 0.6, 1.0, 20.0, 300.0, 0.5
a_fs = alpha_fs * 1e-4
Rf_fs = 1.0 - T_fs
finite_sum = 0.0
for n in range(200):
    path_fwd = x_fs * L_fs + 2.0 * n * L_fs
    weight_fwd = 0.5 * T_fs * (Rf_fs * Rb_fs) ** n
    finite_sum += weight_fwd * exp(-a_fs * path_fwd)
    path_bwd = (2.0 - x_fs) * L_fs + 2.0 * n * L_fs
    weight_bwd = 0.5 * T_fs * Rb_fs * (Rf_fs * Rb_fs) ** n
    finite_sum += weight_bwd * exp(-a_fs * path_bwd)
ck(abs(finite_sum - facet_escape_fraction(T_fs, Rb_fs, alpha_fs, L_fs)) < 1e-10,
   'facet_escape_fraction matches an explicit 200-round-trip path-length sum '
   '(R_back=1, alpha=20/cm, L=300um, T=0.6)')

# item 5: the degenerate denominator (R_back*R_front*prop_rt == 1) must raise
# a clear ValueError rather than silently dividing by zero. alpha_cm=0 makes
# prop_rt=1, so T=0 (R_front=1) with R_back=1 forces the denominator to 0.
try:
    facet_escape_fraction(0.0, 1.0, 0.0, 100.0)
    ck(False, 'degenerate facet cavity (R_back*R_front*prop_rt==1) raises ValueError')
except ValueError:
    ck(True, 'degenerate facet cavity (R_back*R_front*prop_rt==1) raises ValueError')

# item 5: citation and textbook cross-check note are present in the module.
_wg_source = inspect.getsource(wgmod)
ck(_wg_source.count('Coldren, Corzine & Masanovic') >= 2
  and '2nd ed. (2012), ch. 2' in _wg_source
  and 'Coldren & Corzine' not in _wg_source,
   'module cites the three-author textbook edition consistently (old two-author form gone)')
ck('0.9563' in _wg_source and '0.9636' in _wg_source,
   'module keeps the textbook cross-check note (series 0.9563 vs F1 0.9636)')

print(f'{sum(checks)}/{len(checks)} waveguide checks passed')
sys.exit(0 if all(checks) else 1)
