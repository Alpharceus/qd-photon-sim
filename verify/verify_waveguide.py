"""Independent numerical checks for fsim_core.waveguide (published-class values)."""
import sys
from pathlib import Path
from math import pi, tan
import numpy as np
from scipy.optimize import brentq
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fsim_core.waveguide import (Layer, slab_modes, effective_index_ridge, beta_factor,
                                 facet_transmission, na_collection, hkust_ridge_stack,
                                 edge_emission)

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
s=hkust_ridge_stack()
rm=effective_index_ridge(s,668,2000,1200)
g=rm.vertical.gamma_layer('dot')
ck(3.05<rm.n_eff<3.22, 'HKUST n'); ck(.005<g<.05, 'HKUST gamma')
r=edge_emission(s,2000,1200,668,500,.5)
ck(0<r.eta_total<1 and r.n_g>0, 'edge result')

# council review 2026-09-05 item 3: T_facet must be multiplied into eta_total
# exactly once (front stays the pure geometric 0.5 split, no R_back given).
expected_total = r.beta * 0.5 * r.T_facet * r.eta_prop * r.eta_NA
ck(abs(expected_total - r.eta_total) < 1e-9, 'facet transmission included in eta_total')
ck(r.T_facet < 1.0, 'facet transmission is not unity (so omitting it was not a no-op)')

# council review 2026-09-06 item 1: two facet models, chosen by whether
# R_back is given, and T_facet must appear exactly once in eta_total either
# way (the previous bug applied it a second time on top of the escape-rate
# fraction, ~4% flux under-report at the class R_back=0.95 point).
r_back0 = edge_emission(s, 2000, 1200, 668, 500, .5, R_back=0.0)
ck(abs(r_back0.eta_total - r.eta_total) < 1e-9,
   'R_back=0.0 reproduces the uncoated (R_back=None) eta_total exactly')
r_back95 = edge_emission(s, 2000, 1200, 668, 500, .5, R_back=0.95)
expected_escape = r_back95.beta * (r_back95.T_facet / (r_back95.T_facet + (1 - 0.95))) \
    * r_back95.eta_prop * r_back95.eta_NA
ck(abs(expected_escape - r_back95.eta_total) < 1e-9,
   'R_back=0.95 gives the escape-rate facet value exactly (no extra factor of T_facet)')
r_backs = [edge_emission(s, 2000, 1200, 668, 500, .5, R_back=rb).eta_total
           for rb in (0.0, 0.3, 0.6, 0.9, 0.95, 0.99)]
ck(all(a <= b + 1e-12 for a, b in zip(r_backs, r_backs[1:])),
   'eta_total is monotone non-decreasing in R_back (escape-rate model)')

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

print(f'{sum(checks)}/{len(checks)} waveguide checks passed')
sys.exit(0 if all(checks) else 1)
