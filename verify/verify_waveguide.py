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
print(f'{sum(checks)}/{len(checks)} waveguide checks passed')
sys.exit(0 if all(checks) else 1)
