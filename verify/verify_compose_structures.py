"""Checks for the composable-engine structure registry (fsim_core/compose/structures.py, P3).

Closed forms, published numbers and bit-for-bit comparisons against the modules the adapters wrap.
Run: python verify/verify_compose_structures.py   (exit code 0 iff all pass)
"""
from __future__ import annotations

import json
import math
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import numpy as np                                                       # noqa: E402

from fsim_core import dbr, nitride_cavity, waveguide                      # noqa: E402
from fsim_core import nitride_nanowire_photonics as nwp                   # noqa: E402
from fsim_core.spectral import cavity_transmission                        # noqa: E402
from fsim_core.compose import registry as R                               # noqa: E402
from fsim_core.compose import structures as S                             # noqa: E402
from fsim_core.compose.types import EmitterModel, StructureResponse       # noqa: E402

RESULTS = []
H_UEV_PER_GHZ = 4.135667696                # h in ueV/GHz (CODATA 2018, exact)


def check(name, ok, detail=""):
    RESULTS.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))


def make_em(**kw):
    d = dict(E_X_meV=1850.0, lambda_nm=670.0, n_host=3.4, gamma0_ns=1.0, k_nr_ns=0.0, gamma_star_ns=0.1,
             phonon=None, delta_xx_meV=3.0, gamma_xx_ratio=1.0, fss_ueV=0.0, transitions=frozenset({"X"}),
             dipole="x", dEdT_meV_K=-0.3, sigma_sd_ueV=0.0, tau_sd_ns=1.0)
    d.update(kw)
    return EmitterModel(**d)


def resp(kind, vals=None, em=None, T=4.2, NA=0.7):
    return R.get("structure", kind).response(vals or {}, em or make_em(), T, NA)


class V:
    """Accessor mirroring structures._V for calling the module's raw helpers."""
    def __init__(self, kind, vals):
        self.e, self.v = R.get("structure", kind), vals
    __call__ = lambda self, n: self.v.get(n, self.e.param(n).default)           # noqa: E731
    f = lambda self, n: float(self(n))                                            # noqa: E731
    pos = f
    frac = f


em0 = make_em()

# 1. F_P closed form identity (micropillar, measured Q) -- independent re-derivation of L_pen and V
lam, n_h, n_l, d, Q = 670.0, 3.5, 2.95, 2.5, 6124.0
Lpen = lam / (4 * (n_h - n_l)) * n_h * n_l / n_h ** 2             # first_layer='L', n_c = n_h
Leff = (lam / n_h + 2 * Lpen) * 1e-3
Vum = math.pi * d * d / 4 * Leff
F_ref = 3 / (4 * math.pi ** 2) * (lam * 1e-3 / n_h) ** 3 * Q / Vum
r = resp("micropillar", {"Q_measured": Q, "d_um": d, "n_high": n_h, "n_low": n_l})
check("F_P closed form (micropillar, Q 6124, d 2.5 um)", abs(r.F_P / F_ref - 1) < 1e-12, f"{r.F_P:.12g} vs {F_ref:.12g}")
r2 = resp("generic_cavity", {"F_P": 7.3, "Q": 5000.0})
check("generic_cavity returns F_P and kappa = E_mode/Q exactly",
      r2.F_P == 7.3 and abs(r2.kappa_meV - 1850.0 / 5000.0) < 1e-15)
rd = resp("generic_cavity", {"Q": 5000.0, "detuning_ueV": 200.0})
check("detuning convention: E_mode = E_X - detuning (positive = dot above mode)",
      abs(rd.E_mode_meV - (1850.0 - 0.2)) < 1e-12 and abs(rd.kappa_meV - rd.E_mode_meV / 5000.0) < 1e-15)
rg = resp("generic_cavity", {"F_P": 9.0, "kappa_ueV": 100.0})
check("g from F_P = 4 g^2/(kappa hbar Gamma0)",
      abs(S.purcell_from_g(rg.g_meV, rg.kappa_meV, 1.0) / 9.0 - 1) < 1e-12)

# 2. Tomm numbers (open_cavity direct)
em_t = make_em(gamma0_ns=2 * math.pi * 0.30)                      # gamma/2pi = 0.30 GHz [V]
g_ueV, k_ueV = 4.16 * H_UEV_PER_GHZ, 23.3 * H_UEV_PER_GHZ
rt = resp("open_cavity", {"g_ueV": g_ueV, "kappa_ueV": k_ueV}, em=em_t)
check("Tomm F_P = 4g^2/(kappa gamma) = 9.9 +- 0.05", abs(rt.F_P - 9.9) < 0.05, f"{rt.F_P:.3f}")
check("Tomm 1 + F_P = 10.9 vs measured 520/47.5 = 10.95 within 2%", abs((1 + rt.F_P) / (520 / 47.5) - 1) < 0.02,
      f"{1 + rt.F_P:.3f}")
check("open_cavity sets g_meV for the regime test", abs(rt.g_meV - g_ueV * 1e-3) < 1e-15 and rt.kappa_meV == k_ueV * 1e-3)
rgeo = resp("open_cavity", {"mode": "geometry"}, em=em_t)
check("open_cavity geometry mode finite, eta_out = T_top/(T_top + T_bot + loss)",
      math.isfinite(rgeo.F_P) and rgeo.F_P > 0 and abs(rgeo.eta_out - 1e-3 / (1e-3 + 5e-5 + 5e-5)) < 1e-12)

# 3. Somaschi beta
rs = resp("generic_cavity", {"F_P": 7.6})
beta = (rs.beta0 * rs.F_bg + rs.F_P) / (rs.F_bg + rs.F_P)          # rates.combine formula at L = 1
check("Somaschi F_P 7.6, leaky = Gamma0 -> beta 0.884 (paper 0.88)", abs(beta - 0.884) < 0.01, f"{beta:.4f}")
rl = resp("cbg", {"F_P": 4.0, "leak": 1.0})
check("cbg leak: beta = F/(1 + F + leak)",
      abs((rl.beta0 * rl.F_bg + rl.F_P) / (rl.F_bg + rl.F_P) - 4.0 / (1 + 4.0 + 1.0)) < 1e-15)

# 4. Lorentzian detuning
kap, gh = 100.0, 1e-6
check("L(0) = 1 for a narrow line", abs(S.lorentz_L(0.0, gh, kap) - 1.0) < 1e-6)
check("L(kappa/2) = 0.5 for gamma_h << kappa", abs(S.lorentz_L(kap / 2, gh, kap) - 0.5) < 1e-6)
x = np.linspace(-300, 300, 41)
check("lorentz_L == spectral.cavity_transmission bit-for-bit",
      np.array_equal(S.lorentz_L(x, 12.0, kap), cavity_transmission(x, 12.0, kap)))
check("L(0) = kappa/(gamma_h + kappa) for a broad line", abs(S.lorentz_L(0.0, 50.0, 100.0) - 100.0 / 150.0) < 1e-15)

# 5. bulk_planar number
rb = resp("bulk_planar", {"n": 3.5}, NA=0.7)
ref5 = 3 / 8 * 0.2 ** 2 * 4 * 3.5 / 4.5 ** 2
check("bulk_planar eta(n 3.5, NA 0.7) = 0.0104 +-1e-4", abs(rb.eta_out - ref5) < 1e-12 and abs(rb.eta_out - 0.0104) < 1e-4,
      f"{rb.eta_out:.5f}")
rm = resp("bulk_planar", {"n": 3.5, "mirror": True}, NA=0.7)
check("bulk_planar mirror doubles eta", abs(rm.eta_out - 2 * rb.eta_out) < 1e-15)
rh = resp("bulk_planar", {"n": 3.5, "sil": "hemi", "n_s": 2.0}, NA=0.2)
rw = resp("bulk_planar", {"n": 3.5, "sil": "weierstrass", "n_s": 2.0}, NA=0.1)
b02 = S.bulk_planar_eta(0.2, 3.5)
b01 = S.bulk_planar_eta(0.1, 3.5)
check("SIL gains: hemi x n_s^2, Weierstrass x n_s^4 (below the cap)",
      abs(rh.eta_out / b02 - 4.0) < 1e-12 and abs(rw.eta_out / b01 - 16.0) < 1e-12)
rc = resp("bulk_planar", {"n": 3.5, "sil": "weierstrass", "n_s": 3.5}, NA=0.9)
check("SIL cap NA_eff/n <= 1 keeps eta bounded", rc.eta_out == 0.375 * 4 * 3.5 / 4.5 ** 2 and rc.eta_out < 1.0)
check("bulk_planar eta_out_fn(NA) == eta_out", rb.eta_out_fn(0.7) == rb.eta_out)
check("bulk_planar: no cavity (F_P 0, beta0 1)", rb.F_P == 0.0 and rb.beta0 == 1.0 and rb.kappa_meV is None)

# 6. T_NA
na = np.linspace(0.01, 1.0, 60)
tna = S.T_NA(na, 930.0, 1.0)
check("T_NA monotone in NA", bool(np.all(np.diff(tna) > 0)) and S.T_NA(0.0, 930.0, 1.0) == 0.0)
check("T_NA -> 1 at large NA", S.T_NA(1.0, 930.0, 5.0) > 1 - 1e-12)
theta = 0.93 / (math.pi * 1.0)
check("T_NA closed form 1 - exp(-2 NA^2/theta^2)", abs(S.T_NA(0.3, 930.0, 1.0) - (1 - math.exp(-2 * 0.09 / theta ** 2))) < 1e-15)
rp = resp("micropillar", {}, NA=0.7)
check("micropillar eta_out_fn(NA) == eta_out and rises with NA", rp.eta_out_fn(0.7) == rp.eta_out
      and rp.eta_out_fn(0.9) > rp.eta_out_fn(0.4))

# 7. Adapters bit-for-bit
# edge_ridge vs waveguide.edge_emission (cards/edge-inp-gainp-design.yaml: ridge 2000, NA 0.75, R_back 0.95, L 250, alpha 5, 770 nm)
em_e = make_em(lambda_nm=770.0)
re = resp("edge_ridge", {"ridge_width_nm": 2000.0, "L_um": 250.0, "R_back": 0.95, "alpha_cm": 5.0}, em=em_e, NA=0.75)
e = waveguide.edge_emission(waveguide.hkust_ridge_stack(770.0), 2000.0, 1200.0, 770.0, 250.0, 0.75,
                            alpha_cm=5.0, R_back=0.95, coating=None, na_method="gaussian", dot_position=0.5)
check("edge_ridge F_P == edge_emission.F_wg bit-for-bit", re.F_P == e.F_wg)
check("edge_ridge beta and total reproduce edge_emission (beta0 = 0)",
      (re.F_P / (1 + re.F_P)) == e.beta and abs(re.eta_out * (re.F_P / (1 + re.F_P)) / e.eta_total - 1) < 1e-14)
rb_none = resp("edge_ridge", {"R_back": -1.0}, em=em_e, NA=0.5)
e2 = waveguide.edge_emission(waveguide.hkust_ridge_stack(770.0), 2000.0, 1200.0, 770.0, 250.0, 0.5, alpha_cm=5.0, R_back=None)
check("edge_ridge R_back < 0 -> R_back=None (bare facet)", rb_none.F_P == e2.F_wg and abs(rb_none.eta_out * e2.beta / e2.eta_total - 1) < 1e-14)

# nitride_planar_cavity vs nitride_cavity.response (cards/nitride-cavity-pulse-design.yaml: Q 2000, V 2.0, overlap 1, eta_out 0.1, T_track 300, slope -0.04)
em_n = make_em(E_X_meV=2700.0, lambda_nm=459.0, gamma0_ns=1.0, k_nr_ns=0.0, gamma_star_ns=0.0, delta_xx_meV=-10.0)
rn = resp("nitride_planar_cavity", {}, em=em_n, T=300.0)
prm = nitride_cavity.NitrideCavityParams(Q=2000.0, mode_volume_norm=2.0, spatial_overlap=1.0, eta_out=0.1,
                                         T_track=300.0, dEdT_cav_meV_K=-0.04, detuning_offset_meV=0.0, purcell_enabled=True)
gam = S.HBAR_MEV_NS * 1.0
rr = nitride_cavity.response(prm, T_K=300.0, E_X_eV=2.7, E_X_track_eV=2.7, gamma_X_meV=gam, gamma_XX_meV=gam,
                             delta_xx_meV=-10.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
check("nitride_planar_cavity F_P, kappa, E_mode, eta_out == nitride_cavity.response bit-for-bit",
      rn.F_P == rr["Fp_add"] and rn.kappa_meV == rr["kappa_meV"] and rn.E_mode_meV == 1000.0 * rr["E_cav_eV"]
      and rn.eta_out == rr["eta_out"])
check("nitride F_P = (3/(4 pi^2)) Q/V closed form", abs(rn.F_P - 3 / (4 * math.pi ** 2) * 2000.0 / 2.0) < 1e-12)
rn2 = resp("nitride_planar_cavity", {"purcell_enabled": False}, em=em_n, T=300.0)
check("nitride purcell_enabled False -> F_P 0 (cavity filtering kept)", rn2.F_P == 0.0 and rn2.kappa_meV == rn.kappa_meV)

# planar_dbr vs dbr functions
vp = V("planar_dbr", {})
rpd = resp("planar_dbr", {}, em=em0)
top, spacer, bottom = dbr.cavity_stack(3.5, 2.95, 3.5, 20, 30, 670.0)
mode = dbr.cavity_mode(top + [spacer] + bottom, 670.0, n_in=1.0, n_out=3.5)
dc = 670.0 / 3.5
Ft = dbr.planar_total_rate(top, bottom, 3.5, dc, 0.5 * dc, 670.0, n_in=1.0, n_out=3.5, orientation="inplane")
check("planar_dbr E_mode, kappa == dbr.cavity_mode bit-for-bit", rpd.E_mode_meV == mode["E_cav_meV"] and rpd.kappa_meV == mode["kappa_meV"])
check("planar_dbr F_bg == dbr.planar_total_rate bit-for-bit and F_P == 0 (never the 1-D LDOS)",
      rpd.F_bg == Ft and rpd.F_P == 0.0)
Tt = 1 - dbr.qw_peak_reflectivity(1.0, 3.5, 3.5, 2.95, 20)
Tb = 1 - dbr.qw_peak_reflectivity(3.5, 3.5, 3.5, 2.95, 30)
check("planar_dbr eta_out = T_top/(T_top + T_bot) from the stack", rpd.eta_out == Tt / (Tt + Tb))
F1d = dbr.planar_purcell(top, bottom, 3.5, dc, 0.5 * dc, 670.0, n_in=1.0, n_out=3.5)
check("planar_dbr notes quote the 1-D LDOS but F_P stays 0", any(f"{F1d:.3g}" in n for n in rpd.notes) and rpd.F_P == 0.0)

# photonic_nanowire vs nitride_nanowire_photonics.response
rw_ = resp("photonic_nanowire", {"d_nm": 220.0, "n_wire": 3.45, "n_group": 3.6, "taper_T": 0.95, "mirror_R": 0.9}, em=em0, NA=0.75)
pm = nwp.NitrideNanowirePhotonicsParams(family="vertical_photonic", NA=0.75, n_wire=3.45, n_group_override=3.6,
                                         taper_transmission=0.95, bottom_reflectivity=0.9, radiative_rate_factor=1.0)
nr = nwp.response(pm, lambda_nm=670.0, outer_radius_nm=110.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
check("photonic_nanowire beta0 == beta_HE11 bit-for-bit", rw_.beta0 == nr["beta_HE11"] and rw_.F_P == 0.0)
check("photonic_nanowire beta0 x eta_out reproduces eta_collection_X", abs(rw_.beta0 * rw_.eta_out / nr["eta_collection_X"] - 1) < 1e-13)
check("photonic_nanowire eta_out_fn(NA) matches the module at another NA",
      abs(rw_.eta_out_fn(0.4) - nwp.response(nwp.NitrideNanowirePhotonicsParams(
          family="vertical_photonic", NA=0.4, n_wire=3.45, n_group_override=3.6, taper_transmission=0.95,
          bottom_reflectivity=0.9), lambda_nm=670.0, outer_radius_nm=110.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)["eta_collection_X"]
          / nr["beta_HE11"]) < 1e-12)
rwi = resp("photonic_nanowire", {"F_bg": 0.4}, em=em0)
check("photonic_nanowire F_bg < 1 inhibition is passed through", rwi.F_bg == 0.4)

# other physics identities
fw = 3 / (4 * math.pi) * (0.67 / 3.4) ** 2 * (5.0 / 3.4) / 0.03
rpw = resp("phc_waveguide", {"n_g": 5.0, "A_eff_um2": 0.03}, em=em0)
check("phc_waveguide F_wg = (3/(4 pi)) (lambda/n)^2 (n_g/n)/A_eff", abs(rpw.F_P / fw - 1) < 1e-12 and rpw.kappa_meV is None)
rpb = resp("phc_waveguide", {"beta_direct": 0.92, "F_bg": 0.05}, em=em0)
check("phc_waveguide beta_direct: F_wg/(F_wg + F_bg) = 0.92", abs(rpb.F_P / (rpb.F_P + rpb.F_bg) - 0.92) < 1e-12)
rpg = resp("phc_cavity", {"g_ueV": 100.0, "Q": 1e4}, em=em0)
check("phc_cavity g entered: F_P = 4g^2/(kappa hbar Gamma0)", abs(rpg.F_P / S.purcell_from_g(0.1, rpg.kappa_meV, 1.0) - 1) < 1e-12)
rpq = resp("phc_cavity", {"Q": 1e4, "V_lambda3": 0.7}, em=em0)
check("phc_cavity F_P = 3/(4 pi^2) Q/V", abs(rpq.F_P / (3 / (4 * math.pi ** 2) * 1e4 / 0.7) - 1) < 1e-12)
fn = S.eta_table_fn("0.5:0.2, 0.9:0.6")
check("eta table: interpolation, flat above, cone scaling below",
      abs(fn(0.7) - 0.4) < 1e-12 and fn(1.0) == 0.6 and abs(fn(0.25) - 0.05) < 1e-12)
rcb = resp("cbg", {"eta_table": "0.5:0.2, 0.9:0.6"}, NA=0.7)
check("cbg eta_out follows the table at NA", abs(rcb.eta_out - 0.4) < 1e-12 and rcb.eta_out_fn(0.9) == 0.6)
mp = resp("micropillar", {"eta_out": 0.5}, NA=0.7)
check("micropillar eta_out override replaces the mirror asymmetry (T_NA still applied)",
      abs(mp.eta_out - 0.5 * S.T_NA(0.7, 670.0, 2.5 / 2.7)) < 1e-15)
mq = resp("micropillar", {"Q_measured": 3000.0})
check("micropillar Q_measured overrides computed Q (kappa = E/Q)", abs(mq.kappa_meV - 1850.0 / 3000.0) < 1e-12)

# 8. every entry finite at defaults
kinds = R.kinds("structure")
want = {"generic_cavity", "micropillar", "open_cavity", "cbg", "phc_cavity", "phc_waveguide", "photonic_nanowire",
        "planar_dbr", "bulk_planar", "sil_mesa", "mesa_reflector", "edge_ridge", "nitride_planar_cavity"}
check("all 13 structure kinds registered", want <= set(kinds), f"missing {sorted(want - set(kinds))}")
bad = []
for k in sorted(want):
    try:
        rr_ = resp(k, {}, em=em0, T=300.0, NA=0.65)
        fin = all(math.isfinite(x) for x in (rr_.F_P, rr_.F_bg, rr_.beta0, rr_.eta_out, rr_.pol_factor))
        fin &= all(x is None or math.isfinite(x) for x in (rr_.kappa_meV, rr_.E_mode_meV, rr_.g_meV))
        fin &= bool(isinstance(rr_, StructureResponse) and 0.0 <= rr_.eta_out <= 1.0 and rr_.F_P >= 0 and isinstance(rr_.inputs, dict) and rr_.notes)
        if not fin:
            bad.append(k)
    except Exception as ex:                                              # noqa: BLE001
        bad.append(f"{k}: {ex!r}")
check("response() finite at defaults for every entry", not bad, "; ".join(bad))

# 9. registry payload
pl = R.registry_payload()
try:
    json.dumps(pl)
    ser = True
except TypeError:
    ser = False
check("registry_payload JSON-serialisable", ser)
miss = []
for s in pl["structures"]:
    for p in s["params"]:
        if not (p["unit"] is not None and p["tag"] in ("V", "DR", "E", "A") and p["help"] and p["source"]):
            miss.append(f"{s['kind']}.{p['name']}")
check("every ParamSpec has unit / tag / help / source", not miss, ", ".join(miss))
from fsim_core.compose.params import resolve_leaf                            # noqa: E402
errs = []
for k in want:
    for p in R.get("structure", k).params:
        try:
            resolve_leaf(p.default, p, f"structure.{p.name}")
        except ValueError as ex:
            errs.append(f"{k}.{p.name}: {ex}")
check("all ParamSpec defaults resolve as leaves", not errs, "; ".join(errs))

n_ok = sum(RESULTS)
print(f"{n_ok}/{len(RESULTS)} compose_structures checks passed")
sys.exit(0 if n_ok == len(RESULTS) else 1)
