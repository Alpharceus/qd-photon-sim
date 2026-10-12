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

# 1. micropillar F_P: arithmetic regression of V = (pi/4) w0^2 L_eff (w0 = d/2.7), PLUS the independent checks below
lam, n_h, n_l, d, Q = 670.0, 3.5, 2.95, 2.5, 6124.0
Lpen = lam / (4 * (n_h - n_l)) * n_h * n_l / n_h ** 2             # first_layer='L', n_c = n_h
Leff = (lam / n_h + 2 * Lpen) * 1e-3
Vum = math.pi / 4 * (d / 2.7) ** 2 * Leff
F_ref = 3 / (4 * math.pi ** 2) * (lam * 1e-3 / n_h) ** 3 * Q / Vum
r = resp("micropillar", {"Q_measured": Q, "d_um": d, "n_high": n_h, "n_low": n_l})
check("F_P arithmetic (micropillar, Q 6124, d 2.5 um): V = (pi/4) w0^2 L_eff, w0 = d/2.7",
      abs(r.F_P / F_ref - 1) < 1e-12, f"{r.F_P:.12g} vs {F_ref:.12g}")

# 1b. review finding 1/14: the IDEAL F_P must bound the measured one from above. Ding et al., PRL 116, 020401 (2016):
# d 2.5 um, Q 6124, lambda 897 nm, GaAs/AlAs; measured F_P = 6.3 [V] (= ideal x r x cos^2 x spectral overlap < ideal).
em_d = make_em(lambda_nm=897.0, E_X_meV=1239.84198 / 0.897)
rdg = resp("micropillar", {"Q_measured": 6124.0, "d_um": 2.5}, em=em_d)
check("Ding ideal F_P (d 2.5, Q 6124, 897 nm, GaAs/AlAs) >= measured 6.3 (ideal bounds measured from above)",
      rdg.F_P >= 6.3, f"ideal F_P = {rdg.F_P:.3f}")
rdc = resp("micropillar", {"Q_measured": 6124.0, "d_um": 2.5, "mode_area": "cylinder"}, em=em_d)
check("cylinder option is an upper bound on V: F_P(cylinder) < F_P(gaussian)", rdc.F_P < rdg.F_P,
      f"{rdc.F_P:.3f} < {rdg.F_P:.3f}")
rdp = resp("micropillar", {"Q_measured": 6124.0, "d_um": 2.5, "pol_factor": 0.5, "r_factor": 0.8}, em=em_d)
check("pol_factor and r_factor are folded into F_P (response pol_factor = 1)",
      abs(rdp.F_P / rdg.F_P - 0.4) < 1e-12 and rdp.pol_factor == 1.0)


def tmm_half_length_nm(lam_nm, nh, nl, nc, Nt, Nb):
    """Independent TMM (own code, not fsim_core.dbr): int n^2|E|^2 dz/(nc^2 |E_max|^2) over both mirrors plus the
    cavity, with the field started at the spacer antinode (E = 1, H = 0) and propagated outward through each side:
    half spacer lam/(2 nc), then Nt/Nb x (low lam/(4 nl), high lam/(4 nh)) [the order the repo's micropillar uses]."""
    def side(N):
        layers = [(nc, lam_nm / (2 * nc))] + [(nl, lam_nm / (4 * nl)), (nh, lam_nm / (4 * nh))] * N
        E, H, tot = 1.0 + 0j, 0j, 0.0
        for n, dz in layers:
            k = 2 * math.pi * n / lam_nm
            zs = np.linspace(0.0, dz, 400)
            Ez = E * np.cos(k * zs) + 1j * H / n * np.sin(k * zs)
            tot += np.trapezoid(n * n * np.abs(Ez) ** 2, zs)
            E, H = E * math.cos(k * dz) + 1j * H / n * math.sin(k * dz), 1j * n * E * math.sin(k * dz) + H * math.cos(k * dz)
        return tot / (nc * nc)
    return side(Nt) + side(Nb)


cd = S._micropillar_calc(V("micropillar", {"d_um": 2.5}), em_d)
Lz = tmm_half_length_nm(897.0, 3.5, 2.95, 3.5, 25, 40)
check("TMM identity: int n^2|E|^2 dz / (n_c^2 |E_max|^2) = L_eff/2 within 1% (standing wave halves the energy length)",
      abs(Lz / (0.5 * cd["L_eff_um"] * 1e3) - 1) < 0.01, f"{Lz:.1f} nm vs L_eff/2 = {0.5 * cd['L_eff_um'] * 1e3:.1f} nm")
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
rl = resp("cbg", {"F_P": 4.0, "leak": 1.0, "table_mode": "cavity_only"})
check("cbg leak (table_mode cavity_only): beta = F/(1 + F + leak)",
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
rw_ = resp("photonic_nanowire", {"d_nm": 220.0, "n_wire": 3.45, "n_group": 3.6, "taper_T": 0.95, "mirror_R": 0.9,
                                 "mirror_phase_ok": True}, em=em0, NA=0.75)      # FX4: the antinode case is explicit
pm = nwp.NitrideNanowirePhotonicsParams(family="vertical_photonic", NA=0.75, n_wire=3.45, n_group_override=3.6,
                                         taper_transmission=0.95, bottom_reflectivity=0.9, radiative_rate_factor=1.0)
nr = nwp.response(pm, lambda_nm=670.0, outer_radius_nm=110.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
pm0 = nwp.NitrideNanowirePhotonicsParams(family="vertical_photonic", NA=0.75, n_wire=3.45, n_group_override=3.6,
                                          taper_transmission=0.95, bottom_reflectivity=0.0, radiative_rate_factor=1.0)
nr0 = nwp.response(pm0, lambda_nm=670.0, outer_radius_nm=110.0, gamma_X0_ns=1.0, gamma_XX0_ns=1.0)
check("photonic_nanowire eta_out == taper_T x module HE11 acceptance bit-for-bit (F_P 0)",
      rw_.eta_out == 0.95 * nr["diagnostics"]["obj_accept"] and rw_.F_P == 0.0)
check("photonic_nanowire eta_out_fn(NA) matches the module acceptance at another NA",
      abs(rw_.eta_out_fn(0.4) - 0.95 * nwp.response(nwp.NitrideNanowirePhotonicsParams(
          family="vertical_photonic", NA=0.4, n_wire=3.45, n_group_override=3.6, taper_transmission=0.95,
          bottom_reflectivity=0.0), lambda_nm=670.0, outer_radius_nm=110.0, gamma_X0_ns=1.0,
          gamma_XX0_ns=1.0)["diagnostics"]["obj_accept"]) < 1e-12)
rw0 = resp("photonic_nanowire", {"d_nm": 220.0, "n_wire": 3.45, "n_group": 3.6, "taper_T": 0.95, "mirror_R": 0.0}, em=em0, NA=0.75)
beta_w = 2.0 * rw0.beta0                                            # R = 0: upward fraction = beta/2
rr_ = math.sqrt(0.9)
check("photonic_nanowire mirror (Munsch Eq. 1 closed form): beta0 = (beta/2)(1 + |r|)^2/(1 + beta|r|), rate x (1 + beta|r|)",
      abs(rw_.beta0 - 0.5 * beta_w * (1 + rr_) ** 2 / (1 + beta_w * rr_)) < 1e-13
      and abs(rw_.F_bg - (1 + beta_w * rr_)) < 1e-13 and rw0.F_bg == 1.0)
rwi = resp("photonic_nanowire", {"F_bg": 0.4, "mirror_R": 0.0}, em=em0)
check("photonic_nanowire F_bg < 1 inhibition is passed through (no mirror)", rwi.F_bg == 0.4)
rwp = resp("photonic_nanowire", {"d_nm": 220.0, "n_wire": 3.45, "n_group": 3.6, "taper_T": 0.95, "mirror_R": 0.9,
                                 "mirror_phase_ok": False}, em=em0, NA=0.75)
check("photonic_nanowire mirror_phase_ok false: incoherent (beta/2)(1 + R), rate unchanged",
      abs(rwp.beta0 - 0.5 * beta_w * 1.9) < 1e-13 and rwp.F_bg == 1.0)

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
rcb = resp("cbg", {"eta_table": "0.5:0.2, 0.9:0.6", "table_mode": "cavity_only"}, NA=0.7)
bres = 10.0 / 11.0                                                   # beta at resonance, F_P 10, leak 0, beta0 0
check("cbg cavity_only (legacy): eta_out = table/beta_res and beta x eta_out = table at L = 1",
      abs(rcb.eta_out - 0.4 / bres) < 1e-12 and abs(rcb.eta_out_fn(0.9) - 0.6 / bres) < 1e-12
      and abs(bres * rcb.eta_out - 0.4) < 1e-12)
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

# 10. devlib-fx2a: cbg r_factor, nanowire mirror / taper expansion, phc beta_direct
c1 = resp("cbg", {"F_P": 12.0})
c5 = resp("cbg", {"F_P": 12.0, "r_factor": 0.5})
check("cbg r_factor is folded into F_P (F_P x r_factor; default 1 is bit-identical to the old model)",
      c5.F_P == 6.0 and c1.F_P == 12.0 and resp("cbg", {"F_P": 12.0, "r_factor": 1.0}).F_P == c1.F_P)
b0, b5 = c1.F_P / (c1.F_bg + c1.F_P), c5.F_P / (c5.F_bg + c5.F_P)      # F_bg = 1 + leak
check("cbg r_factor lowers beta consistently: beta = F_P r/(1 + leak + F_P r) at leak 0",
      abs(b5 - 6.0 / 7.0) < 1e-12 and abs(b0 - 12.0 / 13.0) < 1e-12)
emn = make_em(E_X_meV=1409.0, lambda_nm=879.9, n_host=3.37, gamma0_ns=0.6)
nw = {"d_nm": 200.0, "n_wire": 3.37}


def up_extraction(**kw):
    rr = resp("photonic_nanowire", {**nw, "mirror_phase_ok": True, **kw}, em=emn, NA=1.0)   # NA 1: HE11 acceptance is exactly 1
    return rr.beta0 * rr.eta_out


x0, x1, x3 = up_extraction(mirror_R=0.0), up_extraction(mirror_R=1.0), up_extraction(mirror_R=0.3)
bw = 2.0 * x0                                                            # NA 1, taper 0.95: x0 = (beta/2) T
bw /= 0.95
check("nanowire mirror (FX3a, coherent): extraction(R)/extraction(0) = (1 + |r|)^2/(1 + beta|r|)  (R 1 and R 0.3)",
      abs(x1 / x0 - 4.0 / (1 + bw)) < 1e-9 and abs(x3 / x0 - (1 + math.sqrt(0.3)) ** 2 / (1 + bw * math.sqrt(0.3))) < 1e-9,
      f"{x1 / x0:.6f} {x3 / x0:.6f}")
rn = resp("photonic_nanowire", {**nw, "mirror_R": 0.3, "mirror_phase_ok": True}, em=emn, NA=0.75)
rt_ = resp("photonic_nanowire", {**nw, "mirror_R": 0.3, "mirror_phase_ok": True, "taper_mfr_nm": 1000.0}, em=emn, NA=0.75)
bare = nwp.response(nwp.NitrideNanowirePhotonicsParams(
    family="vertical_photonic", NA=0.75, n_wire=3.37, taper_transmission=0.95, bottom_reflectivity=0.3,
    radiative_rate_factor=1.0, n_group_override=3.6), lambda_nm=879.9, outer_radius_nm=100.0,
    gamma_X0_ns=0.6, gamma_XX0_ns=0.6)
check("nanowire taper_mfr_nm = 0 (default): eta_out = taper_T x the bare-facet nitride_nanowire_photonics acceptance bit-for-bit",
      rn.eta_out == 0.95 * bare["diagnostics"]["obj_accept"])
x_exp = rt_.beta0 * rt_.eta_out
x_ceil = up_extraction(mirror_R=0.3)                                     # the NA-1 ceiling (acceptance 1)
check("nanowire taper_mfr_nm: a large top-facet mode radius narrows the far field -> NA-0.75 extraction approaches the "
      "NA-1 ceiling from below and exceeds the bare-facet value",
      rn.beta0 * rn.eta_out < x_exp <= x_ceil + 1e-12 and x_exp > 0.99 * x_ceil,
      f"bare {rn.beta0 * rn.eta_out:.4f} < {x_exp:.4f} <= ceiling {x_ceil:.4f}")
from fsim_core.compose.structures import PHC_WAVEGUIDE                       # noqa: E402
from fsim_core.compose import observables as OBS                             # noqa: E402
rb = resp("phc_waveguide", {"beta_direct": 0.9})
check("phc_waveguide beta_direct 0.9: F_wg = beta/(1 - beta) F_bg and observable beta_entered is a passthrough of that leaf",
      abs(rb.F_P - 9.0) < 1e-12 and OBS.get_observable("beta_entered").passthrough_of == "structure.beta_direct"
      and PHC_WAVEGUIDE.param("beta_direct").tag in ("V", "DR", "E", "A"))

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

# 11. devlib-fx3a: cbg total-efficiency table, nanowire coherent mirror + HE11 beta
def first_lens(r, Lval=1.0):
    return (r.beta0 * r.F_bg + r.F_P * Lval) / (r.F_bg + r.F_P * Lval) * r.eta_out


tbl = {"eta_table": "0.7:0.26"}
c_lo, c_lo_old = resp("cbg", {**tbl, "F_P": 0.1}, NA=0.7), resp("cbg", {**tbl, "F_P": 0.1, "table_mode": "cavity_only"}, NA=0.7)
check("cbg total mode, low F_P 0.1 with table 0.26: first lens = 0.26 (legacy cavity_only gave F/(1 + F) = 0.091)",
      abs(first_lens(c_lo) - 0.26) < 1e-12 and abs(first_lens(c_lo_old) - 0.1 / 1.1) < 1e-12,
      f"{first_lens(c_lo):.4f} vs {first_lens(c_lo_old):.4f}")
c_hi, c_hi_old = resp("cbg", {**tbl, "F_P": 20.0}, NA=0.7), resp("cbg", {**tbl, "F_P": 20.0, "table_mode": "cavity_only"}, NA=0.7)
check("cbg high F_P 20: total and cavity_only agree at L = 1 (the table is the device number there)",
      abs(first_lens(c_hi) - 0.26) < 1e-12 and abs(first_lens(c_hi_old) - 0.26) < 1e-12)
Ld = float(S.lorentz_L(3.0, 0.01, 4.0))
check("cbg total mode: detuning (L < 1) changes the rate F_P L only, not the collected fraction beta x eta_out = 0.26 x 1",
      abs(first_lens(c_hi, Ld) - 0.26) < 1e-12 and first_lens(c_hi_old, Ld) < 0.26 - 1e-3 and Ld < 0.9)
check("cbg total mode: F_P and r_factor set only the rate (F_P r_factor still folded, F_bg = 1 + leak)",
      resp("cbg", {**tbl, "F_P": 8.0, "r_factor": 0.5, "leak": 0.5}).F_P == 4.0
      and resp("cbg", {**tbl, "F_P": 8.0, "leak": 0.5}).F_bg == 1.5 and resp("cbg", {**tbl}).beta0 == 1.0)
check("cbg table_mode default is total", S.CBG.param("table_mode").default == "total")

# nanowire: Munsch 2013 Eq. 1 numbers, quasi-static limit, Bleuse peak
check("nanowire_mirror reproduces Munsch Eq. 1: beta 0.96, |r| 0.95, T_Gauss 0.93 -> 0.89 [V]",
      abs(S.nanowire_mirror(0.96, 0.95 ** 2, True)[0] * 0.93 - 0.89) < 0.005, f"{S.nanowire_mirror(0.96, 0.9025, True)[0] * 0.93:.4f}")
check("nanowire_mirror limits: R 0 -> (beta/2, rate 1); |r| 1, beta 1 -> upward fraction 1 (all light up), rate 2",
      S.nanowire_mirror(0.8, 0.0, True) == (0.4, 1.0) and abs(S.nanowire_mirror(1.0, 1.0, True)[0] - 1.0) < 1e-15
      and abs(S.nanowire_mirror(1.0, 1.0, True)[1] - 2.0) < 1e-15)
sq = (2 / (3.45 ** 2 + 1)) ** 2
check("nanowire leaky screening: kappa -> (2/(n^2 + 1))^2 as d -> 0 and -> 1 for d >> lambda/n (closed forms)",
      abs(S.nanowire_leak_inhibition(3.45, 1e-4) - sq) < 1e-6 and abs(S.nanowire_leak_inhibition(3.45, 5.0) - 1.0) < 1e-12)


def beta_wire(dl, lam=930.0, **kw):
    e = make_em(lambda_nm=lam, E_X_meV=1239.84198 / lam * 1e3, n_host=3.45)
    rr = resp("photonic_nanowire", {"d_nm": dl * lam, "n_wire": 3.45, "n_group": 3.6, "taper_T": 1.0, "mirror_R": 0.0, **kw},
              em=e, NA=1.0)
    return 2.0 * rr.beta0


dls = np.linspace(0.15, 0.40, 51)
bcurve = np.array([beta_wire(x) for x in dls])
ipk = int(np.argmax(bcurve))
check("nanowire HE11 beta peak vs Bleuse PRL 106, 103601 (2011) [E]: ~0.95 within 0.02 at d/lambda 0.21-0.26 (GaAs n 3.45)",
      abs(bcurve[ipk] - 0.95) < 0.02 and 0.21 <= dls[ipk] <= 0.26, f"{bcurve[ipk]:.4f} at d/lambda {dls[ipk]:.3f}")
check("nanowire HE11 beta at d/lambda 0.40 ~0.7 within 0.05 (Bleuse 2011 / Claudon 2010 multimode decline) [E]",
      abs(beta_wire(0.40) - 0.70) < 0.05, f"{beta_wire(0.40):.4f}")
check("nanowire HE11 beta -> 0 for a vanishing wire and is below the peak for d/lambda 0.10",
      beta_wire(0.02) < 0.05 and beta_wire(0.10) < bcurve[ipk])
bm = beta_wire(240.0 / 950.0, lam=950.0)
check("nanowire HE11 beta at the Munsch 2013 design point (d 240 nm, lambda 950 nm) vs the paper's 0.96 [V]: within 0.05 "
      "(not tuned; the plateau beta > 0.9)", abs(bm - 0.96) < 0.05, f"{bm:.4f}")

# ---- FX4 (final review 7, 8, 10)
rwd = resp("photonic_nanowire", {"d_nm": 220.0, "n_wire": 3.45, "n_group": 3.6, "taper_T": 0.95, "mirror_R": 0.9}, em=em0, NA=0.75)
check("FX4: photonic_nanowire mirror_phase_ok defaults to False (incoherent power return, rate unchanged, [A])",
      R.get("structure", "photonic_nanowire").param("mirror_phase_ok").default is False and rwd.F_bg == 1.0
      and abs(rwd.beta0 - 0.5 * 2.0 * rw0.beta0 * 1.9) < 1e-13)
emo = make_em(E_X_meV=1450.0, lambda_nm=855.0, n_host=1.0, gamma0_ns=1.0)
ro1 = resp("open_cavity", {}, em=emo)
ro2 = resp("open_cavity", {"pol_dipole": 0.5}, em=emo)
ro3 = resp("open_cavity", {"pol_collect": 0.5}, em=emo)
check("FX4: open_cavity pol_dipole = cos^2 enters g^2 (g x sqrt, F_P x pol_dipole), not the collection factor",
      abs(ro2.g_meV / ro1.g_meV - math.sqrt(0.5)) < 1e-12 and abs(ro2.F_P / ro1.F_P - 0.5) < 1e-12 and ro2.pol_factor == 1.0)
check("FX4: open_cavity pol_collect is brightness only (response pol_factor), g and F_P unchanged",
      ro3.pol_factor == 0.5 and ro3.g_meV == ro1.g_meV and ro3.F_P == ro1.F_P)
rg = resp("open_cavity", {"mode": "geometry", "ROC_um": 20.0, "L_um": 2.0}, em=emo)
w0_ref = math.sqrt((0.855 / math.pi) * math.sqrt(2.0 * 18.0))
check("FX4: structure response carries the open-cavity waist w0^2 = (lambda/pi) sqrt(L (ROC - L)) in extra (scene reads it)",
      abs(rg.extra["waist_um"] - w0_ref) < 1e-12 and rg.extra["stable"] is True
      and resp("open_cavity", {"ROC_um": 2.0, "L_um": 5.0}, em=emo).extra["waist_um"] is None)
rm = resp("micropillar", {}, em=em_d)
check("FX4: micropillar response carries the quarter-wave layer thicknesses lambda/4n (scene reads them)",
      abs(rm.extra["t_high_nm"] - em_d.lambda_nm / (4 * 3.5)) < 1e-9 and abs(rm.extra["t_low_nm"] - em_d.lambda_nm / (4 * 2.95)) < 1e-9
      and abs(rm.extra["spacer_unit_nm"] - em_d.lambda_nm / 3.5) < 1e-9)

n_ok = sum(RESULTS)
print(f"{n_ok}/{len(RESULTS)} compose_structures checks passed")
sys.exit(0 if n_ok == len(RESULTS) else 1)
