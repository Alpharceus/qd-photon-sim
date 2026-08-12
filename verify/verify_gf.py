"""Phase-T2 core regression suite (fsim_core.qd_gf): IBM polaron lineshape
with dot-geometry form factor. Same @check style/registry as verify_sde.py.

Run: python verify/verify_gf.py   (exit code 0 iff all pass)

SECOND METHOD (Yuki standard): the module computes the spectrum as
ZPL + FFT of the sideband correlator on a fixed 5 fs / 82 ps grid; check 9
recomputes the SAME Fourier integral by direct quadrature -- an independently
assembled cumulant phi(t) (its own energy grid and cutoff) and an explicit
cos/sin time trapezoid on its own t-grid, transforming the FULL correlator
C(t) (ZPL included, no ZPL/sideband split) -- and demands 1e-4 relative
agreement at spot frequencies. Closed forms (spherical cutoff, T=0
Huang-Rhys) back the remaining checks.

SIDE CONVENTION under test (checks 10/11): on the module's omega axis
(photon energy minus ZPL energy, blue positive) the phonon-EMISSION (Stokes)
side of the PHOTON emission spectrum is omega < 0 -- the photon comes out
red because the emitted LA phonon carries the difference -- so at low T the
red integrated weight must dominate, with the red/blue ratio falling toward 1
(detailed balance exp(E/kBT) -> 1) at 300 K.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import quad

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core import spectral
from fsim_core.qd_gf import (
    HBAR,
    KB,
    PhononParams,
    coupling_alpha_meV,
    coupling_alpha_ps2,
    effective_gamma_zpl,
    huang_rhys,
    ibm_spectrum,
    ibm_transmission,
    sideband_fraction,
    spectral_density,
    zpl_weight,
    _spectrum_parts,
)

CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


P = PhononParams()
C_NMPS = P.c_s_m_s * 1e-3  # nm/ps


# ============================================================ spectral density

@check("J(E): nonnegative everywhere, J(0) = 0, super-exponential decay at "
       "large E (tail < 1e-8 of peak), single broad peak near the geometric "
       "cutoff")
def _():
    E = np.linspace(-2.0, 40.0, 4001)
    J = spectral_density(E, P)
    assert np.all(J >= 0.0)
    assert spectral_density(0.0, P) == 0.0
    assert np.all(J[E <= 0.0] == 0.0)
    assert J.max() > 0.0
    assert J[-1] < 1e-8 * J.max(), J[-1]
    # peak near the E_c band spanned by the two confinement lengths
    e_pk = E[np.argmax(J)]
    ec_lo = HBAR * np.sqrt(2.0) * C_NMPS / max(P.l_xy_nm, P.l_z_nm)
    ec_hi = HBAR * np.sqrt(2.0) * C_NMPS / min(P.l_xy_nm, P.l_z_nm)
    assert 0.5 * ec_lo < e_pk < 3.0 * ec_hi, (e_pk, ec_lo, ec_hi)


@check("PLAUSIBILITY: alpha from the default InP-class material card lands in "
       "the standard InAs/GaAs literature band 0.01-0.1 ps^2 (order-of-"
       "magnitude anchor, not a validation), and the alpha_ps2 override is "
       "honored exactly")
def _():
    a = coupling_alpha_ps2(P)
    assert 0.01 <= a <= 0.1, a  # default card computes ~0.018 ps^2
    assert coupling_alpha_ps2(PhononParams(alpha_ps2=0.027)) == 0.027
    assert abs(coupling_alpha_meV(P) - a / HBAR**2) < 1e-15


@check("spherical limit: anisotropic angular-averaged form factor with "
       "l_xy = l_z = l matches the closed form alpha E^3 exp(-(E/E_c)^2), "
       "E_c = hbar sqrt(2) c_s / l, to 1e-6 relative")
def _():
    for l in (2.0, 3.0, 5.0):
        ps = PhononParams(l_xy_nm=l, l_z_nm=l)
        e_c = HBAR * np.sqrt(2.0) * C_NMPS / l
        E = np.linspace(0.01, 4.0 * e_c, 400)
        j_num = spectral_density(E, ps)
        j_closed = coupling_alpha_meV(ps) * E**3 * np.exp(-((E / e_c) ** 2))
        assert np.max(np.abs(j_num - j_closed) / j_closed) < 1e-6, l


# ============================================================ Huang-Rhys / ZPL

@check("geometry monotonicity: smaller dot -> larger S_total (l = 2 nm vs "
       "5 nm, spherical, T = 4 K), tracking the closed-form 1/l^2 trend")
def _():
    s2 = huang_rhys(PhononParams(l_xy_nm=2.0, l_z_nm=2.0), 4.0)
    s5 = huang_rhys(PhononParams(l_xy_nm=5.0, l_z_nm=5.0), 4.0)
    assert s2 > s5, (s2, s5)
    # T=0 closed form says exactly (5/2)^2 = 6.25 (checked tight); at 4 K the
    # softer dot (l=5, lower E_c) picks up proportionally more coth
    # enhancement, so the ratio drops below 6.25 (numerically ~4.1) -- demand
    # only that the strong 1/l^2 trend survives (> 3x)
    s2_0 = huang_rhys(PhononParams(l_xy_nm=2.0, l_z_nm=2.0), 0.0)
    s5_0 = huang_rhys(PhononParams(l_xy_nm=5.0, l_z_nm=5.0), 0.0)
    assert abs(s2_0 / s5_0 - 6.25) < 0.01, s2_0 / s5_0
    assert 3.0 < s2 / s5 < 6.25, s2 / s5


@check("S_total(T) monotone increasing on [0, 300] K; Z(T) = exp(-S_total) in "
       "(0, 1]; T = 0 matches an independent scipy.quad of J/E^2 to 5e-5 and "
       "the spherical closed form alpha E_c^2/2 to 1e-4")
def _():
    temps = [0.0, 4.0, 30.0, 77.0, 150.0, 300.0]
    svals = [huang_rhys(P, t) for t in temps]
    for a, b in zip(svals, svals[1:]):
        assert b > a, (a, b)
    for t, s in zip(temps, svals):
        z = zpl_weight(P, t)
        assert 0.0 < z <= 1.0 and abs(z - np.exp(-s)) < 1e-12, (t, z)
    e_hi = 15.0 * HBAR * np.sqrt(2.0) * C_NMPS / min(P.l_xy_nm, P.l_z_nm)
    s0_quad, _err = quad(lambda e: spectral_density(e, P) / e**2, 0.0, e_hi, limit=200)
    assert abs(huang_rhys(P, 0.0) - s0_quad) / s0_quad < 5e-5, (huang_rhys(P, 0.0), s0_quad)
    ps = PhononParams(l_xy_nm=3.0, l_z_nm=3.0)
    e_c = HBAR * np.sqrt(2.0) * C_NMPS / 3.0
    s0_closed = coupling_alpha_meV(ps) * e_c**2 / 2.0
    assert abs(huang_rhys(ps, 0.0) - s0_closed) / s0_closed < 1e-4


# ========================================================== coupling -> 0 limit

P_TINY = PhononParams(alpha_ps2=coupling_alpha_ps2(P) * 1e-6)


@check("coupling -> 0 (alpha x 1e-6): spectrum collapses to the pure "
       "Lorentzian of width gamma_zpl (max abs deviation < 1e-4 of the "
       "Lorentzian peak, same-grid normalization)")
def _():
    grid = np.linspace(-25.0, 25.0, 50001)
    s = ibm_spectrum(grid, P_TINY, 77.0, 0.1)
    lor = spectral.lorentzian(grid, 0.0, 0.1)
    lor = lor / np.trapezoid(lor, grid)
    assert np.max(np.abs(s - lor)) / lor.max() < 1e-4


@check("coupling -> 0: ibm_transmission reproduces spectral.tophat_"
       "transmission AND spectral.cavity_transmission to 1e-4 at several "
       "line-filter offsets (acceptance definitions mirrored exactly)")
def _():
    for d in (0.0, 0.4, 1.2):
        t_top = ibm_transmission(d, P_TINY, 77.0, 0.1, w_meV=1.5)
        assert abs(t_top - spectral.tophat_transmission(d, 0.1, 1.5)) < 1e-4, d
        t_cav = ibm_transmission(d, P_TINY, 77.0, 0.1, kappa_meV=0.5)
        assert abs(t_cav - spectral.cavity_transmission(d, 0.1, 0.5)) < 1e-4, d
    # no-filter branch mirrors spectral.transmission's 1.0
    assert ibm_transmission(0.3, P_TINY, 77.0, 0.1) == 1.0


# ================================================================ spectrum norm

@check("spectrum: unit area on the grid to 1e-6, nonnegative everywhere "
       "(pre-clip FFT ringing < 1e-10 of peak, checked on the raw parts), "
       "and raw ZPL + sideband weight accounts for the whole line to 5e-3")
def _():
    grid = np.linspace(-30.0, 30.0, 60001)
    s = ibm_spectrum(grid, P, 77.0, 0.1)
    assert np.all(s >= 0.0)
    assert abs(np.trapezoid(s, grid) - 1.0) < 1e-6
    omega, s_sb, z = _spectrum_parts(P, 77.0, 0.1)
    peak = z * spectral.lorentzian(0.0, 0.0, 0.1) + s_sb.max()
    assert s_sb.min() > -1e-10 * peak, s_sb.min()  # ringing budget (clipped in module)
    assert abs(z + np.trapezoid(s_sb, omega) - 1.0) < 5e-3


# =============================================================== second method

@check("SECOND METHOD: direct quadrature of the same Fourier integral -- "
       "independently assembled phi(t) (own E-grid) and explicit cos/sin "
       "time trapezoid of the FULL correlator (no ZPL split) -- matches the "
       "module's FFT spectrum to 1e-4 relative at 5 spot frequencies "
       "(T = 77 K, gamma_zpl = 0.8 meV)")
def _():
    T, gam = 77.0, 0.8
    e_hi = 15.0 * HBAR * np.sqrt(2.0) * C_NMPS / min(P.l_xy_nm, P.l_z_nm)
    E = np.linspace(0.0, e_hi, 8001)
    J = spectral_density(E, P)
    g = np.zeros_like(E)
    g[1:] = J[1:] / E[1:] ** 2
    p = np.empty_like(E)
    p[1:] = g[1:] / np.tanh(E[1:] / (2.0 * KB * T))
    p[0] = 2.0 * coupling_alpha_meV(P) * KB * T
    t = np.arange(0.0, 30.0, 0.004)
    phi = np.empty(t.size, dtype=complex)
    for i in range(0, t.size, 200):
        ph = np.multiply.outer(t[i:i + 200], E) / HBAR
        phi[i:i + 200] = (np.trapezoid(p * np.cos(ph), E, axis=1)
                          - 1j * np.trapezoid(g * np.sin(ph), E, axis=1))
    f = np.exp(phi - phi[0]) * np.exp(-gam * t / (2.0 * HBAR))
    spots = np.array([-2.0, -0.75, -0.3, 0.5, 1.5])
    s2 = np.array([
        np.trapezoid(f.real * np.cos(w * t / HBAR) + f.imag * np.sin(w * t / HBAR), t)
        for w in spots
    ]) / (np.pi * HBAR)
    omega, s_sb, z = _spectrum_parts(P, T, gam)
    s1 = z * spectral.lorentzian(spots, 0.0, gam) + np.interp(spots, omega, s_sb)
    rel = np.abs(s1 - s2) / np.abs(s1)
    assert np.max(rel) < 1e-4, rel


# ========================================================== detailed balance

@check("side convention / detailed balance at 4 K: integrated weight on the "
       "RED side (omega < 0, phonon-EMISSION side for photon emission) "
       "exceeds the blue (anti-Stokes) side by > 3x -- the T -> 0 one-sided "
       "sideband, sign convention as documented")
def _():
    grid = np.linspace(-15.0, 15.0, 30001)
    s = ibm_spectrum(grid, P, 4.0, 0.001)
    red = np.trapezoid(s[grid < -0.05], grid[grid < -0.05])
    blue = np.trapezoid(s[grid > 0.05], grid[grid > 0.05])
    assert red > blue
    assert red / blue > 3.0, red / blue


@check("detailed balance at 300 K: red/blue integrated ratio falls to within "
       "tens of percent of 1 (kB T >> E_c washes out the asymmetry) while "
       "staying > 1 (exact KMS ordering)")
def _():
    grid = np.linspace(-15.0, 15.0, 30001)
    s = ibm_spectrum(grid, P, 300.0, 0.001)
    red = np.trapezoid(s[grid < -0.05], grid[grid < -0.05])
    blue = np.trapezoid(s[grid > 0.05], grid[grid > 0.05])
    ratio = red / blue
    assert 1.0 < ratio < 1.4, ratio


# ============================================================== transmissions

@check("wide-open filter: ibm_transmission with w = 1e4 meV collects the "
       "whole line (-> 1 within 1e-4) even at 300 K where the sideband "
       "carries ~90% of the weight")
def _():
    t_open = ibm_transmission(0.7, P, 300.0, 0.1, w_meV=1e4)
    assert abs(t_open - 1.0) < 1e-4, t_open


@check("sideband fraction: 3 nm spherical dot has 1 - Z in (0, 1) at both "
       "4 K and 300 K, larger at 300 K, and equal to 1 - zpl_weight exactly")
def _():
    ps = PhononParams(l_xy_nm=3.0, l_z_nm=3.0)
    sf4 = sideband_fraction(ps, 4.0)
    sf300 = sideband_fraction(ps, 300.0)
    assert 0.0 < sf4 < 1.0 and 0.0 < sf300 < 1.0, (sf4, sf300)
    assert sf300 > sf4, (sf4, sf300)
    assert abs(sf4 - (1.0 - zpl_weight(ps, 4.0))) < 1e-15


# ======================================================== bridge / three-layer

@check("broadening bridge + three-layer rule: effective_gamma_zpl is "
       "bit-identical to the validated spectral.gamma_of_T form, and "
       "fsim_core.qd_gf has no GUI/plotting imports")
def _():
    temps = np.array([0.1, 4.0, 77.0, 200.0, 300.0])
    a = effective_gamma_zpl(temps, 0.002, 1e-4, 5.0, 30.0)
    b = spectral.gamma_of_T(temps, 0.002, 1e-4, 5.0, 30.0)
    assert np.array_equal(a, b)
    src = (Path(__file__).resolve().parents[1] / "fsim_core" / "qd_gf.py").read_text(
        encoding="utf-8")
    for banned in ("tkinter", "PyQt", "PySide", "wx", "matplotlib"):
        assert banned not in src, banned


# ==================================================== device wiring (T2 gate)

@check("device wiring: dot.lineshape default 'lorentzian' keeps evaluate() "
       "BIT-IDENTICAL to legacy (full default-grid curve comparison)")
def _():
    from fsim_core.device import DeviceDesign, evaluate
    ref = evaluate(DeviceDesign())
    d = DeviceDesign()
    assert d.dot.lineshape == "lorentzian"
    same = evaluate(d)
    for key in ("g2", "eps", "rho2"):
        assert np.array_equal(ref["curves"][key], same["curves"][key],
                              equal_nan=True), key


@check("device wiring + T2 GATE: IBM lineshape fattens eps in the "
       "Lorentzian-optimism direction at the staged decision points "
       "(77/120 K, delta=5) and converges toward Lorentzian at 300 K where "
       "Gamma already dominates; both finite everywhere")
def _():
    from fsim_core.device import DeviceDesign, evaluate
    Tg = [77.0, 120.0, 300.0]
    lor = DeviceDesign(); lor.dot.delta_xx = 5.0
    ibm = DeviceDesign(); ibm.dot.delta_xx = 5.0; ibm.dot.lineshape = "ibm"
    el = evaluate(lor, T_grid=Tg)["curves"]["eps"]
    ei = evaluate(ibm, T_grid=Tg)["curves"]["eps"]
    assert np.all(np.isfinite(el)) and np.all(np.isfinite(ei))
    # sidebands leak the XX line into the X window: eps_ibm > eps_lor cold
    assert ei[0] > 1.2 * el[0], (ei[0], el[0])   # 77 K: >20% fattening
    assert ei[1] > 1.2 * el[1], (ei[1], el[1])   # 120 K
    # 300 K: the fitted Gamma dominates; correction shrinks to the % level
    assert abs(ei[2] / el[2] - 1.0) < 0.05, (ei[2], el[2])


@check("device wiring: YAML round-trip preserves lineshape + phonon dict; "
       "phonon geometry override actually changes eps (smaller dot -> "
       "stronger sidebands -> larger eps at 77 K)")
def _():
    import tempfile
    from pathlib import Path as _P
    from fsim_core.device import DeviceDesign, evaluate
    d = DeviceDesign()
    d.dot.lineshape = "ibm"
    d.dot.phonon = {"l_xy_nm": 3.0, "l_z_nm": 1.0}
    with tempfile.TemporaryDirectory() as td:
        p = _P(td) / "card.yaml"
        d.save(p)
        d2 = DeviceDesign.load(p)
    assert d2.dot.lineshape == "ibm" and d2.dot.phonon == d.dot.phonon
    big = DeviceDesign(); big.dot.lineshape = "ibm"
    big.dot.phonon = {"l_xy_nm": 6.0, "l_z_nm": 3.0}
    small = DeviceDesign(); small.dot.lineshape = "ibm"
    small.dot.phonon = {"l_xy_nm": 3.0, "l_z_nm": 1.0}
    eb = evaluate(big, T_grid=[77.0])["curves"]["eps"][0]
    es = evaluate(small, T_grid=[77.0])["curves"]["eps"][0]
    assert es > eb, (es, eb)


@check("V-a GATE (T2): the six Chatzarakis points re-evaluated with the IBM "
       "lineshape at the FITTED Lorentzian parameters shift eps by < 0.03 "
       "per point (inside the V-a acceptance band) -- the validated fit "
       "survives the sideband correction without refitting. Arsenide-class "
       "coupling (Ramsay alpha=0.027 ps^2 [DR]); an eps-level statement at "
       "fixed params, not a full refit (fitting-layer IBM option is the "
       "noted follow-up)")
def _():
    import json
    from fsim_core.spectral import epsilon as _eps
    fp = json.loads((Path(__file__).resolve().parents[1] / "out" / "phase0" /
                     "fit_params.json").read_text())["params"]
    dxx = 5.9   # cards/chatzarakis.yaml delta_xx [DR]
    pp = PhononParams(alpha_ps2=0.027, c_s_m_s=4780.0, l_xy_nm=4.5, l_z_nm=2.0)
    pts = [(78.0, 1.40, -0.37), (120.0, 3.19, -0.70), (150.0, 4.12, -1.27),
           (170.0, 4.74, -1.65), (210.0, 6.69, -2.82), (230.0, 7.66, -3.41)]
    for T, w, dx in pts:
        gam = float(spectral.gamma_of_T(T, fp["gamma0"], fp["a_ac"],
                                        fp["b_lo"], fp["E_lo"]))
        gxx = fp["r_xx"] * gam
        el = _eps(dxx, gam, gxx, w=w, dx=dx).eps
        ei = (float(ibm_transmission(dx - dxx, pp, T, gxx, w_meV=w))
              / float(ibm_transmission(dx, pp, T, gam, w_meV=w)))
        assert ei > el, (T, "sideband correction must fatten eps")
        assert abs(ei - el) < 0.03, (T, el, ei)


def main():
    failed = 0
    for name, fn in CHECKS:
        try:
            fn()
            print(f"  PASS  {name}")
        except AssertionError as e:
            failed += 1
            print(f"* FAIL  {name}  {e}")
    n = len(CHECKS)
    print(f"\n{n - failed}/{n} checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
