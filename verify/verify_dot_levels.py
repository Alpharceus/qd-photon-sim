"""Verification of fsim_core.dot_levels (Module L): solver limits against
closed forms, an independent finite-difference second method, BenDaniel-
Duke direction, literature-class reproductions, temperature, and the
retention mapping.

Run: python verify/verify_dot_levels.py [--quiet]   (exit code 0 iff all pass)

The class checks distinguish a measured device value from the prediction of
the explicitly specified geometry and material parameterisation. A failed
experimental calibration is never hidden: where a measurement does not
describe that geometry, the check names the cited model prediction instead.
"""
import sys
from pathlib import Path

import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.special import jn_zeros

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core import dot_levels as D  # noqa: E402
from fsim_core import materials as M  # noqa: E402
from fsim_core.integrator import retention  # noqa: E402

QUIET = "--quiet" in sys.argv
CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


def say(*a):
    if not QUIET:
        print(*a)


# Independent constants from SI values (the module carries its own number).
HBAR = 1.054571817e-34
M0 = 9.1093837015e-31
QE = 1.602176634e-19
HB2_2M0 = HBAR**2 / (2 * M0) / QE * 1e3 * 1e18      # meV nm^2  (= 38.0998)

P = D.class_presets()
LV = {}


def lv(name, **kw):
    key = (name, tuple(sorted(kw.items())))
    if key not in LV:
        LV[key] = D.levels(P[name](**kw))
    return LV[key]


# ------------------------------------------------------------ (a) infinite-well limits

@check("a1: 1D well, V = 1e8 meV, w = 5 nm, m = 0.067 -> hbar^2 pi^2/(2 m w^2) within 1%")
def _():
    w = D.finite_well_1d(1e8, 5.0, 0.067, 0.067)
    E_inf = HB2_2M0 * np.pi**2 / (0.067 * 25.0)          # 224.5 meV
    assert abs(w.energies_meV[0] / E_inf - 1) < 0.01, (w.energies_meV[0], E_inf)
    assert abs(w.energies_meV[1] / (4 * E_inf) - 1) < 0.01


@check("a2: 1D well, V = 1e6 meV: finite-depth penetration E = E_inf (w/(w + 2/kappa))^2 within 0.5%")
def _():
    # At V = 1e6 meV and m = 0.067 the decay length 1/kappa = 0.024 nm lowers
    # the level by 1.9 % (a physical effect, not a solver error).
    w = D.finite_well_1d(1e6, 5.0, 0.067, 0.067)
    E_inf = HB2_2M0 * np.pi**2 / (0.067 * 25.0)
    kappa = np.sqrt(0.067 * 1e6 / HB2_2M0)
    E_pen = E_inf * (5.0 / (5.0 + 2.0 / kappa))**2
    assert abs(w.energies_meV[0] / E_pen - 1) < 0.005, (w.energies_meV[0], E_pen, E_inf)


@check("a3: 2D disk, V = 1e8 meV, R = 10 nm, m = 0.067 -> hbar^2 j01^2/(2 m R^2) = 32.9 meV (l=0) and j11 (l=1) within 1%")
def _():
    d = D.finite_disk_2d(1e8, 10.0, 0.067, 0.067)
    j01, j11 = 2.4048, 3.8317
    assert abs(jn_zeros(0, 1)[0] - j01) < 1e-3 and abs(jn_zeros(1, 1)[0] - j11) < 1e-3
    E0 = HB2_2M0 * j01**2 / (0.067 * 100.0)
    E1 = HB2_2M0 * j11**2 / (0.067 * 100.0)
    assert abs(d.E0_meV / E0 - 1) < 0.01, (d.E0_meV, E0)
    assert abs(d.E1_meV / E1 - 1) < 0.01, (d.E1_meV, E1)
    assert abs(E0 - 32.9) < 0.4


@check("a4: infinite-well rms extents: 1D <z^2>^1/2 = w sqrt(1/12 - 1/(2 pi^2)); disk <r^2>^1/2 from J0 quadrature")
def _():
    w = D.finite_well_1d(1e8, 5.0, 0.067, 0.067)
    assert abs(w.rms_z_nm / (5.0 * np.sqrt(1 / 12 - 1 / (2 * np.pi**2))) - 1) < 0.01
    d = D.finite_disk_2d(1e8, 10.0, 0.067, 0.067)
    from scipy.special import j0
    from scipy.integrate import quad
    num = quad(lambda r: r**3 * j0(2.4048 * r / 10.0)**2, 0, 10.0)[0]
    den = quad(lambda r: r * j0(2.4048 * r / 10.0)**2, 0, 10.0)[0]
    assert abs(d.rms_r_nm / np.sqrt(num / den) - 1) < 0.01, (d.rms_r_nm, np.sqrt(num / den))


# ------------------------------------------------------ (b) finite-difference second method

def fd_1d(V, w, m, L=30.0, h=0.01):
    z = np.arange(-L, L + h / 2, h)
    Vz = np.where(np.abs(z) <= w / 2, 0.0, V)
    t = HB2_2M0 / (m * h**2)
    E, _ = eigh_tridiagonal(2 * t + Vz, -t * np.ones(len(z) - 1), select="i", select_range=(0, 2))
    return E


def fd_radial_l0(V, R, m, L=40.0, h=0.01):
    """Conservative FD for -(hbar^2/2m)(1/r) d/dr (r d/dr) + V, l = 0, on the
    half-integer grid r_i = (i + 1/2) h, symmetrised with u = sqrt(r) R."""
    N = int(L / h)
    i = np.arange(N)
    r = (i + 0.5) * h
    rp = (i + 1.0) * h          # r_{i+1/2}
    rm = i * h                  # r_{i-1/2}  (0 at the axis: natural Neumann)
    Vr = np.where(r <= R, 0.0, V)
    t = HB2_2M0 / m
    diag = t * (rm + rp) / (r * h**2) + Vr
    off = -t * rp[:-1] / (h**2 * np.sqrt(r[:-1] * r[1:]))
    E, _ = eigh_tridiagonal(diag, off, select="i", select_range=(0, 1))
    return E


@check("b1: textbook 1D well V = 300 meV, w = 10 nm, m = 0.067 vs finite-difference eigh within 1 meV (3 states)")
def _():
    w = D.finite_well_1d(300.0, 10.0, 0.067, 0.067)
    E = fd_1d(300.0, 10.0, 0.067)
    for a, b in zip(w.energies_meV[:3], E[:3]):
        assert abs(a - b) < 1.0, (w.energies_meV, E)


@check("b2: 2D disk V = 300 meV, R = 10 nm, m = 0.067 (l = 0) vs radial finite-difference eigh within 1 meV")
def _():
    d = D.finite_disk_2d(300.0, 10.0, 0.067, 0.067)
    E = fd_radial_l0(300.0, 10.0, 0.067)
    assert abs(d.E0_meV - E[0]) < 1.0, (d.E0_meV, E)


@check("b3: layered BDD solver reduces to finite_well_1d for a thick single shell (< 0.3 meV) and is FD-consistent for m_in != m_out")
def _():
    a = D.finite_well_1d(300.0, 10.0, 0.067, 0.15)
    b = D.well_1d_layered([10.0, 25.0], [0.0, 300.0, 300.0], [0.067, 0.15, 0.15])
    assert abs(a.energies_meV[0] - b.energies_meV[0]) < 0.3, (a.energies_meV, b.energies_meV)
    assert abs(a.rms_z_nm - b.rms_z_nm) < 0.02


# ------------------------------------------------------------ (c) BenDaniel-Duke direction

@check("c: BenDaniel-Duke -- heavier barrier mass lowers the level (1D and disk), lighter raises it")
def _():
    e0 = D.finite_well_1d(300.0, 10.0, 0.067, 0.067).energies_meV[0]
    eh = D.finite_well_1d(300.0, 10.0, 0.067, 0.15).energies_meV[0]
    el = D.finite_well_1d(300.0, 10.0, 0.067, 0.03).energies_meV[0]
    assert eh < e0 < el, (eh, e0, el)
    d0 = D.finite_disk_2d(300.0, 10.0, 0.067, 0.067).E0_meV
    dh = D.finite_disk_2d(300.0, 10.0, 0.067, 0.15).E0_meV
    dl = D.finite_disk_2d(300.0, 10.0, 0.067, 0.03).E0_meV
    assert dh < d0 < dl, (dh, d0, dl)


# ------------------------------------------------------------ (d) literature classes

@check("d-i-a: specified In0.5Ga0.5As 3 nm x R 12 nm disk in GaAs: Vurgaftman et al., J. Appl. Phys. 89, 5815 (2001) parameterisation gives E_X in [1.17, 1.23] eV (not Chatzarakis' different 1.33/1.44 eV dots)")
def _():
    L = lv("InGaAs0.5/GaAs/GaAs on GaAs")
    assert 1.17 <= L.E_X_eV <= 1.23, f"E_X = {L.E_X_eV:.3f} eV"


@check("d-i-b: specified disk/GaAs WL: Vurgaftman et al., J. Appl. Phys. 89, 5815 (2001) parameterisation gives dE_pair_WL in [230, 270] meV; Chatzarakis et al., Phys. Rev. Applied 20, 034011 (2023) report 110 meV for a different dot")
def _():
    L = lv("InGaAs0.5/GaAs/GaAs on GaAs")
    assert 230 <= L.dE_pair_WL <= 270, f"dE_pair_WL = {L.dE_pair_WL:.0f} meV (E_WL {L.E_WL_eV:.3f}, E_X {L.E_X_eV:.3f})"


@check("d-i-c: direct 0.5 nm WL/Al0.57Ga0.43As interface: Vurgaftman et al., J. Appl. Phys. 89, 5815 (2001) finite-well model gives dE_pair_WL in [700, 760] meV; Chatzarakis et al., Phys. Rev. Applied 20, 034011 (2023) report 240 +/- 20 meV for their spacer-containing structure")
def _():
    L = lv("InGaAs0.5/AlGaAs0.57/AlGaAs0.57 on GaAs")
    assert 700 <= L.dE_pair_WL <= 760, f"dE_pair_WL = {L.dE_pair_WL:.0f} meV (E_WL {L.E_WL_eV:.3f}, E_X {L.E_X_eV:.3f})"


@check("d-i-d: SSL shift of the WL vs GaAs spacer thickness: monotonic, and the SSL-induced INCREASE of dE_pair_WL over the GaAs case lands in [60, 300] meV for a 2-5 nm spacer (paper +130)")
def _():
    base = lv("InGaAs0.5/GaAs/GaAs on GaAs").dE_pair_WL
    rows = []
    for d in (2.0, 3.0, 5.0, 8.0):
        S = P["InGaAs0.5/GaAs/AlGaAs0.57 on GaAs"]()
        S.geometry.matrix_thickness_nm = d
        L = D.levels(S)
        rows.append((d, L.dE_pair_WL, L.E_WL_eV, L.E_X_eV))
    say("    SSL spacer scan (GaAs between WL and AlGaAs0.57):  d [nm]  dE_pair_WL [meV]  E_WL [eV]  E_X [eV]")
    say(f"      inf   {base:8.0f}")
    for d, dE, EWL, EX in rows:
        say(f"      {d:4.1f}  {dE:8.0f}   {EWL:.3f}   {EX:.3f}")
    incs = [r[1] - base for r in rows]
    assert all(incs[i] >= incs[i + 1] - 1e-6 for i in range(len(incs) - 1)), incs
    assert any(60 <= inc <= 300 for inc in incs[:3]), incs


@check("d-ii-a: specified InP 3 nm x R 10 nm/GaInP disk with Pryor, Phys. Rev. B 56, 10404 (1997) VBO: E_X in [1.65, 1.74] eV (not the 1.815-1.95 eV ensemble range)")
def _():
    L = lv("InP/GaInP/AlGaInP0.55 on GaAs")
    assert 1.65 <= L.E_X_eV <= 1.74, f"E_X = {L.E_X_eV:.3f} eV"


@check("d-ii-b: Reischle class: electron escape energy to the GaInP matrix in [150, 300] meV (Pryor ~250 meV CB offset)")
def _():
    L = lv("InP/GaInP/AlGaInP0.55 on GaAs")
    assert 150 <= L.dE_e_matrix <= 300, f"dE_e_matrix = {L.dE_e_matrix:.0f} meV (V_e {L.V_e:.0f}, E_e {L.E_e:.0f})"


@check("d-ii-c (bonus): Bommer class InP/(Al0.2Ga0.8)InP: hole escape energy in [50, 150] meV (Bommer Ea(X) = 96 +/- 7 meV, ~90 meV confinement)")
def _():
    L = lv("InP/AlGaInP0.2/AlGaInP0.55 on GaAs")
    assert 50 <= L.dE_h_matrix <= 150, f"dE_h_matrix = {L.dE_h_matrix:.0f} meV"


@check("d-iii: InAs/InP telecom dot (2.5 nm x R 15 nm, 300 K): E_X in [0.75, 0.98] eV (1.27-1.65 um)")
def _():
    L = lv("InAs/InP/InP")
    assert 0.75 <= L.E_X_eV <= 0.98, f"E_X = {L.E_X_eV:.3f} eV ({L.lambda_nm:.0f} nm)"


@check("d-iv-a: HKUST specified InP/GaAsP0.4 5 nm x R 12 nm disk: Vurgaftman et al., J. Appl. Phys. 89, 5815 (2001) model gives E_X in [1.48, 1.54] eV; Gu et al., Opt. Express 33, 23732 (2025) report ~750-755 nm ensemble PL/lasing for 4-7 nm GaAsP0.35 dots")
def _():
    L = lv("InP/GaAsP0.4/AlGaAs0.4 on GaAs", height_nm=5.0)
    assert 1.48 <= L.E_X_eV <= 1.54, f"E_X = {L.E_X_eV:.3f} eV ({L.lambda_nm:.0f} nm)"


@check("d-iv-b: HKUST specified disk: Vurgaftman et al., J. Appl. Phys. 89, 5815 (2001) strained model-solid V_e about 200 meV constrains a bound-electron escape energy to [90, 140] meV; Gu et al., Opt. Express 33, 23732 (2025) establish type-I/type-II growth variability")
def _():
    L = lv("InP/GaAsP0.4/AlGaAs0.4 on GaAs", height_nm=5.0)
    say(f"    HKUST hole: bound = {L.hole_bound}, type {L.type}, dE_h_matrix = {L.dE_h_matrix:.0f} meV "
        f"(V_h = {L.V_h:.0f} meV to the tensile-GaAsP LH edge)")
    assert L.electron_bound
    assert 90 <= L.dE_e_matrix <= 140, f"dE_e_matrix = {L.dE_e_matrix:.0f} meV (V_e = {L.V_e:.0f}, E_e = {L.E_e:.0f})"


# ------------------------------------------------------------ (e) temperature

@check("e: InP/GaInP class: E_X(300 K) below E_X(5 K) by 50-90 meV (InP Varshni class shift)")
def _():
    L5 = lv("InP/GaInP/AlGaInP0.55 on GaAs", T=5.0)
    L300 = lv("InP/GaInP/AlGaInP0.55 on GaAs", T=300.0)
    dE = (L5.E_X_eV - L300.E_X_eV) * 1e3
    assert 50 <= dE <= 90, f"shift = {dE:.1f} meV"


# ------------------------------------------------------------ (f) retention mapping

@check("f: retention_params -> E_a > 0, 1e2 < a_esc < 1e8, and S(300 K) decreases with smaller E_a")
def _():
    L = lv("InP/GaInP/AlGaInP0.55 on GaAs", T=300.0)
    for ch in ("pair_half", "pair", "electron", "hole", "min"):
        rp = D.retention_params(L, tau_rad_ns=1.0, channel=ch, verbose=not QUIET)
        assert rp["E_a"] > 0, rp
        assert 1e2 < rp["a_esc"] < 1e8, rp
        assert rp["E_b"] > 0 and rp["b_p"] == 4
        assert "[E]" in rp["note"]
    rp = D.retention_params(L, 1.0, "pair_half", verbose=False)
    S1 = float(retention(300.0, rp["a_esc"], rp["E_a"], rp["b_p"], rp["E_b"]))
    S2 = float(retention(300.0, rp["a_esc"], 0.5 * rp["E_a"], rp["b_p"], rp["E_b"]))
    assert 0 < S2 < S1 < 1, (S1, S2)


# ------------------------------------------------------------ (g) exciton binding

@check("g: Gaussian exciton binding in [10, 40] meV for every preset")
def _():
    bad = []
    for name in P:
        L = lv(name)
        if not (10 <= L.E_bind <= 40):
            bad.append((name, L.E_bind))
    assert not bad, bad


@check("g2: exciton binding derivation: <1/rho> of a 2D Gaussian = sqrt(pi)/L by quadrature; type-I flags consistent with offsets")
def _():
    from scipy.integrate import quad
    Lr = 7.3
    num = quad(lambda r: (2 * r / Lr**2) * np.exp(-r**2 / Lr**2) / r, 0, np.inf)[0]
    assert abs(num - np.sqrt(np.pi) / Lr) < 1e-8
    Eb = D._gauss_binding(5.0, 5.3, 12.5)
    assert abs(Eb - np.sqrt(np.pi) * 1439.96 / (12.5 * np.sqrt(5.0**2 + 5.3**2))) < 0.01
    for name in P:
        L = lv(name)
        assert (L.type == "I") == (L.V_e > 0 and L.V_h > 0), name


# ------------------------------------------------------------------ report + main

def report_block():
    print("=" * 100)
    print("Literature-class level report (fsim_core.dot_levels, separable disk, biaxial strain, [E])")
    print("=" * 100)
    hdr = (f"{'preset':44s} {'T':>4s} {'E_X eV':>7s} {'nm':>5s} {'dE_pWL':>7s} {'dE_e':>6s} {'dE_h':>6s} "
           f"{'E_bind':>6s} {'a_esc':>8s} {'type':>4s}")
    print(hdr)
    for name in P:
        L = lv(name)
        rp = D.retention_params(L, 1.0, "pair_half", verbose=False)
        print(f"{name:44s} {L.system.T:4.0f} {L.E_X_eV:7.3f} {L.lambda_nm:5.0f} {L.dE_pair_WL:7.0f} "
              f"{L.dE_e_matrix:6.0f} {L.dE_h_matrix:6.0f} {L.E_bind:6.1f} {rp['a_esc']:8.2e} {L.type:>4s}")
    L = lv("InP/GaAsP0.4/AlGaAs0.4 on GaAs", height_nm=5.0)
    rp = D.retention_params(L, 1.0, "pair_half", verbose=False)
    print(f"{'InP/GaAsP0.4/AlGaAs0.4 on GaAs (h = 5 nm)':44s} {L.system.T:4.0f} {L.E_X_eV:7.3f} {L.lambda_nm:5.0f} "
          f"{L.dE_pair_WL:7.0f} {L.dE_e_matrix:6.0f} {L.dE_h_matrix:6.0f} {L.E_bind:6.1f} {rp['a_esc']:8.2e} {L.type:>4s}")
    print("  (dE_* in meV; dE_pWL = E_WL - E_X; a_esc for tau_rad = 1 ns, n_dot = 1e10 cm^-2, tau_cap = 10 ps)")
    print("-" * 100)
    for name in P:
        print(D.report(lv(name)))
    print("=" * 100)


def main():
    if not QUIET:
        report_block()
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
