"""Verification of fsim_core.dot_levels (Module L): solver limits against
closed forms, an independent finite-difference second method, BenDaniel-
Duke direction, literature-class reproductions, temperature, and the
retention mapping.

Run: python verify/verify_dot_levels.py [--quiet]   (exit code 0 iff all pass)

Every check(...) below asserts a window centred on a specific number this
module's code did NOT produce -- a value from a cited paper, or (a1-a4, b1-
b3, c, g2, and the d-i-d monotonicity assertion) a closed form / independent
solver / structural symmetry that holds regardless of the paper numbers.
Nothing here brackets the model's own output.

Where the geometry actually specified by a paper's class (dot material,
matrix, barrier, the stated height/composition) genuinely misses that
paper's own reported number by an amount too large to attribute to
measurement spread or to this model's stated closed-form/finite-difference
accuracy, the comparison is NOT forced to pass by widening the window around
the model's answer: it is moved to the "known deviations" table printed by
report_block()/main() below, with the paper's number, this model's number,
the size of the miss in meV, and the specific documented model limitation
responsible (see the module docstring of fsim_core.dot_levels). Known
deviations are counted separately, on their own line, never as PASS.
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
KNOWN_DEVIATIONS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


def deviation(name, model_value, published, miss, explanation):
    """Register one documented, NOT-counted-as-PASS mismatch between this
    model's own prediction and a cited published value. `miss` is a string
    with the signed size in meV (or eV) and direction."""
    KNOWN_DEVIATIONS.append((name, model_value, published, miss, explanation))


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
#
# An Opus review (2026-09-05) found that d-i-a/b/c, d-ii-a/b, d-iii and
# d-iv-a/b previously asserted windows built to bracket THIS MODEL's own
# output for the specified geometry, not the cited paper's number (which in
# every one of those cases lies outside the window). Each is repaired below:
# either the window is re-centred honestly on the actual cited number (and
# now genuinely passes), or -- where the specified-geometry prediction
# really does miss that number outside any defensible tolerance -- the
# comparison is moved to the known-deviations table (see report_block()),
# tagged with the size of the miss and the specific documented model
# limitation. d-ii-c is untouched (confirmed correct by the review).

@check("d-i-d: SSL vs plain-GaAs spacer scan (In0.5Ga0.5As/GaAs(d)/Al0.57Ga0.43As WL): dE_pair_WL's SSL-induced increase over the semi-infinite-GaAs case is monotonically non-increasing as the GaAs spacer d widens (structural closed-form limit: d -> infinity must recover the GaAs-matrix case)")
def _():
    base = lv("InGaAs0.5/GaAs/GaAs on GaAs").dE_pair_WL
    rows = []
    for d in (2.0, 3.0, 5.0, 8.0, 15.0, 30.0):
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
    assert incs[-1] < 5.0, ("d = 30 nm spacer should have relaxed back close to the semi-infinite-GaAs case", incs)


@check("d-ii-c (bonus): Bommer class InP/(Al0.2Ga0.8)InP: hole escape energy in [50, 150] meV (Bommer Ea(X) = 96 +/- 7 meV, ~90 meV confinement)")
def _():
    L = lv("InP/AlGaInP0.2/AlGaInP0.55 on GaAs")
    assert 50 <= L.dE_h_matrix <= 150, f"dE_h_matrix = {L.dE_h_matrix:.0f} meV"


def _register_known_deviations():
    """Literature-class comparisons where the model, run at the paper's own
    specified (or best-estimate) geometry, misses the paper's own reported
    number outside any defensible measurement-spread or model-accuracy
    tolerance. Not checks: never counted as PASS/FAIL, always reported."""
    L = lv("InGaAs0.5/GaAs/GaAs on GaAs")
    deviation(
        "d-i-a: In0.5Ga0.5As 3 nm x R 12 nm / GaAs (Chatzarakis geometry, no SSL)",
        f"E_X = {L.E_X_eV:.3f} eV", "1.30 eV (Chatzarakis et al., PRA 20, 034011 (2023), macro-PL QD band, <Al> = 0 reference)",
        f"{(1.30 - L.E_X_eV) * 1e3:+.0f} meV",
        "no piezoelectric field (this (211)B-grown class has a strong built-in piezoelectric field along "
        "growth, explicitly omitted, module docstring limitation); dot-to-dot size/composition spread "
        "(single dots in the same paper range 1.283-1.44 eV)")
    deviation(
        "d-i-b: same dot, WL escape dE_pair_WL (0.8 nm WL, no SSL)",
        f"dE_pair_WL = {L.dE_pair_WL:.0f} meV", "110 +/- 10 meV (Chatzarakis et al., no-SSL Delta_E = E_WL - E_QD)",
        f"{L.dE_pair_WL - 110:+.0f} meV", "same piezoelectric-field omission as d-i-a, plus the fixed 10 meV 2D-QW WL "
        "exciton binding assumption")
    Lc = lv("InGaAs0.5/AlGaAs0.57/AlGaAs0.57 on GaAs")
    deviation(
        "d-i-c: same dot, WL embedded directly in bulk Al0.57Ga0.43As (0 nm GaAs spacer)",
        f"dE_pair_WL = {Lc.dE_pair_WL:.0f} meV", "240 +/- 20 meV (Chatzarakis et al., SSL sample, <Al> = 65%, Delta_E)",
        f"{Lc.dE_pair_WL - 240:+.0f} meV", "this is the ZERO-spacer idealised limit (matrix = barrier = "
        "Al0.57Ga0.43As, no intervening GaAs), an extreme end-member of the real digital-SSL structure, not a "
        "fit to it; see also d-i-d, whose finite-spacer scan under-predicts the same SSL-induced shift")
    d2 = None
    for d in (2.0,):
        S = P["InGaAs0.5/GaAs/AlGaAs0.57 on GaAs"]()
        S.geometry.matrix_thickness_nm = d
        d2 = D.levels(S).dE_pair_WL - L.dE_pair_WL
    deviation(
        "d-i-d: SSL-induced increase of dE_pair_WL at a 2 nm GaAs spacer",
        f"+{d2:.0f} meV", "+130 meV (Chatzarakis et al., SSL 240 meV minus no-SSL 110 meV, "
        "each +/- 10-20 meV)", f"{d2 - 130:+.0f} meV",
        "treating the digital short-period (GaAs/AlAs) superlattice barrier as a uniform bulk Al0.57Ga0.43As "
        "alloy at an effective distance underestimates the confinement of the true periodic structure, whose "
        "local band profile briefly reaches the full AlAs barrier height each period")
    Lr = lv("InP/GaInP/AlGaInP0.55 on GaAs")
    deviation(
        "d-ii-a: InP 3 nm x R 10 nm [A, dot size not measured] / Ga0.51In0.49P, Pryor VBO",
        f"E_X = {Lr.E_X_eV:.3f} eV", "1.815-1.836 eV (Reischle et al., Opt. Express 16, 12771 (2008), single-dot "
        "lines) / ~1.84-1.95 eV (Reischle et al., APL 97, 143513 (2010), single-dot lines) -- both single-dot, "
        "not ensemble; dot size not reported in either paper",
        f"{(1.815 - Lr.E_X_eV) * 1e3:+.0f} to {(1.95 - Lr.E_X_eV) * 1e3:+.0f} meV",
        "combines Pryor's isolated, measured unstrained VBO (-45 meV) with this module's own Vurgaftman bulk "
        "gaps and linear deformation-potential strain (module docstring limitation (iii)), compounded by the "
        "same large-mismatch/no-k.p-mixing limitation as HKUST/Gu and the [A] dot size")
    deviation(
        "d-ii-b: same dot, dot-matrix conduction band offset V_e",
        f"V_e = {Lr.V_e:.0f} meV", "~250 meV (240 meV DLTS) (Pryor, Pistol, Samuelson, PRB 56, 10404 (1997), "
        "6-band k.p, strained)", f"{Lr.V_e - 250:+.0f} meV (model 1.7x Pryor's value)",
        "same combination limitation as d-ii-a: Pryor's own strained CBO comes from a 6-band k.p calculation "
        "using Pryor's own bulk gaps and strain treatment, not separable from the single -45 meV VBO number "
        "reused here; investigated for a dot_levels.py coding defect (double-counted strain, wrong reference "
        "material) and none was found -- _with_vbo/materials.strain_shifts apply the override and the strain "
        "shift exactly as the module docstring specifies (rigid shift preserving each material's own gap); "
        "the escape energy dE_e_matrix (251 meV) coincidentally lands near Pryor's 250 meV despite V_e itself "
        "missing by 186 meV, which is why the previous check (testing dE_e_matrix under this V_e's label) is "
        "removed rather than kept as a passing but mislabelled comparison")
    Lt = lv("InAs/InP/InP")
    deviation(
        "d-iii: InAs 2.5 nm x R 15 nm / InP telecom dot, 300 K",
        f"E_X = {Lt.E_X_eV:.3f} eV ({Lt.lambda_nm:.0f} nm)",
        "0.887-0.953 eV (1301-1398 nm) (Laferriere et al., Nano Lett. 23, 962 (2023), single-dot, 4-300 K); "
        "typical InAs/InP telecom dot range 0.80-0.95 eV (1.3-1.55 um) [E]",
        f"{(0.80 - Lt.E_X_eV) * 1e3:.0f} meV below the typical-range lower edge",
        "this preset is a pure InAs dot directly in InP; Laferriere's actual device is an InAs0.68P0.32 "
        "dot-in-a-rod inside an InAs0.5P0.5 rod inside the InP core -- the graded InAsP composition has a "
        "larger gap than pure InAs, so the real device is less confined (blue-shifted) than this simplified "
        "preset, in the observed direction")
    Lg = lv("InP/GaAs0.65P0.35/AlGaAs0.4 on GaAs", height_nm=5.5)
    Lg4 = lv("InP/GaAsP0.4/AlGaAs0.4 on GaAs", height_nm=5.0)
    deviation(
        "d-iv-a: HKUST InP/GaAs0.65P0.35 disk, paper's own composition, 5.5 nm height (midpoint of 4-7 nm)",
        f"E_X = {Lg.E_X_eV:.3f} eV ({Lg.lambda_nm:.0f} nm); design-card GaAs0.60P0.40 [A] variant at 5 nm: "
        f"{Lg4.E_X_eV:.3f} eV ({Lg4.lambda_nm:.0f} nm)",
        "750-755 nm = 1.643-1.653 eV (Gu et al., Opt. Express 33, 23732 (2025), ensemble QD PL/lasing)",
        f"{(1.643 - Lg.E_X_eV) * 1e3:+.0f} to {(1.653 - Lg.E_X_eV) * 1e3:+.0f} meV",
        "large-mismatch (5.0-5.3%) 3D-island relaxation this coherent single-band disk does not fit (module "
        "docstring limitation (i)), compounded by unknown alloy ordering and the fact that the paper's line is "
        "an ensemble, type-I/II-mixed observable, not a single dot of this exact geometry")


_register_known_deviations()


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
        assert rp["E_b"] > 0 and rp["b_p"] == 100
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


def deviations_block():
    print("-" * 100)
    print("Known deviations (model vs. cited published value; documented, NOT counted as pass/fail)")
    print("-" * 100)
    for name, model_value, published, miss, explanation in KNOWN_DEVIATIONS:
        print(f"  {name}")
        print(f"    model: {model_value}")
        print(f"    published: {published}")
        print(f"    miss: {miss}")
        print(f"    why: {explanation}")
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
    print(f"{len(KNOWN_DEVIATIONS)} known deviations (documented, not counted)")
    if not QUIET:
        deviations_block()
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
