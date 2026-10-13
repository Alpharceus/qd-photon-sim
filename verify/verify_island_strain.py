"""Verification of the opt-in buried-island strain model (r2-physics-b Q2, Fable 2026-10-12):
materials.eshelby_shear_factor, strain_shifts/offsets(shear_factor=...), DotGeometry.strain_model.

Run: python verify/verify_island_strain.py     (exit 0 iff all pass; prints N/N island_strain checks passed)

Independent references (none produced by the code under test):
  * the Eshelby oblate-spheroid shear factor is recomputed here from the DEFINING I-integrals
    (Mura, Micromechanics of Defects in Solids (1987), sec. 11) by scipy.quad, at two Poisson ratios,
    and compared with the production closed form;
  * the shape invariance of tr(eps) is checked against the closed form tr = 2 (1-2 nu)/(1-nu) eps_par
    (isotropic elasticity, film and sphere);
  * Chatzarakis et al., PRApplied 20, 034011 (2023): E_X (macro-PL band) 1.30 eV and Delta_E = E_WL - E_QD
    = 110 +/- 10 meV [V] are two observables closed by ONE shear factor s (the one-cause claim).

Also writes .workers/devlib/r2-rtedge-island-sensitivity.md (RT-edge presets under film vs island vs sphere).
The RT-edge default is NOT switched (user decision).
"""
import sys
from dataclasses import replace
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core import dot_levels as D  # noqa: E402
from fsim_core import materials as M  # noqa: E402
from fsim_core.integrator import retention  # noqa: E402

CHECKS = []


def check(name):
    def deco(fn):
        CHECKS.append((name, fn))
        return fn
    return deco


def s_quad(a, nu):
    """Shear factor from the defining integrals (independent of the closed form in materials.py)."""
    a1 = a2 = 1.0
    a3 = a
    pre = 2 * np.pi * a1 * a2 * a3

    def Delta(s):
        return np.sqrt((a1 * a1 + s) * (a2 * a2 + s) * (a3 * a3 + s))

    def I(i):
        return pre * quad(lambda s: 1 / ((i * i + s) * Delta(s)), 0, np.inf, epsabs=1e-13, epsrel=1e-12)[0]

    def Iij(i, j):
        return pre * quad(lambda s: 1 / ((i * i + s) * (j * j + s) * Delta(s)), 0, np.inf,
                          epsabs=1e-13, epsrel=1e-12)[0]

    I1, I3 = I(a1), I(a3)
    I11, I12, I13, I33 = Iij(a1, a1), Iij(a1, a2), Iij(a1, a3), Iij(a3, a3)
    c = 1 / (8 * np.pi * (1 - nu))
    d = (1 - 2 * nu) * c
    e11 = (3 * c * I11 + d * I1) + (c * I12 - d * I1) + (c * a3 ** 2 * I13 - d * I1)
    e33 = 2 * (c * I13 - d * I3) + (3 * c * a3 ** 2 * I33 + d * I3)
    return (e33 - e11) / ((1 + nu) / (1 - nu))


@check("1: eshelby_shear_factor == independent quadrature of Mura's I-integrals at a = 0.05..0.95 (nu 0.31 and 0.357), < 1e-6")
def _():
    for a in (0.05, 0.125, 0.25, 0.3, 0.5, 0.75, 0.9, 0.95):
        for nu in (0.31, 0.357):
            ref = s_quad(a, nu)
            got = M.eshelby_shear_factor(a)
            assert abs(ref - got) < 1e-6, (a, nu, ref, got)


@check("2: s is nu-independent (quadrature at nu 0.2 / 0.31 / 0.45 agree to 1e-8) and monotone decreasing 1 -> 0")
def _():
    for a in (0.25, 0.5):
        v = [s_quad(a, nu) for nu in (0.2, 0.31, 0.45)]
        assert max(v) - min(v) < 1e-8, v
    xs = np.linspace(0.02, 0.98, 25)
    ss = [M.eshelby_shear_factor(x) for x in xs]
    assert all(b < a_ for a_, b in zip(ss, ss[1:])), ss
    assert abs(M.eshelby_shear_factor(1e-9) - 1) < 1e-4 and abs(M.eshelby_shear_factor(1 - 1e-9)) < 1e-6


@check("3: published-probe anchors s(0.25) = 0.555 and s(0.5) = 0.291 (+/- 0.002), s(0.125) = 0.746, s(0.75) = 0.119 (+/- 0.002)")
def _():
    for a, ref in ((0.125, 0.746), (0.25, 0.555), (0.5, 0.291), (0.75, 0.119)):
        assert abs(M.eshelby_shear_factor(a) - ref) <= 0.002, (a, M.eshelby_shear_factor(a))


@check("4: tr(eps) is shape-invariant: strain_shifts(shear_factor=s).tr_eps is s-independent and equals 2 (1-2 nu)/(1-nu) eps_par")
def _():
    for lay, sub in ((M.binary("InP"), M.binary("GaAs")), (M.InGaAs(0.5), M.binary("GaAs")),
                     (M.binary("InAs"), M.binary("InP"))):
        nu = lay.p["C12"] / (lay.p["C11"] + lay.p["C12"])
        e0 = M.strain_shifts(lay, sub, 5.0)
        for s in (0.0, 0.3, 1.0):
            e = M.strain_shifts(lay, sub, 5.0, shear_factor=s)
            assert e.tr_eps == e0.tr_eps and e.eps_par == e0.eps_par and e.E_c == e0.E_c
        assert abs(e0.tr_eps - 2 * (1 - 2 * nu) / (1 - nu) * e0.eps_par) < 1e-12, (lay.label, e0.tr_eps)


@check("5: shear_factor scales exactly Q: E_hh(s) = E_hh(1) + (1-s) Q_film (closed form), E_hh moves DOWN as s falls under compression")
def _():
    lay, sub = M.binary("InP"), M.binary("GaAs")
    e1 = M.strain_shifts(lay, sub, 5.0)
    Q = lay.p["b"] * (e1.eps_perp - e1.eps_par)
    assert Q < 0
    for s in (0.0, 0.291, 0.555):
        e = M.strain_shifts(lay, sub, 5.0, shear_factor=s)
        assert abs(e.E_hh - (e1.E_hh + (1.0 - s) * Q)) < 1e-12, s
        assert e.E_hh < e1.E_hh


def _scan_s(system, s):
    """DotLevels with a forced dot shear factor s (verifier-local override of the model's chooser)."""
    orig = D._shear_factor
    D._shear_factor = lambda g: s
    try:
        return D.levels(system)
    finally:
        D._shear_factor = orig


@check("6: ONE-CAUSE closure (Chatzarakis): the s that puts E_X at the macro-PL 1.30 eV also puts dE_pair_WL "
       "in 110 +/- 10 meV (two observables, one parameter)")
def _():
    S = D.class_presets()["InGaAs0.5/GaAs/GaAs on GaAs"]()
    s_star = brentq(lambda s: _scan_s(S, s).E_X_eV - 1.30, 0.0, 1.0, xtol=1e-6)
    L = _scan_s(S, s_star)
    assert 0.0 < s_star < 1.0, s_star
    assert abs(L.dE_pair_WL - 110.0) <= 10.0, (s_star, L.dE_pair_WL)
    # and the production island model for the preset's own h/R lies between the film and sphere limits
    f, sph = D.levels_bracket(S)
    assert f.E_X_eV < 1.30 < sph.E_X_eV and sph.dE_pair_WL < 110.0 < f.dE_pair_WL


@check("7: levels_bracket contains the Chatzarakis macro-PL E_X (1.30 in [film, sphere]) and Delta_E = 110 in [sphere, film]; "
       "Gu ensemble line 1.648 is NOT contained (documented residual)")
def _():
    f, sph = D.levels_bracket(D.class_presets()["InGaAs0.5/GaAs/GaAs on GaAs"]())
    assert f.E_X_eV <= 1.30 <= sph.E_X_eV
    assert sph.dE_pair_WL <= 110.0 <= f.dE_pair_WL
    fg, sg = D.levels_bracket(D.class_presets()["InP/GaAs0.65P0.35/AlGaAs0.4 on GaAs"](height_nm=5.5))
    assert not (fg.E_X_eV <= 1.648 <= sg.E_X_eV), (fg.E_X_eV, sg.E_X_eV)


@check("8: Hole side: V_h falls with s; Reischle (Pryor VBO) hole becomes unbound (type II) in the sphere limit "
       "(Pryor: holes confined by residual shear only)")
def _():
    f, sph = D.levels_bracket(D.class_presets()["InP/GaInP/AlGaInP0.55 on GaAs"]())
    assert f.V_h > 0 > sph.V_h and f.hole_bound and not sph.hole_bound


@check("9: bad strain_model raises; sensitivity file written")
def _():
    S = D.class_presets()["InAs/InP/InP"]()
    try:
        D.levels(replace(S, geometry=replace(S.geometry, strain_model="pyramid")))
    except ValueError:
        pass
    else:
        raise AssertionError("no ValueError")
    assert write_sensitivity().exists()


def write_sensitivity():
    """RT-edge sensitivity row: Reischle / Bommer InP dot chains at 300 K, film vs island vs sphere."""
    out = ROOT / ".workers" / "devlib" / "r2-rtedge-island-sensitivity.md"
    out.parent.mkdir(parents=True, exist_ok=True)
    P = D.class_presets()
    rows = []
    for name in ("InP/GaInP/AlGaInP0.55 on GaAs", "InP/AlGaInP0.2/AlGaInP0.55 on GaAs"):
        S = P[name](T=300.0)
        for mdl in ("film", "island", "sphere"):
            L = D.levels(replace(S, geometry=replace(S.geometry, strain_model=mdl)))
            rp = D.retention_params(L, 1.0, "pair_half", verbose=False)
            Sret = float(retention(300.0, rp["a_esc"], rp["E_a"], rp["b_p"], rp["E_b"]))
            rows.append((name, mdl, L.shear_factor, L.V_e, L.V_h, L.dE_e_matrix, L.dE_h_matrix, L.dE_pair_WL,
                         L.E_X_eV, rp["E_a"], Sret, L.type))
    lines = [
        "# R2 item 10: RT-edge island-shear sensitivity (strain_model film vs island vs sphere), 300 K",
        "",
        "Generated by verify/verify_island_strain.py (dot_levels.levels + retention_params(tau_rad = 1 ns, "
        "channel pair_half)). Presets: reischle (Pryor VBO -45 meV, h 3 nm R 10 nm) and bommer (transferred "
        "VBO). The RT-edge default stays strain_model='film'; switching it is a USER DECISION because it "
        "collapses the hole escape energy (retention input) -- the S(300 K) column is the Module-E retention "
        "of that single 1-ns-radiative dot.",
        "",
        "| preset | model | s | V_e meV | V_h meV | dE_e_matrix | dE_h_matrix | dE_pair_WL | E_X eV | E_a (pair_half) meV | S(300 K) | type |",
        "|---|---|---|---|---|---|---|---|---|---|---|---|",
    ]
    for r in rows:
        lines.append("| {} | {} | {:.3f} | {:.0f} | {:.0f} | {:.0f} | {:.0f} | {:.0f} | {:.3f} | {:.0f} | {:.3g} | {} |".format(*r))
    lines += ["", "Pins that would move on a default switch: verify_device_rt (legacy npz), rt_edge sweeps, "
              "verify_rt_edge_papers reischle10 row. None are touched by this opt-in change.", ""]
    out.write_text("\n".join(lines), encoding="utf-8")
    return out


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
    print(f"\n{n - failed}/{n} island_strain checks passed")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
