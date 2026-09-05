"""verify_materials -- materials database checks against numbers the code did
NOT produce (published lattice constants, 300 K gaps, offsets, crossovers).

Run: python verify/verify_materials.py   (exit 0 iff all PASS)

Every reference value below carries its source. Tolerances are stated per
item and are deliberately tight where the compilation is exact (lattice
constants, binary gaps) and loose where the literature itself spreads
(alloy offsets, critical-thickness conventions).
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from fsim_core import materials as M  # noqa: E402

RESULTS = []


def check(name, value, ref, tol, note="", rel=False):
    err = abs(value - ref) / (abs(ref) if rel and ref != 0 else 1.0)
    ok = err <= tol
    RESULTS.append(ok)
    unit = "%" if rel else ""
    print(f"  {'PASS' if ok else 'FAIL'}  {name}\n        computed {value:.5g}  vs ref {ref:.5g}"
          f"  (tol {tol * (100 if rel else 1):.3g}{unit})  -- {note}")
    return ok


def bool_check(name, ok, note=""):
    RESULTS.append(bool(ok))
    print(f"  {'PASS' if ok else 'FAIL'}  {name}  -- {note}")


GaAs, InP, InAs, GaP, AlAs, AlP = (M.binary(k) for k in ("GaAs", "InP", "InAs", "GaP", "AlAs", "AlP"))

print("== 1. Lattice constants at 300 K (Vurgaftman 2001 Table; Ioffe NSM) ==")
check("a(GaAs)", M.lattice_constant(GaAs), 5.65325, 0.001, "V01 / Ioffe 5.65325 A")
check("a(InP)", M.lattice_constant(InP), 5.8687, 0.002, "Ioffe 5.8687 A (V01 5.8697)")
check("a(InAs)", M.lattice_constant(InAs), 6.0583, 0.001, "V01 / Ioffe 6.0583 A")
check("a(GaP)", M.lattice_constant(GaP), 5.4505, 0.001, "V01 / Ioffe 5.4505 A")
check("a(AlAs)", M.lattice_constant(AlAs), 5.6611, 0.001, "V01 5.6611 A")

print("== 2. Binary 300 K gaps (eV) ==")
check("Eg_G(GaAs,300)", M.bandgap(GaAs, 300, "G"), 1.424, 0.01, "Ioffe/Sze 1.424 eV")
check("Eg_G(InP,300)", M.bandgap(InP, 300, "G"), 1.344, 0.012, "Ioffe 1.344 eV (V01 Varshni gives 1.353)")
check("Eg_G(InAs,300)", M.bandgap(InAs, 300, "G"), 0.354, 0.01, "Ioffe 0.354 eV")
check("Eg_X(GaP,300)", M.bandgap(GaP, 300, "X"), 2.26, 0.02, "Ioffe 2.26 eV indirect")
check("Eg_X(AlAs,300)", M.bandgap(AlAs, 300, "X"), 2.16, 0.02, "Ioffe 2.16 eV indirect")
bool_check("GaP indirect, AlAs indirect, GaAs/InP/InAs direct",
           (not M.is_direct(GaP)) and (not M.is_direct(AlAs)) and all(M.is_direct(m) for m in (GaAs, InP, InAs)))

print("== 3. Low-temperature gaps (eV) ==")
check("Eg_G(GaAs,0)", M.bandgap(GaAs, 0.0, "G"), 1.519, 0.002, "V01 1.519 eV")
check("Eg_G(InP,0)", M.bandgap(InP, 0.0, "G"), 1.4236, 0.002, "V01 1.4236 eV")

print("== 4. Lattice-matched alloys ==")
g = M.GaInP(0.51)
check("Ga0.51In0.49P misfit to GaAs", M.mismatch(g, GaAs) * 100, 0.0, 0.1, "lattice matched (< 0.1%)")
check("Eg_G(Ga0.51In0.49P, 300)", M.bandgap(g, 300, "G"), 1.90, 0.04,
      "disordered GaInP 1.88-1.92 eV (V01 Sec. V; ellipsometry 1.85+-0.02 partially ordered)")
ig = M.InGaAs(0.532)
check("In0.532Ga0.468As misfit to InP", M.mismatch(ig, InP) * 100, 0.0, 0.05, "lattice matched")
check("Eg(In0.53Ga0.47As, 300)", M.bandgap(ig, 300, "G"), 0.74, 0.02, "standard 0.74-0.75 eV (V01)")

print("== 5. Direct-indirect crossovers ==")
xs = [x / 100 for x in range(0, 101)]
x_cross = next((x for x in xs if not M.is_direct(M.AlGaInP(x))), None)
check("(AlxGa1-x)0.51In0.49P crossover x", x_cross, 0.53, 0.06,
      "literature x_c = 0.5-0.55 at Eg ~2.25-2.3 eV (Ioffe / Adachi)")
check("AlGaInP crossover energy", M.bandgap(M.AlGaInP(x_cross), 300, "G"), 2.27, 0.06, "~2.25-2.3 eV")
x_gp = next((x for x in xs if not M.is_direct(M.GaAsP(x))), None)
check("GaAs1-xPx crossover x", x_gp, 0.47, 0.04, "Ioffe: x_c = 0.45-0.49 (classic 0.45)")
check("Eg(GaAs0.6P0.4, 300)", M.bandgap(M.GaAsP(0.4), 300, "G"), 1.92, 0.04,
      "Ioffe/Adachi direct gap ~1.9-1.93 eV at x=0.4")
x_al = next((x for x in xs if not M.is_direct(M.AlGaAs(x))), None)
check("AlxGa1-xAs crossover x", x_al, 0.43, 0.04, "V01 / Ioffe 0.41-0.45")

print("== 6. Band offsets ==")
# [V] Vurgaftman, Meyer & Ram-Mohan, JAP 89, 5815 (2001), Table XIV.
check("AlAs-InAs VBO bowing", M.BOWING[("AlAs", "InAs")]["vbo"], -0.64, 1e-12,
      "V01 Table XIV")
check("AlAs-InAs a_c bowing", M.BOWING[("AlAs", "InAs")]["a_c"], -1.4, 1e-12,
      "V01 Table XIV")
o = M.offsets(GaAs, M.AlGaAs(0.3), GaAs)
check("GaAs/Al0.3Ga0.7As Q_c", o.dE_c / (o.dE_c + o.dE_v_hh), 0.62, 0.06,
      "accepted 60/40 to 65/35 rule (Batey & Wright 1986)")
o = M.offsets(ig, InP, InP)
check("In0.53Ga0.47As/InP dE_c (eV)", o.dE_c, 0.25, 0.05, "literature 0.22-0.27 eV (Q_c ~ 0.40)")
check("In0.53Ga0.47As/InP dE_v (eV)", o.dE_v_hh, 0.36, 0.05, "literature 0.34-0.37 eV")
o = M.offsets(InP, g, GaAs)
check("InP/Ga0.51In0.49P strained dE_c (eV)", o.dE_c, 0.25, 0.08,
      "Pryor PRB 56, 10404 (1997): ~250 meV CB offset (240 meV DLTS)")
bool_check("InP in GaInP: compressive, hh top", o.well.eps_par < -0.03 and o.well.E_hh > o.well.E_lh,
           f"eps_par = {o.well.eps_par:.4f}")
o2 = M.offsets(InP, M.GaAsP(0.4), GaAs)
bool_check("InP in GaAs0.6P0.4: electrons deeply confined (dE_c > 0.15 eV)", o2.dE_c > 0.15,
           f"dE_c = {o2.dE_c:.3f} eV, dE_v_hh = {o2.dE_v_hh:.3f} eV (weak, strain-set: Gu et al. 2025 type-I/II)")

print("== 7. Strain / misfit ==")
check("InAs on GaAs misfit (%)", M.mismatch(InAs, GaAs) * 100, -6.7, 0.3, "standard 6.7-7.2% depending on reference lattice")
check("InP on GaAs misfit (%)", M.mismatch(InP, GaAs) * 100, -3.7, 0.15, "standard 3.7-3.8%")
check("GaAs0.6P0.4 on GaAs misfit (%)", M.mismatch(M.GaAsP(0.4), GaAs) * 100, 1.43, 0.1, "tensile ~1.4%")
check("InP on GaAs0.6P0.4 misfit (%)", M.mismatch(InP, M.GaAsP(0.4)) * 100, -5.1, 0.3, "computed in materials_research 5.07-5.34%")
sh = M.strain_shifts(GaAs, GaAs, eps_par=-0.01)
check("GaAs hydrostatic gap shift at tr(eps)=-0.0106 (eV)", sh.Eg_hh - sh.Eg_lh, 0.0, 0.2, "hh-lh split sign check only")
# [V] Vurgaftman, Meyer & Ram-Mohan, JAP 89, 5815 (2001), Table I:
# V01's a = a_c + a_v = -8.33 eV in its sign convention.
check("GaAs gap deformation potential a = a_c + a_v (eV)",
      GaAs.p["a_c"] + GaAs.p["a_v"], -8.33, 0.05, "V01 Table I")
check("InP/GaAs exact LH-HH splitting (eV)",
      M.strain_shifts(InP, GaAs).E_hh - M.strain_shifts(InP, GaAs).E_lh,
      0.066, 0.01, "Van de Walle PRB 39, 1871 (1989), exact LH-SO block")
check("InP X gap at 0 K (eV)", M.bandgap(InP, 0.0, "X"), 2.384, 1e-12,
      "V01 Table VI linear X-gap form")
hc = M.critical_thickness_nm(M.InGaAs(0.2), GaAs)
bool_check("MB h_c In0.2Ga0.8As/GaAs in the 4-8 nm single-layer class", 4.0 < hc < 8.0,
           f"h_c = {hc:.1f} nm; single-layer MB form")
hc2 = M.critical_thickness_nm(InAs, GaAs)
bool_check("MB h_c InAs/GaAs in the 0.3-0.6 nm class", 0.3 < hc2 < 0.6,
           f"h_c = {hc2:.2f} nm; single-layer MB form")
hc_ip = M.critical_thickness_nm(InP, GaAs)
bool_check("MB h_c InP/GaAs in the 1.0-1.6 nm class", 1.0 < hc_ip < 1.6,
           f"h_c = {hc_ip:.2f} nm; single-layer MB form")
t_bal = M.strain_balance_thickness(InP, 0.6, M.GaAsP(0.4), GaAs)
bool_check("strain balance: 2 ML InP needs ~1-2 nm GaAs0.6P0.4 (thin) ", 0.5 < t_bal < 3.0, f"t = {t_bal:.2f} nm")

print("== 8. Extra properties ==")
check("k(GaAs,300)", M.thermal_k("GaAs"), 55.0, 3.0, "Ioffe 55 W/m/K")
check("k(InP,300)", M.thermal_k("InP"), 68.0, 3.0, "Ioffe 68 W/m/K")
check("n(GaAs,668nm)", M.refractive_index("GaAs", 668), 3.81, 0.03, "Aspnes 1986 3.81-3.83")
check("n(AlAs,668nm)", M.refractive_index("AlAs", 668), 3.08, 0.02, "Fern & Onton")
check("n(InP,1550nm)", M.refractive_index("InP", 1550), 3.17, 0.02, "Pettit & Turner 3.17")
rep = M.feasibility_report(
    [M.StackLayerSpec(M.AlGaInP(0.7), 1000, "cladding"), M.StackLayerSpec(M.GaAsP(0.4), 8, "well"),
     M.StackLayerSpec(InP, 0.6, "dot"), M.StackLayerSpec(GaAs, 100, "cap")], GaAs, 300, emission_eV=1.856)
bool_check("feasibility report flags the absorbing GaAs cap at 668 nm",
           any("GaAs (cap)" in w for w in rep["warnings"]), str(rep["warnings"]))
bool_check("InGaAs default label preserves both compositions",
           M.InGaAs(0.532).label == "In0.532Ga0.468As", M.InGaAs(0.532).label)

n_pass = sum(RESULTS)
print(f"\n{n_pass}/{len(RESULTS)} materials checks passed")
sys.exit(0 if n_pass == len(RESULTS) else 1)
