"""Independent algebra and policy checks for nitride orientation factors.

Source transcription is separate from the assumed orientation policy.  The
normal-polarization factors are deliberately [A], while electrostatic source
constants are [V] Bernardini et al., PRB 56, R10024 (1997), Table II and
Bernardini & Fiorentini, phys. status solidi (b) 216, 391 (1999), Sec. III.
Schade et al., phys. status solidi (b) (2011), doi:10.1002/pssb.201046350,
motivates the caveat that this is not a rotated-band-structure model.
"""
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from fsim_core.nitride_materials import EPS0_SI, binary, ingaN, orientation_factor, polarization_field

checks = []
def check(label, ok):
    checks.append(bool(ok))
    if not ok:
        print("FAIL", label)

def close(a, b, rel=1e-12, abs_tol=1e-12):
    return math.isclose(a, b, rel_tol=rel, abs_tol=abs_tol)

g = binary("GaN")
i = binary("InN")
check("source transcription B97 spontaneous polarization", g.Psp_Cm2 == -0.029 and i.Psp_Cm2 == -0.032)
check("source transcription BF99 GaN dielectric", g.eps_r == 10.28)

# Independently typed zero-strain VCA constants.  This isolates spontaneous
# polarization from piezoelectric strain and checks SI-to-kV/cm conversion.
x = 0.25
dot = ingaN(x)
p_dot_zero = (1.0-x)*(-0.029) + x*(-0.032)
f_zero_hand = ((-0.029)-p_dot_zero)/(8.8541878128e-12*((1.0-x)*10.28+x*14.61))*1e-5
check("analytical zero-strain alloy field", close(polarization_field(dot, g, 300, strain_fraction=0.0), f_zero_hand, rel=1e-10))
check("zero-strain spontaneous term retained", f_zero_hand != 0.0 and polarization_field(dot, g, 300, strain_fraction=0.0) != 0.0)

print("source transcription checks complete")
check("assumed orientation policy", orientation_factor("c_plane") == 1.0 and orientation_factor("semipolar_11_22") == 0.2 and orientation_factor("m_plane") == 0.0 and orientation_factor("a_plane") == 0.0)
check("semipolar sensitivity factors", orientation_factor("semipolar_11_22", 0.1) == 0.1 and orientation_factor("semipolar_11_22", 0.2) == 0.2 and orientation_factor("semipolar_11_22", 0.3) == 0.3)
check("explicit zero semipolar override", orientation_factor("semipolar_11_22", 0.0) == 0.0)

base = polarization_field(dot, g, 300, strain_fraction=0.0)
for screen in (0.0, 0.5, 1.0):
    for ext in (-100.0, 0.0, 100.0):
        intrinsic = (1.0-screen)*f_zero_hand
        c = polarization_field(dot, g, 300, strain_fraction=0.0, screening_fraction=screen, external_field_kVcm=ext)
        semi = polarization_field(dot, g, 300, strain_fraction=0.0, screening_fraction=screen, external_field_kVcm=ext, orientation="semipolar_11_22")
        m = polarization_field(dot, g, 300, strain_fraction=0.0, screening_fraction=screen, external_field_kVcm=ext, orientation="m_plane")
        a = polarization_field(dot, g, 300, strain_fraction=0.0, screening_fraction=screen, external_field_kVcm=ext, orientation="a_plane")
        check("algebra screen=%s external=%s" % (screen, ext), close(c, intrinsic+ext) and close(semi, 0.2*intrinsic+ext) and close(m, ext) and close(a, ext))
        check("nonpolar external retained screen=%s external=%s" % (screen, ext), m == ext and a == ext)

check("zero-field orientation ratios", close(polarization_field(dot,g,300,strain_fraction=0.0)/base,1.0) and close(polarization_field(dot,g,300,strain_fraction=0.0,orientation="semipolar_11_22")/base,0.2) and polarization_field(dot,g,300,strain_fraction=0.0,orientation="m_plane") == 0.0)
check("m and a model identity", polarization_field(dot,g,300,orientation="m_plane",external_field_kVcm=17.0) == polarization_field(dot,g,300,orientation="a_plane",external_field_kVcm=17.0))

# Legacy equality is deliberately exact, against the pre-change formula with
# independently spelled-out arithmetic, at varied composition and temperature.
for composition, temperature in ((0.15, 230.0), (0.25, 300.0), (0.4, 273.0)):
    d = ingaN(composition)
    # Repeat the documented VCA/strain expression without band_edges.
    a_dot = (1-composition)*3.189 + composition*3.545
    c13 = (1-composition)*106.0 + composition*121.0
    c33 = (1-composition)*398.0 + composition*182.0
    e31 = (1-composition)*(-0.49) + composition*(-0.57)
    e33 = (1-composition)*0.73 + composition*0.97
    psp = (1-composition)*(-0.029) + composition*(-0.032)
    ep = (3.189-a_dot)/a_dot
    ez = -2.0*c13/c33*ep
    p_dot = psp+2.0*e31*ep+e33*ez
    eps = (1-composition)*10.28+composition*14.61
    expected = (1.0-0.5)*((-0.029)-p_dot)/(8.8541878128e-12*eps)*1e-5+12.5
    actual = polarization_field(d,g,temperature,screening_fraction=0.5,external_field_kVcm=12.5)
    check("legacy exact composition=%s temperature=%s" % (composition, temperature), actual == expected)

bad = [
    lambda: orientation_factor("bad"),
    lambda: orientation_factor(True),
    lambda: orientation_factor([]),
    lambda: orientation_factor("semipolar_11_22", True),
    lambda: orientation_factor("semipolar_11_22", float("nan")),
    lambda: orientation_factor("semipolar_11_22", float("inf")),
    lambda: orientation_factor("semipolar_11_22", -0.1),
    lambda: orientation_factor("semipolar_11_22", 1.1),
    lambda: orientation_factor("c_plane", 0.2),
    lambda: orientation_factor("m_plane", 0.2),
    lambda: orientation_factor("a_plane", 1.0),
]
for n, fn in enumerate(bad):
    try:
        fn()
        ok = False
    except ValueError:
        ok = True
    check("reject malformed orientation input %d" % n, ok)

print("analytical limits and assumed orientation policy checks complete")
print("%d/%d nitride orientation checks passed" % (sum(checks), len(checks)))
raise SystemExit(0 if all(checks) else 1)
