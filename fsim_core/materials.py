"""Module M -- III-V material database, alloy interpolation, strain and band
alignment (RT edge-emitter tier, 2026-09-02).

WHY THIS MODULE EXISTS. Before this tier the simulator carried no bandgaps,
lattice constants, band offsets or refractive indices at all (code audit
2026-09-02, section 3): it could not say whether a proposed stack is even
growable, and every thermal-escape energy was an arsenide-class fit proxy.
This module makes those numbers explicit, tagged, and testable.

DATA PROVENANCE. Binary band parameters (bandgaps, Varshni coefficients,
valence-band offsets, deformation potentials, elastic constants, effective
masses, Luttinger parameters, spin-orbit splittings) are the recommended
values of Vurgaftman, Meyer and Ram-Mohan, J. Appl. Phys. 89, 5815 (2001)
["V01"], Tables I-XII, with the ternary bowing parameters from the same
compilation. Thermal conductivities and refractive indices are collected
per material with their own citations in MATERIAL_EXTRA. Every value is
tagged [DR] (derived from a published compilation) unless stated; values
tagged [E] are estimates and are called out in the field comment.

SIGN CONVENTIONS (stated once, used everywhere):
  * Strain: eps_par = (a_sub - a_layer)/a_layer for a pseudomorphic layer on a
    substrate of in-plane lattice constant a_sub (negative = compressive).
    eps_perp = -2 (C12/C11) eps_par. tr(eps) = 2 eps_par + eps_perp.
  * Conduction band: dE_c = a_c tr(eps)   (a_c < 0: compression raises E_c).
  * Valence band average: dE_v,av = a_v_VdW tr(eps) with a_v_VdW = -a_v(V01)
    > 0, i.e. the Van de Walle (PRB 39, 1871 (1989)) convention in which the
    average valence band moves DOWN under compression. The band-gap
    deformation potential is a = a_c - a_v_VdW (GaAs: -7.17 - 1.16 = -8.33 eV,
    the accepted value), which fixes the convention unambiguously.
  * Biaxial (shear) splitting, Chuang "Physics of Optoelectronic Devices"
    ch. 4: Q_eps = b (eps_perp - eps_par) (b < 0 so Q_eps < 0 under
    compression); E_hh = E_v,av + Delta_so/3 - Q_eps  ... we quote band edges
    relative to the unstrained VB top of the same material, so
        E_hh(strained) = E_v0 + a_v_VdW tr(eps) - Q_eps
        E_lh(strained) is evaluated from the exact k=0 LH-SO block below.
    Under compression the heavy hole is the top valence band.
  * Absolute band alignment: E_v0 = VBO (V01 Table, relative to InSb VB = 0),
    E_c0 = E_v0 + E_g^Gamma(T). Alloy VBOs interpolate linearly (V01's own
    recommendation absent a bowing entry, except InGaAs and AlInAs where V01
    gives bowing).

ALLOY INTERPOLATION. Ternary A_xB_(1-x)C: P(x) = x P_A + (1-x) P_B
- x(1-x) C_bow. Quaternaries of the (A_xB_(1-x))_yC_(1-y)D type (AlGaInP,
AlGaInAs on InP) use a weighted average of the two constituent ternaries,
with the remaining A-B pair bowing applied on the cation sublattice. Temperature enters only through the Varshni bandgap
and the lattice-constant expansion coefficient.

WHAT THIS MODULE DOES NOT DO: no k.p band mixing, no alloy ordering (CuPt
ordering in GaInP lowers E_g by up to ~100 meV; we quote the DISORDERED
alloy, which is what MOVPE on misoriented substrates delivers -- flagged in
`ordering_note`), no piezoelectric fields, no temperature dependence of the
deformation potentials.

Every closed form in this module has an independent check in
verify/verify_materials.py (against published numbers the code did not
produce: lattice constants, 300 K gaps, GaAs/AlAs and InP/InGaAs offsets,
the AlGaInP direct-indirect crossover, the classic InP/GaInP mismatch).
"""
from __future__ import annotations

import re
from dataclasses import dataclass, field

import numpy as np

KB_EV = 8.617333262e-5  # eV/K

# ----------------------------------------------------------------------------- binaries

@dataclass(frozen=True)
class Binary:
    """Binary III-V band parameters, V01 recommended values.

    a300      lattice constant at 300 K [A];  da_dT [A/K]
    Eg_G/X/L  bandgap at 0 K for each valley [eV], with Varshni (alpha [eV/K],
              beta [K]) per valley.  GaP's Gamma gap uses V01's Bose form
              (handled in `bandgap`).
    vbo       valence-band offset, V01 Table (InSb VB = 0) [eV]
    d_so      spin-orbit splitting [eV]
    m_e       Gamma-valley electron mass [m0]
    g1,g2,g3  Luttinger parameters
    a_c, a_v  hydrostatic deformation potentials, V01 sign (a_v < 0 in V01's
              table; we convert to Van de Walle sign in `strain_shifts`) [eV]
    b         shear deformation potential [eV]
    C11, C12  elastic constants [GPa]
    """
    name: str
    a300: float
    da_dT: float
    Eg_G: float
    alpha_G: float
    beta_G: float
    Eg_X: float
    alpha_X: float
    beta_X: float
    Eg_L: float
    alpha_L: float
    beta_L: float
    vbo: float
    d_so: float
    m_e: float
    g1: float
    g2: float
    g3: float
    a_c: float
    a_v: float
    b: float
    C11: float
    C12: float
    source: str = "Vurgaftman, Meyer, Ram-Mohan, JAP 89, 5815 (2001)"
    tag: str = "DR"


BINARIES = {
    "GaAs": Binary("GaAs", 5.65325, 3.88e-5,
                   1.519, 0.5405e-3, 204.0,
                   1.981, 0.460e-3, 204.0,
                   1.815, 0.605e-3, 204.0,
                   -0.80, 0.341, 0.067, 6.98, 2.06, 2.93,
                   -7.17, -1.16, -2.0, 1221.0, 566.0),
    "AlAs": Binary("AlAs", 5.6611, 2.90e-5,
                   3.099, 0.885e-3, 530.0,
                   2.24, 0.70e-3, 530.0,
                   2.46, 0.605e-3, 204.0,
                   -1.33, 0.28, 0.15, 3.76, 0.82, 1.42,
                   -5.64, -2.47, -2.3, 1250.0, 534.0),
    "InAs": Binary("InAs", 6.0583, 2.74e-5,
                   0.417, 0.276e-3, 93.0,
                   1.433, 0.276e-3, 93.0,
                   1.133, 0.276e-3, 93.0,
                   -0.59, 0.39, 0.026, 20.0, 8.5, 9.2,
                   -5.08, -1.00, -1.8, 832.9, 452.6),
    "InP": Binary("InP", 5.8697, 2.79e-5,
                  1.4236, 0.363e-3, 162.0,
                  # V01: E_X(T) = 2.384 - 3.7e-4 T  (linear) -> encode as
                  # Varshni with beta = 0 and alpha = 3.7e-4 (exact for beta=0)
                  2.384, 0.37e-3, 0.0,
                  2.014, 0.363e-3, 162.0,
                  -0.94, 0.108, 0.0795, 5.08, 1.60, 2.10,
                  -6.0, -0.6, -2.0, 1011.0, 561.0),
    "GaP": Binary("GaP", 5.4505, 2.92e-5,
                  # Gamma gap: V01 Bose form 2.886 + 0.1081[1 - coth(164/T)];
                  # alpha/beta unused for GaP Gamma (see bandgap()).
                  2.886, 0.0, 0.0,
                  2.35, 0.5771e-3, 372.0,
                  2.72, 0.5771e-3, 372.0,
                  -1.27, 0.08, 0.13, 4.05, 0.49, 2.93,
                  -8.2, -1.7, -1.6, 1405.0, 620.3),
    "AlP": Binary("AlP", 5.4672, 2.92e-5,
                  3.63, 0.5771e-3, 372.0,
                  2.52, 0.318e-3, 588.0,
                  3.57, 0.318e-3, 588.0,
                  -1.74, 0.07, 0.22, 3.35, 0.71, 1.23,
                  -5.7, -3.0, -1.5, 1330.0, 630.0),
}

# Ternary bowing parameters (V01 Tables), keyed by the alloy's two binaries
# in the order (A, B) for A_x B_(1-x). Fields absent -> zero bowing.
# Bowing enters as P = x P_A + (1-x) P_B - x(1-x) C.
BOWING = {
    ("AlAs", "GaAs"): {"Eg_G": "algaas", "Eg_X": 0.055, "Eg_L": 0.0,
                       "note": "V01: Gamma bowing = -0.127 + 1.310 x (x = Al fraction)"},
    ("InAs", "GaAs"): {"Eg_G": 0.477, "Eg_X": 1.4, "Eg_L": 0.33, "d_so": 0.15,
                       "m_e": 0.0091, "vbo": -0.38, "a_c": 2.61},
    ("GaAs", "GaP"):  {"Eg_G": 0.19, "Eg_X": 0.24, "Eg_L": 0.16},
    ("InAs", "InP"):  {"Eg_G": 0.10, "Eg_X": 0.27, "Eg_L": 0.27, "d_so": 0.16},
    ("GaP", "InP"):   {"Eg_G": 0.65, "Eg_X": 0.20, "Eg_L": 1.03, "m_e": 0.051},
    ("AlP", "InP"):   {"Eg_G": -0.48, "Eg_X": 0.38, "d_so": -0.19, "m_e": 0.22},
    ("AlP", "GaP"):   {"Eg_G": 0.0, "Eg_X": 0.13},
    # [V] Vurgaftman, Meyer & Ram-Mohan, JAP 89, 5815 (2001), Table XIV.
    ("AlAs", "InAs"): {"Eg_G": 0.70, "Eg_X": 0.0, "Eg_L": 0.0, "d_so": 0.15,
                         "m_e": 0.049, "vbo": -0.64, "a_c": -1.4},
}

_BOWABLE = ("Eg_G", "Eg_X", "Eg_L", "d_so", "m_e", "vbo", "a_c")
_LINEAR = ("a300", "da_dT", "alpha_G", "beta_G", "alpha_X", "beta_X", "alpha_L",
           "beta_L", "g1", "g2", "g3", "a_v", "b", "C11", "C12")


# --------------------------------------------------------------------- alloy objects

@dataclass
class Material:
    """A binary, ternary or quaternary with all parameters resolved to numbers.
    `label` is a human name; `composition` records how it was built."""
    label: str
    p: dict
    composition: dict = field(default_factory=dict)
    tag: str = "DR"
    ordering_note: str = ""
    # temperature-dependent gap tree: binaries mixed at T (V01 procedure),
    # then ternary bowing terms subtracted: list of (pair, xa, weight)
    bow_terms: list = field(default_factory=list)

    def __getattr__(self, k):
        if k in ("p", "__setstate__", "__deepcopy__"):
            raise AttributeError(k)
        try:
            return self.p[k]
        except KeyError as e:
            raise AttributeError(k) from e


def _bin_dict(b: Binary) -> dict:
    return {k: getattr(b, k) for k in (_BOWABLE + _LINEAR)}


def binary(name: str) -> Material:
    b = BINARIES[name]
    return Material(name, _bin_dict(b), {name: 1.0}, b.tag)


def _binary_gap(name: str, T: float, valley: str) -> float:
    b = BINARIES[name]
    if name == "GaP" and valley == "G":
        T_ = max(float(T), 1.0)
        return b.Eg_G + 0.1081 * (1.0 - 1.0 / np.tanh(164.0 / T_))
    return float(_varshni(getattr(b, f"Eg_{valley}"), getattr(b, f"alpha_{valley}"),
                          getattr(b, f"beta_{valley}"), T))


def _bow(key, pair, x):
    C = BOWING.get(pair, {}).get(key, 0.0)
    if C == "algaas":  # V01 composition-dependent Gamma bowing for AlGaAs
        C = -0.127 + 1.310 * x
    return C


def ternary(A: str, B: str, x: float, label: str | None = None) -> Material:
    """A_x B_(1-x): e.g. ternary("GaP", "InP", 0.51) = Ga0.51In0.49P;
    ternary("InAs","GaAs",0.53) = In0.53Ga0.47As; ternary("AlAs","GaAs",0.45)."""
    if not 0.0 <= x <= 1.0:
        raise ValueError("composition x must be in [0,1]")
    pa, pb = _bin_dict(BINARIES[A]), _bin_dict(BINARIES[B])
    pair = (A, B) if (A, B) in BOWING else (B, A)
    xa = x if pair == (A, B) else 1.0 - x
    out = {}
    for k in _LINEAR:
        out[k] = x * pa[k] + (1 - x) * pb[k]
    for k in _BOWABLE:
        out[k] = x * pa[k] + (1 - x) * pb[k] - xa * (1 - xa) * _bow(k, pair, xa)
    lab = label or f"{A}{x:.2f}{B}{1.0 - x:.2f}"
    return Material(lab, out, {A: x, B: 1 - x}, "DR", bow_terms=[(pair, xa, 1.0)])


def quaternary_from_ternaries(terns: list, weights: list, label: str) -> Material:
    """Weighted average of resolved ternaries (V01 eq. 2.9 spirit). Used for
    (Al_x Ga_1-x)_y In_1-y P and (Al_x Ga_1-x)_y In_1-y As where the two
    Al/Ga ternaries with In are the natural constituents."""
    w = np.asarray(weights, float)
    w = w / w.sum()
    out = {}
    for k in (_BOWABLE + _LINEAR):
        out[k] = float(sum(wi * t.p[k] for wi, t in zip(w, terns)))
    comp = {}
    bows = []
    for wi, t in zip(w, terns):
        for b, f in t.composition.items():
            comp[b] = comp.get(b, 0.0) + wi * f
        bows += [(pair, xa, float(wi) * wt) for (pair, xa, wt) in t.bow_terms]
    return Material(label, out, comp, "DR", bow_terms=bows)


# ------------------------------------------------------------- named material builders

def GaInP(x_ga: float = 0.51) -> Material:
    """Ga_x In_1-x P (x = 0.51 lattice-matched to GaAs at 300 K, V01)."""
    m = ternary("GaP", "InP", x_ga, f"Ga{x_ga:.2f}In{1 - x_ga:.2f}P")
    m.ordering_note = ("disordered-alloy gap; CuPt-ordered GaInP is up to ~0.1 eV "
                       "lower (V01 Sec. V.C) -- MOVPE on 6-10 deg off (100) suppresses ordering")
    return m


def AlInP(x_al: float = 0.52) -> Material:
    return ternary("AlP", "InP", x_al, f"Al{x_al:.2f}In{1 - x_al:.2f}P")


def AlGaInP(x_al: float, y_III: float = 0.51) -> Material:
    """(Al_x Ga_1-x)_y In_1-y P; y = 0.51 lattice-matched to GaAs. Built as the
    x-weighted mix of the AlInP and GaInP ternaries at the same In content
    (exact for the linear parameters; the Al-Ga bowing on the P sublattice,
    V01 AlGaP Eg_X bowing 0.13, is applied as a small correction)."""
    t_al = AlInP(y_III)
    t_ga = GaInP(y_III)
    m = quaternary_from_ternaries([t_al, t_ga], [x_al, 1 - x_al],
                                  f"(Al{x_al:.2f}Ga{1 - x_al:.2f}){y_III:.2f}In{1 - y_III:.2f}P")
    # [V] Vurgaftman, Meyer & Ram-Mohan, JAP 89, 5815 (2001), Table XVII:
    # pairwise Al-Ga bowing is weighted by the squared P-sublattice fraction.
    m.bow_terms.append((("AlP", "GaP"), x_al, y_III ** 2))
    m.ordering_note = t_ga.ordering_note
    return m


def GaAsP(x_p: float) -> Material:
    """GaAs_1-x P_x."""
    return ternary("GaP", "GaAs", x_p, f"GaAs{1 - x_p:.2f}P{x_p:.2f}")


def InGaAs(x_in: float = 0.532) -> Material:
    return ternary("InAs", "GaAs", x_in, f"In{x_in:.3f}Ga{1 - x_in:.3f}As")


def InAsP(x_as: float) -> Material:
    return ternary("InAs", "InP", x_as, f"InAs{x_as:.2f}P{1 - x_as:.2f}")


def AlGaAs(x_al: float) -> Material:
    return ternary("AlAs", "GaAs", x_al, f"Al{x_al:.2f}Ga{1 - x_al:.2f}As")


def AlInAs(x_al: float = 0.48) -> Material:
    return ternary("AlAs", "InAs", x_al, f"Al{x_al:.2f}In{1 - x_al:.2f}As")


def AlGaInAs(x_al: float, y_III: float = 0.47) -> Material:
    """(Al_x Ga_1-x)_y In_1-y As; y = 0.47-0.48 lattice-matched to InP."""
    t_al = AlInAs(y_III)
    t_ga = InGaAs(1 - y_III)
    return quaternary_from_ternaries([t_al, t_ga], [x_al, 1 - x_al],
                                     f"(Al{x_al:.2f}Ga{1 - x_al:.2f}){y_III:.2f}In{1 - y_III:.2f}As")


# ------------------------------------------------------------------- band gaps etc.

def _varshni(E0, alpha, beta, T):
    T = np.asarray(T, float)
    correction = np.divide(alpha * T**2, T + beta,
                           out=np.zeros_like(T, dtype=float),
                           where=(T + beta) != 0)
    return E0 - correction


def bandgap(m: Material, T: float = 300.0, valley: str = "G") -> float:
    """Bandgap [eV] of valley 'G' | 'X' | 'L' at temperature T, or 'min' for
    the lowest of the three. V01 procedure: the binaries' temperature-
    dependent gaps are mixed linearly in composition, then the ternary bowing
    terms (temperature-independent) are subtracted."""
    if valley == "min":
        return min(bandgap(m, T, v) for v in ("G", "X", "L"))
    E = sum(f * _binary_gap(b, T, valley) for b, f in m.composition.items())
    for pair, xa, wt in m.bow_terms:
        E -= wt * xa * (1.0 - xa) * _bow(f"Eg_{valley}", pair, xa)
    return float(E)


def is_direct(m: Material, T: float = 300.0) -> bool:
    return bandgap(m, T, "G") <= min(bandgap(m, T, "X"), bandgap(m, T, "L"))


def lattice_constant(m: Material, T: float = 300.0) -> float:
    """[A]"""
    return float(m.p["a300"] + m.p["da_dT"] * (T - 300.0))


def mismatch(layer: Material, substrate: Material, T: float = 300.0) -> float:
    """In-plane misfit strain eps_par = (a_sub - a_layer)/a_layer (negative =
    compressive layer)."""
    a_l, a_s = lattice_constant(layer, T), lattice_constant(substrate, T)
    return (a_s - a_l) / a_l


def hh_mass_z(m: Material) -> float:
    """Heavy-hole mass along [001] from Luttinger: 1/(g1 - 2 g2)."""
    return 1.0 / (m.p["g1"] - 2.0 * m.p["g2"])


def lh_mass_z(m: Material) -> float:
    return 1.0 / (m.p["g1"] + 2.0 * m.p["g2"])


def hh_mass_inplane(m: Material) -> float:
    """In-plane heavy-hole mass (axial approx.): 1/(g1 + g2)."""
    return 1.0 / (m.p["g1"] + m.p["g2"])


# ------------------------------------------------------------------------ strain

@dataclass
class StrainedEdges:
    eps_par: float
    eps_perp: float
    tr_eps: float
    E_c: float       # absolute CB edge [eV] (InSb VB = 0 scale)
    E_hh: float      # absolute HH edge
    E_lh: float      # absolute LH edge
    E_v0: float      # unstrained VB top (absolute)
    E_c0: float      # unstrained CB edge (absolute)
    Eg_hh: float     # strained gap to HH
    Eg_lh: float
    dE_c: float
    dE_hh: float
    dE_lh: float


def band_edges(m: Material, T: float = 300.0) -> tuple[float, float]:
    """Unstrained absolute (E_v0, E_c0) on the V01 VBO scale."""
    E_v0 = m.p["vbo"]
    return E_v0, E_v0 + bandgap(m, T, "G")


def strain_shifts(layer: Material, substrate: Material, T: float = 300.0,
                  eps_par: float | None = None) -> StrainedEdges:
    """Pseudomorphic biaxial strain of `layer` on `substrate` (or a given
    eps_par), and the resulting absolute CB / HH / LH edges. Conventions in
    the module docstring."""
    if eps_par is None:
        eps_par = mismatch(layer, substrate, T)
    C11, C12 = layer.p["C11"], layer.p["C12"]
    eps_perp = -2.0 * (C12 / C11) * eps_par
    tr = 2.0 * eps_par + eps_perp
    a_c = layer.p["a_c"]
    a_v_vdw = -layer.p["a_v"]            # V01 table sign -> Van de Walle sign
    b = layer.p["b"]
    Q = b * (eps_perp - eps_par)         # < 0 under compression
    E_v0, E_c0 = band_edges(layer, T)
    dE_c = a_c * tr
    dE_av = a_v_vdw * tr
    E_hh = E_v0 + dE_av - Q
    # [V] Van de Walle, PRB 39, 1871 (1989), Eq. (8): exact k=0 LH-SO block.
    # The linear dE_av + Q form fails for strongly strained InP/InAs dots.
    delta = layer.p["d_so"]
    E_lh = E_v0 + dE_av - delta / 3.0 - delta / 6.0 + Q / 2.0 + 0.5 * np.sqrt(
        delta ** 2 + 2.0 * delta * Q + 9.0 * Q ** 2)
    E_c = E_c0 + dE_c
    return StrainedEdges(eps_par, eps_perp, tr, E_c, E_hh, E_lh, E_v0, E_c0,
                         E_c - E_hh, E_c - E_lh, dE_c, E_hh - E_v0, E_lh - E_v0)


def critical_thickness_nm(layer: Material, substrate: Material, T: float = 300.0,
                          nu: float | None = None, b_burgers_A: float | None = None) -> float:
    """Matthews-Blakeslee equilibrium critical thickness [nm] for a single
    strained layer with 60-degree misfit dislocations, solved iteratively:
        h_c = b (1 - nu cos^2 alpha) / (8 pi f (1+nu) cos lambda) [ln(h_c/b) + 1]
    with cos alpha = cos lambda = 1/2 (60-degree dislocations), f = |eps_par|,
    b = a/sqrt(2) the Burgers vector, nu = C12/(C11+C12). Matthews &
    Blakeslee, J. Cryst. Growth 27, 118 (1974), eq. (3) family; cf. People &
    Bean, APL 47, 322 (1985) for the comparison of conventions (single-layer
    vs capped-layer forms differ by a factor of 2). Returns inf for f = 0.
    [V] The single-layer form gives InP/GaAs about 1.3 nm and InAs/GaAs
    about 0.4 nm; the Stranski-Krastanov transition is a separate kinetic
    criterion, not an MB equilibrium thickness.
    [DR] -- equilibrium value; MOVPE/MBE layers can exceed it metastably."""
    f = abs(mismatch(layer, substrate, T))
    if f == 0:
        return float("inf")
    a = lattice_constant(layer, T)
    b = b_burgers_A if b_burgers_A is not None else a / np.sqrt(2.0)
    nu = nu if nu is not None else layer.p["C12"] / (layer.p["C11"] + layer.p["C12"])
    h = 10.0 * b
    for _ in range(200):
        h_new = (b * (1.0 - nu / 4.0) / (8.0 * np.pi * f * (1.0 + nu) * 0.5)) * (np.log(h / b) + 1.0)
        if h_new <= b:
            h_new = b * 1.0001
        if abs(h_new - h) < 1e-6 * h:
            h = h_new
            break
        h = 0.5 * (h + h_new)
    return float(h / 10.0)   # A -> nm


def strain_balance_thickness(layer_a: Material, t_a_nm: float, layer_b: Material,
                             substrate: Material, T: float = 300.0) -> float:
    """Thickness of layer_b [nm] that zero-balances the average in-plane
    stress of (layer_a, t_a) + (layer_b, t_b) on the substrate, using the
    thickness-weighted stress criterion  sum_i t_i A_i eps_i = 0 with
    A_i = C11 + C12 - 2 C12^2/C11 (Ekins-Daukes, Kawaguchi, Zhang, Cryst.
    Growth Des. 2, 287 (2002)). Returns inf if the two layers strain the same
    way (no compensation possible)."""
    def A(m):
        return m.p["C11"] + m.p["C12"] - 2.0 * m.p["C12"]**2 / m.p["C11"]
    ea, eb = mismatch(layer_a, substrate, T), mismatch(layer_b, substrate, T)
    if ea * eb >= 0:
        return float("inf")
    return float(-t_a_nm * A(layer_a) * ea / (A(layer_b) * eb))


# ---------------------------------------------------------- band offsets (heterojunction)

@dataclass
class Offsets:
    dE_c: float     # E_c(barrier) - E_c(well)  [eV]  (> 0: electrons confined in well)
    dE_v_hh: float  # E_hh(well) - E_hh(barrier)      (> 0: heavy holes confined in well)
    dE_v_lh: float
    type: str       # "I" if both > 0 else "II"
    well: StrainedEdges
    barrier: StrainedEdges


def offsets(well: Material, barrier: Material, substrate: Material,
            T: float = 300.0) -> Offsets:
    """Strained type-I/II band offsets of `well` inside `barrier`, both
    pseudomorphic on `substrate` (the barrier is usually lattice-matched, in
    which case its strain terms vanish)."""
    w = strain_shifts(well, substrate, T)
    b = strain_shifts(barrier, substrate, T)
    dEc = b.E_c - w.E_c
    dEv_hh = w.E_hh - b.E_hh
    dEv_lh = w.E_lh - b.E_lh
    typ = "I" if (dEc > 0 and dEv_hh > 0) else "II"
    return Offsets(dEc, dEv_hh, dEv_lh, typ, w, b)


# ---------------------------------------------------- extra per-material properties

# Thermal conductivity k300 [W/m/K] and exponent alpha for k = k300 (300/T)^alpha;
# refractive index samples n(lambda_nm) at the wavelengths this program uses
# (668 nm for the InP-dot red platform, 1310/1550 nm for the InAs/InP variant);
# thermo-optic dn/dT [1/K]. Sources per entry. Alloy k values are [E]: alloy
# scattering suppresses k by 5-10x relative to the binaries (Adachi, JAP 54,
# 1844 (1983), "Lattice thermal resistivity of III-V compound alloys").
MATERIAL_EXTRA = {
    "GaAs": {"k300": 55.0, "alpha_k": 1.25, "tag_k": "DR",
             "src_k": "Adachi 2007 / Ioffe NSM: 55 W/m/K at 300 K",
             "n": {668: 3.83, 930: 3.59, 1310: 3.41, 1550: 3.37}, "tag_n": "DR",
             "src_n": "668: repo dispersion table (Aspnes JAP 60, 754 (1986)); 930-1550: "
                      "Skauli et al., JAP 94, 6447 (2003) Sellmeier",
             "dn_dT": 2.67e-4, "absorbing_below_nm": 873},
    "AlAs": {"k300": 91.0, "alpha_k": 1.2, "tag_k": "DR",
             "src_k": "Adachi 2007: 0.91 W/cm/K",
             "n": {668: 3.08, 930: 2.97, 1310: 2.91}, "tag_n": "DR",
             "src_n": "Fern & Onton JAP 42, 3499 (1971) Sellmeier (repo dispersion table)",
             "dn_dT": 1.43e-4},
    "InP": {"k300": 68.0, "alpha_k": 1.4, "tag_k": "DR",
            "src_k": "Adachi 2007 / Ioffe NSM: 0.68 W/cm/K",
            "n": {668: 3.50, 930: 3.29, 1310: 3.20, 1550: 3.17}, "tag_n": "E",
            "src_n": "[E] Pettit & Turner JAP 36, 2081 (1965); 930 nm is an edge extrapolation",
            "dn_dT": 2.0e-4},
    "InAs": {"k300": 27.0, "alpha_k": 1.2, "tag_k": "DR",
             "src_k": "Adachi 2007: 0.27 W/cm/K",
             "n": {1310: 3.55, 1550: 3.51}, "tag_n": "E",
             "src_n": "Ioffe NSM (bulk InAs, above-gap region at 1.3 um: absorbing)"},
    "GaP": {"k300": 110.0, "alpha_k": 1.4, "tag_k": "DR",
            "src_k": "Adachi 2007: 1.1 W/cm/K",
            "n": {668: 3.30, 930: 3.19}, "tag_n": "DR",
            "src_n": "Bond, J. Appl. Phys. 36, 1674 (1965) / Ioffe NSM: n = 3.31 at 650 nm"},
    "Ga0.51In0.49P": {"k300": 5.2, "alpha_k": 0.5, "tag_k": "E",
                      "src_k": "Adachi 2007 alloy model (~0.05 W/cm/K for x~0.5); "
                               "measured 4.9-5.5 W/m/K class",
                      "n": {668: 3.573, 650: 3.599}, "tag_n": "V",
                      "src_n": "[V] Schubert et al., JAP 77, 3416 (1995): n = 3.573 at 666.74 nm, "
                               "3.599 at 650.15 nm",
                      "dn_dT": 2.5e-4},
    "(Al0.50Ga0.50)0.51In0.49P": {"k300": 6.0, "alpha_k": 0.5, "tag_k": "E",
                                  "src_k": "alloy class, Adachi 2007",
                                  "n": {668: 3.22, 650: 3.24}, "tag_n": "E",
                                  "src_n": "[E] AlxGa1-xInP: Moser et al., APL 64, 235 (1994); "
                                           "Schubert et al., JAP 86, 2025 (1999)",
                                  "dn_dT": 2.0e-4},
    "Al0.52In0.48P": {"k300": 8.0, "alpha_k": 0.5, "tag_k": "E",
                      "src_k": "alloy class, Adachi 2007",
                       "n": {668: 3.05, 650: 3.07}, "tag_n": "E",
                       "src_n": "[E] AlInP: Moser et al., APL 64, 235 (1994); "
                                "Schubert et al., JAP 86, 2025 (1999)",
                      "dn_dT": 1.8e-4},
    "GaAs0.60P0.40": {"k300": 10.0, "alpha_k": 0.5, "tag_k": "E",
                      "src_k": "Adachi 2007 alloy model, GaAsP x~0.4",
                      "n": {668: 3.55}, "tag_n": "E",
                      "src_n": "GaAs-GaP interpolation near the direct edge (E_g ~1.92 eV "
                               "at x = 0.4): the value is within ~0.1 of the true n; MEASURE"},
    "Al0.45Ga0.55As": {"k300": 9.8, "alpha_k": 0.5, "tag_k": "V",
                       "src_k": "[V] Ioffe NSM AlGaAs thermal formula: 0.55-2.12x+2.48x^2 W/cm/K",
                       "n": {668: 3.495}, "tag_n": "DR",
                       "src_n": "Papatryfonos et al., AIP Adv. 11, 025327 (2021) x=0.452 "
                                "(repo dispersion table)",
                       "dn_dT": 2.1e-4},
    "In0.53Ga0.47As": {"k300": 5.0, "alpha_k": 0.5, "tag_k": "E",
                       "src_k": "Adachi 2007: ~0.05 W/cm/K",
                       "n": {1310: 3.60, 1550: 3.56}, "tag_n": "DR",
                       "src_n": "Ioffe NSM / Adachi optical constants (absorbing at 1.31 um)"},
    "InGaAsP-Q1.15": {"k300": 4.5, "alpha_k": 0.5, "tag_k": "E",
                      "src_k": "Adachi 2007 quaternary alloy class",
                      "n": {1310: 3.36, 1550: 3.32}, "tag_n": "DR",
                      "src_n": "Adachi, JAP 53, 5863 (1982) InGaAsP/InP refractive-index model",
                      "dn_dT": 2.5e-4},
    "Si": {"k300": 148.0, "alpha_k": 1.35, "tag_k": "DR", "src_k": "Glassbrenner & Slack 1964",
           "n": {1310: 3.50, 1550: 3.48}, "tag_n": "DR", "src_n": "Li 1980 tabulation"},
}


def extra(label: str) -> dict:
    """Extra properties (thermal, optical) for a named material label; raises
    KeyError with the available labels if unknown (never silently guesses)."""
    try:
        return MATERIAL_EXTRA[label]
    except KeyError as e:
        raise KeyError(f"no MATERIAL_EXTRA entry for {label!r}; known: "
                       f"{sorted(MATERIAL_EXTRA)}") from e


# [E] Composition-aware fallback for two alloy families whose MATERIAL_EXTRA
# table only tabulates a few fixed x (cw-grid-and-algainp-index item 2): a
# requested composition not directly in the table is built by LINEAR
# interpolation in x, at the requested wavelength, between the nearest
# tabulated compositions -- Adachi, "Optical Constants of Crystalline and
# Amorphous Semiconductors"; Kato et al., J. Appl. Phys. 75, 3335 (1994) for
# the AlGaInP dispersion trend; Aspnes et al., JAP 60, 754 (1986) for AlGaAs.
# Anchors are (x, MATERIAL_EXTRA label) pairs sorted by x; the family is
# lattice-matched (Al_x Ga_1-x)_0.51 In_0.49 P / Al_x Ga_1-x As so a single x
# fully determines the composition. Below the direct gap n(x) is monotone
# (checked in verify_materials): interpolating between real, monotonically
# ordered tabulated points preserves that.
_ALGAINP_RE = re.compile(r"^\(Al(?P<xal>\d*\.?\d+)Ga\d*\.?\d+\)(?P<y>\d*\.?\d+)In\d*\.?\d+P$")
_ALGAAS_RE = re.compile(r"^Al(?P<xal>\d*\.?\d+)Ga\d*\.?\d+As$")
_ALGAINP_ANCHORS = ((0.0, "Ga0.51In0.49P"), (0.5, "(Al0.50Ga0.50)0.51In0.49P"), (1.0, "Al0.52In0.48P"))
_ALGAAS_ANCHORS = ((0.0, "GaAs"), (0.45, "Al0.45Ga0.55As"), (1.0, "AlAs"))


def _lookup_n(label: str, lambda_nm: float, extrap_tol_nm: float = 30.0) -> float:
    """Core single-composition wavelength lookup (linear interpolation
    between tabulated points; extrapolation refused beyond extrap_tol_nm of
    a lone tabulated point)."""
    tab = extra(label)["n"]
    lams = np.array(sorted(tab))
    ns = np.array([tab[l] for l in lams])
    if lambda_nm < lams.min() - 1e-9 or lambda_nm > lams.max() + 1e-9:
        if len(lams) == 1 and abs(lambda_nm - lams[0]) < extrap_tol_nm:
            return float(ns[0])
        raise ValueError(f"{label}: n tabulated only at {list(lams)} nm; asked {lambda_nm}")
    return float(np.interp(lambda_nm, lams, ns))


def _interp_n_by_x(x: float, anchors: tuple, lambda_nm: float) -> float:
    """[E] n(x) by linear interpolation between the two tabulated anchors
    bracketing x (clamped at the ends). Anchor lookups use a widened
    single-point extrapolation tolerance (150 nm) so a family stays usable
    across the 600-800 nm red-emission window even where an anchor
    composition has only one measured wavelength (e.g. Al0.45Ga0.55As)."""
    if not 0.0 <= x <= 1.0:
        raise ValueError(f"composition x must be in [0,1], got {x}")
    xs = [a[0] for a in anchors]
    i = 0
    while i < len(xs) - 2 and x > xs[i + 1]:
        i += 1
    (x_lo, lab_lo), (x_hi, lab_hi) = anchors[i], anchors[i + 1]
    n_lo = _lookup_n(lab_lo, lambda_nm, extrap_tol_nm=150.0)
    n_hi = _lookup_n(lab_hi, lambda_nm, extrap_tol_nm=150.0)
    t = (x - x_lo) / (x_hi - x_lo)
    return float(n_lo + t * (n_hi - n_lo))


def refractive_index(label: str, lambda_nm: float) -> float:
    """n at the nearest tabulated wavelength (linear interpolation between
    tabulated points; extrapolation refused).

    `label` may also be a Material (e.g. the object AlGaInP(x) or AlGaAs(x)
    returns) -- its `.label` is used. For a composition not directly in
    MATERIAL_EXTRA, labels of the form "(Al<x>Ga<1-x>)0.51In0.49P" or
    "Al<x>Ga<1-x>As" fall back to [E] linear interpolation in x between the
    nearest tabulated compositions of that family (see _interp_n_by_x)."""
    if not isinstance(label, str):
        label = label.label
    if label not in MATERIAL_EXTRA:
        m = _ALGAINP_RE.match(label)
        if m:
            y = float(m.group("y"))
            if abs(y - 0.51) > 0.02:
                raise ValueError(f"{label}: composition interpolation only covers the "
                                 "y=0.51 lattice-matched (AlxGa1-x)0.51In0.49P family")
            return _interp_n_by_x(float(m.group("xal")), _ALGAINP_ANCHORS, lambda_nm)
        m = _ALGAAS_RE.match(label)
        if m:
            return _interp_n_by_x(float(m.group("xal")), _ALGAAS_ANCHORS, lambda_nm)
    return _lookup_n(label, lambda_nm)


def thermal_k(label: str, T: float = 300.0) -> float:
    e = extra(label)
    return float(e["k300"] * (300.0 / T) ** e["alpha_k"])


# --------------------------------------------------------------- feasibility check

@dataclass
class StackLayerSpec:
    material: Material
    thickness_nm: float
    role: str = ""          # "substrate" | "cladding" | "barrier" | "well" | "dot" | "cap"


def feasibility_report(layers: list[StackLayerSpec], substrate: Material,
                       T: float = 300.0, emission_eV: float | None = None) -> dict:
    """Growth/optics feasibility of an epitaxial stack:
      * misfit and Matthews-Blakeslee critical thickness per layer, with a
        flag when thickness > h_c (relaxation / dislocations expected);
      * direct/indirect gap of each layer;
      * whether each non-active layer is transparent at `emission_eV`
        (E_g > emission energy + 3 kT margin);
      * net strain-thickness product (sign and magnitude) for strain balance.
    Returns a dict of per-layer rows plus overall verdict strings. Pure
    reporting: nothing here modifies the physics chain."""
    rows, warnings = [], []
    net_st = 0.0
    for L in layers:
        m = L.material
        eps = mismatch(m, substrate, T)
        hc = critical_thickness_nm(m, substrate, T)
        Eg = bandgap(m, T, "min")
        direct = is_direct(m, T)
        row = {"material": m.label, "role": L.role, "t_nm": L.thickness_nm,
               "eps_par": eps, "h_c_nm": hc, "Eg_min_eV": Eg, "direct": direct}
        if L.role not in ("dot", "well", "active") and L.thickness_nm > hc:
            warnings.append(f"{m.label} ({L.role}): {L.thickness_nm:.0f} nm exceeds "
                            f"h_c = {hc:.1f} nm at |eps| = {abs(eps):.2%} -> relaxation expected")
        if emission_eV is not None and L.role in ("cladding", "barrier", "cap", "substrate"):
            margin = emission_eV + 3.0 * KB_EV * T
            row["transparent"] = Eg > margin
            if Eg <= margin:
                warnings.append(f"{m.label} ({L.role}): E_g,min = {Eg:.3f} eV absorbs at "
                                f"{emission_eV:.3f} eV emission")
        # [V] Ekins-Daukes, Kawaguchi & Zhang, Cryst. Growth Des. 2,
        # 287 (2002): use the same biaxial average-stress weight as
        # strain_balance_thickness(), rather than an unweighted strain sum.
        A = m.p["C11"] + m.p["C12"] - 2.0 * m.p["C12"] ** 2 / m.p["C11"]
        net_st += A * eps * L.thickness_nm
        rows.append(row)
    return {"rows": rows, "warnings": warnings, "net_strain_thickness_nm": net_st,
            "substrate": substrate.label, "T": T}
