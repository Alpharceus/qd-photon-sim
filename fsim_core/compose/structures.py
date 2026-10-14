"""Structure registry entries of the composable engine (design 2.3, physics brief section 3).

Each entry's `response(vals, em, T_K, NA) -> StructureResponse` maps resolved
parameter values + an EmitterModel to the photonic response consumed by
rates.combine:

    Gamma = Gamma0 (F_bg + F_P L),   beta = (beta0 F_bg + F_P L) / (F_bg + F_P L)

Conventions (devlib-common):
  * F_P is the ADDITIONAL rate ratio (total = 1 + F_P L at F_bg = 1); it is the
    cavity-MODE Purcell factor at zero detuning, ideal position and polarisation
    unless a position factor is stated. L(Delta) is applied by rates.combine.
  * A leaky / background channel is F_bg (in units of Gamma0), never F_P.
  * Cavity detuning is detuning_ueV = E_X - E_mode (positive: the dot is above
    the mode, as in fsim_core.cavity); E_mode_meV = em.E_X_meV - detuning.
  * eta_out already includes the NA-dependent far-field collection (T_NA or the
    structure's own NA model) AT THE `NA` ARGUMENT; `eta_out_fn(NA)` returns the
    same quantity for other NA (sweeps). The engine must not multiply T_NA again.
  * Units: rates 1/ns, energies meV (detuning / kappa / g inputs in ueV), lengths
    as named in the parameter.

Provenance tags: [V] verified against the cited paper, [DR] derived, [E]
estimate / class range, [A] assumption.
"""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Callable, Optional

import numpy as np

from .. import dbr, nitride_cavity, nitride_nanowire_photonics as _nwp, waveguide
from ..spectral import cavity_transmission
from .params import ParamSpec
from .registry import Structure, register
from .types import StructureResponse

# hbar = 6.582119569e-16 eV s (CODATA 2018) = lindblad.HBAR_MEV_NS [V].
HBAR_MEV_NS = 6.582119569e-4
D_OVER_W0 = 2.7          # default micropillar d / w0 [E]: HE11 Gaussian approximation (Gerard et al., PRL 81, 1110 (1998) class); a ParamSpec leaf "d_over_w0"

_REF_PURCELL = "[DR] F_P = (3/(4 pi^2)) (lambda/n)^3 Q/V (Purcell, Phys. Rev. 69, 681 (1946))"
_REF_GERARD = "[E] Gerard et al., Phys. Rev. Lett. 81, 1110 (1998): micropillar Q/V class (V = (pi/4) w0^2 L_eff here; the cylinder pi d^2/4 L_eff is an upper bound)"
_REF_TOMM = "[V] Tomm et al., Nat. Nanotechnol. 16, 399 (2021): g/2pi 4.16, kappa/2pi 23.3, gamma/2pi 0.30 GHz, F_P ~ 10"
_REF_SOMASCHI = "[V] Somaschi et al., Nat. Photonics 10, 340 (2016): F_P 7.6, beta 0.88, eta_out 0.70"
_REF_LECAMP = "[DR] Lecamp, Lalanne & Hugonin, Phys. Rev. Lett. 99, 023902 (2007): F_wg = (3/(4 pi)) (lambda/n)^2 (n_g/n)/A_eff"
_REF_BLEUSE = "[E] Bleuse et al., Phys. Rev. Lett. 106, 103601 (2011): beta_HE11 peaks near d/lambda 0.22-0.25"


# --------------------------------------------------------------------------
# Public helpers
# --------------------------------------------------------------------------
def T_NA(NA, lambda_nm, w0_um):
    """Fraction of a Gaussian far field inside the cone NA [DR]:
    T_NA = 1 - exp(-2 NA^2 / theta_d^2), theta_d = lambda/(pi w0) (paraxial;
    NA is used as the half-angle sine = angle). Monotone in NA, -> 1 at large NA.
    Scalars give a float; arrays an array."""
    theta_d = (np.asarray(lambda_nm, dtype=float) * 1e-3) / (math.pi * np.asarray(w0_um, dtype=float))
    out = 1.0 - np.exp(-2.0 * np.asarray(NA, dtype=float) ** 2 / theta_d ** 2)
    return float(out) if np.ndim(out) == 0 else out


def lorentz_L(det, gamma_h, kappa):
    """Overlap L(Delta) of the emitter line (FWHM gamma_h) with the peak-normalised
    cavity Lorentzian (FWHM kappa): L = kappa (gamma_h + kappa)/4 / (Delta^2 +
    ((gamma_h + kappa)/2)^2) [DR]. Delegates to spectral.cavity_transmission (the
    kernel resonant_source.purcell_rate and nitride_cavity use), so it is the repo
    expression bit-for-bit. L(0) = kappa/(gamma_h + kappa) -> 1 for a narrow line."""
    return cavity_transmission(det, gamma_h, kappa)


def purcell_from_g(g_meV, kappa_meV, gamma0_ns):
    """F_P = 4 g^2/(kappa hbar Gamma0) [DR] (g, kappa as energies in meV)."""
    return 4.0 * g_meV ** 2 / (kappa_meV * HBAR_MEV_NS * gamma0_ns)


def g_from_purcell(F_P, kappa_meV, gamma0_ns):
    """Inverse of purcell_from_g: g = sqrt(F_P kappa hbar Gamma0/4) [DR]."""
    return 0.5 * math.sqrt(max(F_P, 0.0) * kappa_meV * HBAR_MEV_NS * gamma0_ns)


def parse_eta_table(text: str) -> list:
    """'NA:eta, NA:eta, ...' -> sorted [(NA, eta)]."""
    pts = []
    for tok in str(text).replace(";", ",").split(","):
        tok = tok.strip()
        if not tok:
            continue
        try:
            a, b = tok.split(":")
            pts.append((float(a), float(b)))
        except ValueError:
            raise ValueError(f"eta table entry {tok!r} is not 'NA:eta'") from None
    if not pts:
        raise ValueError("eta table is empty; expected 'NA:eta[, NA:eta ...]'")
    for a, b in pts:
        if not (a > 0.0 and 0.0 <= b <= 1.0):
            raise ValueError(f"eta table point {a}:{b} needs NA > 0 and eta in [0, 1]")
    return sorted(pts)


def eta_table_fn(text: str) -> Callable[[float], float]:
    """eta(NA) from a table: piecewise-linear between points, flat above the last
    point, and eta0 (NA/NA0)^2 below the first (cone-area scaling) [A]."""
    pts = parse_eta_table(text)
    xs = [p[0] for p in pts]
    ys = [p[1] for p in pts]

    def fn(NA):
        NA = float(NA)
        if NA <= 0.0:
            return 0.0
        if NA < xs[0]:
            return ys[0] * (NA / xs[0]) ** 2
        return float(np.interp(NA, xs, ys))
    return fn


def bulk_planar_eta(NA, n, mirror=False, sil="none", n_s=2.0, dipole="in_plane"):
    """First-lens collected fraction of a dipole in a planar high-index
    host, small cone [DR]: eta = (3/8) x^2 T_F for an in-plane dipole (azimuth-averaged
    pattern (3/16 pi)(1 + cos^2 th) per sr), (3/16) x^4 T_F for an axial one (pattern
    (3/8 pi) sin^2 th per sr: no light on axis), x = NA_eff/n capped at 1, NA_eff =
    NA x (1 | n_s | n_s^2) for no SIL | hemispherical | Weierstrass, T_F = 4n/(1+n)^2;
    a bottom mirror doubles it (capped at 1)."""
    mult = {"none": 1.0, "hemi": float(n_s), "weierstrass": float(n_s) ** 2}[sil]
    x = min(float(NA) * mult / float(n), 1.0)
    T_F = 4.0 * n / (1.0 + n) ** 2
    eta = (0.1875 * x ** 4 if str(dipole) == "axial" else 0.375 * x * x) * T_F
    if mirror:
        eta *= 2.0
    return float(min(eta, 1.0))


# --------------------------------------------------------------------------
# Dipole orientation (devlib-r2-fix1b, audit finding 10)
# --------------------------------------------------------------------------
# The material `dipole` (in_plane | axial) is the orientation of the transition dipole relative to the growth axis.
# Every cavity / waveguide mode of the entries below is polarised in the plane (HE11 and micropillar fields are
# transverse, photonic-crystal slab modes and TE ridge modes lie in the plane), so the coupling is the cos^2 projection
# of the dipole on the mode field [DR]: 1 for an in-plane dipole, 0 for an exactly axial one. An exactly axial dipole is
# an idealisation (tilt, light-hole admixture), so a small residual projection is kept [A].
AXIAL_RESIDUAL = 0.05          # [A] cos^2 projection of an axial dipole on an in-plane mode field (tilt / LH admixture)


def dipole_projection(em):
    """cos^2 projection of the emitter's transition dipole on an in-plane mode field: 1.0 (in_plane) or AXIAL_RESIDUAL
    (axial) [DR/A]. Applied to the COMPUTED Purcell factors (micropillar, geometry-mode open cavity, Q/V phc_cavity,
    phc_waveguide, edge_ridge, nitride cavity); entered F_P / g numbers (generic_cavity, cbg, direct open cavity, entered
    g) describe the device's own dot and are left as entered."""
    return AXIAL_RESIDUAL if str(getattr(em, "dipole", "in_plane")) == "axial" else 1.0


# --------------------------------------------------------------------------
# Entry machinery
# --------------------------------------------------------------------------
@dataclass(frozen=True)
class _Struct(Structure):
    fn: Optional[Callable] = None

    def response(self, vals, em, T_K, NA):
        return self.fn(self, vals, em, T_K, NA)


# display symbols (KaTeX TeX; Scully-Zubairy / Metcalf-van der Straten notation), keyed by the param id; UI only
_TEX = {
    "F_P": r"F_{\mathrm P}",
    "kappa_ueV": r"\hbar\kappa",
    "Q": r"Q",
    "detuning_ueV": r"\hbar\Delta_{\mathrm c}",
    "beta0": r"\beta_0",
    "eta_out": r"\eta_{\mathrm{out}}",
    "F_bg": r"F_{\mathrm{bg}}",
    "d_um": r"d",
    "N_top": r"N_{\mathrm{top}}",
    "N_bot": r"N_{\mathrm{bot}}",
    "n_high": r"n_{\mathrm H}",
    "n_low": r"n_{\mathrm L}",
    "Q_measured": r"Q_{\mathrm{meas}}",
    "loss_rt": r"\mathcal{L}_{\mathrm{rt}}",
    "Q_sw_1um": r"Q_{\mathrm{sw}}(1\,\mu\mathrm m)",
    "Q_sw_alpha": r"\alpha_{\mathrm{sw}}",
    "mode_area": r"\mathrm{mode\ area}",
    "d_over_w0": r"d/w_0",
    "r_factor": r"\xi_{\mathrm r}",
    "pol_factor": r"\cos^2\theta_{\mathrm{pol}}",
    "mode": r"\mathrm{mode}",
    "g_ueV": r"\hbar g",
    "ROC_um": r"R_{\mathrm{ROC}}",
    "L_um": r"L",
    "L_pen_um": r"L_{\mathrm{pen}}",
    "T_top": r"T_{\mathrm{top}}",
    "T_bot": r"T_{\mathrm{bot}}",
    "loss": r"\mathcal{L}",
    "pol_collect": r"\eta_{\mathrm{pol,c}}",
    "pol_dipole": r"\cos^2\theta_{\mu}",
    "eta_table": r"\eta(\mathrm{NA})",
    "table_mode": r"\mathrm{table\ mode}",
    "leak": r"F_{\mathrm{leak}}",
    "V_lambda3": r"V/(\lambda/n)^3",
    "spatial_overlap": r"\xi",
    "n_g": r"n_{\mathrm g}",
    "A_eff_um2": r"A_{\mathrm{eff}}",
    "beta_direct": r"\beta_{\mathrm{wg}}",
    "pos_factor": r"\xi_{\mathrm r}",
    "eta_coupler": r"\eta_{\mathrm{coup}}",
    "eta_prop": r"\eta_{\mathrm{prop}}",
    "d_nm": r"d",
    "n_wire": r"n_{\mathrm{wire}}",
    "n_group": r"n_{\mathrm{g,wire}}",
    "taper_T": r"T_{\mathrm{taper}}",
    "mirror_R": r"R_{\mathrm m}",
    "mirror_phase_ok": r"\mathrm{antinode}",
    "taper_mfr_nm": r"w_{\mathrm{taper}}",
    "n_cavity": r"n_{\mathrm c}",
    "n_substrate": r"n_{\mathrm{sub}}",
    "z_dot_frac": r"z_{\mathrm{dot}}/d_{\mathrm c}",
    "orientation": r"\hat{\boldsymbol\mu}",
    "n": r"n",
    "mirror": r"\mathrm{mirror}",
    "sil": r"\mathrm{SIL}",
    "n_s": r"n_{\mathrm{SIL}}",
    "rate_factor": r"F_{\mathrm{rate}}",
    "ridge_width_nm": r"w_{\mathrm{ridge}}",
    "etch_depth_nm": r"h_{\mathrm{etch}}",
    "alpha_cm": r"\alpha_{\mathrm{loss}}",
    "R_back": r"R_{\mathrm{back}}",
    "lambda_nm": r"\lambda",
    "na_method": r"\mathrm{NA\ method}",
    "dot_position": r"x_{\mathrm{dot}}/L",
    "mode_volume_norm": r"V/(\lambda/n)^3",
    "T_track": r"T_{\mathrm{track}}",
    "dEdT_cav_meV_K": r"\mathrm dE_{\mathrm c}/\mathrm dT",
    "detuning_offset_meV": r"\delta E_{\mathrm{off}}",
    "purcell_enabled": r"\mathrm{Purcell\ on}",
    "E_X_track_eV": r"E_X(T_{\mathrm{track}})",
}


def _p(name, unit, default, lo=None, hi=None, tag="A", source="", help="", scale="lin",
       choices=None, group="main", visible_if=None, sweepable=True):
    return ParamSpec(name, unit, default, lo=lo, hi=hi, scale=scale, choices=choices, tag=tag,
                     source=source or "assumed value; replace with the device's own number [A]",
                     help=help, group=group, visible_if=visible_if, sweepable=sweepable, tex=_TEX[name])


class _V:
    """Resolved-value accessor with registry defaults."""

    def __init__(self, entry, vals):
        self.e, self.v = entry, vals

    def __call__(self, name):
        x = self.v.get(name, self.e.param(name).default)
        return x

    def f(self, name):
        x = float(self(name))
        if not math.isfinite(x):
            raise ValueError(f"{self.e.kind}.{name} must be finite, got {x!r}")
        return x

    def pos(self, name):
        x = self.f(name)
        if x <= 0.0:
            raise ValueError(f"{self.e.kind}.{name} must be > 0, got {x!r}")
        return x

    def frac(self, name):
        x = self.f(name)
        if not 0.0 <= x <= 1.0:
            raise ValueError(f"{self.e.kind}.{name} must be in [0, 1], got {x!r}")
        return x


def _path(name):
    return f"structure.{name}"


@dataclass(frozen=True)
class StructureResponseX(StructureResponse):
    """StructureResponse plus `extra`: geometry numbers the structure itself derived (open-cavity waist, DBR layer
    thicknesses) so the 3D scene reads them from the response instead of re-deriving them (FX4, final review 10)."""
    extra: dict = field(default_factory=dict)


def _cavity_resp(entry, em, *, F_P, kappa_meV, detuning_meV, g_meV=None, F_bg=1.0, beta0=0.0,
                 eta_out, eta_out_fn=None, pol_factor=1.0, notes=(), inputs=None, tag_floor="", extra=None):
    E_mode = em.E_X_meV - detuning_meV
    if g_meV is None and kappa_meV:
        g_meV = g_from_purcell(F_P, kappa_meV, em.gamma0_ns)
    return StructureResponseX(F_P=float(F_P), kappa_meV=None if kappa_meV is None else float(kappa_meV),
                             E_mode_meV=float(E_mode), g_meV=None if g_meV is None else float(g_meV),
                             F_bg=float(F_bg), beta0=float(beta0), eta_out=float(eta_out),
                             eta_out_fn=eta_out_fn, pol_factor=float(pol_factor),
                             notes=tuple(notes), inputs=dict(inputs or {}), tag_floor=tag_floor,
                             extra=dict(extra or {}))


_COMMON_NOTE_EFFNA = ("eta_out includes the far-field/NA collection at the setup NA; "
                      "eta_out_fn(NA) gives it at other NA (do not multiply T_NA again)")


# --------------------------------------------------------------------------
# generic_cavity
# --------------------------------------------------------------------------
def _generic_cavity(self, vals, em, T_K, NA):
    v = _V(self, vals)
    F_P = v.f("F_P")
    if F_P < 0:
        raise ValueError("generic_cavity.F_P must be >= 0")
    det = v.f("detuning_ueV") * 1e-3
    E_mode = em.E_X_meV - det
    Q = v.f("Q")
    if Q > 0:
        kappa = E_mode / Q
        kin = "Q"
        notes = ["kappa = E_mode/Q from the entered Q [DR]"]
    else:
        kappa = v.pos("kappa_ueV") * 1e-3
        kin = "kappa_ueV"
        notes = ["kappa entered directly (Q = 0 means 'use kappa_ueV')"]
    notes.append("F_P is the user-entered cavity-mode Purcell factor at zero detuning [A]")
    inputs = {"F_P": _path("F_P"), "kappa_meV": _path(kin), "E_mode_meV": _path("detuning_ueV"),
              "F_bg": _path("F_bg"), "beta0": _path("beta0"), "eta_out": _path("eta_out")}
    return _cavity_resp(self, em, F_P=F_P, kappa_meV=kappa, detuning_meV=det, F_bg=v.pos("F_bg"),
                        beta0=v.frac("beta0"), eta_out=v.frac("eta_out"), notes=notes, inputs=inputs)


GENERIC_CAVITY = _Struct(
    kind="generic_cavity", label="Generic cavity (enter F_P, kappa)",
    summary=r"Any single-mode cavity described by its Purcell factor $F_{\mathrm P}$, linewidth $\hbar\kappa$, detuning $\hbar\Delta_{\mathrm c}$ and outcoupling $\eta_{\mathrm{out}}$; the Fock / legacy path.",
    params=(
        _p("F_P", "1", 5.0, 0.0, 1000.0, help=r"Cavity-mode Purcell factor $F_{\mathrm P}$ at zero detuning (ADDITIONAL rate ratio: total $= F_{\mathrm{bg}} + F_{\mathrm P} L$)."),
        _p("kappa_ueV", "ueV", 150.0, 1.0, 5000.0, scale="log", help=r"Cavity linewidth $\hbar\kappa$ (FWHM) in ueV; ignored when $Q > 0$."),
        _p("Q", "1", 0.0, 0.0, 1e6, help=r"Quality factor $Q$; if $> 0$, $\kappa = E_{\mathrm{mode}}/Q$ and $\hbar\kappa$ is ignored (0 = unused).", group="advanced"),
        _p("detuning_ueV", "ueV", 0.0, -5000.0, 5000.0, help=r"Dot minus mode energy $\hbar\Delta_{\mathrm c} = E_X - E_{\mathrm{mode}}$ (positive = dot above the mode)."),
        _p("beta0", "1", 0.0, 0.0, 1.0, help=r"Fraction $\beta_0$ of the background ($F_{\mathrm{bg}}$) emission that is already in the collected mode."),
        _p("eta_out", "1", 0.5, 0.0, 1.0, help=r"Outcoupling efficiency $\eta_{\mathrm{out}}$ from the mode to the first lens, including the far-field collection."),
        _p("F_bg", "1", 1.0, 0.0, 10.0, help=r"Background (non-cavity) decay rate $F_{\mathrm{bg}}$ in units of $\Gamma_0$ (1 = free; $< 1$ inhibited; $> 1$ leaky)."),
    ),
    refs=(_REF_PURCELL,), eta_hint=(0.05, 0.9), fn=_generic_cavity)


# --------------------------------------------------------------------------
# micropillar
# --------------------------------------------------------------------------
def _micropillar_calc(v, em):
    lam = em.lambda_nm
    n_h, n_l = v.pos("n_high"), v.pos("n_low")
    if n_h <= n_l:
        raise ValueError("micropillar needs n_high > n_low")
    n_c = n_h                                          # lambda-cavity spacer = high-index layer [A]
    Nt, Nb = int(v.f("N_top")), int(v.f("N_bot"))
    L_pen = dbr.penetration_depth(n_h, n_l, lam, n_c, "L")      # nm, medium n_c, low-index layer first
    L_eff_um = (lam / n_c + 2.0 * L_pen) * 1e-3
    d = v.pos("d_um")
    w0 = d / v.pos("d_over_w0")
    # V = int eps|E|^2 / (eps|E|^2)_max. Transverse: Gaussian HE11 profile, int exp(-2 r^2/w0^2) 2 pi r dr = pi w0^2/2.
    # Longitudinal: the standing wave halves the energy integral, int n^2|E|^2 dz/(n_c^2|E_max|^2) = L_eff/2
    # (reviewer's TMM of the 25/40-pair GaAs/AlAs stack: 489.1 nm = L_eff/2 = 978.3/2 nm). V = (pi/4) w0^2 L_eff [DR/E].
    # "cylinder" is the OLD pi d^2/4 x L_eff: it neither halves L_eff nor uses the HE11 area (J0 hard-wall area is only
    # 0.2695 x pi d^2/4, Gaussian 0.274x), so it is an UPPER BOUND on V, i.e. a LOWER BOUND on F_P.
    if v("mode_area") == "cylinder":
        V_um3 = math.pi * d * d / 4.0 * L_eff_um
    else:
        V_um3 = (math.pi / 4.0) * w0 * w0 * L_eff_um
    Rt = dbr.qw_peak_reflectivity(1.0, n_c, n_h, n_l, Nt)
    Rb = dbr.qw_peak_reflectivity(n_h, n_c, n_h, n_l, Nb)
    Tt, Tb = 1.0 - Rt, 1.0 - Rb
    loss = v.f("loss_rt")
    Qpl = 4.0 * math.pi * n_c * L_eff_um * 1e3 / (lam * (Tt + Tb + loss))   # = 2 pi (n L_eff/lambda) / ((Tt+Tb+loss)/2)
    Qsw = v.pos("Q_sw_1um") * d ** v.f("Q_sw_alpha")
    Qcalc = 1.0 / (1.0 / Qpl + 1.0 / Qsw)
    Qm = v.f("Q_measured")
    Q = Qm if Qm > 0 else Qcalc
    eta_mirror = Tt / (Tt + Tb + loss)
    eo = v.f("eta_out")
    if eo >= 0:
        eta_mirror = eo
    return dict(t_high_nm=lam / (4.0 * n_h), t_low_nm=lam / (4.0 * n_l), spacer_unit_nm=lam / n_c,
                lam=lam, n_c=n_c, L_pen_nm=L_pen, L_eff_um=L_eff_um, V_um3=V_um3, Tt=Tt, Tb=Tb,
                Q_planar=Qpl, Q_sw=Qsw, Q_calc=Qcalc, Q=Q, eta_mirror=eta_mirror, w0=w0,
                V_norm=V_um3 / (lam * 1e-3 / n_c) ** 3)


def _micropillar(self, vals, em, T_K, NA):
    v = _V(self, vals)
    c = _micropillar_calc(v, em)
    r_factor, pol = v.frac("r_factor"), v.frac("pol_factor")
    dp = dipole_projection(em)
    F_P = (3.0 / (4.0 * math.pi ** 2)) * (c["lam"] * 1e-3 / c["n_c"]) ** 3 * c["Q"] / c["V_um3"] * r_factor * pol * dp
    det = v.f("detuning_ueV") * 1e-3
    kappa = (em.E_X_meV - det) / c["Q"]
    lam, w0, em_ = c["lam"], c["w0"], c["eta_mirror"]
    fn = lambda NA_: em_ * T_NA(NA_, lam, w0)          # noqa: E731
    notes = [_COMMON_NOTE_EFFNA,
             f"V_eff = {'pi d^2/4 x L_eff (UPPER bound on V, lower bound on F_P)' if v('mode_area') == 'cylinder' else '(pi/4) w0^2 L_eff (Gaussian HE11 area pi w0^2/2 x standing-wave L_eff/2; w0 = d/' + format(v.f('d_over_w0'), 'g') + ' [E])'}, "
             f"L_eff = lambda/n + 2 L_pen, L_pen = lambda/(4 dn) n_h n_l/n_c^2 = {c['L_pen_nm']:.0f} nm "
             f"(dbr.penetration_depth, low-index layer first) [DR/E]; V/(lambda/n)^3 = {c['V_norm']:.2f}",
             "the ideal F_P (r_factor = pol_factor = 1, zero detuning) bounds a measured F_P from above; the measured/ideal "
             "ratio is r x cos^2(pol) x spectral overlap (Ding 2016: ideal ~13 vs measured 6.3, ratio ~0.47) [DR/E]",
             f"Q = {c['Q']:.0f} ({'measured, entered by the user' if v.f('Q_measured') > 0 else 'computed: 1/Q = 1/Q_planar + 1/Q_sw'}); "
             f"Q_planar = {c['Q_planar']:.0f} from DBR T_top = {c['Tt']:.2e}, T_bot = {c['Tb']:.2e} [DR], "
             f"Q_sw = {c['Q_sw']:.0f} (A d^alpha, [A])",
             "fractional DBR pair counts are truncated to integers (dbr uses int(n_pairs))",
             "position factor r_factor and polarisation factor pol_factor (cos^2 dipole-field) are both folded into F_P (brief 3); "
             "the response pol_factor is 1",
             f"material dipole {getattr(em, 'dipole', 'in_plane')}: cos^2 projection on the transverse HE11 field = {dp:.4g} folded into F_P [DR/A]"]
    inputs = {"F_P": _path("Q_measured") if v.f("Q_measured") > 0 else _path("d_um"),
              "kappa_meV": _path("detuning_ueV"), "eta_out": _path("N_top"), "pol_factor": _path("pol_factor"),
              "beta0": _path("beta0"), "F_bg": _path("F_bg")}
    return _cavity_resp(self, em, F_P=F_P, kappa_meV=kappa, detuning_meV=det, F_bg=v.pos("F_bg"),
                        beta0=v.frac("beta0"), eta_out=fn(NA), eta_out_fn=fn, pol_factor=1.0,
                        notes=notes, inputs=inputs, tag_floor="E",           # dbr.penetration_depth, n_c = n_high [E]
                        extra={"t_high_nm": c["t_high_nm"], "t_low_nm": c["t_low_nm"],
                               "spacer_unit_nm": c["spacer_unit_nm"], "w0_um": c["w0"], "L_eff_um": c["L_eff_um"]})


MICROPILLAR = _Struct(
    kind="micropillar", label="DBR micropillar",
    summary="Etched GaAs/AlAs DBR micropillar: Q from the mirrors and sidewall scattering (or "
            "measured), V from the pillar diameter, outcoupling from the mirror asymmetry.",
    params=(
        _p("d_um", "um", 2.5, 0.5, 10.0, tag="V", source="Ding et al., PRL 116, 020401 (2016): d = 2.5 um", help=r"Pillar diameter $d$."),
        _p("N_top", "1", 25, 1, 60, help=r"Top DBR pairs $N_{\mathrm{top}}$ (fractional counts are truncated)."),
        _p("N_bot", "1", 40, 1, 60, help=r"Bottom DBR pairs $N_{\mathrm{bot}}$ (fractional counts are truncated)."),
        _p("n_high", "1", 3.5, 2.0, 4.0, tag="E", source="GaAs near 930 nm [E]", help=r"High-index DBR layer and spacer index $n_{\mathrm H}$ (GaAs)."),
        _p("n_low", "1", 2.95, 1.2, 3.5, tag="E", source="AlAs near 930 nm [E]", help=r"Low-index DBR layer index $n_{\mathrm L}$ (AlAs)."),
        _p("Q_measured", "1", 0.0, 0.0, 1e6, tag="A", source="user-entered measured Q", help=r"Measured $Q_{\mathrm{meas}}$; if $> 0$ it overrides the computed $Q$ (0 = compute)."),
        _p("loss_rt", "1", 1e-4, 0.0, 0.1, help=r"Round-trip absorption and scattering loss $\mathcal{L}_{\mathrm{rt}}$ of the planar cavity.", group="advanced"),
        _p("Q_sw_1um", "1", 200.0, 1.0, 1e9, scale="log", help=r"Sidewall-scattering quality factor of a 1 um pillar, $Q_{\mathrm{sw}} = A\,d^{\alpha_{\mathrm{sw}}}$.", group="advanced"),
        _p("Q_sw_alpha", "1", 4.0, 0.0, 8.0, help=r"Diameter exponent $\alpha_{\mathrm{sw}}$ of the sidewall-scattering $Q_{\mathrm{sw}}$.", group="advanced"),
        _p("mode_area", "enum", "gaussian", choices=("gaussian", "cylinder"),
           help=r"Mode volume: Gaussian HE11, $V = (\pi/4)\,w_0^2 L_{\mathrm{eff}}$ (default), or 'cylinder', $V = \pi d^2 L_{\mathrm{eff}}/4$ (an upper bound on $V$, i.e. a lower bound on $F_{\mathrm P}$).", group="advanced"),
        _p("d_over_w0", "1", D_OVER_W0, 1.5, 4.0, tag="E", source="HE11 Gaussian approximation w0 ~ d/2.7 [E] (Gerard et al., PRL 81, 1110 (1998) class)",
           help=r"Pillar diameter over Gaussian waist radius, $w_0 = d/(d/w_0)$.", group="advanced"),
        _p("detuning_ueV", "ueV", 0.0, -2000.0, 2000.0, help=r"Dot minus mode energy $\hbar\Delta_{\mathrm c} = E_X - E_{\mathrm{mode}}$ (positive = dot above the mode)."),
        _p("r_factor", "1", 1.0, 0.0, 1.0, help=r"Position factor $\xi_{\mathrm r} = |E(\mathbf r)|^2/|E_{\max}|^2$ of the dot (folded into $F_{\mathrm P}$)."),
        _p("pol_factor", "1", 1.0, 0.0, 1.0, help=r"Polarisation factor $\cos^2\theta_{\mathrm{pol}}$ between the dipole and the cavity field (folded into $F_{\mathrm P}$)."),
        _p("eta_out", "1", -1.0, -1.0, 1.0, help=r"Mirror outcoupling override $\eta_{\mathrm{out}}$; negative = computed $T_{\mathrm{top}}/(T_{\mathrm{top}} + T_{\mathrm{bot}} + \mathcal{L}_{\mathrm{rt}})$."),
        _p("beta0", "1", 0.0, 0.0, 1.0, help=r"Fraction $\beta_0$ of background emission already in the mode."),
        _p("F_bg", "1", 1.0, 0.0, 10.0, tag="A", source="leaky rate = Gamma0 [A], Somaschi 2016 (beta 0.88 at F_P 7.6) [V]",
           help=r"Background decay rate $F_{\mathrm{bg}}$ in units of $\Gamma_0$."),
    ),
    refs=(_REF_PURCELL, _REF_GERARD, _REF_SOMASCHI), eta_hint=(0.1, 0.8), viz="micropillar", fn=_micropillar)


# --------------------------------------------------------------------------
# open_cavity
# --------------------------------------------------------------------------
def _open_cavity(self, vals, em, T_K, NA):
    v = _V(self, vals)
    lam_um = em.lambda_nm * 1e-3
    det = v.f("detuning_ueV") * 1e-3
    E_mode = em.E_X_meV - det
    if v("mode") == "direct":
        g = v.pos("g_ueV") * 1e-3
        kappa = v.pos("kappa_ueV") * 1e-3
        F_P = purcell_from_g(g, kappa, em.gamma0_ns)
        eta = v.frac("eta_out")
        notes = ["direct (g, kappa) entry, F_P = 4 g^2/(kappa hbar Gamma0) [DR]; the preferred route (Tomm 2021)"]
        inputs = {"F_P": _path("g_ueV"), "kappa_meV": _path("kappa_ueV"), "g_meV": _path("g_ueV"),
                  "eta_out": _path("eta_out"), "E_mode_meV": _path("detuning_ueV")}
    else:
        R, L = v.pos("ROC_um"), v.pos("L_um")
        if L >= R:
            raise ValueError("open_cavity geometry mode needs L_um < ROC_um (stable cavity)")
        Tt, Tb, loss = v.f("T_top"), v.f("T_bot"), v.f("loss")
        w0_2 = (lam_um / math.pi) * math.sqrt(L * (R - L))
        Leff = L + v.f("L_pen_um")
        V = (math.pi / 4.0) * w0_2 * Leff
        Q = 4.0 * math.pi * Leff / (lam_um * (Tt + Tb + loss))
        kappa = E_mode / Q
        F_P = (3.0 / (4.0 * math.pi ** 2)) * lam_um ** 3 * Q / V
        g = g_from_purcell(F_P, kappa, em.gamma0_ns)
        eta = Tt / (Tt + Tb + loss)
        notes = [f"geometry mode: w0^2 = (lambda/pi) sqrt(L (R - L)), V = (pi/4) w0^2 L_eff, "
                 f"Q = 4 pi L_eff/(lambda (T_top + T_bot + loss)) = {Q:.0f} [DR]; air gap n = 1",
                 "eta_out = T_top/(T_top + T_bot + loss) [DR]; fibre / mode-match is not included (enter via eta_out in direct mode)"]
        inputs = {"F_P": _path("ROC_um"), "kappa_meV": _path("L_um"), "eta_out": _path("T_top"),
                  "E_mode_meV": _path("detuning_ueV")}
    R_g, L_g = v.f("ROC_um"), v.f("L_um")
    lam_g = em.lambda_nm * 1e-3
    waist = math.sqrt((lam_g / math.pi) * math.sqrt(L_g * (R_g - L_g))) if 0.0 < L_g < R_g else None
    pol_d = v.frac("pol_dipole")
    dp = dipole_projection(em)
    if v("mode") != "direct" and dp < 1.0:
        pol_d *= dp
        notes.append(f"material dipole {em.dipole}: cos^2 projection {dp:.4g} on the transverse cavity field folded into pol_dipole "
                     "(geometry mode only; a directly entered g is the device's own coupling) [DR/A]")
    if pol_d < 1.0:
        F_P *= pol_d
        g *= math.sqrt(pol_d)
        notes.append(f"pol_dipole = {pol_d:.4g}: dipole cos^2 folded into g^2 (g x sqrt) and F_P [DR]")
    notes.append("pol_collect is a brightness-only collection factor (response pol_factor); the dipole alignment pol_dipole "
                 "is in g^2 [FX4]")
    notes.append("g is returned so the regime test can compare 4g with kappa + Gamma; large g -> strong coupling (JC)")
    inputs.update({"F_bg": _path("F_bg"), "beta0": _path("beta0"), "pol_factor": _path("pol_collect"),
                   "F_P": inputs["F_P"] + "+" + _path("pol_dipole"), "g_meV": inputs.get("g_meV", inputs["F_P"]) + "+" + _path("pol_dipole")})
    return _cavity_resp(self, em, F_P=F_P, kappa_meV=kappa, detuning_meV=det, g_meV=g, F_bg=v.pos("F_bg"),
                        beta0=v.frac("beta0"), eta_out=eta, pol_factor=v.frac("pol_collect"),
                        notes=notes, inputs=inputs,
                        extra={"waist_um": waist, "stable": waist is not None})


OPEN_CAVITY = _Struct(
    kind="open_cavity", label="Open / fibre micro-cavity",
    summary="Tunable open cavity: enter the coupling g and linewidth kappa directly (preferred), or "
            "radius of curvature, length and mirror transmissions.",
    params=(
        _p("mode", "enum", "direct", choices=("direct", "geometry"), help=r"Enter $(g, \kappa)$ directly, or build them from the geometry."),
        _p("g_ueV", "ueV", 17.2, 0.1, 500.0, tag="V", source="Tomm 2021: g/2pi = 4.16 GHz = 17.2 ueV [V]",
           help=r"Dot-cavity coupling $\hbar g$ (not the Rabi splitting $2\hbar g$).", visible_if={"mode": "direct"}),
        _p("kappa_ueV", "ueV", 96.4, 1.0, 5000.0, tag="V", source="Tomm 2021: kappa/2pi = 23.3 GHz = 96.4 ueV [V]",
           help=r"Cavity linewidth $\hbar\kappa$ (FWHM).", visible_if={"mode": "direct"}),
        _p("ROC_um", "um", 20.0, 1.0, 1000.0, help=r"Radius of curvature $R_{\mathrm{ROC}}$ of the concave mirror.", visible_if={"mode": "geometry"}),
        _p("L_um", "um", 2.0, 0.1, 100.0, help=r"Mirror separation $L$ (air gap).", visible_if={"mode": "geometry"}),
        _p("L_pen_um", "um", 0.0, 0.0, 5.0, help=r"Total mirror penetration depth $L_{\mathrm{pen}}$ added to the length (0 = ignore).",
           visible_if={"mode": "geometry"}, group="advanced"),
        _p("T_top", "1", 1e-3, 1e-6, 0.5, scale="log", help=r"Output-mirror power transmission $T_{\mathrm{top}}$.", visible_if={"mode": "geometry"}),
        _p("T_bot", "1", 5e-5, 1e-6, 0.5, scale="log", help=r"Back-mirror power transmission $T_{\mathrm{bot}}$.", visible_if={"mode": "geometry"}),
        _p("loss", "1", 5e-5, 0.0, 0.1, help=r"Round-trip scattering and absorption loss $\mathcal{L}$.", visible_if={"mode": "geometry"}),
        _p("eta_out", "1", 0.5, 0.0, 1.0, help=r"Mode-to-detector outcoupling $\eta_{\mathrm{out}}$ including the fibre match (direct mode).", visible_if={"mode": "direct"}),
        _p("detuning_ueV", "ueV", 0.0, -2000.0, 2000.0, help=r"Dot minus mode energy $\hbar\Delta_{\mathrm c} = E_X - E_{\mathrm{mode}}$ (positive = dot above the mode)."),
        _p("beta0", "1", 0.0, 0.0, 1.0, help=r"Fraction $\beta_0$ of background emission already in the mode."),
        _p("pol_collect", "1", 1.0, 0.0, 1.0,
           help=r"Polarisation collection factor $\eta_{\mathrm{pol,c}}$ (brightness only): fraction of the cavity output reaching the detector in the detected polarisation (mode splitting, analyser); it does not change the dot-cavity coupling."),
        _p("pol_dipole", "1", 1.0, 0.0, 1.0,
           help=r"Dipole alignment $\cos^2\theta_{\mu}$ between the dipole and the cavity field: scales $g^2$ (and $F_{\mathrm P}$) of the dot-cavity coupling, $g \to g\sqrt{\cos^2\theta_{\mu}}$. 1 = aligned.", group="advanced"),
        _p("F_bg", "1", 1.0, 0.0, 10.0, help=r"Background decay rate $F_{\mathrm{bg}}$ in units of $\Gamma_0$."),
    ),
    refs=(_REF_TOMM, _REF_PURCELL), eta_hint=(0.2, 0.7), viz="open_cavity", fn=_open_cavity)


# --------------------------------------------------------------------------
# cbg
# --------------------------------------------------------------------------
def _cbg(self, vals, em, T_K, NA):
    v = _V(self, vals)
    F_P = v.f("F_P")
    if F_P < 0:
        raise ValueError("cbg.F_P must be >= 0")
    r_factor = v.frac("r_factor")
    F_P = F_P * r_factor                        # position factor folded into F_P like the micropillar [DR]
    det = v.f("detuning_ueV") * 1e-3
    Q = v.pos("Q")
    kappa = (em.E_X_meV - det) / Q
    leak = v.f("leak")
    if leak < 0:
        raise ValueError("cbg.leak must be >= 0")
    tab = eta_table_fn(v("eta_table"))
    mode = v("table_mode")
    if mode == "cavity_only":
        # legacy (devlib-fx2a) behaviour: only the cavity-mode channel is collected at table/beta_res (capped at 1),
        # so a low-Purcell CBG collapses to beta = F/(1 + F + leak) [DR]
        beta_res = F_P / (1.0 + leak + F_P) if F_P > 0 else 0.0
        fn = (lambda NA_: min(tab(NA_) / beta_res, 1.0)) if beta_res > 0 else (lambda NA_: 0.0)      # noqa: E731
        beta0 = 0.0
        notes = [_COMMON_NOTE_EFFNA,
                 f"table_mode = cavity_only: eta_table = TOTAL first-lens efficiency (beta x eta_out) only for the cavity "
                 f"channel; eta_out = table/beta_res with beta_res = F_P/(1 + F_P + leak) = {beta_res:.3f} (capped at 1); "
                 f"beta = F_P/(1 + F_P + leak) via F_bg = 1 + leak = {1.0 + leak:.3g} [DR]"]
    else:
        # FX3a: the table is the TOTAL first-lens efficiency of the DEVICE at its design F_P (the FDTD / FEM number).
        # A CBG's broadband extraction is far-field redistribution by the grating, which acts on ALL emitted light, not
        # only on the cavity-mode photons (Davanco, Appl. Phys. Lett. 99, 041102 (2011); Liu, Nat. Nanotechnol. 14, 586
        # (2019): non-Purcell CBG sources of efficiency 0.6-0.8 at F_P ~ 1-3). So every radiative photon enters the
        # collected channel (beta0 = 1) and is collected at eta_total(NA); F_P and the leak only set the RATE [E].
        fn = lambda NA_: min(max(tab(NA_), 0.0), 1.0)                                                  # noqa: E731
        beta0 = 1.0
        notes = [_COMMON_NOTE_EFFNA,
                 "table_mode = total: eta_table is the TOTAL first-lens efficiency eta_total(NA) of ALL radiative emission "
                 "(beta0 = 1, eta_out = eta_total(NA)); F_P and leak set only the decay RATE, Gamma = Gamma0 (1 + leak + "
                 "F_P L) [E] (grating far-field redistribution, Davanco 2011 / Liu 2019 class). Not valid if the entered "
                 "table is a cavity-photons-only number: use table_mode = cavity_only",
                 "consequence: brightness no longer collapses to beta = F/(1 + F) for low-Purcell CBGs; detuning changes "
                 "the rate and lifetime but not the collected fraction [E]"]
    notes += ["no mode solver: F_P, Q and the far-field table are entered numbers [A] (FDTD / COMSOL values)",
              f"position factor r_factor = {r_factor:.3g} (|E(r)|^2/|E_max|^2 of the dot) is folded into F_P like the micropillar: "
              f"F_P used = {F_P:.4g} [DR]",
              "eta_out(NA) from the table: piecewise linear, flat above the last point, eta0 (NA/NA0)^2 below the first [A]"]
    inputs = {"F_P": _path("F_P") + "+r_factor", "kappa_meV": _path("Q"), "E_mode_meV": _path("detuning_ueV"),
              "F_bg": _path("leak"), "eta_out": _path("eta_table") + "+table_mode"}
    return _cavity_resp(self, em, F_P=F_P, kappa_meV=kappa, detuning_meV=det, F_bg=1.0 + leak,
                        beta0=beta0, eta_out=fn(NA), eta_out_fn=fn, notes=notes, inputs=inputs)


CBG = _Struct(
    kind="cbg", label="Circular Bragg grating (bullseye)",
    summary="Broadband bullseye resonator: Purcell factor, low Q, a far-field efficiency table "
            "versus NA and a leaky-channel fraction, all entered numbers.",
    params=(
        _p("F_P", "1", 10.0, 0.0, 100.0, tag="A", source="class value; Wang 2019 PRL F_P 11.3, Kolatschek F_P ~ 4, Gschwandtner F_P ~ 18 [V-abs]",
           help=r"Cavity-mode Purcell factor $F_{\mathrm P}$ at zero detuning (additional rate ratio)."),
        _p("Q", "1", 150.0, 20.0, 1000.0, tag="E", source="CBG Q ~ 100-300 class [E]", help=r"Quality factor $Q$ (broad: $\hbar\kappa \sim$ 5-10 meV)."),
        _p("eta_table", "", "0.7:0.6", tag="A", source="single FDTD-class point (NA 0.7 -> 0.6) [A]; Wang 2019 PRL 0.62 at F 11.3 [V-abs]",
           help=r"TOTAL first-lens efficiency table $\eta(\mathrm{NA})$ 'NA:eta[, NA:eta ...]' of the device at its design $F_{\mathrm P}$ (the FDTD / FEM number, as quoted in the papers)."),
        _p("table_mode", "enum", "total", choices=("total", "cavity_only"), tag="E",
           source="total: grating far-field redistribution acts on all emission (Davanco 2011 / Liu 2019 class) [E]; cavity_only = legacy model",
           help=r"total: the table applies to ALL radiative emission ($\beta_0 = 1$; $F_{\mathrm P}$ only sets the rate). cavity_only: legacy, only the cavity channel is collected ($\eta_{\mathrm{out}} = \eta/\beta_{\mathrm{res}}$).", group="advanced"),
        _p("r_factor", "1", 1.0, 0.0, 1.0, help=r"Position factor $\xi_{\mathrm r} = |E(\mathbf r)|^2/|E_{\max}|^2$ of the dot (folded into $F_{\mathrm P}$ like the micropillar)."),
        _p("leak", "1", 0.0, 0.0, 10.0, help=r"Extra leaky decay $F_{\mathrm{leak}}$ in units of $\Gamma_0$ (cavity_only mode: $\beta = F_{\mathrm P}/(1 + F_{\mathrm P} + F_{\mathrm{leak}})$; total mode: rate only)."),
        _p("detuning_ueV", "ueV", 0.0, -20000.0, 20000.0, help=r"Dot minus mode energy $\hbar\Delta_{\mathrm c} = E_X - E_{\mathrm{mode}}$ (positive = dot above the mode)."),
    ),
    refs=(_REF_PURCELL,), eta_hint=(0.2, 0.7), viz="bullseye", fn=_cbg)


# --------------------------------------------------------------------------
# phc_cavity
# --------------------------------------------------------------------------
def _phc_cavity(self, vals, em, T_K, NA):
    v = _V(self, vals)
    Q, Vn = v.pos("Q"), v.pos("V_lambda3")
    det = v.f("detuning_ueV") * 1e-3
    kappa = (em.E_X_meV - det) / Q
    notes = [f"F_P = 3/(4 pi^2) Q/(V/(lambda/n)^3) x spatial_overlap [DR]; eta_out is an [E] vertical-emission estimate "
             "unless the far field is engineered",
             "usually strong coupling (g ~ 100 ueV vs kappa ~ 50-100 ueV): the regime test flags it and the JC solver is needed"]
    g_in = v.f("g_ueV")
    dp = dipole_projection(em)
    F_P = (3.0 / (4.0 * math.pi ** 2)) * Q / Vn * v.frac("spatial_overlap") * dp
    if dp < 1.0:
        notes.append(f"material dipole {em.dipole}: cos^2 projection {dp:.4g} on the in-plane slab mode folded into the Q/V F_P "
                     "(not into an entered g) [DR/A]")
    inputs = {"F_P": _path("Q") + "+V_lambda3+spatial_overlap", "kappa_meV": _path("Q"), "E_mode_meV": _path("detuning_ueV"),
              "eta_out": _path("eta_out"), "F_bg": _path("F_bg"), "beta0": _path("beta0")}
    if g_in > 0:
        g = g_in * 1e-3
        F_g = purcell_from_g(g, kappa, em.gamma0_ns)
        notes.append(f"g_ueV entered: F_P = 4 g^2/(kappa hbar Gamma0) = {F_g:.3g} replaces the Q/V value {F_P:.3g} [DR]")
        F_P = F_g
        inputs.update({"F_P": _path("g_ueV"), "g_meV": _path("g_ueV")})
    else:
        g = None
    return _cavity_resp(self, em, F_P=F_P, kappa_meV=kappa, detuning_meV=det, g_meV=g, F_bg=v.pos("F_bg"),
                        beta0=v.frac("beta0"), eta_out=v.frac("eta_out"), notes=notes, inputs=inputs)


PHC_CAVITY = _Struct(
    kind="phc_cavity", label="Photonic-crystal cavity (L3 / H1 / nanobeam)",
    summary="High-Q small-volume photonic-crystal cavity: Q and V in cubic wavelengths, optional g, "
            "entered vertical outcoupling.",
    params=(
        _p("Q", "1", 1e4, 100.0, 1e7, scale="log", tag="E", source="L3 / H1 class Q 1e3-1e5 [E]", help=r"Quality factor $Q = \omega_{\mathrm c}/\kappa$."),
        _p("V_lambda3", "1", 0.7, 0.3, 10.0, tag="E", source="L3 class V ~ 0.7 (lambda/n)^3 [E]", help=r"Mode volume $V$ in units of $(\lambda/n)^3$."),
        _p("g_ueV", "ueV", 0.0, 0.0, 500.0, help=r"Coupling $\hbar g$; if $> 0$, $F_{\mathrm P} = 4g^2/(\kappa\Gamma_0)$ replaces the $Q/V$ value (0 = unused)."),
        _p("spatial_overlap", "1", 1.0, 0.0, 1.0, help=r"Position and polarisation overlap factor $\xi = |E(\mathbf r)|^2/|E_{\max}|^2 \times \cos^2\theta$ (folded into $F_{\mathrm P}$)."),
        _p("eta_out", "1", 0.05, 0.0, 1.0, tag="E", source="few % vertical unless far-field engineered [E]", help=r"Mode-to-first-lens outcoupling $\eta_{\mathrm{out}}$."),
        _p("detuning_ueV", "ueV", 0.0, -2000.0, 2000.0, help=r"Dot minus mode energy $\hbar\Delta_{\mathrm c} = E_X - E_{\mathrm{mode}}$ (positive = dot above the mode)."),
        _p("beta0", "1", 0.0, 0.0, 1.0, help=r"Fraction $\beta_0$ of background emission already in the mode."),
        _p("F_bg", "1", 1.0, 0.0, 10.0, help=r"Background decay rate $F_{\mathrm{bg}}$ in units of $\Gamma_0$."),
    ),
    refs=(_REF_PURCELL,), eta_hint=(0.01, 0.2), viz="phc_cavity", fn=_phc_cavity)


# --------------------------------------------------------------------------
# phc_waveguide
# --------------------------------------------------------------------------
def _phc_waveguide(self, vals, em, T_K, NA):
    v = _V(self, vals)
    F_bg = v.pos("F_bg")
    bd = v.f("beta_direct")
    notes = ["broadband (kappa = None): L = 1; F_wg is the guided-mode rate in units of Gamma0 [DR]",
             "beta = F_wg/(F_wg + F_bg) with F_bg the leaky / background rate (< 1 = inhibited) [DR]",
             "reuses waveguide.beta_factor (Lecamp 2007 with the (n_g/n) energy-normalised form) rather than "
             "purcell_mode_width, which needs a field profile",
             "eta_out = eta_coupler x eta_prop; no NA dependence (guided output)"]
    if bd > 0:
        if not 0 < bd < 1:
            raise ValueError("phc_waveguide.beta_direct must be in (0, 1) when used")
        F = bd / (1.0 - bd) * F_bg
        notes.append(f"beta entered directly: F_wg = beta/(1 - beta) F_bg = {F:.3g}")
        inputs = {"F_P": _path("beta_direct")}
    else:
        F, _ = waveguide.beta_factor(v.pos("A_eff_um2"), em.lambda_nm, em.n_host, v.pos("n_g"), v.frac("pos_factor"))
        F *= v.frac("pol_factor") * dipole_projection(em)      # cos^2 dipole-mode overlap folded into F_wg (brief 3)
        if dipole_projection(em) < 1.0:
            notes.append(f"material dipole {em.dipole}: cos^2 projection {dipole_projection(em):.4g} on the in-plane TE mode folded "
                         "into F_wg [DR/A]")
        notes.append("pol_factor (cos^2 dipole-mode) is folded into F_wg; the response pol_factor is 1")
        inputs = {"F_P": _path("A_eff_um2")}
    inputs.update({"F_bg": _path("F_bg"), "eta_out": _path("eta_coupler")})
    return StructureResponse(F_P=float(F), kappa_meV=None, E_mode_meV=None, g_meV=None, F_bg=F_bg, beta0=0.0,
                             eta_out=v.frac("eta_coupler") * v.frac("eta_prop"), eta_out_fn=None,
                             pol_factor=1.0, notes=tuple(notes), inputs=inputs)


PHC_WAVEGUIDE = _Struct(
    kind="phc_waveguide", label="Photonic-crystal waveguide / slow light",
    summary="Broadband guided-mode emitter: slow-light group index and effective area (or beta "
            "directly), coupler and propagation efficiencies.",
    params=(
        _p("n_g", "1", 5.0, 1.0, 100.0, tag="E", source="slow-light group index class value [E]", help=r"Group index $n_{\mathrm g}$ of the guided mode."),
        _p("A_eff_um2", "um2", 0.03, 1e-3, 10.0, scale="log", tag="E", source="energy-normalised mode area class value [E]",
           help=r"Energy-normalised effective mode area $A_{\mathrm{eff}}$ at the emitter."),
        _p("beta_direct", "1", 0.0, 0.0, 0.999, help=r"Waveguide $\beta_{\mathrm{wg}}$ entered directly (0 = use $n_{\mathrm g}$ and $A_{\mathrm{eff}}$)."),
        _p("pos_factor", "1", 1.0, 0.0, 1.0, help=r"Position factor $\xi_{\mathrm r} = |E(\mathbf r)|^2/|E_{\max}|^2$ of the dot."),
        _p("eta_coupler", "1", 0.7, 0.0, 1.0, tag="E", source="grating-coupler class [E]", help=r"End-coupler efficiency $\eta_{\mathrm{coup}}$."),
        _p("eta_prop", "1", 0.95, 0.0, 1.0, tag="E", source="short-waveguide propagation [E]", help=r"Propagation transmission $\eta_{\mathrm{prop}}$ to the coupler."),
        _p("F_bg", "1", 1.0, 0.0, 10.0, help=r"Leaky / background decay rate $F_{\mathrm{bg}}$ in units of $\Gamma_0$ ($< 1$ inhibited)."),
        _p("pol_factor", "1", 1.0, 0.0, 1.0, help=r"Dipole-to-mode polarisation factor $\cos^2\theta_{\mathrm{pol}}$ (folded into $F_{\mathrm{wg}}$; ignored when $\beta_{\mathrm{wg}}$ is used)."),
    ),
    refs=(_REF_LECAMP,), eta_hint=(0.1, 0.8), viz="phc_waveguide", fn=_phc_waveguide)


# --------------------------------------------------------------------------
# photonic_nanowire
# --------------------------------------------------------------------------
# HE11 beta estimator correction (FX3a). The nitride_nanowire_photonics estimator is beta = G/(G + (1 - c)) with
# G = c n_g/n and c the Gaussian-core power fraction: it lets the off-core power radiate at the FREE rate, which for a
# thin high-index wire it does not. Emission into the leaky (radiation) modes of a wire with d < lambda/n is suppressed by
# dielectric screening of the transverse dipole, (2/(n^2 + 1))^2 in the quasi-static limit d -> 0 (Wang, Gudiksen, Duan,
# Cui and Lieber, Science 293, 1455 (2001); Bleuse et al., PRL 106, 103601 (2011) report the same strong inhibition
# of leaky modes) [DR], and the screening fades as the wire grows to ~lambda/n. With the unscreened estimator a GaAs wire
# peaks at beta ~0.83 instead of the ~0.95 of Bleuse 2011 / Claudon 2010 / Munsch 2013 theory. The leaky rate is therefore
# multiplied by kappa = s_q + (1 - s_q)(1 - exp(-(n d/lambda / X0)^2)), s_q = (2/(n^2 + 1))^2. The ONE free constant X0 is
# fixed so the peak of beta_HE11(d/lambda) for n = 3.45 equals the Bleuse 2011 peak ~0.95 (near d/lambda 0.24) [E]; it is
# NOT tuned to the Munsch 2013 card numbers (those are only checked afterwards). The interpolation shape is [A].
_NW_INH_X0 = 1.7                       # [E] calibrated to the Bleuse PRL 106, 103601 (2011) peak beta ~0.95 (n 3.45)


def nanowire_leak_inhibition(n_wire, d_over_lambda):
    """kappa in (s_q, 1): leaky-mode rate relative to the unscreened estimate (see the block comment above)."""
    s_q = (2.0 / (n_wire * n_wire + 1.0)) ** 2
    return s_q + (1.0 - s_q) * (1.0 - math.exp(-(n_wire * d_over_lambda / _NW_INH_X0) ** 2))


def _nanowire_resp(v, em, NA):
    lam = em.lambda_nm
    prm = _nwp.NitrideNanowirePhotonicsParams(
        family="vertical_photonic", NA=min(float(NA), 1.0), n_wire=v.pos("n_wire"),
        taper_transmission=v.frac("taper_T"), bottom_reflectivity=0.0,       # the mirror is handled coherently below
        radiative_rate_factor=v.pos("F_bg"),
        taper_output_mfr_nm=v.pos("taper_mfr_nm") if v.f("taper_mfr_nm") > 0 else None,
        n_group_override=v.pos("n_group") if v.f("n_group") > 0 else v.pos("n_wire"))
    return _nwp.response(prm, lambda_nm=lam, outer_radius_nm=0.5 * v.pos("d_nm"),
                         gamma_X0_ns=em.gamma0_ns, gamma_XX0_ns=em.gamma0_ns)


def _nanowire_beta(v, em, r):
    """beta_HE11 (both directions) with the screened leaky rate; returns (beta, beta_unscreened_module, kappa).

    Dipole orientation (audit finding 10): the HE11 field is transverse, so an in-plane dipole couples as above. An axial
    (growth-axis) dipole on the wire axis sees E_z = 0 of HE11 and radiates almost entirely into radiation modes
    (Bleuse et al., PRL 106, 103601 (2011) and Claudon et al., Nat. Photonics 4, 174 (2010) use in-plane dipoles [E]): the
    guided rate is the in-plane one times the projection AXIAL_RESIDUAL = eps [A] (G -> eps G) and the remaining radiative
    rate is the free one, unscreened (the dielectric screening acts on the TRANSVERSE dipole only):
    beta = eps G/(eps G + 1 - eps c) with kappa = 1 [DR/E]."""
    dg = r["diagnostics"]
    c, gr = dg["confinement_fraction"], dg["group_index_ratio"]
    kappa = nanowire_leak_inhibition(v.pos("n_wire"), v.pos("d_nm") / em.lambda_nm)
    G = c * gr
    if dipole_projection(em) < 1.0:
        eps = dipole_projection(em)
        den = eps * G + (1.0 - eps * c)
        b_sm = eps * G / den if den > 0 else 0.0
        return min(max(b_sm * dg["beta_multimode_penalty"], 0.0), 1.0), r["beta_HE11"], 1.0
    den = G + (1.0 - c) * kappa
    b_sm = G / den if den > 0 else 0.0
    return min(max(b_sm * dg["beta_multimode_penalty"], 0.0), 1.0), r["beta_HE11"], kappa


def _nanowire_eta_up(v, r):
    """Top-facet-to-lens efficiency of the upward guided photon: taper power x HE11 far-field acceptance at the NA."""
    return min(max(v.frac("taper_T") * r["diagnostics"]["obj_accept"], 0.0), 1.0)


def nanowire_mirror(beta, R, phase_ok):
    """(upward_fraction, rate_boost) for a bottom mirror of POWER reflectivity R = |r|^2 under an emitter on a field
    antinode: rate boost 1 + beta|r| (guided rate x (1 + |r|), leaky rate unchanged) and upward fraction
    (beta/2)(1 + |r|)^2/(1 + beta|r|) = Munsch, PRL 110, 177402 (2013), Eq. 1 without T_Gauss [V]. Without the antinode
    (phase_ok False) the phase is unknown: the incoherent power return, rate unchanged and (beta/2)(1 + R) [A]."""
    if phase_ok:
        r_abs = math.sqrt(R)
        return 0.5 * beta * (1.0 + r_abs) ** 2 / (1.0 + beta * r_abs), 1.0 + beta * r_abs
    return 0.5 * beta * (1.0 + R), 1.0


def _photonic_nanowire(self, vals, em, T_K, NA):
    v = _V(self, vals)
    r = _nanowire_resp(v, em, NA)
    beta, beta_mod, kappa = _nanowire_beta(v, em, r)
    R, phase_ok = v.frac("mirror_R"), bool(v("mirror_phase_ok"))
    up, boost = nanowire_mirror(beta, R, phase_ok)
    eta_out = _nanowire_eta_up(v, r)

    def fn(NA_):
        return _nanowire_eta_up(v, _nanowire_resp(v, em, NA_))

    notes = [_COMMON_NOTE_EFFNA,
             "beta0 = UPWARD collected fraction of the emission (the downward half only counts through the mirror), "
             "eta_out = taper_T x HE11 far-field acceptance at NA (nitride_nanowire_photonics obj_accept); F_P = 0 (no "
             "Purcell gain; the mirror's standing-wave boost enters the rate, inhibition via F_bg) [DR/E]",
             ("bottom mirror (emitter on a field antinode, mirror_phase_ok = true): power R = |r|^2 -> |r| = "
              f"{math.sqrt(R):.4f}; guided rate x (1 + |r|), total rate x (1 + beta|r|) = {boost:.4f}; upward fraction "
              f"(beta/2)(1 + |r|)^2/(1 + beta|r|) = {up:.4f} (Munsch 2013 Eq. 1 without T_Gauss) [V/DR]" if phase_ok else
              f"bottom mirror, mirror_phase_ok = false: dot-mirror phase unknown -> incoherent power return, rate unchanged, "
              f"upward fraction (beta/2)(1 + R) = {up:.4f} [A]"),
             "taper: taper_T is a pure POWER transmission; the far-field effect of an adiabatic taper (the guided mode "
             "expands to a larger mode-field radius at the top facet and the emission narrows) enters only through "
             "taper_mfr_nm. With taper_mfr_nm = 0 the acceptance is that of the BARE wire facet (mode radius = Marcuse w, "
             "~0.5 at NA 0.75 for d/lambda 0.23), which under-collects a tapered / trumpet device [DR]",
             (f"material dipole = axial: HE11 coupling x {AXIAL_RESIDUAL:g}, unscreened free leaky rate (beta much smaller than in-plane) [DR/E]"
              if dipole_projection(em) < 1.0 else "material dipole = in_plane (couples to the transverse HE11 field)"),
             f"V = {r['V_number']:.3f}, d/lambda = {2 * r['radius_over_lambda']:.3f}, beta_HE11 = {beta:.3f} (module "
             f"estimator without leaky-mode screening: {beta_mod:.3f}; screening kappa = {kappa:.3f} on the leaky rate, "
             f"X0 = {_NW_INH_X0} calibrated to the Bleuse 2011 peak ~0.95 [E])",
             "the Marcuse / multimode envelopes are calibrated on GaAs/InAs-class wires [A] (Bleuse 2011 class)"]
    if not r["valid"]:
        notes.append("module reports invalid: " + "; ".join(r["invalid_reasons"]))
    if float(NA) > 1.0:
        notes.append("NA > 1 clipped to 1 (the HE11 acceptance integral saturates at the forward hemisphere)")
    inputs = {"beta0": _path("d_nm") + "+n_wire+n_group+mirror_R+mirror_phase_ok", "eta_out": _path("taper_T") + "+taper_mfr_nm",
              "F_bg": _path("F_bg") + "+mirror_R"}
    # geometry numbers for the 3D scene (R2-3D): read from the module's own diagnostics, no new physics
    extra = {"V_number": float(r["V_number"]), "mode_field_radius_nm": float(r["diagnostics"]["mode_field_radius_nm"]),
             "beta_he11": float(beta), "upward_fraction": float(up), "rate_boost": float(boost),
             "mirror_field_r": math.sqrt(R), "phase_ok": bool(phase_ok)}
    return StructureResponseX(F_P=0.0, kappa_meV=None, E_mode_meV=None, g_meV=None, F_bg=v.pos("F_bg") * boost,
                              beta0=float(up), eta_out=float(eta_out), eta_out_fn=fn,
                              pol_factor=1.0, notes=tuple(notes), inputs=inputs,
                              tag_floor="A", extra=extra)       # Marcuse / multimode envelopes calibrated on GaAs/InAs wires [A]


PHOTONIC_NANOWIRE = _Struct(
    kind="photonic_nanowire", label="Photonic nanowire / trumpet",
    summary=r"Single-mode dielectric nanowire with a tapered output: HE11 coupling $\beta$, taper and mirror efficiencies, far-field acceptance at the setup NA. The HE11 field is transverse: an axial dipole couples poorly.",
    params=(
        _p("d_nm", "nm", 220.0, 50.0, 2000.0, tag="E", source="Bleuse 2011 class: beta_HE11 peaks near d/lambda 0.22-0.25 [E]", help=r"Wire diameter $d$."),
        _p("n_wire", "1", 3.45, 1.5, 4.5, tag="E", source="GaAs near 930 nm [E]", help=r"Refractive index $n_{\mathrm{wire}}$ of the wire."),
        _p("n_group", "1", 3.6, 0.0, 5.0, tag="E", source="GaAs group index near 930 nm [E]; 0 = use n_wire (group-index ratio 1, [A])",
           help=r"Group index $n_{\mathrm{g,wire}}$ of the wire material (needed by the beta estimator for a non-GaN wire).", group="advanced"),
        _p("taper_T", "1", 0.95, 0.0, 1.0, tag="E", source="Bleuse 2011 / Claudon 2010 class 0.9-0.95 [E]", help=r"Taper power transmission $T_{\mathrm{taper}}$."),
        _p("mirror_R", "1", 0.9, 0.0, 1.0, tag="E", source="bottom-mirror class [E]; Munsch 2013 |r| = 0.95 -> R = 0.9025 [V]",
           help=r"Bottom-mirror POWER reflectivity $R_{\mathrm m} = |r|^2$ (the field reflectivity $|r| = \sqrt{R_{\mathrm m}}$ sets the coherent boost)."),
        _p("mirror_phase_ok", "bool", False, tag="A", source="dot-mirror phase unknown unless the paper states antinode placement: incoherent power return [A] (FX4; Munsch 2013 designs an antinode, set True for such devices)",
           help=r"True: the dot sits on a standing-wave antinode (rate $\times (1 + \beta|r|)$, upward fraction $(\beta/2)(1 + |r|)^2/(1 + \beta|r|)$). False: phase unknown, incoherent power return $(1 + R_{\mathrm m})/2$.", sweepable=False),
        _p("F_bg", "1", 1.0, 0.0, 10.0, tag="E", source="Gamma_HE11/Gamma_0 ~ 1; < 1 for leaky-mode inhibition [E]",
           help=r"Decay-rate factor $F_{\mathrm{bg}}$ in units of $\Gamma_0$ ($< 1$ = inhibition of leaky modes)."),
        _p("taper_mfr_nm", "nm", 0.0, 0.0, 3000.0, tag="A", source="0 = no mode expansion (bare-wire facet) [A]",
           help=r"Mode-field radius $w_{\mathrm{taper}}$ at the top facet after an adiabatic taper (0 = none: the emission leaves with the bare-wire Marcuse mode radius). A larger radius narrows the far field and raises the NA acceptance; $T_{\mathrm{taper}}$ stays the separate power loss.", group="advanced"),
    ),
    refs=(_REF_BLEUSE,), eta_hint=(0.2, 0.8), fn=_photonic_nanowire)


# --------------------------------------------------------------------------
# planar_dbr
# --------------------------------------------------------------------------
def planar_dbr_calc(v, em):
    """Raw dbr calls for the planar_dbr entry (also used by the verifier)."""
    lam = em.lambda_nm
    n_h, n_l, n_c, n_s = v.pos("n_high"), v.pos("n_low"), v.pos("n_cavity"), v.pos("n_substrate")
    Nt, Nb = int(v.f("N_top")), int(v.f("N_bot"))
    top, spacer, bottom = dbr.cavity_stack(n_h, n_l, n_c, Nt, Nb, lam)
    flat = top + [spacer] + bottom
    mode = dbr.cavity_mode(flat, lam, n_in=1.0, n_out=n_s)
    return top, spacer, bottom, mode, n_c, n_s, lam


def _planar_dbr(self, vals, em, T_K, NA):
    v = _V(self, vals)
    top, spacer, bottom, mode, n_c, n_s, lam = planar_dbr_calc(v, em)
    d_c = lam / n_c
    z = v.frac("z_dot_frac") * d_c
    F_1d = dbr.planar_purcell(top, bottom, n_c, d_c, z, lam, n_in=1.0, n_out=n_s)
    orient = v("orientation")
    if str(getattr(em, "dipole", "in_plane")) == "axial" and orient == "inplane":
        orient = "vertical"                      # material dipole axial = growth-axis dipole [DR]
    F_tot = dbr.planar_total_rate(top, bottom, n_c, d_c, z, lam, n_in=1.0, n_out=n_s, orientation=orient)
    n_h, n_l = v.pos("n_high"), v.pos("n_low")
    Tt = 1.0 - dbr.qw_peak_reflectivity(1.0, n_c, n_h, n_l, int(v.f("N_top")))
    Tb = 1.0 - dbr.qw_peak_reflectivity(n_s, n_c, n_h, n_l, int(v.f("N_bot")))
    kappa = mode["kappa_meV"]
    kappa_ok = kappa is not None and math.isfinite(kappa)
    notes = ["planar DBR microcavity: 1-D transfer matrices, no lateral confinement, no mode volume [E->analytic-1D]",
             f"F_P = 0: dbr.planar_purcell (on-axis 1-D LDOS = {F_1d:.3g}) is NOT a cavity-mode F_P and the total rate is "
             "not F_P either (dbr honesty wall, audit H5); the angle-integrated total rate planar_total_rate = "
             f"{F_tot:.3g} x Gamma0 is returned as F_bg [DR]",
             "kappa, E_mode from dbr.cavity_mode (informational: with F_P = 0 the mode does not change the rate)",
             "eta_out = T_top/(T_top + T_bot) from the quarter-wave stack peak reflectivities; no NA model "
             "(planar emission is not Gaussian), so eta_out is independent of NA",
             "beta0 must be entered: the 1-D model cannot give the collected cavity-cone fraction (0 gives zero brightness)"]
    if orient != v("orientation"):
        notes.append("material dipole = axial: the total rate uses the vertical-dipole LDOS instead of the in-plane one [DR]")
    if not kappa_ok:
        notes.append("cavity_mode found no resolvable resonance (kappa NaN); kappa returned as None")
    inputs = {"F_bg": _path("z_dot_frac"), "E_mode_meV": _path("N_top"), "kappa_meV": _path("N_bot"),
              "eta_out": _path("N_top"), "beta0": _path("beta0")}
    # geometry numbers for the 3D scene (R2-3D): the same quarter-wave stack the transfer matrix used
    extra = {"t_high_nm": float(lam / (4.0 * n_h)), "t_low_nm": float(lam / (4.0 * n_l)), "spacer_nm": float(d_c),
             "z_dot_nm": float(z), "lambda_nm": float(lam)}
    return StructureResponseX(F_P=0.0, kappa_meV=float(kappa) if kappa_ok else None,
                              E_mode_meV=float(mode["E_cav_meV"]), g_meV=None, F_bg=float(F_tot),
                              beta0=v.frac("beta0"), eta_out=float(Tt / (Tt + Tb)), eta_out_fn=None,
                              pol_factor=1.0, notes=tuple(notes), inputs=inputs, tag_floor="E", extra=extra)


PLANAR_DBR = _Struct(
    kind="planar_dbr", label="Planar DBR microcavity",
    summary="Unetched lambda-cavity between two DBRs (1-D transfer matrix): mode position and linewidth, "
            "total emission-rate change, mirror asymmetry.",
    params=(
        _p("N_top", "1", 20, 1, 60, help=r"Top DBR pairs $N_{\mathrm{top}}$."),
        _p("N_bot", "1", 30, 1, 60, help=r"Bottom DBR pairs $N_{\mathrm{bot}}$."),
        _p("n_high", "1", 3.5, 2.0, 4.0, tag="E", source="GaAs near 930 nm [E]", help=r"High-index layer index $n_{\mathrm H}$."),
        _p("n_low", "1", 2.95, 1.2, 3.5, tag="E", source="AlAs near 930 nm [E]", help=r"Low-index layer index $n_{\mathrm L}$."),
        _p("n_cavity", "1", 3.5, 1.5, 4.0, tag="E", source="GaAs spacer [E]", help=r"Spacer index $n_{\mathrm c}$ ($\lambda$ cavity)."),
        _p("n_substrate", "1", 3.5, 1.0, 4.0, tag="E", source="GaAs substrate [E]", help=r"Exit-medium index $n_{\mathrm{sub}}$ below the bottom DBR."),
        _p("z_dot_frac", "1", 0.5, 0.0, 1.0, help=r"Dot depth $z_{\mathrm{dot}}$ in the spacer as a fraction of its thickness $d_{\mathrm c}$ (0.5 = antinode)."),
        _p("orientation", "enum", "inplane", choices=("inplane", "vertical", "isotropic"), help=r"Dipole orientation $\hat{\boldsymbol\mu}$ for the total rate (an axial material dipole overrides in-plane with vertical)."),
        _p("beta0", "1", 0.0, 0.0, 1.0, help=r"Fraction $\beta_0$ of emission collected into the cavity cone (not computable in 1-D)."),
    ),
    refs=("[DR] Bjork et al., Phys. Rev. A 44, 669 (1991); Benisty et al., IEEE JQE 34, 1612 (1998): planar-cavity total rate",),
    eta_hint=(0.01, 0.5), fn=_planar_dbr)


# --------------------------------------------------------------------------
# bulk_planar
# --------------------------------------------------------------------------
def _bulk_planar(self, vals, em, T_K, NA):
    v = _V(self, vals)
    n, n_s = v.pos("n"), v.pos("n_s")
    mirror, sil = bool(v("mirror")), v("sil")
    dip = str(getattr(em, "dipole", "in_plane"))
    fn = lambda NA_: bulk_planar_eta(NA_, n, mirror, sil, n_s, dip)       # noqa: E731
    mult = {"none": 1.0, "hemi": n_s, "weierstrass": n_s ** 2}[sil]
    notes = [_COMMON_NOTE_EFFNA,
             "eta(NA) = (3/8) (NA_eff/n)^2 T_F, T_F = 4n/(1 + n)^2 (in-plane dipole, small cone) [DR]; bottom mirror x2; "
             "SIL: NA_eff = NA x n_s (hemispherical) or n_s^2 (Weierstrass), NA_eff/n capped at 1",
             ("material dipole = axial: sin^2(theta) pattern, no emission on axis, eta = (3/16)(NA_eff/n)^4 T_F [DR]" if dip == "axial"
              else "material dipole = in_plane (the formula above)"),
             "no cavity: F_P = 0, beta0 = 1 and eta_out is the collected fraction of all emission; rate_factor "
             "(planar interference, 0.5-1.5 [E]) enters as F_bg"]
    if NA * mult / n > 1.0:
        notes.append("NA_eff/n capped at 1 (full forward hemisphere of the host)")
    inputs = {"eta_out": _path("n") + "+mirror+sil+n_s", "F_bg": _path("rate_factor")}
    return StructureResponse(F_P=0.0, kappa_meV=None, E_mode_meV=None, g_meV=None, F_bg=v.pos("rate_factor"),
                             beta0=1.0, eta_out=fn(NA), eta_out_fn=fn, pol_factor=1.0,
                             notes=tuple(notes), inputs=inputs)


BULK_PLANAR = _Struct(
    kind="bulk_planar", label="Bulk / planar (optional mirror and SIL)",
    summary="Dot in a planar high-index host with no cavity: collected fraction from the escape cone, "
            "optionally with a back mirror and a solid-immersion lens.",
    params=(
        _p("n", "1", 3.5, 1.5, 4.5, tag="E", source="GaAs near 930 nm [E]", help=r"Host refractive index $n$."),
        _p("mirror", "", False, help=r"Back mirror $\mathrm{mirror}$ (doubles the upward collection)."),
        _p("sil", "enum", "none", choices=("none", "hemi", "weierstrass"), help=r"Solid-immersion lens $\mathrm{SIL}$ type."),
        _p("n_s", "1", 2.0, 1.0, 4.5, tag="A", source="SIL index [A]", help=r"SIL refractive index $n_{\mathrm{SIL}}$.", visible_if={"sil": ["hemi", "weierstrass"]}),
        _p("rate_factor", "1", 1.0, 0.5, 1.5, tag="E", source="planar interference factor 0.5-1.5 [E]", help=r"Total emission-rate factor $F_{\mathrm{rate}}$ (enters as $F_{\mathrm{bg}}$)."),
    ),
    refs=("[DR] small-cone escape fraction (3/8)(NA/n)^2 T_F (Schweickert et al. planar low-Q class, [V-abs])",),
    eta_hint=(0.005, 0.1), viz="mesa", fn=_bulk_planar)


# --------------------------------------------------------------------------
# sil_mesa / mesa_reflector
# --------------------------------------------------------------------------
def _table_struct(self, vals, em, T_K, NA):
    v = _V(self, vals)
    fn = eta_table_fn(v("eta_table"))
    notes = [_COMMON_NOTE_EFFNA,
             "no cavity: F_P = 0, beta0 = 1; eta_out(NA) is an entered table [E] (Vajner-class 10-30 %), piecewise linear, "
             "flat above the last point, eta0 (NA/NA0)^2 below the first",
             "rate_factor enters as F_bg"]
    inputs = {"eta_out": _path("eta_table"), "F_bg": _path("rate_factor"), "beta0": _path("eta_table")}
    return StructureResponse(F_P=0.0, kappa_meV=None, E_mode_meV=None, g_meV=None, F_bg=v.pos("rate_factor"),
                             beta0=1.0, eta_out=fn(NA), eta_out_fn=fn, pol_factor=1.0,
                             notes=tuple(notes), inputs=inputs)


_VAJNER = ("[E] Vajner et al., review of QD single-photon-source structures (class range 10-30 % for mesa + reflector / SIL); "
           "the exact figures are not checked against the paper")
SIL_MESA = _Struct(
    kind="sil_mesa", label="Mesa with solid-immersion lens",
    summary="Etched mesa under a solid-immersion lens: no Purcell gain, collection efficiency from an entered NA table.",
    params=(
        _p("eta_table", "", "0.7:0.3", tag="E", source=_VAJNER, help=r"Collected fraction $\eta(\mathrm{NA})$ versus NA, 'NA:eta[, NA:eta ...]'."),
        _p("rate_factor", "1", 1.0, 0.5, 1.5, tag="E", source="F_P ~ 1 class [E]", help=r"Total emission-rate factor $F_{\mathrm{rate}}$ (enters as $F_{\mathrm{bg}}$)."),
    ),
    refs=(_VAJNER,), eta_hint=(0.1, 0.4), viz="mesa", fn=_table_struct)

MESA_REFLECTOR = _Struct(
    kind="mesa_reflector", label="Mesa with metal reflector",
    summary="Etched mesa on a metal back reflector: no Purcell gain, collection efficiency from an entered NA table.",
    params=(
        _p("eta_table", "", "0.7:0.2", tag="E", source=_VAJNER, help=r"Collected fraction $\eta(\mathrm{NA})$ versus NA, 'NA:eta[, NA:eta ...]'."),
        _p("rate_factor", "1", 1.0, 0.5, 1.5, tag="E", source="F_P ~ 1 class [E]", help=r"Total emission-rate factor $F_{\mathrm{rate}}$ (enters as $F_{\mathrm{bg}}$)."),
    ),
    refs=(_VAJNER,), eta_hint=(0.05, 0.3), viz="mesa", fn=_table_struct)


# --------------------------------------------------------------------------
# edge_ridge
# --------------------------------------------------------------------------
def edge_ridge_stack(v, em):
    """(stack, family_key, family_label, lambda_nm, reason): the ridge platform for the stack wavelength
    (waveguide.ridge_stack_for). stack is None and reason names the incompatibility when no platform covers it."""
    lam = em.lambda_nm if v.f("lambda_nm") <= 0 else v.f("lambda_nm")
    stack, key, label = waveguide.ridge_stack_for(lam)
    if stack is None:
        return None, None, None, lam, label           # `label` carries the reason in this case
    return stack, key, label, lam, None


def edge_ridge_incompatibility(em, vals=None):
    """Reason string when edge_ridge cannot be built for this emitter (greyed in the picker), else None. Pure function
    of the stack wavelength (the structure's lambda_nm override, else the emitter's)."""
    lam = em.lambda_nm
    if vals:
        try:
            lo = float(vals.get("lambda_nm", 0.0))
            lam = lo if lo > 0 else lam
        except (TypeError, ValueError):
            pass
    return waveguide.ridge_family(lam)[2]


def edge_ridge_args(v, em):
    """(stack, ridge_width_nm, etch_depth_nm, lambda_nm, L_um, kwargs) for waveguide.edge_emission. stack is None when
    the wavelength has no ridge platform (see edge_ridge_stack)."""
    stack, _key, _label, lam, _why = edge_ridge_stack(v, em)
    rb = v.f("R_back")
    kw = dict(alpha_cm=v.f("alpha_cm"), R_back=None if rb < 0 else rb, coating=None,
              na_method=v("na_method"), dot_position=v.f("dot_position"))
    return stack, v.pos("ridge_width_nm"), v.pos("etch_depth_nm"), lam, v.f("L_um"), kw


def _edge_ridge_incompatible(reason, lam):
    """Neutral, finite response for a (material, edge_ridge) pair that cannot be built: no Purcell gain, nothing collected.
    `extra["incompatible"]` carries the reason for the UI; the pair is declared incompatible, never an exception."""
    return StructureResponseX(F_P=0.0, kappa_meV=None, E_mode_meV=None, g_meV=None, F_bg=1.0, beta0=0.0, eta_out=0.0,
                              eta_out_fn=None, pol_factor=1.0, notes=("INCOMPATIBLE: " + reason,), inputs={},
                              tag_floor="E", extra={"incompatible": reason, "lambda_nm": float(lam)})


def _edge_ridge(self, vals, em, T_K, NA):
    v = _V(self, vals)
    stack, key, family, lam, why = edge_ridge_stack(v, em)
    if stack is None:
        return _edge_ridge_incompatible(why, lam)
    _, w, etch, lam, L, kw = edge_ridge_args(v, em)
    try:
        e = waveguide.edge_emission(stack, w, etch, lam, L, NA, **kw)
        mode = waveguide.effective_index_ridge(stack, lam, w, etch)
    except ValueError as err:                       # no bound mode / etch cannot be applied: a named incompatibility
        return _edge_ridge_incompatible(f"{family}: the ridge mode solve failed ({err})", lam)
    dp = dipole_projection(em)
    eta_out = e.eta_total / e.beta if e.beta > 0 else 0.0
    notes = ["adapter around waveguide.edge_emission: F_P = F_wg (Lecamp, energy-normalised area), beta0 = 0 so "
             "beta = F_wg/(1 + F_wg); eta_out = eta_total/beta = facet escape fraction x NA collection [DR]",
             f"broadband ridge: kappa = None (L = 1); stack is waveguide.ridge_stack_for({lam:g} nm) = {family} [E] "
             "(red InP: waveguide.hkust_ridge_stack; GaAs / InP / GaN platforms from the materials tables)",
             "NA is applied inside edge_emission (Gaussian or numeric); no eta_out_fn (each NA needs a mode solve)"]
    if dp < 1.0:
        notes.append(f"material dipole {em.dipole}: cos^2 projection {dp:.4g} on the in-plane TE ridge mode folded into F_P [DR/A]")
    notes.extend(e.notes)
    # F_bg = 1 and beta0 = 0 are fixed by the model; F_P / eta_out inputs are the leaves that change them when perturbed
    inputs = {"F_P": "+".join(_path(k) for k in ("ridge_width_nm", "etch_depth_nm", "lambda_nm")),
              "eta_out": "+".join(_path(k) for k in ("ridge_width_nm", "etch_depth_nm", "L_um", "alpha_cm", "R_back",
                                                      "lambda_nm", "dot_position", "na_method"))}
    # geometry numbers for the 3D scene (R2-3D): the stack the mode solve used and the mode's 1/e^2 intensity radii
    extra = {"stack": [{"name": L_.name, "t_nm": float(L_.thickness_nm), "n": float(L_.n), "is_dot": bool(L_.is_dot),
                        "material": L_.material or ""} for L_ in stack],
             "n_eff": float(e.n_eff), "n_g": float(e.n_g), "T_facet": float(e.T_facet),
             "wx_um": float(mode.wx_um), "wy_um": float(mode.wy_um), "lambda_nm": float(lam),
             "stack_family": key, "stack_label": family}
    return StructureResponseX(F_P=float(e.F_wg) * dp, kappa_meV=None, E_mode_meV=None, g_meV=None, F_bg=1.0, beta0=0.0,
                              eta_out=float(eta_out), eta_out_fn=None, pol_factor=1.0, notes=tuple(notes), inputs=inputs,
                              tag_floor="E", extra=extra)       # every ridge stack is an [E] class design


EDGE_RIDGE = _Struct(
    kind="edge_ridge", label="Edge-emitting ridge waveguide",
    summary=r"Dot in a ridge waveguide emitting from the cleaved facet: guided-mode Purcell $F_{\mathrm P}$, facet escape and NA collection. The stack (AlGaInP, GaAs/AlGaAs, InP/InGaAsP, GaN/AlGaN) follows the emission wavelength; other wavelengths are incompatible.",
    params=(
        _p("ridge_width_nm", "nm", 2000.0, 200.0, 10000.0, tag="A", source="class default, edge-emitter design [A]", help=r"Ridge width $w_{\mathrm{ridge}}$."),
        _p("etch_depth_nm", "nm", 1200.0, 50.0, 5000.0, tag="A", source="class default, edge-emitter design [A]", help=r"Ridge etch depth $h_{\mathrm{etch}}$ (from the top)."),
        _p("L_um", "um", 250.0, 10.0, 2000.0, tag="A", source="class default, edge-emitter design [A]", help=r"Mirror separation $L$ (air gap)."),
        _p("alpha_cm", "1/cm", 5.0, 0.0, 100.0, tag="E", source="ridge propagation loss class estimate [E]", help=r"Propagation loss $\alpha_{\mathrm{loss}}$."),
        _p("R_back", "1", 0.95, -1.0, 1.0, tag="A", source="class default, edge-emitter design [A]",
           help=r"Back-facet reflectivity $R_{\mathrm{back}}$; negative = bare cleaved facet (waveguide.edge_emission R_back=None)."),
        _p("lambda_nm", "nm", 0.0, 0.0, 2000.0, help=r"Wavelength $\lambda$ override for the stack; 0 = the emitter's wavelength. The stack platform (AlGaInP, GaAs/AlGaAs, InP/InGaAsP or GaN/AlGaN) follows $\lambda$; outside the covered bands the pair is declared incompatible.", group="advanced"),
        _p("na_method", "enum", "gaussian", choices=("gaussian", "numeric"), help=r"Far-field $\mathrm{NA}$ collection method.", group="advanced"),
        _p("dot_position", "1", 0.5, 0.0, 1.0, help=r"Dot position $x_{\mathrm{dot}}/L$ along the ridge (fraction).", group="advanced"),
    ),
    refs=(_REF_LECAMP, "[DR] Coldren, Corzine & Masanovic, Diode Lasers and Photonic Integrated Circuits, 2nd ed. (2012): facet escape"),
    eta_hint=(0.01, 0.4), fn=_edge_ridge)


# --------------------------------------------------------------------------
# nitride_planar_cavity
# --------------------------------------------------------------------------
def nitride_cavity_args(v, em, T_K):
    prm = nitride_cavity.NitrideCavityParams(
        Q=v.f("Q"), mode_volume_norm=v.f("mode_volume_norm"), spatial_overlap=v.f("spatial_overlap"),
        eta_out=v.f("eta_out"), T_track=v.f("T_track"), dEdT_cav_meV_K=v.f("dEdT_cav_meV_K"),
        detuning_offset_meV=v.f("detuning_offset_meV"), purcell_enabled=bool(v("purcell_enabled")))
    ex = em.E_X_meV * 1e-3
    track = v.f("E_X_track_eV") if v.f("E_X_track_eV") > 0 else ex
    # em.delta_xx_meV: NEGATIVE = bound; nitride_cavity.response wants delta_xx POSITIVE = binding (E_XX = E_X - delta_xx)
    gam = HBAR_MEV_NS * (em.gamma0_ns + em.k_nr_ns + em.gamma_star_ns)
    kw = dict(T_K=T_K, E_X_eV=ex, E_X_track_eV=track, gamma_X_meV=gam, gamma_XX_meV=gam * em.gamma_xx_ratio,
              delta_xx_meV=-em.delta_xx_meV, gamma_X0_ns=em.gamma0_ns, gamma_XX0_ns=em.gamma0_ns)
    return prm, kw


def _nitride_planar_cavity(self, vals, em, T_K, NA):
    v = _V(self, vals)
    prm, kw = nitride_cavity_args(v, em, T_K)
    r = nitride_cavity.response(prm, **kw)
    notes = ["adapter around nitride_cavity.response: F_P = Fp_add (additive cavity channel), kappa = 1000 E_cav/Q, "
             "E_mode = tracked cavity energy (anchored at E_X_track_eV, slope dEdT_cav) [DR/A]",
             "emitter linewidth for the response = hbar (Gamma0 + k_nr + gamma*); the XX line uses em.gamma_xx_ratio and "
             "Gamma0 for its bare rate (only t_X/t_XX use them, not the returned F_P) [A]",
             "eta_out is the card's [A] outcoupling; no NA dependence"]
    inputs = {"F_P": _path("Q"), "kappa_meV": _path("Q"), "E_mode_meV": _path("dEdT_cav_meV_K"), "eta_out": _path("eta_out"),
              "F_bg": _path("purcell_enabled")}
    dp = dipole_projection(em)
    if dp < 1.0:
        notes.append(f"material dipole {em.dipole}: cos^2 projection {dp:.4g} on the in-plane cavity field folded into F_P [DR/A]")
    return _cavity_resp(self, em, F_P=r["Fp_add"] * dp, kappa_meV=r["kappa_meV"], detuning_meV=em.E_X_meV - 1000.0 * r["E_cav_eV"],
                        eta_out=r["eta_out"], notes=notes, inputs=inputs)


NITRIDE_PLANAR_CAVITY = _Struct(
    kind="nitride_planar_cavity", label="Nitride planar / micropillar cavity",
    summary="InGaN/GaN dot in a vertical cavity with a tracked mode (nitride_cavity.response adapter).",
    params=(
        _p("Q", "1", 2000.0, 10.0, 1e5, scale="log", tag="A", source="class default, nitride cavity design [A]", help=r"Quality factor $Q$."),
        _p("mode_volume_norm", "1", 2.0, 0.3, 50.0, tag="A", source="class default, mode volume in (lambda/n)^3 [A]", help=r"Mode volume $V$ in $(\lambda/n)^3$."),
        _p("spatial_overlap", "1", 1.0, 0.0, 1.0, tag="A", source="class default, nitride cavity design [A]", help=r"Dot-field overlap $\xi$."),
        _p("eta_out", "1", 0.1, 0.0, 1.0, tag="A", source="class default, nitride cavity design [A]", help=r"Outcoupling efficiency $\eta_{\mathrm{out}}$."),
        _p("T_track", "K", 300.0, 4.0, 500.0, tag="A", source="class default, nitride cavity design [A]", help=r"Temperature $T_{\mathrm{track}}$ at which the mode is aligned with the dot."),
        _p("dEdT_cav_meV_K", "meV/K", -0.04, -0.2, 0.0, tag="A", source="class default, nitride cavity design [A]", help=r"Cavity-mode temperature slope $\mathrm dE_{\mathrm c}/\mathrm dT$."),
        _p("detuning_offset_meV", "meV", 0.0, -50.0, 50.0, tag="A", source="class default, nitride cavity design [A]", help=r"Static mode offset $\delta E_{\mathrm{off}}$ below the tracked dot energy."),
        _p("purcell_enabled", "", True, help=r"False = cavity filtering without added Purcell ($F_{\mathrm P} = 0$)."),
        _p("E_X_track_eV", "eV", 0.0, 0.0, 6.0, help=r"Dot energy $E_X(T_{\mathrm{track}})$ at $T_{\mathrm{track}}$ (0 = the emitter's $E_X$ at the run temperature).", group="advanced"),
    ),
    refs=("[DR] Purcell, Phys. Rev. 69, 681 (1946); geometry, overlap, eta_out and slope [A]",), eta_hint=(0.05, 0.3), fn=_nitride_planar_cavity)


for _e in (GENERIC_CAVITY, MICROPILLAR, OPEN_CAVITY, CBG, PHC_CAVITY, PHC_WAVEGUIDE, PHOTONIC_NANOWIRE,
           PLANAR_DBR, BULK_PLANAR, SIL_MESA, MESA_REFLECTOR, EDGE_RIDGE, NITRIDE_PLANAR_CAVITY):
    register(_e)


# Declared (material, structure) incompatibilities (audit finding 3): kind -> fn(em, vals) -> reason or None. The reason is
# shown by the Composer picker (greyed); the entry's response() returns a finite, named-incompatible response instead of
# raising (extra["incompatible"]). Every other pair is expected to be finite.
STRUCTURE_COMPAT = {"edge_ridge": edge_ridge_incompatibility}


def incompatibility(kind, em, vals=None):
    """Reason string if the structure `kind` cannot be built for the emitter `em`, else None."""
    fn = STRUCTURE_COMPAT.get(kind)
    return fn(em, vals) if fn else None
