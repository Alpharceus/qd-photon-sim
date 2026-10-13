"""FOUR VALIDATION CLASSES (full definitions and the table sorting every
check into one of them: README.md "Validation record"):
  (N) numerical verification -- one method checked against another under
      the same physical assumptions, with no external data involved.
  (T) source transcription   -- a number, formula, or claim in the code
      checked against the paper that states it.
  (C) parameter calibration  -- free parameters fit to data and then
      evaluated on that same data. Not held-out validation.
  (P) held-out prediction    -- a fitted or class-range model checked
      against data it did not see during fitting or class-range selection.
      Currently empty across this repository.

This module performs a (C) parameter calibration -- the V-a fit below
(Chatzarakis joint fit) is evaluated on the same data it was fit to; it is
not, and is not presented as, held-out validation.

Phase-0/1 fit driver: joint fit of {eps(T), rho(T)} to the validation card.

Residual blocks (each optional block activates when the card carries it):
  1. g2_vs_T      -- the published g2(0) series; rows carry {w, dx} filter
                     geometry and an optional bound:"upper" flag (one-sided).
  2. gamma_vs_T   -- digitized Gamma_X(T) [and Gamma_XX(T)] linewidths (S2);
                     constrains the Gamma model and the XX/X width ratio r_xx.
  3. tau_vs_T     -- published X decay times (S4). tau(T)/tau(T0) is an
                     independent measure of the retention ratio S(T)/S(T0)
                     (assumes Gamma_r const -- assumption ledger).
Anchors/priors: Gamma(250 K) anchor (only when no gamma_vs_T data), E_a prior,
log10 b_p prior. V-a criterion: +-0.03 per g2 point (planning doc section 1).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.interpolate import CubicSpline
from scipy.optimize import least_squares

from .card import Card, Tag, widest
from .integrator import g2_of_T, g2_from, retention
from .loading import f1b_g2
from .spectral import gamma_of_T

_ROOT = Path(__file__).resolve().parents[1]
UPPER_BOUND_ERR_DIVISOR = 3.0  # [A] IBM bundle: a published inequality is weighted with err/3, not +-err (r2-physics-a Q3)

BASE_FREE = ("gamma0", "a_ac", "b_lo", "a_esc", "E_a", "b_p", "b0", "beta")
LOG_FREE = {"a_esc", "b_p", "b0", "beta"}  # span decades: fitted in log10 space

V_A_TOL = 0.03          # per-point absolute tolerance on g2
GAMMA_RELERR = 0.15     # digitization-class relative error on linewidths
TAU_RELERR = 0.10       # relative error on lifetime ratios


@dataclass
class FitResult:
    params: dict                 # full flat parameter dict (fitted + fixed)
    data: list[dict]             # g2 rows: T, g2, err, w, dx, bound
    model_g2: np.ndarray
    residuals: np.ndarray        # g2 block: model - data (0 for satisfied bounds)
    passed: bool                 # all g2 |resid| <= V_A_TOL
    gamma_anchor: tuple          # (model Gamma(T_anchor), anchor, tol, T_anchor)
    tag: Tag
    cost: float = 0.0
    notes: list[str] = field(default_factory=list)
    aux: dict = field(default_factory=dict)   # per-block worst residuals etc.


# ---------------------------------------------------------------- IBM (sideband-aware) lineshape refit
# r2-physics-a Q3 / devlib-r2-pha2 (D1).  The Lorentzian fit absorbs the XX line's sideband leakage into a
# T-independent b0; the IBM line carries it physically.  The measured Gamma(T) (Fig. S2) is a TOTAL FWHM, so the ZPL
# width fed to the IBM line is the root of FWHM_IBM(gamma_zpl, T) = Gamma(T) -- the full width is NOT used as the ZPL
# width (that double-counts the phonon broadening, device.py:2154 before this change).  Opt-in only.

_IBM_GRID_LO, _IBM_GRID_HI, _IBM_GRID_N = 0.01, 14.0, 70   # meV, log-spaced ZPL-width table [A]


def chatz2023_phonon() -> dict:
    """PhononParams kwargs of the chatzarakis2023 device card (geometry + Ramsay-class alpha) [DR/A]."""
    import yaml
    dev = yaml.safe_load((_ROOT / "cards" / "chatzarakis2023.yaml").read_text(encoding="utf-8"))["device"]
    return {**dev["geometry"], "alpha_ps2": float(dev["phonon_alpha_ps2"])}


def ibm_fwhm(pp, T, gamma_zpl_meV):
    """Numeric FWHM [meV] of the IBM emission line (ZPL Lorentzian of FWHM gamma_zpl + sidebands): the contiguous
    region around the peak where S >= peak/2.  Not normalised (the half-maximum is scale free)."""
    from .qd_gf import _spectrum_parts
    omega, s_sb, z = _spectrum_parts(pp, float(T), float(gamma_zpl_meV))
    g = float(gamma_zpl_meV)
    grid = np.arange(-40.0, 40.0 + 1e-9, 0.002)
    s = z * (g / (2.0 * np.pi)) / (grid**2 + (g / 2.0) ** 2) + np.interp(grid, omega, s_sb, left=0.0, right=0.0)
    i0 = int(np.argmax(s))
    half = 0.5 * s[i0]
    lo = i0
    while lo > 0 and s[lo] >= half:
        lo -= 1
    hi = i0
    while hi < len(s) - 1 and s[hi] >= half:
        hi += 1
    # linear interpolation of the two half-maximum crossings (grid step 2 ueV)
    x_lo = grid[lo] + (half - s[lo]) / (s[lo + 1] - s[lo]) * (grid[lo + 1] - grid[lo])
    x_hi = grid[hi - 1] + (s[hi - 1] - half) / (s[hi - 1] - s[hi]) * (grid[hi] - grid[hi - 1])
    return float(x_hi - x_lo)


class IbmLineTables:
    """Per-temperature tables of FWHM_IBM(gamma_zpl) on a log grid (inverted by interpolation) and, per filter row, of
    the IBM transmissions t(gamma_zpl) of the X and XX lines.  Splines in (log gz, log t): the least-squares loop is
    then cheap; the final solution is re-evaluated EXACTLY by `exact_eps` (the spline error is reported, never hidden)."""

    def __init__(self, pp, rows, delta_xx):
        from .qd_gf import ibm_transmission
        self.pp, self.delta_xx = pp, float(delta_xx)
        self.grid = np.geomspace(_IBM_GRID_LO, _IBM_GRID_HI, _IBM_GRID_N)
        lg = np.log(self.grid)
        self._fw, self._tx, self._txx = {}, [], []
        for T in sorted({float(r["T"]) for r in rows}):
            fw = np.array([ibm_fwhm(pp, T, g) for g in self.grid])
            self._fw[T] = np.maximum.accumulate(fw)   # monotone for the inversion
        for r in rows:
            T, w, dx = float(r["T"]), float(r["w"]), float(r["dx"])
            tx = np.array([float(ibm_transmission(dx, pp, T, g, w_meV=w)) for g in self.grid])
            txx = np.array([float(ibm_transmission(dx - self.delta_xx, pp, T, g, w_meV=w)) for g in self.grid])
            self._tx.append(CubicSpline(lg, np.log(tx)))
            self._txx.append(CubicSpline(lg, np.log(txx)))

    def gamma_zpl(self, T, gamma_total):
        """ZPL width with FWHM_IBM = gamma_total, clipped to [0.05, 1] x gamma_total (spec bracket); also returns
        whether the root was bracketed (False -> clamped at 0.05 Gamma, i.e. the sidebands alone exceed Gamma)."""
        fw = self._fw[float(T)]
        if gamma_total <= fw[0]:
            return 0.05 * gamma_total, False
        if gamma_total >= fw[-1]:
            return gamma_total, True
        gz = float(np.exp(np.interp(gamma_total, fw, np.log(self.grid))))
        lo = 0.05 * gamma_total
        return (lo, False) if gz < lo else (min(gz, gamma_total), True)

    def eps(self, i, gz, r_xx):
        """eps = t_XX/t_X of filter row i at ZPL width gz (XX ZPL = r_xx gz); table-interpolated."""
        lgz = np.log(gz)
        lgx = np.log(max(r_xx * gz, _IBM_GRID_LO))
        return float(np.exp(self._txx[i](lgx) - self._tx[i](lgz))),             float(np.exp(self._tx[i](lgz)))


@lru_cache(maxsize=256)
def ibm_zpl_width(pp, T, gamma_total):
    """ZPL FWHM gz [meV] with FWHM_IBM(gz, T) = gamma_total (brentq on [0.05, 1] x gamma_total, clamped to the
    bracket end when the sidebands alone already exceed / the full ZPL stays below gamma_total).  Cached on
    (PhononParams, T, Gamma)."""
    from scipy.optimize import brentq
    lo, hi = 0.05 * gamma_total, gamma_total
    f = lambda g: ibm_fwhm(pp, T, g) - gamma_total
    if f(lo) > 0:
        return float(lo)
    if f(hi) < 0:
        return float(hi)
    return float(brentq(f, lo, hi, xtol=1e-6))


def ibm_exact_eps(pp, T, w, dx, delta_xx, gamma_total, r_xx):
    """Direct (non-tabulated) IBM eps at the narrowed ZPL width (XX ZPL = r_xx gz); returns (eps, t_x, gamma_zpl)."""
    from .qd_gf import ibm_transmission
    gz = ibm_zpl_width(pp, float(T), float(gamma_total))
    t_x = float(ibm_transmission(dx, pp, T, gz, w_meV=w))
    t_xx = float(ibm_transmission(dx - delta_xx, pp, T, r_xx * gz, w_meV=w))
    return t_xx / t_x, t_x, float(gz)


def fit_phase0(card: Card, dataset="g2_vs_T", seed=0, n_starts=40,
               lineshape="lorentzian", phonon=None) -> FitResult:
    """lineshape "lorentzian" (default, the legacy production fit, bit-identical) | "ibm" (opt-in IBM-calibrated
    refit, see the IBM block above; `phonon` = PhononParams kwargs, default the chatzarakis2023 card block)."""
    if lineshape not in ("lorentzian", "ibm"):
        raise ValueError(f"lineshape must be 'lorentzian' or 'ibm', got {lineshape!r}")
    ibm = lineshape == "ibm"
    # free set = every candidate the card declares as a range (swept, never averaged)
    CAND = BASE_FREE + ("r_xx", "E_lo")
    FREE = [k for k in CAND if k in card.params and card[k].is_free]

    def pack(values):
        return np.array([np.log10(values[k]) if k in LOG_FREE else values[k] for k in FREE])

    def unpack(x):
        return {k: (10.0**v if k in LOG_FREE else v) for k, v in zip(FREE, x)}

    ds = card.datasets[dataset]
    rows = sorted((dict(r) for r in ds.rows), key=lambda r: r["T"])
    for r in rows:
        r.setdefault("err", V_A_TOL)
        r.setdefault("bound", "value")

    gam_rows = sorted(card.datasets["gamma_vs_T"].rows, key=lambda r: r["T"]) \
        if "gamma_vs_T" in card.datasets else []
    tau_rows = sorted(card.datasets["tau_vs_T"].rows, key=lambda r: r["T"]) \
        if "tau_vs_T" in card.datasets else []

    fixed = {k: card[k].fixed for k in ("delta_xx", "E_lo", "E_b")
             if k in card.params and not card[k].is_free}
    if "mu_op" in card.params:
        fixed["mu"] = card["mu_op"].fixed  # F1b operating point (Module D)
    T_anchor = card["gamma_anchor_T"].fixed
    g_anchor = card["gamma_anchor"].fixed
    g_anchor_tol = card["gamma_anchor_tol"].fixed
    Ea_prior, Ea_tol = card["E_a_prior"].fixed, card["E_a_prior_tol"].fixed
    bp_prior, bp_tol = card["b_p_prior_log10"].fixed, card["b_p_prior_log10_tol"].fixed

    if ibm:
        from .qd_gf import PhononParams
        pp = PhononParams(**(phonon if phonon is not None else chatz2023_phonon()))
        tables = IbmLineTables(pp, rows, fixed["delta_xx"])

    def model_g2_all(th):
        if ibm:
            # same blocks as integrator.g2_of_T (retention, background, F1b loading); only eps comes from the IBM line
            out = []
            for i, r in enumerate(rows):
                gam = float(gamma_of_T(r["T"], th["gamma0"], th["a_ac"], th["b_lo"],
                                       th.get("E_lo", fixed.get("E_lo"))))
                gz, _ = tables.gamma_zpl(r["T"], gam)
                eps_i, _ = tables.eps(i, gz, th["r_xx"])
                S = float(retention(r["T"], th["a_esc"], th["E_a"], th["b_p"], fixed["E_b"]))
                rho = S / (S + th["b0"] + th["beta"] * (1.0 - S))
                out.append(g2_from(float(f1b_g2(fixed["mu"], eps_i)) if fixed.get("mu") else eps_i, rho))
            return np.array(out)
        return np.array([
            g2_of_T(r["T"], {**fixed, **th, "w": r["w"], "dx": r["dx"]}).g2
            for r in rows
        ])

    def residuals(x):
        th = unpack(x)
        out = []
        m = model_g2_all(th)
        out.extend(
            max(0.0, mi - r["g2"]) / (r["err"] / UPPER_BOUND_ERR_DIVISOR if ibm else r["err"])
            if r["bound"] == "upper"
            else (mi - r["g2"]) / r["err"]
            for mi, r in zip(m, rows))
        gm = lambda T: gamma_of_T(T, th["gamma0"], th["a_ac"], th["b_lo"],
                                  th.get("E_lo", fixed.get("E_lo")))
        for r in gam_rows:
            out.append((gm(r["T"]) - r["gamma_x"]) / (GAMMA_RELERR * r["gamma_x"]))
            if "gamma_xx" in r and "r_xx" in th:
                out.append((th["r_xx"] * gm(r["T"]) - r["gamma_xx"])
                           / (GAMMA_RELERR * r["gamma_xx"]))
        if tau_rows:
            S = lambda T: retention(T, th["a_esc"], th["E_a"], th["b_p"], fixed["E_b"])
            T0, tau0 = tau_rows[0]["T"], tau_rows[0]["tau_ns"]
            S0 = S(T0)
            for r in tau_rows[1:]:
                ratio_d = r["tau_ns"] / tau0
                out.append((S(r["T"]) / S0 - ratio_d) / (TAU_RELERR * ratio_d))
        if not gam_rows:  # anchor only when no real linewidth data
            out.append((gm(T_anchor) - g_anchor) / g_anchor_tol)
        out.append((th["E_a"] - Ea_prior) / Ea_tol)
        out.append((np.log10(th["b_p"]) - bp_prior) / bp_tol)
        return np.array(out)

    lo = pack({k: card[k].bounds[0] for k in FREE})
    hi = pack({k: card[k].bounds[1] for k in FREE})

    rng = np.random.default_rng(seed)
    best = None
    for i in range(n_starts):
        x0 = lo + (hi - lo) * (0.5 if i == 0 else rng.random(len(FREE)))
        try:
            sol = least_squares(residuals, x0, bounds=(lo, hi), method="trf")
        except Exception:
            continue
        if best is None or sol.cost < best.cost:
            best = sol

    th = unpack(best.x)
    p_full = {**fixed, **th}
    m = model_g2_all(th)
    ibm_aux = {}
    if ibm:
        # re-evaluate the solution EXACTLY (no spline): this is what the bundle and the ledger use
        m_tab, m, gzs, eps_x, brk = m, [], [], [], []
        for r in rows:
            gam = float(gamma_of_T(r["T"], th["gamma0"], th["a_ac"], th["b_lo"], th.get("E_lo", fixed.get("E_lo"))))
            e, _, gz = ibm_exact_eps(pp, r["T"], r["w"], r["dx"], fixed["delta_xx"], gam, th["r_xx"])
            S = float(retention(r["T"], th["a_esc"], th["E_a"], th["b_p"], fixed["E_b"]))
            rho = S / (S + th["b0"] + th["beta"] * (1.0 - S))
            m.append(g2_from(float(f1b_g2(fixed["mu"], e)) if fixed.get("mu") else e, rho))
            gzs.append(gz)
            eps_x.append(e)
            brk.append(bool(ibm_fwhm(pp, r["T"], gz) <= gam * 1.02 and gz > 0.05 * gam * 1.001))
        m = np.array(m)
        ibm_aux = {"gamma_zpl_meV": [float(x) for x in gzs], "eps": [float(x) for x in eps_x],
                   "T": [float(r["T"]) for r in rows], "fwhm_root_bracketed": brk,
                   "spline_vs_exact_max_abs_dg2": float(np.max(np.abs(m_tab - m))),
                   "upper_bound_err_divisor": UPPER_BOUND_ERR_DIVISOR,
                   "phonon": dict(phonon if phonon is not None else chatz2023_phonon())}
    resid = np.array([
        max(0.0, mi - r["g2"]) if r["bound"] == "upper" else mi - r["g2"]
        for mi, r in zip(m, rows)
    ])
    E_lo_fit = th.get("E_lo", fixed.get("E_lo"))
    g_at_anchor = float(gamma_of_T(T_anchor, th["gamma0"], th["a_ac"], th["b_lo"], E_lo_fit))

    aux = {}
    if ibm:
        aux["ibm"] = ibm_aux
    if gam_rows:
        gr = [abs(gamma_of_T(r["T"], th["gamma0"], th["a_ac"], th["b_lo"], E_lo_fit)
                  - r["gamma_x"]) / r["gamma_x"] for r in gam_rows]
        aux["gamma_worst_relerr"] = float(max(gr))
    if tau_rows:
        S = lambda T: retention(T, th["a_esc"], th["E_a"], th["b_p"], fixed["E_b"])
        S0, tau0 = S(tau_rows[0]["T"]), tau_rows[0]["tau_ns"]
        tr = [abs(S(r["T"]) / S0 - r["tau_ns"] / tau0) / (r["tau_ns"] / tau0)
              for r in tau_rows[1:]]
        aux["tau_worst_relerr"] = float(max(tr))

    tag = widest(
        ds.tag,
        *(card.datasets[d].tag for d in ("gamma_vs_T", "tau_vs_T") if d in card.datasets),
        *(card[k].tag for k in FREE),
        *(card[k].tag for k in ("delta_xx", "E_lo", "E_b") if k in card.params),
    )

    notes = []
    ph = card.placeholders()
    if ph:
        notes.append(
            "PLACEHOLDER inputs present (%s): fit exercises the machinery only -- "
            "the V-a gate CANNOT close on this run." % ", ".join(ph)
        )

    return FitResult(
        params=p_full, data=rows,
        model_g2=m, residuals=resid,
        passed=bool(np.all(np.abs(resid) <= V_A_TOL)),
        gamma_anchor=(g_at_anchor, g_anchor, g_anchor_tol, T_anchor),
        tag=tag, cost=float(best.cost), notes=notes, aux=aux,
    )
