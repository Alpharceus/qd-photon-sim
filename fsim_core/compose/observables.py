"""Observables chain (design 2.6, physics brief section 5) and tag propagation.

Brightness accounting (brief 5), every factor is reported:

    B_lens  = n_ph x eta_pol x beta x eta_out x T_NA x Z_filt          (first lens, "B_first_lens")
    B_fibre = B_lens x eta_fibre
    counts  = f_rep x B_fibre x eta_optics x eta_det                   (eta_pol appears once, in B_lens)

Cavity filtering of the phonon sideband (review finding 4) [E: Iles-Smith et al., Nat. Photonics 11, 521
(2017); Grange et al., PRL 114, 193601 (2015)]. A structure with a finite kappa transmits only the part of the
emission spectrum inside its Lorentzian T(E) = (kappa/2)^2/((E - E_mode)^2 + (kappa/2)^2). With the unit-area
IBM spectrum S = Z lorentzian(ZPL) + S_sb (qd_gf), the ZPL photons that reach the cavity channel are already
counted by the Purcell overlap L inside beta, so only the sideband needs filtering:

    Z_eff = (Z + int S_sb(E) T(E) dE)/(Z + int S_sb)       Z <= Z_eff <= 1
    q     = Z/Z_eff                                         ZPL fraction of the collected light
    efficiency x Z_eff,   hom_I = I x q^2                   (I_ZPL x (ZPL fraction)^2)

limits: kappa -> inf (T = 1) gives Z_eff = 1, q = Z: the unfiltered result (efficiency x 1, I x Z^2);
kappa -> 0 gives Z_eff = Z, q = 1: the ZPL-filter result (efficiency x Z, I x 1). Broadband structures
(kappa None) keep the unfiltered behaviour; setup.zpl_filter forces Z_eff = Z, q = 1. The sideband is
assumed to see the same cavity Lorentzian as the ZPL (no Purcell re-weighting of the sidebands) and the
background-channel photons are filtered the same way (conservative) [A].

n_ph is the solver's photons per pulse into the radiative channel. It replaces the product
eta_prep x eta_rad of the brief because it equals that product without re-excitation and also contains
re-excitation (both factors are reported separately: `eta_prep`, `eta_rad`, `eta_prep_x_eta_rad`).
eta_prep (FX2c) is the preparation PROBABILITY, as papers define it: P(n >= 1) with n the excitations a pulse
produces from the ground state (photons emitted inside the window + X, XX (TPE) or cavity photon (jc) still
present at its end), computed exactly as 1 - p0 from the no-jump evolution (solvers._window_p0), so 0 <= eta_prep
<= 1. The expected COUNT (FX2b's definition, P_X(end) + emitted, which exceeds 1 under re-excitation) is the
observable `n_exc`; the old end-of-window population is `P_X_end` (informational). Incoherent drives report
p_loaded as eta_prep and no n_exc.
XX channel (TPE cascade only): tau_XX_ps, g2_XX, n_XX, hom_XX (hom_XX carries the same q^2 as hom_I).
T_NA = 1: the structure's response() already received the setup NA, so eta_out contains the NA
collection (T_NA is reported as 1 for the record). Z_filt = Z(T) when the setup filters to the ZPL,
Z_eff (below) for a finite-kappa structure, else 1. Phonon bookkeeping (brief trap 3, tier P0): HOM factor
Z^2 UNFILTERED, efficiency x Z FILTERED, never both; a cavity interpolates between the two.

g2: `g2_0` is intrinsic. `g2_raw` adds uncorrelated Poissonian background of rate B = f_leak S + dark
(S = detected signal rate): g2_raw = 1 - rho^2 (1 - g2), rho = S/(S + B) -> g2 + 2 f_leak for small leak.
Laser leak is counted once: the drive's (plan.inputs["laser_leak"]) when nonzero, else the setup's.

Dark counts for a PULSED drive (devlib-fx2a (a)): the coincidence histogram integrates only a window
`setup.coinc_window_ns` around each peak, so a detector's dark counts reach the central peak only inside that
window. Per pulse and detector, with S_pp = S/f_rep the detected signal probability:

    B_pp = f_leak S_pp + dark_cps x coinc_window_ns x 1e-9          (window capped at the pulse period)
    rho  = S_pp/(S_pp + B_pp),   g2_raw = 1 - rho^2 (1 - g2)         (unchanged functional form)

equivalently a background RATE of f_leak S + dark_cps x min(window x f_rep x 1e-9, 1). Continuous drives keep
B = f_leak S + dark_cps. Check (Muller-like: S 2750/s at 76 MHz, dark 250/s, window 2 ns): the dark share is
0.152 x 250 = 38/s, rho = 0.9863 and g2_raw - g2 = 0.027 (the unwindowed treatment gave 0.16) [DR].

Blinking (devlib-fx2a (c)): `material.eta_blink` (on-fraction, default 1 [A]) multiplies the collected
brightness (B_lens, B_fibre, rate_detected_Hz) and is reported in `eta_rad_eff = eta_blink x eta_rad`; it
never changes the lifetime (blinking removes emission events, not the radiative rate).

HOM (FX2c, replaces the FX1 relation V_raw = ((1 - g2) I - g2)(1 - eps) that double counted g2): `hom_I` = the
solver's Grange overlap I = Int Int |G1(t,t')|^2 dt dt' / N^2, N = the TOTAL mean photons per pulse, x the phonon
factor. G1 is the first-order coherence of the whole photon stream, so I ALREADY contains the multi-photon
(re-excitation) terms: for an ideal emitter that re-emits twice it is far below the single-photon wavepacket
overlap M_s (closed case, verify_compose_engine: two instantaneous pi pulses, n = 0 or 2 photons, I = 1/2 while
every photon is a perfect wavepacket, M_s = 1). The exact zero-delay result for two independent identical sources
on a 50/50 splitter, normalised to the uncorrelated (far) peak as in the experiments, is [DR]
    C_par(0) = (N^2/4)(2 g2 + 2 - 2 I),   C_far = N^2,   V_raw = 1 - 2 C_par/C_far = I - g2
(second-order operator identity <c+ d+ d c> with c,d = (a +- b)/sqrt 2; Ollivier et al., PRL 126, 063602 (2021)
treat the same combination with M_s, V = (1 - g2) M_s - g2, which is the same physics only when I is the
single-photon M_s, so using that form with the solver's I subtracted g2 twice). An uncorrelated background
scales I by rho^2 = (1 - g2_raw)/(1 - g2) (N -> N/rho, Int Int |G1|^2 unchanged) and g2 -> g2_raw:
    hom_V_raw = (rho^2 I - g2_raw)(1 - eps_int),   hom_V_corr = (V_raw + g2_raw)/(1 - g2_raw) = I/(1 - g2) at
eps_int = 0 (the corrected value is the estimate of the single-photon M_s an experiment would quote: Tomm 0.925,
g2 0.021 -> 0.966). eps_int from a setup param "eps_int" if present else 0 [A]. Compare raw to raw and corrected
to corrected, never across.

Gated / binned / IRF observables (devlib-p10e, electrical drives). All [DR]; the closed forms are in the functions below.
  * `tau_EL_ps` = 1e3/(Gamma + k_nr + k_off): the electroluminescence decay of a gated diode. It is REPORTED next to
    tau_ps (= the radiative-lifetime proxy 1e3/(Gamma + k_nr)) and never fed back into Gamma0 / F_P (brief 3, trap 4).
  * `g2_binned` = g2_0 (1 - e^{-gamma T/2}) + e^{-gamma T/2} (pulsed, histogram bin = the pulse period T,
    gamma = the EL decay rate), with the Gaussian detector IRF (drive input irf_ps, or setup.irf_ps if present) leaking
    more of every peak into the neighbours (`binned_g2`); the peak-area g2 is `g2_0`. (Brief 2.3, trap 6.)
  * `hom_I_gated` = the time-filtered HOM overlap of the exciton photon inside a detection gate [t_g, t_g + T]
    (`gated_hom`), optionally averaged over a Gaussian relative emission-time delay (`gated_hom_jitter`); times the
    same ZPL-fraction q^2 as hom_I. Limits: T -> 0 gives 1, T -> inf gives (Gamma + k)/(Gamma + k + gamma*).
  * `g2_irf` (cw): the intrinsic dip g2_0 seen through a Gaussian IRF of FWHM irf_ps, 1 - (1 - g2_0) F with F the
    cw_g2 single-exponential dip factor at the dip recovery time 1/(Gamma + k + r); equals g2_0 when irf_ps = 0.

Spectral leakage (devlib-p10e2, brief 4.2-4.3, nitride high-T): `spectral_overlap` = eps = t_XX/t_X, the XX-line
transmission relative to X through the top-hat window (setup.filter_width_meV, offset filter_offset_meV; both lines
Lorentzian, widths hbar(Gamma + k + gamma*) and xx_width_ratio times that; XX offset E_b = -delta_xx_meV), reusing
`spectral.epsilon`. Conventions: rho = 1/(1 + (I_XX/I_X) eps) (setup.xx_x_ratio) and g2 = 1 - rho^2 (the uncorrelated
Tamariz Eq. 4 form: valid when the XX photon is not correlated with the following X inside the coherence time).

Tags. Each observable lists `deps`: tokens that are either a stage name (material | structure | drive |
setup = the leaves that stage actually READ during the run, tracked by the engine) or an fnmatch pattern
on leaf paths. Tag = widest(DR, tags of the matched leaves); DR is the floor because the model itself is
derived. T_K is an operating condition and does not enter tags [A].
"""
from __future__ import annotations

import fnmatch
import math
from collections.abc import Mapping
from dataclasses import dataclass
from typing import Any

from .. import qd_gf
from .params import widest_tag
from .registry import Observable, register

NAN = float("nan")
STAGES = ("material", "structure", "drive", "setup")
BRIGHTNESS_CONVENTION = ("B_lens = n_ph * eta_pol * beta * eta_out * T_NA * Z_filt (first lens; T_NA is "
                         "inside eta_out); B_fibre = B_lens * eta_fibre; detected = f_rep * B_fibre * "
                         "eta_optics * eta_det. n_ph = photons per pulse from the solver "
                         "(= eta_prep * eta_rad without re-excitation).")


def _read(c, stage, name, default):
    """Untracked-dict value of an optional leaf that is not part of a stage's ParamSpec-driven chain/build
    object; the read is recorded for provenance (observables also list it as a pattern dep)."""
    vals = c.vals.get(stage, {})
    if name in vals:
        c.reads.setdefault(stage, set()).add(name)
        return vals[name]
    return default


def pulsed_dark_rate(dark_cps, coinc_window_ns, f_rep_Hz):
    """Dark-count background RATE that reaches the central coincidence peak for a pulsed drive:
    dark_cps x min(window x f_rep, 1) (the window cannot exceed the pulse period)."""
    return float(dark_cps) * min(float(coinc_window_ns) * 1e-9 * float(f_rep_Hz), 1.0)


def _int_exp(c, lo, hi, shift=0.0):
    """int_lo^hi e^{c t + shift} dt, overflow-safe (anchored at the larger endpoint exponent) and stable for c -> 0."""
    w = hi - lo
    if abs(c * w) < 1e-12:
        return w * math.exp(c * lo + shift)
    if c > 0.0:
        return math.exp(c * hi + shift) * (-math.expm1(-c * w)) / c
    return math.exp(c * lo + shift) * math.expm1(c * w) / c


def _gate_F(Gamma_d, gamma_star, lo, hi):
    """F(lo, hi) = int_lo^hi int_lo^hi e^{-Gamma_d (t + t')} e^{-gamma_star |t - t'|} dt dt' (closed form)."""
    if not hi > lo:
        return 0.0
    a = Gamma_d + gamma_star
    return (2.0 / a) * (_int_exp(-2.0 * Gamma_d, lo, hi) - _int_exp(gamma_star - Gamma_d, lo, hi, shift=-a * hi))


def gated_hom(Gamma, gamma_star, t_gate_ns, T_ns, k_loss=0.0, delay_ns=0.0):
    """Time-filtered HOM overlap of the photon of one excitation inside the detection gate [t_g, t_g + T] [DR].
    Repo convention: population decay Gamma_d = Gamma + k_loss, coherence decay (Gamma_d + gamma*)/2, first-order
    coherence G1(t, t') = Gamma e^{-Gamma_d (t + t')/2} e^{-gamma* |t' - t|/2} (t, t' >= 0). Two photons from
    independent identical emitters, the second emitted `delay_ns` later, detected in the same gate:
        I_T = int int_gate^2 G_A(t, t') G_B(t', t) / (n_A n_B),   n = int_gate Gamma e^{-Gamma_d (t - t_emit)} dt,
    which is the Grange overlap I = 2 S / N^2 restricted to the gate (delay 0). Limits: T -> 0 gives 1;
    T -> inf, t_g = 0 gives Gamma_d/(Gamma_d + gamma*) (= T2/(2 T1) for k = 0); the delay enters as
    e^{Gamma_d delay} before the common start max(t_g, 0, delay) (Gamma_d |delay| for the unwindowed packet:
    I(delay) = Gamma_d/(Gamma_d + gamma*) e^{-Gamma_d |delay|}). Returns dict(n_T, I_T): n_T = the fraction of the
    photons that falls in the gate (per excitation, relative to the radiative yield Gamma/Gamma_d)."""
    Gamma, gamma_star, t_g, T = float(Gamma), float(gamma_star), float(t_gate_ns), float(T_ns)
    Gd = Gamma + float(k_loss)
    d = float(delay_ns)
    E = t_g + T
    if not (T > 0.0 and Gd > 0.0):
        return {"n_T": 0.0, "I_T": NAN}
    lo = max(t_g, 0.0, d)
    n_A = Gamma * _int_exp(-Gd, max(t_g, 0.0), E)
    n_B = Gamma * math.exp(Gd * d) * _int_exp(-Gd, max(t_g, d), E) if E > d else 0.0
    if not (n_A > 0.0 and n_B > 0.0):
        return {"n_T": n_A / (Gamma / Gd) if n_A > 0.0 else 0.0, "I_T": NAN, "num": 0.0, "nn": 0.0}
    num = Gamma ** 2 * math.exp(Gd * d) * _gate_F(Gd, gamma_star, lo, E)
    return {"n_T": float(n_A / (Gamma / Gd)), "I_T": float(num / (n_A * n_B)), "num": float(num),
            "nn": float(n_A * n_B)}


def gated_hom_jitter(Gamma, gamma_star, t_gate_ns, T_ns, sigma_ns, k_loss=0.0, n_nodes=161):
    """Gate-filtered HOM overlap averaged over a Gaussian RELATIVE emission-time delay of std `sigma_ns`
    (Bennett 2008: the Stark sweep jitters the emission time) [DR/A]. The average is over COINCIDENCES,
        I = int P(d) N(d) dd / int P(d) n_A n_B(d) dd,   N = the gate-restricted overlap numerator of gated_hom,
    i.e. each delay is weighted by its detected-pair probability (the brief's unweighted int P I_T lets the vanishing
    detection probability of a late partner inflate I_T -> 1 near d -> T). Trapezoid on +/- 6 sigma; sigma = 0
    returns gated_hom. The jitter is the one [A] input of the Bennett comparison; for two independent sources with
    per-photon jitter s use sigma = sqrt(2) s."""
    import numpy as np
    if not sigma_ns > 0.0:
        return gated_hom(Gamma, gamma_star, t_gate_ns, T_ns, k_loss)["I_T"]
    x = np.linspace(-6.0, 6.0, int(n_nodes))
    w = np.exp(-0.5 * x * x)
    g = [gated_hom(Gamma, gamma_star, t_gate_ns, T_ns, k_loss, delay_ns=float(sigma_ns * xi)) for xi in x]
    num = np.array([v["num"] for v in g])
    nn = np.array([v["nn"] for v in g])
    den = float(np.sum(w * nn))
    return float(np.sum(w * num) / den) if den > 0.0 else NAN


def binned_g2(g2_0, gamma_ns, T_bin_ns, irf_ps=0.0):
    """Central-bin g2 of a pulsed histogram whose bin equals the pulse period T (GHz drives) [DR, brief 2.3]:
        g2_bin = 1 - (1 - g2_0) c0,   c0 = fraction of one coincidence peak inside its own bin,
    two-sided exponential peaks exp(-gamma|tau|) of rate gamma (the EL decay) give c0 = 1 - e^{-gamma T/2}, i.e.
    g2_bin = g2_0 (1 - e^{-gamma T/2}) + e^{-gamma T/2} (IRF-free closed form). A Gaussian detector IRF of FWHM irf_ps
    (the cross-correlation of two detectors has std sqrt(2) sigma_IRF) broadens each peak and lowers c0 further
    (peak density = Laplace convolved with a normal, integrated numerically). The peak-area g2 is g2_0 itself."""
    gamma, T = float(gamma_ns), float(T_bin_ns)
    if not (gamma > 0.0 and T > 0.0):
        return NAN
    if not irf_ps > 0.0:
        c0 = -math.expm1(-gamma * T / 2.0)
    else:
        from scipy.integrate import quad
        from scipy.special import ndtr
        sig = math.sqrt(2.0) * float(irf_ps) * 1e-3 / (2.0 * math.sqrt(2.0 * math.log(2.0)))   # ns, cross-correlation std

        def f(x):                       # Laplace(gamma) * N(0, sig^2) density
            a = gamma * sig * sig
            e1 = math.exp(min(0.5 * gamma * gamma * sig * sig - gamma * x, 700.0)) * ndtr((x - a) / sig)
            e2 = math.exp(min(0.5 * gamma * gamma * sig * sig + gamma * x, 700.0)) * ndtr(-(x + a) / sig)
            return 0.5 * gamma * (e1 + e2)
        c0 = 2.0 * quad(f, 0.0, T / 2.0, epsabs=1e-13, epsrel=1e-12, limit=200)[0]      # f is even
    return float(1.0 - (1.0 - float(g2_0)) * c0)


def irf_dip_g2(g2_0, tau_dip_ns, irf_fwhm_ps):
    """g2(0) of a single-exponential dip 1 - (1 - g2_0) exp(-|tau|/tau_dip) seen through a Gaussian IRF of FWHM
    irf_fwhm_ps: 1 - (1 - g2_0) F (cw_g2._dip_factor, the Reischle 2008 procedure). irf = 0 returns g2_0."""
    if not irf_fwhm_ps > 0.0:
        return float(g2_0)
    from .. import cw_g2
    return float(1.0 - (1.0 - float(g2_0)) * cw_g2._dip_factor(float(tau_dip_ns), float(irf_fwhm_ps), "gaussian"))


def hom_v_raw(I, g2_raw, eps_int=0.0, g2_int=None):
    """Predicted raw visibility V_raw = (rho^2 I - g2_raw)(1 - eps_int) with I the solver's whole-stream overlap
    and rho^2 = (1 - g2_raw)/(1 - g2_int); g2_int=None means no background (rho = 1, V = I - g2). The pre-FX2c
    form ((1 - g2) I - g2) is the Ollivier relation for the single-photon M_s, not for this I (module doc)."""
    rho2 = 1.0 if (g2_int is None or not g2_int < 1.0) else (1.0 - g2_raw) / (1.0 - g2_int)
    return (rho2 * I - g2_raw) * (1.0 - eps_int)


def hom_v_corrected(V_raw, g2):
    """g2-corrected visibility V_corr = (V_raw + g2)/(1 - g2) (Tomm 0.925, 0.021 -> 0.966)."""
    return (V_raw + g2) / (1.0 - g2)


def cavity_sideband_filter(ph, T_K, kappa_meV, detuning_meV, gamma_zpl_meV):
    """(Z, Z_eff) of the cavity-filtered output (module doc): detuning_meV = E_ZPL - E_mode."""
    import numpy as np
    omega, s_sb, z = qd_gf._spectrum_parts(ph, float(T_K), float(gamma_zpl_meV))
    sb_area = float(np.trapezoid(s_sb, omega))
    xc = omega + float(detuning_meV)
    acc = (kappa_meV / 2.0) ** 2 / (xc ** 2 + (kappa_meV / 2.0) ** 2)
    t_sb = float(np.trapezoid(s_sb * acc, omega))
    return float(z), float((z + t_sb) / (z + sb_area))


def g2_with_background(g2, signal_rate, background_rate):
    s, b = signal_rate, background_rate
    if b == 0.0:
        return g2                                  # no background: g2_raw = g2 even at zero signal
    if s <= 0.0:
        return NAN
    rho = s / (s + b)
    return 1.0 - rho * rho * (1.0 - g2)


class TrackedVals(Mapping):
    """Read-tracking view of one stage's parameter values (any access, including `in`, `get` and
    iteration of items/values, counts as a read)."""

    def __init__(self, data: dict, reads: set):
        self._d = dict(data)
        self._r = reads

    def __getitem__(self, k):
        v = self._d[k]
        self._r.add(k)
        return v

    def __iter__(self):
        return iter(self._d)

    def __len__(self):
        return len(self._d)


class Context:
    """Everything an observable needs for one (T) point."""

    def __init__(self, em, sr, rt, plan, res, setup, T_K, vals, reads):
        self.em, self.sr, self.rt, self.plan, self.res, self.setup = em, sr, rt, plan, res, setup
        self.T_K = float(T_K)
        self.vals = vals                      # {stage: {param: value}} (untracked copies)
        self.reads = reads                    # {stage: set(param names)}
        self._d = None

    @property
    def d(self) -> dict:
        if self._d is None:
            self._d = _derive(self)
        return self._d


def _derive(c: Context) -> dict:
    em, sr, rt, plan, res, su = c.em, c.sr, c.rt, c.plan, c.res, c.setup
    d: dict[str, Any] = {}
    why: dict[str, str] = {}
    warnings: list = []
    Gamma, k = float(rt.Gamma_ns), float(rt.k_nr_ns)
    ph = em.phonon
    Z = float(qd_gf.zpl_weight(ph, c.T_K)) if ph is not None else 1.0
    kap = sr.kappa_meV
    if su.zpl_filter:
        Zf, q = Z, 1.0
    elif ph is not None and kap is not None and kap > 0.0 and sr.E_mode_meV is not None:
        _, Zf = cavity_sideband_filter(ph, c.T_K, float(kap), float(em.E_X_meV) - float(sr.E_mode_meV),
                                       lb_hbar_ueV() * 1e-3 * (Gamma + k + float(rt.gamma_star_ns)))
        q = Z / Zf
    else:
        Zf, q = 1.0, Z
    beta = float(res.extra.get("beta", rt.beta))
    eta_pol = float(su.eta_pol) * float(sr.pol_factor)
    n_ph = float(res.n_ph)
    pulsed = plan.timing == "pulsed"
    f_rep = float(plan.f_rep_MHz) * 1e6 if plan.f_rep_MHz else None
    eta_blink = float(_read(c, "material", "eta_blink", 1.0))
    d.update(eta_blink=eta_blink, eta_rad_eff=eta_blink * (Gamma / (Gamma + k)))
    d.update(Z_phonon=Z, Z_filt=Zf, beta=beta, eta_out=float(sr.eta_out), T_NA=1.0, eta_pol=eta_pol,
             eta_rad=Gamma / (Gamma + k), n_ph=n_ph, Gamma_ns=Gamma, gamma_star_ns=float(rt.gamma_star_ns),
             F_P_eff=float(rt.F_P_eff), F_P_peak=float(sr.F_P), lambda_nm=float(em.lambda_nm),
             linewidth_ueV=lb_hbar_ueV() * (Gamma + k + float(rt.gamma_star_ns)),
             detuning_ueV=float(rt.detuning_ueV), hom_zpl_fraction=q)
    d["tau_purcell_ps"] = 1e3 / (Gamma + k)
    # devlib-p10e: the EL decay of a gated diode (tau_ps / Gamma0 stay the radiative values, never refit from EL)
    k_off = float(res.extra.get("k_esc_off_ns", 0.0) or 0.0)
    d["tau_EL_ps"] = float(res.extra["tau_EL_ps"]) if "tau_EL_ps" in res.extra else 1e3 / (Gamma + k + k_off)
    d["tau_ps"] = float(res.extra["tau_jc_ps"]) if "tau_jc_ps" in res.extra else d["tau_purcell_ps"]
    d["tau_1e_ps"] = float(res.extra["tau_1e_ps"]) if "tau_1e_ps" in res.extra else d["tau_ps"]
    d["eta_extract"] = beta * d["eta_out"]
    bd = _read(c, "structure", "beta_direct", 0.0)
    if bd is not None and float(bd) > 0.0:
        d["beta_entered"] = float(bd)
    else:
        d["beta_entered"] = NAN
        why["beta_entered"] = "the structure has no entered beta (structure.beta_direct absent or 0)"
    d["inversion"] = NAN if res.inversion is None else float(res.inversion)
    irf_ps = float(plan.inputs.get("irf_ps", 0.0) or 0.0) or float(_read(c, "setup", "irf_ps", 0.0) or 0.0)
    d["irf_ps"] = irf_ps
    if "eta_prep" in res.extra:
        # FX2c: preparation probability P(n >= 1) = 1 - p0 (the expected count is n_exc)
        d["eta_prep"] = float(res.extra["eta_prep"])
    elif res.inversion is None and "p_loaded" in res.extra:
        d["eta_prep"] = float(res.extra["p_loaded"])
    else:
        d["eta_prep"] = d["inversion"]
    d["n_exc"] = float(res.extra["n_exc"]) if "n_exc" in res.extra else NAN
    if "n_exc" not in res.extra:
        why["n_exc"] = "only the lindblad / jc / polaron coherent solvers report the expected excitation count"
    d["P_X_end"] = float(res.extra["P_X_end"]) if "P_X_end" in res.extra else d["inversion"]
    if not math.isfinite(d["P_X_end"]):
        why["P_X_end"] = "the solver has no end-of-window population (incoherent / cw drive)"
    d["eta_prep_x_eta_rad"] = d["eta_prep"] * d["eta_rad"]
    if pulsed and math.isfinite(n_ph):
        B = n_ph * eta_blink * eta_pol * beta * d["eta_out"] * d["T_NA"] * Zf
        d["B_lens"] = B
        d["B_fibre"] = B * float(su.eta_fibre)
        S = f_rep * d["B_fibre"] * float(su.eta_optics) * float(su.eta_det) if f_rep else NAN
    else:
        d["B_lens"] = d["B_fibre"] = NAN
        why["B_lens"] = why["B_fibre"] = "continuous drive: brightness per pulse undefined (see rate_detected_Hz)"
        flux = res.extra.get("flux_ns")
        S = (float(flux) * 1e9 * eta_blink * eta_pol * beta * d["eta_out"] * Zf * float(su.eta_fibre)
             * float(su.eta_optics) * float(su.eta_det)) if flux is not None else NAN
    d["rate_detected_Hz"] = S
    if not math.isfinite(S):
        why["rate_detected_Hz"] = "no photon flux available from the solver / plan"
    leak = float(plan.inputs.get("laser_leak", 0.0) or 0.0) or float(su.laser_leak)
    g2 = float(res.g2_0)
    dark = float(su.dark_cps)
    if pulsed and f_rep:
        win = float(_read(c, "setup", "coinc_window_ns", 2.0))
        dark = pulsed_dark_rate(dark, win, f_rep)          # darks only inside the coincidence window
        d["coinc_window_ns"] = win
    d["dark_bg_Hz"] = dark
    bg = leak * S + dark if math.isfinite(S) else NAN
    d["laser_leak_used"] = leak
    d["g2_0"] = g2
    if leak == 0 and su.dark_cps == 0:
        d["g2_raw"] = g2                                       # no background: raw = intrinsic, whatever S is
    elif math.isfinite(bg):
        d["g2_raw"] = g2_with_background(g2, S, bg)
        if not math.isfinite(d["g2_raw"]):
            why["g2_raw"] = "no detected signal (S <= 0) with a non-zero background: g2_raw undefined"
    else:
        d["g2_raw"] = NAN
        why["g2_raw"] = "background rate needs the detected signal rate S, which is unavailable"
    # binned (GHz) histogram and IRF-convolved dip
    if pulsed and f_rep:
        d["g2_binned"] = binned_g2(g2, Gamma + k + k_off, 1e9 / f_rep, irf_ps)
    else:
        d["g2_binned"] = NAN
        why["g2_binned"] = "binned histogram g2 needs a pulsed drive (bin = pulse period)"
    if not pulsed and math.isfinite(float(res.extra.get("tau_dip_ns", NAN))):
        d["g2_irf"] = irf_dip_g2(g2, float(res.extra["tau_dip_ns"]), irf_ps)
    else:
        d["g2_irf"] = NAN
        why["g2_irf"] = "IRF-convolved dip needs a cw rate-solver drive (dip recovery time)"
    I1 = float(res.I_single)
    d["hom_I"] = I1 * q * q if math.isfinite(I1) else NAN
    hg = plan.inputs.get("hom_gate")
    if hg:
        sg = float(hg.get("jitter_ns", 0.0))
        gs = float(rt.gamma_star_ns)
        Ig = gated_hom_jitter(Gamma, gs, float(hg["start_ns"]), float(hg["width_ns"]), sg, k + k_off)
        d["hom_I_gated"] = Ig * q * q if math.isfinite(Ig) else NAN
        d["hom_gate_fraction"] = gated_hom(Gamma, gs, float(hg["start_ns"]), float(hg["width_ns"]), k + k_off)["n_T"]
    else:
        d["hom_I_gated"] = NAN
        why["hom_I_gated"] = "no HOM detection gate set (drive.hom_gate_ps = 0)"
    # spectral X/XX leakage through a top-hat window (nitride high-T bookkeeping, brief 4.2)
    d.update(spectral_overlap_terms(c, em))
    # FX2b: XX-line stream of a TPE cascade (cascade solver only)
    for key in ("tau_XX_ps", "g2_XX", "hom_XX", "n_XX"):
        d[key] = float(res.extra[key]) if key in res.extra else NAN
        if key not in res.extra:
            why[key] = "no XX channel: only the TPE (cascade) solver reports the biexciton photon stream"
    if math.isfinite(d["hom_XX"]):
        d["hom_XX"] *= q * q                                   # same ZPL-fraction bookkeeping as hom_I [A]
    if not math.isfinite(d["hom_I"]):
        why["hom_I"] = why["hom_V"] = "solver has no coherence information (rate solver / cw drive)"
    eps = float(_read(c, "setup", "eps_int", 0.0))
    if math.isfinite(d["hom_I"]) and math.isfinite(d["g2_raw"]):
        d["hom_V_raw"] = hom_v_raw(d["hom_I"], d["g2_raw"], eps, g2_int=g2)
        d["hom_V_corr"] = hom_v_corrected(d["hom_V_raw"], d["g2_raw"])
    else:
        d["hom_V_raw"] = d["hom_V_corr"] = NAN
    if bool(plan.inputs.get("cross_pol")) and eta_pol >= 1.0:
        warnings.append("cross-polarised collection with eta_pol = 1: an unpolarised (trion / unpolarised X) emitter loses "
                        "1/2 in a cross-polarised setup unless the structure is polarisation-selective (brief trap 9); "
                        "set setup.eta_pol = 0.5 if that applies [A]")
    d["_warnings"] = warnings
    d["_why"] = why
    return d


def spectral_overlap_terms(c, em) -> dict:
    """eps = t_XX/t_X, rho = 1/(1 + ratio eps), g2 = 1 - rho^2 through the setup's top-hat window [DR]."""
    from .. import spectral
    w = float(_read(c, "setup", "filter_width_meV", 0.0) or 0.0)
    off = float(_read(c, "setup", "filter_offset_meV", 0.0) or 0.0)
    wr = float(_read(c, "setup", "xx_width_ratio", 1.5))
    ratio = float(_read(c, "setup", "xx_x_ratio", 0.0) or 0.0)
    gx = lb_hbar_ueV() * 1e-3 * (float(c.rt.Gamma_ns) + float(c.rt.k_nr_ns) + float(c.rt.gamma_star_ns))
    res = spectral.epsilon(-float(em.delta_xx_meV), gx, wr * gx, w=w if w > 0.0 else None, dx=off)
    eps = float(res.eps)
    rho = 1.0 / (1.0 + ratio * eps)
    return {"spectral_overlap": eps, "rho_spectral": rho, "g2_spectral": 1.0 - rho * rho,
            "t_X_window": float(res.t_x), "t_XX_window": float(res.t_xx)}


def lb_hbar_ueV():
    from .. import lindblad
    return lindblad.HBAR_MEV_NS * 1e3


@dataclass(frozen=True)
class _Obs(Observable):
    key: str = ""
    by_conv: Any = None                      # {convention: derived key}

    def compute(self, ctx, convention="intrinsic") -> float:
        if convention not in self.conventions:
            raise ValueError(f"observable {self.id!r} has no convention {convention!r}; "
                             f"known: {list(self.conventions)}")
        key = (self.by_conv or {}).get(convention, self.key or self.id)
        return float(ctx.d[key])

    def invalid_reason(self, ctx, convention="intrinsic") -> str:
        key = (self.by_conv or {}).get(convention, self.key or self.id)
        return ctx.d["_why"].get(key, ctx.d["_why"].get(self.id, ""))


def _o(id_, label, unit, deps, conventions=("intrinsic",), passthrough_of=None, key="", by_conv=None):
    return register(_Obs(id=id_, label=label, unit=unit, deps=tuple(deps), conventions=tuple(conventions),
                         passthrough_of=passthrough_of, key=key or id_, by_conv=by_conv))


_MS = ("material", "structure")
_ALL = ("material", "structure", "drive")
_CH = ("material", "structure", "drive", "setup", "material.eta_blink", "setup.coinc_window_ns", "setup.eps_int")
_o("lambda_nm", "Emission wavelength", "nm", ("material.lambda_nm", "material.E_X_meV"),
   passthrough_of="material.lambda_nm")
_o("F_P_eff", "Effective Purcell factor F_P L", "1", _MS)
_o("F_P_peak", "Cavity-mode Purcell factor at zero detuning (structure.F_P)", "1", _MS, passthrough_of="structure.F_P")
_o("beta", "Beta factor", "1", _MS, passthrough_of="structure.beta")
_o("beta_entered", "Beta entered directly (structure.beta_direct, e.g. phc_waveguide)", "1",
   ("structure.beta_direct",), passthrough_of="structure.beta_direct")
_o("eta_out", "Extraction into the collected mode", "1", ("structure",))
_o("eta_extract", "Extraction beta x eta_out", "1", _MS)
_o("T_NA", "NA transmission (inside eta_out)", "1", ("setup.NA",))
_o("eta_pol", "Polarisation factor", "1", ("setup.eta_pol", "structure"))
_o("eta_prep", "State-preparation probability P(n >= 1)", "1", _ALL)
_o("n_exc", "Expected excitation count per pulse (P_X(end) + photons emitted in the window; may exceed 1)", "1", _ALL)
_o("eta_rad", "Radiative efficiency Gamma/(Gamma+k)", "1", _MS)
_o("eta_rad_eff", "Effective radiative efficiency eta_blink x Gamma/(Gamma+k)", "1", _MS + ("material.eta_blink",))
_o("eta_blink", "Blinking on-fraction (brightness factor, not lifetime)", "1", ("material.eta_blink",),
   passthrough_of="material.eta_blink")
_o("eta_prep_x_eta_rad", "eta_prep x eta_rad", "1", _ALL)
_o("P_X_end", "Population at the end of the pulse window (informational; includes decay during the window)", "1", _ALL)
_o("tau_XX_ps", "XX lifetime 1e3/(gamma_XX + k_XX) read off the cascade Liouvillian (TPE)", "ps", _ALL)
_o("g2_XX", "g2(0) of the XX-line photon stream (TPE)", "1", _ALL)
_o("hom_XX", "HOM overlap of the XX-line photons, instantaneous preparation (TPE)", "1", _ALL + ("setup",))
_o("n_XX", "XX-line photons per pulse (TPE)", "1", _ALL)
_o("inversion", "Inversion after the pulse", "1", _ALL)
_o("Z_phonon", "ZPL weight Z(T)", "1", ("material",))
_o("Z_filt", "ZPL / cavity sideband filter factor", "1", ("material", "structure", "setup.zpl_filter"))
_o("hom_zpl_fraction", "ZPL fraction q of the collected light (hom_I = I q^2)", "1", ("material", "structure", "setup.zpl_filter"))
_o("gamma_star_ns", "Pure dephasing (repo convention)", "1/ns", ("material",))
_o("Gamma_ns", "Radiative rate", "1/ns", _MS)
_o("n_ph", "Photons per pulse", "1", _ALL)
_o("tau_ps", "Lifetime (slowest Liouvillian decay eigenvalue when solver jc, else Purcell)", "ps", _ALL)
_o("tau_1e_ps", "JC flux 1/e time after the peak (informational; = tau_ps for non-jc)", "ps", _ALL)
_o("tau_purcell_ps", "Purcell lifetime 1e3/(Gamma+k)", "ps", _MS)
_o("tau_EL_ps", "Electroluminescence decay time 1e3/(Gamma+k+k_off) of a gated diode (reported, never fed into Gamma0)", "ps",
   _ALL)
_o("g2_binned", "Central-bin g2 of the pulsed histogram (bin = pulse period; EL decay leaks neighbouring peaks, IRF adds more)",
   "1", _ALL + ("setup",))
_o("spectral_overlap", "X-XX spectral leakage eps = t_XX/t_X through the top-hat window (rho, g2 = 1 - rho^2 conventions)",
   "1", _MS + ("setup.filter_*", "setup.xx_*"), ("intrinsic", "rho", "g2"),
   key="spectral_overlap", by_conv={"intrinsic": "spectral_overlap", "rho": "rho_spectral", "g2": "g2_spectral"})
_o("g2_irf", "g2(0) of the cw dip seen through the detector IRF", "1", _ALL + ("setup",))
_o("hom_I_gated", "HOM overlap inside a detection gate (Stark-gated source), q^2 included", "1", _ALL + ("setup",))
_o("linewidth_ueV", "Emitter linewidth hbar(Gamma+k+gamma*)", "ueV", _MS)
_o("detuning_ueV", "Emitter minus cavity detuning", "ueV", _MS)
_o("g2_0", "g2(0)", "1", _ALL + ("setup.laser_leak", "setup.dark_cps"), ("intrinsic", "raw"),
   by_conv={"intrinsic": "g2_0", "raw": "g2_raw"})
_o("g2_raw", "g2(0) with laser leak and dark counts", "1", _CH, key="g2_raw")
_o("hom_I", "HOM overlap I", "1", _ALL + ("setup",), ("intrinsic", "raw", "corrected"),
   by_conv={"intrinsic": "hom_I", "raw": "hom_V_raw", "corrected": "hom_V_corr"})
_o("hom_V", "HOM visibility (raw|corrected)", "1", _CH, ("raw", "corrected"),
   by_conv={"raw": "hom_V_raw", "corrected": "hom_V_corr"})
_o("hom_V_raw", "HOM visibility, raw", "1", _CH, key="hom_V_raw")
_o("hom_V_corr", "HOM visibility, g2-corrected", "1", _CH, key="hom_V_corr")
_o("B_first_lens", "First-lens brightness", "1", _CH, key="B_lens")
_o("B_lens", "First-lens brightness (alias of B_first_lens)", "1", _CH, key="B_lens")
_o("B_fibre", "Brightness into the fibre", "1", _CH)
_o("rate_detected_Hz", "Detected count rate", "Hz", _CH)


def get_observable(oid):
    from .registry import get
    return get("observable", oid)


def parse_request(spec_str):
    """'id' or 'id:convention' -> (id, convention)."""
    if ":" in spec_str:
        a, b = spec_str.split(":", 1)
        return a, b
    return spec_str, "intrinsic"


def resolve_inputs(obs, leaves: dict, reads: dict) -> list:
    """Leaf paths an observable depends on (module doc)."""
    out: list = []

    def stage_paths(stage):
        return [f"{stage}.{n}" for n in sorted(reads.get(stage, ())) if f"{stage}.{n}" in leaves]

    fallback = None
    for tok in obs.deps:
        if tok in STAGES:
            out += stage_paths(tok)
        else:
            hit = [p for p in leaves if fnmatch.fnmatchcase(p, tok)]
            out += hit
            fallback = fallback or tok.split(".")[0]
    if not out and fallback in STAGES:
        out = stage_paths(fallback)
    seen = set()
    return [p for p in out if not (p in seen or seen.add(p))]


def tag_of(obs, leaves: dict, reads: dict, floor: str = ""):
    """(tag, inputs): widest of DR and the tags of the leaves the observable depends on. `floor` is the
    structure response's tag_floor (constants hidden inside response(), e.g. DBR penetration [E]); it applies
    to observables that depend on the structure."""
    inputs = resolve_inputs(obs, leaves, reads)
    tags = [leaves[p].tag for p in inputs]
    if floor and ("structure" in obs.deps or any(p.startswith("structure.") for p in inputs)):
        return widest_tag("DR", floor, *tags), inputs
    if inputs and not any(t in STAGES for t in obs.deps):      # pure passthrough of entered leaves
        return widest_tag(*tags), inputs
    return widest_tag("DR", *tags), inputs
