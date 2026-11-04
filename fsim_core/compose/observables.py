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
r4 red-team item 1: with setup.filter_fwhm_GHz > 0 (and no setup.zpl_filter, which keeps priority) the HOM ZPL fraction q comes from the
etalon's own Lorentzian window (window_filter_parts, width = the filter FWHM in meV, ZPL-centred [A]): q = Z T_zpl/Zf_window. The
brightness is taxed only when filter_applies_to = g2_hom_brightness, by the sideband loss alone (Z_filt = Z + (Zf_window - Z T_zpl)/T_zpl),
because the ZPL Lorentzian transmission is already in the cascade's filtered n_ph; a narrow filter therefore gives Z x T_zpl.

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
collection (T_NA is reported as 1 for the record). Z_filt = Z(T) when the setup filters to the ZPL, else 1; r4-cavity A2: for
a finite-kappa structure the CAVITY-collected light has Z_filt = 1 (the cavity is the filter and its sideband loss is already in
beta / the entered F_P) or Z/Z_eff with a ZPL filter, and only the leaky beta0 channel keeps Z / 1. Phonon bookkeeping (brief trap 3, tier P0): HOM factor
Z^2 UNFILTERED, efficiency x Z FILTERED, never both; a cavity interpolates between the two.

g2: `g2_0` is intrinsic (it already contains a drive-embedded `bg_frac`, plan.inputs["bg_frac"], where a drive has one:
above_band, p_shell, electrical). `g2_raw` adds uncorrelated Poissonian background of rate B = f_leak S + dark
(S = detected signal rate): g2_raw = 1 - rho^2 (1 - g2), rho = S/(S + B) -> g2 + 2 f_leak for small leak.
Laser leak is counted once: the drive's (plan.inputs["laser_leak"]) when nonzero, else the setup's.
`g2_0:bg_corrected` (devlib-fx3b) is the paper convention g_bgc = (g2_raw - (1 - rho^2))/rho^2 with the BACKGROUND rho =
(1 - bg_frac)/(1 + laser_leak) (drive-embedded fraction and laser leak; dark counts are not corrected, the "dark-count
corrected" value of a paper is a separate number). It is the exact inverse of the mixing law, so it returns the signal
g2 (Sartison, Lazic Eq. 5, Munsch g2*) and is what a background-corrected paper value is compared with.

Photon rates (devlib-fx3b): `B_first_lens_Hz` / `B_fibre_Hz` = the photon rate into the first lens / fibre, before
eta_optics and eta_det: pulsed f_rep x B_first_lens (B_fibre); cw (flux_ns = Gamma n_X,ss from the solver) x the same
chain beta eta_out eta_pol Z_filt (x eta_fibre); `rate_detected_Hz` is that times eta_optics eta_det.

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
  * `hom_I_gated` = the time-filtered HOM overlap of the exciton photon inside a detection gate on the emission-time
    axis (`gated_hom`): a top-hat [t_g, t_g + T] or (devlib-fx3b, gate_shape "gauss") a Gaussian window of intensity
    FWHM T centred at t_g (Bennett's 31 ps window is a Gaussian), optionally averaged over a Gaussian relative timing
    jitter of the two photons' gates (`gated_hom_jitter`: the displaced gates share only part of the photon, so the
    overlap is monotone non-increasing in the jitter and acts at any gate position); times the same ZPL-fraction q^2 as
    hom_I. Limits: T -> 0 gives 1, T -> inf gives (Gamma + k)/(Gamma + k + gamma*). Full-line Gaussian-gate closed form
    (gate far from the emission onset): I_T = e^{x^2} erfc(x), x = gamma* sigma_gate, independent of Gamma, k and the position.
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


_GATE_SHAPES = ("rect", "gauss")


def _gauss_sigma(fwhm_ns):
    return float(fwhm_ns) / (2.0 * math.sqrt(2.0 * math.log(2.0)))


def _gauss_n(Gd, c, sig):
    """int_0^inf exp(-(t - c)^2/(2 sig^2)) e^{-Gd t} dt (the Gaussian-window weight of an exponential emission that cannot
    start before the excitation at t = 0)."""
    from scipy.special import log_ndtr
    return sig * math.sqrt(2.0 * math.pi) * math.exp(-Gd * c + 0.5 * Gd * Gd * sig * sig
                                                     + float(log_ndtr((c - Gd * sig * sig) / sig)))


def _gauss_F(Gd, gamma_star, c, sig):
    """F = int_0^inf int_0^inf W(t) W(t') e^{-Gd (t + t')} e^{-gamma_star |t - t'|} dt dt', W(t) = exp(-(t - c)^2/(2 sig^2)).
    Splitting at t' > t gives 2 int_0^inf W(t) e^{-(Gd - gs) t} [int_t^inf W(t') e^{-(Gd + gs) t'} dt'] dt and the inner
    integral is a normal tail, so the double integral is one smooth 1-D quadrature (log-scaled against overflow)."""
    from scipy.integrate import quad
    from scipy.special import log_ndtr
    a = Gd + gamma_star
    b = Gd - gamma_star
    k0 = -a * c + 0.5 * a * a * sig * sig + math.log(sig * math.sqrt(2.0 * math.pi))

    def f(t):
        return math.exp(-0.5 * ((t - c) / sig) ** 2 - b * t + k0 + float(log_ndtr(-(t - c + a * sig * sig) / sig)))
    hi = c + 14.0 * sig + 2.0 * sig * sig * max(-b, 0.0)
    pk = min(max(c - b * sig * sig, 0.0), hi)
    return 2.0 * quad(f, 0.0, hi, points=[pk] if 0.0 < pk < hi else None, epsabs=0.0, epsrel=1e-11, limit=400)[0]


def gated_hom(Gamma, gamma_star, t_gate_ns, T_ns, k_loss=0.0, delay_ns=0.0, gate_shape="rect", shift_ns=0.0):
    """Time-filtered HOM overlap of the photon of one excitation inside a detection gate [DR].
    Repo convention: population decay Gamma_d = Gamma + k_loss, coherence decay (Gamma_d + gamma*)/2, first-order
    coherence G1(t, t') = Gamma e^{-Gamma_d (t + t')/2} e^{-gamma* |t' - t|/2} (t, t' >= 0). The gate is a window on the
    EMISSION-time axis. Two photons from independent identical emitters (the second emitted `delay_ns` later) are
    detected through their own gates A and B, which are displaced from each other by `shift_ns` (B - A; A is moved by
    -shift/2 and B by +shift/2 about the nominal position, so the mean position is unchanged). Only the part of the
    wavepackets inside BOTH gates interferes (the gate acts on the field amplitude a(t) = sqrt(W(t)), so the
    interfering amplitude is sqrt(W_A W_B)):
        I_T = int int sqrt(W_A W_B)(t) sqrt(W_A W_B)(t') G_A(t, t') G_B(t', t) / (n_A n_B),   n = int W n(t) dt,
    which is the Grange overlap I = 2 S / N^2 restricted to the gate (delay = shift = 0).
      gate_shape "rect"  : W = 1 on [t_gate, t_gate + T] (closed form below).
      gate_shape "gauss" : W = exp(-(t - t_gate)^2/(2 sigma^2)), T = INTENSITY FWHM = 2.3548 sigma, t_gate = the gate
                           CENTRE (Bennett's 31 ps window is a Gaussian). One smooth 1-D quadrature (_gauss_F).
    Limits: T -> 0 gives 1; T -> inf, t_gate = 0 gives Gamma_d/(Gamma_d + gamma*) (= T2/(2 T1) for k = 0); a delay
    enters as e^{Gamma_d delay} before the common start (I(delay) = Gamma_d/(Gamma_d + gamma*) e^{-Gamma_d |delay|}
    for the unwindowed packet). Returns dict(n_T, I_T, num, nn): n_T = the fraction of the photons that falls in
    gate A (per excitation, relative to the radiative yield Gamma/Gamma_d); num, nn the unnormalised overlap numerator
    and n_A n_B (used to weight a jitter average by the detected-pair probability)."""
    Gamma, gamma_star, t_g, T = float(Gamma), float(gamma_star), float(t_gate_ns), float(T_ns)
    if gate_shape not in _GATE_SHAPES:
        raise ValueError(f"gate_shape {gate_shape!r} not in {_GATE_SHAPES}")
    Gd = Gamma + float(k_loss)
    d, h = float(delay_ns), 0.5 * float(shift_ns)
    if not (T > 0.0 and Gd > 0.0):
        return {"n_T": 0.0, "I_T": NAN}
    cA, cB = t_g - h, t_g + h
    if gate_shape == "gauss":
        if d != 0.0:
            raise ValueError("gate_shape 'gauss' supports delay_ns = 0 (use shift_ns for the relative gate timing)")
        sig = _gauss_sigma(T)
        nA = Gamma * _gauss_n(Gd, cA, sig)
        nB = Gamma * _gauss_n(Gd, cB, sig)
        # sqrt(W_A W_B) = exp(-(t - cbar)^2/(2 sig^2)) exp(-(cB - cA)^2/(8 sig^2)); the amplitude factor enters twice
        cbar = 0.5 * (cA + cB)
        num = Gamma ** 2 * math.exp(-((cB - cA) ** 2) / (4.0 * sig * sig)) * _gauss_F(Gd, gamma_star, cbar, sig)
        den = nA * nB
        return {"n_T": float(nA / (Gamma / Gd)), "I_T": float(num / den) if den > 0.0 else NAN, "num": float(num),
                "nn": float(den)}
    a0, a1 = cA, cA + T
    b0, b1 = cB, cB + T
    lo = max(a0, b0, 0.0, d)
    hi = min(a1, b1)
    n_A = Gamma * _int_exp(-Gd, max(a0, 0.0), a1) if a1 > max(a0, 0.0) else 0.0
    n_B = Gamma * math.exp(Gd * d) * _int_exp(-Gd, max(b0, d), b1) if b1 > max(b0, d) else 0.0
    if not (n_A > 0.0 and n_B > 0.0 and hi > lo):
        return {"n_T": n_A / (Gamma / Gd) if n_A > 0.0 else 0.0, "I_T": NAN if not (n_A > 0.0 and n_B > 0.0) else 0.0,
                "num": 0.0, "nn": n_A * n_B}
    num = Gamma ** 2 * math.exp(Gd * d) * _gate_F(Gd, gamma_star, lo, hi)
    return {"n_T": float(n_A / (Gamma / Gd)), "I_T": float(num / (n_A * n_B)), "num": float(num),
            "nn": float(n_A * n_B)}


def gated_hom_jitter(Gamma, gamma_star, t_gate_ns, T_ns, sigma_ns, k_loss=0.0, n_nodes=161, gate_shape="rect"):
    """Gate-filtered HOM overlap averaged over a Gaussian RELATIVE timing jitter of the two photons' detection gates
    (std `sigma_ns` of the displacement between gate A and gate B; Bennett 2008: the Stark sweep jitters the gate
    position) [DR/A]. Both photons are emitted together (pulsed excitation); a displaced gate keeps only the common
    part of the two windows, so I_T(shift) is non-increasing in |shift| and its Gaussian average
        I = int P(s) I_T(s) ds                (the brief's per-pair average; a pair whose windows do not overlap adds 0)
    is MONOTONICALLY non-increasing in sigma, equals gated_hom at sigma = 0, and acts even when the gate opens well
    after the excitation (the earlier emission-delay model of this function did neither). Trapezoid of the normal density
    over +/- 6 sigma, restricted to the shifts where the two gates still overlap (so sigma >> T is resolved).
    The coincidence-weighted average int P N / int P n_A n_B differs by < 0.01 for sigma <= T/2 (Bennett gate) but
    plateaus at large sigma because a gate pushed before the excitation loses its photons (clipping at t = 0), which
    would break the monotone limit; the per-pair form is kept [A]. For two independent sources with per-photon gate
    jitter s use sigma = sqrt(2) s."""
    import numpy as np
    if not sigma_ns > 0.0:
        return gated_hom(Gamma, gamma_star, t_gate_ns, T_ns, k_loss, gate_shape=gate_shape)["I_T"]
    sg = float(sigma_ns)
    # I_T(s) vanishes beyond the gate overlap: |s| < T (top-hat) / e^{-s^2/(4 sigma_g^2)} (Gaussian, cut at 12 sigma_g); the
    # nodes cover min(6 sigma, that support), so a jitter >> T is not aliased onto a few nodes
    support = float(T_ns) if gate_shape == "rect" else 12.0 * _gauss_sigma(T_ns)
    L = min(6.0 * sg, support)
    s_ = np.linspace(-L, L, int(n_nodes))
    pdf = np.exp(-0.5 * (s_ / sg) ** 2) / (sg * math.sqrt(2.0 * math.pi))
    v = np.empty(s_.size)
    for i, si in enumerate(s_):
        it = gated_hom(Gamma, gamma_star, t_gate_ns, T_ns, k_loss, gate_shape=gate_shape, shift_ns=float(si))["I_T"]
        v[i] = it if math.isfinite(it) else 0.0
    return float(np.trapezoid(pdf * v, s_))


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


def single_photon_rate_corrected(rate_Hz, g2_raw):
    """Multi-photon-corrected single-photon rate R_1 = R_det sqrt(1 - g2_raw) (Kolatschek Ref 38, Nowak Ref 31 convention)."""
    return float(rate_Hz) * math.sqrt(1.0 - float(g2_raw))


# ---- r4-incoh: peak-window g2 with a slow single-photon tail (Fable R4 review 3.6, 3.10)
def binned_g2_tail(g2_0, components, T_bin_ns):
    """Central-bin g2 of a pulsed histogram whose bin equals the period T when the single-photon emission-time density is a
    sum of exponentials f(t) = sum_i a_i gamma_i exp(-gamma_i t), `components` = [(a_i, gamma_i/ns)], sum a_i = 1 [DR]:
        g2_bin = 1 - (1 - g2_0) c0,   c0 = 2 sum_ij a_i a_j gamma_i/(gamma_i + gamma_j) (1 - exp(-gamma_j T/2))
    (autocorrelation of f integrated over |tau| < T/2). One component reproduces `binned_g2` (IRF-free); a small slow
    fraction a_s adds g2_bin ~ g2_0 + 2 a_s exp(-gamma_s T/2) (the neighbours' tails leaking into the central window)."""
    T = float(T_bin_ns)
    comp = [(float(a), float(g)) for a, g in components if float(a) > 0.0]
    if not (T > 0.0 and comp and all(g > 0.0 for _, g in comp)):
        return NAN
    tot = sum(a for a, _ in comp)
    c0 = 0.0
    for ai, gi in comp:
        for aj, gj in comp:
            c0 += (ai / tot) * (aj / tot) * gi / (gi + gj) * (-math.expm1(-gj * T / 2.0))
    return float(1.0 - (1.0 - float(g2_0)) * 2.0 * c0)


def _bs_factor(bs_R=0.5):
    """2RT/(R^2 + T^2), T = 1 - R: the fringe-contrast factor of an unbalanced beam splitter (1 at R = 1/2) [DR]. Used only by the
    declared santori2002 estimator; the engine laws below carry the exact g2 dependence of the splitter."""
    R = float(bs_R)
    T = 1.0 - R
    return 2.0 * R * T / (R * R + T * T)


def hom_v_raw(I, g2_raw, eps_int=0.0, g2_int=None, bs_R=0.5):
    """Raw visibility in the CROSS-polarised normalisation the experiments use (r4-hom, Fable 2.3; exact splitter law r4 review G):
    V_raw = 2RT rho^2 I (1 - eps_int)/(R^2 + T^2 + 2RT g2_raw), T = 1 - R, with I the solver's whole-stream overlap and
    rho^2 = (1 - g2_raw)/(1 - g2_int); g2_int=None means no background (rho = 1). Normal-ordered expansion for two independent
    identical sources into an R/T splitter: C_par = N^2 [R^2 + T^2 + 2RT g2 - 2RT I], C_perp = N^2 [R^2 + T^2 + 2RT g2], so
    V = 1 - C_par/C_perp = 2RT I/(R^2 + T^2 + 2RT g2) [DR; agrees with Ollivier et al. 2021 SI Eq. S.2 / S.5, refs/Third pass/
    ollivier2021_arxiv.txt]. R = 1/2 gives I/(1 + g2). The pre-r4 law (rho^2 I - g2_raw) normalised to the far peak; it agrees to O(g2)."""
    R = float(bs_R)
    T = 1.0 - R
    rho2 = 1.0 if (g2_int is None or not g2_int < 1.0) else (1.0 - g2_raw) / (1.0 - g2_int)
    return 2.0 * R * T * rho2 * I * (1.0 - eps_int) / (R * R + T * T + 2.0 * R * T * g2_raw)


def hom_v_corrected(V_raw, g2, eps_int=0.0, bs_R=0.5):
    """Exact inverse of `hom_v_raw` for the single-photon overlap M_s: M = V_raw (R^2 + T^2 + 2RT g2)/(2RT (1 - eps_int)(1 - g2))
    (the stream overlap is I = M_s (1 - g2) + O(g2^2) only for a distinguishable extra photon; a filtered
    stream has a first-order excess, see the cap in _derive). With a background and g2 = g2_raw this equals I/(1 - g2_signal)
    exactly (Tomm 0.925-type values: with eps = 0, R = 1/2 it is the old I/(1 - g2)) [DR]."""
    R = float(bs_R)
    T = 1.0 - R
    return V_raw * (R * R + T * T + 2.0 * R * T * g2) / (2.0 * R * T * (1.0 - eps_int) * (1.0 - g2))


HOM_CORRECTIONS = ("engine", "santori2002", "ollivier2021")


def hom_v_corrected_by(correction, V_raw, g2, eps_int=0.0, bs_R=0.5):
    """The corrected visibility under the convention a paper declares (card item field `correction:`) [A]:
    engine = the exact inverse of the engine's raw law; ollivier2021 = (V + g2)/(1 - g2), Ollivier et al., PRL 126,
    063602 (2021), applied to the raw V as printed (no beam-splitter / interferometer term); santori2002 = that same
    g2 treatment after removing the beam-splitter and interferometer terms, (V/[(1 - eps) 2RT/(R^2 + T^2)] + g2)/(1 - g2),
    the Santori et al. (2002) type with R, T, eps and g2.
    r4 red-team link 7: formula not verified against Santori 2002 (text not on disk). The on-disk refs/Third pass/santori2004.txt
    is "Submicrosecond correlations in photoluminescence from InAs quantum dots" (the blinking paper, PRB 2004), not the 2002
    Nature HOM paper, and contains no visibility correction; it cannot source this convention. The Ollivier Eq. 2 form is the
    Msn = 0 limit of Ollivier Eq. 1 (justified for the in-pulse photon an etalon removes) [A]."""
    if correction == "engine":
        return hom_v_corrected(V_raw, g2, eps_int, bs_R)
    if correction == "ollivier2021":
        return (V_raw + g2) / (1.0 - g2)
    if correction == "santori2002":
        return (V_raw / ((1.0 - eps_int) * _bs_factor(bs_R)) + g2) / (1.0 - g2)
    raise ValueError(f"unknown HOM correction {correction!r}; valid: {HOM_CORRECTIONS}")


def cavity_sideband_parts(ph, T_K, kappa_meV, detuning_meV, gamma_zpl_meV):
    """(Z, t_sb, sb) of the polaron spectrum S = Z Lor + S_sb seen by a Lorentzian cavity (peak-normalised, FWHM kappa) at
    detuning_meV = E_ZPL - E_mode: t_sb = int S_sb T_cav(w + Delta) dw, sb = int S_sb dw. The cavity-channel rate is the overlap
    of the WHOLE spectrum with the cavity response, Gamma_cav = Gamma_0 F_P (Z L + t_sb)/(Z + sb) (r4 review A) [DR]."""
    import numpy as np
    omega, s_sb, z = qd_gf._spectrum_parts(ph, float(T_K), float(gamma_zpl_meV))
    sb_area = float(np.trapezoid(s_sb, omega))
    xc = omega + float(detuning_meV)
    acc = (kappa_meV / 2.0) ** 2 / (xc ** 2 + (kappa_meV / 2.0) ** 2)
    t_sb = float(np.trapezoid(s_sb * acc, omega))
    return float(z), t_sb, sb_area


def cavity_sideband_filter(ph, T_K, kappa_meV, detuning_meV, gamma_zpl_meV):
    """(Z, Z_eff) of the cavity-filtered output (module doc): detuning_meV = E_ZPL - E_mode. Z_eff = (Z + t_sb)/(Z + sb) is the
    ZPL-overlap-free (L = 1) limit of the cavity-channel spectral overlap L_spec = (Z L + t_sb)/(Z + sb)."""
    z, t_sb, sb_area = cavity_sideband_parts(ph, T_K, kappa_meV, detuning_meV, gamma_zpl_meV)
    return z, float((z + t_sb) / (z + sb_area))


# ---- r4-phonon I1: Zf_window (Fable r4 adjudication section 2) -----------------------------------------------------
WINDOW_SHAPES = ("lorentzian", "gaussian", "tophat")


def zpl_window_transmission(w_meV, offset_meV, gamma_meV, shape="lorentzian"):
    """Fraction of a unit-area Lorentzian ZPL (FWHM gamma, centred `offset` from the window centre) that a spectral window of
    FWHM/width `w` passes. Lorentzian window (peak 1): (w (w + g)/4)/(off^2 + ((w + g)/2)^2), = w/(w + g) on centre (convolution of
    two Lorentzians); tophat: [atan((w/2 - off)/(g/2)) + atan((w/2 + off)/(g/2))]/pi; gaussian: Re wofz((off + i g/2)/(s sqrt2)),
    s = w/(2 sqrt(2 ln2)). w = 0 passes 0, w -> inf passes 1."""
    import math
    w, off, g = float(w_meV), float(offset_meV), float(gamma_meV)
    if w <= 0.0:
        return 0.0
    if g <= 0.0:                                                   # delta-function ZPL: the window value at the line
        x = off
        return (1.0 if abs(x) <= w / 2.0 else 0.0) if shape == "tophat" else (
            math.exp(-4.0 * math.log(2.0) * x * x / (w * w)) if shape == "gaussian"
            else (w / 2.0) ** 2 / (x * x + (w / 2.0) ** 2))
    if shape == "lorentzian":
        return (w * (w + g) / 4.0) / (off * off + ((w + g) / 2.0) ** 2)
    if shape == "tophat":
        return (math.atan((w / 2.0 - off) / (g / 2.0)) + math.atan((w / 2.0 + off) / (g / 2.0))) / math.pi
    if shape == "gaussian":
        from scipy.special import wofz
        sg = w / (2.0 * math.sqrt(2.0 * math.log(2.0)))
        return float(wofz(complex(off, g / 2.0) / (sg * math.sqrt(2.0))).real)
    raise ValueError(f"unknown filter shape {shape!r}; valid: {WINDOW_SHAPES}")


def _window_acceptance(x, w, shape):
    import numpy as np
    if shape == "lorentzian":
        return (w / 2.0) ** 2 / (x ** 2 + (w / 2.0) ** 2)
    if shape == "gaussian":
        return np.exp(-4.0 * np.log(2.0) * (x / w) ** 2)
    if shape == "tophat":
        return (np.abs(x) <= w / 2.0).astype(float)
    raise ValueError(f"unknown filter shape {shape!r}; valid: {WINDOW_SHAPES}")


def window_filter_parts(ph, T_K, w_meV, offset_meV, gamma_zpl_meV, shape="lorentzian", kappa_meV=None, detuning_meV=0.0,
                        L_zpl=1.0):
    """Zf of a spectral window (module doc / adjudication I1) on the dot spectrum Z delta_Gamma + S_sb.
    offset_meV = ZPL minus window centre (the setup's filter_offset_meV); the sideband sits at offset + omega.
    Returns dict z, T_zpl, Zf (window-passed share of the spectrum seen by the window, normalised by z + sideband area as in
    cavity_sideband_filter: w -> 0 gives Z, w -> inf gives 1), Zf_cav (cavity-collected light: ZPL passed whole, sideband
    through a Lorentzian cavity kappa at `detuning_meV` = E_ZPL - E_mode, then the window; None without kappa) and Z_eff_cav
    (the cavity-only Z_eff of cavity_sideband_filter). r4 review A: L_zpl = the ZPL-cavity overlap L weights the ZPL term of the
    cavity light, Z_eff_cav = (Z L + b_c)/(Z + sb), Zf_cav = (Z L T_zpl + b_w)/(Z + sb); the default L_zpl = 1 is the pre-review value."""
    import numpy as np
    omega, s_sb, z = qd_gf._spectrum_parts(ph, float(T_K), float(gamma_zpl_meV))
    z = float(z)
    sb_area = float(np.trapezoid(s_sb, omega))
    w = float(w_meV)
    t_zpl = zpl_window_transmission(w, offset_meV, gamma_zpl_meV, shape)
    t_sb = float(np.trapezoid(s_sb * _window_acceptance(omega + float(offset_meV), w, shape), omega)) if w > 0.0 else 0.0
    out = {"z": z, "T_zpl": t_zpl, "Zf": (z * t_zpl + t_sb) / (z + sb_area), "Zf_cav": None, "Z_eff_cav": None}
    if kappa_meV is not None:
        acc = (kappa_meV / 2.0) ** 2 / ((omega + float(detuning_meV)) ** 2 + (kappa_meV / 2.0) ** 2)
        b_c = float(np.trapezoid(s_sb * acc, omega))
        b_w = float(np.trapezoid(s_sb * acc * _window_acceptance(omega + float(offset_meV), w, shape), omega)) if w > 0.0 else 0.0
        out["Z_eff_cav"] = (z * float(L_zpl) + b_c) / (z + sb_area)
        out["Zf_cav"] = (z * float(L_zpl) * t_zpl + b_w) / (z + sb_area)           # window-passed, absolute (unit-area spectrum)
    return out
# ---- end r4-phonon I1 ----------------------------------------------------------------------------------------------


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


_H_MEV_PER_GHZ = 4.135667696e-3                      # h in meV/GHz [V: CODATA]


def _pulse_end_margin_ns(plan) -> float:
    """Time (ns) from the end of the last pulse to the end of the solver window (r4 review E). A Gaussian ends where its field is
    e^-8 of the peak, t_c + 4 sigma (window end t_c + 8 sigma); a rectangle ends with its window; chirped / sinc / explicit-field
    pulses are treated as filling their window (margin 0, conservative). Window as polaron.PulseTrain: W = 2 max(|t0| + half)."""
    from . import polaron as _pl
    from .. import resonant_source as _rs
    ps = plan.pulses
    if not ps or plan.inputs.get("slm_field"):
        return 0.0
    half = [abs(p.t0_ps) * 1e-3 + _pl._half_window(p) for p in ps]
    W = 2.0 * max(half)
    ends = []
    for p in ps:
        tc = 0.5 * W + p.t0_ps * 1e-3
        e = _pl._half_window(p)
        if p.shape == "gauss" and float(p.chirp_meV_ps) == 0.0:
            e = 0.5 * _rs.GAUSS_WINDOW_SIGMA * _rs.gaussian_sigma(p.fwhm_ps * 1e-3)       # 4 sigma
        ends.append(tc + e)
    return max(0.0, W - max(ends))


def _etalon_window(plan, res):
    """(w_meV, taxes_brightness) of the Lorentzian etalon on the detected stream (setup.filter_fwhm_GHz; (0, False) when off
    or not applied by the solver, i.e. res.extra["filter_applied"] is not set: engine.py then lists the filter in
    invalid_reasons). r4 red-team item 1; r4 review fix 1."""
    B = float(plan.inputs.get("filter_fwhm_GHz", 0.0) or 0.0)
    if B <= 0.0 or not res.extra.get("filter_applied"):
        return 0.0, False
    return B * _H_MEV_PER_GHZ, str(plan.inputs.get("filter_applies_to", "g2_hom")) == "g2_hom_brightness"


def _derive(c: Context) -> dict:
    em, sr, rt, plan, res, su = c.em, c.sr, c.rt, c.plan, c.res, c.setup
    d: dict[str, Any] = {}
    why: dict[str, str] = {}
    warnings: list = []
    Gamma, k = float(rt.Gamma_ns), float(rt.k_nr_ns)
    ph = em.phonon
    Z = float(qd_gf.zpl_weight(ph, c.T_K)) if ph is not None else 1.0
    kap = sr.kappa_meV
    if ph is not None and kap is not None and kap > 0.0 and sr.E_mode_meV is not None:
        # r4-cavity A2: no Z_filt double count for CAVITY-collected photons. The cavity is the spectral filter, so the cavity
        # channel's sideband loss is already in beta (rates L_ph) or in the entered/lifetime-derived F_P; its Z_filt is 1
        # without a ZPL filter, and a ZPL filter removes only Z/Z_eff of the cavity light (1 for a narrow pillar). The
        # leaky (beta0) channel is not cavity filtered: it keeps the setup bookkeeping (Z with a ZPL filter, else 1; q = Z).
        # r4 review A: the dot width the cavity sees is w_dot (no Purcell channel), as in rates.combine; the cavity light's ZPL
        # share carries the ZPL-cavity overlap L, q_c = Z L/(Z L + t_sb), L_spec = (Z L + t_sb)/(Z + sb).
        w_dot_meV = float(rt.dot_linewidth_ueV) * 1e-3
        det_meV = float(em.E_X_meV) - float(sr.E_mode_meV)
        Zs, t_sb, sb_ = cavity_sideband_parts(ph, c.T_K, float(kap), det_meV, w_dot_meV)
        Lc = float(rt.L)
        Z_eff = (Zs * Lc + t_sb) / (Zs + sb_)              # L_spec: cavity-channel spectral overlap
        qc0 = Zs * Lc / (Zs * Lc + t_sb) if (Zs * Lc + t_sb) > 0.0 else 1.0
        w_et, bright_et = _etalon_window(plan, res)
        if su.zpl_filter:
            Zc, qc, Zl, ql = qc0, 1.0, Z, 1.0
            P_c, P_l = Zc, Zl
        elif w_et > 0.0:
            # r4 red-team item 1: the etalon (setup.filter_fwhm_GHz, Lorentzian, ZPL-centred [A]) sets q from its own sideband
            # transmission; the ZPL Lorentzian transmission is in the cascade's n_ph, so the brightness gets only the sideband tax.
            wp = window_filter_parts(ph, c.T_K, w_et, 0.0, w_dot_meV, "lorentzian", float(kap), det_meV, Lc)
            Tz = wp["T_zpl"]
            Zc_w = wp["Zf_cav"] / Tz if Tz > 0.0 else Z
            Zl_w = Z + (wp["Zf"] - Z * Tz) / Tz if Tz > 0.0 else Z
            Zc, Zl = (Zc_w / wp["Z_eff_cav"], Zl_w) if bright_et else (1.0, 1.0)
            qc, ql = Z * Lc * Tz / wp["Zf_cav"], Z * Tz / wp["Zf"]
            # r4 review B: the HOM path is ALWAYS filtered, whatever the brightness convention
            P_c, P_l = Zc_w / wp["Z_eff_cav"], Zl_w
        else:
            Zc, qc, Zl, ql = 1.0, qc0, 1.0, Z
            P_c, P_l = 1.0, 1.0
            w_win = float(_read(c, "setup", "filter_width_meV", 0.0) or 0.0)
            if w_win > 0.0:
                # r4-phonon I1 (Zf_window): the setup's spectral window also truncates the ZPL and the sideband. Cavity light:
                # window-passed share relative to the cavity-collected light; leaky channel: the window on the whole spectrum.
                wp = window_filter_parts(ph, c.T_K, w_win, float(_read(c, "setup", "filter_offset_meV", 0.0) or 0.0),
                                         w_dot_meV, str(_read(c, "setup", "filter_shape", "lorentzian")),
                                         float(kap), det_meV, Lc)
                Zc, qc = wp["Zf_cav"] / wp["Z_eff_cav"], Z * Lc * wp["T_zpl"] / wp["Zf_cav"]
                Zl, ql = wp["Zf"], Z * wp["T_zpl"] / wp["Zf"]
                P_c, P_l = Zc, Zl
        G0r = float(rt.Gamma0_ns)
        sh_c = float(rt.F_P_eff) * G0r / Gamma if (Gamma > 0.0 and G0r > 0.0) else 1.0
        sh_b = float(rt.beta0) * float(rt.Gamma_bg_ns) / Gamma if Gamma > 0.0 else 0.0
        w_c = sh_c / (sh_c + sh_b) if (sh_c + sh_b) > 0.0 else 1.0
        Zf = w_c * Zc + (1.0 - w_c) * Zl
        # r4 review B: q weights the channels by the light each PASSES the HOM-path filter, q = sum w P q / sum w P
        qden = w_c * P_c + (1.0 - w_c) * P_l
        q = (w_c * P_c * qc + (1.0 - w_c) * P_l * ql) / qden if qden > 0.0 else w_c * qc + (1.0 - w_c) * ql
    elif su.zpl_filter:
        Zf, q = Z, 1.0
    elif _etalon_window(plan, res)[0] > 0.0:
        w_et, bright_et = _etalon_window(plan, res)           # r4 red-team item 1 (see the cavity branch)
        wp = window_filter_parts(ph, c.T_K, w_et, 0.0, lb_hbar_ueV() * 1e-3 * (Gamma + k + float(rt.gamma_star_ns)), "lorentzian") \
            if ph is not None else None
        if wp is None:
            Zf, q = 1.0, 1.0
        else:
            Tz = wp["T_zpl"]
            Zf = (Z + (wp["Zf"] - Z * Tz) / Tz) if (bright_et and Tz > 0.0) else 1.0
            q = Z * Tz / wp["Zf"]
    else:
        Zf, q = 1.0, Z
        w_win = float(_read(c, "setup", "filter_width_meV", 0.0) or 0.0) if ph is not None else 0.0
        if w_win > 0.0:                                 # r4-phonon I1 (Zf_window): the window truncates the ZPL and the sideband
            wp = window_filter_parts(ph, c.T_K, w_win, float(_read(c, "setup", "filter_offset_meV", 0.0) or 0.0),
                                     lb_hbar_ueV() * 1e-3 * (Gamma + k + float(rt.gamma_star_ns)),
                                     str(_read(c, "setup", "filter_shape", "lorentzian")))
            Zf, q = wp["Zf"], Z * wp["T_zpl"] / wp["Zf"]
    beta = float(res.extra.get("beta", rt.beta))
    eta_pol = float(su.eta_pol) * float(sr.pol_factor)
    pX = float(plan.inputs.get("p_X", 1.0) or 1.0)     # r4-incoh: X share of the loaded pairs (non-resonant pump); 1 = bit-identical
    n_ph = float(res.n_ph) * pX
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
    # R2 FIX-1a (Fable-A Q1): lifetime ratio tau_0/tau. "radiative" = Gamma/Gamma0 = F_bg + F_P L (the additional-rate
    # convention, total = F_bg + F_P L); "intrinsic" = tau_0/tau with tau_0 = 1/(Gamma0 + k) the bulk lifetime and tau the
    # tau_ps observable (JC-aware), i.e. what a TRPL experiment quotes as "F_P = tau_bulk/tau". They coincide for k_nr = 0.
    G0 = float(em.gamma0_ns)
    d["tau_ratio_rad"] = Gamma / G0 if G0 > 0.0 else NAN
    d["tau_ratio"] = (1e3 / (G0 + k)) / d["tau_ps"] if (G0 + k) > 0.0 and d["tau_ps"] > 0.0 else NAN
    # cw resonant drive: MvdS on-resonance saturation parameter s_0 = I/I_sat = 2 Omega^2/Gamma_pop^2 (Gamma_pop = Gamma + k) and
    # the code's s = Omega^2 T1 T2 (T1 = 1/(Gamma + k), T2 = 2/(Gamma + k + Gamma_ph)); they are equal at Gamma_ph = 0 [DR, ch02]
    om_cw = float(plan.omega_cw_ns)
    if plan.family == "coherent" and plan.timing == "cw" and om_cw > 0.0:
        gp = Gamma + k
        d["s_0"] = 2.0 * om_cw ** 2 / gp ** 2
        d["s_sat"] = om_cw ** 2 / (gp * 0.5 * (gp + float(rt.gamma_star_ns)))
    else:
        d["s_0"] = d["s_sat"] = NAN
        why["s_0"] = why["s_sat"] = "saturation parameters exist only for the cw resonant drive"
    d["coherent_fraction"] = float(res.extra["coherent_fraction"]) if "coherent_fraction" in res.extra else NAN
    if "coherent_fraction" not in res.extra:
        why["coherent_fraction"] = "only the cw resonant drive reports the coherently scattered fraction"
    d["eta_extract"] = beta * d["eta_out"]
    bd = _read(c, "structure", "beta_direct", 0.0)
    if bd is not None and float(bd) > 0.0:
        d["beta_entered"] = float(bd)
    else:
        d["beta_entered"] = NAN
        why["beta_entered"] = "the structure has no entered beta (structure.beta_direct absent or 0)"
    d["inversion"] = NAN if res.inversion is None else float(res.inversion)
    # r4-hom (Fable 2.4): `inversion` (convention count = intrinsic) is the expected excitation count n_exc (P_X(end) + photons emitted
    # in the window), which does not depend on how long the window runs after the pulse; the end-of-window population stays as the
    # observable P_X_end.
    if "n_exc" in res.extra and res.inversion is not None and math.isfinite(float(res.extra["n_exc"])):
        d["inversion"] = float(res.extra["n_exc"])
    # r4 review E: convention population = the exciton population right after the pulse, <= 1 by definition (a simulated / theory
    # "final occupation"; n_exc adds the photons emitted during the pulse and the re-excited ones). Once the drive is off the
    # two-level population decays exactly as exp(-(Gamma + k) t), so P_pulse_end = P_X_end exp((Gamma + k) m) with m the time from
    # the end of the last pulse (field e^-8 of the peak = t_c + 4 sigma for a Gaussian) to the window end [DR, Fable r4 review E].
    d["inversion_pop"] = NAN
    if "P_X_end" in res.extra and math.isfinite(float(res.extra["P_X_end"])) and res.inversion is not None:
        if "cascade" in plan.needs or "n_cav_end" in res.extra:
            why["inversion_pop"] = ("the pulse-end population is a closed form for the two-level (lindblad / polaron) solvers; the "
                                    "jc and TPE cascade solvers would need the propagator state at the pulse end")
        else:
            d["inversion_pop"] = float(res.extra["P_X_end"]) * math.exp((Gamma + k) * _pulse_end_margin_ns(plan))
    else:
        why["inversion_pop"] = "the solver has no pulse-end population (incoherent / cw drive)"
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
    d["P_X_end"] = float(res.extra["P_X_end"]) if "P_X_end" in res.extra else (
        NAN if res.inversion is None else float(res.inversion))
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
        S = (float(flux) * pX * 1e9 * eta_blink * eta_pol * beta * d["eta_out"] * Zf * float(su.eta_fibre)
             * float(su.eta_optics) * float(su.eta_det)) if flux is not None else NAN
    d["rate_detected_Hz"] = S
    # FX3b: photon rates into the first lens / fibre (before eta_optics, eta_det), cw (flux = Gamma n_X,ss) and pulsed (f_rep B)
    if pulsed and f_rep and math.isfinite(d["B_lens"]):
        d["B_first_lens_Hz"], d["B_fibre_Hz"] = f_rep * d["B_lens"], f_rep * d["B_fibre"]
    elif not pulsed and res.extra.get("flux_ns") is not None:
        Bl = float(res.extra["flux_ns"]) * pX * 1e9 * eta_blink * eta_pol * beta * d["eta_out"] * d["T_NA"] * Zf
        d["B_first_lens_Hz"], d["B_fibre_Hz"] = Bl, Bl * float(su.eta_fibre)
    else:
        d["B_first_lens_Hz"] = d["B_fibre_Hz"] = NAN
        why["B_first_lens_Hz"] = why["B_fibre_Hz"] = "no photon flux available from the solver / plan"
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
    # FX3b background-corrected g2 = the paper formula (g2_raw - (1 - rho^2))/rho^2 with the BACKGROUND rho = (1 - bg_frac)
    # /(1 + laser_leak): the drive's embedded bg_frac (already inside g2_0) and the laser leak; dark counts are NOT
    # corrected (the "dark-count corrected" value is a separate convention). It is the exact inverse of the mixing law.
    bg_emb = float(plan.inputs.get("bg_frac", 0.0) or 0.0)
    rho_bg = (1.0 - bg_emb) / (1.0 + leak)
    d["rho_background"] = rho_bg
    if math.isfinite(d["g2_raw"]) and rho_bg > 0.0:
        d["g2_bgc"] = (d["g2_raw"] - (1.0 - rho_bg * rho_bg)) / (rho_bg * rho_bg)
    else:
        d["g2_bgc"] = NAN
        why["g2_bgc"] = "background-corrected g2 needs a defined raw g2 and a background rho > 0"
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
    # r4-incoh ----------------------------------------------------------------------------------------------------------
    # g2_far = g2_on / eta_blink (papers that normalise to FAR peaks with visible blinking: two-state emitter, tau_blink >> T_rep)
    # (intrinsic = background free). r4 review F: convention raw = the long-delay value a detector gives with the background,
    # 1 - rho^2 (1 - g2/eta_blink), rho = S/(S + B), S the blink-averaged detected rate (the blink division acts on the intrinsic g2
    # first, the background mixing after) [DR].
    if eta_blink > 0.0:
        d["g2_far"] = g2 / eta_blink
        if leak == 0 and su.dark_cps == 0:
            d["g2_far_raw"] = d["g2_far"]
        elif math.isfinite(bg):
            d["g2_far_raw"] = g2_with_background(d["g2_far"], S, bg)
            if not math.isfinite(d["g2_far_raw"]):
                why["g2_far_raw"] = "no detected signal (S <= 0) with a non-zero background: g2_far_raw undefined"
        else:
            d["g2_far_raw"] = NAN
            why["g2_far_raw"] = "background rate needs the detected signal rate S, which is unavailable"
    else:
        d["g2_far"] = d["g2_far_raw"] = NAN
        why["g2_far"] = why["g2_far_raw"] = "needs a non-zero blinking on-fraction"
    # single-photon rate corrected for the multi-photon share, R_1 = R_det sqrt(1 - g2_raw) (Kolatschek Ref 38 / Nowak Ref 31)
    if math.isfinite(S) and math.isfinite(d["g2_raw"]) and d["g2_raw"] < 1.0:
        d["rate_single_corr_Hz"] = single_photon_rate_corrected(S, d["g2_raw"])
    else:
        d["rate_single_corr_Hz"] = NAN
        why["rate_single_corr_Hz"] = "needs a detected rate and a raw g2 below 1"
    # peak-window g2 with a slow single-photon tail (drive inputs a_slow, tau_slow_ps); no tail = the one-exponential binned g2
    a_sl = float(plan.inputs.get("a_slow", 0.0) or 0.0)
    if pulsed and f_rep:
        comp = [(1.0 - a_sl, Gamma + k + k_off)]
        if a_sl > 0.0:
            comp.append((a_sl, 1e3 / float(plan.inputs["tau_slow_ps"])))
        d["g2_window"] = binned_g2_tail(g2, comp, 1e9 / f_rep)
    else:
        d["g2_window"] = NAN
        why["g2_window"] = "the peak-window g2 needs a pulsed drive (window = pulse period)"
    # fitted X-line decay time of the capacity-2 reservoir cascade (incoherent pulsed drives)
    d["tau_fit_ps"] = d["t_mean_ps"] = NAN
    if plan.family == "incoherent" and pulsed and float(plan.relax_ns) > 0.0 and float(plan.pump_rate_ns) > 0.0:
        try:
            from . import solvers as _sv
            ratio = float(em.gamma_xx_ratio)
            gXX = ratio * float(getattr(rt, "Gamma0_ns", Gamma))
            # r4 review D: a ladder plan (tau_res / q_T > 0) uses its OWN generator (blocked pairs wait and refill the emptied dot,
            # NR as the plan's N_R); the S3d scheme keeps the stand-alone estimator where surplus pairs vanish [A] (= the ladder
            # with tau_res -> 0 for blocked pairs, validated on wijitpatima)
            lcfg = _sv._ladder_cfg(plan)
            t_, I_ = _sv.x_line_decay(float(plan.pump_rate_ns) * float(plan.tau_on_ns), float(plan.relax_ns), Gamma, k,
                                      gXX, k * ratio, float(plan.inputs.get("tau_res_ps", 0.0) or 0.0) * 1e-3,
                                      NR=(lcfg["N_R"] or None) if lcfg else None, full_loss=lcfg is None)
            tf, tm = _sv.tau_fit_from_decay(t_, I_)
            d["tau_fit_ps"], d["t_mean_ps"] = 1e3 * tf, 1e3 * tm
        except (ValueError, RuntimeError, FloatingPointError) as exc:       # pragma: no cover - degenerate inputs
            why["tau_fit_ps"] = why["t_mean_ps"] = f"decay fit failed: {exc}"
    if not math.isfinite(d["tau_fit_ps"]):
        why.setdefault("tau_fit_ps", "needs an incoherent pulsed drive with a capture time (above_band / p_shell); a cw drive or a "
                                     "coherent single-X pulse decays as the single exponential tau_ps, and the TPE cascade X line "
                                     "(two exponentials) is not fitted")
        why.setdefault("t_mean_ps", why["tau_fit_ps"])
    I1 = float(res.I_single)
    d["hom_I"] = I1 * q * q if math.isfinite(I1) else NAN
    # r4-hom (Fable 2.5): two-timescale dephasing. The solver's I carries the FAST dephasing gamma* only; the slow spectral
    # diffusion (material.sigma_sd_ueV, tau_sd_ns) reduces the overlap of two photons setup.hom_delay_ns apart by
    # F_sd = hom_vs_separation(dt)/hom_vs_separation(0) (OU walk; -> 1 for dt << tau_sd). 0 (default) = off, bit-identical.
    d["hom_sd_factor"] = 1.0
    hdly = float(_read(c, "setup", "hom_delay_ns", 0.0) or 0.0)
    if hdly > 0.0 and math.isfinite(d["hom_I"]):
        from .. import resonant_source as _rs
        from .. import lindblad as _lb
        _args = (Gamma, float(res.extra.get("hom_gamma_star_fast_ns", rt.gamma_star_ns)), float(_lb.mev_to_rate(float(em.sigma_sd_ueV) * 1e-3)), float(em.tau_sd_ns))
        _f0 = float(_rs.hom_vs_separation(0.0, *_args))
        _f1 = float(_rs.hom_vs_separation(hdly, *_args))
        if _f0 > 0.0 and math.isfinite(_f1):
            d["hom_sd_factor"] = _f1 / _f0
            d["hom_I"] *= d["hom_sd_factor"]
    hg = plan.inputs.get("hom_gate")
    if hg:
        sg = float(hg.get("jitter_ns", 0.0))
        gs = float(rt.gamma_star_ns)
        shape = str(hg.get("shape", "rect"))
        Ig = gated_hom_jitter(Gamma, gs, float(hg["start_ns"]), float(hg["width_ns"]), sg, k + k_off, gate_shape=shape)
        d["hom_I_gated"] = Ig * q * q if math.isfinite(Ig) else NAN
        d["hom_gate_fraction"] = gated_hom(Gamma, gs, float(hg["start_ns"]), float(hg["width_ns"]), k + k_off,
                                           gate_shape=shape)["n_T"]
    else:
        d["hom_I_gated"] = NAN
        why["hom_I_gated"] = "no HOM gate set; set the HOM gate width (drive.hom_gate_ps > 0)"
        d["hom_gate_fraction"] = 1.0             # no gate: every emitted photon is inside the (whole-period) window
    # R2 dead controls (F10b): the delivered pump Fano factor of an electrical drive, F_dot = 1 - eta_cap (1 - F_pump)
    # (loading.f8b_thin_fano; F_pump from the mechanism: shot-noise rail 1, impedance divider Z_j/(Z_j + Z_d) [E] with the
    # finite-regulation correction, a priced Coulomb-blockade element eps_cyc, a pulsed rail F_rail). REPORTED: the counting
    # solvers do not feed it back into g2 (their notes say so).
    itf = plan.interface
    if itf is not None:
        d["F_dot"] = float(itf.F_dot)
    else:
        d["F_dot"] = NAN
        why["F_dot"] = "only the electrical drives carry a pump interface (Fano factor of the carrier supply)"
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
        why["hom_I"] = why["hom_V"] = ("a cw stream has no per-pulse HOM overlap" if not pulsed
                                       else "solver has no coherence information (rate solver)")
    eps = float(_read(c, "setup", "eps_int", 0.0))
    bsR = float(_read(c, "setup", "bs_R", 0.5))
    corr_keys = {"engine": "hom_V_corr", "santori2002": "hom_V_corr_santori2002", "ollivier2021": "hom_V_corr_ollivier2021"}
    # r4 red-team item 2: filter_applies_to = hom puts the etalon in the HOM path only. The reported g2 stays unfiltered, but the
    # g2 inside the HOM law (the multi-photon share of the stream in the interferometer) is that of the filtered stream.
    g2_hm = g2_hs = None
    if str(plan.inputs.get("filter_applies_to", "")) == "hom" and "g2_filtered" in res.extra:
        g2_hs = float(res.extra["g2_filtered"])
        if leak == 0 and su.dark_cps == 0:
            g2_hm = g2_hs
        else:
            # r4 review C: the stream entering the HOM interferometer is the etalon-filtered signal S_hom = S n_f/n_u; the dark counts
            # arise in the detector behind the filter (unchanged) and the laser leak is a fraction of the detected signal [A]
            nf_, nu_ = res.extra.get("n_ph_filtered"), res.extra.get("n_ph_unfiltered")
            if nf_ is not None and nu_ is not None and math.isfinite(float(nf_)) and float(nu_) > 0.0 and math.isfinite(S):
                S_hom = S * float(nf_) / float(nu_)
                g2_hm = g2_with_background(g2_hs, S_hom, leak * S_hom + dark)
            else:
                g2_hm = NAN
    g2_r = d["g2_raw"] if g2_hm is None else g2_hm
    if math.isfinite(d["hom_I"]) and math.isfinite(g2_r):
        # FX4: the overlap I is the SIGNAL stream's, so rho^2 = (1 - g2_raw)/(1 - g2_signal) must use the signal g2 (before
        # the drive-embedded bg_frac); g2_0 already contains that background and would double-count it (final review 1).
        g2_sig = float(res.extra.get("g2_signal", g2)) if g2_hs is None else g2_hs
        d["hom_V_raw"] = hom_v_raw(d["hom_I"], g2_r, eps, g2_int=g2_sig, bs_R=bsR)
        # r4-hom (Fable 2.3): corrected = the exact inverse of the raw law (eps_int and the splitter term are removed too).
        # r4 Fable filtered-HOM note (.workers/devlib/r4-fable-filtered-hom.md): the stream identity I = M_s (1 - g2) + O(g2^2)
        # holds only for a DISTINGUISHABLE extra photon. In general I = M_s (1 - g2) + g2 [Tr(sigma G_2) - M_s + M_s p0] + O(g2^2);
        # a re-excitation photon that overlaps the single photon (spectral filter: |s|^2 = Gamma/(Gamma + kappa)) gives a
        # FIRST-order excess g2 [2|s|^2/(1 + |s|^2) + p0] [DR]. The only rigorous bound is I <= 1 (G1 is a positive kernel), so
        # each corrected convention is capped at its I = 1 value: above the cap is an inconsistency (NaN); between 1 and the
        # cap is physics and is clipped to 1 with a note (the paper-style estimator is biased high).
        V_raw_max = hom_v_raw(1.0, g2_r, eps, g2_int=g2_sig, bs_R=bsR)
        for cname, key in corr_keys.items():
            v = hom_v_corrected_by(cname, d["hom_V_raw"], g2_r, eps, bsR)
            v_max = hom_v_corrected_by(cname, V_raw_max, g2_r, eps, bsR)
            if v > v_max * (1.0 + 1e-6):
                v = NAN
                why[key] = "corrected visibility exceeds its I = 1 cap: solver overlap I > 1 (numerics)"
            elif v > 1.0:
                warnings.append(f"{key} = {v:.6f} exceeds 1 by {v - 1.0:.2e} <= cap {v_max - 1.0:.2e} = its value at I = 1: the "
                                "multi-photon component overlaps the single-photon mode (first order in g2, "
                                "I = M_s(1 - g2) + g2[2|s|^2/(1+|s|^2) + p0]); the paper estimator is biased high, clipped to 1")
                v = 1.0
            d[key] = v
    else:
        d["hom_V_raw"] = NAN
        for key in corr_keys.values():
            d[key] = NAN
    if bool(plan.inputs.get("cross_pol")) and eta_pol >= 1.0 and "vsys" not in res.extra:
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


# display symbols (KaTeX TeX; Scully-Zubairy / Metcalf-van der Straten notation), keyed by the observable id
_OBS_TEX = {
    "lambda_nm": r"\lambda",
    "F_P_eff": r"F_{\mathrm P}L",
    "F_P_peak": r"F_{\mathrm P}",
    "beta": r"\beta",
    "beta_entered": r"\beta",
    "eta_out": r"\eta_{\mathrm{out}}",
    "eta_extract": r"\beta\,\eta_{\mathrm{out}}",
    "T_NA": r"T_{\mathrm{NA}}",
    "eta_pol": r"\eta_{\mathrm{pol}}",
    "eta_prep": r"\eta_{\mathrm{prep}}",
    "n_exc": r"n_{\mathrm{exc}}",
    "eta_rad": r"\eta_{\mathrm{rad}}",
    "eta_rad_eff": r"\eta_{\mathrm{rad}}^{\mathrm{eff}}",
    "eta_blink": r"\eta_{\mathrm{blink}}",
    "eta_prep_x_eta_rad": r"\eta_{\mathrm{prep}}\,\eta_{\mathrm{rad}}",
    "P_X_end": r"P_X(t_{\mathrm{end}})",
    "tau_XX_ps": r"\tau_{XX}",
    "g2_XX": r"g^{(2)}_{XX}(0)",
    "hom_XX": r"I_{\mathrm{HOM}}^{XX}",
    "n_XX": r"n_{XX}",
    "inversion": r"w",
    "Z_phonon": r"Z(T)",
    "Z_filt": r"Z_{\mathrm{filt}}",
    "hom_zpl_fraction": r"q_{\mathrm{ZPL}}",
    "gamma_star_ns": r"\Gamma_{\mathrm{ph}}",
    "tau_ratio": r"\tau_0/\tau",
    "s_0": r"s_0",
    "s_sat": r"s",
    "coherent_fraction": r"f_{\mathrm{coh}}",
    "Gamma_ns": r"\Gamma",
    "n_ph": r"n_{\mathrm{ph}}",
    "tau_ps": r"\tau",
    "tau_1e_ps": r"\tau_{1/e}",
    "tau_purcell_ps": r"\tau_{\mathrm P}",
    "tau_EL_ps": r"\tau_{\mathrm{EL}}",
    "g2_binned": r"g^{(2)}_{\mathrm{bin}}(0)",
    "g2_window": r"g^{(2)}_{\mathrm{win}}(0)",
    "g2_far": r"g^{(2)}_{\mathrm{far}}(0)",
    "rate_single_corr_Hz": r"R_{1}",
    "tau_fit_ps": r"\tau_{\mathrm{fit}}",
    "t_mean_ps": r"\langle t\rangle",
    "spectral_overlap": r"\varepsilon_{XX}",
    "g2_irf": r"g^{(2)}_{\mathrm{IRF}}(0)",
    "hom_I_gated": r"I_{\mathrm{HOM}}^{\mathrm{gate}}",
    "hom_gate_fraction": r"f_{\mathrm{gate}}",
    "F_dot": r"F_{\mathrm{dot}}",
    "linewidth_ueV": r"\hbar\Gamma_{\mathrm{tot}}",
    "detuning_ueV": r"\delta_{\mathrm{EC}}",
    "g2_0": r"g^{(2)}(0)",
    "g2_raw": r"g^{(2)}_{\mathrm{raw}}(0)",
    "hom_I": r"I_{\mathrm{HOM}}",
    "hom_V": r"V_{\mathrm{HOM}}",
    "hom_V_raw": r"V_{\mathrm{HOM}}^{\mathrm{raw}}",
    "hom_V_corr": r"V_{\mathrm{HOM}}^{\mathrm{corr}}",
    "B_first_lens": r"B_{\mathrm{lens}}",
    "B_lens": r"B_{\mathrm{lens}}",
    "B_fibre": r"B_{\mathrm{fibre}}",
    "rate_detected_Hz": r"R_{\mathrm{det}}",
    "B_first_lens_Hz": r"\dot B_{\mathrm{lens}}",
    "B_fibre_Hz": r"\dot B_{\mathrm{fibre}}",
}


def _o(id_, label, unit, deps, conventions=("intrinsic",), passthrough_of=None, key="", by_conv=None):
    return register(_Obs(id=id_, label=label, unit=unit, deps=tuple(deps), conventions=tuple(conventions),
                         passthrough_of=passthrough_of, key=key or id_, by_conv=by_conv, tex=_OBS_TEX.get(id_)))


_MS = ("material", "structure")
_ALL = ("material", "structure", "drive")
_CH = ("material", "structure", "drive", "setup", "material.eta_blink", "setup.coinc_window_ns", "setup.eps_int", "setup.bs_R")
_o("lambda_nm", "Emission wavelength $\\lambda$", "nm", ("material.lambda_nm", "material.E_X_meV"),
   passthrough_of="material.lambda_nm")
_o("F_P_eff", "Effective Purcell factor $F_{\\mathrm P}L$", "1", _MS)
_o("F_P_peak", "Cavity-mode Purcell factor $F_{\\mathrm P}$ at zero detuning (structure.F_P)", "1", _MS, passthrough_of="structure.F_P")
_o("beta", "Beta factor $\\beta$", "1", _MS, passthrough_of="structure.beta")
_o("beta_entered", "Beta $\\beta$ entered directly (structure.beta_direct, e.g. phc_waveguide)", "1",
   ("structure.beta_direct",), passthrough_of="structure.beta_direct")
_o("eta_out", "Extraction into the collected mode $\\eta_{\\mathrm{out}}$", "1", ("structure",))
_o("eta_extract", "Extraction $\\beta\\,\\eta_{\\mathrm{out}}$", "1", _MS)
_o("T_NA", "NA transmission $T_{\\mathrm{NA}}$ (inside $\\eta_{\\mathrm{out}}$)", "1", ("setup.NA",))
_o("eta_pol", "Polarisation factor $\\eta_{\\mathrm{pol}}$", "1", ("setup.eta_pol", "structure"))
_o("eta_prep", "State-preparation probability $P(n \\geq 1)$", "1", _ALL)
_o("n_exc", "Expected excitation count per pulse ($P_X(t_{\\mathrm{end}})$ + photons emitted in the window; may exceed 1)", "1", _ALL)
_o("eta_rad", "Radiative efficiency $\\eta_{\\mathrm{rad}} = \\Gamma/(\\Gamma + k)$", "1", _MS)
_o("eta_rad_eff", "Effective radiative efficiency $\\eta_{\\mathrm{blink}}\\,\\Gamma/(\\Gamma + k)$", "1", _MS + ("material.eta_blink",))
_o("eta_blink", "Blinking on-fraction $\\eta_{\\mathrm{blink}}$ (brightness factor, not lifetime)", "1", ("material.eta_blink",),
   passthrough_of="material.eta_blink")
_o("eta_prep_x_eta_rad", "$\\eta_{\\mathrm{prep}}\\,\\eta_{\\mathrm{rad}}$", "1", _ALL)
_o("P_X_end", "Population $P_X(t_{\\mathrm{end}})$ at the end of the pulse window (informational; includes decay during the window)", "1", _ALL)
_o("tau_XX_ps", "XX lifetime $\\tau_{XX} = 1/(\\Gamma_{XX} + k_{XX})$ read off the cascade Liouvillian (TPE)", "ps", _ALL)
_o("g2_XX", "$g^{(2)}(0)$ of the XX-line photon stream (TPE)", "1", _ALL)
_o("hom_XX", "HOM overlap of the XX-line photons, instantaneous preparation (TPE)", "1", _ALL + ("setup",))
_o("n_XX", "XX-line photons per pulse (TPE)", "1", _ALL)
_o("inversion", "Inversion $w$ after the pulse (count = expected excitation count $n_{\\mathrm{exc}}$, may exceed 1; population = exciton "
   "population right after the pulse, $\\leq 1$)", "1", _ALL, ("intrinsic", "count", "population"),
   by_conv={"intrinsic": "inversion", "count": "inversion", "population": "inversion_pop"})
_o("Z_phonon", "ZPL weight $Z(T)$", "1", ("material",))
_o("Z_filt", "ZPL / cavity sideband filter factor $Z_{\\mathrm{filt}}$", "1", ("material", "structure", "setup.zpl_filter", "setup.filter_width_meV", "setup.filter_offset_meV", "setup.filter_shape", "setup.filter_fwhm_GHz", "setup.filter_applies_to"))
_o("hom_zpl_fraction", "ZPL fraction $q$ of the collected light ($I_{\\mathrm{HOM}} = I\\,q^2$)", "1", ("material", "structure", "setup.zpl_filter", "setup.filter_width_meV", "setup.filter_offset_meV", "setup.filter_shape", "setup.filter_fwhm_GHz", "setup.filter_applies_to"))
_o("gamma_star_ns", "Pure-dephasing linewidth $\\Gamma_{\\mathrm{ph}} = 2\\gamma_{\\mathrm{ph}} = 2/T_2^{*}$ (code gamma*; coherence decays at $\\gamma_{\\mathrm{ph}}$)", "1/ns", ("material",))
_o("Gamma_ns", "Radiative rate $\\Gamma$", "1/ns", _MS)
_o("n_ph", "Photons per pulse $n_{\\mathrm{ph}}$", "1", _ALL)
_o("tau_ps", "Lifetime $\\tau$ (slowest Liouvillian decay eigenvalue when solver jc, else Purcell)", "ps", _ALL)
_o("tau_ratio", "Lifetime ratio $\\tau_0/\\tau$ (bulk / device; convention radiative: $F_{\\mathrm{bg}} + F_{\\mathrm P}L = \\Gamma/\\Gamma_0$)", "1", _ALL,
   ("intrinsic", "radiative"), key="tau_ratio", by_conv={"intrinsic": "tau_ratio", "radiative": "tau_ratio_rad"})
_o("s_0", "On-resonance saturation parameter $s_0 = 2\\Omega^2/\\Gamma^2 = I/I_{\\mathrm{sat}}$ (cw resonant)", "1", ("material", "structure", "drive"))
_o("s_sat", "Saturation parameter $s = \\Omega^2 T_1 T_2$ (cw resonant; equals $s_0$ at $\\Gamma_{\\mathrm{ph}} = 0$)", "1", _ALL)
_o("coherent_fraction", "Coherently scattered fraction $|\\rho_{GX}|^2/\\rho_{XX}$ (cw resonant)", "1", _ALL)
_o("tau_1e_ps", "JC flux $1/e$ time $\\tau_{1/e}$ after the peak (informational; $= \\tau$ for non-jc)", "ps", _ALL)
_o("tau_purcell_ps", "Purcell lifetime $\\tau_{\\mathrm P} = 1/(\\Gamma + k)$", "ps", _MS)
_o("tau_EL_ps", "Electroluminescence decay time $\\tau_{\\mathrm{EL}} = 1/(\\Gamma + k + k_{\\mathrm{off}})$ of a gated diode (reported, never fed into $\\Gamma_0$)", "ps",
   _ALL)
_o("g2_binned", "Central-bin $g^{(2)}$ of the pulsed histogram (bin = pulse period; EL decay leaks neighbouring peaks, IRF adds more)",
   "1", _ALL + ("setup",))
_o("g2_window", "Central-window $g^{(2)}$ of the pulsed histogram with a slow single-photon tail (window = pulse period; "
   "drive a_slow / tau_slow_ps; no tail = g2_binned without IRF)", "1", _ALL + ("setup",))
_o("g2_far", "$g^{(2)}(0)$ normalised to the far peaks when blinking is visible: $g^{(2)}_{\\mathrm{on}}/\\eta_{\\mathrm{blink}}$ "
   "(intrinsic; raw = with laser leak and dark counts, $1 - \\rho^2(1 - g^{(2)}/\\eta_{\\mathrm{blink}})$)", "1",
   _ALL + ("setup.laser_leak", "setup.dark_cps", "setup.coinc_window_ns", "material.eta_blink"), ("intrinsic", "raw"),
   by_conv={"intrinsic": "g2_far", "raw": "g2_far_raw"})
_o("rate_single_corr_Hz", "Multi-photon-corrected single-photon rate $R_{\\mathrm{det}}\\sqrt{1 - g^{(2)}_{\\mathrm{raw}}}$", "Hz", _CH)
_o("tau_fit_ps", "Single-exponential fit time of the X-line decay of the capacity-2 reservoir cascade (peak, 3 decades, linear least squares)",
   "ps", _ALL)
_o("t_mean_ps", "Mean emission time $\\int t I\\,dt/\\int I\\,dt$ of the X-line decay of the capacity-2 reservoir cascade", "ps", _ALL)
_o("spectral_overlap", "X-XX spectral leakage $\\varepsilon = t_{XX}/t_X$ through the top-hat window ($\\rho$, $g^{(2)} = 1 - \\rho^2$ conventions)",
   "1", _MS + ("setup.filter_*", "setup.xx_*"), ("intrinsic", "rho", "g2"),
   key="spectral_overlap", by_conv={"intrinsic": "spectral_overlap", "rho": "rho_spectral", "g2": "g2_spectral"})
_o("g2_irf", "$g^{(2)}(0)$ of the cw dip seen through the detector IRF", "1", _ALL + ("setup",))
_o("hom_I_gated", "HOM overlap inside a detection gate (Stark-gated source), $q^2$ included", "1", _ALL + ("setup",))
_o("hom_gate_fraction", "Fraction of the emitted photons inside the HOM detection gate (1 = no gate; set by the gate start and width)", "1",
   _ALL)
_o("F_dot", "Delivered pump Fano factor $F_{\\mathrm{dot}} = 1 - \\eta_{\\mathrm{cap}}(1 - F_{\\mathrm{pump}})$ of the electrical carrier supply "
   "(reported; the counting solvers do not feed it back into $g^{(2)}$)", "1", ("drive",))
_o("linewidth_ueV", "Emitter linewidth $\\hbar(\\Gamma + k + \\Gamma_{\\mathrm{ph}})$", "ueV", _MS)
_o("detuning_ueV", "Emitter minus cavity detuning $\\delta_{\\mathrm{EC}}$", "ueV", _MS)
_o("g2_0", "$g^{(2)}(0)$", "1", _ALL + ("setup.laser_leak", "setup.dark_cps"), ("intrinsic", "raw", "bg_corrected"),
   by_conv={"intrinsic": "g2_0", "raw": "g2_raw", "bg_corrected": "g2_bgc"})
_o("g2_raw", "$g^{(2)}(0)$ with laser leak and dark counts", "1", _CH, key="g2_raw")
_o("hom_I", "HOM overlap $I_{\\mathrm{HOM}}$", "1", _ALL + ("setup",),
   ("intrinsic", "raw", "corrected", "corrected_santori2002", "corrected_ollivier2021"),
   by_conv={"intrinsic": "hom_I", "raw": "hom_V_raw", "corrected": "hom_V_corr",
            "corrected_santori2002": "hom_V_corr_santori2002", "corrected_ollivier2021": "hom_V_corr_ollivier2021"})
_o("hom_V", "HOM visibility $V_{\\mathrm{HOM}}$ (raw|corrected)", "1", _CH,
   ("raw", "corrected", "corrected_santori2002", "corrected_ollivier2021"),
   by_conv={"raw": "hom_V_raw", "corrected": "hom_V_corr",
            "corrected_santori2002": "hom_V_corr_santori2002", "corrected_ollivier2021": "hom_V_corr_ollivier2021"})
_o("hom_sd_factor", "Spectral-diffusion factor $F_{\\mathrm{sd}}$ of the HOM overlap at the photon separation setup.hom_delay_ns (1 = off)",
   "1", _CH + ("setup.hom_delay_ns", "material.sigma_sd_ueV", "material.tau_sd_ns"))
_o("hom_V_raw", "HOM visibility, raw", "1", _CH, key="hom_V_raw")
_o("hom_V_corr", "HOM visibility, $g^{(2)}$-corrected", "1", _CH, key="hom_V_corr")
_o("hom_V_corr_santori2002", "HOM visibility corrected by the Santori-2002-type formula (R, T, eps_int, g2)", "1", _CH,
   key="hom_V_corr_santori2002")
_o("hom_V_corr_ollivier2021", "HOM visibility corrected by the Ollivier-2021 formula (V + g2)/(1 - g2)", "1", _CH,
   key="hom_V_corr_ollivier2021")
_o("B_first_lens", "First-lens brightness $B_{\\mathrm{lens}}$", "1", _CH, key="B_lens")
_o("B_lens", "First-lens brightness $B_{\\mathrm{lens}}$ (alias of B_first_lens)", "1", _CH, key="B_lens")
_o("B_fibre", "Brightness into the fibre $B_{\\mathrm{fibre}}$", "1", _CH)
_o("rate_detected_Hz", "Detected count rate $R_{\\mathrm{det}}$", "Hz", _CH)
_o("B_first_lens_Hz", "Photon rate $\\dot B_{\\mathrm{lens}}$ into the first lens (cw: flux x chain; pulsed: $f_{\\mathrm{rep}}B_{\\mathrm{lens}}$)", "Hz", _CH)
_o("B_fibre_Hz", "Photon rate $\\dot B_{\\mathrm{fibre}}$ into the fibre (before $\\eta_{\\mathrm{optics}}$, $\\eta_{\\mathrm{det}}$)", "Hz", _CH)


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
