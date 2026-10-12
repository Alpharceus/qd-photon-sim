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
eta_prep (FX2b) is the excitation probability BEFORE decay: P_X(end of the pulse window) + photons emitted
inside the window (+ <a+a>(end) for jc; P_XX + XX photons for TPE), from the lindblad / jc / polaron solvers;
the old end-of-window population is `P_X_end` (informational). Incoherent drives report p_loaded.
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

HOM: `hom_I` = single-photon wavepacket overlap (intrinsic, g2-free) x phonon factor.
hom_V_raw = ((1 - g2_raw) I - g2_raw) (1 - eps_int) [E: Ollivier et al., PRL 126, 063602 (2021);
Santori 2002 form], eps_int from a setup param "eps_int" if present else 0 [A]; it uses g2_raw, the g2 a
measurement of the same stream would see. hom_V_corr = (V_raw + g2)/(1 - g2) (Tomm: 0.925, 0.021 ->
0.966). Compare raw to raw and corrected to corrected, never across.

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


def hom_v_raw(I, g2, eps_int=0.0):
    """Predicted raw visibility V_raw = ((1 - g2) I - g2)(1 - eps_int)."""
    return ((1.0 - g2) * I - g2) * (1.0 - eps_int)


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
    if "eta_prep" in res.extra:
        # FX2b: excitation probability before decay = P_X(end) + photons emitted inside the window (+ <a+a> for jc)
        d["eta_prep"] = float(res.extra["eta_prep"])
    elif res.inversion is None and "p_loaded" in res.extra:
        d["eta_prep"] = float(res.extra["p_loaded"])
    else:
        d["eta_prep"] = d["inversion"]
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
    I1 = float(res.I_single)
    d["hom_I"] = I1 * q * q if math.isfinite(I1) else NAN
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
        d["hom_V_raw"] = hom_v_raw(d["hom_I"], d["g2_raw"], eps)
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
_o("eta_prep", "State-preparation probability", "1", _ALL)
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
