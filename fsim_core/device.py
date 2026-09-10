"""Device assembly: one editable description of a complete device, evaluated
through the validated F-series chain (Modules A-E).

The device is NOT an arbitrary network -- the physics defines one signal chain
(dot -> barriers/retention -> cavity -> filter -> detector) with drive and
thermal attached. DeviceDesign holds the configurable blocks of that chain;
evaluate() runs the whole thing and returns curves + scalars.

Headless by construction (three-layer rule): fsim_gui/designer.py is a thin
client over this module. Designs round-trip to YAML in cards/ so every GUI
session stays reproducible from its file.

HONESTY: a design evaluated at point values yields conditional numbers, not
predictions -- unmeasured inputs keep their [A] tags and every result carries
the widest tag of its chain. The envelope treatment lives in run_phase3.
"""
from __future__ import annotations

import copy
import itertools
import json
from dataclasses import asdict, dataclass, field
from pathlib import Path

import numpy as np
import yaml
from scipy.integrate import quad

from .cavity import purcell_eff, tracking_detuning
from .integrator import g2_from, retention
from .drive_mech import mech_from_card, mech_set, set_feasibility
from .loading import (
    BackgroundChannel,
    E_SI,
    KB_SI,
    aperture_g2,
    b_injection,
    f1b_g2,
    f8_g2,
    f8b_thin_fano,
    gamma_eff,
    island_radius_nm,
    loading_probs,
    n_window_competitors,
)
from .qd_gf import PhononParams, ibm_purcell_transmission, ibm_transmission
from .spectral import SpectralResult, epsilon, epsilon2, gamma_of_T
from .thermal import Layer, Stack, t_junction
from .linewidth import LinewidthParams, gamma_anchor
from . import cw_g2, dot_levels, materials, pulse_counting, transport, waveguide
from . import nitride_levels, nitride_transport, nitride_cavity, nitride_materials
try:  # Piece 3 is deliberately an optional sibling during staged integration.
    from .nitride_stark import resolve_bias
except ImportError:  # pragma: no cover - exercised only before piece 3 lands.
    def resolve_bias(diode, *, T_j_K, current_uA=None, junction_voltage_V=None,
                     field_polarity=1, external_field_kVcm=0.0):
        """Staged copy of piece 3's resolver; removed from use when available."""
        out = {"T_j_K": T_j_K, "current_uA": current_uA, "V_j": junction_voltage_V,
               "V_terminal": np.nan, "diode_field_kVcm": np.nan,
               "applied_field_kVcm": np.nan, "bias_valid": False, "invalid_reasons": []}
        if ((current_uA is None) == (junction_voltage_V is None) or
                isinstance(field_polarity, bool) or field_polarity not in (-1, 1) or
                not isinstance(T_j_K, (int, float)) or isinstance(T_j_K, bool) or not np.isfinite(T_j_K) or T_j_K <= 0):
            out["invalid_reasons"].append("invalid bias controls")
            return out
        try:
            if current_uA is not None:
                terminal, vj = diode.v_of_i(float(current_uA) * 1e-6, T_j_K)
                out.update(V_j=vj, V_terminal=terminal)
            else:
                current_A = diode.j_of_vj(float(junction_voltage_V), T_j_K) * diode.area_cm2
                out.update(current_uA=current_A * 1e6, V_terminal=float(junction_voltage_V) + current_A * diode.R_s_ohm)
            dep = diode.depletion(out["V_j"], T_j_K)
            out.update(diode_field_kVcm=dep.F_kVcm,
                       applied_field_kVcm=field_polarity * dep.F_kVcm + external_field_kVcm,
                       bias_valid=all(np.isfinite(out[k]) for k in ("current_uA", "V_j", "V_terminal")))
        except (ArithmeticError, ValueError, OverflowError) as exc:
            out["invalid_reasons"].append(str(exc))
        return out

_ROOT = Path(__file__).resolve().parents[1]


def class_proxy_params() -> dict:
    """Arsenide-class Gamma/retention parameters from the V-a joint fit [A as
    proxy for unmeasured platforms]."""
    p = json.loads((_ROOT / "out" / "phase0" / "fit_params.json").read_text())["params"]
    return {k: p[k] for k in ("gamma0", "a_ac", "b_lo", "E_lo", "r_xx",
                              "a_esc", "E_a", "b_p", "E_b", "b0", "beta")}


@dataclass
class DotBlock:
    delta_xx: float = 3.5        # meV [A on InP/GaAsP]
    gamma_scale: float = 1.0     # multiplies the class Gamma(T) proxy [A]
    r_xx: float = 0.72           # Gamma_XX/Gamma_X [DR class]
    lineshape: str = "lorentzian"  # T2 tier: "lorentzian" (legacy, bit-
                                   # identical) | "ibm" (polaron ZPL+sidebands;
                                   # the Lorentzian-optimism correction)
    phonon: dict = field(default_factory=dict)  # qd_gf.PhononParams overrides
                                   # (geometry l_xy_nm/l_z_nm, materials) [A/DR]
    linewidth: str = "class"    # "class" proxy | "anchored" [A; Matsuda, PRB 2001]
    gamma0: float = 0.25         # meV [V] Bommer et al., JAP 110, 063108 (2011)
    a_ac: float = 2.0e-3         # meV/K [V] Ortner et al., PRB 70, 201301 (2004)
    E_LO: float = 43.0           # meV [V] Ioffe NSM, InP LO phonon
    gamma300: float = 12.0       # meV [A] Matsuda et al., PRB 63, 121304 (2001) class anchor


@dataclass
class RetentionBlock:
    # class-proxy Arrhenius; override any field [A]
    a_esc: float = 0.0
    E_a: float = 0.0             # meV; 0 -> use class proxy
    b_p: float = 0.0
    E_b: float = 35.0
    b0: float = 0.0
    beta: float = 0.0
    mode: str = "proxy"         # "proxy" keeps class fit; "confinement" derives levels [E]
    preset: str = ""
    system: dict = field(default_factory=dict)  # inline DotSystem stack (mutually
                                  # exclusive with preset); never carries overrides
    tau_rad_ns: float = 1.0      # ns [E] dot_levels.retention_params convention
    channel: str = "pair_half"  # Gelinas et al., arXiv:0910.0480 [DR]
    overrides: dict = field(default_factory=dict)  # confinement mode only: explicit
                                  # a_esc/E_a/b_p/E_b/b0/beta replacements for the
                                  # retention_params()-derived values [A]; separate
                                  # from `system` so a preset stack can still be
                                  # combined with explicit overrides. An explicit
                                  # 0.0 here always wins (never a truthy-or fallback).
    n_dot_cm2: float | None = None  # cm^-2; peer-review-triage.md finding 1b -- the
                                  # confinement escape prefactor (a_esc, via
                                  # dot_levels.retention_params' states_per_dot =
                                  # n2d_cm2 / n_dot_cm2) otherwise silently used
                                  # retention_params' own 1e10 cm^-2 default instead
                                  # of the SAME dot density device.py's transport
                                  # call already uses (aperture.density_cm2 /
                                  # drive.n_dot_cm2). None (default): _confinement_params
                                  # falls back to the caller-supplied n_dot_cm2 (that
                                  # shared density), then to retention_params' own
                                  # 1e10 default -- bit-identical to today when a card
                                  # sets neither this nor aperture.density_cm2 [DR].
    tau_cap_ps: float | None = None  # ps; None -> retention_params' own 10.0 ps [E]
                                  # capture-time class default (dot_levels.py).
    tau_cap_scales_with_density: bool = False  # [A] nu_esc0 = (1/tau_cap) * N2D/N_dot
                                  # (dot_levels.retention_params): a sparser n_dot_cm2
                                  # already raises a_esc through N2D/N_dot alone (the
                                  # default, "no-cancellation" convention). True opts
                                  # into the alternative "full-cancellation" convention
                                  # -- tau_cap scaled by (1e10 / n_eff) so per-dot
                                  # capture slows in step with the sparser density and
                                  # the two effects cancel out of a_esc. False (default)
                                  # is legacy/bit-identical whenever n_dot_cm2 is also
                                  # unset.


@dataclass
class DriveBlock:
    V: float = 1.9               # V
    I_uA: float = 10.0
    duty: float = 1.0
    mu: float = 0.5              # cap-2 loading per pulse [A]
    b_e: float = 0.02            # total injection background, signal units [A]
    b_e_m: float = 1.5           # current exponent of b_e [A]
    b_e_Eact: float = 100.0      # meV activation of the WL/barrier EL share [A]
    I_ref_uA: float = 10.0
    mode: str = "EL"             # "EL" (electrical: dg_inj + injection background
                                  # both active) | "PL" (optical excitation: both
                                  # disabled) [A]
    dg_inj: float = 0.0          # meV; F5' injection broadening amplitude [A]
    p_inj: float = 1.0           # F5' injection broadening current exponent [A]
    F_p: float = 1.0             # F8 pump Fano factor (1=Poisson) [A]
    eta_capture: float = 1.0     # F8b dot-capture fraction (thinning) [A]
    C_dep_pF: float = 10.0       # pF; F9 depletion capacitance, informational
                                  # (drives granularity_N/island_radius_nm only
                                  # -- never wired into g2) [E]
    mechanism: str = ""           # T1 drive-mechanism library: "" = legacy path
                                  # (bit-identical); else a drive_mech name
                                  # (poisson-rail/quiet-rail/pulsed/set-*/rti)
    mech_params: dict = field(default_factory=dict)  # per-mechanism params [A]
    diode: dict = field(default_factory=dict)  # transport preset plus Diode overrides [E/DR]
    diode_area_um2: float = 0.0  # um^2; council review 2026-09-05 item 3 -- 0 (default)
                                  # means the diode's injection area defaults to the
                                  # APERTURE area pi*(aperture.diameter_um/2)^2 (V_j,
                                  # N_dots and the per-dot share must all use the SAME
                                  # area the current/dots actually occupy); > 0 is an
                                  # explicit [A] override for current crowding (a
                                  # current path narrower than the optical aperture).
                                  # drive.diode's own "area_um2" key (a literal Diode
                                  # field override) still wins over both if given.
    n_dot_cm2: float = 0.0       # cm^-2; 0 uses aperture density [A]
    cw: bool = False              # opt-in CW g2(tau)/IRF diagnostics (docs/
                                  # rt_edge_contract.md); requires mode=
                                  # "EL-transport" (the CW pump rate is
                                  # transport.loading.r_dot, never inferred
                                  # from the dimensionless pulsed mu) [A]
    cw_irf_fwhm_ps: float = 500.0  # detector IRF FWHM, ps; Reischle 2008
                                  # class 0.5 ns [E]
    cw_irf_shape: str = "gaussian"  # cw_g2.IRF_SHAPES
    cw_pump_ratio: float = 1.0    # X->XX secondary pump ratio [A]. Despite
                                  # the "cw_" prefix this field is shared by
                                  # BOTH pump-rate dynamics this module
                                  # drives with cw_g2.generator: drive.cw's
                                  # CW steady state AND drive.finite_pulse's
                                  # per-period moment hierarchy (pulse_
                                  # counting.pulse_g2's own pump_ratio
                                  # argument) -- one X->XX secondary-pump
                                  # assumption for both, not a separate
                                  # pulsed knob (pr-pkg4-fix item 10; kept
                                  # as one field rather than adding
                                  # drive.pulsed_pump_ratio to avoid a
                                  # second [A] the two opt-in diagnostics
                                  # could silently disagree on).
    cw_tau_max_ns: float = 10.0   # CW g2(tau) window half-width, ns [A]
    rep_rate_hz: float = 0.0      # council review 2026-09-05 item 1 -- pulse
                                  # repetition rate, Hz. 0 (default): derive
                                  # from duty/tau_pulse_ns as before (legacy).
                                  # An explicit EL-transport pulsed card MUST
                                  # set drive.diode['tau_pulse_ns'] and either
                                  # duty or this field (evaluate() raises
                                  # otherwise): DriveBlock's own defaults
                                  # (duty=1.0, tau_pulse_ns=1.0) silently read
                                  # as a 1 GHz/100%-duty DC drive [A].
    loading_model: str = "auto"   # peer-review-triage.md finding 2 -- "auto"
                                  # (default, bit-identical): keeps today's
                                  # switch, f8_g2 (moment-matched) whenever
                                  # F_p != 1.0 exactly, else f1b_g2 (capped
                                  # Poisson) -- an infinitesimal Fano change
                                  # at F_p=1 jumps g2 (0.2205 -> 0.4325 at the
                                  # mu=0.8/eps=0.2 reference point) because
                                  # the two are different, equally valid
                                  # cap-2 conventions that are NOT identical
                                  # at finite mu (loading.py module
                                  # docstring) -- they agree only to O(mu^2).
                                  # "capped_poisson" forces f1b_g2 always;
                                  # "moment_matched" forces f8_g2 always
                                  # (continuous across F_p=1 by construction).
                                  # [A] which convention to use where neither
                                  # is independently anchored; validated at
                                  # evaluate() time, ValueError otherwise.
    finite_pulse: bool = False    # peer-review-triage.md finding 1 -- False
                                  # (default, bit-identical): g2_dot above
                                  # comes from a STATIC per-pulse loading
                                  # distribution (f1b_g2/f8_g2), exact only
                                  # in the limit of an instantaneous pulse.
                                  # True opts into fsim_core.pulse_counting.
                                  # pulse_g2, which propagates the exact
                                  # state-resolved factorial-moment hierarchy
                                  # through the pump/dark windows of one
                                  # pulse period -- at the RT edge-emitter
                                  # operating point the escape rates k_X/k_XX
                                  # (cw_g2.escape_rates_from_retention) and
                                  # the ~100 ps pump are fast relative to the
                                  # radiative lifetime, so a dot can be
                                  # re-excited and re-emit repeatedly within
                                  # one pulse (Hanschke et al., npj Quantum
                                  # Inf. 4, 43 (2018)). Requires
                                  # EL-transport (injection is not None) for
                                  # the physical pump rate
                                  # injection.loading.r_dot; a legacy PL/
                                  # EL-fixed-mu card never touches this path.
                                  # [A] rectangular pump waveform; [DR]
                                  # moment hierarchy on the existing cw_g2
                                  # generator (see pulse_counting.py).
                                  # pr-pkg4-fix item 1 (gate-consistent
                                  # counting, Opus review): whenever this is
                                  # True (and evaluable), rho is ALSO always
                                  # recomputed from pulse_counting's own
                                  # gate-restricted counts (see gate_ns
                                  # below) -- the legacy static-loading rho
                                  # (G*S/(G*S+B) above) is never mixed with
                                  # the dynamic g2_dot; a prior round did
                                  # exactly that mixing and measured a 28.8%
                                  # internal inconsistency at the 230 K
                                  # corner. An unconverged periodic steady
                                  # state (pulse_counting's own `converged`
                                  # flag) sets g2_dot and rho to nan and
                                  # marks the row invalid rather than
                                  # reporting a number from a bad fixed
                                  # point (exposed as finite_pulse_converged
                                  # on the evaluation dict).
    gate_ns: float | None = None  # peer-review-triage.md finding 4 -- the
                                  # width (ns) of the one explicit counting
                                  # gate finite_pulse's rho is built from,
                                  # measured from the PUMP ONSET (module
                                  # pulse_counting.py's own convention).
                                  # None (default): the gate is the WHOLE
                                  # pulse period (tau_pulse_ns + the dark
                                  # window) -- bit-identical to gate_ns
                                  # equal to that period, never the legacy
                                  # per-collected-X/CW-time-model mix (that
                                  # mixing was pr-pkg4-fix item 1's bug: see
                                  # finite_pulse above). n_X+n_XX comes from
                                  # pulse_counting's own gate-restricted
                                  # mean_counts (already includes escape and
                                  # the XX line); gate_ns restricts the
                                  # SIGNAL counting window only (pr-pkg4-fix3
                                  # item 1, Opus review: an earlier version
                                  # of this comment falsely claimed n_bg was
                                  # the same injection background rate the
                                  # CW path uses, integrated over that same
                                  # gate width). The background is
                                  # integrated separately, over min(gate_ns,
                                  # tau_pulse_ns) only: background photons
                                  # are counted while injection current
                                  # actually flows, and any afterglow past
                                  # the pulse end is neglected [A] (see
                                  # finite_pulse above and
                                  # pulse_counting.py's module docstring).
                                  # signal = G*(n_X+n_XX); B_fp = b0 +
                                  # beta*(1-S) + n_bg + b_res*n_X, WITHOUT
                                  # the legacy pulsed rho's own G*S factor on
                                  # the background terms (another earlier
                                  # false claim in this comment, that the
                                  # formula used "the legacy pulsed rho's
                                  # own G and b0/beta placement"): n_bg and
                                  # b_res*n_X are already absolute counts on
                                  # the same footing as G*mean_counts, so
                                  # only the signal term carries G; rho =
                                  # signal/(signal+B_fp). [DR] counts, [A]
                                  # gate width and the afterglow-neglect
                                  # above -- recommended gate_ns =
                                  # tau_pulse_ns + 5*d.ret.tau_rad_ns
                                  # (catches >99% of a single-exponential
                                  # decay tail). Neglecting background
                                  # afterglow is one-sided and optimistic:
                                  # rho_pulsed is an UPPER BOUND -- at the
                                  # FAVOURABLE corner (T_hs=230 K,
                                  # gamma300=6, delta_xx=8, NA=0.8,
                                  # R_back=0.95, L=250 um) n_bg is ~17% of
                                  # B_fp, and a ~1 ns background carrier
                                  # lifetime (~11x more background photons
                                  # there) would move rho_pulsed from 0.858
                                  # to ~0.69, while g2_op moves only
                                  # 0.9817 -> 0.98826 [A]. At T_hs=300 K
                                  # instead (same corner otherwise) n_bg/
                                  # B_fp = 0.117152, rho_pulsed = 0.8662227488,
                                  # g2_op = 0.9976578221 -- quote both
                                  # temperatures, they are not
                                  # interchangeable. Exposed as
                                  # finite_pulse_gate_ns_used
                                  # on the evaluation dict (the resolved
                                  # value, whichever of the two applied,
                                  # NaN if finite_pulse was requested but
                                  # the operating point never converged --
                                  # see pr-pkg4-fix3 item 2). Ignored unless
                                  # finite_pulse is also True.
    b_res: float = 0.0            # council review 2026-09-05 item 8 -- residual
                                  # background channel, T-independent:
                                  # background photons in the collection
                                  # window per COLLECTED X photon (same
                                  # normalization as the item-2 fixed b_e/t_X
                                  # ratio). 0.0 (default): legacy, no residual
                                  # channel. Models whatever the injection
                                  # background model does not capture (e.g.
                                  # the 80 K electrical anchor, Reischle et
                                  # al., Optics Express 16, 12771 (2008)
                                  # (DOI 10.1364/OE.16.012771), measured
                                  # rho ~ 0.88 -> b_res = 1/0.88-1 = 0.136);
                                  # a card setting this NON-zero must
                                  # carry its own provenance note for the
                                  # anchor [E/A] -- device.py does not invent
                                  # a value here, only the mechanism.
    # Nitride-only full-cycle loading controls.  They are deliberately
    # independent of loading_model/F8, which remains a legacy Poisson model.
    cycle_loading: str = "rectangular"
    eta_load: float = 1.0
    set_params: dict = field(default_factory=dict)


@dataclass
class ThermalBlock:
    mesa_diameter_um: float = 1.0
    T_hs: float = 77.0
    layers: list = field(default_factory=lambda: [
        {"name": "epi/cladding", "t_um": 1.5, "k300": 15.0, "alpha": 0.5, "spread": False},
    ])
    substrate: dict = field(default_factory=lambda: {"name": "GaAs", "k300": 55.0, "alpha": 1.25})


@dataclass
class CavityBlock:
    enabled: bool = False
    type: str = "planar"         # "planar" (F6 resonant filter/gain) | "sin_waveguide"
    kappa: float = 1.0           # meV [A -> MEEP/COMSOL]; unused when type=sin_waveguide
    T_track: float = 120.0       # tracking-rule target (K)
    E_X0: float = 1.88           # eV cryogenic X
    dEdT_cav: float = -0.04      # meV/K
    F_P: float = 10.0            # [A]; unused when type=sin_waveguide
    G: float = 8.0               # collection gain on signal [A]; unused when type=sin_waveguide
    beta_sin: float = 1.0        # SiN evanescent coupling [A]; Lemma 1: brightness ONLY,
                                  # never eps/rho/g2/T_c (applied inline in evaluate())
    purcell_wire: bool = False   # R2 (roadmap 13), opt-in: wire F_eff into
                                  # (a) retention (radiative rate x rate_mult
                                  # -> escape ratios divided) and (b) the IBM
                                  # ZPL redistribution (Z_eff funneling).
                                  # [E->analytic-1D, confirm: SIM-B]. False =
                                  # legacy bit-identical (T_c conservative).


@dataclass
class FilterBlock:
    enabled: bool = True
    auto_w: bool = True          # w = Gamma(T_j) operating convention [A]
    auto_w_scale: float = 1.0    # council review 2026-09-05 item 9: multiplies
                                 #   the auto_w window (1.0 = legacy). auto_w
                                 #   sets the collection window to the FULL
                                 #   linewidth Gamma(T_j) -- the F-series
                                 #   operating convention, not a bug -- but a
                                 #   sweep must be able to vary how wide that
                                 #   window is relative to the line without
                                 #   abandoning the auto_w convention entirely
                                 #   (auto_w=False, an explicit fixed w) [A].
    w: float = 2.0               # meV, used when auto_w = False
    dx: float = 0.0              # X offset from window center (meV)
    track: str = "mode"          # "mode": slit follows the cavity mode (legacy,
                                 #   bit-identical) | "hold": F6 slit-held --
                                 #   slit centered ON X at each operating T
                                 #   (lab monochromator convention) while the
                                 #   cavity mode walks per dn/dT. The two
                                 #   coincide at T_track by construction.
    track_material: str = ""     # "" = legacy hard-coded GaAs Varshni tracking
                                 #   (cavity.tracking_detuning, bit-identical);
                                 #   "dot" | "matrix" opts into stack-based
                                 #   tracking via materials.bandgap() on the
                                 #   named layer of the shared inline stack
                                 #   (ret.preset / ret.system) -- never GaAs
                                 #   by default for a non-GaAs stack [DR].  With
                                 #   cavity.enabled=False (no cavity mode to net
                                 #   the walk against) a TRACKING filter is
                                 #   centred ON the dot at every T: dx = filter.dx
                                 #   (0 by default), never the material's own
                                 #   Varshni walk -- see hold_window below
                                 #   [DR, council review 2026-09-05 item 4].
    hold_window: bool = False    # cavity-less only: explicit opt-in for a filter
                                 #   window HELD FIXED at its cavity.T_track
                                 #   position (e.g. a physically fixed
                                 #   monochromator slit calibrated at T_track)
                                 #   while the named track_material layer's own
                                 #   Varshni walk moves the line underneath it --
                                 #   dx = 1e3*(bandgap(mat,Tj)-bandgap(mat,T_track)).
                                 #   Deliberately NOT named "tracking": a filter
                                 #   that tracks the line follows it (dx=0
                                 #   above); this is the opposite, explicitly
                                 #   named case [A, council review 2026-09-05
                                 #   item 4 -- round 1 hard-coded this dx as the
                                 #   ONLY cavity-less track_material behaviour,
                                 #   an artefact that collapsed t_X 0.5 ->
                                 #   0.00106 at the worst sweep corner].


@dataclass
class ApertureBlock:
    density_cm2: float | None = None  # cm^-2; None (default): every consumer
                                 # resolves the legacy 7.0e8 [E] class value
                                 # (pr-pkg1-fix2, item 1) -- EXCEPT
                                 # _confinement_params' own n_dot_cm2, which
                                 # instead falls through to its 1e10 default
                                 # (peer-review-triage.md finding 1b): a
                                 # truthy 7.0e8 default here silently fed
                                 # confinement too, so no legacy design could
                                 # ever reach retention_params' real 1e10
                                 # default. A card that sets this explicitly
                                 # (both edge-emitter cards: 3.0e8) forwards
                                 # that SAME value to confinement and
                                 # transport alike.
    diameter_um: float = 1.0
    sigma_inh: float = 40.0
    comp_brightness: float = 0.3
    compose: bool = False        # False (legacy): aperture_g2_penalty stays an
                                 #   informational scalar (integer-rounded
                                 #   competitor count), never folded into g2.
                                 #   True: continuous [A] Poisson-competitor-
                                 #   bath composition into g2/g2_cw0 (docs/
                                 #   rt_edge_contract.md "Aperture assumptions").


@dataclass
class EmissionBlock:
    type: str = "none"           # "none" (legacy) | "edge" opts into
                                 #   waveguide.edge_emission out-coupling
                                 #   (Lemma 1: brightness only, never g2/eps/rho).
    ridge_width_nm: float = 2000.0
    etch_depth_nm: float = 1200.0
    L_um: float = 500.0
    NA: float = 0.5
    alpha_cm: float = 5.0        # ridge propagation loss, class estimate [E]
    R_back: float | None = None  # None: bare cleaved facet (R_back = R_front)
                                 #   in waveguide.facet_escape_fraction's
                                 #   ray-probability escape-fraction model [A]
    coating: dict = field(default_factory=dict)  # facet_transmission() override [A]
    core_half_nm: float = 148.0  # matrix (well) layer thickness each side of
                                 #   the dot plane, hkust_ridge_stack class [E]
    cladding_nm: float = 1000.0  # barrier/cladding thickness each side [E]
    lambda_nm: float = 0.0       # 0 (default): the actual confinement-derived
                                 #   X transition, dot_levels.levels(system).
                                 #   lambda_nm (contract wording). MATERIAL_
                                 #   EXTRA's optical tables are sparse (a
                                 #   handful of tabulated points per alloy), so
                                 #   most physically-derived InP-dot
                                 #   wavelengths land outside any table and
                                 #   are correctly reported ineligible
                                 #   (missing optics evidence) rather than
                                 #   guessed; an explicit >0 value here is a
                                 #   stated [A] design wavelength (e.g. the
                                 #   668 nm point where the phosphide stack's
                                 #   indices are actually tabulated) used
                                 #   instead of the natural transition.


@dataclass
class DeviceDesign:
    name: str = "my-device"
    dot: DotBlock = field(default_factory=DotBlock)
    ret: RetentionBlock = field(default_factory=RetentionBlock)
    drive: DriveBlock = field(default_factory=DriveBlock)
    thermal: ThermalBlock = field(default_factory=ThermalBlock)
    cavity: CavityBlock = field(default_factory=CavityBlock)
    filter: FilterBlock = field(default_factory=FilterBlock)
    aperture: ApertureBlock = field(default_factory=ApertureBlock)
    emission: EmissionBlock = field(default_factory=EmissionBlock)
    platform: str = "legacy"
    # Planar InGaN/GaN inputs.  These are all explicit card inputs on the
    # nitride tier; defaults are only inert schema defaults [A].
    nitride: dict = field(default_factory=dict)
    provenance: dict = field(default_factory=dict)  # persisted card-side evidence notes

    # ---- YAML round-trip (same reproducibility rule as the cards)
    def save(self, path):
        doc = {"meta": {"name": self.name, "role": "device-design",
                        "note": "editable device description; evaluate with "
                                "fsim_core.device.evaluate"},
               "design": asdict(self)}
        Path(path).write_text(yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")

    @staticmethod
    def load(path) -> "DeviceDesign":
        d = yaml.safe_load(Path(path).read_text(encoding="utf-8"))["design"]
        # item 5 (pr-pkg1-fix2): PyYAML's YAML 1.1 float resolver does not
        # recognize an unsigned exponent (e.g. "3.0e8", vs. "3.0e+8") as a
        # float, so it reads back as a str -- a TypeError the first time
        # arithmetic touches it. Coerce the None-defaulted density/time
        # fields explicitly at load (reproduced TypeError without this).
        ret_kw = _coerce_optional_floats(d["ret"], ("n_dot_cm2", "tau_cap_ps"))
        aperture_kw = _coerce_optional_floats(d["aperture"], ("density_cm2",))
        # finding 4 (peer-review-triage.md): drive.gate_ns is the same
        # None-defaulted float class as ret.n_dot_cm2/aperture.density_cm2
        # above -- guard it against the same YAML 1.1 unsigned-exponent bug.
        drive_kw = _coerce_optional_floats(d["drive"], ("gate_ns",))
        # pr-pkg1-fix3 item 9: ret.n_dot_cm2/tau_cap_ps are None-means-
        # "use the caller-supplied/class default" (_confinement_params'
        # is-None resolution chain, see RetentionBlock's own docstring) --
        # an explicit 0.0 is a real value there, not a second "unset", and
        # divides straight into dot_levels.retention_params' states_per_dot
        # (N2D/n_dot_cm2) and its 1/tau_cap escape-attempt rate, so a card
        # that sets either to exactly 0.0 would silently warn/produce
        # inf/nan deep inside the confinement solve instead of failing at
        # load, where the bad field is still named.
        for key in ("n_dot_cm2", "tau_cap_ps"):
            if ret_kw.get(key) is not None and ret_kw[key] <= 0:
                raise ValueError(f"ret.{key} must be > 0 if set (got {ret_kw[key]!r})")
        # pr-pkg1-fix4 item 5: aperture.density_cm2 is the same None-means-
        # "use the legacy/confinement default" class as ret.n_dot_cm2 above
        # (see ApertureBlock.density_cm2's own docstring) -- an explicit
        # 0.0 divides straight into transport's f_qd and confinement's
        # states_per_dot (N2D/n_dot_cm2) the same way, so it must fail here
        # too, not deep inside evaluate().
        if aperture_kw.get("density_cm2") is not None and aperture_kw["density_cm2"] <= 0:
            raise ValueError("aperture.density_cm2 must be > 0 if set "
                             f"(got {aperture_kw['density_cm2']!r})")
        # pr-pkg4-fix item 9: drive.loading_model and drive.gate_ns are both
        # validated again in evaluate() (below), but a card that misspells
        # either should fail at load(), the same as every other malformed
        # card field above, not only once evaluate() is actually called.
        loading_model = drive_kw.get("loading_model", DriveBlock.loading_model)
        if loading_model not in ("auto", "capped_poisson", "moment_matched"):
            raise ValueError(f"unknown drive.loading_model {loading_model!r}")
        if drive_kw.get("gate_ns") is not None and drive_kw["gate_ns"] <= 0:
            raise ValueError(f"drive.gate_ns must be > 0 if set (got {drive_kw['gate_ns']!r})")
        # Nitride opt-in (spec: "Add DeviceDesign.platform ... reject unknown
        # values"): validate the platform string itself, and require the
        # dot's InP-shaped shared fields (gamma0/a_ac/E_LO/gamma300/r_xx/
        # delta_xx) to be EXPLICITLY present in the card's own YAML mapping
        # -- DotBlock's class defaults are InP class numbers [V]/[A] (e.g.
        # E_LO=43 meV, the InP LO phonon), so a nitride card that leaves any
        # of these unset would silently inherit a wrong-material default
        # with no error. This presence check is only possible here, against
        # the raw YAML dict, before DotBlock's own dataclass defaults fill
        # in; evaluate() cannot re-derive "was this explicit" from an
        # already-built in-memory design, so it only checks finiteness
        # there (see _evaluate_nitride) [A, conservative reading of "require
        # explicit ... on cards"].
        platform = d.get("platform", "legacy")
        if platform not in ("legacy", "ingan_gan_planar"):
            raise ValueError(f"unknown platform {platform!r}")
        if platform == "ingan_gan_planar":
            dot_raw = d.get("dot", {})
            required_dot_fields = ("gamma0", "a_ac", "E_LO", "gamma300", "r_xx", "delta_xx")
            missing = [k for k in required_dot_fields if k not in dot_raw]
            if missing:
                raise ValueError("ingan_gan_planar cards must set dot." +
                                 ", dot.".join(missing) + " explicitly")
        return DeviceDesign(
            name=d.get("name", "my-device"),
            dot=DotBlock(**d["dot"]), ret=RetentionBlock(**ret_kw),
            drive=DriveBlock(**drive_kw), thermal=ThermalBlock(**d["thermal"]),
            cavity=CavityBlock(**d["cavity"]), filter=FilterBlock(**d["filter"]),
            aperture=ApertureBlock(**aperture_kw),
            emission=EmissionBlock(**d.get("emission", {})),
            platform=platform, nitride=dict(d.get("nitride", {})),
            provenance=dict(d.get("provenance", {})),
        )


def _coerce_optional_floats(mapping: dict, keys: tuple) -> dict:
    """Shallow copy of `mapping` with each of `keys` coerced through float()
    when present and not None (item 5, pr-pkg1-fix2): guards the None-
    defaulted density/time card fields against PyYAML's YAML 1.1 float
    resolver reading an unsigned-exponent literal like "3.0e8" back as a
    str."""
    out = dict(mapping)
    for key in keys:
        if out.get(key) is not None:
            out[key] = float(out[key])
    return out


def _stack(th: ThermalBlock) -> Stack:
    return Stack(
        layers=[Layer(L["t_um"] * 1e-6, L["k300"], L.get("alpha", 1.25),
                      spread=bool(L.get("spread", False))) for L in th.layers],
        k_sub300=th.substrate["k300"], sub_alpha=th.substrate.get("alpha", 1.25),
    )


def _material_from_mapping(value):
    """Decode a card material without changing the public material API [DR]."""
    if isinstance(value, str):
        return materials.binary(value)
    if not isinstance(value, dict):
        raise ValueError("ret.system material must be a name or mapping")
    kind = value.get("kind", value.get("alloy", "binary"))
    if kind == "binary":
        return materials.binary(value["name"])
    constructors = {"GaAsP": materials.GaAsP, "GaInP": materials.GaInP,
                    "AlGaAs": materials.AlGaAs, "AlGaInP": materials.AlGaInP,
                    "InGaAs": materials.InGaAs}
    if kind not in constructors:
        raise ValueError(f"ret.system unknown material kind {kind!r}")
    kw = {k: v for k, v in value.items() if k not in ("kind", "alloy", "name", "tag", "source")}
    return constructors[kind](**kw)


def _retention_system(ret: RetentionBlock):
    if ret.preset and ret.system:
        raise ValueError("ret.preset and ret.system cannot both be set")
    if ret.preset:
        presets = dot_levels.class_presets()
        if ret.preset not in presets:
            raise ValueError(f"unknown ret.preset {ret.preset!r}")
        # The actual junction temperature, rather than a preset's cryogenic
        # display temperature, is installed by _confinement_params below.
        return presets[ret.preset]()
    if not ret.system:
        raise ValueError("ret.mode='confinement' requires ret.preset or ret.system")
    raw = ret.system
    geometry = dot_levels.DotGeometry(**raw.get("geometry", {}))
    return dot_levels.DotSystem(
        _material_from_mapping(raw["dot"]), _material_from_mapping(raw["matrix"]),
        _material_from_mapping(raw["barrier"]), _material_from_mapping(raw["substrate"]),
        T=float(raw.get("T", 300.0)), geometry=geometry,
        vbo_override_eV=dict(raw.get("vbo_override_eV", {})),
        eps_r=raw.get("eps_r"), name=raw.get("name", "inline-system"),
    )


def _diode_from_drive(drive: DriveBlock, aperture_diameter_um: float):
    raw = dict(drive.diode)
    preset = raw.pop("preset", "")
    for key in ("tau_pulse_ns", "tau_rad_ns", "w_meV", "dE_WL_meV",
                "eta_rad_matrix", "E_urbach_meV", "eta_total"):
        raw.pop(key, None)
    # Inline material spellings are deliberately decoded here: no card is
    # allowed to silently replace its stated stack with a preset [DR].
    for key in ("p_cladding", "n_cladding", "active", "barrier", "substrate"):
        if key in raw:
            raw[key] = _material_from_mapping(raw[key])
    if "area_um2" not in raw:
        # item 3 (council review 2026-09-05): V_j was solved at J = I / a
        # PRESET mesa area (Diode.area_um2's own 0.785 um^2 default) while
        # N_dots and the per-dot share use the APERTURE area -- 6.25x
        # smaller at the edge-emitter cards' 0.4 um aperture -- so the two
        # halves of the same injection disagreed about where the current
        # actually flows.  Default the diode's injection area to the SAME
        # aperture area used everywhere else; drive.diode_area_um2 > 0 is
        # the explicit [A] override for current crowding (a current path
        # narrower than the optical aperture); an explicit "area_um2" key
        # inside drive.diode itself (a literal Diode field override, same
        # mechanism as R_s_ohm above) still wins over both.
        raw["area_um2"] = (drive.diode_area_um2 if drive.diode_area_um2 > 0
                           else float(np.pi * (aperture_diameter_um / 2.0) ** 2))
    if preset == "hkust":
        return transport.hkust_preset(**raw)
    if preset in ("red", "red_diode"):
        return transport.red_diode_preset(**raw)
    raise ValueError("drive.mode='EL-transport' requires diode.preset 'hkust' or 'red'")


def _legacy_density_cm2(density_cm2: float | None) -> float:
    """item 1 (pr-pkg1-fix2): ApertureBlock.density_cm2's pre-fix literal
    7.0e8 [E] default, applied by every TRANSPORT-side consumer (evaluate_
    injection, the aperture competitor-bath composition) when the card
    leaves density_cm2 unset -- legacy designs must stay bit-identical.
    Confinement's own n_dot_cm2 resolution is deliberately NOT routed
    through this (see _confinement_params / ApertureBlock.density_cm2):
    an unset density there falls to retention_params' real 1e10 default
    instead (peer-review-triage.md finding 1b)."""
    return density_cm2 if density_cm2 is not None else 7.0e8


def _confinement_params(ret: RetentionBlock, Tj: float, n_dot_cm2: float | None = None,
                        tau_cap_ps: float | None = None) -> dict:
    """Resolve confinement at the temperature consumed by escape [DR].  Also
    carries the level table's own E_X_eV (resolved transition energy) for
    the item-2 sub-turn-on loading suppression -- the SAME confinement
    solve, never a second one.

    n_dot_cm2 / tau_cap_ps (peer-review-triage.md finding 1b): the escape
    prefactor a_esc = tau_rad * nu_esc0, nu_esc0 = (1/tau_cap) * N2D/N_dot
    (dot_levels.retention_params), used to silently take retention_params'
    own 1e10 cm^-2 / 10 ps defaults no matter what dot density the SAME
    card's transport call (device.py's aperture.density_cm2 / drive.
    n_dot_cm2) already used. Resolution order: the card's own ret.n_dot_cm2
    / ret.tau_cap_ps win; else the caller-supplied n_dot_cm2/tau_cap_ps (the
    shared aperture/drive density this call site passes in); else
    retention_params' own 1e10 cm^-2 / 10 ps [E] class defaults -- so a
    legacy card (neither set, no aperture.density_cm2 forwarded) is
    bit-identical to before this fix. [DR] coupling algebra (the a_esc
    formula itself, unchanged, lives in dot_levels.retention_params).

    item 4 (pr-pkg1-fix2): T_ref=Tj is now forwarded to retention_params so
    the N2D/N_dot thermal state density is evaluated at the SAME junction
    temperature as the level solve above, not retention_params' own 300 K
    default -- a_esc at the 230 K corner was ~30% too large otherwise."""
    # item 5 (pr-pkg1-fix2): explicit is-None tests, not truthy-or chains, so
    # an explicit 0.0 would win (device.py's own convention, see
    # RetentionBlock.overrides above) -- a truthy-or here would silently
    # replace an explicit 0.0 with the next fallback.
    n_eff = ret.n_dot_cm2 if ret.n_dot_cm2 is not None else (
        n_dot_cm2 if n_dot_cm2 is not None else 1e10)
    tau_cap = ret.tau_cap_ps if ret.tau_cap_ps is not None else (
        tau_cap_ps if tau_cap_ps is not None else 10.0)  # [E] dot_levels.py capture-time class default
    if ret.tau_cap_scales_with_density:
        # [A] full-cancellation convention: nu_esc0 = (1/tau_cap) * N2D/N_dot
        # already makes a sparser n_eff raise a_esc through N2D/N_dot alone
        # (the default, no-cancellation convention); this additionally
        # scales tau_cap by (1e10 / n_eff) so per-dot capture slows in step
        # with the sparser density and the two effects cancel -- a_esc then
        # stops moving with n_eff at all.
        tau_cap = tau_cap * (1e10 / n_eff)
    system = _retention_system(ret)
    # DotSystem is mutable, so copy rather than mutating the card-derived
    # object.  This is also why evaluate leaves its input design untouched.
    system = copy.copy(system)
    system.T = float(Tj)
    lv = dot_levels.levels(system)
    params = dot_levels.retention_params(lv, ret.tau_rad_ns, ret.channel,
                                         T_ref=Tj, n_dot_cm2=n_eff, tau_cap_ps=tau_cap,
                                         verbose=False)
    params["E_X_eV"] = lv.E_X_eV
    params["n_dot_cm2_used"] = n_eff
    return params


def _tracked_material(ret: RetentionBlock, layer: str) -> materials.Material:
    """Resolve filter.track_material ("dot" | "matrix") against the SAME
    shared inline stack (ret.preset / ret.system) used by confinement
    retention -- never an independently selected material [DR]."""
    if layer not in ("dot", "matrix"):
        raise ValueError(f"unknown filter.track_material {layer!r} (use 'dot' or 'matrix')")
    system = _retention_system(ret)
    return system.dot if layer == "dot" else system.matrix


def _aperture_lambda(density_cm2, aperture_um2, w_full_meV, sigma_inh_meV,
                     comp_brightness) -> tuple:
    """Continuous (unrounded) competitor count N_w and composed bath mean
    lambda = N_w * comp_brightness (docs/rt_edge_contract.md "Aperture
    assumptions"). NaN in, NaN out (no meaningful operating window)."""
    if not np.isfinite(w_full_meV):
        return float("nan"), float("nan")
    Nw = float(n_window_competitors(density_cm2, aperture_um2, w_full_meV, sigma_inh_meV))
    return Nw, Nw * comp_brightness


def _compose_aperture_g2(g_target, lam):
    """Continuous aperture composition [A]: the competitor bath is an
    aggregate Poisson process of mean lam = N_w * comp_brightness (N_w
    continuous, never rounded), independent of the target and of any other
    loss already folded into g_target (filter/capture/loading are upstream of
    this call, composed exactly once). Factorial-moment derivation, target
    brightness normalized to 1:

        <n(n-1)> = g_target*1^2  +  lam^2  +  2*1*lam     (target + Poisson
                                                            bath + independent
                                                            cross term)
        <n>^2   = (1 + lam)^2

        g_mix = (g_target + 2 lam + lam^2) / (1 + lam)^2

    Continuous in lam (no integer rounding/clipping/additive penalty); lam=0
    (N_w=0 or comp_brightness=0) returns g_target exactly."""
    lam = np.asarray(lam, dtype=float)
    return (np.asarray(g_target, dtype=float) + 2.0 * lam + lam * lam) / (1.0 + lam) ** 2


def _resolve_edge(ret: RetentionBlock, emission: EmissionBlock, Tj: float):
    """Convert the shared inline stack (ret.preset / ret.system) into
    waveguide.Layer objects at the actual emission wavelength (the
    confinement-derived X transition, dot_levels.levels(system).lambda_nm)
    and run waveguide.edge_emission's documented public API. Never uses the
    hkust_ridge_stack() convenience builder (that stack's own composition is
    a separate literature fixture, not inferred from any card's stack) and
    never touches private waveguide internals."""
    if not (ret.preset or ret.system):
        raise ValueError("emission.type='edge' requires ret.preset or ret.system "
                          "(the shared inline stack; no independently selected preset)")
    system = copy.copy(_retention_system(ret))
    system.T = float(Tj)
    lv = dot_levels.levels(system)
    lambda_nm = emission.lambda_nm if emission.lambda_nm > 0 else lv.lambda_nm
    g = system.geometry

    def n_at(mat):
        return materials.refractive_index(mat.label, lambda_nm)

    # council review 2026-09-05 item 5: each layer carries its materials.py
    # label so edge_emission's group-index finite difference can re-resolve
    # it at a shifted wavelength (generalizes beyond the hkust_ridge_stack
    # special case that used to be the only dispersive stack).
    layers = [
        waveguide.Layer("barrier_lower", n_at(system.barrier), emission.cladding_nm,
                        material=system.barrier.label),
        waveguide.Layer("matrix_lower", n_at(system.matrix), emission.core_half_nm,
                        material=system.matrix.label),
        waveguide.Layer("dot", n_at(system.dot), max(g.height_nm, 0.1), True,
                        material=system.dot.label),
        waveguide.Layer("matrix_upper", n_at(system.matrix), emission.core_half_nm,
                        material=system.matrix.label),
        waveguide.Layer("barrier_upper", n_at(system.barrier), emission.cladding_nm,
                        material=system.barrier.label),
    ]
    edge = waveguide.edge_emission(
        layers, emission.ridge_width_nm, emission.etch_depth_nm, lambda_nm,
        emission.L_um, emission.NA, alpha_cm=emission.alpha_cm,
        R_back=emission.R_back, coating=(emission.coating or None))
    return edge, lambda_nm


def _transport_options(drive: DriveBlock, w_meV: float) -> dict:
    """Card-only transport knobs; units are ns and meV [A/E]."""
    raw = drive.diode
    return dict(tau_pulse_ns=float(raw.get("tau_pulse_ns", 1.0)),
                tau_rad_ns=float(raw.get("tau_rad_ns", 1.0)),
                w_meV=float(raw.get("w_meV", w_meV)),
                dE_WL_meV=float(raw.get("dE_WL_meV", 100.0)),
                eta_rad_matrix=float(raw.get("eta_rad_matrix", 0.1)),
                E_urbach_meV=raw.get("E_urbach_meV"),
                eta_total=float(raw.get("eta_total", 0.01)))


def _nitride_flat_background_acceptance(kappa_meV, w_meV, dx_w_meV, detuning_meV,
                                        reservoir_offset_meV=0.0):
    """Average (flat-spectrum) cavity+slit transmission across the SAME
    collection window the X line is filtered through.

    transport.xi_window's own window is "centred dE_WL below the [matrix/WL]
    peak" -- i.e. it IS the X-line collection window itself (dE_WL_meV is
    only the energy distance used inside the Urbach exponential, not a
    second window position). A flat/broadband continuum crossing that same
    window is filtered by the cavity+slit's spectral acceptance function
    averaged (not Lorentzian-line-weighted like spectral.epsilon2's t_X)
    across the reservoir-centred window -- spec: "not t_X" [DR extension of spectral.py's
    transmission2/combined_transmission kernel, flat-weighted instead of
    line-shape-weighted].
    """
    half = w_meV / 2.0
    off = dx_w_meV - detuning_meV  # cavity center, in slit-relative meV (spectral.transmission2 convention)
    val, _ = quad(lambda x: (kappa_meV / 2.0) ** 2 / ((x - off) ** 2 + (kappa_meV / 2.0) ** 2),
                 reservoir_offset_meV-half, reservoir_offset_meV+half, limit=200)
    return val / w_meV


def _nitride_reservoir_energy_eV(dot_kw, T_K, background):
    """Reservoir emission energy independent of the confined dot line.

    An explicit ``reservoir_energy_eV`` is an opt-in card override [A].
    Otherwise a nonzero InGaN wetting layer uses its strained material
    continuum edge from nitride_materials.band_edges; a zero-thickness WL
    uses the GaN barrier Varshni edge.  A 25 meV localization/Urbach
    downshift [A] is applied to the derived continuum, a conservative class
    offset rather than a dot-QCSE-dependent line.  band_edges carries the
    Rinke et al., PRB 77, 075202 (2008) strain and Tsai & Bayram, ACS Omega
    5, 3917 (2020) band-edge provenance [V].

    Documented disagreement (Opus re-review, 2026-09-09, low finding):
    this function reads the material continuum from nitride.dot's OWN
    wl_thickness_nm (0.0 on the shipped headline cards -> GaN barrier
    edge, ~3.41 eV at 300 K), while the diode's SRH background
    (nitride_transport.evaluate_injection) lives in drive.diode's
    separate 0.5-nm InGaN layer (edge ~2.37-3.02 eV over the swept x_in
    range). These two "reservoir" concepts are independently configured
    and currently disagree about which material hosts the background
    carriers; there is no numerical consequence today because this
    function's return value feeds only the flat-spectrum cavity/slit
    ACCEPTANCE of the SRH rate (a spectral filter), never the SRH rate
    itself. Left as a known modeling inconsistency rather than silently
    reconciled -- see out/nitride_cavity/results.md's Limitations section.
    """
    if "reservoir_energy_eV" in background:
        value = float(background["reservoir_energy_eV"])
        if not np.isfinite(value) or value <= 0:
            raise ValueError("invalid nitride.background.reservoir_energy_eV")
        return value
    wl_nm = float(dot_kw.get("wl_thickness_nm", 0.0))
    gan = nitride_materials.binary("GaN")
    if wl_nm > 0:
        mat = nitride_materials.ingaN(float(dot_kw["x_in"]))
        edges = nitride_materials.band_edges(mat, T_K, substrate=gan,
            strain_fraction=float(dot_kw.get("strain_fraction", 1.0)),
            strain_c_fraction=float(dot_kw.get("strain_c_fraction", .7)))
    else:
        edges = nitride_materials.band_edges(gan, T_K, substrate=gan)
    return float(edges["Ec_eV"] - edges["Ev_eV"] - .025)  # [A] 25 meV localization/Urbach offset


def _evaluate_nitride(d: DeviceDesign, T_grid=None) -> dict:
    """Opt-in planar InGaN/GaN cycle evaluator.

    Bias field [Opus fix-round, 2026-09-09]: the field passed to
    nitride_levels as ``external_field_kVcm`` is the diode's OWN
    depletion(vj, Tj).F_kVcm -- the physical p-i-n depletion field, which
    shrinks toward 0 as V_j -> V_bi (forward bias reduces the depletion
    field; never a lumped -vj/d_i estimate, which ignores V_bi entirely and
    grows in the wrong direction with forward bias). Sign convention: this
    diode field is added DIRECTLY (same sign, no extra negation) to
    nitride_levels' own intrinsic polarization field, which follows the
    +c-positive convention documented in nitride_materials.py's module
    docstring (fsim_core/nitride_materials.py `polarization_field`). Whether
    that addition reinforces or partially cancels the intrinsic QCSE field
    depends on the (unmodeled) relative orientation of the p-i-n stack's
    growth axis vs. the dot's own polarization axis -- this direct-addition
    convention is the documented [A] lumped choice; the rejected
    alternative (always negating depletion().F_kVcm to force screening) is
    recorded under this spec's fix-round decisions. It is intentionally
    separate from the legacy cubic-material transport path.
    """
    n = dict(d.nitride)
    allowed = {"dot", "cavity", "k_nr_ns", "background_tau_ns", "eta_background", "tau_rad0_ns", "background", "bias"}
    extra = set(n) - allowed
    if extra: raise ValueError("unknown nitride keys: " + ", ".join(sorted(extra)))
    if d.drive.mode != "EL-transport" or d.drive.diode.get("preset") != "nitride-planar":
        raise ValueError("ingan_gan_planar requires EL-transport and diode.preset='nitride-planar'")
    if not d.drive.finite_pulse or d.drive.mechanism:
        raise ValueError("ingan_gan_planar requires finite_pulse=True and empty mechanism")
    if d.drive.cycle_loading not in ("rectangular", "deterministic_pair"):
        raise ValueError("unknown cycle_loading")
    if d.dot.linewidth != "anchored" or d.dot.lineshape != "lorentzian" or d.ret.mode != "nitride_confinement":
        raise ValueError("nitride requires anchored lorentzian and ret.mode='nitride_confinement'")
    if d.emission.type != "vertical_cavity" or not d.cavity.enabled or d.cavity.type != "nitride_planar":
        raise ValueError("nitride requires enabled nitride_planar vertical_cavity")
    if d.ret.system or d.ret.preset: raise ValueError("nitride does not use legacy ret.system/preset")
    # Legacy cavity.F_P/G would ambiguously compete with nitride.cavity's own
    # Q/mode-volume/eta_out model; reject a card that overrides them away
    # from their inert class defaults instead of silently ignoring them.
    if d.cavity.F_P != CavityBlock.F_P or d.cavity.G != CavityBlock.G:
        raise ValueError("nitride does not use legacy cavity.F_P/G; nitride.cavity owns Q/V/collection")
    dot_kw = dict(n.get("dot", {})); cav_kw = dict(n.get("cavity", {}))
    required_dot = {"height_nm", "radius_nm", "x_in"}
    if not required_dot <= set(dot_kw): raise ValueError("nitride.dot requires height_nm, radius_nm, x_in")
    # PyYAML's YAML 1.1 float resolver does not recognize an unsigned
    # exponent (e.g. "3.0e8" vs "3.0e+8") as a float, so a card written that
    # way reads a dot_kw value back as a str. NitrideDotSystem itself has no
    # __post_init__ validation (nitride_levels._validate only runs later, at
    # evaluation time), so an unvalidated str would otherwise reach
    # math.isfinite() deep inside nitride_levels.levels() and raise
    # TypeError instead of this module's ValueError contract (Opus fix-round
    # finding: "malformed nitride inputs escape the ValueError card-schema
    # contract"). Reject it here, at construction, with the offending key
    # named.
    # Piece 2 adds string-valued geometry/orientation selectors and an
    # optional polarization factor.  Numeric leaves remain finite scalars;
    # NitrideDotSystem validates selector combinations at the solver boundary.
    _dot_numeric = {"height_nm", "radius_nm", "x_in", "wl_thickness_nm",
                    "strain_fraction", "screening_fraction", "external_field_kVcm",
                    "vbo_InN_GaN_eV", "strain_c_fraction", "top_radius_fraction"}
    for _k, _v in dot_kw.items():
        if _k in _dot_numeric and (isinstance(_v, bool) or not isinstance(_v, (int, float)) or not np.isfinite(_v)):
            raise ValueError(f"invalid nitride.dot.{_k}: must be a finite number (got {_v!r})")
    try:
        system0 = nitride_levels.NitrideDotSystem(**dot_kw)
    except TypeError as exc:
        raise ValueError("invalid nitride.dot field(s): " + str(exc)) from exc
    try:
        cav = nitride_cavity.NitrideCavityParams(**cav_kw)
    except TypeError as exc:
        raise ValueError("invalid nitride.cavity field(s): " + str(exc)) from exc
    bias_kw = dict(n.get("bias", {}))
    if set(bias_kw) - {"mode", "field_polarity", "V_j_V", "T_j_K", "cavity_reference_V_j_V"}:
        raise ValueError("unknown nitride.bias keys")
    bias_mode = bias_kw.get("mode", "current")
    if bias_mode not in ("current", "junction_voltage"):
        raise ValueError("nitride.bias.mode must be current or junction_voltage")
    polarity = bias_kw.get("field_polarity", 1)
    if isinstance(polarity, bool) or not isinstance(polarity, int) or polarity not in (-1, 1):
        raise ValueError("nitride.bias.field_polarity must be integer +1 or -1")
    def _bias_number(name):
        value = bias_kw[name]
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not np.isfinite(value):
            raise ValueError("nitride.bias.%s must be finite numeric" % name)
        return float(value)
    if bias_mode == "current":
        if "V_j_V" in bias_kw or "T_j_K" in bias_kw:
            raise ValueError("current bias mode does not accept V_j_V or T_j_K")
    else:
        if "V_j_V" not in bias_kw or "T_j_K" not in bias_kw:
            raise ValueError("junction_voltage bias mode requires V_j_V and T_j_K")
        if _bias_number("V_j_V") < 0 or _bias_number("T_j_K") <= 0:
            raise ValueError("junction_voltage requires V_j_V >= 0 and T_j_K > 0")
    cavity_reference_V = bias_kw.get("cavity_reference_V_j_V")
    if cavity_reference_V is not None and _bias_number("cavity_reference_V_j_V") < 0:
        raise ValueError("nitride.bias.cavity_reference_V_j_V must be >= 0 if set")
    if system0.geometry_type == "qw_fluctuation":
        for _key in ("wl_thickness_nm", "x_in"):
            _value = d.drive.diode.get(_key)
            if isinstance(_value, bool) or not isinstance(_value, (int, float)) or _value != getattr(system0, _key):
                raise ValueError("QW cards require drive.diode.%s == nitride.dot.%s" % (_key, _key))
        if d.drive.diode.get("d_i_nm", 0) < system0.height_nm:
            raise ValueError("QW cards require drive.diode.d_i_nm >= dot height")
        if "reservoir_energy_eV" in dict(n.get("background", {})):
            raise ValueError("QW cards cannot override the resolved reservoir energy")
    tau_on = float(d.drive.diode.get("tau_pulse_ns", 0.0)); rep = float(d.drive.rep_rate_hz)
    if tau_on <= 0 or rep <= 0: raise ValueError("nitride requires positive explicit pulse width and rep_rate_hz")
    period = 1e9 / rep
    if tau_on > period: raise ValueError("pulse width exceeds period")
    # Card-configuration validation (malformed opt-in -> raise unconditionally,
    # same as every other check in this block); deliberately NOT inside the
    # per-T try/except below, which is reserved for physically invalid
    # operating points (unbound dot, thermal runaway), not malformed cards.
    if d.drive.cycle_loading == "deterministic_pair":
        sp0 = dict(d.drive.set_params)
        if not ({"eps_r", "R_T_ohm", "ec_margin"} <= set(sp0)
               and ("radius_nm" in sp0 or "C_sigma_F" in sp0)):
            raise ValueError("SET requires explicit island and feasibility inputs")
    duty = tau_on * rep * 1e-9  # electrical pulse duty cycle (spec worked example: 0.1 ns @ 8e7 Hz -> 0.008)
    density = d.drive.n_dot_cm2 if d.drive.n_dot_cm2 > 0 else d.aperture.density_cm2
    if density is None or density <= 0: raise ValueError("nitride requires explicit positive dot density")
    tau_cap = d.ret.tau_cap_ps
    if tau_cap is None or tau_cap <= 0: raise ValueError("nitride requires explicit positive capture time")
    for key in ("gamma0", "a_ac", "E_LO", "gamma300", "r_xx", "delta_xx"):
        val = getattr(d.dot, key)
        try:
            finite = np.isfinite(val)
        except TypeError as exc:
            raise ValueError(f"invalid dot.{key}: not numeric ({val!r})") from exc
        if not finite: raise ValueError("invalid dot." + key)
    eta_bg = float(n.get("eta_background", .1)); bg_tau = float(n.get("background_tau_ns", 0.0))
    bg_kw = dict(n.get("background", {}))
    if set(bg_kw) - {"reservoir_energy_eV"}:
        raise ValueError("unknown nitride.background keys")
    if not 0 <= eta_bg <= 1 or bg_tau < 0: raise ValueError("invalid nitride background settings")
    if not 0 <= d.drive.eta_load <= 1: raise ValueError("drive.eta_load must be in [0, 1]")
    # cw_pump_ratio (X->XX secondary pump) is forwarded into pulse_counting.
    # pulse_g2 for the rectangular Poisson pump below; deterministic_pair has
    # no continuous pump to apply it to, so a non-default value there is an
    # ambiguous competing override, rejected the same way the other
    # legacy-overlap fields above are -- never silently dropped (Opus
    # fix-round finding: "cw_pump_ratio silently dropped on the nitride
    # branch").
    if d.drive.cycle_loading == "deterministic_pair" and d.drive.cw_pump_ratio != DriveBlock.cw_pump_ratio:
        raise ValueError("drive.cw_pump_ratio is not applicable to deterministic_pair loading")
    tau_rad0_ns = float(n.get("tau_rad0_ns", 1.0))  # [A brief default]; e.g. Deshpande et al.,
    # APL 105, 141109 (2014), DOI 10.1063/1.4897640, abstract 300 K tau=1.3+/-0.3 ns [V abstract-only]
    diode_kw = dict(d.drive.diode); diode_kw.pop("preset", None); diode_kw.pop("tau_pulse_ns", None)
    try:
        diode = nitride_transport.planar_pin(**diode_kw)
    except TypeError as exc:
        raise ValueError("invalid drive.diode field(s): " + str(exc)) from exc
    st, a = _stack(d.thermal), .5*d.thermal.mesa_diameter_um*1e-6
    aperture = float(np.pi*(d.aperture.diameter_um/2)**2)
    if bias_mode == "junction_voltage":
        fixed_T = _bias_number("T_j_K")
        if T_grid is not None:
            requested_ts = np.asarray(T_grid, dtype=float)
            if requested_ts.size != 1 or requested_ts[0] != fixed_T:
                raise ValueError("junction_voltage mode only accepts omitted or singleton T_grid=[T_j_K]")
        ts = np.asarray([fixed_T], dtype=float)
    else:
        ts = np.asarray(T_grid if T_grid is not None else [d.thermal.T_hs], dtype=float)
    # Cavity-tracking reference: the SAME dot, at the SAME supplied current
    # (the design's own bias state), evaluated at T_track instead of the
    # per-row Tj -- a single fixed design anchor, not recomputed per
    # operating-point Tj [A, conservative reading of "same stated
    # bias-control convention"]. Opus fix-round finding 2: the previous
    # version anchored on the UNBIASED dot (system0, external_field_kVcm as
    # given by the card, normally 0) while the operating-point E_X used the
    # bias-resolved dot below, so the two differed by the full Stark shift
    # at every T -- never on-resonance even when T_j == T_track. Resolving
    # vj at T_track through the SAME diode.v_of_i/depletion() path used for
    # the operating point makes the two calls agree on the bias field (and
    # hence ~zero detuning) whenever a row's T_j lands on T_track.
    _track_v = cavity_reference_V
    if _track_v is None and bias_mode == "junction_voltage":
        _track_v = _bias_number("V_j_V")
    _track_bias = resolve_bias(diode, T_j_K=cav.T_track,
                               current_uA=d.drive.I_uA if _track_v is None else None,
                               junction_voltage_V=_track_v, field_polarity=polarity,
                               external_field_kVcm=system0.external_field_kVcm)
    if not _track_bias["bias_valid"]:
        raise ValueError("invalid cavity reference bias: " + "; ".join(_track_bias["invalid_reasons"]))
    _track_dsys = nitride_levels.NitrideDotSystem(**{**dot_kw,
        "external_field_kVcm": _track_bias["applied_field_kVcm"]})
    track0 = nitride_levels.levels(_track_dsys, cav.T_track)
    # The whole-period gate, resolved once (T-independent); used as the
    # reported gate_ns_used on an invalid row too (a requested experimental
    # setting, not a derived physics result).
    requested_gate = min(d.drive.gate_ns if d.drive.gate_ns is not None else period, period)
    rows=[]
    for ths in ts:
        # Electro-thermal self-consistency: diode.junction_power(I_uA, T)
        # depends on T (through v_of_i's own T-dependence), so a one-shot
        # T_j = t_junction(power-at-heat-sink-T, ...) is not self-consistent
        # (Opus fix-round finding: "no self-consistency loop and no
        # convergence flag, unlike the legacy path"). Mirror the legacy
        # transport self-heating loop exactly: up to 12 fixed-point
        # iterations, 1e-10 K tolerance, explicit convergence flag -- so the
        # power_W reported below (from the SAME diode.junction_power call,
        # default eta_total) is the one actually used to solve T_j, not a
        # stale heat-sink-T estimate.
        if bias_mode == "junction_voltage":
            Tj = _bias_number("T_j_K")
            converged = True
        else:
            Tj = t_junction(duty*diode.junction_power(d.drive.I_uA, ths), a, st, ths)
            converged = False
            for _ in range(12):
                if not np.isfinite(Tj):
                    break
                next_Tj = t_junction(duty*diode.junction_power(d.drive.I_uA, Tj), a, st, ths)
                if abs(next_Tj - Tj) < 1e-10:
                    Tj = next_Tj; converged = True; break
                Tj = next_Tj
        reasons=[]
        va = vj = gam = np.nan
        rr = dict(gamma_X0_ns=np.nan, gamma_XX0_ns=np.nan, k_X_ns=np.nan, k_XX_ns=np.nan,
                 E_a_meV=np.nan, S0=np.nan, tau_cap_ps_used=np.nan, valid=False)
        cr = dict(gamma_X_ns=np.nan, gamma_XX_ns=np.nan, t_X=np.nan, t_XX=np.nan,
                 eta_out=cav.eta_out, kappa_meV=np.nan, detuning_meV=np.nan,
                 Fp_add=np.nan, F_eff_X=np.nan, F_eff_XX=np.nan)
        cnt = dict(mean_counts=np.nan, g2=np.nan, gate_ns_used=requested_gate, converged=False)
        r_dot_val = mu_val = power_val = np.nan
        signal = sx = sxx = bg_counts = np.nan
        rho = g2 = np.nan
        one_pair = False; blocked = np.nan
        feas = {}; priced_fp = np.nan
        supplied = False
        bias = None
        # Invalid-row placeholder: an ALL-NaN NitrideLevels, never a real
        # evaluation of some other (e.g. unbiased, heat-sink-T) state (Opus
        # fix-round finding 3: a thermal-runaway row was reporting finite,
        # "plausible" E_X_eV/field_kVcm/overlap_sq/electron_bound=hole_bound
        # =True for a state that was never the operating point). Overwritten
        # below with the REAL bias-resolved levels only once the operating
        # point is actually reached.
        lv = nitride_levels.NitrideLevels(
            E_X_eV=np.nan, lambda_nm=np.nan, electron_bound=False, hole_bound=False,
            overlap_sq=np.nan, field_kVcm=np.nan, E_e_meV=np.nan, E_h_meV=np.nan,
            dE_e_meV=np.nan, dE_h_meV=np.nan, dE_pair_meV=None,
            sp_split_e_meV=np.nan, sp_split_h_meV=np.nan,
            m_e_matrix_xy=np.nan, m_h_matrix_xy=np.nan, valid=False,
            invalid_reasons=("not evaluated: operating point not reached",),
            provenance="[A] placeholder, row not evaluated at operating point")
        try:
            if not converged:
                raise ValueError("transport self-heating did not converge")
            bias = resolve_bias(diode, T_j_K=Tj,
                                current_uA=d.drive.I_uA if bias_mode == "current" else None,
                                junction_voltage_V=None if bias_mode == "current" else _bias_number("V_j_V"),
                                field_polarity=polarity,
                                external_field_kVcm=system0.external_field_kVcm)
            if not bias["bias_valid"]:
                raise ValueError("invalid bias: " + "; ".join(bias["invalid_reasons"]))
            resolved_current_uA = bias["current_uA"]
            va, vj = bias["V_terminal"], bias["V_j"]
            supplied = resolved_current_uA*1e-6*tau_on*1e-9/E_SI >= 1
            # Bias field: the diode's OWN depletion(vj, Tj).F_kVcm (physical
            # p-i-n depletion field, shrinks toward 0 as V_j -> V_bi), added
            # directly to nitride_levels' intrinsic polarization field --
            # see this function's docstring for the sign convention. Never
            # the lumped -vj/d_i estimate (ignores V_bi, grows the wrong way
            # with forward bias).
            dsys = nitride_levels.NitrideDotSystem(**{**dot_kw,
                "external_field_kVcm": bias["applied_field_kVcm"]})
            lv = nitride_levels.levels(dsys, Tj)
            if not lv.valid:
                raise ValueError("unbound dot: " + "; ".join(lv.invalid_reasons))
            # Fix round 2 (2026-09-09, Opus re-review): ret.tau_cap_scales_with_density
            # was accepted by the card schema but never forwarded here, so the
            # axis was inert on the nitride branch (nitride_levels.rates
            # defaults to False regardless of the card). Forward it explicitly.
            rr = nitride_levels.rates(lv, Tj, tau_rad0_ns=tau_rad0_ns, n_dot_cm2=density, tau_cap_ps=tau_cap, tau_cap_scales_with_density=d.ret.tau_cap_scales_with_density, channel=d.ret.channel, k_nr_ns=float(n.get("k_nr_ns",0.0)))
            if not rr["valid"]:
                raise ValueError("invalid escape rates: " + "; ".join(rr["invalid_reasons"]))
            gam = float(gamma_anchor(Tj, LinewidthParams(d.dot.gamma0,d.dot.a_ac,d.dot.E_LO,d.dot.gamma300)))
            w_val = gam if d.filter.auto_w else d.filter.w  # SAME collection window for cavity accept. and transport bg.
            cr = nitride_cavity.response(cav,T_K=Tj,E_X_eV=lv.E_X_eV,E_X_track_eV=track0.E_X_eV,gamma_X_meV=gam,gamma_XX_meV=gam,delta_xx_meV=d.dot.delta_xx,gamma_X0_ns=rr["gamma_X0_ns"],gamma_XX0_ns=rr["gamma_XX0_ns"],w_meV=w_val,dx_w_meV=d.filter.dx)
            # The reservoir is a material continuum, not E_X plus a fixed
            # offset: it must not follow dot-height QCSE.  See
            # _nitride_reservoir_energy_eV for the [V]/[A] convention.
            reservoir_energy_eV = _nitride_reservoir_energy_eV(dot_kw, Tj, bg_kw)
            reservoir_offset_meV = (reservoir_energy_eV-lv.E_X_eV)*1e3
            # tau_rad_ns is the overlap-SCALED radiative lifetime of this
            # actual (biased) dot state -- tau_rad0_ns/overlap_sq, i.e.
            # 1/gamma_X0_ns -- not the bare reference tau_rad0_ns (which
            # implicitly assumes overlap_sq=1); S_dot=rr["S0"] already uses
            # the same gamma_X0_ns internally, so this keeps the ONE
            # radiative lifetime consistent between the two arguments (Opus
            # fix-round finding: "mixing two radiative lifetimes").
            tau_rad_resolved_ns = tau_rad0_ns / lv.overlap_sq
            inj = nitride_transport.evaluate_injection(diode,resolved_current_uA,Tj,density,aperture,tau_on,w_val, max(reservoir_offset_meV,1.),S_dot=rr["S0"],tau_rad_ns=tau_rad_resolved_ns,E_X_eV=lv.E_X_eV)
            r_dot_val, mu_val, power_val = inj.loading.r_dot, inj.mu, duty*inj.P_junction_W
            if d.drive.cycle_loading == "rectangular":
                cnt = pulse_counting.pulse_g2(inj.loading.r_dot/1e9,cr["gamma_X_ns"],cr["gamma_XX_ns"],rr["k_X_ns"],rr["k_XX_ns"],cr["t_X"],cr["t_XX"],tau_on,period-tau_on,pump_ratio=d.drive.cw_pump_ratio,gate_ns=d.drive.gate_ns,split=True)
                cnt["gate_ns_used"] = requested_gate
            else:
                cnt = pulse_counting.deterministic_cycle_g2(cr["gamma_X_ns"],cr["gamma_XX_ns"],rr["k_X_ns"],rr["k_XX_ns"],cr["t_X"],cr["t_XX"],period,eta_load=d.drive.eta_load,gate_ns=d.drive.gate_ns,split=True)
                one_pair=cnt["one_pair_valid"]; blocked=cnt["blocked_load_probability"]
            if not cnt.get("converged", True):
                raise ValueError("photon counting did not converge")
            gate=cnt.get("gate_ns_used",period); bg_window=min(gate,tau_on)
            # Reservoir RATE (photons/s, transport.Background.rate_bg_window --
            # NOT rate_x, which is the dot's own X-channel normalization) is
            # integrated over the injection window plus optional exponential
            # afterglow, then filtered by the cavity/slit's flat-spectrum
            # acceptance (not t_X, spec) and collected once via eta_background [A].
            bg_accept = _nitride_flat_background_acceptance(cr["kappa_meV"], w_val, d.filter.dx, cr["detuning_meV"], reservoir_offset_meV)
            bg_counts=inj.background.rate_bg_window*1e-9*bg_window
            if bg_tau>0 and gate>tau_on: bg_counts += inj.background.rate_bg_window*1e-9*bg_tau*(1-np.exp(-(gate-tau_on)/bg_tau))
            bg_counts *= eta_bg*bg_accept
            signal=cr["eta_out"]*cnt["mean_counts"]; sx=cr["eta_out"]*cnt.get("mean_counts_x",np.nan); sxx=cr["eta_out"]*cnt.get("mean_counts_xx",np.nan)
            # Fix round 2 (2026-09-09, Opus re-review): drive.b_res was read
            # only on the legacy transport path (this module's "item 8"
            # comments below evaluate()), leaving it inert here. Same "per
            # collected X photon" residual-background convention as the
            # legacy branch: b_res*n_X in absolute count units, added
            # directly to bg_counts -- sx is already the eta_out-collected X
            # count, so no further G/S conversion factor applies.
            if np.isfinite(sx): bg_counts += d.drive.b_res*sx
            rho=signal/(signal+bg_counts) if signal+bg_counts>0 else np.nan
            g2=1-rho*rho*(1-cnt["g2"]) if np.isfinite(rho) and np.isfinite(cnt["g2"]) else np.nan
            if not np.isfinite(g2):
                raise ValueError("non-finite g2 at operating point (zero collected signal+background or non-finite counting result)")
            if d.drive.cycle_loading == "deterministic_pair":
                # sp0 (validated as complete above, outside the loop) plus
                # the resolved cycle rate for this row.
                sp=dict(sp0); sp["f_cycle_Hz"]=rep
                feas=set_feasibility(Tj,**sp)
                # ALWAYS priced against the SAME turnstile mechanism (spec); read
                # the delivered F_p back from the mechanism interface rather than
                # re-deriving the feasible/infeasible branch a second time here.
                mech_iface = mech_set("turnstile",Tj,eps_cycle=0.0,**sp)
                priced_fp = float(mech_iface.F_p)
                # Max island radius still satisfying the charging-energy screen
                # (E_C >= ec_margin*kT), inverted through the SAME e^2/C and
                # isolated-sphere convention (loading.island_radius_nm) that
                # set_feasibility itself uses -- not e^2/(2C).
                C_sigma_max = E_SI**2 / (sp["ec_margin"] * KB_SI * Tj)
                feas["radius_max_nm"] = float(island_radius_nm(C_sigma_max, sp["eps_r"]))
        except ValueError as exc:
            if str(exc) not in reasons: reasons.append(str(exc))
        # Unbound dot, non-converged thermal/counting solution or impossible
        # field is an explicit invalid row (never a silent InP-class
        # fallback): gate on convergence explicitly, not only on whatever
        # happened to still be finite downstream. Every branch above that can
        # make `valid` False also raises with its own specific reason, so
        # this generic fallback (deliberately NOT the literal ['invalid
        # row'] the Opus fix-round flagged) should be unreachable; it stays
        # as a defensive net that names itself as such rather than
        # pretending to be a diagnosed physical cause.
        valid=bool(converged and lv.valid and rr["valid"] and cnt.get("converged",True) and np.isfinite(g2))
        if not valid and not reasons:
            reasons.append("invalid row: cause not captured by an explicit check (internal inconsistency)")
        rows.append(dict(T_hs=ths,T_j=Tj,g2=g2,rho=rho,lv=lv,rr=rr,cr=cr,
                         r_dot=r_dot_val,mu=mu_val,power=power_val,
                         cnt=cnt,signal=signal,sx=sx,sxx=sxx,bg=bg_counts,valid=valid,reasons=reasons,
                         feas=feas,priced_fp=priced_fp,supplied=supplied,one_pair=one_pair,blocked=blocked,vj=vj,gam=gam,
                         bias=bias, reservoir_energy_eV=locals().get("reservoir_energy_eV", np.nan),
                         reservoir_offset_meV=locals().get("reservoir_offset_meV", np.nan)))
    keys={"g2":lambda r:r["g2"],"rho2":lambda r:r["rho"]**2,"Tj":lambda r:r["T_j"],"gamma":lambda r:r["gam"],"eps":lambda r:r["cr"]["t_XX"]/r["cr"]["t_X"],"signal_flux":lambda r:r["signal"]*rep}
    curves={k:np.array([f(r) for r in rows]) for k,f in keys.items()}; op=rows[int(np.argmin(abs(ts-d.thermal.T_hs)))]
    x=op; c=x["cr"]; rr=x["rr"]; lv=x["lv"]; feas=x["feas"]
    scalars={"platform":"ingan_gan_planar","cycle_loading":d.drive.cycle_loading,"T_hs":x["T_hs"],"T_j":x["T_j"],"g2_op":x["g2"],"rho_pulsed":x["rho"],"collected_flux_pulsed_s":x["signal"]*rep,"collected_flux_x_s":x["sx"]*rep,"collected_flux_xx_s":x["sxx"]*rep,"background_flux_s":x["bg"]*rep,"total_detected_flux_s":(x["signal"]+x["bg"])*rep,"mean_counts":x["cnt"]["mean_counts"],"mean_counts_x":x["cnt"].get("mean_counts_x",np.nan),"mean_counts_xx":x["cnt"].get("mean_counts_xx",np.nan),"E_X_eV":lv.E_X_eV,"lambda_nm":lv.lambda_nm,"field_kVcm":lv.field_kVcm,"overlap_sq":lv.overlap_sq,"electron_bound":lv.electron_bound,"hole_bound":lv.hole_bound,"E_a_meV":rr["E_a_meV"],"k_X_ns":rr["k_X_ns"],"k_XX_ns":rr["k_XX_ns"],"gamma_X0_ns":rr["gamma_X0_ns"],"gamma_XX0_ns":rr["gamma_XX0_ns"],"gamma_X_ns":c["gamma_X_ns"],"gamma_XX_ns":c["gamma_XX_ns"],"S_X":c["gamma_X_ns"]/(c["gamma_X_ns"]+rr["k_X_ns"]),"S_XX":c["gamma_XX_ns"]/(c["gamma_XX_ns"]+rr["k_XX_ns"]),"Q":cav.Q,"kappa_meV":c["kappa_meV"],"detuning_meV":c["detuning_meV"],"Fp_add":c["Fp_add"],"F_eff_X":c["F_eff_X"],"F_eff_XX":c["F_eff_XX"],"eta_out":c["eta_out"],"gate_ns_used":x["cnt"]["gate_ns_used"],"rep_rate_hz":rep,"r_dot_s":x["r_dot"],"mu_resolved":x["mu"],"n_dot_cm2_used":density,"tau_cap_ps_used":rr["tau_cap_ps_used"],"I_pair_pA":E_SI*rep*1e12,"transport_current_uA":d.drive.I_uA,"resolved_current_uA":x["bias"]["current_uA"] if x["bias"] else np.nan,"V_j":x["vj"],"V_terminal":x["bias"]["V_terminal"] if x["bias"] else np.nan,"power_W":x["power"],"counting_converged":x["cnt"].get("converged",True),"blocked_load_probability":x["blocked"],"one_pair_valid":x["one_pair"],"set_feasible":feas.get("feasible",False),"set_priced_F_p":x["priced_fp"],"set_E_C_meV":feas.get("E_C_meV",np.nan),"set_EC_over_kT":feas.get("EC_over_kT",np.nan),"set_radius_nm":feas.get("radius_nm",np.nan),"set_radius_max_nm":feas.get("radius_max_nm",np.nan),"set_C_sigma_F":feas.get("C_sigma_F",np.nan),"set_R_T_over_RQ":feas.get("R_T_over_RQ",np.nan),"set_f_max_Hz":feas.get("f_max_Hz",np.nan),"pair_supply_possible":x["supplied"],"ideal_load_F_p":0.0 if d.drive.cycle_loading=="deterministic_pair" else np.nan,"valid":x["valid"],"invalid_reasons":x["reasons"],"provenance":{"nitride":"[A/E/DR] planar integration; cavity/transport inputs retain module provenance"}}
    spectroscopy_valid = bool(lv.valid and rr.get("valid", False) and np.isfinite(c["gamma_X_ns"]) and c["gamma_X_ns"] > 0)
    scalars.update({
        "tau_rad_bare_ns": (1.0 / rr["gamma_X0_ns"] if rr["gamma_X0_ns"] > 0 else np.inf),
        "tau_rad_cavity_ns": (1.0 / c["gamma_X_ns"] if c["gamma_X_ns"] > 0 else np.inf),
        "spectroscopy_valid": spectroscopy_valid,
        "spectroscopy_invalid_reasons": [] if spectroscopy_valid else list(x["reasons"]),
        "temperature_mode": "fixed_junction" if bias_mode == "junction_voltage" else "self_consistent",
        "evaluation_kind": "stark_diagnostic" if bias_mode == "junction_voltage" else "source",
        "field_polarity": polarity,
        "diode_field_kVcm": x["bias"]["diode_field_kVcm"] if x["bias"] else np.nan,
        "applied_field_kVcm": x["bias"]["applied_field_kVcm"] if x["bias"] else np.nan,
        "reservoir_energy_eV": x["reservoir_energy_eV"], "reservoir_offset_meV": x["reservoir_offset_meV"],
        "reservoir_kind": lv.reservoir_kind, "geometry_type": system0.geometry_type,
        "effective_height_nm": lv.effective_height_nm, "effective_radius_nm": lv.effective_radius_nm,
        "cavity_reference_V_j_V": _track_bias["V_j"], "cavity_reference_transition_eV": track0.E_X_eV,
        "cavity_reference_convention": "fixed_junction_voltage" if cavity_reference_V is not None or bias_mode == "junction_voltage" else "current_controlled",
    })
    scalars["device_pass"]=bool(bias_mode == "current" and scalars["valid"] and scalars["g2_op"]<.5 and scalars["collected_flux_pulsed_s"]>=1000 and (d.drive.cycle_loading!="deterministic_pair" or (scalars["one_pair_valid"] and scalars["set_feasible"] and scalars["pair_supply_possible"])))
    return {"curves":curves,"scalars":scalars}


def evaluate(design: DeviceDesign, T_grid=None) -> dict:
    """Run the full chain. Returns {'curves': {...}, 'scalars': {...}}."""
    d = design
    if d.platform == "ingan_gan_planar":
        return _evaluate_nitride(d, T_grid)
    if d.platform != "legacy":
        raise ValueError(f"unknown platform {d.platform!r}")
    # cycle_loading's "deterministic_pair" opt-in is nitride-only (spec);
    # a legacy design never reads drive.cycle_loading below, so silently
    # accepting a non-default value here would be misleading, not merely
    # inert -- every legacy card leaves this at its "rectangular" default.
    if d.drive.cycle_loading != "rectangular":
        raise ValueError("drive.cycle_loading is nitride-only (platform='ingan_gan_planar')")
    if d.drive.mode == "PL" and d.drive.diode:
        raise ValueError("drive.mode='PL' cannot be used with a non-empty diode")
    if d.drive.mode == "EL-transport" and d.drive.dg_inj:
        raise ValueError("drive.mode='EL-transport' cannot be used with dg_inj")
    if d.drive.mode == "EL-transport" and d.drive.mechanism:
        raise ValueError("drive.mode='EL-transport' cannot be combined with "
                         "drive.mechanism (conflicting loading resolutions)")
    # pr-pkg4-fix item 6: drive.mechanism resolves its own (mu, F_p, eta),
    # overriding transport's own -- finite_pulse's rates below (r_ns_fp from
    # injection.loading.r_dot, k_X/k_XX from the confinement retention
    # params) are wired straight to the injection/confinement chain and
    # never consult a mechanism's DriveInterface at all, so combining the
    # two would silently ignore the mechanism rather than resolve it into
    # the finite-pulse pump rate.
    if d.drive.finite_pulse and d.drive.mechanism:
        raise ValueError("finite_pulse does not support drive.mechanism overrides")
    if d.dot.linewidth not in ("class", "anchored"):
        raise ValueError(f"unknown dot.linewidth {d.dot.linewidth!r}")
    if d.ret.mode not in ("proxy", "confinement"):
        raise ValueError(f"unknown ret.mode {d.ret.mode!r}")
    if d.drive.mode not in ("EL", "PL", "EL-transport"):
        raise ValueError(f"unknown drive.mode {d.drive.mode!r}")
    if d.drive.loading_model not in ("auto", "capped_poisson", "moment_matched"):
        raise ValueError(f"unknown drive.loading_model {d.drive.loading_model!r}")
    if d.drive.gate_ns is not None and d.drive.gate_ns <= 0:
        raise ValueError(f"drive.gate_ns must be > 0 if set (got {d.drive.gate_ns!r})")
    if d.emission.type not in ("none", "edge"):
        raise ValueError(f"unknown emission.type {d.emission.type!r}")
    if d.filter.track_material not in ("", "dot", "matrix"):
        raise ValueError(f"unknown filter.track_material {d.filter.track_material!r}")
    # pr-pkg6-fix item 1 (pr-pkg6-stale-text item C2): both fields are
    # validated > 0 if set at DeviceDesign.load() (see load()'s own checks),
    # but load() is not the only way to build a DeviceDesign -- a
    # programmatic caller that constructs one directly and sets either field
    # to exactly 0.0 (or negative) skipped that check entirely. The
    # ZeroDivisionError this used to reach was NOT confinement-only: it also
    # fires ~170 lines before the (formerly confinement-only) re-validation
    # even ran, inside transport.evaluate_injection's dot_loading (f_qd:
    # 1e10 / n_dot_cm2) for EL-transport under ANY ret.mode, not just
    # "confinement" (dot_levels.retention_params' states_per_dot division is
    # the confinement-only site). Validated here, before either consumer
    # runs, for every mode that actually reaches one of those two divisions
    # -- ret.mode == "confinement" (retention_params, both fields) or
    # drive.mode == "EL-transport" (dot_loading's f_qd, aperture.density_cm2
    # only). aperture.compose's OWN use of aperture.density_cm2
    # (n_window_competitors: a pure product, n_qd_cm2 * aperture_um2 *
    # w/sigma, never a divisor) is a THIRD, unrelated consumer that is
    # genuinely fine at 0.0 -- N=0 competitors is exactly how "continuous
    # aperture N=0 reproduces the unmixed g2 exactly" (verify_device_rt.py)
    # exercises it under the untouched EL/proxy defaults, so this guard must
    # NOT fire for that design (neither ret.mode == "confinement" nor
    # drive.mode == "EL-transport" there).
    #
    # pr-pkg4-fix4 item 2 (Opus review): a NEGATIVE aperture.density_cm2 is
    # unlike the `== 0.0` case above -- it is not "genuinely fine" for
    # aperture.compose's own product-only use either (n_window_competitors
    # would silently compute a negative Nw, not raise), so it is checked
    # UNCONDITIONALLY, in every mode, ahead of (and independent from) the
    # mode-gated `== 0.0` guard below.
    if d.aperture.density_cm2 is not None and d.aperture.density_cm2 < 0:
        raise ValueError("aperture.density_cm2 must be non-negative if set, "
                         "in every mode (a negative competitor density has "
                         f"no physical meaning) (got {d.aperture.density_cm2!r})")
    if d.ret.mode == "confinement" or d.drive.mode == "EL-transport":
        if d.aperture.density_cm2 is not None and d.aperture.density_cm2 == 0.0:
            raise ValueError("aperture.density_cm2 must be > 0 if set in "
                             "ret.mode='confinement' or drive.mode='EL-transport' "
                             f"(got {d.aperture.density_cm2!r}; 0.0 is only "
                             "valid outside those modes, where aperture.compose's "
                             "own product-only use of it is genuinely fine at 0)")
        if d.ret.n_dot_cm2 is not None and d.ret.n_dot_cm2 <= 0:
            raise ValueError("ret.n_dot_cm2 must be > 0 if set "
                             f"(got {d.ret.n_dot_cm2!r})")
    if d.emission.type == "edge" and d.cavity.enabled:
        # Reject conflicting SiN/edge and resonant-cavity collection (docs/
        # rt_edge_contract.md): coexistence needs an explicit supported
        # leakage accounting this tier does not implement.
        raise ValueError("emission.type='edge' cannot be combined with "
                         "cavity.enabled (no supported cavity-leakage accounting)")
    if d.drive.cw and d.drive.mode != "EL-transport":
        # The CW pump rate is transport.loading.r_dot; never inferred from
        # the dimensionless pulsed mu (docs/rt_edge_contract.md).
        raise ValueError("drive.cw=True requires drive.mode='EL-transport'")
    # council review 2026-09-05 item 1: DriveBlock's own defaults (duty=1.0,
    # drive.diode's tau_pulse_ns default 1.0) silently read as a 1 GHz, 100%
    # duty DC drive under EL-transport pulsed operation (device.py:944-949,
    # pre-fix) -- a card that never touched either field got a fabricated
    # repetition rate with no warning. Pulsed (drive.cw=False) EL-transport
    # now REQUIRES the card to set drive.diode['tau_pulse_ns'] and either
    # drive.duty or drive.rep_rate_hz explicitly (away from their untouched
    # defaults); CW is exempt (no pulsing, no repetition rate to fabricate).
    tau_pulse_ns_val = float((d.drive.diode or {}).get("tau_pulse_ns", 1.0))
    tau_pulse_explicit = "tau_pulse_ns" in (d.drive.diode or {})
    duty_explicit = d.drive.duty != DriveBlock.duty
    rep_explicit = d.drive.rep_rate_hz > 0.0
    # pr-pkg4-fix item 5: drive.finite_pulse needs this SAME explicit pulse
    # period even when drive.cw=True -- the finite-pulse block below reads
    # tau_pulse_ns_val/duty_eff/rep_rate_hz exactly like the legacy pulsed
    # path, so a CW card that opts into finite_pulse without ALSO stating an
    # explicit tau_pulse_ns and duty/rep_rate_hz would otherwise fabricate
    # the same silent 1 GHz/100%-duty DC drive this guard exists to prevent
    # (previously only the drive.cw=False branch below was covered).
    if d.drive.mode == "EL-transport" and (not d.drive.cw or d.drive.finite_pulse):
        if not (tau_pulse_explicit and (duty_explicit or rep_explicit)):
            missing = []
            if not tau_pulse_explicit:
                missing.append("drive.diode['tau_pulse_ns']")
            if not (duty_explicit or rep_explicit):
                missing.append("drive.duty or drive.rep_rate_hz")
            reason = ("pulsed operation (drive.cw=False)" if not d.drive.cw
                     else "drive.finite_pulse=True (needs a real pulse period even "
                          "under drive.cw=True)")
            raise ValueError(
                f"drive.mode='EL-transport' {reason} requires "
                "explicit " + " and ".join(missing) + " -- unset defaults (duty=1.0, "
                "tau_pulse_ns=1.0) silently report a 1 GHz/100%-duty DC drive "
                "(council review 2026-09-05 item 1)")
    # duty_eff feeds the thermal-power calculation below (P, and the
    # transport self-heating loop's own duty*P_junction average) in place of
    # the raw d.drive.duty field: when the card gave rep_rate_hz explicitly
    # instead of duty, the physically consistent duty is tau_pulse*rep_rate,
    # never the untouched duty=1.0 default. Whenever rep_rate_hz is NOT the
    # (sole) explicit choice -- every legacy caller, since rep_rate_hz
    # defaults to 0 -- duty_eff is exactly d.drive.duty (bit-identical).
    duty_eff = (d.drive.rep_rate_hz * tau_pulse_ns_val * 1e-9
               if (rep_explicit and not duty_explicit) else d.drive.duty)
    proxy = class_proxy_params() if (d.dot.linewidth == "class" or d.ret.mode == "proxy") else {}
    gp = ({k: proxy[k] for k in ("gamma0", "a_ac", "b_lo", "E_lo")}
          if d.dot.linewidth == "class" else None)
    rp = ({k: (getattr(d.ret, k) or proxy[k]) for k in ("a_esc", "E_a", "b_p", "b0", "beta")}
          if d.ret.mode == "proxy" else {})
    if d.ret.mode == "proxy":
        rp["E_b"] = d.ret.E_b or proxy["E_b"]
    retention_note = "class-proxy Arrhenius fit [A]"
    diode = (_diode_from_drive(d.drive, d.aperture.diameter_um)
             if d.drive.mode == "EL-transport" else None)

    a = 0.5 * d.thermal.mesa_diameter_um * 1e-6
    st = _stack(d.thermal)
    P = duty_eff * d.drive.I_uA * 1e-6 * d.drive.V
    chan = [BackgroundChannel("injection", A=d.drive.b_e, m=d.drive.b_e_m,
                              E_act=d.drive.b_e_Eact, I_ref=d.drive.I_ref_uA)]
    # SiN evanescent coupling: no resonant line (no kappa acceptance, no G boost);
    # beta_sin is applied to brightness only, below -- never to eps/rho/g2/T_c (Lemma 1).
    sin_mode = d.cavity.enabled and d.cavity.type == "sin_waveguide"
    # council review 2026-09-05 item 7: filter.track_material only actually
    # moves the filter window (a) whenever the cavity tracks it (cavity
    # enabled, non-SiN -- dx is always recomputed from the stack there), or
    # (b) cavity-less/SiN with the explicit hold_window opt-in. Cavity-less
    # with hold_window=False (the default) resolves _tracked_material only to
    # validate the name -- dx stays filter.dx, exactly as if track_material
    # were never set -- so the provenance note below must say so, not claim
    # "materials.bandgap on stack layer ..." as if it were effective.
    cavity_tracks = d.cavity.enabled and not sin_mode
    track_material_effective = bool(d.filter.track_material) and (cavity_tracks or d.filter.hold_window)

    def _nan_fp_rates():
        # finding 1 -- rates pulse_counting.pulse_g2 was (or would have been)
        # called with; all-nan whenever drive.finite_pulse is off or the
        # operating point is invalid before reaching that call.
        return dict(r_ns=np.nan, gamma_X_ns=np.nan, gamma_XX_ns=np.nan,
                    k_X=np.nan, k_XX=np.nan, tau_on_ns=np.nan, tau_dark_ns=np.nan)

    def one(T_hs):
        # The legacy branch deliberately retains its one-shot P=IV calculation.
        # Transport uses the diode's junction power in a short fixed-point loop.
        Tj = t_junction(P, a, st, T_hs)
        if diode is not None:
            if d.drive.I_uA <= 0:
                return dict(Tj=np.nan, gam=np.nan, eps=np.nan, rho=np.nan,
                            g2=np.nan, t_x=np.nan, runaway=False, mu=np.nan,
                            eta_capture=np.nan, b_e=np.nan, injection=None, S=np.nan,
                            g2_cw0=np.nan, g2_cw0_raw=np.nan, cw_r_ns=np.nan,
                            cw_gamma_X_ns=np.nan, cw_rho=np.nan, w=np.nan,
                            n_dot_cm2_used=np.nan,
                            finite_pulse_g2_dot=np.nan, finite_pulse_mean_counts=np.nan,
                            finite_pulse_mean_counts_x=np.nan, finite_pulse_rates=_nan_fp_rates(),
                            finite_pulse_gate_ns_used=np.nan, finite_pulse_converged=False,
                            rho_pulsed=np.nan, retention_params_used=None,
                            invalid_reason="EL-transport requires positive current")
            converged = False
            for _ in range(12):
                trial = transport.evaluate_injection(
                    diode=diode, I_uA=d.drive.I_uA, T=Tj,
                    n_dot_cm2=d.drive.n_dot_cm2 or _legacy_density_cm2(d.aperture.density_cm2),
                    aperture_um2=np.pi * (d.aperture.diameter_um / 2) ** 2,
                    **_transport_options(d.drive, 1.0))
                next_Tj = t_junction(duty_eff * trial.P_junction_W, a, st, T_hs)
                if abs(next_Tj - Tj) < 1e-10:
                    converged = True
                    break
                Tj = next_Tj
            if not converged:
                return dict(Tj=np.nan, gam=np.nan, eps=np.nan, rho=np.nan,
                            g2=np.nan, t_x=np.nan, runaway=False, mu=np.nan,
                            eta_capture=np.nan, b_e=np.nan, injection=None, S=np.nan,
                            g2_cw0=np.nan, g2_cw0_raw=np.nan, cw_r_ns=np.nan,
                            cw_gamma_X_ns=np.nan, cw_rho=np.nan, w=np.nan,
                            n_dot_cm2_used=np.nan,
                            finite_pulse_g2_dot=np.nan, finite_pulse_mean_counts=np.nan,
                            finite_pulse_mean_counts_x=np.nan, finite_pulse_rates=_nan_fp_rates(),
                            finite_pulse_gate_ns_used=np.nan, finite_pulse_converged=False,
                            rho_pulsed=np.nan, retention_params_used=None,
                            invalid_reason="transport self-heating did not converge")
        if not np.isfinite(Tj):
            return dict(Tj=np.inf, gam=np.nan, eps=np.nan, rho=np.nan,
                        g2=np.nan, t_x=np.nan, runaway=True, mu=np.nan,
                        eta_capture=np.nan, b_e=np.nan, injection=None, S=np.nan,
                        g2_cw0=np.nan, g2_cw0_raw=np.nan, cw_r_ns=np.nan,
                        cw_gamma_X_ns=np.nan, cw_rho=np.nan, w=np.nan,
                        n_dot_cm2_used=np.nan,
                        finite_pulse_g2_dot=np.nan, finite_pulse_mean_counts=np.nan,
                        finite_pulse_mean_counts_x=np.nan, finite_pulse_rates=_nan_fp_rates(),
                        finite_pulse_gate_ns_used=np.nan, finite_pulse_converged=False,
                        rho_pulsed=np.nan, retention_params_used=None,
                        invalid_reason="thermal runaway")
        gam_base = (float(gamma_anchor(Tj, LinewidthParams(d.dot.gamma0, d.dot.a_ac,
                    d.dot.E_LO, d.dot.gamma300))) if d.dot.linewidth == "anchored"
                    else float(gamma_of_T(Tj, **gp)))
        gam = d.dot.gamma_scale * gam_base
        # F5' injection broadening (v1.1): EL-mode-only, applied BEFORE the
        # auto-w filter width (below) and BEFORE epsilon.
        if d.drive.mode != "PL" and d.drive.dg_inj:
            gam = gamma_eff(gam, d.drive.dg_inj, d.drive.I_uA / d.drive.I_ref_uA, d.drive.p_inj)
        kappa = d.cavity.kappa if (d.cavity.enabled and not sin_mode) else None
        dx = d.filter.dx
        if d.cavity.enabled and not sin_mode:  # sin waveguide has no mode to track
            if d.filter.track_material:
                # Stack-based tracking [DR]: the E_X0/T=0 anchor cancels
                # between the two bandgap differences below, so only the
                # named layer's OWN Varshni shape enters -- never GaAs unless
                # the stack's own dot/matrix material happens to be GaAs.
                mat = _tracked_material(d.ret, d.filter.track_material)
                dx = (1e3 * (materials.bandgap(mat, Tj) - materials.bandgap(mat, d.cavity.T_track))
                      - d.cavity.dEdT_cav * (Tj - d.cavity.T_track))
            else:
                dx = 1e3 * float(tracking_detuning(Tj, d.cavity.E_X0, d.cavity.T_track,
                                                   "GaAs", d.cavity.dEdT_cav * 1e-3))
        elif d.filter.track_material:
            # item 4, round 2 (council review 2026-09-05): round 1's fix here
            # made track_material do SOMETHING with cavity.enabled=False (it
            # was previously gated behind cavity.enabled and did nothing at
            # all), but the something was wrong: it held the filter window
            # FIXED at its cavity.T_track position and let the dot's own
            # Varshni-walked line drift out of it (t_X 0.5 -> 0.00106 at the
            # worst sweep corner -- the regression's largest single term). A
            # filter that TRACKS the dot is, by definition, centred ON the
            # dot at every T: dx stays d.filter.dx (0 by default; any
            # explicit manual offset still applies), exactly like the
            # track_material == "" path -- there is no cavity mode here for a
            # tracking filter to net against, so nothing computed from the
            # material's Varshni shape belongs in dx by default.
            # filter.hold_window=True is the explicit, differently-named
            # opt-in for the round-1 behaviour (a genuinely fixed window,
            # e.g. a monochromator slit calibrated at T_track and never
            # retuned): _tracked_material is still resolved unconditionally
            # so an unknown track_material name is rejected either way.
            mat = _tracked_material(d.ret, d.filter.track_material)
            if d.filter.hold_window:
                dx = 1e3 * (materials.bandgap(mat, Tj) - materials.bandgap(mat, d.cavity.T_track))
            else:
                dx = d.filter.dx
        w = None
        if d.filter.enabled:
            # item 9 (council review 2026-09-05): auto_w_scale multiplies the
            # auto_w (full-linewidth) convention window; 1.0 (default) is
            # exactly legacy.
            w = gam * d.filter.auto_w_scale if d.filter.auto_w else d.filter.w
        # Slit-held tracking (T-1): slit centered on X (dx_w = user offset
        # only), cavity at the physical mode walk (dx stays the tracking
        # detuning). Under "mode" (legacy) both centers share dx exactly.
        held = (d.filter.track == "hold" and d.filter.enabled
                and d.cavity.enabled and not sin_mode)
        dx_w = d.filter.dx if held else dx
        # R2 opt-in: overlap-penalized Purcell at this operating point,
        # scaled by the cavity Lorentzian at the X detuning (a mode far off
        # the line enhances nothing -- the factor that kills the 77 K +
        # 120 K-tracked wiring honestly).
        wire = (d.cavity.purcell_wire and d.cavity.enabled and not sin_mode)
        if wire:
            det = (kappa / 2.0) ** 2 / (dx**2 + (kappa / 2.0) ** 2)
            F_eff_op = 1.0 + (float(purcell_eff(d.cavity.F_P, kappa, gam))
                              - 1.0) * det
        else:
            F_eff_op = 1.0
        rate_mult = 1.0
        if d.dot.lineshape == "ibm":
            # T2 tier: sideband-aware transmissions. Same acceptance stack and
            # delta conventions as spectral.epsilon; the XX line reuses the X
            # exciton's PhononParams (single-form-factor approximation, [E] --
            # the XX phonon coupling differs at the O(1)-factor level, a
            # refinement that waits for data). ZPL widths carry the standing
            # Gamma(T) model unchanged (see qd_gf.effective_gamma_zpl note).
            pp = PhononParams(**d.dot.phonon)
            if wire:
                # R2: Purcell funnels emission into the ZPL (Z_eff) and
                # multiplies the radiative rate; X line only (the XX line's
                # own Purcell overlap differs -- kept unenhanced, [E]).
                t_x, _z_eff, rate_mult = ibm_purcell_transmission(
                    dx_w, pp, Tj, gam, F_eff_op, w_meV=w, kappa_meV=kappa,
                    delta_c_meV=dx if held else None)
            else:
                t_x = float(ibm_transmission(dx_w, pp, Tj, gam, w_meV=w,
                                             kappa_meV=kappa,
                                             delta_c_meV=dx if held else None))
            t_xx = float(ibm_transmission(dx_w - d.dot.delta_xx, pp, Tj,
                                          d.dot.r_xx * gam, w_meV=w,
                                          kappa_meV=kappa,
                                          delta_c_meV=(dx - d.dot.delta_xx)
                                          if held else None))
            spec = SpectralResult(eps=t_xx / t_x, t_x=t_x, t_xx=t_xx)
        elif held:
            spec = epsilon2(d.dot.delta_xx, gam, d.dot.r_xx * gam, w=w,
                            kappa=kappa, dx_w=dx_w, dx_c=dx)
        else:
            spec = epsilon(d.dot.delta_xx, gam, d.dot.r_xx * gam, w=w, kappa=kappa, dx=dx)
        if d.ret.mode == "confinement":
            # Finding 1b: the SAME EXPLICIT dot density transport.
            # evaluate_injection uses below (device.py's documented aperture/
            # drive convention) is forwarded here -- unlike transport's own
            # call site, an unset (None) aperture.density_cm2 is passed
            # through as-is (not resolved to the 7.0e8 legacy default;
            # item 1), so a legacy card with no explicit density falls to
            # _confinement_params' real 1e10 default instead of silently
            # inheriting transport's unrelated 7.0e8 class value.
            #
            # pr-pkg1-fix3 item 2: drive.n_dot_cm2 is TRANSPORT-only (the EL
            # injection density, RetentionBlock's own docstring and
            # _legacy_density_cm2's callers agree) -- it must never reach
            # confinement. The previous `d.drive.n_dot_cm2 or
            # d.aperture.density_cm2` forwarded it whenever a legacy design
            # set drive.n_dot_cm2 but left aperture.density_cm2 unset,
            # silently changing that design's confinement retention (S
            # 9.637738e-4 -> 1.021141e-4 for the InP/GaAsP0.4/AlGaAs0.4
            # confinement preset at drive.n_dot_cm2=1e9 -- see
            # verify_device_rt.py's dedicated legacy-invariance check).
            # aperture.density_cm2 is forwarded as-is (None passes through
            # unresolved, per the comment above); _confinement_params
            # itself is_not_None-tests it against ret.n_dot_cm2 first, so
            # this value wins over the 1e10 class default whenever
            # ret.n_dot_cm2 itself is unset. pr-pkg1-fix4 item 5: both
            # aperture.density_cm2 and ret.n_dot_cm2 are validated > 0 if set
            # at DeviceDesign.load() (see load()'s own checks) -- but load()
            # is not the only way to build a DeviceDesign, and a 0.0 here
            # feeds dot_levels.retention_params' states_per_dot (N2D/
            # n_dot_cm2) division (this branch only runs for ret.mode ==
            # "confinement", so that condition already holds here).
            # pr-pkg6-fix item 1: that re-validation is no longer done here
            # -- it is hoisted to the top of evaluate() (before ANY use,
            # gated on ret.mode == "confinement" or drive.mode ==
            # "EL-transport", the two consumers that actually divide by
            # these fields) since transport.evaluate_injection's own
            # dot_loading (f_qd, EL-transport, ANY ret.mode) divides by the
            # same aperture.density_cm2 ~170 lines before this branch is
            # even reached, so a confinement-only guard here left that
            # earlier ZeroDivisionError reachable under ret.mode ==
            # "proxy" + drive.mode == "EL-transport" (pr-pkg6-stale-text
            # item C2).
            derived = _confinement_params(d.ret, Tj, n_dot_cm2=d.aperture.density_cm2)
            params = {k: derived[k] for k in ("a_esc", "E_a", "b_p", "E_b")}
            params["b0"], params["beta"] = d.ret.b0, d.ret.beta
            # ret.overrides is an explicit dict (distinct from ret.system, which
            # only ever carries the inline stack): its presence, not its
            # truthiness, decides whether a coefficient is replaced, so an
            # explicitly supplied zero always wins -- unlike legacy's
            # 0="use proxy" shorthand on the RetentionBlock fields themselves.
            for key, value in d.ret.overrides.items():
                if key not in ("a_esc", "E_a", "b_p", "E_b", "b0", "beta"):
                    raise ValueError(f"unsupported ret.overrides key {key!r}")
                params[key] = float(value)
            # item 2: the confinement-derived transition energy feeds
            # transport's sub-turn-on loading suppression below; proxy mode
            # has no confinement level structure, so E_X_eV stays None there
            # (transport.dot_loading's documented legacy-numerics default).
            E_X_eV = derived["E_X_eV"]
            n_dot_cm2_used = derived["n_dot_cm2_used"]
        else:
            params = rp
            E_X_eV = None
            n_dot_cm2_used = float("nan")  # proxy mode: no confinement density used
        # R2: Purcell speeds the radiative channel by rate_mult (photon-
        # weighted; Lorentzian path uses F_eff directly since there is no
        # sideband split), so the escape-to-radiative ratios divide by it.
        # rm == 1.0 whenever `wire` is False (rate_mult/F_eff_op are then
        # left at their 1.0 initializers) -- also reused below by the CW
        # gamma_X_ns "radiative enhancement" (docs/rt_edge_contract.md).
        rm = rate_mult if d.dot.lineshape == "ibm" else F_eff_op
        if wire:
            S = float(retention(Tj, params["a_esc"] / rm, params["E_a"],
                                params["b_p"] / rm, params["E_b"]))
        else:
            S = float(retention(Tj, params["a_esc"], params["E_a"], params["b_p"], params["E_b"]))
        # PL mode: optical excitation, no injection-current background channel.
        injection = None
        if diode is not None:
            w_f = w if w is not None else 1.0
            opts = _transport_options(d.drive, w_f)
            injection = transport.evaluate_injection(
                diode=diode, I_uA=d.drive.I_uA, T=Tj,
                n_dot_cm2=d.drive.n_dot_cm2 or _legacy_density_cm2(d.aperture.density_cm2),
                aperture_um2=np.pi * (d.aperture.diameter_um / 2) ** 2,
                S_dot=S, E_X_eV=E_X_eV, **opts)
            # item 5 (council review 2026-09-05): injection.b_e integrates the
            # matrix/WL continuum over the TRANSPORT window opts['w_meV']
            # (== w_f unless drive.diode.w_meV explicitly overrides it, e.g.
            # to model a collection bandwidth different from the filter's
            # own). A flat continuum's collected rate scales linearly with
            # the actual filter window width -- never by the X line's own
            # transmission (that would be item 5's CW-path bug, below, for a
            # broadband background) -- so rescale by min(1, w_f/w_transport);
            # a no-op (1.0) whenever no such override is given.
            w_transport = opts["w_meV"]
            win_scale = min(1.0, w_f / w_transport) if w_transport > 0 else 0.0
            inj_bg = injection.b_e * win_scale
        else:
            inj_bg = float(b_injection(chan, d.drive.I_uA, Tj)) if d.drive.mode != "PL" else 0.0
        G = d.cavity.G if (d.cavity.enabled and not sin_mode) else 1.0
        if diode:
            # item 2 (council review 2026-09-05, round 3): inj_bg (==
            # injection.b_e * win_scale) is transport's background
            # normalized per UNFILTERED X photon (transport.py's rate_x
            # carries no t_X), while docs/rt_edge_contract.md's rho is "per
            # COLLECTED X photon" -- the CW path already divides by t_X
            # (sig = t_X*I_X + t_XX*I_XX below); the pulsed path did not
            # (rho_op 0.98896 vs the correct 0.97815 == cw_rho_op at the
            # gaasp operating point). Dividing by spec.t_x converts the
            # per-emitted-X ratio into per-collected-X, matching the CW
            # convention exactly.
            bg_per_collected = inj_bg / spec.t_x if spec.t_x > 0 else float("inf")
            # item 8: drive.b_res is already stated per collected X photon
            # (same normalization) -- added once, alongside it.
            B = params["b0"] + params["beta"] * (1.0 - S) + G * S * (bg_per_collected + d.drive.b_res)
        else:
            B = params["b0"] + params["beta"] * (1.0 - S) + inj_bg + d.drive.b_res
        rho = G * S / (G * S + B)
        # T1 mechanism library (opt-in): resolve the card's mechanism into the
        # DriveInterface at THIS junction temperature (SET pricing is T-honest)
        # and use its (mu, F_p, eta). mechanism == "" keeps the legacy fields.
        mu_use, Fp_use, eta_use = d.drive.mu, d.drive.F_p, d.drive.eta_capture
        if injection is not None:
            mu_use, eta_use = injection.mu or 0.0, injection.eta_capture
        if d.drive.mechanism:
            iface = mech_from_card(d.drive.mechanism, d.drive.mech_params,
                                   T_K=Tj, I_uA=d.drive.I_uA,
                                   mu=d.drive.mu, eta=d.drive.eta_capture)
            mu_use, Fp_use, eta_use = iface.mu, iface.F_p, iface.eta_capture
        if mu_use > 0:
            # finding 2 (peer-review-triage.md): drive.loading_model picks the
            # cap-2 loading convention explicitly -- "auto" (default)
            # reproduces today's switch (f8_g2 whenever F_p != 1.0 exactly)
            # bit-identically; "moment_matched" always uses f8_g2 (continuous
            # across F_p=1 by construction); "capped_poisson" always uses
            # f1b_g2. [A]
            use_f8 = (d.drive.loading_model == "moment_matched"
                     or (d.drive.loading_model == "auto" and Fp_use != 1.0))
            if use_f8:
                # F8b thinning: effective pump Fano factor AT the dot, after
                # dot-capture; F8: moment-matched (mu, Fano) cap-2 loading.
                F_eff = f8b_thin_fano(eta_use, Fp_use)
                if d.drive.mechanism:
                    # mechanism path: out-of-domain (mu, F) at this T is a
                    # gap, not a crash (a mechanism may drift out of the
                    # cap-2 domain as Tj moves along the curve)
                    try:
                        g2_dot = float(f8_g2(mu_use, F_eff, spec.eps))
                    except ValueError:
                        g2_dot = np.nan
                else:
                    # legacy F_p path: keep v1.1 semantics EXACTLY -- domain
                    # violations raise; the envelope legs (design_meta) do
                    # their own gap handling and verify4 asserts the raise
                    g2_dot = float(f8_g2(mu_use, F_eff, spec.eps))
            else:
                # loading_model in ("auto" with F_p == 1.0, "capped_poisson"):
                # keep the original f1b_g2 path EXACTLY -- f8_g2(mu, 1, eps)
                # is not bit-identical to f1b_g2(mu, eps) at finite mu
                # (different, equally valid cap-2 conventions; see loading.py
                # module docstring), so every pre-v1.1 result stays
                # reproducible unless loading_model or F_p is explicitly
                # changed from its default.
                g2_dot = float(f1b_g2(mu_use, spec.eps))
        else:
            g2_dot = spec.eps
        # finding 1 (peer-review-triage.md), gate-consistent counting fix
        # (pr-pkg4-fix, Opus review): drive.finite_pulse replaces the static
        # per-pulse g2_dot above with pulse_counting.pulse_g2's exact
        # state-resolved moment-hierarchy result, propagated through the
        # pump ("on") and dark ("off") windows of one pulse period -- see
        # DriveBlock.finite_pulse and pulse_counting.py. Requires EL-
        # transport (injection is not None) for the physical pump rate
        # injection.loading.r_dot; False (default) or no injection leaves
        # g2_dot AND rho untouched (legacy, bit-identical). gamma_X_ns/k_X/
        # k_XX use the EXACT same expressions as the CW branch below
        # (rm-enhanced radiative rate, cw_g2.escape_rates_from_retention at
        # this Tj) so the two opt-in diagnostics stay consistent with each
        # other.
        #
        # Whenever finite_pulse IS evaluable here, rho is ALSO always
        # rebuilt from the SAME gate-restricted counts pulse_g2 just
        # produced -- gate_fp defaults to the whole period when
        # drive.gate_ns is not set (DriveBlock.gate_ns), so signal and
        # background never again live on different windows (a prior round
        # divided a full-period signal by a gate-WIDTH background integral,
        # a 28.8% internal inconsistency at the 230 K corner). signal =
        # G*(n_X+n_XX); B = b0 + beta*(1-S) + n_bg + b_res*n_X. The b0/beta
        # terms sit exactly where the legacy static-loading rho (G*S/
        # (G*S+B) above) has them, but n_bg and b_res*n_X do NOT carry that
        # legacy formula's G*S factor (pr-pkg4-fix2, Opus review item 3: an
        # earlier version of this comment claimed they did -- false. The
        # legacy diode B multiplies (bg_per_collected + b_res) by G*S
        # because bg_per_collected is a RATIO per generated X photon that
        # needs converting to absolute count units via the collected
        # signal G*S; b_res is folded in at the same conversion for
        # convenience.) Here n_bg is already an absolute background photon
        # count integrated directly from transport's background rate (see
        # below), and b_res*n_X multiplies b_res (per-collected-X-photon,
        # same normalization as the legacy path) by pulse_g2's own
        # absolute mean_counts_x -- both are already on the same "absolute
        # count" footing as signal_fp = G*mean_counts, so neither needs a
        # further G or S factor. This instead mirrors the CW branch's own
        # bg formula below (bg = abs_bg_in_window_ns*win_scale +
        # b_res*t_X*I_X -- no G or S factor there either), which is what
        # item 1 (pr-pkg4-fix2) requires: rho_pulsed must reduce to
        # cw_rho_op as tau_dark_ns -> 0, which only happens if B_fp and the
        # CW bg share the same structure. Cards leave the cavity disabled
        # (G == 1.0 whenever d.cavity.enabled is False), so this
        # correction does not move any card's numbers today regardless of
        # which placement had been chosen.
        #
        # pr-pkg4-fix3 item 3 (Opus review); pr-pkg4-fix4 item 1 (Opus
        # review, correcting this same comment): the "reduces exactly to
        # cw_rho_op" claim below holds ONLY when the RESOLVED retention
        # params actually used above -- params["b0"] == params["beta"] ==
        # 0 -- AND the cavity is disabled (true of every shipped card) --
        # with either active the two rho definitions differ BY
        # CONSTRUCTION, not by approximation error: cw_rho_op's own bg (the
        # CW branch further below) never carries a b0/beta term at all, so
        # a nonzero b0/beta only ever inflates B_fp, pulling rho_pulsed
        # below cw_rho_op no matter how small tau_dark_ns gets; and
        # cw_rho_op's own signal (sig = t_X*I_X + t_XX*I_XX, CW branch
        # further below) never carries a G factor at all, so G != 1.0
        # (cavity enabled) scales signal_fp above but leaves cw_rho_op's
        # numerator untouched, again breaking the reduction regardless of
        # tau_dark_ns. The condition is on params["b0"]/params["beta"], NOT
        # the raw design.ret.b0/design.ret.beta fields: in ret.mode=
        # "proxy" (the default) an unset (0.0) ret.b0/ret.beta resolves
        # through the class-proxy Arrhenius fit to b0=0.01276, beta=1e-4 --
        # nonzero -- via `rp` above, so a proxy-mode card can fail this
        # reduction while its own raw RetentionBlock fields still read
        # 0.0; and in ret.mode="confinement" an explicit ret.overrides
        # entry can likewise set a nonzero b0/beta on top of the raw
        # (still-0.0) RetentionBlock fields. Exposed to callers as
        # scalars["retention_params_used"] (b0/beta only) so a guard
        # outside this module (verify_device_rt.py's _cw_reduction_guard)
        # can test the condition that actually applies here rather than
        # the raw fields. Checks E/F below skip (with an explanatory
        # message) whenever the RESOLVED b0/beta is nonzero or the cavity
        # is enabled.
        #
        # n_bg (item 1, pr-pkg4-fix2) is integrated over min(gate_fp,
        # tau_pulse_ns_val), NOT the full gate: background photons are
        # counted only while injection current actually flows, and any
        # afterglow past the pulse end is neglected [A] (see also
        # pulse_counting.py's module docstring). Before this fix the
        # background integrated over the whole gate while signal came only
        # from the ~0.1 ns pump window, giving rho_pulsed = 0.211 against
        # cw_rho_op = 0.858 at the same (230 K, gainp) operating point;
        # with the fix, rho_pulsed -> cw_rho_op exactly as tau_dark_ns ->
        # 0, for any gate_fp that still covers the whole pump window
        # (gate_fp == period or gate_fp == tau_pulse_ns + 5*tau_rad_ns
        # both qualify) -- subject to the b0/beta/cavity caveat just above.
        #
        # pr-pkg4-fix3 item 4 (physics honesty): neglecting background
        # afterglow past the pulse end is one-sided and optimistic --
        # rho_pulsed is therefore an UPPER BOUND on the true pulsed rho,
        # not a central estimate. At the FAVOURABLE corner (T_hs=230 K,
        # gamma300=6, delta_xx=8, NA=0.8, R_back=0.95, L=250 um) n_bg is
        # ~17% of B_fp; a ~1 ns background carrier lifetime (afterglow
        # decaying on that scale rather than being cut off at the pulse
        # end) would carry ~11x more background photons in the counting
        # window and move rho_pulsed from 0.858 to ~0.69, while g2_op
        # moves only from 0.9817 to 0.98826 over the same change [A] (g2 is
        # far less sensitive to the background model than rho is, since
        # g2 depends on the cascade dynamics rather than the signal-to-
        # background ratio directly). At T_hs=300 K instead (same corner
        # otherwise) n_bg/B_fp = 0.117152, rho_pulsed = 0.8662227488,
        # g2_op = 0.9976578221 -- quote both temperatures, they are not
        # interchangeable.
        finite_pulse_g2_dot = float("nan")
        finite_pulse_mean_counts = float("nan")
        finite_pulse_mean_counts_x = float("nan")
        finite_pulse_rates = _nan_fp_rates()
        finite_pulse_gate_ns_used = float("nan")
        finite_pulse_converged = False
        rho_pulsed = float("nan")
        invalid_reason = None
        if d.drive.finite_pulse and injection is not None and mu_use > 0:
            gamma_X_ns_fp = rm / d.ret.tau_rad_ns
            k_X_fp, k_XX_fp = cw_g2.escape_rates_from_retention(
                gamma_X_ns_fp, params["a_esc"], params["E_a"], params["b_p"], params["E_b"], Tj)
            r_ns_fp = injection.loading.r_dot / 1e9
            # rep_rate_hz is resolved identically to the outer-scope pulsed-
            # flux bookkeeping below (duty_eff/tau_pulse_ns_val whenever
            # drive.rep_rate_hz was not itself the explicit choice); recomputed
            # here (rather than threaded in) since `one` is a closure over the
            # same tau_pulse_ns_val/rep_explicit/duty_eff already in scope.
            rep_rate_hz_fp = (d.drive.rep_rate_hz if rep_explicit
                             else (duty_eff / (tau_pulse_ns_val * 1e-9)
                                   if tau_pulse_ns_val > 0 and duty_eff > 0 else float("nan")))
            tau_dark_ns_fp = (1e9 / rep_rate_hz_fp - tau_pulse_ns_val
                             if np.isfinite(rep_rate_hz_fp) and rep_rate_hz_fp > 0 else float("nan"))
            if (gamma_X_ns_fp > 0 and spec.t_x > 0 and r_ns_fp > 0
                    and np.isfinite(tau_dark_ns_fp) and tau_dark_ns_fp >= 0):
                period_fp = tau_pulse_ns_val + tau_dark_ns_fp
                gate_fp = d.drive.gate_ns if d.drive.gate_ns is not None else period_fp
                pc = pulse_counting.pulse_g2(
                    r_ns=r_ns_fp, gamma_X_ns=gamma_X_ns_fp, gamma_XX_ns=2.0 * gamma_X_ns_fp,
                    k_X=k_X_fp, k_XX=k_XX_fp, t_X=spec.t_x, t_XX=spec.eps * spec.t_x,
                    tau_on_ns=tau_pulse_ns_val, tau_dark_ns=tau_dark_ns_fp,
                    pump_ratio=d.drive.cw_pump_ratio, split=True, gate_ns=gate_fp)
                finite_pulse_mean_counts = pc["mean_counts"]
                finite_pulse_mean_counts_x = pc.get("mean_counts_x", float("nan"))
                finite_pulse_rates = dict(r_ns=r_ns_fp, gamma_X_ns=gamma_X_ns_fp,
                                          gamma_XX_ns=2.0 * gamma_X_ns_fp, k_X=k_X_fp, k_XX=k_XX_fp,
                                          tau_on_ns=tau_pulse_ns_val, tau_dark_ns=tau_dark_ns_fp)
                finite_pulse_converged = bool(pc["converged"])
                # finding 7 (Opus review): an unconverged periodic steady
                # state (or the mean_counts underflow pulse_g2 also reports
                # as converged=False) makes the moment-hierarchy result
                # meaningless -- nan g2_dot/rho and mark the row invalid the
                # same way the early-return operating points above do,
                # rather than reporting a number from a bad fixed point.
                # pr-pkg4-fix3 item 2 (Opus review): finite_pulse_gate_ns_used
                # is assigned ONLY here, after the convergence test succeeds,
                # so it stays NaN in the else branch below -- otherwise the
                # "finite_pulse_waveform" provenance fallback ("requested;
                # operating point invalid (<reason>)") never fires for the
                # unconverged case, since np.isfinite(gate_ns_used) would
                # already be True from a value set before this check ran.
                if finite_pulse_converged and np.isfinite(finite_pulse_mean_counts):
                    finite_pulse_gate_ns_used = gate_fp
                    finite_pulse_g2_dot = pc["g2"]
                    g2_dot = finite_pulse_g2_dot
                    n_bg = (injection.background.rate_bg_window * win_scale
                           * min(gate_fp, tau_pulse_ns_val) * 1e-9)
                    signal_fp = G * finite_pulse_mean_counts
                    B_fp = (params["b0"] + params["beta"] * (1.0 - S) + n_bg
                           + d.drive.b_res * finite_pulse_mean_counts_x)
                    rho_pulsed = signal_fp / (signal_fp + B_fp) if (signal_fp + B_fp) > 0 else float("nan")
                    rho = rho_pulsed
                else:
                    finite_pulse_g2_dot = float("nan")
                    g2_dot = float("nan")
                    rho_pulsed = float("nan")
                    rho = float("nan")
                    invalid_reason = ("finite-pulse moment-hierarchy result invalid: "
                                      "pulse_counting.pulse_g2's periodic steady state did "
                                      "not converge (or mean_counts underflowed to 0)")
        # Continuous aperture composition [A] (docs/rt_edge_contract.md
        # "Aperture assumptions"): composed AFTER loading/capture/filter are
        # already folded into g2_dot (spec.eps upstream, mu/F_p above), and
        # BEFORE the background law below -- every loss is applied exactly
        # once. Legacy (compose=False) leaves g2_dot untouched; the rounded
        # aperture_g2_penalty scalar (below, at the operating point only)
        # stays informational-only in that case, unchanged from before.
        # item 9: w_ap is the SAME collection-window concept as the filter's
        # own w above, so it carries the same auto_w_scale (1.0 = legacy).
        w_ap = gam * d.filter.auto_w_scale if (d.filter.auto_w or not d.filter.enabled) else d.filter.w
        if d.aperture.compose:
            _, lam_row = _aperture_lambda(
                _legacy_density_cm2(d.aperture.density_cm2), np.pi * (d.aperture.diameter_um / 2) ** 2,
                w_ap, d.aperture.sigma_inh, d.aperture.comp_brightness)
            g2_dot = float(_compose_aperture_g2(g2_dot, lam_row))
        else:
            lam_row = float("nan")

        # CW diagnostics (opt-in, drive.cw=True; requires EL-transport so the
        # pump rate comes from transport.loading.r_dot, never inferred from
        # the dimensionless pulsed mu). gamma_X_ns carries the same Purcell
        # "radiative enhancement" (rm) as the pulsed retention S above; k_X/
        # k_XX are the matching ABSOLUTE non-radiative rates for that same
        # enhanced gamma_X_ns (integrator.retention-consistent, see
        # cw_g2.escape_rates_from_retention). rho_cw is derived from the
        # TOTAL filtered X+XX signal, not 1/(1+b_e): b_e is a per-X ratio, so
        # naively applying it to the X+XX total would double-discount the XX
        # contribution.
        g2_cw0 = g2_cw0_raw = float("nan")
        cw_r_ns = cw_gamma_X_ns = cw_rho = float("nan")
        if d.drive.cw and injection is not None:
            gamma_X_ns = rm / d.ret.tau_rad_ns
            k_X, k_XX = cw_g2.escape_rates_from_retention(
                gamma_X_ns, params["a_esc"], params["E_a"], params["b_p"], params["E_b"], Tj)
            r_ns = injection.loading.r_dot / 1e9
            t_X, t_XX = spec.t_x, spec.eps * spec.t_x
            if r_ns > 0 and gamma_X_ns > 0 and t_X > 0:
                I_X, I_XX = cw_g2.photon_rates(r_ns, gamma_X_ns, 2.0 * gamma_X_ns,
                                               k_X, k_XX, d.drive.cw_pump_ratio)
                sig = t_X * I_X + t_XX * I_XX
                # item 3 fix (round 1, council review 2026-09-05): the old
                # `inj_bg * t_X * I_X` chained the PULSED ratio inj_bg
                # (denominator = transport's own min(r_dot,1/tau_rad)*S,
                # Background.rate_x) onto the CW rate equation's OWN X rate
                # I_X -- consistent only away from saturation, 13x too large
                # at the card operating point where the two X rates diverge.
                # Carry transport's absolute in-window background rate
                # (photons/s, the SAME physical operating point regardless
                # of pulsed/CW framing) into the CW calculation directly.
                #
                # item 5 fix (round 2, council review 2026-09-05):
                # abs_bg_in_window_ns is a BROADBAND continuum already
                # integrated over the TRANSPORT window w_transport
                # (opts['w_meV']); attenuating it a second time by the X
                # LINE's own spectral transmission t_X (round 1's formula)
                # undercounts a flat continuum by ~940x at the sweep's worst
                # corner (the line's own window transmission collapses while
                # the transport window does not -- unrelated quantities). A
                # flat continuum's collected rate instead scales linearly
                # with the actual collection (filter) window width: reuse
                # the SAME min(1, w_f/w_transport) ratio as the pulsed inj_bg
                # path above (win_scale; a no-op, 1.0, unless
                # drive.diode.w_meV explicitly overrides the transport
                # window away from the filter's own).
                abs_bg_in_window_ns = injection.background.rate_bg_window / 1e9
                # item 8 (council review 2026-09-05, round 3): drive.b_res is
                # stated per collected X photon (t_X*I_X, the same
                # normalization the pulsed path now uses) -- added once,
                # alongside the transport background, so pulsed and CW share
                # the identical residual-channel convention.
                bg = abs_bg_in_window_ns * win_scale + d.drive.b_res * t_X * I_X
                rho_cw = sig / (sig + bg) if (sig + bg) > 0 else float("nan")
                report = cw_g2.cw_report(
                    r_ns, gamma_X_ns, 2.0 * gamma_X_ns, k_X, k_XX, t_X, t_XX, rho_cw,
                    d.drive.cw_irf_fwhm_ps, tau_max_ns=d.drive.cw_tau_max_ns,
                    pump_ratio=d.drive.cw_pump_ratio, irf_shape=d.drive.cw_irf_shape)
                if d.aperture.compose and np.isfinite(lam_row):
                    tau = report["curves"]["tau"]
                    g_mix = _compose_aperture_g2(report["curves"]["g2_dot"], lam_row)
                    g_meas = cw_g2.g2_with_background(g_mix, rho_cw)
                    g_raw = (cw_g2.convolve_irf(tau, g_meas, d.drive.cw_irf_fwhm_ps,
                                                d.drive.cw_irf_shape)
                             if d.drive.cw_irf_fwhm_ps > 0 else g_meas)
                    g2_cw0 = float(np.interp(0.0, tau, g_meas))
                    g2_cw0_raw = float(np.interp(0.0, tau, g_raw))
                else:
                    g2_cw0 = report["g2_meas0"]
                    g2_cw0_raw = report["g2_raw0"]
                cw_r_ns, cw_gamma_X_ns, cw_rho = r_ns, gamma_X_ns, rho_cw
        return dict(Tj=Tj, gam=gam, eps=spec.eps, rho=rho,
                    g2=g2_from(g2_dot, rho), t_x=spec.t_x, runaway=False,
                    mu=mu_use, eta_capture=eta_use, b_e=inj_bg,
                    injection=injection, S=S, g2_cw0=g2_cw0, g2_cw0_raw=g2_cw0_raw,
                    cw_r_ns=cw_r_ns, cw_gamma_X_ns=cw_gamma_X_ns, cw_rho=cw_rho,
                    w=(w if w is not None else float("nan")),
                    n_dot_cm2_used=n_dot_cm2_used,
                    finite_pulse_g2_dot=finite_pulse_g2_dot,
                    finite_pulse_mean_counts=finite_pulse_mean_counts,
                    finite_pulse_mean_counts_x=finite_pulse_mean_counts_x,
                    finite_pulse_rates=finite_pulse_rates,
                    finite_pulse_gate_ns_used=finite_pulse_gate_ns_used,
                    finite_pulse_converged=finite_pulse_converged,
                    rho_pulsed=rho_pulsed,
                    # pr-pkg4-fix4 item 1: the RESOLVED b0/beta this
                    # operating point actually used (params["b0"]/
                    # params["beta"] above) -- NOT design.ret.b0/
                    # design.ret.beta, which can differ from these (proxy
                    # mode's class-fit fallback, confinement mode's
                    # ret.overrides). See the pr-pkg4-fix3 item 3 comment
                    # above for why the cw-reduction identity depends on
                    # THESE values, not the raw RetentionBlock fields.
                    retention_params_used={"b0": params["b0"], "beta": params["beta"]},
                    **({"invalid_reason": invalid_reason} if invalid_reason else {}))

    Ts = np.asarray(T_grid if T_grid is not None else np.linspace(4.0, 350.0, 120))
    rows = [one(T) for T in Ts]
    curves = {
        "T_hs": Ts,
        "Tj": np.array([r["Tj"] for r in rows]),
        "g2": np.array([r["g2"] for r in rows]),
        "eps": np.array([r["eps"] for r in rows]),
        "rho2": np.array([r["rho"] for r in rows]) ** 2,
        "gamma": np.array([r["gam"] for r in rows]),
    }
    if d.drive.cw:
        # T-aligned CW diagnostics (docs/rt_edge_contract.md); g2_op/g2 above
        # remain the pulsed-intrinsic headline metric (Lemma 1 note).
        curves["g2_cw0"] = np.array([r["g2_cw0"] for r in rows])
        curves["g2_cw0_raw"] = np.array([r["g2_cw0_raw"] for r in rows])

    # item 2 (pr-pkg1-fix2, runtime regression): T_grid frequently IS exactly
    # [thermal.T_hs] (e.g. the card verifier's own evaluate() call), which
    # previously recomputed the whole operating point here a second time --
    # self-heating fixed point, confinement, transport, CW cw_report -- for a
    # T already in `rows`. Reuse that row (instead of re-deriving it) for a
    # T already in `rows`; a generic T_grid essentially never lands on T_hs,
    # so this falls through to the original call there.
    #
    # pr-pkg1-fix3 item 10: match with abs(Ts - T_hs) < 1e-9 rather than
    # exact (==) float equality -- a T_grid built by arithmetic (e.g.
    # np.linspace's endpoint, or a caller's own T_hs +/- delta sweep) can
    # land a few ULPs off T_hs even when it is "the same" temperature by
    # construction, and would silently fall through to a second full
    # evaluate() at that point instead of reusing `rows`' matching row --
    # correct either way (within 1e-9 K, numerically equivalent to
    # ~1e-12 relative), but wasteful. pr-pkg1-fix4 item 5: dropped the
    # stronger "bit-identical" claim above -- a few ULPs of T offset does
    # move the self-heating fixed point by that same tiny amount, it just
    # doesn't matter at any precision this module reports.
    _hs_hits = np.flatnonzero(np.abs(Ts - d.thermal.T_hs) < 1e-9)
    op = rows[int(_hs_hits[0])] if _hs_hits.size else one(d.thermal.T_hs)
    # T_c: first heatsink temperature where g2 crosses 0.5 (above the g2 minimum)
    g2c = curves["g2"]
    Tc = np.nan
    if np.isfinite(g2c).any():
        imin = int(np.nanargmin(g2c))
        for i in range(imin, len(Ts) - 1):
            if np.isfinite(g2c[i]) and np.isfinite(g2c[i + 1]) and g2c[i] < 0.5 <= g2c[i + 1]:
                Tc = float(np.interp(0.5, [g2c[i], g2c[i + 1]], [Ts[i], Ts[i + 1]]))
                break

    _, P1, P2 = loading_probs(op["mu"])
    # item 9: same auto_w_scale as the filter's own resolved w (1.0 = legacy).
    w_ap = op["gam"] * d.filter.auto_w_scale if d.filter.auto_w or not d.filter.enabled else d.filter.w
    if np.isfinite(w_ap):
        Nw = float(n_window_competitors(
            _legacy_density_cm2(d.aperture.density_cm2), np.pi * (d.aperture.diameter_um / 2) ** 2,
            w_ap, d.aperture.sigma_inh))
        n_comp = max(int(round(Nw)), 0)
        ap_pen = float(aperture_g2([1.0] + [d.aperture.comp_brightness] * n_comp)) \
            if n_comp else 0.0
    else:
        Nw, ap_pen = float("nan"), 0.0  # runaway: no meaningful operating window
    lam_op = Nw * d.aperture.comp_brightness if d.aperture.compose and np.isfinite(Nw) else float("nan")

    gain_factor = d.cavity.G if (d.cavity.enabled and not sin_mode) else 1.0
    beta_factor = d.cavity.beta_sin if sin_mode else 1.0
    if d.drive.finite_pulse and diode is not None and np.isfinite(op["finite_pulse_mean_counts"]):
        # finding 1: pulse_counting's m1 (mean_counts) is already the
        # detected-photon count per pulse period -- escape (k_X/k_XX enter
        # the generator M directly) and the filter transmission (J =
        # t_X*gamma_X, t_XX*gamma_XX) are BOTH already folded in, so do NOT
        # multiply by t_x or S again here (that would double-count
        # retention/collection). This replaces the legacy static-loading
        # peak-area brightness (P1+P2)*t_x[*S] below exactly when
        # finite_pulse is opted in.
        brightness = float(op["finite_pulse_mean_counts"]) * gain_factor * beta_factor
    else:
        brightness = float((P1 + P2)) * (op["t_x"] if np.isfinite(op["t_x"]) else 0.0) \
            * gain_factor * beta_factor
        if diode is not None:
            # `mu` is captured-dot loading already: multiply retention and the
            # spectral window once, never eta_inj/eta_capture again.
            brightness *= op["S"]

    # emission.type="edge" (Lemma 1: brightness only -- eta_total already
    # contains beta, the NA factor, and ONE facet factor: a ray-probability
    # series with the dot at mid-ridge, eta_facet = 0.5*T*exp(-a*L/2) *
    # (1 + R_back*exp(-a*L)) / (1 - R_back*R_front*exp(-2*a*L)), R_front =
    # 1 - T, R_back=None meaning the bare cleaved facet (R_back=R_front),
    # single-pass propagation already folded into eta_facet [DR] Coldren,
    # Corzine & Masanovic, *Diode Lasers and Photonic Integrated Circuits*,
    # 2nd ed., ch. 2. None of beta/facet/cavity.beta_sin are applied a
    # second time here; see waveguide.edge_emission / facet_escape_fraction
    # for the model in use).
    edge, edge_lambda_nm, edge_err = None, float("nan"), None
    if d.emission.type == "edge" and np.isfinite(op["Tj"]):
        try:
            edge, edge_lambda_nm = _resolve_edge(d.ret, d.emission, op["Tj"])
            brightness *= edge.eta_total
        except ValueError as e:
            edge_err = str(e)
            brightness = float("nan")
    elif d.emission.type == "edge":
        edge_err = "no finite operating temperature for edge emission"
        brightness = float("nan")

    # council review 2026-09-05: injection area actually used (item 3), the
    # background's own sub-turn-on suppression (item 1), and a carrier-budget
    # closure diagnostic (item 2) -- all read off the operating-point
    # InjectionResult, never recomputed.
    inj_op = op["injection"]
    if inj_op is not None:
        injection_area_um2 = diode.area_um2
        injection_f_qfl_bg = inj_op.background.f_qfl_bg
        supply_active = inj_op.leakage.eta_inj * (d.drive.I_uA * 1e-6 / transport.Q_SI)
        ld_op = inj_op.loading
        # item 6 (council review 2026-09-05, round 3): the round-2
        # carrier_budget_closure = (r_captured+r_matrix)/supply_active was
        # identically 1 BY CONSTRUCTION (transport.py's dot_loading defines
        # r_matrix as supply_active - r_captured's complement exactly, see
        # its own docstring item 2) -- it could never fail, so it caught
        # nothing. photon_budget instead sums FOUR independently-read-off
        # channels -- dot radiative (X+XX, the S-survived share of
        # r_captured), dot non-radiative (the (1-S) escape/retention loss),
        # matrix radiative (the eta_rad_matrix share of r_matrix, the FULL
        # emission, not the b_e collection window slice), matrix
        # non-radiative -- each a genuinely separate term rather than a
        # pre-closed pair, so a future change that drops one silently is
        # visible here (verify's own self-test drops one deliberately and
        # confirms the ratio moves off 1).
        eta_rad_matrix_op = float((d.drive.diode or {}).get("eta_rad_matrix", 0.1))
        pb_dot_rad = ld_op.r_captured * op["S"]
        pb_dot_nonrad = ld_op.r_captured * (1.0 - op["S"])
        pb_matrix_rad = ld_op.r_matrix * eta_rad_matrix_op
        pb_matrix_nonrad = ld_op.r_matrix * (1.0 - eta_rad_matrix_op)
        photon_budget = ((pb_dot_rad + pb_dot_nonrad + pb_matrix_rad + pb_matrix_nonrad)
                        / supply_active if supply_active > 0 else float("nan"))
    else:
        injection_area_um2 = float("nan")
        injection_f_qfl_bg = float("nan")
        supply_active = float("nan")
        pb_dot_rad = pb_dot_nonrad = pb_matrix_rad = pb_matrix_nonrad = float("nan")
        photon_budget = float("nan")

    # item 1 (council review 2026-09-05, round 3): the pulsed collected photon
    # flux [A] converts brightness_per_pulse (photons collected per pulse) to
    # a rate via the repetition period tau_pulse_ns/duty_eff -- rep_rate_hz is
    # now either the card's own explicit drive.rep_rate_hz, or derived from
    # duty_eff/tau_pulse_ns_val (both resolved in the validation block above;
    # duty_eff already equals d.drive.duty whenever rep_rate_hz was not the
    # explicit choice, so this reduces to the legacy formula exactly). Only
    # meaningful for the EL-transport diode path, which is the only place
    # tau_pulse_ns is a physical pulse width rather than an unused legacy
    # field.
    if diode is not None:
        tau_pulse_ns = tau_pulse_ns_val
        rep_rate_hz = (d.drive.rep_rate_hz if rep_explicit
                       else (duty_eff / (tau_pulse_ns * 1e-9)
                             if tau_pulse_ns > 0 and duty_eff > 0 else float("nan")))
        collected_flux_pulsed_s = (brightness * rep_rate_hz
                                   if np.isfinite(brightness) and np.isfinite(rep_rate_hz)
                                   else float("nan"))
    else:
        tau_pulse_ns = float("nan")
        rep_rate_hz = float("nan")
        collected_flux_pulsed_s = float("nan")
    flux_measurable = bool(np.isfinite(collected_flux_pulsed_s) and collected_flux_pulsed_s >= 1e3)

    scalars = {
        "T_j_op": op["Tj"], "dT_J": op["Tj"] - d.thermal.T_hs,
        "runaway": bool(op["runaway"]),
        "gamma_op": op["gam"], "eps_op": op["eps"], "rho_op": op["rho"],
        "g2_op": op["g2"], "t_x_op": op["t_x"],
        "brightness_per_pulse": brightness,
        "T_c": Tc,
        "F_eff": float(purcell_eff(d.cavity.F_P, d.cavity.kappa, op["gam"]))
        if (d.cavity.enabled and not sin_mode and np.isfinite(op["gam"])) else np.nan,
        "N_w": Nw, "aperture_g2_penalty": ap_pen,
        "tag_chain": "[A]",  # unmeasured inputs are always in the chain today
        "linewidth_source": d.dot.linewidth,
        "retention_source": d.ret.mode,
        # pr-pkg4-fix4 item 1: the RESOLVED b0/beta this operating point
        # actually used (None whenever `one()` returned via an early-exit
        # branch before params was ever computed -- I_uA<=0, non-
        # convergence, thermal runaway). See the pr-pkg4-fix3 item 3
        # comment in the finite_pulse block above for why callers (e.g.
        # verify_device_rt.py's _cw_reduction_guard) must test THIS, not
        # design.ret.b0/design.ret.beta.
        "retention_params_used": op["retention_params_used"],
        "drive_source": d.drive.mode,
        "mu_resolved": op["mu"],
        "eta_capture_resolved": op["eta_capture"],
        "b_e_resolved": op["b_e"],
        "V_j_op": op["injection"].V_j if op["injection"] is not None else d.drive.V,
        "V_j": op["injection"].V_j if op["injection"] is not None else d.drive.V,
        "V_applied": op["injection"].V_applied if op["injection"] is not None else d.drive.V,
        "V_bi": op["injection"].V_bi if op["injection"] is not None else np.nan,
        "P_junction_W": (op["injection"].P_junction_W
                         if op["injection"] is not None else d.drive.I_uA * 1e-6 * d.drive.V),
        "eta_inj": (op["injection"].leakage.eta_inj
                    if op["injection"] is not None else np.nan),
        "T_j_transport": op["injection"].T if op["injection"] is not None else np.nan,
        "S_resolved": op["S"],
        "loading.r_dot": (op["injection"].loading.r_dot
                          if op["injection"] is not None else np.nan),
        # Pre-collection photon rates staged for integration-b (transport.py
        # native units, photons/s -- not rescaled to ns^-1 here so that these
        # stay exact pass-throughs of the verified Background dataclass).
        "background.rate_bg_window": (op["injection"].background.rate_bg_window
                                      if op["injection"] is not None else np.nan),
        "background.rate_x": (op["injection"].background.rate_x
                              if op["injection"] is not None else np.nan),
        # council review 2026-09-05 round 2: injection area used (item 3),
        # the background's own sub-turn-on suppression (item 1), and the
        # pulsed flux eligibility floor inputs (item 9).
        "injection.area_um2": injection_area_um2,
        "injection.f_qfl_bg": injection_f_qfl_bg,
        # round 3 item 6: photon_budget replaces the tautological
        # carrier_budget_closure; the four summed channels are exposed too so
        # a self-test can recompute the ratio with one deliberately dropped.
        "photon_budget": photon_budget,
        "photon_budget.dot_radiative": pb_dot_rad,
        "photon_budget.dot_nonradiative": pb_dot_nonrad,
        "photon_budget.matrix_radiative": pb_matrix_rad,
        "photon_budget.matrix_nonradiative": pb_matrix_nonrad,
        "photon_budget.supply_active": supply_active,
        # round 3 item 1: the repetition-rate bookkeeping, now explicit.
        "rep_rate_hz": rep_rate_hz,
        "tau_pulse_ns": tau_pulse_ns,
        "duty_resolved": duty_eff,
        "collected_flux_pulsed_s": collected_flux_pulsed_s,
        "flux_measurable": flux_measurable,
        # emission.type="edge" (Lemma 1: reported, never re-multiplied into
        # g2/eps/rho -- see the brightness composition above).
        "emission_type": d.emission.type,
        "edge_lambda_nm": edge_lambda_nm,
        "edge_beta": edge.beta if edge is not None else np.nan,
        "edge_T_facet": edge.T_facet if edge is not None else np.nan,
        "edge_eta_prop": edge.eta_prop if edge is not None else np.nan,
        "edge_eta_NA": edge.eta_NA if edge is not None else np.nan,
        "edge_eta_total": edge.eta_total if edge is not None else np.nan,
        "edge_n_eff": edge.n_eff if edge is not None else np.nan,
        "edge_n_g": edge.n_g if edge is not None else np.nan,
        "edge_Gamma_dot": edge.Gamma_dot if edge is not None else np.nan,
        "edge_A_mode_um2": edge.A_mode_um2 if edge is not None else np.nan,
        # continuous aperture composition (docs/rt_edge_contract.md); N_w
        # above is already the unrounded competitor count.
        "aperture_compose": d.aperture.compose,
        "aperture_lambda_op": lam_op,
        # material tracking (filter.track_material)
        "track_material": d.filter.track_material,
        "track_material_effective": track_material_effective,
        # item 9 (council review 2026-09-05, round 3): the resolved collection
        # window and the fraction of the XX line inside it -- eps == t_xx/t_x
        # by construction (spectral.epsilon), so t_xx = eps_op * t_x_op reads
        # off the SAME transmission the g2/eps chain already used, never a
        # second spectral evaluation.
        "w_resolved": op["w"],
        "xx_in_window": (op["eps"] * op["t_x"]
                         if np.isfinite(op["eps"]) and np.isfinite(op["t_x"]) else np.nan),
        # CW diagnostics (drive.cw)
        "g2_cw0": op["g2_cw0"], "g2_cw0_raw": op["g2_cw0_raw"],
        "cw_r_ns": op["cw_r_ns"], "cw_gamma_X_ns": op["cw_gamma_X_ns"],
        "cw_rho_op": op["cw_rho"],
        # finding 1/4 (peer-review-triage.md): drive.finite_pulse diagnostics
        # -- nan/nan rates whenever finite_pulse is off or not evaluable.
        "finite_pulse_g2_dot": op["finite_pulse_g2_dot"],
        "finite_pulse_mean_counts": op["finite_pulse_mean_counts"],
        "finite_pulse_mean_counts_x": op["finite_pulse_mean_counts_x"],
        "finite_pulse_r_ns": op["finite_pulse_rates"]["r_ns"],
        "finite_pulse_gamma_X_ns": op["finite_pulse_rates"]["gamma_X_ns"],
        "finite_pulse_gamma_XX_ns": op["finite_pulse_rates"]["gamma_XX_ns"],
        "finite_pulse_k_X": op["finite_pulse_rates"]["k_X"],
        "finite_pulse_k_XX": op["finite_pulse_rates"]["k_XX"],
        "finite_pulse_tau_on_ns": op["finite_pulse_rates"]["tau_on_ns"],
        "finite_pulse_tau_dark_ns": op["finite_pulse_rates"]["tau_dark_ns"],
        # pr-pkg4-fix items 1/7: the gate width actually used (None resolves
        # to the whole period -- see DriveBlock.gate_ns) and whether
        # pulse_counting.pulse_g2's periodic steady state converged (False
        # nans g2_dot/rho and marks the row invalid; see "finding 7" above).
        "finite_pulse_gate_ns_used": op["finite_pulse_gate_ns_used"],
        "finite_pulse_converged": op["finite_pulse_converged"],
        "rho_pulsed": op["rho_pulsed"],
        "invalid_reasons": ([] if np.isfinite(op["g2"])
                            else [op.get("invalid_reason", "invalid operating point")])
        + ([edge_err] if edge_err else [])
        + (["CW diagnostics requested but not evaluable at this operating point"]
           if d.drive.cw and not np.isfinite(op["g2_cw0"]) else [])
        # item 2: sub-turn-on warning, informational -- the sweep can see the
        # dot is being loaded almost entirely off the Boltzmann tail of the
        # carrier reservoirs (qV_j far below E_X) without this making the
        # operating point itself invalid (g2_op stays finite).
        + (["sub_turn_on: f_qfl < 0.05 (qV_j well below the dot transition; "
            "loading exponentially suppressed, Sze & Ng ch. 12)"]
           if (op["injection"] is not None and op["injection"].f_qfl < 0.05) else []),
        "provenance": {
            "linewidth": {"tag": "A" if d.dot.linewidth == "anchored" else "A",
                          "note": "anchored linewidth model" if d.dot.linewidth == "anchored"
                                  else "arsenide class proxy"},
            "retention": {"tag": "E" if d.ret.mode == "confinement" else "A",
                          "note": retention_note,
                          # finding 1b: the dot density _confinement_params actually
                          # resolved and forwarded to dot_levels.retention_params
                          # (NaN in proxy mode, where no confinement density is used).
                          "n_dot_cm2_used": op["n_dot_cm2_used"]},
            "drive": {"tag": "E" if diode is not None else "A",
                      "note": "transport Diode injection result" if diode is not None
                              else "legacy drive fields"},
            "mu_resolved": {"tag": "E" if diode is not None else "A",
                            "note": "transport dot loading" if diode is not None else "drive.mu"},
            "b_e_resolved": {"tag": "E" if diode is not None else "A",
                              "note": "in-window background per collected X photon"},
            "emission": {"tag": "E" if edge is not None else "A",
                        "note": ("waveguide.edge_emission on the shared inline stack"
                                 if edge is not None else "no edge collection (legacy)")},
            "aperture": {"tag": "A",
                        "note": ("continuous Poisson competitor-bath composition"
                                 if d.aperture.compose
                                 else "legacy informational aperture_g2_penalty only")},
            # item 7 (council review 2026-09-05, round 3): the note must say
            # "materials.bandgap on stack layer ..." ONLY when track_material
            # actually moves the filter window (track_material_effective,
            # computed above) -- cavity disabled + hold_window=False resolves
            # the material only to validate the name, changes nothing, and
            # must not claim otherwise.
            "tracking": {"tag": "DR" if track_material_effective else "A",
                        "note": (f"materials.bandgap on stack layer {d.filter.track_material!r}"
                                 if track_material_effective
                                 else ("track_material set but inert (cavity disabled, "
                                       "hold_window false)" if d.filter.track_material
                                       else "legacy hard-coded GaAs Varshni"))},
            "cw": {"tag": "A", "note": ("cw_g2 rate-equation model, transport-derived pump rate"
                                        if d.drive.cw else "not requested")},
            # finding 1/4 (peer-review-triage.md), item 11 (pr-pkg4-fix):
            # finite_pulse's moment hierarchy is [DR] on the existing cw_g2
            # generator; the rectangular pump waveform and the gate_ns
            # counting window are separate [A] assumptions -- two entries,
            # not one tag flipped to DR for both. Not requested (legacy
            # static loading) unless opted in.
            "finite_pulse_hierarchy": {
                "tag": "DR" if d.drive.finite_pulse else "A",
                "note": ("pulse_counting factorial-moment hierarchy, exact on the "
                         "existing cw_g2 generator (Hanschke et al., npj Quantum "
                         "Inf. 4, 43 (2018))"
                         if d.drive.finite_pulse else "not requested (legacy static loading)")},
            # item 5 (pr-pkg4-fix2): "not requested" is only true when
            # finite_pulse itself is off -- with finite_pulse on but the
            # operating point invalid (finite_pulse_gate_ns_used stays NaN
            # whenever the block above never reaches a converged result),
            # the fallback must say so was requested, not silently claim
            # legacy static loading was in effect.
            "finite_pulse_waveform": {
                "tag": "A",
                "note": (f"rectangular pump waveform; gate-restricted counting/rho at "
                         f"gate_ns={op['finite_pulse_gate_ns_used']:g}; background "
                         f"afterglow past the pulse end is neglected (one-sided, "
                         f"optimistic) so rho_pulsed is an upper bound -- see "
                         f"DriveBlock.gate_ns"
                         if d.drive.finite_pulse and np.isfinite(op["finite_pulse_gate_ns_used"])
                         else (f"requested; operating point invalid "
                               f"({op.get('invalid_reason', 'finite-pulse block did not run')})"
                               if d.drive.finite_pulse
                               else "not requested (legacy static loading)"))},
            # item 8: b_res is a residual background channel that this module
            # never invents a value for -- 0.0 (legacy) carries no anchor; a
            # non-zero value's own provenance (e.g. the 80 K Reischle 2008
            # electrical anchor) is set by the CARD, not fabricated here.
            "b_res": {"tag": "E" if d.drive.b_res != 0.0 else "A",
                     "note": ("residual background channel, applied once in the pulsed and "
                              "CW rho -- value and anchor must come from the card's own "
                              "provenance" if d.drive.b_res != 0.0
                              else "not requested (0.0, legacy)")},
            # item 9: document the auto_w operating convention explicitly.
            "filter_window": {"tag": "A",
                             "note": (f"auto_w sets the collection window to the full "
                                      f"linewidth Gamma(T_j) (F-series operating convention); "
                                      f"auto_w_scale={d.filter.auto_w_scale:g} scales it "
                                      f"(1.0 = legacy)" if d.filter.auto_w
                                      else "explicit fixed filter.w")},
        },
    }
    return {"curves": curves, "scalars": scalars}


# ------------------------------------------------------------- envelope mode (D2)

_ENVELOPE_CAP = 4096
_ENV_CURVES = ("g2", "eps", "rho2", "Tj", "gamma")
_ENV_SCALARS = ("T_c", "g2_op", "eps_op", "rho_op", "dT_J")


def _set_path(design: DeviceDesign, path: str, value) -> None:
    block_name, field_name = path.split(".", 1)
    setattr(getattr(design, block_name), field_name, value)


def _value_set(spec, mode: str) -> tuple:
    """One parameter's sample set from a ranged spec: a 2-tuple (lo, hi) under
    mode="extremes" becomes the 2-level factorial [lo, hi]; anything else
    (len != 2, or an explicit 2-point grid the caller wants used verbatim) is
    taken as given."""
    seq = tuple(spec)
    if len(seq) == 2:
        if mode != "extremes":
            raise ValueError(f"evaluate_envelope: unsupported mode {mode!r} for a "
                             "(lo, hi) pair -- only mode='extremes' is implemented")
        return seq
    return seq


def _cartesian_eval(design: DeviceDesign, paths: list, value_sets: dict, T_grid,
                    raise_if_all_dropped=True):
    """Evaluate design at every point of the cartesian product of value_sets
    (one value set per path, in `paths` order); returns (curves_list,
    scalars_list), one entry per sample. raise_if_all_dropped=False lets a
    caller (the tornado's collapsed runs) receive an all-NaN result instead
    of an error when the whole collapsed box is outside the F8 domain."""
    combos = list(itertools.product(*(value_sets[p] for p in paths))) if paths else [()]
    curves_list, scalars_list = [], []
    n_dropped = 0
    for values in combos:
        d = copy.deepcopy(design)
        for path, v in zip(paths, values):
            _set_path(d, path, v)
        try:
            res = evaluate(d, T_grid=T_grid)
        except ValueError:
            # F8 domain violation (mu < 1 - F_eff): no cap-2 loading
            # distribution with these moments exists -- same category as
            # thermal runaway, so the sample contributes a gap, not a crash
            # (Opus v1.1 finding). The point-evaluation path keeps raising.
            n_T = len(T_grid) if T_grid is not None else 120
            nanarr = np.full(n_T, np.nan)
            res = {"curves": {k: nanarr.copy() for k in _ENV_CURVES},
                   "scalars": {k: float("nan") for k in _ENV_SCALARS}}
            res["scalars"]["runaway"] = False
            n_dropped += 1
        curves_list.append(res["curves"])
        scalars_list.append(res["scalars"])
    if n_dropped == len(combos) and raise_if_all_dropped:
        raise ValueError(
            "evaluate_envelope: every sample in the ranged box violates the F8 "
            "loading domain (mu >= 1 - F_eff). Narrow the mu range or raise "
            "F_p/eta_capture.")
    return curves_list, scalars_list


def _curve_bands(curves_list: list) -> dict:
    bands = {}
    for name in _ENV_CURVES:
        stacked = np.array([c[name] for c in curves_list], dtype=float)
        # nanmin/nanmax warn (via warnings.warn, not an FP flag) on all-NaN
        # slices -- e.g. a fully-runaway box; guard explicitly (Opus D2 finding)
        any_finite = np.isfinite(stacked).any(axis=0)
        lo = np.full(stacked.shape[1], np.nan)
        hi = np.full(stacked.shape[1], np.nan)
        if any_finite.any():
            lo[any_finite] = np.nanmin(stacked[:, any_finite], axis=0)
            hi[any_finite] = np.nanmax(stacked[:, any_finite], axis=0)
        bands[name] = (lo, hi)
    return bands


def _scalar_bands(scalars_list: list) -> dict:
    out = {}
    for name in _ENV_SCALARS:
        vals = np.array([s[name] for s in scalars_list], dtype=float)
        finite = vals[np.isfinite(vals)]
        if finite.size == 0:
            out[name] = (float("nan"), float("nan"))
        else:
            out[name] = (float(finite.min()), float(finite.max()))
    return out


def evaluate_envelope(design: DeviceDesign, ranged: dict, T_grid=None,
                      mode: str = "extremes") -> dict:
    """Headless envelope evaluation: sweep `ranged` inputs instead of pinning
    them to a point, per the honesty discipline (unmeasured [A] inputs are
    swept, never averaged into a single "prediction").

    ranged: {META path ("dot.delta_xx", ...): (lo, hi) | explicit value list}.
    A 2-tuple is expanded to a 2-level factorial [lo, hi] under the default
    mode="extremes" -- HONEST CAVEAT: this only bounds the true envelope when
    the response is monotone in that parameter across the box; a non-monotone
    response can have interior extrema the 2-level factorial misses. Anything
    else (an explicit list/tuple of length != 2) is used verbatim as that
    parameter's grid, e.g. to include known interior points.

    Samples the cartesian product over every ranged path's value set (capped
    at 4096 evaluations -- raises ValueError above that). Returns:
      bands:        {curve_name: (lo_array, hi_array)} pointwise NaN-safe
                    min/max over samples, for g2/eps/rho2/Tj/gamma.
      mid:          curves dict at the all-midpoint sample (every ranged path
                    pinned to the midpoint of its value set).
      scalar_bands: {scalar_name: (lo, hi)} NaN-safe (all-NaN -> (nan, nan))
                    for T_c/g2_op/eps_op/rho_op/dT_J.
      tornado:      {path: width_reduction} -- for each ranged path, the drop
                    in the baseline scalar's band width (T_c if any sample has
                    a finite T_c, else g2_op) when that one path alone is
                    collapsed to its midpoint (others stay ranged). Larger is
                    a higher measurement priority.
      n_samples:    number of cartesian-product evaluations.

    ranged={} degenerates exactly to evaluate(design, T_grid): one sample, so
    bands/mid/scalar_bands all collapse to that single evaluate() call.
    """
    Ts = np.asarray(T_grid if T_grid is not None else np.linspace(4.0, 350.0, 120))
    paths = list(ranged)
    value_sets = {p: _value_set(ranged[p], mode) for p in paths}

    n_samples = 1
    for p in paths:
        n_samples *= len(value_sets[p])
    if n_samples > _ENVELOPE_CAP:
        raise ValueError(
            f"evaluate_envelope: {len(paths)} ranged parameters -> {n_samples} "
            f"cartesian-product evaluations, over the cap of {_ENVELOPE_CAP}. "
            "Sweep fewer parameters at once, or use coarser/explicit grids."
        )

    curves_list, scalars_list = _cartesian_eval(design, paths, value_sets, Ts)
    bands = _curve_bands(curves_list)
    scalar_bands = _scalar_bands(scalars_list)

    mid_values = {p: 0.5 * (min(value_sets[p]) + max(value_sets[p])) for p in paths}
    mid_design = copy.deepcopy(design)
    for p in paths:
        _set_path(mid_design, p, mid_values[p])
    try:
        mid = evaluate(mid_design, T_grid=Ts)["curves"]
    except ValueError:
        # midpoint itself outside the F8 loading domain (possible for an
        # explicit user range straddling the floor): NaN reference curve,
        # bands remain valid (Opus v1.1 finding, mid-path leg)
        n_T = len(Ts) if Ts is not None else 120
        mid = {k: np.full(n_T, np.nan) for k in _ENV_CURVES}

    baseline = "T_c" if np.isfinite(scalar_bands["T_c"][0]) else "g2_op"
    full_lo, full_hi = scalar_bands[baseline]
    full_width = full_hi - full_lo
    tornado = {}
    for p in paths:
        collapsed = dict(value_sets)
        collapsed[p] = (mid_values[p],)
        _, sc_list = _cartesian_eval(design, paths, collapsed, Ts,
                                     raise_if_all_dropped=False)
        c_lo, c_hi = _scalar_bands(sc_list)[baseline]
        # collapsed midpoint may sit outside the F8 domain -> width undefined
        collapsed_width = c_hi - c_lo
        tornado[p] = full_width - collapsed_width

    return {"bands": bands, "mid": mid, "scalar_bands": scalar_bands,
            "tornado": tornado, "n_samples": n_samples}
