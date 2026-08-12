"""Module D -- drive-mechanism library (tier plan Phase T1, 2026-08-12).

Every electrical driving mechanism is a GENERATOR of one interface (the
separation-theorem wall: drive touches g2 only through these variables):

    DriveInterface = { F_p, mu, eta_capture, timing, notes }
    + the F9 granularity/feasibility numbers where hardware pricing applies.

The F8/F8b/F9 physics itself lives in fsim_core.loading and is IMPORTED here,
never re-derived (three-layer / one-implementation rule). This module adds:

  * mechanism blocks: M-1 Poisson rail, M-2 quiet rail (+ finite regulation
    bandwidth [E]), M-3 pulsed/gated injection, M-4a/b/c the SET family
    (priced, not assumed), M-5a resonant-tunneling injection;
  * SET feasibility pricing (charging energy vs kT, R_Q, cycle rate);
  * WP-M2' re-excitation tier: the driven |0>-|X>-|XX> cascade under a pump
    window, solved EXACTLY for the photon-counting factorial moments via the
    moment-hierarchy ODEs (linear, closed -- no truncation), with a Gillespie
    Monte-Carlo second method in verify/verify_drive.py.

Units: energies meV, T in K, times ps unless a suffix says otherwise; F9
feasibility works in SI (capacitance F, resistance Ohm) like loading.py's F9.

Honesty tags: fano_regulated is [E] (a monotone-limits parameterization of
finite regulation bandwidth, not a derived noise theory -- same epistemic
class as loading.fano_pump). SET pricing uses the isolated-island E_C = e^2/C
convention of F9 (NOT e^2/2C; stated so numbers are comparable to the
granularity_N wall). WP-M2' counts X-line photons under an ideal filter --
binomial thinning by any filter transmission leaves g2 invariant (verified by
MC), so the re-excitation penalty composes with spectral leakage eps
downstream rather than being recomputed per filter.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.integrate import solve_ivp

from .loading import (
    E_SI,
    KB_SI,
    f8_g2_load,
    fano_pump,
    granularity_N,
    island_radius_nm,
)

R_Q_OHM = 25812.807459  # von Klitzing resistance h/e^2 -- tunnel junctions need R_T >> R_Q


# ------------------------------------------------------------------ the interface

@dataclass
class DriveInterface:
    """What any mechanism hands the F-series chain (separation theorem)."""

    mechanism: str
    F_p: float                    # pump Fano factor delivered to the junction
    mu: float                     # operating point (per-window mean loading)
    eta_capture: float = 1.0      # delivery efficiency (F8b thins F_p with this)
    timing: str = "cw"            # cw | pulsed | cycled
    feasible: bool = True         # hardware-priced mechanisms may say no
    notes: dict = field(default_factory=dict)

    def g2_load(self) -> float:
        """Loading-statistics purity factor for this interface (F8 small-mu
        factorization), computed from the DELIVERED Fano factor after F8b
        thinning: F_dot = 1 - eta (1 - F_p).

        Inherits F8's cap-2 domain wall: raises ValueError unless
        mu >= 1 - F_dot (a pump quieter than that at this mu cannot be
        represented by cap-2 moment matching -- the caller must move mu or
        accept the domain edge, never silently clamp)."""
        F_dot = 1.0 - self.eta_capture * (1.0 - self.F_p)
        return float(f8_g2_load(self.mu, F_dot))


# ------------------------------------------------------- M-2: regulation bandwidth [E]

def fano_regulated(F_p0, f_reg_Hz, f_sample_Hz):
    """[E] Finite-regulation-bandwidth correction to a quiet rail's Fano factor.

    The junction/rail feedback suppresses current noise only below its
    regulation bandwidth f_reg ~ 1/(2 pi R C). Emission samples pump noise at
    f_sample (the recapture/emission rate). Monotone-limits parameterization:

        F_eff = 1 - (1 - F_p0) * f_reg / (f_reg + f_sample)

    Limits: f_reg >> f_sample -> F_p0 (full regulation survives);
    f_reg << f_sample -> 1 (regulation too slow, Poisson restored). This is
    the minimal form of the planned pump-regulation-bandwidth extension
    behind the Zhao-magnitude open item -- a parameterization, not a theory."""
    F_p0 = np.asarray(F_p0, dtype=float)
    return 1.0 - (1.0 - F_p0) * f_reg_Hz / (f_reg_Hz + f_sample_Hz)


# ------------------------------------------------------------- SET family pricing (M-4)

def set_feasibility(T_K, C_sigma_F=None, radius_nm=None, eps_r=13.0,
                    R_T_ohm=1e6, ec_margin=10.0, f_cycle_Hz=None):
    """Price a Coulomb-blockade element (M-4a metallic / M-4b gated island /
    M-4c turnstile junction) at temperature T_K. Returns a verdict dict, not
    a design. Conventions: E_C = e^2/C_sigma (matches F9's granularity wall;
    the e^2/2C convention would halve E_C -- margin absorbs it).

    Feasibility requires BOTH:
      * E_C >= ec_margin * kT  (charge quantization vs thermal smearing)
      * R_T >= R_Q             (charge localization on the island)
    Max turnstile cycle rate is RC-limited: f_max ~ 1/(10 R_T C_sigma) [E],
    giving I_max = e f_max for one carrier per cycle."""
    if C_sigma_F is None:
        if radius_nm is None:
            raise ValueError("give C_sigma_F or radius_nm")
        # invert loading.island_radius_nm's isolated-sphere convention
        C_sigma_F = radius_nm / island_radius_nm(1.0, eps_r)
    E_C_J = E_SI**2 / C_sigma_F
    kT_J = KB_SI * T_K
    ec_over_kt = E_C_J / kT_J
    f_max = 1.0 / (10.0 * R_T_ohm * C_sigma_F)
    f_cyc = f_cycle_Hz if f_cycle_Hz is not None else f_max
    return {
        "C_sigma_F": float(C_sigma_F),
        "radius_nm": float(island_radius_nm(C_sigma_F, eps_r)),
        "E_C_meV": float(E_C_J / E_SI * 1e3),
        "EC_over_kT": float(ec_over_kt),
        "R_T_over_RQ": float(R_T_ohm / R_Q_OHM),
        "f_max_Hz": float(f_max),
        "I_max_pA": float(E_SI * f_max * 1e12),
        "N_granularity": float(granularity_N(C_sigma_F, T_K)),
        "feasible": bool(ec_over_kt >= ec_margin and R_T_ohm >= R_Q_OHM
                         and f_cyc <= f_max),
    }


# ------------------------------------------------------------------ mechanism blocks

def mech_poisson_rail(mu, eta_capture=1.0):
    """M-1 -- CW p-i-n on a voltage-class (shot-noise-limited) rail. The
    validated V-c slice: F_p = 1 exactly, g2_load = 1 identically."""
    return DriveInterface("M-1 poisson-rail", F_p=1.0, mu=mu,
                          eta_capture=eta_capture, timing="cw")


def mech_quiet_rail(mu, Z_drive_ohm, Z_junction_ohm, eta_capture=1.0,
                    f_reg_Hz=None, f_sample_Hz=None):
    """M-2 -- high-impedance (current-regulated) rail. F_p from the F9
    impedance divider [E]; optional finite regulation bandwidth [E]."""
    F_p = float(fano_pump(Z_drive_ohm, Z_junction_ohm))
    if f_reg_Hz is not None and f_sample_Hz is not None:
        F_p = float(fano_regulated(F_p, f_reg_Hz, f_sample_Hz))
    return DriveInterface("M-2 quiet-rail", F_p=F_p, mu=mu,
                          eta_capture=eta_capture, timing="cw",
                          notes={"f_reg_Hz": f_reg_Hz})


def mech_pulsed(I_uA, tau_on_ns, eta_capture, F_p_rail=1.0,
                R_drive_ohm=None, C_dep_F=None):
    """M-3 -- pulsed/gated injection. Carriers per window N_c = I tau_on / e;
    mu = eta * N_c (cap-2 truncation happens downstream in F8). If the rail
    RC is given, the usable window shrinks by the rise time t_r ~ R C [E]."""
    tau_eff_ns = tau_on_ns
    t_rise_ns = 0.0
    if R_drive_ohm is not None and C_dep_F is not None:
        t_rise_ns = R_drive_ohm * C_dep_F * 1e9
        tau_eff_ns = max(tau_on_ns - t_rise_ns, 0.0)
    n_c = I_uA * 1e-6 * tau_eff_ns * 1e-9 / E_SI
    return DriveInterface("M-3 pulsed", F_p=F_p_rail,
                          mu=float(eta_capture * n_c),
                          eta_capture=eta_capture, timing="pulsed",
                          notes={"carriers_per_window": float(n_c),
                                 "t_rise_ns": float(t_rise_ns),
                                 "tau_eff_ns": float(tau_eff_ns)})


def mech_set(kind, T_K, eps_cycle=0.01, mu=1.0, **feas_kw):
    """M-4 -- the SET family, priced. kind in {'metallic','gated','turnstile'}.

    All three are the mu -> 1, F_p -> 0 corner of F8; what differs is the
    hardware wall. eps_cycle is the per-cycle error (missed/extra carrier:
    co-tunneling, thermal activation, missed capture) -- for a cycled
    single-carrier source the delivered Fano factor is F_p ~ eps_cycle [E]
    (rare uncorrelated errors on a deterministic stream). The turnstile runs
    mu = 1 by construction (mu argument ignored); series-regulator kinds
    keep mu as the free operating knob it is on any rail."""
    if kind not in ("metallic", "gated", "turnstile"):
        raise ValueError(f"unknown SET kind: {kind}")
    feas = set_feasibility(T_K, **feas_kw)
    mu = 1.0 if kind == "turnstile" else mu
    F_p = float(eps_cycle) if feas["feasible"] else 1.0
    return DriveInterface(f"M-4 SET/{kind}", F_p=F_p, mu=mu,
                          eta_capture=1.0,
                          timing="cycled" if kind == "turnstile" else "cw",
                          feasible=feas["feasible"], notes=feas)


def mech_rti(mu, eta_capture_base, eta_boost=2.0, dg_inj_factor=0.5,
             F_p_rail=1.0):
    """M-5a -- resonant-tunneling injection (Kitamura-class delivery hardware).
    Prices as: capture efficiency boosted (carrier-selective delivery),
    injection broadening reduced (dg_inj scaled by dg_inj_factor < 1), at a
    series-voltage cost carried by SIM-C. Both factors are [A] until the
    EL-vs-PL measurement (#2) lands."""
    eta = min(eta_capture_base * eta_boost, 1.0)
    return DriveInterface("M-5a RTI", F_p=F_p_rail, mu=mu, eta_capture=eta,
                          timing="cw", notes={"dg_inj_factor": dg_inj_factor})


def mech_from_card(name, params=None, *, T_K=300.0, I_uA=10.0, mu=0.5,
                   eta=1.0):
    """Card-facing factory: resolve a mechanism name + params dict (the
    device card's drive.mechanism / drive.mech_params) into a DriveInterface.
    T_K / I_uA / mu / eta come from the evaluation context (junction
    temperature, drive current, card loading/capture) and are overridable
    per-mechanism through params. Unknown names raise (no silent fallback)."""
    p = dict(params or {})
    mu = p.pop("mu", mu)
    eta = p.pop("eta_capture", eta)
    if name in ("m1", "poisson-rail"):
        return mech_poisson_rail(mu=mu, eta_capture=eta)
    if name in ("m2", "quiet-rail"):
        return mech_quiet_rail(mu=mu, eta_capture=eta,
                               Z_drive_ohm=p.pop("Z_drive_ohm", 1e5),
                               Z_junction_ohm=p.pop("Z_junction_ohm", 1e3),
                               f_reg_Hz=p.pop("f_reg_Hz", None),
                               f_sample_Hz=p.pop("f_sample_Hz", None))
    if name in ("m3", "pulsed"):
        return mech_pulsed(I_uA=p.pop("I_uA", I_uA),
                           tau_on_ns=p.pop("tau_on_ns", 1.0),
                           eta_capture=eta,
                           F_p_rail=p.pop("F_p_rail", 1.0),
                           R_drive_ohm=p.pop("R_drive_ohm", None),
                           C_dep_F=p.pop("C_dep_F", None))
    if name in ("m4a", "set-metallic", "m4b", "set-gated", "m4c",
                "set-turnstile"):
        kind = {"m4a": "metallic", "set-metallic": "metallic",
                "m4b": "gated", "set-gated": "gated",
                "m4c": "turnstile", "set-turnstile": "turnstile"}[name]
        if "C_sigma_F" not in p and "radius_nm" not in p:
            p["radius_nm"] = 20.0    # representative few-tens-nm island [A]
        return mech_set(kind, T_K=p.pop("T_K", T_K), mu=mu,
                        eps_cycle=p.pop("eps_cycle", 0.01), **p)
    if name in ("m5a", "rti"):
        return mech_rti(mu=mu, eta_capture_base=eta,
                        eta_boost=p.pop("eta_boost", 2.0),
                        dg_inj_factor=p.pop("dg_inj_factor", 0.5),
                        F_p_rail=p.pop("F_p_rail", 1.0))
    raise ValueError(f"unknown drive mechanism: {name!r}")


def roster():
    """The mechanism roster with representative default instantiations --
    for reporting/GUI enumeration, not for physics."""
    return [
        mech_poisson_rail(mu=0.4),
        mech_quiet_rail(mu=0.7, Z_drive_ohm=3e3, Z_junction_ohm=7e3),
        mech_pulsed(I_uA=10.0, tau_on_ns=1.0, eta_capture=1e-5),
        mech_set("metallic", T_K=300.0, radius_nm=20.0),
        mech_set("gated", T_K=77.0, radius_nm=20.0),
        mech_set("turnstile", T_K=4.0, radius_nm=20.0, R_T_ohm=2e5),
        mech_rti(mu=0.4, eta_capture_base=0.05),
    ]


# ==================================================== WP-M2': re-excitation tier
#
# Driven cascade |0> --p--> |X> --p--> |XX>, decays |XX> --g_B--> |X> (XX
# photon, filtered out) and |X> --g_X--> |0> (X photon, COUNTED). Pump rate
# p(t) = p for t in [0, tau_on], 0 after; evolve to t_end >> lifetimes.
#
# Exact counting moments via the moment hierarchy: with occupation
# probabilities P_s, first-moment densities F_s = E[m 1_s], and second-
# factorial-moment densities S_s = E[m(m-1) 1_s] (s in {0, X, B}), the
# Markov process gives a CLOSED linear ODE system (each X emission increments
# m by exactly 1 as the state jumps X -> 0):
#
#   dP0/dt = -p P0 + gX PX               dPX/dt = p P0 - (gX + p) PX + gB PB
#   dPB/dt = p PX - gB PB
#   dF0/dt = -p F0 + gX FX + gX PX      (E[(m+1) 1_0] inflow = gX (FX + PX))
#   dFX/dt = p F0 - (gX + p) FX + gB FB
#   dFB/dt = p FX - gB FB
#   dS0/dt = -p S0 + gX SX + 2 gX FX    (E[(m+1)m 1_0] inflow = gX (SX + 2 FX))
#   dSX/dt = p S0 - (gX + p) SX + gB SB
#   dSB/dt = p SX - gB SB
#
# g2_reexc = sum(S) / sum(F)^2 at t_end. No truncation anywhere -- this is
# exact for the model. Gillespie MC in verify_drive.py is the second method.

def reexc_g2(p_per_ps, tau_on_ps, tau_x_ps, tau_xx_ps=None, mu0=0.0,
             t_end_factor=12.0, rtol=1e-10, atol=1e-14):
    """WP-M2' re-excitation g2(0) of the X line under a rectangular pump
    window, exact via the counting-moment ODEs (module comment above).

    p_per_ps: pump/capture rate during the window (both 0->X and X->XX).
    tau_on_ps: window length. tau_x_ps / tau_xx_ps: X / XX lifetimes
    (tau_xx defaults to tau_x/2, the standard cascade ratio [DR class]).
    mu0: probability the dot is ALREADY loaded with one exciton at t=0
    (models a reservoir/prior-window handoff; 0 = start empty).

    Returns dict {g2, mean_photons, m2_factorial}. Limits (verify-gated):
    weak pump (p -> 0 at fixed window) gives a FINITE g2 plateau set by the
    window/lifetime geometry -- both <m(m-1)> and <m>^2 scale as p^2, so
    re-excitation g2 is pump-power-independent at weak drive (the design
    lever is the WINDOW, not the power); tau_on -> 0 drives g2 -> 0 (short
    pulses are the fix, the standard pulsed-source result); tau_on = 0 with
    mu0 = 1 gives exactly one photon, g2 = 0; p*tau_on >> 1 with
    gX*tau_on >> 1 approaches Poisson-like g2 ~ 1."""
    if tau_xx_ps is None:
        tau_xx_ps = tau_x_ps / 2.0
    gX, gB = 1.0 / tau_x_ps, 1.0 / tau_xx_ps

    def rhs(t, y, p):
        P0, PX, PB, F0, FX, FB, S0, SX, SB = y
        return [
            -p * P0 + gX * PX,
            p * P0 - (gX + p) * PX + gB * PB,
            p * PX - gB * PB,
            -p * F0 + gX * FX + gX * PX,
            p * F0 - (gX + p) * FX + gB * FB,
            p * FX - gB * FB,
            -p * S0 + gX * SX + 2.0 * gX * FX,
            p * S0 - (gX + p) * SX + gB * SB,
            p * SX - gB * SB,
        ]

    y0 = [1.0 - mu0, mu0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]
    if tau_on_ps > 0:
        sol = solve_ivp(rhs, (0.0, tau_on_ps), y0, args=(p_per_ps,),
                        rtol=rtol, atol=atol, method="LSODA")
        y0 = sol.y[:, -1]
    t_tail = t_end_factor * max(tau_x_ps, tau_xx_ps)
    sol = solve_ivp(rhs, (0.0, t_tail), y0, args=(0.0,),
                    rtol=rtol, atol=atol, method="LSODA")
    y = sol.y[:, -1]
    mean = y[3] + y[4] + y[5]
    m2f = y[6] + y[7] + y[8]
    g2 = float(m2f / mean**2) if mean > 0 else 0.0
    return {"g2": g2, "mean_photons": float(mean),
            "m2_factorial": float(m2f)}


def mc_reexc_g2(p_per_ps, tau_on_ps, tau_x_ps, tau_xx_ps=None, mu0=0.0,
                t_filter=1.0, n_pulses=200_000, seed=0):
    """Gillespie MC second method for reexc_g2, with optional binomial
    thinning of counted X photons by t_filter (checks the thinning-invariance
    claim: g2 is independent of t_filter). Returns (g2, stderr)."""
    if tau_xx_ps is None:
        tau_xx_ps = tau_x_ps / 2.0
    gX, gB = 1.0 / tau_x_ps, 1.0 / tau_xx_ps
    rng = np.random.default_rng(seed)
    counts = np.zeros(n_pulses, dtype=np.int64)
    for i in range(n_pulses):
        t, state, m = 0.0, (1 if rng.random() < mu0 else 0), 0
        # state: 0 = ground, 1 = X, 2 = XX
        while True:
            p = p_per_ps if t < tau_on_ps else 0.0
            rates = ((p, 0.0) if state == 0 else
                     (p, gX) if state == 1 else
                     (0.0, gB))
            total = rates[0] + rates[1]
            if total == 0.0:
                break
            dt = rng.exponential(1.0 / total)
            if t < tau_on_ps <= t + dt and p > 0.0:
                # pump switches off inside this waiting time: redraw from the
                # switch point with the pump-off rates (thinning-correct)
                t = tau_on_ps
                continue
            t += dt
            if rng.random() < rates[0] / total:
                state += 1                       # excitation 0->X or X->XX
            else:
                state -= 1                       # decay XX->X or X->0
                if state == 0 and rng.random() < t_filter:
                    m += 1                       # counted X photon
        counts[i] = m
    mean = counts.mean()
    pairs = counts * (counts - 1)
    g2 = float(pairs.mean() / mean**2) if mean > 0 else 0.0
    se = float(pairs.std() / np.sqrt(n_pulses) / mean**2) if mean > 0 else 0.0
    return g2, se
