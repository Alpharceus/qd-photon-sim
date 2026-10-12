"""Polaron solver of the composable engine (P9): tier P1 of the physics brief 2.2.

A two-level emitter (G, X) driven by one or several pulses (complex Rabi frequency in the frame of pulse 1,
so SUPER and dichromatic drives are the same code) and coupled to bulk LA phonons through the independent
boson model. Tier P1 = weak-coupling phonon relaxation in the INSTANTANEOUS DRESSED BASIS of H(t):

    H(t)  = Delta |X><X| + (Omega_r(t)/2) |X><G| + h.c.,          Omega_r = B Omega,  B = sqrt(Z(T))
    Lambda = sqrt(Delta^2 + |Omega_r|^2),   sin(2 theta) = |Omega_r|/Lambda
    Gamma_down = (pi/2) sin^2(2 theta) J(Lambda) [n(Lambda) + 1]   (upper dressed -> lower dressed)
    Gamma_up   = (pi/2) sin^2(2 theta) J(Lambda)  n(Lambda)
    L(t)  = L_0 - i[H_drive, .] + Gamma_down D[|lo><hi|] + Gamma_up D[|hi><lo|] ,

with |hi>, |lo> the eigenvectors of H(t) RE-ROTATED AT EVERY STEP (brief trap 12; a fixed-basis Lindblad with
constant phonon rates gets the detuning dependence wrong) and L_0 the undriven Liouvillian (radiative Gamma,
photonless k, pure dephasing gamma*, detuning Delta) of resonant_source.TwoLevel [DR].

Derivation of the rate [DR]: golden rule for H_int = |X><X| sum_q g_q (b_q + b_q^+), J(w) = sum_q g_q^2
delta(w - w_q); the matrix element of |X><X| between the dressed states is |<hi|P_X|lo>|^2 = sin^2(2 theta)/4,
so Gamma_down = 2 pi J(Lambda)(n+1) sin^2(2 theta)/4 = (pi/2) sin^2(2 theta) J (n+1). Detailed balance
Gamma_up/Gamma_down = n/(n+1) = exp(-hbar Lambda/kT) holds exactly. On resonance (theta = pi/4) the Rabi
oscillation is damped at (pi/2) J(Omega) coth(hbar Omega/2kT) [E: Ramsay et al., PRL 104, 017402 (2010);
Nazir & McCutcheon, J. Phys. Condens. Matter 28, 103002 (2016)].

Spectral density and units. qd_gf.spectral_density(E) = alpha_meV E^3 F2(E) is in meV on the energy axis, so
J(omega) [rad/ps] = spectral_density(hbar omega)/hbar_ps and the rate in 1/ns is 1e3 x that (the Huang-Rhys
integral int dE J/E^2 = int d(omega) J(omega)/omega^2 is the same number). B = sqrt(Z) uses qd_gf.zpl_weight.

Time grid. Fine propagation grid: dt <= DT_FACTOR/rate_max with rate_max = |Delta| + sum_k |Omega_r,k| +
max(beat + chirp) (a bound on |H|; brief 2.1 asks dt <= 0.05/Lambda_max, Lambda <= rate_max/2 here) and
dt <= sigma/8 for Gaussians (resonant_source). Counting hierarchy (g2, mean photons, periodic fixed point),
inversion, and the Grange HOM integral are the same machinery as resonant_source, written for a generic list
of step propagators (the HOM quadrature groups m fine steps per quadrature step, so stiff pulses do not need
an O(n^2) matrix on the fine grid).

Limitations (reported in SolverResult.reasons):
  * P1 is weak coupling. Quilter-class pulses are not deep in that regime; the full polaron master equation
    (tier P2) is not implemented.
  * Multi-colour drives: the instantaneous basis of H(t) rotates with the beat phase; the P1 rates for a
    fast beat (|delta w| >~ Lambda) are an approximation (a Floquet-Markov treatment would be needed) [A].
  * I_single is the ZPL-coherence overlap of the single-photon wavepacket; the Z^2 / sideband bookkeeping
    stays in observables.py (not applied here, no double counting).
  * The polaron (Franck-Condon) energy shift is absorbed into the line position; detunings are measured
    from the ZPL.

Drive-induced dephasing (devlib-fx3b, Koong-type excitation-induced dephasing) [E]. Optional, default OFF
(plan.inputs["eid_ps"] = 0 is bit-identical): the pure-dephasing rate rises with the instantaneous drive intensity,
    gamma_EID(t) = c_EID |Omega(t)|^2     (c_EID = eid_ps x 1e-3 ns; Omega in rad/ns, so gamma_EID in 1/ns),
added to gamma* as a D[sqrt(gamma_EID) |X><X|] term of L(t) (repo convention: coherence decay gamma_EID/2). It uses the
BARE Rabi frequency (the photon-field intensity), not the phonon-renormalised B Omega. No coefficient comes from a
paper here: the class value is [E] (Ramsay et al., PRL 104, 017402 (2010) report damping that grows with the pulse
intensity, not a universal c), so any card that enables it must say where c comes from.

Optional trion branch (devlib-fx3b, Vajner-type excess-charge Rabi damping) [E]: `apply_trion` scales a result by a
per-pulse probability p_trion(area) = p_trion sin^2(area/2) that the excitation ends in a dark charged state
(default 0 = bit-identical); see its docstring.

Needs the temperature: plan.inputs["T_K"] (the engine adds it).
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Optional

import numpy as np
from scipy.linalg import expm

from .. import lindblad as lb
from .. import qd_gf
from .. import resonant_source as rs
from ..spectral import KB
from .registry import Solver, register
from .types import SolverResult

NAN = float("nan")
DT_FACTOR = 0.1               # [A] fine step <= DT_FACTOR / rate_max (rate_max >= 2 Lambda_max)
HOM_NODE_CAP = 2048           # [A] cap on HOM quadrature nodes (memory O(n^2))
HOM_TOL = 1e-6                # [A] same as solvers._HOM_TOL
HOM_DOUBLINGS = 5
_ALLOWED_NEEDS = {"phonons", "two_colour", "td_detuning"}
_HBAR_PS = qd_gf.HBAR          # meV ps
_SX = np.array([[0, 1], [1, 0]], dtype=complex)
_SY = np.array([[0, -1j], [1j, 0]], dtype=complex)      # basis (G, X)
_LSX = lb.hamiltonian_super(0.5 * _SX)
_LSY = lb.hamiltonian_super(0.5 * _SY)
_LDET = lb.hamiltonian_super(np.array([[0, 0], [0, 1]], dtype=complex))    # d L / d(det_X): H = det |X><X|
_IDX_PX = 3                    # vec index of rho[X, X] (column stacking, dim 2)


_PX = np.array([[0, 0], [0, 1]], dtype=complex)
_DEPH_X = lb.dissipator(_PX)          # D[|X><X|]: coherence decay 1/2 per unit rate (repo gamma* convention)


def trion_probability(p_trion, area_rad):
    """p_trion(area) = p_trion sin^2(area/2) [E]: the chance that a pulse of area `area_rad` leaves the dot in a dark
    trion (excess charge) state, proportional to the ideal two-level inversion sin^2(area/2) it creates; p_trion is the
    value at a pi pulse. No paper gives this function (Vajner 2024 only gives the reduced Rabi-oscillation fidelity
    0.82 vs 0.98 simulated); it is the minimal, parameter-light stand-in."""
    return float(p_trion) * math.sin(0.5 * float(area_rad)) ** 2


def apply_trion(res, p):
    """Dark-trion branch applied to a SolverResult: with probability p per pulse the excitation is lost to a dark charged
    state, so every emission/preparation measure scales with (1 - p) and the normalised g2 = m2/m1^2 with 1/(1 - p)
    (independent pulses, a dark pulse contributes no photon). The HOM overlap I = int int |G1|^2 / N^2 is unchanged
    (|G1|^2 ~ (1 - p)^2, N^2 ~ (1 - p)^2). p = 0 returns `res` itself (bit-identical)."""
    p = float(p)
    if not p > 0.0:
        return res
    if not p < 1.0:
        raise ValueError(f"p_trion must be < 1, got {p!r}")
    q = 1.0 - p
    ex = dict(res.extra)
    for k in ("eta_prep", "n_exc", "P_X_end", "photons_in_window", "n_XX", "P_X_max"):
        if k in ex and isinstance(ex[k], (int, float)):
            ex[k] = float(ex[k]) * q
    ex["p_trion_pulse"] = p
    ex["notes"] = list(ex.get("notes", [])) + [f"dark trion branch: p = {p:.4g} per pulse removes the excitation "
                                               "(n_ph, inversion x (1 - p); g2 / (1 - p); HOM unchanged) [E]"]
    g2 = res.g2_0 / q if math.isfinite(res.g2_0) else res.g2_0
    inv = None if res.inversion is None else float(res.inversion) * q
    return SolverResult(g2_0=g2, n_ph=res.n_ph * q, I_single=res.I_single, inversion=inv, extra=ex,
                        valid=res.valid, reasons=res.reasons)


# ---------------------------------------------------------------- phonon model
class PhononModel:
    """Spectral density, Bose occupation and dressed-state rates for one (PhononParams, T) pair."""

    def __init__(self, ph, T_K: float):
        self.ph = ph
        self.T_K = float(T_K)
        if ph is None:
            self.alpha_ps2 = 0.0
            self.Z = 1.0
        else:
            self.alpha_ps2 = float(qd_gf.coupling_alpha_ps2(ph))
            self.Z = float(qd_gf.zpl_weight(ph, self.T_K)) if self.alpha_ps2 > 0.0 else 1.0
        self.B = math.sqrt(self.Z)

    @property
    def active(self) -> bool:
        return self.ph is not None and self.alpha_ps2 > 0.0

    def J_rate_ns(self, lam_ns):
        """Spectral density J(omega = lam) in 1/ns (rad/ns x ...): spectral_density(hbar lam)/hbar_ps x 1e3."""
        lam = np.asarray(lam_ns, dtype=float)
        if not self.active:
            return np.zeros_like(lam)
        E = lb.HBAR_MEV_NS * lam
        return np.asarray(qd_gf.spectral_density(np.atleast_1d(E), self.ph), dtype=float).reshape(lam.shape) \
            / _HBAR_PS * 1e3

    def occupation(self, lam_ns):
        """Bose occupation n(hbar lam) at T."""
        lam = np.asarray(lam_ns, dtype=float)
        if self.T_K <= 0.0:
            return np.zeros_like(lam)
        x = np.minimum(lb.HBAR_MEV_NS * lam / (KB * self.T_K), 700.0)
        x = np.maximum(x, 1e-300)
        return 1.0 / np.expm1(x)

    def rates(self, lam_ns, mix):
        """(Gamma_down, Gamma_up) in 1/ns: 2 pi mix J (n+1), 2 pi mix J n with mix = |<hi|P_X|lo>|^2."""
        J = self.J_rate_ns(lam_ns)
        n = self.occupation(lam_ns)
        c = 2.0 * math.pi * np.asarray(mix, dtype=float) * J
        return c * (n + 1.0), c * n


# ---------------------------------------------------------------- pulse train
@dataclass(frozen=True)
class _Comp:
    shape: str
    dur_ns: float
    sigma_ns: float
    area: float
    amp: float          # peak Omega (rect: area/dur; gauss: area/(sigma sqrt 2 pi))
    tc: float           # centre (ns)
    dw: float           # laser frequency minus the frame (pulse 0) laser, rad/ns
    phase: float
    kappa: float        # laser chirp d(omega_L)/dt, rad/ns^2
    halfwin: float


class PulseTrain:
    """Complex Rabi frequency Omega(t) in the frame of pulse 0: sum_k Omega_k(t) exp(-i[dw_k (t - tc_k)
    + kappa_k (t - tc_k)^2/2 + phi_k]); window W = 2 max_k(|t0_k| + half window), centred at W/2."""

    def __init__(self, pulses, det_profile=None):
        if not pulses:
            raise ValueError("PulseTrain needs at least one pulse")
        det0 = float(lb.mev_to_rate(pulses[0].detuning_meV))
        # devlib-p10e Stark gate: det_X(t) = det_off + (det_on - det_off) g(t), g = RC-smoothed boxcar centred in the
        # window (width t_gate, time constant t_r; t_r = 0 sharp). det_off (REPO det_X, meV) is the value outside
        # the gate AND after the window (the dark-window Liouvillian L0 uses it); None = constant det0 (bit-identical).
        self.det_profile = None if det_profile is None else tuple(float(x) for x in det_profile)
        self._det_off = self._det_on = self._t_gate = self._t_r = 0.0
        if self.det_profile is not None:
            off, on, tg, tr = self.det_profile
            if not (math.isfinite(sum(self.det_profile)) and tg >= 0.0 and tr >= 0.0):
                raise ValueError("det_profile = (det_off_meV, det_on_meV, t_gate_ps >= 0, t_r_ps >= 0)")
            det0 = float(lb.mev_to_rate(off))
            self._det_off, self._det_on = det0, float(lb.mev_to_rate(on))
            self._t_gate, self._t_r = tg * 1e-3, tr * 1e-3
        half = []
        for p in pulses:
            dur = p.fwhm_ps * 1e-3
            half.append(abs(p.t0_ps) * 1e-3 + 0.5 * rs.pulse_window(p.shape, dur))
        self.W = 2.0 * max(half)
        comps = []
        for p in pulses:
            dur = p.fwhm_ps * 1e-3
            sig = rs.gaussian_sigma(dur)
            amp = p.area_rad / dur if p.shape == "rect" else p.area_rad / (sig * math.sqrt(2 * math.pi))
            comps.append(_Comp(
                shape=p.shape, dur_ns=dur, sigma_ns=sig, area=float(p.area_rad), amp=float(amp),
                tc=0.5 * self.W + p.t0_ps * 1e-3,
                dw=det0 - float(lb.mev_to_rate(p.detuning_meV)), phase=float(p.phase_rad),
                kappa=float(lb.mev_to_rate(p.chirp_meV_ps)) * 1e3,
                halfwin=0.5 * rs.pulse_window(p.shape, dur)))
        self.comps = tuple(comps)
        self.det0_ns = det0

    def delta_t(self, t):
        """det_X(t) (1/ns) at the times t: constant det0 without a det_profile."""
        t = np.asarray(t, dtype=float)
        if self.det_profile is None:
            return np.full(t.shape, self.det0_ns)
        ts, te = 0.5 * self.W - 0.5 * self._t_gate, 0.5 * self.W + 0.5 * self._t_gate
        if self._t_r <= 0.0:
            g = np.where((t >= ts) & (t < te), 1.0, 0.0)
        else:
            ge = 1.0 - math.exp(-self._t_gate / self._t_r)
            g = np.where(t < ts, 0.0, np.where(t < te, 1.0 - np.exp(-(t - ts) / self._t_r),
                                               ge * np.exp(-(t - te) / self._t_r)))
        return self._det_off + (self._det_on - self._det_off) * g

    def delta_cell(self, t_mid, dt, n_sub=8):
        """Cell average of det_X(t) over [t_mid - dt/2, t_mid + dt/2] (edge-weighted for a sharp gate)."""
        t_mid = np.asarray(t_mid, dtype=float)
        off = (np.arange(n_sub) + 0.5) / n_sub - 0.5
        return self.delta_t(t_mid[:, None] + off[None, :] * dt).mean(axis=1)

    def omega(self, t):
        t = np.asarray(t, dtype=float)
        out = np.zeros(t.shape, dtype=complex)
        for c in self.comps:
            x = t - c.tc
            if c.shape == "rect":
                env = np.where(np.abs(x) < 0.5 * c.dur_ns, c.amp, 0.0)
            else:
                env = np.where(np.abs(x) <= rs.GAUSS_WINDOW_SIGMA * c.sigma_ns,
                               c.amp * np.exp(-0.5 * (x / c.sigma_ns) ** 2), 0.0)
            out = out + env * np.exp(-1j * (c.dw * x + 0.5 * c.kappa * x * x + c.phase))
        return out

    def rate_max(self, B: float) -> float:
        """Upper bound on |H| (1/ns): |Delta| + sum |Omega_r| + max (beat + chirp sweep)."""
        sweep = max(abs(c.dw) + abs(c.kappa) * 0.5 * self.W for c in self.comps)
        dmax = abs(self.det0_ns) if self.det_profile is None else max(abs(self._det_off), abs(self._det_on))
        return dmax + B * sum(abs(c.amp) for c in self.comps) + sweep

    def sigma_min(self) -> Optional[float]:
        s = [c.sigma_ns for c in self.comps if c.shape == "gauss" and c.area != 0.0]
        return min(s) if s else None


# ---------------------------------------------------------------- generator
class Generator:
    """L(t) of the phonon-dressed driven two-level emitter on a grid of step midpoints."""

    def __init__(self, train: PulseTrain, tl: rs.TwoLevel, pm: PhononModel, eid_ns: float = 0.0):
        self.train = train
        self.pm = pm
        self.eid_ns = float(eid_ns)          # c_EID in ns: gamma_EID(t) = eid_ns |Omega(t)|^2 (0 = off, bit-identical)
        s0 = tl.system(0.0)
        self.L0 = np.asarray(s0.L, dtype=complex)
        self.J = lb.counting_jump(s0)
        self.A = np.sqrt(tl.gamma_ns) * s0.ops["c_X"]
        self.W = train.W
        self.delta = train.det0_ns
        self._cache = {}

    def dressed(self, t_mid, delta=None):
        """(Omega_r, Lambda, mix, vecs) at the given times. vecs[:, :, 0] = |lo>, [:, :, 1] = |hi>.
        delta = det_X at those times (default: the constant det0, or the det_profile value at t_mid)."""
        om = self.pm.B * self.train.omega(t_mid)
        n = om.size
        H = np.zeros((n, 2, 2), dtype=complex)
        if delta is None:
            delta = self.delta if self.train.det_profile is None else self.train.delta_t(t_mid)
        H[:, 1, 1] = delta
        H[:, 1, 0] = 0.5 * om
        H[:, 0, 1] = 0.5 * np.conj(om)
        ev, vec = np.linalg.eigh(H)
        lam = ev[:, 1] - ev[:, 0]
        mix = np.abs(vec[:, 1, 0] * np.conj(vec[:, 1, 1])) ** 2         # |<hi|P_X|lo>|^2
        return om, lam, mix, vec

    def step_L(self, t_mid, dt=None):
        """Liouvillians (n, 4, 4) at the midpoints t_mid (dt = step width, used to cell-average a det_profile)."""
        t_mid = np.asarray(t_mid, dtype=float)
        if self.train.det_profile is None:
            om, lam, mix, vec = self.dressed(t_mid)
        else:
            dcell = self.train.delta_cell(t_mid, dt) if dt is not None else self.train.delta_t(t_mid)
            om, lam, mix, vec = self.dressed(t_mid, delta=dcell)
        L = np.broadcast_to(self.L0, (om.size, 4, 4)).copy()
        L += om.real[:, None, None] * _LSX + om.imag[:, None, None] * _LSY
        if self.train.det_profile is not None:
            L += (dcell - self.delta)[:, None, None] * _LDET
        if self.eid_ns > 0.0:
            g_eid = self.eid_ns * (np.abs(om) / self.pm.B) ** 2            # bare |Omega|^2 (om = B Omega)
            L += g_eid[:, None, None] * _DEPH_X
        if self.pm.active:
            ok = lam > 1e-9 * max(1.0, abs(self.delta))
            gd = np.zeros(om.size)
            gu = np.zeros(om.size)
            gd[ok], gu[ok] = self.pm.rates(lam[ok], mix[ok])
            for i in np.nonzero(ok & (mix > 1e-30))[0]:
                lo, hi = vec[i, :, 0], vec[i, :, 1]
                down = np.outer(lo, np.conj(hi))
                up = np.outer(hi, np.conj(lo))
                L[i] += gd[i] * lb.dissipator(down) + gu[i] * lb.dissipator(up)
        return L

    def grid(self, n):
        """(edges, L (n,4,4), U (n,4,4)) on n equal steps over the window; cached."""
        if n not in self._cache:
            edges = np.linspace(0.0, self.W, n + 1)
            mid = 0.5 * (edges[:-1] + edges[1:])
            dt = np.diff(edges)
            L = self.step_L(mid, dt[0] if self.train.det_profile is not None else None)
            U = np.array([expm(L[i] * dt[i]) for i in range(n)])
            self._cache[n] = (edges, L, U)
        return self._cache[n]


# ---------------------------------------------------------------- statistics + HOM
def _counting(gen: Generator, n_f, f_rep_MHz):
    edges, L, U = gen.grid(n_f)
    dt = np.diff(edges)
    period = 1.0e3 / float(f_rep_MHz)
    dark = period - gen.W
    if dark <= 0:
        raise ValueError("repetition period must exceed the pulse window")
    P = np.eye(4, dtype=complex)
    for u in U:
        P = u @ P
    v0 = rs._fixed_point(expm(gen.L0 * dark) @ P, 2)
    V = np.zeros(12, dtype=complex)
    V[:4] = v0
    # FX2b: the same propagators started from |G><G| give the single-pulse excitation COUNT
    # n_exc = P_X(end of window) + photons emitted inside the window (linear in the start state).
    # FX2c: eta_prep is the PROBABILITY P(n >= 1) = 1 - p0, p0 = weight of the no-jump (no photon emitted
    # in the window) evolution that ends in |G> (the no-jump generator is L(t) - J), exact and in [0, 1].
    V0 = np.zeros(12, dtype=complex)
    V0[0] = 1.0
    vnj = np.zeros(4, dtype=complex)
    vnj[0] = 1.0
    px = np.empty(n_f + 1)
    px[0] = V[_IDX_PX].real
    J = gen.J
    for i in range(n_f):
        M = expm(rs._aug(L[i], J) * dt[i])
        V = M @ V
        V0 = M @ V0
        vnj = expm((L[i] - J) * dt[i]) @ vnj
        px[i + 1] = V[_IDX_PX].real
    rho_post = lb.unvec(V[:4].copy(), 2)
    n_exc = float(V0[_IDX_PX].real + (lb.trace_row(2) @ V0[4:8]).real)
    eta_prep = 1.0 - float(vnj[0].real)                      # vec index 0 = rho[G, G]
    V = expm(rs._aug(gen.L0, J) * dark) @ V
    tr = lb.trace_row(2)
    m1 = float((tr @ V[4:8]).real)
    m2 = float((tr @ V[8:]).real)
    return {"g2": m2 / m1 ** 2 if m1 > 1e-300 else NAN, "mean_photons": m1, "rho_post": rho_post, "px": px,
            "t": edges, "rho_start": lb.unvec(v0, 2), "eta_prep": eta_prep,
            "n_exc": n_exc}


def _group(U, m):
    """Compose m consecutive step propagators."""
    if m == 1:
        return U
    n = U.shape[0] // m
    out = np.empty((n, 4, 4), dtype=complex)
    for j in range(n):
        P = U[j * m]
        for k in range(1, m):
            P = U[j * m + k] @ P
        out[j] = P
    return out


def _hom_pulse(gen: Generator, Uq, h):
    """Pulse-region integrals on the quadrature grid of step propagators Uq (resonant_source._hom_pulse_once)."""
    dim = 2
    A = gen.A
    SA = lb.sprepost(A, np.eye(dim))
    a = lb.vec(A.conj())
    b = lb.vec((A.conj().T @ A).T)
    v = lb.vec(np.array([[1, 0], [0, 0]], dtype=complex))
    rhos = [v.copy()]
    for u in Uq:
        v = u @ v
        rhos.append(v.copy())
    R = np.array(rhos).T
    npts = R.shape[1]
    wq = np.full(npts, h)
    wq[0] = wq[-1] = 0.5 * h
    N_pulse = float(np.real(np.sum(wq * (b @ R))))
    active = SA @ R
    M = np.zeros((npts, npts))
    for j in range(npts):
        if j > 0:
            active[:, :j] = Uq[j - 1] @ active[:, :j]
        g = a @ active[:, :j + 1]
        M[:j + 1, j] = g.real ** 2 + g.imag ** 2
    F_aa = float(2.0 * (wq @ np.triu(M, 1) @ wq) + (wq ** 2) @ np.diag(M))
    X = (active * wq[None, :]) @ active.conj().T
    return F_aa, X, N_pulse, lb.unvec(R[:, -1], dim), a


def _hom(gen: Generator, nq0, fine_for, tol=HOM_TOL, max_doublings=HOM_DOUBLINGS):
    """Grange HOM overlap with Richardson-doubled quadrature. fine_for(nq) -> number of fine steps (a
    multiple of nq) used to build the quadrature-step propagators. -> (I, info)."""
    L = gen.L0
    P = lb._stationary_projector(L)
    Q = np.eye(L.shape[0]) - P
    Lp = L - P

    def evaluate(nq):
        n_f = fine_for(nq)
        _, _, U = gen.grid(n_f)
        Uq = _group(U, n_f // nq)
        F_aa, X, N_pulse, rho_end, a = _hom_pulse(gen, Uq, gen.W / nq)
        Y = rs._lyap_pair(Lp, Q @ X @ Q.conj().T)
        F_ab = float(np.real(a @ Y @ a.conj()))
        I_c, N_c = lb.indistinguishability(L, 2, rho_end, gen.A, return_counts=True)
        S_c = 0.5 * I_c * N_c ** 2
        Ntot = N_pulse + N_c
        if not Ntot > 1e-300:
            return NAN, Ntot                      # no photon is emitted (zero area): overlap undefined
        return (F_aa + 2.0 * F_ab + 2.0 * S_c) / Ntot ** 2, Ntot

    nq = nq0
    prev, Nprev = evaluate(nq)
    if not math.isfinite(prev):
        return NAN, {"N": float(Nprev), "n_quad": int(nq), "conv_err": NAN, "converged": False}
    out, Ntot, conv, err = prev, Nprev, False, NAN
    for _ in range(max_doublings):
        if 2 * nq > HOM_NODE_CAP:
            break
        nq *= 2
        cur, Ncur = evaluate(nq)
        err = abs(cur - prev) / 3.0
        out = cur + (cur - prev) / 3.0
        Ntot = Ncur + (Ncur - Nprev) / 3.0
        if err < tol * max(abs(out), 1e-300):
            conv = True
            break
        prev, Nprev = cur, Ncur
    return float(out), {"N": float(Ntot), "n_quad": int(nq), "conv_err": float(err), "converged": bool(conv)}


# ---------------------------------------------------------------- solver
class PolaronSolver(Solver):
    name = "polaron"
    label = "Polaron (tier P1, LA phonons)"
    summary = ("Driven two-level emitter with LA-phonon relaxation in the instantaneous dressed basis "
               "(Omega_r = sqrt(Z) Omega, Gamma_down/up = (pi/2) sin^2(2 theta) J(Lambda)[n+1, n]); one or "
               "several pulses (SUPER, dichromatic). Pulsed g2, photon number, inversion, HOM overlap.")
    handles_phonons = True

    def supports(self, plan, em, rt):
        if plan.family != "coherent" or plan.timing != "pulsed":
            return False, "polaron: coherent pulsed drives only"
        extra = sorted(set(plan.needs) - _ALLOWED_NEEDS)
        if extra:
            return False, f"polaron: unsupported needs {extra}"
        if not plan.pulses:
            return False, "polaron: plan has no pulses"
        for p in plan.pulses:
            if p.shape not in ("rect", "gauss"):
                return False, f"polaron: pulse shape {p.shape!r} not supported (rect|gauss)"
            if not (p.fwhm_ps > 0 and math.isfinite(p.area_rad) and p.area_rad >= 0):
                return False, "polaron: pulse needs fwhm_ps > 0 and a finite non-negative area"
        if plan.f_rep_MHz is None or plan.f_rep_MHz <= 0:
            return False, "polaron: pulsed plan needs f_rep_MHz"
        if "T_K" not in plan.inputs:
            return False, "polaron: plan.inputs['T_K'] missing (the engine supplies the temperature)"
        if plan.det_profile is not None and len(plan.det_profile) != 4:
            return False, "polaron: det_profile must be (det_off_meV, det_on_meV, t_gate_ps, t_r_ps)"
        W = PulseTrain(plan.pulses, plan.det_profile).W
        if 1e3 / float(plan.f_rep_MHz) <= W:
            return False, "polaron: repetition period must exceed the pulse window"
        return True, ""

    # --- grid choice
    @staticmethod
    def _grids(train: PulseTrain, B: float):
        rmax = max(train.rate_max(B), 1.0)
        dt = DT_FACTOR / rmax
        smin = train.sigma_min()
        if smin is not None:
            dt = min(dt, smin / rs.GAUSS_STEPS_PER_SIGMA)
        n_f = max(1, int(math.ceil(train.W / dt)))
        if smin is not None:
            nq0 = int(math.ceil(train.W / (smin / rs.GAUSS_STEPS_PER_SIGMA)))
        else:
            nq0 = 32
        nq0 = max(nq0, min(HOM_NODE_CAP // 2, int(math.ceil(train.W * rmax / 4.0))))
        nq0 = min(nq0, HOM_NODE_CAP // 2)
        return n_f, nq0, dt

    def _setup(self, plan, em, rt):
        T = float(plan.inputs["T_K"])
        pm = PhononModel(em.phonon, T)
        train = PulseTrain(plan.pulses, plan.det_profile)
        tl = rs.TwoLevel(float(rt.Gamma_ns), k_ns=float(rt.k_nr_ns), gamma_star_ns=float(rt.gamma_star_ns),
                         delta_ns=train.det0_ns)
        eid = float(plan.inputs.get("eid_ps", 0.0) or 0.0) * 1e-3
        return pm, train, Generator(train, tl, pm, eid_ns=eid)

    def run(self, plan, em, rt, with_hom=True):
        pm, train, gen = self._setup(plan, em, rt)
        _, nq0, dt = self._grids(train, pm.B)

        def fine_for(nq):
            """Fine steps for an nq-node quadrature: m substeps per node so that dt_sub <= dt."""
            m = max(1, int(math.ceil(train.W / nq / dt)))
            return nq * m

        n_f = fine_for(nq0)
        st = _counting(gen, n_f, plan.f_rep_MHz)
        reasons = ["polaron P1: weak-coupling dressed-state rates (no full polaron master equation)"]
        if len({round(c.dw, 9) for c in train.comps}) > 1 and pm.active:
            reasons.append("multi-colour drive with phonons: instantaneous-basis rates for a rotating beat "
                           "phase are an approximation [A]")
        if with_hom:
            I, info = _hom(gen, nq0, fine_for)
            if not info["converged"]:
                reasons.append(f"HOM quadrature not converged (rel err {info['conv_err']:.2e}, "
                               f"{info['n_quad']} nodes)")
        else:
            I, info = NAN, {}
        px = st["px"]
        ipk = int(np.argmax(px))
        _, lam_c, mix_c, _ = gen.dressed(np.array([0.5 * gen.W]))
        gd, gu = pm.rates(lam_c, mix_c)
        extra = {"eta_prep": float(st["eta_prep"]), "n_exc": float(st["n_exc"]), "P_X_end": float(st["rho_post"][1, 1].real),
                 "P_X_max": float(px[ipk]), "t_P_X_max_ps": float(st["t"][ipk] * 1e3),
                 "B_renorm": pm.B, "Z": pm.Z, "T_K": pm.T_K, "alpha_ps2": pm.alpha_ps2,
                 "n_fine_steps": int(n_f), "dt_fine_fs": float(train.W / n_f * 1e6),
                 "window_ps": float(train.W * 1e3), "hom": info, "det_X_meV": float(plan.pulses[0].detuning_meV),
                 "gamma_down_centre_ns": float(gd[0]), "gamma_up_centre_ns": float(gu[0]),
                 "notes": reasons}
        if plan.det_profile is not None:
            extra["lz"] = self._lz(train, pm.B)
        if gen.eid_ns > 0.0:
            extra["eid_ps"] = gen.eid_ns * 1e3
            reasons.append(f"excitation-induced dephasing gamma_EID = {gen.eid_ns * 1e3:.4g} ps x |Omega|^2 [E]")
        res = SolverResult(g2_0=float(st["g2"]), n_ph=float(st["mean_photons"]), I_single=float(I),
                           inversion=float(st["rho_post"][1, 1].real), extra=extra, valid=True,
                           reasons=tuple(reasons))
        pt = float(plan.inputs.get("p_trion", 0.0) or 0.0)
        return apply_trion(res, trion_probability(pt, math.sqrt(sum(float(p.area_rad) ** 2 for p in plan.pulses)))) if pt > 0.0 else res

    @staticmethod
    def _lz(train: PulseTrain, B: float) -> dict:
        """Landau-Zener check of a Stark gate: the rising-edge chirp through resonance is adiabatic when
        Omega^2 >> |d det/dt| (brief 2.1). Reported: Omega (rad/ns), the sweep rate of the RC edge
        |det_on - det_off|/t_r (rad/ns^2; inf for a sharp edge), the ratio Omega^2/sweep (adiabatic if >> 1) and
        the Landau-Zener diabatic probability exp(-pi Omega^2/(2 sweep)) for a linear sweep [DR]."""
        om = B * max(abs(c.amp) for c in train.comps)
        jump = abs(train._det_on - train._det_off)
        sweep = float("inf") if train._t_r <= 0.0 else jump / train._t_r
        if sweep == float("inf"):
            ratio, p_dia = 0.0, 1.0
        elif sweep == 0.0:
            ratio, p_dia = float("inf"), 0.0
        else:
            ratio, p_dia = om * om / sweep, math.exp(-math.pi * om * om / (2.0 * sweep))
        return {"omega_ns": float(om), "sweep_ns2": float(sweep), "adiabaticity": float(ratio),
                "adiabatic": bool(ratio > 10.0), "p_lz_diabatic": float(p_dia)}

    def eta_seconds(self, plan, n_points):
        try:
            train = PulseTrain(plan.pulses, plan.det_profile)
            n_f, _, _ = self._grids(train, 0.95)
        except Exception:                                          # noqa: BLE001
            n_f = 3000
        return (n_f * 4e-4 * 3.0 + 1.0) * max(1, int(n_points))


register(PolaronSolver())
