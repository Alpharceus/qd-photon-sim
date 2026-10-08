"""Resonant / two-photon-excited Fock-state source tier (numpy/scipy only).

Specification: pieces F1-F4 of the quantum-tier spec. Built on the Liouvillian
core ``fsim_core.lindblad`` (Markovian baths). This module must not import
``fsim_core.device`` and no legacy module imports it, so verify_fsim (51/51)
and audit_physics (23/23) are untouched by construction. Regression suite:
verify/verify_resonant_source.py.

Units as in lindblad.py: rates and angular frequencies in 1/ns, times in ns,
energies in meV. Provenance tags: [V] verified against the cited paper, [DR]
derived, [E] estimate or class range, [A] assumption.

Why a separate tier
-------------------
Fock-state boson sampling needs per-photon g2(0) -> 0, pairwise HOM
indistinguishability I -> 1 and high end-to-end efficiency (Renema et al., PRL
120, 220502 (2018); Oszmaniec & Brod, NJP 20, 092002 (2018)). The QD
demonstration with 20 photons (Wang et al., PRL 123, 250503 (2019)) used
resonantly driven QD-micropillars (Somaschi et al., Nat. Photonics 10, 340
(2016); Ding et al., PRL 116, 020401 (2016)). Electrically injected,
incoherently pumped sources do not reach this regime at cryogenic temperature;
the V-a Gamma_0 = 119 ueV is a time-integrated non-resonant linewidth mixing
fast dephasing with spectral diffusion and is NOT used as a homogeneous width
here. Transform-limited
lines under resonant excitation: Kuhlmann et al., Nat. Commun. 6, 8204 (2015).

F1  Resonant pulsed two-level source
------------------------------------
H(t) = Delta |e><e| + Omega(t)/2 (|e><g| + h.c.); dissipators: radiative decay
Gamma (D[sqrt(Gamma)|g><e|], the photon channel), photonless loss k
(D[sqrt(k)|g><e|]), pure dephasing gamma* (D[sqrt(gamma*)|e><e|], the
lindblad.py convention: e-g coherence decays at (Gamma + k + gamma*)/2). A
charged exciton X+ is treated as its own two-level system (ground state = the
resident hole) [A]; the two-level model has no trion fine structure.

Re-excitation. Per-period factorial moments use the counting hierarchy
    dp/dt = L(t) p,  dm1/dt = L(t) m1 + J p,  dm2/dt = L(t) m2 + 2 J m1,
    g2(0) = Tr m2 / (Tr m1)^2.
[DR] For a rectangular pi pulse at small Gamma*tau_p, emission at time t
inside the pulse has probability density Gamma sin^2(Omega t/2), after which
the remaining area pi - Omega t re-excites with probability cos^2(Omega t/2).
Hence P2 ~ Gamma*tau_p <sin^2 u cos^2 u>_{u in [0, pi/2]} = Gamma*tau_p/8 and
    g2(0) ~ Gamma*tau_p / 4.                                           [DR]
g2(0) depends on the pulse only through Gamma*tau_p. Purcell enhancement
(Gamma = Gamma_0 (1 + F_P L)) therefore RAISES re-excitation at fixed pulse
length: a cavity that speeds the decay by F_P must shorten the pulse by the
same factor to keep g2(0). Re-excitation-induced two-photon emission: Fischer
et al., Nat. Phys. 13, 649 (2017). Purcell convention: ``purcell_rate`` below.

HOM integral for the photon emitted during and after one pulse (Grange et al.,
PRL 114, 193601 (2015) Eq. 2; Bylander et al., EPJD 22, 295 (2003)),
    I = [ int int_{[0,inf)^2} |G1(t,t')|^2 dt dt' ] / (int n(t) dt)^2,
    G1(t,t') = <A^+(t') A(t)> = Tr[A^+ U(t',t)(A rho(t))]  (t <= t'),
    n(t) = Tr[A^+ A rho(t)],  A = sqrt(Gamma)|g><e|,
with the emitter in |g> at the pulse start and U the propagator of the
time-dependent L(t) for both t and tau = t' - t. |G1|^2 is symmetric, so the
full-plane integral is split into  (a) both times inside the pulse (2-D
trapezoid on the step grid), (b) t inside, t' after the pulse, and (c) both
after the pulse. (b) and (c) are exact: they are Lyapunov/Sylvester solves
with the stationary mode removed (derivation in
``lindblad.indistinguishability``; no diagonalisability needed). For (b) the
t-integral over the pulse is a trapezoid sum of w_i w_i^+, with w_i the
vector A rho(t_i) propagated to the pulse end. (c) is
``lindblad.indistinguishability`` started from the post-pulse state. The
trapezoid grid is doubled until the Richardson-corrected value converges
(h^2 error: I_R = I_2n + (I_2n - I_n)/3). This is the FIRST-ORDER-coherence
definition; for g2(0) != 0 it ignores the two-photon component, i.e. it is the
HOM visibility of the single-photon wavepacket, which is the work-order
definition.

F2  Two-photon excitation of the cascade  (see ``tpe_*`` functions)
F3  Spectral diffusion and HOM versus photon separation (``hom_*``)
F4  Boson-sampling figure of merit, card loader (``bs_figure_of_merit``)
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
from scipy.linalg import expm, solve, solve_continuous_lyapunov

from fsim_core import lindblad as lb
from fsim_core.spectral import cavity_transmission

#: Gaussian pulses are truncated at +-8 sigma around the centre; the neglected
#: area fraction is erfc(8/sqrt 2) = 1.2e-15. [A] numerical choice.
GAUSS_WINDOW_SIGMA = 8.0
#: Default Gaussian step: dt <= sigma / 8 (midpoint rule on a Gaussian is
#: spectrally accurate; error ~ exp(-2 pi^2 (sigma/dt)^2)). [A] numerical choice.
GAUSS_STEPS_PER_SIGMA = 8


# ------------------------------------------------------------------ pulses

def gaussian_sigma(duration):
    """Amplitude std sigma of a Gaussian pulse whose INTENSITY (Omega^2)
    FWHM is ``duration``: Omega^2 ~ exp(-(t-tc)^2/sigma^2), FWHM =
    2 sigma sqrt(ln 2). [DR]"""
    return float(duration) / (2.0 * np.sqrt(np.log(2.0)))


def pulse_window(shape, duration):
    """Length (ns) of the time window the pulse occupies: ``duration`` for a
    rectangle, 2 * GAUSS_WINDOW_SIGMA * sigma for a Gaussian (centred)."""
    if shape == "rect":
        return float(duration)
    if shape == "gauss":
        return 2.0 * GAUSS_WINDOW_SIGMA * gaussian_sigma(duration)
    raise ValueError("shape must be 'rect' or 'gauss'")


def pulse_envelope(t, shape, area, duration):
    """Rabi frequency Omega(t) (rad/ns) of a pulse of area theta =
    int Omega dt (rad), starting at t = 0.

    shape 'rect' : Omega = area/duration on [0, duration), 0 elsewhere.
    shape 'gauss': Omega = area/(sigma sqrt(2 pi)) exp(-(t-tc)^2/(2 sigma^2))
                   on [0, 2 GAUSS_WINDOW_SIGMA sigma], tc = window/2, with
                   ``duration`` the INTENSITY FWHM (gaussian_sigma)."""
    t = np.asarray(t, dtype=float)
    if shape == "rect":
        return np.where((t >= 0.0) & (t < duration), area / duration, 0.0)
    if shape == "gauss":
        sig = gaussian_sigma(duration)
        win = 2.0 * GAUSS_WINDOW_SIGMA * sig
        tc = 0.5 * win
        env = area / (sig * np.sqrt(2.0 * np.pi)) * np.exp(-0.5 * ((t - tc) / sig) ** 2)
        return np.where((t >= 0.0) & (t <= win), env, 0.0)
    raise ValueError("shape must be 'rect' or 'gauss'")


def _default_steps(shape, duration, n_steps, max_dt=None):
    if n_steps is not None:
        return int(n_steps)
    win = pulse_window(shape, duration)
    if shape == "rect":
        n = 1
    else:
        n = int(np.ceil(win / (gaussian_sigma(duration) / GAUSS_STEPS_PER_SIGMA)))
    if max_dt is not None:
        n = max(n, int(np.ceil(win / max_dt)))
    return n


def purcell_rate(gamma0_ns, F_P, detuning, line_fwhm, kappa):
    """Purcell-enhanced decay rate Gamma = Gamma_0 (1 + F_P L), [DR] with the
    overlap convention of ``spectral.cavity_transmission``: L is the overlap of
    the emitter line (unit-area Lorentzian, FWHM ``line_fwhm``, offset
    ``detuning`` from the cavity) with the peak-normalised cavity Lorentzian
    (FWHM ``kappa``), L = kappa (line_fwhm + kappa)/4 / (detuning^2 +
    ((line_fwhm + kappa)/2)^2), so L(0) = kappa/(line_fwhm + kappa) -> 1 for a
    narrow line and F_P is the usual resonant Purcell factor. detuning,
    line_fwhm and kappa share one unit (any); the result has the unit of
    gamma0_ns."""
    return float(gamma0_ns) * (1.0 + float(F_P) * float(cavity_transmission(detuning, line_fwhm, kappa)))


# ------------------------------------------------------------------ model

@dataclass(frozen=True)
class TwoLevel:
    """Two-level emitter (ground |g>, excited |e>); rates in 1/ns.

    gamma_ns      radiative decay Gamma (the photon channel)
    k_ns          photonless loss k
    gamma_star_ns pure dephasing gamma* (lindblad.py convention)
    delta_ns      transition detuning from the laser, H = delta |e><e|
    label         'X' or 'X+' (X+ = own two-level system, ground = resident hole)
    """
    gamma_ns: float
    k_ns: float = 0.0
    gamma_star_ns: float = 0.0
    delta_ns: float = 0.0
    label: str = "X"

    def system(self, omega_ns=0.0):
        deph = [(self.gamma_star_ns, (0.0, 1.0))] if self.gamma_star_ns > 0 else []
        return lb.build_system(levels=2, gamma_X_ns=self.gamma_ns, k_X=self.k_ns,
                               deph=deph, det_X_ns=self.delta_ns, omega_ns=omega_ns)


@dataclass
class _Pieces:
    """Liouvillian pieces of a driven emitter: L(t) = L0 + Omega(t) L1."""
    L0: np.ndarray
    L1: np.ndarray
    J: np.ndarray
    A: np.ndarray
    ops: dict
    dim: int
    max_dt: float | None = None     # step cap (ns) for a large static detuning (TPE)


def _pieces(tl):
    s0 = tl.system(0.0)
    sd = tl.system(1.0)
    L1 = lb.hamiltonian_super(0.5 * sd.ops["sigma_x_drive"])
    A = np.sqrt(tl.gamma_ns) * s0.ops["c_X"]
    return _Pieces(L0=s0.L, L1=L1, J=lb.counting_jump(s0), A=A, ops=s0.ops, dim=s0.dim)


def _grid(shape, area, duration, n_steps, max_dt=None):
    """Step edges and the Rabi frequency at each step midpoint."""
    n = _default_steps(shape, duration, n_steps, max_dt)
    win = pulse_window(shape, duration)
    edges = np.linspace(0.0, win, n + 1)
    mid = 0.5 * (edges[:-1] + edges[1:])
    return edges, pulse_envelope(mid, shape, area, duration)


def _step_matrices(pc, om, dt, shape):
    """exp((L0 + om_k L1) dt) per step; a rectangle needs one expm."""
    if shape == "rect":
        U = expm((pc.L0 + om[0] * pc.L1) * dt[0])
        return [U] * len(om)
    return [expm((pc.L0 + o * pc.L1) * d) for o, d in zip(om, dt)]


def _ground(dim):
    r = np.zeros((dim, dim), dtype=complex)
    r[0, 0] = 1.0
    return r


def propagate_pulse(tl, shape, area, duration, rho0=None, n_steps=None):
    """Density matrix after the pulse (piecewise-constant midpoint propagation
    of L(t), as lindblad.evolve_td). Returns (rho_end, list of rho at the step
    edges)."""
    pc = _pieces(tl)
    edges, om = _grid(shape, area, duration, n_steps, pc.max_dt)
    dt = np.diff(edges)
    Us = _step_matrices(pc, om, dt, shape)
    v = lb.vec(_ground(pc.dim) if rho0 is None else rho0)
    traj = [lb.unvec(v.copy(), pc.dim)]
    for U in Us:
        v = U @ v
        traj.append(lb.unvec(v.copy(), pc.dim))
    return traj[-1], traj


def _aug(L, J):
    """Generator of the counting hierarchy (p, m1, m2)."""
    N = L.shape[0]
    Z = np.zeros((N, N), dtype=complex)
    return np.block([[L, Z, Z], [J, L, Z], [Z, 2.0 * J, L]])


def _fixed_point(Phi, dim):
    """Periodic state: Phi v = v with Tr = 1 (linear solve; one row of
    Phi - 1 replaced by the trace row)."""
    N = dim * dim
    M = Phi - np.eye(N)
    M[0, :] = lb.trace_row(dim)
    b = np.zeros(N, dtype=complex)
    b[0] = 1.0
    return solve(M, b)


def _stats(pc, shape, area, duration, f_rep_MHz, n_steps, periodic):
    """Counting hierarchy over one repetition period for any pieces ``pc``."""
    period = 1.0e3 / float(f_rep_MHz)
    win = pulse_window(shape, duration)
    if period <= win:
        raise ValueError("repetition period must exceed the pulse window")
    edges, om = _grid(shape, area, duration, n_steps, pc.max_dt)
    dt = np.diff(edges)
    N = pc.dim ** 2
    dark = period - win
    if periodic:
        Us = _step_matrices(pc, om, dt, shape)
        P = np.eye(N, dtype=complex)
        for U in Us:
            P = U @ P
        v_start = _fixed_point(expm(pc.L0 * dark) @ P, pc.dim)
    else:
        v_start = lb.vec(_ground(pc.dim))
    V = np.zeros(3 * N, dtype=complex)
    V[:N] = v_start
    cache = {}
    for o, d in zip(om, dt):
        key = (round(float(o), 14), round(float(d), 14))
        if key not in cache:
            cache[key] = expm(_aug(pc.L0 + o * pc.L1, pc.J) * d)
        V = cache[key] @ V
    rho_post = lb.unvec(V[:N].copy(), pc.dim)
    V = expm(_aug(pc.L0, pc.J) * dark) @ V
    tr = lb.trace_row(pc.dim)
    m1 = float((tr @ V[N:2 * N]).real)
    m2 = float((tr @ V[2 * N:]).real)
    return {"g2": m2 / m1 ** 2 if m1 > 1e-300 else float("nan"),
            "mean_photons": m1, "rho_post_pulse": rho_post,
            "rho_start": lb.unvec(v_start, pc.dim)}


def pulsed_photon_statistics(tl, shape, area, duration, f_rep_MHz=80.0,
                             n_steps=None, periodic=True):
    """Per-period factorial moments of the photon stream of a pulsed two-level
    emitter under a time-dependent L(t) (counting hierarchy in the module
    docstring; piecewise-constant midpoint propagation). Counting runs over
    the whole repetition period T = 1e3/f_rep_MHz ns.

    periodic=True starts each period in the periodic fixed point of the one-
    period map (what lindblad.pulsed_counting does); False starts in |g>.
    Returns dict(g2, mean_photons, rho_post_pulse (state at the pulse end,
    counting included), rho_start, P_e_post)."""
    out = _stats(_pieces(tl), shape, area, duration, f_rep_MHz, n_steps, periodic)
    out["P_e_post"] = float(out["rho_post_pulse"][1, 1].real)
    return out


# ------------------------------------------------------------------ HOM

def _lyap_pair(Lp, X):
    """Solve Lp Y + Y Lp^+ = -X."""
    return solve_continuous_lyapunov(Lp, -X)


def _hom_pulse_once(pc, shape, area, duration, n_quad, rho0=None):
    """Pulse-region full-plane integrals on an n_quad-step trapezoid grid.
    Returns (F_aa + 2 F_ab, N_pulse, rho_end)."""
    dim = pc.dim
    n = dim * dim
    edges, om = _grid(shape, area, duration, n_quad)
    dt = np.diff(edges)
    h = float(dt[0])
    Us = _step_matrices(pc, om, dt, shape)
    SA = lb.sprepost(pc.A, np.eye(dim))
    a = lb.vec(pc.A.conj())               # (A^+)^T: Tr[A^+ X] = a . vec(X)
    b = lb.vec((pc.A.conj().T @ pc.A).T)
    # rho_i at the grid points
    v = lb.vec(_ground(dim) if rho0 is None else rho0)
    rhos = [v.copy()]
    for U in Us:
        v = U @ v
        rhos.append(v.copy())
    R = np.array(rhos).T                              # n x (n_quad+1)
    npts = R.shape[1]
    wq = np.full(npts, h)
    wq[0] = wq[-1] = 0.5 * h
    N_pulse = float(np.real(np.sum(wq * (b @ R))))
    # u_i = A rho_i propagated; G[i, j] = a . u_i(t_j), j >= i
    Ucols = SA @ R                                    # columns u_i
    M = np.zeros((npts, npts))                        # |G1(t_i, t_j)|^2, j >= i (upper triangle)
    active = Ucols.copy()
    for j in range(npts):
        if j > 0:
            active[:, :j] = Us[j - 1] @ active[:, :j]   # advance columns i < j to t_j
        g = a @ active[:, :j + 1]
        M[:j + 1, j] = g.real ** 2 + g.imag ** 2
    # full-square trapezoid of the symmetric |G1|^2 = twice the strict upper triangle + the diagonal
    F_aa = float(2.0 * (wq @ np.triu(M, 1) @ wq) + (wq ** 2) @ np.diag(M))
    # (b): w_i = u_i at the pulse end; X = sum_i wq_i w_i w_i^+
    W = active                                        # all columns at t = T_p
    X = (W * wq[None, :]) @ W.conj().T
    return F_aa, X, N_pulse, lb.unvec(R[:, -1], dim), (a, n)


def pulsed_indistinguishability(tl, shape, area, duration, n_quad=None, tol=1e-8,
                                max_doublings=7, rho0=None, return_info=False):
    """HOM indistinguishability of the photon emitted during and after one
    pulse (definition and method in the module docstring). The emitter starts
    in |g> at the pulse start [A] unless ``rho0`` is given
    (rho0 = |e><e| with area 0 is the exactly solvable free-decay test of the
    three-region assembly: I = Gamma_t/(Gamma_t + gamma*)).

    n_quad : trapezoid steps in the pulse region (default 32 for a
             rectangle, the propagation grid for a Gaussian). The grid is
             doubled until the Richardson-extrapolated value changes by < tol
             (rel) or ``max_doublings`` is reached. return_info=True returns
             (I, info) with info = dict(N, n_quad, conv_err, converged)."""
    pc = _pieces(tl)
    dim = pc.dim
    n0 = int(n_quad) if n_quad is not None else (
        32 if shape == "rect" else _default_steps(shape, duration, None))
    L = np.asarray(pc.L0, dtype=complex)
    P = lb._stationary_projector(L)
    Q = np.eye(L.shape[0]) - P
    Lp = L - P

    def evaluate(nq):
        F_aa, X, N_pulse, rho_end, (a, n) = _hom_pulse_once(pc, shape, area, duration, nq, rho0)
        Xq = Q @ X @ Q.conj().T
        Y = _lyap_pair(Lp, Xq)
        F_ab = float(np.real(a @ Y @ a.conj()))
        I_c, N_c = lb.indistinguishability(L, dim, rho_end, pc.A, return_counts=True)
        S_c = 0.5 * I_c * N_c ** 2
        Ntot = N_pulse + N_c
        return (F_aa + 2.0 * F_ab + 2.0 * S_c) / Ntot ** 2, Ntot

    nq = n0
    prev, Nprev = evaluate(nq)
    out, Ntot, conv, err = prev, Nprev, False, float("nan")
    for _ in range(max_doublings):
        nq *= 2
        cur, Ncur = evaluate(nq)
        err = abs(cur - prev) / 3.0          # estimated error of cur (h^2 rule)
        out = cur + (cur - prev) / 3.0       # Richardson extrapolation
        Ntot = Ncur + (Ncur - Nprev) / 3.0
        if err < tol * max(abs(out), 1e-300):
            conv = True
            break
        prev, Nprev = cur, Ncur
    if return_info:
        return float(out), {"N": float(Ntot), "n_quad": int(nq), "conv_err": float(err),
                            "converged": bool(conv)}
    return float(out)


# ------------------------------------------------------------------ Monte Carlo (second method)

def pulsed_counting_mc(tl, area, duration, f_rep_MHz=80.0, n_traj=2_000_000,
                       seed=20261008, n_batches=50):
    """Gillespie / quantum-jump waiting-time Monte Carlo of the rectangular-
    pulse two-level source, started in |g> (compare ``pulsed_photon_statistics
    (..., periodic=False)``). An independent method: it never builds a
    Liouvillian.

    Between jumps the unnormalised state evolves under
    H_eff = [[0, Om/2], [Om/2, Delta - i (Gamma + k)/2]] (basis g, e); the
    survival S(t) = |psi(t)|^2 is the probability of no jump, so a jump
    occurs when S hits a uniform r. After a jump the emitter is in |g>
    (counted with probability Gamma/(Gamma + k), else photonless loss).
    Inside the pulse S(t) is found by bisection of the exact eigen-
    decomposition; after the pulse the laser is off and S(t) = |c_g|^2 +
    |c_e|^2 exp(-(Gamma + k) t), inverted in closed form -- no time-step
    error. Photons are counted for 0 <= t < T = 1e3/f_rep_MHz. Pure dephasing
    is not supported (raises). Returns dict(g2, g2_se, mean, mean_se, n_traj);
    the standard error is the spread of ``n_batches`` batch estimates.
    Fixed ``seed`` -> bit-reproducible."""
    if tl.gamma_star_ns != 0.0:
        raise ValueError("pulsed_counting_mc supports gamma* = 0 only")
    G, k, D = tl.gamma_ns, tl.k_ns, tl.delta_ns
    Gt = G + k
    om = float(area) / float(duration)
    tp = float(duration)
    period = 1.0e3 / float(f_rep_MHz)
    H = np.array([[0.0, 0.5 * om], [0.5 * om, D - 0.5j * Gt]], dtype=complex)
    lam, V = np.linalg.eig(H)
    c0 = np.linalg.solve(V, np.array([1.0, 0.0], dtype=complex))

    def psi(t):
        ph = np.exp(-1j * np.outer(t, lam)) * c0[None, :]
        return ph @ V.T                                  # (len t, 2): (c_g, c_e)

    rng = np.random.default_rng(seed)
    n = np.zeros(n_traj, dtype=np.int64)
    s = np.zeros(n_traj)
    act = np.arange(n_traj)
    for _ in range(200):
        if act.size == 0:
            break
        r = rng.random(act.size)
        u = rng.random(act.size)                         # photon vs loss
        rem = tp - s[act]
        S_end = np.sum(np.abs(psi(rem)) ** 2, axis=1)
        inside = r >= S_end
        # ---- jump inside the pulse: bisection of S(t') = r on [0, rem]
        ia = np.flatnonzero(inside)
        if ia.size:
            lo = np.zeros(ia.size)
            hi = rem[ia].copy()
            rr = r[ia]
            for _b in range(60):
                mid = 0.5 * (lo + hi)
                Sm = np.sum(np.abs(psi(mid)) ** 2, axis=1)
                go_right = Sm > rr                       # still surviving
                lo = np.where(go_right, mid, lo)
                hi = np.where(go_right, hi, mid)
            tj = 0.5 * (lo + hi)
            idx = act[ia]
            s[idx] += tj
            n[idx] += (u[ia] < G / Gt).astype(np.int64)
        # ---- no jump inside the pulse: free decay afterwards
        oa = np.flatnonzero(~inside)
        if oa.size:
            p = psi(rem[oa])
            cg2 = np.abs(p[:, 0]) ** 2
            ce2 = np.abs(p[:, 1]) ** 2
            rr = r[oa]
            has = rr > cg2
            with np.errstate(divide="ignore", invalid="ignore"):
                t_free = -np.log((rr - cg2) / ce2) / Gt
            t_abs = tp + t_free
            counted = has & (t_abs < period) & (u[oa] < G / Gt)
            n[act[oa]] += counted.astype(np.int64)
        act = act[ia]
    nf = n.astype(float)
    nb = n_traj // n_batches
    m1 = np.array([nf[i * nb:(i + 1) * nb].mean() for i in range(n_batches)])
    m2 = np.array([(nf[i * nb:(i + 1) * nb] * (nf[i * nb:(i + 1) * nb] - 1)).mean()
                   for i in range(n_batches)])
    g2b = m2 / m1 ** 2
    return {"g2": float(m2.mean() / m1.mean() ** 2),
            "g2_se": float(g2b.std(ddof=1) / np.sqrt(n_batches)),
            "mean": float(nf.mean()), "mean_se": float(m1.std(ddof=1) / np.sqrt(n_batches)),
            "n_traj": int(n_traj)}


# ====================================================================== F2
# Two-photon excitation (TPE) of the G-X-XX cascade
#
# Laser frame at the two-photon resonance omega_L = (E_XX - E_G)/(2 hbar):
#     H = delta |X><X| + (Omega(t)/2) (|X><G| + |XX><X| + h.c.),
# the XX level at energy 0 (two-photon resonant), G at 0, X at
# delta = -Delta_XX/2 with Delta_XX = (E_XX - E_X) - (E_X - E_G) the card /
# lindblad.build_system ``det_XX_ns`` convention (Delta_XX < 0 for a bound
# biexciton, so delta > 0). Note: an alternative convention defines
# delta = Delta_XX/(2 hbar); the opposite sign is immaterial for populations
# (H -> -H* with the gauge
# |X> -> -|X> maps one into the other; checked in verify (t1)), and the
# physical X level of a bound biexciton lies ABOVE the two-photon energy.
#
# Exact closed-system solution [DR] (used by the verifier as an independent
# closed form). With bright |B> = (|G> + |XX>)/sqrt 2 and dark |D> =
# (|G> - |XX>)/sqrt 2: H|D> = 0, <X|H|B> = Omega/sqrt 2, so {B, X} is a
# two-level system with eigenvalues lambda_pm = (delta +- R)/2,
# R = sqrt(delta^2 + 2 Omega^2). Starting in |G> = (B + D)/sqrt 2,
#     c_XX(t) = (b(t) - 1)/2,  b(t) = <B|psi(t)> .
# The low branch carries weight ~ 1 - (Omega/delta)^2/2 and phase
# exp(-i lambda_- t) with lambda_- = (delta - R)/2, so
#     P_XX(t) ~ sin^2(Omega_eff t / 2),
#     Omega_eff = (R - delta)/2 = Omega^2/(2 delta) [1 - Omega^2/(2 delta^2) + ...].
# The leading term Omega^2/(2 delta) is the adiabatic-elimination result (the
# AC Stark shifts of G and XX are equal, -Omega^2/(4 delta), and cancel); the
# first correction is RELATIVE -(Omega/delta)^2/2, i.e. -1.25e-3 at
# Omega/delta = 0.05. For a smooth (Gaussian) pulse the B -> low-branch
# following is adiabatic and the final state is (e^{-i phi} |B> + |D>)/sqrt 2,
# phi = int |lambda_-(Omega(t))| dt, so full G -> XX transfer needs phi = pi.

def tpe_detuning_ns(delta_xx_meV):
    """X-level detuning delta (rad/ns) in the laser frame at the two-photon
    resonance, delta = -Delta_XX/(2 hbar) with the card convention (module
    comment above): Delta_XX = -2 meV (binding energy 2 meV) -> delta > 0."""
    return -float(lb.mev_to_rate(delta_xx_meV)) / 2.0


def tpe_exact_rabi(omega_ns, delta_ns):
    """Exact low-branch frequency (sqrt(delta^2 + 2 Omega^2) - |delta|)/2 of
    the closed symmetric ladder (derivation above); -> Omega^2/(2|delta|)."""
    d = abs(float(delta_ns))
    return 0.5 * (np.sqrt(d * d + 2.0 * float(omega_ns) ** 2) - d)


@dataclass(frozen=True)
class Cascade:
    """Ladder |G>, |X>, |XX> in the two-photon laser frame (rates in 1/ns).

    gamma_X_ns, gamma_XX_ns radiative decays X -> G, XX -> X; k_X_ns, k_XX_ns
    photonless losses; deph lindblad.build_system pure-dephasing list over
    (G, X, XX) (e.g. ((gamma*, (0, 1, 0)),) = X-G coherence only);
    delta_ns the X-level detuning in the laser frame (tpe_detuning_ns)."""
    gamma_X_ns: float
    gamma_XX_ns: float = 0.0
    k_X_ns: float = 0.0
    k_XX_ns: float = 0.0
    deph: tuple = ()
    delta_ns: float = 0.0

    def system(self):
        return lb.build_system(levels=3, gamma_X_ns=self.gamma_X_ns, gamma_XX_ns=self.gamma_XX_ns,
                               k_X=self.k_X_ns, k_XX=self.k_XX_ns, deph=list(self.deph),
                               det_X_ns=self.delta_ns, det_XX_ns=-self.delta_ns)

    def closed(self):
        """Same Hamiltonian, no dissipation (the Gamma tau_p -> 0 limit)."""
        return Cascade(0.0, 0.0, 0.0, 0.0, (), self.delta_ns)


def _tpe_drive_op(ops):
    return ops["c_X"] + ops["c_X"].conj().T + ops["c_XX"] + ops["c_XX"].conj().T


def _pieces_cascade(casc, channels=("X",), max_dt_per_delta=0.5):
    s0 = casc.system()
    L1 = lb.hamiltonian_super(0.5 * _tpe_drive_op(s0.ops))
    J = lb.counting_jump(s0, names=list(channels))
    A = np.sqrt(casc.gamma_X_ns) * s0.ops["c_X"]
    max_dt = max_dt_per_delta / abs(casc.delta_ns) if casc.delta_ns else None
    return _Pieces(L0=s0.L, L1=L1, J=J, A=A, ops=s0.ops, dim=s0.dim, max_dt=max_dt)


def tpe_population_series(casc, omega_ns, times):
    """(P_G, P_X, P_XX) at ``times`` for a CONSTANT drive Omega, starting in
    |G> (lindblad.evolve on the cascade Liouvillian)."""
    pc = _pieces_cascade(casc)
    rhos = lb.evolve(pc.L0 + float(omega_ns) * pc.L1, _ground(3), times)
    return np.array([np.real(np.diag(r)) for r in rhos])


def _tpe_ode(delta_ns, shape, area, duration, rtol=1e-13, atol=1e-15):
    """Closed-system amplitudes (c_G, c_X, c_XX) after the pulse, by a
    high-order Runge-Kutta integration of the Schroedinger equation -- an
    integrator independent of the expm propagation used elsewhere."""
    from scipy.integrate import solve_ivp
    win = pulse_window(shape, duration)

    def rhs(t, y):
        om = float(pulse_envelope(t, shape, area, duration))
        cg, cx, cxx = y
        return [-0.5j * om * cx, -1j * (delta_ns * cx + 0.5 * om * (cg + cxx)), -0.5j * om * cx]

    sol = solve_ivp(rhs, (0.0, win), np.array([1, 0, 0], dtype=complex), method="DOP853",
                    rtol=rtol, atol=atol, max_step=win / 200.0)
    return sol.y[:, -1]


def _quad(f, a, b, points=None):
    from scipy.integrate import quad
    return quad(f, a, b, points=points, epsabs=0, epsrel=1e-11, limit=400)[0]


@lru_cache(maxsize=256)
def _tpe_ode_cached(delta_ns, shape, area, duration):
    return tuple(_tpe_ode(delta_ns, shape, area, duration))


@lru_cache(maxsize=256)
def _tpe_area_cached(delta_ns, shape, duration):
    return tpe_pi_area(delta_ns, shape, duration)


def tpe_prepare(casc, shape, area, duration):
    """State after a closed-system (dissipation off) TPE pulse: the
    Gamma tau_p -> 0 limit. Returns dict(psi, rho, P=(P_G, P_X, P_XX))."""
    psi = np.array(_tpe_ode_cached(float(casc.delta_ns), shape, float(area), float(duration)))
    return {"psi": psi, "rho": np.outer(psi, psi.conj()),
            "P": tuple(float(abs(c) ** 2) for c in psi)}


def tpe_pi_area(delta_ns, shape, duration):
    """Pulse area theta = int Omega dt (rad) of the single-photon Rabi
    frequency that transfers G -> XX with unit fidelity, closed system.

    gauss: the adiabatic phase estimate int |lambda_-(Omega(t))| dt = pi
           (brentq) is refined by solving Im(c_G conj(c_XX)) = 0 on the exact
           amplitudes: c_G/c_XX = i cot(phi/2) in the bright/dark picture,
           which crosses zero LINEARLY at phi = pi, so the residual G
           amplitude is driven to ~1e-14 (maximising P_XX would only reach
           sqrt(machine eps)).
    rect : sudden switching excites the upper branch, so exact unit transfer
           is not guaranteed; the area maximising P_XX near the adiabatic
           estimate is returned (bounded Brent)."""
    from scipy.optimize import brentq, minimize_scalar
    delta = abs(float(delta_ns))
    win = pulse_window(shape, duration)

    def phase(theta):
        f = lambda t: tpe_exact_rabi(float(pulse_envelope(t, shape, theta, duration)), delta)
        return _quad(f, 0.0, win, None if shape == "rect" else [0.5 * win])

    th0 = brentq(lambda th: phase(th) - np.pi, 1e-9, 1e3 * max(1.0, delta * duration),
                 xtol=1e-14, rtol=1e-14)
    if shape == "gauss":
        def f(th):
            psi = _tpe_ode(delta, shape, th, duration)
            return float(np.imag(psi[0] * np.conj(psi[2])) / abs(psi[2]) ** 2)
        return float(brentq(f, 0.95 * th0, 1.05 * th0, xtol=1e-15, rtol=8.9e-16))
    res = minimize_scalar(lambda th: -abs(_tpe_ode(delta, shape, th, duration)[2]) ** 2,
                          bounds=(0.9 * th0, 1.1 * th0), method="bounded",
                          options={"xatol": 1e-12})
    return float(res.x)


def tpe_photon_statistics(casc, shape, area, duration, f_rep_MHz=80.0, n_steps=None,
                          channel="X", periodic=True):
    """Per-period factorial moments of the ``channel`` photon stream ('X' =
    the X-line photons, the line a spectral filter selects; 'all' = X + XX)
    under TPE with the dissipators ON during the pulse (finite Gamma tau_p).
    Same hierarchy and midpoint propagation as pulsed_photon_statistics; the
    step is capped at 0.5/|delta| so the static detuning is resolved
    (verify (t4) checks convergence in n_steps). Returns dict(g2,
    mean_photons, rho_post_pulse, rho_start, P_XX_post)."""
    names = ("X", "XX") if channel == "all" else (channel,)
    pc = _pieces_cascade(casc, names)
    out = _stats(pc, shape, area, duration, f_rep_MHz, n_steps, periodic)
    out["P_XX_post"] = float(out["rho_post_pulse"][2, 2].real)
    return out


def tpe_cascade_hom(casc, shape, area, duration):
    """HOM indistinguishability of the X-line photon emitted after a TPE-
    prepared biexciton, in the Gamma tau_p -> 0 limit: the pulse is applied
    with dissipation OFF (tpe_prepare), then lindblad.indistinguishability
    (Lyapunov solve, valid for defective L) runs on the dissipative cascade
    with A = sqrt(gamma_X)|G><X|.

    For an exact |XX> start the cascade identity holds [DR]: the X photon of
    a cascade is delayed by the stochastic XX lifetime, which imprints a
    random start time on the X wavepacket (Schoell et al., PRL 125, 233605
    (2020)); with pure dephasing on the X-G coherence only,
        I_{XX->X} = I_X Gamma_XX/(Gamma_XX + Gamma_X),
        Gamma_X = gamma_X + k_X,  Gamma_XX = gamma_XX + k_XX.
    Returns dict(I, N, P=(P_G, P_X, P_XX) of the prepared state).

    Caveat (finding 6, documented): this prepares with dissipation OFF
    whereas tpe_photon_statistics runs with it ON, so the TPE pulse factor
    omits in-pulse re-excitation / decay (it is the Gamma tau_p -> 0 cap);
    the g2 and mean photons do include it. The factor is therefore optimistic
    by O(Gamma tau_p) [A]."""
    prep = tpe_prepare(casc, shape, area, duration)
    s = casc.system()
    A = np.sqrt(casc.gamma_X_ns) * s.ops["c_X"]
    I, N = lb.indistinguishability(s.L, s.dim, prep["rho"], A, return_counts=True)
    return {"I": float(I), "N": float(N), "P": prep["P"]}


# ====================================================================== F3
# Spectral diffusion and HOM versus photon separation
#
# Split the line into a fast homogeneous part (pure dephasing gamma*,
# Markovian, in L) and slow spectral diffusion: a Gaussian detuning process
# delta(t), standard deviation sigma_SD, correlation exp(-|dt|/tau_SD)
# (Ornstein-Uhlenbeck) [A]. The detuning is frozen over one photon lifetime
# (quasi-static, valid for tau_SD >> 1/Gamma) [A]. Two photons emitted dt
# apart see detunings (delta_1, delta_2), jointly Gaussian with correlation
# exp(-dt/tau_SD); the difference d = delta_1 - delta_2 is Gaussian with
# variance 2 sigma_SD^2 [1 - exp(-dt/tau_SD)]. Two remote sources are
# uncorrelated (variance 2 sigma_SD^2) plus a static residual offset after
# Stark tuning [A].
#
# Cross-HOM of two independent single photons from emitters 1 and 2 (the
# generalisation of the Grange integral to unequal detuning), HOM
# visibility = coincidence suppression of the 50:50 beam splitter [DR]:
#     I = Re int int G1^(1)(t,t') conj(G1^(2)(t,t')) dt dt' / (N1 N2),
#     G1^(j)(t,t+tau) = Tr[A^+ e^{L_j tau}(A rho_j(t))],  N_j = int Tr[A^+A rho_j] dt,
# and the full-plane integral equals twice the t < t' half. With the stationary
# mode removed (lindblad.indistinguishability), L_j' = L_j - P,
#     Y: L_1' Y + Y L_2'^+ = -v_1 v_2^+,   W: L_1' W + W L_2'^+ = -B Y B^+,
#     I = 2 Re(a^T W a^*) / (N_1 N_2)          (Sylvester equations).
# Both detunings enter only through delta_1 - delta_2: a common shift c
# multiplies G1^(j) by exp(i c tau), which cancels in G1^(1) conj(G1^(2)).
# Closed form for exponential wavepackets with instantaneous excitation [DR]
# (coherence decay (Gamma + gamma*)/2, G1 = Gamma e^{-Gamma t} e^{-(Gamma+gamma*)tau/2}
# e^{-i delta tau}, so int_0^inf e^{-(Gamma+gamma*)tau} e^{-i d tau} dtau =
# 1/(Gamma+gamma*+i d)):
#     I(d) = Gamma (Gamma + gamma*) / ((Gamma + gamma*)^2 + d^2),
# which is |<psi_1|psi_2>|^2 = Gamma^2/(Gamma^2 + d^2) at gamma* = 0. For
# Gaussian d of standard deviation s (mean 0), gamma* = 0:
#     I = sqrt(pi/2) (Gamma/s) erfcx(Gamma/(sqrt 2 s)),
# and with gamma* the same with Gamma -> Gamma + gamma* inside and the
# prefactor Gamma/(Gamma + gamma*) outside.

_HOM_GL_ORDER = 24
_HOM_CHUNK = 20000      # batch size of the Sylvester solves [A] memory choice


def _two_level_cross_pieces(Gamma, gamma_star):
    tl = TwoLevel(Gamma, gamma_star_ns=gamma_star)
    s0 = tl.system(0.0)
    L0 = np.asarray(s0.L, dtype=complex)
    ee = np.zeros((2, 2), dtype=complex)
    ee[1, 1] = 1.0
    Ldet = lb.hamiltonian_super(ee)               # d/d(delta) of L
    A = np.sqrt(Gamma) * s0.ops["c_X"]
    P = lb._stationary_projector(L0)
    return L0, Ldet, A, P


def _hom_cross_batch(Gamma, gamma_star, d1, d2):
    """Numerical cross-HOM I for arrays of detunings (d1, d2) (rad/ns) of two
    independent two-level emitters with equal Gamma, gamma* (module comment
    above): batched Kronecker form of the two Sylvester equations."""
    d1, d2 = np.broadcast_arrays(np.atleast_1d(np.asarray(d1, dtype=float)),
                                 np.atleast_1d(np.asarray(d2, dtype=float)))
    d1, d2 = d1.reshape(-1), d2.reshape(-1)
    if d1.size > _HOM_CHUNK:                       # bound the (nb, 16, 16) work arrays
        return np.concatenate([_hom_cross_batch(Gamma, gamma_star, d1[i:i + _HOM_CHUNK], d2[i:i + _HOM_CHUNK])
                               for i in range(0, d1.size, _HOM_CHUNK)])
    L0, Ldet, A, P = _two_level_cross_pieces(Gamma, gamma_star)
    n = 4
    Q = np.eye(n) - P
    I4 = np.eye(n)
    rho0 = np.zeros((2, 2), dtype=complex)
    rho0[1, 1] = 1.0
    v0 = Q @ lb.vec(rho0)
    B = Q @ lb.sprepost(A, np.eye(2)) @ Q
    a = lb.vec(A.conj())
    b = lb.vec((A.conj().T @ A).T)
    nb = len(d1)
    Lp1 = L0[None] - P[None] + d1[:, None, None] * Ldet[None]
    Lp2 = L0[None] - P[None] + d2[:, None, None] * Ldet[None]
    # K vec(W) = vec(L1' W + W L2'^+), column-stacking vec:
    #   vec(L1' W) = (1 kron L1') vec W,  vec(W L2'^+) = (conj(L2') kron 1) vec W
    K = (np.einsum("ij,nkl->nikjl", I4, Lp1).reshape(nb, n * n, n * n)
         + np.einsum("nij,kl->nikjl", Lp2.conj(), I4).reshape(nb, n * n, n * n))

    def vec_batch(M):          # (nb, n, n) -> (nb, n*n), column stacking
        return np.transpose(M, (0, 2, 1)).reshape(nb, n * n)

    def unvec_batch(v):        # inverse
        return v.reshape(nb, n, n).transpose(0, 2, 1)

    rhs1 = np.broadcast_to(-np.outer(v0, v0.conj())[None], (nb, n, n))
    Ym = unvec_batch(np.linalg.solve(K, vec_batch(rhs1)[..., None])[..., 0])
    BYB = B[None] @ Ym @ B.conj().T[None]
    Wm = unvec_batch(np.linalg.solve(K, -vec_batch(BYB)[..., None])[..., 0])
    S = np.einsum("i,nij,j->n", a, Wm, a.conj())
    vb = np.broadcast_to(v0, (nb, n))[..., None]
    N1 = np.linalg.solve(-Lp1, vb)[..., 0] @ b
    N2 = np.linalg.solve(-Lp2, vb)[..., 0] @ b
    return 2.0 * np.real(S) / np.real(N1 * N2)


def hom_pair(Gamma, gamma_star, d_ns):
    """Cross-HOM visibility of two photons from emitters whose transition
    frequencies differ by d_ns (rad/ns), Gamma the (total) population decay
    rate and gamma_star the pure-dephasing rate (lindblad.py convention),
    instantaneous excitation. Numerical (Sylvester); closed form in the
    module comment. Array-valued d_ns returns an array."""
    d = np.asarray(d_ns, dtype=float)
    out = _hom_cross_batch(float(Gamma), float(gamma_star), d.reshape(-1), np.zeros(d.size))
    return float(out[0]) if d.ndim == 0 else out.reshape(d.shape)


def _gauss_panels(mean, std, scale_pts):
    """Panel edges for the composite Gauss-Legendre rule over
    [mean - 12 std, mean + 12 std]: refined geometrically (ratio 2) around
    the Gaussian peak and every point in ``scale_pts`` (Lorentzian centre)."""
    lo, hi = mean - 12.0 * std, mean + 12.0 * std
    base = min(std, scale_pts[1]) / 4.0
    edges = {lo, hi, mean}
    for c in (mean, scale_pts[0]):
        k = 0
        while base * 2 ** k < (hi - lo):
            for sgn in (-1.0, 1.0):
                x = c + sgn * base * 2 ** k
                if lo < x < hi:
                    edges.add(x)
            k += 1
        if lo < c < hi:
            edges.add(c)
    return np.array(sorted(edges))


def hom_gauss_average(Gamma, gamma_star, mean_ns, std_ns):
    """E[ I(d) ] for d ~ N(mean_ns, std_ns^2) (rad/ns), I(d) = hom_pair.

    Composite Gauss-Legendre quadrature (order 24) on panels refined around
    the Gaussian peak and the Lorentzian centre d = 0. Decision: the work
    order names Gauss-Hermite; for Gamma/s = 0.1 the Lorentzian pole sits
    0.07 sigma from the real axis and plain Gauss-Hermite needs ~1e4 nodes
    for 1e-8, so a panelled rule is used and plain Gauss-Hermite is kept as
    an independent cross-check in the verifier (s7)."""
    std = float(std_ns)
    if std <= 0.0:
        return float(hom_pair(Gamma, gamma_star, float(mean_ns)))
    Gt = float(Gamma) + float(gamma_star)
    edges = _gauss_panels(float(mean_ns), std, (0.0, Gt))
    xg, wg = np.polynomial.legendre.leggauss(_HOM_GL_ORDER)
    half = 0.5 * np.diff(edges)
    mid = 0.5 * (edges[1:] + edges[:-1])
    x = (mid[:, None] + half[:, None] * xg[None, :]).reshape(-1)
    w = (half[:, None] * wg[None, :]).reshape(-1)
    pdf = np.exp(-0.5 * ((x - float(mean_ns)) / std) ** 2) / (std * np.sqrt(2.0 * np.pi))
    return float(np.sum(w * pdf * hom_pair(Gamma, gamma_star, x)))


def hom_vs_separation(dt_ns, Gamma, gamma_star, sigma_sd, tau_sd):
    """HOM visibility of two photons emitted dt_ns apart by the same source
    whose detuning performs an Ornstein-Uhlenbeck walk (sigma_sd in rad/ns,
    correlation time tau_sd in ns) [A]: Gaussian average of the cross-HOM over
    d with variance 2 sigma_sd^2 [1 - exp(-dt/tau_sd)]. dt -> 0 gives the
    single-emitter value Gamma/(Gamma + gamma*); dt >> tau_sd gives
    hom_remote(offset = 0). tau_sd = 0 means uncorrelated at any dt > 0."""
    dt = float(dt_ns)
    if dt <= 0.0:
        var = 0.0
    elif tau_sd <= 0.0:
        var = 2.0 * float(sigma_sd) ** 2
    else:
        var = 2.0 * float(sigma_sd) ** 2 * (1.0 - np.exp(-dt / float(tau_sd)))
    return hom_gauss_average(Gamma, gamma_star, 0.0, np.sqrt(var))


def hom_remote(Gamma, gamma_star, sigma_sd, offset=0.0):
    """HOM visibility of two remote sources: uncorrelated detunings (variance
    2 sigma_sd^2 of the difference) plus a static residual offset (rad/ns)
    after Stark tuning [A]."""
    return hom_gauss_average(Gamma, gamma_star, float(offset), np.sqrt(2.0) * float(sigma_sd))


# ====================================================================== F4
# Boson-sampling figure of merit, card loader
#
# Efficiency chain [DR]: eta = eta_source * eta_demux * eta_circuit * eta_det,
#   eta_source = n_ph * beta * eta_out * Z
# with n_ph the mean photons per pulse in the selected line (pulsed
# statistics), beta = (beta0 + F_P L)/(1 + F_P L) the fraction of the emission
# that reaches the collected mode (beta0 = bare-dot collection fraction [A];
# F_P L the Purcell-enhanced cavity channel, L = spectral overlap of
# ``purcell_rate``), eta_out the single-mode extraction/transfer efficiency
# of the collected mode [E/A], and Z = qd_gf.zpl_weight(T) when the card
# selects the ZPL-filtered output (1 otherwise). The unfiltered output keeps the
# sideband photons in eta but multiplies the off-diagonal HOM by Z^2
# (``hom_phonon_factor``); the filtered output has eta x Z and no HOM penalty.
# N-fold coincidence rate for time demultiplexing ONE source into N modes [DR]:
#     R_N = (f_rep / N) eta^N
# (one N-photon event needs N consecutive pulses of the source; every photon
# survives with probability eta).
# Per-photon g2(0) from the pulsed statistics (F1 / F2); the pairwise HOM
# matrix from F3 with dt_ij = |i - j| / f_rep, optionally times the
# single-pulse factor (re-excitation, cascade) -- see ``bs_figure_of_merit``.


@dataclass(frozen=True)
class FockSource:
    """One concrete operating point of the Fock-state source (all numbers are
    inputs from the card; units in the field names). ``provenance`` is a tuple
    of (card entry, tag) pairs."""
    transition: str = "X+"
    excitation: str = "resonant_pi"
    T_K: float = 4.0
    gamma0_ns: float = 1.5
    F_P: float = 5.0
    kappa_ueV: float = 150.0
    cav_detune_ueV: float = 0.0
    pulse_shape: str = "gauss"
    pulse_fwhm_ps: float = 10.0
    hom_width_multiple: float = 1.5
    sigma_sd_ueV: float = 0.3
    tau_sd_ns: float = 100.0
    remote_offset_ueV: float = 0.1
    delta_xx_meV: float = -2.0
    gamma_xx_ratio: float = 2.0
    f_rep_MHz: float = 80.0
    beta0: float = 0.1
    eta_out: float = 0.5
    eta_demux: float = 0.8
    eta_circuit: float = 0.8
    eta_det: float = 0.9
    zpl_output: bool = False
    phonon_alpha_ps2: float = 0.0181
    dephasing_mode: str = "enhanced"
    provenance: tuple = ()

    def __post_init__(self):
        if self.transition not in ("X", "X+"):
            raise ValueError("transition must be 'X' or 'X+'")
        if self.excitation not in ("resonant_pi", "tpe"):
            raise ValueError("excitation must be 'resonant_pi' or 'tpe'")
        if self.excitation == "tpe" and self.transition != "X":
            raise ValueError("two-photon excitation needs the neutral X-XX cascade (transition 'X')")
        if self.pulse_shape not in ("rect", "gauss"):
            raise ValueError("pulse_shape must be 'rect' or 'gauss'")
        if self.dephasing_mode not in ("enhanced", "dot"):
            raise ValueError("dephasing_mode must be 'enhanced' or 'dot'")
        if self.hom_width_multiple < 1.0:
            raise ValueError("hom_width_multiple >= 1 (the transform limit is 1)")


def _rate_from_ueV(e_ueV):
    return float(lb.mev_to_rate(float(e_ueV) * 1e-3))


@lru_cache(maxsize=None)
def _pulse_physics(excitation, shape, Gamma, gamma_star, tau_ns, f_rep_MHz, delta_xx_meV, gxx_ratio, tol,
                   gamma0_ns):
    """g2, mean photons and the single-pulse HOM factor of one pulsed
    configuration (cached; depends on nothing else).

    TPE: the X line (rate ``Gamma``) sits in the cavity but the XX line is
    displaced by |Delta_XX| = 1-4 meV >> kappa = 50-300 ueV, so the biexciton
    decay is NOT Purcell-enhanced: gamma_XX = gxx_ratio * Gamma_0 (bare)
    [A] (Fable review finding 2). The cascade cap is then
    Gamma_XX/(Gamma_XX + Gamma_X) = 2 G0/(2 G0 + G0 (1 + F_P L)) (0.25 at
    F_P L = 5, 0.154 at 10).

    Pulse factor: I_p / inst is the first-order-coherence (single-photon
    wavepacket) ratio; it ignores the two-photon component, which lowers the
    MEASURED two-photon-interference visibility by O(g2) (Fischer et al.,
    Phys. Rev. Lett. 119, 2017; Ollivier et al., Phys. Rev. Lett. 126, 063602
    (2021)) [E]; the F4 verdict therefore uses the conservative I - g2.
    No dephasing is double counted: re-excitation enters I_p, dephasing only
    through gamma_star in the same Liouvillian."""
    inst = Gamma / (Gamma + gamma_star)               # (d1) instantaneous-excitation value
    if excitation == "resonant_pi":
        tl = TwoLevel(Gamma, gamma_star_ns=gamma_star)
        st = pulsed_photon_statistics(tl, shape, np.pi, tau_ns, f_rep_MHz)
        I_p = pulsed_indistinguishability(tl, shape, np.pi, tau_ns, tol=tol, max_doublings=4)
        return st["g2"], st["mean_photons"], I_p / inst
    delta = tpe_detuning_ns(delta_xx_meV)
    casc = Cascade(Gamma, gxx_ratio * gamma0_ns, deph=((gamma_star, (0.0, 1.0, 0.0)),) if gamma_star > 0 else (),
                   delta_ns=delta)
    area = _tpe_area_cached(delta, shape, tau_ns)
    st = tpe_photon_statistics(casc, shape, area, tau_ns, f_rep_MHz)
    hom = tpe_cascade_hom(casc, shape, area, tau_ns)
    return st["g2"], st["mean_photons"], hom["I"] / inst


@lru_cache(maxsize=None)
def _hom_sep(dt_ns, Gamma, gamma_star, sigma, tau_sd):
    return hom_vs_separation(dt_ns, Gamma, gamma_star, sigma, tau_sd)


def source_rates(src):
    """Derived rates of a FockSource (1/ns): Gamma (Purcell-enhanced total
    radiative decay), gamma_star (pure dephasing), L (cavity overlap), Gamma0.

    Homogeneous width = hom_width_multiple m times the transform limit,
    angular Gamma + gamma* = m Gamma_ref [DR]; the V-a Gamma_0 = 119 ueV is
    NOT used. ``dephasing_mode`` [A] selects what Gamma_ref is (sweep both,
    never average; Fable review finding 3):

      'enhanced'  Gamma_ref = Gamma (Purcell-enhanced): gamma* = (m - 1) Gamma.
                  I = Gamma/(Gamma + gamma*) = 1/m is independent of F_P BY
                  CONSTRUCTION (the dephasing is assumed to scale with the
                  rate; there is no Purcell protection). Pessimistic reading
                  of a fixed linewidth multiple.
      'dot'       Gamma_ref = Gamma_0 (bare): gamma* = (m - 1) Gamma_0 is a
                  property of the dot (phonon/charge environment), not of the
                  cavity, so I = Gamma/(Gamma + gamma*) rises with F_P
                  (Grange et al., PRL 114, 193601 (2015); Iles-Smith et al.,
                  Nat. Photonics 11, 521 (2017)).

    The Purcell overlap L uses the SAME line FWHM that the dephasing implies,
    hbar (Gamma + gamma*), self-consistently (finding 8): L = T(w(L)),
    T = spectral.cavity_transmission(detuning, w, kappa) and
      enhanced: w = hbar m Gamma_0 (1 + F_P L),
      dot:      w = hbar Gamma_0 (m + F_P L),
    solved by Brent's method on L in (0, 1] (the root is unique: T decreases
    in w and w increases in L) [DR]."""
    from scipy.optimize import brentq
    h = lb.HBAR_MEV_NS * 1e3                         # ueV * ns
    m, G0, F = src.hom_width_multiple, src.gamma0_ns, src.F_P
    b = m * F if src.dephasing_mode == "enhanced" else F

    def width(L):
        return h * G0 * (m + b * L)

    def resid(L):
        return float(cavity_transmission(src.cav_detune_ueV, width(L), src.kappa_ueV)) - L

    L = float(brentq(resid, 0.0, 1.0, xtol=1e-15, rtol=8.9e-16))
    Gamma = G0 * (1.0 + F * L)
    ref = Gamma if src.dephasing_mode == "enhanced" else G0
    return {"Gamma": Gamma, "gamma_star": (m - 1.0) * ref, "L": L, "Gamma0": G0}


def phonon_Z(src):
    """Debye-Waller ZPL weight Z(T) of the source's phonon coupling
    (qd_gf.zpl_weight, alpha_ps2 from the card) [E], independent of whether the
    output is ZPL-filtered."""
    from fsim_core.qd_gf import PhononParams, zpl_weight
    return float(zpl_weight(PhononParams(alpha_ps2=src.phonon_alpha_ps2), src.T_K))


def zpl_weight_for(src):
    """Efficiency-chain Z: phonon_Z(src) when the card selects the
    ZPL-filtered output (sideband photons are discarded: eta x Z), 1 otherwise
    (sideband photons are kept in eta, but they are distinguishable, see
    ``hom_phonon_factor``)."""
    if not src.zpl_output:
        return 1.0
    return phonon_Z(src)


def hom_phonon_factor(src):
    """Phonon factor on the off-diagonal HOM [DR/E]. Weak-coupling
    independent-boson theory (Iles-Smith et al., Nat. Photonics 11, 521
    (2017); Grange et al., PRL 114, 193601 (2015)): the ZPL photons carry the
    coherence, the phonon-sideband photons are distinguishable, so the
    UNFILTERED two-photon visibility is I ~= Z^2 I_ZPL, while a ZPL filter
    restores I_ZPL at the price eta x Z. Returns Z^2 for zpl_output False and
    1 for True (Fable review finding 1; Z from the same PhononParams(alpha)
    as the efficiency chain)."""
    return 1.0 if src.zpl_output else phonon_Z(src) ** 2


def end_to_end_eta(source, n_ph, eta_chain=None):
    """Per-photon end-to-end efficiency eta = n_ph beta eta_out Z eta_demux
    eta_circuit eta_det (efficiency-chain comment above). Returns (eta,
    components dict)."""
    ch = {"demux": source.eta_demux, "circuit": source.eta_circuit, "det": source.eta_det}
    if eta_chain:
        ch.update({k: float(v) for k, v in eta_chain.items()})
    r = source_rates(source)
    beta = (source.beta0 + source.F_P * r["L"]) / (1.0 + source.F_P * r["L"])
    Z = zpl_weight_for(source)
    eta_source = n_ph * beta * source.eta_out * Z
    eta = eta_source * ch["demux"] * ch["circuit"] * ch["det"]
    return float(eta), {"source": float(eta_source), "demux": ch["demux"], "circuit": ch["circuit"],
                        "det": ch["det"], "n_ph": float(n_ph), "beta": float(beta),
                        "eta_out": source.eta_out, "Z": Z}


def bs_figure_of_merit(source, N, f_rep_MHz=None, eta_chain=None, tol=1e-6):
    """Boson-sampling figure of merit of a time-demultiplexed Fock source.

    source     FockSource
    N          number of photons entering the network
    f_rep_MHz  repetition rate (default source.f_rep_MHz)
    eta_chain  optional dict overriding eta_demux / eta_circuit / eta_det
    Returns dict with
      N, f_rep_MHz, g2 (list of N per-photon g2(0)), eta (end-to-end per
      photon), eta_components, hom_matrix (N x N nested list; symmetric, unit
      diagonal by definition, off-diagonals F3 ``hom_vs_separation`` at
      dt_ij = |i-j| 1e3/f_rep_MHz ns), hom_matrix_pulse_corrected (same
      times the single-pulse factor: re-excitation for resonant pi pulses,
      the cascade cap for TPE; the factorisation is an [A] assumption),
      pulse_factor, hom_phonon_factor (Z^2 unfiltered, 1 ZPL-filtered; already
      in both matrices), mean_pairwise_I, mean_pairwise_I_corrected,
      mean_pairwise_I_conservative (corrected I - g2: the two-photon component
      lowers the measured visibility by O(g2), the verdict quantity), R_N_Hz
      (= f_rep/N * eta^N with f_rep in Hz), provenance."""
    f = float(source.f_rep_MHz if f_rep_MHz is None else f_rep_MHz)
    r = source_rates(source)
    tau_ns = source.pulse_fwhm_ps * 1e-3
    g2, n_ph, pf = _pulse_physics(source.excitation, source.pulse_shape, r["Gamma"], r["gamma_star"],
                                  tau_ns, f, source.delta_xx_meV, source.gamma_xx_ratio, tol, source.gamma0_ns)
    eta, comps = end_to_end_eta(source, n_ph, eta_chain)
    sigma = _rate_from_ueV(source.sigma_sd_ueV)
    dt_unit = 1.0e3 / f
    H = np.ones((N, N))
    for i in range(N):
        for j in range(i + 1, N):
            H[i, j] = H[j, i] = _hom_sep(round((j - i) * dt_unit, 9), r["Gamma"], r["gamma_star"],
                                         sigma, float(source.tau_sd_ns))
    off = ~np.eye(N, dtype=bool)
    zf = hom_phonon_factor(source)
    H[off] = H[off] * zf                               # unfiltered: Z^2; ZPL-filtered: 1
    Hc = H.copy()
    Hc[off] = H[off] * pf
    mean = float(H[off].mean()) if N > 1 else 1.0
    meanc = float(Hc[off].mean()) if N > 1 else 1.0
    return {"N": int(N), "f_rep_MHz": f, "g2": [float(g2)] * int(N), "eta": float(eta),
            "eta_components": comps,
            "hom_matrix": H.tolist(), "hom_matrix_pulse_corrected": Hc.tolist(),
            "pulse_factor": float(pf), "hom_phonon_factor": float(zf), "mean_pairwise_I": mean,
            "mean_pairwise_I_corrected": meanc,
            "mean_pairwise_I_conservative": float(meanc - g2),
            "R_N_Hz": coincidence_rate_Hz(f, N, eta),
            "provenance": dict(source.provenance)}


def coincidence_rate_Hz(f_rep_MHz, N, eta):
    """N-fold coincidence rate of time demultiplexing ONE source,
    R_N = (f_rep/N) eta^N [DR], in Hz (f_rep given in MHz)."""
    return float(f_rep_MHz) * 1.0e6 / int(N) * float(eta) ** int(N)


# ------------------------------------------------------------------ card

_CARD_TAGS = ("V", "DR", "E", "A")


@dataclass
class FockParam:
    """One card entry: exactly one of ``value`` (fixed), ``range`` (swept,
    never averaged) or ``choices`` (+ ``default``); always unit, tag, source."""
    name: str
    unit: str
    tag: str
    source: str
    value: object = None
    range: tuple | None = None
    choices: tuple | None = None
    default: object = None
    scale: str = "lin"

    def points(self, n=3):
        """Sweep points: ``n`` points over a range (lin or log spacing), the
        single value, or every choice."""
        if self.range is not None:
            lo, hi = self.range
            if self.scale == "log":
                return [float(v) for v in np.geomspace(lo, hi, n)]
            return [float(v) for v in np.linspace(lo, hi, n)]
        if self.choices is not None:
            return list(self.choices)
        return [self.value]

    def central(self):
        if self.range is not None:
            lo, hi = self.range
            return float(np.sqrt(lo * hi)) if self.scale == "log" else 0.5 * (lo + hi)
        if self.choices is not None:
            return self.default
        return self.value

    def to_dict(self):
        d = {"unit": self.unit, "tag": self.tag, "source": self.source}
        if self.range is not None:
            d["range"] = [float(self.range[0]), float(self.range[1])]
            if self.scale != "lin":
                d["scale"] = self.scale
        elif self.choices is not None:
            d["choices"] = list(self.choices)
            d["default"] = self.default
        else:
            d["value"] = self.value
        return d


@dataclass
class FockCard:
    meta: dict
    params: dict

    def to_dict(self):
        return {"meta": dict(self.meta), "params": {n: p.to_dict() for n, p in self.params.items()}}

    def __getitem__(self, name):
        return self.params[name]

    def provenance(self):
        return tuple((n, p.tag) for n, p in self.params.items())


def _parse_fock_param(name, raw):
    for key in ("unit", "tag", "source"):
        if key not in raw or raw[key] in (None, ""):
            raise ValueError(f"card entry {name}: missing required field '{key}'")
    if str(raw["tag"]) not in _CARD_TAGS:
        raise ValueError(f"card entry {name}: tag must be one of {_CARD_TAGS}")
    kinds = [k for k in ("value", "range", "choices") if k in raw]
    if len(kinds) != 1:
        raise ValueError(f"card entry {name}: exactly one of value/range/choices required, got {kinds}")
    p = FockParam(name=name, unit=str(raw["unit"]), tag=str(raw["tag"]), source=str(raw["source"]))
    if "range" in raw:
        lo, hi = (float(raw["range"][0]), float(raw["range"][1]))
        if not lo < hi:
            raise ValueError(f"card entry {name}: range must satisfy lo < hi")
        p.range = (lo, hi)
        p.scale = str(raw.get("scale", "lin"))
        if p.scale not in ("lin", "log"):
            raise ValueError(f"card entry {name}: scale must be lin or log")
        if p.scale == "log" and lo <= 0:
            raise ValueError(f"card entry {name}: log scale needs lo > 0")
    elif "choices" in raw:
        p.choices = tuple(raw["choices"])
        if raw.get("default") not in p.choices:
            raise ValueError(f"card entry {name}: default must be one of the choices")
        p.default = raw["default"]
    else:
        p.value = raw["value"]
    return p


def load_fock_card(path):
    """Load a Fock-source card (YAML): {meta: {...}, params: {name: {unit,
    tag, source, value | range [+scale] | choices + default}}}. Every leaf
    must carry unit, tag and source. Ranges are swept, never averaged."""
    import yaml
    from pathlib import Path
    doc = yaml.safe_load(Path(path).read_text(encoding="utf-8"))
    if "meta" not in doc or "name" not in doc["meta"]:
        raise ValueError("card needs meta.name")
    params = {n: _parse_fock_param(n, raw) for n, raw in doc.get("params", {}).items()}
    return FockCard(meta=dict(doc["meta"]), params=params)


def fock_card_from_dict(doc):
    return FockCard(meta=dict(doc["meta"]),
                    params={n: _parse_fock_param(n, raw) for n, raw in doc["params"].items()})


def dump_fock_card(card, path):
    """Write ``card`` as YAML (round-trips through ``load_fock_card``)."""
    import yaml
    from pathlib import Path
    Path(path).write_text(yaml.safe_dump(card.to_dict(), sort_keys=False, allow_unicode=True), encoding="utf-8")


_SOURCE_FIELDS = {
    "T_K": "T_K", "transition": "transition", "excitation": "excitation", "gamma0_ns": "gamma0_ns",
    "F_P": "F_P", "kappa_ueV": "kappa_ueV", "cav_detune_ueV": "cav_detune_ueV",
    "pulse_shape": "pulse_shape", "pulse_fwhm_ps": "pulse_fwhm_ps",
    "hom_width_multiple": "hom_width_multiple", "sigma_sd_ueV": "sigma_sd_ueV",
    "tau_sd_ns": "tau_sd_ns", "remote_offset_ueV": "remote_offset_ueV",
    "delta_xx_meV": "delta_xx_meV", "gamma_xx_ratio": "gamma_xx_ratio", "f_rep_MHz": "f_rep_MHz",
    "beta0": "beta0", "eta_out": "eta_out", "eta_demux": "eta_demux", "eta_circuit": "eta_circuit",
    "eta_det": "eta_det", "zpl_output": "zpl_output", "phonon_alpha_ps2": "phonon_alpha_ps2",
    "dephasing_mode": "dephasing_mode",
}


def source_from_card(card, **overrides):
    """FockSource at the card's central point (value, range midpoint, choice
    default) with ``overrides`` (card entry name -> value) applied. The
    central point is a convenience for a single evaluation; sweeps use
    ``FockParam.points`` and never average a range."""
    kw = {}
    for name, fld in _SOURCE_FIELDS.items():
        if name in overrides:
            kw[fld] = overrides[name]
        elif name in card.params:
            kw[fld] = card.params[name].central()
    unknown = set(overrides) - set(_SOURCE_FIELDS)
    if unknown:
        raise KeyError(f"unknown card entries: {sorted(unknown)}")
    return FockSource(provenance=card.provenance(), **kw)
