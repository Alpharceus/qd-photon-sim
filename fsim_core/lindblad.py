"""Hand-written Lindblad master-equation core (numpy/scipy only, no third-party quantum toolbox).

Optional quantum layer for the qd-photon-sim F-series. Nothing in the legacy
rate-equation paths imports this module, so verify_fsim (51/51) and
audit_physics (23/23) are untouched by construction. Design and seed:
internal quantum audit (Phase A, 2026-09-22). Regression suite:
verify/verify_lindblad.py.

Units
-----
Rates and angular frequencies in 1/ns, times in ns, energies in meV.
E [meV] -> angular frequency E / HBAR_MEV_NS [1/ns]  (helper ``mev_to_rate``).
A line FWHM of F meV corresponds to an amplitude (coherence) decay rate
F / (2 HBAR_MEV_NS) and a Lorentzian FWHM F / HBAR_MEV_NS in angular units.

Equations as implemented
------------------------
Master equation (Lindblad form; Gorini-Kossakowski-Sudarshan, J. Math. Phys.
17, 821 (1976); Lindblad, Commun. Math. Phys. 48, 119 (1976)):

    d rho/dt = L rho = -i [H, rho] + sum_c ( c rho c^+ - 1/2 {c^+ c, rho} ).

Vectorisation is column stacking, vec(A X B) = (B^T kron A) vec(X), so
    spre(A)      = I kron A          (X -> A X)
    spost(A)     = A^T kron I        (X -> X A)
    sprepost(A,B)= B^T kron A        (X -> A X B)
    D[c]         = sprepost(c, c^+) - 1/2 spre(c^+c) - 1/2 spost(c^+c).
Tr[A X] = vec(A^T) . vec(X).

Emitter basis (``levels=3``): |G>=0, |X>=1, |XX>=2 (``levels=2`` drops XX).
Optional charged state |X+> (``charged=True``) is appended as the last
emitter level; it is fed from |G> at ``r_charge_ns`` and relaxes to |G> at
``gamma_charge_ns`` (radiative, a photon channel of its own) -- a minimal
blinking/charged-trion channel [A].
Dissipators (all [DR], standard quantum-optical master-equation terms,
e.g. Breuer & Petruccione, The Theory of Open Quantum Systems, OUP 2002):
    D[sqrt(r)        |X><G| ]      incoherent pump G -> X
    D[sqrt(p r)      |XX><X|]      incoherent pump X -> XX (p = pump_ratio)
    D[sqrt(gamma_X)  |G><X| ]      X radiative (photon channel c_X)
    D[sqrt(gamma_XX) |X><XX|]      XX radiative (photon channel c_XX)
    D[sqrt(k_X)      |G><X| ]      X escape (no photon)
    D[sqrt(k_XX)     |X><XX|]      XX escape (no photon)
    D[sqrt(c_j) sum_i v_ji |i><i|] pure dephasing (``deph`` list of (c_j, v_j));
                                   the i-k coherence decays at an extra
                                   (1/2) sum_j c_j (v_ji - v_jk)^2.
    D[sqrt(kappa) a]               cavity escape (optional cavity mode a)
Hamiltonian (rotating frame of the X transition, so ``det_X`` is the X
transition's detuning from the frame, typically 0):
    H = det_X |X><X| + (det_X + det_XX) |XX><XX|
        + (Omega/2) (|X><G| + |G><X|)                  coherent drive (Rabi Omega)
        + g (|G><X| a^+ + |X><G| a) + det_c a^+ a      Jaynes-Cummings cavity
Cascaded Lorentzian filter (sensor) mode f of FWHM w, driven without
back-action by the emitted field (Gardiner, PRL 70, 2269 (1993);
Carmichael, PRL 70, 2273 (1993)); two-port filter, input and output ports
each w/2, so the peak power transmission is 1 and the filter FWHM is w:
    source field  s = sqrt(gamma_X) c_X + sqrt(gamma_XX) c_XX (radiative only)
    collapse      C = s + sqrt(w/2) f      (replaces D[sqrt(gamma_X) c_X] etc.)
    collapse      sqrt(w/2) f               (detected output port)
    H_casc = (i/2) sqrt(w/2) (s^+ f - f^+ s) + det_f f^+ f
(det_f is the filter centre relative to the frame.) The filtered photon
stream has rate (w/2) <f^+ f> and zero-delay g2_f(0) = <f^+2 f^2>/<f^+ f>^2.
Because <f^+ f> is small (~1e-8 at the device points), the steady state is
solved in a similarity-scaled basis rho_nm -> rho_nm / s^(n+m) (n, m the
sensor photon numbers), as in the audit prototype quantum_filtered_g2.py.

Steady state: L rho = 0 with one row replaced by Tr rho = 1.
Time evolution: exp(L t) through one eigendecomposition (checked against
scipy.linalg.expm, and falling back to expm when the eigenbasis is
ill-conditioned); time-dependent H(t) = H0 + f(t) H1 by piecewise-constant
midpoint products of expm.
Two-time correlations (quantum regression theorem, Lax 1963; Carmichael,
Statistical Methods in Quantum Optics 1, Springer 1999):
    G2_ij(tau) = Tr[ c_j^+ c_j e^{L tau} (c_i rho_ss c_i^+) ],
    g2(tau) = sum_ij w_i w_j G2_ij(tau) / (sum_i w_i <c_i^+ c_i>)^2,
with classical per-photon detection probabilities w_i (Bernoulli thinning,
the cw_g2 convention).
Pulsed photon counting (factorial moments, the pulse_counting hierarchy
with the Liouvillian in place of the rate matrix M):
    d p/dt = L p,  d m1/dt = L m1 + J p,  d m2/dt = L m2 + 2 J m1,
    J rho = sum_i w_i gamma_i c_i rho c_i^+,   g2 = Tr m2 / (Tr m1)^2.
HOM indistinguishability for a single excitation (Grange et al., PRL 114,
193601 (2015), Eq. 2; Bylander, Robert-Philip & Abram, Eur. Phys. J. D 22,
295 (2003)):
    I = 2 int_0^inf dt int_0^inf dtau |<A^+(t+tau) A(t)>|^2 / (int_0^inf <A^+A> dt)^2,
    <A^+(t+tau) A(t)> = Tr[A^+ e^{L tau}(A rho(t))],
evaluated exactly by two Lyapunov (Sylvester) solves on L with the stationary
mode projected out; no diagonalizability needed (defective cascaded filters
are fine). The eigenmode Laplace-sum path is kept as a private reference.

Validity: Markovian baths only (Lorentzian ZPL). Non-Markovian phonon
sidebands (fsim_core.qd_gf) are a separate lineshape layer [A].

Public API
----------
Superoperators : spre, spost, sprepost, vec, unvec, dissipator,
                 liouvillian, jump_super
Model          : build_system(...) -> LindbladSystem, dephasing_for_fwhm(...)
Solvers        : steady_state, evolve, evolve_td, Propagator
Observables    : expect, g2_tau, g2_zero, filtered_g2_zero, filtered_flux,
                 pulsed_counting, deterministic_cycle_counting,
                 indistinguishability, check_physical
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import eig, expm, solve, solve_continuous_lyapunov

# hbar in meV*ns. [V] CODATA 2018 (Tiesinga et al., Rev. Mod. Phys. 93,
# 025010 (2021)): hbar = 6.582119569e-16 eV s (exact in the 2019 SI).
HBAR_MEV_NS = 6.582119569e-4
# Tolerances used by the numerical guards below. [A] numerical choices.
_EIG_COND_MAX = 1e6           # eigenbasis condition number above which expm is used
_EIG_RESID_MAX = 1e-12        # max ||V diag(lam) V^-1 - L|| / ||L|| for the eig path
_ZERO_MODE_REL = 1e-12        # |lambda| below this x max|lambda| is a stationary mode
_PERIODIC_TOL = 1e-13         # fixed-point tolerance of the one-period map
_PERIODIC_MAX_ITER = 20000


def mev_to_rate(e_mev):
    """Energy (meV) -> angular frequency (1/ns). [DR] E/hbar."""
    return np.asarray(e_mev, dtype=float) / HBAR_MEV_NS


# ------------------------------------------------------------ superoperators

def spre(A):
    """Superoperator of X -> A X (column-stacking vec)."""
    return np.kron(np.eye(A.shape[0]), A)


def spost(A):
    """Superoperator of X -> X A."""
    return np.kron(A.T, np.eye(A.shape[0]))


def sprepost(A, B):
    """Superoperator of X -> A X B."""
    return np.kron(B.T, A)


def vec(X):
    """Column-stacked vector of a matrix."""
    return np.asarray(X, dtype=complex).reshape(-1, order="F")


def unvec(v, n):
    """Inverse of vec for an n x n matrix."""
    return np.asarray(v).reshape((n, n), order="F")


def dissipator(c):
    """D[c] = c . c^+ - 1/2 {c^+ c, .}."""
    c = np.asarray(c, dtype=complex)
    cdc = c.conj().T @ c
    return sprepost(c, c.conj().T) - 0.5 * spre(cdc) - 0.5 * spost(cdc)


def liouvillian(H, c_ops):
    """L = -i[H, .] + sum_c D[c]."""
    H = np.asarray(H, dtype=complex)
    L = -1j * (spre(H) - spost(H))
    for c in c_ops:
        L = L + dissipator(c)
    return L


def hamiltonian_super(H):
    """-i[H, .] alone (for time-dependent drive terms)."""
    H = np.asarray(H, dtype=complex)
    return -1j * (spre(H) - spost(H))


def jump_super(c):
    """J rho = c rho c^+."""
    c = np.asarray(c, dtype=complex)
    return sprepost(c, c.conj().T)


def trace_row(n):
    """Row vector t with t . vec(X) = Tr X."""
    return vec(np.eye(n)).conj()


def expect(op, rho):
    """Tr[op rho]."""
    return complex(np.trace(np.asarray(op) @ np.asarray(rho)))


# ------------------------------------------------------------ model building

@dataclass
class LindbladSystem:
    """A built model: Liouvillian plus the operators needed downstream.

    Attributes
    ----------
    L : (d^2, d^2) complex ndarray, the Liouvillian.
    H : (d, d) Hamiltonian (time-independent part).
    c_ops : list of collapse operators.
    dim : Hilbert-space dimension d.
    ops : dict of named operators ("P_G", "P_X", "P_XX", "c_X", "c_XX",
          "a" (cavity), "f" (filter), "sigma_x_drive", ...). Emitter
          operators are embedded in the full space.
    photon : dict name -> (operator c, rate gamma) of the radiative photon
          channels, so that the photon flux of channel i is gamma <c^+ c>.
    sensor_number : (d,) int array, the filter photon number of each basis
          state (all zero without a filter); used for similarity scaling.
    """
    L: np.ndarray
    H: np.ndarray
    c_ops: list
    dim: int
    ops: dict = field(default_factory=dict)
    photon: dict = field(default_factory=dict)
    sensor_number: np.ndarray | None = None
    meta: dict = field(default_factory=dict)


def dephasing_for_fwhm(fwhm_X_ns, fwhm_XX_ns=None, *, r_ns=0.0, gamma_X_ns=0.0,
                       gamma_XX_ns=0.0, k_X=0.0, k_XX=0.0, pump_ratio=1.0,
                       mode="corr", levels=3):
    """Pure-dephasing list reproducing target line FWHMs (angular, 1/ns).

    The population outflows already broaden the lines [DR]:
        X-G  line FWHM_base = r (1 + p) + gamma_X + k_X          (levels=3)
                            = r + gamma_X + k_X                  (levels=2)
        XX-X line FWHM_base = gamma_XX + k_XX + gamma_X + k_X + p r.
    The remainder is assigned to pure dephasing, following the audit
    prototype quantum_filtered_g2.build [A, section 1 of the audit]:
      mode="corr": one fluctuating field, D[sqrt(c) (|X><X| + s |XX><XX|)],
                   c = FWHM_X - base_X, s = 1 + sqrt((FWHM_XX - base_XX)/c);
      mode="indep": D[sqrt(c)|X><X|] and D[sqrt(c')|XX><XX|] with
                   c' = FWHM_XX - base_XX - c (independent level noise).
    Negative remainders are clipped to 0 (lifetime-limited line). Returns a
    list of (rate, diagonal-vector) pairs for build_system(deph=...)."""
    p = pump_ratio
    if levels == 2:
        base_X = r_ns + gamma_X_ns + k_X
        c = max(fwhm_X_ns - base_X, 0.0)
        return [(c, (0.0, 1.0))] if c > 0 else []
    base_X = r_ns * (1.0 + p) + gamma_X_ns + k_X
    base_XX = gamma_XX_ns + k_XX + gamma_X_ns + k_X + p * r_ns
    c = max(fwhm_X_ns - base_X, 0.0)
    if fwhm_XX_ns is None:
        fwhm_XX_ns = fwhm_X_ns
    exXX = max(fwhm_XX_ns - base_XX, 0.0)
    if mode == "corr":
        s = 1.0 + np.sqrt(exXX / c) if c > 0 else 2.0
        return [(c, (0.0, 1.0, s))] if c > 0 else []
    if mode == "indep":
        out = []
        if c > 0:
            out.append((c, (0.0, 1.0, 0.0)))
        cxx = max(exXX - c, 0.0)
        if cxx > 0:
            out.append((cxx, (0.0, 0.0, 1.0)))
        return out
    raise ValueError("mode must be 'corr' or 'indep'")


def _destroy(n):
    return np.diag(np.sqrt(np.arange(1, n)), 1).astype(complex)


def build_system(*, levels=3, r_ns=0.0, gamma_X_ns=1.0, gamma_XX_ns=0.0,
                 k_X=0.0, k_XX=0.0, pump_ratio=1.0, deph=(),
                 omega_ns=0.0, det_X_ns=0.0, det_XX_ns=0.0,
                 charged=False, r_charge_ns=0.0, gamma_charge_ns=0.0,
                 cavity=None, filt=None):
    """Build the emitter [+ cavity] [+ cascaded filter] Liouvillian.

    Parameters (rates 1/ns; detunings angular 1/ns):
      levels      2 (G, X) or 3 (G, X, XX)
      r_ns, pump_ratio, gamma_X_ns, gamma_XX_ns, k_X, k_XX : see module doc
      deph        iterable of (rate, diag_vector) pure-dephasing terms over the
                  emitter levels (see dephasing_for_fwhm)
      omega_ns    coherent Rabi frequency Omega on G<->X (H = Omega/2 sigma_x)
      det_X_ns, det_XX_ns : transition detunings from the frame
      charged     append a charged level |X+> (fed from G at r_charge_ns,
                  radiative back to G at gamma_charge_ns)
      cavity      None or dict(g_ns, kappa_ns, det_ns=0.0, n_max=2): a mode a
                  coupled to the G<->X transition; n_max = number of Fock
                  states kept (0..n_max-1)
      filt        None or dict(fwhm_ns, det_ns=0.0, n_max=3, source="all"):
                  cascaded two-port Lorentzian filter driven by the radiative
                  field ("all" = X + XX (+X+), "X" = X line only)
    """
    if levels not in (2, 3):
        raise ValueError("levels must be 2 or 3")
    ne = levels + (1 if charged else 0)
    nc = int(cavity.get("n_max", 2)) if cavity else 1
    nf = int(filt.get("n_max", 3)) if filt else 1
    dim = ne * nc * nf
    Ic, If = np.eye(nc), np.eye(nf)

    def emb_e(m):
        return np.kron(np.kron(m, Ic), If)

    def e(i, j):
        m = np.zeros((ne, ne), dtype=complex)
        m[i, j] = 1.0
        return emb_e(m)

    G, X = 0, 1
    XX = 2 if levels == 3 else None
    XP = ne - 1 if charged else None

    ops = {"P_G": e(G, G), "P_X": e(X, X)}
    if XX is not None:
        ops["P_XX"] = e(XX, XX)
    if XP is not None:
        ops["P_Xp"] = e(XP, XP)
    ops["c_X"] = e(G, X)
    if XX is not None:
        ops["c_XX"] = e(X, XX)
    if XP is not None:
        ops["c_Xp"] = e(G, XP)

    H = np.zeros((dim, dim), dtype=complex)
    H += det_X_ns * e(X, X)
    if XX is not None:
        H += (det_X_ns + det_XX_ns) * e(XX, XX)
    if omega_ns:
        sx = e(X, G) + e(G, X)
        ops["sigma_x_drive"] = sx
        H += 0.5 * omega_ns * sx

    c_ops = []
    photon = {}
    # incoherent pump
    if r_ns:
        c_ops.append(np.sqrt(r_ns) * e(X, G))
        if XX is not None and pump_ratio * r_ns:
            c_ops.append(np.sqrt(pump_ratio * r_ns) * e(XX, X))
    if XP is not None and r_charge_ns:
        c_ops.append(np.sqrt(r_charge_ns) * e(XP, G))
    # escape (no photon)
    if k_X:
        c_ops.append(np.sqrt(k_X) * e(G, X))
    if XX is not None and k_XX:
        c_ops.append(np.sqrt(k_XX) * e(X, XX))
    # pure dephasing
    for rate, v in deph:
        if rate <= 0:
            continue
        v = np.asarray(v, dtype=float)
        if v.size < ne:
            v = np.concatenate([v, np.zeros(ne - v.size)])
        c_ops.append(np.sqrt(rate) * emb_e(np.diag(v[:ne]).astype(complex)))

    # radiative channels
    rad = [("X", ops["c_X"], gamma_X_ns)]
    if XX is not None:
        rad.append(("XX", ops["c_XX"], gamma_XX_ns))
    if XP is not None:
        rad.append(("Xp", ops["c_Xp"], gamma_charge_ns))
    for name, c, gam in rad:
        photon[name] = (c, gam)

    # cavity
    if cavity:
        a = np.kron(np.kron(np.eye(ne), _destroy(nc)), If)
        ops["a"] = a
        g = float(cavity["g_ns"])
        H += g * (ops["c_X"] @ a.conj().T + ops["c_X"].conj().T @ a)
        H += float(cavity.get("det_ns", 0.0)) * (a.conj().T @ a)
        c_ops.append(np.sqrt(float(cavity["kappa_ns"])) * a)

    sensor_number = np.zeros(dim, dtype=int)
    if filt:
        f = np.kron(np.kron(np.eye(ne), Ic), _destroy(nf))
        ops["f"] = f
        sensor_number = np.tile(np.arange(nf), ne * nc)
        w = float(filt["fwhm_ns"])
        k1 = k2 = 0.5 * w
        src = filt.get("source", "all")
        s = np.zeros((dim, dim), dtype=complex)
        for name, c, gam in rad:
            if gam <= 0:
                continue
            if src == "all" or src == name:
                s = s + np.sqrt(gam) * c
            else:
                c_ops.append(np.sqrt(gam) * c)
        H += 0.5j * np.sqrt(k1) * (s.conj().T @ f - f.conj().T @ s)
        H += float(filt.get("det_ns", 0.0)) * (f.conj().T @ f)
        c_ops.append(s + np.sqrt(k1) * f)
        c_ops.append(np.sqrt(k2) * f)
        ops["f_out"] = np.sqrt(k2) * f
    else:
        for name, c, gam in rad:
            if gam > 0:
                c_ops.append(np.sqrt(gam) * c)

    L = liouvillian(H, c_ops)
    return LindbladSystem(L=L, H=H, c_ops=c_ops, dim=dim, ops=ops, photon=photon,
                          sensor_number=sensor_number,
                          meta=dict(levels=levels, charged=charged,
                                    cavity=cavity, filt=filt))


# ------------------------------------------------------------ solvers

def _scaling_vector(sensor_number, scale):
    """Diagonal similarity T for vec(rho): T[(i,j)] = scale^-(n_i + n_j)."""
    n = np.asarray(sensor_number, dtype=float)
    # column-stack index k = i + j d  (i row, j column)
    return (float(scale) ** -(n[:, None] + n[None, :])).reshape(-1, order="F")


def steady_state(L, dim, *, sensor_number=None, scale=1.0):
    """Stationary rho with Tr rho = 1 (one row of L replaced by the trace row).

    With a filter (sensor_number given and scale < 1) the problem is solved
    in the similarity-scaled basis rho~ = T rho (T = scale^-(n+m)), which
    keeps full relative precision on the ~1e-8 sensor populations."""
    L = np.asarray(L)
    if sensor_number is None or scale == 1.0:
        t = np.ones(dim * dim)
    else:
        t = _scaling_vector(sensor_number, scale)
    Lt = (t[:, None] * L) / t[None, :]
    A = Lt.copy()
    A[0, :] = trace_row(dim) / t
    b = np.zeros(dim * dim, dtype=complex)
    b[0] = 1.0
    v = solve(A, b) / t
    rho = unvec(v, dim)
    rho = 0.5 * (rho + rho.conj().T)
    return rho / np.trace(rho).real


class Propagator:
    """exp(L t) v via one eigendecomposition, falling back to expm when the
    eigenbasis is ill-conditioned (condition number > _EIG_COND_MAX) or does
    not reconstruct L to _EIG_RESID_MAX (e.g. degenerate, purely Hamiltonian
    Liouvillians)."""

    def __init__(self, L):
        self.L = np.asarray(L)
        lam, V = eig(self.L)
        self.lam, self.V = lam, V
        self.ok = False
        try:
            cond = np.linalg.cond(V)
            if np.isfinite(cond) and cond < _EIG_COND_MAX:
                self.Vinv = np.linalg.inv(V)
                nL = np.linalg.norm(self.L)
                resid = np.linalg.norm((V * lam) @ self.Vinv - self.L) / (nL if nL > 0 else 1.0)
                self.ok = bool(resid < _EIG_RESID_MAX)
        except np.linalg.LinAlgError:
            self.ok = False

    def apply(self, v0, taus):
        taus = np.atleast_1d(np.asarray(taus, dtype=float))
        v0 = np.asarray(v0, dtype=complex)
        if self.ok:
            c = self.Vinv @ v0
            out = self.V @ (np.exp(np.outer(self.lam, taus)) * c[:, None])
        else:
            out = np.column_stack([expm(self.L * t) @ v0 for t in taus])
        out[:, taus == 0.0] = v0[:, None]
        return out


def evolve(L, rho0, times, *, method="expm"):
    """rho(t) for each t in `times` (list of dim x dim arrays).
    method="expm" (default; scipy.linalg.expm per time, trace kept to
    ~1e-15 even for stiff L) or "eig" (Propagator; faster for long time
    grids, trace error up to ~1e-11 when L spans >1e5 in rate)."""
    rho0 = np.asarray(rho0, dtype=complex)
    n = rho0.shape[0]
    times = np.atleast_1d(np.asarray(times, dtype=float))
    v0 = vec(rho0)
    if method == "expm":
        return [unvec(expm(np.asarray(L) * t) @ v0, n) for t in times]
    V = Propagator(L).apply(v0, times)
    return [unvec(V[:, k], n) for k in range(times.size)]


def evolve_td(L0, L1, envelope, rho0, t_grid):
    """Time-dependent L(t) = L0 + envelope(t) L1 (e.g. L1 = hamiltonian_super
    of a drive term), piecewise constant on the intervals of t_grid with the
    envelope sampled at each interval midpoint; exact expm per interval.
    Returns the list of rho at every t_grid point (t_grid[0] -> rho0)."""
    rho0 = np.asarray(rho0, dtype=complex)
    n = rho0.shape[0]
    t_grid = np.asarray(t_grid, dtype=float)
    v = vec(rho0)
    out = [unvec(v.copy(), n)]
    for t0, t1 in zip(t_grid[:-1], t_grid[1:]):
        dt = t1 - t0
        f = envelope(0.5 * (t0 + t1))
        v = expm((np.asarray(L0) + f * np.asarray(L1)) * dt) @ v
        out.append(unvec(v.copy(), n))
    return out


# ------------------------------------------------------------ observables

def check_physical(rho):
    """(trace error, min eigenvalue, Hermiticity error) of a density matrix."""
    rho = np.asarray(rho)
    herm = float(np.max(np.abs(rho - rho.conj().T)))
    ev = np.linalg.eigvalsh(0.5 * (rho + rho.conj().T))
    return float(abs(np.trace(rho) - 1.0)), float(ev.min()), herm


def g2_tau(system_or_L, taus, channels=None, weights=None, *, dim=None, rho_ss=None):
    """Stationary g2(tau) of the detected stream (QRT).

    channels : list of (c, gamma) photon channels (default: every radiative
               channel of the system, in the order X, XX, Xp)
    weights  : per-photon detection probabilities (Bernoulli), default 1.
    Returns (g2 array, rho_ss, detected rate I)."""
    if isinstance(system_or_L, LindbladSystem):
        L, dim = system_or_L.L, system_or_L.dim
        if channels is None:
            channels = list(system_or_L.photon.values())
    else:
        L = system_or_L
    if weights is None:
        weights = [1.0] * len(channels)
    if rho_ss is None:
        rho_ss = steady_state(L, dim)
    P = Propagator(L)
    taus = np.abs(np.atleast_1d(np.asarray(taus, dtype=float)))
    I = 0.0
    for (c, gam), w in zip(channels, weights):
        I += w * gam * expect(c.conj().T @ c, rho_ss).real
    G = np.zeros(taus.size)
    for (ci, gi), wi in zip(channels, weights):
        if wi == 0 or gi == 0:
            continue
        vt = P.apply(vec(ci @ rho_ss @ ci.conj().T), taus)
        for (cj, gj), wj in zip(channels, weights):
            if wj == 0 or gj == 0:
                continue
            G += wi * wj * gi * gj * np.real(vec((cj.conj().T @ cj).T) @ vt)
    return G / I ** 2, rho_ss, I


def g2_zero(system_or_L, channels=None, weights=None, *, dim=None, rho_ss=None):
    """g2(0) of the detected stream (QRT at tau = 0)."""
    g2, _, _ = g2_tau(system_or_L, [0.0], channels, weights, dim=dim, rho_ss=rho_ss)
    return float(g2[0])


def _filter_state(system, scale=None, iterations=4):
    f = system.ops["f"]
    if scale is None:
        scale = 1e-3
        for _ in range(iterations):
            rho = steady_state(system.L, system.dim,
                               sensor_number=system.sensor_number, scale=scale)
            n1 = expect(f.conj().T @ f, rho).real
            scale = float(np.clip(np.sqrt(abs(n1)), 1e-12, 1.0))
    return steady_state(system.L, system.dim, sensor_number=system.sensor_number,
                        scale=scale)


def filtered_g2_zero(system, *, scale=None, return_state=False):
    """Zero-delay g2 of the cascaded-filter output, <f^+2 f^2>/<f^+ f>^2."""
    rho = _filter_state(system, scale)
    f = system.ops["f"]
    fd = f.conj().T
    n1 = expect(fd @ f, rho).real
    n2 = expect(fd @ fd @ f @ f, rho).real
    g2 = n2 / n1 ** 2
    return (g2, rho) if return_state else g2


def filtered_flux(system, *, scale=None):
    """Detected output-port photon rate (w/2) <f^+ f> (1/ns)."""
    rho = _filter_state(system, scale)
    fo = system.ops["f_out"]
    return expect(fo.conj().T @ fo, rho).real


def _augmented(L, J):
    N = L.shape[0]
    Z = np.zeros((N, N), dtype=complex)
    return np.block([[L, Z, Z], [J, L, Z], [Z, 2.0 * J, L]])


def counting_jump(system, weights=None, names=None):
    """J rho = sum_i w_i gamma_i c_i rho c_i^+ over the radiative channels."""
    names = list(system.photon) if names is None else names
    weights = [1.0] * len(names) if weights is None else weights
    N = system.dim ** 2
    J = np.zeros((N, N), dtype=complex)
    for nm, w in zip(names, weights):
        c, gam = system.photon[nm]
        J = J + w * gam * jump_super(c)
    return J


def _periodic_fixed_point(Phi, v0):
    v = v0
    for _ in range(_PERIODIC_MAX_ITER):
        vn = Phi @ v
        if np.max(np.abs(vn - v)) < _PERIODIC_TOL:
            return vn, True
        v = vn
    return v, False


def pulsed_counting(L_on, L_off, J, dim, tau_on_ns, tau_dark_ns, *, gate_ns=None,
                    rho_start=None, return_states=False):
    """Per-period factorial moments of the counted stream under a rectangular
    on/off modulation of the Liouvillian (pulse_counting hierarchy, L in place
    of M). Counting runs for 0 <= t < gate_ns from the start of the on window
    (None = whole period). Returns dict(g2, mean_counts, rho_period,
    converged[, states])."""
    period = tau_on_ns + tau_dark_ns
    gate = period if gate_ns is None else min(gate_ns, period)
    N = dim * dim
    Phi = expm(np.asarray(L_off) * tau_dark_ns) @ expm(np.asarray(L_on) * tau_on_ns)
    if rho_start is None:
        rho_start = np.zeros((dim, dim), dtype=complex)
        rho_start[0, 0] = 1.0
    v_ss, conv = _periodic_fixed_point(Phi, vec(rho_start))
    Z = np.zeros((N, N), dtype=complex)
    V = np.zeros(3 * N, dtype=complex)
    V[:N] = v_ss
    states = [unvec(v_ss, dim)]
    windows = []
    on_count = min(gate, tau_on_ns)
    if on_count > 0:
        windows.append((L_on, J, on_count))
    if tau_on_ns - on_count > 0:
        windows.append((L_on, Z, tau_on_ns - on_count))
    dark_count = min(max(gate - tau_on_ns, 0.0), tau_dark_ns)
    if dark_count > 0:
        windows.append((L_off, J, dark_count))
    if tau_dark_ns - dark_count > 0:
        windows.append((L_off, Z, tau_dark_ns - dark_count))
    for Lw, Jw, dt in windows:
        V = expm(_augmented(np.asarray(Lw), Jw) * dt) @ V
        states.append(unvec(V[:N], dim))
    tr = trace_row(dim)
    m1 = (tr @ V[N:2 * N]).real
    m2 = (tr @ V[2 * N:]).real
    out = {"g2": float(m2 / m1 ** 2) if m1 > 1e-300 else float("nan"),
           "mean_counts": float(m1), "rho_period": unvec(V[:N], dim),
           "converged": bool(conv)}
    if return_states:
        out["states"] = states
    return out


def load_kraus(dim_e=3, eta_load=1.0):
    """Kraus operators of the instantaneous one-pair electrical load
    (pulse_counting._deterministic_load_map lifted to a CPTP map [DR]):
        K0 = sqrt(1-eta) I, K1 = sqrt(eta)(|X><G| + |XX><X|), K2 = sqrt(eta)|XX><XX|,
    sum K^+K = I. Emitter-only basis (G, X, XX)."""
    I = np.eye(dim_e, dtype=complex)
    K1 = np.zeros((dim_e, dim_e), dtype=complex)
    K1[1, 0] = 1.0
    K1[2, 1] = 1.0
    K2 = np.zeros((dim_e, dim_e), dtype=complex)
    K2[2, 2] = 1.0
    return [np.sqrt(1.0 - eta_load) * I, np.sqrt(eta_load) * K1, np.sqrt(eta_load) * K2]


def kraus_super(kraus):
    return sum(sprepost(K, K.conj().T) for K in kraus)


def deterministic_cycle_counting(L_off, J, dim, period_ns, kraus, *, gate_ns=None):
    """Factorial moments for a deterministic instantaneous load (Kraus map)
    followed by free evolution L_off for one period, with counting in
    0 <= t < gate_ns (default: whole period). Returns dict(g2, mean_counts,
    rho_before_load, rho_after_load, converged)."""
    N = dim * dim
    K = kraus_super(kraus)
    gate = period_ns if gate_ns is None else min(float(gate_ns), period_ns)
    Phi = expm(np.asarray(L_off) * period_ns) @ K
    rho0 = np.zeros((dim, dim), dtype=complex)
    rho0[0, 0] = 1.0
    v_before, conv = _periodic_fixed_point(Phi, vec(rho0))
    v_after = K @ v_before
    V = np.zeros(3 * N, dtype=complex)
    V[:N] = v_after
    if gate > 0:
        V = expm(_augmented(np.asarray(L_off), J) * gate) @ V
    if period_ns - gate > 0:
        Z = np.zeros((N, N), dtype=complex)
        V = expm(_augmented(np.asarray(L_off), Z) * (period_ns - gate)) @ V
    tr = trace_row(dim)
    m1 = (tr @ V[N:2 * N]).real
    m2 = (tr @ V[2 * N:]).real
    return {"g2": float(m2 / m1 ** 2) if m1 > 1e-300 else float("nan"),
            "mean_counts": float(m1), "rho_before_load": unvec(v_before, dim),
            "rho_after_load": unvec(v_after, dim), "converged": bool(conv)}


def _indistinguishability_eig(L, dim, rho0, A, *, return_counts=False):
    """Eigenmode reference path of ``indistinguishability`` (kept for the
    (d4) cross-check; WRONG at defective L, see the public function).

    HOM indistinguishability of the field A after an instantaneous
    excitation into rho0 (Grange et al., PRL 114, 193601 (2015), Eq. 2):

        I = 2 S / N^2,  S = int dt int dtau |Tr[A^+ e^{L tau}(A rho(t))]|^2,
                        N = int dt Tr[A^+ A rho(t)],

    via analytic Laplace sums over the eigenmodes of L (stationary modes,
    |lambda| < _ZERO_MODE_REL max|lambda|, carry no weight for a decaying
    single excitation). return_counts=True also returns N (photons emitted
    through A per excitation, in units of the A normalisation).

    gamma* convention: for a two-level emitter, I = gamma/(gamma + gamma*)
    (Grange 2015 Eq. 1), with gamma the population decay rate and gamma* the
    rate in D[sqrt(gamma*)|e><e|] (a FWHM contribution: the e-g coherence
    decays at (gamma + gamma*)/2). Equivalently I = Gamma/(Gamma + 2 gamma*_coh)
    with gamma*_coh = gamma*/2 the extra coherence (amplitude) decay rate."""
    L = np.asarray(L)
    lam, V = eig(L)
    Vinv = np.linalg.inv(V)
    Ad = np.asarray(A).conj().T
    c = Vinv @ vec(rho0)
    SB = sprepost(np.asarray(A, dtype=complex), np.eye(dim))   # X -> A X
    W = Vinv @ SB @ V
    a = vec(Ad.T) @ V                      # Tr[A^+ v_j]
    M = (a[:, None] * W) * c[None, :]      # f(t,tau) = sum_jk M_jk e^{l_j tau} e^{l_k t}
    scale = np.max(np.abs(lam))
    small = np.abs(lam) < _ZERO_MODE_REL * scale
    M[small, :] = 0.0
    M[:, small] = 0.0
    den = -(lam[:, None] + lam.conj()[None, :])
    tiny = _ZERO_MODE_REL * scale
    with np.errstate(divide="ignore", invalid="ignore"):
        D = np.where(np.abs(den) > tiny, 1.0 / den, 0.0)
    S = np.einsum("jk,lm,jl,km->", M, M.conj(), D, D)
    b = vec((Ad @ A).T) @ V
    ck = c.copy()
    ck[small] = 0.0
    with np.errstate(divide="ignore", invalid="ignore"):
        inv = np.where(small, 0.0, 1.0 / (-lam))
    N = np.real(np.sum(b * ck * inv))
    I = float(2.0 * np.real(S) / N ** 2)
    return (I, float(N)) if return_counts else I


def _stationary_projector(L):
    """Spectral projector P onto the null space of L along its range.

    Built from the SVD of L (right null vectors R, left null vectors Lt,
    P = R (Lt^+ R)^-1 Lt^+); singular values below _ZERO_MODE_REL x max
    count as zero. L P = P L = 0 and the zero eigenvalue of a Liouvillian is
    semisimple, so P is well defined even when the rest of L is defective."""
    U, sv, Vh = np.linalg.svd(L)
    null = sv < _ZERO_MODE_REL * sv[0]
    R = Vh[null].conj().T
    Ln = U[:, null]
    return R @ np.linalg.solve(Ln.conj().T @ R, Ln.conj().T)


def indistinguishability(L, dim, rho0, A, *, return_counts=False):
    """HOM indistinguishability of the field A after an instantaneous
    excitation into rho0 (Grange et al., PRL 114, 193601 (2015), Eq. 2):

        I = 2 S / N^2,  S = int dt int dtau |Tr[A^+ e^{L tau}(A rho(t))]|^2,
                        N = int dt Tr[A^+ A rho(t)],

    valid for ANY physical L, including defective (non-diagonalizable) ones
    such as a cascaded filter whose rate equals an emitter rate. return_counts
    =True also returns N (photons emitted through A per excitation, in units
    of the A normalisation).

    Method [DR]. Stationary modes are removed explicitly: for a decaying
    single excitation A rho_ss = 0 and Tr[A^+ rho_ss] = 0, so the stationary
    component carries no weight. With P the spectral projector onto null(L)
    (``_stationary_projector``) and Q = 1 - P, the integrand is
        f(t, tau) = a^T e^{L' tau} B e^{L' t} v0,
    with L' = L - P (identical to L on range Q, -1 on null(L), so L' is
    invertible and has the same Jordan structure), v0 = Q vec(rho0),
    B = Q S_A Q (S_A: X -> A X) and a = vec(A^+^T). Then, with the Gramians
        L' Y + Y L'^+ = -v0 v0^+,   L' W + W L'^+ = -B Y B^+,
        S = a^T W a^*   (since |a^T y|^2 = a^T y y^+ a^*),
        N = b^T (-L')^-1 v0,        b = vec((A^+ A)^T).
    This is the doubled-space solve S = (a x a*)^T K^-1 (S_A x S_A*) K^-1
    (v0 x v0*), K = L' x 1 + 1 x L'^*, written as two Sylvester equations
    (Bartels-Stewart on the Schur form) so the n^4 x n^4 matrix is never
    formed; no eigenvector basis is needed. Chosen over a Van Loan block
    exponential because it is exact (no quadrature) and over the doubled
    linear solve for cost; both are algebraically identical.

    gamma* convention: for a two-level emitter, I = gamma/(gamma + gamma*)
    (Grange 2015 Eq. 1), with gamma the population decay rate and gamma* the
    rate in D[sqrt(gamma*)|e><e|] (a FWHM contribution: the e-g coherence
    decays at (gamma + gamma*)/2). Equivalently I = Gamma/(Gamma + 2 gamma*_coh)
    with gamma*_coh = gamma*/2 the extra coherence (amplitude) decay rate."""
    L = np.asarray(L, dtype=complex)
    A = np.asarray(A, dtype=complex)
    Ad = A.conj().T
    n = L.shape[0]
    P = _stationary_projector(L)
    Q = np.eye(n) - P
    Lp = L - P
    v0 = Q @ vec(rho0)
    B = Q @ sprepost(A, np.eye(dim)) @ Q
    a = vec(Ad.T)
    b = vec((Ad @ A).T)
    Y = solve_continuous_lyapunov(Lp, -np.outer(v0, v0.conj()))
    W = solve_continuous_lyapunov(Lp, -(B @ Y @ B.conj().T))
    S = a @ W @ a.conj()
    N = np.real(b @ solve(-Lp, v0))
    I = float(2.0 * np.real(S) / N ** 2)
    return (I, float(N)) if return_counts else I
