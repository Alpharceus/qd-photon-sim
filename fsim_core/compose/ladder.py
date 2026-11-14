"""Four-level exciton-biexciton ladder solver (devlib R3-LADDER; Piccinini et al. 2025, arXiv:2510.21428v2, Methods Eqs 1-4).

Phonon-free Lindblad model on the basis {G, XH, XV, XX}: H-polarised laser (Piccinini Methods [A]), so |G> <-> |XH> <-> |XX>
are driven and |XV> is reached only by the XX -> XV emission. In the frame of laser 1 (the engine's PulseTrain frame, repo
det_X = omega_X - omega_L1; the Rabi frequency Omega(t) already carries the beat of the other pulses):

    H/hbar = e1 |XH><XH| + e2 |XX><XX| + (Omega(t)/2) (|XH><G| + |XX><XH|) + h.c.,
    e1 = det_X + FSS/2,   e2 = 2 det_X + delta_xx   (delta_xx = -E_b, the engine's FockSource sign convention),
    Lindblad  sqrt(Gamma){|G><XH|, |G><XV|}, sqrt(r_xx Gamma){|XH><XX|, |XV><XX|}   (r_xx = Gamma_XX/Gamma_X = 1 [V]).

|XV> carries no coherence (nothing drives it), so the density matrix is propagated as the 3 x 3 block {G, XH, XX} plus the scalar
rho_VV: eight complex numbers per area pair, batched over all pairs; RK4 with the field sampled exactly at the half steps.

Observables (X line only, the line the HBT setup of the paper detects) [DR]:
  N_X    = Int Gamma (rho_HH + rho_VV) dt + (rho_HH + rho_VV + rho_XX)(end): every excitation left at the window end
           ends as exactly one X-line photon (XH -> G, XV -> G, and XX -> X -> G gives one X-line photon);
  g2(0)  = m2 / N_X^2, m2 = 2 Int dt1 Gamma P_X(t1) N_X,after(t1 | restart in |G><G|) (quantum regression of the factorial
           moment <n(n-1)>; the peak-area g2 of a pulse train whose pulses do not overlap). N_X,after(t | G) is the G-G element
           of the Heisenberg-picture counting operator Lambda(t), integrated BACKWARD once (dLambda/dtau = L^dagger Lambda + E,
           tau = T - t, Lambda(0) = projector on the excited states, E = Gamma (|H><H| + |V><V|)): all restart times at once, no
           quadrature over restarted trajectories (a restart-trajectory quadrature needs a spacing << the 1 ps Rabi period).
The 80 MHz period is ~60 lifetimes, so the pulses do not overlap and one pulse is enough.
"""
from __future__ import annotations

import math
from copy import copy
from dataclasses import dataclass, replace
from typing import Optional

import numpy as np
from scipy.integrate import simpson

from .. import lindblad as lb
from .polaron import PulseTrain, comp_norm

DT_PS = 0.02                 # RK4 step [A: numerical]; convergence against 0.005 ps is checked in verify_compose_ladder
CHUNK = 2048                 # area pairs per RK4 batch (cache size) [A]
# state slots of Y (8, B): populations G, XH, XX, XV, cumulative X-line emission N, coherences GH, GXX, HXX
_G, _H, _B, _V, _N, _GH, _GB, _HB = range(8)


@dataclass(frozen=True)
class LadderSetup:
    train: PulseTrain
    n: int                       # number of RK4 steps
    dt: float                    # ns
    series: tuple                # per component: unit-amplitude complex Rabi series at the 2n+1 nodes
    norms: tuple                 # per component: Int env dt for peak amplitude 1 (ns)
    e1: float                    # rad/ns
    e2: float
    gamma: float                 # 1/ns, Gamma_X
    gamma_xx: float              # 1/ns, Gamma_XX (per channel)
    cxx: float                   # 1 = XH <-> XX driven, 0 = XX coupling off (two-level control)
    info: dict


def prepare(plan, em, rt, *, fss_ueV: float = 0.0, xx_coupling: bool = True, dt_ps: float = DT_PS,
            gamma_ns: Optional[float] = None) -> LadderSetup:
    """Build the ladder for a coherent pulsed plan (polaron PulseTrain pulses, any registered shape)."""
    if getattr(rt, "gamma_star_ns", 0.0) or getattr(rt, "k_nr_ns", 0.0):
        raise ValueError("ladder: the paper model has no pure dephasing and no photonless loss (gamma* = k_nr = 0 required)")
    train = PulseTrain(plan.pulses, plan.det_profile, plan.inputs.get("slm_field"))
    if train.field:
        raise ValueError("ladder: the pulse-shaper field (slm_mask) is not supported")
    if plan.det_profile is not None:
        raise ValueError("ladder: a time-dependent detuning (Stark gate, det_profile) is not supported")
    n = max(1, int(math.ceil(train.W / (dt_ps * 1e-3))))
    dt = train.W / n
    tn = np.arange(2 * n + 1) * (0.5 * dt)
    series, norms = [], []
    for j, c in enumerate(train.comps):
        t2 = copy(train)
        t2.comps = tuple(replace(cc, amp=1.0 if k == j else 0.0) for k, cc in enumerate(train.comps))
        series.append(t2.omega(tn))
        norms.append(comp_norm(c))
    det = float(train.det0_ns)
    gam = float(rt.Gamma_ns if gamma_ns is None else gamma_ns)
    return LadderSetup(
        train=train, n=n, dt=dt, series=tuple(series), norms=tuple(norms),
        e1=det + float(lb.mev_to_rate(0.5 * fss_ueV * 1e-3)), e2=2.0 * det + float(lb.mev_to_rate(em.delta_xx_meV)),
        gamma=gam, gamma_xx=gam * float(em.gamma_xx_ratio), cxx=1.0 if xx_coupling else 0.0,
        info={"dt_ps": dt * 1e3, "n_steps": n, "window_ps": train.W * 1e3, "fss_ueV": float(fss_ueV),
              "E_b_meV": -float(em.delta_xx_meV), "Gamma_ns": gam, "Gamma_XX_ns": gam * float(em.gamma_xx_ratio),
              "xx_coupling": bool(xx_coupling)})


def _comm(Y, om, e1, e2, cxx):
    """U = -i[H, Y] on the 3 x 3 block {G, XH, XX} of Y (8, B); the XV / counter slots stay 0. om (B,) complex (rad/ns)."""
    r00, r11, r22 = Y[_G].real, Y[_H].real, Y[_B].real
    r01, r02, r12 = Y[_GH], Y[_GB], Y[_HB]
    c01, c02, c12 = np.conj(r01), np.conj(r02), np.conj(r12)
    oh = 0.5 * om
    h10, h01 = oh, np.conj(oh)
    h21 = cxx * oh
    h12 = np.conj(h21)
    # M = H Y (H tridiagonal: H01 = conj(Om)/2, H10 = Om/2, H12 = conj(H21), H21, diagonal 0, e1, e2); -i[H,Y] = -i (M - M^dagger)
    m00 = h01 * c01
    m01 = h01 * r11
    m10 = h10 * r00 + e1 * c01 + h12 * c02
    m11 = h10 * r01 + e1 * r11 + h12 * c12
    m02 = h01 * r12
    m12 = h10 * r02 + e1 * r12 + h12 * r22
    m20 = h21 * c01 + e2 * c02
    m21 = h21 * r11 + e2 * c12
    m22 = h21 * r12 + e2 * r22
    U = np.zeros_like(Y)
    U[_G] = 2.0 * m00.imag
    U[_H] = 2.0 * m11.imag
    U[_B] = 2.0 * m22.imag
    U[_GH] = -1j * (m01 - np.conj(m10))
    U[_GB] = -1j * (m02 - np.conj(m20))
    U[_HB] = -1j * (m12 - np.conj(m21))
    return U


def _rhs(Y, om, e1, e2, gam, gxx, cxx):
    """Forward Lindblad dY/dt for Y (8, B) complex (slots: rho_GG, rho_HH, rho_XX, rho_VV, N emitted, rho_GH, rho_GXX, rho_HXX)."""
    D = _comm(Y, om, e1, e2, cxx)
    r11, r22, rvv = Y[_H].real, Y[_B].real, Y[_V].real
    D[_G] += gam * (r11 + rvv)
    D[_H] += gxx * r22 - gam * r11
    D[_B] += -2.0 * gxx * r22
    D[_V] += gxx * r22 - gam * rvv
    D[_N] += gam * (r11 + rvv)
    D[_GH] -= 0.5 * gam * Y[_GH]
    D[_GB] -= gxx * Y[_GB]
    D[_HB] -= (0.5 * gam + gxx) * Y[_HB]
    return D


def _rhs_adj(L, om, e1, e2, gam, gxx, cxx):
    """Backward (tau = T - t) Heisenberg-picture counting operator, same slots (slot _V = Lambda_VV): dL/dtau = L^dagger(L) + E."""
    D = -_comm(L, om, e1, e2, cxx)
    l00, l11, l22, lvv = L[_G].real, L[_H].real, L[_B].real, L[_V].real
    D[_H] += gam * (l00 - l11) + gam
    D[_B] += gxx * (l11 + lvv - 2.0 * l22)
    D[_V] += gam * (l00 - lvv) + gam
    D[_GH] -= 0.5 * gam * L[_GH]
    D[_GB] -= gxx * L[_GB]
    D[_HB] -= (0.5 * gam + gxx) * L[_HB]
    return D


def _step(Y, o0, o1, o2, dt, p, rhs=_rhs):
    k1 = rhs(Y, o0, *p)
    k2 = rhs(Y + (0.5 * dt) * k1, o1, *p)
    k3 = rhs(Y + (0.5 * dt) * k2, o1, *p)
    k4 = rhs(Y + dt * k3, o2, *p)
    return Y + (dt / 6.0) * (k1 + 2.0 * k2 + 2.0 * k3 + k4)


def _omega_nodes(S: LadderSetup, areas_rad) -> np.ndarray:
    """(2n+1, B) complex Rabi frequency at the half-step nodes for the per-component areas (each (B,) rad)."""
    om = 0.0
    for k, ser in enumerate(S.series):
        om = om + (np.asarray(areas_rad[k], dtype=float)[None, :] / S.norms[k]) * ser[:, None]
    return om


def _ground(B):
    Y = np.zeros((8, B), dtype=complex)
    Y[_G] = 1.0
    return Y


def scan(S: LadderSetup, areas_pi: list) -> dict:
    """Propagate |G><G| through the pulse train for every area tuple. `areas_pi`: one 1-D array (B,) per component.
    -> {N_X, pX_end, pXX_end} (arrays (B,))."""
    B = len(areas_pi[0])
    out = {k: np.empty(B) for k in ("N_X", "pX_end", "pXX_end")}
    p = (S.e1, S.e2, S.gamma, S.gamma_xx, S.cxx)
    for c0 in range(0, B, CHUNK):
        sl = slice(c0, min(B, c0 + CHUNK))
        om = _omega_nodes(S, [np.asarray(a)[sl] * math.pi for a in areas_pi])
        Y = _ground(om.shape[1])
        for i in range(S.n):
            Y = _step(Y, om[2 * i], om[2 * i + 1], om[2 * i + 2], S.dt, p)
        pX = Y[_H].real + Y[_V].real
        out["N_X"][sl] = Y[_N].real + pX + Y[_B].real
        out["pX_end"][sl] = pX
        out["pXX_end"][sl] = Y[_B].real
    return out


def g2_at(S: LadderSetup, areas_pi: list) -> dict:
    """X-line N_X, g2(0), n_ph and the end-of-window populations at ONE area tuple (see the module docstring)."""
    n = S.n
    p = (S.e1, S.e2, S.gamma, S.gamma_xx, S.cxx)
    om = _omega_nodes(S, [np.array([a * math.pi]) for a in areas_pi])
    Y = _ground(1)
    pxs = np.empty(n + 1)
    pxs[0] = 0.0
    for i in range(n):
        Y = _step(Y, om[2 * i], om[2 * i + 1], om[2 * i + 2], S.dt, p)
        pxs[i + 1] = Y[_H, 0].real + Y[_V, 0].real
    px_end, pxx_end = float(pxs[n]), float(Y[_B, 0].real)
    m1 = float(Y[_N, 0].real) + px_end + pxx_end
    # backward pass: Lambda_GG(t_i) = expected X-line photons after t_i for a dot restarted in |G><G| at t_i
    Lam = np.zeros((8, 1), dtype=complex)
    Lam[_H] = Lam[_B] = Lam[_V] = 1.0
    nafter = np.empty(n + 1)
    nafter[n] = 0.0
    for j in range(n):
        m = 2 * (n - j)
        Lam = _step(Lam, om[m], om[m - 1], om[m - 2], S.dt, p, rhs=_rhs_adj)
        nafter[n - j - 1] = Lam[_G, 0].real
    m2 = 2.0 * S.gamma * float(simpson(pxs * nafter, dx=S.dt))
    return {"N_X": m1, "m2": m2, "g2_0": m2 / m1 ** 2, "n_ph": m1, "pX_end": px_end, "pXX_end": pxx_end}
