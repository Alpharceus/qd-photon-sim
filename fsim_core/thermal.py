"""Module A -- Thermal, analytic first pass (spreading resistance).

Geometry: circular mesa (radius a) on a layer stack over a thick substrate.
  * pillar layers (etched mesa walls, laterally insulated): R = t / (k pi a^2)
  * cone layers (unetched, 45-degree spreading):  R = t / (k pi a (a+t)),
    effective radius grows a -> a + t
  * substrate half-space (isothermal disk):        R = 1 / (4 k a_eff)

k(T) = k300 (300/T)^alpha; the T_j solve iterates k at the mean film
temperature (first-order Kirchhoff correction).

This is the planning-doc "analytic spreading-resistance first pass"; COMSOL
replaces it at Phase-1 exit. The verify suite bounds the cone-approximation
error against an axisymmetric finite-difference solve.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

# k300 (W/m/K), alpha exponent. Alloy/DBR values are effective and [A]-grade:
# override from the card, do not trust silently.
MATERIALS = {
    "GaAs": (55.0, 1.25),
    "Si": (148.0, 1.35),
    "InP": (68.0, 1.4),
    "AlAs": (91.0, 1.2),
    "GaAs_buffer_on_Si": (44.0, 1.25),  # defective heteroepitaxial buffer [E]
}


def k_of_T(k300, alpha, T):
    return k300 * (300.0 / np.asarray(T, dtype=float)) ** alpha


@dataclass
class Layer:
    thickness: float          # m
    k300: float               # W/m/K
    alpha: float = 1.25
    spread: bool = False      # False: pillar (mesa), True: 45-degree cone


@dataclass
class Stack:
    layers: list              # top (junction side) to bottom
    k_sub300: float
    sub_alpha: float = 1.25


def stack_rth(a, stack: Stack, T=300.0):
    """Thermal resistance junction -> heatsink for mesa radius a (m) at
    temperature T (evaluates k(T) of every layer at T)."""
    R, r = 0.0, float(a)
    for L in stack.layers:
        k = k_of_T(L.k300, L.alpha, T)
        if L.spread:
            R += L.thickness / (np.pi * k * r * (r + L.thickness))
            r += L.thickness
        else:
            R += L.thickness / (np.pi * k * r**2)
    R += 1.0 / (4.0 * k_of_T(stack.k_sub300, stack.sub_alpha, T) * r)
    return R


def t_junction(P, a, stack: Stack, T_hs, n_iter=80, T_max=1500.0, relax=0.7):
    """Junction temperature for dissipated power P (W): damped fixed-point
    iteration with k evaluated at the mean film temperature (T_hs + T_j)/2.

    Returns np.inf when the k(T) feedback diverges (T_j would exceed T_max):
    the model predicts thermal runaway at this operating point -- report it,
    do not extrapolate through it."""
    Tj = float(T_hs)
    for _ in range(n_iter):
        Tj_new = T_hs + P * stack_rth(a, stack, T=0.5 * (T_hs + Tj))
        if Tj_new > T_max:
            return np.inf
        if abs(Tj_new - Tj) < 1e-4:
            return Tj_new
        Tj = (1 - relax) * Tj + relax * Tj_new
    return Tj


# ------------------------------------------------------------------ FD second method

def fd_disk_rth(a, k, R_dom=None, Z_dom=None, n=220):
    """Axisymmetric finite-difference solve (constant k): uniform-flux disk of
    radius a on a cylinder (T=0 at bottom and outer wall, insulated top outside
    the disk). Returns R_th = mean disk temperature / P. Second method for the
    spreading closed forms (uniform-flux disk exact: 8/(3 pi^2 k a) ~ 0.270/ka;
    isothermal disk: 0.250/ka)."""
    from scipy.sparse import lil_matrix
    from scipy.sparse.linalg import spsolve

    R_dom = R_dom or 24.0 * a
    Z_dom = Z_dom or 24.0 * a
    nr = nz = n
    dr, dz = R_dom / (nr - 1), Z_dom / (nz - 1)
    rs = np.arange(nr) * dr

    def idx(i, j):
        return i * nz + j

    A = lil_matrix((nr * nz, nr * nz))
    rhs = np.zeros(nr * nz)
    P_tot = 1.0
    q = P_tot / (np.pi * a**2)  # W/m^2 over the disk

    for i in range(nr):
        for j in range(nz):
            m = idx(i, j)
            if i == nr - 1 or j == nz - 1:      # far boundary: heatsink T=0
                A[m, m] = 1.0
                continue
            # conservative FV around node (i,j); face areas ~ 2 pi r dz / dr etc.
            r_e = rs[i] + dr / 2
            r_w = rs[i] - dr / 2 if i > 0 else 0.0
            ae = 2 * np.pi * r_e * dz / dr
            aw = 2 * np.pi * r_w * dz / dr if i > 0 else 0.0
            ring = np.pi * ((rs[i] + dr / 2) ** 2 - (max(rs[i] - dr / 2, 0.0)) ** 2)
            an_ = ring / dz if j > 0 else 0.0
            as_ = ring / dz
            diag = 0.0
            if i > 0:
                A[m, idx(i - 1, j)] = k * aw; diag -= k * aw
            A[m, idx(i + 1, j)] = k * ae; diag -= k * ae
            if j > 0:
                A[m, idx(i, j - 1)] = k * an_; diag -= k * an_
            A[m, idx(i, j + 1)] = k * as_; diag -= k * as_
            A[m, m] = diag
            if j == 0 and rs[i] <= a:           # flux disk at the top surface
                rhs[m] = -q * ring
    T = spsolve(A.tocsr(), rhs)
    disk = [idx(i, 0) for i in range(nr) if rs[i] <= a]
    w = np.array([np.pi * ((rs[i] + dr / 2) ** 2 - (max(rs[i] - dr / 2, 0)) ** 2)
                  for i in range(nr) if rs[i] <= a])
    T_disk = float(np.average(T[disk], weights=w))
    return T_disk / P_tot
