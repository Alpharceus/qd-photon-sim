"""Spectrally filtered and polarisation-gated pulsed photon streams (round 4, numpy/scipy only).

Two opt-in extensions of the pulsed resonant source of ``resonant_source`` (Fable r4 review, pulsed 2.1 / 2.2 and
cavity D). Nothing here is imported by a legacy path: with the filter off (0 GHz) and no fine-structure splitting the
callers never enter this module, so verify_fsim / audit_physics / verify_compose_engine stay bit-identical.
Provenance tags: [V] verified against the cited paper, [DR] derived, [E] estimate, [A] assumption.

Units as in lindblad.py: rates and angular frequencies in 1/ns, times in ns.

1. Spectral filter on the photon stream  (``FilterEmbedding``)
---------------------------------------------------------------
A Lorentzian etalon of FWHM B (Hz) is a filter mode a with decay kappa = 2 pi B (rad/ns), two ports of kappa/2 each
(peak power transmission 1, the convention of ``lindblad.build_system(filt=...)``). The emitter field c_in (the
radiative channel that reaches the detector) drives it through the Gardiner-Carmichael cascade [DR]:
    d rho/dt = L_emit rho + kappa D[a] rho - sqrt(kappa/2) ( [a^+, c_in rho] + [rho c_in^+, a] ),
which equals the joint collapse operator c_in + sqrt(kappa/2) a with the cascade Hamiltonian (i/2) sqrt(kappa/2)
(c_in^+ a - a^+ c_in) plus (kappa/2) D[a]. The transmitted stream is counted with the output port,
J_f rho = (kappa/2) a rho a^+ and A_f = sqrt(kappa/2) a. The Fock space of a is truncated at {0, 1, 2} (two photons can
sit in the filter at once when they are emitted within 1/kappa). Because the cascade is exact, the filtered count
rates, g2 and the Grange HOM overlap come out of the unchanged counting hierarchy of ``resonant_source._stats`` and
``pulsed_indistinguishability`` run on the embedded generator.
Free decay through the filter transmits  T_l = B / (B + (Gamma + gamma*)/(2 pi))  of the photons (area-normalised
emitter Lorentzian against a peak-normalised filter Lorentzian) [DR].
``mode="sensor"`` is the del Valle et al., PRL 109, 183601 (2012) variant: a weakly coupled detector mode,
H = eps (a^+ sigma + sigma^+ a), no back-action; normalised quantities (g2, I) do not depend on eps to O(eps^2).

2. FSS cross-polarised gating  (``VSystem``)
--------------------------------------------
One V-system {G, X, Y} with the fine-structure splitting delta (H = (delta/2)(|X><X| - |Y><Y|)), laser polarisation
|V> = cos(phi)|X> + sin(phi)|Y> (phi to the dot axes), collected polarisation |H> = -sin(phi)|X> + cos(phi)|Y> (cross-
polarised), drive (Omega/2)(|G><V| + h.c.) and photon channels sqrt(Gamma_H)|G><H| (counted) and sqrt(Gamma_V)|G><V|.
Pure dephasing acts on the exciton manifold |X><X| + |Y><Y| and the photonless loss on |G><X|, |G><Y|. A V exciton
prepared by a short pulse emits into H with probability  p_H = (1/2) delta^2 / (Gamma^2 + delta^2)  (Gamma_H = Gamma_V =
Gamma) [DR, Somaschi et al., Nat. Photonics 10, 340 (2016) SI / Fable r4 probe C].
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from . import lindblad as lb
from . import resonant_source as rs

NF = 3                      # filter Fock states {0, 1, 2} [A]
DEFAULT_HOM_TOL = 1e-6      # [A] same as solvers._HOM_TOL
DEFAULT_HOM_DOUBLINGS = 5


def kappa_from_fwhm_GHz(B_GHz: float) -> float:
    """Filter amplitude-decay-rate pair kappa = 2 pi B (rad/ns, total energy decay of the filter mode)."""
    return 2.0 * np.pi * float(B_GHz)


def free_decay_transmission(B_GHz, Gamma_ns, gamma_star_ns=0.0, k_nr_ns=0.0):
    """T_l = B / (B + (Gamma + k_nr + gamma*)/(2 pi)): fraction of a freely decaying Lorentzian photon the etalon passes [DR].
    The line width is the coherence decay of the cascade, (Gamma + k_nr + gamma*)/2 (k_nr = 0 is the pre-r4-review value)."""
    w = (float(Gamma_ns) + float(k_nr_ns) + float(gamma_star_ns)) / (2.0 * np.pi)
    return float(B_GHz) / (float(B_GHz) + w)


def _fock_a(nf):
    return np.diag(np.sqrt(np.arange(1, nf, dtype=float)), 1).astype(complex)


class FilterEmbedding:
    """Emitter (dimension d) x filter-mode (Fock nf) embedding of a Lorentzian filter on the emitter field ``c_in``.

    c_in      (d, d) emitter operator of the radiative channel that feeds the filter (rate included)
    sigma     (d, d) unit-normalised lowering operator of that channel (sensor mode only)
    kappa_ns  filter FWHM in rad/ns (2 pi B)
    mode      "cascade" (exact, default) or "sensor"
    Attributes: D = d nf, L_fc (cascade + filter decay generator), J_out (counting superoperator of the output port),
    A_out (output field operator), lift(S) (emitter superoperator (..., d^2, d^2) -> joint (..., D^2, D^2))."""

    def __init__(self, d, c_in, kappa_ns, mode="cascade", sigma=None, eps_ns=1.0, det_ns=0.0, nf=NF):
        if mode not in ("cascade", "sensor"):
            raise ValueError(f"filter mode must be cascade|sensor, got {mode!r}")
        self.d, self.nf, self.D, self.mode = int(d), int(nf), int(d) * int(nf), mode
        self.kappa = float(kappa_ns)
        D = self.D
        Ie, If = np.eye(self.d, dtype=complex), np.eye(self.nf, dtype=complex)
        a = np.kron(Ie, _fock_a(self.nf))
        ad = a.conj().T
        Lfc = self.kappa * lb.dissipator(a)
        if det_ns:
            Lfc = Lfc + lb.hamiltonian_super(float(det_ns) * (ad @ a))
        if mode == "cascade":
            c = np.kron(np.asarray(c_in, dtype=complex), If)
            kin = 0.5 * self.kappa
            Lfc = Lfc - np.sqrt(kin) * (lb.spre(ad @ c) - lb.sprepost(c, ad) + lb.spost(c.conj().T @ a)
                                        - lb.sprepost(a, c.conj().T))
            kout = 0.5 * self.kappa
        else:
            if sigma is None:
                raise ValueError("sensor mode needs the unit lowering operator sigma")
            s = np.kron(np.asarray(sigma, dtype=complex), If)
            Lfc = Lfc + lb.hamiltonian_super(float(eps_ns) * (ad @ s + s.conj().T @ a))
            kout = self.kappa
        self.L_fc = Lfc
        self.kout = kout
        self.J_out = kout * lb.jump_super(a)
        self.A_out = np.sqrt(kout) * a
        # index map: joint vec index of emitter vec index E and filter pair (i, j); joint basis = kron(emitter, filter)
        d2, nf_ = self.d ** 2, self.nf
        idx = np.empty((d2, nf_, nf_), dtype=int)
        for je in range(self.d):
            for ie in range(self.d):
                for jf in range(nf_):
                    for jf2 in range(nf_):
                        idx[ie + self.d * je, jf2, jf] = (ie * nf_ + jf2) + D * (je * nf_ + jf)
        self._rows = idx[:, None, :, :]
        self._cols = idx[None, :, :, :]

    def lift(self, S):
        """Emitter superoperator S (..., d^2, d^2) acting as S (x) identity on the filter."""
        S = np.asarray(S, dtype=complex)
        out = np.zeros(S.shape[:-2] + (self.D ** 2, self.D ** 2), dtype=complex)
        out[..., self._rows, self._cols] = S[..., :, :, None, None]
        return out

    def pieces(self, pc):
        """``resonant_source._Pieces`` of the filtered output stream built from the emitter pieces ``pc``."""
        return rs._Pieces(L0=self.lift(pc.L0) + self.L_fc, L1=self.lift(pc.L1), J=self.J_out, A=self.A_out,
                          ops={}, dim=self.D, max_dt=pc.max_dt)


# ---------------------------------------------------------------- generic driver on pieces
FINE_STEP_FACTOR = 8      # [A] filtered counting uses 8x the default Gaussian steps: the default sigma/8 grid leaves 1.1e-3 in g2_f


def fine_steps(shape, duration, max_dt=None):
    """Step count of the filtered counting hierarchy: FINE_STEP_FACTOR x the resonant_source default (a rectangle is exact in one step)."""
    return None if shape == "rect" else FINE_STEP_FACTOR * rs._default_steps(shape, duration, None, max_dt)


def pulse_stream(pc, shape, area, duration, f_rep_MHz, hom=True, tol=DEFAULT_HOM_TOL,
                 max_doublings=DEFAULT_HOM_DOUBLINGS, fine=True):
    """(g2, mean photons per period, HOM overlap I, stats dict) of the counted stream of the pieces ``pc``
    (periodic steady state, the resonant_source machinery). fine=True refines the counting grid (see FINE_STEP_FACTOR):
    sigma/8 steps are converged to 1.1e-3 for a filtered stream, 8x finer to 1.5e-5 (convergence run, Ding-like 15 ps pulse)."""
    st = rs._stats(pc, shape, area, duration, f_rep_MHz, fine_steps(shape, duration, pc.max_dt) if fine else None, True)
    I = float("nan")
    if hom:
        I = rs.pulsed_indistinguishability(None, shape, area, duration, tol=tol, max_doublings=max_doublings,
                                           pieces=pc)
    return float(st["g2"]), float(st["mean_photons"]), float(I), st


def filtered_two_level(tl, shape, area, duration, f_rep_MHz, B_GHz, mode="cascade", eps_ns=1.0, hom=True,
                       tol=DEFAULT_HOM_TOL, max_doublings=DEFAULT_HOM_DOUBLINGS):
    """g2, n_ph, I of the pulsed two-level source seen through a Lorentzian filter of FWHM B_GHz.
    Returns dict(g2, n_ph, I, g2_unf, n_unf, I_unf). ``tl`` is a ``resonant_source.TwoLevel``."""
    pc = rs._pieces(tl)
    emb = FilterEmbedding(2, pc.A, kappa_from_fwhm_GHz(B_GHz), mode=mode, eps_ns=eps_ns,
                          sigma=np.array([[0, 1], [0, 0]], dtype=complex))
    g2f, nf, If, _ = pulse_stream(emb.pieces(pc), shape, area, duration, f_rep_MHz, hom, tol, max_doublings)
    st = rs.pulsed_photon_statistics(tl, shape, area, duration, f_rep_MHz)
    out = {"g2": g2f, "n_ph": nf, "I": If, "g2_unf": float(st["g2"]), "n_unf": float(st["mean_photons"])}
    out["I_unf"] = rs.pulsed_indistinguishability(tl, shape, area, duration, tol=tol, max_doublings=max_doublings) \
        if hom else float("nan")
    return out


# ---------------------------------------------------------------- V-system
@dataclass(frozen=True)
class VSystem:
    """FSS V-system {G, X, Y}; rates in 1/ns. fss_ns = delta (rad/ns) = FSS/hbar; phi in rad."""
    Gamma_H_ns: float
    Gamma_V_ns: float
    fss_ns: float
    phi_rad: float
    k_ns: float = 0.0
    gamma_star_ns: float = 0.0
    delta_ns: float = 0.0          # laser detuning of both exciton states, omega_X - omega_L

    def _ket(self, which):
        c, s = np.cos(self.phi_rad), np.sin(self.phi_rad)
        v = np.array([0.0, c, s], dtype=complex) if which == "V" else np.array([0.0, -s, c], dtype=complex)
        return v

    def ops(self):
        G = np.array([1.0, 0.0, 0.0], dtype=complex)
        H, V = self._ket("H"), self._ket("V")
        c_H = np.sqrt(self.Gamma_H_ns) * np.outer(G, H.conj())          # sqrt(Gamma_H) |G><H|
        c_V = np.sqrt(self.Gamma_V_ns) * np.outer(G, V.conj())
        Pexc = np.diag([0.0, 1.0, 1.0]).astype(complex)
        H0 = 0.5 * self.fss_ns * np.diag([0.0, 1.0, -1.0]).astype(complex) + self.delta_ns * Pexc
        H1 = 0.5 * (np.outer(G, V.conj()) + np.outer(V, G.conj()))        # (1/2)(|G><V| + h.c.), Omega = 1
        return dict(c_H=c_H, c_V=c_V, Pexc=Pexc, H0=H0, H1=H1)

    def _L0(self, o):
        cops = [o["c_H"], o["c_V"]]
        if self.k_ns > 0:
            for i in (1, 2):
                m = np.zeros((3, 3), dtype=complex)
                m[0, i] = np.sqrt(self.k_ns)
                cops.append(m)
        if self.gamma_star_ns > 0:
            cops.append(np.sqrt(self.gamma_star_ns) * o["Pexc"])
        return lb.liouvillian(o["H0"], cops)

    def pieces(self, channel="H"):
        """Pieces counting the H photons ('H'), the V photons ('V') or both ('all')."""
        o = self.ops()
        L0, L1 = self._L0(o), lb.hamiltonian_super(o["H1"])
        J = {"H": lb.jump_super(o["c_H"]), "V": lb.jump_super(o["c_V"]),
             "all": lb.jump_super(o["c_H"]) + lb.jump_super(o["c_V"])}[channel]
        A = o["c_H"] if channel == "H" else (o["c_V"] if channel == "V" else None)
        return rs._Pieces(L0=L0, L1=L1, J=J, A=A, ops=o, dim=3)

    def filtered_pieces(self, B_GHz, mode="cascade", eps_ns=1.0):
        pc = self.pieces("H")
        sig = np.outer(np.array([1.0, 0, 0], dtype=complex), self._ket("H").conj())
        # r4 review fix 2: dot-centred filter, det_ns as in the two-level branch (filterhook: det_ns=tl.delta_ns)
        emb = FilterEmbedding(3, pc.A, kappa_from_fwhm_GHz(B_GHz), mode=mode, sigma=sig, eps_ns=eps_ns, det_ns=self.delta_ns)
        return emb.pieces(pc), emb


def h_yield_short_pulse(Gamma_ns, fss_ns, phi_rad=np.pi / 4):
    """Closed form of the H-photon yield of a V exciton prepared instantaneously (phi = 45 deg, Gamma_H = Gamma_V):
    (1/2) delta^2/(Gamma^2 + delta^2) [DR]."""
    if abs(phi_rad - np.pi / 4) > 1e-12:
        raise ValueError("closed form is for phi = 45 deg")
    return 0.5 * fss_ns ** 2 / (Gamma_ns ** 2 + fss_ns ** 2)
