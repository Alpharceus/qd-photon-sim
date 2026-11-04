"""Compose-layer glue of the round-4 spectral filter and FSS polarisation gating (Fable r4 pulsed 2.1 / 2.2, cavity D).

Physics lives in ``fsim_core.filtered_stream``; this module only
  * turns the opt-in leaves (setup.filter_fwhm_GHz / filter_applies_to, material.fss_ueV, drive.pol_angle_deg,
    structure.mode_split_ueV / F_P_ratio_VH) into plan.inputs entries (``attach``, called by engine._prepare),
  * post-processes the lindblad two-level result (``lindblad_postprocess``) and the polaron step-propagator result
    (``polaron_filter``).
Every function returns its input unchanged when the leaves are at their defaults (filter 0 GHz, FSS 0), so the default
numbers are bit-identical. Leaves are read through the tracking view only when they are active, so the provenance of
an unfiltered run is unchanged.
"""
from __future__ import annotations

import math
from dataclasses import replace

import numpy as np
from scipy.linalg import expm

from .. import filtered_stream as fs
from .. import lindblad as lb
from .. import resonant_source as rs

NAN = float("nan")
APPLIES = ("g2_hom", "g2_hom_brightness", "hom")


def _brightness(plan) -> bool:
    return str(plan.inputs.get("filter_applies_to", "g2_hom")) == "g2_hom_brightness"


def _g2_filtered(plan) -> bool:
    """r4 red-team item 2: ``hom`` = the etalon sits only in the HOM path (g2 and brightness unfiltered, the HOM overlap
    is that of the filtered stream)."""
    return str(plan.inputs.get("filter_applies_to", "g2_hom")) != "hom"


# ---------------------------------------------------------------- plan inputs
def mode_ratio_H(rt, split_ueV: float, ratio_VH: float, resonant_mode: str = "V") -> float:
    """F_P(H)/F_P(V-reference) of the polarisation-split cavity, where the V-reference is the mode with the full F_P
    carried by ``rt`` (the mode the dot is resonant with).
    resonant_mode "V" (default): L(split) / (F_P ratio V/H), L = 1/(1 + (2 split/kappa)^2) (Lorentzian cavity mode,
    narrow dot line) [DR, Fable r4 cavity C/D]. 1 = degenerate modes.
    resonant_mode "H" [A]: the dot is resonant with the H (collected) mode, which gets the full F_P: returns 1.0 (the V
    mode then gets L(split) x ratio, see ``mode_ratio_V``)."""
    if str(resonant_mode) == "H":
        return 1.0
    L = 1.0
    if split_ueV > 0.0 and rt.kappa_ns is not None and rt.kappa_ns > 0.0:
        kap_ueV = float(rt.kappa_ns) * lb.HBAR_MEV_NS * 1e3
        L = 1.0 / (1.0 + (2.0 * float(split_ueV) / kap_ueV) ** 2)
    return L / float(ratio_VH)


def mode_ratio_V(rt, split_ueV: float, ratio_VH: float, resonant_mode: str = "V") -> float:
    """F_P(V)/F_P(resonant mode): 1.0 when the dot is resonant with V (default); L(split) x ratio_VH when it is
    resonant with H (the V mode is then the detuned one) [A]."""
    if str(resonant_mode) != "H":
        return 1.0
    L = 1.0
    if split_ueV > 0.0 and rt.kappa_ns is not None and rt.kappa_ns > 0.0:
        kap_ueV = float(rt.kappa_ns) * lb.HBAR_MEV_NS * 1e3
        L = 1.0 / (1.0 + (2.0 * float(split_ueV) / kap_ueV) ** 2)
    return L * float(ratio_VH)


def attach(plan, em, sr, rt, vals, tr):
    """Add the opt-in filter / V-system entries to a coherent pulsed plan (``vals`` = untracked values per axis,
    ``tr(axis)`` = tracking view). Returns the plan unchanged when nothing is active."""
    if plan.family != "coherent" or plan.timing != "pulsed" or not plan.pulses:
        return plan
    add = {}
    B = float(vals.get("setup", {}).get("filter_fwhm_GHz", 0.0) or 0.0)
    if B > 0.0:
        S = tr("setup")
        add["filter_fwhm_GHz"] = float(S["filter_fwhm_GHz"])
        add["filter_applies_to"] = str(S["filter_applies_to"])
    fss = float(getattr(em, "fss_ueV", 0.0) or 0.0)
    if fss > 0.0 and plan.inputs.get("cross_pol"):
        tr("material")["fss_ueV"]
        phi = float(tr("drive")["pol_angle_deg"])
        st = vals.get("structure", {})
        split = float(st.get("mode_split_ueV", 0.0) or 0.0)
        ratio = float(st.get("F_P_ratio_VH", 1.0) or 1.0)
        if split > 0.0:
            tr("structure")["mode_split_ueV"]
        if ratio != 1.0:
            tr("structure")["F_P_ratio_VH"]
        res_mode = str(st.get("resonant_mode", "V") or "V")
        if res_mode != "V":
            tr("structure")["resonant_mode"]
        rH = mode_ratio_H(rt, split, ratio, res_mode)
        rV = mode_ratio_V(rt, split, ratio, res_mode)
        G_bg = float(rt.Gamma_bg_ns)
        Gam_res = float(rt.Gamma_ns)
        Gam_V = Gam_res if rV == 1.0 else G_bg + (Gam_res - G_bg) * rV
        Gam_H = G_bg + (Gam_res - G_bg) * rH
        G0 = float(rt.Gamma0_ns)
        F_bg = G_bg / G0 if G0 > 0.0 else 0.0
        FPL_H = float(rt.F_P_eff) * rH
        beta_H = (float(rt.beta0) * F_bg + FPL_H) / (F_bg + FPL_H) if (F_bg + FPL_H) > 0.0 else float(rt.beta)
        add["vsys"] = {"fss_ueV": fss, "pol_angle_deg": phi, "Gamma_V_ns": Gam_V, "Gamma_H_ns": Gam_H,
                       "mode_ratio_H": rH, "beta_H": beta_H, "mode_split_ueV": split, "F_P_ratio_VH": ratio}
        if res_mode != "V":
            add["vsys"].update(resonant_mode=res_mode, mode_ratio_V=rV)
    if not add:
        return plan
    return replace(plan, inputs={**plan.inputs, **add})


def _vsystem(plan, rt):
    v = plan.inputs["vsys"]
    return fs.VSystem(Gamma_H_ns=float(v["Gamma_H_ns"]), Gamma_V_ns=float(v["Gamma_V_ns"]),
                      fss_ns=float(lb.mev_to_rate(float(v["fss_ueV"]) * 1e-3)),
                      phi_rad=math.radians(float(v["pol_angle_deg"])), k_ns=float(rt.k_nr_ns),
                      gamma_star_ns=float(rt.gamma_star_ns), delta_ns=float(plan.pulses[0].detuning_meV) *
                      float(lb.mev_to_rate(1.0)))


def _filter_note(B, kind):
    return f"{kind}: Lorentzian etalon {B:g} GHz FWHM on the detected stream (exact Gardiner-Carmichael cascade, Fock {{0,1,2}}) [DR]"


# ---------------------------------------------------------------- lindblad (two-level / V-system)
def lindblad_postprocess(res, plan, em, rt, window_counts, window_p0):
    """Filtered / FSS-gated replacement of the unfiltered two-level pulsed result ``res`` (lindblad solver)."""
    B = float(plan.inputs.get("filter_fwhm_GHz", 0.0) or 0.0)
    vs = plan.inputs.get("vsys")
    if not (B > 0.0 or vs):
        return res
    p = plan.pulses[0]
    dur = p.fwhm_ps * 1e-3
    f_rep = float(plan.f_rep_MHz)
    extra = dict(res.extra)
    notes = list(extra.get("notes", []))
    reasons = list(res.reasons)
    if vs:
        V = _vsystem(plan, rt)
        pcH = V.pieces("H")
        pc_all = V.pieces("all")
        g2u, nu, Iu, stu = fs.pulse_stream(pcH, p.shape, p.area_rad, dur, f_rep)
        rho_w, emitted = window_counts(pc_all, p.shape, p.area_rad, dur)
        p0 = window_p0(pc_all, p.shape, p.area_rad, dur, (0,))
        rho_post = stu["rho_post_pulse"]
        extra.update(eta_prep=1.0 - p0, n_exc=float((rho_w[1, 1] + rho_w[2, 2]).real) + emitted, photons_in_window=emitted,
                     P_X_end=float((rho_post[1, 1] + rho_post[2, 2]).real), beta=float(vs["beta_H"]),
                     p_H_collected=float(nu), vsys=dict(vs))
        g2, n, I, inv = g2u, nu, Iu, float((rho_post[1, 1] + rho_post[2, 2]).real)
        notes.append(f"FSS gating: V-system {{G,X,Y}}, FSS {vs['fss_ueV']:g} ueV, drive polarisation {vs['pol_angle_deg']:g} deg to the "
                     f"dot axes, H collection; Gamma_V {vs['Gamma_V_ns']:.4g}/ns, Gamma_H {vs['Gamma_H_ns']:.4g}/ns "
                     f"(H-mode / V-mode F_P ratio {vs['mode_ratio_H']:.3g}); n_ph is the H-photon number per pulse [DR]")
        if B > 0.0:
            pf, emb = V.filtered_pieces(B)
            g2f_, nf, I, _ = fs.pulse_stream(pf, p.shape, p.area_rad, dur, f_rep, max_doublings=4)
            g2 = g2f_ if _g2_filtered(plan) else g2u
            extra.update(g2_unfiltered=g2u, g2_filtered=float(g2f_), I_unfiltered=Iu, n_ph_unfiltered=nu)
            n = nf if _brightness(plan) else nu
            notes.append(_filter_note(B, "filter"))
    else:
        tl = rs.TwoLevel(float(rt.Gamma_ns), k_ns=float(rt.k_nr_ns), gamma_star_ns=float(rt.gamma_star_ns),
                         delta_ns=float(lb.mev_to_rate(float(p.detuning_meV))))
        pc = rs._pieces(tl)
        emb = fs.FilterEmbedding(2, pc.A, fs.kappa_from_fwhm_GHz(B), det_ns=tl.delta_ns)
        g2f_, nf, I, _ = fs.pulse_stream(emb.pieces(pc), p.shape, p.area_rad, dur, f_rep)
        g2 = g2f_ if _g2_filtered(plan) else float(res.g2_0)
        extra.update(g2_unfiltered=float(res.g2_0), g2_filtered=float(g2f_), I_unfiltered=float(res.I_single),
                     n_ph_unfiltered=float(res.n_ph))
        n = nf if _brightness(plan) else float(res.n_ph)
        inv = res.inversion
        notes.append(_filter_note(B, "filter"))
    if B > 0.0:
        extra.update(filter_fwhm_GHz=B, n_ph_filtered=float(n) if _brightness(plan) else NAN,
                     filter_T_free_decay=fs.free_decay_transmission(B, float(rt.Gamma_ns), float(rt.gamma_star_ns)),
                     filter_applied=True)
    extra["notes"] = notes
    return replace(res, g2_0=float(g2), n_ph=float(n), I_single=float(I), inversion=inv, extra=extra,
                   reasons=tuple(reasons))


# ---------------------------------------------------------------- polaron step propagators
class _FilteredGen:
    """Duck-typed ``polaron.Generator`` of the emitter x filter product space (same grid, lifted step propagators)."""

    def __init__(self, gen, emb):
        self.gen, self.emb = gen, emb
        self.dim = emb.D
        self.L0 = emb.lift(gen.L0) + emb.L_fc
        self.J, self.A, self.W = emb.J_out, emb.A_out, gen.W
        self._cache = {}

    def grid(self, n):
        if n not in self._cache:
            edges, L, _ = self.gen.grid(n)
            Lf = self.emb.lift(L) + self.emb.L_fc
            dt = np.diff(edges)
            U = np.array([expm(Lf[i] * dt[i]) for i in range(n)])
            self._cache[n] = (edges, Lf, U)
        return self._cache[n]


def _counting_f(fg, n_f, f_rep_MHz):
    edges, L, U = fg.grid(n_f)
    dt = np.diff(edges)
    D = fg.dim
    N = D * D
    period = 1.0e3 / float(f_rep_MHz)
    dark = period - fg.W
    P = np.eye(N, dtype=complex)
    for u in U:
        P = u @ P
    v0 = rs._fixed_point(expm(fg.L0 * dark) @ P, D)
    V = np.zeros(3 * N, dtype=complex)
    V[:N] = v0
    for i in range(n_f):
        V = expm(rs._aug(L[i], fg.J) * dt[i]) @ V
    V = expm(rs._aug(fg.L0, fg.J) * dark) @ V
    tr = lb.trace_row(D)
    m1, m2 = float((tr @ V[N:2 * N]).real), float((tr @ V[2 * N:]).real)
    return (m2 / m1 ** 2 if m1 > 1e-300 else NAN), m1


def _group_f(U, m):
    if m == 1:
        return U
    n = U.shape[0] // m
    out = np.empty((n,) + U.shape[1:], dtype=complex)
    for j in range(n):
        P = U[j * m]
        for k in range(1, m):
            P = U[j * m + k] @ P
        out[j] = P
    return out


def _hom_pulse_f(fg, Uq, h):
    D = fg.dim
    A = fg.A
    SA = lb.sprepost(A, np.eye(D))
    a = lb.vec(A.conj())
    b = lb.vec((A.conj().T @ A).T)
    g0 = np.zeros((D, D), dtype=complex)
    g0[0, 0] = 1.0
    v = lb.vec(g0)
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
    return F_aa, X, N_pulse, lb.unvec(R[:, -1], D), a


def _hom_f(fg, nq0, fine_for, tol, max_doublings, node_cap):
    L = fg.L0
    D = fg.dim
    P = lb._stationary_projector(L)
    Q = np.eye(L.shape[0]) - P
    Lp = L - P

    def evaluate(nq):
        n_f = fine_for(nq)
        _, _, U = fg.grid(n_f)
        Uq = _group_f(U, n_f // nq)
        F_aa, X, N_pulse, rho_end, a = _hom_pulse_f(fg, Uq, fg.W / nq)
        Y = rs._lyap_pair(Lp, Q @ X @ Q.conj().T)
        F_ab = float(np.real(a @ Y @ a.conj()))
        I_c, N_c = lb.indistinguishability(L, D, rho_end, fg.A, return_counts=True)
        Ntot = N_pulse + N_c
        if not Ntot > 1e-300:
            return NAN, Ntot
        return (F_aa + 2.0 * F_ab + I_c * N_c ** 2) / Ntot ** 2, Ntot

    nq = nq0
    prev, Nprev = evaluate(nq)
    if not math.isfinite(prev):
        return NAN
    out = prev
    for _ in range(max_doublings):
        if 2 * nq > node_cap:
            break
        nq *= 2
        cur, Ncur = evaluate(nq)
        err = abs(cur - prev) / 3.0
        out = cur + (cur - prev) / 3.0
        if err < tol * max(abs(out), 1e-300):
            break
        prev, Nprev = cur, Ncur
    return float(out)


def polaron_filter(gen, plan, st, I, nq0, fine_for, n_f, with_hom, hom_tol, hom_doublings, node_cap):
    """Filtered g2 / n_ph / I for the polaron step propagators ``gen`` (phonon rates and chirp already inside gen.grid).
    Returns (st', I', extra) with st' a copy of the unfiltered counting dict with g2 (and n_ph when the filter applies to
    brightness) replaced."""
    B = float(plan.inputs["filter_fwhm_GHz"])
    emb = fs.FilterEmbedding(2, gen.A, fs.kappa_from_fwhm_GHz(B), det_ns=float(gen.delta))
    fg = _FilteredGen(gen, emb)
    g2f, nf = _counting_f(fg, n_f, plan.f_rep_MHz)
    If = _hom_f(fg, nq0, fine_for, hom_tol, hom_doublings, node_cap) if with_hom else NAN
    st2 = dict(st)
    if _g2_filtered(plan):
        st2["g2"] = g2f
    if _brightness(plan):
        st2["mean_photons"] = nf
    extra = {"filter_fwhm_GHz": B, "filter_applied": True, "g2_unfiltered": float(st["g2"]),
             "I_unfiltered": float(I), "n_ph_unfiltered": float(st["mean_photons"]), "g2_filtered": float(g2f),
             "n_ph_filtered": float(nf) if _brightness(plan) else NAN}
    return st2, If, extra
