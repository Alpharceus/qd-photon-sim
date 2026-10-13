"""Checks for the composable engine (rates, solvers, observables, engine; P5).

Fixture entries (kind prefix "test_") are registered here so the checks do not depend on the real
materials / structures / drives; a final block also runs real registry entries end to end.

Run: python verify/verify_compose_engine.py   (exit 0 iff all pass)
"""
from __future__ import annotations

import json
import math
import sys
from dataclasses import dataclass
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core import lindblad as lb                                   # noqa: E402
from fsim_core import qd_gf                                            # noqa: E402
from fsim_core import resonant_source as rs                            # noqa: E402
from fsim_core.compose import engine, observables as OB, rates, solvers  # noqa: E402
from fsim_core.compose import registry as R                            # noqa: E402
from fsim_core.compose.params import ParamSpec                         # noqa: E402
from fsim_core.compose.types import DrivePlan, EmitterModel, Pulse, StructureResponse  # noqa: E402
from fsim_studio.serialize import dumps                                # noqa: E402

RESULTS = []
HB = lb.HBAR_MEV_NS


def check(name, ok, detail=""):
    RESULTS.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))


def near(a, b, rel=0.0, ab=0.0):
    return abs(a - b) <= max(rel * abs(b), ab)


# ------------------------------------------------------------------ fixtures
def P(name, unit, default, tag="V", **kw):
    return ParamSpec(name, unit, default, tag=tag, source="test fixture", **kw)


@dataclass(frozen=True)
class TMat(R.Material):
    def build(self, vals, T_K):
        alpha = vals["alpha_ps2"]
        return EmitterModel(
            E_X_meV=1380.0, lambda_nm=vals["lambda_nm"], n_host=3.5, gamma0_ns=vals["gamma0_ns"],
            k_nr_ns=vals["k_nr_ns"], gamma_star_ns=vals["gamma_star_ns"],
            phonon=qd_gf.PhononParams(alpha_ps2=alpha) if alpha > 0 else None,
            delta_xx_meV=-2.0, gamma_xx_ratio=2.0, fss_ueV=0.0, transitions=frozenset({"X", "XX"}),
            dipole="in_plane", dEdT_meV_K=0.0, sigma_sd_ueV=0.5, tau_sd_ns=50.0)


@dataclass(frozen=True)
class TStruct(R.Structure):
    def response(self, vals, em, T_K, NA):
        g_ns, k_ns = vals["g_ns"], vals["kappa_ns"]
        kap = k_ns * HB if k_ns > 0 else None
        FP = vals["F_P"]
        if g_ns > 0:
            FP = 4 * g_ns ** 2 / (k_ns * em.gamma0_ns)
        return StructureResponse(F_P=FP, kappa_meV=kap, E_mode_meV=(em.E_X_meV if kap else None),
                                 g_meV=(g_ns * HB if g_ns > 0 else None), F_bg=vals["F_bg"],
                                 beta0=vals["beta0"], eta_out=vals["eta_out"])


@dataclass(frozen=True)
class TDrive(R.Drive):
    flavor: str = "pulse"

    def plan(self, vals, em, rt):
        if self.flavor == "pulse":
            p = Pulse(vals["shape"], vals["fwhm_ps"], vals["area_pi"] * math.pi,
                      detuning_meV=-vals["detuning_paper_meV"])           # repo det_X
            return DrivePlan("coherent", "pulsed", vals["f_rep_MHz"], pulses=(p,),
                             inputs={"laser_leak": 0.0})
        if self.flavor == "cw":
            return DrivePlan("coherent", "cw", None, omega_cw_ns=vals["omega_ns"],
                             inputs={"det_X_meV": 0.0})
        if self.flavor == "incoh_cw":
            return DrivePlan("incoherent", "cw", None, pump_rate_ns=vals["r_ns"], tau_on_ns=0.0, relax_ns=0.0)
        return DrivePlan("incoherent", "pulsed", vals["f_rep_MHz"], pump_rate_ns=vals["r_ns"],
                         tau_on_ns=vals["tau_on_ns"], relax_ns=vals["relax_ns"])


def reg_fixtures():
    R.register(TMat(kind="test_mat", label="t", summary="t", params=(
        P("lambda_nm", "nm", 900.0), P("gamma0_ns", "1/ns", 1.0), P("k_nr_ns", "1/ns", 0.0),
        P("gamma_star_ns", "1/ns", 0.0), P("alpha_ps2", "ps^2", 0.0), P("unused", "1", 1.0, tag="A"))))
    R.register(TStruct(kind="test_struct", label="t", summary="t", params=(
        P("F_P", "1", 0.0, tag="E"), P("g_ns", "1/ns", 0.0), P("kappa_ns", "1/ns", 0.0),
        P("F_bg", "1", 1.0), P("beta0", "1", 1.0), P("eta_out", "1", 1.0))))
    R.register(TDrive(kind="test_pulse", label="t", summary="t", flavor="pulse", params=(
        P("shape", "", "rect"), P("fwhm_ps", "ps", 0.5), P("area_pi", "1", 1.0),
        P("detuning_paper_meV", "meV", 0.0), P("f_rep_MHz", "MHz", 20.0))))
    R.register(TDrive(kind="test_cw", label="t", summary="t", flavor="cw", params=(
        P("omega_ns", "1/ns", 6.0),)))
    R.register(TDrive(kind="test_incoh", label="t", summary="t", flavor="incoh", params=(
        P("r_ns", "1/ns", 5.0), P("tau_on_ns", "ns", 0.05), P("relax_ns", "ns", 0.5),
        P("f_rep_MHz", "MHz", 20.0))))


def spec(mat=None, st=None, dr="test_pulse", drp=None, su=None, solver="auto", T=4.2):
    return {"schema": "fsim-composed/1", "T_K": {"value": T, "unit": "K", "tag": "V", "source": "t"},
            "material": {"kind": "test_mat", "params": mat or {}},
            "structure": {"kind": "test_struct", "params": st or {}},
            "drive": {"kind": dr, "params": drp or {}},
            "setup": {"kind": "lens_fibre_det", "params": su or {}}, "solver": solver}


def leaf(v, tag="V"):
    return {"value": v, "unit": "", "tag": tag, "source": "t"}


def main():
    reg_fixtures()
    exact = lambda v, unit, tag="V": {"value": v, "unit": unit, "tag": tag, "source": "t"}  # noqa: E731

    # 1. rect pi pulse g2 -> Gamma tau_p / 4
    for gt in (0.005, 0.01, 0.02):
        tau_ps = gt * 1e3                     # Gamma = 1/ns
        r = engine.evaluate_spec(spec(dr="test_pulse", drp={"shape": "rect", "fwhm_ps": tau_ps}))
        g2 = r["scalars"]["g2_0"]
        check(f"1. rect pi g2 ~ Gamma tau_p/4 at Gamma tau_p={gt}", near(g2, gt / 4, rel=0.10),
              f"g2={g2:.5g} vs {gt / 4:.5g}")

    # 2. instantaneous excitation I = Gamma/(Gamma+gamma*)
    gs = 0.4
    r = engine.evaluate_spec(spec(mat={"gamma_star_ns": exact(gs, "1/ns")},
                                  drp={"shape": "rect", "fwhm_ps": 0.1}))
    I = r["scalars"]["hom_I"]
    check("2. short-pulse HOM I = Gamma/(Gamma+gamma*_repo) to 1e-3", abs(I - 1.0 / (1.0 + gs)) < 1e-3,
          f"{I:.6f} vs {1 / (1 + gs):.6f}")

    # 3. Rabi
    ok, det = True, ""
    for a in (0.5, 1.0, 1.5):
        r = engine.evaluate_spec(spec(drp={"shape": "rect", "fwhm_ps": 0.2, "area_pi": a}))
        x = r["scalars"]["inversion"]
        ok &= abs(x - math.sin(a * math.pi / 2) ** 2) < 1e-3
        det += f" {x:.5f}"
    check("3. Rabi P_X = sin^2(area/2) for pi/2, pi, 3pi/2 to 1e-3", ok, det)

    # 4. detuned rect
    tau, dpap = 1e-3, 0.5
    r = engine.evaluate_spec(spec(drp={"shape": "rect", "fwhm_ps": 1.0, "area_pi": 1.0,
                                       "detuning_paper_meV": exact(dpap, "meV")}))
    Om, De = math.pi / tau, dpap / HB
    Lam = math.hypot(Om, De)
    want = (Om / Lam) ** 2 * math.sin(Lam * tau / 2) ** 2
    got = r["scalars"]["inversion"]
    check("4. detuned rect P_X = (Omega/Lambda)^2 sin^2(Lambda tau/2) to 1e-3", abs(got - want) < 1e-3,
          f"{got:.6f} vs {want:.6f}")

    # 5. cw RF Mollow
    Gam, Om = 2.0, 6.0
    r = engine.evaluate_spec(spec(mat={"gamma0_ns": exact(Gam, "1/ns")}, dr="test_cw",
                                  drp={"omega_ns": exact(Om, "1/ns")}))
    s = r["series"]["g2_tau"]
    t = np.array(s["x"])
    Op = math.sqrt(Om ** 2 - Gam ** 2 / 16)
    cf = 1 - np.exp(-0.75 * Gam * t) * (np.cos(Op * t) + 3 * Gam / (4 * Op) * np.sin(Op * t))
    err = float(np.max(np.abs(np.array(s["y"]) - cf)))
    check("5a. cw RF g2(tau) = Mollow closed form to 1e-6", err < 1e-6, f"max err {err:.2e}")
    check("5b. cw RF g2(0) = 0", abs(r["scalars"]["g2_0"]) < 1e-9, f"{r['scalars']['g2_0']:.2e}")

    # 6. beta / L / Gamma == resonant_source
    m, F, kap_ueV, det_ueV, b0, G0 = 1.5, 5.0, 150.0, 40.0, 0.1, 1.5
    src = rs.FockSource(gamma0_ns=G0, F_P=F, kappa_ueV=kap_ueV, cav_detune_ueV=det_ueV, beta0=b0,
                        hom_width_multiple=m, dephasing_mode="dot")
    sr_ = rs.source_rates(src)
    _, comps = rs.end_to_end_eta(src, 1.0)
    em = EmitterModel(E_X_meV=1380.0, lambda_nm=900, n_host=3.5, gamma0_ns=G0, k_nr_ns=0.0,
                      gamma_star_ns=sr_["gamma_star"], phonon=None, delta_xx_meV=-2.0, gamma_xx_ratio=2.0,
                      fss_ueV=0, transitions=frozenset({"X"}), dipole="x", dEdT_meV_K=0, sigma_sd_ueV=0,
                      tau_sd_ns=1)
    sr = StructureResponse(F_P=F, kappa_meV=kap_ueV * 1e-3, E_mode_meV=1380.0 - det_ueV * 1e-3, g_meV=None,
                           F_bg=1.0, beta0=b0)
    rt = rates.combine(em, sr, 4.0)
    check("6a. beta == end_to_end_eta beta to 1e-12", abs(rt.beta - comps["beta"]) < 1e-12,
          f"{rt.beta!r} {comps['beta']!r}")
    check("6b. L and Gamma == source_rates to 1e-11 (root conditioning; spec only fixes beta at 1e-12)",
          abs(rt.L - sr_["L"]) < 1e-11 and abs(rt.Gamma_ns - sr_["Gamma"]) < 1e-11,
          f"dL={rt.L - sr_['L']:.2e} dG={rt.Gamma_ns - sr_['Gamma']:.2e}")

    # 7. HOM correction
    v = OB.hom_v_corrected(0.925, 0.021)
    check("7a. V_corr(0.925, 0.021) = 0.966 +- 0.001", abs(v - 0.966) < 1e-3, f"{v:.5f}")
    # FX2c: the solver's I is the whole-stream overlap (multi-photon terms inside), so V_raw = I - g2 and the
    # g2-corrected value is I/(1 - g2) (7b'/7b'' keep the old (1 - g2) I - g2 form as the Ollivier M_s relation only)
    check("7b. FX2c V_raw = I - g2 (no background), V_corr = I/(1 - g2)",
          near(OB.hom_v_raw(0.9, 0.03), 0.9 - 0.03, ab=1e-15)
          and near(OB.hom_v_corrected(OB.hom_v_raw(0.9, 0.03), 0.03), 0.9 / 0.97, ab=1e-12))
    vb = OB.hom_v_raw(0.9, 0.05, 0.1, g2_int=0.03)
    check("7b2. FX2c background: V_raw = (rho^2 I - g2_raw)(1 - eps), rho^2 = (1 - g2_raw)/(1 - g2_int); the corrected "
          "value removes the background: V_corr = I/(1 - g2_int) for eps = 0",
          near(vb, (0.9 * 0.95 / 0.97 - 0.05) * 0.9, ab=1e-12)
          and near(OB.hom_v_corrected(OB.hom_v_raw(0.9, 0.05, 0.0, g2_int=0.03), 0.05), 0.9 / 0.97, ab=1e-12))

    # 7c. cavity sideband filter limits (review finding 4) and g2_raw without background (finding 11)
    from fsim_core import qd_gf as _qg
    _ph = _qg.PhononParams(alpha_ps2=0.0181)
    _Z = _qg.zpl_weight(_ph, 4.0)
    z0, ze_narrow = OB.cavity_sideband_filter(_ph, 4.0, 1e-6, 0.0, 1e-3)
    _, ze_wide = OB.cavity_sideband_filter(_ph, 4.0, 1e5, 0.0, 1e-3)
    check("7c. Z_eff limits: kappa -> 0 gives Z (ZPL filter), kappa -> inf gives 1 (unfiltered); Z <= Z_eff(150 ueV) <= 1",
          abs(ze_narrow - _Z) < 2e-4 and abs(ze_wide - 1.0) < 1e-4 and z0 == _Z
          and _Z <= OB.cavity_sideband_filter(_ph, 4.0, 0.15, 0.0, 1e-3)[1] <= 1.0,
          f"Z {_Z:.5f}, narrow {ze_narrow:.5f}, wide {ze_wide:.5f}")
    check("7d. g2_raw = g2 when the background is zero, also at zero signal",
          OB.g2_with_background(0.03, 0.0, 0.0) == 0.03 and OB.g2_with_background(0.03, 1e5, 0.0) == 0.03)

    # 8. incoherent S3d HOM vs closed form
    ok, det = True, ""
    for relax in (0.05, 2.0):
        gr, Gm, gst = 1 / relax, 1.0, 0.3
        r = engine.evaluate_spec(spec(mat={"gamma_star_ns": exact(gst, "1/ns")}, dr="test_incoh",
                                      drp={"relax_ns": exact(relax, "ns")}))
        got = r["scalars"]["hom_I"]
        want = gr / (gr + Gm) * Gm / (Gm + gst)
        ok &= r["scalars"]["solver"] == "lindblad_incoherent" and near(got, want, rel=0.02)
        det += f" tau_cap={relax}: {got:.5f}/{want:.5f}"
    check("8a. S3d HOM vs I=[g_r/(g_r+G)]G/(G+gamma*) within 2% (fast and slow capture)", ok, det)
    # no reservoir: counting hierarchy == pulse_counting.pulse_g2
    sp = spec(dr="test_incoh", drp={"relax_ns": exact(0.0, "ns")}, solver="lindblad_incoherent")
    g_l = engine.evaluate_spec(sp)["scalars"]["g2_0"]
    g_r = engine.evaluate_spec({**sp, "solver": "rate"})["scalars"]["g2_0"]
    check("8b. lindblad_incoherent g2 (no reservoir) == pulse_g2 to 1e-8", abs(g_l - g_r) < 1e-8,
          f"{g_l!r} {g_r!r}")

    # 9. JC
    G0j, g, kp = 1.0, 5.0, 200.0
    r = engine.evaluate_spec(spec(mat={"gamma0_ns": exact(G0j, "1/ns")},
                                  st={"g_ns": exact(g, "1/ns"), "kappa_ns": exact(kp, "1/ns")},
                                  drp={"fwhm_ps": 1.0}, solver="jc"))
    want = 1e3 / (G0j * (1 + 4 * g * g / (kp * G0j)))
    got = r["scalars"]["tau_ps"]
    check("9a. JC bad-cavity lifetime = 1/(G0(1+4g^2/(kappa G0))) within 2%", near(got, want, rel=0.02),
          f"{got:.2f} vs {want:.2f} ps")
    gT, kT, yT = 2 * math.pi * 4.16, 2 * math.pi * 23.3, 2 * math.pi * 0.30
    r = engine.evaluate_spec(spec(mat={"gamma0_ns": exact(yT, "1/ns")},
                                  st={"g_ns": exact(gT, "1/ns"), "kappa_ns": exact(kT, "1/ns")},
                                  drp={"fwhm_ps": 3.0, "shape": "gauss"}, solver="jc"))
    sc = r["scalars"]
    print(f"      Tomm parameters: JC tau = {sc['tau_ps']:.2f} ps, Purcell tau = {sc['tau_purcell_ps']:.2f} ps, "
          f"F_P_eff = {sc['F_P_eff']:.2f}, 4g/kappa flag = {sc['jc_recommended']}")
    check("9b. Tomm: JC and Purcell lifetimes differ by < 20%",
          abs(sc["tau_ps"] - sc["tau_purcell_ps"]) / sc["tau_purcell_ps"] < 0.20)
    check("9c. Tomm: jc_recommended (4g/kappa = 0.71) and auto picks jc",
          sc["jc_recommended"] and engine.evaluate_spec({**spec(
              mat={"gamma0_ns": exact(yT, "1/ns")}, st={"g_ns": exact(gT, "1/ns"), "kappa_ns": exact(kT, "1/ns")},
              drp={"fwhm_ps": 3.0, "shape": "gauss"})})["scalars"]["solver"] == "jc")

    # 9d. Tomm 47.5 ps (measured, Nat. Nanotechnol. 16, 399 (2021) [V]); tolerance PRE-REGISTERED at +-15% before the
    # run (tol_reason: the paper's own Purcell-model F_P and the dot-to-dot spread; the JC model uses the entered g,
    # kappa, gamma only). Not tuned: the value is reported either way.
    tau_T = sc["tau_ps"]
    print(f"      Tomm: tau_ps (slowest Liouvillian eigenvalue) = {tau_T:.2f} ps, tau_1e_ps = {sc['tau_1e_ps']:.2f} ps, "
          f"measured 47.5 ps, Purcell {sc['tau_purcell_ps']:.2f} ps")
    check("9d. Tomm: JC tau_ps within +-15% of the measured 47.5 ps (pre-registered)", abs(tau_T / 47.5 - 1.0) <= 0.15,
          f"{tau_T:.2f} ps, deviation {100 * (tau_T / 47.5 - 1):+.1f}%")
    # independent: single-excitation amplitude matrix M = [[-G/2, -i g], [-i g, -kappa/2]] (own code, no Lindblad);
    # populations decay at -2 Re(lambda_M); the verifier material has gamma* = 0 (population modes unaffected by dephasing).
    Mmat = np.array([[-yT / 2, -1j * gT], [-1j * gT, -kT / 2]])
    tau_M = 1e3 / (-2.0 * np.linalg.eigvals(Mmat).real.max())
    check("9e. tau_ps == 1/(2|Re lambda_M|) of the 2x2 non-Hermitian JC matrix (independent, rel 1e-3)",
          abs(tau_T / tau_M - 1.0) < 1e-3, f"{tau_T:.3f} vs {tau_M:.3f} ps")
    check("9f. tau_1e_ps (old 1/e-after-peak) is reported and differs from the eigenvalue lifetime",
          math.isfinite(sc["tau_1e_ps"]) and abs(sc["tau_1e_ps"] - tau_T) > 1.0, f"{sc['tau_1e_ps']:.2f} ps")
    check("9g. F_P_peak is the structure's F_P (passthrough) and F_P_eff = F_P_peak x L",
          abs(sc["F_P_eff"] - sc["F_P_peak"] * sc["L"]) < 1e-12 * max(1.0, sc["F_P_peak"]), f"{sc['F_P_peak']:.3f} {sc['F_P_eff']:.3f}")

    # 10. tags and envelope
    sp = spec(mat={"lambda_nm": exact(900.0, "nm", "V"), "gamma0_ns": exact(1.0, "1/ns", "V")},
              st={"F_P": {"range": [2.0, 6.0], "unit": "1", "tag": "E", "source": "t"}},
              drp={"shape": "rect", "fwhm_ps": 2.0})
    r = engine.evaluate_spec(sp)
    pv = r["scalars"]["provenance"]
    check("10a. [V]+[E] inputs -> [E]; passthrough lambda_nm stays [V]",
          pv["tau_ps"]["tag"] == "E" and pv["lambda_nm"]["tag"] == "V", f"{pv['tau_ps']['tag']} {pv['lambda_nm']['tag']}")
    check("10b. unread param does not widen tags; range_note set",
          "material.unused" not in pv["tau_ps"]["inputs"] and r["scalars"]["range_note"] == "midpoint, not a prediction")
    env = engine.envelope(sp, ["tau_ps", "g2_0", "B_lens"])
    ok = True
    for k in ("tau_ps", "g2_0", "B_lens"):
        lo, hi = env[k]["range"]
        ok &= lo <= env[k]["central"] <= hi and lo < hi and "mean" not in env[k]
    lo2 = engine.evaluate_spec({**sp, "structure": {"kind": "test_struct", "params": {"F_P": exact(2.0, "1", "E")}}})
    hi2 = engine.evaluate_spec({**sp, "structure": {"kind": "test_struct", "params": {"F_P": exact(6.0, "1", "E")}}})
    ok &= near(env["tau_ps"]["range"][1], lo2["scalars"]["tau_ps"], rel=1e-9) and \
        near(env["tau_ps"]["range"][0], hi2["scalars"]["tau_ps"], rel=1e-9)
    check("10c. envelope [min,max] contains central, matches endpoint runs, no mean key", ok, str(env["tau_ps"]))
    sr_ = engine.series(sp, "structure.F_P", [0.0, 1.0, 2.0], ["tau_ps"])
    check("10d. series arrays", len(sr_["y"]["tau_ps"]) == 3 and sr_["y"]["tau_ps"][0] > sr_["y"]["tau_ps"][2])

    # phonon bookkeeping: Z^2 unfiltered, Z filtered, never both
    a = 0.027
    Zt = qd_gf.zpl_weight(qd_gf.PhononParams(alpha_ps2=a), 4.2)
    base = engine.evaluate_spec(spec(drp={"shape": "rect", "fwhm_ps": 0.1}))["scalars"]
    un = engine.evaluate_spec(spec(mat={"alpha_ps2": exact(a, "ps^2")}, drp={"shape": "rect", "fwhm_ps": 0.1}))["scalars"]
    fi = engine.evaluate_spec(spec(mat={"alpha_ps2": exact(a, "ps^2")}, su={"zpl_filter": True},
                                   drp={"shape": "rect", "fwhm_ps": 0.1}))["scalars"]
    check("10e. phonons: hom x Z^2 unfiltered; filtered hom x1, B_lens x Z",
          near(un["hom_I"], base["hom_I"] * Zt ** 2, rel=1e-9) and near(un["B_lens"], base["B_lens"], rel=1e-9)
          and near(fi["hom_I"], base["hom_I"], rel=1e-9) and near(fi["B_lens"], base["B_lens"] * Zt, rel=1e-9))
    lk = engine.evaluate_spec(spec(su={"laser_leak": 0.01}, drp={"shape": "rect", "fwhm_ps": 0.1}))["scalars"]
    check("10f. g2_raw ~ g2 + 2 f_leak", near(lk["g2_raw"], lk["g2_0"] + 0.02, ab=3e-4),
          f"{lk['g2_raw']:.5f} vs {lk['g2_0'] + 0.02:.5f}")

    # select: phonons unsupported
    plan = DrivePlan("coherent", "pulsed", 20.0, pulses=(Pulse("gauss", 10.0, 5 * math.pi),),
                     needs=frozenset({"phonons"}))
    sv, iss = solvers.select(plan, em, rt, "auto")
    check("10g. phonon plan -> no solver, issue 'needs P9 polaron solver'",
          sv is None and any("P9 polaron" in i for i in iss), str(iss))

    # 11. contract + JSON
    r = engine.evaluate_spec(spec(drp={"shape": "rect", "fwhm_ps": 2.0}), T_grid=[4.0, 10.0])
    keys = {"curves", "scalars", "series"}
    cv = r["curves"]
    check("11a. contract keys", set(r) == keys and {"T_hs", "Tj", "g2", "hom_I", "B_first_lens", "Gamma_ns"} <= set(cv)
          and len(cv["g2"]) == 2 and {"solver", "regime", "valid", "invalid_reasons", "tag_chain", "provenance"} <= set(r["scalars"]))
    try:
        json.loads(dumps(r))
        okj = True
    except Exception as e:                                                  # noqa: BLE001
        okj = False
        print(e)
    check("11b. JSON-serialisable via fsim_studio.serialize.dumps", okj)
    check("11c. estimate > 0", engine.estimate(spec(), 3) > 0)

    # real registry entries end to end
    try:
        real = {"schema": "fsim-composed/1", "T_K": 4.2,
                "material": {"kind": "inas_gaas_sk", "params": {}},
                "structure": {"kind": "micropillar", "params": {"Q_measured": 6124}},
                "drive": {"kind": "resonant_pulse", "params": {}},
                "setup": {"kind": "lens_fibre_det", "params": {}}, "solver": "auto"}
        rr = engine.evaluate_spec(real)["scalars"]
        fin = all(math.isfinite(rr[k]) for k in ("g2_0", "hom_I", "B_lens", "tau_ps", "inversion"))
        tpe = dict(real, material={"kind": "gaas_algaas_droplet", "params": {}},
                   structure={"kind": "bulk_planar", "params": {}}, drive={"kind": "tpe", "params": {}})
        rt_ = engine.evaluate_spec(tpe)["scalars"]
        check("12a. real inas_gaas_sk + micropillar + resonant_pulse finite", fin, str({k: rr[k] for k in ("solver", "g2_0", "hom_I")}))
        check("12b. real droplet + bulk + tpe runs (lindblad cascade)", rt_["solver"] == "lindblad" and math.isfinite(rt_["g2_0"]),
              str({k: rt_[k] for k in ("solver", "g2_0", "hom_I")}))
    except Exception as e:                                                  # noqa: BLE001
        check("12. real registry integration", False, f"{type(e).__name__}: {e}")

    # 13. devlib-fx2a observables
    drp76 = {"shape": "rect", "fwhm_ps": 0.1, "f_rep_MHz": 76.0}
    S_full = engine.evaluate_spec(spec(drp=drp76))["scalars"]["rate_detected_Hz"]
    su_lo = {"eta_optics": 0.5 * 2750.0 / S_full}                    # Muller-like detected rate S = 2750/s (S is linear in eta_optics)
    base_d = engine.evaluate_spec(spec(drp=drp76, su=su_lo))["scalars"]
    S_det = base_d["rate_detected_Hz"]
    for win, tag_ in ((2.0, "a"), (0.5, "b")):
        dk = engine.evaluate_spec(spec(drp=drp76, su={**su_lo, "dark_cps": 250.0, "coinc_window_ns": win}))["scalars"]
        B_h = 250.0 * win * 1e-9 * 76e6                                   # hand: dark per second inside the window
        rho_h = S_det / (S_det + B_h)
        g_h = 1.0 - rho_h ** 2 * (1.0 - base_d["g2_0"])
        check(f"13{tag_}. pulsed g2_raw: darks only inside the {win} ns window (hand formula)",
              near(dk["g2_raw"], g_h, rel=1e-9) and dk["g2_raw"] - base_d["g2_0"] > 1e-3,
              f"S={S_det:.0f}/s raw {dk['g2_raw']:.6f} vs hand {g_h:.6f}")
    z0 = engine.evaluate_spec(spec(drp=drp76, su={**su_lo, "coinc_window_ns": 5.0}))["scalars"]
    check("13c. zero dark counts: pulsed g2_raw equals the intrinsic g2 whatever the window",
          z0["g2_raw"] == z0["g2_0"])
    hi = engine.evaluate_spec(spec(drp=drp76, su={**su_lo, "dark_cps": 250.0, "coinc_window_ns": 1.0e6}))["scalars"]     # window > period -> capped at the period
    rho_c = S_det / (S_det + 250.0)
    check("13d. window longer than the pulse period is capped at the period (full dark rate)",
          near(hi["g2_raw"], 1.0 - rho_c ** 2 * (1.0 - base_d["g2_0"]), rel=1e-9))
    # Muller-like numbers through the helper: S 2750/s, dark 250/s, 2 ns, 76 MHz
    Sm, Bm = 2750.0, OB.pulsed_dark_rate(250.0, 2.0, 76e6)
    inc_win = OB.g2_with_background(0.0, Sm, Bm)
    inc_old = OB.g2_with_background(0.0, Sm, 250.0)
    check("13e. Muller-like (S 2750/s, dark 250/s, 76 MHz, 2 ns): windowed increase 0.027 << unwindowed 0.16",
          near(inc_win, 0.0272, ab=5e-4) and near(inc_old, 0.1597, ab=5e-4), f"{inc_win:.4f} vs {inc_old:.4f}")
    cw_ = engine.evaluate_spec(spec(dr="test_cw", su={"dark_cps": 250.0, "coinc_window_ns": 0.01}))["scalars"]
    Sc = cw_["rate_detected_Hz"]
    check("13f. cw drive keeps the full dark rate (window ignored; or both undefined without a flux)",
          near(cw_["g2_raw"], OB.g2_with_background(cw_["g2_0"], Sc, 250.0), rel=1e-9) if math.isfinite(Sc)
          else not math.isfinite(cw_["g2_raw"]), f"S={Sc}")
    eps_r = engine.evaluate_spec(spec(drp={"shape": "rect", "fwhm_ps": 0.1}, su={"eps_int": 0.1}))["scalars"]
    e0_r = engine.evaluate_spec(spec(drp={"shape": "rect", "fwhm_ps": 0.1}))["scalars"]
    check("13g. eps_int ParamSpec exists on lens_fibre_det and scales hom_V_raw by (1 - eps_int)",
          R.get("setup", "lens_fibre_det").param("eps_int").default == 0.0
          and near(eps_r["hom_V_raw"], e0_r["hom_V_raw"] * 0.9, rel=1e-9))
    # blinking: real material (has eta_blink); brightness scales, lifetime does not
    def real(blink):
        return {"schema": "fsim-composed/1", "T_K": 4.2,
                "material": {"kind": "inas_gaas_sk", "params": {"eta_blink": {"value": blink, "unit": "1", "tag": "A", "source": "t"}}},
                "structure": {"kind": "micropillar", "params": {"Q_measured": 6124}},
                "drive": {"kind": "resonant_pulse", "params": {}},
                "setup": {"kind": "lens_fibre_det", "params": {}}, "solver": "auto"}
    b1, b6 = engine.evaluate_spec(real(1.0))["scalars"], engine.evaluate_spec(real(0.65))["scalars"]
    check("13h. eta_blink scales B_lens / B_fibre / detected rate, not tau, g2 or n_ph",
          near(b6["B_lens"], 0.65 * b1["B_lens"], rel=1e-9) and near(b6["rate_detected_Hz"], 0.65 * b1["rate_detected_Hz"], rel=1e-9)
          and b6["tau_ps"] == b1["tau_ps"] and b6["g2_0"] == b1["g2_0"])
    sp6 = engine._as_spec(real(0.65))
    ent6 = engine._entries(sp6)
    c6 = engine._run_point(ent6, engine._central_values(engine._leaves(sp6, ent6)), 4.2, "auto")
    check("13i. eta_rad_eff = eta_blink x eta_rad is reported in the chain (eta_rad unchanged)",
          near(c6.d["eta_rad_eff"], 0.65 * c6.d["eta_rad"], rel=1e-12) and 0.0 < c6.d["eta_rad"] <= 1.0)
    check("13j. eta_blink is a registered material param (default 1, 0.5-1) on every material kind",
          all(R.get("material", k).param("eta_blink").default == 1.0 and R.get("material", k).param("eta_blink").lo == 0.5
              for k in R.kinds("material") if not k.startswith("test_")))
    # beta_entered: passthrough, declared, NaN when the structure has no entered beta
    bo = OB.get_observable("beta_entered")
    check("13k. beta_entered declares passthrough_of structure.beta_direct", bo.passthrough_of == "structure.beta_direct")
    phc = dict(real(1.0), structure={"kind": "phc_waveguide", "params": {"beta_direct": {"value": 0.98, "unit": "1", "tag": "V", "source": "t"}}})
    pe = engine.evaluate_spec(phc)
    envp = engine.envelope(phc, ["beta_entered"])
    check("13l. beta_entered returns the entered beta (0.98) with tag V; NaN for a structure without one",
          near(envp["beta_entered"]["central"], 0.98, ab=1e-12) and envp["beta_entered"]["tag"] == "V"
          and not math.isfinite(engine.envelope(real(1.0), ["beta_entered"])["beta_entered"]["central"]))

    # ------------------------------------------------------------------ 14. FX2b: eta_prep, XX channel, refill / background / thermal
    from scipy.integrate import solve_ivp
    from fsim_core.compose.rates import RatesEx as _RatesEx
    from fsim_core.spectral import KB as _KB

    def bloch_prep(shape, fwhm_ps, area_pi, G, k, gs):
        """INDEPENDENT Bloch-equation integral of (P_X + emitted photons) over the pulse window, resonant drive:
        dP = -(G+k) P + Om y, dy = -(G+k+gs)/2 y - (Om/2)(2P - 1), dE = G P (own code, no Lindblad module)."""
        dur = fwhm_ps * 1e-3
        if shape == "rect":
            win, Om = dur, (lambda t: area_pi * math.pi / dur)
        else:
            sg = dur / (2.0 * math.sqrt(math.log(2.0)))
            win = 2.0 * rs.GAUSS_WINDOW_SIGMA * sg
            Om = lambda t: area_pi * math.pi / (sg * math.sqrt(2 * math.pi)) * math.exp(-0.5 * ((t - win / 2) / sg) ** 2)
        gc = 0.5 * (G + k + gs)
        f = lambda t, y: [-(G + k) * y[0] + Om(t) * y[1], -gc * y[1] - 0.5 * Om(t) * (2 * y[0] - 1.0), G * y[0]]
        s = solve_ivp(f, [0.0, win], [0.0, 0.0, 0.0], method="DOP853", rtol=1e-12, atol=1e-14, max_step=win / 200)
        return float(s.y[0, -1]), float(s.y[2, -1])

    def lind_plan(shape, fwhm_ps, area_pi, f=80.0, needs=()):
        return DrivePlan("coherent", "pulsed", f, pulses=(Pulse(shape, fwhm_ps, area_pi * math.pi),),
                         needs=frozenset(needs), inputs={})

    em0 = R.get("material", "test_mat").build({"lambda_nm": 900.0, "gamma0_ns": 1.0, "k_nr_ns": 0.0,
                                               "gamma_star_ns": 0.0, "alpha_ps2": 0.0}, 4.0)

    def rtx(G, gs=0.0, k=0.0, G0=None):
        return _RatesEx(Gamma_ns=G, gamma_star_ns=gs, k_nr_ns=k, L=1.0, F_P_eff=0.0, beta=0.5,
                        Gamma0_ns=G if G0 is None else G0)
    LS = solvers.LindbladSolver()
    ok, det = True, ""
    for shape, fw, G, k, gs in (("rect", 100.0, 10.0, 0.5, 3.0), ("rect", 60.0, 4.0, 0.0, 0.3), ("gauss", 20.0, 8.0, 0.2, 1.0)):
        r = LS.run(lind_plan(shape, fw, 1.0), em0, rtx(G, gs, k))
        P_end, E = bloch_prep(shape, fw, 1.0, G, k, gs)
        tol = 2e-6 if shape == "rect" else 3e-5       # gauss: midpoint rule at 8 steps/sigma (the pre-existing P_X error is 1.1e-5)
        ok &= near(r.extra["n_exc"], P_end + E, ab=tol) and near(r.extra["P_X_end"], P_end, ab=tol)
        ok &= r.extra["n_exc"] > r.extra["P_X_end"] + 1e-3
        det += f" {shape}{fw:g}ps: {r.extra['n_exc']:.6f}/{P_end + E:.6f}"
    check("14a. FX2b/FX2c lindblad n_exc (the expected count P_X(end) + emitted) equals an independent Bloch ODE (rect + "
          "gauss, decay, dephasing, k_nr) to 2e-6 (rect) / 3e-5 (gauss, grid error); P_X_end keeps the old end-of-window "
          "population", ok, det)
    ok, det = True, ""
    for fw in (0.2, 0.5, 1.0):
        G = 2.0
        r = LS.run(lind_plan("rect", fw, 1.0), em0, rtx(G))
        gt = G * fw * 1e-3
        e_new, e_old = abs(1.0 - r.extra["eta_prep"]), abs(1.0 - r.extra["P_X_end"])      # eta_prep = P(n >= 1) now
        ok &= e_new <= 0.5 * gt and e_new < 0.5 * e_old
        det += f" Gtau={gt:.0e}: |1-eta_prep|={e_new:.1e} vs |1-P_X_end|={e_old:.1e};"
    check("14b. FX2b short pi pulse: |1 - eta_prep| <= Gamma tau_p/2 and below half the old |1 - P_X_end|", ok, det)
    r = LS.run(lind_plan("rect", 60.0, 1.0), em0, rtx(4.0))
    check("14c. count conservation: k_nr = 0 -> n_ph = n_exc (every prepared excitation is emitted), and eta_prep = "
          "P(n >= 1) <= min(1, n_exc)",
          near(r.n_ph, r.extra["n_exc"], ab=1e-9) and r.extra["eta_prep"] <= min(1.0, r.extra["n_exc"]) + 1e-12,
          f"n_ph {r.n_ph:.10f} vs n_exc {r.extra['n_exc']:.10f}, eta_prep {r.extra['eta_prep']:.10f}")

    # 14d. jc: bad-cavity limit eta_prep == lindblad with the Purcell rate (validated path); the old P_X_end differs
    gj, kj = 5.0, 200.0
    stj = {"g_ns": exact(gj, "1/ns"), "kappa_ns": exact(kj, "1/ns")}
    aj = engine.evaluate_spec(spec(st=stj, drp={"fwhm_ps": 100.0}, solver="jc"))["scalars"]
    bj = engine.evaluate_spec(spec(st=stj, drp={"fwhm_ps": 100.0}, solver="lindblad"))["scalars"]
    check("14d. FX2b jc eta_prep (P_X + <a+a> + emitted) = Purcell-lindblad eta_prep to 1e-3 in the bad-cavity limit, "
          "while the old jc P_X_end is > 0.05 away",
          abs(aj["eta_prep"] - bj["eta_prep"]) < 1e-3 and abs(aj["P_X_end"] - bj["eta_prep"]) > 0.05,
          f"jc {aj['eta_prep']:.5f}, lindblad {bj['eta_prep']:.5f}, jc P_X_end {aj['P_X_end']:.5f}")
    cj = engine.evaluate_spec(spec(st={"g_ns": exact(2 * math.pi * 4.16, "1/ns"), "kappa_ns": exact(2 * math.pi * 23.3, "1/ns")},
                                   mat={"gamma0_ns": exact(2 * math.pi * 0.30, "1/ns")},
                                   drp={"fwhm_ps": 3.0, "shape": "gauss"}, solver="jc"))["scalars"]
    sp_cj = engine._as_spec(spec(st={"g_ns": exact(2 * math.pi * 4.16, "1/ns"), "kappa_ns": exact(2 * math.pi * 23.3, "1/ns")},
                                 mat={"gamma0_ns": exact(2 * math.pi * 0.30, "1/ns")},
                                 drp={"fwhm_ps": 3.0, "shape": "gauss"}, solver="jc"))
    ent_cj = engine._entries(sp_cj)
    c_cj = engine._run_point(ent_cj, engine._central_values(engine._leaves(sp_cj, ent_cj)), 4.2, "jc")
    check("14e. FX2b/FX2c jc (Tomm-class coupling): P_X_end <= eta_prep <= 1 and n_ph = n_exc (k_nr = 0); the observable "
          "n_exc is exposed",
          cj["P_X_end"] <= cj["eta_prep"] <= 1.0 and near(cj["n_ph"], c_cj.res.extra["n_exc"], ab=1e-6)
          and near(c_cj.d["n_exc"], c_cj.res.extra["n_exc"], ab=1e-15),
          f"eta_prep {cj['eta_prep']:.5f}, P_X_end {cj['P_X_end']:.5f}, n_ph {cj['n_ph']:.5f}, "
          f"n_exc {c_cj.res.extra['n_exc']:.5f}")

    # 14f-j. XX channel of the TPE cascade
    GX, ratio, gsx = 4.0, 2.0, 0.5
    emx = R.get("material", "test_mat").build({"lambda_nm": 900.0, "gamma0_ns": GX, "k_nr_ns": 0.0,
                                               "gamma_star_ns": gsx, "alpha_ps2": 0.0}, 4.0)
    emx = EmitterModel(**{**emx.__dict__, "delta_xx_meV": -4.0, "gamma_xx_ratio": ratio})
    GXX = ratio * GX
    gx = {fw: LS.run(lind_plan("gauss", fw, 1.0, needs=("cascade",)), emx, rtx(GX, gsx, G0=GX)) for fw in (3.0, 12.0)}
    r3 = gx[3.0]
    check("14f. TPE tau_XX_ps = 1e3/Gamma_XX (read off the cascade Liouvillian, 1e-9)",
          near(r3.extra["tau_XX_ps"], 1e3 / GXX, rel=1e-9), f"{r3.extra['tau_XX_ps']:.6f} vs {1e3 / GXX:.6f} ps")
    check("14g. TPE g2_XX -> 0 for short pulses: < 5e-3 at Gamma tau_p = 0.012 and below the 12 ps value",
          r3.extra["g2_XX"] < 5e-3 and r3.extra["g2_XX"] < gx[12.0].extra["g2_XX"],
          f"3 ps {r3.extra['g2_XX']:.2e}, 12 ps {gx[12.0].extra['g2_XX']:.2e}")
    want_xx = GXX / (GXX + GX + gsx)
    check("14h. TPE hom_XX (instantaneous prep) = Gamma_XX/(Gamma_XX + Gamma_X + gamma*) to 1e-6 (derived: the "
          "|X><XX| coherence decays at (Gamma_XX + Gamma_X + gamma*)/2, intensity rate Gamma_XX)",
          near(r3.extra["hom_XX"], want_xx, ab=1e-6), f"{r3.extra['hom_XX']:.8f} vs {want_xx:.8f}")
    emy = EmitterModel(**{**emx.__dict__, "gamma0_ns": 0.01, "gamma_xx_ratio": 800.0})
    ry = LS.run(lind_plan("gauss", 3.0, 1.0, needs=("cascade",)), emy, rtx(0.01, 0.5, G0=0.01))
    check("14i. TPE hom_XX -> Gamma_XX/(Gamma_XX + gamma*) when the intermediate X is long lived (Gamma_X << Gamma_XX): "
          "the spec's form is the Gamma_X -> 0 limit", near(ry.extra["hom_XX"], 8.0 / 8.5, rel=2e-3),
          f"{ry.extra['hom_XX']:.5f} vs {8.0 / 8.5:.5f}")
    tpe_spec = {"schema": "fsim-composed/1", "T_K": {"value": 4.0, "unit": "K", "tag": "V", "source": "t"},
                "material": {"kind": "gaas_algaas_droplet", "params": {}},
                "structure": {"kind": "bulk_planar", "params": {}},
                "drive": {"kind": "tpe", "params": {"shape": {"value": "gauss", "unit": "", "tag": "A", "source": "t"},
                                                    "fwhm_ps": exact(5.0, "ps", "A")}},
                "setup": {"kind": "lens_fibre_det", "params": {}}, "solver": "auto"}
    ts = engine.evaluate_spec(tpe_spec)["scalars"]
    em_t = R.get("material", "gaas_algaas_droplet").build({}, 4.0)
    check("14j. engine surfaces the XX channel (tau_XX_ps, g2_XX, hom_XX, n_XX) for TPE only",
          all(math.isfinite(ts[k]) for k in ("tau_XX_ps", "g2_XX", "hom_XX", "n_XX", "P_X_end", "eta_prep"))
          and near(ts["tau_XX_ps"], 1e3 / (em_t.gamma_xx_ratio * em_t.gamma0_ns), rel=1e-9)
          and "g2_XX" not in engine.evaluate_spec(spec())["scalars"],
          f"tau_XX {ts['tau_XX_ps']:.2f} ps, g2_XX {ts['g2_XX']:.2e}, hom_XX {ts['hom_XX']:.4f}, "
          f"eta_prep {ts['eta_prep']:.4f}, P_XX_end {ts['P_X_end']:.4f}")

    # 14k-r. incoherent drives: refill, background, thermal escape (S3d solver)
    IS = solvers.IncoherentLindbladSolver()

    def inc_plan(P, win_ps, relax_ps, f=10.0, **inp):
        tau = win_ps * 1e-3
        return DrivePlan("incoherent", "pulsed", f, pump_rate_ns=P / tau, tau_on_ns=tau, pump_ratio=P,
                         relax_ns=relax_ps * 1e-3, inputs=dict(inp))
    rt1 = rtx(1.0)
    base = IS.run(inc_plan(30, 1.0, 10.0), None, rt1)
    same = IS.run(inc_plan(30, 1.0, 10.0, p_refill=0.0, bg_frac=0.0, k0_ns=0.0, E_act_meV=25.0, T_K=300.0), None, rt1)
    check("14k. p_refill = bg_frac = k0 = 0 reproduces the original S3d result bit-for-bit (g2, n_ph, I)",
          (base.g2_0, base.n_ph, base.I_single) == (same.g2_0, same.n_ph, same.I_single))
    ok, det = True, ""
    p1 = 1.0 - math.exp(-30.0)
    for p in (0.02, 0.1, 0.3):
        r = IS.run(inc_plan(30, 1.0, 10.0, p_refill=p, tau_refill_ps=100.0), None, rt1)
        ok &= near(r.g2_0, 2 * p / p1, rel=0.01, ab=2e-4) and near(r.n_ph, p1 / (1 - p), rel=2e-3)
        det += f" p={p}: g2 {r.g2_0:.4f}/{2 * p / p1:.4f}, n {r.n_ph:.4f}/{p1 / (1 - p):.4f};"
    check("14l. trap refilling: g2 = 2 p_refill / p_loaded and n_ph = p_loaded/(1 - p_refill) (geometric refill chain, "
          "closed form) within 1% / 0.2%", ok, det)
    ok = True
    for bg in (0.1, 0.3):
        r = IS.run(inc_plan(30, 1.0, 10.0, bg_frac=bg), None, rt1)
        ok &= near(r.g2_0, 1.0 - (1.0 - bg) ** 2 * (1.0 - base.g2_0), ab=1e-12) and near(r.g2_0, 1.0 - (1.0 - bg) ** 2, ab=1e-3)
    check("14m. background fraction: g2 = 1 - rho^2 (1 - g2_signal), rho = 1 - bg_frac (Brouri form), -> 1 - rho^2 for "
          "a perfect source", ok)
    Eth, Tth = 12.0, 77.0
    k_arr = 5.0 * math.exp(-Eth / (_KB * Tth))
    r = IS.run(inc_plan(30, 1.0, 10.0, k0_ns=5.0, E_act_meV=Eth, T_K=Tth), None, rt1)
    check("14n. Arrhenius escape: reported k_esc = k0 exp(-E_act/kT) (independent exp)",
          near(r.extra["k_esc_ns"], k_arr, rel=1e-12), f"{r.extra['k_esc_ns']:.6g} vs {k_arr:.6g} /ns")
    kk, period = 1e-3, 1e3 / 10.0
    r = IS.run(inc_plan(30, 1.0, 10.0, k0_ns=kk, E_act_meV=0.0, T_K=50.0), None, rt1)
    lam = kk * period
    want = (2 * p1 * lam + lam ** 2) / (p1 + lam) ** 2
    check("14o. thermal recapture: pulsed g2 = (2 p1 lam + lam^2)/(p1 + lam)^2 with lam = k_esc T_rep (Poisson background "
          "on a single-photon pulse) within 3% (dead time ~ 1/(Gamma T_rep) + k/Gamma)", near(r.g2_0, want, rel=0.03),
          f"{r.g2_0:.5f} vs {want:.5f}")
    g2T = [IS.run(inc_plan(30, 1.0, 10.0, f=76.0, k0_ns=0.05, E_act_meV=15.0, T_K=T), None, rt1).g2_0 for T in (4.0, 40.0, 77.0)]
    base76 = IS.run(inc_plan(30, 1.0, 10.0, f=76.0), None, rt1)
    check("14p. thermal escape: g2 grows with T (Arrhenius) and the 4 K limit (k_esc ~ 1e-20/ns) reproduces the "
          "T-independent g2 to 1e-9", g2T[0] < g2T[1] < g2T[2] and abs(g2T[0] - base76.g2_0) < 1e-9,
          f"g2(4,40,77 K) = {g2T[0]:.3e}, {g2T[1]:.3e}, {g2T[2]:.3e}; no-thermal {base76.g2_0:.3e}")
    d_ps = R.get("drive", "p_shell")
    d_ab = R.get("drive", "above_band")
    check("14q. drive params: p_shell has p_refill/tau_refill_ps/bg_frac/E_act_meV/k0_ns (default 0 = off), above_band "
          "has refill + background but no thermal params; the plan carries them",
          all(d_ps.param(n).default == 0.0 for n in ("p_refill", "bg_frac", "E_act_meV", "k0_ns"))
          and d_ps.param("tau_refill_ps").default > 0 and all(d_ab.param(n).default == 0.0 for n in ("p_refill", "bg_frac"))
          and not any(p.name in ("E_act_meV", "k0_ns") for p in d_ab.params)
          and d_ps.plan({"p_refill": 0.05, "k0_ns": 3.0, "E_act_meV": 10.0}, em0, rtx(1.0)).inputs["k0_ns"] == 3.0)

    def shell_spec(T, **dp):
        pars = {"pump": 5.0, "window_ps": 3.0, "f_rep_MHz": 76.0, **dp}
        return {"schema": "fsim-composed/1", "T_K": {"value": T, "unit": "K", "tag": "V", "source": "t"},
                "material": {"kind": "inas_gaas_sk", "params": {}},
                "structure": {"kind": "bulk_planar", "params": {}},
                "drive": {"kind": "p_shell", "params": {k: exact(v, d_ps.param(k).unit, "A") for k, v in pars.items()}},
                "setup": {"kind": "lens_fibre_det", "params": {}}, "solver": "auto"}
    g_a = [engine.evaluate_spec(shell_spec(T, k0_ns=0.05, E_act_meV=15.0))["scalars"]["g2_0"] for T in (4.0, 40.0, 77.0)]
    g_0 = [engine.evaluate_spec(shell_spec(T))["scalars"]["g2_0"] for T in (4.0, 77.0)]
    check("14r. engine p_shell: g2(T) rises with the Arrhenius escape and is flat without it (T from the engine)",
          g_a[0] < g_a[1] < g_a[2] and abs(g_0[0] - g_0[1]) < 1e-12 and abs(g_a[0] - g_0[0]) < 1e-9,
          f"with escape {g_a[0]:.3e}/{g_a[1]:.3e}/{g_a[2]:.3e}, without {g_0[0]:.3e}")

    # ------------------------------------------------------------------ 15. FX2c: eta_prep is a probability; HOM relation
    def bloch_p0(shape, fwhm_ps, area_pi, G, k, gs):
        """INDEPENDENT no-jump Bloch ODE (own code): the radiative return G P -> rho_GG is dropped, the nonradiative
        k P -> rho_GG stays. Q = rho_GG(end) is p0, so eta_prep = 1 - Q.
        dP = -(G+k) P + Om y, dQ = k P - Om y, dy = -(G+k+gs)/2 y + (Om/2)(Q - P)."""
        dur = fwhm_ps * 1e-3
        if shape == "rect":
            win, Om = dur, (lambda t: area_pi * math.pi / dur)
        else:
            sg = dur / (2.0 * math.sqrt(math.log(2.0)))
            win = 2.0 * rs.GAUSS_WINDOW_SIGMA * sg
            Om = lambda t: area_pi * math.pi / (sg * math.sqrt(2 * math.pi)) * math.exp(-0.5 * ((t - win / 2) / sg) ** 2)
        gc = 0.5 * (G + k + gs)
        # state order (P, Q, y): dQ = k P - Om y
        f = lambda t, y: [-(G + k) * y[0] + Om(t) * y[2], k * y[0] - Om(t) * y[2],
                          -gc * y[2] + 0.5 * Om(t) * (y[1] - y[0])]
        sol = solve_ivp(f, [0.0, win], [0.0, 1.0, 0.0], method="DOP853", rtol=1e-12, atol=1e-14, max_step=win / 200)
        return 1.0 - float(sol.y[1, -1])

    ok, det = True, ""
    for shape, fw, G, k, gs in (("rect", 100.0, 10.0, 0.5, 3.0), ("rect", 60.0, 4.0, 0.0, 0.3), ("gauss", 20.0, 8.0, 0.2, 1.0)):
        for A in (0.6, 1.0, 2.2):
            r = LS.run(lind_plan(shape, fw, A), em0, rtx(G, gs, k))
            ref = bloch_p0(shape, fw, A, G, k, gs)
            tol = 2e-6 if shape == "rect" else 1e-4       # gauss: 8 steps/sigma midpoint grid error (6.4e-5 at 2.2 pi)
            ok &= near(r.extra["eta_prep"], ref, ab=tol)
            if shape == "gauss":                          # the same code on a 1600-step grid converges onto the ODE
                tl_f = rs.TwoLevel(G, k_ns=k, gamma_star_ns=gs, delta_ns=0.0)
                fine = 1.0 - solvers._window_p0(rs._pieces(tl_f), shape, A * math.pi, fw * 1e-3, (0,), n_steps=1600)
                ok &= near(fine, ref, ab=3e-6)
        det += f" {shape}{fw:g}ps A=2.2pi: {r.extra['eta_prep']:.6f}/{ref:.6f};"
    check("15a. FX2c lindblad eta_prep = P(n >= 1) = 1 - p0 equals an independent no-jump Bloch ODE (rect + gauss, three "
          "areas, decay, dephasing, k_nr) to 2e-6 (rect) / 1e-4 (gauss default grid; 3e-6 on a 1600-step grid)", ok, det)

    ok, det, worst = True, "", 0.0
    for shape in ("rect", "gauss"):
        for gt in (0.02, 0.1, 0.3):                       # Gamma tau_p
            for A in np.linspace(0.25, 3.0, 12):          # (area 0: pulsed HOM is undefined, see 15b0)
                fw = gt / 4.0 * 1e3
                r = LS.run(lind_plan(shape, fw, float(A)), em0, rtx(4.0))
                e = r.extra["eta_prep"]
                ok &= -1e-12 <= e <= 1.0 + 1e-12
                worst = max(worst, e)
    p0_zero = solvers._window_p0(rs._pieces(rs.TwoLevel(4.0)), "rect", 0.0, 0.05, (0,))
    check("15b. FX2c property: 0 <= eta_prep <= 1 for rect and gauss, areas 0.25..3 pi, Gamma tau_p 0.02..0.3, and "
          "eta_prep = 0 at area 0", ok and abs(p0_zero - 1.0) < 1e-12, f"max eta_prep {worst:.6f}")

    ok, det = True, ""
    for A in (0.3, 0.7, 1.0, 1.5):
        gt = 2.0 * 0.5e-3
        r = LS.run(lind_plan("rect", 0.5, A), em0, rtx(2.0))
        want = math.sin(0.5 * A * math.pi) ** 2
        ok &= abs(r.extra["eta_prep"] - want) <= gt
        det += f" {A}pi: {r.extra['eta_prep']:.5f}/{want:.5f};"
    check("15c. FX2c anchor: short pulse eta_prep = sin^2(area/2) within Gamma tau_p (area theorem)", ok, det)

    # 15d. the counting hierarchy gives the same probability as m1 - m2/2 when n <= 2 (state the truncation)
    from scipy.linalg import expm
    def hier_prep(shape, fw, A, G, gs):
        tl = rs.TwoLevel(G, k_ns=0.0, gamma_star_ns=gs, delta_ns=0.0)
        pc = rs._pieces(tl)
        dur = fw * 1e-3
        edges, om = rs._grid(shape, A * math.pi, dur, None, pc.max_dt)
        dt = np.diff(edges)
        V = np.zeros(12, dtype=complex)
        V[0] = 1.0
        for o, d in zip(om, dt):
            V = rs.expm(rs._aug(pc.L0 + o * pc.L1, pc.J) * d) @ V if hasattr(rs, "expm") else \
                __import__("scipy.linalg", fromlist=["expm"]).expm(rs._aug(pc.L0 + o * pc.L1, pc.J) * d) @ V
        tr = lb.trace_row(2)
        PX = np.zeros((2, 2), dtype=complex)
        PX[1, 1] = 1.0
        V1, V2 = V[4:8], V[8:]
        rho = lb.unvec(V[:4].copy(), 2)
        m1 = (tr @ V1).real + rho[1, 1].real                        # n = N_window + 1_X(end)
        EN_X = (lb.vec(PX.T) @ V1).real                              # E[N_window 1_X(end)]
        m2 = (tr @ V2).real + 2.0 * EN_X
        return float(m1 - 0.5 * m2)
    ok, det = True, ""
    for shape, fw, G, gs, A in (("rect", 60.0, 4.0, 0.3, 1.0), ("gauss", 20.0, 8.0, 1.0, 2.2), ("rect", 75.0, 4.0, 0.0, 3.0)):
        r = LS.run(lind_plan(shape, fw, A), em0, rtx(G, gs, 0.0))
        h = hier_prep(shape, fw, A, G, gs)
        gap = r.extra["eta_prep"] - h                    # = P(n=3) + 3 P(n=4) + ... >= 0
        ok &= -1e-9 <= gap <= 1e-2 and 0.0 <= h <= 1.0 + 1e-9
        det += f" {shape}{fw:g}ps A={A}pi G tau_p={G * fw * 1e-3:.2f}: 1-p0 - (m1-m2/2) = {gap:.1e};"
    check("15d. FX2c m1 - m2/2 (window moments + end population) equals P(n >= 1) - P(n=3) - 3 P(n=4) - ...: it never "
          "exceeds 1 - p0 and is within 1e-2 of it up to Gamma tau_p 0.3, 3 pi (the truncation at n <= 2); both in [0, 1]",
          ok, det)

    # 15e. HOM relation: exact two-source coincidence algebra on an independent Fock-space model.
    # Source state sqrt(1-q)|0,0> + sqrt(q)|1_u,1_v> with a random optical phase between independent sources (the
    # two-instantaneous-pi-pulse emitter: photon u before the
    # second pulse, photon v after it; n = 0 or 2, g2 = 1/(2q), whole-stream overlap I = 1/2 while every photon
    # is a perfect wavepacket, M_s = 1).
    def fock_hom(q, cut=5):
        a1 = np.diag(np.sqrt(np.arange(1, cut)), 1)
        eye = np.eye(cut)
        def op(k, A):
            mats = [eye] * 4
            mats[k] = A
            out = mats[0]
            for m in mats[1:]:
                out = np.kron(out, m)
            return out
        uA, vA, uB, vB = (op(k, a1) for k in range(4))
        nC = 0.5 * ((uA + uB).conj().T @ (uA + uB) + (vA + vB).conj().T @ (vA + vB))
        nD = 0.5 * ((uA - uB).conj().T @ (uA - uB) + (vA - vB).conj().T @ (vA - vB))
        vac = np.zeros(cut ** 4)
        vac[0] = 1.0
        # independent sources with random optical phase: the photon-number-dephased product state, i.e. a mixture
        # of |0>, |uA vA>, |uB vB>, |uA vA uB vB> with weights (1-q)^2, q(1-q), q(1-q), q^2
        creA, creB = uA.conj().T @ vA.conj().T, uB.conj().T @ vB.conj().T
        comps = [((1 - q) ** 2, vac), (q * (1 - q), creA @ vac), (q * (1 - q), creB @ vac), (q * q, creA @ (creB @ vac))]
        C_par = sum(w * float(c @ (nC @ (nD @ c))) for w, c in comps)
        Nc = sum(w * float(c @ (nC @ c)) for w, c in comps)
        Nd = sum(w * float(c @ (nD @ c)) for w, c in comps)
        return C_par, Nc, Nd
    ok, det = True, ""
    for q in (0.4, 0.78, 0.95):
        C_par, Nc, Nd = fock_hom(q)
        N = 2.0 * q                                     # mean photons per source
        C_far = Nc * Nd                                 # uncorrelated pulses: product of the singles
        V_exact = 1.0 - 2.0 * C_par / C_far
        g2q, Iq = 1.0 / (2.0 * q), 0.5
        ok &= near(Nc, N, ab=1e-12) and near(Nd, N, ab=1e-12)
        ok &= near(V_exact, OB.hom_v_raw(Iq, g2q), ab=1e-12) and near(V_exact, (q - 1.0) / (2.0 * q), ab=1e-12)
        ok &= abs(((1 - g2q) * Iq - g2q) - V_exact) > 0.02          # the pre-FX2c relation is off
        det += f" q={q}: V {V_exact:.6f} = I - g2 {Iq - g2q:.6f} (old {(1 - g2q) * Iq - g2q:.6f});"
    check("15e. FX2c HOM: exact beam-splitter coincidences of two identical sources (Fock space, n = 0 or 2 photons) give "
          "V = 1 - 2 C_par/C_far = I - g2 for I = Int Int |G1|^2/N^2 = 1/2, g2 = 1/(2q); the old ((1 - g2) I - g2) is "
          "off by > 0.02 (it would be right for the single-photon overlap M_s = 1: 1 - 2 g2)", ok, det)

    sp_h = engine._as_spec(spec())
    ent_h = engine._entries(sp_h)
    c_h = engine._run_point(ent_h, engine._central_values(engine._leaves(sp_h, ent_h)), 4.2, "auto")
    dh = c_h.d
    want_h = (dh["hom_I"] * (1.0 - dh["g2_raw"]) / (1.0 - dh["g2_0"]) - dh["g2_raw"])
    check("15f. engine hom_V_raw = (rho^2 hom_I - g2_raw) with rho^2 = (1 - g2_raw)/(1 - g2_0); hom_V_corr = hom_I/(1 - g2_0)",
          math.isfinite(dh["hom_V_raw"]) and near(dh["hom_V_raw"], want_h, ab=1e-12)
          and near(dh["hom_V_corr"], dh["hom_I"] / (1.0 - dh["g2_0"]), ab=1e-9),
          f"hom_I {dh['hom_I']:.5f}, g2_0 {dh['g2_0']:.5f}, V_raw {dh['hom_V_raw']:.5f}, V_corr {dh['hom_V_corr']:.5f}")

    # 16. devlib-fx3b (real registry entries): M-3 eta_capture double count, bg_corrected g2, cw brightness
    E_Q = 1.602176634e-19

    def rspec(kind, drp, mat=None, su=None):
        return {"schema": "fsim-composed/1", "T_K": 4.2,
                "material": {"kind": "inas_gaas_sk", "params": {"tau_bare_ns": 0.599, "gamma_star_4K_ns": 0.0,
                                                               "a_ac_ns_per_K": 0.0, **(mat or {})}},
                "structure": {"kind": "bulk_planar", "params": {}}, "drive": {"kind": kind, "params": drp},
                "setup": {"kind": "lens_fibre_det", "params": su or {}}, "solver": "rate"}

    ok16, det16 = True, ""
    mu3 = 3.0
    I3 = mu3 * E_Q / 0.2e-9 * 1e6                 # uA for mu = I tau / e = 3 over 200 ps
    for prof in ("rect", "gauss", "rc"):
        runs = {}
        for eta in (1.0, 0.4, 0.05):
            drp = {"eta_capture": eta, "tau_on_ns": 0.2, "I_uA": I3 / eta, "f_rep_MHz": 200.0, "pump_profile": prof,
                   "t_rise_ps": 100.0}
            runs[eta] = engine.evaluate_spec(rspec("electrical_pulsed", drp))["scalars"]
        # mu = eta N_c = I tau/e = 3 for every eta: loaded probability 1 - exp(-3) (independent closed form, Poisson loading)
        ok16 &= all(near(r_["eta_prep"], 1.0 - math.exp(-mu3), ab=1e-12) for r_ in runs.values())
        ok16 &= all(near(runs[e_]["n_ph"], runs[1.0]["n_ph"], rel=1e-12) and near(runs[e_]["g2_0"], runs[1.0]["g2_0"], rel=1e-12)
                    for e_ in (0.4, 0.05))
        det16 += f" {prof}: p_loaded {runs[0.4]['eta_prep']:.6f};"
    check("16a. M-3 pulsed: eta_capture is applied ONCE (mu = eta N_c = I tau/e): with I scaled by 1/eta the loaded probability is "
          "1 - e^{-mu} = 0.9502 and n_ph, g2 are identical for eta = 1, 0.4, 0.05 (rect, Gaussian, RC profiles; the old double "
          "factor gave 1 - e^{-eta mu})", ok16, det16)
    mu_a = engine.evaluate_spec(rspec("electrical_pulsed", {"eta_capture": 0.4, "tau_on_ns": 0.2, "I_uA": I3, "f_rep_MHz": 200.0}))["scalars"]
    check("16a'. and eta_capture does scale mu: the same current at eta = 0.4 loads 1 - e^{-0.4 mu}",
          near(mu_a["eta_prep"], 1.0 - math.exp(-0.4 * mu3), ab=1e-12), f"{mu_a['eta_prep']:.6f}")
    # R2 D3: every mechanism's mu is the pairs loaded INTO THE DOT; eta_capture does not rescale it (cw r = mu (Gamma + k_nr))
    cw_r = engine.evaluate_spec(rspec("electrical_cw", {"mechanism": "M-1 poisson-rail", "mu": 0.5, "eta_capture": 0.2}))["scalars"]
    G16 = float(cw_r["Gamma_ns"])
    n_x = 0.5 / (0.5 + 1.0)                      # r = mu Gamma = 0.5 Gamma (k_nr = 0 at 4 K): n_X = r/(r + Gamma)
    check("16a''. M-1 cw r = mu Gamma, eta_capture = 0.2 does NOT rescale mu (R2 D3): flux_ns/Gamma = r/(r + Gamma) "
          "(independent two-level steady state)", near(float(cw_r["flux_ns"]) / G16, n_x, rel=1e-9),
          f"{float(cw_r['flux_ns']) / G16:.6f} vs {n_x:.6f}")

    # bg_corrected convention
    e_b = engine.envelope(rspec("electrical_cw", {"bg_frac": 0.08}, su={"laser_leak": 0.10}),
                          ["g2_0", "g2_0:raw", "g2_0:bg_corrected"])
    g_sig16 = 0.0                                 # an ideal cw two-level emitter: g2(0) = 0
    g_emb = 1.0 - 0.92 ** 2 * (1.0 - g_sig16)     # drive-embedded background
    g_raw16 = 1.0 - (1.0 / 1.10) ** 2 * (1.0 - g_emb)       # + laser leak 0.10 of the signal
    check("16b. g2_0:bg_corrected = (g2_raw - (1 - rho^2))/rho^2 with rho = (1 - bg_frac)/(1 + laser_leak) is the exact inverse of "
          "the background law: ideal emitter, bg_frac 0.08 + laser leak 0.10 -> intrinsic 0.1536, raw "
          f"{g_raw16:.4f}, bg_corrected 0 (1e-12)",
          near(e_b["g2_0"]["central"], g_emb, ab=1e-12) and near(e_b["g2_0:raw"]["central"], g_raw16, ab=1e-12)
          and abs(e_b["g2_0:bg_corrected"]["central"] - g_sig16) < 1e-12)
    e_d = engine.envelope(rspec("electrical_cw", {"bg_frac": 0.08}, su={"dark_cps": 5000.0}), ["g2_0:raw", "g2_0:bg_corrected"])
    check("16b'. dark counts are NOT removed by the background correction (the dark-count-corrected value is a separate "
          "convention): with darks the corrected value stays above the signal g2 = 0 and below the raw value",
          0.0 < e_d["g2_0:bg_corrected"]["central"] < e_d["g2_0:raw"]["central"],
          f"raw {e_d['g2_0:raw']['central']:.4f}, bgc {e_d['g2_0:bg_corrected']['central']:.4f}")
    e_r = engine.envelope(rspec("electrical_cw", {}), ["g2_0:bg_corrected", "g2_0"])
    check("16b''. no background at all: bg_corrected == intrinsic", near(e_r["g2_0:bg_corrected"]["central"], e_r["g2_0"]["central"], ab=1e-12))

    # cw brightness
    def chain(sc_):
        return float(sc_["eta_extract"]) * float(sc_["eta_pol"]) * float(sc_["Z_filt"])
    flux_hz = G16 * n_x * 1e9
    e_c = engine.envelope(rspec("electrical_cw", {"mechanism": "M-1 poisson-rail", "mu": 0.5, "eta_capture": 0.2}),
                          ["B_first_lens_Hz", "B_fibre_Hz", "rate_detected_Hz"])
    want_lens = flux_hz * chain(cw_r)
    check("16c. cw brightness (electrical, rate solver): photon rate into the first lens = Gamma n_X,ss beta eta_out eta_pol (ns -> Hz), "
          "into the fibre x eta_fibre, detected x eta_optics eta_det = 0.25 (setup defaults); n_X,ss = r/(r + Gamma) is the "
          "independent two-level steady state",
          near(e_c["B_first_lens_Hz"]["central"], want_lens, rel=1e-9) and near(e_c["B_fibre_Hz"]["central"], want_lens, rel=1e-9)
          and near(e_c["rate_detected_Hz"]["central"], want_lens * 0.25, rel=1e-9),
          f"{e_c['B_first_lens_Hz']['central']:.4g} Hz vs {want_lens:.4g}")
    # cw resonance fluorescence on resonance: pump rate R = Omega^2/(2 gamma_2) (gamma_2 = coherence decay), rho_ee = R/(Gamma_t + 2R)
    # = (s/2)/(1 + s) with s = Omega^2 T1 T2 = Omega^2/(gamma_2 Gamma_t) (Loudon: pure radiative s = 2 Omega^2/Gamma^2)
    for s_ in (0.5, 2.0, 10.0):
        sp_f = {**rspec("cw_resonant", {"saturation_s": s_}, mat={"gamma_star_4K_ns": 0.7}), "solver": "auto"}
        e_f = engine.envelope(sp_f, ["B_first_lens_Hz", "rate_detected_Hz"])
        sc_f = engine.evaluate_spec(sp_f)["scalars"]
        rho_ee = (s_ / 2.0) / (1.0 + s_)
        want = float(sc_f["Gamma_ns"]) * rho_ee * 1e9 * chain(sc_f)
        check(f"16d. cw resonant, s = {s_}: first-lens photon rate = Gamma rho_ee chain with rho_ee = (s/2)/(1 + s) (Loudon)",
              near(e_f["B_first_lens_Hz"]["central"], want, rel=2e-4) and near(e_f["rate_detected_Hz"]["central"], want * 0.25, rel=2e-4),
              f"{e_f['B_first_lens_Hz']['central']:.5g} vs {want:.5g}")
    e_p = engine.envelope(rspec("electrical_pulsed", {"I_uA": I3, "tau_on_ns": 0.2, "eta_capture": 1.0, "f_rep_MHz": 200.0}),
                          ["B_first_lens", "B_first_lens_Hz", "B_fibre", "B_fibre_Hz"])
    check("16e. pulsed: B_first_lens_Hz = f_rep B_first_lens and B_fibre_Hz = f_rep B_fibre",
          near(e_p["B_first_lens_Hz"]["central"], 200e6 * e_p["B_first_lens"]["central"], rel=1e-12)
          and near(e_p["B_fibre_Hz"]["central"], 200e6 * e_p["B_fibre"]["central"], rel=1e-12))
    check("16f. the cw drives offer B_first_lens_Hz / B_fibre_Hz / rate_detected_Hz",
          all(o in R.get("drive", k).observables for k in ("electrical_cw", "cw_resonant", "above_band") for o in ("B_first_lens_Hz", "B_fibre_Hz", "rate_detected_Hz")))

    # ---- FX4 (final review 1, 2, 6): real paper devices through the real registry
    import copy as _copy
    import yaml as _yaml

    def card_spec(cid, drive_extra=None):
        raw = _yaml.safe_load((ROOT / "cards" / "papers" / f"{cid}.yaml").read_text(encoding="utf-8"))
        d = _copy.deepcopy(raw["device"])
        for k, v in (drive_extra or {}).items():
            d["drive"]["params"][k] = {"value": v, "unit": "1" if k == "bg_frac" or k == "p_trion" else "ps", "tag": "A",
                                      "source": "FX4 check"}
        return {"schema": "fsim-composed/1", **d}

    I0 = None
    vc, ok_bg = [], True
    for bg_ in (0.0, 0.1, 0.3, 0.5):
        sc_b = engine.evaluate_spec(card_spec("somaschi2016_qd1", {"bg_frac": bg_}))["scalars"]
        I0 = sc_b["hom_I"] if I0 is None else I0
        vc.append(sc_b["hom_V_corr"])
        ok_bg &= near(sc_b["hom_I"], I0, rel=1e-12)
        # independent: V_raw = rho^2 I - g2_raw with rho = 1 - bg (no leak / dark counts on this card)
        want_raw = (1.0 - bg_) ** 2 * sc_b["hom_I"] - sc_b["g2_raw"]
        ok_bg &= near(sc_b["hom_V_raw"], want_raw, ab=1e-9)
    check("FX4-1. HOM with bg_frac: V_corr = I/(1 - g2_signal) is invariant under bg_frac 0..0.5 and V_raw = rho^2 I - g2_raw",
          ok_bg and max(vc) - min(vc) < 1e-9, f"V_corr {[round(v, 6) for v in vc]}")
    sc_t0 = engine.evaluate_spec(card_spec("tomm2021"))["scalars"]
    sc_t1 = engine.evaluate_spec(card_spec("tomm2021", {"p_trion": 1e-6}))["scalars"]
    check("FX4-2. JC regime with p_trion 1e-6: still the JC solver, valid, tau within 1e-4 (relative) of p_trion 0 "
          "(trion = post-processing on the JC result)",
          sc_t1["solver"] == "jc" and sc_t1["valid"] and near(sc_t1["tau_ps"], sc_t0["tau_ps"], rel=1e-4),
          f"{sc_t0['tau_ps']:.5f} vs {sc_t1['tau_ps']:.5f}")
    sc_e = engine.evaluate_spec(card_spec("tomm2021", {"eid_ps": 1.0}))["scalars"]
    check("FX4-2b. JC regime with EID: result invalid with the reason 'JC regime: EID not supported by the JC solver'",
          (not sc_e["valid"]) and any("JC regime: EID not supported by the JC solver" in r for r in sc_e["invalid_reasons"]),
          str(sc_e["invalid_reasons"]))
    sc_w = engine.evaluate_spec(card_spec("wang2019b"))["scalars"]
    check("FX4-6. corrected HOM visibility > 1 -> NaN with the reason 'corrected visibility > 1: ...' (never a pass)",
          not math.isfinite(sc_w["hom_V_corr"]) and "corrected visibility > 1" in sc_w.get("hom_V_corr_invalid_reason", "")
          and math.isfinite(sc_w["hom_V_raw"]))

    # 17. R2: cw incoherent pump P broadens the dot line the cavity filter sees (ch08 eq. dp-lor): w_dot = hbar (Gamma_bg + k + gamma* + P)
    def cav_cw(pump):
        return {"schema": "fsim-composed/1", "T_K": 4.2,
                "material": {"kind": "inas_gaas_sk", "params": {"tau_bare_ns": 0.599, "gamma_star_4K_ns": 0.0,
                                                               "a_ac_ns_per_K": 0.0}},
                "structure": {"kind": "generic_cavity", "params": {"F_P": 5.0, "kappa_ueV": 20.0, "detuning_ueV": 0.0}},
                "drive": {"kind": "above_band", "params": {"cw": True, "pump": pump}},
                "setup": {"kind": "lens_fibre_det", "params": {}}, "solver": "rate"}
    sc0 = engine.evaluate_spec(cav_cw(1e-6))["scalars"]
    sc1 = engine.evaluate_spec(cav_cw(3.0))["scalars"]
    G0_ = 1.0 / 0.599
    hb = lb.HBAR_MEV_NS * 1e3
    kap_ = 20.0
    L0_ = kap_ / (kap_ + hb * G0_)                         # zero-detuning Lorentzian overlap kappa/(kappa + w) [DR]
    Gam0_ = G0_ * (1.0 + 5.0 * L0_)
    r_ = 3.0 * Gam0_                                       # planner: r = P/P_sat (Gamma + k) with the UNBROADENED Gamma
    L1_ = kap_ / (kap_ + hb * (G0_ + r_))
    Gam1_ = G0_ * (1.0 + 5.0 * L1_)
    want_flux = Gam1_ * r_ / (r_ + Gam1_)
    check("17a. cw above_band, cavity: flux = Gamma_P r/(r + Gamma_P) with Gamma_P = Gamma0 (1 + F_P kappa/(kappa + hbar(Gamma0 + r))) "
          "(independent closed form, 1e-6)", abs(float(sc1["flux_ns"]) - want_flux) < 1e-6 * want_flux,
          f"{float(sc1['flux_ns']):.6f} vs {want_flux:.6f}")
    check("17b. the pump broadening lowers Gamma seen by the solver: Gamma_P < Gamma_unbroadened (3 P_sat)", Gam1_ < Gam0_ - 1e-6, f"{Gam1_:.4f} < {Gam0_:.4f}")

    # ------------------------------------------------------------------ 18. R2 FIX-1a
    from scipy.signal import find_peaks
    # 18a. combine(..., pump_ns) hook: pump_ns = 0 is rates.combine bit for bit; the observables see the PUMPED Gamma / beta
    em18 = R.get("material", "inas_gaas_sk").build({}, 4.2)
    sr18 = R.get("structure", "generic_cavity").response({"F_P": 5.0, "kappa_ueV": 20.0, "detuning_ueV": 0.0}, em18, 4.2, 0.65)
    check("18a. engine.combine(pump_ns=0) == rates.combine bit for bit",
          engine.combine(em18, sr18, 4.2) == rates.combine(em18, sr18, 4.2))
    sc_p = engine.evaluate_spec(cav_cw(3.0))["scalars"]
    check("18b. cw pump hook: scalars Gamma_ns / tau_ps / F_P_eff are the PUMPED values (closed form of 17, 1e-9), not the unpumped Gamma0 (1 + F_P L0)",
          near(sc_p["Gamma_ns"], Gam1_, rel=1e-9) and near(sc_p["tau_ps"], 1e3 / Gam1_, rel=1e-9)
          and near(sc_p["F_P_eff"], 5.0 * L1_, rel=1e-9) and sc_p["Gamma_ns"] < Gam0_ - 1e-6,
          f"Gamma {sc_p['Gamma_ns']:.5f} vs {Gam1_:.5f} (unpumped {Gam0_:.5f}); L {sc_p['F_P_eff'] / 5.0:.5f} vs {L1_:.5f}")

    # 18c. tau_ratio = tau_0/tau: radiative convention = F_bg + F_P L (= Gamma/Gamma0), intrinsic = tau_0/tau (bulk 1/(Gamma0 + k))
    mk = {"gamma0_ns": exact(1.0, "1/ns"), "k_nr_ns": exact(0.5, "1/ns")}
    def ctx_of(sp_):
        sp_ = engine._as_spec(sp_)
        ent_ = engine._entries(sp_)
        return engine._run_point(ent_, engine._central_values(engine._leaves(sp_, ent_)), 4.2, sp_.solver)
    ctx18 = ctx_of(spec(mat=mk, st={"F_P": exact(5.0, "1", "E")}, drp={"fwhm_ps": 1.0}))
    t_rad = OB.get_observable("tau_ratio").compute(ctx18, "radiative")
    t_int = OB.get_observable("tau_ratio").compute(ctx18, "intrinsic")
    check("18c. tau_ratio:radiative = F_bg + F_P L = 1 + 5 = 6 and tau_ratio = tau_0/tau = (1/(1 + 0.5))/(1/(6 + 0.5)) = 13/3 "
          "(closed forms; equal for k = 0)", near(t_rad, 6.0, rel=1e-12) and near(t_int, 13.0 / 3.0, rel=1e-9),
          f"{t_rad!r} {t_int!r}")
    rj = engine.evaluate_spec(spec(mat={"gamma0_ns": exact(G0j, "1/ns")}, st={"g_ns": exact(g, "1/ns"), "kappa_ns": exact(kp, "1/ns")},
                                   drp={"fwhm_ps": 1.0}, solver="jc"))["scalars"]
    check("18d. tau_ratio with the jc solver is tau_0/tau_jc = 1 + 4 g^2/(kappa Gamma0) = 1.5 (bad cavity, 2%)",
          near(rj["tau_ratio"], 1.0 + 4 * g * g / (kp * G0j), rel=0.02), f"{rj['tau_ratio']:.4f}")

    # 18e. s_0 / s for cw_resonant: s = Omega^2 T1 T2, s_0 = 2 Omega^2/Gamma^2 (MvdS); equal at gamma* = 0, s_0/s = (Gamma + k + gamma*)/(Gamma + k)
    def cwres(gs_, s_):
        return {"schema": "fsim-composed/1", "T_K": 4.2,
                "material": {"kind": "inas_gaas_sk", "params": {"tau_bare_ns": 0.8, "gamma_star_4K_ns": exact(gs_, "1/ns"),
                                                               "a_ac_ns_per_K": 0.0}},
                "structure": {"kind": "bulk_planar", "params": {}},
                "drive": {"kind": "cw_resonant", "params": {"saturation_s": exact(s_, "1")}},
                "setup": {"kind": "lens_fibre_det", "params": {}}, "solver": "auto"}
    a0 = engine.evaluate_spec(cwres(0.0, 3.0))["scalars"]
    a1 = engine.evaluate_spec(cwres(1.2, 3.0))["scalars"]
    check("18e. cw_resonant: s_0 = 2 Omega^2/Gamma^2 equals the entered s = Omega^2 T1 T2 at gamma* = 0 (1e-9), and s_sat = s always",
          near(a0["s_0"], 3.0, rel=1e-9) and near(a0["s_sat"], 3.0, rel=1e-9) and near(a1["s_sat"], 3.0, rel=1e-9),
          f"{a0['s_0']!r} {a0['s_sat']!r} {a1['s_sat']!r}")
    Gq = a1["Gamma_ns"] + 0.0
    check("18f. cw_resonant with dephasing: s_0/s = (Gamma + k + Gamma_ph)/(Gamma + k) (closed form, 1e-9); the elastic fraction is reported",
          near(a1["s_0"] / a1["s_sat"], (Gq + 1.2) / Gq, rel=1e-9) and 0.0 <= a1["coherent_fraction"] <= 1.0,
          f"{a1['s_0'] / a1['s_sat']!r} vs {(Gq + 1.2) / Gq!r}")

    # 18g. strong-coupling jc paths -------------------------------------------------------------------------------------
    R.register(TDrive(kind="test_incoh_cw", label="t", summary="t", flavor="incoh_cw", params=(P("r_ns", "1/ns", 1.2),)))
    jc_lim = {"g_ns": exact(0.05, "1/ns"), "kappa_ns": exact(100.0, "1/ns"), "beta0": exact(0.0, "1")}      # weak-coupling limit
    # (i) cw resonant: jc -> the bare-emitter Lindblad result as g -> 0 (Mollow flux, inversion, g2(0), elastic fraction)
    mat18 = {"gamma0_ns": exact(2.0, "1/ns")}
    sj = engine.evaluate_spec(spec(mat=mat18, st=jc_lim, dr="test_cw", drp={"omega_ns": exact(6.0, "1/ns")}, solver="jc"))
    sl = engine.evaluate_spec(spec(mat=mat18, st=jc_lim, dr="test_cw", drp={"omega_ns": exact(6.0, "1/ns")}, solver="lindblad"))
    cj, cl = sj["scalars"], sl["scalars"]
    check("18g. jc cw resonant, g -> 0: flux, inversion, elastic fraction equal the bare-emitter Lindblad (1e-3) and g2(0) = 0 (1e-3)",
          cj["solver"] == "jc" and near(cj["flux_ns"], cl["flux_ns"], rel=1e-3) and near(cj["inversion"], cl["inversion"], rel=1e-3)
          and near(cj["coherent_fraction"], cl["coherent_fraction"], rel=1e-3) and abs(cj["g2_0"]) < 1e-3,
          f"flux {cj['flux_ns']:.5f}/{cl['flux_ns']:.5f} inv {cj['inversion']:.5f}/{cl['inversion']:.5f} g2 {cj['g2_0']:.2e}")
    # (ii) pulsed incoherent: jc -> lindblad_incoherent (instantaneous loading) as g -> 0: g2, photons per pulse, HOM
    mat_d = {"gamma0_ns": exact(1.0, "1/ns"), "gamma_star_ns": exact(0.4, "1/ns")}
    ip = {"relax_ns": exact(0.0, "ns")}
    ij = engine.evaluate_spec(spec(mat=mat_d, st=jc_lim, dr="test_incoh", drp=ip, solver="jc"))["scalars"]
    il = engine.evaluate_spec(spec(mat=mat_d, st=jc_lim, dr="test_incoh", drp=ip, solver="lindblad_incoherent"))["scalars"]
    check("18h. jc pulsed incoherent, g -> 0: g2 and photons per pulse equal lindblad_incoherent (1e-3), cavity-output HOM = Gamma/(Gamma + gamma*) (2%)",
          ij["solver"] == "jc" and near(ij["g2_0"], il["g2_0"], rel=1e-3, ab=1e-6) and near(ij["B_lens"], il["B_lens"], rel=2e-3)
          and near(ij["hom_I"], il["hom_I"], rel=0.02),
          f"g2 {ij['g2_0']:.5f}/{il['g2_0']:.5f} B {ij['B_lens']:.5f}/{il['B_lens']:.5f} hom {ij['hom_I']:.4f}/{il['hom_I']:.4f}")
    # (iii) cw incoherent: jc -> the rate solver (flux = Gamma r/(r + Gamma), g2(0) = 0)
    ccj = engine.evaluate_spec(spec(mat=mat18, st=jc_lim, dr="test_incoh_cw", drp={"r_ns": exact(2.0, "1/ns")}, solver="jc"))["scalars"]
    ccr = engine.evaluate_spec(spec(mat=mat18, st=jc_lim, dr="test_incoh_cw", drp={"r_ns": exact(2.0, "1/ns")}, solver="rate"))["scalars"]
    check("18i. jc cw incoherent, g -> 0: flux equals the rate solver and Gamma r/(r + Gamma) (1e-3), g2(0) = 0",
          near(ccj["flux_ns"], ccr["flux_ns"], rel=1e-3) and near(ccj["flux_ns"], 2.0 * 2.0 / 4.0, rel=2e-3) and abs(ccj["g2_0"]) < 1e-3,
          f"{ccj['flux_ns']:.5f} {ccr['flux_ns']:.5f} (closed 1.0)")
    # (iv) Purcell rate 4 g^2/kappa in the bad-cavity limit: cw incoherent flux at r = Gamma_P is Gamma_P/2, the cavity share is F/(1 + F)
    bg_, bk_ = 5.0, 500.0
    FPb = 4 * bg_ ** 2 / (bk_ * 1.0)
    GP = 1.0 + 4 * bg_ ** 2 / bk_
    pu = engine.evaluate_spec(spec(st={"g_ns": exact(bg_, "1/ns"), "kappa_ns": exact(bk_, "1/ns"), "beta0": exact(0.0, "1")},
                                   dr="test_incoh_cw", drp={"r_ns": exact(GP, "1/ns")}, solver="jc"))["scalars"]
    check("18j. jc cw incoherent, bad cavity (kappa 500, g 5): flux at r = Gamma_P is Gamma_P/2 (1%) and beta = F_P/(1 + F_P) = 0.1667 (1%) "
          "[Purcell rate 4 g^2/kappa]", near(pu["flux_ns"], GP / 2, rel=0.01) and near(pu["beta"], FPb / (1 + FPb), rel=0.01),
          f"flux {pu['flux_ns']:.5f} vs {GP / 2:.5f}; beta {pu['beta']:.5f} vs {FPb / (1 + FPb):.5f}")
    # (v) vacuum Rabi splitting 2g in the cw emission spectrum (weak incoherent pump, g >> kappa)
    vr = engine.evaluate_spec(spec(mat={"gamma0_ns": exact(0.02, "1/ns")}, st={"g_ns": exact(5.0, "1/ns"), "kappa_ns": exact(0.5, "1/ns")},
                                   dr="test_incoh_cw", drp={"r_ns": exact(0.01, "1/ns")}, solver="jc"))
    sp_ = vr["series"]["spectrum_cav"]
    x_, y_ = np.array(sp_["x"]), np.array(sp_["y"])
    pk, _ = find_peaks(y_)
    top = sorted(pk, key=lambda i: -y_[i])[:2]
    split = abs(x_[top[0]] - x_[top[1]]) if len(top) == 2 else float("nan")
    check("18k. vacuum Rabi splitting: the two spectral peaks are 2 g = 10/ns apart (1%)", near(split, 10.0, rel=0.01),
          f"split {split:.4f}/ns, peaks at {x_[top[0]]:.3f}, {x_[top[1]]:.3f}")

    # 18l. audit 4: EVERY drive on a strong-coupling phc_cavity is either valid on jc or invalid with a reason that says why, and the
    # reason never recommends the jc solver when jc refuses the run; when it does recommend it, the explicit jc run succeeds
    def sc_spec(drv, solver="auto"):
        return {"schema": "fsim-composed/1", "T_K": 4.2,
                "material": {"kind": "inas_gaas_sk", "params": {}},
                "structure": {"kind": "phc_cavity", "params": {}},
                "drive": {"kind": drv, "params": {}},
                "setup": {"kind": "lens_fibre_det", "params": {}}, "solver": solver}
    bad_all, n_jc, n_inv, det18 = [], 0, 0, []
    for dk in R.kinds("drive"):
        sc_ = engine.evaluate_spec(sc_spec(dk))["scalars"]
        assert sc_["regime"] == "strong", (dk, sc_["regime"])
        if sc_["valid"]:
            n_jc += 1
            if sc_["solver"] != "jc":
                bad_all.append(f"{dk}: valid but solver {sc_['solver']}")
            continue
        n_inv += 1
        why_ = " ".join(sc_["invalid_reasons"])
        recommends = "use solver 'jc'" in why_
        try:
            engine.evaluate_spec(sc_spec(dk, "jc"))
            jc_runs = True
        except Exception:                                   # noqa: BLE001  (SolverUnavailable: jc refuses this plan)
            jc_runs = False
        if recommends and not jc_runs:
            bad_all.append(f"{dk}: recommends jc but jc refuses")
        if (not recommends) and ("not supported" not in why_ or "Jaynes-Cummings" not in why_):
            bad_all.append(f"{dk}: invalid without a clear unsupported reason: {why_[:80]}")
        det18.append(dk)
    check("18l. phc_cavity (strong coupling, inas_gaas_sk): every drive is valid on jc or invalid with an explicit 'not supported' reason; "
          "no reason recommends jc where jc refuses (audit 4)", not bad_all and n_jc >= 6,
          f"jc-valid {n_jc}, invalid-with-reason {n_inv} {det18}; problems {bad_all}")
    sc_tpe = engine.evaluate_spec(sc_spec("tpe"))["scalars"]
    check("18m. TPE in the strong-coupling regime: invalid, the reason names the missing biexciton ladder and gives the remedy (reduce Q / raise V)",
          (not sc_tpe["valid"]) and "biexciton" in " ".join(sc_tpe["invalid_reasons"]) and "weaken the coupling" in " ".join(sc_tpe["invalid_reasons"]),
          " ".join(sc_tpe["invalid_reasons"])[:160])

    # 18n. captions follow the solver / timing (audit 17)
    cap_j = engine.evaluate_spec(sc_spec("resonant_pulse"))["scalars"]["captions"]
    cap_c = engine.evaluate_spec(sc_spec("cw_resonant"))["scalars"]["captions"]
    cap_w = engine.evaluate_spec(spec(drp={"fwhm_ps": 1.0}))["scalars"]["captions"]
    check("18n. captions: jc tau = slowest Liouvillian decay (not 1/(Gamma + k)); weak-coupling tau = 1/(Gamma + k_nr); cw rate has no f_rep; pulsed has",
          "Liouvillian" in cap_j["tau_ps"] and "NOT 1/(Gamma + k_nr)" in cap_j["tau_ps"] and cap_w["tau_ps"].startswith("tau = 1/(Gamma + k_nr)")
          and "f_rep" not in cap_c["rate_detected_Hz"].replace("no f_rep", "") and "f_rep . B_fibre" in cap_w["rate_detected_Hz"],
          f"{cap_j['tau_ps'][:50]} | {cap_c['rate_detected_Hz'][:50]}")

    # 18o. audit 6: cw drives do not advertise hom_I / hom_V, and a requested hom on a cw run is NaN with the reason
    cwk = [k for k in ("cw_resonant", "electrical_cw") if k in R.kinds("drive")]
    adv_ok = all("hom_I" not in R.get("drive", k).observables and "hom_V" not in R.get("drive", k).observables for k in cwk)
    sc_ec = engine.evaluate_spec({"schema": "fsim-composed/1", "T_K": 4.2,
                                  "material": {"kind": "inas_gaas_sk", "params": {}}, "structure": {"kind": "bulk_planar", "params": {}},
                                  "drive": {"kind": "electrical_cw", "params": {}}, "setup": {"kind": "lens_fibre_det", "params": {}},
                                  "solver": "auto"})["scalars"]
    check("18o. cw drives (cw_resonant, electrical_cw) advertise no hom_I / hom_V; a hom request on a cw run is NaN with a cw reason",
          adv_ok and not math.isfinite(sc_ec["hom_I"]) and "cw stream" in sc_ec.get("hom_I_invalid_reason", ""),
          sc_ec.get("hom_I_invalid_reason", ""))

    # 18p. the nitride RT example: g2 sits in the Tamariz 2020 class range (0.23 +- 0.05, GaN/AlN, 300 K [V]); HOM ~ 0 at 300 K
    import yaml
    ex = yaml.safe_load((ROOT / "cards" / "composed" / "example-nitride-rt.yaml").read_text(encoding="utf-8"))
    sc_ex = engine.evaluate_spec(ex["design"]["composed"])["scalars"]
    check("18p. example-nitride-rt: g2(0) within 0.23 +- 0.05 (Tamariz 2020 average at 300 K) and HOM < 0.01 (no longer 2.6e-4)",
          abs(sc_ex["g2_0"] - 0.23) <= 0.05 and sc_ex["hom_I"] < 0.01 and sc_ex["valid"], f"g2 {sc_ex['g2_0']:.4f}, hom {sc_ex['hom_I']:.2e}")

    # 18q. invariant: a valid pulsed run that prepares < 5 % carries a note (the SUPER default no longer does, a wrong area does)
    def super_spec(a1_, a2_):
        p = {"area1_pi": exact(a1_, "pi"), "area2_pi": exact(a2_, "pi")}
        return {"schema": "fsim-composed/1", "T_K": 4.2, "material": {"kind": "inas_gaas_sk", "params": {}},
                "structure": {"kind": "bulk_planar", "params": {}}, "drive": {"kind": "super", "params": p},
                "setup": {"kind": "lens_fibre_det", "params": {}}, "solver": "auto"}
    s_bad = engine.evaluate_spec(super_spec(2.0, 2.0))["scalars"]
    s_def = engine.evaluate_spec({**super_spec(42.0, 40.5), "drive": {"kind": "super", "params": {}}})["scalars"]
    check("18q. SUPER at registry defaults prepares > 0.8 (Bracht 2021 theory ~0.97, Piccinini 2025 ~0.95 [V]); 2 pi / 2 pi areas prepare < 5 % and the run carries the note",
          s_def["eta_prep"] > 0.8 and s_def["inversion"] > 0.8 and s_bad["eta_prep"] < 0.05
          and any("preparation probability" in n for n in s_bad["notes"]) and not any("preparation probability" in n for n in s_def["notes"]),
          f"default eta_prep {s_def['eta_prep']:.3f}, 2/2 pi eta_prep {s_bad['eta_prep']:.4f}")

    # 18r. FIX-1b incompatibility wiring: validate_spec warning, engine invalid reason, registry compat list
    from fsim_core.compose import spec as SPEC
    sp_inc = {"schema": "fsim-composed/1", "T_K": 4.2, "material": {"kind": "gan_aln", "params": {}},
              "structure": {"kind": "edge_ridge", "params": {}}, "drive": {"kind": "above_band", "params": {}},
              "setup": {"kind": "lens_fibre_det", "params": {}}, "solver": "auto"}
    vi = SPEC.validate_spec(sp_inc)
    sc_inc = engine.evaluate_spec(sp_inc)["scalars"]
    cp = R.registry_payload()["compat"]
    check("18r. edge_ridge + gan_aln: validate_spec warns, the engine run is invalid with the named reason, the registry compat list names the pair "
          "(and not inp_gainp_red)",
          any(i["level"] == "warning" and "incompatible" in i["message"] for i in vi["issues"]) and vi["ok"]
          and (not sc_inc["valid"]) and any("incompatible material / structure pair" in r for r in sc_inc["invalid_reasons"])
          and {"material": "gan_aln", "structure": "edge_ridge"}.items() <= next(c for c in cp if c["material"] == "gan_aln").items()
          and not any(c["material"] == "inp_gainp_red" for c in cp), f"{len(cp)} pairs")

    n, p = len(RESULTS), sum(RESULTS)
    print(f"{p}/{n} compose_engine checks passed")
    return 0 if p == n else 1


if __name__ == "__main__":
    sys.exit(main())
