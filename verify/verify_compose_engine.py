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

    n, p = len(RESULTS), sum(RESULTS)
    print(f"{p}/{n} compose_engine checks passed")
    return 0 if p == n else 1


if __name__ == "__main__":
    sys.exit(main())
