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
    check("7b. raw/corrected round trip", abs(OB.hom_v_corrected(OB.hom_v_raw(0.9, 0.03), 0.03) - 0.9) < 1e-12)

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

    n, p = len(RESULTS), sum(RESULTS)
    print(f"{p}/{n} compose_engine checks passed")
    return 0 if p == n else 1


if __name__ == "__main__":
    sys.exit(main())
