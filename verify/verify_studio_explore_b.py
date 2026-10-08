"""FSIM Studio explorers B (fsim_studio/api_explore_b.py): phonon spectra (qd_gf),
transport + nitride Stark, Langevin QD-laser (sde / Zhao).

Backend checks are the physics brief's section 3, 4 and 5 correctness checks. Every number an
endpoint returns must EQUAL (== on floats, no tolerance) the value the same public fsim_core call
returns here; the physics checks themselves (closed forms, limits, detailed balance, estimator
control) are evaluated on the public fsim_core API with the brief's tolerances. Each check fails
on the pre-explorer code: the endpoints do not exist there (404), and the pinned caveats, the
refusal of disable_stim / beta_sp, and the "never a compatibility verdict" rule are asserted on
the served payloads.

Run: python verify/verify_studio_explore_b.py   (exit code 0 iff all pass)
"""
import math
import os
import re
import sys
import tempfile
import time
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok), detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""), flush=True)


def rel(a, b):
    return abs(a - b) / max(abs(b), 1e-300)


def main():
    tmp = tempfile.TemporaryDirectory(prefix="fsim_studio_explore_b_")
    os.environ["FSIM_STUDIO_CACHE"] = tmp.name
    import numpy as np

    from fsim_core import loading, nitride_levels, nitride_stark, nitride_transport, qd_gf, sde, spectral, transport
    from fsim_studio import api_explore_b as eb
    from fsim_studio.server import create_app

    app = create_app()
    c = app.test_client()
    jm = app.config["JOB_MANAGER"]
    H = {"X-FSIM-Studio": "1"}

    def post(url, body):
        return c.post(url, json=body, headers=H)

    def wait(resp, timeout=240):
        """Result of a submit_call response (cached inline, or poll the job)."""
        d = resp.get_json() if hasattr(resp, "get_json") else resp
        if "result" in d:
            return d["result"]
        t0 = time.time()
        while time.time() - t0 < timeout:
            snap = c.get(f"/api/jobs/{d['job_id']}").get_json()
            if snap["state"] == "done":
                return snap["result"]
            if snap["state"] in ("error", "cancelled"):
                raise RuntimeError(f"job {d['job_id']} {snap['state']}: {snap.get('error')}")
            time.sleep(0.2)
        raise TimeoutError(d["job_id"])

    try:
        # ================================================================ phonon (brief section 3)
        meta = c.get("/api/explore/phonon/meta").get_json()
        check("phonon meta: caveat pinned verbatim + not-brightness note",
              meta["caveat"] == eb.PHONON_CAVEAT and "Independent-boson model, bulk LA deformation potential only" in meta["caveat"]
              and "Lengths are wavefunction extents, not dot sizes" in meta["caveat"]
              and "not brightness" in meta["not_brightness"])
        check("phonon meta: default alpha is the computed 0.0181 ps^2 inside 0.01-0.1; geometry tagged A, material DR",
              0.01 <= meta["alpha_default_ps2"] <= 0.1 and abs(meta["alpha_default_ps2"] - qd_gf.coupling_alpha_ps2(qd_gf.PhononParams())) == 0
              and meta["tags"]["geometry"] == "A" and meta["tags"]["material"] == "DR" and meta["tags"]["model"] == "E",
              f"alpha {meta['alpha_default_ps2']:.5f}")

        # Z(T) / S(T) curve (live)
        zt = c.get("/api/explore/phonon/zt?l_xy=4.5&l_z=1.5").get_json()
        P0 = qd_gf.PhononParams(l_xy_nm=4.5, l_z_nm=1.5)
        Ts = np.linspace(1.0, 300.0, 60)
        check("phonon zt: S_total and Z equal huang_rhys / zpl_weight exactly",
              zt["S_total"] == [qd_gf.huang_rhys(P0, float(T)) for T in Ts] and zt["Z"] == [qd_gf.zpl_weight(P0, float(T)) for T in Ts])
        check("phonon zt: Z = exp(-S_total) to 1e-12, Z in (0,1], S monotone in T",
              all(abs(z - math.exp(-s)) < 1e-12 and 0 < z <= 1 for z, s in zip(zt["Z"], zt["S_total"]))
              and all(b > a for a, b in zip(zt["S_total"], zt["S_total"][1:])))
        sph = c.get("/api/explore/phonon/zt?l_xy=3&l_z=3&alpha=0.0181").get_json()
        Ps = qd_gf.PhononParams(l_xy_nm=3.0, l_z_nm=3.0, alpha_ps2=0.0181)
        # independent closed form: S_0 = alpha_meV E_c^2 / 2 with E_c = hbar sqrt2 c_s / l (c_s in nm/ps)
        e_c = qd_gf.HBAR * math.sqrt(2.0) * (4600.0 * 1e-3) / 3.0
        s0_closed = qd_gf.coupling_alpha_meV(Ps) * e_c ** 2 / 2.0
        check("spherical limit: S_total(T=0) = alpha_meV E_c^2 / 2 (rtol 1e-4); the 1 K endpoint value equals huang_rhys exactly",
              rel(qd_gf.huang_rhys(Ps, 0.0), s0_closed) < 1e-4 and sph["S_total"][0] == qd_gf.huang_rhys(Ps, 1.0),
              f"closed {s0_closed:.8f} T=0 {qd_gf.huang_rhys(Ps, 0.0):.8f}")

        # frames job (batch of T in one worker) == direct calls
        args = {"l_xy": 4.5, "l_z": 1.5, "gamma_zpl": 0.5, "w": 2.0, "kappa": 1.0, "F_cav": 1.0, "delta": 0.0}
        t0 = time.time()
        resp = post("/api/explore/phonon/frames", args)
        out = wait(resp)
        t_frames = time.time() - t0
        check("phonon frames: one job returns all T frames", out["T"] == list(eb.PHONON_T_FRAMES) and len(out["frames"]) == len(eb.PHONON_T_FRAMES),
              f"{len(out['frames'])} frames in {t_frames:.1f} s")
        fr = {f["T"]: f for f in out["frames"]}
        f100 = fr[100.0]
        grid = np.asarray(eb.OMEGA_GRID)
        s_direct = qd_gf.ibm_spectrum(grid, P0, 100.0, 0.5)
        check("phonon frame T=100: spectrum equals ibm_spectrum exactly; S_total, Z equal huang_rhys / zpl_weight",
              f100["spectrum"]["S"] == s_direct.tolist() and f100["S_total"] == qd_gf.huang_rhys(P0, 100.0)
              and f100["Z"] == qd_gf.zpl_weight(P0, 100.0))
        area = float(np.trapezoid(np.asarray(f100["spectrum"]["S"]), grid))
        check("phonon frame: unit area on the grid to 1e-6, S >= 0 everywhere",
              abs(area - 1.0) < 1e-6 and min(f100["spectrum"]["S"]) >= 0.0, f"area {area:.9f}")
        t_ibm = qd_gf.ibm_transmission(np.asarray(f100["transmission"]["delta"]), P0, 100.0, 0.5, w_meV=2.0, kappa_meV=1.0)
        t_lor = [float(spectral.transmission(float(d), 0.5, w=2.0, kappa=1.0)) for d in f100["transmission"]["delta"]]
        check("phonon frame: ibm and Lorentzian-only transmission curves equal ibm_transmission / spectral.transmission",
              f100["transmission"]["ibm"] == np.asarray(t_ibm).tolist() and f100["transmission"]["lorentz"] == t_lor)
        pt = qd_gf.ibm_purcell_transmission(0.0, P0, 100.0, 0.5, 1.0, w_meV=2.0, kappa_meV=1.0)
        # Z_eff / rate_mult normalise by the NUMERIC sideband area, which differs from 1 - Z by ~7e-9, so 1e-7 there
        check("F_cav = 1 reproduces ibm_transmission (t to 1e-12; Z_eff == Z and rate_mult == 1 to 1e-7)",
              f100["point"]["t_purcell"] == pt[0] and abs(pt[0] - qd_gf.ibm_transmission(0.0, P0, 100.0, 0.5, w_meV=2.0, kappa_meV=1.0)) < 1e-12
              and abs(f100["point"]["Z_eff"] - f100["Z"]) < 1e-7 and abs(f100["point"]["rate_mult"] - 1.0) < 1e-7,
              f"t {pt[0]:.8f}")
        # Z_eff = Z F / (Z F + 1 - Z)
        lr = c.get("/api/explore/phonon/live?T=100&l_xy=4.5&l_z=1.5&gamma_zpl=0.5&w=2&kappa=1&F_cav=8&delta=0.3")
        live = wait(lr)["frames"][0]
        Z = qd_gf.zpl_weight(P0, 100.0)
        check("Z_eff = Z F / (Z F + 1 - Z) (1e-7: numeric sideband area) and live point equals ibm_purcell_transmission exactly",
              abs(live["point"]["Z_eff"] - Z * 8 / (Z * 8 + 1 - Z)) < 1e-7
              and (live["point"]["t_purcell"], live["point"]["Z_eff"], live["point"]["rate_mult"])
              == qd_gf.ibm_purcell_transmission(0.3, P0, 100.0, 0.5, 8.0, w_meV=2.0, kappa_meV=1.0))
        check("phonon live (one-frame pool job): spectrum / transmission equal direct calls",
              live["spectrum"]["S"] == qd_gf.ibm_spectrum(grid, P0, 100.0, 0.5).tolist()
              and live["point"]["t_ibm"] == float(qd_gf.ibm_transmission(0.3, P0, 100.0, 0.5, w_meV=2.0, kappa_meV=1.0))
              and live["gamma_fit"] == float(qd_gf.effective_gamma_zpl(
                  100.0, 0.11910667926117442, 1.0192382311813207e-16, 11.051021511907562, 18.273934741108743)))
        # narrow filter: t_IBM >= Z t_Lorentzian (w 0.5 meV, 100 K)
        nar = wait(c.get("/api/explore/phonon/live?T=100&gamma_zpl=0.5&w=0.5&kappa_on=0&F_cav=1&delta=0"))["frames"][0]
        check("narrow filter: t_IBM >= Z t_Lorentzian (w 0.5 meV, 100 K)",
              nar["point"]["t_ibm"] >= nar["Z"] * nar["point"]["t_lorentz"],
              f"{nar['point']['t_ibm']:.3f} >= {nar['Z'] * nar['point']['t_lorentz']:.3f}")
        # detailed balance as integrated weights (red/blue), not pointwise exp(-E/kT)
        g = np.linspace(-15.0, 15.0, 30001)

        def red_blue(T):
            s = qd_gf.ibm_spectrum(g, qd_gf.PhononParams(), T, 0.001)
            return np.trapezoid(s[g < -0.05], g[g < -0.05]) / np.trapezoid(s[g > 0.05], g[g > 0.05])
        r4, r300 = red_blue(4.0), red_blue(300.0)
        check("detailed balance (integrated weights): red/blue > 3 at 4 K, 1 < red/blue < 1.4 at 300 K", r4 > 3.0 and 1.0 < r300 < 1.4,
              f"{r4:.2f}, {r300:.3f}")
        # collapse to Lorentzian / wide filter
        P_tiny = qd_gf.PhononParams(alpha_ps2=qd_gf.coupling_alpha_ps2(qd_gf.PhononParams()) * 1e-6)
        gg = np.linspace(-25, 25, 50001)
        lor = spectral.lorentzian(gg, 0.0, 0.1)
        lor = lor / np.trapezoid(lor, gg)
        coll = np.max(np.abs(qd_gf.ibm_spectrum(gg, P_tiny, 77.0, 0.1) - lor)) / lor.max()
        check("alpha x 1e-6 collapses to spectral.lorentzian / tophat_transmission (1e-4)",
              coll < 1e-4 and abs(qd_gf.ibm_transmission(0.4, P_tiny, 77.0, 0.1, w_meV=1.5) - spectral.tophat_transmission(0.4, 0.1, 1.5)) < 1e-4,
              f"{coll:.2e}")
        check("w = 1e4 meV collects 1 within 1e-4 at 300 K", abs(qd_gf.ibm_transmission(0.0, P0, 300.0, 0.5, w_meV=1e4) - 1.0) < 1e-4)
        # review finding 11 (brief 3(e)): T / geometry select precomputed frames; the request thread never integrates
        from fsim_core import qd_gf as _qg
        misses0 = _qg._phi_window.cache_info().misses
        cold = c.get("/api/explore/phonon/live?T=100&l_xy=5.37&l_z=1.93&gamma_zpl=0.5&w=2&kappa=1&F_cav=1&delta=0")
        cj = cold.get_json()
        check("review 11: phonon/live at a new geometry returns a pool job (no top-level spectrum, no phi integration in the "
              "request thread: the Flask process's phi cache did not miss)",
              cold.status_code in (200, 202) and "spectrum" not in cj and ("job_id" in cj)
              and _qg._phi_window.cache_info().misses == misses0, f"{cold.status_code} misses {_qg._phi_window.cache_info().misses - misses0}")
        offgrid = [c.get(f"/api/explore/phonon/live?T={t}").status_code for t in (123, 99.5, 301)]
        check("review 11: phonon/live refuses a T that is not a precomputed frame temperature (400)", offgrid == [400] * 3, str(offgrid))
        wait(cold)
        check("phonon refuses out-of-range inputs with 400",
              c.get("/api/explore/phonon/live?T=400").status_code == 400
              and c.get("/api/explore/phonon/live?l_xy=20").status_code == 400
              and c.get("/api/explore/phonon/zt?alpha=1").status_code == 400
              and post("/api/explore/phonon/frames", {"w": 99}).status_code == 400)

        # ================================================================ transport (brief section 4)
        tm = c.get("/api/explore/transport/meta").get_json()
        check("transport meta: InP and nitride caveats pinned verbatim; Zhang is a labelled guide, non-gating",
              tm["caveat"] == eb.TRANSPORT_CAVEAT and "leak_valleys='G' is a LOWER bound on leakage" in tm["caveat"]
              and tm["nitride_caveat"] == eb.NITRIDE_CAVEAT and "Zhang is 10 K, x=0.15, another device" in tm["nitride_caveat"]
              and tm["zhang"]["slope_meV_per_V"] == -10.0 and "guide only" in tm["zhang"]["label"])
        q = "I_uA=1&T=230&n_dot=2e9&aperture=1&tau_pulse=0.1&w=2&dE_WL=100&tau_rad=1"
        pr = c.get(f"/api/explore/transport/point?preset=red&{q}").get_json()
        red = transport.red_diode_preset()
        d_ref = transport.evaluate_injection(red, 1.0, 230.0, 2e9, 1.0, 0.1, 2.0, 100.0, tau_rad_ns=1.0)
        check("transport point equals evaluate_injection exactly (V_j, V_applied, eta_inj, mu, b_e, P_junction, F, C_dep)",
              (pr["V_j"], pr["V_applied"], pr["eta_inj"], pr["mu"], pr["b_e"], pr["P_junction_W"], pr["F_kVcm"], pr["C_dep_pF"])
              == (d_ref.V_j, d_ref.V_applied, d_ref.leakage.eta_inj, d_ref.mu, d_ref.b_e, d_ref.P_junction_W,
                  d_ref.depletion.F_kVcm, d_ref.depletion.C_dep_pF))
        flux = pr["eta_inj"] * 1e-6 / 1.602176634e-19
        check("r_captured + r_matrix == eta_inj I / q (rtol 1e-12), with and without E_X",
              rel(pr["r_captured"] + pr["r_matrix"], flux) < 1e-12
              and (lambda p2: rel(p2["r_captured"] + p2["r_matrix"], p2["eta_inj"] * 1e-6 / 1.602176634e-19) < 1e-12 and p2["f_qfl"] <= 1.0)(
                  c.get(f"/api/explore/transport/point?preset=red&{q}&E_X_eV=2.2").get_json()),
              f"{pr['r_captured'] + pr['r_matrix']:.15e}")
        U = pr["E_U_meV"]
        assert U == d_ref.background.E_U_meV
        check("xi_window == exp(-dE/E_U) 2 sinh(w/2E_U) for dE >= w/2 (1e-9), xi <= 1",
              rel(pr["xi"], math.exp(-100.0 / U) * 2.0 * math.sinh(2.0 / (2.0 * U))) < 1e-9 and 0 < pr["xi"] <= 1.0, f"xi {pr['xi']:.7e}")
        check("eta_inj in (0, 1] for all three presets",
              all(0 < c.get(f"/api/explore/transport/point?preset={p}&{q}").get_json()["eta_inj"] <= 1.0 for p in eb.PRESETS))
        dep = c.get("/api/explore/transport/depletion?preset=gaas&T=300").get_json()
        g = transport.gaas_homojunction()
        vbi_ref = transport.homojunction_vbi(g.active, 1e17, 1e17, 300.0)
        check("GaAs homojunction V_bi == (kT/q) ln(N_A N_D / n_i^2) (1e-12)", rel(dep["V_bi"], vbi_ref) < 1e-12, f"{dep['V_bi']:.8f}")
        check("depletion: flat band (F = 0) for V_j >= V_bi, W and C_dep continuous at V_bi, curve equals Diode.depletion",
              dep["flat_band"][-1] is True and dep["F_kVcm"][-1] == 0.0 and dep["F_kVcm"][-2] < 1e-6 + 1e-3
              and abs(dep["W_nm"][-1] - dep["W_nm"][-2]) / dep["W_nm"][-1] < 1e-6 and abs(dep["C_dep_pF"][-1] - dep["C_dep_pF"][-2]) / dep["C_dep_pF"][-1] < 1e-6
              and dep["F_kVcm"][10] == g.depletion(dep["V_j"][10], 300.0).F_kVcm,
              f"F(V_bi-) {dep['F_kVcm'][-2]:.3e}")
        sw = wait(post("/api/explore/transport/sweep", {"preset": "red", "T": 230, "n_dot": 2e9, "aperture": 1, "tau_pulse": 0.1,
                                                         "w": 2, "dE_WL": 100, "tau_rad": 1}))
        cols = sw["columns"]
        check("I sweep: 50 points equal evaluate_injection (V_j, mu, b_e exact); V_j monotone in I",
              len(cols["I_uA"]) == 50 and all(b > a for a, b in zip(cols["V_j"], cols["V_j"][1:]))
              and cols["V_j"][20] == transport.evaluate_injection(red, cols["I_uA"][20], 230.0, 2e9, 1.0, 0.1, 2.0, 100.0, tau_rad_ns=1.0).V_j
              and cols["b_e"][37] == transport.evaluate_injection(red, cols["I_uA"][37], 230.0, 2e9, 1.0, 0.1, 2.0, 100.0, tau_rad_ns=1.0).b_e)
        check("legacy b_e fit is shown with its max_rel_err (to_background_channel)", sw["bg_fit"]["max_rel_err"] >= 0.0 and "b_e_Eact" in sw["bg_fit"],
              f"max_rel_err {sw['bg_fit']['max_rel_err']:.3g}")
        unsat = [b for b, s in zip(cols["b_e"], cols["saturated"]) if not s]
        sat = [(i, b) for i, b, s in zip(cols["I_uA"], cols["b_e"], cols["saturated"]) if s]
        flat_ok = len(unsat) < 2 or (max(unsat) / min(unsat) - 1) < 1e-9
        lin_ok = len(sat) < 2 or all(abs((b / i) / (sat[0][1] / sat[0][0]) - 1) < 1e-6 for i, b in sat)
        check("b_e independent of I below saturation (rtol 1e-9) and linear above (docstring item 7)", flat_ok and lin_ok,
              f"{len(unsat)} unsaturated, {len(sat)} saturated points")

        # nitride
        npnt = c.get("/api/explore/transport/nitride_point?I_uA=0.02&T=300&n_dot=2e9&aperture=1&tau_pulse=0.1&w=2&dE_WL=100&tau_rad=1").get_json()
        nd = nitride_transport.NitrideDiode()
        b_ref = nitride_stark.resolve_bias(nd, T_j_K=300.0, current_uA=0.02)
        nr = nitride_transport.evaluate_injection(nd, 0.02, 300.0, 2e9, 1.0, 0.1, 2.0, 100.0, tau_rad_ns=1.0)
        check("nitride point equals nitride_transport.evaluate_injection and resolve_bias exactly",
              (npnt["V_j"], npnt["mu"], npnt["b_e"]) == (nr.V_j, nr.mu, nr.b_e)
              and npnt["bias"]["V_j"] == b_ref["V_j"] and npnt["bias"]["applied_field_kVcm"] == b_ref["applied_field_kVcm"])
        rt = nitride_stark.resolve_bias(nd, T_j_K=300.0, junction_voltage_V=nitride_stark.resolve_bias(nd, T_j_K=300.0, current_uA=1.0)["V_j"])
        check("resolve_bias round trip I -> V_j -> I (1e-12 relative)", abs(rt["current_uA"] / 1.0 - 1.0) < 1e-12, f"{rt['current_uA']!r}")
        # Stark traces (jobs)
        st = post("/api/explore/transport/stark", {"orientation": "c_plane", "polarity": 1, "T": 300}).get_json()
        traces = {t["screening_fraction"]: wait(t) for t in st["traces"]}
        tr0 = traces[0.0]
        rows_ref = []
        for i, V in enumerate(eb.STARK_V):
            b = nitride_stark.resolve_bias(nd, T_j_K=300.0, junction_voltage_V=V, field_polarity=1, external_field_kVcm=0.0)
            lv = nitride_levels.levels(nitride_levels.NitrideDotSystem(orientation="c_plane", screening_fraction=0.0,
                                                                       external_field_kVcm=b["applied_field_kVcm"]), 300.0)
            rows_ref.append({"row_id": i, "V_j": V, "E_X_eV": lv.E_X_eV, "spectroscopy_valid": bool(lv.valid and b["bias_valid"]),
                             "depletion_regime": b["depletion_regime"], "flat_band": b["flat_band"]})
        der_ref = nitride_stark.stark_derivatives(rows_ref)
        check("Stark trace (c-plane, screening 0): E_X(V_j) and dE_X/dV equal levels() + stark_derivatives exactly",
              [r["E_X_eV"] for r in tr0["rows"]] == [r["E_X_eV"] for r in rows_ref]
              and [r["dE_X_dV_meV_per_V"] for r in tr0["rows"]] == [r["dE_X_dV_meV_per_V"] for r in der_ref]
              and len(tr0["rows"]) == 17)
        sl = [r["dE_X_dV_meV_per_V"] for r in tr0["rows"] if r["derivative_valid"]]
        check("c-plane slope is negative at polarity +1 over 0..3.2 V (brief: -11.1 to -43.0 meV/V)", len(sl) >= 10 and max(sl) < 0,
              f"{min(sl):.1f} to {max(sl):.1f} meV/V")
        check("three screening series (0, 0.5, 1) computed; screening changes E_X",
              sorted(traces) == [0.0, 0.5, 1.0] and traces[1.0]["rows"][4]["E_X_eV"] != tr0["rows"][4]["E_X_eV"])
        a_tr = wait(post("/api/explore/transport/stark", {"orientation": "a_plane", "polarity": 1, "T": 300, "screenings": [0.0]}).get_json()["traces"][0])
        check("a-plane level field_kVcm == 0 exactly with no external field", all(r["field_kVcm"] == 0.0 for r in [
            {"field_kVcm": nitride_levels.levels(nitride_levels.NitrideDotSystem(orientation="a_plane", external_field_kVcm=0.0), 300.0).field_kVcm}]),
              "the Studio trace applies the diode field, so its level field is the applied field (checked below)")
        check("a-plane trace: level field equals the applied diode field (no polarization term)",
              all(abs(r["field_kVcm"] - r["applied_field_kVcm"]) < 1e-9 for r in a_tr["rows"] if r["spectroscopy_valid"]))
        flat = eb.stark_trace("c_plane", 0.0, 1, 300.0, 0.0, [3.0, 3.2, 3.39, 3.5])
        fr_rows = {r["V_j"]: r for r in flat["rows"]}
        check("flat_band at V_bi + 0.1 V; derivatives never cross the regime kink",
              fr_rows[3.39]["flat_band"] is True and fr_rows[3.2]["flat_band"] is False and fr_rows[3.2]["derivative_valid"] is False
              and fr_rows[3.39]["derivative_valid"] is False, f"V_bi {nd.vbi(300.0):.4f}")
        blob = repr(tr0) + repr(npnt) + repr(tm)
        src = (ROOT / "fsim_studio" / "web" / "js" / "workspaces" / "explorers" / "transport.js")
        js = src.read_text(encoding="utf-8") if src.exists() else ""
        check("never a 'compatible' verdict on a slope (screening_compatibility stays out of the explorer)",
              "compatib" not in blob.lower() and "compatib" not in js.lower() and "screening_compatibility" not in Path(eb.__file__).read_text(encoding="utf-8").split('"""')[2])

        # ================================================================ sde (brief section 5)
        sm = c.get("/api/explore/sde/meta").get_json()
        check("sde meta: Zhao caveat pinned verbatim (normal 1.7 sigma, quiet does NOT: 5.8 sigma, FAIL both rows)",
              sm["caveat"] == eb.SDE_CAVEAT and "quiet pump does NOT: 1.017 vs 0.982, 5.8 sigma" in sm["caveat"]
              and "not a single-photon source; g2 ~ 1 is correct physics" in sm["caveat"] and "I/I_th, not amperes" in sm["axis_note"])
        card = c.get("/api/explore/sde/card").get_json()
        pts = {p["pump"]: p for p in card["card_points"]}
        check("sde card: 1.0224+-0.003 normal and 0.9823+-0.006 quiet at I/I_th = 4 [V]; comparison verdict FAIL both rows verbatim",
              (pts["normal"]["g2"], pts["normal"]["err"], pts["quiet"]["g2"], pts["quiet"]["err"]) == (1.0224, 0.003, 0.9823, 0.006)
              and all(p["I_over_Ith"] == 4.0 for p in pts.values()) and card["card_tag"] == "V"
              and [r["verdict"] for r in card["comparison"]] == ["FAIL", "FAIL"] and len(card["curve"]) == 3)
        li = wait(post("/api/explore/sde/li", {}))
        P = sde.QDLaserParams()
        ssd = sde.steady_state(P, 4.0)
        k4 = li["I"].index(4.0)
        check("L-I job: steady_state at I = 4 equals the direct call exactly; residual < 1e-9; rho in [0,1]; S monotone",
              li["S"][k4] == ssd["S"] and max(li["residual"]) < 1e-9 and all(0 <= r <= 1 for r in li["rho_ES"] + li["rho_GS"])
              and all(b > a for a, b in zip(li["S"], li["S"][1:])), f"residual {max(li['residual']):.2e}")
        ith = sde.find_threshold(P)
        check("find_threshold: equals the direct call, reproducible to 1e-9, inside (0.1, 5)",
              li["I_th"] == ith and abs(ith - sde.find_threshold(P)) < 1e-9 and 0.1 < ith < 5.0, f"I_th {ith:.4f}")
        fr1 = wait(post("/api/explore/sde/frames", {"I_grid": [4.0], "F_pumps": [1.0], "n_runs": 8, "t_end": 1e-9, "seed": 1}).get_json()["frames"][0])
        d1 = sde.g2_vs_pump(P, 4.0, n_runs=8, t_end=1e-9, dt=1e-13, seed=1, F_pump=1.0)
        check("sde frame equals g2_vs_pump exactly (seeded): g2_mean, g2_se, mean_S; SE finite and positive",
              fr1["g2_mean"] == d1["g2_mean"] and fr1["g2_se"] == d1["g2_se"] and fr1["mean_S"] == d1["mean_S"] and fr1["g2_se"] > 0,
              f"g2 {fr1['g2_mean']:.5f} +- {fr1['g2_se']:.5f}")
        tr = wait(post("/api/explore/sde/trace", {"I": 4.0, "F_pump": 1.0, "t_end": 1e-9, "seed": 1}))
        sim = sde.simulate(P, 4.0, 1e-9, 1e-13, seed=1, F_pump=1.0)
        check("trace job: g2_0 / g2_0_ext equal simulate() exactly; burn-in fraction reported (0.1)",
              tr["g2_0"] == sim["g2_0"] and tr["g2_0_ext"] == sim["g2_0_ext"] and tr["burnin_frac"] == sde.BURNIN_FRAC == 0.1
              and tr["n_burn"] == sim["n_burn"])
        check("thin_stats leaves g2 invariant (1e-9) while the Fano factor follows loading.f8b_thin_fano",
              abs(tr["g2_0"] - tr["g2_0_ext"]) < 1e-9 and rel(tr["fano_ext"], loading.f8b_thin_fano(tr["eta_ext"], tr["fano"])) < 1e-9,
              f"g2 {tr['g2_0']:.5f} = {tr['g2_0_ext']:.5f}")
        # estimator control (verify recipe ONLY: beta_sp 0.35, 2 I_th, disable_stim, n_runs 24, t_end 1e-9)
        Pt = sde.QDLaserParams(beta_sp=0.35)
        ithp = sde.find_threshold(Pt, I_hi=8.0)
        rctl = sde.g2_vs_pump(Pt, 2.0 * ithp, n_runs=24, t_end=1e-9, dt=1e-13, seed=1, disable_stim=True, F_pump=1.0)
        check("estimator control (verify recipe): disable_stim, beta_sp 0.35, 2 I_th -> g2 = 1 within max(4 SE, 0.01)",
              abs(rctl["g2_mean"] - 1.0) < max(4.0 * rctl["g2_se"], 0.01), f"g2 {rctl['g2_mean']:.4f} +- {rctl['g2_se']:.4f}")
        bad = [post("/api/explore/sde/frames", {"disable_stim": True}), post("/api/explore/sde/frames", {"beta_sp": 0.35}),
               post("/api/explore/sde/trace", {"N_dots": 100}), post("/api/explore/sde/li", {"disable_stim": 1})]
        check("disable_stim, beta_sp and the tuned Tier-2 parameters are refused with 400 (never exposed)",
              all(b.status_code == 400 and "not exposed" in b.get_json()["error"] for b in bad))
        # dt halving
        a = sde.g2_vs_pump(P, 4.0, n_runs=8, t_end=1e-9, dt=1e-13, seed=3, F_pump=1.0)
        b = sde.g2_vs_pump(P, 4.0, n_runs=8, t_end=1e-9, dt=5e-14, seed=3, F_pump=1.0)
        se = max(a["g2_se"], b["g2_se"])
        check("dt halving moves g2 by less than 4 SE + 5e-3", abs(a["g2_mean"] - b["g2_mean"]) < 4 * se + 5e-3,
              f"{abs(a['g2_mean'] - b['g2_mean']):.4f} < {4 * se + 5e-3:.4f}")
        wins = 0
        for s in range(6):
            n = eb.sde_frame(4.0, 1.0, 20, 1e-9, 1e-13, 100 + s)
            qd = eb.sde_frame(4.0, 0.08, 20, 1e-9, 1e-13, 100 + s)
            wins += qd["g2_mean"] < n["g2_mean"]
        check("quiet < normal in 6 of 6 seed batches (a sign test only, not a magnitude z-test)", wins == 6, f"{wins}/6")
        check("sde refuses out-of-range I, n_runs, t_end",
              post("/api/explore/sde/frames", {"I_grid": [9.0]}).status_code == 400
              and post("/api/explore/sde/frames", {"n_runs": 500}).status_code == 400
              and post("/api/explore/sde/trace", {"t_end": 1.0}).status_code == 400)

        # ================================================================ review finding 2: transport tags follow the widest rule
        tmeta = c.get("/api/explore/transport/meta").get_json()["tags"]
        # I, T, n_dot, aperture, tau_pulse, tau_rad are [A] controls (transport.js); w, dE_WL [E]. Every readout depends on
        # at least one [A] control, so the widest-tag rule gives A for all of them (the old table said DR / E).
        outs = ("V_j", "V_bi", "F", "C_dep", "eta_inj", "mu", "b_e", "P_junction", "r_dot", "xi")
        check("review 2: transport output tags are the widest of the model line and the controls each depends on (all [A]); "
              "none stays DR/E although I, T, n_dot, aperture are [A]",
              all(tmeta.get(k) == "A" for k in outs), str({k: tmeta.get(k) for k in outs}))
        pmeta = c.get("/api/explore/phonon/meta").get_json()["tags"]
        check("review 2: phonon result tag is the widest of geometry/alpha/material/model/filter/cavity tags",
              pmeta["result"] == "A" and set(pmeta.values()) - {"A"} <= {"DR", "E"}, str(pmeta))
        import re as _re
        jsd = ROOT / "fsim_studio/web/js/workspaces/explorers"
        tsrc = (jsd / "transport.js").read_text(encoding="utf-8")
        psrc = (jsd / "phonon.js").read_text(encoding="utf-8")
        check("review 2: transport.js and phonon.js hover templates / annotations carry no hard-coded tag text",
              not any(_re.search(r"\[(?:DR|E|V|A)\]", seg) for src_ in (tsrc.replace("μ = 0.5 [A]", ""), psrc)
                      for seg in _re.findall(r"hovertemplate: [^<]*<extra>|hovertemplate: \"[^\"]*\"|text: `[^`]*`", src_)), "")

        # ================================================================ review finding 8: malformed bodies are 400, never 500
        bad = [post("/api/explore/sde/frames", {"I_grid": 5}), post("/api/explore/sde/frames", {"F_pumps": 1}),
               post("/api/explore/sde/frames", {"I_grid": ["x"]}), post("/api/explore/sde/frames", {"I_grid": [[1]]}),
               post("/api/explore/sde/li", {"I_grid": ["x"]}), post("/api/explore/sde/li", {"I_grid": 3}),
               post("/api/explore/sde/li", {"I_grid": [None]}),
               post("/api/explore/transport/stark", {"screenings": 1}), post("/api/explore/transport/stark", {"screenings": ["a"]}),
               post("/api/explore/transport/stark", {"screenings": {"a": 1}}), post("/api/explore/transport/stark", {"screenings": [[0]]})]
        codes = [r.status_code for r in bad]
        check("review 8: sde/frames I_grid 5 / F_pumps 1 / ['x'] / [[1]], sde/li ['x'] / 3 / [None], transport/stark screenings 1 / ['a'] / {} / [[0]] "
              "all answer 400 with an error message (no 500)",
              codes == [400] * len(bad) and all("error" in (r.get_json() or {}) for r in bad), str(codes))
    finally:
        jm.shutdown()

    n_pass = sum(1 for _, ok, _ in CHECKS if ok)
    print(f"\n{n_pass}/{len(CHECKS)} studio explorer B checks passed")
    try:
        tmp.cleanup()
    except OSError:
        pass
    return 0 if n_pass == len(CHECKS) else 1


if __name__ == "__main__":
    sys.exit(main())
