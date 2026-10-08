"""FSIM Studio model explorers A (fsim_studio/api_explore_a.py): cw_g2, pulse_counting, lindblad.

Backend checks = the correctness checks of .workers/studio/phase2/physics-brief.md sections 1, 2 and 6.
Every endpoint number must EQUAL (== on floats, no tolerance) the value the same fsim_core call returns
here; the physics anchors (closed forms, published values) are computed in this file independently of
fsim_core. A check that goes through the endpoint fails on a tree without api_explore_a.py (404).
"[direct]" marks a brief anchor that exercises fsim_core only (kept so the brief's list is complete).

Run: python verify/verify_studio_explore_a.py   (exit code 0 iff all pass)
"""
import math
import os
import sys
import tempfile
import warnings
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
warnings.filterwarnings("ignore")

CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok), detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""), flush=True)


def main():
    tmp = tempfile.TemporaryDirectory(prefix="fsim_studio_explore_a_")
    os.environ["FSIM_STUDIO_CACHE"] = tmp.name
    import json

    import numpy as np

    from fsim_core import cw_g2 as cw
    from fsim_core import lindblad as lb
    from fsim_core import pulse_counting as pc
    from fsim_core.loading import f1b_g2
    from fsim_studio import api_explore_a as ea
    from fsim_studio import campaigns
    from fsim_studio.server import create_app

    app = create_app()
    c = app.test_client()

    def get(url):
        r = c.get(url)
        return r.status_code, (r.get_json() if r.is_json else None)

    def ok200(url):
        code, d = get(url)
        assert code == 200, (url, code, d)
        return d

    fit = json.loads((ROOT / "out/phase0/fit_params.json").read_text(encoding="utf-8"))["params"]
    try:
        # ================================================================ controls / tags
        code, ctl = get("/api/explore/a/controls")
        check("controls: endpoint answers", code == 200 and "controls" in ctl)
        cw_c = {s["id"]: s for s in ctl["controls"]["cw_g2"]}
        check("controls cw_g2: defaults and tags are the brief's (r .5 A, gamma_X 1 E, eps .1 A, rho .9 A, IRF 500 V)",
              (cw_c["r"]["default"], cw_c["r"]["tag"]) == (0.5, "A") and (cw_c["gamma_X"]["default"], cw_c["gamma_X"]["tag"]) == (1.0, "E")
              and (cw_c["eps"]["default"], cw_c["eps"]["tag"]) == (0.1, "A") and (cw_c["rho"]["default"], cw_c["rho"]["tag"]) == (0.9, "A")
              and (cw_c["irf_fwhm_ps"]["default"], cw_c["irf_fwhm_ps"]["tag"]) == (500.0, "V"))
        pc_c = {s["id"]: s for s in ctl["controls"]["pulse_counting"]}
        check("controls pulse_counting: tau_on 0.1 [A], period 12.5 [A], range 0.01-2 ns",
              (pc_c["tau_on"]["default"], pc_c["tau_on"]["tag"], pc_c["tau_on"]["min"], pc_c["tau_on"]["max"]) == (0.1, "A", 0.01, 2.0)
              and (pc_c["period"]["default"], pc_c["period"]["tag"]) == (12.5, "A"))
        check("controls: every control carries a tag",
              all(s["tag"] in ("V", "DR", "E", "A") for grp in (ctl["controls"]["cw_g2"], ctl["controls"]["pulse_counting"],
                                                                *ctl["controls"]["lindblad"].values()) for s in grp))

        # ================================================================ 1. cw_g2
        d = ok200("/api/explore/cw_g2")
        rep = cw.cw_report(0.5, 1.0, 2.0, 0.0, 0.0, 1.0, 0.1, 0.9, 500.0, 10.0, 1.0, "exponential")
        q = d["readouts"]
        check("cw_g2 default: readouts equal direct cw_report (== no tolerance)",
              all(q[k] == rep[k] for k in ("g2_dot0", "g2_meas0", "g2_raw0", "tau_dip", "A_dip", "g2_intrinsic_from_raw")))
        check("cw_g2 default: g2_dot0 = 0.1473922902 (brief anchor, 1e-9)", abs(q["g2_dot0"] - 0.1473922902) < 1e-9, f"{q['g2_dot0']:.10f}")
        check("cw_g2 default: curves equal cw_report curves exactly (full grid, not thinned)",
              d["n_tau"] == len(rep["curves"]["tau"]) and not d["decimated"] and d["curves"]["g2_dot"] == rep["curves"]["g2_dot"].tolist()
              and d["curves"]["g2_raw"] == rep["curves"]["g2_raw"].tolist() and d["curves"]["tau"] == rep["curves"]["tau"].tolist())
        check("cw_g2: g2_cw_zero == g2_cw(0) to 1e-10 and both consistency rows ok",
              abs(q["g2_dot0"] - float(cw.g2_cw(np.array([0.0]), 0.5, 1.0, 2.0, 1.0, 0.1)[0])) < 1e-10 and all(x["ok"] for x in d["checks"]))
        check("cw_g2: Plot B = g2_cw_zero on the r grid, limits = g2_cw_zero_low_pump and gamma_X/(eps gamma_XX)",
              d["vs_r"]["g2_dot0"] == cw.g2_cw_zero(np.array(d["vs_r"]["r"]), 1.0, 2.0, 0.1, 0.0, 0.0, 1.0).tolist()
              and d["vs_r"]["low_limit"] == cw.g2_cw_zero_low_pump(1.0, 2.0, 0.1, 0.0, 0.0, 1.0)
              and abs(d["vs_r"]["high_limit"] - 5.0) < 1e-12)
        # two-level limit: pump_ratio 0
        d0 = ok200("/api/explore/cw_g2?pump_ratio=0&r=0.5&gamma_X=1&irf_fwhm_ps=0")
        tau = np.array(d0["curves"]["tau"])
        ref = 1.0 - np.exp(-(0.5 + 1.0) * np.abs(tau))
        check("cw_g2 two-level limit (pump_ratio 0): g2_dot = 1 - exp(-(r+gamma_X+k_X)|tau|) to 1e-6",
              np.max(np.abs(np.array(d0["curves"]["g2_dot"]) - ref)) < 1e-6)
        # g2 -> 1 at tau = 50/min(rate)
        dl = ok200("/api/explore/cw_g2?r=5&gamma_X=5&tau_max=20&irf_fwhm_ps=0")
        rates = [5.0, 5.0, 5.0, 10.0]
        check("cw_g2: g2_dot -> 1 at tau = 50/min(rate) to 1e-6", 50.0 / min(rates) <= 20.0 and abs(dl["curves"]["g2_dot"][-1] - 1.0) < 1e-6)
        # eps = 0
        de = ok200("/api/explore/cw_g2?eps=0")
        mid = len(de["curves"]["tau"]) // 2
        check("cw_g2: eps = 0 gives g2_dot(0) = 0 exactly (closed form and curve centre), no high-pump limit offered",
              de["readouts"]["g2_dot0"] == 0.0 and de["curves"]["tau"][mid] == 0.0 and de["curves"]["g2_dot"][mid] == 0.0
              and de["vs_r"]["high_limit"] is None and all(v == 0.0 for v in de["vs_r"]["g2_dot0"]))
        # low pump limit within 2%
        dlo = ok200("/api/explore/cw_g2?r=0.001")
        check("cw_g2 low-pump limit: r = 1e-3 gamma_X is within 2% of eps p S_XX/S_X",
              abs(dlo["readouts"]["g2_dot0"] / 0.1 - 1.0) < 0.02, f"{dlo['readouts']['g2_dot0']:.5f}")
        # bunching
        db = ok200("/api/explore/cw_g2?r=20&eps=0.1")
        gr = db["vs_r"]["g2_dot0"]
        check("cw_g2 bunching: r=20, eps=0.1 gives 2.456 < 5.0 = gamma_X/(eps gamma_XX), flagged as cascade bunching",
              abs(db["readouts"]["g2_dot0"] - 2.456) < 5e-4 and db["readouts"]["g2_dot0"] < db["vs_r"]["high_limit"] == 5.0 and db["bunching"] is True
              and ea.CAVEAT_CW_BUNCH in db["caveat_bunching"], f"{db['readouts']['g2_dot0']:.4f}")
        check("cw_g2: Plot B is monotone toward the high-pump limit (eps = 0.1)",
              all(b >= a for a, b in zip(gr, gr[1:])) and gr[-1] < 5.0)
        # detailed balance (independent of the module's own steady_state)
        M = cw.generator(0.5, 1.0, 2.0, 0.0, 0.0, 1.0)
        P = cw.steady_state(0.5, 1.0, 2.0, 0.0, 0.0, 1.0)
        check("cw_g2 [direct]: generator @ steady_state == 0 to 1e-12 and P_X/P_0 = r/(gamma_X + k_X)",
              np.max(np.abs(M @ P)) < 1e-12 and abs(P[1] / P[0] - 0.5 / 1.0) < 1e-12)
        # background law
        gm = np.array(d["curves"]["g2_meas"])
        gd = np.array(d["curves"]["g2_dot"])
        check("cw_g2: g2_meas = 1 - rho^2 (1 - g2_dot) pointwise to 1e-12", np.max(np.abs(gm - (1 - 0.9 ** 2 * (1 - gd)))) < 1e-12)
        # evenness
        check("cw_g2: curves even in tau to 1e-12",
              np.max(np.abs(gd - gd[::-1])) < 1e-12 and np.max(np.abs(np.array(d["curves"]["g2_raw"]) - np.array(d["curves"]["g2_raw"])[::-1])) < 1e-12)
        # IRF kernel / deconvolution / Reischle anchor (direct)
        ok_k = all(abs(cw.irf_kernel(0.01, 500.0, sh).sum() - 1.0) < 1e-12 for sh in ("gaussian", "exponential"))
        check("cw_g2 [direct]: IRF kernel sums to 1 (1e-12) for both shapes", ok_k)
        worst = 0.0
        for sh in ("gaussian", "exponential"):
            A = 0.85
            raw = 1.0 - A * cw._dip_factor(1.13, 500.0, sh)
            worst = max(worst, abs(cw.deconvolve_dip(raw, 1.13, 500.0, sh) - A))
        check("cw_g2 [direct]: deconvolve_dip(dip_convolved(A)) = A to 1e-9, both shapes", worst < 1e-9, f"{worst:.1e}")
        reis = cw.dip_convolved_exp(0.85, 1.13, 0.5)
        td = 0.5 * (1 - 0.41) / (0.85 - 0.59)
        check("cw_g2 [direct]: Reischle anchor tau_d = 1.13 ns and dip_convolved_exp(0.85, 1.13, 0.5) = 0.410",
              abs(td - 1.13) < 5e-3 and abs(reis - 0.410) < 5e-3, f"tau_d {td:.3f}, {reis:.4f}")
        raws = [ok200(f"/api/explore/cw_g2?irf_fwhm_ps={w}")["readouts"]["g2_raw0"] for w in (0, 100, 300, 500, 1000)]
        check("cw_g2: g2_raw(0) monotone in IRF FWHM (through the endpoint)", all(b > a for a, b in zip(raws, raws[1:])), " ".join(f"{x:.3f}" for x in raws))
        # escape (T) path
        dT = ok200("/api/explore/cw_g2?T=230")
        kX, kXX = cw.escape_rates_from_retention(1.0, fit["a_esc"], fit["E_a"], fit["b_p"], fit["E_b"], 230.0)
        check("cw_g2: T = 230 K gives k_X = 9.16/ns equal to escape_rates_from_retention with the V-a fit, tagged [DR]-or-wider",
              dT["k_X"] == kX and dT["k_XX"] == kXX and abs(kX - 9.16) < 5e-3 and dT["stiff"] is False and dT["tags"]["k_X"] in ("DR", "E", "A"), f"{kX:.4f}")
        rT = cw.cw_report(0.5, 1.0, 2.0, kX, kXX, 1.0, 0.1, 0.9, 500.0, 10.0, 1.0, "exponential")
        check("cw_g2: T = 230 K readouts equal direct cw_report with those rates", dT["readouts"]["g2_raw0"] == rT["g2_raw0"] and dT["readouts"]["g2_dot0"] == rT["g2_dot0"])
        dS = ok200("/api/explore/cw_g2?T=350&r=100")
        check("cw_g2: stiff escape (T = 350 K) answers, thins the wire arrays and flags stiff",
              dS["stiff"] is True and len(dS["curves"]["tau"]) <= ea.WIRE_MAX and dS["decimated"] is True and 0.0 in dS["curves"]["tau"])
        # tags, caveat, gates, validation
        check("cw_g2: pinned caveat text is exact", d["caveat"] == ea.CAVEAT_CW and "not valid for resonant (Rabi) drive" in d["caveat"]
              and "g2_raw is what a 500 ps HBT measures" in d["caveat"])
        check("cw_g2: every readout has a tag (V/DR/E/A)", all(d["tags"][k] in ("V", "DR", "E", "A") for k in ("g2_dot0", "g2_meas0", "g2_raw0", "tau_dip", "A_dip", "g2_intrinsic_from_raw")))
        check("cw_g2: the 0.5 line is the /api/gates value, not typed", d["gate_g2"]["value"] == campaigns.gates()["g2_ceiling"]["value"])
        codes = [get(u)[0] for u in ("/api/explore/cw_g2?r=1000", "/api/explore/cw_g2?irf_shape=box", "/api/explore/cw_g2?eps=-0.1", "/api/explore/cw_g2?T=1", "/api/explore/cw_g2?r=abc")]
        check("cw_g2: out-of-range or malformed inputs answer 400", codes == [400] * 5, str(codes))

        # ================================================================ 2. pulse_counting
        p = ok200("/api/explore/pulse_counting")
        pt = pc.pulse_g2(10.0, 1.0, 2.0, 0.0, 0.0, 1.0, 0.1, 0.1, 12.4, 1.0, split=True, gate_ns=None, adjacent=True)
        dd = pc.deterministic_cycle_g2(1.0, 2.0, 0.0, 0.0, 1.0, 0.1, 12.5, eta_load=1.0, gate_ns=None, adjacent=True)
        P0 = p["point"]
        check("pulse_counting default: point readouts equal direct pulse_g2 (== no tolerance)",
              all(P0[k] == pt[k] for k in ("g2", "g2_adj", "mean_counts", "mean_counts_x", "mean_counts_xx", "adjacent_peak_factor", "converged")))
        check("pulse_counting default: deterministic readouts equal direct deterministic_cycle_g2",
              all(p["deterministic"][k] == dd[k] for k in ("g2", "g2_adj", "mean_counts", "blocked_load_probability", "one_pair_valid", "converged")))
        rr = np.array(p["vs_r"]["r"])
        direct_g2 = [pc.pulse_g2(float(x), 1.0, 2.0, 0.0, 0.0, 1.0, 0.1, 0.1, 12.4, 1.0, gate_ns=None, adjacent=True)["g2"] for x in rr[::6]]
        check("pulse_counting: Plot A g2 sweep and f1b overlay equal direct calls",
              p["vs_r"]["g2"][::6] == direct_g2 and p["vs_r"]["f1b_g2"] == np.asarray(f1b_g2(rr * 0.1, 0.1), float).tolist())
        per = np.array(p["vs_period"]["period"])
        dg = [pc.deterministic_cycle_g2(1.0, 2.0, 0.0, 0.0, 1.0, 0.1, float(x), eta_load=1.0, gate_ns=None)["g2"] for x in per[::7]]
        check("pulse_counting: Plot C deterministic sweep equals direct calls", p["vs_period"]["g2"][::7] == dg)
        # instantaneous limit (direct; tau_on = 1e-6 is below the slider range)
        inst = pc.pulse_g2(1e6, 1.0, 2.0, 0.0, 0.0, 1.0, 0.1, 1e-6, 12.4, 1.0)["g2"]
        refi = float(f1b_g2(1.0, 0.1))
        check("pulse_counting [direct]: tau_on = 1e-6, r tau_on = 1 reproduces f1b_g2(1, 0.1) to rtol 1e-4", abs(inst / refi - 1) < 1e-4, f"{inst:.7f} vs {refi:.7f}")
        pi = ok200("/api/explore/pulse_counting?r=100&tau_on=0.01&eps=0.1")
        check("pulse_counting: mu = r tau_on = 1 at tau_on = 0.01 sits within 5% of the instantaneous f1b value",
              abs(pi["point"]["g2"] / pi["point"]["f1b_g2_inst"] - 1) < 0.05, f"{pi['point']['g2']:.5f} vs {pi['point']['f1b_g2_inst']:.5f}")
        # Lemma 1 (direct)
        a = pc.pulse_g2(10.0, 1.0, 2.0, 0.5, 1.0, 1.0, 0.1, 0.1, 12.4, 1.0)["g2"]
        b = pc.pulse_g2(10.0, 1.0, 2.0, 0.5, 1.0, 0.37, 0.037, 0.1, 12.4, 1.0)["g2"]
        check("pulse_counting [direct]: Lemma 1, common scaling of t_X, t_XX leaves g2 unchanged to 1e-12", abs(a - b) < 1e-12, f"{abs(a - b):.1e}")
        # two-level limit through the endpoint
        ps = ok200("/api/explore/pulse_counting?pump_ratio=0&eps=0&r=10000&tau_on=0.1&period=50&gamma_X=1")
        lam = 1.0 * 0.1
        ref2 = (2 * lam + lam ** 2) / (1 + lam) ** 2
        check("pulse_counting two-level (pump_ratio 0, eps 0, k 0): r = 1e4, tau_on = 0.1 gives (2 G t + (G t)^2)/(1 + G t)^2 rtol 2e-3",
              abs(ps["point"]["g2"] / ref2 - 1) < 2e-3, f"{ps['point']['g2']:.5f} vs {ref2:.5f}")
        pw = ok200("/api/explore/pulse_counting?pump_ratio=0&eps=0&r=0.01&tau_on=0.01&period=50&gamma_X=1")
        check("pulse_counting two-level weak pump: r = 0.01, tau_on = 0.01 gives G tau / 3 to rtol 3e-2",
              abs(pw["point"]["g2"] / (0.01 / 3.0) - 1) < 3e-2, f"{pw['point']['g2']:.6f} vs {0.01 / 3:.6f}")
        pwd = pc.pulse_g2(1e-3, 1.0, 2.0, 0.0, 0.0, 1.0, 0.0, 1e-3, 200.0, pump_ratio=0.0)["g2"]
        check("pulse_counting [direct]: weak short pump (tau_on = 1e-3) G tau / 3 to rtol 1e-2", abs(pwd / (1e-3 / 3) - 1) < 1e-2)
        # deterministic cycle, ideal
        pdet = ok200("/api/explore/pulse_counting?eps=0&period=12.5&eta_load=1")
        D = pdet["deterministic"]
        check("pulse_counting deterministic: eta_load 1, t_XX 0, period 12.5 gives mean -> 1 (0.99999627), g2 = 0 exactly, blocked -> 0",
              abs(D["mean_counts"] - 0.99999627) < 5e-8 and D["g2"] == 0.0 and D["blocked_load_probability"] < 1e-9 and D["one_pair_valid"] is True,
              f"mean {D['mean_counts']:.8f}")
        # adjacent
        pf = ok200("/api/explore/pulse_counting?period=50&r=10")
        check("pulse_counting: dark window >> relaxation gives g2_adj == g2 (adjacent_peak_factor = 1 to 1e-12)",
              abs(pf["point"]["adjacent_peak_factor"] - 1.0) < 1e-12 and abs(pf["point"]["g2_adj"] / pf["point"]["g2"] - 1) < 1e-12)
        check("pulse_counting: p_period sums to 1 to 1e-12 and converged is true", abs(sum(P0["p_period"]) - 1) < 1e-12 and P0["converged"] is True)
        check("pulse_counting: g2 differs from g2_adj when carry-over is large (default period, k = 0): normalisations are distinct fields",
              P0["adjacent_peak_factor"] != 1.0 and P0["g2"] != P0["g2_adj"])
        # gate
        pg = ok200("/api/explore/pulse_counting?gate=auto&gamma_X=2")
        ptg = pc.pulse_g2(10.0, 2.0, 4.0, 0.0, 0.0, 1.0, 0.1, 0.1, 12.4, 1.0, split=True, gate_ns=0.1 + 5.0 / 2.0, adjacent=True)
        check("pulse_counting: gate auto = tau_on + 5/gamma_X ([A]) equals the direct gated call", pg["point"]["g2"] == ptg["g2"] and pg["inputs"]["gate_ns"] == 0.1 + 2.5)
        pT = ok200("/api/explore/pulse_counting?T=230")
        ptT = pc.pulse_g2(10.0, 1.0, 2.0, kX, kXX, 1.0, 0.1, 0.1, 12.4, 1.0, split=True, adjacent=True)
        check("pulse_counting: T = 230 K point equals direct pulse_g2 with the same escape rates", pT["point"]["g2"] == ptT["g2"] and pT["k_X"] == kX)
        check("pulse_counting: pinned caveat text is exact and the background note is present",
              p["caveat"] == ea.CAVEAT_PULSED and "g2 is normalised to the long-delay peak" in p["caveat"] and "rho is a device.py addition" in p["caveat_extra"])
        check("pulse_counting: gate line is the /api/gates value", p["gate_g2"]["value"] == campaigns.gates()["g2_ceiling"]["value"])
        codes = [get(u)[0] for u in ("/api/explore/pulse_counting?period=0.05&tau_on=0.1", "/api/explore/pulse_counting?gate=maybe",
                                     "/api/explore/pulse_counting?tau_on=5", "/api/explore/pulse_counting?eta_load=2")]
        check("pulse_counting: bad inputs answer 400 (period <= tau_on, bad gate, out of range)", codes == [400] * 4, str(codes))

        # ================================================================ 6. lindblad
        worst = 0.0
        for qs in ("", "&r=20&eps=0.3", "&T=230&r=5", "&pump_ratio=0.4&gamma_X=2.5&eps=0.9"):
            li = ok200("/api/explore/lindblad?panel=incoherent" + qs)
            worst = max(worst, li["readouts"]["max_rel_diff"])
        check("lindblad (a1): incoherent 3-level g2_tau == cw_g2.g2_cw at 4 points (rel <= 1e-10)", worst <= 1e-10, f"{worst:.1e}")
        li = ok200("/api/explore/lindblad?panel=incoherent")
        s = lb.build_system(levels=3, r_ns=0.5, gamma_X_ns=1.0, gamma_XX_ns=2.0, k_X=0.0, k_XX=0.0, pump_ratio=1.0)
        g2d, _, _ = lb.g2_tau(s, np.linspace(0, 10.0, 201), weights=[1.0, 0.1])
        check("lindblad: incoherent curve equals direct g2_tau (== no tolerance)", li["curves"]["g2_lindblad"] == g2d.tolist())
        # (b) resonance fluorescence
        rb = ok200("/api/explore/lindblad?panel=rabi&omega=3&gamma=1&deph=0.5")
        g2c = 1.0 / 2 + 0.5 / 2
        ref_ee = 3 ** 2 / (2 * (1.0 * g2c + 3 ** 2))
        check("lindblad (b1): rho_ee = Omega^2/(2(gamma_1 gamma_2 + Omega^2)) = 0.46153846 at Omega 3, gamma 1, c 0.5 to 1e-10",
              abs(rb["readouts"]["rho_ee"] - ref_ee) < 1e-10 and abs(ref_ee - 0.46153846) < 1e-8, f"{rb['readouts']['rho_ee']:.8f}")
        rb0 = ok200("/api/explore/lindblad?panel=rabi&omega=3&gamma=1&deph=0")
        t = np.array(rb0["curves"]["tau"])
        mu = math.sqrt(9 - 1 / 16)
        cwc = 1 - np.exp(-0.75 * t) * (np.cos(mu * t) + 0.75 / mu * np.sin(mu * t))
        check("lindblad (b3): resonance-fluorescence g2(tau) == Carmichael-Walls closed form to 1e-8, g2(0) = 0",
              np.max(np.abs(np.array(rb0["curves"]["g2"]) - cwc)) < 1e-8 and abs(rb0["readouts"]["g2_0"]) < 1e-9)
        check("lindblad: Rabi ringing at Omega 3, gamma 1 peaks at 1.4545 near 1.06 ns (the code's value; the brief says 1.40)",
              abs(rb0["readouts"]["g2_max"] - 1.4545) < 1e-3 and abs(rb0["readouts"]["tau_at_max"] - 1.06) < 0.02,
              f"{rb0['readouts']['g2_max']:.4f} at {rb0['readouts']['tau_at_max']:.2f}")
        rbp = ok200("/api/explore/lindblad?panel=rabi&omega=20&gamma=0.2&deph=0")
        check("lindblad (b1): saturation limit rho_ee -> 1/2 for Omega >> gamma", abs(rbp["readouts"]["rho_ee"] - 0.5) < 5e-3)
        # (d) HOM
        hm = ok200("/api/explore/lindblad?panel=hom&gamma=1&gstar=10")
        check("lindblad (d1): HOM I = gamma/(gamma + gamma*) = 1/11 = 0.0909091 at gamma* = 10 to 1e-9",
              abs(hm["readouts"]["indistinguishability"] - 1 / 11) < 1e-9, f"{hm['readouts']['indistinguishability']:.7f}")
        gs = np.array(hm["curves"]["gstar"])
        check("lindblad (d1): the swept I(gamma*) curve equals gamma/(gamma + gamma*) to 1e-9 and falls monotonically",
              np.max(np.abs(np.array(hm["curves"]["indistinguishability"]) - 1.0 / (1.0 + gs))) < 1e-9
              and all(b < a for a, b in zip(hm["curves"]["indistinguishability"], hm["curves"]["indistinguishability"][1:])))
        hm0 = ok200("/api/explore/lindblad?panel=hom&gamma=2&gstar=0")
        check("lindblad (d1): gamma* = 0 is perfectly indistinguishable", abs(hm0["readouts"]["indistinguishability"] - 1.0) < 1e-9)
        # (e2) filter
        fl = ok200("/api/explore/lindblad?panel=filter&r=0.5&gamma_X=1&L=7.6&w=100")
        Dd = 0.5 + 1.0
        L_, w_ = 7.6, 100.0
        refg = 2 * Dd * (w_ + L_) / ((Dd + w_) * (3 * w_ + L_))
        check("lindblad (e2): filtered two-level g2_f(0) == 2 D (w+L)/((D+w)(3w+L)) (weak coupling) to rtol 2e-3",
              abs(fl["readouts"]["g2_filtered_0"] / refg - 1) < 2e-3, f"{fl['readouts']['g2_filtered_0']:.6f} vs {refg:.6f}")
        fl2 = ok200("/api/explore/lindblad?panel=filter&r=0.5&gamma_X=1&L=100&w=100")
        check("lindblad (e2): at w = L the filtered g2_f(0) = dip/(dip + w), dip = r + gamma_X + k_X, to rtol 2e-3",
              abs(fl2["readouts"]["g2_filtered_0"] / (Dd / (Dd + 100.0)) - 1) < 2e-3, f"{fl2['readouts']['g2_filtered_0']:.6f}")
        check("lindblad: a line narrower than the lifetime limit is reported as lifetime-limited",
              ok200("/api/explore/lindblad?panel=filter&r=0.5&L=0.1&w=2")["inputs"]["lifetime_limited"] is True and fl["inputs"]["lifetime_limited"] is False)
        # (f1) every state physical
        phys = []
        for pn in ("incoherent", "rabi", "filter"):
            phys.append(ok200(f"/api/explore/lindblad?panel={pn}")["state"])
        phys.append(ok200("/api/explore/lindblad?panel=incoherent&T=230&r=5")["state"])
        check("lindblad (f1): every returned state has |Tr - 1| <= 1e-12, min eig >= -1e-12, Hermitian (<= 1e-12)",
              all(x["ok"] for x in phys), " ".join(f"{x['trace_error']:.0e}" for x in phys))
        check("lindblad: pinned caveat is exact on every panel",
              all(ok200(f"/api/explore/lindblad?panel={pn}")["caveat"] == ea.CAVEAT_LINDBLAD for pn in ("incoherent", "rabi", "hom", "filter"))
              and "Markovian baths, Lorentzian ZPL" in ea.CAVEAT_LINDBLAD)
        codes = [get(u)[0] for u in ("/api/explore/lindblad?panel=nope", "/api/explore/lindblad?panel=rabi&omega=0", "/api/explore/lindblad?panel=hom&gstar=1000")]
        check("lindblad: bad panel or out-of-range input answers 400", codes == [400] * 3, str(codes))
        check("lindblad: every readout group carries a tag", all(ok200(f"/api/explore/lindblad?panel={pn}")["tags"]["readouts"] in ("V", "DR", "E", "A")
                                                                    for pn in ("incoherent", "rabi", "hom", "filter")))
        # ================================================================ review finding 2: widest-tag rule, no hard-coded tags
        rank = {"V": 0, "DR": 1, "E": 2, "A": 3}

        def widest(*t):
            return max((x for x in t if x), key=rank.get)
        cwt = {s_["id"]: s_["tag"] for s_ in ctl["controls"]["cw_g2"]}
        pct = {s_["id"]: s_["tag"] for s_ in ctl["controls"]["pulse_counting"]}
        d_nt = ok200("/api/explore/cw_g2")
        d_t = ok200("/api/explore/cw_g2?T=230")
        want_lim = widest(cwt["gamma_X"], cwt["eps"], cwt["pump_ratio"])
        check("review 2: cw_g2 'limits' tag is the widest of the controls it depends on (gamma_X [E], eps [A], p [A] -> A), "
              "not a hard-coded DR; same with the DR escape rates in",
              d_nt["tags"]["limits"] == want_lim == "A" and d_t["tags"]["limits"] == widest(want_lim, "DR") == "A",
              f"{d_nt['tags']['limits']} {d_t['tags']['limits']}")
        pp = ok200("/api/explore/pulse_counting")
        check("review 2: pulse_counting f1b tag is the widest of DR and its inputs r, tau_on, eps (all [A] -> A)",
              pp["tags"]["f1b"] == widest("DR", pct["r"], pct["tau_on"], pct["eps"]) == "A", pp["tags"]["f1b"])
        src = "".join((ROOT / "fsim_studio/web/js/workspaces/explorers" / f).read_text(encoding="utf-8") for f in ("pulsed.js", "cwg2.js"))
        import re as _re
        check("review 2: pulsed.js and cwg2.js carry no hard-coded [DR]/[E]/[V] tag text in hover templates or annotations",
              not _re.search(r"(?:hovertemplate|text)\s*:[^\n]*\[(?:DR|E|V)\]", src), "")
        # r -> inf limit gamma_X/(eps gamma_XX): independent of escape rates and of p (derivation in studio_support)
        dhi = ok200("/api/explore/cw_g2?T=230&eps=0.2&gamma_X=2&pump_ratio=0.5")
        dz = ok200("/api/explore/cw_g2?eps=0")
        dp0 = ok200("/api/explore/cw_g2?pump_ratio=0")
        check("review 4: r->inf limit = gamma_X/(eps gamma_XX) = 1/(2 eps) (2.5 at eps .2) with escape and p<1; None for eps = 0 or p = 0",
              abs(dhi["vs_r"]["high_limit"] - 2.0 / (0.2 * 4.0)) < 1e-12 and dz["vs_r"]["high_limit"] is None and dp0["vs_r"]["high_limit"] is None,
              str(dhi["vs_r"]["high_limit"]))
        check("review 4: the Flask layer holds no cw_g2 limit formula (api_explore_a imports it from fsim_core.studio_support)",
              "g2_cw_high_pump_limit" in (ROOT / "fsim_studio/api_explore_a.py").read_text(encoding="utf-8")
              and "gamma_X / (eps * gamma_XX)" not in (ROOT / "fsim_studio/api_explore_a.py").read_text(encoding="utf-8").replace("# r -> inf limit gamma_X/(eps gamma_XX)", ""))
    finally:
        try:
            app.config["JOB_MANAGER"].terminate()
        except Exception:
            pass
        tmp.cleanup()

    n_ok = sum(1 for _, ok, _ in CHECKS if ok)
    print(f"\n{n_ok}/{len(CHECKS)} studio explorer A checks passed")
    return 0 if n_ok == len(CHECKS) else 1


if __name__ == "__main__":
    sys.exit(main())
