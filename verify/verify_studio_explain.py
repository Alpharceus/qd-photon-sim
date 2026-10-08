"""FSIM Studio explain endpoints (fsim_studio/api_explain.py): exact equality with
direct fsim_core calls and with the committed out/phase0 files.

The endpoints port the Streamlit app.py explainers; every number must equal the
value the same fsim_core call returns here (== on floats, no tolerance), and
the V-a bundle must equal the CSV/JSON text in out/phase0.

Run: python verify/verify_studio_explain.py   (exit code 0 iff all pass)
"""
import csv
import json
import os
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok), detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""), flush=True)


def main():
    tmp = tempfile.TemporaryDirectory(prefix="fsim_studio_explain_")
    os.environ["FSIM_STUDIO_CACHE"] = tmp.name
    import numpy as np

    from fsim_core.card import load_card
    from fsim_core.integrator import g2_of_T
    from fsim_core.loading import f1b_g2, gamma_eff, loading_probs
    from fsim_core.spectral import epsilon, gamma_of_T, lorentzian
    from fsim_studio.server import create_app

    app = create_app()
    c = app.test_client()
    jm = app.config["JOB_MANAGER"]
    try:
        card = load_card(ROOT / "cards" / "chatzarakis.yaml")
        fit = json.loads((ROOT / "out/phase0/fit_params.json").read_text(encoding="utf-8"))
        p = fit["params"]
        rows = sorted(card.datasets["g2_vs_T"].rows, key=lambda r: r["T"])
        Ts = np.array([r["T"] for r in rows])
        w_of_T = lambda T: float(np.interp(T, Ts, np.array([r["w"] for r in rows])))
        dx_of_T = lambda T: float(np.interp(T, Ts, np.array([r["dx"] for r in rows])))

        # 1-4 spectral, published window, PL
        for T in (78.0, 230.0):
            r = c.get(f"/api/explain/spectral?T={T}&pub=1")
            d = r.get_json()
            g = float(gamma_of_T(T, p["gamma0"], p["a_ac"], p["b_lo"], p["E_lo"]))
            w, dx = w_of_T(T), dx_of_T(T)
            s = epsilon(p["delta_xx"], g, g, w=w, dx=dx)
            check(f"spectral T={T:g} pub: gamma, w, dx equal direct calls",
                  r.status_code == 200 and d["gamma"] == g and d["w"] == w and d["dx"] == dx)
            check(f"spectral T={T:g} pub: eps, t_x, g2_0 equal epsilon()/f1b_g2()",
                  d["eps"] == s.eps and d["t_x"] == s.t_x and d["g2_dot"] == float(f1b_g2(p["mu"], s.eps)),
                  f"eps {d['eps']}")
        # 5-6 spectral, manual window, EL, mu = 0 (F1)
        d = c.get("/api/explain/spectral?T=200&pub=0&w=2.5&dx=0.4&mode=EL&dg_inj=3&I_ratio=2&mu=0").get_json()
        g = float(gamma_of_T(200.0, p["gamma0"], p["a_ac"], p["b_lo"], p["E_lo"]))
        ge = float(gamma_eff(g, 3.0, 2.0, 1.0))
        s = epsilon(p["delta_xx"], ge, ge, w=2.5, dx=0.4)
        check("spectral EL: gamma_eff equals loading.gamma_eff", d["gamma_el"] == ge and d["gamma"] == ge)
        check("spectral mu=0: g2_0 is eps (F1)", d["g2_dot"] == s.eps and d["eps"] == s.eps)
        # 7 curves are lorentzian()/peak on the same grid
        x = np.linspace(-p["delta_xx"] - 6 * ge - 4, 6 * ge + 4, 1200)
        LX = lorentzian(x, 0.0, ge)
        check("spectral curves equal lorentzian()/peak", d["curves"]["x"] == x.tolist()
              and d["curves"]["X"] == (LX / LX.max()).tolist())
        # 8 range validation
        bad = c.get("/api/explain/spectral?T=999")
        check("spectral rejects T outside 4-320 K with 400", bad.status_code == 400 and "error" in bad.get_json())

        # 9-10 cascade
        for T in (150.0, 230.0):
            d = c.get(f"/api/explain/cascade?T={T}").get_json()
            pt = g2_of_T(T, {**p, "w": w_of_T(T), "dx": dx_of_T(T)})
            P0, P1, P2 = loading_probs(p["mu"])
            check(f"cascade T={T:g}: g2, rho, eps, gamma equal integrator.g2_of_T; P equal loading_probs",
                  d["g2"] == pt.g2 and d["rho"] == pt.rho and d["eps"] == pt.eps and d["gamma"] == pt.gamma
                  and (d["P0"], d["P1"], d["P2"]) == (P0, P1, P2) and d["g2_from_check"] == pt.g2,
                  f"g2 {d['g2']}")

        # 11-12 V-a bundle verbatim
        d = c.get("/api/explain/va_fit").get_json()
        with open(ROOT / "out/phase0/fit_data_points.csv", newline="", encoding="utf-8") as f:
            pts = list(csv.DictReader(f))
        ok = len(d["points"]) == len(pts) and all(
            d["points"][i]["T_K"] == float(r["T_K"]) and d["points"][i]["residual"] == float(r["residual"])
            and d["points"][i]["bound"] == r["bound"] for i, r in enumerate(pts))
        check("va_fit points equal out/phase0/fit_data_points.csv", ok, f"{len(pts)} rows")
        with open(ROOT / "out/phase0/fit_model_curve.csv", newline="", encoding="utf-8") as f:
            cur = list(csv.DictReader(f))
        check("va_fit curve equals fit_model_curve.csv; Tc equals fit_params.json",
              d["curve"]["g2_model"] == [float(r["g2_model"]) for r in cur] and d["Tc_K"] == fit["Tc_K"]
              and d["label"] == "calibration (C), not held-out")
    finally:
        jm.shutdown()

    n_pass = sum(1 for _, ok, _ in CHECKS if ok)
    print(f"\n{n_pass}/{len(CHECKS)} studio explain checks passed")
    try:
        tmp.cleanup()
    except OSError:
        pass
    return 0 if n_pass == len(CHECKS) else 1


if __name__ == "__main__":
    sys.exit(main())
