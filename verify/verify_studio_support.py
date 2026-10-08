"""Checks for fsim_core/studio_support.py (the closed forms the FSIM Studio shows next to the literature
cards and the CW explorer) against independent closed forms, published numbers and the committed out/ CSVs.

Every reference here is computed in this file (explicit exponentials, algebra) or read from the committed
CSV/YAML; nothing is checked against the function under test itself.

Run: python verify/verify_studio_support.py   (exit code 0 iff all pass)
"""
import csv
import json
import math
import sys
from pathlib import Path

import numpy as np
import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok), detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""), flush=True)


def read_csv(rel):
    with open(ROOT / rel, newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def main():
    from fsim_core import cw_g2
    from fsim_core import integrator
    from fsim_core import studio_support as ss
    from fsim_core.spectral import epsilon_narrow_filter

    fit = json.loads((ROOT / "out/phase0/fit_params.json").read_text(encoding="utf-8"))["params"]
    P = {k: fit[k] for k in ("a_esc", "E_a", "b_p", "E_b")}
    kB = 8.617333262e-2   # meV/K (CODATA), written out here independently of integrator.KB

    # ---------------------------------------------------------------- tau_anchored
    def S_indep(T):
        return 1.0 / (1.0 + P["a_esc"] * math.exp(-P["E_a"] / (kB * T)) + P["b_p"] * math.exp(-P["E_b"] / (kB * T)))

    card = yaml.safe_load((ROOT / "cards/chatzarakis.yaml").read_text(encoding="utf-8"))
    rows = sorted(card["data"]["tau_vs_T"]["points"], key=lambda r: r["T"])
    T0, tau0 = rows[0]["T"], rows[0]["tau_ns"]
    Ts = np.linspace(T0, rows[-1]["T"], 25)
    got = [float(ss.tau_anchored(T, T0, tau0, **P)) for T in Ts]
    want = [tau0 * S_indep(T) / S_indep(T0) for T in Ts]
    check("tau_anchored: tau(T0) = tau0 exactly, and tau(T) = tau0 S(T)/S(T0) with S written out explicitly (rel 1e-9)",
          ss.tau_anchored(T0, T0, tau0, **P) == tau0 and max(abs(a - b) / b for a, b in zip(got, want)) < 1e-9,
          f"max rel {max(abs(a - b) / b for a, b in zip(got, want)):.2e}")
    check("tau_anchored: the lifetime ratio is the retention ratio of integrator.retention at two temperatures",
          abs(float(ss.tau_anchored(250.0, 100.0, 1.0, **P)) - float(integrator.retention(250.0, **P) / integrator.retention(100.0, **P))) < 1e-12)

    # ---------------------------------------------------------------- inverted_background_g2 (F2 law inverted)
    oe = yaml.safe_load((ROOT / "cards/reischle.yaml").read_text(encoding="utf-8"))["data"]["g2_oe2008"]["points"]
    inv = [float(ss.inverted_background_g2(r["g2_b"], r["rho"])) for r in oe]
    check("inverted_background_g2: reproduces the published g2_corr of Reischle OE 2008 within 0.011 (0.04, 0.03)",
          all(abs(a - r["g2_corr"]) <= 0.011 for a, r in zip(inv, oe)), str([round(x, 4) for x in inv]))
    rt = [integrator.g2_from(g2s, r["rho"]) for g2s, r in zip(inv, oe)]
    check("inverted_background_g2: the F-series background law integrator.g2_from(g2_s, rho) returns the input g2_b (1e-12)",
          all(abs(a - r["g2_b"]) < 1e-12 for a, r in zip(rt, oe)))
    try:
        ss.inverted_background_g2(0.2, 0.0)
        zero_rho = False
    except ValueError:
        zero_rho = True
    check("inverted_background_g2: rho = 0 raises (no division by zero)", zero_rho)

    # ---------------------------------------------------------------- g2_cw_high_pump_limit
    gX, gXX, eps = 2.0, 4.0, 0.2
    ok_num, detail = True, []
    for k_X, k_XX, p in ((0.0, 0.0, 1.0), (9.16, 0.4, 0.5), (3.0, 2.0, 1.0)):
        val = float(cw_g2.g2_cw_zero(np.array([1e8]), gX, gXX, eps, k_X, k_XX, p)[0])
        lim = ss.g2_cw_high_pump_limit(gX, gXX, eps, p)
        ok_num &= lim == gX / (eps * gXX) and abs(val - lim) / lim < 1e-5
        detail.append(f"{val:.6f}")
    check("g2_cw_high_pump_limit = gamma_X/(eps gamma_XX) = 2.5, independent of escape rates and p, and matches the exact "
          "g2_cw_zero at r = 1e8/ns (rel 1e-5)", ok_num, " ".join(detail))
    check("g2_cw_high_pump_limit: None for eps = 0 or p = 0 (g2_dot(0) has no r -> infinity limit there)",
          ss.g2_cw_high_pump_limit(1.0, 2.0, 0.0, 1.0) is None and ss.g2_cw_high_pump_limit(1.0, 2.0, 0.1, 0.0) is None)
    check("g2_cw_high_pump_limit: the brief's anchor gamma_X/(eps gamma_XX) = 5.0 at eps = 0.1, gamma_XX = 2 gamma_X",
          abs(ss.g2_cw_high_pump_limit(1.0, 2.0, 0.1, 1.0) - 5.0) < 1e-12)

    # ---------------------------------------------------------------- implied_dip_recovery_time
    td = [float(ss.implied_dip_recovery_time(r["g2_raw"], r["g2_b"], 0.5)) for r in oe]
    anchor = 0.5 * (1 - 0.41) / (0.85 - 0.59)    # physics brief 1(c): 1.13 ns
    back = [float(cw_g2.dip_convolved_exp(1.0 - r["g2_b"], t, 0.5)) for t, r in zip(td, oe)]
    check("implied_dip_recovery_time: 1.13 ns for the 0.41 -> 0.15 pair (brief anchor) and 1.58 ns for 0.43 -> 0.25",
          abs(td[0] - anchor) < 1e-12 and abs(td[0] - 1.1346) < 1e-3 and abs(td[1] - 0.5 * 0.57 / 0.18) < 1e-12, str([round(x, 4) for x in td]))
    check("implied_dip_recovery_time: cw_g2.dip_convolved_exp(A, tau_d, 0.5) returns the published raw g2(0) (0.410, 0.430) to 1e-12",
          all(abs(b - r["g2_raw"]) < 1e-12 for b, r in zip(back, oe)), str([round(x, 6) for x in back]))
    try:
        ss.implied_dip_recovery_time(0.2, 0.9, 0.5)
        bad = False
    except ValueError:
        bad = True
    check("implied_dip_recovery_time: an intrinsic dip shallower than the raw dip has no finite tau_d (ValueError)", bad)

    # ---------------------------------------------------------------- committed CSVs the Library reads instead of recomputing
    reis = yaml.safe_load((ROOT / "cards/reischle.yaml").read_text(encoding="utf-8"))["data"]["g2_vs_ERR"]["points"]
    g2_100 = next(r["g2"] for r in reis if r["device"] == 1 and r["position"] == 2 and r["ERR_MHz"] == 100.0)
    win = read_csv("out/phase1/vc_rho_window.csv")
    check("out/phase1/vc_rho_window.csv: rho_required = sqrt(1 - g2) of the 100 MHz point (g2 = 0.37 -> 0.79373) on every row, 1e-12",
          len(win) > 5 and all(abs(float(r["rho_required"]) - math.sqrt(1.0 - g2_100)) < 1e-12 for r in win), f"{len(win)} rows")
    pz = yaml.safe_load((ROOT / "cards/qcap-piezo-variant.yaml").read_text(encoding="utf-8"))["params"]
    glo, ghi = pz["gamma_300_lo"]["value"], pz["gamma_300_hi"]["value"]
    mp = read_csv("out/phase3/map_rho_required_300K.csv")
    bad_rows = 0
    for r in mp:
        d = float(r["delta_meV"])
        for g, key in ((glo, "rho_req_gamma6.0"), (ghi, "rho_req_gamma7.0")):
            e = float(epsilon_narrow_filter(d, g, g))
            rr = math.sqrt(0.5 / (1.0 - e)) if e < 1.0 else float("inf")
            want_v = rr if rr <= 1.0 else float("nan")
            got_v = float(r[key])
            if not ((math.isnan(want_v) and math.isnan(got_v)) or abs(got_v - want_v) <= 1e-9 * max(1.0, abs(want_v))):
                bad_rows += 1
    check("out/phase3/map_rho_required_300K.csv: rho >= sqrt(0.5/(1 - eps_min(300 K))) per row and Gamma(300) edge, NaN where rho > 1",
          len(mp) == 61 and bad_rows == 0 and glo == 6.0 and ghi == 7.0, f"{len(mp)} rows, {bad_rows} mismatches")
    zc = read_csv("out/zhao/zhao_fit_comparison.csv")
    check("out/zhao/zhao_fit_comparison.csv: both verdicts FAIL and agree with sigma_ratio > 1 (1.69, 5.82); the Studio shows these words",
          [r["verdict"] for r in zc] == ["FAIL", "FAIL"] and all(float(r["sigma_ratio"]) > 1.0 for r in zc)
          and abs(float(zc[0]["sigma_ratio"]) - 1.69) < 0.01 and abs(float(zc[1]["sigma_ratio"]) - 5.82) < 0.01)

    n_ok = sum(1 for _, ok, _ in CHECKS if ok)
    print(f"\n{n_ok}/{len(CHECKS)} studio support checks passed")
    return 0 if n_ok == len(CHECKS) else 1


if __name__ == "__main__":
    sys.exit(main())
