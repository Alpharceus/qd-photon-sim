"""Verification for deterministic one-pair-per-cycle photon counting.

The Deshpande et al., Applied Physics Letters 105, 141109 (2014), DOI
10.1063/1.4897640 abstract reports tau=1.3 +/- 0.3 ns at 300 K and a maximum
excitation repetition of 200 MHz [V abstract-only].  Its lifetime is used
below only as a source-transcription and analytic-limit fixture; the other
checks are numerical verification [A] budgets, not held-out predictions.
"""
import argparse
import hashlib
import json
import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp
from scipy.linalg import expm

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.device import DeviceDesign, evaluate
from fsim_core import pulse_counting

pulse_g2 = pulse_counting.pulse_g2

BASELINE = ROOT / "verify" / "data" / "deterministic_cycle_legacy_baseline.json"
SOURCE = ROOT / "fsim_core" / "pulse_counting.py"


def _canonical(value):
    """JSON-safe, bit-preserving recursive representation for regression data."""
    if isinstance(value, np.ndarray):
        return {"array": {"dtype": value.dtype.str, "shape": list(value.shape),
                           "bytes": value.tobytes().hex()}}
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float):
        return "nan" if math.isnan(value) else ("inf" if value == math.inf else
               ("-inf" if value == -math.inf else {"float_hex": value.hex()}))
    if isinstance(value, dict):
        return {str(k): _canonical(v) for k, v in sorted(value.items())}
    if isinstance(value, (list, tuple)):
        return [_canonical(v) for v in value]
    if isinstance(value, (str, int, bool)) or value is None:
        return value
    return repr(value)


def _legacy_cases():
    # Covers zero/full/restricted gates, split, escape, and the rectangular
    # path without changing its call order or defaults.
    return {
        "ungated": pulse_g2(2.5, 1.0, 2.0, 0.3, 0.6, 0.8, 0.15, 0.1, 12.4),
        "gate_zero": pulse_g2(2.5, 1.0, 2.0, 0.3, 0.6, 0.8, 0.15, 0.1, 12.4,
                              gate_ns=0.0),
        "gate_restricted_split": pulse_g2(2.5, 1.0, 2.0, 0.3, 0.6, 0.8, 0.15,
                                            0.1, 12.4, split=True, gate_ns=0.35),
        "gate_full_split": pulse_g2(0.6, 1.0, 2.0, 0.0, 0.0, 1.0, 0.05,
                                      0.2, 4.8, split=True, gate_ns=5.0),
    }


def _card_cases():
    """Evaluate both shipped InP cards under both tau_cap density conventions.

    Raw (non-canonical) dict, shared by the capture path and the live
    regression comparison so both go through the same fsim_core.device.evaluate
    calls in the same order.
    """
    cards = {}
    for name in ("edge-inp-gainp-design.yaml", "edge-inp-gaasp-design.yaml"):
        design = DeviceDesign.load(ROOT / "cards" / name)
        # Both shipped retention density conventions are captured explicitly.
        for density_mode in (False, True):
            design.ret.tau_cap_scales_with_density = density_mode
            cards[f"{name}:tau_cap_density={density_mode}"] = evaluate(design)
    return cards


def _capture_legacy(force=False):
    current_hash = hashlib.sha256(SOURCE.read_bytes()).hexdigest()
    if BASELINE.exists() and not force:
        old = json.loads(BASELINE.read_text(encoding="utf-8"))
        old_hash = old.get("source_sha256")
        if old_hash != current_hash:
            print("refusing to overwrite pre-edit baseline: recorded source_sha256 "
                  f"{old_hash} does not match current fsim_core/pulse_counting.py hash "
                  f"{current_hash}; pass --force-recapture to override")
            return 1
    payload = {"purpose": "pre-edit bitwise regression capture; not a literature oracle",
               "source_sha256": current_hash,
               "pulse_g2": _canonical(_legacy_cases()), "cards": _canonical(_card_cases())}
    BASELINE.parent.mkdir(exist_ok=True)
    BASELINE.write_text(json.dumps(payload, sort_keys=True, indent=2) + "\n", encoding="utf-8")
    print("captured deterministic-cycle legacy baseline")
    return 0


def _independent(case):
    """Own load map/generator and solve_ivp moment propagation [A]."""
    gx, gxx, kx, kxx, tx, txx, period, eta, gate = case
    M = np.array([[0.0, gx + kx, 0.0], [0.0, -(gx + kx), gxx + kxx],
                  [0.0, 0.0, -(gxx + kxx)]])
    L1 = np.array([[0.0, 0.0, 0.0], [1.0, 0.0, 0.0], [0.0, 1.0, 1.0]])
    L = (1.0 - eta) * np.eye(3) + eta * L1
    phi = expm(M * period) @ L
    p = np.array([1.0, 0.0, 0.0])
    for _ in range(10000):
        q = phi @ p
        if np.max(np.abs(q - p)) < 1e-12:
            p = q
            break
        p = q
    J = np.zeros((3, 3)); J[0, 1] = tx * gx; J[1, 2] = txx * gxx
    def rhs(counting):
        def f(_t, v):
            pv, m1, m2 = v[:3], v[3:6], v[6:]
            jj = J if counting else np.zeros((3, 3))
            return np.r_[M @ pv, M @ m1 + jj @ pv, M @ m2 + 2.0 * jj @ m1]
        return f
    v = np.zeros(9); v[:3] = L @ p
    if gate > 0:
        v = solve_ivp(rhs(True), (0, gate), v, rtol=1e-11, atol=1e-13).y[:, -1]
    if gate < period:
        v = solve_ivp(rhs(False), (gate, period), v, rtol=1e-11, atol=1e-13).y[:, -1]
    return v


def _close(a, b):
    return abs(a - b) <= 1e-9 + 1e-7 * abs(b)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--capture-legacy", action="store_true")
    parser.add_argument("--force-recapture", action="store_true",
                         help="overwrite the baseline even if the recorded pre-edit "
                              "source_sha256 no longer matches fsim_core/pulse_counting.py")
    args = parser.parse_args()
    if args.capture_legacy:
        return _capture_legacy(force=args.force_recapture)
    checks = []
    def ok(label, condition):
        checks.append(bool(condition)); print(("PASS" if condition else "FAIL") + " " + label)

    # Numerical verification: independent solve_ivp for short/long lifetime and escape.
    for n, case in enumerate(((1.0, 2.0, 0.2, 0.4, 0.8, 0.1, 5.0, 0.7, 0.6),
                              (0.08, 0.16, 0.03, 0.06, 1.0, 0.3, 30.0, 1.0, 30.0))):
        got = pulse_counting.deterministic_cycle_g2(*case[:7], eta_load=case[7], gate_ns=case[8])
        ref = _independent(case)
        ok(f"numerical verification solve_ivp case {n+1}",
           _close(got["mean_counts"], float(ref[3:6].sum())) and
           _close(got["mean_factorial2"], float(ref[6:9].sum())))

    # Source transcription, explicitly abstract-only rather than a performance target:
    # the literal published citation and numbers must appear verbatim in the module
    # source (fsim_core/pulse_counting.py), not just be re-asserted here.
    source_text = SOURCE.read_text(encoding="utf-8")
    ok("source transcription Deshpande APL 2014 abstract-only tau=1.3 ns T=300 K",
       "Applied Physics Letters 105, 141109 (2014)" in source_text and
       "10.1063/1.4897640" in source_text and
       "tau=1.3 +/- 0.3 ns" in source_text and
       "300 K" in source_text and
       "[V abstract-only]" in source_text)
    gx = 1.0 / 1.3
    iso = pulse_counting.deterministic_cycle_g2(gx, 2 * gx, 0, 0, 0.7, 0.0, 100.0, gate_ns=0.7)
    target = 0.7 * gx / gx * (1 - math.exp(-gx * 0.7))
    ok("numerical verification isolated one-X Bernoulli limit", _close(iso["mean_counts"], target)
       and iso["mean_factorial2"] < 1e-12)
    kx = 0.05
    iso_esc = pulse_counting.deterministic_cycle_g2(gx, 2 * gx, kx, 0, 0.7, 0.0, 100.0, gate_ns=0.7)
    target_esc = 0.7 * gx / (gx + kx) * (1 - math.exp(-(gx + kx) * 0.7))
    ok("numerical verification isolated one-X Bernoulli limit with escape",
       _close(iso_esc["mean_counts"], target_esc) and iso_esc["mean_factorial2"] < 1e-12)
    base = pulse_counting.deterministic_cycle_g2(1, 2, .1, .2, 1, .3, 6, eta_load=.6)
    loss = pulse_counting.deterministic_cycle_g2(1, 2, .1, .2, .4, .12, 6, eta_load=.6)
    thin = pulse_counting.deterministic_cycle_g2(1, 2, .1, .2, .5, .15, 6, eta_load=.6)
    ok("numerical verification loss monotonicity", loss["mean_counts"] < base["mean_counts"])
    ok("numerical verification common thinning invariance", abs(base["g2"] - thin["g2"]) < 1e-12)
    split = pulse_counting.deterministic_cycle_g2(1, 2, .1, .2, .8, .1, 6, eta_load=.6, split=True)
    ok("numerical verification split means sum", abs(split["mean_counts"] - split["mean_counts_x"] - split["mean_counts_xx"]) < 1e-12)
    clamp = pulse_counting.deterministic_cycle_g2(1, 2, .1, .2, .8, .1, 6, gate_ns=20)
    full = pulse_counting.deterministic_cycle_g2(1, 2, .1, .2, .8, .1, 6)
    ok("numerical verification gate clamp and full-period state", _close(clamp["mean_counts"], full["mean_counts"])
       and np.allclose(clamp["p_period"], full["p_period"]))
    zero = pulse_counting.deterministic_cycle_g2(1, 2, 0, 0, 1, 0, 4, eta_load=0)
    dark = pulse_counting.deterministic_cycle_g2(0, 0, 0, 0, 1, 1, 4)
    ok("numerical verification eta_load=0", zero["mean_counts"] == 0 and math.isnan(zero["g2"]))
    ok("numerical verification dark source invalid reason", dark["invalid_reason"] == "no_detected_radiative_photons")
    residual = pulse_counting.deterministic_cycle_g2(.02, .04, 0, 0, 1, 1, .1)
    ok("numerical verification residual-X cascade has factorial pairs", residual["mean_factorial2"] > 0)
    slow = pulse_counting.deterministic_cycle_g2(.001, .002, 0, 0, 1, 1, .1)
    ok("numerical verification blocked residual pair invalidates exact-one-pair claim",
       slow["blocked_load_probability"] > 0 and not slow["one_pair_valid"])
    invalid = [(-1, 2, 0, 0, 1, 0, 1), (1, 2, 0, 0, 1.1, 0, 1), (1, 2, 0, 0, 1, 0, 0)]
    ok("numerical verification invalid inputs", all(_raises(x) for x in invalid))

    if not BASELINE.exists():
        ok("regression capture baseline exists", False)
    else:
        old = json.loads(BASELINE.read_text(encoding="utf-8"))
        current = {"pulse_g2": _canonical(_legacy_cases())}
        # Genuinely re-evaluate both cards through fsim_core.device.evaluate
        # (not just re-read the stored snapshot) so the "cards" half of the
        # comparison below is live evidence, not a comparison of the baseline
        # against itself.
        current["cards"] = _canonical(_card_cases())
        ok("regression capture legacy rectangular pulse and card outputs bit-identical",
           current["pulse_g2"] == old["pulse_g2"] and current["cards"] == old["cards"])
    print(f"{sum(checks)}/{len(checks)} deterministic_cycle_counting checks passed")
    return 0 if all(checks) else 1


def _raises(args):
    try:
        pulse_counting.deterministic_cycle_g2(*args)
    except ValueError:
        return True
    return False


if __name__ == "__main__":
    raise SystemExit(main())
