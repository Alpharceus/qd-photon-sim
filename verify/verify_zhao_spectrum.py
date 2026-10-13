"""Verifier for the output-flux noise SPECTRUM estimator of fsim_core.sde (devlib-r2-pha2, Zhao, r2-physics-a Q6).

What is checked, against numbers the code did not produce:
  * the Welch Fano estimator on SYNTHETIC streams with closed-form spectra (Poisson: F = 1; an OU-modulated Poisson
    stream: F(f) = 1 + mu m^2 (1 - a^2)/(1 - 2 a cos(2 pi f dt) + a^2), exact for the discrete AR(1) modulation) [DR];
  * Parseval: the frequency average of F equals the time-domain binned-count Fano factor (same data);
  * agreement with the EXISTING sde.output_counting_stats windowed Fano on the same seed (identical d_loss series);
  * the binomial-detection identity F_det = 1 - eta (1 - F) and the dB conversion;
  * the Yamamoto-Machida-Nilsson asymptote [E]: closed-form values 0.31 (-5.1 dB) at 4 I_th, F_pump 0.08, eta_q 1, and
    1.0 for a Poisson pump; the deterministic eta_q of the model (diagnostic, reported, not a pin);
  * qualitative physics: pump regulation lowers the low-frequency output noise (quiet < normal, > 3 sigma).
NOT checked (and NOT registered as a pass): the paper's measured squeezing (3.1 dB card vs 0.9 +- 0.1 dB digest: the
conflict cannot be resolved from the workspace sources, ledger `unavailable:` zhao_squeezing_dB) and Zhao's g2 values
(zhao_normal / zhao_quiet stay expect-fail in the ledger).

Sequential, OMP_NUM_THREADS-friendly (two SDE runs, ~1 min).  Exit 0 iff all pass; "N/N zhao_spectrum checks passed".
"""
import sys
import time
from pathlib import Path

import numpy as np
from scipy.signal import lfilter

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core import sde  # noqa: E402
from fsim_core.sde import (QDLaserParams, differential_efficiency, find_threshold,  # noqa: E402
                           output_counting_stats, output_noise_spectrum, welch_fano, yamamoto_fano)

RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append(bool(ok))
    print(("PASS " if ok else "FAIL ") + name + (f"  [{detail}]" if detail else ""), flush=True)


rng = np.random.default_rng(12345)

# ---------------------------------------------------------------- 1. estimator on synthetic streams (closed forms)
n_runs, n_bins, n_seg = 8, 20000, 256
mu0 = 50.0
pois = rng.poisson(mu0, size=(n_runs, n_bins)).astype(float)
F_run, _ = welch_fano(pois, n_seg)
F_mean = float(F_run[:, 1:-1].mean())
check("Welch Fano of a Poisson stream is 1 (normalisation, 3 %)", abs(F_mean - 1.0) < 0.03, f"F = {F_mean:.4f}")

dt_bin, tau_c, m, mu1 = 0.01, 0.5, 0.3, 0.2    # ns, ns, modulation depth, mean counts per bin (R = 20 / ns)
a = np.exp(-dt_bin / tau_c)
x = lfilter([np.sqrt(1 - a * a)], [1.0, -a], rng.standard_normal((n_runs, n_bins)), axis=1)   # unit-variance AR(1)
mod = rng.poisson(np.maximum(mu1 * (1.0 + m * x), 0.0)).astype(float)
F_run, _ = welch_fano(mod, n_seg)
Fm = F_run.mean(axis=0)
f_idx = np.arange(n_seg // 2 + 1)
F_cf = 1.0 + mu1 * m**2 * (1 - a * a) / (1 - 2 * a * np.cos(2 * np.pi * f_idx / n_seg) + a * a)
# a Hann-windowed periodogram returns the spectrum convolved with the window's power kernel: weights (1/6, 2/3, 1/6)
# at bin offsets (-1, 0, +1) for the periodic Hann window (even extension at f = 0).  The closed form is smoothed with it.
Fpad = np.concatenate([[F_cf[1]], F_cf, [F_cf[-2]]])
F_cf_s = (2.0 / 3.0) * F_cf + (1.0 / 6.0) * (Fpad[:-2] + Fpad[2:])
sel = slice(1, n_seg // 2)
rel = np.abs(Fm[sel] / F_cf_s[sel] - 1.0)
check("OU-modulated Poisson: Welch F(f) matches the closed-form spectrum (Hann-kernel smoothed; median 5 %, max 20 %)",
      np.median(rel) < 0.05 and rel.max() < 0.20, f"median {np.median(rel):.3f} max {rel.max():.3f}")
check("OU-modulated Poisson: lowest-bin excess F(1/L) - 1 within 10 % of the Hann-smoothed closed form",
      abs((Fm[1] - 1.0) / (F_cf_s[1] - 1.0) - 1.0) < 0.10, f"{Fm[1]:.3f} vs {F_cf_s[1]:.3f} (unsmoothed {F_cf[1]:.3f})")

# ---------------------------------------------------------------- 2. SDE: quiet vs normal, sequential runs
p = QDLaserParams()
I_th = find_threshold(p)
I_op = 4.0 * I_th
T_END, N_RUNS, SEED = 40.0e-9, 16, 7
out = {}
t0 = time.time()
for label, F_pump in (("normal", 1.0), ("quiet", 0.08)):
    out[label] = output_noise_spectrum(p, I_op, n_runs=N_RUNS, t_end=T_END, seed=SEED, F_pump=F_pump)
    print(f"  [{label}: {time.time() - t0:.0f} s]", flush=True)
    o = out[label]
    print(f"  {label}: window {o['window_ns']:.2f} ns, band {o['band_Hz'][0] / 1e9:.3f}-{o['band_Hz'][1] / 1e9:.3f} GHz, "
          f"F_band {o['F_band']:.3f} +- {o['F_band_se']:.3f}, F_low {o['F_low']:.3f} +- {o['F_low_se']:.3f}, "
          f"F_window {o['F_window']:.3f} +- {o['F_window_se']:.3f}, F_bin {o['F_bin_time']:.2f} / parseval {o['F_bin_parseval']:.2f}")

for label in ("normal", "quiet"):
    o = out[label]
    check(f"{label}: Parseval, frequency-average F vs time-domain binned-count Fano (5 %)",
          abs(o["F_bin_parseval"] / o["F_bin_time"] - 1.0) < 0.05,
          f"{o['F_bin_parseval']:.3f} vs {o['F_bin_time']:.3f}")

# the same seed gives the same d_loss in the existing windowed estimator: identical windows -> identical Fano factors
for label, F_pump in (("normal", 1.0), ("quiet", 0.08)):
    o = out[label]
    cs = output_counting_stats(p, I_op, o["window_ns"], n_runs=N_RUNS, t_end=T_END, dt=o["dt"], seed=SEED, F_pump=F_pump)
    check(f"{label}: spectrum module's windowed Fano == existing output_counting_stats on the same seed (1e-9)",
          abs(cs["fano"][0] - o["F_window"]) < 1e-9, f"{cs['fano'][0]:.6f} vs {o['F_window']:.6f}")

n_, q_ = out["normal"], out["quiet"]
sig = np.hypot(n_["F_band_se"], q_["F_band_se"])
check("pump regulation lowers the low-frequency output noise: F_band(quiet) < F_band(normal) by > 3 sigma",
      n_["F_band"] - q_["F_band"] > 3.0 * sig, f"{q_['F_band']:.3f} vs {n_['F_band']:.3f}, sigma {sig:.3f}")
for label in ("normal", "quiet"):
    o = out[label]
    s2 = np.hypot(o["F_low_se"], o["F_window_se"])
    check(f"{label}: lowest-bin F(f -> 0) consistent with the long-window Fano F_T within 3 sigma (same data)",
          abs(o["F_low"] - o["F_window"]) < 3.0 * s2,
          f"{o['F_low']:.3f} vs {o['F_window']:.3f}, sigma {s2:.3f}")

# ---------------------------------------------------------------- 3. detection identity, dB, asymptote
o = sde.output_noise_spectrum(p, I_op, n_runs=4, t_end=8.0e-9, seed=SEED, F_pump=0.08, eta_det=0.5)
check("detection: F_det = 1 - eta (1 - F) and squeezing_dB = 10 log10 F_det (eta 0.5)",
      abs(o["F_det_band"] - (1.0 - 0.5 * (1.0 - o["F_band"]))) < 1e-12
      and abs(o["squeezing_dB"] - 10.0 * np.log10(o["F_det_band"])) < 1e-12)
check("Yamamoto-Machida-Nilsson asymptote [E]: 4 I_th, F_pump 0.08, eta_q 1 -> 0.31 (-5.1 dB); Poisson pump -> 1.0",
      abs(yamamoto_fano(0.08, 1.0, 4.0) - 0.31) < 1e-12 and abs(10 * np.log10(0.31) + 5.08) < 0.05
      and yamamoto_fano(1.0, 0.9, 4.0) == 1.0)
eta_q = differential_efficiency(p, I_op)
check("model differential efficiency eta_q (deterministic L-I slope) is a physical fraction", 0.3 < eta_q <= 1.0,
      f"eta_q = {eta_q:.3f}")
for label, F_pump in (("normal", 1.0), ("quiet", 0.08)):
    asym = yamamoto_fano(F_pump, eta_q, 4.0)
    o = out[label]
    print(f"  DIAGNOSTIC {label}: model lowest-decade F {o['F_band']:.2f} (F_low {o['F_low']:.2f}) vs Yamamoto asymptote "
          f"{asym:.3f} ({10 * np.log10(asym):+.2f} dB); excess low-frequency noise = RO / partition noise (not a pin)")
print("  DIAGNOSTIC card conflict: 0.9 dB -> F_det 0.813 -> eta_total = 0.187/(1 - F(0)); the paper value is UNAVAILABLE "
      "(3.1 vs 0.9 dB unresolved), so no pass/fail is registered.")

n = len(RESULTS)
print(f"{sum(RESULTS)}/{n} zhao_spectrum checks passed")
sys.exit(0 if all(RESULTS) else 1)
