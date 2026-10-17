"""Pulse-shaper (SLM) pulse pair of Karli et al. 2022, Sec. V B (arXiv:2203.00712v1), as a time-domain Rabi field.

The paper builds the SUPER pulses in the frequency domain: ONE broadband Gaussian spectrum is multiplied by a
two-Gaussian FIELD mask and inverse-Fourier transformed; the complex Omega(t) drives the two-level Hamiltonian
H = hbar w0 |x><x| - (hbar/2)(Omega* |g><x| + Omega |x><g|) (their Eq. 1). Quantities (all [V] Sec. V B unless noted):

  broadband   Omega~(w) = alpha~ / (sqrt(2 pi) sigma~) exp(-(w - w_C)^2 / (2 sigma~^2))     (Eq. 2)
              spectral FWHM 2 sqrt(2 ln 2) hbar sigma~ = 17 meV, hbar (w_C - w0) = -8.4 meV, alpha~ = 470 pi
  mask        f(w) = T~1 exp(-(w - w1)^2 / (2 sigma_m^2)) + T~2 exp(-(w - w2)^2 / (2 sigma_m^2))   (Eq. 3)
              hbar sigma_m = 0.2 meV; T~_i is the FIELD transmission, T_i = T~_i^2 is the intensity transmission (the
              y-axis of their Fig. 2); a transmission of 5 % between the two peaks (non-zero SLM extinction)
  detunings   Delta_i = hbar (w_i - w0), paper convention (red = negative)

Conventions fixed here (Fable, r2-karli-pulse-answer.md) [DR]:
  * alpha~ is the SPECTRAL integral over ENERGY in meV (the Reiter/Axt unit system): Int Omega~ dE = alpha~ = 470 pi, and
    (1 / 2 pi) Int Omega~ f exp(-i E t / hbar) dE is hbar Omega(t) in meV; Omega in rad/ps is that over hbar. The envelope
    area of mask line i is T~_i Omega~(E_i) (16.33 pi for line 1 at T1 = 0.5, Delta1 = -4.9 meV; 24.8 pi sqrt(T2) at -10.6).
    The rad/ps reading is short by 1/hbar = 1.519 and cannot reach the paper's 97 %.
  * the 5 % floor is a FIELD weight of 0.05 added to the mask on the interval between the two peak centres [Delta2, Delta1]
    (Eq. 3 plus "a transmission of 5 % between the two peaks"); the alternative (intensity 5 %, field 0.2236) is not used.
  * the spectrum is sampled on a uniform energy grid of step dE (alias period h / dE must exceed the pulse window);
    boundary cells of the floor interval carry their fractional overlap, so the field is continuous in Delta2.
  * the frame is the one of pulse 1 (the engine's convention): carriers exp(-i (E - Delta1) t / hbar).

The field is the sum of three linear parts, Omega = sqrt(T1) F1 + sqrt(T2) F2 + F0 (F0 = the floor), so T1 and T2 scale
amplitudes without recomputing the transforms.
"""
from __future__ import annotations

import math

import numpy as np

from .. import lindblad as lb

HBAR_MEV_PS = lb.HBAR_MEV_NS * 1e3      # meV ps
_CHUNK = 4000


class SlmField:
    def __init__(self, d1_meV, d2_meV, T1, T2, alpha_rad, bb_fwhm_meV, bb_centre_meV, mask_sigma_meV, floor, dE_meV=0.0125):
        self.d1, self.d2 = float(d1_meV), float(d2_meV)
        self.T1, self.T2 = float(T1), float(T2)
        self.alpha, self.floor = float(alpha_rad), float(floor)
        self.sig_m = float(mask_sigma_meV)
        self.dE = float(dE_meV)
        sig_E = float(bb_fwhm_meV) / (2.0 * math.sqrt(2.0 * math.log(2.0)))
        self._sig_E, self._E_c = sig_E, float(bb_centre_meV)
        self._parts = None

    # ---- spectra
    def _bb(self, E):
        return self.alpha / (math.sqrt(2.0 * math.pi) * self._sig_E) * np.exp(-0.5 * ((E - self._E_c) / self._sig_E) ** 2)

    def _grid(self, lo, hi):
        n = max(2, int(math.ceil((hi - lo) / self.dE)))
        return lo + (np.arange(n) + 0.5) * ((hi - lo) / n), (hi - lo) / n

    def _parts_build(self):
        if self._parts is not None:
            return self._parts
        parts = []
        for d in (self.d1, self.d2):                                              # Gaussian mask lines
            E, dE = self._grid(d - 7.0 * self.sig_m, d + 7.0 * self.sig_m)
            w = self._bb(E) * np.exp(-0.5 * ((E - d) / self.sig_m) ** 2)
            parts.append((E, w * (dE / HBAR_MEV_PS) / (2.0 * math.pi)))
        lo, hi = min(self.d1, self.d2), max(self.d1, self.d2)                     # floor interval [Delta2, Delta1]
        E, dE = self._grid(lo - self.dE, hi + self.dE)
        edge_lo, edge_hi = E - 0.5 * dE, E + 0.5 * dE
        frac = np.clip(np.minimum(edge_hi, hi) - np.maximum(edge_lo, lo), 0.0, None) / dE     # fractional overlap with the interval
        parts.append((E, self.floor * frac * self._bb(E) * (dE / HBAR_MEV_PS) / (2.0 * math.pi)))
        self._parts = parts
        return parts

    def series(self, x_ns):
        """(F1, F2, F0) in rad/ns at times x_ns (centred on the pulse maximum), frame of pulse 1."""
        x = np.asarray(x_ns, dtype=float) * 1e3                                   # ps
        out = []
        for E, w in self._parts_build():
            ph = (E - self.d1) / HBAR_MEV_PS                                      # rad/ps
            acc = np.zeros(x.shape, dtype=complex)
            for i in range(0, x.size, _CHUNK):
                acc[i:i + _CHUNK] = np.exp(-1j * np.outer(x[i:i + _CHUNK], ph)) @ w
            out.append(acc * 1e3)                                                 # rad/ps -> rad/ns
        return tuple(out)

    def omega(self, x_ns):
        f1, f2, f0 = self.series(x_ns)
        return math.sqrt(self.T1) * f1 + math.sqrt(self.T2) * f2 + f0

    def peak_rad_ns(self):
        """|Omega(0)|: every part is in phase at t = 0, so this is the maximum (upper bound otherwise) [DR]."""
        s = [float(np.sum(w)) * 1e3 for _, w in self._parts_build()]
        return math.sqrt(self.T1) * s[0] + math.sqrt(self.T2) * s[1] + s[2]

    def line_area_rad(self, which):
        """Envelope area (rad) of mask line 1 or 2 at unit field transmission: the broadband FT at the line centre."""
        d = self.d1 if which == 1 else self.d2
        return float(self._bb(np.array([d]))[0])

    @staticmethod
    def mask_fwhm_ps(sig_m_meV):
        """Intensity FWHM (ps) of the transform-limited Gaussian of a field mask of sigma sig_m: 2 sqrt(ln 2) hbar / sigma_m."""
        return 2.0 * math.sqrt(math.log(2.0)) * HBAR_MEV_PS / float(sig_m_meV)
