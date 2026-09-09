"""Additive cavity response for visible nitride-dot design studies.

This is a factorized rate-times-filter model, not a self-consistent LDOS or
radiation-pattern calculation.  Purcell feeding and the collected spectrum
can therefore be correlated in a real device.  ``eta_out`` is the bounded
first-lens collection after the one modeled spectral acceptance.

The default Q, normalized mode volume, placement overlap, collection, and
cavity temperature slope are independent exploratory design inputs [A], not
a jointly measured cavity.  Q=167 is a reported example (Taylor et al.,
Nanoscale Res. Lett. 5, 608, 2010, DOI 10.1007/s11671-009-9514-4) [V]; its
lack of observed Purcell enhancement must not be read as a limit on all Q.
"""
from __future__ import annotations

from dataclasses import dataclass
from math import isfinite, pi

from .dbr import dbr_reflectivity, power_RT, quarter_wave_stack
from .spectral import cavity_transmission, epsilon2


class NitrideCavityParameterError(ValueError):
    """A cavity input is non-finite or outside this model's physical domain."""


class NitrideCavityEnergyError(ValueError):
    """The tracked cavity energy is not positive."""


@dataclass(frozen=True)
class NitrideCavityParams:
    """Explicit planar/micropillar design inputs.

    ``mode_volume_norm`` is V/(lambda/n)^3.  The ideal single-mode expression
    below is from Purcell, Phys. Rev. 69, 681 (1946) [DR].  Geometry, spatial
    overlap, eta_out, and dEdT_cav_meV_K are assumptions [A].
    """
    Q: float = 2000.0
    mode_volume_norm: float = 2.0
    spatial_overlap: float = 1.0
    eta_out: float = 0.1
    T_track: float = 300.0
    dEdT_cav_meV_K: float = -0.04
    detuning_offset_meV: float = 0.0
    purcell_enabled: bool = True

    def __post_init__(self):
        for name in ("Q", "mode_volume_norm", "spatial_overlap", "eta_out",
                     "T_track", "dEdT_cav_meV_K", "detuning_offset_meV"):
            value = getattr(self, name)
            if not isinstance(value, (int, float)) or not isfinite(value):
                raise NitrideCavityParameterError(f"{name} must be finite")
        if self.Q <= 0.0:
            raise NitrideCavityParameterError("Q must be positive")
        if self.mode_volume_norm <= 0.0:
            raise NitrideCavityParameterError("mode_volume_norm must be positive")
        for name in ("spatial_overlap", "eta_out"):
            if not 0.0 <= getattr(self, name) <= 1.0:
                raise NitrideCavityParameterError(f"{name} must be in [0, 1]")


def _finite_positive(name, value):
    if not isinstance(value, (int, float)) or not isfinite(value) or value <= 0.0:
        raise NitrideCavityParameterError(f"{name} must be finite and positive")
    return float(value)


def varshni_shift_eV(T_K, alpha_eV_K=8.3e-4, beta_K=825.0):
    """E(T)-E(0) in eV for an optional nitride tracking diagnostic.

    Deshpande et al., Nat. Commun. 4, 1675 (2013), Fig. 2b report
    alpha=(8.3 +/- 0.4)e-4 eV/K and beta=825 +/- 16.7 K [V].  It is not
    imposed on composition-specific confinement energies supplied to response.
    """
    T = float(T_K)
    alpha = _finite_positive("alpha_eV_K", alpha_eV_K)
    beta = _finite_positive("beta_K", beta_K)
    if not isfinite(T) or T < 0.0:
        raise NitrideCavityParameterError("T_K must be finite and non-negative")
    return -alpha * T * T / (T + beta)


def dbr_diagnostic(lambda_nm, *, n_high=2.1, n_low=1.46, pairs=10, n_in=1.0,
                   n_out=1.0):
    """Ideal lossless quarter-wave DBR diagnostic, not manufacturing evidence.

    The SiO2/Ta2O5-like numeric indices and ten pairs are exploratory [E/A].
    The transfer-matrix calculation follows Born and Wolf, *Principles of
    Optics*, ch. 1.6 [DR]; it excludes absorption, sidewalls, contacts,
    cracks, electrical access, collection, and 3-D mode volume.

    ``n_in``/``n_out`` default to free space on both sides; they are NOT a
    GaN-substrate assumption.  A real substrate (n_out ~ 2.4 for GaN) must be
    supplied explicitly by the caller -- the module carries no default for
    it, and the free-space default understates R relative to a GaN-backed
    stack (e.g. R=0.99722 free space vs R=0.99884 with n_out=2.4 at the same
    10-pair 2.1/1.46 design [E], since the higher-index exit medium raises
    the admittance contrast the stack presents).  ``transmittance`` is the
    transfer-matrix power T from ``dbr.power_RT`` on this same stack, not a
    lossless-medium ``1-R`` shortcut.
    """
    lam = _finite_positive("lambda_nm", lambda_nm)
    nh, nl = _finite_positive("n_high", n_high), _finite_positive("n_low", n_low)
    n_in = _finite_positive("n_in", n_in)
    n_out = _finite_positive("n_out", n_out)
    if nh <= nl:
        raise NitrideCavityParameterError("n_high must exceed n_low")
    if not isinstance(pairs, int) or pairs <= 0:
        raise NitrideCavityParameterError("pairs must be a positive integer")
    stack = dbr_reflectivity(nh, nl, pairs, lam, n_in=n_in, n_out=n_out)
    R = float(stack["R_at_lambda0"])
    layers = quarter_wave_stack(nh, nl, pairs, lam)
    _, T = power_RT(layers, lam, n_in, n_out)
    return {
        "reflectance": R,
        "transmittance": float(T),
        "d_high_nm": float(layers[0].d_nm),
        "d_low_nm": float(layers[1].d_nm),
        "stopband_nm": tuple(float(x) for x in stack["stopband_analytic_nm"]),
        "provenance": "[DR] Born and Wolf, Principles of Optics, ch. 1.6; indices/pairs [E/A]",
    }


def _line_overlap(line_detuning_meV, gamma_meV, kappa_meV):
    """Broadened single-line spectral overlap O_i, delegating to
    ``spectral.cavity_transmission`` (the peak-normalized Lorentzian-cavity
    kernel) instead of retyping its algebra: at line_detuning=0 this is
    kappa/(kappa+gamma), matching the on-resonance overlap in
    ``cavity.purcell_eff`` [DR within the analytic model]."""
    return cavity_transmission(line_detuning_meV, gamma_meV, kappa_meV)


def response(params, *, T_K, E_X_eV, E_X_track_eV, gamma_X_meV,
             gamma_XX_meV, delta_xx_meV, gamma_X0_ns, gamma_XX0_ns,
             w_meV=None, dx_w_meV=0.0):
    """Return the additive, broadened cavity response at one temperature.

    Positive detuning is emitter above cavity.  E_XX=E_X-delta_xx, so an
    antibinding (negative) delta_xx retains its physical sign.  Line rates
    are enhanced once; t_X/t_XX contain exactly one slit/cavity acceptance.

    ``purcell_enabled=False`` zeroes only the ADDED Purcell channel
    (Fp_add=0, F_eff_*=1); the cavity's own Lorentzian spectral acceptance
    still narrows t_X/t_XX via epsilon2's kappa argument.  This flag is
    "cavity without added Purcell", not a no-cavity baseline -- collection
    stays cavity-filtered either way.

    Every out-of-domain input (non-finite, non-positive energy/Q/volume,
    an out-of-range efficiency, or a non-positive tracked cavity energy)
    raises a named error above rather than being reported through
    ``valid``/``invalid_reasons``; those two fields are therefore always
    ``True``/``[]`` in a dict this function returns -- a construction that
    only completes once every input has already been validated.
    """
    if not isinstance(params, NitrideCavityParams):
        raise NitrideCavityParameterError("params must be NitrideCavityParams")
    T = _finite_positive("T_K", T_K)
    ex = _finite_positive("E_X_eV", E_X_eV)
    ex_track = _finite_positive("E_X_track_eV", E_X_track_eV)
    gx = _finite_positive("gamma_X_meV", gamma_X_meV)
    gxx = _finite_positive("gamma_XX_meV", gamma_XX_meV)
    dx = float(delta_xx_meV)
    dxw = float(dx_w_meV)
    if not isfinite(dx) or not isfinite(dxw):
        raise NitrideCavityParameterError("detunings must be finite")
    rate_x = _finite_positive("gamma_X0_ns", gamma_X0_ns)
    rate_xx = _finite_positive("gamma_XX0_ns", gamma_XX0_ns)
    if w_meV is not None:
        w_meV = _finite_positive("w_meV", w_meV)

    # Anchored once per geometry, not re-tuned at each temperature. [A slope]
    ecav = (ex_track + params.dEdT_cav_meV_K * (T - params.T_track) / 1000.0
            - params.detuning_offset_meV / 1000.0)
    if not isfinite(ecav) or ecav <= 0.0:
        raise NitrideCavityEnergyError("tracked cavity energy must be positive")
    kappa = 1000.0 * ecav / params.Q
    detuning = 1000.0 * (ex - ecav)
    detuning_xx = detuning - dx
    fp = ((3.0 / (4.0 * pi * pi)) * params.Q / params.mode_volume_norm
          * params.spatial_overlap) if params.purcell_enabled else 0.0

    fx = 1.0 + fp * _line_overlap(detuning, gx, kappa)
    fxx = 1.0 + fp * _line_overlap(detuning_xx, gxx, kappa)
    accepted = epsilon2(dx, gx, gxx, w=w_meV, kappa=kappa,
                        dx_w=dxw, dx_c=detuning)
    return {
        "E_cav_eV": ecav, "kappa_meV": kappa, "detuning_meV": detuning,
        "Fp_add": fp, "F_eff_X": fx, "F_eff_XX": fxx,
        "gamma_X_ns": rate_x * fx, "gamma_XX_ns": rate_xx * fxx,
        "t_X": accepted.t_x, "t_XX": accepted.t_xx,
        "eta_out": params.eta_out, "valid": True, "invalid_reasons": [],
        "provenance": {
            "purcell": "[DR] Purcell, Phys. Rev. 69, 681 (1946); V and overlap [A]",
            "overlap": "[DR] Lorentzian convolution extension of cavity.purcell_eff",
            "tracking": "[A] cavity slope/design tracking; emitter energy caller supplied",
            "collection": "[A] factorized rate-times-filter, eta_out after acceptance",
        },
    }
