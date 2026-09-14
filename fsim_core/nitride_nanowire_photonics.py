"""Independent-prediction photonics envelope for the opt-in nitride-nanowire
tier (piece 3; depends on piece 1's frozen contract/evidence ledger).

Two families, kept structurally separate and never blended, per
`docs/nitride_nanowire_contract.md`:

  * ``horizontal_as_built`` -- the dispersed, end-contacted wire lying on
    SiO2/Si that Deshpande et al. actually measured.  A subwavelength
    dielectric wire supports no guided mode; collection is modeled as an
    oriented point-dipole far field above a two-interface (air/SiO2/Si)
    thin-film stack, collected by a normal-incidence objective of the given
    NA.  The reported ~70% axial linear polarization (Deshpande et al., Nat.
    Commun. 4, 1675 (2013), p.5) is attributed there to a dielectric-antenna
    effect of the thin wire, not to intrinsic dot anisotropy; this module
    treats it as an input weighting knob (``dipole_weights``), never as an
    automatic collection efficiency.

  * ``vertical_photonic`` -- a designed, tapered, top-collected wire in the
    HE11-guiding regime.  No new vector-mode eigensolver is implemented (and
    none of `fsim_core.waveguide`'s slab/ridge solver is reused for a
    cylinder): confinement/beta come from a documented, closed-form,
    dimensionless-V-number surrogate (Marcuse's Gaussian mode-field-radius
    fit) that is smooth and cutoff-free down to V=0, exactly as an isolated
    dielectric-cylinder HE11 mode must be.  Maslov and Ning, Opt. Lett. 29,
    572 (2004) and Claudon et al., Nat. Photon. 4, 174 (2010) are both
    ``evidence_status: missing`` in `verify/data/nitride_nanowire_anchors.yaml`
    (piece 1) -- no numeric confinement or extraction anchor from either
    paper is used anywhere in this module.

Both families share: gamma_X_ns = gamma_X0_ns * radiative_rate_factor (a
Purcell/rate envelope, baseline 1.0 [A]) -- collection efficiency NEVER
multiplies the decay rate, and a high beta or high first-lens number never
stands in for gamma.  A card selects exactly one family; the other family's
knobs must sit at their neutral/off default or `response` raises (Section
"family cross-talk" below), so a horizontal card can never be silently
affected by a taper/mirror/top-contact setting and vice versa.

Refractive-index provenance actually used here:
  * GaN, ordinary ray: Barker and Ilegems, Phys. Rev. B 7, 743 (1973)
    Sellmeier n^2 = 1 + 2.60 + 1.75 L^2/(L^2-0.256^2) + 4.1 L^2/(L^2-17.86^2)
    (L in micrometers), stated valid 0.35-10 um [V] (formula and citation
    read from refractiveindex.info's Barker-o transcription of the paper;
    the numeric evaluation at any lambda_nm is this module's own).
  * SiO2 (thermal oxide / fused silica): Malitson, J. Opt. Soc. Am. 55, 1205
    (1965) Sellmeier, the standard three-term fused-silica dispersion [V].
  * Si (absorbing at these visible wavelengths): a two-point complex-index
    anchor table at 450 nm (n=4.676, k=0.091) and 630 nm (n=3.879, k=0.016),
    linearly interpolated; class values consistent with Aspnes and Studna,
    Phys. Rev. B 27, 985 (1983) -- [E], not independently re-verified
    against that paper's primary table in this pass.  Extrapolation outside
    [450, 630] nm is refused (raise), matching `materials._lookup_n`'s
    refuse-to-extrapolate convention; a caller needing another wavelength
    must supply an explicit `n_substrate` override.
  * Thin-film (single-layer) interference and normal Fresnel reflection:
    Born and Wolf, *Principles of Optics*, ch. 1.6 [V, standard formula].
  * Oriented-dipole-above-an-interface s/p decomposition and image-style
    interference: Lukosz and Kunz, J. Opt. Soc. Am. 67, 1607 (1977) [DR,
    standard two-ray formalism]; the well known h=0 perfect-mirror limits
    (a dipole parallel to the mirror radiates nothing; a dipole
    perpendicular to it radiates 4x the free-space intensity) are used as
    independent checks in `verify/verify_nitride_nanowire_photonics.py`.

Everything this module omits is stated, not silently dropped: no substrate
near-field/nonradiative LDOS (so a mirror can raise the modeled "collection"
above what a real device would show at small heights -- flagged, and capped
at 1, rather than hidden); no cavity Q, DBR, or spectral linewidth filter
(owned by device integration); no dipole-orientation dependence of the
guided-mode beta.
"""
from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np
from scipy.special import jn_zeros

FAMILIES = frozenset({"horizontal_as_built", "vertical_photonic"})

# [V] Standard step-index/cylinder V-number single-mode (LP01/LP11) cutoff:
# the first zero of the J0 Bessel function, independently computed here
# (not imported from any sibling nanowire module).
V_CUTOFF_LP11 = float(jn_zeros(0, 1)[0])

# [V] Total power radiated by a unit point dipole into the full 4*pi solid
# angle, in the same normalization as the angular pattern used below
# (integral of sin^2(gamma) over 4*pi, a standard closed-form identity,
# independent of the dipole's orientation).
_FREE_SPACE_TOTAL_POWER = 8.0 * math.pi / 3.0

# [A] Marcuse (1977) fit is stated valid for 0.8 <= V <= 2.5; outside that
# window `approximation_error` grows from 0 as a simple linear distance in
# V-units.  This is a bookkeeping heuristic, not a formal error bound.
_MARCUSE_V_LO, _MARCUSE_V_HI = 0.8, 2.5

# [E] Two-point complex Si index anchor (see module docstring); linear
# interpolation between them, extrapolation refused outside [450, 630] nm.
_SI_INDEX_ANCHORS_NM = (450.0, 630.0)
_SI_INDEX_N = (4.676, 3.879)
_SI_INDEX_K = (0.091, 0.016)

_VERTICAL_ONLY_DEFAULTS = {
    "bottom_reflectivity": 0.0,
    "taper_transmission": 1.0,
    "top_contact_transmission": 1.0,
    "propagation_transmission": 1.0,
    "unguided_collection_scale": 0.0,
    "beta_scale": 1.0,
}
_HORIZONTAL_ONLY_DEFAULTS = {
    "collection_scale": 1.0,
}


class NitrideNanowirePhotonicsError(ValueError):
    """A photonics input is malformed, non-finite, or outside this model's
    declared domain (never raised for a merely unfavorable but well-formed
    physical point -- those come back with ``valid=False``)."""


def _finite(name, value):
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise NitrideNanowirePhotonicsError(f"{name} must be a finite number")
    return float(value)


def _finite_positive(name, value):
    v = _finite(name, value)
    if v <= 0.0:
        raise NitrideNanowirePhotonicsError(f"{name} must be positive")
    return v


def _fraction(name, value):
    v = _finite(name, value)
    if not 0.0 <= v <= 1.0:
        raise NitrideNanowirePhotonicsError(f"{name} must be in [0, 1]")
    return v


@dataclass(frozen=True)
class NitrideNanowirePhotonicsParams:
    """Card-level optical inputs, frozen and normalized separately per
    family (see the module docstring for provenance of each default).

    ``dipole_weights`` is ``(along_wire, transverse_inplane, vertical)``,
    fractions of an assumed EFFECTIVE dipole-orientation mixture used only
    by the horizontal far-field integral; it is a free [A] input, never
    derived from the reported polarization fraction.  ``n_wire``/``n_oxide``/
    ``n_substrate`` default to ``None``, which resolves to the built-in
    dispersion (GaN/SiO2/Si respectively) at the call's ``lambda_nm``; a
    numeric override bypasses that dispersion entirely (no wavelength
    restriction then applies to the override).
    """
    family: str
    NA: float = 0.5                      # [A] objective numerical aperture
    n_wire: float | None = None          # [V] GaN Sellmeier if None
    n_ambient: float = 1.0               # [A] air
    n_oxide: float | None = None         # [V] SiO2 Sellmeier if None
    oxide_thickness_nm: float = 100.0    # [A] design brief, Deshpande 2013 p.3; not a frozen ledger anchor
    n_substrate: complex | float | None = None  # [E] Si anchor table if None
    dipole_weights: tuple = (1.0 / 3.0, 1.0 / 3.0, 1.0 / 3.0)  # [A] isotropic default
    emitter_height_nm: float | None = None      # [A] defaults to outer_radius_nm
    collection_scale: float = 1.0        # [A] horizontal-only omitted-physics envelope
    radiative_rate_factor: float = 1.0   # [A] baseline; independent of beta/eta
    beta_scale: float = 1.0              # [A] confinement-to-beta mapping, vertical-only
    taper_transmission: float = 1.0      # [A] vertical-only, neutral = lossless
    bottom_reflectivity: float = 0.0     # [A] vertical-only, neutral = no mirror
    top_contact_transmission: float = 1.0  # [A] vertical-only, neutral = no contact loss
    propagation_transmission: float = 1.0  # [A] vertical-only, neutral = lossless guiding
    unguided_collection_scale: float = 0.0  # [A] vertical-only, neutral = excluded

    def __post_init__(self):
        if self.family not in FAMILIES:
            raise NitrideNanowirePhotonicsError(f"family must be one of {sorted(FAMILIES)}")
        _fraction("NA", self.NA)
        if self.n_wire is not None:
            _finite_positive("n_wire", self.n_wire)
        _finite_positive("n_ambient", self.n_ambient)
        if self.n_oxide is not None:
            _finite_positive("n_oxide", self.n_oxide)
        thickness = _finite("oxide_thickness_nm", self.oxide_thickness_nm)
        if thickness < 0.0:
            raise NitrideNanowirePhotonicsError("oxide_thickness_nm must be >= 0")
        if self.n_substrate is not None:
            nsub = complex(self.n_substrate)
            if (not math.isfinite(nsub.real) or not math.isfinite(nsub.imag)
                    or nsub.real <= 0.0 or nsub.imag < 0.0):
                raise NitrideNanowirePhotonicsError(
                    "n_substrate must have a positive finite real part and a "
                    "non-negative finite imaginary part")
        w = self.dipole_weights
        if not (isinstance(w, tuple) and len(w) == 3):
            raise NitrideNanowirePhotonicsError(
                "dipole_weights must be a 3-tuple (along_wire, transverse_inplane, vertical)")
        for label, value in zip(("along_wire", "transverse_inplane", "vertical"), w):
            _fraction(f"dipole_weights.{label}", value)
        if not math.isclose(sum(w), 1.0, rel_tol=0.0, abs_tol=1e-9):
            raise NitrideNanowirePhotonicsError("dipole_weights must sum to 1")
        if self.emitter_height_nm is not None:
            _finite_positive("emitter_height_nm", self.emitter_height_nm)
        _fraction("collection_scale", self.collection_scale)
        _finite_positive("radiative_rate_factor", self.radiative_rate_factor)
        _finite_positive("beta_scale", self.beta_scale)
        _fraction("taper_transmission", self.taper_transmission)
        _fraction("bottom_reflectivity", self.bottom_reflectivity)
        _fraction("top_contact_transmission", self.top_contact_transmission)
        _fraction("propagation_transmission", self.propagation_transmission)
        _fraction("unguided_collection_scale", self.unguided_collection_scale)


def _check_family_activation(params: NitrideNanowirePhotonicsParams) -> None:
    """Reject (raise) a card that sets the OTHER family's knobs away from
    their neutral default -- "insensitive or rejects explicitly" is met by
    rejecting explicitly, so a config mistake cannot silently do nothing."""
    if params.family == "horizontal_as_built":
        active = [n for n, d in _VERTICAL_ONLY_DEFAULTS.items() if getattr(params, n) != d]
        if active:
            raise NitrideNanowirePhotonicsError(
                "horizontal_as_built rejects non-default vertical-only parameters: "
                + ", ".join(active))
    else:
        active = [n for n, d in _HORIZONTAL_ONLY_DEFAULTS.items() if getattr(params, n) != d]
        if active:
            raise NitrideNanowirePhotonicsError(
                "vertical_photonic rejects non-default horizontal-only parameters: "
                + ", ".join(active))


# --------------------------------------------------------------- dispersion

def gan_ordinary_index(lambda_nm: float) -> float:
    """[V] Barker and Ilegems, Phys. Rev. B 7, 743 (1973), ordinary-ray
    Sellmeier for wurtzite GaN, valid 0.35-10 um (see module docstring)."""
    lam = _finite_positive("lambda_nm", lambda_nm)
    if not 350.0 <= lam <= 10000.0:
        raise NitrideNanowirePhotonicsError(
            "gan_ordinary_index: lambda_nm outside the cited Sellmeier's 0.35-10 um range")
    L = lam / 1000.0
    n2 = 1.0 + 2.60 + 1.75 * L * L / (L * L - 0.256 ** 2) + 4.1 * L * L / (L * L - 17.86 ** 2)
    return math.sqrt(n2)


def sio2_index(lambda_nm: float) -> float:
    """[V] Malitson, J. Opt. Soc. Am. 55, 1205 (1965) fused-silica Sellmeier."""
    lam = _finite_positive("lambda_nm", lambda_nm)
    if not 210.0 <= lam <= 3710.0:
        raise NitrideNanowirePhotonicsError(
            "sio2_index: lambda_nm outside Malitson's measured 0.21-3.71 um range")
    L = lam / 1000.0
    L2 = L * L
    n2 = (1.0 + 0.6961663 * L2 / (L2 - 0.0684043 ** 2)
          + 0.4079426 * L2 / (L2 - 0.1162414 ** 2)
          + 0.8974794 * L2 / (L2 - 9.896161 ** 2))
    return math.sqrt(n2)


def si_complex_index(lambda_nm: float) -> complex:
    """[E] Two-point linear-interpolated crystalline-Si complex index (see
    module docstring); refuses extrapolation outside [450, 630] nm."""
    lam = _finite_positive("lambda_nm", lambda_nm)
    lo, hi = _SI_INDEX_ANCHORS_NM
    if not lo - 1e-9 <= lam <= hi + 1e-9:
        raise NitrideNanowirePhotonicsError(
            f"si_complex_index: lambda_nm={lam} outside the tabulated {lo}-{hi} nm "
            "anchor range; supply an explicit n_substrate override instead")
    t = (lam - lo) / (hi - lo)
    n = _SI_INDEX_N[0] + t * (_SI_INDEX_N[1] - _SI_INDEX_N[0])
    k = _SI_INDEX_K[0] + t * (_SI_INDEX_K[1] - _SI_INDEX_K[0])
    return complex(n, k)


def v_number(outer_radius_nm: float, lambda_nm: float, n_wire: float, n_ambient: float) -> float:
    """[V] Standard step-index/cylinder dimensionless V-number (Snyder and
    Love, *Optical Waveguide Theory*, 1983): V = 2*pi*R/lambda*sqrt(n1^2-n2^2).
    Uses RADIUS (never diameter); returns 0 rather than raising when the
    wire has no index contrast (n_wire <= n_ambient), the correct V->0
    no-guiding limit."""
    r = _finite_positive("outer_radius_nm", outer_radius_nm)
    lam = _finite_positive("lambda_nm", lambda_nm)
    contrast = max(0.0, n_wire * n_wire - n_ambient * n_ambient)
    return 2.0 * math.pi * r / lam * math.sqrt(contrast)


# ------------------------------------------------------- vertical surrogate

def _marcuse_w_over_a(V: float) -> float:
    """[V] Marcuse, Bell Syst. Tech. J. 56, 703 (1977) Gaussian mode-field
    radius w/a (coefficients as reproduced in standard optical-waveguide
    texts, e.g. Okamoto, *Fundamentals of Optical Waveguides*, Eq. 3.30),
    stated accurate for 0.8<=V<=2.5.  The SAME closed form stays smooth,
    positive, and monotonically decreasing in V for every V>0 (checked in
    the verifier), so it is used unmodified outside that window too, at the
    cost of unquantified (but bookkept) accuracy -- see `approximation_error`.
    """
    if V <= 0.0:
        return math.inf
    return 0.65 + 1.619 * V ** -1.5 + 2.879 * V ** -6.0


def _gaussian_core_fraction(w_over_a: float) -> float:
    """[DR] Fraction of a circular Gaussian beam's power (intensity profile
    exp(-2 r^2/w^2)) falling inside radius a: 1 - exp(-2 a^2/w^2), a direct
    integral of the same Gaussian intensity profile whose spot size
    `_marcuse_w_over_a` supplies. w_over_a=inf (V=0) gives exactly 0."""
    if not math.isfinite(w_over_a):
        return 0.0
    return 1.0 - math.exp(-2.0 / (w_over_a * w_over_a))


def _extrapolation_metric(V: float) -> float:
    return float(max(0.0, _MARCUSE_V_LO - V) + max(0.0, V - _MARCUSE_V_HI))


def _diffraction_half_angle_rad(mode_radius_nm: float, lambda_nm: float) -> float:
    """[A] theta_div ~ lambda/(pi*w): the guided mode's own size (taken as
    the wire's outer radius) sets the diffraction-limited divergence of the
    light emerging from an idealized taper, in the same small-angle Gaussian
    form used for the objective-acceptance integral below."""
    return lambda_nm / (math.pi * mode_radius_nm)


def _objective_acceptance(theta_max_rad: float, theta_div_rad: float) -> float:
    """[DR] Fraction of a circular Gaussian beam of divergence half-angle
    theta_div falling within a cone of half-angle theta_max: the same
    1-exp(-2 x^2) integral as `_gaussian_core_fraction`, in angle space."""
    if theta_div_rad <= 0.0:
        return 0.0
    return 1.0 - math.exp(-2.0 * theta_max_rad * theta_max_rad / (theta_div_rad * theta_div_rad))


def _vertical_collection(params: NitrideNanowirePhotonicsParams, lambda_nm: float,
                          outer_radius_nm: float, V: float):
    w_over_a = _marcuse_w_over_a(V)
    confinement = _gaussian_core_fraction(w_over_a)
    beta_raw = confinement * params.beta_scale
    beta_clipped = beta_raw > 1.0 or beta_raw < 0.0
    beta = min(max(beta_raw, 0.0), 1.0)

    # [DR] beta is BOTH propagation directions (contract); split evenly by
    # the disc's up/down mirror symmetry within the wire.
    beta_up = beta / 2.0
    beta_down = beta / 2.0
    top_up = beta_up * params.propagation_transmission
    # [A] disc-near-base geometry: the down-then-reflect-then-up path
    # crosses the full guided length once (same propagation_transmission),
    # the short disc-to-bottom hop treated as lossless.
    top_down = beta_down * params.bottom_reflectivity * params.propagation_transmission
    guided_at_top = top_up + top_down
    through_taper = guided_at_top * params.taper_transmission * params.top_contact_transmission

    theta_max = math.asin(min(params.NA, 1.0))
    theta_div = _diffraction_half_angle_rad(outer_radius_nm, lambda_nm)
    obj_accept = _objective_acceptance(theta_max, theta_div)
    guided_collected = through_taper * obj_accept

    unguided_fraction = max(0.0, 1.0 - beta)
    unguided_collected = unguided_fraction * params.unguided_collection_scale

    eta_raw = guided_collected + unguided_collected
    diag = {
        "confinement_fraction": confinement, "beta_up": beta_up, "beta_down": beta_down,
        "top_up": top_up, "top_down": top_down, "guided_at_top": guided_at_top,
        "through_taper": through_taper, "theta_max_rad": theta_max,
        "theta_div_rad": theta_div, "obj_accept": obj_accept,
        "unguided_fraction": unguided_fraction, "unguided_collected": unguided_collected,
        "beta_clipped": beta_clipped, "additional_modes_possible": V > V_CUTOFF_LP11,
        "V_cutoff_LP11": V_CUTOFF_LP11,
    }
    return beta, eta_raw, diag


# ---------------------------------------------------- horizontal dipole-cone

_ORIENTATIONS = (
    ("along_wire", (1.0, 0.0, 0.0)),
    ("transverse_inplane", (0.0, 1.0, 0.0)),
    ("vertical", (0.0, 0.0, 1.0)),
)


def _snell_cos_t(n_i, n_t, cos_i):
    cos_i = np.asarray(cos_i, dtype=complex)
    sin_i2 = 1.0 - cos_i * cos_i
    sin_t2 = (n_i / n_t) ** 2 * sin_i2
    return np.sqrt(1.0 - sin_t2)


def _fresnel_rs(n_i, n_t, cos_i, cos_t):
    return (n_i * cos_i - n_t * cos_t) / (n_i * cos_i + n_t * cos_t)


def _fresnel_rp(n_i, n_t, cos_i, cos_t):
    return (n_t * cos_i - n_i * cos_t) / (n_t * cos_i + n_i * cos_t)


def stack_reflection(n0, n1, n2, thickness_nm, lambda_nm, cos0, pol):
    """[V] Two-interface (air/layer/substrate) thin-film power-reflection
    amplitude, Born and Wolf, *Principles of Optics*, ch. 1.6.  ``cos0`` may
    be a scalar or a numpy array; ``n2`` may be complex (absorbing)."""
    cos0c = np.asarray(cos0, dtype=complex)
    cos1 = _snell_cos_t(n0, n1, cos0c)
    cos2 = _snell_cos_t(n1, n2, cos1)
    if pol == "s":
        r01 = _fresnel_rs(n0, n1, cos0c, cos1)
        r12 = _fresnel_rs(n1, n2, cos1, cos2)
    elif pol == "p":
        r01 = _fresnel_rp(n0, n1, cos0c, cos1)
        r12 = _fresnel_rp(n1, n2, cos1, cos2)
    else:
        raise NitrideNanowirePhotonicsError("pol must be 's' or 'p'")
    beta = 2.0 * math.pi / lambda_nm * n1 * thickness_nm * cos1
    phase = np.exp(2j * beta)
    return (r01 + r12 * phase) / (1.0 + r01 * r12 * phase)


def _hemisphere_nodes(theta_max, n_theta, n_phi):
    nodes, weights = np.polynomial.legendre.leggauss(n_theta)
    u_lo, u_hi = math.cos(theta_max), 1.0
    u = 0.5 * (u_hi - u_lo) * nodes + 0.5 * (u_hi + u_lo)
    w_u = 0.5 * (u_hi - u_lo) * weights
    phi = np.linspace(0.0, 2.0 * math.pi, n_phi, endpoint=False)
    w_phi = np.full(n_phi, 2.0 * math.pi / n_phi)
    return u, w_u, phi, w_phi


def dipole_collection_fraction(p, NA, n_ambient, n_oxide, n_substrate,
                                oxide_thickness_nm, height_nm, lambda_nm,
                                n_theta=48, n_phi=96):
    """[DR] Coherent direct+substrate-reflected far-field power of a point
    dipole oriented along unit vector ``p``, collected within a normal-
    incidence objective of the given NA, normalized by the FIXED free-space
    4*pi total (`_FREE_SPACE_TOTAL_POWER`) -- see the module docstring for
    why this can legitimately exceed 1 near a strongly reflecting interface
    at small height (an omitted LDOS/rate effect, not a bug)."""
    if NA <= 0.0:
        return 0.0
    theta_max = math.asin(min(NA, 1.0))
    u, w_u, phi, w_phi = _hemisphere_nodes(theta_max, n_theta, n_phi)
    U = u[:, None]
    sin_t = np.sqrt(np.clip(1.0 - U * U, 0.0, None))
    cos_p = np.cos(phi)[None, :]
    sin_p = np.sin(phi)[None, :]

    s_x, s_y = -np.sin(phi), np.cos(phi)
    Es = (p[0] * s_x + p[1] * s_y)[None, :]

    tm_up_z = -sin_t
    Ep_up = p[0] * (U * cos_p) + p[1] * (U * sin_p) + p[2] * tm_up_z
    Ep_down = p[0] * (-U * cos_p) + p[1] * (-U * sin_p) + p[2] * tm_up_z

    r_s = stack_reflection(n_ambient, n_oxide, n_substrate, oxide_thickness_nm,
                            lambda_nm, u, "s")[:, None]
    r_p = stack_reflection(n_ambient, n_oxide, n_substrate, oxide_thickness_nm,
                            lambda_nm, u, "p")[:, None]

    k0 = 2.0 * math.pi / lambda_nm
    phase = np.exp(2j * k0 * n_ambient * height_nm * U)

    Es_total = Es * (1.0 + r_s * phase)
    Ep_total = Ep_up + r_p * phase * Ep_down

    intensity = np.abs(Es_total) ** 2 + np.abs(Ep_total) ** 2
    collected = np.sum(intensity * w_u[:, None] * w_phi[None, :])
    return float(collected / _FREE_SPACE_TOTAL_POWER)


def _horizontal_collection(params: NitrideNanowirePhotonicsParams, lambda_nm: float,
                            outer_radius_nm: float):
    n_oxide = params.n_oxide if params.n_oxide is not None else sio2_index(lambda_nm)
    n_sub = params.n_substrate if params.n_substrate is not None else si_complex_index(lambda_nm)
    height = params.emitter_height_nm if params.emitter_height_nm is not None else outer_radius_nm
    per_orientation = {}
    total = 0.0
    for (label, vec), weight in zip(_ORIENTATIONS, params.dipole_weights):
        eta = dipole_collection_fraction(vec, params.NA, params.n_ambient, n_oxide, n_sub,
                                          params.oxide_thickness_nm, height, lambda_nm)
        per_orientation[label] = eta
        total += weight * eta
    diag = {
        "n_oxide_used": n_oxide, "n_substrate_used": n_sub,
        "emitter_height_nm_used": height, "eta_by_orientation": per_orientation,
    }
    return total, diag


# --------------------------------------------------------------- top level

def response(params: NitrideNanowirePhotonicsParams, *, lambda_nm: float,
             outer_radius_nm: float, gamma_X0_ns: float, gamma_XX0_ns: float) -> dict:
    """One family's optical response at one wavelength and outer radius.

    X and XX are evaluated at the SAME ``lambda_nm`` (the signature admits
    only one wavelength): ``eta_collection_X`` and ``eta_collection_XX`` are
    therefore numerically identical by construction here -- an explicit [A]
    simplification (the X/XX wavelength split, ~meV scale, is neglected in
    this geometric-optics collection model), not a coincidence to rely on
    once X/XX detuning is modeled elsewhere.
    """
    if not isinstance(params, NitrideNanowirePhotonicsParams):
        raise TypeError("params must be a NitrideNanowirePhotonicsParams")
    lam = _finite_positive("lambda_nm", lambda_nm)
    radius = _finite_positive("outer_radius_nm", outer_radius_nm)
    g_x0 = _finite_positive("gamma_X0_ns", gamma_X0_ns)
    g_xx0 = _finite_positive("gamma_XX0_ns", gamma_XX0_ns)
    _check_family_activation(params)

    n_wire = params.n_wire if params.n_wire is not None else gan_ordinary_index(lam)
    V = v_number(radius, lam, n_wire, params.n_ambient)
    radius_over_lambda = radius / lam

    gamma_X_ns = g_x0 * params.radiative_rate_factor
    gamma_XX_ns = g_xx0 * params.radiative_rate_factor

    notes: list[str] = []
    invalid_reasons: list[str] = []
    approximation_error = 0.0
    diagnostics: dict = {"n_wire_used": n_wire}

    if params.family == "horizontal_as_built":
        beta_HE11 = "not_applicable"
        eta_raw, diag = _horizontal_collection(params, lam, radius)
        eta_x_raw = eta_raw * params.collection_scale
        diagnostics.update(diag)
        diagnostics["eta_geom_raw"] = eta_raw
        notes.append(
            "no substrate near-field/nonradiative LDOS is modeled; this is a "
            "collection envelope, and collection_scale in [0, ~0.3-1.0] is the "
            "assumed range for that omitted physics")
    else:
        if n_wire <= params.n_ambient:
            invalid_reasons.append(
                "n_wire must exceed n_ambient for a designed vertical_photonic wire to guide")
        beta_HE11, eta_x_raw, diag = _vertical_collection(params, lam, radius, V)
        diagnostics.update(diag)
        approximation_error = _extrapolation_metric(V)
        if diag["beta_clipped"]:
            notes.append("beta_HE11 clipped to [0, 1] after applying beta_scale")
        if diag["additional_modes_possible"]:
            notes.append(
                f"V_number={V:.4f} exceeds the single-mode (LP11, V={V_CUTOFF_LP11:.6f}) "
                "cutoff; higher-order guided content is not modeled by this "
                "single-mode surrogate")
        if approximation_error > 0.0:
            notes.append(
                f"V_number={V:.4f} is outside the Marcuse (1977) calibration "
                f"window [{_MARCUSE_V_LO}, {_MARCUSE_V_HI}]; beta_HE11 is an [A] "
                "smooth extrapolation of that fit, not independently validated there")

    if eta_x_raw > 1.0 + 1e-9:
        notes.append(
            "raw collected probability exceeded 1 before capping (expected near a "
            "strongly reflecting interface at small height, since this factorized "
            "model does not move that enhancement into gamma); capped at 1.0")
    eta_x = float(min(max(eta_x_raw, 0.0), 1.0))
    eta_xx = eta_x

    valid = len(invalid_reasons) == 0

    provenance = {
        "V_number": "[V] Snyder and Love, Optical Waveguide Theory (1983), standard "
                    "step-index V-number definition; RADIUS (not diameter), see v_number()",
        "n_wire": "[V] Barker and Ilegems, Phys. Rev. B 7, 743 (1973), GaN ordinary "
                  "Sellmeier" if params.n_wire is None else "[A] caller-supplied override",
        "radiative_rate_factor": "[A] baseline 1.0, independent Purcell/rate envelope; "
                                  "never multiplied by beta or eta_collection",
        "maslov_claudon_status": "Maslov and Ning (2004) and Claudon et al. (2010) are "
                                  "evidence_status=missing in verify/data/"
                                  "nitride_nanowire_anchors.yaml; no numeric anchor from "
                                  "either paper is used in this module",
    }
    if params.family == "horizontal_as_built":
        provenance["n_oxide"] = ("[V] Malitson, J. Opt. Soc. Am. 55, 1205 (1965)"
                                  if params.n_oxide is None else "[A] caller-supplied override")
        provenance["n_substrate"] = ("[E] two-point Si n,k anchor table, class values "
                                      "consistent with Aspnes and Studna (1983); not "
                                      "independently re-verified against that paper's table"
                                      if params.n_substrate is None else "[A] caller-supplied override")
        provenance["dipole_model"] = ("[DR] point-dipole s/p far field above a thin-film "
                                       "stack; Born and Wolf ch. 1.6 (Fresnel/thin-film); "
                                       "Lukosz and Kunz, JOSA 67, 1607 (1977) (two-ray "
                                       "interference formalism)")
        provenance["dipole_weights"] = ("[A] assumed effective dipole-orientation mixture; "
                                         "NOT derived from the reported ~70% axial "
                                         "polarization (Deshpande et al. 2013, p.5)")
    else:
        provenance["confinement_surrogate"] = (
            "[V]/[A] Marcuse (1977) Gaussian mode-field-radius fit, calibrated "
            f"{_MARCUSE_V_LO}<=V<={_MARCUSE_V_HI}; smoothly extrapolated outside that "
            "window (approximation_error quantifies the extrapolation distance)")
        provenance["confinement_to_beta_mapping"] = (
            "[A] confinement fraction is mapped onto beta_HE11 (times beta_scale); "
            "this mapping itself is an assumption, not an independently derived beta")

    return {
        "family": params.family,
        "lambda_nm": lam,
        "radius_over_lambda": radius_over_lambda,
        "V_number": V,
        "beta_HE11": beta_HE11,
        "eta_collection_X": eta_x,
        "eta_collection_XX": eta_xx,
        "radiative_rate_factor": params.radiative_rate_factor,
        "gamma_X_ns": gamma_X_ns,
        "gamma_XX_ns": gamma_XX_ns,
        "valid": valid,
        "invalid_reasons": invalid_reasons,
        "approximation_error": approximation_error,
        "provenance": provenance,
        "notes": notes,
        "diagnostics": diagnostics,
    }
