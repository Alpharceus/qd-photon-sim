"""Scalar effective-index model for epitaxial ridge-waveguide edge emission.

The vertical solver is the scalar TE/TM slab eigenproblem (the TM result uses
the same weak-guidance scalar approximation).  This is appropriate for the
first design pass, but is not a replacement for a full-vector FEM calculation.
The guided-mode Purcell factor follows Lecamp, Lalanne & Hugonin, PRL 99,
023902 (2007) [DR]; the high-beta comparison point is Arcari *et al.*, PRL
113, 093603 (2014) [V].  Index values used by :func:`hkust_ridge_stack` come
from ``materials.MATERIAL_EXTRA`` (Schubert *et al.*, JAP 77, 3416 (1995))
[DR].  Ridge loss 5 /cm and Gaussian far-field collection are class/design
estimates [E].
"""
from __future__ import annotations

from dataclasses import dataclass, field
from math import exp, pi
import numpy as np
from scipy.linalg import eigh_tridiagonal
from scipy.special import erf

from .materials import refractive_index


@dataclass(frozen=True)
class Layer:
    """A vertical layer.  First and last layers also define the exterior [A].

    `material` (council review 2026-09-05 item 5) is an OPTIONAL materials.py
    label ("InP", "Al0.52In0.48P", ...) naming the source of `n`.  It is never
    consulted for the mode solve itself (that always uses the numeric `n`
    supplied); it exists so edge_emission's group-index finite difference can
    re-resolve a wavelength-aware layer's index at lambda +/- 5 nm via
    materials.refractive_index.  None (default; every pre-existing call site)
    means "numeric-only, non-dispersive" -- legacy behaviour is unaffected."""
    name: str
    n: float
    thickness_nm: float
    is_dot: bool = False
    material: str | None = None


@dataclass
class SlabMode:
    n_eff: float
    z_nm: np.ndarray
    field: np.ndarray
    layers: list[Layer]
    confinement: dict[str, float]
    mode_width_um: float

    def gamma_layer(self, name: str) -> float:
        return self.confinement[name]


@dataclass
class RidgeMode:
    n_eff: float
    vertical: SlabMode
    lateral: SlabMode
    n_ridge: float
    n_outside: float
    A_mode_um2: float
    # Gaussian-equivalent 1/e^2 intensity radii (not 1/e radii), in um.
    wx_um: float
    wy_um: float


@dataclass
class EdgeResult:
    n_eff: float
    n_g: float
    Gamma_dot: float
    A_mode_um2: float
    F_wg: float
    beta: float
    T_facet: float
    eta_prop: float
    eta_NA: float
    eta_total: float
    notes: list[str] = field(default_factory=list)


def _grid(layers: list[Layer], max_step_nm: float = 2.0):
    if len(layers) < 3:
        raise ValueError("a slab needs at least lower cladding, core, and upper cladding")
    if any(x.thickness_nm <= 0 or x.n <= 0 for x in layers):
        raise ValueError("layer thicknesses and refractive indices must be positive")
    # The supplied claddings are physical finite layers; Dirichlet ends are
    # safe once their evanescent tails have decayed.  A 2-nm grid makes the
    # independently checked symmetric slab dispersion accurate to <1e-4 [DR].
    edges = np.r_[0.0, np.cumsum([x.thickness_nm for x in layers])]
    total = edges[-1]
    # Resolve ultrathin active sheets explicitly; otherwise 2 nm is ample.
    step = min(max_step_nm, min(x.thickness_nm for x in layers) / 2.0)
    count = max(301, int(np.ceil(total / step)) + 1)
    z = np.linspace(0.0, total, count)
    n = np.empty_like(z)
    owner = np.empty(z.size, dtype=int)
    for i in range(len(layers)):
        mask = (z >= edges[i]) & (z <= edges[i + 1] if i == len(layers)-1 else z < edges[i + 1])
        n[mask], owner[mask] = layers[i].n, i
    return z, n, owner


def slab_modes(layers, lambda_nm, pol="TE", n_modes=1):
    """Return bound scalar slab modes, ordered by effective index.

    ``pol`` is accepted as ``TE`` or ``TM``.  TM is the weak-guidance scalar
    approximation [A]; TE is the usual scalar slab equation [DR].
    """
    if pol.upper() not in ("TE", "TM"):
        raise ValueError("pol must be 'TE' or 'TM'")
    layers = list(layers)
    z, n, owner = _grid(layers)
    dz = (z[1] - z[0]) * 1e-3  # um
    k0 = 2 * pi / (lambda_nm * 1e-3)
    ni = n[1:-1]
    N = ni.size
    # The operator is tridiagonal, so the specialised solver is markedly
    # faster and avoids iterative convergence sensitivity. [DR]
    want = min(max(n_modes + 4, 6), N - 1)
    vals, vecs = eigh_tridiagonal(k0**2*ni**2 - 2/dz**2,
                                  np.ones(N-1)/dz**2, select="i",
                                  select_range=(N-want, N-1))
    order = np.argsort(vals)[::-1]
    nclad = max(layers[0].n, layers[-1].n)
    result = []
    for j in order:
        if vals[j] <= (k0*nclad)**2 * (1 + 1e-9):
            continue
        neff = float(np.sqrt(vals[j]) / k0)
        f = np.zeros_like(z); f[1:-1] = vecs[:, j]
        # Intensity-normalize: integral |E|^2 dz = 1.
        f /= np.sqrt(np.trapezoid(f*f, z))
        denom = np.trapezoid(f*f, z)
        bounds = np.r_[0.0, np.cumsum([L.thickness_nm for L in layers])]
        # Include both interface samples in each adjacent layer integral; the
        # field is continuous, and this avoids losing a thin sheet's endpoint.
        conf = {L.name: float(np.trapezoid(f[(z >= bounds[i]) & (z <= bounds[i+1])]**2,
                                           z[(z >= bounds[i]) & (z <= bounds[i+1])]) / denom)
                for i, L in enumerate(layers)}
        i2 = np.trapezoid(f**4, z)
        width = float(denom*denom / i2 * 1e-3)  # intensity effective width, um
        result.append(SlabMode(neff, z, f, layers, conf, width))
        if len(result) == n_modes:
            break
    if not result:
        raise ValueError("no bound mode: increase index contrast or cladding thickness")
    return result


def _etched(layers, etch_depth_nm):
    """Remove material from the top and terminate the etched stack in air [A]."""
    remain = float(etch_depth_nm)
    out = list(layers)
    while remain > 0 and out:
        top = out[-1]
        if remain >= top.thickness_nm:
            remain -= top.thickness_nm
            out.pop()
        else:
            out[-1] = Layer(top.name, top.n, top.thickness_nm-remain, top.is_dot)
            remain = 0
    if remain > 1e-9 or not out:
        raise ValueError("etch depth cannot remove the lower cladding")
    return out + [Layer("etched_air", 1.0, max(1000.0, layers[-1].thickness_nm))]


def effective_index_ridge(layers, lambda_nm, ridge_width_nm, etch_depth_nm, pol="TE"):
    """Effective-index ridge mode and separable intensity-effective area [DR].

    ``wx_um`` and ``wy_um`` are Gaussian-equivalent 1/e^2 *intensity*
    radii.  Thus a Gaussian has ``A_mode = pi * wx_um * wy_um``.
    """
    if ridge_width_nm <= 0:
        raise ValueError("ridge_width_nm must be positive")
    vertical = slab_modes(layers, lambda_nm, pol, 1)[0]
    outside_layers = _etched(list(layers), etch_depth_nm)
    try:
        outside = slab_modes(outside_layers, lambda_nm, pol, 1)[0]
        n_outside = outside.n_eff
    except ValueError:
        # A fully etched exterior need not itself guide.  Its scalar EIM index
        # is then the highest remaining slab material index [A].
        n_outside = max(outside_layers[0].n, outside_layers[-1].n)
        outside = SlabMode(n_outside, np.array([0., 1.]), np.array([1., 1.]),
                           outside_layers, {}, 1.0)
    # Lateral exterior is padded so its evanescent field reaches negligible size.
    pad = max(2000.0, ridge_width_nm)
    lateral_layers = [Layer("outside_left", n_outside, pad),
                      Layer("ridge", vertical.n_eff, ridge_width_nm),
                      Layer("outside_right", n_outside, pad)]
    lateral = slab_modes(lateral_layers, lambda_nm, pol, 1)[0]
    # Aeff=(int I)^2/int I^2; separability gives product of the two widths [DR].
    A = vertical.mode_width_um * lateral.mode_width_um
    # For I=exp(-2 x^2/w^2), (int I)^2/int I^2=sqrt(pi)*w.
    # Hence these are 1/e^2 intensity radii and Aeff=pi*wx*wy. [DR]
    wx = lateral.mode_width_um / np.sqrt(pi)
    wy = vertical.mode_width_um / np.sqrt(pi)
    return RidgeMode(lateral.n_eff, vertical, lateral, vertical.n_eff, n_outside,
                     float(A), float(wx), float(wy))


def beta_factor(A_mode_um2, lambda_nm, n_dot, n_g, position_factor=1.0):
    """``(F_wg, beta)`` for one guided direction pair, Lecamp 2007 [DR]."""
    if A_mode_um2 <= 0 or n_dot <= 0 or n_g <= 0 or not 0 <= position_factor <= 1:
        raise ValueError("invalid mode area, indices, or position_factor")
    lam = lambda_nm * 1e-3
    F = (3/(4*pi)) * (lam/n_dot)**2 * (n_g/n_dot) / A_mode_um2 * position_factor
    return float(F), float(F/(1+F))


def facet_transmission(n_eff, coating=None):
    """Normal-incidence Fresnel power transmission, or supplied coating override.

    A numeric ``coating`` or ``{'T_facet': value}`` is a measured/design power
    transmission [A], rather than a coating-transfer-matrix calculation.
    """
    if coating is not None:
        value = coating.get("T_facet", coating.get("transmission")) if isinstance(coating, dict) else coating
        if value is None or not 0 <= float(value) <= 1:
            raise ValueError("coating must be a transmission in [0, 1]")
        return float(value)
    R = ((n_eff-1)/(n_eff+1))**2
    return float(1-R)


def na_collection(wx_um, wy_um, lambda_nm, NA):
    """Circular-NA collection of an elliptical Gaussian far field [DR].

    ``wx_um`` and ``wy_um`` are 1/e^2 intensity radii.  The paraxial
    far-field intensity is ``exp(-2 tx^2/theta_x^2 - 2 ty^2/theta_y^2)``,
    with ``theta_i=lambda/(pi*w_i)``.  This integrates that distribution over
    the circular lens cone ``sqrt(tx^2 + ty^2) <= asin(NA)``.  Equal widths
    use the analytic circular result.
    """
    if wx_um <= 0 or wy_um <= 0 or not 0 <= NA <= 1:
        raise ValueError("widths must be positive and NA must be in [0, 1]")
    if NA == 0:
        return 0.0
    if NA == 1:
        return 1.0
    lam = lambda_nm * 1e-3
    theta_x, theta_y = lam/(pi*wx_um), lam/(pi*wy_um)
    theta_max = np.arcsin(NA)
    if np.isclose(theta_x, theta_y, rtol=1e-12, atol=0.0):
        return float(1.0 - np.exp(-2.0 * theta_max**2 / theta_x**2))

    # Integrate the normalized 2-D Gaussian over a disk.  Conditional on x,
    # the y integral is an erf; Gauss-Legendre quadrature avoids an arbitrary
    # Cartesian cutoff while retaining the round-aperture geometry. [DR]
    nodes, weights = np.polynomial.legendre.leggauss(96)
    x = theta_max * nodes
    y_extent = theta_max * np.sqrt(1.0 - nodes**2)
    px = (np.sqrt(2.0 / pi) / theta_x) * np.exp(-2.0 * x**2 / theta_x**2)
    integral = theta_max * np.sum(weights * px * erf(np.sqrt(2.0) * y_extent / theta_y))
    return float(np.clip(integral, 0.0, 1.0))


def na_collection_numeric(field_x, field_y, dx, dy, lambda_nm, NA):
    """Fourier-plane circular-NA collection for a separable scalar mode [DR].

    ``field_x``/``field_y`` are uniformly sampled complex field amplitudes;
    ``dx``/``dy`` are their sample spacings in um.  The result is normalized
    to the propagating forward hemisphere, so ``NA=1`` is exactly one.
    """
    if dx <= 0 or dy <= 0 or lambda_nm <= 0 or not 0 <= NA <= 1:
        raise ValueError("sample spacings and wavelength must be positive and NA must be in [0, 1]")
    fx_field = np.asarray(field_x)
    fy_field = np.asarray(field_y)
    if fx_field.ndim != 1 or fy_field.ndim != 1 or fx_field.size < 2 or fy_field.size < 2:
        raise ValueError("field_x and field_y must be one-dimensional sampled fields")
    if NA == 0:
        return 0.0
    spectrum = np.abs(np.outer(np.fft.fftshift(np.fft.fft(fy_field)),
                               np.fft.fftshift(np.fft.fft(fx_field))))**2
    fx = np.fft.fftshift(np.fft.fftfreq(fx_field.size, d=dx))
    fy = np.fft.fftshift(np.fft.fftfreq(fy_field.size, d=dy))
    transverse_na = (lambda_nm * 1e-3) * np.hypot(fy[:, None], fx[None, :])
    propagating = transverse_na <= 1.0
    denominator = spectrum[propagating].sum()
    if denominator == 0:
        raise ValueError("field spectra contain no propagating power")
    return float(spectrum[propagating & (transverse_na <= NA)].sum() / denominator)


def hkust_ridge_stack(lambda_nm=668.0):
    """AlInP / AlGaInP core with an InP dot in a GaAsP well [E].

    The one-nm well barriers are a thin effective active-region proxy [A].
    Every Layer carries its materials.py label (council review 2026-09-05
    item 5) so edge_emission's group-index finite difference can re-resolve
    it at a shifted wavelength; this is the sole reason hkust_ridge_stack no
    longer needs its own special-cased identity check there.
    """
    cl_label, core_label = "Al0.52In0.48P", "(Al0.50Ga0.50)0.51In0.49P"
    well_label, dot_label = "GaAs0.60P0.40", "InP"
    cl = refractive_index(cl_label, lambda_nm)
    core = refractive_index(core_label, lambda_nm)
    well = refractive_index(well_label, lambda_nm)
    dot = refractive_index(dot_label, lambda_nm)
    return [Layer("lower_cladding", cl, 1000, material=cl_label),
            Layer("core_lower", core, 148, material=core_label),
            Layer("well_lower", well, 1, material=well_label),
            Layer("dot", dot, 2, True, material=dot_label),
            Layer("well_upper", well, 1, material=well_label),
            Layer("core_upper", core, 148, material=core_label),
            Layer("upper_cladding", cl, 1000, material=cl_label)]


def _redispersed(stack, lambda_nm):
    """Rebuild `stack` with each wavelength-aware Layer's n re-evaluated at
    lambda_nm via materials.refractive_index (council review 2026-09-05 item
    5). A numeric-only layer (material is None) keeps its original n
    unchanged, and so does a labelled layer whose OWN material table does not
    reach lambda_nm (graceful per-layer degrade -- one thin layer's sparser
    table (e.g. a 1-2 nm well/dot sheet tabulated only at a few points)
    should not block re-solving the layers whose tables DO reach the shifted
    wavelength; whether the stack is dispersive AT ALL is decided by the
    caller from the un-shifted Layer.material fields, not from this
    best-effort per-layer attempt)."""
    out = []
    for x in stack:
        if x.material is not None:
            try:
                out.append(Layer(x.name, refractive_index(x.material, lambda_nm),
                                 x.thickness_nm, x.is_dot, x.material))
                continue
            except ValueError:
                pass
        out.append(x)
    return out


def facet_escape_fraction(T, R_back, alpha_cm, L_um, dot_position=0.5):
    """Ray-probability escape fraction through the front facet, single
    source of truth for ``edge_emission``'s facet term (peer-review pkg2
    facet fix, 2026-09-07, `.workers/specs/pr-pkg2-facet-fix.md`, correcting
    `.workers/specs/pr-pkg2-facet.md`'s checkpoint). Every independent
    self-check that needs to reproduce the facet factor (verify_waveguide.py,
    verify_device_rt.py, scripts/run_rt_edge.py's forward check) calls THIS
    function, or writes its own hand formula against these same numbers --
    never a formula copy-pasted out of ``edge_emission``'s source.

    Guided-mode emission from the dot splits 50/50 forward/backward into the
    mode [A]. The dot sits at fraction ``dot_position`` of ``L_um`` measured
    from the FRONT facet -- ``dot_position=0.5`` is a single-point
    mid-ridge assumption [A], NOT a uniform average over dot position along
    the ridge (a uniform average would weight terms by
    ``(1 - exp(-alpha_cm*1e-4*L_um)) / (alpha_cm*1e-4*L_um)`` instead).

    The forward half travels ``dot_position*L_um`` to the front facet and
    escapes with probability ``T`` or reflects with ``R_front = 1 - T``. The
    backward half travels ``(1-dot_position)*L_um`` to the back facet,
    reflects with ``R_back``, and then crosses the FULL ridge length
    ``L_um`` to reach the front facet -- not back to the dot's own
    position -- so for a mid-ridge dot (``dot_position=0.5``) the backward
    path to the front facet is ``L_um/2 + L_um``, one half-pass plus one
    full pass, not "one additional round trip of length 2*L_um" as an
    earlier version of this model mis-stated (a back-facet-to-front-facet
    leg is a single pass of ``L_um``, never a round trip). Either way the
    photon then either escapes (``T``) or reflects (``R_front``) for a
    further round trip of length ``2*L_um``. Summing the geometric series
    of round trips, with ``a = alpha_cm * 1e-4`` (per um) and
    ``x = dot_position``::

        prop_rt = exp(-2*a*L_um)
        eta_facet = 0.5*T*exp(-a*x*L_um) * (1 + R_back*exp(-2*a*(1-x)*L_um))
                    / (1 - R_back*R_front*prop_rt)

    At the mid-point ``x=0.5`` the returning first-arrival factor is
    ``exp(-a*L_um) == sqrt(prop_rt)``, NOT ``prop_rt`` itself -- the
    checkpoint version of this formula used ``prop_rt`` there, which
    silently gave the backward half an extra round trip of attenuation it
    had not yet taken.

    [DR] Coldren, Corzine & Masanovic, *Diode Lasers and Photonic
    Integrated Circuits*, 2nd ed. (2012), ch. 2, mirror-loss/escape-fraction
    treatment, extended here to a finite, possibly-uncoated back facet, a
    dot position other than the antinode mid-point, and finite ridge loss
    rather than a lossless HR mirror. Textbook cross-check at alpha_cm=0,
    R_back=0.95, T=0.719371: this series gives 0.9563, vs. Coldren, Corzine
    & Masanovic's lossless front-facet fraction
    ``F1 = (1-R_front)*sqrt(R_back) / ((1-R_front)*sqrt(R_back) +
    (1-R_back)*sqrt(R_front))`` = 0.9636 -- the two are close but not
    identical (different physical quantities: this is a ray-probability
    escape fraction for spontaneous emission, F1 is a lasing-cavity
    front-facet output-coupling fraction); report both, do not conflate
    them.
    """
    if not 0.0 <= T <= 1.0:
        raise ValueError("T must be in [0, 1]")
    if not 0.0 <= R_back <= 1.0:
        raise ValueError("R_back must be in [0, 1]")
    if not 0.0 <= dot_position <= 1.0:
        raise ValueError("dot_position must be in [0, 1]")
    R_front = 1.0 - T
    a = alpha_cm * 1e-4  # per um
    x = dot_position
    prop_rt = exp(-2.0 * a * L_um)
    denom = 1.0 - R_back * R_front * prop_rt
    if denom <= 1e-15:
        if R_back == 1.0 and alpha_cm * L_um == 0.0:
            # Lossless perfect-mirror limit: with zero ridge loss and a
            # fully reflective back facet, the front facet is the ONLY
            # loss channel in the cavity, so by power conservation every
            # photon eventually escapes forward after enough round trips,
            # regardless of T -- this is the well-defined analytic limit
            # of the ray series above (numerator and denominator both -> 0
            # together, ratio -> 1), not a fallback for an undefined case.
            return 1.0
        raise ValueError(
            "degenerate facet cavity: R_back * R_front * exp(-2*alpha_cm*1e-4*L_um) "
            ">= 1, the round-trip geometric series does not converge")
    return (0.5 * T * exp(-a * x * L_um)
            * (1.0 + R_back * exp(-2.0 * a * (1.0 - x) * L_um)) / denom)


def edge_emission(stack, ridge_width_nm, etch_depth_nm, lambda_nm, L_um, NA,
                  alpha_cm=5.0, R_back=None, coating=None, na_method="gaussian",
                  dot_position=0.5):
    """Calculate the collected forward edge-emission probability.

    Facet model (peer-review finding 3, `.workers/review/peer-review-triage.md`,
    2026-09-07; fix pass `.workers/specs/pr-pkg2-facet-fix.md`, 2026-09-07): a
    single ray-probability model, continuous in ``R_back``, replaces the old
    two-branch geometric/escape-rate switch. That switch stepped by +16.33%
    at T_facet=0.719371 crossing R_back=0 (0.5*T=0.3596855 vs the
    escape-rate limit T/(T+1)=0.4183916), and its comment claimed T/(T+1)
    was "strictly below" 0.5*T when in fact
    T/(T+1) - 0.5*T = T(1-T)/(2(T+1)) > 0, i.e. strictly ABOVE.

    The facet factor itself is computed by :func:`facet_escape_fraction`
    (single source of truth); see its docstring for the derivation,
    ``dot_position`` semantics, and the closed-form series. ``total`` below
    no longer applies a further single-pass ``prop`` factor on top of it --
    that would double-count loss already folded into the series. The
    ``eta_prop`` field of the returned :class:`EdgeResult` still reports the
    old single-full-length diagnostic ``exp(-alpha_cm * L_um * 1e-4)`` for
    continuity with existing callers/reports, but it is informational only
    -- it is not multiplied into ``eta_total``.

    ``R_back=None`` means "no HR coating, cleaved back facet": it resolves
    to the UNCOATED Fresnel reflectivity of the bare facet (the same
    formula ``facet_transmission`` uses with no ``coating`` override), not
    to ``1 - T_facet`` of whatever ``coating``-overridden front transmission
    is in effect -- a coated front facet does not change what the back
    facet's own, separately uncoated, Fresnel reflectivity is. When no
    ``coating`` override is given at all, the uncoated reflectivity and
    ``1 - T_facet`` are the same number, so this is a strict generalization
    of the previous (uncoated-only) behaviour, not a change to it. With this
    resolution the model is continuous through ``R_back -> 0`` by
    construction (checked in verify/verify_waveguide.py).
    """
    if L_um < 0 or alpha_cm < 0 or R_back is not None and not 0 <= R_back <= 1:
        raise ValueError("invalid length, loss, or back reflectance")
    if na_method not in ("gaussian", "numeric"):
        raise ValueError("na_method must be 'gaussian' or 'numeric'")
    mode = effective_index_ridge(stack, lambda_nm, ridge_width_nm, etch_depth_nm)
    notes = ["[DR] scalar effective-index mode",
             "[DR] mode widths are Gaussian-equivalent 1/e^2 intensity radii; A_mode = pi wx wy",
             f"[DR] NA collection method: {na_method}", "[A] antinode dipole"]
    # Group index (council review 2026-09-05 item 5): a 5 nm central finite
    # difference on n_eff, re-solving every WAVELENGTH-AWARE layer (one that
    # carries a materials.py label) at lambda +/- 5 nm -- generalizes the
    # former hkust_ridge_stack-only special case to any stack built with
    # dispersive Layers (device.py's _resolve_edge included). A stack with no
    # labelled layer (numeric-only n, e.g. verify's synthetic slabs) keeps the
    # n_g = n_eff fallback and says so in notes, exactly as before.
    ng = mode.n_eff
    dispersive = any(x.material is not None for x in stack)
    if dispersive:
        mminus = mplus = None
        try:
            stack_minus = _redispersed(stack, lambda_nm - 5.0)
            mminus = effective_index_ridge(stack_minus, lambda_nm - 5.0,
                                           ridge_width_nm, etch_depth_nm)
        except ValueError:  # no bound mode at the shifted wavelength [A]
            pass
        try:
            stack_plus = _redispersed(stack, lambda_nm + 5.0)
            mplus = effective_index_ridge(stack_plus, lambda_nm + 5.0,
                                          ridge_width_nm, etch_depth_nm)
        except ValueError:  # no bound mode at the shifted wavelength [A]
            pass
        if mminus is not None and mplus is not None:
            slope = (mplus.n_eff - mminus.n_eff) / 10.0
            notes.append("[DR] n_g: 5 nm central finite difference on wavelength-aware layers")
            ng = mode.n_eff - lambda_nm * slope
        elif mplus is not None:  # one-sided forward difference [DR]
            slope = (mplus.n_eff - mode.n_eff) / 5.0
            notes.append("[DR] n_g: 5 nm forward finite difference (blue endpoint out of table)")
            ng = mode.n_eff - lambda_nm * slope
        elif mminus is not None:  # one-sided backward difference [DR]
            slope = (mode.n_eff - mminus.n_eff) / 5.0
            notes.append("[DR] n_g: 5 nm backward finite difference (red endpoint out of table)")
            ng = mode.n_eff - lambda_nm * slope
        else:
            notes.append("[A] n_g fallback to n_eff: both dispersion endpoints out of table")
    else:
        notes.append("[A] n_g fallback to n_eff: no wavelength-aware (materials-labelled) layer")
    dots = [x for x in stack if x.is_dot]
    if not dots:
        raise ValueError("stack must mark one Layer as is_dot=True")
    gamma = sum(mode.vertical.confinement.get(x.name, 0.0) for x in dots)
    n_dot = dots[0].n
    # Position factor (council review 2026-09-05 item 4): |E(z_dot)|^2 /
    # max|E|^2 evaluated at the DOT LAYER'S OWN centre z (thickness-weighted
    # mean centre if more than one dot layer), not at wherever the field
    # happens to peak -- the latter is identically 1 by construction and
    # cannot express a dot placed off the vertical antinode.
    bounds = np.r_[0.0, np.cumsum([x.thickness_nm for x in stack])]
    centres = [0.5 * (bounds[i] + bounds[i + 1]) for i, x in enumerate(stack) if x.is_dot]
    z_dot = float(np.average(centres, weights=[x.thickness_nm for x in dots]))
    field_at_dot = float(np.interp(z_dot, mode.vertical.z_nm, mode.vertical.field))
    pos = min(1.0, float(field_at_dot ** 2 / np.max(mode.vertical.field ** 2)))
    F, beta = beta_factor(mode.A_mode_um2, lambda_nm, n_dot, ng, pos)
    T = facet_transmission(mode.n_eff, coating)
    # Peer-review pkg2 fix (2026-09-07, item 4): R_back=None resolves to the
    # UNCOATED Fresnel reflectivity of the bare back facet, computed the
    # same way facet_transmission computes the uncoated T -- NOT to
    # 1-T_facet of a coating-overridden front, which would wrongly let a
    # front-facet coating change the back facet's own physics.
    if R_back is None:
        T_uncoated = facet_transmission(mode.n_eff, None)
        R_back_eff = 1.0 - T_uncoated
        notes.append("[A] R_back=None resolved to the uncoated Fresnel "
                     "reflectivity of the bare back facet (no HR coating: "
                     "cleaved back facet, same semiconductor/air interface "
                     "physics as an uncoated front, independent of any "
                     "`coating` override on the front transmission)")
    else:
        R_back_eff = R_back
    facet_factor = facet_escape_fraction(T, R_back_eff, alpha_cm, L_um, dot_position)
    notes.append(
        "[DR] ray-probability escape fraction (facet_escape_fraction, "
        "convention 'ray-series-midpoint' at the dot_position default), "
        "continuous in R_back: eta_facet = 0.5 * T_facet * "
        "exp(-a*dot_position*L_um) * (1 + R_back * "
        "exp(-2*a*(1-dot_position)*L_um)) / (1 - R_back * R_front * "
        "exp(-2*a*L_um)), a = alpha_cm*1e-4, R_front = 1 - T_facet "
        "(Coldren, Corzine & Masanovic, Diode Lasers and Photonic "
        "Integrated Circuits, 2nd ed. (2012), ch. 2, mirror-loss/"
        "escape-fraction; propagation is folded entirely into this series, "
        "so it is NOT applied again as a separate factor below)")
    prop = exp(-alpha_cm * (L_um * 1e-4))  # diagnostic only; see docstring
    if na_method == "gaussian":
        eta_na = na_collection(mode.wx_um, mode.wy_um, lambda_nm, NA)
    else:
        dx = (mode.lateral.z_nm[1] - mode.lateral.z_nm[0]) * 1e-3
        dy = (mode.vertical.z_nm[1] - mode.vertical.z_nm[0]) * 1e-3
        eta_na = na_collection_numeric(mode.lateral.field, mode.vertical.field,
                                       dx, dy, lambda_nm, NA)
    total = beta * facet_factor * eta_na
    return EdgeResult(mode.n_eff, float(ng), float(gamma), mode.A_mode_um2, F, beta, T,
                      prop, eta_na, float(total), notes)
