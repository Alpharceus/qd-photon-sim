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
from math import erf, exp, pi, sin
import numpy as np
from scipy.linalg import eigh_tridiagonal

from .materials import refractive_index


@dataclass(frozen=True)
class Layer:
    """A vertical layer.  First and last layers also define the exterior [A]."""
    name: str
    n: float
    thickness_nm: float
    is_dot: bool = False


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
    """Effective-index ridge mode and separable intensity-effective area [DR]."""
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
    # Gaussian-equivalent 1/e^2 intensity radii: Aeff=sqrt(2*pi)*w per axis.
    wx = lateral.mode_width_um / np.sqrt(2*pi)
    wy = vertical.mode_width_um / np.sqrt(2*pi)
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
    """Gaussian far-field lens collection estimate [E]."""
    if wx_um <= 0 or wy_um <= 0 or not 0 <= NA <= 1:
        raise ValueError("widths must be positive and NA must be in [0, 1]")
    if NA == 0:
        return 0.0
    if NA == 1:
        return 1.0
    lam = lambda_nm * 1e-3
    tx, ty = lam/(pi*wx_um), lam/(pi*wy_um)
    # Clamp angles: the paraxial Gaussian expression becomes a full hemisphere.
    return float(erf(NA/(np.sqrt(2)*sin(min(tx, pi/2)))) *
                 erf(NA/(np.sqrt(2)*sin(min(ty, pi/2)))))


def hkust_ridge_stack(lambda_nm=668.0):
    """AlInP / AlGaInP core with an InP dot in a GaAsP well [E].

    The one-nm well barriers are a thin effective active-region proxy [A].
    """
    cl = refractive_index("Al0.52In0.48P", lambda_nm)
    core = refractive_index("(Al0.50Ga0.50)0.51In0.49P", lambda_nm)
    well = refractive_index("GaAs0.60P0.40", lambda_nm)
    dot = refractive_index("InP", lambda_nm)
    return [Layer("lower_cladding", cl, 1000), Layer("core_lower", core, 148),
            Layer("well_lower", well, 1), Layer("dot", dot, 2, True),
            Layer("well_upper", well, 1), Layer("core_upper", core, 148),
            Layer("upper_cladding", cl, 1000)]


def edge_emission(stack, ridge_width_nm, etch_depth_nm, lambda_nm, L_um, NA,
                  alpha_cm=5.0, R_back=None, coating=None):
    """Calculate the collected forward edge-emission probability.

    Two facets share emission equally unless ``R_back`` is supplied; then the
    stated escape-rate approximation is ``(1-Rf)/[(1-Rf)+(1-Rb)]`` [A].
    """
    if L_um < 0 or alpha_cm < 0 or R_back is not None and not 0 <= R_back <= 1:
        raise ValueError("invalid length, loss, or back reflectance")
    mode = effective_index_ridge(stack, lambda_nm, ridge_width_nm, etch_depth_nm)
    # Material tabulations are sparse; use a 1-nm symmetric difference only
    # when this is the supplied named HKUST stack, otherwise constant-n [A].
    ng = mode.n_eff
    if [x.name for x in stack] == [x.name for x in hkust_ridge_stack(lambda_nm)]:
        try:
            mminus = effective_index_ridge(hkust_ridge_stack(lambda_nm-1), lambda_nm-1,
                                            ridge_width_nm, etch_depth_nm)
            try:
                mplus = effective_index_ridge(hkust_ridge_stack(lambda_nm+1), lambda_nm+1,
                                               ridge_width_nm, etch_depth_nm)
                slope = (mplus.n_eff-mminus.n_eff)/2
            except ValueError:  # use a one-sided endpoint difference [DR]
                slope = mode.n_eff-mminus.n_eff
            ng = mode.n_eff - lambda_nm*slope
        except ValueError:  # endpoint outside an input material table [A]
            pass
    dots = [x for x in stack if x.is_dot]
    if not dots:
        raise ValueError("stack must mark one Layer as is_dot=True")
    gamma = sum(mode.vertical.confinement.get(x.name, 0.0) for x in dots)
    n_dot = dots[0].n
    # Dot is assumed at the ridge centre and vertical antinode [A].
    pos = float((mode.vertical.field[np.argmax(np.abs(mode.vertical.field))] ** 2) /
                np.max(mode.vertical.field**2))
    F, beta = beta_factor(mode.A_mode_um2, lambda_nm, n_dot, ng, pos)
    T = facet_transmission(mode.n_eff, coating)
    front = 0.5 if R_back is None else T / (T + (1-R_back))
    prop = exp(-alpha_cm * (L_um*1e-4))
    eta_na = na_collection(mode.wx_um, mode.wy_um, lambda_nm, NA)
    total = beta * front * prop * eta_na
    return EdgeResult(mode.n_eff, float(ng), float(gamma), mode.A_mode_um2, F, beta, T,
                      prop, eta_na, float(total), ["[DR] scalar effective-index mode", "[E] Gaussian NA", "[A] antinode dipole"])
