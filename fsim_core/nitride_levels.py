"""Finite-well InGaN/GaN dot levels.

This is a single-band, anisotropic-effective-mass, separable disk model [E],
with thick GaN reservoirs and a real-energy (not tunnelling) retention
treatment [A].  Consequently a bound state in a tilted well may still
tunnel; Arrhenius retention is optimistic in that situation.  Material
response is supplied by :mod:`nitride_materials`.

FIELD GEOMETRY.  Fixed-displacement polarization produces a uniform field
inside an isolated dot and zero field outside it, but electrostatic potential
remains continuous: the exterior on each side is held at the value reached at
that face.  The two flat plateaus therefore differ by F*h in the unscreened
parallel-plate limit [DR], Bernardini & Fiorentini, phys. stat. sol. (b) 216,
391 (1999), Eqs. 7-8 [V].  Escape is judged against the lower plateau on the
side toward which each carrier drifts.  `screening_fraction` is an explicit
[A] parameter and its default remains zero.  That unscreened isolated-dot
limit is a pessimistic confinement bound: electrically injected structures
can partially screen the polarization field, as indicated by the <2 meV
current shift in Deshpande et al., Nat. Commun. 4, 1675 (2013), p.5 [V], and
the bias-dependent QCSE in Zhang et al., APL 108, 153102 (2016), Fig. 5 [V].

The opt-in ``qw_fluctuation`` geometry is a same-composition local thickness
fluctuation.  Its lateral depth is the difference between independently
solved local-column and surrounding-QW subband edges [A], not the GaN band
offset.  Carrier edge energies are measured upward from the strained InGaN
band edge; ``reservoir_energy_eV`` is the free electron-hole continuum edge
and deliberately contains no dot Coulomb correction.  The fluctuation model
also assumes the local dot column and the surrounding well are centered on
the SAME z=0 midplane [A]; in the unscreened tilted limit (large intrinsic
or external field, small confinement offset) that placement assumption
becomes the entire predicted lateral depth, which collapses to the
mass-independent geometric estimate |F|*(H-w)/2 (F in eV/nm, H and w in nm)
rather than a quantum-confinement difference -- a one-sided (asymmetric)
fluctuation would give a very different value (0 up to |F|*(H-w)).  This is
disclosed per-row in ``approximation_error``.

ORIENTATION-DEPENDENT MASSES.  For c_plane and semipolar_11_22 growth the
growth (z) axis is taken as the crystal c-axis, so the z-confinement mass is
me_z/mh_z and the lateral (radial) mass is me_xy/mh_xy, as before.  For
m_plane/a_plane (nonpolar) growth the growth axis is perpendicular to c, so
the z-confinement mass is instead the crystal's PERPENDICULAR (xy)
component and the c-axis component is used for the lateral direction [A;
the true nonpolar lateral direction is itself anisotropic -- one in-plane
direction lies along c (mass mh_z) and the other is again perpendicular
(mass mh_xy) -- this single-band radial model cannot represent that split
and uses the c-axis value for both].  Rinke et al., PRB 77, 075202 (2008),
Table V.  At orientation=a_plane, height_nm=3 (otherwise-default system)
this swap moves E_X by about +32 meV and the hole escape depth by about
-24 meV relative to (incorrectly) using c-plane masses at every
orientation.  semipolar_11_22 keeps c-plane masses [A; orientation-
dependent valence reordering, Schade et al., phys. status solidi (b)
(2011), is not modelled].

RESERVOIR KIND AND ENERGY.  For an isolated dot (no wl_thickness_nm)
``reservoir_kind`` is always ``'gan_barrier'`` and ``reservoir_energy_eV``
is the bulk GaN Varshni gap at the row's T_K minus the same 25 meV [A]
localization/Urbach offset that ``fsim_core.device._nitride_reservoir_energy_eV``
applies on its own zero-wetting-layer branch, so the two independently
computed values agree numerically (about 3.41 eV at 300 K); the field is
therefore orientation/height/radius/screening independent, as a bulk
material edge must be.  device.py does NOT read this field for isolated-dot
rows -- it keeps its own bulk-edge helper as the value it reports -- so
this field is provided for documentation/consistency (and for any future or
external consumer), not as the number device.py's isolated-dot rows
currently forward.  For a qw_fluctuation dot each carrier independently
selects its own escape channel, min(surrounding-QW subband edge, GaN
plateau); ``reservoir_electron_edge_meV``/``reservoir_hole_edge_meV`` report
that SELECTED per-carrier edge (not the raw surrounding-well edge).
``reservoir_kind`` is ``'ingan_qw'`` only when BOTH carriers select the
well, ``'gan_barrier'`` only when BOTH select the GaN plateau, and
``'mixed'`` when the two carriers select different channels (e.g. electron
confined by the well, hole escaping into the wider GaN barrier) -- a
spatially indirect situation no single reservoir label describes exactly;
``reservoir_energy_eV`` still reports gap+eth+hth in that case (device.py
only reads reservoir_kind as a label and reservoir_energy_eV as a number,
never branches on the kind string, so a 'mixed' row evaluates normally).
"""
from dataclasses import dataclass
from functools import lru_cache
import math
import numpy as np
from scipy.linalg import eigh_tridiagonal
from .nitride_materials import binary, ingaN, band_edges, polarization_field, orientation_factor, KB_EV
from .dot_levels import finite_disk_2d

_HBAR2_2M0 = 0.0380998212 # eV nm2 [DR] CODATA 2018 constants
_HBAR2_2M0_MEV = 38.0998212 # meV nm2, same constant [DR] CODATA 2018
_HBAR_SI = 1.054571817e-34 # J s [V] CODATA 2018
_M0 = 9.1093837015e-31 # kg [V] CODATA 2018
_KB_SI = 1.380649e-23 # J K-1 [V] SI 2019
_HC_EV_NM = 1239.841984 # h c [V] CODATA 2018, eV nm
_E2_4PIEPS0_EV_NM = 1.439964 # e^2/(4 pi eps0) [DR] CODATA 2018, eV nm
TAIL_CONVERGENCE_TOL_MEV = 1e-3 # [E] padding-doubling bound-state convergence threshold

# [E] Documented, NOT runtime-computed, default-grid interior discretization
# bound (hardening round after commit 17a333f, Opus re-review finding).
# electron_padding_delta_meV/hole_padding_delta_meV below only re-solve with
# exterior_nm doubled at the SAME z_points, so they measure padding
# (exterior-domain) sufficiency, not interior mesh resolution, and read
# close to 0 meV at these heights even though a doubling of BOTH z_points
# and exterior_nm together (the check the verifier's convergence matrix
# actually performs) moves E_X well past the 0.5 meV convergence budget at
# the DEFAULT grid (z_points=1201, exterior_nm=45). Measured once, offline,
# at commit 17a333f (c-plane, R=17.5 nm, screening_fraction=0, T=300 K);
# doubling both controls at construction cost for every default-settings
# levels() call is explicitly NOT an acceptable way to close this finding.
# See verify_nitride_geometry.py's own (h, R, orientation, geometry_type,
# screening) convergence matrix for the general sweep this table summarizes.
DEFAULT_GRID_EX_BOUND_MEV_17A333F = {5.: 0.438, 6.: 0.525, 7.: 0.612} # h_nm -> meV [E]
DEFAULT_GRID_EX_RESIDUAL_MEV_17A333F_H7 = 1.2 # [E] residual at h=7nm after one further doubling
_GRID_NOTE = ('default grid (z_points=1201, exterior_nm=45) has a height-dependent interior-'
    'resolution discretization error NOT captured by padding_delta (exterior-only doubling): '
    'measured (not recomputed per call) at commit 17a333f, c-plane, R=17.5 nm, '
    'screening_fraction=0, doubling z_points AND exterior_nm together moves E_X by '
    + repr(DEFAULT_GRID_EX_BOUND_MEV_17A333F) + ' meV (h_nm -> meV) with a ~'
    + repr(DEFAULT_GRID_EX_RESIDUAL_MEV_17A333F_H7)
    + ' meV residual at h=7 nm even after that further doubling; pass explicit '
    'z_points/exterior_nm for a tighter bound')

# Source-transcription targets: named module constants (compared, in the
# verifier, to INDEPENDENTLY typed literals so the check is a real
# transcription cross-check rather than a literal-equals-itself tautology).
DESHPANDE2013_HEIGHT_NM = 2.0         # [V] Deshpande et al., Nat. Commun. 4, 1675 (2013), p.2, Fig. 1d
DESHPANDE2013_DIAMETER_NM = 25.0      # [V] ibid., p.2, Fig. 1c (25 +/- 5 nm)
DESHPANDE2013_X_IN = 0.25             # [V] ibid., p.2 (single-wire EDX)
DESHPANDE2013_X_NM = 436.56           # [V] ibid., p.2-3, Fig. 3c (single-dot EL, X line)
DESHPANDE2013_XX_ANTIBIND_MEV = -10.0 # [V] ibid., p.2-3, Fig. 3c (XX above X: antibinding)
ZHANG2013_HEIGHT_NM = 3.0             # [V] Zhang et al., APL 103, 192114 (2013), Fig. 2
ZHANG2013_DIAMETER_NM = 29.0          # [V] ibid., Fig. 2 (D = 29 nm pillar)
ZHANG2013_X_IN = 0.15                 # [V] ibid., Fig. 2 (In0.15Ga0.85N)
ZHANG2013_ZPL_EV = 2.95               # [V] ibid., Fig. 2 (10 K)
ZHANG2013_LIFETIME_NS = 3.40          # [V] ibid., Fig. 2 (mono-exp lifetime, 10 K)
ZHANG2016_SLOPE_MEV_PER_V = -10.0     # [V] Zhang et al., APL 108, 153102 (2016), Fig. 5 (below 2 V)
INGAN_DOT_LIFETIME_RANGE_NS = (1.0, 10.0) # [E] digest class range, nitride_digests.md Sec.11 item 5
TAU_RAD0_DEFAULT_NS = 1.3              # [E] Deshpande et al., APL 105, 141109 (2014), abstract:
                                        # RT electrically-driven recombination lifetime 1.3+/-0.3 ns;
                                        # transferred here as a field-free radiative-lifetime DEFAULT,
                                        # not a verified measurement of this planar design's own
                                        # baseline (their number mixes radiative and any non-radiative
                                        # loss at RT) -- [E], overridable by callers.

@dataclass(frozen=True)
class NitrideDotSystem:
    height_nm: float = 3.0; radius_nm: float = 10.0; x_in: float = .25
    wl_thickness_nm: float = 0.0; strain_fraction: float = 1.0
    screening_fraction: float = 0.0; external_field_kVcm: float = 0.0
    vbo_InN_GaN_eV: float = 1.15; strain_c_fraction: float = .7
    orientation: str = "c_plane"; polarization_factor: object = None
    shape: str = "disc"; top_radius_fraction: float = 1.0
    geometry_type: str = "isolated_dot"; shape_height_fraction: object = None

@dataclass(frozen=True)
class NitrideLevels:
    """Resolved nitride levels and reservoirs.

    ``reservoir_*`` fields describe the selected per-carrier escape channels.
    ``optical_reservoir_energy_eV``/``optical_reservoir_kind`` instead describe
    the continuum used for flat optical-background acceptance: the surrounding
    InGaN QW ground subbands for QW fluctuations, or GaN minus 25 meV [A] for
    isolated dots.  The escape-rate DOS is a bulk-channel approximation [A].
    """
    E_X_eV: float; lambda_nm: float; electron_bound: bool; hole_bound: bool
    overlap_sq: float; field_kVcm: float; E_e_meV: float; E_h_meV: float
    dE_e_meV: float; dE_h_meV: float; dE_pair_meV: object
    sp_split_e_meV: float; sp_split_h_meV: float; m_e_matrix_xy: float
    m_h_matrix_xy: float; valid: bool; invalid_reasons: tuple; provenance: str
    electron_in_dot_probability: float = float('nan')
    hole_in_dot_probability: float = float('nan')
    electron_padding_delta_meV: float = float('nan')
    hole_padding_delta_meV: float = float('nan')
    electron_exterior_left_eV: float = float('nan')
    electron_exterior_right_eV: float = float('nan')
    hole_exterior_left_eV: float = float('nan')
    hole_exterior_right_eV: float = float('nan')
    reservoir_kind: str = 'gan_barrier'
    reservoir_energy_eV: float = float('nan')
    reservoir_electron_edge_meV: float = float('nan')
    reservoir_hole_edge_meV: float = float('nan')
    effective_height_nm: float = float('nan')
    effective_radius_nm: float = float('nan')
    shape_volume_nm3: float = float('nan')
    geometry_approximation: str = ''
    approximation_error: str = ''
    polarization_factor_used: float = float('nan')
    physical_height_nm: float = float('nan')
    physical_radius_nm: float = float('nan')
    aspect_ratio: float = float('nan')
    optical_reservoir_energy_eV: float = float('nan')
    optical_reservoir_kind: str = 'gan_barrier'
    escape_e_matrix_xy: float = float('nan')
    escape_h_matrix_xy: float = float('nan')

def _validate(s):
    bad=[]
    for n in ('height_nm','radius_nm','x_in','wl_thickness_nm','strain_fraction',
              'screening_fraction','strain_c_fraction','external_field_kVcm','vbo_InN_GaN_eV'):
        v=getattr(s,n)
        if not math.isfinite(v): bad.append(n+' is not finite')
    if s.height_nm<=0: bad.append('height_nm must be positive')
    if s.radius_nm<=0: bad.append('radius_nm must be positive')
    if not 0<=s.x_in<=1: bad.append('x_in must be in [0,1]')
    if not 0<=s.wl_thickness_nm<s.height_nm: bad.append('wl_thickness_nm must be >=0 and < height_nm')
    if not 0<=s.strain_fraction<=1: bad.append('strain_fraction must be in [0,1]')
    if not 0<=s.screening_fraction<=1: bad.append('screening_fraction must be in [0,1]')
    if s.shape not in ('disc','lens','truncated_cone'): bad.append('shape is invalid')
    if s.geometry_type not in ('isolated_dot','qw_fluctuation'): bad.append('geometry_type is invalid')
    if not isinstance(s.top_radius_fraction,(int,float)) or isinstance(s.top_radius_fraction,bool) or not math.isfinite(s.top_radius_fraction) or not 0<=s.top_radius_fraction<=1:
        bad.append('top_radius_fraction must be finite in [0,1]')
    elif s.shape != 'truncated_cone' and s.top_radius_fraction != 1.:
        bad.append('top_radius_fraction is only supported for truncated_cone')
    if s.shape_height_fraction is not None:
        q=s.shape_height_fraction
        if not isinstance(q,(int,float)) or isinstance(q,bool) or not math.isfinite(q):
            bad.append('shape_height_fraction must be finite')
        else:
            native = .5 if s.shape == 'lens' else ((1+s.top_radius_fraction+s.top_radius_fraction*s.top_radius_fraction)/3. if s.shape == 'truncated_cone' else 1.)
            if s.shape == 'disc' or not native<=q<=1.: bad.append('shape_height_fraction is unsupported or outside [native,1]')
    if s.geometry_type == 'qw_fluctuation':
        if s.shape != 'disc' or s.shape_height_fraction is not None: bad.append('qw_fluctuation requires disc shape and no shape_height_fraction')
        if not 0 < s.wl_thickness_nm < s.height_nm: bad.append('qw_fluctuation requires 0 < wl_thickness_nm < height_nm')
    elif s.wl_thickness_nm != 0: bad.append('wl_thickness_nm is only supported for qw_fluctuation')
    try: orientation_factor(s.orientation,s.polarization_factor)
    except ValueError as exc: bad.append(str(exc))
    return bad

def _geometry(s):
    """Effective separable confinement mapping [A/DR elementary volume integration].

    shape='lens' is specifically a PARABOLOID of revolution (h_eff=H/2 by
    volume-preserving construction); a spherical-cap lens of the same H,R
    has h_eff = H/2 + H**3/(6*R*R) instead (e.g. 3.25 vs 3.00 nm at H=6,
    R=12 nm) [A elementary volume integration]. That profile choice --
    paraboloid vs spherical cap, and disc/cone vs any real 3-D solve -- is
    part of the UNQUANTIFIED shape systematic error disclosed in
    `approximation_error`, not a rigorous bound on it.
    """
    if s.shape == 'disc': native=1.; label='full-height disc'
    elif s.shape == 'lens': native=.5; label='paraboloidal lens volume mapping'
    else:
        t=s.top_radius_fraction; native=(1+t+t*t)/3.; label='truncated-cone volume mapping'
    frac=native if s.shape_height_fraction is None else s.shape_height_fraction
    return s.height_nm*frac, s.radius_nm, math.pi*s.radius_nm*s.radius_nm*s.height_nm*native, label

def _growth_masses(mat, orientation):
    """Return (me_growth, me_lateral, mh_growth, mh_lateral) for `mat` along
    the actual growth axis of `orientation`. See the module docstring
    ("ORIENTATION-DEPENDENT MASSES") for the physical justification and
    magnitude of the m/a-plane swap; semipolar_11_22 keeps c-plane masses."""
    if orientation in ('m_plane', 'a_plane'):
        return mat.me_xy, mat.me_z, mat.mh_xy, mat.mh_z
    return mat.me_z, mat.me_xy, mat.mh_z, mat.mh_xy

def _z_grid(half, pad, target_n):
    """Node grid with +/-half landing EXACTLY on a node (finite-volume cell
    face), so the mass/potential step at the dot interface is not smeared
    to first order by an arbitrary offset between a linspace node and the
    physical edge. `target_n` sets the INTERIOR half-well resolution
    (n_half = target_n/2, floored at 200): the well is nm-scale while
    `pad` (exterior_nm) is tens of nm, so splitting a single uniform dz
    proportionally across the whole padded domain (the previous scheme)
    starves the interior of points for a thin dot -- e.g. z_points=1201,
    exterior_nm=45 gave only ~7 half-well nodes at height_nm=1, a first-
    order interface error far above the 0.5 meV budget. Fixing dz from the
    interior requirement and extending the SAME dz into the pad instead
    keeps the interface well resolved regardless of pad size."""
    n_half = max(200, int(round(target_n / 2.0)))
    dz = half / n_half
    n_pad = max(1, int(round(pad / dz)))
    idx = np.arange(-(n_half + n_pad), n_half + n_pad + 1)
    return idx * dz, dz, n_half

def _z_potential(height, barrier, field, sign, z):
    """Continuous tilted-well potential in eV on coordinates `z` in nm.

    The electrostatic term is linear inside the dot and is clipped to its
    interface value outside.  Thus the field vanishes in each exterior while
    the face potentials retain the F*h offset. [DR] Bernardini & Fiorentini,
    pss(b) 216, 391 (1999), Eqs. 7-8, isolated-slab parallel-plate limit.
    """
    half = height / 2.
    z = np.asarray(z, dtype=float)
    inside = np.abs(z) <= half + 1e-12
    tilt = sign * field * 1e-4 * np.clip(z, -half, half)
    return np.where(inside, 0., barrier) + tilt


def _z_state(height, barrier, md, mb, field, sign, n=1201, pad=45.):
    """BDD finite-volume state; field is kV/cm, sign is carrier charge sign.

    `field` tilts the potential inside the dot; each exterior is flat at the
    adjacent face value, so V is electrostatically continuous apart from the
    material band offset itself.  The exterior padding is only a numerical
    convergence control.  The external-field contribution is folded into
    the same clipped potential [A]; modeling a diode-scale depletion ramp is
    outside this nm-scale solver.  BenDaniel & Duke, PR 152, 683 (1966) [V]
    supplies the flux-matching boundary condition (mass-weighted hopping).
    """
    half = height / 2.
    z, dz, n_half = _z_grid(half, pad, n)
    n_pts = len(z)
    node_idx = np.round(z / dz).astype(int)
    inside = np.abs(node_idx) <= n_half
    mass = np.where(inside, md, mb)
    # e*(kV/cm)*(nm) = 1e-4 eV; sign gives opposite electron/hole tilt.
    V = _z_potential(height, barrier, field, sign, z)
    invface=2/(mass[:-1]+mass[1:]); a=_HBAR2_2M0*invface/dz**2
    diag=np.empty(n_pts-2); off=-a[1:-1]
    diag[:] = a[:-1]+a[1:]+V[1:-1]
    vals, vecs=eigh_tridiagonal(diag,off,select='i',select_range=(0,1))
    psi=np.zeros(n_pts); psi[1:-1]=vecs[:,0]; psi/=math.sqrt(np.trapezoid(psi*psi,z))
    cont=min(V[0],V[-1]); loc=float(np.trapezoid(psi[inside]**2,z[inside]))
    e1 = vals[1] if len(vals)>1 else float('nan')
    return vals[0], e1, cont, loc, z, psi

def _radial(remaining_offset, radius, md, mb):
    """Finite circular well (BenDaniel-Duke matching, barrier mass mb),
    reusing dot_levels' material-independent disk kernel [E]; `remaining_offset`
    is the barrier headroom left after the z-confinement energy is removed
    (adiabatic separable decoupling: the disk sees V - E_z, not the full
    offset). Returns (E0_eV, E1_eV or nan, bound, p_bound, rms_r_nm or None)."""
    if remaining_offset <= 0.:
        return float('nan'), float('nan'), False, False, None
    d = finite_disk_2d(remaining_offset * 1000., radius, md, mb)
    if not d.bound:
        return float('nan'), float('nan'), False, False, None
    e0 = d.E0_meV / 1000.
    e1 = d.E1_meV / 1000. if d.p_bound else float('nan')
    return e0, e1, True, d.p_bound, d.rms_r_nm

def _coulomb_binding_eV(l_e_xy, l_h_xy, z_sep_nm, eps_r):
    """Screened e-h Coulomb binding, frozen-orbital Gaussian envelopes [E]
    (same in-plane method as dot_levels._gauss_binding, Bernardini/BenDaniel
    geometry aside) with the field-driven z separation added in quadrature
    to the relative-coordinate width, so QCSE separation SUPPRESSES (never
    hard-zeroes) the binding as the carriers are pulled apart:
        E_b = e^2/(4 pi eps0 eps_r) * sqrt(pi) / sqrt(l_e^2 + l_h^2 + 2 z_sep^2)
    Evaluated from the normalized envelopes actually solved above, not
    fitted to any target wavelength."""
    L2 = l_e_xy*l_e_xy + l_h_xy*l_h_xy + 2.*z_sep_nm*z_sep_nm
    if L2 <= 0.: return float('nan')
    return math.sqrt(math.pi) * _E2_4PIEPS0_EV_NM / (eps_r * math.sqrt(L2))

def _qw_levels(s,d,m,de,Ve,Vh,F,ee,eh,ce,ch,le,lh,ze,pe,zh,ph,de_pad,dh_pad,n,pad,h_eff,r_eff,volume,geometry_label,
               ee1,eh1,d_ez,d_exy,d_hz,d_hxy,m_ez,m_exy,m_hz,m_hxy):
    """Same-composition QW fluctuation, adiabatic local-column model [A].
    `d_ez`/`d_exy`/`d_hz`/`d_hxy` and `m_ez`/`m_exy`/`m_hz`/`m_hxy` are the
    orientation-swapped growth-axis/lateral masses from `_growth_masses`,
    shared with the dot-column solve the caller already performed for
    `ee`/`eh`/`ce`/`ch` so both the local dot column and the surrounding
    well use the SAME (orientation-correct) masses."""
    er,_,cer,_,_,_=_z_state(s.wl_thickness_nm,Ve,d_ez,m_ez,F,-1,n,pad)
    hr,_,chr,_,_,_=_z_state(s.wl_thickness_nm,Vh,d_hz,m_hz,F,+1,n,pad)
    erp,_,_,_,_,_=_z_state(s.wl_thickness_nm,Ve,d_ez,m_ez,F,-1,n,2.*pad)
    hrp,_,_,_,_,_=_z_state(s.wl_thickness_nm,Vh,d_hz,m_hz,F,+1,n,2.*pad)
    bad=[]
    if not (ee < ce and eh < ch and er < cer and hr < chr): bad.append('dot-column or surrounding-QW vertical reservoir is unbound')
    if abs(erp-er)*1000. > TAIL_CONVERGENCE_TOL_MEV or abs(hrp-hr)*1000. > TAIL_CONVERGENCE_TOL_MEV:
        bad.append('surrounding-QW vertical reservoir padding-unconverged')
    de_lat=er-ee; dh_lat=hr-eh
    if de_lat <= 0: bad.append('electron QW lateral depth is nonpositive: no localized dot')
    if dh_lat <= 0: bad.append('hole QW lateral depth is nonpositive: no localized dot')
    re,rpe,oke,pbe,rmse=_radial(de_lat,r_eff,d_exy,d_exy)
    rh,rph,okh,pbh,rmsh=_radial(dh_lat,r_eff,d_hxy,d_hxy)
    ebe=ee+(re if oke else float('nan')); hbe=eh+(rh if okh else float('nan'))
    eth=min(er,ce); hth=min(hr,ch)
    ede=(eth-ebe)*1000. if oke else float('nan'); hde=(hth-hbe)*1000. if okh else float('nan')
    if not (oke and ede>0 and de_pad<=TAIL_CONVERGENCE_TOL_MEV): bad.append('electron fluctuation state is unbound, padding-unconverged, or laterally exhausted')
    if not (okh and hde>0 and dh_pad<=TAIL_CONVERGENCE_TOL_MEV): bad.append('hole fluctuation state is unbound, padding-unconverged, or laterally exhausted')
    if bad: return _invalid(s,bad,F)
    ov=float(np.trapezoid(pe*ph,ze)**2); ov=max(0.,min(1.,ov))
    zsep=abs(float(np.trapezoid(ze*pe*pe,ze))-float(np.trapezoid(zh*ph*ph,zh)))
    lee=rmse if rmse else r_eff; lhh=rmsh if rmsh else r_eff
    coul=_coulomb_binding_eV(lee,lhh,zsep,d.eps_r)
    ex=(de['Ec_eV']-de['Ev_eV'])+ebe+hbe-coul
    # Validity gate (Opus fix-round finding 1): a nonphysical (nonpositive
    # or non-finite) transition energy must never be reported valid=True.
    if not (math.isfinite(ex) and ex>0.):
        return _invalid(s,['nonphysical E_X'],F)
    # First excited state = min(z, radial) for the DOT COLUMN, same rule as
    # the isolated branch (Opus fix-round finding: ee1/eh1 were computed by
    # the caller but discarded here, hardwiring zg_e=zg_h=inf). Admission is
    # against THIS branch's own threshold eth/hth=min(er,ce)/min(hr,ch), not
    # the isolated-dot continuum ce/ch (hardening-round finding: ee1<ce
    # admitted a z-excited state above the branch's own selected continuum,
    # e.g. sp_split_e_meV=88.69 meV reported for a state 33 meV above the
    # surrounding-well edge at h=7, R=5, w=3.5, screening=1).
    zg_e=(ee1-ee) if (math.isfinite(ee1) and ee1<eth) else float('inf'); rg_e=(rpe-re) if (oke and math.isfinite(rpe) and ee+rpe<eth) else float('inf')
    zg_h=(eh1-eh) if (math.isfinite(eh1) and eh1<hth) else float('inf'); rg_h=(rph-rh) if (okh and math.isfinite(rph) and eh+rph<hth) else float('inf')
    eleft=float(_z_potential(h_eff,Ve,F,-1,np.array([-h_eff/2.-1.]))[0]); eright=float(_z_potential(h_eff,Ve,F,-1,np.array([h_eff/2.+1.]))[0])
    hleft=float(_z_potential(h_eff,Vh,F,+1,np.array([-h_eff/2.-1.]))[0]); hright=float(_z_potential(h_eff,Vh,F,+1,np.array([h_eff/2.+1.]))[0])
    # reservoir_kind/reservoir_energy_eV (Opus fix-round finding 4, refined by
    # the hardening round after commit 17a333f): report the channel each
    # carrier ACTUALLY escapes into (min(er,ce), min(hr,ch)) rather than
    # hardcoding the QW edge. 'ingan_qw' only when BOTH carriers select the
    # surrounding well, 'gan_barrier' only when BOTH select the wider GaN
    # plateau, and 'mixed' when the two carriers select DIFFERENT channels
    # (a spatially indirect situation -- electron confined by the well while
    # the hole escapes to the plateau, or vice versa -- that neither single
    # label describes exactly; the prior both-carriers rule mislabelled this
    # case 'gan_barrier'). reservoir_energy_eV tracks the selected (eth,hth)
    # pair in every case, which reduces exactly to the QW edge (gap+er+hr)
    # when both carriers select the QW and to gap+ce+ch when both select the
    # barrier; device.py only reads reservoir_kind as a label and
    # reservoir_energy_eV as a number (see module docstring) and never
    # branches on the kind string, so 'mixed' rows evaluate normally.
    # reservoir_electron_edge_meV/reservoir_hole_edge_meV report the SAME
    # selected per-carrier edge (eth, hth), not the raw surrounding-well
    # edge (er, hr) previously emitted here regardless of which channel was
    # actually selected.
    e_is_qw = er <= ce; h_is_qw = hr <= ch
    if e_is_qw and h_is_qw: reservoir_kind = 'ingan_qw'
    elif e_is_qw or h_is_qw: reservoir_kind = 'mixed'
    else: reservoir_kind = 'gan_barrier'
    reservoir_energy_eV = (de['Ec_eV']-de['Ev_eV'])+eth+hth
    # The optical acceptance sees the QW e-h continuum, irrespective of
    # which lower per-carrier escape plateaux control retention. [A]
    optical_reservoir_energy_eV = (de['Ec_eV']-de['Ev_eV'])+er+hr
    escape_e_mass = d_exy if e_is_qw else m_exy
    escape_h_mass = d_hxy if h_is_qw else m_hxy
    return NitrideLevels(ex,_HC_EV_NM/ex if ex>0 else float('nan'),True,True,ov,F,ebe*1000.,hbe*1000.,ede,hde,None,
        min(zg_e,rg_e)*1000. if math.isfinite(min(zg_e,rg_e)) else float('nan'),min(zg_h,rg_h)*1000. if math.isfinite(min(zg_h,rg_h)) else float('nan'),
        d_exy,d_hxy,True,(),'[V] BenDaniel & Duke, PR 152, 683 (1966); [A] same-composition adiabatic local-column QW fluctuation, both columns centered on the same z=0 midplane; lateral interface electrostatics neglected',
        le,lh,de_pad,dh_pad,eleft,eright,hleft,hright,reservoir_kind,reservoir_energy_eV,eth*1000.,hth*1000.,h_eff,r_eff,volume,
        geometry_label+'; [A] local-column thickness fluctuation with effective field length',
        'numerical padding/discretization reported separately; lateral interface electrostatics, shape and nonpolar strain/valence systematic error UNQUANTIFIED; '
        'symmetric-fluctuation placement assumption: in the unscreened tilted limit the reported lateral depth collapses to the mass-independent geometric '
        'estimate |F|*(H-w)/2 rather than a quantum-confinement difference (see module docstring); ' + _GRID_NOTE,
        orientation_factor(s.orientation,s.polarization_factor),s.height_nm,s.radius_nm,s.height_nm/(2.*s.radius_nm),
        optical_reservoir_energy_eV=optical_reservoir_energy_eV, optical_reservoir_kind='ingan_qw',
        escape_e_matrix_xy=escape_e_mass, escape_h_matrix_xy=escape_h_mass)

@lru_cache(maxsize=256)
def _levels_cached(s,T_K,n,pad):
    bad=_validate(s)
    if not math.isfinite(T_K) or T_K<=0: bad.append('T_K must be positive and finite')
    if bad: return _invalid(s,bad)
    d,m=ingaN(s.x_in),binary('GaN')
    d_ez,d_exy,d_hz,d_hxy=_growth_masses(d,s.orientation)
    m_ez,m_exy,m_hz,m_hxy=_growth_masses(m,s.orientation)
    h_eff,r_eff,volume,geometry_label=_geometry(s)
    de=band_edges(d,T_K,substrate=m,strain_fraction=s.strain_fraction,vbo_InN_GaN_eV=s.vbo_InN_GaN_eV,strain_c_fraction=s.strain_c_fraction)
    be=band_edges(m,T_K,substrate=m)
    Ve=be['Ec_eV']-de['Ec_eV']; Vh=de['Ev_eV']-be['Ev_eV']
    F=polarization_field(d,m,T_K,strain_fraction=s.strain_fraction,screening_fraction=s.screening_fraction,external_field_kVcm=s.external_field_kVcm,orientation=s.orientation,polarization_factor=s.polarization_factor)
    if Ve<=0: bad.append('electron offset is nonpositive')
    if Vh<=0: bad.append('hole offset is nonpositive')
    if bad: return _invalid(s,bad,F)
    ee,ee1,ce,le,ze,pe=_z_state(h_eff,Ve,d_ez,m_ez,F,-1,n,pad)
    eh,eh1,ch,lh,zh,ph=_z_state(h_eff,Vh,d_hz,m_hz,F,+1,n,pad)
    # A true discrete state has an exponentially decaying exterior tail, so
    # its eigenenergy is insensitive to doubling an already-large padding.
    # This replaces the arbitrary in-dot-probability > 0.5 validity gate [E].
    ee_pad,_,_,_,_,_=_z_state(h_eff,Ve,d_ez,m_ez,F,-1,n,2.*pad)
    eh_pad,_,_,_,_,_=_z_state(h_eff,Vh,d_hz,m_hz,F,+1,n,2.*pad)
    de_pad=abs(ee_pad-ee)*1000.
    dh_pad=abs(eh_pad-eh)*1000.
    if s.geometry_type == 'qw_fluctuation':
        return _qw_levels(s,d,m,de,Ve,Vh,F,ee,eh,ce,ch,le,lh,ze,pe,zh,ph,de_pad,dh_pad,n,pad,h_eff,r_eff,volume,geometry_label,
                           ee1,eh1,d_ez,d_exy,d_hz,d_hxy,m_ez,m_exy,m_hz,m_hxy)
    re,rpe,ok_e,pb_e,rms_e=_radial(ce-ee,r_eff,d_exy,m_exy)
    rh,rph,ok_h,pb_h,rms_h=_radial(ch-eh,r_eff,d_hxy,m_hxy)
    eb=ee+(re if ok_e else float('nan')); hb=eh+(rh if ok_h else float('nan'))
    if not (ok_e and eb<ce and de_pad<=TAIL_CONVERGENCE_TOL_MEV):
        bad.append('electron unbound, padding-unconverged, or laterally exhausted offset')
    if not (ok_h and hb<ch and dh_pad<=TAIL_CONVERGENCE_TOL_MEV):
        bad.append('hole unbound, padding-unconverged, or laterally exhausted offset')
    if bad: return _invalid(s,bad,F)
    # Envelopes are separable; normalized radial ground envelopes cancel in
    # overlap ratio under the common-disk approximation [E].
    ov=float(np.trapezoid(pe*ph,ze)**2); ov=max(0.,min(1.,ov))
    eps=(d.eps_r+m.eps_r)/2
    z_sep=abs(float(np.trapezoid(ze*pe*pe,ze))-float(np.trapezoid(zh*ph*ph,zh)))
    l_e_xy = rms_e if rms_e else s.radius_nm; l_h_xy = rms_h if rms_h else s.radius_nm
    coul=_coulomb_binding_eV(l_e_xy,l_h_xy,z_sep,eps) # eV, screened Gaussian-envelope proxy [E]
    ex=(de['Ec_eV']-de['Ev_eV'])+eb+hb-coul
    # Validity gate (Opus fix-round finding 1): a nonphysical (nonpositive
    # or non-finite) transition energy must never be reported valid=True.
    if not (math.isfinite(ex) and ex>0.):
        return _invalid(s,['nonphysical E_X'],F)
    # First excited state = min(z-excitation, radial p-shell excitation),
    # each admitted only if it is itself bound below the local continuum.
    zg_e = (ee1-ee) if (math.isfinite(ee1) and ee1<ce) else float('inf')
    rg_e = (rpe-re) if (ok_e and math.isfinite(rpe)) else float('inf')
    sp_e = min(zg_e,rg_e)
    zg_h = (eh1-eh) if (math.isfinite(eh1) and eh1<ch) else float('inf')
    rg_h = (rph-rh) if (ok_h and math.isfinite(rph)) else float('inf')
    sp_h = min(zg_h,rg_h)
    valid=True
    half=h_eff/2.
    eleft=float(_z_potential(h_eff,Ve,F,-1,np.array([-half-1.]))[0])
    eright=float(_z_potential(h_eff,Ve,F,-1,np.array([half+1.]))[0])
    hleft=float(_z_potential(h_eff,Vh,F,+1,np.array([-half-1.]))[0])
    hright=float(_z_potential(h_eff,Vh,F,+1,np.array([half+1.]))[0])
    return NitrideLevels(ex,_HC_EV_NM/ex if ex>0 else float('nan'),True,True,ov,F,
                          eb*1000,hb*1000,(ce-eb)*1000,(ch-hb)*1000,None,
                          sp_e*1000 if math.isfinite(sp_e) else float('nan'),
                          sp_h*1000 if math.isfinite(sp_h) else float('nan'),
                          m_exy,m_hxy,valid,tuple(bad),
                          '[V] Rinke et al., PRB 77, 075202 (2008), Table V; '
                          '[V] Bernardini et al., PRB 56, R10024 (1997); '
                          '[V/DR] Bernardini & Fiorentini, pss(b) 216, 391 (1999), '
                          'Eqs. 7-8, continuous isolated-slab field geometry; '
                          '[V] BenDaniel & Duke, PR 152, 683 (1966); '
                          '[E] separable disk, finite-barrier radial confinement, and screened '
                          'Gaussian-envelope Coulomb approximation; '
                          '[A] partial screening parameter and real-energy retention',
                          le,lh,de_pad,dh_pad,eleft,eright,hleft,hright,
                          # reservoir_energy_eV (hardening round after commit 17a333f): the bulk
                          # GaN Varshni edge at this T_K minus the SAME 25 meV [A] localization/
                          # Urbach offset fsim_core.device._nitride_reservoir_energy_eV applies on
                          # its own zero-wetting-layer branch (both read the identical band_edges(
                          # GaN, T_K, substrate=GaN) call: `be` above), not the dot's OWN strained
                          # InGaN gap (de['Ec_eV']-de['Ev_eV']) previously reported here, which was
                          # numerically BELOW this row's own E_X and disagreed with the GaN-barrier
                          # reservoir_kind already reported alongside it. device.py does not read
                          # this field for isolated-dot rows (it keeps its own bulk-edge helper's
                          # value); this field exists for documentation/consistency, see module
                          # docstring "RESERVOIR KIND AND ENERGY".
                          'gan_barrier',(be['Ec_eV']-be['Ev_eV'])-.025,0.,0.,h_eff,r_eff,volume,
                          geometry_label+'; [A] effective field length used consistently',
                          'numerical padding/discretization reported separately; shape and nonpolar strain/valence systematic error UNQUANTIFIED; effective-height spread is sensitivity only; ' + _GRID_NOTE,
                          orientation_factor(s.orientation,s.polarization_factor),s.height_nm,s.radius_nm,s.height_nm/(2.*s.radius_nm),
                          optical_reservoir_energy_eV=(be['Ec_eV']-be['Ev_eV'])-.025,
                          optical_reservoir_kind='gan_barrier', escape_e_matrix_xy=m_exy,
                          escape_h_matrix_xy=m_hxy)

def _invalid(s,reasons,F=0.):
    return NitrideLevels(float('nan'),float('nan'),False,False,0.,F,float('nan'),float('nan'),float('nan'),float('nan'),None,float('nan'),float('nan'),float('nan'),float('nan'),False,tuple(reasons),'[A] invalid geometry/offset rejected before model evaluation')

def levels(system, T_K=300.0, *, z_points=1201, exterior_nm=45.0):
    """Return immutable cached levels. Numerical controls are cache keys."""
    if not isinstance(system,NitrideDotSystem): raise TypeError('system must be NitrideDotSystem')
    if system.geometry_type == 'isolated_dot' and math.isfinite(system.wl_thickness_nm) and system.wl_thickness_nm>0:
        raise ValueError('nonzero wl_thickness_nm is unsupported: no wetting-layer continuum is implemented')
    return _levels_cached(system,float(T_K),int(z_points),float(exterior_nm))

def rates(lv,T_K,*,tau_rad0_ns=TAU_RAD0_DEFAULT_NS,n_dot_cm2=1e10,tau_cap_ps=10.,tau_cap_scales_with_density=False,channel='min',k_nr_ns=0.):
    """Absolute detailed-balance escape rates in 1/ns; no cavity/Purcell input.

    tau_rad0_ns default is TAU_RAD0_DEFAULT_NS [E], see module docstring
    constants; callers (fsim_core.device) always pass an explicit value.
    The N2D prefactor is a bulk-channel DOS approximation [A]: for QW rows
    it uses the in-plane mass of each carrier's selected escape channel
    (surrounding InGaN well or GaN plateau), not necessarily the dot mass.
    Pair channels remain unsupported [A]."""
    bad=list(lv.invalid_reasons)
    if not lv.valid: bad.append('levels are invalid')
    # Opus fix-round finding 3: an overlap below the numerical noise floor
    # (e.g. 1e-21) is not a resolved lifetime prediction, just underflow in
    # the separable-envelope integral; flag it instead of reporting an
    # astronomical (numerically meaningless) radiative lifetime.
    elif lv.overlap_sq < 1e-12: bad.append('overlap_unresolved')
    if T_K<=0 or tau_rad0_ns<=0 or n_dot_cm2<=0 or tau_cap_ps<=0 or k_nr_ns<0: bad.append('invalid rate input')
    if channel not in ('min','electron','hole','pair','pair_half'): bad.append('unknown channel')
    if channel in ('pair','pair_half'): bad.append('pair channel requires a bound wetting-layer continuum')
    if bad: return dict(gamma_X0_ns=float('nan'),gamma_XX0_ns=float('nan'),k_X_ns=float('nan'),k_XX_ns=float('nan'),escape_prefactor_ns=float('nan'),E_a_meV=float('nan'),S0=float('nan'),tau_cap_ps_used=float('nan'),valid=False,invalid_reasons=tuple(bad),provenance='[A] invalid rate request')
    tc=tau_cap_ps*(1e10/n_dot_cm2) if tau_cap_scales_with_density else tau_cap_ps
    e_mass = lv.escape_e_matrix_xy if math.isfinite(lv.escape_e_matrix_xy) else lv.m_e_matrix_xy
    h_mass = lv.escape_h_matrix_xy if math.isfinite(lv.escape_h_matrix_xy) else lv.m_h_matrix_xy
    choices={'electron':(lv.dE_e_meV,e_mass),'hole':(lv.dE_h_meV,h_mass)}
    name=min(choices,key=lambda q:choices[q][0]) if channel=='min' else channel
    ea,mass=choices[name]
    if ea<0: return dict(gamma_X0_ns=float('nan'),gamma_XX0_ns=float('nan'),k_X_ns=float('nan'),k_XX_ns=float('nan'),escape_prefactor_ns=float('nan'),E_a_meV=ea,S0=float('nan'),tau_cap_ps_used=tc,valid=False,invalid_reasons=('negative escape barrier',),provenance='[A] negative barrier is not floored')
    n2d=mass*_M0*_KB_SI*T_K/(math.pi*_HBAR_SI**2)/1e4
    attempt=(1000/tc)*(n2d/n_dot_cm2)
    k=attempt*math.exp(-(ea/1000)/(KB_EV*T_K))+k_nr_ns
    gx=lv.overlap_sq/tau_rad0_ns; gxx=2*gx
    return dict(gamma_X0_ns=gx,gamma_XX0_ns=gxx,k_X_ns=k,k_XX_ns=2*k,escape_prefactor_ns=attempt,E_a_meV=ea,S0=gx/(gx+k),tau_cap_ps_used=tc,valid=True,invalid_reasons=(),provenance='[A] 2D-reservoir detailed balance N2D=m kT/(pi hbar2); [A] XX=2X cascade transfer, Reischle et al., OE 16, 12771 (2008); selected '+name+' channel')
