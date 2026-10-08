"""Checks for fsim_core.nitride_levels.z_profile (FSIM Studio band view wrapper).

z_profile is a read-only wrapper over the private _z_potential/_z_state solve
that levels() already performs.  These checks confirm it (a) reproduces the
private functions bit-for-bit, (b) reports the same overlap_sq as levels(),
(c) obeys the documented field geometry (tilt F*h inside, flat exterior;
e*(kV/cm)*nm = 1e-4 eV), and (d) on the shipped nitride-cavity-pulse-design
card reproduces the EVALUATED overlap_sq/field_kVcm at the operating point
within 1e-9.  The operating point is captured by wrapping nitride_levels.levels
during evaluate(), and separately rebuilt from the evaluated scalars
(applied_field_kVcm, T_j) the way fsim_studio/scene.py does.
"""
import sys
import dataclasses
from pathlib import Path
import numpy as np
sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import fsim_core.nitride_levels as nl
from fsim_core.nitride_materials import ingaN, binary, band_edges, polarization_field
from fsim_core.device import DeviceDesign, evaluate

checks = []
def ck(ok, name):
    checks.append(bool(ok))
    if not ok: print("FAIL", name)

ROOT = Path(__file__).resolve().parents[1]

def private_reference(s, T, n=1201, pad=45.):
    """Independent replay of the private call sequence in _levels_cached."""
    d, m = ingaN(s.x_in), binary('GaN')
    dez, _, dhz, _ = nl._growth_masses(d, s.orientation)
    mez, _, mhz, _ = nl._growth_masses(m, s.orientation)
    h, _, _, _ = nl._geometry(s)
    de = band_edges(d, T, substrate=m, strain_fraction=s.strain_fraction,
                    vbo_InN_GaN_eV=s.vbo_InN_GaN_eV, strain_c_fraction=s.strain_c_fraction)
    be = band_edges(m, T, substrate=m)
    Ve = be['Ec_eV'] - de['Ec_eV']; Vh = de['Ev_eV'] - be['Ev_eV']
    F = polarization_field(d, m, T, strain_fraction=s.strain_fraction,
                           screening_fraction=s.screening_fraction,
                           external_field_kVcm=s.external_field_kVcm,
                           orientation=s.orientation, polarization_factor=s.polarization_factor)
    e = nl._z_state(h, Ve, dez, mez, F, -1, n, pad)
    hh = nl._z_state(h, Vh, dhz, mhz, F, +1, n, pad)
    return dict(e=e, h=hh, F=F, h_eff=h, de=de,
                pot_e=nl._z_potential(h, Ve, F, -1, e[4]),
                pot_h=nl._z_potential(h, Vh, F, +1, hh[4]))

systems = {
    'default c-plane': (nl.NitrideDotSystem(), 300.),
    'card dot h3 r10 x.25': (nl.NitrideDotSystem(height_nm=3., radius_nm=10., x_in=.25), 300.),
    'screened 0.6, 250 K': (nl.NitrideDotSystem(height_nm=2., radius_nm=12.5, x_in=.4, screening_fraction=.6), 250.),
    'a-plane': (nl.NitrideDotSystem(orientation='a_plane', height_nm=3.), 300.),
    'lens shape': (nl.NitrideDotSystem(shape='lens', height_nm=4., radius_nm=12.), 273.),
    'external field +500': (nl.NitrideDotSystem(external_field_kVcm=500.), 300.),
}

for label, (s, T) in systems.items():
    zp = nl.z_profile(s, T)
    ref = private_reference(s, T)
    ze, pe = ref['e'][4], ref['e'][5]
    zh, ph = ref['h'][4], ref['h'][5]
    ck(np.array_equal(zp['z_nm'], ze) and np.array_equal(ze, zh), f'{label}: z grid identical to _z_state')
    ck(np.array_equal(zp['psi_e'], pe) and np.array_equal(zp['psi_h'], ph), f'{label}: psi_e/psi_h identical to _z_state')
    ck(np.array_equal(zp['cb_eV'], ref['de']['Ec_eV'] + ref['pot_e'])
       and np.array_equal(zp['vb_eV'], ref['de']['Ev_eV'] - ref['pot_h']), f'{label}: CB/VB identical to _z_potential')
    ck(zp['E_e_z_meV'] == ref['e'][0] * 1000. and zp['E_h_z_meV'] == ref['h'][0] * 1000., f'{label}: z levels identical')
    ck(zp['field_kVcm'] == ref['F'], f'{label}: field identical to polarization_field')
    lv = nl.levels(s, T)
    if lv.valid:
        ck(abs(zp['overlap_sq'] - lv.overlap_sq) <= 1e-15, f'{label}: overlap_sq equals levels().overlap_sq')
    for psi in (zp['psi_e'], zp['psi_h']):
        ck(abs(np.trapezoid(psi * psi, zp['z_nm']) - 1.) < 1e-12, f'{label}: psi normalized')
    # Field geometry [DR]: e*(kV/cm)*(nm) = 1e-4 eV; CB(+h/2)-CB(-h/2) = -F*h*1e-4.
    h = zp['effective_height_nm']; z = zp['z_nm']; cb = zp['cb_eV']; vb = zp['vb_eV']
    i_lo = int(np.argmin(np.abs(z + h / 2))); i_hi = int(np.argmin(np.abs(z - h / 2)))
    tilt = cb[i_hi] - cb[i_lo]
    ck(abs(tilt - (-zp['field_kVcm'] * 1e-4 * (z[i_hi] - z[i_lo]))) < 1e-12, f'{label}: CB tilt = -F*h*1e-4 eV')
    ck(abs((vb[i_hi] - vb[i_lo]) - tilt) < 1e-12, f'{label}: VB tilt parallel to CB')
    out_l = z < -h / 2 - 1e-9; out_r = z > h / 2 + 1e-9
    ck(np.ptp(cb[out_l]) < 1e-12 and np.ptp(cb[out_r]) < 1e-12, f'{label}: exterior flat (zero field outside)')
    ck(abs((cb[i_lo + 0] - zp['Ec_dot_eV']) - (-zp['field_kVcm'] * 1e-4 * z[i_lo])) < 1e-12, f'{label}: CB referenced to Ec_dot at z=0')

# Invalid input raises.
try:
    nl.z_profile(nl.NitrideDotSystem(height_nm=-1.), 300.); ck(False, 'negative height raises')
except ValueError:
    ck(True, 'negative height raises')
ck(nl.z_profile(nl.NitrideDotSystem())['T_K'] == 300.0, 'T=None means 300 K (levels default)')

# Card operating point: capture the (system, Tj) evaluate() hands to levels().
captured = []
_orig = nl.levels
def _spy(system, T_K=300.0, **kw):
    captured.append((system, float(T_K)))
    return _orig(system, T_K, **kw)
nl.levels = _spy
try:
    d = DeviceDesign.load(ROOT / 'cards' / 'nitride-cavity-pulse-design.yaml')
    res = evaluate(d, T_grid=[300.0])
finally:
    nl.levels = _orig
sc = res['scalars']
op = [c for c in captured if abs(c[1] - float(sc['T_j'])) < 1e-12
      and abs(c[0].external_field_kVcm - float(sc['applied_field_kVcm'])) < 1e-12]
ck(len(op) >= 1, 'operating-point levels() call captured')
if op:
    s_op, T_op = op[-1]
    zp = nl.z_profile(s_op, T_op)
    ck(abs(zp['overlap_sq'] - float(sc['overlap_sq'])) < 1e-9,
       f"card overlap_sq: z_profile {zp['overlap_sq']:.12g} vs evaluated {float(sc['overlap_sq']):.12g} (1e-9)")
    ck(abs(zp['field_kVcm'] - float(sc['field_kVcm'])) < 1e-9, 'card field_kVcm matches evaluated (1e-9)')
    # scene.py reconstruction: card dot block + applied_field_kVcm at T_j.
    dot_kw = dict(d.nitride['dot'])
    s_rb = nl.NitrideDotSystem(**{**dot_kw, 'external_field_kVcm': float(sc['applied_field_kVcm'])})
    ck(dataclasses.astuple(s_rb) == dataclasses.astuple(s_op), 'scene reconstruction equals captured operating system')
    ck(abs(nl.z_profile(s_rb, float(sc['T_j']))['overlap_sq'] - float(sc['overlap_sq'])) < 1e-9,
       'scene reconstruction overlap_sq matches evaluated (1e-9)')
    # 03-3d.md design note: about 1.17 eV of band tilt across the 3 nm dot.
    ck(abs(abs(zp['field_kVcm']) * 1e-4 * zp['effective_height_nm'] - 1.17) < 0.01, 'card tilt across dot ~1.17 eV')

print(f'{sum(checks)}/{len(checks)} nitride z_profile checks passed')
sys.exit(0 if all(checks) else 1)
