"""Helper functions for presentation2/sections/03_driving.json's
`repo_numbers` `how` expressions (SCHEMA.md form (b)).

Every `how` expression must stay a single Python expression; these thin
wrappers exist only so that expression can call the real fsim_core API
(fsim_core.transport, fsim_core.device) with the literal, documented
overrides a slide describes, instead of repeating a multi-line DeviceDesign
construction inline in the JSON. No independent physics lives here -- each
function does exactly what its slide's file/how note says: load a card
through fsim_core.device.DeviceDesign.load, apply the stated overrides, and
call fsim_core.device.evaluate (or, for the two transport-only numbers,
construct a bare fsim_core.transport.Diode and call its own methods).

Not a figure script itself (writes no PNG); imported by the JSON's `how`
expressions as `presentation2.figures.helpers_03` with the repository root
on sys.path (the validator/test-command convention), and importable the
same way interactively.
"""
from __future__ import annotations

import math
import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[2]
if str(_ROOT) not in sys.path:
    sys.path.insert(0, str(_ROOT))

from fsim_core import materials, transport  # noqa: E402
from fsim_core.device import DeviceDesign, evaluate  # noqa: E402

# The favourable diagnostic corner (out/rt_edge/verdict.md's best pulsed-g2
# row for the gainp card): dot.delta_xx=8 meV, dot.gamma300=6 meV,
# emission.NA=0.8, emission.R_back=0.95, emission.L_um=250, pulsed drive at
# 80 MHz / 100 ps -- see presentation2/figures/fig_03_22_efficiency_chain.py.
FAVORABLE_OVERRIDES = {
    "dot.delta_xx": 8.0,
    "dot.gamma300": 6.0,
    "emission.NA": 0.8,
    "emission.R_back": 0.95,
    "emission.L_um": 250.0,
}
FAVORABLE_PULSE_WIDTH_NS = 0.1
FAVORABLE_REP_RATE_HZ = 80.0e6
# CW raw diagnostic min (verdict.md) is the sweep's low-irf_ps grid endpoint
# (scripts/run_rt_edge.py RANGE_BOUNDS["irf_ps"] = (50.0, 200.0)); g2_cw0
# itself does not depend on irf_ps.
FAVORABLE_CW_IRF_PS = 50.0


def _default_diode() -> transport.Diode:
    """A transport.Diode with every field at its dataclass default. f_qd and
    area_um2 do not depend on the material choice, so any valid Material
    fills the five required positional layers."""
    gaas = materials.binary("GaAs")
    return transport.Diode(gaas, gaas, gaas, gaas, gaas)


def f_qd(n_dot_cm2: float) -> float:
    """fsim_core.transport.Diode's default f_qd(n_dot_cm2)."""
    return _default_diode().f_qd(n_dot_cm2)


def mesa_over_aperture_area_ratio(aperture_diameter_um: float = 0.4) -> float:
    """Diode's default mesa area_um2 divided by the aperture area implied by
    aperture_diameter_um (0.4 um: both edge-emitter cards' aperture.diameter_um)."""
    return _default_diode().area_um2 / (math.pi * (aperture_diameter_um / 2.0) ** 2)


def _apply_overrides(design: DeviceDesign, overrides: dict) -> None:
    for dotted, value in overrides.items():
        obj = design
        *parents, leaf = dotted.split(".")
        for name in parents:
            obj = getattr(obj, name)
        setattr(obj, leaf, value)


def load_design(card_path: str, favorable: bool = False, cw: bool = False,
                 T_hs: float | None = None, irf_ps: float | None = None) -> DeviceDesign:
    """Load `card_path` (relative to the repository root) fresh and apply the
    overrides a slide describes: `favorable` for the favourable diagnostic
    corner (pulsed by default unless `cw`), `cw` to force CW drive (with
    `irf_ps` overriding drive.cw_irf_fwhm_ps), `T_hs` to override
    thermal.T_hs."""
    design = DeviceDesign.load(_ROOT / card_path)
    if favorable:
        _apply_overrides(design, FAVORABLE_OVERRIDES)
        if not cw:
            design.drive.cw = False
            design.drive.duty = FAVORABLE_PULSE_WIDTH_NS * 1e-9 * FAVORABLE_REP_RATE_HZ
            design.drive.diode["tau_pulse_ns"] = FAVORABLE_PULSE_WIDTH_NS
    if cw:
        design.drive.cw = True
        if irf_ps is not None:
            design.drive.cw_irf_fwhm_ps = float(irf_ps)
    if T_hs is not None:
        design.thermal.T_hs = float(T_hs)
    return design


def evaluate_card(card_path: str, key: str | None = None, favorable: bool = False,
                   cw: bool = False, T_hs: float | None = None,
                   irf_ps: float | None = None):
    """Fresh DeviceDesign.load(card_path) -> evaluate() at thermal.T_hs
    (after any override), returning scalars[key] if given, else the full
    scalars dict."""
    design = load_design(card_path, favorable=favorable, cw=cw, T_hs=T_hs, irf_ps=irf_ps)
    scalars = evaluate(design, T_grid=[design.thermal.T_hs])["scalars"]
    return scalars[key] if key is not None else scalars


def n_dots_expected(card_path: str) -> float:
    """aperture.density_cm2 * pi*(aperture.diameter_um/2)**2 * 1e-8 -- the
    Poisson-expected dot count under the card's own aperture, no evaluate()
    needed (aperture geometry, not a device-chain output)."""
    design = DeviceDesign.load(_ROOT / card_path)
    ap = design.aperture
    return ap.density_cm2 * math.pi * (ap.diameter_um / 2.0) ** 2 * 1e-8


def drive_field(card_path: str, field: str) -> float:
    """A literal DriveBlock field as the card actually stores it (e.g.
    'I_uA'), read fresh rather than copied from the card text."""
    design = DeviceDesign.load(_ROOT / card_path)
    return float(getattr(design.drive, field))


def gate_pass(card_path: str, threshold: float = 0.5, favorable: bool = False,
              T_hs: float | None = None) -> float:
    """1.0 if scalars['g2_op'] < threshold else 0.0, at the given corner --
    the section's pass/fail gate, as a float for the how-expression contract."""
    g2 = evaluate_card(card_path, "g2_op", favorable=favorable, T_hs=T_hs)
    return 1.0 if g2 < threshold else 0.0
