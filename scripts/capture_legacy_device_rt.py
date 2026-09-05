"""Capture the pre-RT device regression fixture.

The key list is intentionally the legacy public result surface.  Keeping it
fixed means this script remains a meaningful replay check after opt-in RT
fields are added.
"""
from __future__ import annotations

from dataclasses import asdict
from pathlib import Path
import sys

import numpy as np


ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from fsim_core.device import DeviceDesign, evaluate
OUT = ROOT / "verify" / "data" / "legacy_device_rt.npz"
T_GRID = np.array([4.0, 77.0, 120.0, 220.0, 300.0, 350.0])
CURVES = ("T_hs", "Tj", "g2", "eps", "rho2", "gamma")
SCALARS = ("T_j_op", "dT_J", "runaway", "gamma_op", "eps_op", "rho_op",
           "g2_op", "t_x_op", "brightness_per_pulse", "T_c", "F_eff", "N_w",
           "aperture_g2_penalty", "tag_chain")


def designs():
    staged = DeviceDesign.load(ROOT / "cards" / "staged-device-design.yaml")
    cavity = DeviceDesign(name="legacy-cavity")
    cavity.cavity.enabled = True
    held = DeviceDesign(name="legacy-held")
    held.cavity.enabled = True
    held.filter.track = "hold"
    ibm = DeviceDesign(name="legacy-ibm")
    ibm.dot.lineshape = "ibm"
    pl = DeviceDesign(name="legacy-pl")
    pl.drive.mode = "PL"
    el = DeviceDesign(name="legacy-el")
    el.drive.mode = "EL"
    return [staged, cavity, held, ibm, pl, el]


def main():
    payload = {"T_grid": T_GRID, "names": np.array([d.name for d in designs()])}
    for i, design in enumerate(designs()):
        result = evaluate(design, T_GRID)
        for name in CURVES:
            payload[f"{i}/curve/{name}"] = result["curves"][name]
        for name in SCALARS:
            payload[f"{i}/scalar/{name}"] = np.array(result["scalars"][name])
        # Keep the original dataclass surface, so regeneration is byte-for-byte
        # stable after optional RT fields have been appended.
        legacy = asdict(design)
        legacy.pop("provenance", None)
        for block, fields in {
            "dot": ("linewidth", "gamma0", "a_ac", "E_LO", "gamma300"),
            "ret": ("mode", "preset", "system", "tau_rad_ns", "channel"),
            "drive": ("diode", "n_dot_cm2"),
        }.items():
            for name in fields:
                legacy[block].pop(name, None)
        payload[f"{i}/design"] = np.array(repr(legacy))
    np.savez(OUT, **payload)


if __name__ == "__main__":
    main()
