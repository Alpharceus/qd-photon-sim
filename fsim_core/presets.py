"""Preset registry for the device designer (phase D1): named fragments for the
dot, thermal template, cavity, and drive blocks, plus a composer that builds a
full DeviceDesign from four preset keys.

Headless by construction (three-layer rule): this module returns/mutates
DeviceDesign objects only -- no GUI, no plotting. The designer combo boxes are
a thin client that pick keys out of the dicts below.

Every preset here is a point choice out of a wider, honestly-tagged band (see
design_meta.META) -- picking a preset does not upgrade its provenance tag.
"""
from __future__ import annotations

from .device import CavityBlock, DeviceDesign, DotBlock, DriveBlock

DOT_PRESETS = {
    "chatzarakis-class": {
        "delta_xx": 5.9, "gamma_scale": 1.0, "r_xx": 0.72,
        "note": "measured (211)B InAs/GaAs cavity dot [DR]",
    },
    "staged-inp-gaasp": {
        "delta_xx": 3.5,
        "note": "[A] midpoint, unmeasured platform",
    },
    "piezo-variant": {
        "delta_xx": 8.0,
        "note": "[A] selection tail of the (211)B 4-13 meV distribution",
    },
}

TEMPLATE_PRESETS = {
    "GaAs": {
        "layers": [
            {"name": "epi/cladding", "t_um": 1.5, "k300": 15.0, "alpha": 0.5, "spread": False},
        ],
        "substrate": {"name": "GaAs", "k300": 55.0, "alpha": 1.25},
    },
    "GaAs-Si": {
        "layers": [
            {"name": "epi/cladding", "t_um": 1.5, "k300": 15.0, "alpha": 0.5, "spread": False},
            {"name": "GaAs buffer", "t_um": 2.0, "k300": 44.0, "alpha": 1.25, "spread": True},
        ],
        "substrate": {"name": "Si", "k300": 148.0, "alpha": 1.35},
    },
}

CAVITY_PRESETS = {
    "none": {"enabled": False},
    "planar-lambda": {"enabled": True, "type": "planar", "kappa": 2.0, "G": 8.0, "F_P": 5.0},
    "micropillar": {"enabled": True, "type": "planar", "kappa": 0.5, "G": 12.0, "F_P": 15.0},
    "sin-waveguide": {"enabled": True, "type": "sin_waveguide", "beta_sin": 0.3},
}

DRIVE_PRESETS = {
    "cw-electrical": {"V": 1.9, "I_uA": 10.0, "duty": 1.0, "mu": 0.5, "b_e": 0.02},
    "pulsed-electrical": {"V": 1.9, "I_uA": 50.0, "duty": 0.05, "mu": 0.5, "b_e": 0.02},
    "optical": {"V": 0.0, "I_uA": 0.0, "duty": 0.0, "mu": 0.33, "b_e": 0.0},  # no injection, no junction heating
}


def apply_dot_preset(design: DeviceDesign, key: str) -> DeviceDesign:
    p = DOT_PRESETS[key]
    b = DotBlock()
    design.dot.delta_xx = p.get("delta_xx", b.delta_xx)
    design.dot.gamma_scale = p.get("gamma_scale", b.gamma_scale)
    design.dot.r_xx = p.get("r_xx", b.r_xx)
    return design


def apply_template_preset(design: DeviceDesign, key: str) -> DeviceDesign:
    p = TEMPLATE_PRESETS[key]
    design.thermal.layers = [dict(L) for L in p["layers"]]
    design.thermal.substrate = dict(p["substrate"])
    return design


def apply_cavity_preset(design: DeviceDesign, key: str) -> DeviceDesign:
    p = CAVITY_PRESETS[key]
    b = CavityBlock()
    design.cavity.enabled = bool(p["enabled"])
    design.cavity.type = p.get("type", b.type)
    design.cavity.kappa = p.get("kappa", b.kappa)
    design.cavity.G = p.get("G", b.G)
    design.cavity.F_P = p.get("F_P", b.F_P)
    design.cavity.beta_sin = p.get("beta_sin", b.beta_sin)
    return design


def apply_drive_preset(design: DeviceDesign, key: str) -> DeviceDesign:
    p = DRIVE_PRESETS[key]
    b = DriveBlock()
    design.drive.V = p.get("V", b.V)
    design.drive.I_uA = p.get("I_uA", b.I_uA)
    design.drive.duty = p.get("duty", b.duty)
    design.drive.mu = p.get("mu", b.mu)
    design.drive.b_e = p.get("b_e", b.b_e)
    return design


def preset_device(dot: str, template: str, cavity: str, drive: str,
                  name: str = "preset-device") -> DeviceDesign:
    """Compose one DeviceDesign from four preset keys (dot/template/cavity/drive)."""
    d = DeviceDesign(name=name)
    apply_dot_preset(d, dot)
    apply_template_preset(d, template)
    apply_cavity_preset(d, cavity)
    apply_drive_preset(d, drive)
    return d
