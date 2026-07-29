"""Provenance metadata for every editable scalar field of DeviceDesign.

One entry per dotted field path ("block.field") holding {unit, tag, source,
lo, hi}. lo/hi are the widest *known-physical* band for the quantity -- they
drive GUI warnings only, never clamping (an out-of-band point is still a
legal, honestly-tagged conditional number).

Headless by construction (three-layer rule): no GUI imports here. The
designer looks up META by path to draw its tag bullets/tooltips/warnings.
"""
from __future__ import annotations

import dataclasses

from .device import DeviceDesign

# tag in {"V", "DR", "E", "A"} -- see fsim_core.card.Tag for the ordering rule.
META = {
    # ---- dot
    "dot.delta_xx": {"unit": "meV", "tag": "A", "lo": 2.0, "hi": 14.0,
                      "source": "(211)B InAs/GaAs 4-13 meV distribution; platform-dependent"},
    "dot.gamma_scale": {"unit": "1", "tag": "A", "lo": 0.1, "hi": 5.0,
                         "source": "unitless scale on the class Gamma(T) proxy"},
    "dot.r_xx": {"unit": "1", "tag": "DR", "lo": 0.3, "hi": 1.5,
                 "source": "Gamma_XX/Gamma_X, class-proxy V-a joint fit"},

    # ---- drive
    "drive.V": {"unit": "V", "tag": "E", "lo": 0.0, "hi": 5.0,
                "source": "typical QD diode forward bias"},
    "drive.I_uA": {"unit": "uA", "tag": "E", "lo": 0.0, "hi": 500.0,
                   "source": "typical injection current range"},
    "drive.duty": {"unit": "1", "tag": "A", "lo": 0.0, "hi": 1.0,
                   "source": "duty cycle, device-dependent"},
    "drive.mu": {"unit": "1", "tag": "A", "lo": 0.0, "hi": 5.0,
                 "source": "cap-2 mean loading per pulse"},
    "drive.b_e": {"unit": "1", "tag": "A", "lo": 0.0, "hi": 1.0,
                  "source": "injection background amplitude (signal units)"},
    "drive.b_e_m": {"unit": "1", "tag": "A", "lo": 0.5, "hi": 4.0,
                    "source": "injection background current exponent"},
    "drive.b_e_Eact": {"unit": "meV", "tag": "A", "lo": 0.0, "hi": 300.0,
                        "source": "WL/barrier EL activation energy"},
    "drive.I_ref_uA": {"unit": "uA", "tag": "A", "lo": 0.1, "hi": 1000.0,
                        "source": "reference current for the injection-background law"},
    "drive.mode": {"unit": "enum", "tag": "A", "lo": 0.0, "hi": 1.0,
                   "source": "EL (electrical injection) | PL (optical excitation, "
                             "disables dg_inj and the injection background)"},
    "drive.dg_inj": {"unit": "meV", "tag": "A", "lo": 0.0, "hi": 20.0,
                      "source": "Kitamura EL~4xPL anchor; unmeasured on our platform"},
    "drive.p_inj": {"unit": "1", "tag": "A", "lo": 0.5, "hi": 3.0,
                     "source": "F5' injection broadening current exponent"},
    "drive.F_p": {"unit": "1", "tag": "A", "lo": 0.0, "hi": 2.0,
                  "source": "pump Fano factor; 1=Poisson, <1 quiet pump (Zhao RT), "
                            "0=turnstile"},
    "drive.eta_capture": {"unit": "1", "tag": "A", "lo": 0.0, "hi": 1.0,
                           "source": "dot capture fraction; F8b thinning"},
    "drive.C_dep_pF": {"unit": "pF", "tag": "E", "lo": 0.1, "hi": 100.0,
                        "source": "Zhao: down to ~3.5 pF"},

    # ---- retention (class-proxy Arrhenius overrides; 0 -> use class proxy)
    "ret.a_esc": {"unit": "1", "tag": "A", "lo": 0.0, "hi": 1e7,
                  "source": "class-proxy Arrhenius prefactor override"},
    "ret.E_a": {"unit": "meV", "tag": "A", "lo": 0.0, "hi": 400.0,
                "source": "class-proxy Arrhenius activation override"},
    "ret.b_p": {"unit": "1", "tag": "A", "lo": 0.0, "hi": 1e5,
                "source": "class-proxy phonon-escape prefactor override"},
    "ret.E_b": {"unit": "meV", "tag": "A", "lo": 0.0, "hi": 100.0,
                "source": "class-proxy phonon-escape activation override"},
    "ret.b0": {"unit": "1", "tag": "A", "lo": 0.0, "hi": 1.0,
               "source": "T-independent background floor override"},
    "ret.beta": {"unit": "1", "tag": "A", "lo": 0.0, "hi": 20.0,
                 "source": "escape-coupled background weight override"},

    # ---- thermal (layers/substrate excluded -- own stack editor)
    "thermal.mesa_diameter_um": {"unit": "um", "tag": "E", "lo": 0.1, "hi": 10.0,
                                  "source": "typical processed mesa diameter range"},
    "thermal.T_hs": {"unit": "K", "tag": "E", "lo": 4.0, "hi": 350.0,
                      "source": "cryostat heatsink operating range"},

    # ---- cavity
    "cavity.enabled": {"unit": "bool", "tag": "A", "lo": 0.0, "hi": 1.0,
                        "source": "design choice: cavity/waveguide present"},
    "cavity.type": {"unit": "enum", "tag": "A", "lo": 0.0, "hi": 1.0,
                     "source": "planar (F6 resonant) vs sin_waveguide (evanescent)"},
    "cavity.kappa": {"unit": "meV", "tag": "A", "lo": 0.1, "hi": 20.0,
                      "source": "MEEP/COMSOL cavity linewidth estimate"},
    "cavity.T_track": {"unit": "K", "tag": "A", "lo": 4.0, "hi": 350.0,
                        "source": "mode-tracking design target (F6 iv)"},
    "cavity.E_X0": {"unit": "eV", "tag": "E", "lo": 1.2, "hi": 2.2,
                     "source": "cryogenic X transition energy"},
    "cavity.dEdT_cav": {"unit": "meV/K", "tag": "E", "lo": -0.1, "hi": 0.0,
                         "source": "cavity mode refractive-index drift dn/dT"},
    "cavity.F_P": {"unit": "1", "tag": "A", "lo": 1.0, "hi": 30.0,
                   "source": "target Purcell factor"},
    "cavity.G": {"unit": "1", "tag": "A", "lo": 1.0, "hi": 30.0,
                 "source": "collection/rate gain on signal"},
    "cavity.beta_sin": {"unit": "1", "tag": "A", "lo": 0.0, "hi": 5.0,
                         "source": "SiN evanescent coupling coefficient (Lemma 1: brightness only)"},

    # ---- filter
    "filter.enabled": {"unit": "bool", "tag": "A", "lo": 0.0, "hi": 1.0,
                        "source": "design choice: slit filter present"},
    "filter.auto_w": {"unit": "bool", "tag": "A", "lo": 0.0, "hi": 1.0,
                       "source": "auto width = Gamma(T_j) operating convention"},
    "filter.w": {"unit": "meV", "tag": "A", "lo": 0.1, "hi": 20.0,
                 "source": "manual filter width when auto_w=False"},
    "filter.dx": {"unit": "meV", "tag": "A", "lo": -10.0, "hi": 10.0,
                  "source": "filter offset from the X window center"},

    # ---- aperture / ensemble
    "aperture.density_cm2": {"unit": "cm^-2", "tag": "E", "lo": 1e6, "hi": 1e10,
                              "source": "areal QD density estimate"},
    "aperture.diameter_um": {"unit": "um", "tag": "A", "lo": 0.1, "hi": 10.0,
                              "source": "aperture/collection diameter"},
    "aperture.sigma_inh": {"unit": "meV", "tag": "E", "lo": 5.0, "hi": 100.0,
                            "source": "inhomogeneous linewidth estimate"},
    "aperture.comp_brightness": {"unit": "1", "tag": "A", "lo": 0.0, "hi": 1.0,
                                  "source": "relative brightness of in-window competitor dots"},
}

# ---- envelope defaults (phase D2): default sweep ranges for [A]-tagged
# parameters, taken from the program's prediction cards (qcap-staged.yaml /
# qcap-piezo-variant.yaml). Program decision: [A] inputs default to RANGED
# on in the designer, rendering bands instead of single-point lines.
ENV_DEFAULTS = {
    "dot.delta_xx": (1.5, 5.0),
    "dot.gamma_scale": (0.7, 1.3),
    "drive.b_e": (1e-3, 0.3),
    "drive.mu": (0.1, 1.0),
    "drive.F_p": (0.3, 1.0),
    "drive.eta_capture": (0.1, 1.0),
    "cavity.kappa": (0.3, 3.0),
    "cavity.G": (5.0, 15.0),
    "cavity.F_P": (5.0, 20.0),
    "aperture.density_cm2": (1e8, 2e10),
}


_DEFAULT_RANGED_EXCLUDED = {"aperture.density_cm2", "drive.F_p", "drive.eta_capture"}


def default_ranged(design: DeviceDesign) -> dict:
    """The ENV_DEFAULTS subset applicable to `design`: cavity.* only makes
    sense (and is only wired into evaluate()) when the cavity is enabled and
    not a SiN waveguide (kappa/G/F_P are unused for cavity.type=sin_waveguide,
    see CavityBlock/evaluate). aperture.density_cm2 is excluded by default --
    it only feeds the F5 aperture-penalty scalar, not the curves, so ranging
    it by default would band nothing a user is looking at. drive.F_p and
    drive.eta_capture (F8/F8b, v1.1) are also excluded by default: F8's
    domain (mu >= 1 - f8b_thin_fano(eta_capture, F_p)) couples them to
    drive.mu, and ranging all three together by default risks a ValueError
    on any design whose default mu doesn't span the coupled band -- F_p/
    eta_capture ranging is opt-in (still available via ENV_DEFAULTS for a
    caller who has chosen a compatible mu), not default-on."""
    out = {}
    for path, bounds in ENV_DEFAULTS.items():
        if path in _DEFAULT_RANGED_EXCLUDED:
            continue
        block = path.split(".", 1)[0]
        if block == "cavity" and not (design.cavity.enabled
                                      and design.cavity.type != "sin_waveguide"):
            continue
        if path == "drive.mu" and design.drive.F_p != 1.0:
            # F8 domain floor (Opus v1.1 finding): with a fixed quiet pump,
            # mu below 1 - F_eff has no realizable cap-2 loading distribution.
            # Clamp the default mu range to the physical domain; if the whole
            # range is below the floor, don't range mu at all.
            F_eff = 1.0 - design.drive.eta_capture * (1.0 - design.drive.F_p)
            floor = (1.0 - F_eff) + 1e-9
            lo, hi = bounds
            if floor >= hi:
                continue
            bounds = (max(lo, floor), hi)
        out[path] = bounds
    return out


_EXCLUDED = {("thermal", "layers"), ("thermal", "substrate")}


def missing_meta() -> list:
    """Editable numeric/bool/str field paths of DeviceDesign's blocks that lack
    a META entry (thermal.layers/substrate excluded -- they have their own
    stack editor, not a scalar-parameter one)."""
    d = DeviceDesign()
    missing = []
    for f in dataclasses.fields(d):
        block = getattr(d, f.name)
        if not dataclasses.is_dataclass(block):
            continue  # DeviceDesign.name: not a block
        for bf in dataclasses.fields(block):
            if (f.name, bf.name) in _EXCLUDED:
                continue
            val = getattr(block, bf.name)
            if isinstance(val, (int, float, bool, str)):
                path = f"{f.name}.{bf.name}"
                if path not in META:
                    missing.append(path)
    return missing
