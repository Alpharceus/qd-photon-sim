"""Tier passthroughs of the composable engine (design 0, 1.4 "tier passthrough?", 4.2, P6).

The existing device tiers are exposed inside the composed world WITHOUT being
re-implemented:

  legacy            V-a class device (staged-device-design.yaml)
  rt_edge           InP edge-emitting room-temperature device (emission.type "edge")
  nitride_planar    InGaN/GaN planar vertical-cavity device (platform ingan_gan_planar)
  nitride_nanowire  InGaN/GaN nanowire, horizontal / vertical (platform ingan_gan_nanowire)
  fock              4 K resonant Fock-state source card (resonant_source.load_fock_card)

For the four DeviceDesign tiers `TIERS[id]["run"](card, overrides, T_grid)` loads
the DeviceDesign, applies the overrides through `device._set_path` (the same
META dotted paths the Designer edits) and returns `device.evaluate(...)`
UNCHANGED (bit-identical). The Fock tier's `run` returns
`resonant_source.bs_figure_of_merit` at the card's central point.

`fock_card_to_spec(path)` maps the Fock card onto a composed spec
(inas_gaas_sk emitter, generic_cavity, resonant_pulse / tpe, lens_fibre_det) so
that the generic engine can be checked against `bs_figure_of_merit`.
`spec_for_tier_card(card)` returns the {"tier", "card"} handoff descriptor the
Studio preset row uses. `tiers_payload()` is the JSON description P10 merges into
`registry_payload()["tiers"]`.

Only this module imports fsim_core.device, and only inside functions.
"""
from __future__ import annotations

import copy
from pathlib import Path

import yaml

from .params import widest_tag
from . import registry as _registry

ROOT = Path(__file__).resolve().parents[2]
CARDS_DIR = ROOT / "cards"
FOCK_CARD = "fock-source-4K-design"
SCHEMA = "fsim-composed/1"


# --------------------------------------------------------------------------
# card discovery / classification
# --------------------------------------------------------------------------
def _card_file(card) -> Path:
    """A card given as a path, a file name, or a bare card name in cards/."""
    p = Path(card)
    if p.suffix in (".yaml", ".yml"):
        if p.is_file():
            return p
        if (CARDS_DIR / p.name).is_file():
            return CARDS_DIR / p.name
        raise FileNotFoundError(f"card {str(card)!r} not found (looked at {p} and {CARDS_DIR / p.name})")
    q = CARDS_DIR / (p.name + ".yaml")
    if q.is_file():
        return q
    raise FileNotFoundError(f"card {str(card)!r} not found in {CARDS_DIR}")


def _tier_of_raw(doc: dict):
    """-> (tier id, variant | None) from the raw YAML mapping, or (None, None)
    for cards that are neither a DeviceDesign nor the Fock card."""
    if not isinstance(doc, dict):
        return None, None
    if (doc.get("meta") or {}).get("role") == "fock-source-design":
        return "fock", None
    d = doc.get("design")
    if not isinstance(d, dict):
        return None, None
    plat = d.get("platform", "legacy")
    if plat == "ingan_gan_planar":
        return "nitride_planar", None
    if plat == "ingan_gan_nanowire":
        fam = str(((d.get("nitride") or {}).get("nanowire") or {}).get("family", ""))
        return "nitride_nanowire", ("horizontal" if fam.startswith("horizontal")
                                    else "vertical" if fam.startswith("vertical") else None)
    if plat == "legacy":
        if (d.get("emission") or {}).get("type") == "edge":
            return "rt_edge", None
        return "legacy", None
    return None, None                      # composed or unknown: not a tier card


def _scan_cards() -> dict:
    """{card name: (tier id, variant)} for every cards/*.yaml that belongs to a tier."""
    out = {}
    for f in sorted(CARDS_DIR.glob("*.yaml")):
        try:
            doc = yaml.safe_load(f.read_text(encoding="utf-8"))
        except Exception:
            continue
        tier, variant = _tier_of_raw(doc)
        if tier is not None:
            out[f.stem] = (tier, variant)
    return out


def spec_for_tier_card(card):
    """UI preset handoff: {"tier": id, "card": name} (+ "variant" for the nanowire
    tier) for a shipped card path / name / loaded DeviceDesign, or None when the
    card is not a tier card (validation cards, prediction cards, composed cards)."""
    if hasattr(card, "platform") and hasattr(card, "name") and not isinstance(card, (str, Path)):
        return None if getattr(card, "platform", "") == "composed" else _desc_from_design(card)
    try:
        f = _card_file(card)
    except FileNotFoundError:
        return None
    tier, variant = _tier_of_raw(yaml.safe_load(f.read_text(encoding="utf-8")))
    if tier is None:
        return None
    out = {"tier": tier, "card": f.stem}
    if variant:
        out["variant"] = variant
    return out


def _desc_from_design(design):
    """Descriptor for an in-memory DeviceDesign (card name = design.name)."""
    doc = {"design": {"platform": design.platform, "emission": {"type": design.emission.type},
                      "nitride": design.nitride}}
    tier, variant = _tier_of_raw(doc)
    if tier is None:
        return None
    out = {"tier": tier, "card": design.name}
    if variant:
        out["variant"] = variant
    return out


# --------------------------------------------------------------------------
# passthrough runners
# --------------------------------------------------------------------------
def _load_design(card_path_or_design):
    from ..device import DeviceDesign
    if isinstance(card_path_or_design, DeviceDesign):
        return copy.deepcopy(card_path_or_design)
    return DeviceDesign.load(_card_file(card_path_or_design))


def _apply_overrides(design, overrides):
    from ..device import _set_path
    for path, value in (overrides or {}).items():
        if not isinstance(path, str) or "." not in path:
            raise KeyError(f"override path {path!r} is not a META dotted path 'block.field'")
        block, field = path.split(".", 1)
        blk = getattr(design, block, None)
        if blk is None or isinstance(blk, (dict, str, int, float, bool)) or not hasattr(blk, field):
            raise KeyError(f"override path {path!r} is not an editable DeviceDesign field")
        _set_path(design, path, value)
    return design


def run_design(card_path_or_design, overrides=None, T_grid=None) -> dict:
    """Load the DeviceDesign, apply META-path overrides, return
    `device.evaluate(design, T_grid)` unchanged. The input design object is never
    mutated."""
    from ..device import evaluate
    design = _apply_overrides(_load_design(card_path_or_design), overrides)
    return evaluate(design, T_grid)


def _make_design_runner(tier_id):
    def run(card_path_or_design, overrides=None, T_grid=None):
        if not hasattr(card_path_or_design, "platform"):
            doc = yaml.safe_load(_card_file(card_path_or_design).read_text(encoding="utf-8"))
            got, _ = _tier_of_raw(doc)
            if got != tier_id:
                raise ValueError(f"card {str(card_path_or_design)!r} belongs to tier {got!r}, not {tier_id!r}")
        return run_design(card_path_or_design, overrides, T_grid)
    run.__name__ = f"run_{tier_id}"
    return run


def run_fock(card_path=None, overrides=None, T_grid=None, N: int = 2) -> dict:
    """Fock tier: `resonant_source.bs_figure_of_merit` at the card's central point
    with `overrides` (card entry name -> value). T_grid is ignored (the card has
    one cryostat temperature, T_K)."""
    from .. import resonant_source as rs
    card = rs.load_fock_card(_card_file(card_path or FOCK_CARD))
    return rs.bs_figure_of_merit(rs.source_from_card(card, **(overrides or {})), N)


def _meta_paths() -> int:
    from ..design_meta import META
    return len(META)


_DESIGN_TIERS = (
    ("legacy", "Legacy V-a class",
     "Electrically injected InP-class quantum dot, incoherent pump, class-proxy confinement and "
     "transport (the V-a validated chain). Runs the whole legacy device evaluate()."),
    ("rt_edge", "InP edge RT",
     "Room-temperature electrically driven InP dot in an edge-emitting ridge "
     "(waveguide.edge_emission out-coupling; brightness only, never g2)."),
    ("nitride_planar", "Planar nitride",
     "InGaN/GaN dot in a planar vertical cavity: polar/nonpolar, QW-fluctuation and "
     "SET-injection variants."),
    ("nitride_nanowire", "Nanowire nitride (H/V)",
     "InGaN/GaN nanowire dot, horizontal as-built or vertical photonic, with wire thermal "
     "and injector models."),
)


def _build_tiers() -> dict:
    scan = _scan_cards()
    tiers = {}
    for tid, label, summary in _DESIGN_TIERS:
        cards = [c for c, (t, _v) in scan.items() if t == tid]
        entry = {"id": tid, "label": label, "summary": summary, "kind": "device_design",
                 "cards": cards, "handoff": "design", "run": _make_design_runner(tid)}
        if tid == "nitride_nanowire":
            entry["variants"] = {v: [c for c in cards if scan[c][1] == v]
                                 for v in ("horizontal", "vertical")}
        tiers[tid] = entry
    tiers["fock"] = {
        "id": "fock", "label": "Fock 4 K",
        "summary": "Cryogenic resonant (pi pulse) or two-photon Fock-state source card: boson-sampling "
                   "figure of merit from resonant_source. Maps onto a composed spec "
                   "(fock_card_to_spec) so the generic engine can reproduce it.",
        "kind": "fock_card", "cards": [c for c, (t, _v) in scan.items() if t == "fock"],
        "handoff": "compose", "run": run_fock}
    return tiers


TIERS = _build_tiers()


def tiers_payload() -> list:
    """JSON-serialisable tier list for registry_payload()["tiers"] (P10 merges it):
    [{id, label, summary, kind, cards, handoff, meta_paths, variants?}]. `meta_paths`
    is the number of editable DeviceDesign META paths for DeviceDesign tiers and
    the number of card entries for the Fock tier."""
    scan_fock = None
    out = []
    for tid, t in TIERS.items():
        if tid == "fock":
            from .. import resonant_source as rs
            if scan_fock is None:
                scan_fock = len(rs.load_fock_card(_card_file(FOCK_CARD)).params)
            n = scan_fock
        else:
            n = _meta_paths()
        d = {k: copy.deepcopy(t[k]) for k in ("id", "label", "summary", "kind", "cards", "handoff")}
        d["meta_paths"] = n
        if "variants" in t:
            d["variants"] = copy.deepcopy(t["variants"])
        out.append(d)
    return out


# --------------------------------------------------------------------------
# Fock card -> composed spec
# --------------------------------------------------------------------------
def _unit(axis, kind, name):
    return _registry.get(axis, kind).param(name).unit


def _leaf(axis, kind, name, value, tag, source):
    """Registry-checked leaf mapping {value, unit?, tag, source}."""
    unit = _unit(axis, kind, name)            # KeyError if the param does not exist
    out = {"value": value}
    if unit not in ("", "enum"):
        out["unit"] = unit
    out["tag"] = tag
    out["source"] = source
    return out


def _card_central(card, name):
    """(central value, tag, source text) of a card entry; ranges give the
    declared central point and say so (never a prediction)."""
    p = card.params[name]
    v = p.central()
    src = p.source
    if p.range is not None:
        how = "geometric mean" if p.scale == "log" else "midpoint"
        src = (f"{src} [Fock card range {list(p.range)} {p.unit} ({p.scale}); {how} used here, "
               f"a central point and not a prediction]")
    if isinstance(v, bool):
        pass
    elif isinstance(v, (int, float)):
        v = float(v)
    return v, p.tag, src


def fock_card_to_spec(path=None, overrides=None) -> dict:
    """Map the resonant Fock card onto a composed spec dict (schema fsim-composed/1).

    The card's ranges are evaluated at their declared central point (as
    `resonant_source.source_from_card`), so each spec leaf is a value that keeps
    the card's tag and source, with the range recorded in the source text.

      material  inas_gaas_sk (no generic two-level kind exists in materials.py):
                tau_bare = 1/gamma0 [DR], E_b = -delta_xx [DR sign convention],
                gamma* from resonant_source.source_rates [DR], a_ac = 0 (the card has
                no T dependence), FSS = 0 (not modelled by the card), alpha pinned
      structure generic_cavity: F_P, kappa, detuning, beta0, eta_out direct; F_bg = 1
      drive     resonant_pulse (target = transition) or tpe, per `excitation`
      setup     lens_fibre_det: eta_optics = eta_demux x eta_circuit, eta_det, zpl_filter

    `overrides` (card entry name -> value) as in source_from_card."""
    from .. import resonant_source as rs
    card = rs.load_fock_card(_card_file(path or FOCK_CARD))
    ov = dict(overrides or {})
    src = rs.source_from_card(card, **ov)

    def C(name):
        if name in ov:
            p = card.params[name]
            return ov[name], p.tag, f"{p.source} [overridden to {ov[name]!r}]"
        return _card_central(card, name)

    tpe = src.excitation == "tpe"
    rates = rs.source_rates(src)
    from ..qd_gf import PhononParams
    _PH = PhononParams()
    _PH_SRC = ("[DR] qd_gf.PhononParams default, the phonon cutoff geometry that resonant_source.phonon_Z "
               "uses for the Fock card (only alpha_ps2 is a card entry); pinned so Z(T) matches")

    g0, g0_tag, g0_src = C("gamma0_ns")
    dxx, dxx_tag, dxx_src = C("delta_xx_meV")
    m_tag = C("hom_width_multiple")[1]
    gs_tags = [m_tag, g0_tag, C("dephasing_mode")[1]]
    if src.dephasing_mode == "enhanced":
        gs_tags += [C("F_P")[1], C("kappa_ueV")[1], C("cav_detune_ueV")[1]]
    mat_kind = "inas_gaas_sk"
    material = {"kind": mat_kind, "params": {
        "tau_bare_ns": _leaf("material", mat_kind, "tau_bare_ns", 1.0 / float(src.gamma0_ns), g0_tag,
                             f"[DR] tau_bare = 1/gamma0 of the Fock card entry gamma0_ns. {g0_src}"),
        "gamma_xx_ratio": _leaf("material", mat_kind, "gamma_xx_ratio", float(src.gamma_xx_ratio),
                                *C("gamma_xx_ratio")[1:]),
        "E_b_meV": _leaf("material", mat_kind, "E_b_meV", -float(src.delta_xx_meV), dxx_tag,
                         f"[DR] E_b = -delta_xx_meV (card: negative = bound; material: positive = "
                         f"binding). {dxx_src}"),
        "gamma_star_4K_ns": _leaf(
            "material", mat_kind, "gamma_star_4K_ns", float(rates["gamma_star"]), widest_tag(*gs_tags),
            f"[DR] pure dephasing gamma* = (m-1) Gamma_ref from resonant_source.source_rates at the card "
            f"central point (dephasing_mode {src.dephasing_mode!r}, m = {src.hom_width_multiple:.6g}); "
            "the fixed point of the engine's self-consistent linewidth is the card's own"),
        "a_ac_ns_per_K": _leaf("material", mat_kind, "a_ac_ns_per_K", 0.0, "A",
                               "Fock card has no temperature dependence of gamma*; slope set to 0"),
        "fss_ueV": _leaf("material", mat_kind, "fss_ueV", 0.0, "A",
                         "Fock card does not model a fine-structure splitting"),
        "alpha_ps2": _leaf("material", mat_kind, "alpha_ps2", float(src.phonon_alpha_ps2),
                           *C("phonon_alpha_ps2")[1:]),
        # Z(T) of the card comes from qd_gf.PhononParams(alpha_ps2=alpha) with the OTHER phonon
        # fields at their qd_gf defaults; the material's own defaults differ (c_LA 5110, l_xy 4.0), so
        # the cutoff, and with it Z, would not match unless these are pinned to the same numbers.
        "c_LA_m_s": _leaf("material", mat_kind, "c_LA_m_s", float(_PH.c_s_m_s), "DR", _PH_SRC),
        "l_xy_nm": _leaf("material", mat_kind, "l_xy_nm", float(_PH.l_xy_nm), "DR", _PH_SRC),
        "l_z_nm": _leaf("material", mat_kind, "l_z_nm", float(_PH.l_z_nm), "DR", _PH_SRC),
        "sigma_sd_ueV": _leaf("material", mat_kind, "sigma_sd_ueV", float(src.sigma_sd_ueV),
                              *C("sigma_sd_ueV")[1:]),
        "tau_sd_ns": _leaf("material", mat_kind, "tau_sd_ns", float(src.tau_sd_ns), *C("tau_sd_ns")[1:]),
        "k_nr_ns": _leaf("material", mat_kind, "k_nr_ns", 0.0, "A",
                         "Fock card has no non-radiative loss channel"),
        "charge_state": _leaf("material", mat_kind, "charge_state", src.transition, *C("transition")[1:]),
    }}

    st_kind = "generic_cavity"
    structure = {"kind": st_kind, "params": {
        "F_P": _leaf("structure", st_kind, "F_P", float(src.F_P), *C("F_P")[1:]),
        "kappa_ueV": _leaf("structure", st_kind, "kappa_ueV", float(src.kappa_ueV), *C("kappa_ueV")[1:]),
        "Q": _leaf("structure", st_kind, "Q", 0.0, "A", "Q = 0: use the card's kappa directly"),
        "detuning_ueV": _leaf("structure", st_kind, "detuning_ueV", float(src.cav_detune_ueV),
                              *C("cav_detune_ueV")[1:]),
        "beta0": _leaf("structure", st_kind, "beta0", float(src.beta0), *C("beta0")[1:]),
        "eta_out": _leaf("structure", st_kind, "eta_out", float(src.eta_out), *C("eta_out")[1:]),
        "F_bg": _leaf("structure", st_kind, "F_bg", 1.0, "A",
                      "Fock card: no background Purcell suppression (total rate = Gamma0 (1 + F_P L))"),
    }}

    dr_kind = "tpe" if tpe else "resonant_pulse"
    dparams = {
        "shape": _leaf("drive", dr_kind, "shape", src.pulse_shape, *C("pulse_shape")[1:]),
        "fwhm_ps": _leaf("drive", dr_kind, "fwhm_ps", float(src.pulse_fwhm_ps), *C("pulse_fwhm_ps")[1:]),
        "area_pi": _leaf("drive", dr_kind, "area_pi", 1.0, "A",
                         "Fock card drives a pi pulse (resonant) or the TPE pi-equivalent area"),
        "detuning_meV": _leaf("drive", dr_kind, "detuning_meV", 0.0, "A",
                              "Fock card: laser on resonance (for TPE: on the two-photon resonance)"),
        "f_rep_MHz": _leaf("drive", dr_kind, "f_rep_MHz", float(src.f_rep_MHz), *C("f_rep_MHz")[1:]),
    }
    if not tpe:
        dparams["target"] = _leaf("drive", dr_kind, "target", src.transition, *C("transition")[1:])
    drive = {"kind": dr_kind, "params": dparams}

    se_kind = "lens_fibre_det"
    d_tag = widest_tag(C("eta_demux")[1], C("eta_circuit")[1])
    setup = {"kind": se_kind, "params": {
        "eta_optics": _leaf("setup", se_kind, "eta_optics", float(src.eta_demux) * float(src.eta_circuit),
                            d_tag, f"[DR] eta_demux x eta_circuit = {src.eta_demux:.6g} x "
                                   f"{src.eta_circuit:.6g} of the Fock card (central points)"),
        "eta_det": _leaf("setup", se_kind, "eta_det", float(src.eta_det), *C("eta_det")[1:]),
        "eta_pol": _leaf("setup", se_kind, "eta_pol", 1.0, "A", "Fock card efficiency chain has no polariser"),
        "eta_fibre": _leaf("setup", se_kind, "eta_fibre", 1.0, "A",
                           "Fock card efficiency chain has no separate fibre-coupling factor"),
        "zpl_filter": _leaf("setup", se_kind, "zpl_filter", bool(src.zpl_output), *C("zpl_output")[1:]),
    }}

    T_val, T_tag, T_src = C("T_K")
    return {"schema": SCHEMA,
            "T_K": {"value": float(T_val), "unit": "K", "tag": T_tag, "source": T_src},
            "material": material, "structure": structure, "drive": drive, "setup": setup,
            "solver": "auto"}


def fock_source_from_spec(spec) -> dict:
    """Inverse of fock_card_to_spec for round-trip checks: the card entries
    (`FockSource` field names) recovered from a mapped spec's leaf values, as a dict.
    eta_demux and eta_circuit are not separable (their product is `eta_optics`);
    hom_width_multiple is not stored (only the derived gamma*)."""
    s = spec.to_dict() if hasattr(spec, "to_dict") else spec

    def v(axis, name):
        leaf = s[axis]["params"][name]
        return leaf["value"] if isinstance(leaf, dict) else leaf

    tpe = s["drive"]["kind"] == "tpe"
    g0 = 1.0 / v("material", "tau_bare_ns")
    out = {
        "T_K": s["T_K"]["value"] if isinstance(s["T_K"], dict) else s["T_K"],
        "transition": v("material", "charge_state"),
        "excitation": "tpe" if tpe else "resonant_pi",
        "gamma0_ns": g0,
        "F_P": v("structure", "F_P"), "kappa_ueV": v("structure", "kappa_ueV"),
        "cav_detune_ueV": v("structure", "detuning_ueV"),
        "beta0": v("structure", "beta0"), "eta_out": v("structure", "eta_out"),
        "pulse_shape": v("drive", "shape"), "pulse_fwhm_ps": v("drive", "fwhm_ps"),
        "f_rep_MHz": v("drive", "f_rep_MHz"),
        "delta_xx_meV": -v("material", "E_b_meV"),
        "gamma_xx_ratio": v("material", "gamma_xx_ratio"),
        "sigma_sd_ueV": v("material", "sigma_sd_ueV"), "tau_sd_ns": v("material", "tau_sd_ns"),
        "phonon_alpha_ps2": v("material", "alpha_ps2"),
        "eta_det": v("setup", "eta_det"), "zpl_output": v("setup", "zpl_filter"),
        "eta_optics": v("setup", "eta_optics"),            # = eta_demux x eta_circuit (not separable)
        "gamma_star_ns": v("material", "gamma_star_4K_ns"),
    }
    return out
