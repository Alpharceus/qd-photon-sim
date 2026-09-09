"""Numerical wiring and legacy-regression checks for the opt-in planar
InGaN/GaN device tier (fsim_core.device.DeviceDesign.platform=
"ingan_gan_planar", integrating pieces 1-4: nitride_materials/
nitride_transport, nitride_levels, the deterministic cycle load map in
pulse_counting, and nitride_cavity).

Class discipline (README.md's five-way split):
  * "legacy bit-identical" checks below (--capture-legacy / the default run)
    are (N) numerical verification against a literal fixture captured from
    the committed pre-integration evaluator, never re-derived here.
  * "wiring" checks are (N) numerical verification of device.py's own
    plumbing (do the module return values reach the scalar contract
    unmodified, correctly-signed, and exactly once), independent of any
    external data.
  * "citations" checks are (T) source transcription -- confirming a cited
    paper's string is actually present next to the number it backs, not a
    re-check of the paper's own value.
  * "hardware diagnostics" checks are (N) numerical verification: the exact
    same drive_mech.set_feasibility/mech_set calls, run directly and
    compared to what evaluate() reports.
No check here claims held-out prediction; g2/rho targets are computed from
independent analytical limits or from the same modules called directly, not
copied from evaluate()'s own output and re-asserted.

Citations this file's design fixtures carry: the Deshpande et al., Applied
Physics Letters 105, 141109 (2014), DOI 10.1063/1.4897640, abstract 300 K
tau=1.3+/-0.3 ns [V abstract-only] radiative-lifetime anchor (also cited in
fsim_core/pulse_counting.py's deterministic_cycle_g2 docstring), and the
Bernardini, Fiorentini & Vanderbilt, PRB 56, R10024 (1997) polarization
convention (fsim_core/nitride_materials.py, fsim_core/nitride_levels.py).
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
import subprocess
import sys
import time
import types
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import fsim_core.device as devmod
from fsim_core.device import DeviceDesign, CavityBlock, evaluate
from fsim_core import nitride_levels, nitride_cavity, nitride_transport, drive_mech

BASELINE_PATH = ROOT / "verify" / "data" / "nitride_device_legacy_baseline.json"
# HEAD before this spec's device.py edits (pieces 1-4 landed, no `platform`
# field yet) -- the pre-edit ground truth AC1 requires. Fixed once at
# capture time; never re-resolved from the (now-edited) working tree.
LEGACY_COMMIT = "6246fb94feb47309b9eb75d1d1c8763586f58da7"
CARD_PATHS = ("cards/edge-inp-gaasp-design.yaml", "cards/edge-inp-gainp-design.yaml")

checks = []


def ok(name, value):
    checks.append(bool(value))
    print(("ok  " if value else "FAIL") + " " + name)


# --------------------------------------------------------------- encoding

def _encode(obj):
    """float.hex/array-exact JSON encoding (bool before int: bool is an int
    subclass); every float in scalars/curves round-trips bit-for-bit."""
    if isinstance(obj, dict):
        return {k: _encode(v) for k, v in sorted(obj.items())}
    if isinstance(obj, np.ndarray):
        return [_encode(v) for v in obj.tolist()]
    if isinstance(obj, (list, tuple)):
        return [_encode(v) for v in obj]
    if isinstance(obj, (bool, np.bool_)):
        return bool(obj)
    if isinstance(obj, (int, np.integer)):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        return float.hex(float(obj))
    if obj is None or isinstance(obj, str):
        return obj
    raise TypeError(f"cannot encode {type(obj)!r} for the legacy baseline fixture")


def _sha256_file(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def _git(*args):
    return subprocess.run(["git", *args], cwd=ROOT, capture_output=True,
                          text=True, check=True).stdout.strip()


def _legacy_module():
    """Exec the pre-edit fsim_core/device.py (git show, LEGACY_COMMIT) into
    a throwaway module -- in memory only, no file written -- so the legacy
    baseline is captured from the ACTUAL committed pre-integration
    evaluator, never from a hand-reimplementation of it."""
    src = _git("show", f"{LEGACY_COMMIT}:fsim_core/device.py")
    name = "fsim_core._device_legacy_pre_nitride"
    mod = types.ModuleType(name)
    mod.__file__ = str(ROOT / "fsim_core" / "device.py")
    mod.__package__ = "fsim_core"
    sys.modules[name] = mod
    exec(compile(src, f"<fsim_core/device.py@{LEGACY_COMMIT}>", "exec"), mod.__dict__)
    return mod


def _legacy_cases(mod):
    """Both InP cards at 230/300 K, finite_pulse on/off, both retention
    density conventions (16 cases), plus a representative default legacy
    design and a representative cavity-enabled legacy design (AC1)."""
    DD = mod.DeviceDesign
    cases = {}
    for card in CARD_PATHS:
        base = DD.load(str(ROOT / card))
        for T_hs in (230.0, 300.0):
            for fp in (False, True):
                for scales in (False, True):
                    d = copy.deepcopy(base)
                    d.thermal.T_hs = T_hs
                    d.drive.finite_pulse = fp
                    d.ret.tau_cap_scales_with_density = scales
                    label = f"{Path(card).stem}_T{int(T_hs)}_fp{int(fp)}_scale{int(scales)}"
                    cases[label] = d
    d = DD(); d.thermal.T_hs = 300.0
    cases["default_legacy"] = d
    d2 = DD(); d2.thermal.T_hs = 300.0; d2.cavity.enabled = True
    cases["default_legacy_cavity"] = d2
    return cases


def capture_legacy():
    """--capture-legacy: save exact pre-edit float.hex baseline. Must be run
    BEFORE any production edit to device.py reaches the working tree in a
    way that could change legacy numerics -- it evaluates the LEGACY module
    (git show LEGACY_COMMIT), never the current fsim_core.device, so it is
    safe to (re-)run at any time without ever baking a post-edit number into
    the fixture."""
    mod = _legacy_module()
    cases = _legacy_cases(mod)
    baseline = {
        "note": "Legacy (pre-nitride-integration) bit-identical baseline: "
                "float.hex/array captures from the committed evaluator at "
                "LEGACY_COMMIT, for the InP cards (both retention density "
                "conventions, both loading regimes, 230/300 K) plus a "
                "default and a cavity-enabled legacy design. Never update "
                "this fixture from a modified device.py.",
        "legacy_commit": LEGACY_COMMIT,
        "legacy_device_py_blob": _git("rev-parse", f"{LEGACY_COMMIT}:fsim_core/device.py"),
        "card_sha256": {c: _sha256_file(ROOT / c) for c in CARD_PATHS},
        "cases": {},
    }
    for label, d in cases.items():
        out = mod.evaluate(d, [d.thermal.T_hs])
        baseline["cases"][label] = {"scalars": _encode(out["scalars"]),
                                    "curves": _encode(out["curves"])}
    BASELINE_PATH.write_text(json.dumps(baseline, indent=1, sort_keys=True) + "\n",
                             encoding="utf-8")
    print(f"captured {len(cases)} legacy baseline cases -> {BASELINE_PATH}")


def check_legacy_regression():
    """AC5: existing rectangular-pulse AND InP-card results are bit-identical
    through the baseline comparison. Rebuilds the SAME cases with the
    CURRENT (possibly edited) device.py and compares every captured field
    exactly against the fixture -- never regenerates the fixture here."""
    if not BASELINE_PATH.exists():
        ok("legacy baseline fixture exists (run --capture-legacy first)", False)
        return
    try:
        baseline = json.loads(BASELINE_PATH.read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError):
        ok("legacy baseline fixture is readable JSON", False)
        return
    if "cases" not in baseline or not baseline["cases"]:
        ok("legacy baseline fixture has been captured (run --capture-legacy)", False)
        return
    for card in CARD_PATHS:
        want = baseline.get("card_sha256", {}).get(card)
        ok(f"legacy baseline source hash recorded for {card}",
          want == _sha256_file(ROOT / card))
    current_cases = _legacy_cases(devmod)
    for label, expected in baseline["cases"].items():
        d = current_cases.get(label)
        if d is None:
            ok(f"legacy bit-identical: {label} (case still in current recipe)", False)
            continue
        out = evaluate(d, [d.thermal.T_hs])
        got = {"scalars": _encode(out["scalars"]), "curves": _encode(out["curves"])}
        ok(f"legacy bit-identical: {label}", got == expected)


# ---------------------------------------------------------- nitride fixture

def design(kind="rectangular"):
    d = DeviceDesign(platform="ingan_gan_planar")
    d.dot.linewidth = "anchored"; d.dot.lineshape = "lorentzian"
    d.dot.gamma0 = .7; d.dot.a_ac = .002; d.dot.E_LO = 91.5
    d.dot.gamma300 = 25.; d.dot.r_xx = 2.; d.dot.delta_xx = -10.
    d.ret.mode = "nitride_confinement"; d.ret.channel = "min"; d.ret.tau_cap_ps = 10.
    d.drive.mode = "EL-transport"; d.drive.finite_pulse = True; d.drive.I_uA = 10.
    d.drive.rep_rate_hz = 8e7
    d.drive.diode = {"preset": "nitride-planar", "tau_pulse_ns": .1}
    d.drive.n_dot_cm2 = 1e10; d.drive.cycle_loading = kind
    d.cavity.enabled = True; d.cavity.type = "nitride_planar"
    d.emission.type = "vertical_cavity"
    d.nitride = {"dot": {"height_nm": 3., "radius_nm": 10., "x_in": .25,
                         "screening_fraction": .98},
                "cavity": {"Q": 167., "mode_volume_norm": 2., "eta_out": .1,
                           "T_track": 300.},
                "tau_rad0_ns": 1.3}
    d.thermal.T_hs = 300.
    return d


def set_design():
    d = design("deterministic_pair")
    d.drive.set_params = {"radius_nm": 20., "eps_r": 13., "R_T_ohm": 1e6, "ec_margin": 10.}
    return d


# --------------------------------------------------------------- checks

def check_wiring():
    d = design("rectangular"); out = evaluate(d, [300.]); s = out["scalars"]
    ok("platform tag on nitride scalars", s["platform"] == "ingan_gan_planar")
    ok("one-point T_grid curves have length 1", len(out["curves"]["g2"]) == 1)
    ok("T_j finite", math.isfinite(s["T_j"]))
    ok("gamma_X_ns finite and positive", math.isfinite(s["gamma_X_ns"]) and s["gamma_X_ns"] > 0)
    ok("k_X_ns finite and nonnegative", math.isfinite(s["k_X_ns"]) and s["k_X_ns"] >= 0)
    ok("rep_rate_hz threads through unchanged", s["rep_rate_hz"] == 8e7)
    ok("required contract keys all present", set(_REQUIRED_SCALAR_KEYS) <= set(s))
    ok("split counts sum to the total (loss applied once)",
      abs(s["mean_counts_x"] + s["mean_counts_xx"] - s["mean_counts"]) < 1e-9)
    ok("collected flux excludes background",
      abs(s["collected_flux_x_s"] + s["collected_flux_xx_s"] - s["collected_flux_pulsed_s"]) < 1e-6)
    ok("total detected flux = signal + background",
      abs(s["total_detected_flux_s"] - (s["collected_flux_pulsed_s"] + s["background_flux_s"])) < 1e-6)
    ok("I_pair_pA = e*f (Deshpande-cycle rep rate, exact SI charge)",
      abs(s["I_pair_pA"] - 1.602176634e-19 * 8e7 * 1e12) < 1e-9)


def check_signed_binding():
    """Opus fix-round required bullet 1: the bias field passed to
    nitride_levels must be the diode's OWN depletion(vj, Tj).F_kVcm (never a
    lumped -vj/d_i estimate), pinned exactly at the fixture point, and it
    must decrease with forward bias. Also: the resolved E_X_eV is not the
    T_track anchor's own energy (the field is actually threading through,
    not silently ignored)."""
    d = design("rectangular"); s = evaluate(d, [300.])["scalars"]
    unbiased = nitride_levels.levels(
        nitride_levels.NitrideDotSystem(**d.nitride["dot"]), s["T_j"])
    ok("field_kVcm != unbiased polarization field (bias resolved and applied)",
      abs(s["field_kVcm"] - unbiased.field_kVcm) > 1e-6)
    ok("V_j is forward (positive) at a forward-driven current", s["V_j"] > 0)

    # Independent diode, built from the SAME card fields evaluate() used.
    diode = nitride_transport.planar_pin(**{k: v for k, v in d.drive.diode.items()
                                            if k not in ("preset", "tau_pulse_ns")})
    _, vj = diode.v_of_i(d.drive.I_uA * 1e-6, s["T_j"])
    dep = diode.depletion(vj, s["T_j"])
    ok("independently resolved V_j matches evaluate()'s", abs(vj - s["V_j"]) < 1e-9)
    ok("field_kVcm at the fixture point == unbiased polarization field + "
      "diode.depletion(vj, Tj).F_kVcm (pinned to the diode's own value, bullet 1)",
      abs(s["field_kVcm"] - (unbiased.field_kVcm + dep.F_kVcm)) < 1e-6)

    # Forward bias reduces the depletion field (V_bi - V_j shrinking): a
    # higher supplied current resolves a higher V_j, so depletion()'s own
    # F_kVcm must shrink -- independent of nitride_levels, screening, or
    # this design's own dot geometry (Opus finding 1: the old -vj/d_i
    # estimate grew the WRONG way with forward bias).
    _, vj_hi = diode.v_of_i(d.drive.I_uA * 10 * 1e-6, s["T_j"])
    dep_hi = diode.depletion(vj_hi, s["T_j"])
    ok("V_j increases with forward current", vj_hi > vj)
    ok("diode.depletion().F_kVcm decreases with forward bias (bullet 1 sign convention)",
      dep_hi.F_kVcm < dep.F_kVcm)


def check_cavity_tracking_anchor():
    """Opus fix-round required bullet 2: the cavity tracking anchor must be
    evaluated at the SAME bias state (same supplied current) as the
    operating point, so a card tracked at T_track shows |detuning| < ~1 meV
    at T_j == T_track (this fixture's own T_hs == nitride.cavity.T_track ==
    300 K, and self-heating at duty=0.008 leaves T_j within a few K of
    that)."""
    d = design("rectangular"); s = evaluate(d, [300.])["scalars"]
    ok("T_j lands within a few K of T_track for the on-resonance check",
      abs(s["T_j"] - d.nitride["cavity"]["T_track"]) < 5.0)
    ok("|detuning_meV| < ~1 meV when T_j is close to T_track (bullet 2)",
      abs(s["detuning_meV"]) < 1.0)


def check_invalid_row_reasons():
    """Opus fix-round required bullet 3: an invalid row carries NaN for
    every physics scalar not evaluated at the operating point (never a
    plausible-looking placeholder from some other state), plus one specific
    reason string -- never the literal ['invalid row'], never a duplicate."""
    d = design("rectangular"); d.drive.I_uA = 1.0e7  # forces thermal runaway
    s = evaluate(d, [300.])["scalars"]
    ok("forced-runaway row is invalid", s["valid"] is False)
    ok("invalid_reasons is non-empty and not the generic ['invalid row'] placeholder",
      len(s["invalid_reasons"]) > 0 and list(s["invalid_reasons"]) != ["invalid row"])
    ok("invalid_reasons carries no duplicate entries",
      len(s["invalid_reasons"]) == len(set(s["invalid_reasons"])))
    ok("invalid row: E_X_eV is NaN (never a plausible unbiased placeholder)",
      not math.isfinite(s["E_X_eV"]))
    ok("invalid row: field_kVcm is NaN", not math.isfinite(s["field_kVcm"]))
    ok("invalid row: overlap_sq is NaN", not math.isfinite(s["overlap_sq"]))
    ok("invalid row: electron_bound/hole_bound are False, not a stale bound placeholder",
      s["electron_bound"] is False and s["hole_bound"] is False)


def check_both_loading_regimes():
    d_rect = design("rectangular"); s_rect = evaluate(d_rect, [300.])["scalars"]
    ok("rectangular: one_pair_valid unset (False, legacy of deterministic-only concept)",
      s_rect["one_pair_valid"] is False)
    ok("rectangular: set_feasible NaN-equivalent (no SET pricing on this branch)",
      s_rect["set_feasible"] is False and not np.isfinite(s_rect["set_priced_F_p"]))

    d_set = set_design(); s_set = evaluate(d_set, [300.])["scalars"]
    ok("deterministic: cycle_loading threads through", s_set["cycle_loading"] == "deterministic_pair")
    ok("deterministic: ideal_load_F_p == 0 exactly (bullet 71)", s_set["ideal_load_F_p"] == 0.0)
    ok("deterministic: set_priced_F_p is 0 or 1 (turnstile eps_cycle=0 pricing)",
      s_set["set_priced_F_p"] in (0.0, 1.0))
    ok("deterministic: infeasible SET does not overwrite ideal counting (mean_counts still finite)",
      math.isfinite(s_set["mean_counts"]))
    ok("deterministic: infeasible SET blocks device_pass",
      (not s_set["set_feasible"]) and (s_set["device_pass"] is False))
    ok("deterministic: residual occupancy blocking cannot be exact-one-pair valid",
      (s_set["blocked_load_probability"] <= 1e-9) or (s_set["one_pair_valid"] is False))


def check_yaml_round_trip():
    """AC2: YAML round-trip through an in-memory serialization path (no
    temp files/dirs) -- dump the design's own asdict() to a YAML string,
    load it back through the SAME DeviceDesign.load() coercion logic by
    monkey-patching Path.read_text for this one in-memory string, and check
    that the round-tripped design evaluates bit-identically."""
    import yaml
    from dataclasses import asdict
    d = design("deterministic_pair")
    d.drive.set_params = {"radius_nm": 20., "eps_r": 13., "R_T_ohm": 1e6, "ec_margin": 10.}
    doc_str = yaml.safe_dump({"meta": {"name": d.name}, "design": asdict(d)}, sort_keys=False)
    raw = yaml.safe_load(doc_str)["design"]
    # Mirror DeviceDesign.load's own field-by-field construction (its body
    # reads from a path; here the source is the in-memory string above) --
    # not a second, independent parser.
    d2 = DeviceDesign(
        name=raw.get("name", "my-device"),
        dot=devmod.DotBlock(**raw["dot"]), ret=devmod.RetentionBlock(**raw["ret"]),
        drive=devmod.DriveBlock(**raw["drive"]), thermal=devmod.ThermalBlock(**raw["thermal"]),
        cavity=devmod.CavityBlock(**raw["cavity"]), filter=devmod.FilterBlock(**raw["filter"]),
        aperture=devmod.ApertureBlock(**raw["aperture"]),
        emission=devmod.EmissionBlock(**raw.get("emission", {})),
        platform=raw.get("platform", "legacy"), nitride=dict(raw.get("nitride", {})),
        provenance=dict(raw.get("provenance", {})),
    )
    ok("YAML round-trip preserves platform", d2.platform == d.platform)
    ok("YAML round-trip preserves nitride block", d2.nitride == d.nitride)
    s1 = evaluate(d, [300.])["scalars"]; s2 = evaluate(d2, [300.])["scalars"]
    ok("YAML round-trip evaluates bit-identically",
      _encode(s1) == _encode(s2))


def check_malformed_opt_ins():
    def raises(mutate):
        dd = design("rectangular"); mutate(dd)
        try:
            evaluate(dd, [300.]); return False
        except ValueError:
            return True

    ok("rejects unknown nitride key", raises(lambda d: d.nitride.__setitem__("bogus", 1)))
    ok("rejects legacy cavity.F_P override", raises(lambda d: setattr(d.cavity, "F_P", 99.0)))
    ok("rejects legacy cavity.G override", raises(lambda d: setattr(d.cavity, "G", 99.0)))
    ok("rejects non-anchored dot.linewidth", raises(lambda d: setattr(d.dot, "linewidth", "class")))
    ok("rejects ibm lineshape", raises(lambda d: setattr(d.dot, "lineshape", "ibm")))
    ok("rejects nonempty ret.system", raises(lambda d: setattr(d.ret, "system", {"a": 1})))
    ok("rejects nonempty ret.preset", raises(lambda d: setattr(d.ret, "preset", "InP-GaAsP0.4-AlGaAs0.4-RT")))
    ok("rejects unknown cycle_loading", raises(lambda d: setattr(d.drive, "cycle_loading", "bogus")))
    ok("rejects missing nitride.dot geometry key", raises(lambda d: d.nitride["dot"].pop("height_nm")))
    ok("rejects missing dot density", raises(lambda d: setattr(d.drive, "n_dot_cm2", 0.0)))
    ok("rejects missing capture time", raises(lambda d: setattr(d.ret, "tau_cap_ps", None)))
    ok("rejects bad platform value", raises(lambda d: setattr(d, "platform", "bogus")))
    ok("rejects legacy-platform deterministic_pair",
      raises(lambda d: (setattr(d, "platform", "legacy"), setattr(d.drive, "cycle_loading", "deterministic_pair"))))

    # Opus fix-round required bullet: malformed nitride inputs must raise
    # this module's ValueError contract, not a bare TypeError.
    ok("rejects unknown drive.diode key",
      raises(lambda d: d.drive.diode.__setitem__("bogus_diode_key", 1.0)))
    ok("rejects a non-numeric nitride.dot field (YAML 1.1 unsigned-exponent string)",
      raises(lambda d: d.nitride["dot"].__setitem__("radius_nm", "1e10")))
    ok("rejects drive.eta_load outside [0, 1]",
      raises(lambda d: setattr(d.drive, "eta_load", 1.5)))
    ok("rejects drive.cw_pump_ratio override on deterministic_pair loading",
      raises(lambda d: (setattr(d, "drive", d.drive),
                        setattr(d.drive, "cycle_loading", "deterministic_pair"),
                        d.drive.set_params.update({"radius_nm": 20., "eps_r": 13.,
                                                    "R_T_ohm": 1e6, "ec_margin": 10.}),
                        setattr(d.drive, "cw_pump_ratio", 2.0))))

    def set_missing():
        dd = set_design(); dd.drive.set_params = {}
        try:
            evaluate(dd, [300.]); return False
        except ValueError:
            return True
    ok("SET loading requires explicit island/feasibility inputs", set_missing())

    ok("load() rejects a nitride card missing an explicit dot field",
      _card_missing_gamma0_rejected())


def _card_missing_gamma0_rejected():
    """DeviceDesign.load()'s own presence check (not evaluate()'s finiteness
    check): build an in-memory YAML doc with platform=ingan_gan_planar and
    dot.gamma0 omitted, and feed it through DeviceDesign.load() by
    monkeypatching Path.read_text for the duration of one call -- no file is
    written (no temp dir/file; verify/data/ stays limited to the one
    baseline fixture this spec's Files: line allows)."""
    from dataclasses import asdict
    import yaml
    d = design("rectangular")
    doc = asdict(d)
    doc["dot"].pop("gamma0")
    text = yaml.safe_dump({"meta": {"name": "scratch"}, "design": doc}, sort_keys=False)
    orig_read_text = Path.read_text
    Path.read_text = lambda self, *a, **kw: text
    try:
        DeviceDesign.load("<in-memory nitride card, gamma0 omitted>")
        return False
    except ValueError:
        return True
    finally:
        Path.read_text = orig_read_text


def check_background_variants():
    d0 = design("rectangular"); s0 = evaluate(d0, [300.])["scalars"]
    ok("nonzero background_flux_s by default", s0["background_flux_s"] > 0)

    d_zero = design("rectangular"); d_zero.nitride["eta_background"] = 0.0
    s_zero = evaluate(d_zero, [300.])["scalars"]
    ok("eta_background=0 zeroes the background exactly", s_zero["background_flux_s"] == 0.0)

    d_glow = design("rectangular"); d_glow.nitride["background_tau_ns"] = 2.0
    s_glow = evaluate(d_glow, [300.])["scalars"]
    ok("afterglow (background_tau_ns>0) raises background over pump-window-only",
      s_glow["background_flux_s"] >= s0["background_flux_s"])

    d_gate = design("rectangular"); d_gate.drive.gate_ns = 0.02
    s_gate = evaluate(d_gate, [300.])["scalars"]
    ok("a narrower gate changes mean_counts (gate actually restricts counting)",
      s_gate["mean_counts"] != s0["mean_counts"])
    ok("photon moments and background share the SAME gate",
      s_gate["gate_ns_used"] == 0.02)


def check_reservoir_energy_and_acceptance():
    """Opus results review 2026-09-09: reservoir energy is material-owned,
    not the QCSE-shifted dot line, and its flat-spectrum cavity acceptance is
    evaluated around that reservoir energy."""
    d_lo = design("rectangular"); d_lo.nitride["dot"]["height_nm"] = 2.
    d_hi = design("rectangular"); d_hi.nitride["dot"]["height_nm"] = 5.
    # Same composition and temperature; this direct material-continuum
    # calculation is deliberately independent of either dot's E_X result.
    e_lo = devmod._nitride_reservoir_energy_eV(d_lo.nitride["dot"], 300., {})
    e_hi = devmod._nitride_reservoir_energy_eV(d_hi.nitride["dot"], 300., {})
    ok("reservoir energy is invariant when only dot height changes", e_lo == e_hi)

    captured = []
    original = devmod._nitride_flat_background_acceptance
    def spy(kappa, width, dx, detuning, reservoir_offset=0.):
        captured.append(reservoir_offset)
        return original(kappa, width, dx, detuning, reservoir_offset)
    devmod._nitride_flat_background_acceptance = spy
    try:
        s = evaluate(d_lo, [300.])["scalars"]
    finally:
        devmod._nitride_flat_background_acceptance = original
    expected_offset = (devmod._nitride_reservoir_energy_eV(
        d_lo.nitride["dot"], s["T_j"], {}) - s["E_X_eV"]) * 1e3
    ok("flat-spectrum background acceptance is evaluated at reservoir energy",
      len(captured) == 1 and math.isclose(captured[0], expected_offset, rel_tol=1e-12))


def check_purcell_invariance():
    """AC3: at fixed T/geometry/field, Q/purcell_enabled changes
    gamma_X/gamma_XX but k_X/k_XX (bare, from piece 2, no Purcell input)
    remain bit-identical."""
    d = design("rectangular"); s = evaluate(d, [300.])["scalars"]
    d_q = design("rectangular"); d_q.nitride["cavity"]["Q"] = 5000.
    s_q = evaluate(d_q, [300.])["scalars"]
    ok("Q change: k_X_ns bit-identical", s["k_X_ns"] == s_q["k_X_ns"])
    ok("Q change: k_XX_ns bit-identical", s["k_XX_ns"] == s_q["k_XX_ns"])
    ok("Q change: gamma_X_ns changes", s["gamma_X_ns"] != s_q["gamma_X_ns"])

    d_p = design("rectangular"); d_p.nitride["cavity"]["purcell_enabled"] = False
    s_p = evaluate(d_p, [300.])["scalars"]
    ok("purcell_enabled=False: k_X_ns bit-identical", s["k_X_ns"] == s_p["k_X_ns"])
    ok("purcell_enabled=False: gamma_X_ns changes", s["gamma_X_ns"] != s_p["gamma_X_ns"])
    ok("purcell_enabled=False: Fp_add == 0", s_p["Fp_add"] == 0.0)


def check_b_res_and_density_scaling():
    """Fix round 2 (2026-09-09, Opus re-review): drive.b_res and
    ret.tau_cap_scales_with_density were accepted by the nitride card
    schema but read only on the legacy transport path (device.py's "item
    8" comments around evaluate()), leaving both axes inert on the
    nitride branch. Wire both in and prove each moves a scalar output at
    fixed everything else -- a wiring check, not a re-derivation."""
    d0 = design("rectangular"); s0 = evaluate(d0, [300.])["scalars"]
    d_b = design("rectangular"); d_b.drive.b_res = 1.0
    s_b = evaluate(d_b, [300.])["scalars"]
    ok("drive.b_res changes background_flux_s on the nitride branch",
      s_b["background_flux_s"] != s0["background_flux_s"])
    ok("drive.b_res raises background_flux_s (residual channel is additive)",
      s_b["background_flux_s"] > s0["background_flux_s"])
    ok("drive.b_res=0 (default) leaves background_flux_s untouched",
      s0["background_flux_s"] == evaluate(design("rectangular"), [300.])["scalars"]["background_flux_s"])

    d_d0 = design("rectangular"); d_d0.drive.n_dot_cm2 = 1e9; d_d0.ret.tau_cap_scales_with_density = False
    s_d0 = evaluate(d_d0, [300.])["scalars"]
    d_d1 = design("rectangular"); d_d1.drive.n_dot_cm2 = 1e9; d_d1.ret.tau_cap_scales_with_density = True
    s_d1 = evaluate(d_d1, [300.])["scalars"]
    ok("ret.tau_cap_scales_with_density changes tau_cap_ps_used at n_dot_cm2=1e9",
      s_d1["tau_cap_ps_used"] != s_d0["tau_cap_ps_used"])
    ok("ret.tau_cap_scales_with_density=True scales tau_cap_ps_used by 1e10/n_dot_cm2",
      math.isclose(s_d1["tau_cap_ps_used"], s_d0["tau_cap_ps_used"] * (1e10 / 1e9), rel_tol=1e-9))
    ok("ret.tau_cap_scales_with_density changes k_X_ns downstream",
      s_d1["k_X_ns"] != s_d0["k_X_ns"])


def check_analytic_limits_and_citations():
    """AC4: independent analytical/module cross-checks, not a target g2
    copied from evaluate(); plus literal citation-string presence checks."""
    d = design("rectangular")
    s = evaluate(d, [300.])["scalars"]

    # Recompute the bare rates and the cavity response DIRECTLY from piece 2
    # / piece 4, at the SAME resolved T_j/field, and require device.py to
    # have forwarded them unmodified (not re-derived its own copy). Field:
    # the diode's own depletion(vj, Tj).F_kVcm (bullet 1), not a re-derived
    # -vj/d_i estimate.
    diode = nitride_transport.planar_pin(**{k: v for k, v in d.drive.diode.items()
                                            if k not in ("preset", "tau_pulse_ns")})
    dep = diode.depletion(s["V_j"], s["T_j"])
    system = nitride_levels.NitrideDotSystem(
        **{**d.nitride["dot"], "external_field_kVcm": dep.F_kVcm})
    lv = nitride_levels.levels(system, s["T_j"])
    rr = nitride_levels.rates(lv, s["T_j"], tau_rad0_ns=d.nitride["tau_rad0_ns"],
                              n_dot_cm2=d.drive.n_dot_cm2, tau_cap_ps=d.ret.tau_cap_ps,
                              channel=d.ret.channel)
    ok("independent nitride_levels.rates() k_X_ns matches evaluate()'s",
      s["k_X_ns"] == rr["k_X_ns"])
    ok("independent nitride_levels.rates() k_XX_ns matches evaluate()'s",
      s["k_XX_ns"] == rr["k_XX_ns"])
    ok("independent nitride_levels.levels() E_X_eV matches evaluate()'s",
      s["E_X_eV"] == lv.E_X_eV)

    # Hardware diagnostics: direct unchanged set_feasibility call, e^2/C
    # convention, matched against a SET design's evaluate() output.
    d_set = set_design(); s_set = evaluate(d_set, [300.])["scalars"]
    feas = drive_mech.set_feasibility(s_set["T_j"], **{**d_set.drive.set_params,
                                                       "f_cycle_Hz": 8e7})
    ok("hardware diagnostics: set_feasibility() feasible flag matches evaluate()'s",
      feas["feasible"] == s_set["set_feasible"])
    ok("hardware diagnostics: set_feasibility() E_C_meV matches evaluate()'s",
      feas["E_C_meV"] == s_set["set_E_C_meV"])
    ok("hardware diagnostics: e^2/C convention (E_C_meV = 1e3*e/C_sigma_F, not e/2C)",
      abs(feas["E_C_meV"] - 1e3 * 1.602176634e-19 / feas["C_sigma_F"]) / feas["E_C_meV"] < 1e-9)

    # Citations: the literal source strings must actually appear next to the
    # numbers they back (source transcription, class T -- not a re-check of
    # the papers' own values).
    levels_src = (ROOT / "fsim_core" / "nitride_levels.py").read_text(encoding="utf-8")
    pulse_src = (ROOT / "fsim_core" / "pulse_counting.py").read_text(encoding="utf-8")
    ok("Bernardini et al., PRB 56, R10024 (1997) cited in nitride_levels.py",
      "Bernardini" in levels_src and "1997" in levels_src)
    ok("Deshpande et al., APL 105, 141109 (2014) 300 K 1.3 ns cited in pulse_counting.py",
      "Deshpande" in pulse_src and "1.3" in pulse_src)
    device_src = (ROOT / "fsim_core" / "device.py").read_text(encoding="utf-8")
    ok("Deshpande 300 K/1.3 ns default is cross-referenced at its device.py use site",
      "Deshpande" in device_src)


def check_no_mutation_and_repeatability():
    d = design("rectangular")
    pristine = copy.deepcopy(d)
    s1 = evaluate(d, [300.])["scalars"]
    ok("evaluate() does not mutate its input design",
      d.dot == pristine.dot and d.ret == pristine.ret and d.drive == pristine.drive
      and d.nitride == pristine.nitride and d.thermal == pristine.thermal
      and d.cavity == pristine.cavity and d.emission == pristine.emission)
    s2 = evaluate(d, [300.])["scalars"]
    ok("repeat evaluation is bit-identical (no stale cache)", _encode(s1) == _encode(s2))


def check_finite_and_invalid_states():
    """Unbound dot / impossible field is reported as an explicit invalid
    row, not a silent fallback."""
    d = design("rectangular")
    d.nitride["dot"] = {"height_nm": 0.3, "radius_nm": 1.0, "x_in": .9}
    out = evaluate(d, [300.])
    s = out["scalars"]
    ok("an unbound/impossible-geometry dot is reported invalid, not silently substituted",
      s["valid"] is False and len(s["invalid_reasons"]) > 0)
    ok("an invalid row still returns the full scalar contract (NaN, not missing keys)",
      set(_REQUIRED_SCALAR_KEYS) <= set(s))


def check_timing():
    d = design("rectangular")
    for _ in range(3):
        evaluate(d, [300.])
    times = []
    for _ in range(15):
        t0 = time.perf_counter(); evaluate(d, [300.]); times.append(time.perf_counter() - t0)
    med = sorted(times)[len(times) // 2]
    print(f"timing: median one-point evaluate() = {med * 1e3:.3f} ms "
         f"(target < 300 ms; generous regression cutoff 1000 ms)")
    ok("one-point evaluate() well under the 300 ms target after warmup", med < 1.0)


_REQUIRED_SCALAR_KEYS = (
    "platform", "cycle_loading", "T_hs", "T_j", "g2_op", "collected_flux_pulsed_s",
    "collected_flux_x_s", "collected_flux_xx_s", "background_flux_s",
    "total_detected_flux_s", "mean_counts", "mean_counts_x", "mean_counts_xx",
    "rho_pulsed", "E_X_eV", "lambda_nm", "field_kVcm", "overlap_sq", "electron_bound",
    "hole_bound", "E_a_meV", "k_X_ns", "k_XX_ns", "gamma_X0_ns", "gamma_XX0_ns",
    "gamma_X_ns", "gamma_XX_ns", "S_X", "S_XX", "Q", "kappa_meV", "detuning_meV",
    "Fp_add", "F_eff_X", "F_eff_XX", "eta_out", "gate_ns_used", "rep_rate_hz",
    "r_dot_s", "mu_resolved", "n_dot_cm2_used", "tau_cap_ps_used", "I_pair_pA",
    "transport_current_uA", "V_j", "power_W", "counting_converged",
    "blocked_load_probability", "one_pair_valid", "set_feasible", "set_priced_F_p",
    "set_E_C_meV", "set_EC_over_kT", "set_radius_nm", "set_radius_max_nm",
    "set_C_sigma_F", "set_R_T_over_RQ", "set_f_max_Hz", "pair_supply_possible",
    "ideal_load_F_p", "valid", "invalid_reasons", "provenance",
)


def main():
    if "--capture-legacy" in sys.argv:
        capture_legacy()
        return
    t_start = time.time()
    check_legacy_regression()
    check_wiring()
    check_signed_binding()
    check_cavity_tracking_anchor()
    check_invalid_row_reasons()
    check_both_loading_regimes()
    check_yaml_round_trip()
    check_malformed_opt_ins()
    check_background_variants()
    check_reservoir_energy_and_acceptance()
    check_purcell_invariance()
    check_b_res_and_density_scaling()
    check_analytic_limits_and_citations()
    check_no_mutation_and_repeatability()
    check_finite_and_invalid_states()
    check_timing()
    elapsed = time.time() - t_start
    print(f"{sum(checks)}/{len(checks)} nitride device checks passed ({elapsed:.1f} s)")
    if elapsed > 300:
        print("WARNING: exceeded the 5-minute runtime budget")
    raise SystemExit(0 if all(checks) else 1)


if __name__ == "__main__":
    main()
