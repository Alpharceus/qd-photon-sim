"""Source-transcription and numerical checks for nitride_nanowire_surface.

Class discipline (README.md five-way split): the S_cm_s=1e3 secondary
attribution and the 0.52 ensemble PL ratio are read from
verify/data/nitride_nanowire_anchors.yaml as two DISTINCT source
transcriptions -- this file fails if either is missing, or if 0.52 is
ever treated as an absolute single-dot IQE or baked into the module as a
fitted literal.  The unit-conversion and channel-independence checks are
independently derived arithmetic, computed here with plain literals, never
by calling surface_rates()/yield_ratio() to generate their own expected
value.  The yield_ratio-vs-0.52 comparison at the end is printed as a
non-gating observation only; no assertion requires it to match.
"""
from __future__ import annotations

import ast
import math
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.nitride_nanowire_surface import (
    NitrideNanowireSurfaceError,
    NitrideNanowireSurfaceParams,
    surface_rates,
    yield_ratio,
)

LEDGER = ROOT / "verify" / "data" / "nitride_nanowire_anchors.yaml"
MODULE_SOURCE = (ROOT / "fsim_core" / "nitride_nanowire_surface.py").read_text(encoding="utf-8")
# Collapsed to a single line so a documentation phrase that happens to wrap
# across source lines still matches a plain substring search below.
MODULE_SOURCE_FLAT = " ".join(MODULE_SOURCE.split())


def close(a, b, rtol=1e-9, atol=1e-12):
    return math.isclose(float(a), float(b), rel_tol=rtol, abs_tol=atol)


def expect(exc, fn):
    try:
        fn()
    except exc:
        return True
    return False


def main() -> int:
    checks = []

    # --- Acceptance 1: velocity attribution and PL ratio as DISTINCT ---
    # source transcriptions; 0.52 must never be an absolute IQE or a fit.
    doc = yaml.safe_load(LEDGER.read_text(encoding="utf-8"))
    # Ledger layout after the piece-1 fix round (Opus finding, 2026-09-14):
    # the velocity lives in its own [E]-tagged anchor (secondary attribution,
    # ref. 35 never read); the PL ratio stays in the [V] thermal/PL anchor.
    vel = doc["anchors"]["deshpande2013_surface_velocity"]
    anchor = doc["anchors"]["deshpande2013_thermal_and_pl"]
    checks.append(vel["value"]["S_cm_s"] == 1000.0)
    checks.append(anchor["value"]["ensemble_PL_300K_over_10K"] == 0.52)
    checks.append(vel["tag"] == "E" and anchor["tag"] == "V")
    # The two numbers are transcribed as separate fields of one anchor row,
    # not merged into a single quantity, and the ledger's own transfer note
    # already refuses the IQE/fit reading -- this line pins that text so a
    # future edit cannot quietly delete the caveat.
    checks.append("never a fitted S, activation energy, or lifetime" in anchor["transfer_notes"])
    checks.append(anchor["value"]["ensemble_PL_300K_over_10K"] != vel["value"]["S_cm_s"])
    # The module may DOCUMENT 0.52 (it does, as a disclaimed reference in
    # the docstrings), but it must never appear as an actual numeric literal
    # in the code (a fitted parameter or gating threshold).  Parse the AST
    # so a docstring/comment mention of "0.52" (just text inside a string)
    # cannot trip this check the way a plain substring search would.
    tree = ast.parse(MODULE_SOURCE)
    numeric_literals = {n.value for n in ast.walk(tree)
                         if isinstance(n, ast.Constant) and isinstance(n.value, (int, float))
                         and not isinstance(n.value, bool)}
    checks.append(0.52 not in numeric_literals)
    # The module must say, in its own words, that 0.52 is not a single-dot
    # IQE -- guards against a future edit relabelling the ensemble ratio.
    checks.append("not a single-dot internal quantum efficiency" in MODULE_SOURCE_FLAT)
    checks.append("never gates on it" in MODULE_SOURCE_FLAT)

    # --- Acceptance 2: independent unit conversion ---
    # k_side = 2*S/(R_cm) in s^-1, converted to ns^-1; computed here from
    # plain literals, not via surface_rates().
    S, R_nm = 1000.0, 12.5
    expected_k_side = 2.0 * S / (R_nm * 1e-7) * 1e-9
    checks.append(close(expected_k_side, 1.6))
    p_none = NitrideNanowireSurfaceParams(shell="none")
    r = surface_rates(p_none, core_radius_nm=R_nm, T_K=300.0)
    checks.append(close(r["k_side_ns"], expected_k_side))
    checks.append(close(r["k_side_ns"], 1.6))

    r_2R = surface_rates(p_none, core_radius_nm=2.0 * R_nm, T_K=300.0)
    checks.append(close(r_2R["k_side_ns"], r["k_side_ns"] / 2.0))
    checks.append(close(r_2R["k_surface_reservoir_ns"], r["k_surface_reservoir_ns"] / 2.0))
    checks.append(close(r_2R["k_surface_X_ns"], r["k_surface_X_ns"] / 2.0))
    checks.append(close(r_2R["k_surface_XX_ns"], r["k_surface_XX_ns"] / 2.0))

    p_2S = NitrideNanowireSurfaceParams(S_cm_s=2.0 * S, shell="none")
    r_2S = surface_rates(p_2S, core_radius_nm=R_nm, T_K=300.0)
    checks.append(close(r_2S["k_side_ns"], r["k_side_ns"] * 2.0))

    p_zero_S = NitrideNanowireSurfaceParams(S_cm_s=0.0, shell="none")
    r_zero_S = surface_rates(p_zero_S, core_radius_nm=R_nm, T_K=300.0)
    checks.append(r_zero_S["k_side_ns"] == 0.0)
    checks.append(r_zero_S["k_surface_reservoir_ns"] == 0.0)
    checks.append(r_zero_S["k_surface_X_ns"] == 0.0)
    checks.append(r_zero_S["k_surface_XX_ns"] == 0.0)

    p_zero_access = NitrideNanowireSurfaceParams(reservoir_access=0.0, occupied_dot_access=0.0, shell="none")
    r_zero_access = surface_rates(p_zero_access, core_radius_nm=R_nm, T_K=300.0)
    checks.append(r_zero_access["k_surface_reservoir_ns"] == 0.0)
    checks.append(r_zero_access["k_surface_X_ns"] == 0.0)
    checks.append(r_zero_access["k_surface_XX_ns"] == 0.0)
    checks.append(r_zero_access["k_side_ns"] > 0.0)  # the uncombined uniform rate is unaffected

    # --- Acceptance 3: channel independence and shell-multiplier semantics ---
    base = NitrideNanowireSurfaceParams(shell="none", reservoir_access=0.4, occupied_dot_access=0.9)
    r_base = surface_rates(base, core_radius_nm=R_nm, T_K=300.0)
    only_res_changed = NitrideNanowireSurfaceParams(shell="none", reservoir_access=0.9, occupied_dot_access=0.9)
    r_res = surface_rates(only_res_changed, core_radius_nm=R_nm, T_K=300.0)
    checks.append(close(r_res["k_surface_X_ns"], r_base["k_surface_X_ns"]))
    checks.append(close(r_res["k_surface_XX_ns"], r_base["k_surface_XX_ns"]))
    checks.append(not close(r_res["k_surface_reservoir_ns"], r_base["k_surface_reservoir_ns"]))

    only_dot_changed = NitrideNanowireSurfaceParams(shell="none", reservoir_access=0.4, occupied_dot_access=0.2)
    r_dot = surface_rates(only_dot_changed, core_radius_nm=R_nm, T_K=300.0)
    checks.append(close(r_dot["k_surface_reservoir_ns"], r_base["k_surface_reservoir_ns"]))
    checks.append(not close(r_dot["k_surface_X_ns"], r_base["k_surface_X_ns"]))

    # k_surface_XX = 2 * k_surface_X exactly.
    checks.append(close(r_base["k_surface_XX_ns"], 2.0 * r_base["k_surface_X_ns"]))

    # AlGaN multiplier applies exactly once: shell_multiplier_used equals the
    # raw shell_multiplier field (not squared/compounded), and the resulting
    # channel rates are exactly the shell=none rates times that one factor.
    p_algan = NitrideNanowireSurfaceParams(shell="AlGaN", shell_multiplier=0.1, reservoir_access=0.4, occupied_dot_access=0.9)
    r_algan = surface_rates(p_algan, core_radius_nm=R_nm, T_K=300.0)
    checks.append(r_algan["shell_multiplier_used"] == 0.1)
    checks.append(close(r_algan["k_surface_reservoir_ns"], r_base["k_surface_reservoir_ns"] * 0.1))
    checks.append(close(r_algan["k_surface_X_ns"], r_base["k_surface_X_ns"] * 0.1))
    # k_side_ns is the bare uniform-cylinder rate (shell-independent by
    # construction); only the channel-specific rates below it carry the
    # shell factor, applied exactly once.
    checks.append(close(r_algan["k_side_ns"], r_base["k_side_ns"]))

    # shell='none' retains unpassivated rates even if shell_multiplier is
    # left at a non-1.0 value: the AlGaN factor never leaks into shell=none.
    p_none_custom_mult = NitrideNanowireSurfaceParams(shell="none", shell_multiplier=0.1)
    r_none_custom = surface_rates(p_none_custom_mult, core_radius_nm=R_nm, T_K=300.0)
    checks.append(r_none_custom["shell_multiplier_used"] == 1.0)
    checks.append(close(r_none_custom["k_side_ns"], expected_k_side))

    # --- Acceptance 4: invalid/NaN rejection ---
    checks.append(expect(NitrideNanowireSurfaceError, lambda: NitrideNanowireSurfaceParams(S_cm_s=-1.0)))
    checks.append(expect(NitrideNanowireSurfaceError, lambda: NitrideNanowireSurfaceParams(S_cm_s=float("nan"))))
    checks.append(expect(NitrideNanowireSurfaceError, lambda: NitrideNanowireSurfaceParams(shell="sio2")))
    checks.append(expect(NitrideNanowireSurfaceError, lambda: NitrideNanowireSurfaceParams(shell_multiplier=1.5)))
    checks.append(expect(NitrideNanowireSurfaceError, lambda: NitrideNanowireSurfaceParams(shell_multiplier=-0.1)))
    checks.append(expect(NitrideNanowireSurfaceError, lambda: NitrideNanowireSurfaceParams(reservoir_access=-0.01)))
    checks.append(expect(NitrideNanowireSurfaceError, lambda: NitrideNanowireSurfaceParams(occupied_dot_access=1.01)))
    checks.append(expect(NitrideNanowireSurfaceError, lambda: surface_rates(p_none, core_radius_nm=0.0, T_K=300.0)))
    checks.append(expect(NitrideNanowireSurfaceError, lambda: surface_rates(p_none, core_radius_nm=-5.0, T_K=300.0)))
    checks.append(expect(NitrideNanowireSurfaceError, lambda: surface_rates(p_none, core_radius_nm=float("nan"), T_K=300.0)))
    checks.append(expect(NitrideNanowireSurfaceError, lambda: surface_rates(p_none, core_radius_nm=R_nm, T_K=0.0)))
    checks.append(expect(NitrideNanowireSurfaceError, lambda: surface_rates(p_none, core_radius_nm=R_nm, T_K=-10.0)))
    checks.append(expect(NitrideNanowireSurfaceError, lambda: surface_rates(p_none, core_radius_nm=R_nm, T_K=float("nan"))))
    checks.append(expect(NitrideNanowireSurfaceError, lambda: surface_rates("not-params", core_radius_nm=R_nm, T_K=300.0)))
    checks.append(expect(NitrideNanowireSurfaceError, lambda: yield_ratio(gamma_300_ns=-1.0, loss_300_ns=1.0, gamma_10_ns=1.0, loss_10_ns=1.0)))
    checks.append(expect(NitrideNanowireSurfaceError, lambda: yield_ratio(gamma_300_ns=1.0, loss_300_ns=1.0, gamma_10_ns=1.0, loss_10_ns=float("nan"))))

    # Independent competing-rate/quantum-yield limits, computed with plain
    # literals rather than by calling yield_ratio() to generate its own
    # expected value.
    g300, l300, g10, l10 = 0.5, 1.5, 1.0, 0.0
    y300 = g300 / (g300 + l300)
    y10 = g10 / (g10 + l10)
    expected_ratio = y300 / y10
    yr = yield_ratio(gamma_300_ns=g300, loss_300_ns=l300, gamma_10_ns=g10, loss_10_ns=l10)
    checks.append(yr["defined"] is True)
    checks.append(close(yr["yield_300K"], y300))
    checks.append(close(yr["yield_10K"], y10))
    checks.append(close(yr["ratio"], expected_ratio))

    # Zero reference (10 K) radiative yield -> undefined with a reason,
    # not a ZeroDivisionError and not silently 0/0.
    yr_zero_ref = yield_ratio(gamma_300_ns=0.5, loss_300_ns=0.5, gamma_10_ns=0.0, loss_10_ns=1.0)
    checks.append(yr_zero_ref["defined"] is False)
    checks.append(yr_zero_ref["ratio"] is None)
    checks.append(isinstance(yr_zero_ref["reason"], str) and len(yr_zero_ref["reason"]) > 0)

    # Zero total rate at either condition is also undefined, not a crash.
    yr_zero_total_10 = yield_ratio(gamma_300_ns=0.5, loss_300_ns=0.5, gamma_10_ns=0.0, loss_10_ns=0.0)
    checks.append(yr_zero_total_10["defined"] is False)
    yr_zero_total_300 = yield_ratio(gamma_300_ns=0.0, loss_300_ns=0.0, gamma_10_ns=0.5, loss_10_ns=0.5)
    checks.append(yr_zero_total_300["defined"] is False)

    passed = sum(checks)
    total = len(checks)

    # Non-gating observation only: this module's own uniform-cylinder rate
    # has no built-in temperature law (by design -- see module docstring),
    # so a proxy built ONLY from k_surface_X_ns at two temperatures is
    # trivially 1.0; it is printed to show the helper composes correctly,
    # not as a stand-in for the full 300K/10K device model (other pieces
    # add the temperature-dependent capture/intrinsic channels).  It is
    # printed next to, never asserted against, the 0.52 ensemble anchor.
    p_headline = NitrideNanowireSurfaceParams(shell="none")
    rr_300 = surface_rates(p_headline, core_radius_nm=R_nm, T_K=300.0)
    rr_10 = surface_rates(p_headline, core_radius_nm=R_nm, T_K=10.0)
    proxy = yield_ratio(gamma_300_ns=1.0, loss_300_ns=rr_300["k_surface_X_ns"],
                         gamma_10_ns=1.0, loss_10_ns=rr_10["k_surface_X_ns"])
    print("non-gating observation: surface-only proxy ratio at R=12.5 nm, shell=none "
          "(no T-dependence in this module by design): " + repr(proxy["ratio"]) +
          "; measured ensemble I(300K)/I(10K)=0.52 (Deshpande 2013, separate non-gating anchor, not compared here)")

    print(f"{passed}/{total} nitride nanowire surface checks passed")
    return 0 if passed == total else 1


if __name__ == "__main__":
    raise SystemExit(main())
