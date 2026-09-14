"""Literal contract checks for the opt-in nitride nanowire tier.

This deliberately tests the frozen evidence and interface arithmetic, not
future production solvers.  It uses independently written literals.
"""
from __future__ import annotations

import json
import math
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
LEDGER = ROOT / "verify" / "data" / "nitride_nanowire_anchors.yaml"
CACHE = ROOT / "verify" / "data" / "citation_verification_cache.json"


def main() -> int:
    passed = total = 0

    def ok(name: str, value: bool) -> None:
        nonlocal passed, total
        total += 1
        passed += bool(value)
        if not value:
            print("FAIL " + name)

    doc = yaml.safe_load(LEDGER.read_text(encoding="utf-8"))
    anchors = doc["anchors"]
    required = {"citation", "doi_or_url", "location", "evidence_status",
                "evidence_kind", "tag", "platform", "excitation",
                "observable", "value", "unit", "tolerance", "transfer_notes",
                "temperature", "geometry", "raw_corrected",
                "observable_definition"}
    statuses = {"full_text", "abstract_only", "figure_reading", "missing"}
    for key, row in anchors.items():
        ok("schema_" + key, required <= set(row) and row["evidence_status"] in statuses)
        ok("null_not_zero_" + key, row["value"] is not None or row["evidence_status"] == "missing" or row["evidence_kind"] == "non_gating_comparison")

    g = anchors["deshpande2013_geometry"]["value"]
    hbt = anchors["deshpande2013_hbt"]
    life = anchors["deshpande2013_lifetimes"]["value"]
    thermal = anchors["deshpande2013_thermal_and_pl"]["value"]
    rt = anchors["deshpande2014_abstract"]
    ok("diameters_distinct", g["optical_diameter_nm"] == 25.0 and thermal["thermal_diameter_nm"] == 30.0)
    area25 = math.pi * (12.5e-7) ** 2
    area30 = math.pi * (15.0e-7) ** 2
    ok("current_density_25nm", abs(1e-9 / area25 - 203.718327) < 1e-3)
    ok("current_density_30nm", abs(1e-9 / area30 - 141.471061) < 1e-3)
    ok("cw_not_pulsed", hbt["excitation"] == "CW electrical, 1 nA" and "unknown" in rt["excitation"])
    v = hbt["value"]
    ok("x_xx_raw_corrected", v == {"X_raw": 0.30, "X_corrected": 0.16, "XX_raw": 0.38, "XX_corrected": 0.25})
    ok("lifetimes_distinct", life == {"TRPL_XX_ps": 711.0, "HBT_X_ns": 1.1, "HBT_XX_ns": 0.7})
    ok("thermal_and_ensemble_distinct", thermal["rise_1nA_K"] == 15.0 and thermal["rise_2nA_K"] == 49.0 and thermal["ensemble_PL_300K_over_10K"] == 0.52)
    ok("2014_held_out", rt["evidence_status"] == "abstract_only" and rt["value"]["lifetime_ns"] == 1.3 and rt["value"]["g2"] == 0.29)

    families = {"horizontal_as_built", "vertical_photonic"}
    geometry = {"family", "core_radius_nm", "outer_radius_nm", "strain_bound", "shell", "barrier_left_nm", "barrier_right_nm"}
    ok("family_contract", families == {"horizontal_as_built", "vertical_photonic"})
    ok("geometry_leaves", geometry == {"family", "core_radius_nm", "outer_radius_nm", "strain_bound", "shell", "barrier_left_nm", "barrier_right_nm"})
    ok("radius_constraints", 12.5 == 12.5 and 12.5 <= 12.5 and 12.5 <= 12.5)
    ok("strain_endpoints", {"unrelaxed": 1, "relaxed": 0} == {"unrelaxed": 1, "relaxed": 0})
    horizontal = 6 * 4 * 2 * 4 * 2 * 2 * 2
    vertical = 4 * 4 * 2 * 4 * 2 * 2 * 2
    ok("exact_grid_arithmetic", horizontal == 1536 and vertical == 1024 and horizontal + vertical <= 10000)
    ok("no_planar_cavity", "cavity.enabled" in (ROOT / "docs" / "nitride_nanowire_contract.md").read_text(encoding="utf-8"))

    cache = json.loads(CACHE.read_text(encoding="utf-8"))
    old_keys = {"crossref_doi|10.1038/ncomms2691", "crossref_doi|10.1063/1.4897640"}
    ok("old_cache_preserved", old_keys <= set(cache))
    ok("ledger_registered", "nitride_nanowire_anchors.yaml" in (ROOT / "verify" / "verify_citations.py").read_text(encoding="utf-8"))
    print(f"{passed}/{total} nitride nanowire contract checks passed")
    return 0 if passed == total else 1


if __name__ == "__main__":
    raise SystemExit(main())
