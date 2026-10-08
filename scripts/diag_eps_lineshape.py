#!/usr/bin/env python3
"""Diagnostic script: biexciton-leak lineshape sensitivity at RT edge points.

This script recomputes epsilon (XX/X spectral transmission ratio) for the
edge-emitter cards at cryogenic and room temperatures, under both Lorentzian
and IBM lineshapes, sweeping the spectral filter width. The Lorentzian rows
must reproduce the RT edge-emitter review table to three decimals.

CARDS USED (in-memory copies; no files modified):
  - cards/edge-inp-gaasp-design.yaml (GaAsP, 300 K)
  - cards/edge-inp-gainp-design.yaml (GaInP, 230 K and 300 K)
  - cards/nitride-cavity-pulse-design.yaml (nitride, 300 K)

CRITICAL NOTE: The 300 K epsilon values are measurements of total
transmission through a spectral filter, not decomposed homogeneous-width
numbers. A measured 300 K spectrum is the decisive input for trusting any
300 K lineshape model. This script exists to document sensitivity to the
choice of lineshape (Lorentzian vs IBM) over the range of filter widths
that the system can probe.

IBM ROWS (lineshape == "ibm"): the code's own independent-boson path, exactly
as fsim_core/device.py does when dot.lineshape == "ibm" (device.py ~1959-1978):
  * pp = qd_gf.PhononParams(**dot.phonon). None of the three cards carries a
    dot.phonon block, so pp is the DEFAULT PhononParams() on every card (InP-
    class deformation potentials / density / sound speed, l_xy = 4.5 nm,
    l_z = 1.5 nm, [DR] class proxy, [A, null] until AFM/TEM). For the nitride
    card this is an InP-class proxy, NOT InGaN/GaN constants, and it is not a
    bound in either direction (InGaN/GaN phonon coupling may be stronger or
    weaker): read its ibm rows as an InP-class proxy, not a sensitivity bound
    and not a prediction. CAVEAT (review finding 7, 2026-10-08): the card's
    Gamma(T) is a TOTAL width, yet it is used as the ZPL width and the
    qd_gf phonon sideband is added on top, i.e. the sideband is counted twice
    (the Q4 'total' double count); the ibm rows overstate the broadening.
  * t_X  = qd_gf.ibm_transmission(0,        pp, T, Gamma,        w_meV=w)
    t_XX = qd_gf.ibm_transmission(-d_xx,    pp, T, r_xx * Gamma, w_meV=w)
    eps = t_XX / t_X (same top-hat acceptance and delta convention as
    spectral.epsilon; no cavity). ZPL widths carry the card's Gamma(T) unchanged.
    The sign of delta_xx is kept (the sideband is asymmetric).
  * Z_T = qd_gf.zpl_weight(pp, T) is the ZPL weight (1.0 for Lorentzian rows).
NARROW-FILTER DEFINITION (both lineshapes): the w -> 0 limit of eps, i.e. the
ratio of unit-area spectral densities at the filter center,
  Lorentzian: spectral.epsilon_narrow_filter (closed form),
  IBM: S_XX(+d_xx) / S_X(0), S = [Z*Lorentz(omega; Gamma_zpl) + S_sb(omega)]
       / (Z + int S_sb), the exact normalization inside ibm_transmission
       (omega = photon energy relative to the line's ZPL, blue positive).
Cross-check (run at every execution): the direct qd_gf eps equals
fsim_core.device.evaluate(in-memory ibm copy, auto_w=False, w=2 meV).eps_op at
the evaluated (Tj, Gamma) to 1e-12, and the narrow definition matches
ibm_transmission at w = 0.5 meV to 1% (smaller w is dominated by the 6 ueV FFT-grid mask quantisation) (consistency of the w -> 0 definition).

Author: Claude Haiku 4.5
Date: 2026-10-08
"""
from __future__ import annotations

import sys
import csv
import copy
import warnings
from pathlib import Path
import numpy as np

# Add repo root to path
ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import fsim_core.device as device
import fsim_core.linewidth as linewidth
import fsim_core.spectral as spectral
import fsim_core.qd_gf as qd_gf


def load_card_copy(card_name: str) -> device.DeviceDesign:
    """Load a card and return an independent in-memory copy."""
    return device.DeviceDesign.load(ROOT / "cards" / card_name)


def ibm_eps(delta_xx, pp, T, gam_x, gam_xx, w):
    """eps = t_XX/t_X via qd_gf.ibm_transmission, device.py conventions
    (X at filter center, XX offset by -delta_xx, top-hat full width w, no
    cavity). The sign of delta_xx is kept (nitride: -10 meV)."""
    t_x = float(qd_gf.ibm_transmission(0.0, pp, T, gam_x, w_meV=w))
    t_xx = float(qd_gf.ibm_transmission(0.0 - delta_xx, pp, T, gam_xx, w_meV=w))
    return t_xx / t_x


def _ibm_density(omega0, pp, T, gamma_zpl):
    """Unit-area spectral density at omega0 with ibm_transmission's own
    normalization [Z*Lorentz + S_sb] / (Z + int S_sb)."""
    omega, s_sb, z = qd_gf._spectrum_parts(pp, float(T), float(gamma_zpl))
    sb_area = float(np.trapezoid(s_sb, omega))
    s = z * float(qd_gf.lorentzian(np.array([omega0]), 0.0, gamma_zpl)[0]) \
        + float(np.interp(omega0, omega, s_sb, left=0.0, right=0.0))
    return s / (z + sb_area)


def ibm_eps_narrow(delta_xx, pp, T, gam_x, gam_xx):
    """w -> 0 limit: X density at the filter center (omega = 0) over the XX
    density at omega = +delta_xx (ibm_transmission's delta = -d_xx gives
    x = omega - d_xx = 0)."""
    return _ibm_density(delta_xx, pp, T, gam_xx) / _ibm_density(0.0, pp, T, gam_x)


def crosscheck_against_evaluate(card_filename):
    """Direct qd_gf eps vs fsim_core.device.evaluate on an in-memory ibm copy
    (auto_w=False, w=2 meV) at the evaluated (Tj, Gamma); asserts 1e-12.
    Also checks the narrow-filter definition against ibm_transmission."""
    d = copy.deepcopy(load_card_copy(card_filename))
    d.dot.lineshape = "ibm"
    d.filter.auto_w = False
    d.filter.w = 2.0
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        sc = device.evaluate(d)["scalars"]
    pp = qd_gf.PhononParams(**d.dot.phonon)
    tj, gam = float(sc["T_j_op"]), float(sc["gamma_op"])
    direct = ibm_eps(d.dot.delta_xx, pp, tj, gam, d.dot.r_xx * gam, 2.0)
    diff = abs(direct - float(sc["eps_op"]))
    assert diff < 1e-12, f"direct qd_gf eps {direct!r} != evaluate {sc['eps_op']!r}"
    num = ibm_eps(d.dot.delta_xx, pp, tj, gam, d.dot.r_xx * gam, 0.5)
    nar = ibm_eps_narrow(d.dot.delta_xx, pp, tj, gam, d.dot.r_xx * gam)
    assert abs(num - nar) < 1e-2 * nar, f"narrow definition {nar} vs w=0.5: {num}"
    return diff, direct, tj, gam


def compute_epsilon_row(
    card_name: str,
    card: device.DeviceDesign,
    temperature_K: float,
    lineshape: str,
    filter_width_meV: float | None,
) -> dict:
    """Compute one epsilon result row.

    Args:
        card_name: e.g., "edge-inp-gainp-design"
        card: DeviceDesign object
        temperature_K: e.g., 230.0 or 300.0
        lineshape: "lorentzian" or "ibm"
        filter_width_meV: filter width in meV, or None for narrow-filter limit

    Returns:
        Dictionary with keys: card, T_K, lineshape, filter_width_meV,
        gamma_x_meV, gamma_xx_meV, delta_xx_meV, eps_narrow, eps_at_w,
        provenance
    """

    # Extract dot parameters
    dot_params = card.dot

    # Determine which linewidth model to use
    if lineshape == "lorentzian":
        # Use the LinewidthParams from the card
        lp = linewidth.LinewidthParams(
            gamma0=getattr(dot_params, "gamma0", 0.25),
            a_ac=getattr(dot_params, "a_ac", 2.0e-3),
            E_LO=getattr(dot_params, "E_LO", 43.0),
            gamma300=getattr(dot_params, "gamma300", 12.0),
        )
        gam_x = linewidth.gamma_anchor(temperature_K, lp)
    elif lineshape == "ibm":
        # Same Gamma(T) anchor as the Lorentzian rows (device.py uses the
        # card's Gamma(T) for the ZPL widths unchanged).
        lp = linewidth.LinewidthParams(
            gamma0=getattr(dot_params, "gamma0", 0.25),
            a_ac=getattr(dot_params, "a_ac", 2.0e-3),
            E_LO=getattr(dot_params, "E_LO", 43.0),
            gamma300=getattr(dot_params, "gamma300", 12.0),
        )
        gam_x = linewidth.gamma_anchor(temperature_K, lp)
    else:
        raise ValueError(f"Unknown lineshape: {lineshape}")

    # XX linewidth is r_xx * X linewidth
    r_xx = getattr(dot_params, "r_xx", 0.72)
    gam_xx = r_xx * gam_x

    # Delta XX (should be negative for nitride, positive for edge)
    delta_xx = getattr(dot_params, "delta_xx", 7.0)

    z_t = 1.0
    if lineshape == "lorentzian":
        eps_narrow = spectral.epsilon_narrow_filter(
            abs(delta_xx), gam_x, gam_xx
        )
        if filter_width_meV is None:
            eps_at_w = eps_narrow
            w_label = "w->0"
        else:
            result = spectral.epsilon(
                abs(delta_xx), gam_x, gam_xx, w=filter_width_meV
            )
            eps_at_w = result.eps
            w_label = f"w={filter_width_meV:.3f}"
        prov = f"[DR] recomputed from {card_name} using spectral.epsilon"
    else:
        pp = qd_gf.PhononParams(**dot_params.phonon)
        z_t = qd_gf.zpl_weight(pp, temperature_K)
        eps_narrow = ibm_eps_narrow(delta_xx, pp, temperature_K, gam_x, gam_xx)
        if filter_width_meV is None:
            eps_at_w = eps_narrow
            w_label = "w->0"
        else:
            eps_at_w = ibm_eps(delta_xx, pp, temperature_K, gam_x, gam_xx,
                               filter_width_meV)
            w_label = f"w={filter_width_meV:.3f}"
        prov = (f"[DR] qd_gf.ibm_transmission (device.py ibm path), default "
                f"PhononParams (InP-class proxy, not a bound; card Gamma(T) is a total width "
                f"used as ZPL width + sideband = Q4 double count), from {card_name}")

    return {
        "card": card_name,
        "T_K": temperature_K,
        "lineshape": lineshape,
        "filter_width_meV": filter_width_meV if filter_width_meV is not None else float("nan"),
        "w_label": w_label,
        "gamma_x_meV": gam_x,
        "gamma_xx_meV": gam_xx,
        "delta_xx_meV": delta_xx,
        "eps_narrow": eps_narrow,
        "eps_at_w": eps_at_w,
        "Z_T": z_t,
        "provenance": prov,
    }


def main():
    """Run the diagnostic sweep."""

    # Create output directory
    out_dir = ROOT / "out" / "diag"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Cards to process
    cards_to_load = [
        ("edge-inp-gaasp-design.yaml", [(300.0, "300 K")]),
        ("edge-inp-gainp-design.yaml", [(230.0, "230 K"), (300.0, "300 K")]),
        ("nitride-cavity-pulse-design.yaml", [(300.0, "300 K")]),
    ]

    rows = []

    diff, direct, tj, gam = crosscheck_against_evaluate("edge-inp-gaasp-design.yaml")
    print(f"cross-check (gaasp, ibm, w=2 meV, Tj={tj:.9f}, gamma={gam:.9f}): "
          f"|direct - evaluate.eps_op| = {diff:.2e} (eps={direct!r})")

    # Process each card
    for card_filename, temps in cards_to_load:
        print(f"Loading {card_filename}...")
        card = load_card_copy(card_filename)

        # Get simplified card name
        card_name = card_filename.replace(".yaml", "")

        # Process each temperature
        for temp_K, temp_label in temps:
            print(f"  Processing {temp_label}...")

            # Process each lineshape
            for lineshape in ["lorentzian", "ibm"]:
                # Narrow-filter limit (w -> 0)
                row_narrow = compute_epsilon_row(
                    card_name, card, temp_K, lineshape, filter_width_meV=None
                )
                rows.append(row_narrow)
                print(
                    f"    {lineshape}: eps_narrow = {row_narrow['eps_narrow']:.4f}"
                )

                # At w = gamma_x (= Gamma linewidth)
                gam_x = row_narrow["gamma_x_meV"]
                row_at_w = compute_epsilon_row(
                    card_name, card, temp_K, lineshape, filter_width_meV=gam_x
                )
                rows.append(row_at_w)
                print(
                    f"    {lineshape}: eps(w=gamma_x={gam_x:.3f}) = {row_at_w['eps_at_w']:.4f}"
                )

                # Over a range of filter widths
                filter_widths = [0.5, 1.0, 2.0, 5.0, 10.0, 20.0]
                for w in filter_widths:
                    row_w = compute_epsilon_row(
                        card_name, card, temp_K, lineshape, filter_width_meV=w
                    )
                    rows.append(row_w)

    # Write CSV
    csv_path = out_dir / "eps_lineshape.csv"
    print(f"\nWriting results to {csv_path}...")

    with csv_path.open("w", newline="", encoding="utf-8") as f:
        fieldnames = [
            "card",
            "T_K",
            "lineshape",
            "filter_width_meV",
            "w_label",
            "gamma_x_meV",
            "gamma_xx_meV",
            "delta_xx_meV",
            "eps_narrow",
            "eps_at_w",
            "Z_T",
            "provenance",
        ]
        writer = csv.DictWriter(f, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Done! Wrote {len(rows)} rows to {csv_path}")

    # Print summary of Lorentzian rows vs. expected table
    print("\n=== LORENTZIAN ROWS (expected reference values) ===")
    print(
        "Reference table values:\n"
        "  GaAsP, 300 K: eps_narrow=0.594, eps_at_w=0.907\n"
        "  GaInP, 300 K: eps_narrow=0.383, eps_at_w=0.651\n"
        "  nitride, 300 K: eps_narrow=0.754, eps_at_w=0.897\n"
        "  GaInP, 230 K: eps_narrow=0.154, eps_at_w=0.237\n"
    )

    lorentzian_rows = [r for r in rows if r["lineshape"] == "lorentzian"]
    for row in lorentzian_rows:
        if np.isnan(row["filter_width_meV"]):  # Narrow filter limit
            print(
                f"{row['card']:30s} T={row['T_K']:3.0f} K: "
                f"eps_narrow={row['eps_narrow']:.4f}"
            )
        elif row["filter_width_meV"] == row["gamma_x_meV"]:  # At w=gamma_x
            print(
                f"{row['card']:30s} T={row['T_K']:3.0f} K: "
                f"eps(w=gamma_x)={row['eps_at_w']:.4f} (gamma_x={row['gamma_x_meV']:.3f} meV)"
            )


if __name__ == "__main__":
    main()
