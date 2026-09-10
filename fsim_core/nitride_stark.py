"""Bias and Stark diagnostics for planar nitride traces.

These functions reduce rows made by an evaluator; they deliberately do not
solve optical states or invert a measured slope into a screening fraction.
The abrupt-depletion relations follow the depletion approximation [DR Sze &
Ng, Physics of Semiconductor Devices, 3rd ed. (2007)].  The Zhang reference
is a non-gating transcription [V Zhang et al., Appl. Phys. Lett. 108,
153102 (2016), Fig. 5].
"""
from __future__ import annotations

import math


ZHANG2016_SLOPE_MEV_PER_V = -10.0
ZHANG2016_COMPARISON = {
    "slope_meV_per_V": ZHANG2016_SLOPE_MEV_PER_V,
    "citation": "Zhang et al., Appl. Phys. Lett. 108, 153102 (2016), Fig. 5",
    "url": "https://arxiv.org/abs/1602.02325",
    "measurement": "bias-dependent PL below 2 V",
    "temperature": "unknown in supplied digest",
    "excitation": "unknown in supplied digest",
    "use": "non-gating comparison; not a planar-card fitted target",
}


def _finite(value):
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(value)


def resolve_bias(diode, *, T_j_K, current_uA=None, junction_voltage_V=None,
                 field_polarity=1, external_field_kVcm=0.0):
    """Resolve one forward-bias control into electrical and applied fields.

    Forward flat-band clipping is inherited unchanged from ``diode.depletion``.
    Negative junction bias, reverse leakage, and breakdown are intentionally
    outside this compact diagnostic's domain and are reported invalid.
    """
    reasons = []
    one_control = (current_uA is None) != (junction_voltage_V is None)
    if not one_control:
        reasons.append("supply exactly one of current_uA or junction_voltage_V")
    if not _finite(T_j_K) or T_j_K <= 0:
        reasons.append("T_j_K must be finite and positive")
    if isinstance(field_polarity, bool) or not isinstance(field_polarity, int) or field_polarity not in (-1, 1):
        reasons.append("field_polarity must be integer +1 or -1")
    if not _finite(external_field_kVcm):
        reasons.append("external_field_kVcm must be finite")
    if current_uA is not None and (not _finite(current_uA) or current_uA < 0):
        reasons.append("current_uA must be finite and nonnegative")
    if junction_voltage_V is not None and (not _finite(junction_voltage_V) or junction_voltage_V < 0):
        reasons.append("junction_voltage_V must be finite and nonnegative; reverse bias unsupported")
    result = {"T_j_K": T_j_K, "current_uA": current_uA, "V_j": junction_voltage_V,
              "V_terminal": float("nan"), "diode_field_kVcm": float("nan"),
              "applied_field_kVcm": float("nan"), "bias_valid": False,
              "invalid_reasons": reasons}
    if reasons:
        return result
    try:
        if current_uA is not None:
            terminal, vj = diode.v_of_i(current_uA * 1e-6, T_j_K)
            result["V_j"] = vj
            result["V_terminal"] = terminal
        else:
            vj = junction_voltage_V
            current_A = diode.j_of_vj(vj, T_j_K) * diode.area_cm2
            if not _finite(current_A):
                raise OverflowError("junction current overflow")
            result["current_uA"] = current_A * 1e6
            result["V_terminal"] = vj + current_A * diode.R_s_ohm
        if not all(_finite(result[key]) for key in ("V_j", "V_terminal", "current_uA")):
            raise OverflowError("unsupported operating point")
        dep = diode.depletion(result["V_j"], T_j_K)
        field = dep.F_kVcm
        if not _finite(field):
            raise OverflowError("depletion field is nonfinite")
        result["diode_field_kVcm"] = field
        result["applied_field_kVcm"] = field_polarity * field + external_field_kVcm
        result["bias_valid"] = True
    except (ArithmeticError, OverflowError, ValueError) as exc:
        result["invalid_reasons"].append("unsupported operating point: %s" % exc)
    return result


def _row_id(row, index):
    return row.get("row_id", row.get("id", index))


def _regime_marker(row):
    """Caller-declared discontinuities that must not be differentiated across."""
    return tuple((key, row.get(key)) for key in
                 ("depletion_regime", "depletion_clipped", "flat_band", "trace_branch") if key in row)


def _three_point_derivative(xs, ys, at):
    # Derivative of the Lagrange interpolating quadratic, valid on nonuniform x.
    total = 0.0
    for j in range(3):
        other = [k for k in range(3) if k != j]
        a, b = other
        total += ys[j] * ((2.0 * xs[at] - xs[a] - xs[b]) /
                           ((xs[j] - xs[a]) * (xs[j] - xs[b])))
    return total


def stark_derivatives(rows, *, voltage_key="V_j", energy_key="E_X_eV",
                      valid_key="spectroscopy_valid"):
    """Return same-order local dE/dV diagnostics using adjacent triplets.

    Policy: a valid interior point uses itself and its two adjacent physical
    samples; endpoints use the nearest three. Invalid samples and declared
    depletion-regime changes break the trace.  Callers must pre-partition
    other kinks and trace branches.
    """
    rows = list(rows)
    voltages = []
    for row in rows:
        v = row.get(voltage_key)
        if not _finite(v):
            raise ValueError("trace voltage must be finite")
        voltages.append(float(v))
    if any(voltages[i + 1] <= voltages[i] for i in range(len(voltages) - 1)):
        raise ValueError("trace voltages must be strictly increasing and distinct")
    out = []
    for i, row in enumerate(rows):
        rec = dict(row)
        rec.update(dE_X_dV_meV_per_V=float("nan"), derivative_valid=False,
                   derivative_row_ids=[] , derivative_policy="adjacent_three_point_nonuniform")
        valid = bool(row.get(valid_key, False)) and _finite(row.get(energy_key))
        if valid and len(rows) >= 3:
            indices = (0, 1, 2) if i == 0 else (len(rows)-3, len(rows)-2, len(rows)-1) if i == len(rows)-1 else (i-1, i, i+1)
            selected = [rows[j] for j in indices]
            if (all(bool(x.get(valid_key, False)) and _finite(x.get(energy_key)) for x in selected)
                    and len({_regime_marker(x) for x in selected}) == 1):
                xs = [voltages[j] for j in indices]
                ys = [float(rows[j][energy_key]) for j in indices]
                rec["dE_X_dV_meV_per_V"] = 1000.0 * _three_point_derivative(xs, ys, indices.index(i))
                rec["derivative_valid"] = True
                rec["derivative_row_ids"] = [_row_id(rows[j], j) for j in indices]
        out.append(rec)
    return out


_COORDINATE_KEYS = {"row_id", "id", "V_j", "V_j_V", "junction_voltage_V", "current_uA", "T_j_K",
                    "E_X_eV", "spectroscopy_valid", "spectroscopy_invalid_reasons", "valid",
                    "tau_rad_bare_ns", "tau_rad_cavity_ns", "dE_X_dV_meV_per_V", "derivative_valid",
                    "derivative_row_ids", "invalid_reasons", "screening", "screening_fraction"}
_TRACE_KEYS = ("geometry", "shape", "composition", "surrounding_well", "regime",
               "field_polarity", "applied_static_field_kVcm", "temperature_mode",
               "controlled_temperature", "T_hs", "cavity_reference_V_j_V")


def _freeze(value):
    if isinstance(value, dict): return tuple(sorted((k, _freeze(v)) for k, v in value.items()))
    if isinstance(value, (list, tuple)): return tuple(_freeze(v) for v in value)
    return value


def _group_key(row):
    # Explicit physics keys are retained when present; remaining fixed scalar
    # nuisance settings are included too, so unlike geometries never pool.
    keys = set(_TRACE_KEYS) | {k for k in row if k not in _COORDINATE_KEYS and k not in {"screening", "screening_fraction"}}
    return tuple(sorted((k, _freeze(row.get(k))) for k in keys if k in row))


def _screening(row):
    return row.get("screening_fraction", row.get("screening"))


def _fit(samples, voltage_key, energy_key):
    x = [float(r[voltage_key]) for r in samples]
    y = [1000.0 * float(r[energy_key]) for r in samples]
    xm, ym = sum(x)/len(x), sum(y)/len(y)
    den = sum((a-xm)**2 for a in x)
    slope = sum((a-xm)*(b-ym) for a, b in zip(x, y))/den
    intercept = ym-slope*xm
    rms = math.sqrt(sum((b-(intercept+slope*a))**2 for a,b in zip(x,y))/len(x))
    return slope, rms


def screening_compatibility(rows, *, slope_range_meV_per_V, voltage_window_V,
                            lifetime_range_ns=None):
    """Fit same-window slopes per supplied screening hypothesis.

    A result is conditional on all fixed nuisance settings in its group.  It
    reports compatibility only, never an interpolated or unique screening.
    """
    try: lo, hi = slope_range_meV_per_V
    except (TypeError, ValueError): raise ValueError("slope_range_meV_per_V must have two values")
    try: vlo, vhi = voltage_window_V
    except (TypeError, ValueError): raise ValueError("voltage_window_V must have two values")
    if not all(_finite(x) for x in (lo, hi, vlo, vhi)) or lo > hi or vlo >= vhi:
        raise ValueError("invalid slope range or voltage window")
    if lifetime_range_ns is not None:
        required = {"voltage_V", "min_ns", "max_ns", "kind"}
        if set(lifetime_range_ns) != required or lifetime_range_ns["kind"] not in ("bare", "cavity"):
            raise ValueError("lifetime_range_ns must exactly specify voltage_V, min_ns, max_ns, kind")
        if not all(_finite(lifetime_range_ns[k]) and lifetime_range_ns[k] > 0 for k in ("voltage_V", "min_ns", "max_ns")) or lifetime_range_ns["min_ns"] > lifetime_range_ns["max_ns"]:
            raise ValueError("invalid lifetime range")
    grouped = {}
    for index, row in enumerate(rows):
        grouped.setdefault((_group_key(row), _screening(row)), []).append((index, row))
    records = []
    for (gkey, hypothesis), numbered in grouped.items():
        # Preserve physical sequence; split on invalid gaps and voltage turns.
        branches, current, direction = [], [], 0
        for index, row in numbered:
            v = row.get("V_j", row.get("V_j_V", row.get("junction_voltage_V")))
            good = bool(row.get("spectroscopy_valid", row.get("valid", False))) and _finite(v) and _finite(row.get("E_X_eV"))
            if not good:
                if current: branches.append(current); current = []; direction = 0
                continue
            if current:
                delta = float(v) - float(current[-1][1].get("V_j", current[-1][1].get("V_j_V", current[-1][1].get("junction_voltage_V"))))
                newdir = 1 if delta > 0 else -1 if delta < 0 else 0
                if newdir == 0 or (direction and newdir != direction):
                    if current: branches.append(current)
                    current = [(index, row)]; direction = 0
                    continue
                direction = newdir
            current.append((index, row))
        if current: branches.append(current)
        made = False
        for branch_id, branch in enumerate(branches):
            window = [(i, r) for i, r in branch if vlo <= float(r.get("V_j", r.get("V_j_V", r.get("junction_voltage_V")))) <= vhi]
            if len(window) < 3:
                continue
            vkey = "V_j" if "V_j" in window[0][1] else "V_j_V" if "V_j_V" in window[0][1] else "junction_voltage_V"
            slope, rms = _fit([r for _, r in window], vkey, "E_X_eV")
            lifetime_ok = True
            lifetime_value = None
            if lifetime_range_ns is not None:
                key = "tau_rad_bare_ns" if lifetime_range_ns["kind"] == "bare" else "tau_rad_cavity_ns"
                matches = [r for _, r in branch if float(r.get(vkey)) == float(lifetime_range_ns["voltage_V"])]
                lifetime_value = matches[0].get(key) if matches else None
                lifetime_ok = (_finite(lifetime_value) and lifetime_range_ns["min_ns"] <= lifetime_value <= lifetime_range_ns["max_ns"])
            # A tiny arithmetic allowance preserves inclusive measurement
            # endpoints when an exactly linear fixture is accumulated in
            # binary floating point; it is not a model uncertainty.
            endpoint_eps = 1e-12 * max(1.0, abs(lo), abs(hi), abs(slope))
            compatible = lo - endpoint_eps <= slope <= hi + endpoint_eps and lifetime_ok
            records.append({"group": gkey, "screening": hypothesis, "branch_id": branch_id,
                            "voltage_window_V": (vlo, vhi), "fitted_slope_meV_per_V": slope,
                            "slope_interval_meV_per_V": (slope, slope), "residual_rms_meV": rms,
                            "nonlinearity_diagnostic": "rms_residual_meV", "row_ids": [_row_id(r, i) for i,r in window],
                            "compatible": compatible, "lifetime_value_ns": lifetime_value,
                            "identification_status": "compatible" if compatible else "incompatible"})
            made = True
        if not made:
            records.append({"group": gkey, "screening": hypothesis, "branch_id": None,
                            "voltage_window_V": (vlo, vhi), "fitted_slope_meV_per_V": float("nan"),
                            "slope_interval_meV_per_V": None, "residual_rms_meV": float("nan"),
                            "nonlinearity_diagnostic": "incomplete_model_coverage", "row_ids": [],
                            "compatible": False, "lifetime_value_ns": None,
                            "identification_status": "incomplete_model_coverage"})
    # A shared nonpolar curve cannot identify screening.  Equality is tested
    # from the explicitly fitted same-window predictions, never arbitrary order.
    by_group = {}
    for rec in records: by_group.setdefault((rec["group"], rec["branch_id"]), []).append(rec)
    for members in by_group.values():
        usable = [r for r in members if _finite(r["fitted_slope_meV_per_V"])]
        if len(usable) >= 2 and len({round(r["fitted_slope_meV_per_V"], 12) for r in usable}) == 1:
            shared = usable[0]["compatible"]
            for rec in usable:
                rec["identification_status"] = "screening_unidentifiable"
                rec["compatible"] = shared
    return records
