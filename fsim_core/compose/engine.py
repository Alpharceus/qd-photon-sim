"""Composed-engine entry points. STUB until P5 (rates, solvers, observables).

Validation already works; a spec that passes it raises NotImplementedError.
"""
from __future__ import annotations

from .spec import ComposedSpec, validate_spec


def evaluate_spec(spec, T_grid=None) -> dict:
    """spec: ComposedSpec or its dict form. Returns the legacy contract
    {'curves', 'scalars', 'series'} once P5 lands."""
    if isinstance(spec, dict):
        spec = ComposedSpec.from_dict(spec)
    res = validate_spec(spec)
    if not res["ok"]:
        lines = [f"{i['level']} {i['path']}: {i['message']}" for i in res["issues"]]
        raise ValueError("invalid composed spec:\n" + "\n".join(lines))
    raise NotImplementedError("compose engine: P5")


def evaluate_composed(design, T_grid=None) -> dict:
    """device.evaluate() arm for platform == 'composed'."""
    return evaluate_spec(ComposedSpec.from_dict(design.composed), T_grid)
