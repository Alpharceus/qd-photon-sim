"""Composed device spec (schema "fsim-composed/1", design 1.3) and its validator.

`ComposedSpec.from_dict` is structural validation only (required keys, schema
string, solver choice, block shapes). `validate_spec` checks the spec against
the registry: kinds exist, params are known, leaves resolve (units, tags,
choices) and requires/provides are compatible.
"""
from __future__ import annotations

import copy
from dataclasses import dataclass, replace
from typing import Optional

from . import registry as _registry
from .params import ParamSpec, resolve_leaf

SCHEMA = "fsim-composed/1"
SOLVERS = ("auto", "rate", "lindblad", "lindblad_incoherent", "polaron", "jc")
AXIS_KEYS = ("material", "structure", "drive", "setup")
_TOP_KEYS = ("schema", "T_K") + AXIS_KEYS + ("solver",)
T_K_SPEC = ParamSpec("T_K", "K", 4.2, lo=0.1, hi=400.0, tag="A",
                     source="operating temperature of the cryostat / heat sink",
                     help="Sample temperature. A range is an envelope, not an average.")


def _check_block(axis: str, blk) -> dict:
    if not isinstance(blk, dict):
        raise ValueError(f"composed.{axis} must be a mapping with kind/params, got {type(blk).__name__}")
    extra = set(blk) - {"kind", "params"}
    if extra:
        raise ValueError(f"composed.{axis}: unknown keys {sorted(extra)} (allowed: kind, params)")
    if not isinstance(blk.get("kind"), str) or not blk["kind"]:
        raise ValueError(f"composed.{axis}.kind must be a non-empty string")
    params = blk.get("params", {})
    if params is None:
        params = {}
    if not isinstance(params, dict) or not all(isinstance(k, str) for k in params):
        raise ValueError(f"composed.{axis}.params must be a mapping of parameter name -> leaf")
    return {"kind": blk["kind"], "params": copy.deepcopy(params)}


@dataclass(frozen=True)
class ComposedSpec:
    T_K: object                      # raw leaf: scalar or {value|range, unit, tag, source|where}
    material: dict
    structure: dict
    drive: dict
    setup: dict
    solver: str = "auto"
    schema: str = SCHEMA

    @staticmethod
    def from_dict(d) -> "ComposedSpec":
        if not isinstance(d, dict):
            raise ValueError(f"composed spec must be a mapping, got {type(d).__name__}")
        unknown = set(d) - set(_TOP_KEYS)
        if unknown:
            raise ValueError(f"unknown top-level key(s) in composed spec: {sorted(unknown)}; "
                             f"allowed: {list(_TOP_KEYS)}")
        missing = [k for k in ("schema", "T_K") + AXIS_KEYS if k not in d]
        if missing:
            raise ValueError(f"composed spec missing required key(s): {missing}")
        if d["schema"] != SCHEMA:
            raise ValueError(f"composed spec schema must be {SCHEMA!r}, got {d['schema']!r}")
        solver = d.get("solver", "auto")
        if solver not in SOLVERS:
            raise ValueError(f"composed.solver must be one of {list(SOLVERS)}, got {solver!r}")
        return ComposedSpec(T_K=copy.deepcopy(d["T_K"]),
                            solver=solver, schema=d["schema"],
                            **{a: _check_block(a, d[a]) for a in AXIS_KEYS})

    def to_dict(self) -> dict:
        out = {"schema": self.schema, "T_K": copy.deepcopy(self.T_K)}
        for a in AXIS_KEYS:
            blk = getattr(self, a)
            out[a] = {"kind": blk["kind"], "params": copy.deepcopy(blk["params"])}
        out["solver"] = self.solver
        return out


def _issue(level: str, path: str, message: str) -> dict:
    return {"level": level, "path": path, "message": message}


def _resolved_entry(leaf) -> dict:
    r = {"unit": leaf.unit, "tag": leaf.tag}
    if leaf.kind == "value":
        r["value"] = leaf.value
    elif leaf.kind == "range":
        r["range"] = list(leaf.range)
        r["scale"] = leaf.scale
    else:
        r["choice"] = list(leaf.choice)
    if leaf.defaulted:
        r["default"] = True
    return r


def _central_vals(entry, blk, path_prefix):
    from .params import central
    out = {}
    for p in entry.params:
        leaf = resolve_leaf(blk["params"][p.name] if p.name in blk["params"] else p.default, p, f"{path_prefix}.{p.name}")
        out[p.name] = central(leaf)
    return out


def _pair_incompatibility(spec, entries):
    """Reason string if structures.incompatibility() declares the (material, structure) pair unbuildable, else None.
    Never raises (a failing build is reported by the other validation paths)."""
    try:
        from . import structures as _st
        from .params import central
        T = float(central(resolve_leaf(spec.T_K, T_K_SPEC, "T_K")))
        em = entries["material"].build(_central_vals(entries["material"], spec.material, "material"), T)
        return _st.incompatibility(entries["structure"].kind, em, _central_vals(entries["structure"], spec.structure, "structure"))
    except Exception:                                   # noqa: BLE001
        return None


def validate_spec(spec, registry: Optional[_registry.Registry] = None) -> dict:
    """-> {ok, solver, regime, issues:[{level, path, message}], resolved:{path:{value|range, unit, tag}}}.
    ok is False iff any issue has level 'error'. `regime` is only known after
    rates.combine (P5), so it is "unknown" here."""
    issues: list = []
    resolved: dict = {}
    if isinstance(spec, dict):
        try:
            spec = ComposedSpec.from_dict(spec)
        except ValueError as e:
            return {"ok": False, "solver": None, "regime": "unknown",
                    "issues": [_issue("error", "", str(e))], "resolved": {}}
    reg = _registry._reg(registry)

    try:
        resolved["T_K"] = _resolved_entry(resolve_leaf(spec.T_K, T_K_SPEC, "T_K"))
    except ValueError as e:
        issues.append(_issue("error", "T_K", str(e)))

    entries = {}
    for axis in AXIS_KEYS:
        blk = getattr(spec, axis)
        kinds = reg.kinds(axis)
        if not kinds:
            issues.append(_issue("error", f"{axis}.kind", f"no {axis} registered yet"))
            continue
        if blk["kind"] not in kinds:
            issues.append(_issue("error", f"{axis}.kind",
                                 f"unknown {axis} kind {blk['kind']!r}; registered: {kinds}"))
            continue
        entry = reg.get(axis, blk["kind"])
        entries[axis] = entry
        if entry.status == "needs_physics":
            issues.append(_issue("error", f"{axis}.kind",
                                 f"{axis} kind {entry.kind!r} is registered but awaiting physics"))
        known = {p.name for p in entry.params}
        for name in blk["params"]:
            if name not in known:
                issues.append(_issue("error", f"{axis}.params.{name}",
                                     f"unknown param {name!r} for {axis} kind {entry.kind!r}; "
                                     f"known: {sorted(known)}"))
        for p in entry.params:
            path = f"{axis}.{p.name}"
            try:
                if p.name in blk["params"]:
                    leaf = resolve_leaf(blk["params"][p.name], p, path)
                else:
                    leaf = resolve_leaf(p.default, p, path)
                    leaf = replace(leaf, defaulted=True)
            except ValueError as e:
                issues.append(_issue("error", f"{axis}.params.{p.name}", str(e)))
                continue
            resolved[path] = _resolved_entry(leaf)

    # requires / provides: every requirement of an entry must be provided by
    # some OTHER chosen entry.
    for axis, entry in entries.items():
        if not entry.requires:
            continue
        have = set()
        for other_axis, other in entries.items():
            if other_axis != axis:
                have |= set(other.provides)
        lacking = sorted(set(entry.requires) - have)
        if lacking:
            issues.append(_issue("error", f"{axis}.kind",
                                 f"{axis} kind {entry.kind!r} requires {lacking}, not provided by "
                                 f"the other chosen entries "
                                 f"({ {a: sorted(e.provides) for a, e in entries.items() if a != axis} })"))

    # FIX-1b/1a: declared (material, structure) incompatibility (structures.STRUCTURE_COMPAT); a WARNING so the engine still
    # reports valid=False with the reason (engine._validity) and the Composer picker greys the pair
    if "material" in entries and "structure" in entries:
        why = _pair_incompatibility(spec, entries)
        if why:
            issues.append(_issue("warning", "structure.kind", f"incompatible with the chosen material: {why}"))

    solvers = reg.kinds("solver")
    if spec.solver == "auto":
        if not solvers:
            issues.append(_issue("warning", "solver", "no solver registered yet"))
    elif not solvers:
        issues.append(_issue("warning", "solver", "no solver registered yet"))
    elif spec.solver not in solvers:
        issues.append(_issue("error", "solver",
                             f"unknown solver {spec.solver!r}; registered: {solvers}"))

    ok = not any(i["level"] == "error" for i in issues)
    return {"ok": ok, "solver": spec.solver, "regime": "unknown", "issues": issues,
            "resolved": resolved}
