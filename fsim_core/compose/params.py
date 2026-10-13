"""Parameter specs and provenance-tagged leaves for the composable engine.

A ParamSpec is the self-description of one registry parameter (design 2.1).
A Leaf is one resolved spec value: a value, a range (+ scale) or a discrete
choice list, with unit, tag in {V, DR, E, A} and source / where.

Ranges are never averaged (card rule): point runs use `central`, labelled by
the caller "midpoint, not a prediction"; envelopes use `corners` / `lhs`.
Units are fixed per parameter and never converted.
"""
from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from typing import Any, Optional

import numpy as np

from ..card import Tag, widest

TAGS = ("V", "DR", "E", "A")
_SCALES = ("lin", "log")
_GROUPS = ("main", "advanced")
_LEAF_KEYS = {"value", "range", "choice", "scale", "unit", "tag", "source", "where"}


@dataclass(frozen=True)
class ParamSpec:
    name: str
    unit: str
    default: Any
    lo: Optional[float] = None       # physical band: UI warning only, never clamps (as META)
    hi: Optional[float] = None
    scale: str = "lin"
    choices: Optional[tuple] = None
    tag: str = "A"                   # provenance of the DEFAULT
    source: str = ""
    help: str = ""
    group: str = "main"              # main | advanced
    visible_if: Optional[dict] = None
    sweepable: bool = True
    tex: Optional[str] = None        # display symbol in KaTeX TeX (UI only; never used by the physics)

    def __post_init__(self):
        if self.tag not in TAGS:
            raise ValueError(f"{self.name}: tag must be one of {TAGS}, got {self.tag!r}")
        if self.scale not in _SCALES:
            raise ValueError(f"{self.name}: scale must be lin|log, got {self.scale!r}")
        if self.group not in _GROUPS:
            raise ValueError(f"{self.name}: group must be main|advanced, got {self.group!r}")
        if self.choices is not None:
            object.__setattr__(self, "choices", tuple(self.choices))

    def to_meta(self) -> dict:
        """Same dict shape as an entry of fsim_core.design_meta.META:
        {unit, tag, source} plus lo/hi (numeric params) and/or choices, each
        present only when set (META omits them rather than writing null)."""
        m = {"unit": self.unit, "tag": self.tag, "source": self.source}
        if self.lo is not None:
            m["lo"] = self.lo
        if self.hi is not None:
            m["hi"] = self.hi
        if self.choices is not None:
            m["choices"] = list(self.choices)
        if self.tex is not None:
            m["tex"] = self.tex
        return m


@dataclass(frozen=True)
class Leaf:
    """kind is 'value' | 'range' | 'choice'. `range` is (lo, hi); `choice` is
    a tuple of alternatives (a discrete scan; central = first)."""
    path: str
    kind: str
    unit: str
    tag: str
    source: str = ""
    where: str = ""
    value: Any = None
    range: Optional[tuple] = None
    scale: str = "lin"
    choice: Optional[tuple] = None
    defaulted: bool = False          # True when the spec omitted it and the registry default was used

    @property
    def is_free(self) -> bool:
        return self.kind != "value"


def widest_tag(*tags) -> str:
    """Order V < DR < E < A (fsim_core.card.Tag); returns the widest as a string."""
    if not tags:
        raise ValueError("widest_tag needs at least one tag")
    for t in tags:
        if str(t) not in TAGS:
            raise ValueError(f"unknown tag {t!r}; expected one of {TAGS}")
    return widest(*(Tag[str(t)] for t in tags)).name


def _is_num(x) -> bool:
    return isinstance(x, (int, float)) and not isinstance(x, bool) and math.isfinite(x)


def _check_scalar(v, spec: ParamSpec, path: str):
    d = spec.default
    if isinstance(d, bool):
        ok = isinstance(v, bool)
    elif isinstance(d, str):
        ok = isinstance(v, str)
    else:
        ok = _is_num(v)
    if not ok:
        raise ValueError(f"{path}: expected {type(d).__name__} for {spec.name}, got {v!r}")
    if spec.choices is not None and v not in spec.choices:
        raise ValueError(f"{path}: {v!r} not in choices {list(spec.choices)}")


def resolve_leaf(raw, spec: ParamSpec, path: Optional[str] = None) -> Leaf:
    """Design 1.3 leaf rules. bare scalar -> value with the spec's default
    tag/source; mapping -> exactly one of value|range|choice plus unit (must
    equal spec.unit), tag in V|DR|E|A, source and/or where."""
    path = path or spec.name
    if not isinstance(raw, dict):
        if isinstance(raw, (list, tuple)):
            raise ValueError(f"{path}: a bare list is not a leaf; use a mapping with range: or choice:")
        _check_scalar(raw, spec, path)
        return Leaf(path=path, kind="value", unit=spec.unit, tag=spec.tag,
                    source=spec.source, value=raw)
    extra = set(raw) - _LEAF_KEYS
    if extra:
        raise ValueError(f"{path}: unknown leaf keys {sorted(extra)}")
    present = [k for k in ("value", "range", "choice") if k in raw]
    if len(present) != 1:
        raise ValueError(f"{path}: leaf must hold exactly one of value/range/choice, got {present or 'none'}")
    if "unit" in raw:
        if raw["unit"] != spec.unit:
            raise ValueError(f"{path}: unit {raw['unit']!r} given, {spec.unit!r} expected "
                             "(no unit conversion)")
    elif spec.unit not in ("", "enum"):
        raise ValueError(f"{path}: unit missing, {spec.unit!r} expected")
    tag = raw.get("tag", spec.tag)
    if tag not in TAGS:
        raise ValueError(f"{path}: tag {tag!r} not in {TAGS}")
    source = raw.get("source", spec.source if "where" not in raw else "")
    where = raw.get("where", "")
    scale = raw.get("scale", spec.scale)
    if scale not in _SCALES:
        raise ValueError(f"{path}: scale must be lin|log, got {scale!r}")
    kw = dict(path=path, unit=spec.unit, tag=tag, source=source, where=where, scale=scale)
    kind = present[0]
    if kind == "value":
        _check_scalar(raw["value"], spec, path)
        return Leaf(kind="value", value=raw["value"], **kw)
    if kind == "range":
        r = raw["range"]
        if (not isinstance(r, (list, tuple)) or len(r) != 2 or not all(_is_num(x) for x in r)
                or isinstance(spec.default, (str, bool))):
            raise ValueError(f"{path}: range must be [lo, hi] of two numbers for a numeric parameter")
        if r[0] > r[1]:
            raise ValueError(f"{path}: range must be [lo, hi] with lo <= hi, got {list(r)}")
        if scale == "log" and r[0] <= 0:
            raise ValueError(f"{path}: log-scale range needs lo > 0, got {list(r)}")
        return Leaf(kind="range", range=(float(r[0]), float(r[1])), **kw)
    c = raw["choice"]
    if not isinstance(c, (list, tuple)) or len(c) < 1:
        raise ValueError(f"{path}: choice must be a non-empty list of alternatives")
    for v in c:
        _check_scalar(v, spec, path)
    return Leaf(kind="choice", choice=tuple(c), **kw)


def central(leaf: Leaf):
    """Declared central value: value; lin midpoint / log geometric mean of a
    range; first alternative of a choice. Not a prediction for ranges."""
    if leaf.kind == "value":
        return leaf.value
    if leaf.kind == "range":
        lo, hi = leaf.range
        return math.sqrt(lo * hi) if leaf.scale == "log" else 0.5 * (lo + hi)
    return leaf.choice[0]


def _items(leaves):
    return list(leaves.items()) if isinstance(leaves, dict) else [(l.path, l) for l in leaves]


def corners(leaves, max_n: int = 16) -> list:
    """Corner points of the range box, as a list of {path: value} dicts.
    Range leaves take lo/hi, every other leaf its central value. Always starts
    with all-lo and all-hi. If 2^k > max_n the set is truncated
    deterministically to: all-lo, all-hi, then single flips from each."""
    items = _items(leaves)
    free = [(p, l) for p, l in items if l.kind == "range"]
    base = {p: central(l) for p, l in items}
    if not free:
        return [dict(base)]
    k = len(free)
    lo = {p: l.range[0] for p, l in free}
    hi = {p: l.range[1] for p, l in free}
    all_lo, all_hi = {**base, **lo}, {**base, **hi}
    if 2 ** k <= max_n:
        pts = [all_lo, all_hi]
        for bits in itertools.product((0, 1), repeat=k):
            pt = {**base, **{p: (hi if b else lo)[p] for (p, _), b in zip(free, bits)}}
            if pt != all_lo and pt != all_hi:
                pts.append(pt)
        return pts[:max_n]
    pts = [all_lo, all_hi]
    for p, _ in free:
        pts.append({**base, **lo, p: hi[p]})
        pts.append({**base, **hi, p: lo[p]})
    return pts[:max_n]


def lhs(leaves, n: int = 32, seed: int = 0) -> list:
    """Latin-hypercube samples over the range leaves (log leaves sampled in
    log space), list of {path: value}; choice leaves are drawn uniformly,
    value leaves held fixed."""
    items = _items(leaves)
    rng = np.random.default_rng(seed)
    out = [{p: central(l) for p, l in items} for _ in range(n)]
    for p, l in items:
        if l.kind == "range":
            u = (rng.permutation(n) + rng.random(n)) / n
            lo, hi = l.range
            if l.scale == "log":
                vals = np.exp(math.log(lo) + u * (math.log(hi) - math.log(lo)))
            else:
                vals = lo + u * (hi - lo)
            for i in range(n):
                out[i][p] = float(vals[i])
        elif l.kind == "choice":
            idx = rng.integers(0, len(l.choice), size=n)
            for i in range(n):
                out[i][p] = l.choice[int(idx[i])]
    return out
