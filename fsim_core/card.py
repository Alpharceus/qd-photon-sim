"""Parameter-card schema (FSIM planning doc, section 3).

Every card entry = {value | range, unit, tag in [V/DR/E/A], source}.
Ranges are swept or used as fit bounds -- never averaged (council rule).
Derived quantities inherit the *widest* tag in their input chain.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from enum import IntEnum
from pathlib import Path
from typing import Optional

import yaml


class Tag(IntEnum):
    """Provenance tags, ordered from strongest to widest (weakest)."""

    V = 0    # verified / published measurement on this device
    DR = 1   # derived from published data (digitized, re-fit, restated)
    E = 2    # estimated from comparable systems or literature
    A = 3    # assumed; no measurement backs this number

    @property
    def label(self) -> str:
        return f"[{self.name}]"


def widest(*tags: Tag) -> Tag:
    """Tag propagation rule: a derived number is only as good as its weakest input."""
    return Tag(max(int(t) for t in tags))


@dataclass
class Param:
    name: str
    unit: str
    tag: Tag
    source: str
    value: Optional[float] = None
    range: Optional[tuple[float, float]] = None

    def __post_init__(self):
        if (self.value is None) == (self.range is None):
            raise ValueError(f"{self.name}: exactly one of value/range required")
        if self.range is not None and self.range[0] > self.range[1]:
            raise ValueError(f"{self.name}: range must be (lo, hi)")

    @property
    def is_free(self) -> bool:
        """Range-valued params are sweep/fit intervals, never point estimates."""
        return self.range is not None

    @property
    def fixed(self) -> float:
        if self.value is None:
            raise ValueError(
                f"{self.name} is range-valued; ranges are swept or fit, never averaged"
            )
        return self.value

    @property
    def bounds(self) -> tuple[float, float]:
        return self.range if self.range is not None else (self.value, self.value)


@dataclass
class DataSet:
    name: str
    tag: Tag
    source: str
    columns: list[str]
    rows: list[dict]


@dataclass
class Card:
    name: str
    meta: dict
    params: dict[str, Param]
    datasets: dict[str, DataSet] = field(default_factory=dict)
    path: Optional[Path] = None

    def __getitem__(self, key: str) -> Param:
        return self.params[key]

    def tag_chain(self, *names: str) -> Tag:
        return widest(*(self.params[n].tag for n in names))

    def placeholders(self) -> list[str]:
        """Entries whose source marks them as placeholder -- these block any V-a claim."""
        out = [n for n, p in self.params.items() if "PLACEHOLDER" in p.source.upper()]
        out += [n for n, d in self.datasets.items() if "PLACEHOLDER" in d.source.upper()]
        return out


def _parse_param(name: str, raw: dict) -> Param:
    for key in ("unit", "tag", "source"):
        if key not in raw:
            raise ValueError(f"param {name}: missing required field '{key}'")
    rng = raw.get("range")
    return Param(
        name=name,
        unit=str(raw["unit"]),
        tag=Tag[raw["tag"]],
        source=str(raw["source"]),
        value=None if "value" not in raw else float(raw["value"]),
        range=None if rng is None else (float(rng[0]), float(rng[1])),
    )


def load_card(path: str | Path) -> Card:
    path = Path(path)
    doc = yaml.safe_load(path.read_text(encoding="utf-8"))
    params = {n: _parse_param(n, raw) for n, raw in doc.get("params", {}).items()}
    datasets = {}
    for n, raw in doc.get("data", {}).items():
        rows = [dict(r) for r in raw["points"]]
        cols = sorted({k for r in rows for k in r})
        datasets[n] = DataSet(
            name=n, tag=Tag[raw["tag"]], source=str(raw["source"]),
            columns=cols, rows=rows,
        )
    return Card(
        name=doc["meta"]["name"], meta=doc["meta"],
        params=params, datasets=datasets, path=path,
    )
