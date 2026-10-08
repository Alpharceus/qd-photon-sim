"""Provenance-tag bookkeeping (no physics): parse [V]/[DR]/[E]/[A] tags out
of fsim_core's tag_chain / provenance payloads and pick the widest.

Order (narrow -> wide): V < DR < E < A. The widest tag of a chain is the
least certain link, which is what a readout must show (CONTRACT rule 2).
"""
from __future__ import annotations

import re

ORDER = ("V", "DR", "E", "A")
_RANK = {t: i for i, t in enumerate(ORDER)}
_BRACKET = re.compile(r"\[([^\[\]]{1,80})\]")
_SPLIT = re.compile(r"[/,;\s]+")

SRC_CORE = "fsim_core tag_chain"
SRC_PARSED = "parsed from provenance"
SRC_ASSUMED = "no tag returned; assumed widest"


def widest(tags):
    tags = [t for t in tags if t in _RANK]
    return max(tags, key=_RANK.get) if tags else None


def tags_in_text(text: str) -> set:
    """Every V/DR/E/A token inside [...] brackets of a string."""
    out = set()
    for m in _BRACKET.finditer(str(text)):
        for tok in _SPLIT.split(m.group(1).strip()):
            if tok.upper() in _RANK and tok == tok.upper():
                out.add(tok)
    return out


def tags_in(obj) -> set:
    """Walk a provenance payload (str / dict / list, nested): bracket tags in
    every string plus the value of every "tag" key."""
    out = set()
    if obj is None:
        return out
    if isinstance(obj, str):
        return tags_in_text(obj)
    if isinstance(obj, dict):
        for k, v in obj.items():
            if k == "tag" and isinstance(v, str):
                t = v.strip("[] ").upper()
                if t in _RANK:
                    out.add(t)
                else:
                    out |= tags_in_text(v)
            else:
                out |= tags_in(v)
        return out
    if isinstance(obj, (list, tuple)):
        for v in obj:
            out |= tags_in(v)
        return out
    return out


def chain_of(scalars: dict) -> tuple:
    """(tag, source) for a run's scalars: fsim_core's tag_chain when it
    returns one, else the widest tag parsed from scalars["provenance"], else
    "A" with the source saying it was assumed (interface I1)."""
    tc = (scalars or {}).get("tag_chain")
    if isinstance(tc, str) and tc.strip():
        t = widest(tags_in_text(tc) | ({tc.strip("[] ").upper()} & set(ORDER)))
        if t:
            return t, SRC_CORE
    t = widest(tags_in((scalars or {}).get("provenance")))
    if t:
        return t, SRC_PARSED
    return "A", SRC_ASSUMED
