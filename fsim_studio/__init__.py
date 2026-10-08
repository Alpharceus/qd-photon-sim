"""FSIM Studio: local web front end over fsim_core.

Three-layer rule: this package computes no physics. It loads cards, calls
fsim_core.device.evaluate / evaluate_envelope, reads committed out/ files and
serves the results as JSON to the static front end in fsim_studio/web/.
"""
from __future__ import annotations

import hashlib
import os
import threading
from pathlib import Path

__version__ = "0.1.0"

ROOT = Path(__file__).resolve().parents[1]
PKG = Path(__file__).resolve().parent
WEB = PKG / "web"
CARDS = ROOT / "cards"
OUT = ROOT / "out"
# FSIM_STUDIO_CACHE overrides the cache location (verify scripts use a temp dir)
CACHE_DIR = Path(os.environ.get("FSIM_STUDIO_CACHE") or (ROOT / ".cache" / "studio"))

_HASH_LOCK = threading.Lock()
_HASH_MEMO = {"sig": None, "hash": None}
_SEP = bytes(1)  # a single NUL separator between name and contents


def _core_signature() -> tuple:
    return tuple((p.name, st.st_mtime_ns, st.st_size)
                 for p in sorted((ROOT / "fsim_core").glob("*.py")) for st in (p.stat(),))


def fsim_core_hash() -> str:
    """sha256 over every fsim_core/*.py source (sorted by name). Part of every
    cache key, so a physics change invalidates cached results automatically.
    Memoized on the (name, mtime, size) signature of the sources, so an edit
    while the server runs produces the new hash (H14)."""
    sig = _core_signature()
    with _HASH_LOCK:
        if _HASH_MEMO["sig"] == sig:
            return _HASH_MEMO["hash"]
    h = hashlib.sha256()
    for p in sorted((ROOT / "fsim_core").glob("*.py")):
        h.update(p.name.encode("utf-8"))
        h.update(_SEP)
        h.update(p.read_bytes())
        h.update(_SEP)
    digest = h.hexdigest()
    with _HASH_LOCK:
        _HASH_MEMO["sig"], _HASH_MEMO["hash"] = sig, digest
    return digest
