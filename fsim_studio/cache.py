"""Disk cache .cache/studio/<sha256>.json (gitignored). The key covers the
canonical request (design, mode, ranged, T grid) plus the fsim_core source
hash, so a physics change invalidates every cached result automatically."""
from __future__ import annotations

import hashlib
import json
import os
import threading

from . import CACHE_DIR, fsim_core_hash
from .serialize import canonical, dumps

_LOCK = threading.Lock()
# Version of the Studio's result shaping (tag_chain normalisation, labels,
# ...). Part of every key, so a change to how results are shaped invalidates
# cached results even when fsim_core did not change (H14).
STUDIO_SCHEMA = "fsim-studio-result/2"


def core_hash() -> str:
    """Indirection so tests can monkeypatch the source hash."""
    return fsim_core_hash()


def make_key(payload: dict) -> str:
    text = canonical({"payload": payload, "fsim_core_hash": core_hash(), "studio_schema": STUDIO_SCHEMA})
    return hashlib.sha256(text.encode("utf-8")).hexdigest()


def path_for(key: str):
    return CACHE_DIR / f"{key}.json"


def get(key: str):
    p = path_for(key)
    if not p.is_file():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def put(key: str, record: dict) -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    p = path_for(key)
    tmp = p.with_suffix(f".{os.getpid()}.{threading.get_ident()}.tmp")
    with _LOCK:
        tmp.write_text(dumps(record), encoding="utf-8")
        os.replace(tmp, p)


def find_run(run_id: str):
    """Cached record for a run_id ("FS-<card>-<hex6>"): the 6-hex tail is the
    cache-key prefix; the record's own run_id must match exactly."""
    prefix = run_id.rsplit("-", 1)[-1]
    if not CACHE_DIR.is_dir() or len(prefix) != 6 or not all(c in "0123456789abcdef" for c in prefix):
        return None
    hits = sorted(CACHE_DIR.glob(f"{prefix}*.json"), key=lambda p: p.stat().st_mtime, reverse=True)
    for p in hits:
        try:
            rec = json.loads(p.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        if rec.get("run_id") == run_id:
            return rec
    return None
