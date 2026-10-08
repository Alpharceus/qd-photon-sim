"""numpy/NaN-safe JSON for the Studio API.

NaN and +-inf become null. When a dict carries a NaN value under key k and a
non-empty explanation under "k_invalid_reason" (the fsim_core convention,
e.g. g2_op / g2_op_invalid_reason), the serialized dict also gets
"k__nan_reason" so the UI can print the reason instead of "n/a"."""
from __future__ import annotations

import dataclasses
import json
import math
from pathlib import Path

import numpy as np


def _is_bad_float(v) -> bool:
    return isinstance(v, float) and not math.isfinite(v)


def to_jsonable(obj):
    if obj is None or isinstance(obj, (bool, str)):
        return obj
    if isinstance(obj, np.bool_):
        return bool(obj)
    if isinstance(obj, (int, np.integer)) and not isinstance(obj, bool):
        return int(obj)
    if isinstance(obj, (float, np.floating)):
        f = float(obj)
        return f if math.isfinite(f) else None
    if isinstance(obj, np.ndarray):
        return [to_jsonable(x) for x in obj.tolist()]
    if isinstance(obj, dict):
        out = {}
        for k, v in obj.items():
            out[str(k)] = to_jsonable(v)
        for k, v in obj.items():
            if isinstance(k, str) and (_is_bad_float(v) or
                                       (isinstance(v, np.floating) and not np.isfinite(v))):
                reason = obj.get(f"{k}_invalid_reason")
                if isinstance(reason, str) and reason:
                    out[f"{k}__nan_reason"] = reason
        return out
    if isinstance(obj, (list, tuple, set)):
        return [to_jsonable(x) for x in obj]
    if dataclasses.is_dataclass(obj) and not isinstance(obj, type):
        return to_jsonable(dataclasses.asdict(obj))
    if isinstance(obj, Path):
        return str(obj)
    return str(obj)


def dumps(obj) -> str:
    return json.dumps(to_jsonable(obj), allow_nan=False, separators=(",", ":"))


def canonical(obj) -> str:
    """Stable JSON text for hashing (sorted keys, NaN-safe)."""
    return json.dumps(to_jsonable(obj), allow_nan=False, sort_keys=True, separators=(",", ":"))


def nan_from_null(values) -> list:
    """Inverse helper for arrays coming back from JSON: null -> nan."""
    return [float("nan") if v is None else v for v in values]
