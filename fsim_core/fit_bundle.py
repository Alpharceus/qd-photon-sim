"""Single resolver for the V-a joint fit bundle (phase0 fit_params.json).

out/ is gitignored, so a fresh clone has no out/phase0/fit_params.json.  The
tracked byte copy fsim_core/data/phase0_fit_params.json is used in that case.
Every reader of the bundle goes through this module so the two sources cannot
drift apart in code.
"""
import json
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
_OUT_PATH = _ROOT / "out" / "phase0" / "fit_params.json"
_SNAPSHOT_PATH = _ROOT / "fsim_core" / "data" / "phase0_fit_params.json"


def fit_params_path() -> Path:
    """The generated bundle if it exists, else the tracked snapshot."""
    return _OUT_PATH if _OUT_PATH.is_file() else _SNAPSHOT_PATH


def fit_params_source() -> str:
    """Provenance label matching fit_params_path()."""
    if fit_params_path() == _OUT_PATH:
        return "out/phase0/fit_params.json"
    return "fsim_core/data/phase0_fit_params.json (tracked snapshot)"


def load_fit_params() -> dict:
    """Full bundle JSON (keys include "params")."""
    return json.loads(fit_params_path().read_text(encoding="utf-8"))
