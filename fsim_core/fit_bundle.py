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


# D1 (devlib-r2-pha2): opt-in IBM-calibrated refit bundle (fitting.fit_phase0(..., lineshape="ibm"),
# written by scripts/fit_ibm_bundle.py).  Tracked file only; the Lorentzian bundle above is unchanged.
_IBM_PATH = _ROOT / "fsim_core" / "data" / "phase0_fit_params_ibm.json"


def load_fit_params(lineshape: str = "lorentzian") -> dict:
    """Full bundle JSON (keys include "params").  lineshape "lorentzian" (default, legacy, bit-identical)
    or "ibm" (opt-in IBM-calibrated refit, fsim_core/data/phase0_fit_params_ibm.json)."""
    if lineshape == "ibm":
        return json.loads(_IBM_PATH.read_text(encoding="utf-8"))
    if lineshape != "lorentzian":
        raise ValueError(f"lineshape must be 'lorentzian' or 'ibm', got {lineshape!r}")
    return json.loads(fit_params_path().read_text(encoding="utf-8"))
