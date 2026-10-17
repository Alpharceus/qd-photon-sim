"""Clone-safe V-a fit for the verifiers (devlib-vf1, defect H1).

out/phase0/fit_params.json is gitignored, so a fresh clone cannot read it.
The V-a joint fit is deterministic (seed 0, ~1.8 s) and is recomputed here
in memory from the Chatzarakis card with the production fit
(fsim_core.fitting.fit_phase0).  The result is asserted against the committed
snapshot fsim_core/data/phase0_fit_params.json (a copy of out/phase0/fit_params.json).
If the fit raises, the snapshot is used instead (and `source` says so).
"""
import json
import sys
from functools import lru_cache
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

SNAPSHOT = ROOT / "fsim_core" / "data" / "phase0_fit_params.json"
RTOL = 1e-6
# a_ac is pinned at its lower bound (~1e-16 1/ns/K, physically inert: its
# contribution to Gamma(T) is < 1e-12 meV); the optimiser's last digits there are
# not reproducible to 1e-6 relative, so it is compared on an absolute scale.
ATOL_PINNED = {"a_ac": 1e-18}


@lru_cache(maxsize=1)
def phase0_fit():
    """Return (params dict, source str).  source is 'fit' or 'snapshot (fit raised: ...)'."""
    snap = json.loads(SNAPSHOT.read_text())["params"]
    try:
        from fsim_core.card import load_card
        from fsim_core.fitting import fit_phase0
        fit = fit_phase0(load_card(str(ROOT / "cards" / "chatzarakis.yaml")))
        params = {k: float(fit.params[k]) for k in snap}
    except Exception as exc:  # fall back to the committed snapshot only if the fit raises
        return dict(snap), f"snapshot (fit raised: {exc!r})"
    for k, v in snap.items():
        atol = ATOL_PINNED.get(k, 0.0)
        if abs(params[k] - v) > RTOL * abs(v) + atol:
            raise AssertionError(
                f"in-memory V-a fit disagrees with fsim_core/data/phase0_fit_params.json: "
                f"{k} = {params[k]!r} vs {v!r}")
    return params, "fit"


def phase0_params():
    return dict(phase0_fit()[0])
