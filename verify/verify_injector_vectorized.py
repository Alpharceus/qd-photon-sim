"""Equality checks for the studio-p2c CPU speed-up (no physics change).

1. The vectorized nanowire tunnel-injector engine
   (fsim_core.nitride_nanowire_injector._transmission_reflection_vec) against
   the retained scalar reference (_transmission_reflection_scalar) on dense
   energy grids covering the full domain the cards use (both carriers, every
   nanowire card, the card fields plus zero field and a bias sweep, dense
   zooms on every resonance, the evanescent/propagating transitions at every
   tilted band edge, and the window edges): <= 1e-12 relative on T and R.
2. Whole-card regression: evaluate() of every nanowire and planar nitride card
   (and the two edge cards, which exercise dot_levels) at 3 temperatures,
   compared leaf-by-leaf with verify/data/perf_baseline.json, which was
   captured from the pre-speed-up code by `--capture`: <= 1e-10 relative on
   every scalar and curve, exact equality on every non-numeric leaf.

Exits 0 iff all checks pass and prints "N/N ... passed".
"""
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

BASELINE = ROOT / "verify" / "data" / "perf_baseline.json"
CARDS = [
    "nitride-nanowire-vertical-pulse-design",
    "nitride-nanowire-vertical-set-design",
    "nitride-nanowire-horizontal-pulse-design",
    "nitride-nanowire-horizontal-set-design",
    "nitride-cavity-pulse-design",
    "nitride-cavity-set-design",
    "nitride-deshpande2014-comparison-design",
    "nitride-nonpolar-pulse-design",
    "nitride-nonpolar-set-design",
    "nitride-qw-fluctuation-pulse-design",
    "nitride-qw-fluctuation-set-design",
    "edge-inp-gainp-design",
    "edge-inp-gaasp-design",
]
TEMPS = (230.0, 265.0, 300.0)
CARD_RTOL = 1e-10
ENGINE_RTOL = 1e-12

c = []


def check(label, value):
    c.append(bool(value))
    if not value:
        print("FAIL", label)


def _to_json(x):
    """Recursive, lossless-enough JSON form of an evaluate() result: floats
    stay floats (json writes repr, which round-trips exactly), arrays become
    lists, NaN/inf are written as Python's json extensions."""
    if isinstance(x, dict):
        return {str(k): _to_json(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [_to_json(v) for v in x]
    if isinstance(x, np.ndarray):
        return [_to_json(v) for v in x.tolist()]
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (int, np.integer)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        return float(x)
    if x is None or isinstance(x, str):
        return x
    return repr(x)


def _evaluate(card, T):
    from fsim_core.device import DeviceDesign, evaluate
    d = DeviceDesign.load(ROOT / "cards" / f"{card}.yaml")
    return _to_json(evaluate(d, T_grid=[T]))


def capture():
    out = {"temps": list(TEMPS), "cards": {}}
    for card in CARDS:
        out["cards"][card] = {}
        for T in TEMPS:
            t0 = time.perf_counter()
            out["cards"][card][repr(T)] = _evaluate(card, T)
            print(f"captured {card} T={T}: {time.perf_counter() - t0:.2f} s", flush=True)
    BASELINE.parent.mkdir(parents=True, exist_ok=True)
    BASELINE.write_text(json.dumps(out, indent=1, sort_keys=True))
    print("wrote", BASELINE)


def _num_close(a, b, rtol):
    if isinstance(a, bool) or isinstance(b, bool):
        return a is b
    a, b = float(a), float(b)
    if math.isnan(a) or math.isnan(b):
        return math.isnan(a) and math.isnan(b)
    if math.isinf(a) or math.isinf(b):
        return a == b
    return abs(a - b) <= rtol * max(abs(a), abs(b))


def _compare(path, a, b, rtol, worst, mism):
    """Leaf-by-leaf comparison; records the worst relative deviation."""
    if isinstance(a, dict) and isinstance(b, dict):
        if set(a) != set(b):
            mism.append(f"{path}: key sets differ {sorted(set(a) ^ set(b))}")
        for k in sorted(set(a) & set(b)):
            _compare(f"{path}.{k}", a[k], b[k], rtol, worst, mism)
        return
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            mism.append(f"{path}: length {len(a)} != {len(b)}")
            return
        for i, (x, y) in enumerate(zip(a, b)):
            _compare(f"{path}[{i}]", x, y, rtol, worst, mism)
        return
    numeric = (int, float)
    if (isinstance(a, numeric) and not isinstance(a, bool)
            and isinstance(b, numeric) and not isinstance(b, bool)):
        if not _num_close(a, b, rtol):
            mism.append(f"{path}: {a!r} != {b!r}")
        fa, fb = float(a), float(b)
        if math.isfinite(fa) and math.isfinite(fb) and max(abs(fa), abs(fb)) > 0:
            rel = abs(fa - fb) / max(abs(fa), abs(fb))
            if rel > worst[0]:
                worst[0], worst[1] = rel, path
        return
    if a != b:
        mism.append(f"{path}: {a!r} != {b!r}")


class _Recorder:
    """Wraps the injector's per-energy and grid entry points during the card
    runs to collect every (stack, energy) the cards actually evaluate."""

    def __init__(self, inj):
        self.inj = inj
        self.calls = {}   # (params, bias, field, carrier) -> list of energies
        self._orig_s = inj._transmission_scalar
        self._orig_v = inj._transmission_vec

    def __enter__(self):
        inj, calls = self.inj, self.calls

        def s(params, E, bias_V, field_kVcm, carrier):
            calls.setdefault((params, float(bias_V), float(field_kVcm), carrier), []).append(float(E))
            return self._orig_s(params, E, bias_V, field_kVcm, carrier)

        def v(params, Es, bias_V, field_kVcm, carrier):
            calls.setdefault((params, float(bias_V), float(field_kVcm), carrier), []).extend(
                np.atleast_1d(np.asarray(Es, dtype=float)).tolist())
            return self._orig_v(params, Es, bias_V, field_kVcm, carrier)

        inj._transmission_scalar, inj._transmission_vec = s, v
        return self

    def __exit__(self, *exc):
        self.inj._transmission_scalar = self._orig_s
        self.inj._transmission_vec = self._orig_v


def _test_energies(energies):
    """Base test grid for one stack: a 4001-point uniform grid over the
    full energy domain the cards used (plus 5 % headroom), the domain edges,
    and an evenly strided subsample (about 3000) of the energies the cards
    actually evaluated.  engine_checks adds the resonance zooms and the
    band-edge (evanescent/propagating) transitions."""
    E = np.unique(np.asarray(energies, dtype=float))
    lo, hi = max(float(E.min()), 1e-9), float(E.max()) * 1.05
    dense = np.linspace(lo, hi, 4001)
    out = [dense, E[:: max(1, E.size // 3000)], np.array([lo, hi, 1e-9])]
    return dense, out


def engine_checks(recorded):
    from fsim_core import nitride_nanowire_injector as inj
    stacks = {}
    for (params, bias_V, field_kVcm, carrier), Es in recorded.items():
        key = inj._engine_inputs(params, bias_V, field_kVcm, carrier)
        stacks.setdefault(key, (params, bias_V, field_kVcm, carrier, []))[4].extend(Es)
    worst_T = worst_R = 0.0
    n_pts = n_exact = 0
    cached_exact = True
    for key, (params, bias_V, field_kVcm, carrier, Es) in stacks.items():
        segs, m_well, tilt, v_right_eV, pol, sl, ms, slices = key
        dense, parts = _test_energies(Es)
        T_d = inj._transmission_vec(params, dense, bias_V, field_kVcm, carrier)[0]
        # resonance zooms: every interior local maximum of the dense grid
        step = dense[1] - dense[0]
        for j in range(1, dense.size - 1):
            if T_d[j] >= T_d[j - 1] and T_d[j] >= T_d[j + 1] and T_d[j] > 0.0:
                parts.append(np.linspace(dense[j] - 2 * step, dense[j] + 2 * step, 401))
        # evanescent/propagating transitions at every distinct band edge
        edges = sorted({float(V) for (_dl, V, _m) in slices} | {0.0, float(v_right_eV)})
        trans = []
        for V in edges:
            for d in (-1e-6, -1e-9, -1e-12, 0.0, 1e-12, 1e-9, 1e-6):
                x = V + d * max(abs(V), 1e-3)
                if x > 0.0:
                    trans.append(x)
            if V > 0.0:
                trans.extend([float(np.nextafter(V, np.inf)), float(np.nextafter(V, -np.inf))])
        parts.append(np.asarray(trans, dtype=float))
        grid = np.unique(np.concatenate(parts))
        grid = grid[grid > 0.0]
        T_v, R_v = inj._transmission_vec(params, grid, bias_V, field_kVcm, carrier)
        for j, E in enumerate(grid.tolist()):
            T_s, R_s = inj._transmission_scalar_reference(params, E, bias_V, field_kVcm, carrier)
            n_pts += 1
            if T_s == T_v[j] and R_s == R_v[j]:
                n_exact += 1
            for ref, got, kind in ((T_s, T_v[j], "T"), (R_s, R_v[j], "R")):
                if math.isnan(ref) or math.isnan(got):
                    rel = 0.0 if (math.isnan(ref) and math.isnan(got)) else math.inf
                else:
                    scale = max(abs(ref), abs(got))
                    rel = abs(ref - got) / scale if scale > 0 else 0.0
                if kind == "T":
                    worst_T = max(worst_T, rel)
                else:
                    worst_R = max(worst_R, rel)
            if j % 37 == 0:   # cached scalar path vs reference, bitwise
                if inj._transmission_scalar(params, E, bias_V, field_kVcm, carrier) != (T_s, R_s):
                    cached_exact = False
    print(f"engine: {len(stacks)} distinct stacks, {n_pts} energies, "
          f"{n_exact} bit-identical; max rel dev T {worst_T:.3e}, R {worst_R:.3e}")
    check("engine stacks recorded from the nanowire cards", len(stacks) > 0)
    check(f"vectorized T within {ENGINE_RTOL} relative of scalar", worst_T <= ENGINE_RTOL)
    check(f"vectorized R within {ENGINE_RTOL} relative of scalar", worst_R <= ENGINE_RTOL)
    check("cached scalar path bit-identical to the reference", cached_exact)
    return worst_T, worst_R


def card_checks(baseline):
    from fsim_core import nitride_nanowire_injector as inj
    per_card = {}
    with _Recorder(inj) as rec:
        for card in CARDS:
            worst = [0.0, ""]
            mism = []
            for T in TEMPS:
                got = _evaluate(card, T)
                _compare(f"{card}@{T}", baseline["cards"][card][repr(T)], got, CARD_RTOL, worst, mism)
            for m in mism[:5]:
                print("  mismatch", m)
            per_card[card] = worst
            print(f"card {card}: max rel dev {worst[0]:.3e}" + (f" at {worst[1]}" if worst[0] else ""))
            check(f"{card} matches baseline to {CARD_RTOL}", not mism)
    return rec.calls


if __name__ == "__main__":
    if "--capture" in sys.argv:
        capture()
        sys.exit(0)
    import warnings
    warnings.simplefilter("ignore")
    baseline = json.loads(BASELINE.read_text())
    check("baseline temperatures unchanged", tuple(baseline["temps"]) == TEMPS)
    recorded = card_checks(baseline)
    engine_checks(recorded)
    n_ok = sum(c)
    print(f"{n_ok}/{len(c)} injector-vectorization checks passed")
    sys.exit(0 if n_ok == len(c) else 1)
