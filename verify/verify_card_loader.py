"""fsim_core.card.load_card loader checks.

  1. every card in cards/*.yaml loads (24), including the two device-tier cards
     (chatzarakis2023, laferriere2023) whose datasets carry no per-dataset `source` and whose
     numbers live in a `device:` block (both raised KeyError('source') before the fix);
  2. for every card the pre-fix loader (reproduced verbatim below as LEGACY_LOAD) could read,
     the current loader returns an identical object: same meta, params (value, range, unit,
     tag, source), dataset rows/columns/tag/source, name and path;
  3. a dataset without `source` takes meta.source, a dataset with one keeps its own, and the
     device-tier cards get empty params and the datasets the YAML names.

Run: python verify/verify_card_loader.py   (exit code 0 iff all pass)
"""
import sys
import tempfile
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.card import Card, DataSet, Tag, _parse_param, load_card  # noqa: E402

KNOWN = {"chatzarakis", "chatzarakis2023", "edge-inp-gaasp-design", "edge-inp-gainp-design", "kitamura",
         "laferriere", "laferriere2023", "nitride-cavity-pulse-design", "nitride-cavity-set-design",
         "nitride-deshpande2014-comparison-design", "nitride-nanowire-horizontal-pulse-design",
         "nitride-nanowire-horizontal-set-design", "nitride-nanowire-vertical-pulse-design",
         "nitride-nanowire-vertical-set-design", "nitride-nonpolar-pulse-design", "nitride-nonpolar-set-design",
         "nitride-qw-fluctuation-pulse-design", "nitride-qw-fluctuation-set-design", "qcap-cavity",
         "qcap-piezo-variant", "qcap-staged", "reischle", "staged-device-design", "zhao"}
DEVICE_TIER = ("chatzarakis2023", "laferriere2023")
CHECKS = []


def check(name, ok, detail=""):
    CHECKS.append((name, bool(ok), detail))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""), flush=True)


def LEGACY_LOAD(path):
    """fsim_core/card.py load_card as it was before the device-tier fix (verbatim)."""
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


def main():
    # The 24 cards that exist for fsim_core.card. Later cards with their own loader (e.g.
    # fock-source-4K-design, loaded by resonant_source.load_fock_card) are out of scope here.
    paths = [p for p in sorted((ROOT / "cards").glob("*.yaml")) if p.stem in KNOWN]
    loaded, failed = {}, {}
    for p in paths:
        try:
            loaded[p.stem] = load_card(p)
        except Exception as exc:  # noqa: BLE001
            failed[p.stem] = f"{type(exc).__name__}: {exc}"
    check(f"all {len(paths)} cards in cards/ load", not failed and len(loaded) == len(paths) == 24,
          f"{len(loaded)} loaded; failed: {failed}")

    legacy_ok, same, diffs = [], [], []
    for p in paths:
        try:
            old = LEGACY_LOAD(p)
        except Exception:  # noqa: BLE001
            continue
        legacy_ok.append(p.stem)
        new = loaded.get(p.stem)
        if new is not None and old == new:  # dataclass equality: every field
            same.append(p.stem)
        else:
            diffs.append(p.stem)
    check("the cards the legacy loader could read are exactly the 22 non-device-tier cards",
          len(legacy_ok) == len(paths) - len(DEVICE_TIER) and not set(legacy_ok) & set(DEVICE_TIER),
          f"{len(legacy_ok)} legacy-loadable")
    check("every legacy-loadable card loads to an identical object (name, meta, params, datasets, path)",
          not diffs and len(same) == len(legacy_ok), f"{len(same)}/{len(legacy_ok)}; differing: {diffs}")

    for n in DEVICE_TIER:
        c = loaded.get(n)
        raw = yaml.safe_load((ROOT / "cards" / f"{n}.yaml").read_text(encoding="utf-8"))
        ok = (c is not None and c.params == {} and set(c.datasets) == set(raw["data"])
              and all(ds.source == str(raw["meta"]["source"]) for ds in c.datasets.values())
              and all(len(ds.rows) == len(raw["data"][k]["points"]) for k, ds in c.datasets.items()))
        check(f"{n}: loads with no params, its datasets, source = meta.source", ok)

    # a dataset with its own source keeps it; one without takes meta.source
    doc = {"meta": {"name": "t", "source": "META"},
           "data": {"a": {"tag": "V", "source": "OWN", "points": [{"x": 1.0}]},
                    "b": {"tag": "E", "points": [{"x": 2.0}]}}}
    with tempfile.TemporaryDirectory() as td:
        f = Path(td) / "t.yaml"
        f.write_text(yaml.safe_dump(doc), encoding="utf-8")
        c = load_card(f)
    check("per-dataset source is kept when present, falls back to meta.source when absent",
          c.datasets["a"].source == "OWN" and c.datasets["b"].source == "META" and c.params == {})

    n_pass = sum(ok for _, ok, _ in CHECKS)
    print(f"\n{n_pass}/{len(CHECKS)} card loader checks passed")
    return 0 if n_pass == len(CHECKS) else 1


if __name__ == "__main__":
    sys.exit(main())
