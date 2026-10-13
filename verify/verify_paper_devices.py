"""Paper-device gate (design 3.2): schema + class rules + reproduce + expectation ledger.

Run: python verify/verify_paper_devices.py [--only id] [--schema-only] [--self-test] [--jobs N]

For every card in cards/papers (files starting with "_" skipped) it validates the schema and class rules,
runs the composed engine at the paper's device and prints one line per ledger item. Exit 0 iff every
card passes schema + class rules AND every outcome equals the card's `expect` (known gaps are recorded
as `expect: fail|unavailable` with a reason; any flip in either direction breaks the run for review).
Only class M items count as tests; the honest count is "M tests (central within tolerance): x/y" -- a paper value
that is only inside the input envelope (outcome pass-envelope) is NOT a pass, and class R (paper theory /
simulation reference) is never a test. Weak passes (`weak: true` in the ledger) are counted separately.
"""
from __future__ import annotations

import argparse
import copy
import math
import sys
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core.compose import papers as PP                      # noqa: E402

CARDS = ROOT / "cards" / "papers"

# R2 audit 7: an M item that a NULL model (g2 -> 0, HOM overlap/visibility -> 1, efficiency/preparation -> ideal 1) would
# also pass cannot discriminate; it must carry `weak: true`. KNOWN_OPEN = (card id, item id) pairs that currently pass the
# null model without being flagged weak; the cards are flagged in the next wave, then delete the entries.
KNOWN_OPEN_NULL = {   # R2 audit 7; measured 2026-10-12: each passes the null model and is not flagged weak (tolerance wider than the effect)
    ("ding2016", "eta_prep_96"),
    ("ding2016", "hom_corr_2p1ns"),
    ("he2013", "hom_corr"),
    ("karli2022", "g2_corr"),
    ("muller2014", "g2_X_darkcorr"),
    ("muller2014", "g2_X_raw"),
    ("nawrath2022", "g2_analyt_76"),
    ("piccinini2025_tpe", "g2_TPE"),
    ("vajner2024", "g2_XX"),
    ("wang2019a", "hom_corr_pillar"),
}


def null_value(observable: str):
    """Value of the trivial model for an observable, or None when no null model is defined."""
    if observable.startswith("g2"):
        return 0.0
    if observable.startswith("hom_"):
        return 1.0
    if observable.startswith("eta_"):
        return 1.0
    return None


def _num(v) -> str:
    return "nan" if v is None or (isinstance(v, float) and math.isnan(v)) else f"{v:.4g}"


def item_line(pid: str, r) -> str:
    if r.kind == "scalar":
        paper = _num(r.paper) + (f"+-{_num(r.err)}" if r.err is not None else "")
        model = f"model {_num(r.model)} [{_num(r.lo)},{_num(r.hi)}]"
        what = f"paper {paper} {model}"
    elif r.kind == "series":
        what = f"series rms_rel {_num(r.model)}"
    else:
        what = f"prediction model {_num(r.model)}"
    tol = f" tol {r.tol_str}" if r.tol_str else ""
    why = f" -- {r.note}" if r.note and r.outcome in ("unavailable", "error") else ""
    weak = f" [WEAK: {r.weak_reason}]" if r.weak and r.outcome == "pass" else ""
    return (f"{pid} {r.id} {r.cls} {what}{tol} {r.outcome.upper()}{weak} "
            f"(expect {r.expect}) {'MATCH' if r.match else 'MISMATCH'}{why}")


def exit_code(reports, schema_ok: bool = True) -> int:
    """0 iff the schema/class gate passed and every outcome equals its expectation."""
    return 0 if schema_ok and all(rep.matched for rep in reports) else 1


def _run_card(path: str, envelope: bool = True):
    p = PP.load_paper(path)
    return PP.reproduce(p, envelope=envelope)


def run_cards(only=None, schema_only=False, jobs=1, directory=CARDS) -> int:
    files = sorted(f for f in Path(directory).glob("*.yaml") if not f.name.startswith("_")) \
        if Path(directory).is_dir() else []
    if only:
        files = [f for f in files if f.stem == only]
        if not files:
            print(f"no paper card with id {only!r}")
            return 1
    if not files:
        print("0 paper cards")
        print("0/0 paper-device checks matched ledger")
        return 0
    total = matched = 0
    schema_fail = 0
    loaded = []
    for f in files:                                           # schema + class rules first
        total += 1
        try:
            p = PP.load_paper(f)
        except PP.PaperSchemaError as e:
            schema_fail += 1
            print(f"{f.stem} SCHEMA/CLASS FAIL")
            for pr in e.problems:
                print(f"    {pr}")
            continue
        except Exception as e:                                # noqa: BLE001  (yaml errors etc.)
            schema_fail += 1
            print(f"{f.stem} LOAD FAIL: {type(e).__name__}: {e}")
            continue
        matched += 1
        loaded.append(p)
        print(f"{p.id} schema+class OK ({len(p.reproduce)} ledger items, "
              f"{len(p.inputs_from_paper)} paper-sourced inputs)")
    if schema_only:
        print(f"{matched}/{total} paper-device checks matched ledger (schema only)")
        return 0 if schema_fail == 0 else 1

    reports = []
    if jobs and jobs > 1 and len(loaded) > 1:
        from concurrent.futures import ProcessPoolExecutor
        with ProcessPoolExecutor(max_workers=jobs) as ex:
            reports = list(ex.map(_run_card, [p.path for p in loaded]))
    else:
        reports = [PP.reproduce(p) for p in loaded]
    m_n = m_pass = m_env = m_weak = 0
    r_tot = {o: 0 for o in PP.OUTCOMES}
    c_tot = {c: 0 for c in PP.CLASSES}
    for p, rep in zip(loaded, reports):
        for r in rep.items:
            total += 1
            matched += int(r.match)
            print(item_line(p.id, r))
        if not rep.valid:
            print(f"{p.id} note: central run invalid/unavailable: {'; '.join(rep.invalid_reasons)}")
        summ = rep.class_summary()
        print(f"{p.id} classes: " + ", ".join(
            f"{c}={s['n']} (pass {s['pass']}, envelope {s['pass-envelope']}, fail {s['fail']}, "
            f"unavail {s['unavailable']}, err {s['error']})"
            for c, s in summ.items() if s["n"]))
        m_n += summ["M"]["n"]
        m_pass += summ["M"]["pass"]
        m_env += summ["M"]["pass-envelope"]
        m_weak += sum(1 for r in rep.items if r.cls == "M" and r.outcome == "pass" and r.weak)
        for c in PP.CLASSES:
            c_tot[c] += summ[c]["n"]
        for o in PP.OUTCOMES:
            r_tot[o] += summ["R"][o]
    # ---- null-model discrimination (R2 audit 7)
    null_rows, null_bad, null_known = [], [], []
    for p, rep in zip(loaded, reports):
        for r in rep.items:
            if r.cls != "M" or r.kind != "scalar" or r.outcome != "pass" or r.expect != "pass":
                continue
            nv = null_value(r.observable)
            if nv is None or r.paper is None:
                continue
            if PP._scalar_outcome(r.paper, r.err, nv, 0.0, 0.0, r.tol) == "pass":
                null_rows.append((p.id, r.id, r.observable, r.paper, nv, r.weak))
                if not r.weak:
                    (null_known if (p.id, r.id) in KNOWN_OPEN_NULL else null_bad).append((p.id, r.id))
    print(f"null-model report: {len(null_rows)} M passes would also pass the null model "
          f"({sum(1 for x in null_rows if x[5])} flagged weak, {len(null_known)} KNOWN-OPEN, {len(null_bad)} unflagged)")
    for pid, iid, obs, pv, nv, weak in null_rows:
        tag = "weak" if weak else ("KNOWN-OPEN" if (pid, iid) in KNOWN_OPEN_NULL else "UNFLAGGED")
        print(f"  null passes {pid} {iid} ({obs}: paper {pv:.4g}, null {nv:g}) [{tag}]")
    stale = sorted(KNOWN_OPEN_NULL - {(a, b) for a, b, *_ in null_rows if True})
    for pid, iid in stale:
        print(f"  [FIXED?] KNOWN_OPEN_NULL entry ({pid}, {iid}) no longer passes the null model unflagged: delete it")
    total += 1
    matched += int(not null_bad)
    print(f"{'PASS' if not null_bad else 'FAIL'} every M pass that the null model also passes is flagged weak"
          + (f"  ({len(null_bad)} unflagged: {', '.join(f'{a}/{b}' for a, b in null_bad)})" if null_bad else ""))
    print(f"M discriminating passes: {m_pass - len(null_rows)}/{m_n} (M tests passing minus those the null model also passes; headline)")
    print(f"M tests (central within tolerance): {m_pass}/{m_n}")
    print(f"M pass-envelope (paper value inside the input envelope only, not a pass): {m_env}")
    print(f"M weak passes (flagged weak in the ledger, included in the {m_pass} above): {m_weak}")
    print(f"R paper theory/simulation references (not tests): {c_tot['R']} (pass {r_tot['pass']}, "
          f"envelope {r_tot['pass-envelope']}, fail {r_tot['fail']}, unavail {r_tot['unavailable']})")
    print("per-class item counts: " + ", ".join(f"{c}={c_tot[c]}" for c in PP.CLASSES))
    print(f"{matched}/{total} paper-device checks matched ledger")
    return exit_code(reports, schema_ok=schema_fail == 0 and not null_bad)


# ------------------------------------------------------------------ self-test
def self_test() -> int:
    res = []

    def check(name, ok, detail=""):
        res.append(bool(ok))
        print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail else ""))

    path = CARDS / "_fixture.yaml"
    raw = yaml.safe_load(path.read_text(encoding="utf-8"))
    try:
        paper = PP.load_paper(path)
        check("fixture loads and passes schema + class rules", True,
              f"{len(paper.reproduce)} items, inputs {sorted(paper.inputs_from_paper)}")
    except PP.PaperSchemaError as e:
        check("fixture loads and passes schema + class rules", False, "; ".join(e.problems))
        return 1
    check("fixture id is _fixture", paper.id == "_fixture")
    check("list_papers skips files starting with '_'",
          all(not p.path.split("/")[-1].startswith("_") for p in PP.list_papers(CARDS)))

    def rejected(name, mutate, expect_text):
        d = copy.deepcopy(raw)
        mutate(d)
        try:
            PP.paper_from_dict(d, stem="_fixture")
            check(name, False, "accepted a violating card")
        except PP.PaperSchemaError as e:
            hit = any(expect_text in pr for pr in e.problems)
            print(f"    rejected: {e.problems[0]}")
            check(name, hit, f"expected a problem containing {expect_text!r}")

    def lam_m(d):
        for e in d["reproduce"]:
            if e["item"] == "lambda":
                e["class"] = "M"
    rejected("(a) passthrough of a paper input labelled M is rejected", lam_m, "must be T")

    def fit_m(d):
        for e in d["reproduce"]:
            if e["item"] == "tau_vs_FP":
                e["class"] = "M"
    rejected("fitted item labelled M is rejected", fit_m, "must be C")

    def p_with_value(d):
        for e in d["reproduce"]:
            if e["item"] == "g2_prediction":
                e["item"] = "g2"
                e["class"] = "P"
        d["reproduce"] = [e for e in d["reproduce"] if not (e["item"] == "g2" and e["class"] == "M")]
    rejected("class P on an item that has a reported value is rejected", p_with_value, "no reported value")

    def no_where(d):
        del d["device"]["material"]["params"]["lambda_nm"]["where"]
    rejected("device leaf tagged V without where is rejected", no_where, "needs")

    def bad_conv(d):
        d["reported"][1]["convention"] = "corrected"
    rejected("convention invalid for the observable is rejected", bad_conv, "convention")

    def no_reason(d):
        for e in d["reproduce"]:
            if e["item"] == "tau_scalar":
                del e["reason"]
    rejected("expect != pass without a reason is rejected", no_reason, "reason required")

    def undigit(d):
        d["series"][0]["tag"] = "V"
    rejected("digitized series not tagged DR is rejected", undigit, "digitized")

    def no_tol_reason(d):
        del d["reproduce"][1]["tol_reason"]
    rejected("missing tol_reason is rejected", no_tol_reason, "tol_reason")

    rep = PP.reproduce(paper)
    by = {r.id: r for r in rep.items}
    check("fixture reproduces: every outcome equals its expectation", rep.matched,
          ", ".join(f"{r.id}={r.outcome}" for r in rep.items))
    check("T item model equals the transcribed input", by["lambda"].model == 900.0
          and by["lambda"].inputs == ["material.lambda_nm"])
    check("ranged input: paper brightness lies inside the model [min,max]",
          by["brightness"].lo < by["brightness"].paper < by["brightness"].hi
          and by["brightness"].lo < by["brightness"].hi)
    be = by["brightness_env"]
    check("FX4: central misses the tolerance but the paper value is inside [min,max] -> pass-envelope (not pass)",
          be.outcome == "pass-envelope" and be.lo < be.paper < be.hi and abs(be.model - be.paper) > 2 * be.err
          and be.match)
    summ = rep.class_summary()
    check("FX4: pass-envelope is not counted in the M pass count", summ["M"]["pass"] == 2 and summ["M"]["pass-envelope"] == 1
          and summ["M"]["n"] == 4, str(summ["M"]))
    check("FX4: class R item is reported, counted as R and not as an M test",
          by["g2_theory"].cls == "R" and by["g2_theory"].outcome == "pass" and summ["R"]["n"] == 1)

    def weak_no_reason(d):
        for e in d["reproduce"]:
            if e["item"] == "g2":
                e["weak"] = True
    rejected("FX4: weak: true without weak_reason is rejected", weak_no_reason, "weak_reason")

    def env_expect_pass(d):
        for e in d["reproduce"]:
            if e["item"] == "brightness_env":
                e["expect"] = "pass"
                e.pop("reason", None)
    d_e = copy.deepcopy(raw)
    env_expect_pass(d_e)
    rep_e = PP.reproduce(PP.paper_from_dict(d_e, stem="_fixture"))
    check("FX4: ledger expecting pass on an envelope-only item mismatches (an envelope is not a pass)",
          not rep_e.matched and exit_code([rep_e]) == 1)

    wk = copy.deepcopy(raw)
    for e in wk["reproduce"]:
        if e["item"] == "g2":
            e["weak"], e["weak_reason"] = True, "synthetic weak pass"
    rep_w = PP.reproduce(PP.paper_from_dict(wk, stem="_fixture"))
    gw = {r.id: r for r in rep_w.items}["g2"]
    check("FX4: weak flag and reason carried on the result", gw.weak and gw.weak_reason == "synthetic weak pass"
          and "[WEAK: synthetic weak pass]" in item_line("_fixture", gw))
    check("known gap is recorded as expect fail and is a fail", by["tau_scalar"].outcome == "fail"
          and by["tau_scalar"].match)
    check("P item reports a model number and no comparison",
          by["g2_prediction"].paper is None and math.isfinite(by["g2_prediction"].model))
    check("exit-code logic: matched ledger -> 0", exit_code([rep]) == 0)

    flipped = copy.deepcopy(raw)
    for e in flipped["reproduce"]:
        if e["item"] == "tau_scalar":
            e["expect"] = "pass"
            e.pop("reason", None)
    rep_f = PP.reproduce(PP.paper_from_dict(flipped, stem="_fixture"))
    check("(b) flipped expect breaks the run: exit-code logic returns 1",
          exit_code([rep_f]) == 1 and not rep_f.matched,
          "; ".join(f"{r.id} {r.outcome} vs expect {r.expect}" for r in rep_f.items if not r.match))
    check("schema failure forces exit 1", exit_code([rep], schema_ok=False) == 1)

    unavail = copy.deepcopy(raw)
    # FX2c: arp is ready now; register a throwaway needs_physics drive (restored afterwards, never shipped)
    from fsim_core.compose import drives as _D
    from fsim_core.compose import registry as _R
    thr = _D._placeholder("test_needs_physics", "throwaway", "self-test only", "self-test placeholder",
                          (_D._GAUSS,), ("none",), _D.OBS_COHERENT_PULSED)
    _R.register(thr)
    unavail["device"]["drive"] = {"kind": "test_needs_physics", "params": {}}
    for e in unavail["reproduce"]:
        e["expect"], e["reason"] = "unavailable", "drive test_needs_physics is awaiting physics"
    try:
        rep_u = PP.reproduce(PP.paper_from_dict(unavail, stem="_fixture"))
    finally:
        _R.unregister("drive", "test_needs_physics")
    check("needs_physics drive -> every item unavailable, ledger matched",
          all(r.outcome == "unavailable" for r in rep_u.items) and rep_u.matched,
          rep_u.items[0].note)

    n_ok = sum(res)
    print(f"{n_ok}/{len(res)} paper-device self-test checks passed")
    return 0 if n_ok == len(res) else 1


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--only", help="check one card id")
    ap.add_argument("--schema-only", action="store_true", help="schema and class rules only")
    ap.add_argument("--self-test", action="store_true", help="run on cards/papers/_fixture.yaml")
    ap.add_argument("--jobs", type=int, default=1, help="process-parallel cards")
    a = ap.parse_args(argv)
    if a.self_test:
        return self_test()
    return run_cards(a.only, a.schema_only, a.jobs)


if __name__ == "__main__":
    sys.exit(main())
