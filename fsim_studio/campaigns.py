"""out/ scanner: campaign manifests, VERDICT / BEST line parsing, typed
sweep.csv tables and filtered queries, figure listings.

Verdicts are only ever parsed from VERDICT/BEST lines in the committed
results.md / verdict.md (never recomputed). Values are kept as the exact text
printed in the file (no rounding); the UI rounds. pandas is not installed:
csv module + plain lists.
"""
from __future__ import annotations

import csv
import json
import math
import re
import subprocess
import threading
from functools import lru_cache
from pathlib import Path

from . import OUT, ROOT

EXCLUDED_DIRS = {"c6_scratch_lw_quick", "quick", "mplconfig", "presentation", "presentation2",
                 "presentation2_selftest", "_scratch"}
LEGACY_DIRS = ("phase0", "phase1", "phase2", "phase3", "spec", "validation", "zhao")
COLLECTION_DIRS = ("rt_campaign", "tier_geometry", "tier_device")
PRE_AUDIT = {
    "tier_geometry": ["corner_frontier.csv"],
    "tier_device": None,  # every file
}
PRE_AUDIT_NOTE = "pre-audit (2026-08-12), may carry H5 1-D Purcell"
NULL_TOKENS = {"", "nan", "NaN", "None", "not_applicable"}
FIGURE_EXT = {".png", ".svg", ".pdf"}
HEADLINE_FILTER = {"rt_edge": {"model_finite_pulse": True, "model_tau_cap_density": False}}

KEY_RE = re.compile(r"^[A-Za-z_][\w.]*=")


class CampaignError(Exception):
    def __init__(self, status: int, message: str):
        super().__init__(message)
        self.status = status
        self.message = message


# ----------------------------------------------------------- line parsing

def parse_kv_line(line: str) -> dict:
    """'VERDICT: FAIL k=v k=v ...' or 'BEST_PASSING_FLUX k=v ...' -> {raw, word,
    fields, note}. Tolerates prose: bare tokens after a key=value continue
    that value (e.g. 'model=finite electrical pulse'); a parenthesised run
    '( ... )' goes to note; bare tokens before any key form the word."""
    raw = line.strip()
    body = raw
    head = None
    if body.startswith("VERDICT:"):
        head = "VERDICT"
        body = body[len("VERDICT:"):].strip()
    else:
        m = re.match(r"^(BEST[\w]*)\s*:?\s*(.*)$", body)
        if m:
            head = m.group(1)
            body = m.group(2)
    fields, word_parts, notes = {}, [], []
    last = None
    in_paren = False
    for tok in body.split():
        if in_paren or tok.startswith("("):
            notes.append(tok)
            in_paren = not tok.endswith(")")
            continue
        if KEY_RE.match(tok):
            k, _, v = tok.partition("=")
            fields[k] = v
            last = k
        elif last is None:
            word_parts.append(tok)
        else:
            fields[last] = f"{fields[last]} {tok}"
    word = " ".join(word_parts) if word_parts else fields.get("idealized_status")
    out = {"raw": raw, "kind": head, "word": word, "fields": fields}
    if notes:
        out["note"] = " ".join(notes).strip("()")
    return out


def _md_file(cdir: Path) -> Path | None:
    for name in ("verdict.md", "results.md"):
        p = cdir / name
        if p.is_file():
            return p
    return None


def parse_md(path: Path) -> tuple[list, list]:
    verdicts, best = [], []
    if path is None or not path.is_file():
        return verdicts, best
    for i, line in enumerate(path.read_text(encoding="utf-8").splitlines(), start=1):
        s = line.strip()
        if s.startswith("VERDICT:"):
            v = parse_kv_line(s)
            v["line"] = i
            verdicts.append(v)
        elif re.match(r"^BEST[A-Z_0-9]*\s", s):
            b = parse_kv_line(s)
            b["line"] = i
            best.append(b)
    return verdicts, best


# ----------------------------------------------------------- discovery

def _rel(p: Path) -> str:
    return p.relative_to(OUT).as_posix()


@lru_cache(maxsize=64)
def _git_last(rel_dir: str) -> tuple:
    try:
        out = subprocess.run(["git", "log", "-1", "--format=%h|%cI", "--", f"out/{rel_dir}"],
                             cwd=ROOT, capture_output=True, text=True, timeout=10)
        if out.returncode == 0 and "|" in out.stdout:
            h, d = out.stdout.strip().split("|", 1)
            return h, d
    except (OSError, subprocess.SubprocessError):
        pass
    return None, None


def _manifest_dirs() -> list:
    dirs = []
    for p in sorted(OUT.glob("*/manifest.json")) + sorted(OUT.glob("*/*/manifest.json")):
        cdir = p.parent
        if any(part in EXCLUDED_DIRS for part in cdir.relative_to(OUT).parts):
            continue
        dirs.append(cdir)
    return dirs


def _titleize(cid: str) -> str:
    return cid.replace("/", " / ").replace("_", " ")


def _status(verdicts: list) -> tuple:
    words = {}
    for v in verdicts:
        if v["word"]:
            words[v["word"]] = words.get(v["word"], 0) + 1
    if not words:
        return None, words
    return (next(iter(words)) if len(words) == 1 else "mixed"), words


def _load_manifest(cdir: Path) -> dict | None:
    p = cdir / "manifest.json"
    if not p.is_file():
        return None
    try:
        return json.loads(p.read_text(encoding="utf-8"))
    except (OSError, ValueError):
        return None


def _commit_and_date(cid: str, manifest: dict | None) -> tuple:
    commit = None
    generated = None
    if manifest:
        commit = manifest.get("report_only_generation_commit") or manifest.get("generation_commit")
        generated = manifest.get("generated_utc")
    gh, gd = _git_last(cid)
    if commit:
        return commit[:7], (generated or gd), "manifest"
    return gh, (generated or gd), ("git log (last commit touching the folder)" if gh else None)


def _stale(cid: str) -> tuple:
    if cid in PRE_AUDIT:
        files = PRE_AUDIT[cid]
        return True, files
    return False, None


def list_campaigns() -> list:
    out = []
    for cdir in _manifest_dirs():
        cid = _rel(cdir)
        manifest = _load_manifest(cdir)
        verdicts, best = parse_md(_md_file(cdir))
        word, words = _status(verdicts)
        commit, generated, csrc = _commit_and_date(cid, manifest)
        out.append({"id": cid, "title": _titleize(cid), "kind": "campaign",
                    "path": f"out/{cid}", "has_manifest": True,
                    "verdicts": verdicts, "best": best, "status_word": word,
                    "status_words": words, "generated": generated, "commit": commit,
                    "commit_source": csrc,
                    "stale": False, "stale_files": None})
    for name in COLLECTION_DIRS + LEGACY_DIRS:
        cdir = OUT / name
        if not cdir.is_dir():
            continue
        stale, stale_files = _stale(name)
        commit, generated, csrc = _commit_and_date(name, None)
        legacy = name in LEGACY_DIRS
        out.append({"id": name, "title": _titleize(name),
                    "kind": "legacy" if legacy else "collection",
                    "collection_label": "legacy F-series" if legacy else "figures + CSV",
                    "path": f"out/{name}", "has_manifest": False, "verdicts": [], "best": [],
                    "status_word": None, "status_words": {}, "generated": generated,
                    "commit": commit, "commit_source": csrc, "stale": stale, "stale_files": stale_files,
                    "stale_note": PRE_AUDIT_NOTE if stale else None})
    return out


def campaign_dir(cid: str) -> Path:
    if not cid or ".." in cid.split("/") or cid.startswith("/"):
        raise CampaignError(400, f"invalid campaign id {cid!r}")
    known = {c["id"] for c in list_campaigns()}
    if cid not in known:
        raise CampaignError(404, f"no campaign {cid!r}")
    return OUT / cid


def _files(cdir: Path) -> dict:
    csvs, figs, mds, other = [], [], [], []
    for p in sorted(cdir.iterdir()):
        if not p.is_file():
            continue
        rel = _rel(p)
        ext = p.suffix.lower()
        if ext == ".csv":
            csvs.append(rel)
        elif ext in FIGURE_EXT:
            figs.append(rel)
        elif ext == ".md":
            mds.append(rel)
        elif ext in (".json", ".html", ".yaml"):
            other.append(rel)
    return {"csv": csvs, "figures": figs, "md": mds, "other": other}


# ----------------------------------------------------------- typed tables

def _typed(text: str):
    if text in NULL_TOKENS:
        return None
    if text == "True":
        return True
    if text == "False":
        return False
    try:
        f = float(text)
    except ValueError:
        return text
    return f if math.isfinite(f) else None


class Table:
    def __init__(self, path: Path):
        self.path = path
        self.mtime = path.stat().st_mtime
        with path.open(newline="", encoding="utf-8") as fh:
            reader = csv.reader(fh)
            self.columns = next(reader, [])
            self.rows = [[_typed(c) for c in row] for row in reader]
        self.index = {c: i for i, c in enumerate(self.columns)}
        self._kinds = None

    def kinds(self) -> dict:
        if self._kinds is None:
            kinds = {}
            for i, c in enumerate(self.columns):
                seen = set()
                for r in self.rows:
                    v = r[i] if i < len(r) else None
                    if v is None:
                        continue
                    seen.add("bool" if isinstance(v, bool) else
                             "number" if isinstance(v, float) else "string")
                    if len(seen) > 1:
                        break
                kinds[c] = (seen.pop() if len(seen) == 1 else
                            "mixed" if seen else "empty")
            self._kinds = kinds
        return self._kinds

    def cardinality(self, col: str, cap: int = 9) -> int:
        i = self.index[col]
        vals = set()
        for r in self.rows:
            vals.add(r[i] if i < len(r) else None)
            if len(vals) >= cap:
                break
        return len(vals)


_TABLES: dict = {}
_TABLES_LOCK = threading.Lock()


def get_table(path: Path) -> Table:
    key = str(path)
    with _TABLES_LOCK:
        t = _TABLES.get(key)
        if t is not None and t.mtime == path.stat().st_mtime:
            return t
    t = Table(path)
    with _TABLES_LOCK:
        _TABLES[key] = t
    return t


_VERDICT_HINTS = ("eligible", "pass", "valid", "qualified", "feasible", "invalid_reason",
                  "runaway", "converged", "measurable", "bound_reversal", "compatible")
_INPUT_HINTS = {"card_id", "card_class", "config_id", "card_file", "row_id", "row_kind",
                "family", "regime", "orientation", "strain_bound", "bound_role", "T_hs",
                "T_hs_K", "T_hs_requested", "cavity_tracking", "screening_fraction",
                "height_nm", "radius_nm", "core_radius_nm", "x_in", "Q", "current_uA",
                "rep_rate_hz", "occupied_dot_access", "model_finite_pulse",
                "model_tau_cap_density", "emission_NA", "emission_R_back", "emission_L_um",
                "emission_alpha_cm", "irf_ps", "delta_xx_meV", "gamma300_meV", "b_res",
                "bias_mode", "V_j", "group", "screening", "voltage_window_V",
                "sensitivity_axis", "sensitivity_value"}


def column_roles(table: Table, manifest: dict | None) -> list:
    axes = set()
    if manifest:
        for k in ("grid", "axes"):
            v = manifest.get(k)
            if isinstance(v, dict):
                for name in v:
                    axes.add(name)
                    axes.add(name.split(".")[-1])
    kinds = table.kinds()
    out = []
    for c in table.columns:
        kind = kinds[c]
        low = c.lower()
        if c in _INPUT_HINTS or c in axes or c.endswith("_tag"):
            role = "input"
        elif kind == "bool" or any(h in low for h in _VERDICT_HINTS):
            role = "verdict"
        else:
            role = "output"
        out.append({"name": c, "role": role, "kind": kind})
    return out


def _primary_csv(cdir: Path) -> Path | None:
    p = cdir / "sweep.csv"
    if p.is_file():
        return p
    csvs = sorted(cdir.glob("*.csv"))
    return csvs[0] if csvs else None


def get_campaign(cid: str) -> dict:
    cdir = campaign_dir(cid)
    manifest = _load_manifest(cdir)
    md = _md_file(cdir)
    verdicts, best = parse_md(md)
    word, words = _status(verdicts)
    commit, generated, csrc = _commit_and_date(cid, manifest)
    primary = _primary_csv(cdir)
    columns = column_roles(get_table(primary), manifest) if primary else []
    stale, stale_files = _stale(cid.split("/")[0])
    return {"id": cid, "manifest": manifest, "verdicts": verdicts, "best": best,
            "status_word": word, "status_words": words, "commit": commit,
            "commit_source": csrc, "generated": generated, "md": _rel(md) if md else None,
            "files": _files(cdir), "primary_csv": _rel(primary) if primary else None,
            "columns": columns, "stale": stale, "stale_files": stale_files,
            "stale_note": PRE_AUDIT_NOTE if stale else None,
            "headline_filter": HEADLINE_FILTER.get(cid)}


# ----------------------------------------------------------- queries

def _norm_filter_value(v):
    if isinstance(v, str):
        if v == "True":
            return True
        if v == "False":
            return False
        if v in NULL_TOKENS:
            return None
        try:
            return float(v)
        except ValueError:
            return v
    if isinstance(v, bool) or v is None:
        return v
    if isinstance(v, (int, float)):
        return float(v)
    return v


def _match_fn(spec):
    if isinstance(spec, dict) and ("min" in spec or "max" in spec):
        lo = spec.get("min")
        hi = spec.get("max")

        def f(x):
            if x is None or isinstance(x, (bool, str)):
                return False
            return (lo is None or x >= lo) and (hi is None or x <= hi)
        return f
    if isinstance(spec, list):
        allowed = [_norm_filter_value(v) for v in spec]

        def f(x):
            return any(x == a and type(x) is type(a) for a in allowed)
        return f
    target = _norm_filter_value(spec)
    return lambda x: x == target and type(x) is type(target)


def query_sweep(cid: str, file: str | None = None, filt: dict | None = None,
                cols: list | None = None, limit: int = 5000, headline: bool = False,
                offset: int = 0) -> dict:
    cdir = campaign_dir(cid)
    if file:
        if "/" in file or "\\" in file or ".." in file or not file.endswith(".csv"):
            raise CampaignError(400, f"invalid file {file!r}")
        path = cdir / file
    else:
        path = _primary_csv(cdir)
    if path is None or not path.is_file():
        raise CampaignError(404, f"no csv {file!r} in campaign {cid!r}")
    t = get_table(path)
    conds = {}
    if headline:
        hf = HEADLINE_FILTER.get(cid)
        if hf is None:
            raise CampaignError(400, f"campaign {cid!r} has no headline-model filter")
        conds.update(hf)
    if filt:
        if not isinstance(filt, dict):
            raise CampaignError(400, "filter must be a JSON object")
        if headline:
            # H12: a user filter may not override the headline-model pin while
            # the response still says headline:true.
            clash = sorted(k for k, v in filt.items()
                           if k in conds and _norm_filter_value(v) != _norm_filter_value(conds[k]))
            if clash:
                raise CampaignError(400, f"filter keys {clash} conflict with the headline-model pin "
                                         f"{HEADLINE_FILTER.get(cid)}; drop headline=1 to query other models")
        conds.update(filt)
    checks = []
    for col, spec in conds.items():
        if col not in t.index:
            raise CampaignError(400, f"unknown column {col!r}")
        checks.append((t.index[col], _match_fn(spec)))
    if cols:
        missing = [c for c in cols if c not in t.index]
        if missing:
            raise CampaignError(400, f"unknown columns {missing}")
        out_cols = list(cols)
    else:
        out_cols = list(t.columns)
    idx = [t.index[c] for c in out_cols]
    ncols = len(t.columns)
    matched = []
    for r in t.rows:
        if len(r) < ncols:
            r = r + [None] * (ncols - len(r))
        if all(fn(r[i]) for i, fn in checks):
            matched.append(r)
    page = matched[offset:offset + max(0, int(limit))] if limit is not None else matched[offset:]
    rows = [[r[i] for i in idx] for r in page]
    return {"file": _rel(path), "columns": out_cols, "rows": rows, "total": len(t.rows),
            "filtered": len(matched), "returned": len(rows), "offset": offset,
            "headline": bool(headline), "filter": conds}


# ----------------------------------------------------------- gates (I3)

GATE_FILES = {"verdict": "rt_edge/verdict.md", "contract": ROOT / "docs" / "rt_edge_contract.md"}
_G2_GATE_RE = re.compile(r"pulsed intrinsic g2\(0\) < ([0-9.]+)")
_G2_CONTRACT_RE = re.compile(r"pulsed intrinsic `g2\(0\)` is below ([0-9.]+)")
_FLUX_FLOOR_RE = re.compile(r"eligibility floor \[([A-Z/]+)\]: ([0-9.eE+]+) (photons/s)")


def _first_match(path: Path, rx):
    if not path.is_file():
        return None, None
    text = path.read_text(encoding="utf-8")
    m = rx.search(text)
    if not m:
        return None, None
    return m, text.count("\n", 0, m.start()) + 1


def gates() -> dict:
    """The two acceptance gates, parsed (never typed in): the g2(0) ceiling
    from out/rt_edge/verdict.md, cross-checked against the gate definition
    in docs/rt_edge_contract.md, and the collected-flux floor with the tag
    the verdict prints next to it. A gate the files do not state comes back
    with value null and a note."""
    vpath = OUT / GATE_FILES["verdict"]
    vrel = f"out/{GATE_FILES['verdict']}"
    crel = GATE_FILES["contract"].relative_to(ROOT).as_posix()
    m, ln = _first_match(vpath, _G2_GATE_RE)
    mc, lnc = _first_match(GATE_FILES["contract"], _G2_CONTRACT_RE)
    g2 = {"value": None, "tag": "A", "unit": "",
          "tag_source": "declared: acceptance-gate convention (no bracket tag printed next to it)",
          "source": None, "definition": None}
    if m:
        g2["value"] = float(m.group(1))
        g2["source"] = f"{vrel}:{ln}"
        g2["definition"] = "pulsed intrinsic g2(0) < " + m.group(1)
    if mc:
        g2["contract_source"] = f"{crel}:{lnc}"
        if m and float(mc.group(1)) != g2["value"]:
            g2["note"] = (f"{vrel} says {m.group(1)} but {crel} says {mc.group(1)}; value withheld")
            g2["value"] = None
    if not m:
        g2["note"] = f"no 'pulsed intrinsic g2(0) < x' line in {vrel}"
    mf, lnf = _first_match(vpath, _FLUX_FLOOR_RE)
    flux = {"value": None, "unit": "photons/s", "tag": "A", "source": None}
    if mf:
        tags = [t for t in mf.group(1).split("/") if t in ("V", "DR", "E", "A")]
        order = ("V", "DR", "E", "A")
        flux.update({"value": float(mf.group(2)), "unit": mf.group(3),
                     "tag": max(tags, key=order.index) if tags else "A",
                     "tag_source": f"bracket tag printed in {vrel}",
                     "source": f"{vrel}:{lnf}", "definition": "pulsed collected-flux eligibility floor"})
    else:
        flux["note"] = f"no 'eligibility floor [tag]: x photons/s' line in {vrel}"
    return {"g2_ceiling": g2, "flux_floor": flux}


# ----------------------------------------------------- card -> verdict (I6)

_HEADLINE_MODEL_TOKEN = "finite_pulse:true,tau_cap_density:false"


def _same_number(text, value) -> bool:
    try:
        return float(text) == float(value)
    except (TypeError, ValueError):
        return False


def card_verdict(name: str) -> dict:
    """The committed VERDICT line(s) that speak for one design card. Matching
    is exact on what the card itself carries -- regime (drive.cycle_loading),
    family / orientation (nanowire.family; planar dot.orientation, nonpolar
    a_plane -> family a_plane), strain bound and rep rate (nanowire) -- and
    never picks among several matches: match is "exact" (one line),
    "ambiguous" (N lines, all returned) or "none" (reason says why)."""
    from . import api_cards
    from .api_cards import load_design
    # A Studio-made design (the unsaved draft route, or a card saved into cards/studio/) is not a
    # card any campaign swept: matching its numbers to a committed VERDICT line would be a claim
    # about a design that was never run there. Refuse with a reason, never map it.
    if name in api_cards.RESERVED_NAMES:
        return {"campaign": None, "card": name, "lines": [], "match": "none", "criteria": {}, "refused": "unsaved",
                "reason": "An unsaved Studio design: no committed VERDICT line speaks for it "
                          "(the campaigns swept the shipped cards, not this one)."}
    if api_cards.NAME_RE.match(name or "") and not name.endswith(".yaml") \
            and not (api_cards.CARDS / f"{name}.yaml").is_file() \
            and (api_cards.USER_CARDS / f"{name}.yaml").is_file():
        return {"campaign": None, "card": name, "lines": [], "match": "none", "criteria": {}, "refused": "saved",
                "reason": "A design saved in the Studio: the committed VERDICT lines were computed for the "
                          "shipped cards, not for this one."}
    d = load_design(name)
    regime = getattr(d.drive, "cycle_loading", None)
    crit: dict = {}
    reason = ""
    if d.platform == "legacy":
        if getattr(d.emission, "type", "none") != "edge":
            return {"campaign": None, "lines": [], "match": "none", "criteria": {},
                    "reason": "legacy (non-edge) cards are not swept by any committed campaign with VERDICT lines"}
        cid = "rt_edge"
        crit = {"model": _HEADLINE_MODEL_TOKEN}
        reason = ("rt_edge sweeps both edge cards; its VERDICT line is the pooled headline-model "
                  "verdict (per-card detail is in the CARD lines of verdict.md)")
    elif d.platform == "ingan_gan_planar":
        dot = (d.nitride or {}).get("dot", {}) if isinstance(d.nitride, dict) else {}
        if dot.get("geometry_type") == "qw_fluctuation":
            return {"campaign": "nitride_geometry_stark", "lines": [], "match": "none",
                    "criteria": {"geometry_type": "qw_fluctuation"},
                    "reason": "QW-fluctuation geometry is a supplementary sweep (geometry_supplement.csv) "
                              "with no VERDICT line of its own"}
        if (d.provenance or {}).get("comparison_status") or "comparison" in name:
            return {"campaign": None, "lines": [], "match": "none", "criteria": {},
                    "reason": "non-gating comparison card (design.provenance.comparison_status); "
                              "no campaign VERDICT line speaks for it"}
        orient = dot.get("orientation") or "c_plane"
        orient = "a_plane" if orient in ("nonpolar", "a_plane") else orient
        if orient == "c_plane":
            cid = "nitride_cavity"
            crit = {"regime": regime}
            reason = "c-plane planar cavity cards are the nitride_cavity campaign's subject (regime match)"
        else:
            cid = "nitride_geometry_stark"
            crit = {"family": orient, "regime": regime}
            reason = f"{orient} orientation: nitride_geometry_stark family {orient}, regime match"
    elif d.platform == "ingan_gan_nanowire":
        nw = (d.nitride or {}).get("nanowire", {}) if isinstance(d.nitride, dict) else {}
        cid = "nitride_nanowire/full"
        crit = {"family": nw.get("family"), "regime": regime,
                "strain_bound": nw.get("strain_bound") or "relaxed",
                "rep_rate_hz": float(d.drive.rep_rate_hz)}
        reason = "nanowire: exact family, regime, strain bound and rep rate from the card"
    else:
        return {"campaign": None, "lines": [], "match": "none", "criteria": {},
                "reason": f"no campaign for platform {d.platform!r}"}
    md = _md_file(OUT / cid)
    verdicts, _ = parse_md(md)
    hits = []
    for v in verdicts:
        f = v["fields"]
        ok = True
        for k, want in crit.items():
            have = f.get(k)
            if isinstance(want, float):
                ok = _same_number(have, want)
            else:
                ok = have == want
            if not ok:
                break
        if ok:
            hits.append(v)
    n = len(hits)
    match = "exact" if n == 1 else ("ambiguous" if n > 1 else "none")
    if n > 1:
        reason += f"; ambiguous: {n} lines match (not picked)"
        differ = sorted({k for v in hits for k in v["fields"]
                         if len({h["fields"].get(k) for h in hits}) > 1})
        reason += f"; they differ in {differ}"
    elif n == 0:
        reason += f"; no VERDICT line in out/{cid} matches {crit}"
    return {"campaign": cid, "card": name, "md": _rel(md) if md else None, "criteria": crit,
            "lines": hits, "match": match, "n_matches": n, "reason": reason}
