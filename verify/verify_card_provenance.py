#!/usr/bin/env python
"""Deterministic provenance checker for the paper cards (no-hallucination audit step H1).

For every cards/papers/<id>.yaml (files starting with `_` are skipped) each leaf / reported value / series
that carries tag V (verbatim) or DR (derived from paper values) and a `where:` string is checked against the
text of the paper(s) mapped to that card in verify/data/provenance_sources.yaml:

  (a) NUMBER   V only: the value (or the value as printed in the `where:` text) must occur in the paper text,
               allowing format variants (0.95 / .95 / 95 % ; 6.5e6 / 6.5 x 10^6 ; ns<->ps, GHz<->MHz, ueV<->meV ;
               the pypdf artefact "¼" for "=" ; +- / ±) and unit-scaled forms (x1e3, x1e6, x100, ...).
  (b) QUOTE    every quoted span (single or double quotes) of >= 6 words in `where:` must occur in the text after
               normalisation; fuzzy sliding-window ratio >= 0.9 (difflib), best ratio reported. Text after a colon
               is checked too but is informational unless --strict-colon (it is usually a paraphrase).
  (c) DR       every number written in the `where:` derivation string (page / figure refs, small integers, unit
               constants and the item's own value removed) must occur in the paper text. Arithmetic is NOT redone.

Tags E and A are skipped and counted. Verdicts: FOUND / NOT_FOUND / NO_SOURCE / NOT_CHECKABLE.

Also checks verify/data/reproduce_ledger.yaml `paper:` blocks and the [V]/[DR] comments of
fsim_core/sde.py::QDLaserParams.zhao_table_s1() against the Zhao supplement.

Usage:
    python verify/verify_card_provenance.py                 # full audit, writes the report, exit 0
    python verify/verify_card_provenance.py --self-test     # planted-error self-test, exit 1 on failure
    python verify/verify_card_provenance.py --only cao2015 [--no-ledger] [--no-zhao] [--strict-colon]
"""
from __future__ import annotations

import argparse
import bisect
import difflib
import glob
import math
import os
import re
import sys
import unicodedata
from collections import Counter, OrderedDict

import yaml

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CARDS_DIR = os.path.join(ROOT, "cards", "papers")
MAP_PATH = os.path.join(ROOT, "verify", "data", "provenance_sources.yaml")
LEDGER_PATH = os.path.join(ROOT, "verify", "data", "reproduce_ledger.yaml")
SDE_PATH = os.path.join(ROOT, "fsim_core", "sde.py")
REPORT_PATH = os.path.join(ROOT, ".workers", "devlib", "r4-provenance-report.md")
PDF_CACHE = os.path.join(ROOT, ".workers", "devlib", "r4-provenance-txt")

QUOTE_MIN_WORDS = 6
QUOTE_THRESHOLD = 0.9
SCALES = [1e3, 1e-3, 1e6, 1e-6, 1e9, 1e-9, 1e12, 1e-12, 100.0, 0.01]   # unit / percent rescalings
DECADES = [10.0 ** k for k in range(-15, 16)]                            # mantissa-only matching (Zhao SI)
CONSTANTS = [1239.84, 1239.842, 1239.8, 12398.4, 1.2398, 4.1357, 4.135667696, 0.6582, 0.658212, 2.99792458,
             1.380649, 8.617333, 8.617, 1.054571817, 6.62607015, 1.602176634, 9.1093837015, 3.14159, 6.022]


# ---------------------------------------------------------------------------------------------- text prep
_SUP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺", "0123456789-+")
_PAGE_RE = re.compile(r"^=====PAGE (\d+)=====\s*$")


def prep(text: str) -> str:
    """Normalise raw pypdf text for matching. Newlines are preserved one-for-one (line numbers stay valid)."""
    text = re.sub(r"[⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+", lambda m: "^" + m.group(0).translate(_SUP), text)
    text = text.replace("¼", "=").replace("\u00ad", "")          # pypdf artefact for "=", soft hyphen
    text = unicodedata.normalize("NFKC", text)                   # ligatures fi fl ffi, fullwidth, mu sign
    for ch in "\u2212\u2013\u2014\u2010\u2011\u2012\u2015":
        text = text.replace(ch, "-")
    for ch in "\u2018\u2019\u201a\u201b\u2032\u00b4`":
        text = text.replace(ch, "'")
    for ch in "\u201c\u201d\u201e\u201f\u2033":
        text = text.replace(ch, '"')
    text = text.replace("\u03bc", "u").replace("\u00b5", "u")    # mu -> u (ueV)
    text = text.replace("\u223c", "~").replace("\u2217", "*")
    text = text.replace("\u00b1", "+-")
    return text


def tokenise(s: str) -> list:
    return re.findall(r"[a-z0-9]+", s.lower())


# ---------------------------------------------------------------------------------------------- numbers
NUM_RE = re.compile(
    r"""(?<![A-Za-z0-9_.])
    (?:
      (?P<sci>(?P<m>\d+(?:\.\d+)?|\.\d+)\s*(?:[x\u00d7\u00b7*]|\\times|\\cdot)\s*10\s*\^?\s*\{?\s*(?P<e>[-+]?\d+)\s*\}?)
     |(?P<e2>(?:\d+(?:\.\d+)?|\.\d+)[eE][-+]?\d+)
     |(?P<num>\d{1,3}(?:,\d{3})+(?:\.\d+)?|\d+(?:\.\d+)?|\.\d+)
    )""", re.X)


def sig_digits(tok: str) -> int:
    t = tok.replace(",", "").lstrip("0.").replace(".", "")
    if "." not in tok:
        t = t.rstrip("0") or "0"
    return len(t.strip("0")) if "." not in tok else len(t)


def key6(x: float) -> str:
    x = abs(x)
    return "0" if x == 0 else f"{x:.6g}"


def mant_key(x: float) -> str:
    x = abs(x)
    if x == 0:
        return "0"
    return f"{x / 10 ** math.floor(math.log10(x)):.5g}"


def parse_numbers(s: str, with_tokens=False):
    """All numbers in s (sign dropped). Returns list of (value, raw_token, start)."""
    out = []
    for m in NUM_RE.finditer(s):
        try:
            if m.group("sci"):
                v = float(m.group("m")) * 10.0 ** int(m.group("e"))
                out.append((v, m.group(0), m.start()))
                mm = m.group("m")
                out.append((float(mm), mm, m.start()))       # the bare mantissa too
            elif m.group("e2"):
                out.append((float(m.group("e2")), m.group("e2"), m.start()))
            else:
                tok = m.group("num")
                out.append((float(tok.replace(",", "")), tok, m.start()))
        except ValueError:
            pass
    return out


_REF_RES = [
    re.compile(r"\bpp?\.\s*\d+(?:\s*[-/,&]\s*\d+)*", re.I),
    re.compile(r"\bpages?\s+\d+(?:\s*[-/,&]\s*\d+)*", re.I),
    re.compile(r"\b(?:Fig(?:ure)?s?|Tables?|Tab|Sec(?:tion)?s?|Eqs?|Refs?|SI|Supplement\w*|Note|Panel)\.?\s*S?\d+[A-Za-z]?"
               r"(?:\([a-z]\))?(?:\s*[-/,&]\s*(?:and\s*)?S?\d+[A-Za-z]?(?:\([a-z]\))?)*", re.I),
    re.compile(r"arXiv:\s*\d+\.\d+(?:v\d)?", re.I),
    re.compile(r"doi:\s*\S+", re.I),
    re.compile(r"\[\d+(?:[,-]\d+)*\]"),
    re.compile(r"\bv\d\b"),
]


def strip_refs(s: str) -> str:
    for r in _REF_RES:
        s = r.sub(" ", s)
    return s


def where_numbers(where: str) -> list:
    """Numbers written in a where: string, minus page/figure references and year-like integers."""
    s = prep(strip_refs(where))
    out = []
    for v, tok, _ in parse_numbers(s):
        if re.fullmatch(r"(?:19|20)\d\d", tok):
            continue
        out.append((v, tok))
    return out


def value_numbers(v) -> list:
    """Floats carried by a card value (number, range list, string with numbers)."""
    if isinstance(v, bool) or v is None:
        return []
    if isinstance(v, (int, float)):
        return [float(v)]
    if isinstance(v, str):
        return [x[0] for x in parse_numbers(prep(v))]
    if isinstance(v, (list, tuple)):
        out = []
        for e in v:
            out.extend(value_numbers(e))
        return out
    return []


# ---------------------------------------------------------------------------------------------- source text
class Source:
    """One paper text file: number index (with page/line) and token index for quote matching."""

    def __init__(self, label: str, text: str):
        self.label = label
        self.raw = text
        self.P = prep(text)
        lines = self.P.split("\n")
        self.line_page = []
        pg = 0
        for ln in lines:
            m = _PAGE_RE.match(ln.strip())
            if m:
                pg = int(m.group(1))
            self.line_page.append(pg)
        self.nl = [m.start() for m in re.finditer("\n", self.P)]
        self.nums = {}                      # key6 -> [(value, token, page, line)]
        self.mants = {}                     # mant_key -> [(value, token, page, line)]
        self.allvals = []                   # [(value, token, page, line)]
        for v, tok, st in parse_numbers(self.P):
            ln = bisect.bisect_right(self.nl, st)               # 0-based line index
            rec = (v, tok, self.line_page[min(ln, len(lines) - 1)], ln + 1)
            self.nums.setdefault(key6(v), []).append(rec)
            self.mants.setdefault(mant_key(v), []).append(rec)
            self.allvals.append(rec)
        # token variants for quotes: V1 joins line-break hyphenation, V2 splits it
        self.variants = []
        segs = []                           # (page, text)
        cur_pg, buf = 0, []
        for ln in lines:
            m = _PAGE_RE.match(ln.strip())
            if m:
                segs.append((cur_pg, "\n".join(buf)))
                cur_pg, buf = int(m.group(1)), []
            else:
                buf.append(ln)
        segs.append((cur_pg, "\n".join(buf)))
        for join in (True, False):
            toks, pages = [], []
            for pg, t in segs:
                t = t.lower()
                t = re.sub(r"-[ \t]*\n[ \t]*(?=[a-z0-9])", "" if join else " ", t)
                tk = tokenise(t)
                toks.extend(tk)
                pages.extend([pg] * len(tk))
            compact = "".join(toks)
            cum, acc = [], 0
            for tk in toks:
                cum.append(acc)
                acc += len(tk)
            self.variants.append({"toks": toks, "pages": pages, "compact": compact, "cum": cum,
                                  "tokset": None})

    # ---- numbers
    def find_number(self, x: float, scaled=True, decade=False):
        """Return (record, scale, kind) for the first occurrence of x (exact, then unit-scaled), else None."""
        if decade:
            recs = self.mants.get(mant_key(x))
            if recs:
                return recs[0], None, "decade"
            return None
        recs = self.nums.get(key6(x))
        if recs:
            return recs[0], 1.0, "exact"
        if scaled:
            for s in SCALES:
                recs = [r for r in self.nums.get(key6(x * s), []) if sig_digits(r[1]) >= 2 or "10" in r[1]]
                if recs:
                    return recs[0], s, "scaled"
        return None

    def nearest(self, x: float):
        best = None
        for rec in self.allvals:
            if rec[0] == 0 or x == 0:
                continue
            d = abs(math.log(abs(rec[0]) / abs(x)))
            if best is None or d < best[0]:
                best = (d, rec)
        return best[1] if best else None

    # ---- quotes
    def best_quote(self, q: str):
        """Best fuzzy ratio of quote q against the text. Returns (ratio, page, snippet)."""
        qt = tokenise(prep(q))
        if not qt:
            return 0.0, 0, ""
        qcomp = "".join(qt)
        best = (0.0, 0, "")
        for var in self.variants:
            idx = var["compact"].find(qcomp)
            if idx >= 0:
                ti = max(0, bisect.bisect_right(var["cum"], idx) - 1)
                return 1.0, var["pages"][ti], " ".join(var["toks"][ti:ti + len(qt)])
            r = self._fuzzy(var, qt)
            if r[0] > best[0]:
                best = r
        return best

    @staticmethod
    def _fuzzy(var, qt):
        toks, n = var["toks"], len(qt)
        N = len(toks)
        if N < 2:
            return (0.0, 0, "")
        qset = set(qt)
        pre = [0]
        for t in toks:
            pre.append(pre[-1] + (1 if t in qset else 0))
        w = min(n + 2, N)
        cands = []
        for s in range(0, N - w + 1):
            sc = pre[s + w] - pre[s]
            if sc >= 0.45 * n:
                cands.append((sc, s))
        cands.sort(reverse=True)
        qs = " ".join(qt)
        best = (0.0, 0, "")
        seen = set()
        for sc, s in cands[:10]:
            for s0 in (s - 1, s, s + 1, s + 2):
                for ln in (n - 1, n, n + 1, n + 2):
                    if s0 < 0 or ln < 2 or s0 + ln > N or (s0, ln) in seen:
                        continue
                    seen.add((s0, ln))
                    cand = " ".join(toks[s0:s0 + ln])
                    sm = difflib.SequenceMatcher(None, qs, cand, autojunk=False)
                    if sm.real_quick_ratio() < best[0] or sm.quick_ratio() < best[0]:
                        continue
                    r = sm.ratio()
                    if r > best[0]:
                        best = (r, var["pages"][s0], cand)
        return best


_SRC_CACHE = {}


def resolve_text_path(rel: str):
    """Path (repo-relative or absolute) of a text source; a .pdf is extracted with pypdf into PDF_CACHE."""
    p = rel if os.path.isabs(rel) else os.path.normpath(os.path.join(ROOT, rel))
    if p.lower().endswith(".pdf"):
        os.makedirs(PDF_CACHE, exist_ok=True)
        out = os.path.join(PDF_CACHE, re.sub(r"[^A-Za-z0-9_.-]+", "_", os.path.basename(p)[:-4]) + ".txt")
        if not os.path.exists(out):
            from pypdf import PdfReader
            rd = PdfReader(p)
            parts = []
            for i, pg in enumerate(rd.pages, 1):
                parts.append(f"=====PAGE {i}=====\n" + (pg.extract_text() or ""))
            with open(out, "w", encoding="utf-8") as f:
                f.write("\n".join(parts))
        return out
    return p


def load_source(rel: str):
    if rel in _SRC_CACHE:
        return _SRC_CACHE[rel]
    try:
        p = resolve_text_path(rel)
        with open(p, encoding="utf-8", errors="replace") as f:
            text = f.read()
        src = Source(os.path.basename(p), text) if text.strip() else None
    except Exception as e:                                       # missing / unreadable
        print(f"  [warn] cannot load source {rel}: {e}", file=sys.stderr)
        src = None
    _SRC_CACHE[rel] = src
    return src


_ALL_CORPUS = None


def all_corpus():
    """Every candidate paper text (used only for 'found elsewhere' hints on NOT_FOUND numbers)."""
    global _ALL_CORPUS
    if _ALL_CORPUS is None:
        paths = sorted(glob.glob(os.path.join(ROOT, ".workers", "devlib", "papers", "raw", "*", "*.txt")) +
                       glob.glob(os.path.join(ROOT, "..", "refs", "**", "*.txt"), recursive=True))
        _ALL_CORPUS = [s for s in (load_source(os.path.normpath(p)) for p in paths) if s]
    return _ALL_CORPUS


# ---------------------------------------------------------------------------------------------- item checks
class Result:
    def __init__(self, path, tag, value, where):
        self.path, self.tag, self.value, self.where = path, tag, value, where
        self.verdict = "NOT_CHECKABLE"
        self.num = "n/a"                # text for the number column
        self.quote = "-"
        self.notes = []
        self.fails = []
        self.best = ""                  # best match text for NOT_FOUND reporting
        self.n_colon_low = 0


def extract_quotes(where: str):
    """[(text, kind)] with kind 'quote' (inside '...' / "...") or 'colon' (text after a colon)."""
    w = prep(where)
    spans = []
    rest = w
    for rx in (re.compile(r'"([^"]{12,}?)"'), re.compile(r"(?<![A-Za-z0-9])'([^']{12,}?)'(?![A-Za-z0-9])")):
        for m in rx.finditer(w):
            spans.append((m.group(1).strip(), "quote"))
        rest = rx.sub(" ", rest)
    out = [(t, k) for t, k in spans if len(re.findall(r"[A-Za-z]+", t)) >= 1 and len(t.split()) >= QUOTE_MIN_WORDS]
    if ":" in rest:
        seg = rest.split(":", 1)[1].split(";")[0].strip()
        if len(re.findall(r"[A-Za-z]{2,}", seg)) >= QUOTE_MIN_WORDS and len(seg.split()) >= QUOTE_MIN_WORDS:
            out.append((seg, "colon"))
    return out


def _equal_under_scale(w: float, v: float) -> bool:
    if key6(w) == key6(v):
        return True
    return any(key6(w * s) == key6(v) for s in SCALES)


def find_any(srcs, x: float, decade=False):
    """Exact occurrence in any source first, then unit-scaled. Returns (src, rec, scale, kind) or None."""
    if decade:
        for s in srcs:
            r = s.find_number(x, decade=True)
            if r:
                return s, r[0], None, "decade"
        return None
    for s in srcs:
        r = s.find_number(x, scaled=False)
        if r:
            return s, r[0], 1.0, "exact"
    for s in srcs:
        r = s.find_number(x, scaled=True)
        if r:
            return s, r[0], r[1], "scaled"
    return None


def loc(h) -> str:
    s, rec, sc, kind = h
    base = f"{s.label} p{rec[2]} L{rec[3]} '{rec[1].strip()}'"
    if kind == "scaled":
        base += f" (x{sc:g})"
    if kind == "decade":
        base += " (decade)"
    if kind == "exact" and rec[0] != 0 and sig_digits(rec[1]) <= 1:
        base += " [weak: 1 digit]"
    return base


def check_item(it: dict, srcs: list, strict_colon=False, decade=False) -> Result:
    """it: {path, tag, value(list of floats), where, digitized, table, label}. Returns Result."""
    tag = it["tag"]
    res = Result(it["path"], tag, it.get("label", ""), it.get("where", ""))
    if not srcs:
        res.verdict = "NO_SOURCE"
        res.num = "no source"
        return res
    where = it.get("where") or ""
    checked = False
    vals = it.get("values") or []
    own = list(vals)

    # (a) NUMBER, V only
    if tag == "V" and vals and not it.get("digitized") and not it.get("table"):
        wn = where_numbers(where)
        locs, miss = [], []
        for v in vals:
            printed = [w for w, _ in wn if _equal_under_scale(w, v)]
            hit = None
            for c in [v] + printed:
                hit = find_any(srcs, c, decade)
                if hit:
                    break
            if hit:
                locs.append(loc(hit))
            else:
                miss.append(v)
        checked = True
        if miss:
            near = [srcs[0].nearest(m) for m in miss]
            res.num = "NO"
            nr = "; ".join(f"{m:g} nearest {n[1]} p{n[2]}" for m, n in zip(miss, near) if n)
            res.notes.append("missing " + ", ".join(f"{m:g}" for m in miss) + (f" ({nr})" if nr else ""))
            res.fails.append("number")
            for m in miss:                                     # hint: present in a different paper?
                for cs in all_corpus():
                    if cs not in srcs and cs.find_number(m, scaled=False) and sig_digits(f"{m:g}") >= 2:
                        res.notes.append(f"{m:g} occurs in OTHER paper {cs.label}")
                        break
            # keep the first nearest as 'best match'
            if near and near[0]:
                res.best = f"nearest number {near[0][1]} (p{near[0][2]} L{near[0][3]})"
        else:
            res.num = "yes: " + " | ".join(locs[:3])
    elif tag == "V" and (it.get("digitized") or it.get("table")):
        res.num = "n/a (table/digitised)"
    elif tag == "DR":
        wn = where_numbers(where)
        keep = []
        for w, tok in wn:
            if "." not in tok and w == int(w) and abs(w) <= 10:
                continue
            if any(mant_key(w) == mant_key(c) for c in CONSTANTS):
                continue
            if any(_equal_under_scale(w, o) or key6(w) == key6(o) for o in own):
                continue
            keep.append((w, tok))
        found, missing = [], []
        seen = set()
        for w, tok in keep:
            if key6(w) in seen:
                continue
            seen.add(key6(w))
            h = find_any(srcs, w, decade)
            (found if h else missing).append((w, tok, h))
        if keep:
            checked = True
            res.num = f"DR nums {len(found)}/{len(found) + len(missing)}"
            if found:
                res.notes.append("found " + ", ".join(f"{t.strip()}@{loc(h).split(' ')[1]}" for _, t, h in found[:6]))
            if missing:
                res.notes.append("not found " + ", ".join(t.strip() for _, t, _ in missing))
                res.fails.append("DR-numbers")
                res.best = "DR numbers not in text: " + ", ".join(t.strip() for _, t, _ in missing)
        else:
            res.num = "DR (no paper numbers in where)"

    # (b) QUOTES
    qs = extract_quotes(where)
    ratios = []
    for q, kind in qs:
        best = (0.0, 0, "")
        for s in srcs:
            r = s.best_quote(q)
            if r[0] > best[0]:
                best = (r[0], r[1], r[2], s.label) if len(r) == 3 else r
        r, pg, snip = best[0], best[1], best[2]
        ratios.append((kind, r, pg))
        if kind == "quote":
            checked = True
            if r < QUOTE_THRESHOLD:
                res.fails.append("quote")
                res.notes.append(f"quote ratio {r:.2f}: '{q[:70]}'")
                res.best = (res.best + " | " if res.best else "") + f"best {r:.2f} p{pg}: '{snip[:90]}'"
        else:
            if r < QUOTE_THRESHOLD:
                res.n_colon_low += 1
                if strict_colon:
                    checked = True
                    res.fails.append("colon-span")
                    res.best = (res.best + " | " if res.best else "") + f"colon best {r:.2f} p{pg}: '{snip[:90]}'"
            elif strict_colon:
                checked = True
    if ratios:
        res.quote = " ".join(("q" if k == "quote" else "c") + f"{r:.2f}" for k, r, _ in ratios)

    if res.fails:
        res.verdict = "NOT_FOUND"
    elif checked:
        res.verdict = "FOUND"
    else:
        res.verdict = "NOT_CHECKABLE"
    return res


# ---------------------------------------------------------------------------------------------- card walking
def walk_card(card: dict):
    """Yield item dicts for every node with a tag, plus counts of skipped E / A."""
    items, skipped = [], Counter()

    def label_of(d):
        if "value" in d:
            v = d["value"]
        elif "range" in d:
            v = d["range"]
        elif "points" in d:
            v = f"<{len(d['points'])} points>"
        else:
            v = ""
        s = str(v)
        if "err" in d:
            s += f" ± {d['err']}"
        return s if len(s) <= 40 else s[:37] + "..."

    def rec(node, path):
        if isinstance(node, dict):
            if "tag" in node and isinstance(node.get("tag"), str) and node["tag"] in ("V", "DR", "E", "A"):
                tag = node["tag"]
                if tag in ("E", "A"):
                    skipped[tag] += 1
                else:
                    where = node.get("where")
                    note = ""
                    if not where:
                        where = node.get("source", "") or ""
                        note = "no where:, used source"
                    base = {"path": path, "tag": tag, "where": where, "label": label_of(node),
                            "digitized": bool(node.get("digitized")), "note": note}
                    table = "points" in node or (isinstance(node.get("value"), (list, tuple)) and
                                                 len(value_numbers(node.get("value"))) > 12)
                    vals = []
                    if "value" in node:
                        vals += value_numbers(node["value"])
                    if "range" in node:
                        vals += value_numbers(node["range"])
                    items.append(dict(base, values=vals, table=table))
                    if "err" in node and node.get("err") is not None and not table:
                        ev = value_numbers(node["err"])
                        if ev and tag == "V":
                            items.append(dict(base, path=path + ":err", values=ev, table=False,
                                              label=f"err {node['err']}"))
            for k, v in node.items():
                if k in ("history", "meta"):
                    continue
                if isinstance(v, (dict, list)):
                    rec(v, f"{path}.{k}" if path else k)
        elif isinstance(node, list):
            for i, v in enumerate(node):
                nm = v.get("id") if isinstance(v, dict) and isinstance(v.get("id"), str) else str(i)
                rec(v, f"{path}[{nm}]")

    rec(card, "")
    return items, skipped


# ---------------------------------------------------------------------------------------------- reporting
def esc(s) -> str:
    return str(s).replace("|", "\\|").replace("\n", " ")


def card_table(rows) -> str:
    out = ["| item path | tag | value | number found? (where: page/line) | quote ratio | verdict | notes |",
           "|---|---|---|---|---|---|---|"]
    for r in rows:
        out.append(f"| {esc(r.path)} | {r.tag} | {esc(r.value)} | {esc(r.num)} | {esc(r.quote)} | {r.verdict} | "
                   f"{esc('; '.join(r.notes))[:240]} |")
    return "\n".join(out)


def load_map():
    with open(MAP_PATH, encoding="utf-8") as f:
        return yaml.safe_load(f) or {}


def card_files():
    return sorted(f for f in glob.glob(os.path.join(CARDS_DIR, "*.yaml")) if not os.path.basename(f).startswith("_"))


def audit_cards(only=None, strict_colon=False):
    m = load_map()
    cmap = m.get("cards", {}) or {}
    per_card, totals, unmapped, skipped_total = OrderedDict(), Counter(), [], Counter()
    colon_low = 0
    for f in card_files():
        cid = os.path.basename(f)[:-5]
        if only and cid != only:
            continue
        with open(f, encoding="utf-8") as fh:
            card = yaml.safe_load(fh)
        items, skipped = walk_card(card)
        skipped_total.update(skipped)
        paths = cmap.get(cid) or []
        srcs = [s for s in (load_source(p) for p in paths) if s]
        if not srcs:
            unmapped.append(cid)
        rows = []
        for it in items:
            r = check_item(it, srcs, strict_colon)
            if it.get("note"):
                r.notes.append(it["note"])
            rows.append(r)
            totals[r.verdict] += 1
            colon_low += r.n_colon_low
        per_card[cid] = (rows, skipped, [s.label for s in srcs])
    return per_card, totals, unmapped, skipped_total, colon_low


# ---------------------------------------------------------------------------------------------- ledger / Zhao
def audit_ledger(strict_colon=False):
    m = load_map()
    kws = m.get("ledger_keywords", []) or []
    with open(LEDGER_PATH, encoding="utf-8") as f:
        led = yaml.safe_load(f)
    rows, nopaper = [], 0
    for it in led.get("items", []):
        p = it.get("paper")
        if not isinstance(p, dict):
            nopaper += 1
            continue
        tag = p.get("tag")
        if tag not in ("V", "DR"):
            continue
        where = p.get("where", "") or ""
        files = []
        for ent in kws:
            if re.search(ent["match"], where + " " + str(it.get("desc", "")), re.I):
                files = ent.get("files") or []
                break
        srcs = [s for s in (load_source(x) for x in files) if s]
        vals = value_numbers(p.get("value"))
        d = {"path": f"ledger:{it['id']}", "tag": tag, "where": where, "label": str(p.get("value")),
             "values": vals, "digitized": False, "table": False}
        rows.append(check_item(d, srcs, strict_colon))
        if p.get("err") is not None and tag == "V":
            d2 = dict(d, path=f"ledger:{it['id']}:err", values=value_numbers(p["err"]), label=f"err {p['err']}")
            rows.append(check_item(d2, srcs, strict_colon))
    return rows, nopaper


def audit_zhao():
    m = load_map()
    zpath = m.get("zhao_supplement")
    src = load_source(zpath) if zpath else None
    with open(SDE_PATH, encoding="utf-8") as f:
        lines = f.read().split("\n")
    start = next(i for i, l in enumerate(lines) if "def zhao_table_s1" in l)
    end = next(i for i in range(start + 1, len(lines)) if lines[i].startswith("# ====="))
    rows = []
    for i in range(start, end):
        l = lines[i]
        mt = re.search(r"#\s*\[(V|DR)\b[^\]]*\](.*)$", l)
        if not mt or "=" not in l.split("#")[0] and "(" not in l.split("#")[0] and ":" not in l.split("#")[0]:
            continue
        code = l.split("#")[0]
        if not code.strip() or code.strip().startswith(("def ", '"""')):
            continue
        tag = mt.group(1)
        lit = [(v, t) for v, t, _ in parse_numbers(code) if not re.fullmatch(r"\d", t)] + \
              [(float(t), t) for t in re.findall(r"(?<![\w.])\d(?![\w.])", code)]
        name = re.sub(r"\s+", " ", code.strip())[:60]
        label = name
        if not lit:
            rows.append((i + 1, name, tag, "-", "NOT_CHECKABLE", "no numeric literal"))
            continue
        for v, t in lit:
            if src is None:
                rows.append((i + 1, name, tag, t, "NO_SOURCE", "supplement text unavailable"))
                continue
            h = find_any([src], v, decade=True)
            weak = " [weak]" if (sig_digits(t) <= 1) else ""
            if tag == "V":
                rows.append((i + 1, name, tag, t, "FOUND" if h else "NOT_FOUND",
                             (loc(h) + weak) if h else "no mantissa match in supplement"))
            else:
                rows.append((i + 1, name, tag, t, "FOUND" if h else "NOT_FOUND",
                             (loc(h) + weak) if h else "DR literal not printed in supplement (derived value)"))
    return rows


# ---------------------------------------------------------------------------------------------- self-test
def _first_text_words(raw: str, start: int, before: int, after: int) -> str:
    words = list(re.finditer(r"\S+", raw))
    ix = next(i for i, w in enumerate(words) if w.start() <= start < w.end() or w.start() > start)
    return " ".join(w.group(0) for w in words[max(0, ix - before): ix + after])


def run_self_test() -> int:
    raw = lambda cid: os.path.join(ROOT, ".workers", "devlib", "papers", "raw", cid, cid + ".txt")
    A = load_source(raw("cao2015"))
    B = load_source(raw("ellis2008"))
    assert A and B, "self-test texts missing"
    items = []                         # (item dict, srcs, expect, kind)

    def mk(path, tag, vals, where, label=None):
        return {"path": path, "tag": tag, "values": vals, "where": where,
                "label": label if label is not None else (f"{vals[0]:g}" if vals else ""), "digitized": False,
                "table": False}

    # ---- (iv) two correct items (percent form + a real quote broken across lines and with a ligature)
    items.append((mk("ok.percent_value", "V", [0.028], "Abstract (multiple photon emission probability reduced to 2.8%)"),
                  [A], "FOUND", "correct"))
    items.append((mk("ok.real_quote", "V", [1.324],
                     "text p.1: 'applying high frequency electrical pulses to a single quantum dot in a p-i-n diode "
                     "under resonant laser excitation' at E = 1.324 eV"), [A], "FOUND", "correct"))
    # ---- (i) three numbers not in the text: a real value with one digit changed
    for nm, real in (("fake.num_gate", 300.0), ("fake.num_efield", 1.324), ("fake.num_rep", 200.0)):
        s = f"{real:g}"
        bad = None
        for pos in range(len(s) - 1, -1, -1):
            if not s[pos].isdigit():
                continue
            for d in "97531":
                cand = s[:pos] + d + s[pos + 1:]
                if cand != s and cand[0] != "0" and not find_any([A], float(cand)):
                    bad = float(cand)
                    break
            if bad:
                break
        assert bad is not None
        items.append((mk(nm, "V", [bad], f"text p.7 ({bad:g} as printed)"), [A], "NOT_FOUND", "fake number"))
    # ---- (ii) three fabricated quotes with paper-like wording (number in each item IS in the text)
    fab = ["the multi photon emission probability was suppressed by a factor of four under pulsed resonant excitation of the exciton",
           "we find that the emission linewidth is limited by the radiative lifetime when the gate pulse is shorter than the Rabi period",
           "the resonance fluorescence count rate increased linearly with the repetition rate of the applied electrical pulse train"]
    for k, q in enumerate(fab):
        assert A.best_quote(q)[0] < 0.8, "fabricated quote accidentally present"
        items.append((mk(f"fake.quote_{k}", "V", [300.0], f'text p.8 / Fig 3c ("{q}")'), [A], "NOT_FOUND", "fake quote"))
    # ---- (iii) two real quotes from a DIFFERENT paper (ellis2008)
    bt = B.variants[0]["toks"]
    picked = []
    for frac in (0.30, 0.62):
        i = int(len(bt) * frac)
        for j in range(i, len(bt) - 14):
            win = bt[j:j + 12]
            if all(t.isalpha() and len(t) > 2 for t in win) and A.best_quote(" ".join(win))[0] < 0.8:
                picked.append(" ".join(win))
                break
    assert len(picked) == 2
    for k, q in enumerate(picked):
        assert B.best_quote(q)[0] >= 0.99
        items.append((mk(f"fake.otherpaper_{k}", "V", [300.0], f"text p.3 ('{q}')"), [A], "NOT_FOUND", "quote from other paper"))
    # ---- three true positives that need format normalisation, each against a real paper that prints the form
    corpus = [(c, load_source(raw(c))) for c in
              ("nowak2014", "nawrath2022", "schweickert2018", "tomm2021", "tamariz2020", "lazic2017")]
    tp = {}
    for cid, s in corpus:
        r = s.raw
        if "ueV" not in tp:
            m = re.search(r"\d[\d.]*\s*[μµ]eV", r)
            if m:
                span = _first_text_words(r, m.start(), 4, 6)
                tp["ueV"] = (cid, s, span)
        if "x10" not in tp:
            m = next((x for x in NUM_RE.finditer(s.P) if x.group("sci")), None)
            if m:
                tp["x10"] = (cid, s, (float(m.group("m")) * 10.0 ** int(m.group("e")), m.group(0)))
    for cid in ("nowak2014", "nawrath2022", "tomm2021", "karli2022", "wang2019", "uppu2020", "he2013", "thomas2021"):
        s = load_source(raw(cid))
        m = re.search(r"\S+\s+\S+\s+\S+\s+\S+\s+\S+\s+¼\s*\d[\d.]*\s*\S*\s*\S*", s.raw) if s else None
        if m:
            tp["eq"] = (cid, s, m.group(0))
            break
    assert {"ueV", "x10", "eq"} <= set(tp), f"format examples missing: {sorted(tp)}"
    cid, s, span = tp["ueV"]
    nval = float(re.search(r"(\d[\d.]*)\s*[μµ]eV", span).group(1))
    items.append((mk("tp.mu_eV", "V", [nval], f"text p.1: '{span.replace('μ', 'u').replace('µ', 'u')}'"),
                  [s], "FOUND", f"true positive mu-eV ({cid})"))
    cid, s, (v, tok) = tp["x10"]
    items.append((mk("tp.times10", "V", [v], f"text p.1 ({v:.4g} printed as '{tok.strip()}')"),
                  [s], "FOUND", f"true positive x10^ ({cid}: {tok.strip()})"))
    cid, s, span = tp["eq"]
    num = float(re.search(r"¼\s*(\d[\d.]*)", span).group(1))
    items.append((mk("tp.quarter_eq", "V", [num], "text p.1: '" + span.replace("¼", "=") + "'"),
                  [s], "FOUND", f"true positive ¼ for = ({cid})"))

    bad = 0
    print(f"self-test: {len(items)} planted items "
          f"(3 fake numbers, 3 fabricated quotes, 2 other-paper quotes, 2 correct, 3 format true positives)")
    for it, srcs, expect, kind in items:
        r = check_item(it, srcs)
        ok = r.verdict == expect
        bad += (not ok)
        print(f"  [{'ok ' if ok else 'BAD'}] {it['path']:<22} expect {expect:<9} got {r.verdict:<9} {kind}"
              f" | num={r.num[:50]} | quote={r.quote} | {'; '.join(r.notes)[:100]}")
    print(f"self-test {'PASSED' if not bad else 'FAILED'} ({len(items) - bad}/{len(items)} as expected)")
    return 1 if bad else 0


# ---------------------------------------------------------------------------------------------- main
def main(argv=None) -> int:
    try:
        sys.stdout.reconfigure(encoding="utf-8", errors="replace")
    except Exception:
        pass
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--self-test", action="store_true")
    ap.add_argument("--only", help="check a single card id")
    ap.add_argument("--no-ledger", action="store_true")
    ap.add_argument("--no-zhao", action="store_true")
    ap.add_argument("--strict-colon", action="store_true", help="count failing colon-spans as NOT_FOUND")
    ap.add_argument("--report", default=REPORT_PATH)
    a = ap.parse_args(argv)

    if a.self_test:
        return run_self_test()

    per_card, totals, unmapped, skipped_total, colon_low = audit_cards(a.only, a.strict_colon)
    led_rows, led_nopaper, zrows = [], 0, []
    led_tot, z_tot = Counter(), Counter()
    if not a.only and not a.no_ledger:
        led_rows, led_nopaper = audit_ledger(a.strict_colon)
        led_tot.update(r.verdict for r in led_rows)
    if not a.only and not a.no_zhao:
        try:
            zrows = audit_zhao()
            z_tot.update(r[4] for r in zrows)
        except Exception as e:
            print(f"[warn] Zhao check failed: {e}", file=sys.stderr)

    # ---- report
    L = ["# R4 provenance report (verify/verify_card_provenance.py)", "",
         "Deterministic check of every V/DR card value and `where:` quote against the paper text. Tags E and A are "
         "skipped. Verdicts: FOUND / NOT_FOUND / NO_SOURCE / NOT_CHECKABLE. Source map: "
         "`verify/data/provenance_sources.yaml`. Quote threshold 0.9 (difflib sliding window), colon-spans "
         f"informational{' (STRICT)' if a.strict_colon else ''}. Exit code 0 (audit tool, not yet a gate).", ""]
    L += ["## Totals (cards)", "", "| verdict | items |", "|---|---|"]
    for v in ("FOUND", "NOT_FOUND", "NO_SOURCE", "NOT_CHECKABLE"):
        L.append(f"| {v} | {totals.get(v, 0)} |")
    L.append(f"| skipped E | {skipped_total.get('E', 0)} |")
    L.append(f"| skipped A | {skipped_total.get('A', 0)} |")
    L += ["", f"Colon-spans below 0.9 (informational, mostly paraphrase): {colon_low}", ""]
    L += [f"Cards with no source text (NO_SOURCE): {', '.join(unmapped) if unmapped else 'none'}", ""]
    for cid, (rows, skipped, labels) in per_card.items():
        c = Counter(r.verdict for r in rows)
        L += [f"## {cid}", "", f"sources: {', '.join(labels) if labels else 'NONE'} | "
              + " ".join(f"{k}={v}" for k, v in sorted(c.items())) + f" | skipped E={skipped.get('E', 0)} A={skipped.get('A', 0)}", ""]
        L += [card_table(rows), ""]
    if led_rows:
        L += ["## reproduce_ledger.yaml `paper:` blocks", "",
              f"{len(led_rows)} V/DR rows checked, {led_nopaper} items without a `paper:` block. "
              + " ".join(f"{k}={v}" for k, v in sorted(led_tot.items())), "", card_table(led_rows), ""]
    if zrows:
        L += ["## fsim_core/sde.py zhao_table_s1() [V]/[DR] literals vs Zhao supplement (mantissa/decade match)", "",
              " ".join(f"{k}={v}" for k, v in sorted(z_tot.items())), "",
              "| sde.py line | statement | tag | literal | verdict | found at |", "|---|---|---|---|---|---|"]
        for ln, nm, tg, t, vd, where in zrows:
            L.append(f"| {ln} | {esc(nm)} | {tg} | {esc(t)} | {vd} | {esc(where)} |")
        L.append("")
    nf = [(cid, r) for cid, (rows, _, _) in per_card.items() for r in rows if r.verdict == "NOT_FOUND"]
    L += ["## NOT_FOUND list (cards)", "",
          "| card | item | value | where | best match |", "|---|---|---|---|---|"]
    for cid, r in nf:
        L.append(f"| {cid} | {esc(r.path)} | {esc(r.value)} | {esc(r.where)[:300]} | {esc(r.best or '; '.join(r.notes))[:300]} |")
    nf2 = [r for r in led_rows if r.verdict == "NOT_FOUND"]
    if nf2:
        L += ["", "## NOT_FOUND list (ledger)", "", "| item | value | where | best match |", "|---|---|---|---|"]
        for r in nf2:
            L.append(f"| {esc(r.path)} | {esc(r.value)} | {esc(r.where)[:200]} | {esc(r.best or '; '.join(r.notes))[:240]} |")
    status = ("STATUS: " + " ".join(f"{k}={totals.get(k, 0)}" for k in ("FOUND", "NOT_FOUND", "NO_SOURCE", "NOT_CHECKABLE"))
              + f" skippedE={skipped_total.get('E', 0)} skippedA={skipped_total.get('A', 0)}"
              + f" unmapped_cards={len(unmapped)} ledger[{' '.join(f'{k}={v}' for k, v in sorted(led_tot.items()))}]"
              + f" zhao[{' '.join(f'{k}={v}' for k, v in sorted(z_tot.items()))}]")
    L += ["", "## STATUS", "", "```", status, "```", ""]
    os.makedirs(os.path.dirname(a.report), exist_ok=True)
    with open(a.report, "w", encoding="utf-8") as f:
        f.write("\n".join(L))

    print(f"cards checked: {len(per_card)}  report: {a.report}")
    for v in ("FOUND", "NOT_FOUND", "NO_SOURCE", "NOT_CHECKABLE"):
        print(f"  {v:<14}{totals.get(v, 0)}")
    print(f"  skipped E={skipped_total.get('E', 0)} A={skipped_total.get('A', 0)}   colon-spans<0.9 (info)={colon_low}")
    print(f"  unmapped cards: {', '.join(unmapped) if unmapped else 'none'}")
    if led_rows:
        print(f"  ledger: {dict(led_tot)} ({led_nopaper} items without paper block)")
    if zrows:
        print(f"  zhao_table_s1: {dict(z_tot)}")
    print("NOT_FOUND:")
    for cid, r in nf:
        print(f"  {cid} | {r.path} | {r.value} | {r.where[:110]} | {(r.best or '; '.join(r.notes))[:140]}")
    for r in nf2:
        print(f"  ledger | {r.path} | {r.value} | {r.where[:90]} | {(r.best or '; '.join(r.notes))[:120]}")
    print(status)
    return 0


if __name__ == "__main__":
    sys.exit(main())
