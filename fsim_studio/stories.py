"""Story mode: load / save / bake story JSON files (fsim_studio/stories/).

A story is {name, title, scenes: [...]}; each scene carries {id, title,
claim, numbers: [{name, label, unit, tag_declared, source}], visual, caveat,
notes, catalog_id}. Numbers are NEVER typed in: each number's `source` names
an out/ file location or a live fsim_core call, and bake() resolves it,
writing value / raw / resolved_from (and the resolved tag) next to the
source. Values are copied exactly as the file prints them (no rounding; the
UI rounds).

Prose (title, claim, caveat, notes) is a template: `{name}` placeholders
refer to a number's `name` (story-wide namespace) and are interpolated from
the BAKED values on GET (interface I5); the stored file keeps the templates
and GET also returns them as `<field>_template`. Placeholder forms:
  {name}          value (verdict/best/md_regex: the printed text; else a
                  compact number; lists joined with ", ")
  {name[i]}       list element (negative i from the end)
  {name:.3g}      Python format spec applied to the numeric value
  {name.tag}      the number's resolved provenance tag
  {name.T}        the temperature a live / csv_cell source was evaluated at
  {name.say}      a status word in room language (STATUS_PHRASES, e.g.
                  pass_hardware_infeasible -> "passes optically, fails in
                  hardware"); raw tokens stay for the status chips

Tags are resolved at bake, from the source itself where the source carries
one: live -> the run's tag_chain (fsim_studio.tags.chain_of); md_regex ->
bracket tags inside the matched text; json -> the file's top-level
"tag_chain"; csv_cell / csv_agg -> a "<col>_tag" column. Where the source
prints no tag (VERDICT/BEST lines, data CSVs, commits) the author's
`tag_declared` stands and `tag_source` says so.

Source kinds (file paths are relative to the repo root and must lie under
out/ -- H13):
  verdict   {file, match?: {field: text}, index?: int, field | word: true}
  best      {file, best_kind?: "BEST_PASSING_FLUX", match?: {...}, field}
  md_regex  {file, pattern (<= 300 chars, no nested quantifiers), group?:
             int (default 1), flags?: "m"}
  json      {file, key: "a.b.c"}
  csv_cell  {file, match: {col: value}, col}
  csv_agg   {file, filter?: {col: value|[..]|{min,max}}, col, agg: count|min|max|values|unique}
  commit    {dir: "rt_edge"}            campaign data commit (manifest or git log)
  live      {card, mode: point|headline, overrides?: {path: value}, T_grid?: [..], scalar}

Staleness (H4): baked_hash covers the fsim_core source hash, every
out/**/manifest.json, every out/ file a number reads, the card file of every
live source, and every number's resolved fields (value, raw, resolved_from,
tag, tag_source). A value typed into the file, a moved out/ file or an edited
card all mark the story stale. PUT never accepts resolved fields: they are
stripped from the body and restored from the stored story for numbers whose
source did not change.
"""
from __future__ import annotations

import copy
import hashlib
import json
import math
import re
from pathlib import Path

from . import CARDS, OUT, PKG, ROOT, fsim_core_hash
from .serialize import to_jsonable
from .tags import chain_of, tags_in_text, widest

STORIES = PKG / "stories"
NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,80}$")
PROSE_FIELDS = ("title", "claim", "caveat", "notes")
RESOLVED_FIELDS = ("value", "raw", "resolved_from", "error", "live_skipped", "tag", "tag_source", "at_T")
STORY_STATE_FIELDS = ("current_hash", "stale", "baked", "baked_hash", "baked_fsim_core_hash", "baked_at",
                      "bake_errors", "prose_errors")
MAX_PATTERN = 300
# nested quantifier, e.g. (a+)+ / (a*)* / (a+){2,}: catastrophic-backtracking shape
_NESTED_QUANT = re.compile(r"\((?:[^()\\]|\\.)*[+*}](?:[^()\\]|\\.)*\)\s*[+*{]")
PLACEHOLDER = re.compile(r"\{([A-Za-z_]\w*)(?:\[(-?\d+)\])?(?:\.(tag|T|say))?(?::([^{}]*))?\}")
# Room language for printed status words ({name.say}); the raw token stays in the
# number's value/raw and renders only inside status chips. Wording, not physics.
STATUS_PHRASES = {
    "pass_hardware_infeasible": "passes optically, fails in hardware",
    "no_idealized_pass": "has no optical pass",
    "FAIL": "fails",
    "PASS": "passes",
}
RAW_TEXT_KINDS = ("verdict", "best", "md_regex", "commit")


class StoryError(Exception):
    def __init__(self, status: int, message: str):
        super().__init__(message)
        self.status = status
        self.message = message


def _path(name: str) -> Path:
    if not NAME_RE.match(name or ""):
        raise StoryError(400, f"invalid story name {name!r}")
    return STORIES / f"{name}.json"


def _repo_file(rel: str) -> Path:
    """A bake source file: must resolve inside out/ (H13: no .git, no
    sources, no cards through a story)."""
    p = (ROOT / str(rel)).resolve()
    if OUT.resolve() not in p.parents:
        raise ValueError(f"story sources must lie under out/: {rel!r}")
    if not p.is_file():
        raise FileNotFoundError(rel)
    return p


def _num_or_text(text):
    if isinstance(text, (int, float, bool)) or text is None:
        return text
    t = str(text).strip()
    try:
        f = float(t)
    except ValueError:
        return t
    return f if math.isfinite(f) else t


# ----------------------------------------------------------- resolvers
# Each returns (value, raw, where, extra); extra may carry tag / tag_source / at_T.

def _r_verdict(src):
    from .campaigns import parse_md
    verdicts, _ = parse_md(_repo_file(src["file"]))
    cands = verdicts
    for k, v in (src.get("match") or {}).items():
        cands = [x for x in cands if x["fields"].get(k) == str(v)]
    idx = src.get("index", 0)
    if not cands or idx >= len(cands):
        raise LookupError(f"no VERDICT line matching {src.get('match')} in {src['file']}")
    line = cands[idx]
    extra = _bracket_tag(line["raw"], f"bracket tag on {src['file']}:{line['line']}")
    if src.get("word"):
        return line["word"], line["word"], f"{src['file']}:{line['line']} (VERDICT word)", extra
    raw = line["fields"][src["field"]]
    return _num_or_text(raw), raw, f"{src['file']}:{line['line']} VERDICT {src['field']}", extra


def _r_best(src):
    from .campaigns import parse_md
    _, best = parse_md(_repo_file(src["file"]))
    cands = [b for b in best if b["kind"] == src.get("best_kind", "BEST_PASSING_FLUX")]
    for k, v in (src.get("match") or {}).items():
        cands = [x for x in cands if x["fields"].get(k) == str(v)]
    if not cands:
        raise LookupError(f"no BEST line matching {src.get('match')} in {src['file']}")
    line = cands[src.get("index", 0)]
    raw = line["fields"][src["field"]]
    extra = _bracket_tag(line["raw"], f"bracket tag on {src['file']}:{line['line']}")
    return _num_or_text(raw), raw, f"{src['file']}:{line['line']} {line['kind']} {src['field']}", extra


def _check_pattern(pattern) -> str:
    if not isinstance(pattern, str) or not pattern:
        raise ValueError("md_regex pattern must be a non-empty string")
    if len(pattern) > MAX_PATTERN:
        raise ValueError(f"md_regex pattern longer than {MAX_PATTERN} characters")
    if _NESTED_QUANT.search(pattern) or re.search(r"\\[1-9]", pattern):
        raise ValueError("md_regex pattern rejected: nested quantifier or backreference")
    re.compile(pattern)
    return pattern


def _r_md_regex(src):
    pattern = _check_pattern(src.get("pattern"))
    text = _repo_file(src["file"]).read_text(encoding="utf-8")
    flags = re.MULTILINE if "m" in src.get("flags", "m") else 0
    m = re.search(pattern, text, flags)
    if not m:
        raise LookupError(f"pattern not found in {src['file']}: {pattern}")
    raw = m.group(src.get("group", 1)).strip()
    line_no = text.count("\n", 0, m.start()) + 1
    extra = _bracket_tag(m.group(0), f"bracket tag inside the matched text ({src['file']}:{line_no})")
    return _num_or_text(raw), raw, f"{src['file']}:{line_no}", extra


def _r_json(src):
    doc = json.loads(_repo_file(src["file"]).read_text(encoding="utf-8"))
    obj = doc
    if src["key"] in obj:  # literal key first ("passed_pm0.03")
        obj = obj[src["key"]]
    else:
        for part in src["key"].split("."):
            obj = obj[part]
    extra = {}
    if isinstance(doc, dict) and isinstance(doc.get("tag_chain"), str):
        t = widest(tags_in_text(doc["tag_chain"]) | {doc["tag_chain"].strip("[] ")})
        if t:
            extra = {"tag": t, "tag_source": f"{src['file']} tag_chain"}
    return obj, json.dumps(obj), f"{src['file']} [{src['key']}]", extra


def _csv_rows(rel):
    import csv
    with _repo_file(rel).open(newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh))


def _cell_match(cell: str, want) -> bool:
    if isinstance(want, dict):
        try:
            x = float(cell)
        except (TypeError, ValueError):
            return False
        return (want.get("min") is None or x >= want["min"]) and \
               (want.get("max") is None or x <= want["max"])
    if isinstance(want, list):
        return any(_cell_match(cell, w) for w in want)
    if isinstance(want, bool):
        return cell == str(want)
    if isinstance(want, (int, float)):
        try:
            return float(cell) == float(want)
        except (TypeError, ValueError):
            return False
    return cell == str(want)


def _sort_key(s: str):
    v = _num_or_text(s)
    return (0, v, "") if isinstance(v, float) else (1, 0.0, str(s))


def _filter_rows(rows, filt):
    for k, v in (filt or {}).items():
        rows = [r for r in rows if _cell_match(r.get(k), v)]
    return rows


def _csv_tag(rows, col, rel) -> dict:
    tcol = f"{col}_tag"
    if not rows or tcol not in rows[0]:
        return {}
    tags = set()
    for r in rows:
        tags |= tags_in_text(r.get(tcol) or "") | ({(r.get(tcol) or "").strip("[] ")} & {"V", "DR", "E", "A"})
    t = widest(tags)
    return {"tag": t, "tag_source": f"{rel} column {tcol}"} if t else {}


_T_KEYS = ("T_K", "T_hs", "T_hs_K")


def _r_csv_cell(src):
    rows = _filter_rows(_csv_rows(src["file"]), src["match"])
    if len(rows) != 1:
        raise LookupError(f"{src['file']}: {len(rows)} rows match {src['match']} (need 1)")
    raw = rows[0][src["col"]]
    extra = _csv_tag(rows, src["col"], src["file"])
    for k in _T_KEYS:
        if k in (src.get("match") or {}):
            extra["at_T"] = float(src["match"][k])
    return _num_or_text(raw), raw, f"{src['file']} row {src['match']} col {src['col']}", extra


def _r_csv_agg(src):
    rows = _filter_rows(_csv_rows(src["file"]), src.get("filter"))
    agg = src["agg"]
    where = f"{src['file']} {agg}({src.get('col', '*')}) where {src.get('filter') or {}}"
    if agg == "count":
        return len(rows), str(len(rows)), where, {}
    extra = _csv_tag(rows, src["col"], src["file"])
    vals = [r[src["col"]] for r in rows]
    if agg in ("values", "unique"):
        seq = vals if agg == "values" else sorted(set(vals), key=_sort_key)
        return [_num_or_text(v) for v in seq], ",".join(seq), where, extra
    nums = [(float(v), v) for v in vals if v not in ("", "nan", "None") and math.isfinite(float(v))]
    if not nums:
        raise LookupError(f"no finite values for {where}")
    pick = min(nums) if agg == "min" else max(nums)
    return pick[0], pick[1], where, extra


def _r_commit(src):
    from .campaigns import _commit_and_date, _load_manifest
    if not re.match(r"^[A-Za-z0-9_][A-Za-z0-9_/.-]*$", src["dir"]) or ".." in src["dir"]:
        raise ValueError(f"invalid campaign dir {src['dir']!r}")
    commit, generated, csrc = _commit_and_date(src["dir"], _load_manifest(OUT / src["dir"]))
    return commit, commit, f"out/{src['dir']} ({csrc})", {}


def _apply_override(d, path: str, value) -> None:
    parts = path.split(".")
    if parts[0] == "quantum" and len(parts) == 1:
        d.quantum = dict(value)
        return
    obj = getattr(d, parts[0])
    if len(parts) == 2:
        setattr(obj, parts[1], value)
        return
    # block.dictfield.key (e.g. drive.diode.tau_pulse_ns): copy-on-write dict
    container = dict(getattr(obj, parts[1]) or {})
    container[parts[2]] = value
    setattr(obj, parts[1], container)


def _r_live(src):
    from fsim_core.device import evaluate

    from .api_cards import load_design
    from .jobs import headline_design
    d = load_design(src["card"])
    for path, value in (src.get("overrides") or {}).items():
        _apply_override(d, path, value)
    if src.get("mode", "point") == "headline":
        d = headline_design(d)
    T_grid = src.get("T_grid") or [float(d.thermal.T_hs)]
    scalars = evaluate(d, T_grid=T_grid)["scalars"]
    v = to_jsonable(scalars[src["scalar"]])
    where = (f"live fsim_core.device.evaluate(cards/{src['card']}.yaml, mode={src.get('mode', 'point')}, "
             f"overrides={src.get('overrides') or {}}, T_grid={T_grid})['scalars']['{src['scalar']}']")
    tag, tsrc = chain_of(scalars)
    extra = {"tag": tag, "tag_source": f"live run tag_chain ({tsrc})"}
    if len(T_grid) == 1:
        extra["at_T"] = float(T_grid[0])
    return v, repr(v), where, extra


RESOLVERS = {"verdict": _r_verdict, "best": _r_best, "md_regex": _r_md_regex, "json": _r_json,
             "csv_cell": _r_csv_cell, "csv_agg": _r_csv_agg, "commit": _r_commit, "live": _r_live}


def _bracket_tag(text, where) -> dict:
    t = widest(tags_in_text(text))
    return {"tag": t, "tag_source": where} if t else {}


# ----------------------------------------------------------- hashing

def _referenced_files(story: dict) -> list:
    files = set()
    for sc in story.get("scenes", []):
        for n in sc.get("numbers", []):
            src = n.get("source") or {}
            if src.get("file"):
                files.add(src["file"])
            if src.get("kind") == "live" and src.get("card"):
                files.add(f"cards/{src['card']}.yaml")
    return sorted(files)


def _file_digest(rel: str) -> bytes:
    p = (ROOT / rel).resolve()
    allowed = (OUT.resolve(), CARDS.resolve())
    if not any(a in p.parents for a in allowed) or not p.is_file():
        return b"missing"
    return hashlib.sha256(p.read_bytes()).digest()


def _numbers_digest(story: dict) -> bytes:
    rows = []
    for sc in story.get("scenes", []):
        for n in sc.get("numbers", []):
            rows.append([sc.get("id"), n.get("name"), n.get("label"), n.get("source"),
                         {k: n.get(k) for k in ("value", "raw", "resolved_from", "tag", "tag_source", "at_T")}])
    return hashlib.sha256(json.dumps(rows, sort_keys=True, allow_nan=True, default=str).encode()).digest()


def current_hash(story: dict) -> str:
    """fsim_core source hash + every out/**/manifest.json + every out/ file and
    live-source card a number reads + every number's resolved fields. Any
    change marks a baked story stale."""
    h = hashlib.sha256()
    h.update(fsim_core_hash().encode())
    out = ROOT / "out"
    for p in sorted(out.glob("**/manifest.json")):
        if any(part.startswith("presentation") for part in p.relative_to(out).parts):
            continue
        h.update(p.relative_to(ROOT).as_posix().encode())
        h.update(hashlib.sha256(p.read_bytes()).digest())
    for rel in _referenced_files(story):
        h.update(rel.encode())
        h.update(_file_digest(rel))
    h.update(b"numbers")
    h.update(_numbers_digest(story))
    return h.hexdigest()


# ----------------------------------------------------------- prose

def _compact(v):
    if v is None:
        return "n/a"
    if isinstance(v, bool):
        return "true" if v else "false"
    if isinstance(v, int):
        return str(v)
    if isinstance(v, float):
        if not math.isfinite(v):
            return "n/a"
        if v.is_integer() and abs(v) < 1e15:
            return str(int(v))
        return f"{v:.0f}" if abs(v) >= 1000 else f"{v:.4g}"
    if isinstance(v, (list, tuple)):
        return ", ".join(_compact(x) for x in v)
    return str(v)


def _render_number(n: dict, index, attr, spec):
    if attr == "tag":
        return n.get("tag") or n.get("tag_declared") or "?"
    if attr == "T":
        return _compact(n.get("at_T"))
    if attr == "say":
        word = str(n.get("raw") if isinstance(n.get("raw"), str) and n.get("raw") else _compact(n.get("value")))
        return STATUS_PHRASES.get(word, word.replace("_", " "))
    v = n.get("value")
    if index is not None:
        if not isinstance(v, list):
            raise LookupError(f"{n.get('name')} is not a list")
        v = v[int(index)]
    if spec:
        return format(v, spec)
    kind = (n.get("source") or {}).get("kind")
    if index is None and kind in RAW_TEXT_KINDS and isinstance(n.get("raw"), str) and n["raw"]:
        return n["raw"]
    return _compact(v)


def _number_index(story: dict) -> tuple:
    idx, dup = {}, []
    for sc in story.get("scenes", []):
        for n in sc.get("numbers", []):
            nm = n.get("name")
            if not nm:
                continue
            if nm in idx:
                dup.append(nm)
            idx[nm] = n
    return idx, dup


def interpolate(text, numbers: dict, errors: list, where: str):
    if not isinstance(text, str):
        return text

    def sub(m):
        name, index, attr, spec = m.group(1), m.group(2), m.group(3), m.group(4)
        n = numbers.get(name)
        if n is None:
            errors.append({"where": where, "placeholder": m.group(0), "error": "unknown number name"})
            return m.group(0)
        if n.get("value") is None and attr is None:
            errors.append({"where": where, "placeholder": m.group(0), "error": "number not baked"})
            return "n/a"
        try:
            return _render_number(n, index, attr, spec)
        except (LookupError, ValueError, TypeError, IndexError) as exc:
            errors.append({"where": where, "placeholder": m.group(0), "error": str(exc)})
            return m.group(0)

    return PLACEHOLDER.sub(sub, text)


def render(story: dict) -> dict:
    """Interpolate every scene's prose from the baked numbers; keep the
    templates as <field>_template."""
    numbers, dup = _number_index(story)
    errors = [{"where": "numbers", "placeholder": d, "error": "duplicate number name"} for d in sorted(set(dup))]
    for sc in story.get("scenes", []):
        for f in PROSE_FIELDS:
            if f in sc:
                sc[f"{f}_template"] = sc[f]
                sc[f] = interpolate(sc[f], numbers, errors, f"{sc.get('id')}.{f}")
    story["prose_errors"] = errors
    return story


# ----------------------------------------------------------- API

def list_stories() -> list:
    out = []
    for p in sorted(STORIES.glob("*.json")):
        try:
            s = json.loads(p.read_text(encoding="utf-8"))
        except (OSError, ValueError):
            continue
        out.append({"name": p.stem, "title": s.get("title"), "n_scenes": len(s.get("scenes", [])),
                    "baked_at": s.get("baked_at")})
    return out


def load_story(name: str) -> dict:
    p = _path(name)
    if not p.is_file():
        raise StoryError(404, f"no story {name!r}")
    return json.loads(p.read_text(encoding="utf-8"))


def get_story(name: str) -> dict:
    story = load_story(name)
    cur = current_hash(story)
    story["current_hash"] = cur
    story["baked"] = bool(story.get("baked_hash"))
    story["stale"] = story.get("baked_hash") != cur
    return render(story)


def _write(name: str, story: dict) -> None:
    p = _path(name)
    STORIES.mkdir(parents=True, exist_ok=True)
    p.write_text(json.dumps(story, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
                 encoding="utf-8")


def _number_key(sc_id, n):
    return (sc_id, n.get("name") or n.get("label"), json.dumps(n.get("source"), sort_keys=True))


def save_story(name: str, story: dict) -> dict:
    """PUT: authoring fields only. Resolved number fields and bake state are
    stripped from the body (H4) and restored from the stored story for
    numbers whose (scene, name, source) is unchanged; rendered prose is
    replaced by its <field>_template when the body carries one."""
    if not isinstance(story.get("scenes"), list):
        raise StoryError(400, "story must carry scenes: [...]")
    story = copy.deepcopy(story)
    try:
        stored = load_story(name)
    except StoryError:
        stored = {}
    prev = {}
    for sc in stored.get("scenes", []):
        for n in sc.get("numbers", []):
            prev[_number_key(sc.get("id"), n)] = n
    for k in STORY_STATE_FIELDS:
        story.pop(k, None)
    for sc in story["scenes"]:
        if not isinstance(sc, dict):
            raise StoryError(400, "each scene must be an object")
        for f in PROSE_FIELDS:
            tpl = sc.pop(f"{f}_template", None)
            if tpl is not None:
                sc[f] = tpl
        nums = []
        for n in sc.get("numbers", []) or []:
            if not isinstance(n, dict):
                raise StoryError(400, "each number must be an object")
            n = {k: v for k, v in n.items() if k not in RESOLVED_FIELDS}
            old = prev.get(_number_key(sc.get("id"), n))
            if old is not None:
                for k in RESOLVED_FIELDS:
                    if k in old:
                        n[k] = old[k]
            nums.append(n)
        sc["numbers"] = nums
    for k in ("baked_hash", "baked_fsim_core_hash", "baked_at", "bake_errors"):
        if k in stored:
            story[k] = stored[k]
    _write(name, story)
    return get_story(name)


def bake(name: str, run_live: bool = True) -> dict:
    """Resolve every number's source and store value / raw / resolved_from /
    tag, plus baked_hash. Live sources run fsim_core inline (seconds); with
    run_live=False they keep their previous value and are marked."""
    import datetime as _dt
    story = load_story(name)
    for k in ("current_hash", "stale", "baked", "prose_errors"):
        story.pop(k, None)
    errors = []
    for sc in story.get("scenes", []):
        for n in sc.get("numbers", []):
            if "tag_declared" not in n and "tag" in n:
                n["tag_declared"] = n["tag"]  # the author's tag, kept for comparison
            src = n.get("source") or {}
            kind = src.get("kind")
            if kind == "live" and not run_live:
                n["live_skipped"] = True
                continue
            n.pop("live_skipped", None)
            try:
                if kind not in RESOLVERS:
                    raise ValueError(f"unknown source kind {kind!r}")
                value, raw, where, extra = RESOLVERS[kind](copy.deepcopy(src))
                n["value"] = to_jsonable(value)
                n["raw"] = raw
                n["resolved_from"] = where
                if extra.get("tag"):
                    n["tag"], n["tag_source"] = extra["tag"], extra["tag_source"]
                else:
                    n["tag"] = n.get("tag_declared") or "A"
                    n["tag_source"] = ("declared in the story (the source prints no tag)"
                                       if n.get("tag_declared") else "no tag in source; assumed widest")
                if "at_T" in extra:
                    n["at_T"] = extra["at_T"]
                else:
                    n.pop("at_T", None)
                n.pop("error", None)
            except Exception as exc:  # noqa: BLE001 -- recorded per number, never typed in
                n["value"] = None
                n["raw"] = None
                n["error"] = f"{type(exc).__name__}: {exc}"
                errors.append({"scene": sc.get("id"), "label": n.get("label"), "error": n["error"]})
    numbers, dup = _number_index(story)
    for d in sorted(set(dup)):
        errors.append({"scene": None, "label": d, "error": "duplicate number name"})
    prose_errs = []
    for sc in story.get("scenes", []):
        for f in PROSE_FIELDS:
            interpolate(sc.get(f), numbers, prose_errs, f"{sc.get('id')}.{f}")
    errors += [{"scene": e["where"], "label": e["placeholder"], "error": e["error"]} for e in prose_errs]
    story["baked_hash"] = current_hash(story)
    story["baked_fsim_core_hash"] = fsim_core_hash()
    story["baked_at"] = _dt.datetime.now(_dt.timezone.utc).isoformat(timespec="seconds")
    story["bake_errors"] = errors
    _write(name, story)
    return get_story(name)
