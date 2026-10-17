"""Display-symbol checks of the composable registries (devlib-r2-math).

Every ParamSpec of the material, drive and setup axes, and every observable, carries a `tex` symbol; every `$...$`
segment of a help / label / summary string and every `tex` parses in KaTeX (the vendored 0.16.11, run under node with
throwOnError and trust off), so the Studio never falls back to raw text. `tex` is UI-only: it is absent from
`ParamSpec.to_meta()` when unset, so existing metadata dicts are unchanged. Structure params carry symbols too (devlib-r2-fix1b);
the repo gamma* wording (GAMMA_STAR_CONVENTION, the dephasing material help) is the pure-dephasing linewidth Gamma_ph.

Run from the repo root: python verify/verify_param_tex.py   (needs node on PATH)
Prints one line per check and "N/N param_tex checks passed"; exits 0 iff all pass.
"""
from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
import warnings

warnings.filterwarnings("ignore")
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

from fsim_core.compose import registry  # noqa: E402
from fsim_core.compose.params import ParamSpec  # noqa: E402

KATEX = os.path.join(ROOT, "fsim_studio", "web", "vendor", "katex", "katex.min.js")
IN_SCOPE_AXES = ("material", "drive", "setup", "structure")
DOLLAR = re.compile(r"(?<!\\)\$\$?(.+?)(?<!\\)\$\$?", re.S)

RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append(bool(ok))
    print(f"{'PASS' if ok else 'FAIL'}  {name}{('  (' + detail + ')') if detail else ''}")


def segments(text):
    """The TeX bodies of the `$...$` segments of a string (an escaped \\$ is a literal dollar)."""
    t = text.replace("\\$", "")
    return [m.group(1) for m in DOLLAR.finditer(t)]


def katex_parse(sources):
    """-> list of error strings ('' = parsed) from the vendored KaTeX under node."""
    node = shutil.which("node")
    if not node:
        raise SystemExit("node not found on PATH")
    js = (
        "const k=require(process.argv[1]);const src=JSON.parse(require('fs').readFileSync(0,'utf8'));"
        "const out=src.map(s=>{try{k.renderToString(s,{throwOnError:true,trust:false,strict:'ignore'});return ''}"
        "catch(e){return String(e.message)}});process.stdout.write(JSON.stringify(out));"
    )
    r = subprocess.run([node, "-e", js, KATEX], input=json.dumps(sources), capture_output=True, text=True, timeout=120)
    if r.returncode != 0:
        raise SystemExit(f"node failed: {r.stderr[:400]}")
    return json.loads(r.stdout)


def main():
    registry.ensure_loaded()
    R = registry.REGISTRY
    params = []   # (axis, kind, ParamSpec)
    for axis in registry.AXES:
        if axis in ("solver", "observable"):
            continue
        for k in R.kinds(axis):
            for p in R.get(axis, k).params:
                params.append((axis, k, p))
    obs = [R.get("observable", k) for k in R.kinds("observable")]

    check("katex.min.js is vendored (offline)", os.path.isfile(KATEX), KATEX)
    check("vendored KaTeX fonts are present (woff2)", len([f for f in os.listdir(os.path.join(os.path.dirname(KATEX), "fonts")) if f.endswith(".woff2")]) >= 19)

    # ParamSpec.tex is optional metadata
    plain = ParamSpec("x", "1", 0.0)
    check("ParamSpec.tex defaults to None and is absent from to_meta()", plain.tex is None and "tex" not in plain.to_meta()
          and set(plain.to_meta()) == {"unit", "tag", "source"})
    withtex = ParamSpec("x", "1", 0.0, tex=r"\eta")
    check("ParamSpec.tex appears in to_meta() when set, other keys unchanged",
          withtex.to_meta().get("tex") == r"\eta" and {k: v for k, v in withtex.to_meta().items() if k != "tex"} == plain.to_meta())

    inscope = [(a, k, p) for a, k, p in params if a in IN_SCOPE_AXES]
    missing = [f"{a}:{k}.{p.name}" for a, k, p in inscope if not p.tex]
    check(f"every material / drive / setup param has a tex symbol ({len(inscope)} params)", not missing, ", ".join(missing[:6]))
    obs_missing = [o.id for o in obs if not o.tex]
    check(f"every observable has a tex symbol ({len(obs)} observables)", not obs_missing, ", ".join(obs_missing[:6]))
    struct = [(a, k, p) for a, k, p in params if a == "structure"]
    struct_missing = [f"{k}.{p.name}" for _a, k, p in struct if not p.tex]
    check(f"every structure param has a tex symbol ({len(struct)} params)", not struct_missing, ", ".join(struct_missing[:6]))
    # help texts of the structures are marked up: every structure param help has at least one $...$ symbol
    plain_help = [f"{k}.{p.name}" for _a, k, p in struct if not segments(p.help or "")]
    check("every structure help text carries $...$ markup", not plain_help, ", ".join(plain_help[:6]))

    # consistency: one id never carries two different symbols within an axis
    by_name = {}
    for a, k, p in inscope:
        if p.tex:
            by_name.setdefault((a, p.name), set()).add(p.tex)
    clash = [f"{a}.{n}: {sorted(v)}" for (a, n), v in by_name.items() if len(v) > 1]
    check("an id carries the same symbol on every kind of its axis", not clash, "; ".join(clash[:3]))

    # no ASCII snake_case leaked into a symbol
    snake = [f"{a}:{k}.{p.name}={p.tex}" for a, k, p in inscope if p.tex and re.search(r"[A-Za-z]{2,}_[A-Za-z]{2,}", p.tex.replace(r"\_", ""))]
    check("no tex symbol contains a raw snake_case word", not snake, "; ".join(snake[:3]))

    # every tex and every $...$ segment of the help / label / summary strings parses
    items = []   # (where, tex)
    for a, k, p in params:
        if p.tex:
            items.append((f"{a}:{k}.{p.name}.tex", p.tex))
        for i, s in enumerate(segments(p.help or "")):
            items.append((f"{a}:{k}.{p.name}.help[{i}]", s))
    for o in obs:
        if o.tex:
            items.append((f"observable:{o.id}.tex", o.tex))
        for i, s in enumerate(segments(o.label or "")):
            items.append((f"observable:{o.id}.label[{i}]", s))
    for axis in ("material", "drive", "setup", "structure"):
        for k in R.kinds(axis):
            e = R.get(axis, k)
            for fld in ("summary", "label"):
                for i, s in enumerate(segments(getattr(e, fld, "") or "")):
                    items.append((f"{axis}:{k}.{fld}[{i}]", s))
    errs = katex_parse([t for _w, t in items])
    bad = [f"{w}: {e[:80]}" for (w, _t), e in zip(items, errs) if e]
    check(f"all {len(items)} TeX strings parse in KaTeX (throwOnError, trust off)", not bad, " | ".join(bad[:3]))
    print(f"INFO  {len([1 for w, _ in items if '.help[' in w])} help segments, {len([1 for w, _ in items if w.endswith('.tex')])} symbols checked")

    from fsim_core.compose import rates as _rates
    gs_items = segments(_rates.GAMMA_STAR_CONVENTION)
    gs_err = [e for e in katex_parse(gs_items) if e]
    check("GAMMA_STAR_CONVENTION is $...$ markup naming Gamma_ph and parses", bool(gs_items) and not gs_err
          and r"\Gamma_{\mathrm{ph}}" in _rates.GAMMA_STAR_CONVENTION, "; ".join(gs_err[:2]))
    repo_gs = [f"{a}:{k}.{p.name}" for a, k, p in params if "repo gamma" in (p.help or "").lower() or "repo convention" in (p.help or "").lower()]
    check("no material help still says 'repo gamma*' / 'repo convention'", not repo_gs, ", ".join(repo_gs[:4]))

    # help strings keep balanced dollars (an odd count would swallow the prose after it)
    unbalanced = []
    for a, k, p in params:
        h = (p.help or "").replace("\\$", "")
        if h.count("$") % 2:
            unbalanced.append(f"{a}:{k}.{p.name}")
    for o in obs:
        if (o.label or "").replace("\\$", "").count("$") % 2:
            unbalanced.append(f"observable:{o.id}")
    check("every help / observable label has an even number of $", not unbalanced, ", ".join(unbalanced[:4]))

    # the registry payload (what the Studio reads) carries tex on params and observables
    pay = registry.registry_payload()
    pp = [q for ent in pay["drives"] + pay["materials"] + pay["setups"] for q in ent["params"]]
    check("registry payload: every in-scope param dict has its tex", all(q.get("tex") for q in pp), f"{len(pp)} params")
    check("registry payload: observables carry tex", all(o.get("tex") for o in pay["observables"]), f"{len(pay['observables'])}")

    # the help of a few well-known params uses the expected markup (spot check against the literal)
    g = R.get("drive", "resonant_pulse").param("fwhm_ps")
    check("spot check: resonant_pulse fwhm_ps symbol and help", g.tex == r"\tau_{\mathrm{FWHM}}" and r"$\Omega^2$" in g.help)
    m = R.get("material", sorted(R.kinds("material"))[0]).param("gamma_star_4K_ns")
    check("spot check: gamma_star_4K_ns symbol is Gamma_ph(4 K)", m.tex == r"\Gamma_{\mathrm{ph}}(4\,\mathrm K)" and r"\Gamma_{\mathrm{ph}}" in m.help)
    o = R.get("observable", "g2_0")
    check("spot check: g2_0 observable symbol and label", o.tex == "g^{(2)}(0)" and o.label == "$g^{(2)}(0)$")

    n = len(RESULTS)
    ok = sum(RESULTS)
    print(f"{ok}/{n} param_tex checks passed")
    return 0 if ok == n else 1


if __name__ == "__main__":
    sys.exit(main())
