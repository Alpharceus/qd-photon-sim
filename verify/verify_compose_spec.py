"""Checks for the composable-engine skeleton (fsim_core/compose, P1).

Covers params (ParamSpec / leaf rules / central / corners / lhs / widest_tag), the registry,
ComposedSpec round trip and validation, the additive `platform: composed` arm in
fsim_core/device.py, legacy save() byte-identity against the HEAD version of device.py for every
shipped design card, and the fsim_core_hash glob.

Run: python verify/verify_compose_spec.py   (exit code 0 iff all pass)
"""
from __future__ import annotations

import copy
import importlib.util
import json
import math
import subprocess
import sys
import tempfile
from dataclasses import dataclass
from pathlib import Path

import yaml

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from fsim_core import design_meta                                      # noqa: E402
from fsim_core.compose import params as P                              # noqa: E402
from fsim_core.compose import registry as R                            # noqa: E402
from fsim_core.compose import spec as S                                # noqa: E402
from fsim_core.compose.types import EmitterModel, StructureResponse, DrivePlan  # noqa: E402
from fsim_core.device import DeviceDesign, evaluate                    # noqa: E402

RESULTS = []


def check(name, ok, detail=""):
    RESULTS.append(bool(ok))
    print(f"[{'PASS' if ok else 'FAIL'}] {name}" + (f"  ({detail})" if detail and not ok else ""))


def raises(exc, fn, *a, contains=(), **kw):
    try:
        fn(*a, **kw)
    except exc as e:
        msg = str(e)
        return all(c in msg for c in contains), msg
    except Exception as e:  # wrong exception type
        return False, f"{type(e).__name__}: {e}"
    return False, "no exception"


# ---------------------------------------------------------------- params
num = P.ParamSpec("Q", "1", 5000.0, lo=10.0, hi=1e6, scale="log", tag="E", source="class")
txt = P.ParamSpec("shape", "", "gauss", choices=("gauss", "rect"), tag="A", source="conv")

m = num.to_meta()
numeric_meta = design_meta.META["dot.delta_xx"]            # unit, tag, lo, hi, source
choice_meta = design_meta.META["dot.lineshape"]            # unit, tag, choices, source
check("to_meta keys of a numeric ParamSpec == numeric META entry keys",
      set(m) == set(numeric_meta), f"{sorted(m)} vs {sorted(numeric_meta)}")
check("to_meta keys of a choice ParamSpec == choice META entry keys",
      set(txt.to_meta()) == {"unit", "tag", "source", "choices"} == set(choice_meta) | set(),
      f"{sorted(txt.to_meta())} vs {sorted(choice_meta)}")
check("every META entry key set is a subset of {unit,tag,lo,hi,source,choices}",
      all(set(v) <= {"unit", "tag", "lo", "hi", "source", "choices"} for v in design_meta.META.values()))
check("to_meta values carried", m["unit"] == "1" and m["tag"] == "E" and m["lo"] == 10.0 and m["hi"] == 1e6)

lf = P.resolve_leaf(6124, num)
check("bare scalar -> value leaf with default tag/source",
      lf.kind == "value" and lf.value == 6124 and lf.tag == "E" and lf.source == "class")
lf = P.resolve_leaf({"value": 2.0, "unit": "1", "tag": "V", "where": "p.2"}, num)
check("mapping value leaf", lf.kind == "value" and lf.tag == "V" and lf.where == "p.2")
lf = P.resolve_leaf({"range": [1e3, 1e5], "unit": "1", "tag": "E", "source": "x", "scale": "log"}, num)
check("mapping range leaf", lf.kind == "range" and lf.range == (1e3, 1e5) and lf.scale == "log")
lf = P.resolve_leaf({"choice": ["gauss", "rect"], "unit": ""}, txt)
check("mapping choice leaf", lf.kind == "choice" and lf.choice == ("gauss", "rect"))
ok, msg = raises(ValueError, P.resolve_leaf, {"value": 1, "unit": "ns"}, num, "structure.Q",
                 contains=("structure.Q", "'ns'", "'1'"))
check("bad unit -> ValueError naming path, given and expected unit", ok, msg)
check("two of value/range -> ValueError",
      raises(ValueError, P.resolve_leaf, {"value": 1, "range": [1, 2], "unit": "1"}, num)[0])
check("none of value/range/choice -> ValueError",
      raises(ValueError, P.resolve_leaf, {"unit": "1"}, num)[0])
check("bad tag -> ValueError",
      raises(ValueError, P.resolve_leaf, {"value": 1, "unit": "1", "tag": "Z"}, num)[0])
check("reversed range -> ValueError",
      raises(ValueError, P.resolve_leaf, {"range": [5, 1], "unit": "1"}, num)[0])
check("choice outside spec.choices -> ValueError",
      raises(ValueError, P.resolve_leaf, "sech", txt)[0])
check("non-numeric value for numeric param -> ValueError",
      raises(ValueError, P.resolve_leaf, "abc", num)[0])

lin = P.resolve_leaf({"range": [2.0, 8.0], "unit": "1", "scale": "lin"}, num)
log = P.resolve_leaf({"range": [2.0, 8.0], "unit": "1", "scale": "log"}, num)
check("central: lin midpoint", P.central(lin) == 5.0)
check("central: log geometric mean", abs(P.central(log) - 4.0) < 1e-12)

leaves = {f"a{i}": P.resolve_leaf({"range": [1.0 + i, 2.0 + i], "unit": "1", "scale": "lin"}, num)
          for i in range(3)}
leaves["fixed"] = P.resolve_leaf(7.0, num)
c = P.corners(leaves)
check("corners: 2^3 = 8 points <= 16, extremes present, fixed leaf constant",
      len(c) == 8 and {"a0": 1.0, "a1": 2.0, "a2": 3.0, "fixed": 7.0} in c
      and {"a0": 2.0, "a1": 3.0, "a2": 4.0, "fixed": 7.0} in c and all(p["fixed"] == 7.0 for p in c))
many = {f"b{i}": P.resolve_leaf({"range": [1.0, 2.0], "unit": "1"}, num) for i in range(8)}
c = P.corners(many)
check("corners: 8 ranges (256) truncated to <= 16 and keep both extremes",
      len(c) <= 16 and {k: 1.0 for k in many} in c and {k: 2.0 for k in many} in c)
s = P.lhs(leaves, n=32, seed=1)
inside = all(1.0 <= p["a0"] <= 2.0 and 2.0 <= p["a1"] <= 3.0 and p["fixed"] == 7.0 for p in s)
strata = sorted(int((p["a0"] - 1.0) * 32) for p in s) == list(range(32))
check("lhs: 32 samples, inside ranges, one per stratum, reproducible",
      len(s) == 32 and inside and strata and s == P.lhs(leaves, n=32, seed=1))
check("widest_tag order V < DR < E < A",
      P.widest_tag("V") == "V" and P.widest_tag("V", "DR") == "DR" and P.widest_tag("DR", "E", "V") == "E"
      and P.widest_tag("A", "V") == "A" and raises(ValueError, P.widest_tag, "Q")[0])

# ---------------------------------------------------------------- registry
@dataclass(frozen=True)
class _Mat(R.Material):
    def build(self, vals, T_K):
        return EmitterModel(E_X_meV=1380.0, lambda_nm=897.0, n_host=3.5, gamma0_ns=1.5, k_nr_ns=0.0,
                            gamma_star_ns=0.0, phonon=None, delta_xx_meV=2.0, gamma_xx_ratio=1.0,
                            fss_ueV=0.0, transitions=frozenset({"X"}), dipole="linear",
                            dEdT_meV_K=0.0, sigma_sd_ueV=0.0, tau_sd_ns=1.0)


@dataclass(frozen=True)
class _Str(R.Structure):
    def response(self, vals, em, T_K, NA):
        return StructureResponse(F_P=1.0, kappa_meV=None, E_mode_meV=None, g_meV=None)


@dataclass(frozen=True)
class _Drv(R.Drive):
    def plan(self, vals, em, rt):
        return DrivePlan(family="resonant", timing="pulsed", f_rep_MHz=80.0)


def mk(cls, kind, **kw):
    return cls(kind=kind, label=kind, summary="throwaway", **kw)


t_mat = mk(_Mat, "t_mat", params=(P.ParamSpec("lambda_nm", "nm", 900.0, lo=300.0, hi=1600.0, tag="E", source="t"),),
           provides=frozenset({"X"}))
t_str = mk(_Str, "t_str", params=(P.ParamSpec("Q", "1", 5000.0, lo=1.0, hi=1e6, scale="log", tag="E", source="t"),))
t_drv = mk(_Drv, "t_drv", params=(num, txt), requires=frozenset({"X"}))
t_drv_xx = mk(_Drv, "t_drv_xx", requires=frozenset({"XX"}))

fresh = R.Registry()
check("a fresh Registry is empty", all(fresh.kinds(a) == [] for a in R.AXES))
fresh.register(t_mat)
check("register + get + kinds", fresh.get("material", "t_mat") is t_mat and fresh.kinds("material") == ["t_mat"])
ok, msg = raises(ValueError, fresh.register, t_mat, contains=("duplicate", "t_mat"))
check("register duplicate kind -> ValueError", ok, msg)
ok, msg = raises(KeyError, R.get, "material", "nope", fresh, contains=("t_mat",))
check("get unknown kind -> KeyError listing registered kinds", ok, msg)
check("get unknown axis -> KeyError", raises(KeyError, fresh.get, "nope", "x")[0])
check("Entry rejects duplicate param names",
      raises(ValueError, mk, _Mat, "dup", params=(num, num))[0])
check("abstract Material cannot be instantiated directly", raises(TypeError, mk, R.Material, "abs")[0])
check("built-in setup lens_fibre_det registered on the default registry",
      "lens_fibre_det" in R.kinds("setup"))
chain = R.get("setup", "lens_fibre_det").chain({"NA": 0.7})
check("Setup.chain -> SetupChain with defaults filled",
      chain.NA == 0.7 and chain.eta_det == 0.5 and chain.zpl_filter is False)

fresh.register(t_str)
fresh.register(t_drv)
pay = R.registry_payload(fresh)
txt_payload = json.dumps(pay)
check("registry_payload JSON-serialisable with the design 2.1 top-level keys",
      set(pay) == {"schema", "materials", "structures", "drives", "setups", "solvers", "observables",
                   "tiers", "compat"} and pay["schema"] == "fsim-compose-registry/1" and len(txt_payload) > 10)
pm = pay["materials"][0]["params"][0]
check("payload params carry to_meta fields, help, group, visible_if, choices, default",
      {"unit", "tag", "lo", "hi", "source", "choices", "help", "group", "visible_if", "default"} <= set(pm)
      and pm["default"] == 900.0 and pm["unit"] == "nm")
check("default payload serialises too (lens_fibre_det present)",
      any(e["kind"] == "lens_fibre_det" for e in json.loads(json.dumps(R.registry_payload()))["setups"]))


def ensure_loaded_probe():
    """Missing lazy modules are skipped; an ImportError INSIDE an existing module propagates.
    The broken module lives in a temp dir appended to the package __path__ (nothing is written
    into the repo)."""
    import fsim_core.compose as pkg
    R.ensure_loaded(force=True)
    absent = [m for m in R._LAZY_MODULES if not (ROOT / "fsim_core" / "compose" / f"{m}.py").exists()]
    if not absent:
        return None            # every lazy module has shipped; nothing left to probe with
    name = absent[-1]
    with tempfile.TemporaryDirectory() as tdir:
        (Path(tdir) / f"{name}.py").write_text("import definitely_not_a_module_xyz\n", encoding="utf-8")
        pkg.__path__.append(tdir)
        sys.modules.pop(f"fsim_core.compose.{name}", None)
        try:
            R.ensure_loaded(force=True)
            return False
        except ModuleNotFoundError as e:
            return e.name == "definitely_not_a_module_xyz"
        finally:
            pkg.__path__.remove(tdir)
            sys.modules.pop(f"fsim_core.compose.{name}", None)
            R._STATE["loaded"] = False


r = ensure_loaded_probe()
check("ensure_loaded skips absent modules and propagates other import errors", r in (None, True))
R.ensure_loaded(force=True)

# ---------------------------------------------------------------- spec
GOOD = {
    "schema": "fsim-composed/1",
    "T_K": {"value": 4.2, "unit": "K", "tag": "V", "source": "paper"},
    "material": {"kind": "t_mat", "params": {"lambda_nm": {"value": 897.4, "unit": "nm", "tag": "V", "where": "p.2"}}},
    "structure": {"kind": "t_str", "params": {"Q": 6124}},
    "drive": {"kind": "t_drv", "params": {"shape": "gauss", "Q": {"range": [1e3, 1e4], "unit": "1", "scale": "log"}}},
    "setup": {"kind": "lens_fibre_det", "params": {"NA": 0.68, "zpl_filter": False}},
    "solver": "auto",
}
cs = S.ComposedSpec.from_dict(GOOD)
check("ComposedSpec round trip identical", cs.to_dict() == GOOD and S.ComposedSpec.from_dict(cs.to_dict()).to_dict() == GOOD)
bad = dict(GOOD, bogus=1)
ok, msg = raises(ValueError, S.ComposedSpec.from_dict, bad, contains=("bogus",))
check("unknown top-level key rejected", ok, msg)
check("wrong schema string rejected", raises(ValueError, S.ComposedSpec.from_dict, dict(GOOD, schema="x/2"))[0])
check("bad solver rejected", raises(ValueError, S.ComposedSpec.from_dict, dict(GOOD, solver="magic"))[0])
check("missing required key rejected",
      raises(ValueError, S.ComposedSpec.from_dict, {k: v for k, v in GOOD.items() if k != "drive"})[0])
check("unknown key inside a block rejected",
      raises(ValueError, S.ComposedSpec.from_dict, dict(GOOD, drive={"kind": "t_drv", "oops": 1}))[0])

# validate_spec against the DEFAULT registry (throwaway entries cleaned up in the finally below)
R.register(t_mat)
R.register(t_str)
R.register(t_drv)
R.register(t_drv_xx)
try:
    v = S.validate_spec(cs)
    check("validate_spec ok on a good spec; keys ok/solver/regime/issues/resolved",
          v["ok"] and set(v) == {"ok", "solver", "regime", "issues", "resolved"}
          and not [i for i in v["issues"] if i["level"] == "error"], v["issues"])
    rs = v["resolved"]
    check("resolved carries value|range, unit, tag per path",
          rs["material.lambda_nm"] == {"value": 897.4, "unit": "nm", "tag": "V"}
          and rs["drive.Q"]["range"] == [1e3, 1e4] and rs["drive.Q"]["tag"] == "E"
          and rs["structure.Q"]["value"] == 6124 and rs["T_K"]["tag"] == "V"
          and rs["setup.NA"]["value"] == 0.68 and rs["setup.eta_det"]["default"] is True)
    json.dumps(v)

    d2 = copy.deepcopy(GOOD)
    d2["drive"]["kind"] = "nope"
    v2 = S.validate_spec(S.ComposedSpec.from_dict(d2))
    iss = [i for i in v2["issues"] if i["path"] == "drive.kind"]
    check("unknown kind -> issue listing registered kinds", not v2["ok"] and iss and "t_drv" in iss[0]["message"]
          and "t_drv_xx" in iss[0]["message"], v2["issues"])

    d3 = copy.deepcopy(GOOD)
    d3["structure"]["params"]["not_a_param"] = 1
    v3 = S.validate_spec(S.ComposedSpec.from_dict(d3))
    check("unknown param -> issue", not v3["ok"] and any(i["path"] == "structure.params.not_a_param"
                                                          for i in v3["issues"]), v3["issues"])

    d4 = copy.deepcopy(GOOD)
    d4["structure"]["params"]["Q"] = {"value": 1, "unit": "ns"}
    v4 = S.validate_spec(S.ComposedSpec.from_dict(d4))
    check("bad unit on a leaf -> issue naming the path",
          not v4["ok"] and any("structure.Q" in i["message"] for i in v4["issues"]), v4["issues"])

    d5 = copy.deepcopy(GOOD)
    d5["drive"]["kind"] = "t_drv_xx"
    d5["drive"]["params"] = {}
    v5 = S.validate_spec(S.ComposedSpec.from_dict(d5))
    check("drive requires XX, material provides only X -> compat issue",
          not v5["ok"] and any("requires ['XX']" in i["message"] for i in v5["issues"]), v5["issues"])

    # axis with nothing registered -> "no <axis> registered yet"
    v6 = S.validate_spec(cs, registry=R.Registry())
    check("empty axis -> 'no <axis> registered yet'",
          not v6["ok"] and any(i["message"] == "no material registered yet" for i in v6["issues"]), v6["issues"])

    # evaluate() on a composed design: validation passes -> NotImplementedError("compose engine: P5")
    dd = DeviceDesign(name="t", platform="composed", composed=copy.deepcopy(GOOD))
    ok, msg = raises(NotImplementedError, evaluate, dd, contains=("compose engine: P5",))
    check("device.evaluate(composed) -> NotImplementedError('compose engine: P5')", ok, msg)
    dd2 = DeviceDesign(name="t", platform="composed", composed=d2)
    ok, msg = raises(ValueError, evaluate, dd2, contains=("drive.kind",))
    check("device.evaluate(composed, invalid) -> ValueError with the issues", ok, msg)
finally:
    for ax, k in (("material", "t_mat"), ("structure", "t_str"), ("drive", "t_drv"), ("drive", "t_drv_xx")):
        R.unregister(ax, k)

check("throwaway entries cleaned up", R.kinds("material") == [] or "t_mat" not in R.kinds("material"))

# ---------------------------------------------------------------- device load / save
with tempfile.TemporaryDirectory() as td:
    td = Path(td)
    doc = {"meta": {"name": "c1", "role": "device-design"},
           "design": {"name": "c1", "platform": "composed", "composed": copy.deepcopy(GOOD), "provenance": {}}}
    p1 = td / "c1.yaml"
    p1.write_text(yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")
    d = DeviceDesign.load(p1)
    check("device.load of a composed YAML works", d.platform == "composed" and d.composed == GOOD)
    p2, p3 = td / "c2.yaml", td / "c3.yaml"
    d.save(p2)
    DeviceDesign.load(p2).save(p3)
    t2, t3 = p2.read_text(encoding="utf-8"), p3.read_text(encoding="utf-8")
    saved = yaml.safe_load(t2)
    check("composed save/load byte-stable", t2 == t3)
    check("composed save writes only name/platform/composed/provenance",
          set(saved["design"]) == {"name", "platform", "composed", "provenance"}
          and saved["meta"]["role"] == "device-design")

    for blk in ("dot", "ret", "drive", "cavity", "filter", "aperture", "emission", "nitride"):
        bd = copy.deepcopy(doc)
        bd["design"][blk] = {}
        pb = td / f"bad_{blk}.yaml"
        pb.write_text(yaml.safe_dump(bd, sort_keys=False), encoding="utf-8")
        ok, msg = raises(ValueError, DeviceDesign.load, pb, contains=("inert",))
        check(f"composed + legacy block '{blk}' -> error containing 'inert'", ok, msg)
    be = copy.deepcopy(doc)
    be["design"]["composed"] = {}
    pe = td / "empty.yaml"
    pe.write_text(yaml.safe_dump(be, sort_keys=False), encoding="utf-8")
    check("empty composed block rejected", raises(ValueError, DeviceDesign.load, pe)[0])
    bu = copy.deepcopy(doc)
    bu["design"]["composed"]["bogus"] = 1
    pu = td / "unk.yaml"
    pu.write_text(yaml.safe_dump(bu, sort_keys=False), encoding="utf-8")
    check("composed block failing ComposedSpec.from_dict rejected at load",
          raises(ValueError, DeviceDesign.load, pu, contains=("bogus",))[0])

    # legacy save() byte-identity against the HEAD device.py
    head_src = None
    try:
        head_src = subprocess.run(["git", "show", "HEAD:fsim_core/device.py"], cwd=ROOT, check=True,
                                  capture_output=True).stdout.decode("utf-8")
    except Exception as e:  # noqa: BLE001
        print(f"note: git show HEAD failed ({e}); falling back to dict-level comparison")
    head_mod = None
    if head_src is not None:
        hp = td / "device_head.py"
        hp.write_text(head_src, encoding="utf-8")
        spec_ = importlib.util.spec_from_file_location("fsim_core._device_head", hp)
        head_mod = importlib.util.module_from_spec(spec_)
        sys.modules["fsim_core._device_head"] = head_mod
        spec_.loader.exec_module(head_mod)
    cards = sorted(p for p in (ROOT / "cards").glob("*.yaml")
                   if isinstance(yaml.safe_load(p.read_text(encoding="utf-8")), dict)
                   and "design" in yaml.safe_load(p.read_text(encoding="utf-8")))
    check("14 design cards found", len(cards) == 14, str(len(cards)))
    same = []
    for cp in cards:
        new = DeviceDesign.load(cp)
        a = td / "new.yaml"
        new.save(a)
        ta = a.read_text(encoding="utf-8")
        if head_mod is not None:
            b = td / "head.yaml"
            head_mod.DeviceDesign.load(cp).save(b)
            same.append((cp.name, ta == b.read_text(encoding="utf-8")))
        else:
            same.append((cp.name, "composed" not in yaml.safe_load(ta)["design"]))
    check(("save() text identical to HEAD device.py" if head_mod else "save() has no composed key")
          + f" for all {len(same)} design cards", all(ok for _, ok in same),
          str([n for n, ok in same if not ok]))

# ---------------------------------------------------------------- fsim_core_hash
import fsim_studio  # noqa: E402

h0 = fsim_studio.fsim_core_hash()
comp = ROOT / "fsim_core" / "compose"
probe_py = comp / "_hash_probe.py"
probe_yaml = comp / "data" / "_hash_probe.yaml"
try:
    probe_py.write_text("X = 1\n", encoding="utf-8")
    h1 = fsim_studio.fsim_core_hash()
finally:
    probe_py.unlink(missing_ok=True)
h2 = fsim_studio.fsim_core_hash()
try:
    probe_yaml.write_text("a: 1\n", encoding="utf-8")
    h3 = fsim_studio.fsim_core_hash()
finally:
    probe_yaml.unlink(missing_ok=True)
h4 = fsim_studio.fsim_core_hash()
check("fsim_core_hash changes when a fsim_core/compose/*.py file appears and restores after",
      h1 != h0 and h2 == h0)
check("fsim_core_hash changes when a fsim_core/compose/data/*.yaml file appears and restores after",
      h3 != h0 and h4 == h0)

n, total = sum(RESULTS), len(RESULTS)
print(f"{n}/{total} compose spec checks passed")
sys.exit(0 if n == total else 1)
