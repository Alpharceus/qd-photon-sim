"""Cards, META, presets, validation and save. Configuration only: every value
comes from cards/*.yaml, fsim_core.design_meta or fsim_core.presets."""
from __future__ import annotations

import dataclasses
import os
import re
import tempfile
from pathlib import Path

import yaml

from . import CARDS

NAME_RE = re.compile(r"^[A-Za-z0-9][A-Za-z0-9._-]{0,120}$")
# Saved designs (studio-p2a "Save as card"): lower-case slug, kept apart from the shipped
# cards/*.yaml (which scripts, gates and fsim_gui glob) in cards/studio/; FSIM_STUDIO_USER_CARDS
# overrides the folder (verify scripts use a temp dir).
SLUG_RE = re.compile(r"^[a-z0-9]+(?:-[a-z0-9]+)*$")
SLUG_MAX = 80
RESERVED_NAMES = frozenset({"unsaved", "studio"})
SAVED_BY = "fsim-studio"
USER_CARDS = Path(os.environ.get("FSIM_STUDIO_USER_CARDS") or (CARDS / "studio"))


class CardError(Exception):
    def __init__(self, status: int, message: str):
        super().__init__(message)
        self.status = status
        self.message = message


def _card_path(name: str) -> Path:
    if not NAME_RE.match(name or "") or name.endswith(".yaml"):
        raise CardError(400, f"invalid card name {name!r}")
    return CARDS / f"{name}.yaml"


def _read_yaml(path: Path) -> dict:
    return yaml.safe_load(path.read_text(encoding="utf-8")) or {}


def _find_card(name: str) -> Path:
    """A shipped card (cards/<name>.yaml) wins; otherwise a saved one (USER_CARDS)."""
    path = _card_path(name)
    if path.is_file():
        return path
    user = USER_CARDS / f"{name}.yaml"
    return user if user.is_file() else path


def _card_files():
    for p in sorted(CARDS.glob("*.yaml")):
        yield p, False
    if USER_CARDS.is_dir() and USER_CARDS.resolve() != CARDS.resolve():
        shipped = {p.stem for p in CARDS.glob("*.yaml")}
        for p in sorted(USER_CARDS.glob("*.yaml")):
            if p.stem not in shipped and NAME_RE.match(p.stem):
                yield p, True


def list_cards() -> list:
    out = []
    for p, saved in _card_files():
        try:
            doc = _read_yaml(p)
        except (OSError, yaml.YAMLError):
            continue
        meta = doc.get("meta") or {}
        is_design = isinstance(doc.get("design"), dict)
        platform = (doc["design"].get("platform", "legacy") if is_design else None)
        device = meta.get("device")
        out.append({
            "name": p.stem,
            "kind": "design" if is_design else "param",
            "platform": platform,
            "device": " ".join(str(device).split()) if device else None,
            "title": meta.get("name", p.stem),
            "role": meta.get("role"),
            **({"saved": True, "presets": meta.get("presets")} if saved else {}),
        })
    return out


def load_design(name: str):
    from fsim_core.device import DeviceDesign
    path = _find_card(name)
    if not path.is_file():
        raise CardError(404, f"no card {name!r}")
    doc = _read_yaml(path)
    if not isinstance(doc.get("design"), dict):
        raise CardError(422, f"card {name!r} is a parameter card, not a design card")
    try:
        return DeviceDesign.load(path)
    except (ValueError, TypeError) as exc:
        raise CardError(422, f"card {name!r} failed to load: {exc}") from exc


def design_to_dict(design) -> dict:
    return dataclasses.asdict(design)


def design_from_dict(dct: dict):
    """DeviceDesign from a JSON dict through DeviceDesign.load itself (the same
    validation and coercion as a card file): the dict is written as a
    one-off YAML document in the system temp dir and loaded."""
    from fsim_core.device import DeviceDesign
    if not isinstance(dct, dict):
        raise CardError(400, "design must be an object")
    fd, tmp = tempfile.mkstemp(prefix="fsim_studio_", suffix=".yaml")
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            yaml.safe_dump({"design": dct}, fh, sort_keys=False)
        try:
            return DeviceDesign.load(tmp)
        except (ValueError, TypeError, KeyError) as exc:
            raise CardError(422, f"invalid design: {exc}") from exc
    finally:
        try:
            os.unlink(tmp)
        except OSError:
            pass


def _meta_for(path: str) -> dict | None:
    from fsim_core import design_meta
    m = design_meta.META.get(path)
    if m is None:
        return None
    band = [m["lo"], m["hi"]] if "lo" in m and "hi" in m else None
    out = {"tag": m.get("tag"), "unit": m.get("unit"), "band": band,
           "source": m.get("source"), "label": path.split(".", 1)[1]}
    if "choices" in m:
        out["choices"] = list(m["choices"])
    return out


def card_meta(design_dict: dict) -> dict:
    """META annotations for every dotted path present in the design dict."""
    out = {}
    for block, fields in design_dict.items():
        if not isinstance(fields, dict):
            continue
        for field_name in fields:
            path = f"{block}.{field_name}"
            m = _meta_for(path)
            if m is not None:
                out[path] = m
    return out


def get_card(name: str) -> dict:
    d = load_design(name)
    dd = design_to_dict(d)
    out = {"name": name, "design": dd, "meta": card_meta(dd), "platform": d.platform}
    if _find_card(name).parent != CARDS:
        out["saved"] = True
    return out


def _shown_path(path: Path) -> str:
    try:
        return str(path.resolve().relative_to(CARDS.parent.resolve())).replace("\\", "/")
    except ValueError:
        return str(path)


def save_card(name: str, design_dict: dict, presets_used: dict | None = None) -> dict:
    """Save as card (studio-p2a): a lower-case slug; a shipped card name is refused with 409
    (new names are allowed, and a design saved here earlier may be overwritten). The file goes
    to USER_CARDS through DeviceDesign.save, then is re-read through DeviceDesign.load (the
    round-trip must reproduce the design) before the Studio marker is added to its meta block."""
    if not isinstance(name, str) or not SLUG_RE.match(name) or len(name) > SLUG_MAX:
        raise CardError(400, f"invalid card name {name!r}: use a lower-case slug "
                             f"(letters, digits, single hyphens; at most {SLUG_MAX} characters)")
    if name in RESERVED_NAMES:
        raise CardError(400, f"card name {name!r} is reserved")
    if (CARDS / f"{name}.yaml").exists():
        raise CardError(409, f"refusing to overwrite shipped card {name!r}; save under a new name")
    path = USER_CARDS / f"{name}.yaml"
    if path.exists():
        try:
            prior = (_read_yaml(path).get("meta") or {}).get("saved_by")
        except (OSError, yaml.YAMLError):
            prior = None
        if prior != SAVED_BY:
            raise CardError(409, f"refusing to overwrite {name!r}: not a design saved by FSIM Studio")
    if presets_used is not None and not isinstance(presets_used, dict):
        raise CardError(400, "presets must be an object {axis: preset name}")
    d = design_from_dict(design_dict)
    d.name = name
    USER_CARDS.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(f".{os.getpid()}.tmp")
    try:
        d.save(tmp)
        from fsim_core.device import DeviceDesign
        back = DeviceDesign.load(tmp)
        if design_to_dict(back) != design_to_dict(d):
            raise CardError(422, f"design {name!r} does not round-trip through DeviceDesign.save/load")
        doc = _read_yaml(tmp)
        doc.setdefault("meta", {})["saved_by"] = SAVED_BY
        if presets_used:
            doc["meta"]["presets"] = {str(k): str(v) for k, v in presets_used.items() if v}
        tmp.write_text(yaml.safe_dump(doc, sort_keys=False), encoding="utf-8")
        os.replace(tmp, path)
    finally:
        tmp.unlink(missing_ok=True)   # no <name>.<pid>.tmp is ever left behind (also after os.replace: no-op)
    return {"name": name, "path": _shown_path(path), "saved": True}


def meta_payload() -> dict:
    from fsim_core import design_meta, presets
    return {
        "META": design_meta.META,
        "ENV_DEFAULTS": {k: list(v) for k, v in design_meta.ENV_DEFAULTS.items()},
        "presets": {
            "dot": presets.DOT_PRESETS,
            "template": presets.TEMPLATE_PRESETS,
            "cavity": presets.CAVITY_PRESETS,
            "drive": presets.DRIVE_PRESETS,
            "injection": presets.INJECTION_PRESETS,
        },
    }


def apply_presets(body: dict) -> dict:
    """Same call order as the designer's apply_presets(): dot, template,
    cavity, drive presets via preset_device, then the injection preset."""
    from fsim_core import presets
    try:
        d = presets.preset_device(body["dot"], body["template"], body["cavity"], body["drive"],
                                  name=body.get("name", "preset-device"))
        inj = body.get("injection")
        if inj:
            presets.apply_injection_preset(d, inj)
    except KeyError as exc:
        raise CardError(400, f"unknown or missing preset key: {exc}") from exc
    dd = design_to_dict(d)
    return {"design": dd, "meta": card_meta(dd)}


def f8_floor(design) -> float:
    """F8 loading-domain floor mu >= 1 - F_eff at the design's own drive
    fields, F_eff from fsim_core.loading.f8b_thin_fano."""
    from fsim_core.loading import f8b_thin_fano
    F_eff = f8b_thin_fano(design.drive.eta_capture, design.drive.F_p)
    return 1.0 - F_eff


def validate(design_dict: dict) -> dict:
    from fsim_core import design_meta
    d = design_from_dict(design_dict)
    violations = []
    for path, m in design_meta.META.items():
        if "lo" not in m or "hi" not in m or m.get("unit") in ("bool", "enum"):
            continue
        block, field_name = path.split(".", 1)
        value = getattr(getattr(d, block, None), field_name, None)
        if isinstance(value, bool) or not isinstance(value, (int, float)):
            continue
        if value < m["lo"] or value > m["hi"]:
            violations.append({"path": path, "value": value, "band": [m["lo"], m["hi"]]})
    return {"violations": violations,
            "default_ranged": {k: list(v) for k, v in design_meta.default_ranged(d).items()},
            "f8_floor": f8_floor(d)}
