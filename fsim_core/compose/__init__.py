"""Composable device engine (dot x structure x drive x setup). See
.workers/devlib/p1-design.md. P1 skeleton: types, registries, spec validation
and an engine stub; later pieces fill the registries and the engine.

Importing this package never imports fsim_core.device.
"""
from .engine import evaluate_composed, evaluate_spec
from .params import Leaf, ParamSpec, central, corners, lhs, resolve_leaf, widest_tag
from .registry import (Drive, Entry, Material, Observable, Registry, Setup, Solver, Structure,
                       ensure_loaded, get, kinds, register, registry_payload)
from .spec import ComposedSpec, validate_spec
from .types import (DrivePlan, EmitterModel, Pulse, Rates, SetupChain, SolverResult,
                    StructureResponse)

__all__ = [
    "evaluate_composed", "evaluate_spec", "Leaf", "ParamSpec", "central", "corners", "lhs",
    "resolve_leaf", "widest_tag", "Drive", "Entry", "Material", "Observable", "Registry", "Setup",
    "Solver", "Structure", "ensure_loaded", "get", "kinds", "register", "registry_payload",
    "ComposedSpec", "validate_spec", "DrivePlan", "EmitterModel", "Pulse", "Rates", "SetupChain",
    "SolverResult", "StructureResponse",
]
