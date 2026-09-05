# qd-photon-sim

F-series quantum-dot single-photon-source simulator (Python 3, numpy/scipy). Branch `rt-edge-emitter` adds the room-temperature, electrically driven, edge-emitting InP-dot tier (see `../_goal/GOAL_PROMPT.md`, `../_goal/PROGRESS.md`, `../_goal/materials_research.md`, `../_goal/code_audit.md`).

## Conventions

- Every literature number in code carries a provenance tag in its comment or docstring: `[V]` verified against the cited paper, `[DR]` derived, `[E]` estimate/class range, `[A]` assumption. Cite the paper (author, journal, year).
- Legacy paths must stay bit-identical: `verify/verify_fsim.py` (51/51) and `verify/audit_physics.py` (23/23) must keep passing after any change to `fsim_core/`.
- New modules ship with a `verify/verify_<module>.py` that checks against numbers the code did not produce (published values), exits 0 iff all pass, and prints `N/N ... passed`.
- No pushes to the remote. No commits by workers.

## Sandbox note for workers

`verify/verify_fsim.py` creates a `tempfile.TemporaryDirectory()`; inside a read-only or workspace-write sandbox this can raise `PermissionError`. If that is the ONLY failure, report `TESTS: pass (verify_fsim skipped: sandbox temp dir)` and `STATUS: DONE`; the orchestrator re-runs the full test command outside the sandbox.

## Test command

```
python verify/verify_fsim.py && python verify/audit_physics.py && python verify/verify_materials.py && python verify/verify_cw_g2.py && python verify/verify_dot_levels.py && python verify/verify_waveguide.py && python verify/verify_transport.py
```
