# qd-photon-sim

F-series quantum-dot single-photon-source simulator (Python 3, numpy/scipy). Branch `rt-edge-emitter` adds the room-temperature, electrically driven, edge-emitting InP-dot tier (see `../_goal/GOAL_PROMPT.md`, `../_goal/PROGRESS.md`, `../_goal/materials_research.md`, `../_goal/code_audit.md`).

## Conventions

- Every literature number in code carries a provenance tag in its comment or docstring: `[V]` verified against the cited paper, `[DR]` derived, `[E]` estimate/class range, `[A]` assumption. Cite the paper (author, journal, year).
- Legacy paths must stay bit-identical: `verify/verify_fsim.py` (51/51) and `verify/audit_physics.py` (23/23) must keep passing after any change to `fsim_core/`.
- New modules ship with a `verify/verify_<module>.py` that checks against numbers the code did not produce (published values), exits 0 iff all pass, and prints `N/N ... passed`.
- No pushes to the remote. No commits by workers.

## Test command

```
python verify/verify_fsim.py && python verify/audit_physics.py && python verify/verify_materials.py && python verify/verify_cw_g2.py
```

## Shell note for codex workers (2026-09-09)

On this machine the default PowerShell of codex's exec tool cannot start inside the sandbox: any command that lets the tool pick its own shell fails with `CreateProcessAsUserW failed: -1073283067` (it resolves `pwsh.exe` to the MSIX package path, which the sandbox token cannot execute). Whatever exec tool you have (`shell` or `exec_command`), ALWAYS name cmd.exe explicitly as the program, e.g. `C:\Windows\system32\cmd.exe /c type .workers\specsoo.md`, `C:\Windows\system32\cmd.exe /c python verifyerify_x.py`, `C:\Windows\system32\cmd.exe /c git diff 0871997 -- path`. If a command fails with a process-creation error, re-issue it through cmd.exe as above and continue; never stop the task because of it, and never conclude that no shell is available.
