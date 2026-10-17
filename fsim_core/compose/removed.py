"""Parameters that were REMOVED from a registry kind because they moved no output (R2 audit finding 10, 2026-10-13).

Rule of the product: a control in the Composer either changes an output in some context or is not shown. A saved spec or
a paper card that still carries a removed key keeps loading: validate_spec reports it as a WARNING (never an error) and
the engine ignores it (the engine reads only the parameters a kind declares).

{(axis, kind): {param: why}}; kind "*" = every kind of the axis.
"""
from __future__ import annotations

REMOVED_PARAMS = {
    ("material", "*"): {
        "fss_ueV": "the exciton fine-structure splitting is a polarisation quantum beat; the engine has no polarisation-"
                   "resolved or entanglement observable that depends on it (removed, R2 F10)",
    },
    ("drive", "electrical_pulsed"): {
        "mu": "every pulsed mechanism ignores mu (M-3 computes mu = eta I tau / e from the current, the M-4 turnstile fixes "
              "mu = 1); mu is the cw operating point (removed, R2 F10b)",
        "dg_inj_factor": "the injection-broadening factor of resonant-tunnelling injection has no consumer: no solver "
                         "models the injection linewidth (removed, R2 F10b)",
    },
    ("drive", "electrical_cw"): {
        "dg_inj_factor": "the injection-broadening factor of resonant-tunnelling injection has no consumer: no solver "
                         "models the injection linewidth (removed, R2 F10b)",
    },
}


def removed_reason(axis: str, kind: str, name: str):
    """Why `name` was removed from (axis, kind), or None if it was not."""
    for key in ((axis, kind), (axis, "*")):
        why = REMOVED_PARAMS.get(key, {}).get(name)
        if why:
            return why
    return None
