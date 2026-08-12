#!/usr/bin/env bash
# FSIM full verification pass -- all suites, serial (house rule), from anywhere.
# Exit 0 iff every suite passes; prints one summary line per suite.
cd "$(dirname "$0")"
PY=.venv/bin/python
if [ ! -x "$PY" ]; then
    echo "error: $PY not found -- create the venv first:" >&2
    echo "  python3 -m venv .venv && .venv/bin/pip install -r requirements.txt" >&2
    exit 1
fi
suites=(verify_fsim verify4 verify_sde verify_drive verify_gf verify_dbr
        audit_physics gate_v11_gui gate_d1 gate_d2 gate_d3 gate_spec)
fail=0
for s in "${suites[@]}"; do
    printf '%-14s ' "$s"
    if out=$("$PY" "verify/$s.py" 2>&1); then
        echo "$out" | tail -1
    else
        fail=1
        echo "FAILED"
        echo "$out" | tail -15
    fi
done
exit $fail
