#!/usr/bin/env bash
# FSIM Studio (local web app over fsim_core) -- launch from anywhere
#   ./run-studio.sh            native window (pywebview)
#   ./run-studio.sh --browser  default browser instead
cd "$(dirname "$0")"
exec .venv/bin/python -m fsim_studio "$@"
