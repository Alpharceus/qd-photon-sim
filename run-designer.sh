#!/usr/bin/env bash
# FSIM device designer (Dear PyGui) -- launch from anywhere
cd "$(dirname "$0")"
exec .venv/bin/python fsim_gui/designer.py "$@"
