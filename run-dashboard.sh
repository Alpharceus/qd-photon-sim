#!/usr/bin/env bash
# FSIM validation dashboard (Streamlit) -- launch from anywhere
cd "$(dirname "$0")"
exec .venv/bin/python -m streamlit run fsim_gui/app.py "$@"
