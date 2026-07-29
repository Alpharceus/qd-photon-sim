@echo off
rem FSIM validation dashboard (Streamlit) -- launch from anywhere
cd /d "%~dp0"
streamlit run fsim_gui/app.py
