@echo off
rem FSIM Studio (local web app over fsim_core) -- launch from anywhere
rem   run-studio.bat            native window (pywebview)
rem   run-studio.bat --browser  default browser instead
cd /d "%~dp0"
python -m fsim_studio %*
