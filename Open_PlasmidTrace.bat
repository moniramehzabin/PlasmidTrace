@echo off
cd /d "%~dp0"
py -3 browser_app.py
if errorlevel 1 (
 echo.
 echo PlasmidTrace could not start. Keep the complete extracted app folder together.
 echo Windows Python 3.10 or newer is required.
 pause
)
