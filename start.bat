@echo off
REM ============================================================
REM  AptaRank - double-click this file to start.
REM
REM  Deliberately does not call "python": on many machines that
REM  resolves to an interpreter without the scientific packages,
REM  and the resulting error is impossible to interpret. The
REM  launcher below finds the right interpreter itself.
REM ============================================================
setlocal
cd /d "%~dp0"

if exist ".venv\Scripts\python.exe" (
    ".venv\Scripts\python.exe" "scripts\start.py" %*
    goto :done
)

where py >nul 2>nul
if %errorlevel%==0 (
    py -3 "scripts\start.py" %*
    goto :done
)

echo.
echo   AptaRank needs Python 3.10 or newer, and could not find it.
echo.
echo   Install it from https://www.python.org/downloads/
echo   During setup, tick "Add python.exe to PATH".
echo   Then double-click this file again.
echo.

:done
echo.
pause
