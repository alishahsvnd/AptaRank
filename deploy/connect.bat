@echo off
REM ============================================================
REM  Open AptaRank (running on the lab server).
REM
REM  Double-click this file. It opens a secure connection to the
REM  server and then opens AptaRank in your browser.
REM
REM  Leave the window that appears open while you are using it.
REM  Closing it disconnects you (any analysis already running on
REM  the server keeps going).
REM ============================================================
setlocal

set REMOTE=H200
set PORT=8510

echo.
echo   Connecting to the AptaRank server...
echo.

REM -N: no remote shell, just the tunnel. The server's SSH keys are the
REM authentication, so there is no extra password to remember or share.
start "AptaRank connection - keep this window open" cmd /c ^
    "ssh -N -L %PORT%:127.0.0.1:%PORT% %REMOTE% || (echo. & echo   Could not connect to %REMOTE%. & echo   Check that you are on the lab network and that your SSH key is set up. & echo. & pause)"

echo   Waiting for the connection...
timeout /t 4 /nobreak >nul

start http://localhost:%PORT%

echo.
echo   AptaRank should now be open in your browser at:
echo       http://localhost:%PORT%
echo.
echo   If the page does not load, wait a few seconds and refresh.
echo.
timeout /t 6 /nobreak >nul
