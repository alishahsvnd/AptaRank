<#
.SYNOPSIS
    Deploy AptaRank to the lab server from this Windows workstation.

.DESCRIPTION
    Pushes the committed code over SSH to a bare repository on the server,
    checks it out, installs or updates the environment, and restarts the
    dashboard.

    Deploys from git rather than copying files, so the server always runs an
    identifiable commit — the same one the paper cites — and an unfinished
    edit on this machine can never reach it.

.EXAMPLE
    .\deploy\deploy.ps1
    .\deploy\deploy.ps1 -Host H200 -SkipRestart
#>
[CmdletBinding()]
param(
    [string]$RemoteHost = "H200",
    [string]$AppDir     = "aptarank",
    [string]$DataDir    = "aptarank-data",
    [string]$Branch     = "main",
    [switch]$SkipRestart,
    [switch]$SkipInstall
)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
Set-Location $repo

function Step($message) { Write-Host "`n==> $message" -ForegroundColor Cyan }
function Ok($message)   { Write-Host "    ok  $message" -ForegroundColor Green }
function Warn($message) { Write-Host "    !!  $message" -ForegroundColor Yellow }

Step "Checking the working tree"
$dirty = git status --porcelain
if ($dirty) {
    Warn "Uncommitted changes will NOT be deployed:"
    $dirty -split "`n" | Select-Object -First 10 | ForEach-Object { Write-Host "        $_" }
    if ((Read-Host "    Continue anyway? [y/N]") -notmatch '^[Yy]') { exit 1 }
}
$commit = (git rev-parse --short HEAD).Trim()
Ok "deploying $Branch @ $commit"

Step "Preparing the server"
$prepare = @"
set -e
mkdir -p ~/$AppDir ~/$DataDir
if [ ! -d ~/$AppDir.git ]; then
    git init --bare -q ~/$AppDir.git
    echo '    created bare repository'
fi
"@
ssh $RemoteHost $prepare
Ok "remote repository ready"

Step "Pushing code"
$remoteUrl = "${RemoteHost}:$AppDir.git"
git push --force $remoteUrl "${Branch}:${Branch}" 2>&1 | ForEach-Object { Write-Host "    $_" }

$checkout = @"
set -e
cd ~/$AppDir
if [ ! -d .git ]; then
    git init -q .
    git remote add origin ~/$AppDir.git 2>/dev/null || true
fi
git remote set-url origin ~/$AppDir.git
git fetch -q origin $Branch
git checkout -q -B $Branch origin/$Branch
git rev-parse --short HEAD
"@
$deployed = (ssh $RemoteHost $checkout | Select-Object -Last 1).Trim()
Ok "server now at $deployed"
if ($deployed -ne $commit) { Warn "server commit differs from local ($commit)" }

if (-not $SkipInstall) {
    Step "Installing / updating the environment (this can take a few minutes)"
    ssh $RemoteHost "chmod +x ~/$AppDir/deploy/*.sh; APTARANK_APP_DIR=~/$AppDir APTARANK_DATA_DIR=~/$DataDir bash ~/$AppDir/deploy/server_install.sh"
}

if (-not $SkipRestart) {
    Step "Restarting the dashboard"
    ssh $RemoteHost "APTARANK_APP_DIR=~/$AppDir APTARANK_DATA_DIR=~/$DataDir bash ~/$AppDir/deploy/aptarank.sh restart"

    Step "Health check"
    $health = ssh $RemoteHost "sleep 3; curl -sS -o /dev/null -w '%{http_code}' http://127.0.0.1:8501/ || echo failed"
    if ($health -match "200") {
        Ok "dashboard responding"
    } else {
        Warn "health check returned '$health' — check: ssh $RemoteHost 'tail -40 ~/$DataDir/logs/dashboard.log'"
        exit 1
    }
}

Write-Host ""
Write-Host "Deployed." -ForegroundColor Green
Write-Host ""
Write-Host "  To use it from this machine, double-click:   deploy\connect.bat"
Write-Host "  or run:                                      ssh -N -L 8501:127.0.0.1:8501 $RemoteHost"
Write-Host "  then open:                                   http://localhost:8501"
Write-Host ""
