<#
.SYNOPSIS
    Deploy AptaRank to the lab server from this Windows workstation.

.DESCRIPTION
    Pushes the committed code over SSH to a bare repository on the server,
    checks it out, installs or updates the environment, and restarts the
    dashboard.

    Deploys from git rather than copying files, so the server always runs an
    identifiable commit (the same one the paper cites), and an unfinished edit
    on this machine can never reach it.

    Kept to plain ASCII on purpose: Windows PowerShell 5.1 reads a BOM-less
    UTF-8 script as ANSI, and a single stray character mangles the parse with
    an error that points at the wrong line.

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
    [int]$Port          = 8510,
    [switch]$SkipRestart,
    [switch]$SkipInstall,
    [switch]$Force        # deploy despite uncommitted changes, without prompting
)

$ErrorActionPreference = "Stop"
$repo = Split-Path -Parent $PSScriptRoot
Set-Location $repo

function Step($message) { Write-Host "`n==> $message" -ForegroundColor Cyan }
function Ok($message)   { Write-Host "    ok  $message" -ForegroundColor Green }
function Warn($message) { Write-Host "    !!  $message" -ForegroundColor Yellow }

function Invoke-Remote {
    <#
      Run a script on the server.

      PowerShell here-strings carry CRLF line endings; bash receives the
      trailing CR as part of each command and fails in ways that look like the
      script is fine ("set -e" reported as an invalid option, directories
      created with an invisible CR in the name). Normalise before sending.
    #>
    param([string]$Script)
    ssh $RemoteHost ($Script -replace "`r`n", "`n")
}

Step "Checking the working tree"
$dirty = git status --porcelain
if ($dirty) {
    Warn "Uncommitted changes will NOT be deployed:"
    $dirty -split "`n" | Select-Object -First 10 | ForEach-Object { Write-Host "        $_" }
    if ($Force) {
        Warn "-Force given: deploying the last commit regardless"
    } elseif ([Environment]::UserInteractive -and -not [Console]::IsInputRedirected) {
        if ((Read-Host "    Continue anyway? [y/N]") -notmatch '^[Yy]') { exit 1 }
    } else {
        Warn "Refusing to deploy with a dirty tree in a non-interactive shell."
        Warn "Commit the changes, or re-run with -Force."
        exit 1
    }
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
Invoke-Remote $prepare
Ok "remote repository ready"

Step "Pushing code"
$remoteUrl = "${RemoteHost}:$AppDir.git"
git push --force $remoteUrl "${Branch}:${Branch}" 2>&1 | ForEach-Object { Write-Host "    $_" }

# reset --hard, not checkout: the working tree is disposable code, and the
# server has no business carrying local edits. core.fileMode=false stops the
# executable bit (set when the install script runs chmod) from registering as a
# modification and blocking the next deploy. Untracked paths - .venv and the
# generated configs/server.yaml - are deliberately left alone.
$checkout = @"
set -e
cd ~/$AppDir
if [ ! -d .git ]; then
    git init -q .
    git remote add origin ~/$AppDir.git 2>/dev/null || true
fi
git config core.fileMode false
git remote set-url origin ~/$AppDir.git
git fetch -q origin $Branch
git reset -q --hard origin/$Branch
git rev-parse --short HEAD
"@
$deployed = (Invoke-Remote $checkout | Select-Object -Last 1).Trim()
Ok "server now at $deployed"
if ($deployed -ne $commit) { Warn "server commit differs from local ($commit)" }

if (-not $SkipInstall) {
    Step "Installing / updating the environment (this can take a few minutes)"
    ssh $RemoteHost "chmod +x ~/$AppDir/deploy/*.sh; APTARANK_APP_DIR=~/$AppDir APTARANK_DATA_DIR=~/$DataDir APTARANK_PORT=$Port bash ~/$AppDir/deploy/server_install.sh"
}

if (-not $SkipRestart) {
    Step "Restarting the dashboard"
    ssh $RemoteHost "APTARANK_APP_DIR=~/$AppDir APTARANK_DATA_DIR=~/$DataDir APTARANK_PORT=$Port bash ~/$AppDir/deploy/aptarank.sh restart"

    Step "Health check"
    # Check that OUR process is serving the port, not merely that something is.
    # Another user's Streamlit was already on 8501, and an HTTP-only check
    # happily reported their app as a successful deployment of ours.
    $probe = @"
pid=`$(cat ~/$DataDir/aptarank.pid 2>/dev/null || echo none)
kill -0 "`$pid" 2>/dev/null || { echo "dead:`$pid"; exit 0; }
owner=`$(ss -ltnp "sport = :$Port" 2>/dev/null | grep -o 'pid=[0-9]*' | head -1 | cut -d= -f2)
code=`$(curl -sS -o /dev/null -w '%{http_code}' http://127.0.0.1:$Port/ 2>/dev/null || echo 000)
echo "pid=`$pid owner=`$owner code=`$code"
"@
    $result = (Invoke-Remote $probe | Select-Object -Last 1)
    if ($result -match 'pid=(\d+) owner=(\d*) code=(\d+)' -and
        $Matches[3] -eq "200" -and $Matches[1] -eq $Matches[2]) {
        Ok "dashboard responding on port $Port (pid $($Matches[1]))"
    } else {
        Warn "health check failed: $result"
        Warn "Check: ssh $RemoteHost 'tail -40 ~/$DataDir/logs/dashboard.log'"
        exit 1
    }
}

Write-Host ""
Write-Host "Deployed." -ForegroundColor Green
Write-Host ""
Write-Host "  To use it from this machine, double-click:   deploy\connect.bat"
Write-Host "  or run:                                      ssh -N -L $Port`:127.0.0.1:$Port $RemoteHost"
Write-Host "  then open:                                   http://localhost:$Port"
Write-Host ""
