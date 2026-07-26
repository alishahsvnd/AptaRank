#!/usr/bin/env bash
# Install or update AptaRank on a Linux server, without root.
#
# The target machine is shared and we do not have sudo, so everything lands
# under $HOME: the Python environment, the fpocket binary built from source,
# and all mutable state. Nothing outside the two directories below is touched.
#
#   ~/aptarank        the code checkout (replaced on every deploy)
#   ~/aptarank-data   runs, uploads, caches, reference libraries (never replaced)
#
# Keeping data outside the checkout is what makes a redeploy safe: the code can
# be thrown away and restored from git; a biologist's uploaded sequences cannot.
#
# Idempotent — safe to run repeatedly.

set -euo pipefail

APP_DIR="${APTARANK_APP_DIR:-$HOME/aptarank}"
DATA_DIR="${APTARANK_DATA_DIR:-$HOME/aptarank-data}"
VENV_DIR="$APP_DIR/.venv"
LOCAL_BIN="$HOME/.local/bin"
PYTHON="${APTARANK_PYTHON:-python3}"

log()  { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
warn() { printf '\033[33m    !! %s\033[0m\n' "$*"; }
ok()   { printf '\033[32m    ok  %s\033[0m\n' "$*"; }

log "AptaRank server install"
echo "    app        $APP_DIR"
echo "    data       $DATA_DIR"
echo "    python     $($PYTHON --version 2>&1)"
echo "    host       $(hostname), $(nproc) cores"

# -- 1. persistent data directories -------------------------------------

log "Preparing data directories"
for sub in runs uploads cache/corpus cache/calibration cache/targets cache/structures \
           data/corpus data/libraries logs; do
    mkdir -p "$DATA_DIR/$sub"
done
chmod 700 "$DATA_DIR"          # uploaded sequences are unpublished data
ok "$DATA_DIR"

# -- 2. python environment ----------------------------------------------

log "Python environment"
if [ ! -x "$VENV_DIR/bin/python" ]; then
    "$PYTHON" -m venv "$VENV_DIR"
    ok "created $VENV_DIR"
else
    ok "reusing $VENV_DIR"
fi
PY="$VENV_DIR/bin/python"
"$PY" -m pip install --quiet --upgrade pip setuptools wheel cython
ok "build tools"

# -- 3. ushuffle FIRST, because it never installs cleanly ----------------
#
# It is a hard dependency of AptaRank, so installing the package first would
# make pip try to build it under build isolation — where Cython is absent, the
# shipped pre-generated C is used, and the build fails on any Python >= 3.9
# (tp_print on 3.9-3.11, longintrepr.h on 3.12). Building it here, with
# --no-build-isolation, lets setup.py re-cythonize the .pyx against the current
# interpreter; the package install below then sees the requirement satisfied.

if ! "$PY" -c "import ushuffle" 2>/dev/null; then
    log "Building ushuffle from source"
    work="$(mktemp -d)"
    "$PY" -m pip download ushuffle --no-binary :all: --no-deps -d "$work" --quiet
    tar -xzf "$work"/ushuffle-*.tar.gz -C "$work"
    if "$PY" -m pip install --quiet --no-build-isolation "$work"/ushuffle-*/; then
        ok "ushuffle built"
    else
        warn "ushuffle build failed — the shuffled-control check will be unavailable"
    fi
    rm -rf "$work"
else
    ok "ushuffle present"
fi

log "Installing AptaRank"
"$PY" -m pip install --quiet -e "$APP_DIR[dashboard]"
ok "AptaRank installed"

# -- 4. fpocket, so Tier 2 works on real structures ----------------------

if ! command -v fpocket >/dev/null 2>&1 && [ ! -x "$LOCAL_BIN/fpocket" ]; then
    log "Building fpocket from source"
    mkdir -p "$LOCAL_BIN"
    work="$(mktemp -d)"
    git clone --depth 1 https://github.com/Discngine/fpocket.git "$work" >/dev/null 2>&1
    # Serial build: fpocket's makefile races on -j and fails compiling qhull
    # before its headers are in place.
    if (cd "$work" && make >"$DATA_DIR/logs/fpocket_build.log" 2>&1); then
        install -m 0755 "$work"/bin/fpocket "$LOCAL_BIN/fpocket"
        install -m 0755 "$work"/bin/dpocket "$LOCAL_BIN/dpocket" 2>/dev/null || true
        ok "fpocket -> $LOCAL_BIN/fpocket ($(cd "$work" && git rev-parse --short HEAD))"
    else
        warn "fpocket build failed; see $DATA_DIR/logs/fpocket_build.log"
        warn "Tier 2 will only work with target bundles built elsewhere."
    fi
    rm -rf "$work"
else
    ok "fpocket present: $(command -v fpocket || echo "$LOCAL_BIN/fpocket")"
fi

# -- 5. server configuration ---------------------------------------------

log "Server configuration"
cat > "$APP_DIR/configs/server.yaml" <<YAML
# Written by deploy/server_install.sh — edit deploy settings there, not here.
#
# Resource policy for a shared machine. AptaRank is CPU-bound RNA folding and
# will happily saturate every core it is given; on a box other people are using,
# that is the difference between a useful service and an antisocial one.

corpus:
  cache_dir: $DATA_DIR/cache/corpus

tier1:
  parallel:
    workers: ${APTARANK_WORKERS_PER_JOB:-16}

tier2:
  bundle_dir: $DATA_DIR/cache/targets
  structure_cache_dir: $DATA_DIR/cache/structures
  calibration:
    cache_dir: $DATA_DIR/cache/calibration

output:
  dir: $DATA_DIR/runs
YAML
ok "configs/server.yaml"

# -- 6. verify ------------------------------------------------------------

log "Checking the installation"
"$PY" - <<'PYCODE'
import importlib, shutil, sys
missing = []
for name in ("aptarank", "RNA", "forgi", "ushuffle", "Bio", "streamlit", "altair"):
    try:
        importlib.import_module(name)
    except Exception as exc:
        missing.append(f"{name} ({exc.__class__.__name__})")
print("    fpocket   ", shutil.which("fpocket") or "not on PATH")
if missing:
    print("    MISSING   ", ", ".join(missing))
    sys.exit(1)
print("    imports    all present")
PYCODE

if [ -d "$APP_DIR/tests" ]; then
    if "$PY" -m pytest -q "$APP_DIR/tests" >"$DATA_DIR/logs/tests.log" 2>&1; then
        ok "test suite passed ($(grep -oE '[0-9]+ passed' "$DATA_DIR/logs/tests.log" | tail -1))"
    else
        warn "tests failed — see $DATA_DIR/logs/tests.log"
        warn "not starting a service in this state"
        exit 1
    fi
fi

log "Done"
echo "    Start it with:  $APP_DIR/deploy/aptarank.sh start"
