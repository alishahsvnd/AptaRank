#!/usr/bin/env bash
# Start, stop and inspect the AptaRank dashboard on a server without root.
#
#   ./deploy/aptarank.sh start | stop | restart | status | logs | tunnel-help
#
# The dashboard binds to 127.0.0.1 only. Streamlit has no authentication of its
# own, and this service accepts file uploads and launches subprocesses, so it is
# not something to expose on a network — even an internal one. Users reach it
# through an SSH tunnel, which means the server's existing SSH keys *are* the
# authentication, per person, with no shared password to leak.
#
# See deploy/README.md for the reverse-proxy upgrade path if browser access
# without a tunnel is ever needed.

set -euo pipefail

APP_DIR="${APTARANK_APP_DIR:-$HOME/aptarank}"
DATA_DIR="${APTARANK_DATA_DIR:-$HOME/aptarank-data}"
# 8501 and 8502 are Streamlit's defaults and were already taken on this machine
# by someone else's app - which also meant an early health check was cheerfully
# reporting *their* service as ours. Pick something unlikely to collide.
PORT="${APTARANK_PORT:-8510}"
BIND="${APTARANK_BIND:-127.0.0.1}"

PID_FILE="$DATA_DIR/aptarank.pid"
LOG_FILE="$DATA_DIR/logs/dashboard.log"
PY="$APP_DIR/.venv/bin/python"

# Shared-machine resource policy, read by dashboard/jobs.py.
export APTARANK_MAX_CONCURRENT_JOBS="${APTARANK_MAX_CONCURRENT_JOBS:-2}"
export APTARANK_WORKERS_PER_JOB="${APTARANK_WORKERS_PER_JOB:-16}"
export APTARANK_DATA_DIR="$DATA_DIR"
# Layered over configs/default.yaml so every run writes into the data
# directory rather than into the checkout that the next deploy replaces.
export APTARANK_CONFIG="${APTARANK_CONFIG:-$APP_DIR/configs/server.yaml}"
# AptaRank is CPU-only. Never let it probe or claim a GPU on a box whose GPUs
# other people are using.
export CUDA_VISIBLE_DEVICES=""
export PATH="$HOME/.local/bin:$PATH"

running_pid() {
    [ -f "$PID_FILE" ] || return 1
    local pid
    pid="$(cat "$PID_FILE" 2>/dev/null || true)"
    [ -n "$pid" ] || return 1
    kill -0 "$pid" 2>/dev/null || return 1
    echo "$pid"
}

case "${1:-status}" in

start)
    if pid="$(running_pid)"; then
        echo "Already running (pid $pid) on http://$BIND:$PORT"
        exit 0
    fi
    [ -x "$PY" ] || { echo "Not installed: run deploy/server_install.sh first"; exit 1; }
    if ss -ltn "sport = :$PORT" 2>/dev/null | grep -q LISTEN; then
        echo "Port $PORT is already in use by something else on this machine."
        echo "Set APTARANK_PORT to a free port and try again."
        exit 1
    fi
    mkdir -p "$DATA_DIR/logs"

    # nice: this is a background service on someone else's compute box.
    nohup nice -n 10 "$PY" -m streamlit run "$APP_DIR/dashboard/streamlit_app.py" \
        --server.address "$BIND" \
        --server.port "$PORT" \
        --server.headless true \
        --server.maxUploadSize 200 \
        --server.enableXsrfProtection true \
        --browser.gatherUsageStats false \
        -- --runs-dir "$DATA_DIR/runs" \
        >>"$LOG_FILE" 2>&1 &

    echo $! > "$PID_FILE"
    sleep 4
    if pid="$(running_pid)"; then
        echo "Started (pid $pid)"
        echo "  listening on   http://$BIND:$PORT   (loopback only, by design)"
        echo "  job slots      $APTARANK_MAX_CONCURRENT_JOBS × $APTARANK_WORKERS_PER_JOB workers"
        echo "  data           $DATA_DIR"
        echo "  log            $LOG_FILE"
        echo
        echo "From a workstation:  ssh -N -L $PORT:127.0.0.1:$PORT <this-host>"
        echo "then open            http://localhost:$PORT"
    else
        echo "Failed to start. Last lines of $LOG_FILE:"
        tail -20 "$LOG_FILE"
        exit 1
    fi
    ;;

stop)
    if pid="$(running_pid)"; then
        # The dashboard's child analyses are deliberately left alone: they are
        # detached, they write their own results, and killing a half-finished
        # scientific run to restart a web page would be the wrong trade.
        kill "$pid"
        rm -f "$PID_FILE"
        echo "Stopped (pid $pid). Any analyses already running will finish."
    else
        echo "Not running."
    fi
    ;;

restart)
    "$0" stop || true
    sleep 2
    "$0" start
    ;;

status)
    if pid="$(running_pid)"; then
        echo "running   pid $pid   http://$BIND:$PORT"
    else
        echo "stopped"
    fi
    echo "app       $APP_DIR ($(cd "$APP_DIR" && git rev-parse --short HEAD 2>/dev/null || echo 'not a git checkout'))"
    echo "data      $DATA_DIR"
    if [ -d "$DATA_DIR/runs/jobs" ]; then
        echo "analyses  $(find "$DATA_DIR/runs/jobs" -maxdepth 1 -mindepth 1 -type d | wc -l) on disk"
    fi
    echo "fpocket   $(command -v fpocket || echo 'not available — Tier 2 needs prepared bundles')"
    ;;

logs)
    tail -n "${2:-60}" -f "$LOG_FILE"
    ;;

tunnel-help)
    cat <<EOF

To use AptaRank from your own computer:

  1. Open a terminal (PowerShell on Windows) and run:

       ssh -N -L $PORT:127.0.0.1:$PORT $(whoami)@$(hostname -I 2>/dev/null | awk '{print $1}')

     Leave that window open.

  2. Open a browser at:

       http://localhost:$PORT

Windows users can double-click deploy/connect.bat instead, which does both.

EOF
    ;;

*)
    echo "usage: $0 {start|stop|restart|status|logs|tunnel-help}"
    exit 2
    ;;
esac
