#!/usr/bin/env bash
# Give colleagues browser access to AptaRank, with a login, without root.
#
#   ./deploy/enable_access.sh add <username>      create or reset a login
#   ./deploy/enable_access.sh start               serve on the lab network
#   ./deploy/enable_access.sh stop
#   ./deploy/enable_access.sh status
#   ./deploy/enable_access.sh users
#
# What this sets up:
#
#   biologist's browser --HTTPS--> Caddy :8443 --HTTP--> Streamlit 127.0.0.1:8510
#                                  (per-person login)     (still loopback only)
#
# Streamlit itself stays bound to loopback and gains no authentication; the
# proxy in front of it is what makes access safe to hand out. Anyone who can
# reach the port must present a password, the traffic is encrypted, and each
# person has their own credential, so one can be revoked without disturbing the
# others. Handing out a single shared password would undo most of that.
#
# The certificate comes from Caddy's own internal authority, so the first visit
# shows a browser warning. That is the honest cost of not having an institutional
# certificate; see the note printed by `start`, and prefer a real certificate if
# your institution issues them.

set -euo pipefail

APP_DIR="${APTARANK_APP_DIR:-$HOME/aptarank}"
DATA_DIR="${APTARANK_DATA_DIR:-$HOME/aptarank-data}"
UPSTREAM_PORT="${APTARANK_PORT:-8510}"
LISTEN_PORT="${APTARANK_PUBLIC_PORT:-8443}"

CADDY="$HOME/.local/bin/caddy"
CADDYFILE="$DATA_DIR/Caddyfile"
USERS_FILE="$DATA_DIR/users.conf"          # username:bcrypt-hash, one per line
PID_FILE="$DATA_DIR/caddy.pid"
LOG_FILE="$DATA_DIR/logs/caddy.log"
CADDY_VERSION="${CADDY_VERSION:-2.8.4}"

log()  { printf '\n\033[1m==> %s\033[0m\n' "$*"; }
warn() { printf '\033[33m    !! %s\033[0m\n' "$*"; }
ok()   { printf '\033[32m    ok  %s\033[0m\n' "$*"; }

ensure_caddy() {
    [ -x "$CADDY" ] && return 0
    log "Fetching Caddy $CADDY_VERSION"
    mkdir -p "$HOME/.local/bin"
    local url="https://github.com/caddyserver/caddy/releases/download/v${CADDY_VERSION}/caddy_${CADDY_VERSION}_linux_amd64.tar.gz"
    local work
    work="$(mktemp -d)"
    curl -fsSL "$url" -o "$work/caddy.tar.gz"
    tar -xzf "$work/caddy.tar.gz" -C "$work" caddy
    install -m 0755 "$work/caddy" "$CADDY"
    rm -rf "$work"
    ok "$($CADDY version | head -1)"
}

write_caddyfile() {
    mkdir -p "$DATA_DIR/logs"
    touch "$USERS_FILE"
    chmod 600 "$USERS_FILE"

    local block=""
    while IFS=: read -r user hash; do
        [ -n "${user:-}" ] || continue
        block+="        $user $hash"$'\n'
    done < "$USERS_FILE"

    if [ -z "$block" ]; then
        warn "no logins defined yet — run: $0 add <username>"
        return 1
    fi

    cat > "$CADDYFILE" <<EOF
# Written by deploy/enable_access.sh. Edit logins with that script, not here.
{
    admin off
    storage file_system $DATA_DIR/caddy-storage
}

:$LISTEN_PORT {
    tls internal

    basic_auth {
$block    }

    # Streamlit talks over websockets and streams progress for the length of an
    # analysis; without generous timeouts the page goes blank mid-run.
    reverse_proxy 127.0.0.1:$UPSTREAM_PORT {
        transport http {
            read_timeout 24h
            write_timeout 24h
        }
    }

    request_body {
        max_size 200MB
    }

    # Deliberately minimal logging: request paths on this service can carry
    # nothing sensitive, but uploaded sequences are unpublished data and have no
    # business in a log file.
    log {
        output file $DATA_DIR/logs/caddy_access.log
        format console
    }
}
EOF
}

case "${1:-status}" in

add)
    user="${2:-}"
    [ -n "$user" ] || { echo "usage: $0 add <username>"; exit 2; }
    ensure_caddy
    mkdir -p "$DATA_DIR"
    touch "$USERS_FILE"; chmod 600 "$USERS_FILE"

    # A generated password, not one chosen under pressure at a keyboard.
    password="$(tr -dc 'A-Za-z0-9' </dev/urandom | head -c 16)"
    hash="$("$CADDY" hash-password --plaintext "$password")"

    grep -v "^$user:" "$USERS_FILE" > "$USERS_FILE.tmp" 2>/dev/null || true
    mv "$USERS_FILE.tmp" "$USERS_FILE"
    echo "$user:$hash" >> "$USERS_FILE"
    chmod 600 "$USERS_FILE"

    echo
    echo "  Login created. Send these to $user (and only to $user):"
    echo
    echo "      address    https://$(hostname -I | awk '{print $1}'):$LISTEN_PORT"
    echo "      username   $user"
    echo "      password   $password"
    echo
    echo "  The password is not stored anywhere in readable form. If it is lost,"
    echo "  run this command again to issue a new one."
    echo
    if [ -f "$PID_FILE" ]; then
        echo "  Apply it with:  $0 restart"
    else
        echo "  Start serving with:  $0 start"
    fi
    ;;

users)
    if [ -s "$USERS_FILE" ]; then
        echo "Logins:"
        cut -d: -f1 "$USERS_FILE" | sed 's/^/    /'
    else
        echo "No logins defined. Create one with: $0 add <username>"
    fi
    ;;

remove)
    user="${2:-}"
    [ -n "$user" ] || { echo "usage: $0 remove <username>"; exit 2; }
    grep -v "^$user:" "$USERS_FILE" > "$USERS_FILE.tmp" 2>/dev/null || true
    mv "$USERS_FILE.tmp" "$USERS_FILE"
    chmod 600 "$USERS_FILE"
    echo "Removed $user. Apply with: $0 restart"
    ;;

start)
    ensure_caddy
    if [ -f "$PID_FILE" ] && kill -0 "$(cat "$PID_FILE")" 2>/dev/null; then
        echo "Already serving on port $LISTEN_PORT (pid $(cat "$PID_FILE"))"
        exit 0
    fi
    write_caddyfile || exit 1

    if ! ss -ltn "sport = :$UPSTREAM_PORT" 2>/dev/null | grep -q LISTEN; then
        warn "AptaRank is not running on $UPSTREAM_PORT — start it first:"
        warn "  $APP_DIR/deploy/aptarank.sh start"
        exit 1
    fi

    nohup "$CADDY" run --config "$CADDYFILE" --adapter caddyfile \
        >>"$LOG_FILE" 2>&1 &
    echo $! > "$PID_FILE"
    sleep 3

    if kill -0 "$(cat "$PID_FILE")" 2>/dev/null; then
        address="https://$(hostname -I | awk '{print $1}'):$LISTEN_PORT"
        echo
        ok "serving on $address"
        echo
        echo "    logins       $(cut -d: -f1 "$USERS_FILE" | tr '\n' ' ')"
        echo "    upstream     127.0.0.1:$UPSTREAM_PORT (still loopback only)"
        echo "    log          $LOG_FILE"
        echo
        echo "  First visit shows a certificate warning: the certificate is issued"
        echo "  by Caddy's own authority rather than a public one. Colleagues can"
        echo "  click through it once, or you can install the authority's"
        echo "  certificate on their machines to remove the warning:"
        echo
        echo "      $DATA_DIR/caddy-storage/pki/authorities/local/root.crt"
        echo
        echo "  If your institution issues certificates, use one instead of"
        echo "  'tls internal' in $CADDYFILE."
        echo
    else
        warn "Caddy failed to start. Last lines of $LOG_FILE:"
        tail -20 "$LOG_FILE"
        exit 1
    fi
    ;;

stop)
    if [ -f "$PID_FILE" ] && kill -0 "$(cat "$PID_FILE")" 2>/dev/null; then
        kill "$(cat "$PID_FILE")"
        rm -f "$PID_FILE"
        echo "Stopped. AptaRank itself is still running on 127.0.0.1:$UPSTREAM_PORT"
        echo "and remains reachable through an SSH tunnel."
    else
        echo "Not serving."
    fi
    ;;

restart)
    "$0" stop || true
    sleep 1
    "$0" start
    ;;

status)
    if [ -f "$PID_FILE" ] && kill -0 "$(cat "$PID_FILE")" 2>/dev/null; then
        echo "serving   https://$(hostname -I | awk '{print $1}'):$LISTEN_PORT (pid $(cat "$PID_FILE"))"
        echo "logins    $(cut -d: -f1 "$USERS_FILE" 2>/dev/null | tr '\n' ' ')"
    else
        echo "not serving (AptaRank reachable only via SSH tunnel)"
    fi
    ;;

*)
    echo "usage: $0 {add <user>|remove <user>|users|start|stop|restart|status}"
    exit 2
    ;;
esac
