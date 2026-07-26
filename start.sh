#!/usr/bin/env bash
# AptaRank — run ./start.sh to start the dashboard.
#
# Deliberately avoids a bare `python`: the launcher resolves a suitable
# interpreter itself and explains what to install if it cannot find one.
set -euo pipefail
cd "$(dirname "$0")"

if [ -x ".venv/bin/python" ]; then
    exec .venv/bin/python scripts/start.py "$@"
fi

for candidate in python3.12 python3.11 python3.10 python3; do
    if command -v "$candidate" >/dev/null 2>&1; then
        exec "$candidate" scripts/start.py "$@"
    fi
done

cat <<'EOF'

  AptaRank needs Python 3.10 or newer, and could not find it.

    macOS:  brew install python@3.11
    Ubuntu: sudo apt install python3.11 python3.11-venv build-essential

  Then run ./start.sh again.

EOF
exit 2
