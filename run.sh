#!/usr/bin/env bash
# ProbeK launcher. On first run, prompts to set up a local Python environment
# (creates .venv/ in this folder and installs ProbeK + its dependencies --
# nothing system-wide). Every run after that skips straight to launching
# ProbeK, which handles its own first-run BLAST+/reference-data setup.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

VENV_DIR=".venv"
PROBEK_BIN="$VENV_DIR/bin/probek"

is_non_interactive=false
for arg in "$@"; do
    if [ "$arg" = "--non-interactive" ]; then
        is_non_interactive=true
    fi
done

if [ ! -x "$PROBEK_BIN" ]; then
    if ! command -v python3 >/dev/null 2>&1; then
        echo "ERROR: python3 was not found on PATH." >&2
        echo "Install Python 3.10+ (e.g. 'sudo apt install python3') and re-run ./run.sh." >&2
        exit 1
    fi

    if ! python3 -c 'import sys; sys.exit(0 if sys.version_info[:2] >= (3, 10) else 1)'; then
        found_version=$(python3 -c 'import sys; print("%d.%d" % sys.version_info[:2])')
        echo "ERROR: ProbeK needs Python 3.10+; found $found_version." >&2
        exit 1
    fi

    if $is_non_interactive; then
        echo "ERROR: Python environment isn't set up yet, and --non-interactive was given." >&2
        echo "Run ./run.sh once interactively first to complete first-time setup." >&2
        exit 1
    fi

    echo
    echo "ProbeK needs to set up a local Python environment before it can run."
    echo "This creates '$VENV_DIR/' in this folder and installs a few Python"
    echo "packages into it (pandas, numpy, requests, pyarrow) -- nothing is"
    echo "installed system-wide."
    read -r -p "Proceed? [y/N]: " reply
    case "$reply" in
        [yY]|[yY][eE][sS]) ;;
        *) echo "Setup declined. Nothing was installed."; exit 1 ;;
    esac

    if [ ! -d "$VENV_DIR" ]; then
        echo "Creating virtual environment..."
        venv_err=$(python3 -m venv "$VENV_DIR" 2>&1) || {
            if echo "$venv_err" | grep -qi "ensurepip"; then
                echo "python3's built-in pip bootstrap isn't available on this system;" \
                     "falling back to a manual pip install..."
                rm -rf "$VENV_DIR"
                python3 -m venv "$VENV_DIR" --without-pip
                curl -sS https://bootstrap.pypa.io/get-pip.py | "$VENV_DIR/bin/python"
            else
                echo "$venv_err" >&2
                exit 1
            fi
        }
    fi

    echo "Installing ProbeK and its dependencies..."
    "$VENV_DIR/bin/python" -m pip install --upgrade pip -q
    "$VENV_DIR/bin/python" -m pip install -e . -q
    echo "Setup complete."
    echo
fi

exec "$PROBEK_BIN" "$@"
