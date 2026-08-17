#!/usr/bin/env bash
# ProbeK updater. Pulls the latest code from GitHub and reinstalls
# dependencies. Never touches input_csvs/, reference_data/, results/, or
# .venv/ -- those aren't tracked by git, so a pull leaves them alone.
set -euo pipefail

cd "$(dirname "${BASH_SOURCE[0]}")"

if [ ! -d .git ]; then
    echo "ERROR: This folder isn't a git checkout, so it can't be updated automatically." >&2
    echo "This usually means ProbeK was installed via GitHub's 'Download ZIP' button" >&2
    echo "instead of 'git clone'." >&2
    echo >&2
    echo "To switch to an updatable install:" >&2
    echo "  1. Install git if needed (e.g. 'sudo apt install git')." >&2
    echo "  2. Clone a fresh copy elsewhere:" >&2
    echo "       git clone https://github.com/Sirawyn7/ProbeK.git ProbeK-new" >&2
    echo "  3. Copy your data into it (skips re-downloading 1GB+ of reference data):" >&2
    echo "       cp -r input_csvs/*.csv ProbeK-new/input_csvs/" >&2
    echo "       cp -r reference_data ProbeK-new/" >&2
    echo "  4. Use ProbeK-new/ from now on -- 'bash update.sh' will work there." >&2
    exit 1
fi

if ! command -v git >/dev/null 2>&1; then
    echo "ERROR: git was not found on PATH." >&2
    exit 1
fi

if [ -n "$(git status --porcelain --untracked-files=no)" ]; then
    echo "ERROR: This folder has local changes to tracked files, so it's not safe to auto-update." >&2
    git status --short --untracked-files=no >&2
    echo "Resolve or discard these changes first (see 'git status'), then try again." >&2
    exit 1
fi

echo "Checking for updates..."
git fetch --quiet origin

branch=$(git rev-parse --abbrev-ref HEAD)
upstream="origin/$branch"
if ! git rev-parse --verify "$upstream" >/dev/null 2>&1; then
    echo "ERROR: No '$upstream' branch found on the remote -- can't check for updates." >&2
    exit 1
fi

behind=$(git rev-list --count "HEAD..$upstream")
if [ "$behind" -eq 0 ]; then
    echo "ProbeK is already up to date."
    exit 0
fi

echo
echo "A new version of ProbeK is available."
read -r -p "Update now? [y/N]: " reply
case "$reply" in
    [yY]|[yY][eE][sS]) ;;
    *) echo "Update skipped."; exit 0 ;;
esac

git pull --ff-only origin "$branch"

VENV_DIR=".venv"
if [ -x "$VENV_DIR/bin/python" ]; then
    echo "Updating installed dependencies..."
    "$VENV_DIR/bin/python" -m pip install -e .
else
    echo "(No local Python environment yet -- run 'bash run.sh' next to set one up.)"
fi

echo
echo "Update complete."
