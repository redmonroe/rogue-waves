#!/bin/bash
# One-time setup on a new machine (macOS/Linux). Run from the repo root.
set -e
cd "$(dirname "$0")"

# 1. Rogue wave paper code (Häfner et al. 2023, MIT license); only needed with --data
[ -d code ] || git clone --depth 1 https://github.com/dionhaefner/rogue-wave-discovery.git code

# 2. Lean toolchain + prebuilt Mathlib
command -v elan >/dev/null || curl -sSf https://raw.githubusercontent.com/leanprover/elan/master/elan-init.sh | sh -s -- -y --default-toolchain none
export PATH="$HOME/.elan/bin:$PATH"
(cd RogueWaves && lake exe cache get && lake build)

# 3. Python env for reading the data
[ -d .venv ] || python3 -m venv .venv
.venv/bin/pip install -q pyarrow sympy mpmath numpy scipy

# 4. Rogue wave data (optional): 157 parquet files, 6.2 GB, from erda.ku.dk (resumable)
if [ "$1" = "--data" ]; then bash fetch_data.sh; fi
