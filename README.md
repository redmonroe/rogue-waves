# rogue-waves

Testing data-driven rogue wave probability models against public buoy data, with a Lean 4 + Mathlib project for any formal work.

## Layout

- `RogueWaves/` — Lean 4 project (Mathlib, toolchain pinned in `lean-toolchain`)
- `fetch_data.sh` — downloads the preprocessed FOWD-CDIP dataset into `data/` (157 parquet files, 6.2 GB, 14,985,244 rows, 106 columns)
- `setup.sh` — one-time setup on a new machine: paper code, Lean, Mathlib cache, Python env, data
- `research_notes/`, `reports/` — background research on whether this problem is worth pursuing
- `code/` (not tracked) — clone of github.com/dionhaefner/rogue-wave-discovery
- `data/` (not tracked)

## Setup

    bash setup.sh

## Sources

- Häfner, Gemmrich, Jochum. "Machine-guided discovery of a real-world rogue wave model." PNAS 2023. doi:10.1073/pnas.2306275120
- Data: https://erda.ku.dk/archives/ee6b452c1907fbd48271b071c3cee10e/published-archive.html

## Notes

- The paper's full training run needs 64 GB RAM and about 2 days on 16 cores. Scoring the published formula against the data needs far less.
- `lakefile.toml` pulls Mathlib by git URL rather than through the Lean package index.
