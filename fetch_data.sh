#!/bin/bash
# Download the preprocessed FOWD-CDIP rogue wave inputs (Häfner et al. 2023).
# File list and base URL come from the paper's own code/get_rogue_data.sh. Resumable.
cd "$(dirname "$0")"
mkdir -p data
BASE_URL=$(grep '^BASE_URL=' code/get_rogue_data.sh | cut -d'"' -f2)
grep -oE 'fowd_cdip_[0-9a-z]+_filtered-preprocessed-agg\.parquet' code/get_rogue_data.sh | sort -u > data/filelist.txt
xargs -a data/filelist.txt -P 4 -I{} curl -sS -f -C - --retry 5 --retry-delay 5 -o data/{} "$BASE_URL/{}"
echo "exit=$? files=$(ls data/*.parquet 2>/dev/null | wc -l) size=$(du -sh data | cut -f1)"
