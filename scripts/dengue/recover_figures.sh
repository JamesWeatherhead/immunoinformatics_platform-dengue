#!/usr/bin/env bash
# Recovery: ensure all 6 figures exist. The generate_figures.py script
# has placeholder fallback for missing data, so this is mostly about
# making sure pip + matplotlib are available.

set -uo pipefail
REPO=/data/james/ica-dengue/immunoinformatics_platform-dengue
LOG=/data/james/ica-dengue/logs/recover_figures.log
FIG_OUT=/data/james/ica-dengue/outputs/figures

log() { printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$LOG"; }

log "figures recovery starting"
cd "$REPO"
git pull --rebase --autostash origin master 2>&1 | tail -3 | tee -a "$LOG"

mkdir -p "$FIG_OUT"

# Ensure deps
python3 -c "import matplotlib, numpy" 2>/dev/null || pip3 install --user matplotlib numpy 2>&1 | tail -3 | tee -a "$LOG"

python3 "$REPO/scripts/dengue/generate_figures.py" \
    --table4_dir /data/james/ica-dengue/outputs/table4_dengue \
    --phase3_dir /data/james/ica-dengue/outputs/phase3_retrospective \
    --bcell_dir /data/james/ica-dengue/outputs/dengue_bcell \
    --tcell_dir /data/james/ica-dengue/outputs/dengue_tcell \
    --out_dir "$FIG_OUT" \
    >> "$LOG" 2>&1 || log "  generate_figures FAILED; placeholder PNGs will be missing"

log "figures recovery done; rendered $(ls $FIG_OUT/*.png 2>/dev/null | wc -l) figures"
