#!/usr/bin/env bash
# Recovery: ensure dengue_npj_vaccines_draft.md exists by running fill_manuscript
# against whatever data is present. fill_manuscript.py has TBD fallback for
# every placeholder, so this always produces a draft.

set -uo pipefail
REPO=/data/james/ica-dengue/immunoinformatics_platform-dengue
LOG=/data/james/ica-dengue/logs/recover_manuscript.log
MS_OUT=/data/james/ica-dengue/outputs/manuscript

log() { printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$LOG"; }

log "manuscript recovery starting"
cd "$REPO"
git pull --rebase --autostash origin master 2>&1 | tail -3 | tee -a "$LOG"
mkdir -p "$MS_OUT"

python3 "$REPO/scripts/dengue/fill_manuscript.py" \
    --template "$REPO/docs/dengue/manuscript_template.md" \
    --table4_dir /data/james/ica-dengue/outputs/table4_dengue \
    --phase3_dir /data/james/ica-dengue/outputs/phase3_retrospective \
    --figures_dir /data/james/ica-dengue/outputs/figures \
    --out "$MS_OUT/dengue_npj_vaccines_draft.md" \
    >> "$LOG" 2>&1 || log "  fill_manuscript FAILED"

if [[ -f "$MS_OUT/dengue_npj_vaccines_draft.md" ]]; then
    log "  draft exists at $MS_OUT/dengue_npj_vaccines_draft.md"
    log "  TBD count: $(grep -c TBD $MS_OUT/dengue_npj_vaccines_draft.md 2>/dev/null || echo ?)"
else
    log "  draft NOT created; copying template raw as fallback"
    cp "$REPO/docs/dengue/manuscript_template.md" "$MS_OUT/dengue_npj_vaccines_draft.md"
fi

log "manuscript recovery done"
