#!/usr/bin/env bash
# Recovery: aggressive disk cleanup when /data exceeds 92% used.
# Targets: docker layers, old work_*, old tmp files, old logs.

set -uo pipefail
LOG=/data/james/ica-dengue/logs/recover_disk.log
log() { printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$LOG"; }

log "disk recovery starting; df: $(df -h /data | tail -1)"

# 1. Docker prune
docker system prune -af --volumes 2>&1 | tail -3 | tee -a "$LOG"

# 2. Old nextflow work dirs (keep dengue current ones)
find /data/james/ica-dengue -maxdepth 4 -type d -name "work_*" -mtime +1 \
    -exec rm -rf {} + 2>/dev/null

# 3. Old /tmp
find /tmp -maxdepth 2 -mtime +3 -delete 2>/dev/null

# 4. Old AF MSA intermediate files (jackhmmer leaves big tmp files)
# (only purge older than 1 day so we don't kill running AF)
find /data/james/ica-dengue/outputs/alphafold_e_dimers -name "msas" -type d -mtime +1 \
    -exec rm -rf {} + 2>/dev/null

log "disk recovery done; df: $(df -h /data | tail -1)"
