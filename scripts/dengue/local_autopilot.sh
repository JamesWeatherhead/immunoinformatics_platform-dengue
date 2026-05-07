#!/usr/bin/env bash
# Mac-side autopilot: detects v1.0-dengue-results tag landing, then
# triggers paperclip enrichment + headless Claude rewrite of the manuscript.
#
# Runs hourly via launchd (com.dengue.autopilot.plist). Safe to invoke
# manually for testing. Idempotent: marks done with .autopilot_done so
# repeat runs don't redo the rewrite.
#
# Usage (one-shot test):
#   bash scripts/dengue/local_autopilot.sh
#
# Usage (launchd-driven; auto):
#   launchctl load ~/Library/LaunchAgents/com.dengue.autopilot.plist

set -uo pipefail

REPO="$HOME/Desktop/dengue-fork"
LOG="$HOME/Library/Logs/dengue_autopilot.log"
DONE_MARKER="$REPO/outputs/manuscript/.autopilot_done"
mkdir -p "$(dirname "$LOG")"

log() { printf '[%s] %s\n' "$(date -u +%Y-%m-%dT%H:%M:%SZ)" "$*" | tee -a "$LOG"; }

# Source secrets for paperclip + claude API key
if [[ -f "$HOME/.claude/secrets.env" ]]; then
    source "$HOME/.claude/secrets.env"
fi

log "autopilot iteration starting"

if [[ ! -d "$REPO" ]]; then
    log "ERROR: repo $REPO not found; cannot proceed"
    exit 1
fi

cd "$REPO"

# ----------------------------------------------------------------------------
# Step 1: pull latest from GitHub (cluster pushes outputs every 30 min)
# ----------------------------------------------------------------------------
log "  pulling latest from origin/master"
git pull --rebase --autostash origin master 2>&1 | tail -3 | tee -a "$LOG"

# Sync remote tags
git fetch --tags origin 2>&1 | tail -3 | tee -a "$LOG"

# ----------------------------------------------------------------------------
# Step 2: check for v1.0-dengue-results tag
# ----------------------------------------------------------------------------
if ! git tag -l | grep -q "^v1.0-dengue-results$"; then
    log "  v1.0-dengue-results not yet tagged; checking cluster status via SSH"

    # Snapshot cluster state for posterity
    if [[ -x /tmp/ssh_cluster.exp ]] && [[ -n "${SSH_CLUSTER_PASS:-}" ]]; then
        export SSH_CLUSTER_PASS
        /tmp/ssh_cluster.exp "bash $REPO/scripts/dengue/morning_status.sh 2>&1 | head -80" \
            >> "$HOME/Library/Logs/dengue_morning_snapshots.log" 2>&1 || \
            log "  cluster SSH check failed (cluster offline?)"
    fi
    log "autopilot iteration done (waiting for tag)"
    exit 0
fi

log "  v1.0-dengue-results TAG DETECTED"

# ----------------------------------------------------------------------------
# Step 3: idempotency - skip if already enriched
# ----------------------------------------------------------------------------
if [[ -f "$DONE_MARKER" ]]; then
    log "  autopilot already completed (marker $DONE_MARKER exists); skipping"
    exit 0
fi

# Stay on master. Tag is just a "done" signal; outputs are already on master HEAD.
log "  staying on master branch (tag is a completion signal, outputs already at HEAD)"

# ----------------------------------------------------------------------------
# Step 4: paperclip enrichment (citation graph)
# ----------------------------------------------------------------------------
log "  running paperclip enrichment"
if command -v paperclip >/dev/null 2>&1; then
    bash "$REPO/scripts/dengue/run_paperclip_enrichment.sh" 2>&1 | tail -40 | tee -a "$LOG"
else
    log "  WARN: paperclip CLI not found; skipping enrichment, going straight to Claude rewrite"
fi

# ----------------------------------------------------------------------------
# Step 5: dispatch headless Claude (Opus 4.7) to do paperclip-grounded rewrite
# ----------------------------------------------------------------------------
log "  dispatching headless Claude Opus 4.7 for paperclip-grounded manuscript rewrite"

if ! command -v claude >/dev/null 2>&1; then
    log "  ERROR: claude CLI not found; cannot do automated rewrite"
    log "  manual fallback: open $REPO/outputs/manuscript/dengue_npj_vaccines_draft.md and rewrite with paperclip"
    exit 2
fi

# Three Claude invocations with progressively scoped tasks. Each writes to a
# specific output file so they don't collide. The first one is the heavy
# paperclip-grounding pass; subsequent ones polish.

DRAFT="$REPO/outputs/manuscript/dengue_npj_vaccines_draft.md"
ROUND1="$REPO/outputs/manuscript/dengue_round1_paperclip.md"
ROUND2="$REPO/outputs/manuscript/dengue_round2_lancet_concerns.md"
ROUND3="$REPO/outputs/manuscript/dengue_round3_polish.md"
FINAL="$REPO/outputs/manuscript/dengue_npj_vaccines_FINAL.md"

# Round 1: paperclip-grounded citation enrichment
log "    Round 1: paperclip-grounded citation enrichment"
claude --model opus -p \
    --permission-mode bypassPermissions \
    --max-budget-usd 8 \
    --add-dir "$REPO" \
    "Read the dengue Perspective draft at $DRAFT.

Use paperclip to verify and enrich every citation. Specifically:

1. The McCaffrey paper (Versiani et al. 2026 Sci Adv 12:eaeb2066, immunoinformatics_platform) — pull this into paperclip and confirm we cite it correctly for the pipeline pedigree in Section 2.1 and the limitations in Section 4.

2. The Estofolete et al. 2026 npj Vaccines paper (composite immune correlate, Table 4) — this is THE paper our pipeline implements. Pull it into paperclip and confirm Section 1 framing accurately represents Tier A and Tier B definitions.

3. Each Phase 3 vaccine reference (Sridhar 2018 NEJM, Tricou 2024 NEJM, Kallas 2024 NEJM) — verify exact efficacy point estimates against published numbers.

4. Bhatt 2013 Nature for global burden — verify 390M / 96M numbers.

5. Reviewer 3's accusation: 9 citations were cited for HLA inequity but none measured HLA. Use paperclip to FIND papers that DID measure HLA in dengue-endemic regions, and replace any vague citations.

Write the enriched manuscript to $ROUND1. Keep the structure of $DRAFT but substitute proper citations everywhere a claim was vague. Preserve all pipeline numbers from fill_manuscript.py exactly.

Use Agent subagents in parallel (general-purpose) to research each citation block independently. Each subagent should return paperclip doc IDs + line refs.

Spawn no more than 6 subagents at once. Total budget for this round: 30 min wall time.

Important: write your final result to $ROUND1 using the Write tool before exiting." \
    --output-format text \
    > "$REPO/outputs/manuscript/.claude_round1_log.txt" 2>&1 &
ROUND1_PID=$!
log "    Round 1 PID: $ROUND1_PID; will wait up to 60 min"

# Wait with timeout
wait_with_timeout() {
    local pid=$1
    local timeout=$2
    local elapsed=0
    while kill -0 "$pid" 2>/dev/null; do
        if (( elapsed > timeout )); then
            log "    TIMEOUT after ${timeout}s; killing PID $pid"
            kill "$pid" 2>/dev/null || true
            return 1
        fi
        sleep 60
        elapsed=$((elapsed + 60))
    done
    return 0
}
wait_with_timeout "$ROUND1_PID" 3600

# Round 2: address Lancet ID reviewer concerns
log "    Round 2: address Lancet ID reviewer concerns explicitly"
claude --model opus -p \
    --permission-mode bypassPermissions \
    --max-budget-usd 4 \
    --add-dir "$REPO" \
    "Read the manuscript at ${ROUND1} (or $DRAFT if Round 1 didn't write it).

The prior version of this work (THELANCETID-D-26-00582) was rejected at Lancet Infectious Diseases on three grounds:

R1: Misrepresentation of cited literature.
R2: Speculative writing register.
R3: 9 references for HLA inequity, none of which actually measured HLA.

Section 4.3 of the draft addresses each, but lazily. Rewrite Section 4.3 so each subsection:
- Names the specific paper(s) we previously cited that triggered the concern (use paperclip if needed)
- States exactly what we changed
- Points to the specific section/figure of THIS Perspective that fixes it

Save the result to $ROUND2 using the Write tool. Do NOT change Sections 1-3 or 4.1-4.2 unless directly tied to a reviewer concern.

Budget: 20 min." \
    --output-format text \
    > "$REPO/outputs/manuscript/.claude_round2_log.txt" 2>&1 &
ROUND2_PID=$!
wait_with_timeout "$ROUND2_PID" 2400

# Round 3: polish for npj Vaccines register
log "    Round 3: language polish for npj Vaccines register"
claude --model opus -p \
    --permission-mode bypassPermissions \
    --max-budget-usd 3 \
    --add-dir "$REPO" \
    "Read the manuscript at ${ROUND2} (or fall back to $ROUND1 or $DRAFT).

Polish the language for npj Vaccines Perspective register:
- Active voice where possible
- No em-dashes (use commas, semicolons, colons)
- No overclaiming (we have r=~0.X across n=3 vaccine programs; don't say it 'predicts' efficacy, say it 'rank-orders')
- Preserve every pipeline number EXACTLY
- Keep word count under 5000

Save the polished result to $FINAL using the Write tool.

Budget: 15 min." \
    --output-format text \
    > "$REPO/outputs/manuscript/.claude_round3_log.txt" 2>&1 &
ROUND3_PID=$!
wait_with_timeout "$ROUND3_PID" 1800

# ----------------------------------------------------------------------------
# Step 6: commit + push final manuscript
# ----------------------------------------------------------------------------
cd "$REPO"
git add outputs/manuscript/ docs/dengue/ 2>>"$LOG" || true
if git diff --cached --quiet; then
    log "  no manuscript changes to commit (Claude rounds may have failed)"
else
    git commit -m "manuscript: paperclip-grounded rewrite via 3-round headless Claude Opus 4.7" 2>&1 | tail -3 | tee -a "$LOG"
    git push origin master 2>&1 | tail -3 | tee -a "$LOG"
fi

touch "$DONE_MARKER"
log "autopilot complete; manuscript at $FINAL"
log "open the manuscript: open '$FINAL'"
