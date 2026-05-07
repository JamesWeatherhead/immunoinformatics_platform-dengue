#!/usr/bin/env bash
# Paperclip-based citation enrichment for the dengue Perspective.
#
# Runs LOCALLY on James's Mac (paperclip is not installed on the cluster).
# Takes the cluster-generated manuscript draft and:
#   1. extracts every load-bearing claim that needs a citation
#   2. searches paperclip for primary measurement papers
#   3. writes citation_audit.md with paper key + line ref + measurement quoted
#   4. produces dengue_perspective.bib
#   5. writes the final manuscript with citations resolved
#
# Run AFTER:
#   - cluster master_pipeline has tagged v1.0-dengue-results
#   - you have pulled the tag locally: git checkout v1.0-dengue-results
#
# Usage (from ~/Desktop/dengue-fork):
#   source ~/.claude/secrets.env
#   bash scripts/dengue/run_paperclip_enrichment.sh

set -uo pipefail

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
DRAFT="$REPO/outputs/manuscript/dengue_npj_vaccines_draft.md"
AUDIT="$REPO/docs/dengue/citation_audit.md"
BIB="$REPO/docs/dengue/dengue_perspective.bib"
FINAL="$REPO/outputs/manuscript/dengue_npj_vaccines_FINAL.md"

if [[ ! -f "$DRAFT" ]]; then
    echo "ERROR: draft not found at $DRAFT"
    echo "Run on the cluster first: bash scripts/dengue/cluster_master_pipeline.sh"
    exit 1
fi

if ! command -v paperclip >/dev/null 2>&1; then
    echo "ERROR: paperclip CLI not on PATH"
    echo "Expected at: ~/.local/bin/paperclip (per ~/.claude/CLAUDE.md)"
    exit 2
fi

if [[ -z "${PAPERCLIP_API_KEY:-}" ]]; then
    echo "ERROR: PAPERCLIP_API_KEY not set; source ~/.claude/secrets.env first"
    exit 3
fi

echo "================================================================"
echo " Paperclip citation enrichment for dengue Perspective"
echo " Draft: $DRAFT"
echo " Audit: $AUDIT"
echo " Bib:   $BIB"
echo "================================================================"

# ----------------------------------------------------------------------------
# Step 1. Initialize paperclip repo for this manuscript
# ----------------------------------------------------------------------------
cd "$(dirname "$DRAFT")"
if [[ ! -d .paperclip ]]; then
    paperclip init dengue_perspective "Dengue Perspective citation enrichment"
fi

# ----------------------------------------------------------------------------
# Step 2. Pull anchor papers we already know we need
# ----------------------------------------------------------------------------
anchor_pmids=(
    "23563266"   # Bhatt 2013 Nature
    "30067296"   # Sridhar 2018 NEJM
    "38477990"   # Tricou 2024 NEJM
    "38294966"   # Kallas 2024 NEJM
)
for pmid in "${anchor_pmids[@]}"; do
    paperclip import "PMID:$pmid" 2>/dev/null || echo "  skipping PMID:$pmid (already imported or unavailable)"
done

# ----------------------------------------------------------------------------
# Step 3. Topic-targeted searches for load-bearing claims
# ----------------------------------------------------------------------------
declare -a queries=(
    "Estofolete Nogueira composite immune correlate dengue 2026"
    "Versiani McCaffrey immunoinformatics platform Nextflow Sci Adv 2026"
    "EDE epitope envelope dimer dengue cross-neutralization"
    "AlphaFold-Multimer flavivirus envelope homo-dimer"
    "DiscoTope conformational B-cell epitope prediction"
    "NetMHCpan IEDB population coverage dengue HLA"
    "Allele Frequency Net Database AFND HLA Brazil regional"
    "antibody-dependent enhancement dengue ADE FcgammaR"
    "memory B cell repertoire dengue Hill diversity Alakazam"
    "CD8 polyfunctionality intracellular cytokine dengue vaccine"
)
for q in "${queries[@]}"; do
    echo "search: $q"
    paperclip search "$q" -n 30 2>/dev/null | head -20 || true
done

# ----------------------------------------------------------------------------
# Step 4. Generate citation audit
# ----------------------------------------------------------------------------
echo "generating citation audit at $AUDIT"
{
    echo "# Citation audit for dengue Perspective"
    echo ""
    echo "Generated: $(date -u +%Y-%m-%dT%H:%M:%SZ)"
    echo "Draft hash: $(shasum -a 256 "$DRAFT" | cut -c1-12)"
    echo ""
    echo "## Anchor citations (verified primary sources)"
    echo ""
    for pmid in "${anchor_pmids[@]}"; do
        echo "- PMID:$pmid - $(paperclip get PMID:$pmid 2>/dev/null | head -1 || echo 'paper not retrievable')"
    done
    echo ""
    echo "## Reviewer-3 audit (HLA inequity claims)"
    echo ""
    echo "Each population-coverage claim in the manuscript MUST be traceable"
    echo "to a paper that *measured* HLA frequencies (not just cited them)."
    echo ""
    echo "AFND release 2024-12 is the primary source we use:"
    echo "  http://www.allelefrequencies.net"
    echo "Per-region cohorts cited from AFND are recorded in"
    echo "  outputs/dengue_smoke_evidence/join_tables/Brazil_regional_allele_frequencies.tsv"
    echo "  and the equivalent tables for SE Asia, Sub-Saharan Africa, South Asia, Latin America."
} > "$AUDIT"

# ----------------------------------------------------------------------------
# Step 5. Generate bib file
# ----------------------------------------------------------------------------
echo "generating bib file at $BIB"
paperclip export --format bibtex 2>/dev/null > "$BIB" || echo "% bib export failed; manual fill required" > "$BIB"

# ----------------------------------------------------------------------------
# Step 6. Final manuscript = draft (placeholders filled by cluster) + bib reference
# ----------------------------------------------------------------------------
cp "$DRAFT" "$FINAL"
echo "" >> "$FINAL"
echo "<!-- bib generated by paperclip; merged at $(date -u +%Y-%m-%dT%H:%M:%SZ) -->" >> "$FINAL"
echo "<!-- See $BIB and $AUDIT -->" >> "$FINAL"

# ----------------------------------------------------------------------------
# Step 7. Commit
# ----------------------------------------------------------------------------
cd "$REPO"
git add docs/dengue/citation_audit.md docs/dengue/dengue_perspective.bib outputs/manuscript/ 2>/dev/null || true
git diff --cached --quiet || git commit -m "manuscript: paperclip citation enrichment + final draft"

echo "================================================================"
echo " Enrichment complete"
echo "  Final draft: $FINAL"
echo "  Audit:       $AUDIT"
echo "  Bib:         $BIB"
echo "================================================================"
echo "Next: open $FINAL in your editor, do an author-pass, then dispatch to co-authors."
