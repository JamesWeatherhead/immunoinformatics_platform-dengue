#!/usr/bin/env bash
# Fetch dengue reference and Phase 3 vaccine parent strain protein sequences
# for the dengue extension of the immunoinformatics_platform pipeline.
#
# Outputs (relative to repo root):
#   ref_fastas/dengue/DENV-{1,2,3,4}_polyprotein_*.fasta            - UniProt SwissProt references
#   ref_fastas/dengue/dengue_polyproteins_multiseq.fasta            - concatenated input for Nextflow
#   ref_fastas/dengue_phase3_vaccines/{CYD-TDV,TAK-003,Butantan-DV}_DENV{1-4}_*.fasta
#   ref_fastas/dengue_phase3_vaccines/dengue_phase3_vaccines_multiseq.fasta
#   ref_fastas/dengue_phase3_vaccines/PROVENANCE.md
#
# Usage:
#   bash scripts/dengue/fetch_dengue_data.sh

set -euo pipefail
cd "$(dirname "$0")/../.."

mkdir -p ref_fastas/dengue ref_fastas/dengue_phase3_vaccines

# DENV 1-4 reference polyproteins from UniProt SwissProt
echo "[1/3] DENV 1-4 reference polyproteins from UniProt"
declare -a REFS=(
  "DENV-1 P17763"
  "DENV-2 P29990"
  "DENV-3 P27915"
  "DENV-4 P09866"
)
for line in "${REFS[@]}"; do
  read -r sero acc <<< "$line"
  outfile="ref_fastas/dengue/${sero}_polyprotein_${acc}.fasta"
  if [[ -s "$outfile" ]]; then
    echo "  exists: $outfile"
    continue
  fi
  echo "  fetch: $sero ($acc)"
  curl -fsSL "https://rest.uniprot.org/uniprotkb/${acc}.fasta" > "$outfile"
done

# Combined multi-FASTA for pipeline input
cat ref_fastas/dengue/DENV-*.fasta > ref_fastas/dengue/dengue_polyproteins_multiseq.fasta
echo "  combined: ref_fastas/dengue/dengue_polyproteins_multiseq.fasta ($(grep -c '^>' ref_fastas/dengue/dengue_polyproteins_multiseq.fasta) sequences)"

# Phase 3 vaccine parent strain GenBank entries
echo "[2/3] Phase 3 vaccine parent strain protein sequences from GenBank"
EFETCH="https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"

# Format: VACCINE_NAME SEROTYPE GENBANK_ACCESSION PARENT_STRAIN_LABEL
declare -a CONSTRUCTS=(
  # CYD-TDV (Sanofi Dengvaxia) parent strains; YF17D chimeric backbone
  "CYD-TDV     1 M87512   PUO-359-82_DENV-1"
  "CYD-TDV     2 M29095   PUO-218_DENV-2"
  "CYD-TDV     3 M93130   PaH881-88_DENV-3"
  "CYD-TDV     4 M14931   1228_DENV-4"
  # TAK-003 (Takeda Qdenga) parent strains; PDK-53 attenuated DENV-2 backbone
  "TAK-003     1 EU081230 16007_PDK-13_DENV-1"
  "TAK-003     2 M20558   16681_PDK-53_DENV-2"
  "TAK-003     3 EU081240 16562_PGMK-30_DENV-3"
  "TAK-003     4 EU081245 1036_PDK-48_DENV-4"
  # Butantan-DV / TV003 (NIAID/Butantan) rDEN attenuated tetravalent
  "Butantan-DV 1 EU848545 WestPac-74_rDEN1d30_DENV-1"
  "Butantan-DV 2 GU289914 NewGuineaC_rDEN2-4d30_DENV-2"
  "Butantan-DV 3 EF629369 Sleman-78_rDEN3d30-31_DENV-3"
  "Butantan-DV 4 AF326825 Dominica-814669_rDEN4d30_DENV-4"
)

for line in "${CONSTRUCTS[@]}"; do
  read -r vaccine sero acc strain <<< "$line"
  outfile="ref_fastas/dengue_phase3_vaccines/${vaccine}_DENV${sero}_${acc}.fasta"
  if [[ -s "$outfile" ]]; then
    echo "  exists: $outfile"
    continue
  fi
  echo "  fetch: $vaccine DENV-$sero ($acc, $strain)"
  # Try fasta_cds_aa first (preferred: gives translated proteins).
  # If that returns empty (some older GenBank records lack CDS_AA), fall
  # back to fetching the nucleotide and translating the longest ORF.
  curl -fsS --retry 3 --retry-delay 2 \
      "${EFETCH}?db=nuccore&id=${acc}&rettype=fasta_cds_aa&retmode=text" \
      -o "$outfile.tmp" 2>/dev/null || true
  if [[ ! -s "$outfile.tmp" ]] || [[ $(wc -c < "$outfile.tmp") -lt 50 ]]; then
    echo "    fasta_cds_aa empty for $acc, falling back to nucleotide + translate"
    rm -f "$outfile.tmp"
    nuc_tmp="$outfile.nuc.tmp"
    if curl -fsS --retry 3 --retry-delay 2 \
        "${EFETCH}?db=nuccore&id=${acc}&rettype=fasta&retmode=text" \
        -o "$nuc_tmp"; then
      python3 - <<PYNUC
import sys
from pathlib import Path
nuc_path = Path("$nuc_tmp")
out_path = Path("$outfile.tmp")
text = nuc_path.read_text().strip()
if not text or not text.startswith(">"):
    sys.exit(0)
header, *seq_lines = text.splitlines()
seq = "".join(seq_lines).upper().replace("U", "T")
# Translate all 3 forward frames; take longest ORF (M..*) per frame
def translate(s):
    table = {
        'TTT':'F','TTC':'F','TTA':'L','TTG':'L','CTT':'L','CTC':'L','CTA':'L','CTG':'L',
        'ATT':'I','ATC':'I','ATA':'I','ATG':'M','GTT':'V','GTC':'V','GTA':'V','GTG':'V',
        'TCT':'S','TCC':'S','TCA':'S','TCG':'S','CCT':'P','CCC':'P','CCA':'P','CCG':'P',
        'ACT':'T','ACC':'T','ACA':'T','ACG':'T','GCT':'A','GCC':'A','GCA':'A','GCG':'A',
        'TAT':'Y','TAC':'Y','TAA':'*','TAG':'*','CAT':'H','CAC':'H','CAA':'Q','CAG':'Q',
        'AAT':'N','AAC':'N','AAA':'K','AAG':'K','GAT':'D','GAC':'D','GAA':'E','GAG':'E',
        'TGT':'C','TGC':'C','TGA':'*','TGG':'W','CGT':'R','CGC':'R','CGA':'R','CGG':'R',
        'AGT':'S','AGC':'S','AGA':'R','AGG':'R','GGT':'G','GGC':'G','GGA':'G','GGG':'G',
    }
    return "".join(table.get(s[i:i+3], "X") for i in range(0, len(s)-2, 3))
best = ""
for fr in range(3):
    p = translate(seq[fr:])
    # find longest M..stop subsequence; flaviviruses have one big polyprotein
    cur = ""
    for ch in p:
        if cur:
            if ch == "*":
                if len(cur) > len(best):
                    best = cur
                cur = ""
            else:
                cur += ch
        elif ch == "M":
            cur = "M"
    if cur and len(cur) > len(best):
        best = cur
if best:
    out_path.write_text(f">${vaccine}_DENV${sero}_${acc}_translated_longest_ORF\n{best}\n")
    print(f"    translated nucleotide -> {len(best)} aa")
else:
    print("    translation produced no ORF")
PYNUC
      rm -f "$nuc_tmp"
    fi
  fi
  if [[ -s "$outfile.tmp" ]]; then
    # Rewrite header to make construct provenance explicit, keep all CDS records
    python3 - <<PY
from pathlib import Path
src = Path("$outfile.tmp").read_text()
records = src.split('\n>')
out = []
for i, rec in enumerate(records):
    if not rec.strip():
        continue
    rec = rec if rec.startswith('>') or i == 0 else '>' + rec
    if not rec.startswith('>'):
        rec = '>' + rec
    head, _, body = rec.partition('\n')
    head = f">${vaccine}_DENV${sero}_${acc}_CDS{i+1} | parent_strain=${strain} | upstream_header={head[1:].strip()[:120]}"
    out.append(head + '\n' + body.rstrip())
Path("$outfile").write_text('\n'.join(out) + '\n')
PY
    rm -f "$outfile.tmp"
  else
    echo "    WARN: still empty after fallback for $acc; will need manual resolution"
    rm -f "$outfile.tmp"
  fi
  sleep 1
done

# Combined Phase 3 multi-FASTA for pipeline input (only includes successful fetches)
shopt -s nullglob
phase3_fastas=(ref_fastas/dengue_phase3_vaccines/*_DENV*.fasta)
if (( ${#phase3_fastas[@]} > 0 )); then
  cat "${phase3_fastas[@]}" > ref_fastas/dengue_phase3_vaccines/dengue_phase3_vaccines_multiseq.fasta
  echo "  combined: ref_fastas/dengue_phase3_vaccines/dengue_phase3_vaccines_multiseq.fasta ($(grep -c '^>' ref_fastas/dengue_phase3_vaccines/dengue_phase3_vaccines_multiseq.fasta) sequences)"
fi

# Provenance documentation
echo "[3/3] writing provenance documentation"
cat > ref_fastas/dengue_phase3_vaccines/PROVENANCE.md <<'EOF'
# Phase 3 dengue vaccine parent strain provenance

This directory holds GenBank protein records for the parent strains of the three
licensed Phase 3 dengue vaccine constructs. They are used by the retrospective
validation script (`scripts/dengue/retrospective_phase3_validation.py`) to
score CYD-TDV, TAK-003, and Butantan-DV through the pipeline and compare to
published clinical efficacy.

## Important: parent strain != licensed construct

The licensed vaccines carry attenuation mutations (3' UTR deletions in the
Butantan-DV components, PDK passage adaptations in TAK-003, YF17D chimerism
in CYD-TDV). For the Geometry and Equity axes that depend on the E protein
primary sequence, these mutations are not expected to substantively change
the predicted scores. For a stricter analysis, replace these FASTAs with
construct-specific sequences from each program's FDA Biologics License
Application supplement.

## Construct summary

| Vaccine | Sponsor | Backbone strategy | Parent strains (DENV-1/2/3/4) | Pivotal trial |
|---|---|---|---|---|
| CYD-TDV (Dengvaxia) | Sanofi Pasteur | YF17D chimeric | PUO-359/82, PUO-218, PaH881/88, 1228 | Capeding 2014, Hadinegoro 2015, Sridhar 2018 |
| TAK-003 (Qdenga) | Takeda | DENV-2 PDK-53 | 16007 PDK-13, 16681 PDK-53, 16562 PGMK-30, 1036 PDK-48 | Biswal 2019, Tricou 2024 |
| Butantan-DV / TV003 | NIAID / Butantan | rDEN attenuated tetravalent | WestPac-74, NewGuineaC (chimeric DENV-2/4), Sleman-78, Dominica-814669 | Kirkpatrick 2016, Kallas 2024, Nogueira 2024 |

## Citations

- Capeding MR et al. Lancet 2014;384:1358 (CYD-TDV Phase 3 Asia)
- Hadinegoro SR et al. NEJM 2015;373:1195 (CYD-TDV long-term)
- Sridhar S et al. NEJM 2018;379:327 (CYD-TDV serostatus stratified)
- Huang CY-H et al. J Virol 2013;87:7036 (TAK-003 backbone)
- Biswal S et al. NEJM 2019;381:2009 (TAK-003 Phase 3 1y)
- Tricou V et al. Lancet Glob Health 2024;12:e257 (TAK-003 4.5y)
- Whitehead SS et al. PLoS NTD 2017;11:e0005584 (TV003 Phase 1)
- Kirkpatrick BD et al. Sci Transl Med 2016;8:330ra36 (TV003 challenge)
- Kallas EG et al. NEJM 2024;390:397 (Butantan-DV Phase 3 single-dose)
- Nogueira ML et al. Lancet Infect Dis 2024;24 (Butantan-DV extended follow-up)
- Estofolete CF, Saivish MV, Nogueira ML, Vasilakis N. npj Vaccines 2026;11:68
EOF

echo
echo "DONE."
echo
echo "Inputs ready at:"
echo "  ref_fastas/dengue/dengue_polyproteins_multiseq.fasta            (4 sequences)"
echo "  ref_fastas/dengue_phase3_vaccines/dengue_phase3_vaccines_multiseq.fasta"
