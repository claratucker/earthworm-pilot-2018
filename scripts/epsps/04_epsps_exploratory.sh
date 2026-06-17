#!/usr/bin/env bash
# 04_epsps_exploratory.sh - EPSPS sensitivity-class layer.
#
# Method: PICRUSt2 predicts which taxa carry the EPSPS gene (aroA, K00800)
# from 16S phylogenetic placement. It does not observe EPSPS directly.
# EPSPSClass then classifies reference EPSPS proteins for those taxa into
# sensitivity classes I-IV. Both steps are predictions, not measurements,
# so results should be reported as hypotheses, not findings. Measuring
# EPSPS class directly requires shotgun metagenomics or targeted EPSPS
# amplicons.
#
# Citations: Leino et al. 2021 Environ Int 149:106334 (classification
# framework), epspsclass repo for this implementation (note the 40% identity
# gate), Douglas et al. 2020 Nat Biotechnol 38:685 (PICRUSt2).
set -euo pipefail

PROJECT="$HOME/pilot2018"
BUCKET="s3://earthworm-pilot-2018"
PIC="$PROJECT/results/picrust2_retry"
EOUT="$PROJECT/results/epsps"; mkdir -p "$EOUT"
REFDIR="$PROJECT/refs/epsps"; mkdir -p "$REFDIR"
cd "$PROJECT"

if ! command -v epspsclass >/dev/null 2>&1; then
  python3 -m venv "$PROJECT/.epsps-venv"
  source "$PROJECT/.epsps-venv/bin/activate"
  pip install --quiet --upgrade pip
  pip install --quiet "git+https://github.com/claratucker/epspsclass.git"
else
  echo "epspsclass already on PATH"
fi
epspsclass validate-markers | head -20 || true

KO_PRED=$(find "$PIC" -maxdepth 1 -name 'combined_KO_predicted.tsv*' | head -1 || true)
if [ -z "$KO_PRED" ]; then
  KO_PRED=$(find "$PIC" -maxdepth 1 -name 'KO_predicted.tsv*' | head -1 || true)
fi
BAC_KO=$(find "$PIC" -maxdepth 1 -name 'bac_KO_predicted.tsv*' | head -1 || true)
ARC_KO=$(find "$PIC" -maxdepth 1 -name 'arc_KO_predicted.tsv*' | head -1 || true)

if [ -z "$KO_PRED" ]; then
  echo "No combined_KO_predicted.tsv* or KO_predicted.tsv* found under $PIC."
  echo "Checking for separate bac_KO_predicted.tsv.gz / arc_KO_predicted.tsv.gz."
fi

if [ -z "$KO_PRED" ] && { [ -n "$BAC_KO" ] || [ -n "$ARC_KO" ]; }; then
  echo "Merging bac/arc KO tables: $BAC_KO + $ARC_KO"
  python3 - "$BAC_KO" "$ARC_KO" "$EOUT/asv_aroA_predicted.tsv" <<'PY'
import sys, gzip, csv
bac_path, arc_path, out = sys.argv[1], sys.argv[2], sys.argv[3]

def read_k00800(path):
    rows = {}
    if not path:
        return rows
    op = gzip.open if path.endswith('.gz') else open
    with op(path, 'rt') as fh:
        r = csv.reader(fh, delimiter='\t'); header = next(r)
        candidates = [c for c in header if c == 'K00800' or c == 'ko:K00800']
        if not candidates:
            return rows
        j = header.index(candidates[0])
        for row in r:
            rows[row[0]] = row[j]
    return rows

bac = read_k00800(bac_path)
arc = read_k00800(arc_path)
all_asvs = set(bac) | set(arc)
with open(out, 'w') as o:
    o.write("asv\tK00800_predicted_copies\n")
    for asv in sorted(all_asvs):
        val = bac.get(asv, arc.get(asv, "0"))
        o.write(f"{asv}\t{val}\n")
print("Wrote merged per-ASV aroA prediction:", out, "(", len(all_asvs), "ASVs )")
PY
elif [ -n "$KO_PRED" ]; then
  echo "Found per-ASV KO predictions: $KO_PRED"
  python3 - "$KO_PRED" "$EOUT/asv_aroA_predicted.tsv" <<'PY'
import sys, gzip, csv
src, out = sys.argv[1], sys.argv[2]
op = gzip.open if src.endswith('.gz') else open
with op(src, 'rt') as fh:
    r = csv.reader(fh, delimiter='\t'); header = next(r)
    candidates = [c for c in header if c == 'K00800' or c == 'ko:K00800']
    if not candidates:
        print("K00800 not in KO table; no EPSPS predicted at all."); sys.exit(0)
    j = header.index(candidates[0])
    asv_col = 0
    with open(out, 'w') as o:
        o.write("asv\tK00800_predicted_copies\n")
        for row in r:
            o.write(f"{row[asv_col]}\t{row[j]}\n")
print("Wrote per-ASV aroA prediction:", out)
PY
else
  echo "No KO predictions found under $PIC. Check the picrust2 output."
fi

REF_FAA="$REFDIR/reference_epsps.faa"
if [ ! -s "$REF_FAA" ]; then
  cat <<'NOTE'
Build the reference EPSPS protein FASTA before running classification.
  File: refs/epsps/reference_epsps.faa
  One EPSPS (aroA) protein per genus flagged in step 1.
  Header: >Genus|accession   example: >Pseudomonas|WP_003114564.1
  Source: UniProt (query "EPSP synthase <genus>") or RefSeq protein (aroA).
Run scripts/r/05_epsps_overlay.R first for the exact genus list.
NOTE
  echo "No reference FASTA yet. Skipping classification."
  exit 0
fi

epspsclass classify \
  --input "$REF_FAA" \
  --output "$EOUT/epsps_classified.tsv" \
  --threshold 40 --summary 2> "$EOUT/epsps_summary.txt"

echo "EPSPSClass summary (reference proteins for predicted genera):"
cat "$EOUT/epsps_summary.txt"
echo "Per-sequence classes: $EOUT/epsps_classified.tsv"

aws s3 sync "$PROJECT/results/epsps/" "$BUCKET/results/epsps/"
echo "Done. Run scripts/r/05_epsps_overlay.R to join classes to genera and treatment."
