#!/usr/bin/env bash
# 01_diagnose.sh -- INSPECT the raw data before committing any processing
# parameters. This is the most important step for a legacy dataset: it tells us
# whether the reads are single-end, already-merged, or multiplexed, what the
# read-length distribution is, and whether primers/barcodes are still inline.
#
# WHY THIS MATTERS: the 2026 pipeline hardcodes a PAIRED-end DADA2 step with V4
# truncation tuned for ~253 bp paired reads. This 2018 run is single-end ~272 bp
# (per the MR DNA count file), with the original Caporaso 515F primer
# (GTGCCAGCMGCCGCGGTAA), not the 515F-Parada (GTGYCAGC...) used in 2026. If we
# applied the 2026 parameters blindly, denoising would fail silently. So we
# measure first, then set parameters in 02.
set -euo pipefail

PROJECT="$HOME/pilot2018"
RAW="$PROJECT/data/raw"
DIAG="$PROJECT/results/diagnostics"
mkdir -p "$DIAG"
cd "$PROJECT"

echo "=== Files present in raw/ ==="
ls -lhR "$RAW"

# Identify candidate read files (handles .fastq and .fastq.gz)
mapfile -t FASTQS < <(find "$RAW" -type f \( -name '*.fastq' -o -name '*.fastq.gz' \) | sort)
echo "=== Detected ${#FASTQS[@]} fastq file(s) ==="
printf '%s\n' "${FASTQS[@]}"

# Helper: cat or zcat depending on extension
read_cmd () { case "$1" in *.gz) zcat "$1";; *) cat "$1";; esac; }

for f in "${FASTQS[@]}"; do
  base="$(basename "$f")"
  echo "---------------------------------------------"
  echo "FILE: $base"

  # Read count (lines / 4)
  n=$(read_cmd "$f" | wc -l); echo "  reads: $(( n / 4 ))"

  # Read-length distribution (first 50k reads is plenty to characterize)
  read_cmd "$f" | awk 'NR%4==2{print length($0)}' | head -50000 \
    | sort -n | uniq -c | sort -rn | head -15 \
    > "$DIAG/${base}.lengthdist.txt"
  echo "  top read lengths (count length):"
  sed 's/^/    /' "$DIAG/${base}.lengthdist.txt" | head -5

  # Is the forward primer still inline at the 5' end? Check both 515F variants.
  echo "  primer scan (first 50k reads):"
  CAP=$(read_cmd "$f" | awk 'NR%4==2' | head -50000 | grep -c '^.\{0,12\}GTGCCAGC' || true)
  PAR=$(read_cmd "$f" | awk 'NR%4==2' | head -50000 | grep -c '^.\{0,12\}GTGYCAGC\|^.\{0,12\}GTG[CT]CAGC' || true)
  echo "    Caporaso 515F (GTGCCAGC...) leading matches: $CAP"
  echo "    515F-Parada    (GTGYCAGC...) leading matches: $PAR"

  # Is the reverse primer present (suggests merged/long reads spanning the amplicon)?
  REV=$(read_cmd "$f" | awk 'NR%4==2' | head -50000 | grep -c 'GGACTAC' || true)
  echo "    reverse-primer (GGACTAC...) anywhere matches: $REV"
done

echo "=============================================="
echo "INTERPRETATION GUIDE:"
echo " * One fastq, ~270 bp reads, forward primer inline, reverse primer also"
echo "   present  => single multiplexed OR merged single-end product spanning V4."
echo "             Use the SINGLE-END path in 02 (denoise-single)."
echo " * Many fastqs (one per sample), no barcodes inline => already demultiplexed;"
echo "   build a single-end manifest pointing at each file (02 covers this)."
echo " * Forward primer matches concentrate on Caporaso (GTGCCAGC), confirming"
echo "   the mapping file: trim THAT sequence in 02, not the 2026 Parada primer."
echo ""
echo "Record which case you are in, then run scripts/qiime2/02_qiime2_upstream.sh"
