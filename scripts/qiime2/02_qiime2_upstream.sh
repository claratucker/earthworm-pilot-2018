#!/usr/bin/env bash
# 02_qiime2_upstream.sh -- single-end V4 processing for the 2018 pilot.
#
# Matches the 2026 toolchain (QIIME 2 / DADA2 / SILVA 138) so the two datasets
# are comparable, but uses the SINGLE-END path and the original Caporaso 515F
# primer that 01_diagnose.sh confirmed. Long steps run inside `screen`.
#
# Bolyen et al. 2019 Nat. Biotechnol. 37:852 (QIIME 2)
# Caporaso et al. 2011 PNAS 108:4516 (515F/806R primers)
# Martin 2011 EMBnet J. 17:10 (cutadapt)
# Callahan et al. 2016 Nat. Methods 13:581 (DADA2)
# Quast et al. 2013 NAR 41:D590; Bokulich et al. 2018 Microbiome 6:90;
#   Robeson et al. 2021 PLoS Comput Biol 17:e1009581 (SILVA 138 + classifier)
# Katoh & Standley 2013 MBE 30:772 (MAFFT); Price et al. 2010 PLoS ONE 5:e9490 (FastTree)
set -euo pipefail

PROJECT="$HOME/pilot2018"
BUCKET="s3://earthworm-pilot-2018"
THREADS="${THREADS:-4}"
Q="$PROJECT/results/qiime2"
cd "$PROJECT"
source /home/ubuntu/miniforge3/etc/profile.d/conda.sh && conda activate qiime2-amplicon-2026.1

# Primers CONFIRMED from the 2018 mapping file:
FWD_PRIMER="GTGCCAGCMGCCGCGGTAA"      # Caporaso 515F (note: C at pos 4, not Y)
REV_PRIMER="GGACTACVSGGGTATCTAAT"     # 806R variant from mapping

# =============================================================================
# IMPORT + DEMULTIPLEX
# CONFIRMED from the raw reads: this is ONE multiplexed single-end fastq with the
# 8 bp barcode PREPENDED to each read, immediately followed by the 515F primer:
#   [8bp barcode][GTGCCAGCMGCCGCGGTAA][16S...]
# e.g. read 1 = GAGTTCTG + GTGCCAGC... = sample C2.Soil.
# So we import the single multiplexed file, then use cutadapt to demultiplex on
# the inline 5' barcodes AND strip the primer in the same workflow.
#
# Step 1: import the single multiplexed fastq as MultiplexedSingleEndBarcodeInSequence
cp "$PROJECT/data/raw/122117ND515F2-full.fastq" "$Q/forward.fastq"
gzip -f "$Q/forward.fastq"           # QIIME importer expects .fastq.gz
qiime tools import \
  --type MultiplexedSingleEndBarcodeInSequence \
  --input-path "$Q/forward.fastq.gz" \
  --output-path "$Q/multiplexed-seqs.qza"

# Step 2: build the barcode file QIIME expects (sample-id + barcode columns).
#         metadata.tsv already has a 'barcode' column, so reuse it directly.
# Step 3: cutadapt demux on inline 5' barcodes. This removes the barcode; we then
#         remove the primer in a second cutadapt pass for clarity/logging.
qiime cutadapt demux-single \
  --i-seqs "$Q/multiplexed-seqs.qza" \
  --m-barcodes-file "$PROJECT/data/metadata.tsv" \
  --m-barcodes-column barcode \
  --p-error-rate 0 \
  --o-per-sample-sequences "$Q/demux.qza" \
  --o-untrimmed-sequences "$Q/demux-untrimmed.qza" \
  --verbose 2>&1 | tee "$Q/demux-log.txt"

qiime demux summarize --i-data "$Q/demux.qza" \
  --o-visualization "$Q/demux-summary.qzv"

# =============================================================================
# PRIMER REMOVAL (cutadapt, single-end). The M degeneracy in 515F resolves to
# both A and C in the reads (GTGCCAGCA... and GTGCCAGCC...); cutadapt handles the
# IUPAC code automatically. --discard-untrimmed drops anything without the primer.
# =============================================================================
qiime cutadapt trim-single \
  --i-demultiplexed-sequences "$Q/demux.qza" \
  --p-front "$FWD_PRIMER" \
  --p-discard-untrimmed --p-cores "$THREADS" \
  --o-trimmed-sequences "$Q/trimmed.qza" \
  --verbose 2>&1 | tee "$Q/cutadapt-log.txt"

# =============================================================================
# DADA2 DENOISE (single-end). Set --p-trunc-len from the demux-summary quality
# plot. Reads are ~272 bp incl. primer; after trimming the ~19 bp primer the
# insert is ~253 bp. A trunc-len of 0 (no truncation) is acceptable if quality
# holds to the 3' end; otherwise truncate where median quality drops below ~Q25.
# START at 0, inspect denoising-stats, raise truncation only if retention is low.
# =============================================================================
TRUNC_LEN="${TRUNC_LEN:-0}"
qiime dada2 denoise-single \
  --i-demultiplexed-seqs "$Q/trimmed.qza" \
  --p-trunc-len "$TRUNC_LEN" \
  --p-max-ee 2 \
  --p-n-threads "$THREADS" \
  --o-table "$Q/table.qza" \
  --o-representative-sequences "$Q/rep-seqs.qza" \
  --o-denoising-stats "$Q/denoising-stats.qza" \
  --o-base-transition-stats "$Q/base-transition-stats.qza" \

qiime metadata tabulate \
  --m-input-file "$Q/denoising-stats.qza" \
  --o-visualization "$Q/denoising-stats.qzv"

# =============================================================================
# TAXONOMY (SILVA 138 -- SAME classifier artifact as the 2026 run)
# =============================================================================
qiime feature-classifier classify-sklearn \
  --i-classifier "$PROJECT/refs/silva-138-99-nb-classifier.qza" \
  --i-reads "$Q/rep-seqs.qza" \
  --p-n-jobs "$THREADS" \
  --o-classification "$Q/taxonomy.qza"

# Remove mitochondria/chloroplast. CRITICAL for gut samples: earthworm host
# mitochondrial 16S co-amplifies with universal primers and otherwise appears
# as a spurious dominant "taxon". (Drake & Horn 2007 Annu Rev Microbiol 61:169
# for the gut-as-filter context.)
qiime taxa filter-table \
  --i-table "$Q/table.qza" --i-taxonomy "$Q/taxonomy.qza" \
  --p-exclude mitochondria,chloroplast --p-mode contains \
  --o-filtered-table "$Q/table-nohost.qza"

# Low-frequency feature filter (lenient -- this is a small pilot)
qiime feature-table filter-features \
  --i-table "$Q/table-nohost.qza" \
  --p-min-frequency 5 --p-min-samples 1 \
  --o-filtered-table "$Q/table-final.qza"

# =============================================================================
# PICRUSt2 -- predicts gene-family (KO) content per ASV from phylogenetic
# placement. Needed by the EPSPS exploratory module (aroA = K00800). This is
# PREDICTION, not measurement: PICRUSt2 infers gene content from 16S, it does
# not observe EPSPS. (Douglas et al. 2020 Nat Biotechnol 38:685.)
#
# IMPORTANT (lesson from the 2026 run): the q2-picrust2 PLUGIN often will not
# load inside the qiime2-amplicon-2026.1 env due to a version mismatch (the
# plugin is built against 2024.5). The reliable route, validated on the 2026
# run, is to EXPORT the table + rep-seqs and run STANDALONE picrust2_pipeline.py
# from a dedicated picrust2 conda env. That is the default below; the plugin
# call is kept commented as a fallback if your env has a matching plugin.
# =============================================================================
PIC_OUT="$PROJECT/results/picrust2"
PIC_EXPORT="$PROJECT/results/picrust2_input"; mkdir -p "$PIC_EXPORT"

# Export the inputs PICRUSt2 needs (biom table + rep-seq fasta)
qiime tools export --input-path "$Q/table-final.qza" --output-path "$PIC_EXPORT"
qiime tools export --input-path "$Q/rep-seqs.qza"    --output-path "$PIC_EXPORT"
# -> $PIC_EXPORT/feature-table.biom and $PIC_EXPORT/dna-sequences.fasta

# Standalone PICRUSt2 (run in its own env; do NOT nest conda activate inside the
# qiime env without deactivating first). Open a fresh shell or use `conda run`.
cat <<EOF

>>> Run PICRUSt2 in its dedicated env (matches the 2026 working route):
    conda deactivate
    conda activate picrust2            # or: q2-picrust2-qiime2-2024.5
    picrust2_pipeline.py \\
      -i "$PIC_EXPORT/feature-table.biom" \\
      -s "$PIC_EXPORT/dna-sequences.fasta" \\
      -o "$PIC_OUT" \\
      -p "$THREADS" --verbose
    conda deactivate && source /home/ubuntu/miniforge3/etc/profile.d/conda.sh && conda activate qiime2-amplicon-2026.1
<<< Then re-run this script's remaining steps, or just continue manually.

EOF

# --- Fallback: the plugin form (only if your env has a compatible q2-picrust2) ---
# qiime picrust2 full-pipeline \
#   --i-table "$Q/table-final.qza" --i-seq "$Q/rep-seqs.qza" \
#   --p-threads "$THREADS" --p-hsp-method mp --p-max-nsti 2.0 \
#   --output-dir "$PIC_OUT" --verbose

# Standalone PICRUSt2 writes KO predictions to:
#   $PIC_OUT/KO_metagenome_out/pred_metagenome_unstrat.tsv.gz   (per-sample KO)
#   $PIC_OUT/KO_predicted.tsv.gz                                (per-ASV KO; aroA)
#   $PIC_OUT/marker_predicted_and_nsti.tsv.gz                   (per-ASV NSTI)
# The EPSPS module (04) looks for KO_predicted.tsv* under results/picrust2.

# =============================================================================
# PHYLOGENY (for UniFrac / Faith's PD). Produces BOTH unrooted and rooted .qza
# trees. phyloseq needs the ROOTED tree -- the export below pulls rooted-tree.qza
# and `qiime tools export` writes it as a plain Newick file named tree.nwk.
# (This is the "where does tree.nwk come from" answer: it is the EXPORTED form of
#  rooted-tree.qza, not a separate phylogeny output. Always export the ROOTED
#  tree, never unrooted-tree.qza, or UniFrac/Faith's PD will be wrong.)
# =============================================================================
qiime phylogeny align-to-tree-mafft-fasttree \
  --i-sequences "$Q/rep-seqs.qza" --p-n-threads "$THREADS" \
  --o-alignment "$Q/aligned-rep-seqs.qza" \
  --o-masked-alignment "$Q/masked-aligned-rep-seqs.qza" \
  --o-tree "$Q/unrooted-tree.qza" \
  --o-rooted-tree "$Q/rooted-tree.qza"

# =============================================================================
# EXPORT for R
# =============================================================================
EXP="$PROJECT/results/r/exported"; mkdir -p "$EXP"
qiime tools export --input-path "$Q/table-final.qza" --output-path "$EXP"
biom convert -i "$EXP/feature-table.biom" -o "$EXP/feature-table.tsv" --to-tsv
qiime tools export --input-path "$Q/taxonomy.qza"     --output-path "$EXP"
qiime tools export --input-path "$Q/rooted-tree.qza"  --output-path "$EXP"
# -> $EXP/tree.nwk  (this is the rooted tree, exported to Newick; 03/05 read it)
test -f "$EXP/tree.nwk" || echo "WARNING: tree.nwk not found; check rooted-tree.qza export"

aws s3 sync "$PROJECT/results/" "$BUCKET/results/qiime2-artifacts/"
echo "Upstream done. Next: Rscript scripts/r/03_pilot_analysis.R"
