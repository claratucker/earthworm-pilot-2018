#!/usr/bin/env bash
# 00_setup_and_stage.sh -- one-time setup for the 2018 earthworm gut/soil pilot.
# Creates the S3 bucket, the project tree on the EC2 instance, and pulls the
# raw data + mapping file down from S3 onto the instance.
#
# RUN CONTEXT:
#   - Launch an EC2 instance from your existing AMI: qiime2-earthworm-base
#     (the same image used for the 2026 dose-gradient run; it already has the
#      QIIME 2 conda env, awscli, and R toolchain baked in).
#   - Suggested instance type: t3.xlarge (4 vCPU, 16 GB RAM). The dataset is
#     small (~1.6M reads total), but the SILVA classifier step needs ~12-16 GB
#     RAM regardless of read count, so 16 GB is the floor. A t3.large (8 GB)
#     will OOM during classify-sklearn.
#   - This dataset is single-end V4, so there is no DADA2 paired-merge step and
#     the run is fast; you will NOT need an m7i.2xlarge like the 2026 run.
#
# BEFORE YOU RUN:
#   On your LAPTOP, after downloading the unzipped files out of Google Drive:
#     aws s3 mb s3://earthworm-pilot-2018
#     aws s3 cp 122117ND515F2-full.fastq      s3://earthworm-pilot-2018/raw-data/
#     aws s3 cp 122117ND515F4-mapping.txt     s3://earthworm-pilot-2018/raw-data/
#   (Upload whatever the real read file turns out to be -- see 01_diagnose.sh.
#    If there are per-sample fastqs instead of one multiplexed file, sync the
#    whole folder:  aws s3 sync ./per_sample_fastqs s3://earthworm-pilot-2018/raw-data/per_sample/ )
set -euo pipefail

BUCKET="s3://earthworm-pilot-2018"
PROJECT="$HOME/pilot2018"

# ---- Project tree (mirrors the 2026 layout: data / results/qiime2 / results/r)
mkdir -p "$PROJECT"/{data,results/qiime2,results/r,results/diagnostics,refs}
cd "$PROJECT"

# ---- Pull raw data + mapping down from S3 ------------------------------------
aws s3 sync "$BUCKET/raw-data/" "$PROJECT/data/raw/"

echo "Staged raw data into $PROJECT/data/raw/ :"
ls -lh "$PROJECT/data/raw/"

# ---- Reference classifier ----------------------------------------------------
# Reuse the SAME SILVA 138 classifier artifact as the 2026 run so the two
# datasets are taxonomically comparable. If it is already baked into the AMI,
# copy it into refs/; otherwise pull from your 2026 bucket.
if [ -f "$HOME/project/data/silva-138-99-nb-classifier.qza" ]; then
  cp "$HOME/project/data/silva-138-99-nb-classifier.qza" "$PROJECT/refs/"
else
  echo "NOTE: SILVA classifier not found at the 2026 path."
  echo "      Pull it from your 2026 bucket, e.g.:"
  echo "      aws s3 cp s3://earthworm-microbiome-2026/refs/silva-138-99-nb-classifier.qza $PROJECT/refs/"
fi

echo "Setup complete. Next: bash scripts/qiime2/01_diagnose.sh"
