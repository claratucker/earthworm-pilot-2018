#!/usr/bin/env bash
# 99_init_github_repo.sh -- turn this folder into a GitHub repo.
#
# Run this LOCALLY (on the instance or wherever the scripts live). It does NOT
# commit data: the .gitignore keeps fastq/qza/results out, so only code, the
# metadata table, and the procedure doc are tracked. That is exactly what you
# want in version control; the data stays in s3://earthworm-pilot-2018.
set -euo pipefail

REPO_NAME="${REPO_NAME:-earthworm-2018-pilot}"
cd "$(dirname "$0")/.."   # repo root (the pilot2018 folder)

git init -b main
git add .gitignore README.md AWS_SETUP_STEP_BY_STEP.md data/metadata.tsv scripts/ *.docx
git status   # review what is staged BEFORE committing
git commit -m "Initial commit: 2018 earthworm gut/soil glyphosate pilot pipeline"

cat <<EOF

Next steps (choose ONE):

A) With the GitHub CLI (easiest):
     gh repo create $REPO_NAME --public --source=. --remote=origin --push

B) Manually:
     1. Create an empty repo named "$REPO_NAME" at https://github.com/new
        (do NOT initialize it with a README/license -- you already have them)
     2. git remote add origin git@github.com:<your-username>/$REPO_NAME.git
     3. git push -u origin main

To confirm data is NOT being pushed, check that 'git status' above showed no
.fastq / .qza / results/ files staged.
EOF
