#!/usr/bin/env bash
# 00b_drive_to_s3_via_rclone.sh -- move ONLY the needed files from Google Drive
# to S3, running entirely ON THE EC2 INSTANCE so your laptop storage is never
# used. Drive -> instance disk -> S3, then delete the local copy.
#
# WHY rclone: it streams from Drive using your Google login (OAuth), so large
# files never pass through your laptop. The instance has the disk and bandwidth.
#
# ONE-TIME rclone SETUP (interactive, ~2 min):
#   sudo apt-get update && sudo apt-get install -y rclone   # (or: curl https://rclone.org/install.sh | sudo bash)
#   rclone config
#     n) new remote
#     name> gdrive
#     storage> drive            (Google Drive)
#     client_id> (leave blank or use your own for speed)
#     scope> 1                  (full access; read-only "2" also fine here)
#     Use auto config? > N      (you are on a headless EC2 box)
#       -> it prints a URL; open it on ANY device, log in, paste the token back
#     Configure as team drive? > N
#   That creates the "gdrive:" remote.
#
# TIP: rather than hunt folder names, the EASIEST path is to put the 3 needed
# files into ONE Drive folder (e.g. "earthworm_pilot_min") on your end first,
# then point this script at that folder.
set -euo pipefail

BUCKET="s3://earthworm-pilot-2018"
STAGE="$HOME/pilot2018/data/raw"
mkdir -p "$STAGE"

# ---- Option A: you made a single clean Drive folder with just the 3 files ----
# Replace with your folder name:
DRIVE_FOLDER="${DRIVE_FOLDER:-earthworm_pilot_min}"

echo "Copying from gdrive:$DRIVE_FOLDER -> $STAGE"
rclone copy "gdrive:$DRIVE_FOLDER" "$STAGE" --drive-acknowledge-abuse -P

# ---- Option B: copy specific files by name from wherever they live ----------
# If you did NOT make a clean folder, copy individual files instead. rclone can
# match by name across your drive with --drive-shared-with-me / a path. Example:
#   rclone copy "gdrive:path/to/122117ND515F2-full.fastq" "$STAGE" -P
#   rclone copy "gdrive:path/to/122117ND515F4-mapping.txt" "$STAGE" -P

echo "Files staged locally:"
ls -lh "$STAGE"

# ---- Push to S3, then free the instance disk --------------------------------
aws s3 sync "$STAGE" "$BUCKET/raw-data/"
echo "Uploaded to $BUCKET/raw-data/. Verifying:"
aws s3 ls "$BUCKET/raw-data/"

# Optional: reclaim space on the instance once it is safely in S3
# rm -f "$STAGE"/*.fastq
echo "Done. The data is in S3 and never touched your laptop."
