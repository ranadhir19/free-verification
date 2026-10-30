#!/bin/sh
# Reassemble data/chembl36_activities.csv.gz from its two byte parts.
# The single file (28,256,850 bytes) exceeds GitHub's 25 MB web-upload
# limit, so the GitHub copy stores it split. The reassembled file is
# byte-identical to the deposited artifact; its SHA-256 is verified here.
set -e
cd "$(dirname "$0")"
cat data/chembl36_activities.csv.gz.part00 data/chembl36_activities.csv.gz.part01 > data/chembl36_activities.csv.gz
echo "1583f55cbb72676ff2a8f7a7083fcc27c3f2f5825162170d237445f81be45af7  data/chembl36_activities.csv.gz" | sha256sum -c -
echo "join_data: data/chembl36_activities.csv.gz reassembled and verified."
