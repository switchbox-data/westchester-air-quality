#!/usr/bin/env bash
# Download the Aclima mobile-monitoring data from the public S3 bucket.
#
#   scripts/fetch_data.sh                  # Westchester pipeline input only (157 MB)
#   scripts/fetch_data.sh segments tvoc    # add statewide folders (see README)
#   scripts/fetch_data.sh all              # everything (about 2.7 GB)
#
# No AWS account is needed. Each file is checked against
# scripts/data_manifest.sha256; files that are already present and correct are
# skipped. Westchester files go to Data/; statewide files go to Data/aclima/.
set -euo pipefail

BASE_URL="https://switchbox-ny-air-monitoring.s3.us-west-2.amazonaws.com"
ROOT="$(cd "$(dirname "$0")/.." && pwd)"
MANIFEST="$ROOT/scripts/data_manifest.sha256"

if command -v sha256sum >/dev/null; then sha() { sha256sum "$1" | cut -d' ' -f1; }
else sha() { shasum -a 256 "$1" | cut -d' ' -f1; }; fi

folders=("westchester" "$@")
[[ " $* " == *" all "* ]] && folders=(westchester segments tvoc raw_1s)
for f in "${folders[@]}"; do
  case "$f" in westchester|segments|tvoc|raw_1s) ;; *)
    echo "unknown folder: $f (use segments, tvoc, raw_1s or all)" >&2; exit 1 ;;
  esac
done

while read -r hash key; do
  folder=${key%%/*}
  [[ " ${folders[*]} " == *" $folder "* ]] || continue
  if [[ $folder == westchester ]]; then dest="$ROOT/Data/${key#*/}"; else dest="$ROOT/Data/aclima/$key"; fi
  if [[ -f $dest && $(sha "$dest") == "$hash" ]]; then echo "ok       $key"; continue; fi
  echo "fetch    $key"
  mkdir -p "$(dirname "$dest")"
  curl --fail --silent --show-error --location --retry 3 -o "$dest.part" "$BASE_URL/$key"
  if [[ $(sha "$dest.part") != "$hash" ]]; then
    echo "checksum mismatch: $key" >&2; rm -f "$dest.part"; exit 1
  fi
  mv "$dest.part" "$dest"
done < "$MANIFEST"
