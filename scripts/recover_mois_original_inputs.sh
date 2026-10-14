#!/usr/bin/env bash
# Recover ORIGINAL FCP full source inputs for the independent mois macro manuscript.
# This is download + identity verification ONLY, not automatic data redistribution,
# a source-model rerun, DOI deposit, or permission to publish third-party photographs.
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT="${MOIS_SOURCE_DIR:-$ROOT/raw_inputs}"
mkdir -p "$OUT"

if ! command -v gh >/dev/null 2>&1; then
  echo "ERROR: install GitHub CLI (gh) and authenticate with rights to read the FCP workflow artifacts" >&2
  exit 2
fi
if ! command -v sha256sum >/dev/null 2>&1; then
  echo "ERROR: sha256sum is needed for exact original-archive identity verification" >&2
  exit 2
fi

retrieve() {
  local artifact_id="$1" name="$2" digest="$3"
  local target="$OUT/$name" temp="$OUT/.$name.partial"
  if [[ ! -f "$target" ]]; then
    rm -f "$temp"
    echo "Fetching original FCP artifact $artifact_id -> $target"
    gh api -H 'Accept: application/vnd.github+json' \
      "/repos/zuizui0223/fcp/actions/artifacts/$artifact_id/zip" > "$temp"
    printf '%s  %s\n' "$digest" "$temp" | sha256sum -c -
    mv "$temp" "$target"
  fi
  printf '%s  %s\n' "$digest" "$target" | sha256sum -c -
  if command -v unzip >/dev/null 2>&1; then unzip -tq "$target" >/dev/null; fi
  echo "PASS: $name exact original FCP SHA256 and ZIP integrity"
}

retrieve 11587115220 "original_100543_photo_sites.zip" "2d6097ba1b0434632e13bbe79dd78180608a287ff5dabebc770b5c0a2cfb52a0"
retrieve 11653183273 "original_100543_photo_sites_all36_environment.zip" "6e287294b3a05eaa383be3de644551825915bfecca9c4bde6950b4890f5e16ca"

cat <<EOF
COMPLETE: two original full source ZIPs recovered and source SHA256 verified.
Files remain LOCAL in $OUT, excluded by repository .gitignore.
They were NOT redistributed, uploaded to Dryad/Zenodo, or published by this script.
Before journal peer review, independently check photo licensing, permissions,
code and linked data integrity, then create an author-approved anonymous stable
reviewer-access archive and record its private link in the blinded manuscript.
EOF
