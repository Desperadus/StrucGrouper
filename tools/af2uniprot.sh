#!/usr/bin/env bash
# af2uniprot.sh
# Usage:
#   ./af2uniprot.sh input.tsv > output.tsv
#   # or pipe:
#   cat input.tsv | ./af2uniprot.sh > output.tsv

awk -v FS='\t' -v OFS='\t' '
NR==1 { print; next }  # keep header as-is
{
  # Match AF-<anything>-F<number>-model_v<number>
  if ($1 ~ /^AF-.*-F[0-9]+-model_v[0-9]+$/) {
    id = $1
    sub(/^AF-/, "", id)                          # drop leading "AF-"
    sub(/-F[0-9]+-model_v[0-9]+$/, "", id)       # drop trailing "-F*-model_v*"
    $1 = id
  }
  print
}' "${1:-/dev/stdin}"
