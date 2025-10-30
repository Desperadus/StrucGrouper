#!/usr/bin/env bash
# af2uniprot.sh
# Usage:
#   ./af2uniprot.sh input.tsv > output.tsv
#   # or pipe:
#   cat input.tsv | ./af2uniprot.sh > output.tsv

awk -F'\t' 'BEGIN{OFS=FS}
NR==1 { print; next }  # keep header as-is
{
  # Only transform when it matches the AlphaFold pattern
  if ($1 ~ /^AF-[^-]+-F[0-9]+-model_v[0-9]+$/) {
    # AF-<UniProt>-F*-model_v*  → take the 2nd dash-separated token
    split($1, a, "-");
    $1 = a[2];
  }
  print
}' "${1:-/dev/stdin}"
