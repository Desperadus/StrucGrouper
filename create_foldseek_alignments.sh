#!/usr/bin/env bash
# Build a Foldseek DB from a directory of structures and run an all-vs-all top-K search.
# Produces a TSV (.m8) with per-query top-K neighbors (preferably with alntmscore).
#
# Usage:
#   ./create_foldseek_alignments.sh -i afdb_pdb [-q afdb_pdb] [-d afdbDB] [-o knn_raw.m8] \
#                         [-T foldseek_tmp] [-k 50] [-t 32] [-f]
#
# Options:
#   -i  Input folder with PDB/mmCIF files (required)
#   -q  Query folder (default: same as -i; use to do queries vs DB)
#   -d  Output DB name/prefix (default: afdbDB)
#   -o  Output TSV path (default: knn_raw.m8)
#   -T  Temporary/index directory (default: foldseek_tmp)
#   -k  Top-K hits to keep per query (default: 50)
#   -t  Threads (default: number of CPUs)
#   -f  Force rebuild if DB exists
#
# Notes:
#   - Tries to run TM-align mode (--alignment-type 1).
#   - First attempts a format with 'alntmscore'; on failure, falls back to a
#     format without it, so older builds still work.

set -euo pipefail

INPUT=""
QUERY=""
DB="afdbDB"
OUT="knn_raw.m8"
TMPDIR="foldseek_tmp"
K=50
THREADS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)"
FORCE=0

usage() {
  sed -n '2,999p' "$0" | sed 's/^# \{0,1\}//' | awk '/^Usage:/,0' | sed '1,1!d'
}

while getopts ":i:q:d:o:T:k:t:fh" opt; do
  case "${opt}" in
  i) INPUT="${OPTARG}" ;;
  q) QUERY="${OPTARG}" ;;
  d) DB="${OPTARG}" ;;
  o) OUT="${OPTARG}" ;;
  T) TMPDIR="${OPTARG}" ;;
  k) K="${OPTARG}" ;;
  t) THREADS="${OPTARG}" ;;
  f) FORCE=1 ;;
  h)
    usage
    exit 0
    ;;
  \?)
    echo "Unknown option: -$OPTARG" >&2
    usage
    exit 1
    ;;
  :)
    echo "Option -$OPTARG requires an argument." >&2
    usage
    exit 1
    ;;
  esac
done

if [[ -z "${INPUT}" ]]; then
  echo "Error: -i <input_folder> is required." >&2
  usage
  exit 1
fi
[[ -z "${QUERY}" ]] && QUERY="${INPUT}"

if ! command -v foldseek >/dev/null 2>&1; then
  echo "Error: 'foldseek' not found in PATH. Activate your conda env first:" >&2
  echo "  conda activate foldseek-umap" >&2
  exit 1
fi

if [[ ! -d "${INPUT}" ]]; then
  echo "Error: input folder '${INPUT}' does not exist." >&2
  exit 1
fi
if [[ ! -d "${QUERY}" ]]; then
  echo "Error: query folder '${QUERY}' does not exist." >&2
  exit 1
fi

echo "=== Foldseek: Build & Search ==="
echo "Input dir : ${INPUT}"
echo "Query dir : ${QUERY}"
echo "DB name   : ${DB}"
echo "Output    : ${OUT}"
echo "Tmp/index : ${TMPDIR}"
echo "Top-K     : ${K}"
echo "Threads   : ${THREADS}"
echo "Force     : ${FORCE}"

mkdir -p "${TMPDIR}"

# If DB already exists and not forcing, bail out
if [[ -e "${DB}.lookup" || -e "${DB}.dbtype" || -d "${DB}_h" || -d "${DB}_ss" || -d "${DB}_ca" ]]; then
  if [[ "${FORCE}" -eq 0 ]]; then
    echo "Database '${DB}' already exists. Use -f to overwrite." >&2
  else
    echo "Forcing DB rebuild of '${DB}'..."
    rm -f "${DB}.lookup" "${DB}.dbtype" || true
    rm -rf "${DB}_h" "${DB}_ss" "${DB}_ca" 2>/dev/null || true
  fi
fi

# Build DB and index
if [[ ! -e "${DB}.dbtype" ]]; then
  echo "[1/3] foldseek createdb '${INPUT}' '${DB}'"
  foldseek createdb "${INPUT}" "${DB}"
else
  echo "[1/3] DB files exist, skipping createdb."
fi

echo "[2/3] foldseek createindex '${DB}' '${TMPDIR}' --threads ${THREADS}"
foldseek createindex "${DB}" "${TMPDIR}" --threads "${THREADS}"

# Run easy-search with TM-align mode and try preferred output format with alntmscore
echo "[3/3] foldseek easy-search (TM-align mode, top-K=${K})"
echo "      Writing: ${OUT}"

set +e
foldseek easy-search "${QUERY}" "${DB}" "${OUT}" "${TMPDIR}" \
  --threads "${THREADS}" \
  --max-seqs "${K}" \
  --alignment-type 1 \
  --format-output "query,target,alnlen,qcov,tcov,evalue,bits,alntmscore"
rc=$?
set -e

if [[ $rc -ne 0 ]]; then
  echo "   > 'alntmscore' not supported by this Foldseek build. Falling back..."
  foldseek easy-search "${QUERY}" "${DB}" "${OUT}" "${TMPDIR}" \
    --threads "${THREADS}" \
    --max-seqs "${K}" \
    --alignment-type 1 \
    --format-output "query,target,alnlen,qcov,tcov,evalue,score"
fi

echo "Done."
echo "TSV written to: ${OUT}"
echo "Next step: run your Python UMAP script on '${OUT}'."
