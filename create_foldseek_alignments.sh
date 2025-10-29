#!/usr/bin/env bash
# Create a Foldseek structure database and index from a directory of PDB/mmCIF files.
# Usage:
#   ./make_foldseek_db.sh -i afdb_pdb -d afdbDB -T foldseek_tmp -t 32
#
# Required:
#   -i  Path to input folder with structures (e.g., PDB/mmCIF files)
#
# Optional:
#   -d  Output database prefix/name (default: afdbDB)
#   -T  Temporary/index directory (default: foldseek_tmp)
#   -t  Threads (default: number of available CPUs)
#   -f  Force overwrite if DB already exists (default: off)

set -euo pipefail

# Defaults
INPUT=""
DB="afdbDB"
TMPDIR="foldseek_tmp"
THREADS="$(getconf _NPROCESSORS_ONLN 2>/dev/null || echo 1)"
FORCE=0

usage() {
  sed -n '2,25p' "$0" | sed 's/^# \{0,1\}//'
}

# Parse args
while getopts ":i:d:T:t:fh" opt; do
  case "${opt}" in
  i) INPUT="${OPTARG}" ;;
  d) DB="${OPTARG}" ;;
  T) TMPDIR="${OPTARG}" ;;
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

if ! command -v foldseek >/dev/null 2>&1; then
  echo "Error: 'foldseek' not found in PATH. Activate your conda env first:" >&2
  echo "  conda activate foldseek-umap" >&2
  exit 1
fi

if [[ ! -d "${INPUT}" ]]; then
  echo "Error: input folder '${INPUT}' does not exist." >&2
  exit 1
fi

echo "=== Foldseek DB builder ==="
echo "Input dir : ${INPUT}"
echo "DB name   : ${DB}"
echo "Tmp/index : ${TMPDIR}"
echo "Threads   : ${THREADS}"
echo "Force     : ${FORCE}"

# Create tmp dir
mkdir -p "${TMPDIR}"

# If DB already exists and not forcing, bail out
if [[ -e "${DB}.lookup" || -e "${DB}_h" || -e "${DB}_ss" || -e "${DB}.dbtype" ]]; then
  if [[ "${FORCE}" -eq 0 ]]; then
    echo "Database '${DB}' seems to already exist. Use -f to overwrite." >&2
    exit 1
  else
    echo "Forcing rebuild of existing DB '${DB}'..."
    rm -f "${DB}.lookup" "${DB}.dbtype"
    rm -rf "${DB}_h" "${DB}_ss" "${DB}_ca" 2>/dev/null || true
  fi
fi

echo "[1/2] Running: foldseek createdb '${INPUT}' '${DB}'"
foldseek createdb "${INPUT}" "${DB}"

echo "[2/2] Running: foldseek createindex '${DB}' '${TMPDIR}' --threads ${THREADS}"
foldseek createindex "${DB}" "${TMPDIR}" --threads "${THREADS}"

echo "Done. Database prefix: ${DB}"
echo "You can now run 'foldseek easy-search' against '${DB}'."
