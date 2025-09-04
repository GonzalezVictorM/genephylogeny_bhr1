#!/usr/bin/env bash
set -euo pipefail
shopt -s nullglob
export LC_ALL=C

# === Define directories ===
DATA_DIR="local_data"
SET_DIR="$DATA_DIR/speciestree"
IN_DIR="$SET_DIR/seq_alignments"
OUT_DIR="$SET_DIR/trimmed_alignments"
LOGS_DIR="$DATA_DIR/logs/trim_array"

mkdir -p "$OUT_DIR" "$LOGS_DIR"

# Resolve to absolute paths (safer for SLURM nodes)
IN_DIR="$(readlink -f "$IN_DIR")"
OUT_DIR="$(readlink -f "$OUT_DIR")"
SET_DIR="$(readlink -f "$SET_DIR")"
LOGS_DIR="$(readlink -f "$LOGS_DIR")"

SET_FILE="$SET_DIR/untrimmed_files.txt"

# === Cleanup: Delete any empty alignment files ===
echo "Cleaning up empty files in $OUT_DIR..."
find "$OUT_DIR" -type f -empty -delete || true
echo "Cleanup complete."

# === Compute pending files ===
mapfile -t in_basenames  < <(find "$IN_DIR"  -maxdepth 1 -type f -name '*.fa'      -printf '%f\n' | sed -E 's/\.fa$//'      | sort -u)
mapfile -t out_basenames < <(find "$OUT_DIR" -maxdepth 1 -type f -name '*_trim.fa' -printf '%f\n' | sed -E 's/_trim\.fa$//' | sort -u)

in_tmp=$(mktemp); out_tmp=$(mktemp)
printf "%s\n" "${in_basenames[@]:-}"  > "$in_tmp"
printf "%s\n" "${out_basenames[@]:-}" > "$out_tmp"

# Write missing names (with .fa restored) to the list file
comm -23 "$in_tmp" "$out_tmp" | sed 's/$/.fa/' > "$SET_FILE"

rm -f "$in_tmp" "$out_tmp"

if [[ -s "$SET_FILE" ]]; then
  npending=$(wc -l < "$SET_FILE")
  echo "There are $npending unaligned files."

  cap=${MAX_ARRAY_SIZE:-380}
  (( npending > cap )) && { echo "Capping array size to $cap to respect submit limit."; npending=$cap; }
  echo "Submitting SLURM array for $npending unaligned files…"
  SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
  jobid=$(sbatch --parsable --array=1-"$npending"%100 "$SCRIPT_DIR/trim-par.sh" "$SET_FILE" "$IN_DIR" "$OUT_DIR")
  echo "Submitted job $jobid"
else
  echo "All files in $IN_DIR have corresponding alignments in $OUT_DIR."
fi

