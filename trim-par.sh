#!/bin/bash -l
#SBATCH --account=project_2002833
#SBATCH --job-name=trim_array
#SBATCH --output=local_data/logs/trim_array/%x_%A_%a.out
#SBATCH --error=local_data/logs/trim_array/%x_%A_%a.stderr
#SBATCH --time=00:10:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1 
#SBATCH --mem-per-cpu=1G
#SBATCH --partition=small

set -euo pipefail

echo "=== Job started at $(date) ==="
echo "SLURM job ID: $SLURM_JOB_ID"
echo "Working dir: $(pwd)"

# Threads from Slurm:
THREADS="${SLURM_CPUS_PER_TASK:-1}"

# === Define paths to executables ===
TRIMAL_BIN="/scratch/project_2002833/VG/software/trimal-1.5.0/source/trimal"

echo "trimAl version:"
"$TRIMAL_BIN" --version

# === Establish paths ===
set_file="$1"
in_dir="$2"
out_dir="$3"

# Get the line that matches this task ID (basename with .fa)
fname=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$set_file")
in_path="$in_dir/$fname"
base="${fname%.fa}"
out_path="$out_dir/${base}_trim.fa"

# skip if already aligned
[[ -s "$out_path" ]] && { echo "Exists: $out_path — skipping."; exit 0; }

echo "Running $TRIMAL_BIN -in $in_path -gt 0.8 -cons 10 -out $out_path"
"$TRIMAL_BIN" -in "$in_path" -gt 0.8 -cons 10 -out "$out_path"

echo "Job completed!"
echo "=== Job ended at $(date) ==="
