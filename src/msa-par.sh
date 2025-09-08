#!/bin/bash -l
#SBATCH --account=project_2002833
#SBATCH --job-name=msa_array
#SBATCH --output=local_data/logs/msa_array/%x_%A_%a.out
#SBATCH --error=local_data/logs/msa_array/%x_%A_%a.stderr
#SBATCH --time=1:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8 
#SBATCH --mem-per-cpu=8G
#SBATCH --partition=small

set -euo pipefail

echo "=== Job started at $(date) ==="
echo "SLURM job ID: $SLURM_JOB_ID"
echo "Working dir: $(pwd)"

# Threads from Slurm:
THREADS="${SLURM_CPUS_PER_TASK:-1}"

# === Modules ===
module load mafft

echo "MAFFT version:"
mafft --version

# === Establish paths ===
set_file="$1"
in_dir="$2"
out_dir="$3"

# Get the line that matches this task ID (basename with .fa)
fname=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$set_file")
in_path="$in_dir/$fname"
base="${fname%.fa}"
out_path="$out_dir/${base}_mafft.fa"

# skip if already aligned
[[ -s "$out_path" ]] && { echo "Exists: $out_path — skipping."; exit 0; }

echo "Running: mafft --thread $THREADS --retree 2 --maxiterate 1000 $in_path > $out_path"
mafft --thread "$THREADS" --retree 2 --maxiterate 1000 "$in_path" > "$out_path"

echo "Job completed!"
echo "=== Job ended at $(date) ==="
