#!/bin/bash -l
#SBATCH --account=project_2002833
#SBATCH --job-name=iqtree_array
#SBATCH --output=local_data/logs/iqtree_array/%x_%A_%a.out
#SBATCH --error=local_data/logs/iqtree_array/%x_%A_%a.stderr
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=16
#SBATCH --mem-per-cpu=2G
#SBATCH --partition=small

set -euo pipefail

echo "=== Job started at $(date) ==="
echo "SLURM job ID: $SLURM_JOB_ID"
echo "Working dir: $(pwd)"

# Threads from Slurm:
THREADS="${SLURM_CPUS_PER_TASK:-1}"

# === Define paths to executables ===
IQTREE_BIN="/scratch/project_2002833/VG/software/iqtree-3.0.1-Linux/bin/iqtree3"

echo "IQTREE version:"
"$IQTREE_BIN" --version

# === Establish paths ===
set_file="$1"
in_dir="$2"
out_dir="$3"

# Get the line that matches this task ID (basename with .fa)
fname=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$set_file")
in_path="$in_dir/$fname"
base="${fname%.fa}"
out_path="$out_dir/${base}.treefile"

# skip if already aligned
[[ -s "$out_path" ]] && { echo "Exists: $out_path — skipping."; exit 0; }

echo "Running $IQTREE_BIN -s $in_path -m TEST -T $THREADS -B 1000 -alrt 1000 -pre $out_dir/$base"
"$IQTREE_BIN" -s "$in_path" -m TEST -T $THREADS -B 1000 -alrt 1000 -pre "$out_dir/$base"


echo "Job completed!"
echo "=== Job ended at $(date) ==="
