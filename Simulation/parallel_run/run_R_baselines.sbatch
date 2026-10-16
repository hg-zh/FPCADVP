#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=570G
#SBATCH -J para_R2_baselines
#SBATCH -o para_R2_baselines_mse_t.out
#SBATCH -t 120:00:00

# Compatibility launcher for the published IN benchmarks in both settings.
# Submit from the repository root, Code root, or this script's directory.
set -euo pipefail
submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"
code_dir="${SIMULATION_CODE_DIR:-}"
if [[ -z "$code_dir" ]]; then
  for candidate in "$submit_dir/Code/Simulation/parallel_run" "$submit_dir/Simulation/parallel_run" "$submit_dir"; do
    if [[ -f "$candidate/para_R2_mse_t_refined.R" ]]; then code_dir="$candidate"; break; fi
  done
fi
if [[ -z "$code_dir" || ! -f "$code_dir/para_R2_mse_t_refined.R" ]]; then
  echo "Cannot locate simulation code; set SIMULATION_CODE_DIR." >&2
  exit 2
fi
cd "$code_dir"
export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1 NUMEXPR_NUM_THREADS=1
# Load a site-specific module only when requested; otherwise use Rscript on PATH.
if [[ -n "${R_MODULE:-}" ]]; then module load "$R_MODULE"; fi
command -v Rscript >/dev/null || { echo "Rscript is unavailable" >&2; exit 2; }
bash "run_benchmarks_IN.sbatch"
