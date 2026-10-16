#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=300G
#SBATCH --time=120:00:00
#SBATCH --job-name=fig3_benchIN
#SBATCH --output=benchmarks_IN_%j.out
#SBATCH --error=benchmarks_IN_%j.err

# Published LQD and Tangent benchmarks use fdapace IN.
# Keep both data generators, evaluation grids, and seeds fixed.
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
result_root="${RESULTS_ROOT:-results/paper_comparison}/benchmarks_IN"
for setting in 1 2; do
  if [[ -d "$result_root/setting$setting" ]]; then
    echo "Refusing to overwrite: $result_root/setting$setting" >&2
    exit 1
  fi
done
for setting in 1 2; do
  if [[ "$setting" == 1 ]]; then
    beta_a=0.6,3
    beta_b=0.6,3
  else
    beta_a=3,4
    beta_b=1,2
  fi
  echo "Starting IN benchmarks, setting $setting, at $(date -Is)"
  Rscript para_baselines_mse_t.R \
    outdir="$result_root/setting$setting" reps=200 seed_base=2000 t_design=random \
    n_vec=20,50,100,200,400 N_vec=3,5,10,20 m_vec=200 \
    methodXi=IN decayR=1 Ku=10 basis_type=bspline bs_degree=3 \
    beta_a="$beta_a" beta_b="$beta_b" t_grid_len=51 \
    save_each_rep=TRUE save_mse_t=TRUE plot_mse_t=FALSE
  Rscript -e 'r <- readRDS(commandArgs(TRUE)[1]); stopifnot(nrow(r)==4000L, all(r$ok), all(is.finite(r$imse_sq_lqd)), all(is.finite(r$imse_sq_tangent))); cat("PASS: all 4000 benchmark records succeeded\n")' \
    "$result_root/setting$setting/raw_results.rds"
done
echo "Both IN benchmark settings completed at $(date -Is)"
