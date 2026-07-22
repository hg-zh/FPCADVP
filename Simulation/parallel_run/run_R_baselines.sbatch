#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=570G
#SBATCH -J para_R2_baselines
#SBATCH -o para_R2_baselines_mse_t.out
#SBATCH -t 120:00:00

# Run from the submission directory
cd "$SLURM_SUBMIT_DIR"

# Avoid thread over-subscription if BLAS/OpenMP are multithreaded
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

module load r/4.5.0



Rscript para_baselines_mse_t.R \
  outdir=Res_baselines_mse_t_b075_beta063 \
  save_mse_t=TRUE \
  plot_mse_t=FALSE \
  t_grid_len=51 \
  methodXi=CE \
  decayR=1.5 \
  Ku=10 basis_type=bspline bs_degree=3 \
  beta_a=0.6,3 \
  beta_b=0.6,3

