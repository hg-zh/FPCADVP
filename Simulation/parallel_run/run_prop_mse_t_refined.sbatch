#!/bin/bash
#SBATCH --mail-type=END,FAIL
#SBATCH -N 1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=50
#SBATCH --mem=570G
#SBATCH -J para_R2_prop
#SBATCH -o para_R2_prop.out
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


Rscript para_R2_mse_t_refined.R \
  outdir=Res_LatentF_CE_095_8_b075_beta063 \
  reps=200 seed_base=2000 t_design=random \
  latent=F methodXi=CE \
  decayR=0.75 \
  save_mse_t=TRUE save_mse_plots=FALSE t_grid_len=51 mse_rep_only=NA \
  eps_u=0.05 u_cap=0.95 Z_cap=8 \
  "beta_a=0.6,3" \
  "beta_b=0.6,3"