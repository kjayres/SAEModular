#!/usr/bin/env bash
# Three patients, one context, <=2500 sealed calls; no MCMC bank or sampler.
#SBATCH --job-name=sab-ceiling
#SBATCH --mail-user=reub0840@stats.ox.ac.uk
#SBATCH --mail-type=FAIL
#SBATCH --cluster=srf_cpu_01
#SBATCH --partition=standard-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:20:00
#SBATCH --chdir=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes
#SBATCH --output=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/ceiling-%j.out
#SBATCH --error=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/ceiling-%j.err
set -Eeuo pipefail
umask 027
WORKSPACE_ROOT="/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0"
export SAB_PROJECT_ROOT="$WORKSPACE_ROOT/projects/saem_anchored_bayes"
RSCRIPT_BIN="$WORKSPACE_ROOT/software/miniforge3/envs/model_fusion_r/bin/Rscript"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$RSCRIPT_BIN" --vanilla -e 'testthat::test_dir("tests/testthat", filter="residual_ceiling", reporter="summary", stop_on_failure=TRUE)'
exec "$RSCRIPT_BIN" --vanilla experiments/system_a_residual_ceiling.R \
  "$SAB_PROJECT_ROOT/outputs/residual_ceiling_${SLURM_JOB_ID}"
