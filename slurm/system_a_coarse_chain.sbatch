#!/usr/bin/env bash
# Three patients x three fixed contexts x four chains x two methods.
#SBATCH --job-name=sab-coarse-chain
#SBATCH --mail-user=reub0840@stats.ox.ac.uk
#SBATCH --mail-type=FAIL
#SBATCH --cluster=srf_cpu_01
#SBATCH --partition=standard-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --time=02:00:00
#SBATCH --chdir=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes
#SBATCH --output=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/coarse-chain-%j.out
#SBATCH --error=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/coarse-chain-%j.err
set -Eeuo pipefail
umask 027
WORKSPACE_ROOT="/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0"
export SAB_PROJECT_ROOT="$WORKSPACE_ROOT/projects/saem_anchored_bayes"
RSCRIPT_BIN="$WORKSPACE_ROOT/software/miniforge3/envs/model_fusion_r/bin/Rscript"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 SAB_SYSTEM_A_INTEGRATION=true
"$RSCRIPT_BIN" --vanilla -e 'testthat::test_dir("tests/testthat", reporter="summary", stop_on_failure=TRUE)'
exec "$RSCRIPT_BIN" --vanilla experiments/system_a_coarse_chain.R \
  "$SAB_PROJECT_ROOT/outputs/coarse_chain_${SLURM_JOB_ID}"
