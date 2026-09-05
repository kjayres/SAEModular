#!/usr/bin/env bash
# Three patients, three shared contexts; no retained chain unless separately gated.
#SBATCH --job-name=sab-coarse
#SBATCH --mail-user=reub0840@stats.ox.ac.uk
#SBATCH --mail-type=FAIL
#SBATCH --cluster=srf_cpu_01
#SBATCH --partition=standard-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --chdir=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes
#SBATCH --output=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/coarse-%j.out
#SBATCH --error=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/coarse-%j.err
set -Eeuo pipefail
umask 027
WORKSPACE_ROOT="/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0"
export SAB_PROJECT_ROOT="$WORKSPACE_ROOT/projects/saem_anchored_bayes"
RSCRIPT_BIN="$WORKSPACE_ROOT/software/miniforge3/envs/model_fusion_r/bin/Rscript"
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 SAB_SYSTEM_A_INTEGRATION=true
"$RSCRIPT_BIN" --vanilla -e 'testthat::test_dir("tests/testthat", filter="coarse_guide", reporter="summary", stop_on_failure=TRUE)'
exec "$RSCRIPT_BIN" --vanilla experiments/system_a_coarse_guide.R \
  "$SAB_PROJECT_ROOT/outputs/coarse_guide_${SLURM_JOB_ID}"
