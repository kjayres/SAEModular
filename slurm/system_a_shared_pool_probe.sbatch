#!/usr/bin/env bash
# Two bounded offline questions; no new sampler or SAEM fit.
#SBATCH --job-name=sab-shared-probe
#SBATCH --mail-user=reub0840@stats.ox.ac.uk
#SBATCH --mail-type=FAIL
#SBATCH --cluster=srf_cpu_01
#SBATCH --partition=standard-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=6G
#SBATCH --time=01:00:00
#SBATCH --chdir=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes
#SBATCH --output=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/shared-probe-%j.out
#SBATCH --error=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/shared-probe-%j.err
set -Eeuo pipefail
umask 027
WORKSPACE_ROOT="/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0"
export SAB_PROJECT_ROOT="$WORKSPACE_ROOT/projects/saem_anchored_bayes"
RSCRIPT_BIN="$WORKSPACE_ROOT/software/miniforge3/envs/model_fusion_r/bin/Rscript"
export SAB_SYSTEM_A_INTEGRATION=true
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
"$RSCRIPT_BIN" --vanilla -e 'invisible(lapply(c(list.files("R", "[.]R$", full.names=TRUE), list.files("experiments", "[.]R$", full.names=TRUE)), function(p) parse(file=p))); testthat::test_dir("tests/testthat", reporter="summary", stop_on_failure=TRUE)'
if [[ -n "${SAB_REPLAY_SOURCE:-}" ]]; then
  exec "$RSCRIPT_BIN" --vanilla experiments/system_a_shared_pool_probe.R \
    --replay-saved "$SAB_REPLAY_SOURCE" \
    "$SAB_PROJECT_ROOT/outputs/shared_pool_probe_${SLURM_JOB_ID}"
fi
exec "$RSCRIPT_BIN" --vanilla experiments/system_a_shared_pool_probe.R \
  "$SAB_PROJECT_ROOT/outputs/shared_pool_probe_${SLURM_JOB_ID}"
