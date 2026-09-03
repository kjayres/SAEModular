#!/usr/bin/env bash
# Convert and validate one completed production SAEM run against the full
# sealed System A VODE-BDF target. This job does not launch patient banks.
# Required: SAB_SAEM_RUN_DIR. Optional: SAB_ANCHOR_OUTPUT_DIR.
#SBATCH --job-name=sab-sysa-prep
#SBATCH --mail-user=reub0840@stats.ox.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --cluster=srf_cpu_01
#SBATCH --partition=standard-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=02:00:00
#SBATCH --chdir=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes
#SBATCH --output=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/system-a-prepare-anchor-%j.out
#SBATCH --error=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/system-a-prepare-anchor-%j.err

set -Eeuo pipefail
umask 027

WORKSPACE_ROOT="/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0"
PROJECT_ROOT="$WORKSPACE_ROOT/projects/saem_anchored_bayes"
RSCRIPT_BIN="$WORKSPACE_ROOT/software/miniforge3/envs/model_fusion_r/bin/Rscript"
SAEM_RUN_DIR="${SAB_SAEM_RUN_DIR:?Set SAB_SAEM_RUN_DIR to a completed production SAEM run}"
OUTPUT_DIR="${SAB_ANCHOR_OUTPUT_DIR:-$PROJECT_ROOT/outputs/system_a_anchor_${SLURM_JOB_ID}}"

[[ -x "$RSCRIPT_BIN" ]] || {
  printf 'Rscript is not executable: %s\n' "$RSCRIPT_BIN" >&2
  exit 2
}
[[ -d "$SAEM_RUN_DIR" ]] || {
  printf 'Completed SAEM run directory is absent: %s\n' "$SAEM_RUN_DIR" >&2
  exit 2
}
[[ ! -e "$OUTPUT_DIR" ]] || {
  printf 'Refusing to overwrite anchor output: %s\n' "$OUTPUT_DIR" >&2
  exit 2
}
[[ -f "$PROJECT_ROOT/R/system_a_adapter.R" ]] || {
  printf 'System A adapter module is absent.\n' >&2
  exit 2
}
[[ -f "$PROJECT_ROOT/R/system_a_saem_anchor.R" ]] || {
  printf 'System A anchor-conversion module is absent.\n' >&2
  exit 2
}

export SAB_WORKSPACE_ROOT="$WORKSPACE_ROOT"
export SAB_PROJECT_ROOT="$PROJECT_ROOT"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

exec "$RSCRIPT_BIN" --vanilla \
  experiments/system_a_prepare_anchor.R "$SAEM_RUN_DIR" "$OUTPUT_DIR"
