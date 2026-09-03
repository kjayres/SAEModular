#!/usr/bin/env bash
# Exact-target fixed-psi conditional MCMC banks for the two components of the
# defensive reference: g_SAEM and g_priorcentral.  Two endpoints times twelve
# patient targets give 24 array tasks.  This launcher must not be used until
# the canonical SAEM anchor, anchor banks, and endpoint plan are complete.
# Required: SAB_ANCHOR_RDS and SAB_ENDPOINT_PLAN. Optional bank controls match
# slurm/system_a_anchor_banks.sbatch.
#SBATCH --job-name=sab-sysa-candidate-bank
#SBATCH --mail-user=reub0840@stats.ox.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --cluster=srf_cpu_01
#SBATCH --partition=standard-cpu
#SBATCH --array=0-23%12
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --chdir=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes
#SBATCH --output=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/system-a-candidate-bank-%A_%a.out
#SBATCH --error=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/system-a-candidate-bank-%A_%a.err

set -Eeuo pipefail
umask 027

WORKSPACE_ROOT="/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0"
PROJECT_ROOT="$WORKSPACE_ROOT/projects/saem_anchored_bayes"
RSCRIPT_BIN="$WORKSPACE_ROOT/software/miniforge3/envs/model_fusion_r/bin/Rscript"
ANCHOR_RDS="${SAB_ANCHOR_RDS:?Set SAB_ANCHOR_RDS to the canonical anchor RDS}"
ENDPOINT_PLAN="${SAB_ENDPOINT_PLAN:?Set SAB_ENDPOINT_PLAN to endpoint_plan.rds}"
OUTPUT_ROOT="${SAB_BANK_OUTPUT_DIR:-$PROJECT_ROOT/outputs/system_a_candidate_banks_${SLURM_ARRAY_JOB_ID}}"

[[ -x "$RSCRIPT_BIN" ]] || {
  printf 'Rscript is not executable: %s\n' "$RSCRIPT_BIN" >&2
  exit 2
}
[[ -f "$ANCHOR_RDS" ]] || {
  printf 'Canonical SAEM anchor is absent: %s\n' "$ANCHOR_RDS" >&2
  exit 2
}
[[ -f "$ENDPOINT_PLAN" ]] || {
  printf 'Endpoint plan is absent: %s\n' "$ENDPOINT_PLAN" >&2
  exit 2
}

mkdir -p "$OUTPUT_ROOT"

export SAB_WORKSPACE_ROOT="$WORKSPACE_ROOT"
export SAB_PROJECT_ROOT="$PROJECT_ROOT"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

exec "$RSCRIPT_BIN" --vanilla \
  experiments/system_a_anchor_bank_patient.R \
  "$ANCHOR_RDS" "$OUTPUT_ROOT" "${SLURM_ARRAY_TASK_ID:?Missing array index}" \
  "$ENDPOINT_PLAN"
