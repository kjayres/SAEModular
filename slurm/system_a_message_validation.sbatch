#!/usr/bin/env bash
# ODE-free planning or final raw-message/bridge assessment. Required variables:
# SAB_MESSAGE_MODE=plan|assess, SAB_ANCHOR_RDS, SAB_ANCHOR_BANK_DIR,
# SAB_MESSAGE_OUTPUT_DIR. Assess mode additionally requires SAB_ENDPOINT_PLAN
# and SAB_CANDIDATE_BANK_DIR.
#SBATCH --job-name=sab-sysa-message
#SBATCH --mail-user=reub0840@stats.ox.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --cluster=srf_cpu_01
#SBATCH --partition=standard-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=4G
#SBATCH --time=00:30:00
#SBATCH --chdir=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes
#SBATCH --output=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/system-a-message-%j.out
#SBATCH --error=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/system-a-message-%j.err

set -Eeuo pipefail
umask 027

WORKSPACE_ROOT="/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0"
PROJECT_ROOT="$WORKSPACE_ROOT/projects/saem_anchored_bayes"
RSCRIPT_BIN="$WORKSPACE_ROOT/software/miniforge3/envs/model_fusion_r/bin/Rscript"
MODE="${SAB_MESSAGE_MODE:?Set SAB_MESSAGE_MODE to plan or assess}"
ANCHOR_RDS="${SAB_ANCHOR_RDS:?Set SAB_ANCHOR_RDS}"
ANCHOR_BANK_DIR="${SAB_ANCHOR_BANK_DIR:?Set SAB_ANCHOR_BANK_DIR}"
OUTPUT_DIR="${SAB_MESSAGE_OUTPUT_DIR:?Set SAB_MESSAGE_OUTPUT_DIR}"

[[ "$MODE" == "plan" || "$MODE" == "assess" ]] || {
  printf 'SAB_MESSAGE_MODE must be plan or assess.\n' >&2
  exit 2
}
[[ -x "$RSCRIPT_BIN" && -f "$ANCHOR_RDS" && -d "$ANCHOR_BANK_DIR" ]] || {
  printf 'Rscript, anchor, or anchor-bank input is missing.\n' >&2
  exit 2
}

export SAB_WORKSPACE_ROOT="$WORKSPACE_ROOT"
export SAB_PROJECT_ROOT="$PROJECT_ROOT"
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1

if [[ "$MODE" == "plan" ]]; then
  exec "$RSCRIPT_BIN" --vanilla \
    experiments/system_a_message_validation.R \
    plan "$ANCHOR_RDS" "$ANCHOR_BANK_DIR" "$OUTPUT_DIR"
fi

ENDPOINT_PLAN="${SAB_ENDPOINT_PLAN:?Set SAB_ENDPOINT_PLAN for assess mode}"
CANDIDATE_BANK_DIR="${SAB_CANDIDATE_BANK_DIR:?Set SAB_CANDIDATE_BANK_DIR for assess mode}"
[[ -f "$ENDPOINT_PLAN" && -d "$CANDIDATE_BANK_DIR" ]] || {
  printf 'Endpoint plan or candidate-bank directory is missing.\n' >&2
  exit 2
}
exec "$RSCRIPT_BIN" --vanilla \
  experiments/system_a_message_validation.R \
  assess "$ANCHOR_RDS" "$ANCHOR_BANK_DIR" "$ENDPOINT_PLAN" \
  "$CANDIDATE_BANK_DIR" "$OUTPUT_DIR"
