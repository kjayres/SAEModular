#!/usr/bin/env bash
# Final bounded System A experiment: pure-SAEM fixed-psi reference banks for
# the twelve predeclared patients.  Each task runs four sequential pCN chains.
# The scientific budget is frozen here and must not be extended after seeing
# diagnostics. Required: SAB_ANCHOR_RDS, SAB_BANK_OUTPUT_DIR,
# SAB_BANK_BASE_SEED, SAB_EXPECTED_ANCHOR_SHA256.
#SBATCH --job-name=sab-sysa-pure-ref
#SBATCH --mail-user=reub0840@stats.ox.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --cluster=srf_cpu_01
#SBATCH --partition=standard-cpu
#SBATCH --array=0-11%12
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=01:00:00
#SBATCH --chdir=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes
#SBATCH --output=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/system-a-pure-ref-%A_%a.out
#SBATCH --error=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/system-a-pure-ref-%A_%a.err

set -Eeuo pipefail
umask 027

WORKSPACE_ROOT="/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0"
PROJECT_ROOT="$WORKSPACE_ROOT/projects/saem_anchored_bayes"
RSCRIPT_BIN="$WORKSPACE_ROOT/software/miniforge3/envs/model_fusion_r/bin/Rscript"
ANCHOR_RDS="${SAB_ANCHOR_RDS:?Set SAB_ANCHOR_RDS to a corrected canonical anchor}"
OUTPUT_ROOT="${SAB_BANK_OUTPUT_DIR:?Set SAB_BANK_OUTPUT_DIR to a new branch-specific directory}"
BASE_SEED="${SAB_BANK_BASE_SEED:?Set a predeclared branch-specific base seed}"
EXPECTED_ANCHOR_SHA256="${SAB_EXPECTED_ANCHOR_SHA256:?Pin the corrected canonical anchor SHA-256}"

[[ -x "$RSCRIPT_BIN" && -f "$ANCHOR_RDS" ]] || {
  printf 'Rscript or canonical anchor is missing.\n' >&2
  exit 2
}
OBSERVED_ANCHOR_SHA256="$(sha256sum -- "$ANCHOR_RDS" | awk '{print $1}')"
[[ "$OBSERVED_ANCHOR_SHA256" == "$EXPECTED_ANCHOR_SHA256" ]] || {
  printf 'Canonical anchor hash does not match the predeclared branch.\n' >&2
  exit 2
}
[[ ! -e "$OUTPUT_ROOT" || -d "$OUTPUT_ROOT" ]] || {
  printf 'Bank output root exists but is not a directory: %s\n' "$OUTPUT_ROOT" >&2
  exit 2
}
mkdir -p "$OUTPUT_ROOT"

export SAB_WORKSPACE_ROOT="$WORKSPACE_ROOT"
export SAB_PROJECT_ROOT="$PROJECT_ROOT"
export SAB_BANK_CHAINS=4
export SAB_BANK_WARMUP=250
export SAB_BANK_DRAWS=500
export SAB_BANK_THIN=1
export SAB_BANK_INITIAL_BETA=0.4
export SAB_BANK_TARGET_ACCEPTANCE=0.3
export SAB_BANK_ADAPTATION_BLOCK=50
export SAB_BANK_BASE_SEED="$BASE_SEED"
export SAB_BANK_START_OFFSET_SD=1
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

exec "$RSCRIPT_BIN" --vanilla \
  experiments/system_a_anchor_bank_patient.R \
  "$ANCHOR_RDS" "$OUTPUT_ROOT" "${SLURM_ARRAY_TASK_ID:?Missing array index}" \
  pure_saem_reference
