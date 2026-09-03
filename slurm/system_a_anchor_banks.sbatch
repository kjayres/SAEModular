#!/usr/bin/env bash
# Exact-target fixed-psi conditional MCMC banks for the twelve predeclared
# System A patients under h_i = 0.5*g_SAEM,i + 0.5*g_priorcentral,i. Each array
# task owns one patient and runs independent mixture-proposal chains
# sequentially on one CPU. This launcher does not run message calculations.
# Required: SAB_ANCHOR_RDS. Optional controls are SAB_BANK_OUTPUT_DIR,
# SAB_BANK_CHAINS, SAB_BANK_WARMUP, SAB_BANK_DRAWS, SAB_BANK_THIN,
# SAB_BANK_INITIAL_BETA, SAB_BANK_TARGET_ACCEPTANCE,
# SAB_BANK_ADAPTATION_BLOCK, SAB_BANK_BASE_SEED, and
# SAB_BANK_START_OFFSET_SD.
#SBATCH --job-name=sab-sysa-bank
#SBATCH --mail-user=reub0840@stats.ox.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --cluster=srf_cpu_01
#SBATCH --partition=standard-cpu
#SBATCH --array=0-11%12
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --time=24:00:00
#SBATCH --chdir=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes
#SBATCH --output=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/system-a-bank-%A_%a.out
#SBATCH --error=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/system-a-bank-%A_%a.err

set -Eeuo pipefail
umask 027

WORKSPACE_ROOT="/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0"
PROJECT_ROOT="$WORKSPACE_ROOT/projects/saem_anchored_bayes"
RSCRIPT_BIN="$WORKSPACE_ROOT/software/miniforge3/envs/model_fusion_r/bin/Rscript"
ANCHOR_RDS="${SAB_ANCHOR_RDS:?Set SAB_ANCHOR_RDS to a canonical anchor RDS}"
OUTPUT_ROOT="${SAB_BANK_OUTPUT_DIR:-$PROJECT_ROOT/outputs/system_a_anchor_banks_${SLURM_ARRAY_JOB_ID}}"

[[ -x "$RSCRIPT_BIN" ]] || {
  printf 'Rscript is not executable: %s\n' "$RSCRIPT_BIN" >&2
  exit 2
}
[[ -f "$ANCHOR_RDS" ]] || {
  printf 'Canonical SAEM anchor is absent: %s\n' "$ANCHOR_RDS" >&2
  exit 2
}
[[ -f "$PROJECT_ROOT/R/system_a_adapter.R" ]] || {
  printf 'System A adapter module is absent.\n' >&2
  exit 2
}
[[ -f "$PROJECT_ROOT/R/system_a_patient_bank.R" ]] || {
  printf 'System A patient-bank module is absent.\n' >&2
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
  "$ANCHOR_RDS" "$OUTPUT_ROOT" "${SLURM_ARRAY_TASK_ID:?Missing array index}"
