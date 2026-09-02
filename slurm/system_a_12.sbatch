#!/usr/bin/env bash
# ORACLE mechanism pilot only: no retained posterior chain and no bridge step.
# All ODE work runs on the allocated compute node and is parallelised over the
# twelve predeclared patients.
#SBATCH --job-name=sab-sysa12
#SBATCH --mail-user=reub0840@stats.ox.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --cluster=srf_cpu_01
#SBATCH --partition=standard-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --time=12:00:00
#SBATCH --chdir=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes
#SBATCH --output=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/system-a-12-%j.out
#SBATCH --error=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/system-a-12-%j.err

set -Eeuo pipefail
umask 027

WORKSPACE_ROOT="/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0"
PROJECT_ROOT="$WORKSPACE_ROOT/projects/saem_anchored_bayes"
RSCRIPT_BIN="$WORKSPACE_ROOT/software/miniforge3/envs/model_fusion_r/bin/Rscript"
OUTPUT_DIR="${SAB_SYSTEM_A_OUTPUT_DIR:-$PROJECT_ROOT/outputs/system_a_12_oracle_v1}"

[[ -x "$RSCRIPT_BIN" ]] || {
  printf 'Rscript is not executable: %s\n' "$RSCRIPT_BIN" >&2
  exit 2
}
[[ -d "$WORKSPACE_ROOT/projects/.worktrees/modular_bayes_system_a_factor_v1" ]] || {
  printf 'Pinned System A worktree is absent.\n' >&2
  exit 2
}

export SAB_WORKSPACE_ROOT="$WORKSPACE_ROOT"
export SAB_SYSTEM_A_CORES="${SLURM_CPUS_PER_TASK:?Missing Slurm CPU allocation}"
export SAB_SYSTEM_A_INTEGRATION=true
export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

mkdir -p "$OUTPUT_DIR"

"$RSCRIPT_BIN" --vanilla -e \
  'testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)'

exec "$RSCRIPT_BIN" --vanilla experiments/system_a_12.R "$OUTPUT_DIR"
