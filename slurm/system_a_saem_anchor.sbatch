#!/usr/bin/env bash
# Anchor-only all-patient SAEM fit. The subsequent patient-message banks use
# the separately pinned, exact System A VODE-BDF target.
#SBATCH --job-name=sab-sysa-saem
#SBATCH --mail-user=reub0840@stats.ox.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --cluster=srf_cpu_01
#SBATCH --partition=standard-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --chdir=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes
#SBATCH --output=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/system-a-saem-%j.out
#SBATCH --error=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/system-a-saem-%j.err

set -Eeuo pipefail
umask 027

WORKSPACE_ROOT="/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0"
PROJECT_ROOT="$WORKSPACE_ROOT/projects/saem_anchored_bayes"
SAEM_SOURCE_ROOT="$WORKSPACE_ROOT/projects/saem_model"
RSCRIPT_BIN="$WORKSPACE_ROOT/software/miniforge3/envs/model_fusion_r/bin/Rscript"
RUN_TAG="${SAB_SAEM_RUN_TAG:-saemodular_system_a_anchor_${SLURM_JOB_ID}}"
RUN_DIR="$SAEM_SOURCE_ROOT/outputs/$RUN_TAG"

[[ -x "$RSCRIPT_BIN" ]] || {
  printf 'Rscript is not executable: %s\n' "$RSCRIPT_BIN" >&2
  exit 2
}
[[ ! -e "$RUN_DIR" ]] || {
  printf 'Refusing to overwrite an existing SAEM output: %s\n' "$RUN_DIR" >&2
  exit 2
}

sha256sum --check "$PROJECT_ROOT/config/system_a_saem_sources.sha256"

export SAEM_PROJECT_ROOT="$SAEM_SOURCE_ROOT"
export SAEM_WORKSPACE_ROOT="$WORKSPACE_ROOT"
export SAEM_RUN_TAG="$RUN_TAG"
export HIV_BENCHMARK_DATA="$SAEM_SOURCE_ROOT/data/hiv_benchmark_canonical.csv"
export SAEM_SEED="${SAB_SAEM_SEED:-8}"
export SAEM_K1="${SAB_SAEM_K1:-300}"
export SAEM_K2="${SAB_SAEM_K2:-100}"
export SAEM_NB_MCMC="${SAB_SAEM_NB_MCMC:-2,2,2,0}"
export SAEM_NB_SA="${SAB_SAEM_NB_SA:-150}"
export SAEM_NB_BURN="${SAB_SAEM_NB_BURN:-5}"
export SAEM_NB_MAP="${SAB_SAEM_NB_MAP:-5}"
export SAEM_ENABLE_MAP="${SAB_SAEM_ENABLE_MAP:-true}"
export SAEM_ENABLE_FIM="${SAB_SAEM_ENABLE_FIM:-false}"
export SAEM_ENABLE_LLIS="${SAB_SAEM_ENABLE_LLIS:-false}"

# This is deliberately the established SAEM likelihood configuration. Its
# LSODA result is an anchor, never the exact retained target. Numerical target
# evaluation and every conditional bank remain pinned to VODE-BDF.
export SAEM_ODE_METHOD="lsoda"
export SAEM_ODE_RTOL="1e-8"
export SAEM_ODE_ATOL="1e-10"
export SAEM_ODE_MAXSTEPS="100000"
export SAEM_TIME_EPS="1e-6"

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

printf 'SAEM anchor output: %s\n' "$RUN_DIR"
exec "$RSCRIPT_BIN" --vanilla "$SAEM_SOURCE_ROOT/scripts/run_hiv_saem_fit.R"
