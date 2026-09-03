#!/usr/bin/env bash
# Corrected all-patient System A SAEM anchor.  The upstream source is hash
# pinned; this job executes the minimal audited omega-initialisation overlay.
# Scientific fits use the original 300 + 100 SAEM schedule unless explicitly
# overridden for a smoke run.
#SBATCH --job-name=sab-sysa-saem2
#SBATCH --mail-user=reub0840@stats.ox.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --cluster=srf_cpu_01
#SBATCH --partition=standard-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --time=48:00:00
#SBATCH --chdir=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes
#SBATCH --output=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/system-a-saem-corrected-%j.out
#SBATCH --error=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/system-a-saem-corrected-%j.err

set -Eeuo pipefail
umask 027

WORKSPACE_ROOT="/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0"
PROJECT_ROOT="$WORKSPACE_ROOT/projects/saem_anchored_bayes"
SAEM_SOURCE_ROOT="$WORKSPACE_ROOT/projects/saem_model"
RSCRIPT_BIN="$WORKSPACE_ROOT/software/miniforge3/envs/model_fusion_r/bin/Rscript"
RUN_TAG="${SAB_SAEM_RUN_TAG:-saemodular_system_a_corrected_${SLURM_JOB_ID}}"
RUN_DIR="$PROJECT_ROOT/outputs/$RUN_TAG"

[[ -x "$RSCRIPT_BIN" ]] || {
  printf 'Rscript is not executable: %s\n' "$RSCRIPT_BIN" >&2
  exit 2
}
[[ ! -e "$RUN_DIR" ]] || {
  printf 'Refusing to overwrite an existing corrected SAEM output: %s\n' \
    "$RUN_DIR" >&2
  exit 2
}

sha256sum --check "$PROJECT_ROOT/config/system_a_saem_sources.sha256"
sha256sum \
  "$PROJECT_ROOT/R/system_a_saem_omega_overlay.R" \
  "$PROJECT_ROOT/experiments/system_a_saem_corrected.R" \
  "$PROJECT_ROOT/slurm/system_a_saem_corrected.sbatch"

export SAEMODULAR_PROJECT_ROOT="$PROJECT_ROOT"
export SAEM_PROJECT_ROOT="$SAEM_SOURCE_ROOT"
export SAEM_OUTPUT_ROOT="$PROJECT_ROOT/outputs"
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

# Anchor construction uses the established SAEM LSODA likelihood.  Exact
# message-bank evaluation remains the separately pinned VODE-BDF target.
export SAEM_ODE_METHOD="lsoda"
export SAEM_ODE_RTOL="1e-8"
export SAEM_ODE_ATOL="1e-10"
export SAEM_ODE_MAXSTEPS="100000"
export SAEM_TIME_EPS="1e-6"

export OPENBLAS_NUM_THREADS=1
export OMP_NUM_THREADS=1
export MKL_NUM_THREADS=1

printf 'Corrected SAEM anchor output: %s\n' "$RUN_DIR"
printf 'Corrected SAEM seed: %s; schedule: K1=%s K2=%s\n' \
  "$SAEM_SEED" "$SAEM_K1" "$SAEM_K2"
exec "$RSCRIPT_BIN" --vanilla \
  "$PROJECT_ROOT/experiments/system_a_saem_corrected.R"
