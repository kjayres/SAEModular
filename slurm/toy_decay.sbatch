#!/usr/bin/env bash
#SBATCH --job-name=sab-toy-decay
#SBATCH --mail-user=reub0840@stats.ox.ac.uk
#SBATCH --mail-type=FAIL
#SBATCH --cluster=srf_cpu_01
#SBATCH --partition=standard-cpu
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=2G
#SBATCH --time=00:30:00
#SBATCH --chdir=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes
#SBATCH --output=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/toy-decay-%j.out
#SBATCH --error=/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0/projects/saem_anchored_bayes/logs/toy-decay-%j.err

set -Eeuo pipefail
umask 027

WORKSPACE_ROOT="/ceph/volumes/homes_rstudent/reub0840/f341ced6-637f-43b2-a7a8-b9e1a63a82f0"
PROJECT_ROOT="$WORKSPACE_ROOT/projects/saem_anchored_bayes"
RSCRIPT_BIN="$WORKSPACE_ROOT/software/miniforge3/envs/model_fusion_r/bin/Rscript"
OUTPUT_DIR="${SAB_TOY_OUTPUT_DIR:-$PROJECT_ROOT/outputs/toy_decay_v1}"

[[ -x "$RSCRIPT_BIN" ]] || {
  printf 'Rscript is not executable: %s\n' "$RSCRIPT_BIN" >&2
  exit 2
}

export OMP_NUM_THREADS=1
export OPENBLAS_NUM_THREADS=1
export MKL_NUM_THREADS=1
export VECLIB_MAXIMUM_THREADS=1
export NUMEXPR_NUM_THREADS=1

mkdir -p "$OUTPUT_DIR"

"$RSCRIPT_BIN" --vanilla -e \
  'testthat::test_dir("tests/testthat", reporter = "summary", stop_on_failure = TRUE)'

exec "$RSCRIPT_BIN" --vanilla experiments/toy_decay.R "$OUTPUT_DIR"
