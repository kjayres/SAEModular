#!/usr/bin/env Rscript

# Execute the hash-pinned upstream SAEM driver with only four audited changes:
# load the omega overlay after the upstream model, keep outputs in this project,
# enable saemix warnings, and verify/save the iteration-0 allpar row.

project_root <- Sys.getenv("SAEMODULAR_PROJECT_ROOT", unset = "")
upstream_root <- Sys.getenv("SAEM_PROJECT_ROOT", unset = "")
if (!nzchar(project_root) || !dir.exists(project_root)) {
  stop("SAEMODULAR_PROJECT_ROOT must name the SAEModular checkout.")
}
if (!nzchar(upstream_root) || !dir.exists(upstream_root)) {
  stop("SAEM_PROJECT_ROOT must name the pinned upstream SAEM checkout.")
}

overlay_path <- file.path(project_root, "R", "system_a_saem_omega_overlay.R")
runner_path <- file.path(upstream_root, "scripts", "run_hiv_saem_fit.R")
if (!file.exists(overlay_path) || !file.exists(runner_path)) {
  stop("The correction overlay or pinned upstream runner is missing.")
}

source(overlay_path, local = FALSE)
upstream_lines <- readLines(runner_path, warn = TRUE)
effective_lines <- sab_build_corrected_saem_runner(upstream_lines)

# Parse before starting any work.  Exact insertion counts are checked by the
# transformer, so upstream source drift fails instead of changing semantics.
effective_expression <- parse(text = effective_lines, keep.source = TRUE)
run_tag <- Sys.getenv("SAEM_RUN_TAG", unset = "")
output_root <- Sys.getenv("SAEM_OUTPUT_ROOT", unset = "")
if (!nzchar(run_tag) || !nzchar(output_root) || !dir.exists(output_root)) {
  stop("SAEM_RUN_TAG and an existing SAEM_OUTPUT_ROOT are required.")
}
run_directory <- file.path(output_root, run_tag)
if (dir.exists(run_directory) || !dir.create(run_directory)) {
  stop("Refusing to overwrite or unable to create: ", run_directory)
}
writeLines(
  effective_lines,
  file.path(run_directory, "effective_corrected_saem_runner.R"),
  useBytes = TRUE
)
eval(effective_expression, envir = .GlobalEnv)
