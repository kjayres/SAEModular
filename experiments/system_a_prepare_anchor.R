#!/usr/bin/env Rscript

# Convert one completed production SAEM export into the canonical System A
# anchor used by the message experiment.  This is a fail-closed conversion and
# validation step: it performs no fitting and never substitutes SAEM's LSODA
# likelihood for the sealed VODE-BDF target.

sab_prepare_stop <- function(...) {
  stop(paste0(...), call. = FALSE)
}

sab_prepare_env <- function(name) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) sab_prepare_stop("Missing environment variable ", name, ".")
  value
}

sab_prepare_sha256 <- function(path) {
  binary <- Sys.which("sha256sum")
  if (!nzchar(binary)) sab_prepare_stop("sha256sum is required.")
  path <- normalizePath(path, mustWork = TRUE)
  output <- suppressWarnings(system2(
    binary, c("--", shQuote(path)), stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  if ((!is.null(status) && status != 0L) || length(output) != 1L) {
    sab_prepare_stop("Could not calculate SHA-256 for ", path, ".")
  }
  digest <- tolower(strsplit(trimws(output), "[[:space:]]+")[[1L]][[1L]])
  if (!grepl("^[0-9a-f]{64}$", digest)) {
    sab_prepare_stop("sha256sum returned an invalid digest for ", path, ".")
  }
  digest
}

sab_prepare_atomic_file <- function(path, writer) {
  if (file.exists(path)) sab_prepare_stop("Refusing to overwrite ", path, ".")
  partial <- paste0(path, ".partial")
  if (file.exists(partial)) {
    sab_prepare_stop("Refusing to overwrite stale staging file ", partial, ".")
  }
  writer(partial)
  if (!file.exists(partial) || is.na(file.info(partial)$size) ||
      file.info(partial)$size <= 0) {
    sab_prepare_stop("Writer did not create a non-empty file for ", path, ".")
  }
  if (!file.rename(partial, path)) {
    sab_prepare_stop("Atomic file publication failed for ", path, ".")
  }
  invisible(path)
}

sab_prepare_validate_production <- function(anchor) {
  source <- anchor$source
  audit <- if (is.list(source)) source$audit else NULL
  run_config <- if (is.list(source)) source$run_config else NULL
  fit <- if (is.list(run_config)) run_config$fit_control else NULL
  valid <- is.list(source) && is.list(audit) &&
    identical(audit$status, "completed") &&
    identical(audit$n_patients, 115L) &&
    is.list(fit) && is.integer(fit$k1) && length(fit$k1) == 1L &&
    is.integer(fit$k2) && length(fit$k2) == 1L &&
    fit$k1 >= 300L && fit$k2 >= 100L &&
    isTRUE(fit$map)
  if (!isTRUE(valid)) {
    sab_prepare_stop(
      "Production anchor preparation requires a completed 115-patient SAEM ",
      "run with MAP output, K1 >= 300, and K2 >= 100."
    )
  }
  invisible(TRUE)
}

sab_prepare_validate_exact <- function(anchor) {
  exact <- anchor$validation$exact_target
  summary <- if (is.list(exact)) exact$summary else NULL
  ledger <- if (is.list(exact)) exact$ledger else NULL
  valid <- inherits(exact, "sab_system_a_saem_exact_validation") &&
    is.data.frame(summary) && nrow(summary) == 1L &&
    is.data.frame(ledger) && nrow(ledger) == 115L &&
    identical(as.character(ledger$patient_id), anchor$patient_ids) &&
    all(ledger$prediction_calls == 1L) &&
    all(ledger$likelihood_calls == 1L) &&
    all(!is.na(ledger$ode_integrations)) &&
    all(ledger$ode_integrations >= 0L) &&
    all(!ledger$callback_error) &&
    isTRUE(exact$completed_without_callback_error) &&
    summary$patients[[1L]] == 115L &&
    summary$prediction_calls[[1L]] == 115L &&
    summary$likelihood_calls[[1L]] == 115L &&
    summary$unknown_ode_counts[[1L]] == 0L &&
    summary$callback_errors[[1L]] == 0L &&
    isTRUE(summary$completed_without_callback_error[[1L]])
  if (!isTRUE(valid)) {
    sab_prepare_stop(
      "All 115 SAEM states were not cleanly evaluated by the sealed exact ",
      "System A callbacks."
    )
  }
  # An SAEM MAP state may be an ordinary rejection under the more stringent
  # exact target.  Such rows remain explicit as -Inf and are not callback
  # failures; downstream banks must start those patients from a valid fallback.
  invisible(TRUE)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 2L) {
  sab_prepare_stop(
    "Usage: system_a_prepare_anchor.R <completed-SAEM-run-dir> ",
    "<new-output-dir>."
  )
}
command_line <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_line, value = TRUE)
if (length(file_argument) != 1L) {
  sab_prepare_stop("Could not identify this preparation script's path.")
}

workspace_root <- normalizePath(
  sab_prepare_env("SAB_WORKSPACE_ROOT"), mustWork = TRUE
)
project_root <- normalizePath(
  sab_prepare_env("SAB_PROJECT_ROOT"), mustWork = TRUE
)
run_directory <- normalizePath(arguments[[1L]], mustWork = TRUE)
output_argument <- arguments[[2L]]
if (!nzchar(output_argument) || basename(output_argument) %in% c(".", "..")) {
  sab_prepare_stop("The output directory must have a non-empty final name.")
}
output_parent <- dirname(output_argument)
if (!dir.exists(output_parent) &&
    !dir.create(output_parent, recursive = TRUE, showWarnings = FALSE)) {
  sab_prepare_stop("Could not create output parent directory ", output_parent, ".")
}
output_directory <- file.path(
  normalizePath(output_parent, mustWork = TRUE), basename(output_argument)
)
if (file.exists(output_directory) || dir.exists(output_directory)) {
  sab_prepare_stop("Refusing to overwrite output ", output_directory, ".")
}

adapter_path <- file.path(project_root, "R", "system_a_adapter.R")
converter_path <- file.path(project_root, "R", "system_a_saem_anchor.R")
script_path <- normalizePath(sub("^--file=", "", file_argument), mustWork = TRUE)
source_paths <- c(
  "R/system_a_adapter.R" = adapter_path,
  "R/system_a_saem_anchor.R" = converter_path,
  "experiments/system_a_prepare_anchor.R" = script_path
)
if (any(!file.exists(source_paths))) {
  sab_prepare_stop("Required anchor preparation source is absent.")
}
source_sha256_before <- vapply(source_paths, sab_prepare_sha256, character(1L))

sys.source(adapter_path, envir = globalenv(), keep.source = FALSE)
sys.source(converter_path, envir = globalenv(), keep.source = FALSE)
adapter <- sab_load_system_a_adapter(
  workspace_root = workspace_root, patient_ids = NULL
)
if (!identical(adapter$numerical_target, "sealed_deSolve_VODE_BDF") ||
    length(adapter$patient_ids) != 115L || anyDuplicated(adapter$patient_ids)) {
  sab_prepare_stop("The loaded adapter is not the full sealed System A target.")
}

message("Reading completed production SAEM export: ", run_directory)
message("Validating all 115 SAEM states through sealed VODE-BDF callbacks.")
anchor <- sab_read_system_a_saem_anchor(
  run_directory = run_directory,
  adapter = adapter,
  validate_exact = TRUE
)
sab_prepare_validate_production(anchor)
sab_prepare_validate_exact(anchor)
sab_validate_system_a_saem_anchor(anchor, adapter)

source_sha256_after <- vapply(source_paths, sab_prepare_sha256, character(1L))
if (!identical(source_sha256_after, source_sha256_before)) {
  sab_prepare_stop("Preparation source changed while the anchor was built.")
}
target_source_after <- sab_validate_system_a_sources(workspace_root)
target_source_fields <- c(
  "upstream_commit", "source_sha256", "canonical_source_sha256",
  "oracle_certificate_sha256", "target_fingerprint"
)
if (!identical(
    unclass(target_source_after[target_source_fields]),
    unclass(adapter$source_validation[target_source_fields])
  )) {
  sab_prepare_stop("The sealed System A source changed during validation.")
}

stage_directory <- file.path(
  dirname(output_directory),
  paste0(".", basename(output_directory), ".job_",
         Sys.getenv("SLURM_JOB_ID", unset = "unknown"), ".tmp")
)
if (file.exists(stage_directory) || dir.exists(stage_directory)) {
  sab_prepare_stop("Refusing to overwrite stale staging directory ",
                   stage_directory, ".")
}
if (!dir.create(stage_directory, recursive = FALSE, showWarnings = FALSE)) {
  sab_prepare_stop("Could not create staging directory ", stage_directory, ".")
}
published <- FALSE
on.exit({
  if (!published && dir.exists(stage_directory)) {
    unlink(stage_directory, recursive = TRUE, force = FALSE)
  }
}, add = TRUE)

anchor_path <- file.path(stage_directory, "canonical_anchor.rds")
ledger_path <- file.path(stage_directory, "exact_validation_ledger.csv")
summary_path <- file.path(stage_directory, "exact_validation_summary.csv")
manifest_path <- file.path(stage_directory, "manifest.json")
sab_prepare_atomic_file(anchor_path, function(path) {
  saveRDS(anchor, path, version = 3L, compress = "gzip")
})
sab_prepare_atomic_file(ledger_path, function(path) {
  utils::write.csv(
    anchor$validation$exact_target$ledger, path,
    row.names = FALSE, quote = TRUE, na = "NA"
  )
})
sab_prepare_atomic_file(summary_path, function(path) {
  utils::write.csv(
    anchor$validation$exact_target$summary, path,
    row.names = FALSE, quote = TRUE, na = "NA"
  )
})

output_paths <- c(
  canonical_anchor = anchor_path,
  exact_validation_ledger = ledger_path,
  exact_validation_summary = summary_path
)
output_sha256 <- vapply(output_paths, sab_prepare_sha256, character(1L))
manifest <- list(
  schema_version = "sab_system_a_prepared_anchor_manifest_v1",
  purpose = "SAEM-anchored fixed-psi patient-message experiment",
  production_gate = list(
    passed = TRUE,
    minimum_k1 = 300L,
    minimum_k2 = 100L,
    observed_k1 = anchor$source$run_config$fit_control$k1,
    observed_k2 = anchor$source$run_config$fit_control$k2,
    saem_status = anchor$source$audit$status
  ),
  saem = list(
    run_directory = run_directory,
    run_tag = anchor$source$run_config$run_tag,
    patients = length(anchor$patient_ids),
    artifact_sha256 = as.list(anchor$source$artifact_sha256)
  ),
  exact_target = list(
    numerical_target = adapter$numerical_target,
    target_fingerprint = anchor$target_fingerprint,
    likelihood_signature = adapter$likelihood_signature,
    oracle_audit = adapter$oracle_audit,
    source_validation = unclass(adapter$source_validation),
    validation_schema = anchor$validation$exact_target$schema_version,
    validation = as.list(anchor$validation$exact_target$summary[1L, ])
  ),
  anchor_coordinates = list(
    eta = as.list(anchor$eta),
    psi = as.list(anchor$psi)
  ),
  inputs = list(
    source_paths = as.list(source_paths),
    source_sha256 = as.list(source_sha256_after)
  ),
  outputs = list(
    files = as.list(stats::setNames(basename(output_paths), names(output_paths))),
    sha256 = as.list(output_sha256)
  ),
  execution = list(
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = NA_character_),
    hostname = Sys.info()[["nodename"]],
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    r_version = R.version.string,
    session_info = utils::capture.output(utils::sessionInfo())
  )
)
sab_prepare_atomic_file(manifest_path, function(path) {
  jsonlite::write_json(
    manifest, path, auto_unbox = TRUE, pretty = TRUE,
    null = "null", na = "null", digits = NA
  )
})

# Recheck every long-lived input immediately before publishing the directory.
if (!identical(
    vapply(source_paths, sab_prepare_sha256, character(1L)),
    source_sha256_before
  ) || !identical(
    vapply(anchor$source$artifact_files, sab_prepare_sha256, character(1L)),
    anchor$source$artifact_sha256
  )) {
  sab_prepare_stop("Source or SAEM artifacts changed before publication.")
}
if (!file.rename(stage_directory, output_directory)) {
  sab_prepare_stop("Atomic anchor-directory publication failed for ",
                   output_directory, ".")
}
published <- TRUE

message("Published canonical production anchor: ",
        file.path(output_directory, "canonical_anchor.rds"))
message(
  "Exact validation: ",
  anchor$validation$exact_target$summary$patients[[1L]], " patients, ",
  anchor$validation$exact_target$summary$callback_errors[[1L]],
  " callback errors, ",
  anchor$validation$exact_target$summary$ordinary_target_rejections[[1L]],
  " ordinary exact-target rejections."
)
