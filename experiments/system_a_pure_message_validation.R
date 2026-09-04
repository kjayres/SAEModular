#!/usr/bin/env Rscript

# ODE-free planning and assessment for the final pure-SAEM System A test.
# Patient ODE work is performed only by separately submitted bank workers.
# Both modes publish a complete diagnostic directory before returning nonzero
# when a gate fails, allowing Slurm afterok dependencies to enforce the stop.

sab_pure_cli_stop <- function(...) stop(paste0(...), call. = FALSE)

sab_pure_cli_sha256 <- function(path) {
  binary <- Sys.which("sha256sum")
  if (!nzchar(binary)) sab_pure_cli_stop("sha256sum is required.")
  output <- suppressWarnings(system2(
    binary, c("--", shQuote(normalizePath(path, mustWork = TRUE))),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  if ((!is.null(status) && status != 0L) || length(output) != 1L) {
    sab_pure_cli_stop("Could not hash ", path, ".")
  }
  digest <- tolower(strsplit(trimws(output), "[[:space:]]+")[[1L]][[1L]])
  if (!grepl("^[0-9a-f]{64}$", digest)) {
    sab_pure_cli_stop("Invalid SHA-256 output for ", path, ".")
  }
  digest
}

sab_pure_cli_load_reference_banks <- function(root, patient_ids) {
  paths <- stats::setNames(file.path(
    root, sprintf("patient_%03d", as.integer(patient_ids)),
    "patient_banks.rds"
  ), patient_ids)
  missing <- names(paths)[!file.exists(paths)]
  if (length(missing)) {
    sab_pure_cli_stop("Missing pure reference banks for patients: ",
                      paste(missing, collapse = ", "), ".")
  }
  paths <- vapply(paths, normalizePath, character(1L), mustWork = TRUE)
  list(
    paths = paths,
    sha256 = vapply(paths, sab_pure_cli_sha256, character(1L)),
    artifacts = lapply(paths, readRDS)
  )
}

sab_pure_cli_load_candidate_banks <- function(root, endpoints, patient_ids) {
  artifacts <- stats::setNames(vector("list", nrow(endpoints)),
                               endpoints$endpoint_id)
  paths <- hashes <- artifacts
  for (endpoint_id in endpoints$endpoint_id) {
    endpoint_paths <- stats::setNames(file.path(
      root, endpoint_id, sprintf("patient_%03d", as.integer(patient_ids)),
      "patient_banks.rds"
    ), patient_ids)
    missing <- names(endpoint_paths)[!file.exists(endpoint_paths)]
    if (length(missing)) {
      sab_pure_cli_stop("Missing candidate banks for ", endpoint_id,
                        ", patients: ", paste(missing, collapse = ", "), ".")
    }
    endpoint_paths <- vapply(
      endpoint_paths, normalizePath, character(1L), mustWork = TRUE
    )
    paths[[endpoint_id]] <- endpoint_paths
    hashes[[endpoint_id]] <- vapply(
      endpoint_paths, sab_pure_cli_sha256, character(1L)
    )
    artifacts[[endpoint_id]] <- lapply(endpoint_paths, readRDS)
  }
  list(paths = paths, sha256 = hashes, artifacts = artifacts)
}

sab_pure_cli_publish <- function(output_directory, files, report) {
  output_directory <- path.expand(output_directory)
  if (file.exists(output_directory) || dir.exists(output_directory)) {
    sab_pure_cli_stop("Refusing to overwrite output directory ",
                      output_directory, ".")
  }
  parent <- dirname(output_directory)
  if (!dir.exists(parent) &&
      !dir.create(parent, recursive = TRUE, showWarnings = FALSE) &&
      !dir.exists(parent)) {
    sab_pure_cli_stop("Could not create output parent ", parent, ".")
  }
  stage <- paste0(output_directory, ".partial")
  if (file.exists(stage) || dir.exists(stage) ||
      !dir.create(stage, recursive = FALSE, mode = "0750")) {
    sab_pure_cli_stop("Could not create clean staging directory ", stage, ".")
  }
  published <- FALSE
  on.exit({
    if (!published && dir.exists(stage)) {
      unlink(stage, recursive = TRUE, force = FALSE)
    }
  }, add = TRUE)
  for (name in names(files)) {
    value <- files[[name]]
    path <- file.path(stage, name)
    if (grepl("\\.rds$", name)) {
      saveRDS(value, path, version = 3L, compress = "gzip")
    } else if (grepl("\\.csv$", name)) {
      utils::write.csv(
        value, path, row.names = FALSE, quote = TRUE, na = "NA"
      )
    } else {
      writeLines(as.character(value), path, useBytes = TRUE)
    }
  }
  writeLines(report, file.path(stage, "REPORT.txt"), useBytes = TRUE)
  manifest_paths <- list.files(stage, full.names = TRUE)
  manifest <- data.frame(
    file = basename(manifest_paths),
    sha256 = vapply(manifest_paths, sab_pure_cli_sha256, character(1L)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(
    manifest, file.path(stage, "manifest.csv"),
    row.names = FALSE, quote = TRUE
  )
  writeLines("completed", file.path(stage, "COMPLETED"), useBytes = TRUE)
  if (!file.rename(stage, output_directory)) {
    sab_pure_cli_stop("Atomic publication failed for ", output_directory, ".")
  }
  published <- TRUE
  normalizePath(output_directory, mustWork = TRUE)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (!length(arguments) %in% c(4L, 6L) ||
    !arguments[[1L]] %in% c("plan", "assess")) {
  sab_pure_cli_stop(
    "Usage: system_a_pure_message_validation.R plan <anchor.rds> ",
    "<pure-reference-bank-root> <output-dir> OR assess <anchor.rds> ",
    "<pure-reference-bank-root> <endpoint-plan.rds> ",
    "<candidate-bank-root> <output-dir>."
  )
}
mode <- arguments[[1L]]
if ((mode == "plan" && length(arguments) != 4L) ||
    (mode == "assess" && length(arguments) != 6L)) {
  sab_pure_cli_stop("The selected mode has the wrong number of arguments.")
}
workspace_root <- normalizePath(
  Sys.getenv("SAB_WORKSPACE_ROOT", unset = ""), mustWork = TRUE
)
project_root <- normalizePath(
  Sys.getenv("SAB_PROJECT_ROOT", unset = ""), mustWork = TRUE
)
anchor_path <- normalizePath(arguments[[2L]], mustWork = TRUE)
reference_bank_root <- normalizePath(arguments[[3L]], mustWork = TRUE)
modules <- file.path(project_root, "R", c(
  "system_a_adapter.R", "system_a_saem_anchor.R", "patient_messages.R",
  "system_a_message_validation.R", "system_a_pure_message_validation.R"
))
if (any(!file.exists(modules))) {
  sab_pure_cli_stop("Required pure-message R modules are absent.")
}
invisible(lapply(modules, sys.source, envir = globalenv(), keep.source = FALSE))

adapter <- sab_load_system_a_adapter(workspace_root)
canonical_anchor <- readRDS(anchor_path)
sab_validate_system_a_saem_anchor(canonical_anchor, adapter)
design <- sab_system_a_pure_message_design()
reference_banks <- sab_pure_cli_load_reference_banks(
  reference_bank_root, design$patient_ids
)
anchor_sha256 <- sab_pure_cli_sha256(anchor_path)
module_sha256 <- stats::setNames(
  vapply(modules, sab_pure_cli_sha256, character(1L)), basename(modules)
)

if (mode == "plan") {
  output_directory <- arguments[[4L]]
  result <- sab_system_a_plan_pure_message_endpoints(
    adapter = adapter,
    canonical_anchor = canonical_anchor,
    anchor_bank_artifacts = reference_banks$artifacts,
    anchor_path = anchor_path,
    anchor_sha256 = anchor_sha256,
    anchor_bank_paths = reference_banks$paths,
    anchor_bank_sha256 = reference_banks$sha256
  )
  result$plan$provenance <- list(
    module_sha256 = module_sha256,
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = NA_character_),
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
  )
  report <- c(
    "System A final pure-SAEM fixed-psi pilot",
    "",
    "The reference is exactly g_SAEM; no prior-centred mixture is used.",
    "Chains 01-02 alone decide this Stage-1 pilot gate. Chains 03-04 remain",
    "held out for assessment if every predeclared endpoint passes.",
    "Treatment endpoints move beta_nelf by one anchor omega_eta_pi SD.",
    "Scale endpoints multiply omega_lambda by 1/1.25 or 1.25.",
    "The diagnostic_fallback field is retained only for worker-schema",
    "compatibility: TRUE marks a failed endpoint; no fallback is run.",
    "The projected 115-patient MCSE is an extrapolation from the fixed",
    "12-patient panel, not a random-sample standard error or proof for 115.",
    "",
    paste("Pilot passed:", result$pilot_passed),
    "",
    capture.output(print(result$selection, row.names = FALSE))
  )
  published <- sab_pure_cli_publish(
    output_directory,
    list(
      "endpoint_plan.rds" = result$plan,
      "pilot_diagnostics.csv" = result$pilot_diagnostics,
      "pilot_mcmc_diagnostics.csv" = result$pilot_mcmc_diagnostics,
      "pilot_patient_mcse.csv" = result$pilot_patient_mcse,
      "selection.csv" = result$selection,
      "endpoints.csv" = result$plan$endpoints,
      "reference_bank_performance.csv" =
        result$reference_bank_performance
    ),
    report
  )
  message("Published pure-SAEM endpoint plan: ", published)
  if (!isTRUE(result$pilot_passed)) {
    message("Primary Stage-1 gate failed; candidate afterok jobs must not run.")
    quit(save = "no", status = 3L, runLast = FALSE)
  }
  quit(save = "no", status = 0L, runLast = FALSE)
}

plan_path <- normalizePath(arguments[[4L]], mustWork = TRUE)
candidate_bank_root <- normalizePath(arguments[[5L]], mustWork = TRUE)
output_directory <- arguments[[6L]]
plan <- readRDS(plan_path)
plan_sha256 <- sab_pure_cli_sha256(plan_path)
if (!identical(plan$anchor$path, anchor_path) ||
    !identical(plan$anchor$sha256, anchor_sha256) ||
    !identical(plan$anchor_banks$paths, reference_banks$paths) ||
    !identical(plan$anchor_banks$sha256, reference_banks$sha256) ||
    !identical(plan$provenance$module_sha256, module_sha256)) {
  sab_pure_cli_stop("Anchor or held-out banks differ from the pure pilot plan.")
}
if (!all(plan$endpoints$pilot_gate_passed)) {
  sab_pure_cli_stop("Assessment is forbidden because the pilot gate failed.")
}
candidate_banks <- sab_pure_cli_load_candidate_banks(
  candidate_bank_root, plan$endpoints, design$patient_ids
)
result <- sab_system_a_assess_pure_messages(
  adapter = adapter,
  canonical_anchor = canonical_anchor,
  plan = plan,
  plan_sha256 = plan_sha256,
  anchor_bank_artifacts = reference_banks$artifacts,
  candidate_bank_artifacts = candidate_banks$artifacts
)
result$provenance <- list(
  anchor = c(path = anchor_path, sha256 = anchor_sha256),
  endpoint_plan = c(path = plan_path, sha256 = plan_sha256),
  reference_bank_paths = reference_banks$paths,
  reference_bank_sha256 = reference_banks$sha256,
  candidate_bank_paths = candidate_banks$paths,
  candidate_bank_sha256 = candidate_banks$sha256,
  module_sha256 = module_sha256,
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = NA_character_),
  completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE)
)
report <- c(
  "System A final pure-SAEM 12-patient fixed-psi assessment",
  "",
  "Forward estimates use only held-out pure-reference chains 03-04.",
  "Reverse estimates use candidate chains 01-02. Bridge estimates use",
  "reference chains 01-02 (the Stage-1 pilot chains) and independent",
  "candidate chains 03-04; no chain used by either raw estimator enters",
  "the bridge on that same side.",
  "Every failed bridge or unavailable required MCMC diagnostic fails closed.",
  "The full-115 MCSE figures extrapolate the panel's 8 controls and 4",
  "treated patients separately to the full cohort's 80/35 composition.",
  "This is not a 115-patient posterior, dynamic-psi inference, exact frozen-",
  "bank posterior, or evidence about the other corrected SAEM branch.",
  "",
  paste("Status:", result$status),
  "",
  capture.output(print(result$endpoint_summary, row.names = FALSE))
)
published <- sab_pure_cli_publish(
  output_directory,
  list(
    "result.rds" = result,
    "endpoint_summary.csv" = result$endpoint_summary,
    "patient_summary.csv" = result$patient_summary,
    "directional_diagnostics.csv" = result$directional_diagnostics,
    "mcmc_diagnostics.csv" = result$mcmc_diagnostics,
    "bridge_chain_pairs.csv" = result$bridge_chain_pairs,
    "reference_bank_performance.csv" = result$reference_bank_performance,
    "cost_summary.csv" = result$cost_summary
  ),
  report
)
message("Published pure-SAEM assessment: ", published)
if (!isTRUE(result$falsification_passed)) {
  message("Final pure-SAEM falsification gate failed.")
  quit(save = "no", status = 4L, runLast = FALSE)
}
quit(save = "no", status = 0L, runLast = FALSE)
