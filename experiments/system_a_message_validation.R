#!/usr/bin/env Rscript

# ODE-free endpoint planning and post-bank validation for the raw System A
# patient-message experiment.  Heavy conditional sampling is performed only by
# the separate Slurm bank arrays.

sab_cli_stop <- function(...) stop(paste0(...), call. = FALSE)

sab_cli_sha256 <- function(path) {
  binary <- Sys.which("sha256sum")
  if (!nzchar(binary)) sab_cli_stop("sha256sum is required.")
  output <- suppressWarnings(system2(
    binary, c("--", shQuote(normalizePath(path, mustWork = TRUE))),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  if ((!is.null(status) && status != 0L) || length(output) != 1L) {
    sab_cli_stop("Could not hash ", path, ".")
  }
  digest <- tolower(strsplit(trimws(output), "[[:space:]]+")[[1L]][[1L]])
  if (!grepl("^[0-9a-f]{64}$", digest)) {
    sab_cli_stop("Invalid SHA-256 output for ", path, ".")
  }
  digest
}

sab_cli_load_anchor_banks <- function(root, patient_ids) {
  paths <- stats::setNames(file.path(
    root, sprintf("patient_%03d", as.integer(patient_ids)),
    "patient_banks.rds"
  ), patient_ids)
  missing <- names(paths)[!file.exists(paths)]
  if (length(missing)) {
    sab_cli_stop("Missing anchor banks for patients: ",
                 paste(missing, collapse = ", "), ".")
  }
  paths <- vapply(paths, normalizePath, character(1L), mustWork = TRUE)
  list(
    paths = paths,
    sha256 = vapply(paths, sab_cli_sha256, character(1L)),
    artifacts = lapply(paths, readRDS)
  )
}

sab_cli_load_candidate_banks <- function(root, endpoints, patient_ids) {
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
      sab_cli_stop("Missing candidate banks for ", endpoint_id,
                   ", patients: ", paste(missing, collapse = ", "), ".")
    }
    endpoint_paths <- vapply(
      endpoint_paths, normalizePath, character(1L), mustWork = TRUE
    )
    paths[[endpoint_id]] <- endpoint_paths
    hashes[[endpoint_id]] <- vapply(
      endpoint_paths, sab_cli_sha256, character(1L)
    )
    artifacts[[endpoint_id]] <- lapply(endpoint_paths, readRDS)
  }
  list(paths = paths, sha256 = hashes, artifacts = artifacts)
}

sab_cli_publish <- function(output_directory, files, report) {
  output_directory <- path.expand(output_directory)
  if (file.exists(output_directory) || dir.exists(output_directory)) {
    sab_cli_stop("Refusing to overwrite output directory ",
                 output_directory, ".")
  }
  parent <- dirname(output_directory)
  if (!dir.exists(parent) && !dir.create(parent, recursive = TRUE) &&
      !dir.exists(parent)) {
    sab_cli_stop("Could not create output parent ", parent, ".")
  }
  stage <- paste0(output_directory, ".partial")
  if (file.exists(stage) || dir.exists(stage) ||
      !dir.create(stage, recursive = FALSE, mode = "0750")) {
    sab_cli_stop("Could not create clean staging directory ", stage, ".")
  }
  for (name in names(files)) {
    value <- files[[name]]
    path <- file.path(stage, name)
    if (grepl("\\.rds$", name)) {
      saveRDS(value, path, version = 3L, compress = "gzip")
    } else if (grepl("\\.csv$", name)) {
      utils::write.csv(value, path, row.names = FALSE, quote = TRUE, na = "NA")
    } else {
      writeLines(as.character(value), path, useBytes = TRUE)
    }
  }
  writeLines(report, file.path(stage, "REPORT.txt"), useBytes = TRUE)
  manifest_paths <- list.files(stage, full.names = TRUE)
  manifest <- data.frame(
    file = basename(manifest_paths),
    sha256 = vapply(manifest_paths, sab_cli_sha256, character(1L)),
    stringsAsFactors = FALSE
  )
  utils::write.csv(manifest, file.path(stage, "manifest.csv"),
                   row.names = FALSE, quote = TRUE)
  writeLines("completed", file.path(stage, "COMPLETED"), useBytes = TRUE)
  if (!file.rename(stage, output_directory)) {
    sab_cli_stop("Atomic publication failed for ", output_directory, ".")
  }
  normalizePath(output_directory, mustWork = TRUE)
}

arguments <- commandArgs(trailingOnly = TRUE)
if (!length(arguments) %in% c(4L, 6L) ||
    !arguments[[1L]] %in% c("plan", "assess")) {
  sab_cli_stop(
    "Usage: system_a_message_validation.R plan <anchor.rds> ",
    "<anchor-bank-root> <output-dir> OR assess <anchor.rds> ",
    "<anchor-bank-root> <endpoint-plan.rds> <candidate-bank-root> ",
    "<output-dir>."
  )
}
mode <- arguments[[1L]]
if ((mode == "plan" && length(arguments) != 4L) ||
    (mode == "assess" && length(arguments) != 6L)) {
  sab_cli_stop("The selected mode has the wrong number of arguments.")
}
workspace_root <- normalizePath(
  Sys.getenv("SAB_WORKSPACE_ROOT", unset = ""), mustWork = TRUE
)
project_root <- normalizePath(
  Sys.getenv("SAB_PROJECT_ROOT", unset = ""), mustWork = TRUE
)
anchor_path <- normalizePath(arguments[[2L]], mustWork = TRUE)
anchor_bank_root <- normalizePath(arguments[[3L]], mustWork = TRUE)

modules <- file.path(project_root, "R", c(
  "system_a_adapter.R", "system_a_saem_anchor.R", "patient_messages.R",
  "system_a_message_validation.R"
))
if (any(!file.exists(modules))) sab_cli_stop("Required R modules are absent.")
invisible(lapply(modules, sys.source, envir = globalenv(), keep.source = FALSE))

adapter <- sab_load_system_a_adapter(workspace_root)
canonical_anchor <- readRDS(anchor_path)
sab_validate_system_a_saem_anchor(canonical_anchor, adapter)
design <- sab_system_a_message_design()
anchor_banks <- sab_cli_load_anchor_banks(
  anchor_bank_root, design$patient_ids
)
anchor_sha256 <- sab_cli_sha256(anchor_path)

if (mode == "plan") {
  output_directory <- arguments[[4L]]
  result <- sab_system_a_plan_message_endpoints(
    adapter = adapter,
    canonical_anchor = canonical_anchor,
    anchor_bank_artifacts = anchor_banks$artifacts,
    anchor_path = anchor_path,
    anchor_sha256 = anchor_sha256,
    anchor_bank_paths = anchor_banks$paths,
    anchor_bank_sha256 = anchor_banks$sha256
  )
  report <- c(
    "System A fixed-psi raw-message endpoint plan",
    "",
    "This ODE-free pilot used only chain_01 and chain_02 from each anchor bank.",
    "chain_03 and chain_04 remain held out for the final comparison.",
    "The two endpoints are the separated SAEM and prior-centred components",
    "of h = 0.5*g_SAEM + 0.5*g_priorcentral. Failed pilot gates cannot later",
    "count as validation.",
    "",
    capture.output(print(result$selection, row.names = FALSE))
  )
  published <- sab_cli_publish(
    output_directory,
    list(
      "endpoint_plan.rds" = result$plan,
      "pilot_diagnostics.csv" = result$pilot_diagnostics,
      "component_identity.csv" = result$component_identity,
      "reference_bank_performance.csv" =
        result$reference_bank_performance,
      "selection.csv" = result$selection,
      "endpoints.csv" = result$plan$endpoints
    ),
    report
  )
  message("Published endpoint plan: ", published)
} else {
  plan_path <- normalizePath(arguments[[4L]], mustWork = TRUE)
  candidate_bank_root <- normalizePath(arguments[[5L]], mustWork = TRUE)
  output_directory <- arguments[[6L]]
  plan <- readRDS(plan_path)
  plan_sha256 <- sab_cli_sha256(plan_path)
  if (!identical(plan$anchor$path, anchor_path) ||
      !identical(plan$anchor$sha256, anchor_sha256) ||
      !identical(plan$anchor_banks$paths, anchor_banks$paths) ||
      !identical(plan$anchor_banks$sha256, anchor_banks$sha256)) {
    sab_cli_stop("Anchor or held-out bank artifacts differ from the pilot plan.")
  }
  candidate_banks <- sab_cli_load_candidate_banks(
    candidate_bank_root, plan$endpoints, design$patient_ids
  )
  result <- sab_system_a_assess_messages(
    adapter = adapter,
    canonical_anchor = canonical_anchor,
    plan = plan,
    plan_sha256 = plan_sha256,
    anchor_bank_artifacts = anchor_banks$artifacts,
    candidate_bank_artifacts = candidate_banks$artifacts
  )
  result$provenance <- list(
    anchor = c(path = anchor_path, sha256 = anchor_sha256),
    endpoint_plan = c(path = plan_path, sha256 = plan_sha256),
    anchor_bank_paths = anchor_banks$paths,
    anchor_bank_sha256 = anchor_banks$sha256,
    candidate_bank_paths = candidate_banks$paths,
    candidate_bank_sha256 = candidate_banks$sha256,
    module_sha256 = stats::setNames(
      vapply(modules, sab_cli_sha256, character(1L)), basename(modules)
    )
  )
  report <- c(
    "System A 12-patient fixed-psi raw-message validation",
    "",
    "Forward raw importance uses only held-out anchor chains 03-04.",
    "Reverse importance uses separately simulated candidate-target banks.",
    "The candidate-augmented bridge estimates reuse the held-out h chains,",
    "so raw-versus-bridge agreement is a diagnostic, not an independent",
    "replicate comparison. All population densities are normalized callback",
    "evaluations with the observed treatment covariate retained.",
    "This panel validates only the selected fixed-psi endpoints, not the full",
    "115-patient posterior or uncertainty in shared dynamic parameters.",
    "",
    paste("Status:", result$status),
    "",
    capture.output(print(result$endpoint_summary, row.names = FALSE))
  )
  published <- sab_cli_publish(
    output_directory,
    list(
      "result.rds" = result,
      "endpoint_summary.csv" = result$endpoint_summary,
      "patient_summary.csv" = result$patient_summary,
      "component_identity.csv" = result$component_identity,
      "bridge_component_identity.csv" = result$bridge_component_identity,
      "mcmc_diagnostics.csv" = result$mcmc_diagnostics,
      "reference_bank_performance.csv" =
        result$reference_bank_performance,
      "directional_diagnostics.csv" = result$directional_diagnostics,
      "bridge_chain_pairs.csv" = result$bridge_chain_pairs
    ),
    report
  )
  message("Published message validation: ", published)
}
