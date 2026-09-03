#!/usr/bin/env Rscript

# Build independent exact-target conditional MCMC banks for one predeclared
# System A patient.  With three command-line arguments the target fixes psi at
# the canonical SAEM anchor and uses the declared defensive population
# reference.  An optional endpoint-plan argument instead targets one of the
# population endpoints selected by the separate ODE-free pilot.  Slurm
# assigns patient/endpoint pairs across array tasks; chains for one pair run
# sequentially inside its one-CPU task.
# This worker does not evaluate messages, choose population endpoints, fit a
# transport map, or use comparator output.  The anchor-bank target uses the
# fixed defensive reference h_i = 0.5 g_SAEM,i + 0.5 g_priorcentral,i and an
# independence proposal from h_i.  Candidate endpoint banks retain the
# diagonal-Gaussian pCN kernel.  These System A samplers place no family
# restriction on population densities evaluated by the generic message code.

sab_worker_stop <- function(...) {
  stop(paste0(...), call. = FALSE)
}

sab_worker_integer <- function(value, label, minimum = 0L,
                               maximum = .Machine$integer.max) {
  parsed <- suppressWarnings(as.numeric(value))
  if (length(parsed) != 1L || !is.finite(parsed) || parsed != floor(parsed) ||
      parsed < minimum || parsed > maximum) {
    sab_worker_stop(label, " must be one integer in [", minimum, ", ",
                    maximum, "].")
  }
  as.integer(parsed)
}

sab_worker_number <- function(value, label, lower, upper,
                              lower_open = FALSE, upper_open = FALSE) {
  parsed <- suppressWarnings(as.numeric(value))
  lower_bad <- if (lower_open) parsed <= lower else parsed < lower
  upper_bad <- if (upper_open) parsed >= upper else parsed > upper
  if (length(parsed) != 1L || !is.finite(parsed) || lower_bad || upper_bad) {
    left <- if (lower_open) "(" else "["
    right <- if (upper_open) ")" else "]"
    sab_worker_stop(label, " must be one finite number in ", left, lower,
                    ", ", upper, right, ".")
  }
  parsed
}

sab_worker_env <- function(name, default = NULL) {
  value <- Sys.getenv(name, unset = "")
  if (!nzchar(value)) {
    if (is.null(default)) sab_worker_stop("Missing environment variable ", name, ".")
    return(default)
  }
  value
}

sab_worker_sha256 <- function(path) {
  binary <- Sys.which("sha256sum")
  if (!nzchar(binary)) sab_worker_stop("sha256sum is required.")
  output <- suppressWarnings(system2(
    binary, c("--", shQuote(normalizePath(path, mustWork = TRUE))),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  if ((!is.null(status) && status != 0L) || length(output) != 1L) {
    sab_worker_stop("Could not calculate SHA-256 for ", path, ".")
  }
  digest <- tolower(strsplit(trimws(output), "[[:space:]]+")[[1L]][[1L]])
  if (!grepl("^[0-9a-f]{64}$", digest)) {
    sab_worker_stop("sha256sum returned an invalid digest for ", path, ".")
  }
  digest
}

sab_worker_validate_production_source <- function(source) {
  required_artifacts <- c(
    "population", "omega", "individuals", "run_config", "input_data",
    "audit"
  )
  files <- if (is.list(source)) source$artifact_files else NULL
  hashes <- if (is.list(source)) source$artifact_sha256 else NULL
  audit <- if (is.list(source)) source$audit else NULL
  run_config <- if (is.list(source)) source$run_config else NULL
  fit <- if (is.list(run_config)) run_config$fit_control else NULL
  row_order <- if (is.list(source)) source$row_order_validation else NULL
  row_error_fields <- c(
    "max_abs_map_csv_error", "max_abs_fixed_csv_error",
    "max_abs_omega_csv_error"
  )
  valid <-
    is.character(source$run_directory) &&
    length(source$run_directory) == 1L &&
    !is.na(source$run_directory) && dir.exists(source$run_directory) &&
    is.character(files) && identical(names(files), required_artifacts) &&
    all(file.exists(files)) &&
    is.character(hashes) && identical(names(hashes), required_artifacts) &&
    all(grepl("^[0-9a-f]{64}$", hashes)) &&
    is.list(audit) && identical(audit$status, "completed") &&
    identical(audit$path, unname(files[["audit"]])) &&
    identical(audit$sha256, unname(hashes[["audit"]])) &&
    identical(audit$n_patients, 115L) &&
    is.list(run_config) &&
    identical(run_config$path, unname(files[["run_config"]])) &&
    identical(run_config$sha256, unname(hashes[["run_config"]])) &&
    is.character(run_config$run_tag) && length(run_config$run_tag) == 1L &&
    !is.na(run_config$run_tag) && nzchar(run_config$run_tag) &&
    is.list(fit) && is.integer(fit$k1) && length(fit$k1) == 1L &&
    is.integer(fit$k2) && length(fit$k2) == 1L &&
    fit$k1 >= 300L && fit$k2 >= 100L &&
    is.character(run_config$hiv_data_path) &&
    length(run_config$hiv_data_path) == 1L &&
    file.exists(run_config$hiv_data_path) &&
    is.character(run_config$hiv_data_sha256) &&
    length(run_config$hiv_data_sha256) == 1L &&
    grepl("^[0-9a-f]{64}$", run_config$hiv_data_sha256) &&
    is.list(row_order) && identical(
      row_order$schema_version,
      "sab_system_a_saem_row_order_validation_v1"
    ) &&
    isTRUE(row_order$validated) &&
    identical(row_order$method, "saemix_internal_index_binding") &&
    is.character(source$input_patient_ids) &&
    length(source$input_patient_ids) == 115L &&
    !anyNA(source$input_patient_ids) &&
    !anyDuplicated(source$input_patient_ids) &&
    identical(row_order$patient_ids, source$input_patient_ids) &&
    is.character(row_order$fit_path) && length(row_order$fit_path) == 1L &&
    !is.na(row_order$fit_path) && file.exists(row_order$fit_path) &&
    identical(
      normalizePath(row_order$fit_path, mustWork = TRUE),
      normalizePath(file.path(
        source$run_directory, paste0(run_config$run_tag, ".fit.rds")
      ), mustWork = TRUE)
    ) &&
    is.character(row_order$fit_sha256) &&
    length(row_order$fit_sha256) == 1L &&
    grepl("^[0-9a-f]{64}$", row_order$fit_sha256) &&
    all(vapply(row_order[row_error_fields], function(value) {
      is.numeric(value) && length(value) == 1L && is.finite(value) &&
        value >= 0
    }, logical(1L)))
  if (!isTRUE(valid)) {
    sab_worker_stop(
      "Banks require a completed production SAEM anchor with K1 >= 300 ",
      "and K2 >= 100; smoke or malformed run provenance is forbidden."
    )
  }
  observed <- vapply(files, sab_worker_sha256, character(1L))
  if (!identical(observed, hashes) ||
      !identical(
        sab_worker_sha256(run_config$hiv_data_path),
        run_config$hiv_data_sha256
      ) ||
      !identical(
        sab_worker_sha256(row_order$fit_path), row_order$fit_sha256
      )) {
    sab_worker_stop("SAEM source artifacts or input data changed after conversion.")
  }
  list(
    audit_status = audit$status,
    k1 = fit$k1,
    k2 = fit$k2,
    paths = c(
      files, hiv_source_data = run_config$hiv_data_path,
      saem_fit = row_order$fit_path
    ),
    sha256 = c(
      hashes, hiv_source_data = run_config$hiv_data_sha256,
      saem_fit = row_order$fit_sha256
    )
  )
}

sab_worker_atomic_save_rds <- function(object, path) {
  if (file.exists(path)) sab_worker_stop("Refusing to overwrite ", path, ".")
  temporary <- paste0(path, ".partial")
  if (file.exists(temporary)) {
    sab_worker_stop("Refusing to overwrite stale staging file ", temporary, ".")
  }
  saveRDS(object, temporary, version = 3L, compress = "gzip")
  if (!file.rename(temporary, path)) {
    sab_worker_stop("Atomic RDS publication failed for ", path, ".")
  }
  invisible(path)
}

sab_worker_atomic_write_csv <- function(value, path) {
  if (file.exists(path)) sab_worker_stop("Refusing to overwrite ", path, ".")
  temporary <- paste0(path, ".partial")
  if (file.exists(temporary)) {
    sab_worker_stop("Refusing to overwrite stale staging file ", temporary, ".")
  }
  utils::write.csv(value, temporary, row.names = FALSE, quote = TRUE,
                   na = "NA")
  if (!file.rename(temporary, path)) {
    sab_worker_stop("Atomic CSV publication failed for ", path, ".")
  }
  invisible(path)
}

sab_worker_validate_anchor <- function(anchor, adapter, patient_id) {
  required <- c(
    "schema_version", "target_fingerprint", "source", "coordinate_names",
    "patient_ids", "eta", "psi", "patient_states", "validation"
  )
  if (!inherits(anchor, "sab_system_a_saem_anchor") ||
      !is.list(anchor) || length(setdiff(required, names(anchor))) ||
      !identical(anchor$schema_version, "sab_system_a_saem_anchor_v1")) {
    sab_worker_stop("The anchor is not a canonical System A SAEM anchor.")
  }
  if (!identical(anchor$target_fingerprint, adapter$target_fingerprint)) {
    sab_worker_stop("Anchor and exact target fingerprints differ.")
  }
  expected_coordinates <- adapter$coordinate_names[c(
    "local", "population", "global"
  )]
  if (!is.list(anchor$coordinate_names) ||
      !identical(anchor$coordinate_names, expected_coordinates)) {
    sab_worker_stop("Anchor coordinate names/order differ from the exact target.")
  }
  if (!is.character(anchor$patient_ids) || length(anchor$patient_ids) != 115L ||
      anyNA(anchor$patient_ids) || anyDuplicated(anchor$patient_ids) ||
      !patient_id %in% anchor$patient_ids) {
    sab_worker_stop("Anchor patient IDs are malformed or incomplete.")
  }
  if (!is.numeric(anchor$eta) ||
      !identical(names(anchor$eta), expected_coordinates$population) ||
      any(!is.finite(anchor$eta)) || !isTRUE(adapter$eta_in_domain(anchor$eta))) {
    sab_worker_stop("Anchor eta is not finite, canonical, and in target support.")
  }
  if (!is.numeric(anchor$psi) ||
      !identical(names(anchor$psi), expected_coordinates$global) ||
      any(!is.finite(anchor$psi)) || !isTRUE(adapter$psi_in_domain(anchor$psi))) {
    sab_worker_stop("Anchor psi is not finite, canonical, and in target support.")
  }
  states <- anchor$patient_states
  if (!is.matrix(states) || !is.numeric(states) ||
      !identical(dim(states), c(115L, length(expected_coordinates$local))) ||
      !identical(rownames(states), anchor$patient_ids) ||
      !identical(colnames(states), expected_coordinates$local) ||
      any(!is.finite(states))) {
    sab_worker_stop("Anchor patient states are not a finite canonical 115 x 8 matrix.")
  }
  population_validation <- if (is.list(anchor$validation)) {
    anchor$validation$population_log_density
  } else {
    NULL
  }
  if (!is.numeric(population_validation) ||
      !identical(names(population_validation), anchor$patient_ids) ||
      any(!is.finite(population_validation))) {
    sab_worker_stop("Anchor population-density validation is malformed.")
  }
  exact_validation <- if (is.list(anchor$validation)) {
    anchor$validation$exact_target
  } else {
    NULL
  }
  exact_summary <- if (is.list(exact_validation)) {
    exact_validation$summary
  } else {
    NULL
  }
  exact_ledger <- if (is.list(exact_validation)) {
    exact_validation$ledger
  } else {
    NULL
  }
  ledger_columns <- c(
    "patient_id", "population_log_density", "prediction_calls",
    "ode_integrations", "prediction_ok", "prediction_reason",
    "callback_error", "likelihood_calls", "loglik", "finite_loglik",
    "ordinary_target_rejection", "passed"
  )
  summary_columns <- c(
    "patients", "prediction_calls", "ode_integrations",
    "unknown_ode_counts", "prediction_failures", "likelihood_calls",
    "nonfinite_loglik", "ordinary_target_rejections", "callback_errors",
    "completed_without_callback_error", "passed"
  )
  allowed_reasons <- c(
    "ok", "nonfinite_natural_parameter", "invalid_equilibrium",
    "ode_failure", "invalid_ode_output", "invalid_ode_times"
  )
  exact_structure_ok <-
    is.list(anchor$source) && is.list(anchor$validation) &&
    inherits(exact_validation, "sab_system_a_saem_exact_validation") &&
    identical(
      exact_validation$schema_version,
      "sab_system_a_saem_exact_validation_v1"
    ) &&
    identical(exact_validation$target_fingerprint,
              anchor$target_fingerprint) &&
    is.data.frame(exact_ledger) &&
    identical(names(exact_ledger), ledger_columns) &&
    nrow(exact_ledger) == 115L &&
    identical(as.character(exact_ledger$patient_id), anchor$patient_ids) &&
    is.data.frame(exact_summary) &&
    identical(names(exact_summary), summary_columns) &&
    nrow(exact_summary) == 1L &&
    is.logical(exact_validation$completed_without_callback_error) &&
    length(exact_validation$completed_without_callback_error) == 1L &&
    !is.na(exact_validation$completed_without_callback_error) &&
    is.logical(exact_validation$passed) &&
    length(exact_validation$passed) == 1L &&
    !is.na(exact_validation$passed)
  if (!isTRUE(exact_structure_ok)) {
    sab_worker_stop(
      "Anchor provenance or all-patient exact-target validation is malformed."
    )
  }
  ledger_ok <-
    is.numeric(exact_ledger$population_log_density) &&
    all(is.finite(exact_ledger$population_log_density)) &&
    all(exact_ledger$prediction_calls == 1L) &&
    is.numeric(exact_ledger$ode_integrations) &&
    all(is.finite(exact_ledger$ode_integrations)) &&
    all(exact_ledger$ode_integrations >= 0) &&
    is.logical(exact_ledger$prediction_ok) &&
    !anyNA(exact_ledger$prediction_ok) &&
    is.character(exact_ledger$prediction_reason) &&
    !anyNA(exact_ledger$prediction_reason) &&
    all(exact_ledger$prediction_reason %in% allowed_reasons) &&
    is.logical(exact_ledger$callback_error) &&
    !anyNA(exact_ledger$callback_error) &&
    !any(exact_ledger$callback_error) &&
    all(exact_ledger$likelihood_calls == 1L) &&
    is.numeric(exact_ledger$loglik) &&
    !anyNA(exact_ledger$loglik) &&
    is.logical(exact_ledger$finite_loglik) &&
    !anyNA(exact_ledger$finite_loglik) &&
    is.logical(exact_ledger$ordinary_target_rejection) &&
    !anyNA(exact_ledger$ordinary_target_rejection) &&
    is.logical(exact_ledger$passed) && !anyNA(exact_ledger$passed) &&
    identical(exact_ledger$prediction_ok,
              exact_ledger$prediction_reason == "ok") &&
    identical(exact_ledger$finite_loglik, is.finite(exact_ledger$loglik)) &&
    identical(
      exact_ledger$passed,
      exact_ledger$prediction_ok & exact_ledger$finite_loglik
    ) &&
    identical(exact_ledger$ordinary_target_rejection,
              !exact_ledger$passed) &&
    all(exact_ledger$loglik[!exact_ledger$passed] == -Inf)
  expected_passed <- all(exact_ledger$passed)
  summary_ok <-
    exact_summary$patients[[1L]] == 115L &&
    exact_summary$prediction_calls[[1L]] == 115L &&
    exact_summary$ode_integrations[[1L]] ==
      sum(exact_ledger$ode_integrations) &&
    exact_summary$unknown_ode_counts[[1L]] == 0L &&
    exact_summary$prediction_failures[[1L]] ==
      sum(!exact_ledger$prediction_ok) &&
    exact_summary$likelihood_calls[[1L]] == 115L &&
    exact_summary$nonfinite_loglik[[1L]] ==
      sum(!exact_ledger$finite_loglik) &&
    exact_summary$ordinary_target_rejections[[1L]] ==
      sum(exact_ledger$ordinary_target_rejection) &&
    exact_summary$callback_errors[[1L]] == 0L &&
    isTRUE(exact_summary$completed_without_callback_error[[1L]]) &&
    isTRUE(exact_validation$completed_without_callback_error) &&
    identical(exact_summary$passed[[1L]], expected_passed) &&
    identical(exact_validation$passed, expected_passed)
  if (!isTRUE(ledger_ok) || !isTRUE(summary_ok)) {
    sab_worker_stop(
      "Anchor exact-target validation contains callback/contract errors."
    )
  }
  recalculated_population_density <- adapter$log_population_density(
    patient_id, states[patient_id, ], anchor$eta
  )
  if (!is.numeric(recalculated_population_density) ||
      length(recalculated_population_density) != 1L ||
      !is.finite(recalculated_population_density) ||
      !isTRUE(all.equal(
        as.numeric(recalculated_population_density),
        as.numeric(population_validation[[patient_id]]),
        tolerance = 1e-12
      )) ||
      !isTRUE(all.equal(
        as.numeric(recalculated_population_density),
        exact_ledger$population_log_density[
          exact_ledger$patient_id == patient_id
        ],
        tolerance = 1e-12
      ))) {
    sab_worker_stop("Selected anchor population-density validation is stale.")
  }
  production_source <- sab_worker_validate_production_source(anchor$source)
  invisible(production_source)
}

sab_worker_prior_directions <- function(local_names) {
  # Rows of a fixed order-eight Hadamard design.  Dense directions disperse
  # starts across the local prior without a data-dependent search.
  if (length(local_names) != 8L) {
    sab_worker_stop("The predeclared System A start design requires 8 local coordinates.")
  }
  directions <- rbind(
    c(1, 1, 1, 1, 1, 1, 1, 1),
    c(1, -1, 1, -1, 1, -1, 1, -1),
    c(1, 1, -1, -1, 1, 1, -1, -1),
    c(1, -1, -1, 1, 1, -1, -1, 1),
    c(1, 1, 1, 1, -1, -1, -1, -1),
    c(1, -1, 1, -1, -1, 1, -1, 1),
    c(1, 1, -1, -1, -1, -1, 1, 1),
    c(1, -1, -1, 1, -1, 1, 1, -1)
  )
  colnames(directions) <- local_names
  directions
}

sab_worker_candidates <- function(saem_x, population_mean, population_sd,
                                  offset_sd, chain_index) {
  local_names <- names(population_mean)
  directions <- sab_worker_prior_directions(local_names)
  offsets <- do.call(rbind, lapply(seq_len(nrow(directions)), function(index) {
    rbind(
      population_mean + offset_sd * population_sd * directions[index, ],
      population_mean - offset_sd * population_sd * directions[index, ]
    )
  }))
  rownames(offsets) <- unlist(lapply(seq_len(nrow(directions)), function(index) {
    c(sprintf("prior_offset_%02d_plus", index),
      sprintf("prior_offset_%02d_minus", index))
  }))
  candidates <- rbind(saem_x = saem_x, population_mean = population_mean,
                      offsets)
  desired <- if (chain_index == 1L) {
    1L
  } else if (chain_index == 2L) {
    2L
  } else {
    2L + ((chain_index - 3L) %% nrow(offsets)) + 1L
  }
  candidates[c(desired, setdiff(seq_len(nrow(candidates)), desired)), ,
             drop = FALSE]
}

sab_worker_defensive_candidates <- function(
    saem_x, saem_mean, saem_sd, prior_mean, prior_sd,
    offset_sd, chain_index) {
  local_names <- names(saem_mean)
  directions <- sab_worker_prior_directions(local_names)
  direction_index <- (chain_index - 1L) %% nrow(directions) + 1L
  direction <- directions[direction_index, ]
  preferred <- switch(
    as.character((chain_index - 1L) %% 4L + 1L),
    "1" = saem_x,
    "2" = prior_mean,
    "3" = saem_mean + offset_sd * saem_sd * direction,
    "4" = prior_mean + offset_sd * prior_sd * direction
  )
  names(preferred) <- local_names
  candidates <- rbind(
    requested = preferred,
    saem_state = saem_x,
    saem_component_mean = saem_mean,
    priorcentral_component_mean = prior_mean,
    saem_opposite_offset = saem_mean - offset_sd * saem_sd * direction,
    priorcentral_opposite_offset =
      prior_mean - offset_sd * prior_sd * direction
  )
  candidates[!duplicated(as.data.frame(candidates)), , drop = FALSE]
}

arguments <- commandArgs(trailingOnly = TRUE)
if (!length(arguments) %in% c(3L, 4L)) {
  sab_worker_stop(
    "Usage: system_a_anchor_bank_patient.R <anchor.rds> <output-root> ",
    "<zero-based-array-index> [endpoint-plan.rds]."
  )
}
command_line <- commandArgs(trailingOnly = FALSE)
file_argument <- grep("^--file=", command_line, value = TRUE)
if (length(file_argument) != 1L) {
  sab_worker_stop("Could not identify this worker's source path.")
}
worker_script_path <- normalizePath(
  sub("^--file=", "", file_argument), mustWork = TRUE
)

workspace_root <- normalizePath(
  sab_worker_env("SAB_WORKSPACE_ROOT"), mustWork = TRUE
)
project_root <- normalizePath(
  sab_worker_env("SAB_PROJECT_ROOT"), mustWork = TRUE
)
anchor_path <- normalizePath(arguments[[1L]], mustWork = TRUE)
output_root <- normalizePath(arguments[[2L]], mustWork = TRUE)
endpoint_plan_path <- if (length(arguments) == 4L) {
  normalizePath(arguments[[4L]], mustWork = TRUE)
} else {
  NULL
}
endpoint_plan <- if (is.null(endpoint_plan_path)) NULL else
  readRDS(endpoint_plan_path)
candidate_mode <- !is.null(endpoint_plan)
if (candidate_mode) {
  endpoint_shape_ok <-
    is.list(endpoint_plan) &&
    identical(
      endpoint_plan$schema_version,
      "sab_system_a_message_endpoint_plan_v1"
    ) &&
    is.data.frame(endpoint_plan$endpoints) &&
    nrow(endpoint_plan$endpoints) == 2L &&
    identical(
      names(endpoint_plan$endpoints),
      c(
        "endpoint_id", "stage", "kind", "axis", "parameter", "sign", "step",
        "pilot_gate_passed", "diagnostic_fallback"
      )
    ) &&
    is.matrix(endpoint_plan$endpoint_eta) &&
    is.numeric(endpoint_plan$endpoint_eta) &&
    nrow(endpoint_plan$endpoint_eta) == 2L &&
    identical(
      rownames(endpoint_plan$endpoint_eta),
      endpoint_plan$endpoints$endpoint_id
    )
  if (!isTRUE(endpoint_shape_ok)) {
    sab_worker_stop("The endpoint plan has an invalid schema or shape.")
  }
}
task_maximum <- if (candidate_mode) {
  12L * nrow(endpoint_plan$endpoints) - 1L
} else {
  11L
}
task_index <- sab_worker_integer(
  arguments[[3L]], "array task index", 0L, task_maximum
)
slurm_task <- Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "")
if (nzchar(slurm_task) &&
    sab_worker_integer(
      slurm_task, "SLURM_ARRAY_TASK_ID", 0L, task_maximum
    ) !=
      task_index) {
  sab_worker_stop("Command-line and Slurm array task indices differ.")
}

audited_patient_ids <- c(
  "3", "6", "16", "20", "55", "71", "74", "88", "105", "111",
  "117", "122"
)
patient_index <- if (candidate_mode) task_index %% 12L else task_index
endpoint_index <- if (candidate_mode) task_index %/% 12L + 1L else NA_integer_
patient_id <- audited_patient_ids[[patient_index + 1L]]

configuration <- list(
  chains = sab_worker_integer(
    sab_worker_env("SAB_BANK_CHAINS", "4"), "SAB_BANK_CHAINS", 2L, 64L
  ),
  warmup = sab_worker_integer(
    sab_worker_env("SAB_BANK_WARMUP", "1000"), "SAB_BANK_WARMUP", 0L
  ),
  draws = sab_worker_integer(
    sab_worker_env("SAB_BANK_DRAWS", "1000"), "SAB_BANK_DRAWS", 1L
  ),
  thin = sab_worker_integer(
    sab_worker_env("SAB_BANK_THIN", "1"), "SAB_BANK_THIN", 1L
  ),
  initial_beta = sab_worker_number(
    sab_worker_env("SAB_BANK_INITIAL_BETA", "0.2"),
    "SAB_BANK_INITIAL_BETA", 0.005, 0.95
  ),
  target_acceptance = sab_worker_number(
    sab_worker_env("SAB_BANK_TARGET_ACCEPTANCE", "0.3"),
    "SAB_BANK_TARGET_ACCEPTANCE", 0, 1, TRUE, TRUE
  ),
  adaptation_block = sab_worker_integer(
    sab_worker_env("SAB_BANK_ADAPTATION_BLOCK", "50"),
    "SAB_BANK_ADAPTATION_BLOCK", 1L
  ),
  base_seed = sab_worker_integer(
    sab_worker_env("SAB_BANK_BASE_SEED", "710003"),
    "SAB_BANK_BASE_SEED", 0L, .Machine$integer.max - 1000000L
  ),
  start_offset_sd = sab_worker_number(
    sab_worker_env("SAB_BANK_START_OFFSET_SD", "1"),
    "SAB_BANK_START_OFFSET_SD", 0.1, 3
  )
)

adapter_path <- file.path(project_root, "R", "system_a_adapter.R")
bank_path <- file.path(project_root, "R", "system_a_patient_bank.R")
if (!file.exists(adapter_path) || !file.exists(bank_path)) {
  sab_worker_stop("Required System A modules are absent.")
}
input_sha256 <- c(
  anchor = sab_worker_sha256(anchor_path),
  "R/system_a_adapter.R" = sab_worker_sha256(adapter_path),
  "R/system_a_patient_bank.R" = sab_worker_sha256(bank_path),
  "experiments/system_a_anchor_bank_patient.R" =
    sab_worker_sha256(worker_script_path)
)
if (candidate_mode) {
  input_sha256 <- c(
    input_sha256,
    endpoint_plan = sab_worker_sha256(endpoint_plan_path)
  )
}
sys.source(adapter_path, envir = globalenv(), keep.source = FALSE)
sys.source(bank_path, envir = globalenv(), keep.source = FALSE)

adapter <- sab_load_system_a_adapter(
  workspace_root = workspace_root, patient_ids = patient_id
)
if (!identical(adapter$numerical_target, "sealed_deSolve_VODE_BDF")) {
  sab_worker_stop("Patient banks require the sealed exact VODE-BDF target.")
}
anchor <- readRDS(anchor_path)
production_source <- sab_worker_validate_anchor(anchor, adapter, patient_id)
if (candidate_mode) {
  plan_ok <-
    identical(endpoint_plan$target_fingerprint, adapter$target_fingerprint) &&
    identical(endpoint_plan$anchor$path, anchor_path) &&
    identical(endpoint_plan$anchor$sha256,
              unname(input_sha256[["anchor"]])) &&
    identical(endpoint_plan$patient_ids, audited_patient_ids) &&
    identical(colnames(endpoint_plan$endpoint_eta),
              adapter$coordinate_names$population) &&
    identical(endpoint_plan$psi, anchor$psi) &&
    all(is.finite(endpoint_plan$endpoint_eta))
  if (!isTRUE(plan_ok)) {
    sab_worker_stop(
      "Endpoint plan, canonical anchor, patients, or target coordinates differ."
    )
  }
}
source_hash_names <- paste0("saem_source/", names(production_source$sha256))
input_sha256 <- c(
  input_sha256,
  stats::setNames(production_source$sha256, source_hash_names)
)
saem_start_validation <- anchor$validation$exact_target$ledger[
  anchor$validation$exact_target$ledger$patient_id == patient_id, ,
  drop = FALSE
]
if (nrow(saem_start_validation) != 1L) {
  sab_worker_stop("Exact validation does not identify one SAEM patient start.")
}
saem_start_rejection_reason <- if (
  saem_start_validation$ordinary_target_rejection
) {
  if (identical(saem_start_validation$prediction_reason, "ok")) {
    "nonfinite_loglik"
  } else {
    saem_start_validation$prediction_reason
  }
} else {
  NA_character_
}

eta <- if (candidate_mode) {
  value <- endpoint_plan$endpoint_eta[endpoint_index, ]
  names(value) <- colnames(endpoint_plan$endpoint_eta)
  value
} else {
  anchor$eta
}
psi <- anchor$psi
if (!isTRUE(adapter$eta_in_domain(eta))) {
  sab_worker_stop("The selected conditional-bank endpoint is outside support.")
}
endpoint <- if (candidate_mode) {
  endpoint_plan$endpoints[endpoint_index, , drop = FALSE]
} else {
  data.frame(
    endpoint_id = "defensive_reference", stage = "reference",
    kind = "equal_defensive_mixture", axis = "reference",
    parameter = NA_character_, sign = 0L, step = 0,
    pilot_gate_passed = TRUE,
    diagnostic_fallback = FALSE, stringsAsFactors = FALSE
  )
}
local_names <- adapter$coordinate_names$local
saem_x <- anchor$patient_states[patient_id, ]
names(saem_x) <- local_names
population_mean <- adapter$population_mean(patient_id, eta)
if (!is.numeric(population_mean) ||
    !identical(names(population_mean), local_names) ||
    any(!is.finite(population_mean))) {
  sab_worker_stop("Exact adapter returned a malformed patient population mean.")
}
scale_names <- grep(
  "^log_omega_", adapter$coordinate_names$population, value = TRUE
)
if (length(scale_names) != length(local_names)) {
  sab_worker_stop("System A population scales do not match local coordinates.")
}
population_sd <- exp(eta[scale_names])
names(population_sd) <- local_names
if (any(!is.finite(population_sd)) || any(population_sd <= 0)) {
  sab_worker_stop("Anchor patient population scales are invalid.")
}
priorcentral_eta <- adapter$prior_reference$eta
if (!is.numeric(priorcentral_eta) ||
    !identical(names(priorcentral_eta), adapter$coordinate_names$population) ||
    any(!is.finite(priorcentral_eta)) ||
    !isTRUE(adapter$eta_in_domain(priorcentral_eta))) {
  sab_worker_stop("The pinned prior-centred population reference is malformed.")
}
priorcentral_mean <- adapter$population_mean(patient_id, priorcentral_eta)
priorcentral_sd <- exp(priorcentral_eta[scale_names])
names(priorcentral_sd) <- local_names
if (!is.numeric(priorcentral_mean) ||
    !identical(names(priorcentral_mean), local_names) ||
    any(!is.finite(priorcentral_mean)) ||
    any(!is.finite(priorcentral_sd)) || any(priorcentral_sd <= 0)) {
  sab_worker_stop("The prior-centred patient population component is invalid.")
}

seed_namespace <- if (candidate_mode) 500000L else 0L
chain_seeds <- configuration$base_seed + seed_namespace + 10000L * task_index +
  seq_len(configuration$chains)
if (any(chain_seeds > .Machine$integer.max)) {
  sab_worker_stop("The deterministic chain seeds exceed the R integer range.")
}
RNGkind("L'Ecuyer-CMRG")
chains <- setNames(vector("list", configuration$chains),
                   sprintf("chain_%02d", seq_len(configuration$chains)))
start_design <- setNames(vector("list", configuration$chains), names(chains))

for (chain_index in seq_len(configuration$chains)) {
  candidates <- if (candidate_mode) {
    sab_worker_candidates(
      saem_x = saem_x,
      population_mean = population_mean,
      population_sd = population_sd,
      offset_sd = configuration$start_offset_sd,
      chain_index = chain_index
    )
  } else {
    sab_worker_defensive_candidates(
      saem_x = saem_x,
      saem_mean = population_mean,
      saem_sd = population_sd,
      prior_mean = priorcentral_mean,
      prior_sd = priorcentral_sd,
      offset_sd = configuration$start_offset_sd,
      chain_index = chain_index
    )
  }
  seed <- as.integer(chain_seeds[[chain_index]])
  message("Patient ", patient_id, ", chain ", chain_index, "/",
          configuration$chains, ", seed ", seed, ".")
  chains[[chain_index]] <- sab_build_system_a_patient_bank(
    adapter = adapter,
    patient_id = patient_id,
    eta = eta,
    psi = psi,
    initial_candidates = candidates,
    warmup = configuration$warmup,
    draws = configuration$draws,
    thin = configuration$thin,
    initial_beta = configuration$initial_beta,
    target_acceptance = configuration$target_acceptance,
    adaptation_block = configuration$adaptation_block,
    seed = seed,
    proposal_mode = if (candidate_mode) {
      "pcn"
    } else {
      "defensive_independence"
    },
    defensive_eta = if (candidate_mode) NULL else priorcentral_eta
  )
  start_design[[chain_index]] <- list(
    requested = rownames(candidates)[[1L]],
    candidate_order = rownames(candidates),
    accepted_initial_candidate = rownames(candidates)[[
      chains[[chain_index]]$initial_candidate_index
    ]]
  )
}

summary_rows <- do.call(rbind, lapply(seq_along(chains), function(index) {
  bank <- chains[[index]]
  total <- bank$ledger[bank$ledger$phase == "total", , drop = FALSE]
  components <- bank$proposal$component_ledger
  component_value <- function(component, field) {
    if (is.null(components)) return(NA_real_)
    row <- components[
      components$phase == "sampling" & components$component == component,
      , drop = FALSE
    ]
    proposals <- row$proposals[[1L]]
    if (field == "acceptance") {
      if (proposals == 0L) return(NA_real_)
      return(row$accepted[[1L]] / proposals)
    }
    if (field == "proposal_fraction") {
      total_proposals <- sum(
        components$proposals[components$phase == "sampling"]
      )
      if (total_proposals == 0L) return(NA_real_)
      return(proposals / total_proposals)
    }
    if (field == "prediction_failure_rate") {
      if (proposals == 0L) return(NA_real_)
      return(row$prediction_failures[[1L]] / proposals)
    }
    if (field == "invalid_rate") {
      if (proposals == 0L) return(NA_real_)
      return(row$nonfinite_loglik[[1L]] / proposals)
    }
    as.numeric(row[[field]][[1L]])
  }
  data.frame(
    patient_id = patient_id,
    chain = names(chains)[[index]],
    seed = chain_seeds[[index]],
    requested_start = start_design[[index]]$requested,
    accepted_initial_candidate =
      start_design[[index]]$accepted_initial_candidate,
    saem_start_exact_passed = saem_start_validation$passed,
    saem_start_rejection_reason = saem_start_rejection_reason,
    warmup_acceptance = unname(bank$acceptance[["warmup"]]),
    sampling_acceptance = unname(bank$acceptance[["sampling"]]),
    final_beta = bank$beta$final,
    min_ess_x = min(bank$ess$x),
    min_ess_x_squared = min(bank$ess$x_squared),
    min_ess_x_per_ode = if (total$ode_integrations[[1L]] > 0) {
      min(bank$ess$x) / total$ode_integrations[[1L]]
    } else {
      NA_real_
    },
    saem_component_proposals = component_value("saem", "proposals"),
    saem_component_proposal_fraction =
      component_value("saem", "proposal_fraction"),
    saem_component_acceptance = component_value("saem", "acceptance"),
    saem_component_prediction_failure_rate =
      component_value("saem", "prediction_failure_rate"),
    saem_component_invalid_rate = component_value("saem", "invalid_rate"),
    saem_component_prediction_failures =
      component_value("saem", "prediction_failures"),
    saem_component_nonfinite_loglik =
      component_value("saem", "nonfinite_loglik"),
    priorcentral_component_proposals =
      component_value("priorcentral", "proposals"),
    priorcentral_component_proposal_fraction =
      component_value("priorcentral", "proposal_fraction"),
    priorcentral_component_acceptance =
      component_value("priorcentral", "acceptance"),
    priorcentral_component_prediction_failure_rate =
      component_value("priorcentral", "prediction_failure_rate"),
    priorcentral_component_invalid_rate =
      component_value("priorcentral", "invalid_rate"),
    priorcentral_component_prediction_failures =
      component_value("priorcentral", "prediction_failures"),
    priorcentral_component_nonfinite_loglik =
      component_value("priorcentral", "nonfinite_loglik"),
    exact_prediction_calls = total$prediction_calls,
    exact_ode_integrations = total$ode_integrations,
    population_density_rejections = total$population_density_rejections,
    prediction_failures = total$prediction_failures,
    nonfinite_loglik = total$nonfinite_loglik,
    stringsAsFactors = FALSE
  )
}))
rownames(summary_rows) <- NULL

artifact <- structure(list(
  schema_version = if (candidate_mode) {
    "sab_system_a_endpoint_patient_banks_v1"
  } else {
    "sab_system_a_reference_patient_banks_v1"
  },
  bank_role = if (candidate_mode) {
    "candidate_endpoint"
  } else {
    "defensive_reference"
  },
  patient_id = patient_id,
  patient_context = list(
    population_covariates = c(
      treat_nelf = unname(adapter$treatment[[patient_id]])
    ),
    manifest = adapter$patient_manifest[1L, , drop = FALSE]
  ),
  audited_patient_ids = audited_patient_ids,
  array_task_index = task_index,
  conditional_target = list(
    endpoint = endpoint,
    eta = eta,
    psi = psi,
    reference = if (candidate_mode) NULL else list(
      kind = "equal_defensive_mixture",
      weights = c(saem = 0.5, priorcentral = 0.5),
      eta_components = list(
        saem = anchor$eta,
        priorcentral = priorcentral_eta
      )
    ),
    endpoint_plan = if (candidate_mode) {
      list(path = endpoint_plan_path,
           sha256 = unname(input_sha256[["endpoint_plan"]]))
    } else {
      NULL
    }
  ),
  target = list(
    numerical_target = adapter$numerical_target,
    target_fingerprint = adapter$target_fingerprint,
    likelihood_signature = adapter$likelihood_signature,
    oracle_audit = adapter$oracle_audit
  ),
  anchor = list(
    path = anchor_path,
    sha256 = unname(input_sha256[["anchor"]]),
    schema_version = anchor$schema_version,
    source = anchor$source,
    production_gate = production_source[c("audit_status", "k1", "k2")],
    validation = anchor$validation,
    eta = anchor$eta,
    psi = psi,
    saem_patient_state = saem_x,
    saem_patient_state_validation = saem_start_validation
  ),
  configuration = configuration,
  chain_seeds = setNames(as.integer(chain_seeds), names(chains)),
  start_design = start_design,
  chains = chains,
  summary = summary_rows,
  provenance = list(
    source_sha256 = input_sha256[names(input_sha256) != "anchor"],
    slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = NA_character_),
    slurm_array_job_id = Sys.getenv(
      "SLURM_ARRAY_JOB_ID", unset = NA_character_
    ),
    slurm_array_task_id = Sys.getenv(
      "SLURM_ARRAY_TASK_ID", unset = NA_character_
    ),
    hostname = Sys.info()[["nodename"]],
    completed_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
    r_version = R.version.string,
    session_info = utils::capture.output(utils::sessionInfo())
  )
), class = c("sab_system_a_patient_banks_artifact", "list"))

final_input_sha256 <- c(
  anchor = sab_worker_sha256(anchor_path),
  "R/system_a_adapter.R" = sab_worker_sha256(adapter_path),
  "R/system_a_patient_bank.R" = sab_worker_sha256(bank_path),
  "experiments/system_a_anchor_bank_patient.R" =
    sab_worker_sha256(worker_script_path),
  stats::setNames(
    vapply(production_source$paths, sab_worker_sha256, character(1L)),
    source_hash_names
  )
)
if (candidate_mode) {
  final_input_sha256 <- append(
    final_input_sha256,
    c(endpoint_plan = sab_worker_sha256(endpoint_plan_path)),
    after = 4L
  )
}
if (!identical(final_input_sha256, input_sha256)) {
  sab_worker_stop("Anchor or source files changed while the bank task ran.")
}

patient_label <- sprintf("patient_%03d", as.integer(patient_id))
relative_directory <- if (candidate_mode) {
  file.path(endpoint$endpoint_id[[1L]], patient_label)
} else {
  patient_label
}
final_directory <- file.path(output_root, relative_directory)
if (file.exists(final_directory) || dir.exists(final_directory)) {
  sab_worker_stop("Refusing to overwrite completed output ", final_directory, ".")
}
stage_directory <- file.path(
  output_root,
  paste0(".", endpoint$endpoint_id[[1L]], "_", patient_label, ".job_",
         Sys.getenv("SLURM_JOB_ID", unset = "unknown"), ".tmp")
)
if (!dir.create(stage_directory, recursive = FALSE, mode = "0750")) {
  sab_worker_stop("Could not create unique staging directory ",
                  stage_directory, ".")
}

rds_path <- file.path(stage_directory, "patient_banks.rds")
summary_path <- file.path(stage_directory, "summary.csv")
manifest_path <- file.path(stage_directory, "manifest.csv")
sab_worker_atomic_save_rds(artifact, rds_path)
sab_worker_atomic_write_csv(summary_rows, summary_path)

manifest <- data.frame(
  field = c(
    "schema_version", "bank_role", "endpoint_id", "axis", "parameter",
    "stage", "kind", "sign", "step", "pilot_gate_passed",
    "diagnostic_fallback",
    "patient_id", "array_task_index", "target_fingerprint",
    "numerical_target", "anchor_path", "anchor_sha256", "chains", "warmup",
    "draws_per_chain", "thin", "base_seed", "saem_k1", "saem_k2",
    "saem_audit_status", "saem_start_exact_passed",
    "saem_start_rejection_reason", "rds_sha256", "summary_sha256"
  ),
  value = c(
    artifact$schema_version, artifact$bank_role, endpoint$endpoint_id,
    endpoint$axis, endpoint$parameter, endpoint$stage, endpoint$kind,
    endpoint$sign, endpoint$step,
    endpoint$pilot_gate_passed, endpoint$diagnostic_fallback,
    patient_id, task_index,
    adapter$target_fingerprint, adapter$numerical_target, anchor_path,
    artifact$anchor$sha256, configuration$chains, configuration$warmup,
    configuration$draws, configuration$thin, configuration$base_seed,
    production_source$k1, production_source$k2,
    production_source$audit_status,
    saem_start_validation$passed,
    saem_start_rejection_reason,
    sab_worker_sha256(rds_path), sab_worker_sha256(summary_path)
  ),
  stringsAsFactors = FALSE
)
sab_worker_atomic_write_csv(manifest, manifest_path)

final_parent <- dirname(final_directory)
if (!dir.exists(final_parent) &&
    !dir.create(final_parent, recursive = TRUE, mode = "0750") &&
    !dir.exists(final_parent)) {
  sab_worker_stop("Could not create endpoint output directory ",
                  final_parent, ".")
}
if (file.exists(final_directory) || dir.exists(final_directory) ||
    !file.rename(stage_directory, final_directory)) {
  sab_worker_stop("Atomic patient-directory publication failed for ",
                  final_directory, ".")
}
message("Published exact-target fixed-parameter MCMC banks: ", final_directory)
