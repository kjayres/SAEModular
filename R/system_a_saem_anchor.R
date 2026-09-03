# Strict conversion of HIV SAEM CSV artifacts into a System A anchor.
#
# This module contains no model implementation and performs no fitting.  It
# translates one completed SAEM export into the exact coordinate convention
# used by the sealed System A adapter.  Every mapping is enumerated below; no
# positional or fuzzy name matching is permitted.

#' Artifact and coordinate contract for a System A SAEM anchor
#'
#' @return An immutable-by-convention list of filenames and exact mappings.
#' @export
sab_system_a_saem_anchor_contract <- function() {
  local_names <- c(
    "u_log_lambda", "u_log_mu_t", "u_log_mu_a", "u_log_p",
    "u_log_alpha_l", "u_pi", "u_eta_rti", "u_eta_pi"
  )
  population_names <- c(
    "mu_log_lambda", "mu_log_mu_t", "mu_log_mu_a", "mu_log_p",
    "mu_log_alpha_l", "mu_u_pi", "mu_u_eta_rti", "mu_u_eta_pi",
    "beta_nelf", "log_omega_lambda", "log_omega_mu_t",
    "log_omega_mu_a", "log_omega_p", "log_omega_alpha_l",
    "log_omega_pi", "log_omega_eta_rti", "log_omega_eta_pi"
  )
  global_names <- c(
    "log_gamma_pop", "log_mu_l_pop", "log_sigma_v", "log_sigma_t"
  )
  population_location_map <- c(
    mu_log_lambda = "u_log_lambda",
    mu_log_mu_t = "u_log_mu_t",
    mu_log_mu_a = "u_log_mu_a",
    mu_log_p = "u_log_p",
    mu_log_alpha_l = "u_log_alpha_l",
    mu_u_pi = "u_pi",
    mu_u_eta_rti = "u_eta_rti",
    mu_u_eta_pi = "u_eta_pi",
    beta_nelf = "beta_treat_nelf(u_eta_pi)"
  )
  global_map <- c(
    log_gamma_pop = "u_log_gamma",
    log_mu_l_pop = "u_log_mu_l",
    log_sigma_v = "u_log_sigma_v",
    log_sigma_t = "u_log_sigma_t"
  )
  omega_map <- c(
    log_omega_lambda = "omega_u_log_lambda",
    log_omega_mu_t = "omega_u_log_mu_t",
    log_omega_mu_a = "omega_u_log_mu_a",
    log_omega_p = "omega_u_log_p",
    log_omega_alpha_l = "omega_u_log_alpha_l",
    log_omega_pi = "omega_u_pi",
    log_omega_eta_rti = "omega_u_eta_rti",
    log_omega_eta_pi = "omega_u_eta_pi"
  )
  structural_zero_omega <- c(
    "omega_u_log_gamma", "omega_u_log_mu_l",
    "omega_u_log_sigma_v", "omega_u_log_sigma_t"
  )
  list(
    schema_version = "sab_system_a_saem_anchor_contract_v1",
    anchor_schema_version = "sab_system_a_saem_anchor_v1",
    target_fingerprint =
      "123fc00a2c0151215f9c550613819d4563a02c5bbe789eee5f2aee2c20844925",
    canonical_hiv_data_sha256 =
      "9397f5f6abb4a934f5885ef8305493c2dd5532b85be62738c7b78f3863d3f1f2",
    expected_patient_count = 115L,
    files = c(
      population = "population_estimates.csv",
      omega = "omega_estimates.csv",
      individuals = "individual_estimates.csv",
      run_config = "run_config.json",
      input_data = "saem_input_data.csv"
    ),
    audit_pattern = "\\.audit\\.json$",
    coordinate_names = list(
      local = local_names,
      population = population_names,
      global = global_names
    ),
    population_location_map = population_location_map,
    global_map = global_map,
    omega_map = omega_map,
    structural_zero_omega = structural_zero_omega,
    individual_map = stats::setNames(local_names, local_names),
    expected_population_parameters = c(
      unname(population_location_map), unname(global_map)
    ),
    expected_omega_parameters = c(
      "omega_u_log_lambda", "omega_u_log_gamma", "omega_u_log_mu_t",
      "omega_u_log_mu_l", "omega_u_log_mu_a", "omega_u_log_p",
      "omega_u_log_alpha_l", "omega_u_pi", "omega_u_eta_rti",
      "omega_u_eta_pi", "omega_u_log_sigma_v", "omega_u_log_sigma_t"
    ),
    expected_individual_latent_columns = c(
      "u_log_lambda", "u_log_gamma", "u_log_mu_t", "u_log_mu_l",
      "u_log_mu_a", "u_log_p", "u_log_alpha_l", "u_pi",
      "u_eta_rti", "u_eta_pi", "u_log_sigma_v", "u_log_sigma_t"
    ),
    expected_input_columns = c(
      "ID", "TIME", "DV", "DV_NCENS", "cens", "YTYPE", "outcome",
      "TREATMENT", "treatment_label", "treat_nelf"
    ),
    ordinary_prediction_reasons = c(
      "ok", "nonfinite_natural_parameter", "invalid_equilibrium",
      "ode_failure", "invalid_ode_output", "invalid_ode_times"
    )
  )
}

#' Read a completed HIV SAEM export as a canonical System A anchor
#'
#' The SAEM estimates are an initialization/proposal anchor only.  Reading
#' them does not alter the sealed System A posterior target.  The exported
#' individual rows are reordered by their explicit IDs into the adapter's
#' canonical patient order.
#'
#' @param run_directory Directory containing the contracted CSV, run-config,
#'   audit, and fit artifacts.
#' @param adapter A full, certified System A adapter (or a test double with the
#'   same callback contract).
#' @param validate_exact If `TRUE`, additionally call the sealed prediction and
#'   likelihood callbacks once for every SAEM patient state and attach the
#'   resulting ledger.  Ordinary target failures are recorded, not hidden.
#'
#' @return A `sab_system_a_saem_anchor` with canonical `eta`, `psi`, and a
#'   115-by-8 `patient_states` matrix.
#' @export
sab_read_system_a_saem_anchor <- function(run_directory, adapter,
                                           validate_exact = FALSE) {
  contract <- sab_system_a_saem_anchor_contract()
  if (!is.logical(validate_exact) || length(validate_exact) != 1L ||
      is.na(validate_exact)) {
    stop("validate_exact must be TRUE or FALSE.", call. = FALSE)
  }
  .sab_saem_validate_adapter(adapter, contract, require_exact = validate_exact)
  if (!is.character(run_directory) || length(run_directory) != 1L ||
      is.na(run_directory) || !nzchar(run_directory) ||
      !dir.exists(run_directory)) {
    stop("run_directory must name an existing SAEM output directory.",
         call. = FALSE)
  }
  run_directory <- normalizePath(run_directory, mustWork = TRUE)
  run_provenance <- .sab_saem_read_run_provenance(run_directory, contract)
  paths <- file.path(run_directory, unname(contract$files))
  names(paths) <- names(contract$files)
  missing <- names(paths)[
    !file.exists(paths) | is.na(file.info(paths)$size) |
      file.info(paths)$size <= 0
  ]
  if (length(missing)) {
    stop(
      "Missing or empty completed SAEM artifact(s): ",
      paste(contract$files[missing], collapse = ", "), ".",
      call. = FALSE
    )
  }

  population <- .sab_saem_read_csv(
    paths[["population"]], c("parameter", "estimate"), "population"
  )
  omega <- .sab_saem_read_csv(
    paths[["omega"]], c("parameter", "variance", "sd"), "omega"
  )
  individuals <- .sab_saem_read_csv(
    paths[["individuals"]], NULL, "individual"
  )
  if (!identical(.sab_saem_sha256(paths[["input_data"]]),
                 contract$canonical_hiv_data_sha256)) {
    stop("SAEM input data do not match the pinned canonical HIV data.",
         call. = FALSE)
  }
  input_data <- .sab_saem_read_csv(
    paths[["input_data"]], contract$expected_input_columns, "input-data"
  )
  .sab_saem_require_parameter_rows(
    population, contract$expected_population_parameters, "population"
  )
  .sab_saem_require_parameter_rows(
    omega, contract$expected_omega_parameters, "omega"
  )
  .sab_saem_validate_numeric_column(population$estimate, "population estimate")
  .sab_saem_validate_numeric_column(omega$variance, "omega variance")
  .sab_saem_validate_numeric_column(omega$sd, "omega sd")
  if (any(omega$variance < 0) || any(omega$sd < 0)) {
    stop("SAEM omega variances and SDs must be non-negative.", call. = FALSE)
  }
  omega_expected_variance <- omega$sd^2
  omega_tolerance <- 64 * .Machine$double.eps * pmax(
    abs(omega$variance), abs(omega_expected_variance),
    .Machine$double.xmin
  )
  if (any(abs(omega$variance - omega_expected_variance) >
          omega_tolerance)) {
    stop("SAEM omega variance and SD columns are inconsistent.",
         call. = FALSE)
  }

  population_values <- stats::setNames(
    as.numeric(population$estimate), population$parameter
  )
  omega_sd <- stats::setNames(as.numeric(omega$sd), omega$parameter)
  mapped_sd <- omega_sd[unname(contract$omega_map)]
  if (any(!is.finite(mapped_sd)) || any(mapped_sd <= 0)) {
    bad <- names(contract$omega_map)[!is.finite(mapped_sd) | mapped_sd <= 0]
    stop(
      "Mapped System A omega SDs must be finite and strictly positive: ",
      paste(bad, collapse = ", "), ".",
      call. = FALSE
    )
  }
  fixed_sd <- omega_sd[contract$structural_zero_omega]
  omega_variance <- stats::setNames(
    as.numeric(omega$variance), omega$parameter
  )
  fixed_variance <- omega_variance[contract$structural_zero_omega]
  if (any(!is.finite(fixed_sd)) || any(!is.finite(fixed_variance)) ||
      any(fixed_sd != 0) || any(fixed_variance != 0)) {
    stop(
      "Structurally non-random SAEM omega rows must be exactly zero: ",
      paste(contract$structural_zero_omega, collapse = ", "), ".",
      call. = FALSE
    )
  }

  eta <- c(
    population_values[unname(contract$population_location_map)],
    log(mapped_sd)
  )
  names(eta) <- c(
    names(contract$population_location_map), names(contract$omega_map)
  )
  eta <- eta[contract$coordinate_names$population]
  psi <- population_values[unname(contract$global_map)]
  names(psi) <- names(contract$global_map)
  psi <- psi[contract$coordinate_names$global]

  .sab_saem_validate_individual_columns(individuals, contract)
  input_ids <- .sab_saem_patient_ids(individuals$ID)
  if (length(input_ids) != contract$expected_patient_count) {
    stop(
      "The SAEM individual artifact must contain exactly ",
      contract$expected_patient_count, " patients; observed ",
      length(input_ids), ".", call. = FALSE
    )
  }
  if (anyDuplicated(input_ids)) {
    stop("The SAEM individual artifact contains duplicate patient IDs.",
         call. = FALSE)
  }
  expected_ids <- adapter$patient_ids
  if (!setequal(input_ids, expected_ids)) {
    stop(
      "SAEM individual IDs do not match the full System A adapter IDs. ",
      "Missing: ", paste(setdiff(expected_ids, input_ids), collapse = ","),
      "; unexpected: ",
      paste(setdiff(input_ids, expected_ids), collapse = ","), ".",
      call. = FALSE
    )
  }
  data_ids <- unique(.sab_saem_patient_ids(input_data$ID))
  if (length(data_ids) != contract$expected_patient_count ||
      !setequal(data_ids, expected_ids)) {
    stop(
      "SAEM input-data IDs do not match the full System A adapter IDs.",
      call. = FALSE
    )
  }
  if (!identical(run_provenance$audit$n_rows, nrow(input_data)) ||
      !identical(run_provenance$audit$n_patients,
                 length(data_ids))) {
    stop("SAEM audit counts do not match saem_input_data.csv.",
         call. = FALSE)
  }
  row_order_validation <- if (inherits(adapter, "sab_system_a_adapter")) {
    .sab_saem_validate_fit_binding(
      run_directory, run_provenance, population, omega, individuals,
      input_ids, contract
    )
  } else {
    list(
      schema_version = "sab_system_a_saem_row_order_validation_v1",
      validated = TRUE,
      method = "synthetic_test_fixture",
      fit_path = NA_character_,
      fit_sha256 = NA_character_,
      patient_ids = input_ids,
      max_abs_map_csv_error = 0,
      max_abs_fixed_csv_error = 0,
      max_abs_omega_csv_error = 0
    )
  }
  patient_rows <- match(expected_ids, input_ids)
  patient_states <- as.matrix(individuals[
    patient_rows, unname(contract$individual_map), drop = FALSE
  ])
  storage.mode(patient_states) <- "double"
  colnames(patient_states) <- names(contract$individual_map)
  rownames(patient_states) <- expected_ids
  if (any(!is.finite(patient_states))) {
    stop("SAEM patient states must all be finite.", call. = FALSE)
  }

  # SAEMIX repeats fixed shared coordinates in map.psi.  Checking those
  # columns guards against combining population and individual files from
  # different runs while allowing only negligible CSV round-off.
  fixed_columns <- unname(contract$global_map)
  fixed_matrix <- as.matrix(individuals[
    patient_rows, fixed_columns, drop = FALSE
  ])
  storage.mode(fixed_matrix) <- "double"
  fixed_reference <- population_values[fixed_columns]
  fixed_error <- sweep(fixed_matrix, 2L, fixed_reference, "-")
  fixed_tolerance <- 1e-10 * pmax(1, abs(fixed_reference))
  if (any(!is.finite(fixed_matrix)) ||
      any(abs(fixed_error) > matrix(
        fixed_tolerance, nrow = nrow(fixed_matrix),
        ncol = ncol(fixed_matrix), byrow = TRUE
      ))) {
    stop(
      "Individual fixed shared coordinates disagree with the SAEM ",
      "population artifact.", call. = FALSE
    )
  }

  if (!isTRUE(adapter$eta_in_domain(eta))) {
    stop("The converted SAEM eta is outside the System A target domain.",
         call. = FALSE)
  }
  if (!isTRUE(adapter$psi_in_domain(psi))) {
    stop("The converted SAEM psi is outside the System A target domain.",
         call. = FALSE)
  }
  population_log_density <- vapply(
    expected_ids,
    function(patient_id) {
      value <- adapter$log_population_density(
        patient_id, patient_states[patient_id, ], eta
      )
      if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
        stop(
          "Non-finite System A population density for SAEM patient ",
          patient_id, ".", call. = FALSE
        )
      }
      as.numeric(value)
    },
    numeric(1L)
  )

  source_paths <- c(paths, audit = run_provenance$audit$path)
  anchor <- structure(
    list(
      schema_version = contract$anchor_schema_version,
      target_fingerprint = contract$target_fingerprint,
      source = list(
        run_directory = run_directory,
        artifact_files = source_paths,
        artifact_sha256 = vapply(
          source_paths, .sab_saem_sha256, character(1L)
        ),
        input_patient_ids = input_ids,
        reordered_to_canonical_ids = !identical(input_ids, expected_ids),
        audit = run_provenance$audit,
        run_config = run_provenance$run_config,
        row_order_validation = row_order_validation
      ),
      coordinate_names = contract$coordinate_names,
      patient_ids = expected_ids,
      eta = eta,
      psi = psi,
      patient_states = patient_states,
      validation = list(
        population_log_density = population_log_density,
        exact_target = NULL
      )
    ),
    class = c("sab_system_a_saem_anchor", "list")
  )
  sab_validate_system_a_saem_anchor(anchor, adapter)
  if (validate_exact) {
    anchor$validation$exact_target <-
      sab_validate_system_a_saem_anchor_exact(anchor, adapter)
  }
  anchor
}

#' Validate a canonical SAEM anchor against its System A adapter
#'
#' @param anchor Object returned by [sab_read_system_a_saem_anchor()].
#' @param adapter Full System A adapter used to interpret the anchor.
#'
#' @return `anchor`, invisibly.
#' @export
sab_validate_system_a_saem_anchor <- function(anchor, adapter) {
  contract <- sab_system_a_saem_anchor_contract()
  .sab_saem_validate_adapter(adapter, contract, require_exact = FALSE)
  .sab_saem_validate_stored_provenance(
    anchor$source, contract,
    require_fit_binding = inherits(adapter, "sab_system_a_adapter")
  )
  valid <- inherits(anchor, "sab_system_a_saem_anchor") &&
    is.list(anchor) &&
    identical(anchor$schema_version, contract$anchor_schema_version) &&
    identical(anchor$target_fingerprint, contract$target_fingerprint) &&
    identical(anchor$target_fingerprint, adapter$target_fingerprint) &&
    identical(anchor$coordinate_names, contract$coordinate_names) &&
    identical(anchor$patient_ids, adapter$patient_ids) &&
    is.numeric(anchor$eta) &&
    identical(names(anchor$eta), contract$coordinate_names$population) &&
    all(is.finite(anchor$eta)) &&
    is.numeric(anchor$psi) &&
    identical(names(anchor$psi), contract$coordinate_names$global) &&
    all(is.finite(anchor$psi)) &&
    is.matrix(anchor$patient_states) && is.numeric(anchor$patient_states) &&
    identical(dim(anchor$patient_states),
              c(contract$expected_patient_count,
                length(contract$coordinate_names$local))) &&
    identical(rownames(anchor$patient_states), adapter$patient_ids) &&
    identical(colnames(anchor$patient_states),
              contract$coordinate_names$local) &&
    all(is.finite(anchor$patient_states)) &&
    is.list(anchor$source) &&
    is.list(anchor$validation) &&
    is.numeric(anchor$validation$population_log_density) &&
    identical(names(anchor$validation$population_log_density),
              adapter$patient_ids) &&
    all(is.finite(anchor$validation$population_log_density))
  if (!isTRUE(valid)) {
    stop("Malformed canonical System A SAEM anchor.", call. = FALSE)
  }
  if (!isTRUE(adapter$eta_in_domain(anchor$eta)) ||
      !isTRUE(adapter$psi_in_domain(anchor$psi))) {
    stop("Canonical System A SAEM anchor is outside the target domain.",
         call. = FALSE)
  }
  recalculated <- vapply(
    adapter$patient_ids,
    function(patient_id) {
      value <- adapter$log_population_density(
        patient_id, anchor$patient_states[patient_id, ], anchor$eta
      )
      if (!is.numeric(value) || length(value) != 1L || !is.finite(value)) {
        stop(
          "Non-finite System A population density for canonical SAEM ",
          "patient ", patient_id, ".", call. = FALSE
        )
      }
      as.numeric(value)
    },
    numeric(1L)
  )
  if (!isTRUE(all.equal(
      unname(recalculated),
      unname(anchor$validation$population_log_density),
      tolerance = 1e-12, check.attributes = FALSE
  ))) {
    stop("Stored SAEM-anchor population densities are stale or malformed.",
         call. = FALSE)
  }
  invisible(anchor)
}

#' Evaluate every SAEM patient state through the sealed exact target
#'
#' This diagnostic deliberately calls `solve_prediction()` first and then
#' `loglik_from_prediction()`, including for ordinary failed prediction
#' objects.  It never substitutes the SAEM/LSODA likelihood.
#'
#' @param anchor Canonical System A SAEM anchor.
#' @param adapter The same full System A adapter used for conversion.
#'
#' @return A list containing a patient-level ledger, aggregate counts, and a
#'   `passed` flag.  `completed_without_callback_error` distinguishes ordinary
#'   exact-target rejection of an SAEM MAP state from a broken callback or
#'   callback contract; neither is ever converted into a finite likelihood.
#' @export
sab_validate_system_a_saem_anchor_exact <- function(anchor, adapter) {
  contract <- sab_system_a_saem_anchor_contract()
  .sab_saem_validate_adapter(adapter, contract, require_exact = TRUE)
  sab_validate_system_a_saem_anchor(anchor, adapter)

  rows <- lapply(anchor$patient_ids, function(patient_id) {
    x <- anchor$patient_states[patient_id, ]
    prediction_error <- NULL
    prediction <- tryCatch(
      adapter$solve_prediction(patient_id, x, anchor$psi),
      error = function(error) {
        prediction_error <<- conditionMessage(error)
        NULL
      }
    )
    if (!is.null(prediction_error)) {
      return(.sab_saem_exact_row(
        patient_id, anchor$validation$population_log_density[[patient_id]],
        ode_integrations = NA_integer_, prediction_ok = FALSE,
        prediction_reason = paste0("prediction_callback_error: ",
                                   prediction_error),
        likelihood_calls = 0L, loglik = NA_real_, callback_error = TRUE
      ))
    }
    prediction_check <- .sab_saem_prediction_summary(
      prediction, contract$ordinary_prediction_reasons
    )
    likelihood_error <- NULL
    loglik <- tryCatch(
      adapter$loglik_from_prediction(prediction, anchor$psi),
      error = function(error) {
        likelihood_error <<- conditionMessage(error)
        NA_real_
      }
    )
    if (!is.null(likelihood_error)) {
      prediction_check$reason <- paste0(
        prediction_check$reason, "; likelihood_callback_error: ",
        likelihood_error
      )
    }
    invalid_loglik <- !is.numeric(loglik) || length(loglik) != 1L ||
      is.na(loglik) || is.nan(loglik) ||
      (is.infinite(loglik) && loglik > 0)
    .sab_saem_exact_row(
      patient_id, anchor$validation$population_log_density[[patient_id]],
      ode_integrations = prediction_check$ode_integrations,
      prediction_ok = prediction_check$ok,
      prediction_reason = prediction_check$reason,
      likelihood_calls = 1L, loglik = loglik,
      callback_error = prediction_check$contract_error ||
        !is.null(likelihood_error) || invalid_loglik
    )
  })
  ledger <- do.call(rbind, rows)
  rownames(ledger) <- NULL
  passed <- all(ledger$passed)
  completed_without_callback_error <- !any(ledger$callback_error)
  structure(
    list(
      schema_version = "sab_system_a_saem_exact_validation_v1",
      target_fingerprint = anchor$target_fingerprint,
      ledger = ledger,
      summary = data.frame(
        patients = nrow(ledger),
        prediction_calls = sum(ledger$prediction_calls),
        ode_integrations = if (anyNA(ledger$ode_integrations)) {
          NA_integer_
        } else {
          as.integer(sum(ledger$ode_integrations))
        },
        unknown_ode_counts = sum(is.na(ledger$ode_integrations)),
        prediction_failures = sum(!ledger$prediction_ok),
        likelihood_calls = sum(ledger$likelihood_calls),
        nonfinite_loglik = sum(!ledger$finite_loglik),
        ordinary_target_rejections =
          sum(ledger$ordinary_target_rejection),
        callback_errors = sum(ledger$callback_error),
        completed_without_callback_error =
          completed_without_callback_error,
        passed = passed
      ),
      completed_without_callback_error = completed_without_callback_error,
      passed = passed
    ),
    class = c("sab_system_a_saem_exact_validation", "list")
  )
}

.sab_saem_read_run_provenance <- function(run_directory, contract) {
  if (!requireNamespace("jsonlite", quietly = TRUE)) {
    stop("jsonlite is required to verify SAEM run provenance.",
         call. = FALSE)
  }
  audit_paths <- list.files(
    run_directory, pattern = contract$audit_pattern,
    full.names = TRUE, recursive = FALSE
  )
  audit_paths <- audit_paths[
    file.exists(audit_paths) & !is.na(file.info(audit_paths)$size) &
      file.info(audit_paths)$size > 0
  ]
  if (length(audit_paths) != 1L) {
    stop(
      "A completed SAEM run directory must contain exactly one non-empty ",
      "*.audit.json file; observed ", length(audit_paths), ".",
      call. = FALSE
    )
  }
  config_path <- file.path(
    run_directory, contract$files[["run_config"]]
  )
  if (!file.exists(config_path) || is.na(file.info(config_path)$size) ||
      file.info(config_path)$size <= 0) {
    stop("Missing or empty completed SAEM run_config.json.", call. = FALSE)
  }
  audit <- .sab_saem_read_json(audit_paths[[1L]], "audit")
  config <- .sab_saem_read_json(config_path, "run-config")

  status <- .sab_saem_json_character(audit$status, "audit status")
  if (!identical(status, "completed")) {
    stop("SAEM audit status must be exactly `completed`.", call. = FALSE)
  }
  started_at <- .sab_saem_json_character(
    audit$started_at, "audit started_at"
  )
  updated_at <- .sab_saem_json_character(
    audit$updated_at, "audit updated_at"
  )
  events <- audit$events
  if (!is.list(events) || !length(events) ||
      any(!vapply(events, is.list, logical(1L)))) {
    stop("SAEM audit events are malformed.", call. = FALSE)
  }
  event_names <- vapply(events, function(event) {
    tryCatch(
      .sab_saem_json_character(event$event, "audit event"),
      error = function(error) NA_character_
    )
  }, character(1L))
  if (anyNA(event_names) ||
      sum(event_names == "hiv_benchmark_loaded") != 1L ||
      sum(event_names == "saem_fit_completed") != 1L) {
    stop(
      "SAEM audit must contain exactly one HIV-data load and one completed ",
      "fit event.", call. = FALSE
    )
  }
  data_event <- events[[which(event_names == "hiv_benchmark_loaded")]]
  fit_event <- events[[which(event_names == "saem_fit_completed")]]
  n_rows <- .sab_saem_json_integer(
    data_event$n_rows, "audit n_rows", length = 1L, minimum = 1L
  )
  n_patients <- .sab_saem_json_integer(
    data_event$n_patients, "audit n_patients", length = 1L, minimum = 1L
  )
  if (n_patients != contract$expected_patient_count) {
    stop(
      "SAEM audit must report exactly ", contract$expected_patient_count,
      " patients.", call. = FALSE
    )
  }
  ode_failure_count <- .sab_saem_json_integer(
    fit_event$ode_failure_count, "audit ODE failure count",
    length = 1L, minimum = 0L
  )

  run_tag <- .sab_saem_json_character(config$run_tag, "run_config run_tag")
  expected_audit_name <- paste0(run_tag, ".audit.json")
  if (!identical(basename(audit_paths[[1L]]), expected_audit_name)) {
    stop("SAEM audit filename does not match run_config run_tag.",
         call. = FALSE)
  }
  configured_run_directory <- .sab_saem_json_character(
    config$run_dir, "run_config run_dir"
  )
  if (!dir.exists(configured_run_directory) ||
      !identical(normalizePath(configured_run_directory, mustWork = TRUE),
                 run_directory)) {
    stop("run_config run_dir does not identify the supplied run directory.",
         call. = FALSE)
  }
  hiv_data_path <- .sab_saem_json_character(
    config$hiv_data_path, "run_config hiv_data_path"
  )
  if (!file.exists(hiv_data_path) || is.na(file.info(hiv_data_path)$size) ||
      file.info(hiv_data_path)$size <= 0) {
    stop("run_config HIV source data do not exist or are empty.",
         call. = FALSE)
  }
  hiv_data_path <- normalizePath(hiv_data_path, mustWork = TRUE)
  hiv_data_sha256 <- .sab_saem_sha256(hiv_data_path)
  if (!identical(hiv_data_sha256, contract$canonical_hiv_data_sha256)) {
    stop("run_config HIV source data fail the pinned canonical hash.",
         call. = FALSE)
  }

  fit <- config$fit_control
  if (!is.list(fit)) stop("run_config fit_control is malformed.",
                          call. = FALSE)
  saem_iterations <- .sab_saem_json_integer(
    fit[["nbiter.saemix"]], "fit_control nbiter.saemix",
    length = 2L, minimum = 0L
  )
  if (saem_iterations[[1L]] < 1L) {
    stop("The first SAEM phase must contain at least one iteration.",
         call. = FALSE)
  }
  mcmc_iterations <- .sab_saem_json_integer(
    fit[["nbiter.mcmc"]], "fit_control nbiter.mcmc",
    length = 4L, minimum = 0L
  )
  fit_control <- list(
    k1 = saem_iterations[[1L]],
    k2 = saem_iterations[[2L]],
    nbiter_mcmc = mcmc_iterations,
    nbiter_sa = .sab_saem_json_integer(
      fit[["nbiter.sa"]], "fit_control nbiter.sa",
      length = 1L, minimum = 0L
    ),
    nbiter_burn = .sab_saem_json_integer(
      fit[["nbiter.burn"]], "fit_control nbiter.burn",
      length = 1L, minimum = 0L
    ),
    nbiter_map = .sab_saem_json_integer(
      fit[["nbiter.map"]], "fit_control nbiter.map",
      length = 1L, minimum = 0L
    ),
    map = .sab_saem_json_logical(fit$map, "fit_control map"),
    fim = .sab_saem_json_logical(fit$fim, "fit_control fim"),
    seed = .sab_saem_json_integer(
      fit$seed, "fit_control seed", length = 1L, minimum = 0L
    ),
    nb_chains = .sab_saem_json_integer(
      fit[["nb.chains"]], "fit_control nb.chains",
      length = 1L, minimum = 1L
    )
  )
  ode <- config$ode_control
  if (!is.list(ode)) stop("run_config ode_control is malformed.",
                          call. = FALSE)
  ode_control <- list(
    method = .sab_saem_json_character(ode$method, "ode_control method"),
    rtol = .sab_saem_json_positive_number(ode$rtol, "ode_control rtol"),
    atol = .sab_saem_json_positive_number(ode$atol, "ode_control atol"),
    maxsteps = .sab_saem_json_integer(
      ode$maxsteps, "ode_control maxsteps", length = 1L, minimum = 1L
    ),
    time_eps = .sab_saem_json_positive_number(
      ode$time_eps, "ode_control time_eps"
    )
  )
  list(
    audit = list(
      path = normalizePath(audit_paths[[1L]], mustWork = TRUE),
      sha256 = .sab_saem_sha256(audit_paths[[1L]]),
      status = status,
      started_at = started_at,
      updated_at = updated_at,
      n_rows = n_rows,
      n_patients = n_patients,
      ode_failure_count = ode_failure_count
    ),
    run_config = list(
      path = normalizePath(config_path, mustWork = TRUE),
      sha256 = .sab_saem_sha256(config_path),
      run_tag = run_tag,
      hiv_data_path = hiv_data_path,
      hiv_data_sha256 = hiv_data_sha256,
      fit_control = fit_control,
      ode_control = ode_control
    )
  )
}

.sab_saem_validate_fit_binding <- function(
    run_directory, run_provenance, population, omega, individuals,
    input_ids, contract) {
  fit_paths <- list.files(
    run_directory, pattern = "\\.fit\\.rds$", full.names = TRUE,
    recursive = FALSE
  )
  fit_paths <- fit_paths[
    file.exists(fit_paths) & !is.na(file.info(fit_paths)$size) &
      file.info(fit_paths)$size > 0
  ]
  if (length(fit_paths) != 1L) {
    stop(
      "A real System A anchor requires exactly one non-empty SAEM fit RDS ",
      "for patient-row validation; observed ", length(fit_paths), ".",
      call. = FALSE
    )
  }
  expected_name <- paste0(
    run_provenance$run_config$run_tag, ".fit.rds"
  )
  if (!identical(basename(fit_paths[[1L]]), expected_name)) {
    stop("SAEM fit filename does not match run_config run_tag.",
         call. = FALSE)
  }
  fit <- tryCatch(
    readRDS(fit_paths[[1L]]),
    error = function(error) {
      stop("Could not read SAEM fit for row-order validation: ",
           conditionMessage(error), call. = FALSE)
    }
  )
  if (!isS4(fit) || !all(c("results", "data") %in% methods::slotNames(fit))) {
    stop("SAEM fit does not expose results and data slots.", call. = FALSE)
  }
  results <- methods::slot(fit, "results")
  data_object <- methods::slot(fit, "data")
  if (!isS4(results) || !isS4(data_object) ||
      !all(c("map.psi", "fixed.effects", "name.fixed", "omega") %in%
           methods::slotNames(results)) ||
      !"data" %in% methods::slotNames(data_object)) {
    stop("SAEM fit omits fields required for artifact binding.",
         call. = FALSE)
  }

  map <- as.data.frame(methods::slot(results, "map.psi"),
                       check.names = FALSE)
  fit_data <- as.data.frame(methods::slot(data_object, "data"),
                            check.names = FALSE)
  if (nrow(map) != contract$expected_patient_count ||
      !identical(names(map), contract$expected_individual_latent_columns) ||
      !identical(rownames(map),
                 as.character(seq_len(contract$expected_patient_count))) ||
      !all(c("index", "ID") %in% names(fit_data))) {
    stop("SAEM fit has no canonical internal patient-row index.",
         call. = FALSE)
  }
  index <- fit_data$index
  if ((!is.numeric(index) && !is.integer(index)) || any(!is.finite(index)) ||
      any(index != floor(index))) {
    stop("SAEM fit patient indices are malformed.", call. = FALSE)
  }
  index_id <- unique(fit_data[c("index", "ID")])
  index_id <- index_id[order(index_id$index), , drop = FALSE]
  if (nrow(index_id) != contract$expected_patient_count ||
      !identical(as.integer(index_id$index),
                 seq_len(contract$expected_patient_count))) {
    stop("SAEM fit patient indices are not one-to-one and contiguous.",
         call. = FALSE)
  }
  internal_ids <- .sab_saem_patient_ids(index_id$ID)
  if (!identical(internal_ids, input_ids)) {
    stop(
      "SAEM fit internal index-to-ID map disagrees with individual CSV rows.",
      call. = FALSE
    )
  }

  map_matrix <- as.matrix(map)
  csv_map_matrix <- as.matrix(individuals[
    , contract$expected_individual_latent_columns, drop = FALSE
  ])
  storage.mode(map_matrix) <- "double"
  storage.mode(csv_map_matrix) <- "double"
  map_error <- .sab_saem_max_abs_error(
    map_matrix, csv_map_matrix, "SAEM patient MAP"
  )

  fixed <- as.numeric(methods::slot(results, "fixed.effects"))
  fixed_names <- as.character(methods::slot(results, "name.fixed"))
  if (length(fixed) != length(fixed_names) || anyDuplicated(fixed_names) ||
      any(!is.finite(fixed)) ||
      !setequal(fixed_names, contract$expected_population_parameters)) {
    stop("SAEM fit fixed effects do not match the population CSV contract.",
         call. = FALSE)
  }
  fit_fixed <- stats::setNames(fixed, fixed_names)[population$parameter]
  fixed_error <- .sab_saem_max_abs_error(
    fit_fixed, population$estimate, "SAEM fixed effect"
  )

  fit_omega <- as.matrix(methods::slot(results, "omega"))
  if (!is.numeric(fit_omega) || nrow(fit_omega) != ncol(fit_omega) ||
      nrow(fit_omega) != length(contract$expected_omega_parameters) ||
      any(!is.finite(fit_omega))) {
    stop("SAEM fit omega matrix does not match the omega CSV contract.",
         call. = FALSE)
  }
  fit_variance <- stats::setNames(
    diag(fit_omega), contract$expected_omega_parameters
  )[omega$parameter]
  omega_error <- .sab_saem_max_abs_error(
    fit_variance, omega$variance, "SAEM omega"
  )

  list(
    schema_version = "sab_system_a_saem_row_order_validation_v1",
    validated = TRUE,
    method = "saemix_internal_index_binding",
    fit_path = normalizePath(fit_paths[[1L]], mustWork = TRUE),
    fit_sha256 = .sab_saem_sha256(fit_paths[[1L]]),
    patient_ids = internal_ids,
    max_abs_map_csv_error = map_error,
    max_abs_fixed_csv_error = fixed_error,
    max_abs_omega_csv_error = omega_error
  )
}

.sab_saem_max_abs_error <- function(observed, expected, label) {
  if (!is.numeric(observed) || !is.numeric(expected) ||
      !identical(dim(observed), dim(expected)) ||
      length(observed) != length(expected) ||
      any(!is.finite(observed)) || any(!is.finite(expected))) {
    stop(label, " artifact comparison is malformed.", call. = FALSE)
  }
  scale <- pmax(1, abs(observed), abs(expected))
  difference <- abs(observed - expected)
  if (any(difference > 128 * .Machine$double.eps * scale)) {
    stop(label, " values disagree between fit and CSV artifact.",
         call. = FALSE)
  }
  if (length(difference)) max(difference) else 0
}

.sab_saem_read_json <- function(path, label) {
  value <- tryCatch(
    jsonlite::fromJSON(path, simplifyVector = FALSE),
    error = function(error) {
      stop("Could not parse SAEM ", label, " JSON: ",
           conditionMessage(error), call. = FALSE)
    }
  )
  if (!is.list(value)) {
    stop("SAEM ", label, " JSON must contain one object.", call. = FALSE)
  }
  value
}

.sab_saem_json_character <- function(value, label) {
  if (!is.character(value) || length(value) != 1L || is.na(value) ||
      !nzchar(value)) {
    stop(label, " must be one non-empty string.", call. = FALSE)
  }
  value
}

.sab_saem_json_logical <- function(value, label) {
  if (!is.logical(value) || length(value) != 1L || is.na(value)) {
    stop(label, " must be one JSON boolean.", call. = FALSE)
  }
  value
}

.sab_saem_json_integer <- function(value, label, length, minimum) {
  value <- unlist(value, recursive = FALSE, use.names = FALSE)
  if ((!is.numeric(value) && !is.integer(value)) ||
      base::length(value) != length || any(!is.finite(value)) ||
      any(value != floor(value)) || any(value < minimum) ||
      any(value > .Machine$integer.max)) {
    stop(label, " must contain exactly ", length,
         " integer value(s) >= ", minimum, ".", call. = FALSE)
  }
  as.integer(value)
}

.sab_saem_json_positive_number <- function(value, label) {
  value <- unlist(value, recursive = FALSE, use.names = FALSE)
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value <= 0) {
    stop(label, " must be one finite positive number.", call. = FALSE)
  }
  as.numeric(value)
}

.sab_saem_validate_stored_provenance <- function(
    source, contract, require_fit_binding) {
  expected_path_names <- c(names(contract$files), "audit")
  if (!is.list(source) ||
      !is.character(source$run_directory) ||
      length(source$run_directory) != 1L ||
      is.na(source$run_directory) || !dir.exists(source$run_directory) ||
      !is.character(source$artifact_files) ||
      !identical(names(source$artifact_files), expected_path_names) ||
      anyNA(source$artifact_files) || any(!nzchar(source$artifact_files)) ||
      !is.character(source$artifact_sha256) ||
      !identical(names(source$artifact_sha256), expected_path_names) ||
      anyNA(source$artifact_sha256) ||
      any(!grepl("^[0-9a-f]{64}$", source$artifact_sha256)) ||
      !is.character(source$input_patient_ids) ||
      length(source$input_patient_ids) != contract$expected_patient_count ||
      anyNA(source$input_patient_ids) ||
      any(!nzchar(source$input_patient_ids)) ||
      anyDuplicated(source$input_patient_ids) ||
      !is.logical(source$reordered_to_canonical_ids) ||
      length(source$reordered_to_canonical_ids) != 1L ||
      is.na(source$reordered_to_canonical_ids) ||
      !is.list(source$audit) ||
      !is.character(source$audit$path) ||
      length(source$audit$path) != 1L || is.na(source$audit$path) ||
      !nzchar(source$audit$path) || !is.list(source$run_config) ||
      !is.list(source$row_order_validation)) {
    stop("Malformed System A SAEM run provenance.", call. = FALSE)
  }
  run_directory <- normalizePath(source$run_directory, mustWork = TRUE)
  expected_paths <- c(
    file.path(run_directory, unname(contract$files)),
    audit = source$audit$path
  )
  names(expected_paths) <- expected_path_names
  if (any(!file.exists(source$artifact_files)) ||
      !identical(
        unname(normalizePath(source$artifact_files, mustWork = TRUE)),
        unname(normalizePath(expected_paths, mustWork = TRUE))
      )) {
    stop("Stored System A SAEM artifact paths are missing or inconsistent.",
         call. = FALSE)
  }
  observed_hashes <- vapply(
    source$artifact_files, .sab_saem_sha256, character(1L)
  )
  if (!identical(unname(observed_hashes),
                 unname(source$artifact_sha256))) {
    stop("Stored System A SAEM artifact hashes no longer match their files.",
         call. = FALSE)
  }

  audit <- source$audit
  audit_valid <- is.character(audit$path) && length(audit$path) == 1L &&
    identical(audit$path, source$artifact_files[["audit"]]) &&
    is.character(audit$sha256) && length(audit$sha256) == 1L &&
    identical(audit$sha256, source$artifact_sha256[["audit"]]) &&
    identical(audit$status, "completed") &&
    is.character(audit$started_at) && length(audit$started_at) == 1L &&
    !is.na(audit$started_at) && nzchar(audit$started_at) &&
    is.character(audit$updated_at) && length(audit$updated_at) == 1L &&
    !is.na(audit$updated_at) && nzchar(audit$updated_at) &&
    is.integer(audit$n_rows) && length(audit$n_rows) == 1L &&
    !is.na(audit$n_rows) && audit$n_rows >= 1L &&
    identical(audit$n_patients, contract$expected_patient_count) &&
    is.integer(audit$ode_failure_count) &&
    length(audit$ode_failure_count) == 1L &&
    !is.na(audit$ode_failure_count) && audit$ode_failure_count >= 0L
  if (!isTRUE(audit_valid)) {
    stop("Malformed or incomplete System A SAEM audit provenance.",
         call. = FALSE)
  }

  run_config <- source$run_config
  fit <- run_config$fit_control
  ode <- run_config$ode_control
  fit_valid <- is.list(fit) && identical(
    names(fit),
    c(
      "k1", "k2", "nbiter_mcmc", "nbiter_sa", "nbiter_burn",
      "nbiter_map", "map", "fim", "seed", "nb_chains"
    )
  ) &&
    is.integer(fit$k1) && length(fit$k1) == 1L && fit$k1 >= 1L &&
    is.integer(fit$k2) && length(fit$k2) == 1L && fit$k2 >= 0L &&
    is.integer(fit$nbiter_mcmc) && length(fit$nbiter_mcmc) == 4L &&
    all(fit$nbiter_mcmc >= 0L) &&
    all(vapply(
      fit[c("nbiter_sa", "nbiter_burn", "nbiter_map", "seed")],
      function(value) {
        is.integer(value) && length(value) == 1L &&
          !is.na(value) && value >= 0L
      }, logical(1L)
    )) &&
    is.integer(fit$nb_chains) && length(fit$nb_chains) == 1L &&
    fit$nb_chains >= 1L &&
    is.logical(fit$map) && length(fit$map) == 1L && !is.na(fit$map) &&
    is.logical(fit$fim) && length(fit$fim) == 1L && !is.na(fit$fim)
  ode_valid <- is.list(ode) && identical(
    names(ode), c("method", "rtol", "atol", "maxsteps", "time_eps")
  ) &&
    is.character(ode$method) && length(ode$method) == 1L &&
    !is.na(ode$method) && nzchar(ode$method) &&
    all(vapply(ode[c("rtol", "atol", "time_eps")], function(value) {
      is.numeric(value) && length(value) == 1L &&
        is.finite(value) && value > 0
    }, logical(1L))) &&
    is.integer(ode$maxsteps) && length(ode$maxsteps) == 1L &&
    ode$maxsteps >= 1L
  run_valid <- is.character(run_config$path) &&
    length(run_config$path) == 1L &&
    identical(run_config$path, source$artifact_files[["run_config"]]) &&
    is.character(run_config$sha256) && length(run_config$sha256) == 1L &&
    identical(run_config$sha256,
              source$artifact_sha256[["run_config"]]) &&
    is.character(run_config$run_tag) && length(run_config$run_tag) == 1L &&
    !is.na(run_config$run_tag) && nzchar(run_config$run_tag) &&
    identical(basename(audit$path),
              paste0(run_config$run_tag, ".audit.json")) &&
    is.character(run_config$hiv_data_path) &&
    length(run_config$hiv_data_path) == 1L &&
    file.exists(run_config$hiv_data_path) &&
    is.character(run_config$hiv_data_sha256) &&
    length(run_config$hiv_data_sha256) == 1L &&
    identical(run_config$hiv_data_sha256,
              contract$canonical_hiv_data_sha256) &&
    identical(source$artifact_sha256[["input_data"]],
              contract$canonical_hiv_data_sha256) &&
    identical(
      .sab_saem_sha256(run_config$hiv_data_path),
      run_config$hiv_data_sha256
    ) && isTRUE(fit_valid) && isTRUE(ode_valid)
  if (!isTRUE(run_valid)) {
    stop("Malformed System A SAEM run-config provenance.", call. = FALSE)
  }

  row_order <- source$row_order_validation
  row_order_valid <- identical(
    row_order$schema_version,
    "sab_system_a_saem_row_order_validation_v1"
  ) && isTRUE(row_order$validated) &&
    is.character(row_order$method) && length(row_order$method) == 1L &&
    !is.na(row_order$method) && nzchar(row_order$method) &&
    is.character(row_order$patient_ids) &&
    identical(row_order$patient_ids, source$input_patient_ids) &&
    all(vapply(
      row_order[c(
        "max_abs_map_csv_error", "max_abs_fixed_csv_error",
        "max_abs_omega_csv_error"
      )],
      function(value) {
        is.numeric(value) && length(value) == 1L &&
          is.finite(value) && value >= 0
      }, logical(1L)
    ))
  if (isTRUE(require_fit_binding)) {
    row_order_valid <- row_order_valid &&
      identical(row_order$method, "saemix_internal_index_binding") &&
      is.character(row_order$fit_path) &&
      length(row_order$fit_path) == 1L && !is.na(row_order$fit_path) &&
      file.exists(row_order$fit_path) &&
      is.character(row_order$fit_sha256) &&
      length(row_order$fit_sha256) == 1L &&
      grepl("^[0-9a-f]{64}$", row_order$fit_sha256) &&
      identical(.sab_saem_sha256(row_order$fit_path),
                row_order$fit_sha256)
  } else {
    row_order_valid <- row_order_valid &&
      identical(row_order$method, "synthetic_test_fixture")
  }
  if (!isTRUE(row_order_valid)) {
    stop("Malformed or unverified System A SAEM patient-row binding.",
         call. = FALSE)
  }
  invisible(source)
}

.sab_saem_validate_adapter <- function(adapter, contract, require_exact) {
  callbacks <- c(
    "eta_in_domain", "psi_in_domain", "log_population_density"
  )
  if (require_exact) {
    callbacks <- c(callbacks, "solve_prediction", "loglik_from_prediction")
  }
  real_adapter <- inherits(adapter, "sab_system_a_adapter")
  mock_adapter <- inherits(adapter, "sab_saem_anchor_mock_adapter")
  valid <- (real_adapter || mock_adapter) && is.list(adapter) &&
    identical(adapter$target_fingerprint, contract$target_fingerprint) &&
    is.character(adapter$patient_ids) &&
    length(adapter$patient_ids) == contract$expected_patient_count &&
    !anyNA(adapter$patient_ids) && !any(!nzchar(adapter$patient_ids)) &&
    !anyDuplicated(adapter$patient_ids) &&
    identical(adapter$coordinate_names$local,
              contract$coordinate_names$local) &&
    identical(adapter$coordinate_names$population,
              contract$coordinate_names$population) &&
    identical(adapter$coordinate_names$global,
              contract$coordinate_names$global) &&
    all(vapply(callbacks, function(name) {
      is.function(adapter[[name]])
    }, logical(1L)))
  if (!isTRUE(valid)) {
    stop("A full canonical System A adapter is required for SAEM conversion.",
         call. = FALSE)
  }
  if (real_adapter &&
      !exists("sab_validate_system_a_adapter", mode = "function")) {
    stop("Load the certified System A adapter module before conversion.",
         call. = FALSE)
  }
  if (real_adapter) {
    sab_validate_system_a_adapter(adapter)
  }
  invisible(adapter)
}

.sab_saem_read_csv <- function(path, exact_columns, label) {
  value <- tryCatch(
    utils::read.csv(
      path, stringsAsFactors = FALSE, check.names = FALSE,
      na.strings = c("NA", "NaN", "Inf", "-Inf")
    ),
    error = function(error) {
      stop("Could not parse SAEM ", label, " artifact: ",
           conditionMessage(error), call. = FALSE)
    }
  )
  if (!is.data.frame(value) || nrow(value) < 1L ||
      anyNA(names(value)) || any(!nzchar(names(value))) ||
      anyDuplicated(names(value))) {
    stop("Malformed SAEM ", label, " artifact.", call. = FALSE)
  }
  if (!is.null(exact_columns) && !identical(names(value), exact_columns)) {
    stop(
      "SAEM ", label, " artifact columns must be exactly: ",
      paste(exact_columns, collapse = ", "), ".",
      call. = FALSE
    )
  }
  value
}

.sab_saem_require_parameter_rows <- function(value, expected, label) {
  parameter <- value$parameter
  if (!is.character(parameter) || anyNA(parameter) ||
      any(!nzchar(parameter)) || anyDuplicated(parameter)) {
    stop("SAEM ", label, " parameter names must be non-empty and unique.",
         call. = FALSE)
  }
  if (!setequal(parameter, expected) || length(parameter) != length(expected)) {
    stop(
      "SAEM ", label, " parameter rows do not match the exact HIV export ",
      "contract. Missing: ", paste(setdiff(expected, parameter), collapse = ","),
      "; unexpected: ", paste(setdiff(parameter, expected), collapse = ","),
      ".", call. = FALSE
    )
  }
  invisible(value)
}

.sab_saem_validate_numeric_column <- function(value, label) {
  if (!is.numeric(value) || any(!is.finite(value))) {
    stop("Every ", label, " must be finite numeric data.", call. = FALSE)
  }
  invisible(value)
}

.sab_saem_validate_individual_columns <- function(individuals, contract) {
  required <- c("ID", contract$expected_individual_latent_columns)
  missing <- setdiff(required, names(individuals))
  if (length(missing)) {
    stop(
      "SAEM individual artifact is missing latent column(s): ",
      paste(missing, collapse = ", "), ".", call. = FALSE
    )
  }
  for (name in contract$expected_individual_latent_columns) {
    .sab_saem_validate_numeric_column(
      individuals[[name]], paste0("individual `", name, "`")
    )
  }
  invisible(individuals)
}

.sab_saem_patient_ids <- function(value) {
  if (is.factor(value)) value <- as.character(value)
  if (is.numeric(value) || is.integer(value)) {
    if (any(!is.finite(value)) || any(value != floor(value))) {
      stop("Numeric SAEM patient IDs must be finite integers.", call. = FALSE)
    }
    value <- format(value, scientific = FALSE, trim = TRUE)
  }
  if (!is.character(value) || anyNA(value) || any(!nzchar(value)) ||
      any(grepl("^\\s|\\s$", value))) {
    stop("SAEM patient IDs must be non-empty canonical labels.",
         call. = FALSE)
  }
  value
}

.sab_saem_sha256 <- function(path) {
  executable <- Sys.which("sha256sum")
  if (!nzchar(executable)) {
    stop("sha256sum is required to fingerprint SAEM artifacts.",
         call. = FALSE)
  }
  output <- suppressWarnings(system2(
    executable, args = c("--", shQuote(normalizePath(path, mustWork = TRUE))),
    stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  if ((!is.null(status) && status != 0L) || length(output) != 1L) {
    stop("Could not calculate SAEM artifact SHA-256.", call. = FALSE)
  }
  hash <- sub("[[:space:]].*$", "", output[[1L]])
  if (!grepl("^[0-9a-f]{64}$", hash)) {
    stop("Malformed SHA-256 output for SAEM artifact.", call. = FALSE)
  }
  hash
}

.sab_saem_prediction_summary <- function(prediction, ordinary_reasons) {
  malformed <- !is.list(prediction) ||
    !is.logical(prediction$ok) || length(prediction$ok) != 1L ||
    is.na(prediction$ok)
  if (malformed) {
    return(list(ok = FALSE, reason = "malformed_prediction",
                ode_integrations = NA_integer_, contract_error = TRUE))
  }
  reason <- if (isTRUE(prediction$ok)) {
    "ok"
  } else if (is.character(prediction$reason) &&
             length(prediction$reason) == 1L &&
             !is.na(prediction$reason) && nzchar(prediction$reason)) {
    prediction$reason
  } else {
    "malformed_prediction_failure_reason"
  }
  ode_integrations <- .sab_saem_prediction_ode_count(prediction)
  known_reason <- if (isTRUE(prediction$ok)) {
    identical(reason, "ok")
  } else {
    reason %in% setdiff(ordinary_reasons, "ok")
  }
  list(
    ok = isTRUE(prediction$ok), reason = reason,
    ode_integrations = ode_integrations,
    contract_error = !known_reason || is.na(ode_integrations)
  )
}

.sab_saem_prediction_ode_count <- function(prediction) {
  if (!is.null(prediction$ode_integrations)) {
    value <- prediction$ode_integrations
    if ((!is.numeric(value) && !is.integer(value)) || length(value) != 1L ||
        !is.finite(value) || value != floor(value) || value < 0L) {
      return(NA_integer_)
    }
    return(as.integer(value))
  }
  if (isTRUE(prediction$ok)) {
    times <- prediction$adjusted_positive_times
    if (!is.numeric(times) || any(!is.finite(times))) return(NA_integer_)
    return(as.integer(length(times) > 0L))
  }
  if (is.character(prediction$reason) && length(prediction$reason) == 1L &&
      !is.na(prediction$reason) && prediction$reason %in%
      c("ode_failure", "invalid_ode_output", "invalid_ode_times")) {
    return(1L)
  }
  0L
}

.sab_saem_exact_row <- function(patient_id, population_log_density,
                                ode_integrations, prediction_ok,
                                prediction_reason, likelihood_calls,
                                loglik, callback_error) {
  finite_loglik <- is.numeric(loglik) && length(loglik) == 1L &&
    is.finite(loglik)
  data.frame(
    patient_id = patient_id,
    population_log_density = as.numeric(population_log_density),
    prediction_calls = 1L,
    ode_integrations = as.integer(ode_integrations),
    prediction_ok = isTRUE(prediction_ok),
    prediction_reason = as.character(prediction_reason),
    callback_error = isTRUE(callback_error),
    likelihood_calls = as.integer(likelihood_calls),
    loglik = if (is.numeric(loglik) && length(loglik) == 1L) {
      as.numeric(loglik)
    } else {
      NA_real_
    },
    finite_loglik = finite_loglik,
    ordinary_target_rejection = !isTRUE(callback_error) &&
      (!isTRUE(prediction_ok) || !finite_loglik),
    passed = !isTRUE(callback_error) && isTRUE(prediction_ok) &&
      finite_loglik,
    stringsAsFactors = FALSE
  )
}
