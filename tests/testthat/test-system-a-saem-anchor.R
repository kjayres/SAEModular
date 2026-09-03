if (!exists("sab_system_a_saem_anchor_contract", mode = "function")) {
  source(file.path(core_root, "R", "system_a_saem_anchor.R"), local = FALSE)
}

.sab_saem_test_data_path <- function() {
  normalizePath(
    file.path(
      core_root, "..", "saem_model", "data",
      "hiv_benchmark_canonical.csv"
    ),
    mustWork = TRUE
  )
}

.sab_saem_test_patient_ids <- function() {
  data <- utils::read.csv(
    .sab_saem_test_data_path(), stringsAsFactors = FALSE,
    check.names = FALSE
  )
  as.character(unique(data$ID))
}

.sab_saem_test_adapter <- function(failed_patient = NULL,
                                   callback_error_patient = NULL,
                                   nonfinite_density_patient = NULL) {
  contract <- sab_system_a_saem_anchor_contract()
  calls <- new.env(parent = emptyenv())
  calls$prediction <- 0L
  calls$likelihood <- 0L
  adapter <- list(
    target_fingerprint = contract$target_fingerprint,
    patient_ids = .sab_saem_test_patient_ids(),
    coordinate_names = contract$coordinate_names,
    eta_in_domain = function(eta) {
      is.numeric(eta) && identical(names(eta),
                                   contract$coordinate_names$population) &&
        all(is.finite(eta))
    },
    psi_in_domain = function(psi) {
      is.numeric(psi) && identical(names(psi),
                                   contract$coordinate_names$global) &&
        all(is.finite(psi))
    },
    log_population_density = function(patient_id, x, eta) {
      if (identical(patient_id, nonfinite_density_patient)) return(-Inf)
      -0.5 * sum(x^2) - sum(eta[grep("^log_omega_", names(eta))])
    },
    solve_prediction = function(patient_id, x, psi) {
      calls$prediction <- calls$prediction + 1L
      if (identical(patient_id, callback_error_patient)) {
        stop("deliberate prediction error")
      }
      if (identical(patient_id, failed_patient)) {
        return(list(
          ok = FALSE, reason = "ode_failure", ode_integrations = 1L
        ))
      }
      list(ok = TRUE, state = x, ode_integrations = 1L)
    },
    loglik_from_prediction = function(prediction, psi) {
      calls$likelihood <- calls$likelihood + 1L
      if (!isTRUE(prediction$ok)) return(-Inf)
      -sum(prediction$state^2) - sum(psi^2)
    },
    test_calls = calls
  )
  class(adapter) <- c("sab_saem_anchor_mock_adapter", "list")
  adapter
}

.sab_saem_test_fixture <- function() {
  contract <- sab_system_a_saem_anchor_contract()
  directory <- tempfile("sab_saem_anchor_")
  dir.create(directory)

  population_values <- stats::setNames(
    seq_along(contract$expected_population_parameters) / 10,
    contract$expected_population_parameters
  )
  population_values[c(
    "u_log_gamma", "u_log_mu_l", "u_log_sigma_v", "u_log_sigma_t"
  )] <- log(c(0.0021, 0.0092, 0.464, 0.254))
  population <- data.frame(
    parameter = rev(names(population_values)),
    estimate = unname(rev(population_values)),
    stringsAsFactors = FALSE
  )

  omega_sd <- stats::setNames(
    numeric(length(contract$expected_omega_parameters)),
    contract$expected_omega_parameters
  )
  omega_sd[unname(contract$omega_map)] <- seq(0.35, 1.05,
                                               length.out = 8L)
  omega <- data.frame(
    parameter = rev(names(omega_sd)),
    variance = unname(rev(omega_sd^2)),
    sd = unname(rev(omega_sd)),
    stringsAsFactors = FALSE
  )

  ids <- rev(as.integer(.sab_saem_test_patient_ids()))
  latent <- matrix(
    0, nrow = length(ids),
    ncol = length(contract$expected_individual_latent_columns),
    dimnames = list(NULL, contract$expected_individual_latent_columns)
  )
  for (index in seq_along(contract$coordinate_names$local)) {
    latent[, contract$coordinate_names$local[[index]]] <-
      ids / 1000 + index / 10
  }
  for (name in unname(contract$global_map)) {
    latent[, name] <- population_values[[name]]
  }
  individuals <- as.data.frame(latent, check.names = FALSE)
  individuals$ID <- ids
  individuals$unused_natural_scale_export <- exp(latent[, "u_log_lambda"])
  input_data <- utils::read.csv(
    .sab_saem_test_data_path(), stringsAsFactors = FALSE,
    check.names = FALSE
  )

  utils::write.csv(
    population, file.path(directory, contract$files[["population"]]),
    row.names = FALSE, quote = FALSE
  )
  utils::write.csv(
    omega, file.path(directory, contract$files[["omega"]]),
    row.names = FALSE, quote = FALSE
  )
  utils::write.csv(
    individuals, file.path(directory, contract$files[["individuals"]]),
    row.names = FALSE, quote = FALSE
  )
  copied <- file.copy(
    .sab_saem_test_data_path(),
    file.path(directory, contract$files[["input_data"]]),
    overwrite = FALSE, copy.mode = FALSE, copy.date = FALSE
  )
  if (!isTRUE(copied)) stop("Could not construct pinned input-data fixture.")
  run_tag <- basename(directory)
  jsonlite::write_json(
    list(
      project_root = dirname(dirname(directory)),
      workspace_root = dirname(dirname(dirname(directory))),
      run_tag = run_tag,
      run_dir = normalizePath(directory, mustWork = TRUE),
      hiv_data_path = .sab_saem_test_data_path(),
      fit_control = list(
        map = TRUE, fim = FALSE,
        `nbiter.saemix` = c(5L, 2L),
        `nbiter.sa` = 2L, `nbiter.burn` = 1L, `nbiter.map` = 1L,
        `nb.chains` = 1L, seed = 8L,
        `nbiter.mcmc` = c(2L, 2L, 2L, 0L)
      ),
      ode_control = list(
        method = "lsoda", rtol = 1e-8, atol = 1e-10,
        maxsteps = 100000L, time_eps = 1e-6
      )
    ),
    file.path(directory, contract$files[["run_config"]]),
    auto_unbox = TRUE, pretty = TRUE
  )
  audit_path <- file.path(directory, paste0(run_tag, ".audit.json"))
  jsonlite::write_json(
    list(
      status = "completed",
      started_at = "2026-01-01T00:00:00+0000",
      events = list(
        list(
          event = "hiv_benchmark_loaded", phase = "io",
          n_rows = nrow(input_data),
          n_patients = contract$expected_patient_count
        ),
        list(
          event = "saem_fit_completed", phase = "fit",
          ode_failure_count = 0L
        )
      ),
      updated_at = "2026-01-01T00:01:00+0000"
    ),
    audit_path, auto_unbox = TRUE, pretty = TRUE
  )
  list(
    directory = directory, contract = contract,
    population_values = population_values, omega_sd = omega_sd,
    population = population, omega = omega, individuals = individuals,
    input_data = input_data, audit_path = audit_path
  )
}

.sab_saem_test_rewrite <- function(fixture, artifact, value) {
  utils::write.csv(
    value,
    file.path(fixture$directory, fixture$contract$files[[artifact]]),
    row.names = FALSE, quote = FALSE
  )
}

testthat::test_that("the SAEM mapping contract is complete and explicit", {
  contract <- sab_system_a_saem_anchor_contract()
  testthat::expect_identical(contract$expected_patient_count, 115L)
  testthat::expect_length(contract$coordinate_names$local, 8L)
  testthat::expect_length(contract$coordinate_names$population, 17L)
  testthat::expect_length(contract$coordinate_names$global, 4L)
  testthat::expect_length(contract$population_location_map, 9L)
  testthat::expect_length(contract$omega_map, 8L)
  testthat::expect_identical(
    contract$population_location_map,
    c(
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
  )
  testthat::expect_identical(
    contract$global_map,
    c(
      log_gamma_pop = "u_log_gamma",
      log_mu_l_pop = "u_log_mu_l",
      log_sigma_v = "u_log_sigma_v",
      log_sigma_t = "u_log_sigma_t"
    )
  )
  testthat::expect_identical(
    contract$omega_map,
    c(
      log_omega_lambda = "omega_u_log_lambda",
      log_omega_mu_t = "omega_u_log_mu_t",
      log_omega_mu_a = "omega_u_log_mu_a",
      log_omega_p = "omega_u_log_p",
      log_omega_alpha_l = "omega_u_log_alpha_l",
      log_omega_pi = "omega_u_pi",
      log_omega_eta_rti = "omega_u_eta_rti",
      log_omega_eta_pi = "omega_u_eta_pi"
    )
  )
  testthat::expect_identical(
    anyDuplicated(contract$expected_population_parameters), 0L
  )
  testthat::expect_identical(
    anyDuplicated(contract$expected_omega_parameters), 0L
  )
})

testthat::test_that("completed SAEM artifacts map into canonical coordinates", {
  fixture <- .sab_saem_test_fixture()
  on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
  adapter <- .sab_saem_test_adapter()

  anchor <- sab_read_system_a_saem_anchor(fixture$directory, adapter)
  testthat::expect_s3_class(anchor, "sab_system_a_saem_anchor")
  testthat::expect_invisible(
    sab_validate_system_a_saem_anchor(anchor, adapter)
  )
  testthat::expect_identical(anchor$patient_ids, adapter$patient_ids)
  testthat::expect_identical(
    dim(anchor$patient_states), c(115L, 8L)
  )
  testthat::expect_identical(
    rownames(anchor$patient_states), adapter$patient_ids
  )
  testthat::expect_identical(
    colnames(anchor$patient_states), fixture$contract$coordinate_names$local
  )
  expected_eta <- c(
    mu_log_lambda = fixture$population_values[["u_log_lambda"]],
    mu_log_mu_t = fixture$population_values[["u_log_mu_t"]],
    mu_log_mu_a = fixture$population_values[["u_log_mu_a"]],
    mu_log_p = fixture$population_values[["u_log_p"]],
    mu_log_alpha_l = fixture$population_values[["u_log_alpha_l"]],
    mu_u_pi = fixture$population_values[["u_pi"]],
    mu_u_eta_rti = fixture$population_values[["u_eta_rti"]],
    mu_u_eta_pi = fixture$population_values[["u_eta_pi"]],
    beta_nelf =
      fixture$population_values[["beta_treat_nelf(u_eta_pi)"]],
    log_omega_lambda = log(fixture$omega_sd[["omega_u_log_lambda"]]),
    log_omega_mu_t = log(fixture$omega_sd[["omega_u_log_mu_t"]]),
    log_omega_mu_a = log(fixture$omega_sd[["omega_u_log_mu_a"]]),
    log_omega_p = log(fixture$omega_sd[["omega_u_log_p"]]),
    log_omega_alpha_l =
      log(fixture$omega_sd[["omega_u_log_alpha_l"]]),
    log_omega_pi = log(fixture$omega_sd[["omega_u_pi"]]),
    log_omega_eta_rti = log(fixture$omega_sd[["omega_u_eta_rti"]]),
    log_omega_eta_pi = log(fixture$omega_sd[["omega_u_eta_pi"]])
  )
  expected_psi <- c(
    log_gamma_pop = fixture$population_values[["u_log_gamma"]],
    log_mu_l_pop = fixture$population_values[["u_log_mu_l"]],
    log_sigma_v = fixture$population_values[["u_log_sigma_v"]],
    log_sigma_t = fixture$population_values[["u_log_sigma_t"]]
  )
  testthat::expect_equal(anchor$eta, expected_eta, tolerance = 1e-14)
  testthat::expect_equal(anchor$psi, expected_psi, tolerance = 1e-14)
  testthat::expect_equal(
    anchor$patient_states["2", ],
    stats::setNames(seq_len(8L) / 10 + 0.002,
                    fixture$contract$coordinate_names$local),
    tolerance = 1e-14
  )
  testthat::expect_true(anchor$source$reordered_to_canonical_ids)
  testthat::expect_identical(
    names(anchor$source$artifact_sha256),
    c(names(fixture$contract$files), "audit")
  )
  testthat::expect_true(all(grepl(
    "^[0-9a-f]{64}$", anchor$source$artifact_sha256
  )))
  testthat::expect_identical(anchor$source$audit$status, "completed")
  testthat::expect_identical(
    anchor$source$run_config$fit_control$k1, 5L
  )
  testthat::expect_identical(
    anchor$source$run_config$fit_control$k2, 2L
  )
  testthat::expect_identical(
    anchor$source$run_config$fit_control$nbiter_mcmc,
    c(2L, 2L, 2L, 0L)
  )
  testthat::expect_true(anchor$source$run_config$fit_control$map)
  testthat::expect_identical(
    anchor$source$run_config$hiv_data_sha256,
    fixture$contract$canonical_hiv_data_sha256
  )
})

testthat::test_that("completed-run audit and configuration are mandatory", {
  adapter <- .sab_saem_test_adapter()

  fixture <- .sab_saem_test_fixture()
  on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
  audit <- jsonlite::fromJSON(fixture$audit_path, simplifyVector = FALSE)
  audit$status <- "running"
  jsonlite::write_json(
    audit, fixture$audit_path, auto_unbox = TRUE, pretty = TRUE
  )
  testthat::expect_error(
    sab_read_system_a_saem_anchor(fixture$directory, adapter),
    "must be exactly `completed`", fixed = TRUE
  )

  fixture2 <- .sab_saem_test_fixture()
  on.exit(unlink(fixture2$directory, recursive = TRUE), add = TRUE)
  file.copy(
    fixture2$audit_path,
    file.path(fixture2$directory, "second.audit.json")
  )
  testthat::expect_error(
    sab_read_system_a_saem_anchor(fixture2$directory, adapter),
    "exactly one", fixed = TRUE
  )

  fixture3 <- .sab_saem_test_fixture()
  on.exit(unlink(fixture3$directory, recursive = TRUE), add = TRUE)
  tampered_data <- fixture3$input_data
  tampered_data$DV[[1L]] <- tampered_data$DV[[1L]] + 1
  .sab_saem_test_rewrite(fixture3, "input_data", tampered_data)
  testthat::expect_error(
    sab_read_system_a_saem_anchor(fixture3$directory, adapter),
    "pinned canonical HIV data", fixed = TRUE
  )
})

testthat::test_that("SAEM conversion fails closed on bad parameter artifacts", {
  adapter <- .sab_saem_test_adapter()

  fixture <- .sab_saem_test_fixture()
  on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
  duplicate <- fixture$population
  duplicate$parameter[[2L]] <- duplicate$parameter[[1L]]
  .sab_saem_test_rewrite(fixture, "population", duplicate)
  testthat::expect_error(
    sab_read_system_a_saem_anchor(fixture$directory, adapter),
    "non-empty and unique", fixed = TRUE
  )

  fixture2 <- .sab_saem_test_fixture()
  on.exit(unlink(fixture2$directory, recursive = TRUE), add = TRUE)
  missing <- fixture2$population[-1L, ]
  .sab_saem_test_rewrite(fixture2, "population", missing)
  testthat::expect_error(
    sab_read_system_a_saem_anchor(fixture2$directory, adapter),
    "do not match the exact HIV export contract", fixed = TRUE
  )

  fixture3 <- .sab_saem_test_fixture()
  on.exit(unlink(fixture3$directory, recursive = TRUE), add = TRUE)
  inconsistent <- fixture3$omega
  inconsistent$variance[[1L]] <- inconsistent$variance[[1L]] + 0.1
  .sab_saem_test_rewrite(fixture3, "omega", inconsistent)
  testthat::expect_error(
    sab_read_system_a_saem_anchor(fixture3$directory, adapter),
    "variance and SD columns are inconsistent", fixed = TRUE
  )
})

testthat::test_that("mapped zero omega and nonzero fixed omega are rejected", {
  adapter <- .sab_saem_test_adapter()
  fixture <- .sab_saem_test_fixture()
  on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
  row <- match("omega_u_log_lambda", fixture$omega$parameter)
  fixture$omega$sd[[row]] <- 0
  fixture$omega$variance[[row]] <- 0
  .sab_saem_test_rewrite(fixture, "omega", fixture$omega)
  testthat::expect_error(
    sab_read_system_a_saem_anchor(fixture$directory, adapter),
    "strictly positive", fixed = TRUE
  )

  fixture2 <- .sab_saem_test_fixture()
  on.exit(unlink(fixture2$directory, recursive = TRUE), add = TRUE)
  row <- match("omega_u_log_gamma", fixture2$omega$parameter)
  fixture2$omega$sd[[row]] <- 0.01
  fixture2$omega$variance[[row]] <- 0.0001
  .sab_saem_test_rewrite(fixture2, "omega", fixture2$omega)
  testthat::expect_error(
    sab_read_system_a_saem_anchor(fixture2$directory, adapter),
    "must be exactly zero", fixed = TRUE
  )

  fixture3 <- .sab_saem_test_fixture()
  on.exit(unlink(fixture3$directory, recursive = TRUE), add = TRUE)
  row <- match("omega_u_log_lambda", fixture3$omega$parameter)
  fixture3$omega$sd[[row]] <- 1e-6
  fixture3$omega$variance[[row]] <- 0
  .sab_saem_test_rewrite(fixture3, "omega", fixture3$omega)
  testthat::expect_error(
    sab_read_system_a_saem_anchor(fixture3$directory, adapter),
    "variance and SD columns are inconsistent", fixed = TRUE
  )

  fixture4 <- .sab_saem_test_fixture()
  on.exit(unlink(fixture4$directory, recursive = TRUE), add = TRUE)
  row <- match("omega_u_log_gamma", fixture4$omega$parameter)
  fixture4$omega$variance[[row]] <- 1e-12
  .sab_saem_test_rewrite(fixture4, "omega", fixture4$omega)
  testthat::expect_error(
    sab_read_system_a_saem_anchor(fixture4$directory, adapter),
    "variance and SD columns are inconsistent", fixed = TRUE
  )
})

testthat::test_that("patient IDs, states, and fixed columns are fail-closed", {
  adapter <- .sab_saem_test_adapter()

  fixture <- .sab_saem_test_fixture()
  on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
  fixture$individuals$ID[[1L]] <- fixture$individuals$ID[[2L]]
  .sab_saem_test_rewrite(fixture, "individuals", fixture$individuals)
  testthat::expect_error(
    sab_read_system_a_saem_anchor(fixture$directory, adapter),
    "duplicate patient IDs", fixed = TRUE
  )

  fixture2 <- .sab_saem_test_fixture()
  on.exit(unlink(fixture2$directory, recursive = TRUE), add = TRUE)
  fixture2$individuals$ID[[1L]] <- 999L
  .sab_saem_test_rewrite(fixture2, "individuals", fixture2$individuals)
  testthat::expect_error(
    sab_read_system_a_saem_anchor(fixture2$directory, adapter),
    "do not match", fixed = TRUE
  )

  fixture3 <- .sab_saem_test_fixture()
  on.exit(unlink(fixture3$directory, recursive = TRUE), add = TRUE)
  fixture3$individuals$u_log_lambda[[1L]] <- NA_real_
  .sab_saem_test_rewrite(fixture3, "individuals", fixture3$individuals)
  testthat::expect_error(
    sab_read_system_a_saem_anchor(fixture3$directory, adapter),
    "must be finite numeric data", fixed = TRUE
  )

  fixture4 <- .sab_saem_test_fixture()
  on.exit(unlink(fixture4$directory, recursive = TRUE), add = TRUE)
  fixture4$individuals$u_log_gamma[[1L]] <-
    fixture4$individuals$u_log_gamma[[1L]] + 0.01
  .sab_saem_test_rewrite(fixture4, "individuals", fixture4$individuals)
  testthat::expect_error(
    sab_read_system_a_saem_anchor(fixture4$directory, adapter),
    "disagree with the SAEM population artifact", fixed = TRUE
  )
})

testthat::test_that("target domains and population densities are validated", {
  fixture <- .sab_saem_test_fixture()
  on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)

  bad_domain <- .sab_saem_test_adapter()
  bad_domain$eta_in_domain <- function(eta) FALSE
  testthat::expect_error(
    sab_read_system_a_saem_anchor(fixture$directory, bad_domain),
    "eta is outside", fixed = TRUE
  )

  bad_density <- .sab_saem_test_adapter(nonfinite_density_patient = "19")
  testthat::expect_error(
    sab_read_system_a_saem_anchor(fixture$directory, bad_density),
    "population density", fixed = TRUE
  )
})

testthat::test_that("sealed exact validation records a complete ODE ledger", {
  fixture <- .sab_saem_test_fixture()
  on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
  adapter <- .sab_saem_test_adapter()

  anchor <- sab_read_system_a_saem_anchor(
    fixture$directory, adapter, validate_exact = TRUE
  )
  exact <- anchor$validation$exact_target
  testthat::expect_s3_class(
    exact, "sab_system_a_saem_exact_validation"
  )
  testthat::expect_true(exact$passed)
  testthat::expect_true(exact$completed_without_callback_error)
  testthat::expect_true(all(exact$ledger$passed))
  testthat::expect_false(any(exact$ledger$callback_error))
  testthat::expect_identical(exact$summary$patients, 115L)
  testthat::expect_identical(exact$summary$prediction_calls, 115L)
  testthat::expect_identical(exact$summary$ode_integrations, 115L)
  testthat::expect_identical(exact$summary$likelihood_calls, 115L)
  testthat::expect_identical(adapter$test_calls$prediction, 115L)
  testthat::expect_identical(adapter$test_calls$likelihood, 115L)
})

testthat::test_that("ordinary and callback exact-target failures stay visible", {
  fixture <- .sab_saem_test_fixture()
  on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
  adapter <- .sab_saem_test_adapter(
    failed_patient = "7", callback_error_patient = "9"
  )
  anchor <- sab_read_system_a_saem_anchor(fixture$directory, adapter)
  exact <- sab_validate_system_a_saem_anchor_exact(anchor, adapter)

  testthat::expect_false(exact$passed)
  testthat::expect_false(exact$completed_without_callback_error)
  testthat::expect_false(exact$ledger$passed[exact$ledger$patient_id == "7"])
  testthat::expect_true(
    exact$ledger$ordinary_target_rejection[
      exact$ledger$patient_id == "7"
    ]
  )
  testthat::expect_false(
    exact$ledger$callback_error[exact$ledger$patient_id == "7"]
  )
  testthat::expect_identical(
    exact$ledger$prediction_reason[exact$ledger$patient_id == "7"],
    "ode_failure"
  )
  testthat::expect_match(
    exact$ledger$prediction_reason[exact$ledger$patient_id == "9"],
    "prediction_callback_error"
  )
  testthat::expect_true(
    exact$ledger$callback_error[exact$ledger$patient_id == "9"]
  )
  testthat::expect_identical(exact$summary$unknown_ode_counts, 1L)
  testthat::expect_identical(exact$summary$nonfinite_loglik, 2L)
  testthat::expect_identical(exact$summary$likelihood_calls, 114L)
  testthat::expect_identical(exact$summary$ordinary_target_rejections, 1L)
  testthat::expect_identical(exact$summary$callback_errors, 1L)
  testthat::expect_false(exact$summary$completed_without_callback_error)
})

testthat::test_that("ordinary target rejection does not imply callback failure", {
  fixture <- .sab_saem_test_fixture()
  on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
  adapter <- .sab_saem_test_adapter(failed_patient = "7")
  anchor <- sab_read_system_a_saem_anchor(fixture$directory, adapter)
  exact <- sab_validate_system_a_saem_anchor_exact(anchor, adapter)

  testthat::expect_false(exact$passed)
  testthat::expect_true(exact$completed_without_callback_error)
  testthat::expect_identical(exact$summary$ordinary_target_rejections, 1L)
  testthat::expect_identical(exact$summary$callback_errors, 0L)
})

testthat::test_that("unknown or contradictory prediction reasons are errors", {
  fixture <- .sab_saem_test_fixture()
  on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)

  adapter <- .sab_saem_test_adapter()
  original_solve <- adapter$solve_prediction
  adapter$solve_prediction <- function(patient_id, x, psi) {
    if (identical(patient_id, "7")) {
      return(list(ok = FALSE, reason = "programming_bug",
                  ode_integrations = 0L))
    }
    original_solve(patient_id, x, psi)
  }
  anchor <- sab_read_system_a_saem_anchor(fixture$directory, adapter)
  exact <- sab_validate_system_a_saem_anchor_exact(anchor, adapter)
  testthat::expect_false(exact$completed_without_callback_error)
  testthat::expect_true(
    exact$ledger$callback_error[exact$ledger$patient_id == "7"]
  )

  adapter2 <- .sab_saem_test_adapter()
  original_solve2 <- adapter2$solve_prediction
  adapter2$solve_prediction <- function(patient_id, x, psi) {
    if (identical(patient_id, "7")) {
      return(list(ok = FALSE, reason = "ok", ode_integrations = 0L))
    }
    original_solve2(patient_id, x, psi)
  }
  anchor2 <- sab_read_system_a_saem_anchor(fixture$directory, adapter2)
  exact2 <- sab_validate_system_a_saem_anchor_exact(anchor2, adapter2)
  testthat::expect_false(exact2$completed_without_callback_error)
  testthat::expect_true(
    exact2$ledger$callback_error[exact2$ledger$patient_id == "7"]
  )
})

testthat::test_that("anchor validation detects stale population-density data", {
  fixture <- .sab_saem_test_fixture()
  on.exit(unlink(fixture$directory, recursive = TRUE), add = TRUE)
  adapter <- .sab_saem_test_adapter()
  anchor <- sab_read_system_a_saem_anchor(fixture$directory, adapter)
  anchor$patient_states["2", "u_log_lambda"] <-
    anchor$patient_states["2", "u_log_lambda"] + 0.5

  testthat::expect_error(
    sab_validate_system_a_saem_anchor(anchor, adapter),
    "stale or malformed", fixed = TRUE
  )
})

testthat::test_that("completed smoke artifacts bind SAEM internal rows to IDs", {
  testthat::skip_if_not(
    identical(tolower(Sys.getenv("SAB_SYSTEM_A_INTEGRATION")), "true"),
    "Set SAB_SYSTEM_A_INTEGRATION=true in an HPC test job."
  )
  workspace_root <- normalizePath(
    file.path(core_root, "..", ".."), mustWork = TRUE
  )
  run_directory <- file.path(
    workspace_root, "projects", "saem_model", "outputs",
    "saemodular_system_a_anchor_121939"
  )
  testthat::skip_if_not(dir.exists(run_directory))
  adapter <- sab_load_system_a_adapter(workspace_root)
  anchor <- sab_read_system_a_saem_anchor(run_directory, adapter)

  testthat::expect_true(anchor$source$row_order_validation$validated)
  testthat::expect_identical(
    anchor$source$row_order_validation$method,
    "saemix_internal_index_binding"
  )
  testthat::expect_identical(
    anchor$source$row_order_validation$patient_ids,
    anchor$source$input_patient_ids
  )
  testthat::expect_true(file.exists(
    anchor$source$row_order_validation$fit_path
  ))
  testthat::expect_match(
    anchor$source$row_order_validation$fit_sha256, "^[0-9a-f]{64}$"
  )
  # This run proves the binding machinery only.  Its K1=5, K2=2 schedule is
  # intentionally below the separate production-anchor gate.
  testthat::expect_identical(anchor$source$run_config$fit_control$k1, 5L)
  testthat::expect_identical(anchor$source$run_config$fit_control$k2, 2L)
})
