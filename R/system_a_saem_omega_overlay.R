# Fail-closed correction for the System A SAEM initial covariance.
#
# saemix initially simulates all 12 individual coordinates, including the four
# whose IIV is structurally disabled by covariance.model.  A singular
# omega.init therefore makes initialiseMainAlgo replace the complete matrix by
# the identity before iteration 0.  This overlay keeps the eight published IIV
# starts and assigns variance 1e-6 (SD 0.001 on the latent scale) to the four
# structurally fixed coordinates.  This remains numerically positive while its
# transient initialization perturbation is negligible; it is never estimated
# and disappears when saemix removes non-IIV coordinates after initialization.

.sab_saem_model_coordinates <- c(
  "u_log_lambda", "u_log_gamma", "u_log_mu_t", "u_log_mu_l",
  "u_log_mu_a", "u_log_p", "u_log_alpha_l", "u_pi",
  "u_eta_rti", "u_eta_pi", "u_log_sigma_v", "u_log_sigma_t"
)

.sab_saem_active_iiv_coordinates <- c(
  "u_log_lambda", "u_log_mu_t", "u_log_mu_a", "u_log_p",
  "u_log_alpha_l", "u_pi", "u_eta_rti", "u_eta_pi"
)

.sab_saem_fixed_iiv_coordinates <- setdiff(
  .sab_saem_model_coordinates, .sab_saem_active_iiv_coordinates
)

.sab_saem_initial_sd_fields <- c(
  u_log_lambda = "omega_lambda",
  u_log_mu_t = "omega_mu_t",
  u_log_mu_a = "omega_mu_a",
  u_log_p = "omega_p",
  u_log_alpha_l = "omega_alpha_l",
  u_pi = "omega_pi",
  u_eta_rti = "omega_eta_rti",
  u_eta_pi = "omega_eta_pi"
)

sab_corrected_saem_omega_contract <- function(initials,
                                               dummy_fixed_variance = 1e-6) {
  if (!is.list(initials) || !is.list(initials$omega)) {
    stop("SAEM initials must contain an omega list.", call. = FALSE)
  }
  if (length(dummy_fixed_variance) != 1L ||
      !is.finite(dummy_fixed_variance) || dummy_fixed_variance <= 0) {
    stop("dummy_fixed_variance must be one finite positive number.",
         call. = FALSE)
  }

  missing_fields <- setdiff(
    unname(.sab_saem_initial_sd_fields), names(initials$omega)
  )
  if (length(missing_fields)) {
    stop(
      "Published SAEM omega fields are missing: ",
      paste(missing_fields, collapse = ", "), ".",
      call. = FALSE
    )
  }
  published_sd <- vapply(
    unname(.sab_saem_initial_sd_fields),
    function(field) as.numeric(initials$omega[[field]]),
    numeric(1)
  )
  names(published_sd) <- names(.sab_saem_initial_sd_fields)
  if (any(!is.finite(published_sd)) || any(published_sd <= 0)) {
    stop("All eight published IIV standard deviations must be positive.",
         call. = FALSE)
  }

  variance <- stats::setNames(
    rep(as.numeric(dummy_fixed_variance),
        length(.sab_saem_model_coordinates)),
    .sab_saem_model_coordinates
  )
  variance[names(published_sd)] <- published_sd^2
  omega_init <- diag(unname(variance), nrow = length(variance))
  dimnames(omega_init) <- list(names(variance), names(variance))

  estimated <- .sab_saem_model_coordinates %in%
    .sab_saem_active_iiv_coordinates
  covariance_model <- diag(as.integer(estimated), nrow = length(estimated))
  dimnames(covariance_model) <- dimnames(omega_init)

  list(
    schema_version = "sab_system_a_corrected_omega_v1",
    dummy_fixed_variance = as.numeric(dummy_fixed_variance),
    coordinates = .sab_saem_model_coordinates,
    active_coordinates = .sab_saem_active_iiv_coordinates,
    fixed_coordinates = .sab_saem_fixed_iiv_coordinates,
    published_sd = published_sd,
    omega_init = omega_init,
    covariance_model = covariance_model
  )
}

.sab_assert_numeric_matrix <- function(x, expected_dim, label) {
  if (!is.matrix(x) || !is.numeric(x) ||
      !identical(dim(x), as.integer(expected_dim)) || any(!is.finite(x))) {
    stop(label, " has the wrong dimensions or non-finite entries.",
         call. = FALSE)
  }
  invisible(TRUE)
}

sab_assert_corrected_saem_matrices <- function(model_names,
                                                covariance_model,
                                                omega_init,
                                                contract,
                                                tolerance = 1e-14) {
  if (!identical(as.character(model_names), contract$coordinates)) {
    stop("SAEM model coordinate order differs from the correction contract.",
         call. = FALSE)
  }
  expected_dim <- c(length(contract$coordinates),
                    length(contract$coordinates))
  .sab_assert_numeric_matrix(covariance_model, expected_dim,
                             "covariance.model")
  .sab_assert_numeric_matrix(omega_init, expected_dim, "omega.init")
  if (!isSymmetric(omega_init, tol = tolerance)) {
    stop("Corrected omega.init is not symmetric.", call. = FALSE)
  }
  if (inherits(try(chol(omega_init), silent = TRUE), "try-error")) {
    stop("Corrected omega.init is not positive definite.", call. = FALSE)
  }
  reciprocal_condition <- rcond(omega_init)
  if (!is.finite(reciprocal_condition) || reciprocal_condition < 1e-9) {
    stop("Corrected omega.init is too ill-conditioned for initialization.",
         call. = FALSE)
  }
  if (!isTRUE(all.equal(
    unname(covariance_model), unname(contract$covariance_model),
    tolerance = tolerance, check.attributes = FALSE
  ))) {
    stop("covariance.model no longer fixes exactly the four contracted IIVs.",
         call. = FALSE)
  }
  if (!isTRUE(all.equal(
    unname(omega_init), unname(contract$omega_init),
    tolerance = tolerance, check.attributes = FALSE
  ))) {
    stop("Constructed omega.init differs from the corrected start contract.",
         call. = FALSE)
  }
  invisible(TRUE)
}

sab_assert_fitted_structural_omega <- function(fitted_omega, contract,
                                                tolerance = 1e-14) {
  expected_dim <- c(length(contract$coordinates),
                    length(contract$coordinates))
  .sab_assert_numeric_matrix(fitted_omega, expected_dim,
                             "fitted omega")
  fitted_variance <- diag(fitted_omega)
  names(fitted_variance) <- contract$coordinates
  fixed_variance <- fitted_variance[contract$fixed_coordinates]
  if (any(abs(fixed_variance) > tolerance)) {
    stop(
      "The fitted omega does not retain zero variance for all four ",
      "structurally fixed IIV coordinates.", call. = FALSE
    )
  }
  active_variance <- fitted_variance[contract$active_coordinates]
  if (any(active_variance <= 0)) {
    stop("At least one fitted active IIV variance is not positive.",
         call. = FALSE)
  }
  invisible(TRUE)
}

.sab_corrected_saem_run_directory <- function() {
  output_root <- Sys.getenv("SAEM_OUTPUT_ROOT", unset = "")
  run_tag <- Sys.getenv("SAEM_RUN_TAG", unset = "")
  if (!nzchar(output_root) || !nzchar(run_tag)) return(NA_character_)
  file.path(output_root, run_tag)
}

.sab_atomic_write_csv <- function(x, path) {
  directory <- dirname(path)
  if (!dir.exists(directory)) {
    stop("Audit output directory does not exist: ", directory,
         call. = FALSE)
  }
  temporary <- tempfile(pattern = paste0(basename(path), "."),
                        tmpdir = directory)
  on.exit(unlink(temporary), add = TRUE)
  utils::write.csv(x, temporary, row.names = FALSE, quote = TRUE)
  if (!file.rename(temporary, path)) {
    stop("Could not atomically publish ", path, ".", call. = FALSE)
  }
  invisible(path)
}

sab_write_constructed_omega_audit <- function(contract, run_directory) {
  if (!is.character(run_directory) || length(run_directory) != 1L ||
      is.na(run_directory) || !nzchar(run_directory)) {
    stop("A concrete SAEM run directory is required for the omega audit.",
         call. = FALSE)
  }
  variance <- diag(contract$omega_init)
  table <- data.frame(
    coordinate = contract$coordinates,
    covariance_model_diagonal = diag(contract$covariance_model),
    omega_initial_variance = unname(variance),
    initialization_role = ifelse(
      contract$coordinates %in% contract$active_coordinates,
      "published_iiv_variance", "positive_dummy_for_structurally_fixed_iiv"
    ),
    stringsAsFactors = FALSE
  )
  .sab_atomic_write_csv(
    table, file.path(run_directory, "constructed_omega_initialization.csv")
  )
}

sab_install_corrected_saem_model_builder <- function(
    environment = parent.frame(), dummy_fixed_variance = 1e-6) {
  if (!exists("build_hiv_saemix_model_object", envir = environment,
              inherits = FALSE, mode = "function")) {
    stop("Upstream build_hiv_saemix_model_object must be loaded first.",
         call. = FALSE)
  }
  if (exists(".sab_upstream_hiv_saemix_model_builder", envir = environment,
             inherits = FALSE)) {
    stop("The corrected SAEM model overlay was installed more than once.",
         call. = FALSE)
  }

  upstream_builder <- get(
    "build_hiv_saemix_model_object", envir = environment, inherits = FALSE
  )
  assign(".sab_upstream_hiv_saemix_model_builder", upstream_builder,
         envir = environment)

  corrected_builder <- function(initials, ode_control = list()) {
    model <- .sab_upstream_hiv_saemix_model_builder(
      initials, ode_control = ode_control
    )
    contract <- sab_corrected_saem_omega_contract(
      initials, dummy_fixed_variance = dummy_fixed_variance
    )

    # Use the package setter so S4 validity is checked after replacement.
    model["omega.init"] <- contract$omega_init
    sab_assert_corrected_saem_matrices(
      model["name.modpar"], model["covariance.model"],
      model["omega.init"], contract
    )
    run_directory <- .sab_corrected_saem_run_directory()
    if (is.na(run_directory)) {
      stop("SAEM_OUTPUT_ROOT and SAEM_RUN_TAG must identify the audit output.",
           call. = FALSE)
    }
    sab_write_constructed_omega_audit(contract, run_directory)
    model
  }
  assign("build_hiv_saemix_model_object", corrected_builder,
         envir = environment)
  invisible(TRUE)
}

sab_assert_fit_iteration0 <- function(fit, initials, run_directory,
                                       dummy_fixed_variance = 1e-6,
                                       tolerance = 1e-12) {
  contract <- sab_corrected_saem_omega_contract(
    initials, dummy_fixed_variance = dummy_fixed_variance
  )
  if (!methods::is(fit, "SaemixObject")) {
    stop("The fitted object is not a SaemixObject.", call. = FALSE)
  }
  fitted_model <- methods::slot(fit, "model")
  sab_assert_corrected_saem_matrices(
    fitted_model["name.modpar"], fitted_model["covariance.model"],
    fitted_model["omega.init"], contract, tolerance = tolerance
  )

  results <- methods::slot(fit, "results")
  sab_assert_fitted_structural_omega(
    as.matrix(methods::slot(results, "omega")), contract,
    tolerance = tolerance
  )
  allpar <- methods::slot(results, "allpar")
  random_names <- as.character(methods::slot(results, "name.random"))
  if (!is.matrix(allpar) || !is.numeric(allpar) || nrow(allpar) < 1L ||
      any(!is.finite(allpar[1L, ]))) {
    stop("SAEM results do not contain a finite iteration-0 allpar row.",
         call. = FALSE)
  }
  if (length(random_names) != length(contract$active_coordinates) ||
      ncol(allpar) < length(random_names)) {
    stop("SAEM allpar has an unexpected random-effect block.",
         call. = FALSE)
  }
  random_columns <- tail(seq_len(ncol(allpar)), length(random_names))
  if (!identical(colnames(allpar)[random_columns], random_names)) {
    stop("SAEM allpar random-effect columns are not the final named block.",
         call. = FALSE)
  }

  expected <- diag(contract$omega_init)[contract$active_coordinates]
  observed <- as.numeric(allpar[1L, random_columns])
  if (!isTRUE(all.equal(
    observed, unname(expected), tolerance = tolerance,
    check.attributes = FALSE
  ))) {
    stop(
      "Iteration-0 allpar does not contain the eight published IIV variances; ",
      "saemix may have replaced omega.init.", call. = FALSE
    )
  }

  table <- data.frame(
    allpar_column = colnames(allpar),
    iteration_0_value = as.numeric(allpar[1L, ]),
    block = ifelse(seq_len(ncol(allpar)) %in% random_columns,
                   "iiv_variance", "fixed_effect"),
    expected_published_variance = NA_real_,
    stringsAsFactors = FALSE
  )
  table$expected_published_variance[random_columns] <- unname(expected)
  .sab_atomic_write_csv(
    table, file.path(run_directory, "iteration_0_allpar.csv")
  )
  invisible(TRUE)
}

sab_build_corrected_saem_runner <- function(upstream_lines) {
  if (!is.character(upstream_lines) || !length(upstream_lines)) {
    stop("upstream_lines must contain the pinned runner source.",
         call. = FALSE)
  }
  replace_once <- function(lines, pattern, replacement, label, fixed = TRUE) {
    indices <- grep(pattern, lines, fixed = fixed)
    if (length(indices) != 1L) {
      stop("Expected exactly one ", label, " insertion point; found ",
           length(indices), ".", call. = FALSE)
    }
    index <- indices[[1L]]
    before <- if (index > 1L) lines[seq_len(index - 1L)] else character()
    after <- if (index < length(lines)) {
      lines[seq.int(index + 1L, length(lines))]
    } else {
      character()
    }
    c(before, replacement(lines[index]), after)
  }

  lines <- replace_once(
    upstream_lines,
    'source(file.path(project_root, "R", "hiv_latent_saem_model.R"))',
    function(line) c(
      line,
      'source(file.path(Sys.getenv("SAEMODULAR_PROJECT_ROOT"), "R",',
      '                 "system_a_saem_omega_overlay.R"))',
      'sab_install_corrected_saem_model_builder(environment = environment())'
    ),
    "model-overlay"
  )
  lines <- replace_once(
    lines,
    'run_dir <- file.path(project_root, "outputs", run_tag)',
    function(line) c(
      'saem_output_root <- Sys.getenv("SAEM_OUTPUT_ROOT",',
      '                               unset = file.path(project_root, "outputs"))',
      'run_dir <- file.path(saem_output_root, run_tag)'
    ),
    "output-root"
  )
  lines <- replace_once(
    lines,
    "    print = TRUE,",
    function(line) c(line, "    warnings = TRUE,"),
    "warnings-control"
  )
  lines <- replace_once(
    lines,
    "      fit <- saemix::saemix(saem_model, saem_data, control = fit_control)",
    function(line) c(
      line,
      "      sab_assert_fit_iteration0(fit, initials, run_dir)"
    ),
    "iteration-0 assertion"
  )
  lines
}
