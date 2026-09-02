#!/usr/bin/env Rscript

# ORACLE System A mechanism pilot.
#
# This is deliberately not a deployable SAEM workflow and not a retained
# posterior chain.  The target-matched full-joint comparator is used only to
# choose a fixed global anchor, patient starting values, and the dynamic-psi
# scale.  At that anchor, independent exact patient-conditional pCN chains fit
# frozen affine maps.  A separate conditional bank then compares fixed-x and
# affine endpoints using the certified VODE-BDF likelihood.

options(warn = 1)

sab_stop <- function(...) stop(..., call. = FALSE)

sab_integer_env <- function(name, default, minimum = 1L) {
  text <- Sys.getenv(name, unset = "")
  if (!nzchar(text)) return(as.integer(default))
  value <- suppressWarnings(as.numeric(text))
  if (length(value) != 1L || !is.finite(value) || value != floor(value) ||
      value < minimum) {
    sab_stop(name, " must be an integer >= ", minimum, ".")
  }
  as.integer(value)
}

sab_numeric_env <- function(name, default, lower = -Inf, upper = Inf) {
  text <- Sys.getenv(name, unset = "")
  if (!nzchar(text)) return(as.numeric(default))
  value <- suppressWarnings(as.numeric(text))
  if (length(value) != 1L || !is.finite(value) ||
      value < lower || value > upper) {
    sab_stop(name, " must be finite and in [", lower, ", ", upper, "].")
  }
  value
}

sab_sha256_file <- function(path) {
  path <- normalizePath(path, mustWork = TRUE)
  executable <- Sys.which("sha256sum")
  arguments <- c("--", shQuote(path))
  if (!nzchar(executable)) {
    executable <- Sys.which("shasum")
    arguments <- c("-a", "256", shQuote(path))
  }
  if (!nzchar(executable)) sab_stop("No SHA-256 executable is available.")
  output <- suppressWarnings(system2(
    executable, arguments, stdout = TRUE, stderr = TRUE
  ))
  status <- attr(output, "status")
  if ((!is.null(status) && status != 0L) || length(output) != 1L) {
    sab_stop("Could not hash ", path, ".")
  }
  digest <- tolower(strsplit(trimws(output), "[[:space:]]+")[[1L]][[1L]])
  if (!grepl("^[0-9a-f]{64}$", digest)) {
    sab_stop("Invalid SHA-256 output for ", path, ".")
  }
  digest
}

sab_assert_named_vector <- function(value, expected_names, label) {
  original_names <- names(value)
  if ((!is.numeric(value) && !is.integer(value)) ||
      length(value) != length(expected_names) ||
      !identical(original_names, expected_names) ||
      any(!is.finite(value))) {
    sab_stop(label, " must be finite and named in canonical order.")
  }
  value <- as.numeric(value)
  names(value) <- expected_names
  value
}

sab_log_mean_exp <- function(value) {
  if (!is.numeric(value) || !length(value) || anyNA(value) ||
      any(is.nan(value))) {
    return(NA_real_)
  }
  if (any(value == Inf)) return(Inf)
  maximum <- max(value)
  if (maximum == -Inf) return(-Inf)
  maximum + log(mean(exp(value - maximum)))
}

sab_effective_sample_size <- function(value) {
  value <- as.numeric(value)
  n <- length(value)
  if (n < 4L || any(!is.finite(value)) || stats::var(value) <= 0) return(0)
  maximum_lag <- min(n - 1L, max(20L, floor(n / 2L)))
  correlation <- as.numeric(stats::acf(
    value, lag.max = maximum_lag, plot = FALSE,
    type = "correlation", demean = TRUE
  )$acf)[-1L]
  if (length(correlation) < 2L) return(n)
  correlation <- correlation[seq_len(2L * floor(length(correlation) / 2L))]
  pair_sums <- correlation[seq(1L, length(correlation), by = 2L)] +
    correlation[seq(2L, length(correlation), by = 2L)]
  first_nonpositive <- which(!is.finite(pair_sums) | pair_sums <= 0)[1L]
  if (!is.na(first_nonpositive)) {
    pair_sums <- head(pair_sums, first_nonpositive - 1L)
  }
  if (length(pair_sums) > 1L) {
    pair_sums <- cummin(pair_sums)
  }
  integrated_time <- max(1, 1 + 2 * sum(pair_sums))
  min(n, n / integrated_time)
}

sab_read_cmdstan_retained <- function(path, expected_sha256, selected,
                                      warmup, retained) {
  if (!file.exists(path)) sab_stop("Missing comparator CSV: ", path)
  observed <- sab_sha256_file(path)
  if (!identical(observed, expected_sha256)) {
    sab_stop("Comparator hash mismatch for ", path, ".")
  }
  connection <- file(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  header <- NULL
  repeat {
    line <- readLines(connection, n = 1L, warn = FALSE)
    if (!length(line)) break
    if (nzchar(line) && substr(line, 1L, 1L) != "#") {
      header <- line
      break
    }
  }
  close(connection)
  on.exit(NULL, add = FALSE)
  if (is.null(header)) sab_stop("Comparator has no CSV header: ", path)
  fields <- strsplit(header, ",", fixed = TRUE)[[1L]]
  if (anyDuplicated(fields)) sab_stop("Comparator CSV has duplicate fields.")
  absent <- setdiff(selected, fields)
  if (length(absent)) {
    sab_stop("Comparator omits selected fields: ", paste(absent, collapse = ", "))
  }
  classes <- rep("NULL", length(fields))
  classes[match(selected, fields)] <- "numeric"
  draws <- utils::read.csv(
    path,
    comment.char = "#",
    check.names = FALSE,
    colClasses = classes,
    stringsAsFactors = FALSE
  )
  if (nrow(draws) != warmup + retained ||
      !setequal(names(draws), selected) ||
      any(!is.finite(as.matrix(draws)))) {
    sab_stop("Comparator row count, selected field order, or values are invalid.")
  }
  draws <- draws[, selected, drop = FALSE]
  draws[seq.int(warmup + 1L, warmup + retained), , drop = FALSE]
}

sab_oracle_anchor <- function(workspace_root, adapter, panel) {
  comparator_directory <- file.path(
    workspace_root, "projects", "full_joint_model", "stan", "outputs",
    "hiv_current_1000w1000s_20260506_120805"
  )
  comparator_sha256 <- c(
    "18b96caccb71b4dd5391286889286451e7b5c075622993355bddfe6e283e54d7",
    "ec57dcde7563b9a5ecfb165cdfccc8808d08b59dc2b01ff48af32fbc15e83ee2",
    "46240f79d182078869e9af802e214f2ae97327b6baa0f32b00b92a61ef41a198",
    "905a501289cb65c28ea086593d8282fc26f3a152719b7e71fd0280e933319e59"
  )
  eta_location <- c(
    "mu_log_lambda", "mu_log_mu_t", "mu_log_mu_a", "mu_log_p",
    "mu_log_alpha_l", "mu_u_pi", "mu_u_eta_rti", "mu_u_eta_pi"
  )
  eta_scale_natural <- c(
    "omega_lambda", "omega_mu_t", "omega_mu_a", "omega_p",
    "omega_alpha_l", "omega_pi", "omega_eta_rti", "omega_eta_pi"
  )
  psi_natural <- c("gamma_pop", "mu_l_pop", "sigma_v", "sigma_t")
  local_natural <- c(
    "lambda", "mu_t", "mu_a", "p", "alpha_l", "pi", "eta_rti",
    "eta_pi"
  )
  local_fields <- unlist(lapply(panel$patient_index, function(index) {
    paste0(local_natural, ".", index)
  }), use.names = FALSE)
  selected <- c(
    eta_location, "beta_nelf", eta_scale_natural, psi_natural, local_fields
  )
  selected <- unique(selected)

  eta_by_chain <- vector("list", 4L)
  psi_by_chain <- vector("list", 4L)
  x_by_chain <- vector("list", 4L)
  comparator_paths <- character(4L)
  for (chain in seq_len(4L)) {
    path <- file.path(
      comparator_directory,
      sprintf("hiv_current_1000w1000s_20260506_120805_chain%02d-1.csv", chain)
    )
    comparator_paths[[chain]] <- path
    draw <- sab_read_cmdstan_retained(
      path, comparator_sha256[[chain]], selected,
      warmup = 1000L, retained = 1000L
    )
    if (any(as.matrix(draw[, c(eta_scale_natural, psi_natural)]) <= 0)) {
      sab_stop("Comparator contains non-positive natural scale parameters.")
    }
    eta <- cbind(
      as.matrix(draw[, eta_location, drop = FALSE]),
      beta_nelf = draw$beta_nelf,
      log(as.matrix(draw[, eta_scale_natural, drop = FALSE]))
    )
    colnames(eta) <- adapter$coordinate_names$population
    psi <- log(as.matrix(draw[, psi_natural, drop = FALSE]))
    colnames(psi) <- adapter$coordinate_names$global
    patient_x <- setNames(lapply(seq_len(nrow(panel)), function(row) {
      fields <- paste0(local_natural, ".", panel$patient_index[[row]])
      natural <- as.matrix(draw[, fields, drop = FALSE])
      if (any(natural[, seq_len(5L), drop = FALSE] <= 0) ||
          any(natural[, 6:8, drop = FALSE] <= 0) ||
          any(natural[, 6:8, drop = FALSE] >= 1)) {
        sab_stop("Comparator local state is outside the working-scale domain.")
      }
      x <- cbind(log(natural[, seq_len(5L), drop = FALSE]),
                 stats::qlogis(natural[, 6:8, drop = FALSE]))
      colnames(x) <- adapter$coordinate_names$local
      x
    }), panel$patient_id)
    eta_by_chain[[chain]] <- eta
    psi_by_chain[[chain]] <- psi
    x_by_chain[[chain]] <- patient_x
  }

  eta_draws <- do.call(rbind, eta_by_chain)
  psi_draws <- do.call(rbind, psi_by_chain)
  x_draws <- setNames(lapply(panel$patient_id, function(patient_id) {
    do.call(rbind, lapply(x_by_chain, `[[`, patient_id))
  }), panel$patient_id)
  eta_center <- colMeans(eta_draws)
  psi_center <- colMeans(psi_draws)
  x_center <- lapply(x_draws, colMeans)
  eta_center <- sab_assert_named_vector(
    eta_center, adapter$coordinate_names$population, "oracle eta center"
  )
  psi_center <- sab_assert_named_vector(
    psi_center, adapter$coordinate_names$global, "oracle psi center"
  )
  x_center <- lapply(x_center, sab_assert_named_vector,
                     expected_names = adapter$coordinate_names$local,
                     label = "oracle patient center")

  dynamic_names <- adapter$coordinate_names$dynamic_global
  dynamic_covariance <- stats::cov(psi_draws[, dynamic_names, drop = FALSE])
  dynamic_chol <- t(chol(dynamic_covariance))
  dimnames(dynamic_chol) <- list(dynamic_names, dynamic_names)
  # This published audit value catches accidental inclusion of warmup or use
  # of natural rather than log coordinates without defining the estimator.
  audited_chol <- rbind(
    c(0.3705365474, 0),
    c(0.00878601058, 0.1169029412)
  )
  dimnames(audited_chol) <- dimnames(dynamic_chol)
  audit_difference <- max(abs(dynamic_chol - audited_chol))
  if (!is.finite(audit_difference) || audit_difference > 5e-6) {
    sab_stop(
      "Derived dynamic-psi Cholesky disagrees with the pinned comparator audit: ",
      format(audit_difference, scientific = TRUE), "."
    )
  }
  list(
    schema_version = "sab_system_a_oracle_anchor_v1",
    declaration = paste(
      "ORACLE comparator anchor/scale only; not a SAEM estimate and not",
      "deployable evidence"
    ),
    eta = eta_center,
    psi = psi_center,
    x = x_center,
    dynamic_covariance = dynamic_covariance,
    dynamic_chol = dynamic_chol,
    audited_chol_max_abs_difference = audit_difference,
    comparator_paths = comparator_paths,
    comparator_sha256 = comparator_sha256,
    retained_draws = nrow(psi_draws),
    # Kept only in memory to find a finite deterministic start if a nonlinear
    # average is outside the equilibrium domain. These are stripped from the
    # final anchor artifact below.
    candidate_dynamic = psi_draws[, dynamic_names, drop = FALSE],
    candidate_x = x_draws
  )
}

sab_initial_candidates <- function(anchor, patient_id, psi, maximum = 50L) {
  dynamic_names <- colnames(anchor$candidate_dynamic)
  delta <- sweep(
    anchor$candidate_dynamic, 2L, psi[dynamic_names], FUN = "-"
  )
  whitened <- t(forwardsolve(anchor$dynamic_chol, t(delta)))
  nearest <- order(rowSums(whitened^2))
  nearest <- head(nearest, maximum)
  candidates <- rbind(anchor$x[[patient_id]],
                      anchor$candidate_x[[patient_id]][nearest, , drop = FALSE])
  colnames(candidates) <- names(anchor$x[[patient_id]])
  candidates
}

sab_evaluate_patient <- function(adapter, patient_id, x, psi) {
  prediction <- adapter$solve_prediction(patient_id, x, psi)
  if (!isTRUE(prediction$ok)) {
    return(list(
      loglik = -Inf,
      solver_failure = TRUE,
      reason = as.character(prediction$reason %||% "unknown_failure")
    ))
  }
  loglik <- adapter$loglik_from_prediction(prediction, psi)
  if (length(loglik) != 1L || is.na(loglik) || is.nan(loglik) || loglik == Inf) {
    sab_stop("Certified likelihood returned an invalid value.")
  }
  list(
    loglik = as.numeric(loglik),
    solver_failure = FALSE,
    reason = if (is.finite(loglik)) "ok" else "nonfinite_loglik"
  )
}

`%||%` <- function(left, right) {
  if (is.null(left) || !length(left) || is.na(left[[1L]])) right else left
}

sab_run_patient_conditional <- function(adapter, patient_id, eta, psi,
                                        initial_candidates, warmup, draws, thin,
                                        initial_beta, target_acceptance,
                                        adaptation_block, seed) {
  set.seed(seed)
  local_names <- adapter$coordinate_names$local
  if (is.numeric(initial_candidates) && is.null(dim(initial_candidates))) {
    initial_candidates <- matrix(
      initial_candidates, nrow = 1L,
      dimnames = list(NULL, names(initial_candidates))
    )
  }
  if (!is.matrix(initial_candidates) || !is.numeric(initial_candidates) ||
      nrow(initial_candidates) < 1L ||
      !identical(colnames(initial_candidates), local_names) ||
      any(!is.finite(initial_candidates))) {
    sab_stop("Initial patient candidates are malformed.")
  }
  population_mean <- adapter$population_mean(patient_id, eta)
  population_mean <- sab_assert_named_vector(
    population_mean, local_names, "patient population mean"
  )
  log_scale_names <- grep(
    "^log_omega_", adapter$coordinate_names$population, value = TRUE
  )
  if (length(log_scale_names) != length(local_names)) {
    sab_stop("System A population scales do not match local coordinates.")
  }
  population_sd <- exp(eta[log_scale_names])
  names(population_sd) <- local_names
  if (any(!is.finite(population_sd)) || any(population_sd <= 0)) {
    sab_stop("Patient population scale is not representable.")
  }
  current <- NULL
  x <- NULL
  initialization_ode_calls <- 0L
  initialization_solver_failures <- 0L
  for (candidate_index in seq_len(nrow(initial_candidates))) {
    candidate <- initial_candidates[candidate_index, ]
    names(candidate) <- local_names
    evaluated <- sab_evaluate_patient(adapter, patient_id, candidate, psi)
    initialization_ode_calls <- initialization_ode_calls + 1L
    initialization_solver_failures <- initialization_solver_failures +
      as.integer(evaluated$solver_failure)
    if (is.finite(evaluated$loglik)) {
      x <- candidate
      current <- evaluated
      break
    }
  }
  if (is.null(current)) {
    sab_stop("No pinned-comparator candidate has finite conditional density: ",
             patient_id, ".")
  }

  beta <- initial_beta
  beta_logit <- stats::qlogis(beta)
  total_sampling <- draws * thin
  total_iterations <- warmup + total_sampling
  kept <- matrix(NA_real_, nrow = draws, ncol = length(local_names),
                 dimnames = list(NULL, local_names))
  kept_loglik <- numeric(draws)
  warmup_accepted <- 0L
  sampling_accepted <- 0L
  block_accepted <- 0L
  solver_failures <- 0L
  nonfinite_loglik <- 0L
  kept_index <- 0L
  adaptation_index <- 0L

  for (iteration in seq_len(total_iterations)) {
    persistence <- sqrt(1 - beta^2)
    proposed_x <- population_mean +
      persistence * (x - population_mean) +
      beta * population_sd * stats::rnorm(length(x))
    names(proposed_x) <- local_names
    proposed <- sab_evaluate_patient(adapter, patient_id, proposed_x, psi)
    solver_failures <- solver_failures + as.integer(proposed$solver_failure)
    nonfinite_loglik <- nonfinite_loglik + as.integer(!is.finite(proposed$loglik))
    accepted <- log(stats::runif(1L)) < proposed$loglik - current$loglik
    if (accepted) {
      x <- proposed_x
      current <- proposed
    }

    if (iteration <= warmup) {
      warmup_accepted <- warmup_accepted + as.integer(accepted)
      block_accepted <- block_accepted + as.integer(accepted)
      if (iteration %% adaptation_block == 0L) {
        adaptation_index <- adaptation_index + 1L
        observed <- block_accepted / adaptation_block
        gain <- 0.8 / sqrt(adaptation_index)
        beta_logit <- beta_logit + gain * (observed - target_acceptance)
        beta <- min(0.95, max(0.005, stats::plogis(beta_logit)))
        beta_logit <- stats::qlogis(beta)
        block_accepted <- 0L
      }
    } else {
      sampling_accepted <- sampling_accepted + as.integer(accepted)
      if ((iteration - warmup) %% thin == 0L) {
        kept_index <- kept_index + 1L
        kept[kept_index, ] <- x
        kept_loglik[[kept_index]] <- current$loglik
      }
    }
  }
  if (kept_index != draws || any(!is.finite(kept)) ||
      any(!is.finite(kept_loglik))) {
    sab_stop("Conditional bank was not filled with finite states.")
  }
  coordinate_ess <- apply(kept, 2L, sab_effective_sample_size)
  list(
    draws = kept,
    loglik = kept_loglik,
    final_x = x,
    final_beta = beta,
    warmup_acceptance = warmup_accepted / warmup,
    sampling_acceptance = sampling_accepted / total_sampling,
    coordinate_ess = coordinate_ess,
    solver_failures = solver_failures,
    nonfinite_loglik = nonfinite_loglik,
    initialization_ode_calls = initialization_ode_calls,
    initialization_solver_failures = initialization_solver_failures,
    exact_ode_calls = total_iterations + initialization_ode_calls,
    seed = seed
  )
}

sab_repaired_cholesky <- function(draws, relative_floor) {
  covariance <- stats::cov(draws)
  decomposition <- eigen(covariance, symmetric = TRUE)
  maximum <- max(decomposition$values)
  if (!is.finite(maximum) || maximum <= 0) {
    sab_stop("Conditional covariance has no positive eigenvalue.")
  }
  floor <- max(1e-12, relative_floor * maximum)
  repaired_values <- pmax(decomposition$values, floor)
  repaired <- decomposition$vectors %*%
    (repaired_values * t(decomposition$vectors))
  repaired <- (repaired + t(repaired)) / 2
  list(
    chol = t(chol(repaired)),
    raw_min_eigenvalue = min(decomposition$values),
    applied_floor = floor,
    regularized = any(decomposition$values < floor)
  )
}

sab_fit_patient_map <- function(patient_id, banks, design_z, anchor,
                                adapter, covariance_floor) {
  local_names <- adapter$coordinate_names$local
  expected_design <- c(
    "center", "plus_axis_1", "minus_axis_1", "plus_axis_2", "minus_axis_2"
  )
  if (!identical(names(banks), expected_design) ||
      !identical(rownames(design_z), expected_design)) {
    sab_stop("Patient map banks do not follow the frozen axial design.")
  }
  means <- do.call(rbind, lapply(banks, function(bank) {
    colMeans(bank$draws)
  }))
  radius <- abs(design_z["plus_axis_1", "z1"])
  if (!is.finite(radius) || radius <= 0 ||
      design_z["minus_axis_1", "z1"] != -radius ||
      design_z["plus_axis_2", "z2"] != radius ||
      design_z["minus_axis_2", "z2"] != -radius) {
    sab_stop("Patient map design is not symmetric around its anchor.")
  }
  mean_coefficients <- cbind(
    intercept = means["center", ],
    z1 = (means["plus_axis_1", ] - means["minus_axis_1", ]) /
      (2 * radius),
    z2 = (means["plus_axis_2", ] - means["minus_axis_2", ]) /
      (2 * radius)
  )
  rownames(mean_coefficients) <- local_names

  cholesky <- lapply(banks, function(bank) {
    sab_repaired_cholesky(bank$draws, covariance_floor)
  })
  p <- length(local_names)
  lower_index <- which(lower.tri(matrix(0, p, p), diag = TRUE), arr.ind = TRUE)
  features <- t(vapply(cholesky, function(item) {
    values <- item$chol[lower_index]
    diagonal <- lower_index[, 1L] == lower_index[, 2L]
    values[diagonal] <- log(values[diagonal])
    values
  }, numeric(nrow(lower_index))))
  chol_coefficients <- cbind(
    intercept = features["center", ],
    z1 = (features["plus_axis_1", ] - features["minus_axis_1", ]) /
      (2 * radius),
    z2 = (features["plus_axis_2", ] - features["minus_axis_2", ]) /
      (2 * radius)
  )

  dynamic_names <- adapter$coordinate_names$dynamic_global
  dynamic_center <- anchor$psi[dynamic_names]
  dynamic_chol <- anchor$dynamic_chol
  map <- local({
    mean_coef <- mean_coefficients
    chol_coef <- chol_coefficients
    index <- lower_index
    center <- dynamic_center
    global_chol <- dynamic_chol
    coordinate_names <- local_names
    id <- patient_id
    sab_new_affine_map(
      mean_fn = function(global) {
        global <- sab_assert_named_vector(global, dynamic_names,
                                          "dynamic global state")
        z <- as.numeric(forwardsolve(global_chol, global - center))
        value <- as.numeric(mean_coef %*% c(1, z))
        names(value) <- coordinate_names
        value
      },
      chol_fn = function(global) {
        global <- sab_assert_named_vector(global, dynamic_names,
                                          "dynamic global state")
        z <- as.numeric(forwardsolve(global_chol, global - center))
        values <- as.numeric(chol_coef %*% c(1, z))
        diagonal <- index[, 1L] == index[, 2L]
        values[diagonal] <- exp(values[diagonal])
        value <- matrix(0, p, p)
        value[index] <- values
        value
      },
      name = paste0("system_a_patient_", id),
      check_at = dynamic_center
    )
  })
  list(
    map = map,
    coefficients = list(
      mean = mean_coefficients,
      chol = chol_coefficients,
      lower_index = lower_index
    ),
    covariance_diagnostics = data.frame(
      design = rownames(design_z),
      raw_min_eigenvalue = vapply(
        cholesky, `[[`, numeric(1L), "raw_min_eigenvalue"
      ),
      applied_floor = vapply(cholesky, `[[`, numeric(1L), "applied_floor"),
      regularized = vapply(cholesky, `[[`, logical(1L), "regularized"),
      stringsAsFactors = FALSE
    )
  )
}

sab_bank_diagnostic_row <- function(patient_id, bank_name, purpose, bank) {
  data.frame(
    patient_id = patient_id,
    bank = bank_name,
    purpose = purpose,
    warmup_acceptance = bank$warmup_acceptance,
    sampling_acceptance = bank$sampling_acceptance,
    final_beta = bank$final_beta,
    minimum_coordinate_ess = min(bank$coordinate_ess),
    median_coordinate_ess = stats::median(bank$coordinate_ess),
    solver_failure_rate = bank$solver_failures /
      (bank$exact_ode_calls - bank$initialization_ode_calls),
    nonfinite_loglik_rate = bank$nonfinite_loglik /
      (bank$exact_ode_calls - bank$initialization_ode_calls),
    exact_ode_calls = bank$exact_ode_calls,
    initialization_ode_calls = bank$initialization_ode_calls,
    initialization_solver_failures = bank$initialization_solver_failures,
    seed = bank$seed,
    stringsAsFactors = FALSE
  )
}

sab_patient_endpoints <- function(adapter, patient_id, eta, psi_center,
                                  dynamic_center, moves, heldout, map) {
  n_draw <- nrow(heldout$draws)
  n_move <- nrow(moves)
  methods <- c("fixed_x", "affine")
  output <- setNames(lapply(methods, function(method) {
    list(
      work = matrix(NA_real_, nrow = n_draw, ncol = n_move),
      solver_failure = matrix(FALSE, nrow = n_draw, ncol = n_move),
      nonfinite_loglik = matrix(FALSE, nrow = n_draw, ncol = n_move)
    )
  }), methods)
  current_population <- apply(heldout$draws, 1L, function(x) {
    names(x) <- adapter$coordinate_names$local
    adapter$log_population_density(patient_id, x, eta)
  })
  if (any(!is.finite(current_population))) {
    sab_stop("Held-out patient bank has non-finite population density.")
  }
  current_logdet <- sab_affine_log_abs_det(map, dynamic_center)

  for (move_index in seq_len(n_move)) {
    proposed_dynamic <- c(
      log_gamma_pop = moves$log_gamma_pop[[move_index]],
      log_mu_l_pop = moves$log_mu_l_pop[[move_index]]
    )
    proposed_psi <- psi_center
    proposed_psi[names(proposed_dynamic)] <- proposed_dynamic
    proposed_logdet <- sab_affine_log_abs_det(map, proposed_dynamic)
    for (draw_index in seq_len(n_draw)) {
      current_x <- heldout$draws[draw_index, ]
      names(current_x) <- adapter$coordinate_names$local
      for (method in methods) {
        proposed_x <- if (method == "fixed_x") {
          current_x
        } else {
          sab_affine_transport(
            map, current_x, dynamic_center, proposed_dynamic
          )
        }
        evaluation <- sab_evaluate_patient(
          adapter, patient_id, proposed_x, proposed_psi
        )
        output[[method]]$solver_failure[draw_index, move_index] <-
          evaluation$solver_failure
        output[[method]]$nonfinite_loglik[draw_index, move_index] <-
          !is.finite(evaluation$loglik)
        if (method == "fixed_x") {
          output[[method]]$work[draw_index, move_index] <-
            evaluation$loglik - heldout$loglik[[draw_index]]
        } else {
          proposed_population <- adapter$log_population_density(
            patient_id, proposed_x, eta
          )
          output[[method]]$work[draw_index, move_index] <-
            evaluation$loglik + proposed_population + proposed_logdet -
            heldout$loglik[[draw_index]] -
            current_population[[draw_index]] - current_logdet
        }
      }
    }
  }
  output
}

sab_work_diagnostics <- function(patient_results, moves) {
  rows <- list()
  row_index <- 0L
  for (patient_id in names(patient_results)) {
    for (method in names(patient_results[[patient_id]])) {
      work_matrix <- patient_results[[patient_id]][[method]]$work
      for (move_index in seq_len(nrow(moves))) {
        work <- work_matrix[, move_index]
        finite <- is.finite(work)
        log_normalizer <- sab_log_mean_exp(work)
        forward_d2 <- sab_log_mean_exp(2 * work) - 2 * log_normalizer
        # Change of measure gives D2(s_c || s_c') using the same independent
        # centre bank: log E_c exp(-W) + log E_c exp(W).
        reverse_d2 <- sab_log_mean_exp(-work) + log_normalizer
        midpoint <- floor(length(work) / 2L)
        halves <- list(
          first = seq_len(midpoint),
          second = seq.int(midpoint + 1L, length(work))
        )
        half_forward <- vapply(halves, function(index) {
          sab_log_mean_exp(2 * work[index]) -
            2 * sab_log_mean_exp(work[index])
        }, numeric(1L))
        half_reverse <- vapply(halves, function(index) {
          sab_log_mean_exp(-work[index]) + sab_log_mean_exp(work[index])
        }, numeric(1L))
        row_index <- row_index + 1L
        rows[[row_index]] <- data.frame(
          patient_id = patient_id,
          method = method,
          move_id = moves$move_id[[move_index]],
          scale = moves$scale[[move_index]],
          axis = moves$axis[[move_index]],
          sign = moves$sign[[move_index]],
          finite_work_fraction = mean(finite),
          work_variance_finite = if (sum(finite) >= 2L) {
            stats::var(work[finite])
          } else {
            NA_real_
          },
          d2_forward = forward_d2,
          d2_reverse_importance = reverse_d2,
          d2_forward_half_difference = abs(diff(half_forward)),
          d2_reverse_half_difference = abs(diff(half_reverse)),
          stringsAsFactors = FALSE
        )
      }
    }
  }
  do.call(rbind, rows)
}

sab_aggregate_work_diagnostics <- function(patient_diagnostics) {
  groups <- split(
    patient_diagnostics,
    interaction(
      patient_diagnostics$method,
      patient_diagnostics$move_id,
      drop = TRUE
    )
  )
  output <- do.call(rbind, lapply(groups, function(group) {
    all_patient_work_finite <- all(group$finite_work_fraction == 1)
    variance_finite <- all(is.finite(group$work_variance_finite))
    forward_d2_finite <- all(is.finite(group$d2_forward))
    reverse_d2_finite <- all(is.finite(group$d2_reverse_importance))
    total_variance <- if (variance_finite) {
      sum(group$work_variance_finite)
    } else {
      NA_real_
    }
    total_forward_d2 <- if (forward_d2_finite) {
      sum(group$d2_forward)
    } else {
      NA_real_
    }
    total_reverse_d2 <- if (reverse_d2_finite) {
      sum(group$d2_reverse_importance)
    } else {
      NA_real_
    }
    data.frame(
      method = group$method[[1L]],
      move_id = group$move_id[[1L]],
      scale = group$scale[[1L]],
      axis = group$axis[[1L]],
      sign = group$sign[[1L]],
      all_patient_work_finite = all_patient_work_finite,
      all_diagnostics_finite = variance_finite && forward_d2_finite &&
        reverse_d2_finite,
      total_work_variance = total_variance,
      total_d2_forward = total_forward_d2,
      total_d2_reverse_importance = total_reverse_d2,
      forward_weight_relative_ess = if (is.finite(total_forward_d2)) {
        exp(-total_forward_d2)
      } else {
        NA_real_
      },
      reverse_weight_relative_ess = if (is.finite(total_reverse_d2)) {
        exp(-total_reverse_d2)
      } else {
        NA_real_
      },
      maximum_forward_half_difference =
        max(group$d2_forward_half_difference),
      maximum_reverse_half_difference =
        max(group$d2_reverse_half_difference),
      stringsAsFactors = FALSE
    )
  }))
  rownames(output) <- NULL
  output
}

sab_replay_summary <- function(adapter, eta, psi_center, patient_results,
                               moves, seed) {
  methods <- names(patient_results[[1L]])
  n_draw <- nrow(patient_results[[1L]][[1L]]$work)
  n_patient <- length(patient_results)
  set.seed(seed)
  common_log_uniform <- matrix(
    log(stats::runif(n_draw * nrow(moves))),
    nrow = n_draw,
    ncol = nrow(moves)
  )
  current_hyperprior <- adapter$log_hyperprior(eta, psi_center)
  if (!is.finite(current_hyperprior)) {
    sab_stop("Oracle anchor has non-finite System A hyperprior density.")
  }
  rows <- list()
  row_index <- 0L
  draw_rows <- list()
  draw_row_index <- 0L
  for (method in methods) {
    for (move_index in seq_len(nrow(moves))) {
      proposed_psi <- psi_center
      proposed_psi[c("log_gamma_pop", "log_mu_l_pop")] <- c(
        moves$log_gamma_pop[[move_index]],
        moves$log_mu_l_pop[[move_index]]
      )
      hyperprior_difference <-
        adapter$log_hyperprior(eta, proposed_psi) - current_hyperprior
      total_work <- Reduce(`+`, lapply(patient_results, function(patient) {
        patient[[method]]$work[, move_index]
      }))
      log_ratio <- hyperprior_difference + total_work
      log_acceptance <- pmin(0, log_ratio)
      acceptance_probability <- ifelse(
        is.finite(log_acceptance), exp(log_acceptance), 0
      )
      accepted <- common_log_uniform[, move_index] < log_acceptance
      solver_failure <- Reduce(`+`, lapply(patient_results, function(patient) {
        as.integer(patient[[method]]$solver_failure[, move_index])
      }))
      nonfinite <- Reduce(`+`, lapply(patient_results, function(patient) {
        as.integer(patient[[method]]$nonfinite_loglik[, move_index])
      }))
      squared_jump <- moves$scale[[move_index]]^2
      row_index <- row_index + 1L
      rows[[row_index]] <- data.frame(
        method = method,
        move_id = moves$move_id[[move_index]],
        scale = moves$scale[[move_index]],
        axis = moves$axis[[move_index]],
        sign = moves$sign[[move_index]],
        expected_acceptance = mean(acceptance_probability),
        realised_acceptance = mean(accepted),
        expected_whitened_esjd = squared_jump * mean(acceptance_probability),
        realised_whitened_esjd = squared_jump * mean(accepted),
        expected_esjd_per_ode = squared_jump *
          mean(acceptance_probability) / n_patient,
        realised_esjd_per_ode = squared_jump * mean(accepted) / n_patient,
        proposed_exact_ode_calls = n_draw * n_patient,
        solver_failure_rate = sum(solver_failure) / (n_draw * n_patient),
        nonfinite_loglik_rate = sum(nonfinite) / (n_draw * n_patient),
        stringsAsFactors = FALSE
      )
      draw_row_index <- draw_row_index + 1L
      draw_rows[[draw_row_index]] <- data.frame(
        method = method,
        move_id = moves$move_id[[move_index]],
        draw = seq_len(n_draw),
        log_ratio = log_ratio,
        acceptance_probability = acceptance_probability,
        accepted = accepted,
        patient_solver_failures = solver_failure,
        patient_nonfinite_loglik = nonfinite,
        stringsAsFactors = FALSE
      )
    }
  }
  by_move <- do.call(rbind, rows)
  by_scale <- do.call(rbind, lapply(split(
    by_move, interaction(by_move$method, by_move$scale, drop = TRUE)
  ), function(group) {
    data.frame(
      method = group$method[[1L]],
      scale = group$scale[[1L]],
      expected_acceptance = mean(group$expected_acceptance),
      realised_acceptance = mean(group$realised_acceptance),
      expected_whitened_esjd = mean(group$expected_whitened_esjd),
      realised_whitened_esjd = mean(group$realised_whitened_esjd),
      expected_esjd_per_ode = mean(group$expected_esjd_per_ode),
      realised_esjd_per_ode = mean(group$realised_esjd_per_ode),
      solver_failure_rate = mean(group$solver_failure_rate),
      nonfinite_loglik_rate = mean(group$nonfinite_loglik_rate),
      proposed_exact_ode_calls = sum(group$proposed_exact_ode_calls),
      stringsAsFactors = FALSE
    )
  }))
  rownames(by_scale) <- NULL
  list(
    by_move = by_move,
    by_scale = by_scale,
    draws = do.call(rbind, draw_rows)
  )
}

command_args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(command_args) == 1L) command_args[[1L]] else {
  if (length(command_args) > 1L) sab_stop("Expected at most one output path.")
  file.path("outputs", "system_a_12_oracle_v1")
}
workspace_root <- Sys.getenv("SAB_WORKSPACE_ROOT", unset = "")
if (!nzchar(workspace_root) || !dir.exists(workspace_root)) {
  sab_stop("SAB_WORKSPACE_ROOT must name the workspace containing projects/.")
}
workspace_root <- normalizePath(workspace_root, mustWork = TRUE)
if (dir.exists(output_dir) && length(list.files(
    output_dir, all.files = TRUE, no.. = TRUE
))) {
  sab_stop("Output directory is not empty: ", output_dir)
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path("R", "affine_map.R"), local = FALSE)
source(file.path("R", "transport_mh.R"), local = FALSE)
source(file.path("R", "system_a_adapter.R"), local = FALSE)

config <- list(
  schema_version = "sab_system_a_12_oracle_config_v1",
  declaration = paste(
    "ORACLE 12-patient endpoint replay; comparator anchor/scale;",
    "not SAEM, not a retained chain, not deployable proof"
  ),
  map_radius = sab_numeric_env("SAB_SYSA_MAP_RADIUS", 1.0, 0.1, 2.0),
  map_warmup = sab_integer_env("SAB_SYSA_MAP_WARMUP", 750L, 50L),
  map_draws = sab_integer_env("SAB_SYSA_MAP_DRAWS", 1000L, 100L),
  map_thin = sab_integer_env("SAB_SYSA_MAP_THIN", 1L, 1L),
  heldout_warmup = sab_integer_env("SAB_SYSA_HELDOUT_WARMUP", 750L, 50L),
  heldout_draws = sab_integer_env("SAB_SYSA_HELDOUT_DRAWS", 500L, 100L),
  heldout_thin = sab_integer_env("SAB_SYSA_HELDOUT_THIN", 1L, 1L),
  initial_pcn_beta = sab_numeric_env(
    "SAB_SYSA_INITIAL_PCN_BETA", 0.20, 0.005, 0.95
  ),
  pcn_target_acceptance = sab_numeric_env(
    "SAB_SYSA_PCN_TARGET", 0.30, 0.05, 0.95
  ),
  adaptation_block = sab_integer_env("SAB_SYSA_ADAPT_BLOCK", 25L, 5L),
  covariance_relative_floor = sab_numeric_env(
    "SAB_SYSA_COVARIANCE_FLOOR", 1e-6, 1e-12, 1e-2
  ),
  proposal_scales = c(0.5, 1.0),
  minimum_bank_ess = sab_numeric_env("SAB_SYSA_MIN_BANK_ESS", 25, 1, Inf),
  minimum_efficiency_ratio = sab_numeric_env(
    "SAB_SYSA_MIN_EFFICIENCY_RATIO", 1.5, 1, Inf
  ),
  minimum_affine_acceptance = sab_numeric_env(
    "SAB_SYSA_MIN_AFFINE_ACCEPTANCE", 0.15, 0, 1
  ),
  maximum_solver_failure_rate = sab_numeric_env(
    "SAB_SYSA_MAX_SOLVER_FAILURE", 0.05, 0, 1
  ),
  map_seed = 202609210L,
  heldout_seed = 202609310L,
  replay_seed = 202609410L,
  cores = sab_integer_env("SAB_SYSTEM_A_CORES", 1L, 1L)
)
if (config$cores > 1L && .Platform$OS.type != "unix") {
  sab_stop("Parallel System A pilot requires Unix fork support.")
}
if (config$map_radius < max(config$proposal_scales)) {
  sab_stop("Map anchors must cover every replay proposal scale.")
}

panel_path <- file.path(
  workspace_root, "projects", "active", "modular_bayes",
  "system_a_failure_audit_v1_20260902_c73810b", "failure_audit_panel.tsv"
)
panel_sha256 <- "7ead58001da0c3e58db64fdb1acaee2e2d17bbbbd2ae388f5eb3fa848eacabf7"
if (!file.exists(panel_path) || !identical(sab_sha256_file(panel_path), panel_sha256)) {
  sab_stop("The predeclared System A audit panel is absent or changed.")
}
panel_original <- utils::read.delim(
  panel_path, stringsAsFactors = FALSE, check.names = FALSE
)
required_panel_fields <- c(
  "patient_index", "patient_id", "treat_nelf", "selection_role"
)
if (length(setdiff(required_panel_fields, names(panel_original))) ||
    nrow(panel_original) != 12L || anyDuplicated(panel_original$patient_index) ||
    anyDuplicated(panel_original$patient_id)) {
  sab_stop("The predeclared System A audit panel is malformed.")
}
panel_original$patient_id <- as.character(panel_original$patient_id)
panel <- panel_original[order(panel_original$patient_index), , drop = FALSE]
rownames(panel) <- NULL
expected_canonical_ids <- c(
  "3", "6", "16", "20", "55", "71", "74", "88", "105", "111",
  "117", "122"
)
if (!identical(panel$patient_id, expected_canonical_ids)) {
  sab_stop("Canonical panel IDs differ from the predeclared audit panel.")
}

adapter <- sab_load_system_a_adapter(
  workspace_root = workspace_root,
  patient_ids = panel$patient_id
)
expected_treatment <- setNames(as.integer(panel$treat_nelf), panel$patient_id)
if (!identical(adapter$treatment, expected_treatment)) {
  sab_stop("Panel treatment coding differs from the certified target.")
}
anchor <- sab_oracle_anchor(workspace_root, adapter, panel)
message("Validated the certified System A target and pinned ORACLE anchor.")
dynamic_names <- adapter$coordinate_names$dynamic_global
dynamic_center <- anchor$psi[dynamic_names]

design_z <- rbind(
  center = c(0, 0),
  plus_axis_1 = c(config$map_radius, 0),
  minus_axis_1 = c(-config$map_radius, 0),
  plus_axis_2 = c(0, config$map_radius),
  minus_axis_2 = c(0, -config$map_radius)
)
colnames(design_z) <- c("z1", "z2")
design_psi <- t(apply(design_z, 1L, function(z) {
  value <- anchor$psi
  value[dynamic_names] <- dynamic_center + anchor$dynamic_chol %*% z
  value
}))
colnames(design_psi) <- adapter$coordinate_names$global

patient_tasks <- parallel::mclapply(
  seq_len(nrow(panel)),
  function(patient_row) {
    patient_id <- panel$patient_id[[patient_row]]
    message("Building conditional banks for patient ", patient_id, ".")
    fit_banks <- setNames(lapply(seq_len(nrow(design_z)), function(index) {
      psi <- design_psi[index, ]
      names(psi) <- adapter$coordinate_names$global
      sab_run_patient_conditional(
        adapter = adapter,
        patient_id = patient_id,
        eta = anchor$eta,
        psi = psi,
        initial_candidates = sab_initial_candidates(
          anchor, patient_id, psi
        ),
        warmup = config$map_warmup,
        draws = config$map_draws,
        thin = config$map_thin,
        initial_beta = config$initial_pcn_beta,
        target_acceptance = config$pcn_target_acceptance,
        adaptation_block = config$adaptation_block,
        seed = config$map_seed + 1000L * patient_row + index
      )
    }), rownames(design_z))
    heldout <- sab_run_patient_conditional(
      adapter = adapter,
      patient_id = patient_id,
      eta = anchor$eta,
      psi = anchor$psi,
      initial_candidates = sab_initial_candidates(
        anchor, patient_id, anchor$psi
      ),
      warmup = config$heldout_warmup,
      draws = config$heldout_draws,
      thin = config$heldout_thin,
      initial_beta = config$initial_pcn_beta,
      target_acceptance = config$pcn_target_acceptance,
      adaptation_block = config$adaptation_block,
      seed = config$heldout_seed + 1000L * patient_row
    )
    message("Completed conditional banks for patient ", patient_id, ".")
    list(patient_id = patient_id, fit_banks = fit_banks, heldout = heldout)
  },
  mc.cores = min(config$cores, nrow(panel)),
  mc.preschedule = TRUE,
  mc.set.seed = FALSE
)
if (any(vapply(patient_tasks, inherits, logical(1L), what = "try-error"))) {
  failures <- vapply(patient_tasks, function(value) {
    if (inherits(value, "try-error")) as.character(value) else ""
  }, character(1L))
  sab_stop("Parallel patient-bank task failed: ",
           paste(failures[nzchar(failures)], collapse = " | "))
}
names(patient_tasks) <- panel$patient_id

fitted <- setNames(lapply(patient_tasks, function(task) {
  sab_fit_patient_map(
    task$patient_id, task$fit_banks, design_z, anchor, adapter,
    config$covariance_relative_floor
  )
}), panel$patient_id)
maps <- lapply(fitted, `[[`, "map")
message("Fitted and froze all twelve patient maps.")

bank_diagnostics <- do.call(rbind, lapply(patient_tasks, function(task) {
  rbind(
    do.call(rbind, lapply(names(task$fit_banks), function(name) {
      sab_bank_diagnostic_row(task$patient_id, name, "map_fit",
                              task$fit_banks[[name]])
    })),
    sab_bank_diagnostic_row(task$patient_id, "center_independent",
                            "heldout_replay", task$heldout)
  )
}))
rownames(bank_diagnostics) <- NULL
covariance_diagnostics <- do.call(rbind, lapply(names(fitted), function(id) {
  cbind(patient_id = id, fitted[[id]]$covariance_diagnostics,
        stringsAsFactors = FALSE)
}))
rownames(covariance_diagnostics) <- NULL

moves <- do.call(rbind, lapply(config$proposal_scales, function(scale) {
  do.call(rbind, lapply(seq_len(2L), function(axis) {
    do.call(rbind, lapply(c(-1, 1), function(sign) {
      z <- numeric(2L)
      z[[axis]] <- sign * scale
      proposed_dynamic <- dynamic_center + anchor$dynamic_chol %*% z
      data.frame(
        move_id = sprintf("z%.2f_axis%d_%s", scale, axis,
                          if (sign < 0) "minus" else "plus"),
        scale = scale,
        axis = axis,
        sign = sign,
        z1 = z[[1L]],
        z2 = z[[2L]],
        log_gamma_pop = proposed_dynamic[[1L]],
        log_mu_l_pop = proposed_dynamic[[2L]],
        stringsAsFactors = FALSE
      )
    }))
  }))
}))
rownames(moves) <- NULL

endpoint_results <- parallel::mclapply(
  seq_len(nrow(panel)),
  function(patient_row) {
    patient_id <- panel$patient_id[[patient_row]]
    message("Replaying exact endpoints for patient ", patient_id, ".")
    value <- sab_patient_endpoints(
      adapter = adapter,
      patient_id = patient_id,
      eta = anchor$eta,
      psi_center = anchor$psi,
      dynamic_center = dynamic_center,
      moves = moves,
      heldout = patient_tasks[[patient_id]]$heldout,
      map = maps[[patient_id]]
    )
    message("Completed exact endpoint replay for patient ", patient_id, ".")
    value
  },
  mc.cores = min(config$cores, nrow(panel)),
  mc.preschedule = TRUE,
  mc.set.seed = FALSE
)
if (any(vapply(endpoint_results, inherits, logical(1L), what = "try-error"))) {
  failures <- vapply(endpoint_results, function(value) {
    if (inherits(value, "try-error")) as.character(value) else ""
  }, character(1L))
  sab_stop("Parallel endpoint task failed: ",
           paste(failures[nzchar(failures)], collapse = " | "))
}
names(endpoint_results) <- panel$patient_id

work_diagnostics <- sab_work_diagnostics(endpoint_results, moves)
aggregate_work_diagnostics <- sab_aggregate_work_diagnostics(work_diagnostics)
replay <- sab_replay_summary(
  adapter, anchor$eta, anchor$psi, endpoint_results, moves,
  seed = config$replay_seed
)

scale_one <- replay$by_scale[replay$by_scale$scale == 1, , drop = FALSE]
fixed_efficiency <- scale_one$expected_esjd_per_ode[
  scale_one$method == "fixed_x"
]
affine_efficiency <- scale_one$expected_esjd_per_ode[
  scale_one$method == "affine"
]
affine_acceptance <- scale_one$expected_acceptance[
  scale_one$method == "affine"
]
affine_failure <- scale_one$solver_failure_rate[
  scale_one$method == "affine"
]
if (length(fixed_efficiency) != 1L || length(affine_efficiency) != 1L ||
    length(affine_acceptance) != 1L || length(affine_failure) != 1L) {
  sab_stop("One-standard-deviation replay summary is incomplete.")
}
efficiency_ratio <- if (fixed_efficiency > 0) {
  affine_efficiency / fixed_efficiency
} else if (affine_efficiency > 0) {
  Inf
} else {
  NA_real_
}
bank_adequate <- all(
  bank_diagnostics$minimum_coordinate_ess >= config$minimum_bank_ess
) && all(bank_diagnostics$sampling_acceptance > 0.01) &&
  all(bank_diagnostics$sampling_acceptance < 0.99)
mechanism_pass <- isTRUE(bank_adequate) && !is.na(efficiency_ratio) &&
  efficiency_ratio >= config$minimum_efficiency_ratio &&
  affine_acceptance >= config$minimum_affine_acceptance &&
  affine_failure <= config$maximum_solver_failure_rate

map_artifact <- list(
  schema_version = "sab_system_a_frozen_map_v1",
  declaration = config$declaration,
  target_fingerprint = adapter$target_fingerprint,
  patient_ids = panel$patient_id,
  dynamic_names = dynamic_names,
  dynamic_center = dynamic_center,
  dynamic_chol = anchor$dynamic_chol,
  fixed_eta = anchor$eta,
  fixed_observation_psi = anchor$psi[
    adapter$coordinate_names$observation_global
  ],
  design_z = design_z,
  coefficients = lapply(fitted, `[[`, "coefficients"),
  covariance_relative_floor = config$covariance_relative_floor
)
anchor_artifact <- anchor
anchor_artifact$candidate_dynamic <- NULL
anchor_artifact$candidate_x <- NULL
result <- list(
  schema_version = "sab_system_a_12_oracle_result_v1",
  completed_at_utc = format(Sys.time(), tz = "UTC", usetz = TRUE),
  slurm_job_id = Sys.getenv("SLURM_JOB_ID", unset = NA_character_),
  session_info = utils::capture.output(sessionInfo()),
  project_source_sha256 = vapply(c(
    "R/affine_map.R", "R/transport_mh.R", "R/system_a_adapter.R",
    "experiments/system_a_12.R"
  ), sab_sha256_file, character(1L)),
  declaration = config$declaration,
  mechanism_pass = mechanism_pass,
  bank_adequate = bank_adequate,
  efficiency_ratio_at_scale_one = efficiency_ratio,
  config = config,
  panel_sha256 = panel_sha256,
  panel_original_order = panel_original,
  panel_canonical_order = panel,
  target = list(
    fingerprint = adapter$target_fingerprint,
    likelihood_signature = adapter$likelihood_signature,
    numerical_target = adapter$numerical_target,
    source_validation = adapter$source_validation,
    oracle_audit = adapter$oracle_audit
  ),
  oracle_anchor = anchor_artifact,
  map = map_artifact,
  bank_diagnostics = bank_diagnostics,
  covariance_diagnostics = covariance_diagnostics,
  moves = moves,
  replay_by_move = replay$by_move,
  replay_by_scale = replay$by_scale,
  replay_draws = replay$draws,
  patient_endpoint_results = endpoint_results,
  work_diagnostics = work_diagnostics,
  aggregate_work_diagnostics = aggregate_work_diagnostics,
  bank_exact_ode_calls = sum(bank_diagnostics$exact_ode_calls),
  replay_proposed_exact_ode_calls = sum(replay$by_move$proposed_exact_ode_calls)
)

utils::write.csv(bank_diagnostics, file.path(output_dir, "bank_diagnostics.csv"),
                 row.names = FALSE)
utils::write.csv(covariance_diagnostics,
                 file.path(output_dir, "covariance_diagnostics.csv"),
                 row.names = FALSE)
utils::write.csv(replay$by_move, file.path(output_dir, "replay_by_move.csv"),
                 row.names = FALSE)
utils::write.csv(replay$by_scale, file.path(output_dir, "replay_by_scale.csv"),
                 row.names = FALSE)
utils::write.csv(work_diagnostics,
                 file.path(output_dir, "patient_work_diagnostics.csv"),
                 row.names = FALSE)
utils::write.csv(
  aggregate_work_diagnostics,
  file.path(output_dir, "aggregate_work_diagnostics.csv"),
  row.names = FALSE
)
saveRDS(map_artifact, file.path(output_dir, "frozen_map.rds"), version = 3)
saveRDS(result, file.path(output_dir, "result.rds"), version = 3)

report <- c(
  "System A affine transport: 12-patient ORACLE mechanism pilot",
  paste0("mechanism_pass: ", toupper(as.character(mechanism_pass))),
  paste0("bank_adequate: ", toupper(as.character(bank_adequate))),
  paste0("declaration: ", config$declaration),
  paste0("target fingerprint: ", adapter$target_fingerprint),
  paste0("upstream commit: ", adapter$source_validation$upstream_commit),
  paste0("canonical patient IDs: ", paste(panel$patient_id, collapse = ",")),
  paste0("panel-role patient IDs: ",
         paste(panel_original$patient_id, collapse = ",")),
  sprintf("minimum map/heldout coordinate ESS: %.1f",
          min(bank_diagnostics$minimum_coordinate_ess)),
  sprintf("sampling acceptance range: %.3f--%.3f",
          min(bank_diagnostics$sampling_acceptance),
          max(bank_diagnostics$sampling_acceptance)),
  sprintf("affine/fixed expected ESJD-per-ODE ratio at zeta=1: %.3f",
          efficiency_ratio),
  sprintf("affine expected acceptance at zeta=1: %.3f", affine_acceptance),
  sprintf("affine solver-failure rate at zeta=1: %.5f", affine_failure),
  sprintf(
    "affine total work-variance range at zeta=1: %.3f--%.3f",
    min(aggregate_work_diagnostics$total_work_variance[
      aggregate_work_diagnostics$method == "affine" &
        aggregate_work_diagnostics$scale == 1
    ]),
    max(aggregate_work_diagnostics$total_work_variance[
      aggregate_work_diagnostics$method == "affine" &
        aggregate_work_diagnostics$scale == 1
    ])
  ),
  sprintf(
    "affine total forward-D2 range at zeta=1: %.3f--%.3f",
    min(aggregate_work_diagnostics$total_d2_forward[
      aggregate_work_diagnostics$method == "affine" &
        aggregate_work_diagnostics$scale == 1
    ]),
    max(aggregate_work_diagnostics$total_d2_forward[
      aggregate_work_diagnostics$method == "affine" &
        aggregate_work_diagnostics$scale == 1
    ])
  ),
  paste0("bank exact ODE calls: ", result$bank_exact_ode_calls),
  paste0("replay proposed exact ODE calls: ",
         result$replay_proposed_exact_ode_calls),
  paste(
    "Interpretation: this pass/fail concerns only frozen-map endpoint",
    "mechanics at an oracle anchor. It is not evidence for a SAEM handoff,",
    "posterior correctness of a retained chain, or full-cohort scaling."
  ),
  paste(
    "Current likelihood values are treated as cached during replay; each",
    "endpoint is charged exactly one VODE-BDF solve per patient."
  ),
  paste(
    "Acceptance ratios assume the audited symmetric additive dynamic-psi",
    "proposal. Zeta uses the full-115 comparator scale, although this replay",
    "contains only the twelve audited patients."
  )
)
writeLines(report, file.path(output_dir, "REPORT.txt"))
writeLines("complete", file.path(output_dir, "COMPLETED"))
cat(paste(report, collapse = "\n"), "\n")
