# Experimental numerical target: one VODE solve at the union of each patient's
# ALREADY adjusted positive times. This never replaces the sealed callbacks.
# A union grid can alter adaptive solver paths and numerical failure boundaries.

sab_make_shared_trajectory_probe <- function(adapter) {
  sab_validate_system_a_adapter(adapter)
  owner <- environment(adapter$solve_prediction)
  upstream <- get("upstream", envir = owner, inherits = FALSE)
  exact <- get("exact_likelihood", envir = owner, inherits = FALSE)
  patients <- exact$patients[adapter$patient_ids]
  controls <- patients[[1L]]$controls
  if (!all(vapply(patients, function(p) identical(p$controls, controls),
                  logical(1L)))) stop("Shared trajectory requires common controls.")
  positive_times <- lapply(patients, function(p) {
    upstream$sysa_stan_positive_times(p$obs_time, p$controls$time_eps)
  })
  union_times <- sort(unique(as.numeric(unlist(positive_times))))
  target_name <- "experimental_union_grid_vode_bdf_v1"
  fingerprint <- upstream$sysam_sha256_object(list(
    numerical_target = target_name,
    reference_target = adapter$target_fingerprint,
    reference_likelihood = adapter$likelihood_signature,
    patients = adapter$patient_ids, controls = controls,
    positive_times = positive_times, union_times = union_times,
    solve_body = body(.sab_shared_solve),
    prediction_body = body(.sab_shared_prediction)
  ))
  structure(list(
    numerical_target = target_name, target_fingerprint = fingerprint,
    reference_target_fingerprint = adapter$target_fingerprint,
    adapter = adapter, upstream = upstream, patients = patients,
    controls = controls, positive_times = positive_times,
    union_times = union_times
  ), class = "sab_shared_trajectory_probe")
}

.sab_shared_solve <- function(probe, x, psi, ode_function = deSolve::ode) {
  upstream <- probe$upstream
  psi <- upstream$sysa_named_finite_vector(
    psi, upstream$sysa_psi_names(), "psi"
  )
  subject <- upstream$sysa_local_to_natural(x)
  shared <- upstream$sysa_psi_to_natural(psi)
  failed <- function(reason, equilibrium = NULL, calls = 0L, message = NULL) {
    list(ok = FALSE, reason = reason, equilibrium = equilibrium,
         ode_integrations = calls, message = message)
  }
  if (any(!is.finite(subject)) || any(!is.finite(shared))) {
    return(failed("nonfinite_natural_parameter"))
  }
  y0 <- upstream$sysa_stan_equilibrium(subject, shared, probe$controls$mu_v)
  if (is.null(y0)) return(failed("invalid_equilibrium"))
  solution <- NULL
  calls <- 0L
  if (length(probe$union_times)) {
    parameters <- c(subject, gamma = shared[["gamma"]],
                    mu_l = shared[["mu_l"]], mu_v = probe$controls$mu_v)
    requested_times <- c(0, probe$union_times)
    expected_failure <- NULL
    calls <- 1L
    solution <- withCallingHandlers(tryCatch(
      ode_function(
        y = y0, times = requested_times, func = upstream$sysa_stan_rhs,
        parms = parameters, method = "vode", mf = 22,
        rtol = probe$controls$rel_tol, atol = probe$controls$abs_tol,
        maxsteps = probe$controls$max_num_steps
      ), error = function(condition) {
        if (!upstream$sysa_is_expected_ode_condition(condition)) stop(condition)
        expected_failure <<- conditionMessage(condition)
        NULL
      }), warning = function(condition) {
        if (!upstream$sysa_is_expected_ode_condition(condition)) {
          stop("Unexpected warning from union-grid solve: ",
               conditionMessage(condition), call. = FALSE)
        }
        expected_failure <<- conditionMessage(condition)
        invokeRestart("muffleWarning")
      })
    if (!is.null(expected_failure)) {
      return(failed("ode_failure", y0, calls, expected_failure))
    }
    checked <- upstream$sysa_validate_solver_output(solution, requested_times)
    if (!isTRUE(checked$ok)) return(failed(checked$reason, y0, calls))
    solution <- checked$solution
  }
  list(ok = TRUE, equilibrium = y0, solution = solution,
       ode_integrations = calls)
}

.sab_shared_prediction <- function(probe, solved, patient_id, psi) {
  patient <- probe$patients[[patient_id]]
  positive <- probe$positive_times[[patient_id]]
  # An integration failure must not invalidate patients observed only at t<=0.
  if (!isTRUE(solved$ok) &&
      (is.null(solved$equilibrium) || length(positive))) {
    return(list(ok = FALSE, reason = solved$reason))
  }
  rows <- match(positive, probe$union_times)
  if (anyNA(rows)) stop("Patient effective time is missing from the union grid.")
  y0 <- solved$equilibrium
  means <- numeric(length(patient$y))
  index <- 0L
  for (j in seq_along(means)) {
    state <- y0
    if (patient$obs_time[[j]] > 0) {
      index <- index + 1L
      state <- solved$solution[rows[[index]] + 1L, 2:6]
    }
    means[[j]] <- if (patient$ytype[[j]] == 1L) {
      log10(max(1000 * (state[[4L]] + state[[5L]]), 1e-30))
    } else max(state[[1L]] + state[[2L]] + state[[3L]], 1e-12)
  }
  list(
    ok = TRUE, patient_id = patient_id, observation_mean = means,
    y = patient$y, cens = patient$cens, ytype = patient$ytype,
    adjusted_positive_times = positive, equilibrium = y0,
    dynamic_psi = psi[c("log_gamma_pop", "log_mu_l_pop")],
    solver = probe$numerical_target
  )
}

# Counts original integrations from the sealed solver's exhaustive outcomes:
# natural-parameter/equilibrium rejections occur before its sole integration.
.sab_shared_original_calls <- function(prediction) {
  if (isTRUE(prediction$ok)) {
    return(as.integer(length(prediction$adjusted_positive_times) > 0L))
  }
  if (prediction$reason %in% c("nonfinite_natural_parameter",
                              "invalid_equilibrium")) return(0L)
  if (prediction$reason %in% c("ode_failure", "invalid_ode_output",
                              "invalid_ode_times")) return(1L)
  stop("Unexpected original prediction failure: ", prediction$reason)
}

# Both paths include reconstruction and observation likelihood in elapsed time.
# The common path always integrates the full frozen union, even for a subset of
# observations; callers must not vary the grid based on x, psi or solver success.
sab_shared_trajectory_evaluate <- function(
    probe, x, psi, compare_original = TRUE, common_first = FALSE) {
  if (!inherits(probe, "sab_shared_trajectory_probe")) {
    stop("probe must come from sab_make_shared_trajectory_probe().")
  }
  ids <- names(probe$patients)
  evaluate <- function(common) {
    started <- proc.time()[["elapsed"]]
    solved <- if (common) .sab_shared_solve(probe, x, psi) else NULL
    predictions <- setNames(lapply(ids, function(id) {
      if (common) .sab_shared_prediction(probe, solved, id, psi) else {
        probe$adapter$solve_prediction(id, x, psi)
      }
    }), ids)
    loglik <- vapply(predictions, function(p) {
      probe$adapter$loglik_from_prediction(p, psi)
    }, numeric(1L))
    calls <- if (common) solved$ode_integrations else {
      sum(vapply(predictions, .sab_shared_original_calls, integer(1L)))
    }
    ode_failed <- function(p) {
      !isTRUE(p$ok) && p$reason %in% c(
        "ode_failure", "invalid_ode_output", "invalid_ode_times"
      )
    }
    list(predictions = predictions, loglik = loglik, ledger = data.frame(
      method = if (common) "common" else "original",
      prediction_calls = if (common) 1L else length(ids),
      observation_likelihood_calls = length(ids), ode_integrations = calls,
      ode_failures = if (common) as.integer(ode_failed(solved)) else {
        sum(vapply(predictions, ode_failed, logical(1L)))
      },
      prediction_failures = sum(!vapply(predictions, function(p) isTRUE(p$ok),
                                       logical(1L))),
      elapsed_seconds = proc.time()[["elapsed"]] - started
    ))
  }
  original <- common <- NULL
  if (isTRUE(common_first)) common <- evaluate(TRUE)
  if (isTRUE(compare_original)) original <- evaluate(FALSE)
  if (is.null(common)) common <- evaluate(TRUE)
  comparison <- NULL
  if (!is.null(original)) {
    both_finite <- is.finite(original$loglik) & is.finite(common$loglik)
    difference <- common$loglik - original$loglik
    difference[!both_finite] <- NA_real_
    comparison <- data.frame(
      patient_id = ids, original_loglik = unname(original$loglik),
      common_loglik = unname(common$loglik),
      loglik_difference = unname(difference),
      finite_support_agrees = is.finite(original$loglik) == is.finite(common$loglik),
      prediction_success_agrees =
        vapply(original$predictions, function(p) isTRUE(p$ok), logical(1L)) ==
        vapply(common$predictions, function(p) isTRUE(p$ok), logical(1L)),
      row.names = NULL
    )
  }
  list(
    numerical_target = probe$numerical_target,
    target_fingerprint = probe$target_fingerprint,
    reference_target_fingerprint = probe$reference_target_fingerprint,
    common_predictions = common$predictions,
    original_predictions = original$predictions,
    common_loglik = common$loglik, original_loglik = original$loglik,
    comparison = comparison,
    ledger = rbind(original$ledger, common$ledger)
  )
}
