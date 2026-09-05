# ODE-free System A population updates. The location conditional matches
# modular_bayes/research/system_a_modular_v0/R/stage2_exact_workflow.R:800.
# All scientific priors are supplied by the pinned target's caller.

.sab_population_names <- function() list(
  local = c("u_log_lambda", "u_log_mu_t", "u_log_mu_a", "u_log_p",
            "u_log_alpha_l", "u_pi", "u_eta_rti", "u_eta_pi"),
  eta = c("mu_log_lambda", "mu_log_mu_t", "mu_log_mu_a", "mu_log_p",
          "mu_log_alpha_l", "mu_u_pi", "mu_u_eta_rti", "mu_u_eta_pi",
          "beta_nelf", "log_omega_lambda", "log_omega_mu_t", "log_omega_mu_a",
          "log_omega_p", "log_omega_alpha_l", "log_omega_pi",
          "log_omega_eta_rti", "log_omega_eta_pi")
)

.sab_validate_population_inputs <- function(x, eta, treatment, prior_mean, prior_sd) {
  expected <- .sab_population_names()
  if (!is.matrix(x) || !is.numeric(x) || nrow(x) < 1L || ncol(x) != 8L ||
      any(!is.finite(x)) || !identical(colnames(x), expected$local) ||
      is.null(rownames(x)) || anyNA(rownames(x)) ||
      any(!nzchar(rownames(x))) || anyDuplicated(rownames(x))) {
    stop("x requires finite values, canonical columns and unique patient row names.")
  }
  for (value in list(eta, prior_mean, prior_sd)) {
    if (!is.numeric(value) || !is.null(dim(value)) || length(value) != 17L ||
        !identical(names(value), expected$eta) || any(!is.finite(value))) {
      stop("eta and prior vectors require finite canonical named coordinates.")
    }
  }
  if (any(prior_sd <= 0)) stop("prior_sd must be positive.")
  if (!is.numeric(treatment) || !is.null(dim(treatment)) ||
      length(treatment) != nrow(x) || anyNA(treatment) ||
      any(!treatment %in% c(0, 1)) ||
      (!is.null(names(treatment)) && !identical(names(treatment), rownames(x)))) {
    stop("treatment must be aligned with patient rows and contain only 0 or 1.")
  }
  invisible(TRUE)
}

sab_system_a_location_conditional <- function(x, eta, treatment, prior_mean, prior_sd) {
  .sab_validate_population_inputs(x, eta, treatment, prior_mean, prior_sd)
  weight <- exp(-2 * eta[10:17])
  prior_precision <- 1 / prior_sd[1:9]^2
  if (any(!is.finite(weight)) || any(!is.finite(prior_precision)) ||
      any(prior_precision <= 0)) stop("Population precision is not representable.")
  precision <- prior_precision[1:7] + nrow(x) * weight[1:7]
  means <- (prior_precision[1:7] * prior_mean[1:7] +
              colSums(x[, 1:7, drop = FALSE]) * weight[1:7]) / precision
  design <- cbind(1, as.numeric(treatment))
  pair_precision <- diag(as.numeric(prior_precision[8:9])) +
    crossprod(design) * weight[[8L]]
  rhs <- prior_precision[8:9] * prior_mean[8:9] +
    as.numeric(crossprod(design, x[, 8L])) * weight[[8L]]
  pair_mean <- setNames(as.numeric(solve(pair_precision, rhs)), names(eta)[8:9])
  if (any(!is.finite(c(means, precision, pair_mean, pair_precision)))) {
    stop("Population location conditional is not representable.")
  }
  list(independent_mean = setNames(as.numeric(means), names(eta)[1:7]),
       independent_sd = setNames(1 / sqrt(as.numeric(precision)), names(eta)[1:7]),
       pair_mean = pair_mean, pair_precision = pair_precision,
       pair_chol_precision = chol(pair_precision))
}

# Conditional log density with respect to d(log omega). The Gaussian prior
# is already on log omega, so there is NO additional exp(s) Jacobian here.
sab_system_a_scale_logdensity <- function(s, n, sse, prior_mean, prior_sd) {
  values <- list(s, n, sse, prior_mean, prior_sd)
  if (any(vapply(values, function(v) !is.numeric(v) || length(v) != 1L ||
                 !is.finite(v), logical(1L))) || n < 1 || n != floor(n) ||
      sse < 0 || prior_sd <= 0) stop("Invalid scale conditional arguments.")
  quadratic <- if (sse == 0) 0 else {
    log_quadratic <- log(sse) - 2 * s - log(2)
    if (log_quadratic > log(.Machine$double.xmax)) Inf else exp(log_quadratic)
  }
  value <- -n * s - quadratic - .5 * ((s - prior_mean) / prior_sd)^2
  if (is.finite(value)) value else -Inf
}

# Unbounded step-out and shrinkage. A cap is a fatal diagnostic, never a silent
# rejected transition; completing a chain after such truncation would need a
# separate invariance argument. Production callers always use the defaults.
.sab_population_slice <- function(current, log_density, width = .5, max_evaluations = 1000L) {
  if (!is.numeric(current) || length(current) != 1L || !is.finite(current) ||
      !is.function(log_density) || length(width) != 1L || !is.finite(width) ||
      width <= 0 || length(max_evaluations) != 1L || !is.finite(max_evaluations) ||
      max_evaluations < 1 || max_evaluations != floor(max_evaluations)) {
    stop("Invalid scalar slice arguments.")
  }
  evaluations <- 0L
  evaluate <- function(value) {
    evaluations <<- evaluations + 1L
    if (evaluations > max_evaluations) stop("Slice evaluation cap exceeded.")
    result <- log_density(value)
    if (!is.numeric(result) || length(result) != 1L || is.na(result) || result == Inf) {
      stop("Slice log density returned an invalid value.")
    }
    result
  }
  initial <- evaluate(current)
  if (!is.finite(initial)) stop("Slice current state must have finite log density.")
  level <- initial - stats::rexp(1L)
  left <- current - stats::runif(1L) * width
  right <- left + width
  while (evaluate(left) > level) left <- left - width
  while (evaluate(right) > level) right <- right + width
  repeat {
    proposal <- stats::runif(1L, left, right)
    if (evaluate(proposal) >= level) return(proposal)
    if (proposal < current) left <- proposal else right <- proposal
    if (!(left < right)) stop("Slice bracket collapsed numerically.")
  }
}

sab_system_a_population_update <- function(x, eta, treatment, prior_mean, prior_sd) {
  conditional <- sab_system_a_location_conditional(x, eta, treatment, prior_mean, prior_sd)
  updated <- eta
  updated[1:7] <- stats::rnorm(7L, conditional$independent_mean, conditional$independent_sd)
  updated[8:9] <- conditional$pair_mean +
    backsolve(conditional$pair_chol_precision, stats::rnorm(2L))
  residual <- sweep(x, 2L, updated[1:8], "-")
  residual[, 8L] <- residual[, 8L] - treatment * updated[[9L]]
  sse <- colSums(residual^2)
  if (any(!is.finite(sse))) stop("Population residual sum of squares overflowed.")
  for (j in seq_len(8L)) {
    coordinate <- j + 9L
    updated[[coordinate]] <- .sab_population_slice(updated[[coordinate]], function(s) {
      sab_system_a_scale_logdensity(s, nrow(x), sse[[j]],
                                   prior_mean[[coordinate]], prior_sd[[coordinate]])
    })
  }
  updated
}
