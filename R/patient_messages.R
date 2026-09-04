# Stable estimators for SAEM-anchored patient messages.
#
# If X_s is sampled from the normalized anchor conditional
#
#   pi_0(x) = q_0(x) / Z_0,
#
# and `log_weights[s] = log q_1(X_s) - log q_0(X_s)`, then
#
#   log(Z_1 / Z_0) = log E_pi0[exp(log_weights)].
#
# For population messages, q_1 / q_0 reduces to
#
#   g_(i,M)(x | eta, z_i) / g_(i,A)(x | eta_A, z_i)
#
# when the observation model and shared parameters are fixed.  The population
# densities may have any evaluable normalized form and may depend on model
# structure and known patient covariates.  A product of patient messages assumes
# that these densities factor by patient conditional on the global state.
# These functions operate only on supplied log densities or log weights.  They
# never evaluate an ODE and make no Gaussian, diagonal, or transport-map
# approximation.

.sab_pm_validate_log_vector <- function(value,
                                        argument,
                                        minimum_length = 1L,
                                        allow_negative_infinity = TRUE,
                                        require_positive_weight = FALSE) {
  if (!is.numeric(value) || is.object(value) || !is.null(dim(value)) ||
      length(value) < minimum_length) {
    stop(
      "`", argument, "` must be a numeric vector of length at least ",
      minimum_length, ".",
      call. = FALSE
    )
  }
  if (anyNA(value) || any(is.nan(value)) || any(value == Inf)) {
    stop(
      "`", argument, "` must not contain NA, NaN, or positive infinity.",
      call. = FALSE
    )
  }
  if (!allow_negative_infinity && any(value == -Inf)) {
    stop("`", argument, "` must contain only finite values.", call. = FALSE)
  }
  if (require_positive_weight && !any(is.finite(value))) {
    stop(
      "`", argument, "` must contain at least one finite log weight.",
      call. = FALSE
    )
  }
  as.numeric(value)
}

.sab_pm_validate_positive_scalar <- function(value, argument) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value <= 0) {
    stop("`", argument, "` must be a finite positive scalar.", call. = FALSE)
  }
  as.numeric(value)
}

.sab_pm_validate_count <- function(value, argument, minimum = 1L) {
  if (!is.numeric(value) || length(value) != 1L || !is.finite(value) ||
      value != floor(value) || value < minimum ||
      value > .Machine$integer.max) {
    stop(
      "`", argument, "` must be an integer at least ", minimum, ".",
      call. = FALSE
    )
  }
  as.integer(value)
}

#' Classical split-Rhat and multi-chain autocorrelation ESS
#'
#' These diagnose Markov-chain persistence and between-chain disagreement.
#' They are deliberately separate from importance-weight ESS, which measures
#' weight concentration and ignores draw order.  Rhat here is the classical
#' split statistic (not rank-normalized); tail behaviour must be examined with
#' the weight and split-message diagnostics as well.
#'
#' @param chains Named or unnamed list of at least two equal-length finite
#'   numeric chain vectors.  Each chain must contain at least eight draws.
#'
#' @return Split-Rhat, Geyer initial-positive-sequence MCMC ESS, chain count,
#'   split length, and the range of unsplit chain means.
#' @export
sab_mcmc_scalar_diagnostics <- function(chains) {
  if (!is.list(chains) || length(chains) < 2L ||
      any(!vapply(chains, is.numeric, logical(1L)))) {
    stop("`chains` must be a list of at least two numeric vectors.",
         call. = FALSE)
  }
  lengths <- vapply(chains, length, integer(1L))
  if (length(unique(lengths)) != 1L || lengths[[1L]] < 8L ||
      any(!vapply(chains, function(value) all(is.finite(value)), logical(1L)))) {
    stop("MCMC chains must be finite, equal-length, and contain >= 8 draws.",
         call. = FALSE)
  }
  half <- lengths[[1L]] %/% 2L
  split_chains <- unlist(lapply(chains, function(value) {
    list(value[seq_len(half)], tail(value, half))
  }), recursive = FALSE)
  draws <- do.call(cbind, split_chains)
  n <- nrow(draws)
  m <- ncol(draws)
  chain_means <- colMeans(draws)
  chain_variances <- apply(draws, 2L, stats::var)
  within <- mean(chain_variances)
  between <- n * stats::var(chain_means)
  variance_plus <- (n - 1) / n * within + between / n
  all_constant <- all(draws == draws[[1L]])
  if (within == 0) {
    split_rhat <- if (all_constant) 1 else Inf
    mcmc_ess <- if (all_constant) n * m else 1
  } else {
    # Sampling variation can make the classical finite-sample estimate dip
    # just below one.  Report the conventional lower-bounded value.
    split_rhat <- max(1, sqrt(variance_plus / within))
    rho <- numeric(max(0L, n - 1L))
    if (length(rho)) {
      for (lag in seq_along(rho)) {
        autocovariance <- mean(vapply(seq_len(m), function(chain) {
          centered <- draws[, chain] - chain_means[[chain]]
          sum(
            centered[seq_len(n - lag)] *
              centered[seq.int(lag + 1L, n)]
          ) / n
        }, numeric(1L)))
        rho[[lag]] <- 1 - (within - autocovariance) / variance_plus
      }
    }
    paired_sum <- numeric()
    if (length(rho) >= 2L) {
      for (start in seq.int(1L, length(rho) - 1L, by = 2L)) {
        pair <- rho[[start]] + rho[[start + 1L]]
        if (!is.finite(pair) || pair < 0) break
        paired_sum <- c(paired_sum, pair)
      }
    }
    integrated_time <- max(1, 1 + 2 * sum(paired_sum))
    mcmc_ess <- min(n * m, max(1, n * m / integrated_time))
  }
  structure(list(
    split_rhat = as.numeric(split_rhat),
    mcmc_ess = as.numeric(mcmc_ess),
    relative_mcmc_ess = as.numeric(mcmc_ess / (n * m)),
    original_chains = length(chains),
    split_chains = m,
    draws_per_split = n,
    chain_mean_range = diff(range(vapply(chains, mean, numeric(1L))))
  ), class = c("sab_mcmc_scalar_diagnostics", "list"))
}

#' Form distribution-free patient population log weights
#'
#' This is the minimal model-facing numerical interface for a
#' fixed-shared-parameter raw message.  The caller evaluates the normalized
#' anchor and candidate
#' patient population densities on the same stored states.  Those evaluators
#' may depend on patient covariates and model structure and need not be
#' Gaussian.  Their normalizing terms must not be omitted.  This numerical
#' helper validates the supplied values only; absolute continuity and tail
#' overlap away from the observed bank are obligations of the caller.
#'
#' @param log_g_candidate Values of
#'   `log g_(i,M)(x_s | eta, z_i)` on anchor-bank draws.  Negative infinity is
#'   allowed where the candidate density is zero.
#' @param log_g_anchor Values of
#'   `log g_(i,A)(x_s | eta_A, z_i)` on the same draws and in the same order.
#'   These must be finite: an anchor-conditional draw cannot lie outside the
#'   support of its own population density.
#'
#' @return The ordered log importance weights for
#'   [sab_raw_message_diagnostics()].
#' @export
sab_population_log_weights <- function(log_g_candidate, log_g_anchor) {
  log_g_candidate <- .sab_pm_validate_log_vector(
    log_g_candidate, "log_g_candidate"
  )
  log_g_anchor <- .sab_pm_validate_log_vector(
    log_g_anchor, "log_g_anchor", allow_negative_infinity = FALSE
  )
  if (length(log_g_candidate) != length(log_g_anchor)) {
    stop(
      "Candidate and anchor log population densities must have equal length.",
      call. = FALSE
    )
  }
  log_weights <- log_g_candidate - log_g_anchor
  .sab_pm_validate_log_vector(log_weights, "population log weights")
}

#' Evaluate a population-message density ratio on stored patient states
#'
#' This helper deliberately knows nothing about the population family.  The
#' supplied callback may implement a correlated, non-Gaussian, mixture, or
#' covariate-dependent normalized density.  It is called with the same patient
#' context for the anchor and candidate models, so observed covariates cannot
#' be silently dropped between the two evaluations.
#'
#' @param draws Finite numeric matrix.  Rows are stored patient states and
#'   columns are the declared local coordinates.
#' @param anchor_parameter,candidate_parameter Arbitrary objects understood by
#'   `log_population_density`; these may include model structure as well as
#'   continuous population parameters.
#' @param patient_context Fixed patient metadata, including any covariates,
#'   passed unchanged to both density evaluations.
#' @param log_population_density Function of `(x, parameter, patient_context)`
#'   returning one normalized log density with respect to the declared local
#'   state measure.  Negative infinity is allowed for zero candidate density.
#'
#' @return A list containing anchor and candidate log densities and their
#'   ordered log ratio.
#' @export
sab_evaluate_population_message <- function(
    draws, anchor_parameter, candidate_parameter, patient_context,
    log_population_density) {
  if (!is.matrix(draws) || !is.numeric(draws) || nrow(draws) < 1L ||
      ncol(draws) < 1L || any(!is.finite(draws)) ||
      is.null(colnames(draws)) || anyNA(colnames(draws)) ||
      any(!nzchar(colnames(draws))) || anyDuplicated(colnames(draws))) {
    stop(
      "`draws` must be a finite numeric matrix with unique coordinate names.",
      call. = FALSE
    )
  }
  if (!is.function(log_population_density)) {
    stop("`log_population_density` must be a function.", call. = FALSE)
  }
  evaluate <- function(parameter, label) {
    values <- vapply(seq_len(nrow(draws)), function(index) {
      value <- log_population_density(
        draws[index, ], parameter, patient_context
      )
      if (!is.numeric(value) || length(value) != 1L || is.na(value) ||
          is.nan(value) || value == Inf) {
        stop(
          "The population-density callback returned an invalid value for ",
          label, " at draw ", index, ".", call. = FALSE
        )
      }
      as.numeric(value)
    }, numeric(1L))
    .sab_pm_validate_log_vector(values, label)
  }

  log_g_anchor <- evaluate(anchor_parameter, "anchor log population density")
  log_g_candidate <- evaluate(
    candidate_parameter, "candidate log population density"
  )
  log_weights <- sab_population_log_weights(
    log_g_candidate = log_g_candidate,
    log_g_anchor = log_g_anchor
  )
  structure(
    list(
      log_g_anchor = log_g_anchor,
      log_g_candidate = log_g_candidate,
      log_weights = log_weights
    ),
    class = c("sab_population_message_evaluation", "list")
  )
}

.sab_pm_log_sum_exp <- function(value) {
  largest <- max(value)
  if (largest == -Inf) {
    return(-Inf)
  }
  largest + log(sum(exp(value - largest)))
}

.sab_pm_log_add_exp <- function(first, second) {
  largest <- pmax(first, second)
  largest + log(exp(first - largest) + exp(second - largest))
}

#' Stable logarithm of a mean of exponentials
#'
#' @param log_values Numeric log values.  Negative infinity is allowed, but
#'   positive infinity and missing values are rejected.
#'
#' @return `log(mean(exp(log_values)))`, evaluated without exponentiating the
#'   largest input.
#' @export
sab_log_mean_exp <- function(log_values) {
  log_values <- .sab_pm_validate_log_vector(log_values, "log_values")
  .sab_pm_log_sum_exp(log_values) - log(length(log_values))
}

.sab_pm_weight_summary <- function(log_weights, argument = "log_weights") {
  log_weights <- .sab_pm_validate_log_vector(
    log_weights,
    argument,
    require_positive_weight = TRUE
  )
  n_draws <- length(log_weights)
  log_sum_weights <- .sab_pm_log_sum_exp(log_weights)
  log_normalized_weights <- log_weights - log_sum_weights
  log_sum_squared_weights <- .sab_pm_log_sum_exp(
    2 * log_normalized_weights
  )

  weight_ess <- exp(-log_sum_squared_weights)
  # Protect exact identities from harmless floating-point excursions.
  weight_ess <- min(as.numeric(n_draws), max(1, weight_ess))
  d2 <- log(n_draws) + log_sum_squared_weights
  d2 <- max(0, d2)

  list(
    n_draws = n_draws,
    log_ratio = log_sum_weights - log(n_draws),
    weight_ess = weight_ess,
    relative_weight_ess = weight_ess / n_draws,
    max_normalized_weight = exp(max(log_normalized_weights)),
    d2 = d2
  )
}

#' Contiguous batch-means MCSE for a log message estimate
#'
#' The batches are contiguous and differ in size by at most one.  If `B` is the
#' number of batches, `n_b` is a batch size, `M_b` is its mean raw weight, and
#' `M` is the overall mean raw weight, the returned log-scale MCSE is
#'
#'   sqrt(sum_b n_b (M_b / M - 1)^2 / ((B - 1) n)).
#'
#' This reduces to the usual equal-size non-overlapping batch-means estimator,
#' followed by the delta method for `log(M)`.  The calculation uses log batch
#' means, so a common extreme scale in the weights does not cause overflow.  A
#' caller must choose batches long relative to the chain's autocorrelation
#' range and should check stability under coarser batching.
#'
#' @param log_weights Log importance weights in Markov-chain order.
#' @param n_batches Number of contiguous batches.  Every batch must contain at
#'   least two draws.
#'
#' @return A list containing the batch boundaries, batch log estimates, and
#'   the log-scale MCSE.
#' @export
sab_batch_log_mcse <- function(log_weights, n_batches) {
  log_weights <- .sab_pm_validate_log_vector(
    log_weights,
    "log_weights",
    minimum_length = 4L,
    require_positive_weight = TRUE
  )
  n_batches <- .sab_pm_validate_count(n_batches, "n_batches", minimum = 2L)
  n_draws <- length(log_weights)
  if (2L * n_batches > n_draws) {
    stop("Every contiguous batch must contain at least two draws.", call. = FALSE)
  }

  base_size <- n_draws %/% n_batches
  remainder <- n_draws %% n_batches
  batch_sizes <- rep.int(base_size, n_batches)
  if (remainder > 0L) {
    batch_sizes[seq_len(remainder)] <- batch_sizes[seq_len(remainder)] + 1L
  }
  batch_ends <- cumsum(batch_sizes)
  batch_starts <- c(1L, head(batch_ends, -1L) + 1L)
  batch_log_ratios <- vapply(
    seq_len(n_batches),
    function(batch) {
      sab_log_mean_exp(log_weights[batch_starts[batch]:batch_ends[batch]])
    },
    numeric(1)
  )

  overall_log_ratio <- sab_log_mean_exp(log_weights)
  relative_batch_means <- exp(batch_log_ratios - overall_log_ratio)
  log_scale_variance <- sum(
    batch_sizes * (relative_batch_means - 1)^2
  ) / ((n_batches - 1) * n_draws)

  structure(
    list(
      n_batches = n_batches,
      batch_sizes = batch_sizes,
      batch_starts = batch_starts,
      batch_ends = batch_ends,
      batch_log_ratios = batch_log_ratios,
      relative_batch_means = relative_batch_means,
      log_scale_mcse = sqrt(log_scale_variance)
    ),
    class = "sab_batch_log_mcse"
  )
}

#' Split-bank stability diagnostics for a raw patient message
#'
#' @param log_weights Log importance weights in Markov-chain order.  At least
#'   four values are required.  Both contiguous halves must contain a positive
#'   empirical weight.
#'
#' @return First- and second-half log ratios, weight diagnostics, and absolute
#'   differences in the log ratio and empirical D2 divergence.
#' @export
sab_split_message_diagnostics <- function(log_weights) {
  log_weights <- .sab_pm_validate_log_vector(
    log_weights,
    "log_weights",
    minimum_length = 4L,
    require_positive_weight = TRUE
  )
  first_end <- length(log_weights) %/% 2L
  first <- .sab_pm_weight_summary(
    log_weights[seq_len(first_end)],
    argument = "first-half log_weights"
  )
  second <- .sab_pm_weight_summary(
    log_weights[seq.int(first_end + 1L, length(log_weights))],
    argument = "second-half log_weights"
  )

  structure(
    list(
      first = first,
      second = second,
      signed_log_ratio_difference = first$log_ratio - second$log_ratio,
      absolute_log_ratio_difference = abs(first$log_ratio - second$log_ratio),
      absolute_d2_difference = abs(first$d2 - second$d2)
    ),
    class = "sab_split_message_diagnostics"
  )
}

#' Estimate a raw SAEM-anchored patient message and its reliability
#'
#' @param log_weights Ordered log population-density ratios returned by
#'   [sab_population_log_weights()], or more general ratios
#'   `log q_new(x_s) - log q_anchor(x_s)`, on draws from the normalized anchor
#'   conditional.
#' @param n_batches Number of contiguous batches used for the batch-means
#'   log-scale MCSE.  The default uses at most ten batches and leaves at least
#'   two draws in each batch.
#'
#' @return A list with the log message ratio, raw-weight ESS, maximum normalized
#'   weight, empirical order-2 Renyi divergence (`d2`), batch MCSE, and split
#'   diagnostics.  Here `d2 = log(n / ESS)` exactly for the empirical weights.
#'   Weight ESS measures importance-weight concentration only; it does not
#'   account for MCMC autocorrelation.
#' @export
sab_raw_message_diagnostics <- function(log_weights, n_batches = NULL) {
  log_weights <- .sab_pm_validate_log_vector(
    log_weights,
    "log_weights",
    minimum_length = 4L,
    require_positive_weight = TRUE
  )
  if (is.null(n_batches)) {
    n_batches <- min(10L, length(log_weights) %/% 2L)
  }

  summary <- .sab_pm_weight_summary(log_weights)
  summary$batch <- sab_batch_log_mcse(log_weights, n_batches)
  summary$split <- sab_split_message_diagnostics(log_weights)
  structure(summary, class = "sab_raw_message_diagnostics")
}

#' Bidirectional optimal bridge estimate of a normalizing-constant ratio
#'
#' Let `q0` and `q1` be unnormalized densities with constants `Z0` and `Z1`.
#' This function receives evaluations of both densities on independent draws
#' from `q0 / Z0` and `q1 / Z1`.  It solves the optimal bridge fixed-point
#' equation for `log(Z1 / Z0)`, using sample-size fractions in the bridge
#' function.  Neither density has to be Gaussian.  Density evaluations must be
#' finite, which deliberately fails closed when the supplied banks have
#' unresolved support mismatch.
#'
#' @param log_q0_on_0,log_q1_on_0 Log densities on draws from `q0 / Z0`.
#' @param log_q0_on_1,log_q1_on_1 Log densities on draws from `q1 / Z1`.
#' @param tolerance Absolute convergence tolerance on the log ratio.
#' @param max_iterations Maximum fixed-point iterations.
#' @param initial_log_ratio Optional finite starting value.  By default, the
#'   midpoint of the forward and reverse raw-importance estimates is used.
#'
#' @return A list containing `log_ratio`, convergence history, iteration count,
#'   and the two directional raw-importance estimates.
#' @export
sab_bridge_log_ratio <- function(log_q0_on_0,
                                 log_q1_on_0,
                                 log_q0_on_1,
                                 log_q1_on_1,
                                 tolerance = 1e-10,
                                 max_iterations = 10000L,
                                 initial_log_ratio = NULL) {
  log_q0_on_0 <- .sab_pm_validate_log_vector(
    log_q0_on_0, "log_q0_on_0", minimum_length = 2L,
    allow_negative_infinity = FALSE
  )
  log_q1_on_0 <- .sab_pm_validate_log_vector(
    log_q1_on_0, "log_q1_on_0", minimum_length = 2L,
    allow_negative_infinity = FALSE
  )
  log_q0_on_1 <- .sab_pm_validate_log_vector(
    log_q0_on_1, "log_q0_on_1", minimum_length = 2L,
    allow_negative_infinity = FALSE
  )
  log_q1_on_1 <- .sab_pm_validate_log_vector(
    log_q1_on_1, "log_q1_on_1", minimum_length = 2L,
    allow_negative_infinity = FALSE
  )
  if (length(log_q0_on_0) != length(log_q1_on_0)) {
    stop("The two log-density vectors on q0 draws must have equal length.",
         call. = FALSE)
  }
  if (length(log_q0_on_1) != length(log_q1_on_1)) {
    stop("The two log-density vectors on q1 draws must have equal length.",
         call. = FALSE)
  }
  tolerance <- .sab_pm_validate_positive_scalar(tolerance, "tolerance")
  max_iterations <- .sab_pm_validate_count(
    max_iterations, "max_iterations", minimum = 1L
  )

  log_density_ratio_on_0 <- log_q1_on_0 - log_q0_on_0
  log_density_ratio_on_1 <- log_q1_on_1 - log_q0_on_1
  if (any(!is.finite(log_density_ratio_on_0)) ||
      any(!is.finite(log_density_ratio_on_1))) {
    stop("Log-density ratios must be finite.", call. = FALSE)
  }

  forward_log_ratio <- sab_log_mean_exp(log_density_ratio_on_0)
  reverse_log_ratio <- -sab_log_mean_exp(-log_density_ratio_on_1)
  if (is.null(initial_log_ratio)) {
    log_ratio <- (forward_log_ratio + reverse_log_ratio) / 2
  } else {
    if (!is.numeric(initial_log_ratio) || length(initial_log_ratio) != 1L ||
        !is.finite(initial_log_ratio)) {
      stop("`initial_log_ratio` must be a finite numeric scalar.", call. = FALSE)
    }
    log_ratio <- as.numeric(initial_log_ratio)
  }

  n0 <- length(log_density_ratio_on_0)
  n1 <- length(log_density_ratio_on_1)
  log_s0 <- log(n0) - log(n0 + n1)
  log_s1 <- log(n1) - log(n0 + n1)
  history <- c(log_ratio)

  for (iteration in seq_len(max_iterations)) {
    log_denominator_on_0 <- .sab_pm_log_add_exp(
      log_s0 + log_ratio,
      log_s1 + log_density_ratio_on_0
    )
    log_denominator_on_1 <- .sab_pm_log_add_exp(
      log_s0 + log_ratio,
      log_s1 + log_density_ratio_on_1
    )
    updated_log_ratio <-
      sab_log_mean_exp(
        log_density_ratio_on_0 - log_denominator_on_0
      ) -
      sab_log_mean_exp(-log_denominator_on_1)

    if (!is.finite(updated_log_ratio)) {
      stop(.sab_pm_bridge_nonconvergence(
        message = paste0(
          "Bridge iteration produced a non-finite estimate; ",
          "check bank overlap."
        ),
        failure_code = "nonfinite_iteration",
        log_ratio = log_ratio,
        iterations = iteration,
        final_increment = NA_real_,
        history = history,
        n0 = n0,
        n1 = n1,
        forward_log_ratio = forward_log_ratio,
        reverse_log_ratio = reverse_log_ratio
      ))
    }
    history <- c(history, updated_log_ratio)
    increment <- abs(updated_log_ratio - log_ratio)
    log_ratio <- updated_log_ratio
    if (increment <= tolerance) {
      return(structure(
        list(
          log_ratio = log_ratio,
          iterations = iteration,
          converged = TRUE,
          final_increment = increment,
          history = history,
          n0 = n0,
          n1 = n1,
          forward_log_ratio = forward_log_ratio,
          reverse_log_ratio = reverse_log_ratio
        ),
        class = "sab_bridge_log_ratio"
      ))
    }
  }

  stop(.sab_pm_bridge_nonconvergence(
    message = "Bridge iteration did not converge within `max_iterations`.",
    failure_code = "max_iterations",
    log_ratio = log_ratio,
    iterations = max_iterations,
    final_increment = increment,
    history = history,
    n0 = n0,
    n1 = n1,
    forward_log_ratio = forward_log_ratio,
    reverse_log_ratio = reverse_log_ratio
  ))
}

.sab_pm_bridge_nonconvergence <- function(message, failure_code, log_ratio,
                                          iterations, final_increment, history,
                                          n0, n1, forward_log_ratio,
                                          reverse_log_ratio) {
  estimate <- structure(list(
    log_ratio = as.numeric(log_ratio),
    iterations = as.integer(iterations),
    converged = FALSE,
    final_increment = as.numeric(final_increment),
    history = as.numeric(history),
    n0 = as.integer(n0),
    n1 = as.integer(n1),
    forward_log_ratio = as.numeric(forward_log_ratio),
    reverse_log_ratio = as.numeric(reverse_log_ratio)
  ), class = "sab_bridge_log_ratio")
  structure(
    list(
      message = message,
      call = NULL,
      failure_code = failure_code,
      estimate = estimate
    ),
    class = c("sab_bridge_nonconvergence", "error", "condition")
  )
}
