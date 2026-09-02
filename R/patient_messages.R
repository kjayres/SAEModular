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
# These functions operate only on supplied log densities or log weights.  They
# never evaluate an ODE and make no transport-map approximation.

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
#' @param log_weights Log ratios `log q_new(x_s) - log q_anchor(x_s)` for
#'   ordered draws from the normalized anchor conditional.
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
#' function.  Density evaluations must be finite, which deliberately fails
#' closed when the supplied banks have unresolved support mismatch.
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
      stop(
        "Bridge iteration produced a non-finite estimate; check bank overlap.",
        call. = FALSE
      )
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

  stop(
    "Bridge iteration did not converge within `max_iterations`.",
    call. = FALSE
  )
}
