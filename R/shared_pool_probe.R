# Offline feasibility calculations for a live pool of continuous patient states.
#
# For a fixed normalized reserve density q, an unused slot has density q and
# patient i's occupied slot has unnormalized density f_i. Holding other occupied
# slots fixed, a uniformly selected reserve is accepted with
#   min(1, (f_i(reserve)/q(reserve)) / (f_i(current)/q(current))).
# Gibbs probabilities over current + unused slots are proportional to f_i/q.
#
# This module evaluates those formulas on supplied states only. It does not run
# an assignment chain, estimate posterior ESS, or imply that different patients
# may simultaneously occupy the same slot. No Gaussian family is assumed.

.sab_pool_log_values <- function(x, argument, finite = FALSE) {
  if (!is.numeric(x) || !length(x) || anyNA(x) || any(x == Inf) ||
      (finite && any(!is.finite(x)))) {
    stop("`", argument, "` must contain ",
         if (finite) "finite log densities." else "finite values or -Inf.",
         call. = FALSE)
  }
  invisible(x)
}

#' Paired reserve proposal diagnostics for a fixed common q
#'
#' Rows are patients; columns are reserve draws sampled directly from q. Every
#' reserve, including solver failures with log f = -Inf, remains in the reported
#' denominators. Current states must have positive finite target/proposal density.
#' Patientwise unknown constants in f_i cancel, but q must be evaluated as the
#' same density at both current and reserve states.
#'
#' `squared_jump`, when supplied, is a nonnegative patient-by-reserve matrix in a
#' fixed, declared coordinate scaling. Its acceptance-weighted mean is a proposal
#' movement proxy, not posterior ESS or population ESJD.
#'
#' The usefulness threshold is descriptive and must be declared before viewing
#' results. Coverage of two patients means alternative uses of a trajectory; it
#' never means simultaneous occupancy. These summaries alone are not a go gate.
sab_shared_pool_probe <- function(log_f_reserve,
                                  log_q_reserve,
                                  log_f_current,
                                  log_q_current,
                                  useful_acceptance = 0.1,
                                  squared_jump = NULL) {
  if (!is.matrix(log_f_reserve) || any(dim(log_f_reserve) < 1L)) {
    stop("`log_f_reserve` must be a nonempty patient-by-reserve matrix.",
         call. = FALSE)
  }
  .sab_pool_log_values(log_f_reserve, "log_f_reserve")
  n_patient <- nrow(log_f_reserve)
  n_reserve <- ncol(log_f_reserve)
  vectors <- list(log_q_reserve = log_q_reserve,
                  log_f_current = log_f_current,
                  log_q_current = log_q_current)
  expected_lengths <- c(n_reserve, n_patient, n_patient)
  for (index in seq_along(vectors)) {
    value <- vectors[[index]]
    argument <- names(vectors)[[index]]
    if (!is.null(dim(value)) || length(value) != expected_lengths[[index]]) {
      stop("`", argument, "` has the wrong vector length.", call. = FALSE)
    }
    .sab_pool_log_values(value, argument, finite = TRUE)
  }
  if (!is.numeric(useful_acceptance) || length(useful_acceptance) != 1L ||
      !is.finite(useful_acceptance) || useful_acceptance <= 0 ||
      useful_acceptance > 1) {
    stop("`useful_acceptance` must be in (0, 1].", call. = FALSE)
  }
  if (!is.null(squared_jump) &&
      (!is.matrix(squared_jump) || !is.numeric(squared_jump) ||
       !identical(dim(squared_jump), dim(log_f_reserve)) ||
       any(!is.finite(squared_jump)) || any(squared_jump < 0))) {
    stop("`squared_jump` must be a finite nonnegative matrix matching log_f_reserve.",
         call. = FALSE)
  }

  log_w_reserve <- sweep(log_f_reserve, 2L, log_q_reserve, "-")
  log_w_current <- log_f_current - log_q_current
  .sab_pool_log_values(log_w_reserve, "reserve log f/q")
  .sab_pool_log_values(log_w_current, "current log f/q", finite = TRUE)
  log_ratio <- sweep(log_w_reserve, 1L, log_w_current, "-")
  acceptance <- exp(pmin(log_ratio, 0))
  gibbs <- t(vapply(seq_len(n_patient), function(i) {
    log_weights <- c(log_w_current[[i]], log_w_reserve[i, ])
    weights <- exp(log_weights - max(log_weights))
    weights / sum(weights)
  }, numeric(n_reserve + 1L)))

  patient_ids <- rownames(log_f_reserve)
  if (is.null(patient_ids)) patient_ids <- as.character(seq_len(n_patient))
  reserve_ids <- colnames(log_f_reserve)
  if (is.null(reserve_ids)) reserve_ids <- as.character(seq_len(n_reserve))
  dimnames(acceptance) <- dimnames(log_w_reserve) <-
    list(patient_ids, reserve_ids)
  dimnames(gibbs) <- list(patient_ids, c("current", reserve_ids))
  useful_counts <- colSums(acceptance >= useful_acceptance)
  movement <- if (is.null(squared_jump)) rep(NA_real_, n_patient) else
    rowMeans(acceptance * squared_jump)
  patient <- data.frame(
    patient_id = patient_ids,
    reserve_count = n_reserve,
    finite_target_fraction = rowMeans(is.finite(log_f_reserve)),
    mean_reserve_mh_acceptance = rowMeans(acceptance),
    useful_reserve_fraction = rowMeans(acceptance >= useful_acceptance),
    gibbs_move_probability = rowSums(gibbs[, -1L, drop = FALSE]),
    expected_accepted_squared_jump = movement,
    row.names = NULL
  )
  reserve <- data.frame(
    reserve_id = reserve_ids,
    finite_target_patients = colSums(is.finite(log_f_reserve)),
    useful_patients = useful_counts,
    sum_pairwise_acceptance = colSums(acceptance),
    row.names = NULL
  )
  # Mean second-largest alpha is the probability of at least two successful
  # *hypothetical* tests using one common U~Uniform(0,1), conditional on these
  # current states and uniformly choosing one reserve. No assignment occurs.
  shared_uniform_two <- if (n_patient < 2L) 0 else mean(apply(
    acceptance, 2L, function(a) sort(a, decreasing = TRUE)[[2L]]
  ))
  list(
    schema = "sab_shared_pool_probe_v1",
    useful_acceptance = useful_acceptance,
    log_w_current = log_w_current,
    log_w_reserve = log_w_reserve,
    acceptance = acceptance,
    gibbs_probabilities = gibbs,
    patient = patient,
    reserve = reserve,
    summary = data.frame(
      patient_count = n_patient,
      reserve_count = n_reserve,
      mean_reserve_mh_acceptance = mean(acceptance),
      fraction_reserves_useful_to_any_patient = mean(useful_counts >= 1L),
      fraction_reserves_useful_to_multiple_patients = mean(useful_counts >= 2L),
      hypothetical_common_uniform_two_patient_probability = shared_uniform_two,
      expected_accepted_squared_jump = if (is.null(squared_jump)) NA_real_ else
        mean(movement)
    )
  )
}
