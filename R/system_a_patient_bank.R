# Exact-target fixed-anchor patient MCMC banks for System A.
#
# This file only defines functions.  In particular, sourcing it never loads the
# System A target and never starts a chain.  The caller must supply a validated
# adapter, an anchor, and explicit run settings.

#' Draw one exact-target System A patient-conditional MCMC bank
#'
#' The proposal is pCN with respect to System A's diagonal Gaussian anchor
#' population density `g_i(x | eta, z_i)`.  Consequently `g_i` and the pCN
#' proposal densities cancel exactly, leaving only the sealed
#' patient-likelihood ratio in the Metropolis correction.  `beta` is adapted
#' during warm-up blocks and then frozen.  This diagonal-Gaussian requirement
#' belongs only to this System A bank kernel; raw candidate messages can use any
#' normalized evaluable population density on the same local-state space.
#' Compatibility between the adapter density and the pCN reference is checked
#' at every evaluated state; a mismatch is a contract error, not an MH
#' rejection.  Retained states are finite-chain MCMC output, not iid exact
#' conditional draws, so convergence must be assessed across chains.
#'
#' @param adapter A loaded System A adapter.  Tests may supply an object with the
#'   same callback contract.
#' @param patient_id One patient exposed by `adapter`.
#' @param eta,psi Finite anchor vectors in canonical adapter order.
#' @param initial_candidates A canonically named vector or numeric matrix whose
#'   rows are dispersed starting candidates.  Candidates are tried in order
#'   until one has finite exact likelihood.
#' @param warmup Number of warm-up proposals.
#' @param draws Number of retained draws.
#' @param thin Number of post-warm-up proposals per retained draw.
#' @param initial_beta Initial pCN innovation scale in `[0.005, 0.95]`.
#' @param target_acceptance Warm-up adaptation target in `(0, 1)`.
#' @param adaptation_block Positive warm-up block length.
#' @param seed Integer random seed.
#' @param proposal_mode Either `"pcn"` for the System A diagonal-Gaussian
#'   population law at `eta`, or `"defensive_independence"` for an equal
#'   mixture of that law and the law at `defensive_eta`.  The latter proposes
#'   independently from the normalized mixture, so the mixture cancels and the
#'   exact MH correction is again the sealed likelihood ratio.
#' @param defensive_eta Canonically named second population parameter vector,
#'   required only for `"defensive_independence"`.
#'
#' @return A `sab_system_a_patient_bank` containing draws, cached likelihoods,
#'   first- and second-moment ESS diagnostics, tuning history, and an exact-call
#'   ledger.
#' @export
sab_build_system_a_patient_bank <- function(
    adapter, patient_id, eta, psi, initial_candidates,
    warmup, draws, thin = 1L, initial_beta = 0.2,
    target_acceptance = 0.3, adaptation_block = 50L, seed,
    proposal_mode = c("pcn", "defensive_independence"),
    defensive_eta = NULL) {
  .sab_patient_bank_validate_adapter(adapter)
  local_names <- adapter$coordinate_names$local
  eta <- .sab_patient_bank_named_vector(
    eta, adapter$coordinate_names$population, "eta"
  )
  proposal_mode <- match.arg(proposal_mode)
  if (proposal_mode == "defensive_independence") {
    defensive_eta <- .sab_patient_bank_named_vector(
      defensive_eta, adapter$coordinate_names$population, "defensive_eta"
    )
  } else if (!is.null(defensive_eta)) {
    stop("defensive_eta is only valid for defensive_independence mode.",
         call. = FALSE)
  }
  psi <- .sab_patient_bank_named_vector(
    psi, adapter$coordinate_names$global, "psi"
  )
  patient_id <- as.character(patient_id)
  if (length(patient_id) != 1L || is.na(patient_id) ||
      !patient_id %in% adapter$patient_ids) {
    stop("patient_id is not in the supplied System A adapter.", call. = FALSE)
  }
  if (!isTRUE(adapter$eta_in_domain(eta)) ||
      (proposal_mode == "defensive_independence" &&
       !isTRUE(adapter$eta_in_domain(defensive_eta))) ||
      !isTRUE(adapter$psi_in_domain(psi))) {
    stop("The fixed patient-bank anchor is outside the target domain.",
         call. = FALSE)
  }
  warmup <- .sab_patient_bank_integer(warmup, "warmup", minimum = 0L)
  draws <- .sab_patient_bank_integer(draws, "draws", minimum = 1L)
  thin <- .sab_patient_bank_integer(thin, "thin", minimum = 1L)
  adaptation_block <- .sab_patient_bank_integer(
    adaptation_block, "adaptation_block", minimum = 1L
  )
  seed <- .sab_patient_bank_integer(seed, "seed", minimum = 0L)
  total_sampling <- as.double(draws) * thin
  total_iterations <- warmup + total_sampling
  if (!is.finite(total_iterations) ||
      total_iterations > .Machine$integer.max) {
    stop("The requested patient bank has too many iterations.", call. = FALSE)
  }
  if (!is.numeric(initial_beta) || length(initial_beta) != 1L ||
      !is.finite(initial_beta) || initial_beta < 0.005 ||
      initial_beta > 0.95) {
    stop("initial_beta must be finite and in [0.005, 0.95].",
         call. = FALSE)
  }
  if (!is.numeric(target_acceptance) || length(target_acceptance) != 1L ||
      !is.finite(target_acceptance) || target_acceptance <= 0 ||
      target_acceptance >= 1) {
    stop("target_acceptance must be finite and strictly between zero and one.",
         call. = FALSE)
  }

  initial_candidates <- .sab_patient_bank_candidates(
    initial_candidates, local_names
  )
  population_mean <- .sab_patient_bank_named_vector(
    adapter$population_mean(patient_id, eta),
    local_names,
    "patient population mean"
  )
  log_scale_names <- grep(
    "^log_omega_", adapter$coordinate_names$population, value = TRUE
  )
  if (length(log_scale_names) != length(local_names)) {
    stop("System A population scales do not match its local coordinates.",
         call. = FALSE)
  }
  population_sd <- exp(eta[log_scale_names])
  names(population_sd) <- local_names
  if (any(!is.finite(population_sd)) || any(population_sd <= 0)) {
    stop("Patient population scales are not representable.", call. = FALSE)
  }
  defensive_mean <- defensive_sd <- NULL
  if (proposal_mode == "defensive_independence") {
    defensive_mean <- .sab_patient_bank_named_vector(
      adapter$population_mean(patient_id, defensive_eta),
      local_names, "defensive patient population mean"
    )
    defensive_sd <- exp(defensive_eta[log_scale_names])
    names(defensive_sd) <- local_names
    if (any(!is.finite(defensive_sd)) || any(defensive_sd <= 0)) {
      stop("Defensive patient population scales are not representable.",
           call. = FALSE)
    }
  }

  phases <- c("initialization", "warmup", "sampling")
  ledger <- matrix(
    0, nrow = length(phases), ncol = 7L,
    dimnames = list(
      phases,
      c(
        "prediction_calls", "ode_integrations", "mh_proposals",
        "mh_accepted", "population_density_rejections",
        "prediction_failures", "nonfinite_loglik"
      )
    )
  )
  failure_reasons <- setNames(vector("list", length(phases)), phases)
  current <- NULL
  x <- NULL
  initial_candidate_index <- NA_integer_
  for (candidate_index in seq_len(nrow(initial_candidates))) {
    candidate <- initial_candidates[candidate_index, ]
    names(candidate) <- local_names
    evaluated <- .sab_patient_bank_evaluate(
      adapter, patient_id, candidate, eta, psi,
      population_mean, population_sd,
      reject_zero_population = proposal_mode == "pcn"
    )
    if (proposal_mode == "defensive_independence") {
      evaluated <- .sab_patient_bank_add_defensive_reference(
        evaluated, adapter, patient_id, candidate, defensive_eta,
        defensive_mean, defensive_sd
      )
    }
    ledger["initialization", "prediction_calls"] <-
      ledger["initialization", "prediction_calls"] +
      as.integer(evaluated$prediction_called)
    ledger["initialization", "ode_integrations"] <-
      ledger["initialization", "ode_integrations"] +
      evaluated$ode_integrations
    ledger["initialization", "population_density_rejections"] <-
      ledger["initialization", "population_density_rejections"] +
      as.integer(evaluated$population_density_rejection)
    ledger["initialization", "prediction_failures"] <-
      ledger["initialization", "prediction_failures"] +
      as.integer(evaluated$prediction_failure)
    ledger["initialization", "nonfinite_loglik"] <-
      ledger["initialization", "nonfinite_loglik"] +
      as.integer(!is.finite(evaluated$loglik))
    if (!identical(evaluated$reason, "ok")) {
      failure_reasons[["initialization"]] <- .sab_patient_bank_increment(
        failure_reasons[["initialization"]], evaluated$reason
      )
    }
    if (is.finite(evaluated$loglik)) {
      x <- candidate
      current <- evaluated
      initial_candidate_index <- candidate_index
      break
    }
  }
  if (is.null(current)) {
    stop("No initial patient candidate has finite exact likelihood: ",
         patient_id, ".", call. = FALSE)
  }

  set.seed(seed)
  beta <- as.numeric(initial_beta)
  beta_logit <- stats::qlogis(beta)
  block_accepted <- 0L
  block_proposals <- 0L
  adaptation_index <- 0L
  adaptation <- data.frame(
    iteration = integer(), block_acceptance = numeric(), beta = numeric()
  )
  kept <- matrix(
    NA_real_, nrow = draws, ncol = length(local_names),
    dimnames = list(NULL, local_names)
  )
  kept_loglik <- numeric(draws)
  kept_index <- 0L
  proposal_components <- c(saem = 0L, priorcentral = 0L)
  component_ledger <- expand.grid(
    phase = c("warmup", "sampling"),
    component = names(proposal_components),
    stringsAsFactors = FALSE
  )
  component_ledger$proposals <- 0L
  component_ledger$accepted <- 0L
  component_ledger$prediction_failures <- 0L
  component_ledger$nonfinite_loglik <- 0L

  for (iteration in seq_len(as.integer(total_iterations))) {
    phase <- if (iteration <= warmup) "warmup" else "sampling"
    component_name <- NA_character_
    if (proposal_mode == "pcn") {
      persistence <- sqrt(1 - beta^2)
      proposed_x <- population_mean +
        persistence * (x - population_mean) +
        beta * population_sd * stats::rnorm(length(x))
    } else {
      component <- if (stats::runif(1L) < 0.5) 1L else 2L
      component_name <- names(proposal_components)[[component]]
      proposal_components[[component_name]] <-
        proposal_components[[component_name]] + 1L
      component_mean <- if (component == 1L) population_mean else defensive_mean
      component_sd <- if (component == 1L) population_sd else defensive_sd
      proposed_x <- component_mean + component_sd * stats::rnorm(length(x))
    }
    names(proposed_x) <- local_names
    proposed <- .sab_patient_bank_evaluate(
      adapter, patient_id, proposed_x, eta, psi,
      population_mean, population_sd,
      reject_zero_population = proposal_mode == "pcn"
    )
    if (proposal_mode == "defensive_independence") {
      proposed <- .sab_patient_bank_add_defensive_reference(
        proposed, adapter, patient_id, proposed_x, defensive_eta,
        defensive_mean, defensive_sd
      )
    }
    component_row <- if (proposal_mode == "defensive_independence") {
      which(
        component_ledger$phase == phase &
          component_ledger$component == component_name
      )
    } else {
      integer()
    }
    if (length(component_row)) {
      component_ledger$proposals[[component_row]] <-
        component_ledger$proposals[[component_row]] + 1L
      component_ledger$prediction_failures[[component_row]] <-
        component_ledger$prediction_failures[[component_row]] +
        as.integer(proposed$prediction_failure)
      component_ledger$nonfinite_loglik[[component_row]] <-
        component_ledger$nonfinite_loglik[[component_row]] +
        as.integer(!is.finite(proposed$loglik))
    }
    ledger[phase, "prediction_calls"] <-
      ledger[phase, "prediction_calls"] +
      as.integer(proposed$prediction_called)
    ledger[phase, "ode_integrations"] <-
      ledger[phase, "ode_integrations"] + proposed$ode_integrations
    ledger[phase, "mh_proposals"] <- ledger[phase, "mh_proposals"] + 1L
    ledger[phase, "population_density_rejections"] <-
      ledger[phase, "population_density_rejections"] +
      as.integer(proposed$population_density_rejection)
    ledger[phase, "prediction_failures"] <-
      ledger[phase, "prediction_failures"] +
      as.integer(proposed$prediction_failure)
    ledger[phase, "nonfinite_loglik"] <-
      ledger[phase, "nonfinite_loglik"] +
      as.integer(!is.finite(proposed$loglik))
    if (!identical(proposed$reason, "ok")) {
      failure_reasons[[phase]] <- .sab_patient_bank_increment(
        failure_reasons[[phase]], proposed$reason
      )
    }
    accepted <- log(stats::runif(1L)) < proposed$loglik - current$loglik
    if (accepted) {
      x <- proposed_x
      current <- proposed
      ledger[phase, "mh_accepted"] <-
        ledger[phase, "mh_accepted"] + 1L
      if (length(component_row)) {
        component_ledger$accepted[[component_row]] <-
          component_ledger$accepted[[component_row]] + 1L
      }
    }

    if (phase == "warmup") {
      if (proposal_mode == "pcn") {
        block_accepted <- block_accepted + as.integer(accepted)
        block_proposals <- block_proposals + 1L
        if (iteration %% adaptation_block == 0L || iteration == warmup) {
          adaptation_index <- adaptation_index + 1L
          observed <- block_accepted / block_proposals
          gain <- 0.8 / sqrt(adaptation_index)
          beta_logit <- beta_logit + gain *
            (observed - target_acceptance)
          beta <- min(0.95, max(0.005, stats::plogis(beta_logit)))
          beta_logit <- stats::qlogis(beta)
          adaptation <- rbind(
            adaptation,
            data.frame(
              iteration = iteration,
              block_acceptance = observed,
              beta = beta
            )
          )
          block_accepted <- 0L
          block_proposals <- 0L
        }
      }
    } else if ((iteration - warmup) %% thin == 0L) {
      kept_index <- kept_index + 1L
      kept[kept_index, ] <- x
      kept_loglik[[kept_index]] <- current$loglik
    }
  }
  if (kept_index != draws || any(!is.finite(kept)) ||
      any(!is.finite(kept_loglik))) {
    stop("The patient bank was not filled with finite retained states.",
         call. = FALSE)
  }

  ledger <- as.data.frame(ledger)
  ledger$phase <- rownames(ledger)
  rownames(ledger) <- NULL
  ledger <- ledger[, c(
    "phase", "prediction_calls", "ode_integrations", "mh_proposals",
    "mh_accepted", "population_density_rejections",
    "prediction_failures", "nonfinite_loglik"
  )]
  total_row <- data.frame(
    phase = "total",
    prediction_calls = sum(ledger$prediction_calls),
    ode_integrations = sum(ledger$ode_integrations),
    mh_proposals = sum(ledger$mh_proposals),
    mh_accepted = sum(ledger$mh_accepted),
    population_density_rejections =
      sum(ledger$population_density_rejections),
    prediction_failures = sum(ledger$prediction_failures),
    nonfinite_loglik = sum(ledger$nonfinite_loglik)
  )
  ledger <- rbind(ledger, total_row)
  failure_reasons <- .sab_patient_bank_reason_table(failure_reasons)
  ess_x <- apply(kept, 2L, .sab_patient_bank_ess)
  ess_x_squared <- apply(kept^2, 2L, .sab_patient_bank_ess)

  structure(
    list(
      schema_version = "sab_system_a_patient_bank_v1",
      patient_id = patient_id,
      target_fingerprint = .sab_patient_bank_scalar_character(
        adapter$target_fingerprint
      ),
      anchor = list(
        eta = eta, psi = psi,
        population_mean = population_mean,
        population_sd = population_sd,
        population_reference = proposal_mode,
        pcn_reference = if (proposal_mode == "pcn") {
          "diagonal_gaussian"
        } else {
          NA_character_
        },
        pcn_reference_checked = proposal_mode == "pcn",
        defensive_eta = defensive_eta,
        defensive_population_mean = defensive_mean,
        defensive_population_sd = defensive_sd,
        defensive_mixture_weights = if (
          proposal_mode == "defensive_independence"
        ) c(saem = 0.5, priorcentral = 0.5) else NULL,
        defensive_reference_checked =
          proposal_mode == "defensive_independence"
      ),
      draws = kept,
      loglik = kept_loglik,
      initial_x = initial_candidates[initial_candidate_index, ],
      final_x = x,
      final_loglik = current$loglik,
      initial_candidate_index = initial_candidate_index,
      beta = list(
        initial = if (proposal_mode == "pcn") {
          as.numeric(initial_beta)
        } else NA_real_,
        final = if (proposal_mode == "pcn") beta else NA_real_,
        target_acceptance = as.numeric(target_acceptance),
        adaptation_block = adaptation_block,
        history = adaptation
      ),
      acceptance = c(
        warmup = if (warmup > 0L) {
          ledger$mh_accepted[ledger$phase == "warmup"] / warmup
        } else NA_real_,
        sampling = ledger$mh_accepted[ledger$phase == "sampling"] /
          total_sampling
      ),
      ess = list(x = ess_x, x_squared = ess_x_squared),
      ledger = ledger,
      exact_prediction_calls =
        ledger$prediction_calls[ledger$phase == "total"],
      exact_ode_integrations =
        ledger$ode_integrations[ledger$phase == "total"],
      rejection_reasons = failure_reasons,
      proposal = list(
        mode = proposal_mode,
        component_counts = if (
          proposal_mode == "defensive_independence"
        ) proposal_components else NULL,
        component_ledger = if (
          proposal_mode == "defensive_independence"
        ) component_ledger else NULL
      ),
      run = list(
        warmup = warmup, draws = draws, thin = thin, seed = seed,
        rng_kind = RNGkind(), r_version = R.version.string
      )
    ),
    class = c("sab_system_a_patient_bank", "list")
  )
}

.sab_patient_bank_validate_adapter <- function(adapter) {
  callbacks <- c(
    "solve_prediction", "loglik_from_prediction", "population_mean",
    "log_population_density", "eta_in_domain", "psi_in_domain"
  )
  coordinates <- c("local", "population", "global")
  valid_class <- inherits(adapter, "sab_system_a_adapter") ||
    inherits(adapter, "sab_patient_bank_mock_adapter")
  valid <- valid_class && is.list(adapter) &&
    is.character(adapter$patient_ids) && length(adapter$patient_ids) >= 1L &&
    !anyNA(adapter$patient_ids) && !anyDuplicated(adapter$patient_ids) &&
    is.list(adapter$coordinate_names) &&
    all(coordinates %in% names(adapter$coordinate_names)) &&
    all(vapply(adapter$coordinate_names[coordinates], function(value) {
      is.character(value) && length(value) >= 1L && !anyNA(value) &&
        !anyDuplicated(value)
    }, logical(1L))) &&
    all(vapply(callbacks, function(name) {
      is.function(adapter[[name]])
    }, logical(1L)))
  if (!isTRUE(valid)) {
    stop("Malformed System A patient-bank adapter.", call. = FALSE)
  }
  if (inherits(adapter, "sab_system_a_adapter")) {
    if (!exists("sab_validate_system_a_adapter", mode = "function")) {
      stop("Load the System A adapter module before building patient banks.",
           call. = FALSE)
    }
    sab_validate_system_a_adapter(adapter)
  }
  invisible(adapter)
}

.sab_patient_bank_named_vector <- function(value, expected_names, label) {
  original_names <- names(value)
  if ((!is.numeric(value) && !is.integer(value)) ||
      length(value) != length(expected_names) ||
      !identical(original_names, expected_names) ||
      any(!is.finite(value))) {
    stop(label, " must be finite and named in canonical order.",
         call. = FALSE)
  }
  value <- as.numeric(value)
  names(value) <- expected_names
  value
}

.sab_patient_bank_integer <- function(value, label, minimum) {
  if ((!is.numeric(value) && !is.integer(value)) || length(value) != 1L ||
      !is.finite(value) || value != floor(value) || value < minimum ||
      value > .Machine$integer.max) {
    stop(label, " must be one integer >= ", minimum, ".", call. = FALSE)
  }
  as.integer(value)
}

.sab_patient_bank_candidates <- function(value, local_names) {
  if (is.numeric(value) && is.null(dim(value))) {
    value <- matrix(value, nrow = 1L, dimnames = list(NULL, names(value)))
  }
  if (!is.matrix(value) || !is.numeric(value) || nrow(value) < 1L ||
      !identical(colnames(value), local_names) || any(!is.finite(value))) {
    stop("initial_candidates must be finite and in canonical local order.",
         call. = FALSE)
  }
  value
}

.sab_patient_bank_evaluate <- function(adapter, patient_id, x, eta, psi,
                                       population_mean, population_sd,
                                       reject_zero_population = TRUE) {
  population_log_density <- adapter$log_population_density(
    patient_id, x, eta
  )
  if (!is.numeric(population_log_density) ||
      length(population_log_density) != 1L ||
      is.na(population_log_density) || is.nan(population_log_density) ||
      population_log_density == Inf) {
    stop("The patient population density returned an invalid value.",
         call. = FALSE)
  }
  expected_log_density <- sum(stats::dnorm(
    x, mean = population_mean, sd = population_sd, log = TRUE
  ))
  both_negative_infinity <-
    population_log_density == -Inf && expected_log_density == -Inf
  comparison_scale <- max(
    1, abs(population_log_density), abs(expected_log_density)
  )
  reference_matches <- both_negative_infinity || (
    is.finite(population_log_density) && is.finite(expected_log_density) &&
      abs(population_log_density - expected_log_density) <=
        4096 * .Machine$double.eps * comparison_scale
  )
  if (!isTRUE(reference_matches)) {
    stop(
      "The adapter population density is incompatible with the diagonal-",
      "Gaussian pCN reference.", call. = FALSE
    )
  }
  if (population_log_density == -Inf && isTRUE(reject_zero_population)) {
    return(list(
      loglik = -Inf,
      population_log_density = -Inf,
      population_density_rejection = TRUE,
      prediction_called = FALSE,
      prediction_failure = FALSE,
      ode_integrations = 0L,
      reason = "nonfinite_population_density"
    ))
  }

  prediction <- adapter$solve_prediction(patient_id, x, psi)
  # Always route prediction states through the sealed callback.  It owns the
  # whitelist of ordinary target rejections and must distinguish those from
  # malformed/programming-error states.
  loglik <- adapter$loglik_from_prediction(prediction, psi)
  if (!is.list(prediction) || !is.logical(prediction$ok) ||
      length(prediction$ok) != 1L || is.na(prediction$ok)) {
    stop("The exact patient prediction returned a malformed state.",
         call. = FALSE)
  }
  if (!is.numeric(loglik) || length(loglik) != 1L || is.na(loglik) ||
      is.nan(loglik) || loglik == Inf) {
    stop("The exact patient likelihood returned an invalid value.",
         call. = FALSE)
  }
  prediction_failure <- !isTRUE(prediction$ok)
  reason <- if (prediction_failure) {
    if (!is.character(prediction$reason) || length(prediction$reason) != 1L ||
        is.na(prediction$reason) || !nzchar(prediction$reason)) {
      stop("The exact patient prediction has no valid failure reason.",
           call. = FALSE)
    }
    prediction$reason
  } else if (is.finite(loglik)) {
    "ok"
  } else {
    "nonfinite_loglik"
  }
  list(
    loglik = as.numeric(loglik),
    population_log_density = as.numeric(population_log_density),
    population_density_rejection = FALSE,
    prediction_called = TRUE,
    prediction_failure = prediction_failure,
    ode_integrations = .sab_patient_bank_ode_integrations(prediction),
    reason = reason
  )
}

.sab_patient_bank_add_defensive_reference <- function(
    evaluated, adapter, patient_id, x, defensive_eta,
    defensive_mean, defensive_sd) {
  defensive_log_density <- adapter$log_population_density(
    patient_id, x, defensive_eta
  )
  if (!is.numeric(defensive_log_density) ||
      length(defensive_log_density) != 1L ||
      is.na(defensive_log_density) || is.nan(defensive_log_density) ||
      defensive_log_density == Inf) {
    stop("The defensive population density returned an invalid value.",
         call. = FALSE)
  }
  expected <- sum(stats::dnorm(
    x, mean = defensive_mean, sd = defensive_sd, log = TRUE
  ))
  both_negative_infinity <-
    defensive_log_density == -Inf && expected == -Inf
  comparison_scale <- max(1, abs(defensive_log_density), abs(expected))
  matches <- both_negative_infinity || (
    is.finite(defensive_log_density) && is.finite(expected) &&
      abs(defensive_log_density - expected) <=
        4096 * .Machine$double.eps * comparison_scale
  )
  if (!isTRUE(matches)) {
    stop(
      "The defensive adapter density is incompatible with its diagonal-",
      "Gaussian independence component.", call. = FALSE
    )
  }
  first <- evaluated$population_log_density + log(0.5)
  second <- defensive_log_density + log(0.5)
  largest <- max(first, second)
  mixture <- if (largest == -Inf) {
    -Inf
  } else {
    largest + log(exp(first - largest) + exp(second - largest))
  }
  if (!is.finite(mixture)) {
    stop("The defensive reference mixture is non-finite.", call. = FALSE)
  }
  evaluated$population_log_density <- mixture
  evaluated
}

.sab_patient_bank_ode_integrations <- function(prediction) {
  if (!is.null(prediction$ode_integrations)) {
    value <- prediction$ode_integrations
    if ((!is.numeric(value) && !is.integer(value)) || length(value) != 1L ||
        !is.finite(value) || value != floor(value) || value < 0L) {
      stop("Prediction `ode_integrations` must be one non-negative integer.",
           call. = FALSE)
    }
    return(as.integer(value))
  }
  if (isTRUE(prediction$ok)) {
    times <- prediction$adjusted_positive_times
    if (!is.numeric(times) || any(!is.finite(times))) {
      stop("Successful prediction omits valid adjusted positive times.",
           call. = FALSE)
    }
    return(as.integer(length(times) > 0L))
  }
  if (prediction$reason %in%
      c("ode_failure", "invalid_ode_output", "invalid_ode_times")) {
    return(1L)
  }
  0L
}

.sab_patient_bank_increment <- function(counts, label) {
  if (is.null(counts)) counts <- integer()
  if (!label %in% names(counts)) counts[[label]] <- 0L
  counts[[label]] <- counts[[label]] + 1L
  counts
}

.sab_patient_bank_reason_table <- function(by_phase) {
  rows <- lapply(names(by_phase), function(phase) {
    counts <- by_phase[[phase]]
    if (!length(counts)) return(NULL)
    data.frame(
      phase = rep(phase, length(counts)),
      reason = names(counts), count = as.integer(counts),
      stringsAsFactors = FALSE
    )
  })
  rows <- Filter(Negate(is.null), rows)
  if (!length(rows)) {
    return(data.frame(
      phase = character(), reason = character(), count = integer()
    ))
  }
  do.call(rbind, rows)
}

.sab_patient_bank_ess <- function(value) {
  value <- as.numeric(value)
  n <- length(value)
  variance <- if (n >= 2L) stats::var(value) else NA_real_
  if (n < 4L || any(!is.finite(value)) || !is.finite(variance) ||
      variance <= 0) return(0)
  maximum_lag <- min(n - 1L, max(20L, floor(n / 2L)))
  correlation <- as.numeric(stats::acf(
    value, lag.max = maximum_lag, plot = FALSE,
    type = "correlation", demean = TRUE
  )$acf)[-1L]
  if (length(correlation) < 2L) return(as.numeric(n))
  correlation <- correlation[
    seq_len(2L * floor(length(correlation) / 2L))
  ]
  pair_sums <- correlation[seq(1L, length(correlation), by = 2L)] +
    correlation[seq(2L, length(correlation), by = 2L)]
  first_nonpositive <- which(!is.finite(pair_sums) | pair_sums <= 0)[1L]
  if (!is.na(first_nonpositive)) {
    pair_sums <- head(pair_sums, first_nonpositive - 1L)
  }
  if (length(pair_sums) > 1L) pair_sums <- cummin(pair_sums)
  integrated_time <- max(1, 1 + 2 * sum(pair_sums))
  min(n, n / integrated_time)
}

.sab_patient_bank_scalar_character <- function(value) {
  if (is.character(value) && length(value) == 1L && !is.na(value)) {
    value
  } else {
    NA_character_
  }
}
