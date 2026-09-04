# Final, bounded pure-SAEM fixed-psi message falsification for System A.
#
# This module deliberately does not know how banks are simulated.  It consumes
# frozen exact-target patient-bank artifacts, evaluates normalized population
# density ratios without ODE calls, and compares forward importance estimates
# with independent reverse banks and bidirectional bridge estimates.  The
# dynamic shared parameters remain fixed at each SAEM branch's own anchor.

.sab_pure_required_functions <- function() {
  required <- c(
    "sab_evaluate_population_message", "sab_raw_message_diagnostics",
    "sab_mcmc_scalar_diagnostics", "sab_batch_log_mcse",
    "sab_log_mean_exp", ".sab_me_bridge",
    ".sab_me_bridge_scores", ".sab_me_rescaled_exp_chains"
  )
  missing <- required[!vapply(required, exists, logical(1L), mode = "function")]
  if (length(missing)) {
    stop(
      "Pure-message validation dependencies are not loaded: ",
      paste(missing, collapse = ", "), ".",
      call. = FALSE
    )
  }
  invisible(TRUE)
}

#' Frozen design for the final pure-SAEM System A falsification
#'
#' The treatment-effect displacement is resolved at each anchor to one fitted
#' population SD of `u_eta_pi`.  The scale displacement changes omega_lambda by
#' a factor of 1.25.  This is below the Gaussian raw-weight second-moment
#' boundary sqrt(2), unlike the previously suggested exp(0.35) widening.
#'
#' @return A fixed design and strict go/no-go thresholds.
#' @export
sab_system_a_pure_message_design <- function() {
  list(
    schema_version = "sab_system_a_pure_message_design_v1",
    patient_ids = c(
      "3", "6", "16", "20", "55", "71", "74", "88", "105",
      "111", "117", "122"
    ),
    pilot_chains = c("chain_01", "chain_02"),
    # These sets are deliberately disjoint within the final assessment.
    # Forward and reverse raw importance therefore provide an independent
    # two-bank comparison, while the bridge uses neither raw estimator's
    # chains.  Pilot chains are reused by the bridge only after the pilot has
    # passed; the endpoints themselves are fixed below rather than selected by
    # that pilot.
    forward_chains = c("chain_03", "chain_04"),
    reverse_chains = c("chain_01", "chain_02"),
    bridge_reference_chains = c("chain_01", "chain_02"),
    bridge_candidate_chains = c("chain_03", "chain_04"),
    bank_contract = list(
      chains = 4L,
      chain_names = sprintf("chain_%02d", 1:4),
      warmup = 250L,
      draws = 500L,
      thin = 1L,
      initial_beta = 0.4,
      target_acceptance = 0.3,
      adaptation_block = 50L,
      start_offset_sd = 1,
      maximum_initial_candidates_per_chain = 18L,
      maximum_prediction_calls_per_patient_target = 3072
    ),
    endpoint_template = data.frame(
      endpoint_id = c(
        "treatment_minus", "treatment_plus",
        "lambda_scale_minus", "lambda_scale_plus"
      ),
      stage = rep("final_pure_saem_falsification", 4L),
      kind = rep("population_shift", 4L),
      axis = c(
        "treatment_effect", "treatment_effect",
        "lambda_population_scale", "lambda_population_scale"
      ),
      parameter = c(
        "beta_nelf", "beta_nelf",
        "log_omega_lambda", "log_omega_lambda"
      ),
      sign = c(-1L, 1L, -1L, 1L),
      step_rule = c(
        "one_anchor_omega_eta_pi", "one_anchor_omega_eta_pi",
        "log_1.25", "log_1.25"
      ),
      stringsAsFactors = FALSE
    ),
    projection_patient_count = 115L,
    full_treatment_counts = c(`0` = 80L, `1` = 35L),
    invariant_control_log_tolerance = 1e-8,
    pilot_gates = list(
      minimum_relative_weight_ess = 0.20,
      maximum_normalized_weight = 0.05,
      maximum_split_log_ratio_difference = 0.25,
      maximum_cohort_chain_range = 0.50,
      maximum_split_rhat = 1.01,
      minimum_relevant_mcmc_ess = 400,
      minimum_relevant_tail_ess = 200,
      maximum_projected_115_log_mcse = 0.50,
      maximum_projected_variance_share = c(
        treatment_effect = 0.50,
        lambda_population_scale = 0.25
      ),
      maximum_batch_mcse_ratio = 2
    ),
    evaluation_gates = list(
      minimum_relative_weight_ess = 0.20,
      maximum_normalized_weight = 0.05,
      maximum_split_log_ratio_difference = 0.25,
      maximum_cohort_chain_range = 0.50,
      minimum_agreement_tolerance = 0.25,
      minimum_cohort_agreement_tolerance = 0.50,
      maximum_split_rhat = 1.01,
      minimum_relevant_mcmc_ess = 400,
      minimum_relevant_tail_ess = 200,
      maximum_projected_115_log_mcse = 0.50,
      maximum_projected_variance_share = c(
        treatment_effect = 0.50,
        lambda_population_scale = 0.25
      ),
      maximum_batch_mcse_ratio = 2
    )
  )
}

.sab_pure_reference <- function(canonical_anchor) {
  list(kind = "population", eta = canonical_anchor$eta)
}

.sab_pure_density <- function(adapter) {
  force(adapter)
  function(x, parameter, patient_context) {
    if (!is.list(parameter) || !identical(parameter$kind, "population") ||
        !is.numeric(parameter$eta) ||
        !identical(names(parameter$eta),
                   adapter$coordinate_names$population)) {
      stop("Pure-message population parameter is malformed.", call. = FALSE)
    }
    adapter$log_population_density(
      patient_context$patient_id, x, parameter$eta
    )
  }
}

.sab_pure_endpoints <- function(adapter, canonical_anchor) {
  eta <- canonical_anchor$eta
  required <- c("beta_nelf", "log_omega_lambda", "log_omega_eta_pi")
  if (!is.numeric(eta) ||
      !identical(names(eta), adapter$coordinate_names$population) ||
      any(!is.finite(eta)) || length(setdiff(required, names(eta)))) {
    stop("The canonical anchor cannot resolve pure-message endpoints.",
         call. = FALSE)
  }
  beta_step <- exp(eta[["log_omega_eta_pi"]])
  scale_step <- log(1.25)
  if (!is.finite(beta_step) || beta_step <= 0 ||
      !is.finite(scale_step) || scale_step <= 0 ||
      exp(scale_step) >= sqrt(2)) {
    stop("The resolved pure-message endpoint steps are invalid.",
         call. = FALSE)
  }
  template <- sab_system_a_pure_message_design()$endpoint_template
  endpoints <- template[, c(
    "endpoint_id", "stage", "kind", "axis", "parameter", "sign"
  )]
  endpoints$step <- c(beta_step, beta_step, scale_step, scale_step)
  endpoint_eta <- matrix(
    rep(eta, times = nrow(endpoints)), nrow = nrow(endpoints), byrow = TRUE,
    dimnames = list(endpoints$endpoint_id, names(eta))
  )
  for (index in seq_len(nrow(endpoints))) {
    parameter <- endpoints$parameter[[index]]
    endpoint_eta[index, parameter] <- eta[[parameter]] +
      endpoints$sign[[index]] * endpoints$step[[index]]
    candidate <- endpoint_eta[index, ]
    names(candidate) <- colnames(endpoint_eta)
    if (!isTRUE(adapter$eta_in_domain(candidate))) {
      stop("Resolved endpoint is outside the System A population domain: ",
           endpoints$endpoint_id[[index]], ".", call. = FALSE)
    }
  }
  endpoints$pilot_gate_passed <- FALSE
  endpoints$diagnostic_fallback <- FALSE
  list(endpoints = endpoints, endpoint_eta = endpoint_eta)
}

.sab_pure_patient_context <- function(artifact, adapter, patient_id) {
  context <- artifact$patient_context$population_covariates
  valid <- is.numeric(context) && identical(names(context), "treat_nelf") &&
    length(context) == 1L && context[[1L]] %in% c(0, 1) &&
    identical(unname(context[[1L]]),
              unname(adapter$treatment[[patient_id]]))
  if (!isTRUE(valid)) {
    stop("Patient bank has malformed or stale treatment context: ",
         patient_id, ".", call. = FALSE)
  }
  list(patient_id = patient_id, treat_nelf = unname(context[[1L]]))
}

.sab_pure_validate_chain <- function(chain, adapter, patient_id, eta, psi,
                                     bank_contract) {
  valid <- is.list(chain) && inherits(chain, "sab_system_a_patient_bank") &&
    identical(chain$patient_id, patient_id) &&
    is.matrix(chain$draws) && is.numeric(chain$draws) &&
    nrow(chain$draws) >= 8L &&
    identical(colnames(chain$draws), adapter$coordinate_names$local) &&
    all(is.finite(chain$draws)) &&
    is.numeric(chain$loglik) && length(chain$loglik) == nrow(chain$draws) &&
    all(is.finite(chain$loglik)) &&
    identical(chain$anchor$eta, eta) && identical(chain$anchor$psi, psi) &&
    identical(chain$anchor$population_reference, "pcn") &&
    identical(chain$anchor$pcn_reference, "diagonal_gaussian") &&
    isTRUE(chain$anchor$pcn_reference_checked) &&
    identical(chain$proposal$mode, "pcn") &&
    is.list(chain$beta) &&
    isTRUE(all.equal(chain$beta$initial, bank_contract$initial_beta)) &&
    isTRUE(all.equal(
      chain$beta$target_acceptance, bank_contract$target_acceptance
    )) &&
    identical(chain$beta$adaptation_block, bank_contract$adaptation_block) &&
    is.list(chain$run) &&
    identical(chain$run$warmup, bank_contract$warmup) &&
    identical(chain$run$draws, bank_contract$draws) &&
    identical(chain$run$thin, bank_contract$thin) &&
    nrow(chain$draws) == bank_contract$draws
  if (!isTRUE(valid)) {
    stop("Malformed pure pCN chain for patient ", patient_id, ".",
         call. = FALSE)
  }
  invisible(TRUE)
}

.sab_pure_bank_budget_status <- function(artifact, design) {
  contract <- design$bank_contract
  configuration <- artifact$configuration
  summary <- artifact$summary
  configuration_ok <- is.list(configuration) &&
    identical(configuration$chains, contract$chains) &&
    identical(configuration$warmup, contract$warmup) &&
    identical(configuration$draws, contract$draws) &&
    identical(configuration$thin, contract$thin) &&
    isTRUE(all.equal(configuration$initial_beta, contract$initial_beta)) &&
    isTRUE(all.equal(
      configuration$target_acceptance, contract$target_acceptance
    )) &&
    identical(configuration$adaptation_block, contract$adaptation_block) &&
    isTRUE(all.equal(
      configuration$start_offset_sd, contract$start_offset_sd
    )) &&
    identical(configuration$maximum_initial_candidates_per_chain,
              contract$maximum_initial_candidates_per_chain) &&
    is.numeric(configuration$exact_prediction_call_cap) &&
    length(configuration$exact_prediction_call_cap) == 1L &&
    is.finite(configuration$exact_prediction_call_cap) &&
    configuration$exact_prediction_call_cap ==
      contract$maximum_prediction_calls_per_patient_target
  summary_ok <- is.data.frame(summary) && nrow(summary) == contract$chains &&
    identical(as.character(summary$chain), contract$chain_names) &&
    all(c("exact_prediction_calls", "exact_ode_integrations") %in%
          names(summary)) &&
    is.numeric(summary$exact_prediction_calls) &&
    all(is.finite(summary$exact_prediction_calls)) &&
    all(summary$exact_prediction_calls >=
          contract$warmup + contract$draws * contract$thin + 1L) &&
    is.numeric(summary$exact_ode_integrations) &&
    all(is.finite(summary$exact_ode_integrations)) &&
    all(summary$exact_ode_integrations >= 0)
  observed_calls <- if (summary_ok) {
    sum(summary$exact_prediction_calls)
  } else {
    Inf
  }
  observed_ode_integrations <- if (summary_ok) {
    sum(summary$exact_ode_integrations)
  } else {
    Inf
  }
  list(
    configuration_ok = configuration_ok,
    summary_ok = summary_ok,
    observed_prediction_calls = observed_calls,
    observed_ode_integrations = observed_ode_integrations,
    maximum_prediction_calls =
      contract$maximum_prediction_calls_per_patient_target,
    within_budget = configuration_ok && summary_ok &&
      observed_calls <= contract$maximum_prediction_calls_per_patient_target &&
      observed_ode_integrations <=
        contract$maximum_prediction_calls_per_patient_target
  )
}

.sab_pure_validate_reference_artifacts <- function(
    adapter, canonical_anchor, artifacts, design) {
  reference <- .sab_pure_reference(canonical_anchor)
  if (!is.list(artifacts) || !identical(names(artifacts), design$patient_ids)) {
    stop("Pure-SAEM reference artifacts are not in audited patient order.",
         call. = FALSE)
  }
  required_chains <- design$bank_contract$chain_names
  for (patient_id in design$patient_ids) {
    artifact <- artifacts[[patient_id]]
    target <- artifact$conditional_target
    valid <- is.list(artifact) &&
      identical(artifact$schema_version,
                "sab_system_a_pure_reference_patient_banks_v1") &&
      identical(artifact$bank_role, "pure_saem_reference") &&
      identical(artifact$patient_id, patient_id) &&
      identical(artifact$target$target_fingerprint,
                adapter$target_fingerprint) &&
      identical(target$eta, canonical_anchor$eta) &&
      identical(target$psi, canonical_anchor$psi) &&
      identical(target$reference, reference) &&
      is.list(artifact$chains) &&
      identical(names(artifact$chains), required_chains)
    if (!isTRUE(valid)) {
      stop("Malformed or incompatible pure-SAEM reference artifact: ",
           patient_id, ".", call. = FALSE)
    }
    .sab_pure_patient_context(artifact, adapter, patient_id)
    for (chain_name in required_chains) {
      .sab_pure_validate_chain(
        artifact$chains[[chain_name]], adapter, patient_id,
        canonical_anchor$eta, canonical_anchor$psi, design$bank_contract
      )
    }
  }
  invisible(TRUE)
}

.sab_pure_validate_candidate_artifact <- function(
    artifact, adapter, canonical_anchor, plan, plan_sha256,
    endpoint, candidate_eta, patient_id) {
  target <- artifact$conditional_target
  valid <- is.list(artifact) &&
    identical(artifact$schema_version,
              "sab_system_a_endpoint_patient_banks_v1") &&
    identical(artifact$bank_role, "candidate_endpoint") &&
    identical(artifact$patient_id, patient_id) &&
    identical(artifact$target$target_fingerprint, adapter$target_fingerprint) &&
    identical(target$endpoint$endpoint_id, endpoint$endpoint_id) &&
    identical(target$eta, candidate_eta) &&
    identical(target$psi, canonical_anchor$psi) &&
    identical(target$endpoint_plan$sha256, plan_sha256) &&
    is.list(artifact$chains) &&
    identical(names(artifact$chains),
              sab_system_a_pure_message_design()$bank_contract$chain_names)
  if (!isTRUE(valid)) {
    stop("Malformed candidate artifact for ", endpoint$endpoint_id,
         ", patient ", patient_id, ".", call. = FALSE)
  }
  .sab_pure_patient_context(artifact, adapter, patient_id)
  for (chain in artifact$chains) {
    .sab_pure_validate_chain(
      chain, adapter, patient_id, candidate_eta, canonical_anchor$psi,
      sab_system_a_pure_message_design()$bank_contract
    )
  }
  invisible(TRUE)
}

.sab_pure_diagnostic_row <- function(endpoint, patient_id, treatment, chain,
                                     direction, diagnostic, log_weights) {
  if (!is.numeric(log_weights) || !length(log_weights) ||
      any(!is.finite(log_weights))) {
    stop("Message log weights are malformed.", call. = FALSE)
  }
  data.frame(
    endpoint_id = endpoint$endpoint_id,
    axis = endpoint$axis,
    parameter = endpoint$parameter,
    step = endpoint$step,
    sign = endpoint$sign,
    patient_id = patient_id,
    treatment = treatment,
    chain = chain,
    direction = direction,
    raw_log_ratio = diagnostic$log_ratio,
    oriented_log_ratio = if (direction == "reverse") {
      -diagnostic$log_ratio
    } else {
      diagnostic$log_ratio
    },
    relative_weight_ess = diagnostic$relative_weight_ess,
    max_normalized_weight = diagnostic$max_normalized_weight,
    d2 = diagnostic$d2,
    batch_log_mcse = diagnostic$batch$log_scale_mcse,
    split_log_ratio_difference =
      diagnostic$split$absolute_log_ratio_difference,
    split_d2_difference = diagnostic$split$absolute_d2_difference,
    maximum_absolute_log_weight = max(abs(log_weights)),
    stringsAsFactors = FALSE
  )
}

.sab_pure_relevant_state_quantities <- function(endpoint, context, evaluations,
                                                 prefix) {
  parameter <- endpoint$parameter[[1L]]
  if (parameter == "beta_nelf" && context$treat_nelf == 0L) {
    return(list())
  }
  if (parameter == "beta_nelf") {
    coordinate <- "u_eta_pi"
    return(stats::setNames(list(lapply(
      evaluations, function(value) value$chain$draws[, coordinate]
    )), paste0(prefix, "_x/", coordinate)))
  }
  if (parameter == "log_omega_lambda") {
    coordinate <- "u_log_lambda"
    raw <- lapply(evaluations, function(value) value$chain$draws[, coordinate])
    return(stats::setNames(
      list(raw, lapply(raw, function(value) value^2)),
      c(paste0(prefix, "_x/", coordinate),
        paste0(prefix, "_x_squared/", coordinate))
    ))
  }
  stop("Unknown pure-message endpoint parameter.", call. = FALSE)
}

.sab_pure_rank_normalize_chains <- function(chains, folded = FALSE) {
  lengths <- vapply(chains, length, integer(1L))
  if (!length(lengths) || length(unique(lengths)) != 1L ||
      any(!vapply(chains, function(value) {
        is.numeric(value) && all(is.finite(value))
      }, logical(1L)))) {
    stop("Rank diagnostics require equal-length finite chains.", call. = FALSE)
  }
  pooled <- unlist(chains, use.names = FALSE)
  if (folded) pooled <- abs(pooled - stats::median(pooled))
  ranks <- rank(pooled, ties.method = "average")
  normal_scores <- stats::qnorm(
    (ranks - 3 / 8) / (length(ranks) + 1 / 4)
  )
  starts <- cumsum(c(1L, head(lengths, -1L)))
  result <- lapply(seq_along(chains), function(index) {
    normal_scores[starts[[index]]:(starts[[index]] + lengths[[index]] - 1L)]
  })
  names(result) <- names(chains)
  result
}

.sab_pure_rank_tail_diagnostics <- function(chains) {
  rank_diagnostic <- sab_mcmc_scalar_diagnostics(
    .sab_pure_rank_normalize_chains(chains)
  )
  folded_diagnostic <- sab_mcmc_scalar_diagnostics(
    .sab_pure_rank_normalize_chains(chains, folded = TRUE)
  )
  pooled <- unlist(chains, use.names = FALSE)
  quantiles <- stats::quantile(
    pooled, probs = c(0.05, 0.95), names = FALSE, type = 8
  )
  lower <- lapply(chains, function(value) as.numeric(value <= quantiles[[1L]]))
  upper <- lapply(chains, function(value) as.numeric(value >= quantiles[[2L]]))
  tail_ess <- min(
    sab_mcmc_scalar_diagnostics(lower)$mcmc_ess,
    sab_mcmc_scalar_diagnostics(upper)$mcmc_ess
  )
  list(
    split_rhat = max(rank_diagnostic$split_rhat,
                     folded_diagnostic$split_rhat),
    rank_normalized_split_rhat = rank_diagnostic$split_rhat,
    folded_split_rhat = folded_diagnostic$split_rhat,
    bulk_mcmc_ess = rank_diagnostic$mcmc_ess,
    tail_mcmc_ess = tail_ess,
    relative_bulk_mcmc_ess = rank_diagnostic$relative_mcmc_ess,
    chain_mean_range = diff(range(vapply(chains, mean, numeric(1L))))
  )
}

.sab_pure_mcmc_rows <- function(quantities, endpoint_id, patient_id,
                                required_names) {
  do.call(rbind, lapply(names(quantities), function(quantity) {
    diagnostic <- .sab_pure_rank_tail_diagnostics(quantities[[quantity]])
    data.frame(
      endpoint_id = endpoint_id,
      patient_id = patient_id,
      quantity = quantity,
      required_for_gate = quantity %in% required_names,
      split_rhat = diagnostic$split_rhat,
      rank_normalized_split_rhat = diagnostic$rank_normalized_split_rhat,
      folded_split_rhat = diagnostic$folded_split_rhat,
      mcmc_ess = diagnostic$bulk_mcmc_ess,
      tail_mcmc_ess = diagnostic$tail_mcmc_ess,
      relative_mcmc_ess = diagnostic$relative_bulk_mcmc_ess,
      chain_mean_range = diagnostic$chain_mean_range,
      diagnostic_available = TRUE,
      failure_code = NA_character_,
      failure_message = NA_character_,
      stringsAsFactors = FALSE
    )
  }))
}

.sab_pure_patient_log_mcse <- function(evaluations) {
  combine <- function(n_batches) {
    values <- vapply(evaluations, function(value) {
      sab_batch_log_mcse(
        value$population$log_weights, n_batches = n_batches
      )$log_scale_mcse
    }, numeric(1L))
    sqrt(sum(values^2)) / length(values)
  }
  mcse_5 <- combine(5L)
  mcse_10 <- combine(10L)
  largest <- max(mcse_5, mcse_10)
  smallest <- min(mcse_5, mcse_10)
  ratio <- if (largest == 0) {
    1
  } else if (smallest == 0) {
    Inf
  } else {
    largest / smallest
  }
  list(
    conservative = largest,
    batch_5 = mcse_5,
    batch_10 = mcse_10,
    batch_ratio = ratio
  )
}

.sab_pure_projection_summary <- function(patient_mcse, treatment, design) {
  if (!is.numeric(patient_mcse) || length(patient_mcse) !=
      length(design$patient_ids) || any(!is.finite(patient_mcse)) ||
      any(patient_mcse < 0) || !is.numeric(treatment) ||
      length(treatment) != length(patient_mcse) || anyNA(treatment) ||
      !all(treatment %in% c(0, 1)) ||
      !identical(names(design$full_treatment_counts), c("0", "1")) ||
      sum(design$full_treatment_counts) != design$projection_patient_count) {
    stop("Stratified patient MCSE inputs are malformed.", call. = FALSE)
  }
  projected_contributions <- numeric(length(patient_mcse))
  for (stratum in names(design$full_treatment_counts)) {
    selected <- which(treatment == as.numeric(stratum))
    if (!length(selected)) {
      stop("The 12-patient panel omits a treatment stratum.", call. = FALSE)
    }
    projected_contributions[selected] <-
      design$full_treatment_counts[[stratum]] / length(selected) *
      patient_mcse[selected]^2
  }
  projected_variance <- sum(projected_contributions)
  maximum_share <- if (projected_variance == 0) {
    0
  } else {
    max(projected_contributions) / projected_variance
  }
  list(
    projected_log_mcse = sqrt(projected_variance),
    projected_variance = projected_variance,
    maximum_projected_variance_share = maximum_share,
    worst_observed_patient_envelope =
      sqrt(design$projection_patient_count) * max(patient_mcse),
    delta_log_bias_proxy = 0.5 * projected_variance,
    projected_contributions = projected_contributions
  )
}

.sab_pure_variance_share_limit <- function(endpoint, gates) {
  axis <- as.character(endpoint$axis[[1L]])
  limits <- gates$maximum_projected_variance_share
  if (!is.numeric(limits) || is.null(names(limits)) ||
      !axis %in% names(limits) || !is.finite(limits[[axis]]) ||
      limits[[axis]] <= 0 || limits[[axis]] > 1) {
    stop("The endpoint variance-concentration gate is malformed.",
         call. = FALSE)
  }
  unname(limits[[axis]])
}

.sab_pure_reference_performance <- function(artifacts) {
  rows <- do.call(rbind, lapply(artifacts, function(artifact) artifact$summary))
  rownames(rows) <- NULL
  required <- c(
    "patient_id", "chain", "sampling_acceptance", "min_ess_x",
    "min_ess_x_squared", "exact_prediction_calls", "exact_ode_integrations"
  )
  if (!is.data.frame(rows) || length(setdiff(required, names(rows))) ||
      any(!is.finite(rows$exact_prediction_calls)) ||
      any(!is.finite(rows$exact_ode_integrations))) {
    stop("Pure-SAEM bank summaries omit required cost diagnostics.",
         call. = FALSE)
  }
  rows
}

#' Plan the final pure-SAEM endpoints from pilot reference chains
#'
#' @return A frozen four-endpoint plan and ODE-free pilot diagnostics.
#' @export
sab_system_a_plan_pure_message_endpoints <- function(
    adapter, canonical_anchor, anchor_bank_artifacts,
    anchor_path, anchor_sha256, anchor_bank_paths, anchor_bank_sha256) {
  .sab_pure_required_functions()
  design <- sab_system_a_pure_message_design()
  .sab_pure_validate_reference_artifacts(
    adapter, canonical_anchor, anchor_bank_artifacts, design
  )
  if (!identical(names(anchor_bank_paths), design$patient_ids) ||
      !identical(names(anchor_bank_sha256), design$patient_ids) ||
      !is.character(anchor_path) || length(anchor_path) != 1L ||
      !is.character(anchor_sha256) || length(anchor_sha256) != 1L ||
      !grepl("^[0-9a-f]{64}$", anchor_sha256) ||
      any(!grepl("^[0-9a-f]{64}$", anchor_bank_sha256))) {
    stop("Pure-SAEM anchor paths or hashes are malformed.", call. = FALSE)
  }
  resolved <- .sab_pure_endpoints(adapter, canonical_anchor)
  endpoints <- resolved$endpoints
  endpoint_eta <- resolved$endpoint_eta
  density <- .sab_pure_density(adapter)
  reference <- .sab_pure_reference(canonical_anchor)
  diagnostic_rows <- list()
  mcmc_rows <- list()
  patient_mcse_rows <- list()
  diagnostic_index <- mcmc_index <- patient_index <- 0L

  for (endpoint_index in seq_len(nrow(endpoints))) {
    endpoint <- endpoints[endpoint_index, , drop = FALSE]
    candidate_eta <- endpoint_eta[endpoint$endpoint_id, ]
    names(candidate_eta) <- colnames(endpoint_eta)
    candidate <- list(kind = "population", eta = candidate_eta)
    for (patient_id in design$patient_ids) {
      artifact <- anchor_bank_artifacts[[patient_id]]
      context <- .sab_pure_patient_context(artifact, adapter, patient_id)
      evaluations <- lapply(design$pilot_chains, function(chain_name) {
        chain <- artifact$chains[[chain_name]]
        population <- sab_evaluate_population_message(
          chain$draws, reference, candidate, context, density
        )
        diagnostic <- sab_raw_message_diagnostics(population$log_weights)
        diagnostic_index <<- diagnostic_index + 1L
        diagnostic_rows[[diagnostic_index]] <<- .sab_pure_diagnostic_row(
          endpoint, patient_id, context$treat_nelf, chain_name, "forward",
          diagnostic, population$log_weights
        )
        list(chain = chain, population = population, diagnostic = diagnostic)
      })
      names(evaluations) <- design$pilot_chains
      quantities <- c(
        list(forward_raw_weight_contribution =
               .sab_me_rescaled_exp_chains(lapply(
                 evaluations, function(value) value$population$log_weights
               ))),
        .sab_pure_relevant_state_quantities(
          endpoint, context, evaluations, "reference"
        )
      )
      required <- names(quantities)
      rows <- .sab_pure_mcmc_rows(
        quantities, endpoint$endpoint_id, patient_id, required
      )
      mcmc_index <- mcmc_index + 1L
      mcmc_rows[[mcmc_index]] <- rows
      mcse <- .sab_pure_patient_log_mcse(evaluations)
      patient_index <- patient_index + 1L
      patient_mcse_rows[[patient_index]] <- data.frame(
        endpoint_id = endpoint$endpoint_id,
        patient_id = patient_id,
        treatment = context$treat_nelf,
        forward_log_mcse = mcse$conservative,
        forward_log_mcse_batch_5 = mcse$batch_5,
        forward_log_mcse_batch_10 = mcse$batch_10,
        forward_log_mcse_batch_ratio = mcse$batch_ratio,
        stringsAsFactors = FALSE
      )
    }
  }
  diagnostics <- do.call(rbind, diagnostic_rows)
  mcmc <- do.call(rbind, mcmc_rows)
  patient_mcse <- do.call(rbind, patient_mcse_rows)
  rownames(diagnostics) <- rownames(mcmc) <- rownames(patient_mcse) <- NULL

  gates <- design$pilot_gates
  budget_status <- lapply(
    anchor_bank_artifacts, .sab_pure_bank_budget_status, design = design
  )
  reference_cost_passed <- all(vapply(
    budget_status, function(value) value$within_budget, logical(1L)
  ))
  maximum_reference_calls <- max(vapply(
    budget_status,
    function(value) value$observed_prediction_calls,
    numeric(1L)
  ))
  selection <- do.call(rbind, lapply(endpoints$endpoint_id, function(id) {
    endpoint <- endpoints[endpoints$endpoint_id == id, , drop = FALSE]
    variance_share_limit <- .sab_pure_variance_share_limit(endpoint, gates)
    selected <- diagnostics[diagnostics$endpoint_id == id, , drop = FALSE]
    selected_mcmc <- mcmc[
      mcmc$endpoint_id == id & mcmc$required_for_gate, , drop = FALSE
    ]
    selected_mcse <- patient_mcse[patient_mcse$endpoint_id == id, ]
    cohort <- stats::aggregate(raw_log_ratio ~ chain, selected, sum)
    projection <- .sab_pure_projection_summary(
      selected_mcse$forward_log_mcse, selected_mcse$treatment, design
    )
    invariant_control <- selected$parameter == "beta_nelf" &
      selected$treatment == 0L
    invariant_control_passed <- all(
      !invariant_control |
        selected$maximum_absolute_log_weight <=
          design$invariant_control_log_tolerance
    )
    pass <-
      all(selected$relative_weight_ess >=
            gates$minimum_relative_weight_ess) &&
      all(selected$max_normalized_weight <=
            gates$maximum_normalized_weight) &&
      all(selected$split_log_ratio_difference <=
            gates$maximum_split_log_ratio_difference) &&
      diff(range(cohort$raw_log_ratio)) <=
        gates$maximum_cohort_chain_range &&
      nrow(selected_mcmc) > 0L &&
      all(selected_mcmc$split_rhat <= gates$maximum_split_rhat) &&
      all(selected_mcmc$mcmc_ess >= gates$minimum_relevant_mcmc_ess) &&
      all(selected_mcmc$tail_mcmc_ess >=
            gates$minimum_relevant_tail_ess) &&
      max(selected_mcse$forward_log_mcse_batch_ratio) <=
        gates$maximum_batch_mcse_ratio &&
      projection$projected_log_mcse <=
        gates$maximum_projected_115_log_mcse &&
      projection$maximum_projected_variance_share <=
        variance_share_limit + 1e-12 &&
      invariant_control_passed &&
      reference_cost_passed
    data.frame(
      endpoint_id = id,
      minimum_relative_weight_ess = min(selected$relative_weight_ess),
      maximum_normalized_weight = max(selected$max_normalized_weight),
      maximum_d2 = max(selected$d2),
      maximum_split_log_ratio_difference =
        max(selected$split_log_ratio_difference),
      cohort_chain_range = diff(range(cohort$raw_log_ratio)),
      maximum_relevant_split_rhat = max(selected_mcmc$split_rhat),
      minimum_relevant_mcmc_ess = min(selected_mcmc$mcmc_ess),
      minimum_relevant_tail_ess = min(selected_mcmc$tail_mcmc_ess),
      maximum_forward_batch_mcse_ratio =
        max(selected_mcse$forward_log_mcse_batch_ratio),
      projected_115_forward_log_mcse = projection$projected_log_mcse,
      maximum_forward_projected_variance_share =
        projection$maximum_projected_variance_share,
      maximum_projected_variance_share_limit = variance_share_limit,
      invariant_control_identity_passed = invariant_control_passed,
      forward_worst_observed_patient_envelope =
        projection$worst_observed_patient_envelope,
      forward_delta_log_bias_proxy = projection$delta_log_bias_proxy,
      maximum_reference_prediction_calls = maximum_reference_calls,
      reference_cost_passed = reference_cost_passed,
      passed = pass,
      stringsAsFactors = FALSE
    )
  }))
  rownames(selection) <- NULL
  endpoints$pilot_gate_passed <- selection$passed[
    match(endpoints$endpoint_id, selection$endpoint_id)
  ]
  endpoints$diagnostic_fallback <- !endpoints$pilot_gate_passed

  plan <- structure(list(
    schema_version = "sab_system_a_pure_message_endpoint_plan_v1",
    design_schema_version = design$schema_version,
    target_fingerprint = adapter$target_fingerprint,
    numerical_target = adapter$numerical_target,
    fixed_shared_parameters = TRUE,
    psi = canonical_anchor$psi,
    reference = reference,
    patient_ids = design$patient_ids,
    chain_roles = list(
      pilot = design$pilot_chains,
      forward = design$forward_chains,
      reverse = design$reverse_chains,
      bridge_reference = design$bridge_reference_chains,
      bridge_candidate = design$bridge_candidate_chains
    ),
    bank_contract = design$bank_contract,
    endpoints = endpoints,
    endpoint_eta = endpoint_eta,
    projection_patient_count = design$projection_patient_count,
    full_treatment_counts = design$full_treatment_counts,
    anchor = list(
      path = anchor_path, sha256 = anchor_sha256,
      eta = canonical_anchor$eta
    ),
    anchor_banks = list(
      paths = anchor_bank_paths, sha256 = anchor_bank_sha256
    ),
    pilot_selection = selection
  ), class = c("sab_system_a_pure_message_endpoint_plan", "list"))

  list(
    plan = plan,
    pilot_diagnostics = diagnostics,
    pilot_mcmc_diagnostics = mcmc,
    pilot_patient_mcse = patient_mcse,
    selection = selection,
    reference_bank_performance =
      .sab_pure_reference_performance(anchor_bank_artifacts),
    pilot_passed = all(selection$passed)
  )
}

.sab_pure_validate_plan <- function(plan, adapter, canonical_anchor, design) {
  expected <- .sab_pure_endpoints(adapter, canonical_anchor)
  endpoint_names <- c(
    "endpoint_id", "stage", "kind", "axis", "parameter", "sign", "step",
    "pilot_gate_passed", "diagnostic_fallback"
  )
  valid <- inherits(plan, "sab_system_a_pure_message_endpoint_plan") &&
    identical(plan$schema_version,
              "sab_system_a_pure_message_endpoint_plan_v1") &&
    identical(plan$design_schema_version, design$schema_version) &&
    identical(plan$target_fingerprint, adapter$target_fingerprint) &&
    isTRUE(plan$fixed_shared_parameters) &&
    identical(plan$psi, canonical_anchor$psi) &&
    identical(plan$reference, .sab_pure_reference(canonical_anchor)) &&
    identical(plan$patient_ids, design$patient_ids) &&
    identical(plan$chain_roles, list(
      pilot = design$pilot_chains,
      forward = design$forward_chains,
      reverse = design$reverse_chains,
      bridge_reference = design$bridge_reference_chains,
      bridge_candidate = design$bridge_candidate_chains
    )) &&
    identical(plan$bank_contract, design$bank_contract) &&
    identical(plan$anchor$eta, canonical_anchor$eta) &&
    is.data.frame(plan$endpoints) && nrow(plan$endpoints) == 4L &&
    identical(names(plan$endpoints), endpoint_names) &&
    identical(plan$endpoints[, setdiff(endpoint_names,
                                       c("pilot_gate_passed",
                                         "diagnostic_fallback"))],
              expected$endpoints[, setdiff(endpoint_names,
                                           c("pilot_gate_passed",
                                             "diagnostic_fallback"))]) &&
    is.matrix(plan$endpoint_eta) && is.numeric(plan$endpoint_eta) &&
    identical(dim(plan$endpoint_eta), c(4L, length(canonical_anchor$eta))) &&
    identical(rownames(plan$endpoint_eta), plan$endpoints$endpoint_id) &&
    identical(colnames(plan$endpoint_eta), names(canonical_anchor$eta)) &&
    isTRUE(all.equal(unname(plan$endpoint_eta),
                     unname(expected$endpoint_eta), tolerance = 0)) &&
    is.data.frame(plan$pilot_selection) &&
    identical(as.character(plan$pilot_selection$endpoint_id),
              as.character(plan$endpoints$endpoint_id)) &&
    identical(as.logical(plan$pilot_selection$passed),
              as.logical(plan$endpoints$pilot_gate_passed)) &&
    identical(plan$projection_patient_count, 115L) &&
    identical(plan$full_treatment_counts, design$full_treatment_counts)
  if (!isTRUE(valid)) {
    stop("Malformed or incompatible pure-SAEM endpoint plan.", call. = FALSE)
  }
  invisible(TRUE)
}

.sab_pure_unavailable_mcmc_rows <- function(endpoint_id, patient_id,
                                            quantities, bridge) {
  do.call(rbind, lapply(quantities, function(quantity) data.frame(
    endpoint_id = endpoint_id,
    patient_id = patient_id,
    quantity = quantity,
    required_for_gate = TRUE,
    split_rhat = NA_real_,
    rank_normalized_split_rhat = NA_real_,
    folded_split_rhat = NA_real_,
    mcmc_ess = NA_real_,
    tail_mcmc_ess = NA_real_,
    relative_mcmc_ess = NA_real_,
    chain_mean_range = NA_real_,
    diagnostic_available = FALSE,
    failure_code = bridge$failure_code,
    failure_message = bridge$failure_message,
    stringsAsFactors = FALSE
  )))
}

.sab_pure_cost_summary <- function(reference_artifacts, candidate_artifacts) {
  design <- sab_system_a_pure_message_design()
  summarize <- function(artifacts, role, endpoint_id) {
    rows <- do.call(rbind, lapply(artifacts, function(artifact) artifact$summary))
    statuses <- lapply(
      artifacts, .sab_pure_bank_budget_status, design = design
    )
    data.frame(
      role = role,
      endpoint_id = endpoint_id,
      bank_tasks = length(unique(rows$patient_id)),
      chains = nrow(rows),
      exact_prediction_calls = sum(rows$exact_prediction_calls),
      exact_ode_integrations = sum(rows$exact_ode_integrations),
      maximum_patient_target_prediction_calls = max(vapply(
        statuses, function(value) value$observed_prediction_calls, numeric(1L)
      )),
      maximum_patient_target_ode_integrations = max(vapply(
        statuses, function(value) value$observed_ode_integrations, numeric(1L)
      )),
      prediction_call_cap_per_patient_target =
        design$bank_contract$maximum_prediction_calls_per_patient_target,
      within_frozen_budget = all(vapply(
        statuses, function(value) value$within_budget, logical(1L)
      )),
      stringsAsFactors = FALSE
    )
  }
  rows <- list(summarize(
    reference_artifacts, "deployable_reference", NA_character_
  ))
  index <- 1L
  for (endpoint_id in names(candidate_artifacts)) {
    artifacts <- candidate_artifacts[[endpoint_id]]
    candidate <- do.call(rbind, lapply(artifacts, function(artifact) {
      artifact$summary
    }))
    required <- c("patient_id", "exact_prediction_calls",
                  "exact_ode_integrations")
    if (!is.data.frame(candidate) ||
        length(setdiff(required, names(candidate))) ||
        any(!is.finite(candidate$exact_prediction_calls)) ||
        any(!is.finite(candidate$exact_ode_integrations))) {
      stop("Candidate summaries omit required cost diagnostics.",
           call. = FALSE)
    }
    index <- index + 1L
    rows[[index]] <- summarize(
      artifacts, "validation_candidate", endpoint_id
    )
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

#' Assess the final pure-SAEM fixed-psi message experiment
#'
#' @return Held-out forward, independent reverse, bridge, mixing, projected
#'   cohort-error, and exact-cost diagnostics with a strict pass/fail status.
#' @export
sab_system_a_assess_pure_messages <- function(
    adapter, canonical_anchor, plan, plan_sha256, anchor_bank_artifacts,
    candidate_bank_artifacts, bridge_tolerance = 1e-10,
    bridge_max_iterations = 10000L) {
  .sab_pure_required_functions()
  design <- sab_system_a_pure_message_design()
  .sab_pure_validate_plan(plan, adapter, canonical_anchor, design)
  .sab_pure_validate_reference_artifacts(
    adapter, canonical_anchor, anchor_bank_artifacts, design
  )
  if (!is.character(plan_sha256) || length(plan_sha256) != 1L ||
      !grepl("^[0-9a-f]{64}$", plan_sha256) ||
      !is.list(candidate_bank_artifacts) ||
      !identical(names(candidate_bank_artifacts), plan$endpoints$endpoint_id)) {
    stop("Candidate-bank collection or endpoint-plan hash is malformed.",
         call. = FALSE)
  }
  density <- .sab_pure_density(adapter)
  reference <- plan$reference
  directional_rows <- bridge_rows <- patient_rows <- mcmc_rows <- list()
  directional_index <- bridge_index <- patient_index <- mcmc_index <- 0L

  for (endpoint_index in seq_len(nrow(plan$endpoints))) {
    endpoint <- plan$endpoints[endpoint_index, , drop = FALSE]
    endpoint_id <- endpoint$endpoint_id[[1L]]
    candidate_eta <- plan$endpoint_eta[endpoint_id, ]
    names(candidate_eta) <- colnames(plan$endpoint_eta)
    candidate <- list(kind = "population", eta = candidate_eta)
    candidates <- candidate_bank_artifacts[[endpoint_id]]
    if (!is.list(candidates) ||
        !identical(names(candidates), design$patient_ids)) {
      stop("Candidate banks are incomplete for ", endpoint_id, ".",
           call. = FALSE)
    }

    for (patient_id in design$patient_ids) {
      anchor_artifact <- anchor_bank_artifacts[[patient_id]]
      candidate_artifact <- candidates[[patient_id]]
      .sab_pure_validate_candidate_artifact(
        candidate_artifact, adapter, canonical_anchor, plan, plan_sha256,
        endpoint, candidate_eta, patient_id
      )
      context <- .sab_pure_patient_context(
        anchor_artifact, adapter, patient_id
      )
      evaluate_anchor <- function(chain_names, record_direction = FALSE) {
        result <- lapply(chain_names, function(chain_name) {
          chain <- anchor_artifact$chains[[chain_name]]
          population <- sab_evaluate_population_message(
            chain$draws, reference, candidate, context, density
          )
          diagnostic <- sab_raw_message_diagnostics(population$log_weights)
          if (record_direction) {
            directional_index <<- directional_index + 1L
            directional_rows[[directional_index]] <<- .sab_pure_diagnostic_row(
              endpoint, patient_id, context$treat_nelf, chain_name,
              "forward", diagnostic, population$log_weights
            )
          }
          list(chain = chain, population = population, diagnostic = diagnostic)
        })
        names(result) <- chain_names
        result
      }
      evaluate_candidate <- function(chain_names, record_direction = FALSE) {
        result <- lapply(chain_names, function(chain_name) {
        chain <- candidate_artifact$chains[[chain_name]]
        population <- sab_evaluate_population_message(
          chain$draws, candidate, reference, context, density
        )
        diagnostic <- sab_raw_message_diagnostics(population$log_weights)
        if (record_direction) {
          directional_index <<- directional_index + 1L
          directional_rows[[directional_index]] <<- .sab_pure_diagnostic_row(
            endpoint, patient_id, context$treat_nelf, chain_name,
            "reverse", diagnostic, population$log_weights
          )
        }
        list(chain = chain, population = population, diagnostic = diagnostic)
        })
        names(result) <- chain_names
        result
      }
      forward_evaluations <- evaluate_anchor(
        design$forward_chains, record_direction = TRUE
      )
      reverse_evaluations <- evaluate_candidate(
        design$reverse_chains, record_direction = TRUE
      )
      bridge_anchor_evaluations <- evaluate_anchor(
        design$bridge_reference_chains
      )
      bridge_candidate_evaluations <- evaluate_candidate(
        design$bridge_candidate_chains
      )

      forward_weights <- unlist(lapply(
        forward_evaluations, function(value) value$population$log_weights
      ), use.names = FALSE)
      reverse_weights <- unlist(lapply(
        reverse_evaluations, function(value) value$population$log_weights
      ), use.names = FALSE)
      pooled_bridge <- .sab_me_bridge(
        bridge_anchor_evaluations, bridge_candidate_evaluations,
        endpoint_id = endpoint_id, patient_id = patient_id,
        bridge_scope = "pooled", tolerance = bridge_tolerance,
        max_iterations = bridge_max_iterations
      )

      quantities <- c(
        list(
          forward_raw_weight_contribution =
            .sab_me_rescaled_exp_chains(lapply(
              forward_evaluations,
              function(value) value$population$log_weights
            )),
          reverse_raw_weight_contribution =
            .sab_me_rescaled_exp_chains(lapply(
              reverse_evaluations,
              function(value) value$population$log_weights
            ))
        ),
        .sab_pure_relevant_state_quantities(
          endpoint, context, forward_evaluations, "reference_forward"
        ),
        .sab_pure_relevant_state_quantities(
          endpoint, context, reverse_evaluations, "candidate_reverse"
        ),
        .sab_pure_relevant_state_quantities(
          endpoint, context, bridge_anchor_evaluations, "reference_bridge"
        ),
        .sab_pure_relevant_state_quantities(
          endpoint, context, bridge_candidate_evaluations, "candidate_bridge"
        )
      )
      required_names <- names(quantities)
      # All coordinates remain visible as secondary diagnostics, but only the
      # endpoint-relevant sufficient statistics decide this narrow gate.
      for (coordinate in adapter$coordinate_names$local) {
        quantities[[paste0("reference_x_all/", coordinate)]] <- lapply(
          forward_evaluations,
          function(value) value$chain$draws[, coordinate]
        )
        quantities[[paste0("candidate_x_all/", coordinate)]] <- lapply(
          reverse_evaluations,
          function(value) value$chain$draws[, coordinate]
        )
      }
      if (pooled_bridge$converged) {
        bridge_scores <- .sab_me_bridge_scores(
          bridge_anchor_evaluations, bridge_candidate_evaluations,
          pooled_bridge$log_ratio
        )
        quantities$reference_bridge_score <- bridge_scores$reference
        quantities$candidate_bridge_score <- bridge_scores$candidate
        required_names <- c(
          required_names, "reference_bridge_score", "candidate_bridge_score"
        )
      }
      rows <- .sab_pure_mcmc_rows(
        quantities, endpoint_id, patient_id, required_names
      )
      mcmc_index <- mcmc_index + 1L
      mcmc_rows[[mcmc_index]] <- rows
      if (!pooled_bridge$converged) {
        mcmc_index <- mcmc_index + 1L
        mcmc_rows[[mcmc_index]] <- .sab_pure_unavailable_mcmc_rows(
          endpoint_id, patient_id,
          c("reference_bridge_score", "candidate_bridge_score"),
          pooled_bridge
        )
      }

      pair_bridges <- numeric()
      pair_converged <- logical()
      for (anchor_name in names(bridge_anchor_evaluations)) {
        for (candidate_name in names(bridge_candidate_evaluations)) {
          if (pooled_bridge$converged) {
            bridge <- .sab_me_bridge(
              bridge_anchor_evaluations[anchor_name],
              bridge_candidate_evaluations[candidate_name],
              endpoint_id = endpoint_id, patient_id = patient_id,
              bridge_scope = "chain_pair", anchor_chain = anchor_name,
              candidate_chain = candidate_name,
              tolerance = bridge_tolerance,
              max_iterations = bridge_max_iterations
            )
          } else {
            bridge <- list(
              log_ratio = NA_real_, last_log_ratio = NA_real_,
              iterations = 0L, converged = FALSE,
              final_increment = NA_real_,
              failure_code = "not_attempted_after_pooled_failure",
              failure_message = paste0(
                "Chain-pair bridge was not attempted because the pooled ",
                "bridge did not converge."
              )
            )
          }
          pair_name <- paste(anchor_name, candidate_name, sep = "__")
          pair_bridges[[pair_name]] <- bridge$log_ratio
          pair_converged[[pair_name]] <- bridge$converged
          bridge_index <- bridge_index + 1L
          bridge_rows[[bridge_index]] <- data.frame(
            endpoint_id = endpoint_id,
            patient_id = patient_id,
            anchor_chain = anchor_name,
            candidate_chain = candidate_name,
            log_ratio = bridge$log_ratio,
            iterations = bridge$iterations,
            converged = bridge$converged,
            final_increment = bridge$final_increment,
            last_log_ratio = bridge$last_log_ratio,
            failure_code = bridge$failure_code,
            failure_message = bridge$failure_message,
            stringsAsFactors = FALSE
          )
        }
      }
      if (all(pair_converged)) {
        bridge_matrix <- matrix(
          pair_bridges, nrow = length(bridge_anchor_evaluations),
          ncol = length(bridge_candidate_evaluations), byrow = TRUE
        )
        bridge_instability <- sqrt(
          stats::var(rowMeans(bridge_matrix)) / nrow(bridge_matrix) +
            stats::var(colMeans(bridge_matrix)) / ncol(bridge_matrix)
        )
        bridge_pair_range <- diff(range(pair_bridges))
      } else {
        bridge_instability <- NA_real_
        bridge_pair_range <- NA_real_
      }
      forward_mcse <- .sab_pure_patient_log_mcse(forward_evaluations)
      reverse_mcse <- .sab_pure_patient_log_mcse(reverse_evaluations)
      forward_ratio <- sab_log_mean_exp(forward_weights)
      reverse_ratio <- -sab_log_mean_exp(reverse_weights)
      bridge_ratio <- pooled_bridge$log_ratio
      bridge_available <- pooled_bridge$converged && all(pair_converged) &&
        is.finite(bridge_instability)
      invariant_control_required <-
        endpoint$parameter[[1L]] == "beta_nelf" && context$treat_nelf == 0L
      invariant_control_identity_passed <- !invariant_control_required ||
        (
          max(abs(forward_weights)) <=
            design$invariant_control_log_tolerance &&
          max(abs(reverse_weights)) <=
            design$invariant_control_log_tolerance &&
          bridge_available &&
          abs(bridge_ratio) <= design$invariant_control_log_tolerance &&
          all(abs(pair_bridges) <= design$invariant_control_log_tolerance)
        )
      if (bridge_available) {
        forward_tolerance <- max(
          design$evaluation_gates$minimum_agreement_tolerance,
          2 * sqrt(forward_mcse$conservative^2 + bridge_instability^2)
        )
        reverse_tolerance <- max(
          design$evaluation_gates$minimum_agreement_tolerance,
          2 * sqrt(reverse_mcse$conservative^2 + bridge_instability^2)
        )
      } else {
        forward_tolerance <- reverse_tolerance <- NA_real_
      }
      forward_reverse_tolerance <- max(
        design$evaluation_gates$minimum_agreement_tolerance,
        2 * sqrt(forward_mcse$conservative^2 +
                   reverse_mcse$conservative^2)
      )
      patient_index <- patient_index + 1L
      patient_rows[[patient_index]] <- data.frame(
        endpoint_id = endpoint_id,
        axis = endpoint$axis,
        parameter = endpoint$parameter,
        sign = endpoint$sign,
        step = endpoint$step,
        patient_id = patient_id,
        treatment = context$treat_nelf,
        invariant_control_required = invariant_control_required,
        invariant_control_identity_passed =
          invariant_control_identity_passed,
        forward_log_ratio = forward_ratio,
        reverse_log_ratio = reverse_ratio,
        bridge_log_ratio = bridge_ratio,
        bridge_last_log_ratio = pooled_bridge$last_log_ratio,
        bridge_converged = pooled_bridge$converged,
        bridge_iterations = pooled_bridge$iterations,
        bridge_final_increment = pooled_bridge$final_increment,
        bridge_failure_code = pooled_bridge$failure_code,
        bridge_failure_message = pooled_bridge$failure_message,
        all_bridge_pairs_converged = all(pair_converged),
        bridge_pair_failure_count = sum(!pair_converged),
        bridge_chain_pair_instability = bridge_instability,
        bridge_pair_range = bridge_pair_range,
        forward_log_mcse = forward_mcse$conservative,
        forward_log_mcse_batch_5 = forward_mcse$batch_5,
        forward_log_mcse_batch_10 = forward_mcse$batch_10,
        forward_log_mcse_batch_ratio = forward_mcse$batch_ratio,
        reverse_log_mcse = reverse_mcse$conservative,
        reverse_log_mcse_batch_5 = reverse_mcse$batch_5,
        reverse_log_mcse_batch_10 = reverse_mcse$batch_10,
        reverse_log_mcse_batch_ratio = reverse_mcse$batch_ratio,
        forward_reverse_difference = abs(forward_ratio - reverse_ratio),
        forward_reverse_tolerance = forward_reverse_tolerance,
        forward_reverse_passed =
          abs(forward_ratio - reverse_ratio) <= forward_reverse_tolerance,
        forward_bridge_difference = abs(forward_ratio - bridge_ratio),
        reverse_bridge_difference = abs(reverse_ratio - bridge_ratio),
        forward_bridge_tolerance = forward_tolerance,
        reverse_bridge_tolerance = reverse_tolerance,
        forward_bridge_passed = bridge_available &&
          abs(forward_ratio - bridge_ratio) <= forward_tolerance,
        reverse_bridge_passed = bridge_available &&
          abs(reverse_ratio - bridge_ratio) <= reverse_tolerance,
        stringsAsFactors = FALSE
      )
    }
  }

  directional <- do.call(rbind, directional_rows)
  bridges <- do.call(rbind, bridge_rows)
  patients <- do.call(rbind, patient_rows)
  mcmc <- do.call(rbind, mcmc_rows)
  rownames(directional) <- rownames(bridges) <- rownames(patients) <-
    rownames(mcmc) <- NULL
  gates <- design$evaluation_gates
  cost_summary <- .sab_pure_cost_summary(
    anchor_bank_artifacts, candidate_bank_artifacts
  )

  endpoint_summary <- do.call(rbind, lapply(
    plan$endpoints$endpoint_id, function(endpoint_id) {
      endpoint <- plan$endpoints[
        plan$endpoints$endpoint_id == endpoint_id, , drop = FALSE
      ]
      diagnostic <- directional[
        directional$endpoint_id == endpoint_id, , drop = FALSE
      ]
      patient <- patients[patients$endpoint_id == endpoint_id, , drop = FALSE]
      bridge <- bridges[bridges$endpoint_id == endpoint_id, , drop = FALSE]
      endpoint_mcmc <- mcmc[
        mcmc$endpoint_id == endpoint_id & mcmc$required_for_gate,
        , drop = FALSE
      ]
      forward <- diagnostic[diagnostic$direction == "forward", , drop = FALSE]
      reverse <- diagnostic[diagnostic$direction == "reverse", , drop = FALSE]
      forward_totals <- stats::aggregate(
        oriented_log_ratio ~ chain, forward, sum
      )$oriented_log_ratio
      reverse_totals <- stats::aggregate(
        oriented_log_ratio ~ chain, reverse, sum
      )$oriented_log_ratio
      bridge_groups <- split(
        bridge,
        interaction(bridge$anchor_chain, bridge$candidate_chain, drop = TRUE)
      )
      bridge_totals <- vapply(bridge_groups, function(value) {
        if (nrow(value) != nrow(patient) || !all(value$converged) ||
            any(!is.finite(value$log_ratio))) return(NA_real_)
        sum(value$log_ratio)
      }, numeric(1L))
      bridge_converged <- nrow(bridge) > 0L &&
        all(patient$bridge_converged) &&
        all(patient$all_bridge_pairs_converged) &&
        all(is.finite(bridge_totals))
      projected_forward <- .sab_pure_projection_summary(
        patient$forward_log_mcse, patient$treatment, design
      )
      projected_reverse <- .sab_pure_projection_summary(
        patient$reverse_log_mcse, patient$treatment, design
      )
      variance_share_limit <- .sab_pure_variance_share_limit(endpoint, gates)
      forward_total <- sum(patient$forward_log_ratio)
      reverse_total <- sum(patient$reverse_log_ratio)
      bridge_total <- if (bridge_converged) {
        sum(patient$bridge_log_ratio)
      } else {
        NA_real_
      }
      forward_cohort_mcse <- sqrt(sum(patient$forward_log_mcse^2))
      reverse_cohort_mcse <- sqrt(sum(patient$reverse_log_mcse^2))
      bridge_cohort_instability <- if (bridge_converged) {
        stats::sd(bridge_totals)
      } else {
        NA_real_
      }
      forward_reverse_tolerance <- max(
        gates$minimum_cohort_agreement_tolerance,
        2 * sqrt(forward_cohort_mcse^2 + reverse_cohort_mcse^2)
      )
      forward_bridge_tolerance <- if (bridge_converged) {
        max(
          gates$minimum_cohort_agreement_tolerance,
          2 * sqrt(forward_cohort_mcse^2 + bridge_cohort_instability^2)
        )
      } else {
        NA_real_
      }
      reverse_bridge_tolerance <- if (bridge_converged) {
        max(
          gates$minimum_cohort_agreement_tolerance,
          2 * sqrt(reverse_cohort_mcse^2 + bridge_cohort_instability^2)
        )
      } else {
        NA_real_
      }
      overlap_passed <-
        min(diagnostic$relative_weight_ess) >=
          gates$minimum_relative_weight_ess &&
        max(diagnostic$max_normalized_weight) <=
          gates$maximum_normalized_weight &&
        max(diagnostic$split_log_ratio_difference) <=
          gates$maximum_split_log_ratio_difference
      replication_passed <- bridge_converged &&
        diff(range(forward_totals)) <= gates$maximum_cohort_chain_range &&
        diff(range(reverse_totals)) <= gates$maximum_cohort_chain_range &&
        diff(range(bridge_totals)) <= gates$maximum_cohort_chain_range
      agreement_passed <- bridge_converged &&
        all(patient$forward_bridge_passed) &&
        all(patient$reverse_bridge_passed) &&
        all(patient$forward_reverse_passed) &&
        abs(forward_total - reverse_total) <= forward_reverse_tolerance &&
        abs(forward_total - bridge_total) <= forward_bridge_tolerance &&
        abs(reverse_total - bridge_total) <= reverse_bridge_tolerance
      mcmc_passed <- nrow(endpoint_mcmc) > 0L &&
        all(endpoint_mcmc$diagnostic_available) &&
        all(is.finite(endpoint_mcmc$split_rhat)) &&
        max(endpoint_mcmc$split_rhat) <= gates$maximum_split_rhat &&
        min(endpoint_mcmc$mcmc_ess) >= gates$minimum_relevant_mcmc_ess &&
        all(is.finite(endpoint_mcmc$tail_mcmc_ess)) &&
        min(endpoint_mcmc$tail_mcmc_ess) >=
          gates$minimum_relevant_tail_ess
      projected_mcse_passed <-
        projected_forward$projected_log_mcse <=
          gates$maximum_projected_115_log_mcse &&
        projected_reverse$projected_log_mcse <=
          gates$maximum_projected_115_log_mcse &&
        projected_forward$maximum_projected_variance_share <=
          variance_share_limit + 1e-12 &&
        projected_reverse$maximum_projected_variance_share <=
          variance_share_limit + 1e-12 &&
        max(patient$forward_log_mcse_batch_ratio) <=
          gates$maximum_batch_mcse_ratio &&
        max(patient$reverse_log_mcse_batch_ratio) <=
          gates$maximum_batch_mcse_ratio
      invariant_control_identity_passed <-
        all(patient$invariant_control_identity_passed)
      endpoint_cost <- cost_summary[
        is.na(cost_summary$endpoint_id) |
          cost_summary$endpoint_id == endpoint_id, , drop = FALSE
      ]
      cost_passed <- nrow(endpoint_cost) == 2L &&
        all(endpoint_cost$within_frozen_budget)
      validated <- endpoint$pilot_gate_passed && overlap_passed &&
        replication_passed && agreement_passed && mcmc_passed &&
        projected_mcse_passed && invariant_control_identity_passed &&
        cost_passed
      data.frame(
        endpoint_id = endpoint_id,
        axis = endpoint$axis,
        parameter = endpoint$parameter,
        sign = endpoint$sign,
        step = endpoint$step,
        pilot_gate_passed = endpoint$pilot_gate_passed,
        forward_log_ratio = forward_total,
        reverse_log_ratio = reverse_total,
        bridge_log_ratio = bridge_total,
        forward_reverse_difference = abs(forward_total - reverse_total),
        forward_reverse_tolerance = forward_reverse_tolerance,
        forward_bridge_difference = abs(forward_total - bridge_total),
        forward_bridge_tolerance = forward_bridge_tolerance,
        reverse_bridge_difference = abs(reverse_total - bridge_total),
        reverse_bridge_tolerance = reverse_bridge_tolerance,
        minimum_relative_weight_ess = min(diagnostic$relative_weight_ess),
        maximum_normalized_weight = max(diagnostic$max_normalized_weight),
        maximum_d2 = max(diagnostic$d2),
        maximum_split_log_ratio_difference =
          max(diagnostic$split_log_ratio_difference),
        forward_chain_range = diff(range(forward_totals)),
        reverse_chain_range = diff(range(reverse_totals)),
        bridge_pair_range = if (bridge_converged) {
          diff(range(bridge_totals))
        } else {
          NA_real_
        },
        projected_115_forward_log_mcse =
          projected_forward$projected_log_mcse,
        projected_115_reverse_log_mcse =
          projected_reverse$projected_log_mcse,
        maximum_forward_projected_variance_share =
          projected_forward$maximum_projected_variance_share,
        maximum_reverse_projected_variance_share =
          projected_reverse$maximum_projected_variance_share,
        maximum_projected_variance_share_limit = variance_share_limit,
        invariant_control_identity_passed =
          invariant_control_identity_passed,
        forward_worst_observed_patient_envelope =
          projected_forward$worst_observed_patient_envelope,
        reverse_worst_observed_patient_envelope =
          projected_reverse$worst_observed_patient_envelope,
        forward_delta_log_bias_proxy =
          projected_forward$delta_log_bias_proxy,
        reverse_delta_log_bias_proxy =
          projected_reverse$delta_log_bias_proxy,
        maximum_forward_batch_mcse_ratio =
          max(patient$forward_log_mcse_batch_ratio),
        maximum_reverse_batch_mcse_ratio =
          max(patient$reverse_log_mcse_batch_ratio),
        maximum_relevant_split_rhat = if (nrow(endpoint_mcmc)) {
          max(endpoint_mcmc$split_rhat, na.rm = TRUE)
        } else {
          NA_real_
        },
        minimum_relevant_mcmc_ess = if (nrow(endpoint_mcmc)) {
          min(endpoint_mcmc$mcmc_ess, na.rm = TRUE)
        } else {
          NA_real_
        },
        minimum_relevant_tail_ess = if (nrow(endpoint_mcmc)) {
          min(endpoint_mcmc$tail_mcmc_ess, na.rm = TRUE)
        } else {
          NA_real_
        },
        pooled_bridge_failures = sum(!patient$bridge_converged),
        chain_pair_bridge_failures = sum(!bridge$converged),
        overlap_passed = overlap_passed,
        replication_passed = replication_passed,
        agreement_passed = agreement_passed,
        mcmc_diagnostics_passed = mcmc_passed,
        projected_mcse_passed = projected_mcse_passed,
        cost_passed = cost_passed,
        validated = validated,
        stringsAsFactors = FALSE
      )
    }
  ))
  rownames(endpoint_summary) <- NULL
  passed <- all(endpoint_summary$validated)
  list(
    endpoint_summary = endpoint_summary,
    patient_summary = patients,
    directional_diagnostics = directional,
    mcmc_diagnostics = mcmc,
    bridge_chain_pairs = bridges,
    reference_bank_performance =
      .sab_pure_reference_performance(anchor_bank_artifacts),
    cost_summary = cost_summary,
    falsification_passed = passed,
    status = if (passed) {
      "pure_saem_falsification_passed"
    } else {
      "pure_saem_falsification_failed"
    },
    scope = paste0(
      "Twelve-patient fixed-psi local population-message test at one ",
      "predeclared SAEM branch; not a 115-patient posterior or dynamic-psi ",
      "validation."
    )
  )
}
