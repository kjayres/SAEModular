# Minimal fixed-psi raw-message validation for the twelve audited System A
# patients.  The population density is always evaluated through a callback;
# none of the message calculations below assumes a Gaussian population law.

sab_system_a_message_design <- function() {
  list(
    schema_version = "sab_system_a_message_design_v1",
    patient_ids = c(
      "3", "6", "16", "20", "55", "71", "74", "88", "105",
      "111", "117", "122"
    ),
    pilot_chains = c("chain_01", "chain_02"),
    evaluation_chains = c("chain_03", "chain_04"),
    reference = list(
      kind = "equal_defensive_mixture",
      weights = c(saem = 0.5, priorcentral = 0.5)
    ),
    core_endpoints = data.frame(
      endpoint_id = c("saem_component", "priorcentral_component"),
      stage = rep("core_component", 2L),
      kind = rep("population_component", 2L),
      axis = rep("reference_component", 2L),
      parameter = rep(NA_character_, 2L),
      sign = rep(0L, 2L),
      step = rep(0, 2L),
      stringsAsFactors = FALSE
    ),
    # These local-stress directions are intentionally not part of the first
    # candidate-bank run.  They are activated only after both separated
    # reference components pass raw/bridge validation.
    deferred_local_stress = list(
      beta_nelf = c(0.1, 0.2, 0.4, 0.8),
      log_omega_lambda = c(0.04375, 0.0875, 0.175, 0.35)
    ),
    pilot_gates = list(
      minimum_relative_weight_ess = 0.20,
      maximum_normalized_weight = 0.05,
      maximum_split_log_ratio_difference = 0.25,
      maximum_cohort_chain_range = 0.50,
      maximum_component_probability_chain_range = 0.10,
      maximum_mixture_identity_error = 1e-12,
      maximum_split_rhat = 1.01,
      minimum_mcmc_ess = 400
    ),
    evaluation_gates = list(
      minimum_relative_weight_ess = 0.10,
      maximum_normalized_weight = 0.05,
      maximum_cohort_chain_range = 0.50,
      minimum_agreement_tolerance = 0.25,
      maximum_component_probability_chain_range = 0.10,
      maximum_mixture_identity_error = 1e-12,
      maximum_split_rhat = 1.01,
      minimum_mcmc_ess = 400
    )
  )
}

sab_system_a_plan_message_endpoints <- function(
    adapter, canonical_anchor, anchor_bank_artifacts,
    anchor_path, anchor_sha256, anchor_bank_paths, anchor_bank_sha256) {
  design <- sab_system_a_message_design()
  .sab_me_validate_inputs(
    adapter, canonical_anchor, anchor_bank_artifacts, design,
    expected_role = "defensive_reference",
    expected_reference = .sab_me_reference(canonical_anchor, adapter)
  )
  if (!identical(names(anchor_bank_paths), design$patient_ids) ||
      !identical(names(anchor_bank_sha256), design$patient_ids)) {
    stop("Anchor-bank paths and hashes must be named in audited patient order.",
         call. = FALSE)
  }

  density <- .sab_me_system_a_density(adapter)
  reference <- .sab_me_reference(canonical_anchor, adapter)
  endpoint_eta <- rbind(
    saem_component = canonical_anchor$eta,
    priorcentral_component = adapter$prior_reference$eta
  )
  colnames(endpoint_eta) <- names(canonical_anchor$eta)
  pilot_rows <- list()
  row_index <- 0L
  for (endpoint_id in rownames(endpoint_eta)) {
    candidate <- list(kind = "population", eta = endpoint_eta[endpoint_id, ])
    names(candidate$eta) <- colnames(endpoint_eta)
    endpoint <- design$core_endpoints[
      design$core_endpoints$endpoint_id == endpoint_id, , drop = FALSE
    ]
    for (patient_id in design$patient_ids) {
      artifact <- anchor_bank_artifacts[[patient_id]]
      context <- .sab_me_patient_context(artifact, patient_id)
      for (chain_name in design$pilot_chains) {
        chain <- artifact$chains[[chain_name]]
        evaluated <- sab_evaluate_population_message(
          chain$draws,
          anchor_parameter = reference,
          candidate_parameter = candidate,
          patient_context = context,
          log_population_density = density
        )
        diagnostic <- sab_raw_message_diagnostics(evaluated$log_weights)
        row_index <- row_index + 1L
        pilot_rows[[row_index]] <- .sab_me_diagnostic_row(
          endpoint$axis, endpoint$parameter, endpoint$step, endpoint$sign,
          patient_id, chain_name, "forward", diagnostic,
          endpoint_id = endpoint_id
        )
      }
    }
  }
  pilot <- do.call(rbind, pilot_rows)
  rownames(pilot) <- NULL
  identity <- .sab_me_component_identity(
    adapter, canonical_anchor, anchor_bank_artifacts,
    design$pilot_chains, density
  )
  pilot_gate <- .sab_me_component_pilot_gate(pilot, identity, design)
  endpoints <- design$core_endpoints
  endpoints$pilot_gate_passed <- vapply(
    endpoints$endpoint_id,
    function(endpoint_id) pilot_gate$passed[pilot_gate$endpoint_id == endpoint_id],
    logical(1L)
  )
  endpoints$diagnostic_fallback <- !endpoints$pilot_gate_passed

  plan <- structure(list(
    schema_version = "sab_system_a_message_endpoint_plan_v1",
    target_fingerprint = adapter$target_fingerprint,
    numerical_target = adapter$numerical_target,
    fixed_shared_parameters = TRUE,
    psi = canonical_anchor$psi,
    reference = reference,
    patient_ids = design$patient_ids,
    design = design,
    anchor = list(
      path = normalizePath(anchor_path, mustWork = TRUE),
      sha256 = anchor_sha256,
      eta = canonical_anchor$eta,
      psi = canonical_anchor$psi
    ),
    anchor_banks = list(
      paths = anchor_bank_paths,
      sha256 = anchor_bank_sha256
    ),
    selection = pilot_gate,
    endpoints = endpoints,
    endpoint_eta = endpoint_eta,
    pilot_diagnostics = pilot,
    pilot_component_identity = identity
  ), class = c("sab_system_a_message_endpoint_plan", "list"))
  performance <- .sab_me_reference_bank_performance(anchor_bank_artifacts)
  list(
    plan = plan, pilot_diagnostics = pilot, selection = pilot_gate,
    component_identity = identity, reference_bank_performance = performance
  )
}

sab_system_a_assess_messages <- function(
    adapter, canonical_anchor, plan, plan_sha256, anchor_bank_artifacts,
    candidate_bank_artifacts) {
  design <- sab_system_a_message_design()
  .sab_me_validate_plan(plan, adapter, canonical_anchor, design)
  .sab_me_validate_inputs(
    adapter, canonical_anchor, anchor_bank_artifacts, design,
    expected_role = "defensive_reference",
    expected_reference = plan$reference
  )
  density <- .sab_me_system_a_density(adapter)
  reference <- plan$reference
  directional_rows <- list()
  bridge_rows <- list()
  patient_rows <- list()
  mcmc_rows <- list()
  directional_index <- bridge_index <- patient_index <- 0L
  mcmc_index <- 0L

  for (endpoint_index in seq_len(nrow(plan$endpoints))) {
    endpoint <- plan$endpoints[endpoint_index, , drop = FALSE]
    endpoint_id <- endpoint$endpoint_id[[1L]]
    candidate_eta <- plan$endpoint_eta[endpoint_id, ]
    names(candidate_eta) <- colnames(plan$endpoint_eta)
    candidate <- list(kind = "population", eta = candidate_eta)
    endpoint_candidates <- candidate_bank_artifacts[[endpoint_id]]
    if (!is.list(endpoint_candidates) ||
        !identical(names(endpoint_candidates), design$patient_ids)) {
      stop("Candidate banks are incomplete for endpoint ", endpoint_id, ".",
           call. = FALSE)
    }

    for (patient_id in design$patient_ids) {
      anchor_artifact <- anchor_bank_artifacts[[patient_id]]
      candidate_artifact <- endpoint_candidates[[patient_id]]
      .sab_me_validate_candidate(
        candidate_artifact, adapter, canonical_anchor, plan, endpoint,
        plan_sha256, patient_id, candidate_eta
      )
      context <- .sab_me_patient_context(anchor_artifact, patient_id)

      anchor_evaluations <- lapply(design$evaluation_chains, function(chain_name) {
        chain <- anchor_artifact$chains[[chain_name]]
        population <- sab_evaluate_population_message(
          chain$draws, reference, candidate, context, density
        )
        diagnostic <- sab_raw_message_diagnostics(population$log_weights)
        directional_index <<- directional_index + 1L
        directional_rows[[directional_index]] <<- .sab_me_diagnostic_row(
          endpoint$axis, endpoint$parameter, endpoint$step, endpoint$sign,
          patient_id, chain_name, "forward", diagnostic,
          endpoint_id = endpoint_id
        )
        list(chain = chain, population = population, diagnostic = diagnostic)
      })
      names(anchor_evaluations) <- design$evaluation_chains

      candidate_chain_names <- names(candidate_artifact$chains)
      candidate_evaluations <- lapply(candidate_chain_names, function(chain_name) {
        chain <- candidate_artifact$chains[[chain_name]]
        # Swap anchor/candidate arguments so the returned weights estimate
        # Z_anchor/Z_candidate on candidate-target draws.
        population <- sab_evaluate_population_message(
          chain$draws, candidate, reference, context, density
        )
        diagnostic <- sab_raw_message_diagnostics(population$log_weights)
        directional_index <<- directional_index + 1L
        row <- .sab_me_diagnostic_row(
          endpoint$axis, endpoint$parameter, endpoint$step, endpoint$sign,
          patient_id, chain_name, "reverse", diagnostic,
          endpoint_id = endpoint_id
        )
        row$oriented_log_ratio <- -row$raw_log_ratio
        directional_rows[[directional_index]] <<- row
        list(chain = chain, population = population, diagnostic = diagnostic)
      })
      names(candidate_evaluations) <- candidate_chain_names

      # The message estimators average exp(log weight), so their MCMC ESS must
      # diagnose those raw contributions rather than the log weights.  A
      # single common rescaling per direction prevents overflow and preserves
      # autocorrelation and split-Rhat.
      mcmc_quantities <- list(
        forward_raw_weight_contribution = .sab_me_rescaled_exp_chains(lapply(
          anchor_evaluations, function(value) value$population$log_weights
        )),
        reverse_raw_weight_contribution = .sab_me_rescaled_exp_chains(lapply(
          candidate_evaluations, function(value) value$population$log_weights
        ))
      )
      for (coordinate in adapter$coordinate_names$local) {
        mcmc_quantities[[paste0("reference_x/", coordinate)]] <- lapply(
          anchor_evaluations,
          function(value) value$chain$draws[, coordinate]
        )
        mcmc_quantities[[paste0("candidate_x/", coordinate)]] <- lapply(
          candidate_evaluations,
          function(value) value$chain$draws[, coordinate]
        )
      }

      forward_weights <- unlist(lapply(
        anchor_evaluations, function(value) value$population$log_weights
      ), use.names = FALSE)
      reverse_weights <- unlist(lapply(
        candidate_evaluations, function(value) value$population$log_weights
      ), use.names = FALSE)
      pooled_bridge <- .sab_me_bridge(
        anchor_evaluations, candidate_evaluations
      )
      bridge_scores <- .sab_me_bridge_scores(
        anchor_evaluations, candidate_evaluations, pooled_bridge$log_ratio
      )
      mcmc_quantities$reference_bridge_score <- bridge_scores$reference
      mcmc_quantities$candidate_bridge_score <- bridge_scores$candidate
      for (quantity in names(mcmc_quantities)) {
        diagnostic <- sab_mcmc_scalar_diagnostics(mcmc_quantities[[quantity]])
        mcmc_index <- mcmc_index + 1L
        mcmc_rows[[mcmc_index]] <- data.frame(
          endpoint_id = endpoint_id, patient_id = patient_id,
          quantity = quantity, split_rhat = diagnostic$split_rhat,
          mcmc_ess = diagnostic$mcmc_ess,
          relative_mcmc_ess = diagnostic$relative_mcmc_ess,
          chain_mean_range = diagnostic$chain_mean_range,
          stringsAsFactors = FALSE
        )
      }
      pair_bridges <- numeric()
      for (anchor_name in names(anchor_evaluations)) {
        for (candidate_name in names(candidate_evaluations)) {
          value <- .sab_me_bridge(
            anchor_evaluations[anchor_name],
            candidate_evaluations[candidate_name]
          )
          pair_name <- paste(anchor_name, candidate_name, sep = "__")
          pair_bridges[[pair_name]] <- value$log_ratio
          bridge_index <- bridge_index + 1L
          bridge_rows[[bridge_index]] <- data.frame(
            endpoint_id = endpoint_id, patient_id = patient_id,
            anchor_chain = anchor_name, candidate_chain = candidate_name,
            log_ratio = value$log_ratio, iterations = value$iterations,
            stringsAsFactors = FALSE
          )
        }
      }
      bridge_matrix <- matrix(
        pair_bridges,
        nrow = length(anchor_evaluations),
        ncol = length(candidate_evaluations),
        byrow = TRUE
      )
      # These crossed chain-pair estimates share draws.  This is an empirical
      # instability scale, not an independent replicate MCSE.
      bridge_chain_pair_instability <- sqrt(
        stats::var(rowMeans(bridge_matrix)) / nrow(bridge_matrix) +
          stats::var(colMeans(bridge_matrix)) / ncol(bridge_matrix)
      )
      forward_chain <- vapply(
        anchor_evaluations,
        function(value) value$diagnostic$log_ratio,
        numeric(1L)
      )
      reverse_chain <- -vapply(
        candidate_evaluations,
        function(value) value$diagnostic$log_ratio,
        numeric(1L)
      )
      forward_mcse <- sqrt(mean(vapply(
        anchor_evaluations,
        function(value) value$diagnostic$batch$log_scale_mcse^2,
        numeric(1L)
      )) / length(anchor_evaluations))
      reverse_mcse <- sqrt(mean(vapply(
        candidate_evaluations,
        function(value) value$diagnostic$batch$log_scale_mcse^2,
        numeric(1L)
      )) / length(candidate_evaluations))
      bridge_pair_sd <- stats::sd(pair_bridges)
      forward_ratio <- sab_log_mean_exp(forward_weights)
      reverse_ratio <- -sab_log_mean_exp(reverse_weights)
      bridge_ratio <- pooled_bridge$log_ratio
      forward_tolerance <- max(
        design$evaluation_gates$minimum_agreement_tolerance,
        2 * sqrt(forward_mcse^2 + bridge_chain_pair_instability^2)
      )
      reverse_tolerance <- max(
        design$evaluation_gates$minimum_agreement_tolerance,
        2 * sqrt(reverse_mcse^2 + bridge_chain_pair_instability^2)
      )
      patient_index <- patient_index + 1L
      patient_rows[[patient_index]] <- data.frame(
        endpoint_id = endpoint_id,
        axis = endpoint$axis,
        patient_id = patient_id,
        treatment = context$treat_nelf,
        forward_log_ratio = forward_ratio,
        reverse_log_ratio = reverse_ratio,
        bridge_log_ratio = bridge_ratio,
        forward_bridge_difference = abs(forward_ratio - bridge_ratio),
        reverse_bridge_difference = abs(reverse_ratio - bridge_ratio),
        forward_bridge_tolerance = forward_tolerance,
        reverse_bridge_tolerance = reverse_tolerance,
        forward_bridge_passed =
          abs(forward_ratio - bridge_ratio) <= forward_tolerance,
        reverse_bridge_passed =
          abs(reverse_ratio - bridge_ratio) <= reverse_tolerance,
        forward_log_mcse = forward_mcse,
        reverse_log_mcse = reverse_mcse,
        bridge_pair_sd = bridge_pair_sd,
        bridge_chain_pair_instability = bridge_chain_pair_instability,
        forward_chain_range = diff(range(forward_chain)),
        reverse_chain_range = diff(range(reverse_chain)),
        bridge_pair_range = diff(range(pair_bridges)),
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
  identity <- .sab_me_component_identity(
    adapter, canonical_anchor, anchor_bank_artifacts,
    design$evaluation_chains, density
  )
  bridge_identity <- .sab_me_bridge_identity(patients)
  endpoints <- .sab_me_endpoint_summary(
    plan, directional, bridges, patients, identity, bridge_identity,
    mcmc, design
  )
  calibration_passed <- all(endpoints$validated) &&
    all(bridge_identity$passed)
  list(
    endpoint_summary = endpoints,
    patient_summary = patients,
    component_identity = identity,
    bridge_component_identity = bridge_identity,
    mcmc_diagnostics = mcmc,
    reference_bank_performance =
      .sab_me_reference_bank_performance(anchor_bank_artifacts),
    directional_diagnostics = directional,
    bridge_chain_pairs = bridges,
    calibration_passed = calibration_passed,
    status = if (calibration_passed) {
      "calibration_passed_ready_for_local_stress"
    } else {
      "inconclusive_calibration_or_mixing_failure"
    }
  )
}

.sab_me_reference_bank_performance <- function(artifacts) {
  rows <- do.call(rbind, lapply(artifacts, function(artifact) {
    artifact$summary
  }))
  rownames(rows) <- NULL
  required <- c(
    "patient_id", "chain", "sampling_acceptance", "min_ess_x",
    "min_ess_x_squared", "min_ess_x_per_ode", "exact_prediction_calls",
    "exact_ode_integrations", "saem_component_proposals",
    "saem_component_proposal_fraction", "saem_component_acceptance",
    "saem_component_prediction_failure_rate", "saem_component_invalid_rate",
    "saem_component_prediction_failures",
    "saem_component_nonfinite_loglik", "priorcentral_component_proposals",
    "priorcentral_component_proposal_fraction",
    "priorcentral_component_acceptance",
    "priorcentral_component_prediction_failure_rate",
    "priorcentral_component_invalid_rate",
    "priorcentral_component_prediction_failures",
    "priorcentral_component_nonfinite_loglik"
  )
  if (length(setdiff(required, names(rows)))) {
    stop("Reference bank summaries omit component performance fields.",
         call. = FALSE)
  }
  rows
}

.sab_me_system_a_density <- function(adapter) {
  force(adapter)
  function(x, parameter, patient_context) {
    if (!is.list(parameter) || !is.character(parameter$kind) ||
        length(parameter$kind) != 1L) {
      stop("Population-density parameter object is malformed.", call. = FALSE)
    }
    if (parameter$kind == "population") {
      return(adapter$log_population_density(
        patient_context$patient_id, x, parameter$eta
      ))
    }
    if (parameter$kind != "equal_defensive_mixture" ||
        !identical(parameter$weights,
                   c(saem = 0.5, priorcentral = 0.5)) ||
        !is.list(parameter$eta_components) ||
        !identical(names(parameter$eta_components),
                   c("saem", "priorcentral"))) {
      stop("Unknown or malformed population reference.", call. = FALSE)
    }
    terms <- c(
      saem = log(0.5) + adapter$log_population_density(
        patient_context$patient_id, x, parameter$eta_components$saem
      ),
      priorcentral = log(0.5) + adapter$log_population_density(
        patient_context$patient_id, x,
        parameter$eta_components$priorcentral
      )
    )
    largest <- max(terms)
    if (largest == -Inf) return(-Inf)
    largest + log(sum(exp(terms - largest)))
  }
}

.sab_me_rescaled_exp_chains <- function(log_chains) {
  if (!is.list(log_chains) || !length(log_chains) ||
      any(!vapply(log_chains, is.numeric, logical(1L)))) {
    stop("Log-weight chains must be a nonempty list of numeric vectors.",
         call. = FALSE)
  }
  all_values <- unlist(log_chains, use.names = FALSE)
  if (!length(all_values) || anyNA(all_values) || any(is.nan(all_values)) ||
      any(all_values == Inf) || !any(is.finite(all_values))) {
    stop("Log-weight chains must contain valid finite contributions.",
         call. = FALSE)
  }
  common_scale <- max(all_values)
  lapply(log_chains, function(value) exp(value - common_scale))
}

.sab_me_reference <- function(canonical_anchor, adapter) {
  priorcentral <- adapter$prior_reference$eta
  if (!is.numeric(priorcentral) ||
      !identical(names(priorcentral), names(canonical_anchor$eta)) ||
      any(!is.finite(priorcentral)) ||
      !isTRUE(adapter$eta_in_domain(priorcentral))) {
    stop("Pinned prior-centred eta is malformed.", call. = FALSE)
  }
  list(
    kind = "equal_defensive_mixture",
    weights = c(saem = 0.5, priorcentral = 0.5),
    eta_components = list(
      saem = canonical_anchor$eta,
      priorcentral = priorcentral
    )
  )
}

.sab_me_patient_context <- function(artifact, patient_id) {
  context <- artifact$patient_context$population_covariates
  if (!is.numeric(context) || !identical(names(context), "treat_nelf") ||
      !context[["treat_nelf"]] %in% c(0, 1)) {
    stop("Patient bank has malformed covariate context: ", patient_id, ".",
         call. = FALSE)
  }
  list(patient_id = patient_id, treat_nelf = unname(context[["treat_nelf"]]))
}

.sab_me_validate_inputs <- function(adapter, canonical_anchor, artifacts,
                                    design, expected_role,
                                    expected_reference) {
  if (!is.list(artifacts) || !identical(names(artifacts), design$patient_ids)) {
    stop("Patient-bank artifacts are not in the predeclared audited order.",
         call. = FALSE)
  }
  for (patient_id in design$patient_ids) {
    artifact <- artifacts[[patient_id]]
    valid <- is.list(artifact) &&
      identical(artifact$schema_version,
                "sab_system_a_reference_patient_banks_v1") &&
      identical(artifact$patient_id, patient_id) &&
      identical(artifact$bank_role, expected_role) &&
      identical(artifact$target$target_fingerprint,
                adapter$target_fingerprint) &&
      identical(artifact$conditional_target$reference,
                expected_reference) &&
      identical(artifact$conditional_target$psi, canonical_anchor$psi) &&
      identical(
        unname(artifact$patient_context$population_covariates[["treat_nelf"]]),
        unname(adapter$treatment[[patient_id]])
      ) &&
      identical(names(artifact$chains),
                sprintf("chain_%02d", seq_along(artifact$chains))) &&
      all(c(design$pilot_chains, design$evaluation_chains) %in%
            names(artifact$chains))
    if (!isTRUE(valid)) {
      stop("Malformed or incompatible patient-bank artifact: ", patient_id,
           ".", call. = FALSE)
    }
    .sab_me_validate_chains(
      artifact$chains, adapter$coordinate_names$local,
      expected_reference, canonical_anchor$psi, patient_id,
      mode = "defensive_independence"
    )
  }
  invisible(TRUE)
}

.sab_me_validate_plan <- function(plan, adapter, canonical_anchor, design) {
  valid <- inherits(plan, "sab_system_a_message_endpoint_plan") &&
    identical(plan$schema_version, "sab_system_a_message_endpoint_plan_v1") &&
    identical(plan$target_fingerprint, adapter$target_fingerprint) &&
    isTRUE(plan$fixed_shared_parameters) &&
    identical(plan$psi, canonical_anchor$psi) &&
    identical(plan$reference, .sab_me_reference(canonical_anchor, adapter)) &&
    identical(plan$patient_ids, design$patient_ids) &&
    identical(plan$anchor$eta, canonical_anchor$eta) &&
    is.data.frame(plan$endpoints) && nrow(plan$endpoints) == 2L &&
    identical(names(plan$endpoints), c(
      "endpoint_id", "stage", "kind", "axis", "parameter", "sign",
      "step", "pilot_gate_passed", "diagnostic_fallback"
    )) &&
    is.matrix(plan$endpoint_eta) && nrow(plan$endpoint_eta) == 2L &&
    identical(rownames(plan$endpoint_eta), plan$endpoints$endpoint_id) &&
    identical(colnames(plan$endpoint_eta), names(canonical_anchor$eta))
  if (!isTRUE(valid)) stop("Malformed or incompatible endpoint plan.", call. = FALSE)
  invisible(TRUE)
}

.sab_me_validate_candidate <- function(artifact, adapter, canonical_anchor,
                                       plan, endpoint, plan_sha256,
                                       patient_id, candidate_eta) {
  target <- artifact$conditional_target
  valid <- is.list(artifact) &&
    identical(artifact$schema_version,
              "sab_system_a_endpoint_patient_banks_v1") &&
    identical(artifact$bank_role, "candidate_endpoint") &&
    identical(artifact$patient_id, patient_id) &&
    identical(artifact$target$target_fingerprint,
              adapter$target_fingerprint) &&
    identical(target$endpoint$endpoint_id, endpoint$endpoint_id) &&
    identical(target$eta, candidate_eta) &&
    identical(target$psi, canonical_anchor$psi) &&
    identical(
      unname(artifact$patient_context$population_covariates[["treat_nelf"]]),
      unname(adapter$treatment[[patient_id]])
    ) &&
    identical(target$endpoint_plan$sha256, plan_sha256) &&
    length(artifact$chains) >= 2L
  if (!isTRUE(valid)) {
    stop("Malformed candidate bank for ", endpoint$endpoint_id, ", patient ",
         patient_id, ".", call. = FALSE)
  }
  .sab_me_validate_chains(
    artifact$chains, adapter$coordinate_names$local,
    list(eta_components = list(saem = candidate_eta)),
    canonical_anchor$psi, patient_id, mode = "pcn"
  )
  invisible(TRUE)
}

.sab_me_validate_chains <- function(chains, local_names, expected_reference,
                                    expected_psi, patient_id,
                                    mode = c("pcn", "defensive_independence")) {
  mode <- match.arg(mode)
  valid <- is.list(chains) && length(chains) >= 2L &&
    !is.null(names(chains)) && !anyNA(names(chains)) &&
    !anyDuplicated(names(chains)) &&
    all(vapply(chains, function(chain) {
      is.list(chain) &&
        inherits(chain, "sab_system_a_patient_bank") &&
        identical(chain$patient_id, patient_id) &&
        is.matrix(chain$draws) && is.numeric(chain$draws) &&
        nrow(chain$draws) >= 4L &&
        identical(colnames(chain$draws), local_names) &&
        all(is.finite(chain$draws)) &&
        is.numeric(chain$loglik) &&
        length(chain$loglik) == nrow(chain$draws) &&
        all(is.finite(chain$loglik)) &&
        identical(chain$anchor$eta,
                  expected_reference$eta_components$saem) &&
        identical(chain$anchor$psi, expected_psi) &&
        if (mode == "pcn") {
          identical(chain$anchor$population_reference, "pcn") &&
            identical(chain$anchor$pcn_reference, "diagonal_gaussian") &&
            isTRUE(chain$anchor$pcn_reference_checked) &&
            identical(chain$proposal$mode, "pcn")
        } else {
          identical(chain$anchor$population_reference,
                    "defensive_independence") &&
            identical(chain$anchor$defensive_eta,
                      expected_reference$eta_components$priorcentral) &&
            identical(chain$anchor$defensive_mixture_weights,
                      expected_reference$weights) &&
            isTRUE(chain$anchor$defensive_reference_checked) &&
            identical(chain$proposal$mode, "defensive_independence")
        }
    }, logical(1L)))
  if (!isTRUE(valid)) {
    stop("Malformed conditional MCMC chains for patient ", patient_id, ".",
         call. = FALSE)
  }
  invisible(TRUE)
}

.sab_me_diagnostic_row <- function(axis, parameter, step, sign, patient_id,
                                   chain, direction, diagnostic,
                                   endpoint_id = NA_character_) {
  data.frame(
    endpoint_id = endpoint_id, axis = as.character(axis),
    parameter = as.character(parameter), step = as.numeric(step),
    sign = as.integer(sign), patient_id = as.character(patient_id),
    chain = as.character(chain), direction = direction,
    raw_log_ratio = diagnostic$log_ratio,
    oriented_log_ratio = diagnostic$log_ratio,
    relative_weight_ess = diagnostic$relative_weight_ess,
    max_normalized_weight = diagnostic$max_normalized_weight,
    d2 = diagnostic$d2,
    batch_log_mcse = diagnostic$batch$log_scale_mcse,
    split_log_ratio_difference =
      diagnostic$split$absolute_log_ratio_difference,
    stringsAsFactors = FALSE
  )
}

.sab_me_component_identity <- function(
    adapter, canonical_anchor, artifacts, chain_names, density) {
  design <- sab_system_a_message_design()
  reference <- .sab_me_reference(canonical_anchor, adapter)
  components <- list(
    saem = list(kind = "population", eta = canonical_anchor$eta),
    priorcentral = list(
      kind = "population", eta = adapter$prior_reference$eta
    )
  )
  rows <- list()
  index <- 0L
  for (patient_id in design$patient_ids) {
    artifact <- artifacts[[patient_id]]
    context <- .sab_me_patient_context(artifact, patient_id)
    patient_start <- index + 1L
    responsibilities <- list()
    saem_log_weights <- prior_log_weights <- list()
    for (chain_name in chain_names) {
      chain <- artifact$chains[[chain_name]]
      saem <- sab_evaluate_population_message(
        chain$draws, reference, components$saem, context, density
      )
      prior <- sab_evaluate_population_message(
        chain$draws, reference, components$priorcentral, context, density
      )
      probability_saem <- 0.5 * exp(saem$log_weights)
      probability_priorcentral <- 0.5 * exp(prior$log_weights)
      responsibilities[[chain_name]] <- probability_saem
      saem_log_weights[[chain_name]] <- saem$log_weights
      prior_log_weights[[chain_name]] <- prior$log_weights
      identity <- probability_saem + probability_priorcentral
      index <- index + 1L
      rows[[index]] <- data.frame(
        patient_id = patient_id,
        treatment = context$treat_nelf,
        chain = chain_name,
        mean_probability_saem = mean(probability_saem),
        mean_probability_priorcentral = mean(probability_priorcentral),
        maximum_pointwise_identity_error = max(abs(identity - 1)),
        stringsAsFactors = FALSE
      )
    }
    responsibility_diagnostic <- sab_mcmc_scalar_diagnostics(responsibilities)
    saem_weight_diagnostic <- sab_mcmc_scalar_diagnostics(
      .sab_me_rescaled_exp_chains(saem_log_weights)
    )
    prior_weight_diagnostic <- sab_mcmc_scalar_diagnostics(
      .sab_me_rescaled_exp_chains(prior_log_weights)
    )
    patient_rows <- patient_start:index
    for (row in patient_rows) {
      rows[[row]]$responsibility_split_rhat <-
        responsibility_diagnostic$split_rhat
      rows[[row]]$responsibility_mcmc_ess <-
        responsibility_diagnostic$mcmc_ess
      rows[[row]]$responsibility_relative_mcmc_ess <-
        responsibility_diagnostic$relative_mcmc_ess
      rows[[row]]$saem_weight_split_rhat <-
        saem_weight_diagnostic$split_rhat
      rows[[row]]$saem_weight_mcmc_ess <-
        saem_weight_diagnostic$mcmc_ess
      rows[[row]]$prior_weight_split_rhat <-
        prior_weight_diagnostic$split_rhat
      rows[[row]]$prior_weight_mcmc_ess <-
        prior_weight_diagnostic$mcmc_ess
    }
  }
  result <- do.call(rbind, rows)
  rownames(result) <- NULL
  result
}

.sab_me_component_pilot_gate <- function(pilot, identity, design) {
  gates <- design$pilot_gates
  probability_ranges <- vapply(
    split(identity$mean_probability_saem, identity$patient_id),
    function(value) diff(range(value)), numeric(1L)
  )
  reference_passed <-
    max(identity$maximum_pointwise_identity_error) <=
      gates$maximum_mixture_identity_error &&
    max(probability_ranges) <=
      gates$maximum_component_probability_chain_range &&
    max(c(identity$responsibility_split_rhat,
          identity$saem_weight_split_rhat,
          identity$prior_weight_split_rhat)) <=
      gates$maximum_split_rhat &&
    min(c(identity$responsibility_mcmc_ess,
          identity$saem_weight_mcmc_ess,
          identity$prior_weight_mcmc_ess)) >= gates$minimum_mcmc_ess
  do.call(rbind, lapply(design$core_endpoints$endpoint_id, function(endpoint_id) {
    selected <- pilot[pilot$endpoint_id == endpoint_id, ]
    cohort <- aggregate(raw_log_ratio ~ chain, selected, sum)
    pass <- reference_passed &&
      all(selected$relative_weight_ess >=
            gates$minimum_relative_weight_ess) &&
      all(selected$max_normalized_weight <=
            gates$maximum_normalized_weight) &&
      all(selected$split_log_ratio_difference <=
            gates$maximum_split_log_ratio_difference) &&
      diff(range(cohort$raw_log_ratio)) <=
        gates$maximum_cohort_chain_range
    data.frame(
      endpoint_id = endpoint_id,
      minimum_relative_weight_ess = min(selected$relative_weight_ess),
      maximum_normalized_weight = max(selected$max_normalized_weight),
      maximum_split_log_ratio_difference =
        max(selected$split_log_ratio_difference),
      cohort_chain_range = diff(range(cohort$raw_log_ratio)),
      maximum_component_probability_chain_range = max(probability_ranges),
      maximum_mixture_identity_error =
        max(identity$maximum_pointwise_identity_error),
      maximum_responsibility_split_rhat =
        max(identity$responsibility_split_rhat),
      minimum_responsibility_mcmc_ess =
        min(identity$responsibility_mcmc_ess),
      maximum_component_weight_split_rhat = max(c(
        identity$saem_weight_split_rhat,
        identity$prior_weight_split_rhat
      )),
      minimum_component_weight_mcmc_ess = min(c(
        identity$saem_weight_mcmc_ess,
        identity$prior_weight_mcmc_ess
      )),
      passed = pass,
      stringsAsFactors = FALSE
    )
  }))
}

.sab_me_bridge <- function(anchor_evaluations, candidate_evaluations) {
  q0_on_0 <- q1_on_0 <- q0_on_1 <- q1_on_1 <- numeric()
  for (value in anchor_evaluations) {
    q0_on_0 <- c(q0_on_0,
                 value$chain$loglik + value$population$log_g_anchor)
    q1_on_0 <- c(q1_on_0,
                 value$chain$loglik + value$population$log_g_candidate)
  }
  for (value in candidate_evaluations) {
    # Candidate evaluations have swapped labels: `anchor` is g_candidate and
    # `candidate` is g_SAEM-anchor.
    q0_on_1 <- c(q0_on_1,
                 value$chain$loglik + value$population$log_g_candidate)
    q1_on_1 <- c(q1_on_1,
                 value$chain$loglik + value$population$log_g_anchor)
  }
  sab_bridge_log_ratio(q0_on_0, q1_on_0, q0_on_1, q1_on_1)
}

.sab_me_bridge_scores <- function(anchor_evaluations, candidate_evaluations,
                                  log_ratio) {
  n0 <- sum(vapply(
    anchor_evaluations, function(value) length(value$population$log_weights),
    integer(1L)
  ))
  n1 <- sum(vapply(
    candidate_evaluations,
    function(value) length(value$population$log_weights), integer(1L)
  ))
  log_s0 <- log(n0) - log(n0 + n1)
  log_s1 <- log(n1) - log(n0 + n1)
  log_score <- function(log_density_ratio, side) {
    first <- log_s0 + log_ratio
    second <- log_s1 + log_density_ratio
    largest <- pmax(first, second)
    log_denominator <- largest +
      log(exp(first - largest) + exp(second - largest))
    if (side == "reference") {
      log_density_ratio - log_denominator
    } else {
      -log_denominator
    }
  }
  reference <- lapply(anchor_evaluations, function(value) {
    log_score(value$population$log_weights, "reference")
  })
  candidate <- lapply(candidate_evaluations, function(value) {
    log_score(-value$population$log_weights, "candidate")
  })
  # A common side-specific scale does not change autocorrelation or split-Rhat
  # and avoids numerical overflow when exponentiating bridge contributions.
  reference_scale <- max(unlist(reference, use.names = FALSE))
  candidate_scale <- max(unlist(candidate, use.names = FALSE))
  list(
    reference = lapply(reference, function(value) exp(value - reference_scale)),
    candidate = lapply(candidate, function(value) exp(value - candidate_scale))
  )
}

.sab_me_bridge_identity <- function(patients) {
  saem <- patients[patients$endpoint_id == "saem_component", ]
  prior <- patients[patients$endpoint_id == "priorcentral_component", ]
  prior <- prior[match(saem$patient_id, prior$patient_id), ]
  if (nrow(saem) == 0L || nrow(prior) != nrow(saem) ||
      anyNA(prior$patient_id)) {
    stop("Cannot pair component bridge estimates by patient.", call. = FALSE)
  }
  closure <- 0.5 * exp(saem$bridge_log_ratio) +
    0.5 * exp(prior$bridge_log_ratio)
  conservative_instability <-
    0.5 * exp(saem$bridge_log_ratio) *
      saem$bridge_chain_pair_instability +
    0.5 * exp(prior$bridge_log_ratio) *
      prior$bridge_chain_pair_instability
  tolerance <- pmax(0.02, 2 * conservative_instability)
  data.frame(
    patient_id = saem$patient_id,
    treatment = saem$treatment,
    bridge_saem_log_ratio = saem$bridge_log_ratio,
    bridge_priorcentral_log_ratio = prior$bridge_log_ratio,
    mixture_identity = closure,
    absolute_identity_error = abs(closure - 1),
    conservative_identity_instability = conservative_instability,
    tolerance = tolerance,
    passed = abs(closure - 1) <= tolerance,
    stringsAsFactors = FALSE
  )
}

.sab_me_endpoint_summary <- function(plan, directional, bridges, patients,
                                     identity, bridge_identity, mcmc, design) {
  gates <- design$evaluation_gates
  do.call(rbind, lapply(plan$endpoints$endpoint_id, function(endpoint_id) {
    endpoint <- plan$endpoints[plan$endpoints$endpoint_id == endpoint_id, ]
    patient <- patients[patients$endpoint_id == endpoint_id, ]
    diagnostic <- directional[directional$endpoint_id == endpoint_id, ]
    bridge <- bridges[bridges$endpoint_id == endpoint_id, ]
    forward <- diagnostic[diagnostic$direction == "forward", ]
    reverse <- diagnostic[diagnostic$direction == "reverse", ]
    forward_replicates <- aggregate(
      oriented_log_ratio ~ chain, forward, sum
    )$oriented_log_ratio
    reverse_replicates <- aggregate(
      oriented_log_ratio ~ chain, reverse, sum
    )$oriented_log_ratio
    bridge_chain_pair_totals <- aggregate(
      log_ratio ~ anchor_chain + candidate_chain, bridge, sum
    )$log_ratio
    forward_mcse <- sqrt(sum(vapply(
      split(forward, forward$patient_id),
      function(value) mean(value$batch_log_mcse^2) / nrow(value),
      numeric(1L)
    )))
    reverse_mcse <- sqrt(sum(vapply(
      split(reverse, reverse$patient_id),
      function(value) mean(value$batch_log_mcse^2) / nrow(value),
      numeric(1L)
    )))
    bridge_spread <- stats::sd(bridge_chain_pair_totals)
    forward_ratio <- sum(patient$forward_log_ratio)
    reverse_ratio <- sum(patient$reverse_log_ratio)
    bridge_ratio <- sum(patient$bridge_log_ratio)
    forward_tolerance <- max(
      gates$minimum_agreement_tolerance,
      2 * sqrt(forward_mcse^2 + bridge_spread^2)
    )
    reverse_tolerance <- max(
      gates$minimum_agreement_tolerance,
      2 * sqrt(reverse_mcse^2 + bridge_spread^2)
    )
    overlap_passed <-
      min(diagnostic$relative_weight_ess) >=
        gates$minimum_relative_weight_ess &&
      max(diagnostic$max_normalized_weight) <=
        gates$maximum_normalized_weight
    replication_passed <-
      diff(range(forward_replicates)) <=
        gates$maximum_cohort_chain_range &&
      diff(range(reverse_replicates)) <=
        gates$maximum_cohort_chain_range &&
      diff(range(bridge_chain_pair_totals)) <=
        gates$maximum_cohort_chain_range
    agreement_passed <-
      abs(forward_ratio - bridge_ratio) <= forward_tolerance &&
      abs(reverse_ratio - bridge_ratio) <= reverse_tolerance &&
      all(patient$forward_bridge_passed) &&
      all(patient$reverse_bridge_passed)
    component_ranges <- vapply(
      split(identity$mean_probability_saem, identity$patient_id),
      function(value) diff(range(value)), numeric(1L)
    )
    reference_passed <-
      max(identity$maximum_pointwise_identity_error) <=
        gates$maximum_mixture_identity_error &&
      max(component_ranges) <=
        gates$maximum_component_probability_chain_range &&
      max(c(identity$responsibility_split_rhat,
            identity$saem_weight_split_rhat,
            identity$prior_weight_split_rhat)) <=
        gates$maximum_split_rhat &&
      min(c(identity$responsibility_mcmc_ess,
            identity$saem_weight_mcmc_ess,
            identity$prior_weight_mcmc_ess)) >= gates$minimum_mcmc_ess
    endpoint_mcmc <- mcmc[mcmc$endpoint_id == endpoint_id, ]
    mcmc_passed <-
      all(is.finite(endpoint_mcmc$split_rhat)) &&
      max(endpoint_mcmc$split_rhat) <= gates$maximum_split_rhat &&
      min(endpoint_mcmc$mcmc_ess) >= gates$minimum_mcmc_ess
    candidate_augmented_identity_passed <- all(bridge_identity$passed)
    data.frame(
      endpoint_id = endpoint_id, axis = endpoint$axis,
      sign = endpoint$sign, step = endpoint$step,
      pilot_gate_passed = endpoint$pilot_gate_passed,
      diagnostic_fallback = endpoint$diagnostic_fallback,
      forward_log_ratio = forward_ratio,
      reverse_log_ratio = reverse_ratio,
      bridge_log_ratio = bridge_ratio,
      forward_bridge_difference = abs(forward_ratio - bridge_ratio),
      reverse_bridge_difference = abs(reverse_ratio - bridge_ratio),
      forward_agreement_tolerance = forward_tolerance,
      reverse_agreement_tolerance = reverse_tolerance,
      minimum_relative_weight_ess =
        min(diagnostic$relative_weight_ess),
      maximum_normalized_weight =
        max(diagnostic$max_normalized_weight),
      forward_chain_range = diff(range(forward_replicates)),
      reverse_chain_range = diff(range(reverse_replicates)),
      bridge_pair_range = diff(range(bridge_chain_pair_totals)),
      overlap_passed = overlap_passed,
      replication_passed = replication_passed,
      agreement_passed = agreement_passed,
      maximum_component_probability_chain_range = max(component_ranges),
      maximum_mixture_identity_error =
        max(identity$maximum_pointwise_identity_error),
      reference_mixing_passed = reference_passed,
      mcmc_diagnostics_passed = mcmc_passed,
      candidate_augmented_bridge_identity_passed =
        candidate_augmented_identity_passed,
      validated = endpoint$pilot_gate_passed && overlap_passed &&
        replication_passed && agreement_passed && reference_passed &&
        mcmc_passed && candidate_augmented_identity_passed,
      stringsAsFactors = FALSE
    )
  }))
}
