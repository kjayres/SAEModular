if (!exists("sab_system_a_plan_message_endpoints", mode = "function")) {
  source(file.path(core_root, "R", "patient_messages.R"), local = FALSE)
  source(
    file.path(core_root, "R", "system_a_message_validation.R"),
    local = FALSE
  )
}

.sab_message_validation_fixture <- function(separation = 3) {
  design <- sab_system_a_message_design()
  eta_saem <- c(mu = separation, log_omega = 0)
  eta_prior <- c(mu = -separation, log_omega = 0)
  psi <- c(psi = 0)
  reference <- list(
    kind = "equal_defensive_mixture",
    weights = c(saem = 0.5, priorcentral = 0.5),
    eta_components = list(saem = eta_saem, priorcentral = eta_prior)
  )
  adapter <- list(
    target_fingerprint = "synthetic-message-target",
    numerical_target = "synthetic_exact_target",
    patient_ids = design$patient_ids,
    treatment = stats::setNames(
      rep(c(0L, 1L), length.out = length(design$patient_ids)),
      design$patient_ids
    ),
    coordinate_names = list(local = "x", population = names(eta_saem)),
    prior_reference = list(eta = eta_prior),
    eta_in_domain = function(eta) {
      is.numeric(eta) && identical(names(eta), c("mu", "log_omega")) &&
        all(is.finite(eta))
    },
    log_population_density = function(patient_id, x, eta) {
      stats::dnorm(x[["x"]], eta[["mu"]], exp(eta[["log_omega"]]),
                   log = TRUE)
    }
  )
  canonical_anchor <- list(eta = eta_saem, psi = psi)
  # Alternating draws make both halves and all four chains cover both mixture
  # components, while the two component centres remain six SD apart.
  x <- rep(c(-separation, separation), 500L)
  draws <- matrix(x, ncol = 1L, dimnames = list(NULL, "x"))
  make_chain <- function(patient_id) {
    structure(list(
      patient_id = patient_id,
      draws = draws,
      loglik = rep(0, nrow(draws)),
      anchor = list(
        eta = eta_saem, psi = psi,
        population_reference = "defensive_independence",
        pcn_reference = NA_character_, pcn_reference_checked = FALSE,
        defensive_eta = eta_prior,
        defensive_mixture_weights = c(saem = 0.5, priorcentral = 0.5),
        defensive_reference_checked = TRUE
      ),
      proposal = list(mode = "defensive_independence")
    ), class = c("sab_system_a_patient_bank", "list"))
  }
  artifacts <- stats::setNames(lapply(design$patient_ids, function(patient_id) {
    summary <- data.frame(
      patient_id = patient_id,
      chain = sprintf("chain_%02d", 1:4),
      sampling_acceptance = 1,
      min_ess_x = 1000,
      min_ess_x_squared = 1000,
      min_ess_x_per_ode = 1,
      exact_prediction_calls = 1000,
      exact_ode_integrations = 1000,
      saem_component_proposals = 500,
      saem_component_proposal_fraction = 0.5,
      saem_component_acceptance = 1,
      saem_component_prediction_failure_rate = 0,
      saem_component_invalid_rate = 0,
      saem_component_prediction_failures = 0,
      saem_component_nonfinite_loglik = 0,
      priorcentral_component_proposals = 500,
      priorcentral_component_proposal_fraction = 0.5,
      priorcentral_component_acceptance = 1,
      priorcentral_component_prediction_failure_rate = 0,
      priorcentral_component_invalid_rate = 0,
      priorcentral_component_prediction_failures = 0,
      priorcentral_component_nonfinite_loglik = 0
    )
    list(
      schema_version = "sab_system_a_reference_patient_banks_v1",
      patient_id = patient_id,
      bank_role = "defensive_reference",
      target = list(target_fingerprint = adapter$target_fingerprint),
      patient_context = list(population_covariates = c(
        treat_nelf = unname(adapter$treatment[[patient_id]])
      )),
      conditional_target = list(reference = reference, psi = psi),
      summary = summary,
      chains = stats::setNames(
        replicate(4L, make_chain(patient_id), simplify = FALSE),
        sprintf("chain_%02d", 1:4)
      )
    )
  }), design$patient_ids)
  list(
    design = design, adapter = adapter, anchor = canonical_anchor,
    reference = reference, artifacts = artifacts
  )
}

testthat::test_that("defensive reference plan uses separated components", {
  fixture <- .sab_message_validation_fixture()
  temporary <- tempfile(fileext = ".rds")
  saveRDS(fixture$anchor, temporary)
  on.exit(unlink(temporary), add = TRUE)
  paths <- stats::setNames(
    rep(normalizePath(temporary), length(fixture$design$patient_ids)),
    fixture$design$patient_ids
  )
  hashes <- stats::setNames(
    rep(paste(rep("a", 64L), collapse = ""), length(paths)), names(paths)
  )
  result <- sab_system_a_plan_message_endpoints(
    fixture$adapter, fixture$anchor, fixture$artifacts,
    temporary, paste(rep("b", 64L), collapse = ""), paths, hashes
  )

  testthat::expect_s3_class(
    result$plan, "sab_system_a_message_endpoint_plan"
  )
  testthat::expect_identical(
    result$plan$endpoints$endpoint_id,
    c("saem_component", "priorcentral_component")
  )
  testthat::expect_identical(
    result$plan$reference$weights,
    c(saem = 0.5, priorcentral = 0.5)
  )
  testthat::expect_true(all(result$selection$passed))
  testthat::expect_lt(
    max(result$component_identity$maximum_pointwise_identity_error),
    1e-12
  )
  testthat::expect_equal(
    result$component_identity$mean_probability_saem +
      result$component_identity$mean_probability_priorcentral,
    rep(1, nrow(result$component_identity)), tolerance = 1e-14
  )
})

testthat::test_that("component gate detects chains trapped in different modes", {
  fixture <- .sab_message_validation_fixture()
  identity <- expand.grid(
    patient_id = fixture$design$patient_ids,
    chain = fixture$design$pilot_chains,
    stringsAsFactors = FALSE
  )
  identity$treatment <- 0
  identity$mean_probability_saem <- 0.5
  identity$mean_probability_priorcentral <- 0.5
  identity$maximum_pointwise_identity_error <- 0
  identity$responsibility_split_rhat <- 1
  identity$responsibility_mcmc_ess <- 1000
  identity$responsibility_relative_mcmc_ess <- 1
  identity$saem_weight_split_rhat <- 1
  identity$saem_weight_mcmc_ess <- 1000
  identity$prior_weight_split_rhat <- 1
  identity$prior_weight_mcmc_ess <- 1000
  identity$mean_probability_saem[
    identity$patient_id == fixture$design$patient_ids[[1L]] &
      identity$chain == "chain_01"
  ] <- 0.95
  identity$mean_probability_saem[
    identity$patient_id == fixture$design$patient_ids[[1L]] &
      identity$chain == "chain_02"
  ] <- 0.05
  pilot <- do.call(rbind, lapply(
    fixture$design$core_endpoints$endpoint_id,
    function(endpoint_id) do.call(rbind, lapply(
      fixture$design$pilot_chains,
      function(chain) data.frame(
        endpoint_id = endpoint_id,
        patient_id = fixture$design$patient_ids,
        chain = chain,
        raw_log_ratio = 0,
        relative_weight_ess = 1,
        max_normalized_weight = 0.01,
        split_log_ratio_difference = 0
      )
    ))
  ))

  gate <- .sab_me_component_pilot_gate(pilot, identity, fixture$design)
  testthat::expect_false(any(gate$passed))
  testthat::expect_equal(
    gate$maximum_component_probability_chain_range,
    rep(0.9, nrow(gate))
  )
})

testthat::test_that("evaluation gates contain every component threshold", {
  gates <- sab_system_a_message_design()$evaluation_gates
  testthat::expect_true(all(c(
    "maximum_component_probability_chain_range",
    "maximum_mixture_identity_error", "maximum_split_rhat",
    "minimum_mcmc_ess"
  ) %in% names(gates)))
})

testthat::test_that("message MCMC contributions exponentiate on one stable scale", {
  log_chains <- list(
    chain_01 = c(1000, 999, -Inf, 998),
    chain_02 = c(997, 996, 995, 994)
  )
  contributions <- .sab_me_rescaled_exp_chains(log_chains)
  testthat::expect_equal(max(unlist(contributions)), 1)
  testthat::expect_equal(
    contributions$chain_01[[2L]] / contributions$chain_02[[1L]],
    exp(2), tolerance = 1e-14
  )
  testthat::expect_equal(contributions$chain_01[[3L]], 0)
})

testthat::test_that("stored treatment must match the pinned adapter", {
  fixture <- .sab_message_validation_fixture()
  patient_id <- fixture$design$patient_ids[[1L]]
  fixture$artifacts[[patient_id]]$patient_context$
    population_covariates[["treat_nelf"]] <-
      1 - fixture$adapter$treatment[[patient_id]]
  temporary <- tempfile(fileext = ".rds")
  saveRDS(fixture$anchor, temporary)
  on.exit(unlink(temporary), add = TRUE)
  paths <- stats::setNames(rep(temporary, length(fixture$design$patient_ids)),
                           fixture$design$patient_ids)
  hashes <- stats::setNames(rep(strrep("a", 64L), length(paths)), names(paths))
  testthat::expect_error(
    sab_system_a_plan_message_endpoints(
      fixture$adapter, fixture$anchor, fixture$artifacts,
      temporary, strrep("b", 64L), paths, hashes
    ),
    "Malformed or incompatible patient-bank artifact"
  )
})

.sab_message_candidate_artifacts <- function(fixture, plan, plan_sha256) {
  make_chain <- function(patient_id, eta) {
    draws <- matrix(
      rep(eta[["mu"]], 1000L), ncol = 1L,
      dimnames = list(NULL, "x")
    )
    structure(list(
      patient_id = patient_id,
      draws = draws,
      loglik = rep(0, nrow(draws)),
      anchor = list(
        eta = eta, psi = fixture$anchor$psi,
        population_reference = "pcn",
        pcn_reference = "diagonal_gaussian",
        pcn_reference_checked = TRUE
      ),
      proposal = list(mode = "pcn")
    ), class = c("sab_system_a_patient_bank", "list"))
  }
  stats::setNames(lapply(plan$endpoints$endpoint_id, function(endpoint_id) {
    endpoint <- plan$endpoints[
      plan$endpoints$endpoint_id == endpoint_id, , drop = FALSE
    ]
    eta <- plan$endpoint_eta[endpoint_id, ]
    names(eta) <- colnames(plan$endpoint_eta)
    stats::setNames(lapply(fixture$design$patient_ids, function(patient_id) {
      list(
        schema_version = "sab_system_a_endpoint_patient_banks_v1",
        bank_role = "candidate_endpoint",
        patient_id = patient_id,
        target = list(
          target_fingerprint = fixture$adapter$target_fingerprint
        ),
        patient_context = list(population_covariates = c(
          treat_nelf = unname(fixture$adapter$treatment[[patient_id]])
        )),
        conditional_target = list(
          endpoint = endpoint, eta = eta, psi = fixture$anchor$psi,
          endpoint_plan = list(sha256 = plan_sha256)
        ),
        chains = stats::setNames(
          replicate(2L, make_chain(patient_id, eta), simplify = FALSE),
          c("chain_01", "chain_02")
        )
      )
    }), fixture$design$patient_ids)
  }), plan$endpoints$endpoint_id)
}

testthat::test_that("full component assessment is coherent on an exact fixture", {
  # Coincident components make every finite-bank identity exact.  The earlier
  # separated-component test exercises the nontrivial defensive-reference
  # calculation; this fixture isolates assessment schemas and gate wiring.
  fixture <- .sab_message_validation_fixture(separation = 0)
  temporary <- tempfile(fileext = ".rds")
  saveRDS(fixture$anchor, temporary)
  on.exit(unlink(temporary), add = TRUE)
  paths <- stats::setNames(
    rep(normalizePath(temporary), length(fixture$design$patient_ids)),
    fixture$design$patient_ids
  )
  hashes <- stats::setNames(rep(strrep("a", 64L), length(paths)), names(paths))
  planned <- sab_system_a_plan_message_endpoints(
    fixture$adapter, fixture$anchor, fixture$artifacts,
    temporary, strrep("b", 64L), paths, hashes
  )
  plan_sha256 <- strrep("c", 64L)
  candidates <- .sab_message_candidate_artifacts(
    fixture, planned$plan, plan_sha256
  )
  assessed <- sab_system_a_assess_messages(
    fixture$adapter, fixture$anchor, planned$plan, plan_sha256,
    fixture$artifacts, candidates
  )

  testthat::expect_true(assessed$calibration_passed)
  testthat::expect_identical(
    assessed$status, "calibration_passed_ready_for_local_stress"
  )
  testthat::expect_equal(assessed$endpoint_summary$forward_log_ratio, c(0, 0))
  testthat::expect_equal(assessed$endpoint_summary$bridge_log_ratio, c(0, 0))
  testthat::expect_true(all(
    assessed$bridge_component_identity$mixture_identity == 1
  ))
  testthat::expect_true(all(c(
    "forward_raw_weight_contribution",
    "reverse_raw_weight_contribution", "reference_bridge_score",
    "candidate_bridge_score"
  ) %in% assessed$mcmc_diagnostics$quantity))
  testthat::expect_true(
    "candidate_augmented_bridge_identity_passed" %in%
      names(assessed$endpoint_summary)
  )
  testthat::expect_true(
    "bridge_chain_pair_instability" %in% names(assessed$patient_summary)
  )
})
