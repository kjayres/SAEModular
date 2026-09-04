source(file.path(core_root, "R", "patient_messages.R"), local = FALSE)
source(file.path(core_root, "R", "system_a_message_validation.R"),
       local = FALSE)
source(file.path(core_root, "R", "system_a_pure_message_validation.R"),
       local = FALSE)

.sab_pure_test_fixture <- function(draws_per_chain = 500L) {
  design <- sab_system_a_pure_message_design()
  local_names <- c(
    "u_log_lambda", "u_log_mu_t", "u_log_mu_a", "u_log_p",
    "u_log_alpha_l", "u_pi", "u_eta_rti", "u_eta_pi"
  )
  eta_names <- c(
    "mu_log_lambda", "mu_log_mu_t", "mu_log_mu_a", "mu_log_p",
    "mu_log_alpha_l", "mu_u_pi", "mu_u_eta_rti", "mu_u_eta_pi",
    "beta_nelf", "log_omega_lambda", "log_omega_mu_t",
    "log_omega_mu_a", "log_omega_p", "log_omega_alpha_l",
    "log_omega_pi", "log_omega_eta_rti", "log_omega_eta_pi"
  )
  eta <- stats::setNames(rep(0, length(eta_names)), eta_names)
  eta[["log_omega_eta_pi"]] <- log(0.6)
  psi <- c(psi = 0)
  treatment <- stats::setNames(
    c(rep(0L, 8L), rep(1L, 4L)),
    design$patient_ids
  )
  population_mean <- function(patient_id, parameter) {
    value <- parameter[seq_len(8L)]
    names(value) <- local_names
    value[["u_eta_pi"]] <- value[["u_eta_pi"]] +
      treatment[[patient_id]] * parameter[["beta_nelf"]]
    value
  }
  adapter <- list(
    target_fingerprint = "synthetic-pure-message-target",
    numerical_target = "synthetic_exact_target",
    patient_ids = design$patient_ids,
    treatment = treatment,
    coordinate_names = list(local = local_names, population = eta_names),
    eta_in_domain = function(parameter) {
      is.numeric(parameter) && identical(names(parameter), eta_names) &&
        all(is.finite(parameter))
    },
    log_population_density = function(patient_id, x, parameter) {
      mean <- population_mean(patient_id, parameter)
      sd <- exp(parameter[10:17])
      sum(stats::dnorm(x, mean = mean, sd = sd, log = TRUE))
    }
  )
  anchor <- list(eta = eta, psi = psi)
  make_chain <- function(patient_id, chain_index, chain_eta = eta,
                         seed_namespace = 0L) {
    set.seed(
      81000L + seed_namespace + 100L * treatment[[patient_id]] + chain_index
    )
    mean <- population_mean(patient_id, chain_eta)
    sd <- exp(chain_eta[10:17])
    marginal <- stats::qnorm(
      (seq_len(draws_per_chain) - 0.5) / draws_per_chain
    )
    draws <- vapply(seq_along(local_names), function(index) {
      sample(marginal, length(marginal), replace = FALSE)
    }, numeric(draws_per_chain))
    draws <- sweep(draws, 2L, sd, "*")
    draws <- sweep(draws, 2L, mean, "+")
    colnames(draws) <- local_names
    structure(list(
      patient_id = patient_id,
      draws = draws,
      loglik = rep(0, draws_per_chain),
      anchor = list(
        eta = chain_eta, psi = psi,
        population_reference = "pcn",
        pcn_reference = "diagonal_gaussian",
        pcn_reference_checked = TRUE
      ),
      proposal = list(mode = "pcn"),
      beta = list(
        initial = 0.4, target_acceptance = 0.3,
        adaptation_block = 50L
      ),
      run = list(warmup = 250L, draws = 500L, thin = 1L)
    ), class = c("sab_system_a_patient_bank", "list"))
  }
  reference <- list(kind = "population", eta = eta)
  artifacts <- stats::setNames(lapply(design$patient_ids, function(patient_id) {
    chains <- stats::setNames(
      lapply(seq_len(4L), function(index) make_chain(patient_id, index)),
      sprintf("chain_%02d", seq_len(4L))
    )
    summary <- data.frame(
      patient_id = patient_id,
      chain = names(chains),
      sampling_acceptance = 0.3,
      min_ess_x = draws_per_chain,
      min_ess_x_squared = draws_per_chain,
      exact_prediction_calls = draws_per_chain + 251L,
      exact_ode_integrations = draws_per_chain + 250L,
      stringsAsFactors = FALSE
    )
    list(
      schema_version = "sab_system_a_pure_reference_patient_banks_v1",
      bank_role = "pure_saem_reference",
      patient_id = patient_id,
      target = list(target_fingerprint = adapter$target_fingerprint),
      patient_context = list(population_covariates = c(
        treat_nelf = unname(treatment[[patient_id]])
      )),
      conditional_target = list(
        eta = eta, psi = psi, reference = reference
      ),
      configuration = list(
        chains = 4L, warmup = 250L, draws = 500L, thin = 1L,
        initial_beta = 0.4, target_acceptance = 0.3,
        adaptation_block = 50L, start_offset_sd = 1,
        maximum_initial_candidates_per_chain = 18L,
        exact_prediction_call_cap = 3072
      ),
      summary = summary,
      chains = chains
    )
  }), design$patient_ids)
  temporary <- tempfile(fileext = ".rds")
  saveRDS(anchor, temporary)
  paths <- stats::setNames(
    rep(normalizePath(temporary), length(design$patient_ids)),
    design$patient_ids
  )
  hashes <- stats::setNames(rep(strrep("a", 64L), length(paths)), names(paths))
  list(
    design = design, adapter = adapter, anchor = anchor,
    artifacts = artifacts, temporary = temporary,
    paths = paths, hashes = hashes, make_chain = make_chain
  )
}

.sab_pure_test_candidates <- function(fixture, plan, plan_sha256) {
  endpoints <- plan$endpoints
  result <- stats::setNames(vector("list", nrow(endpoints)),
                            endpoints$endpoint_id)
  for (endpoint_index in seq_len(nrow(endpoints))) {
    endpoint <- endpoints[endpoint_index, , drop = FALSE]
    endpoint_id <- endpoint$endpoint_id[[1L]]
    eta <- plan$endpoint_eta[endpoint_id, ]
    names(eta) <- colnames(plan$endpoint_eta)
    result[[endpoint_id]] <- stats::setNames(lapply(
      fixture$design$patient_ids, function(patient_id) {
        chains <- stats::setNames(lapply(seq_len(4L), function(chain_index) {
          fixture$make_chain(
            patient_id, chain_index, chain_eta = eta,
            seed_namespace = 1000L * endpoint_index
          )
        }), sprintf("chain_%02d", seq_len(4L)))
        summary <- data.frame(
          patient_id = patient_id,
          chain = names(chains),
          sampling_acceptance = 0.3,
          min_ess_x = 500,
          min_ess_x_squared = 500,
          exact_prediction_calls = 751,
          exact_ode_integrations = 750,
          stringsAsFactors = FALSE
        )
        list(
          schema_version = "sab_system_a_endpoint_patient_banks_v1",
          bank_role = "candidate_endpoint",
          patient_id = patient_id,
          target = list(target_fingerprint = fixture$adapter$target_fingerprint),
          patient_context = list(population_covariates = c(
            treat_nelf = unname(fixture$adapter$treatment[[patient_id]])
          )),
          conditional_target = list(
            endpoint = endpoint, eta = eta, psi = fixture$anchor$psi,
            endpoint_plan = list(sha256 = plan_sha256)
          ),
          configuration = list(
            chains = 4L, warmup = 250L, draws = 500L, thin = 1L,
            initial_beta = 0.4, target_acceptance = 0.3,
            adaptation_block = 50L, start_offset_sd = 1,
            maximum_initial_candidates_per_chain = 18L,
            exact_prediction_call_cap = 3072
          ),
          summary = summary,
          chains = chains
        )
      }
    ), fixture$design$patient_ids)
  }
  result
}

.sab_pure_force_pilot_pass <- function(plan) {
  plan$endpoints$pilot_gate_passed[] <- TRUE
  plan$endpoints$diagnostic_fallback[] <- FALSE
  plan$pilot_selection$passed[] <- TRUE
  plan
}

testthat::test_that("pure endpoints obey the finite-variance contract", {
  fixture <- .sab_pure_test_fixture(20L)
  on.exit(unlink(fixture$temporary), add = TRUE)
  resolved <- .sab_pure_endpoints(fixture$adapter, fixture$anchor)
  endpoints <- resolved$endpoints

  testthat::expect_identical(
    endpoints$endpoint_id,
    c("treatment_minus", "treatment_plus",
      "lambda_scale_minus", "lambda_scale_plus")
  )
  testthat::expect_equal(endpoints$step[1:2], rep(0.6, 2L))
  testthat::expect_equal(endpoints$step[3:4], rep(log(1.25), 2L))
  testthat::expect_lt(exp(endpoints$step[[4L]]), sqrt(2))
  testthat::expect_equal(
    resolved$endpoint_eta["treatment_minus", "beta_nelf"], -0.6
  )
  testthat::expect_equal(
    resolved$endpoint_eta["treatment_plus", "beta_nelf"], 0.6
  )
  testthat::expect_equal(
    exp(resolved$endpoint_eta["lambda_scale_plus", "log_omega_lambda"]),
    1.25
  )
  testthat::expect_equal(
    .sab_pure_variance_share_limit(
      endpoints[1L, , drop = FALSE],
      sab_system_a_pure_message_design()$pilot_gates
    ),
    0.50
  )
  testthat::expect_equal(
    .sab_pure_variance_share_limit(
      endpoints[3L, , drop = FALSE],
      sab_system_a_pure_message_design()$pilot_gates
    ),
    0.25
  )
})

testthat::test_that("115-patient MCSE projection is explicit and exact", {
  design <- sab_system_a_pure_message_design()
  patient_mcse <- rep(0.047, length(design$patient_ids))
  treatment <- c(rep(0, 8L), rep(1, 4L))
  expected <- sqrt(
    80 / 8 * sum(patient_mcse[treatment == 0]^2) +
      35 / 4 * sum(patient_mcse[treatment == 1]^2)
  )
  result <- .sab_pure_projection_summary(patient_mcse, treatment, design)
  testthat::expect_equal(
    result$projected_log_mcse, expected
  )
  testthat::expect_equal(result$delta_log_bias_proxy, expected^2 / 2)
  testthat::expect_equal(
    result$maximum_projected_variance_share, 10 / 115
  )
  testthat::expect_error(
    .sab_pure_projection_summary(patient_mcse[-1L], treatment, design),
    "malformed"
  )
})

testthat::test_that("pure plan is distinct and binds all bank hashes", {
  fixture <- .sab_pure_test_fixture()
  on.exit(unlink(fixture$temporary), add = TRUE)
  result <- sab_system_a_plan_pure_message_endpoints(
    fixture$adapter, fixture$anchor, fixture$artifacts,
    fixture$temporary, strrep("b", 64L), fixture$paths, fixture$hashes
  )

  testthat::expect_s3_class(
    result$plan, "sab_system_a_pure_message_endpoint_plan"
  )
  testthat::expect_identical(
    result$plan$schema_version,
    "sab_system_a_pure_message_endpoint_plan_v1"
  )
  testthat::expect_identical(
    result$plan$reference,
    list(kind = "population", eta = fixture$anchor$eta)
  )
  testthat::expect_identical(
    result$plan$anchor_banks$sha256, fixture$hashes
  )
  testthat::expect_equal(nrow(result$pilot_diagnostics), 4L * 12L * 2L)
  testthat::expect_true(all(c(
    "projected_115_forward_log_mcse", "minimum_relevant_mcmc_ess",
    "maximum_relevant_split_rhat", "passed"
  ) %in% names(result$selection)))

  defensive <- fixture$artifacts
  defensive[[1L]]$schema_version <-
    "sab_system_a_reference_patient_banks_v1"
  defensive[[1L]]$bank_role <- "defensive_reference"
  testthat::expect_error(
    sab_system_a_plan_pure_message_endpoints(
      fixture$adapter, fixture$anchor, defensive,
      fixture$temporary, strrep("b", 64L), fixture$paths, fixture$hashes
    ),
    "pure-SAEM reference artifact"
  )
})

testthat::test_that("rank, folded, and tail diagnostics expose chain mismatch", {
  set.seed(991L)
  matched <- list(stats::rnorm(500), stats::rnorm(500))
  good <- .sab_pure_rank_tail_diagnostics(matched)
  shifted <- matched
  shifted[[2L]] <- shifted[[2L]] + 3
  bad <- .sab_pure_rank_tail_diagnostics(shifted)

  testthat::expect_true(is.finite(good$tail_mcmc_ess))
  testthat::expect_gt(bad$split_rhat, 1.01)
  testthat::expect_true(
    bad$rank_normalized_split_rhat > 1.01 || bad$folded_split_rhat > 1.01
  )
})

testthat::test_that("reference plan enforces the frozen chain and cost contract", {
  fixture <- .sab_pure_test_fixture()
  on.exit(unlink(fixture$temporary), add = TRUE)

  wrong_run <- fixture$artifacts
  wrong_run[[1L]]$chains$chain_01$run$warmup <- 251L
  testthat::expect_error(
    sab_system_a_plan_pure_message_endpoints(
      fixture$adapter, fixture$anchor, wrong_run,
      fixture$temporary, strrep("b", 64L), fixture$paths, fixture$hashes
    ),
    "Malformed pure pCN chain"
  )

  costly <- fixture$artifacts
  costly[[1L]]$summary$exact_prediction_calls[[1L]] <- 1000
  planned <- sab_system_a_plan_pure_message_endpoints(
    fixture$adapter, fixture$anchor, costly,
    fixture$temporary, strrep("b", 64L), fixture$paths, fixture$hashes
  )
  testthat::expect_false(planned$pilot_passed)
  testthat::expect_true(all(!planned$selection$reference_cost_passed))

  excessive_ode <- fixture$artifacts
  excessive_ode[[1L]]$summary$exact_ode_integrations[[1L]] <- 1000
  planned_ode <- sab_system_a_plan_pure_message_endpoints(
    fixture$adapter, fixture$anchor, excessive_ode,
    fixture$temporary, strrep("b", 64L), fixture$paths, fixture$hashes
  )
  testthat::expect_false(planned_ode$pilot_passed)
  testthat::expect_true(all(!planned_ode$selection$reference_cost_passed))

  wrong_configuration <- fixture$artifacts
  wrong_configuration[[1L]]$configuration$initial_beta <- 0.2
  planned_configuration <- sab_system_a_plan_pure_message_endpoints(
    fixture$adapter, fixture$anchor, wrong_configuration,
    fixture$temporary, strrep("b", 64L), fixture$paths, fixture$hashes
  )
  testthat::expect_false(planned_configuration$pilot_passed)
  testthat::expect_true(all(
    !planned_configuration$selection$reference_cost_passed
  ))

  configuration_mutations <- list(
    initial_beta = 0.2,
    target_acceptance = 0.4,
    adaptation_block = 25L,
    start_offset_sd = 1.5
  )
  for (field in names(configuration_mutations)) {
    artifact <- fixture$artifacts[[1L]]
    artifact$configuration[[field]] <- configuration_mutations[[field]]
    testthat::expect_false(
      .sab_pure_bank_budget_status(artifact, fixture$design)$configuration_ok,
      info = paste("configuration field", field)
    )
  }

  beta_mutations <- list(
    initial = 0.2,
    target_acceptance = 0.4,
    adaptation_block = 25L
  )
  for (field in names(beta_mutations)) {
    chain <- fixture$artifacts[[1L]]$chains$chain_01
    chain$beta[[field]] <- beta_mutations[[field]]
    testthat::expect_error(
      .sab_pure_validate_chain(
        chain, fixture$adapter, fixture$design$patient_ids[[1L]],
        fixture$anchor$eta, fixture$anchor$psi,
        fixture$design$bank_contract
      ),
      "Malformed pure pCN chain",
      info = paste("chain beta field", field)
    )
  }
})

testthat::test_that("assessment uses disjoint banks and can pass exact normals", {
  fixture <- .sab_pure_test_fixture()
  on.exit(unlink(fixture$temporary), add = TRUE)
  planned <- sab_system_a_plan_pure_message_endpoints(
    fixture$adapter, fixture$anchor, fixture$artifacts,
    fixture$temporary, strrep("b", 64L), fixture$paths, fixture$hashes
  )
  plan <- .sab_pure_force_pilot_pass(planned$plan)
  plan_sha256 <- strrep("c", 64L)
  candidates <- .sab_pure_test_candidates(fixture, plan, plan_sha256)
  assessed <- sab_system_a_assess_pure_messages(
    fixture$adapter, fixture$anchor, plan, plan_sha256,
    fixture$artifacts, candidates
  )

  testthat::expect_true(assessed$falsification_passed)
  testthat::expect_true(all(assessed$endpoint_summary$agreement_passed))
  testthat::expect_true(all(assessed$endpoint_summary$cost_passed))
  testthat::expect_true(all(assessed$cost_summary$within_frozen_budget))
  testthat::expect_identical(
    unique(assessed$directional_diagnostics$chain[
      assessed$directional_diagnostics$direction == "forward"
    ]), c("chain_03", "chain_04")
  )
  testthat::expect_identical(
    unique(assessed$directional_diagnostics$chain[
      assessed$directional_diagnostics$direction == "reverse"
    ]), c("chain_01", "chain_02")
  )
  testthat::expect_identical(
    sort(unique(assessed$bridge_chain_pairs$anchor_chain)),
    c("chain_01", "chain_02")
  )
  testthat::expect_identical(
    sort(unique(assessed$bridge_chain_pairs$candidate_chain)),
    c("chain_03", "chain_04")
  )
})

testthat::test_that("direct disagreement and excessive cost fail closed", {
  fixture <- .sab_pure_test_fixture()
  on.exit(unlink(fixture$temporary), add = TRUE)
  planned <- sab_system_a_plan_pure_message_endpoints(
    fixture$adapter, fixture$anchor, fixture$artifacts,
    fixture$temporary, strrep("b", 64L), fixture$paths, fixture$hashes
  )
  plan <- .sab_pure_force_pilot_pass(planned$plan)
  plan_sha256 <- strrep("d", 64L)
  candidates <- .sab_pure_test_candidates(fixture, plan, plan_sha256)

  costly <- candidates
  costly[["lambda_scale_plus"]][["3"]]$summary$exact_prediction_calls[[1L]] <-
    1000
  cost_result <- sab_system_a_assess_pure_messages(
    fixture$adapter, fixture$anchor, plan, plan_sha256,
    fixture$artifacts, costly
  )
  testthat::expect_false(cost_result$falsification_passed)
  testthat::expect_false(cost_result$endpoint_summary$cost_passed[
    cost_result$endpoint_summary$endpoint_id == "lambda_scale_plus"
  ])

  disagreeing <- candidates
  for (chain_name in c("chain_01", "chain_02")) {
    chain <- disagreeing[["lambda_scale_plus"]][["3"]]$chains[[chain_name]]
    chain$draws[, "u_log_lambda"] <- chain$draws[, "u_log_lambda"] + 6
    disagreeing[["lambda_scale_plus"]][["3"]]$chains[[chain_name]] <- chain
  }
  disagreement_result <- sab_system_a_assess_pure_messages(
    fixture$adapter, fixture$anchor, plan, plan_sha256,
    fixture$artifacts, disagreeing
  )
  testthat::expect_false(disagreement_result$falsification_passed)
  testthat::expect_false(disagreement_result$endpoint_summary$agreement_passed[
    disagreement_result$endpoint_summary$endpoint_id == "lambda_scale_plus"
  ])
})

testthat::test_that("bridge failure and substituted chain counts fail closed", {
  fixture <- .sab_pure_test_fixture()
  on.exit(unlink(fixture$temporary), add = TRUE)
  planned <- sab_system_a_plan_pure_message_endpoints(
    fixture$adapter, fixture$anchor, fixture$artifacts,
    fixture$temporary, strrep("b", 64L), fixture$paths, fixture$hashes
  )
  plan <- .sab_pure_force_pilot_pass(planned$plan)
  plan_sha256 <- strrep("e", 64L)
  candidates <- .sab_pure_test_candidates(fixture, plan, plan_sha256)

  failed_bridge <- sab_system_a_assess_pure_messages(
    fixture$adapter, fixture$anchor, plan, plan_sha256,
    fixture$artifacts, candidates,
    bridge_tolerance = 1e-16, bridge_max_iterations = 1L
  )
  testthat::expect_false(failed_bridge$falsification_passed)
  testthat::expect_true(any(
    failed_bridge$endpoint_summary$pooled_bridge_failures > 0L
  ))

  candidates[[1L]][[1L]]$chains$chain_05 <-
    candidates[[1L]][[1L]]$chains$chain_04
  testthat::expect_error(
    sab_system_a_assess_pure_messages(
      fixture$adapter, fixture$anchor, plan, plan_sha256,
      fixture$artifacts, candidates
    ),
    "Malformed candidate artifact"
  )
})
