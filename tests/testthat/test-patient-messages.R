source(file.path(core_root, "R", "patient_messages.R"), local = FALSE)

testthat::test_that("log message calculations are stable at extreme scales", {
  constant <- sab_raw_message_diagnostics(rep(10000, 20), n_batches = 4)

  testthat::expect_equal(constant$log_ratio, 10000, tolerance = 1e-12)
  testthat::expect_equal(constant$weight_ess, 20, tolerance = 1e-12)
  testthat::expect_equal(constant$relative_weight_ess, 1, tolerance = 1e-12)
  testthat::expect_equal(constant$max_normalized_weight, 0.05, tolerance = 1e-12)
  testthat::expect_equal(constant$d2, 0, tolerance = 1e-12)
  testthat::expect_equal(constant$batch$log_scale_mcse, 0, tolerance = 1e-12)
  testthat::expect_equal(
    constant$split$absolute_log_ratio_difference,
    0,
    tolerance = 1e-12
  )

  expected <- -10000 + log((1 + exp(-1)) / 2)
  testthat::expect_equal(
    sab_log_mean_exp(c(-10000, -10001)),
    expected,
    tolerance = 1e-12
  )
  testthat::expect_identical(sab_log_mean_exp(c(-Inf, -Inf)), -Inf)
})

testthat::test_that("raw anchor weights recover an analytic Gaussian message", {
  set.seed(20260902)
  n_draws <- 100000L
  eta_anchor <- -0.1
  eta_new <- 0.35
  y <- 0.4
  tau <- 0.8
  sigma <- 0.5
  conditional_variance <- 1 / (1 / tau^2 + 1 / sigma^2)
  conditional_mean <- conditional_variance *
    (eta_anchor / tau^2 + y / sigma^2)
  patient_draws <- stats::rnorm(
    n_draws,
    mean = conditional_mean,
    sd = sqrt(conditional_variance)
  )
  log_weights <-
    stats::dnorm(patient_draws, eta_new, tau, log = TRUE) -
    stats::dnorm(patient_draws, eta_anchor, tau, log = TRUE)

  estimate <- sab_raw_message_diagnostics(log_weights, n_batches = 20)
  marginal_sd <- sqrt(tau^2 + sigma^2)
  truth <-
    stats::dnorm(y, eta_new, marginal_sd, log = TRUE) -
    stats::dnorm(y, eta_anchor, marginal_sd, log = TRUE)

  testthat::expect_equal(estimate$log_ratio, truth, tolerance = 0.01)
  testthat::expect_gt(estimate$weight_ess, 0)
  testthat::expect_lte(estimate$weight_ess, n_draws)
  testthat::expect_equal(
    estimate$relative_weight_ess,
    exp(-estimate$d2),
    tolerance = 1e-12
  )
  testthat::expect_gt(estimate$batch$log_scale_mcse, 0)
  testthat::expect_lt(estimate$split$absolute_log_ratio_difference, 0.03)
})

testthat::test_that("population messages accept non-Gaussian covariate models", {
  # Four deterministic states have exactly the anchor Bernoulli frequencies.
  # The candidate probability is a covariate-dependent alternative model.
  states <- c(0L, 1L, 1L, 1L)
  anchor_probability <- 0.75
  patient_covariate <- 2
  candidate_coefficient <- 0.2
  candidate_probability <- stats::plogis(
    candidate_coefficient * patient_covariate
  )
  log_probability <- function(state, probability) {
    ifelse(state == 1L, log(probability), log1p(-probability))
  }

  log_weights <- sab_population_log_weights(
    log_g_candidate = log_probability(states, candidate_probability),
    log_g_anchor = log_probability(states, anchor_probability)
  )
  estimate <- sab_raw_message_diagnostics(log_weights, n_batches = 2L)

  testthat::expect_equal(estimate$log_ratio, 0, tolerance = 1e-14)
  testthat::expect_length(log_weights, length(states))
})

testthat::test_that("density callbacks retain covariates without family assumptions", {
  draws <- matrix(
    c(-1.5, -0.25, 0.5, 2), ncol = 1L,
    dimnames = list(NULL, "x")
  )
  context <- list(treatment = 1, dose = 2.5)
  anchor <- list(intercept = -0.2, treatment = 0.4, scale = 1.1, df = 5)
  candidate <- list(intercept = 0.1, treatment = -0.15, scale = 1.4, df = 3)
  log_student_population <- function(x, parameter, patient_context) {
    location <- parameter$intercept +
      parameter$treatment * patient_context$treatment +
      0.05 * patient_context$dose
    stats::dt(
      (x[["x"]] - location) / parameter$scale,
      df = parameter$df, log = TRUE
    ) - log(parameter$scale)
  }

  evaluated <- sab_evaluate_population_message(
    draws, anchor, candidate, context, log_student_population
  )
  expected_anchor <- apply(
    draws, 1L, log_student_population,
    parameter = anchor, patient_context = context
  )
  expected_candidate <- apply(
    draws, 1L, log_student_population,
    parameter = candidate, patient_context = context
  )

  testthat::expect_s3_class(
    evaluated, "sab_population_message_evaluation"
  )
  testthat::expect_equal(evaluated$log_g_anchor, expected_anchor)
  testthat::expect_equal(evaluated$log_g_candidate, expected_candidate)
  testthat::expect_equal(
    evaluated$log_weights, expected_candidate - expected_anchor
  )
})

testthat::test_that("density callback evaluation fails closed", {
  draws <- matrix(0, nrow = 2L, ncol = 1L,
                  dimnames = list(NULL, "x"))
  invalid <- function(x, parameter, patient_context) {
    if (identical(parameter, "candidate")) NaN else 0
  }
  testthat::expect_error(
    sab_evaluate_population_message(
      draws, "anchor", "candidate", list(), invalid
    ),
    "candidate log population density at draw 1",
    fixed = TRUE
  )
  testthat::expect_error(
    sab_evaluate_population_message(
      unname(draws), "anchor", "candidate", list(), invalid
    ),
    "unique coordinate names",
    fixed = TRUE
  )
})

testthat::test_that("equal defensive-component messages obey their identity", {
  states <- matrix(
    c(-8, -1, 0, 1, 9), ncol = 1L,
    dimnames = list(NULL, "x")
  )
  log_g_saem <- stats::dnorm(states[, 1L], -3, 0.8, log = TRUE)
  log_g_prior <- stats::dt(states[, 1L] / 2, df = 4, log = TRUE) - log(2)
  largest <- pmax(log_g_saem, log_g_prior)
  log_h <- log(0.5) + largest +
    log(exp(log_g_saem - largest) + exp(log_g_prior - largest))
  weight_saem <- exp(log_g_saem - log_h)
  weight_prior <- exp(log_g_prior - log_h)

  testthat::expect_equal(
    0.5 * weight_saem + 0.5 * weight_prior,
    rep(1, nrow(states)), tolerance = 2e-15
  )
  testthat::expect_lte(max(weight_saem), 2)
  testthat::expect_lte(max(weight_prior), 2)
})

testthat::test_that("MCMC diagnostics are distinct from weight concentration", {
  constant <- sab_mcmc_scalar_diagnostics(list(rep(2, 20), rep(2, 20)))
  testthat::expect_equal(constant$split_rhat, 1)
  testthat::expect_equal(constant$mcmc_ess, 40)

  separated <- sab_mcmc_scalar_diagnostics(list(rep(-1, 20), rep(1, 20)))
  testthat::expect_true(is.infinite(separated$split_rhat))
  testthat::expect_equal(separated$mcmc_ess, 1)

  persistent <- list(cumsum(rep(c(1, -0.9), 20)),
                     cumsum(rep(c(0.9, -1), 20)))
  diagnostic <- sab_mcmc_scalar_diagnostics(persistent)
  testthat::expect_true(is.finite(diagnostic$split_rhat))
  testthat::expect_gte(diagnostic$mcmc_ess, 1)
  testthat::expect_lte(diagnostic$mcmc_ess, 80)
})

testthat::test_that("population log weights validate observed density arrays", {
  testthat::expect_equal(
    sab_population_log_weights(c(0, -Inf), c(0, 0)),
    c(0, -Inf)
  )
  testthat::expect_error(
    sab_population_log_weights(c(0, 0), c(0, -Inf)),
    "finite values",
    fixed = TRUE
  )
  testthat::expect_error(
    sab_population_log_weights(c(0, 0), 0),
    "equal length",
    fixed = TRUE
  )
})

testthat::test_that("batch log MCSE matches a hand calculation", {
  raw_weights <- c(1, 3, 2, 4, 5, 7)
  estimate <- sab_batch_log_mcse(log(raw_weights), n_batches = 2)
  batch_means <- c(2, 16 / 3)
  overall_mean <- mean(raw_weights)
  expected <- sqrt(stats::var(batch_means) / 2) / overall_mean

  testthat::expect_equal(estimate$batch_sizes, c(3L, 3L))
  testthat::expect_equal(
    estimate$batch_log_ratios,
    log(batch_means),
    tolerance = 1e-12
  )
  testthat::expect_equal(estimate$log_scale_mcse, expected, tolerance = 1e-12)
})

testthat::test_that("nonconstant weight ESS and D2 match direct arithmetic", {
  raw_weights <- c(1, 2, 4, 8)
  estimate <- sab_raw_message_diagnostics(log(raw_weights), n_batches = 2)
  direct_ess <- sum(raw_weights)^2 / sum(raw_weights^2)

  testthat::expect_equal(estimate$weight_ess, direct_ess, tolerance = 1e-12)
  testthat::expect_equal(
    estimate$d2,
    log(length(raw_weights) / direct_ess),
    tolerance = 1e-12
  )
})

testthat::test_that("bidirectional bridge sampling recovers known normal constants", {
  set.seed(76123)
  n0 <- 60000L
  n1 <- 50000L
  shifted_mean <- 0.45
  shifted_sd <- 1.15
  log_amplitude <- 0.3
  draws_0 <- stats::rnorm(n0)
  draws_1 <- stats::rnorm(n1, shifted_mean, shifted_sd)
  log_q0 <- function(x) -0.5 * x^2
  log_q1 <- function(x) {
    log_amplitude - 0.5 * ((x - shifted_mean) / shifted_sd)^2
  }

  estimate <- sab_bridge_log_ratio(
    log_q0_on_0 = log_q0(draws_0),
    log_q1_on_0 = log_q1(draws_0),
    log_q0_on_1 = log_q0(draws_1),
    log_q1_on_1 = log_q1(draws_1)
  )
  truth <- log_amplitude + log(shifted_sd)

  testthat::expect_true(estimate$converged)
  testthat::expect_lt(estimate$iterations, 100L)
  testthat::expect_equal(estimate$log_ratio, truth, tolerance = 0.01)
})

testthat::test_that("bridge sampling is exact for proportional densities", {
  log_q0_on_0 <- c(-3, -1, 0, -2)
  log_q0_on_1 <- c(-0.5, -4, -1.5, -2.5, -3.5)
  log_constant <- 700
  estimate <- sab_bridge_log_ratio(
    log_q0_on_0,
    log_q0_on_0 + log_constant,
    log_q0_on_1,
    log_q0_on_1 + log_constant
  )

  testthat::expect_equal(estimate$log_ratio, log_constant, tolerance = 1e-12)

  reversed <- sab_bridge_log_ratio(
    log_q0_on_1 + log_constant,
    log_q0_on_1,
    log_q0_on_0 + log_constant,
    log_q0_on_0
  )
  testthat::expect_equal(
    estimate$log_ratio,
    -reversed$log_ratio,
    tolerance = 1e-12
  )
})

testthat::test_that("message estimators fail closed on malformed input", {
  testthat::expect_error(sab_log_mean_exp(numeric()), "length at least")
  testthat::expect_error(sab_log_mean_exp(c(0, NA_real_)), "must not contain")
  testthat::expect_error(sab_log_mean_exp(c(0, Inf)), "positive infinity")
  testthat::expect_error(sab_log_mean_exp(matrix(1:4, 2, 2)), "numeric vector")
  testthat::expect_error(
    sab_raw_message_diagnostics(rep(-Inf, 8), n_batches = 2),
    "at least one finite"
  )
  testthat::expect_error(
    sab_raw_message_diagnostics(1:3, n_batches = 2),
    "length at least 4"
  )
  testthat::expect_error(
    sab_batch_log_mcse(1:8, n_batches = 5),
    "at least two draws"
  )
  testthat::expect_error(
    sab_raw_message_diagnostics(c(rep(0, 4), rep(-Inf, 4)), n_batches = 2),
    "second-half"
  )

  testthat::expect_error(
    sab_bridge_log_ratio(1:3, 1:2, 1:3, 1:3),
    "q0 draws must have equal length"
  )
  testthat::expect_error(
    sab_bridge_log_ratio(c(0, -Inf), c(0, 0), c(0, 0), c(0, 0)),
    "only finite"
  )
  testthat::expect_error(
    sab_bridge_log_ratio(1:3, 1:3, 1:3, 1:3, tolerance = 0),
    "finite positive"
  )
  testthat::expect_error(
    sab_bridge_log_ratio(1:3, 1:3, 1:3, 1:3, max_iterations = 0),
    "integer at least 1"
  )
  testthat::expect_error(
    sab_bridge_log_ratio(
      c(0, 0), c(-1, 0), c(0, 0), c(0, 1),
      tolerance = 1e-16,
      max_iterations = 1,
      initial_log_ratio = 100
    ),
    "did not converge"
  )
})
