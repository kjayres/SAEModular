if (!exists("sab_build_system_a_patient_bank", mode = "function")) {
  source(file.path(core_root, "R", "system_a_patient_bank.R"), local = FALSE)
}

.sab_make_patient_bank_mock <- function(
    likelihood = function(x, psi) 0,
    solver = function(x, psi) {
      list(ok = TRUE, x = x, ode_integrations = 1L)
    },
    population_log_density = function(x, eta) {
      sum(stats::dnorm(
        x,
        mean = unname(eta[c("mu_a", "mu_b")]),
        sd = exp(unname(eta[c("log_omega_a", "log_omega_b")])),
        log = TRUE
      ))
    }) {
  counts <- new.env(parent = emptyenv())
  counts$solve <- 0L
  counts$loglik <- 0L
  adapter <- list(
    target_fingerprint = "synthetic-target",
    patient_ids = "patient_1",
    coordinate_names = list(
      local = c("u_a", "u_b"),
      population = c("mu_a", "mu_b", "log_omega_a", "log_omega_b"),
      global = "psi"
    ),
    population_mean = function(patient_id, eta) {
      stats::setNames(unname(eta[c("mu_a", "mu_b")]), c("u_a", "u_b"))
    },
    solve_prediction = function(patient_id, x, psi) {
      counts$solve <- counts$solve + 1L
      solver(x, psi)
    },
    loglik_from_prediction = function(prediction, psi) {
      counts$loglik <- counts$loglik + 1L
      if (!is.list(prediction) || !is.logical(prediction$ok) ||
          length(prediction$ok) != 1L || is.na(prediction$ok)) {
        stop("synthetic malformed prediction", call. = FALSE)
      }
      if (!isTRUE(prediction$ok)) {
        if (!is.character(prediction$reason) ||
            length(prediction$reason) != 1L ||
            !prediction$reason %in% c("synthetic_equilibrium_failure")) {
          stop("synthetic unexpected prediction failure", call. = FALSE)
        }
        return(-Inf)
      }
      likelihood(prediction$x, psi)
    },
    log_population_density = function(patient_id, x, eta) {
      population_log_density(x, eta)
    },
    eta_in_domain = function(eta) all(is.finite(eta)),
    psi_in_domain = function(psi) all(is.finite(psi))
  )
  class(adapter) <- c("sab_patient_bank_mock_adapter", "list")
  list(adapter = adapter, counts = counts)
}

.sab_patient_bank_anchor <- function() {
  list(
    eta = c(mu_a = 1, mu_b = -2, log_omega_a = log(2),
            log_omega_b = log(0.5)),
    psi = c(psi = 0.25)
  )
}

testthat::test_that("constant likelihood exposes the exact prior-reversible pCN path", {
  mock <- .sab_make_patient_bank_mock()
  anchor <- .sab_patient_bank_anchor()
  initial <- c(u_a = 5, u_b = -4)
  beta <- 0.4
  seed <- 402L
  n_draws <- 6L

  bank <- sab_build_system_a_patient_bank(
    adapter = mock$adapter,
    patient_id = "patient_1",
    eta = anchor$eta,
    psi = anchor$psi,
    initial_candidates = initial,
    warmup = 0L,
    draws = n_draws,
    thin = 1L,
    initial_beta = beta,
    target_acceptance = 0.3,
    adaptation_block = 2L,
    seed = seed
  )

  population_mean <- c(u_a = 1, u_b = -2)
  population_sd <- c(u_a = 2, u_b = 0.5)
  expected <- matrix(
    NA_real_, nrow = n_draws, ncol = 2L,
    dimnames = list(NULL, names(initial))
  )
  set.seed(seed)
  current <- initial
  for (iteration in seq_len(n_draws)) {
    current <- population_mean + sqrt(1 - beta^2) *
      (current - population_mean) +
      beta * population_sd * stats::rnorm(2L)
    stats::runif(1L)
    expected[iteration, ] <- current
  }

  testthat::expect_s3_class(bank, "sab_system_a_patient_bank")
  testthat::expect_equal(bank$draws, expected, tolerance = 1e-14)
  testthat::expect_identical(colnames(bank$draws), c("u_a", "u_b"))
  testthat::expect_equal(bank$loglik, rep(0, n_draws))
  testthat::expect_equal(unname(bank$acceptance[["sampling"]]), 1)
  testthat::expect_true(is.na(bank$acceptance[["warmup"]]))
  testthat::expect_identical(names(bank$ess$x), c("u_a", "u_b"))
  testthat::expect_identical(names(bank$ess$x_squared), c("u_a", "u_b"))
  testthat::expect_equal(mock$counts$solve, 1L + n_draws)
  testthat::expect_equal(mock$counts$loglik, 1L + n_draws)
  testthat::expect_equal(
    bank$ledger$prediction_calls[bank$ledger$phase == "total"],
    1 + n_draws
  )
  testthat::expect_equal(
    bank$ledger$ode_integrations[bank$ledger$phase == "total"],
    1 + n_draws
  )
})

testthat::test_that("adaptation is warmup-only and current likelihood is cached", {
  likelihood <- function(x, psi) -0.5 * (x[[1L]]^2 + 2 * x[[2L]]^2)
  anchor <- .sab_patient_bank_anchor()
  build <- function(draws) {
    mock <- .sab_make_patient_bank_mock(likelihood = likelihood)
    bank <- sab_build_system_a_patient_bank(
      adapter = mock$adapter,
      patient_id = "patient_1",
      eta = anchor$eta,
      psi = anchor$psi,
      initial_candidates = c(u_a = 0, u_b = 0),
      warmup = 4L,
      draws = draws,
      thin = 2L,
      initial_beta = 0.3,
      target_acceptance = 0.4,
      adaptation_block = 2L,
      seed = 91L
    )
    list(bank = bank, counts = mock$counts)
  }
  short <- build(3L)
  long <- build(9L)

  testthat::expect_equal(
    short$bank$beta$history$iteration, c(2L, 4L)
  )
  testthat::expect_equal(
    short$bank$beta$history, long$bank$beta$history
  )
  testthat::expect_equal(short$bank$beta$final, long$bank$beta$final)
  testthat::expect_equal(short$counts$solve, 1L + 4L + 2L * 3L)
  testthat::expect_equal(short$counts$loglik, short$counts$solve)
  expected_loglik <- apply(short$bank$draws, 1L, likelihood, psi = anchor$psi)
  testthat::expect_equal(short$bank$loglik, expected_loglik)
  testthat::expect_equal(
    short$bank$final_loglik,
    likelihood(short$bank$final_x, anchor$psi)
  )
  testthat::expect_equal(
    short$bank$ledger$prediction_calls[
      short$bank$ledger$phase == "warmup"
    ],
    4
  )
  testthat::expect_equal(
    short$bank$ledger$prediction_calls[
      short$bank$ledger$phase == "sampling"
    ],
    6
  )
})

testthat::test_that("dispersed candidates fail over with an auditable ledger", {
  solver <- function(x, psi) {
    if (x[[1L]] == 99) {
      return(list(ok = FALSE, reason = "synthetic_equilibrium_failure"))
    }
    list(ok = TRUE, x = x, ode_integrations = 1L)
  }
  likelihood <- function(x, psi) {
    if (x[[1L]] == 88) -Inf else 0
  }
  mock <- .sab_make_patient_bank_mock(likelihood, solver)
  anchor <- .sab_patient_bank_anchor()
  candidates <- rbind(
    c(u_a = 99, u_b = 0),
    c(u_a = 88, u_b = 0),
    c(u_a = 0, u_b = 0)
  )

  bank <- sab_build_system_a_patient_bank(
    adapter = mock$adapter,
    patient_id = "patient_1",
    eta = anchor$eta,
    psi = anchor$psi,
    initial_candidates = candidates,
    warmup = 0L,
    draws = 4L,
    thin = 1L,
    initial_beta = 0.1,
    target_acceptance = 0.3,
    adaptation_block = 2L,
    seed = 7L
  )

  initialization <- bank$ledger[bank$ledger$phase == "initialization", ]
  testthat::expect_equal(bank$initial_candidate_index, 3L)
  testthat::expect_equal(initialization$prediction_calls, 3)
  testthat::expect_equal(initialization$ode_integrations, 2)
  testthat::expect_equal(initialization$prediction_failures, 1)
  testthat::expect_equal(initialization$nonfinite_loglik, 2)
  testthat::expect_equal(
    bank$rejection_reasons,
    data.frame(
      phase = rep("initialization", 2L),
      reason = c("synthetic_equilibrium_failure", "nonfinite_loglik"),
      count = c(1L, 1L),
      stringsAsFactors = FALSE
    )
  )
})

testthat::test_that("population-zero states are rejected before prediction", {
  population_log_density <- function(x, eta) {
    if (x[[1L]] == 10) -Inf else 0
  }
  mock <- .sab_make_patient_bank_mock(
    population_log_density = population_log_density
  )
  anchor <- .sab_patient_bank_anchor()
  candidates <- rbind(
    c(u_a = 10, u_b = 0),
    c(u_a = 0, u_b = 0)
  )
  bank <- sab_build_system_a_patient_bank(
    adapter = mock$adapter,
    patient_id = "patient_1",
    eta = anchor$eta,
    psi = anchor$psi,
    initial_candidates = candidates,
    warmup = 0L,
    draws = 4L,
    thin = 1L,
    initial_beta = 0.1,
    target_acceptance = 0.3,
    adaptation_block = 2L,
    seed = 7L
  )

  initialization <- bank$ledger[bank$ledger$phase == "initialization", ]
  testthat::expect_equal(bank$initial_candidate_index, 2L)
  testthat::expect_equal(initialization$prediction_calls, 1)
  testthat::expect_equal(initialization$ode_integrations, 1)
  testthat::expect_equal(initialization$population_density_rejections, 1)
  testthat::expect_equal(mock$counts$solve, 1L + 4L)
})

testthat::test_that("sealed failure validation cannot be bypassed", {
  anchor <- .sab_patient_bank_anchor()
  malformed <- .sab_make_patient_bank_mock(
    solver = function(x, psi) "not a prediction"
  )
  testthat::expect_error(
    sab_build_system_a_patient_bank(
      malformed$adapter, "patient_1", anchor$eta, anchor$psi,
      c(u_a = 0, u_b = 0), warmup = 0L, draws = 2L,
      initial_beta = 0.1, target_acceptance = 0.3,
      adaptation_block = 2L, seed = 1L
    ),
    "synthetic malformed prediction",
    fixed = TRUE
  )

  unknown_failure <- .sab_make_patient_bank_mock(
    solver = function(x, psi) {
      list(ok = FALSE, reason = "unrecognised_failure", ode_integrations = 0L)
    }
  )
  testthat::expect_error(
    sab_build_system_a_patient_bank(
      unknown_failure$adapter, "patient_1", anchor$eta, anchor$psi,
      c(u_a = 0, u_b = 0), warmup = 0L, draws = 2L,
      initial_beta = 0.1, target_acceptance = 0.3,
      adaptation_block = 2L, seed = 1L
    ),
    "synthetic unexpected prediction failure",
    fixed = TRUE
  )
})

testthat::test_that("rejected proposals cannot replace the cached state", {
  likelihood <- function(x, psi) {
    if (all(x == 0)) 0 else -Inf
  }
  mock <- .sab_make_patient_bank_mock(likelihood = likelihood)
  anchor <- .sab_patient_bank_anchor()
  initial <- c(u_a = 0, u_b = 0)
  bank <- sab_build_system_a_patient_bank(
    mock$adapter, "patient_1", anchor$eta, anchor$psi, initial,
    warmup = 0L, draws = 5L, initial_beta = 0.2,
    target_acceptance = 0.3, adaptation_block = 2L, seed = 14L
  )

  testthat::expect_equal(bank$draws, matrix(
    rep(initial, 5L), nrow = 5L, byrow = TRUE,
    dimnames = list(NULL, names(initial))
  ))
  testthat::expect_equal(bank$final_x, initial)
  testthat::expect_equal(bank$final_loglik, 0)
  testthat::expect_equal(unname(bank$acceptance[["sampling"]]), 0)
  testthat::expect_equal(mock$counts$loglik, 6L)
})

testthat::test_that("patient-bank inputs cannot silently reorder coordinates", {
  mock <- .sab_make_patient_bank_mock()
  anchor <- .sab_patient_bank_anchor()
  arguments <- list(
    adapter = mock$adapter,
    patient_id = "patient_1",
    eta = anchor$eta,
    psi = anchor$psi,
    initial_candidates = c(u_a = 0, u_b = 0),
    warmup = 0L,
    draws = 4L,
    thin = 1L,
    initial_beta = 0.2,
    target_acceptance = 0.3,
    adaptation_block = 2L,
    seed = 3L
  )

  reordered_candidates <- arguments
  reordered_candidates$initial_candidates <- c(u_b = 0, u_a = 0)
  testthat::expect_error(
    do.call(sab_build_system_a_patient_bank, reordered_candidates),
    "canonical local order",
    fixed = TRUE
  )

  reordered_eta <- arguments
  reordered_eta$eta <- rev(anchor$eta)
  testthat::expect_error(
    do.call(sab_build_system_a_patient_bank, reordered_eta),
    "eta must be finite and named in canonical order",
    fixed = TRUE
  )
})
