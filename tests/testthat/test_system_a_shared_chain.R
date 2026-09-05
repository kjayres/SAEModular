# Synthetic algebraic forward maps only. The test launcher runs these on Slurm;
# no sealed ODE, comparator data, or scientific inference is needed here.
shared_chain_fixture <- function(dynamic_bound = Inf) {
  engine <- new.env(parent = globalenv())
  for (file in c("system_a_population_updates.R", "shared_pool_kernel.R",
                 "system_a_shared_chain.R")) {
    sys.source(file.path(core_root, "R", file), envir = engine)
  }
  coordinates <- engine$.sab_population_names()
  ids <- c("a", "b", "equilibrium")
  treatment <- setNames(c(0L, 1L, 0L), ids)
  prior_mean <- setNames(rep(0, 17), coordinates$eta)
  prior_sd <- setNames(rep(0.8, 17), coordinates$eta)
  psi_names <- c("log_gamma_pop", "log_mu_l_pop", "log_sigma_v", "log_sigma_t")
  initial <- list(
    x = matrix(seq(-0.3, 0.3, length.out = 24), 3L, 8L,
                dimnames = list(ids, coordinates$local)),
    eta = prior_mean,
    psi = setNames(rep(0, 4), psi_names)
  )
  recorder <- new.env(parent = emptyenv())
  recorder$exact_calls <- 0L
  recorder$common_calls <- 0L
  recorder$likelihood_calls <- 0L
  recorder$stale_calls <- 0L
  recorder$common_dynamic_states <- list()
  observations <- list(a = c(-0.5, 0.1), b = c(0.3, -0.2), equilibrium = c(0.2, 0.4))
  predict <- function(id, x, psi) list(
    ok = TRUE, patient_id = id, x = x,
    adjusted_positive_times = if (id == "equilibrium") numeric() else 1,
    dynamic_psi = psi[1:2],
    observation_mean = c(x[[1L]] + psi[[1L]], x[[2L]] + psi[[2L]])
  )
  adapter <- list(
    patient_ids = ids, treatment = treatment,
    coordinate_names = list(local = coordinates$local, population = coordinates$eta,
                            global = psi_names),
    target_fingerprint = "synthetic_exact_target",
    population_mean = function(id, eta) {
      mu <- setNames(as.numeric(eta[1:8]), coordinates$local)
      mu[[8L]] <- mu[[8L]] + eta[[9L]] * treatment[[id]]
      mu
    },
    eta_in_domain = function(eta) all(is.finite(eta)) && all(is.finite(exp(eta[10:17]))),
    psi_in_domain = function(psi) all(is.finite(psi)) && all(is.finite(exp(psi[3:4]))),
    log_hyperprior = function(eta, psi) {
      if (any(abs(psi[1:2]) > dynamic_bound)) return(-Inf)
      sum(stats::dnorm(eta, prior_mean, prior_sd, log = TRUE)) +
        sum(stats::dnorm(psi, 0, 1, log = TRUE))
    },
    solve_prediction = function(id, x, psi) {
      recorder$exact_calls <- recorder$exact_calls + 1L
      predict(id, x, psi)
    },
    loglik_from_prediction = function(prediction, psi) {
      recorder$likelihood_calls <- recorder$likelihood_calls + 1L
      if (!identical(prediction$dynamic_psi, psi[1:2])) {
        recorder$stale_calls <- recorder$stale_calls + 1L
        stop("Synthetic likelihood received a stale dynamic-psi cache.")
      }
      sum(stats::dnorm(observations[[prediction$patient_id]],
                      prediction$observation_mean, exp(psi[3:4]), log = TRUE))
    }
  )
  engine$.sab_shared_original_calls <- function(prediction) {
    as.integer(length(prediction$adjusted_positive_times) > 0L)
  }
  engine$.sab_shared_solve <- function(probe, x, psi) {
    recorder$common_calls <- recorder$common_calls + 1L
    recorder$common_dynamic_states[[recorder$common_calls]] <- psi[1:2]
    list(ok = TRUE, x = x, dynamic_psi = psi[1:2], ode_integrations = 1L)
  }
  engine$.sab_shared_prediction <- function(probe, solved, id, psi) {
    if (!identical(solved$dynamic_psi, psi[1:2])) stop("Stale synthetic common solve.")
    predict(id, solved$x, psi)
  }
  probe <- list(patients = list(a = list(), b = list()),
                positive_times = list(a = 1, b = 1),
                reference_target_fingerprint = "synthetic_exact_target",
                target_fingerprint = "synthetic_common_proposal")
  second <- initial
  second$eta[[1L]] <- 0.4
  second$eta[[9L]] <- -0.7
  second$eta[10:17] <- log(seq(0.7, 1.4, length.out = 8))
  components <- engine$sab_shared_pool_reference_components(adapter, list(initial, second))
  list(engine = engine, adapter = adapter, probe = probe, initial = initial,
       components = components, prior_mean = prior_mean, prior_sd = prior_sd,
       recorder = recorder)
}

shared_chain_run_fixture <- function(f, method, budget = 1000L, max_sweeps = 18L) {
  f$engine$sab_run_system_a_shared_chain(
    f$adapter, f$probe, f$initial, f$components, f$prior_mean, f$prior_sd,
    method = method, seed = 472L,
    config = list(budget = budget, warmup = 3L, dynamic_every = 1L,
                  reserve_refresh = 1L, pool_size = 6L, max_sweeps = max_sweeps,
                  initial_dynamic_sd = c(0.02, 0.02), checkpoint_every = 100L))
}

testthat::test_that("vectorized population density and equal-mixture q match direct formulas", {
  f <- shared_chain_fixture()
  x <- rbind(f$initial$x, -2 * f$initial$x)
  mu <- seq(-0.4, 0.3, length.out = 8)
  sd <- seq(0.6, 1.3, length.out = 8)
  expected <- apply(x, 1L, function(value) sum(stats::dnorm(value, mu, sd, log = TRUE)))
  testthat::expect_equal(unname(f$engine$.sab_chain_normal_logdensity(x, mu, sd)),
                         unname(expected))
  testthat::expect_length(f$components, 4L)
  mixture <- apply(x, 1L, function(value) {
    log(mean(vapply(f$components, function(component) {
      prod(stats::dnorm(value, component$mean, component$sd))
    }, numeric(1L))))
  })
  testthat::expect_equal(unname(f$engine$.sab_chain_log_q(x, f$components)),
                         unname(mixture), tolerance = 1e-12)
  testthat::expect_equal(unname(f$engine$.sab_chain_log_q(x[1L, , drop = FALSE], f$components)),
                         unname(mixture[[1L]]), tolerance = 1e-12)
})

testthat::test_that("both complete short chains retain globals and honor the forward budget", {
  for (method in c("baseline", "recycling")) {
    f <- shared_chain_fixture()
    result <- shared_chain_run_fixture(f, method, budget = 60L, max_sweeps = 100L)
    testthat::expect_gt(nrow(result$draws), result$warmup)
    testthat::expect_identical(ncol(result$draws), 22L)
    testthat::expect_identical(colnames(result$draws),
      c(names(f$initial$eta), names(f$initial$psi), "treated_logit_location"))
    testthat::expect_true(all(is.finite(result$draws)))
    testthat::expect_equal(unname(result$draws[, "treated_logit_location"]),
      unname(result$draws[, "mu_u_eta_pi"] + result$draws[, "beta_nelf"]))
    testthat::expect_lte(sum(result$ledger$budget_units), 60)
    testthat::expect_equal(tail(result$costs, 1L), sum(result$ledger$budget_units))
    testthat::expect_equal(sum(result$ledger$ode_integrations), sum(result$ledger$budget_units))
    testthat::expect_equal(sum(result$ledger$attempts[result$ledger$category != "common"]),
                           f$recorder$exact_calls)
    testthat::expect_equal(result$ledger$attempts[result$ledger$category == "common"],
                           f$recorder$common_calls)
    testthat::expect_identical(f$recorder$stale_calls, 0L)
    testthat::expect_gt(result$counters[["dynamic_accepts"]], 0)
    testthat::expect_true(all(is.finite(result$final_state$x)))
    testthat::expect_equal(unname(result$local_attempts), rep(nrow(result$draws), 3L))
    if (method == "recycling") {
      testthat::expect_equal(result$counters[["assignment_attempts"]], 2 * nrow(result$draws))
      testthat::expect_equal(result$counters[["reserve_refreshes"]], nrow(result$draws))
      testthat::expect_gt(f$recorder$common_calls, 6L)
      testthat::expect_equal(result$counters[["assignment_moves"]], sum(result$assignment_moves))
    } else testthat::expect_identical(f$recorder$common_calls, 0L)
  }
})

testthat::test_that("rejected dynamic proposals keep all active prediction versions valid", {
  # A proper prior truncated to this narrow interval forces the proposed
  # 0.02-scale dynamic moves to reject for the frozen test seed.
  f <- shared_chain_fixture(dynamic_bound = 1e-12)
  result <- shared_chain_run_fixture(f, "recycling", max_sweeps = 8L)
  testthat::expect_equal(result$counters[["dynamic_attempts"]], 8)
  testthat::expect_equal(result$counters[["dynamic_accepts"]], 0)
  testthat::expect_equal(result$final_state$psi[1:2], f$initial$psi[1:2])
  testthat::expect_true(all(result$draws[, names(f$initial$psi)[1:2]] == 0))
  testthat::expect_identical(f$recorder$stale_calls, 0L)
  testthat::expect_true(all(vapply(f$recorder$common_dynamic_states,
    function(value) identical(value, f$initial$psi[1:2]), logical(1L))))
})
