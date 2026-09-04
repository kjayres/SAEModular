source(file.path(core_root, "R", "shared_trajectory_probe.R"), local = FALSE)

# Use pinned pure helpers, with a synthetic integrator: no ODE integration or
# inference is performed in these tests.
shared_probe_fixture <- function() {
  upstream <- new.env(parent = globalenv())
  source_root <- file.path(core_root, "..", "..", "projects", ".worktrees",
                          "modular_bayes_system_a_factor_v1", "research",
                          "system_a_modular_v0", "R")
  testthat::skip_if_not(file.exists(file.path(source_root, "exact_likelihood.R")))
  for (file in c("target_manifest.R", "gme_adapter.R", "exact_likelihood.R")) {
    sys.source(file.path(source_root, file), envir = upstream)
  }
  controls <- list(mu_v = 30, rel_tol = 1e-6, abs_tol = 1e-8,
                   max_num_steps = 100000L, time_eps = 1e-6)
  patient <- function(id, times, types) list(
    patient_id = id, obs_time = times, ytype = types,
    y = rep(1, length(times)), cens = rep(0L, length(times)), controls = controls
  )
  patients <- list(
    a = patient("a", c(-1, 2, 2), c(2L, 1L, 2L)),
    b = patient("b", c(0, 2, 3), c(1L, 2L, 1L)),
    baseline = patient("baseline", c(-1, 0), c(1L, 2L))
  )
  positive <- lapply(patients, function(p) {
    upstream$sysa_stan_positive_times(p$obs_time, controls$time_eps)
  })
  structure(list(
    upstream = upstream, patients = patients, controls = controls,
    positive_times = positive,
    union_times = sort(unique(as.numeric(unlist(positive)))),
    numerical_target = "experimental_union_grid_vode_bdf_v1"
  ), class = "sab_shared_trajectory_probe")
}

shared_probe_state <- function(probe) {
  list(x = setNames(c(log(100), log(.1), log(.2), log(100), log(.01), 0, 1, 1),
                    probe$upstream$sysa_local_coordinate_names()),
       psi = setNames(log(c(.003, .015, .3, .2)),
                      probe$upstream$sysa_psi_names()))
}

testthat::test_that("union output maps original jitter and outcome order exactly", {
  probe <- shared_probe_fixture()
  state <- shared_probe_state(probe)
  calls <- 0L
  fake_ode <- function(y, times, ...) {
    calls <<- calls + 1L
    testthat::expect_equal(times, c(0, 2, 2 + 1e-6, 3), tolerance = 0)
    cbind(time = times, outer(times, seq_len(5), "+"))
  }
  solved <- .sab_shared_solve(probe, state$x, state$psi, fake_ode)
  testthat::expect_identical(calls, 1L)
  testthat::expect_identical(solved$ode_integrations, 1L)
  a <- .sab_shared_prediction(probe, solved, "a", state$psi)
  b <- .sab_shared_prediction(probe, solved, "b", state$psi)
  testthat::expect_equal(a$observation_mean[[2L]], log10(13000))
  testthat::expect_equal(a$observation_mean[[3L]], 12 + 3e-6)
  testthat::expect_equal(b$observation_mean[[2L]], 12)
  testthat::expect_equal(a$adjusted_positive_times, c(2, 2 + 1e-6), tolerance = 0)
  testthat::expect_identical(a$y, probe$patients$a$y)
  testthat::expect_identical(a$solver, probe$numerical_target)
  testthat::expect_true(is.finite(probe$upstream$sysa_loglik_from_prediction(
    a, state$psi
  )))
  changed <- state$psi
  changed[[1L]] <- changed[[1L]] + .01
  testthat::expect_error(probe$upstream$sysa_loglik_from_prediction(a, changed),
                         "different dynamic psi")
})

testthat::test_that("union integration failure preserves equilibrium-only patients", {
  probe <- shared_probe_fixture()
  state <- shared_probe_state(probe)
  solved <- .sab_shared_solve(probe, state$x, state$psi, function(...) {
    warning("DVODE- excessive work done")
    NULL
  })
  testthat::expect_false(solved$ok)
  testthat::expect_identical(solved$ode_integrations, 1L)
  testthat::expect_identical(solved$reason, "ode_failure")
  testthat::expect_false(.sab_shared_prediction(probe, solved, "a", state$psi)$ok)
  baseline <- .sab_shared_prediction(probe, solved, "baseline", state$psi)
  testthat::expect_true(baseline$ok)
  testthat::expect_true(is.finite(probe$upstream$sysa_loglik_from_prediction(
    baseline, state$psi
  )))
  testthat::expect_error(.sab_shared_solve(probe, state$x, state$psi,
    function(...) stop("synthetic programming defect")), "programming defect")
  testthat::expect_error(.sab_shared_solve(probe, state$x, state$psi,
    function(...) warning("synthetic programming warning")), "Unexpected warning")
})

testthat::test_that("invalid equilibria and baseline-only design spend no integration", {
  probe <- shared_probe_fixture()
  state <- shared_probe_state(probe)
  should_not_run <- function(...) stop("Integrator should not be called")
  invalid <- state$x
  invalid[["u_log_lambda"]] <- -100
  solved <- .sab_shared_solve(probe, invalid, state$psi, should_not_run)
  testthat::expect_identical(solved$reason, "invalid_equilibrium")
  testthat::expect_identical(solved$ode_integrations, 0L)
  testthat::expect_false(.sab_shared_prediction(probe, solved,
                                               "baseline", state$psi)$ok)
  probe$union_times <- numeric()
  probe$patients <- probe$patients["baseline"]
  solved <- .sab_shared_solve(probe, state$x, state$psi, should_not_run)
  testthat::expect_true(solved$ok)
  testthat::expect_identical(solved$ode_integrations, 0L)
  prediction <- .sab_shared_prediction(probe, solved, "baseline", state$psi)
  testthat::expect_identical(.sab_shared_original_calls(prediction), 0L)
  testthat::expect_identical(.sab_shared_original_calls(
    list(ok = FALSE, reason = "ode_failure")), 1L)
  testthat::expect_error(.sab_shared_original_calls(
    list(ok = FALSE, reason = "unknown_patient")), "Unexpected original")
})

testthat::test_that("malformed union solver output is rejected, never interpolated", {
  probe <- shared_probe_fixture()
  state <- shared_probe_state(probe)
  solved <- .sab_shared_solve(probe, state$x, state$psi, function(...) matrix(1, 2, 6))
  testthat::expect_false(solved$ok)
  testthat::expect_identical(solved$reason, "invalid_ode_output")
  testthat::expect_identical(solved$ode_integrations, 1L)
})
