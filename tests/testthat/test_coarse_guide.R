source(file.path(core_root, "R", "system_a_coarse_guide.R"), local = FALSE)

coarse_uniforms <- function(values) {
  index <- 0L
  function() { index <<- index + 1L; values[[index]] }
}

testthat::test_that("identical targets give unit macro correction and atomic accepted caches", {
  target <- function(x) list(log_target = -.5 * x^2, cache = paste0("at_", x))
  current <- list(x = 0, exact = target(0), cheap = target(0))
  result <- sab_surrogate_macro_step(current,
    function(x) list(x = x + .1, log_reverse_minus_forward = 0), target, target,
    steps = 3, uniform = coarse_uniforms(rep(.01, 4)))
  testthat::expect_equal(result$current$x, .3)
  testthat::expect_equal(result$log_correction, 0)
  testthat::expect_equal(result$inner_accepts, 3)
  testthat::expect_equal(result$exact_calls, 1)
  testthat::expect_identical(result$current$exact$cache, "at_0.3")
  testthat::expect_identical(current$exact$cache, "at_0")
})

testthat::test_that("rejection and impossible exact endpoints preserve all active caches", {
  cheap <- function(x) list(log_target = 0, cache = x)
  current <- list(x = 0, exact = list(log_target = 0, cache = "original"), cheap = cheap(0))
  for (value in c(-100, -Inf)) {
    result <- sab_surrogate_macro_step(current,
      function(x) list(x = 1, log_reverse_minus_forward = 0), cheap,
      function(x) list(log_target = value, cache = "must_not_be_retained"), 1,
      uniform = coarse_uniforms(c(.1, .9)))
    testthat::expect_false(result$accepted)
    testthat::expect_identical(result$current, current)
  }
})

testthat::test_that("cheap proposal asymmetry and self transitions are handled", {
  current <- list(x = 0, exact = list(log_target = 0), cheap = list(log_target = 0))
  result <- sab_surrogate_macro_step(current,
    function(x) list(x = 1, log_reverse_minus_forward = -10),
    function(x) list(log_target = 0), function(x) stop("Exact call not needed"), 2,
    uniform = coarse_uniforms(c(.1, .1)))
  testthat::expect_false(result$moved)
  testthat::expect_equal(result$exact_calls, 0)
  testthat::expect_identical(result$current, current)
  testthat::expect_error(sab_surrogate_macro_step(current,
    function(x) list(x = 1, log_reverse_minus_forward = 0),
    function(x) list(log_target = NA_real_), function(x) list(log_target = 0), 1), "cheap density")
})

testthat::test_that("non-Gaussian finite targets obey detailed balance for macro lengths", {
  exact <- c(.1, .4, .2, .3); cheap <- c(.6, .1, .1, .2); proposal <- c(.1, .2, .5, .2)
  kernel <- matrix(0, 4, 4)
  for (i in 1:4) for (j in 1:4) {
    probability <- proposal[[j]] * min(1, cheap[[j]] * proposal[[i]] /
                                      (cheap[[i]] * proposal[[j]]))
    kernel[i, j] <- kernel[i, j] + probability
    kernel[i, i] <- kernel[i, i] + proposal[[j]] - probability
  }
  for (steps in c(1L, 3L, 8L)) {
    macro <- diag(4)
    for (k in seq_len(steps)) macro <- macro %*% kernel
    corrected <- matrix(0, 4, 4)
    for (i in 1:4) for (j in 1:4) {
      probability <- macro[i, j] * min(1, exact[[j]] / cheap[[j]] * cheap[[i]] / exact[[i]])
      corrected[i, j] <- corrected[i, j] + probability
      corrected[i, i] <- corrected[i, i] + macro[i, j] - probability
    }
    flow <- corrected * exact
    testthat::expect_equal(rowSums(corrected), rep(1, 4), tolerance = 1e-12)
    testthat::expect_equal(flow, t(flow), tolerance = 1e-12)
    testthat::expect_equal(as.numeric(exact %*% corrected), exact, tolerance = 1e-12)
  }
})

testthat::test_that("retained-chain mean and quantile diagnostics are available", {
  set.seed(6105)
  values <- matrix(stats::rnorm(4000), 1000, 4)
  testthat::expect_true(is.finite(posterior::ess_mean(values)))
  testthat::expect_length(posterior::mcse_quantile(values, probs = c(.05, .5, .95)), 3L)
})

testthat::test_that("proposal-only guide cannot change the sealed solver or false-zero its support", {
  testthat::skip_if_not(identical(Sys.getenv("SAB_SYSTEM_A_INTEGRATION"), "true"))
  source(file.path(core_root, "R", "shared_trajectory_probe.R"), local = FALSE)
  workspace <- normalizePath(file.path(core_root, "..", ".."))
  adapter <- sab_load_system_a_adapter(workspace, "3")
  owner <- environment(adapter$solve_prediction)
  before <- get("exact_likelihood", owner, inherits = FALSE)$patients[["3"]]$controls
  guide <- sab_make_system_a_coarse_guide(adapter, before$rel_tol, before$abs_tol)
  x <- setNames(c(log(100), log(.1), log(.2), log(100), log(.01), 0, 1, 1), adapter$coordinate_names$local)
  psi <- setNames(log(c(.003, .015, .3, .2)), adapter$coordinate_names$global)
  expected <- adapter$log_likelihood("3", x, psi)
  value <- guide$evaluate("3", x, psi)
  testthat::expect_true(is.finite(expected))
  testthat::expect_equal(value$raw_loglik, expected, tolerance = 0)
  testthat::expect_identical(value$source, "approximate")
  loose <- sab_make_system_a_coarse_guide(adapter, .1, .001)
  after <- get("exact_likelihood", owner, inherits = FALSE)$patients[["3"]]$controls
  testthat::expect_identical(before, after)
  invalid <- x; invalid[[1L]] <- -100
  invalid_value <- loose$evaluate("3", invalid, psi)
  testthat::expect_identical(invalid_value$raw_loglik, -Inf)
  testthat::expect_equal(invalid_value$loglik, -1e6)
  testthat::expect_equal(invalid_value$ode_integrations, 0)
})
