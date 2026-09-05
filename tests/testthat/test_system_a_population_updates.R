source(file.path(core_root, "R", "system_a_population_updates.R"), local = FALSE)

population_fixture <- function() {
  coords <- .sab_population_names()
  list(x = matrix(seq(-1, 1, length.out = 24), 3, 8,
                  dimnames = list(c("a", "b", "c"), coords$local)),
       eta = setNames(c(rep(.1, 8), -.3, rep(log(.7), 8)), coords$eta),
       treatment = c(a = 0L, b = 1L, c = 1L),
       prior_mean = setNames(seq(-.2, .2, length.out = 17), coords$eta),
       prior_sd = setNames(seq(.5, 1.5, length.out = 17), coords$eta))
}

population_joint_logdensity <- function(eta, f) {
  sum(stats::dnorm(eta, f$prior_mean, f$prior_sd, log = TRUE)) + sum(vapply(
    seq_len(nrow(f$x)), function(i) {
      mu <- eta[1:8]; mu[[8L]] <- mu[[8L]] + eta[[9L]] * f$treatment[[i]]
      sum(stats::dnorm(f$x[i, ], mu, exp(eta[10:17]), log = TRUE))
    }, numeric(1L)))
}

testthat::test_that("location conditional reproduces joint density ratios and covariance", {
  f <- population_fixture()
  conditional <- do.call(sab_system_a_location_conditional, f)
  density <- function(eta) {
    delta <- eta[8:9] - conditional$pair_mean
    sum(stats::dnorm(eta[1:7], conditional$independent_mean,
                    conditional$independent_sd, log = TRUE)) -
      .5 * drop(crossprod(delta, conditional$pair_precision %*% delta))
  }
  other <- f$eta; other[1:9] <- other[1:9] + seq(-.5, .5, length.out = 9)
  testthat::expect_equal(density(other) - density(f$eta),
    population_joint_logdensity(other, f) - population_joint_logdensity(f$eta, f),
    tolerance = 1e-12)
  transform <- backsolve(conditional$pair_chol_precision, diag(2))
  testthat::expect_equal(tcrossprod(transform), solve(conditional$pair_precision))
  testthat::expect_lt(solve(conditional$pair_precision)[1, 2], 0)
})

testthat::test_that("scale conditional has correct normalization and no extra Jacobian", {
  f <- population_fixture()
  for (j in c(1L, 8L)) {
    mu <- f$eta[[j]] + if (j == 8L) f$treatment * f$eta[[9L]] else 0
    sse <- sum((f$x[, j] - mu)^2)
    density <- function(s) sab_system_a_scale_logdensity(
      s, 3L, sse, f$prior_mean[[j + 9L]], f$prior_sd[[j + 9L]])
    other <- f$eta; other[[j + 9L]] <- other[[j + 9L]] + .4
    testthat::expect_equal(density(other[[j + 9L]]) - density(f$eta[[j + 9L]]),
      population_joint_logdensity(other, f) - population_joint_logdensity(f$eta, f),
      tolerance = 1e-12)
  }
  testthat::expect_identical(sab_system_a_scale_logdensity(-1000, 3, 1, 0, 1), -Inf)
  testthat::expect_true(is.finite(sab_system_a_scale_logdensity(-1000, 3, 0, 0, 1)))
})

testthat::test_that("slice cap fails explicitly and successful update preserves schema", {
  testthat::expect_error(.sab_population_slice(0, function(s) 0, max_evaluations = 8L),
                         "cap exceeded")
  testthat::expect_error(.sab_population_slice(0, function(s) NaN), "invalid value")
  set.seed(52)
  f <- population_fixture()
  result <- do.call(sab_system_a_population_update, f)
  testthat::expect_identical(names(result), names(f$eta))
  testthat::expect_true(all(is.finite(result)))
  testthat::expect_false(identical(result, f$eta))
  f$treatment <- f$treatment[c(2, 1, 3)]
  testthat::expect_error(do.call(sab_system_a_population_update, f), "aligned")
  f <- population_fixture(); f$eta <- rev(f$eta)
  testthat::expect_error(do.call(sab_system_a_population_update, f), "canonical")
  f <- population_fixture(); f$prior_sd[[1L]] <- 0
  testthat::expect_error(do.call(sab_system_a_population_update, f), "positive")
})
