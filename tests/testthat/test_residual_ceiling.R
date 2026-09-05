source(file.path(core_root, "R", "residual_ceiling.R"), local = FALSE)

testthat::test_that("a perfectly matched component is not falsely ruled out", {
  result <- sab_residual_ceiling(rep(0, 1000), 0, caps = c(1, 2, 10))
  testthat::expect_equal(result$mean_capped_ratio, rep(1, 3))
  testthat::expect_equal(result$collapse_upper, rep(1, 3))
  shifted <- sab_residual_ceiling(rep(-10000, 1000), -10000, caps = c(1, 2, 10))
  testthat::expect_equal(shifted, result)
})

testthat::test_that("known shape mismatch gives conservative arithmetic ceilings", {
  # F(x)=1+9*x on [0,1], q=h=Uniform[0,1], c_probe=c*=1, M=5.5.
  # Midpoints test the arithmetic against analytic integrals, not IID coverage.
  x <- (seq_len(4096) - 0.5) / 4096
  result <- sab_residual_ceiling(log1p(9 * x), 0, alpha_per_cap = 0.01)
  testthat::expect_true(all(result$collapse_upper >= 1 / 5.5))
  testthat::expect_lt(min(result$collapse_upper), 0.2)
  testthat::expect_equal(result$mean_capped_ratio[2], 5.5)
  expected <- pmax(0, result$mean_capped_ratio -
                     result$cap * sqrt(log(100) / (2 * length(x))))
  testthat::expect_equal(result$mass_ratio_lower, expected)
})

testthat::test_that("partial support can hide residual mass but broader h can expose it", {
  # q=Uniform[0,1]. F=1 on [0,1], F=99 on (1,2], hence c*=1 and M=100.
  compact <- sab_residual_ceiling(rep(0, 4096), 0)
  testthat::expect_equal(compact$collapse_upper, rep(1, 3))
  # For h=Uniform[0,2], F/h is 2 or 198 with equal probability.
  broad <- sab_residual_ceiling(log(rep(c(2, 198), each = 2048)), 0)
  testthat::expect_gte(min(broad$collapse_upper), 0.01)
  testthat::expect_lt(min(broad$collapse_upper), 0.05)
})

testthat::test_that("exact zeros, clipping and malformed inputs are handled explicitly", {
  zero <- sab_residual_ceiling(numeric(), -Inf)
  testthat::expect_equal(zero$collapse_upper, rep(0, 3))
  clipped <- sab_residual_ceiling(c(-Inf, 10000), 0, caps = 10)
  testthat::expect_equal(clipped$mean_capped_ratio, 5)
  testthat::expect_error(sab_residual_ceiling(c(0, NA), 0), "log_ratio")
  testthat::expect_error(sab_residual_ceiling(Inf, 0), "log_ratio")
  testthat::expect_error(sab_residual_ceiling(0, Inf), "log_c_probe")
  testthat::expect_error(sab_residual_ceiling(numeric(), 0), "log_ratio")
  testthat::expect_error(sab_residual_ceiling(0, 0, caps = 0), "caps")
  testthat::expect_error(sab_residual_ceiling(0, 0, alpha_per_cap = 1), "alpha_per_cap")
})

testthat::test_that("one-dimensional truncated Gaussian integrates to one", {
  location <- 1.2; scale <- 2.3; mass <- 0.8
  factor <- matrix(scale, 1, 1)
  radius <- sqrt(stats::qchisq(mass, df = 1))
  density <- function(x) exp(sab_truncated_gaussian_logdensity(
    matrix(x, ncol = 1), location, factor, mass))
  integral <- stats::integrate(density, location - scale * radius,
                              location + scale * radius)$value
  testthat::expect_equal(integral, 1, tolerance = 1e-8)
  testthat::expect_identical(sab_truncated_gaussian_logdensity(
    location + scale * (radius + 0.01), location, factor, mass), -Inf)
})

testthat::test_that("eight-dimensional draws have the normalized ball radial law", {
  set.seed(8104)
  d <- 8L; n <- 4000L; mass <- 0.8; location <- rep(0.5, d)
  factor <- diag(seq(1, 2, length.out = d)); factor[2, 1] <- 0.3
  x <- sab_truncated_gaussian_draw(n, location, factor, mass)
  z <- t(forwardsolve(factor, t(sweep(x, 2, location, "-"))))
  radial_cdf <- stats::pchisq(rowSums(z^2), df = d)
  testthat::expect_identical(dim(x), c(n, d))
  testthat::expect_true(all(radial_cdf < mass))
  testthat::expect_equal(mean(radial_cdf), mass / 2, tolerance = 0.015)
  testthat::expect_equal(colMeans(z), rep(0, d), tolerance = 0.06)
  expected <- -0.5 * (d * log(2 * pi) + rowSums(z^2)) -
    sum(log(diag(factor))) - log(mass)
  testthat::expect_equal(sab_truncated_gaussian_logdensity(
    x, location, factor, mass), expected, tolerance = 1e-12)
})

testthat::test_that("Gaussian shape and support inputs fail closed", {
  testthat::expect_error(sab_truncated_gaussian_draw(1.5, 0, matrix(1), 0.8), "n")
  testthat::expect_error(sab_truncated_gaussian_draw(1, 0, matrix(0), 0.8), "Cholesky")
  testthat::expect_error(sab_truncated_gaussian_draw(1, 0, matrix(1), 1), "mass")
  upper <- diag(2); upper[1, 2] <- 0.1
  testthat::expect_error(sab_truncated_gaussian_draw(1, c(0, 0), upper, 0.8), "Cholesky")
  testthat::expect_error(sab_truncated_gaussian_logdensity(Inf, 0, matrix(1), 0.8), "x")
})
