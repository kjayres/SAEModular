testthat::test_that("affine transport preserves u and has the stated inverse", {
  mean_fn <- function(global) {
    c(1 + global[[1]], -2 + 0.5 * global[[2]])
  }
  chol_fn <- function(global) {
    matrix(
      c(exp(0.2 * global[[1]]),
        0,
        0.1 * global[[2]], exp(-0.1 * global[[2]])),
      nrow = 2,
      byrow = TRUE
    )
  }
  map <- sab_new_affine_map(mean_fn, chol_fn, check_at = c(0, 0))

  from <- c(-0.4, 0.2)
  to <- c(0.7, -0.5)
  u <- c(1.2, -0.8)
  x <- sab_affine_unstandardise(map, u, from)
  transported <- sab_affine_transport(map, x, from, to)

  testthat::expect_equal(sab_affine_standardise(map, x, from), u, tolerance = 1e-12)
  testthat::expect_equal(
    sab_affine_standardise(map, transported, to),
    u,
    tolerance = 1e-12
  )
  testthat::expect_equal(
    sab_affine_inverse_transport(map, transported, from, to),
    x,
    tolerance = 1e-12
  )
})

testthat::test_that("reported affine Jacobian equals the matrix determinant", {
  map <- sab_new_affine_map(
    mean_fn = function(global) c(global, -global),
    chol_fn = function(global) {
      matrix(c(exp(global), 0, 0.25, exp(-0.3 * global)), 2, 2, byrow = TRUE)
    }
  )
  from <- -0.2
  to <- 0.45
  from_chol <- sab_affine_map_components(map, from)$chol
  to_chol <- sab_affine_map_components(map, to)$chol
  direct_matrix <- to_chol %*% solve(from_chol)

  testthat::expect_equal(
    sab_affine_transport_log_jacobian(map, from, to),
    as.numeric(determinant(direct_matrix, logarithm = TRUE)$modulus),
    tolerance = 1e-12
  )
  testthat::expect_equal(
    sab_affine_transport_log_jacobian(map, from, to),
    -sab_affine_transport_log_jacobian(map, to, from),
    tolerance = 1e-12
  )
})

testthat::test_that("transport preserves auditable coordinate names", {
  coordinate_names <- c("first", "second")
  map <- sab_new_affine_map(
    mean_fn = function(global) {
      stats::setNames(c(global, -global), coordinate_names)
    },
    chol_fn = function(global) diag(c(exp(global), exp(-global))),
    check_at = 0
  )
  x <- stats::setNames(c(0.2, -0.4), coordinate_names)

  transported <- sab_affine_transport(map, x, -0.1, 0.3)
  standardised <- sab_affine_standardise(map, x, -0.1)

  testthat::expect_identical(names(transported), coordinate_names)
  testthat::expect_identical(names(standardised), coordinate_names)
  testthat::expect_identical(
    names(sab_affine_inverse_transport(map, transported, -0.1, 0.3)),
    coordinate_names
  )
  testthat::expect_error(
    sab_affine_transport(map, rev(x), -0.1, 0.3),
    "must match the affine-map coordinates in order",
    fixed = TRUE
  )
})

testthat::test_that("a map cannot silently change coordinate names", {
  map <- sab_new_affine_map(
    mean_fn = function(global) {
      if (global < 0) c(first = 0, second = 0) else c(second = 0, first = 0)
    },
    chol_fn = function(global) diag(2),
    check_at = -1
  )

  testthat::expect_error(
    sab_affine_transport(map, c(first = 0, second = 0), -1, 1),
    "changed coordinate names",
    fixed = TRUE
  )
})

testthat::test_that("invalid affine maps fail before entering an MCMC kernel", {
  testthat::expect_error(
    sab_new_affine_map(function(global) c(0, 0), function(global) diag(c(1, -1)),
                   check_at = 0),
    "strictly positive"
  )
  testthat::expect_error(
    sab_new_affine_map(
      function(global) c(0, 0),
      function(global) matrix(c(1, 0.2, 0, 1), 2, 2, byrow = TRUE),
      check_at = 0
    ),
    "lower triangular"
  )
  testthat::expect_error(
    sab_new_affine_map(function(global) c(0, 0), function(global) diag(1),
                   check_at = 0),
    "invalid Cholesky"
  )
})
