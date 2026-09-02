testthat::test_that("patient transports sum their log Jacobians", {
  map_1 <- sab_new_affine_map(
    function(global) global,
    function(global) matrix(exp(0.2 * global), 1, 1),
    name = "patient_1"
  )
  map_2 <- sab_new_affine_map(
    function(global) c(global, -global),
    function(global) diag(c(exp(global), exp(-0.5 * global))),
    name = "patient_2"
  )
  maps <- list(patient_1 = map_1, patient_2 = map_2)
  states <- list(patient_1 = 0.3, patient_2 = c(-0.2, 0.7))

  result <- sab_transport_local_states(maps, states, -0.1, 0.4)
  expected <- sum(vapply(
    maps,
    sab_affine_transport_log_jacobian,
    numeric(1),
    from_global = -0.1,
    to_global = 0.4
  ))

  testthat::expect_equal(result$log_abs_det_jacobian, expected, tolerance = 1e-12)
  testthat::expect_equal(
    sab_transport_local_states(maps, result$local_states, 0.4, -0.1)$local_states,
    states,
    tolerance = 1e-12
  )
})

testthat::test_that("an exact Gaussian conditional map cancels local dependence", {
  map <- sab_new_affine_map(
    mean_fn = function(global) 1.5 * global,
    chol_fn = function(global) matrix(exp(0.25 * global), 1, 1)
  )
  exact_joint <- function(global, locals) {
    components <- sab_affine_map_components(map, global)
    u <- sab_affine_standardise(map, locals[[1]], global)
    stats::dnorm(global, log = TRUE) +
      stats::dnorm(u, log = TRUE) - components$log_abs_det
  }
  symmetric_proposal <- function(to, from) 0
  current_global <- -0.6
  proposed_global <- 0.35
  current_u <- 1.1
  current_x <- sab_affine_unstandardise(map, current_u, current_global)

  proposal <- sab_transport_mh_proposal(
    current_global = current_global,
    proposed_global = proposed_global,
    current_locals = list(current_x),
    maps = list(map),
    log_target = exact_joint,
    log_global_proposal = symmetric_proposal
  )

  expected <- stats::dnorm(proposed_global, log = TRUE) -
    stats::dnorm(current_global, log = TRUE)
  testthat::expect_equal(proposal$log_ratio, expected, tolerance = 1e-12)
  testthat::expect_equal(
    sab_affine_standardise(map, proposal$proposed_state$locals[[1]], proposed_global),
    current_u,
    tolerance = 1e-12
  )
})

testthat::test_that("a marginal independence proposal gives unit Gaussian ratio", {
  map <- sab_new_affine_map(
    mean_fn = function(global) global,
    chol_fn = function(global) matrix(exp(0.1 * global), 1, 1)
  )
  exact_joint <- function(global, locals) {
    u <- sab_affine_standardise(map, locals[[1]], global)
    stats::dnorm(global, log = TRUE) + stats::dnorm(u, log = TRUE) -
      sab_affine_log_abs_det(map, global)
  }
  marginal_independence <- function(to, from) stats::dnorm(to, log = TRUE)

  proposal <- sab_transport_mh_proposal(
    current_global = -0.2,
    proposed_global = 0.8,
    current_locals = list(sab_affine_unstandardise(map, -1.3, -0.2)),
    maps = list(map),
    log_target = exact_joint,
    log_global_proposal = marginal_independence
  )

  testthat::expect_equal(proposal$log_ratio, 0, tolerance = 1e-12)
  testthat::expect_equal(proposal$log_acceptance_probability, 0)
})

testthat::test_that("acceptance selects target evaluations atomically", {
  map <- sab_new_affine_map(
    mean_fn = function(global) global,
    chol_fn = function(global) matrix(1, 1, 1)
  )
  evaluator <- function(global, locals) {
    list(
      log_density = -0.5 * global^2 - 0.5 * locals[[1]]^2,
      cache = paste0("cache_at_", global)
    )
  }
  proposal <- sab_transport_mh_proposal(
    current_global = 0,
    proposed_global = 1,
    current_locals = list(0),
    maps = list(map),
    log_target = evaluator,
    log_global_proposal = function(to, from) 0
  )

  accepted <- sab_decide_transport_mh(proposal, log_uniform = -100)
  rejected <- sab_decide_transport_mh(proposal, log_uniform = -0.01)

  testthat::expect_true(accepted$accepted)
  testthat::expect_identical(accepted$evaluation$cache, "cache_at_1")
  testthat::expect_false(rejected$accepted)
  testthat::expect_identical(rejected$evaluation$cache, "cache_at_0")
  testthat::expect_equal(rejected$state$global, 0)
})

testthat::test_that("an invalid proposed target is rejected without NaN ratios", {
  map <- sab_new_affine_map(
    mean_fn = function(global) 0,
    chol_fn = function(global) matrix(1, 1, 1)
  )
  target <- function(global, locals) {
    if (global > 0) {
      return(list(log_density = -Inf, failure = "outside support"))
    }
    list(log_density = -0.5 * locals[[1]]^2, cache = "current")
  }
  proposal <- sab_transport_mh_proposal(
    current_global = -1,
    proposed_global = 1,
    current_locals = list(0),
    maps = list(map),
    log_target = target,
    log_global_proposal = function(to, from) 0
  )
  decision <- sab_decide_transport_mh(proposal, log_uniform = -Inf)

  testthat::expect_identical(proposal$log_ratio, -Inf)
  testthat::expect_false(decision$accepted)
  testthat::expect_identical(decision$evaluation$cache, "current")
})

testthat::test_that("patient names cannot be silently reordered", {
  map <- sab_new_affine_map(
    mean_fn = function(global) 0,
    chol_fn = function(global) matrix(1, 1, 1)
  )
  testthat::expect_error(
    sab_transport_local_states(
      maps = list(a = map, b = map),
      local_states = list(b = 0, a = 0),
      from_global = 0,
      to_global = 1
    ),
    "identical"
  )
})
