source(file.path(core_root, "R", "shared_pool_kernel.R"), local = FALSE)

testthat::test_that("surrogate Gibbs excludes other occupied slots and divides by q", {
  p <- sab_shared_pool_assignment_probabilities(
    assigned = c(1L, 2L), patient_index = 1L,
    log_approx_likelihood = log(c(2, 1000, 1)),
    log_population = log(c(0.2, 0.4, 0.4)), log_q = log(c(0.2, 0.3, 0.5)))
  testthat::expect_identical(p$slots, c(1L, 3L))
  testthat::expect_equal(p$probabilities, c(2, 0.8) / 2.8)
  testthat::expect_identical(sab_shared_pool_reserve_slots(c(1L, 2L), 3L), 3L)
})

testthat::test_that("finite auxiliary target obeys detailed balance for nonuniform q", {
  # Six injective assignments of two patients to three fixed distinct slots.
  # Non-Gaussian finite densities and an inaccurate surrogate test both q and
  # exact-correction terms against the independently written enlarged target.
  states <- rbind(c(1L, 2L), c(1L, 3L), c(2L, 1L),
                  c(2L, 3L), c(3L, 1L), c(3L, 2L))
  q <- c(0.2, 0.3, 0.5)
  g <- rbind(c(0.1, 0.6, 0.3), c(0.5, 0.2, 0.3))
  exact <- rbind(c(0.9, 0.2, 0.7), c(0.3, 0.8, 0.4))
  approx <- rbind(c(0.2, 0.9, 0.3), c(0.8, 0.2, 0.5))
  target <- apply(states, 1L, function(a) {
    prod(exact[cbind(1:2, a)] * g[cbind(1:2, a)]) * q[setdiff(1:3, a)]
  })
  target <- target / sum(target)
  transition <- matrix(0, 6L, 6L)
  for (row in 1:6) for (i in 1:2) {
    a <- states[row, ]
    proposal <- sab_shared_pool_assignment_probabilities(
      a, i, log(approx[i, ]), log(g[i, ]), log(q))
    boundaries <- c(0, cumsum(proposal$probabilities))
    for (j in seq_along(proposal$slots)) {
      result <- sab_shared_pool_corrected_assignment(
        a, i, log(approx[i, ]), log(g[i, ]), log(q),
        exact_current = list(loglik = log(exact[i, a[[i]]])),
        evaluate_exact = function(slot) list(loglik = log(exact[i, slot])),
        selection_uniform = mean(boundaries[c(j, j + 1L)]),
        acceptance_uniform = 0.5)
      testthat::expect_identical(result$selected_slot, proposal$slots[[j]])
      b <- a; b[[i]] <- result$selected_slot
      next_row <- which(apply(states, 1L, function(s) identical(s, b)))
      mass <- 0.5 * proposal$probabilities[[j]]
      transition[row, next_row] <- transition[row, next_row] +
        mass * result$acceptance_probability
      transition[row, row] <- transition[row, row] +
        mass * (1 - result$acceptance_probability)
    }
  }
  flow <- sweep(transition, 1L, target, "*")
  testthat::expect_equal(rowSums(transition), rep(1, 6), tolerance = 1e-14)
  testthat::expect_equal(flow, t(flow), tolerance = 1e-14)
  testthat::expect_equal(as.numeric(target %*% transition), target, tolerance = 1e-14)
})

testthat::test_that("rejection retains active cache and acceptance replaces it atomically", {
  assigned <- c(1L, 2L)
  current <- list(loglik = 0, prediction = list(value = 17, token = "old"))
  candidate <- list(loglik = log(0.01), prediction = list(value = 29, token = "new"))
  call <- function(u) sab_shared_pool_corrected_assignment(
    assigned, 1L, rep(0, 3), rep(0, 3), rep(0, 3), current,
    function(slot) { stopifnot(slot == 3L); candidate },
    selection_uniform = 0.75, acceptance_uniform = u)
  rejected <- call(0.5)
  testthat::expect_false(rejected$accepted)
  testthat::expect_identical(rejected$assigned, assigned)
  testthat::expect_identical(rejected$exact_current, current)
  testthat::expect_identical(rejected$evaluated_exact, candidate)
  testthat::expect_equal(rejected$acceptance_probability, 0.01)
  accepted <- call(0.001)
  testthat::expect_true(accepted$moved)
  testthat::expect_identical(accepted$assigned, c(3L, 2L))
  testthat::expect_identical(accepted$exact_current, candidate)
  testthat::expect_identical(current$prediction$token, "old")
})

testthat::test_that("self assignments spend no exact callback and zero likelihood rejects", {
  never <- function(slot) stop("Unexpected exact evaluation")
  current <- list(loglik = 0, prediction = "kept")
  self <- sab_shared_pool_corrected_assignment(
    c(1L, 2L), 1L, rep(0, 3), rep(0, 3), rep(0, 3), current, never,
    selection_uniform = 0.1, acceptance_uniform = 0.5)
  testthat::expect_false(self$moved)
  testthat::expect_identical(self$exact_callback_calls, 0L)
  testthat::expect_identical(self$exact_current, current)
  zero <- sab_shared_pool_corrected_assignment(
    c(1L, 2L), 1L, rep(0, 3), rep(0, 3), rep(0, 3), current,
    function(slot) list(loglik = -Inf, prediction = "failed sealed solve"),
    selection_uniform = 0.75, acceptance_uniform = 0)
  testthat::expect_false(zero$moved)
  testthat::expect_equal(zero$acceptance_probability, 0)
  testthat::expect_identical(zero$exact_current, current)
})

testthat::test_that("fixed surrogate floor retains support and ratios remain stable", {
  floored <- sab_shared_pool_positive_loglik(c(-Inf, -10000, 2), -500)
  testthat::expect_equal(floored, c(-500, -500, 2))
  original <- sab_shared_pool_assignment_probabilities(
    c(1L, 2L), 1L, c(0, 100, -2), rep(0, 3), rep(0, 3))
  shifted <- sab_shared_pool_assignment_probabilities(
    c(1L, 2L), 1L, c(0, 100, -2) + 10000, rep(0, 3), rep(0, 3))
  testthat::expect_equal(shifted$probabilities, original$probabilities)
  same_target <- sab_shared_pool_corrected_assignment(
    c(1L, 2L), 1L, c(10000, 0, 9998), rep(0, 3), rep(0, 3),
    list(loglik = 10003), function(slot) list(loglik = 10001),
    selection_uniform = 0.99, acceptance_uniform = 0.999)
  testthat::expect_true(same_target$moved)
  testthat::expect_equal(same_target$log_acceptance_ratio, 0)
})

testthat::test_that("malformed assignments, q support and callback errors fail closed", {
  testthat::expect_error(sab_shared_pool_reserve_slots(c(1L, 1L), 3L), "unique")
  testthat::expect_error(sab_shared_pool_reserve_slots(c(1L, 4L), 3L), "valid")
  testthat::expect_error(sab_shared_pool_assignment_probabilities(
    c(1L, 2L), 1L, rep(0, 3), rep(0, 3), c(0, 0, -Inf)), "log_q")
  testthat::expect_error(sab_shared_pool_assignment_probabilities(
    c(1L, 2L), 1L, c(0, 0, -Inf), rep(0, 3), rep(0, 3)), "log_approx")
  testthat::expect_error(sab_shared_pool_corrected_assignment(
    c(1L, 2L), 1L, rep(0, 3), rep(0, 3), rep(0, 3), list(loglik = 0),
    function(slot) list(loglik = NaN), 0.75, 0.5), "evaluate_exact")
  testthat::expect_error(sab_shared_pool_corrected_assignment(
    c(1L, 2L), 1L, rep(0, 3), rep(0, 3), rep(0, 3), list(loglik = 0),
    function(slot) stop("unexpected solver defect"), 0.75, 0.5), "unexpected solver defect")
})
