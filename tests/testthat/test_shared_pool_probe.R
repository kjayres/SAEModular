source(file.path(core_root, "R", "shared_pool_probe.R"), local = FALSE)

testthat::test_that("reserve MH and Gibbs use the common q at both endpoints", {
  # Non-Gaussian/discrete positive densities are enough to test the algebra.
  result <- sab_shared_pool_probe(
    log_f_reserve = log(matrix(c(2, 1, 1, 4), nrow = 2L, byrow = TRUE)),
    log_q_reserve = log(c(0.5, 0.25)),
    log_f_current = log(c(2, 3)), log_q_current = log(c(0.25, 0.5)),
    useful_acceptance = 0.75,
    squared_jump = matrix(c(4, 1, 9, 4), nrow = 2L, byrow = TRUE)
  )
  expected_alpha <- matrix(c(0.5, 0.5, 1 / 3, 1), 2L, byrow = TRUE)
  expected_gibbs <- rbind(c(8, 4, 4) / 16, c(6, 2, 16) / 24)
  testthat::expect_equal(unname(result$acceptance), expected_alpha)
  testthat::expect_equal(unname(result$gibbs_probabilities), expected_gibbs)
  testthat::expect_equal(result$patient$expected_accepted_squared_jump,
                         c(1.25, 3.5))
  testthat::expect_equal(result$reserve$useful_patients, c(0, 1))
  testthat::expect_equal(
    result$summary$hypothetical_common_uniform_two_patient_probability,
    (1 / 3 + 0.5) / 2
  )
})

testthat::test_that("zero likelihood reserves remain in every denominator", {
  result <- sab_shared_pool_probe(
    matrix(c(0, -Inf, 0, -Inf), nrow = 2L, byrow = TRUE),
    c(0, 0), c(0, 0), c(0, 0)
  )
  testthat::expect_equal(unname(result$acceptance),
                         matrix(c(1, 0, 1, 0), 2L, byrow = TRUE))
  testthat::expect_equal(result$patient$reserve_count, c(2, 2))
  testthat::expect_equal(result$patient$mean_reserve_mh_acceptance, c(0.5, 0.5))
  testthat::expect_equal(result$patient$gibbs_move_probability, c(0.5, 0.5))
  testthat::expect_equal(result$summary$fraction_reserves_useful_to_multiple_patients,
                         0.5)
  failed <- sab_shared_pool_probe(matrix(-Inf, 1L, 3L), rep(0, 3), 0, 0)
  testthat::expect_equal(unname(failed$gibbs_probabilities),
                         matrix(c(1, 0, 0, 0), 1L))
  testthat::expect_equal(failed$patient$mean_reserve_mh_acceptance, 0)
})

testthat::test_that("unknown patient normalizers cancel at extreme log scales", {
  log_f <- matrix(c(1, 3, -Inf, -4, 0, 5), 2L, byrow = TRUE)
  current <- c(0, 2)
  original <- sab_shared_pool_probe(log_f, c(0, 1, -2), current, c(1, -1))
  shifted <- sab_shared_pool_probe(
    sweep(log_f, 1L, c(10000, -10000), "+"),
    c(0, 1, -2), current + c(10000, -10000), c(1, -1)
  )
  testthat::expect_equal(shifted$acceptance, original$acceptance)
  testthat::expect_equal(shifted$gibbs_probabilities, original$gibbs_probabilities)
  testthat::expect_equal(rowSums(shifted$gibbs_probabilities), c(`1` = 1, `2` = 1))
})

testthat::test_that("single reserve and patient retain matrix dimensions", {
  result <- sab_shared_pool_probe(matrix(log(2), 1L, 1L), 0, 0, 0)
  testthat::expect_equal(dim(result$acceptance), c(1L, 1L))
  testthat::expect_equal(unname(result$gibbs_probabilities),
                         matrix(c(1 / 3, 2 / 3), 1L))
  testthat::expect_equal(result$summary$fraction_reserves_useful_to_multiple_patients,
                         0)
  testthat::expect_equal(
    result$summary$hypothetical_common_uniform_two_patient_probability, 0
  )
})

testthat::test_that("support violations and invalid calculations fail closed", {
  testthat::expect_error(sab_shared_pool_probe(matrix(0, 1L), -Inf, 0, 0),
                         "log_q_reserve")
  testthat::expect_error(sab_shared_pool_probe(matrix(0, 1L), 0, -Inf, 0),
                         "log_f_current")
  testthat::expect_error(sab_shared_pool_probe(matrix(NaN, 1L), 0, 0, 0),
                         "log_f_reserve")
  testthat::expect_error(sab_shared_pool_probe(matrix(Inf, 1L), 0, 0, 0),
                         "log_f_reserve")
  testthat::expect_error(sab_shared_pool_probe(matrix(0, 1L), c(0, 0), 0, 0),
                         "wrong vector length")
  testthat::expect_error(sab_shared_pool_probe(matrix(0, 1L), 0, 0, 0,
                                             useful_acceptance = 0),
                         "useful_acceptance")
  testthat::expect_error(sab_shared_pool_probe(matrix(0, 1L), 0, 0, 0,
                                             squared_jump = matrix(-1, 1L)),
                         "squared_jump")
})

testthat::test_that("candidate likelihood alone would give a different answer", {
  # Equal likelihood/population f values do not cancel a nonuniform reserve q.
  result <- sab_shared_pool_probe(matrix(0, 1L), log(0.8), 0, log(0.2))
  testthat::expect_equal(unname(result$acceptance), matrix(0.25, 1L))
  testthat::expect_equal(unname(result$gibbs_probabilities),
                         matrix(c(0.8, 0.2), 1L))
})
