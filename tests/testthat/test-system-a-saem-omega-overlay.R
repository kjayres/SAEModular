if (!exists("sab_corrected_saem_omega_contract", mode = "function")) {
  source(
    file.path(core_root, "R", "system_a_saem_omega_overlay.R"),
    local = FALSE
  )
}

.sab_overlay_test_initials <- function() {
  list(omega = list(
    omega_lambda = 0.55,
    omega_mu_t = 0.44,
    omega_mu_a = 0.399,
    omega_p = 0.9,
    omega_alpha_l = 0.678,
    omega_pi = 0.45,
    omega_eta_rti = 2.93,
    omega_eta_pi = 3.19
  ))
}

testthat::test_that("corrected omega preserves published active starts", {
  contract <- sab_corrected_saem_omega_contract(
    .sab_overlay_test_initials()
  )
  expected <- c(
    u_log_lambda = 0.55^2,
    u_log_mu_t = 0.44^2,
    u_log_mu_a = 0.399^2,
    u_log_p = 0.9^2,
    u_log_alpha_l = 0.678^2,
    u_pi = 0.45^2,
    u_eta_rti = 2.93^2,
    u_eta_pi = 3.19^2
  )
  testthat::expect_equal(
    diag(contract$omega_init)[names(expected)], expected,
    tolerance = 0
  )
  testthat::expect_equal(
    diag(contract$omega_init)[contract$fixed_coordinates],
    stats::setNames(rep(1e-6, 4), contract$fixed_coordinates),
    tolerance = 0
  )
  testthat::expect_silent(chol(contract$omega_init))
  testthat::expect_gt(rcond(contract$omega_init), 1e-9)
  testthat::expect_lt(kappa(contract$omega_init, exact = TRUE), 2e7)
  testthat::expect_equal(
    unname(diag(contract$covariance_model)[contract$active_coordinates]),
    rep(1L, 8L)
  )
  testthat::expect_equal(
    unname(diag(contract$covariance_model)[contract$fixed_coordinates]),
    rep(0L, 4L)
  )
})

testthat::test_that("fitted structural omega must return to exact zeros", {
  contract <- sab_corrected_saem_omega_contract(
    .sab_overlay_test_initials()
  )
  fitted <- diag(rep(0, length(contract$coordinates)))
  dimnames(fitted) <- list(contract$coordinates, contract$coordinates)
  active_index <- match(contract$active_coordinates, contract$coordinates)
  fitted[cbind(active_index, active_index)] <- seq(0.1, 0.8,
                                                   length.out = 8)
  testthat::expect_silent(
    sab_assert_fitted_structural_omega(fitted, contract)
  )
  fitted[2, 2] <- 1e-6
  testthat::expect_error(
    sab_assert_fitted_structural_omega(fitted, contract),
    "structurally fixed"
  )
})

testthat::test_that("matrix assertion fails closed on altered inputs", {
  contract <- sab_corrected_saem_omega_contract(
    .sab_overlay_test_initials()
  )
  testthat::expect_silent(sab_assert_corrected_saem_matrices(
    contract$coordinates, contract$covariance_model,
    contract$omega_init, contract
  ))

  wrong <- contract$omega_init
  wrong[1, 1] <- 1
  testthat::expect_error(
    sab_assert_corrected_saem_matrices(
      contract$coordinates, contract$covariance_model, wrong, contract
    ),
    "differs"
  )
  singular <- contract$omega_init
  singular[2, 2] <- 0
  testthat::expect_error(
    sab_assert_corrected_saem_matrices(
      contract$coordinates, contract$covariance_model, singular, contract
    ),
    "positive definite"
  )
})

testthat::test_that("runner transformation is exact and fail-closed", {
  fixture <- c(
    'source(file.path(project_root, "R", "hiv_latent_saem_model.R"))',
    'run_dir <- file.path(project_root, "outputs", run_tag)',
    "    print = TRUE,",
    "      fit <- saemix::saemix(saem_model, saem_data, control = fit_control)"
  )
  transformed <- sab_build_corrected_saem_runner(fixture)
  testthat::expect_true(any(grepl(
    "system_a_saem_omega_overlay.R", transformed, fixed = TRUE
  )))
  testthat::expect_true(any(grepl(
    "SAEM_OUTPUT_ROOT", transformed, fixed = TRUE
  )))
  testthat::expect_true(any(grepl(
    "warnings = TRUE", transformed, fixed = TRUE
  )))
  testthat::expect_true(any(grepl(
    "sab_assert_fit_iteration0", transformed, fixed = TRUE
  )))
  testthat::expect_error(
    sab_build_corrected_saem_runner(fixture[-1]), "model-overlay"
  )
  testthat::expect_error(
    sab_build_corrected_saem_runner(c(fixture, fixture[1])),
    "model-overlay"
  )
})

testthat::test_that("the pinned upstream runner still matches the overlay", {
  runner <- normalizePath(
    file.path(core_root, "..", "saem_model", "scripts",
              "run_hiv_saem_fit.R"),
    mustWork = TRUE
  )
  testthat::expect_silent(
    parse(text = sab_build_corrected_saem_runner(
      readLines(runner, warn = TRUE)
    ))
  )
})
