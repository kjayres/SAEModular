if (!exists("sab_system_a_contract", mode = "function")) {
  source(file.path(core_root, "R", "system_a_adapter.R"), local = FALSE)
}

testthat::test_that("System A provenance pins are complete and immutable", {
  contract <- sab_system_a_contract()

  testthat::expect_identical(
    contract$upstream_commit,
    "a3e0367b06c82ad4d07280c99f10d4f9bac69978"
  )
  testthat::expect_identical(
    contract$target_fingerprint,
    "123fc00a2c0151215f9c550613819d4563a02c5bbe789eee5f2aee2c20844925"
  )
  testthat::expect_identical(
    names(contract$source_sha256),
    c(
      "R/target_manifest.R",
      "R/gme_adapter.R",
      "R/exact_likelihood.R",
      "R/oracle_audit_contract.R",
      "R/oracle_audit_cases.R",
      "scripts/oracle_audit_run.R",
      "stan/patient_loglik_oracle.stan"
    )
  )
  testthat::expect_true(all(grepl(
    "^[0-9a-f]{64}$", unname(contract$source_sha256)
  )))
  testthat::expect_match(
    contract$oracle_certificate_sha256, "^[0-9a-f]{64}$"
  )
})

testthat::test_that("System A patient selection rejects ambiguous ordering", {
  canonical <- as.character(1:4)

  testthat::expect_identical(
    .sab_system_a_validate_patient_selection(NULL, canonical), canonical
  )
  testthat::expect_identical(
    .sab_system_a_validate_patient_selection(c("1", "3"), canonical),
    c("1", "3")
  )
  testthat::expect_error(
    .sab_system_a_validate_patient_selection(c("3", "1"), canonical),
    "canonical System A order",
    fixed = TRUE
  )
  testthat::expect_error(
    .sab_system_a_validate_patient_selection(c("1", "1"), canonical),
    "unique subset",
    fixed = TRUE
  )
  testthat::expect_error(
    .sab_system_a_validate_patient_selection("5", canonical),
    "unique subset",
    fixed = TRUE
  )
})

testthat::test_that("System A source validation fails closed", {
  empty_root <- tempfile("sab_missing_workspace_")
  dir.create(empty_root)
  on.exit(unlink(empty_root, recursive = TRUE), add = TRUE)

  testthat::expect_error(
    sab_validate_system_a_sources(empty_root),
    "pinned System A worktree does not exist",
    fixed = TRUE
  )
})

testthat::test_that("pinned System A sources and certificate validate", {
  workspace_root <- normalizePath(
    file.path(core_root, "..", ".."), mustWork = TRUE
  )
  expected_worktree <- file.path(
    workspace_root,
    sab_system_a_contract()$upstream_repository
  )
  testthat::skip_if_not(dir.exists(expected_worktree))

  validation <- sab_validate_system_a_sources(workspace_root)
  testthat::expect_s3_class(
    validation, "sab_system_a_source_validation"
  )
  testthat::expect_true(validation$upstream_clean)
  testthat::expect_identical(
    validation$upstream_commit,
    sab_system_a_contract()$upstream_commit
  )
  testthat::expect_identical(
    validation$source_sha256,
    sab_system_a_contract()$source_sha256
  )
  testthat::expect_identical(
    validation$canonical_source_sha256,
    sab_system_a_contract()$source_sha256
  )
})

testthat::test_that("certified System A callbacks load for a patient subset", {
  testthat::skip_if_not(
    identical(tolower(Sys.getenv("SAB_SYSTEM_A_INTEGRATION")), "true"),
    "Set SAB_SYSTEM_A_INTEGRATION=true in an HPC test job."
  )
  workspace_root <- normalizePath(
    file.path(core_root, "..", ".."), mustWork = TRUE
  )
  # These are the first twelve IDs in the hash-pinned patient manifest.  This
  # is an adapter smoke subset, not the eventual scientific pilot selection.
  pilot_ids <- as.character(2:13)
  adapter <- sab_load_system_a_adapter(workspace_root, pilot_ids)

  testthat::expect_s3_class(adapter, "sab_system_a_adapter")
  testthat::expect_invisible(sab_validate_system_a_adapter(adapter))
  testthat::expect_identical(adapter$patient_ids, pilot_ids)
  testthat::expect_length(adapter$coordinate_names$local, 8L)
  testthat::expect_length(adapter$coordinate_names$population, 17L)
  testthat::expect_identical(
    adapter$coordinate_names$dynamic_global,
    c("log_gamma_pop", "log_mu_l_pop")
  )
  testthat::expect_true(adapter$eta_in_domain(adapter$prior_reference$eta))
  testthat::expect_true(adapter$psi_in_domain(adapter$prior_reference$psi))

  patient_id <- adapter$patient_ids[[1L]]
  local_state <- adapter$population_mean(
    patient_id, adapter$prior_reference$eta
  )
  testthat::expect_true(is.finite(adapter$log_population_density(
    patient_id, local_state, adapter$prior_reference$eta
  )))
  testthat::expect_true(is.finite(adapter$log_hyperprior(
    adapter$prior_reference$eta, adapter$prior_reference$psi
  )))
  testthat::expect_error(
    adapter$log_likelihood("not-a-patient", local_state,
                           adapter$prior_reference$psi),
    "not in this System A adapter",
    fixed = TRUE
  )
})
