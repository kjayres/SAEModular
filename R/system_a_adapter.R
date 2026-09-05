# Read-only, fail-closed bridge to the certified System A numerical target.
#
# This file contains no copy of the System A model.  It sources a hash-pinned
# implementation from a clean, commit-pinned worktree and lets that
# implementation reconstruct the immutable target from its existing manifest
# and oracle-audit certificate.  The canonical `projects/modular_bayes` mirror
# is checked as well because the certified constructor deliberately resolves
# its target artifacts and audit provenance there.

#' Frozen provenance contract for the System A adapter
#'
#' @return A list of commit, source, target, and certificate pins.
#' @export
sab_system_a_contract <- function() {
  source_sha256 <- c(
    "R/target_manifest.R" =
      "d4d2221245f7c215107d844b5faf4186436c5e53be0a4a3d4593dc9f294b1816",
    "R/gme_adapter.R" =
      "a6d61d490fe42463a2181e87a6f25e75a95f1dea09452e0ea88981fc697e5fef",
    "R/exact_likelihood.R" =
      "ea7931622ab3f888fba491b4d859c6678b2906d107a9d3c92423ad2c51072b80",
    "R/oracle_audit_contract.R" =
      "1e0766171c91b6ffb0f1021e081b6f07cd5da83349f61d0b390025e3590258ca",
    "R/oracle_audit_cases.R" =
      "f9bfff5ad732cfb94a5796ccf62c087eb0b969f715b5d1ec8db66cbbf96564d6",
    "scripts/oracle_audit_run.R" =
      "e785ca4cf6c78711968d4e6a7fdb5fddbe37960cfa445548bfdcd62c6dec3eec",
    "stan/patient_loglik_oracle.stan" =
      "d06d9e504446dc69711cc76aa8e7fc1cb58269bb58079a51d5d4f1ebe0652c24"
  )
  list(
    schema_version = "sab_system_a_contract_v1",
    upstream_commit =
      "a3e0367b06c82ad4d07280c99f10d4f9bac69978",
    upstream_repository = file.path(
      "projects", ".worktrees", "modular_bayes_system_a_factor_v1"
    ),
    system_a_source_root = file.path(
      "research", "system_a_modular_v0"
    ),
    canonical_source_root = file.path(
      "projects", "modular_bayes", "research", "system_a_modular_v0"
    ),
    sourced_files = c(
      "R/target_manifest.R",
      "R/gme_adapter.R",
      "R/exact_likelihood.R",
      "R/oracle_audit_contract.R"
    ),
    source_sha256 = source_sha256,
    target_fingerprint =
      "123fc00a2c0151215f9c550613819d4563a02c5bbe789eee5f2aee2c20844925",
    oracle_directory = file.path(
      "projects", "modular_bayes", "research", "system_a_modular_v0",
      "outputs", "production_v3", "oracle_audit"
    ),
    oracle_certificate_sha256 =
      "db341c06d16d891132f18d3f25d6c8f765d7e97f1b84a062e317c27d7c323deb"
  )
}

#' Validate every external source used by the System A adapter
#'
#' Validation is deliberately performed before any upstream file is sourced.
#' The complete upstream worktree must be at the pinned commit and clean.  Both
#' its System A sources and the canonical source mirror used by the certificate
#' must have the declared byte-level identities.
#'
#' @param workspace_root Explicit root containing `projects/`.
#'
#' @return A validated provenance record.
#' @export
sab_validate_system_a_sources <- function(workspace_root) {
  contract <- sab_system_a_contract()
  workspace_root <- .sab_system_a_normalize_directory(
    workspace_root, "workspace_root"
  )
  upstream_repository <- file.path(
    workspace_root, contract$upstream_repository
  )
  upstream_repository <- .sab_system_a_normalize_directory(
    upstream_repository, "pinned System A worktree"
  )

  git <- Sys.which("git")
  if (!nzchar(git)) {
    stop("git is required to validate the pinned System A worktree.",
         call. = FALSE)
  }
  head <- .sab_system_a_git(
    git, upstream_repository, c("rev-parse", "--verify", "HEAD")
  )
  if (length(head) != 1L || !identical(head, contract$upstream_commit)) {
    stop(
      "System A worktree is at commit ", paste(head, collapse = " "),
      "; expected ", contract$upstream_commit, ".",
      call. = FALSE
    )
  }
  status <- .sab_system_a_git(
    git,
    upstream_repository,
    c("status", "--porcelain=v1", "--untracked-files=all")
  )
  if (length(status)) {
    stop(
      "Pinned System A worktree is not clean; refusing to load it: ",
      paste(status, collapse = " | "),
      call. = FALSE
    )
  }

  upstream_source_root <- file.path(
    upstream_repository, contract$system_a_source_root
  )
  upstream_source_root <- .sab_system_a_normalize_directory(
    upstream_source_root, "pinned System A source root"
  )
  canonical_source_root <- file.path(
    workspace_root, contract$canonical_source_root
  )
  canonical_source_root <- .sab_system_a_normalize_directory(
    canonical_source_root, "canonical System A source root"
  )

  upstream_hashes <- .sab_system_a_check_hashes(
    upstream_source_root,
    contract$source_sha256,
    "pinned worktree"
  )
  canonical_hashes <- .sab_system_a_check_hashes(
    canonical_source_root,
    contract$source_sha256,
    "canonical target mirror"
  )

  oracle_directory <- file.path(workspace_root, contract$oracle_directory)
  oracle_directory <- .sab_system_a_normalize_directory(
    oracle_directory, "certified System A oracle directory"
  )
  certificate_path <- file.path(oracle_directory, "certificate.rds")
  certificate_hash <- .sab_system_a_sha256_file(certificate_path)
  if (!identical(
      certificate_hash, contract$oracle_certificate_sha256
  )) {
    stop(
      "System A oracle certificate hash mismatch; expected ",
      contract$oracle_certificate_sha256, ", observed ", certificate_hash,
      ".",
      call. = FALSE
    )
  }

  structure(
    list(
      schema_version = "sab_system_a_source_validation_v1",
      workspace_root = workspace_root,
      upstream_repository = upstream_repository,
      upstream_commit = head,
      upstream_clean = TRUE,
      upstream_source_root = upstream_source_root,
      canonical_source_root = canonical_source_root,
      source_sha256 = upstream_hashes,
      canonical_source_sha256 = canonical_hashes,
      oracle_directory = oracle_directory,
      oracle_certificate_sha256 = certificate_hash,
      target_fingerprint = contract$target_fingerprint
    ),
    class = c("sab_system_a_source_validation", "list")
  )
}

#' Load the certified System A callbacks for an exact pilot
#'
#' The full System A manifest and VODE-BDF likelihood remain owned by the
#' upstream project.  This function exposes only a selected, canonically
#' ordered patient view and the exact callbacks required by the affine pilot.
#' No target artifact or upstream source is written.
#'
#' @param workspace_root Explicit root containing `projects/`.
#' @param patient_ids Optional character vector of patients in canonical order.
#'   `NULL` exposes all 115 patients.  A 12-patient pilot should pass its
#'   predeclared IDs explicitly.
#'
#' @return A validated object of class `sab_system_a_adapter`.
#' @export
sab_load_system_a_adapter <- function(workspace_root, patient_ids = NULL) {
  validation <- sab_validate_system_a_sources(workspace_root)
  contract <- sab_system_a_contract()

  upstream <- new.env(parent = globalenv())
  for (relative_path in contract$sourced_files) {
    sys.source(
      file.path(validation$upstream_source_root, relative_path),
      envir = upstream,
      keep.source = FALSE
    )
  }
  required_functions <- c(
    "sysa_target_manifest",
    "sysa_oracle_audit_load_attached_likelihood",
    "sysa_validate_exact_likelihood_object",
    "sysa_local_coordinate_names",
    "sysa_eta_names",
    "sysa_psi_names",
    "sysa_population_mean",
    "sysa_log_population_density",
    "sysa_log_hyperprior",
    "sysa_eta_in_domain",
    "sysa_psi_in_domain",
    "sysa_default_psi_reference"
  )
  missing_functions <- required_functions[!vapply(
    required_functions,
    exists,
    logical(1L),
    envir = upstream,
    mode = "function",
    inherits = FALSE
  )]
  if (length(missing_functions)) {
    stop(
      "Pinned System A sources omit required functions: ",
      paste(missing_functions, collapse = ", "), ".",
      call. = FALSE
    )
  }

  manifest <- upstream$sysa_target_manifest(validation$workspace_root)
  if (!identical(manifest$target_fingerprint, contract$target_fingerprint)) {
    stop("Loaded System A manifest has the wrong target fingerprint.",
         call. = FALSE)
  }
  exact_likelihood <- upstream$sysa_oracle_audit_load_attached_likelihood(
    manifest = manifest,
    output_dir = validation$oracle_directory,
    solver = "vode_bdf"
  )
  upstream$sysa_validate_exact_likelihood_object(
    exact_likelihood,
    manifest = manifest,
    require_numerical_audit = TRUE
  )
  lockEnvironment(upstream, bindings = TRUE)

  canonical_ids <- as.character(manifest$patient_manifest$patient_id)
  selected_ids <- .sab_system_a_validate_patient_selection(
    patient_ids, canonical_ids
  )
  patient_rows <- match(selected_ids, canonical_ids)
  patient_manifest <- manifest$patient_manifest[
    patient_rows,
    c("patient_id", "treatment", "treatment_label", "treat_nelf", "n_obs"),
    drop = FALSE
  ]
  patient_manifest$patient_id <- as.character(patient_manifest$patient_id)
  rownames(patient_manifest) <- NULL
  treatment <- setNames(
    as.integer(patient_manifest$treat_nelf),
    patient_manifest$patient_id
  )

  callback_environment <- new.env(parent = baseenv())
  callback_environment$exact_likelihood <- exact_likelihood
  callback_environment$upstream <- upstream
  callback_environment$selected_ids <- selected_ids
  callback_environment$treatment <- treatment
  callback_environment$stan_data <- manifest$processed_data
  callbacks <- evalq({
    validate_patient_id <- function(patient_id) {
      patient_id <- as.character(patient_id)
      if (length(patient_id) != 1L || is.na(patient_id) ||
          !patient_id %in% selected_ids) {
        stop("patient_id is not in this System A adapter.", call. = FALSE)
      }
      patient_id
    }
    list(
      log_likelihood = function(patient_id, x, psi) {
        patient_id <- validate_patient_id(patient_id)
        exact_likelihood$log_likelihood(patient_id, x, psi)
      },
      solve_prediction = function(patient_id, x, psi) {
        patient_id <- validate_patient_id(patient_id)
        exact_likelihood$solve_prediction(patient_id, x, psi)
      },
      loglik_from_prediction = function(prediction, psi) {
        exact_likelihood$loglik_from_prediction(prediction, psi)
      },
      population_mean = function(patient_id, eta) {
        patient_id <- validate_patient_id(patient_id)
        upstream$sysa_population_mean(eta, treatment[[patient_id]])
      },
      log_population_density = function(patient_id, x, eta) {
        patient_id <- validate_patient_id(patient_id)
        upstream$sysa_log_population_density(
          x, eta, treatment[[patient_id]]
        )
      },
      log_hyperprior = function(eta, psi) {
        upstream$sysa_log_hyperprior(eta, psi, stan_data)
      },
      eta_in_domain = function(eta) {
        upstream$sysa_eta_in_domain(eta, stan_data)
      },
      psi_in_domain = function(psi) {
        upstream$sysa_psi_in_domain(psi, stan_data)
      }
    )
  }, envir = callback_environment)
  # These callbacks close over the immutable resolved Stan data in a private
  # locked environment, so callers cannot rebind the target objects.
  lockEnvironment(callback_environment, bindings = TRUE)

  certificate <- exact_likelihood$numerical_audit_certificate
  adapter <- structure(
    c(
      list(
        schema_version = "sab_system_a_adapter_v1",
        target_fingerprint = manifest$target_fingerprint,
        likelihood_signature = exact_likelihood$likelihood_signature,
        numerical_target = "sealed_deSolve_VODE_BDF",
        patient_ids = selected_ids,
        patient_manifest = patient_manifest,
        treatment = treatment,
        coordinate_names = list(
          local = upstream$sysa_local_coordinate_names(),
          population = upstream$sysa_eta_names(),
          global = upstream$sysa_psi_names(),
          dynamic_global = exact_likelihood$dynamic_psi_names,
          observation_global = exact_likelihood$observation_psi_names
        ),
        prior_reference = list(
          eta = manifest$priors$eta_normal_mean,
          eta_sd = manifest$priors$eta_normal_sd,
          psi = upstream$sysa_default_psi_reference(
            manifest$processed_data
          )
        ),
        oracle_audit = list(
          status = exact_likelihood$numerical_equivalence_status,
          certificate_sha256 = validation$oracle_certificate_sha256,
          audit_contract_fingerprint =
            certificate$audit_contract_fingerprint,
          max_abs_loglik_error = certificate$max_abs_loglik_error,
          declared_max_abs_loglik_error =
            certificate$declared_max_abs_loglik_error,
          desolve_version = certificate$desolve_version,
          r_version = certificate$r_version,
          platform = certificate$platform
        ),
        source_validation = validation
      ),
      callbacks
    ),
    class = c("sab_system_a_adapter", "list")
  )
  sab_validate_system_a_adapter(adapter)
  adapter
}

#' Validate a loaded System A adapter
#'
#' @param adapter Object returned by [sab_load_system_a_adapter()].
#'
#' @return `adapter`, invisibly.
#' @export
sab_validate_system_a_adapter <- function(adapter) {
  required_callbacks <- c(
    "log_likelihood", "solve_prediction", "loglik_from_prediction",
    "population_mean", "log_population_density", "log_hyperprior",
    "eta_in_domain", "psi_in_domain"
  )
  required_coordinates <- c(
    "local", "population", "global", "dynamic_global", "observation_global"
  )
  valid <- inherits(adapter, "sab_system_a_adapter") &&
    identical(adapter$schema_version, "sab_system_a_adapter_v1") &&
    identical(
      adapter$target_fingerprint,
      sab_system_a_contract()$target_fingerprint
    ) &&
    is.character(adapter$patient_ids) && length(adapter$patient_ids) >= 1L &&
    !anyNA(adapter$patient_ids) && !anyDuplicated(adapter$patient_ids) &&
    is.data.frame(adapter$patient_manifest) &&
    identical(adapter$patient_manifest$patient_id, adapter$patient_ids) &&
    is.integer(adapter$treatment) &&
    identical(names(adapter$treatment), adapter$patient_ids) &&
    all(adapter$treatment %in% c(0L, 1L)) &&
    is.list(adapter$coordinate_names) &&
    identical(names(adapter$coordinate_names), required_coordinates) &&
    identical(adapter$coordinate_names$dynamic_global,
              c("log_gamma_pop", "log_mu_l_pop")) &&
    identical(adapter$coordinate_names$observation_global,
              c("log_sigma_v", "log_sigma_t")) &&
    is.list(adapter$oracle_audit) &&
    identical(
      adapter$oracle_audit$status, "stan_bdf_oracle_audit_passed"
    ) &&
    is.list(adapter$source_validation) &&
    inherits(
      adapter$source_validation, "sab_system_a_source_validation"
    ) &&
    all(vapply(required_callbacks, function(name) {
      is.function(adapter[[name]])
    }, logical(1L)))
  if (!isTRUE(valid)) {
    stop("Malformed or uncertified System A adapter.", call. = FALSE)
  }
  invisible(adapter)
}

.sab_system_a_validate_patient_selection <- function(patient_ids,
                                                      canonical_ids) {
  if (is.null(patient_ids)) return(canonical_ids)
  if (!is.character(patient_ids) || !length(patient_ids) ||
      anyNA(patient_ids) || any(!nzchar(patient_ids)) ||
      anyDuplicated(patient_ids) || length(setdiff(patient_ids, canonical_ids))) {
    stop("patient_ids must be a non-empty unique subset of System A IDs.",
         call. = FALSE)
  }
  canonical_selection <- canonical_ids[canonical_ids %in% patient_ids]
  if (!identical(patient_ids, canonical_selection)) {
    stop("patient_ids must be supplied in canonical System A order.",
         call. = FALSE)
  }
  patient_ids
}

.sab_system_a_normalize_directory <- function(path, label) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !nzchar(path) || !dir.exists(path)) {
    stop(label, " does not exist or is not one directory.", call. = FALSE)
  }
  normalizePath(path, mustWork = TRUE)
}

.sab_system_a_git <- function(git, repository, arguments) {
  output <- suppressWarnings(system2(
    git,
    args = c("-C", shQuote(repository), arguments),
    stdout = TRUE,
    stderr = TRUE
  ))
  exit_status <- attr(output, "status")
  if (!is.null(exit_status) && exit_status != 0L) {
    stop(
      "git failed while validating the System A worktree: ",
      paste(output, collapse = " | "),
      call. = FALSE
    )
  }
  output <- trimws(output)
  output[nzchar(output)]
}

.sab_system_a_check_hashes <- function(root, expected, label) {
  observed <- vapply(names(expected), function(relative_path) {
    path <- file.path(root, relative_path)
    if (!file.exists(path)) {
      stop(label, " source is absent: ", path, call. = FALSE)
    }
    .sab_system_a_sha256_file(path)
  }, character(1L), USE.NAMES = TRUE)
  mismatch <- names(expected)[observed != expected]
  if (length(mismatch)) {
    details <- vapply(mismatch, function(relative_path) {
      paste0(
        relative_path, " expected ", expected[[relative_path]],
        " observed ", observed[[relative_path]]
      )
    }, character(1L))
    stop(
      label, " source hash mismatch: ", paste(details, collapse = "; "),
      call. = FALSE
    )
  }
  observed
}

.sab_system_a_sha256_file <- function(path) {
  if (!is.character(path) || length(path) != 1L || is.na(path) ||
      !file.exists(path)) {
    stop("Cannot hash missing file: ", path, call. = FALSE)
  }
  path <- normalizePath(path, mustWork = TRUE)
  sha256sum <- Sys.which("sha256sum")
  if (nzchar(sha256sum)) {
    output <- suppressWarnings(system2(
      sha256sum,
      args = c("--", shQuote(path)),
      stdout = TRUE,
      stderr = TRUE
    ))
  } else {
    shasum <- Sys.which("shasum")
    if (!nzchar(shasum)) {
      stop("Neither sha256sum nor shasum is available.", call. = FALSE)
    }
    output <- suppressWarnings(system2(
      shasum,
      args = c("-a", "256", shQuote(path)),
      stdout = TRUE,
      stderr = TRUE
    ))
  }
  exit_status <- attr(output, "status")
  if ((!is.null(exit_status) && exit_status != 0L) || length(output) != 1L) {
    stop("SHA-256 calculation failed for ", path, ".", call. = FALSE)
  }
  digest <- tolower(strsplit(trimws(output), "[[:space:]]+")[[1L]][[1L]])
  if (!grepl("^[0-9a-f]{64}$", digest)) {
    stop("SHA-256 tool returned an invalid digest for ", path, ".",
         call. = FALSE)
  }
  digest
}
