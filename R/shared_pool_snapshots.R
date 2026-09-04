# Two predeclared, coherent joint draws from the Bayesian System A comparator.
# This supplies an oracle feasibility design, not deployable initialization or
# independent conditional banks. The caller must reevaluate the sealed target.

sab_shared_pool_snapshots <- function(workspace_root, adapter) {
  sab_validate_system_a_adapter(adapter)
  if (length(adapter$patient_ids) != 115L) {
    stop("Shared-pool snapshots require the full 115-patient adapter.",
         call. = FALSE)
  }
  run <- "hiv_current_1000w1000s_20260506_120805"
  directory <- file.path(workspace_root, "projects", "full_joint_model",
                         "stan", "outputs", run)
  chains <- c(1L, 3L)
  hashes <- c(
    "18b96caccb71b4dd5391286889286451e7b5c075622993355bddfe6e283e54d7",
    "46240f79d182078869e9af802e214f2ae97327b6baa0f32b00b92a61ef41a198"
  )
  patient_hash <-
    "f3e773f466844a5ed9742a97c26b8b9c3389fbc2371bd2c48f2b04e8b017fab5"
  locations <- c("mu_log_lambda", "mu_log_mu_t", "mu_log_mu_a", "mu_log_p",
                 "mu_log_alpha_l", "mu_u_pi", "mu_u_eta_rti", "mu_u_eta_pi")
  scales <- c("omega_lambda", "omega_mu_t", "omega_mu_a", "omega_p",
              "omega_alpha_l", "omega_pi", "omega_eta_rti", "omega_eta_pi")
  globals <- c("gamma_pop", "mu_l_pop", "sigma_v", "sigma_t")
  locals <- c("lambda", "mu_t", "mu_a", "p", "alpha_l", "pi",
              "eta_rti", "eta_pi")
  snapshots <- lapply(seq_along(chains), function(k) {
    stem <- sprintf("%s_chain%02d", run, chains[[k]])
    patient_path <- file.path(directory, paste0(stem, "_patients.csv"))
    path <- file.path(directory, paste0(stem, "-1.csv"))
    if (!identical(.sab_system_a_sha256_file(patient_path), patient_hash) ||
        !identical(.sab_system_a_sha256_file(path), hashes[[k]])) {
      stop("Shared-pool comparator artifact hash mismatch.", call. = FALSE)
    }
    patients <- utils::read.csv(patient_path, stringsAsFactors = FALSE)
    patients$patient_id <- as.character(patients$patient_id)
    manifest <- adapter$patient_manifest
    fields <- c("patient_id", "treatment", "treatment_label", "treat_nelf",
                "n_obs")
    if (nrow(patients) != 115L || anyDuplicated(patients$patient_id) ||
        !all(fields %in% names(patients)) ||
        !setequal(patients$patient_id, adapter$patient_ids)) {
      stop("Comparator patient manifest is malformed.", call. = FALSE)
    }
    # Stan's suffix is row index in its patient manifest, never patient ID.
    index <- match(adapter$patient_ids, patients$patient_id)
    for (field in fields) {
      if (!identical(as.character(patients[[field]][index]),
                     as.character(manifest[[field]]))) {
        stop("Comparator/target patient manifest mismatch: ", field,
             call. = FALSE)
      }
    }
    local_fields <- unlist(lapply(index, function(i) paste0(locals, ".", i)),
                           use.names = FALSE)
    selected <- c(locations, "beta_nelf", scales, globals, local_fields)
    draw <- .sab_shared_pool_read_joint_row(path, selected, 1750L, 2000L)
    if (any(draw[c(scales, globals)] <= 0)) {
      stop("Comparator has nonpositive population/shared scales.",
           call. = FALSE)
    }
    eta <- setNames(c(draw[locations], draw["beta_nelf"], log(draw[scales])),
                    adapter$coordinate_names$population)
    psi <- setNames(log(draw[globals]), adapter$coordinate_names$global)
    natural <- matrix(draw[local_fields], nrow = 115L, byrow = TRUE)
    if (any(natural[, 1:5] <= 0) || any(natural[, 6:8] <= 0) ||
        any(natural[, 6:8] >= 1)) {
      stop("Comparator local draw is outside its transformed domain.",
           call. = FALSE)
    }
    x <- cbind(log(natural[, 1:5]), stats::qlogis(natural[, 6:8]))
    dimnames(x) <- list(adapter$patient_ids, adapter$coordinate_names$local)
    if (!isTRUE(adapter$eta_in_domain(eta)) ||
        !isTRUE(adapter$psi_in_domain(psi)) || any(!is.finite(x))) {
      stop("Comparator snapshot is outside the System A domain.",
           call. = FALSE)
    }
    list(
      id = sprintf("comparator_chain%02d_retained0750", chains[[k]]),
      eta = eta, psi = psi, x = x,
      source = list(path = normalizePath(path), sha256 = hashes[[k]],
                    chain = chains[[k]], retained_row = 750L,
                    warmup_rows = 1000L, numeric_row = 1750L,
                    patient_manifest_path = normalizePath(patient_path),
                    patient_manifest_sha256 = patient_hash,
                    patient_index = setNames(index, adapter$patient_ids),
                    target_fingerprint = adapter$target_fingerprint)
    )
  })
  names(snapshots) <- vapply(snapshots, `[[`, character(1L), "id")
  snapshots
}

.sab_shared_pool_read_joint_row <- function(path, selected, row, expected_rows) {
  connection <- file(path, open = "rt")
  on.exit(close(connection), add = TRUE)
  header <- NULL
  count <- 0L
  result <- NULL
  repeat {
    lines <- readLines(connection, n = 256L, warn = FALSE)
    if (!length(lines)) break
    lines <- lines[nzchar(lines) & !startsWith(lines, "#")]
    for (line in lines) {
      fields <- strsplit(line, ",", fixed = TRUE)[[1L]]
      if (is.null(header)) {
        header <- fields
        if (anyDuplicated(header) || !all(selected %in% header)) {
          stop("Comparator CSV header is incompatible.", call. = FALSE)
        }
      } else {
        count <- count + 1L
        if (length(fields) != length(header)) {
          stop("Comparator CSV row width changed.", call. = FALSE)
        }
        if (count == row) {
          result <- setNames(as.numeric(fields[match(selected, header)]),
                             selected)
        }
      }
    }
  }
  if (count != expected_rows || is.null(result) || any(!is.finite(result))) {
    stop("Comparator retained row/count is invalid.", call. = FALSE)
  }
  result
}
