#!/usr/bin/env Rscript
# Bounded oracle feasibility only. Run through the associated Slurm launcher.
options(warn = 1)
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Usage: system_a_shared_pool_probe.R <new-output-dir>")
project <- normalizePath(Sys.getenv("SAB_PROJECT_ROOT"), mustWork = TRUE)
workspace <- normalizePath(file.path(project, "..", ".."), mustWork = TRUE)
output <- args[[1L]]
if (file.exists(output) || !dir.create(output, recursive = TRUE)) {
  stop("Refusing existing or uncreatable output directory: ", output)
}
modules <- file.path(project, "R", c("system_a_adapter.R", "shared_pool_snapshots.R",
                                      "shared_trajectory_probe.R", "shared_pool_probe.R"))
invisible(lapply(modules, sys.source, envir = globalenv()))
started <- proc.time()[["elapsed"]]
ids <- c("3", "6", "16", "20", "55", "71", "74", "88", "105", "111", "117", "122")
full <- sab_load_system_a_adapter(workspace)
panel <- sab_load_system_a_adapter(workspace, ids)
snapshots <- sab_shared_pool_snapshots(workspace, full)
probes <- list(panel12 = sab_make_shared_trajectory_probe(panel),
               grid115 = sab_make_shared_trajectory_probe(full))
cores <- min(4L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
stopifnot(is.finite(cores), cores >= 1L)
representatives <- vapply(0:1, function(t) full$patient_ids[which(full$treatment == t)[1]], "")
components <- unlist(lapply(snapshots, function(s) lapply(representatives, function(id) {
  list(eta = s$eta, patient_id = id, mean = full$population_mean(id, s$eta),
       sd = setNames(exp(s$eta[10:17]), full$coordinate_names$local))
})), recursive = FALSE)
log_q <- function(x) {
  values <- vapply(components, function(component) full$log_population_density(
    component$patient_id, x, component$eta), numeric(1L))
  maximum <- max(values)
  maximum + log(mean(exp(values - maximum)))
}
set.seed(20260904)
reserve <- t(vapply(seq_len(128L), function(j) {
  component <- components[[sample.int(length(components), 1L)]]
  component$mean + component$sd * rnorm(8L)
}, numeric(8L)))
dimnames(reserve) <- list(sprintf("reserve_%03d", 1:128), full$coordinate_names$local)
source_paths <- c(modules, file.path(project, "experiments", "system_a_shared_pool_probe.R"),
                  file.path(project, "SHARED_POOL_CONTRACT.md"))
write.csv(data.frame(path = source_paths, sha256 = vapply(
  source_paths, .sab_system_a_sha256_file, "")), file.path(output, "sources.csv"), row.names = FALSE)
saveRDS(list(snapshots = snapshots, reserve = reserve, components = components,
             reference_target = full$target_fingerprint,
             experimental_targets = lapply(probes, function(p) p[c(
               "target_fingerprint", "positive_times", "union_times", "controls")]),
             contract = "SHARED_POOL_CONTRACT.md", seed = 20260904),
        file.path(output, "design.rds"))
write.csv(do.call(rbind, lapply(snapshots, function(s) data.frame(
  snapshot = s$id, parameter = c(names(s$eta), names(s$psi)),
  value = c(s$eta, s$psi)))), file.path(output, "contexts.csv"), row.names = FALSE)

evaluate_vectors <- function(probe, vectors, snapshot, grid) {
  cat("Evaluating", snapshot$id, grid, nrow(vectors), "vectors on", cores, "workers\n")
  flush.console()
  values <- parallel::mclapply(seq_len(nrow(vectors)), function(j) {
    result <- sab_shared_trajectory_evaluate(
      probe, vectors[j, ], snapshot$psi, common_first = (j %% 2L == 0L))
    result$comparison$snapshot <- snapshot$id
    result$comparison$grid <- grid
    result$comparison$vector_id <- rownames(vectors)[[j]]
    result$ledger$snapshot <- snapshot$id
    result$ledger$grid <- grid
    result$ledger$vector_id <- rownames(vectors)[[j]]
    list(original = result$original_loglik, common = result$common_loglik,
         comparison = result$comparison, ledger = result$ledger)
  }, mc.cores = cores, mc.preschedule = TRUE, mc.set.seed = FALSE)
  errors <- vapply(values, inherits, logical(1L), "try-error")
  if (any(errors)) stop("Vector worker failed: ", paste(values[errors], collapse = "\n"))
  list(original = do.call(cbind, lapply(values, `[[`, "original")),
       common = do.call(cbind, lapply(values, `[[`, "common")),
       comparison = do.call(rbind, lapply(values, `[[`, "comparison")),
       ledger = do.call(rbind, lapply(values, `[[`, "ledger")))
}
wilson <- function(successes, n) {
  z <- qnorm(.975); p <- successes/n
  centre <- (p + z*z/(2*n))/(1 + z*z/n)
  half <- z * sqrt(p*(1-p)/n + z*z/(4*n*n))/(1 + z*z/n)
  c(lower = centre-half, upper = centre+half)
}
maximum_finite <- function(x) if (any(is.finite(x))) max(abs(x[is.finite(x)])) else NA_real_
all_comparisons <- all_ledgers <- all_patients <- all_reserves <- all_moves <- list()
decisions <- timings <- list()
for (k in seq_along(snapshots)) {
  snapshot <- snapshots[[k]]
  donor <- snapshot$x[ids, , drop = FALSE]
  rownames(donor) <- paste0("patient_", ids)
  vectors <- rbind(reserve, donor)
  paired <- evaluate_vectors(probes$panel12, vectors, snapshot, "panel12")
  selected <- c(1:4, 128L + match(c("3", "20", "74", "122"), ids))
  large <- evaluate_vectors(probes$grid115, vectors[selected, , drop = FALSE],
                            snapshot, "grid115")
  saveRDS(list(panel = paired, grid115 = large),
          file.path(output, paste0(snapshot$id, "_paired.rds")))
  # Current likelihoods are the diagonal of the paired current-patient block.
  current_columns <- 128L + seq_along(ids)
  current_ll <- paired$original[cbind(seq_along(ids), current_columns)]
  current_common <- paired$common[cbind(seq_along(ids), current_columns)]
  if (any(!is.finite(current_ll)) || any(!is.finite(current_common))) {
    stop("A predeclared joint current state is invalid under a tested interface; see paired.rds.")
  }
  population <- t(vapply(ids, function(id) apply(vectors, 1L, function(x) {
    panel$log_population_density(id, x, snapshot$eta)
  }), numeric(nrow(vectors))))
  logq <- apply(vectors, 1L, log_q)
  current_g <- population[cbind(seq_along(ids), current_columns)]
  current_q <- logq[current_columns]
  scale <- exp(snapshot$eta[10:17])
  jump <- t(vapply(seq_along(ids), function(i) rowSums(sweep(
    sweep(vectors, 2L, donor[i, ], "-"), 2L, scale, "/")^2), numeric(nrow(vectors))))
  result <- sab_shared_pool_probe(
    paired$original + population, logq, current_ll + current_g, current_q,
    squared_jump = jump)
  common_result <- sab_shared_pool_probe(
    paired$common + population, logq, current_common + current_g, current_q,
    squared_jump = jump)
  # Primary evidence uses ONLY actual direct q draws, not posterior donor states.
  direct <- sab_shared_pool_probe(
    (paired$original + population)[, 1:128, drop = FALSE], logq[1:128],
    current_ll + current_g, current_q, squared_jump = jump[, 1:128, drop = FALSE])
  alpha <- result$acceptance
  alpha_common <- common_result$acceptance
  colnames(alpha) <- colnames(alpha_common) <- rownames(vectors)
  rownames(alpha) <- rownames(alpha_common) <- ids
  patient <- direct$patient; patient$patient_id <- ids; patient$snapshot <- snapshot$id
  patient$treat_nelf <- unname(panel$treatment)
  patient$n_obs <- panel$patient_manifest$n_obs
  reserve_summary <- direct$reserve; reserve_summary$snapshot <- snapshot$id
  reserve_summary$reserve_id <- rownames(reserve)
  moves <- expand.grid(patient_id = ids, vector_id = rownames(vectors), stringsAsFactors = FALSE)
  moves$snapshot <- snapshot$id
  moves$source <- rep(c(rep("direct_q", 128), rep("occupied_patient", 12)), each = 12)
  moves$acceptance <- as.numeric(alpha)
  moves$common_acceptance <- as.numeric(alpha_common)
  moves$scaled_squared_jump <- as.numeric(jump)
  moves$log_f_over_q <- as.numeric(result$log_w_reserve)
  moves$self_comparison <- moves$vector_id == paste0("patient_", moves$patient_id)
  combined <- rbind(paired$comparison, large$comparison)
  ledger <- rbind(paired$ledger, large$ledger)
  timing <- do.call(rbind, lapply(c("panel12", "grid115"), function(grid) {
    original <- subset(ledger, ledger$grid == grid & method == "original")
    common <- subset(ledger, ledger$grid == grid & method == "common")
    stopifnot(identical(original$vector_id, common$vector_id))
    # Integration count >0 and absence of failures distinguishes actual valid
    # integration cases from cheap equilibrium rejections.
    valid <- common$ode_integrations > 0 & common$ode_failures == 0 &
      original$ode_integrations > 0 & original$ode_failures == 0
    data.frame(snapshot = snapshot$id, grid = grid,
      original_seconds = sum(original$elapsed_seconds),
      common_seconds = sum(common$elapsed_seconds),
      elapsed_speedup = sum(original$elapsed_seconds)/sum(common$elapsed_seconds),
      valid_integration_speedup = sum(original$elapsed_seconds[valid])/
        sum(common$elapsed_seconds[valid]),
      valid_common_integrations = sum(valid),
      original_integrations = sum(original$ode_integrations),
      common_integrations = sum(common$ode_integrations))
  }))
  sharing <- mean(reserve_summary$useful_patients >= 2L)
  interval <- wilson(sum(reserve_summary$useful_patients >= 2L), 128L)
  error <- maximum_finite(combined$loglik_difference)
  alpha_error <- max(abs(alpha[, 1:128] - alpha_common[, 1:128]))
  peer <- moves[moves$source == "occupied_patient" & !moves$self_comparison, ]
  decision <- data.frame(snapshot = snapshot$id,
    direct_mean_acceptance = mean(direct$acceptance),
    fraction_shared_reserves = sharing, sharing_lower95 = interval[[1L]],
    sharing_upper95 = interval[[2L]],
    patients_mean_acceptance_ge_002 = sum(patient$mean_reserve_mh_acceptance >= .02),
    peer_mean_acceptance_excluding_self = mean(peer$acceptance),
    max_abs_loglik_difference = error, max_abs_reserve_acceptance_difference = alpha_error,
    support_disagreements = sum(!combined$finite_support_agrees),
    prediction_outcome_disagreements = sum(!combined$prediction_success_agrees),
    usefulness_pass = sharing >= .10 & sum(patient$mean_reserve_mh_acceptance >= .02) >= 6L,
    numerical_pass = all(combined$finite_support_agrees) &
      all(combined$prediction_success_agrees) & is.finite(error) & error <= .01 & alpha_error <= .01,
    speed_pass = all(is.finite(timing$elapsed_speedup) & timing$elapsed_speedup >= 3) &
      all(timing$valid_common_integrations >= 4L))
  decision$continue_to_discuss_reduced_sampler <- with(decision,
    usefulness_pass & numerical_pass & speed_pass)
  all_comparisons[[k]] <- combined; all_ledgers[[k]] <- ledger
  all_patients[[k]] <- patient; all_reserves[[k]] <- reserve_summary
  all_moves[[k]] <- moves; decisions[[k]] <- decision; timings[[k]] <- timing
  # Persist each completed context so scheduler/network interruptions lose no work.
  tables <- list(comparison = all_comparisons, ledger = all_ledgers, patient = all_patients,
                 reserve = all_reserves, moves = all_moves, decision = decisions, timing = timings)
  for (name in names(tables)) write.csv(do.call(rbind, tables[[name]]),
    file.path(output, paste0(name, ".csv")), row.names = FALSE)
  print(decision)
}
ledger <- do.call(rbind, all_ledgers)
stopifnot(sum(ledger$prediction_calls) <= 5496L)
decisions <- do.call(rbind, decisions)
report <- c(
  "Shared-trajectory feasibility screen; not a retained posterior sampler.",
  paste("Both contexts pass:", all(decisions$continue_to_discuss_reduced_sampler)),
  paste("Prediction-interface calls:", sum(ledger$prediction_calls)),
  paste("Actual integration attempts:", sum(ledger$ode_integrations)),
  paste("Integration failures:", sum(ledger$ode_failures)),
  paste("Elapsed wall seconds (including loading):", proc.time()[["elapsed"]] - started),
  "All direct-q draws including failed solves remain counted. No SAEM rerun.",
  "Common-grid target is experimental; this screen is not a full numerical certification.",
  "No global ESS, dynamic-psi cache lifetime, model-reuse or posterior agreement measured.",
  "A failed gate stops this design, not every conceivable sampler or reserve density.")
writeLines(report, file.path(output, "REPORT.txt"))
writeLines(capture.output(sessionInfo()), file.path(output, "session.txt"))
paths <- list.files(output, full.names = TRUE)
write.csv(data.frame(file = basename(paths), sha256 = vapply(paths, .sab_system_a_sha256_file, "")),
          file.path(output, "manifest.csv"), row.names = FALSE)
writeLines("completed", file.path(output, "COMPLETED"))
cat(paste(report, collapse = "\n"), "\n")
