#!/usr/bin/env Rscript
# Cheap numerical guide feasibility only; a pass authorises a frozen chain pilot.
options(warn = 1)
if (!nzchar(Sys.getenv("SLURM_JOB_ID"))) stop("Run through Slurm.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply a new output directory.")
project <- normalizePath(Sys.getenv("SAB_PROJECT_ROOT"), mustWork = TRUE)
workspace <- normalizePath(file.path(project, "..", ".."), mustWork = TRUE)
output <- args[[1L]]
if (file.exists(output) || !dir.create(output, recursive = TRUE)) stop("Refusing output directory.")
modules <- file.path(project, "R", c("system_a_adapter.R", "shared_trajectory_probe.R", "system_a_coarse_guide.R"))
invisible(lapply(modules, sys.source, envir = globalenv()))
source_paths <- file.path(project, c("outputs/residual_ceiling_122250/design_and_states.rds",
  "outputs/residual_ceiling_122250/ledger.csv", "outputs/system_a_12_oracle_rescue_v2/frozen_map.rds"))
hashes <- c("f350728c91de7b568ff50e44ece79408d9c7682c534bde473b5716d8c1125c4c",
  "80882179f21135c3d15ee825f46c4e368e63eadb35dd26e503bf60cd1a60a8e1",
  "ec4b29a4efc7162da43f2295c0530d1f7a02e72e1c785df1cd0a7e3b6b590a78")
stopifnot(identical(unname(vapply(source_paths, .sab_system_a_sha256_file, "")), hashes))
saved <- readRDS(source_paths[[1L]]); old_ledger <- read.csv(source_paths[[2L]])
map <- readRDS(source_paths[[3L]])
ids <- c("3", "74", "122")
adapter <- sab_load_system_a_adapter(workspace, ids)
stopifnot(identical(saved$target_fingerprint, adapter$target_fingerprint),
          identical(map$target_fingerprint, adapter$target_fingerprint))
local_names <- adapter$coordinate_names$local
contexts <- list(centre = saved$psi, dynamic_1_plus = saved$psi, dynamic_2_plus = saved$psi)
contexts$dynamic_1_plus[1:2] <- contexts$dynamic_1_plus[1:2] + map$dynamic_chol[, 1L]
contexts$dynamic_2_plus[1:2] <- contexts$dynamic_2_plus[1:2] + map$dynamic_chol[, 2L]
profiles <- list(rtol_1e4 = c(1e-4, 1e-6), rtol_1e2 = c(1e-2, 1e-4), rtol_1e1 = c(1e-1, 1e-3))
guides <- lapply(profiles, function(tol) sab_make_system_a_coarse_guide(adapter, tol[[1L]], tol[[2L]]))
evaluate_exact <- function(id, x, psi) {
  p <- adapter$solve_prediction(id, x, psi)
  ll <- adapter$loglik_from_prediction(p, psi)
  list(loglik = ll, raw_loglik = ll,
       ok = isTRUE(p$ok), reason = if (isTRUE(p$ok)) "ok" else p$reason,
       ode_integrations = .sab_shared_original_calls(p))
}
ledger <- list(); n <- 0L
started <- proc.time()
for (id in ids) {
  selected <- head(which(as.character(old_ledger$patient_id) == id & old_ledger$phase == "validation"), 96L)
  stopifnot(length(selected) == 96L)
  states <- saved$states[selected, paste0("original_", local_names), drop = FALSE]
  colnames(states) <- local_names
  for (context_name in names(contexts)) {
    psi <- contexts[[context_name]]
    for (j in seq_len(nrow(states))) for (repeat_id in 1:2) {
      names_all <- c("exact", names(profiles))
      order <- ((seq_along(names_all) + j + repeat_id - 2L) %% length(names_all)) + 1L
      for (profile in names_all[order]) {
        begin <- proc.time()
        value <- if (profile == "exact") evaluate_exact(id, states[j, ], psi) else
          guides[[profile]]$evaluate(id, states[j, ], psi)
        cost <- proc.time() - begin
        n <- n + 1L
        if (n > 6912L) stop("Predeclared call cap exceeded.")
        ledger[[n]] <- data.frame(patient_id = id, context = context_name, state = j,
          repeat_id = repeat_id, profile = profile, loglik = value$loglik,
          raw_loglik = value$raw_loglik, ok = value$ok, reason = value$reason,
          attempts = 1L, ode_integrations = value$ode_integrations,
          cpu_sec = sum(cost[c("user.self", "sys.self")]), elapsed_sec = cost[["elapsed"]])
      }
    }
    cat("Finished", id, context_name, "calls", n, "\n"); flush.console()
  }
}
ledger <- do.call(rbind, ledger)
summaries <- list()
for (id in ids) for (context in names(contexts)) {
  group <- ledger[ledger$patient_id == id & ledger$context == context, ]
  exact <- group[group$profile == "exact" & group$repeat_id == 1L, ]
  exact <- exact[order(exact$state), ]
  exact_repeat <- group[group$profile == "exact" & group$repeat_id == 2L, ]
  exact_repeat <- exact_repeat[order(exact_repeat$state), ]
  stopifnot(identical(exact$loglik, exact_repeat$loglik))
  for (profile in names(profiles)) {
    approx <- group[group$profile == profile & group$repeat_id == 1L, ]
    approx <- approx[order(approx$state), ]
    finite <- is.finite(exact$loglik)
    residual <- exact$loglik - approx$loglik
    differences <- residual[seq(2, 96, 2)] - residual[seq(1, 95, 2)]
    differences <- differences[is.finite(differences)]
    repeated <- group[group$profile == profile & group$repeat_id == 2L, ]
    repeated <- repeated[order(repeated$state), ]
    stopifnot(identical(approx$loglik, repeated$loglik))
    exact_time <- sum(group$cpu_sec[group$profile == "exact"])
    approx_time <- sum(group$cpu_sec[group$profile == profile])
    speedup <- if (approx_time > 0) exact_time / approx_time else NA_real_
    error95 <- if (length(differences) >= 24L) unname(quantile(abs(differences), .95)) else Inf
    summaries[[length(summaries) + 1L]] <- data.frame(patient_id = id, context = context,
      profile = profile, exact_cpu_sec = exact_time, guide_cpu_sec = approx_time,
      cpu_speedup = speedup, exact_finite_fraction = mean(finite),
      false_negative_guide_failures = sum(finite & !approx$ok),
      absolute_residual_difference95 = error95, finite_design_pairs = length(differences),
      case_pass = is.finite(speedup) && speedup >= 1.5 && mean(finite) >= .9 &&
        !any(finite & !approx$ok) && error95 <= 1)
  }
}
summary <- do.call(rbind, summaries)
decision <- do.call(rbind, lapply(split(summary, summary$profile), function(g) {
  data.frame(profile = g$profile[[1L]], median_cpu_speedup = median(g$cpu_speedup),
    minimum_cpu_speedup = min(g$cpu_speedup), worst_residual_difference95 = max(g$absolute_residual_difference95),
    passing_cases = sum(g$case_pass), cases = nrow(g),
    chain_pilot_go = all(g$case_pass) && median(g$cpu_speedup) >= 2)
}))
write.csv(ledger, file.path(output, "ledger.csv"), row.names = FALSE)
write.csv(summary, file.path(output, "summary.csv"), row.names = FALSE)
write.csv(decision, file.path(output, "decision.csv"), row.names = FALSE)
saveRDS(list(contexts = contexts, profiles = profiles, selected_source_rows = lapply(ids, function(id)
  head(which(as.character(old_ledger$patient_id) == id & old_ledger$phase == "validation"), 96L)),
  target_fingerprint = adapter$target_fingerprint), file.path(output, "design.rds"))
source_paths <- c(source_paths, modules, file.path(project, c("experiments/system_a_coarse_guide.R", "COARSE_GUIDE_CONTRACT.md")))
write.csv(data.frame(path = source_paths, sha256 = vapply(source_paths, .sab_system_a_sha256_file, "")),
  file.path(output, "sources.csv"), row.names = FALSE)
report <- c("Coarse-tolerance guide preflight, not retained Bayesian sampling.",
  paste("Exact calls:", sum(ledger$attempts[ledger$profile == "exact"])),
  paste("Guide calls:", sum(ledger$attempts[ledger$profile != "exact"])),
  paste("Qualifying profiles:", sum(decision$chain_pilot_go)),
  paste("Elapsed seconds:", (proc.time() - started)[["elapsed"]]),
  "Residual differences are evaluated on fixed proposal-design states, not conditional posterior draws.",
  "No modularity, new-population-model reuse or full uncertainty claim follows from this preflight.")
writeLines(report, file.path(output, "REPORT.txt"))
writeLines(capture.output(sessionInfo()), file.path(output, "session.txt"))
paths <- list.files(output, full.names = TRUE)
write.csv(data.frame(file = basename(paths), sha256 = vapply(paths, .sab_system_a_sha256_file, "")),
  file.path(output, "manifest.csv"), row.names = FALSE)
writeLines("completed", file.path(output, "COMPLETED"))
cat(paste(report, collapse = "\n"), "\n")
