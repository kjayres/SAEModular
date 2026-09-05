#!/usr/bin/env Rscript
# Actual reduced posterior comparison. Scientific work runs only under Slurm.
options(warn = 1)
if (!nzchar(Sys.getenv("SLURM_JOB_ID"))) stop("Run this experiment via Slurm.")
arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 1L) stop("Usage: system_a_shared_pool_chain.R <new-output-dir>")
started <- proc.time()
project <- normalizePath(Sys.getenv("SAB_PROJECT_ROOT"), mustWork = TRUE)
workspace <- normalizePath(file.path(project, "..", ".."), mustWork = TRUE)
output <- arguments[[1L]]
if (file.exists(output) || !dir.create(output, recursive = TRUE)) {
  stop("Refusing existing or uncreatable output directory: ", output)
}
modules <- file.path(project, "R", c("system_a_adapter.R", "shared_pool_snapshots.R",
  "shared_trajectory_probe.R", "shared_pool_kernel.R", "system_a_population_updates.R",
  "system_a_shared_chain.R"))
invisible(lapply(modules, sys.source, envir = globalenv()))
if (!requireNamespace("posterior", quietly = TRUE)) stop("posterior is required.")
smoke <- identical(Sys.getenv("SMOKE", "0"), "1")
config <- list(budget = 100000L, warmup = 500L, dynamic_every = 5L,
  reserve_refresh = 4L, pool_size = 36L, direct_probability = .2,
  initial_pcn_beta = sqrt(.19), initial_dynamic_sd = c(.04, .04),
  initial_noise_sd = .06, log_floor = -1e6, checkpoint_every = 250L,
  max_sweeps = 20000L)
if (smoke) config[c("budget", "warmup", "max_sweeps")] <- list(1500L, 10L, 60L)
ids <- c("3", "6", "16", "20", "55", "71", "74", "88", "105", "111", "117", "122")
full <- sab_load_system_a_adapter(workspace)
adapter <- sab_load_system_a_adapter(workspace, ids)
snapshots <- sab_shared_pool_snapshots(workspace, full, chains = 1:4)
components <- sab_shared_pool_reference_components(full, snapshots[c(1L, 3L)])
panel_probe <- sab_make_shared_trajectory_probe(adapter)
ode_ids <- names(panel_probe$positive_times)[lengths(panel_probe$positive_times) > 0L]
stopifnot(length(ode_ids) == 9L)
probe <- sab_make_shared_trajectory_probe(sab_load_system_a_adapter(workspace, ode_ids))
prior_mean <- adapter$prior_reference$eta
prior_sd <- adapter$prior_reference$eta_sd
stopifnot(identical(names(prior_mean), adapter$coordinate_names$population),
          identical(names(prior_sd), names(prior_mean)), all(prior_sd > 0))
cores <- if (smoke) 1L else min(4L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
stopifnot(is.finite(cores), cores >= 1L)
jobs <- expand.grid(method = c("baseline", "recycling"),
                    chain = seq_len(if (smoke) 1L else 4L), stringsAsFactors = FALSE)
jobs$seed <- 202609050L + 100L * jobs$chain + match(jobs$method, c("baseline", "recycling"))
source_paths <- c(modules, file.path(project, "experiments", "system_a_shared_pool_chain.R"),
                  file.path(project, "SHARED_POOL_CHAIN_CONTRACT.md"))
write.csv(data.frame(path = source_paths, sha256 = vapply(
  source_paths, .sab_system_a_sha256_file, "")), file.path(output, "sources.csv"), row.names = FALSE)
saveRDS(list(config = config, jobs = jobs, initial_snapshots = snapshots,
             reference_components = components, patient_ids = ids,
             exact_target = adapter$target_fingerprint,
             proposal_target = probe$target_fingerprint, smoke = smoke),
        file.path(output, "design.rds"))
preparation <- proc.time() - started
preparation_cpu <- unname(sum(preparation[c("user.self", "sys.self")]))
results <- parallel::mclapply(seq_len(nrow(jobs)), function(j) {
  job <- jobs[j, ]
  label <- sprintf("%s_chain%02d", job$method, job$chain)
  initial <- snapshots[[job$chain]]
  initial$x <- initial$x[ids, , drop = FALSE]
  tryCatch({
    value <- sab_run_system_a_shared_chain(adapter, probe, initial, components,
      prior_mean, prior_sd, method = job$method, seed = job$seed, config = config,
      checkpoint_path = file.path(output, paste0(label, "_checkpoint.rds")))
    value$chain <- job$chain
    saveRDS(value, file.path(output, paste0(label, ".rds")))
    value
  }, error = function(e) {
    writeLines(conditionMessage(e), file.path(output, paste0(label, "_FAILED.txt")))
    list(error = conditionMessage(e), method = job$method, chain = job$chain)
  })
}, mc.cores = cores, mc.preschedule = FALSE, mc.set.seed = FALSE)
failed <- vapply(results, function(x) !is.list(x) || !is.null(x$error) ||
                   !is.matrix(x$draws) || !is.data.frame(x$ledger), logical(1L))
if (any(failed)) stop("Chain worker failed; completed chains and failure messages were saved.")
ledger <- do.call(rbind, lapply(results, function(x) {
  value <- x$ledger; value$method <- x$method; value$chain <- x$chain; value
}))
write.csv(ledger, file.path(output, "ledger.csv"), row.names = FALSE)
counters <- do.call(rbind, lapply(results, function(x) data.frame(
  method = x$method, chain = x$chain, counter = names(x$counters), value = unname(x$counters))))
write.csv(counters, file.path(output, "counters.csv"), row.names = FALSE)
cost <- do.call(rbind, lapply(results, function(x) data.frame(method = x$method,
  chain = x$chain, saved_sweeps = nrow(x$draws), warmup = x$warmup,
  retained_sweeps = max(0L, nrow(x$draws) - x$warmup), cpu_sec = x$cpu_sec,
  elapsed_sec = x$elapsed_sec, budget_units = sum(x$ledger$budget_units))))
write.csv(cost, file.path(output, "cost.csv"), row.names = FALSE)
stopifnot(all(cost$budget_units <= config$budget), all(is.finite(cost$cpu_sec)),
          all(cost$cpu_sec > 0), all(cost$saved_sweeps <= config$max_sweeps))
summaries <- lengths <- list()
for (method in c("baseline", "recycling")) {
  chains <- Filter(function(x) identical(x$method, method), results)
  retained <- lapply(chains, function(x) x$draws[seq_len(nrow(x$draws)) > x$warmup, , drop = FALSE])
  n <- min(vapply(retained, nrow, integer(1L)))
  variables <- colnames(retained[[1L]])
  stopifnot(all(vapply(retained, function(x) identical(colnames(x), variables), logical(1L))))
  cpu <- sum(vapply(chains, `[[`, numeric(1L), "cpu_sec")) + preparation_cpu / 2
  elapsed <- sum(vapply(chains, `[[`, numeric(1L), "elapsed_sec")) + unname(preparation[["elapsed"]]) / 2
  units <- sum(vapply(chains, function(x) sum(x$ledger$budget_units), numeric(1L)))
  integrations <- sum(vapply(chains, function(x) sum(x$ledger$ode_integrations), numeric(1L)))
  lengths[[method]] <- data.frame(method = method, chains = length(chains),
    retained_prefix = n, unused_retained_sweeps = sum(vapply(retained, nrow, integer(1L))) - n * length(chains),
    charged_cpu_sec = cpu, charged_sum_chain_wall_sec = elapsed,
    charged_budget_units = units, charged_ode_integrations = integrations)
  summaries[[method]] <- do.call(rbind, lapply(variables, function(variable) {
    matrix <- vapply(retained, function(x) x[seq_len(n), variable], numeric(n))
    if (n >= 20L && all(is.finite(matrix))) {
      uncertainty <- c(posterior::mcse_mean(matrix), posterior::mcse_sd(matrix),
                       posterior::mcse_quantile(matrix, probs = c(.05, .5, .95)))
      diagnostic <- c(posterior::rhat(matrix), posterior::ess_bulk(matrix), posterior::ess_tail(matrix))
    } else {
      uncertainty <- rep(NA_real_, 5L); diagnostic <- rep(NA_real_, 3L)
    }
    estimates <- if (n > 0L) c(mean(matrix), stats::sd(as.numeric(matrix)),
                         stats::quantile(matrix, c(.05, .5, .95), names = FALSE)) else rep(NA_real_, 5L)
    data.frame(method = method, variable = variable, mean = estimates[1L], sd = estimates[2L],
      q05 = estimates[3L], q50 = estimates[4L], q95 = estimates[5L],
      mcse_mean = uncertainty[1L], mcse_sd = uncertainty[2L], mcse_q05 = uncertainty[3L],
      mcse_q50 = uncertainty[4L], mcse_q95 = uncertainty[5L], rhat = diagnostic[1L],
      bulk_ess = diagnostic[2L], tail_ess = diagnostic[3L],
      bulk_ess_per_cpu_sec = diagnostic[2L] / cpu, tail_ess_per_cpu_sec = diagnostic[3L] / cpu,
      bulk_ess_per_sum_wall_sec = diagnostic[2L] / elapsed,
      bulk_ess_per_1000_budget_units = 1000 * diagnostic[2L] / units,
      tail_ess_per_1000_budget_units = 1000 * diagnostic[3L] / units,
      bulk_ess_per_1000_ode_integrations = if (integrations > 0) 1000 * diagnostic[2L] / integrations else NA_real_,
      tail_ess_per_1000_ode_integrations = if (integrations > 0) 1000 * diagnostic[3L] / integrations else NA_real_,
      diagnostic_passed = !smoke && length(chains) == 4L && all(is.finite(diagnostic)) &&
        diagnostic[1L] <= 1.01 && diagnostic[2L] >= 400 && diagnostic[3L] >= 200)
  }))
}
summary <- do.call(rbind, summaries)
write.csv(summary, file.path(output, "posterior_summary.csv"), row.names = FALSE)
write.csv(do.call(rbind, lengths), file.path(output, "analysis_prefixes.csv"), row.names = FALSE)
baseline <- summaries$baseline; recycling <- summaries$recycling
stopifnot(identical(baseline$variable, recycling$variable))
statistics <- c("mean", "sd", "q05", "q50", "q95")
z <- stats::qnorm(1 - .05 / (2 * nrow(baseline) * length(statistics)))
agreement <- do.call(rbind, lapply(statistics, function(statistic) {
  mcse <- paste0("mcse_", statistic)
  error <- sqrt(baseline[[mcse]]^2 + recycling[[mcse]]^2)
  difference <- recycling[[statistic]] - baseline[[statistic]]
  data.frame(variable = baseline$variable, statistic = statistic,
    baseline = baseline[[statistic]], recycling = recycling[[statistic]],
    difference = difference, combined_mcse = error, simultaneous_z = z,
    tolerance = z * error, passed = is.finite(error) & is.finite(difference) & abs(difference) <= z * error)
}))
write.csv(agreement, file.path(output, "posterior_agreement.csv"), row.names = FALSE)
efficiency <- data.frame(variable = baseline$variable,
  bulk_cpu_ratio = recycling$bulk_ess_per_cpu_sec / baseline$bulk_ess_per_cpu_sec,
  tail_cpu_ratio = recycling$tail_ess_per_cpu_sec / baseline$tail_ess_per_cpu_sec,
  bulk_wall_ratio = recycling$bulk_ess_per_sum_wall_sec / baseline$bulk_ess_per_sum_wall_sec,
  bulk_budget_ratio = recycling$bulk_ess_per_1000_budget_units / baseline$bulk_ess_per_1000_budget_units,
  bulk_ode_ratio = recycling$bulk_ess_per_1000_ode_integrations / baseline$bulk_ess_per_1000_ode_integrations)
write.csv(efficiency, file.path(output, "efficiency.csv"), row.names = FALSE)
primary <- c("beta_nelf", "treated_logit_location", "log_omega_eta_pi", "mu_log_lambda", "log_omega_lambda")
stopifnot(all(c(primary, adapter$coordinate_names$dynamic_global) %in% efficiency$variable))
population_ratios <- efficiency$bulk_cpu_ratio[match(primary, efficiency$variable)]
dynamic_ratios <- efficiency$bulk_cpu_ratio[match(adapter$coordinate_names$dynamic_global, efficiency$variable)]
decision <- data.frame(smoke = smoke, all_global_diagnostics_pass = all(summary$diagnostic_passed),
  posterior_agreement_pass = all(agreement$passed), primary_median_bulk_cpu_ratio = median(population_ratios),
  primary_minimum_bulk_cpu_ratio = min(population_ratios), dynamic_minimum_bulk_cpu_ratio = min(dynamic_ratios))
decision$efficiency_pass <- all(is.finite(c(population_ratios, dynamic_ratios))) &&
  median(population_ratios) >= 1.5 && min(population_ratios) >= 1 && min(dynamic_ratios) >= .8
decision$continue_to_discuss <- !smoke && decision$all_global_diagnostics_pass &&
  decision$posterior_agreement_pass && decision$efficiency_pass
write.csv(decision, file.path(output, "decision.csv"), row.names = FALSE)
report <- c("Reduced 12-patient exact hierarchical comparison; original priors and sealed likelihood.",
  paste("Implementation smoke:", smoke), paste("All global convergence gates:", decision$all_global_diagnostics_pass),
  paste("Independent-arm posterior agreement:", decision$posterior_agreement_pass),
  paste("Efficiency gate:", decision$efficiency_pass), paste("Continue to discuss:", decision$continue_to_discuss),
  paste("Preparation CPU seconds:", preparation_cpu),
  paste("Driver elapsed seconds:", proc.time()[["elapsed"]] - started[["elapsed"]]),
  "All costs include warmup and unused retained tails; historical comparator fitting is excluded.",
  "Nonconvergence means an inconclusive posterior comparison, not a general impossibility result.",
  "No 115-patient inference, deployable SAEM cost or model-structure reuse is established.")
writeLines(report, file.path(output, "REPORT.txt"))
writeLines(capture.output(sessionInfo()), file.path(output, "session.txt"))
paths <- list.files(output, full.names = TRUE)
write.csv(data.frame(file = basename(paths), sha256 = vapply(paths, .sab_system_a_sha256_file, "")),
          file.path(output, "manifest.csv"), row.names = FALSE)
writeLines("completed", file.path(output, "COMPLETED"))
cat(paste(report, collapse = "\n"), "\n")
