#!/usr/bin/env Rscript
# Small actual conditional-chain comparison, not a joint or modular posterior.
options(warn = 1)
if (!nzchar(Sys.getenv("SLURM_JOB_ID"))) stop("Run through Slurm.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply a new output directory.")
started <- proc.time()
project <- normalizePath(Sys.getenv("SAB_PROJECT_ROOT"), mustWork = TRUE)
workspace <- normalizePath(file.path(project, "..", ".."), mustWork = TRUE)
output <- args[[1L]]
if (file.exists(output) || !dir.create(output, recursive = TRUE)) stop("Refusing output directory.")
modules <- file.path(project, "R", c("system_a_adapter.R", "shared_trajectory_probe.R", "system_a_coarse_guide.R"))
invisible(lapply(modules, sys.source, envir = globalenv()))
if (!requireNamespace("posterior", quietly = TRUE)) stop("posterior is required.")
paths <- file.path(project, c("outputs/residual_ceiling_122250/design_and_states.rds",
  "outputs/residual_ceiling_122250/ledger.csv", "outputs/system_a_12_oracle_rescue_v2/frozen_map.rds",
  "outputs/coarse_guide_122252/decision.csv", "outputs/coarse_guide_122252/design.rds"))
hashes <- c("f350728c91de7b568ff50e44ece79408d9c7682c534bde473b5716d8c1125c4c",
  "80882179f21135c3d15ee825f46c4e368e63eadb35dd26e503bf60cd1a60a8e1",
  "ec4b29a4efc7162da43f2295c0530d1f7a02e72e1c785df1cd0a7e3b6b590a78",
  "8a9979ca489137c764b7b6a588fa12b82e096eb00bc1b493d6bcbcd0b52b1f60",
  "6553bec3c7be7b0cf3307a4bd7452b4c7c0aa3d6991db0faa945410477865003")
stopifnot(identical(unname(vapply(paths, .sab_system_a_sha256_file, "")), hashes))
saved <- readRDS(paths[[1L]]); source_ledger <- read.csv(paths[[2L]])
map <- readRDS(paths[[3L]]); preflight <- read.csv(paths[[4L]]); design <- readRDS(paths[[5L]])
eligible <- preflight[preflight$chain_pilot_go, ]
if (!nrow(eligible)) stop("No guide passed the preflight.")
profile <- eligible$profile[[which.max(eligible$median_cpu_speedup)]]
tol <- design$profiles[[profile]]
stopifnot(identical(profile, "rtol_1e1"))
smoke <- identical(Sys.getenv("SMOKE", "0"), "1")
ids <- if (smoke) "3" else c("3", "74", "122")
contexts <- if (smoke) design$contexts["centre"] else design$contexts
adapter <- sab_load_system_a_adapter(workspace, ids)
stopifnot(identical(saved$target_fingerprint, adapter$target_fingerprint))
guide <- sab_make_system_a_coarse_guide(adapter, tol[[1L]], tol[[2L]])
local_names <- adapter$coordinate_names$local
config <- list(iterations = c(baseline = 6000L, macro = 1500L),
  warmup = c(baseline = 500L, macro = 125L), inner_steps = 8L,
  preflight_cpu_sec = 40.182, preflight_wall_sec = 46, checkpoint_every = 500L)
if (smoke) { config$iterations <- c(baseline = 24L, macro = 8L); config$warmup <- c(baseline = 4L, macro = 2L) }
jobs <- expand.grid(method = c("baseline", "macro"), chain = seq_len(if (smoke) 1L else 4L),
                    patient_id = ids, context = names(contexts), stringsAsFactors = FALSE)
jobs$seed <- 2026090600L + seq_len(nrow(jobs))
sources <- c(paths, modules, file.path(project, c("experiments/system_a_coarse_chain.R", "COARSE_GUIDE_CHAIN_CONTRACT.md")))
write.csv(data.frame(path = sources, sha256 = vapply(sources, .sab_system_a_sha256_file, "")),
  file.path(output, "sources.csv"), row.names = FALSE)
saveRDS(list(jobs = jobs, contexts = contexts, eta = saved$eta, config = config,
  profile = profile, target_fingerprint = adapter$target_fingerprint, smoke = smoke), file.path(output, "design.rds"))

run <- function(job) {
  set.seed(job$seed)
  begin_chain <- proc.time()
  id <- job$patient_id; psi <- contexts[[job$context]]; eta <- saved$eta
  label <- paste(id, job$context, job$method, job$chain, sep = "_")
  n <- config$iterations[[job$method]]; warmup <- config$warmup[[job$method]]
  budget <- matrix(0, 3L, 4L, dimnames = list(c("initial", "exact", "guide"),
    c("attempts", "ode_integrations", "failures", "cpu_sec")))
  exact_calls <- cheap_calls <- 0L
  evaluate_exact <- function(x, category = "exact") {
    exact_calls <<- exact_calls + 1L
    if (exact_calls > n + 16L) stop("Chain exact-call cap exceeded.")
    t <- proc.time()
    prediction <- adapter$solve_prediction(id, x, psi)
    ll <- adapter$loglik_from_prediction(prediction, psi)
    value <- list(log_target = ll + adapter$log_population_density(id, x, eta),
                  loglik = ll, prediction = prediction)
    dt <- proc.time() - t
    budget[category, ] <<- budget[category, ] + c(1, .sab_shared_original_calls(prediction),
      !isTRUE(prediction$ok), sum(dt[c("user.self", "sys.self")]))
    value
  }
  evaluate_cheap <- function(x) {
    cheap_calls <<- cheap_calls + 1L
    if (cheap_calls > n * config$inner_steps + 1L) stop("Chain cheap-call cap exceeded.")
    t <- proc.time(); value <- guide$evaluate(id, x, psi); dt <- proc.time() - t
    budget["guide", ] <<- budget["guide", ] + c(1, value$ode_integrations, !value$ok,
      sum(dt[c("user.self", "sys.self")]))
    list(log_target = value$loglik + adapter$log_population_density(id, x, eta), loglik = value$loglik)
  }
  rows <- which(as.character(source_ledger$patient_id) == id & source_ledger$phase == "validation")
  start <- c(193L, 257L, 321L, 385L)[[job$chain]]
  current <- NULL
  for (index in start + 0:15) {
    x <- setNames(as.numeric(saved$states[rows[[index]], paste0("original_", local_names)]), local_names)
    exact <- evaluate_exact(x, "initial")
    if (is.finite(exact$log_target)) { current <- list(x = x, exact = exact); break }
  }
  if (is.null(current)) stop("No finite predeclared initial state.")
  if (job$method == "macro") current$cheap <- evaluate_cheap(current$x)
  coef <- map$coefficients[[id]]; index <- coef$lower_index
  entries <- coef$chol[, "intercept"]; diag_entry <- index[, 1L] == index[, 2L]
  entries[diag_entry] <- exp(entries[diag_entry])
  scale <- matrix(0, 8L, 8L); scale[index] <- entries
  log_step <- log(2.38 / sqrt(8))
  propose <- function(x) list(x = setNames(x + exp(log_step) * as.numeric(scale %*% stats::rnorm(8L)), local_names),
                              log_reverse_minus_forward = 0)
  draws <- matrix(NA_real_, n, 9L, dimnames = list(NULL, c(local_names, "log_likelihood")))
  attempts <- accepts <- inner_attempts <- inner_accepts <- moves <- 0L
  snapshot <- function(iteration) {
    dt <- proc.time() - begin_chain
    list(job = job, draws = draws[seq_len(iteration), , drop = FALSE], warmup = warmup,
      ledger = data.frame(category = rownames(budget), budget, row.names = NULL),
      tuning = exp(log_step), attempts = attempts, accepts = accepts, moves = moves,
      inner_attempts = inner_attempts, inner_accepts = inner_accepts,
      exact_calls = exact_calls, cheap_calls = cheap_calls,
      cpu_sec = unname(sum(dt[c("user.self", "sys.self")])), elapsed_sec = unname(dt[["elapsed"]]))
  }
  for (iteration in seq_len(n)) {
    if (job$method == "baseline") {
      candidate <- propose(current$x); value <- evaluate_exact(candidate$x)
      accepted <- log(stats::runif(1L)) < min(0, value$log_target - current$exact$log_target)
      moved <- accepted; inner_rate <- as.numeric(accepted)
      if (accepted) current <- list(x = candidate$x, exact = value)
      inner_attempts <- inner_attempts + 1L; inner_accepts <- inner_accepts + accepted
    } else {
      value <- sab_surrogate_macro_step(current, propose, evaluate_cheap, evaluate_exact, config$inner_steps)
      current <- value$current; accepted <- value$accepted; moved <- value$moved
      inner_rate <- value$inner_accepts / config$inner_steps
      inner_attempts <- inner_attempts + config$inner_steps; inner_accepts <- inner_accepts + value$inner_accepts
    }
    attempts <- attempts + 1L; accepts <- accepts + accepted; moves <- moves + moved
    if (iteration <= warmup) log_step <- max(log(.05), min(log(3),
      log_step + (iteration + 10)^(-.6) * (inner_rate - .234)))
    draws[iteration, ] <- c(current$x, current$exact$loglik)
    if (iteration %% config$checkpoint_every == 0L) saveRDS(snapshot(iteration),
      file.path(output, paste0(label, "_checkpoint.rds")))
  }
  result <- snapshot(n)
  saveRDS(result, file.path(output, paste0(label, ".rds")))
  cat("Finished", label, "CPU", result$cpu_sec, "exact", exact_calls, "guide", cheap_calls, "\n")
  flush.console()
  result
}
preparation <- proc.time() - started
cores <- if (smoke) 1L else min(4L, as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", "1")))
results <- parallel::mclapply(seq_len(nrow(jobs)), function(j) {
  tryCatch(run(jobs[j, ]), error = function(e) {
    writeLines(conditionMessage(e), file.path(output, paste0("job_", j, "_FAILED.txt")))
    list(error = conditionMessage(e))
  })
}, mc.cores = cores, mc.preschedule = FALSE, mc.set.seed = FALSE)
if (any(vapply(results, function(x) !is.list(x) || !is.null(x$error), logical(1L)))) stop("A chain failed; see saved failure records.")
cost <- do.call(rbind, lapply(results, function(x) cbind(x$job, data.frame(cpu_sec = x$cpu_sec,
  elapsed_sec = x$elapsed_sec, exact_calls = x$exact_calls, guide_calls = x$cheap_calls,
  accepted = x$accepts, moves = x$moves, outer_attempts = x$attempts,
  inner_attempts = x$inner_attempts, inner_accepts = x$inner_accepts, proposal_scale = x$tuning))))
write.csv(cost, file.path(output, "cost.csv"), row.names = FALSE)
ledger <- do.call(rbind, lapply(results, function(x) cbind(x$job[rep(1, nrow(x$ledger)), ], x$ledger)))
write.csv(ledger, file.path(output, "ledger.csv"), row.names = FALSE)
summary <- list()
for (id in ids) for (context in names(contexts)) for (method in c("baseline", "macro")) {
  chains <- Filter(function(x) x$job$patient_id == id && x$job$context == context && x$job$method == method, results)
  retained <- lapply(chains, function(x) {
    draw <- x$draws[seq_len(nrow(x$draws)) > x$warmup, , drop = FALSE]
    z <- sweep(sweep(draw[, local_names, drop = FALSE], 2L, adapter$population_mean(id, saved$eta), "-"),
      2L, exp(saved$eta[10:17]), "/")
    squares <- z^2 - 1; colnames(squares) <- paste0("variance_score_", seq_len(8L))
    cbind(draw, squares)
  })
  n <- nrow(retained[[1L]])
  shared_cost <- sum(preparation[c("user.self", "sys.self")]) / (2 * length(ids) * length(contexts))
  cpu <- sum(vapply(chains, `[[`, numeric(1L), "cpu_sec")) + shared_cost +
    if (method == "macro") config$preflight_cpu_sec / (length(ids) * length(contexts)) else 0
  wall <- sum(vapply(chains, `[[`, numeric(1L), "elapsed_sec")) +
    preparation[["elapsed"]] / (2 * length(ids) * length(contexts)) +
    if (method == "macro") config$preflight_wall_sec / (length(ids) * length(contexts)) else 0
  exact_count <- sum(vapply(chains, `[[`, numeric(1L), "exact_calls"))
  for (variable in colnames(retained[[1L]])) {
    m <- vapply(retained, function(x) x[, variable], numeric(n))
    valid <- !smoke && ncol(m) == 4L && all(is.finite(m))
    diag <- if (valid) c(posterior::rhat(m), posterior::ess_bulk(m), posterior::ess_tail(m)) else rep(NA_real_, 3)
    mcse <- if (valid) c(posterior::mcse_mean(m), posterior::mcse_sd(m),
      posterior::mcse_quantile(m, probs = c(.05, .5, .95))) else rep(NA_real_, 5)
    mean_ess <- if (valid) posterior::ess_mean(m) else NA_real_
    estimates <- c(mean(m), stats::sd(as.numeric(m)), stats::quantile(m, c(.05, .5, .95), names = FALSE))
    lag_one <- mean(apply(m, 2L, function(v) if (length(v) > 2 && stats::sd(v) > 0) stats::cor(head(v, -1), tail(v, -1)) else NA_real_))
    summary[[length(summary) + 1L]] <- data.frame(patient_id = id, context = context, method = method, variable = variable,
      mean = estimates[[1L]], sd = estimates[[2L]], q05 = estimates[[3L]], q50 = estimates[[4L]], q95 = estimates[[5L]],
      mcse_mean = mcse[[1L]], mcse_sd = mcse[[2L]], mcse_q05 = mcse[[3L]], mcse_q50 = mcse[[4L]], mcse_q95 = mcse[[5L]],
      rhat = diag[[1L]], bulk_ess = diag[[2L]], tail_ess = diag[[3L]],
      bulk_iact = length(m) / diag[[2L]], mean_iact = length(m) / mean_ess,
      lag_one = lag_one, charged_cpu_sec = cpu, charged_sum_wall_sec = wall,
      bulk_ess_per_cpu = diag[[2L]] / cpu, tail_ess_per_cpu = diag[[3L]] / cpu,
      bulk_ess_per_sum_wall_sec = diag[[2L]] / wall,
      bulk_ess_per_1000_exact_calls = 1000 * diag[[2L]] / exact_count,
      convergence_pass = valid && all(is.finite(diag)) && diag[[1L]] <= 1.01 && diag[[2L]] >= 400 && diag[[3L]] >= 200)
  }
}
summary <- do.call(rbind, summary)
write.csv(summary, file.path(output, "posterior_summary.csv"), row.names = FALSE)
baseline <- summary[summary$method == "baseline", ]; macro <- summary[summary$method == "macro", ]
key <- function(x) paste(x$patient_id, x$context, x$variable)
stopifnot(identical(key(baseline), key(macro)))
efficiency <- baseline[, c("patient_id", "context", "variable")]
efficiency$bulk_cpu_ratio <- macro$bulk_ess_per_cpu / baseline$bulk_ess_per_cpu
efficiency$tail_cpu_ratio <- macro$tail_ess_per_cpu / baseline$tail_ess_per_cpu
efficiency$bulk_wall_ratio <- macro$bulk_ess_per_sum_wall_sec / baseline$bulk_ess_per_sum_wall_sec
efficiency$bulk_exact_call_ratio <- macro$bulk_ess_per_1000_exact_calls / baseline$bulk_ess_per_1000_exact_calls
write.csv(efficiency, file.path(output, "efficiency.csv"), row.names = FALSE)
statistics <- c("mean", "sd", "q05", "q50", "q95")
z <- stats::qnorm(1 - .05 / (2 * nrow(baseline) * length(statistics)))
agreement <- do.call(rbind, lapply(statistics, function(statistic) {
  error <- sqrt(baseline[[paste0("mcse_", statistic)]]^2 + macro[[paste0("mcse_", statistic)]]^2)
  difference <- macro[[statistic]] - baseline[[statistic]]
  cbind(baseline[, c("patient_id", "context", "variable")], data.frame(statistic = statistic,
    difference = difference, combined_mcse = error, simultaneous_z = z,
    pass = is.finite(error) & is.finite(difference) & abs(difference) <= z * error))
}))
write.csv(agreement, file.path(output, "posterior_agreement.csv"), row.names = FALSE)
score <- efficiency[efficiency$variable != "log_likelihood", ]
group_medians <- vapply(split(score$bulk_cpu_ratio, interaction(score$patient_id, score$context)), median, 0)
decision <- data.frame(smoke = smoke, all_convergence_pass = all(summary$convergence_pass),
  all_agreement_pass = all(agreement$pass), median_score_bulk_cpu_ratio = median(score$bulk_cpu_ratio),
  minimum_case_median_ratio = min(group_medians), minimum_score_ratio = min(score$bulk_cpu_ratio))
decision$efficiency_pass <- !smoke && all(is.finite(score$bulk_cpu_ratio)) &&
  decision$median_score_bulk_cpu_ratio >= 1.5 && decision$minimum_case_median_ratio >= 1 && decision$minimum_score_ratio >= .8
decision$go_to_modular_reuse_test <- !smoke && decision$all_convergence_pass && decision$all_agreement_pass && decision$efficiency_pass
write.csv(decision, file.path(output, "decision.csv"), row.names = FALSE)
report <- c("Exact conditional patient sampling, not a joint population or modular posterior.",
  paste("Implementation smoke:", smoke), paste("All convergence gates:", decision$all_convergence_pass),
  paste("All posterior agreement gates:", decision$all_agreement_pass),
  paste("Median score ESS/CPU ratio:", decision$median_score_bulk_cpu_ratio),
  paste("Go to modular reuse test:", decision$go_to_modular_reuse_test),
  paste("Exact calls:", sum(cost$exact_calls)), paste("Guide calls:", sum(cost$guide_calls)),
  paste("Elapsed driver seconds:", (proc.time() - started)[["elapsed"]]),
  "All chain, warmup and initialization CPU is charged; full preflight CPU is charged to the guide.",
  "Historical oracle proposal covariance fitting is excluded from both arms; no deployable cost claim.",
  "Nonconvergence makes posterior comparison inconclusive, not proof of bias.",
  "Modular reuse across population structures and full shared uncertainty remain untested.")
writeLines(report, file.path(output, "REPORT.txt"))
writeLines(capture.output(sessionInfo()), file.path(output, "session.txt"))
paths <- list.files(output, full.names = TRUE)
write.csv(data.frame(file = basename(paths), sha256 = vapply(paths, .sab_system_a_sha256_file, "")),
  file.path(output, "manifest.csv"), row.names = FALSE)
writeLines("completed", file.path(output, "COMPLETED"))
cat(paste(report, collapse = "\n"), "\n")
