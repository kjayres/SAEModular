#!/usr/bin/env Rscript
# Structural checks of the supplied bound recipe, not a new sampler.
options(warn = 1)
if (!nzchar(Sys.getenv("SLURM_JOB_ID"))) stop("Run via Slurm, not the login node.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply a new output directory.")
project <- normalizePath(Sys.getenv("SAB_PROJECT_ROOT"), mustWork = TRUE)
workspace <- normalizePath(file.path(project, "..", ".."), mustWork = TRUE)
output <- args[[1L]]
if (file.exists(output) || !dir.create(output, recursive = TRUE)) stop("Refusing output directory.")
modules <- file.path(project, "R", c("system_a_adapter.R", "shared_trajectory_probe.R"))
invisible(lapply(modules, sys.source, envir = globalenv()))
source_path <- file.path(project, "outputs/shared_pool_probe_122228/design.rds")
source_hash <- "51fc4d28614e42563f9b0f2f460a77f9b630e3a96f5920bc17933d4189ce3fd8"
stopifnot(identical(.sab_system_a_sha256_file(source_path), source_hash))
design <- readRDS(source_path)
adapter <- sab_load_system_a_adapter(workspace, "3")
stopifnot(identical(design$reference_target, adapter$target_fingerprint))
owner <- environment(adapter$solve_prediction)
upstream <- get("upstream", owner, inherits = FALSE)
patient <- get("exact_likelihood", owner, inherits = FALSE)$patients[["3"]]
stopifnot(any(patient$obs_time > 0), any(patient$obs_time <= 0 & patient$ytype == 2L))
baseline_row <- which(patient$obs_time <= 0 & patient$ytype == 2L)[[1L]]
started <- proc.time()[["elapsed"]]
ledger <- support <- tails <- regions <- list()
counter <- 0L
evaluate <- function(x, psi, context, category, label) {
  counter <<- counter + 1L
  if (counter > 40L) stop("Predeclared prediction-call cap exceeded.")
  begin <- proc.time()[["elapsed"]]
  prediction <- adapter$solve_prediction("3", x, psi)
  likelihood <- adapter$loglik_from_prediction(prediction, psi)
  ledger[[counter]] <<- data.frame(context = context, category = category, label = label,
    prediction_calls = 1L, ode_integrations = .sab_shared_original_calls(prediction),
    ok = isTRUE(prediction$ok), reason = if (isTRUE(prediction$ok)) "ok" else prediction$reason,
    log_likelihood = likelihood, elapsed_sec = proc.time()[["elapsed"]] - begin)
  list(prediction = prediction, loglik = likelihood)
}
set.seed(2026090501)
for (s in design$snapshots) {
  x <- s$x["3", ]; eta <- s$eta; psi <- s$psi
  natural <- upstream$sysa_local_to_natural(x)
  shared <- upstream$sysa_psi_to_natural(psi)
  equilibrium <- upstream$sysa_stan_equilibrium(natural, shared, patient$controls$mu_v)
  stopifnot(!is.null(equilibrium))
  anchor <- evaluate(x, psi, s$id, "anchor", "current")
  stopifnot(is.finite(anchor$loglik))
  critical_lambda <- natural[["mu_t"]] * equilibrium[["T"]]
  coefficient <- (1 - natural[["pi"]]) / (natural[["alpha_l"]] + shared[["mu_l"]]) +
    (natural[["alpha_l"]] + natural[["pi"]] * shared[["mu_l"]]) /
      (natural[["mu_a"]] * (natural[["alpha_l"]] + shared[["mu_l"]]))
  anchor_cd4 <- equilibrium[["T"]] + coefficient * (natural[["lambda"]] - critical_lambda)
  stopifnot(isTRUE(all.equal(anchor_cd4, anchor$prediction$observation_mean[[baseline_row]],
                            tolerance = 1e-11)))
  for (factor in c(.5, .9)) {
    candidate <- x; candidate[[1L]] <- log(critical_lambda * factor)
    value <- evaluate(candidate, psi, s$id, "invalid_support", as.character(factor))
    log_g <- adapter$log_population_density("3", candidate, eta)
    stopifnot(identical(value$prediction$reason, "invalid_equilibrium"),
              identical(value$loglik, -Inf), is.finite(log_g))
    support[[length(support) + 1L]] <- data.frame(context = s$id, lambda_fraction = factor,
      log_lambda = candidate[[1L]], log_critical_lambda = log(critical_lambda),
      log_likelihood = value$loglik, log_population_density = log_g,
      conditional_invalid_probability_given_other_anchor_coordinates = stats::pnorm(
        log(critical_lambda), eta[[1L]], exp(eta[[10L]])))
  }
  for (d in c(.1, .25, .5, 1, 2, 4)) {
    candidate <- x; candidate[[1L]] <- x[[1L]] + d
    value <- evaluate(candidate, psi, s$id, "affine_tail", as.character(d))
    analytic <- equilibrium[["T"]] + coefficient * (exp(candidate[[1L]]) - critical_lambda)
    tangent <- anchor_cd4 + coefficient * natural[["lambda"]] * d
    residual <- coefficient * natural[["lambda"]] * (expm1(d) - d)
    observed <- if (isTRUE(value$prediction$ok)) value$prediction$observation_mean[[baseline_row]] else NA_real_
    stopifnot(isTRUE(all.equal(analytic - tangent, residual, tolerance = 1e-10)))
    if (is.finite(observed)) stopifnot(isTRUE(all.equal(analytic, observed, tolerance = 1e-10)))
    tails[[length(tails) + 1L]] <- data.frame(context = s$id, log_lambda_displacement = d,
      analytic_baseline_cd4 = analytic, sealed_baseline_cd4 = observed,
      anchor_tangent_cd4 = tangent, absolute_affine_error = residual,
      error_in_fixed_anchor_cd4_sd = residual / (exp(psi[[4L]]) * anchor_cd4),
      solver_ok = isTRUE(value$prediction$ok), log_likelihood = value$loglik)
  }
  # Cheap population draws check the support algebra, not posterior inference.
  n_draws <- 4096L
  mean <- adapter$population_mean("3", eta); sd <- exp(eta[10:17])
  bank <- sweep(sweep(matrix(stats::rnorm(n_draws * 8L), n_draws, 8L), 2L, sd, "*"), 2L, mean, "+")
  colnames(bank) <- adapter$coordinate_names$local
  valid <- apply(bank, 1L, function(v) !is.null(upstream$sysa_stan_equilibrium(
    upstream$sysa_local_to_natural(v), shared, patient$controls$mu_v)))
  a <- bank[, 1L] + bank[, 4L] - bank[, 2L] - bank[, 3L] + psi[[1L]] - log(patient$controls$mu_v)
  a_mean <- mean[[1L]] + mean[[4L]] - mean[[2L]] - mean[[3L]] + psi[[1L]] - log(patient$controls$mu_v)
  a_sd <- sqrt(sum(sd[1:4]^2))
  mass <- numeric(17L); members <- vector("list", 17L)
  for (j in seq_len(17L)) {
    k <- 2^(-9L + j)
    alpha_threshold <- psi[[2L]] + log(k); a_threshold <- log1p(1 / k)
    mass[[j]] <- stats::pnorm(alpha_threshold, mean[[5L]], sd[[5L]], lower.tail = FALSE) *
      stats::pnorm(a_threshold, a_mean, a_sd, lower.tail = FALSE)
    inside <- bank[, 5L] >= alpha_threshold & a >= a_threshold
    stopifnot(!any(inside & !valid))
    members[[j]] <- which(inside)
    regions[[length(regions) + 1L]] <- data.frame(context = s$id, k = k,
      gaussian_population_mass = mass[[j]], empirical_population_mass = base::mean(inside),
      membership_mcse = sqrt(mass[[j]] * (1 - mass[[j]]) / n_draws),
      population_draws = n_draws, physical_equilibrium_valid_fraction = base::mean(valid),
      sampled_inner_region_violations = sum(inside & !valid))
  }
  best <- which.max(mass)
  selected <- head(members[[best]], 4L)
  for (j in selected) evaluate(bank[j, ], psi, s$id, "support_region_solver",
                              paste0("k=", 2^(-9L + best), ";draw=", j))
}
write.csv(do.call(rbind, support), file.path(output, "support_witnesses.csv"), row.names = FALSE)
write.csv(do.call(rbind, tails), file.path(output, "affine_tail.csv"), row.names = FALSE)
write.csv(do.call(rbind, regions), file.path(output, "support_regions.csv"), row.names = FALSE)
write.csv(do.call(rbind, ledger), file.path(output, "ledger.csv"), row.names = FALSE)
write.csv(data.frame(patient_id = "3", observations = length(patient$y),
  positive_time_rows = sum(patient$obs_time > 0),
  baseline_cd4_rows = sum(patient$obs_time <= 0 & patient$ytype == 2L),
  censored_viral_rows = sum(patient$ytype == 1L & patient$cens == 1L),
  cd4_sd_depends_on_prediction = TRUE), file.path(output, "observation_model.csv"), row.names = FALSE)
write.csv(data.frame(positive_full_support_gaussian_lower_bound_ruled_out = TRUE,
  global_affine_forward_fixed_precision_error_bound_ruled_out = TRUE,
  cheap_gaussian_physical_support_regions_available = TRUE,
  useful_lower_likelihood_component_constructed = FALSE,
  sealed_numerical_likelihood_certificate_available = FALSE,
  residual_collapse_method_as_a_whole_falsified = FALSE, full_sampler_go = FALSE),
  file.path(output, "decision.csv"), row.names = FALSE)
source_paths <- c(modules, source_path, file.path(project,
  c("experiments/system_a_residual_bound_probe.R", "RESIDUAL_BOUND_FEASIBILITY.md")))
write.csv(data.frame(path = source_paths, sha256 = vapply(source_paths, .sab_system_a_sha256_file, "")),
          file.path(output, "sources.csv"), row.names = FALSE)
writeLines(capture.output(sessionInfo()), file.path(output, "session.txt"))
report <- c("System A residual-bound structural probe completed; no sampler fitted.",
  paste("Sealed prediction calls:", counter),
  paste("Elapsed seconds:", proc.time()[["elapsed"]] - started),
  "Unmodified positive Gaussian and global affine-error recipes are inapplicable.",
  "Simple Gaussian-integrable regions repair mathematical equilibrium support only.",
  "No useful likelihood lower bound, residual-gap estimate or full sampler go is claimed.")
writeLines(report, file.path(output, "REPORT.txt"))
paths <- list.files(output, full.names = TRUE)
write.csv(data.frame(file = basename(paths), sha256 = vapply(paths, .sab_system_a_sha256_file, "")),
          file.path(output, "manifest.csv"), row.names = FALSE)
writeLines("completed", file.path(output, "COMPLETED"))
cat(paste(report, collapse = "\n"), "\n")
