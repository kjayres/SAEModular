#!/usr/bin/env Rscript
# A finite-cost falsification screen; no sampler or certified lower bound.
options(warn = 1)
if (!nzchar(Sys.getenv("SLURM_JOB_ID"))) stop("Run through Slurm.")
args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1L) stop("Supply a new output directory.")
project <- normalizePath(Sys.getenv("SAB_PROJECT_ROOT"), mustWork = TRUE)
workspace <- normalizePath(file.path(project, "..", ".."), mustWork = TRUE)
output <- args[[1L]]
if (file.exists(output) || !dir.create(output, recursive = TRUE)) stop("Refusing output directory.")
modules <- file.path(project, "R", c("system_a_adapter.R", "shared_trajectory_probe.R", "residual_ceiling.R"))
invisible(lapply(modules, sys.source, envir = globalenv()))
artifact_path <- file.path(project, "outputs/system_a_12_oracle_rescue_v2/frozen_map.rds")
stopifnot(identical(.sab_system_a_sha256_file(artifact_path),
  "ec4b29a4efc7162da43f2295c0530d1f7a02e72e1c785df1cd0a7e3b6b590a78"))
artifact <- readRDS(artifact_path)
ids <- c("3", "74", "122")
adapter <- sab_load_system_a_adapter(workspace, ids)
stopifnot(identical(artifact$target_fingerprint, adapter$target_fingerprint))
eta <- artifact$fixed_eta
psi <- c(artifact$dynamic_center, artifact$fixed_observation_psi)[adapter$coordinate_names$global]
local_names <- adapter$coordinate_names$local
chart_names <- c("log_equilibrium_excess", local_names[-1L])
owner <- environment(adapter$solve_prediction)
exact <- get("exact_likelihood", owner, inherits = FALSE)
upstream <- get("upstream", owner, inherits = FALSE)

softplus <- function(v) pmax(v, 0) + log1p(exp(-abs(v)))
logadd <- function(a, b) max(a, b) + log1p(exp(-abs(a - b)))
critical <- function(x, mu_v) {
  log_pi <- stats::plogis(x[[6L]], log.p = TRUE)
  x[[2L]] + x[[3L]] + log(mu_v) - psi[[1L]] - x[[4L]] +
    logadd(x[[5L]], psi[[2L]]) - logadd(x[[5L]], log_pi + psi[[2L]])
}
to_chart <- function(x, mu_v) {
  delta <- x[[1L]] - critical(x, mu_v)
  if (!is.finite(delta) || delta <= 0) stop("Saved component centre is infeasible.")
  y <- x; y[[1L]] <- delta + log1p(-exp(-delta))
  setNames(y, chart_names)
}
from_chart <- function(y, mu_v) {
  x <- setNames(y, local_names)
  x[[1L]] <- critical(x, mu_v) + softplus(y[[1L]])
  x
}
chart_jacobian <- function(x, y) {
  pi <- stats::plogis(x[[6L]])
  grad <- numeric(8L); grad[2:4] <- c(1, 1, -1)
  grad[[5L]] <- stats::plogis(x[[5L]] - psi[[2L]]) -
    stats::plogis(x[[5L]] - log(pi) - psi[[2L]])
  grad[[6L]] <- -(1 - pi) * stats::plogis(log(pi) + psi[[2L]] - x[[5L]])
  jac <- diag(8L); jac[1L, ] <- (jac[1L, ] - grad) / stats::plogis(y[[1L]])
  jac
}

masses <- c(.5, .8, .95)
caps <- c(2, 10, 50)
alpha <- .05 / (length(ids) * length(masses) * length(caps))
ledger <- rows <- states <- components <- list()
call_count <- 0L
started <- proc.time()
evaluate <- function(patient_id, y, mu_v, phase, mass, index) {
  call_count <<- call_count + 1L
  if (call_count > 2500L) stop("Predeclared prediction-call cap exceeded.")
  x <- from_chart(y, mu_v)
  begin <- proc.time()[["elapsed"]]
  prediction <- adapter$solve_prediction(patient_id, x, psi)
  ll <- adapter$loglik_from_prediction(prediction, psi)
  log_g <- adapter$log_population_density(patient_id, x, eta)
  log_jac <- -softplus(-y[[1L]])
  log_f <- ll + log_g + log_jac
  if (is.na(log_f) || identical(log_f, Inf)) stop("Invalid complete patient factor.")
  ledger[[call_count]] <<- data.frame(patient_id = patient_id, phase = phase,
    gaussian_mass = mass, index = index, prediction_calls = 1L,
    ode_integrations = .sab_shared_original_calls(prediction),
    ok = isTRUE(prediction$ok), reason = if (isTRUE(prediction$ok)) "ok" else prediction$reason,
    log_likelihood = ll, log_population = log_g, log_jacobian = log_jac,
    log_factor = log_f, elapsed_sec = proc.time()[["elapsed"]] - begin)
  states[[call_count]] <<- c(y, setNames(x, paste0("original_", names(x))))
  log_f
}

for (i in seq_along(ids)) {
  id <- ids[[i]]
  mu_v <- exact$patients[[id]]$controls$mu_v
  coef <- artifact$coefficients[[id]]
  x_mean <- setNames(coef$mean[, "intercept"], local_names)
  index <- coef$lower_index
  entries <- coef$chol[, "intercept"]
  diagonal <- index[, 1L] == index[, 2L]
  entries[diagonal] <- exp(entries[diagonal])
  x_chol <- matrix(0, 8L, 8L); x_chol[index] <- entries
  centre <- to_chart(x_mean, mu_v)
  jac <- chart_jacobian(x_mean, centre)
  propagated <- jac %*% x_chol
  scale <- t(chol(tcrossprod(propagated)))
  stopifnot(max(abs(from_chart(centre, mu_v) - x_mean)) < 1e-10)
  # Independent finite differences check the complete triangular chart Jacobian.
  step <- 1e-5
  fd <- vapply(seq_len(8L), function(j) {
    plus <- minus <- x_mean; plus[[j]] <- plus[[j]] + step; minus[[j]] <- minus[[j]] - step
    (to_chart(plus, mu_v) - to_chart(minus, mu_v)) / (2 * step)
  }, numeric(8L))
  stopifnot(max(abs(fd - jac)) < 1e-5,
    abs(as.numeric(determinant(jac, logarithm = TRUE)$modulus) - softplus(-centre[[1L]])) < 1e-10)
  natural <- upstream$sysa_local_to_natural(x_mean)
  shared <- upstream$sysa_psi_to_natural(psi)
  equilibrium <- upstream$sysa_stan_equilibrium(natural, shared, mu_v)
  stopifnot(!is.null(equilibrium), isTRUE(all.equal(equilibrium[["VI"]],
    natural[["mu_t"]] / shared[["gamma"]] * exp(centre[[1L]]), tolerance = 1e-9)))
  components[[id]] <- list(mean = centre, chol = scale, original_mean = x_mean,
    original_chol = x_chol, chart_jacobian = jac, finite_difference_error = max(abs(fd - jac)))

  # Freeze all probe multipliers before generating any validation state.
  log_probe <- numeric(length(masses))
  for (j in seq_along(masses)) {
    mass <- masses[[j]]
    set.seed(2026090500L + 100L * i + j)
    random <- sab_truncated_gaussian_draw(64L, centre, scale, mass)
    radius <- sqrt(stats::qchisq(mass, 8L))
    axes <- rbind(diag(8L), -diag(8L)) * (.999 * radius)
    probes <- rbind(centre, sweep(axes %*% t(scale), 2L, centre, "+"), random)
    log_q <- sab_truncated_gaussian_logdensity(probes, centre, scale, mass)
    stopifnot(all(is.finite(log_q)))
    log_f <- vapply(seq_len(nrow(probes)), function(k) {
      evaluate(id, probes[k, ], mu_v, "probe", mass, k)
    }, numeric(1L))
    log_probe[[j]] <- min(log_f - log_q)
  }
  set.seed(2026099500L + i)
  standard <- matrix(stats::rnorm(512L * 8L), 512L, 8L)
  validation <- sweep(standard %*% t(scale), 2L, centre, "+")
  log_h <- -.5 * (8L * log(2 * pi) + rowSums(standard^2)) - sum(log(diag(scale)))
  # 'pi' above is the mathematical constant, not the patient latent probability.
  log_f <- vapply(seq_len(nrow(validation)), function(k) {
    evaluate(id, validation[k, ], mu_v, "validation", NA_real_, k)
  }, numeric(1L))
  for (j in seq_along(masses)) {
    bound <- sab_residual_ceiling(log_f - log_h, log_probe[[j]], caps, alpha)
    bound$patient_id <- id
    bound$gaussian_mass <- masses[[j]]
    bound$radius <- sqrt(stats::qchisq(masses[[j]], 8L))
    bound$log_c_probe <- log_probe[[j]]
    bound$n_validation <- 512L
    bound$alpha_per_cap <- alpha
    rows[[length(rows) + 1L]] <- bound
  }
  cat("Completed patient", id, "prediction calls so far", call_count, "\n")
  flush.console()
}
bound_table <- do.call(rbind, rows)
decision <- do.call(rbind, lapply(split(bound_table,
  interaction(bound_table$patient_id, bound_table$gaussian_mass, drop = TRUE)), function(group) {
    upper <- min(group$collapse_upper)
    data.frame(patient_id = group$patient_id[[1L]], gaussian_mass = group$gaussian_mass[[1L]],
      simultaneous95_collapse_upper = upper, rules_out_80_percent = upper < .8,
      certified_lower_component = FALSE, sampler_go = FALSE)
  }))
write.csv(bound_table, file.path(output, "bounds.csv"), row.names = FALSE)
write.csv(decision, file.path(output, "decision.csv"), row.names = FALSE)
write.csv(do.call(rbind, ledger), file.path(output, "ledger.csv"), row.names = FALSE)
saveRDS(list(components = components, eta = eta, psi = psi, masses = masses, caps = caps,
  alpha_per_cap = alpha, artifact_sha256 = .sab_system_a_sha256_file(artifact_path),
  target_fingerprint = adapter$target_fingerprint, states = do.call(rbind, states)),
  file.path(output, "design_and_states.rds"))
source_paths <- c(modules, artifact_path, file.path(project,
  c("experiments/system_a_residual_ceiling.R", "RESIDUAL_CEILING_CONTRACT.md")))
write.csv(data.frame(path = source_paths, sha256 = vapply(source_paths, .sab_system_a_sha256_file, "")),
  file.path(output, "sources.csv"), row.names = FALSE)
elapsed <- proc.time() - started
report <- c("Bank-free full-eight-coordinate local-component rejection screen.",
  paste("Prediction calls:", call_count), paste("ODE integrations:", sum(vapply(ledger, function(r) r$ode_integrations, 0))),
  paste("Components ruled out for 80% collapse:", sum(decision$rules_out_80_percent), "of", nrow(decision)),
  paste("CPU seconds:", sum(elapsed[c("user.self", "sys.self")])),
  paste("Elapsed seconds:", elapsed[["elapsed"]]),
  "Simultaneous 95% upper confidence bounds; no deterministic likelihood certificate.",
  "High ceilings are inconclusive. Low ceilings reject only these frozen components.",
  "Oracle moment construction costs are historical and excluded; no deployable cost claim.",
  "No sampler, posterior approximation, SAEM, or full-cohort experiment was constructed.")
writeLines(report, file.path(output, "REPORT.txt"))
writeLines(capture.output(sessionInfo()), file.path(output, "session.txt"))
paths <- list.files(output, full.names = TRUE)
write.csv(data.frame(file = basename(paths), sha256 = vapply(paths, .sab_system_a_sha256_file, "")),
  file.path(output, "manifest.csv"), row.names = FALSE)
writeLines("completed", file.path(output, "COMPLETED"))
cat(paste(report, collapse = "\n"), "\n")
