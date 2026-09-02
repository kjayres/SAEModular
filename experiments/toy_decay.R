#!/usr/bin/env Rscript

# Minimal falsification experiment for SAEM-anchored affine transport.
#
# Model:
#   x_i ~ N(0, tau^2)
#   kappa ~ N(kappa_prior_mean, kappa_prior_sd^2)
#   y_ij | x_i, kappa ~ N(exp(x_i - exp(kappa) * time_j), sigma^2)
#
# The two distinct observation times make the patient conditional nonlinear
# and prevent any scalar affine transformation from preserving both fitted
# observations.  The affine proposal is nevertheless exact because every
# global proposal is corrected with the original joint density and its
# Jacobian.

options(warn = 1)

command_args <- commandArgs(trailingOnly = TRUE)
output_dir <- if (length(command_args) >= 1L) command_args[[1L]] else {
  file.path("outputs", "toy_decay_v1")
}
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path("R", "affine_map.R"), local = FALSE)
source(file.path("R", "transport_mh.R"), local = FALSE)

config <- list(
  schema_version = "sab_toy_decay_config_v1",
  data_seed = 20260902L,
  n_patient = 24L,
  times = c(0.5, 1.5),
  tau = 0.6,
  sigma = 0.4,
  kappa_true = log(0.8),
  kappa_prior_mean = log(0.8),
  kappa_prior_sd = 0.7,
  saem_seed = 61341L,
  saem_iterations = 5000L,
  saem_exploration = 1000L,
  saem_gain_exponent = 0.7,
  saem_pcn_beta = 0.30,
  map_seed = 73001L,
  map_half_width = 0.15,
  map_warmup = 1000L,
  map_draws = 3000L,
  map_thin = 2L,
  map_pcn_beta = 0.30,
  # The first run used 0.5, 1.0, 1.5. Both methods selected the upper
  # boundary, so the frozen robustness rerun extends the same grid rather than
  # treating a boundary optimum as tuned.
  proposal_multipliers = c(0.5, 1.0, 1.5, 2.0, 2.5),
  tuning_seed = 84001L,
  tuning_iterations = 3500L,
  tuning_burn = 750L,
  evaluation_seed = 95001L,
  evaluation_chains = 4L,
  evaluation_iterations = 12000L,
  evaluation_burn = 3000L,
  local_pcn_beta = 0.30,
  quadrature_orders = c(80L, 160L),
  reference_grid = seq(-2.5, 1.5, length.out = 2401L)
)

sab_log_sum_exp <- function(value) {
  maximum <- max(value)
  if (!is.finite(maximum)) return(maximum)
  maximum + log(sum(exp(value - maximum)))
}

sab_toy_patient_loglik <- function(x, kappa, data, config) {
  if (length(x) != nrow(data) || !is.finite(kappa) ||
      any(!is.finite(x))) {
    return(rep(-Inf, nrow(data)))
  }
  rate <- exp(kappa)
  amplitude <- exp(x)
  if (!is.finite(rate) || any(!is.finite(amplitude))) {
    return(rep(-Inf, nrow(data)))
  }
  decay <- exp(-rate * config$times)
  prediction <- amplitude %o% decay
  contributions <- matrix(
    stats::dnorm(
      as.numeric(data),
      mean = as.numeric(prediction),
      sd = config$sigma,
      log = TRUE
    ),
    nrow = nrow(data),
    ncol = ncol(data)
  )
  value <- rowSums(contributions)
  value[!is.finite(value)] <- -Inf
  value
}

sab_toy_log_target <- function(global, locals, data, config) {
  kappa <- as.numeric(global)
  if (length(kappa) != 1L || !is.finite(kappa)) {
    return(list(log_density = -Inf, patient_loglik = rep(-Inf, nrow(data))))
  }
  x <- vapply(locals, function(value) as.numeric(value)[[1L]], numeric(1L))
  if (length(x) != nrow(data) || any(!is.finite(x))) {
    return(list(log_density = -Inf, patient_loglik = rep(-Inf, nrow(data))))
  }
  patient_loglik <- sab_toy_patient_loglik(x, kappa, data, config)
  sab_toy_cached_evaluation(kappa, x, patient_loglik, config)
}

sab_toy_cached_evaluation <- function(kappa, x, patient_loglik, config) {
  log_density <- stats::dnorm(
    kappa,
    config$kappa_prior_mean,
    config$kappa_prior_sd,
    log = TRUE
  ) + sum(stats::dnorm(x, 0, config$tau, log = TRUE)) +
    sum(patient_loglik)
  if (!is.finite(log_density)) log_density <- -Inf
  list(log_density = log_density, patient_loglik = patient_loglik)
}

sab_gauss_hermite <- function(order) {
  order <- as.integer(order)
  if (length(order) != 1L || is.na(order) || order < 2L) {
    stop("Gauss-Hermite order must be an integer of at least two.",
         call. = FALSE)
  }
  jacobi <- matrix(0, nrow = order, ncol = order)
  off_diagonal <- sqrt(seq_len(order - 1L) / 2)
  jacobi[cbind(seq_len(order - 1L), 2:order)] <- off_diagonal
  jacobi[cbind(2:order, seq_len(order - 1L))] <- off_diagonal
  decomposition <- eigen(jacobi, symmetric = TRUE)
  list(
    nodes = sqrt(2) * decomposition$values,
    weights = decomposition$vectors[1L, ]^2
  )
}

sab_reference_posterior <- function(order, grid, data, config) {
  rule <- sab_gauss_hermite(order)
  x_nodes <- config$tau * rule$nodes
  log_weights <- log(rule$weights)
  log_posterior <- vapply(grid, function(kappa) {
    prediction <- exp(x_nodes) %o% exp(-exp(kappa) * config$times)
    patient_terms <- vapply(seq_len(nrow(data)), function(patient) {
      node_loglik <- rowSums(matrix(
        stats::dnorm(
          rep(data[patient, ], each = length(x_nodes)),
          mean = as.numeric(prediction),
          sd = config$sigma,
          log = TRUE
        ),
        nrow = length(x_nodes),
        ncol = ncol(data)
      ))
      sab_log_sum_exp(log_weights + node_loglik)
    }, numeric(1L))
    stats::dnorm(
      kappa,
      config$kappa_prior_mean,
      config$kappa_prior_sd,
      log = TRUE
    ) + sum(patient_terms)
  }, numeric(1L))
  weights <- exp(log_posterior - max(log_posterior))
  weights <- weights / sum(weights)
  mean <- sum(grid * weights)
  variance <- sum((grid - mean)^2 * weights)
  cumulative <- cumsum(weights)
  quantile <- function(probability) {
    stats::approx(
      x = cumulative,
      y = grid,
      xout = probability,
      ties = "ordered",
      rule = 2
    )$y
  }
  list(
    order = order,
    grid = grid,
    log_posterior = log_posterior,
    weights = weights,
    mean = mean,
    sd = sqrt(variance),
    mode = grid[[which.max(log_posterior)]],
    quantiles = setNames(
      vapply(c(0.025, 0.1, 0.5, 0.9, 0.975), quantile, numeric(1L)),
      c("q025", "q10", "q50", "q90", "q975")
    ),
    edge_mass = sum(weights[c(seq_len(10L), length(weights) - 0:9L)])
  )
}

sab_pcn_patient_step <- function(x, current_loglik, kappa, data, config,
                                 beta) {
  proposed_x <- sqrt(1 - beta^2) * x +
    beta * stats::rnorm(length(x), 0, config$tau)
  proposed_loglik <- sab_toy_patient_loglik(
    proposed_x, kappa, data, config
  )
  accepted <- log(stats::runif(length(x))) <
    proposed_loglik - current_loglik
  x[accepted] <- proposed_x[accepted]
  current_loglik[accepted] <- proposed_loglik[accepted]
  list(x = x, loglik = current_loglik, accepted = accepted)
}

sab_run_saem <- function(data, config) {
  set.seed(config$saem_seed)
  kappa <- config$kappa_prior_mean
  initial_level <- pmax(rowMeans(data), 0.05)
  x <- log(initial_level) + exp(kappa) * mean(config$times)
  x <- pmax(pmin(x, 3 * config$tau), -3 * config$tau)
  current_loglik <- sab_toy_patient_loglik(x, kappa, data, config)
  statistic_y_exp_x <- numeric(length(config$times))
  statistic_exp_2x <- 0
  trace <- numeric(config$saem_iterations)
  accepted <- 0L

  for (iteration in seq_len(config$saem_iterations)) {
    step <- sab_pcn_patient_step(
      x, current_loglik, kappa, data, config, config$saem_pcn_beta
    )
    x <- step$x
    current_loglik <- step$loglik
    accepted <- accepted + sum(step$accepted)

    realised_y_exp_x <- colSums(data * exp(x))
    realised_exp_2x <- sum(exp(2 * x))
    gain <- if (iteration <= config$saem_exploration) {
      1
    } else {
      (iteration - config$saem_exploration)^(-config$saem_gain_exponent)
    }
    statistic_y_exp_x <- statistic_y_exp_x +
      gain * (realised_y_exp_x - statistic_y_exp_x)
    statistic_exp_2x <- statistic_exp_2x +
      gain * (realised_exp_2x - statistic_exp_2x)

    objective <- function(candidate_kappa) {
      decay <- exp(-exp(candidate_kappa) * config$times)
      expected_log_joint <-
        (sum(decay * statistic_y_exp_x) -
           0.5 * statistic_exp_2x * sum(decay^2)) / config$sigma^2 +
        stats::dnorm(
          candidate_kappa,
          config$kappa_prior_mean,
          config$kappa_prior_sd,
          log = TRUE
        )
      -expected_log_joint
    }
    kappa <- stats::optimize(objective, interval = c(-2.5, 1.5))$minimum
    # The next patient transition targets a different conditional after the
    # M-step.  Refresh the cached current likelihood under that new kappa;
    # comparing it with a proposal evaluated at the new kappa would otherwise
    # be an invalid Metropolis ratio.
    current_loglik <- sab_toy_patient_loglik(x, kappa, data, config)
    trace[[iteration]] <- kappa
  }

  list(
    kappa = kappa,
    x = x,
    trace = trace,
    patient_acceptance = accepted /
      (config$n_patient * config$saem_iterations),
    patient_likelihood_evaluations = config$n_patient *
      (2L * config$saem_iterations + 1L)
  )
}

sab_conditional_moments <- function(kappa, initial_x, data, config, seed) {
  set.seed(seed)
  x <- initial_x
  current_loglik <- sab_toy_patient_loglik(x, kappa, data, config)
  kept <- matrix(
    NA_real_, nrow = config$n_patient, ncol = config$map_draws
  )
  accepted <- 0L
  total_iterations <- config$map_warmup +
    config$map_draws * config$map_thin
  kept_index <- 0L
  for (iteration in seq_len(total_iterations)) {
    step <- sab_pcn_patient_step(
      x, current_loglik, kappa, data, config, config$map_pcn_beta
    )
    x <- step$x
    current_loglik <- step$loglik
    accepted <- accepted + sum(step$accepted)
    after_warmup <- iteration > config$map_warmup
    on_thinning_grid <- after_warmup &&
      ((iteration - config$map_warmup) %% config$map_thin == 0L)
    if (on_thinning_grid) {
      kept_index <- kept_index + 1L
      kept[, kept_index] <- x
    }
  }
  if (kept_index != config$map_draws || any(!is.finite(kept))) {
    stop("Conditional moment bank was not filled with finite draws.",
         call. = FALSE)
  }
  list(
    mean = rowMeans(kept),
    sd = apply(kept, 1L, stats::sd),
    acceptance = accepted / (config$n_patient * total_iterations),
    last_x = x,
    patient_likelihood_evaluations = config$n_patient *
      (total_iterations + 1L)
  )
}

sab_fit_frozen_maps <- function(anchor, initial_x, data, config) {
  design <- anchor + c(-config$map_half_width, 0, config$map_half_width)
  banks <- vector("list", length(design))
  for (index in seq_along(design)) {
    banks[[index]] <- sab_conditional_moments(
      kappa = design[[index]],
      initial_x = initial_x,
      data = data,
      config = config,
      seed = config$map_seed + index
    )
  }
  mean_at_design <- do.call(cbind, lapply(banks, `[[`, "mean"))
  sd_at_design <- do.call(cbind, lapply(banks, `[[`, "sd"))
  if (any(!is.finite(sd_at_design)) || any(sd_at_design <= 0)) {
    stop("Conditional map scales must be finite and positive.", call. = FALSE)
  }
  design_matrix <- cbind(1, design - anchor)
  mean_coefficients <- t(stats::lm.fit(
    x = design_matrix,
    y = t(mean_at_design)
  )$coefficients)
  log_sd_coefficients <- t(stats::lm.fit(
    x = design_matrix,
    y = t(log(sd_at_design))
  )$coefficients)

  patient_ids <- sprintf("patient_%03d", seq_len(config$n_patient))
  maps <- setNames(lapply(seq_len(config$n_patient), function(patient) {
    force(patient)
    sab_new_affine_map(
      mean_fn = function(global) {
        delta <- as.numeric(global)[[1L]] - anchor
        delta <- max(-config$map_half_width,
                     min(config$map_half_width, delta))
        mean_coefficients[patient, 1L] +
          mean_coefficients[patient, 2L] * delta
      },
      chol_fn = function(global) {
        delta <- as.numeric(global)[[1L]] - anchor
        delta <- max(-config$map_half_width,
                     min(config$map_half_width, delta))
        value <- exp(
          log_sd_coefficients[patient, 1L] +
            log_sd_coefficients[patient, 2L] * delta
        )
        matrix(value, nrow = 1L, ncol = 1L)
      },
      name = patient_ids[[patient]],
      check_at = c(kappa = anchor)
    )
  }), patient_ids)

  list(
    schema_version = "sab_frozen_toy_map_v1",
    anchor = anchor,
    design = design,
    mean_coefficients = mean_coefficients,
    log_sd_coefficients = log_sd_coefficients,
    maps = maps,
    bank_acceptance = vapply(banks, `[[`, numeric(1L), "acceptance"),
    patient_likelihood_evaluations = sum(vapply(
      banks, `[[`, numeric(1L), "patient_likelihood_evaluations"
    ))
  )
}

sab_fixed_maps <- function(n_patient) {
  patient_ids <- sprintf("patient_%03d", seq_len(n_patient))
  setNames(lapply(patient_ids, function(patient_id) {
    sab_new_affine_map(
      mean_fn = function(global) 0,
      chol_fn = function(global) matrix(1, 1L, 1L),
      name = paste0(patient_id, "_fixed"),
      check_at = c(kappa = 0)
    )
  }), patient_ids)
}

sab_initial_state <- function(kappa, x, data, config) {
  patient_ids <- sprintf("patient_%03d", seq_len(config$n_patient))
  locals <- setNames(lapply(x, function(value) value), patient_ids)
  global <- c(kappa = kappa)
  evaluation <- sab_toy_log_target(global, locals, data, config)
  if (!is.finite(evaluation$log_density)) {
    stop("Initial toy state has non-finite density.", call. = FALSE)
  }
  list(global = global, locals = locals, evaluation = evaluation)
}

sab_run_exact_chain <- function(method, maps, proposal_sd, initial_kappa,
                                initial_x, iterations, seed, data, config) {
  if (!method %in% c("fixed", "affine")) {
    stop("Unknown toy method.", call. = FALSE)
  }
  set.seed(seed)
  state <- sab_initial_state(initial_kappa, initial_x, data, config)
  trace <- numeric(iterations)
  global_accepted <- logical(iterations)
  local_acceptance <- numeric(iterations)
  patient_ids <- names(state$locals)
  target <- function(global, locals) {
    sab_toy_log_target(global, locals, data, config)
  }
  symmetric_proposal <- function(to, from) 0

  for (iteration in seq_len(iterations)) {
    current_x <- vapply(
      state$locals, function(value) as.numeric(value)[[1L]], numeric(1L)
    )
    local_step <- sab_pcn_patient_step(
      x = current_x,
      current_loglik = state$evaluation$patient_loglik,
      kappa = as.numeric(state$global),
      data = data,
      config = config,
      beta = config$local_pcn_beta
    )
    state$locals <- setNames(
      lapply(local_step$x, function(value) value), patient_ids
    )
    state$evaluation <- sab_toy_cached_evaluation(
      kappa = as.numeric(state$global),
      x = local_step$x,
      patient_loglik = local_step$loglik,
      config = config
    )
    local_acceptance[[iteration]] <- mean(local_step$accepted)

    proposed_global <- c(
      kappa = as.numeric(state$global) + stats::rnorm(1L, 0, proposal_sd)
    )
    proposal <- sab_transport_mh_proposal(
      current_global = state$global,
      proposed_global = proposed_global,
      current_locals = state$locals,
      maps = maps,
      log_target = target,
      log_global_proposal = symmetric_proposal,
      current_evaluation = state$evaluation
    )
    decision <- sab_decide_transport_mh(proposal)
    state$global <- decision$state$global
    state$locals <- decision$state$locals
    state$evaluation <- decision$evaluation
    global_accepted[[iteration]] <- decision$accepted
    trace[[iteration]] <- as.numeric(state$global)
  }

  list(
    method = method,
    trace = trace,
    global_acceptance = mean(global_accepted),
    local_acceptance = mean(local_acceptance),
    patient_likelihood_evaluations = config$n_patient *
      (1L + 2L * iterations),
    patient_likelihood_evaluations_per_iteration =
      2L * config$n_patient,
    final_state = state
  )
}

sab_initial_positive_ess <- function(value) {
  value <- as.numeric(value)
  n <- length(value)
  if (n < 4L || !is.finite(stats::var(value)) || stats::var(value) == 0) {
    return(NA_real_)
  }
  lag_max <- min(n - 1L, max(100L, floor(n / 2L)))
  correlation <- as.numeric(stats::acf(
    value, lag.max = lag_max, plot = FALSE, demean = TRUE
  )$acf)[-1L]
  if (length(correlation) %% 2L == 1L) {
    correlation <- correlation[-length(correlation)]
  }
  pair_sums <- correlation[seq.int(1L, length(correlation), by = 2L)] +
    correlation[seq.int(2L, length(correlation), by = 2L)]
  first_nonpositive <- which(!is.finite(pair_sums) | pair_sums <= 0)[1L]
  if (is.na(first_nonpositive)) {
    retained_pairs <- seq_along(pair_sums)
  } else if (first_nonpositive == 1L) {
    retained_pairs <- integer()
  } else {
    retained_pairs <- seq_len(first_nonpositive - 1L)
  }
  tau <- 1 + 2 * sum(pair_sums[retained_pairs])
  min(n, n / max(tau, 1))
}

sab_split_rhat <- function(chains) {
  chain_length <- min(vapply(chains, length, integer(1L)))
  half <- floor(chain_length / 2L)
  if (half < 2L) return(NA_real_)
  split <- c(
    lapply(chains, function(value) value[seq_len(half)]),
    lapply(chains, function(value) tail(value, half))
  )
  matrix_value <- do.call(cbind, split)
  within <- mean(apply(matrix_value, 2L, stats::var))
  between <- half * stats::var(colMeans(matrix_value))
  variance_hat <- (half - 1) / half * within + between / half
  sqrt(variance_hat / within)
}

sab_chain_summary <- function(chains, burn, reference) {
  retained <- lapply(chains, function(chain) {
    chain$trace[seq.int(burn + 1L, length(chain$trace))]
  })
  pooled <- unlist(retained, use.names = FALSE)
  ess_by_chain <- vapply(retained, sab_initial_positive_ess, numeric(1L))
  squared_jumps <- lapply(retained, function(value) {
    (diff(value) / reference$sd)^2
  })
  posterior_scaled_esjd <- mean(unlist(squared_jumps, use.names = FALSE))
  total_posterior_scaled_squared_jump <- sum(vapply(
    squared_jumps, sum, numeric(1L)
  ))
  total_ess <- sum(ess_by_chain)
  total_calls <- sum(vapply(
    chains, `[[`, numeric(1L), "patient_likelihood_evaluations"
  ))
  retained_sampling_calls <- sum(vapply(chains, function(chain) {
    chain$patient_likelihood_evaluations_per_iteration *
      (length(chain$trace) - burn)
  }, numeric(1L)))
  list(
    mean = mean(pooled),
    sd = stats::sd(pooled),
    quantiles = stats::quantile(
      pooled, probs = c(0.025, 0.5, 0.975), names = FALSE
    ),
    split_rhat = sab_split_rhat(retained),
    ess_by_chain = ess_by_chain,
    total_ess = total_ess,
    total_calls = total_calls,
    retained_sampling_calls = retained_sampling_calls,
    ess_per_million_calls = 1e6 * total_ess / total_calls,
    retained_ess_per_million_calls =
      1e6 * total_ess / retained_sampling_calls,
    posterior_scaled_esjd = posterior_scaled_esjd,
    total_posterior_scaled_squared_jump =
      total_posterior_scaled_squared_jump,
    esjd_per_million_calls =
      1e6 * total_posterior_scaled_squared_jump / total_calls,
    mean_error_in_reference_sd =
      abs(mean(pooled) - reference$mean) / reference$sd,
    sd_relative_error = abs(stats::sd(pooled) / reference$sd - 1),
    mean_mcse = stats::sd(pooled) / sqrt(total_ess),
    global_acceptance = mean(vapply(
      chains, `[[`, numeric(1L), "global_acceptance"
    )),
    local_acceptance = mean(vapply(
      chains, `[[`, numeric(1L), "local_acceptance"
    ))
  )
}

set.seed(config$data_seed)
latent_truth <- stats::rnorm(config$n_patient, 0, config$tau)
truth_prediction <- exp(latent_truth) %o%
  exp(-exp(config$kappa_true) * config$times)
data <- matrix(
  stats::rnorm(
    length(truth_prediction),
    mean = as.numeric(truth_prediction),
    sd = config$sigma
  ),
  nrow = config$n_patient,
  ncol = length(config$times)
)

references <- lapply(config$quadrature_orders, function(order) {
  sab_reference_posterior(
    order = order,
    grid = config$reference_grid,
    data = data,
    config = config
  )
})
reference <- references[[length(references)]]
quadrature_mean_change <- abs(references[[1L]]$mean - reference$mean) /
  reference$sd
quadrature_sd_change <- abs(references[[1L]]$sd / reference$sd - 1)
quadrature_stable <- quadrature_mean_change < 0.01 &&
  quadrature_sd_change < 0.01 && reference$edge_mass < 1e-6

saem <- sab_run_saem(data, config)
saem_error_in_reference_sd <- abs(saem$kappa - reference$mode) / reference$sd

frozen_map <- sab_fit_frozen_maps(
  anchor = saem$kappa,
  initial_x = saem$x,
  data = data,
  config = config
)
map_artifact <- frozen_map
map_artifact$maps <- NULL
saveRDS(map_artifact, file.path(output_dir, "frozen_map.rds"), version = 3)

maps_by_method <- list(
  fixed = sab_fixed_maps(config$n_patient),
  affine = frozen_map$maps
)

tuning_rows <- list()
selected_multiplier <- numeric(length(maps_by_method))
names(selected_multiplier) <- names(maps_by_method)
for (method in names(maps_by_method)) {
  method_rows <- vector("list", length(config$proposal_multipliers))
  for (scale_index in seq_along(config$proposal_multipliers)) {
    multiplier <- config$proposal_multipliers[[scale_index]]
    chain <- sab_run_exact_chain(
      method = method,
      maps = maps_by_method[[method]],
      proposal_sd = multiplier * reference$sd,
      initial_kappa = saem$kappa,
      initial_x = saem$x,
      iterations = config$tuning_iterations,
      seed = config$tuning_seed +
        match(method, names(maps_by_method)) * 100L + scale_index,
      data = data,
      config = config
    )
    retained_trace <- chain$trace[seq.int(
      config$tuning_burn + 1L, length(chain$trace)
    )]
    standardised_jump <- diff(retained_trace) / reference$sd
    method_rows[[scale_index]] <- data.frame(
      method = method,
      multiplier = multiplier,
      acceptance = chain$global_acceptance,
      posterior_scaled_esjd = mean(standardised_jump^2),
      stringsAsFactors = FALSE
    )
  }
  method_table <- do.call(rbind, method_rows)
  selected_multiplier[[method]] <- method_table$multiplier[
    which.max(method_table$posterior_scaled_esjd)
  ]
  tuning_rows[[method]] <- method_table
}
tuning <- do.call(rbind, tuning_rows)

evaluation <- list()
initial_probabilities <- seq(
  0.1, 0.9, length.out = config$evaluation_chains
)
initial_kappas <- vapply(initial_probabilities, function(probability) {
  stats::approx(
    x = cumsum(reference$weights),
    y = reference$grid,
    xout = probability,
    ties = "ordered",
    rule = 2
  )$y
}, numeric(1L))

for (method in names(maps_by_method)) {
  chains <- vector("list", config$evaluation_chains)
  for (chain_index in seq_len(config$evaluation_chains)) {
    initial_kappa <- initial_kappas[[chain_index]]
    initial_x <- vapply(maps_by_method$affine, function(map) {
      sab_affine_map_components(map, c(kappa = initial_kappa))$mean[[1L]]
    }, numeric(1L))
    chains[[chain_index]] <- sab_run_exact_chain(
      method = method,
      maps = maps_by_method[[method]],
      proposal_sd = selected_multiplier[[method]] * reference$sd,
      initial_kappa = initial_kappa,
      initial_x = initial_x,
      iterations = config$evaluation_iterations,
      seed = config$evaluation_seed +
        match(method, names(maps_by_method)) * 1000L + chain_index,
      data = data,
      config = config
    )
  }
  evaluation[[method]] <- list(
    chains = chains,
    summary = sab_chain_summary(chains, config$evaluation_burn, reference)
  )
}

efficiency_ratio <- evaluation$affine$summary$ess_per_million_calls /
  evaluation$fixed$summary$ess_per_million_calls
esjd_ratio <- evaluation$affine$summary$posterior_scaled_esjd /
  evaluation$fixed$summary$posterior_scaled_esjd
affine_correct <- isTRUE(
  evaluation$affine$summary$split_rhat <= 1.01 &&
  evaluation$affine$summary$total_ess >= 400 &&
  evaluation$affine$summary$mean_error_in_reference_sd <= 0.10 &&
  evaluation$affine$summary$sd_relative_error <= 0.10
)
fixed_adequate <- isTRUE(
  evaluation$fixed$summary$split_rhat <= 1.01 &&
  evaluation$fixed$summary$total_ess >= 400 &&
  evaluation$fixed$summary$mean_error_in_reference_sd <= 0.10 &&
  evaluation$fixed$summary$sd_relative_error <= 0.10
)
mechanism_pass <- isTRUE(quadrature_stable &&
  saem_error_in_reference_sd <= 0.50 &&
  affine_correct && fixed_adequate &&
  efficiency_ratio >= 1.50 && esjd_ratio >= 1.25)

map_cost <- frozen_map$patient_likelihood_evaluations
fixed_rate <- evaluation$fixed$summary$total_ess /
  evaluation$fixed$summary$retained_sampling_calls
affine_rate <- evaluation$affine$summary$total_ess /
  evaluation$affine$summary$retained_sampling_calls
sab_break_even <- function(startup_calls) if (
    is.finite(affine_rate) && is.finite(fixed_rate) &&
    affine_rate > fixed_rate && fixed_rate > 0) {
  startup_calls / (1 / fixed_rate - 1 / affine_rate)
} else {
  Inf
}
map_break_even_effective_samples <- sab_break_even(map_cost)
total_startup_break_even_effective_samples <- sab_break_even(
  saem$patient_likelihood_evaluations + map_cost
)

summary_table <- do.call(rbind, lapply(names(evaluation), function(method) {
  value <- evaluation[[method]]$summary
  data.frame(
    method = method,
    selected_proposal_multiplier = selected_multiplier[[method]],
    global_acceptance = value$global_acceptance,
    local_acceptance = value$local_acceptance,
    mean = value$mean,
    sd = value$sd,
    split_rhat = value$split_rhat,
    total_ess = value$total_ess,
    ess_per_million_patient_likelihood_calls = value$ess_per_million_calls,
    retained_ess_per_million_patient_likelihood_calls =
      value$retained_ess_per_million_calls,
    posterior_scaled_esjd = value$posterior_scaled_esjd,
    esjd_per_million_patient_likelihood_calls = value$esjd_per_million_calls,
    mean_error_in_reference_sd = value$mean_error_in_reference_sd,
    sd_relative_error = value$sd_relative_error,
    stringsAsFactors = FALSE
  )
}))

write.csv(tuning, file.path(output_dir, "tuning.csv"), row.names = FALSE)
write.csv(summary_table, file.path(output_dir, "chain_summary.csv"),
          row.names = FALSE)
saveRDS(
  list(
    schema_version = "sab_toy_decay_result_v1",
    config = config,
    data = data,
    latent_truth = latent_truth,
    references = references,
    reference = reference,
    saem = saem,
    frozen_map = map_artifact,
    tuning = tuning,
    selected_multiplier = selected_multiplier,
    evaluation = evaluation,
    efficiency_ratio = efficiency_ratio,
    esjd_ratio = esjd_ratio,
    map_break_even_effective_samples = map_break_even_effective_samples,
    total_startup_break_even_effective_samples =
      total_startup_break_even_effective_samples,
    gates = list(
      quadrature_stable = quadrature_stable,
      affine_correct = affine_correct,
      fixed_adequate = fixed_adequate,
      mechanism_pass = mechanism_pass
    )
  ),
  file.path(output_dir, "result.rds"),
  version = 3
)

report <- c(
  "SAEM-anchored affine transport: nonlinear decay toy",
  sprintf("mechanism_pass: %s", mechanism_pass),
  sprintf("quadrature_stable: %s", quadrature_stable),
  sprintf("reference mean / sd: %.6f / %.6f", reference$mean, reference$sd),
  sprintf("SAEM anchor: %.6f (%.3f reference SD from mode)",
          saem$kappa, saem_error_in_reference_sd),
  sprintf("SAEM patient acceptance: %.3f", saem$patient_acceptance),
  sprintf("map-bank acceptance at -h,0,+h: %s",
          paste(sprintf("%.3f", frozen_map$bank_acceptance), collapse = ", ")),
  sprintf("fixed ESS/million calls: %.1f",
          evaluation$fixed$summary$ess_per_million_calls),
  sprintf("affine ESS/million calls: %.1f",
          evaluation$affine$summary$ess_per_million_calls),
  sprintf("affine/fixed efficiency ratio: %.3f", efficiency_ratio),
  sprintf("affine/fixed posterior-scaled ESJD ratio: %.3f", esjd_ratio),
  sprintf("affine split-Rhat: %.4f", evaluation$affine$summary$split_rhat),
  sprintf("affine mean error (reference SD): %.4f",
          evaluation$affine$summary$mean_error_in_reference_sd),
  sprintf("affine SD relative error: %.4f",
          evaluation$affine$summary$sd_relative_error),
  sprintf("one-off SAEM calls: %d", saem$patient_likelihood_evaluations),
  sprintf("one-off map calls: %d", map_cost),
  sprintf("map-only break-even effective samples (SAEM already available): %s",
          if (is.finite(map_break_even_effective_samples)) {
            sprintf("%.1f", map_break_even_effective_samples)
          } else {
            "never at measured rates"
          }),
  sprintf("total-startup break-even effective samples (versus no anchor cost): %s",
          if (is.finite(total_startup_break_even_effective_samples)) {
            sprintf("%.1f", total_startup_break_even_effective_samples)
          } else {
            "never at measured rates"
          })
)
writeLines(report, file.path(output_dir, "REPORT.txt"))
writeLines(report)
