# Reduced exact System A hierarchy, with an optional live auxiliary pool.
# The union-grid solver supplies proposals ONLY; active likelihoods always use
# the pinned adapter. Fixed q, distinct assignments and continuous refresh are
# essential. This is not a frozen-bank likelihood approximation.

sab_shared_pool_reference_components <- function(adapter, snapshots) {
  representatives <- vapply(0:1, function(t) {
    adapter$patient_ids[which(adapter$treatment == t)[1L]]
  }, character(1L))
  unlist(lapply(snapshots, function(s) lapply(representatives, function(id) {
    list(eta = s$eta, patient_id = id, mean = adapter$population_mean(id, s$eta),
         sd = setNames(exp(s$eta[10:17]), adapter$coordinate_names$local))
  })), recursive = FALSE)
}

.sab_chain_normal_logdensity <- function(x, mean, sd) {
  rowSums(matrix(stats::dnorm(as.numeric(x),
    mean = rep(as.numeric(mean), each = nrow(x)),
    sd = rep(as.numeric(sd), each = nrow(x)), log = TRUE), nrow = nrow(x)))
}

.sab_chain_log_q <- function(x, components) {
  values <- vapply(components, function(c) {
    .sab_chain_normal_logdensity(x, c$mean, c$sd)
  }, numeric(nrow(x)))
  if (is.null(dim(values))) values <- matrix(values, nrow = nrow(x))
  maximum <- apply(values, 1L, max)
  result <- maximum + log(rowMeans(exp(values - maximum)))
  if (any(!is.finite(result))) stop("Fixed reserve density is not representable.")
  result
}

sab_run_system_a_shared_chain <- function(
    adapter, probe, initial, reference_components, prior_mean, prior_sd,
    method = c("baseline", "recycling"), seed, config, checkpoint_path = NULL) {
  method <- match.arg(method)
  defaults <- list(budget = 100000L, warmup = 500L, dynamic_every = 5L,
    reserve_refresh = 4L, pool_size = 36L, direct_probability = .2,
    initial_pcn_beta = sqrt(.19), initial_dynamic_sd = c(.04, .04),
    initial_noise_sd = .06, log_floor = -1e6, checkpoint_every = 250L,
    max_sweeps = 20000L)
  if (length(setdiff(names(config), names(defaults)))) stop("Unknown chain configuration.")
  config <- utils::modifyList(defaults, config)
  for (key in c("budget", "warmup", "dynamic_every", "reserve_refresh",
                "pool_size", "checkpoint_every", "max_sweeps")) {
    v <- config[[key]]
    if (length(v) != 1L || !is.finite(v) || v < 1 || v != floor(v)) {
      stop("Invalid integer configuration: ", key)
    }
  }
  stopifnot(config$direct_probability > 0, config$direct_probability < 1,
    config$initial_pcn_beta > 0, config$initial_pcn_beta < 1,
    length(config$initial_dynamic_sd) == 2L, all(config$initial_dynamic_sd > 0),
    config$initial_noise_sd > 0, is.finite(config$log_floor))
  set.seed(seed)
  started <- proc.time()
  if (!is.character(adapter$target_fingerprint) || length(adapter$target_fingerprint) != 1L ||
      !identical(probe$reference_target_fingerprint, adapter$target_fingerprint)) {
    stop("Common proposal and exact adapter must share the pinned source target.")
  }
  ids <- adapter$patient_ids
  pool_ids <- names(probe$patients)[lengths(probe$positive_times) > 0L]
  if (!identical(names(probe$patients), pool_ids) ||
      !all(pool_ids %in% ids) || !length(pool_ids)) {
    stop("The common probe must contain only this panel's ODE patients.")
  }
  n <- length(ids); m <- length(pool_ids); K <- config$pool_size
  if (K <= m || config$reserve_refresh > K - m) stop("Too few unused pool slots.")
  if (config$budget < m) stop("Budget cannot cover exact initialization.")
  pool_index <- match(pool_ids, ids)
  ode_patient <- ids %in% pool_ids
  x <- initial$x[ids, , drop = FALSE]; eta <- initial$eta; psi <- initial$psi
  .sab_validate_population_inputs(x, eta, adapter$treatment, prior_mean, prior_sd)
  if (!adapter$psi_in_domain(psi)) stop("Invalid initial shared parameters.")
  categories <- c("initial", "local", "dynamic", "assignment", "common")
  ledger <- matrix(0, length(categories), 6L, dimnames = list(categories,
    c("attempts", "ode_integrations", "prediction_failures", "ode_failures",
      "elapsed_sec", "budget_units")))
  counters <- c(local_attempts = 0, local_accepts = 0, dynamic_attempts = 0,
    dynamic_accepts = 0, noise_attempts = 0, noise_accepts = 0,
    assignment_attempts = 0, assignment_moves = 0, assignment_self = 0,
    correction_rejects = 0, exact_cache_hits = 0, reserve_refreshes = 0)
  local_attempts <- local_accepts <- setNames(integer(n), ids)
  assignment_moves <- setNames(integer(m), pool_ids)
  note <- function(category, prediction, integrations, elapsed, units) {
    failed <- !isTRUE(prediction$ok)
    ode_failed <- failed && prediction$reason %in%
      c("ode_failure", "invalid_ode_output", "invalid_ode_times")
    ledger[category, ] <<- ledger[category, ] +
      c(1, integrations, failed, ode_failed, elapsed, units)
  }
  evaluate <- function(i, value, proposed_psi, category) {
    start <- proc.time()[["elapsed"]]
    prediction <- adapter$solve_prediction(ids[[i]], value, proposed_psi)
    likelihood <- adapter$loglik_from_prediction(prediction, proposed_psi)
    note(category, prediction, .sab_shared_original_calls(prediction),
         proc.time()[["elapsed"]] - start, as.integer(ode_patient[[i]]))
    list(loglik = likelihood, prediction = prediction)
  }
  active <- lapply(seq_len(n), function(i) evaluate(i, x[i, ], psi, "initial"))
  ll <- vapply(active, `[[`, numeric(1L), "loglik")
  if (any(!is.finite(ll))) stop("Initial sealed likelihood is not finite.")
  log_prior <- adapter$log_hyperprior(eta, psi)
  if (!is.finite(log_prior)) stop("Initial hyperprior is not finite.")
  draw_q <- function() {
    component <- reference_components[[sample.int(length(reference_components), 1L)]]
    component$mean + component$sd * stats::rnorm(8L)
  }
  pool <- new.env(parent = emptyenv())
  if (method == "recycling") {
    pool$assigned <- seq_len(m)
    pool$z <- rbind(x[pool_index, , drop = FALSE],
                   t(vapply(seq_len(K - m), function(j) draw_q(), numeric(8L))))
    pool$common <- vector("list", K)
    pool$exact <- vector("list", K * m)
    pool$log_q <- .sab_chain_log_q(pool$z, reference_components)
  }
  exact_index <- function(i, j) i + (j - 1L) * m
  replace_slot <- function(j, value) {
    pool$z[j, ] <- value
    pool$common[j] <- list(NULL)
    pool$exact[((j - 1L) * m + 1L):(j * m)] <- rep(list(NULL), m)
    pool$log_q[[j]] <- .sab_chain_log_q(pool$z[j, , drop = FALSE], reference_components)
  }
  populate_exact_active <- function() {
    for (p in seq_len(m)) {
      pool$exact[[exact_index(p, pool$assigned[[p]])]] <- active[[pool_index[[p]]]]$prediction
    }
  }
  beta <- rep(config$initial_pcn_beta, n)
  dynamic_sd <- config$initial_dynamic_sd
  noise_sd <- rep(config$initial_noise_sd, 2L)
  draw_names <- c(names(eta), names(psi), "treated_logit_location")
  draws <- matrix(NA_real_, config$max_sweeps, length(draw_names),
                  dimnames = list(NULL, draw_names))
  costs <- numeric(config$max_sweeps)
  sweep_number <- 0L
  # Leave enough budget for a complete worst-case next sweep. Dynamic calls,
  # local calls, all K common entries and m exact corrections are bounded.
  max_sweep_units <- 2L * m + if (method == "recycling") K + m else 0L
  snapshot <- function() {
    elapsed <- proc.time() - started
    list(method = method, seed = seed, config = config,
      draws = draws[seq_len(sweep_number), , drop = FALSE],
      warmup = min(config$warmup, sweep_number),
      costs = costs[seq_len(sweep_number)],
      ledger = data.frame(category = rownames(ledger), ledger, row.names = NULL),
      counters = counters, local_attempts = local_attempts, local_accepts = local_accepts,
      assignment_moves = assignment_moves, elapsed_sec = unname(elapsed[["elapsed"]]),
      cpu_sec = unname(sum(elapsed[c("user.self", "sys.self")])),
      tuning = list(pcn_beta = beta, dynamic_sd = dynamic_sd, noise_sd = noise_sd),
      final_state = list(x = x, eta = eta, psi = psi),
      patient_ids = ids, ode_patient_ids = pool_ids,
      target_fingerprint = adapter$target_fingerprint,
      proposal_fingerprint = probe$target_fingerprint)
  }
  repeat {
    if (sweep_number >= config$max_sweeps ||
        sum(ledger[, "budget_units"]) + max_sweep_units > config$budget) break
    sweep_number <- sweep_number + 1L
    warmup <- sweep_number <= config$warmup
    gain <- (sweep_number + 10)^(-.6)
    for (i in seq_len(n)) {
      mean <- adapter$population_mean(ids[[i]], eta)
      sd <- exp(eta[10:17])
      direct <- stats::runif(1L) < config$direct_probability
      proposed <- if (direct) mean + sd * stats::rnorm(8L) else
        mean + sqrt(1 - beta[[i]]^2) * (x[i, ] - mean) +
          beta[[i]] * sd * stats::rnorm(8L)
      names(proposed) <- colnames(x)
      candidate <- evaluate(i, proposed, psi, "local")
      accepted <- log(stats::runif(1L)) < min(0, candidate$loglik - ll[[i]])
      counters[["local_attempts"]] <- counters[["local_attempts"]] + 1L
      counters[["local_accepts"]] <- counters[["local_accepts"]] + accepted
      local_attempts[[i]] <- local_attempts[[i]] + 1L
      local_accepts[[i]] <- local_accepts[[i]] + accepted
      if (accepted) {
        x[i, ] <- proposed; active[[i]] <- candidate; ll[[i]] <- candidate$loglik
        p <- match(i, pool_index)
        if (method == "recycling" && !is.na(p)) replace_slot(pool$assigned[[p]], proposed)
      }
      if (warmup && !direct) beta[[i]] <- min(.95, max(.02,
        stats::plogis(stats::qlogis(beta[[i]]) + gain * (accepted - .3))))
    }
    eta <- sab_system_a_population_update(x, eta, adapter$treatment, prior_mean, prior_sd)
    if (!adapter$eta_in_domain(eta)) stop("Population update left representable domain.")
    log_prior <- adapter$log_hyperprior(eta, psi)
    if (!is.finite(log_prior)) stop("Updated population hyperprior is not representable.")
    for (j in 1:2) {
      proposed_psi <- psi
      proposed_psi[[j + 2L]] <- psi[[j + 2L]] + stats::rnorm(1L, sd = noise_sd[[j]])
      accepted <- FALSE
      if (adapter$psi_in_domain(proposed_psi)) {
        candidate_ll <- vapply(active, function(a) {
          adapter$loglik_from_prediction(a$prediction, proposed_psi)
        }, numeric(1L))
        candidate_prior <- adapter$log_hyperprior(eta, proposed_psi)
        accepted <- log(stats::runif(1L)) < min(0,
          sum(candidate_ll - ll) + candidate_prior - log_prior)
        if (accepted) { psi <- proposed_psi; ll <- candidate_ll; log_prior <- candidate_prior }
      }
      counters[["noise_attempts"]] <- counters[["noise_attempts"]] + 1L
      counters[["noise_accepts"]] <- counters[["noise_accepts"]] + accepted
      if (warmup) noise_sd[[j]] <- min(1, max(.001,
        noise_sd[[j]] * exp(gain * (accepted - .44))))
    }
    if (sweep_number %% config$dynamic_every == 0L) {
      proposed_psi <- psi; proposed_psi[1:2] <- psi[1:2] + stats::rnorm(2L, sd = dynamic_sd)
      accepted <- FALSE
      if (adapter$psi_in_domain(proposed_psi)) {
        candidate <- lapply(seq_len(n), function(i) evaluate(i, x[i, ], proposed_psi, "dynamic"))
        candidate_ll <- vapply(candidate, `[[`, numeric(1L), "loglik")
        candidate_prior <- adapter$log_hyperprior(eta, proposed_psi)
        accepted <- log(stats::runif(1L)) < min(0,
          sum(candidate_ll - ll) + candidate_prior - log_prior)
        if (accepted) {
          psi <- proposed_psi; active <- candidate; ll <- candidate_ll; log_prior <- candidate_prior
          if (method == "recycling") {
            pool$common <- vector("list", K); pool$exact <- vector("list", m * K)
          }
        }
      }
      counters[["dynamic_attempts"]] <- counters[["dynamic_attempts"]] + 1L
      counters[["dynamic_accepts"]] <- counters[["dynamic_accepts"]] + accepted
      if (warmup) dynamic_sd <- pmin(1, pmax(.001, dynamic_sd *
        exp((sweep_number / config$dynamic_every + 10)^(-.6) * (accepted - .234))))
    }
    # Noise changes alter likelihoods, but predictions remain valid. Always
    # synchronize active loglik before passing an endpoint to the correction.
    for (i in seq_len(n)) active[[i]]$loglik <- ll[[i]]
    if (method == "recycling") {
      populate_exact_active()
      unused <- sab_shared_pool_reserve_slots(pool$assigned, K)
      refreshed <- unused[sample.int(length(unused), config$reserve_refresh)]
      for (j in refreshed) replace_slot(j, draw_q())
      counters[["reserve_refreshes"]] <- counters[["reserve_refreshes"]] + length(refreshed)
      for (j in seq_len(K)) if (is.null(pool$common[[j]])) {
        start <- proc.time()[["elapsed"]]
        solved <- .sab_shared_solve(probe, pool$z[j, ], psi)
        pool$common[[j]] <- lapply(pool_ids, function(id) .sab_shared_prediction(probe, solved, id, psi))
        note("common", solved, solved$ode_integrations, proc.time()[["elapsed"]] - start, 1L)
      }
      approximate <- vapply(pool$common, function(predictions) {
        sab_shared_pool_positive_loglik(vapply(predictions, function(prediction) {
          adapter$loglik_from_prediction(prediction, psi)
        }, numeric(1L)), config$log_floor)
      }, numeric(m))
      dim(approximate) <- c(m, K)
      for (p in seq_len(m)) {
        i <- pool_index[[p]]
        pop <- .sab_chain_normal_logdensity(pool$z,
          adapter$population_mean(ids[[i]], eta), exp(eta[10:17]))
        result <- sab_shared_pool_corrected_assignment(pool$assigned, p,
          approximate[p, ], pop, pool$log_q, active[[i]], function(j) {
            key <- exact_index(p, j)
            prediction <- pool$exact[[key]]
            if (is.null(prediction)) {
              value <- evaluate(i, pool$z[j, ], psi, "assignment")
              pool$exact[[key]] <- value$prediction
              value
            } else {
              counters[["exact_cache_hits"]] <<- counters[["exact_cache_hits"]] + 1L
              list(prediction = prediction, loglik = adapter$loglik_from_prediction(prediction, psi))
            }
          })
        counters[["assignment_attempts"]] <- counters[["assignment_attempts"]] + 1L
        counters[["assignment_self"]] <- counters[["assignment_self"]] + (result$exact_callback_calls == 0L)
        counters[["correction_rejects"]] <- counters[["correction_rejects"]] + !result$accepted
        counters[["assignment_moves"]] <- counters[["assignment_moves"]] + result$moved
        if (result$moved) {
          pool$assigned <- result$assigned; x[i, ] <- pool$z[result$selected_slot, ]
          active[[i]] <- result$exact_current; ll[[i]] <- active[[i]]$loglik
          assignment_moves[[p]] <- assignment_moves[[p]] + 1L
        }
      }
    }
    draws[sweep_number, ] <- c(eta, psi, eta[[8L]] + eta[[9L]])
    costs[[sweep_number]] <- sum(ledger[, "budget_units"])
    if (!is.null(checkpoint_path) && sweep_number %% config$checkpoint_every == 0L) {
      saveRDS(snapshot(), checkpoint_path)
      cat(method, seed, "sweep", sweep_number, "forward budget", costs[[sweep_number]], "\n")
      flush.console()
    }
  }
  snapshot()
}
