# Exact correction of a surrogate Gibbs assignment in a live auxiliary pool.
#
# Target: p(c) prod_i [L_i(z[a_i];c) g_i(z[a_i];c)] prod_unused_j q(z_j),
# with injective a and FIXED normalized q. Other patients' assignments and c
# are held fixed during one call. The eligible set (current + unused slots) is
# unchanged by the forward/reverse assignment, so the surrogate Gibbs kernel
# is reversible with weights tilde_L_i(z_j) g_i(z_j)/q(z_j). Consequently its
# exact MH correction is [L_i(new)/tilde_L_i(new)]/[L_i(old)/tilde_L_i(old)].
#
# No Gaussian assumption, ODE solver, global update or adaptive surrogate is
# implemented here. The caller owns slot values and caches. Refresh an unused
# slot directly from q, without likelihood screening; invalidate caches whenever
# its state, shared dynamics or numerical-target fingerprint changes.

.sab_pool_kernel_log_vector <- function(x, argument, finite = FALSE) {
  if (!is.numeric(x) || !is.null(dim(x)) || !length(x) || anyNA(x) ||
      any(x == Inf) || (finite && any(!is.finite(x)))) {
    stop("`", argument, "` must be a numeric log-density vector containing ",
         if (finite) "finite values." else "finite values or -Inf.",
         call. = FALSE)
  }
  invisible(x)
}

# The same deterministic floor must be used at both endpoints. This preserves
# surrogate support even when a common-grid solve fails but the sealed one does
# not. A floor is a proposal device only and never replaces an exact likelihood.
sab_shared_pool_positive_loglik <- function(raw_loglik, log_floor) {
  .sab_pool_kernel_log_vector(raw_loglik, "raw_loglik")
  if (!is.numeric(log_floor) || length(log_floor) != 1L || !is.finite(log_floor)) {
    stop("`log_floor` must be a finite scalar fixed before sampling.", call. = FALSE)
  }
  pmax(raw_loglik, log_floor)
}

sab_shared_pool_reserve_slots <- function(assigned, pool_size) {
  if (!is.numeric(pool_size) || length(pool_size) != 1L ||
      !is.finite(pool_size) || pool_size < 1 || pool_size != floor(pool_size) ||
      pool_size > .Machine$integer.max || !is.numeric(assigned) ||
      !is.null(dim(assigned)) || !length(assigned) || anyNA(assigned) ||
      any(!is.finite(assigned)) || any(assigned != floor(assigned)) ||
      any(assigned < 1 | assigned > pool_size) || anyDuplicated(assigned)) {
    stop("Assignments must be unique valid integer slot indices.", call. = FALSE)
  }
  setdiff(seq_len(pool_size), as.integer(assigned))
}

sab_shared_pool_assignment_probabilities <- function(
    assigned, patient_index, log_approx_likelihood, log_population, log_q) {
  .sab_pool_kernel_log_vector(log_q, "log_q", finite = TRUE)
  .sab_pool_kernel_log_vector(log_approx_likelihood, "log_approx_likelihood",
                              finite = TRUE)
  .sab_pool_kernel_log_vector(log_population, "log_population")
  pool_size <- length(log_q)
  if (length(log_approx_likelihood) != pool_size ||
      length(log_population) != pool_size) {
    stop("All slot log-density vectors must have the same length.", call. = FALSE)
  }
  reserves <- sab_shared_pool_reserve_slots(assigned, pool_size)
  if (!is.numeric(patient_index) || length(patient_index) != 1L ||
      !is.finite(patient_index) || patient_index != floor(patient_index) ||
      patient_index < 1 || patient_index > length(assigned)) {
    stop("`patient_index` must identify an assigned patient.", call. = FALSE)
  }
  current <- assigned[[patient_index]]
  if (!is.finite(log_population[[current]])) {
    stop("Current population density must be positive and finite.", call. = FALSE)
  }
  slots <- c(as.integer(current), reserves)
  log_weights <- log_approx_likelihood[slots] + log_population[slots] - log_q[slots]
  .sab_pool_kernel_log_vector(log_weights, "eligible surrogate log weights")
  if (!is.finite(log_weights[[1L]])) {
    stop("Current surrogate weight overflowed.", call. = FALSE)
  }
  weights <- exp(log_weights - max(log_weights))
  list(slots = slots, probabilities = unname(weights / sum(weights)),
       log_weights = unname(log_weights), current_slot = as.integer(current))
}

.sab_pool_kernel_exact_cache <- function(cache, argument, finite) {
  if (!is.list(cache) || !is.numeric(cache$loglik) || length(cache$loglik) != 1L) {
    stop("`", argument, "` must be a list with scalar numeric loglik.", call. = FALSE)
  }
  .sab_pool_kernel_log_vector(cache$loglik, argument, finite = finite)
  invisible(cache)
}

# Returns a proposed assignment decision without mutating caller state. On
# rejection, assigned and exact_current are returned unchanged. evaluated_exact
# is a separate optional memoization result for the candidate, not active state.
# evaluate_exact(slot) must return list(loglik=..., ...) using the sealed target;
# ordinary impossible states have -Inf, while unexpected callback errors abort.
sab_shared_pool_corrected_assignment <- function(
    assigned, patient_index, log_approx_likelihood, log_population, log_q,
    exact_current, evaluate_exact,
    selection_uniform = stats::runif(1L), acceptance_uniform = stats::runif(1L)) {
  proposal <- sab_shared_pool_assignment_probabilities(
    assigned, patient_index, log_approx_likelihood, log_population, log_q)
  .sab_pool_kernel_exact_cache(exact_current, "exact_current", finite = TRUE)
  if (!is.function(evaluate_exact)) stop("`evaluate_exact` must be a function.")
  for (u in list(selection_uniform, acceptance_uniform)) {
    if (!is.numeric(u) || length(u) != 1L || !is.finite(u) || u < 0 || u >= 1) {
      stop("Selection and acceptance uniforms must lie in [0,1).", call. = FALSE)
    }
  }
  positive <- which(proposal$probabilities > 0)
  cumulative <- cumsum(proposal$probabilities[positive])
  cumulative[[length(cumulative)]] <- 1
  position <- positive[[which(selection_uniform < cumulative)[[1L]]]]
  selected <- proposal$slots[[position]]
  current <- proposal$current_slot
  self <- selected == current
  evaluated <- NULL
  log_ratio <- 0
  if (!self) {
    evaluated <- evaluate_exact(selected)
    .sab_pool_kernel_exact_cache(evaluated, "evaluate_exact result", finite = FALSE)
    log_ratio <- (evaluated$loglik - log_approx_likelihood[[selected]]) -
      (exact_current$loglik - log_approx_likelihood[[current]])
    if (is.na(log_ratio)) stop("Exact correction produced an undefined log ratio.")
  }
  probability <- exp(min(0, log_ratio))
  accepted <- self || log(acceptance_uniform) < min(0, log_ratio)
  retained_assignment <- assigned
  retained_exact <- exact_current
  if (accepted && !self) {
    retained_assignment[[patient_index]] <- selected
    retained_exact <- evaluated
  }
  list(
    assigned = retained_assignment, exact_current = retained_exact,
    accepted = accepted, moved = accepted && !self, selected_slot = selected,
    previous_slot = current, selection_probability = proposal$probabilities[[position]],
    log_acceptance_ratio = log_ratio, acceptance_probability = probability,
    exact_callback_calls = as.integer(!self), evaluated_exact = evaluated,
    proposal = proposal
  )
}
