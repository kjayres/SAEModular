# Proposal-only numerical guide. The sealed adapter and its caches are never
# modified. This object contains no population density or patient posterior fit.
sab_make_system_a_coarse_guide <- function(adapter, relative_tolerance,
                                          absolute_tolerance, log_floor = -1e6) {
  sab_validate_system_a_adapter(adapter)
  controls <- c(relative_tolerance, absolute_tolerance, log_floor)
  if (length(controls) != 3L || any(!is.finite(controls)) ||
      relative_tolerance <= 0 || absolute_tolerance <= 0) stop("Invalid guide controls.")
  owner <- environment(adapter$solve_prediction)
  upstream <- get("upstream", owner, inherits = FALSE)
  original <- get("exact_likelihood", owner, inherits = FALSE)
  patients <- original$patients[adapter$patient_ids]
  for (id in names(patients)) {
    patients[[id]]$controls$rel_tol <- relative_tolerance
    patients[[id]]$controls$abs_tol <- absolute_tolerance
  }
  signature <- paste("proposal_only_vode_bdf", format(relative_tolerance, digits = 17),
                     format(absolute_tolerance, digits = 17), sep = "/")
  evaluate <- function(patient_id, x, psi) {
    if (length(patient_id) != 1L || !patient_id %in% names(patients)) stop("Unknown patient.")
    prediction <- upstream$sysa_solve_prediction(patients[[patient_id]], x, psi,
                                                solver = "vode_bdf")
    raw <- upstream$sysa_loglik_from_prediction(prediction, psi)
    # Failure of the approximate solver must not remove exact target support.
    # A fixed likelihood floor is only a proposal device, never an exact value.
    list(loglik = max(raw, log_floor), raw_loglik = raw,
         ok = isTRUE(prediction$ok),
         reason = if (isTRUE(prediction$ok)) "ok" else prediction$reason,
         ode_integrations = .sab_shared_original_calls(prediction),
         source = "approximate", signature = signature)
  }
  list(evaluate = evaluate, signature = signature, log_floor = log_floor,
       relative_tolerance = relative_tolerance, absolute_tolerance = absolute_tolerance,
       reference_target_fingerprint = adapter$target_fingerprint)
}

# Fixed-context Metropolis macro proposal: k repetitions of the SAME reversible
# cheap MH kernel, then an exact endpoint correction. Targets include the actual
# population density; no Gaussian population assumption enters this identity.
# propose(x) returns x and log_reverse_minus_forward. Adaptation/changes of
# globals or guide require a new, consistent current cache outside this call.
sab_surrogate_macro_step <- function(current, propose, evaluate_cheap,
                                     evaluate_exact, steps,
                                     uniform = function() stats::runif(1L)) {
  scalar <- function(v, finite = FALSE) {
    is.numeric(v) && length(v) == 1L && !is.na(v) && v != Inf && (!finite || is.finite(v))
  }
  if (!is.list(current) || !scalar(current$exact$log_target, TRUE) ||
      !scalar(current$cheap$log_target, TRUE) || !is.function(propose) ||
      !is.function(evaluate_cheap) || !is.function(evaluate_exact) ||
      length(steps) != 1L || !is.finite(steps) || steps < 1 || steps != floor(steps)) {
    stop("Invalid macro-kernel state or controls.")
  }
  log_uniform <- function() {
    u <- uniform()
    if (!scalar(u, TRUE) || u < 0 || u >= 1) stop("Uniform outside [0,1).")
    log(u)
  }
  x <- current$x; cheap <- current$cheap; inner_accepts <- 0L
  for (j in seq_len(steps)) {
    candidate <- propose(x)
    if (!is.list(candidate) || !scalar(candidate$log_reverse_minus_forward, TRUE)) {
      stop("Malformed proposal.")
    }
    proposed_cheap <- evaluate_cheap(candidate$x)
    if (!is.list(proposed_cheap) || !scalar(proposed_cheap$log_target)) stop("Invalid cheap density.")
    ratio <- proposed_cheap$log_target - cheap$log_target + candidate$log_reverse_minus_forward
    if (log_uniform() < min(0, ratio)) {
      x <- candidate$x; cheap <- proposed_cheap; inner_accepts <- inner_accepts + 1L
    }
  }
  if (identical(x, current$x)) return(list(current = current, accepted = TRUE,
    moved = FALSE, exact_calls = 0L, cheap_calls = steps, inner_accepts = inner_accepts,
    log_correction = 0))
  exact <- evaluate_exact(x)
  if (!is.list(exact) || !scalar(exact$log_target)) stop("Invalid exact density.")
  ratio <- (exact$log_target - cheap$log_target) -
    (current$exact$log_target - current$cheap$log_target)
  accepted <- log_uniform() < min(0, ratio)
  retained <- if (accepted) list(x = x, exact = exact, cheap = cheap) else current
  list(current = retained, accepted = accepted, moved = accepted,
       exact_calls = 1L, cheap_calls = steps, inner_accepts = inner_accepts,
       log_correction = ratio)
}
