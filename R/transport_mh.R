#' Transport all patient states for a proposed global move
#'
#' Each patient map preserves its own standardised coordinate.  Named inputs
#' are accepted only when map and state names are complete and identical; this
#' prevents silent patient reordering.
#'
#' @param maps A `sab_affine_map` or a list of patient-specific maps.
#' @param local_states A list of patient state vectors in the same order as
#'   `maps`.
#' @param from_global Current global parameter value.
#' @param to_global Proposed global parameter value.
#'
#' @return A list containing transported `local_states`, standardised
#'   coordinates, patient-level log Jacobians, and their sum.
#' @export
sab_transport_local_states <- function(maps,
                                       local_states,
                                       from_global,
                                       to_global) {
  if (inherits(maps, "sab_affine_map")) {
    maps <- list(maps)
  }
  if (!is.list(maps) || length(maps) < 1L ||
      !all(vapply(maps, inherits, logical(1), what = "sab_affine_map"))) {
    stop("`maps` must be a non-empty list of `sab_affine_map` objects.",
         call. = FALSE)
  }
  if (!is.list(local_states) || length(local_states) != length(maps)) {
    stop("`local_states` must be a list with one state per map.",
         call. = FALSE)
  }
  .sab_validate_patient_names(maps, local_states)

  n_patients <- length(maps)
  proposed <- vector("list", n_patients)
  standardised <- vector("list", n_patients)
  log_jacobian <- numeric(n_patients)

  for (i in seq_len(n_patients)) {
    from <- sab_affine_map_components(maps[[i]], from_global)
    to <- sab_affine_map_components(maps[[i]], to_global)
    if (length(from$mean) != length(to$mean)) {
      stop("Patient map ", i, " changed dimension between global states.",
           call. = FALSE)
    }
    if (!identical(names(from$mean), names(to$mean))) {
      stop("Patient map ", i,
           " changed coordinate names between global states.",
           call. = FALSE)
    }
    x <- .sab_validate_affine_vector(
      local_states[[i]],
      expected_length = length(from$mean),
      argument = paste0("local_states[[", i, "]]")
    )
    .sab_match_coordinate_names(
      x, from$mean, paste0("local_states[[", i, "]]" )
    )
    u <- as.numeric(forwardsolve(from$chol, x - from$mean))
    names(u) <- if (!is.null(names(x))) names(x) else names(from$mean)
    standardised[[i]] <- u
    proposed[[i]] <- as.numeric(to$mean + to$chol %*% u)
    names(proposed[[i]]) <- if (!is.null(names(x))) names(x) else names(to$mean)
    log_jacobian[[i]] <- to$log_abs_det - from$log_abs_det
  }

  names(proposed) <- names(local_states)
  names(standardised) <- names(local_states)
  names(log_jacobian) <- names(local_states)

  structure(
    list(
      local_states = proposed,
      standardised = standardised,
      patient_log_jacobian = log_jacobian,
      log_abs_det_jacobian = sum(log_jacobian)
    ),
    class = "sab_transport_result"
  )
}

#' Build an exact deterministic-transport Metropolis proposal
#'
#' This implements a transformation-based Metropolis-Hastings proposal.  A
#' global proposal `global -> global_prime` deterministically induces all local
#' proposals.  If `pi(global, locals)` is the exact joint density with respect
#' to the original local coordinates, the returned ratio is
#'
#' `log pi(proposed) - log pi(current) + log q(current | proposed) -
#' log q(proposed | current) + log |det d locals_proposed / d locals|`.
#'
#' A target evaluator may return either a numeric log density or a list with a
#' numeric `log_density` element and arbitrary additional fields, such as a
#' proposed prediction cache.  Keeping that complete evaluation in the proposal
#' object permits atomic cache replacement after acceptance.
#'
#' This move preserves every standardised local coordinate.  It is therefore
#' not irreducible by itself and must be composed with exact local-coordinate
#' refresh kernels in a retained chain.
#'
#' @param current_global Current global parameter value.
#' @param proposed_global Proposed global parameter value drawn from the global
#'   proposal.
#' @param current_locals List of current patient state vectors.
#' @param maps Patient-specific affine maps; see [sab_transport_local_states()].
#' @param log_target Function `log_target(global, locals)` evaluating the exact
#'   joint target.  It must have no side effects on the current chain state.
#' @param log_global_proposal Function `log_global_proposal(to, from)` returning
#'   the global proposal log density.  Use `function(to, from) 0` for a
#'   symmetric proposal.
#' @param current_evaluation Optional cached current target evaluation, in the
#'   same numeric-or-list form accepted from `log_target`.
#'
#' @return An object of class `sab_transport_mh_proposal` containing both
#'   states, both target evaluations, the Jacobian, and the exact log ratio.
#' @export
sab_transport_mh_proposal <- function(current_global,
                                      proposed_global,
                                      current_locals,
                                      maps,
                                      log_target,
                                      log_global_proposal,
                                      current_evaluation = NULL) {
  if (!is.function(log_target)) {
    stop("`log_target` must be a function.", call. = FALSE)
  }
  if (!is.function(log_global_proposal)) {
    stop("`log_global_proposal` must be a function.", call. = FALSE)
  }

  transport <- sab_transport_local_states(
    maps = maps,
    local_states = current_locals,
    from_global = current_global,
    to_global = proposed_global
  )

  if (is.null(current_evaluation)) {
    current_evaluation <- log_target(current_global, current_locals)
  }
  current_evaluation <- .sab_normalise_target_evaluation(
    current_evaluation,
    label = "current target"
  )
  if (!is.finite(current_evaluation$log_density)) {
    stop("The current state must have finite exact log density.",
         call. = FALSE)
  }

  proposed_evaluation <- .sab_normalise_target_evaluation(
    log_target(proposed_global, transport$local_states),
    label = "proposed target"
  )

  log_q_forward <- .sab_validate_log_value(
    log_global_proposal(proposed_global, current_global),
    label = "forward global proposal",
    allow_negative_infinity = FALSE
  )
  log_q_reverse <- .sab_validate_log_value(
    log_global_proposal(current_global, proposed_global),
    label = "reverse global proposal",
    allow_negative_infinity = TRUE
  )

  if (is.infinite(proposed_evaluation$log_density) ||
      is.infinite(log_q_reverse)) {
    log_ratio <- -Inf
  } else {
    log_ratio <- proposed_evaluation$log_density -
      current_evaluation$log_density +
      log_q_reverse - log_q_forward +
      transport$log_abs_det_jacobian
  }

  structure(
    list(
      current_state = list(
        global = current_global,
        locals = current_locals
      ),
      proposed_state = list(
        global = proposed_global,
        locals = transport$local_states
      ),
      current_evaluation = current_evaluation,
      proposed_evaluation = proposed_evaluation,
      standardised = transport$standardised,
      patient_log_jacobian = transport$patient_log_jacobian,
      log_abs_det_jacobian = transport$log_abs_det_jacobian,
      log_q_forward = log_q_forward,
      log_q_reverse = log_q_reverse,
      log_ratio = log_ratio,
      log_acceptance_probability = min(0, log_ratio)
    ),
    class = "sab_transport_mh_proposal"
  )
}

#' Accept or reject a deterministic-transport Metropolis proposal
#'
#' Supplying `log_uniform` makes the decision reproducible and allows exact
#' boundary tests.  On rejection this function returns the complete current
#' target evaluation; it never exposes the proposed cache as the selected one.
#'
#' @param proposal A `sab_transport_mh_proposal`.
#' @param log_uniform Logarithm of a draw from `Uniform(0, 1)`.  By default a
#'   fresh draw is generated.
#'
#' @return A list with `accepted`, selected `state`, selected `evaluation`, and
#'   the proposal diagnostics.
#' @export
sab_decide_transport_mh <- function(proposal,
                                    log_uniform = log(stats::runif(1L))) {
  if (!inherits(proposal, "sab_transport_mh_proposal")) {
    stop("`proposal` must be a `sab_transport_mh_proposal`.",
         call. = FALSE)
  }
  if (!is.numeric(log_uniform) || length(log_uniform) != 1L ||
      is.na(log_uniform) || is.nan(log_uniform) || log_uniform > 0) {
    stop("`log_uniform` must be a non-positive scalar.", call. = FALSE)
  }

  accepted <- log_uniform < proposal$log_acceptance_probability
  if (accepted) {
    state <- proposal$proposed_state
    evaluation <- proposal$proposed_evaluation
  } else {
    state <- proposal$current_state
    evaluation <- proposal$current_evaluation
  }

  structure(
    list(
      accepted = accepted,
      state = state,
      evaluation = evaluation,
      log_uniform = log_uniform,
      log_ratio = proposal$log_ratio,
      log_acceptance_probability = proposal$log_acceptance_probability,
      patient_log_jacobian = proposal$patient_log_jacobian,
      log_abs_det_jacobian = proposal$log_abs_det_jacobian
    ),
    class = "sab_transport_mh_decision"
  )
}

.sab_validate_patient_names <- function(maps, local_states) {
  map_names <- names(maps)
  state_names <- names(local_states)
  maps_are_named <- !is.null(map_names) && all(nzchar(map_names))
  states_are_named <- !is.null(state_names) && all(nzchar(state_names))

  if (maps_are_named || states_are_named) {
    if (!maps_are_named || !states_are_named ||
        anyDuplicated(map_names) || anyDuplicated(state_names) ||
        !identical(map_names, state_names)) {
      stop(
        "Named `maps` and `local_states` must have identical, unique names ",
        "in identical order.",
        call. = FALSE
      )
    }
  } else if ((!is.null(map_names) && any(nzchar(map_names))) ||
             (!is.null(state_names) && any(nzchar(state_names)))) {
    stop("Patient names must be either complete or absent.", call. = FALSE)
  }
  invisible(TRUE)
}

.sab_normalise_target_evaluation <- function(value, label) {
  if (is.list(value)) {
    if (is.null(value$log_density)) {
      stop("The ", label, " evaluation lacks `log_density`.",
           call. = FALSE)
    }
    value$log_density <- .sab_validate_log_value(
      value$log_density,
      label = label,
      allow_negative_infinity = TRUE
    )
    return(value)
  }

  list(
    log_density = .sab_validate_log_value(
      value,
      label = label,
      allow_negative_infinity = TRUE
    )
  )
}

.sab_validate_log_value <- function(value,
                                    label,
                                    allow_negative_infinity) {
  valid <- is.numeric(value) && length(value) == 1L &&
    !is.na(value) && !is.nan(value) && value != Inf
  if (!allow_negative_infinity) {
    valid <- valid && is.finite(value)
  }
  if (!valid) {
    permitted <- if (allow_negative_infinity) "finite or -Inf" else "finite"
    stop("The ", label, " log density must be ", permitted, ".",
         call. = FALSE)
  }
  as.numeric(value)
}
