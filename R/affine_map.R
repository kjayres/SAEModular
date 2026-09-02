#' Construct a frozen conditional affine map
#'
#' The map represents a patient state as
#' `x = mean(global) + chol(global) %*% u`.  The supplied functions must be
#' deterministic during retained MCMC.  This constructor deliberately accepts
#' functions rather than fitted-model objects so that proposal fitting is kept
#' separate from the exact transition kernel.
#'
#' @param mean_fn Function of `global` returning a finite numeric vector.
#' @param chol_fn Function of `global` returning a finite lower-triangular
#'   matrix with strictly positive diagonal.
#' @param name Short identifier used in validation errors.
#' @param check_at Optional global state at which to validate both functions.
#' @param triangular_tolerance Relative tolerance used when checking that the
#'   upper triangle is zero.
#'
#' @return An object of class `sab_affine_map`.
#' @export
sab_new_affine_map <- function(mean_fn,
                               chol_fn,
                               name = "affine_map",
                               check_at = NULL,
                               triangular_tolerance = sqrt(.Machine$double.eps)) {
  if (!is.function(mean_fn)) {
    stop("`mean_fn` must be a function.", call. = FALSE)
  }
  if (!is.function(chol_fn)) {
    stop("`chol_fn` must be a function.", call. = FALSE)
  }
  if (!is.character(name) || length(name) != 1L || is.na(name) || !nzchar(name)) {
    stop("`name` must be one non-empty string.", call. = FALSE)
  }
  if (!is.numeric(triangular_tolerance) ||
      length(triangular_tolerance) != 1L ||
      !is.finite(triangular_tolerance) ||
      triangular_tolerance < 0) {
    stop("`triangular_tolerance` must be one finite non-negative number.",
         call. = FALSE)
  }

  map <- structure(
    list(
      mean_fn = mean_fn,
      chol_fn = chol_fn,
      name = name,
      triangular_tolerance = triangular_tolerance
    ),
    class = "sab_affine_map"
  )

  if (!is.null(check_at)) {
    sab_affine_map_components(map, check_at)
  }
  map
}

#' Evaluate and validate an affine map
#'
#' @param map A `sab_affine_map`.
#' @param global Global parameter value passed to the map functions.
#'
#' @return A list containing `mean`, `chol`, and `log_abs_det`.
#' @export
sab_affine_map_components <- function(map, global) {
  .sab_assert_affine_map(map)

  mean <- map$mean_fn(global)
  chol <- map$chol_fn(global)
  label <- paste0("Affine map `", map$name, "`")

  if (!is.numeric(mean) || !is.null(dim(mean)) || length(mean) < 1L ||
      anyNA(mean) || any(!is.finite(mean))) {
    stop(label, " returned an invalid mean vector.", call. = FALSE)
  }
  if (!is.matrix(chol) || !is.numeric(chol) ||
      nrow(chol) != length(mean) || ncol(chol) != length(mean) ||
      anyNA(chol) || any(!is.finite(chol))) {
    stop(label, " returned an invalid Cholesky matrix.", call. = FALSE)
  }

  diagonal <- diag(chol)
  if (any(diagonal <= 0)) {
    stop(label, " must have a strictly positive Cholesky diagonal.",
         call. = FALSE)
  }

  upper <- chol[upper.tri(chol)]
  matrix_scale <- max(1, max(abs(chol)))
  if (length(upper) > 0L &&
      any(abs(upper) > map$triangular_tolerance * matrix_scale)) {
    stop(label, " must be lower triangular.", call. = FALSE)
  }

  # Remove tolerated numerical noise so that forward substitution and the
  # reported Jacobian refer to exactly the same transformation.
  chol[upper.tri(chol)] <- 0

  mean_names <- names(mean)
  mean <- as.numeric(mean)
  names(mean) <- mean_names
  .sab_validate_coordinate_names(mean_names, label)
  list(
    mean = mean,
    chol = chol,
    log_abs_det = sum(log(diagonal))
  )
}

#' Convert a patient state to its standardised coordinate
#'
#' @param map A `sab_affine_map`.
#' @param x Patient state in the original coordinates.
#' @param global Global parameter value.
#'
#' @return The standardised patient coordinate `u`.
#' @export
sab_affine_standardise <- function(map, x, global) {
  components <- sab_affine_map_components(map, global)
  x <- .sab_validate_affine_vector(x, length(components$mean), "x")
  .sab_match_coordinate_names(x, components$mean, "x")
  value <- as.numeric(forwardsolve(components$chol, x - components$mean))
  names(value) <- if (!is.null(names(x))) names(x) else names(components$mean)
  value
}

#' Convert a standardised coordinate to a patient state
#'
#' @param map A `sab_affine_map`.
#' @param u Standardised patient coordinate.
#' @param global Global parameter value.
#'
#' @return The patient state `x`.
#' @export
sab_affine_unstandardise <- function(map, u, global) {
  components <- sab_affine_map_components(map, global)
  u <- .sab_validate_affine_vector(u, length(components$mean), "u")
  .sab_match_coordinate_names(u, components$mean, "u")
  value <- as.numeric(components$mean + components$chol %*% u)
  names(value) <- if (!is.null(names(components$mean))) {
    names(components$mean)
  } else {
    names(u)
  }
  value
}

#' Transport a patient state between two global parameter values
#'
#' The standardised coordinate is held fixed.  Thus the transformation is
#' `x_to = mean(to) + chol(to) %*% solve(chol(from), x_from - mean(from))`.
#'
#' @param map A `sab_affine_map`.
#' @param x Patient state at `from_global`.
#' @param from_global Current global parameter value.
#' @param to_global Proposed global parameter value.
#'
#' @return The transported patient state.
#' @export
sab_affine_transport <- function(map, x, from_global, to_global) {
  from <- sab_affine_map_components(map, from_global)
  to <- sab_affine_map_components(map, to_global)
  if (length(from$mean) != length(to$mean)) {
    stop("An affine map changed dimension between global states.",
         call. = FALSE)
  }
  if (!identical(names(from$mean), names(to$mean))) {
    stop("An affine map changed coordinate names between global states.",
         call. = FALSE)
  }
  x <- .sab_validate_affine_vector(x, length(from$mean), "x")
  .sab_match_coordinate_names(x, from$mean, "x")
  u <- as.numeric(forwardsolve(from$chol, x - from$mean))
  value <- as.numeric(to$mean + to$chol %*% u)
  names(value) <- if (!is.null(names(x))) names(x) else names(to$mean)
  value
}

#' Invert a previously defined affine transport
#'
#' @param map A `sab_affine_map`.
#' @param x_to Patient state at `to_global`.
#' @param from_global Original global parameter value.
#' @param to_global Destination global parameter value.
#'
#' @return The corresponding state at `from_global`.
#' @export
sab_affine_inverse_transport <- function(map, x_to, from_global, to_global) {
  sab_affine_transport(
    map = map,
    x = x_to,
    from_global = to_global,
    to_global = from_global
  )
}

#' Log absolute determinant of a map's scale matrix
#'
#' @param map A `sab_affine_map`.
#' @param global Global parameter value.
#'
#' @return `log(abs(det(chol(global))))`.
#' @export
sab_affine_log_abs_det <- function(map, global) {
  sab_affine_map_components(map, global)$log_abs_det
}

#' Log Jacobian of an affine transport
#'
#' @param map A `sab_affine_map`.
#' @param from_global Current global parameter value.
#' @param to_global Proposed global parameter value.
#'
#' @return The log absolute Jacobian determinant for `x_from -> x_to`.
#' @export
sab_affine_transport_log_jacobian <- function(map, from_global, to_global) {
  sab_affine_log_abs_det(map, to_global) -
    sab_affine_log_abs_det(map, from_global)
}

.sab_assert_affine_map <- function(map) {
  if (!inherits(map, "sab_affine_map")) {
    stop("`map` must be a `sab_affine_map`.", call. = FALSE)
  }
  invisible(TRUE)
}

.sab_validate_affine_vector <- function(value, expected_length, argument) {
  if (!is.numeric(value) || !is.null(dim(value)) ||
      length(value) != expected_length || anyNA(value) ||
      any(!is.finite(value))) {
    stop(
      "`", argument, "` must be a finite numeric vector of length ",
      expected_length, ".",
      call. = FALSE
    )
  }
  value_names <- names(value)
  value <- as.numeric(value)
  names(value) <- value_names
  value
}

.sab_validate_coordinate_names <- function(value_names, label) {
  if (!is.null(value_names) &&
      (length(value_names) < 1L || anyNA(value_names) ||
       any(!nzchar(value_names)) || anyDuplicated(value_names))) {
    stop(label, " returned invalid coordinate names.", call. = FALSE)
  }
  invisible(TRUE)
}

.sab_match_coordinate_names <- function(value, reference, argument) {
  value_names <- names(value)
  reference_names <- names(reference)
  if (!is.null(value_names) && !is.null(reference_names) &&
      !identical(value_names, reference_names)) {
    stop("Named `", argument,
         "` coordinates must match the affine-map coordinates in order.",
         call. = FALSE)
  }
  invisible(TRUE)
}
