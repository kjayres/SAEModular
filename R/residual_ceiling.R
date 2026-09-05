# One-sided falsification only: no estimated posterior and no lower certificate.
# Freeze a normalized component q and independent probes giving c_probe >= c
# for every pointwise-valid c*q <= F. For independent IID draws from ANY frozen
# normalized h, supply log_ratio = log(F/h). With Z=min(F/(c_probe*h), cap),
# E_h[Z] <= M/c_probe, even if h has only partial support. Hoeffding therefore
# gives c/M <= min(1, 1/lower) with error probability alpha_per_cap. Account for
# all caps, contexts, candidate components and interim looks outside this helper.
# A high bound proves neither coverage, a useful certificate, nor good mixing.
sab_residual_ceiling <- function(log_ratio, log_c_probe, caps = c(2, 10, 50),
                                 alpha_per_cap = 0.05) {
  if (!is.numeric(log_c_probe) || length(log_c_probe) != 1L ||
      is.na(log_c_probe) || log_c_probe == Inf) {
    stop("`log_c_probe` must be finite or -Inf.", call. = FALSE)
  }
  if (!is.numeric(log_ratio) || !is.null(dim(log_ratio)) || anyNA(log_ratio) ||
      any(log_ratio == Inf) || (!length(log_ratio) && is.finite(log_c_probe))) {
    stop("`log_ratio` must contain finite values or -Inf, with IID draws.",
         call. = FALSE)
  }
  if (!is.numeric(caps) || !is.null(dim(caps)) || !length(caps) ||
      any(!is.finite(caps)) || any(caps <= 0) ||
      !is.numeric(alpha_per_cap) || length(alpha_per_cap) != 1L ||
      !is.finite(alpha_per_cap) || alpha_per_cap <= 0 || alpha_per_cap >= 1) {
    stop("Positive finite `caps` and scalar `alpha_per_cap` in (0, 1) required.",
         call. = FALSE)
  }
  if (log_c_probe == -Inf) {
    # An exact target zero at a point where q>0 forces any pointwise c to zero.
    # A numerically underflowed positive density must NOT be treated as a zero.
    return(data.frame(cap = caps, mean_capped_ratio = NA_real_,
                      mass_ratio_lower = Inf, collapse_upper = 0))
  }
  means <- vapply(caps, function(cap) {
    mean(exp(pmin(log_ratio - log_c_probe, log(cap))))
  }, numeric(1))
  lower <- pmax(0, means - caps * sqrt(-log(alpha_per_cap) /
                                      (2 * length(log_ratio))))
  data.frame(cap = caps, mean_capped_ratio = means, mass_ratio_lower = lower,
             collapse_upper = pmin(1, 1 / lower))
}

.sab_residual_gaussian_check <- function(mean, chol, mass) {
  d <- length(mean)
  if (!is.numeric(mean) || !is.null(dim(mean)) || d < 1L ||
      any(!is.finite(mean)) || !is.matrix(chol) || !is.numeric(chol) ||
      !identical(dim(chol), c(d, d)) || any(!is.finite(chol)) ||
      any(chol[upper.tri(chol)] != 0) || any(diag(chol) <= 0) ||
      !is.numeric(mass) || length(mass) != 1L || !is.finite(mass) ||
      mass <= 0 || mass >= 1) {
    stop("Finite mean, nonsingular lower Cholesky factor and mass in (0, 1) required.",
         call. = FALSE)
  }
  invisible(d)
}

# chol is LOWER triangular: covariance = chol %*% t(chol).
# Truncate to the Gaussian Mahalanobis ball containing probability `mass`.
sab_truncated_gaussian_draw <- function(n, mean, chol, mass) {
  d <- .sab_residual_gaussian_check(mean, chol, mass)
  if (!is.numeric(n) || length(n) != 1L || !is.finite(n) || n < 1 ||
      n != floor(n) || n > .Machine$integer.max) {
    stop("`n` must be a positive integer.", call. = FALSE)
  }
  direction <- matrix(stats::rnorm(n * d), nrow = n, ncol = d)
  norms <- sqrt(rowSums(direction^2))
  if (any(norms == 0)) stop("Zero Gaussian direction encountered.", call. = FALSE)
  radius <- sqrt(stats::qchisq(stats::runif(n) * mass, df = d))
  z <- direction * (radius / norms)
  draws <- sweep(z %*% t(chol), 2L, mean, "+")
  if (any(!is.finite(draws))) stop("Gaussian draws overflowed.", call. = FALSE)
  draws
}

sab_truncated_gaussian_logdensity <- function(x, mean, chol, mass) {
  d <- .sab_residual_gaussian_check(mean, chol, mass)
  if (is.numeric(x) && is.null(dim(x)) && length(x) == d) x <- matrix(x, nrow = 1L)
  if (!is.matrix(x) || !is.numeric(x) || ncol(x) != d || nrow(x) < 1L ||
      any(!is.finite(x))) {
    stop("`x` must have finite rows in the Gaussian coordinate order.", call. = FALSE)
  }
  z <- t(forwardsolve(chol, t(sweep(x, 2L, mean, "-"))))
  radius_sq <- rowSums(z^2)
  log_q <- -0.5 * (d * log(2 * pi) + radius_sq) -
    sum(log(diag(chol))) - log(mass)
  log_q[radius_sq > stats::qchisq(mass, df = d)] <- -Inf
  unname(log_q)
}
