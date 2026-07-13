# Phase-space embedding utilities
#
# Family: embed
# Shared reconstruction utilities consumed by the lyapunov and rqa families
# (and optionally fractal, for box-counting on reconstructed attractors).
# Kept as its own family rather than duplicated internal helpers, since
# multiple downstream families depend on the same reconstruction step.
#
# No reference implementation was specified for this family (no preference
# given) -- all three functions below are clean-room implementations of the
# standard published methods. See inst/COPYRIGHTS.

#' Time-Delay (Takens) Embedding
#'
#' Reconstructs a phase-space trajectory from a scalar time series using
#' Takens' (1981) delay-coordinate method.
#'
#' @param x Numeric vector. The time series to embed.
#' @param m Integer >= 1. Embedding dimension.
#' @param tau Integer >= 1. Time delay (in samples).
#'
#' @return A numeric matrix with `length(x) - (m - 1) * tau` rows and `m`
#'   columns; row `i` is `(x[i], x[i + tau], ..., x[i + (m - 1) * tau])`.
#'
#' @references
#' Takens F. Detecting strange attractors in turbulence. In: Dynamical
#' Systems and Turbulence, Warwick 1980. Lecture Notes in Mathematics,
#' vol 898. Springer; 1981:366-381.
#'
#' @examples
#' embed_time_series(sin(seq(0, 20, length.out = 200)), m = 3, tau = 5)
#'
#' @export
embed_time_series <- function(x, m = 2L, tau = 1L) {
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  m <- as.integer(m)
  tau <- as.integer(tau)
  if (m < 1L) stop("`m` must be >= 1.", call. = FALSE)
  if (tau < 1L) stop("`tau` must be >= 1.", call. = FALSE)

  N <- length(x)
  n_vec <- N - (m - 1L) * tau
  if (n_vec < 1L) {
    stop("`x` is too short for the requested `m`/`tau`.", call. = FALSE)
  }

  traj <- matrix(NA_real_, nrow = n_vec, ncol = m)
  for (d in seq_len(m)) {
    traj[, d] <- x[seq_len(n_vec) + (d - 1L) * tau]
  }
  traj
}

#' Estimate Embedding Time Delay
#'
#' Estimates a suitable time delay `tau` for phase-space reconstruction,
#' via the first local minimum of the average mutual information (Fraser &
#' Swinney 1986) -- the standard choice, more robust to nonlinearity than
#' autocorrelation. An autocorrelation-based fallback (first zero crossing)
#' is also available.
#'
#' @param x Numeric vector. The time series to analyse.
#' @param max_lag Integer or `NULL`. Largest lag considered. `NULL`
#'   (default) uses `min(length(x) %/% 10, 100)`.
#' @param n_bins Integer. Number of equal-width bins used to discretise `x`
#'   for the mutual-information histogram. Default `16`.
#' @param method `"mutual_information"` (default) or `"acf"`.
#'
#' @return Integer: the estimated delay `tau`.
#'
#' @references
#' Fraser AM, Swinney HL. Independent coordinates for strange attractors
#' from mutual information. Phys Rev A 1986;33:1134-1140.
#'
#' @examples
#' set.seed(1)
#' estimate_delay(sin(seq(0, 40, length.out = 1000)) + rnorm(1000, sd = 0.05))
#'
#' @export
estimate_delay <- function(x, max_lag = NULL, n_bins = 16L,
                            method = c("mutual_information", "acf")) {
  method <- match.arg(method)
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  N <- length(x)
  if (is.null(max_lag)) max_lag <- min(N %/% 10L, 100L)
  max_lag <- as.integer(max_lag)
  if (max_lag < 1L) stop("`max_lag` must be >= 1.", call. = FALSE)

  if (method == "acf") {
    ac <- stats::acf(x, lag.max = max_lag, plot = FALSE)$acf[-1]
    below <- which(ac <= 0)
    return(if (length(below) > 0L) below[1] else max_lag)
  }

  breaks <- seq(min(x), max(x), length.out = as.integer(n_bins) + 1L)
  breaks[length(breaks)] <- breaks[length(breaks)] + 1e-9
  bin <- findInterval(x, breaks, all.inside = TRUE)

  ami <- numeric(max_lag)
  for (lag in seq_len(max_lag)) {
    b1 <- bin[seq_len(N - lag)]
    b2 <- bin[(lag + 1L):N]
    joint <- table(b1, b2)
    pij <- joint / sum(joint)
    pi_ <- rowSums(pij)
    pj_ <- colSums(pij)
    nz <- pij > 0
    denom <- outer(pi_, pj_)[nz]
    ami[lag] <- sum(pij[nz] * log(pij[nz] / denom))
  }

  d <- diff(ami)
  first_min <- which(d > 0)[1]
  if (is.na(first_min)) which.min(ami) else first_min
}

#' Estimate Embedding Dimension via False Nearest Neighbours
#'
#' Estimates a suitable embedding dimension `m` for phase-space
#' reconstruction using Kennel, Brown & Abarbanel's (1992) false-nearest-
#' neighbours (FNN) method: for each candidate dimension, the fraction of
#' points whose nearest neighbour turns out not to be a true dynamical
#' neighbour (rather than an artefact of projecting a higher-dimensional
#' attractor down too far) is computed; `m` is chosen where this fraction
#' first drops below `fnn_threshold`.
#'
#' This search is O(N^2) per candidate dimension (brute-force nearest
#' neighbour), so a warning is issued above 5000 points.
#'
#' @param x Numeric vector. The time series to analyse.
#' @param tau Integer. Time delay, e.g. from [estimate_delay()]. Default `1`.
#' @param m_max Integer. Largest embedding dimension to test. Default `10`.
#' @param rtol Numeric. Threshold for the FNN distance-ratio criterion.
#'   Default `15`, as in the original paper.
#' @param atol Numeric. Threshold for the FNN attractor-size criterion.
#'   Default `2`, as in the original paper.
#' @param theiler Integer. Temporal exclusion window: candidate neighbours
#'   within `theiler` samples of a point in time are ignored, to avoid
#'   picking a temporally- rather than dynamically-close neighbour. Default
#'   `0` (no exclusion).
#' @param fnn_threshold Numeric in `[0, 1]`. `m` is chosen at the first
#'   dimension whose FNN fraction is at or below this value. Default `0.01`.
#'
#' @return A list with:
#'   \item{m}{The estimated embedding dimension (`NA` if the FNN fraction
#'     never drops to `fnn_threshold` within `m_max`).}
#'   \item{dim}{The dimensions tested, `1:m_max`.}
#'   \item{fnn_fraction}{FNN fraction at each tested dimension.}
#'
#' @references
#' Kennel MB, Brown R, Abarbanel HDI. Determining embedding dimension for
#' phase-space reconstruction using a geometrical construction. Phys Rev A
#' 1992;45:3403-3411.
#'
#' @examples
#' set.seed(1)
#' x <- sin(seq(0, 40, length.out = 800)) + rnorm(800, sd = 0.02)
#' estimate_embed_dim(x, m_max = 6)
#'
#' @export
estimate_embed_dim <- function(x, tau = 1L, m_max = 10L, rtol = 15, atol = 2,
                                theiler = 0L, fnn_threshold = 0.01) {
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  N <- length(x)
  if (N > 5000L) {
    warning(
      "False-nearest-neighbours search is O(N^2) per dimension; ", N,
      " points may be slow. Consider subsampling `x` first.", call. = FALSE
    )
  }

  dims <- seq_len(as.integer(m_max))
  frac <- fnn_fraction_cpp(
    as.double(x), as.integer(tau), as.integer(dims),
    as.double(rtol), as.double(atol), as.integer(theiler)
  )

  below <- which(frac <= fnn_threshold)
  m <- if (length(below) > 0L) dims[below[1]] else NA_integer_

  list(m = m, dim = dims, fnn_fraction = frac)
}
