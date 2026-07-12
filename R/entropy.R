# Entropy measures
#
# Family: entropy
# Signal-agnostic (any numeric vector, regardless of source — EEG,
# actigraphy, BOLD, HRV, etc.).

#' Permutation Entropy
#'
#' Estimates the complexity of a time series using permutation entropy
#' (Bandt & Pompe 2002): the Shannon entropy of the distribution of
#' ordinal patterns ("motifs") of length `order` found in the series.
#' Ported from Lucas Franca's own `mrpheus` package (part of its AASM
#' sleep-staging feature pipeline), itself validated for exact parity
#' against the `antropy` Python library; re-validated here directly
#' against `antropy` 0.2.2 on synthetic test data (exact match to
#' displayed precision). See `inst/COPYRIGHTS`.
#'
#' @param x Numeric vector. The time series to analyse.
#' @param order Integer >= 2. Length of the ordinal pattern (embedding
#'   dimension). Default `3`.
#' @param delay Integer >= 1. Time delay between pattern elements.
#'   Default `1`.
#' @param normalize Logical. If `TRUE` (default), divide by `log(order!)`
#'   so the result falls in `[0, 1]`. If `FALSE`, return the raw Shannon
#'   entropy in nats (natural log). Note this differs from the `antropy`
#'   Python library, which computes the raw value in bits (log base 2);
#'   the normalized value is base-independent and matches `antropy`
#'   exactly, but the raw values are not directly comparable between the
#'   two. See `inst/COPYRIGHTS`.
#'
#' @return A length-1 numeric: the permutation entropy.
#'
#' @references
#' Bandt C, Pompe B. Permutation entropy: a natural complexity measure for
#' time series. Phys Rev Lett 2002;88:174102.
#'
#' @examples
#' set.seed(1)
#' perm_entropy(rnorm(1000))
#'
#' @export
perm_entropy <- function(x, order = 3L, delay = 1L, normalize = TRUE) {
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  order <- as.integer(order)
  delay <- as.integer(delay)
  if (order < 2L) stop("`order` must be >= 2.", call. = FALSE)
  if (delay < 1L) stop("`delay` must be >= 1.", call. = FALSE)

  h <- perm_entropy_raw_cpp(as.double(x), order, delay)
  if (isTRUE(normalize)) h <- h / sum(log(seq_len(order)))
  h
}

#' Sample Entropy (SampEn)
#'
#' Estimates the complexity of a time series using sample entropy
#' (Richman & Moorman 2000): the negative log ratio of template matches
#' of length `m + 1` to template matches of length `m`, using a fixed
#' Chebyshev-distance tolerance. Direct C++ port of the counting core in
#' PhysioNet's reference `mse.c` (Costa), validated to reproduce the
#' compiled reference binary's output exactly (to its own displayed
#' precision) on synthetic test data. See `inst/COPYRIGHTS`.
#'
#' This is also the building block `multiscale_entropy()` applies at each
#' coarse-grained scale -- see `R/multiscale.R`.
#'
#' @param x Numeric vector. The time series to analyse.
#' @param m Integer >= 1. Template length. Default `2` (the standard MSE
#'   convention, per Costa et al.).
#' @param r Numeric > 0. Tolerance, as a fraction of `x`'s own standard
#'   deviation (i.e. the actual Chebyshev-distance tolerance used is
#'   `sd(x) * r`). Default `0.15` (the standard MSE convention).
#'
#' @return A length-1 numeric: the sample entropy. If there are too few
#'   valid template pairs to compare (`length(x) - m < 2`), returns `NA`.
#'   If there are zero matches at either template length (undefined ratio),
#'   returns the conventional fallback `-log(1 / (N*(N-1)))` used by the
#'   reference implementation, where `N = length(x) - m`.
#'
#' @references
#' Richman JS, Moorman JR. Physiological time-series analysis using
#' approximate entropy and sample entropy. Am J Physiol Heart Circ
#' Physiol 2000;278(6):H2039-H2049.
#'
#' @examples
#' set.seed(1)
#' sample_entropy(rnorm(1000))
#'
#' @export
sample_entropy <- function(x, m = 2L, r = 0.15) {
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  m <- as.integer(m)
  if (m < 1L) stop("`m` must be >= 1.", call. = FALSE)
  if (r <= 0) stop("`r` must be > 0.", call. = FALSE)

  x <- as.double(x)
  tolerance <- r * stats::sd(x)
  .sample_entropy_core(x, m, tolerance)
}

#' @keywords internal
.sample_entropy_core <- function(x, m, tolerance) {
  n <- length(x)
  nlin_m <- n - m
  if (nlin_m < 2L) return(NA_real_)

  counts <- sample_entropy_counts_cpp(x, m, tolerance)
  cont_m <- counts[m]
  cont_m1 <- counts[m + 1]

  if (cont_m == 0 || cont_m1 == 0) {
    return(-log(1 / (nlin_m * (nlin_m - 1))))
  }
  -log(cont_m1 / cont_m)
}

# ---- Planned, not yet implemented -----------------------------------------
#
#   - approx_entropy()       Approximate entropy (ApEn), Pincus 1991
#   - shannon_entropy()      Shannon entropy of a discretised/binned signal
#   - renyi_entropy()        Renyi entropy, parametrised by order q
#
# Native counterpart(s) expected in: src/entropy.cpp
