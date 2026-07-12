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
#'   entropy in nats.
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

# ---- Planned, not yet implemented -----------------------------------------
#
#   - sample_entropy()       Sample entropy (SampEn), Richman & Moorman 2000
#   - approx_entropy()       Approximate entropy (ApEn), Pincus 1991
#   - shannon_entropy()      Shannon entropy of a discretised/binned signal
#   - renyi_entropy()        Renyi entropy, parametrised by order q
#
# Native counterpart(s) expected in: src/entropy.cpp
