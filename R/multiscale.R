# Multiscale metrics
#
# Family: multiscale
# Complexity measures computed across a range of temporal scales
# (coarse-graining). Signal-agnostic (see entropy.R note). Built on the
# entropy family (entropy.R) applied at each coarse-grained scale, rather
# than a separate algorithmic base.

#' Multiscale Entropy (MSE)
#'
#' Estimates the complexity of a time series across a range of temporal
#' scales (Costa, Goldberger & Peng 2002): the series is coarse-grained
#' (non-overlapping block-averaged) at each scale factor, and
#' [sample_entropy()] is computed on each coarse-grained series, using a
#' tolerance held *fixed* relative to the original (not the coarse-grained)
#' series' standard deviation -- this is the standard MSE convention and
#' essential for entropy values to be comparable across scales. Direct
#' C++ port of the coarse-graining + sample-entropy core in PhysioNet's
#' reference `mse.c` (Costa), validated to reproduce the compiled
#' reference binary's output exactly (to its own displayed precision) on
#' synthetic test data. See `inst/COPYRIGHTS`.
#'
#' @param x Numeric vector. The time series to analyse.
#' @param scale_max Integer >= 1. Largest scale factor to evaluate; MSE is
#'   computed at every integer scale from `1` to `scale_max`. Default `20`
#'   (the standard MSE convention).
#' @param m Integer >= 1. Template length, passed to [sample_entropy()].
#'   Default `2` (the standard MSE convention).
#' @param r Numeric > 0. Tolerance as a fraction of the *original* series'
#'   standard deviation, passed to [sample_entropy()]. Default `0.15` (the
#'   standard MSE convention).
#'
#' @return A list with:
#'   \item{scale}{The scale factors evaluated, `1:scale_max`.}
#'   \item{mse}{Sample entropy at each scale (may contain `NA` at large
#'     scales, where the coarse-grained series becomes too short to
#'     estimate reliably).}
#'   \item{m, r}{The parameters used, echoed back for reference.}
#'
#' @references
#' Costa M, Goldberger AL, Peng CK. Multiscale entropy analysis of complex
#' physiologic time series. Phys Rev Lett 2002;89:068102.
#'
#' @examples
#' set.seed(1)
#' res <- multiscale_entropy(rnorm(2000), scale_max = 10)
#' plot(res$scale, res$mse, type = "b")
#'
#' @export
multiscale_entropy <- function(x, scale_max = 20L, m = 2L, r = 0.15) {
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  scale_max <- as.integer(scale_max)
  m <- as.integer(m)
  if (scale_max < 1L) stop("`scale_max` must be >= 1.", call. = FALSE)
  if (m < 1L) stop("`m` must be >= 1.", call. = FALSE)
  if (r <= 0) stop("`r` must be > 0.", call. = FALSE)

  x <- as.double(x)
  tolerance <- r * stats::sd(x) # fixed across scales, per Costa et al.'s convention

  scales <- seq_len(scale_max)
  mse_values <- vapply(scales, function(scale) {
    y <- .coarse_grain(x, scale)
    .sample_entropy_core(y, m, tolerance)
  }, numeric(1))

  list(scale = scales, mse = mse_values, m = m, r = r)
}

#' @keywords internal
.coarse_grain <- function(x, scale) {
  n <- length(x)
  n_out <- n %/% scale
  if (n_out < 1L) return(numeric(0))
  vapply(seq_len(n_out), function(i) {
    mean(x[((i - 1L) * scale + 1L):(i * scale)])
  }, numeric(1))
}

# ---- Planned, not yet implemented -----------------------------------------
#
#   - rcmse()                Refined composite multiscale entropy,
#                            Wu et al. 2014
#   - mse_complexity_index() Area-under-curve summary of an MSE profile
