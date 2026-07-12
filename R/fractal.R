# Fractal and multifractal analysis
#
# Family: fractal
# Wraps established reference implementations for fractal dimension and
# scaling-exponent estimation. Signal-agnostic (any numeric vector,
# regardless of source — EEG, actigraphy, BOLD, HRV, etc.).

#' Detrended Fluctuation Analysis (DFA)
#'
#' Estimates the scaling exponent alpha of a time series using Detrended
#' Fluctuation Analysis (Peng et al. 1994). This is a direct C++ port of
#' the reference `dfa.c` implementation distributed by PhysioNet (Mietus,
#' Peng & Moody), validated to reproduce the original compiled binary's
#' output exactly on synthetic test data — see `inst/COPYRIGHTS`.
#'
#' @param x Numeric vector. The time series to analyse.
#' @param order Integer >= 1. Order of the polynomial detrending fit
#'   (`1` = linear detrending, the DFA default; `2` = quadratic, etc.).
#'   Corresponds to `nfit - 1` in the original `dfa.c`.
#' @param min_box,max_box Integer or `NULL`. Smallest/largest box size to
#'   evaluate. Defaults (as in the original): `min_box = 2 * (order + 1)`,
#'   `max_box = length(x) / 4`.
#' @param integrate Logical. If `TRUE` (default), `x` is treated as an
#'   increment series and cumulatively summed before analysis (the usual
#'   DFA convention — e.g. pass RR-interval deviations, not a cumulative
#'   profile). Set `FALSE` if `x` is already an integrated/cumulative
#'   profile.
#' @param sliding_window Logical. Use overlapping (sliding) windows instead
#'   of non-overlapping boxes. Default `FALSE`.
#'
#' @return A list with:
#'   \item{n}{Box sizes evaluated (integer vector).}
#'   \item{F}{RMS fluctuation at each box size (numeric vector).}
#'   \item{alpha}{The DFA scaling exponent: the slope of `log10(F)` on
#'     `log10(n)`.}
#'
#' @references
#' Peng C-K, Buldyrev SV, Havlin S, Simons M, Stanley HE, Goldberger AL.
#' Mosaic organization of DNA nucleotides. Phys Rev E 1994;49:1685-1689.
#'
#' @examples
#' set.seed(1)
#' dfa(rnorm(2000))
#'
#' @export
dfa <- function(x, order = 1, min_box = NULL, max_box = NULL,
                 integrate = TRUE, sliding_window = FALSE) {
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  if (order < 1) stop("`order` must be >= 1.", call. = FALSE)

  x <- as.double(x)
  nfit <- as.integer(order + 1)
  minbox <- if (is.null(min_box)) 0L else as.integer(min_box)
  maxbox <- if (is.null(max_box)) 0L else as.integer(max_box)

  res <- dfa_cpp(
    x, nfit = nfit, minbox = minbox, maxbox = maxbox,
    integrate = integrate, sliding_window = sliding_window
  )

  fit <- stats::lm(log10(res$F) ~ log10(res$n))
  alpha <- unname(stats::coef(fit)[2])

  list(n = res$n, F = res$F, alpha = alpha)
}

#' Higuchi Fractal Dimension
#'
#' Estimates the fractal dimension of a time series using Higuchi's (1988)
#' curve-length algorithm. A clean-room C++ reimplementation, validated
#' against a MATLAB reference implementation on synthetic test data — see
#' `inst/COPYRIGHTS`.
#'
#' @param x Numeric vector. The time series to analyse.
#' @param k_max Integer (>= 2). Maximum sub-series interval `k`. Choose
#'   based on where the log-log curve plateaus for your signal (see Doyle
#'   et al. for guidance in postural-sway / physiological contexts).
#'
#' @return A list with:
#'   \item{k}{The `k` values evaluated, `1:k_max`.}
#'   \item{L}{Average curve length at each `k`.}
#'   \item{hfd}{The Higuchi Fractal Dimension: negative slope of `log(L)`
#'     on `log(1/k)`.}
#'
#' @references
#' Higuchi T. Approach to an irregular time series on the basis of the
#' fractal theory. Physica D 1988;31(2):277-283.
#'
#' @examples
#' set.seed(1)
#' higuchi_fd(rnorm(1000), k_max = 10)
#'
#' @export
higuchi_fd <- function(x, k_max = 10L) {
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  k_max <- as.integer(k_max)
  if (k_max < 2L) stop("`k_max` must be >= 2.", call. = FALSE)

  L <- higuchi_length_cpp(as.double(x), k_max)
  k <- seq_len(k_max)

  fit <- stats::lm(log(L) ~ log(1 / k))
  hfd <- unname(stats::coef(fit)[2])

  list(k = k, L = L, hfd = hfd)
}

# ---- Planned, not yet implemented -----------------------------------------
#
#   - box_counting_fd()      Box-counting fractal dimension
#   - mfdma()                Multifractal detrending moving average
#                            (reference: MFDMA_1D.m — no license header
#                            present; candidate for clean-room
#                            reimplementation)
#   - chhabra_jensen()       Multifractal spectrum via the Chhabra-Jensen
#                            method (reference: ChhabraJensen_Yuj_w0.m,
#                            co-authored by L. França — no license header
#                            present; candidate for clean-room
#                            reimplementation)
#
# Native counterpart(s) expected in: src/fractal.cpp (for mfdma/chhabra_jensen)
