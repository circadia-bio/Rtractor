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

#' Multifractal Detrending Moving Average (MFDMA)
#'
#' Estimates the multifractal scaling properties of a time series using
#' the detrending moving average algorithm (Gu & Zhou 2010). Clean-room
#' C++ reimplementation from the published algorithm (the reference
#' MATLAB implementation consulted, `MFDMA_1D.m`, had no license header;
#' see `inst/COPYRIGHTS`). The segment-fluctuation core was validated
#' against a Python transliteration of that reference on synthetic test
#' data (exact match to displayed precision).
#'
#' @param x Numeric vector. The time series to analyse.
#' @param n_min,n_max Integer. Lower/upper bound of the segment size `n`.
#'   Following the reference implementation's guidance: `n_min` around
#'   `10`; `n_max` around 10% of `length(x)`. Defaults: `n_min = 10`,
#'   `n_max = round(length(x) / 10)`.
#' @param n_scales Integer. Number of segment sizes to evaluate
#'   (log-spaced between `n_min` and `n_max`). Default `30`.
#' @param theta Numeric in `[0, 1]`. Position of the moving-average window:
#'   `0` (default, recommended) = backward MFDMA, `0.5` = centered,
#'   `1` = forward.
#' @param q Numeric vector. Multifractal orders to evaluate. Default
#'   `seq(-4, 4, by = 0.1)`.
#'
#' @return A list with:
#'   \item{n}{Segment sizes evaluated.}
#'   \item{Fq}{Matrix of the q-th order fluctuation function (segment size
#'     x q).}
#'   \item{tau}{Multifractal scaling exponent tau(q).}
#'   \item{alpha}{Singularity strength alpha(q) (trimmed at both ends by
#'     the local-slope smoothing window; shorter than `q`).}
#'   \item{f}{Multifractal spectrum f(alpha).}
#'   \item{q}{The `q` values corresponding to `alpha`/`f` (trimmed to match).}
#'
#' @references
#' Gu GF, Zhou WX. Detrending moving average algorithm for multifractals.
#' Phys Rev E 2010;82:011136.
#'
#' @examples
#' set.seed(1)
#' x <- rnorm(4000)
#' res <- mfdma(x, n_min = 10, n_max = 400, n_scales = 20)
#' plot(res$alpha, res$f, type = "b")
#'
#' @export
mfdma <- function(x, n_min = 10L, n_max = NULL, n_scales = 30L, theta = 0,
                   q = seq(-4, 4, by = 0.1)) {
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  if (theta < 0 || theta > 1) stop("`theta` must be in [0, 1].", call. = FALSE)

  x <- as.double(x)
  M <- length(x)
  if (is.null(n_max)) n_max <- round(M / 10)
  n_min <- as.integer(n_min)
  n_max <- as.integer(n_max)
  if (n_min < 2L) stop("`n_min` must be >= 2.", call. = FALSE)
  if (n_max <= n_min) stop("`n_max` must be greater than `n_min`.", call. = FALSE)

  scales <- unique(as.integer(round(
    10^seq(log10(n_min), log10(n_max), length.out = n_scales)
  )))

  y <- cumsum(x)
  F_list <- mfdma_fluctuations_cpp(y, scales, theta)

  n_q <- length(q)
  n_scales_actual <- length(scales)
  Fq <- matrix(NA_real_, nrow = n_scales_actual, ncol = n_q)
  for (j in seq_len(n_scales_actual)) {
    f <- F_list[[j]]
    f <- f[is.finite(f) & f > 0]
    if (length(f) == 0L) next
    for (i in seq_len(n_q)) {
      qi <- q[i]
      Fq[j, i] <- if (qi == 0) {
        exp(0.5 * mean(log(f^2)))
      } else {
        (mean(f^qi))^(1 / qi)
      }
    }
  }

  logn <- log(scales)
  h <- vapply(seq_len(n_q), function(i) {
    fq_col <- Fq[, i]
    ok <- is.finite(fq_col) & fq_col > 0
    if (sum(ok) < 2L) return(NA_real_)
    unname(stats::coef(stats::lm(log(fq_col[ok]) ~ logn[ok]))[2])
  }, numeric(1))
  tau <- h * q - 1

  dx <- (7L - 1L) %/% 2L
  n_tau <- length(tau)
  if (n_tau <= 2L * dx) {
    stop("`q` must have more than ", 2L * dx, " values for the alpha/f smoothing window.", call. = FALSE)
  }
  alpha <- rep(NA_real_, n_tau)
  for (i in (dx + 1L):(n_tau - dx)) {
    idx <- (i - dx):(i + dx)
    alpha[i] <- unname(stats::coef(stats::lm(tau[idx] ~ q[idx]))[2])
  }
  keep <- (dx + 1L):(n_tau - dx)
  alpha_trim <- alpha[keep]
  tau_trim <- tau[keep]
  q_trim <- q[keep]
  f_spectrum <- q_trim * alpha_trim - tau_trim

  list(n = scales, Fq = Fq, tau = tau, alpha = alpha_trim, f = f_spectrum, q = q_trim)
}

#' Multifractal Spectrum via the Chhabra-Jensen Method
#'
#' Estimates the multifractal spectrum f(alpha) and generalised dimension
#' spectrum D(q) of a strictly positive time series using the direct
#' box-counting method of Chhabra & Jensen (1989). Clean-room C++
#' reimplementation from the published algorithm (the reference MATLAB
#' implementation consulted, `ChhabraJensen_Yuj_w0.m`, co-authored by
#' L. Franca, had no license header; see `inst/COPYRIGHTS`). The moments
#' core was validated against a Python transliteration of that reference
#' on synthetic test data (exact match to displayed precision).
#'
#' @param x Numeric vector, strictly positive (treat it as a measure —
#'   e.g. apply a sigmoid transform first if your data can be negative).
#'   `length(x)` must be evenly divisible by `2^scale` for every value in
#'   `scales`, i.e. dyadic lengths (powers of two) work best.
#' @param q_values Numeric vector of multifractal orders. Per the original
#'   author's guidance, values strictly between 0 and 1 (exclusive) are
#'   numerically unstable for this method and best avoided -- a warning is
#'   issued if any are supplied. Default skips that range:
#'   `c(seq(-10, -0.1, by = 0.1), seq(1, 10, by = 0.1))`.
#' @param scales Integer vector of box-counting scale exponents; the box
#'   size at each scale is `2^scale`. Default `1:floor(log2(length(x)) - 2)`,
#'   which keeps at least 4 points per box at the coarsest scale.
#'
#' @return A list with:
#'   \item{alpha}{Singularity strength alpha(q).}
#'   \item{falpha}{Multifractal spectrum f(alpha(q)).}
#'   \item{Dq}{Generalised dimension spectrum D(q).}
#'   \item{r_squared_alpha, r_squared_falpha, r_squared_Dq}{R-squared of
#'     the linear fit underlying each of the above, per `q` -- inspect
#'     these before trusting a given q value's estimate.}
#'   \item{q}{The `q_values` used.}
#'   \item{mu_scale, Ma, Mf, Md}{The underlying regression inputs, included
#'     for completeness.}
#'
#' @references
#' Chhabra A, Jensen RV. Direct determination of the f(alpha) singularity
#' spectrum. Phys Rev Lett 1989;62:1327-1330.
#'
#' Franca LGS, Miranda JGV, Leite M, Sharma NK, Walker MC, Lemieux L,
#' Wang Y. Fractal and multifractal properties of electrographic
#' recordings of human brain activity: toward its use as a signal feature
#' for machine learning in clinical applications. Front Physiol
#' 2018;9:1767. Compares MF-DFA, MF-DMA, and Chhabra-Jensen on simulated
#' and human intracranial EEG data and finds Chhabra-Jensen the most
#' stable/reliable of the three -- the basis for implementing this method
#' (rather than MF-DFA) alongside `mfdma()` in Rtractor.
#'
#' @examples
#' set.seed(1)
#' x <- abs(rnorm(1024)) + 0.01
#' res <- chhabra_jensen(x, scales = 1:6)
#' plot(res$alpha, res$falpha, type = "b")
#'
#' @export
chhabra_jensen <- function(x,
                           q_values = c(seq(-10, -0.1, by = 0.1), seq(1, 10, by = 0.1)),
                           scales = NULL) {
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  if (any(x <= 0)) {
    stop("`x` must be strictly positive; consider a sigmoid transform first.", call. = FALSE)
  }
  if (any(q_values > 0 & q_values < 1)) {
    warning("Some `q_values` fall strictly between 0 and 1; this method is numerically unstable there (per the original author's guidance).", call. = FALSE)
  }

  x <- as.double(x)
  L <- length(x)
  if (is.null(scales)) scales <- seq_len(max(1L, floor(log2(L)) - 2L))
  scales <- as.integer(scales)

  for (s in scales) {
    window <- 2^s
    if (L %% window != 0) {
      stop(sprintf(
        "length(x) (%d) is not divisible by 2^%d = %d; choose `scales` that evenly divide length(x).",
        L, s, window
      ), call. = FALSE)
    }
  }

  moments <- chhabra_jensen_moments_cpp(x, as.double(q_values), scales)
  Ma <- moments$Ma
  Mf <- moments$Mf
  Md <- moments$Md

  mu_scale <- -log10(2^scales)
  nq <- length(q_values)

  alpha <- numeric(nq)
  falpha <- numeric(nq)
  Dq <- numeric(nq)
  r2_alpha <- numeric(nq)
  r2_falpha <- numeric(nq)
  r2_Dq <- numeric(nq)

  for (i in seq_len(nq)) {
    fit_a <- stats::lm(Ma[i, ] ~ mu_scale)
    fit_f <- stats::lm(Mf[i, ] ~ mu_scale)
    fit_d <- stats::lm(Md[i, ] ~ mu_scale)

    alpha[i]  <- unname(stats::coef(fit_a)[2])
    falpha[i] <- unname(stats::coef(fit_f)[2])
    b_md      <- unname(stats::coef(fit_d)[2])

    q <- q_values[i]
    Dq[i] <- if (q > 0 && q <= 1) b_md else b_md / (q - 1)

    # Computed directly rather than via summary.lm()$r.squared: these
    # regressions are frequently a near-perfect fit (the data closely
    # follows the theoretical power law), which makes summary.lm() emit
    # "essentially perfect fit: summary may be unreliable" -- a base R
    # warning about its own internal numerics, not a problem with the R^2
    # value itself. Computing R^2 from residuals directly avoids that
    # code path (and the warning) entirely.
    r2_alpha[i]  <- .r_squared(fit_a, Ma[i, ])
    r2_falpha[i] <- .r_squared(fit_f, Mf[i, ])
    r2_Dq[i]     <- .r_squared(fit_d, Md[i, ])
  }

  list(
    alpha = alpha, falpha = falpha, Dq = Dq,
    r_squared_alpha = r2_alpha, r_squared_falpha = r2_falpha, r_squared_Dq = r2_Dq,
    q = q_values, mu_scale = mu_scale, Ma = Ma, Mf = Mf, Md = Md
  )
}

#' @keywords internal
.r_squared <- function(fit, y) {
  yhat <- stats::fitted(fit)
  ss_res <- sum((y - yhat)^2)
  ss_tot <- sum((y - mean(y))^2)
  if (ss_tot == 0) return(1)
  1 - ss_res / ss_tot
}

#' Petrosian Fractal Dimension
#'
#' Estimates fractal dimension from the rate of sign changes in a time
#' series' first difference (a fast proxy for signal irregularity).
#' Ported from Lucas Franca's own `mrpheus` package (AASM staging feature
#' pipeline), itself validated for exact parity against the `antropy`
#' Python library; re-validated here directly against `antropy` 0.2.2 on
#' synthetic test data (exact match). See `inst/COPYRIGHTS`.
#'
#' @param x Numeric vector. The time series to analyse.
#'
#' @return A length-1 numeric: the Petrosian fractal dimension.
#'
#' @references
#' Petrosian A. Kolmogorov complexity of finite sequences and recognition
#' of different preictal EEG patterns. Proceedings of the Eighth IEEE
#' Symposium on Computer-Based Medical Systems 1995:212-217.
#'
#' @examples
#' set.seed(1)
#' petrosian_fd(rnorm(1000))
#'
#' @export
petrosian_fd <- function(x) {
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  petrosian_fd_cpp(as.double(x))
}

#' Hjorth Mobility and Complexity
#'
#' Computes the Hjorth parameters (Hjorth 1970): mobility, a proxy for mean
#' frequency derived from the variance ratio of the first difference to
#' the signal itself; and complexity, a proxy for bandwidth/irregularity
#' derived from the second difference. Ported from Lucas Franca's own
#' `mrpheus` package (AASM staging feature pipeline). Uses Bessel-corrected
#' (`ddof = 1`) variance throughout, deliberately matching R's `var()`
#' convention -- this was mrpheus's original design choice, and differs
#' slightly from the `antropy` Python library's population-variance
#' (`ddof = 0`) convention. The two converge as `length(x)` grows but
#' differ meaningfully for short series; see `inst/COPYRIGHTS`.
#'
#' @param x Numeric vector. The time series to analyse.
#'
#' @return A named list with `mobility` and `complexity`.
#'
#' @references
#' Hjorth B. EEG analysis based on time domain properties.
#' Electroencephalogr Clin Neurophysiol 1970;29(3):306-310.
#'
#' @examples
#' set.seed(1)
#' hjorth_parameters(rnorm(1000))
#'
#' @export
hjorth_parameters <- function(x) {
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  res <- hjorth_cpp(as.double(x))
  list(mobility = unname(res["mobility"]), complexity = unname(res["complexity"]))
}

#' Number of Zero Crossings
#'
#' Counts sign changes in a time series -- a simple time-domain proxy for
#' dominant frequency / signal roughness, commonly reported alongside
#' Hjorth parameters and fractal dimension in EEG complexity work. Ported
#' from Lucas Franca's own `mrpheus` package (AASM staging feature
#' pipeline), itself validated for exact parity against the `antropy`
#' Python library; re-validated here directly against `antropy` 0.2.2 on
#' synthetic test data (exact match). See `inst/COPYRIGHTS`.
#'
#' @param x Numeric vector. The time series to analyse.
#'
#' @return A length-1 integer: the number of zero crossings.
#'
#' @examples
#' set.seed(1)
#' num_zerocross(rnorm(1000))
#'
#' @export
num_zerocross <- function(x) {
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  num_zerocross_cpp(as.double(x))
}

#' Hurst Exponent via Standard Deviation Analysis (SDA)
#'
#' Estimates the Hurst exponent of a time series from the scaling of its
#' Root-Mean-Square deviation from the local mean across a range of window
#' sizes -- no detrending (Russ 1994). This is the method used to compute
#' the Hurst exponent of body-sway time series in da Silva Costa et al.
#' (2017), a fibromyalgia gait/balance study co-authored by Lucas Franca.
#'
#' Ported from Franca's own legacy `mrug` C++ tool (used to produce that
#' paper's figures), with three bugs fixed rather than reproduced -- see
#' `inst/COPYRIGHTS` for the full write-up. Treat this as a *corrected*
#' re-implementation of the method, not a bit-exact match to the 2017
#' paper's published numbers (the exact code snapshot used for that paper
#' could not be identified).
#'
#' @param x Numeric vector. The time series to analyse, used directly (not
#'   as increments) -- e.g. a displacement/position series, as in the
#'   original body-sway application.
#' @param scale_min Integer. Smallest window size. Default `5`, matching
#'   the original implementation.
#' @param scale_max Integer or `NULL`. Largest window size evaluated.
#'   `NULL` (default) uses `floor(length(x) / 2)`, matching the original
#'   tool's convention of scanning up to half the series length.
#' @param growth Numeric > 1. Geometric growth factor applied between
#'   successive window sizes. Default `1.1`, matching the original.
#' @param force_below Integer. Below this window size, the window size is
#'   incremented by 1 rather than by `growth`, to avoid repeated integer
#'   window sizes at small scales (geometric growth by 1.1 often doesn't
#'   move an integer forward until it's larger). Default `20`, matching
#'   the original.
#' @param fit_min,fit_max Integer. Window-size range used for the log-log
#'   linear fit that yields the Hurst exponent. `fit_min` defaults to
#'   `scale_min`; `fit_max` defaults to `floor(length(x) / 4)`, matching
#'   the convention used in the CLI tool that produced the 2017 paper's
#'   results.
#'
#' @return A list with:
#'   \item{k}{Window sizes evaluated (integer vector).}
#'   \item{F}{RMS fluctuation at each window size (numeric vector; `NA`
#'     for any window size with no valid non-constant windows).}
#'   \item{H}{The Hurst exponent: the slope of `log(F)` on `log(k)`, fit
#'     over `[fit_min, fit_max]`.}
#'   \item{se}{Standard error of `H`.}
#'   \item{r_squared}{R-squared of the log-log fit.}
#'   \item{intercept}{Intercept of the log-log fit.}
#'   \item{n_fit}{Number of window sizes used in the fit.}
#'
#' @references
#' Russ JC. Fractal Surfaces. New York, NY: Plenum; 1994.
#'
#' da Silva Costa I, Gamundi A, Miranda JGV, Franca LGS, De Santana CN,
#' Montoya P. Altered functional performance in patients with
#' fibromyalgia. Front Hum Neurosci 2017;11:14.
#'
#' @examples
#' set.seed(1)
#' sda(cumsum(rnorm(2000)), fit_max = 200)
#'
#' @export
sda <- function(x, scale_min = 5L, scale_max = NULL, growth = 1.1,
                force_below = 20L, fit_min = scale_min, fit_max = NULL) {
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  if (growth <= 1) stop("`growth` must be > 1.", call. = FALSE)

  x <- as.double(x)
  N <- length(x)
  scale_min <- as.integer(scale_min)
  if (is.null(scale_max)) scale_max <- N %/% 2L
  scale_max <- as.integer(scale_max)
  if (scale_max < scale_min) {
    stop("`scale_max` must be >= `scale_min`.", call. = FALSE)
  }
  if (is.null(fit_max)) fit_max <- floor(N / 4)
  fit_min <- as.integer(fit_min)
  fit_max <- as.integer(fit_max)

  scales <- .sda_scales(scale_min, scale_max, growth, as.integer(force_below))
  Fk <- sda_fluctuation_cpp(x, scales)

  keep <- scales >= fit_min & scales <= fit_max & is.finite(Fk) & Fk > 0
  n_fit <- sum(keep)
  if (n_fit < 2L) {
    stop(
      "Fewer than 2 valid window sizes fall in [fit_min, fit_max] = [",
      fit_min, ", ", fit_max, "]; widen the range or check `x`.",
      call. = FALSE
    )
  }

  fit <- stats::lm(log(Fk[keep]) ~ log(scales[keep]))
  co <- unname(stats::coef(fit))
  se <- unname(sqrt(diag(stats::vcov(fit))))

  list(
    k = scales, F = Fk,
    H = co[2], se = se[2],
    r_squared = .r_squared(fit, log(Fk[keep])),
    intercept = co[1],
    n_fit = n_fit
  )
}

#' @keywords internal
.sda_scales <- function(scale_min, scale_max, growth, force_below) {
  k <- scale_min
  kaux <- force_below
  scales <- integer(0)
  while (k <= scale_max) {
    scales <- c(scales, k)
    k <- as.integer(growth * k)
    if (k <= kaux) {
      k <- k + 1L
      kaux <- k
    }
  }
  scales
}

# ---- Planned, not yet implemented -----------------------------------------
#
#   - box_counting_fd()      Box-counting fractal dimension
#
# Native counterpart(s) expected in: src/fractal.cpp (for box_counting_fd)
