# Synthetic multifractal/monofractal test-signal generators
#
# Family: simulate
# Reference-generator utilities, distinct from the measurement families
# (entropy, fractal, lyapunov, multiscale, rqa): these produce synthetic
# time series with known/controllable fractal properties, useful for
# testing and demonstrating the measurement functions above (e.g. feeding
# pmodel() output into mfdma() or chhabra_jensen()).

#' Simulate a Multifractal Time Series Using the p-Model
#'
#' Generates a multiplicative binomial cascade ("p-model") time series
#' (Meneveau & Sreenivasan 1987), optionally fractionally integrated in
#' Fourier space to a prescribed power-spectrum slope (Davis et al. 1997)
#' to produce a continuous, nonstationary series superposed on the
#' multifractal cascade.
#'
#' Clean-room R reimplementation of `pmodel.m`, written by the late Victor
#' Venema (no license header present in the source file; see
#' `inst/COPYRIGHTS`). The underlying random sign draws use R's own RNG
#' rather than attempting to replicate MATLAB's specific generator
#' bit-for-bit (not feasible across languages, and not meaningful for a
#' synthetic-data generator characterised by its cascade *rule* rather
#' than particular RNG bytes) -- see the `seed` argument below for how the
#' original's reseed-every-level behaviour is preserved. Validated against
#' two RNG-independent structural properties that hold for *any* random
#' sign sequence: at `p = 0.5` the cascade is provably an exact constant
#' vector of `1`s, and total mass is conserved exactly at every cascade
#' level (`sum(y)` doubles per level) for any `p` -- see
#' `tests/testthat/test-pmodel.R`.
#'
#' @param n_values Integer. Desired length of the output series. Internally
#'   rounded up to the next power of two for the cascade construction, then
#'   truncated back to `n_values`. Default `256`.
#' @param p Numeric in `(0, 1)`. The p-model parameter: values near `0` or
#'   `1` produce a very peaked (strongly multifractal) series; values near
#'   `0.5` produce a calm series (`p = 0.5` gives an exactly constant unit
#'   vector -- no multifractality at all). Default `0.375`.
#' @param slope Numeric or `NULL`. If supplied, the p-model series is
#'   fractionally integrated in Fourier space to have this power-spectrum
#'   slope, producing a continuous, nonstationary series (slopes between
#'   -1 and -3 give stationary increments; steeper than -3 is
#'   differentiable) -- see Davis et al. (1997) for the regimes. If `NULL`
#'   (default), no fractional integration is applied.
#' @param seed Integer or `NULL`. Seeds the random sign draws at *every*
#'   cascade level. The original MATLAB implementation reseeds to a fixed
#'   seed (`42`) at every level unconditionally, making the cascade fully
#'   deterministic by design (each level's random signs are literally a
#'   prefix of the next level's, since the RNG stream restarts from the
#'   same point every time). This port preserves that same
#'   reseed-every-level structure, with `seed` (default `42`, matching the
#'   original) controlling it. Pass `seed = NULL` to draw fresh randomness
#'   at every level instead (no reseeding) -- unlike the original, which
#'   has no way to disable this.
#'
#' @return If `slope` is `NULL`: a numeric vector of length `n_values`,
#'   the p-model series itself. If `slope` is supplied: a list with
#'   \item{x}{The fractionally integrated series.}
#'   \item{y}{The underlying p-model series before integration.}
#'
#' @references
#' Meneveau C, Sreenivasan KR. Simple multifractal cascade model for
#' fully developed turbulence. Phys Rev Lett 1987;59:1424-1427.
#'
#' Davis A, Marshak A, Cahalan R, Wiscombe W. The Landsat scale break in
#' stratocumulus as a three-dimensional radiative transfer effect:
#' implications for cloud remote sensing. J Atmos Sci 1997;54(2):241-260.
#'
#' Venema V, Bachner S, Rust HW, Simmer C. Statistical characteristics of
#' surrogate data based on geophysical measurements. Nonlin Processes
#' Geophys 2006;13:449-466. (Requested citation for this code by its
#' original author, Victor Venema.)
#'
#' @examples
#' y <- pmodel(1024, p = 0.3)
#' range(y)
#'
#' # p = 0.5 always gives a constant series, regardless of seed:
#' all(pmodel(64, p = 0.5) == 1)
#'
#' @export
pmodel <- function(n_values = 256L, p = 0.375, slope = NULL, seed = 42L) {
  if (!is.numeric(p) || length(p) != 1L || p <= 0 || p >= 1) {
    stop("`p` must be a single numeric value in (0, 1).", call. = FALSE)
  }
  n_values <- as.integer(n_values)
  if (n_values < 1L) stop("`n_values` must be >= 1.", call. = FALSE)

  no_orders <- ceiling(log2(n_values))
  y <- 1
  for (i in seq_len(no_orders)) {
    y <- .pmodel_next_step(y, p, seed = seed)
  }
  y <- y[seq_len(n_values)]

  if (is.null(slope)) {
    return(y)
  }

  x <- .pmodel_fractional_integration(y, slope)
  list(x = x, y = y)
}

#' @keywords internal
.pmodel_next_step <- function(y, p, seed = NULL) {
  len <- length(y)
  if (!is.null(seed)) {
    has_seed <- exists(".Random.seed", envir = .GlobalEnv)
    old_seed <- if (has_seed) .Random.seed else NULL
    on.exit({
      if (has_seed) assign(".Random.seed", old_seed, envir = .GlobalEnv)
      else if (exists(".Random.seed", envir = .GlobalEnv)) rm(".Random.seed", envir = .GlobalEnv)
    }, add = TRUE)
    set.seed(seed)
  }
  sign <- sample(c(-1, 1), len, replace = TRUE)

  y2 <- numeric(len * 2)
  y2[seq(1, 2 * len - 1, by = 2)] <- y + sign * (1 - 2 * p) * y
  y2[seq(2, 2 * len, by = 2)]     <- y - sign * (1 - 2 * p) * y
  y2
}

#' @keywords internal
.pmodel_fractal_spectrum <- function(n_values, slope) {
  half <- n_values / 2
  a <- numeric(n_values)
  for (t2 in 1:(half + 1)) {
    index <- t2 - 1
    t4 <- 2 + n_values - t2
    if (t4 > n_values) t4 <- t2
    coeff <- index^slope
    a[t2] <- coeff
    a[t4] <- coeff
  }
  a[1] <- 0
  a
}

#' @keywords internal
.pmodel_fractional_integration <- function(y, slope) {
  n <- length(y)
  fourier_coeff <- .pmodel_fractal_spectrum(n, slope / 2)
  mean_val <- mean(y)
  std_y <- stats::sd(y)

  # R's fft(inverse = TRUE) is unnormalised, unlike MATLAB's ifft(); divide
  # by n to match MATLAB's convention exactly.
  xf <- fft(y - mean_val, inverse = TRUE) / n
  phase <- Arg(xf)
  xf2 <- fourier_coeff * exp(1i * phase)
  x <- Re(fft(xf2))
  x <- x * std_y / stats::sd(x)
  x + mean_val
}
