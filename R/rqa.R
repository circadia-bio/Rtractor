# Recurrence quantification analysis (RQA)
#
# Family: rqa
# Wraps established reference implementations for recurrence plots and their
# derived quantification measures. Signal-agnostic (any numeric vector,
# regardless of source — EEG, actigraphy, BOLD, HRV, etc.).

#' Recurrence Microstates Entropy (maximum-entropy threshold search)
#'
#' Estimates the complexity of a (windowed) time series using the entropy
#' of recurrence microstates (Corso et al. 2018), searching over the
#' vicinity threshold `epsilon` to find the value that *maximises* the
#' Shannon entropy of the microstate distribution — a parameter-free
#' alternative to fixing `epsilon` by hand, related to the optimal-threshold
#' method of Prado et al. (2018).
#'
#' Ported from `circadia-bio/maxEntropy` (`Entropy.jl`, MIT licensed), the
#' script used to produce the entropy estimates in Ferre et al. (2024,
#' PLOS ONE). The deterministic core loop was validated against a Python
#' transliteration of the original on synthetic test data (exact match on
#' `S_max` and `eps_max` given identical sampled index pairs). Only the
#' reusable `Max_Entropy()` routine was ported — the original script's
#' `main()` (file I/O over named subject directories) is study-specific
#' scaffolding, not part of the general-purpose method.
#'
#' @param x Numeric vector. The (already windowed) time series segment to
#'   analyse. In the original study, `x` was min-max normalised to `[0, 1]`
#'   per window before calling this function — do the same if you want
#'   comparable `eps_max` values across windows/subjects.
#' @param block Integer. Side length of the square recurrence microstate
#'   block; the microstate space has `2^(block^2)` possible patterns.
#'   Default `3` (matching the reference implementation).
#' @param n_samples Integer. Number of random point-pairs sampled to
#'   estimate the recurrence microstate distribution. Default `18000`.
#' @param eps_min,eps_max Numeric. Search range for the vicinity threshold.
#'   Defaults match the reference implementation.
#' @param frac,frac2 Integer. Number of steps in the coarse search pass and
#'   number of refinement iterations, respectively. Defaults match the
#'   reference implementation.
#' @param seed Optional integer. Seeds the random pair sampling for
#'   reproducibility. The original Julia script reseeds from system entropy
#'   on every run (`Random.seed!()` with no argument) and is therefore not
#'   reproducible run-to-run by design; pass a `seed` here if you need
#'   reproducible results, or leave `NULL` to match the original's
#'   fresh-randomness-every-call behaviour.
#'
#' @return A list with:
#'   \item{microstate_probs}{Probability of each of the `2^(block^2)`
#'     microstates at the entropy-maximising threshold.}
#'   \item{entropy_max}{The maximum Shannon entropy found (`S_max`).}
#'   \item{eps_max}{The vicinity threshold at which `entropy_max` was
#'     achieved.}
#'
#' @references
#' Corso G, Prado TL, dos Santos Lima GZ, Kurths J, Lopes SR. Quantifying
#' entropy using recurrence matrix microstates. Chaos 2018;28(8):083108.
#'
#' Prado TL, dos Santos Lima GZ, Lobao-Soares B, do Nascimento GC, Corso G,
#' Fontenele-Araujo J, Kurths J, Lopes SR. Optimizing the detection of
#' nonstationary signals by using recurrence analysis. Chaos
#' 2018;28(8):085703.
#'
#' Ferre IBS, Corso G, dos Santos Lima GZ, Lopes SR, Leocadio-Miguel MA,
#' Franca LGS, de Lima Prado T, Araujo JF. Cycling reduces the entropy of
#' neuronal activity in the human adult cortex. PLOS ONE 2024;19(10):e0298703.
#'
#' @examples
#' set.seed(1)
#' x <- (sin(seq(0, 20 * pi, length.out = 300)) + 1) / 2
#' recurrence_microstate_entropy(x, seed = 1)
#'
#' @export
recurrence_microstate_entropy <- function(x, block = 3L, n_samples = 18000L,
                                           eps_min = 1e-9, eps_max = 0.499999999,
                                           frac = 10L, frac2 = 10L, seed = NULL) {
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  window_size <- length(x)
  block <- as.integer(block)
  if (window_size <= block) {
    stop("`x` must be longer than `block`.", call. = FALSE)
  }

  if (!is.null(seed)) {
    if (exists(".Random.seed", envir = .GlobalEnv)) {
      old_seed <- .Random.seed
      on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
    }
    set.seed(seed)
  }

  hi <- window_size - block
  x_idx <- as.integer(round(stats::runif(n_samples, min = 0, max = hi)))
  y_idx <- as.integer(round(stats::runif(n_samples, min = 0, max = hi)))

  res <- recurrence_microstate_entropy_cpp(
    as.double(x), x_idx, y_idx, block,
    as.double(eps_min), as.double(eps_max), as.integer(frac), as.integer(frac2)
  )

  list(
    microstate_probs = res$microstate_probs,
    entropy_max = res$entropy_max,
    eps_max = res$eps_max
  )
}

# ---- Planned, not yet implemented -----------------------------------------
#
#   - plot_recurrence()      Recurrence plot visualisation (Rtractor palette)
#
# Native counterpart(s) expected in: src/rqa.cpp (for the above)

#' Recurrence Matrix
#'
#' Builds a binary recurrence matrix from a time series: `R[i, j] = 1` if
#' the (optionally embedded) trajectory points `i` and `j` are within a
#' distance `radius` of each other, `0` otherwise (Eckmann, Kamphorst &
#' Ruelle 1987). No reference implementation was specified for this port
#' (no preference given) -- implemented directly from the standard
#' formulas in Marwan et al. (2007). See `inst/COPYRIGHTS`.
#'
#' @param x Numeric vector. The time series to analyse.
#' @param m Integer. Embedding dimension. Default `1` (no embedding --
#'   recurrence of the raw scalar series). Use [estimate_embed_dim()] to
#'   choose this for a proper phase-space reconstruction.
#' @param tau Integer. Time delay used if `m > 1`. Default `1`. Use
#'   [estimate_delay()] to choose this.
#' @param radius Numeric or `NULL`. Fixed vicinity threshold. Exactly one
#'   of `radius`/`rr` should be given (or neither, see below).
#' @param rr Numeric in `(0, 1)` or `NULL`. Target recurrence rate: `radius`
#'   is set to the empirical quantile of pairwise distances (outside the
#'   Theiler band) that yields this recurrence rate. If both `radius` and
#'   `rr` are `NULL`, `rr = 0.1` is used (a common default for comparable
#'   RQA measures across signals).
#' @param norm Distance norm: `"euclidean"` (default), `"maximum"`
#'   (Chebyshev), or `"manhattan"`.
#' @param theiler_window Integer. Points within this many samples of the
#'   main diagonal are excluded from `rr`/distance-quantile calculations
#'   (they are trivially close in time, not dynamically informative).
#'   Default `1` (excludes just the main diagonal `i == j`).
#'
#' @return A list with:
#'   \item{matrix}{Logical `N x N` recurrence matrix.}
#'   \item{distance}{Numeric `N x N` distance matrix (before thresholding).}
#'   \item{radius}{The threshold used.}
#'   \item{rr}{The achieved recurrence rate (Theiler band excluded).}
#'   \item{m, tau, theiler_window, norm}{Echoed parameters, for
#'     [rqa_measures()] and downstream use.}
#'
#' This is O(N^2) in both time and memory (a dense `N x N` matrix); a
#' warning is issued above 5000 points.
#'
#' @references
#' Eckmann JP, Kamphorst SO, Ruelle D. Recurrence plots of dynamical
#' systems. Europhys Lett 1987;4(9):973-977.
#'
#' Marwan N, Romano MC, Thiel M, Kurths J. Recurrence plots for the
#' analysis of complex systems. Phys Rep 2007;438(5-6):237-329.
#'
#' @examples
#' set.seed(1)
#' x <- sin(seq(0, 40, length.out = 400)) + rnorm(400, sd = 0.05)
#' rec <- recurrence_matrix(x, m = 3, tau = 5, rr = 0.1)
#' rec$rr
#'
#' @export
recurrence_matrix <- function(x, m = 1L, tau = 1L, radius = NULL, rr = NULL,
                               norm = c("euclidean", "maximum", "manhattan"),
                               theiler_window = 1L) {
  norm <- match.arg(norm)
  if (!is.numeric(x)) stop("`x` must be numeric.", call. = FALSE)
  if (!is.null(radius) && !is.null(rr)) {
    stop("Specify only one of `radius` or `rr`, not both.", call. = FALSE)
  }

  traj <- if (m > 1L) embed_time_series(x, m, tau) else matrix(as.double(x), ncol = 1L)
  N <- nrow(traj)
  if (N > 5000L) {
    warning(
      "Recurrence matrix is O(N^2) in memory/time; ", N,
      " points may be slow/large. Consider subsampling `x` first.",
      call. = FALSE
    )
  }

  norm_code <- switch(norm, euclidean = 1L, maximum = 2L, manhattan = 3L)
  D <- distance_matrix_cpp(traj, norm_code)

  theiler_mask <- abs(row(D) - col(D)) > theiler_window
  if (is.null(radius) && is.null(rr)) rr <- 0.1
  if (!is.null(rr)) {
    if (rr <= 0 || rr >= 1) stop("`rr` must be in (0, 1).", call. = FALSE)
    radius <- stats::quantile(D[theiler_mask], probs = rr, names = FALSE)
  }

  R <- D <= radius
  rr_actual <- mean(R[theiler_mask])

  list(
    matrix = R, distance = D, radius = radius, rr = rr_actual,
    m = m, tau = tau, theiler_window = theiler_window, norm = norm
  )
}

#' Recurrence Quantification Analysis (RQA) Measures
#'
#' Computes the standard RQA measures from a recurrence matrix (Zbilut &
#' Webber 1992; Marwan et al. 2007): recurrence rate, determinism,
#' diagonal-line statistics, laminarity, and trapping time. No reference
#' implementation was specified for this port (no preference given) --
#' implemented directly from the standard formulas. See `inst/COPYRIGHTS`.
#'
#' @param rec Either the list returned by [recurrence_matrix()], or a plain
#'   logical `N x N` matrix (in which case `theiler_window = 1` is assumed).
#' @param lmin Integer. Minimum diagonal line length counted toward `DET`,
#'   `L_mean`, `L_max`, and `ENTR`. Default `2`.
#' @param vmin Integer. Minimum vertical line length counted toward `LAM`
#'   and `TT`. Default `2`.
#'
#' @return A list with:
#'   \item{RR}{Recurrence rate: fraction of recurrent points (Theiler band
#'     excluded).}
#'   \item{DET}{Determinism: fraction of recurrent points forming diagonal
#'     lines of length `>= lmin`.}
#'   \item{L_mean}{Mean diagonal line length (lines `>= lmin`).}
#'   \item{L_max}{Longest diagonal line length.}
#'   \item{ENTR}{Shannon entropy (natural log) of the diagonal line length
#'     distribution (lines `>= lmin`).}
#'   \item{LAM}{Laminarity: fraction of recurrent points forming vertical
#'     lines of length `>= vmin`.}
#'   \item{TT}{Trapping time: mean vertical line length (lines `>= vmin`).}
#'   \item{lmin, vmin}{Echoed parameters.}
#'
#' @references
#' Zbilut JP, Webber CL Jr. Embeddings and delays as derived from
#' quantification of recurrence plots. Phys Lett A 1992;171(3-4):199-203.
#'
#' Marwan N, Romano MC, Thiel M, Kurths J. Recurrence plots for the
#' analysis of complex systems. Phys Rep 2007;438(5-6):237-329.
#'
#' @examples
#' set.seed(1)
#' x <- sin(seq(0, 40, length.out = 400)) + rnorm(400, sd = 0.05)
#' rec <- recurrence_matrix(x, m = 3, tau = 5, rr = 0.1)
#' rqa_measures(rec)
#'
#' @export
rqa_measures <- function(rec, lmin = 2L, vmin = 2L) {
  if (is.list(rec) && !is.null(rec$matrix)) {
    R <- rec$matrix
    theiler_window <- rec$theiler_window
  } else if (is.matrix(rec)) {
    R <- rec
    theiler_window <- 1L
  } else {
    stop("`rec` must be a recurrence_matrix() result or a logical matrix.", call. = FALSE)
  }

  raw <- rqa_line_stats_cpp(R, as.integer(theiler_window))
  total <- raw$total_points
  diag_len <- raw$diag_lengths
  vert_len <- raw$vert_lengths
  N <- raw$N

  rr <- if ((N * N - raw$n_excluded) > 0) total / (N * N - raw$n_excluded) else NA_real_

  det_mask <- diag_len >= lmin
  det <- if (total > 0) sum(diag_len[det_mask]) / total else NA_real_
  l_mean <- if (any(det_mask)) mean(diag_len[det_mask]) else NA_real_
  l_max <- if (length(diag_len) > 0) max(diag_len) else 0L
  p_l <- table(diag_len[det_mask])
  p_l <- p_l / sum(p_l)
  entr <- if (length(p_l) > 0) -sum(p_l * log(p_l)) else NA_real_

  lam_mask <- vert_len >= vmin
  lam <- if (total > 0) sum(vert_len[lam_mask]) / total else NA_real_
  tt <- if (any(lam_mask)) mean(vert_len[lam_mask]) else NA_real_

  list(
    RR = rr, DET = det, L_mean = l_mean, L_max = l_max, ENTR = entr,
    LAM = lam, TT = tt, lmin = lmin, vmin = vmin
  )
}
