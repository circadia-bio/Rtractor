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
#   - recurrence_matrix()    Binary/distance recurrence matrix from a
#                            reconstructed trajectory
#   - rqa_measures()         Determinism, laminarity, recurrence rate,
#                            trapping time, mean/max diagonal line length,
#                            entropy of diagonal line lengths, etc.
#   - plot_recurrence()      Recurrence plot visualisation (Rtractor palette)
#
# Requires phase-space reconstruction first (embedding dimension + time
# delay) — see embed.R.
#
# Native counterpart(s) expected in: src/rqa.cpp (for the above)
