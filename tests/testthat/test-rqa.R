test_that("recurrence_matrix() achieves the requested recurrence rate", {
  set.seed(1)
  x <- sin(seq(0, 40 * pi, length.out = 500)) + rnorm(500, sd = 0.05)
  rec <- recurrence_matrix(x, m = 3, tau = 5, rr = 0.1)

  expect_named(rec, c("matrix", "distance", "radius", "rr", "m", "tau", "theiler_window", "norm"))
  expect_equal(dim(rec$matrix), dim(rec$distance))
  # Empirical-quantile threshold hits the target rr closely (exact modulo
  # ties in the distance distribution).
  expect_equal(rec$rr, 0.1, tolerance = 0.02)
})

test_that("recurrence_matrix() validates inputs", {
  expect_error(recurrence_matrix("not numeric"), "must be numeric")
  x <- rnorm(50)
  expect_error(recurrence_matrix(x, radius = 1, rr = 0.1), "only one of")
  expect_error(recurrence_matrix(x, rr = 1.5), "must be in \\(0, 1\\)")
})

test_that("rqa_measures() gives higher determinism for a periodic signal than for noise", {
  set.seed(1)
  t <- seq(0, 40 * pi, length.out = 800)
  periodic <- sin(t)
  noise <- rnorm(800)

  rec_p <- recurrence_matrix(periodic, m = 3, tau = 10, rr = 0.15)
  rec_n <- recurrence_matrix(noise, m = 3, tau = 10, rr = 0.15)

  m_p <- rqa_measures(rec_p)
  m_n <- rqa_measures(rec_n)

  expect_named(m_p, c("RR", "DET", "L_mean", "L_max", "ENTR", "LAM", "TT", "lmin", "vmin"))
  # Both are thresholded to ~the same recurrence rate, so DET is the fair
  # comparison: a periodic signal's recurrences fall on long diagonal
  # lines, an uncorrelated signal's don't.
  expect_gt(m_p$DET, m_n$DET)
  expect_gt(m_p$L_max, m_n$L_max)
})

test_that("rqa_measures() accepts a plain logical matrix", {
  set.seed(1)
  x <- rnorm(100)
  rec <- recurrence_matrix(x, rr = 0.1)
  res <- rqa_measures(rec$matrix)
  expect_true(is.list(res))
  expect_true(is.finite(res$RR))
})

test_that("rqa_measures() validates inputs", {
  expect_error(rqa_measures(list(foo = 1)), "recurrence_matrix\\(\\) result")
  expect_error(rqa_measures("not a matrix"), "recurrence_matrix\\(\\) result")
})
