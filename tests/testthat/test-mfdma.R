test_that("mfdma() returns expected structure and monofractal-ish scaling for white noise", {
  set.seed(1)
  x <- rnorm(4000)
  res <- mfdma(x, n_min = 10, n_max = 400, n_scales = 20)

  expect_named(res, c("n", "Fq", "tau", "alpha", "f", "q"))
  expect_true(all(diff(res$n) > 0))
  expect_equal(nrow(res$Fq), length(res$n))
  expect_equal(ncol(res$Fq), length(seq(-4, 4, by = 0.1)))
  # White noise is close to monofractal: alpha should cluster near 0.5
  expect_true(all(res$alpha > 0.3 & res$alpha < 0.8))
  # f(alpha) should peak near 1 for a monofractal-ish signal
  expect_lt(abs(max(res$f) - 1), 0.2)
})

test_that("mfdma() validates inputs", {
  expect_error(mfdma("not numeric"), "must be numeric")
  expect_error(mfdma(rnorm(100), n_min = 1), "n_min")
  expect_error(mfdma(rnorm(100), n_min = 20, n_max = 10), "n_max")
  # Use a longer series so the default n_max doesn't collide with n_min
  # (length 100 with defaults gives n_max = 10 = n_min, which would
  # trigger the n_max check before theta is ever validated).
  expect_error(mfdma(rnorm(2000), theta = 1.5), "theta")
})

test_that("mfdma() respects theta argument without erroring", {
  set.seed(2)
  x <- rnorm(2000)
  res_backward <- mfdma(x, n_min = 10, n_max = 200, n_scales = 15, theta = 0)
  res_centered <- mfdma(x, n_min = 10, n_max = 200, n_scales = 15, theta = 0.5)
  res_forward  <- mfdma(x, n_min = 10, n_max = 200, n_scales = 15, theta = 1)

  expect_true(all(is.finite(res_backward$tau)))
  expect_true(all(is.finite(res_centered$tau)))
  expect_true(all(is.finite(res_forward$tau)))
})

test_that("mfdma() correctly detects increasing multifractality in known p-model ground truth", {
  # Round-trip check against pmodel(): the p-model's multifractality is
  # controlled directly by its `p` parameter (p near 0.5 -> near-
  # monofractal; p far from 0.5 -> strongly multifractal), so mfdma()'s
  # estimated spectrum width should increase monotonically as p moves
  # away from 0.5. Note: mfdma() expects the raw p-model output directly
  # (it does its own internal integration) -- do NOT log-difference it
  # first, which destroys the multifractal structure before mfdma() ever
  # sees it.
  y_mono   <- pmodel(8192, p = 0.48, seed = 1)
  y_multi  <- pmodel(8192, p = 0.15, seed = 1)
  y_strong <- pmodel(8192, p = 0.05, seed = 1)

  res_mono   <- mfdma(y_mono,   n_min = 10, n_max = 400, n_scales = 25)
  res_multi  <- mfdma(y_multi,  n_min = 10, n_max = 400, n_scales = 25)
  res_strong <- mfdma(y_strong, n_min = 10, n_max = 400, n_scales = 25)

  width_mono   <- diff(range(res_mono$alpha))
  width_multi  <- diff(range(res_multi$alpha))
  width_strong <- diff(range(res_strong$alpha))

  # Near p = 0.5, the spectrum should collapse to near-zero width
  expect_lt(width_mono, 0.1)

  # Width must increase monotonically as p moves away from 0.5
  expect_lt(width_mono, width_multi)
  expect_lt(width_multi, width_strong)
})
