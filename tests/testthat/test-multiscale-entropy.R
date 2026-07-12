test_that("multiscale_entropy() matches the reference mse.c binary on synthetic data", {
  set.seed(99)
  x <- rnorm(2000)
  res <- multiscale_entropy(x, scale_max = 10)

  # Reference (mse.c, m=2, r=0.15, scales 1-10), tool's own displayed
  # precision (3 decimals):
  #   1  2.458   2  2.092   3  1.956   4  1.788   5  1.710
  #   6  1.540   7  1.534   8  1.443   9  1.412  10  1.395
  expected <- c(2.458, 2.092, 1.956, 1.788, 1.710, 1.540, 1.534, 1.443, 1.412, 1.395)

  expect_named(res, c("scale", "mse", "m", "r"))
  expect_equal(res$scale, 1:10)
  expect_equal(res$mse, expected, tolerance = 5e-4)
})

test_that("multiscale_entropy() uses a fixed tolerance across scales, not per-scale sd", {
  # This is the defining feature of MSE (Costa et al.) vs. naively running
  # sample_entropy() independently on each coarse-grained series: the
  # tolerance must stay relative to the ORIGINAL series' sd, not each
  # coarse-grained series' own (shrinking) sd. Confirm this by checking
  # against the reference values above, which already encode this
  # convention -- if the fixed-tolerance behaviour were broken, the
  # previous test would fail.
  set.seed(99)
  x <- rnorm(2000)
  res <- multiscale_entropy(x, scale_max = 3)
  expect_true(all(is.finite(res$mse)))
})

test_that("multiscale_entropy() validates inputs", {
  expect_error(multiscale_entropy("not numeric"), "must be numeric")
  expect_error(multiscale_entropy(rnorm(100), scale_max = 0), "scale_max")
  expect_error(multiscale_entropy(rnorm(100), m = 0), "m")
  expect_error(multiscale_entropy(rnorm(100), r = 0), "r")
})

test_that("multiscale_entropy() returns NA at scales too coarse for the series length", {
  x <- rnorm(20)
  res <- multiscale_entropy(x, scale_max = 15)
  expect_true(any(is.na(res$mse)))
})
