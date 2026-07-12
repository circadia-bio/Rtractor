test_that("sample_entropy() matches the reference mse.c binary on synthetic data (to displayed precision)", {
  set.seed(99)
  x <- rnorm(2000)
  # Reference (mse.c, m=2, r=0.15) at scale 1 (no coarse-graining), run
  # directly on this exact R-generated sequence (set.seed(99); rnorm(2000)):
  # 1  2.504  (tool's own displayed precision, 3 decimals)
  res <- sample_entropy(x, m = 2, r = 0.15)
  expect_equal(res, 2.504, tolerance = 5e-4)
})

test_that("sample_entropy() decreases (roughly) for more regular signals", {
  set.seed(1)
  x_noise <- rnorm(1000)
  x_smooth <- sin(seq(0, 40 * pi, length.out = 1000))
  expect_lt(sample_entropy(x_smooth), sample_entropy(x_noise))
})

test_that("sample_entropy() validates inputs", {
  expect_error(sample_entropy("not numeric"), "must be numeric")
  expect_error(sample_entropy(rnorm(100), m = 0), "m")
  expect_error(sample_entropy(rnorm(100), r = 0), "r")
  expect_error(sample_entropy(rnorm(100), r = -0.1), "r")
})

test_that("sample_entropy() returns NA for too-short series", {
  expect_true(is.na(sample_entropy(c(1, 2), m = 2)))
})
