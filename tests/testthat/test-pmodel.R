test_that("pmodel() at p = 0.5 gives an exact constant vector of 1s", {
  # sign * (1 - 2*0.5) = sign * 0 = 0 always, regardless of the random
  # sign draws -- this holds for ANY RNG, making it a strong,
  # RNG-independent correctness check for the cascade recursion itself.
  y <- pmodel(1024, p = 0.5)
  expect_true(all(y == 1))

  # Also true with fresh randomness at every level (no reseeding)
  y2 <- pmodel(1024, p = 0.5, seed = NULL)
  expect_true(all(y2 == 1))
})

test_that("pmodel() conserves total mass at every cascade level for any p", {
  # child1 + child2 = 2*y always, regardless of the random sign -- so
  # sum(y) must double at every level, for any p. Check on the full
  # (untruncated) cascade output at a power-of-two length.
  for (test_p in c(0.1, 0.3, 0.375, 0.45)) {
    n <- 2048L
    no_orders <- ceiling(log2(n))
    y <- 1
    for (i in seq_len(no_orders)) y <- Rtractor:::.pmodel_next_step(y, test_p, seed = 42)
    expect_equal(sum(y), 2^no_orders, tolerance = 1e-8)
  }
})

test_that("pmodel() returns a numeric vector of the requested length", {
  y <- pmodel(1000, p = 0.3)
  expect_length(y, 1000)
  expect_true(all(is.finite(y)))
  expect_true(all(y > 0))
})

test_that("pmodel() is more peaked (wider range) for p far from 0.5", {
  y_peaked <- pmodel(2048, p = 0.1, seed = 1)
  y_calm   <- pmodel(2048, p = 0.45, seed = 1)
  expect_gt(diff(range(y_peaked)), diff(range(y_calm)))
})

test_that("pmodel() with a slope performs fractional integration", {
  res <- pmodel(1024, p = 0.3, slope = -2)
  expect_named(res, c("x", "y"))
  expect_length(res$x, 1024)
  expect_length(res$y, 1024)
  expect_true(all(is.finite(res$x)))
  # Rescaling preserves mean and sd of the original p-model series
  expect_equal(mean(res$x), mean(res$y), tolerance = 1e-8)
  expect_equal(stats::sd(res$x), stats::sd(res$y), tolerance = 1e-8)
})

test_that("pmodel() validates inputs", {
  expect_error(pmodel(p = 0), "p")
  expect_error(pmodel(p = 1), "p")
  expect_error(pmodel(p = -0.1), "p")
  expect_error(pmodel(n_values = 0), "n_values")
})

test_that("pmodel() reseed-per-level default is reproducible", {
  y1 <- pmodel(512, p = 0.3, seed = 7)
  y2 <- pmodel(512, p = 0.3, seed = 7)
  expect_identical(y1, y2)
})
