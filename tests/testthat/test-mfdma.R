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
  expect_error(mfdma(rnorm(100), theta = 1.5), "theta")
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
