test_that("chhabra_jensen() returns expected structure", {
  set.seed(1)
  x <- abs(rnorm(1024)) + 0.01
  res <- chhabra_jensen(x, scales = 1:6)

  expect_named(res, c(
    "alpha", "falpha", "Dq",
    "r_squared_alpha", "r_squared_falpha", "r_squared_Dq",
    "q", "mu_scale", "Ma", "Mf", "Md"
  ))
  expect_length(res$alpha, length(res$q))
  expect_length(res$Dq, length(res$q))
  expect_true(all(res$r_squared_Dq >= -1e-8 & res$r_squared_Dq <= 1 + 1e-8))
})

test_that("chhabra_jensen() validates inputs", {
  expect_error(chhabra_jensen("not numeric"), "must be numeric")
  expect_error(chhabra_jensen(c(-1, 2, 3)), "strictly positive")
  expect_error(chhabra_jensen(abs(rnorm(100)) + 0.1, scales = 1:5), "not divisible")
})

test_that("chhabra_jensen() warns on q values strictly between 0 and 1", {
  set.seed(1)
  x <- abs(rnorm(1024)) + 0.01
  expect_warning(
    chhabra_jensen(x, q_values = c(-1, 0.5, 2), scales = 1:6),
    "numerically unstable"
  )
})

test_that("chhabra_jensen() gives Dq close to 1 for near-uniform positive noise", {
  set.seed(5)
  x <- abs(rnorm(1024)) + 0.01
  res <- chhabra_jensen(x, q_values = c(-2, 2), scales = 1:6)
  # i.i.d. positive weights approach a generalized dimension close to 1
  expect_true(all(abs(res$Dq - 1) < 0.3))
})

test_that("chhabra_jensen() correctly detects increasing multifractality in known p-model ground truth", {
  # Round-trip check against pmodel(): the p-model's multifractality is
  # controlled directly by its `p` parameter (p near 0.5 -> near-
  # monofractal; p far from 0.5 -> strongly multifractal), so
  # chhabra_jensen()'s estimated spectrum width should increase
  # monotonically as p moves away from 0.5. p-model output is already
  # strictly positive, so no sigmoid transform is needed here.
  y_mono   <- pmodel(8192, p = 0.48, seed = 1)
  y_multi  <- pmodel(8192, p = 0.15, seed = 1)
  y_strong <- pmodel(8192, p = 0.05, seed = 1)

  res_mono   <- chhabra_jensen(y_mono,   scales = 1:11)
  res_multi  <- chhabra_jensen(y_multi,  scales = 1:11)
  res_strong <- chhabra_jensen(y_strong, scales = 1:11)

  width_mono   <- diff(range(res_mono$alpha))
  width_multi  <- diff(range(res_multi$alpha))
  width_strong <- diff(range(res_strong$alpha))

  # Near p = 0.5, the spectrum should collapse to near-zero width
  expect_lt(width_mono, 0.1)

  # Width must increase monotonically as p moves away from 0.5
  expect_lt(width_mono, width_multi)
  expect_lt(width_multi, width_strong)

  # Theoretical constraint: the spectrum's peak height f(alpha) should be
  # close to 1 (the fractal dimension of the 1D support) in all cases
  expect_equal(max(res_mono$falpha), 1, tolerance = 0.01)
  expect_equal(max(res_multi$falpha), 1, tolerance = 0.01)
  expect_equal(max(res_strong$falpha), 1, tolerance = 0.02)
})
