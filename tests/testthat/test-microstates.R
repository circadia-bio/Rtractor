test_that("recurrence_microstate_entropy() returns valid output structure", {
  set.seed(1)
  x <- (sin(seq(0, 20 * pi, length.out = 300)) + 1) / 2
  res <- recurrence_microstate_entropy(x, n_samples = 2000, frac = 5, frac2 = 5, seed = 1)

  expect_named(res, c("microstate_probs", "entropy_max", "eps_max"))
  expect_length(res$microstate_probs, 2^(3*3))
  expect_equal(sum(res$microstate_probs), 1, tolerance = 1e-8)
  expect_gt(res$entropy_max, 0)
  expect_gt(res$eps_max, 0)
  expect_lt(res$eps_max, 0.5)
})

test_that("recurrence_microstate_entropy() is reproducible with a seed", {
  x <- (sin(seq(0, 20 * pi, length.out = 300)) + 1) / 2
  res1 <- recurrence_microstate_entropy(x, n_samples = 2000, frac = 5, frac2 = 5, seed = 42)
  res2 <- recurrence_microstate_entropy(x, n_samples = 2000, frac = 5, frac2 = 5, seed = 42)

  expect_equal(res1$entropy_max, res2$entropy_max)
  expect_equal(res1$eps_max, res2$eps_max)
})

test_that("recurrence_microstate_entropy() gives lower entropy for a constant signal", {
  # A constant signal has only one possible microstate (all-recurrent),
  # so entropy should collapse toward 0 regardless of threshold.
  x_const <- rep(0.5, 300)
  x_noisy <- runif(300)

  res_const <- recurrence_microstate_entropy(x_const, n_samples = 2000, frac = 5, frac2 = 5, seed = 1)
  res_noisy <- recurrence_microstate_entropy(x_noisy, n_samples = 2000, frac = 5, frac2 = 5, seed = 1)

  expect_lt(res_const$entropy_max, res_noisy$entropy_max)
})

test_that("recurrence_microstate_entropy() validates inputs", {
  expect_error(recurrence_microstate_entropy("not numeric"), "must be numeric")
  expect_error(recurrence_microstate_entropy(rnorm(2), block = 3), "longer than")
})
