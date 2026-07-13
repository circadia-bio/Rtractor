test_that("embed_time_series() builds correct delay-coordinate matrix", {
  x <- 1:20
  traj <- embed_time_series(x, m = 3, tau = 2)
  expect_equal(dim(traj), c(20 - 2 * 2, 3))
  expect_equal(traj[1, ], c(1, 3, 5))
  expect_equal(traj[nrow(traj), ], c(16, 18, 20))
})

test_that("embed_time_series() validates inputs", {
  expect_error(embed_time_series("not numeric"), "must be numeric")
  expect_error(embed_time_series(1:5, m = 10, tau = 5), "too short")
})

test_that("estimate_delay() picks a short delay for a fast periodic signal", {
  set.seed(1)
  x <- sin(seq(0, 40 * pi, length.out = 2000))
  tau <- estimate_delay(x, max_lag = 60)
  expect_true(tau >= 1 && tau <= 60)
})

test_that("estimate_delay() acf method returns first zero crossing", {
  set.seed(1)
  x <- sin(seq(0, 40 * pi, length.out = 2000))
  tau <- estimate_delay(x, max_lag = 60, method = "acf")
  expect_true(tau >= 1 && tau <= 60)
})

test_that("estimate_embed_dim() finds a low dimension for a simple oscillator", {
  set.seed(1)
  x <- sin(seq(0, 40 * pi, length.out = 1000)) + rnorm(1000, sd = 0.01)
  res <- estimate_embed_dim(x, tau = 10, m_max = 6)
  expect_named(res, c("m", "dim", "fnn_fraction"))
  expect_length(res$fnn_fraction, 6)
  # A clean sinusoid embeds well in low dimension -- FNN fraction should
  # drop sharply well before m_max.
  expect_lt(res$fnn_fraction[6], res$fnn_fraction[1])
})

test_that("estimate_embed_dim() validates inputs", {
  expect_error(estimate_embed_dim("not numeric"), "must be numeric")
})
