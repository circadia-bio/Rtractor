test_that("dfa() recovers expected scaling for white noise", {
  set.seed(1)
  x <- rnorm(3000)
  res <- dfa(x)

  expect_named(res, c("n", "F", "alpha"))
  expect_true(all(res$F > 0))
  expect_true(all(diff(res$n) > 0))
  # White noise: alpha ~ 0.5
  expect_gt(res$alpha, 0.35)
  expect_lt(res$alpha, 0.65)
})

test_that("dfa() recovers expected scaling for a random walk", {
  set.seed(1)
  x <- rnorm(3000)
  # A random walk fed through the standard DFA pipeline (default
  # integrate = TRUE) amounts to double integration -- the classic
  # textbook benchmark that gives alpha ~ 1.5.
  walk <- cumsum(x)
  res <- dfa(walk)

  expect_gt(res$alpha, 1.3)
  expect_lt(res$alpha, 1.7)
})

test_that("dfa() with integrate = FALSE analyses the profile directly", {
  set.seed(1)
  x <- rnorm(3000)
  walk <- cumsum(x)
  # Skipping the internal integration and analysing the random walk
  # trajectory directly (as its own profile) recovers its own Hurst
  # exponent, ~0.5 for a standard random walk -- verified against the
  # reference dfa.c binary run with its -i flag.
  res <- dfa(walk, integrate = FALSE)

  expect_gt(res$alpha, 0.35)
  expect_lt(res$alpha, 0.65)
})

test_that("dfa() validates inputs", {
  expect_error(dfa("not numeric"), "must be numeric")
  expect_error(dfa(rnorm(100), order = 0), "order")
})

test_that("dfa() respects order/min_box/max_box arguments", {
  set.seed(2)
  x <- rnorm(2000)
  res <- dfa(x, order = 2, min_box = 8, max_box = 200)
  expect_true(min(res$n) >= 8)
  expect_true(max(res$n) <= 200)
})
