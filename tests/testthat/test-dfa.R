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
  # already-integrated random walk; analyse without re-integrating
  walk <- cumsum(x)
  res <- dfa(walk, integrate = FALSE)

  # Random walk: alpha ~ 1.5
  expect_gt(res$alpha, 1.3)
  expect_lt(res$alpha, 1.7)
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
