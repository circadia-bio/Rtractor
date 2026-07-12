test_that("higuchi_fd() returns HFD near 2 for white noise", {
  set.seed(1)
  x <- rnorm(1000)
  res <- higuchi_fd(x, k_max = 10)

  expect_named(res, c("k", "L", "hfd"))
  expect_length(res$L, 10)
  expect_true(all(diff(res$L) < 0))  # curve length decreases with k
  # White noise is space-filling: HFD ~ 2
  expect_gt(res$hfd, 1.85)
  expect_lt(res$hfd, 2.05)
})

test_that("higuchi_fd() returns HFD near 1 for a smooth signal", {
  x <- sin(seq(0, 20 * pi, length.out = 2000))
  res <- higuchi_fd(x, k_max = 10)

  # Smooth periodic signal: HFD close to 1 (line-like)
  expect_lt(res$hfd, 1.3)
})

test_that("higuchi_fd() validates inputs", {
  expect_error(higuchi_fd("not numeric"), "must be numeric")
  expect_error(higuchi_fd(rnorm(100), k_max = 1), "k_max")
})
