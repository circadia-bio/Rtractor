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
