test_that("perm_entropy() matches antropy's normalized value on a fixed reference sequence", {
  # Fixed literal values (not RNG-derived) so the reference is reproducible
  # across R and the antropy Python library used to validate this port.
  x <- c(0.12, -0.45, 1.3, 0.67, -1.1, 0.05, 0.88, -0.33, 1.5, -0.6,
         0.2, 0.9, -0.15, 1.1, -0.8, 0.4, -1.2, 0.75, -0.05, 1.05)

  # antropy 0.2.2 reference (normalize = TRUE). Note: the *unnormalized*
  # raw entropy is NOT directly comparable between antropy and this port
  # -- antropy computes it in log base 2 (bits) internally, this port
  # uses natural log (nats). The normalized ratio (h / log(order!)) is
  # base-independent, which is why only that value is checked here.
  expect_equal(perm_entropy(x, normalize = TRUE), 0.9259640358613577, tolerance = 1e-8)
})

test_that("perm_entropy() normalize argument behaves as documented", {
  set.seed(1)
  x <- rnorm(500)
  normalized   <- perm_entropy(x, normalize = TRUE)
  unnormalized <- perm_entropy(x, normalize = FALSE)
  expect_true(normalized >= 0 && normalized <= 1)
  expect_equal(normalized, unnormalized / sum(log(1:3)))
})

test_that("perm_entropy() validates inputs", {
  expect_error(perm_entropy("not numeric"), "must be numeric")
  expect_error(perm_entropy(rnorm(100), order = 1), "order")
  expect_error(perm_entropy(rnorm(100), delay = 0), "delay")
})
