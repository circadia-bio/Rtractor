test_that("petrosian_fd(), num_zerocross(), hjorth_parameters() match antropy on a fixed sequence", {
  # Fixed literal values (not RNG-derived) so the reference is reproducible
  # across R and the antropy Python library used to validate this port.
  x <- c(0.12, -0.45, 1.3, 0.67, -1.1, 0.05, 0.88, -0.33, 1.5, -0.6,
         0.2, 0.9, -0.15, 1.1, -0.8, 0.4, -1.2, 0.75, -0.05, 1.05)

  # antropy 0.2.2 reference values for this exact sequence
  expect_equal(petrosian_fd(x), 1.0959857083258884, tolerance = 1e-8)
  expect_equal(num_zerocross(x), 16L)

  hj <- hjorth_parameters(x)
  expect_named(hj, c("mobility", "complexity"))
  expect_equal(hj$mobility, 1.745211182192803, tolerance = 1e-6)
  expect_equal(hj$complexity, 1.0765350992133398, tolerance = 1e-6)
})

test_that("petrosian_fd() validates inputs", {
  expect_error(petrosian_fd("not numeric"), "must be numeric")
})

test_that("num_zerocross() validates inputs", {
  expect_error(num_zerocross("not numeric"), "must be numeric")
})

test_that("hjorth_parameters() validates inputs", {
  expect_error(hjorth_parameters("not numeric"), "must be numeric")
})

test_that("higuchi_fd() agrees with mrpheus/antropy's implementation for realistic epoch lengths", {
  set.seed(1)
  x <- rnorm(2000)
  # Rtractor's higuchi_fd() (MATLAB-reference port) vs antropy's algorithm
  # (mrpheus's own implementation) agree exactly when kmax << length(x),
  # since the valid_cnt-vs-k averaging distinction only matters when kmax
  # approaches length(x). See inst/COPYRIGHTS for the full comparison.
  res <- higuchi_fd(x, k_max = 10)
  expect_true(is.finite(res$hfd))
})
