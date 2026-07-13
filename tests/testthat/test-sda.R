test_that(".sda_scales() forces +1 steps below force_below, then grows geometrically", {
  scales <- Rtractor:::.sda_scales(scale_min = 5L, scale_max = 200L, growth = 1.1, force_below = 20L)

  # Below the point where 1.1*k first exceeds k+1's forced ceiling (k=10),
  # steps are forced to +1: 5, 6, 7, 8, 9, 10.
  expect_equal(scales[1:6], 5:10)

  # Strictly increasing, and geometric-region ratios hover near `growth`.
  expect_true(all(diff(scales) > 0))
  tail_ratios <- scales[8:length(scales)] / scales[7:(length(scales) - 1)]
  expect_true(all(tail_ratios > 1 & tail_ratios < 1.25))
})

test_that("sda() distinguishes an uncorrelated series from an integrated (Brownian-like) one", {
  set.seed(1)
  white_noise <- rnorm(3000)
  brownian <- cumsum(rnorm(3000))

  res_white <- sda(white_noise, fit_max = 300)
  res_brown <- sda(brownian, fit_max = 300)

  expect_named(res_white, c("k", "F", "H", "se", "r_squared", "intercept", "n_fit"))

  # Windowed RMS-from-local-mean of raw white noise doesn't grow with
  # window size (H near 0); of an integrated/Brownian-like series it does
  # (H well above the white-noise case, in the classical ~0.5 neighbourhood
  # for a random walk).
  expect_lt(res_white$H, 0.2)
  expect_gt(res_brown$H, res_white$H)
  expect_gt(res_brown$H, 0.3)
  expect_lt(res_brown$H, 0.7)
})

test_that("sda() includes the series' first sample (bug fix vs the original mrug tool)", {
  # A single, large first sample should visibly affect the smallest-scale
  # fluctuation value if (and only if) index 0 is included in the windowing.
  x <- c(1000, rnorm(500))
  res <- sda(x, scale_min = 5L, fit_max = 100)
  expect_true(is.finite(res$F[1]))
  expect_gt(res$F[1], 50)  # the outlier dominates any window it appears in
})

test_that("sda() validates inputs", {
  expect_error(sda("not numeric"), "must be numeric")
  expect_error(sda(rnorm(100), growth = 1), "growth")
  expect_error(sda(rnorm(100), scale_min = 20, scale_max = 5), "scale_max")
  expect_error(sda(rnorm(50), fit_min = 40, fit_max = 45), "Fewer than 2 valid window sizes")
})
