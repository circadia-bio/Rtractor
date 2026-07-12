test_that("rtractor_palette() returns known palettes", {
  main <- rtractor_palette()
  expect_named(main, c("coral", "cream", "sage", "steel_blue", "ink"))
  expect_length(main, 5)

  core <- rtractor_palette("core")
  expect_length(core, 4)
  expect_false("ink" %in% names(core))
})

test_that("rtractor_palette() respects n and reverse", {
  pal <- rtractor_palette("core", n = 2)
  expect_length(pal, 2)

  pal_rev <- rtractor_palette("core", reverse = TRUE)
  expect_identical(unname(pal_rev), rev(unname(rtractor_palette("core"))))
})

test_that("rtractor_palette() errors on unknown palette or oversized n", {
  expect_error(rtractor_palette("not_a_palette"), "Unknown palette")
  expect_error(rtractor_palette("core", n = 10), "only 4 colours")
})

test_that("rtractor_palettes() returns sizes invisibly", {
  output <- capture.output(sizes <- rtractor_palettes())
  expect_true(length(output) > 0)
  expect_type(sizes, "integer")
  expect_true(all(c("main", "core", "diverging") %in% names(sizes)))
})
