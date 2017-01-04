context("mPRM")

dat <- format_prm(
  year = 1:10,
  catch = 1:10,
  species_cat = "test")

test_that("format_data returns correct values", {
  expect_equal(dat$years_back, 10:1)
  expect_equal(dat$catch, 1:10)
  expect_equal(dat$max_catch, rep(10, 10))
  expect_equal(dat$scaled_catch, seq(0.1, 1, 0.1))
  expect_equal(dat$mean_scaled_catch, rep(mean(seq(0.1, 1, 0.1)), 10))
  expect_equal(dat$scaled_catch1, c(NA, seq(0.1, 0.9, 0.1)))
  expect_equal(dat$scaled_catch2, c(NA, NA, seq(0.1, 0.8, 0.1)))
  expect_equal(dat$catch_to_rolling_max, rep(1, 10))
  expect_equal(dat$time_to_max, rep(10, 10))
  expect_equal(dat$initial_slope, rep(0.1, 10))
})
