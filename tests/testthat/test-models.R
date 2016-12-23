context("Models run")

test_that("CMSY runs", {
  x <- cmsy(blue_gren$yr, ct = blue_gren$ct, reps = 2e4)
  expect_equal(class(x$bbmsy), "data.frame")
})

test_that("COMSIR runs", {
  skip_on_cran()

  x <- comsir(blue_gren$yr, ct = blue_gren$ct, nsim = 1e6)
  expect_equal(class(x$bbmsy), "data.frame")

  expect_warning(comsir(blue_gren$yr, ct = blue_gren$ct, nsim = 100)) # too few samples
})

## test_that("SSCOM runs", {
##   skip_on_cran()
##   skip_on_travis()
##   skip_on_appveyor()
##
##   expect_warning({
##     x <- sscom(blue_gren$yr, blue_gren$ct,
##       NburninPrelim = 1e3,
##       NiterPrelim = 2e3,
##       NthinPrelim = 1,
##       NchainsPrelim = 20,
##       NburninJags = 1e3,
##       NiterJags = 2e3,
##       NthinJags = 1,
##       Nchains = 2, return_jags = TRUE)
##   })
##
##   head(x$bbmsy)
##   x$jags
##   expect_equal(class(x$bbmsy), "data.frame")
##   expect_equal(class(x$jags), "rjags")
## })

test_that("mPRM runs", {
  d <- subset(ram_prm_dat, stockid == "BGRDRSE")
  x <- predict_prm(d, ci = TRUE)
  expect_equal(class(x), "data.frame")
})
