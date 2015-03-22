#' Schaefer production model
#'
#' @param r Intrinsic population growth rate
#' @param k Carrying capacity
#' @param sigR Standard deviation of process noise
#' @param startbt A vector of possible starting biomasses to loop over
#' @param yr Numeric vector of years
#' @param ct Numeric vector of catch
#' @param interyr Interim year within time series for which biomass estimate is
#'   available. Set to yr[2] if no estimates are available.
#' @param interbio_lim A numeric vector of length 2 that gives the lower and
#'   upper biomass limits in the interim year
#' @param prior_log_mean Prior log mean
#' @param prior_log_sd Prior log sd
#'
#' @useDynLib datalimited
#' @importFrom Rcpp sourceCpp
#'
#' @rdname cmsy
#'
#' @examples
#' schaefer_cmsy_r(0.1, 100, sigR = 0.1, startbt = seq(0.2, 0.6, by = 0.05),
#'   yr = 1:10, ct = rnorm(10), prior_log_mean = 0.3, prior_log_sd = 0.1,
#'   interyr = 2, interbio_lim = c(0, 1))
#'
#' # Benchmark the R version vs. the C++ version:
#' r_lim <- c(0.1, 0.4)
#' k_lim <- c(80, 120)
#' N <- 200L
#'
#' r_test <- function() {
#'   r <- runif(N, r_lim[1], r_lim[2])
#'   k <-  exp(runif(N, log(k_lim[1]), log(k_lim[2])))
#'   pars <- cbind(r, k)
#'   apply(pars, 1, function(x) schaefer_cmsy_r(r = x[["r"]], k = x[["k"]],
#'     sigR = 0.1, startbt = seq(0.2, 0.6, by = 0.05),
#'     yr = 1:10, ct = rnorm(10), prior_log_mean = 0.3, prior_log_sd = 0.1,
#'     interyr = 2, interbio_lim = c(0, 1)))
#' }
#'
#' cpp_test <- function() {
#'   schaefer_cmsy(r_lim, k_lim, sigR = 0.1, startbt = seq(0.2, 0.6, by = 0.05),
#'     yr = 1:10, ct = rnorm(10), prior_log_mean = 0.3, prior_log_sd = 0.1,
#'     interyr_index = 2, interbio_lim = c(0, 1), N = N)
#' }
#'
#' op <- microbenchmark::microbenchmark(
#'   r     = r_test(),
#'   cpp   = cpp_test(),
#'   times = 100L)
#' print(op)


# Steps to make this orders of magnitude faster:
# - calculate constants once at beginning (e.g. prior_log_mean - log(2))
# - index things that need to be indexed many times once as a constant (e.g.
#   interbio_lim[1], and [2]
# - generate random variables as a matrix before (rnorm(1,0,sigR), runif(1,0,1))
# - set full vectors before and fill them
# - extract the loopy bits and C++ them

schaefer_cmsy_r <- function(r, k, sigR, startbt, nyr, yr, ct, interyr = 2,
  prior_log_mean, prior_log_sd, interbio_lim = c(0, 1)) {

  nyr <- length(yr)
  bt <- 0
  ell <- 0  # initialize ell
  J <- 0
  for (j in startbt) {
    # print(j)
    # print(ell)
    if (ell == 0) {
      bt[1] <- j * k * exp(rnorm(1, 0, sigR)) # set biomass in first year
      # bt[1] <- j * k * exp(-0.2) # set biomass in first year
      for(i in 1:nyr) { # for all years in the time series
        xt <- rnorm(1,0, sigR)
        # xt <- 0.1
        # calculate biomass as function of previous year's biomass plus net
        # production minus catch:
        bt[i + 1] <- (bt[i] + r * bt[i] * (1 - bt[i]/k) - ct[i]) * exp(xt)
      }
      # Bernoulli likelihood, assign 0 or 1 to each combination of r and k
      ell <- 0
      # New posterior predictive prior on final year biomass
      current.bio.ratio <- bt[nyr + 1]/k
      tmp <- runif(1, 0, 1)
      # tmp <- 0.05
      # print(current.bio.ratio)
      test <- (dlnorm(current.bio.ratio,
        meanlog = prior_log_mean - log(2),
        sdlog = prior_log_sd)) /
        dlnorm(exp(prior_log_mean - log(2)),
          meanlog = prior_log_mean - log(2),
          sdlog = prior_log_sd)
      if (tmp < test &&
          min(bt) > 0 &&
          max(bt) <= k &&
          bt[which(yr == interyr)]/k >= interbio_lim[1] &&
          bt[which(yr == interyr)]/k <= interbio_lim[2]) {
        ell <- 1
      }
      J <- j
    }
  }
  list(ell = ell, J = J)
}
