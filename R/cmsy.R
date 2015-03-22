#' Catch-MSY method
#'
#' Apply the catch-MSY method from Martell and Froese (2013).
#'
#' @details
#' \code{cmsy} is a wrapper for \code{shaefer_cmsy} that implements suggested
#' argument values and translates resiliency categories into ranges of intrinsic
#' growth rate. \code{shaefer_cmsy} is also an exported function so it can be
#' used to enter customized ranges of intrinsic growth rate. The
#' \code{shaefer_cmsy} function is written in C++ with the \pkg{Rcpp} package
#' for speed.
#'
#' @param yr Numeric vector of years
#' @param ct Numeric vector of catch
#' @param sig_r Standard deviation of process noise
#' @param bio_step Step size between lower and upper biomass limits in the iterim year
#' @param resilience A character value designating the resilience of the stock
#' @param interyr_index Index of the interim year within time series for which
#'   biomass estimate is available
#' @param interbio A numeric vector that gives the lower and upper biomass
#'   limits in the interim year
#' @param prior_log_mean Prior log mean
#' @param prior_log_sd Prior log sd
#' @param reps Number of repititions to run the sampling
#' @param remove_ell0 Should sampled rows with \code{ell == 0} be removed
#'   internally?
#' @return A data frame containing columns for intrinsic population growth rates
#'   (\code{r}), carrying capacity (\code{k}), log likelihood (\code{ell}), and
#'   biomass (\code{biomass}). Each row contains an iteration for a total length
#'   of \code{reps}.
#' @references
#' Martell, S., & Froese, R. (2013). A simple method for estimating MSY from
#' catch and resilience. Fish and Fisheries, 14(4), 504-514.
#' \url{http://doi.org/10.1111/j.1467-2979.2012.00485.x}
#' @name cmsy
NULL

#' @rdname cmsy
#' @export
#' @examples
#'
#' # An example using cmsy() with a stock from the RAM Legacy database:
#' # The stock is "BGRDRSE" (Blue Grenadier Southeast Australia)
#'
#' d <- structure(list(
#'   yr = 1979:2011,
#'   ct = c(512, 865, 511, 829, 935, 1390, 1260, 2240, 3370, 2770, 3050, 3290,
#'     4540, 3300, 3500, 3190, 2880, 3490, 5670, 6310, 9550, 8700, 9130, 9160,
#'     8490, 6400, 4420, 3680, 3190, 3960, 3290, 4220, 4220)),
#'   .Names = c("yr", "ct"),
#'   class = "data.frame",
#'   row.names = 1:33)
#' x <- cmsy(d$yr, ct = d$ct, prior_log_mean = 0.035, prior_log_sd = 0.68)
#' head(x)
#' par(mfrow = c(1, 3))
#' plot(d, type = "o", xlab = "Year", ylab = "Catch (t)")
#' plot(density(x$biomass))
#' plot(x$r, x$k)
cmsy <- function(
  yr,
  ct,
  prior_log_mean,
  prior_log_sd,
  interyr_index    = 2L,
  interbio         = c(0, 1),
  bio_step         = 0.05,
  resilience       = c(NA, "Very low", "Low", "Medium", "High"),
  start_k          = c(max(ct), 50 * max(ct)),
  startbio         = if (ct[1] / max(ct) < 0.2) c(0.5, 0.9) else c(0.2, 0.6),
  sig_r            = 0.05,
  reps             = 1e4,
  remove_ell0      = TRUE) {

  if (!identical(length(interbio), 2L))
    stop("interbio must be a vector of length 2")
  if (!identical(length(start_k), 2L))
    stop("start_k must be a vector of length 2")
  if (!identical(length(yr), length(ct)))
    stop("yr and ct must be the same length")

  start_r <- switch(resilience[1],
    "Very low"     = c(0.015, 0.1),
    "Low"          = c(0.050, 0.5),
    "Medium"       = c(0.200, 1.0),
    "High"         = c(0.600, 1.5),
    "NA"           = c(0.015, 1.5))

  out <- schaefer_cmsy(
    r_lim          = start_r,
    k_lim          = start_k,
    sig_r          = sig_r,
    startbio       = seq(startbio[1], startbio[2], by = bio_step),
    yr             = yr,
    ct             = ct,
    interyr_index  = interyr_index,
    prior_log_mean = prior_log_mean,
    prior_log_sd   = prior_log_sd,
    interbio       = interbio,
    reps           = reps)

  if (remove_ell0) out <- out[out$ell == 1, ]

  out
}
