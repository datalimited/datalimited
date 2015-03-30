#' Catch-MSY method
#'
#' Apply the catch-MSY method from Martell and Froese (2013).
#'
#' @param yr Numeric vector of years
#' @param ct Numeric vector of catch
#' @param sig_r Standard deviation of process noise
#' @param bio_step Step size between lower and upper biomass limits in the iterim year
#' @param start_r A numeric vector of length 2 giving the lower and upper
#'   bounds on the population growth rate parameter. This can either be
#'   specified manually or by translating resiliency categories via the function
#'   \code{\link{resilience}}
#' @param startbio Starting biomass
#' @param start_k Numeric vector of length 2 giving the lower and upper starting
#'   bounds on stock biomass at carrying capacity
#' @param interyr_index Index of the interim year within time series for which
#'   biomass estimate is available
#' @param interbio A numeric vector that gives the lower and upper biomass
#'   limits in the interim year
#' @param prior_log_mean Prior log mean
#' @param prior_log_sd Prior log sd
#' @param reps Number of repititions to run the sampling
#' @param remove_ell0 Should sampled rows with \code{ell == 0} be removed
#'   internally?
#' @return A list containing a matrix \code{biomass} and a data frame
#'   \code{schaefer}. The matrix has rows for each iteration and columns for
#'   years. The data frame contains columns for intrinsic population growth
#'   rates (\code{r}), carrying capacity (\code{k}), log likelihood
#'   (\code{ell}), and biomass (\code{biomass}). Each row contains an iteration
#'   for a total length of \code{reps}.
#' @useDynLib datalimited
#' @importFrom Rcpp sourceCpp
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
#' x <- cmsy(blue_gren$yr, ct = blue_gren$ct, prior_log_mean = 0.035,
#'   prior_log_sd = 0.68)
#' head(x$schaefer)
#' par(mfrow = c(2, 2))
#' plot(blue_gren$yr, blue_gren$ct, type = "o", xlab = "Year", ylab = "Catch (t)")
#' plot(blue_gren$yr,  apply(x$biomass, 2, median)[-1], type = "o",
#'   ylab = "Estimated biomass", xlab = "Year")
#' hist(x$schaefer$bmsy)
#' plot(x$schaefer$r, x$schaefer$k)

cmsy <- function(
  yr,
  ct,
  prior_log_mean,
  prior_log_sd,
  interyr_index    = 2L,
  interbio         = c(0, 1),
  bio_step         = 0.05,
  start_r          = resilience(NA),
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

  schaefer_out <- schaefer_cmsy(
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

  schaefer_out <- schaefer_out[schaefer_out$ell == 1, ]
  schaefer_out$bmsy <- schaefer_out$k * 0.5
  schaefer_out$msy  <- schaefer_out$r * schaefer_out$k / 4
  schaefer_out$mean_ln_msy <- mean(log(schaefer_out$msy))

  biomass_out <- get_cmsy_biomass(
    r = schaefer_out$r,
    k = schaefer_out$k,
    j = schaefer_out$J,
    sigR = sig_r,
    nyr = length(yr),
    ct = ct)
  biomass_out <- t(biomass_out)

  list(biomass = biomass_out, schaefer = schaefer_out)
}
