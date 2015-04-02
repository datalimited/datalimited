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
#' @param revise_bounds Should the bounds on r and k be revised after fitting
#'   the algorithm once? The algorithm will then fit a second time with the
#'   revised bounds.
#' @return A list containing a matrix \code{biomass} and a data frame
#'   \code{quantities}. The matrix has rows for each iteration and columns for
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
#' head(x$quantities)
#' par(mfrow = c(2, 2))
#' plot(blue_gren$yr, blue_gren$ct, type = "o", xlab = "Year", ylab = "Catch (t)")
#' plot(blue_gren$yr,  apply(x$biomass, 2, median)[-1], type = "o",
#'   ylab = "Estimated biomass", xlab = "Year")
#' hist(x$quantities$bmsy)
#' plot(x$quantities$r, x$quantities$k)

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
  revise_bounds     = TRUE) {

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

  if (revise_bounds) {
    # Repeat with revised bounds
    # use for batch processing   ****THIS IS AN UPDATE FROM RAINER--Oct 12, 2012:
    # finalbio <- if(ct[length(yr)]/max(ct) > 0.5) {c(0.3,0.7)} else {c(0.01,0.4)}
    # parbound <- list(r = start_r, k = start_k, lambda = finalbio, sig_r=sig_r)
    parbound <- list(r = start_r, k = start_k)
    r1 	<- schaefer_out$r[schaefer_out$ell==1]
    k1 	<- schaefer_out$k[schaefer_out$ell==1]
    j1  <- schaefer_out$J[schaefer_out$ell==1]
    msy1 <- r1*k1/4
    mean_msy1 <- exp(mean(log(msy1)))
    max_k1a <- min(k1[r1<1.1*parbound$r[1]]) ## smallest k1 near initial lower bound of r
    max_k1b <- max(k1[r1*k1/4<mean_msy1]) ## largest k1 that gives mean MSY
    max_k1 <- if(max_k1a < max_k1b) {max_k1a} else {max_k1b}
    if (length(r1)<10) {
      warning(paste0("Too few (", length(r1), ") possible r-k combinations, ",
        "check input parameters"))
      return(NULL)
    } else {
      # set new upper bound of r to 1.2 max r1
      parbound$r[2] <- 1.2 * max(r1)
      # set new lower bound for k to 0.9 min k1 and upper bound to max_k1
      parbound$k <- c(0.9 * min(k1), max_k1)
      # Repeat analysis with new r-k bounds
      schaefer_out <- schaefer_cmsy(
        r_lim          = parbound$r,
        k_lim          = parbound$k,
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
    }
  }

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

  list(biomass = biomass_out, quantities = schaefer_out)
}
