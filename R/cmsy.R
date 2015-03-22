#' Catch-msy method
#'
#' @param yr
#' @param ct
#' @param prior_log_mean
#' @param prior_log_sd
#' @param interyr_index
#' @param interbio
#' @param finalbio
#' @param bio_step
#' @param resilience
#' @param start_k
#' @param startbio
#' @param sig_r
#' @param reps
#'
#' @return A data frame containing columns for intrinsic population growth rates
#'   (\code{r}), carrying capacity (\code{k}), log likelihood (\code{ell}), and
#'   biomass (\code{biomass}). Each row contains an iteration for a total length
#'   of \code{reps}.
#'
#' @seealso \code{\link{shaefer_cmsy}}
#' @export
#' @examples
#'
#' # "BGRDRSE" (Blue Grenadier Southeast Australia) from the RAM Legacy database:
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
#' x1 <- x[x$ell == 1, ]
#' par(mfrow = c(1, 3))
#' plot(d, type = "o")
#' plot(density(x1$biomass))
#' plot(x1$r, x1$k)
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
  reps             = 1e4) {

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

  schaefer_cmsy(
    r_lim          = start_r,
    k_lim          = start_k,
    sig_r          = sig_r,
    startbt        = seq(startbio[1], startbio[2], by = bio_step),
    yr             = yr,
    ct             = ct,
    interyr_index  = interyr_index,
    prior_log_mean = prior_log_mean,
    prior_log_sd   = prior_log_sd,
    interbio_lim   = interbio,
    reps           = reps)
}
