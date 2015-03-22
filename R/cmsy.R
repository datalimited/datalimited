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
#' @seealso \code{\link{shaefer_cmsy}}
#' @export
#' @examples
#' x <- cmsy(1:30, ct = rlnorm(30), prior_log_mean = 1, prior_log_sd = 0.5)
#' head(x)
cmsy <- function(
  yr,
  ct,
  prior_log_mean,
  prior_log_sd,
  interyr_index = 2,
  interbio = c(0, 1),
  bio_step = 0.05,
  resilience = c(NA, "Very low", "Low", "Medium", "High"),
  start_k = c(max(ct), 50 * max(ct)),
  startbio = if (ct[1] / max(ct) < 0.2) c(0.5, 0.9) else c(0.2, 0.6),
  sig_r = 0.05,
  reps = 1e4) {

  start_r <- switch(resilience[1],
    "Very low" = c(0.015, 0.1),
    "Low"      = c(0.050, 0.5),
    "Medium"   = c(0.200, 1.0),
    "High"     = c(0.600, 1.5),
    "NA"       = c(0.015, 1.5))

  schaefer_cmsy(r_lim = start_r,
    k_lim = start_k,
    sig_r = sig_r,
    startbt = seq(startbio[1], startbio[2], by = bio_step),
    yr = yr,
    ct = ct,
    interyr_index = interyr_index,
    prior_log_mean = prior_log_mean,
    prior_log_sd = prior_log_sd,
    interbio_lim = interbio,
    reps = reps)
}