#' Summarize B/B_MSY values
#'
#' @param bbmsy A numeric matrix of B/B_MSY values. Columns should represent
#'   samples of B/B_MSY and rows should represent years.
#' @param probs A numeric vector of quantile probabilities.
#' @param log Logical: should the mean and standard deviation be calculated on
#'   the log scale and then exponentiated at the end?
#' @param ... Other parameters to pass to \code{quantile}, \code{sd}, and
#'   \code{mean}. For example, \code{na.rm = TRUE}.
#' @return A data frame: rows are years and columns represent quantile, mean,
#'   and standard deviation values of B/B_MSY.
#' @export
summarize_bbmsy <- function(bbmsy, probs = c(0.025, 0.25, 0.5, 0.75, 0.975),
  log = TRUE, ...) {
  q <- apply(bbmsy, 2, quantile, probs = probs, ...)
  q <- as.data.frame(t(q))
  names(q) <- paste0("bbmsy_q", names(q))
  names(q) <- gsub("%", "", names(q))
  if (log) {
    x_sd <- apply(log(bbmsy), 2, sd, ...)
    x_mean <- apply(log(bbmsy), 2, mean, ...)
  } else {
    x_sd <- apply(bbmsy, 2, sd, ...)
    x_mean <- apply(bbmsy, 2, mean, ...)
  }
  data.frame(q, bbmsy_sd = x_sd, bbmsy_mean = x_mean)
}
