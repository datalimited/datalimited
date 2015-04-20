#' Format a series of predictors for a panel regression
#'
#' @param year A numeric vector of years
#' @param catch A numeric vector of catches
#' @param bbmsy A numeric vector of B/B_MSY
# @param von_bert_k A single value (or a numeric vector) of Von Bertalanffy
#   growth rate parameter
# @param temperature A single value (or a numeric vector) of the mean preferred
#   temperature of the species
# @param max_length A single value (or a numeric vector) of maximum recorded
#   length for the species
# @param age_maturity A single value (or a numeric vector) of age at which 50%
#   of the individuals are sexually mature for the species
#' @param species_cat A single character value of a
#'   species category. Can be from any set of species categories as long as the
#'   same set are used in model fitting and extrapolation.
#' @export
#' @references Costello, C., D. Ovando, R. Hilborn, S. D. Gaines, O. Deschenes,
#' and S. E. Lester. 2012. Status and Solutions for the World's Unassessed
#' Fisheries. Science 338:517-520.
#' @return A data frame formatted for use with \code{\link{fit_prm}} or
#'   \code{\link{predict_prm}}
#' @family prm
#' @examples
#' d <- dplyr::inner_join(ramts, spp_categories, by = "scientificname")
#' ram_prm_dat <- plyr::ddply(d, "stockid", function(x) {
#'   format_prm(year = x$year, catch = x$catch, bbmsy = x$bbmsy_ram,
#'     species_cat = x$spp_category[1L])
#' })

format_prm <- function(year, catch, bbmsy, species_cat) {

  stopifnot(identical(length(catch), length(bbmsy)))
  stopifnot(identical(length(species_cat), 1L))

  .n <- length(catch)
  max_catch <- max(catch)
  scaled_catch <- catch / max_catch
  scaled_catch1 <- lag(scaled_catch, n = 1L)
  scaled_catch2 <- lag(scaled_catch, n = 2L)
  scaled_catch3 <- lag(scaled_catch, n = 3L)
  scaled_catch4 <- lag(scaled_catch, n = 4L)
  catch_to_rolling_max <- scaled_catch / cummax(scaled_catch)
  mean_scaled_catch <- mean(scaled_catch)
  time_to_max <- min(seq_along(catch)[catch == max(catch)])
  years_back <- seq(length(catch), 1L)
  initial_slope <- coef(lm(scaled_catch[seq_len(6L)] ~ seq_len(6L)))[[2L]]

  data.frame(year, bbmsy, years_back, catch,
    max_catch = rep(max_catch, .n),
    scaled_catch,
    mean_scaled_catch = rep(mean_scaled_catch, .n),
    scaled_catch1, scaled_catch2, scaled_catch3, scaled_catch4,
    catch_to_rolling_max,
    time_to_max = rep(time_to_max, .n),
    initial_slope = rep(initial_slope, .n),
    species_cat = rep(species_cat, .n))
}

#' Fit panel regression
#'
#' @param dat A data frame created by \code{\link{format_prm}}
#' @export
#' @return A linear model
#' @seealso \code{\link{format_prm}}, \code{\link{predict_prm}}
#' @references Costello, C., D. Ovando, R. Hilborn, S. D. Gaines, O. Deschenes,
#' and S. E. Lester. 2012. Status and Solutions for the World's Unassessed
#' Fisheries. Science 338:517-520.
#' @family prm
#' @examples
#' m <- fit_prm(ram_prm_dat)

fit_prm <- function(dat) {
  lm(log(bbmsy) ~
      species_cat +
      max_catch +
      mean_scaled_catch +
      scaled_catch +
      scaled_catch1 +
      scaled_catch2 +
      scaled_catch3 +
      scaled_catch4 +
      catch_to_rolling_max +
      time_to_max +
      years_back +
      initial_slope - 1,
    data = dat)
}

#' Predict from a panel regression model
#'
#' @param newdata A data frame to predict on that has been formatted with
#'   \code{\link{format_prm}}.
#' @param model A linear model to predict from built from \code{\link{fit_prm}}.
#'   Defaults to a cached version of a model fit to the RAM database,
#'   \code{\link{ram_prm_model}}.
#' @param ci Should confidence intervals on B/B_MSY be returned?
#' @param level Confidence interval level
#' @return If \code{ci = FALSE}, a vector of predictions of B/B_MSY. If \code{ci
#'   = TRUE}, a data frame with columns for B/B_MSY \code{fit} and lower and
#'   upper confidence intervals \code{lower}, \code{upper}. Note that the model
#'   is fitted to log(B/B_MSY) and the output from \code{predict_prm} is
#'   exponentiated so the prediction represents an estimate of median B/B_MSY.
#' @export
#' @references Costello, C., D. Ovando, R. Hilborn, S. D. Gaines, O. Deschenes,
#' and S. E. Lester. 2012. Status and Solutions for the World's Unassessed
#' Fisheries. Science 338:517-520.
#' @family prm
#' @examples
#' d <- subset(ram_prm_dat, stockid == "BGRDRSE")
#' x <- predict_prm(d)
#' plot(x)
#'
# with confidence intervals:
#' library("ggplot2")
#' x <- predict_prm(d, ci = TRUE)
#' x$year <- rev(d$years_back)
#' ggplot(x, aes(year, fit)) +
#'   geom_ribbon(aes(ymin = lower, ymax = upper), fill = "#00000040") +
#'   geom_line() +
#'   ylab(expression(B/B[MSY]))

predict_prm <- function(newdata, model = datalimited::ram_prm_model,
  ci = FALSE, level = 0.95) {

  p <- predict(model, newdata = newdata, se = ci)

  if (!ci) {
    return(exp(p))
  } else {
    q <- -qnorm((1 - level)/2)
    data.frame(fit = exp(p$fit),
      lower = exp(p$fit - q * p$se.fit),
      upper = exp(p$fit + q * p$se.fit))
  }
}
