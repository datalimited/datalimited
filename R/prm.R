#' Modified panel regression model
#'
#' Format for, fit, and predict from a panel regression model. The
#' \pkg{datalimited} package includes cached data that can be used with
#' \code{fit_prm} and a cached model that can be used with \code{predict_prm}.
#' These functions can also be used with new data.
#'
#' @name prm
#' @references Costello, C., D. Ovando, R. Hilborn, S. D. Gaines, O. Deschenes,
#' and S. E. Lester. 2012. Status and Solutions for the World's Unassessed
#' Fisheries. Science 338:517-520.
#' @examples
#' # combine two built in datasets:
#' # d <- dplyr::inner_join(ram_ts, spp_categories, by = "scientificname")
#'
#' # ram_prm_dat is built into the package; it can be created with:
#' # ram_prm_dat <- plyr::ddply(d, "stockid", function(x) {
#' #   format_prm(year = x$year, catch = x$catch, bbmsy = x$bbmsy_ram,
#' #     species_cat = x$spp_category[1L])
#' # })
#'
#' # ram_prm_model is built into the package, it is created in the same manner
#' # as this:
#' # m <- fit_prm(ram_prm_dat)
#'
#' # now predict B/Bmsy:
#' d <- subset(ram_prm_dat, stockid == "BGRDRSE")
#' x <- predict_prm(d) # use the built in model
#' plot(x)
#'
#' # with confidence intervals:
#' library("ggplot2")
#' x <- predict_prm(d, ci = TRUE)
#' ggplot(x$bbmsy, aes(year, bbmsy_q50)) +
#'   geom_ribbon(aes(ymin = bbmsy_q2.5, ymax = bbmsy_q97.5), fill = "#00000040") +
#'   geom_line() +
#'   ylab(expression(B/B[MSY]))
NULL

#' @description \code{format_prm}: Format a series of predictors for a panel
#'   regression
#'
#' @param year A numeric vector of years
#' @param catch A numeric vector of catches
#' @param bbmsy A numeric vector of B/B_MSY
#' @param species_cat A single character value of a
#'   species category. Can be from any set of species categories as long as the
#'   same set are used in model fitting and extrapolation.
#' @export
#' @return \code{format_prm}: A data frame formatted for use with \code{fit_prm} or
#'   \code{predict_prm}
#' @rdname prm

format_prm <- function(year, catch, bbmsy, species_cat) {

  stopifnot(identical(length(catch), length(bbmsy)))
  stopifnot(identical(length(species_cat), 1L))

  .n <- length(catch)
  max_catch <- max(catch)
  scaled_catch <- catch / max_catch
  scaled_catch1 <- dplyr::lag(scaled_catch, n = 1L)
  scaled_catch2 <- dplyr::lag(scaled_catch, n = 2L)
  scaled_catch3 <- dplyr::lag(scaled_catch, n = 3L)
  scaled_catch4 <- dplyr::lag(scaled_catch, n = 4L)
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

#' @description \code{fit_prm}: Fit a panel regression
#'
#' @param dat A data frame created by \code{\link{format_prm}}
#' @param eqn A formula describing the regression
#' @param type Either linear model (the original method) or a GBM model
#' @param ... Anything extra to pass to \code{\link[gbm]{gbm}}
#' @export
#' @return \code{fit_prm}: A linear model for use with \code{predict_prm}
#' @rdname prm

fit_prm <- function(dat,
  eqn = log(bbmsy) ~
    max_catch +
    mean_scaled_catch +
    scaled_catch +
    scaled_catch1 +
    scaled_catch2 +
    scaled_catch3 +
    scaled_catch4 +
    species_cat +
    catch_to_rolling_max +
    time_to_max +
    years_back +
    initial_slope - 1,
    type = c("lm", "gbm"), ...) {

  switch(type[1],
    lm = lm(eqn, data = dat),
    gbm = gbm::gbm(eqn, data = dat, distribution = "gaussian", ...))
}

#' \code{predict_prm}: Predict from a panel regression model
#'
#' @param newdata A data frame to predict on that has been formatted with
#'   \code{\link{format_prm}}.
#' @param model A linear model to predict from built from \code{\link{fit_prm}}.
#'   Defaults to a cached version of a model fit to the RAM database,
#'   \code{\link{ram_prm_model}}.
#' @param ci Should confidence intervals on B/B_MSY be returned?
#' @param level Confidence interval level
#' @return \code{predict_prm}: If \code{ci = FALSE}, a vector of predictions of
#'   B/B_MSY. If \code{ci = TRUE}, a data frame with columns for B/B_MSY
#'   \code{fit} and lower and upper confidence intervals \code{lower},
#'   \code{upper}. Note that the model is fitted to log(B/B_MSY) and the output
#'   from \code{predict_prm} is exponentiated so the prediction represents an
#'   estimate of median B/B_MSY.
#' @export
#' @rdname prm


predict_prm <- function(newdata, model = datalimited::ram_prm_model,
  ci = FALSE, level = 0.95) {

  if (class(model) == "lm") {
    p <- predict(model, newdata = newdata, se = ci)
  }
  if (class(model) == "gbm") {
    p <- exp(gbm::predict.gbm(model, newdata = newdata, n.trees = model$n.trees))
  }

  if (!ci) {
    out <- exp(p)
  } else {
    out_main <- prm_ci(p, level = level)

    out1 <- prm_ci(p, level = 0.95)
    out2 <- prm_ci(p, level = 0.50)
    out1 <- dplyr::rename_(out1, bbmsy_q2.5 = "lower", bbmsy_q50 = "fit",
      bbmsy_q97.5 = "upper")
    out2 <- dplyr::rename_(out2, bbmsy_q25 = "lower", bbmsy_q75 = "upper")
    out2$fit <- NULL
    out1$median <- out1$fit
    bbmsy <- data.frame(year = newdata$year, catch = newdata$catch, out1, out2)
    bbmsy <- bbmsy[,c("year", "catch", "bbmsy_q2.5", "bbmsy_q25", "bbmsy_q50",
      "bbmsy_q75", "bbmsy_q97.5")]
    out <- list()
    out$prm <- out_main
    out$bbmsy <- bbmsy
  }
  out
}

prm_ci <- function(p, level = 0.95) {
  q <- -qnorm((1 - level)/2)
  data.frame(fit = exp(p$fit),
    lower = exp(p$fit - q * p$se.fit),
    upper = exp(p$fit + q * p$se.fit))
}
