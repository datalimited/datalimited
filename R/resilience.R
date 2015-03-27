#' Translate resilience categories into population growth rate ranges
#'
#' @param x A character value designating the resilience of the stock. For
#'   convenience, this is case insensitive. Set to any value other than the
#'   named categories to return a default large-ranging value for population
#'   growth rate (\code{c(0.015, 1.5)}).
#' @export
#' @examples
#' resilience("Very low")
#' resilience("medium")
#' resilience(NA)
#' resilience("asdf")
resilience <- function(x = c("very low", "low", "medium", "high")) {
  switch(tolower(as.character(x[1])),
    "very low"     = c(0.015, 0.1),
    "low"          = c(0.050, 0.5),
    "medium"       = c(0.200, 1.0),
    "high"         = c(0.600, 1.5),
                     c(0.015, 1.5))
}
