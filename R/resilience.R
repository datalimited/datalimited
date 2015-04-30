#' Translate resilience categories into population growth rate ranges
#'
#' @param x A character value designating the resilience of the stock. For
#'   convenience, this is case insensitive. Set to any value other than the
#'   named categories to return a default large-ranging value for population
#'   growth rate (\code{c(0.015, 1.5)}).
#' @param unknown_equals_medium Logical: should an unknown designation be
#'   treated the same as the medium category or encompass the full range of
#'   popopulation growth rate (0.015, 1.5)
#' @export
#' @examples
#' resilience("Very low")
#' resilience("medium")
#' resilience(NA)
#' resilience("asdf")
resilience <- function(x = c("very low", "low", "medium", "high"),
  unknown_equals_medium = TRUE) {

  if(unknown_equals_medium)
    unknown_range <- c(0.2, 1)
  else
    unknown_range <- c(0.015, 1.5)

  switch(tolower(as.character(x[1])),
    "very low"     = c(0.015, 0.1),
    "low"          = c(0.050, 0.5),
    "medium"       = c(0.200, 1.0),
    "high"         = c(0.600, 1.5),
                     unknown_range)
}
