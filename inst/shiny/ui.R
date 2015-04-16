library("shiny")
library("dplyr")

source("helper.R")

shinyUI(fluidPage(
  titlePanel("Catch-MSY"),

  sidebarLayout(
    sidebarPanel(

      selectInput(inputId = "stock",
        label = "Stock",
        choices = stocks,
        selected = NULL),

      selectInput("resilience",
        label = "Resilience",
        choices = c("Unknown", "Very low", "Low", "Medium", "High")),

      sliderInput("sig_r",
        label = "Recruitment variability",
        min = 0,
        max = 0.8,
        value = 0.05),

      sliderInput("reps",
        label = "Number of samples",
        min = 500,
        max = 20000,
        value = 3000,
        step = 500),

      sliderInput("interyr_index",
        label = "Biomass depletion reference year",
        min = 1L,
        max = 10L,
        value = 2L,
        step = 1L),

      sliderInput("interbio",
        label = "Limits on biomass depletion in ref. year",
        min = 0,
        max = 1,
        value = c(0, 1)),

       checkboxInput("revise_bounds",
         label = "Revise r and K bounds",
         value = TRUE),

       sliderInput("prior_mean",
         label = "Mean depletion prior",
         min = 0.1,
         max = 4.0,
         value = 0.8),

       sliderInput("prior_sd",
         label = "SD depletion prior",
         min = 0.5,
         max = 4,
         value = 2.5)

#       checkboxGroupInput("status",
#         label = "Choose IUCN statuses",
#         choices = statuses,
#         selected = statuses),
#
#       checkboxGroupInput("habitat",
#         label = "Choose habitats",
#         choices = habitats,
#         selected = habitats)

    ),

    mainPanel(
      # plotOutput("plot_comsir"),
      plotOutput("plot_cmsy")
    )
  )
))
