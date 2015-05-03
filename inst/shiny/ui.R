library("shiny")
library("dplyr")

source("helper.R")

shinyUI(fluidPage(
  titlePanel("Data-limited assessments"),

  sidebarLayout(
    sidebarPanel(

      h2("General options"),

      selectInput(inputId = "stock",
        label = "Stock",
        choices = stocks,
        selected = "Albacore tuna Indian Ocean"),

      selectInput("resilience",
        label = "Resilience",
        choices = c("Unknown", "Very low", "Low", "Medium", "High")),



      h2("COM-SIR options"),

       sliderInput("comsir_cv",
         label = "CV",
         min = 0.05,
         max = 2.0,
         value = c(0.1, 1)),

      sliderInput("comsir_reps",
        label = "Number of initial samples",
        min = 1000,
        max = 500000,
        value = 50000,
        step = 2000),

      sliderInput("comsir_n_posterior",
        label = "Number of posterior samples",
        min = 500,
        max = 20000,
        value = 5000,
        step = 500),

       sliderInput("comsir_a_bounds",
         label = "a range",
         min = 0.001,
         max = 3,
         value = c(0.001, 2)),

       sliderInput("comsir_effort_bounds",
         label = "effort range",
         min = -3,
         max = 3, step = 0.25,
         value = c(-1.0, 1.0)),

      sliderInput("comsir_x_bounds",
        label = "x range",
        min = 0.001,
        max = 3,
        value = c(0.001, 2)),

      sliderInput("comsir_r_bounds",
        label = "r range",
        min = 0.01,
        max = 3,
        value = c(0.2, 1.5)),

      checkboxInput("comsir_obs",
        label = "obs",
        value = FALSE),

      checkboxInput("comsir_dampen",
        label = "dampen",
        value = FALSE),

      h2("Catch-MSY options"),

      sliderInput("sig_r",
        label = "Recruitment variability",
        min = 0,
        max = 0.8,
        value = 0.05),

      sliderInput("cmsy_reps",
        label = "Number of samples",
        min = 500,
        max = 100000,
        value = 2000,
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
        value = TRUE)

#      sliderInput("prior_mean",
#        label = "Mean depletion prior",
#        min = 0.1,
#        max = 4.0,
#        value = 0.8),
#
#      sliderInput("prior_sd",
#        label = "SD depletion prior",
#        min = 0.5,
#        max = 4,
#        value = 2.5)

    ),

    mainPanel(

      tabsetPanel(

        tabPanel("Catch-MSY", plotOutput("plot_cmsy", height = "600px")),
        tabPanel("COMSIR", plotOutput("plot_comsir", height = "1000px")),
        tabPanel("Costello et al. PRM", plotOutput("plot_prm", height = "600px"))

      )

    )
  )
))
