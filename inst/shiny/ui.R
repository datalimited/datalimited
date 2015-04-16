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
        selected = NULL),

      selectInput("resilience",
        label = "Resilience",
        choices = c("Unknown", "Very low", "Low", "Medium", "High")),

      h2("Catch-MSY options"),

      sliderInput("sig_r",
        label = "Recruitment variability",
        min = 0,
        max = 0.8,
        value = 0.05),

      sliderInput("cmsy_reps",
        label = "Number of samples",
        min = 500,
        max = 20000,
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
         value = 2.5),

      h2("COM-SIR options"),

       sliderInput("comsir_a",
         label = "Fraction of K at bioeconomic equilibrium",
         min = 0.05,
         max = 1.0,
         value = 0.8),

       sliderInput("comsir_k",
         label = "Stock size at carrying capacity",
         min = 10,
         max = 3000,
         value = 800),

       sliderInput("comsir_x",
         label = "Intrinsic rate of increase in effort",
         min = 0.05,
         max = 1.0,
         value = 0.5),

       sliderInput("comsir_cv",
         label = "CV",
         min = 0.05,
         max = 2.0,
         value = 0.4),

      sliderInput("comsir_reps",
        label = "Number of initial samples",
        min = 1000,
        max = 100000,
        value = 50000,
        step = 2000),

      sliderInput("comsir_n_posterior",
        label = "Number of posterior samples",
        min = 500,
        max = 20000,
        value = 5000,
        step = 500),

       checkboxInput("comsir_logk",
         label = "Log K",
         value = TRUE),

       checkboxInput("comsir_norm_r",
         label = "Normal r",
         value = FALSE),

       checkboxInput("comsir_norm_k",
         label = "Normal k",
         value = FALSE),

       checkboxInput("comsir_norm_x",
         label = "Normal x",
         value = FALSE),

       checkboxInput("comsir_norm_a",
         label = "Normal a",
         value = FALSE),

       checkboxInput("comsir_normal_like",
         label = "Normal likelihood",
         value = FALSE),

       checkboxInput("comsir_logistic",
         label = "Logistic",
         value = TRUE),

       checkboxInput("comsir_obs",
         label = "Obs??",
         value = FALSE)


    ),

    mainPanel(



      tabsetPanel(
        tabPanel("Catch-MSY", plotOutput("plot_cmsy")),
        tabPanel("COMSIR", plotOutput("plot_comsir"))
      )

    )
  )
))
