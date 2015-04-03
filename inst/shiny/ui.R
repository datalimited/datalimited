library("shiny")

# source("helper.R")

shinyUI(fluidPage(
  titlePanel("COM-SIR"),

  sidebarLayout(
    sidebarPanel(

      sliderInput(
        inputId = "start_r",
        label = "Starting r",
        value = c(0.2, 1.0),
        min = 0.05,
        max = 3.0,
        animate = TRUE)

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
      plotOutput("plot")
    )
  )
))
