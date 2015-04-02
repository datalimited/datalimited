library("shiny")

# source("helper.R")

shinyUI(fluidPage(
  titlePanel("COM-SIR"),

  sidebarLayout(
    sidebarPanel(

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
      dataTableOutput("table")
    )
  )
))
