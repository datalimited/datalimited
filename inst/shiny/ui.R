library("shiny")
library("dplyr")

# source("helper.R")

ramts <- readRDS("data/ramts.rds")
ramts <- ramts %>% filter(!is.na(c_touse)) # TODO won't this create gaps in years?
# stocks <- sort(unique(ramts$stocklong))
stocks <- c(
  "Black Grouper Gulf of Mexico",
  "Black oreo West end of Chatham Rise",
  "Dover sole Gulf of Alaska")

shinyUI(fluidPage(
  titlePanel("Data-limited methods"),

  sidebarLayout(
    sidebarPanel(

      selectInput(inputId = "stock",
        label = "Stock",
        choices = stocks,
        selected = NULL),

      sliderInput(
        inputId = "start_r",
        label = "Starting r",
        value = c(0.2, 1.0),
        min = 0.02,
        max = 2.0,
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
      plotOutput("plot_comsir"),
      plotOutput("plot_cmsy")
    )
  )
))
