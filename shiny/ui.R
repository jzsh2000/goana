#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(DT)

# Define UI for application that draws a histogram
shinyUI(fluidPage(

  # Application title
  titlePanel("goana"),

  # Sidebar with a slider input for number of bins
  sidebarLayout(
    sidebarPanel(
        radioButtons(
            inputId = 'species',
            label = 'Species',
            choices = c('Homo sapiens' = 'human',
                        'Mus musculus' = 'mouse'),
            selected = "human"
        ),
        selectizeInput(
            inputId = 'p_type',
            label = 'P-values adjustment method',
            choices = c('Raw P-values' = 'none',
                        'Benjamini & Hochberg' = 'fdr',
                        'Bonferroni' = 'bonferroni'),
            selected = 'fdr'
        ),
        numericInput(
            inputId = 'threshold',
            label = 'P-values threshold',
            value = 0.01,
            min = 0,
            max = 0.5
        ),
        fluidRow(
            column(width = 9,
                   textAreaInput(
                       inputId = 'gene',
                       label = "Gene List",
                       height = '200px',
                       placeholder = 'Your awesome gene list'
                   )),
            column(width = 3,
                   actionLink(inputId = 'clear', label = 'clear'),
                   hr(),
                   actionLink(inputId = 'example',
                              label = tags$span('example',
                                                class = 'red-box')))
        )
    ),

    # Show a plot of the generated distribution
    mainPanel(
        dataTableOutput('goana_res')
    )
  )
))
