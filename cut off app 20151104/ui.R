library(shiny)
fluidPage(
  fluidRow(
    column(5,
           plotOutput(outputId = "graph")
    ),
    column(4,
           plotOutput(outputId = "ROC")
    )
  ),
  fluidRow(
    column(3,
           h4("Sensitive Distribution"),
           sliderInput(inputId = "senmu",
                       label = "Mean half-life",
                       value = 3, min = 1, max = 6
           ),
           sliderInput(inputId = "sensd",
                       label = "SD",
                       value = 1.45, min = 1, max = 2.1
           )
    ),
    column(3,
           h4("Resistant Distribution"),
           sliderInput(inputId = "resmu",
                       label = "Mean half-life",
                       value = 6.5, min = 5, max = 10
           ),
           sliderInput(inputId = "ressd",
                       label = "SD",
                       value = 1.22, min = 1, max = 2.1
           ),
           sliderInput(inputId = "prop_resist",
                       label = "Proportion resistant",
                       value = .1, min = 0, max = 1
           )
    ),
    column(3,
           numericInput(inputId = "nn",
                        label = "Sample Size:",
                        value = 200
           ),
           
           
           sliderInput(inputId = "cutoff",
                       label = "Cut-off half-life value",
                       value = 5, min = 0, max = 10
           )
    )
    
  ),
  
  
  
  
  verbatimTextOutput(outputId = "check")
  
)