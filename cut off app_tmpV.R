#Shiny app for cut off value assessment tool for clearance half-life

library(shiny)
ui <- fluidPage(
  
  numericInput(inputId = "nn",
               label = "Sample Size:",
               value = 200
               ),
  sliderInput(inputId = "senmu",
              label = "Mean half-life of sensitive distribution",
              value = 3, min = 1, max = 6
  ),
  sliderInput(inputId = "sensd",
              label = "SD of sensitive distribution",
              value = 1.45, min = .1, max = 3
  ),
  sliderInput(inputId = "resmu",
              label = "Mean half-life of resistant distribution",
              value = 6.5, min = 5, max = 10
  ),
  sliderInput(inputId = "ressd",
              label = "SD of resistant distribution",
              value = 1.22, min = .1, max = 3
  ),
  sliderInput(inputId = "cutoff",
              label = "Cut-off half-life value",
              value = 5, min = 0, max = 10
  ),
  sliderInput(inputId = "prop_resist",
              label = "Proportion resistant",
              value = .1, min = 0, max = 1
  ),
  
  plotOutput(outputId = "graph"),
  textOutput(outputId = "cutoffexp")
  
)

server <- function(input, output) {
  senmuR <- reactive({input$senmu})
  sensdR <- reactive({input$sensd})
  resmuR <- reactive({input$resmu})
  ressdR <- reactive({input$ressd})
  
  genData <- reactive({
    senmu <- log(as.numeric(senmuR()))
    sensd <- log(as.numeric(sensdR()))
    resmu <- log(as.numeric(resmuR()))
    ressd <- log(as.numeric(ressdR()))
    
    nn <- input$nn

    prop_resist <- input$prop_resist

    
    sen_pop <- rlnorm(nn*(1-prop_resist),senmu,sensd) #sensitive population
    res_pop <- rlnorm(nn*prop_resist,resmu,ressd) #resistant population
    c(sen_pop,res_pop)
    #total_pop <- c(sen_pop,res_pop)
  })
  
  cutoffR <- reactive({input$cutoff})
  
  output$graph <- renderPlot({
    hist(genData(), freq=FALSE,col="grey",lwd=2,ps=20,breaks=as.numeric(floor(min(genData())):ceiling(max(genData()))))
    lines(density(genData()),lwd=5, col="red")
    abline(v=cutoffR, lwd=3, col="blue")
  })
  

}

shinyApp(ui = ui, server = server)