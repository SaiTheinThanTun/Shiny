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
  
  plotOutput(outputId = "graph")
)

server <- function(input, output) {
  
  
  output$graph <- renderPlot({
    nn <- input$nn
    senmu <- log(input$senmu)
    sensd <- log(input$sensd)
    
    prop_resist <- input$prop_resist
    resmu <- log(input$resmu)
    ressd <- log(input$ressd)
    cutoff <- input$cutoff
    
    
    sen_pop <- rlnorm(nn*(1-prop_resist),senmu,sensd) #sensitive population
    res_pop <- rlnorm(nn*prop_resist,resmu,ressd) #resistant population
    total_pop <- c(sen_pop,res_pop)
    hist(total_pop, freq=FALSE,col="grey",lwd=2,ps=20,breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12))
    lines(density(total_pop),lwd=5, col="red")
    abline(v=cutoff, untf = TRUE, lwd=3, col="blue")
  })
  
}

shinyApp(ui = ui, server = server)