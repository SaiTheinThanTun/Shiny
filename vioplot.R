library(shiny)
library(boot)
library(MASS)
library(segmented)
library(sm)
library(mixtools)
library(vioplot)

ui <- fluidPage(
  fileInput(inputId = "file", label = "Select your input file:"),
  plotOutput(outputId = "violin")
  )

server <- function(input, output) {
  output$violin <- renderPlot({
    mixdat <- input$file
    vioplot(na.omit(mixdat[,1]),na.omit(mixdat[,2]), na.omit(mixdat[,3]), na.omit(mixdat[,4]), na.omit(mixdat[,5]), na.omit(mixdat[,6]), na.omit(mixdat[,7]),na.omit(mixdat[,8]), na.omit(mixdat[,9]), na.omit(mixdat[,10]), na.omit(mixdat[,11]), na.omit(mixdat[,12]),col="grey",pchMed=30)
  })
}

shinyApp(ui = ui, server = server)
