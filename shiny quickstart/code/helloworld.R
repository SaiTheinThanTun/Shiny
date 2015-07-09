library(shiny)
ui <- fluidPage(
  textInput(inputId = "f",
              label = "Input formula in terms of x (eg. x^2)",
              value = "sin(x)-x/2"
            ),
  sliderInput(inputId = "xRange",
              label = "Range of x",
              value = c(0,10), min = -50, max = 50
              ),
  plotOutput(outputId = "graph")
  )

server <- function(input, output) {
  output$graph <- renderPlot({
    #f <- function(x) {as.function(input$f)}
    #f <- function(x) {x^2}
    f <- function(x, y=parse(text=input$f)){eval(y)}
    range1 <- input$xRange[1]
    range2 <- input$xRange[2]
    x <- seq(range1,range2,by=.001)  
    plot(x,f(x))
  })
    
}

shinyApp(ui = ui, server = server)