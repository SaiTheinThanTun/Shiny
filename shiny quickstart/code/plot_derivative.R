library(shiny)
ui <- fluidPage(
  textInput(inputId = "f",
              label = "Input formula in terms of x (eg. x^2)",
              value = "sin(x)-x/25"
            ),
  sliderInput(inputId = "xRange",
              label = "Range of x",
              value = c(-50,50), min = -100, max = 100
              ),
  plotOutput(outputId = "graph"),
  plotOutput(outputId = "differential")
  )

server <- function(input, output) {
  output$graph <- renderPlot({
    f <- function(x, y=parse(text=input$f)){eval(y)}
    range1 <- input$xRange[1]
    range2 <- input$xRange[2]
    x <- seq(range1,range2,by=.001)  
    plot(x,f(x), type="l", main = "Graph of the function")
  })
  output$differential <- renderPlot({
    g <- function(x, y=parse(text=input$f)) {eval(D(y,"x"))}
    range1 <- input$xRange[1]
    range2 <- input$xRange[2]
    x <- seq(range1,range2,by=.001)
    plot(x,g(x), type="l", main = "Graph of the derivatives of the function")
  })
}

shinyApp(ui = ui, server = server)