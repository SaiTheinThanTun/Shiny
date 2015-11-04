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
              value = 1.45, min = 1, max = 2.1
  ),
  sliderInput(inputId = "resmu",
              label = "Mean half-life of resistant distribution",
              value = 6.5, min = 5, max = 10
  ),
  sliderInput(inputId = "ressd",
              label = "SD of resistant distribution",
              value = 1.22, min = 1, max = 2.1
  ),
  sliderInput(inputId = "cutoff",
              label = "Cut-off half-life value",
              value = 5, min = 0, max = 10
  ),
  sliderInput(inputId = "prop_resist",
              label = "Proportion resistant",
              value = .1, min = 0, max = 1
  ),
  
  plotOutput(outputId = "ROC"),
  plotOutput(outputId = "graph"),
  verbatimTextOutput(outputId = "check")
  
)

server <- function(input, output) {
  senmuR <- reactive({log(input$senmu)})
  sensdR <- reactive({log(input$sensd)})
  resmuR <- reactive({log(input$resmu)})
  ressdR <- reactive({log(input$ressd)})
  
  sen_popR <- reactive({rlnorm(input$nn*(1-input$prop_resist),senmuR(),sensdR())})
  res_popR <- reactive({rlnorm(input$nn*input$prop_resist,resmuR(),ressdR())})
  
  
  genData <- reactive({
    #nn <- input$nn
    senmu <- senmuR()
    sensd <- sensdR()
    
    #prop_resist <- input$prop_resist
    resmu <- resmuR()
    ressd <- ressdR()
    
    sen_pop <- sen_popR() #sensitive population
    res_pop <- res_popR() #resistant population
    c(sen_pop,res_pop)
    #total_pop <- c(sen_pop,res_pop)
  })
  output$graph <- renderPlot({
    
    hist(genData(), freq=FALSE,col="grey",lwd=2,ps=20,breaks=as.numeric(floor(min(genData())):ceiling(max(genData()))))
    lines(density(genData()),lwd=5, col="red")
    abline(v=input$cutoff, lwd=3, col="blue")
    
  })
  output$check <- renderPrint({
    true_res <- plnorm(input$cutoff,resmuR(),ressdR(), lower.tail=FALSE) #% identified as resistant from truly resistant pop
    fal_res <- plnorm(input$cutoff,senmuR(),sensdR(), lower.tail=FALSE) #% identified as resistant from sensitive pop
    
    fal_sen <- plnorm(input$cutoff, resmuR(),ressdR()) #% wrongly identified as sensitive from truly resistant pop
    true_sen <- plnorm(input$cutoff, senmuR(), sensdR()) #% identified as sensitive from truely sensitive pop
    cat(round(100*true_res,2),"% true resistant, ", round(100*fal_res,2),"% false resistant \n",round(100*fal_sen,2),"% false sensitive, ",round(100*true_sen,2),"% true sensitive")
  })
  
  output$ROC <- renderPlot({
    popDF <- cbind(genData(),c(rep(0,length(sen_popR())),rep(1,length(res_popR()))))
    
    TPR <- sum(res_popR()>=input$cutoff)/length(res_popR())
    FPR <- sum(sen_popR()>=input$cutoff)/length(sen_popR())
    
    roc(popDF[,2], popDF[,1],  partial.auc.correct=TRUE, partial.auc.focus="sens",ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE, show.thres=TRUE)
    points((1-FPR),TPR, col="red", pch=19)
  })
}

shinyApp(ui = ui, server = server)