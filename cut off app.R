#Shiny app for cut off value assessment tool for clearance half-life

library(shiny)
library(pROC)
library(ggplot2)

ui <- fluidPage(
  fluidRow(
    column(5,
           plotOutput(outputId = "graph")
          ),
    column(5,
           plotOutput(outputId = "ROC")
           ),
    column(5,
           plotOutput(outputId = "stacked")
           ),
    column(5,
           plotOutput(outputId = "densityplot")
           )
  ),
  fluidRow(
    column(3,
      h4("Sensitive Distribution"),
      sliderInput(inputId = "senmu",
                  label = "Mean half-life",
                  value = 3, min = 1, max = 6.5, step = .5
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
                  value = 6.5, min = 5, max = 10, step = .5
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
  
  genData.DF <- reactive({
    cbind(genData(),c(rep(0,length(sen_popR())),rep(1,length(res_popR()))))
  })
  
  output$graph <- renderPlot({
    
    hist(genData(), freq=FALSE,col="grey",lwd=2,ps=20,breaks=as.numeric(floor(min(genData())):ceiling(max(genData()))), main="Histogram of Simulated Half-Lives", xlab="Half-life (hours)")
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
    popDF <- genData.DF()
    
    TPR <- sum(res_popR()>=input$cutoff)/length(res_popR())
    FPR <- sum(sen_popR()>=input$cutoff)/length(sen_popR())
    
    roc(popDF[,2], popDF[,1],  partial.auc.correct=TRUE, partial.auc.focus="sens",ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE, show.thres=TRUE, main="Receiver Operating Characteristic (ROC) Curve")
    points((1-FPR),TPR, col="red", pch=19)
  })
  
  output$stacked <- renderPlot({
    popDF2 <- genData.DF() 
    
    popDF2[popDF2[,2]==0,2] <- "Sensitive"
    popDF2[popDF2[,2]==1,2] <- "Resistant"
    
    popDF2 <- as.data.frame(popDF2)
    popDF2[,1] <- as.numeric(as.character(popDF2[,1]))
    names(popDF2) <- c("Half-life (hours)","Sensitivity")
    
    #ggplot(popDF2, aes(x=`Half-life (hours)` , y= ..density.. , fill=Sensitivity)) +
    #ggplot(popDF2, aes(x=`Half-life (hours)`, fill=Sensitivity)) +
    ggplot(popDF2, aes(x=`Half-life (hours)`, fill=Sensitivity)) +
      geom_histogram(position="stack")  +
      scale_x_continuous(breaks=as.numeric(floor(min(genData())):ceiling(max(genData()))))
    
  })
  output$densityplot <- renderPlot({
    popDF2 <- genData.DF() 
    
    #mxmdl <- normalmixEM(popDF2)
    #plot(mxmdl, which=2)
    popDF2[popDF2[,2]==0,2] <- "Sensitive"
    popDF2[popDF2[,2]==1,2] <- "Resistant"
    
    popDF2 <- as.data.frame(popDF2)
    popDF2[,1] <- as.numeric(as.character(popDF2[,1]))
    names(popDF2) <- c("Half-life (hours)","Sensitivity")
    
    ggplot(popDF2, aes(x=`Half-life (hours)`, fill=Sensitivity, colour= Sensitivity)) +
      geom_histogram(aes(y=(..count..)/sum(..count..)), alpha=.6, position="stack",breaks=as.numeric(floor(min(genData())):ceiling(max(genData())))) +
      scale_x_continuous(breaks=as.numeric(floor(min(genData())):ceiling(max(genData()))))
  })
}

shinyApp(ui = ui, server = server)