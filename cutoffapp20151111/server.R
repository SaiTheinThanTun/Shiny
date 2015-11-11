library(shiny)
library(pROC)
library(ggplot2)

function(input, output) {
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
  output$ROC <- renderPlot({
    popDF <- genData.DF()
    
    TPR <- sum(res_popR()>=input$cutoff)/length(res_popR())
    FPR <- sum(sen_popR()>=input$cutoff)/length(sen_popR())
    
    true_res <- sum(res_popR()>=input$cutoff)
    fal_res <- sum(sen_popR()>=input$cutoff)
    fal_sen <- sum(res_popR()<input$cutoff)
    true_sen <- sum(sen_popR()<input$cutoff)
    overlay <- paste(true_res," truly resistant, ", fal_res, " falsely resistant \n",fal_sen," falsely sensitive, ",true_sen," truly sensitive")
    
    
    roc(popDF[,2], popDF[,1],  partial.auc.correct=TRUE, partial.auc.focus="sens",ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE, plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE, show.thres=TRUE, main="Receiver Operating Characteristic (ROC) Curve")
    points((1-FPR),TPR, col="red", pch=19)
    text(.5,.5,overlay, col="red")
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
    
    if(input$bDensity){
      if(input$bStacked){
        ggplot(popDF2, aes(x=`Half-life (hours)`, fill=Sensitivity, colour= Sensitivity)) + theme_bw() +
          geom_histogram(aes(y=(..count..)/sum(..count..)), alpha=.8, position="stack",breaks=as.numeric(floor(min(genData())):ceiling(max(genData())))) +
          geom_vline(xintercept= input$cutoff, colour="red", size=1) + ylab("Density") + ggtitle("Stacked Histogram of Simulated Half-Lives")+
          theme(plot.title= element_text(face="bold")) +
          scale_x_continuous(breaks=as.numeric(floor(min(genData())):ceiling(max(genData()))))
      }
      
      else {
        ggplot(popDF2, aes(x=`Half-life (hours)`)) + theme_bw() +
          geom_histogram(aes(y=(..count..)/sum(..count..), fill=Sensitivity, colour= Sensitivity), alpha=.4, position="identity",breaks=as.numeric(floor(min(genData())):ceiling(max(genData())))) +
          geom_vline(xintercept= input$cutoff, colour="red", size=1) + ylab("Density") + ggtitle("Histogram of Simulated Half-Lives")+
          theme(plot.title= element_text(face="bold")) +
          scale_x_continuous(breaks=as.numeric(floor(min(genData())):ceiling(max(genData()))))
      }
    }
    
    else{
      if(input$bStacked){
        ggplot(popDF2, aes(x=`Half-life (hours)`, fill=Sensitivity, colour= Sensitivity)) + theme_bw() +
          geom_histogram(alpha=.8, position="stack",breaks=as.numeric(floor(min(genData())):ceiling(max(genData())))) +
          geom_vline(xintercept= input$cutoff, colour="red", size=1) + ggtitle("Stacked Histogram of Simulated Half-Lives")+
          theme(plot.title= element_text(face="bold")) +
          scale_x_continuous(breaks=as.numeric(floor(min(genData())):ceiling(max(genData()))))
      }
      
      else {
        ggplot(popDF2, aes(x=`Half-life (hours)`, fill=Sensitivity, colour= Sensitivity)) + theme_bw() +
          geom_histogram( alpha=.4, position="identity",breaks=as.numeric(floor(min(genData())):ceiling(max(genData())))) +
          geom_vline(xintercept= input$cutoff, colour="red", size=1) + ggtitle("Histogram of Simulated Half-Lives")+
          theme(plot.title= element_text(face="bold")) +
          scale_x_continuous(breaks=as.numeric(floor(min(genData())):ceiling(max(genData()))))
      }
    }
    
    
    
  })
}