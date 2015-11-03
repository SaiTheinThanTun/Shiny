library(shiny)
library(boot)
library(MASS)
library(segmented)
library(sm)
library(mixtools)
library(vioplot)

ui <- fluidPage(
  titlePanel("Mixture modelling of Plasmodium falciparum half lives"),
  fileInput(inputId = "file", label = "Select your input file: (simulated_cloneHLdata_SMRUbyyear.csv in this case)"),
  plotOutput(outputId = "histoplot"),
 
  
  textOutput(outputId = "proportions"),
  textOutput(outputId = "results")
  )

server <- function(input, output) {
        
    output$histoplot <- renderPlot({
    inFile <- input$file
    mixdat <- read.csv(inFile$datapath)
    
    N<-ncol(mixdat)-1
    M <- 5
    
    
    
    
    
    pval<-0.1
    nboot<-100 # number of iterations for bootstrap
    nsim<-1000 # number of iterations for creating probaility of resistance vs HL graphs
    P<-2 # use P or more samples to get geometric means and discard all other samples for permutation analysis
    T<-100 # number of permutations for permutation analysis
    smax=5000
    
    # create output matrices
    output.mu <- matrix(NA,nrow=M,ncol=N)
    output.sigma <- matrix(NA,nrow=M,ncol=N)
    output.lambda <- matrix(NA,nrow=M,ncol=N)
    output.loglik <- matrix(NA,nrow=M,ncol=N)
    output.mu.se <- matrix(NA,nrow=M,ncol=N)
    output.sigma.se <- matrix(NA,nrow=M,ncol=N)
    output.lambda.se <- matrix(NA,nrow=M,ncol=N)
    AIC<-matrix(0,nrow=M,ncol=N)
    AICdelta<-matrix(0,nrow=M,ncol=N)
    
    nb<-na.omit(mixdat[,N+1])
    
    # fit single component model
    for (i in 1:N){
            # 1 COMPONENT LOG NORMAL
            nmixdat<-na.omit(mixdat[,i])
            lmixdat<- log(nmixdat)
            xll<-fitdistr(lmixdat,"normal")
            output.loglik[1,i]<- xll$loglik
            output.mu[1,i]<-xll$estimate[1]
            output.lambda[1,i]<-1
            output.sigma[1,i]<-xll$estimate[2]
            output.mu.se[1,i]<-xll$sd[1]
            output.sigma.se[1,i]<-xll$sd[2]
            output.lambda.se[1,i]<-0
            AIC[1,i]<-2*(3*1-1)-2*output.loglik[1,i]
            AICdelta[1,i]<-0
    }
    
    # fit multiple component models sequntially
    for (i in 1:N){
      nmixdat<-na.omit(mixdat[,i])
      lmixdat<- log(nmixdat)
      # >=2 COMPONENTS LOG NORMAL
      j<-1
      # stop if j-component model is more parsimonious than (j-1)-compnent model
      while((j<=M-1) && AICdelta[j,i]<=pval){
        j<-j+1
        res <- normalmixEM(lmixdat, lambda = matrix((1/j),nrow=1,ncol=j), mu = 2*(1:j)/j, sigma = 0.3*matrix(1,nrow=1,ncol=j))
        resboot <- boot.se(res, B = nboot)
        resboot[c("lambda.se", "mu.se", "sigma.se","loglik.se")]	
        output.loglik[j,i]<-res$loglik
        AIC[j,i]<-2*(3*j-1)-2*output.loglik[j,i]
        AICdelta[j,i]<-exp(-(AIC[j-1,i]-AIC[j,i])/2)
        if(AICdelta[j,i]<=pval){
          output.mu[1:j,i]<-res$mu
          output.sigma[1:j,i]<-res$sigma
          output.lambda[1:j,i]<-res$lambda
          output.mu.se[1:j,i]<-resboot$mu.se
          output.sigma.se[1:j,i]<-resboot$sigma.se
          output.lambda.se[1:j,i]<-resboot$lambda.se		
          
          }
      }
    }
    
    Sys.sleep(0.02)
    for (ds in 1:N){
            Sys.sleep(0.02)
            nmixdat<-na.omit(mixdat[,ds])
            plam<-na.omit(output.lambda[,ds])
            pmu<-na.omit(output.mu[,ds])
            psig<-na.omit(output.sigma[,ds])
            hist(nmixdat,freq=FALSE,main = paste("Northwestern Thai-Myanmar border",2000+ds),xlab = "Clearance half-life (hours)",ylim=c(0,0.6),col="grey",lwd=2,ps=20,breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12))
            x <- seq(0.1, max(nmixdat), length=1000)
            hx<-plam[1]*dlnorm(x,meanlog=(pmu[1]),sdlog=psig[1])
            if(length(plam)>1){
                    for(k in 2:length(plam)){
                            hx<-hx+plam[k]*dlnorm(x,meanlog=(pmu[k]),sdlog=psig[k])
                    }
        }
            lines(x,hx,col="red", lwd=5)
            Sys.sleep(0.02)
    }
    
  })

  output$proportions <- renderText({
          M
  })
  
  output$results <- renderText({"is the estimated proportion of patients with artemesinin resistance"})
  


  
  
       
 
  }

shinyApp(ui = ui, server = server)
