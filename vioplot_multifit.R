library(shiny)
library(boot)
library(MASS)
library(segmented)
library(sm)
library(mixtools)
library(vioplot)

ui <- fluidPage(
  fileInput(inputId = "file", label = "Select your input file: (simulated_cloneHLdata_SMRUbyyear.csv in this case)"),
  plotOutput(outputId = "vioplot")
  )

server <- function(input, output) {
  output$vioplot <- renderPlot({
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
    
    vioplot(na.omit(mixdat[,1]),na.omit(mixdat[,2]), na.omit(mixdat[,3]), na.omit(mixdat[,4]), na.omit(mixdat[,5]), na.omit(mixdat[,6]), na.omit(mixdat[,7]),na.omit(mixdat[,8]), na.omit(mixdat[,9]), na.omit(mixdat[,10]), na.omit(mixdat[,11]), na.omit(mixdat[,12]),col="grey",pchMed=30)
    
    mev1<-matrix(1,nrow=N*M,ncol=1)
    mev2<-matrix(1,nrow=N*M,ncol=1)
    mev3<-matrix(1,nrow=N*M,ncol=1)
    for (a in 1:N){
      for (b in 1:M){
        mev1[M*(a-1)+b]<-2000+a+b/1000
        mev2[M*(a-1)+b]<-exp(output.mu[b,a])
        mev3[M*(a-1)+b]<-(output.lambda[b,a])^0.5
      }
    }
    dfx = data.frame(ev1=mev1, ev2=mev2, ev3=mev3)
    symbols(x=dfx$ev1-2000, y=dfx$ev2, circles=dfx$ev3, inches=1/7, ann=F, bg="red", fg=NULL,xaxp = c(2001, 2015, 14),ylim=c(2,8),add=TRUE)
    title(xlab = "Time (years)",ylab="Clearance half-life (hours)",ps=20)
    
  })
}

shinyApp(ui = ui, server = server)
