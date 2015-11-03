# Here we're trying to print the probabilities of each individual belonging to each distribution all in one go. 

library(shiny)
library(boot)
library(MASS)
library(segmented)
library(sm)
library(mixtools)
library(vioplot)

ui <- fluidPage(
        titlePanel("Mixture modelling of Plasmodium falciparum half lives"),
        
        p("This shiny app allows you to use a model to help identify the presence of artemesinin resistance
          amongst patients with falciparum malaria. The model uses is described in White et al. 2015 PUT IN REFERENCE. 
          The model takes, as its input, the half lives of falciparum malaria parasites in the blood or patients after 
          their starting artemisinin combination therapy. If you click the button below, you will be able to upload
          half lives from your patients of interest. Please have the halve lives in a single column of an excel spread 
          sheet, with one half life recorded for each patient's infection."),
        
        fileInput(inputId = "file", label = "Select your input file: (simulated_cloneHLdata_SMRUbyyear.csv in this case)"),
        plotOutput(outputId = "histoplot"),
        textOutput(outputId = "explanation1"),
        br(),
        textOutput(outputId = "explanation2"),
        br(),
        textOutput(outputId = "explanation3"),
        br(),
        textOutput(outputId = "results1"),
        verbatimTextOutput(outputId = "geometric_means_and_proportions"),
        textOutput(outputId = "full_table"),
        textOutput(outputId = "geometric_means"),
        textOutput(outputId = "results2"),
        textOutput(outputId = "proportions"),
        textOutput(outputId = "conclusions"),
        textInput("text", label = h3("Parasite clearance half life for an individual patient"), value = "Enter text..."),
        textOutput(outputId = "Prob_of_resistance_in_given_patient")
        
        )

server <- function(input, output) {
        
        output$histoplot <- renderPlot({
                inFile <- input$file
                mixdat <- read.csv(inFile$datapath)
                
                N<-ncol(mixdat)
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
                
                nb<-na.omit(mixdat[,N])
                
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
                
                # fit multiple component models sequentially
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
                                        
                                        
                                        output$explanation1 <- renderText({"Below are two graphs. The one on the left 
                                                represents aggregate data from White et al. 2015. On the right is a graph made
                                                from your data."})
                                        
                                        output$explanation2 <- renderText({"The graph on the left depicts two distributions of
                                                half lives of geometric means SOMETHING AND SOMETHING ELSE respecitvely. 
                                                The distribution with geometric mean half life of SOMETHING was intepreted
                                                as representing patients with parasites sensitive to artemisinin. The distribution
                                                on with geometric mean half life of SOMETHING ELSE was interpreted as representing
                                                patients with parasites resistant artemisinin. With this information, you may be abline
                                                to interpret the graph on the right, which represents your own data."})
                                        
                                        output$explanation3 <- renderText({"Below the graphs you can see some statistics from 
                                                the graph representing your data. There is also a list of the probabilities of 
                                                each patient belonging to each of the component distributions depicted in the graph."})
                                        
                                        output$results1 <- renderText({
                                                
                                                paste0("The model predicts ", j, " component geometric mean half lives (hours):")
                                                
                                        })
                                        
                                        
                                        
                                        
                                        
                                        output$geometric_means_and_proportions <- renderPrint({
                                                for (a in 1:j) {
                                                        cat("Distribution",a,"\n",
                                                            
                                                            "Geometric mean = ", res$mu[a],
                                                            
                                                            "\n", 
                                                            
                                                            "Contribution to composite distribution = ", res$lambda[a],
                                                            
                                                            "\n"
                                                            
                                                        )
                                                }
                                                
                                                
                                        })
                                        
                                        
                                        
                                        nmixdat<-na.omit(mixdat)
                                        
                                        prob <- matrix(nrow = length(nmixdat),ncol = length(res$mu))
                                        
                                        output$full_table <- renderPrint({
                                                
                                                
                                                        
                                                        for(d in 1:length(res$mu)){
                                                          for (h in 1:(length(nmixdat))){
                                                                prob[h,d] <-  ((dnorm((nmixdat[h]), mean = exp(res$mu)[d], sd = exp(res$sigma)[d]))*res$lambda[d])       
                                                                
                                                        }
                                                        
                                                }
                                                
                                                print(prob) 
                                                
                                        })
                                        
                                        
                                        
                                        
                                        
                                        
                                        
                                        output$geometric_means <- renderText({
                                                
                                                round(exp(res$mu), digits = 2) 
                                                
                                        })
                                        
                                        
                                        output$results2 <- renderText({"Respective weightings of each component distribution:"})
                                        
                                        
                                        output$proportions <- renderText({
                                                
                                                round(res$lambda, digits = 2)  
                                        })
                                        
                                        
                                        if ((length(res$mu)>=2) & (sum(exp(res$mu)>=6.5)>=1))
                                                
                                        {
                                                a <- "is present, and that the proportion of patients with artemisinin resistant parasites is "
                                                
                                        }
                                        else {
                                                a <- "is not present."   
                                                
                                        }
                                        
                                        output$conclusions <- renderText({
                                                
                                                paste0("The model predicts that artemisinin resistance ", a, round(exp(res$lambda[3]), digits = 2) 
                                                )
                                                
                                        })
                                        
                                        output$value <- renderPrint({ input$text })
                                        
                                        #Find the index of the component distribution with geometric mean closest to 6.5. Assume this is the one corresponding to resistance.
                                        calculate_t_half_6.5ish <- (order((res$mu-6.5)^2))
                                        index <- which(calculate_t_half_6.5ish == 1)[[1]]
                                        t_half_6.5ish <- res$mu[index]
                                        
                                        #Given an the value of an individual patient's half life typed into the shiny page, what is the probability that it occurs in the 'resistant distribution'? Ignore distributions with half lives longer than the 'resistant distribution'.
                                        
                                        unconditional_prob_res <- ((dnorm((as.numeric(input$text)), mean = exp(res$mu)[index], sd = exp(res$sigma)[index]))*res$lambda[index])
                                        
                                        vector_unconditional_prob_not_res <- NA
                                        
                                        for (i in 1:(index-1)) {
                                                
                                                vector_unconditional_prob_not_res[i] <-  ((dnorm((as.numeric(input$text)), mean = exp(res$mu)[i], sd = exp(res$sigma)[i]))*res$lambda[i])       
                                                
                                        }
                                        unconditional_prob_not_res <- sum(vector_unconditional_prob_not_res)
                                        
                                        
                                        conditional_prob_res <- (unconditional_prob_res /(unconditional_prob_res + unconditional_prob_not_res))
                                        
                                        output$Prob_of_resistance_in_given_patient <- renderText({
                                                
                                                
                                                paste0("The model predicts that the probability of this individual having resistant parasites is ", round(conditional_prob_res, digits = 5) 
                                                )
                                                
                                                
                                                
                                                
                                                
                                                
                                        })
                                        
                                        
                                        
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
                        hist(nmixdat,freq=FALSE,main = paste("Distribution of parasite clearance half lives"),xlab = "Clearance half-life (hours)",ylim=c(0,0.6),col="grey",lwd=2,ps=20,breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12))
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
        
        }

shinyApp(ui = ui, server = server)