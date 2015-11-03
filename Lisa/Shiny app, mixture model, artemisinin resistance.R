library(shiny)
library(boot)
library(MASS)
library(segmented)
library(sm)
library(mixtools)
library(vioplot)

ui <- fluidPage(
        titlePanel("Mixture modelling of Plasmodium falciparum half lives"),
        br(),
        br(),
        p("This shiny app allows you to use a model to help identify the presence of artemesinin resistance
          amongst patients with falciparum malaria. The model used is described in", span(a("White et al. 2015", href = "http://journals.plos.org/plosmedicine/article?id=10.1371/journal.pmed.1001823"),tags$sup(1)),". At 
          at the bottom of this webpage can be found some limitations of the applicability of the model."),
        p("The model takes, as its input, the half lives of falciparum malaria parasites in the blood of patients after 
          their starting artemisinin combination therapy. If you click the button below, you will be able to upload
          half lives from your patients of interest. Please have the half lives in a single column of an excel spread 
          sheet, with one half life recorded for each patient's infection."),
        br(),
        br(),
        fileInput(inputId = "file", label = "Select your input file: (simulated_cloneHLdata_SMRUbyyear.csv in this case)"),
        br(),
        br(),
        textOutput(outputId = "explanation1"),
        br(),
        textOutput(outputId = "explanation2"),
        br(),
        
        fluidRow(
                column(width = 6,
                       plotOutput(outputId = "histoplot1")
                ),
                column(width = 6,
                       plotOutput(outputId = "histoplot2")
                ) 
        ),
        
        fluidPage(
                br(),
                textOutput(outputId = "explanation3"),
                br(),
                verbatimTextOutput(outputId = "geometric_means_and_proportions"),
                tableOutput(outputId = "full_table"),
                textOutput(outputId = "geometric_means"),
                textOutput(outputId = "results2"),
                textOutput(outputId = "proportions"),
                textOutput(outputId = "conclusions"),
                br(),
                br(),
                strong("Parasite clearance half-life for an indiviual patient"),
                br(),
                br(),
                p("You can enter the half-life of an individual patient into the box below to determine the probabilities
                  of their belonging to each of the respective component distributions of the data set entered above. Note that
                  if the individual patient is from a different population, such as in terms of age, or location, it may not
                  be useful to determine these probabilities."),
                textInput("INPUT", label = NULL, value = "Enter half-life..."),
                verbatimTextOutput(outputId = "Prob_of_resistance_in_given_patient"),
                br(),
                br(),
                strong("Limitations of the applicability of the model"),
                br(),
                br(),
                p("This model was generated in the low transmission setting of South East Asia. It may not
                  be applicable to high transmission settings, such as in parts of Sub-Saharan Africa. Furthermore, 
                  multi clonal infections were excluded from the data used to generate the model, therefore, the model
                  may not be as applicable to users including multiclonal infections in their data, particularly
                  in high transmission settings where multiclonal infections are more prevalent."),
                br(),
                br(),
                strong("Reference"),
                br(),
                br(),
                em("1 White LJ, Flegg JA, Phyo AP, Wiladpai-ngern JH, Bethell D, Plowe C, et al. Defining the In Vivo Phenotype of Artemisinin-Resistant Falciparum Malaria: A Modelling Approach. 2015 [cited 2015 Jul 1]; Available from: http://dx.plos.org/10.1371/journal.pmed.1001823"),
                br(),
                br(),
                br()
                ))

server <- function(input, output) {
        
        #Generate the histogram (via barplot) of the aggregate data in White et al. 2015. 
        
        output$histoplot1 <- renderPlot({
                
                xval = c(0:12)
                count = c(0.000656598818122127, 0.0899540380827314, 0.278069599474721, 0.192711753118844, 0.137229152987525, 0.10177281680893, 0.0971766250820748, 0.0604070912672357, 0.0275771503611294, 0.0124753775443204, 0.00196979645436638, 0)
                df = as.data.frame(cbind(xval, count))
                mp <- barplot(df$count, space=0,ylim=c(0,0.6), xlim=c(0,12), main = paste("Distribution of parasite clearance half lives","\n", "from White et al. 2015"), xlab = "Clearance half-life (hours)", ylab = "Density", axes=FALSE) 
                axis(side=2, pos=-0.5)
                axis(side=1, at =mp-0.5, labels=df$xval)
                
                # Draw red line on histogram. Need to get values for the data and to increase the number of intervals used to make the curve smooth. 
                x<-c((1:80)*12/80);         
                emp_dens <- c(2.02585057383938e-14, 7.77173757739607e-09, 2.78940113823811e-06, 8.70671331737266e-05, 0.000825882146448662, 0.00395529148563944, 0.0122884693613028, 0.0284849695997866, 0.0536290713696855, 0.0866277262154361, 0.124595870701594, 0.163783977451365, 0.200526991552525, 0.231906529288011, 0.25605189599341, 0.272146977880007, 0.280257646559789, 0.281084902536158, 0.275717402391711, 0.265424150506261, 0.251502660144673, 0.235181347692743, 0.217565936906641, 0.199616244680149, 0.182140174875239, 0.165794643739927, 0.151087232901058, 0.138376458928438, 0.127871731987167, 0.119635803880625, 0.113592694813787, 0.109543062351541, 0.107187299917257, 0.106154927937539, 0.106037540950885, 0.106421952743341, 0.106920268337932, 0.107194265142465, 0.106972454104453, 0.106059270825135, 0.10433681498884, 0.101760284689963, 0.0983486865381704, 0.0941725519498531, 0.0893403059939812, 0.0839846890946096, 0.0782502977104948, 0.072282951899071, 0.0662212633520916, 0.0601904972848442, 0.0542986097015175, 0.0486341997535655, 0.0432660383494602, 0.0382438072817113, 0.0335996945757871, 0.0293505285970653, 0.0255001844494674, 0.0220420524749678, 0.0189614137894603, 0.0162376175975855, 0.0138459972260963, 0.0117594955492012, 0.00994999590669894, 0.00838937253017513, 0.00705028602685715, 0.00590675587678106, 0.0049345443980872, 0.00411138632708848, 0.00341709597563706, 0.00283358061704542, 0.00234478488623254, 0.00193658697570432, 0.00159666355009287, 0.00131433676602087, 0.001080413667746, 0.000887025570021812, 0.000727472834767737, 0.000596078668184522, 0.000488054165556397, 0.000399375760908199);
                dat_dens <- as.matrix(cbind(x,emp_dens));
                lines(dat_dens,lty=1,col="red",lwd=5);
        })
        
        #For the users data, run the mixture model and draw the histogram. 
        
        output$histoplot2 <- renderPlot({
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
                                        
                                        output$explanation1 <- renderText({"Below are two graphs. The graph on the left
                                                represents aggregate data from White et al. 2015. The graph on the right is a graph made from your data."})
                                        
                                        output$explanation2 <- renderText({"The graph on th left depicts two half life distributions
                                                with geometric means SOMETHING AND SOMETHING ELSE respectively. 
                                                The distribution with a geometric mean half life of SOMETHING was intepreted
                                                as representing patients with parasites sensitive to artemisinin. The distribution
                                                with a geometric mean half life of SOMETHING ELSE was interpreted as representing
                                                patients with parasites resistant artemisinin. With this information, you may be able
                                                to interpret the graph on the rigth, which represents your own data."})
                                        
                                        output$explanation3 <- renderText({"Below are some statistics from 
                                                from the graph representing your data. There is also a list of the probabilities of 
                                                each patient belonging to each of the component distributions depicted in the graph."})
                                        
                                        output$geometric_means_and_proportions <- renderPrint({
                                                
                                                cat("The model predicts ", j, " component geometric mean half lives (hours):",
                                                    
                                                    "\n\n")
                                                
                                                for (a in 1:j) {
                                                        cat("Distribution",a,"\n",
                                                            
                                                            "Geometric mean = ", exp(res$mu[a]),
                                                            
                                                            "\n", 
                                                            
                                                            "Contribution to composite distribution = ", res$lambda[a],
                                                            
                                                            "\n\n"
                                                            
                                                        )       
                                                }
                                        })
                                        
                                        
                                        
                                        }
                                
                                
                                Sys.sleep(0.02)
                                for (ds in 1:N){
                                        Sys.sleep(0.02)
                                        nmixdat<-na.omit(mixdat[,ds])
                                        plam<-na.omit(output.lambda[,ds])
                                        pmu<-na.omit(output.mu[,ds])
                                        psig<-na.omit(output.sigma[,ds])
                                        hist(nmixdat,freq=FALSE,main = paste("Distribution of parasite clearance half lives","\n", "from your data"),xlab = "Clearance half-life (hours)",ylim=c(0,0.6),col="grey",lwd=2,ps=20,breaks=c(0,1,2,3,4,5,6,7,8,9,10,11,12))
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
                                
                                        }
                        
                        }
                
                #Make the table containing the probabilities of each patient belonging to each component distribution. 
                
                output$full_table <- renderTable({ 
                        
                        nmixdat<-na.omit(mixdat)
                        
                        prob <- matrix(,nrow = length((nmixdat)[[1]]),ncol = length(res$mu))                      
                        
                        relative_prob <- matrix(,nrow(prob), ncol(prob))   
                        
                        for (h in 1:(length((nmixdat)[[1]]))){
                                
                                for(d in 1:length(res$mu)){
                                        
                                        prob[h,d] <-  ((dnorm((nmixdat[[1]][h]), exp(res$mu)[d], sd = exp(res$sigma)[d]))*res$lambda[d])       
                                        
                                        
                                        for (i in 1:nrow(prob)){
                                                
                                                row_total <- sum(prob[i,])
                                                
                                                for (a in 1:ncol(prob)){
                                                        
                                                        relative_prob[i,a] <- (prob[i,a]/row_total)
                                                        
                                                }
                                                
                                        }        
                                }
                        }
                        
                        
                        rounded_relative_prob <- (round((relative_prob), digits = 2))
                        
                        A <- as.character((1:(nrow(rounded_relative_prob))))
                        B <- rep("Patient", times = (nrow(rounded_relative_prob)))
                        H <- (cbind(B,A))
                        rownames(rounded_relative_prob) <- str_c(H[,1],H[,2], sep = " ")
                        
                        G <- as.character((1:(ncol(rounded_relative_prob))))
                        L <- rep("Probability of being in distribution", times = (ncol(rounded_relative_prob)))
                        W <- (cbind(L,G))
                        colnames(rounded_relative_prob) <- str_c(W[,1],W[,2], sep = " ")
                        
                        print(rounded_relative_prob)
                        
                })      
                
                
                #Calculate the probabilities of an individual patient, typed into the box by the user, belong to each component distribution.
                # Print the results. 
                
                means <- res$mu
                spreads <- exp(res$sigma)
                proportions <- res$lambda
                
                output$Prob_of_resistance_in_given_patient <- renderPrint({
                        
                        absolute_probabilities <- NA
                        
                        for (i in 1:length(means)) {
                                
                                absolute_probabilities[i] <-  ((dnorm((as.numeric(input$INPUT)), mean = means[i], sd = spreads[i]))*proportions[i])       
                                
                        }
                        
                        probabilities <- (absolute_probabilities /(sum(absolute_probabilities)))
                        
                        for (x in 1:length(probabilities)){
                                
                                cat('The probability of this patient being from distribution', x, 'is', round(probabilities[x], digits = 3),
                                    
                                    "\n\n"
                                )
                        }
                })
})
        
        }

shinyApp(ui = ui, server = server)