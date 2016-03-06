library(deSolve)
library(shiny)

ui <- fluidPage(
  fluidRow(
    column(4,
           h3("Transmission parameters"),
           sliderInput(inputId = "mui",
                       label = "birth rate of mosquitos",
                       value = .05, min = .001, max = .1 #10 days survival= 20 half-days survival, therefore 1/20=.05
           ),
           sliderInput(inputId = "muo",
                       label = "death rate of mosquitos",
                       value = .05, min = .001, max = .1
           ),
           sliderInput(inputId = "b",
                       label = "biting rate",
                       value = .1, min = .01, max = .3
           ),
           sliderInput(inputId= "beta",
                       label = "probability of disease transmission per bite for mosquitos",
                       value = .3, min = .01, max = 1
           ),
           sliderInput(inputId= "beta_h",
                       label = "probability of disease transmission per bite for humans",
                       value = .3, min = .01, max = 1
           )
           ),
    column(4,
           h3("Population parameters"),
           numericInput(inputId = "initPh",
                        label = "Human Population",
                        value = 80
           ),
           numericInput(inputId = "initIh",
                        label = "Infected Human Population",
                        value = 30
           ),
           numericInput(inputId = "initP",
                        label = "Initial Mosquito Population",
                        value = 800
           ),
           sliderInput(inputId = "initI",
                       label = "Initial infected Mosquito Population",
                       value = 200, min = 1, max=800
           ),
           numericInput(inputId= "weeks",
                        label= "No. of weeks for the model",
                        value = 52
                        )
           )
  ),
  fluidRow(
    #column(5,
    #       h3("plot of inc2"),
    #       plotOutput(outputId = "inc2")
    #       ),
    column(5,
           h3("human pop"),
           plotOutput(outputId="human_pop")
           ),
    column(5,
           h3("plot of everything_mosq"),
           plotOutput(outputId = "everything_mosq")
           ),
    column(1,h3("lam"),
           textOutput(outputId = "lam")
           )
  )
  
  
)

server <- function(input, output) {
  parameters_r <- reactive({
    c(mui=input$mui,    # birth #lifespan of mosquito 10 days
      muo=input$muo,    # death
      #beta= #per capita effective contact with infected human per unit time
      #ce = (.3*.01), #probability of disease transmission per bite * biting rate
      b = input$b, #biting rate
      beta = input$beta, #probability of disease transmission per bite for mosquitos
      beta_h = input$beta_h #probability of disease transmission per bite for human
        )})
  
  output$lam <- renderPrint({
    parameter <- parameters_r()
    lam <- parameter[3]*parameter[4]*(input$initIh/input$initPh)
    lam
  })
  
  ode_out <- reactive({
    # define the number of weeks to run the model
    times <- seq(0, input$weeks, by = (1/14)) #previously week was 52, 7 days/14 = 1/2 day timestep
    
    #MODEL PARAMETERS
    parameters <- parameters_r()
    
    # MODEL INITIAL CONDITIONS
    initPh <- input$initPh #total human population
    initIh <- input$initIh #no of infected people, ideally this should be changing either from the IBM or the ODE model itself
    #page 75 of A biologist's guide to mathematical modelling
    initRh <- 0 #no of
    initSh <- initPh-(initIh+initRh) #no of suscepitables
    
    initP<- input$initP #209100 
    initI<-input$initI
    initS<-initP-initI
    initD <- 0
    
    state <- c(S = initS, I = initI, D = initD, P = initP, Sh= initSh, Ih = initIh, Rh= initRh, Y=0)
    
    
    # set up a function to solve the model
    mosQ<-function(t, state, parameters) 
    {
      with(as.list(c(state, parameters)),
           {
             
             # define variables
             P <- (S+I)
             Ph <- (Sh+Ih+Rh)
             #seas<-1+amp*cos(2*pi*(Y-phi)/52)
             #beta<-R0*(muo+nui)*gamma/(muo+gamma)
             lam <- b*beta*(Ih/Ph)
             lam_h <- b*beta_h*(I/P)
             
             # rate of change for mosquitos
             dS <- mui*P-muo*S-lam*S
             dI <- -muo*I+lam*S
             dD <- muo*S+muo*I
             dP <- 0
             
             # rate of change for humans
             dSh <- -lam_h*Sh
             dIh <- lam_h*Sh
             dRh <- 0
             dY <- 1
             
             # return the rate of change
             list(c(dS, dI, dD,dP, dSh, dIh, dRh, dY))
           }
      ) 
      
    }
    ode(y = state, times = times, func = mosQ, parms = parameters)
    
    #out <- ode(y = state, times = times, func = mosQ, parms = parameters)
    #out[,]
  })
  

  
  output$everything_mosq <- renderPlot({
    out <- ode_out()
    
    plot(out[,1],out[,2], type="l", col="blue", axes=FALSE, xlab="", ylab="", main="mosq_pop")
    axis(2, ylim=c(0,17),col="blue") 
    mtext("Susceptible mosquitoes",side=2,line=2.5) 
    box()
    par(new=TRUE)
    plot(out[,1],out[,3], type="l", col="red", axes=FALSE, xlab="", ylab="")
    axis(4, ylim=c(0,17),col="red") 
    mtext("Infected mosquitoes",side=4,line=2.5)
    
    axis(1,pretty(range(out[,1]),10))
    mtext("Time (Weeks)",side=1,col="black",line=2.5)
    
    legend("top",legend=c("Susceptibles","Infected"),
           text.col=c("blue","red"),pch= "__", col=c("blue","red"))
  })
  
  output$human_pop <- renderPlot({
    out <- ode_out()
    
    plot(out[,1],out[,6], type="l", col="blue", axes=FALSE, xlab="", ylab="", main="human_pop")
    axis(2, ylim=c(0,17),col="blue") 
    mtext("Susceptible humans",side=2,line=2.5) 
    box()
    par(new=TRUE)
    plot(out[,1],out[,7], type="l", col="red", axes=FALSE, xlab="", ylab="")
    axis(4, ylim=c(0,17),col="red") 
    mtext("Infected humans",side=4,line=2.5)
    
    axis(1,pretty(range(out[,1]),10))
    mtext("Time (Weeks)",side=1,col="black",line=2.5)
    
    legend("top",legend=c("Susceptibles","Infected"),
           text.col=c("blue","red"),pch= "__", col=c("blue","red"))
  })
}

shinyApp(ui = ui, server = server)

#output$inc2 <- renderPlot({
#  parameters <- parameters_r()
#  out <- ode_out()
#  #pop <- out[,2]+out[,3]
  
#  #inc <- parameters[3]*parameters[4]*(out[,7]/out[,5])*out[,3]*out[,2]/pop
#  inc2 <- parameters[3]*parameters[4]*(out[,7]/(out[,6]+out[,7]))*out[,3]
#  #plot(inc2, type="l")
#  plot(out[,3])
#})