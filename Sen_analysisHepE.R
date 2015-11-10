# Practical session II: Economic evaluation of latrine to prevent hepatitis E virus
# Cost-Utility analysis

library(deSolve)
library(shiny)

ui <- fluidPage(
  titlePanel("Pratical session 2: Sensitivity analysis/ EE of latrine to prevent hep E virus"),
  fileInput(inputId="file", label= "Input hepEdata.csv file:"),
  sliderInput(inputId="l_eff", label= "Effectiveness of latrine", min=0,max=1, step=.1, value=.7),
  sliderInput(inputId="C_latrine", label="Cost of latrine per unit $US",min=20,max=400, step=50, value=200),
  sliderInput(inputId="R0", label="Basic reproduction number",min=1,max=5,step=.1, value=3.5),
  tableOutput(outputId="CEA")
)

server <- function(input, output){
  output$CEA <- renderTable({
    
    inFile <- input$file
    cases <- read.csv(inFile$datapath)
    #cases <- read.csv('hepEdata.csv')
    # define number of weeks as the length of the first column of the data matrix
    nw<-length(cases[,1])
    cases<-cases[1:nw-1,]
    nw<-nw-1
    # create the calendar dates from the weeks
    dates_data<-as.Date('31/12/2012',format='%d/%m/%Y')+7*(38+(1:nw))
    
    # we know the epidemic was over by May 2014 so we include one extra data point
    cases<-rbind(cases,c(20,0))
    dates_data<-c(dates_data,as.Date('31/05/2014',format='%d/%m/%Y'))
    
    # define the number of weeks to run the model
    times <- seq(0, 52, by = (1/7))
    
    ############################
    #Basecase (no intervention)
    ############################
    #MODEL PARAMETERS
    base_parameters <- c(mui=(1/(50*52)),    # birth
                         muo=(1/(50*52)),    # death
                         R0=input$R0,               # basic reproduction number
                         omega=(1/(10*52)),  # rate of loss of immunity = 1/(average duration of immunity)
                         gamma=1/4,          # rate of movement from latent to infectious stage = 1/(average latent period)
                         nui=1/4,            # rate of recovery = 1/(average duration of infection)
                         report=1/7,         # proportion of all infections that are reported
                         amp=0,              # relative amplitude of seasonal forcing
                         phi=0,              # week of peak in seasonal forcing
                         lagepiweeks=13,     # number of weeks before the first case report that the outbreak began
                         week_interv=15,     # number of weeks after first case report when the intervention begins
                         l_cov=0,            # coverage of latrines (reported as 0.24)
                         l_eff=0.7           # redcution in transmission given full latrine coverage (reports have suggested that poor hygeine triples the risk of infection)
    )
    
    ############################
    #Latrines (Intervention)
    ############################
    #MODEL PARAMETERS
    lat_parameters <- c(mui=(1/(50*52)),    # birth
                        muo=(1/(50*52)),    # death
                        R0=input$R0,               # basic reproduction number
                        omega=(1/(10*52)),  # rate of loss of immunity = 1/(average duration of immunity)
                        gamma=1/4,          # rate of movement from latent to infectious stage = 1/(average latent period)
                        nui=1/4,            # rate of recovery = 1/(average duration of infection)
                        report=1/7,         # proportion of all infections that are reported
                        amp=0,              # relative amplitude of seasonal forcing
                        phi=0,              # week of peak in seasonal forcing
                        lagepiweeks=13,     # number of weeks before the first case report that the outbreak began
                        week_interv=15,     # number of weeks after first case report when the intervention begins
                        l_cov=0.7,          # coverage of latrines
                        l_eff=input$l_eff           # redcution in transmission given full latrine coverage (reports have suggested that poor hygeine triples the risk of infection)
    )
    
    
    ### Costs ###
    cost_parameters <- c(
      C_latrine <- input$C_latrine,                   # Cost of latrine per unit ($US)
      Num_household <- 7000/5,            # Number of household in the population (average 5 people per household)
      C_hospital <- 30,                   # Cost of hospitalisation in case of fulminant per day ($US)
      hosp_day <- 14,                     # number of hospitalisation day per fulminant case (days)
      l_cov<-0.7                          # Latrine coverage
    )
    
    ## At the sensitivity analysis section of the practical remove the "#" from the below two lines ##
    
    CEA_fun <- function(base_parameters, lat_parameters, 
                        cost_parameters){
      
      # MODEL INITIAL CONDITIONS
      initP<-7000
      initE<-1
      initI<-0
      initR<-0
      initS<-initP-initE-initI-initR
      
      state <- c(S = initS, E=initE, I = initI,R = initR, Y=0)
      times <- seq(0, 52, by = (1/7))
      
      # set up a function to solve the model
      HepE<-function(t, state, parms) 
      {
        with(as.list(c(state, parms)),
             {
               
               # define intervention variables
               latrine<-(Y<=(week_interv+lagepiweeks))+(1-l_cov*l_eff)*(Y>(week_interv+lagepiweeks))
               
               # define variables
               P <- (S+E+I+R)
               seas<-1+amp*cos(2*pi*(Y-phi)/52)
               beta<-R0*(muo+nui)*gamma/(muo+gamma)
               # include redcution in transmission due to latrines
               lam <- latrine*beta*seas*I/P
               
               # rate of change
               dS <- mui*P-muo*S-lam*S+omega*R
               dE <- -muo*E+lam*S-gamma*E
               dI <- -muo*I+gamma*E-nui*I
               dR <- -muo*R+nui*I-omega*R
               dY <- 1
               
               # return the rate of change
               list(c(dS, dE, dI, dR, dY))
               
             }
        ) 
        
      }
      
      #Basecase Output
      out_base <- ode(y = state, times = times, func = HepE, parms = base_parameters)
      pop_base <-out_base[,2]+out_base[,3]+out_base[,4]+out_base[,5] #to check model output
      
      #define the dates for the model output
      time<-out_base[,1]
      dates_model<-dates_data[1]+(time)*7-base_parameters[10]*7
      
      #Total cases from model prediction
      cumulativereports_base<- base_parameters[7]*out_base[,5]
      
      Total_case_base <- cumulativereports_base[365]
      
      #Intervention Output
      out_latrine <- ode(y = state, times = times, func = HepE, parms = lat_parameters)
      pop_latrine <-out_latrine[,2]+out_latrine[,3]+out_latrine[,4]+out_latrine[,5] #to check model output
      
      # define the dates for the model output
      time<-out_latrine[,1]
      dates_model<-dates_data[1]+(time)*7-lat_parameters[10]*7
      
      #Total cases from model prediction
      cumulativereports_latrine<- lat_parameters[7]*out_latrine[,5]
      Total_case_latrine <- cumulativereports_latrine[365]
      
      ############################
      # Other relevant outcomes
      ############################
      #Convert outputs to from the transmission model to use in economic analysis
      
      #Parameters; no of HepE cases, no of acute jaundice, no of death
      risk_fulminant <- 0.05             #Risk of fulminant amongst Hepatitis E virus cases (about 5% of HepE)
      risk_death_fulminant <- 0.08       #Risk of death amongst Hepatitis E cases with fulminant (about 8% of fulminant)
      risk_mildill <- 1-risk_fulminant
      
      num_fulminant_base<-Total_case_base*risk_fulminant           #number of fulminant hep cases = number of cases*risk of fulminant
      num_fulminant_latrine<-Total_case_latrine*risk_fulminant
      
      num_mildill_base<-Total_case_base*risk_mildill               #number of mild illness (HepE case without fulminant hepatitis) = number of cases*(1-risk of fulminant)
      num_mildill_latrine<-Total_case_latrine*risk_mildill
      
      num_death_base<-num_fulminant_base*risk_death_fulminant      #number of death due to fulminant hepatitis) = number of cases*risk of death
      num_death_latrine<-num_fulminant_latrine*risk_death_fulminant
      
      
      #####################
      #Economic Evaluation 
      #####################
      
      C_latrine <- cost_parameters[1]            # Cost of latrine per unit ($US)
      Num_household <- cost_parameters[2]        # Number of household
      C_hospital <- cost_parameters[3]           # Cost of hospitalisation in case of fulminant per day ($US)
      hosp_day <- cost_parameters[4]             # number of hospitalisation day per fulminant case (days)
      l_cov<-cost_parameters[5]                  # Latrine coverage (%)
      
      Total_C_base <- C_hospital*hosp_day*num_fulminant_base
      Total_C_latrine <- C_latrine*Num_household*l_cov+C_hospital*hosp_day*num_fulminant_latrine
      
      #Check total cost of each strategy
      Total_C_base
      Total_C_latrine
      
      #Outcomes; DALY=LYDs+LYLs or life year with disablitiy (due to mild illness and fulminant) + life year lost (due to early death)
      #LYDs
      dw_mildill <- 0.1               #DALY weight due to mild illness from Hep E  (assume as 0.1)
      num_mildill_day <- 7            #Number of mild illness day from Hep E (assume as 1 week)
      
      dw_fulminant <- 0.3               #DALY weight due to fulminant  (assume as 0.3)
      num_fulminant_day <- 14           #Number of hospitalisation day for Hep E with fulminant case (assume as 2 weeks)
      
      mildill_day_base <- num_mildill_base*num_mildill_day             #number of mild illness-day
      mildill_day_latrine <- num_mildill_latrine*num_mildill_day       #number of mild illness-day with latrine
      
      fulminant_day_base <- num_fulminant_base*num_fulminant_day             #number of fulminant-day
      fulminant_day_latrine <- num_fulminant_latrine*num_fulminant_day       #number of fulminant-day with latrine
      
      LYDis_base <- (dw_mildill*mildill_day_base+dw_fulminant*fulminant_day_base)/365     #number of lifeyear with disability (Basecase)
      LYDis_latrine <- (dw_mildill*mildill_day_latrine+dw_fulminant*fulminant_day_latrine)/365  #number of lifeyear with disability (with latrine)
      
      #LYLs
      LYlost<-35        #Average life year lost per one death
      
      LYlost_base <- num_death_base*LYlost            #number of life year lost due to fulminant hepatitis death
      LYlost_latrine <- num_death_latrine*LYlost      
      
      
      #Total DALY (combined LYDs and LYLs) of each strategy
      Total_DALY_base <- LYDis_base+LYlost_base
      Total_DALY_latrine <- LYDis_latrine+LYlost_latrine
      
      #Check total DALY lost of each strategy
      Total_DALY_base
      Total_DALY_latrine
      
      ############################
      #Cost-Utility analysis 
      ############################
      # No intervention vs Latrine
      Net_cost <- Total_C_latrine-Total_C_base              # Incremental cost due to the intervention
      
      DALY_averted <- Total_DALY_base-Total_DALY_latrine    # Incremental benefit provided by the intervention (Number of DALY averted)
      ICER <-Net_cost/DALY_averted
      
      Res <- cbind(Net_cost,DALY_averted, ICER)
      Res
      
      # For the sensitivity analysis section delete the below "#"
    }
    CEA_fun(base_parameters, lat_parameters, cost_parameters)
  })
}

#Plot to see results on the ICER plane
#plot(DALY_averted, Net_cost, xlim=c(-50, 100), ylim=c(-80000,200000), ylab="Incremental costs", xlab="DALY averted")
#abline(v = 0, h = 0)
#abline(a=0, b=675,lty = 3)


shinyApp(ui = ui, server = server)