library(shiny)
library(boot)
library(MASS)
library(segmented)
library(sm)
library(mixtools)


shinyUI(fluidPage(
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
    )))