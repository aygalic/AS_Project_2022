#
# This is the user-interface definition of a Shiny web application. You can
# run the application by clicking 'Run App' above.
#
# Find out more about building applications with Shiny here:
#
#    http://shiny.rstudio.com/
#

library(shiny)
library(shinythemes)
library(DT)
library("zoo")
library(mvtnorm)
library(rgl)
library(car)
library(plotly)
library(shiny)
library(factoextra)
library(MASS)

### LOAD DATA
#setwd("/Users/lucamainini/Documents/GitHub/AS_Project_2022")
#load(file.path("Dataset","breast_auc_data.Rdata"))
load("breast_auc_data.Rdata")
load("selected_genes_data.Rdata")
breast_auc = t(breast_data_treatment_auc)
data_5 = na.aggregate(breast_auc)

# per mostrare dati
# vchoices <- 1:ncol(data_5)
# names(vchoices) <- names(data_5)


### MEANS
cell_means = apply(data_5,1,mean)

### PCA
cov2 = cov(data_5)
decomp2 = eigen(cov2)
P3 = as.matrix(decomp2$vectors[,1:3])
reduced_M_ = data_5%*%P3
colnames(reduced_M_) <- c("v1","v2", "v3")

data_plot = data.frame(reduced_M_)
data_plot$cell_means= cell_means

# Define UI for application that draws a histogram
shinyUI(fluidPage(theme = shinytheme("united"), #cerulean, united, flatly
                  titlePanel("Clusters on treatment efficacy (AUC score)"),
                  
                
                  
                  navbarPage("Type of clusters",
                             tabPanel(icon("home"),
                                      scroller::use_scroller(), # add use_scroller() in the UI
                                      #fluidRow(column(
                                      tags$head(tags$style( # per modificare dimensione foto in base dispositivo
                                        type="text/css",
                                        "#d img {max-width: 95%; max-height: 95%; width: 95%; height: auto}"
                                      )),
                                        imageOutput("d"),
                                        #tags$img(src="project_scope.png",width="480px",height="270px")
                                       # ,width=6)),
                                       # 
                                      #br(),
                                      #br(),
                                      #br(),  

                                      p("This interactive app was developed to visualise the clusters obtained on AUC values. 
                                        The data was taken from Cancer Cell Line Encyclopedia. 
                                        After being preprocessed, you can view the AUC dataset here", 
                                        #strong("But do not worry!"), "you will find alternatives to learn all these technical aspects independently.",
                                        style="text-align:justify;color:black;background-color:papayawhip;padding:15px;border-radius:10px"), #color:papayawhip,lavender
                                      #br(),
                                      selectInput("select", "Select 1 or more treatments to display", colnames(data_5), multiple = TRUE),
                                      
                                      fluidRow(column(DT::dataTableOutput("mytable"),
                                                      width = 12)),
                                      
                                      tags$head(tags$style( # per modificare dimensione foto in base dispositivo
                                        type="text/css",
                                        "#step1 img {max-width: 90%; width: 90%; height: auto}"
                                      )),
                                      
                      
                                      #fluidRow(column(
                                        imageOutput("step1"),
                                        #tags$img(src="STEP1.png",width="480px",height="270px"),
                                        #width=6)),
                                      #br(),
                                      #br(), 
                                        
                                      p("The aim of this step is to divide the patients according to the effectiveness of the treatments. 
                                      For example, we could try to obtain a group where all treatments are effective and one where treatments generally do not work. 
                                      Remember we define the drug performance as in-vitro efficacy (measured by AUC scores).
                                      In the other sections of this app, you can view the", strong("clusters"), "obtained with a", strong("hierarchical approach"), "or with a", strong("non-hierarchical one")),
                                      p("The first approach shows poor results, whereas k-means divides the clusters quite similarly to a cluster based on the average AUC efficacy.
                                      "),
                                      p("In the second step, key genes are found."),
                                      tags$head(tags$style( # per modificare dimensione foto in base dispositivo
                                        type="text/css",
                                        "#step2 img {max-width: 90%; width: 90%; height: auto}"
                                      )),
                                      #fluidRow(column(
                                        imageOutput("step2"),
                                        
                                        #width=6)),
                                      #br(),
                                      #br(), 
                                      
                                      p("At the end of each section, these genes are used to check which group patients belong to (LDA classification). 
                                      An estimate of the prediction error is provided."),
                                      a("Scroll to top", type = "button", class = "btn btn-danger", href = "#d"), #.btn
                                      br(),   
                                      
                             ),
                             
                             tabPanel("Non-Hierarchical",
                                      fluidRow(column(3,
                                                      selectInput('meth', 'Choose a distance', choices = list(
                                                        "K-mean" = 'kmeans', "Average AUC" = 'AUC'),
                                                        selected = 'kmeans')
                                      ),
                                      
                                      column(3,
                                             sliderInput("k_n", 'Choose number of clusters', min = 2,  max = 10, value = 3)
                                      ),
                                      
                                      
                                      ),
                                      
                                      
                                      fluidRow(
                                        # column(6,plotOutput("plot_silhouette")),
                                        column(6,plotOutput("plot_wss_n"))
                                        
                                      ),
                                      p("Clusters are obtained using all treatments. 
                                        They are projected onto the space of the first three PCs for visualisation purposes only."),
                                      
                                      fluidRow(
                                        column(5,plotlyOutput("plot_mean_n")),
                                        column(5,plotlyOutput("plot_meth"))
                                      ),
                                      h2("Results of LDA Classifation using influential gene "),
                                      fluidRow(
                                        
                                        column(5,plotlyOutput("plot_classified_n")),
                                        column(5,plotlyOutput("plot_meth_2"))
                                      ),
                                      
                                      fluidRow(
                                        
                                        column(5,
                                               p("Prediction matrix"),
                                               verbatimTextOutput("table_n")),
                                        column(3,
                                               p("APER"),
                                               verbatimTextOutput("APER_n")),
                                      ),
                                      h3("Cross Validation "),
                                      fluidRow(
                                        
                                        column(5,
                                               p("CV - Prediction matrix"),
                                               verbatimTextOutput("table_CV_n")),
                                        column(3,
                                               p("CV - APER"),
                                               verbatimTextOutput("AERCV_n"))
                                      ),
                                      br(),       
                             ), 
                             tabPanel("Hierarchical",
                                      fluidRow(column(3,
                                                      selectInput('dist', 'Choose a distance', choices = list(
                                                        "Eucledian" = 'euclidean', "Manhattan" = 'manhattan', "Canberra" = 'canberra'),
                                                        selected = 'euclidean')
                                      ),
                                      column(3,
                                             selectInput('linkage', 'Choose a Linkage', choices = list(
                                               "Single" = 'single', "Average" = 'average', "Complete" = 'complete', "Ward"='ward.D2'),
                                               selected = 'average')
                                      ),
                                      
                                      column(3,
                                             sliderInput("k", 'Choose number of clusters', min = 2,  max = 10, value = 3)
                                      ),
                                      
                                      column(3, 
                                             p("Cophenetic Correlation Coefficient:"),
                                             verbatimTextOutput("cophern_value"))
                                      
                                      
                                      ),
                                      
                                      
                                      fluidRow(
                                        column(5,plotOutput("plot_dendogram")),
                                        column(5,plotlyOutput("plot_dist"))
                                        
                                      ),
                                      
                                      p("Clusters are obtained using all treatments. 
                                        They are projected onto the space of the first three PCs for visualisation purposes only."),
                                      fluidRow(
                                        column(5,plotlyOutput("plot_mean")),
                                        column(5,plotlyOutput("plot_dist_2"))
                                      ),
                                      
                                      h2("Results of LDA Classifation using influential gene "),
                                      fluidRow(
                                        
                                        column(5,plotlyOutput("plot_classified")),
                                        column(5,plotlyOutput("plot_dist_3"))
                                      ),
                                      
                                      fluidRow(
                                        
                                        column(5,
                                               p("Prediction matrix"),
                                               verbatimTextOutput("table")),
                                        column(3,
                                               p("APER"),
                                               verbatimTextOutput("APER")),
                                      ),
                                      h3("Cross Validation "),
                                    
                                      fluidRow(
                                        
                                        column(5,
                                               p("CV - Prediction matrix"),
                                               verbatimTextOutput("table_CV")),
                                        column(3,
                                               p("CV - APER"),
                                               verbatimTextOutput("AERCV"))
                                        
                                      ),
                                      br(),         
                             ),
                  )
))