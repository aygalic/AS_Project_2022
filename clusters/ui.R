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
                                      fluidRow(column(tags$img(src="project_scope.png",width="480px",height="270px"),width=2)),
                                      fluidRow(column(tags$img(src="STEP1.png",width="480px",height="270px"),width=2))
                                 
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
                                      
                                      column(3, verbatimTextOutput("cophern_value"))
                                      
                                      
                                      ),
                                      
                                      
                                      fluidRow(
                                        column(5,plotOutput("plot_dendogram")),
                                        column(5,plotlyOutput("plot_dist"))
                                        
                                      ),
                                      
                                      
                                      fluidRow(
                                        column(5,plotlyOutput("plot_mean")),
                                        column(5,plotlyOutput("plot_dist_2"))
                                      ),
                                      
                                      
                                      fluidRow(
                                        
                                        column(5,plotlyOutput("plot_classified")),
                                        column(5,plotlyOutput("plot_dist_3"))
                                      ),
                                      
                                      fluidRow(
                                        
                                        column(5,verbatimTextOutput("table")),
                                        column(3,verbatimTextOutput("APER")),
                                      ),
                                      
                                      fluidRow(
                                        
                                        column(5,verbatimTextOutput("table_CV")),
                                        column(3,verbatimTextOutput("AERCV"))
                                      ),
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
                                      
                                      
                                      fluidRow(
                                        column(5,plotlyOutput("plot_mean_n")),
                                        column(5,plotlyOutput("plot_meth"))
                                      ),
                                      
                                      fluidRow(
                                        
                                        column(5,plotlyOutput("plot_classified_n")),
                                        column(5,plotlyOutput("plot_meth_2"))
                                      ),
                                      
                                      fluidRow(
                                        
                                        column(5,verbatimTextOutput("table_n")),
                                        column(3,verbatimTextOutput("APER_n")),
                                      ),
                                      
                                      fluidRow(
                                        
                                        column(5,verbatimTextOutput("table_CV_n")),
                                        column(3,verbatimTextOutput("AERCV_n"))
                                      ),
                                      br(),       
                             )
                  )
))