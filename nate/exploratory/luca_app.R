#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.

library("zoo")
library(mvtnorm)
library(rgl)
library(car)
library(plotly)
library(shiny)
library(factoextra)

source("../utils/nate_utils.R")
source("../utils/exp_base_script.R")

### LOAD DATA
#setwd("/Users/lucamainini/Documents/GitHub/AS_Project_2022")
#load(file.path("Dataset","breast_auc_data.Rdata"))
load("breast_cns.Rda")
data_5 = as.matrix(tosave)

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
ui <- fluidPage(
  
  # Application title
  titlePanel("Can we found any good custer?"),
  
  # Sidebar with a slider input for number of bins 
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
  
  
  fluidRow(column(3,
                  selectInput('cancer', 'Choose a cancer type', choices = cancer_types,
                    selected = 'BREAST', multiple=TRUE)
  )
  
  
  ),
  
  
  fluidRow(
    column(5,plotOutput("plot_dendogram")),
    column(5,plotlyOutput("plot_dist"))
    
  ),
  
  
  fluidRow(
    column(5,plotlyOutput("plot_mean")),
    column(5,plotlyOutput("plot_dist_2"))
  )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  data_5 = reactive({
    block_dat(input$cancer, auc)
  })
  
  dist_mat <- reactive({
    dist(data_5, method=input$dist)
  })
  data <- reactive({
    hclust(dist_mat(), method=input$linkage)
  })
  
  cluster_group <- reactive({
    cutree(data(), k=input$k)
  }) 
  
  title <- reactive({
    paste("Dendogram", input$dist,input$linkage, sep = " - ")
  })
  
  cophern <- reactive({
    cor(dist_mat(), cophenetic(data()))
  })
  
  output$cophern_value <- renderPrint({ cophern()})
  
  output$plot_dendogram <- renderPlot({  
    plot(data(), main=title(), hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
    rect.hclust(data(), k=input$k, border = 2:5)
  })
  
  #fviz_nbclust(data_5, FUN = kmeans, method = "silhouette") 
  #fviz_nbclust(data_5, FUN = kmeans, method = "wss")
  
  output$plot_mean <- renderPlotly({
    plot_ly(data = data_plot, x = ~v1, y = ~v2, z = ~v3,
            mode   = 'markers',
            color = ~cell_means,
            type="scatter3d"
            #colorscale='earth'
    ) %>% layout(title = 'Visualization of cells on AUC first 3 PCs',
                 legend = list(title=list(text='average of treatment efficacy'))
    )#colors based on treatment efficacy average  
  })
  
  output$plot_dist <- renderPlotly({
    name = input$dist
    # generate 3d plot based on the name of clusters
    plot_ly(data = data_plot, x = ~v1, y = ~v2, z = ~v3,
            mode   = 'markers',
            type="scatter3d",
            color = as.character(cluster_group())
            #colorscale='earth'
    ) %>% layout(title = 'Visualization of clusters on first 3 PCs'
    )
    #fig <- fig %>% add_markers()
  })
  
  output$plot_dist_2 <- renderPlotly({
    name = input$dist
    # generate 3d plot based on the name of clusters
    plot_ly(data = data_plot, x = ~v1, y = ~v2, z = ~v3,
            mode   = 'markers',
            type="scatter3d",
            color = as.character(cluster_group())
            #colorscale='earth'
    ) %>% layout(title = 'Visualization of clusters on first 3 PCs'
    )
    #fig <- fig %>% add_markers()
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
