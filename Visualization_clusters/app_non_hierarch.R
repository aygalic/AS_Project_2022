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

### LOAD DATA
#setwd("/Users/lucamainini/Documents/GitHub/AS_Project_2022")
#load(file.path("Dataset","breast_auc_data.Rdata"))
load("breast_auc_data.Rdata")
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

### COMPUTATION OF DISTANCES 

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Can we found any good custer?"),

    # Sidebar with a slider input for number of bins 
    fluidRow(column(3,
              selectInput('dist', 'Choose a distance', choices = list(
                "K-mean" = 'kmeans', "K-mean" = 'kmeans'),
                selected = 'euclidean')
    ),
    
            column(3,
                   sliderInput("k", 'Choose number of clusters', min = 2,  max = 10, value = 3)
            ),
    
    
    ),
  
    
    fluidRow(
    column(6,plotOutput("plot_silhouette")),
    #column(6,plotlyOutput("plot_wss"))
    
    ),
    
    
    fluidRow(
      column(5,plotlyOutput("plot_mean")),
      column(5,plotlyOutput("plot_dist"))
    )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  cluster_group <- reactive({
     result.k <- kmeans(data_5, centers=input$k)
   }) 

  output$plot_silhouette <- renderPlot({  
    fviz_nbclust(data_5, FUN = kmeans, method = "silhouette")
  })

  output$plot_wss <- renderPlot({  
    fviz_nbclust(data_5, FUN = kmeans, method = "wss")
  })

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
            color = as.character(cluster_group()$cluster)
            #colorscale='earth'
    ) %>% layout(title = 'Visualization of clusters on first 3 PCs'
    )
    #fig <- fig %>% add_markers()
  })
  
}

# Run the application 
shinyApp(ui = ui, server = server)
