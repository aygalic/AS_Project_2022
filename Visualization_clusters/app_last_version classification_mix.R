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

### COMPUTATION OF DISTANCES 
# d.e <- dist(data_5, method='euclidean')
# d.m <- dist(data_5, method='manhattan')
# d.c <- dist(data_5, method='canberra')
# 
# d.es <- hclust(d.e, method='single')
# d.ea <- hclust(d.e, method='average')
# d.ec <- hclust(d.e, method='complete')
# d.w <- hclust(d.e, method='ward.D2')
# 
# d.ms <- hclust(d.m, method='single')
# d.ma <- hclust(d.m, method='average')
# d.mc <- hclust(d.m, method='complete')
# d.mw <- hclust(d.m, method='ward.D2')
# 
# d.cs <- hclust(d.c, method='single')
# d.ca <- hclust(d.c, method='average')
# d.cc <- hclust(d.c, method='complete')
# d.cw <- hclust(d.c, method='ward.D2')


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
)
)


# Define server logic required to draw a histogram
server <- function(input, output) {
  dist_mat <- reactive({
    dist(data_5, method=input$dist)
  })
  data_c <- reactive({
    hclust(dist_mat(), method=input$linkage)
  })
  
  cluster_group <- reactive({
    return(as.factor(cutree(data_c(), k=input$k)))
  }) 
  
  title <- reactive({
    paste("Dendogram", input$dist,input$linkage, sep = " - ")
  })
  
  cophern <- reactive({
    cor(dist_mat(), cophenetic(data_c()))
    })
  
  output$cophern_value <- renderPrint({ cophern()})

  output$plot_dendogram <- renderPlot({  
    plot(data_c(), main=title(), hang=-0.1, xlab='', labels=F, cex=0.6, sub='')
    rect.hclust(data_c(), k=input$k, border = 2:5)
  })
  
  lda.fit <- reactive({lda(cluster_group() ~., data= data[,1:13])}) 
  Lda.pred <- reactive({predict(lda.fit(), as.data.frame(data[,1:13]))}) 
  
  lda.fit <- reactive({lda(cluster_group() ~., data= data[,1:13])}) 
  Lda.pred <- reactive({predict(lda.fit(), as.data.frame(data[,1:13]))}) 
  
  n       <- reactive({ length(cluster_group()) })      # total number of observations
  ng      <- reactive({ table(cluster_group())  })      # number of obs. in each group
  group_names   <- reactive({ levels(cluster_group())  })     # name of groups
  g       <- reactive({ length(group_names)})
  
  
  misc <- reactive({table(class.true=cluster_group(), class.assigned=Lda.pred()$class)})
  #print(misc) #CONFUSION MATRIX
  errors <- reactive({(Lda.pred()$class != cluster_group())})
  APER <- reactive({
    APER=0
    for(gi in 1:g()){
      APER <- APER + sum(misc()[gi,-gi])/sum(misc()[gi,]) * lda.fit()$prior[gi]}
    APER
  })
  
  LdaCV.aut <- reactive({lda(cluster_group() ~., data= data[,1:13], CV=TRUE) })  
  misc_cv <- reactive({table(class.true=cluster_group(), class.assigned=LdaCV.aut()$class)})
  errorsCV <- reactive({(LdaCV.aut()$class != cluster_group())})
  AERCV <- reactive({
    AERCV=0
    for(gi in 1:g()){
      AERCV <- AERCV + sum(misc_cv()[gi,-gi])/sum(misc_cv()[gi,]) * lda.fit()$prior[gi]}
    AERCV
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
  
  output$plot_dist_3 <- renderPlotly({
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
  
  output$plot_classified <- renderPlotly({
    name = input$dist
    # generate 3d plot based on the name of clusters
    plot_ly(data = data_plot, x = ~v1, y = ~v2, z = ~v3,
            mode   = 'markers',
            type="scatter3d",
            color = as.character(Lda.pred()$class)
            #colorscale='earth'
    ) %>% layout(title = 'Visualization of classified on first 3 PCs'
    )
    #fig <- fig %>% add_markers()
  })
  
  output$table <- renderPrint({ misc()})
  output$table_CV <- renderPrint({ misc_cv()})
  
  output$APER <- renderPrint({ 
    "APER:"
    APER()})
  
  output$AERCV <- renderPrint({ 
    AERCV()})
  
}

# Run the application 
shinyApp(ui = ui, server = server)
