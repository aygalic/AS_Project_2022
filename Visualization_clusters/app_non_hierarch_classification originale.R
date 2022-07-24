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

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Can we found any good custer?"),

    # Sidebar with a slider input for number of bins 
    fluidRow(column(3,
              selectInput('dist', 'Choose a distance', choices = list(
                "K-mean" = 'kmeans', "K-mean" = 'kmeans'),
                selected = 'kmeans')
    ),
    
            column(3,
                   sliderInput("k", 'Choose number of clusters', min = 2,  max = 10, value = 3)
            ),
    
    
    ),
  
    
    fluidRow(
    # column(6,plotOutput("plot_silhouette")),
      column(6,plotOutput("plot_wss"))
    
    ),
    
    
    fluidRow(
      column(5,plotlyOutput("plot_mean")),
      column(5,plotlyOutput("plot_dist"))
    ),
    
    fluidRow(
      
      column(5,plotlyOutput("plot_classified")),
      column(5,plotlyOutput("plot_dist_2"))
    ),
    
    fluidRow(
      
      column(5,verbatimTextOutput("table")),
      column(3,verbatimTextOutput("APER")),
    ),
    
    fluidRow(
      
      column(5,verbatimTextOutput("table_CV")),
      column(3,verbatimTextOutput("AERCV"))
    ),
    # 
    # fluidRow(
    #   
    #   column(6,plotOutput("box_plot")),
    # )
)


# Define server logic required to draw a histogram
server <- function(input, output) {
  
  cluster_group <- reactive({
     result.k <- kmeans(data_5, centers=input$k)
   }) 
  
  # data$y <- reactive({
  #   cluster_group()$cluster
  # }) 
  
  lda.fit <- reactive({lda(cluster_group()$cluster ~., data= data[,1:13])}) 
  Lda.pred <- reactive({predict(lda.fit(), as.data.frame(data[,1:13]))}) 
  
  n       <- reactive({ length(cluster_group()$cluster) })      # total number of observations
  ng      <- reactive({ table(cluster_group()$cluster)  })      # number of obs. in each group
  group_names   <- reactive({ levels(cluster_group()$cluster)  })     # name of groups
  g       <- reactive({ length(group_names)})
  
  
  misc <- reactive({table(class.true=cluster_group()$cluster, class.assigned=Lda.pred()$class)})
  #print(misc) #CONFUSION MATRIX
  errors <- reactive({(Lda.pred()$class != cluster_group()$cluster)})
  APER <- reactive({
    APER=0
  for(gi in 1:g()){
    APER <- APER + sum(misc()[gi,-gi])/sum(misc()[gi,]) * lda.fit()$prior[gi]}
    APER
    })
  
  LdaCV.aut <- reactive({lda(cluster_group()$cluster ~., data= data[,1:13], CV=TRUE) })  
  misc_cv <- reactive({table(class.true=cluster_group()$cluster, class.assigned=LdaCV.aut()$class)})
  errorsCV <- reactive({(LdaCV.aut()$class != cluster_group()$cluster)})
  AERCV <- reactive({
    AERCV=0
    for(gi in 1:g()){
      AERCV <- AERCV + sum(misc_cv()[gi,-gi])/sum(misc_cv()[gi,]) * lda.fit()$prior[gi]}
    AERCV
  })
  
  
  # df1 <- reactive({pivot_longer(as.data.frame(t(colMeans(data_pc[cluster_group()$cluster==1,]))), cols=1:13, names_to = "Gene", values_to = "Expression")
  # df1$Group <- rep('1',13)})
  # 
  # df2 <- reactive({pivot_longer(as.data.frame(t(colMeans(data_pc[cluster_group()$cluster==2,]))), cols=1:13, names_to = "Gene", values_to = "Expression")
  # df2$Group <- rep('2',13)})
  # 
  # 
  # df3 <- reactive({pivot_longer(as.data.frame(t(colMeans(data_pc[cluster_group()$cluster==3,]))), cols=1:13, names_to = "Gene", values_to = "Expression")
  # df3$Group <- rep('3',13)})
  
  # df4 <- pivot_longer(as.data.frame(t(colMeans(data_pc[group==4,]))), cols=1:13, names_to = "Gene", values_to = "Expression")
  # df4$Group <- rep('4',13)
  
  total <- reactive({rbind(df1(),df2(),df3())})
  
  # ggplot(total, aes(fill=Group, y=Expression, x=Gene)) + 
  #   geom_bar(position="dodge", stat="identity")

  # output$plot_silhouette <- renderPlot({  
  #   fviz_nbclust(data_5, FUN = kmeans, method = "silhouette")
  # })

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

  output$plot_dist_2 <- renderPlotly({
    name = input$dist
    # generate 3d plot based on the name of clusters
    plot_ly(data = data_plot, x = ~v1, y = ~v2, z = ~v3,
            mode   = 'markers',
            type="scatter3d",
            color = as.character(cluster_group()$cluster)
            #colorscale='earth'
    ) %>% layout(title = 'Visualization of clusters on first 3 PCs'
    )#  %>% add_trace(
    #   data = data_plot[names,]
    #   , x = ~v1
    #   , y = ~v2
    #   , z = ~v3
    #   , color= "k"
    #   , mode = "markers"
    #   , type = "scatter3d"
    #   , marker = list(size = 5))
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
    ) })
    
    # output$box_plot <- renderPlot({
    #   ggplot(total(), aes(fill=Group, y=Expression, x=Gene)) + geom_bar(position="dodge", stat="identity")
    # })
    
    
    #fig <- fig %>% add_markers()

  
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
