#
# This is the server logic of a Shiny web application. You can run the
# application by clicking 'Run App' above.
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

# Define server logic required to draw a histogram
shinyServer(function(input, output) {

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
  
  
  ## NON-HIERCARICAL
  

  result.k <- reactive({
    kmeans(data_5, centers=input$k_n)
  })
  
  group.k <- reactive({
    as.factor(result.k()$cluster)
  })
  
  cluster_group_n <- reactive({
    if (input$meth == 'kmeans'){
      return(group.k())
    }
    if (input$meth == 'AUC'){
      group <- ifelse(rowMeans(data_5)>0.90,3,2)
      group[which(rowMeans(data_5)<0.82)]=1
      group = as.factor(group)
      return(group)
    }
  }) 
  
  
  lda.fit_n <- reactive({lda(cluster_group_n() ~., data= data[,1:13])}) 
  Lda.pred_n <- reactive({predict(lda.fit_n(), as.data.frame(data[,1:13]))}) 
  
  n_n       <- reactive({ length(cluster_group_n()) })      # total number of observations
  ng_n      <- reactive({ table(cluster_group_n())  })      # number of obs. in each group
  group_names_n   <- reactive({ levels(cluster_group_n())  })     # name of groups
  g_n       <- reactive({ length(group_names_n)})
  
  
  misc_n <- reactive({table(class.true=cluster_group_n(), class.assigned=Lda.pred_n()$class)})
  #print(misc_n) #CONFUSION MATRIX
  errors <- reactive({(Lda.pred_n()$class != cluster_group_n())})
  APER_n <- reactive({
    APER_n=0
    for(gi in 1:g_n()){
      APER_n <- APER_n + sum(misc_n()[gi,-gi])/sum(misc_n()[gi,]) * lda.fit_n()$prior[gi]}
    APER_n
  })
  
  LdaCV.aut_n <- reactive({lda(cluster_group_n() ~., data= data[,1:13], CV=TRUE) })  
  misc_cv_n <- reactive({table(class.true=cluster_group_n(), class.assigned=LdaCV.aut_n()$class)})
  errorsCV <- reactive({(LdaCV.aut_n()$class != cluster_group_n())})
  AERCV_n <- reactive({
    AERCV_n=0
    for(gi in 1:g_n()){
      AERCV_n <- AERCV_n + sum(misc_cv_n()[gi,-gi])/sum(misc_cv_n()[gi,]) * lda.fit_n()$prior[gi]}
    AERCV_n
  })
  
  
  # df1 <- reactive({pivot_longer(as.data.frame(t(colMeans(data_pc[cluster_group_n()==1,]))), cols=1:13, names_to = "Gene", values_to = "Expression")
  # df1$Group <- rep('1',13)})
  # 
  # df2 <- reactive({pivot_longer(as.data.frame(t(colMeans(data_pc[cluster_group_n()==2,]))), cols=1:13, names_to = "Gene", values_to = "Expression")
  # df2$Group <- rep('2',13)})
  # 
  # 
  # df3 <- reactive({pivot_longer(as.data.frame(t(colMeans(data_pc[cluster_group_n()==3,]))), cols=1:13, names_to = "Gene", values_to = "Expression")
  # df3$Group <- rep('3',13)})
  
  # df4 <- pivot_longer(as.data.frame(t(colMeans(data_pc[group==4,]))), cols=1:13, names_to = "Gene", values_to = "Expression")
  # df4$Group <- rep('4',13)
  
  # total <- reactive({rbind(df1(),df2(),df3())})
  
  # ggplot(total, aes(fill=Group, y=Expression, x=Gene)) + 
  #   geom_bar(position="dodge", stat="identity")
  
  # output$plot_silhouette <- renderPlot({  
  #   fviz_nbclust(data_5, FUN = kmeans, method = "silhouette")
  # })
  
  output$plot_wss_n <- renderPlot({  
    if (input$meth == 'kmeans'){
      fviz_nbclust(data_5, FUN = kmeans, method = "wss")}
  })
  
  output$plot_mean_n <- renderPlotly({
    plot_ly(data = data_plot, x = ~v1, y = ~v2, z = ~v3,
            mode   = 'markers',
            color = ~cell_means,
            type="scatter3d"
            #colorscale='earth'
    ) %>% layout(title = 'Visualization of cells on AUC first 3 PCs',
                 legend = list(title=list(text='average of treatment efficacy'))
    )#colors based on treatment efficacy average  
  })
  
  output$plot_meth <- renderPlotly({
    name = input$meth
    # generate 3d plot based on the name of clusters
    plot_ly(data = data_plot, x = ~v1, y = ~v2, z = ~v3,
            mode   = 'markers',
            type="scatter3d",
            color = as.character(cluster_group_n())
            #colorscale='earth'
    ) %>% layout(title = 'Visualization of clusters on first 3 PCs'
    )
    #fig <- fig %>% add_markers()
  })
  
  output$plot_meth_2 <- renderPlotly({
    name = input$meth
    # generate 3d plot based on the name of clusters
    plot_ly(data = data_plot, x = ~v1, y = ~v2, z = ~v3,
            mode   = 'markers',
            type="scatter3d",
            color = as.character(cluster_group_n())
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
  
  output$plot_classified_n <- renderPlotly({
    name = input$meth
    # generate 3d plot based on the name of clusters
    plot_ly(data = data_plot, x = ~v1, y = ~v2, z = ~v3,
            mode   = 'markers',
            type="scatter3d",
            color = as.character(Lda.pred_n()$class)
            #colorscale='earth'
    ) %>% layout(title = 'Visualization of classified on first 3 PCs'
    ) })
  
  # output$box_plot <- renderPlot({
  #   ggplot(total(), aes(fill=Group, y=Expression, x=Gene)) + geom_bar(position="dodge", stat="identity")
  # })
  
  
  #fig <- fig %>% add_markers()
  
  
  output$table_n <- renderPrint({ misc_n()})
  output$table_CV_n <- renderPrint({ misc_cv_n()})
  
  output$APER_n <- renderPrint({ 
    "APER:"
    APER_n()})
  
  output$AERCV_n <- renderPrint({ 
    AERCV_n()})
  
  

})
