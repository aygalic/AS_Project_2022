#
# This is a Shiny web application. You can run the application by clicking
# the 'Run App' button above.

library(shiny)

# Define UI for application that draws a histogram
ui <- fluidPage(

    # Application title
    titlePanel("Can we found any good custer?"),

    # Sidebar with a slider input for number of bins 
    fluidRow(column(3,
      selectizeInput('dist', 'Choose a typer of clustering', choices = list(
        Eucledian = c(`Single Link` = 'cluster.es', `Average Link` = 'cluster.ea', `Complete Link` = 'cluster.ec', `Ward Link` = 'cluster.ew'),
        Camperra = c(`Single Link` = 'cluster.cs', `Average Link` = 'cluster.ca', `Complete Link` = 'cluster.cc', `Ward Link` = 'cluster.cw'),
        Manhattan = c(`Single Link` = 'cluster.ms', `Average Link` = 'cluster.ma', `Complete Link` = 'cluster.mc', `Ward Link` = 'cluster.mw'),
        K_means = c('K-means with k=3' = 'cluster.km_3')
      )
      )
    )),
    
    fluidRow(
      column(6,plotlyOutput("plot_mean")),
      column(6,plotlyOutput("plot2"))
)
)

# Define server logic required to draw a histogram
server <- function(input, output) {

    output$plot_mean <- renderPlotly({
      plot_ly(data = data_plot, x = ~v1, y = ~v2, z = ~v3,
              mode   = 'markers',
              color = ~cell_means,
              type="scatter3d"
              #colorscale='earth'
              ) %>% layout(title = 'Visualization of mean of treatment efficacy by cell on first 3 PCs')
    })
    
    output$plot2 <- renderPlotly({
      name = input$dist
      # generate 3d plot based on the name of clusters
      plot_ly(data = data_plot, x = ~v1, y = ~v2, z = ~v3,
              mode   = 'markers',
              type="scatter3d",
              color = as.character(data_plot[,name])
              #colorscale='earth'
              ) %>% layout(title = 'Visualization of clusters on first 3 PCs'
              )
      #fig <- fig %>% add_markers()
    })
    
    
}

# Run the application 
shinyApp(ui = ui, server = server)
