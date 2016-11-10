server <- function(input, output, session) {

  mgc <- reactive({
    mg <- readRDS('~/tmp/schmidt/multiGSEA-EP-uber_hKO_tWT-hWT_tWT.rds')
    MultiGSEAResultContainer(mg)
  })

  callModule(mgVolcano, 'volcano', mgc, )
  output$distPlot <- renderPlot({
    hist(rnorm(input$obs), col = 'darkgray', border = 'white')
  })
}

ui <- fluidPage(
  fluidRow(
    mgVolcanoUI("volcano")
  ))

shinyApp(ui=ui, server=server)
