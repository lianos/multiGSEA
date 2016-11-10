library(shiny)

library(rprojroot)
devtools::load_all(find_root(is_r_package))

server <- function(input, output, session) {

  mgc <- reactive({
    mg <- readRDS('~/tmp/schmidt/multiGSEA-EP-uber_hWT_tKO-hWT_tWT.rds')
    MultiGSEAResultContainer(mg)
  })

  callModule(mgVolcano, 'volcano', mgc)
}

ui <- fluidPage(
  fluidRow(
    mgVolcanoUI("volcano", hexbin=TRUE)
  ))

shinyApp(ui=ui, server=server)
