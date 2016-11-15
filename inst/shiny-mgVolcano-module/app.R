library(shiny)
library(DT)
library(dplyr)
library(shinyBS)

library(rprojroot)
devtools::load_all(find_root(is_r_package))

server <- function(input, output, session) {

  mgc <- reactive({
    mg <- readRDS('~/tmp/schmidt/multiGSEA-EP-uber_hWT_tKO-hWT_tWT.rds')
    MultiGSEAResultContainer(mg)
  })

  volcano <- callModule(mgVolcano, 'volcano', mgc)

  output$brushed <- DT::renderDataTable({
    req(volcano()) %>%
      select(symbol, featureId, logFC, pval, padj) %>%
      datatable
  })

  observeEvent(volcano(), {
    brushed <- volcano()
    if (is.null(brushed)) {
      msg("No brush here")
    } else {
      msg("Brushed", nrow(brushed), "genes")
    }
  })
}

ui <- fluidPage(
  fluidRow(
    column(4, mgVolcanoUI("volcano", hexbin=TRUE)),
    column(8, DT::dataTableOutput('brushed'))
  ))

shinyApp(ui=ui, server=server)
