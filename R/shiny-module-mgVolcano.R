##' Creates a volcano plot from a source input
##'
##' @export
##' @rdname mgVolcano
##' @importFrom shiny NS tagList tags sliderInput
##' @importFrom shinyjs useShinyjs hidden
##' @param id the shiny id of the output widget
##' @param x the object to build a volcano plot from
##' @param stats the stats to pull from \code{x} (if necessary) to build the
##'   volcano for.
##' @param x the name of the column from the stats table to use on x axis
##' @param y the name of the column from the stats table to use on y axis
mgVolcanoUI <- function(id, x, stats='dge', xaxis='logFC', yaxis='padj',
                        idx='idx', hexbin=TRUE) {
  ns <- NS(id)
  if (hexbin) {
    out <- tagList(
      useShinyjs(),
      tags$a(id=ns('settings'), icon("wrench")),
      # rbokehOutput(ns("plot")),
      plotlyOutput(ns("plot")),
      hidden(
        tags$div(
          id=ns('widgets'),
          sliderInput(
            ns("xhex"), 'x filter', min=0, max=5, step=0.25, value=1),
          sliderInput(
            ns("yhex"), 'y filter', min=0, max=1, step=0.025, value=0.10))))
  } else {
    out <- rbokehOutput(ns("plot"))
  }
  out
}

##' @section Module Return:
##' This module returns a data.frame with the information of the selected
##' features. If no genes are selected in a volcano, it will return \code{NULL}.
##'
##' @rdname mgVolcano
##' @export
##' @importFrom shiny reactive observeEvent updateSliderInput
##' @importFrom shinyjs onclick toggle
##' @rdname mgVolcano
mgVolcano <- function(input, output, session,
                      x, stats='dge', xaxis='logFC', yaxis='pval', idx='idx',
                      tools=c('box_select', 'reset', 'save')) {
  onclick("settings", toggle(id="widgets", anim=TRUE))
  if (missing(idx)) {
    if (stats == 'dge') idx <- 'featureId'
  }

  ## Extract the data used in the volcano to keep it handy
  dat <- reactive({
    req(x())
    volcano.stats.table(x(), stats, xaxis, yaxis, idx)
  })

  ## If UI is showing the hexbin sliders, fix ranges and labels when dat()
  ## is initialized
  observeEvent(dat(), {
    if (!is.null(input$xhex)) {
      updateSliderInput(session, 'yhex', sprintf('%s filter', yaxis),
                        min=0, max=1, step=0.025, value=0.10)
      # max.x <- ceiling(max(abs(dat()[['xaxis']]))) - 0.5
      max.x <- ceiling(max(abs(dat()[['.xvt']]))) - 0.5
      updateSliderInput(session, 'xhex', sprintf('%s filter', xaxis),
                        min=0, max=max.x, step=0.25, value=1)
    }
  })

  plt <- reactive({
    req(x())
    ns <- session$ns
    xhex <- input$xhex
    yhex <- input$yhex
    p <- volcano_plot(x(), stats, xaxis, yaxis, idx, xhex=xhex, yhex=yhex,
                      tools=tools, shiny_source='mgvolcano',
                      width=400, heigh=350)
    p
  })

  output$plot <- renderPlotly({
    req(plt())
  })

  ## This module returns a data.frame containing info genes that are brushed
  ## by the user
  vals <- reactive({
    dat <- req(plt()) %>% plotly_data
    event <- event_data('plotly_selected', source='mgvolcano')
    if (!is.null(event)) {
      # dat <- isolate(plt()) %>% plotly_data
      # selected <- subset(dat, featureId %in% event$key)
      out <- subset(dat, featureId %in% event$key)
    } else {
      out <- NULL
    }
    out
  })

  return(vals)
}
