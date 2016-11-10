##' Creates a volcano plot from a source input
##' @export
##' @param id the shiny id of the output widget
##' @param x the object to build a volcano plot from
##' @param stats the stats to pull from \code{x} (if necessary) to build the
##'   volcano for.
##' @param x the name of the column from the stats table to use on x axis
##' @param y the name of the column from the stats table to use on y axis
mgVolcanoUI <- function(id, x, stats='dge', xaxis='logFC', yaxis='padj',
                        idx='idx', hexbin=TRUE) {
  stopifnot(requireNamespace('shiny'))
  ns <- shiny::NS(id)
  if (hexbin) {
    out <- shiny::tagList(
      sliderInput(ns("xhex"), 'x filter',
                  min=0, max=5, step=0.5, value=1),
      sliderInput(ns("yhex"), 'y filter',
                  min=0, max=1, step=0.025, value=0.10),
      rbokehOutput(ns("volcano_plot")))
  } else {
    out <- rbokehOutput(ns("volcano_plot"))
  }

  out
}

##' @export
mgVolcano <- function(input, output, session,
                      x, stats='dge', xaxis='logFC', yaxis='pval', idx='idx') {
  stopifnot(requireNamespace('shiny'))
  if (FALSE) {
    x <- mg <- readRDS('~/tmp/schmidt/multiGSEA-EP-uber_hKO_tWT-hWT_tWT.rds')
    stats='dge'; xaxis='logFC'; yaxis='padj'; idx='idx';
    xgran=NULL; ygran=NULL
  }

  ## Reset the labels and ranges in the "hexbin" sliders, if the UI element
  ## has hexbin enabled.
  observeEvent(x(), {
    req(input$xhex)
    dat <- volcano.stats.table(x(), stats, xaxis, yaxis, idx)
    updateSliderInput(session, 'yhex', sprintf('%s filter', yaxis),
                      min=0, max=1, step=0.025, value=0.10)
    max.x <- ceiling(max(abs(dat[['xaxis']]))) - 0.5
    updateSliderInput(session, 'xhex', sprintf('%s filter', xaxis),
                      min=0, max=max.x, step=0.5, value=1)
  })

  output$volcano_plot <- renderRbokeh({
    xhex <- input$xhex
    yhex <- input$yhex
    if (!is.null(xhex)) {
      # dat <- volcano.stats.table(x(), stats, xaxis, yaxis, idx)
      # yhex <- approx.target.from.transformed(input$yhex, dat$pval, dat$padj)
    } else {
      yhex <- NULL
    }
    volcano_plot(x(), stats, xaxis, yaxis, idx, xhex=xhex, yhex=yhex)
  })
}
