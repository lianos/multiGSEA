##' Creates a volcano plot from a source input
##' @export
##' @param id the shiny id of the output widget
##' @param x the object to build a volcano plot from
##' @param stats the stats to pull from \code{x} (if necessary) to build the
##'   volcano for.
##' @param x the name of the column from the stats table to use on x axis
##' @param y the name of the column from the stats table to use on y axis
mgVolcanoUI <- function(id, x, stats='dge', xaxis='logFC', yaxis='padj',
                        idx='idx', hex_granularity=TRUE) {
  requireNamespace('shiny')
  ns <- shiny::NS(id)
  dat <- volcano.stats.table(x, stats, xaxis, yaxis, idx)

  xvals <- dat[['xaxis']]
  yvals <- dat[['yaxis']]

  if (hex_granularity) {
    xrange <- range(xvals)
    yrane <- range(yvals)
    out <- shiny::tagList(
      sliderInput(ns("xgran"), sprintf('%s cloud', xaxis),
                  min=xrange[1], max=xrange[2], value=quantile()),
      sliderInput(ns("ygran"), sprintf('%s cloud', yaxis),
                  min=yrange[1], max=yrange[2], value=quantile()),
      rbokehOutput(ns("volcano_plot")))
  } else {
    out <- rbokehOutput(ns("volcano_plot"))
  }

  out
}

mgVolcano <- function(input, output, session,
                      x, stats='dge', xaxis='logFC', yaxis='padj', idx='idx',
                      hex_granularity=TRUE) {
  requireNamespace('shiny')

  if (FALSE) {
    x <- mg <- readRDS('~/tmp/schmidt/multiGSEA-EP-uber_hKO_tWT-hWT_tWT.rds')
    stats='dge'; xaxis='logFC'; yaxis='padj'; idx='idx'; hex_granularity=TRUE
  }

  dat <- volcano.stats.table(x, stats, xaxis, yaxis, idx)

  output$volcano_plot <- rbokehOutput({
    hex_volcano_plot(dat, xaxis, yaxis, idx,
                     xgran=input$xgran, ygran=input$ygran)
  })
}
