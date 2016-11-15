##' A module to explore the GSEA statistics generated from a multiGESA Run.
##'
##' @export
##' @param id the shiny id of the UI module
##' @param mg The \code{MultiGSEAResultContainer} object
##' @param server whether or not certain things (options, datatables) be
##'   "server" side or not. This isn't thoroughly pushed throughout everything
##' @rdname mgTableBrowserModule
mgTableBrowserUI <- function(id, mg=NULL, server=TRUE) {
  stopifnot(requireNamespace('shiny'),
            requireNamespace("DT"))

  ns <- shiny::NS(id)

  shiny::tagList(
    shiny::uiOutput(ns("resultTableMessage")),
    DT::dataTableOutput(ns("gseaResultTable")))
}

##' @export
##' @param input shiny server input object
##' @param output shiny server output object
##' @param session shiny server session object
##' @param mgc the \code{MultiGSEAResultContainer} object
##' @param a reactive that determines which method to explore from this result
##' @param a reactive that gives the fdr threshold to filter results in the
##'   table by.
##' @rdname mgTableBrowserModule
mgTableBrowser <- function(input, output, session, mgc, method, fdr,
                           server=TRUE) {
  stopifnot(requireNamespace('shiny'), requireNamespace('DT'))

  ## under the FDR threshold
  gsea.result.table <- shiny::reactive({
    # browser()
    mg <- shiny::req(mgc()$mg)
    if (is.null(method()) || method() == "") {
      msg("... gseaMethod not selected yet")
      return(NULL)
    }
    ## MultiGSEResult object, method, and FDR thersholds all set, now fetch
    ## the data that corresponds to this criteria
    constructGseaResultTable(mg, method(), fdr())
  })

  output$resultTableMessage <- shiny::renderUI({
    gst <- shiny::req(gsea.result.table())
    if (!is(gst, 'data.frame')) {
      tmsg <- ''
    } else if (nrow(gst) == 0) {
      tmsg <- sprintf('No results at FDR cutoff of %.2f', fdr())
    } else {
      tmsg <- sprintf('Showing %d genesets at FDR cutoff of %.2f',
                      nrow(gst), fdr())
    }
    shiny::tags$h5(tmsg)
  })

  output$gseaResultTable <- DT::renderDataTable({
    shiny::req(gsea.result.table(), mgc())
    renderGseaResultTableDataTable(gsea.result.table(), method(),
                                   mgc()$mg)
  }, server=server)

  shiny::reactive({
    if (!is.null(input$gseaResultTable_row_last_clicked)) {
      idx <- input$gseaResultTable_row_last_clicked
      xcol <- as.character(gsea.result.table()$collection[idx])
      xname <- as.character(gsea.result.table()$name[idx])
      selected <- paste(xcol, xname, sep='_::_')
      msg("Selected: ", selected)
    } else {
      selected <- NULL
    }

    list(stats=gsea.result.table, selected=selected)
  })
}
