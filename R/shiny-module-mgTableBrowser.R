##' A module to explore the GSEA statistics generated from a multiGESA Run.
##'
##' @rdname mgTableBrowserModule
##' @export
##' @importFrom shiny tagList uiOutput NS
##' @importFrom DT dataTableOutput
##' @param id the shiny id of the UI module
##' @param mg The \code{MultiGSEAResultContainer} object
##' @param server whether or not certain things (options, datatables) be
##'   "server" side or not. This isn't thoroughly pushed throughout everything
mgTableBrowserUI <- function(id, mg=NULL, server=TRUE) {
  ns <- NS(id)

  tagList(
    uiOutput(ns("resultTableMessage")),
    DT::dataTableOutput(ns("gseaResultTable")))
}

##' @rdname mgTableBrowserModule
##'
##' @export
##' @importFrom shiny reactive req renderUI tags
##' @importFrom DT renderDataTable
##' @param input shiny server input object
##' @param output shiny server output object
##' @param session shiny server session object
##' @param mgc the \code{MultiGSEAResultContainer} object
##' @param a reactive that determines which method to explore from this result
##' @param a reactive that gives the fdr threshold to filter results in the
##'   table by.
mgTableBrowser <- function(input, output, session, mgc, method, fdr,
                           server=TRUE) {

  ## under the FDR threshold
  gsea.result.table <- reactive({
    # browser()
    mg <- req(mgc()$mg)
    if (is.null(method()) || method() == "") {
      # msg("... gseaMethod not selected yet")
      return(NULL)
    }
    ## MultiGSEResult object, method, and FDR thersholds all set, now fetch
    ## the data that corresponds to this criteria
    constructGseaResultTable(mg, method(), fdr())
  })

  output$resultTableMessage <- renderUI({
    gst <- req(gsea.result.table())
    if (!is(gst, 'data.frame')) {
      tmsg <- ''
    } else if (nrow(gst) == 0) {
      tmsg <- sprintf('No results at FDR cutoff of %.2f', fdr())
    } else {
      tmsg <- sprintf('Showing %d genesets at FDR cutoff of %.2f',
                      nrow(gst), fdr())
    }
    tags$h5(tmsg)
  })

  output$gseaResultTable <- DT::renderDataTable({
    req(gsea.result.table(), mgc())
    renderGseaResultTableDataTable(gsea.result.table(), method(),
                                   mgc()$mg)
  }, server=server)

  reactive({
    if (!is.null(input$gseaResultTable_row_last_clicked)) {
      idx <- input$gseaResultTable_row_last_clicked
      xcol <- as.character(gsea.result.table()$collection[idx])
      xname <- as.character(gsea.result.table()$name[idx])
      selected <- paste(xcol, xname, sep='_::_')
    } else {
      selected <- NULL
    }

    list(stats=gsea.result.table, selected=selected)
  })
}
