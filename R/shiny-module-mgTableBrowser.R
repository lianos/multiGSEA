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
    shiny::fluidRow(
      shiny::column(
        3,
        shiny::selectInput(ns("gseaMethod"), "GSEA Methods", "")),
      column(
        6,
        shiny::sliderInput(ns("gseaReportFDR"), "FDR Cutoff", min=0,
                           max=1, value=0.2, step=0.05)),
      shiny::column(
        3,
        shiny::tags$div(
          style="padding-top: 3.5em;",
          shiny::downloadButton(ns("gseaDownloadStats"), 'Download')))),
    shiny::uiOutput(ns("resultTableMessage")),
    DT::dataTableOutput(ns("gseaResultTable")))
}

##' @export
##' @param input shiny server input object
##' @param output shiny server output object
##' @param session shiny server session object
##' @param mgc the \code{MultiGSEAResultContainer} object
##' @rdname mgTableBrowserModule
mgTableBrowser <- function(input, output, session, mgc, server=TRUE) {
  stopifnot(requireNamespace('shiny'), requireNamespace('DT'))
  ## When the MultiGSEAResult changes, we want to update different aspects of
  ## the application
  shiny::observeEvent(mgc(), {
    shiny::req(mgc())
    shiny::updateSelectInput(session, "gseaMethod",
                             choices=mgc()$methods,
                             selected=mgc()$methods[1L])
  })

  ## under the FDR threshold
  gsea.result.table <- shiny::reactive({
    shiny::req(mgc())
    if (input$gseaMethod == "") {
      msg("... gseaMethod not selected yet")
      return(NULL)
    }
    ## MultiGSEResult object, method, and FDR thersholds all set, now fetch
    ## the data that corresponds to this criteria
    constructGseaResultTable(mgc()$mg, input$gseaMethod, input$gseaReportFDR)
  })

  output$resultTableMessage <- shiny::renderUI({
    gst <- shiny::req(gsea.result.table())
    if (!is(gst, 'data.frame')) {
      tmsg <- ''
    } else if (nrow(gst) == 0) {
      tmsg <- sprintf('No results at FDR cutoff of %.2f',
                      input$gseaReportFDR)
    } else {
      tmsg <- sprintf('Showing %d genesets at FDR cutoff of %.2f',
                      nrow(gst), input$gseaReportFDR)
    }
    shiny::tags$h5(tmsg)
  })

  output$gseaResultTable <- DT::renderDataTable({
    shiny::req(gsea.result.table(), mgc())
    renderGseaResultTableDataTable(gsea.result.table(), input$gseaMethod,
                                   mgc()$mg)
  }, server=server)

  output$gseaDownloadStats <- shiny::downloadHandler(
    filename=function() {
      sprintf('multiGSEA-gsea-statistics-%s.csv', input$gseaMethod)
    },
    content=function(file) {
      write.csv(result(mgc()$mg, input$gseaMethod), file, row.names=FALSE)
    }
  )

  shiny::reactive({
    if (!is.null(input$gseaResultTable_row_last_clicked)) {
      idx <- input$gseaResultTable_row_last_clicked
      xcol <- as.character(gsea.result.table()$collection[idx])
      xname <- as.character(gsea.result.table()$name[idx])
      selected <- paste(xcol, xname, sep='_::_')
    } else {
      selected <- NULL
    }

    list(stats=gsea.result.table,
         method=input$gseaMethod,
         fdr=input$gseaReportFDR,
         selected=selected)
  })
}
