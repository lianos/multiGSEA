##' Shiny module to pick GSEA method and fdr params used for display
##' 
##' @export
##' @rdname mgResultFilter
mgResultFilterUI <- function(id, mg=NULL) {
  stopifnot(requireNamespace('shiny'))
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
          style="padding-top: 2em;",
          shiny::downloadButton(ns("gseaDownloadStats"), 'Download'))))
  )
}

##' @export
##' @rdname mgResultFilter
mgResultFilter <- function(input, output, session, mgc) {
  stopifnot(requireNamespace('shiny'))

  ## When the MultiGSEAResult changes, we want to update different aspects of
  ## the application
  shiny::observeEvent(mgc(), {
    shiny::req(mgc())
    shiny::updateSelectInput(session, "gseaMethod",
                             choices=mgc()$methods,
                             selected=mgc()$methods[1L])
  })

  output$gseaDownloadStats <- shiny::downloadHandler(
    filename=function() {
      sprintf('multiGSEA-gsea-statistics-%s.csv', isolate(input$gseaMethod))
    },
    content=function(file) {
      write.csv(result(mgc()$mg, isolate(input$gseaMethod)), file,
                row.names=FALSE)
    }
  )

  shiny::reactive({
    list(method=reactive(input$gseaMethod), fdr=reactive(input$gseaReportFDR))
  })
}
