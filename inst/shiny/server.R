
shinyServer(function(input, output, session) {
  ## Construct the active MultiGSEAResultContainer object. This is used to feed
  ## different things that want to interact / visualze a MultiGSEAResult
  mgc <- reactive({
    req(input$mgresult)
    failWith(NULL, MultiGSEAResultContainer(input$mgresult$datapath))
  })

  gs_table_browser <- callModule(mgTableBrowser, 'mg_table_browser', mgc,
                                 server=TRUE)

  gs_viewer <- callModule(geneSetContrastView, 'geneset_viewer',
                          mgc, maxOptions=500, server=TRUE)

  observeEvent(gs_table_browser()$selected, {
    if (!is.null(gs_table_browser()$selected)) {
      sval <- gs_table_browser()$selected
      updateGeneSetContrastViewGeneSet(session, 'geneset_viewer',
                                       choices=isolate(mgc()$choices),
                                       selected=gs_table_browser()$selected,
                                       server=TRUE)
    }
  })

  output$gseaMethodSummary <- renderUI({
    obj <- failWith(NULL, expr=mgc(), silent=TRUE)
    if (!is(obj, 'MultiGSEAResultContainer')) {
      tags$p(style="font-weight: bold; color: red",
             "Upload the MultiGSEAResult object to initialize the application")
    } else {
      tagList(
        tags$h4("GSEA Analyses Overview"),
        summaryHTMLTable.multiGSEA(mgc()$mg, mgc()$methods,
                                   gs_table_browser()$fdr,
                                   p.col='padj.by.collection')
      )
    }
  })

})
