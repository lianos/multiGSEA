shinyServer(function(input, output, session) {
  ## If this application was invoked via explore(MultiGSEAResult), then
  ## getOption(EXPLORE_MULTIGSEA_RESULT='path/to/result.rds') was set that
  ## we can load, otherwise this will respond to a user upload.
  mgc <- reactive({
    ## Are we here because the user uploaded something, or did the user ask
    ## to `explore(MultiGSEAResult)`?
    # msg("wiring up mgc")
    # browser()
    if (is.null(input$mgresult)) {
      mg <- getOption('EXPLORE_MULTIGSEA_RESULT', NULL)
      res <- failWith(NULL, MultiGSEAResultContainer(mg), silent=TRUE)
      return(res)
    }
    ## User uploaded a file
    return(failWith(NULL, MultiGSEAResultContainer(input$mgresult$datapath)))
  })

  lfc <- reactive({
    req(mgc()$mg) %>%
      logFC %>%
      arrange(desc(logFC))
  })

  gs_result_filter <- callModule(mgResultFilter, 'mg_result_filter', mgc)

  ## Overview Tab ==============================================================
  output$gseaMethodSummary <- renderUI({
    obj <- failWith(NULL, expr=mgc(), silent=TRUE)
    if (!is(obj, 'MultiGSEAResultContainer')) {
      tags$p(style="font-weight: bold; color: red",
             "Upload the MultiGSEAResult object to initialize the application")
    } else {
      tagList(
        tags$h4("GSEA Analyses Overview"),
        summaryHTMLTable.multiGSEA(mgc()$mg, mgc()$methods,
                                   gs_result_filter()$fdr(),
                                   p.col='padj.by.collection')
      )
    }
  })

  ## GSEA Results Tab ==========================================================
  gs_viewer <- callModule(geneSetContrastView, 'geneset_viewer',
                          mgc, maxOptions=500, server=TRUE)

  ## A table of GSEA statistics/results for the given method and fdr threshold
  ## The table is wired to the gs_viewer so that row clicks can signal updates
  ## to the contrast viewer
  gs_table_browser <- callModule(mgTableBrowser, 'mg_table_browser', mgc,
                                 method=gs_result_filter()$method,
                                 fdr=gs_result_filter()$fdr,
                                 server=TRUE)
  ## clicks on gsea result table update the contrast view
  observeEvent(gs_table_browser$selected(), {
    .mgc <- req(mgc())
    geneset <- req(gs_table_browser$selected())
    updateActiveGeneSetInContrastView(session, gs_viewer, geneset, .mgc)
  })

  ## A table of other genesets that brushed genes in the contrast viewer
  ## belong to. This table is also wired to the contrast viewer, so that
  ## a click on a row of the table will update the contrast view, too.
  other_genesets_gsea <- callModule(mgGeneSetSummaryByGene,
                                    'other_genesets_gsea',
                                    mgc, features=gs_viewer()$selected,
                                    method=gs_result_filter()$method,
                                    fdr=gs_result_filter()$fdr)
  ## DEBUG: Can we add a DT row click listner to the `other_genesets_gsea` so
  ## that it updates the `gs_viewer`? My first shot at doing sends the
  ## application into a tailspin, my best guess is because the selection is
  ## still active in the rbokeh boxp/density plot.

  ## Differential Gene Expression Tab ==========================================
  gene.volcano <- callModule(mgVolcano, 'dge_volcano', mgc)

  output$dge_volcano_genestats <- DT::renderDataTable({
    res <- req(lfc()) %>%
      select(symbol, featureId, logFC, pval, padj)

    selected <- gene.volcano()
    if (!is.null(selected)) {
      res <- filter(res, featureId %in% selected$featureId)
    }

    renderFeatureStatsDataTable(res, filter='top', feature.link.fn=ncbi.entrez.link)
  })

  ## A table of other genesets that brushed genes in the contrast viewer
  ## belong to. This table is also wired to the contrast viewer, so that
  ## a click on a row of the table will update the contrast view, too.
  other_genesets_volcano <- callModule(mgGeneSetSummaryByGene,
                                       'other_genesets_volcano',
                                       mgc, features=gene.volcano,
                                       method=gs_result_filter()$method,
                                       fdr=gs_result_filter()$fdr)
})
