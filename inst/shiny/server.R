
shinyServer(function(input, output, session) {
  ## Construct the active MultiGSEAResultContainer object. This is used to feed
  ## different things that want to interact / visualze a MultiGSEAResult
  mgc <- reactive({
    req(input$mgresult)
    failWith(NULL, MultiGSEAResultContainer(input$mgresult$datapath))
  })

  ## genesetView ===============================================================
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

  ## geneView ==================================================================
  gene.volcano <- callModule(mgVolcano, 'dge_volcano', mgc)
  output$dge_volcano_stats <- DT::renderDataTable({
    genes <- setDF(req(gene.volcano()))
    ## Fetch the genesets that have `genes` in them and show them here.
    mg <- mgc()$mg
    gdb <- geneSetDb(mg)
    browser()
    xgdb <- subsetByFeatures(gdb, genes$featureId)
    sub.gs <- geneSets(xgdb) %>%
      select(collection, name, n, n.sig, logFC=mean.logFC.trim)
    lfc <- logFC(mg) %>%
      setDF %>%
      semi_join(genes, by='featureId') %>%
      select(featureId, symbol, logFC)
    lfc.wide <- t(setNames(lfc$logFC, lfc$symbol))[rep(1, nrow(sub.gs)),,drop=FALSE]
    lfc.wide <- lfc.wide[, order(colnames(lfc.wide)),drop=FALSE]
    out <- bind_cols(sub.gs, as.data.frame(lfc.wide))
    datatable(out) %>% roundDT
  })
})

