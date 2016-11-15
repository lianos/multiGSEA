
shinyServer(function(input, output, session) {
  ## Construct the active MultiGSEAResultContainer object. This is used to feed
  ## different things that want to interact / visualze a MultiGSEAResult
  mgc <- reactive({
    req(input$mgresult)
    failWith(NULL, MultiGSEAResultContainer(input$mgresult$datapath))
  })

  species <- reactive({
    ## if symbols == touppser(symbol) this is human, else it's a mouse (GNE only)
    req(mgc())
    mg <- mgc()$mg
    symbols <- logFC(mg)$symbol
    if (is.character(symbols)) {
      symbols <- symbols[!is.na(symbols)]
      tests <- sample(symbols, min(20, length(symbols)))
      ans <- if (mean(tests == toupper(tests)) > 0.7) 'human' else 'mouse'
    } else {
      ans <- NULL
    }
    ans
  })

  gs_result_filter <- callModule(mgResultFilter, 'mg_result_filter', mgc)

  ## genesetView ===============================================================
  gs_table_browser <- callModule(mgTableBrowser, 'mg_table_browser', mgc,
                                 method=gs_result_filter()$method,
                                 fdr=gs_result_filter()$fdr,
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
                                   gs_result_filter()$fdr(),
                                   p.col='padj.by.collection')
      )
    }
  })

  ## geneView ==================================================================
  gene.volcano <- callModule(mgVolcano, 'dge_volcano', mgc)

  output$dge_volcano_genestats <- DT::renderDataTable({
    genes <- setDF(req(gene.volcano()))
    mg <- mgc()$mg
    res <- logFC(mg) %>%
      filter(featureId %in% genes$featureId) %>%
      select(symbol, featureId, pval, padj)
    if (!is.null(species())) {
      url <- sprintf('http://research.gene.com/genehub/#/summary/gene/%s/%s',
                     res$featureId, species())
      html <- sprintf('<a href="%s" target="_blank">%s</a>', url,
                      res$symbol)
      res$symbol <- html
    }
    datatable(res, filter='top', escape=FALSE) %>% roundDT
  })

  output$dge_volcano_genesetstats <- DT::renderDataTable({
    genes <- setDF(req(gene.volcano()))
    mg <- mgc()$mg

    if (input$dge_genesets_sigonly) {
      method <- gs_result_filter()$method()
      max.p <- gs_result_filter()$fdr()
    } else {
      method <- NULL
      max.p <- NULL
    }
    out <- geneSetSummaryByGenes(mg, genes$featureId, feature.rename='symbol',
                                 method=method, max.p=max.p, .external=FALSE)
    out <- round.dt(out)
    datatable(out, filter='top')
  })
})

if (FALSE) {
  genes <- tibble(
    featureId=c("102633704", "12565", "242384", "19220", "72780",
                "102636809", "319555"))
  feature.rename <- 'symbol'
  method <- 'camera'
  max.p <- 0.20
  mg <- readRDS("/Users/lianogls/tmp/schmidt/multiGSEA-EP-uber_hWT_tKO-hWT_tWT.rds")
  out <- geneSetSummaryByGenes(mg, genes$featureId, feature.rename='symbol',
                               method=method, max.p=max.p)

}
