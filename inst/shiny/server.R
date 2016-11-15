
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
  observeEvent(gs_table_browser()$selected, {
    sval <- gs_table_browser()$selected
    if (!is.null(sval)) {
      updateGeneSetContrastViewGeneSet(session, 'geneset_viewer',
                                       choices=isolate(mgc()$choices),
                                       selected=sval, server=TRUE)
    }
  })

  ## A table of other genesets that brushed genes in the contrast viewer
  ## belong to. This table is also wired to the contrast viewer, so that
  ## a click on a row of the table will update the contrast view, too.
  other_genesets_gsea <- callModule(mgGeneSetSummaryByGene,
                                    'other_genesets_gsea',
                                    mgc, features=gs_viewer()$selected,
                                    method=gs_result_filter()$method,
                                    fdr=gs_result_filter()$fdr)
  # observeEvent(other_genesets_gsea()$selected, {
  #   ## This sends us into a tailspin of updating the plot, then updating
  #   ## the selected genes and repupdating plot, etc ... it shouldn't, but ...
  #   sval <- other_genesets_gsea()$selected
  #   if (!is.null(sval)) {
  #     updateGeneSetContrastViewGeneSet(session, 'geneset_viewer',
  #                                      choices=isolate(mgc()$choices),
  #                                      selected=sval, server=TRUE)
  #   }
  # })

  ## Differential Gene Expression Tab ==========================================
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

  ## A table of other genesets that brushed genes in the contrast viewer
  ## belong to. This table is also wired to the contrast viewer, so that
  ## a click on a row of the table will update the contrast view, too.
  other_genesets_volcano <- callModule(mgGeneSetSummaryByGene,
                                       'other_genesets_volcano',
                                       mgc, features=gene.volcano,
                                       method=gs_result_filter()$method,
                                       fdr=gs_result_filter()$fdr)
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
