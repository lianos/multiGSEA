
shinyServer(function(input, output, session) {
  vals <- reactiveValues(
    geneset=list(),
    gs.name=character(),
    gsea.result.table=NULL,
    genes=character())

  ## The active MultiGSEAResult object
  mg <- reactive({
    infile <- input$mgresult
    if (is.null(infile)) {
      message("NULL mgresult")
      return(NULL)
    }

    out <- readRDS(infile$datapath)
    if (!is(out, 'MultiGSEAResult')) {
      message("File is not a MultiGSEAResult")
      return(NULL)
    }

    methods <- resultNames(out)
    if (length(methods) == 0L) {
      message("No GSEA methods found in MultiGSEAResult")
      return(NULL)
    }
    method <- if ('camera' %in% methods) 'camera' else methods[1L]
    updateSelectInput(session, "gseaMethod", choices=methods, selected=method)
    out
  })

  gsea.result.table <- reactive({
    method <- input$gseaMethod
    if (is.null(mg())) {
      message("NULL mg() in reactive gsea.result.table")
      return(NULL)
    }
    if (method == "") {
      message("... gseaMethod not selected yet")
      return(NULL)
    }

    ## Method is selected: let's act
    method <- input$gseaMethod
    message("... updating result.table with method: ", method)
    output$gseaResultName <- renderUI(h4(sprintf("%s summary", method)))
    subset(result(mg(), method), padj.by.collection <= input$gseaReportFDR)
  })

  observeEvent(gsea.result.table(), {
    if (is.null(gsea.result.table())) {
      message("NULL gsea.result.table() in observeEvent(gsea.result.table())")
      return(NULL)
    }
    output$gseaResultTable <- render.dt(gsea.result.table(), mg())
  })

  output$gseaMethodSummary <- renderUI({
    summaryHTMLTable.multiGSEA(mg(), resultNames(mg()), input$gseaReportFDR,
                               p.col='padj.by.collection')

  })

  observeEvent(input$gseaResultTable_row_last_clicked, {
    idx <- input$gseaResultTable_row_last_clicked
    xcol <- gsea.result.table()$collection[idx]
    xname <- gsea.result.table()$name[idx]
    vals$geneset <- list(collection=xcol, name=xname, stats=geneSet(mg(), xcol, xname))
  })

  ## Update Plot
  ## Update logFC table
  observeEvent(vals$geneset, {
    if (!all(c('collection', 'name', 'stats') %in% names(vals$geneset))) {
      return(NULL)
    }
    xcol <- vals$geneset$collection[1L]
    xname <- vals$geneset$name[1L]

    ## Plot
    iplt <- iplot(mg(), xcol, xname, type=input$genesetView_plot_type)
    output$genesetView_gseaPlot <- renderPlotly(iplt)

    ## Update members table
    gcols <- c('symbol', 'featureId', 'logFC', 'padj')
    gs <- vals$geneset$stats[, gcols, with=FALSE]
    gs <- round.dt(gs)
    output$genesetView_genesetMembers <- DT::renderDataTable(
      gs, filter='top', rownames=FALSE, selection='none',
      options=list(
        pageLength=5,
        lengthMenu=c(5, 10, 15, 20),
        dom='ltipr'
      ))
  })
})
