
shinyServer(function(input, output, session) {
  vals <- reactiveValues(
    geneset=list(),
    gs.name=character(),
    gsea.result.table=NULL,
    genes=character())

  ## Construct the active MultiGSEAResult object. Currently this can only be
  ## set by the user uploading the MultiGSEAResult, but I want to update this
  ## such that the shiny app can be launched with the MultiGSEAResult object
  ## in hand, either because:
  ##   (1) The workspace that launched it already has it available, ie. the user
  ##       called explore(mg); or
  ##   (2) The app is pointed to a path on the filesystem where the
  ##       MultiGSEAResult object is already saved.
  mg <- reactive({
    infile <- input$mgresult
    if (is.null(infile)) {
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
    out
  })

  ## Don't put side effects in your reactives, or Joe Cheng will kill you
  ## When a new MultiGSEAResult (mg) object is set, we need to update/reset some
  ## parts of the UI, including:
  ## (1) Picking which GSEA method is the default to show results for first. I
  ##     prefer to pick `camera` as the default method if it was run.
  observeEvent(mg(), {
    if (!is(mg(), 'MultiGSEAResult')) {
      return(NULL)
    }
    methods <- resultNames(mg())
    method <- if ('camera' %in% methods) 'camera' else methods[1L]
    updateSelectInput(session, "gseaMethod", choices=methods, selected=method)
    # reset any previously picked GSEA result
    vals$genelist <- list()
    vals$gs.name <- character()
  })

  ## Holds the current table of "active" GSEA statistics that are browseable.
  ## These are determined by:
  ##   (1) The GSEA method the user has chosen to browse via the gseaMethod
  ##       selectInput; and
  ##   (2) The maximum FDR selected to show
  gsea.result.table <- reactive({
    if (!is(mg(), 'MultiGSEAResult')) {
      message("mg() is not MultiGSEAResult in reactive gsea.result.table")
      return(NULL)
    }
    if (input$gseaMethod == "") {
      message("... gseaMethod not selected yet")
      return(NULL)
    }

    ## MultiGSEResult object, method, and FDR thersholds all set, now fetch
    ## the data that corresponds to this criteria
    constructGseaResultTable(mg(), input$gseaMethod, input$gseaReportFDR)
  })

  output$gseaResultLabel <- renderUI({
    if (is.null(input$gseaMethod) || input$gseaMethod == '') {
      ''
    } else {
      input$gseaMethod
    }
  })

  output$resultTableMessage <- renderUI({
    req(gst <- gsea.result.table())
    if (!is(gst, 'data.frame')) {
      msg <- ''
    } else if (nrow(gst) == 0) {
      msg <- sprintf('No results at FDR cutoff of %.2f',
                     input$gseaReportFDR)
    } else {
      msg <- sprintf('Showing %d genesets at FDR cutoff of %.2f',
                     nrow(gst), input$gseaReportFDR)
    }
    h5(msg)
  })

  output$gseaResultTable <- DT::renderDataTable({
    req(gsea.result.table())
    req(mg())
    dtargs <- prepareRenderGseaResultTable(gsea.result.table(), mg())
    do.call(datatable, dtargs)
  })

  output$gseaMethodSummary <- renderUI({
    if (!is(mg(), 'MultiGSEAResult')) {
      return(NULL)
    }
    summaryHTMLTable.multiGSEA(mg(), resultNames(mg()), input$gseaReportFDR,
                               p.col='padj.by.collection')
  })

  ## When the user clicks a row in the geneset statistics datatable, we want
  ## to update the current geneset to the row that was clicked.
  observeEvent(input$gseaResultTable_row_last_clicked, {
    idx <- input$gseaResultTable_row_last_clicked
    xcol <- as.character(gsea.result.table()$collection[idx])
    xname <- gsea.result.table()$name[idx]
    xstats <- arrange_(geneSet(mg(), xcol, xname), ~ -logFC)
    vals$geneset <- list(collection=xcol, name=xname, stats=xstats)
  })

  ## Lists the differential expression statistics for the active geneset
  output$genesetView_genesetMembers <- DT::renderDataTable({
    req(mg())
    if (length(vals$geneset) == 0) return(NULL)

    gcols <- c('symbol', 'featureId', 'logFC', 'padj')
    gs <- vals$geneset$stats[, gcols, with=FALSE]
    gs <- round.dt(gs)
    setnames(gs, 'padj', 'FDR')
    if (nrow(gs) < 10) {
      length.opts <- nrow(gs)
    } else {
      length.opts <- c(5, 10, 15, 50, 100)
      length.opts <- length.opts[length.opts < nrow(gs)]
      length.opts <- c(length.opts, nrow(gs))
    }

    dt.opts <- list(
      pageLength=length.opts[1L],
      lengthMenu=length.opts,
      dom='lBtipr',
      buttons=c('copy', 'csv', 'excel'))

    datatable(gs, filter='top', extension='Buttons', rownames=FALSE,
              selection='none', options=dt.opts)
  })

  ## Creates the geneset enrichment plot
  output$genesetView_gseaPlot <- renderPlotly({
    req(mg())
    if (length(vals$geneset) == 0) return(NULL)
    xcol <- vals$geneset$collection[1L]
    xname <- vals$geneset$name[1L]
    iplt <- iplot(mg(), xcol, xname,
                  value=input$gensetView_plot_statistic,
                  type=input$genesetView_plot_type)
  })

})
