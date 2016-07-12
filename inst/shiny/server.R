
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
    if (!is(mg(), 'MultiGSEAResult')) {
      message("mg() is not MultiGSEAResult in reactive gsea.result.table")
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
    out <- subset(result(mg(), method), padj.by.collection <= input$gseaReportFDR)

    if (nrow(out)) {
      has.hallmark <- 'h' %in% out$collection
      if (has.hallmark) {
        lvls <- c('h', sort(setdiff(out$collection, 'h')))
        out[, collection := factor(collection, lvls, ordered=TRUE)]
        out <- arrange_(out, ~collection, ~ -mean.logFC.trim)
      } else {
        lvls <- sort(unique(out$collection))
        out[, collection := factor(collection, lvls)]
        out <- arrange_(out, ~ -mean.logFC.trim)
      }
      # Trigger refresh of "active" geneset to update graphs and members
      xcol <- as.character(out$collection[1L])
      xname <- out$name[1L]
      xstats <- arrange_(geneSet(mg(), xcol, xname), ~ -logFC)
      vals$geneset <- list(collection=xcol, name=xname, stats=xstats)
    } else {
      out <- NULL
    }

    out
  })

  observeEvent(gsea.result.table(), {
    if (is.null(gsea.result.table())) {
      message("NULL gsea.result.table() in observeEvent(gsea.result.table())")
      return(NULL)
    }
    output$gseaResultTable <- render.dt(gsea.result.table(), mg())
  })

  output$gseaMethodSummary <- renderUI({
    if (!is(mg(), 'MultiGSEAResult')) {
      return(NULL)
    }
    summaryHTMLTable.multiGSEA(mg(), resultNames(mg()), input$gseaReportFDR,
                               p.col='padj.by.collection')
  })

  observeEvent(input$gseaResultTable_row_last_clicked, {
    idx <- input$gseaResultTable_row_last_clicked
    xcol <- as.character(gsea.result.table()$collection[idx])
    xname <- gsea.result.table()$name[idx]
    xstats <- arrange_(geneSet(mg(), xcol, xname), ~ -logFC)
    vals$geneset <- list(collection=xcol, name=xname, stats=xstats)
  })

  observeEvent(vals$geneset, {
    ## Whenever our 'active geneset' is updated, we want to update the
    ## geneset plot as well as the geneset table with dge statistics
    if (!all(c('collection', 'name', 'stats') %in% names(vals$geneset))) {
      return(NULL)
    }
    xcol <- vals$geneset$collection[1L]
    xname <- vals$geneset$name[1L]

    ## Update GSEA Plot
    iplt <- iplot(mg(), xcol, xname, value=input$gensetView_plot_statistic,
                  type=input$genesetView_plot_type)
    output$genesetView_gseaPlot <- renderPlotly(iplt)

    ## Update members table
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
    output$genesetView_genesetMembers <- DT::renderDataTable(
      gs, filter='top', extension='Buttons', rownames=FALSE, selection='none',
      options=list(
        pageLength=length.opts[1L],
        lengthMenu=length.opts,
        dom='lBtipr',
        buttons=c('copy', 'csv', 'excel')))
  })

  observeEvent(list(input$genesetView_plot_type, input$gensetView_plot_statistic), {
    if (!is(mg(), 'MultiGSEAResult')) {
      return(NULL)
    }
    xcol <- vals$geneset$collection[1L]
    xname <- vals$geneset$name[1L]
    iplt <- iplot(mg(), xcol, xname,
                  value=input$gensetView_plot_statistic,
                  type=input$genesetView_plot_type)
    output$genesetView_gseaPlot <- renderPlotly(iplt)
  })
})
