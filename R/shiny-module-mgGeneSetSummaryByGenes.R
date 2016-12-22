
##' Module that displays gene sets related to (by membership) a set of genes.
##' @export
##' @importFrom shiny NS tagList checkboxInput uiOutput
##' @importFrom DT dataTableOutput
##' @rdname mgGeneSetSummaryByGene
mgGeneSetSummaryByGeneUI <- function(id, mg=NULL) {
  ns <- NS(id)

  tagList(
    checkboxInput(ns('genesets_sigonly'),
                  'Show membership for significant gene sets only',
                  value=FALSE),
    uiOutput(ns('selected_message')),
    DT::dataTableOutput(ns("other_genesets")))
}

##' @rdname mgGeneSetSummaryByGene
##' @export
##' @importFrom shiny reactive req renderUI tags
##' @importFrom DT renderDataTable datatable
mgGeneSetSummaryByGene <- function(input, output, session, mgc,
                                   features, method, fdr) {
  genesets <- reactive({
    if (input$genesets_sigonly) {
      method <- method()
      max.p <- fdr()
    } else {
      method <- NULL
      max.p <- NULL
    }
    req(mgc())
    fids <- req(features())
    if (is(fids, 'data.frame')) {
      fids <- fids$featureId
    }
    geneSetSummaryByGenes(mgc()$mg, fids, feature.rename='symbol',
                          method=method, max.p=max.p, .external=FALSE)
  })

  output$selected_message <- renderUI({
    fids <- features()
    if (is.null(fids)) {
      n <- 0L
      ngs <- 0L
    } else {
      n <- if (is.vector(fids)) length(fids) else nrow(fids)
      ngs <- nrow(genesets())
    }
    tags$p(sprintf('%d features selected across %d genesets', n, ngs))
  })

  output$other_genesets <- DT::renderDataTable({
    out <- copy(req(genesets()))
    mg <- req(mgc()$mg)
    out[, collection := factor(collection)]
    out[, active := NULL]
    out[, name := {
      url <- geneSetURL(mg, as.character(collection), name)
      xname <- gsub('_', ' ', name)
      html <- '<a href="%s" target="_blank">%s</a>'
      ifelse(is.na(url), xname, sprintf(html, url, xname))
    }]

    out <- round.dt(out)
    DT::datatable(out, filter='top', escape=FALSE,
                  selection=list(mode='single', selected=NA, target='row'),
                  rownames=FALSE)
  })

  ## Return the selected geneset
  reactive({
    idx <- input$other_genesets_row_last_clicked
    if (!is.null(idx)) {
      others <- genesets()
      xcol <- as.character(others$collection[idx])
      xname <- as.character(others$name[idx])
      selected <- paste(xcol, xname, sep='_::_')
      msg("Selected: ", selected)
    } else {
      selected <- NULL
    }
    list(others=genesets, selected=selected)
  })
}
