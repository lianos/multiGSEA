
##' @export
mgGeneSetSummaryByGeneUI <- function(id, mg=NULL) {
  stopifnot(requireNamespace('shiny'))
  ns <- shiny::NS(id)

  shiny::tagList(
    checkboxInput(ns('genesets_sigonly'),
                  'Show only significant gene set membership',
                  value=FALSE),
    DT::dataTableOutput(ns("other_genesets")))
}

##' @export
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

    shiny::req(mgc())
    fids <- shiny::req(features())
    if (is(fids, 'data.frame')) {
      fids <- fids$featureId
    }

    geneSetSummaryByGenes(mgc()$mg, fids, feature.rename='symbol',
                          method=method, max.p=max.p, .external=FALSE)
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
    datatable(out, filter='top', escape=FALSE,
              selection=list(mode='single', selected=NA, target='row'),
              rownames=FALSE)
  })

  ## Return the selected geneset
  shiny::reactive({
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
