##' Interactively explore a MultiGSEAResult via a shiny app
##'
##' @description
##' This will launch a shiny application that enables exploratory analysis of
##' the results from the different GSEA methods run within a
##' \code{MultiGSEAResult}.
##'
##' Reference the "shiny-multiGSEA" vignette for more detailed documentation of
##' the functionality provided by this application.
##'
##' @export
##' @importFrom shiny runApp
##' @param x A \code{MultiGSEAResult} object
##' @examples
##' \dontrun{
##' vm <- exampleExpressionSet()
##' gdb <- exampleGeneSetDb()
##' mg <- multiGSEA(gdb, vm, vm$design, methods=c('camera', 'fry'))
##' explore(mg)
##' }
explore <- function(x) {
  stopifnot(is(x, 'MultiGSEAResult'))
  options(EXPLORE_MULTIGSEA_RESULT=x)
  on.exit(options(EXPLORE_MULTIGSEA_RESULT=NULL))
  runApp(system.file('shiny', package='multiGSEA'))
}

## Gene Set Level Table Helpers ================================================

##' Builds the table of GSEA statistics to present to the user
##'
##' @description
##' The application will present the set of gene sets that pass a given
##' \code{fdr} for a given \code{method} as a central piece of the UI. This
##' function accepts those to arguments and prepares the statistics generated
##' from the analysis for display.
##'
##' @export
##' @param mg \code{MultiGSEAResult} object
##' @param method the method to show statistics for
##' @param fdr the FDR cut off to present statistics for
##' @param prioritize the preffered collections to put at the top of the
##'   list. The collection column of the table is turned into a factor and for
##'   more usful display with datatable's filter. Specifcying collections
##'   here will put those collections at the front of the factor levels and
##'   therofre prioritize their display in the select dropdown for the filter
##' @return a data.table of the statistics that match the filtering criteria.
##'   A 0-row data.table is returned if nothing passes.
constructGseaResultTable <- function(mg, method, fdr, prioritize=c('h')) {
  out <- result(mg, method, .external=FALSE)
  out <- out[padj.by.collection <= fdr]
  if (nrow(out)) {
    colls <- sort(unique(out$collection))
    priority <- intersect(prioritize, colls)
    lvls <- c(priority, setdiff(colls, priority))
    out[, collection := factor(collection, lvls, ordered=TRUE)]
    out <- out[order(mean.logFC.trim, decreasing=TRUE)]
  }
  out
}

##' Creates a DT::datatalbe of geneset level GSEA results for use in shiny bits
##'
##' @export
##' @importFrom DT datatable
##' @param x The set of GSEA statistics generated from from
##'   \code{\link{constructGseaResultTable}}
##' @param method the GSEA method being used fo rdisplay
##' @param mg The \code{MultiGSEAResult} object. This is used swap in the
##' URL links for genesets using \code{\link{geneSetURL}}.
##' @return a DT::DataTable
renderGseaResultTableDataTable <- function(x, method, mg, digits=3) {
  stopifnot(is(x, 'data.table'))
  stopifnot(is.character(method) && length(method) == 1L)
  stopifnot(is(mg, 'MultiGSEAResult'))

  rcols <- c('collection'='collection', 'name'='name', 'n'='n',
             ## 'padj'='padj', 'padj.by.collection'='padjByColl', 'pval'='pval',
             'padj.by.collection'='FDR',
             'mean.logFC.trim'='logFC', 'n.sig.up'='nSigUp',
             'n.sig.down'='nSigDown', 'n.up'='nUp', 'n.down'='nDown')

  res <- x[, names(rcols), with=FALSE]
  setnames(res, names(rcols), rcols)

  res[, name := {
    url <- geneSetURL(mg, as.character(collection), name)
    xname <- gsub('_', ' ', name)
    html <- '<a href="%s" target="_blank">%s</a>'
    ifelse(is.na(url), xname, sprintf(html, url, xname))
  }]

  dt.order <- list()
  lfc.col <- which(colnames(res) == 'logFC') - 1L
  dt.order[[length(dt.order) + 1L]] <- list(lfc.col, 'desc')

  length.opts <- c(10, 20, 50, 100, 250)
  length.opts <- length.opts[length.opts < nrow(res)]
  length.opts <- c(length.opts, nrow(res))

  dt.opts <- list(
    dom='ltpir',
    order=dt.order,
    scrollX=TRUE,
    pageLength=length.opts[1L],
    lengthMenu=length.opts)
  dtargs <- list(data=res, filter='top',
                 selection=list(mode='single', selected=NA, target='row'),
                 # extensions='Buttons',
                 escape=FALSE, rownames=FALSE,
                 options=dt.opts)
  do.call(DT::datatable, dtargs) %>% roundDT(digits=digits)
}

## Gene evel Table Helpers =====================================================

##' Transforms a column in feature table to an external link for that feature.
##'
##' @description
##' When listing features in an interactive table, it's often useful to link
##' the feature to an external webpage that has more information about that
##' feature. Functions to genes to their NCBI or GeneCards webpage via their
##' \code{featureId} are provided via \code{ncbi.entrez.link} and
##' \code{genecards.entrez.link}. The column used to transform into a link
##' is specified by \code{link.col}.
##'
##' If \code{link.col} is not found in the data.frame \code{x} then the provided
##' functions are NO-OPS, ie. the same data.frame is simply returned.
##'
##' @rdname feature-link-functions
##' @export
##' @param x a data.frame from \code{logFC(MultiGSEAResult)}
##' @param link.col the column in \code{x} that should be transformed to a link
##' @return a modified \code{x} with an html link in \code{link.col}
ncbi.entrez.link <- function(x, link.col='symbol') {
  if (is.character(link.col) && is.character(x[[link.col]])) {
    url <- sprintf('https://www.ncbi.nlm.nih.gov/gene/%s', x$featureId)
    html <- sprintf('<a href="%s" target="_blank">%s</a>', url, x$symbol)
    x[[link.col]] <- html
  }
  x
}


##' @rdname feature-link-functions
##' @export
genecards.entrez.link <- function(x, link.col='symbol') {
  if (is.character(link.col) && is.character(x[[link.col]])) {
    url <- sprintf('http://www.genecards.org/cgi-bin/carddisp.pl?gene=%s',
                   x$featureId)
    html <- sprintf('<a href="%s" target="_blank">%s</a>', url, x$symbol)
    x[[link.col]] <- html
  }
  x
}

##' Creates a DT::datatable of feature level statistics for use in shiny bits
##'
##' We often want to display an interactive table of feature level statistics
##' for all features. This function is a convenience wrapper to do that.
##'
##' @export
##' @param x A \code{MultiGSEAResult} or \code{data.frame} of feature level
##'   statistics. When \code{x} is a \code{\link{MultiGSEAResult}}, the
##'   \code{\link{logFC}} feature level statistics will be extracted for
##'   display.
##' @param features A character vector that specifies the subset of
##'   \code{featureId}'s to display from \code{x}. If \code{NULL} (default),
##'   all of \code{x} will be used.
##' @param digits number of digits to round the numeric columns to
##' @param columns the columns from \code{x} to use. If \code{missing}, then
##'   only \code{c('symbol', 'featureId', 'logFC', 'pval', 'padj', order.by)}
##'   will be used. If explicitly set to \code{NULL} all columns will be used.
##'
##' @param feature.link.fn A funcion that receives the data.frame of statistics
##'   to be rendered and transforms one of its columns into a hyperlink for
##'   further reference. Refer to the \code{\link{ncbi.entrez.link}} function
##'   as an example
renderFeatureStatsDataTable <- function(x, features=NULL, digits=3,
                                        columns=NULL, feature.link.fn=NULL,
                                        order.by='logFC',
                                        order.dir=c('desc', 'asc'),
                                        filter='none',
                                        length.opts=c(10, 25, 50, 100, 250)) {
  if (is(x, 'MultiGSEAResult')) {
    x <- copy(logFC(x, .external=FALSE))
  }
  stopifnot(is(x, 'data.frame'), !is.null(x$featureId))
  setDT(x)
  order.dir <- match.arg(order.dir)
  if (is.character(features)) {
    x <- subset(x, featureId %in% features)
  }

  ## Figure out what columns to keep in the outgoing datatable
  if (missing(columns)) {
    columns <- c('symbol', 'featureId', 'logFC', 'pval', 'padj')
  } else if (is.null(columns)) {
    columns <- colnames(x)
  }
  if (!is.null(order.by)) {
    stopifnot(is.character(order.by),
              length(order.by) == 1L,
              !is.null(x[[order.by]]))
    columns <- c(columns, order.by)
  }
  bad.cols <- setdiff(columns, colnames(x))
  if (length(bad.cols)) {
    warning("The following columns not found: ", paste(bad.cols, collapse=','))
  }
  columns <- intersect(columns, colnames(x))
  x <- x[, columns, with=FALSE]

  if (is.function(feature.link.fn)) {
    x <- feature.link.fn(x)
  }
  if (!is.null(order.by)) {
    x <- setorderv(x, order.by, order=if (order.dir == 'asc') 1L else -1L)
  }

  ## Tweak length.opts
  if (nrow(x) <= 10) {
    length.opts <- nrow(x)
  } else {
    length.opts <- length.opts[length.opts <= nrow(x)]
    if (tail(length.opts, 1) > nrow(x)) {
      length.opts <- c(head(length.opts, -1L), nrow(x))
    }
  }

  dt.opts <- list(
    pageLength=length.opts[1L],
    lengthMenu=length.opts,
    dom='ltipr')
  out <- DT::datatable(setDF(x), selection='none', escape=FALSE, rownames=FALSE,
                       options=dt.opts, filter=filter)
  roundDT(out)
}

##' Creates an HTML-ized version of \code{tabuleResults}
##'
##' The table produced here is broken into two sections (left and right). The
##' left provides meta information about the geneset collections tested, ie.
##' their names and number of genesets the contain. The right contains columns
##' of results
##'
##' @export
##' @importFrom shiny tags
##' @inheritParams tabulateResults
##' @return a \code{tagList} version of an HTML table for use in a shiny app
summaryHTMLTable.multiGSEA <- function(x, names=resultNames(x), max.p, p.col) {
  stopifnot(is(x, 'MultiGSEAResult'))
  s <- tabulateResults(x, names, max.p, p.col)

  ## The header of this table is two rows
  thead <- tags$thead(
    ## Super headers
    tags$tr(class='super-header',
            tags$th("Gene Sets", colspan="2"),
            tags$th("Analysis Summary", colspan=length(names))),
    ## Sub headers
    do.call(
      tags$tr,
      c(list(class='sub-header',
             tags$th("Collection", class="multiGSEA-summary-meta"),
             tags$th("Count", class="multiGSEA-summary-meta")),
        lapply(names, function(name) {
          target.dom.id <- paste0('multiGSEA-result-', name)
          ## I wanted the headers that had the method names to link to the
          ## tab with the results, but that takes a bit more tweaking
          ## tags$th(tags$a(href=paste0('#', target.dom.id), name))
          tags$th(name)
        })
      )))

  ## The body of the table one row per collection
  collections <- unique(s$collection)
  tbody.rows <- lapply(collections, function(col) {
    sc <- subset(s, collection == col)
    stopifnot(all(names %in% sc$method))
    stopifnot(length(unique(sc$geneset_count)) == 1L)

    mres <- lapply(names, function(name) {
      with(subset(sc, method == name), {
        tags$td(sprintf("%d (%d up; %d down)", sig_count, sig_up, sig_down))
      })
    })

    do.call(tags$tr,
            c(list(tags$td(sc$collection[1L], class='multiGSEA-summary-meta'),
                   tags$td(sc$geneset_count[1L], class='multiGSEA-summary-meta')),
              mres))
  })

  tbody <- tags$tbody(tbody.rows)
  html <- tags$div(class='multiGSEA-summary-table', tags$table(thead, tbody))
}

##' Round the numeric columns of a DT
##'
##' @export
##' @importFrom DT formatRound datatable
##' @param x a DT::datatable
##' @param digits the number of digits to round. If \code{NA}, then no rounding
##'   is performed
##' @return a rounded DT::datatable
##' @examples
##' \dontrun{
##' df <- data.frame(a=rnorm(10), b=sample(letters, 10), c=rnorm(10))
##' roundDT(datatable(df),  digits=2)
##' }
roundDT <- function(x, digits=3) {
  stopifnot(is(x, "datatables"))
  if (is.na(digits)) {
    return(x)
  }
  round.me <- sapply(x$x$data, function(x) {
    is.numeric(x) && any(as.integer(x) != x)
  })
  if (any(round.me)) {
    x <- formatRound(x, round.me, digits=digits)
  }
  x
}
