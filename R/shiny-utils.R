##' Interactively explore a MultiGSEAResult via a shiny app
##'
##' I hope you're sitting down
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
  runApp(system.file('shiny', package='multiGSEA'))
}

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
  out <- result(mg, method) %>%
    dplyr::filter(padj.by.collection <= fdr)
  if (nrow(out)) {
    colls <- sort(unique(out$collection))
    priority <- intersect(prioritize, colls)
    lvls <- c(priority, setdiff(colls, priority))
    out <- out %>%
      mutate(collection=factor(collection, lvls, ordered=TRUE)) %>%
      arrange_(~ - mean.logFC.trim)
  }
  out
}

##' Prepares the datatable of GSEA result statistics for presentation
##'
##' @export
##' @importFrom DT datatable
##' @param x The set of GSEA statistics generated from from
##'   \code{\link{constructGseaResultTable}}
##' @param mg The \code{MultiGSEAResult} object
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

##' Transforms a feature table to have the symbols link to NCBI webpage
##'
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

genecards.entrez.link <- function(x, link.col='symbol') {
  stop("TODO: create genecards link")
}

##' @export
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
    x <- filter(x, featureId %in% features)
  }

  if (is.null(columns)) {
    columns <- union(c('symbol', 'featureId', 'logFC', 'pval', 'padj'), order.by)
  }
  columns <- intersect(columns, colnames(x))
  x <- x[, columns, with=FALSE]

  if (is.function(feature.link.fn)) {
    x <- feature.link.fn(x)
  }
  if (is.character(order.by) && !is.null(x[[order.by]])) {
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

##' @export
##' @importFrom shiny tags
summaryHTMLTable.multiGSEA <- function(x, names=resultNames(x),
                                       max.p, p.col) {
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
##' @importFrom DT formatRound
##' @param x a DT::datatable
##' @param digits the number of digits to round. If \code{NA}, then no rounding
##'   is performed
##' @return a rounded DT::datatable
##' @examples
##' dt <- tibble(a=rnorm(10), b=sample(letters, 10), c=rnorm(10))
##' datatable(dt) %>% roundDT(digits=2)
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
