##' Builds the subset of the GSEA statistics to present to the user
##'
##' This function sets the collection to a factor and puts the MSigDB Hallmark
##' collection ("h") at the top of the totem pole, if it is included in
##' \code{mg}.
##'
##' @export
##' @param mg \code{MultiGSEAResult} object
##' @param method the method to show statistics for
##' @param the FDR cut off to present statistics for
##' @return a data.table of the statistics that match the filtering criteria.
##'   A 0-row data.table is returned if nothing passes.
constructGseaResultTable <- function(mg, method, fdr) {
  out <- result(mg, method) %>%
    dplyr::filter(padj.by.collection <= fdr)

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
  }

  out
}

##' Prepares the datatable arguments to use for the gsea.result.table display
##' @export
##' @param x The \code{gsea.result.table} \code{data.table} object
##' @param mg The \code{MultiGSEAResult} object
##' @return a DT::DataTable
renderGseaResultTableDataTable <- function(x, method, mg, digits=3) {
  stopifnot(requireNamespace("DT"))
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

  ## If MSigDB hallmark collections are include, always put those first
  ## regardless of padj sorting
  dt.order <- list()
  lfc.col <- which(colnames(res) == 'logFC') - 1L
  dt.order[[length(dt.order) + 1]] <- list(lfc.col, 'desc')

  length.opts <- c(10, 20, 50, 100, 250)
  length.opts <- length.opts[length.opts < nrow(res)]
  length.opts <- c(length.opts, nrow(res))

  ## Download the GSEA statistics is now handled through a shiny::downloadButton
  ## so that we can display this table using `server=TRUE` to make it more
  ## performant.
  ## Configure options for Button so that the downloaded file has a propper
  ## filename
  ## for button options cf: https://github.com/rstudio/DT/issues/213
  ## bfn <- paste0('multiGSEA-', method)
  ## btn.opts <- list(
  ##   # 'colvis',
  ##   list(
  ##     extend='collection',
  ##     buttons=list(
  ##       list(extend='copy'),
  ##       list(extend='csv', filename=bfn),
  ##       list(extend='excel', filename=bfn)),
  ##     text='Export'))
  dt.opts <- list(
    # dom='lBtpir',
    dom='ltpir',
    # order=dt.order,
    scrollX=TRUE,
    pageLength=length.opts[1L],
    # buttons=btn.opts,
    lengthMenu=length.opts)
  dtargs <- list(data=res, filter='top',
                 selection=list(mode='single', selected=NA, target='row'),
                 # extensions='Buttons',
                 escape=FALSE, rownames=FALSE,
                 options=dt.opts)
  do.call(DT::datatable, dtargs) %>% roundDT(digits=digits)
}

##' Prepares the datatable to show the geneset members and their logFC stats
##'
##' @export
##' @return a list of arguments that you can \code{do.call, datatable}
renderGeneSetStatsDataTable <- function(gstats, name, digits=3) {
  requireNamespace('DT')
  gcols <- c('symbol', 'featureId', 'logFC', 'padj')
  gs <- gstats[, gcols, with=FALSE]
  setnames(gs, 'padj', 'FDR')
  if (nrow(gs) < 10) {
    length.opts <- nrow(gs)
  } else {
    length.opts <- c(6, 15, 50, 100)
    length.opts <- length.opts[length.opts < nrow(gs)]
    length.opts <- c(length.opts, nrow(gs))
  }

  gs.name <- gsub('[^A-Za-z0-9]', '_', name)
  bfn <- sprintf('%s-multiGSEA-gene-statistics', gs.name)
  btn.opts <- list(
    # 'colvis',
    list(
      extend='collection',
      buttons=list(
        list(extend='copy'),
        list(extend='csv', filename=bfn),
        list(extend='excel', filename=bfn)),
      text='Export'))

  dt.opts <- list(
    pageLength=length.opts[1L],
    lengthMenu=length.opts,
    # buttons=btn.opts,
    # dom='lBtipr',
    dom='ltipr')

  dtargs <- list(data=gs, selection='none', extension='Buttons', escape=FALSE,
                 rownames=FALSE, options=dt.opts)
  do.call(DT::datatable, dtargs) %>% roundDT(digits=digits)
}

##' @export
summaryHTMLTable.multiGSEA <- function(x, names=resultNames(x),
                                       max.p, p.col) {
  requireNamespace('shiny')
  stopifnot(is(x, 'MultiGSEAResult'))
  s <- tabulateResults(x, names, max.p, p.col)

  ## The header of this table is two rows
  thead <- shiny::tags$thead(
    ## Super headers
    shiny::tags$tr(class='super-header',
                   shiny::tags$th("Gene Sets", colspan="2"),
                   shiny::tags$th("Analysis Summary", colspan=length(names))),
    ## Sub headers
    do.call(
      shiny::tags$tr,
      c(list(class='sub-header',
             shiny::tags$th("Collection", class="multiGSEA-summary-meta"),
             shiny::tags$th("Count", class="multiGSEA-summary-meta")),
        lapply(names, function(name) {
          target.dom.id <- paste0('multiGSEA-result-', name)
          ## I wanted the headers that had the method names to link to the
          ## tab with the results, but that takes a bit more tweaking
          ## tags$th(tags$a(href=paste0('#', target.dom.id), name))
          shiny::tags$th(name)
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
        shiny::tags$td(sprintf("%d (%d up; %d down)", sig_count, sig_up, sig_down))
      })
    })

    do.call(shiny::tags$tr,
            c(list(shiny::tags$td(sc$collection[1L], class='multiGSEA-summary-meta'),
                   shiny::tags$td(sc$geneset_count[1L], class='multiGSEA-summary-meta')),
              mres))
  })

  tbody <- shiny::tags$tbody(tbody.rows)
  html <- shiny::tags$div(class='multiGSEA-summary-table', shiny::tags$table(thead, tbody))
}

##' Round the numeric columns of a DT
##'
##' @export
##' @param a DT::datatable
##' @param digits the number of digits to round. If \code{NA}, then no rounding
##'   is performed
##' @return a rounded DT::datatable
roundDT <- function(x, digits=3) {
  requireNamespace('DT')
  stopifnot(is(x, "datatables"))
  if (is.na(digits)) {
    return(x)
  }
  round.me <- sapply(x$x$data, function(x) {
    is.numeric(x) && any(as.integer(x) != x)
  })
  if (any(round.me)) {
    x <- DT::formatRound(x, round.me, digits=digits)
  }
  x
}
