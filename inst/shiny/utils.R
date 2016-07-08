## This code should be pushed back into the package
round.dt <- function(x, digits=2) {
  stopifnot(is(x, 'data.table'))
  for (cname in names(x)[sapply(x, is.numeric)]) {
    vals <- x[[cname]]
    if (!is.integer(vals)) {
      x[, (cname) := round(vals, digits=2)]
    }
  }
  x
}

render.dt <- function(x, mg) {
  if (FALSE) {
    mg <- xmg
    xresult <- 'camera'
    cnames <- lapply(resultNames(mg), function(n) colnames(result(x, n)))
    Reduce(intersect, cnames[-1], cnames[[1]])
    x <- result(x, 'camera')
  }
  stopifnot(is(x, 'data.table'))

  rcols <- c('collection'='collection', 'name'='name', 'n'='n', 'pval'='pval',
             ## 'padj'='padj', 'padj.by.collection'='padjByColl',
             'padj.by.collection'='FDR',
             'mean.logFC.trim'='logFC', 'n.sig.up'='nSigUp',
             'n.sig.down'='nSigDown', 'n.up'='nUp', 'n.down'='nDown')

  res <- x[, names(rcols), with=FALSE]
  setnames(res, names(rcols), rcols)
  res <- round.dt(res)

  res[, name := {
    url <- geneSetURL(mg, as.character(collection), name)
    xname <- gsub('_', ' ', name)
    html <- '<a href="%s" target="_blank">%s</a>'
    ifelse(is.na(url), xname, sprintf(html, url, xname))
  }]

  ## If MSigDB hallmark collections are include, always put those first
  ## regardless of padj sorting
  dt.order <- list()
  # if ('h' %in% res$collection) {
  #   lvls <- c('h', sort(setdiff(res$collection, 'h')))
  #   ## dt.order <- list(list(0, 'asc'))
  # } else {
  #   lvls <- sort(unique(res$collection))
  # }
  # res[, collection := factor(collection, lvls, ordered=TRUE)]

  res <- res[, rcols, with=FALSE]
  lfc.col <- which(colnames(res) == 'logFC') - 1L
  dt.order[[length(dt.order) + 1]] <- list(lfc.col, 'desc')

  length.opts <- c(10, 20, 50, 100, 250)
  length.opts <- length.opts[length.opts < nrow(res)]
  length.opts <- c(length.opts, nrow(res))

  # if ('h'%in% res$collection) {
  #   res <- arrange(res, collection, -logFC)
  # } else {
  #   res <- arrange(res, -logFC)
  # }

  dt <- datatable(
    res, filter='top',
    selection=list(mode='single', selected=1, target='row'),
    extensions='Buttons',
    escape=FALSE, rownames=FALSE,
    options=list(
      dom='lBtpir',
      # order=dt.order,
      pageLength=length.opts[1L],
      lengthMenu=length.opts,
      buttons=c('copy', 'csv', 'excel')))

  DT::renderDataTable(dt)
}

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
