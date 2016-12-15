setMethod("updateObject", "MultiGSEAResult",
function(object, ..., verbose=FALSE) {
  ## Unfortunately we've internally generated GeneSetDb objects that weren't
  ## entirely properly keyed.
  object@gsd <- updateObject(object@gsd)

  proto <- new("MultiGSEAResult")
  ekeys <- key(proto@logFC)
  xkeys <- key(object@logFC)
  if (!setequal(ekeys, xkeys) || !all.equal(ekeys, xkeys)) {
    setkeyv(object@logFC, ekeys)
  }
  object
})


##' Fetches the GeneSetDb from MultiGSEAResult
##'
##' @export
##'
##' @rdname MultiGSEAResult-utilities
##' @param x \code{MultiGSEAResult}
##' @return The \code{GeneSetDb}
##'
##' @examples
##' vm <- exampleExpressionSet(do.voom=TRUE)
##' gdb <- exampleGeneSetDb()
##' mg <- multiGSEA(gdb, vm, vm$design, 'tumor', methods=NULL)
##' geneSetDb(mg)
geneSetDb <- function(x) {
  stopifnot(is(x, 'MultiGSEAResult'))
  updateObject(x@gsd)
}

setMethod("geneSet", c(x="MultiGSEAResult"),
function(x, i, j, active.only=TRUE, fetch.all=FALSE,
         with.feature.map=FALSE, ..., .external=TRUE) {
  if (!isTRUE(active.only)) {
    warning("active.only ignored on geneSet,MultiGSEAResult")
  }
  if (isTRUE(fetch.all)) {
    warning("fetch.all must be `FALSE` when called on MultiGSEAResult")
  }
  x <- updateObject(x)
  gdb <- geneSetDb(x)
  gs <- geneSet(gdb, i, j, active.only=active.only, fetch.all=fetch.all,
                with.feature.map=with.feature.map, ..., .external=FALSE)
  lfc <- logFC(x, .external=FALSE)[J(gs$featureId)]
  keep <- setdiff(colnames(lfc), colnames(gs))
  out <- cbind(gs, lfc)
  out
})

setMethod("geneSets", c(x="MultiGSEAResult"),
function(x, ..., .external=TRUE) {
  ret.df(geneSets(geneSetDb(x), active.only=TRUE, .external=FALSE),
         .external=.external)
})

## Here's your chance to vectorize this. Once that's done push this code down
## into geneSetUrl,GeneSetDb
setMethod("geneSetURL", c(x="MultiGSEAResult"), function(x, i, j, ...) {
  geneSetURL(geneSetDb(x), i, j, ...)
})

setMethod("geneSetCollectionURLfunction", "MultiGSEAResult",
function(x, i, ...) {
  geneSetCollectionURLfunction(geneSetDb(x), i, ...)
})

setMethod("featureIds", c(x="MultiGSEAResult"),
function(x, i, j, value=c('x.id', 'featureId'), ...) {
  value <- match.arg(value)
  featureIds(geneSetDb(x), i, j, value=value, fetch.all=FALSE,
             active.only=TRUE)
})

##' Fetch the logFC statistics for a geneset
##'
##' @export
##'
##' @param x \code{MultiGSEAResult}
##' @param i The collection name
##' @param j The gene set name
##' @return data.table with the stats for the features in the geneset
geneSetFeatureStats <- function(x, i, j, .external=TRUE) {
  stopifnot(is(x, 'MultiGSEAResult'))
  fids <- featureIds(x@gsd, i, j)
  ret.df(subset(logFC(x, .external=FALSE), featureId %in% fids),
         .external=.external)
}

##' Summarizes useful statistics per gene set from a MultiGSEAResult
##'
##' This function calculates the number of genes that move up/down for the
##' given contrasts, as well as mean and trimmed mean of the logFC and
##' t-statistics. Note that the statistics calculated and returned here are
##' purely a function of the statistics generated at the gene-level stage
##' of the analysis.
##'
##' @export
##'
##' @param x A \code{multiGSEAResult} object
##' @param feature.min.logFC used with \code{feature.max.padj} to identify
##'   the individual features that are to be considered differentially
##'   expressed.
##' @param feature.max.padj used with \code{feature.min.logFC} to identify
##'   the individual features that are to be considered differentially
##'   expressed.
##' @param trim The amount to trim when calculated trimmed \code{t} and
##'   \code{logFC} statistics for each geneset.
##'
##' @return A data.table with statistics on the effect size shift for the gene
##'   set as well as numbers of members that shift up or down.
##'
##' @examples
##'
##' vm <- exampleExpressionSet(do.voom=TRUE)
##' gdb <- exampleGeneSetDb()
##' mg <- multiGSEA(gdb, vm, vm$design, 'tumor', methods=NULL)
##' head(geneSetsStats(mg))
geneSetsStats <- function(x, feature.min.logFC=1,
                          feature.max.padj=0.10, trim=0.10, .external=TRUE) {
  stopifnot(is(x, 'MultiGSEAResult'))
  lfc <- logFC(x, .external=FALSE)

  annotate.lfc <- !missing(feature.min.logFC) ||
    !missing(feature.max.padj) ||
    !all(c('significant', 'direction') %in% names(lfc))
  if (annotate.lfc) {
    lfc <- within(lfc, {
      significant <- abs(logFC) >= feature.min.logFC &
        !is.na(padj) &
        padj <= feature.max.padj
      direction <- ifelse(logFC > 0, 'up', 'down')
    })
  }

  gs <- geneSets(x, .external=FALSE)
  do.by <- key(gs)

  out <- gs[, {
    fids <- featureIds(x, .BY[[1L]], .BY[[2L]])
    stats <- lfc[fids]
    up <- stats$direction == 'up'
    down <- !up
    is.sig <- stats$significant
    t.nona <- stats$t[!is.na(stats$t)]
    list(n.sig=sum(is.sig), n.neutral=sum(!is.sig),
         n.up=sum(up), n.down=sum(down),
         n.sig.up=sum(up & is.sig), n.sig.down=(sum(down & is.sig)),
         mean.logFC=mean(stats$logFC, na.rm=TRUE),
         mean.logFC.trim=mean(stats$logFC, na.rm=TRUE, trim=trim),
         mean.t=mean(stats$t, na.rm=TRUE),
         mean.t.trim=mean(stats$t, na.rm=TRUE, trim=trim))
  }, by=do.by]
  setkeyv(out, do.by)
  ret.df(out, .external=.external)
}

##' Extract the individual fold changes statistics for elements in the
##' expression object.
##'
##' @export
##' @param x A \code{MultiGSEAResult}
##' @return The log fold change \code{data.table}
##'
##' @examples
##' vm <- exampleExpressionSet(do.voom=TRUE)
##' gdb <- exampleGeneSetDb()
##' mg <- multiGSEA(gdb, vm, vm$design, 'tumor', methods=NULL)
##' lfc <- logFC(mg)
logFC <- function(x, .external=TRUE) {
  stopifnot(is(x, 'MultiGSEAResult'))
  ret.df(x@logFC, .external=.external)
}

##' Fetch names of GSEA methods that were run from a \code{MultiGSEAResult}
##'
##' @export
##' @rdname results
##' @param x A MultiGSEAResult
resultNames <- function(x) {
  stopifnot(is(x, 'MultiGSEAResult'))
  names(x@results)
}

##' Returns method names that were not run on a MultiGSEAResult
invalidMethods <- function(x, names, as.error=FALSE) {
  stopifnot(is(x, 'MultiGSEAResult'))
  stopifnot(is.character(names))
  if (length(names) == 0 && FALSE) {
    warning("No mehod `names` passed to invalidMethods", immediate.=TRUE)
  }
  bad.names <- setdiff(names, resultNames(x))
  if (length(bad.names) && as.error) {
    stop("Illegal result names queried: ", paste(bad.names, collapse=','))
  }
  bad.names
}


##' Extract the result from a given enrichment test from a MultiGSEAResult
##'
##' TODO: Enable a way for caller to get a subset of the genesets tested,
##' perhaps this will happen at the "collectoin" level. In this case, it's
##' not clear if the pval correction should be redone to only account for
##' the genesets that are returned from.
##'
##' @export
##' @rdname results
##'
##' @param x MultiGSEAResult
##' @param name the names of the results desired
##' @param stats.only logical, set to \code{FALSE} if you want to return all
##'   (column-wise) data for each result. By default only the pvalues,
##'   adjusted pvalues, and rank are returned.
##' @param rank.by the statistic to use to append a \code{rank} column for the
##'   geneset result. By default we rank by pvalue calculated by the GSEA
##'   method. You can rank the results based on the trimmed mean of the logFC's
##'   calculated for all of the features in the geneset (\code{"logFC"}), the
##'   trimmed t-statistics of the these features (\code{"t"}), or the "J-G"
##'   statistic of the geneset.
##' @param add.suffix If \code{TRUE}, adds \code{.name} as a suffix to the
##'   columns of the \code{method}-specific statistics returned, ie. the
##'   \code{pval} column from the \code{camera} result will be turned to
##'   \code{pval.camera}.
##'
##' @return a data.table with the results from the requested method.
result <- function(x, name, stats.only=FALSE,
                   rank.by=c('pval', 't', 'logFC'),
                   add.suffix=FALSE, .external=TRUE) {
  stopifnot(is(x, 'MultiGSEAResult'))
  if (is.null(resultNames(x)) || length(resultNames(x)) == 0) {
    if (missing(name)) name <- NULL
    return(results(x, name, .external=.external))
  }
  if (length(resultNames(x)) == 1L) {
    name <- resultNames(x)
  }
  ## no methods run?
  stopifnot(isSingleCharacter(name))
  invalidMethods(x, name, as.error=TRUE)
  stopifnot(isSingleLogical(stats.only))
  rank.by <- match.arg(rank.by)
  stopifnot(isSingleLogical(add.suffix))

  out <- copy(geneSets(x, .external=FALSE))

  pval.col <- 'pval'
  rank.col <- 'rank'
  if (add.suffix) {
    pval.col <- paste0('pval.', name)
    rank.col <- paste0('rank.', name)
  }

  res <- local({
    r <- x@results[[name]]
    if (!isTRUE(all.equal(out[, key(out), with=FALSE], r[, key(out), with=FALSE])[1])) {
      stop("Unexpected geneset ordering in `", name, "` result")
    }
    if (stats.only) {
      res.cols <- c('pval', 'padj', 'padj.by.collection',
                    ## roast has these, presumably other methods will have
                    ## others
                    'pval.mixed', 'padj.mixed')
      ## Select any column that starts with pval or padj
      res.cols <- names(r)[grepl('^(pval\\.?|padj\\.?)', names(r))]
    } else {
      res.cols <- setdiff(names(r), names(out))
    }
    r <- r[, res.cols, with=FALSE]
    if (add.suffix) {
      setnames(r, res.cols, paste(res.cols, name, sep='.'))
    }
    r
  })

  for (col in names(res)) {
    out[, (col) := res[[col]]]
  }

  missing.pvals <- is.na(out[[pval.col]])
  n.missing <- sum(missing.pvals)
  if (any(missing.pvals)) {
    msg <- sprintf("%d missing pvalues for actige genesets from %s",
                   n.missing, name)
    warning(msg, immediate.=TRUE)
  }

  ranks <- switch(rank.by,
                  logFC=rank(-abs(out$mean.logFC.trim), ties.method="min"),
                  t=rank(-abs(out$mean.t.trim), ties.method="min"),
                  pval=rank(out[[pval.col]], ties.method="min"))
  out[, (rank.col) := ranks]

  ret.df(out, .external=.external)
}

##' A summary of one/some/all of the results that were run.
##'
##' If the result from more than one method is requested, the column names for
##' these results will be suffixed with the method name.
##'
##' @export
##' @rdname results
##'
##' @param x \code{MultiGSEAResult}
##' @param names The results you want to cbind together for the output
##' @param stats.only logical, set to \code{FALSE} if you want to return all
##'   (column-wise) data for each result. By default only the pvalues,
##'   adjusted pvalues, and rank are returned.
##' @param add.suffix If \code{TRUE}, adds \code{.name} as a suffix to the
##'   columns of the \code{method}-specific statistics returned. This is forced
##'   to \code{TRUE} if \code{length(names) > 1L}
##'
##' @return a data.table with all of the results. The results will optionally
##'   be suffixed with the method name that generated them when
##'   \code{length(names) > 1L}
results <- function(x, names=resultNames(x), stats.only=TRUE,
                    rank.by=c('pval', 'logFC', 't'),
                    add.suffix=length(names) > 1L, .external=TRUE) {
  stopifnot(is(x, 'MultiGSEAResult'))
  if (is.null(resultNames(x)) || length(resultNames(x)) == 0L) {
    ## No methods were run, you can only return geneset stats
    if (!is.null(names) || length(names) > 0L) {
      ## User is asking for something that is not there
      stop("No GSEA methods were run, you can only get set statistics")
    }
    message("No GSEA methods were run, only geneset statistics have ",
            "been returned.")
    return(ret.df(geneSets(x, .external=FALSE)))
  }
  invalidMethods(x, names)
  stopifnot(isSingleLogical(stats.only))
  rank.by <- match.arg(rank.by)
  if (length(names) > 1L && add.suffix == FALSE) {
    add.suffix <- TRUE
    warning("Forcing method suffix to generated result columns",
            immediate.=TRUE)
  }

  ## Idiomatic data.table, non-idiomatic R
  out <- copy(geneSets(x, .external=FALSE))
  for (name in names) {
    res <- result(x, name, stats.only, rank.by, add.suffix, .external=FALSE)
    for (col in setdiff(names(res), names(out))) {
      out[, (col) := res[[col]]]
    }
  }

  ret.df(out, .external=.external)
}

##' Create a summary table to indicate number of significant genesets per
##' collection, per method
##'
##' @export
##' @rdname results
##'
##' @param x \code{MultiGSEAResult}
##' @param names The results you want to cbind together for the output
##' @param max.p The maximum padj value to consider a result significant
##' @param p.col use padj or padj.by.collection?
##'
##' @return a data.table that summarizes the significant results per method
##'   per collection for the GSEA that was run
tabulateResults <- function(x, names=resultNames(x), max.p=0.30,
                            p.col=c('padj', 'padj.by.collection', 'pval')) {
  stopifnot(is(x, 'MultiGSEAResult'))
  invalidMethods(x, names)
  stopifnot(isSingleNumeric(max.p))
  p.col <- match.arg(p.col)
  res <- lapply(names, function(wut) {
    r <- result(x, wut, .external=FALSE)
    ## some results (like goseq) don't have just "padj" or "pval" columns,
    ## because it has pval.over and pval.under, so let's just grab the first
    ## pval or padj "hit"
    idx <- grep(p.col, names(r))[1L]
    if (is.na(idx)) {
      r$pcol <- rep(NA, nrow(r))
    } else {
      r$pcol <- r[[idx]]
    }
    ## r$pcol <- r[[p.col]]
    r[, {
      list(method=wut, geneset_count=length(pcol),
           sig_count=sum(pcol <= max.p, na.rm=TRUE),
           sig_up=sum(pcol <= max.p & mean.logFC.trim > 0, na.rm=TRUE),
           sig_down=sum(pcol <= max.p & mean.logFC.trim < 0, na.rm=TRUE))
    }, by='collection']
  })
  ret.df(rbindlist(res))
}

##' Subset MultiGSEAResult to include include results for specified genesets
##'
##' This isn't exported yet because I don't like its implementation
##' @param x MultiGSEAResult
##' @param keep logical vector as long as there are numbers of results
##' @return a \code{MultiGSEAResult} that has only the results for the specified
##'   genesets.
subset.MultiGSEAResult <- function(x, keep) {
  stopifnot(is(x, 'MultiGSEAResult'))
  did.gsea <- length(x@results) > 0
  ## length of x@results is 0 if no methods were run (only stats calc'd)
  nr <- nrow(if (did.gsea) x@results[[1]] else geneSets(x, .external=FALSE))
  if (!is.logical(keep) && lenght(keep) != nr) {
    stop("The `keep` vector is FUBAR'd")
  }
  if (did.gsea) {
    new.res <- lapply(x@results, function(x) subset(x, keep))
    x@results <- new.res
    gsets <- with(new.res[[1]], paste(collection, name, sep=':'))
    x@gsd@table <- x@gsd@table[paste(collection, name, sep=':') %in% gsets]
  } else {
    x@gsd@table <- x@gsd@table[keep]
  }

  x
}

setMethod("show", "MultiGSEAResult", function(object) {
  msg <- paste("multiGSEA result (max FDR by collection set to 30%)",
               "---------------------------------------------------", sep='\n')
  cat(msg, "\n")
  if (length(resultNames(object)) == 0) {
    cat("No GSEA methods were run, only geneset level statistics calculated")
  } else {
    base::print.data.frame(as.data.frame(tabulateResults(object, max.p=0.30)))
  }
  cat("\n")
})

##' Assembles a matrix of nominal or adjusted pvalues from a multiGSEA result
##'
##' @export
##' @param x The output from multiGSEA
##' @param pval Are we testing pvalues or adjusted pvalues?
##' @return A matrix of the desired pvalues for all genesets
p.matrix <- function(x, names=resultNames(x),
                     pcol=c('padj', 'padj.by.collection', 'pval')) {
  stopifnot(is(x, 'MultiGSEAResult'))
  invalidMethods(x, names)
  pcol <- match.arg(pcol)
  res <- results(x, names, add.suffix=TRUE)
  regex <- sprintf('^%s\\.', pcol)
  col.idx <- grep(regex, names(res))
  if (pcol != 'padj.by.collection') {
    col.idx <- setdiff(col.idx, grep('padj.by.collection', names(res)))
  }
  out <- as.matrix(res[, col.idx, with=FALSE])
  colnames(out) <- sub(regex, '', colnames(out))
  out
}
