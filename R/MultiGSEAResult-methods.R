#' Makes MultiGSEAResult objects compatible with old GNE multiGSEA Explorer App
#'
#' This is not documented or exported on purpose and will be removed in the near
#' future without warning, so: don't use it.
#' @noRd
downgradeObject <- function(x) {
  stopifnot(is(x, 'MultiGSEAResult'))
  x@results <- sapply(names(x@results), function(res) {
    out <- result(x, res, as.dt=TRUE)
    setattr(out, 'rawresult', FALSE)
  }, simplify=FALSE)
  x
}

#' Combines two MultiGSEAResult objects together.
#'
#' This would be useful when you want to add a GSEA result to an already
#' existing one. `append` would be more appropriate, but ...
#'
#' When would you want to do that? Imagine a shiny app that drives multiGSEA.
#' You might want to present the results of each analysis as they come "online",
#' so you would run them independently and make them available to the user
#' immediately after they each finish (ie. in combination with the promises
#' package).
#'
#' @param x A `MultiGSEAResult` object
#' @param y A `MultiGSEAResult` object
#' @param rename.x A named vector that used to match resultNames(x) and remane
#'   them to something different. `names(rename.x)` should match whatever you
#'   want to change in `resultNames(x)`, and the values are the new names of
#'   the result.
#' @param rename.y Same as `rename.x`, but for the results in `y`.
#' @param ... more things
#'
#' @importMethodsFrom BiocGenerics combine
#' @exportMethod combine
#' @return A combined `MultiGSEAResult` object
#' @examples
#' mg1 <- exampleMultiGSEAResult()
#' mg2 <- exampleMultiGSEAResult()
#' mgc <- combine(mg1, mg2)
setMethod("combine", c(x = "MultiGSEAResult", y = "MultiGSEAResult"),
function(x, y, rename.x = NULL, rename.y = NULL, ...) {
  if (!isTRUE(all.equal(geneSetDb(x), geneSetDb(y)))) {
    stop("Can't combine MultiGSEAResults with different internal GeneSetDb's")
  }
  if (!isTRUE(all.equal(logFC(x), logFC(y)))) {
    stop("Can't combine MultiGSEAResults with different internal logFC stats")
  }

  xnames <- .replace(resultNames(x), rename.x)
  ynames <- .replace(resultNames(y), rename.y)
  rnames <- make.unique(c(xnames, ynames))

  out <- x
  out@results <- setNames(c(x@results, y@results), rnames)
  out
})

#' Fetches the GeneSetDb from MultiGSEAResult
#'
#' @export
#'
#' @rdname MultiGSEAResult-utilities
#' @param x \code{MultiGSEAResult}
#' @return The \code{GeneSetDb}
#'
#' @examples
#' vm <- exampleExpressionSet(do.voom=TRUE)
#' gdb <- exampleGeneSetDb()
#' mg <- multiGSEA(gdb, vm, vm$design, 'tumor', methods=NULL)
#' geneSetDb(mg)
geneSetDb <- function(x) {
  stopifnot(is(x, 'MultiGSEAResult'))
  x@gsd
}

#' @rdname geneSet
setMethod("geneSet", c(x="MultiGSEAResult"),
function(x, i, j, active.only=TRUE, with.feature.map=FALSE, ...,
         collection = NULL, name = NULL, as.dt=FALSE) {
  if (!isTRUE(active.only)) {
    warning("active.only set to TRUE for geneSet,MultiGSEAResult")
    active.only <- TRUE
  }
  gdb <- geneSetDb(x)
  gs <- geneSet(gdb, i, j, active.only=TRUE, with.feature.map=with.feature.map,
                ..., collection = collection, name = name, as.dt=TRUE)
  lfc <- logFC(x, as.dt=TRUE)[list(gs$feature_id), on='feature_id']
  stopifnot(all(lfc$feature_id == gs$feature_id))
  lfc$feature_id <- NULL
  # Columns that appear in the geneset and also appear in the `logFC` data.frame
  # will be renamed to avoid collision. The duplicated columns used to be
  # removed with a preference to keep the columns in `gs`, however we might
  # have a `logFC` column in the geneset database, which would then override
  # the `logFC` column from the differential expression analysis.
  rename <- intersect(colnames(lfc), colnames(gs))
  if (length(rename)) setnames(gs, rename, paste0(rename, ".gs"))
  out <- cbind(gs, lfc)
  if (!as.dt) setDF(out)
  out
})

setMethod("geneSets", c(x="MultiGSEAResult"),
function(x, ..., as.dt=FALSE) {
  geneSets(geneSetDb(x), active.only=TRUE, as.dt=as.dt)
})

# Here's your chance to vectorize this. Once that's done push this code down
# into geneSetUrl,GeneSetDb
setMethod("geneSetURL", c(x="MultiGSEAResult"), function(x, i, j, ...) {
  geneSetURL(geneSetDb(x), i, j, ...)
})

setMethod("geneSetCollectionURLfunction", "MultiGSEAResult",
function(x, i, ...) {
  geneSetCollectionURLfunction(geneSetDb(x), i, ...)
})

#' @rdname featureIds
setMethod("featureIds", c(x="MultiGSEAResult"),
function(x, i, j, value=c('feature_id', 'x.id', 'x.idx'),
         active.only=TRUE, ...) {
  value <- match.arg(value)
  if (!isTRUE(active.only)) {
    warning("The featureIds() accessor for a MultiGSEAResult enforces ",
            "`active.only` to be set to TRUE")
  }
  featureIds(geneSetDb(x), i, j, value=value, active.only=TRUE, ...)
})

#' Summarizes useful statistics per gene set from a MultiGSEAResult
#'
#' This function calculates the number of genes that move up/down for the
#' given contrasts, as well as mean and trimmed mean of the logFC and
#' t-statistics. Note that the statistics calculated and returned here are
#' purely a function of the statistics generated at the gene-level stage
#' of the analysis.
#'
#' @export
#'
#' @param x A `MultiGSEAResult` object
#' @param feature.min.logFC used with `feature.max.padj` to identify
#'   the individual features that are to be considered differentially
#'   expressed.
#' @param feature.max.padj used with `feature.min.logFC` to identify
#'   the individual features that are to be considered differentially
#'   expressed.
#' @param trim The amount to trim when calculated trimmed `t` and
#'   `logFC` statistics for each geneset.
#'
#' @return A data.table with statistics at the gene set level across the
#'   prescribed contrast run on `x`. These statistics are independent
#'   of any particular GSEA method, but rather summarize aggregate shifts
#'   of the gene sets individual features. The columns included in the output
#'   are summarized below:
#'
#' * `n.sig`: The number of individual features whose `abs(logFC)` and padj
#'    thersholds satisfy the criteria of the `feature.min.logFC` and
#'    `feature.max.padj` parameters of the original [multiGSEA()] call
#' * `n.neutral`: The number of individual features whose abs(logFC) and padj
#'    thersholds do not satisfy the `feature.*` criteria named above.
#' * `n.up, n.down`: The number of individual features with `logFC > 0` or
#'   `logFC < 0`, respectively, irrespective of the `feature.*` thresholds
#'    referenced above.
#' * `n.sig.up, n.sig.down`: The number of individual features that pass the
#'   `feature.*` thresholds and have logFC > 0 or logFC < 0, respectively.
#' * `mean.logFC, mean.logFC.trim`: The mean (or trimmed mean) of the individual
#'   logFC estimates for the features in the gene set. The amount of trim is
#'   specified in the `trim` parameter of the [multiGSEA()] call.
#' * `mean.t, mean.t.trim`: The mean (or trimmed mean) of the individual
#'   t-statistics for the features in the gene sets. These are `NA` if the input
#'   expression object was a `DGEList`.
#'
#' @examples
#' vm <- exampleExpressionSet(do.voom=TRUE)
#' gdb <- exampleGeneSetDb()
#' mg <- multiGSEA(gdb, vm, vm$design, 'tumor')
#' head(geneSetsStats(mg))
geneSetsStats <- function(x, feature.min.logFC=1, feature.max.padj=0.10,
                          trim=0.10, reannotate.significance = FALSE,
                          as.dt=FALSE) {
  stopifnot(is(x, 'MultiGSEAResult'))
  lfc <- logFC(x, as.dt=TRUE)

  # reannotate.significance was added to better accomodate annotations that
  # were passed in the xmeta. data.frame through the multiGSEA call. You
  # should revisit this when you revisit how data.frame input support superseds
  # the xmeta. hacks we have in place now.
  annotate.lfc <- !missing(feature.min.logFC) ||
    !missing(feature.max.padj)
  annotate.lfc <- annotate.lfc && reannotate.significance
  annotate.lfc <- annotate.lfc ||
    !all(c('significant', 'direction') %in% names(lfc))


  is.ttest <- "logFC" %in% names(lfc)
  if (annotate.lfc) {
    lfc <- within(lfc, {
      if (is.ttest) {
        significant <- abs(logFC) >= feature.min.logFC &
          !is.na(padj) &
          padj <= feature.max.padj
        direction <- ifelse(logFC > 0, 'up', 'down')
      } else {
        # This is an ANOVA
        significant <- !is.na(padj) & padj <= feature.max.padj
        direction <- rep("ambiguous", length(padj))
      }
    })
  }

  gs <- geneSets(x, as.dt=TRUE)
  do.by <- key(gs)

  out <- gs[, {
    fids <- featureIds(x, .BY[[1L]], .BY[[2L]])
    stats <- lfc[fids, on='feature_id']
    up <- stats$direction == 'up'
    down <- stats$direction == "down"
    is.sig <- stats$significant
    t.nona <- stats$t[!is.na(stats$t)]
    if (is.ttest) {
      mean.logFC <- mean(stats$logFC, na.rm=TRUE)
      mean.logFC.trim <- mean(stats$logFC, na.rm=TRUE, trim=trim)
      mean.t <- mean(stats$t, na.rm=TRUE)
      mean.t.trim <- mean(stats$t, na.rm=TRUE, trim=trim)
    } else {
      mean.logFC <- rep(NA_real_, .N)
      mean.logFC.trim <- rep(NA_real_, .N)
      mean.t <- rep(NA_real_, .N)
      mean.t.trim <- rep(NA_real_, .N)
    }
    list(n.sig=sum(is.sig), n.neutral=sum(!is.sig),
         n.up=sum(up), n.down=sum(down),
         n.sig.up=sum(up & is.sig), n.sig.down=(sum(down & is.sig)),
         mean.logFC=mean.logFC,
         mean.logFC.trim=mean.logFC.trim,
         mean.t=mean.t,
         mean.t.trim=mean.t.trim)
  }, by=do.by]
  setkeyv(out, do.by)

  if (!as.dt) setDF(out)
  out
}

#' Extract the individual fold changes statistics for elements in the
#' expression object.
#'
#' @export
#' @param x A [MultiGSEAResult()]
#' @template asdt-param
#' @return The log fold change `data.table``
#'
#' @examples
#' vm <- exampleExpressionSet(do.voom=TRUE)
#' gdb <- exampleGeneSetDb()
#' mg <- multiGSEA(gdb, vm, vm$design, 'tumor', methods=NULL)
#' lfc <- logFC(mg)
logFC <- function(x, as.dt=FALSE) {
  stopifnot(is(x, 'MultiGSEAResult'))
  out <- x@logFC
  if (!as.dt) out <- setDF(copy(out))
  out
}

#' Helper funtion: returns method names that were not run on a MultiGSEAResult
#'
#' @noRd
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

#' Interrogate the results of a multiGSEA analysis stored in a MultiGSEAResult
#'
#' @description
#' The `resultNames`, `result`, and `results` functions enable
#' you to explore the results of the analysis run with \code{\link{multiGSEA}}.
#'
#' The results that are stored within a `MultiGSEAResult` object have a
#' more or less 1:1 mapping with the values passed as `methods`, parameter
#' of the [multiGSEA()] call.
#'
#' @details
#' The product of an indivdual GSEA is consumed by the corresponding
#' `do.<METHOD>` function and converted into a data.table of results that
#' is internally stored.
#'
#' Use the `resultNames()` function to identify which results are available
#' for interrogation. The `result()` function returns the statistics of
#' one individual result, and the `results()` function combines the results
#' from the specified methods into an arbitrarily wide data.table with
#' `method`-suffixed column names.
#'
#' Use the `tabulateResults()` function to create a summary table that
#' tallies the number of significant genesets per collection, per method at
#' the specified FDR thresholds.
#'
#' @export
#' @rdname results
#'
#' @examples
#' ## Refer to the examples in ?multiGSEA
resultNames <- function(x) {
  stopifnot(is(x, 'MultiGSEAResult'))
  names(x@results)
}

#' @export
#' @rdname results
result <- function(x, ...) {
  UseMethod("result", x)
}

#' @export
#' @rdname results
#' @param x MultiGSEAResult
#' @param name the names of the results desired
#' @param stats.only logical, set to `FALSE` if you want to return all
#'   (column-wise) data for each result. By default only the pvalues,
#'   adjusted pvalues, and rank are returned.
#' @param rank.by the statistic to use to append a `rank` column for the
#'   geneset result. By default we rank by pvalue calculated by the GSEA
#'   method. You can rank the results based on the trimmed mean of the logFC's
#'   calculated for all of the features in the geneset (`"logFC"`), or the
#'   trimmed t-statistics of the these features (`"t"`).
#' @param add.suffix If `TRUE`, adds `.name` as a suffix to the
#'   columns of the `method`-specific statistics returned, ie. the
#'   `pval` column from the `"camera"` result will be turned to
#'   `pval.camera`.
#'
#' @return a data.table with the results from the requested method.
result.MultiGSEAResult <- function(x, name = NULL, stats.only=FALSE,
                                   rank.by=c('pval', 't', 'logFC'),
                                   add.suffix=FALSE, as.dt=FALSE, ...) {
  stopifnot(is(x, 'MultiGSEAResult'))
  if (is.null(resultNames(x)) || length(resultNames(x)) == 0) {
    return(results(x, name, as.dt=as.dt))
  }
  if (is.null(name)) name <- resultNames(x)[1L]
  ## no methods run?
  stopifnot(isSingleCharacter(name))
  invalidMethods(x, name, as.error=TRUE)
  stopifnot(isSingleLogical(stats.only))
  rank.by <- match.arg(rank.by)
  stopifnot(isSingleLogical(add.suffix))

  gs <- geneSets(x, as.dt=TRUE)
  base.colnames <- colnames(gs)
  out <- copy(gs)

  res <- local({
    # The <name> of the GSEA methods should not have a "." in them. If there
    # is a ".", we take this to mean (for now) some modification of the
    # method named up to that dot, for instatnce "goseq.up" is "goseq" method.
    # We can have other "."s appear in the future. For instance, someone might
    # want to run camera twice with different values for inter.gene.cor. In
    # this case, one "name" could be "camera.igc01" or something. This part
    # isn't formal, and will be more formalized (and changed) in a future
    # release.
    fnname <- sub("\\..*$", "", name)
    fnname <- paste0('mgres.', fnname)
    fnfetch <- getFunction(fnname)
    r <- fnfetch(x@results[[name]], geneSetDb(x))
    kosher <- isTRUE(
      all.equal(
        out[, key(out), with = FALSE],
        r[, key(out), with = FALSE],
        check.attributes = FALSE))
    if (!isTRUE(kosher)) {
      stop("Unexpected geneset ordering in `", name, "` result")
    }
    ## All mgres.* function must return a pval and padj column
    if (!all(c('pval', 'padj') %in% colnames(r))) {
      stop(fnname, " does not provide pval and/or padj column")
    }

    if (stats.only) {
      ## Select any column that starts with pval or padj
      res.cols <- names(r)[grepl('^(pval\\.?|padj\\.?)', names(r))]
    } else {
      res.cols <- setdiff(names(r), names(out))
    }
    r <- r[, res.cols, with=FALSE]
    r
  })

  for (col in names(res)) {
    out[, (col) := res[[col]]]
  }

  missing.pvals <- is.na(out[['pval']])
  n.missing <- sum(missing.pvals)
  if (any(missing.pvals)) {
    msg <- sprintf("%d missing pvalues for active genesets from %s",
                   n.missing, name)
    warning(msg, immediate.=TRUE)
  }

  out[, padj.by.collection := p.adjust(pval, 'BH'), by='collection']

  ranks <- switch(rank.by,
                  logFC=rank(-abs(out$mean.logFC.trim), ties.method="min"),
                  t=rank(-abs(out$mean.t.trim), ties.method="min"),
                  pval=rank(out[['pval']], ties.method="min"))
  out[, rank := ranks]

  if (add.suffix) {
    rename.cols <- setdiff(colnames(out), base.colnames)
    setnames(out, rename.cols, paste(rename.cols, name, sep='.'))
  }

  # let's get the pvalue related columns up front
  upfront <- c(
    "collection", "name", "active", "N", "n",
    names(out)[grepl('^(pval\\.?|padj\\.?)', names(out))])
  upfront <- intersect(upfront, names(out))
  setcolorder(out, upfront)

  if (!as.dt) setDF(out)
  out
}

#' @export
#' @rdname results
#'
#' @param names The results you want to cbind together for the output, will
#'   default to all of them.
results <- function(x, names=resultNames(x), stats.only=TRUE,
                    rank.by=c('pval', 'logFC', 't'),
                    add.suffix=length(names) > 1L, as.dt=FALSE) {
  stopifnot(is(x, 'MultiGSEAResult'))
  if (is.null(resultNames(x)) || length(resultNames(x)) == 0L) {
    ## No methods were run, you can only return geneset stats
    if (!is.null(names) || length(names) > 0L) {
      ## User is asking for something that is not there
      stop("No GSEA methods were run, you can only get set statistics")
    }
    warning("No GSEA methods were run, only geneset statistics have ",
            "been returned.", immediate.=TRUE)
    out <- geneSets(x, as.dt=as.dt)
  } else {
    invalidMethods(x, names)
    stopifnot(isSingleLogical(stats.only))
    rank.by <- match.arg(rank.by)
    if (length(names) > 1L && add.suffix == FALSE) {
      add.suffix <- TRUE
      warning("Forcing method suffix to generated result columns",
              immediate.=TRUE)
    }

    ## Idiomatic data.table, non-idiomatic R
    out <- copy(geneSets(x, as.dt=TRUE))
    for (name in names) {
      res <- result(x, name, stats.only, rank.by, add.suffix, as.dt=TRUE)
      for (col in setdiff(names(res), names(out))) {
        out[, (col) := res[[col]]]
      }
    }

    # let's get the pvalue related columns up front
    upfront <- c(
      "collection", "name", "active", "N", "n",
      names(out)[grepl('^(pval\\.?|padj\\.?)', names(out))])
    upfront <- intersect(upfront, names(out))
    setcolorder(out, upfront)

    if (!as.dt) setDF(out)
  }
  out
}

#' Summary of geneset level results at a specified FDR
#'
#' Generates a table to indicate the number of genesets per collection that
#' pass a given FDR. The table provides separate groups of rows for each of
#' the `methods` run in the [multiGSEA()] call that generated that
#' generated `x`.
#'
#' @export
#' @rdname results
#'
#' @param x A [MultiGSEAResult()] object.
#' @param names the names of the GSEA methods to be reported. By default,
#'   this function will display results for all methods.
#' @param max.p The maximum padj value to consider a result significant
#' @param p.col use padj or padj.by.collection?
#'
#' @return a data.table that summarizes the significant results per method
#'   per collection for the GSEA that was run
tabulateResults <- function(x, names=resultNames(x), max.p=0.20,
                            p.col=c('padj', 'padj.by.collection', 'pval'),
                            as.dt=FALSE) {
  stopifnot(is(x, 'MultiGSEAResult'))
  invalidMethods(x, names)
  stopifnot(isSingleNumeric(max.p))
  p.col <- match.arg(p.col)
  res <- lapply(names, function(wut) {
    r <- result(x, wut, as.dt=TRUE)
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
    pcol <- NULL # silence R CMD check NOTEs
    r[, {
      list(method=wut, geneset_count=length(pcol),
           sig_count=sum(pcol <= max.p, na.rm=TRUE),
           sig_up=sum(pcol <= max.p & mean.logFC.trim > 0, na.rm=TRUE),
           sig_down=sum(pcol <= max.p & mean.logFC.trim < 0, na.rm=TRUE))
    }, by='collection']
  })
  out <- rbindlist(res)
  if (!as.dt) setDF(out)
  out
}


setMethod("show", "MultiGSEAResult", function(object) {
  msg <- paste("multiGSEA result (max FDR by collection set to 20%)",
               "---------------------------------------------------", sep='\n')
  cat(msg, "\n")
  if (length(resultNames(object)) == 0) {
    cat("No GSEA methods were run, only geneset level statistics calculated")
  } else {
    base::print.data.frame(as.data.frame(tabulateResults(object, max.p=0.20)))
  }
  cat("\n")
})

#' Assembles a matrix of nominal or adjusted pvalues from a multiGSEA result
#'
#' You might want a matrix of pvalues (or FDRs) for the gene sets across all
#' GSEA methods you tried. I think I did, once, so here it is.
#'
#' @export
#' @param x A [MultiGSEAResult()] object.
#' @param names the entries from `resultNames(x)` that you want to include
#'   in the matrix. By default we take all of them.
#' @param pval Are we testing pvalues or adjusted pvalues?
#' @return A matrix of the desired pvalues for all genesets
#'
#' @examples
#' # vm <- exampleExpressionSet(do.voom=TRUE)
#' # gdb <- exampleGeneSetDb()
#' # mg <- multiGSEA(gdb, vm, vm$design, 'tumor', methods=c('cameraPR'))
#' mg <- exampleMultiGSEAResult()
#' pm <- p.matrix(mg)
p.matrix <- function(x, names=resultNames(x),
                     pcol=c('padj', 'padj.by.collection', 'pval')) {
  stopifnot(is(x, 'MultiGSEAResult'))
  invalidMethods(x, names)
  pcol <- match.arg(pcol)
  res <- results(x, names, add.suffix=TRUE, as.dt=TRUE)
  regex <- sprintf('^%s\\.', pcol)
  col.idx <- grep(regex, names(res))
  if (pcol != 'padj.by.collection') {
    col.idx <- setdiff(col.idx, grep('padj.by.collection', names(res)))
  }
  out <- as.matrix(res[, col.idx, with=FALSE])
  colnames(out) <- sub(regex, '', colnames(out))
  out
}
