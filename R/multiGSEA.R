##' Performs several GSEA analysis over an experiment
##'
##' multiGSEA is wrapper function which delegates GSEA analyses to different
##' "workers" that implement each of the GSEA methods specified in the
##' \code{methods} argument. By default, this only gathers geneset level
##' statistics without running any of the GSEA methods in particular.
##'
##'
##' Currently this primarily wraps the GSEA functions from limma, but others
##' will be added over time. Refer to the \code{method} documentation to see
##' which GSEA methods are supported.
##'
##' logFC statistics are always calcualted via a limma pipeline, so if \code{x}
##' is a \code{DGEList}, it is first "voom"ed, then processed.
##'
##' @details This function will write several results to an output directory if
##'   a valid value for \code{outdir} is provided. Furthermore, the results of
##'   these analysis will be cached so that time is saved on subsequent calls.
##'
##' @export
##' @import limma
##'
##' @seealso
##'   \code{\link[limma]{camera}},
##'   \code{\link[limma]{roast}},
##'   \code{\link[limma]{geneSetTest}}
##'
##' @param gsd The \code{\link{GeneSetDb}} to run enrichment over.
##' @param x An ExpressoinSet-like object
##' @param design A design matrix for the study
##' @param contrast The contrast of interest to analyze
##' @param methods A character vector indicating the GSEA methods you want to
##'   run. This can be any combination of:
##'   \enumerate{
##'     \item{camera} {A \code{\link[limma]{camera}} analysis}
##'     \item{roast} {A \code{\link[limma]{roast}} analysis}
##'     \item{gsd} {A \code{\link[limma]{geneSetTest} analysis}}
##'     \item{hyperGeometricTest}{
##'       Tests for enrichment of differentially expressed features in each
##'       geneset. Thresholds of differential expression for the individual
##'       features of a geneset are defind by the \code{feature.min.logFC} and
##'       \code{feature.max.padj} parameters.
##'     }
##'   }
##'   If no methods are specified, then only geneset level statistics, such as
##'   the JG score and the mean logFCs t-statistics, of the features in each
##'   geneset are calculated.
##' @param feature.min.logFC The minimum logFC required for an individual
##'   feature (not geneset) to be considered differentialy expressed. Used in
##'   conjunction with \code{feature.max.padj} primarily for summarization
##'   of genesets (by \code{\link{geneSetsStats}}), but can also be
##'   used by GSEA methods that require differential expression calls at the
##'   individual feature level, like \code{hyperGeometricTest}.
##' @param feature.max.padj The maximum adjusted pvalue used to consider an
##'   individual feature (not geneset) to be differentially expressed. Used in
##'   conjunction with \code{feature.min.logFC}.
##' @param trim The amount to trim when calculated trimmed \code{t} and
##'   \code{logFC} statistics for each geneset. This is passed down to the
##'   \code{\link{geneSetsStats}} function.
##' @param ... The arguments are passed down into the various geneset analysis
##'   functions.
##'
##' @return A \code{MultiGSEAResult}
multiGSEA <- function(gsd, x, design=NULL, contrast=NULL,
                      methods=NULL, feature.min.logFC=1, feature.max.padj=0.10,
                      trim=0.10, ...) {
  if (!is(gsd, 'GeneSetDb')) {
    if (is(gsd, 'GeneSetCollection') || is(gsd, 'GeneSet')) {
      stop("A GeneSetDb is required. GeneSetCollections can be can be ",
           "converted to a GeneSetDb via the GeneSetDb constructor, or a call ",
           "to as(gene_set_collection, 'GeneSetDb). See ?GeneSetDb for more ",
           "info")
    }
    stop("GeneSetDb required. Please see `?GeneSetDb` for ways to turn your ",
         "gene sets into a GeneSetDb")
  }

  if (missing(methods) || length(methods) == 0) {
    methods <- 'logFC'
  }
  ## We know all goes well when on.exit() sees that this variable is set to TRUE
  finished <- FALSE
  on.exit({
    if (!finished) {
      warning("An error in `multiGSEA` stopped it from finishing ...",
              immediate.=TRUE)
    }
  })

  ## ---------------------------------------------------------------------------
  ## Argument sanity checking

  ## error-out if illegal methods were specified here
  .unsupportedGSEAmethods(methods)

  if (length(methods) == 0) {
    stop("0 GSEA methods provided")
  }

  ## ---------------------------------------------------------------------------
  ## Sanitize / cache inputs
  args <- list(gsd=gsd, x=x, design=design, contrast=contrast)
  inputs <- validateInputs(x, design, contrast, methods,
                           require.x.rownames=TRUE)

  x <- inputs$x
  design <- inputs$design
  contrast <- inputs$contrast

  if (is(x, 'DGEList')) {
    vm <- voom(x, design, plot=FALSE)
  } else {
    vm <- x
  }

  if (!is.conformed(gsd, x)) {
    gsd <- conform(gsd, x)
  }

  ## ---------------------------------------------------------------------------
  ## Run the analyses
  logFC <- calculateIndividualLogFC(vm, design, contrast, ...)
  logFC <- within(logFC, {
    significant <- abs(logFC) >= feature.min.logFC & padj <= feature.max.padj
    direction <- ifelse(logFC > 0, 'up', 'down')
  })

  ## the 'logFC' method is just a pass through -- we don't call it if it was
  ## provided
  methods <- setdiff(methods, 'logFC')
  results <- sapply(methods, function(method) {
    fn <- getFunction(paste0('do.', method))
    tryCatch({
      dt <- fn(gsd, x, design, contrast, logFC=logFC,
               feature.min.logFC=feature.min.logFC,
               feature.max.padj=feature.max.padj,
               vm=vm, ...)
      dt[, padj.by.collection := p.adjust(pval, 'BH'), by='collection']
    }, error=function(e) NULL)
  }, simplify=FALSE)
  failed <- sapply(results, is.null)
  if (any(failed)) {
    warning("The following GSEA methods failed: ",
            paste(names(results)[failed], collapse=','))
  }

  setkeyv(logFC, 'featureId')

  out <- .MultiGSEAResult(gsd=gsd, results=results, logFC=logFC)
  gs.stats <- geneSetsStats(out, feature.min.logFC=feature.min.logFC,
                            feature.max.padj=feature.max.padj,
                            trim=trim)
  out@gsd@table <- merge(out@gsd@table, gs.stats, by=key(out@gsd@table))
  finished <- TRUE
  out
}

##' Get the GeneSetDb from MultiGSEAResult
##'
##' @export
##' @rdname geneSetDb-accessor
##'
##' @param x \code{MultiGSEAResult}
##' @return The \code{GeneSetDb}
geneSetDb <- function(x) {
  stopifnot(is(x, 'MultiGSEAResult'))
  x@gsd
}

setMethod("geneSets", c(x="MultiGSEAResult"),
function(x, ...) {
  geneSets(geneSetDb(x), active.only=TRUE)
})

## Here's your chance to vectorize this. Once that's done push this code down
## into geneSetUrl,GeneSetDb
setMethod("geneSetURL", c(x="MultiGSEAResult"), function(x, i, j, ...) {
  geneSetURL(geneSetDb(x), i, j, ...)
})

setMethod("geneSetCollectionURLfunction", "MultiGSEAResult", function(x, i, ...) {
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
geneSetFeatureStats <- function(x, i, j) {
  stopifnot(is(x, 'MultiGSEAResult'))
  fids <- featureIds(x@gsd, i, j)
  subset(logFC(x), featureId %in% fids)
}

##' Summarize useful feature-level statistics across gene sets.
##'
##' This function calculates the mean and trimmed mean of the logFC and
##' t-statistics, as well as the J-G statistic
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
geneSetsStats <- function(x, feature.min.logFC=1,
                          feature.max.padj=0.10, trim=0.10) {
  stopifnot(is(x, 'MultiGSEAResult'))
  lfc <- logFC(x)

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

  gs <- geneSets(x)
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
         JG=sum(t.nona) / sqrt(length(t.nona)),
         mean.logFC=mean(stats$logFC, na.rm=TRUE),
         mean.logFC.trim=mean(stats$logFC, na.rm=TRUE, trim=trim),
         mean.t=mean(stats$t, na.rm=TRUE),
         mean.t.trim=mean(stats$t, na.rm=TRUE, trim=trim))
  }, by=do.by]
  setkeyv(out, do.by)
}

##' Extract the individual fold changes statistics for elements in the
##' expression object.
##'
##' @export
##' @param x A \code{MultiGSEAResult}
##' @return The log fold change \code{data.table}
logFC <- function(x) {
  stopifnot(is(x, 'MultiGSEAResult'))
  x@logFC
}

##' Fetch names of GSEA methods that were run from a \code{MultiGSEAResult}
##'
##' @export
##'
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
                   rank.by=c('pval', 't', 'logFC', 'JG'),
                   add.suffix=FALSE) {
  stopifnot(is(x, 'MultiGSEAResult'))
  ## no methods run?
  if (length(resultNames(x)) == 0) {
    return(results(x))
  }
  stopifnot(isSingleCharacter(name))
  invalidMethods(x, name, as.error=TRUE)
  stopifnot(isSingleLogical(stats.only))
  rank.by <- match.arg(rank.by)
  stopifnot(isSingleLogical(add.suffix))

  out <- copy(geneSets(x))

  pval.col <- 'pval'
  rank.col <- 'rank'
  if (add.suffix) {
    pval.col <- paste0('pval.', name)
    rank.col <- paste0('rank.', name)
  }

  res <- local({
    r <- x@results[[name]]
    if (!all.equal(out[, key(out), with=FALSE], r[, key(out), with=FALSE])) {
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
                  JG=rank(-abs(out$JG), ties.method="min"),
                  logFC=rank(-abs(out$mean.logFC.trim), ties.method="min"),
                  t=rank(-abs(out$mean.t.trim), ties.method="min"),
                  pval=rank(out[[pval.col]], ties.method="min"))
  out[, (rank.col) := ranks]

  out
}

##' A summary of one/some/all of the results that were run.
##'
##' If the result from more than one method is requested, the column names for
##' these results will be suffixed with the method name.
##'
##' @export
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
                    rank.by=c('pval', 'logFC', 't', 'JG'),
                    add.suffix=length(names) > 1L) {
  stopifnot(is(x, 'MultiGSEAResult'))
  invalidMethods(x, names)
  stopifnot(isSingleLogical(stats.only))
  rank.by <- match.arg(rank.by)
  if (length(names) > 1L && add.suffix == FALSE) {
    add.suffix <- TRUE
    warning("Forcing method suffix to generated result columns",
            immediate.=TRUE)
  }

  ## Idiomatic data.table, non-idiomatic R
  out <- copy(geneSets(x))
  for (name in names) {
    res <- result(x, name, stats.only, rank.by, add.suffix)
    for (col in setdiff(names(res), names(out))) {
      out[, (col) := res[[col]]]
    }
  }

  out
}

##' Create a summary table to indicate number of significant genesets per
##' collection, per method
##'
##' @export
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
    r <- result(x, wut)
    r$pcol <- r[[p.col]]
    r[, {
      list(method=wut, geneset_count=length(pcol),
           sig_count=sum(pcol <= max.p, na.rm=TRUE),
           sig_up=sum(pcol <= max.p & mean.logFC.trim > 0, na.rm=TRUE),
           sig_down=sum(pcol <= max.p & mean.logFC.trim < 0, na.rm=TRUE))
    }, by='collection']
  })
  rbindlist(res)
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
  nr <- nrow(if (did.gsea) x@results[[1]] else geneSets(x))
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

if (FALSE) {
##' @export
##' @importFrom BiocGenerics subset
setMethod("subset", "MultiGSEAResult",
function(x, subject, select, drop=FALSE, ...) {

})
}

setMethod("show", "MultiGSEAResult", function(object) {
  msg <- paste("multiGSEA result (max FDR by collection set to 30%)",
               "---------------------------------------------------", sep='\n')
  cat(msg, "\n")
  if (length(resultNames(object)) == 0) {
    cat("No GSEA methods were run, only geneset level statistics calculated")
  } else {
    data.table:::print.data.table(tabulateResults(object, max.p=0.30))
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
