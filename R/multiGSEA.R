##' Performs several GSEA analysis over an experiment
##'
##' multiGSEA is wrapper function which delegates GSEA analyses to several
##' different "workers," which are implemented in different packages. Currently
##' this primarily wraps the GSEA functinos from limma, but I intend to add
##' \code{npGSEA} as a supported anlysis soon. Refer to the \code{method}
##' documentation to see which GSEA methods are supported.
##'
##' This function will write several results to an output directory if a valid
##' value for \code{outdir} is provided. Furthermore, the results of these
##' analysis will be cached so that time is saved on subsequent calls.
##'
##' @section Caching:
##'
##'   (NOTE: this documentation needs further work)
##'
##'   Running a GSEA can be a costly propositin, especially when using methods
##'   that rely on permutation, such as (m)roast. When \code{outdir} is set, the
##'   inputs and intermediate results from each GSEA method are saved in a
##'   directory tree that is rooted at \code{outdir}. This parameter works
##'   intimately with \code{use.cache} and \code{plots.generate}.
##'
##'   While this can be helpful for manually checking results as they are
##'   generated, you can also use it to speed up future re-analyses of the data
##'   -- perhaps you want to add one more method to the GSEA comparison that you
##'   didn't include originally. The result of each analysis method is stored
##'   in the \code{cache} directory encoded with the contrast that was tested.
##'   If that method is run again, and \code{use.cache} is set to \code{TRUE},
##'   then the result will be loaded if it is found instead of recalculating
##'   the result again.
##'
##'   If \code{force.reeval} is set to \code{TRUE} then the results of a method
##'   or plot will be redone even if they are found to already exist in the
##'   cache.
##'
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
##'   }
##' @param outdir A character vector indicating the directory to use inorder to
##'   save/cache the results from the analysis. This is required if
##'   \code{use.cache=TRUE} or if you want to save intermediate results as they
##'   are being generated.
##' @param use.cache Logical to indicate whether or not to take advantage of the
##'   poor-man's caching the multiGSEA does. See the "Caching" section below for
##'   further details, but in short the function will look for a saved result
##'   for a given \code{method} before computing it.
##' @param ... The arguments are passed down into the various geneset analysis
##'   functions.
##'
##' @return A \code{MultiGSEAResult}
multiGSEA <- function(gsd, x, design=NULL, contrast=NULL,
                      methods=c('camera'), outdir=NULL, use.cache=TRUE,
                      cache.inputs=TRUE, keep.outdir.onerror=TRUE,
                      feature.min.logFC=1, feature.max.padj=0.10,
                      trim=0.10, ...) {
  if (!is(gsd, 'GeneSetDb')) {
    stop("GeneSetDb required")
  }

  ## We know all goes well when on.exit() sees that this variable is set to TRUE
  finished <- FALSE

  ## Perhaps we should use some `.call <- match.call()` mojo for some automated
  ## cached filename generation in the future
  ## ---------------------------------------------------------------------------
  ## Argument sanity checking
  .unsupportedGSEAmethods(methods)
  if (is.null(outdir)) {
    outdir.created <- FALSE
    noutdir <- character()
    if (missing(use.cache)) {
      use.cache <- FALSE
    } else if (use.cache) {
      stop("Can't use.cache without an `outdir`")
    }
    if (missing(cache.inputs)) {
      cache.inputs <- FALSE
    } else if (cache.inputs) {
      stop("Can't cache.inputs without an `outdir`")
    }
  } else {
    if (!is.character(outdir)) {
      stop("`outdir` must be path to a directory used ot save results")
    }
    outdir.created <- .initOutDir(outdir)
    noutdir <- normalizePath(outdir)
  }

  if (length(methods) == 0) {
    stop("0 GSEA methods provided")
  }

  on.exit({
    if (!finished) {
      warning("An error in `multiGSEA` stopped it from finishing ...",
              immediate.=TRUE)
    }
    ## if (outdir.created && dir.exists(outdir) && !keep.outdir.onerror) {
    ##   unlink(outdir, recursive=TRUE)
    ## }
  })

  ## ---------------------------------------------------------------------------
  ## Sanitize / cache inputs
  args <- list(gsd=gsd, x=x, design=design, contrast=contrast)
  inputs <- validateInputs(x, design, contrast, methods,
                           require.x.rownames=TRUE)
  x <- inputs$x
  design <- inputs$design
  contrast <- inputs$contrast
  gsd <- conform(gsd, x)

  if (cache.inputs) {
    inputs.fn <- file.path(noutdir, 'cache', 'inputs.rda')
    validated <- list(gsd=gsd, x=x, design=design, contrast=contrast)
    in.list <- c(input=args, validated=inputs)
    save(in.list, file=inputs.fn)
  }

  ## ---------------------------------------------------------------------------
  ## Run the analyses
  logFC <- calculateIndividualLogFC(x, design, contrast, ...)
  logFC <- within(logFC, {
    significant <- abs(logFC) >= feature.min.logFC & padj <= feature.max.padj
    direction <- ifelse(logFC > 0, 'up', 'down')
  })

  results <- sapply(methods, function(method) {
    fn <- getFunction(paste0('do.', method))
    tryCatch({
      dt <- fn(gsd, x, design, contrast, outdir=noutdir, use.cache=use.cache,
               logFC=logFC, feature.min.logFC=feature.min.logFC,
               feature.max.padj=feature.max.padj, ...)
      dt[, padj.by.collection := p.adjust(pval, 'BH'), by='collection']
    }, error=function(e) NULL)
  }, simplify=FALSE)
  failed <- sapply(results, is.null)
  if (any(failed)) {
    warning("The following GSEA methods failed: ",
            paste(names(results)[failed], collapse=','))
  }

  setkeyv(logFC, 'featureId')

  out <- .MultiGSEAResult(gsd=gsd, results=results, logFC=logFC, outdir=noutdir)
  gs.stats <- geneSetFeatureStatistics(out, ...)
  out@gsd@table <- merge(out@gsd@table, gs.stats, by=key(out@gsd@table))
  finished <- TRUE
  out
}

##' Get the GeneSetDb from MultiGSEAResult
##'
##' @export
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

setMethod("featureIds", c(x="MultiGSEAResult"),
function(x, i, j, value=c('x.id', 'featureId'), ...) {
  value <- match.arg(value)
  featureIds(geneSetDb(x), i, j, value=value, fetch.all=FALSE,
             active.only=TRUE)
})

##' Summarize the statistics for each feature at the geneset level
##'
##' @param x A \code{multiGSEAResult} object
##' @param feature.min.logFC used with \code{feature.max.padj} to identify
##'   the individual features that have a significant change.
##' @param feature.max.padj used with \code{feature.min.logFC} to identify
##'   the individual features that have a significant change.
##' @param trim passed to \code{geneSetFeatureStatistics} for trimmed mean
##'   calculations.
##'
##' @return A data.table with statistics on the effect size shift for the gene
##'   set as well as numbers of members that shift up or down.
geneSetFeatureStatistics <- function(x, feature.min.logFC=1,
                                     feature.max.padj=0.10, trim=0.10,
                                     ...) {
  stopifnot(is(x, 'MultiGSEAResult'))
  lfc <- logFC(x)
  annotate.lfc <- !missing(feature.min.logFC) ||
    !missing(feature.max.padj) ||
    !all(c('significant', 'direction') %in% names(lfc))
  if (annotate.lfc) {
    lfc <- within(lfc, {
      significant <- abs(logFC) >= feature.min.logFC & padj <= feature.max.padj
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
    list(n.sig=sum(is.sig), n.neutral=sum(!is.sig),
         n.up=sum(up), n.down=sum(down),
         n.sig.up=sum(up & is.sig), n.sig.down=(sum(down & is.sig)),
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

##' Fetch names of GSEA methods taht were run
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
  stopifnot(is.character(names) && length(names) >= 1L)
  bad.names <- setdiff(names, resultNames(x))
  if (length(bad.names) && as.error) {
    stop("Illegal result names queried: ", paste(bad.names, collapse=','))
  }
  bad.names
}


##' Extract the result from a given enrichment test from a MultiGSEAResult
##'
##' @export
##'
##' @param x MultiGSEAResult
##' @param name the names of the results desired
##' @param stats.only logical, set to \code{FALSE} if you want to return all
##'   (column-wise) data for each result. By default only the pvalues,
##'   adjusted pvalues, and rank are returned.
##' @param rank.by the statistic to use to append a \code{rank} column for the
##'   geneset result
##' @param add.suffix If \code{TRUE}, adds \code{.name} as a suffix to the
##'   columns of the \code{method}-specific statistics returned.
##'
##' @return a data.table with the results from the requested method.
result <- function(x, name, stats.only=FALSE, rank.by=c('pval', 'logFC', 't'),
                   add.suffix=FALSE) {
  stopifnot(is(x, 'MultiGSEAResult'))
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
      res.cols <- c('pval', 'padj', 'padj.by.collection')
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
                  logFC=rank(-abs(out$mean.logFC.trim)),
                  t=rank(-abs(out$mean.t.trim)),
                  pval=rank(out[[pval.col]]))
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
                    rank.by=c('pval', 'logFC', 't'),
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
summarizeResults <- function(x, names=resultNames(x), max.p=0.30,
                             p.col=c('padj', 'padj.by.collection', 'pval')) {
  stopifnot(is(x, 'MultiGSEAResult'))
  invalidMethods(x, names)
  stopifnot(isSingleNumeric(max.p))
  p.col <- match.arg(p.col)
  res <- lapply(resultNames(x), function(wut) {
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

setMethod("show", "MultiGSEAResult", function(object) {
  msg <- paste("multiGSEA result (max FDR by collection set to 30%)",
               "---------------------------------------------------", sep='\n')
  cat(msg, "\n")
  data.table:::print.data.table(summarizeResults(object, max.p=0.30))
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

##' Initializes the output directory for the multiGSEA call.
##'
##' If outdir is NULL, then no directory is checked/created. This also implies
##' that creating plots is not possible.
##'
##' @param A character vector pointing to a directory to check/create
##' @return TRUE if the output directory was created, otherwise FALSE (it might
##' already exist).
.initOutDir <- function(outdir) {
  sub.dirs <- c('data', 'cache')
  if (is.null(outdir)) {
    return(FALSE)
  }
  if (!is.character(outdir) && length(outdir) != 1) {
    stop("character required for `outdir`")
  }
  outdir.created <- FALSE
  if (dir.exists(outdir)) {
    if (!dir.writable(outdir)) {
      stop("Can't write to output directory: ", outdir)
    }
  } else {
    pdir <- dirname(outdir)
    if (!dir.exists(pdir)) {
      stop("Path to outdir does not exist: ", pdir)
    }
    if (!dir.writable(pdir)) {
      stop("Can't create output directory in: ", pdir)
    }
    dir.create(outdir)
    outdir.created <- TRUE
  }
  for (sd in sub.dirs) {
    dir.create(file.path(outdir, sd))
  }

  outdir.created
}

