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
##' @importFrom parallel mclapply
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
                      methods=NULL, use.treat=FALSE,
                      feature.min.logFC=if (use.treat) log2(1.25) else 1,
                      feature.max.padj=0.10,
                      trim=0.10, verbose=FALSE, ..., do.parallel=FALSE) {
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
                           require.x.rownames=TRUE, ...)

  x <- inputs$x
  design <- inputs$design
  contrast <- inputs$contrast

  if (!is.conformed(gsd, x)) {
    gsd <- conform(gsd, x)
  }

  ## ---------------------------------------------------------------------------
  ## Run the analyses
  logFC <- calculateIndividualLogFC(x, design, contrast, use.treat=use.treat,
                                    treat.lfc=feature.min.logFC,
                                    verbose=verbose, ..., .external=FALSE)
  logFC <- within(logFC, {
    significant <- if (use.treat) {
      padj <= feature.max.padj
    } else {
      abs(logFC) >= feature.min.logFC & padj <= feature.max.padj
    }
  })

  ## the 'logFC' method is just a pass through -- we don't call it if it was
  ## provided
  methods <- setdiff(methods, 'logFC')

  ## Many methods create a geneset to rowname/index vector. Let's run it once
  ## here and pass it along
  gs.idxs <- as.expression.indexes(gsd, value='x.idx')
  if (length(methods) > 0L) {
    if (do.parallel) {
      res1 <- mclapply(methods, function(method) {
        tryCatch(mg.run(method, gsd, x, design, contrast, logFC, use.treat,
                        feature.min.logFC, feature.max.padj, verbose=verbose,
                        gs.idxs=gs.idxs, ...),
                 error=function(e) list(NULL))
      })
    } else {
      res1 <- lapply(methods, function(method) {
        tryCatch(mg.run(method, gsd, x, design, contrast, logFC, use.treat,
                        feature.min.logFC, feature.max.padj, verbose=verbose,
                        ...),
                 error=function(e) list(NULL))
      })
    }
    names(res1) <- methods
    results <- unlist(res1, recursive=FALSE)
    names(results) <- sub('\\.all$', '', names(results))
    failed <- sapply(results, is.null)
    if (any(failed)) {
      warning("The following GSEA methods failed: ",
              paste(names(results)[failed], collapse=','))
    }
  } else {
    results <- list()
  }

  setkeyv(logFC, 'featureId')

  out <- .MultiGSEAResult(gsd=gsd, results=results, logFC=logFC)
  gs.stats <- geneSetsStats(out, feature.min.logFC=feature.min.logFC,
                            feature.max.padj=feature.max.padj,
                            trim=trim, .external=FALSE)
  out@gsd@table <- merge(out@gsd@table, gs.stats, by=key(out@gsd@table))
  finished <- TRUE
  out
}

mg.run <- function(method, gsd, x, design, contrast, logFC=NULL,
                   use.treat=TRUE, feature.min.logFC=log2(1.25),
                   feature.max.padj=0.10, verbose=FALSE, ...) {
  fn.name <- paste0('do.', method)
  if (verbose) {
    message("... calling: ", fn.name)
  }

  fn <- getFunction(fn.name)
  res <- fn(gsd, x, design, contrast, logFC=logFC,
            use.treat=use.treat, feature.min.logFC=feature.min.logFC,
            feature.max.padj=feature.max.padj, verbose=verbose, ...)
  if (is.data.table(res)) {
    res <- list(all=res)
  }
  out <- lapply(res, function(dt) {
    pcols <- grepl('^pval\\.?', names(dt))
    if (any(pcols)) {
      for (i in which(pcols)) {
        pcol <- names(dt)[i]
        pname <- paste0(sub('pval', 'padj', pcol), '.by.collection')
        padjs <- p.adjust(dt[[pcol]], 'BH')
        dt[, (pname) := padjs]
      }
    }
    dt
  })
  if (verbose) {
    message("... ", fn.name, " finishd without error.")
  }
  out
}
