##' Performs a plethora of GSEA analyses over a contrast of interest.
##'
##' multiGSEA is wrapper function that delegates GSEA analyses to different
##' "workers", each of which implements the flavor of GSEA of your choosing.
##' The particular analyses that are performed are specified by the
##' \code{methods} argument, and these methods are fine tuned by passing their
##' arguments down through the \code{...} of this wrapper function.
##'
##' Currently this primarily wraps many of the GSEA methods available in
##' edgeR/limma as well as goseq and hypergeometric testing functionality.
##' Others will be added over time. Refer to the \emph{GSEA Methods} and
##' \emph{GSEA Method Parameterization} sections for more details.
##'
##' In addition to performing GSEA, this function also internally orchestrates
##' a differential expression analysis, which can be tweaked by idnetifying
##' the parameters in the \code{\link{calculateIndividualLogFC}} function, and
##' passing them down through \code{...} here. The results of the differential
##' expression analysis are accessible by calling the \code{logFC} function on
##' the \code{\link{MultiGSEAResult}} object returned from this function call.
##'
##' \emph{Please Note}: be sure to cite the original GSEA method when using
##' results generated from it.
##'
##' @section GSEA Methods:
##'
##' You can choose the methods you want to run over your contrast using the
##' \code{methods}, parameter. The currently available methods are:
##'
##' \enumerate{
##'   \item \code{\link[limma]{camera}(*)}
##'   \item \code{\link[limma]{roast}(*)}
##'   \item \code{\link[limma]{fry}(*)}
##'   \item \code{\link[limma]{romer}(*)}
##'   \item \code{\link[limma]{geneSetTest}}
##'   \item \code{\link{goseq}}
##'   \item \code{\link{hyperGeometricTest}}
##' }
##'
##' Methods annotated with a \code{(*)} indicate that these methods require
##' a complete expression object, a valid design matrix, and a contrast
##' specification in order to run. These are all of the same things you need to
##' provide when performing a vanilla differential gene expression analysis.
##'
##' Methods missing a \code{(*)} can be run on a feature-named input vector
##' of gene level statistics which will be used for ranking (ie. a named vector
##' of logFC's or t-statistics for genes).
##'
##' The worker functions that execute these GSEA methods are functions named
##' \code{do.METHOD}. These functions aren't executed, but you can find help
##' on them via the console, via \code{?do.goseq}, or \code{?do.camera},
##' for example.
##'
##' @section GSEA Method Parameterization:
##'
##' Each GSEA method can be tweaked via a custom set of parameters. We leave the
##' documentation of these parameters and how they affect their respective GSEA
##' methods to the documentation available in the packages where they are
##' defined. The caller simply has to pass these parameters down
##' into the \code{...} parameters here. The \code{multiGSEA} function will then
##' pass these along to their worker functions.
##'
##' \emph{What happens when two different GSEA methods have parameters with the
##' same name?}
##'
##' Unfortunately you currently cannot provide different values for these
##' parameters. An upcoming version version of multiGSEA will support this
##' feature via slightly different calling semantics. This will also allow the
##' caller to call the same GSEA method with different parameterizations so that
##' even these can be compared against each other.
##'
##' @section Differential Gene Expression:
##'
##' When the \code{multiGSEA} call is given an expression matrix, design, and
##' contrast, it will also internally orchestrate a gene level differential
##' expression analysis. Depending on the type of expression object passed in
##' via \code{x}, this function will guess on the best method to use for this
##' analysis.
##'
##' If \code{x} is a \code{DGEList}, then ensure that you have already called
##' \code{\link[edgeR]{estimateDisp}} on \code{x} and edgeR's quasilikelihood
##' framework will be used, otherwise we'll use limma (note that \code{x} can
##' be a "voom"d \code{EList}).
##'
##' The parameters of this differential expression analysis can also be
##' customized. Please refer to the \code{\link{calculateIndividualLogFC}}
##' function for more information. The multiGSEA \code{use.treat},
##' \code{feature.min.logFC}, \code{feature.max.padj}, as well as the \code{...}
##' parameters are passed down to that funciton.
##'
##' @export
##' @importFrom BiocParallel bplapply SerialParam
##'
##' @seealso
##'   \code{\link[limma]{camera}},
##'   \code{\link[limma]{roast}},
##'   \code{\link[limma]{geneSetTest}}
##'
##' @param gsd The \code{\link{GeneSetDb}} that defines the gene sets of
##'   interest.
##' @param x An ExpressoinSet-like object
##' @param design A design matrix for the study
##' @param contrast The contrast of interest to analyze. This can be a column
##'   name of \code{design}, or a contrast vector which performs "coefficient
##'   arithmetic" over the columns of \code{design}. The \code{design} and
##'   \code{contrast} parameters are interpreted in \emph{exactly} the same
##'   way as the same parameters in limma's \code{\link[limma]{camera}} and
##'   \code{\link[limma]{roast}} methods.
##' @param methods A character vector indicating the GSEA methods you want to
##'   run. Refer to the GSEA Methods section for more details.
##'   If no methods are specified, only differential gene expressino and geneset
##'   level statistics for the contrast are computed.
##' @param use.treat should we use limma/edgeR's "treat" functionality for the
##'   gene-level differential expression analysis?
##' @param feature.min.logFC The minimum logFC required for an individual
##'   feature (not geneset) to be considered differentialy expressed. Used in
##'   conjunction with \code{feature.max.padj} primarily for summarization
##'   of genesets (by \code{\link{geneSetsStats}}), but can also be
##'   used by GSEA methods that require differential expression calls at the
##'   individual feature level, like \code{\link{goseq}}
##' @param feature.max.padj The maximum adjusted pvalue used to consider an
##'   individual feature (not geneset) to be differentially expressed. Used in
##'   conjunction with \code{feature.min.logFC}.
##' @param trim The amount to trim when calculated trimmed \code{t} and
##'   \code{logFC} statistics for each geneset. This is passed down to the
##'   \code{\link{geneSetsStats}} function.
##' @param verbose make some noise during execution?
##' @param ... The arguments are passed down into the various geneset analysis
##'   functions.
##' @param .parallel by default, \code{.parallel=FALSE} runs each GSEA in a
##'   serial manner. If \code{.parallel=TRUE}, the GSEA execution loop is
##'   parallelized using \code{BiocParallel}. Note that you might want to remove
##'   unnecessary large objects from your workspace when this is \code{TRUE}
##'   because R will likely want to copy them down into your worker threads.
##' @param BPPARAM a \code{\link[BiocParallel]{BatchJobsParam}} passed down to
##'   \code{\link[BiocParallel]{bplapply}}. If not specified, the already
##'   registerd \code{\link[BiocParallel]{BatchJobsParam}} will be used.
##'   If \code{.parallel=FALSE}, this parameter is explicitly ignored and
##'   replaced with a \code{\link[BiocParallel]{SerialParam}} object.
##' @return A \code{\link{MultiGSEAResult}} which holds the results of all the
##'   analyses specified in \code{methods}.
##'
##' @examples
##'
##' vm <- exampleExpressionSet(do.voom=TRUE)
##' gdb <- exampleGeneSetDb()
##' mg <- multiGSEA(gdb, vm, vm$design, 'tumor', methods=c('camera', 'fry'))
##' resultNames(mg)
##' res.camera <- result(mg, 'camera')
##' res.fry <- result(mg, 'fry')
##' res.all <- results(mg)
multiGSEA <- function(gsd, x, design=NULL, contrast=NULL,
                      methods=NULL, use.treat=FALSE,
                      feature.min.logFC=if (use.treat) log2(1.25) else 1,
                      feature.max.padj=0.10, trim=0.10, verbose=FALSE, ...,
                      .parallel=FALSE, BPPARAM=bpparam()) {
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
  treat.lfc <- if (use.treat) feature.min.logFC else NULL
  logFC <- calculateIndividualLogFC(x, design, contrast, use.treat=use.treat,
                                    treat.lfc=treat.lfc, verbose=verbose, ...,
                                    .external=FALSE)
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
  gs.idxs <- as.list(gsd, active.only=TRUE, value='x.idx')

  ## Let's do this!
  if (length(methods) > 0L) {
    if (!.parallel) {
      BPPARAM <- SerialParam(stop.on.error=FALSE)
    }
    stopifnot(is(BPPARAM, 'BiocParallelParam'))

    res1 <- bplapply(methods, function(method) {
      tryCatch(mg.run(method, gsd, x, design, contrast, logFC, use.treat,
                      feature.min.logFC, feature.max.padj, verbose=verbose,
                      gs.idxs=gs.idxs, ...),
               error=function(e) list(NULL))
    }, BPPARAM=BPPARAM)
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

  out <- .MultiGSEAResult(gsd=gsd, results=results, logFC=logFC)
  gs.stats <- geneSetsStats(out, feature.min.logFC=feature.min.logFC,
                            feature.max.padj=feature.max.padj,
                            trim=trim, .external=FALSE)
  out@gsd@table <- merge(out@gsd@table, gs.stats, by=key(out@gsd@table))
  finished <- TRUE
  out
}

## helper function that runs a single GSEA method
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
