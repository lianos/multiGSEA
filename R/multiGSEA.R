#' Performs a plethora of GSEA analyses over a contrast of interest.
#'
#' multiGSEA is wrapper function that delegates GSEA analyses to different
#' "workers", each of which implements the flavor of GSEA of your choosing.
#' The particular analyses that are performed are specified by the
#' `methods` argument, and these methods are fine tuned by passing their
#' arguments down through the `...` of this wrapper function.
#'
#' The bulk of the GSEA methods currently available in this package come from
#' edgeR/limma, however others are included (and are being added), as well.
#' *GSEA Methods* and *GSEA Method Parameterization* sections for more details.
#'
#' In addition to performing GSEA, this function also internally orchestrates
#' a differential expression analysis, which can be tweaked by identifying
#' the parameters in the [calculateIndividualLogFC()] function, and
#' passing them down through `...` here. The results of the differential
#' expression analysis (ie. the [limma::topTable()]) are accessible by calling
#' the [logFC()] function on the [MultiGSEAResult()] object returned from this
#' function call.
#'
#' **Please Note**: be sure to cite the original GSEA method when using
#' results generated from this function.
#'
#' @section GSEA Methods:
#' You can choose the methods you want to run over your contrast using the
#' \code{methods}, parameter. The currently available methods are:
#'
#' - `"camera"`: from [limma::camera()] (*)
#' - `"cameraPR"`: from [limma::cameraPR()]
#' - `"roast"`: from [limma::roast()] (*)
#' - `"fry"`: from [limma::fry()] (*)
#' - `"romer"`: from [limma::romer()] (*)
#' - `"geneSetTest"`: from [limma::geneSetTest()]
#' - `"goseq"`: from [goseq::goseq()]
#' - `"hyperGeometricTest"`
#' - `"fgsea"`: from [fgsea::fgsea()]
#'
#' Methods annotated with a `(*)` indicate that these methods require
#' a complete expression object, a valid design matrix, and a contrast
#' specification in order to run. These are all of the same things you need to
#' provide when performing a vanilla differential gene expression analysis.
#'
#' Methods missing a `(*)` can be run on a feature-named input vector
#' of gene level statistics which will be used for ranking (ie. a named vector
#' of logFC's or t-statistics for genes). They can also be run by providing
#' an expression, design, and contrast vector, and the appropriate statistics
#' vector will be generated internally from the t-statistics, p-values, or
#' log-fold-changes, depending on the value provided in the `score.by`
#' parameter.
#'
#' The worker functions that execute these GSEA methods are functions named
#' `do.METHOD` within this package. These functions are not meant to be executed
#' directly by the user, and are therefore not exported. Look at the respective
#' method's help page (ie. if you are running `"camera"`, look at the
#' [limma::camera()] help page for full details. The formal parameters that
#' these methods take can be passed to them via the `...` in this `multiGSEA()`
#' function.
#'
#' @section GSEA Method Parameterization:
#'
#' Each GSEA method can be tweaked via a custom set of parameters. We leave the
#' documentation of these parameters and how they affect their respective GSEA
#' methods to the documentation available in the packages where they are
#' defined. The \code{multiGSEA} call simply has to pass these parameters down
#' into the \code{...} parameters here. The \code{multiGSEA} function will then
#' pass these along to their worker functions.
#'
#' **What happens when two different GSEA methods have parameters with the
#' same name?**
#'
#' Unfortunately you currently cannot provide different values for these
#' parameters. An upcoming version version of multiGSEA will support this
#' feature via slightly different calling semantics. This will also allow the
#' caller to call the same GSEA method with different parameterizations so that
#' even these can be compared against each other.
#'
#' @section Differential Gene Expression:
#'
#' When the `multiGSEA()` call is given an expression matrix, design, and
#' contrast, it will also internally orchestrate a gene level differential
#' expression analysis. Depending on the type of expression object passed in
#' via `x`, this function will guess on the best method to use for this
#' analysis.
#'
#' If `x` is a \code{DGEList}, then ensure that you have already called
#' [edgeR::estimateDisp()] on `x` and edgeR's quasilikelihood framework will be
#' used, otherwise we'll use limma (note that `x` can be an `EList` run through
#' `voom()`, `voomWithQuailityWeights()`, or when where you have leveraged
#' limma's [limma::duplicateCorrelation()] functionality, even.
#'
#' The parameters of this differential expression analysis can also be
#' customized. Please refer to the [calculateIndividualLogFC()] function for
#' more information. The multiGSEA `use.treat`, `feature.min.logFC`,
#' `feature.max.padj`, as well as the `...` parameters are passed down to that
#' funciton.
#'
#' @export
#' @importFrom BiocParallel bplapply SerialParam bpparam
#'
#' @param gsd The [GeneSetDb()] that defines the gene sets of interest.
#' @param x An ExpressoinSet-like object
#' @param design A design matrix for the study
#' @param contrast The contrast of interest to analyze. This can be a column
#'   name of `design`, or a contrast vector which performs "coefficient
#'   arithmetic" over the columns of `design`. The `design` and `contrast`
#'   parameters are interpreted in *exactly* the same way as the same parameters
#'   in limma's [limma::camera()] and [limma::roast()] methods.
#' @param methods A character vector indicating the GSEA methods you want to
#'   run. Refer to the GSEA Methods section for more details.
#'   If no methods are specified, only differential gene expression and geneset
#'   level statistics for the contrast are computed.
#' @param use.treat should we use limma/edgeR's "treat" functionality for the
#'   gene-level differential expression analysis?
#' @param feature.min.logFC The minimum logFC required for an individual
#'   feature (not geneset) to be considered differentialy expressed. Used in
#'   conjunction with `feature.max.padj` primarily for summarization
#'   of genesets (by [geneSetsStats()], but can also be used by GSEA methods
#'   that require differential expression calls at the individual feature level,
#'   like [goseq()].
#' @param feature.max.padj The maximum adjusted pvalue used to consider an
#'   individual feature (not geneset) to be differentially expressed. Used in
#'   conjunction with `feature.min.logFC`.
#' @param trim The amount to trim when calculated trimmed `t` and
#'   `logFC` statistics for each geneset. This is passed down to the
#'   [geneSetsStats()] function.
#' @param verbose make some noise during execution?
#' @param ... The arguments are passed down into
#'   [calculateIndividualLogFC()] and the various geneset analysis functions.
#' @param .parallel by default, `.parallel=FALSE` runs each GSEA in a
#'   serial manner. If `.parallel=TRUE`, the GSEA execution loop is
#'   parallelized using the *BiocParallel* package. Note that you might want to
#'   remove unnecessary large objects from your workspace when this is `TRUE`
#'   because R will likely want to copy them down into your worker threads.
#' @param BPPARAM a *BiocParallel* parameter definition, like one generated from
#'   [BiocParallel::MulticoreParam()()], or [BiocParallel::BatchtoolsParam()],
#'   for instance, which is passed down to [BiocParallel]::bplapply()]. If not
#'   specified and `.parallel = TRUE`, then the [BiocParallel::bpparam()] object
#'   will be used. If `.parallel = FALSE`, this parameter is explicitly ignored
#'   and replaced with a [BiocParallel]::SerialParam()] object.
#' @return A [MultiGSEAResult()] which holds the results of all the analyses
#'   specified in the `methods` parameter.
#'
#' @examples
#' vm <- exampleExpressionSet(do.voom=TRUE)
#' gdb <- exampleGeneSetDb()
#' mg <- multiGSEA(gdb, vm, vm$design, 'tumor',
#'                 methods=c('camera', 'fry'),
#'                 ## customzie camera parameter:
#'                 inter.gene.cor=0.04)
#' resultNames(mg)
#' res.camera <- result(mg, 'camera')
#' res.fry <- result(mg, 'fry')
#' res.all <- results(mg)
multiGSEA <- function(gsd, x, design=NULL, contrast=NULL,
                      methods=NULL, use.treat=FALSE,
                      feature.min.logFC=if (use.treat) log2(1.25) else 1,
                      feature.max.padj=0.10, trim=0.10, verbose=FALSE, ...,
                      .parallel=FALSE, BPPARAM=bpparam()) {
  if (!is(gsd, "GeneSetDb")) {
    if (is(gsd, "GeneSetCollection") || is(gsd, "GeneSet")) {
      stop("A GeneSetDb is required. GeneSetCollections can be converted to a ",
           "GeneSetDb via the GeneSetDb constructor, or a call to ",
           "`as(gene_set_collection, 'GeneSetDb')`. See ?GeneSetDb for more ",
           "details.")
    }
    stop("GeneSetDb required. Please see `?GeneSetDb` for ways to turn your ",
         "gene sets into a GeneSetDb")
  }
  if (missing(methods) || length(methods) == 0) {
    methods <- "logFC"
  }

  stopifnot(
    is.numeric(feature.min.logFC) && length(feature.min.logFC) == 1L,
    is.numeric(feature.max.padj) && length(feature.max.padj) == 1L,
    is.numeric(trim) && length(trim) == 1L)

  # ----------------------------------------------------------------------------
  # Argument sanity checking and input sanitization
  args <- list(gsd=gsd, x=x, design=design, contrast=contrast)
  inputs <- validateInputs(x, design, contrast, methods,
                           require.x.rownames=TRUE, ...)

  x <- inputs$x
  design <- inputs$design
  contrast <- inputs$contrast

  if (!is.conformed(gsd, x)) {
    gsd <- conform(gsd, x, ...)
  }

  # ----------------------------------------------------------------------------
  # Run the analyses

  # First calculate differential expression statistics, or wrap a pre-ranked
  # vector into a data.frame returned by an internal dge analysis
  treat.lfc <- if (use.treat) feature.min.logFC else NULL
  logFC <- calculateIndividualLogFC(x, design, contrast, use.treat=use.treat,
                                    treat.lfc=treat.lfc, verbose=verbose, ...,
                                    as.dt=TRUE)
  test_type <- attr(logFC, "test_type")

  logFC[, significant := {
    if (test_type == "anova") {
      padj <= feature.max.padj
    } else {
      padj <= feature.max.padj & abs(logFC) >= feature.min.logFC
    }
  }]

  # the 'logFC' method is just a pass through -- we don't call it if it was
  # provided
  methods <- setdiff(methods, "logFC")

  # Let's do this!
  results <- list()
  if (length(methods) > 0L) {
    # I'm being too clever here. The loop that calls the GSEA methods catches
    # errors thrown during iteration. I'm putting some code here to eat those
    # error so that the other GSEA methods that can finish.
    finished <- FALSE
    on.exit({
      if (!finished) {
        warning("An error in `multiGSEA` stopped it from finishing ...",
                immediate.=TRUE)
      }
    })

    ## Many methods create a geneset to rowname/index vector. Let's run it once
    ## here and pass it along
    gs.idxs <- as.list(gsd, active.only = TRUE, value = "x.idx")

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

    failed <- sapply(res1, function(res) is.null(res[[1L]]))
    if (any(failed)) {
      warning("The following GSEA methods failed and are removed from the ",
              "downstream result: ",
              paste(names(res1)[failed], collapse=','), "\n")
      res1 <- res1[!failed]
    }

    results <- unlist(res1, recursive=FALSE)
    names(results) <- sub('\\.all$', '', names(results))
  }

  out <- .MultiGSEAResult(gsd=gsd, results=results, logFC=logFC)
  gs.stats <- geneSetsStats(out, feature.min.logFC=feature.min.logFC,
                            feature.max.padj=feature.max.padj,
                            trim=trim, as.dt=TRUE)
  axe.gsd.cols <- setdiff(colnames(gs.stats), c('collection', 'name'))
  axe.gsd.cols <- intersect(axe.gsd.cols, colnames(out@gsd@table))
  new.table <- copy(out@gsd@table)
  ## Remove any columns in gs.stats that are already in the GeneSetDb@table
  ## (ie. if we got a GeneSetDb from a previous MultiGSEAResult and we don't
  ## remove these columns, you will get thigns like mean.logFC.x and
  ## mean.logFC.y
  if (length(axe.gsd.cols)) {
    name <- NULL
    for (name in axe.gsd.cols) new.table[, c(name) := NULL]
  }
  out@gsd@table <- merge(new.table, gs.stats, by=key(new.table))
  finished <- TRUE
  out
}

#' Helper function that runs a single GSEA method
#'
#' this method insures that all results (even for single data.frames) are
#' returned in a list object. this allows for the goseq.all, goseq.up,
#' goseq.down hacks from a single goseq call
#'
#' @noRd
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
  if (!isTRUE(attr(res, 'mgunlist', TRUE))) {
    res <- list(all=res)
  }
  out <- res
  if (verbose) {
    message("... ", fn.name, " finishd without error.")
  }
  out
}
