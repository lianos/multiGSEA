##' Runs several GSEA methods over a dataset for a single contrast.
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
##' @export
##'
##' @import limma
##' @importFrom parallel mclapply
##'
##' @param x An ExpressoinSet-like object
##' @param gene.sets This can be several things: (a) a character vector that
##'   lists the MSigDB set id's (in conjunction with the \code{species}
##'   parameter); (b) a list of logical vectors into the rows of \code{x}. This
##'   should be a list of lists, where the top-level lists identified the
##'   "group" its child elements (genesets) belong to (like the c1, c2, etc.
##'   grouping) in the MSigDB lists; ie. the provenance of the gene set
##'   (c) An already rigged up GeneSetTable
##' @param design A design matrix for the study
##' @param contrast The contrast of interest to analyze
##' @param methods A character vector indicating the GSEA methods you want to
##'   run. This can be any combination of:
##' \enumerate{
##'   \item{camera}{A \code{limma::camera} analysis}
##'   \item{roast}{A \code{limma::(m)roast} analysis}
##'   \item{gst}{A \code{limma::geneSetTest analysis}}
##' }
##' @param outdir A character vector indicating the directory to use inorder to
##'   save/cache the results from the analysis. This is required if you want to
##'   generate plots from this analysis
##' @param plots.generate A logical indicating whether or not you want to
##'   generate "distribution shift" plots for the GSEA results. If \code{TRUE},
##'   a propoer \code{outdir} parameter is required.
##' @param plots.padj.treshold If \code{plots.generate == TRUE}, this value
##'   determines what FDR cutoff is required for a plot to be generated for a
##'   particular geneset.
##' @param species This is only required if you are utilizing the automatic
##'   construction of gene sets from MSigDB gene set IDs, ie. if
##'   \code{gene.sets} is something like \code{c('c2', 'c7')}. This can only be
##'    either "human" or "mouse".
##' @param mc.cores GSEA analyses can be parallelized **over `methods` only**
##'   by setting this value to > 1. This means that you can run a camera and
##'   roast analysis simultaneously. Unfortunately this is not ||-izing each
##'   analysis over a large number of genesets.
##' @param use.cache Logical to indicate whether or not to take advantage of
##'   the poor-man's caching the multiGSEA does.
##' @param force.reeval If \code{TRUE} no cached results will be re-used if
##'   they are found.
##' @param keep.outdir.onerror If \code{FALSE}, the \code{outdir} is removed
##'   if this function *created* it.
##' @param ... The arguments are passed down into the various geneset analysis
##'   functions.
##'
##' @return An uber data.table with the results from the GSEA \code{methods}
##'   collated into separate columns.
multiGSEA <- function(x, gene.sets, design=NULL, contrast=NULL,
                      methods=c('camera'), outdir=NULL,
                      plots.generate=FALSE, plots.padj.threshold=1,
                      species=NULL, mc.cores=1L, use.cache=FALSE,
                      force.reeval=FALSE, keep.outdir.onerror=TRUE,
                      score.by=c('logFC', 't'), ...) {
  ## Perhaps we should use some `.call <- match.call()` mojo for some automated
  ## cached filename generation in the future
  ## ---------------------------------------------------------------------------
  ## Argument sanity checking
  score.by <- match.arg(score.by)
  if (is.null(outdir)) {
    outdir.created <- FALSE
    if (missing(plots.generate)) {
      plots.generate <- FALSE
    } else if (plots.generate) {
      stop("Can't set `plots.generate=TRUE` without and `outdir`")
    }
    if (missing(use.cache)) {
      use.cache <- FALSE
    } else if (use.cache) {
      stop("Can't use.cache without an `outdir`")
    }
  } else {
    if (!is.character(outdir)) {
      stop("`outdir` must be path to a directory used ot save results")
    }
    outdir.created <- .initOutDir(outdir)
  }

  species <- match.species(species)
  if (length(methods) == 0) {
    stop("0 GSEA methods provided")
  }
  .unsupportedGSEAmethods(methods)
  finished <- FALSE

  on.exit({
    if (!finished) {
      warning("An error in `multiGSEA` stopped it from finishing ...",
              immediate.=TRUE)
    }
    if (outdir.created && dir.exists(outdir) && !keep.outdir.onerror) {
      unlink(outdir, recursive=TRUE)
    }
  })

  ## ---------------------------------------------------------------------------
  ## Sanitize / load inputs
  inputs.fn <- file.path(outdir, 'cache', 'inputs.rda')
  if (use.cache && file.exists(inputs.fn)) {
    load(inputs.fn)
  } else {
    req.x.rownames <- !inherits(gene.sets, 'GeneSetTable')
    args <- list(x=x, gene.sets=gene.sets, design=design, contrast=contrast)
    inputs <- validateInputs(x, design, contrast, methods,
                             require.x.rownames=req.x.rownames)
    x <- inputs$x
    design <- inputs$design
    contrast <- inputs$contrast

    if (!is(gene.sets, 'GeneSetTable')) {
      gst <- GeneSetTable(gene.sets, x)
    } else {
      gst <- conform(gene.sets, x)
    }

    if (!is.null(outdir)) {
      save(x, design, contrast, gst, file=inputs.fn)
    }
  }

  ## ---------------------------------------------------------------------------
  ## Run the analyses
  results <- mclapply(methods, function(method) {
    fn <- getFunction(paste0('do.', method))
    fn(x, gst, design, contrast, outdir=outdir, use.cache=use.cache, ...)
  }, mc.cores=mc.cores)
  names(results) <- methods

  ## calculate the log fold changes for all the elements in x for the given
  ## contrast we are running GSEA over
  istats <- calculateIndividualLogFC(x, design, contrast, provide=score.by, ...)

  out <- summarizeResults(x, design, contrast, gst, results,
                          logFC=istats)

  if (plots.generate) {
    ii <- generate.GSEA.plots(x, design, contrast, out, outdir,
                              use.cache=use.cache, logFC=istats,
                              padj.threshold=plots.padj.threshold,
                              score.by=score.by, ...)
    out[, img.path := ifelse(ii$fn.exists.post, ii$fn, NA_character_)]
  }

  finished <- TRUE
  out
}

##' Creates a summary table from all the methods that were run
summarizeResults <- function(x, design, contrast, gst, results, logFC=NULL,
                             score.by=c('logFC', 't')) {
  def.take <- c('group', 'id', 'pval', 'padj', 'padj.by.group')

  ## These are the anaylsis-specific columns we want to extract from each
  ## individual result so that it is included in the outoing data.table
  add.take <- list(camera=c('Correlation', 'Direction'),
                   roast=c('PropDown', 'PropUp', 'Direction'))
                           ## 'pval.mixed', 'padj.mixed'))

  ## Scores to summarize "effect size" of the gene set.
  gs.scores <- do.geneSetScores(x, design, contrast, gst, logFC.stats=logFC,
                                score.by=score.by)

  meta.cols <- c('group', 'id', 'N', 'n') #, 'membership', 'feature.id')
  meta <- results[[1]][, meta.cols, with=FALSE]

  out <- merge(meta, gs.scores, by=c('group', 'id'))

  for (rn in names(results)) {
    take <- c(def.take, add.take[[rn]])
    res <- results[[rn]][, take, with=FALSE]
    rename <- setdiff(take, c('group', 'id'))
    if (length(rename)) {
      setnames(res, rename, paste(rename, rn, sep='.'))
    }
    out <- merge(out, res, by=c('group', 'id'))
  }

  ## Add the memberhsip and feature.id lists
  more <- results[[1]][, list(group, id, membership, feature.id)]
  m.key <- paste(more$group, more$id, sep='.')
  o.key <- paste(out$group, out$id, sep='.')
  xref <- match(o.key, m.key)
  if (any(is.na(xref))) {
    stop("Cross referencing output back to membership info failed")
  }
  if (any(duplicated(xref))) {
    stop("Duplicated keys when cross referencing membership and output")
  }
  out[, membership := more$membership[xref]]
  out[, feature.id := more$feature.id[xref]]
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
  sub.dirs <- c('images', 'data', 'cache')
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

