##' Runs several GSEA methods over a dataset for a single contrast.
##'
##' Currenty supports the following methods:
##'
##' \enumerate{
##'   \item \code{limma::camera}
##'   \item \code{limma::mroast}
##'   \item \code{limma::geneSetTest}
##' }
##'
##' This function will write several results to an output directory
##'
##' @export
##'
##' @import limma
##' @importFrom parallel mclapply
multiGSEA <- function(x, gene.sets, design=NULL, contrast=NULL,
                      methods=c('camera'), outdir=NULL,
                      plots.generate=FALSE, plots.padj.threshold=1,
                      species=NULL, mc.cores=1L, use.cache=FALSE,
                      force.reeval=FALSE, keep.outdir.onerror=TRUE, ...) {
  ## ---------------------------------------------------------------------------
  ## Argument sanity checking
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
      stop("Can't use.cacue without an `outdir`")
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
  logFC <- calculateIndividualLogFC(x, design, contrast, ...)

  out <- summarizeResults(x, design, contrast, gst, results, out.basem, logFC)

  if (plots.generate) {
    img.info <- generate.GSEA.plots(x, design, contrast, out, outdir,
                                    use.cache=use.cache, logFC=logFC,
                                    padj.threshold=plots.padj.threshold)
    if (!is.null(img.info) && nrow(img.info)) {
      xref <- match(paste(out$group, out$id, sep='.'),
                    paste(img.info$group, img.info$id, sep='.'))
      out[, img.path := img.info$fn[xref]]
    } else {
      out[, img.path := NA_character_]
    }
  }

  finished <- TRUE
  out
}

##' Creates a summary table from all the methods that were run
summarizeResults <- function(x, design, contrast, gst, results, out.base,
                             logFC=NULL) {
  def.take <- c('group', 'id', 'pval', 'padj', 'padj.by.group')

  ## These are the anaylsis-specific columns we want to extract from each
  ## individual result so that it is included in the outoing data.table
  add.take <- list(camera=c('Correlation', 'Direction'),
                   roast=c('PropDown', 'PropUp', 'Direction'))
                           ## 'pval.mixed', 'padj.mixed'))

  ## Scores to summarize "effect size" of the gene set.
  gs.scores <- do.geneSetScores(x, design, contrast, gst, logFC.stats=logFC)

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

