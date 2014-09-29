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
                      methods=c('camera', 'gst'), outdir=NULL,
                      plots.generate=TRUE, plots.padj.threshold=0.3,
                      cleanup.on.fail=TRUE, species=NULL,
                      mc.cores=1L, ...) {
  species <- match.species(species)
  .unsupportedGSEAmethods(methods)
  outdir.created <- .initOutDir(outdir)
  finished <- FALSE

  on.exit({
    if (!finished && outdir.created && dir.exists(outdir)) {
      warning("An error in `mulitGSEA stopped it from finishing ...",
              immediate.=TRUE)
      unlink(outdir, recursive=TRUE)
    }
  })

  args <- list(x=x, gene.sets=gene.sets, design=design, contrast=contrast)
  req.x.rownames <- !inherits(gene.sets, 'GeneSetTable')

  inputs <- validateInputs(x, design, contrast, methods,
                           require.x.rownames=req.x.rownames)
  x <- inputs$x
  design <- inputs$design
  contrast <- inputs$contrast

  if (!is(gene.sets, 'GeneSetTable')) {
    gst <- GeneSetTable(x, gene.sets)
  } else {
    gst <- conform(gene.sets, x)
  }

  results <- mclapply(methods, function(method) {
    fn <- getFunction(paste0('do.', method))
    fn(x, gst, design, contrast, ...)
  }, mc.cores=mc.cores)
  names(results) <- methods

  gs.scores <- do.geneSetScores(x, gst, design, contrast)

  meta.cols <- c('group', 'id', 'N', 'n') #, 'membership', 'feature.id')
  meta <- results[[1]][, meta.cols, with=FALSE]

  out <- merge(meta, gs.scores, by=c('group', 'id'))

  def.take <- c('group', 'id', 'pval', 'padj', 'padj.by.group')
  add.take <- list(camera=c('Correlation', 'Direction'))

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

  finished <- TRUE
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
  sub.dirs <- c('images', 'data')
  for (sd in sub.dirs) {
    dir.create(file.path(outdir, sd))
  }

  outdir.created
}

##' Runs several GSEA methods over a dataset.
##'
##' @param result.name The name/directory used to store the results for this
##' (unique) analysis.
rmd.multi.GSEA <- function(result.name, x, gene.sets, design=NULL,
                           contrast=NULL, ...) {

}

