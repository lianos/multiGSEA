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
##' importFrom parallel mclapply
multiGSEA <- function(x, gene.sets, design=NULL, contrast=NULL,
                      methods=c('camera', 'roast', 'gst'), outdir=NULL,
                      plots.generate=TRUE, plots.padj.threshold=0.3,
                      cleanup.on.fail=TRUE, species=.wehi.msigdb.species,
                      mc.cores=1L, ...) {
  species <- match.arg(species)
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

  gs.table <- preprocessGeneSets(gene.sets, x, species)

  results <- mclapply(methods, function(method) {
    fn <- getFunction(paste0('do.', method))
    fn(x, gs.table, design, contrast, ...)
  }, mc.cores=mc.cores)

  finished <- TRUE
  out
}

.initOutDir <- function(outdir) {
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

