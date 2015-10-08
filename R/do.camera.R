##' @include validateInputs.R
NULL

validate.inputs.camera <- .validate.inputs.full.design
validate.x.camera <- validate.X

##' Runs camera on the given experiment.
##'
##' @section Gene Set Enrichment Methods:
##'
##' camera tends to be very conservative, especially when trying to estimate
##' residual correlation of gene sets from experiments with small N. For this
##' reason,  Gordon has added the ability provide a prespecified
##' \code{inter.gene.correlation} value.
##' \href{https://support.bioconductor.org/p/70005/#70195}{
##' He suggests to try a small positive number (0.05)}
do.camera <- function(gsd, x, design, contrast=ncol(design), ...) {
  stopifnot(is.conformed(gsd, x))

  args <- list(...)
  if (!is.null(args$inter.gene.cor)) {
    ## Preset inter.gene.cor values were implemented in limma v3.24.14
    if (packageVersion('limma') < '3.24.14') {
      warning("inter.gene.cor values for camera require limma >= 3.24.14, ",
              "the parameter is ignored.", immediate.=TRUE)
    }
  }

  call.args <- as.list(formals(limma::camera.default))
  for (arg in intersect(names(args), names(call.args))) {
    call.args[[arg]] <- args[[arg]]
  }

  gs.idxs <- as.expression.indexes(gsd, value='x.idx')
  call.args[['y']] <- x
  call.args[['index']] <- gs.idxs
  call.args[['design']] <- design
  call.args[['contrast']] <- contrast
  call.args[['sort']] <- FALSE
  call.args[['...']] <- NULL

  res <- do.call(camera, call.args)

  out <- cbind(geneSets(gsd)[, list(collection, name)], as.data.table(res))
  out[, NGenes := NULL]
  setnames(out, c('PValue', 'FDR'), c('pval', 'padj'))
}
