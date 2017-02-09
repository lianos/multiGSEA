##' @include validateInputs.R
NULL

validate.inputs.roast <- .validate.inputs.full.design
validate.x.roast <- validate.X

##' Worker function to run roast from within a multiGSEA pipeline
##'
##' \strong{This function is not meant to be called directly, it should only be
##' called internally within \code{multiGSEA}}
do.roast <- function(gsd, x, design, contrast=ncol(design),
                     gs.idxs=as.list(gsd, active.only=TRUE, value='x.idx'), ...) {
  stopifnot(is.conformed(gsd, x))
  args <- list(...)
  call.args <- as.list(formals(limma::mroast.default))
  for (arg in intersect(names(args), names(call.args))) {
    call.args[[arg]] <- args[[arg]]
  }

  call.args[['y']] <- x
  call.args[['index']] <- gs.idxs
  call.args[['design']] <- design
  call.args[['contrast']] <- contrast
  call.args[['sort']] <- 'none'
  call.args[['...']] <- NULL
  ## earlier versions of edgeR::mroast was double-passing var.prior
  if (is(x, 'DGEList')) {
    ## var.prior and df.prior are set internally in edgeR::mroast.DGEList
    ## if we don't nuke them here, they will be passed in twice to limma::mroast
    ## and the call will error out
    call.args[['var.prior']] <- NULL
    call.args[['df.prior']] <- NULL
  }
  res <- do.call(mroast, call.args)

  out <- cbind(geneSets(gsd, .external=FALSE)[, list(collection, name)],
               as.data.table(res))
  out[, NGenes := NULL]

  setnames(out,
           c('PValue', 'FDR', 'PValue.Mixed', 'FDR.Mixed'),
           c('pval', 'padj', 'pval.mixed', 'padj.mixed'))
}
