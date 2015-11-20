##' @include validateInputs.R
NULL

validate.inputs.roast <- .validate.inputs.full.design
validate.x.roast <- validate.X

do.roast <- function(gsd, x, design, contrast=ncol(design), ...) {
  stopifnot(is.conformed(gsd, x))
  args <- list(...)
  call.args <- as.list(formals(limma::mroast.default))
  for (arg in intersect(names(args), names(call.args))) {
    call.args[[arg]] <- args[[arg]]
  }

  gs.idxs <- as.expression.indexes(gsd, value='x.idx')
  call.args[['y']] <- x
  call.args[['index']] <- gs.idxs
  call.args[['design']] <- design
  call.args[['contrast']] <- contrast
  call.args[['sort']] <- 'none'
  call.args[['...']] <- NULL
  ## earlier versions of edgeR::mroast was double-passing var.prior
  if (packageVersion('edgeR') < '3.11') {
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
