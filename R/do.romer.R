##' @include validateInputs.R
NULL

validate.inputs.romer <- .validate.inputs.full.design
validate.x.romer <- validate.DGEList

do.romer <- function(gsd, x, design, contrast=ncol(design), ...) {
  stopifnot(is.conformed(gsd, x))
  args <- list(...)
  call.args <- as.list(formals(limma::romer.default))
  for (arg in intersect(names(args), names(call.args))) {
    call.args[[arg]] <- args[[arg]]
  }

  gs.idxs <- as.expression.indexes(gsd, value='x.idx')
  call.args[['y']] <- x
  call.args[['index']] <- gs.idxs
  call.args[['design']] <- design
  call.args[['contrast']] <- contrast
  call.args[['...']] <- NULL

  res <- do.call(romer, call.args)

  out <- cbind(geneSets(gsd, .external=FALSE)[, list(collection, name)],
               as.data.table(res))
  out[, NGenes := NULL]

  setnames(out,
           c('Up', 'Down', 'Mixed'),
           c('pval.up', 'pval.down', 'pval'))

  out[, padj := p.adjust(pval)]
  out[, padj.up := p.adjust(pval.up)]
  out[, padj.down := p.adjust(pval.down)]

  out
}

