##' @include validateInputs.R
NULL

validate.inputs.romer <- .validate.inputs.full.design

## validate.x.romer <- validate.DGEList

## romer works on normal EList or DGEList. It has not been adapted
## to work with a $weights matrix (that you get from voom)
validate.x.romer <- function(x) {
  if (isTRUE(is(x, 'DGEList'))) {
    if (!is.numeric(x$common.dispersion)) {
      return("dispersion is not estimated, minimally call estimateDisp(x)")
    }
  } else if (is(x, 'EList')) {
    if (is.matrix(x$weights)) {
      ## This is coming from voom?
      warning("x has $weights. romer hasn't been implemented for voom",
              immediate.=TRUE)
    }
  } else if (!is.matrix(x)) {
    return("romer only works with DGEList, EList, or matrix for x")
  }
  return(TRUE)
}

do.romer <- function(gsd, x, design, contrast=ncol(design),
                     gs.idxs=as.list(gsd, value='x.idx'), ...) {
  stopifnot(is.conformed(gsd, x))
  args <- list(...)
  call.args <- as.list(formals(limma::romer.default))
  for (arg in intersect(names(args), names(call.args))) {
    call.args[[arg]] <- args[[arg]]
  }

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
