##' @include validateInputs.R
NULL

validate.inputs.fry <- .validate.inputs.full.design
validate.x.fry <- function(x) {
  ## This is not defined for DGEList yet
  validate.X(x) && !isTRUE(is(x, 'DGEList'))
}

do.fry <- function(gsd, x, design, contrast=ncol(design), ...) {
  ## This function was defined in limma v3.23.13 in April 2015
  ## ------------------------------------------------------------------------
  ## r102075 | smyth@wehi.edu.au | 2015-04-07 22:56:02 -0700 (Tue, 07 Apr 2015)
  if (packageVersion('limma') < '3.23.13') {
    stop("fry not defined until limma 3.23.13")
  }
  stopifnot(is.conformed(gsd, x))
  args <- list(...)
  call.args <- as.list(formals(limma::fry.default))
  for (arg in intersect(names(args), names(call.args))) {
    call.args[[arg]] <- args[[arg]]
  }

  gs.idxs <- as.expression.indexes(gsd, value='x.idx')
  add.dummy <- length(gs.idxs) == 1
  dummy.gs <- 'zYxWvUt'
  if (add.dummy) {
    gs.idxs[[dummy.gs]] <- sample(nrow(x), 3, replace=TRUE)
  }

  call.args[['y']] <- x
  call.args[['index']] <- gs.idxs
  call.args[['design']] <- design
  call.args[['contrast']] <- contrast
  call.args[['sort']] <- FALSE
  call.args[['...']] <- NULL
  res <- do.call(limma::fry, call.args)

  out <- cbind(geneSets(gsd)[, list(collection, name)], as.data.table(res))
  out[, NGenes := NULL]
  out <- subset(out, name != dummy.gs)

  rcols <- c(PValue='pval', FDR='padj', PValue.Mixed='pval.mixed',
             FDR.Mixed='padj.mixed')
  rcols <- rcols[names(rcols) %in% names(out)]
  setnames(out, names(rcols), rcols)
}
