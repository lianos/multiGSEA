##' @include validateInputs.R
NULL

validate.inputs.fry <- .validate.inputs.full.design
validate.x.fry <- function(x) {
  ## This is not defined for DGEList yet
  if (is(x, 'DGEList') && packageVersion('edgeR') < '3.14') {
    warning("fry is not implemented for a DGEList yet", immediate.=TRUE)
    return(FALSE)
  }
  validate.X(x)
}

## NOTE: Upcoming changes to fry for Bioc 3.3
## 1. edgeR::fry.DGEList method added
## 2. sort argument will be same type as roast. In bioc-3.2, fry,sort is
##    TRUE/FALSE, but roast uses character ('directional'). You will have to
##    set this to 'none'.
do.fry <- function(gsd, x, design, contrast=ncol(design),
                   gs.idxs=as.list(gsd, value='x.idx'), ...) {
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

  ## There is a bug in edgeR::fry.DGEList that double passes the standardize
  ## argument if it is passed in via `...`, so remove call.args$standardize
  ## here. The bug exists in v3.17.1, and was reported to Gordon on 11/3/2016
  if (is(x, 'DGEList')) {
    call.args[['standardize']] <- NULL
  }

  res <- do.call(limma::fry, call.args)

  out <- cbind(geneSets(gsd, .external=FALSE)[, list(collection, name)],
               as.data.table(res))
  out[, NGenes := NULL]
  out <- subset(out, name != dummy.gs)

  rcols <- c(PValue='pval', FDR='padj', PValue.Mixed='pval.mixed',
             FDR.Mixed='padj.mixed')
  rcols <- rcols[names(rcols) %in% names(out)]
  setnames(out, names(rcols), rcols)
}
