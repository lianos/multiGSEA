##' @include validateInputs.R
NULL

validate.inputs.svdGeneSetTest <- .validate.inputs.full.design
validate.x.svdGeneSetTest <- validate.X

##' Transforms gene x sample to geneset x sample and run differential expression
##'
do.svdGeneSetTest <- function(gsd, x, design, contrast=ncol(design),
                              gs.idxs=as.expression.indexes(gsd, value='x.idx'),
                              trend.eBayes=FALSE, ...) {
  stopifnot(is.conformed(gsd, x))
  X <- scoreSingleSamples(gsd, x, 'svd', scale=TRUE, center=TRUE,
                          unscale=TRUE, uncenter=TRUE,
                          as.matrix=TRUE)
  stopifnot(all(colnames(x) == colnames(X)))

  if (is(x, 'DGEList')) {
    X <- cpm(x, prior.count=5, log=TRUE)
    trend.eBayes <- TRUE
  } else if (is(x, 'EList') && is.matrix(x$weights)) {
    trend.eBayes <- TRUE
  }

  res <- calculateIndividualLogFC(X, design, contrast, trend.eBayes=trend.eBayes, ...)

  out <- cbind(geneSets(gsd, .external=FALSE)[, list(collection, name)],
               as.data.table(res))
  out[, NGenes := NULL]
}

