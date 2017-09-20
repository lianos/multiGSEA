##' @include validateInputs.R
NULL

validate.inputs.svdGeneSetTest <- .validate.inputs.full.design
validate.x.svdGeneSetTest <- validate.X

##' @title
##' Worker function summarizes gene set activity by SVD score and t-test across
##' contrast
##'
##' @description
##' The idea is to summarize gene set level activity per sample by their
##' svd (GSDecon) score to create a geneset x sample matrix. Then uses this
##' matrix in a limma based framework to assess their "differential expression"
##'
##' This method seems like it should produce sane results, but I haven't really
##' tested it yet and is therefore experimental. The fact that the
##' "svdGeneSetTest" method isn't mentioned anywhere in the documentation yet is
##' by design. Still, you have found that it exists and you can try to use it
##' if you like (or, please contact me to tell me why it's a bad idea!)
##'
##' \strong{This function is not meant to be called directly, it should only be
##' called internally within \code{multiGSEA}}
do.svdGeneSetTest <- function(gsd, x, design, contrast=ncol(design),
                              gs.idxs=as.list(gsd, active.only=TRUE, value='x.idx'),
                              trend.eBayes=FALSE, ...) {
  stop("TODO: Implement svdGeneSetTest")
  stopifnot(is.conformed(gsd, x))
  X <- scoreSingleSamples(gsd, x, 'ewm', scale=TRUE, center=TRUE,
                          unscale=TRUE, uncenter=TRUE,
                          as.matrix=TRUE)
  stopifnot(all(colnames(x) == colnames(X)))

  if (missing(trend.eBayes)) {
    if (is(x, 'DGEList')) {
      X <- cpm(x, prior.count=5, log=TRUE)
      trend.eBayes <- TRUE
    } else if (is(x, 'EList') && is.matrix(x$weights)) {
      trend.eBayes <- TRUE
    }
  }

  res <- calculateIndividualLogFC(X, design, contrast,
                                  trend.eBayes=trend.eBayes, ...)

  out <- cbind(geneSets(gsd, as.dt=TRUE)[, list(collection, name)], setDT(res))
}

