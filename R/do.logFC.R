#' @include validateInputs.R
NULL

validate.inputs.logFC <- .validate.inputs.logFC.only
validate.x.logFC <- validate.X

#' Worker function that is never really called within a multiGSEA pipeline
#'
#' @description
#' This is a "pass through" function. The caller asks for \code{method="logFC"}
#' if she only wants to calculate the geneSet statistics for each of the
#' geneSets, ie. mean.logFC, mean.t for each geneset.
#'
#' **This function is not meant to be called directly.** It should only be
#' called internally within [multiGSEA()].
#'
#' @noRd
do.logFC <- function(gsd, x, design, contrast=ncol(design), outdir=NULL,
                     use.cache=TRUE, logFC=NULL, robust.fit=FALSE,
                     robust.eBayes=FALSE, ...) {
  ## This function is actually never called.
  res <- transform(geneSets(gsd, as.dt=TRUE),
                   logFC=NA_real_, pval=NA_real_, padj=NA_real_)
  setattr(res, 'rawresult', TRUE)
}

mgres.logFC <- function(res, gsd, ...) res
