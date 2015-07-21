##' @include validateInputs.R
NULL

validate.inputs.logFC <- .validate.inputs.logFC.only
validate.x.logFC <- validate.X

##' A "pass through" function. The caller asks for "logFC" if she only wants
##' to calculate the geneSet statistics for each of the geneSets, ie. their
##' JG, mean.logFC, mean.t for each geneset.
##'
##' The user can get all the relenvant stats for this by calling either
##' \code{logFC(MultiGSEAResult} to get the individual logFCs, or
##' \code{geneSets(MultiGSEAResult} to get the geneSet level statistics like
##' the JG score.
do.logFC <- function(gsd, x, design, contrast=ncol(design), outdir=NULL,
                     use.cache=TRUE, logFC=NULL, robust.fit=FALSE,
                     robust.eBayes=FALSE, ...) {
  ## This function is actually never called.
  transform(geneSets(gsd), xlogFC=NA_real_, pval=NA_real_, padj=NA_real_)
}
