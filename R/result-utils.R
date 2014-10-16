## This functions work over the large data.table that gets returned from the
## multiGSEA call.

##' Assembles a matrix of nominal or adjusted pvalues from a multiGSEA result
##'
##' @param x The output from multiGSEA
##' @param pval Are we testing pvalues or adjusted pvalues?
##'
##' @return A matrix of all of the
p.matrix <- function(x, pval=c('adjusted', 'nominal')) {
  if (!is.data.table(x)) {
    stop("The input should (minimally(!)) be a data.table")
  }
  pval <- switch(match.arg(pval), nominal='pval', adjusted='padj')
  regex <- sprintf('^%s\\.', pval)
  as.matrix(x[, grep(regex, names(x)), with=FALSE])
}

##' Identifies which geneset has a significant result in >= 1 GSEAs
##'
##' @export
##'
##' @param x The output from multiGSEA
##' @param pval Are we testing pvalues or adjusted pvalues?
##' @param threshold The threshold that \code{pval} has to be less than
##' @param n.sig The number of pval|padj columns that need to hit significance
##'   for the geneset to be flagged as such
##'
##' @return a logical vector as long as \code{nrow(x)} with \code{TRUE}
##'   indicates the gene sets that are significant given the criteria
significantGeneSets <- function(x, pval=c('adjusted', 'nominal'),
                                threshold=0.10, n.sig=1L) {
  rowSums(p.matrix(x, pval) <= threshold) >= n.sig
}
