##' Runs a vanilla limma differential expression analysis
##'
##' @export
##'
##' @param x The expression object. This can be 1 column matrix if you are not
##' running any analysis, and this function essentially is just a "pass through"
##' @param design The design matrix for the experiment
##' @param contrast The contrast you want to test and provide stats for
##' @param provide Specifies what you want to get out of this, just the
##'   \code{logFC} or \code{t} statistics, or the whole \code{topTable}.
##' @param robust.fit If \code{TRUE}, sets \code{method='robust'} in the
##'   \code{lmFit} call. Defaults to \code{FALSE}
##' @param robust.eBayes The value of the \code{robust} parameter in the
##'   \code{eBayes} call.
##' @param ... parameters passed through to the \code{lmFit} call.
##'
##' @return A vector of \code{logFC} or \code{t} statistics, or the entire
##'   \code{topTable} \code{data.frame}.
calculateIndividualLogFC <- function(x, design, contrast,
                                     provide=c('logFC', 't', 'topTable'),
                                     robust.fit=FALSE, robust.eBayes=FALSE,
                                     with.fit=FALSE, ...) {
  if (ncol(x) > 1) {
    provide <- match.arg(provide)
    ## Do limma fits on this thing and let it rip
    fit <- lmFit(x, design, method=if (robust.fit) 'robust' else 'ls', ...)
    if (length(contrast) > 1) {
      if (length(contrast) != ncol(design)) {
        stop("Invalid contrast vector, must be as long as columns in design")
      }
      fit <- contrasts.fit(fit, contrast)
      contrast <- 1L
    }
    fit  <- eBayes(fit, robust=robust.eBayes)
    tt <- topTable(fit, contrast, number=Inf, sort.by='none')
    if (provide == 'topTable') {
      out <- tt
    } else {
      out <- tt[[provide]]
      names(out) <- rownames(x)
    }
  } else {
    out <- x[, 1L]
    names(out) <- rownames(x)
    fit <- NULL
  }

  if (with.fit) list(result=out, fit=fit) else out
}
