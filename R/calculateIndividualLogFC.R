##' Runs a vanilla limma differential expression analysis
##'
##' @export
##' @importFrom data.table setnames
##'
##' @param x The expression object. This can be 1 column matrix if you are not
##' running any analysis, and this function essentially is just a "pass through"
##' @param design The design matrix for the experiment
##' @param contrast The contrast you want to test and provide stats for
##' @param provide Specifies what you want to get out of this, just the
##'   \code{logFC} or \code{t} statistics, or the whole result \code{table}.
##'   TODO: remove this parameter and just always return the table, if one is
##'   generated as a result of the function call (ie. \code{ncol(x) > 1})
##' @param robust.fit If \code{TRUE}, sets \code{method='robust'} in the
##'   \code{lmFit} call. Defaults to \code{FALSE}
##' @param robust.eBayes The value of the \code{robust} parameter in the
##'   \code{eBayes} call.
##' @param ... parameters passed through to the \code{lmFit} call.
##'
##' @return If \code{provde == 'table'}: a data.frame of results that resembles
##'   the output of limma::topTable with P.Value and adj.P.Val columns renamed
##'   to pval, padj. Even if the \code{x} input is a vector, we reconstruct a
##'   "dummy" \code{data.frame}. Otherwise, the column vector that corresponds
##'   to \code{provide}.
calculateIndividualLogFC <- function(x, design, contrast,
                                     provide=c('table', 'logFC', 't'),
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
    setnames(tt, c('P.Value', 'adj.P.Val'), c('pval', 'padj'))
    out <- tt
  } else {
    out <- data.frame(logFC=x[, 1L], t=x[, 1L], pval=NA_real_, padj=NA_real_)
    rownames(out) <- rownames(x)
    fit <- NULL
  }

  if (provide != 'table') {
    out <- setNames(out[[provide]], rownames(out))
  }

  if (with.fit) list(result=out, fit=fit) else out
}
