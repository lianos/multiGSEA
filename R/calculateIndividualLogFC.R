##' Runs a vanilla limma differential expression analysis
##'
##' Note that differential expression analysis is always run through limma,
##' so if \code{x} is a \code{DGEList}, it will be voomd first and then
##' processed "as usual".
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
##' @return A data.table with the logFC results per row of \code{x}. The
##'   \code{$featureId} column is the rownames (ids) of x.
calculateIndividualLogFC <- function(x, design, contrast=ncol(design),
                                     robust.fit=FALSE, robust.eBayes=FALSE,
                                     with.fit=FALSE, ...) {
  if (is(x, 'DGEList')) {
    message("vooming DGEList to calculate logFCs")
    x <- voom(x, design, plot=FALSE)
  }
  if (ncol(x) > 1) {
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
    tt <- transform(as.data.table(tt), featureId=rownames(x))
    setnames(tt, c('P.Value', 'adj.P.Val'), c('pval', 'padj'))
    out <- tt
  } else {
    out <- data.table(logFC=x[, 1L], AveExpr=NA_real_, t=x[, 1L],
                      pval=NA_real_, padj=NA_real_,
                      featureId=rownames(x))
    fit <- NULL
  }
  out[, x.idx := 1:nrow(x)]
  if ('ID' %in% names(out)) {
    out[, ID := NULL]
  }
  if (with.fit) list(result=out, fit=fit) else out
}

##' Checks that a provided table is "similar enough" the the result generated
##' from calculateIndividualLogFC
##'
##' @param logFC the table to check
##' @param x An expression-like object to further test against.
is.logFC.like <- function(logFC, x, as.error=FALSE) {
  ref.dt <- data.table(logFC=numeric(), t=numeric(), pval=numeric(),
                       padj=numeric(), featureId=character())
  ref.check <- check.dt(logFC, ref.dt)
  if (isTRUE(ref.check)) {
    if (!missing(x)) {
      if (nrow(logFC) != nrow(x)) {
        ref.check <- "nrow(logFC) != nrow(x)"
      } else {
        missed.features <- setdiff(rownames(x), logFC$featureId)
        if (length(missed.features)) {
          ref.check <- sprintf("%d features missing from featureId",
                               length(missed.features))
        }
      }
    }
  }

  if (as.error && is.character(ref.check)) {
    stop("Provided logFC is not a valid table:\n  ",
         paste(ref.check, collapse="\n  "))
  }

  isTRUE(ref.check)
}
