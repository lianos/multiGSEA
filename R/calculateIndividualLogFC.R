##' Utility function to run limma differential expression analysis
##'
##' @details
##' If \code{x} is a \code{DGEList} we require that \code{estimateDisp} has
##' already been called and this will run a quasi-likelihood based differential
##' expression analysis, otherwise a "vanilla" limma run will be run.
##'
##' Lastly, if \code{x} is simply a single column matrix, we assume that we are
##' just passing a single pre-ranked vector of statistics through multiGSEA's
##' analysis pipelines, and a logFC-like data.frame is constructed with these
##' statistics in the \code{logFC} and \code{t} columns.
##'
##' @export
##'
##' @param x The expression object. This can be 1 column matrix if you are not
##' running any analysis, and this function essentially is just a "pass through"
##' @param design The design matrix for the experiment
##' @param contrast The contrast you want to test and provide stats for
##' @param robust.fit If \code{TRUE}, sets \code{method='robust'} in the
##'   \code{lmFit} call. Defaults to \code{FALSE}
##' @param robust.eBayes The value of the \code{robust} parameter in the
##'   \code{eBayes} call.
##' @param treat.lfc If this is numeric, this activates limma's "treat"
##'   functionality and tests for differential expression against this
##'   specified log fold change threshold. This defaults to \code{NULL}.
##' @param confint add confidence intervals to \code{topTable} output (default
##'   \code{TRUE})? Ignored if \code{x} is a \code{DGEList}.
##' @param with.fit If \code{TRUE}, this function returns the fit in addition
##'   to the logFC statistics.
##' @param ... parameters passed through to the \code{lmFit} call.
##' @param .external If called w/in multiGSEA package, we don't check to convert
##'   outgoing data.table to user's preferred type. \code{.external} is by
##'   default \code{TRUE}, as if someone on the outside is calling.
##' @return If \code{with.fit == FALSE} (the default) a \code{data.table} of
##'   logFC statistics for the contrast under test. Otherwise, a list is
##'   returned with \code{$result} containing the logFC statistics, and
##'   \code{$fit} has the limma fit for the data/design/contrast under test.
calculateIndividualLogFC <- function(x, design, contrast=ncol(design),
                                     robust.fit=FALSE, robust.eBayes=FALSE,
                                     trend.eBayes=FALSE, treat.lfc=NULL,
                                     confint=TRUE, with.fit=FALSE, ...,
                                     .external=TRUE) {
  do.contrast <- !is.vector(x) &&
    ncol(x) > 1L &&
    !is.null(design) &&
    length(contrast) > 1L
  if (do.contrast) {
    if (length(contrast) != ncol(design)) {
      stop("Invalid contrast vector, must be as long as columns in design")
    }
  } else if (!is.null(design) && !is.null(contrast) &&
             length(contrast) != 1 &&
             contrast > 0 && contrast <= ncol(design)) {
    stop("Illegal coefficient to test in design")
  }

  use.treat <- FALSE
  if (is.numeric(treat.lfc)) {
    stopifnot(length(treat.lfc) == 1L, treat.lfc > 0)
    use.treat <- TRUE
  }

  if (is(x, 'DGEList')) {
    ## We default to using the quasi-likelihood piepline with edgeR with
    ## robust fitting. Setting robust=TRUE here should be roughly equivalent
    ## to eBayes(..., robust=TRUE) in the limma world.
    if (!disp.estimated(x)) {
      stop("Dispersions not estimated, need to run estimateDisp first")
    }
    fit <- glmQLFit(x, design, robust=TRUE)
    if (use.treat) {
      if (do.contrast) {
        tt <- glmTreat(fit, contrast=contrast, lfc=treat.lfc)
      } else {
        tt <- glmTreat(fit, coef=contrast, lfc=treat.lfc)
      }
    } else {
      if (do.contrast) {
        tt <- glmQLFTest(fit, contrast=contrast)
      } else {
        tt <- glmQLFTest(fit, coef=contrast)
      }
    }
    tt <- as.data.frame(topTags(tt, Inf, sort.by='none'))
    tt <- transform(setDT(tt), t=NA_real_, featureId=rownames(x))
    setnames(tt, c('logCPM', 'PValue', 'FDR'), c('AveExpr', 'pval', 'padj'))
    out <- tt
  } else if (ncol(x) > 1L) {
    ## If x is matrix-like but not a DGEList, we assume you are OK to run the
    ## limma pipeline.
    fit <- lmFit(x, design, method=if (robust.fit) 'robust' else 'ls', ...)
    if (do.contrast) {
      fit <- contrasts.fit(fit, contrast)
      contrast <- 1L
    }
    if (use.treat) {
      fit <- treat(fit, lfc=treat.lfc, robust=robust.eBayes, trend=trend.eBayes)
      tt <- topTreat(fit, contrast, number=Inf, sort.by='none', confint=confint)
    } else {
      fit <- eBayes(fit, robust=robust.eBayes, trend=trend.eBayes)
      tt <- topTable(fit, contrast, number=Inf, sort.by='none', confint=confint)
    }
    tt <- transform(setDT(tt), featureId=rownames(x))
    setnames(tt, c('P.Value', 'adj.P.Val'), c('pval', 'padj'))
    out <- tt
  } else {
    ## The user passed in a vector of statistics. Only a very few number of
    ## GSEA methods support this, but ...
    out <- data.table(logFC=x[, 1L], AveExpr=NA_real_, t=x[, 1L],
                      pval=NA_real_, padj=NA_real_, confint=NA_real_,
                      featureId=rownames(x))
    fit <- NULL
  }

  # x.idx <- ID <- NULL # silence R CMD check NOTEs
  out[, x.idx := 1:nrow(x)]
  if ('ID' %in% names(out)) {
    out[, ID := NULL]
  }

  out <- ret.df(out, .external=.external)
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
