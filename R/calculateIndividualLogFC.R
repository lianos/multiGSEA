##' Utility function to run limma differential expression analysis
##'
##' @details
##' This function fits linear modles (or glms) to perform differential
##' expression analyses. If the \code{x} object is a \code{DGEList} the
##' analysis will be performed using edgeR's quasi-likelihood framework,
##' otherwise limma will be used for all other scenarios.
##'
##' If \code{x} is a \code{DGEList} we require that \code{estimateDisp} has
##' already been called. If you prefer to analyze rnaseq data using voom,
##' be sure that \code{x} is the object that has been returned from a
##' call to \code{\link[limma]{voom}} (or
##' \code{\link[limma]{voomWithQualityWeights}}.
##'
##' The documentation here is speaking the language of a "limma" analysis,
##' however for each parameter, there is an analagous function/parameter that
##' will be delegated to.
##'
##' Lastly, if \code{x} is simply a single column matrix, we assume that we are
##' just passing a single pre-ranked vector of statistics through multiGSEA's
##' analysis pipelines (for use in methods like "fgsea", "cameraPR", etc.), and
##' a logFC-like data.frame is constructed with these statistics in the
##' \code{logFC} and \code{t} columns.
##'
##' @export
##'
##' @param x The expression object. This can be 1 column matrix if you are not
##' running any analysis, and this function essentially is just a "pass through"
##' @param design The design matrix for the experiment
##' @param contrast The contrast you want to test and provide stats for. By
##'   default this tests the last column of the \code{design} matrix. If you
##'   want to test a custom contrast, this can be a contrast vector, which
##'   means that it should be as long as \code{ncol(design)} and it most-often
##'   sum to one. In the future, the user will be able to specify a range of
##'   coefficients over \code{design} to perform an ANOVA analysis.
##' @param robust.fit The value of the \code{robust} parameter to pass down to
##'   the \code{\link[limma]{lmFit}} function. Defaults to \code{FALSE}.
##' @param robust.eBayes The value of the \code{robust} parameter to pass down
##'   to the \code{\link[limma]{eBayes}} function.
##' @param trend.eBayes The value of the \code{trend} parameter to pass down to
##'   the \code{\link[limma]{eBayes}} function.
##' @param treat.lfc If this is numeric, this activates limma's "treat"
##'   functionality and tests for differential expression against this
##'   specified log fold change threshold. This defaults to \code{NULL}.
##' @param confint add confidence intervals to \code{topTable} output (default
##'   \code{TRUE})? Ignored if \code{x} is a \code{DGEList}.
##' @param with.fit If \code{TRUE}, this function returns a list object with
##'   both the fit and the table of logFC statistics, otherwise just the
##'   logFC statistics table is returned.
##' @param use.qlf If \code{TRUE} (default), will use edgeR's quasilikelihood
##'   framework for analysis, otherwise uses glmFit/glmTest.
##' @param ... parameters passed down into the relevant limma/edgeR based
##'   functions.
##' @param as.dt Return the result as a \code{data.table}? Defaults to
##'   \code{FALSE}.
##' @return If \code{with.fit == FALSE} (the default) a \code{data.table} of
##'   logFC statistics for the contrast under test. Otherwise, a list is
##'   returned with \code{$result} containing the logFC statistics, and
##'   \code{$fit} has the limma fit for the data/design/contrast under test.
calculateIndividualLogFC <- function(x, design, contrast=ncol(design),
                                     robust.fit=FALSE, robust.eBayes=FALSE,
                                     trend.eBayes=FALSE, treat.lfc=NULL,
                                     confint=TRUE, with.fit=FALSE,
                                     use.qlf = TRUE, ..., as.dt=FALSE) {
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
    if (use.qlf) {
      fit <- glmQLFit(x, design, robust=TRUE)
    } else {
      fit <- glmFit(x, design)
    }
    if (use.treat) {
      if (do.contrast) {
        tt <- glmTreat(fit, contrast=contrast, lfc=treat.lfc)
      } else {
        tt <- glmTreat(fit, coef=contrast, lfc=treat.lfc)
      }
    } else {
      etest <- if (use.qlf) glmQLFTest else glmLRT
      if (do.contrast) {
        tt <- etest(fit, contrast = contrast)
      } else {
        tt <- etest(fit, coef=contrast)
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

  if (!as.dt) setDF(out)
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
