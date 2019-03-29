#' Utility function to run limma differential expression analysis
#'
#' @details
#' This function fits linear modles (or glms) to perform differential
#' expression analyses. If the `x` object is a `DGEList` the
#' analysis will be performed using edgeR's quasi-likelihood framework,
#' otherwise limma will be used for all other scenarios.
#'
#' If `x` is a [edgeR::DGEList()] we require that [edgeR::estimateDisp()] has
#' already been called. If you prefer to analyze rnaseq data using voom, be sure
#' that `x` is the object that has been returned from a call to [limma::voom()]
#' (or [limma::voomWithQualityWeights()].
#'
#' The documentation here is speaking the language of a "limma" analysis,
#' however for each parameter, there is an analagous function/parameter that
#' will be delegated to.
#'
#' Lastly, if `x` is simply a single column matrix, we assume that we are
#' just passing a single pre-ranked vector of statistics through multiGSEA's
#' analysis pipelines (for use in methods like "fgsea", "cameraPR", etc.), and
#' a logFC-like data.frame is constructed with these statistics in the
#' `logFC` and `t` columns.
#'
#' @export
#'
#' @param x The expression object. This can be 1 column matrix if you are not
#'   running any analysis, and this function essentially is just a
#'   "pass through"
#' @param design The design matrix for the experiment
#' @param contrast The contrast you want to test and provide stats for. By
#'   default this tests the last column of the `design` matrix. If you
#'   want to test a custom contrast, this can be a contrast vector, which
#'   means that it should be as long as `ncol(design)` and it most-often sum to
#'   one. In the future, the user will be able to specify a range of
#'   coefficients over `design` to perform an ANOVA analysis, cf.
#'   Issue #11 (https://github.com/lianos/multiGSEA/issues/11).
#' @param robust.fit The value of the `robust` parameter to pass down to the
#'   [limma::lmFit()] function. Defaults to `FALSE`.
#' @param robust.eBayes The value of the `robust` parameter to pass down to
#'   the limma::eBayes()] function.
#' @param trend.eBayes The value of the `trend` parameter to pass down to the
#'   [limma::eBayes()] function.
#' @param treat.lfc If this is numeric, this activates limma's "treat"
#'   functionality and tests for differential expression against this
#'   specified log fold change threshold. This defaults to `NULL`.
#' @param confint add confidence intervals to `topTable` output (default
#'   `TRUE`)? Ignored if `x` is a `DGEList`.
#' @param with.fit If `TRUE`, this function returns a list object with
#'   both the fit and the table of logFC statistics, otherwise just the
#'   logFC statistics table is returned.
#' @param use.qlf If `TRUE` (default), will use edgeR's quasilikelihood
#'   framework for analysis, otherwise uses glmFit/glmTest.
#' @param ... parameters passed down into the relevant limma/edgeR based
#'   functions.
#' @template asdt-param
#' @return If `with.fit == FALSE` (the default) a `data.table` of
#'   logFC statistics for the contrast under test. Otherwise, a list is
#'   returned with `$result` containing the logFC statistics, and
#'   `$fit` has the limma fit for the data/design/contrast under test.
calculateIndividualLogFC <- function(x, design, contrast=ncol(design),
                                     robust.fit=FALSE, robust.eBayes=FALSE,
                                     trend.eBayes=FALSE, treat.lfc=NULL,
                                     confint=TRUE, with.fit=FALSE,
                                     use.qlf = TRUE, ..., as.dt=FALSE) {
  do.contrast <- !is.vector(x) &&
    ncol(x) > 1L &&
    !is.null(design) &&
    length(contrast) == ncol(design)

  if (is.vector(x)) {
    # this is a preranked vector
    x <- matrix(x, ncol=1L, dimnames=list(names(x), NULL))
  }

  if (ncol(x) == 1L) {
    # The user passed in a vector of statistics. Only a very few number of
    # GSEA methods support this, but they are useful (fgsea and cameraPR). We
    # make the output of this look like a "DGE" data.table
    out <- data.table(logFC=x[, 1L], AveExpr=NA_real_, t=x[, 1L],
                      pval=NA_real_, padj=NA_real_, confint=NA_real_,
                      featureId=rownames(x))
    fit <- NULL
    test_type <- "preranked"
  } else if (do.contrast) {
    test_type <- "ttest"
    assert_matrix(design, nrows = ncol(x))
    assert_numeric(contrast, len = ncol(design))
  } else {
    # Testing a coefficient, or multiple coefficients (anova)
    assert_matrix(design, nrows = ncol(x))
    if (is.character(contrast)) {
      assert_choice(contrast, colnames(design))
      contrast <- match(contrast, colnames(design))
    }
    assert_integerish(contrast, lower = 1, upper = ncol(design),
                      min.len = 1L, max.len = ncol(design) - 1L)
    test_type <- if (length(contrast) == 1L) "ttest" else "anova"
  }

  use.treat <- FALSE
  if (is.numeric(treat.lfc)) {
    stopifnot(length(treat.lfc) == 1L, treat.lfc > 0)
    use.treat <- TRUE
  }

  if (is(x, 'DGEList')) {
    # We default to using the quasi-likelihood piepline with edgeR with
    # robust fitting. Setting robust=TRUE here should be roughly equivalent
    # to eBayes(..., robust=TRUE) in the limma world.
    if (!disp.estimated(x)) {
      stop("Dispersions not estimated, need to run estimateDisp first")
    }
    if (use.qlf) {
      fit <- glmQLFit(x, design, robust=TRUE)
    } else {
      fit <- glmFit(x, design)
    }
    if (use.treat && test_type == "ttest") {
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
    # If x is matrix-like but not a DGEList, we assume you are OK to run the
    # limma pipeline.
    fit <- lmFit(x, design, method = if (robust.fit) 'robust' else 'ls', ...)
    if (do.contrast) {
      fit <- contrasts.fit(fit, contrast)
      contrast <- 1L
    }
    if (use.treat && test_type == "ttest") {
      fit <- treat(fit, lfc = treat.lfc, robust = robust.eBayes,
                   trend = trend.eBayes)
      tt <- topTreat(fit, contrast, number = Inf, sort.by = "none",
                     confint = confint)
    } else {
      fit <- eBayes(fit, robust = robust.eBayes, trend = trend.eBayes)
      tt <- topTable(fit, contrast, number = Inf, sort.by = "none",
                     confint = confint)
    }
    tt <- transform(setDT(tt), featureId = rownames(x))
    setnames(tt, c("P.Value", "adj.P.Val"), c("pval", "padj"))
    out <- tt
  }

  out[, x.idx := 1:nrow(x)]
  if ('ID' %in% names(out)) {
    out[, ID := NULL]
  }

  if (!as.dt) out <- setDF(out)
  if (with.fit) out <- list(result=out, fit=fit)
  setattr(out, "test_type", test_type)
  out
}

#' Checks that a provided table is "similar enough" the the result generated
#' from calculateIndividualLogFC
#'
#' @noRd
#'
#' @param logFC the table to check
#' @param x An expression-like object to further test against.
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
