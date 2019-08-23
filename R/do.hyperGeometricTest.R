#' @include validateInputs.R
NULL

validate.x.hyperGeometricTest <- validate.X
validate.inputs.hyperGeometricTest <- function(x, design, contrast,
                                               xmeta. = NULL, ...) {
  if (!is.data.frame(xmeta.)) {
    default <- .validate.inputs.full.design(x, design, contrast)
    if (length(default)) {
      return(default)
    }
  }
  errs <- list()
  return(errs)
}

#' Performs a hypergeometric test for significance of gene set membership.
#'
#' Genes are selected for testing against each geneset by virture of them
#' passing a maximum FDR and minimum log fold change as perscribed by the
#' `min.logFC` and `max.padj` parameters, respectfully.
#'
#' Note that we are intentionally adding a hyperG.selected column by reference
#' so that this information is kicked back to the caller multiGSEA function
#' and included in downstream reporting.
#'
#' @param gsd The [GeneSetDb()] for analysis
#' @param x The expression object
#' @param design Experimental design
#' @param contrast The contrast to test
#' @param direction Same as direction in `GOstats`
#' @param logFC The logFC data.table from `calculateIndividualLogFC`
#' @param ... arguments to pass down into `calculateIndividualLogFC`
do.hyperGeometricTest <- function(gsd, x, design, contrast=ncol(design),
                                  direction=c('over', 'under'),
                                  use.treat=FALSE,
                                  feature.min.logFC=1, feature.max.padj=0.10,
                                  logFC=NULL, ...) {
  stopifnot(is.conformed(gsd, x))
  direction <- match.arg(direction)

  if (is.null(logFC)) {
    treat.lfc <- if (use.treat) feature.min.logFC else NULL
    logFC <- calculateIndividualLogFC(x, design, contrast, treat.lfc=treat.lfc,
                                      ..., as.dt=TRUE)
  }
  is.logFC.like(logFC, x, as.error=TRUE)
  logFC <- setDT(copy(logFC))
  if (is.null(logFC$significant)) {
    logFC[, significant := {
      padj <= feature.max.padj & abs(logFC) >= feature.min.logFC
    }]
  }

  drawn <- logFC[significant == TRUE]$featureId
  out <- hyperGeometricTest(gsd, drawn, rownames(x), direction, do.conform=FALSE,
                            as.dt=TRUE)
  setattr(out, 'rawresult', TRUE)
}

mgres.hyperGeometricTest <- function(res, gsd, ...) res

#' Perform Hypergeometric Enrichment tests across a GeneSetDb.
#'
#' @export
#'
#' @param gsd The [GeneSetDb()] object to run tests against
#' @param selected The ids of the selected features
#' @param universe The ids of the universe
#' @param direction Same as direction in `GOstats`
#' @param do.conform By default `TRUE`: does some gymnastics to conform
#'   the `gsd` to the `universe` vector. This should neber be set to `FALSE`,
#'   but this parameter is here so that when this function is called from the
#'   [multiGSEA()] codepath, we do not have to reconform the [GeneSetDb()]
#'   object, because it has already been done.
#' @template asdt-param
#' @return A `data.frame` of results
hyperGeometricTest <- function(gsd, selected, universe,
                               direction=c('over', 'under'),
                               do.conform=TRUE, as.dt=FALSE) {
  stopifnot(is(gsd, 'GeneSetDb'))
  stopifnot(is.character(selected))
  stopifnot(is.character(universe))
  selected <- unique(selected)
  universe <- unique(universe)

  ## This needs to be conformed to work
  if (!do.conform) {
    gsd <- conform(gsd, universe)
  }
  direction <- match.arg(direction)
  dir.over <- direction == 'over'

  selected <- intersect(selected, universe)
  numDrawn <- length(selected)

  # silence R CMD check NOTEs
  odds <- expected <- NULL
  if (numDrawn == 0) {
    warning("No selected genes in hyperGeometricTest", immediate.=TRUE)
    out <- geneSets(gsd, as.dt=TRUE)[, list(collection, id)]
    out[, pval := NA_real_]
    out[, odds := NA_real_]
    out[, expected := NA_real_]
    out[, padj := NA_real_]
    out[, padj.by.collection := NA_real_]
  } else {
    numB <- length(universe) ## number of genes in universe
    out <- geneSets(gsd, as.dt=TRUE)[, {
      ids <- featureIds(gsd, .BY[[1]], .BY[[2]])
      numW <- length(ids)
      Wdrawn <- intersect(ids, selected)
      numWdrawn <- length(Wdrawn)
      ## numW: number of genes in GO category
      ## numB: size of universe
      ## numDrawn: number of differentially expressed genes
      ## numWdrawn: the number of genes differentially expressed in category
      hg <- .doHyperGInternal(numW, numB, numDrawn, numWdrawn, dir.over)
      c(list(n.in.set=n, n.drawn=numWdrawn), hg)
    }, by=c('collection', 'name')]
    setnames(out, 'p', 'pval')
    out[, padj := p.adjust(pval, 'BH')]
  }

  if (!as.dt) setDF(out)
  out
}

## -----------------------------------------------------------------------------
## These were copied from the Bioconductor collection package

## We envision the test as follows:
##
## The urn contains genes from the gene universe.  Genes annotated at a
## given collection term are white and the rest black.
##
## The number drawn is the size of the selected gene list.  The
## number of white drawn is the size of the intersection of the
## selected list and the genes annotated at the collection.
## Here's a diagram based on using GO as the collection:
##
##          inGO    notGO
##          White   Black
## selected  n11     n12
## not       n21     n22
##
## numW: number of genes in GO category
## numB: size of universe
## numDrawn: number of differentially expressed genes
## numWdrawn: the number of genes differentially expressed in category
.doHyperGInternal <- function(numW, numB, numDrawn, numWdrawn, over) {
  n21 <- numW - numWdrawn
  n12 <- numDrawn - numWdrawn
  n22 <- numB - n12

  odds_ratio <-  (numWdrawn * n22) / (n12 * n21)

  expected <- (numWdrawn + n12) * (numWdrawn + n21)
  expected <- expected / (numWdrawn + n12 + n21 + n22)

  if (over) {
    ## take the -1 because we want evidence for as extreme or more
    pvals <- phyper(numWdrawn - 1L, numW, numB,
                    numDrawn, lower.tail=FALSE)
  } else {
    pvals <- phyper(numWdrawn, numW, numB,
                    numDrawn, lower.tail=TRUE)
  }

  list(p=pvals, odds=odds_ratio, expected=expected)
}

