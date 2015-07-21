##' @include validateInputs.R
NULL

validate.inputs.hyperGeometricTest <- .validate.inputs.full.design
validate.x.hyperGeometricTest <- validate.X

##' Performs a hypergeometric test for significance of gene set membership.
##'
##' Genes are selected for testing against each geneset by virture of them
##' passing a maximum FDR and minimum log fold change as perscribed by the
##' \code{min.logFC} and \code{max.padj} parameters, respectfully.
##'
##' Note that we are intentionally adding a hyperG.selected column by reference
##' so that this information is kicked back to the caller multiGSEA function
##' and included in downstream reporting.
##'
##' @param gsd The \code{\link{GeneSetDb}} for analysis
##' @param x The expression object
##' @param design Experimental design
##' @param contrast The contrast to test
##' @param direction Same as direction in \code{GOstats}
##' @param logFC The logFC data.table from \code{calculateIndividualLogFC}
##' @param ... arguments to pass down into \code{calculateIndividualLogFC}
do.hyperGeometricTest <- function(gsd, x, design, contrast=ncol(design),
                                  direction='over', logFC=NULL,
                                  feature.min.logFC=1,
                                  feature.max.padj=0.10, vm=x, ...) {
  if (is(x, "DGEList")) {
    x <- vm
  }
  stopifnot(is.conformed(gsd, x))

  if (is.null(logFC)) {
    logFC <- calculateIndividualLogFC(x, design, contrast, ...)
  } else {
    is.logFC.like(logFC, x, as.error=TRUE)
  }

  dir.over <- direction == 'over'

  if (is.null(logFC$hyperG.selected)) {
    logFC[, hyperG.selected := {
      padj <= feature.max.padj & abs(logFC) >= feature.min.logFC
    }]
  }

  ## These are the IDs we are "drawing" out of the bucket
  drawn <- logFC[hyperG.selected == TRUE]$featureId
  numDrawn <- length(drawn)     ## number of genes selected to be "interesting"
  if (numDrawn == 0) {
    warning("No selected genes in hyperGeometricTest", immediate.=TRUE)
    out <- geneSets(gsd)[, list(collection, id)]
    out[, pval := NA_real_]
    out[, odds := NA_real_]
    out[, expected := NA_real_]
    out[, padj := NA_real_]
    out[, padj.by.collection := NA_real_]
  } else {
    numB <- nrow(x)              ## number of genes in universe
    out <- geneSets(gsd)[, {
      ids <- featureIds(gsd, .BY[[1]], .BY[[2]])
      numW <- length(ids)
      Wdrawn <- intersect(ids, drawn)
      numWdrawn <- length(Wdrawn)
      hg <- .doHyperGInternal(numW, numB, numDrawn, numWdrawn, dir.over)
      c(list(n.in.set=n, n.drawn=numWdrawn), hg)
    }, by=c('collection', 'name')]
    setnames(out, 'p', 'pval')
    out[, padj := p.adjust(pval, 'BH')]
  }

  out
}

##' Identify the features that were "selected" for a given geneset for the
##' hyperGeometricTest for enrichment.
##'
##' @export
##'
##' @param x a \code{MultiGSEAResult}
##' @param i the collection of the gene set
##' @param j the id of the geneset
##'
##' @return A character vector of features in the expression object that were
##'   identified as part of the geneset that was "selected"
selectedByHyperG <- function(x, i, j) {
  stopifnot(is(x, 'MultiGSEAResult'))
  if (!'hyperGeometricTest' %in% resultNames(x)) {
    stop("A hypergeometric test was not performed")
  }
  fids <- featureIds(x@gsd, i, j, value='x.id')
  intersect(fids, subset(x@logFC, hyperG.selected == TRUE)$featureId)
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

