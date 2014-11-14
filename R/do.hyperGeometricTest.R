##' Performs a hypergeometric test for significance of gene set membership.
##'
##' Genes are selected for testing against each geneset by virture of them
##' passing a maximum FDR and minimum log fold change as perscribed by the
##' \code{min.logFC} and \code{max.padj} parameters, respectfully.
##'
##' @param x The expression object
##' @param gs.table The \code{\link{GeneSetTable}} for analysis
##' @param design Experimental design
##' @param contrast The contrast to test
##' @param direction Same as direction in \code{GOstats}
##' @param logFC.stats The logFC stats table for the given contrast
##' @param logFC.stats The topTable from limma
do.hyperGeometricTest <- function(x, gs.table, design, contrast,
                                  direction='over', logFC.stats=NULL,
                                  robust.fit=FALSE, robust.eBayes=FALSE,
                                  min.logFC=1, max.padj=0.10, ...) {
  if (!is.data.frame(logFC.stats)) {
    logFC.stats <- calculateIndividualLogFC(x, design, contrast, robust.fit,
                                            robust.eBayes, provide='table',
                                            ...)
  }
  req  <- c('logFC', 'pval', 'padj')
  if (length(missed <- setdiff(names(logFC.stats), req))) {
    stop("Requried columns in logFC.stats are missing: ",
         paste(missed, collapse=','))
  }
  if (!all.equal(rownames(x), rownames(logFC.stats))) {
    stop("rownames for logFC.stats data.frame do not look right.")
  }

  dir.over <- direction == 'over'

  ## Logical vector to indicate which results are significant
  significant <- with(logFC.stats, padj <= max.padj & abs(logFC) >= min.logFC)
  numB <- nrow(x)              ## number of genes in universe
  numDrawn <- sum(significant) ## number of genes selected for differential expr

  res <- lapply(gs.table@table$membership, function(idx) {
    numW <- sum(idx)                  ## Number of genes in geneset
    numWdrawn <- sum(selected & idx)  ## Number of these we selected by dge
    ## n12 <- numW - n11
    ## n21 <- gs.size - n11
    ## n22 <- numB - gs.size
    as.data.table(.doHyperGInternal(numW, numB, numDrawn, numWdrawn, dir.over))
  })
  res <- rbindlist(res)
  setnames(res, 'p', 'pval')

  out <- cbind(gs.table@table[, list(group, id)], hg.resuts)
  out[, padj := p.adjust(pval, 'BH')]
  out[, padj.by.group := p.adjust(pval, 'BH'), by='group']
  out[, feature.id := lapply(membership, function(m) f.ids[m])]
}

## -----------------------------------------------------------------------------
## These were copied from the Bioconductor Category package

## We envision the test as follows:
##
## The urn contains genes from the gene universe.  Genes annotated at a
## given category term are white and the rest black.
##
## The number drawn is the size of the selected gene list.  The
## number of white drawn is the size of the intersection of the
## selected list and the genes annotated at the category.
## Here's a diagram based on using GO as the category:
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

