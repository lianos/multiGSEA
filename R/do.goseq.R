## Note: I won't implement limma::goana since they do not accept their own
##       universe. The methodology in goseq is equivalent except for minor
##       detail when low number of DE genes (cf ?goana)
##' @include validateInputs.R
NULL


validate.inputs.goseq <- .validate.inputs.full.design
validate.x.goseq <- validate.X

##' Performs goseq analysis significance of gene set membership.
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
##' @param feature.bias a vector as long as \code{nrow(x)} that has the "bias"
##'   information for the features/genes tested (ie. vector of gene lengths).
##' @param method The method to use to calculate the unbiased category
##'   enrichment scores
##' @param repcnt Number of random samples to be calculated when random sampling
##'   is used. Ignored unless \code{method="Sampling"}.
##' @param use_genes_without_cat A boolean to indicate whether genes without a
##'   categorie should still be used. For example, a large number of gene may
##'   have no GO term annotated. If this option is set to FALSE, those genes
##'   will be ignored in the calculation of p-values (default behaviour). If
##'   this option is set to TRUE, then these genes will count towards the total
##'   number of genes outside the category being tested.
##' @param direction Same as direction in \code{GOstats}
##' @param plot.fit To plot (or not) the bias in selected genes vs.
##'   \code{feature.bias}.
##' @param logFC The logFC data.table from \code{calculateIndividualLogFC}
##' @param vm A voomed object, if \code{multiGSEA} was original called with
##'   a \code{DGEList}.
##' @param ... arguments to pass down into \code{calculateIndividualLogFC}
##' @return A data.table of goseq results. The "pval" column here refers to
##'   pval.over, for simplicity in other places.
do.goseq <- function(gsd, x, design, contrast=ncol(design),
                     feature.bias=create.glength.vector(x),
                     method="Wallenius",
                     repcnt=2000, use_genes_without_cat=TRUE,
                     direction=c('over', 'under'),
                     plot.fit=FALSE,
                     logFC=NULL,
                     feature.min.logFC=1,
                     feature.max.padj=0.10, vm=x, ...) {
  if (is(x, "DGEList")) {
    x <- vm
  }
  stopifnot(is.conformed(gsd, x))
  direction <- match.arg(direction)

  if (is.null(logFC)) {
    logFC <- calculateIndividualLogFC(x, design, contrast, ...)
  } else {
    is.logFC.like(logFC, x, as.error=TRUE)
  }

  if (is.null(logFC$hyperG.selected)) {
    logFC[, hyperG.selected := {
      padj <= feature.max.padj & abs(logFC) >= feature.min.logFC
    }]
  }

  drawn <- logFC[hyperG.selected == TRUE]$featureId
  res <- multiGSEA::goseq(gsd, drawn, rownames(x), feature.bias, method,
                          repcnt, use_genes_without_cat, plot.fit=plot.fit,
                          do.conform=FALSE)
  setnames(res, c('over_represented_pvalue', 'under_represented_pvalue'),
           c('pval', 'pval.under'))
  res[, padj := p.adjust(pval, 'BH')]
  res[, padj.under := p.adjust(pval.under, 'BH')]
}

##' Perform goseq Enrichment tests across a GeneSetDb.
##'
##' @export
##' @import goseq
##'
##' @param gsd The \code{GeneSetDb} object to run tests against
##' @param selected The ids of the selected features
##' @param universe The ids of the universe
##' @param feature.bias a vector as long as \code{nrow(x)} that has the "bias"
##'   information for the features/genes tested (ie. vector of gene lengths).
##' @param method The method to use to calculate the unbiased category
##'   enrichment scores
##' @param repcnt Number of random samples to be calculated when random sampling
##'   is used. Ignored unless \code{method="Sampling"}.
##' @param use_genes_without_cat A boolean to indicate whether genes without a
##'   categorie should still be used. For example, a large number of gene may
##'   have no GO term annotated. If this option is set to FALSE, those genes
##'   will be ignored in the calculation of p-values (default behaviour). If
##'   this option is set to TRUE, then these genes will count towards the total
##'   number of genes outside the category being tested.
##' @param do.conform By default \code{TRUE}: does some gymnastics to conform
##'   the \code{gsd} to the \code{universe} vector. This should neber be set
##'   to \code{FALSE}, but this parameter is here so that when this function
##'   is called from the \code{\link{multiGSEA}} codepath, we do not have to
##'   reconform the \code{GeneSetDb} object, because it has already been done.
##' @param active.only If \code{TRUE}, only "active" genesets are used
##' @param value The featureId types to extract from \code{gsd}
##' @return A \code{data.table} of results, similar to goseq output.
goseq <- function(gsd, selected, universe,
                  feature.bias=NULL,
                  method=c("Wallenius", "Sampling", "Hypergeometric"),
                  repcnt=2000, use_genes_without_cat=TRUE,
                  plot.fit=TRUE, do.conform=TRUE) {
  stopifnot(is(gsd, 'GeneSetDb'))
  stopifnot(is.character(selected))
  stopifnot(is.character(universe))
  selected <- unique(selected)
  universe <- unique(universe)
  method <- match.arg(method)

  if (missing(feature.bias)) {
    gmat <- matrix(rep(1, length(universe)),
                   nrow=length(universe),
                   dimnames=list(universe, NULL))
    feature.bias <- create.glength.vector(gmat)
  }

  ## This needs to be conformed to work
  if (!is.conformed(gsd) || !do.conform) {
    gsd <- conform(gsd, universe)
  }
  stopifnot(is.conformed(gsd, universe))

  gs <- geneSets(gsd, active.only=TRUE)
  gs[, category := paste(collection, name, sep=';;')]
  g2c <- as.data.frame(gsd, active.only=TRUE, value='x.id')
  g2c <- transform(g2c, category=paste(collection, name, sep=';;'))
  g2c <- g2c[, c('category', 'featureId')]
  selected <- intersect(selected, universe)
  de.genes <- setNames(integer(length(universe)), universe)
  de.genes[selected] <- 1L

  pwf <- goseq::nullp(de.genes, bias.data=feature.bias, plot.fit=plot.fit)
  res <- goseq::goseq(pwf, gene2cat=g2c, method=method, repcnt=repcnt,
                      use_genes_without_cat=use_genes_without_cat)

  ## resort output to match gs ordering
  res <- res[match(gs$category, res$category),]
  ##            !names(res) %in% c('category', 'numInCat')]
  ## out <- cbind(gs[, list(collection, name, n, n.drawn=res$numDEInCat)], res)
  out <- gs[, list(collection, name, n, n.drawn=res$numDEInCat)]
  rcols <- setdiff(names(res), c('category', 'numInCat', 'numDEInCat'))
  for (rcol in rcols) {
    out[, (rcol) := res[rcol]]
  }
  attr(out, 'pwf') <- pwf
  out
}

##' Takes a conformed GeneSetDb and creates a data.frame for goseq to test over
##'
##' @param gsd GeneSetDb to prep
##' @return \code{data.frame} for goseq's \code{gene2cat}
goseq.gene2cat <- function(gsd) {
  stopifnot(is(gsd, 'GeneSetDb'))
  stopifnot(is.conformed(gsd))

  gs <- geneSets(gsd)
  gs[, category := paste(collection, name, sep=';;')]

  gene2cat <- merge(gsd@db, featureIdMap(gsd), by='featureId')
  ## remove features/collections that are inactive
  gene2cat <- subset(gene2cat, !(is.na(x.id) | is.na(x.idx)))
  gene2cat[, category := paste(collection, name, sep=';;')]
  gene2cat <- subset(gene2cat, category %in% gs$category)
  as.data.frame(gene2cat[, list(category, featureId)], stringsAsFactors=FALSE)
}
