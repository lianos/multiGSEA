##' @include validateInputs.R
NULL

validate.inputs.fgsea <- .validate.inputs.full.design
validate.x.fgsea <- validate.X

##' Runs GSEA on a pre-ranked list of differential expression statistcis with fgsea
##'
##' Techinically this doesn't require a full design to run, but it proves
##' expedient for a quick turnaround.
##'
##' Note that fgsea isn't added to import or suggests because rescomp's
##' compiler can't handle this yet.
##'
##' @param gsd The \code{\link{GeneSetDb}} for analysis
##' @inheritParams calculateIndividualLogFC
##' @param ... arguments to pass down into \code{calculateIndividualLogFC}
##' @return A data.table of fgsea results. The "pval" column here refers to
##'   pval.over, for simplicity in other places. If \code{split.updown=TRUE},
##'   a list of data.table's are returned named 'goseq', 'goseq.up', and
##'   'goseq.down' which are the results of running goseq three independent
##'   times.
do.fgsea <- function(gsd, x, design, contrast=ncol(design),
                     minSize=15, maxSize=15, nperm=10000, gseaParam=1,
                     score.by=c('t', 'logFC', 'pval'), use.treat=FALSE,
                     feature.min.logFC=if (use.treat) log2(1.25) else 1,
                     feature.max.padj=0.10, gs.idxs=NULL, logFC=NULL,
                     ...) {
  if (!requireNamespace('fgsea')) {
    stop("The Bioconductor fgsea is required for this functionality")
  }
  if (!is.conformed(gsd, x)) {
    gsd <- conform(gsd, x, min.gs.size=minSize, max.gs.size=maxSize)
    gs.idxs <- NULL
  }
  if (is.null(gsd.idxs) || is.numeric(gsd.idxs[[1]])) {
    gs.idxs <- as.list(gsd, value='x.id')
  }

  if (is.null(logFC)) {
    logFC <- calculateIndividualLogFC(x, design, contrast, use.treat=use.treat,
                                      treat.lfc=feature.min.logFC, ...,
                                      .external=FALSE)
    if (use.treat) {
      logFC[, significant := padj <= feature.max.padj]
    }
  } else {
    is.logFC.like(logFC, x, as.error=TRUE)
  }

  ranks <- logFC[[score.by]]
  if (is.null(ranks) || any(is.na(ranks))) {
    ranks <- logFC[['logFC']]
  }
  names(ranks) <- logFC$featureId

  gs <- subset(geneSets(gsd), active)
  minSize <- min(gs$n, na.rm=TRUE)
  maxSize <- max(gs$n, na.rm=TRUE)

  ## wants stats to be a named vector of values whose names are found in
  ## pathwas (gs.idxs). The code suggest that pathways can be numeric and
  ## this would be kosher, let's see.
  res <- fgsea(gs.idxs, ranks, nperm, minSize, maxSize, gseaParam=gseaParam)
  ## TODO: transform fgsea result to legit multiGSEA result
}
