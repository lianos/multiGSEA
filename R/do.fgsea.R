##' @include validateInputs.R
NULL

validate.inputs.fgsea <- .validate.inputs.preranked
validate.x.fgsea <- validate.X

##' Runs GSEA on a pre-ranked list of differential expression statistcis with fgsea
##'
##' Techinically this doesn't require a full design to run, but it proves
##' expedient for a quick turnaround.
##'
##' Note that fgsea isn't added to import or suggests because rescomp's
##' compiler can't handle this yet.
##'
##' fgsea expects the pathways to be describes as a list of character vectors
##' where the vectors are gene ids, and the names of the pre-ranked vector
##' correspond to those IDs
##'
##' @param gsd The \code{\link{GeneSetDb}} for analysis
##' @inheritParams calculateIndividualLogFC
##' @param ... arguments to pass down into \code{calculateIndividualLogFC}
##' @return A data.table of fgsea results.
do.fgsea <- function(gsd, x, design, contrast=ncol(design),
                     minSize=15, maxSize=15, nperm=10000, gseaParam=1,
                     score.by=c('t', 'logFC', 'pval'), use.treat=FALSE,
                     feature.min.logFC=if (use.treat) log2(1.25) else 1,
                     feature.max.padj=0.10,
                     gs.idxs=as.list(gsd, active.only=TRUE, value='x.idx'),
                     logFC=NULL,
                     ...) {
  score.by <- match.arg(score.by)
  if (!requireNamespace('fgsea')) {
    stop("The Bioconductor fgsea is required for this functionality")
  }
  if (!is.conformed(gsd, x)) {
    gsd <- conform(gsd, x, min.gs.size=minSize, max.gs.size=maxSize)
    gs.idxs <- as.list(gsd, active.only=TRUE, value='x.idx')
  }

  ## The call to calculateIndividualLogFC in multiGSEA puts "the right stuff"
  ## in the logFC data.table, so we just need that.
  ranks <- setNames(logFC[[score.by]], logFC[['featureId']])
  gs <- subset(geneSets(gsd), active)
  gs.range <- range(gs$n, na.rm=TRUE)

  ## fgsea function wans a list of gene identifiers for pathway definition
  pathways <- lapply(gs.idxs, function(idxs) names(ranks)[idxs])

  res <- fgsea::fgsea(pathways, ranks, nperm, minSize=gs.range[1L],
                      maxSize=gs.range[2L], gseaParam=gseaParam)

  # fgsea will sometimes return 0 hits for a pathway. Need to reconstruct a fix!
  # missed <- setdiff(names(gs.idxs), res$pathway)
  # if (length(missed)) {
  #   warning(length(missed), " pathways missed in fgsea!")
  #   lists <- replicate(length(missed), list())
  #   addme <- data.table(pathway=missed, leadingEdge=lists)
  #   res <- rbindlist(list(res, addme), fill=TRUE, use.names=TRUE)
  # }
  # xref <- match(names(gs.idxs), res$pathway)
  # res <- res[xref]
  gs <- geneSets(gsd, .external=FALSE)[, list(collection, name)]
  stopifnot(all.equal(res$pathway, paste(gs$collection, gs$name, sep=';;')))
  out <- cbind(gs, res[, -1, with=FALSE])
  out
}

