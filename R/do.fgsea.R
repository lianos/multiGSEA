#' @include validateInputs.R
NULL

validate.inputs.fgsea <- .validate.inputs.preranked
validate.x.fgsea <- validate.X

#' Runs GSEA on a pre-ranked list of differential expression statistcis with fgsea
#'
#' Techinically this doesn't require a full design to run, but it proves
#' expedient for a quick turnaround.
#'
#' Note that fgsea isn't added to import or suggests because rescomp's
#' compiler can't handle this yet.
#'
#' fgsea expects the pathways to be describes as a list of character vectors
#' where the vectors are gene ids, and the names of the pre-ranked vector
#' correspond to those IDs
#'
#' **This function is not meant to be called directly.** It should only be
#' called internally within [multiGSEA()].
#'
#' @param gsd The [GeneSetDb()] for analysis
#' @inheritParams calculateIndividualLogFC
#' @param ... arguments to pass down into [calculateIndividualLogFC()]
#' @return A data.table of fgsea results.
do.fgsea <- function(gsd, x, design, contrast=ncol(design),
                     minSize=15, maxSize=15, nperm=10000, gseaParam=1,
                     score.by=c('t', 'logFC', 'pval'), logFC=NULL,
                     gs.idxs=as.list(gsd, active.only=TRUE, value='x.idx'),
                     ...) {
  score.by <- match.arg(score.by)
  if (!requireNamespace('fgsea', quietly=TRUE)) {
    stop("The Bioconductor fgsea package is required for this functionality")
  }
  stopifnot(is.conformed(gsd, x))

  stats <- extract_preranked_stats(x, design, contrast, score.by=score.by,
                                   logFC=logFC, ...)

  ## fgsea call needs minSize and maxSize params, we set this to whatever
  ## the min/max active geneset sizes are, since this was already specified
  ## in the internally called conform(gsd, x, min.gs.size, max.gs.size) call
  gs.size <- range(subset(geneSets(gsd), active)$n)

  ## fgsea function wans a list of gene identifiers for pathway definition
  pathways <- lapply(gs.idxs, function(idxs) names(stats)[idxs])

  res <- fgsea::fgsea(pathways, stats, nperm, minSize=gs.size[1L],
                      maxSize=gs.size[2L], gseaParam=gseaParam)
  setattr(res, 'rawresult', TRUE)
}

mgres.fgsea <- function(res, gsd, ...) {
  if (!isTRUE(attr(res, 'rawresult'))) return(res)
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
  gs <- geneSets(gsd, as.dt=TRUE)[, list(collection, name)]
  stopifnot(all.equal(res$pathway, encode_gskey(gs)))
  out <- cbind(gs, as.data.table(res[, -1L, drop=FALSE]))
  out
}
