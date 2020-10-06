#' @include validateInputs.R
NULL

validate.inputs.fgsea <- .validate.inputs.preranked
validate.x.fgsea <- validate.X

#' Runs GSEA on a pre-ranked list of differential expression statistcis with fgsea
#'
#' fgsea expects the pathways to be describes as a list of character vectors
#' where the vectors are gene ids, and the names of the pre-ranked vector
#' correspond to those IDs
#'
#' Deafults fgsea functionality was changed to [fgsea::fgseaMultilevel()] in
#' 1.13.2 and the default parameters here reflect that. If you want to use the
#' original "fgseaSimple" you have to pass in `use.fgsea.simple = TRUE` in your
#' `multiGSEA()` call.
#'
#' `minSize` and `maxSize` are already set by the `conform` logic that was
#' specified in the call to [multiGSEA()] via the `min.gs.size` and
#' `max.gs.size` parameters.
#'
#' **This function is not meant to be called directly.** It should only be
#' called internally within [multiGSEA()].
#'
#' @param gsd The [GeneSetDb()] for analysis
#' @inheritParams calculateIndividualLogFC
#' @param ... arguments to pass down into [calculateIndividualLogFC()]
#' @return A data.table of fgsea results.
do.fgsea <- function(gsd, x, design, contrast = ncol(design),
                     sampleSize = 101, eps = 1e-10,
                     scoreType = c("std", "pos", "neg"),
                     nproc = 0, gseaParam = 1, nPermSimple = 1000,
                     absEps = NULL, use.fgsea.simple = FALSE,
                     score.by = c('t', 'logFC', 'pval'), logFC = NULL,
                     gs.idxs = as.list(gsd, active.only=TRUE, value='x.idx'),
                     ...) {

  scoreType <- match.arg(scoreType)
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

  if (isTRUE(use.fgsea.simple)) {
    res <- fgsea::fgseaSimple(
      pathways, stats, nperm = nPermSimple,
      minSize = gs.size[1L], maxSize = gs.size[2L],
      scoreType = scoreType, nproc = nproc,
      gseaParam = gseaParam)
  } else {
    res <- fgsea::fgseaMultilevel(
      pathways, stats, sampleSize = sampleSize,
      minSize = gs.size[1L], maxSize = gs.size[2L],
      eps = eps, scoreType = scoreType, nproc = nproc,
      gseaParam = gseaParam,
      # nPermSimple = nPermSimple, # TODO: uncomment for next version of bioc
      absEps = absEps)
  }
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
