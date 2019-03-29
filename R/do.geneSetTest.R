##' @include validateInputs.R
NULL

validate.inputs.geneSetTest <- .validate.inputs.logFC.only
validate.x.geneSetTest <- validate.X

#' Worker function to run geneSetTest from within a multiGSEA pipeline
#'
#' **This function is not meant to be called directly.** It should only be
#' called internally within [multiGSEA()].
#'
#' @noRd
do.geneSetTest <- function(gsd, x, design, contrast = ncol(design),
                           score.by = c("t", "logFC", "pval"), logFC = NULL,
                           gs.idxs = NULL, ...) {
  score.by <- match.arg(score.by)
  stopifnot(is.conformed(gsd, x))
  if (!missing(design) && missing(contrast)) {
    contrast <- ncol(design)
  }

  stats <- extract_preranked_stats(x, design, contrast, score.by = score.by,
                                   logFC = logFC, ...)

  args <- list(...)
  call.args <- as.list(formals(limma::geneSetTest))
  for (arg in intersect(names(args), names(call.args))) {
    call.args[[arg]] <- args[[arg]]
  }

  call.args[['statistics']] <- stats

  if (is.null(gs.idxs)) {
    gs.idxs <- as.list(gsd, active.only = TRUE, value = "x.idx")
  }

  pvals <- sapply(gs.idxs, function(idx) {
    xargs <- call.args
    xargs[['index']] <- idx
    do.call(limma::geneSetTest, xargs)
  })

  out <- geneSets(gsd, as.dt = TRUE)[, list(collection, name)]
  kosher <- .gsdlist_conforms_to_gsd(gs.idxs, gsd, active.only = TRUE)
  stopifnot(kosher)
  out[, pval := pvals]
  out[, padj := p.adjust(pval, 'BH')]
  setattr(out, 'rawresult', TRUE)
}

mgres.geneSetTest <- function(res, gsd, ...) res

