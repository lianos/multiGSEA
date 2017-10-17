##' @include validateInputs.R
NULL

validate.inputs.geneSetTest <- .validate.inputs.logFC.only
validate.x.geneSetTest <- validate.X

##' Worker function to run geneSetTest from within a multiGSEA pipeline
##'
##' @description
##'
##' \strong{This function is not meant to be called directly, it should only be
##' called internally within \code{multiGSEA}}
do.geneSetTest <- function(gsd, x, design, contrast,
                           score.by=c('t', 'logFC', 'pval'), logFC=NULL,
                           gs.idxs=as.list(gsd, active.only=TRUE, value='x.idx'),
                           ...) {
  score.by <- match.arg(score.by)
  stopifnot(is.conformed(gsd, x))
  if (!missing(design) && missing(contrast)) {
    contrast <- ncol(design)
  }

  stats <- extract_preranked_stats(x, design, contrast, score.by=score.by,
                                   logFC=logFC, ...)

  args <- list(...)
  call.args <- as.list(formals(limma::geneSetTest))
  for (arg in intersect(names(args), names(call.args))) {
    call.args[[arg]] <- args[[arg]]
  }

  call.args[['statistics']] <- stats

  pvals <- sapply(gs.idxs, function(idx) {
    xargs <- call.args
    xargs[['index']] <- idx
    do.call(limma::geneSetTest, xargs)
  })

  out <- geneSets(gsd, as.dt=TRUE)[, list(collection, name)]
  out[, pval := pvals]
  out[, padj := p.adjust(pval, 'BH')]
  setattr(out, 'rawresult', TRUE)
}

mgres.geneSetTest <- function(res, gsd, ...) res

