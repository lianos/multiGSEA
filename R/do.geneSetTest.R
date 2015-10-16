##' @include validateInputs.R
NULL

validate.inputs.geneSetTest <- .validate.inputs.logFC.only
validate.x.geneSetTest <- validate.X

do.geneSetTest <- function(gsd, x, design, contrast, outdir=NULL,
                           robust.fit=FALSE, robust.eBayes=FALSE,
                           logFC=NULL, score.by=c('t', 'logFC'), ...) {
  score.by <- match.arg(score.by)
  stopifnot(is.conformed(gsd, x))
  if (!missing(design) && missing(contrast)) {
    contrast <- ncol(design)
  }

  if (ncol(x) > 1) {
    if (is.null(logFC)) {
      logFC <- calculateIndividualLogFC(x, design, contrast, robust.fit,
                                        robust.eBayes, ...)
    } else {
      is.logFC.like(logFC, x, as.error=TRUE)
    }
    stats <- setNames(logFC[[score.by]], logFC$featureId)
  } else {
    ## This already a column matrix of precomputed things (logFC, perhaps)
    ## to rank
    stats <- setNames(x[, 1L], rownames(x))
  }

  args <- list(...)
  call.args <- as.list(formals(limma::geneSetTest))
  for (arg in intersect(names(args), names(call.args))) {
    call.args[[arg]] <- args[[arg]]
  }

  gs.idxs <- as.expression.indexes(gsd, value='x.idx')
  stats <- stats[rownames(x)]

  ## Double checking that feature <-> expression IDs match up here. This is
  ## supposed to be redundant, but measure twice & cut once.
  if (any(is.na(stats))
      || length(stats) != nrow(x)
      || any(rownames(x) != names(stats))) {
    stop("individual stats do not match up with rownames of x")
  }

  call.args[['statistics']] <- stats

  pvals <- sapply(gs.idxs, function(idx) {
    xargs <- call.args
    xargs[['index']] <- idx
    do.call(limma::geneSetTest, xargs)
  })

  out <- geneSets(gsd)[, list(collection, name)]
  out[, pval := pvals]
  out[, padj := p.adjust(pval, 'BH')]
}
