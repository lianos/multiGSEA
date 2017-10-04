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
do.geneSetTest <- function(gsd, x, design, contrast, outdir=NULL,
                           robust.fit=FALSE, robust.eBayes=FALSE,
                           logFC=NULL, score.by=c('t', 'logFC', 'pval'),
                           gs.idxs=as.list(gsd, active.only=TRUE, value='x.idx'),
                           ...) {
  score.by <- match.arg(score.by)
  stopifnot(is.conformed(gsd, x))
  if (!missing(design) && missing(contrast)) {
    contrast <- ncol(design)
  }

  if (ncol(x) > 1) {
    if (is.null(logFC)) {
      logFC <- calculateIndividualLogFC(x, design, contrast, robust.fit,
                                        robust.eBayes, ..., as.dt=TRUE)
    } else {
      is.logFC.like(logFC, x, as.error=TRUE)
    }
    ## t will be NA if statistics were computed using edgeR from a DGEList
    if (score.by == 't' && any(is.na(logFC[['t']]))) {
      warning("t statistics not found in dge results, switching to logFC",
              immediate.=TRUE)
      score.by <- 'logFC'
    }
    stats <- setNames(logFC[[score.by]], logFC$featureId)
  } else {
    ## This already a column matrix of precomputed things (logFC, perhaps)
    ## to rank
    stats <- setNames(x[, 1L], rownames(x))
  }

  if (any(is.na(stats))) {
    stop("NA values are found in the stats vector for geneSetTest")
  }

  args <- list(...)
  call.args <- as.list(formals(limma::geneSetTest))
  for (arg in intersect(names(args), names(call.args))) {
    call.args[[arg]] <- args[[arg]]
  }

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

  out <- geneSets(gsd, as.dt=TRUE)[, list(collection, name)]
  out[, pval := pvals]
  out[, padj := p.adjust(pval, 'BH')]
  setattr(out, 'rawresult', TRUE)
}

mgres.geneSetTest <- function(res, gsd, ...) res

