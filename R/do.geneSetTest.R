##' @include validateInputs.R
NULL

validate.inputs.geneSetTest <- .validate.inputs.logFC.only

do.geneSetTest <- function(gsd, x, design, contrast, outdir=NULL,
                           use.cache=TRUE,
                           alternative="mixed", type="auto", ranks.only=TRUE,
                           nsim=9999, robust.fit=FALSE, robust.eBayes=FALSE,
                           logFC=NULL, score.by=c('t', 'logFC'), mc.cores=1L,
                           ...) {
  score.by <- match.arg(score.by)
  stopifnot(is.conformed(gsd, x))
  extra.args <- c('alternative', 'type', 'ranks.only')
  if (!missing(design) && missing(contrast)) {
    contrast <- ncol(design)
  }
  cache.fn <- cache.data.fn('geneSetTest', design, contrast, extra.args,
                            outdir=outdir, ext='rds')
  if (file.exists(cache.fn) && use.cache) {
    pvals <- readRDS(cache.fn)
    ## TODO: check the genesets returned from pvals to ensure that they match
    ##       the active genesets
  } else {
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

    gs.idxs <- as.expression.indexes(gsd, value='x.idx')
    stats <- stats[rownames(x)]

    ## Double checking that feature <-> expression IDs match up here. This is
    ## supposed to be redundant, but measure twice & cut once.
    if (any(is.na(stats))
        || length(stats) != nrow(x)
        || any(rownames(x) != names(stats))) {
      stop("individual stats do not match up with rownames of x")
    }
    pvals <- sapply(gs.idxs, function(idx) {
      geneSetTest(idx, stats, alternative, type, ranks.only, nsim)
    })

    if (is.character(outdir) && isTRUE(file.exists(outdir))) {
      saveRDS(pvals, cache.fn)
    }
  }

  out <- geneSets(gsd)[, list(collection, name)]
  out[, pval := pvals]
  out[, padj := p.adjust(pval, 'BH')]
}
