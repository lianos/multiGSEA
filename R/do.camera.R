##' @include validateInputs.R
NULL

validate.inputs.camera <- .validate.inputs.full.design
validate.x.camera <- validate.X

##' Runs camera on the given experiment.
##'
do.camera <- function(gsd, x, design, contrast=ncol(design), outdir=NULL,
                      use.cache=TRUE, weights=NULL, use.ranks=FALSE,
                      allow.neg.cor=TRUE, trend.var=FALSE, sort=FALSE,
                      ...) {
  stopifnot(is.conformed(gsd, x))
  extra.args <- c('use.ranks', 'allow.neg.cor', 'trend.var', 'sort')
  cache.fn <- cache.data.fn('camera', design, contrast, extra.args,
                            outdir=outdir, ext='rds')
  gs.idxs <- as.expression.indexes(gsd, value='x.idx')

  if (file.exists(cache.fn) && use.cache) {
    res <- readRDS(cache.fn)
    ## TODO: check the genesets returned from pvals to ensure that they match
    ##       the active genesets
  } else {
    res <- camera(x, gs.idxs, design, contrast, weights=weights,
                  use.ranks=use.ranks, allow.neg.cor=allow.neg.cor,
                  trend.var=trend.var, sort=FALSE, ...)
    if (is.character(outdir) && isTRUE(file.exists(outdir))) {
      saveRDS(res, cache.fn)
    }
  }

  out <- cbind(geneSets(gsd)[, list(collection, name)], as.data.table(res))
  out[, NGenes := NULL]
  setnames(out, c('PValue', 'FDR'), c('pval', 'padj'))
}
