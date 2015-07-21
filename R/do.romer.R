##' @include validateInputs.R
NULL

validate.inputs.romer <- .validate.inputs.full.design
validate.x.romer <- validate.DGEList

do.romer <- function(gsd, x, design, contrast=ncol(design), outdir=NULL,
                     use.cache=TRUE, array.weights=NULL, block=NULL,
                     correlation=NULL, set.statistic="mean", nrot=10000,
                     shrink.resid=TRUE, ...) {
  stopifnot(is.conformed(gsd, x))
  extra.args <- c('array.weights', 'block', 'correlation', 'set.statistic',
                  'weights', 'nrot', 'shrink.resid')
  extra.args <- c('set.statistic', 'nrot', 'shrink.resid')
  cache.fn <- cache.data.fn('romer', design, contrast, extra.args,
                            outdir=outdir, ext='rds')
  gs.idxs <- as.expression.indexes(gsd, value='x.idx')
  if (file.exists(cache.fn) && use.cache) {
    res <- readRDS(cache.fn)
    ## TODO: check the genesets returned from pvals to ensure that they match
    ##       the active genesets
  } else {
    res <- romer(x, gs.idxs, design, contrast,
                 array.weights=array.weights,
                 block=block,
                 correlation=correlation,
                 set.statistic=set.statistic,
                 nrot=nrot,
                 shrink.resid=shrink.resid, ...)
    if (is.character(outdir) && isTRUE(file.exists(outdir))) {
      saveRDS(res, cache.fn)
    }
  }

  out <- cbind(geneSets(gsd)[, list(collection, name)], as.data.table(res))
  out[, NGenes := NULL]

  setnames(out,
           c('Up', 'Down', 'Mixed'),
           c('pval.up', 'pval.down', 'pval'))

  out[, padj := p.adjust(pval)]
  out[, padj.up := p.adjust(pval.up)]
  out[, padj.down := p.adjust(pval.down)]

  out
}

