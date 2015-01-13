##' @include validateInputs.R
NULL

validate.inputs.roast <- .validate.inputs.full.design

do.roast <- function(gsd, x, design, contrast=ncol(design), outdir=NULL,
                     use.cache=TRUE, set.statistic='mean', gene.weights=NULL,
                     array.weights=NULL, weights=NULL, block=NULL, correlation,
                     var.prior=NULL, df.prior=NULL, trend.var=FALSE, nrot=10000,
                     approx.zscore=TRUE, ...) {
  stopifnot(is.conformed(gsd, x))
  extra.args <- c('set.statistic', 'gene.weights', 'array.weights', 'weights',
                  'block', 'correlation', 'var.prior', 'df.prior', 'trend.var',
                  'nrot', 'approx.zscore')
  extra.args <- c('set.statistic', 'nrot', 'approx.zscore')
  cache.fn <- cache.data.fn('roast', design, contrast, extra.args,
                            outdir=outdir, ext='rds')
  gs.idxs <- as.expression.indexes(gsd, value='x.idx')
  if (file.exists(cache.fn) && use.cache) {
    res <- readRDS(cache.fn)
    ## TODO: check the genesets returned from pvals to ensure that they match
    ##       the active genesets
  } else {
    res <- mroast(x, gs.idxs, design, contrast, set.statistic=set.statistic,
                  gene.weights=gene.weights, array.weights=array.weights,
                  weights=weights, block=block, correlation=correlation,
                  var.prior=var.prior, df.prior=df.prior, trend.var=trend.var,
                  nrot=nrot, approx.zscore=approx.zscore, sort='none', ...)
    if (is.character(outdir) && isTRUE(file.exists(outdir))) {
      saveRDS(res, cache.fn)
    }
  }

  out <- cbind(geneSets(gsd)[, list(collection, name)], as.data.table(res))
  out[, NGenes := NULL]

  setnames(out,
           c('PValue', 'FDR', 'PValue.Mixed', 'FDR.Mixed'),
           c('pval', 'padj', 'pval.mixed', 'padj.mixed'))
  out[, rank := rank(pval)]
  out[, rank.up := Inf]
  out[, rank.down := Inf]
  out[Direction == 'Up', rank.up := rank(pval)]
  out[Direction == 'Down', rank.down := rank(pval)]
}

